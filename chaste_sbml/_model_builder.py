"""Translation of a loaded SBML model into the internal representation used by the templates.

ModelBuilder walks the libsbml model -- rules, ODEs, species, parameters, reactions, functions
and events -- and populates the collections (state variables, parameters, derived quantities,
equations, ...) that the code templates consume. It uses a NameManager for collision-free
synthetic names and validates the result for residual identifier conflicts.
"""

import re
from math import isnan
from typing import TYPE_CHECKING, Optional

from libsbml import (
    AST_RELATIONAL_EQ,
    AST_RELATIONAL_GEQ,
    AST_RELATIONAL_GT,
    AST_RELATIONAL_LEQ,
    AST_RELATIONAL_LT,
    AST_RELATIONAL_NEQ,
    SBML_ALGEBRAIC_RULE,
    SBML_ASSIGNMENT_RULE,
    SBML_RATE_RULE,
    SBMLTransforms,
    formulaToL3String,
    parseL3Formula,
)

from ._config import (
    AMOUNT_PREFIX,
    CHASTE_PREFIX,
    CONCENTRATION_PREFIX,
    DERIVATIVE_PREFIX,
    DERIVATIVE_SUFFIX,
    INITIAL_ASSIGNMENT_PREFIX,
    NON_DIM_UNITS,
    PLACEHOLDER_STATE_ID,
    PREFIX_SEP,
    DerivedQuantityKind,
    EquationType,
    EventType,
    VarType,
)
from ._expressions import collect_ast_names, formula_to_string, substitute_ast_names
from ._names import CHASTE_RESERVED_NAMES, NameConflictError, find_name_conflicts
from ._records import (
    DerivedQuantity,
    Equation,
    Event,
    EventAssignment,
    Function,
    InitialAssignment,
    LocalParameter,
    Parameter,
    Reaction,
    Rule,
    StateVariable,
)
from ._sbml_reader import get_compartment_size, get_function_definition_arguments

if TYPE_CHECKING:
    from typing import Any

    from libsbml import ASTNode
    from libsbml import Reaction as SbmlReaction
    from libsbml import Species, SpeciesReference


class ModelBuilder:
    """Builds the internal model representation from a loaded SBML model."""

    def __init__(self, sbml_model, name_manager, model_name: str) -> None:
        """Translate a loaded SBML model into the internal representation for the templates.

        :param sbml_model: The libsbml Model to translate.
        :param name_manager: The NameManager providing collision-free synthetic names.
        :param model_name: The model name, used in conflict error messages.
        """
        self._sbml_model = sbml_model
        self._sbml_compartments = sbml_model.getListOfCompartments()
        self._sbml_events = sbml_model.getListOfEvents()
        self._sbml_function_definitions = sbml_model.getListOfFunctionDefinitions()
        self._sbml_initial_assignments = sbml_model.getListOfInitialAssignments()
        self._sbml_parameters = sbml_model.getListOfParameters()
        self._sbml_reactions = sbml_model.getListOfReactions()
        self._sbml_rules = sbml_model.getListOfRules()
        self._sbml_species = sbml_model.getListOfSpecies()
        self._sbml_unit_definitions = sbml_model.getListOfUnitDefinitions()
        self._names = name_manager
        self._model_name = model_name

    def template_data(self) -> dict[str, "Any"]:
        """Return the built collections as a template-variables dict."""
        return dict(
            assignment_rules=self._assignment_rules,
            derived_quantities=self._derived_quantities,
            equations=self._equations,
            events=self._events,
            functions=self._functions,
            parameters=self._parameters,
            reactions=self._reactions,
            state_variables=self._state_variables,
            stoichiometry_variables=self._stoichiometry_variables,
        )

    def _add_amount(self, species: "Species") -> None:
        """Add an amount derived quantity variable to the template variables.

        :param species: The related Species.
        """
        species_id = species.getId()
        compartment_id = species.getCompartment()

        amt_id = self._names.reserve(AMOUNT_PREFIX + PREFIX_SEP + species_id)
        amt_label = f"Amount of {species_id}"
        amt_units = NON_DIM_UNITS  # TODO: Use correct units
        amt_rhs = f"{species_id} * {compartment_id}"
        amt_math = parseL3Formula(amt_rhs)

        self._add_derived_quantity(amt_id, amt_label, None, amt_units, kind=DerivedQuantityKind.CONVERSION)
        self._add_equation(var=amt_id, math=amt_math, eq_type=EquationType.CONVERSION)

    def _add_concentration(self, species: "Species") -> None:
        """Add a concentration derived quantity variable to the template variables.

        For a species stored as an amount (hasOnlySubstanceUnits), the concentration is the
        amount divided by its compartment size. Mirrors _add_amount and reuses the same
        conversion machinery (computed in ComputeDerivedQuantities, not stored as a member).

        :param species: The related Species.
        """
        species_id = species.getId()
        compartment_id = species.getCompartment()

        conc_id = self._names.reserve(CONCENTRATION_PREFIX + PREFIX_SEP + species_id)
        conc_label = f"Concentration of {species_id}"
        conc_units = NON_DIM_UNITS  # TODO: Use correct units
        conc_rhs = f"{species_id} / {compartment_id}"
        conc_math = parseL3Formula(conc_rhs)

        self._add_derived_quantity(conc_id, conc_label, None, conc_units, kind=DerivedQuantityKind.CONVERSION)
        self._add_equation(var=conc_id, math=conc_math, eq_type=EquationType.CONVERSION)

    def _add_assignment_rule(self, id_: str, label: str, var: str, math: "ASTNode") -> None:
        """Add an assignment rule to the template variables.

        :param id_: The rule ID.
        :param label: The rule description.
        :param var: The variable being defined.
        :param math: The right-hand side of the rule as an ASTNode.
        """
        self._assignment_rules.append(Rule(id=id_, label=label, var=var, math=math))

    def _add_derived_quantity(
        self,
        id_: str,
        label: str,
        initial_value: Optional[float],
        units: str = NON_DIM_UNITS,
        kind: DerivedQuantityKind = DerivedQuantityKind.NORMAL,
    ) -> None:
        """Add a derived quantity to the template variables.

        :param id_: The variable ID.
        :param label: The variable description.
        :param initial_value: The variable initial value.
        :param units: The variable units.
        :param kind: The kind of derived quantity. CONVERSION is an amount/concentration conversion
            (computed in ComputeDerivedQuantities rather than stored as a member); REACTION is a
            reaction flux exposed as an output -- its member is already declared by the reaction
            loop, so it is not re-declared and its variable type stays VarType.REACTION.
        """
        self._derived_quantities.append(
            DerivedQuantity(
                id=id_,
                label=label,
                index=len(self._derived_quantities),
                initial_value=initial_value,
                units=units,
                kind=kind,
                name=self._names.sbml_name(id_),
            )
        )

    def _add_equation(
        self,
        var: str,
        math: "ASTNode",
        eq_type: EquationType = EquationType.UNKNOWN,
        local_parameters: Optional[list["LocalParameter"]] = None,
    ) -> None:
        """Add an equation to the template variables.

        These are expected to be equations of the form `var = rhs`.

        :param var: The variable being defined.
        :param math: The right-hand side of the equation as an AST.
        :param eq_type: The type of the equation.
        :param local_parameters: The kinetic-law LocalParameter records in scope for the equation.
        """
        self._equations.append(Equation(var=var, math=math, type=eq_type, local_parameters=local_parameters))

    def _add_event(
        self,
        label: str,
        trigger: str,
        assignments: list["EventAssignment"],
        distance: str,
        event_type: EventType,
        initial_satisfied: bool = True,
        priority: Optional[str] = None,
    ) -> None:
        """Add an event to the template variables.

        :param label: The event description.
        :param trigger: The event trigger formula.
        :param assignments: The event assignments.
        :param distance: The distance for the event trigger.
        :param event_type: The type of the event (e.g., cell division).
        :param initial_satisfied: Whether the event starts as satisfied (from SBML trigger initialValue).
        :param priority: The event priority formula, or None if the event has no priority.
        """
        self._events.append(
            Event(
                label=label,
                index=len(self._events),
                trigger=trigger,
                assignments=assignments,
                distance=distance,
                type=event_type,
                initial_satisfied=initial_satisfied,
                priority=priority,
            )
        )

    def _add_function(self, id_: str, label: str, args: str, body: str) -> None:
        """Add a function to the template variables.

        :param id_: The function ID.
        :param label: The function description.
        :param args: The function arguments.
        :param body: The function body.
        """
        self._functions.append(Function(id=id_, label=label, index=len(self._functions), args=args, body=body))

    def _add_initial_assignment(self, id_: str, label: str, var: str, math: Optional["ASTNode"] = None) -> None:
        """Add an initial assignment to the template variables.

        :param id_: The rule ID.
        :param label: The rule description.
        :param var: The variable being assigned.
        :param math: The right-hand side of the rule as an ASTNode.
        """
        self._initial_assignments.append(InitialAssignment(id=id_, label=label, var=var, math=math))

    def _add_parameter(
        self,
        id_: str,
        label: str,
        initial_value: Optional[float] = None,
        units: str = NON_DIM_UNITS,
        is_const: bool = False,
    ) -> None:
        """Add a parameter variable to the template variables.

        :param id_: The variable ID.
        :param label: The variable description.
        :param initial_value: The variable initial value.
        :param units: The variable units.
        :param is_const: True if variable is constant.
        """
        self._parameters.append(
            Parameter(
                index=len(self._parameters),
                id=id_,
                is_const=is_const,
                label=label,
                initial_value=initial_value,
                units=units,
                name=self._names.sbml_name(id_),
            )
        )

    def _add_rate_rule(self, id_: str, label: str, var: str, math: "ASTNode") -> None:
        """Add a rate rule to the template variables.

        :param id_: The rule ID.
        :param label: The rule description.
        :param var: The variable being defined.
        :param math: The right-hand side of the rule as an ASTNode.
        """
        self._rate_rules.append(Rule(id=id_, label=label, var=var, math=math))

    def _add_reaction(self, id_: str, label: str) -> None:
        """Add a reaction variable to the template variables.

        :param id_: The variable ID.
        :param label: The variable description.
        """
        self._reactions.append(Reaction(index=len(self._reactions), id=id_, label=label))

    def _add_state_variable(
        self, id_: str, label: str, initial_value: Optional[float], units: str = NON_DIM_UNITS
    ) -> "StateVariable":
        """Add a state variable to the template variables.

        :param id_: The variable ID.
        :param label: The variable description.
        :param initial_value: The variable initial value.
        :param units: The variable units.
        :return: The created StateVariable record (callers read its derivative_id).
        """
        state_var = StateVariable(
            index=len(self._state_variables),
            id=id_,
            derivative_id=self._names.reserve(f"{DERIVATIVE_PREFIX}{id_}{DERIVATIVE_SUFFIX}"),
            label=label,
            initial_value=initial_value,
            units=units,
            name=self._names.sbml_name(id_),
        )

        self._state_variables.append(state_var)

        return state_var

    def _add_placeholder_state_variable(self) -> None:
        """Add a synthetic state variable with a zero derivative for a model with no ODEs.

        Chaste's ODE system integrates a vector of state variables, so a model with no continuous
        dynamics still needs at least one. The placeholder is constant (``dy/dt = 0``), starts at
        zero and is not an SBML output, so it never appears in the expected results; it exists only
        to give the solver a trivial state to advance in time while the model's real outputs are
        recomputed each step.
        """
        placeholder_id = self._names.reserve(PREFIX_SEP.join([CHASTE_PREFIX, PLACEHOLDER_STATE_ID]))
        state_var = self._add_state_variable(placeholder_id, "Placeholder state variable (model has no ODEs)", 0.0)
        # Initialise the member to zero
        self._add_equation(var=placeholder_id, math=parseL3Formula("0.0"), eq_type=EquationType.INITIAL_VALUE)
        self._add_equation(var=state_var.derivative_id, math=parseL3Formula("0.0"), eq_type=EquationType.DERIVATIVE)

    def _add_stoichiometry_variable(self, species_reference: "SpeciesReference") -> None:
        """Add a stoichiometry variable to the template variables.

        :param species_reference: The related SpeciesReference.
        """
        id_ = species_reference.getId()
        label = species_reference.getName().strip()

        # A named speciesReference's stoichiometry defaults to 1 when not explicitly set, so use
        # getStoichiometry() (which returns that default) as the initial value even when
        # isSetStoichiometry() is False - as long as it is a real number. Without this a defaulted
        # stoichiometry stays uninitialised (e.g. case 1800). An L3 variable stoichiometry left
        # unset reads as NaN; it is initialised by its rule or assignment instead.
        initial_value = None
        stoichiometry = species_reference.getStoichiometry()
        if not isnan(stoichiometry):
            initial_value = stoichiometry
            math = parseL3Formula(f"{initial_value}")
            self._add_equation(var=id_, math=math, eq_type=EquationType.INITIAL_VALUE)

        # A speciesReference's stoichiometry may be driven by a rate rule or an
        # assignment rule rather than a constant `stoichiometry` attribute. Note
        # these are not in _sbml_parameters, so _format_parameters does not handle them.
        rate_rules = {rule.var: rule.math for rule in self._rate_rules}
        assignment_rules = {rule.var: rule.math for rule in self._assignment_rules}
        if id_ in rate_rules:
            # A rate-rule stoichiometry is integrated as a state variable
            state_var = self._add_state_variable(id_, label, initial_value, NON_DIM_UNITS)  # TODO: Use correct units
            self._add_equation(var=state_var.derivative_id, math=rate_rules[id_], eq_type=EquationType.DERIVATIVE)
        elif id_ in assignment_rules:
            # An assignment-rule stoichiometry is a derived quantity recomputed each step.
            self._add_derived_quantity(id_, label, initial_value, NON_DIM_UNITS)  # TODO: Use correct units
            self._add_equation(var=id_, math=assignment_rules[id_], eq_type=EquationType.ASSIGNMENT_RULE)
        else:
            # A constant stoichiometry parameter.
            self._add_parameter(id_, label, initial_value, NON_DIM_UNITS)  # TODO: Use correct units

        if species_reference.isSetStoichiometryMath():
            math = species_reference.getStoichiometryMath().getMath()
            self._add_equation(var=id_, math=math, eq_type=EquationType.ASSIGNMENT_RULE)

    def _extract_odes(self) -> None:
        """Extract the ODEs equations for each species."""
        # Fast reactions require solving the fast subsystem to equilibrium as an algebraic
        # constraint (a DAE / quasi-steady-state problem) rather than as an ordinary ODE term.
        # This is not implemented, so reject such models instead of emitting incorrect ODEs.
        if any(r.isSetFast() and r.getFast() for r in self._sbml_reactions):
            raise NotImplementedError("Fast reactions are not supported.")

        # Note: rules must be processed before ODE extraction
        odes = {}

        def _update_ode(species_ref: "SpeciesReference", rxn: "SbmlReaction", add: bool) -> None:
            """Update the ODE for a species based on a species reference in a reaction.

            Each ODE will essentially be the sum of the products minus the
            sum of the reactants divided by the compartment volume

            :param species_ref: The species reference.
            :param rxn: The reaction.
            :param add: True if the species reference is a product, False if a reactant.
            """
            # A boundary-condition species' amount is not changed by reactions (its value is
            # held constant or set by a rule), so reactions contribute nothing to its ODE.
            species = self._sbml_model.getSpecies(species_ref.getSpecies())
            if species is not None and species.getBoundaryCondition():
                return

            rhs = rxn.getId()

            # Account for stoichiometry
            if species_ref.isSetId():
                self._add_stoichiometry_variable(species_ref)
                sto_id = species_ref.getId()
                rhs = f"({sto_id} * {rhs})"

            elif species_ref.isSetStoichiometryMath():
                sto_math = species_ref.getStoichiometryMath().getMath()
                sto_formula = self._formula_to_string(sto_math)
                rhs = f"({sto_formula} * {rhs})"

            elif species_ref.isSetStoichiometry():
                sto_value = species_ref.getStoichiometry()
                if sto_value != 1:
                    rhs = f"({sto_value} * {rhs})"

            # Update the ODE
            species_id = species_ref.getSpecies()
            if species_id in odes:
                if add:
                    odes[species_id] += f" + {rhs}"
                else:
                    odes[species_id] += f" - {rhs}"
            else:
                if add:
                    odes[species_id] = rhs
                else:
                    odes[species_id] = f"-{rhs}"

        # Extract ODEs from reactions
        for reaction in self._sbml_reactions:
            products = reaction.getListOfProducts()
            reactants = reaction.getListOfReactants()

            # Decompose reaction into sum of products minus sum of reactants
            for product in products:
                _update_ode(product, reaction, add=True)

            for reactant in reactants:
                _update_ode(reactant, reaction, add=False)

        # Convert string formulas to ASTs
        for var, rhs in odes.items():
            odes[var] = parseL3Formula(f"({rhs})")

        # Extract ODEs from rate rules. A rate rule fully determines its variable's derivative,
        # taking precedence over any reaction contribution.
        for rule in self._rate_rules:
            var = rule.var
            math = rule.math
            odes[var] = math

        # An empty odes dict is allowed: a model with no reactions or rate rules has no
        # continuous dynamics, and build() synthesises a placeholder state variable so the
        # generated ODE system still has something for the solver to integrate.
        self._odes = odes

    def _total_time_derivative(self, ast_node: "ASTNode") -> str:
        """Compute the total time derivative of an AST expression.

        Returns a formula string for d(expr)/dt computed via the chain rule:
          d(expr)/dt = sum_i (∂expr/∂s_i) * (ds_i/dt)
        over all state variables s_i with ODEs.  Returns an empty string when
        the expression has no time-varying dependencies on any state variable.
        """
        terms = []
        for var_id, ode_ast in self._odes.items():
            partial = ast_node.derivative(var_id)
            if partial is None:
                continue
            partial_str = formulaToL3String(partial)
            if partial_str == "0":
                continue
            ode_str = formulaToL3String(ode_ast)
            terms.append(f"({partial_str}) * ({ode_str})")
        return " + ".join(terms)

    def _format_compartments(self) -> None:
        """Add compartments to template variables."""
        # Note: rules must be processed before compartments
        assignment_rules = {rule.var: rule.math for rule in self._assignment_rules}
        initial_assignments = {ia.var: ia.math for ia in self._initial_assignments}

        for compartment in self._sbml_compartments:
            id_ = compartment.getId()
            label = compartment.getName().strip()

            size = get_compartment_size(compartment)
            math = parseL3Formula(f"{size}")
            self._add_equation(
                var=id_,
                math=math,
                eq_type=EquationType.INITIAL_VALUE,
            )
            units = compartment.getUnits() if compartment.isSetUnits() else NON_DIM_UNITS

            is_const = compartment.isSetConstant() and compartment.getConstant()

            if id_ in initial_assignments:
                math = initial_assignments[id_]
                self._add_equation(var=id_, math=math, eq_type=EquationType.INITIAL_ASSIGNMENT)

            if is_const:
                # Derived quantity
                self._add_derived_quantity(id_, label, size, units)
            elif id_ in self._odes:
                # State variable
                math = self._odes[id_]
                state_var = self._add_state_variable(id_, label, size, units)
                self._add_equation(
                    var=state_var.derivative_id,
                    math=math,
                    eq_type=EquationType.DERIVATIVE,
                )
            elif id_ in assignment_rules:
                # Derived quantity
                math = assignment_rules[id_]
                self._add_derived_quantity(id_, label, size, units)
                self._add_equation(var=id_, math=math, eq_type=EquationType.ASSIGNMENT_RULE)
            else:
                # Variable parameter
                self._add_parameter(id_, label, size, units)

    def _format_equations(self) -> None:
        """Convert and sort equations."""
        for eq in self._equations:
            if eq.math:
                eq.rhs = self._formula_to_string(eq.math, eq.local_parameters)

        self._sort_equations()

    def _format_events(self) -> None:
        """Add events to template variables."""
        for event in self._sbml_events:
            trigger = event.getTrigger()
            if trigger is None or trigger.getMath() is None:
                # An event with no trigger, or a trigger with no MathML, can never fire.
                continue
            self._reject_unsupported_delay(event)
            label = event.getName().strip()
            event_type = self._guess_event_type(label)
            # Compute the trigger formula first: it mutates the trigger AST in place (stripping
            # units, renaming the avogadro csymbol), and the distance below deep-copies that AST.
            trigger_formula = self._formula_to_string(trigger.getMath())
            trigger_distance = self._event_trigger_distance(event)
            assignments = self._event_assignments(event)
            self._compensate_compartment_resizes(event, assignments)
            # SBML trigger initialValue="false" means the trigger is treated as false just before
            # t=0, allowing the event to fire immediately if the condition is true at t=0.
            initial_satisfied = event.getTrigger().getInitialValue()
            priority = self._event_priority(event)
            self._add_event(
                label,
                trigger_formula,
                assignments,
                trigger_distance,
                event_type,
                initial_satisfied,
                priority,
            )

    def _reject_unsupported_packages(self) -> None:
        """Raise if the model uses an SBML package the generator cannot translate.

        :raises NotImplementedError: if the model uses the flux-balance-constraints (fbc) package.
        """
        if self._sbml_model.getPlugin("fbc") is not None:
            raise NotImplementedError("Flux balance constraint models (SBML fbc package) are not supported.")

    def _reject_unsupported_delay(self, event) -> None:
        """Raise if the event has a non-zero delay (delays are not supported)."""
        if not event.isSetDelay():
            return
        delay_math = event.getDelay().getMath()
        if delay_math is None:
            # A delay element with no MathML is treated as no delay.
            return
        math = self._formula_to_string(delay_math)
        try:
            delay = float(math)
        except ValueError:
            delay = 9999

        if delay != 0.0:
            # Delay of zero is equivalent to no delay
            raise NotImplementedError("Events with delays are not supported.")

    @staticmethod
    def _guess_event_type(label: str) -> EventType:
        """Guess the event type from its label (cell division vs unknown)."""
        cell_division_terms = ["cell division", "cytokinesis", "mitosis", "meiosis"]
        for term in cell_division_terms:
            if all(word in label.lower() for word in term.split()):
                return EventType.CELL_DIVISION
        return EventType.UNKNOWN

    def _event_trigger_distance(self, event) -> str:
        """Signed distance to the trigger boundary, for CVODE root-finding.

        Must run after the trigger formula: it deep-copies the (already unit-stripped) trigger
        AST and expands function-definition calls so the relational operator is visible. The
        distance is negative when the condition is false, zero at the crossing, positive when true.
        """
        distance_math = event.getTrigger().getMath().deepCopy()
        SBMLTransforms.replaceFD(distance_math, self._sbml_model.getListOfFunctionDefinitions())

        trigger_distance = "1.0"
        node_type = distance_math.getType()
        if node_type in [
            AST_RELATIONAL_LT,
            AST_RELATIONAL_GT,
            AST_RELATIONAL_EQ,
            AST_RELATIONAL_LEQ,
            AST_RELATIONAL_GEQ,
            AST_RELATIONAL_NEQ,
        ]:
            lc = self._formula_to_string(distance_math.getLeftChild())
            rc = self._formula_to_string(distance_math.getRightChild())

            # Distance is negative when the condition is false,
            # zero at the point where the condition switches from false to true,
            # and positive when the condition is true.
            if node_type == AST_RELATIONAL_GT:
                # gt(4.5    , 5.0    ) -> condition=false, dist=-0.5-eps
                # gt(5.0    , 5.0+eps) -> condition=false, dist=-eps-eps
                # gt(5.0    , 5.0    ) -> condition=false, dist=-eps
                # gt(5.0+eps, 5.0    ) -> condition=true, dist=0.0
                # gt(5.5    , 5.0    ) -> condition=true, dist=0.5-eps
                trigger_distance = f"({lc}) - ({rc}) - std::numeric_limits<double>::epsilon()"
            elif node_type == AST_RELATIONAL_GEQ:
                # geq(4.5    , 5.0    ) -> condition=false, dist=-0.5-eps
                # geq(5.0    , 5.0+eps) -> condition=false, dist=-eps-eps
                # geq(5.0    , 5.0    ) -> condition=true, dist=-eps
                # geq(5.0+eps, 5.0    ) -> condition=true, dist=0.0
                # geq(5.5    , 5.0    ) -> condition=true, dist=0.5-eps
                trigger_distance = f"({lc}) - ({rc}) - std::numeric_limits<double>::epsilon()"
            elif node_type == AST_RELATIONAL_LT:
                # lt(5.5    , 5.0    ) -> condition=false, dist=-0.5-eps
                # lt(5.0+eps, 5.0    ) -> condition=false, dist=-eps-eps
                # lt(5.0    , 5.0    ) -> condition=false, dist=-eps
                # lt(5.0    , 5.0+eps) -> condition=true, dist=0.0
                # lt(4.5    , 5.0    ) -> condition=true, dist=0.5-eps
                trigger_distance = f"({rc}) - ({lc}) - std::numeric_limits<double>::epsilon()"
            elif node_type == AST_RELATIONAL_LEQ:
                # leq(5.5    , 5.0    ) -> condition=false, dist=-0.5-eps
                # leq(5.0+eps, 5.0    ) -> condition=false, dist=-eps-eps
                # leq(5.0    , 5.0    ) -> condition=true, dist=-eps
                # leq(5.0    , 5.0+eps) -> condition=true, dist=0.0
                # leq(4.5    , 5.0    ) -> condition=true, dist=0.5-eps
                trigger_distance = f"({rc}) - ({lc}) - std::numeric_limits<double>::epsilon()"
            elif node_type == AST_RELATIONAL_EQ:
                # eq(4.5    , 5.0    ) -> condition=false, dist=-0.5
                # eq(5.0    , 5.0+eps) -> condition=false, dist=-eps
                # eq(5.0    , 5.0    ) -> condition=true, dist=0.0
                # eq(5.0+eps, 5.0    ) -> condition=false, dist=-eps
                # eq(5.5    , 5.0    ) -> condition=false, dist=-0.5
                trigger_distance = f"-std::abs(({lc}) - ({rc}))"
            else:  # AST_RELATIONAL_NEQ
                # neq(4.5    , 5.0    ) -> condition=true, dist=0.5-eps
                # neq(5.0    , 5.0+eps) -> condition=true, dist=0.0
                # neq(5.0    , 5.0    ) -> condition=false, dist=-eps
                # neq(5.0+eps, 5.0    ) -> condition=true, dist=0.0
                # neq(5.5    , 5.0    ) -> condition=true, dist=0.5-eps
                trigger_distance = f"std::abs(({lc}) - ({rc})) - std::numeric_limits<double>::epsilon()"

            # TODO: Distance calculation assumes two operands. Extend to more operands?
            # e.g. trigger: geq(3.0, 6.0, 7.0, 9.0) -> condition=false, dist=max(3.0, 1.0, 2.0)=3.0
            # TODO: Distance calculation assumes a single relational operator. Handle compound expressions?
            # e.g. trigger: and(geq(3.0, 6.0), geq(7.0, 9.0)) -> condition=false, dist=max(3.0, 2.0)=3.0
        return trigger_distance

    def _event_assignments(self, event) -> list["EventAssignment"]:
        """Build the list of assignment records the event applies when it fires."""
        assignments = []
        for assignment in event.getListOfEventAssignments():
            # An event assignment with no MathML assigns nothing, leaving its target
            # unchanged, so skip it rather than dereferencing a null math node.
            if assignment.getMath() is None:
                continue

            lhs = assignment.getVariable()
            type_ = self._get_variable_type(lhs)
            index = self._get_variable_index(lhs)
            rhs = self._formula_to_string(assignment.getMath())

            assignments.append(EventAssignment(index=index, lhs=lhs, rhs=rhs, type=type_))
        return assignments

    def _compensate_compartment_resizes(self, event, assignments) -> None:
        """Append/adjust assignments so a compartment resize conserves species amounts."""
        # An event that resizes a compartment conserves the amount of each concentration
        # species in it (SBML semantics), so the species' concentration must be rescaled by
        # old_size / new_size. The codegen tracks such species as concentration, so without
        # this their amount (concentration * size) would jump with the compartment. Append
        # compensating assignments; like all event assignments they are evaluated from the
        # pre-event state, so the compartment id still reads its old size here.
        compartment_ids = {compartment.getId() for compartment in self._sbml_compartments}
        assignment_by_lhs = {assignment.lhs: assignment for assignment in assignments}
        resized = {a.lhs: a.rhs for a in assignments if a.lhs in compartment_ids}

        # A compartment may also be resized indirectly: an assignment rule sets it from a
        # variable the event assigns (e.g. C = fakeC with the event assigning fakeC). Its new
        # size is the rule expression with the event's assignments substituted in.
        event_assignment_math = {
            ea.getVariable(): ea.getMath() for ea in event.getListOfEventAssignments() if ea.getMath() is not None
        }
        assignment_rule_math = {ar.var: ar.math for ar in self._assignment_rules}
        for compartment_id, rule_math in assignment_rule_math.items():
            if compartment_id not in compartment_ids or compartment_id in resized:
                continue
            referenced = set()
            collect_ast_names(rule_math, referenced)
            if referenced & set(event_assignment_math):
                substituted = substitute_ast_names(rule_math, event_assignment_math)
                resized[compartment_id] = self._formula_to_string(substituted)

        for compartment_id, new_size in resized.items():
            for species in self._sbml_species:
                species_id = species.getId()
                if species.getCompartment() != compartment_id:
                    continue
                has_only_substance_units = species.isSetHasOnlySubstanceUnits() and species.getHasOnlySubstanceUnits()
                species_type = self._get_variable_type(species_id)
                # Amount-tracked species need no rescale (their amount is unchanged, and the
                # concentration = amount / size follows automatically). Assignment-rule species
                # are recomputed each step, so they cannot be assigned here.
                if has_only_substance_units or species_type not in (VarType.STATE_VARIABLE, VarType.PARAMETER):
                    continue

                scale = f"{compartment_id} / ({new_size})"
                explicit = assignment_by_lhs.get(species_id)
                if explicit is not None:
                    # The event already assigns this species a concentration. That value is
                    # taken at the old compartment size, so scale it by old_size / new_size to
                    # carry the resulting amount across the simultaneous resize.
                    explicit.rhs = f"({explicit.rhs}) * {scale}"
                else:
                    # The species is otherwise unchanged by the event, so conserve its amount.
                    assignments.append(
                        EventAssignment(
                            index=self._get_variable_index(species_id),
                            lhs=species_id,
                            rhs=f"{species_id} * {scale}",
                            type=species_type,
                        )
                    )

    def _event_priority(self, event) -> Optional[str]:
        """Return the event priority expression, or None if it has none."""
        if event.isSetPriority() and event.getPriority().getMath() is not None:
            return self._formula_to_string(event.getPriority().getMath())
        return None

    def _format_function_definitions(self) -> None:
        """Add function definitions to template variables."""
        for fd in self._sbml_function_definitions:
            if fd.getBody() is None:
                # A function definition with no MathML body cannot be emitted.
                continue
            fd_id = fd.getId()
            label = fd.getName().strip()
            arg_list = get_function_definition_arguments(fd)
            args = ", ".join(map(lambda x: f"double {x}", arg_list))
            body = self._formula_to_string(fd.getBody())

            self._add_function(
                fd_id,
                label,
                args,
                body,
            )

    def _format_initial_assignments(self) -> None:
        """Add initial assignments to template variables."""
        sbml_param_ids = {p.getId() for p in self._sbml_parameters}
        for ia in self._sbml_initial_assignments:
            if ia.getMath() is not None:
                id_ = ia.getId()
                label = ia.getName().strip()
                var = ia.getSymbol()
                math = ia.getMath()
                self._add_initial_assignment(id_, label, var, math)

                # Species reference stoichiometry variables are not in _sbml_parameters, so
                # _format_parameters won't apply their initial assignments. Do it here.
                # A variable stoichiometry (driven by a rate rule) is a state variable; its
                # initial assignment overrides the speciesReference stoichiometry attribute,
                # and Initialise's SetDefaultInitialCondition then carries it into the ICs.
                if var not in sbml_param_ids:
                    param_ids = {p.id for p in self._parameters}
                    state_var_ids = {s.id for s in self._state_variables}
                    if var in param_ids or var in state_var_ids:
                        self._add_equation(var=var, math=math, eq_type=EquationType.INITIAL_ASSIGNMENT)

    def _format_parameters(self) -> None:
        """Add parameters to template variables."""
        # Note: rules must be processed before parameters
        assignment_rules = {rule.var: rule.math for rule in self._assignment_rules}
        initial_assignments = {ia.var: ia.math for ia in self._initial_assignments}

        for param in self._sbml_parameters:
            param_id = param.getId()
            label = param.getName().strip()

            value = None
            if param.isSetValue():
                value = param.getValue()
                self._add_equation(
                    var=param_id,
                    math=parseL3Formula(f"{value}"),
                    eq_type=EquationType.INITIAL_VALUE,
                )

            units = NON_DIM_UNITS

            is_const = param.isSetConstant() and param.getConstant()

            if param.isSetUnits():
                units = param.getUnits()

            if param_id in initial_assignments:
                math = initial_assignments[param_id]
                self._add_equation(var=param_id, math=math, eq_type=EquationType.INITIAL_ASSIGNMENT)

            if param_id in self._odes:
                # State variable
                math = self._odes[param_id]
                state_var = self._add_state_variable(param_id, label, value, units)
                d_param_id = state_var.derivative_id
                self._add_equation(var=d_param_id, math=math, eq_type=EquationType.DERIVATIVE)

            elif param_id in assignment_rules:
                # Derived quantity
                math = assignment_rules[param_id]
                self._add_derived_quantity(param_id, label, value, units)
                self._add_equation(var=param_id, math=math, eq_type=EquationType.ASSIGNMENT_RULE)

            else:
                # Fixed / variable parameter
                self._add_parameter(param_id, label, value, units, is_const=is_const)

    def _format_reactions(self) -> None:
        """Add reactions to template variables."""
        for reaction in self._sbml_reactions:
            id_ = reaction.getId()
            label = reaction.getName().strip()
            self._add_reaction(id_, label)

            # A reaction's ID denotes its rate (flux), an observable that rules, events and
            # test outputs can read. Expose it as a derived quantity. The flux member is
            # already declared and computed by the reaction machinery, so it keeps its
            # VarType.REACTION type and is not re-declared (kind=REACTION).
            self._add_derived_quantity(id_, label, None, NON_DIM_UNITS, kind=DerivedQuantityKind.REACTION)

            kinetic_law = reaction.getKineticLaw()
            if kinetic_law is None:
                species = [r.getSpecies() for r in reaction.getListOfReactants()]
                rhs = " * ".join(species)
                self._add_equation(
                    var=id_,
                    math=parseL3Formula(rhs),
                    eq_type=EquationType.REACTION,
                )

            else:
                math = kinetic_law.getMath()

                parameters = kinetic_law.getListOfParameters()
                local_parameters = []
                for param in parameters:
                    param_id = param.getId()
                    if not param.isSetValue():
                        raise ValueError(f"Local parameter {param_id} in reaction {id_} has no value.")
                    local_parameters.append(
                        LocalParameter(id=param_id, label=param.getName().strip(), value=str(param.getValue()))
                    )
                self._add_equation(
                    var=id_,
                    math=math,
                    local_parameters=local_parameters,
                    eq_type=EquationType.REACTION,
                )

    def _format_rules(self) -> None:
        """Add rules to template variables."""
        # Algebraic rules are not implemented
        if any(r.getTypeCode() == SBML_ALGEBRAIC_RULE for r in self._sbml_rules):
            raise NotImplementedError("Algebraic rules are not supported.")

        # Get assignment and rate rules
        for rule in self._sbml_rules:
            math = rule.getMath()
            if math is not None:
                id_ = rule.getId()
                label = rule.getName().strip()
                var = rule.getVariable()
                if rule.getTypeCode() == SBML_ASSIGNMENT_RULE:
                    self._add_assignment_rule(id_, label, var, math)
                elif rule.getTypeCode() == SBML_RATE_RULE:
                    self._add_rate_rule(id_, label, var, math)

    def _format_species(self) -> None:
        """Add species to template variables (rules must be processed first)."""
        assignment_rules = {ar.var: ar.math for ar in self._assignment_rules}
        rate_rules = {rr.var: rr.math for rr in self._rate_rules}
        initial_assignments = {ia.var: ia.math for ia in self._initial_assignments}

        for species in self._sbml_species:
            conversion_factor = self._species_conversion_factor(species)

            # If there's a compartment we'll normalise the ODEs, so declare it as non-dimensional
            compartment = self._sbml_compartments.get(species.getCompartment())
            units = NON_DIM_UNITS if compartment else species.getSubstanceUnits()
            has_only_substance_units = species.isSetHasOnlySubstanceUnits() and species.getHasOnlySubstanceUnits()

            initial_value = self._add_species_initial_value(
                species, has_only_substance_units, assignment_rules, initial_assignments
            )
            self._add_species_dynamics(
                species,
                units,
                has_only_substance_units,
                conversion_factor,
                initial_value,
                assignment_rules,
                rate_rules,
                initial_assignments,
            )

            if not has_only_substance_units:
                # Concentration species: add an "amount" derived quantity (amount = conc * volume)
                self._add_amount(species)
            else:
                # Amount species: add a "concentration" derived quantity (conc = amount / volume)
                self._add_concentration(species)

    def _species_conversion_factor(self, species) -> Optional[str]:
        """Return the species' conversion factor, falling back to the model's, else None."""
        if species.isSetConversionFactor():
            return species.getConversionFactor()
        if self._sbml_model.isSetConversionFactor():
            return self._sbml_model.getConversionFactor()
        return None

    def _compartment_time_derivative(self, compartment_id, assignment_rules) -> str:
        """Return the compartment's dC/dt as an infix string, or "" if it is constant.

        A species in concentration in a time-varying compartment needs a dilution term
        (ds/dt gains -s*(dC/dt)/C) to conserve its amount; this supplies that dC/dt from the
        compartment's own ODE or assignment rule.
        """
        if compartment_id in self._odes:
            return formulaToL3String(self._odes[compartment_id])
        if compartment_id in assignment_rules:
            return self._total_time_derivative(assignment_rules[compartment_id])
        return ""

    def _add_species_unit_conversion(self, species_id, description, rhs) -> None:
        """Add a synthetic initial assignment converting a species' initial value between units."""
        ia_id = PREFIX_SEP.join([CHASTE_PREFIX, INITIAL_ASSIGNMENT_PREFIX, species_id])
        ia_math = parseL3Formula(rhs)
        self._add_initial_assignment(ia_id, f"Convert {species_id} {description}", species_id, ia_math)
        self._add_equation(var=species_id, math=ia_math, eq_type=EquationType.INITIAL_ASSIGNMENT)

    def _add_species_initial_value(
        self, species, has_only_substance_units, assignment_rules, initial_assignments
    ) -> Optional[float]:
        """Emit the species' initial-value equation and any concentration<->amount conversion.

        :return: the initial concentration or amount, or None if neither is set.
        """
        species_id = species.getId()
        compartment_id = species.getCompartment()

        if species.isSetInitialConcentration():
            initial_value = species.getInitialConcentration()
            self._add_equation(
                var=species_id, math=parseL3Formula(f"{initial_value}"), eq_type=EquationType.INITIAL_VALUE
            )
            if (
                has_only_substance_units
                and (species_id not in assignment_rules)
                and (species_id not in initial_assignments)
            ):
                # Convert initial concentration to amount via a custom initial assignment
                self._add_species_unit_conversion(
                    species_id, "concentration to amount", f"{species_id} * {compartment_id}"
                )
            return initial_value

        if species.isSetInitialAmount():
            initial_value = species.getInitialAmount()
            self._add_equation(species_id, parseL3Formula(f"{initial_value}"), eq_type=EquationType.INITIAL_VALUE)
            if (
                (not has_only_substance_units)
                and (species_id not in assignment_rules)
                and (species_id not in initial_assignments)
            ):
                # Convert initial amount to concentration via a custom initial assignment
                self._add_species_unit_conversion(
                    species_id, "amount to concentration", f"{species_id} / {compartment_id}"
                )
            return initial_value

        return None

    def _add_species_dynamics(
        self,
        species,
        units,
        has_only_substance_units,
        conversion_factor,
        initial_value,
        assignment_rules,
        rate_rules,
        initial_assignments,
    ) -> None:
        """Classify a species (derived quantity, state variable or parameter) and emit its equation."""
        species_id = species.getId()
        label = species.getName().strip()
        compartment_id = species.getCompartment()

        is_bc = species.isSetBoundaryCondition() and species.getBoundaryCondition()

        if species_id in initial_assignments:
            math = initial_assignments[species_id]
            self._add_equation(var=species_id, math=math, eq_type=EquationType.INITIAL_ASSIGNMENT)

        if species_id in assignment_rules:
            # Derived quantity (includes boundary conditions with assignment rules)
            math = assignment_rules[species_id]
            self._add_derived_quantity(species_id, label, initial_value, units)
            self._add_equation(var=species_id, math=math, eq_type=EquationType.ASSIGNMENT_RULE)

        elif is_bc and (species_id not in rate_rules):
            compartment_ddt = (
                "" if has_only_substance_units else self._compartment_time_derivative(compartment_id, assignment_rules)
            )
            if compartment_ddt:
                # State variable: boundary condition in changing compartment
                rhs = f"(-{species_id} * ({compartment_ddt})) / {compartment_id}"
                if conversion_factor is not None:
                    rhs = f"({rhs}) * {conversion_factor}"
                state_var = self._add_state_variable(species_id, label, initial_value, units)
                self._add_equation(
                    var=state_var.derivative_id, math=parseL3Formula(rhs), eq_type=EquationType.DERIVATIVE
                )
            else:
                # Constant boundary species (no rule, non-time-varying compartment): its value
                # is fixed except when an event changes it. Model it as a (variable) parameter
                # rather than a derived quantity so the value is stored per step (time-resolved
                # in the recorded solution) and event assignments go through the deferred
                # parameter mechanism instead of mutating a member that is never recorded.
                self._add_parameter(species_id, label, initial_value, units)

        elif species_id in self._odes:
            # State variable
            rhs = formulaToL3String(self._odes[species_id])

            # Add parentheses if there are multiple terms
            if "+" in rhs or "-" in rhs[1:]:
                rhs = f"({rhs})"

            # The conversion factor relates reaction extent to the change in species amount,
            # so it multiplies the reaction flux only - not a rate rule, which gives the
            # species' rate of change directly. Apply it before the dilution term below: that
            # term conserves amount in a time-varying compartment and is independent of the
            # conversion factor.
            if conversion_factor is not None and species_id not in rate_rules:
                rhs = f"({rhs}) * {conversion_factor}"

            # Add compartment scaling if defined by a reaction
            if not (has_only_substance_units or species_id in rate_rules):
                compartment_ddt = self._compartment_time_derivative(compartment_id, assignment_rules)
                if compartment_ddt:
                    rhs = f"({rhs} - {species_id} * ({compartment_ddt}))"
                # Scale by compartment volume
                rhs = f"{rhs} / {compartment_id}"

            # TODO: Handle time scaling
            # time_multiplier = self._get_timescale_multiplier()
            # if time_multiplier != 1.0:
            #     f"rDY[{state_variable_index}] *= {time_multiplier};"

            state_var = self._add_state_variable(species_id, label, initial_value, units)
            self._add_equation(var=state_var.derivative_id, math=parseL3Formula(rhs), eq_type=EquationType.DERIVATIVE)

        else:
            # Check whether the compartment is time-varying: a species in concentration with no
            # reactions or rules still needs a dilution ODE (ds/dt = -s*(dC/dt)/C) to conserve
            # its amount.
            compartment_ddt = (
                "" if has_only_substance_units else self._compartment_time_derivative(compartment_id, assignment_rules)
            )
            if compartment_ddt:
                rhs = f"(-{species_id} * ({compartment_ddt})) / {compartment_id}"
                if conversion_factor is not None:
                    rhs = f"({rhs}) * {conversion_factor}"
                state_var = self._add_state_variable(species_id, label, initial_value, units)
                self._add_equation(
                    var=state_var.derivative_id, math=parseL3Formula(rhs), eq_type=EquationType.DERIVATIVE
                )
            elif species_id in self._event_assigned_ids:
                # Otherwise constant but modified by an event: model as a (variable) parameter
                self._add_parameter(species_id, label, initial_value, units)
            else:
                # Truly constant: no reactions, no rules, no events, constant compartment
                self._add_derived_quantity(species_id, label, initial_value, units)

    def _formula_to_string(self, math: "ASTNode", local_parameters: Optional[list["LocalParameter"]] = None) -> str:
        """Convert an AST math formula to an equivalent C++ string.

        Thin wrapper over :func:`chaste_sbml._expressions.formula_to_string`, supplying the
        model's variable-type map and state variables.

        :param math: The AST math formula.
        :param local_parameters: Local parameters in scope for the formula.
        :return: The equivalent C++ string.
        """
        return formula_to_string(math, self._variable_type_map(), self._state_variables, local_parameters)

    def _get_timescale_multiplier(self) -> float:
        """Get the timescale multiplier.

        SBML uses seconds by default and Chaste uses hours.

        :return: The timescale multiplier.
        """
        for unit_def in self._sbml_unit_definitions:
            u_id = unit_def.getId()
            if u_id.lower() == "time":  # Do people ever call this something different?
                timescale = unit_def.getName().strip().lower()
                if "minute" in timescale:
                    return 60.0
                elif "hour" in timescale:
                    return 1.0
        return 3600.0

    def _get_variable_index(self, id_: str) -> int:
        """Get the index of a variable within its own collection.

        :param id_: The variable ID.
        :return: The variable index.
        :raises ValueError: if ``id_`` is not a recognised variable.
        """
        if self._index_of is None:
            self._index_of = self._build_variable_index()
        if id_ not in self._index_of:
            raise ValueError(f"ID '{id_}' is not a recognized variable.")
        return self._index_of[id_]

    def _build_variable_index(self) -> dict[str, int]:
        """Map each variable id to its index within its own collection.

        Built once and cached; only valid once the collections are complete, which holds because
        _get_variable_index is only called while formatting events, after every variable has been
        added. Reaction fluxes exposed as derived quantities are skipped so a reaction's index is
        its position among reactions, matching how it is declared.
        """
        index_of: dict[str, int] = {}
        for collection in (self._parameters, self._reactions, self._state_variables, self._functions):
            for record in collection:
                index_of[record.id] = record.index
        for derived_quantity in self._derived_quantities:
            if derived_quantity.kind != DerivedQuantityKind.REACTION:
                index_of[derived_quantity.id] = derived_quantity.index
        return index_of

    def _variable_type_map(self) -> dict[str, VarType]:
        """Map each variable id to its VarType, derived from the built record collections.

        The type of a variable is implied by which collection its record belongs to, so no separate
        mapping is kept. A reaction flux is both a Reaction (typed VarType.REACTION) and an
        REACTION-kind DerivedQuantity; the latter is skipped so the reaction id keeps its type.
        """
        types: dict[str, VarType] = {}
        for derived_quantity in self._derived_quantities:
            if derived_quantity.kind != DerivedQuantityKind.REACTION:
                types[derived_quantity.id] = VarType.DERIVED_QUANTITY
        for function in self._functions:
            types[function.id] = VarType.FUNCTION
        for parameter in self._parameters:
            types[parameter.id] = VarType.PARAMETER
        for reaction in self._reactions:
            types[reaction.id] = VarType.REACTION
        for state_variable in self._state_variables:
            types[state_variable.id] = VarType.STATE_VARIABLE
        return types

    def _get_variable_type(self, var_id: str) -> VarType:
        """Get the type of a variable based on its ID.

        :param var_id: The variable ID.
        :return: The variable's VarType, or VarType.UNKNOWN if it is not a known variable.
        """
        return self._variable_type_map().get(var_id, VarType.UNKNOWN)

    def _order_derived_quantities(self) -> None:
        """Order derived quantities: normal quantities, then reactions, then conversions.

        Reaction fluxes and the amount/concentration conversions (amt__/conc__) are added
        while processing reactions and species, interleaved with the other derived quantities.
        Grouping them after the model-intrinsic quantities keeps the latter at stable,
        contiguous indices. The sort is stable, so the relative order within each group is
        preserved; indices are then renumbered to match.
        """

        def group(dq: "DerivedQuantity") -> int:
            if dq.kind == DerivedQuantityKind.CONVERSION:
                return 2
            if dq.kind == DerivedQuantityKind.REACTION:
                return 1
            return 0

        self._derived_quantities.sort(key=group)
        for index, dq in enumerate(self._derived_quantities):
            dq.index = index

    def build(self) -> None:
        """Process the SBML model to set up the formatted variables for templates."""
        self._reject_unsupported_packages()

        self._index_of = None  # id -> index cache, built lazily by _get_variable_index
        self._odes = {}

        self._assignment_rules = []
        self._rate_rules = []
        self._initial_assignments = []

        self._state_variables = []
        self._derived_quantities = []
        self._parameters = []
        self._stoichiometry_variables = []
        self._equations = []

        self._reactions = []
        self._events = []
        self._functions = []

        # Ids assigned by some event, so species classification can model an otherwise-constant
        # event-modified species as a (variable) parameter.
        self._event_assigned_ids = {
            ea.getVariable() for event in self._sbml_events for ea in event.getListOfEventAssignments()
        }

        # Names already claimed by real SBML entities and the Chaste base classes. Synthetic
        # identifiers (derivatives, amount/concentration conversions, initial-assignment
        # intermediates) are allocated against this so they never collide.
        self._names.reset()

        # TODO: enforce processing order e.g. rules must be processed first
        self._format_rules()
        self._extract_odes()

        self._format_initial_assignments()

        self._format_compartments()
        self._format_species()
        self._format_parameters()

        if not self._state_variables:
            # The model has no continuous dynamics (no reactions, ODEs or rate rules). Add a
            # placeholder state variable with a zero derivative so the generated ODE system still
            # has something for the solver to integrate; its outputs (constants, assignment rules
            # of time, event-driven changes) are recomputed each step like any other model.
            self._add_placeholder_state_variable()

        self._format_reactions()
        self._format_function_definitions()
        self._format_events()

        self._format_equations()
        self._order_derived_quantities()

        self._check_name_conflicts()

    def _check_name_conflicts(self) -> None:
        """Fail if any generated C++ identifier clashes with another or a reserved name.

        Detect conflicts and raise rather than emit silently incorrect
        C++. Gathers every identifier the templates turn into a C++ name -- parameters, state
        variables and their derivatives, derived quantities (including amount/concentration
        conversions), stoichiometry variables, reactions and model functions -- and checks them
        for duplicates, C++ keywords, reserved Chaste base-class names and invalid identifiers.
        Reaction flux outputs are excluded from the derived quantities here as they are the same
        entities already counted under reactions.

        Kinetic-law local parameters (emitted as ``double <id> = ...;`` inside a reaction block)
        are checked too, against C++ keywords and reserved Chaste names but not for duplication:
        they are scoped per reaction, so an id that repeats across reactions or shadows a global is
        legal, but one equal to a reserved name -- e.g. ``time``, the emitted simulation-time
        parameter -- would shadow it and silently corrupt the reaction, so it is flagged.

        :raises NameConflictError: if any conflict is found.
        """
        identifiers = []
        identifiers += [(p.id, "parameter") for p in self._parameters]
        for var in self._state_variables:
            identifiers.append((var.id, "state variable"))
            identifiers.append((var.derivative_id, "state-variable derivative"))
        for dq in self._derived_quantities:
            if dq.kind == DerivedQuantityKind.REACTION:
                continue  # Declared under reactions; counted there.
            kind_label = (
                "amount/concentration conversion" if dq.kind == DerivedQuantityKind.CONVERSION else "derived quantity"
            )
            identifiers.append((dq.id, kind_label))
        identifiers += [(s.id, "stoichiometry variable") for s in self._stoichiometry_variables]
        identifiers += [(r.id, "reaction") for r in self._reactions]
        identifiers += [(f.id, "function") for f in self._functions]
        # Initial assignments are deliberately excluded: an initial assignment's id is the
        # existing variable it assigns to (emitted as `var = <expr>`), not a new declaration,
        # so counting it would double-count that variable.

        conflicts = find_name_conflicts(identifiers)
        for eq in self._equations:
            for lp in eq.local_parameters or ():
                conflicts.extend(find_name_conflicts([(lp.id, "local parameter")], reserved=CHASTE_RESERVED_NAMES))
        conflicts = sorted(conflicts)
        if conflicts:
            raise NameConflictError(
                f"Cannot generate '{self._model_name}': C++ identifier conflicts detected:\n  - "
                + "\n  - ".join(conflicts)
            )

    def _sort_equations(self) -> None:
        """Sort equations, ordering by groups, then by dependencies.

        Group order (ignored if there are dependencies that require a different order):
        1. Initial value equations e.g. `S0 = 1.0`.
        2. Initial assignment equations e.g. `S1 = 2 * S0`.
        3. Assignment rules e.g. `k0 = k0 + 1.0`.
        4. Reactions e.g. `J0 = k0 * S0 + k1 * S1`.
        5. Derivatives e.g. `d_S0_dt = J0 - J1`.
        6. Amounts e.g. `amt_S0 = S0 * compartment`.

        Dependency order: if equation A depends on B (i.e. A -> B), then B should
        come before A. It is assumed that dependencies are acyclic. This cannot
        sort cyclic dependencies such as A -> B -> A or A -> B -> C -> A.
        """

        def _depends(eq: "Equation", var: str) -> bool:
            """Check if a variable appears in the rhs of an equation.

            :param eq: The equation to check.
            :param var: The variable to check for.
            :return: True if the variable appears in the rhs of the equation, False otherwise.
            """
            # False if there is a local parameter with the same name as the variable.
            if eq.local_parameters:
                if any(var == param.id for param in eq.local_parameters):
                    return False

            return bool(re.search(rf"\b{var}\b", eq.rhs))

        # Sort equations by group
        group_order = [
            EquationType.INITIAL_VALUE,
            EquationType.INITIAL_ASSIGNMENT,
            EquationType.ASSIGNMENT_RULE,
            EquationType.REACTION,
            EquationType.DERIVATIVE,
            EquationType.CONVERSION,
            EquationType.UNKNOWN,
        ]

        sorted_equations = [eq for eq_type in group_order for eq in self._equations if eq.type == eq_type]

        # Sort equations by dependencies
        for _ in range(len(sorted_equations)):  # Max iterations for worst case
            re_sorted_equations = []
            for eq in sorted_equations:
                for i, other_eq in enumerate(re_sorted_equations):
                    # Insert eq before the first other equation that depends on it.
                    if _depends(other_eq, eq.var):
                        re_sorted_equations.insert(i, eq)
                        break
                else:
                    # No other equations depend on this one.
                    re_sorted_equations.append(eq)

            if sorted_equations == re_sorted_equations:
                # No changes were made.
                break

            sorted_equations = re_sorted_equations

        self._equations = sorted_equations
