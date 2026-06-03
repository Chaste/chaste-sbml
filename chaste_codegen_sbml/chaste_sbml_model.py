"""Module for the ChasteSbmlModel class."""

import abc
import os
import pathlib
import re
import shutil
import subprocess
import sys
from typing import TYPE_CHECKING, Optional

from jinja2 import Environment, PackageLoader, select_autoescape
from libsbml import (
    AST_FUNCTION_DELAY,
    AST_RELATIONAL_EQ,
    AST_RELATIONAL_GEQ,
    AST_RELATIONAL_GT,
    AST_RELATIONAL_LEQ,
    AST_RELATIONAL_LT,
    AST_RELATIONAL_NEQ,
    LIBSBML_OPERATION_SUCCESS,
    SBML_ALGEBRAIC_RULE,
    SBML_ASSIGNMENT_RULE,
    SBML_RATE_RULE,
    ConversionProperties,
    SBMLReader,
    formulaToL3String,
    parseL3Formula,
)

from ._config import (
    AMOUNT_PREFIX,
    CHASTE_PREFIX,
    DERIVATIVE_PREFIX,
    DERIVATIVE_SUFFIX,
    INITIAL_ASSIGNMENT_PREFIX,
    NON_DIM_UNITS,
    PREFIX_SEP,
    ROOT_DIR,
    EquationType,
    EventType,
    ModelType,
    VarType,
)
from ._utils import (
    generate_header_guard,
    get_compartment_size,
    get_function_definition_arguments,
    search_ast_type,
    to_camel_case,
    to_cpp_name,
)

if TYPE_CHECKING:
    from typing import Any

    from jinja2.environment import Template
    from libsbml import ASTNode, Reaction, Species, SpeciesReference


class ChasteSbmlModel:
    """Holds information about an SBML model for which Chaste code is to be generated."""

    __metaclass__ = abc.ABCMeta

    _jinja_env = Environment(
        loader=PackageLoader("chaste_codegen_sbml"),
        autoescape=select_autoescape(),
        trim_blocks=True,
        lstrip_blocks=True,
    )
    _jinja_env.globals["AMOUNT_PREFIX"] = AMOUNT_PREFIX
    _jinja_env.globals["EquationType"] = EquationType
    _jinja_env.globals["EventType"] = EventType
    _jinja_env.globals["PREFIX_SEP"] = PREFIX_SEP
    _jinja_env.globals["VarType"] = VarType

    # -- PUBLIC --------------------------------------

    def __init__(self, sbml_file: str, model_name: str = "", model_type: ModelType = ModelType.GENERIC) -> None:
        """Initialise the ChasteSbmlModel.

        :param sbml: The SBML file.
        :param name: The model name.
        :param model_type: The model type e.g. ModelType.SRN.
        """
        self._sbml_file = os.path.abspath(sbml_file)
        if not os.path.isfile(self._sbml_file):
            raise FileNotFoundError(f"Could not find SBML file: {self._sbml_file}")

        if model_name:
            self._model_name = model_name
        else:
            filename = os.path.splitext(os.path.basename(self._sbml_file))[0]
            model_name = to_camel_case(to_cpp_name(filename)) + "Sbml"
            self._model_name = model_name[0].upper() + model_name[1:]

        self._model_type = model_type

        self._ode_class_name = self._model_name + "OdeSystem"
        self._ode_hpp_filename = f"{self._ode_class_name}.hpp"
        self._ode_cpp_filename = f"{self._ode_class_name}.cpp"

        self._srn_class_name = ""
        self._srn_hpp_filename = ""
        self._srn_cpp_filename = ""
        if self._model_type == ModelType.SRN:
            self._srn_class_name = self._model_name + "SrnModel"
            self._srn_hpp_filename = f"{self._srn_class_name}.hpp"
            self._srn_cpp_filename = f"{self._srn_class_name}.cpp"

        self._cell_cycle_class_name = ""
        self._cell_cycle_hpp_filename = ""
        self._cell_cycle_cpp_filename = ""
        if self._model_type == ModelType.CELL_CYCLE:
            self._cell_cycle_class_name = self._model_name + "CellCycleModel"
            self._cell_cycle_hpp_filename = f"{self._cell_cycle_class_name}.hpp"
            self._cell_cycle_cpp_filename = f"{self._cell_cycle_class_name}.cpp"

        # Read the SBML model
        reader = SBMLReader()
        doc = reader.readSBMLFromFile(self._sbml_file)

        if doc.getNumErrors() > 0:
            doc.printErrors()
            raise ValueError(f"Errors found while reading SBML file: {self._sbml_file}")

        # Run required conversions
        config = ConversionProperties()
        # Sort assignment rules in order of dependence.
        config.addOption("sortRules")
        config.addOption("removeUnusedUnits")
        # Convert initial assignments to initial values where possible
        config.addOption("expandInitialAssignments")
        # config.addOption('replaceReactions')
        # config.addOption('expandFunctionDefinitions')

        status = doc.convert(config)
        if status != LIBSBML_OPERATION_SUCCESS:
            doc.printErrors()
            raise ValueError("Errors during conversion")

        self._sbml_model = doc.getModel()
        self._sbml_compartments = self._sbml_model.getListOfCompartments()
        self._sbml_events = self._sbml_model.getListOfEvents()
        self._sbml_function_definitions = self._sbml_model.getListOfFunctionDefinitions()
        self._sbml_initial_assignments = self._sbml_model.getListOfInitialAssignments()
        self._sbml_parameters = self._sbml_model.getListOfParameters()
        self._sbml_reactions = self._sbml_model.getListOfReactions()
        self._sbml_rules = self._sbml_model.getListOfRules()
        self._sbml_species = self._sbml_model.getListOfSpecies()
        self._sbml_unit_definitions = self._sbml_model.getListOfUnitDefinitions()

        self._variable_types = {}  # { id: VarType }
        self._odes = {}  # { id: str }

        self._assignment_rules = []  # [ { id: str, label: str, ... } ]
        self._rate_rules = []  # [ { id: str, label: str, ... } ]
        self._state_variables = []  # [ { id: str, label: str, ... } ]
        self._derived_quantities = []  # [ { id: str, label: str, ... } ]
        self._amounts = []  # [ { id: str, label: str, ... } ]
        self._parameters = []  # [ { id: str, label: str, ... } ]
        self._initial_assignments = []  # [ { id: str, label: str, ... } ]
        self._reactions = []  # [ { id: str, label: str, ... } ]
        self._events = []  # [ { name: str, trigger: str, ... } ]
        self._functions = []  # [ { name: str, args: [str], body: str } ]
        self._stoichiometry_variables = []  # [ { id: str, label: str, ... } ]
        self._equations = []  # [ { lhs: str, rhs: str, label: str } ]

        self._template_vars = {}  # type: dict[str, Any]

        self._outputs = {}  # { filename: code }

        self._process_model()

    @property
    def outputs(self) -> dict[str, str]:
        """Get the generated code outputs.

        :return: A dictionary of filename and code pairs.
        """
        return self._outputs

    def write(self, output_directory=None):
        """Generate Chaste code and write to file.

        :param output_directory: The output directory. Defaults to the current directory.
        """
        # Generate the code
        self._generate_outputs()

        # Write the code to file
        if output_directory:
            root_dir = pathlib.Path(output_directory)
        else:
            root_dir = pathlib.Path().cwd()

        for filename, code in self._outputs.items():
            file_path = root_dir / filename
            with open(file_path, "w") as f:
                f.write(code)

        # Format with clang-format — search the current interpreter's bin/ first
        # so the clang-format installed as a Python dependency is always found.
        python_bin = os.path.dirname(sys.executable)
        search_path = os.pathsep.join([python_bin, os.environ.get("PATH", "")])
        clang_format = shutil.which("clang-format", path=search_path)
        if clang_format is not None:
            for filename in self._outputs:
                file_path = str(root_dir / filename)
                subprocess.run(
                    [
                        clang_format,
                        "-i",
                        f"-style=file:{ROOT_DIR}/.clang-format",
                        str(file_path),
                    ],
                    check=True,
                )

    # -- PRIVATE ---------------------------------------

    def _add_amount(self, species: "Species") -> None:
        """Add an amount derived quantity variable to the template variables.

        :param species: The related Species.
        """
        species_id = species.getId()
        compartment_id = species.getCompartment()

        amt_id = AMOUNT_PREFIX + PREFIX_SEP + species_id
        amt_label = f"Amount of {species_id}"
        amt_units = NON_DIM_UNITS  # TODO: Use correct units
        amt_rhs = f"{species_id} * {compartment_id}"
        amt_math = parseL3Formula(amt_rhs)

        self._add_derived_quantity(amt_id, amt_label, None, amt_units, is_amount=True)
        self._add_equation(var=amt_id, math=amt_math, eq_type=EquationType.AMOUNT)

    def _add_assignment_rule(self, id_: str, label: str, var: str, math: "ASTNode") -> None:
        """Add an assignment rule to the template variables.

        :param id_: The rule ID.
        :param label: The rule description.
        :param var: The variable being defined.
        :param math: The right-hand side of the rule as an ASTNode.
        """
        self._assignment_rules.append(
            {
                "id": id_,
                "label": label,
                "var": var,
                "math": math,
            }
        )

    def _add_derived_quantity(
        self,
        id_: str,
        label: str,
        initial_value: Optional[float],
        units: str = NON_DIM_UNITS,
        is_amount: bool = False,
    ) -> None:
        """Add a derived quantity to the template variables.

        :param id_: The variable ID.
        :param label: The variable description.
        :param initial_value: The variable initial value.
        :param units: The variable units.
        :param is_amount: True if the derived quantity is an amount conversion.
        """
        self._derived_quantities.append(
            {
                "id": id_,
                "label": label,
                "index": len(self._derived_quantities),
                "initial_value": initial_value,
                "is_amount": is_amount,
                "units": units,
            }
        )
        self._variable_types[id_] = VarType.DERIVED_QUANTITY

    def _add_equation(
        self,
        var: str,
        math: "ASTNode",
        eq_type: EquationType = EquationType.UNKNOWN,
        local_parameters: Optional[list[dict[str, str]]] = None,
    ) -> None:
        """Add an equation to the template variables.

        These are expected to be equations of the form `var = rhs`.

        :param var: The variable being defined.
        :param math: The right-hand side of the equation as an AST.
        :param eq_type: The type of the equation.
        :param local_parameters: A list of local parameter names used in the equation.
        """
        eq = {
            "var": var,
            "rhs": "",
            "math": math,
            "type": eq_type,
            "local_parameters": local_parameters,
        }
        self._equations.append(eq)

    def _add_event(
        self,
        label: str,
        trigger: str,
        assignments: list[dict[str, "Any"]],
        distance: str,
        event_type: EventType,
    ) -> None:
        """Add an event to the template variables.

        :param label: The event description.
        :param trigger: The event trigger formula.
        :param assignments: The event assignments.
        :param distance: The distance for the event trigger.
        :param event_type: The type of the event (e.g., cell division).
        """
        self._events.append(
            {
                "label": label,
                "index": len(self._events),
                "trigger": trigger,
                "assignments": assignments,
                "distance": distance,
                "type": event_type,
            }
        )

    def _add_function(self, id_: str, label: str, args: str, body: str) -> None:
        """Add a function to the template variables.

        :param id_: The function ID.
        :param label: The function description.
        :param args: The function arguments.
        :param body: The function body.
        """
        self._functions.append(
            {
                "id": id_,
                "label": label,
                "index": len(self._functions),
                "args": args,
                "body": body,
            }
        )
        self._variable_types[id_] = VarType.FUNCTION

    def _add_initial_assignment(self, id_: str, label: str, var: str, math: Optional["ASTNode"] = None) -> None:
        """Add an initial assignment to the template variables.

        :param id_: The rule ID.
        :param label: The rule description.
        :param var: The variable being assigned.
        :param math: The right-hand side of the rule as an ASTNode.
        """
        self._initial_assignments.append({"id": id_, "label": label, "var": var, "math": math})

    def _add_output(self, filename: str, code: str) -> None:
        """Add generated code to the outputs dictionary.

        :param filename: The filename.
        :param code: The code.
        """
        self._outputs[filename] = code

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
            {
                "index": len(self._parameters),
                "id": id_,
                "is_const": is_const,
                "label": label,
                "initial_value": initial_value,
                "units": units,
            }
        )
        self._variable_types[id_] = VarType.PARAMETER

    def _add_rate_rule(self, id_: str, label: str, var: str, math: "ASTNode") -> None:
        """Add a rate rule to the template variables.

        :param id_: The rule ID.
        :param label: The rule description.
        :param var: The variable being defined.
        :param math: The right-hand side of the rule as an ASTNode.
        """
        self._rate_rules.append(
            {
                "id": id_,
                "label": label,
                "var": var,
                "math": math,
            }
        )

    def _add_reaction(self, id_: str, label: str) -> None:
        """Add a reaction variable to the template variables.

        :param id_: The variable ID.
        :param label: The variable description.
        """
        self._reactions.append(
            {
                "index": len(self._reactions),
                "id": id_,
                "label": label,
            }
        )
        self._variable_types[id_] = VarType.REACTION

    def _add_state_variable(
        self, id_: str, label: str, initial_value: Optional[float], units: str = NON_DIM_UNITS
    ) -> dict:
        """Add a state variable to the template variables.

        :param id_: The variable ID.
        :param label: The variable description.
        :param initial_value: The variable initial value.
        :param units: The variable units.
        """
        # TODO: Check for name clashes with derivative_id's
        state_var = {
            "index": len(self._state_variables),
            "id": id_,
            "derivative_id": f"{DERIVATIVE_PREFIX}{id_}{DERIVATIVE_SUFFIX}",
            "label": label,
            "initial_value": initial_value,
            "units": units,
        }

        self._state_variables.append(state_var)
        self._variable_types[id_] = VarType.STATE_VARIABLE

        return state_var

    def _add_stoichiometry_variable(self, species_reference: "SpeciesReference") -> None:
        """Add a stoichiometry variable to the template variables.

        :param species_reference: The related SpeciesReference.
        """
        id_ = species_reference.getId()
        label = species_reference.getName().strip()

        initial_value = None
        if species_reference.isSetStoichiometry():
            initial_value = species_reference.getStoichiometry()
            math = parseL3Formula(f"{initial_value}")
            self._add_equation(var=id_, math=math, eq_type=EquationType.INITIAL_VALUE)

        self._add_parameter(id_, label, initial_value, NON_DIM_UNITS)  # TODO: Use correct units

        if species_reference.isSetStoichiometryMath():
            math = species_reference.getStoichiometryMath().getMath()
            self._add_equation(var=id_, math=math, eq_type=EquationType.ASSIGNMENT_RULE)

    def _extract_odes(self) -> None:
        """Extract the ODEs equations for each species."""
        # Note: rules must be processed before ODE extraction
        odes = {}

        def _update_ode(species_ref: "SpeciesReference", rxn: "Reaction", add: bool) -> None:
            """Update the ODE for a species based on a species reference in a reaction.

            Each ODE will essentially be the sum of the products minus the
            sum of the reactants divided by the compartment volume

            :param species_ref: The species reference.
            :param rxn: The reaction.
            :param add: True if the species reference is a product, False if a reactant.
            """
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

        # Extract ODEs from rate rules
        for rule in self._rate_rules:
            var = rule["var"]
            math = rule["math"]
            if var in odes:
                raise ValueError(f"{var} has more than one rate rule and/or reaction.")
            odes[var] = math

        if not odes:
            raise NotImplementedError("Models without ODEs are not supported.")

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
        assignment_rules = {rule["var"]: rule["math"] for rule in self._assignment_rules}
        initial_assignments = {ia["var"]: ia["math"] for ia in self._initial_assignments}

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
                    var=state_var["derivative_id"],
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
            if eq["math"]:
                eq["rhs"] = self._formula_to_string(eq["math"])

        self._sort_equations()

    def _format_events(self) -> None:
        """Add events to template variables."""
        # TODO: Add priority

        for event in self._sbml_events:
            if event.isSetDelay():
                math = self._formula_to_string(event.getDelay().getMath())
                try:
                    delay = float(math)
                except ValueError:
                    delay = 9999

                if delay != 0.0:
                    # Delay of zero is equivalent to no delay
                    raise NotImplementedError("Events with delays are not supported.")

            label = event.getName().strip()

            # Try to guess the event type
            event_type = EventType.UNKNOWN

            cell_division_terms = ["cell division", "cytokinesis", "mitosis", "meiosis"]
            for term in cell_division_terms:
                if all(word in label.lower() for word in term.split()):
                    event_type = EventType.CELL_DIVISION
                    break

            trigger_math = event.getTrigger().getMath()
            trigger_formula = self._formula_to_string(trigger_math)

            trigger_distance = "1.0"
            node_type = trigger_math.getType()
            if node_type in [
                AST_RELATIONAL_LT,
                AST_RELATIONAL_GT,
                AST_RELATIONAL_EQ,
                AST_RELATIONAL_LEQ,
                AST_RELATIONAL_GEQ,
                AST_RELATIONAL_NEQ,
            ]:
                lc = self._formula_to_string(trigger_math.getLeftChild())
                rc = self._formula_to_string(trigger_math.getRightChild())

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
                    # geq(4.5    , 5.0    ) -> condition=false, dist=-0.5
                    # geq(5.0    , 5.0+eps) -> condition=false, dist=-eps
                    # geq(5.0    , 5.0    ) -> condition=true, dist=0.0
                    # geq(5.0+eps, 5.0    ) -> condition=true, dist=eps
                    # geq(5.5    , 5.0    ) -> condition=true, dist=0.5
                    trigger_distance = f"({lc}) - ({rc})"
                elif node_type == AST_RELATIONAL_LT:
                    # lt(5.5    , 5.0    ) -> condition=false, dist=-0.5-eps
                    # lt(5.0+eps, 5.0    ) -> condition=false, dist=-eps-eps
                    # lt(5.0    , 5.0    ) -> condition=false, dist=-eps
                    # lt(5.0    , 5.0+eps) -> condition=true, dist=0.0
                    # lt(4.5    , 5.0    ) -> condition=true, dist=0.5-eps
                    trigger_distance = f"({rc}) - ({lc}) - std::numeric_limits<double>::epsilon()"
                elif node_type == AST_RELATIONAL_LEQ:
                    # leq(5.5    , 5.0    ) -> condition=false, dist=-0.5
                    # leq(5.0+eps, 5.0    ) -> condition=false, dist=-eps
                    # leq(5.0    , 5.0    ) -> condition=true, dist=0.0
                    # leq(5.0    , 5.0+eps) -> condition=true, dist=eps
                    # leq(4.5    , 5.0    ) -> condition=true, dist=0.5
                    trigger_distance = f"({rc}) - ({lc})"
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
                # e.g. trigger: geq(3.0, 6.0, 7.0, 9.0) -> condition=false, dist=min(3.0, 1.0, 2.0)=1.0

            assignments = []
            for assignment in event.getListOfEventAssignments():
                lhs = assignment.getVariable()
                type_ = self._get_variable_type(lhs)
                index = self._get_variable_index(lhs)
                rhs = self._formula_to_string(assignment.getMath())

                assignments.append(
                    {
                        "index": index,
                        "lhs": lhs,
                        "rhs": rhs,
                        "type": type_,
                    }
                )

            self._add_event(label, trigger_formula, assignments, trigger_distance, event_type)

    def _format_function_definitions(self) -> None:
        """Add function definitions to template variables."""
        for fd in self._sbml_function_definitions:
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
        for ia in self._sbml_initial_assignments:
            if ia.getMath() is not None:
                id_ = ia.getId()
                label = ia.getName().strip()
                var = ia.getSymbol()
                math = ia.getMath()
                self._add_initial_assignment(id_, label, var, math)

    def _format_parameters(self) -> None:
        """Add parameters to template variables."""
        # Note: rules must be processed before parameters
        assignment_rules = {rule["var"]: rule["math"] for rule in self._assignment_rules}
        initial_assignments = {ia["var"]: ia["math"] for ia in self._initial_assignments}

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
                d_param_id = state_var["derivative_id"]
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
                        {
                            "id": param_id,
                            "label": param.getName().strip(),
                            "value": str(param.getValue()),
                        }
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
        """Add species to template variables."""
        # Note: Rules must be processed before species
        assignment_rules = {ar["var"]: ar["math"] for ar in self._assignment_rules}
        rate_rules = {rr["var"]: rr["math"] for rr in self._rate_rules}
        initial_assignments = {ia["var"]: ia["math"] for ia in self._initial_assignments}

        model_conversion_factor = None
        if self._sbml_model.isSetConversionFactor():
            model_conversion_factor = self._sbml_model.getConversionFactor()

        for species in self._sbml_species:
            species_id = species.getId()
            label = species.getName().strip()

            conversion_factor = None
            if species.isSetConversionFactor():
                conversion_factor = species.getConversionFactor()
            elif model_conversion_factor is not None:
                conversion_factor = model_conversion_factor

            # If there's a compartment we'll normalise the ODEs, so declare it as non-dimensional
            compartment_id = species.getCompartment()
            compartment = self._sbml_compartments.get(compartment_id)

            units = NON_DIM_UNITS if compartment else species.getSubstanceUnits()
            has_only_substance_units = species.isSetHasOnlySubstanceUnits() and species.getHasOnlySubstanceUnits()

            initial_value = None
            if species.isSetInitialConcentration():
                initial_value = species.getInitialConcentration()
                math = parseL3Formula(f"{initial_value}")
                self._add_equation(
                    var=species_id,
                    math=math,
                    eq_type=EquationType.INITIAL_VALUE,
                )

                if (
                    has_only_substance_units
                    and (species_id not in assignment_rules)
                    and (species_id not in initial_assignments)
                ):
                    # Convert initial concentration to amount via a custom initial assignment
                    ia_id = PREFIX_SEP.join([CHASTE_PREFIX, INITIAL_ASSIGNMENT_PREFIX, species_id])
                    ia_label = f"Convert {species_id} concentration to amount"
                    ia_var = species_id
                    ia_rhs = f"{species_id} * {compartment_id}"
                    ia_math = parseL3Formula(ia_rhs)
                    self._add_initial_assignment(ia_id, ia_label, ia_var, ia_math)
                    self._add_equation(var=ia_var, math=ia_math, eq_type=EquationType.INITIAL_ASSIGNMENT)

            elif species.isSetInitialAmount():
                initial_value = species.getInitialAmount()
                math = parseL3Formula(f"{initial_value}")
                self._add_equation(species_id, math, eq_type=EquationType.INITIAL_VALUE)

                if (
                    not has_only_substance_units
                    and (species_id not in assignment_rules)
                    and (species_id not in initial_assignments)
                ):
                    # Convert initial amount to concentration via a custom initial assignment
                    ia_id = PREFIX_SEP.join([CHASTE_PREFIX, INITIAL_ASSIGNMENT_PREFIX, species_id])
                    ia_label = f"Convert {species_id} amount to concentration"
                    ia_var = species_id
                    ia_rhs = f"{species_id} / {compartment_id}"
                    ia_math = parseL3Formula(ia_rhs)
                    self._add_initial_assignment(ia_id, ia_label, ia_var, ia_math)
                    self._add_equation(var=ia_var, math=ia_math, eq_type=EquationType.INITIAL_ASSIGNMENT)

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
                if not has_only_substance_units:
                    if compartment_id in self._odes:
                        compartment_rhs = formulaToL3String(self._odes[compartment_id])
                    elif compartment_id in assignment_rules:
                        compartment_rhs = self._total_time_derivative(assignment_rules[compartment_id])
                    else:
                        compartment_rhs = ""

                if not has_only_substance_units and compartment_rhs:
                    # State variable: boundary condition in changing compartment
                    rhs = f"(-{species_id} * ({compartment_rhs})) / {compartment_id}"
                    if conversion_factor is not None:
                        rhs = f"({rhs}) * {conversion_factor}"
                    math = parseL3Formula(rhs)
                    state_var = self._add_state_variable(species_id, label, initial_value, units)
                    self._add_equation(var=state_var["derivative_id"], math=math, eq_type=EquationType.DERIVATIVE)
                else:
                    # Derived quantity
                    self._add_derived_quantity(species_id, label, initial_value, units)

            elif species_id in self._odes:
                # State variable
                math = self._odes[species_id]
                rhs = formulaToL3String(math)

                # Add parentheses if there are multiple terms
                if "+" in rhs or "-" in rhs[1:]:
                    rhs = f"({rhs})"

                # Add compartment scaling if defined by a reaction
                if not (has_only_substance_units or species_id in rate_rules):
                    # Determine dC/dt for the dilution correction
                    if compartment_id in self._odes:
                        compartment_ddt = formulaToL3String(self._odes[compartment_id])
                    elif compartment_id in assignment_rules:
                        compartment_ddt = self._total_time_derivative(assignment_rules[compartment_id])
                    else:
                        compartment_ddt = ""

                    if compartment_ddt:
                        rhs = f"({rhs} - {species_id} * ({compartment_ddt}))"

                    # Scale by compartment volume
                    rhs = f"{rhs} / {compartment_id}"

                # TODO: Handle time scaling
                # time_multiplier = self._get_timescale_multiplier()
                # if time_multiplier != 1.0:
                #     f"rDY[{state_variable_index}] *= {time_multiplier};"

                if conversion_factor is not None:
                    rhs = f"({rhs}) * {conversion_factor}"

                state_var = self._add_state_variable(species_id, label, initial_value, units)
                self._add_equation(
                    var=state_var["derivative_id"],
                    math=parseL3Formula(rhs),
                    eq_type=EquationType.DERIVATIVE,
                )

            else:
                # Derived quantity (constant value)
                rhs = str(initial_value)
                self._add_derived_quantity(species_id, label, initial_value, units)

            if not has_only_substance_units:
                # Add an extra "amount" derived quantity for the species
                self._add_amount(species)

    def _convert_infix_operator_to_function_syntax(self, formula: str, operator: str, function_name: str) -> str:
        """Convert infix operator expressions to function syntax.

        Example: with operator='^' and function_name='pow', rewrites
        ``a ^ b`` to ``pow(a, b)``.

        This parser handles parenthesized operands (including nested parentheses)
        and simple symbolic/numeric tokens.

        :param formula: Formula text in SBML infix style.
        :param operator: Infix operator token to rewrite.
        :param function_name: Function name used for the replacement.
        :return: Formula text with infix operators rewritten as function calls.
        """
        if not operator:
            raise ValueError("operator must be a non-empty string")
        if not function_name:
            raise ValueError("function_name must be a non-empty string")

        token_chars = set("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_:.eE")

        def read_left_operand(s: str, operator_index: int) -> tuple[int, int] | None:
            i = operator_index - 1
            while i >= 0 and s[i].isspace():
                i -= 1
            if i < 0:
                return None

            if s[i] == ")":
                depth = 1
                j = i - 1
                while j >= 0:
                    if s[j] == ")":
                        depth += 1
                    elif s[j] == "(":
                        depth -= 1
                        if depth == 0:
                            return (j, i + 1)
                    j -= 1
                return None

            j = i
            while j >= 0 and s[j] in token_chars:
                j -= 1

            start = j + 1
            end = i + 1
            if start >= end:
                return None
            return (start, end)

        def read_right_operand(s: str, operator_index: int) -> tuple[int, int] | None:
            i = operator_index + len(operator)
            n = len(s)
            while i < n and s[i].isspace():
                i += 1
            if i >= n:
                return None

            if s[i] == "(":
                depth = 1
                j = i + 1
                while j < n:
                    if s[j] == "(":
                        depth += 1
                    elif s[j] == ")":
                        depth -= 1
                        if depth == 0:
                            return (i, j + 1)
                    j += 1
                return None

            j = i
            if s[j] in "+-":
                j += 1

            while j < n and s[j] in token_chars:
                j += 1

            if j <= i:
                return None
            return (i, j)

        max_rewrites = 500
        rewrites = 0
        search_start = 0
        while rewrites < max_rewrites:
            operator_index = formula.find(operator, search_start)
            if operator_index < 0:
                break

            left = read_left_operand(formula, operator_index)
            right = read_right_operand(formula, operator_index)
            if left is None or right is None:
                search_start = operator_index + len(operator)
                continue

            l_start, l_end = left
            r_start, r_end = right

            lhs = formula[l_start:l_end].strip()
            rhs = formula[r_start:r_end].strip()
            replacement = f"{function_name}({lhs}, {rhs})"

            formula = formula[:l_start] + replacement + formula[r_end:]
            search_start = l_start + len(replacement)
            rewrites += 1

        return formula

    def _convert_infix_power_to_function_syntax(self, formula: str) -> str:
        """Convert infix power expressions (a ^ b) to function syntax pow(a, b)."""
        return self._convert_infix_operator_to_function_syntax(
            formula=formula,
            operator="^",
            function_name="pow",
        )

    def _formula_to_string(self, math: "ASTNode") -> str:
        """Convert an AST math formula to an equivalent C++ string.

        :param math: The AST math formula.
        :return: The equivalent C++ string.
        """
        unsupported_functions = ["delay"]
        for func in unsupported_functions:
            if search_ast_type(math, AST_FUNCTION_DELAY):
                raise NotImplementedError(f"SBML function not supported: '{func}'.")

        formula = formulaToL3String(math)

        # Convert all integer literals to doubles to fix integer division.
        # TODO: Perhaps instead of regex, traverse AST and convert some AST_INTEGER
        # nodes to AST_REAL. This should only need to apply to numbers used in a division.
        formula = re.sub(r"(?<!\.)(?<!e-|E-)\b[0-9]+\b(?!\.)", lambda x: f"{x[0]}.0", formula)

        formula = self._convert_infix_power_to_function_syntax(formula)

        # SBML contants to be replaced with C++ equivalents
        constants = {
            "avogadro": "sm::AVOGADRO",
            "exponentiale": "M_E",
            "inf": "std::numeric_limits<double>::infinity()",
            "infinity": "std::numeric_limits<double>::infinity()",
            "nan": "NAN",
            "notanumber": "NAN",
            "pi": "M_PI",
            "time": "time",
            "t": "time",
            "s": "time",
        }
        # skip: "true", "false"

        # SBML functions with same name as C++ equivalents
        unchanged_functions = {
            "acos",
            "acosh",
            "asin",
            "asinh",
            "atan",
            "atanh",
            "ceil",
            "cos",
            "cosh",
            "exp",
            "floor",
            "pow",
            "sin",
            "sinh",
            "sqrt",
            "tan",
            "tanh",
        }

        # SBML functions with different names in C++
        renamed_functions = {
            "abs": "fabs",
            "arccos": "acos",
            "arccosh": "acosh",
            "arcsin": "asin",
            "arcsinh": "asinh",
            "arctan": "atan",
            "arctanh": "atanh",
            "ceiling": "ceil",
            "ln": "log",
            "power": "pow",
            "rem": "fmod",
        }

        # SBML functions with custom implementations
        custom_functions = {
            "and": "and_",
            "acot": "acot",
            "acoth": "acoth",
            "acsc": "acsc",
            "acsch": "acsch",
            "asec": "asec",
            "asech": "asech",
            "arccot": "acot",
            "arccoth": "acoth",
            "arccsc": "acsc",
            "arccsch": "acsch",
            "arcsec": "asec",
            "arcsech": "asech",
            "cot": "cot",
            "coth": "coth",
            "csc": "csc",
            "csch": "csch",
            "eq": "eq",
            "factorial": "factorial",
            "geq": "geq",
            "gt": "gt",
            "leq": "leq",
            "log": "log",
            "lt": "lt",
            "max": "max",
            "min": "min",
            "neq": "neq",
            "not": "not_",
            "or": "or_",
            "piecewise": "piecewise",
            "quotient": "quotient",
            "root": "root",
            "sec": "sec",
            "sech": "sech",
            "sqr": "sqr",
            "xor": "xor_",
        }

        # TODO: From SBML Level 3 upwards, log defaults to base 10.
        # SBML versions lower than 3 default to base e.
        # See https://sbml.org/software/libsbml/5.18.0/docs/formatted/python-api/namespacelibsbml.html#a8e96a5a70569ae32655c6302638f6dc3  # noqa: B950

        tokens = re.findall(r"\w+|\W+", formula)

        cpp_tokens = []
        for token in tokens:
            cpp_token = token

            # Replace function names and constants.
            if token in constants:
                cpp_token = f"{constants[token]}"

            elif token in unchanged_functions:
                cpp_token = f"std::{token}"

            elif token in renamed_functions:
                cpp_token = f"std::{renamed_functions[token]}"

            elif token in custom_functions:
                cpp_token = f"sm::{custom_functions[token]}"

            cpp_tokens.append(cpp_token)
        cpp_formula = "".join(cpp_tokens)

        results = re.findall(r"rateOf\(([^)]+)\)", cpp_formula)
        if results:
            for var in results:
                rate = "0.0"
                var_type = self._get_variable_type(var)
                if var_type == VarType.STATE_VARIABLE:
                    i = self._get_variable_index(var)
                    state_var = self._state_variables[i]
                    rate = state_var["derivative_id"]
                cpp_formula = cpp_formula.replace(f"rateOf({var})", rate)
        return cpp_formula

    def _generate_output(self, template_path, filename) -> None:
        """Generate a single output file from a template.

        Generated code is stored in the outputs list.

        :param template_path: The path to the template.
        :param filename: The output filename.
        """
        template = self._get_template(template_path)
        code = template.render(self._template_vars)
        self._add_output(filename, code)

    def _generate_outputs(self) -> None:
        """Generate Chaste code for the model."""
        # Generate code for the OdeSystem
        self._generate_output("ode/ode.hpp", self._ode_hpp_filename)
        self._generate_output("ode/ode.cpp", self._ode_cpp_filename)

        # Generate code for the SRN or Cell-Cycle model
        if self._model_type == ModelType.SRN:
            self._generate_output("srn/srn.hpp", self._srn_hpp_filename)
            self._generate_output("srn/srn.cpp", self._srn_cpp_filename)

        elif self._model_type == ModelType.CELL_CYCLE:
            self._generate_output("cell_cycle/cell_cycle.hpp", self._cell_cycle_hpp_filename)
            self._generate_output("cell_cycle/cell_cycle.cpp", self._cell_cycle_cpp_filename)

    def _get_template(self, name: str) -> "Template":
        """Get a Jinja2 template.

        :param name: The template name.
        :return: The template object.
        """
        return self._jinja_env.get_template(name)

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
        """Get the index of a variable.

        :param id_: The variable ID.
        :return: The variable index.
        """
        # TODO: Make this more generic
        var_type = self._get_variable_type(id_)

        if var_type == VarType.AMOUNT:
            for amount in self._amounts:
                if amount["id"] == id_:
                    return amount["index"]

        elif var_type == VarType.DERIVED_QUANTITY:
            for dq in self._derived_quantities:
                if dq["id"] == id_:
                    return dq["index"]

        elif var_type == VarType.FUNCTION:
            for func in self._functions:
                if func["id"] == id_:
                    return func["index"]

        elif var_type == VarType.PARAMETER:
            for param in self._parameters:
                if param["id"] == id_:
                    return param["index"]

        elif var_type == VarType.REACTION:
            for reaction in self._reactions:
                if reaction["id"] == id_:
                    return reaction["index"]

        elif var_type == VarType.STATE_VARIABLE:
            for state_variable in self._state_variables:
                if state_variable["id"] == id_:
                    return state_variable["index"]

        raise ValueError(f"ID '{id_}' is not a recognized variable.")

    def _get_variable_type(self, var_id: str) -> bool:
        """Get the type of a variable based on its ID.

        :param var_id: The variable ID.
        """
        return self._variable_types.get(var_id, VarType.UNKNOWN)

    def _populate_template_vars(self) -> None:
        """Populate the template variables for generating C++ code."""
        template_vars: dict[str, "Any"] = dict(
            amounts=self._amounts,
            assignment_rules=self._assignment_rules,
            derived_quantities=self._derived_quantities,
            equations=self._equations,
            events=self._events,
            functions=self._functions,
            ode_class_name=self._ode_class_name,
            ode_header_guard=generate_header_guard(self._ode_hpp_filename),
            ode_hpp_file=self._ode_hpp_filename,
            parameters=self._parameters,
            reactions=self._reactions,
            state_variables=self._state_variables,
            stoichiometry_variables=self._stoichiometry_variables,
        )

        if self._model_type == ModelType.SRN:
            template_vars.update(
                dict(
                    srn_class_name=self._srn_class_name,
                    srn_header_guard=generate_header_guard(self._srn_hpp_filename),
                    srn_hpp_file=self._srn_hpp_filename,
                )
            )
        elif self._model_type == ModelType.CELL_CYCLE:
            template_vars.update(
                dict(
                    cell_cycle_class_name=self._cell_cycle_class_name,
                    cell_cycle_header_guard=generate_header_guard(self._cell_cycle_hpp_filename),
                    cell_cycle_hpp_file=self._cell_cycle_hpp_filename,
                )
            )
        self._template_vars = template_vars

    def _process_model(self) -> None:
        """Process the SBML model to set up the formatted variables for templates."""
        self._variable_types = {}
        self._odes = {}

        self._assignment_rules = []
        self._rate_rules = []
        self._initial_assignments = []

        self._state_variables = []
        self._derived_quantities = []
        self._amounts = []
        self._parameters = []
        self._stoichiometry_variables = []
        self._equations = []

        self._reactions = []
        self._events = []
        self._functions = []

        # TODO: enforce processing order e.g. rules must be processed first
        self._format_rules()
        self._extract_odes()

        self._format_initial_assignments()

        self._format_compartments()
        self._format_species()
        self._format_parameters()

        if not self._state_variables:
            raise NotImplementedError("Models without state variables are not supported.")

        self._format_reactions()
        self._format_function_definitions()
        self._format_events()

        self._format_equations()

        self._populate_template_vars()

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

        def _depends(eq: dict, var: str) -> bool:
            """Check if a variable appears in the rhs of an equation.

            :param eq: The equation to check.
            :param var: The variable to check for.
            :return: True if the variable appears in the rhs of the equation, False otherwise.
            """
            # False if there is a local parameter with the same name as the variable.
            if eq["local_parameters"]:
                if any(var == param["id"] for param in eq["local_parameters"]):
                    return False

            return bool(re.search(rf"\b{var}\b", eq["rhs"]))

        # Sort equations by group
        group_order = [
            EquationType.INITIAL_VALUE,
            EquationType.INITIAL_ASSIGNMENT,
            EquationType.ASSIGNMENT_RULE,
            EquationType.REACTION,
            EquationType.DERIVATIVE,
            EquationType.AMOUNT,
            EquationType.UNKNOWN,
        ]

        sorted_equations = [eq for eq_type in group_order for eq in self._equations if eq["type"] == eq_type]

        # Sort equations by dependencies
        for _ in range(len(sorted_equations)):  # Max iterations for worst case
            re_sorted_equations = []
            for eq in sorted_equations:
                for i, other_eq in enumerate(re_sorted_equations):
                    # Insert eq before the first other equation that depends on it.
                    if _depends(other_eq, eq["var"]):
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
