"""Module for the ChasteSbmlModel class."""

import abc
import os
import pathlib
import subprocess
from typing import TYPE_CHECKING

from jinja2 import Environment, PackageLoader, select_autoescape
from libsbml import (  # AST_INTEGER,; AST_RATIONAL,; AST_REAL,; AST_REAL_E,
    AST_RELATIONAL_EQ,
    AST_RELATIONAL_GEQ,
    AST_RELATIONAL_GT,
    AST_RELATIONAL_LEQ,
    AST_RELATIONAL_LT,
    AST_RELATIONAL_NEQ,
    SBML_ALGEBRAIC_RULE,
    SBML_ASSIGNMENT_RULE,
    SBML_RATE_RULE,
    SBMLReader,
)

from ._config import (
    AMOUNT_PREFIX,
    CHASTE_PREFIX,
    INITIAL_ASSIGNMENT_PREFIX,
    NON_DIM_UNITS,
    PREFIX_SEP,
    ROOT_DIR,
    EventType,
    ModelType,
    VarType,
)
from ._utils import (
    convert_ast_formula,
    generate_header_guard,
    get_compartment_size,
    get_function_definition_arguments,
    sort_formulas,
    to_camel_case,
    to_cpp_name,
)

if TYPE_CHECKING:
    from typing import Any

    from jinja2.environment import Template
    from libsbml import Reaction, SpeciesReference


class ChasteSbmlModel:
    """Holds information about an SBML model for which Chaste code is to be generated."""

    __metaclass__ = abc.ABCMeta

    _jinja_env = Environment(
        loader=PackageLoader("chaste_codegen_sbml"),
        autoescape=select_autoescape(),
        trim_blocks=True,
        lstrip_blocks=True,
    )
    _jinja_env.globals["VarType"] = VarType
    _jinja_env.globals["EventType"] = EventType
    _jinja_env.globals["AMOUNT_PREFIX"] = AMOUNT_PREFIX
    _jinja_env.globals["PREFIX_SEP"] = PREFIX_SEP

    # -- PUBLIC --------------------------------------

    def __init__(
        self, sbml_file: str, model_name: str = "", model_type: ModelType = ModelType.GENERIC
    ) -> None:
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

        self._srn_class_name = None
        self._srn_hpp_filename = None
        self._srn_cpp_filename = None
        if self._model_type == ModelType.SRN:
            self._srn_class_name = self._model_name + "SrnModel"
            self._srn_hpp_filename = f"{self._srn_class_name}.hpp"
            self._srn_cpp_filename = f"{self._srn_class_name}.cpp"

        self._cell_cycle_class_name = None
        self._cell_cycle_hpp_filename = None
        self._cell_cycle_cpp_filename = None
        if self._model_type == ModelType.CELL_CYCLE:
            self._cell_cycle_class_name = self._model_name + "CellCycleModel"
            self._cell_cycle_hpp_filename = f"{self._cell_cycle_class_name}.hpp"
            self._cell_cycle_cpp_filename = f"{self._cell_cycle_class_name}.cpp"

        self._sbml_model = SBMLReader().readSBMLFromFile(self._sbml_file).getModel()
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
        self._variable_parameters = []  # [ { id: str, label: str, ... } ]
        self._constant_parameters = []  # [ { id: str, label: str, ... } ]
        self._initial_assignments = []  # [ { id: str, label: str, ... } ]
        self._reactions = []  # [ { id: str, label: str, ... } ]
        self._events = []  # [ { name: str, trigger: str, ... } ]
        self._functions = []  # [ { name: str, args: [str], body: str } ]
        self._stoichiometry_variables = []  # [ { id: str, label: str, ... } ]

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

        # Format with clang-format
        for filename in self._outputs:
            file_path = str(root_dir / filename)
            subprocess.run(
                [
                    "clang-format",
                    "-i",
                    f"-style=file:{ROOT_DIR}/.clang-format",
                    str(file_path),
                ],
                check=True,
            )

    # -- PRIVATE ---------------------------------------

    def _add_amount(self, id_: str, units: str, rhs: str) -> None:
        """Add an amount derived quantity to the template variables.

        :param id_: The variable ID.
        :param label: The variable description.
        :param initial_value: The variable initial value.
        :param units: The variable units.
        :param rhs: The variable formula.
        """
        amount_id = AMOUNT_PREFIX + PREFIX_SEP + id_
        self._amounts.append(
            {
                "id": amount_id,
                "index": len(self._amounts),
                "rhs": rhs,
                "units": units,
            }
        )
        self._variable_types[amount_id] = VarType.AMOUNT

    def _add_assignment_rule(self, id_: str, label: str, lhs: str, rhs: str) -> None:
        """Add an assignment rule to the template variables.

        :param id_: The rule ID.
        :param label: The rule description.
        :param lhs: The left-hand side of the rule.
        :param rhs: The right-hand side of the rule.
        """
        self._assignment_rules.append(
            {
                "id": id_,
                "label": label,
                "index": len(self._assignment_rules),
                "lhs": lhs,
                "rhs": rhs,
            }
        )
        self._variable_types[id_] = VarType.ASSIGNMENT_RULE

    def _add_constant_parameter(
        self, id_: str, label: str, value: float, units: str, rhs: str
    ) -> None:
        """Add a constant parameter to the template variables.

        :param id_: The parameter ID.
        :param label: The parameter description.
        :param value: The parameter value.
        :param units: The parameter units.
        """
        self._constant_parameters.append(
            {
                "id": id_,
                "label": label,
                "index": len(self._constant_parameters),
                "rhs": rhs,
                "value": value,
                "units": units,
            }
        )
        self._variable_types[id_] = VarType.CONSTANT_PARAMETER

    def _add_derived_quantity(
        self, id_: str, label: str, initial_value: float, units: str, rhs: str
    ) -> None:
        """Add a derived quantity to the template variables.

        :param id_: The variable ID.
        :param label: The variable description.
        :param initial_value: The variable initial value.
        :param units: The variable units.
        :param rhs: The variable formula.
        """
        self._derived_quantities.append(
            {
                "id": id_,
                "label": label,
                "index": len(self._derived_quantities),
                "initial_value": initial_value,
                "rhs": rhs,
                "units": units,
            }
        )
        self._variable_types[id_] = VarType.DERIVED_QUANTITY

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

    def _add_initial_assignment(
        self, id_: str, label: str, lhs: str, rhs: str, custom: bool = False
    ) -> None:
        """Add an initial assignment to the template variables.

        :param id_: The rule ID.
        :param label: The rule description.
        :param lhs: The left-hand side of the rule.
        :param rhs: The right-hand side of the rule.
        :param custom: True if this is a custom initial assignment added during processing.
        """
        self._initial_assignments.append(
            {
                "id": id_,
                "label": label,
                "index": len(self._initial_assignments),
                "lhs": lhs,
                "rhs": rhs,
                "custom": custom,
            }
        )
        self._variable_types[id_] = VarType.INITIAL_ASSIGNMENT

    def _add_output(self, filename: str, code: str) -> None:
        """Add generated code to the outputs dictionary.

        :param filename: The filename.
        :param code: The code.
        """
        self._outputs[filename] = code

    def _add_rate_rule(self, id_: str, label: str, lhs: str, rhs: str) -> None:
        """Add a rate rule to the template variables.

        :param id_: The rule ID.
        :param label: The rule description.
        :param lhs: The left-hand side of the rule.
        :param rhs: The right-hand side of the rule.
        """
        self._rate_rules.append(
            {
                "id": id_,
                "label": label,
                "index": len(self._rate_rules),
                "lhs": lhs,
                "rhs": rhs,
            }
        )
        self._variable_types[id_] = VarType.RATE_RULE

    def _add_reaction(
        self, id_: str, label: str, rhs: str, parameters: list[dict[str, "Any"]]
    ) -> None:
        """Add a reaction to the template variables.

        :param id_: The reaction ID.
        :param label: The reaction description.
        :param rhs: The reaction formula.
        :param parameters: The reaction parameters.
        """
        self._reactions.append(
            {
                "id": id_,
                "label": label,
                "index": len(self._reactions),
                "rhs": rhs,
                "parameters": parameters,
            }
        )
        self._variable_types[id_] = VarType.REACTION

    def _add_state_variable(
        self, id_: str, label: str, initial_value: float, units: str, rhs: str
    ) -> None:
        """Add a state variable to the template variables.

        :param id_: The variable ID.
        :param label: The variable description.
        :param initial_value: The variable initial value.
        :param units: The variable units.
        :param rhs: The variable formula.
        """
        self._state_variables.append(
            {
                "id": id_,
                "label": label,
                "index": len(self._state_variables),
                "initial_value": initial_value,
                "rhs": rhs,
                "units": units,
            }
        )
        self._variable_types[id_] = VarType.STATE_VARIABLE

    def _add_stoichiometry_variable(
        self, id_: str, label: str, initial_value: float, rhs: str
    ) -> None:
        """Add a stoichiometry variable to the template variables.

        :param id_: The ID of the stoichiometry variable.
        :param label: The description of the stoichiometry variable.
        """
        self._stoichiometry_variables.append(
            {
                "id": id_,
                "index": len(self._stoichiometry_variables),
                "label": label,
                "initial_value": initial_value,
                "rhs": rhs,
            }
        )
        self._variable_types[id_] = VarType.STOICHIOMETRY_VARIABLE

    def _add_variable_parameter(
        self, id_: str, label: str, initial_value: float, units: str = NON_DIM_UNITS
    ) -> None:
        """Add a variable parameter to the template variables.

        :param id_: The parameter ID.
        :param label: The parameter description.
        :param initial_value: The parameter initial value.
        :param units: The parameter units.
        """
        self._variable_parameters.append(
            {
                "id": id_,
                "label": label,
                "index": len(self._variable_parameters),
                "initial_value": initial_value,
                "units": units,
            }
        )
        self._variable_types[id_] = VarType.VARIABLE_PARAMETER

    def _extract_odes(self) -> None:
        """Extract the ODEs equations for each species."""
        # Note: rules must be processed before ODE extraction
        self._odes = {}

        def _update_ode(species_ref: "SpeciesReference", rxn: "Reaction", is_product: bool) -> None:
            """Update the ODE for a species based on a species reference in a reaction.

            :param species_ref: The species reference.
            :param rxn: The reaction.
            :param is_product: True if the species reference is a product, False if a reactant.
            """
            rhs = rxn.getId()

            # Account for stoichiometry
            sto_value = None
            if species_ref.isSetStoichiometry():
                sto_value = species_ref.getStoichiometry()

            sto_formula = None
            if species_ref.isSetStoichiometryMath():
                sto_math = species_ref.getStoichiometryMath().getMath()
                sto_formula = convert_ast_formula(sto_math)

            if species_ref.isSetId():
                sto_id = species_ref.getId()
                sto_label = species_ref.getName().strip()
                self._add_stoichiometry_variable(sto_id, sto_label, sto_value, sto_formula)

                rhs = f"({sto_id} * {rhs})"

            elif sto_formula is not None:
                rhs = f"({sto_formula} * {rhs})"

            elif sto_value is not None and sto_value != 1:
                rhs = f"({sto_value} * {rhs})"

            # Update the ODE
            species_id = species_ref.getSpecies()
            if species_id in self._odes:
                if is_product:
                    self._odes[species_id] += f" + {rhs}"
                else:
                    self._odes[species_id] += f" - {rhs}"
            else:
                if is_product:
                    self._odes[species_id] = rhs
                else:
                    self._odes[species_id] = f"-{rhs}"

        # Extract ODEs from reactions:
        # each ODE will essentially be the sum of the products minus the
        # sum of the reactants divided by the compartment volume
        for reaction in self._sbml_reactions:

            # Decompose reaction into sum of products minus sum of reactants
            for product in reaction.getListOfProducts():
                _update_ode(product, reaction, is_product=True)

            for reactant in reaction.getListOfReactants():
                _update_ode(reactant, reaction, is_product=False)

        # Extract ODEs from rate rules
        for rule in self._rate_rules:
            lhs = rule["lhs"]
            rhs = rule["rhs"]
            if lhs in self._odes:
                raise ValueError(f"{lhs} has more than one rate rule and/or reaction.")
            self._odes[lhs] = rhs

        # Extract ODEs from assignment rules
        # assignment_rules = [r for r in self._sbml_rules if r.getTypeCode() == SBML_ASSIGNMENT_RULE]
        # compartments = {c.getId() for c in self._sbml_compartments}

        # for rule in assignment_rules:
        #     lhs = rule.getVariable()
        #     if lhs in self._odes:
        #         raise ValueError(f"{lhs} has both an assignment rule and a rate rule/reaction.")

        #     # Check for compartments that change size
        #     if lhs in compartments:
        #         math = rule.getMath()
        #         if math.getType() not in [AST_INTEGER, AST_REAL, AST_REAL_E, AST_RATIONAL]:
        #             rhs = convert_ast_formula(math)
        #             rhs = f"({rhs}) - {lhs}"
        #             self._odes[lhs] = rhs

        if not self._odes:
            raise NotImplementedError("Models without ODEs are not supported.")

    def _format_compartments(self) -> None:
        """Add compartments to template variables."""
        # Note: rules must be processed before compartments
        assignment_rules = {rule["lhs"]: rule["rhs"] for rule in self._assignment_rules}

        for compartment in self._sbml_compartments:
            compartment_id = compartment.getId()
            label = compartment.getName().strip()

            value = get_compartment_size(compartment)
            units = compartment.getUnits() if compartment.isSetUnits() else NON_DIM_UNITS

            if compartment.isSetConstant() and compartment.getConstant():
                # Derived quantity
                self._add_derived_quantity(compartment_id, label, value, units, rhs=None)
            elif compartment_id in self._odes:
                # State variable
                rhs = self._odes[compartment_id]
                self._add_state_variable(compartment_id, label, value, units, rhs)
            elif compartment_id in assignment_rules:
                # Derived quantity
                rhs = assignment_rules[compartment_id]
                self._add_derived_quantity(compartment_id, label, value, units, rhs)
            else:
                # Variable parameter
                self._add_variable_parameter(compartment_id, label, value, units)

    def _format_events(self) -> None:
        """Add events to template variables."""
        # TODO: Add priority

        for event in self._sbml_events:
            if event.isSetDelay():
                math = convert_ast_formula(event.getDelay().getMath())
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
            trigger_formula = convert_ast_formula(trigger_math)

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
                lc = convert_ast_formula(trigger_math.getLeftChild())
                rc = convert_ast_formula(trigger_math.getRightChild())

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
                    trigger_distance = (
                        f"std::abs(({lc}) - ({rc})) - std::numeric_limits<double>::epsilon()"
                    )

                # TODO: Distance calculation assumes two operands. Extend to more operands?
                # e.g. trigger: geq(3.0, 6.0, 7.0, 9.0) -> condition=false, dist=min(3.0, 1.0, 2.0)=1.0

            assignments = []
            for assignment in event.getListOfEventAssignments():
                lhs = assignment.getVariable()
                type_ = self._get_variable_type(lhs)
                index = self._get_variable_index(lhs)
                rhs = convert_ast_formula(assignment.getMath())

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
            body = convert_ast_formula(fd.getBody())

            self._add_function(
                fd_id,
                label,
                args,
                body,
            )

    def _format_initial_assignments(self) -> None:
        """Add initial assignments to template variables."""
        # Sort initial assignments - variables on rhs must be defined before they are used
        initial_assignments = [
            ia for ia in self._sbml_initial_assignments if ia.getMath() is not None
        ]

        # TODO: Don't repeat this code from _format_rules
        assignment_rules = [
            r
            for r in self._sbml_rules
            if r.getTypeCode() == SBML_ASSIGNMENT_RULE and r.getMath() is not None
        ]

        assignments = initial_assignments + assignment_rules

        formulas = [(ia.getSymbol(), ia.getMath()) for ia in initial_assignments]
        formulas += [(r.getVariable(), r.getMath()) for r in assignment_rules]

        sort_index = sort_formulas(formulas)
        sorted_assignments = [assignments[i] for i in sort_index]
        sorted_formulas = [formulas[i] for i in sort_index]

        for i, assignment in enumerate(sorted_assignments):
            id_ = assignment.getId()
            label = assignment.getName().strip()
            lhs = sorted_formulas[i][0]
            rhs = convert_ast_formula(sorted_formulas[i][1])
            self._add_initial_assignment(id_, label, lhs, rhs, custom=False)

    def _format_parameters(self) -> None:
        """Add parameters to template variables."""
        # Note: rules must be processed before parameters
        assignment_rules = {rule["lhs"]: rule["rhs"] for rule in self._assignment_rules}
        initial_assignments = {ia["lhs"]: ia["rhs"] for ia in self._initial_assignments}

        for param in self._sbml_parameters:
            param_id = param.getId()
            label = param.getName().strip()

            value = None
            if param.isSetValue():
                value = param.getValue()

            units = NON_DIM_UNITS
            if param.isSetUnits():
                units = param.getUnits()

            if param_id in self._odes:
                # State variable
                rhs = self._odes[param_id]
                self._add_state_variable(param_id, label, value, units, rhs)
            elif param_id in assignment_rules:
                # Derived quantity
                rhs = assignment_rules[param_id]
                self._add_derived_quantity(param_id, label, value, units, rhs)
            elif param.isSetConstant() and param.getConstant():
                # Constant parameter
                if param_id in initial_assignments:
                    rhs = initial_assignments[param_id]
                else:
                    rhs = None
                self._add_constant_parameter(param_id, label, value, units, rhs)
            else:
                # Variable parameter
                self._add_variable_parameter(param_id, label, value, units)

    def _format_reactions(self) -> None:
        """Add reactions to template variables."""
        for reaction in self._sbml_reactions:
            reaction_id = reaction.getId()
            label = reaction.getName().strip()

            rhs = ""
            reaction_parameters = []

            kinetic_law = reaction.getKineticLaw()
            if kinetic_law is None:
                species = [r.getSpecies() for r in reaction.getListOfReactants()]
                rhs = " * ".join(species)

            else:
                rhs = convert_ast_formula(kinetic_law.getMath())

                sbml_parameters = kinetic_law.getListOfParameters()
                for param in sbml_parameters:
                    param_id = param.getId()
                    param_descr = param.getName().strip()
                    param_value = param.getValue()  # TODO: What if not set?
                    reaction_parameters.append(
                        {
                            "id": param_id,
                            "label": param_descr,
                            "value": param_value,
                        }
                    )

            self._add_reaction(reaction_id, label, rhs, reaction_parameters)

    def _format_rules(self) -> None:
        """Add rules to template variables."""
        # Algebraic rules are not implemented
        if any(r.getTypeCode() == SBML_ALGEBRAIC_RULE for r in self._sbml_rules):
            raise NotImplementedError("Algebraic rules are not supported.")

        # Sort assignment rules - variables on rhs must be defined before they are used
        assignment_rules = [r for r in self._sbml_rules if r.getTypeCode() == SBML_ASSIGNMENT_RULE]
        formulas = [(r.getVariable(), r.getMath()) for r in assignment_rules]
        sort_index = sort_formulas(formulas)
        sorted_assignments = [assignment_rules[i] for i in sort_index]

        for rule in sorted_assignments:
            rule_id = rule.getId()
            label = rule.getName().strip()
            lhs = rule.getVariable()

            # if lhs not in self._odes: # Exclude assignments converted to ODEs
            math = rule.getMath()
            if math is not None:
                rhs = convert_ast_formula(math)
                self._add_assignment_rule(rule_id, label, lhs, rhs)

        # Get rate rules
        rate_rules = [r for r in self._sbml_rules if r.getTypeCode() == SBML_RATE_RULE]
        for rule in rate_rules:
            rule_id = rule.getId()
            label = rule.getName().strip()
            lhs = rule.getVariable()

            math = rule.getMath()
            if math is not None:
                rhs = convert_ast_formula(math)
                self._add_rate_rule(rule_id, label, lhs, rhs)

    def _format_species(self) -> None:
        """Add species to template variables."""
        # Note: Rules must be processed before species
        assignment_rules = {r["lhs"]: r["rhs"] for r in self._assignment_rules}
        rate_rules = {r["lhs"]: r["rhs"] for r in self._rate_rules}
        initial_assignments = {r["lhs"]: r["rhs"] for r in self._initial_assignments}

        for species in self._sbml_species:
            species_id = species.getId()
            label = species.getName().strip()

            # If there's a compartment we'll normalise the ODEs, so declare it as non-dimensional
            compartment_id = species.getCompartment()
            compartment = self._sbml_compartments.get(compartment_id)
            units = NON_DIM_UNITS if compartment else species.getSubstanceUnits()

            initial_value = None
            if species.isSetInitialConcentration():
                initial_value = species.getInitialConcentration()

            elif species.isSetInitialAmount():
                initial_value = species.getInitialAmount()

                if (species_id not in assignment_rules) and (species_id not in initial_assignments):
                    # Convert initial amount to concentration via a custom initial assignment
                    ia_id = PREFIX_SEP.join([CHASTE_PREFIX, INITIAL_ASSIGNMENT_PREFIX, species_id])
                    ia_label = f"Convert {species_id} amount to concentration"
                    ia_lhs = species_id
                    ia_rhs = f"{species_id} / {compartment_id}"
                    self._add_initial_assignment(ia_id, ia_label, ia_lhs, ia_rhs, custom=True)

            is_bc = species.isSetBoundaryCondition() and species.getBoundaryCondition()

            if species_id in assignment_rules:
                # Derived quantity (includes boundary conditions with assignment rules)
                rhs = assignment_rules[species_id]
                self._add_derived_quantity(species_id, label, initial_value, units, rhs)

            elif is_bc and species_id not in rate_rules:
                # Derived quantity (boundary condition)
                rhs = str(initial_value)
                self._add_derived_quantity(species_id, label, initial_value, units, rhs)

            elif species_id in self._odes:
                # State variable
                rhs = self._odes[species_id]

                # Add compartment ODE
                if compartment_id in self._odes:
                    compartment_rhs = self._odes[compartment_id]
                    rhs = f"{rhs} - {species_id} * ({compartment_rhs})"

                # Add parentheses if there are multiple terms
                if "+" in rhs or "-" in rhs[1:]:
                    rhs = f"({rhs})"

                # Add compartment scaling if defined by a reaction
                if species_id not in rate_rules:
                    rhs = f"{rhs} / {compartment_id}"

                # TODO: Handle time scaling
                # time_multiplier = self._get_timescale_multiplier()
                # if time_multiplier != 1.0:
                #     # This does not include species defined in algebraic rules
                #     f"rDY[{state_variable_index}] *= {time_multiplier};"

                self._add_state_variable(species_id, label, initial_value, units, rhs)

            else:
                # Derived quantity (constant value)
                rhs = str(initial_value)
                self._add_derived_quantity(species_id, label, initial_value, units, rhs)

            # Add an extra "amount" derived quantity for the species
            rhs = f"{species_id} * {compartment_id}"
            self._add_amount(species_id, units, rhs)

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

        elif var_type == VarType.ASSIGNMENT_RULE:
            for rule in self._assignment_rules:
                if rule["id"] == id_:
                    return rule["index"]

        elif var_type == VarType.CONSTANT_PARAMETER:
            for param in self._constant_parameters:
                if param["id"] == id_:
                    return param["index"]

        elif var_type == VarType.DERIVED_QUANTITY:
            for dq in self._derived_quantities:
                if dq["id"] == id_:
                    return dq["index"]

        elif var_type == VarType.FUNCTION:
            for func in self._functions:
                if func["id"] == id_:
                    return func["index"]

        elif var_type == VarType.INITIAL_ASSIGNMENT:
            for initial_assignment in self._initial_assignments:
                if initial_assignment["id"] == id_:
                    return initial_assignment["index"]

        elif var_type == VarType.STOICHIOMETRY_VARIABLE:
            for stoichiometry_variable in self._stoichiometry_variables:
                if stoichiometry_variable["id"] == id_:
                    return stoichiometry_variable["index"]

        elif var_type == VarType.RATE_RULE:
            for rate_rule in self._rate_rules:
                if rate_rule["id"] == id_:
                    return rate_rule["index"]

        elif var_type == VarType.REACTION:
            for reaction in self._reactions:
                if reaction["id"] == id_:
                    return reaction["index"]
                    return rate_rule["index"]

        elif var_type == VarType.STATE_VARIABLE:
            for state_variable in self._state_variables:
                if state_variable["id"] == id_:
                    return state_variable["index"]

        elif var_type == VarType.VARIABLE_PARAMETER:
            for param in self._variable_parameters:
                if param["id"] == id_:
                    return param["index"]

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
            constant_parameters=self._constant_parameters,
            derived_quantities=self._derived_quantities,
            events=self._events,
            functions=self._functions,
            initial_assignments=self._initial_assignments,
            ode_class_name=self._ode_class_name,
            ode_header_guard=generate_header_guard(self._ode_hpp_filename),
            ode_hpp_file=self._ode_hpp_filename,
            reactions=self._reactions,
            state_variables=self._state_variables,
            stoichiometry_variables=self._stoichiometry_variables,
            variable_parameters=self._variable_parameters,
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
        self._variable_parameters = []
        self._constant_parameters = []
        self._stoichiometry_variables = []

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

        self._format_reactions()
        self._format_function_definitions()
        self._format_events()

        if not self._state_variables:
            raise NotImplementedError("Models without state variables are not supported.")

        self._populate_template_vars()
