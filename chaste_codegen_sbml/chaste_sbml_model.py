"""Module for the ChasteSbmlModel class."""

import abc
import os
import pathlib
import re
import subprocess
from typing import TYPE_CHECKING

from jinja2 import Environment, PackageLoader, select_autoescape
from libsbml import (
    AST_NAME,
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
    formulaToString,
)

from ._config import NON_DIM_UNITS, ROOT_DIR, EventType, ModelType, VarType
from ._utils import (
    generate_header_guard,
    get_function_definition_arguments,
    get_species_concentration,
    to_camel_case,
    to_cpp_name,
)

if TYPE_CHECKING:
    from typing import Any

    from jinja2.environment import Template
    from libsbml import ASTNode, Rule


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
        self._sbml_parameters = self._sbml_model.getListOfParameters()
        self._sbml_reactions = self._sbml_model.getListOfReactions()
        self._sbml_rules = self._sbml_model.getListOfRules()
        self._sbml_species = self._sbml_model.getListOfSpecies()
        self._sbml_unit_definitions = self._sbml_model.getListOfUnitDefinitions()

        self._variable_types = {}  # { id: VarType }
        self._odes = {}  # { id: str }

        self._assignment_rules = []  # [ { id: str, label: str, ... } ]
        self._state_variables = []  # [ { id: str, label: str, ... } ]
        self._derived_quantities = []  # [ { id: str, label: str, ... } ]
        self._variable_parameters = []  # [ { id: str, label: str, ... } ]
        self._constant_parameters = []  # [ { id: str, label: str, ... } ]
        self._rule_based_parameters = []  # [ { id: str, label: str, ... } ]
        self._reactions = []  # [ { id: str, label: str, ... } ]
        self._events = []  # [ { name: str, trigger: str, ... } ]
        self._functions = []  # [ { name: str, args: [str], body: str } ]

        self._template_vars = {}  # type: dict[str, Any]

        self._outputs = {}  # { filename: code }

        self._process_model()

    @property
    def outputs(self) -> dict[str, str]:
        return self._outputs

    def write(self, output_directory=None):
        """Generate Chaste code and write to file."""
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

    def _add_assignment_rule(self, id_: str, label: str, lhs: str, rhs: str) -> None:
        """Add an assignment rule to the template variables."""
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

    def _add_constant_parameter(self, id_: str, label: str, value: float, units: str) -> None:
        """Add a constant parameter to the template variables."""
        self._constant_parameters.append(
            {
                "id": id_,
                "label": label,
                "index": len(self._constant_parameters),
                "value": value,
                "units": units,
            }
        )
        self._variable_types[id_] = VarType.CONSTANT_PARAMETER

    def _add_derived_quantity(
        self, id_: str, label: str, initial_value: float, units: str, rhs: str
    ) -> None:
        """Add a derived quantity to the template variables."""
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

    def _add_output(self, filename: str, code: str) -> None:
        """Add generated code to the outputs dictionary.

        :param filename: The filename.
        :param code: The code.
        """
        self._outputs[filename] = code

    def _add_reaction(
        self, id_: str, label: str, rhs: str, parameters: list[dict[str, "Any"]]
    ) -> None:
        """Add a reaction to the template variables."""
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

    def _add_rule_based_parameter(
        self, id_: str, label: str, initial_value: float, units: str = NON_DIM_UNITS
    ) -> None:
        """Add a rule parameter to the template variables."""
        self._rule_based_parameters.append(
            {
                "id": id_,
                "label": label,
                "index": len(self._rule_based_parameters),
                "initial_value": initial_value,
                "units": units,
            }
        )
        self._variable_types[id_] = VarType.RULE_BASED_PARAMETER

    def _add_state_variable(
        self, id_: str, label: str, initial_value: float, units: str, rhs: str
    ) -> None:
        """Add a state variable to the template variables."""
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

    def _add_variable_parameter(
        self, id_: str, label: str, initial_value: float, units: str = NON_DIM_UNITS
    ) -> None:
        """Add a variable parameter to the template variables."""
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

    def _convert_ast_formula(self, ast_formula: "ASTNode") -> str:
        """Convert SBML AST formula to equivalent C++ string.

        :param ast: The AST formula.
        :param convert_names: Whether to convert state variable and parameter names.
        :return: The equivalent C++ string.
        """
        return self._convert_str_formula(formulaToString(ast_formula))

    def _convert_str_formula(self, formula: str) -> str:
        """Convert SBML string formula to equivalent C++ string.

        :param ast: The string formula.
        :return: The equivalent C++ string.
        """
        # Convert all integers to doubles
        # TODO: Instead of regex, traverse AST and convert AST_INTEGER nodes to AST_REAL
        # This has an adverse effect on literals like 4e-6 i.e. 4e-6.0 is invalid
        formula = re.sub(r"(?<!\.)\b[0-9]+\b(?!\.)", lambda x: f"{x[0]}.0", formula)

        # TODO: implies, lambda, delay

        # SBML contants to be replaced with C++ equivalents
        constants = {
            "avogadro": "sm::AVOGADRO",
            "exponentiale": "M_E",
            "inf": "std::numeric_limits<double>::infinity()",
            "infinity": "std::numeric_limits<double>::infinity()",
            "nan": "NAN",
            "notanumber": "NAN",
            "pi": "M_PI",
            "time": "SimulationTime::Instance()->GetTimeStep()",
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
        return cpp_formula

    def _extract_odes(self) -> None:
        """Extract the ODEs equations for each species."""
        self._odes = {}

        # Extract ODEs from reactions:
        # each ODE will essentially be the sum of the products minus the
        # sum of the reactants divided by the compartment volume
        for reaction in self._sbml_reactions:
            reaction_id = reaction.getId()

            # Decompose reaction into sum of products minus sum of reactants
            products = reaction.getListOfProducts()
            for product in products:
                species_id = product.getSpecies()

                stoichiometry = product.getStoichiometry()
                if float(stoichiometry) == 1.0:
                    rhs = f"{reaction_id}"
                else:
                    rhs = f"({stoichiometry} * {reaction_id})"

                if species_id in self._odes:
                    self._odes[species_id] += f" + {rhs}"
                else:
                    self._odes[species_id] = f"{rhs}"

            reactants = reaction.getListOfReactants()
            for reactant in reactants:
                species_id = reactant.getSpecies()

                stoichiometry = reactant.getStoichiometry()
                if float(stoichiometry) == 1.0:
                    rhs = f"{reaction_id}"
                else:
                    rhs = f"({stoichiometry} * {reaction_id})"

                if species_id in self._odes:
                    self._odes[species_id] += f" - {rhs}"
                else:
                    self._odes[species_id] = f"-{rhs}"

        # Extract ODEs from rate rules
        rate_rules = [r for r in self._sbml_rules if r.getTypeCode() == SBML_RATE_RULE]
        for rule in rate_rules:
            lhs = rule.getVariable()
            if lhs in self._odes:
                raise ValueError(f"{lhs} has both a rate rule and a reaction.")

            rhs = self._convert_str_formula(rule.getFormula())
            self._odes[lhs] = rhs

    def _format_compartments(self) -> None:
        """Add compartments to template variables."""
        for compartment in self._sbml_compartments:
            id_ = compartment.getId()
            label = compartment.getName().strip()
            value = compartment.getSize()
            self._add_variable_parameter(id_, label, value)

    def _format_events(self) -> None:
        """Add events to template variables."""
        # TODO: Add priority

        for event in self._sbml_events:
            label = event.getName().strip()

            # Try to guess the event type
            event_type = EventType.UNKNOWN

            cell_division_terms = ["cell division", "cytokinesis", "mitosis", "meiosis"]
            for term in cell_division_terms:
                if all(word in label.lower() for word in term.split()):
                    event_type = EventType.CELL_DIVISION
                    break

            trigger_math = event.getTrigger().getMath()
            trigger_formula = self._convert_ast_formula(trigger_math)

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
                lc = self._convert_ast_formula(trigger_math.getLeftChild())
                rc = self._convert_ast_formula(trigger_math.getRightChild())

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
                rhs = self._convert_ast_formula(assignment.getMath())

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
            body = self._convert_ast_formula(fd.getBody())

            self._add_function(
                fd_id,
                label,
                args,
                body,
            )

    def _format_parameters(self) -> None:
        """Add parameters to template variables."""
        # Note: rules must be processed before parameters
        if not self._assignment_rules:
            if any(r.getTypeCode() == SBML_ASSIGNMENT_RULE for r in self._sbml_rules):
                raise RuntimeError("Please process rules before parameters.")
        assignment_rules = {r["lhs"]: r["rhs"] for r in self._assignment_rules}

        for param in self._sbml_parameters:
            param_id = param.getId()
            label = param.getName().strip()

            value = param.getValue() if param.isSetValue() else 0.0
            units = param.getUnits() if param.isSetUnits() else NON_DIM_UNITS

            if param_id in self._odes:
                # State variable
                rhs = self._odes[param_id]
                self._add_state_variable(param_id, label, value, units, rhs)
            elif param_id in assignment_rules:
                # Rule-based parameter
                self._add_rule_based_parameter(param_id, label, value, units)
            elif param.isSetConstant() and not param.getConstant():
                # Variable parameter
                self._add_variable_parameter(param_id, label, value, units)
            else:
                # Constant parameter
                self._add_constant_parameter(param_id, label, value, units)

    def _format_reactions(self) -> None:
        """Add reactions to template variables."""
        for reaction in self._sbml_reactions:
            reaction_id = reaction.getId()
            label = reaction.getName().strip()

            kinetic_law = reaction.getKineticLaw()
            rhs = self._convert_str_formula(kinetic_law.getFormula())

            reaction_parameters = []
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
        # Sort assignment rules - variables on rhs must be defined before they are used
        assignment_rules = [r for r in self._sbml_rules if r.getTypeCode() == SBML_ASSIGNMENT_RULE]
        assignment_rules = self._sort_rules(assignment_rules)

        for rule in assignment_rules:
            rule_id = rule.getId()
            label = rule.getName().strip()
            lhs = rule.getVariable()
            rhs = self._convert_str_formula(rule.getFormula())
            self._add_assignment_rule(rule_id, label, lhs, rhs)

        # Algebraic rules are not implemented
        if any(r.getTypeCode() == SBML_ALGEBRAIC_RULE for r in self._sbml_rules):
            raise NotImplementedError("Algebraic rules are not yet supported.")

        # Note: Rate rules are handled during ODE extraction

    def _format_species(self) -> None:
        """Add species to template variables."""
        # Note: rules must be processed before species
        if not self._assignment_rules:
            if any(r.getTypeCode() == SBML_ASSIGNMENT_RULE for r in self._sbml_rules):
                raise RuntimeError("Please process rules before species.")
        assignment_rules = {r["lhs"]: r["rhs"] for r in self._assignment_rules}

        for species in self._sbml_species:
            species_id = species.getId()
            label = species.getName().strip()
            initial_value = get_species_concentration(species)
            boundary_condition = species.getBoundaryCondition()

            # If there's a compartment we'll normalise the ODEs, so declare it as non-dimensional
            compartment_id = species.getCompartment()
            compartment = self._sbml_compartments.get(compartment_id)
            units = NON_DIM_UNITS if compartment else species.getSubstanceUnits()

            if species_id in assignment_rules:
                # Derived quantity (includes boundary conditions with assignment rules)
                rhs = assignment_rules[species_id]
                self._add_derived_quantity(species_id, label, initial_value, units, rhs)

            elif boundary_condition:
                # Derived quantity (boundary condition with no assignment rule)
                rhs = str(initial_value)
                self._add_derived_quantity(species_id, label, initial_value, units, rhs)

            elif species_id in self._odes:
                # State variable
                # Normalised ODE
                rhs = self._odes[species_id]
                if len(rhs[1:].replace("+", "\t").replace("-", "\t").split("\t")) > 1:
                    # Add parentheses if there are multiple terms
                    rhs = f"({rhs})"
                rhs = f"{rhs} / {compartment_id}"
                self._add_state_variable(species_id, label, initial_value, units, rhs)

                # TODO: Handle time scaling
                # time_multiplier = self._get_timescale_multiplier()
                # if time_multiplier != 1.0:
                #     # This does not include species defined in algebraic rules
                #     f"rDY[{state_variable_index}] *= {time_multiplier};"

            else:
                # Variable parameter
                self._add_variable_parameter(species_id, label, initial_value, units)

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
        """Get the index of a variable."""
        var_type = self._get_variable_type(id_)

        if var_type == VarType.STATE_VARIABLE:
            for state_variable in self._state_variables:
                if state_variable["id"] == id_:
                    return state_variable["index"]

        elif var_type == VarType.DERIVED_QUANTITY:
            for dq in self._derived_quantities:
                if dq["id"] == id_:
                    return dq["index"]

        elif var_type == VarType.CONSTANT_PARAMETER:
            for param in self._constant_parameters:
                if param["id"] == id_:
                    return param["index"]

        elif var_type == VarType.VARIABLE_PARAMETER:
            for param in self._variable_parameters:
                if param["id"] == id_:
                    return param["index"]

        elif var_type == VarType.RULE_BASED_PARAMETER:
            for param in self._rule_based_parameters:
                if param["id"] == id_:
                    return param["index"]

        elif var_type == VarType.ASSIGNMENT_RULE:
            for rule in self._assignment_rules:
                if rule["id"] == id_:
                    return rule["index"]

        elif var_type == VarType.FUNCTION:
            for func in self._functions:
                if func["id"] == id_:
                    return func["index"]

        elif var_type == VarType.REACTION:
            for reaction in self._reactions:
                if reaction["id"] == id_:
                    return reaction["index"]

        raise ValueError(f"ID '{id_}' is not a recognized variable.")

    def _get_variable_type(self, var_id: str) -> bool:
        """Get the type of a variable based on its ID."""
        return self._variable_types.get(var_id, VarType.UNKNOWN)

    def _populate_template_vars(self) -> None:
        """Populate the template variables for generating C++ code."""
        template_vars: dict[str, "Any"] = dict(
            assignment_rules=self._assignment_rules,
            constant_parameters=self._constant_parameters,
            derived_quantities=self._derived_quantities,
            events=self._events,
            functions=self._functions,
            ode_class_name=self._ode_class_name,
            ode_header_guard=generate_header_guard(self._ode_hpp_filename),
            ode_hpp_file=self._ode_hpp_filename,
            reactions=self._reactions,
            rule_based_parameters=self._rule_based_parameters,
            state_variables=self._state_variables,
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
        self._state_variables = []
        self._derived_quantities = []
        self._variable_parameters = []
        self._constant_parameters = []
        self._rule_based_parameters = []

        self._reactions = []
        self._events = []
        self._functions = []

        self._extract_odes()

        self._format_rules()
        self._format_compartments()
        self._format_species()
        self._format_parameters()

        self._format_reactions()
        self._format_events()
        self._format_function_definitions()

        self._populate_template_vars()

    def _sort_rules(self, rules: list["Rule"]) -> list["Rule"]:
        """Sort rules based on their dependency.

        Rules are sorted such that if rule A depends on B (A -> B), then B comes
        before A. It is assumed that the input rules are acyclic. This function
        can't sort cyclic dependencies such as A -> B -> A, or A -> B -> C -> A.

        :param rules: The list of rules to sort.
        :return: The sorted list of rules.
        """

        def _search_formula(node: "ASTNode", name: str) -> bool:
            """Recursively search for a variable name in the AST formula."""
            if node is None:
                return False
            if node.getType() == AST_NAME and node.getName() == name:
                return True

            for i in range(node.getNumChildren()):
                child = node.getChild(i)
                if _search_formula(child, name):
                    return True
            return False

        _compare_cache = dict()

        def _compare_rules(rule_a: "Rule", rule_b: "Rule") -> int:
            """Compare two rules based on their dependency.

            :param rule_a: The first rule (A).
            :param rule_b: The second rule (B).

            :return: An integer indicating the order of the rules.
                -1 if A < B (A comes before B)
                1 if A > B (A comes after B)
                0 if the order doesn't matter
            """
            id_a = rule_a.getId()
            id_b = rule_b.getId()

            order = _compare_cache.get((id_a, id_b), None)
            if order is not None:
                return order

            var_a = rule_a.getVariable()
            var_b = rule_b.getVariable()

            rhs_a = rule_a.getMath()
            rhs_b = rule_b.getMath()

            # Check if var_a is in rhs_b
            if _search_formula(rhs_b, var_a):
                # var_a is used in rhs_b: rule_a comes before rule_b
                _compare_cache[(id_a, id_b)] = -1
                _compare_cache[(id_b, id_a)] = 1
                return -1

            # Check if var_b is in rhs_a
            if _search_formula(rhs_a, var_b):
                # var_b is used in rhs_a: rule_b comes before rule_a
                _compare_cache[(id_a, id_b)] = 1
                _compare_cache[(id_b, id_a)] = -1
                return 1

            # Order doesn't matter
            _compare_cache[(id_a, id_b)] = 0
            _compare_cache[(id_b, id_a)] = 0
            return 0

        def _insertion_sort(rules: list["Rule"]) -> list["Rule"]:
            """Sort rules using insertion sort based on their dependency."""
            # We need to compare each rule to all the others until we find a
            # non-zero comparison i.e. a +1 or -1 match (or until we exhaust
            # all options) because of cases such as:
            # Initial order: (a, b, c)
            # a == b (order doesn't matter, comparison is 0);
            # b == c (order doesn't matter, comparison is 0);
            # a > c (a should come after c, comparison is 1);
            # If we only compare (a, b) and (b, c), no changes will be made.
            sorted_rules = []

            for rule_a in rules:
                for i, rule_b in enumerate(sorted_rules):
                    if _compare_rules(rule_a, rule_b) < 0:  # rule_a < rule_b
                        sorted_rules.insert(i, rule_a)
                        break
                else:
                    # rule_a comes after everything already in sorted_rules
                    sorted_rules.append(rule_a)

            return sorted_rules

        return _insertion_sort(rules)
