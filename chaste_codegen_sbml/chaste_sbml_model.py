"""Module for the ChasteSbmlModel class."""

import abc
import os
import pathlib
import re
from typing import TYPE_CHECKING

from jinja2 import Environment, PackageLoader, select_autoescape
from libsbml import (
    AST_RELATIONAL_EQ,
    AST_RELATIONAL_GEQ,
    AST_RELATIONAL_GT,
    AST_RELATIONAL_LEQ,
    AST_RELATIONAL_LT,
    AST_RELATIONAL_NEQ,
    SBML_ASSIGNMENT_RULE,
    SBML_ALGEBRAIC_RULE,
    SBML_RATE_RULE,
    SBMLReader,
    formulaToString,
)

from ._config import ODE_SUFFIX, SHORT_NAME_LEN
from ._utils import (
    get_function_definition_arguments,
    get_species_concentration,
    varname_sanitize,
    varname_staggercase,
)

if TYPE_CHECKING:
    from typing import Any

    from jinja2.environment import Template
    from libsbml import ASTNode, SBase


class ChasteSbmlModel:
    """Holds information about an SBML model for which Chaste code is to be generated."""

    __metaclass__ = abc.ABCMeta

    _jinja_env = Environment(
        loader=PackageLoader("chaste_codegen_sbml"),
        autoescape=select_autoescape(),
        trim_blocks=True,
        lstrip_blocks=True,
    )

    STATE_VARIABLE = 0
    DERIVED_QUANTITY = 1
    VARIABLE_PARAMETER = 2
    CONSTANT_PARAMETER = 3

    # -- PUBLIC ---------------------------------------

    def __init__(self, sbml_file: str, model_name: str = None, model_suffix: str = None) -> None:
        """Initialise the ChasteSbmlModel.

        :param sbml: The SBML file.
        :param name: The model name.
        :param suffix: The model type suffix e.g. "CellCycle".
        """
        self._sbml_file = os.path.abspath(sbml_file)
        if not os.path.isfile(self._sbml_file):
            raise FileNotFoundError(f"Could not find SBML file: {self._sbml_file}")

        if model_name:
            self._model_name = model_name
        else:
            filename = os.path.splitext(os.path.basename(self._sbml_file))[0]
            model_name = varname_staggercase(filename) + "Sbml"
            self._model_name = model_name[0].upper() + model_name[1:]

        self._model_suffix = model_suffix
        self._ode_class_name = self._model_name + ODE_SUFFIX
        self._model_class_name = f"{self._model_name}{self._model_suffix}Model"
        self._wrapper_class_name = f"Sbml{self._model_suffix}WrapperModel"

        self._model = SBMLReader().readSBMLFromFile(self._sbml_file).getModel()
        self._compartments = self._model.getListOfCompartments()
        self._events = self._model.getListOfEvents()
        self._function_definitions = self._model.getListOfFunctionDefinitions()
        self._parameters = self._model.getListOfParameters()
        self._reactions = self._model.getListOfReactions()
        self._rules = self._model.getListOfRules()
        self._species = self._model.getListOfSpecies()
        self._unit_definitions = self._model.getListOfUnitDefinitions()

        self._rule_eqs = None

        self._state_variables = None
        self._derived_quantities = None
        self._variable_parameters = None
        self._constant_parameters = None

        self._outputs = {}  # { filename: code }

        self._process_model()

    def write(self, output_directory=None):
        """Generate Chaste code and write to file."""
        # Generate the code
        self._generate()

        # Write the code to file
        if output_directory:
            root_dir = pathlib.Path(output_directory)
        else:
            root_dir = pathlib.Path().cwd()

        for filename, code in self._outputs.items():
            file_path = root_dir / filename
            with open(file_path, "w") as f:
                f.write(code)

    # -- PRIVATE ---------------------------------------

    @abc.abstractmethod
    def _generate(self) -> None:
        """Generate Chaste code for the model.
        This method should be implemented by subclasses.
        Generated code should be stored in the outputs using _add_output.
        """
        return

    def _add_constant_parameter(self, obj_id: "str") -> None:
        """Add a constant parameter to the model data.

        :param obj_id: The id of the Parameter.
        """
        if obj_id not in self._constant_parameters:
            index = len(self._constant_parameters)
            self._constant_parameters[obj_id] = index

    def _add_derived_quantity(self, obj_id: "str") -> None:
        """Add a derived quantity to the model data.

        :param obj_id: The id of the Species.
        """
        if obj_id not in self._derived_quantities:
            index = len(self._derived_quantities)
            self._derived_quantities[obj_id] = index

    def _add_variable_parameter(self, obj_id: "str") -> None:
        """Add a variable parameter to the model data.

        :param obj_id: The id of the Parameter / Species.
        """
        if obj_id not in self._variable_parameters:
            index = len(self._variable_parameters)
            self._variable_parameters[obj_id] = index

    def _add_state_variable(self, obj_id: "str") -> None:
        """Add a state variable to the model data.

        :param obj_id: The id of the Species.
        """
        if obj_id not in self._state_variables:
            index = len(self._state_variables)
            self._state_variables[obj_id] = index

    def _add_output(self, filename: str, code: str) -> None:
        """Add generated code to the outputs dictionary.

        :param filename: The filename.
        :param code: The code.
        """
        self._outputs[filename] = code

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

    def _format_compartments(self, variable_parameters: list[dict[str, "Any"]]) -> None:
        """Add compartments to template variables."""
        for compartment in self._compartments:
            variable_parameters.append(
                {
                    "id": compartment.getId(),
                    "index": len(variable_parameters),
                    "name": compartment.getName(),
                    "value": compartment.getSize(),
                }
            )

    def _format_events(self) -> list[dict[str, "Any"]]:
        """Get a list of event dictionaries for the model.

        :return: A list of event dictionaries.
        """
        event_dicts = []
        for event in self._events:
            name = event.getName()
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

            assignment_formulas = []
            for assignment in event.getListOfEventAssignments():
                rhs = self._convert_ast_formula(assignment.getMath())

                lhs = assignment.getVariable()
                if self._is_state_variable(lhs):
                    state_variable_index = self._get_state_variable_index(lhs)
                    assignment_formulas.append(f"SetStateVariable({state_variable_index}, {rhs})")
                    assignment_formulas.append(
                        f"SetDefaultInitialCondition({state_variable_index}, {rhs})"
                    )
                elif self._is_variable_parameter(lhs):
                    assignment_formulas.append(f'SetParameter("{lhs}", {rhs})')
                else:
                    assignment_formulas.append(f"{lhs} = {rhs}")

            event_dicts.append(
                {
                    "assignments": assignment_formulas,
                    "distance": trigger_distance,
                    "name": name,
                    "trigger": trigger_formula,
                }
            )
        return event_dicts

    def _format_function_definitions(self) -> list[dict[str, "Any"]]:
        """Get a list of function definition dictionaries for the model.

        :return: A list of function definition dictionaries.
        """
        function_definition_dicts = []
        for fd in self._function_definitions:
            fd_id = fd.getId()
            arg_list = get_function_definition_arguments(fd)
            args = ", ".join(map(lambda x: f"double {x}", arg_list))
            body = self._convert_ast_formula(fd.getBody())

            function_definition_dicts.append(
                {
                    "id": fd_id,
                    "args": args,
                    "body": body,
                }
            )
        return function_definition_dicts

    def _format_header_guard(self, filename: str) -> str:
        """Get the header guard for a file.

        :param filename: The filename.
        :return: The header guard.
        """
        return filename.upper().replace(".", "_") + "_"

    def _format_parameters(self, constant_parameters, variable_parameters) -> None:
        """Add parameters to template variables."""

        for param in self._parameters:
            param_id = param.getId()
            name = param.getName()

            value = param.getValue() if param.isSetValue() else 0.0
            units = param.getUnits() if param.isSetUnits() else "non-dim"

            is_variable = param.isSetConstant() and not param.getConstant()
            index = len(variable_parameters) if is_variable else len(constant_parameters)

            param_dict = {
                "id": param_id,
                "index": index,
                "name": name,
                "units": units,
                "value": value,
            }

            if is_variable:
                variable_parameters.append(param_dict)
            else:
                constant_parameters.append(param_dict)

    def _format_reactions(self) -> list[dict[str, "Any"]]:
        """Get a list of reaction dictionaries for the model.

        :return: A list of reaction dictionaries.
        """
        reaction_dicts = []
        for reaction in self._reactions:
            reaction_id = reaction.getId()
            name = reaction.getName()

            kinetic_law = reaction.getKineticLaw()
            rhs = self._convert_str_formula(kinetic_law.getFormula())

            reaction_dict = {
                "id": reaction_id,
                "name": name,
                "parameters": [],
                "rhs": rhs,
            }

            parameters = kinetic_law.getListOfParameters()
            for param in parameters:
                param_id = param.getId()
                param_name = self._get_name(param)
                param_value = param.getValue()  # TODO: What if not set?

                reaction_dict["parameters"].append(
                    {
                        "id": param_id,
                        "name": param_name,
                        "value": param_value,
                    }
                )

            reaction_dicts.append(reaction_dict)
        return reaction_dicts

    def _format_rules(self) -> list[dict[str, "Any"]]:
        """Get a list of rule dictionaries for the model.

        :return: A list of rule dictionaries.
        """
        rule_dicts = []

        for rule in self._rules:
            rule_id = rule.getId()
            name = rule.getName()
            lhs = rule.getVariable()
            rhs = self._convert_str_formula(rule.getFormula())
            is_state_variable = self._is_state_variable(lhs)
            is_variable_parameter = self._is_variable_parameter(lhs)

            rule_dict = {
                "id": rule_id,
                "name": name,
                "rhs": rhs,
                "lhs": lhs,
                "is_variable_parameter": is_variable_parameter,
                "is_state_variable": is_state_variable,
            }
            rule_dicts.append(rule_dict)
        return rule_dicts

    def _format_species(self, state_variables, derived_quantities, variable_parameters) -> None:
        """Add species to template variables."""
        """Get a list of species dictionaries for the model.

        :return: A list of species dictionaries.
        """
        odes = self._process_odes()

        for species in self._species:
            species_id = species.getId()
            name = species.getName()
            initial_value = get_species_concentration(species)

            comp_id = species.getCompartment()
            compartment = self._compartments.get(comp_id)
            # If there's a compartment we'll normalise the ODE, so declare it as non-dimensional
            units = "non-dim" if compartment else species.getSubstanceUnits()

            index = None
            rhs = None
            if species_id in odes:
                # State variable
                index = len(state_variables)
                rhs = f"({odes[species_id]}) / {comp_id}"

                # TODO: Handle time scaling
                # time_multiplier = self._get_timescale_multiplier()
                # if time_multiplier != 1.0:
                #     # This does not include species defined in algebraic rules
                #     f"rDY[{state_variable_index}] *= {time_multiplier};"

                # TODO: Do something different for boundary conditions?
                # if species.getBoundaryCondition():

            elif species_id in self._rules_dict:
                # Derived quantity
                index = len(derived_quantities)
                rhs = self._rules_dict[species_id]

            else:
                # Variable parameter
                index = len(variable_parameters)

            species_dict = {
                "id": species_id,
                "index": index,
                "initial_value": initial_value,
                "name": name,
                "rhs": rhs,
                "units": units,
            }

    def _get_name(self, obj: "SBase") -> str:
        """Get the name of a libSBML object, or the ID if it doesn't have one.

        :param obj: The object.
        :return: The object name, or ID.
        """
        obj_name = obj.getName().strip()
        if obj_name:
            return obj_name
        return obj.getId()

    def _get_constant_parameter_index(self, obj_id: "str") -> int:
        """Get the index of a constant parameter.

        :param obj_id: The id of the Parameter.
        :return: The constant parameter index.
        """
        return self._constant_parameters.get(obj_id)

    def _get_derived_quantity_index(self, obj_id: "str") -> int:
        """Get the index of a derived quantity.

        :param obj_id: The id of the Species.
        :return: The derived quantity index.
        """
        return self._derived_quantities.get(obj_id)

    def _get_variable_parameter_index(self, obj_id: "str") -> int:
        """Get the index of a variable parameter.

        :param obj_id: The id of the Parameter / Species.
        :return: The variable parameter index.
        """
        return self._variable_parameters.get(obj_id)

    def _get_state_variable_index(self, obj_id: "str") -> int:
        """Get the index of a state variable.

        :param obj_id: The id of the Species.
        :return: The state variable index.
        """
        return self._state_variables.get(obj_id)

    def _get_template(self, name: str) -> "Template":
        """Get a Jinja2 template.

        :param name: The template name.
        :return: The template object.
        """
        return self._jinja_env.get_template(name)

    def _get_template_vars(self, hpp_filename: str) -> dict[str, "Any"]:
        """Generate the template variables for the model C++ files.

        :param hpp_filename: The hpp filename for the model.
        :return: The template variables for the model C++ files.
        """
        return dict(
            compartments=self._formatted_compartments,
            events=self._formatted_events,
            function_definitions=self._formatted_function_definitions,
            header_guard=self._format_header_guard(hpp_filename),
            model_hpp_file=hpp_filename,
            model_class_name=self._model_class_name,
            ode_class_name=self._ode_class_name,
            parameters=self._formatted_parameters,
            reactions=self._formatted_reactions,
            rules=self._formatted_rules,
            num_state_vars=self._num_state_vars,
            species=self._formatted_species,
            wrapper_class_name=self._wrapper_class_name,
        )

    def _get_timescale_multiplier(self) -> float:
        """Get the timescale multiplier.

        SBML uses seconds by default and Chaste uses hours.

        :return: The timescale multiplier.
        """
        for unit_def in self._unit_definitions:
            u_id = unit_def.getId()
            if u_id.lower() == "time":  # Do people ever call this something different?
                timescale = unit_def.getName().lower()
                if "minute" in timescale:
                    return 60.0
                elif "hour" in timescale:
                    return 1.0
        return 3600.0

    def _get_varname(self, obj: "SBase") -> str:
        """Get a suitable C++ variable name for a libSBML object.

        :param obj: The object.
        :return: The variable name.
        """
        # TODO: id is not unique for all objects e.g. reaction parameters
        obj_id = obj.getId()
        if obj_id in self._varnames:
            return self._varnames[obj_id]

        # Prefer the name if it is reasonably short, or shorter than the ID
        obj_name = varname_sanitize(self._get_name(obj))
        if 0 < len(obj_name) <= max(SHORT_NAME_LEN, len(obj_id)):
            var = obj_name
        else:
            var = obj_id

        # Check that all generated variable names are unique
        if var in self._varnames.values():
            i = 0
            while f"{var}_{i}" in self._varnames.values():
                i += 1
            var = f"{var}_{i}"

        self._varnames[obj_id] = var
        return var

    def _is_parameter(self, obj_id: str) -> bool:
        """Check if ID belongs to a parameter.

        :param obj_id: The ID to check.
        :return: True if the ID belongs to a parameter.
        """
        return any(obj_id == p.getId() for p in self._parameters)

    def _is_species(self, obj_id: str) -> bool:
        """Check if ID belongs to a species.

        :param obj_id: The ID to check.
        :return: True if the ID belongs to a species.
        """
        return any(obj_id == p.getId() for p in self._species)

    def _is_constant_parameter(self, obj_id: str) -> bool:
        """Check if ID belongs to a constant parameter.

        :param obj_id: The ID to check.
        :return: True if the ID belongs to a constant parameter.
        """
        return obj_id in self._constant_parameters

    def _is_derived_quantity(self, obj_id: str) -> bool:
        """Check if a Species is defined as a derived quantity for Chaste.

        :param obj_id: The Species ID to check.
        :return: True if defined as a derived quantity, False otherwise.
        """
        return obj_id in self._derived_quantities

    def _is_variable_parameter(self, obj_id: str) -> bool:
        """Check if a Species or Parameter is defined as a variable parameter for Chaste.

        :param obj_id: The Species or Parameter ID to check.
        :return: True if defined as a variable parameter, False otherwise.
        """
        return obj_id in self._variable_parameters

    def _is_state_variable(self, species_id: str) -> bool:
        """Check if a Species is defined as a state variable for Chaste.

        :param species_id: The Species ID to check.
        :return: True if defined as a state variable, False otherwise.
        """
        return species_id in self._state_variables

    def _process_model(self) -> None:
        """Process the SBML model to set up the ODEs and rules dictionaries."""

        state_variables = []
        derived_quantities = []
        variable_parameters = []
        constant_parameters = []

        self._format_compartments(variable_parameters)
        self._format_species(state_variables, derived_quantities, variable_parameters)
        self._format_parameters(constant_parameters, variable_parameters)

        # Process the rules dictionary
        self._process_rules()

        # Process variable types
        parameter_types = {}
        for parameter in variable_parameters:
            _id = parameter["id"]
            parameter_types[_id] = VARIABLE_PARAMETER

        for state_variable in state_variables:
            _id = state_variable["id"]
            parameter_types[_id] = STATE_VARIABLE

        for derived_quantity in derived_quantities:
            _id = derived_quantity["id"]
            parameter_types[_id] = DERIVED_QUANTITY

        for constant_parameter in constant_parameters:
            _id = constant_parameter["id"]
            parameter_types[_id] = CONSTANT_PARAMETER

        # Set the number of state variables
        self._num_state_vars = len(self._state_variables)

    def _process_odes(self) -> dict[str, str]:
        """Get the ODEs equations for each species.

        Each ODE will essentially be the sum of the products minus the sum of
        the reactants divided by the compartment volume
        """
        odes = {}
        for reaction in self._reactions:
            reaction_id = reaction.getId()

            # Decompose reaction into sum of products minus sum of reactants
            products = reaction.getListOfProducts()
            for product in products:
                # Get the species concerning the product
                species_id = product.getSpecies()

                # TODO: Do we need to do something special with boundary conditions?
                # species = self._model.getSpecies(species_id)
                # if species.isSetBoundaryCondition() and not species.getBoundaryCondition():

                if species_id in odes:
                    odes[species_id] += " + " + reaction_id
                else:
                    odes[species_id] = reaction_id

            reactants = reaction.getListOfReactants()
            for reactant in reactants:
                species_id = reactant.getSpecies()

                # TODO: Do we need to do something special with boundary conditions?
                # species = self._model.getSpecies(species_id)
                # if species.isSetBoundaryCondition() and not species.getBoundaryCondition():

                if species_id in odes:
                    odes[species_id] += " - " + reaction_id
                else:
                    odes[species_id] = "-" + reaction_id

            return odes

    def _process_rules(self) -> None:
        """Get the rules for each species."""
        self._rules = {}
        for rule in self._rules:
            type_code = rule.getTypeCode()
            if type_code == SBML_ASSIGNMENT_RULE:
                lhs = rule.getVariable()
                rhs = self._convert_str_formula(rule.getFormula())
                self._rules[lhs] = rhs
            elif type_code == SBML_ALGEBRAIC_RULE:
                # Not implemented
                raise NotImplementedError("Algebraic rules are not yet supported.")
            elif type_code == SBML_RATE_RULE:
                # Not implemented
                raise NotImplementedError("Rate rules are not yet supported.")

    def _process_variable_types(self) -> None:
        """Set the variable types for the model."""
        self._state_variables = {}
        self._derived_quantities = {}
        self._variable_parameters = {}
        self._constant_parameters = {}

        # Species defined by rate rules (ODEs) are state variables.
        # Species defined by assignment rules are derived quantities.
        # All other species are variable parameters.
        for species in self._species:
            species_id = species.getId()
            if species_id in self._odes_dict:
                self._add_state_variable(species_id)
            elif species_id in self._rules_dict:
                self._add_derived_quantity(species_id)
            else:
                self._add_variable_parameter(species_id)

        # Non-constant parameters are variable parameters.
        for param in self._parameters:
            if param.isSetConstant() and not param.getConstant():
                self._add_variable_parameter(param.getId())
            else:
                self._add_constant_parameter(param.getId())
