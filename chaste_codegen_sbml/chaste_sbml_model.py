"""Module for the ChasteSbmlModel class."""

import abc
import os
import pathlib
import re
from typing import TYPE_CHECKING

from jinja2 import Environment, PackageLoader, select_autoescape
from libsbml import Parameter, SBMLReader, formulaToString

from ._config import ODE_SUFFIX, SHORT_NAME_LEN
from ._utils import (
    get_function_definition_arguments,
    get_species_concentration,
    varname_camelcase,
    varname_sanitize,
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

    SPECIAL_PARAMETER_NAMES = ["wnt", "gamma", "ComplexTransit"]

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
            self._model_name = varname_camelcase(filename).title()

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

        self._varnames = {}

        self._odes_dict = {}
        self._rules_dict = {}
        self._state_parameters = {}
        self._state_variables = {}

        self._update_odes_dict()
        self._update_rules_dict()
        self._update_state()

        self._num_state_vars = len(self._state_variables)

        self._outputs = {}  # { filename: code }

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

    def _add_state_parameter(self, obj_id: "str") -> None:
        """Add a state parameter to the model data.

        :param obj_id: The id of the Parameter / Species.
        """
        if obj_id not in self._state_parameters:
            index = len(self._state_parameters)
            self._state_parameters[obj_id] = index

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

    def _convert_ast_formula(self, ast_formula: "ASTNode", convert_names=True) -> str:
        """Convert SBML AST formula to equivalent C++ string.

        :param ast: The AST formula.
        :param convert_names: Whether to convert state variable and parameter names.
        :return: The equivalent C++ string.
        """
        return self._convert_str_formula(formulaToString(ast_formula), convert_names)

    def _convert_str_formula(self, formula: str, convert_names=True) -> str:
        """Convert SBML string formula to equivalent C++ string.

        :param ast: The string formula.
        :param convert_names: Whether to convert state variable and parameter names.
        :return: The equivalent C++ string.
        """

        # TODO: implies, lambda, delay

        # SBML contants to be replaced with C++ equivalents
        constants = {
            "avogadro": "sbmlmath::SM_AVOGADRO",
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
            "and": "sm_and",
            "acot": "sm_acot",
            "acoth": "sm_acoth",
            "acsc": "sm_acsc",
            "acsch": "sm_acsch",
            "asec": "sm_asec",
            "asech": "sm_asech",
            "arccot": "sm_acot",
            "arccoth": "sm_acoth",
            "arccsc": "sm_acsc",
            "arccsch": "sm_acsch",
            "arcsec": "sm_asec",
            "arcsech": "sm_asech",
            "cot": "sm_cot",
            "coth": "sm_coth",
            "csc": "sm_csc",
            "csch": "sm_csch",
            "eq": "sm_eq",
            "factorial": "sm_factorial",
            "geq": "sm_geq",
            "gt": "sm_gt",
            "leq": "sm_leq",
            "log": "sm_log",
            "lt": "sm_lt",
            "max": "sm_max",
            "min": "sm_min",
            "neq": "sm_neq",
            "not": "sm_not",
            "or": "sm_or",
            "piecewise": "sm_piecewise",
            "quotient": "sm_quotient",
            "root": "sm_root",
            "sec": "sm_sec",
            "sech": "sm_sech",
            "sqr": "sm_sqr",
            "xor": "sm_xor",
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
                cpp_token = f"sbmlmath::{custom_functions[token]}"

            # Replace state parameter and variable names.
            elif convert_names:
                if self._is_state_variable(token):
                    index = self._get_state_variable_index(token)
                    cpp_token = f"rY[{index}]"

                elif self._is_state_parameter(token):
                    index = self._get_state_parameter_index(token)
                    cpp_token = f"this->mParameters[{index}]"

            cpp_tokens.append(cpp_token)
        cpp_formula = "".join(cpp_tokens)
        return cpp_formula

    def _format_compartments(self) -> list[dict[str, "Any"]]:
        """Get a list of compartment dictionaries for the model.

        :return: A list of compartment dictionaries.
        """
        return [
            {
                "id": c.getId(),
                "name": c.getName(),
                "varname": self._get_varname(c),
                "size": c.getSize(),
            }
            for c in self._compartments
        ]

    def _format_events(self) -> list[dict[str, "Any"]]:
        """Get a list of event dictionaries for the model.

        :return: A list of event dictionaries.
        """
        event_dicts = []
        for event in self._events:
            trigger = event.getTrigger()
            trigger_formula = self._convert_ast_formula(trigger.getMath())

            assignment_formulas = []
            for assignment in event.getListOfEventAssignments():
                rhs_math = self._convert_ast_formula(assignment.getMath())
                rhs = f"static_cast<double>({rhs_math})"

                lhs = assignment.getVariable()
                if self._is_state_variable(lhs):
                    assignment_formulas.append(f'this->SetStateVariable("{lhs}", {rhs});')
                elif self._is_state_parameter(lhs):
                    assignment_formulas.append(f'this->SetParameter("{lhs}", {rhs});')
                else:
                    assignment_formulas.append(f"{lhs} = {rhs};")

            event_dicts.append(
                {
                    "trigger": trigger_formula,
                    "assignments": assignment_formulas,
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

    def _format_parameters(self) -> list[dict[str, "Any"]]:
        """Get a list of parameter dictionaries for the model.

        :return: A list of parameter dictionaries.
        """
        parameter_dicts = []
        for param in self._parameters:
            param_id = param.getId()
            name = param.getName()
            varname = self._get_varname(param)

            default = 0.0
            is_special = self._is_special_parameter(param)
            if is_special:
                # Special strings in the parameter name
                # TODO: Review how to handle these special parameters
                if any(x in name for x in ["gamma", "ComplexTransit"]):
                    default = 1.0
                # 0.0 for wnt and everything else
            elif param.isSetValue():
                default = param.getValue()

            value = param.getValue()  # TODO: if param.isSetValue() else p_default ?
            units = param.getUnits() if param.isSetUnits() else "non-dim"

            is_defined = (param.isSetValue() or param_id in self._rules_dict) and not is_special

            is_state_parameter = self._is_state_parameter(param_id)
            state_parameter_index = self._get_state_parameter_index(param_id)

            parameter_dicts.append(
                {
                    "default": default,
                    "id": param_id,
                    "is_defined": is_defined,
                    "is_state_parameter": is_state_parameter,
                    "name": name,
                    "state_parameter_index": state_parameter_index,
                    "units": units,
                    "value": value,
                    "varname": varname,
                }
            )
        return parameter_dicts

    def _format_reactions(self) -> list[dict[str, "Any"]]:
        """Get a list of reaction dictionaries for the model.

        :return: A list of reaction dictionaries.
        """
        reaction_dicts = []
        for reaction in self._reactions:
            reaction_id = reaction.getId()
            name = reaction.getName()
            varname = self._get_varname(reaction)

            kinetic_law = reaction.getKineticLaw()
            rhs = self._convert_str_formula(kinetic_law.getFormula())

            reaction_dict = {
                "id": reaction_id,
                "name": name,
                "parameters": [],
                "rhs": rhs,
                "varname": varname,
            }

            parameters = kinetic_law.getListOfParameters()
            for param in parameters:
                param_id = param.getId()
                param_name = self._get_name(param)
                param_varname = self._get_varname(param)
                param_value = param.getValue()  # TODO: What if not set?

                reaction_dict["parameters"].append(
                    {
                        "id": param_id,
                        "name": param_name,
                        "varname": param_varname,
                        "value": param_value,
                    }
                )

            reaction_dicts.append(reaction_dict)
        return reaction_dicts

    def _format_rules(self) -> list[dict[str, "Any"]]:
        """Get a list of rule dictionaries for the model.

        :return: A list of rule dictionaries.
        """
        return [
            {
                "id": r.getId(),
                "name": r.getName(),
                "formula": self._convert_str_formula(r.getFormula()),
            }
            for r in self._rules
        ]

    def _format_species(self) -> list[dict[str, "Any"]]:
        """Get a list of species dictionaries for the model.

        :return: A list of species dictionaries.
        """
        species_dicts = []
        # time_multiplier = self._get_timescale_multiplier()
        ode_index = 0

        for species in self._species:
            species_id = species.getId()
            name = species.getName()
            varname = self._get_varname(species)
            concentration = get_species_concentration(species)

            is_state_parameter = self._is_state_parameter(species_id)
            state_parameter_index = self._get_state_parameter_index(species_id)

            is_state_variable = self._is_state_variable(species_id)
            state_variable_index = self._get_state_variable_index(species_id)

            comp_id = species.getCompartment()
            comp = self._compartments.get(comp_id)
            comp_varname = self._get_varname(comp)

            # If there's a compartment we'll normalise the ODE so declare it as non-dimensional
            units = "non-dim" if comp else species.getSubstanceUnits()

            species_dict = {
                "concentration": concentration,
                "id": species_id,
                "is_state_parameter": is_state_parameter,
                "is_state_variable": is_state_variable,
                "name": name,
                "state_parameter_index": state_parameter_index,
                "state_variable_index": state_variable_index,
                "units": units,
                "varname": varname,
            }

            # ODE system
            lhs = None
            rhs = None
            if species_id in self._odes_dict:
                if not species.getBoundaryCondition():
                    index = ode_index
                    lhs = f"rDY[{index}]"
                    rhs = f"({self._odes_dict[species_id]}) / {comp_varname}"
                    ode_index += 1

                # TODO: Handle time scaling
                # if time_multiplier != 1.0:
                #     # This does not include species defined in algebraic rules
                #     f"rDY[{index}] *= {time_multiplier};"

            elif is_state_variable and (species_id == "drag" or varname == "drag"):
                index = ode_index
                lhs = f"rDY[{index}]"
                rhs = f"(drag - rY[{state_variable_index}]) / {comp_varname}"
                ode_index += 1

            elif species_id in self._rules_dict:
                # Species defined by algebraic rules are not in odes_dict
                # Assuming these are assignments where variables are added together to represent a total
                formula = self._rules_dict[species_id]
                tokens = formula.split(" ")
                rhs_tokens = []
                for token in tokens:
                    state_var_index = self._get_state_variable_index(token)
                    if state_var_index is not None:
                        rhs_tokens.append(f"rDY[{state_var_index}]")

                if rhs_tokens:
                    index = ode_index
                    lhs = f"rDY[{index}]"
                    rhs = f"({' + '.join(rhs_tokens)}) / {comp_varname}"
                    ode_index += 1

            # TODO: include other rules

            if lhs and rhs:
                species_dict["ode"] = {"lhs": lhs, "rhs": rhs}

            species_dicts.append(species_dict)
        return species_dicts

    def _get_name(self, obj: "SBase") -> str:
        """Get the name of a libSBML object, or the ID if it doesn't have one.

        :param obj: The object.
        :return: The object name, or ID.
        """
        obj_name = obj.getName().strip()
        if obj_name:
            return obj_name
        return obj.getId()

    def _get_state_parameter_index(self, obj_id: "str") -> int:
        """Get the index of a state parameter.

        :param obj_id: The id of the Parameter / Species.
        :return: The state parameter index.
        """
        return self._state_parameters.get(obj_id)

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

    def _get_template_vars(self, hpp_filename: str) -> dict[str, str]:
        """Generate the template variables for the model cpp file.

        :param hpp_filename: The hpp filename for the model.
        return: The template variables for the model cpp file.
        """
        return dict(
            compartments=self._format_compartments(),
            events=self._format_events(),
            function_definitions=self._format_function_definitions(),
            header_guard=self._format_header_guard(hpp_filename),
            model_hpp_file=hpp_filename,
            model_class_name=self._model_class_name,
            ode_class_name=self._ode_class_name,
            parameters=self._format_parameters(),
            reactions=self._format_reactions(),
            rules=self._format_rules(),
            num_state_vars=self._num_state_vars,
            species=self._format_species(),
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

    def _is_special_parameter(self, parameter: "Parameter") -> bool:
        """Check if a parameter has a special name.

        :param parameter: The Parameter.
        :return: True if special, False otherwise.
        """
        name = parameter.getName()
        return name and any(s in name for s in self.SPECIAL_PARAMETER_NAMES)

    def _is_state_parameter(self, obj_id: str) -> bool:
        """Check if a Species or Parameter is defined as a state parameter for Chaste.

        :param obj_id: The Species or Parameter ID to check.
        :return: True if defined as a state parameter, False otherwise.
        """
        return obj_id in self._state_parameters

    def _is_state_variable(self, species_id: str) -> bool:
        """Check if a Species is defined as a state variable for Chaste.

        :param species_id: The Species ID to check.
        :return: True if defined as a state variable, False otherwise.
        """
        return species_id in self._state_variables

    def _update_state(self) -> None:
        """Set the state parameters and state variables for the model."""
        self._state_parameters = {}
        self._state_variables = {}

        for species in self._species:
            species_id = species.getId()
            if (species_id in self._odes_dict) or (species_id in self._rules_dict):
                # Any species defined by an ODE or rule is a state variable
                self._add_state_variable(species_id)
            else:
                # All other species are defined as state parameters
                self._add_state_parameter(species_id)

        for param in self._parameters:
            param_id = param.getId()
            if self._is_special_parameter(param):
                # Parameters with special strings in their name are defined as state parameters.
                self._add_state_parameter(param_id)

            elif not (param.isSetValue() or (param_id in self._rules_dict)):
                # Parameters with unset values are also defined as state parameters.
                self._add_state_parameter(param_id)

    def _update_odes_dict(self) -> None:
        """Set the ODEs dictionary of equations corresponding to each species.

        Each ODE will essentially be the sum of the products minus the sum of
        the reactants divided by the compartment volume
        """
        self._odes_dict = {}
        for reaction in self._reactions:
            reaction_var = self._get_varname(reaction)

            # Decompose reaction into sum of products minus sum of reactants
            products = reaction.getListOfProducts()
            for product in products:
                # Get the species concerning the product
                species_id = product.getSpecies()

                # TODO: Do we need to do something special with boundary conditions?
                # species = self._model.getSpecies(species_id)
                # if species.isSetBoundaryCondition() and not species.getBoundaryCondition():

                if species_id in self._odes_dict:
                    self._odes_dict[species_id] += " + " + reaction_var
                else:
                    self._odes_dict[species_id] = reaction_var

            reactants = reaction.getListOfReactants()
            for reactant in reactants:
                species_id = reactant.getSpecies()

                # TODO: Do we need to do something special with boundary conditions?
                # species = self._model.getSpecies(species_id)
                # if species.isSetBoundaryCondition() and not species.getBoundaryCondition():

                if species_id in self._odes_dict:
                    self._odes_dict[species_id] += " - " + reaction_var
                else:
                    self._odes_dict[species_id] = "-" + reaction_var

    def _update_rules_dict(self) -> None:
        """Set the dictionary of species defined by reaction rules."""
        self._rules_dict = {
            r.getId(): self._convert_str_formula(r.getFormula()) for r in self._rules
        }
