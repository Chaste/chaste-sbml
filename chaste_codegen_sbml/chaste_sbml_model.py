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
    SBML_ALGEBRAIC_RULE,
    SBML_ASSIGNMENT_RULE,
    SBML_RATE_RULE,
    SBMLReader,
    formulaToString,
)

from ._config import NON_DIM_UNITS, ODE_SUFFIX, VarType
from ._utils import (
    get_function_definition_arguments,
    get_species_concentration,
    varname_staggercase,
)

if TYPE_CHECKING:
    from typing import Any

    from jinja2.environment import Template
    from libsbml import ASTNode


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

    # -- PUBLIC --------------------------------------

    def __init__(self, sbml_file: str, model_name: str = "", model_suffix: str = "") -> None:
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

        self._sbml_model = SBMLReader().readSBMLFromFile(self._sbml_file).getModel()
        self._sbml_compartments = self._sbml_model.getListOfCompartments()
        self._sbml_events = self._sbml_model.getListOfEvents()
        self._sbml_function_definitions = self._sbml_model.getListOfFunctionDefinitions()
        self._sbml_parameters = self._sbml_model.getListOfParameters()
        self._sbml_reactions = self._sbml_model.getListOfReactions()
        self._sbml_rules = self._sbml_model.getListOfRules()
        self._sbml_species = self._sbml_model.getListOfSpecies()
        self._sbml_unit_definitions = self._sbml_model.getListOfUnitDefinitions()

        self._assignment_rules = []  # [ { id: str, label: str, ... } ]
        self._state_variables = []  # [ { id: str, label: str, ... } ]
        self._derived_quantities = []  # [ { id: str, label: str, ... } ]
        self._variable_parameters = []  # [ { id: str, label: str, ... } ]
        self._constant_parameters = []  # [ { id: str, label: str, ... } ]
        self._rule_based_parameters = []  # [ { id: str, label: str, ... } ]
        self._reactions = []  # [ { id: str, label: str, ... } ]
        self._events = []  # [ { name: str, trigger: str, ... } ]
        self._functions = []  # [ { name: str, args: [str], body: str } ]

        self._variable_types = {}  # { id: VarType }

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
        Generated code should be stored in the outputs using `_add_output`.

        Example:
        cpp_code = ...
        hpp_code = ...
        self._add_output(cpp_filename, cpp_code)
        self._add_output(hpp_filename, hpp_code)
        """
        return

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

    def _add_derived_quantity(self, id_: str, label: str, units: str, rhs: str) -> None:
        """Add a derived quantity to the template variables."""
        self._derived_quantities.append(
            {
                "id": id_,
                "label": label,
                "index": len(self._derived_quantities),
                "rhs": rhs,
                "units": units,
            }
        )
        self._variable_types[id_] = VarType.DERIVED_QUANTITY

    def _add_event(
        self, label: str, trigger: str, assignments: list[dict[str, "Any"]], distance: str
    ) -> None:
        """Add an event to the template variables.

        :param label: The event description.
        :param trigger: The event trigger formula.
        :param assignments: The event assignments.
        :param distance: The distance for the event trigger.
        """
        self._events.append(
            {
                "label": label,
                "index": len(self._events),
                "trigger": trigger,
                "assignments": assignments,
                "distance": distance,
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

            self._add_event(label, trigger_formula, assignments, trigger_distance)

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

    def _format_header_guard(self, filename: str) -> str:
        """Get the header guard for a file.

        :param filename: The filename.
        :return: The header guard.
        """
        return filename.upper().replace(".", "_") + "_"

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

            if param_id in assignment_rules:
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
        for rule in self._sbml_rules:
            rule_id = rule.getId()
            label = rule.getName().strip()

            type_code = rule.getTypeCode()
            if type_code == SBML_ASSIGNMENT_RULE:
                lhs = rule.getVariable()
                rhs = self._convert_str_formula(rule.getFormula())
                self._add_assignment_rule(rule_id, label, lhs, rhs)

            elif type_code == SBML_ALGEBRAIC_RULE:
                # Not implemented
                raise NotImplementedError("Algebraic rules are not yet supported.")

            elif type_code == SBML_RATE_RULE:
                # Not implemented
                raise NotImplementedError("Rate rules are not yet supported.")

    def _format_species(self) -> None:
        """Add species to template variables."""
        odes = self._get_odes()

        # Note: rules must be processed before species
        if not self._assignment_rules:
            if any(r.getTypeCode() == SBML_ASSIGNMENT_RULE for r in self._sbml_rules):
                raise RuntimeError("Please process rules before species.")
        assignment_rules = {r["lhs"]: r["rhs"] for r in self._assignment_rules}

        for species in self._sbml_species:
            species_id = species.getId()
            label = species.getName().strip()
            initial_value = get_species_concentration(species)

            # If there's a compartment we'll normalise the ODEs, so declare it as non-dimensional
            compartment_id = species.getCompartment()
            compartment = self._sbml_compartments.get(compartment_id)
            units = NON_DIM_UNITS if compartment else species.getSubstanceUnits()

            if species_id in assignment_rules:
                # Derived quantity (includes boundary conditions)
                rhs = assignment_rules[species_id]
                self._add_derived_quantity(species_id, label, units, rhs)

            elif species_id in odes:
                # State variable
                rhs = f"({odes[species_id]}) / {compartment_id}"  # Normalised ODE
                self._add_state_variable(species_id, label, initial_value, units, rhs)

                # TODO: Handle time scaling
                # time_multiplier = self._get_timescale_multiplier()
                # if time_multiplier != 1.0:
                #     # This does not include species defined in algebraic rules
                #     f"rDY[{state_variable_index}] *= {time_multiplier};"

            else:
                # Variable parameter
                self._add_variable_parameter(species_id, label, initial_value, units)

    def _get_odes(self) -> dict[str, str]:
        """Get the ODEs equations for each species.

        Each ODE will essentially be the sum of the products minus the sum of
        the reactants divided by the compartment volume
        """
        odes = {}
        for reaction in self._sbml_reactions:
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
            assignment_rules=self._assignment_rules,
            constant_parameters=self._constant_parameters,
            derived_quantities=self._derived_quantities,
            events=self._events,
            functions=self._functions,
            header_guard=self._format_header_guard(hpp_filename),
            model_class_name=self._model_class_name,
            model_hpp_file=hpp_filename,
            ode_class_name=self._ode_class_name,
            reactions=self._reactions,
            rule_based_parameters=self._rule_based_parameters,
            state_variables=self._state_variables,
            variable_parameters=self._variable_parameters,
            wrapper_class_name=self._wrapper_class_name,
        )

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
        """Get the index of a variable"""
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

    def _process_model(self) -> None:
        """Process the SBML model to set up the formatted variables for templates."""

        self._assignment_rules = []

        self._state_variables = []
        self._derived_quantities = []
        self._variable_parameters = []
        self._constant_parameters = []
        self._rule_based_parameters = []

        self._reactions = []
        self._events = []
        self._functions = []

        self._variable_types = {}

        self._format_rules()
        self._format_compartments()
        self._format_species()
        self._format_parameters()

        self._format_reactions()
        self._format_events()
        self._format_function_definitions()
