import os
import pathlib
from typing import TYPE_CHECKING

from libsbml import SBMLReader, formulaToString

from ._utils import varname_camel, varname_sanitize

if TYPE_CHECKING:
    from libsbml import (
        ASTNode,
        FunctionDefinition,
        Parameter,
        SBase,
        Species,
    )


class ChasteModel:
    """Holds information about the SBML model for which code is to be generated."""

    def __init__(self, sbml_file: str, model_name: str = None) -> None:
        self.sbml_file = sbml_file
        self.sbml_model = SBMLReader().readSBMLFromFile(sbml_file).getModel()

        if model_name:
            self.model_name = model_name
        else:
            filename = os.path.splitext(os.path.basename(sbml_file))[0]
            self.model_name = varname_camel(filename)

        self.ode_system_name = f"{self.model_name}OdeSystem"

        self.hpp_filename = f"{self.model_name}.hpp"
        self.cpp_filename = f"{self.model_name}.cpp"

        self.hpp_source = ""
        self.cpp_source = ""

        self.odes_dict = self.get_ode_dict()
        self.rules_dict = self.get_rules_dict()

    @staticmethod
    def convert_formula(formula: str) -> str:
        """Convert a formula to its C++ equivalent.

        :param formula: The formula.
        :return: The C++ equivalent of the formula.
        """

        method_map = {
            "abs": "fabs",
            "arccos": "acos",
            "arccsc": "acsc",
            "arccsch": "acsch",
            "arcsec": "asec",
            "arcsech": "asech",
            "arcsinh": "asinh",
            "arcsin": "asin",
            "arctan": "atan",
            "arctanh": "atanh",
            "max": "fmax",
            "min": "fmin",
        }

        # TODO: Add more method mappings as needed
        # https://sbml.org/software/libsbml/5.18.0/docs/formatted/python-api/namespacelibsbml.html#a8e96a5a70569ae32655c6302638f6dc3

        # Replace method names with the C++ equivalents
        for old, new in method_map.items():
            formula = formula.replace(old, new)

        return formula

    def convert_function_body(self, fn_body: "ASTNode") -> str:
        """Convert a function body to its C++ equivalent.

        :param fn_body: The function body.
        :return: The C++ equivalent of the function body.
        """

        # Get the C++ conversion for each node function in the body
        node_list = self.get_sorted_nodes(fn_body)

        formula_mapping = {}
        for node in node_list:
            if node.isFunction():
                continue

            node_formula = formulaToString(node)
            node_formula_cpp = ""

            if "root" in node_formula:
                # If the function is root, we have to deal with it in an annoying way...

                # Split the root by the comma and rearrange
                split_formula = node_formula.split(",", 1)
                first_part = split_formula[0]
                second_part = split_formula[1]

                # Get the exponent of the root
                index = first_part.find("(")
                exponent = first_part[index + 1 : len(first_part)]

                # Get the base
                base = second_part[0 : len(second_part) - 1]

                node_formula_cpp = f"pow({base}, 1.0 / {exponent})"

            else:
                # The rest can be done by simple string replacement
                node_formula_cpp = node_formula
                node_formula_cpp = self.convert_formula(node_formula_cpp)

            formula_mapping[node_formula] = node_formula_cpp

        # Replace node formulas in function body with C++ equivalents
        formula = formulaToString(fn_body)
        for old, new in formula_mapping.items():
            formula = formula.replace(old, new)

        return formula

    def get_function_definition_arguments(
        self, fn_def: "FunctionDefinition"
    ) -> list[str]:
        """Get the list of arguments in a given function definition.

        :param fn_def: The function definition
        :return: List of arguments in the function definition
        """

        return [
            formulaToString(fn_def.getArgument(i))
            for i in range(fn_def.getNumArguments())
        ]

    def get_parameter_value(self, parameter: "Parameter") -> float:
        """Get initial parameter value.

        :return: The parameter value.
        """

        name = parameter.getName()
        if name and ("gamma" in name or "ComplexTransit" in name):
            return 1.0
        # Default for wnt and all other cases not specified above
        return 0.0

    def get_sorted_nodes(
        self, node: "ASTNode", node_list: list["ASTNode"] = []
    ) -> list["ASTNode"]:
        """Traverse an ASTNode tree and return an ordered list of nodes.

        :param node: The current ASTNode.
        :param node_list: A growing list of nodes in traversal order.
        :return: The node list with the current node and its sub-tree added.
        """

        left_node = node.getLeftChild()
        if left_node:
            self.get_sorted_nodes(left_node, node_list)

        node_list.append(node)

        right_node = node.getRightChild()
        if right_node:
            self.get_sorted_nodes(right_node, node_list)

        return node_list

    def get_species_concentration(self, species: "Species") -> float:
        """Get a initial species concentration.

        :return: The initial species concentration.
        """
        if species.isSetInitialAmount():
            return species.getInitialAmount()
        return species.getInitialConcentration()

    def get_ode_dict(self) -> dict[str, str]:
        """Get the ODEs as a dictionary of equations corresponding to each species.

        Each ODE will essentially be the sum of the products minus the sum of
        the reactants divided by the compartment volume

        :return: A dictionary of species and their corresponding ODEs
        """
        ode_dict = {}
        reactions = self.model.getListOfReactions()
        for reaction in reactions:
            reaction_id = reaction.getId()

            # Decompose reaction into sum of products minus sum of reactants
            products = reaction.getListOfProducts()
            for _, product in products.items():
                # Get the species concerning the product
                species_id = product.getSpecies()

                # TODO: Do we need to do something special with boundary conditions?
                # species = self.model.getSpecies(species_id)
                # if species.isSetBoundaryCondition() and not species.getBoundaryCondition():

                if species_id in ode_dict:
                    ode_dict[species_id] += " + " + reaction_id
                else:
                    ode_dict[species_id] = reaction_id

            reactants = reaction.getListOfReactants()
            for _, reactant in reactants.items():
                species_id = reactant.getSpecies()

                # TODO: Do we need to do something special with boundary conditions?
                # species = self.model.getSpecies(species_id)
                # if species.isSetBoundaryCondition() and not species.getBoundaryCondition():

                if species_id in ode_dict:
                    ode_dict[species_id] += " - " + reaction_id
                else:
                    ode_dict[species_id] = " - " + reaction_id

        return ode_dict

    def get_rules_dict(self) -> dict[str, str]:
        """Get a dictionary of species defined by reaction rules.

        :return: A dictionary of species and their corresponding rules
        """
        rules_dict = {
            r.getId(): self.convert_formula(r.getFormula())
            for r in self.model.getListOfRules()
        }
        return rules_dict

    def get_varname(self, obj: "SBase") -> str:
        """Get suitable C++ variable name for an object.

        :return: The variable name.
        """
        # TODO: Check that all custom variable names are unique

        name = varname_sanitize(obj.getName())
        if name:
            return name
        return obj.getId()

    def is_cc_model(self) -> bool:
        """Determine if the model is a Cell Cycle model.

        Models with events are assumed to be Cell Cycle models.
        Models with no events are assumed to be SRN models.

        :return: True if the model is a Cell Cycle model, False otherwise.
        """
        return not self.is_srn_model()

    def is_parameter_constant(self, i: int) -> bool:
        """Check if parameter is constant

        :param i: The index of the parameter
        :return: True if the parameter is constant, False otherwise
        """
        return self.model.getParameter(i).getConstant()

    def is_parameter_value_set(self, i: int) -> bool:
        """Check if parameter has a pre-set value

        :param i: The index of the parameter
        :return: True if the parameter has a pre-set value, False otherwise
        """
        return self.model.getParameter(i).isSetValue()

    def is_state_parameter(self, obj: "Species" | "Parameter") -> bool:
        """Check if a species or parameter is defined as a state parameter for Chaste.

        :param obj: The species or parameter to check.
        :return: True if defined as a parameter, False otherwise.
        """

        # Any species not defined by an ODE or rule is set as a state parameter
        if isinstance(obj, Species):
            species_id = obj.getId()
            return (species_id not in self.odes_dict) and (
                species_id not in self.rules_dict
            )

        if isinstance(obj, Parameter):
            # Also parameters with special strings in their name are state parameters.
            parameter_name = obj.getName()
            match_strings = ["wnt", "gamma", "ComplexTransit"]
            if parameter_name and any(s in parameter_name for s in match_strings):
                return True

            # Also parameters with unset values are state parameters.
            parameter_id = obj.getId()
            return (not obj.isSetValue()) and (parameter_id not in self.rules_dict)

        return False

    def is_srn_model(self) -> bool:
        """Determine if the model is a subcellular reaction network (SRN) model.

        Models with no events are assumed to be SRN models.
        Models with events are assumed to be Cell Cycle models.

        :return: True if the model is an SRN model, False otherwise.
        """
        return self.model.getNumEvents() == 0

    def get_timescale_multiplier(self) -> float:
        """Get the timescale multiplier.

        SBML uses seconds by default and Chaste uses hours.

        :return: The timescale multiplier.
        """

        unit_defs = self.model.getListOfUnitDefinitions()
        for unit_def in unit_defs:
            u_id = unit_def.getId()
            if u_id.lower() == "time":  # Do people ever call this something different?
                timescale = unit_def.getName().lower()
                if "minute" in timescale:
                    return 60.0
                elif "hour" in timescale:
                    return 1.0
        return 3600.0

    def write_code(self, output_directory=None):
        if output_directory:
            root_dir = pathlib.Path(output_directory)
        else:
            root_dir = pathlib.Path().cwd()

        hpp_file_path = root_dir / self.hpp_filename
        cpp_file_path = root_dir / self.cpp_filename

        with open(hpp_file_path, "w") as f:
            f.write(self.hpp_source)

        with open(cpp_file_path, "w") as f:
            f.write(self.cpp_source)
