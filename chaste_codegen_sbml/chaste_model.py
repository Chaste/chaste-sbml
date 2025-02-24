import os
import pathlib
from typing import TYPE_CHECKING

from libsbml import SBMLReader

if TYPE_CHECKING:
    from libsbml import ASTNode, Compartment, FunctionDefinition, Parameter, Species


ODE_SUFFIX = "OdeSystem"
TAB = "    "


class ChasteModel:
    """Holds information about the SBML model for which code is to be generated."""

    def __init__(self, sbml_file: str, model_name: str = None) -> None:
        self.sbml_file = sbml_file
        self.sbml_model = SBMLReader().readSBMLFromFile(sbml_file).getModel()

        if model_name:
            self.model_name = model_name
        else:
            self.model_name = os.path.splitext(os.path.basename(sbml_file))[0]

        self.ode_system_name = f"{self.model_name}OdeSystem"

        self.hpp_filename = f"{self.model_name}.hpp"
        self.cpp_filename = f"{self.model_name}.cpp"

        self.hpp_source = ""
        self.cpp_source = ""

    def convert_function_body(self, fn_body: "ASTNode") -> str:
        """Convert a function body to its C++ equivalent.

        :param fn_body: The function body.
        :return: The C++ equivalent of the function body.
        """

        op_mapping = {
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

        # Get the C++ conversion for each node function in the body
        node_list = self.get_sorted_nodes(fn_body)

        formula_mapping = {}
        for node in node_list:
            if node.isFunction():
                continue

            node_formula = self.formula_to_string(node)
            cpp_formula = ""

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

                cpp_formula = f"pow({base}, 1.0 / {exponent})"

            else:
                # The rest can be done by simple string replacement
                cpp_formula = node_formula

                for key, val in op_mapping.items():
                    cpp_formula = cpp_formula.replace(key, val)

            formula_mapping[node_formula] = cpp_formula

        # Replace node formulas in function body with C++ equivalents
        formula = self.formula_to_string(fn_body)
        for key, val in formula_mapping.items():
            formula = formula.replace(key, val)

        return formula

    # TODO: Implement formula_to_string
    def formula_to_string(self, node: ASTNode) -> str:
        """Convert an ASTNode formula to a string.

        :param node: The ASTNode
        :return: The formula as a string
        """
        pass

    def get_function_definition_arguments(
        self, fn_def: "FunctionDefinition"
    ) -> list[str]:
        """Get the list of arguments in a given function definition.

        :param fn_def: The function definition
        :return: List of arguments in the function definition
        """

        return [
            self.formula_to_string(fn_def.getArgument(i))
            for i in range(fn_def.getNumArguments())
        ]

    def get_identifier(self, obj: "Compartment" | "Parameter") -> str:
        """Get suitable C++ variable name for a compartment or parameter.

        :return: The object's C++ identifier.
        """

        name = obj.getName()
        if name:
            return name
        return obj.getId()

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

    def get_odes(self) -> dict[str, str]:
        """Get the ODEs as a dictionary of equations corresponding to each species.

        Each ODE will essentially be the sum of the products minus the sum of
        the reactants divided by the compartment volume

        :return: A dictionary of species and their corresponding ODEs
        """
        ode_dict = {}
        n = self.model.getNumReactions()
        for i in range(n):
            reaction = self.model.getReaction(i)
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

    def get_rules(self) -> dict[str, str]:
        """Get a dictionary of species defined by reaction rules.

        :return: A dictionary of species and their corresponding rules
        """
        rules_dict = {}
        n = self.model.getNumRules()
        for i in range(n):
            rule = self.model.getRule(i)
            rule_id = rule.getId()
            rule_formula = rule.getFormula()
            rules_dict[rule_id] = rule_formula
        return rules_dict

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

    def is_srn_model(self) -> bool:
        """Determine if the model is a subcellular reaction network (SRN) model.

        Models with no events are assumed to be SRN models.
        Models with events are assumed to be Cell Cycle models.

        :return: True if the model is an SRN model, False otherwise.
        """
        return self.model.getNumEvents() == 0

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
