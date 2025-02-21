import os
import shutil
import pathlib

from libsbml import SBMLReader

from . import filewriters


def get_node_index(node_list: type, node: type) -> int:
    """Get the corresponding node index from an ASTNodeList.

    :param node_list: The ASTNodeList
    :param node: The ASTNode
    :return: The index of the node in the list
    """
    # ASTNodeLists are not enumerable, nor are they Python lists, so we have to do this.
    for i in range(node_list.getSize()):
        if node_list.get(i) == node:
            return i
    return None


def get_traversal_order(node, node_list, order=[]) -> list[int]:
    """Traverse an ASTNode tree in-order and return the order of nodes.

    Assumes the tree is an ASTNode, as defined in SBML.

    :param node: The current ASTNode
    :param node_list: Unordered collection of all nodes for indexing
    :param order: A growing list of node indices in traversal order
    :return: The list of node indices in traversal order with this node's subtree added
    """

    l_node = node.getLeftChild()
    if l_node:
        get_traversal_order(l_node, node_list, order)

    i = get_node_index(node_list, node)
    order.append(i)

    r_node = node.getRightChild()
    if r_node:
        get_traversal_order(r_node, node_list, order)

    return order


def get_sub_function_equivalents(fn_body):
    """Check if any functions used don't match up with their C++ counterparts."""
    functs = {
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

    subfunct_equivalents = {}

    node_list = fn_body.getListOfNodes()
    node_order = get_traversal_order(fn_body, node_list)

    for i in node_order:
        node = node_list.get(i)
        if node.isFunction():
            node_formula = formula_to_string(node)  # TODO: Implement formula_to_string

            # If the function is root, we have to deal with it in an annoying way...
            if "root" in node_formula:
                # Split the root by the comma and rearrange
                split_formula = node_formula.split(",", 1)
                first_part = split_formula[0]
                second_part = split_formula[1]

                # Get the exponent of the root
                index = first_part.find("(")
                exponent = first_part[index + 1 : len(first_part)]

                # Get the base
                base = second_part[0 : len(second_part) - 1]
                subfunct_equivalents[node_formula] = (
                    "pow(" + base + ", " + "1.0 /" + exponent + ")"
                )

            # The rest can be done by simple string replacement
            else:
                new_formula = node_formula

                for key in functs:
                    new_formula = new_formula.replace(key, functs[key])

                subfunct_equivalents[node_formula] = new_formula

    return subfunct_equivalents


def convert_function_definition_body(self, fn_def: type) -> str:
    """Convrt the body of a given function definition to C++.

    :param fn_def: The function definition
    :return: The body of the function definition in C++.
    """
    fn_body = fn_def.getBody()
    subfunct_equivalents = GetSubfunctionEquivalents(fn_body)

    # Now replace all instances of those sub-functions with their C++ equivalents
    funct_body_formula = formulaToString(funct_def_body)

    for key in subfunct_equivalents:
        funct_body_formula = funct_body_formula.replace(key, subfunct_equivalents[key])

    return funct_body_formula


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

        self.hpp_file = f"{self.model_name}.hpp"
        self.cpp_file = f"{self.model_name}.cpp"

    def get_compartment_id(self, i: int) -> str:
        """Get the ID of the compartment

        :param i: The index of the compartment
        :return: The ID of the compartment
        """
        return self.model.getCompartment(i).getId()

    def get_compartment_name(self, i: int) -> str:
        """Get the name of the compartment

        :param i: The index of the compartment
        :return: The name of the compartment
        """
        return self.model.getCompartment(i).getName()

    def get_compartment_size(self, i: int) -> float:
        """Get the size of the compartment

        :param i: The index of the compartment
        :return: The size of the compartment
        """
        return self.model.getCompartment(i).getSize()

    def get_compartments(self) -> list[str]:
        """Get the list of compartment names

        :return: The list of compartments names
        """
        compartments = []
        n = self.model.getNumCompartments()
        for i in range(n):
            name = self.get_compartment_name(i)
            if not name:
                name = self.get_compartment_id(i)
            compartments.append(name)
        return compartments

    def get_function_definition_arguments(self, fn_def: type) -> list[str]:
        """Get the list of arguments in a given function definition.

        :param fn_def: The function definition
        :return: List of arguments in the function definition
        """
        args = []
        n = fn_def.getNumArguments()
        for i in range(n):
            arg = formula_to_string(
                fn_def.getArgument(i)
            )  # TODO: Implement formula_to_string
            args.append(arg)
        return args

    def get_function_definitions(self) -> dict[str, list[str]]:
        """Get a dictionary of function definitions in the model.

        :return: Dictionary of function IDs and corresponding arguments
        """
        fn_defs = {}
        n = self.model.getNumFunctionDefinitions()
        for i in range(n):
            fn_def = self.model.getFunctionDefinition(i)
            fn_id = fn_def.getId()
            fn_args = self.get_function_definition_arguments(fn_def)
            fn_defs[fn_id] = fn_args
        return fn_defs

    def get_events(self):
        n = self.model.getNumEvents()
        if n == 0:
            return []

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

    def get_parameter_id(self, i: int) -> str:
        """Get the ID of the parameter

        :param i: The index of the parameter
        :return: The ID of the parameter
        """
        return self.model.getParameter(i).getId()

    def get_parameter_name(self, i: int) -> str:
        """Get the name of the parameter

        :param i: The index of the parameter
        :return: The name of the parameter
        """
        return self.model.getParameter(i).getName()

    def get_parameter_value(self, i: int) -> float:
        """Get the value of the parameter

        :param i: The index of the parameter
        :return: The value of the parameter
        """
        return self.model.getParameter(i).getValue()

    def get_parameters(self) -> dict[str, float]:
        """Get a dictionary of parameter names and values.

        :return: A dictionary of parameter names and values.
        """
        parameters = {}
        n = self.model.getNumParameters()
        for i in range(n):
            name = self.get_parameter_name(i)
            value = 0.0  # Default for wnt and all cases not specified below
            if name:
                if "gamma" in name or "ComplexTransit" in name:
                    value = 1.0
            else:
                name = self.get_parameter_id(i)
            parameters[name] = value
        return parameters

    def get_species(self) -> dict[str, float]:
        """Get a dictionary of species IDs and initial concentrations.

        :return: A dictionary of species IDs and initial concentrations
        """
        species_dict = {}
        n = self.model.getNumSpecies()
        for i in range(n):
            species = self.model.getSpecies(i)
            if species.isSetInitialAmount():
                init_conc = species.getInitialAmount()
            else:
                init_conc = species.getInitialConcentration()
            species_dict[i] = init_conc
        return species

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

    def generate_code(self, output_directory=None):
        if self.is_srn_model():
            filewriters.WriteSrnModelToFile(self.model_name, self.model)
        else:
            filewriters.WriteCcmModelToFile(self.model_name, self.model)

        if output_directory:
            current_dir = pathlib.Path().cwd()
            out_dir = pathlib.Path(output_directory)

            shutil.move(current_dir / self.hpp_file, out_dir / self.hpp_file)
            shutil.move(current_dir / self.cpp_file, out_dir / self.cpp_file)
