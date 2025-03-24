"""Utility functions for code generation."""

import re
from typing import TYPE_CHECKING

from libsbml import formulaToString

if TYPE_CHECKING:
    from libsbml import ASTNode, FunctionDefinition, ListOf, SBase, Species


def convert_formula(formula: str) -> str:
    """Convert a formula string to its C++ equivalent.

    :param formula: The formula string.
    :return: The C++ equivalent of the formula.
    """

    if not formula:
        return ""

    # Sorted by length to avoid partial replacements
    token_map = {
        "arccsch": "acsch",
        "arcsech": "asech",
        "arcsinh": "asinh",
        "arctanh": "atanh",
        "arccos": "acos",
        "arccsc": "acsc",
        "arcsec": "asec",
        "arcsin": "asin",
        "arctan": "atan",
        "abs": "fabs",
        "max": "fmax",
        "min": "fmin",
        "and": "&&",
        "geq": ">=",
        "leq": "<=",
        "neq": "!=",
        "not": "!",
        "eq": "==",
        "gt": ">",
        "lt": "<",
        "or": "||",
    }

    # TODO: Add more token mappings as needed
    # https://sbml.org/software/libsbml/5.18.0/docs/formatted/python-api/namespacelibsbml.html#a8e96a5a70569ae32655c6302638f6dc3

    # Replace token names with the C++ equivalents
    for sbml_token, cpp_token in token_map.items():
        formula = re.sub(rf"\b{sbml_token}\b", cpp_token, formula)

    return formula


def convert_function_body(fn_body: "ASTNode") -> str:
    """Convert a function body to its C++ equivalent.

    :param fn_body: The function body.
    :return: The C++ equivalent of the function body.
    """
    # Get the C++ conversion for each node function in the body
    node_list = sort_nodes(fn_body)

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
            node_formula_cpp = convert_formula(node_formula_cpp)

        formula_mapping[node_formula] = node_formula_cpp

    # Replace node formulas in function body with C++ equivalents
    formula = formulaToString(fn_body)
    for sbml_formula, cpp_formula in formula_mapping.items():
        formula = re.sub(rf"\b{sbml_formula}\b", cpp_formula, formula)

    return formula


def formula_to_cpp(formula: "ASTNode") -> str:
    """Convert an ASTNode formula to an equivalent C++ string.

    :param formula: The ASTNode.
    :return: The equivalent C++ formula string.
    """
    pass


def get_function_definition_arguments(fn_def: "FunctionDefinition") -> list[str]:
    """Get the list of arguments in a given function definition.

    :param fn_def: The function definition
    :return: List of arguments in the function definition
    """
    n = fn_def.getNumArguments()
    return [formulaToString(fn_def.getArgument(i)) for i in range(n)]


def get_index_by_obj(obj: "SBase", listof: "ListOf") -> int:
    """Return the index of an object in a libsbml.ListOf.

    :param o: The object.
    :return: The index of the object in the ListOf.
    """
    for i, o in enumerate(listof):
        if obj == o:
            return i
    return None


def get_index_by_id(obj_id: str, listof: "ListOf") -> int:
    """Return the index of an object in a libsbml.ListOf by its id.

    :param o: The object.
    :return: The index of the object in the ListOf.
    """
    for i, o in enumerate(listof):
        if obj_id == o.getId():
            return i
    return None


def get_species_concentration(species: "Species") -> float:
    """Get a initial species concentration.

    :return: The initial species concentration.
    """
    if species.isSetInitialAmount():
        return species.getInitialAmount()
    return species.getInitialConcentration()


def sort_nodes(node: "ASTNode", node_list: list["ASTNode"] = None) -> list["ASTNode"]:
    """Traverse an ASTNode tree and return an ordered list of nodes.

    :param node: The current ASTNode.
    :param node_list: A growing list of nodes in traversal order.
    :return: The list of nodes in traversal order.
    """
    if node_list is None:
        node_list = []

    left_node = node.getLeftChild()
    if left_node:
        sort_nodes(left_node, node_list)

    node_list.append(node)

    right_node = node.getRightChild()
    if right_node:
        sort_nodes(right_node, node_list)

    return node_list


def varname_camelcase(name: str) -> str:
    """Convert an input string to a C++ compatible alphanumeric string in camel case.

    :param name: The variable name.
    :return: The variable name in camel case.
    """
    camel_name = []

    next_caps = False
    for char in name:
        if char.isalpha():
            if next_caps:
                camel_name.append(char.upper())
                next_caps = False
            else:
                camel_name.append(char)
        elif char.isdigit():
            camel_name.append(char)
            next_caps = True
        else:
            next_caps = True

    return "".join(camel_name)


def varname_sanitize(name: str) -> str:
    """Convert an input string to a C++ compatible alphanumeric string.

    :param name: The variable name.
    :return: The variable name in C++ alphanumeric.
    """
    var_name = []

    name = name.strip()

    # Prefix with "_" if name starts with a number
    if name and name[0].isdigit():
        var_name.append("_")

    skip_underscores = False
    for char in name:
        if char.isalpha() or char.isdigit() or char == "_":
            var_name.append(char)
            skip_underscores = False
        else:
            # Replace non-alphanumeric chars with "_"; merging successive "_"s
            if not skip_underscores:
                var_name.append("_")
                skip_underscores = True

    return "".join(var_name)
