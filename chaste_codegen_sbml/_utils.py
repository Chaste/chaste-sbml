"""Utility functions for code generation."""

from typing import TYPE_CHECKING

from libsbml import formulaToString

if TYPE_CHECKING:
    from libsbml import ASTNode, FunctionDefinition, ListOf, SBase, Species


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
