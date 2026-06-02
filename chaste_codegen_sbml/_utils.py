"""Utility functions for code generation."""

from typing import TYPE_CHECKING

from libsbml import Compartment, formulaToString

if TYPE_CHECKING:
    from libsbml import ASTNode, FunctionDefinition


def get_compartment_size(compartment: "Compartment") -> float:
    """Get a compartment size.

    :return: The compartment size.
    """
    if compartment.isSetSize():
        return compartment.getSize()
    return 1.0


def generate_header_guard(filename: str) -> str:
    """Generate a C++ header guard from a filename.

    :param filename: The filename.
    :return: The header guard.
    """
    name = to_cpp_name(filename)
    if not name:
        return ""

    prev = name[0]
    guard = [prev.upper()]

    for char in name[1:]:
        if char.isupper():
            # Add _ before uppercase chars if not sequence of uppercase chars
            if not prev.isupper():
                guard.append("_")
            guard.append(char)
        elif char.islower():
            guard.append(char.upper())
        elif char.isdigit():
            # Add _ between lowercase char and digit
            if prev.islower():
                guard.append("_")
            guard.append(char)
        elif guard[-1] != "_":
            # Replace non-alphanumeric chars with "_"; merging successive "_"s
            guard.append("_")
        prev = char

    return "".join(guard) + "_"


def get_function_definition_arguments(fn_def: "FunctionDefinition") -> list[str]:
    """Get the list of arguments in a given function definition.

    :param fn_def: The function definition
    :return: List of arguments in the function definition
    """
    n = fn_def.getNumArguments()
    return [formulaToString(fn_def.getArgument(i)) for i in range(n)]


def search_ast_type(root: "ASTNode", node_type: int) -> bool:
    """Recursively search the AST for a node of a certain type.

    :param root: The root node of the AST.
    :param node_type: The type of node to search for.

    :return: True if a node matching the spec is found, False otherwise.
    """
    if root is None:
        return False

    if root.getType() == node_type:
        return True

    for i in range(root.getNumChildren()):
        child = root.getChild(i)
        if search_ast_type(child, node_type):
            return True

    return False


def to_camel_case(name: str) -> str:
    """Convert an input name to an alphanumeric name in camel case.

    :param name: The variable name.
    :return: The variable name in camel case.
    """
    camel = []

    caps = False  # True: capitalize next letter
    for char in name.strip():
        if char.isalpha():
            if caps:
                camel.append(char.upper())
                caps = False
            else:
                camel.append(char)
        elif char.isdigit():
            camel.append(char)
            caps = True
        else:
            # Skip other chars
            caps = True

    return "".join(camel)


def to_cpp_name(name: str) -> str:
    """Sanitize an input name to a C++ compatible alphanumeric name.

    :param name: The variable name.
    :return: The sanitized variable name.
    """
    name_ = name.strip()
    if not name_:
        return ""

    cpp_name = []
    # Prefix with "_" if name doesn't start with a letter or "_"
    if name_[0] != "_" and not name_[0].isalpha():
        cpp_name.append("_")

    for char in name_:
        if char.isalpha() or char.isdigit() or char == "_":
            cpp_name.append(char)
        else:
            # Replace other chars with "_"
            cpp_name.append("_")

    # TODO: Check for C++ keywords

    return "".join(cpp_name)
