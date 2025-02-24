"""Utility functions for code generation."""

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from libsbml import ListOf, SBase


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


def varname_camel(name: str, keep_underscores: bool = False) -> str:
    """Convert an input string to a C++ compatible alphanumeric string in camel case.

    :param name: The variable name.
    :param keep_underscores: Whether to keep underscores in the variable name.
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
        elif keep_underscores and char == "_":
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
    for char in name:
        if char.isalpha() or char.isdigit():
            var_name.append(char)
        else:
            var_name.append("_")
    return "".join(var_name)
