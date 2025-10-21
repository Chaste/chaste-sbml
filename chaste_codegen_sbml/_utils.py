"""Utility functions for code generation."""

import re
from typing import TYPE_CHECKING

from libsbml import AST_NAME, formulaToString

if TYPE_CHECKING:
    from libsbml import ASTNode, FunctionDefinition, ListOf, SBase, Species


def convert_ast_formula(ast_formula: "ASTNode") -> str:
    """Convert SBML AST formula to equivalent C++ string.

    :param ast: The AST formula.
    :param convert_names: Whether to convert state variable and parameter names.
    :return: The equivalent C++ string.
    """
    return convert_str_formula(formulaToString(ast_formula))


def convert_str_formula(formula: str) -> str:
    """Convert SBML string formula to equivalent C++ string.

    :param ast: The string formula.
    :return: The equivalent C++ string.
    """
    # Convert all integer literals to doubles
    # TODO: Instead of regex, traverse AST and convert AST_INTEGER nodes to AST_REAL
    # This shouldn't apply to numbers encoded as <cn type="integer">.
    # Fix integer division directly in the AST instead.
    formula = re.sub(r"(?<!\.)(?<!e-|E-)\b[0-9]+\b(?!\.)", lambda x: f"{x[0]}.0", formula)

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

    # TODO: From SBML Level 3 upwards, log defaults to base 10.
    # SBML versions lower than 3 default to base e.
    # See https://sbml.org/software/libsbml/5.18.0/docs/formatted/python-api/namespacelibsbml.html#a8e96a5a70569ae32655c6302638f6dc3  # noqa: B950

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


def sort_formulas(formulas: list[tuple[str, "ASTNode"]]) -> list[int]:
    """Sort formulas based on their dependency.

    Formulas are sorted such that if formula A depends on B (i.e. A -> B),
    then B comes before A. It is assumed that the input formulas are acyclic.
    This function can't sort cyclic dependencies such as A -> B -> C -> A.

    :param formulas: A dictionary of (lhs, rhs) tuples.
    :return: A list of sorted formula indices.
    """

    def _search_formula(node: "ASTNode", name: str) -> bool:
        """Recursively search for a variable name in an AST formula.

        :param node: The AST node to search.
        :param name: The variable name to search for.

        :return: True if the variable name is found, False otherwise.
        """
        if node is None:
            return False
        if node.getType() == AST_NAME and node.getName() == name:
            return True

        for i in range(node.getNumChildren()):
            child = node.getChild(i)
            if _search_formula(child, name):
                return True
        return False

    _compare_cache = dict()

    def _compare_formulas(index_a: int, index_b: int) -> int:
        """Compare two formulas based on their dependency.

        :param index_a: The index of the first formula (A).
        :param index_b: The index of the second formula (B).

        :return: An integer indicating the order of the formulas.
            -1 if A < B (A comes before B)
            1 if A > B (A comes after B)
            0 if the order doesn't matter
        """
        order = _compare_cache.get((index_a, index_b), None)
        if order is not None:
            return order

        var_a = formulas[index_a][0]
        var_b = formulas[index_b][0]

        rhs_a = formulas[index_a][1]
        rhs_b = formulas[index_b][1]

        # Check if var_a is in rhs_b
        if _search_formula(rhs_b, var_a):
            # var_a is used in rhs_b: rule_a comes before rule_b
            _compare_cache[(index_a, index_b)] = -1
            _compare_cache[(index_b, index_a)] = 1
            return -1

        # Check if var_b is in rhs_a
        if _search_formula(rhs_a, var_b):
            # var_b is used in rhs_a: rule_b comes before rule_a
            _compare_cache[(index_a, index_b)] = 1
            _compare_cache[(index_b, index_a)] = -1
            return 1

        # Order doesn't matter
        _compare_cache[(index_a, index_b)] = 0
        _compare_cache[(index_b, index_a)] = 0
        return 0

    def _insertion_sort() -> list[int]:
        """Sort formulas using insertion sort.

        :return: A list of sorted formula indices.
        """
        # We need to compare each formula to all the others until we find a
        # non-zero comparison i.e. a +1 or -1 match (or until we exhaust
        # all options) because of cases such as:
        # Initial order: (a, b, c)
        # a == b (order doesn't matter, comparison is 0);
        # b == c (order doesn't matter, comparison is 0);
        # a > c (a should come after c, comparison is 1);
        # If we only compare (a, b) and (b, c), no changes will be made.
        sorted_indices = []

        for index_a in range(len(formulas)):
            for i, index_b in enumerate(sorted_indices):
                if _compare_formulas(index_a, index_b) < 0:  # rule_a < rule_b
                    sorted_indices.insert(i, index_a)
                    break
            else:
                # rule_a comes after everything already in sorted_rules
                sorted_indices.append(index_a)

        return sorted_indices

    return _insertion_sort()


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
