"""Utility functions for code generation."""

import re
from typing import TYPE_CHECKING

from libsbml import AST_FUNCTION_DELAY, Compartment, formulaToString

from chaste_codegen_sbml._config import DERIVATIVE_PREFIX, DERIVATIVE_SUFFIX

if TYPE_CHECKING:
    from libsbml import ASTNode, FunctionDefinition


def convert_ast_formula(ast_formula: "ASTNode") -> str:
    """Convert SBML AST formula to equivalent C++ string.

    :param ast: The AST formula.
    :param convert_names: Whether to convert state variable and parameter names.
    :return: The equivalent C++ string.
    """
    unsupported_functions = ["delay"]
    for func in unsupported_functions:
        if search_ast_type(ast_formula, AST_FUNCTION_DELAY):
            raise NotImplementedError(f"SBML function not supported: '{func}'.")

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
        "time": "time",
        "t": "time",
        "s": "time",
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

    results = re.findall(r"rateOf\(([^)]+)\)", cpp_formula)
    if results:
        for var in results:
            cpp_formula = cpp_formula.replace(f"rateOf({var})", get_derivative_var(var))
    return cpp_formula


def get_compartment_size(compartment: "Compartment") -> float:
    """Get a compartment size.

    :return: The compartment size.
    """
    if compartment.isSetSize():
        return compartment.getSize()
    return 1.0


def get_derivative_var(var: str) -> str:
    """Get the derivative variable name for a given variable.

    :param var: The variable name.
    :return: The derivative variable name.
    """
    return f"{DERIVATIVE_PREFIX}{var}{DERIVATIVE_SUFFIX}"


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
