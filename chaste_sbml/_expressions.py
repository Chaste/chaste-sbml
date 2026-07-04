"""Translation of SBML MathML/AST expressions into C++ formula strings.

Extracted from ChasteSbmlModel so the expression-formatting logic lives on its own. The
functions are pure apart from the libsbml AST objects they operate on: ``formula_to_string`` is
given the model's variable-type map and state-variable list as arguments rather than reaching
into the model.
"""

import re
from typing import TYPE_CHECKING, Optional

from libsbml import AST_FUNCTION_DELAY, AST_NAME, AST_NAME_AVOGADRO, formulaToL3String

from ._config import CHASTE_PREFIX, PREFIX_SEP, VarType
from ._utils import search_ast_type

if TYPE_CHECKING:
    from libsbml import ASTNode


def collect_ast_names(node: "ASTNode", names: set) -> None:
    """Recursively collect the identifiers referenced by name in an AST.

    :param node: The root AST node.
    :param names: A set to accumulate referenced identifiers into.
    """
    if node is None:
        return
    if node.isName():
        names.add(node.getName())
    for i in range(node.getNumChildren()):
        collect_ast_names(node.getChild(i), names)


def substitute_ast_names(node: "ASTNode", replacements: dict) -> "ASTNode":
    """Return a copy of an AST with each named node replaced by a replacement AST.

    :param node: The AST to copy and substitute into.
    :param replacements: A mapping of identifier to replacement ASTNode.
    :return: A new ASTNode with the substitutions applied.
    """
    if node.isName() and node.getName() in replacements:
        return replacements[node.getName()].deepCopy()
    result = node.deepCopy()
    stack = [result]
    while stack:
        current = stack.pop()
        for i in range(current.getNumChildren()):
            child = current.getChild(i)
            if child.isName() and child.getName() in replacements:
                current.replaceChild(i, replacements[child.getName()].deepCopy())
            else:
                stack.append(child)
    return result


def replace_avogadro_csymbol(node: "ASTNode", placeholder: str) -> None:
    """Recursively rename avogadro csymbol nodes to a placeholder identifier.

    The avogadro csymbol and a parameter both named 'avogadro' are distinct AST nodes that
    ``formulaToL3String`` renders identically. Renaming the csymbol to a placeholder keeps it
    distinct so it can be mapped to ``sm::AVOGADRO`` while the parameter keeps its own name.

    :param node: The root AST node.
    :param placeholder: The identifier to rename avogadro csymbol nodes to.
    """
    if node.getType() == AST_NAME_AVOGADRO:
        node.setType(AST_NAME)
        node.setName(placeholder)
    for i in range(node.getNumChildren()):
        replace_avogadro_csymbol(node.getChild(i), placeholder)


def strip_ast_units(node: "ASTNode") -> None:
    """Recursively strip units annotations from AST nodes.

    SBML allows numeric literals to carry units annotations (e.g. ``<cn sbml:units="mole">``).
    ``formulaToL3String`` includes these in its output (e.g. ``0.00015 mole``), which is not
    valid C++. Stripping them here is safe because the units carry no mathematical information.

    :param node: The root AST node.
    """
    if node.isSetUnits():
        node.unsetUnits()
    for i in range(node.getNumChildren()):
        strip_ast_units(node.getChild(i))


def convert_infix_operator_to_function_syntax(formula: str, operator: str, function_name: str) -> str:
    """Convert infix operator expressions to function syntax.

    Example: with operator='^' and function_name='pow', rewrites
    ``a ^ b`` to ``pow(a, b)``.

    This parser handles parenthesized operands (including nested parentheses)
    and simple symbolic/numeric tokens.

    :param formula: Formula text in SBML infix style.
    :param operator: Infix operator token to rewrite.
    :param function_name: Function name used for the replacement.
    :return: Formula text with infix operators rewritten as function calls.
    """
    if not operator:
        raise ValueError("operator must be a non-empty string")
    if not function_name:
        raise ValueError("function_name must be a non-empty string")

    token_chars = set("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_:.eE")

    def read_left_operand(s: str, operator_index: int) -> tuple[int, int] | None:
        i = operator_index - 1
        while i >= 0 and s[i].isspace():
            i -= 1
        if i < 0:
            return None

        if s[i] == ")":
            depth = 1
            j = i - 1
            while j >= 0:
                if s[j] == ")":
                    depth += 1
                elif s[j] == "(":
                    depth -= 1
                    if depth == 0:
                        break
                j -= 1
            if depth != 0:
                return None
            # j points to the matching '('; also include any function name before it
            k = j - 1
            while k >= 0 and s[k] in token_chars:
                k -= 1
            return (k + 1, i + 1)

        j = i
        while j >= 0 and s[j] in token_chars:
            j -= 1

        start = j + 1
        end = i + 1
        if start >= end:
            return None
        return (start, end)

    def read_right_operand(s: str, operator_index: int) -> tuple[int, int] | None:
        i = operator_index + len(operator)
        n = len(s)
        while i < n and s[i].isspace():
            i += 1
        if i >= n:
            return None

        if s[i] == "(":
            depth = 1
            j = i + 1
            while j < n:
                if s[j] == "(":
                    depth += 1
                elif s[j] == ")":
                    depth -= 1
                    if depth == 0:
                        return (i, j + 1)
                j += 1
            return None

        j = i
        if s[j] in "+-":
            j += 1

        while j < n and s[j] in token_chars:
            j += 1

        if j <= i:
            return None
        return (i, j)

    max_rewrites = 500
    rewrites = 0
    search_start = 0
    while rewrites < max_rewrites:
        operator_index = formula.find(operator, search_start)
        if operator_index < 0:
            break

        left = read_left_operand(formula, operator_index)
        right = read_right_operand(formula, operator_index)
        if left is None or right is None:
            search_start = operator_index + len(operator)
            continue

        l_start, l_end = left
        r_start, r_end = right

        lhs = formula[l_start:l_end].strip()
        rhs = formula[r_start:r_end].strip()
        replacement = f"{function_name}({lhs}, {rhs})"

        formula = formula[:l_start] + replacement + formula[r_end:]
        search_start = l_start + len(replacement)
        rewrites += 1

    return formula


def convert_infix_power_to_function_syntax(formula: str) -> str:
    """Convert infix power expressions (a ^ b) to function syntax pow(a, b)."""
    return convert_infix_operator_to_function_syntax(formula=formula, operator="^", function_name="pow")


def formula_to_string(
    math: "ASTNode",
    variable_types: dict,
    state_variables: list,
    local_parameters: Optional[list[dict[str, str]]] = None,
) -> str:
    """Convert an AST math formula to an equivalent C++ string.

    :param math: The AST math formula.
    :param variable_types: Mapping of model id to :class:`VarType`, used to tell model variables
        apart from SBML constants/functions of the same spelling.
    :param state_variables: The model's state-variable dicts, used to resolve ``rateOf`` to a
        state variable's derivative id.
    :param local_parameters: Local parameters in scope (e.g. a reaction's kinetic-law
        parameters). These shadow global symbols of the same name and are constant, so
        ``rateOf`` applied to one is zero.
    :return: The equivalent C++ string.
    """
    unsupported_functions = ["delay"]
    for func in unsupported_functions:
        if search_ast_type(math, AST_FUNCTION_DELAY):
            raise NotImplementedError(f"SBML function not supported: '{func}'.")

    strip_ast_units(math)
    # The avogadro csymbol and a parameter both named 'avogadro' are distinct AST nodes that
    # render identically. Rename the csymbol to a placeholder (mapped to sm::AVOGADRO in the
    # constants below) so it stays distinct from a same-named parameter.
    avogadro_placeholder = f"{CHASTE_PREFIX}{PREFIX_SEP}avogadro"
    replace_avogadro_csymbol(math, avogadro_placeholder)
    formula = formulaToL3String(math)

    # Convert all integer literals to doubles to fix integer division.
    # TODO: Perhaps instead of regex, traverse AST and convert some AST_INTEGER
    # nodes to AST_REAL. This should only need to apply to numbers used in a division.
    formula = re.sub(r"(?<!\.)(?<!e-|E-)\b[0-9]+\b(?!\.)", lambda x: f"{x[0]}.0", formula)

    formula = convert_infix_power_to_function_syntax(formula)

    # SBML contants to be replaced with C++ equivalents
    constants = {
        "avogadro": "sm::AVOGADRO",
        avogadro_placeholder: "sm::AVOGADRO",
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
        "plus": "plus",
        "quotient": "quotient",
        "root": "root",
        "sec": "sec",
        "sech": "sech",
        "sqr": "sqr",
        "times": "times",
        "xor": "xor_",
    }

    # TODO: From SBML Level 3 upwards, log defaults to base 10.
    # SBML versions lower than 3 default to base e.
    # See https://sbml.org/software/libsbml/5.18.0/docs/formatted/python-api/namespacelibsbml.html#a8e96a5a70569ae32655c6302638f6dc3  # noqa: B950

    tokens = re.findall(r"\w+|\W+", formula)

    local_param_ids = {param["id"] for param in (local_parameters or [])}

    cpp_tokens = []
    for token in tokens:
        cpp_token = token

        # Replace function names and constants, but only when the token is not an actual
        # model variable or a local parameter (e.g. a species named "s" or "t" must not be
        # replaced with the SBML time symbol, and a local parameter named "avogadro" must
        # stay the local, distinct from the avogadro csymbol handled above).
        if token in constants and token not in variable_types and token not in local_param_ids:
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
            rate = "0.0"
            # A local parameter shadows any global symbol of the same name and is
            # constant, so its rate of change is zero. Only fall back to the global
            # variable when the name is not a local parameter.
            if var not in local_param_ids:
                if variable_types.get(var, VarType.UNKNOWN) == VarType.STATE_VARIABLE:
                    for state_var in state_variables:
                        if state_var["id"] == var:
                            rate = state_var["derivative_id"]
                            break
            cpp_formula = cpp_formula.replace(f"rateOf({var})", rate)
    return cpp_formula
