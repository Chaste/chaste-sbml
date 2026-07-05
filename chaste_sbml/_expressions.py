"""Translation of SBML MathML/AST expressions into C++ formula strings.

Extracted from ChasteSbmlModel so the expression-formatting logic lives on its own. The
functions are pure apart from the libsbml AST objects they operate on: ``formula_to_string`` is
given the model's variable-type map and state-variable list as arguments rather than reaching
into the model.
"""

import re
from typing import TYPE_CHECKING, Optional

from libsbml import (
    AST_FUNCTION,
    AST_FUNCTION_DELAY,
    AST_FUNCTION_POWER,
    AST_NAME,
    AST_NAME_AVOGADRO,
    AST_POWER,
    AST_RELATIONAL_EQ,
    AST_RELATIONAL_GEQ,
    AST_RELATIONAL_GT,
    AST_RELATIONAL_LEQ,
    AST_RELATIONAL_LT,
    formulaToL3String,
)

from ._config import CHASTE_PREFIX, PREFIX_SEP, VarType

if TYPE_CHECKING:
    from libsbml import ASTNode

    from ._records import LocalParameter, StateVariable


def collect_ast_names(node: Optional["ASTNode"], names: set) -> None:
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


# SBML MathML relational operators, and the SbmlMath (``sm::``) helper each maps to.
_NARY_RELATIONAL_NAMES = {
    AST_RELATIONAL_LT: "lt",
    AST_RELATIONAL_GT: "gt",
    AST_RELATIONAL_LEQ: "leq",
    AST_RELATIONAL_GEQ: "geq",
    AST_RELATIONAL_EQ: "eq",
}


def rewrite_nary_relational(node: Optional["ASTNode"]) -> None:
    """Rewrite n-ary (>2 operand) relational AST nodes into function-call form, in place.

    SBML MathML relationals are n-ary -- ``lt(a, b, c)`` means ``a < b < c`` (i.e. ``a < b`` and
    ``b < c``), and ``eq(a, b, c)`` means all equal -- but ``formulaToL3String`` renders them as a
    chained infix expression (``a < b < c``), which C++ mis-evaluates left-to-right as
    ``(a < b) < c``. Converting each such node to a function call (rendered ``lt(a, b, c)``, later
    mapped to the variadic ``sm::lt`` helper that implements the correct semantics) fixes this.
    Binary relationals are left as infix operators, so existing generated code is unchanged.
    ``neq`` is binary-only in MathML, so it is never n-ary and is not rewritten here.

    :param node: The root AST node to rewrite.
    """
    if node is None:
        return
    name = _NARY_RELATIONAL_NAMES.get(node.getType())
    if name is not None and node.getNumChildren() > 2:
        node.setType(AST_FUNCTION)
        node.setName(name)
    for i in range(node.getNumChildren()):
        rewrite_nary_relational(node.getChild(i))


def search_ast_type(root: Optional["ASTNode"], node_type: int) -> bool:
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


def rewrite_power(node: Optional["ASTNode"]) -> None:
    """Rewrite power-operator (``a ^ b``) AST nodes into ``pow()`` function-call form, in place.

    ``formulaToL3String`` renders a power as the infix ``a ^ b``, which is not valid C++ (``^`` is
    bitwise xor there). Converting each power node to a ``pow`` call before stringifying (rendered
    ``pow(a, b)``, later mapped to ``std::pow`` by the function-name token map) handles nesting
    correctly -- a string-level rewrite of ``a ^ b ^ c`` would convert only the outer operator and
    leave an inner ``^`` behind. Both power node types are handled: ``AST_POWER`` (from ``^``/``pow``
    in L3 infix, and from MathML ``<power>`` that libsbml normalises) and ``AST_FUNCTION_POWER``
    (from other MathML ``<power>`` forms) -- both render as ``a ^ b``.

    :param node: The root AST node to rewrite.
    """
    if node is None:
        return
    if node.getType() in (AST_POWER, AST_FUNCTION_POWER):
        node.setType(AST_FUNCTION)
        node.setName("pow")
    for i in range(node.getNumChildren()):
        rewrite_power(node.getChild(i))


# Placeholder id the avogadro csymbol is renamed to (see replace_avogadro_csymbol) so it stays
# distinct from a same-named parameter; mapped to sm::AVOGADRO in _CONSTANTS below.
AVOGADRO_PLACEHOLDER = f"{CHASTE_PREFIX}{PREFIX_SEP}avogadro"

# SBML symbolic constants replaced with their C++ equivalents. (true/false are deliberately absent.)
_CONSTANTS = {
    "avogadro": "sm::AVOGADRO",
    AVOGADRO_PLACEHOLDER: "sm::AVOGADRO",
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

# SBML functions whose name matches the C++ (std::) equivalent.
_UNCHANGED_FUNCTIONS = {
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

# SBML functions with a different name in C++ (std::).
_RENAMED_FUNCTIONS = {
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

# SBML functions with custom implementations in SbmlMath (sm::).
_CUSTOM_FUNCTIONS = {
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


def formula_to_string(
    math: "ASTNode",
    variable_types: dict,
    state_variables: list["StateVariable"],
    local_parameters: Optional[list["LocalParameter"]] = None,
) -> str:
    """Convert an AST math formula to an equivalent C++ string.

    :param math: The AST math formula. It is not mutated (a copy is normalised internally).
    :param variable_types: Mapping of model id to :class:`VarType`, used to tell model variables
        apart from SBML constants/functions of the same spelling.
    :param state_variables: The model's state variables, used to resolve ``rateOf`` to a state
        variable's derivative id.
    :param local_parameters: Local parameters in scope (e.g. a reaction's kinetic-law
        parameters). These shadow global symbols of the same name and are constant, so
        ``rateOf`` applied to one is zero.
    :return: The equivalent C++ string.
    """
    if search_ast_type(math, AST_FUNCTION_DELAY):
        raise NotImplementedError("SBML function not supported: 'delay'.")

    # Normalise a copy so the caller's AST -- and the libsbml model it belongs to -- is not mutated.
    math = math.deepCopy()

    strip_ast_units(math)
    # The avogadro csymbol and a parameter both named 'avogadro' are distinct AST nodes that
    # render identically. Rename the csymbol to a placeholder (mapped to sm::AVOGADRO in _CONSTANTS)
    # so it stays distinct from a same-named parameter.
    replace_avogadro_csymbol(math, AVOGADRO_PLACEHOLDER)
    # Route operators that formulaToL3String would render as C++-invalid infix to function-call
    # form before stringifying: n-ary relationals (a < b < c) to the sm:: helpers, and the power
    # operator (a ^ b) to pow(). Binary relationals stay as infix operators.
    rewrite_nary_relational(math)
    rewrite_power(math)
    formula = formulaToL3String(math)

    # Convert all integer literals to doubles to fix integer division.
    # TODO: Perhaps instead of regex, traverse AST and convert some AST_INTEGER
    # nodes to AST_REAL. This should only need to apply to numbers used in a division.
    formula = re.sub(r"(?<!\.)(?<!e-|E-)\b[0-9]+\b(?!\.)", lambda x: f"{x[0]}.0", formula)

    tokens = re.findall(r"\w+|\W+", formula)

    local_param_ids = {param.id for param in (local_parameters or [])}

    cpp_tokens = []
    for token in tokens:
        cpp_token = token

        # Replace function names and constants, but only when the token is not an actual
        # model variable or a local parameter (e.g. a species named "s" or "t" must not be
        # replaced with the SBML time symbol, and a local parameter named "avogadro" must
        # stay the local, distinct from the avogadro csymbol handled above).
        if token in _CONSTANTS and token not in variable_types and token not in local_param_ids:
            cpp_token = f"{_CONSTANTS[token]}"

        elif token in _UNCHANGED_FUNCTIONS:
            cpp_token = f"std::{token}"

        elif token in _RENAMED_FUNCTIONS:
            cpp_token = f"std::{_RENAMED_FUNCTIONS[token]}"

        elif token in _CUSTOM_FUNCTIONS:
            cpp_token = f"sm::{_CUSTOM_FUNCTIONS[token]}"

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
                        if state_var.id == var:
                            rate = state_var.derivative_id
                            break
            cpp_formula = cpp_formula.replace(f"rateOf({var})", rate)
    return cpp_formula
