"""Translation of SBML MathML/AST expressions into C++ formula strings.

Extracted from ChasteSbmlModel so the expression-formatting logic lives on its own. The
functions are pure apart from the libsbml AST objects they operate on: ``formula_to_string`` is
given the model's variable-type map and state-variable list as arguments rather than reaching
into the model.
"""

import math
import re
from typing import TYPE_CHECKING, Optional

from libsbml import (
    AST_CONSTANT_E,
    AST_CONSTANT_PI,
    AST_FUNCTION,
    AST_FUNCTION_DELAY,
    AST_FUNCTION_POWER,
    AST_FUNCTION_RATE_OF,
    AST_NAME,
    AST_NAME_AVOGADRO,
    AST_NAME_TIME,
    AST_POWER,
    AST_REAL,
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


# Placeholder ids the SBML math symbols are renamed to. A parameter may share a symbol's spelling
# (e.g. a parameter named 'pi' or 'time'), and ``formulaToL3String`` renders the symbol and the
# parameter identically; renaming the symbol node to a unique placeholder keeps them distinct, so
# the symbol maps to its C++ value (below) while the parameter keeps its own name.
AVOGADRO_PLACEHOLDER = f"{CHASTE_PREFIX}{PREFIX_SEP}avogadro"
PI_PLACEHOLDER = f"{CHASTE_PREFIX}{PREFIX_SEP}pi"
E_PLACEHOLDER = f"{CHASTE_PREFIX}{PREFIX_SEP}exponentiale"
TIME_PLACEHOLDER = f"{CHASTE_PREFIX}{PREFIX_SEP}time"
INFINITY_PLACEHOLDER = f"{CHASTE_PREFIX}{PREFIX_SEP}infinity"
NEG_INFINITY_PLACEHOLDER = f"{CHASTE_PREFIX}{PREFIX_SEP}neg_infinity"
NAN_PLACEHOLDER = f"{CHASTE_PREFIX}{PREFIX_SEP}nan"

_SYMBOL_PLACEHOLDERS = {
    AST_NAME_AVOGADRO: AVOGADRO_PLACEHOLDER,
    AST_CONSTANT_PI: PI_PLACEHOLDER,
    AST_CONSTANT_E: E_PLACEHOLDER,
    AST_NAME_TIME: TIME_PLACEHOLDER,
}


def replace_constant_symbols(node: "ASTNode") -> None:
    """Recursively rename SBML math-symbol nodes (avogadro, pi, exponentiale, time, inf, NaN).

    Each of these symbols renders identically to a same-named parameter, so it is renamed to a
    unique placeholder that maps to its C++ value, keeping it distinct from a parameter of the same
    spelling. Handling this at the AST level (rather than string-matching the rendered name) means a
    real variable named ``t``, ``pi`` etc. is never mistaken for the symbol. The MathML constants
    ``<infinity/>`` and ``<notanumber/>`` are especially treacherous: libsbml stores them as ordinary
    real nodes that render as ``INF``/``NaN``, indistinguishable from a reference to a parameter so
    named (see test-suite cases 1811/1813), so they are matched here by their non-finite value.

    :param node: The root AST node.
    """
    placeholder = _SYMBOL_PLACEHOLDERS.get(node.getType())
    if placeholder is not None:
        node.setType(AST_NAME)
        node.setName(placeholder)
    elif node.getType() == AST_REAL:
        value = node.getReal()
        if math.isnan(value):
            node.setType(AST_NAME)
            node.setName(NAN_PLACEHOLDER)
        elif math.isinf(value):
            node.setType(AST_NAME)
            node.setName(NEG_INFINITY_PLACEHOLDER if value < 0 else INFINITY_PLACEHOLDER)
    for i in range(node.getNumChildren()):
        replace_constant_symbols(node.getChild(i))


def resolve_rate_of(
    node: "ASTNode",
    variable_types: dict,
    state_variables: list["StateVariable"],
    local_param_ids: set,
) -> None:
    """Recursively resolve ``rateOf(x)`` nodes to a state variable's derivative id, or ``0``.

    ``rateOf`` of a state variable is that variable's time derivative (its ``d_<id>_dt`` id); of
    anything constant -- a parameter, or a kinetic-law local parameter that shadows a global -- it
    is zero. Resolving on the AST (rather than by regex on the rendered string) is robust to the
    argument's spelling and to surrounding parentheses.

    :param node: The root AST node.
    :param variable_types: Mapping of model id to :class:`VarType`.
    :param state_variables: The model's state variables (to look up a derivative id).
    :param local_param_ids: Ids of local parameters in scope (constant, so their rate is zero).
    """
    if node.getType() == AST_FUNCTION_RATE_OF and node.getNumChildren() >= 1:
        var = node.getChild(0).getName() or ""
        derivative = None
        if var not in local_param_ids and variable_types.get(var, VarType.UNKNOWN) == VarType.STATE_VARIABLE:
            derivative = next((sv.derivative_id for sv in state_variables if sv.id == var), None)
        while node.getNumChildren() > 0:
            node.removeChild(0)
        if derivative is not None:
            node.setType(AST_NAME)
            node.setName(derivative)
        else:
            node.setType(AST_REAL)
            node.setValue(0.0)
        return
    for i in range(node.getNumChildren()):
        resolve_rate_of(node.getChild(i), variable_types, state_variables, local_param_ids)


def strip_ast_units(node: "ASTNode") -> None:
    """Recursively strip units annotations from AST nodes.

    SBML allows numeric literals to carry units annotations (e.g. ``<cn sbml:units="mole">``).
    ``formulaToL3String`` includes these in its output (e.g. ``0.00015 mole``), which is not
    valid C++. Stripping them here is safe because the units are not used.

    :param node: The root AST node.
    """
    # TODO: Handle units properly
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

    SBML MathML relationals are n-ary: ``lt(a, b, c)`` means ``a < b < c``
    (i.e. ``a < b`` and ``b < c``), and ``eq(a, b, c)`` means all equal. However,
    ``formulaToL3String`` renders them as a chained infix expression (``a < b < c``),
    which C++ mis-evaluates left-to-right as  ``(a < b) < c``. Converting each
    such node to a function call (rendered ``lt(a, b, c)``, later mapped to the
    variadic ``sm::lt`` helper that implements the correct semantics) fixes this.
    Binary relationals are left as infix operators, so existing generated code
    is unchanged. ``neq`` is binary-only in MathML, so it is never n-ary and is
    not rewritten here.

    :param node: The root AST node to rewrite.
    """
    if node is None:  # pragma: no cover - defensive base case; callers pass a real AST root
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

    ``formulaToL3String`` renders a power as the infix ``a ^ b``, which is not
    valid C++ (``^`` is bitwise xor). Converting each power node to a ``pow``
    call before stringifying (rendered ``pow(a, b)``, later mapped to ``std::pow``
    by the function-name token map) handles nesting correctly. Both power node
    types are handled: ``AST_POWER`` (from ``^``/``pow`` in L3 infix, and from
    MathML ``<power>`` that libsbml normalises) and ``AST_FUNCTION_POWER``
    (from other MathML ``<power>`` forms).

    :param node: The root AST node to rewrite.
    """
    if node is None:  # pragma: no cover - defensive base case; callers pass a real AST root
        return
    if node.getType() in (AST_POWER, AST_FUNCTION_POWER):
        node.setType(AST_FUNCTION)
        node.setName("pow")
    for i in range(node.getNumChildren()):
        rewrite_power(node.getChild(i))


# SBML symbolic constants replaced with their C++ equivalents, keyed by the placeholder each symbol
# is renamed to (see replace_constant_symbols). Every one is matched at the AST level and renamed to
# its placeholder, so a same-named parameter (e.g. one literally called pi, inf or NaN) is never
# mistaken for the constant. (true/false deliberately absent.)
_CONSTANTS = {
    AVOGADRO_PLACEHOLDER: "sm::AVOGADRO",
    PI_PLACEHOLDER: "M_PI",
    E_PLACEHOLDER: "M_E",
    TIME_PLACEHOLDER: "time",
    INFINITY_PLACEHOLDER: "std::numeric_limits<double>::infinity()",
    NEG_INFINITY_PLACEHOLDER: "-std::numeric_limits<double>::infinity()",
    NAN_PLACEHOLDER: "std::numeric_limits<double>::quiet_NaN()",
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
    "implies": "implies",
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

    :param math: The AST math formula. It is not mutated.
    :param variable_types: Mapping of model id to :class:`VarType`, used to tell
        model variables apart from SBML constants/functions of the same spelling.
    :param state_variables: The model's state variables, used to resolve ``rateOf``
        to a state variable's derivative id.
    :param local_parameters: Local parameters in scope (e.g. a reaction's
        kinetic-law parameters). These shadow global symbols of the same name
        and are constant, so ``rateOf`` applied to one is zero.
    :return: The equivalent C++ string.
    """
    if search_ast_type(math, AST_FUNCTION_DELAY):
        raise NotImplementedError("SBML function not supported: 'delay'.")

    # Make a copy so the AST is not mutated.
    math = math.deepCopy()
    local_param_ids = {param.id for param in (local_parameters or [])}

    strip_ast_units(math)
    # Modify the AST, renaming the avogadro/pi/exponentiale/time symbol nodes to
    # placeholders (so each stays distinct from a same-named parameter), and
    # resolve rateOf to a derivative id (or 0).
    replace_constant_symbols(math)
    resolve_rate_of(math, variable_types, state_variables, local_param_ids)
    # Route operators that formulaToL3String would render as C++-invalid infix
    # to function-call form: n-ary relationals (a < b < c) to the SBMLMath
    # helpers (sm::lt(a, b, c)), and the power operator (a ^ b) to std::pow().
    # Binary relationals stay as infix operators (a < b).
    rewrite_nary_relational(math)
    rewrite_power(math)
    formula = formulaToL3String(math)

    # Convert all integer literals to doubles to fix integer division. The lookbehinds
    # keep the digits of a scientific-notation exponent (e.g. 1e-5, 1e+5) intact.
    formula = re.sub(r"(?<!\.)(?<![eE][-+])\b[0-9]+\b(?!\.)", lambda x: f"{x[0]}.0", formula)

    tokens = re.findall(r"\w+|\W+", formula)

    cpp_tokens = []
    for token in tokens:
        cpp_token = token

        # Replace function names and constants, but only when the token is not
        # an actual model variable or a local parameter (e.g. a local parameter
        # named "inf" must stay local).
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
    return cpp_formula
