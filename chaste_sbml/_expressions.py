"""Translation of SBML MathML/AST expressions into C++ formula strings.

The functions are pure apart from the libsbml AST objects they operate on:
``formula_to_string`` is given the model's variable-type map and state-variable
list as arguments rather than reaching into the model.

Rendering strategy: libsbml's ``formulaToL3String`` does all the structural work
(operator precedence, parenthesization and numeric-literal formatting) and we do
not reimplement any of it. Instead the mapping from SBML symbols/functions to
their C++ spellings is applied *on the AST* (``resolve_constants``,
``map_functions``) by renaming nodes to their final C++ names, which
``formulaToL3String`` then renders verbatim. Working on the AST, where a node's
kind is known, means a variable named ``sin`` or ``pi`` is never confused with
the function/constant of the same spelling, unlike a string search of the
rendered output. The only string post-processing is turning integer literals
into doubles (to avoid integer division).
"""

import math
import re
from typing import TYPE_CHECKING, Optional

from libsbml import (
    AST_CONSTANT_E,
    AST_CONSTANT_PI,
    AST_FUNCTION,
    AST_FUNCTION_DELAY,
    AST_FUNCTION_LOG,
    AST_FUNCTION_POWER,
    AST_FUNCTION_RATE_OF,
    AST_FUNCTION_ROOT,
    AST_INTEGER,
    AST_NAME,
    AST_NAME_AVOGADRO,
    AST_NAME_TIME,
    AST_PLUS,
    AST_POWER,
    AST_REAL,
    AST_RELATIONAL_EQ,
    AST_RELATIONAL_GEQ,
    AST_RELATIONAL_GT,
    AST_RELATIONAL_LEQ,
    AST_RELATIONAL_LT,
    AST_TIMES,
    formulaToL3String,
)

from ._config import VarType

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


# Non-finite literals rendered as C++ (matched by value, see resolve_constants).
_INFINITY_CPP = "std::numeric_limits<double>::infinity()"
_NEG_INFINITY_CPP = "-" + _INFINITY_CPP
_NAN_CPP = "std::numeric_limits<double>::quiet_NaN()"

# SBML math-symbol node types renamed directly to their C++ value. Matched by
# node *type*, so a parameter that happens to share a symbol's spelling (e.g.
# one named ``pi`` or ``time``) is an ordinary AST_NAME node and is never
# mistaken for the symbol.
_SYMBOL_CPP = {
    AST_NAME_AVOGADRO: "sm::AVOGADRO",
    AST_CONSTANT_PI: "M_PI",
    AST_CONSTANT_E: "M_E",
    AST_NAME_TIME: "time",
}


def resolve_constants(node: "ASTNode") -> None:
    """Recursively rename SBML math-symbol nodes to their C++ value, in place.

    Handles avogadro, pi, exponentiale and the time csymbol (renamed by node
    type), and the MathML constants ``<infinity/>`` and ``<notanumber/>``. The
    latter are especially treacherous: libsbml stores them as ordinary real
    nodes that render as ``INF``/``NaN``, indistinguishable from a reference to
    a parameter so named (see test-suite cases 1811/1813), so they are matched
    here by their non-finite value. Each node is renamed to a plain name holding
    the C++ text, which ``formulaToL3String`` then renders verbatim; a
    same-named parameter is an ordinary name node and is left untouched.

    :param node: The root AST node.
    """
    cpp = _SYMBOL_CPP.get(node.getType())
    if cpp is not None:
        node.setType(AST_NAME)
        node.setName(cpp)
    elif node.getType() == AST_REAL:
        value = node.getReal()
        if math.isnan(value):
            node.setType(AST_NAME)
            node.setName(_NAN_CPP)
        elif math.isinf(value):
            node.setType(AST_NAME)
            node.setName(_NEG_INFINITY_CPP if value < 0 else _INFINITY_CPP)
    for i in range(node.getNumChildren()):
        resolve_constants(node.getChild(i))


def resolve_rate_of(
    node: "ASTNode",
    variable_types: dict,
    state_variables: list["StateVariable"],
    local_param_ids: set,
) -> None:
    """Recursively resolve ``rateOf(x)`` nodes to a state variable's derivative id, or ``0``.

    ``rateOf`` of a state variable is that variable's time derivative (``d_<id>_dt``).
    ``rateOf`` of anything constant (e.g. a parameter) is zero.
    Resolving on the AST (rather than on the rendered string) is robust to the
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

    SBML allows numeric literals to carry units annotations (e.g.
    ``<cn sbml:units="mole">``). ``formulaToL3String`` includes these in its
    output (e.g. ``0.00015 mole``), which is not valid C++. Stripping them here
    is safe because the units are not used.

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

# SBML functions with custom implementations in SbmlMath (sm::), rendered as
# calls. The logical operators (and/or/not/xor/implies), n-ary relationals
# (lt/gt/leq/geq/eq) and plus/times are handled by their own branches in
# map_functions, so they are not listed here.
_CUSTOM_FUNCTIONS = {
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
    "factorial": "factorial",
    "log": "log",
    "max": "max",
    "min": "min",
    "piecewise": "piecewise",
    "quotient": "quotient",
    "root": "root",
    "sec": "sec",
    "sech": "sech",
    "sqr": "sqr",
}

# SBML logical operators mapped to their C++ spelling. libsbml renders these
# infix (&& || !) when they have enough operands and as a call otherwise, so
# they are renamed in place (keeping the logical node type) to preserve that
# arity-based choice and libsbml's parenthesisation of the operands.
_LOGICAL_CPP = {
    "and": "sm::and_",
    "or": "sm::or_",
    "not": "sm::not_",
    "xor": "sm::xor_",
    "implies": "sm::implies",
}

# TODO: From SBML Level 3 upwards, log defaults to base 10.
# SBML versions lower than 3 default to base e.
# See https://sbml.org/software/libsbml/5.18.0/docs/formatted/python-api/namespacelibsbml.html#a8e96a5a70569ae32655c6302638f6dc3  # noqa: B950


def _rename_call(node: "ASTNode", cpp_name: str) -> None:
    """Turn a node into a plain function call rendered with the given C++ name.

    :param node: The node to rewrite in place.
    :param cpp_name: The C++ function name (e.g. ``std::sin``), rendered verbatim by libsbml.
    """
    node.setType(AST_FUNCTION)
    node.setName(cpp_name)


def _map_root(node: "ASTNode") -> None:
    """Map an SBML ``root`` node to ``std::sqrt`` (square roots) or ``sm::root`` (otherwise).

    libsbml renders a two-argument root with an integer degree of 2 as ``sqrt``;
    every other form (a general degree, a real ``2.0``, or a lone radicand) is a
    general n-th root. We reproduce that split so the output matches libsbml's
    own, dropping the redundant degree child for ``sqrt``.

    :param node: The AST_FUNCTION_ROOT node to rewrite in place.
    """
    degree = node.getChild(0) if node.getNumChildren() >= 2 else None
    if degree is not None and degree.getType() == AST_INTEGER and degree.getInteger() == 2:
        node.removeChild(0)
        _rename_call(node, "std::sqrt")
    else:
        _rename_call(node, "sm::root")


def _map_log(node: "ASTNode") -> None:
    """Map an SBML ``log`` node to ``log10`` (base 10) or ``sm::log`` (a general base).

    libsbml renders a base-10 log as ``log10``; ``log`` with any other explicit
    base is a two-argument general log. We reproduce that split, dropping the
    redundant base-10 child.

    :param node: The AST_FUNCTION_LOG node to rewrite in place.
    """
    base = node.getChild(0) if node.getNumChildren() >= 2 else None
    if base is None or (base.getType() == AST_INTEGER and base.getInteger() == 10):
        if base is not None:
            node.removeChild(0)
        _rename_call(node, "log10")
    else:
        _rename_call(node, "sm::log")


def map_functions(node: "ASTNode") -> None:
    """Recursively rename SBML function/operator nodes to their C++ spelling, in place.

    Most operators libsbml already renders as infix C++ (``+ - * /``, ``&& || !``,
    binary relationals) and are left alone. Handled here:

    - power (``^``, MathML ``<power>``, ``sqr``) -> ``std::pow``
    - ``root``/``log`` -> ``std::sqrt``/``log10`` or ``sm::root``/``sm::log`` (see helpers)
    - degenerate ``plus``/``times`` (fewer than two operands, which libsbml renders as a call rather
      than infix) -> ``sm::plus``/``sm::times``
    - n-ary (>2 operand) relationals -> the variadic ``sm::`` helper (``a < b < c`` -> ``sm::lt``);
      binary relationals stay infix
    - logical operators -> renamed in place to ``sm::and_``/``or_``/``not_``/``xor_``/``implies``,
      keeping the logical node type so libsbml still picks infix vs call by arity (and parenthesises
      the operands the same way)
    - other functions -> ``std::``/``sm::`` per the name maps; an unknown name (a model's own
      function definition) is left as-is

    :param node: The root AST node to rewrite in place.
    """
    node_type = node.getType()
    if node_type in (AST_POWER, AST_FUNCTION_POWER):
        _rename_call(node, "std::pow")
    elif node_type == AST_FUNCTION_ROOT:
        _map_root(node)
    elif node_type == AST_FUNCTION_LOG:
        _map_log(node)
    elif node_type in (AST_PLUS, AST_TIMES) and node.getNumChildren() < 2:
        _rename_call(node, "sm::plus" if node_type == AST_PLUS else "sm::times")
    elif node.isRelational() and node.getNumChildren() > 2:
        _rename_call(node, f"sm::{_NARY_RELATIONAL_NAMES[node_type]}")
    elif node.isLogical() and node.getName() in _LOGICAL_CPP:
        node.setName(_LOGICAL_CPP[node.getName()])
    elif node.isFunction():
        name = node.getName()
        if name in _UNCHANGED_FUNCTIONS:
            _rename_call(node, f"std::{name}")
        elif name in _RENAMED_FUNCTIONS:
            _rename_call(node, f"std::{_RENAMED_FUNCTIONS[name]}")
        elif name in _CUSTOM_FUNCTIONS:
            _rename_call(node, f"sm::{_CUSTOM_FUNCTIONS[name]}")
        # An unknown name is a model's own function definition; leave it unchanged.

    for i in range(node.getNumChildren()):
        map_functions(node.getChild(i))


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
    # Rename symbol and function nodes to their final C++ spelling on the AST, so formulaToL3String
    # renders valid C++ directly (and a same-named variable is never mistaken for a symbol/function).
    resolve_constants(math)
    resolve_rate_of(math, variable_types, state_variables, local_param_ids)
    map_functions(math)

    formula = formulaToL3String(math)

    # Convert all integer literals to doubles to fix integer division. The lookbehinds
    # keep the digits of a scientific-notation exponent (e.g. 1e-5, 1e+5) intact.
    formula = re.sub(r"(?<!\.)(?<![eE][-+])\b[0-9]+\b(?!\.)", lambda x: f"{x[0]}.0", formula)

    return formula
