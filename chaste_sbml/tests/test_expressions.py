"""Tests for SBML expression/AST -> C++ formula translation."""

import libsbml
import pytest
from libsbml import (
    AST_FUNCTION_COS,
    AST_FUNCTION_SIN,
    AST_NAME,
    AST_NAME_AVOGADRO,
    AST_PLUS,
    AST_REAL,
    ASTNode,
    formulaToL3String,
    parseL3Formula,
)

from chaste_sbml._config import VarType
from chaste_sbml._expressions import (
    collect_ast_names,
    formula_to_string,
    map_functions,
    resolve_constants,
    search_ast_type,
    strip_ast_units,
    substitute_ast_names,
)
from chaste_sbml._records import LocalParameter, StateVariable


def _cpp(formula, variable_types=None, state_variables=None, local_parameters=None):
    """Render an L3 formula string to C++ via ``formula_to_string`` (test helper)."""
    return formula_to_string(
        parseL3Formula(formula),
        variable_types or {},
        state_variables=state_variables or [],
        local_parameters=local_parameters,
    )


# A characterisation battery captured from the pre-refactor implementation: renaming SBML symbols and
# functions on the AST then rendering with formulaToL3String must reproduce the old string output
# byte-for-byte. Covers arithmetic grouping, relationals, logical ops, power/root/log and numbers.
_PARAMS = {c: VarType.PARAMETER for c in "abcdxyzkSn"}
_SNAPSHOT = [
    ("a+b", "a + b"),
    ("a-b", "a - b"),
    ("a*b", "a * b"),
    ("a/b", "a / b"),
    ("-a", "-a"),
    ("-(a+b)", "-(a + b)"),
    ("a-b-c", "a - b - c"),
    ("a-b+c", "(a - b) + c"),
    ("a/b*c", "(a / b) * c"),
    ("a*b/c", "a * b / c"),
    ("a+b*c", "a + b * c"),
    ("(a+b)*c", "(a + b) * c"),
    ("a<b", "a < b"),
    ("a>=b", "a >= b"),
    ("a==b", "a == b"),
    ("a!=b", "a != b"),
    ("lt(a,b,c)", "sm::lt(a, b, c)"),
    ("geq(a,b,c)", "sm::geq(a, b, c)"),
    ("eq(a,b,c,d)", "sm::eq(a, b, c, d)"),
    ("and(a,b)", "a && b"),
    ("or(a,b,c)", "a || b || c"),
    ("not(a)", "!a"),
    ("xor(a,b)", "sm::xor_(a, b)"),
    ("xor(a>4, a<1, a<2)", "sm::xor_((a > 4.0), (a < 1.0), (a < 2.0))"),
    ("implies(a,b)", "sm::implies(a, b)"),
    # A logical operator keeps libsbml's parenthesisation of relational operands.
    ("implies(a>b, c<d)", "sm::implies((a > b), (c < d))"),
    ("arcsec(a)", "sm::asec(a)"),
    ("arcsech(a)", "sm::asech(a)"),
    ("piecewise(a,b,c)", "sm::piecewise(a, b, c)"),
    ("piecewise(a,b,c,d,x)", "sm::piecewise(a, b, c, d, x)"),
    ("min(a,b)", "sm::min(a, b)"),
    ("max(a,b,c)", "sm::max(a, b, c)"),
    ("quotient(a,b)", "sm::quotient(a, b)"),
    ("rem(a,b)", "std::fmod(a, b)"),
    ("factorial(n)", "sm::factorial(n)"),
    ("abs(a)", "std::fabs(a)"),
    ("ceiling(a)", "std::ceil(a)"),
    ("ln(a)", "std::log(a)"),
    ("sin(a)", "std::sin(a)"),
    ("cos(a)", "std::cos(a)"),
    ("exp(a)", "std::exp(a)"),
    ("sqrt(a)", "std::sqrt(a)"),
    ("cot(a)", "sm::cot(a)"),
    ("csc(a)", "sm::csc(a)"),
    ("sech(a)", "sm::sech(a)"),
    ("arcsin(a)", "std::asin(a)"),
    ("arccot(a)", "sm::acot(a)"),
    ("a^b", "std::pow(a, b)"),
    ("a^b^c", "std::pow(a, std::pow(b, c))"),
    ("(a+b)^2", "std::pow(a + b, 2.0)"),
    ("pow(a,b)", "std::pow(a, b)"),
    ("sqr(a)", "std::pow(a, 2.0)"),
    ("root(a)", "sm::root(a)"),
    ("root(2,a)", "std::sqrt(a)"),
    ("root(3,a)", "sm::root(3.0, a)"),
    ("root(2.0,a)", "sm::root(2.0, a)"),
    ("log(a)", "log10(a)"),
    ("log(2,a)", "sm::log(2.0, a)"),
    ("log(10,a)", "log10(a)"),
    ("k*sin(a+b)+2", "k * std::sin(a + b) + 2.0"),
    ("2", "2.0"),
    ("3.5", "3.5"),
    (".5", "0.5"),
    ("1e-5", "1e-5"),
    ("1e5", "1e5"),
    ("100000", "100000.0"),
    ("2^10", "std::pow(2.0, 10.0)"),
    ("sin(a)+abs(b)", "std::sin(a) + std::fabs(b)"),
    ("a/b + c*d - e", "a / b + c * d - e"),
    ("1/3", "1.0 / 3.0"),
]


@pytest.mark.parametrize(("formula", "expected"), _SNAPSHOT)
def test_formula_to_string_matches_reference(formula, expected):
    """Every battery formula renders to its captured C++ string (byte-for-byte parity)."""
    assert _cpp(formula, _PARAMS) == expected


@pytest.mark.parametrize(
    ("mathml", "expected"),
    [
        # A base-less MathML log/root exercises the "no explicit degree/base" branch.
        ("<apply><log/><ci>a</ci></apply>", "log10(a)"),
        ("<apply><root/><ci>a</ci></apply>", "std::sqrt(a)"),
        ("<apply><power/><ci>a</ci><cn>2</cn></apply>", "std::pow(a, 2.0)"),
    ],
)
def test_formula_to_string_mathml_forms(mathml, expected):
    """MathML power/root/log forms (distinct AST from the L3 parser) still render correctly."""
    math = libsbml.readMathMLFromString(f'<math xmlns="http://www.w3.org/1998/Math/MathML">{mathml}</math>')
    assert formula_to_string(math, {"a": VarType.PARAMETER}, state_variables=[]) == expected


@pytest.mark.parametrize(
    ("formula", "expected"),
    [
        ("a ^ b", "std::pow(a, b)"),
        ("2 ^ 3", "std::pow(2.0, 3.0)"),
        ("(a + b) ^ 2", "std::pow(a + b, 2.0)"),
        # Nested power fully converts (the old string rewrite left an inner ^ behind).
        ("a ^ b ^ c", "std::pow(a, std::pow(b, c))"),
    ],
)
def test_formula_to_string_power_uses_std_pow(formula, expected):
    """Power operators (including nested) become std::pow calls."""
    variable_types = {name: VarType.PARAMETER for name in ("a", "b", "c")}
    assert _cpp(formula, variable_types) == expected


def test_formula_to_string_functions_constants_and_power():
    """Renders SBML functions/constants and infix power as C++."""
    out = _cpp("k * pow(S, 2) + pi", {"k": VarType.PARAMETER, "S": VarType.STATE_VARIABLE})
    assert out == "k * std::pow(S, 2.0) + M_PI"


def test_formula_to_string_implies_maps_to_sbmlmath():
    """The MathML 'implies' logical operator maps to the sm::implies helper."""
    assert _cpp("implies(a, b)") == "sm::implies(a, b)"


@pytest.mark.parametrize(
    ("formula", "expected"),
    [
        # n-ary relationals (>2 operands) become sm:: calls with the correct semantics...
        ("lt(1, 2, 3)", "sm::lt(1.0, 2.0, 3.0)"),
        ("gt(2, 1, 2)", "sm::gt(2.0, 1.0, 2.0)"),
        ("leq(1, 2, 3, 4)", "sm::leq(1.0, 2.0, 3.0, 4.0)"),
        ("geq(2, 1, 2)", "sm::geq(2.0, 1.0, 2.0)"),
        ("eq(1, 1, 2)", "sm::eq(1.0, 1.0, 2.0)"),
        # ...while binary relationals stay as infix operators (unchanged output).
        ("a < b", "a < b"),
        ("a >= b", "a >= b"),
    ],
)
def test_formula_to_string_nary_relational_uses_sbmlmath(formula, expected):
    """n-ary relationals route to the variadic sm:: helpers; binary ones stay infix."""
    variable_types = {"a": VarType.PARAMETER, "b": VarType.PARAMETER}
    assert _cpp(formula, variable_types) == expected


@pytest.mark.parametrize(
    ("mathml", "expected"),
    [
        # Degenerate (fewer-operand) operators libsbml renders as calls, not infix.
        ("<apply><plus/></apply>", "sm::plus()"),
        ("<apply><times/></apply>", "sm::times()"),
        ("<apply><plus/><cn>1</cn></apply>", "sm::plus(1.0)"),
        ("<apply><and/></apply>", "sm::and_()"),
        ("<apply><or/></apply>", "sm::or_()"),
        ("<apply><not/></apply>", "sm::not_()"),
    ],
)
def test_formula_to_string_degenerate_operators_render_as_calls(mathml, expected):
    """A plus/times/and/or/not with too few operands maps to its sm:: helper call."""
    math = libsbml.readMathMLFromString(f'<math xmlns="http://www.w3.org/1998/Math/MathML">{mathml}</math>')
    assert formula_to_string(math, {}, state_variables=[]) == expected


def test_formula_to_string_unknown_function_is_left_unchanged():
    """A call to the model's own function definition keeps its name (not std::/sm:: mapped)."""
    assert _cpp("myFunc(a, b)", {"a": VarType.PARAMETER, "b": VarType.PARAMETER}) == "myFunc(a, b)"


def test_formula_to_string_variable_named_like_function_is_distinct():
    """A parameter spelled like a function stays a bare name; a real call is qualified.

    The old string token-loop mis-handled this: a parameter named ``sin`` suppressed the ``std::``
    prefix on a genuine ``sin()`` call too. Mapping on the AST (by node kind) keeps them independent.
    """
    out = _cpp("sin + sin(x)", {"sin": VarType.PARAMETER, "x": VarType.PARAMETER})
    assert out == "sin + std::sin(x)"


@pytest.mark.parametrize(
    ("name", "args"),
    [
        ("sin", "<ci>x</ci>"),  # spelled like an _UNCHANGED_FUNCTIONS built-in (std::sin)
        ("min", "<ci>a</ci><ci>b</ci>"),  # spelled like a _CUSTOM_FUNCTIONS built-in (sm::min)
        ("log", "<ci>x</ci>"),  # spelled like the log built-in (log10/sm::log)
    ],
)
def test_formula_to_string_model_function_spelled_like_builtin_is_left_unchanged(name, args):
    """A call to a model's own function definition spelled like a built-in keeps its own name.

    Such a call is a generic ``AST_FUNCTION`` node (see ``_names.py``); only libsbml's *typed*
    built-in function nodes are name-mapped, so the call resolves to the model's function, not the
    C++ built-in of the same spelling.
    """
    mathml = f'<math xmlns="http://www.w3.org/1998/Math/MathML"><apply><ci>{name}</ci>{args}</apply></math>'
    math = libsbml.readMathMLFromString(mathml)
    rendered = args.replace("<ci>", "").replace("</ci>", ", ").rstrip(", ")
    assert formula_to_string(math, {}, state_variables=[]) == f"{name}({rendered})"


def test_formula_to_string_delay_is_unsupported():
    """A delay function raises rather than emitting invalid C++."""
    with pytest.raises(NotImplementedError, match="delay"):
        _cpp("delay(x, tau)", {"x": VarType.PARAMETER, "tau": VarType.PARAMETER})


def test_formula_to_string_rateof_state_variable_uses_derivative():
    """rateOf(state var) resolves to that variable's derivative id."""
    out = formula_to_string(
        parseL3Formula("rateOf(X)"),
        variable_types={"X": VarType.STATE_VARIABLE},
        state_variables=[
            StateVariable(index=0, id="X", derivative_id="d_X_dt", label="X", initial_value=None, units="")
        ],
    )
    assert out == "d_X_dt"


def test_formula_to_string_rateof_local_parameter_is_zero():
    """rateOf of a kinetic-law local parameter (constant) is zero, not a global's derivative."""
    out = formula_to_string(
        parseL3Formula("rateOf(k)"),
        variable_types={"k": VarType.STATE_VARIABLE},
        state_variables=[
            StateVariable(index=0, id="k", derivative_id="d_k_dt", label="", initial_value=None, units="")
        ],
        local_parameters=[LocalParameter(id="k", label="", value="1")],
    )
    assert out == "0.0"


def test_collect_ast_names_gathers_variable_identifiers():
    """collect_ast_names accumulates the variable names, excluding function names."""
    names = set()
    collect_ast_names(parseL3Formula("a + b * sin(c)"), names)
    assert names == {"a", "b", "c"}


def test_collect_ast_names_ignores_none():
    """A None node contributes nothing to the accumulator."""
    names = {"pre_existing"}
    collect_ast_names(None, names)
    assert names == {"pre_existing"}


def test_substitute_ast_names_replaces_named_nodes():
    """Each named node is replaced by a copy of its replacement AST."""
    result = substitute_ast_names(parseL3Formula("a + b"), {"a": parseL3Formula("x * 2")})
    assert formulaToL3String(result) == "x * 2 + b"


def test_substitute_ast_names_replaces_a_bare_name():
    """A replacement applied to the named target node itself returns the replacement AST."""
    result = substitute_ast_names(parseL3Formula("a"), {"a": parseL3Formula("z + 1")})
    assert formulaToL3String(result) == "z + 1"


def test_substitute_ast_names_leaves_original_unchanged():
    """Substitution returns a new AST, leaving the input untouched."""
    original = parseL3Formula("a + b")
    substitute_ast_names(original, {"a": parseL3Formula("z")})
    assert formulaToL3String(original) == "a + b"


def test_resolve_constants_renames_nested_symbol():
    """An avogadro csymbol anywhere in the tree becomes a plain name with its C++ value."""
    one = ASTNode(AST_REAL)
    one.setValue(1.0)
    avogadro = ASTNode(AST_NAME_AVOGADRO)
    expr = ASTNode(AST_PLUS)
    expr.addChild(one)
    expr.addChild(avogadro)

    resolve_constants(expr)

    renamed = expr.getChild(1)
    assert renamed.getType() == AST_NAME
    assert renamed.getName() == "sm::AVOGADRO"


def test_resolve_constants_leaves_other_names():
    """A non-symbol name node is left untouched."""
    name = ASTNode(AST_NAME)
    name.setName("x")
    resolve_constants(name)
    assert name.getName() == "x"


def test_map_functions_leaves_plain_operators():
    """Arithmetic operators are left for the formatter (map_functions renames only calls)."""
    node = parseL3Formula("a + b * c")
    map_functions(node)
    assert formulaToL3String(node) == "a + b * c"


@pytest.mark.parametrize(
    ("mathml_symbol", "expected"),
    [
        ("<pi/>", "M_PI"),
        ("<exponentiale/>", "M_E"),
        ("<infinity/>", "std::numeric_limits<double>::infinity()"),
        ("<notanumber/>", "std::numeric_limits<double>::quiet_NaN()"),
    ],
)
def test_formula_to_string_math_constants(mathml_symbol, expected):
    """SBML math constants map to their C++ form."""
    math = libsbml.readMathMLFromString(f'<math xmlns="http://www.w3.org/1998/Math/MathML">{mathml_symbol}</math>')
    assert formula_to_string(math, {}, state_variables=[]) == expected


def test_formula_to_string_negative_infinity():
    """A negative-infinity literal maps to the negated C++ infinity."""
    math = libsbml.readMathMLFromString(
        '<math xmlns="http://www.w3.org/1998/Math/MathML"><apply><minus/><infinity/></apply></math>'
    )
    assert formula_to_string(math, {}, state_variables=[]) == "-std::numeric_limits<double>::infinity()"


@pytest.mark.parametrize(
    ("mathml_symbol", "param", "expected"),
    [
        ("<infinity/>", "INF", "std::numeric_limits<double>::infinity()"),
        ("<notanumber/>", "NaN", "std::numeric_limits<double>::quiet_NaN()"),
    ],
)
def test_formula_to_string_infinity_nan_distinct_from_parameter(mathml_symbol, param, expected):
    """The <infinity/>/<notanumber/> constants stay distinct from a parameter spelled INF/NaN.

    libsbml stores these constants as real nodes that render as ``INF``/``NaN`` -- the same text as
    a reference to a parameter so named -- so the constant must be resolved from the AST, not its
    rendered spelling (test-suite cases 1811/1813).
    """
    math = libsbml.readMathMLFromString(f'<math xmlns="http://www.w3.org/1998/Math/MathML">{mathml_symbol}</math>')
    assert formula_to_string(math, {param: VarType.PARAMETER}, state_variables=[]) == expected


def test_formula_to_string_pi_constant_and_parameter_are_distinct():
    """The <pi/> constant maps to M_PI while a same-named parameter keeps its own name."""
    math = libsbml.readMathMLFromString(
        '<math xmlns="http://www.w3.org/1998/Math/MathML"><apply><plus/><ci>pi</ci><pi/></apply></math>'
    )
    assert formula_to_string(math, {"pi": VarType.PARAMETER}, state_variables=[]) == "pi + M_PI"


def test_formula_to_string_time_symbol_and_named_variable_are_distinct():
    """The time csymbol maps to `time`; a plain variable named `t` is left as-is (not time)."""
    assert _cpp("time") == "time"
    assert _cpp("t", {"t": VarType.PARAMETER}) == "t"


def test_strip_ast_units_removes_units_recursively():
    """Units annotations are stripped from a node and its children."""
    node = parseL3Formula("2 * 3")
    node.getChild(0).setUnits("mole")
    strip_ast_units(node)
    assert not node.getChild(0).isSetUnits()


def test_search_ast_type_finds_nested_node():
    """search_ast_type reports whether a node type occurs anywhere in the tree."""
    math = parseL3Formula("k + sin(x)")
    assert search_ast_type(math, AST_FUNCTION_SIN) is True
    assert search_ast_type(math, AST_FUNCTION_COS) is False


def test_search_ast_type_handles_none():
    """A None root yields False."""
    assert search_ast_type(None, AST_FUNCTION_SIN) is False
