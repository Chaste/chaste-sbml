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
    AVOGADRO_PLACEHOLDER,
    collect_ast_names,
    formula_to_string,
    replace_constant_symbols,
    rewrite_nary_relational,
    rewrite_power,
    search_ast_type,
    strip_ast_units,
    substitute_ast_names,
)
from chaste_sbml._records import StateVariable


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
    assert formula_to_string(parseL3Formula(formula), variable_types, state_variables=[]) == expected


def test_rewrite_power_converts_only_power_nodes():
    """rewrite_power converts AST_POWER nodes to pow() and leaves other nodes alone."""
    node = parseL3Formula("a ^ b ^ c")
    rewrite_power(node)
    assert formulaToL3String(node) == "pow(a, pow(b, c))"

    unaffected = parseL3Formula("a * b + c")
    rewrite_power(unaffected)
    assert formulaToL3String(unaffected) == "a * b + c"


def test_rewrite_power_converts_mathml_function_power():
    """A MathML <power> (AST_FUNCTION_POWER, distinct from L3 AST_POWER) is also converted.

    Guards against only handling the L3-parsed power node type: SBML MathML powers parse to a
    different AST type that ``formulaToL3String`` still renders as ``^``.
    """
    mathml = '<math xmlns="http://www.w3.org/1998/Math/MathML">' "<apply><power/><ci>a</ci><cn>2</cn></apply></math>"
    node = libsbml.readMathMLFromString(mathml)
    rewrite_power(node)
    assert formulaToL3String(node) == "pow(a, 2)"


def test_formula_to_string_functions_constants_and_power():
    """Renders SBML functions/constants and infix power as C++."""
    math = parseL3Formula("k * pow(S, 2) + pi")
    out = formula_to_string(
        math, variable_types={"k": VarType.PARAMETER, "S": VarType.STATE_VARIABLE}, state_variables=[]
    )
    assert out == "k * std::pow(S, 2.0) + M_PI"


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
    assert formula_to_string(parseL3Formula(formula), variable_types, state_variables=[]) == expected


def test_rewrite_nary_relational_leaves_binary_untouched():
    """rewrite_nary_relational converts only >2-operand relationals, leaving binary ones as-is."""
    binary = parseL3Formula("a < b")
    rewrite_nary_relational(binary)
    assert formulaToL3String(binary) == "a < b"

    nary = parseL3Formula("lt(a, b, c)")
    rewrite_nary_relational(nary)
    assert formulaToL3String(nary) == "lt(a, b, c)"


def test_formula_to_string_rateof_state_variable_uses_derivative():
    """rateOf(state var) resolves to that variable's derivative id."""
    math = parseL3Formula("rateOf(X)")
    out = formula_to_string(
        math,
        variable_types={"X": VarType.STATE_VARIABLE},
        state_variables=[
            StateVariable(index=0, id="X", derivative_id="d_X_dt", label="X", initial_value=None, units="")
        ],
    )
    assert out == "d_X_dt"


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


def test_replace_constant_symbols_renames_nested_node():
    """An avogadro csymbol anywhere in the tree becomes a plain name with its placeholder id."""
    one = ASTNode(AST_REAL)
    one.setValue(1.0)
    avogadro = ASTNode(AST_NAME_AVOGADRO)
    expr = ASTNode(AST_PLUS)
    expr.addChild(one)
    expr.addChild(avogadro)

    replace_constant_symbols(expr)

    renamed = expr.getChild(1)
    assert renamed.getType() == AST_NAME
    assert renamed.getName() == AVOGADRO_PLACEHOLDER


def test_replace_constant_symbols_leaves_other_names():
    """A non-symbol name node is left untouched."""
    name = ASTNode(AST_NAME)
    name.setName("x")
    replace_constant_symbols(name)
    assert name.getName() == "x"


@pytest.mark.parametrize(
    ("mathml_symbol", "expected"),
    [
        ("<pi/>", "M_PI"),
        ("<exponentiale/>", "M_E"),
        ("<infinity/>", "std::numeric_limits<double>::infinity()"),
        ("<notanumber/>", "NAN"),
    ],
)
def test_formula_to_string_math_constants(mathml_symbol, expected):
    """SBML math constants (incl. the capitalised INF/NaN libsbml emits) map to their C++ form."""
    math = libsbml.readMathMLFromString(f'<math xmlns="http://www.w3.org/1998/Math/MathML">{mathml_symbol}</math>')
    assert formula_to_string(math, {}, state_variables=[]) == expected


def test_formula_to_string_pi_constant_and_parameter_are_distinct():
    """The <pi/> constant maps to M_PI while a same-named parameter keeps its own name."""
    math = libsbml.readMathMLFromString(
        '<math xmlns="http://www.w3.org/1998/Math/MathML"><apply><plus/><ci>pi</ci><pi/></apply></math>'
    )
    assert formula_to_string(math, {"pi": VarType.PARAMETER}, state_variables=[]) == "pi + M_PI"


def test_formula_to_string_time_symbol_and_named_variable_are_distinct():
    """The time csymbol maps to `time`; a plain variable named `t` is left as-is (not time)."""
    assert formula_to_string(parseL3Formula("time"), {}, state_variables=[]) == "time"
    assert formula_to_string(parseL3Formula("t"), {"t": VarType.PARAMETER}, state_variables=[]) == "t"


def test_formula_to_string_rateof_local_parameter_is_zero():
    """rateOf of a kinetic-law local parameter (constant) is zero, not a global's derivative."""
    from chaste_sbml._records import LocalParameter

    out = formula_to_string(
        parseL3Formula("rateOf(k)"),
        variable_types={"k": VarType.STATE_VARIABLE},
        state_variables=[
            StateVariable(index=0, id="k", derivative_id="d_k_dt", label="", initial_value=None, units="")
        ],
        local_parameters=[LocalParameter(id="k", label="", value="1")],
    )
    assert out == "0.0"


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
