"""Tests for SBML expression/AST -> C++ formula translation."""

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
    convert_infix_operator_to_function_syntax,
    formula_to_string,
    replace_avogadro_csymbol,
    search_ast_type,
    strip_ast_units,
    substitute_ast_names,
)
from chaste_sbml._records import StateVariable


def test_convert_infix_operator_to_function_syntax_power():
    """Converts infix power expressions to pow calls."""
    converted = convert_infix_operator_to_function_syntax(
        formula="(a + b) ^ (1.0 / 2.0)",
        operator="^",
        function_name="pow",
    )

    assert converted == "pow((a + b), (1.0 / 2.0))"


def test_convert_infix_operator_to_function_syntax_custom_pair():
    """Supports arbitrary operator/function conversion pairs."""
    converted = convert_infix_operator_to_function_syntax(
        formula="alpha @@ beta", operator="@@", function_name="combine"
    )

    assert converted == "combine(alpha, beta)"


def test_convert_infix_operator_to_function_syntax_nested_and_repeated():
    """Converts nested and repeated operator expressions in one formula."""
    converted = convert_infix_operator_to_function_syntax(formula="(x ^ y) ^ z", operator="^", function_name="pow")

    assert converted == "pow((pow(x, y)), z)"


@pytest.mark.parametrize(
    ("operator", "function_name", "error_message"),
    [
        ("", "pow", "operator must be a non-empty string"),
        ("^", "", "function_name must be a non-empty string"),
    ],
)
def test_convert_infix_operator_to_function_syntax_invalid_args(operator: str, function_name: str, error_message: str):
    """Validates required arguments for operator/function conversion."""
    with pytest.raises(ValueError, match=error_message):
        convert_infix_operator_to_function_syntax(formula="a ^ b", operator=operator, function_name=function_name)


def test_formula_to_string_functions_constants_and_power():
    """Renders SBML functions/constants and infix power as C++."""
    math = parseL3Formula("k * pow(S, 2) + pi")
    out = formula_to_string(
        math, variable_types={"k": VarType.PARAMETER, "S": VarType.STATE_VARIABLE}, state_variables=[]
    )
    assert out == "k * std::pow(S, 2.0) + M_PI"


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


def test_replace_avogadro_csymbol_renames_nested_node():
    """An avogadro csymbol anywhere in the tree becomes a plain name with the placeholder id."""
    one = ASTNode(AST_REAL)
    one.setValue(1.0)
    avogadro = ASTNode(AST_NAME_AVOGADRO)
    expr = ASTNode(AST_PLUS)
    expr.addChild(one)
    expr.addChild(avogadro)

    replace_avogadro_csymbol(expr, "AVO")

    renamed = expr.getChild(1)
    assert renamed.getType() == AST_NAME
    assert renamed.getName() == "AVO"


def test_replace_avogadro_csymbol_leaves_other_names():
    """A non-avogadro name node is left untouched."""
    name = ASTNode(AST_NAME)
    name.setName("x")
    replace_avogadro_csymbol(name, "AVO")
    assert name.getName() == "x"


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
