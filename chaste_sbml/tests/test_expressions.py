"""Tests for SBML expression/AST -> C++ formula translation."""

import pytest
from libsbml import parseL3Formula

from chaste_sbml._config import VarType
from chaste_sbml._expressions import convert_infix_operator_to_function_syntax, formula_to_string
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
