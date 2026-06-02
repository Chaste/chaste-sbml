"""Unit tests for ChasteSbmlModel internals."""

import pytest

from chaste_codegen_sbml.chaste_sbml_model import ChasteSbmlModel


def _model_without_init() -> ChasteSbmlModel:
    """Create an instance without running __init__ for pure string helper tests."""
    return object.__new__(ChasteSbmlModel)


def test_convert_infix_operator_to_function_syntax_power():
    """Converts infix power expressions to pow calls."""
    model = _model_without_init()

    converted = model._convert_infix_operator_to_function_syntax(
        formula="(a + b) ^ (1.0 / 2.0)",
        operator="^",
        function_name="pow",
    )

    assert converted == "pow((a + b), (1.0 / 2.0))"


def test_convert_infix_operator_to_function_syntax_custom_pair():
    """Supports arbitrary operator/function conversion pairs."""
    model = _model_without_init()

    converted = model._convert_infix_operator_to_function_syntax(
        formula="alpha @@ beta", operator="@@", function_name="combine"
    )

    assert converted == "combine(alpha, beta)"


def test_convert_infix_operator_to_function_syntax_nested_and_repeated():
    """Converts nested and repeated operator expressions in one formula."""
    model = _model_without_init()

    converted = model._convert_infix_operator_to_function_syntax(
        formula="(x ^ y) ^ z", operator="^", function_name="pow"
    )

    assert converted == "pow((pow(x, y)), z)"


@pytest.mark.parametrize(
    ("operator", "function_name", "error_message"),
    [
        ("", "pow", "operator must be a non-empty string"),
        ("^", "", "function_name must be a non-empty string"),
    ],
)
def test_convert_infix_operator_to_function_syntax_invalid_args(
    operator: str, function_name: str, error_message: str
):
    """Validates required arguments for operator/function conversion."""
    model = _model_without_init()

    with pytest.raises(ValueError, match=error_message):
        model._convert_infix_operator_to_function_syntax(
            formula="a ^ b", operator=operator, function_name=function_name
        )
