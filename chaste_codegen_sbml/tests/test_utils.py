import logging

from libsbml import parseFormula

from chaste_codegen_sbml._utils import varname_camelcase, varname_sanitize

logger = logging.getLogger(__name__)


def donttest_formula_to_cpp():
    """
    Test formula conversion
    """

    test_cases = [
        # basic cases
        ("", ""),
        ("x", "x"),
        # constants and math symbols
        ("true", "true"),
        ("false", "false"),
        ("time", "time"),  # TODO: Fix this
        ("avogadro", ""),  # TODO: Fix this
        ("pi", "M_PI"),
        ("inf", "INFINITY"),
        ("infinity", "INFINITY"),
        ("nan", "NAN"),
        ("notanumber", "NAN"),
        # arithmetic operators
        ("x / y", "x / y"),
        ("x - y", "x - y"),
        ("x + y", "x + y"),
        ("x * y", "x * y"),
        ("x % y", "x % y"),  # TODO: convert to fmod(x, y) instead?
        # arithmentic functions
        ("divide(x, y)", "x / y"),
        ("minus(x, y)", "x - y"),
        ("plus(x, y, z)", "x + y + z"),
        ("times(x, y, z)", "x * y * z"),
        # power operator
        ("x ^ y", "pow(x, y)"),
        # power functions
        ("pow(x, y)", "pow(x, y)"),
        ("power(x, y)", "pow(x, y)"),
        ("root(b, x)", "pow(x, 1.0 / b)"),
        ("sqr(x)", "pow(x, 2)"),
        ("sqrt(x)", "sqrt(x)"),
        # logical operators
        ("!x", "!x"),
        ("x && y", "x && y"),
        ("x || y", "x || y"),
        # logical functions
        ("not(x)", "!(x)"),
        ("and()", "true"),
        ("and(x, y, z)", "((x) && (y) && (z))"),
        ("or()", "false"),
        ("or(x, y, z)", "((x) || (y) || (z))"),
        ("xor()", "false"),
        # TODO: check if this is correct
        ("xor(x, y, z)", "(((x) != 0) ^ ((y) != 0) ^ ((z) != 0))"),
        # comparison operators
        ("x == y", "x == y"),
        ("x >= y", "x >= y"),
        ("x > y", "x > y"),
        ("x <= y", "x <= y"),
        ("x < y", "x < y"),
        ("x != y", "x != y"),
        # comparison functions
        ("eq(x, y, z)", "(((x) == (y)) && ((y) == (z)))"),
        ("geq(x, y, z)", "(((x) >= (y)) && ((y) >= (z)))"),
        ("gt(x, y, z)", "(((x) > (y)) && ((y) > (z)))"),
        ("leq(x, y, z)", "(((x) <= (y)) && ((y) <= (z)))"),
        ("lt(x, y, z)", "(((x) < (y)) && ((y) < (z)))"),
        ("neq(x, y)", "((x) != (y))"),
        # trigonometric functions - no change
        ("acos(x)", "acos(x)"),
        ("acosh(x)", "acosh(x)"),
        ("acsc(x)", "acsc(x)"),
        ("acsch(x)", "acsch(x)"),
        ("asec(x)", "asec(x)"),
        ("asech(x)", "asech(x)"),
        ("asin(x)", "asin(x)"),
        ("asinh(x)", "asinh(x)"),
        ("atan(x)", "atan(x)"),
        ("atanh(x)", "atanh(x)"),
        ("cos(x)", "cos(x)"),
        ("cosh(x)", "cosh(x)"),
        ("sin(x)", "sin(x)"),
        ("sinh(x)", "sinh(x)"),
        ("tan(x)", "tan(x)"),
        ("tanh(x)", "tanh(x)"),
        # trigonometric functions - requires name change
        ("arccos(x)", "acos(x)"),
        ("arcosh(x)", "acosh(x)"),
        ("arccsc(x)", "acsc(x)"),
        ("arccsch(x)", "acsch(x)"),
        ("arcsec(x)", "asec(x)"),
        ("arcsech(x)", "asech(x)"),
        ("arcsin(x)", "asin(x)"),
        ("arcsinh(x)", "asinh(x)"),
        ("arctan(x)", "atan(x)"),
        ("arctanh(x)", "atanh(x)"),
        # trigonometric functions - requires conversion
        ("acot(x)", "atan(1.0 / x)"),
        ("arccot(x)", "atan(1.0 / x)"),
        ("acoth(x)", "atanh(1.0 / x)"),
        ("arccoth(x)", "atanh(1.0 / x)"),
        ("cot(x)", "(1.0 / tan(x))"),
        ("coth(x)", "(1.0 / tanh(x))"),
        ("csc(x)", "(1.0 / sin(x))"),
        ("csch(x)", "(1.0 / sinh(x))"),
        ("sec(x)", "(1.0 / cos(x))"),
        ("sech(x)", "(1.0 / cosh(x))"),
        # other functions - no change
        ("ceil(x)", "ceil(x)"),
        ("exp(x)", "exp(x)"),
        ("floor(x)", "floor(x)"),
        ("log10(x)", "log10(x)"),
        # other functions - requires name change
        ("abs(x)", "fabs(x)"),
        ("ceiling(x)", "ceil(x)"),
        ("log(x)", "log10(x)"),
        ("ln(x)", "log(x)"),
        ("max(x, y)", "fmax(x, y)"),
        ("min(x, y)", "fmin(x, y)"),
        # other functions - requires conversion
        ("log(x, y)", "(log(y) / log(x))"),
        ("piecewise(x1, y1, x2, y2, z)", "(y1 ? x1 : (y2 ? x2 : z))"),
        # other functions - unimplemented
        ("delay(x, y)", "delay(x, y)"),  # TODO: Fix this
        ("factorial(n)", "factorial(n)"),  # TODO: Fix this
    ]

    for formula_string, expected_cpp in test_cases:
        formula_tree = parseFormula(formula_string)
        assert formula_to_cpp(formula_tree) == expected_cpp


def test_varname_camelcase():
    """
    Test camel case variable name creation.
    """

    test_cases = [
        ("", ""),
        (" ", ""),
        ("foo", "foo"),
        ("foo_bar", "fooBar"),
        ("foo bar", "fooBar"),
        ("foo.bar", "fooBar"),
        ("abstract_foo_bar", "abstractFooBar"),
        ("foo1bar", "foo1Bar"),
        ("foo100bar", "foo100Bar"),
    ]

    for test_input, expected_output in test_cases:
        assert varname_camelcase(test_input) == expected_output


def test_varname_sanitize():
    """
    Test variable name sanitization.
    """

    test_cases = [
        ("", ""),
        (" ", ""),
        ("foo", "foo"),
        ("foo_bar", "foo_bar"),
        ("foo___bar", "foo___bar"),
        ("foo bar", "foo_bar"),
        ("foo_ _bar", "foo___bar"),
        ("foo_ . _bar", "foo___bar"),
        ("foo1bar", "foo1bar"),
        ("1foo", "_1foo"),
        ("foo100bar", "foo100bar"),
    ]

    for test_input, expected_output in test_cases:
        assert varname_sanitize(test_input) == expected_output
