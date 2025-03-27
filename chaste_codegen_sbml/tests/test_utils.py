import logging

from chaste_codegen_sbml._utils import varname_camelcase, varname_sanitize

logger = logging.getLogger(__name__)


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
