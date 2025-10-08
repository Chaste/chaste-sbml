"""Tests for utility functions."""

import logging

from chaste_codegen_sbml._utils import generate_header_guard, to_camel_case, to_cpp_name

logger = logging.getLogger(__name__)


def test_generate_header_guard():
    """Test header guard generation."""
    test_cases = [
        ("", ""),
        (" ", ""),
        ("foo.hpp", "FOO_HPP_"),
        ("foo_bar.hpp", "FOO_BAR_HPP_"),
        ("fooBar.hpp", "FOO_BAR_HPP_"),
        ("FooBar.hpp", "FOO_BAR_HPP_"),
        ("Foo100bar.hpp", "FOO_100BAR_HPP_"),
        ("Foo200Bar.hpp", "FOO_200_BAR_HPP_"),
        ("FooB300ar.hpp", "FOO_B300AR_HPP_"),
        ("400FooBar.hpp", "_400_FOO_BAR_HPP_"),
        ("TysonNovak2001.hpp", "TYSON_NOVAK_2001_HPP_"),
        ("TestSemantic00001L2V5Sbml.hpp", "TEST_SEMANTIC_00001_L2_V5_SBML_HPP_"),
    ]

    for test_input, expected_output in test_cases:
        assert generate_header_guard(test_input) == expected_output


def test_to_camel_case():
    """Test variable name case staggering."""
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
        ("TysonNovak2001", "TysonNovak2001"),
    ]

    for test_input, expected_output in test_cases:
        assert to_camel_case(test_input) == expected_output


def test_to_cpp_name():
    """Test variable name sanitization."""
    test_cases = [
        ("", ""),
        (" ", ""),
        ("foo", "foo"),
        ("foo_bar", "foo_bar"),
        ("foo___bar", "foo___bar"),
        ("foo bar", "foo_bar"),
        ("foo_ _bar", "foo___bar"),
        ("foo_ . _bar", "foo_____bar"),
        ("foo1bar", "foo1bar"),
        ("1foo", "_1foo"),
        ("foo100bar", "foo100bar"),
    ]

    for test_input, expected_output in test_cases:
        assert to_cpp_name(test_input) == expected_output
