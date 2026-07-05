"""Tests for C++ identifier conflict detection and naming."""

import libsbml
import pytest

from chaste_sbml._names import (
    CHASTE_RESERVED_NAMES,
    CPP_KEYWORDS,
    NameManager,
    find_name_conflicts,
    generate_header_guard,
    is_valid_cpp_identifier,
    resolve_cpp_name,
    to_camel_case,
    to_cpp_name,
    unique_name,
)


@pytest.mark.parametrize("name", ["C", "amt__X", "d_C_dt", "_leading", "x2", "chaste__avogadro"])
def test_valid_cpp_identifiers(name):
    """Accepts syntactically valid C++ identifiers."""
    assert is_valid_cpp_identifier(name)


@pytest.mark.parametrize("name", ["2x", "a-b", "has space", "", "a.b", "π"])
def test_invalid_cpp_identifiers(name):
    """Rejects names that are not valid C++ identifiers."""
    assert not is_valid_cpp_identifier(name)


def test_no_conflicts_returns_empty():
    """A set of distinct, valid, non-reserved names has no conflicts."""
    assert find_name_conflicts([("C", "state variable"), ("d_C_dt", "derivative"), ("VM1", "parameter")]) == []


def test_detects_duplicate():
    """A name emitted by two entities is flagged, naming both kinds."""
    conflicts = find_name_conflicts([("amt__X", "amount/concentration conversion"), ("amt__X", "state variable")])
    assert conflicts == ["'amt__X' is emitted twice: as amount/concentration conversion and as state variable"]


@pytest.mark.parametrize("keyword", ["double", "int", "class", "new", "template", "this"])
def test_detects_cpp_keyword(keyword):
    """C++ keywords used as identifiers are flagged."""
    assert find_name_conflicts([(keyword, "parameter")]) == [f"parameter '{keyword}' is a C++ keyword"]


@pytest.mark.parametrize("reserved", ["mParameters", "GetParameter", "mStateVariables", "time"])
def test_detects_reserved_chaste_name(reserved):
    """Names inherited from the Chaste base classes are flagged."""
    assert find_name_conflicts([(reserved, "state variable")]) == [
        f"state variable '{reserved}' clashes with a reserved Chaste name"
    ]


def test_detects_invalid_identifier():
    """Non-identifier names are flagged."""
    assert find_name_conflicts([("2x", "parameter")]) == ["parameter '2x' is not a valid C++ identifier"]


def test_reports_multiple_sorted():
    """All conflicts in one batch are reported, sorted for stable output."""
    conflicts = find_name_conflicts([("int", "parameter"), ("X", "state variable"), ("X", "reaction")])
    assert conflicts == [
        "'X' is emitted twice: as state variable and as reaction",
        "parameter 'int' is a C++ keyword",
    ]


def test_keyword_and_reserved_sets_disjoint_from_ordinary_names():
    """Sanity: representative model ids are neither keywords nor reserved."""
    for name in ("C", "M", "X", "VM1", "cell", "reaction1", "amt__C"):
        assert name not in CPP_KEYWORDS
        assert name not in CHASTE_RESERVED_NAMES


def test_unique_name_returns_base_when_free():
    """A free base name is returned unchanged, so names stay clean."""
    assert unique_name("amt__X", set()) == "amt__X"
    assert unique_name("d_C_dt", {"C", "M", "X"}) == "d_C_dt"


def test_unique_name_appends_suffix_on_collision():
    """A taken base name is escaped with the smallest free numeric suffix."""
    assert unique_name("amt__X", {"amt__X"}) == "amt__X_2"
    assert unique_name("amt__X", {"amt__X", "amt__X_2"}) == "amt__X_3"


def test_unique_name_skips_taken_suffixes():
    """Already-taken suffixed names are skipped."""
    assert unique_name("v", {"v", "v_2", "v_3"}) == "v_4"


def test_resolve_cpp_name_leaves_safe_id_unchanged():
    """A real id that is neither keyword nor reserved and is free is returned unchanged."""
    assert resolve_cpp_name("cell", {"C", "M"}) == "cell"


@pytest.mark.parametrize(("keyword", "expected"), [("default", "default_"), ("int", "int_"), ("new", "new_")])
def test_resolve_cpp_name_escapes_keywords(keyword, expected):
    """A keyword id is escaped by suffixing an underscore."""
    assert resolve_cpp_name(keyword, set()) == expected


def test_resolve_cpp_name_escapes_reserved_chaste_name():
    """A reserved Chaste member name is escaped."""
    assert resolve_cpp_name("mParameters", set()) == "mParameters_"


def test_resolve_cpp_name_escape_then_uniquify():
    """If the escaped name is already taken, a numeric suffix is added."""
    assert resolve_cpp_name("default", {"default_"}) == "default__2"


def test_name_manager_reserve_keeps_clean_name_when_free():
    """A synthetic name that does not collide is returned unchanged and then reserved."""
    manager = NameManager(None)
    manager._taken = {"X", "cell"}
    assert manager.reserve("amt__X") == "amt__X"
    assert "amt__X" in manager._taken


def test_name_manager_reserve_escapes_collision_with_real_id():
    """A synthetic name equal to a real id is escaped, and repeats take further suffixes."""
    manager = NameManager(None)
    manager._taken = {"amt__X"}  # a real species literally named amt__X
    assert manager.reserve("amt__X") == "amt__X_2"
    assert manager.reserve("amt__X") == "amt__X_3"


def test_name_manager_resolve_renames_keyword_compartment():
    """A compartment whose id is a C++ keyword is renamed, and its references updated."""
    doc = libsbml.SBMLDocument(3, 2)
    sbml_model = doc.createModel()
    compartment = sbml_model.createCompartment()
    compartment.setId("default")  # 'default' is a C++ keyword
    compartment.setConstant(True)
    compartment.setSize(1.0)
    species = sbml_model.createSpecies()
    species.setId("S")
    species.setCompartment("default")
    species.setConstant(False)
    species.setBoundaryCondition(False)
    species.setHasOnlySubstanceUnits(False)
    species.setInitialAmount(1.0)

    NameManager(sbml_model).resolve_real_id_conflicts()

    assert sbml_model.getElementBySId("default") is None
    assert sbml_model.getElementBySId("default_") is not None
    # The species' compartment reference was rewritten too.
    assert sbml_model.getSpecies("S").getCompartment() == "default_"


def test_name_manager_resolve_leaves_safe_ids_untouched():
    """A model with only safe ids is not modified."""
    doc = libsbml.SBMLDocument(3, 2)
    sbml_model = doc.createModel()
    compartment = sbml_model.createCompartment()
    compartment.setId("cell")
    compartment.setConstant(True)

    NameManager(sbml_model).resolve_real_id_conflicts()

    assert sbml_model.getElementBySId("cell") is not None


def test_name_manager_sbml_name_recovers_original_id():
    """After a keyword id is renamed, sbml_name maps the C++ id back to the original SBML id.

    A variable must be reported to Chaste under its real SBML id even though the emitted C++
    identifier had to be escaped, so it can be looked up by that id.
    """
    doc = libsbml.SBMLDocument(3, 2)
    sbml_model = doc.createModel()
    parameter = sbml_model.createParameter()
    parameter.setId("true")  # 'true' is a C++ keyword, so it is renamed to 'true_'
    parameter.setConstant(True)
    parameter.setValue(1.0)

    manager = NameManager(sbml_model)
    manager.resolve_real_id_conflicts()

    assert sbml_model.getParameter("true_") is not None
    assert manager.sbml_name("true_") == "true"  # renamed id maps back to the real SBML id
    assert manager.sbml_name("k") == "k"  # an un-renamed id is returned unchanged


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
