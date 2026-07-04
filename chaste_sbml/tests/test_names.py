"""Tests for C++ identifier conflict detection (issue #35, phase A)."""

import pytest

from chaste_sbml._names import (
    CHASTE_RESERVED_NAMES,
    CPP_KEYWORDS,
    find_name_conflicts,
    is_valid_cpp_identifier,
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
