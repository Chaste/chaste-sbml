"""Tests for loading and reading SBML with libsbml."""

import logging

import libsbml
import pytest
from libsbml import parseL3Formula

from chaste_sbml._config import ROOT_DIR
from chaste_sbml._sbml_reader import (
    _warn_unsupported_level,
    get_compartment_size,
    get_function_definition_arguments,
    load_sbml_model,
)

GOLDBETER = ROOT_DIR / "SbmlRefModels" / "src" / "reference" / "Goldbeter1991" / "Goldbeter1991.xml"


def test_load_sbml_model_reads_a_valid_file():
    """A valid SBML file loads into a libsbml Model, with its declared time unit and SBML level."""
    model, declared_time_unit, sbml_level = load_sbml_model(str(GOLDBETER))
    assert model.getNumSpecies() == 3
    # Goldbeter1991 is SBML Level 2 and declares no time unit.
    assert sbml_level == 2
    assert declared_time_unit is None


def test_load_sbml_model_raises_on_read_errors(tmp_path):
    """A malformed file raises ValueError rather than returning a broken model."""
    bad = tmp_path / "bad.xml"
    bad.write_text("<not well formed <<<")
    with pytest.raises(ValueError, match="Errors found while reading"):
        load_sbml_model(str(bad))


@pytest.mark.parametrize("level", [1, 4])
def test_warn_unsupported_level_warns(level, caplog):
    """SBML levels outside 2-3 (e.g. Level 1, or a future Level 4) trigger a warning."""
    with caplog.at_level(logging.WARNING):
        _warn_unsupported_level(level)
    assert f"SBML Level {level} is not supported" in caplog.text


@pytest.mark.parametrize("level", [2, 3])
def test_warn_unsupported_level_quiet_for_supported(level, caplog):
    """Supported SBML levels (2 and 3) produce no warning."""
    with caplog.at_level(logging.WARNING):
        _warn_unsupported_level(level)
    assert "not supported" not in caplog.text


def test_get_compartment_size_returns_set_size():
    """An explicitly set compartment size is returned."""
    doc = libsbml.SBMLDocument(3, 2)
    compartment = doc.createModel().createCompartment()
    compartment.setSize(2.5)
    assert get_compartment_size(compartment) == 2.5


def test_get_compartment_size_defaults_to_one():
    """A compartment with no size defaults to 1.0."""
    doc = libsbml.SBMLDocument(3, 2)
    compartment = doc.createModel().createCompartment()
    assert get_compartment_size(compartment) == 1.0


def test_get_function_definition_arguments_lists_bound_variables():
    """The bound variables of a lambda are returned in order, excluding the body."""
    doc = libsbml.SBMLDocument(3, 2)
    fn_def = doc.createModel().createFunctionDefinition()
    fn_def.setId("f")
    fn_def.setMath(parseL3Formula("lambda(x, y, x + y)"))
    assert get_function_definition_arguments(fn_def) == ["x", "y"]
