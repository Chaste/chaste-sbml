"""Tests for SBML→Chaste time-unit detection, resolution and scaling."""

import logging

import libsbml
import pytest

from chaste_sbml import ChasteSbmlModel
from chaste_sbml._config import ROOT_DIR, ModelType, TimeUnit
from chaste_sbml._sbml_reader import _builtin_time_unit, _seconds_factor, detect_time_unit, load_sbml_model

GOLDBETER = ROOT_DIR / "SbmlRefModels" / "src" / "reference" / "Goldbeter1991" / "Goldbeter1991.xml"
TYSON = ROOT_DIR / "SbmlRefModels" / "src" / "reference" / "TysonNovak2001" / "TysonNovak2001.xml"


# --- TimeUnit enum -------------------------------------------------------------------------------


def test_multiplier_literals_are_plain_cpp_doubles():
    """Every unit's multiplier renders as an exact C++ double literal (no scientific notation)."""
    assert TimeUnit.SECOND.multiplier_literal == "3600.0"
    assert TimeUnit.MINUTE.multiplier_literal == "60.0"
    assert TimeUnit.HOUR.multiplier_literal == "1.0"
    assert TimeUnit.MILLISECOND.multiplier_literal == "3600000.0"
    assert TimeUnit.NONE.multiplier_literal == "1.0"


def test_from_cli_maps_every_token():
    """Each --timescale token maps to the matching unit."""
    assert TimeUnit.from_cli("ms") is TimeUnit.MILLISECOND
    assert TimeUnit.from_cli("s") is TimeUnit.SECOND
    assert TimeUnit.from_cli("m") is TimeUnit.MINUTE
    assert TimeUnit.from_cli("h") is TimeUnit.HOUR


def test_from_seconds_factor_matches_and_rejects():
    """A recognised seconds factor maps to a unit; an unknown one returns None."""
    assert TimeUnit.from_seconds_factor(60.0) is TimeUnit.MINUTE
    assert TimeUnit.from_seconds_factor(3600.0) is TimeUnit.HOUR
    assert TimeUnit.from_seconds_factor(1.0) is TimeUnit.SECOND
    assert TimeUnit.from_seconds_factor(1.0e-3) is TimeUnit.MILLISECOND
    assert TimeUnit.from_seconds_factor(7.0) is None


# --- detection helpers ---------------------------------------------------------------------------


def _l2_time_unit_def(multiplier=1.0, scale=0, exponent=1, kind=libsbml.UNIT_KIND_SECOND, extra_units=0):
    """Build a Level 2 model with a <unitDefinition id="time"> and return the model."""
    doc = libsbml.SBMLDocument(2, 4)
    model = doc.createModel()
    unit_def = model.createUnitDefinition()
    unit_def.setId("time")
    unit = unit_def.createUnit()
    unit.setKind(kind)
    unit.setExponent(exponent)
    unit.setScale(scale)
    unit.setMultiplier(multiplier)
    for _ in range(extra_units):
        pad = unit_def.createUnit()
        pad.setKind(libsbml.UNIT_KIND_DIMENSIONLESS)
        pad.setExponent(1)
    return doc, model


@pytest.mark.parametrize(
    ("multiplier", "scale", "expected"),
    [
        (60.0, 0, TimeUnit.MINUTE),
        (3600.0, 0, TimeUnit.HOUR),
        (1.0, 0, TimeUnit.SECOND),
        (1.0e-3, 0, TimeUnit.MILLISECOND),
        (1.0, -3, TimeUnit.MILLISECOND),  # scale contributes 10^-3
    ],
)
def test_detect_l2_unit_definition(multiplier, scale, expected):
    """A Level 2 <unitDefinition id="time"> resolves structurally from its <unit>."""
    _doc, model = _l2_time_unit_def(multiplier=multiplier, scale=scale)
    assert detect_time_unit(model) is expected


def test_detect_returns_none_for_unrecognised_factor():
    """A second-based unit whose factor matches no known unit is undeterminable."""
    _doc, model = _l2_time_unit_def(multiplier=7.0)
    assert detect_time_unit(model) is None


def test_detect_returns_none_for_composite_unit():
    """A composite (multi-<unit>) time definition is not a scalar multiple of seconds."""
    _doc, model = _l2_time_unit_def(multiplier=60.0, extra_units=1)
    assert detect_time_unit(model) is None


def test_detect_returns_none_for_non_second_kind():
    """A time definition whose kind is not second is undeterminable."""
    _doc, model = _l2_time_unit_def(kind=libsbml.UNIT_KIND_MOLE)
    assert detect_time_unit(model) is None


def test_detect_returns_none_for_exponent_not_one():
    """A second unit with exponent != 1 is undeterminable."""
    _doc, model = _l2_time_unit_def(exponent=2)
    assert detect_time_unit(model) is None


def test_detect_returns_none_when_no_units_declared():
    """A model with no time unit declares nothing."""
    doc = libsbml.SBMLDocument(2, 4)
    model = doc.createModel()
    assert detect_time_unit(model) is None


def test_detect_l3_time_units_builtin_second():
    """An L3 timeUnits attribute naming the base unit 'second' resolves to seconds."""
    doc = libsbml.SBMLDocument(3, 2)
    model = doc.createModel()
    model.setTimeUnits("second")
    assert detect_time_unit(model) is TimeUnit.SECOND


def test_detect_l3_time_units_dimensionless_is_no_conversion():
    """An L3 timeUnits of 'dimensionless' means no conversion."""
    doc = libsbml.SBMLDocument(3, 2)
    model = doc.createModel()
    model.setTimeUnits("dimensionless")
    assert detect_time_unit(model) is TimeUnit.NONE


def test_detect_l3_time_units_referencing_unit_definition():
    """An L3 timeUnits attribute pointing at a UnitDefinition resolves from that definition."""
    doc = libsbml.SBMLDocument(3, 2)
    model = doc.createModel()
    unit_def = model.createUnitDefinition()
    unit_def.setId("myTime")
    unit = unit_def.createUnit()
    unit.setKind(libsbml.UNIT_KIND_SECOND)
    unit.setExponent(1)
    unit.setScale(0)
    unit.setMultiplier(60.0)
    model.setTimeUnits("myTime")
    assert detect_time_unit(model) is TimeUnit.MINUTE


def test_detect_l3_unset_time_units_is_none():
    """An L3 model with unset timeUnits and no id='time' declares nothing."""
    doc = libsbml.SBMLDocument(3, 2)
    model = doc.createModel()
    assert detect_time_unit(model) is None


def test_builtin_time_unit_unknown_keyword():
    """A non base-unit keyword is not a builtin time unit."""
    assert _builtin_time_unit("myTime") is None


def test_seconds_factor_applies_scale():
    """_seconds_factor combines multiplier and scale (seconds = multiplier * 10**scale)."""
    _doc, model = _l2_time_unit_def(multiplier=6.0, scale=1)  # 6 * 10^1 = 60 s
    assert _seconds_factor(model.getUnitDefinition("time")) == 60.0


# --- detection through load_sbml_model (survives removeUnusedUnits) --------------------------------


def test_detection_through_load_reads_unreferenced_time_unit(tmp_path):
    """An unreferenced <unitDefinition id='time'> is detected via load_sbml_model.

    Detection runs before the removeUnusedUnits conversion, so an otherwise-unreferenced time unit
    is captured regardless of whether that pass would strip it.
    """
    sbml = """<?xml version="1.0" encoding="UTF-8"?>
<sbml xmlns="http://www.sbml.org/sbml/level2/version4" level="2" version="4">
  <model id="m">
    <listOfUnitDefinitions>
      <unitDefinition id="time" name="minutes">
        <listOfUnits><unit kind="second" multiplier="60"/></listOfUnits>
      </unitDefinition>
    </listOfUnitDefinitions>
    <listOfCompartments><compartment id="c" size="1"/></listOfCompartments>
    <listOfSpecies>
      <species id="S" compartment="c" initialConcentration="1"/>
    </listOfSpecies>
  </model>
</sbml>
"""
    path = tmp_path / "unref_time.xml"
    path.write_text(sbml)

    _model, declared, level = load_sbml_model(str(path))
    assert level == 2
    assert declared is TimeUnit.MINUTE


# --- resolution rule -----------------------------------------------------------------------------


@pytest.fixture(scope="module")
def any_model():
    """A constructed model, used only to reach the pure _resolve_time_unit method."""
    return ChasteSbmlModel(str(GOLDBETER))


def test_resolve_override_wins_and_warns_on_mismatch(any_model, caplog):
    """An explicit override overrides a declared unit and logs a mismatch warning."""
    with caplog.at_level(logging.WARNING):
        resolved = any_model._resolve_time_unit(TimeUnit.MINUTE, TimeUnit.HOUR, 2)
    assert resolved is TimeUnit.MINUTE
    assert "differs" in caplog.text


def test_resolve_override_without_declared_does_not_warn(any_model, caplog):
    """An override with no declared unit is used silently."""
    with caplog.at_level(logging.WARNING):
        resolved = any_model._resolve_time_unit(TimeUnit.SECOND, None, 2)
    assert resolved is TimeUnit.SECOND
    assert "differs" not in caplog.text


def test_resolve_override_matching_declared_does_not_warn(any_model, caplog):
    """An override equal to the declared unit does not warn."""
    with caplog.at_level(logging.WARNING):
        resolved = any_model._resolve_time_unit(TimeUnit.MINUTE, TimeUnit.MINUTE, 2)
    assert resolved is TimeUnit.MINUTE
    assert "differs" not in caplog.text


def test_resolve_uses_declared_when_no_override(any_model):
    """A declared unit is used when there is no override."""
    assert any_model._resolve_time_unit(None, TimeUnit.MINUTE, 2) is TimeUnit.MINUTE


def test_resolve_defaults_level2_to_seconds(any_model):
    """An undeclared Level 2 model defaults to seconds (SBML L2 predefines time = second)."""
    assert any_model._resolve_time_unit(None, None, 2) is TimeUnit.SECOND


def test_resolve_defaults_level3_to_no_conversion(any_model):
    """An undeclared Level 3 model leaves time undefined (no conversion)."""
    assert any_model._resolve_time_unit(None, None, 3) is TimeUnit.NONE


# --- end-to-end resolution + rendering ------------------------------------------------------------


def test_model_detects_minutes_and_scales(tmp_path):
    """A minutes model resolves to minutes and emits the ×60 scaling in the ODE system."""
    model = ChasteSbmlModel(str(TYSON), model_type=ModelType.GENERIC)
    assert model._time_unit is TimeUnit.MINUTE
    model.write(str(tmp_path))
    cpp = (tmp_path / "TysonNovak2001SbmlOdeSystem.cpp").read_text()
    assert "time *= 60.0;" in cpp
    assert "rDY[i] = 60.0 * derivatives[i];" in cpp


def test_timescale_override_forces_conversion(tmp_path):
    """An explicit --timescale-style override converts an otherwise-undeclared model."""
    # Goldbeter1991 declares no unit; forcing minutes must emit ×60 scaling.
    model = ChasteSbmlModel(str(GOLDBETER), time_unit=TimeUnit.MINUTE)
    assert model._time_unit is TimeUnit.MINUTE
    model.write(str(tmp_path))
    cpp = (tmp_path / "Goldbeter1991SbmlOdeSystem.cpp").read_text()
    assert "time *= 60.0;" in cpp


def test_timescale_none_override_suppresses_conversion(tmp_path):
    """Forcing TimeUnit.NONE leaves an otherwise-seconds model unscaled."""
    model = ChasteSbmlModel(str(GOLDBETER), time_unit=TimeUnit.NONE)
    assert model._time_unit is TimeUnit.NONE
    model.write(str(tmp_path))
    cpp = (tmp_path / "Goldbeter1991SbmlOdeSystem.cpp").read_text()
    assert "time *=" not in cpp
