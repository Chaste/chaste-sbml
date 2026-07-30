"""Loading, normalising and reading an SBML model with libsbml.

Reads an SBML file, flattens hierarchical models and runs the conversions the
generator depends on, returning a libsbml ``Model`` ready to be processed.
python-libsbml keeps the owning document alive through the returned model, so
callers need only hold the model.
Also provides small readers over individual SBML elements (compartment size,
function-definition arguments) used while building the internal representation.
"""

import logging
from typing import TYPE_CHECKING, Optional

from libsbml import (
    LIBSBML_OPERATION_SUCCESS,
    UNIT_KIND_SECOND,
    ConversionProperties,
    SBMLReader,
    formulaToString,
)

from ._config import TimeUnit

if TYPE_CHECKING:
    from libsbml import Compartment, FunctionDefinition, Model, UnitDefinition

logger = logging.getLogger(__name__)

# SBML levels the generator supports.
MIN_SBML_LEVEL = 2
MAX_SBML_LEVEL = 3


def load_sbml_model(sbml_file: str) -> tuple["Model", Optional[TimeUnit], int]:
    """Read, comp-flatten and convert an SBML file into a ready-to-process libsbml Model.

    :param sbml_file: Path to the SBML file.
    :return: A tuple of the loaded libsbml ``Model``, its declared time unit (or ``None`` if the
        model declares no determinable time unit), and its SBML level.
    :raises ValueError: if the file has read errors, or flattening/conversion fails.
    """
    reader = SBMLReader()
    doc = reader.readSBMLFromFile(sbml_file)

    if doc.getNumErrors() > 0:
        doc.printErrors()
        raise ValueError(f"Errors found while reading SBML file: {sbml_file}")

    # Flatten hierarchical (comp package) models, composing submodels into
    # a flattened model with submodelId__element names.
    if doc.getPlugin("comp") is not None:
        flatten_props = ConversionProperties()
        flatten_props.addOption("flatten comp", True, "flatten comp")
        flatten_props.addOption("leave_ports", False)
        if doc.convert(flatten_props) != LIBSBML_OPERATION_SUCCESS:  # pragma: no cover - libsbml failure path
            doc.printErrors()
            raise ValueError("Errors during comp flattening")

    # Detect the declared time unit before the conversions below: 'removeUnusedUnits' can strip an
    # unreferenced <unitDefinition id="time">, which would hide it from the builder.
    sbml_level = doc.getLevel()
    _warn_unsupported_level(sbml_level)
    declared_time_unit = detect_time_unit(doc.getModel())

    # Run required conversions
    config = ConversionProperties()
    # Sort assignment rules in order of dependence.
    config.addOption("sortRules")
    config.addOption("removeUnusedUnits")
    # Convert initial assignments to initial values where possible
    config.addOption("expandInitialAssignments")
    # config.addOption('replaceReactions')
    # config.addOption('expandFunctionDefinitions')

    status = doc.convert(config)
    if status != LIBSBML_OPERATION_SUCCESS:  # pragma: no cover - libsbml failure path
        doc.printErrors()
        raise ValueError("Errors during conversion")

    return doc.getModel(), declared_time_unit, sbml_level


def _warn_unsupported_level(sbml_level: int) -> None:
    """Warn if the SBML level is outside the supported range (Levels 2 and 3).

    Level 1 predates much of the structure the generator relies on, and no level above 3 exists yet;
    either way the generated code may be incorrect, so warn rather than fail.

    :param sbml_level: The model's SBML level.
    """
    if sbml_level < MIN_SBML_LEVEL or sbml_level > MAX_SBML_LEVEL:
        logger.warning(
            "SBML Level %d is not supported; only Levels %d and %d are supported. "
            "Generated code may be incorrect.",
            sbml_level,
            MIN_SBML_LEVEL,
            MAX_SBML_LEVEL,
        )


def detect_time_unit(model: "Model") -> Optional[TimeUnit]:
    """Determine the time unit a model *declares*, or ``None`` if it declares none we can resolve.

    Precedence: the SBML Level 3 ``model.timeUnits`` attribute (a base-unit keyword such as
    ``second``, or a ``UnitDefinition`` id), then a ``<unitDefinition id="time">`` (the Level 2
    convention). The unit is resolved structurally from the ``<unit>`` list (see ``_seconds_factor``),
    never from the human-readable unit name.

    :param model: The libsbml model.
    :return: The declared :class:`TimeUnit`, or ``None`` if undeclared/undeterminable.
    """
    unit_id = model.getTimeUnits() if model.isSetTimeUnits() else ""

    # An L3 timeUnits attribute may name a base unit directly rather than a UnitDefinition.
    builtin = _builtin_time_unit(unit_id)
    if builtin is not None:
        return builtin

    unit_def = model.getUnitDefinition(unit_id) if unit_id else None
    if unit_def is None:
        unit_def = model.getUnitDefinition("time")
    if unit_def is None:
        return None

    factor = _seconds_factor(unit_def)
    if factor is None:
        return None
    return TimeUnit.from_seconds_factor(factor)


def _builtin_time_unit(unit_id: str) -> Optional[TimeUnit]:
    """Map an SBML base-unit keyword used as a time unit to a TimeUnit.

    :param unit_id: A unit id, which may be a base-unit keyword (e.g. ``second``).
    :return: The matching :class:`TimeUnit`, or ``None`` if it is not a recognised base-unit keyword.
    """
    if unit_id == "second":
        return TimeUnit.SECOND
    if unit_id == "dimensionless":
        return TimeUnit.NONE
    return None


def _seconds_factor(unit_def: "UnitDefinition") -> Optional[float]:
    """Return the size in seconds of a UnitDefinition that is a scalar multiple of ``second^1``.

    :param unit_def: The unit definition.
    :return: ``multiplier * 10**scale`` seconds, or ``None`` if the definition is not exactly one
        second-kind unit with exponent 1 (composite units, non-second kinds and exponents other than
        1 are treated as undeterminable).
    """
    if unit_def.getNumUnits() != 1:
        return None
    unit = unit_def.getUnit(0)
    if unit.getKind() != UNIT_KIND_SECOND:
        return None
    if unit.getExponentAsDouble() != 1.0:
        return None
    # seconds = multiplier * 10^scale (offset is 0 for time units).
    return unit.getMultiplier() * (10.0**unit.getScale())


def get_compartment_size(compartment: "Compartment") -> float:
    """Get a compartment size.

    :return: The compartment size.
    """
    if compartment.isSetSize():
        return compartment.getSize()
    return 1.0


def get_function_definition_arguments(fn_def: "FunctionDefinition") -> list[str]:
    """Get the list of arguments in a given function definition.

    :param fn_def: The function definition
    :return: List of arguments in the function definition
    """
    n = fn_def.getNumArguments()
    return [formulaToString(fn_def.getArgument(i)) for i in range(n)]
