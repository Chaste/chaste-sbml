"""Configuration settings and constants for code generation."""

import math
import pathlib
from enum import Enum
from typing import Optional

ROOT_DIR = pathlib.Path(__file__).parent.absolute()

NON_DIM_UNITS = "non-dim"

PREFIX_SEP = "__"
AMOUNT_PREFIX = "amt"
CONCENTRATION_PREFIX = "conc"
CHASTE_PREFIX = "chaste"
INITIAL_ASSIGNMENT_PREFIX = "ia"

DERIVATIVE_PREFIX = "d_"
DERIVATIVE_SUFFIX = "_dt"

# Base id for the synthetic zero-derivative state variable added to models with no ODEs, so the
# generated ODE system always has at least one variable for the solver to integrate.
PLACEHOLDER_STATE_ID = "placeholder"


class EquationType(Enum):
    """Enumeration of equation types for code generation."""

    ASSIGNMENT_RULE = 0
    CONVERSION = 1
    DERIVATIVE = 2
    INITIAL_ASSIGNMENT = 3
    INITIAL_VALUE = 4
    REACTION = 5
    UNKNOWN = 6


class EventType(Enum):
    """Enumeration of event types in SBML models."""

    CELL_DIVISION = 1
    UNKNOWN = 2


class VarType(Enum):
    """Enumeration of variable types in SBML models."""

    DERIVED_QUANTITY = 0
    FUNCTION = 1
    PARAMETER = 2
    REACTION = 3
    STATE_VARIABLE = 4
    UNKNOWN = 5


class DerivedQuantityKind(Enum):
    """Kind of derived quantity, which determines how it is declared and computed.

    NORMAL is declared as a member and computed from an equation; REACTION is a reaction flux
    (declared and computed by the reaction machinery, only exposed here as an output); CONVERSION
    is an amount/concentration conversion computed in ComputeDerivedQuantities.
    """

    NORMAL = 0
    REACTION = 1
    CONVERSION = 2


class ModelType(Enum):
    """Enumeration of model types for code generation."""

    GENERIC = 0
    SRN = 1
    CELL_CYCLE = 2


class TimeUnit(Enum):
    """SBML time units relevant to Chaste, which simulates in hours.

    Each member carries its size in seconds (used to match an SBML unit definition), the multiplier
    that converts a derivative from the native unit to per-hour (native units per hour), and a display
    name for logging. ``NONE`` means "no conversion" (multiplier 1.0); it is distinct from ``HOUR``
    (also multiplier 1.0) so logging can tell "already in hours" from "left unscaled".
    """

    #        (seconds_per_unit, units_per_hour, display)
    NONE = (None, 1.0, "native/unknown")
    MILLISECOND = (1.0e-3, 3.6e6, "milliseconds")
    SECOND = (1.0, 3600.0, "seconds")
    MINUTE = (60.0, 60.0, "minutes")
    HOUR = (3600.0, 1.0, "hours")

    def __init__(self, seconds_per_unit: Optional[float], multiplier: float, display: str) -> None:
        self.seconds_per_unit = seconds_per_unit
        self.multiplier = multiplier  # native units per hour; applied to derivatives (dY/d(hours))
        self.display = display

    @property
    def multiplier_literal(self) -> str:
        """C++ double literal for the multiplier, e.g. ``'60.0'``, ``'3600.0'``, ``'3600000.0'``."""
        return repr(float(self.multiplier))

    @classmethod
    def from_cli(cls, token: str) -> "TimeUnit":
        """Map a ``--timescale`` token (``ms``/``s``/``m``/``h``) to a TimeUnit."""
        return {"ms": cls.MILLISECOND, "s": cls.SECOND, "m": cls.MINUTE, "h": cls.HOUR}[token]

    @classmethod
    def from_seconds_factor(cls, factor: float) -> Optional["TimeUnit"]:
        """Map a seconds-per-unit factor to a TimeUnit, or None if it matches no known unit."""
        for unit in (cls.MILLISECOND, cls.SECOND, cls.MINUTE, cls.HOUR):
            if math.isclose(factor, unit.seconds_per_unit, rel_tol=1e-9):
                return unit
        return None
