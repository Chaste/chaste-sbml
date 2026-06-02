"""Configuration settings and constants for code generation."""

import pathlib
from enum import Enum

ROOT_DIR = pathlib.Path(__file__).parent.absolute()

NON_DIM_UNITS = "non-dim"

PREFIX_SEP = "__"
AMOUNT_PREFIX = "amt"
CHASTE_PREFIX = "chaste"
INITIAL_ASSIGNMENT_PREFIX = "ia"

DERIVATIVE_PREFIX = "d_"
DERIVATIVE_SUFFIX = "_dt"


class EquationType(Enum):
    """Enumeration of equation types for code generation."""

    AMOUNT = 0
    ASSIGNMENT_RULE = 1
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

    AMOUNT = 0
    DERIVED_QUANTITY = 1
    FUNCTION = 2
    PARAMETER = 3
    REACTION = 4
    STATE_VARIABLE = 5
    UNKNOWN = 6


class ModelType(Enum):
    """Enumeration of model types for code generation."""

    GENERIC = 0
    SRN = 1
    CELL_CYCLE = 2
