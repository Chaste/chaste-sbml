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

    ASSIGNMENT_RULE = 0
    DERIVATIVE = 1
    INITIAL_ASSIGNMENT = 2
    INITIAL_VALUE = 3
    REACTION = 4
    UNKNOWN = 5


class EventType(Enum):
    """Enumeration of event types in SBML models."""

    CELL_DIVISION = 1
    UNKNOWN = 2


class VarType(Enum):
    """Enumeration of variable types in SBML models."""

    AMOUNT = 0
    ASSIGNMENT_RULE = 1
    DERIVED_QUANTITY = 2
    FUNCTION = 3
    INITIAL_ASSIGNMENT = 4
    PARAMETER = 5
    RATE_RULE = 6
    REACTION = 7
    STATE_VARIABLE = 8
    STOICHIOMETRY_VARIABLE = 9
    UNKNOWN = 10


class ModelType(Enum):
    """Enumeration of model types for code generation."""

    GENERIC = 0
    SRN = 1
    CELL_CYCLE = 2
