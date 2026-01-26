"""Configuration settings and constants for code generation."""

import pathlib
from enum import Enum

ROOT_DIR = pathlib.Path(__file__).parent.absolute()

NON_DIM_UNITS = "non-dim"

PREFIX_SEP = "__"
AMOUNT_PREFIX = "amt"
CHASTE_PREFIX = "chaste"
INITIAL_ASSIGNMENT_PREFIX = "ia"


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


class EventType(Enum):
    """Enumeration of event types in SBML models."""

    CELL_DIVISION = 1
    UNKNOWN = 2


class ModelType(Enum):
    """Enumeration of model types for code generation."""

    GENERIC = 0
    SRN = 1
    CELL_CYCLE = 2
