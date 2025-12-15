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
    CONSTANT_PARAMETER = 2
    DERIVED_QUANTITY = 3
    FUNCTION = 4
    INITIAL_ASSIGNMENT = 5
    RATE_RULE = 6
    REACTION = 7
    REFERENCE_VARIABLE = 8
    STATE_VARIABLE = 9
    VARIABLE_PARAMETER = 10
    UNKNOWN = 11


class EventType(Enum):
    """Enumeration of event types in SBML models."""

    CELL_DIVISION = 1
    UNKNOWN = 2


class ModelType(Enum):
    """Enumeration of model types for code generation."""

    GENERIC = 0
    SRN = 1
    CELL_CYCLE = 2
