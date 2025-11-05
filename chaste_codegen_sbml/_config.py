"""Configuration settings and constants for code generation."""

import pathlib
from enum import Enum

ROOT_DIR = pathlib.Path(__file__).parent.absolute()

NON_DIM_UNITS = "non-dim"

PREFIX_SEP = "__"
AMOUNT_PREFIX = "amt"


class VarType(Enum):
    """Enumeration of variable types in SBML models."""

    STATE_VARIABLE = 0
    DERIVED_QUANTITY = 1
    AMOUNT = 2
    VARIABLE_PARAMETER = 3
    CONSTANT_PARAMETER = 4
    ASSIGNMENT_RULE = 5
    INITIAL_ASSIGNMENT = 6
    FUNCTION = 7
    REACTION = 8
    UNKNOWN = 9


class EventType(Enum):
    """Enumeration of event types in SBML models."""

    CELL_DIVISION = 1
    UNKNOWN = 2


class ModelType(Enum):
    """Enumeration of model types for code generation."""

    GENERIC = 0
    SRN = 1
    CELL_CYCLE = 2
