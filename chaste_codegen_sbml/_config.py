"""Configuration settings and constants for code generation."""

import pathlib
from enum import Enum

ROOT_DIR = pathlib.Path(__file__).parent.absolute()

NON_DIM_UNITS = "non-dim"


class VarType(Enum):
    STATE_VARIABLE = 1
    DERIVED_QUANTITY = 2
    VARIABLE_PARAMETER = 3
    CONSTANT_PARAMETER = 4
    RULE_BASED_PARAMETER = 5
    ASSIGNMENT_RULE = 6
    FUNCTION = 7
    REACTION = 8
    UNKNOWN = 9


class EventType(Enum):
    CELL_DIVISION = 1
    UNKNOWN = 2


class ModelType(Enum):
    SRN = 1
    CELL_CYCLE = 2
    UNKNOWN = 3
