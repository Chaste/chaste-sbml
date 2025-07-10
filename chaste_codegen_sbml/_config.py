"""Configuration settings and constants for code generation."""

import pathlib
from enum import Enum

ROOT_DIR = pathlib.Path(__file__).parent.absolute()

CCM_SUFFIX = "CellCycle"
ODE_SUFFIX = "OdeSystem"
SRN_SUFFIX = "Srn"

NON_DIM_UNITS = "non-dim"


class VarType(Enum):
    STATE_VARIABLE = 1
    DERIVED_QUANTITY = 2
    VARIABLE_PARAMETER = 3
    CONSTANT_PARAMETER = 4
    ASSIGNMENT_RULE = 5
    FUNCTION = 6
    REACTION = 7
    UNKNOWN = 8
