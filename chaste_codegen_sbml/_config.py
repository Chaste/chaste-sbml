"""Configuration settings and constants for code generation."""

import pathlib

ROOT_DIR = pathlib.Path(__file__).parent.absolute()

TAB = "    "

ODE_SYSTEM_SUFFIX = "OdeSystem"

SRN_HEADER_GUARD_SUFFIX = "ODESYSTEMANDSRNMODEL_HPP_"
SRN_MODEL_SUFFIX = "SrnModel"

CC_HEADER_GUARD_SUFFIX = "ODESYSTEMANDCELLCYCLEMODEL_HPP_"
