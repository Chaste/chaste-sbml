"""Main module for the Chaste code generation from SBML package."""

from ._version import __version__
from .chaste_sbml_model import ChasteSbmlModel

__all__ = [
    "__version__",
    "ChasteSbmlModel",
]
