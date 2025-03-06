"""Module for the ChasteCellCycleModel class."""

from .chaste_model import ChasteModel


class ChasteCellCycleModel(ChasteModel):
    """Class for generating Chaste code for a Cell Cycle model from SBML data."""

    def __init__(self, sbml_file: str, model_name: str = None) -> None:
        """Initialise the ChasteCellCycleModel."""
        pass

    def generate_hpp(self):
        """Generate the header file for the Cell Cycle model."""
        return ""

    def generate_cpp(self):
        """Generate the source file for the Cell Cycle model."""
        return ""
