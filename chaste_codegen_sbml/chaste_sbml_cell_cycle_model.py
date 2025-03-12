"""Module for the ChasteCellCycleModel class."""

import pathlib

from ._config import CCM_SUFFIX, ODE_SUFFIX
from .chaste_sbml_model import ChasteSbmlModel


class ChasteSbmlCellCycleModel(ChasteSbmlModel):
    """Class for generating Chaste code for a Cell Cycle model from SBML data."""

    # -- PUBLIC ---------------------------------------

    def __init__(self, sbml_file: str, model_name: str = None, **kwargs) -> None:
        """Initialise the ChasteSbmlCellCycleModel."""
        super().__init__(sbml_file, model_name, **kwargs)

    @property
    def ccm_cpp_filename(self) -> str:
        """Get the output {cell_cycle_model}.cpp filename."""
        return pathlib.Path(self.ccm_hpp_filename).with_suffix(".cpp")

    @property
    def ccm_hpp_filename(self) -> str:
        """Get the output {cell_cycle_model}.hpp filename."""
        return f"{self._model_name}{ODE_SUFFIX}And{CCM_SUFFIX}Model.hpp"

    # -- PRIVATE ---------------------------------------

    def _generate(self) -> None:
        """Generate Chaste code for the Cell Cycle model."""
        self._generate_ccm_hpp()
        self._generate_ccm_cpp()

    def _generate_ccm_cpp(self) -> None:
        """Generate the cpp file code for the Cell Cycle model."""
        pass

    def _generate_ccm_hpp(self) -> None:
        """Generate the hpp file code for the Cell Cycle model."""
        template = self._get_template("ccm.hpp")
        vars = self._get_hpp_vars(CCM_SUFFIX, self.ccm_hpp_filename)
        code = template.render(vars)
        self._outputs[self.ccm_hpp_filename] = code
