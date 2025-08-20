"""Module for the ChasteCellCycleModel class."""

import pathlib

from ._config import CCM_SUFFIX, ODE_SUFFIX
from .chaste_sbml_model import ChasteSbmlModel


class ChasteSbmlCellCycleModel(ChasteSbmlModel):
    """Class for generating Chaste code for a Cell Cycle model from SBML data."""

    # -- PUBLIC ---------------------------------------

    def __init__(self, sbml_file: str, model_name: str = "", **kwargs) -> None:
        """Initialise the ChasteSbmlCellCycleModel."""
        super().__init__(sbml_file, model_name, CCM_SUFFIX, **kwargs)

    @property
    def ccm_cpp_filename(self) -> str:
        """Get the output {cell_cycle_model}.cpp filename."""
        return str(pathlib.Path(self.ccm_hpp_filename).with_suffix(".cpp"))

    @property
    def ccm_hpp_filename(self) -> str:
        """Get the output {cell_cycle_model}.hpp filename."""
        return f"{self._model_name}{ODE_SUFFIX}And{CCM_SUFFIX}Model.hpp"

    # -- PRIVATE ---------------------------------------

    def _generate(self) -> None:
        """Generate Chaste code for the Cell Cycle model.
        Overrides the parent method.
        """
        # Get the template variables
        template_vars = self._get_template_vars(self.ccm_hpp_filename)

        # Generate the cpp code
        cpp_template = self._get_template("ccm.cpp")
        cpp_code = cpp_template.render(template_vars)

        # Generate the hpp code
        hpp_template = self._get_template("ccm.hpp")
        hpp_code = hpp_template.render(template_vars)

        # Add the output files
        self._add_output(self.ccm_cpp_filename, cpp_code)
        self._add_output(self.ccm_hpp_filename, hpp_code)
