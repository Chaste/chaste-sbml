"""Module for the ChasteSbmlSrnModel class."""

import pathlib

from ._config import ODE_SUFFIX, SRN_SUFFIX
from .chaste_sbml_model import ChasteSbmlModel


class ChasteSbmlSrnModel(ChasteSbmlModel):
    """Class for generating Chaste code for an SRN model from SBML data."""

    # -- PUBLIC ---------------------------------------
    def __init__(self, sbml_file: str, model_name: str = "", **kwargs) -> None:
        """Initialise the ChasteSbmlSrnModel."""
        super().__init__(sbml_file, model_name, SRN_SUFFIX, **kwargs)

    @property
    def srn_cpp_filename(self) -> str:
        """Get the output {srn_model}.cpp filename."""
        return str(pathlib.Path(self.srn_hpp_filename).with_suffix(".cpp"))

    @property
    def srn_hpp_filename(self) -> str:
        """Get the output {srn_model}.hpp filename."""
        return f"{self._model_name}{ODE_SUFFIX}And{SRN_SUFFIX}Model.hpp"

    # -- PRIVATE ---------------------------------------

    def _generate(self) -> None:
        """Generate Chaste code for the SRN model.
        Overrides the parent method.
        """
        # Get the template variables
        template_vars = self._get_template_vars(self.srn_hpp_filename)

        # Generate the cpp code
        cpp_template = self._get_template("srn.cpp")
        cpp_code = cpp_template.render(template_vars)

        # Generate the hpp code
        hpp_template = self._get_template("srn.hpp")
        hpp_code = hpp_template.render(template_vars)

        # Add the output files
        self._add_output(self.srn_cpp_filename, cpp_code)
        self._add_output(self.srn_hpp_filename, hpp_code)
