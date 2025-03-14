"""Module for the ChasteSbmlSrnModel class."""

import pathlib

from ._config import ODE_SUFFIX, SRN_SUFFIX
from .chaste_sbml_model import ChasteSbmlModel


class ChasteSbmlSrnModel(ChasteSbmlModel):
    """Class for generating Chaste code for an SRN model from SBML data."""

    # -- PUBLIC ---------------------------------------
    def __init__(self, sbml: str, name: str = None, **kwargs) -> None:
        """Initialise the ChasteSbmlSrnModel."""
        super().__init__(sbml, name, SRN_SUFFIX, **kwargs)

    @property
    def srn_cpp_filename(self) -> str:
        """Get the output {srn_model}.cpp filename."""
        return pathlib.Path(self.srn_hpp_filename).with_suffix(".cpp")

    @property
    def srn_hpp_filename(self) -> str:
        """Get the output {srn_model}.hpp filename."""
        return f"{self._model_name}{ODE_SUFFIX}And{SRN_SUFFIX}Model.hpp"

    # -- PRIVATE ---------------------------------------

    def _generate(self) -> None:
        """Generate Chaste code for the SRN model.
        Overrides the parent method.
        """
        self._generate_srn_hpp()
        self._generate_srn_cpp()

    def _generate_srn_cpp(self) -> None:
        """Generate the cpp file code for the SRN model."""
        template = self._get_template("srn.cpp")
        vars = self._get_srn_cpp_vars()
        code = template.render(vars)
        self._add_output(self.srn_cpp_filename, code)

    def _generate_srn_hpp(self) -> None:
        """Generate the hpp file code for the SRN model."""
        template = self._get_template("srn.hpp")
        vars = self._get_hpp_vars(self.srn_hpp_filename)
        code = template.render(vars)
        self._add_output(self.srn_hpp_filename, code)

    def _get_srn_cpp_vars(self) -> dict[str, str]:
        """Generate the template variables for the SRN model cpp file.

        return: The generated source file as a string.
        """
        # Events
        num_events = self._model.getNumEvents()
        event_vector_init_str = ""
        if num_events > 0:
            event_vector_init_str = f"eventsSatisfied.resize({num_events}, false);"

        # Create inputs for the source file template
        cpp_vars = dict(
            compartments=self._format_compartments(),
            event_vector_init=event_vector_init_str,
            function_definitions=self._format_function_definitions(),
            model_hpp_file=self.srn_hpp_filename,
            model_class_name=self._model_class_name,
            ode_class_name=self._ode_class_name,
            parameters=self._format_parameters(),
            reactions=self._format_reactions(),
            rules=self._format_rules(),
            num_state_vars=self._num_state_vars,
            num_events=self._model.getNumEvents(),
            species=self._format_species(),
            wrapper_class_name=self._wrapper_class_name,
        )

        return cpp_vars
