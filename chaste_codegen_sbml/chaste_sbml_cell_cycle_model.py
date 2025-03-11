"""Module for the ChasteCellCycleModel class."""

from ._config import CCM_HEADER_GUARD_SUFFIX, ODE_SYSTEM_SUFFIX, CCM_MODEL_SUFFIX, TAB
from ._utils import (
    get_function_definition_arguments,
)
from .chaste_sbml_model import ChasteSbmlModel


class ChasteSbmlCellCycleModel(ChasteSbmlModel):
    """Class for generating Chaste code for a Cell Cycle model from SBML data."""

    # -- PUBLIC ---------------------------------------

    def __init__(self, sbml_file: str, model_name: str = None, **kwargs) -> None:
        """Initialise the ChasteSbmlCellCycleModel."""
        super().__init__(sbml_file, model_name, **kwargs)

        self._ccm_name = self._model_name + CCM_MODEL_SUFFIX

        ccm_filename = f"{self._model_name}{ODE_SYSTEM_SUFFIX}And{CCM_MODEL_SUFFIX}"
        self._ccm_hpp_filename = f"{ccm_filename}.hpp"
        self._ccm_cpp_filename = f"{ccm_filename}.cpp"

        ccm_hpp_template = self._get_template("ccm.hpp")
        ccm_hpp_vars = self._get_ccm_hpp_vars()
        ccm_hpp_code = ccm_hpp_template.render(ccm_hpp_vars)
        self._outputs[self._ccm_hpp_filename] = ccm_hpp_code

    @property
    def ccm_cpp_filename(self) -> str:
        """Get the output {cell_cycle_model}.cpp filename."""
        return self._ccm_cpp_filename

    @property
    def ccm_hpp_filename(self) -> str:
        """Get the output {cell_cycle_model}.hpp filename."""
        return self._ccm_hpp_filename

    # -- PRIVATE ---------------------------------------

    def _get_ccm_hpp_vars(self) -> dict[str, str]:
        """Generate the template variables for the Cell Cycle model hpp file.

        return: The generated header file as a string.
        """
        # Get inputs for the header file template
        header_guard_str = self._model_name.upper() + CCM_HEADER_GUARD_SUFFIX

        compartments = self._compartments
        compartment_decls_str = f"\n{TAB}".join(
            [f"double {c.getId()}; // {self._get_varname(c)}" for c in compartments]
        )

        parameters = self._parameters
        parameter_decls_str = f"\n{TAB}".join(
            [f"double {p.getId()}; // {self._get_varname(p)}" for p in parameters]
        )

        function_decls_list = []
        for fn in self._function_definitions:
            fn_id = fn.getId()
            fn_args = get_function_definition_arguments(fn)

            function_decls_list.append(f"double {fn_id}({', '.join(fn_args)});")

        function_decls_str = f"\n{TAB}".join(function_decls_list)

        # Create inputs for the header file template
        hpp_vars = dict(
            header_guard=header_guard_str,
            ode_system_name=self._ode_system_name,
            cell_cycle_model_name=self._ccm_name,
            num_state_vars=self._num_state_vars,
            compartment_decls=compartment_decls_str,
            parameter_decls=parameter_decls_str,
            function_decls=function_decls_str,
        )

        return hpp_vars
