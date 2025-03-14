"""Module for the ChasteCellCycleModel class."""

import pathlib

from libsbml import formulaToString

from ._config import CCM_SUFFIX, ODE_SUFFIX
from ._utils import (
    convert_formula,
    sort_nodes,
)
from .chaste_sbml_model import ChasteSbmlModel


class ChasteSbmlCellCycleModel(ChasteSbmlModel):
    """Class for generating Chaste code for a Cell Cycle model from SBML data."""

    # -- PUBLIC ---------------------------------------

    def __init__(self, sbml: str, name: str = None, **kwargs) -> None:
        """Initialise the ChasteSbmlCellCycleModel."""
        super().__init__(sbml, name, CCM_SUFFIX, **kwargs)

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
        """Generate Chaste code for the Cell Cycle model.
        Overrides the parent method.
        """
        self._generate_ccm_hpp()
        self._generate_ccm_cpp()

    def _generate_ccm_cpp(self) -> None:
        """Generate the cpp file code for the Cell Cycle model."""
        template = self._get_template("ccm.cpp")
        vars = self._get_ccm_cpp_vars()
        code = template.render(vars)
        self._add_output(self.ccm_cpp_filename, code)

    def _generate_ccm_hpp(self) -> None:
        """Generate the hpp file code for the Cell Cycle model."""
        template = self._get_template("ccm.hpp")
        vars = self._get_hpp_vars(self.ccm_hpp_filename)
        code = template.render(vars)
        self._add_output(self.ccm_hpp_filename, code)

    def _get_ccm_cpp_vars(self) -> dict[str, str]:
        """Generate the template variables for the Cell Cycle model cpp file.

        return: The generated source file as a string.
        """
        # Events
        num_events = self._model.getNumEvents()
        event_vector_init_str = ""
        if num_events > 0:
            event_vector_init_str = f"eventsSatisfied.resize({num_events}, false);"

        event_def_tpl = """
    if ({trigger})
    {{
        {event_assignment}
        eventsSatisfied[{i}] = true;
    }}
"""

        event_defs = []
        for i, event in enumerate(self._events):
            trigger = event.getTrigger()
            trigger_nodes = sort_nodes(trigger.getMath())
            tokens = []
            for node in trigger_nodes:
                if node.isNumber():
                    token = convert_formula(str(node.getValue()))
                else:
                    token = convert_formula(node.getName())
                    if node.isName():
                        # Replace species variable name with Chaste equivalent.
                        for j, s_id in enumerate(self._odes_dict):
                            if token == s_id:
                                token = f"rY[{j}]"
                                break
                tokens.append(token)
            trigger_def_str = " ".join(tokens)

            assignment_defs = []
            for assignment in event.getListOfEventAssignments():
                # Replace species variable name with Chaste equivalent.
                variable = assignment.getVariable()
                for j, s_id in enumerate(self._odes_dict):
                    if variable == s_id:
                        variable = f"this->rGetStateVariables()[{j}]"
                        break

                formula = convert_formula(formulaToString(assignment.getMath()))
                formula_tokens = formula.split(" ")
                tokens = []
                for token in formula_tokens:
                    for j, s_id in enumerate(self._odes_dict):
                        if token == s_id:
                            token = f"rY[{j}]"
                            break
                    tokens.append(token)
                value_str = " ".join(tokens)

                assignment_defs.append(f"{variable} = double({value_str});")
            assignment_def_str = "\n".join(assignment_defs)

            event_defs.append(
                event_def_tpl.format(
                    trigger=trigger_def_str, event_assignment=assignment_def_str, i=i
                )
            )

        event_def_str = "\n".join(event_defs)

        # Create inputs for the source file template
        cpp_vars = dict(
            compartments=self._format_compartments(),
            event_defs=event_def_str,
            event_vector_init=event_vector_init_str,
            function_definitions=self._format_function_definitions(),
            model_hpp_file=self.ccm_hpp_filename,
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
