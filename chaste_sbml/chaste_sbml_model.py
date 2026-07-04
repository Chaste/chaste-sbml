"""Module for the ChasteSbmlModel class."""

import abc
import os
from typing import TYPE_CHECKING

from ._config import ModelType
from ._model_builder import ModelBuilder
from ._names import NameManager
from ._rendering import CodeRenderer
from ._sbml_loader import load_sbml_model
from ._utils import generate_header_guard, to_camel_case, to_cpp_name

if TYPE_CHECKING:
    from typing import Any


class ChasteSbmlModel:
    """Generates Chaste C++ code from an SBML model.

    Orchestrates the pipeline: load the SBML (``_sbml_loader``), resolve identifier conflicts
    and reserve synthetic names (``NameManager``), build the internal representation
    (``ModelBuilder``), then render and write the C++ (``CodeRenderer``).
    """

    __metaclass__ = abc.ABCMeta

    # -- PUBLIC --------------------------------------

    def __init__(self, sbml_file: str, model_name: str = "", model_type: ModelType = ModelType.GENERIC) -> None:
        """Initialise the ChasteSbmlModel.

        :param sbml_file: The SBML file to generate code from.
        :param model_name: The model name; derived from the filename when not given.
        :param model_type: The model type e.g. ModelType.SRN.
        """
        self._sbml_file = os.path.abspath(sbml_file)
        if not os.path.isfile(self._sbml_file):
            raise FileNotFoundError(f"Could not find SBML file: {self._sbml_file}")

        if model_name:
            self._model_name = model_name
        else:
            filename = os.path.splitext(os.path.basename(self._sbml_file))[0]
            model_name = to_camel_case(to_cpp_name(filename)) + "Sbml"
            self._model_name = model_name[0].upper() + model_name[1:]

        self._model_type = model_type

        self._ode_class_name = self._model_name + "OdeSystem"
        self._ode_hpp_filename = f"{self._ode_class_name}.hpp"
        self._ode_cpp_filename = f"{self._ode_class_name}.cpp"

        self._srn_class_name = ""
        self._srn_hpp_filename = ""
        self._srn_cpp_filename = ""
        if self._model_type == ModelType.SRN:
            self._srn_class_name = self._model_name + "SrnModel"
            self._srn_hpp_filename = f"{self._srn_class_name}.hpp"
            self._srn_cpp_filename = f"{self._srn_class_name}.cpp"

        self._cell_cycle_class_name = ""
        self._cell_cycle_hpp_filename = ""
        self._cell_cycle_cpp_filename = ""
        if self._model_type == ModelType.CELL_CYCLE:
            self._cell_cycle_class_name = self._model_name + "CellCycleModel"
            self._cell_cycle_hpp_filename = f"{self._cell_cycle_class_name}.hpp"
            self._cell_cycle_cpp_filename = f"{self._cell_cycle_class_name}.cpp"

        # Read, flatten and convert the SBML model
        self._sbml_model = load_sbml_model(self._sbml_file)

        self._names = NameManager(self._sbml_model)

        self._template_vars = {}  # type: dict[str, Any]

        self._outputs = {}  # { filename: code }
        self._renderer = CodeRenderer()

        self._names.resolve_real_id_conflicts()

        self._builder = ModelBuilder(self._sbml_model, self._names, self._model_name)
        self._builder.build()
        self._populate_template_vars()

    @property
    def outputs(self) -> dict[str, str]:
        """Get the generated code outputs.

        :return: A dictionary of filename and code pairs.
        """
        return self._outputs

    def write(self, output_directory=None):
        """Generate Chaste code and write to file.

        :param output_directory: The output directory. Defaults to the current directory.
        """
        # Generate the code
        self._generate_outputs()

        # Write the code to file (formatted with clang-format)
        self._renderer.write(self._outputs, output_directory)

    # -- PRIVATE ---------------------------------------

    def _add_output(self, filename: str, code: str) -> None:
        """Add generated code to the outputs dictionary.

        :param filename: The filename.
        :param code: The code.
        """
        self._outputs[filename] = code

    def _generate_output(self, template_path, filename) -> None:
        """Generate a single output file from a template.

        Generated code is stored in the outputs list.

        :param template_path: The path to the template.
        :param filename: The output filename.
        """
        code = self._renderer.render(template_path, self._template_vars)
        self._add_output(filename, code)

    def _generate_outputs(self) -> None:
        """Generate Chaste code for the model."""
        # Generate code for the OdeSystem
        self._generate_output("ode/ode.hpp", self._ode_hpp_filename)
        self._generate_output("ode/ode.cpp", self._ode_cpp_filename)

        # Generate code for the SRN or Cell-Cycle model
        if self._model_type == ModelType.SRN:
            self._generate_output("srn/srn.hpp", self._srn_hpp_filename)
            self._generate_output("srn/srn.cpp", self._srn_cpp_filename)

        elif self._model_type == ModelType.CELL_CYCLE:
            self._generate_output("cell_cycle/cell_cycle.hpp", self._cell_cycle_hpp_filename)
            self._generate_output("cell_cycle/cell_cycle.cpp", self._cell_cycle_cpp_filename)

    def _populate_template_vars(self) -> None:
        """Populate the template variables for generating C++ code."""
        template_vars: dict[str, "Any"] = dict(
            **self._builder.template_data(),
            ode_class_name=self._ode_class_name,
            ode_header_guard=generate_header_guard(self._ode_hpp_filename),
            ode_hpp_file=self._ode_hpp_filename,
        )

        if self._model_type == ModelType.SRN:
            template_vars.update(
                dict(
                    srn_class_name=self._srn_class_name,
                    srn_header_guard=generate_header_guard(self._srn_hpp_filename),
                    srn_hpp_file=self._srn_hpp_filename,
                )
            )
        elif self._model_type == ModelType.CELL_CYCLE:
            template_vars.update(
                dict(
                    cell_cycle_class_name=self._cell_cycle_class_name,
                    cell_cycle_header_guard=generate_header_guard(self._cell_cycle_hpp_filename),
                    cell_cycle_hpp_file=self._cell_cycle_hpp_filename,
                )
            )
        self._template_vars = template_vars
