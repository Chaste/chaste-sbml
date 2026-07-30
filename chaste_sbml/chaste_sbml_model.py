"""Module for the ChasteSbmlModel class."""

import abc
import logging
import os
from typing import TYPE_CHECKING, Optional

from ._config import ModelType, TimeUnit
from ._model_builder import ModelBuilder
from ._names import NameManager, generate_header_guard, to_camel_case, to_cpp_name
from ._rendering import CodeRenderer
from ._sbml_reader import load_sbml_model

if TYPE_CHECKING:
    from typing import Any

logger = logging.getLogger(__name__)


class ChasteSbmlModel:
    """Generates Chaste C++ code from an SBML model.

    Orchestrates the pipeline: load the SBML (``_sbml_reader``), resolve identifier conflicts
    and reserve synthetic names (``NameManager``), build the internal representation
    (``ModelBuilder``), then render and write the C++ (``CodeRenderer``).
    """

    __metaclass__ = abc.ABCMeta

    def __init__(
        self,
        sbml_file: str,
        model_name: str = "",
        model_type: ModelType = ModelType.GENERIC,
        generate_tests: bool = True,
        time_unit: Optional[TimeUnit] = None,
    ) -> None:
        """Initialise the ChasteSbmlModel.

        :param sbml_file: The SBML file to generate code from.
        :param model_name: The model name; derived from the filename when not given.
        :param model_type: The model type e.g. ModelType.SRN.
        :param generate_tests: Whether to generate a placeholder test for the model.
        :param time_unit: The model's native time unit, overriding auto-detection. ``None``
            (the default) auto-detects from the SBML; pass ``TimeUnit.NONE`` to force no conversion.
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

        self._generate_tests = generate_tests
        self._test_hpp_filename = f"Test{self._model_name}.hpp"

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
        self._sbml_model, declared_time_unit, sbml_level = load_sbml_model(self._sbml_file)

        # Resolve the time unit used to scale derivatives to Chaste's hours (see _resolve_time_unit).
        self._time_unit = self._resolve_time_unit(time_unit, declared_time_unit, sbml_level)

        self._names = NameManager(self._sbml_model)

        self._template_vars = {}  # type: dict[str, Any]

        self._outputs = {}  # { filename: code }
        self._test_outputs = {}  # { filename: code } for generated placeholder tests
        self._renderer = CodeRenderer()

        self._names.resolve_real_id_conflicts()

        self._builder = ModelBuilder(self._sbml_model, self._names, self._model_name)
        self._builder.build()
        self._populate_template_vars()

    def _resolve_time_unit(
        self, override: Optional[TimeUnit], declared: Optional[TimeUnit], sbml_level: int
    ) -> TimeUnit:
        """Resolve the time unit to convert derivatives from, into Chaste's hours.

        Precedence: an explicit override wins (warning if it contradicts a declared unit); otherwise a
        declared unit is used; otherwise the SBML default applies -- Level 2 predefines ``time`` as
        seconds, while Level 3 leaves an unset time unit undefined (no conversion). When no unit is
        declared and no override is given, a warning hints at the ``--timescale`` option.

        :param override: An explicit time unit (from ``--timescale``), or ``None`` to auto-detect.
        :param declared: The time unit declared by the model, or ``None`` if none was declared.
        :param sbml_level: The model's SBML level.
        :return: The resolved :class:`TimeUnit`.
        """
        if override is not None:
            if declared is not None and declared is not override:
                logger.warning(
                    "Time unit override (%s) differs from the unit declared by the model (%s); "
                    "using the override.",
                    override.display,
                    declared.display,
                )
            resolved = override
        elif declared is not None:
            resolved = declared
        elif sbml_level == 2:
            # SBML Level 2 predefines the time unit as seconds when none is declared.
            resolved = TimeUnit.SECOND
            logger.warning(
                "The model declares no time unit; assuming seconds per the SBML Level 2 default and "
                "scaling derivatives by %s to convert to Chaste hours. Pass --timescale ms|s|m|h to "
                "set the model's time unit explicitly.",
                resolved.multiplier,
            )
        else:
            # SBML Level 3 leaves an unset time unit undefined, so apply no conversion.
            resolved = TimeUnit.NONE
            logger.warning(
                "The model declares no time unit (SBML Level 3 leaves it undefined); no time "
                "conversion is applied. Pass --timescale ms|s|m|h to convert the model to Chaste "
                "hours.",
            )

        logger.info("Using time unit %s (derivative multiplier %s).", resolved.display, resolved.multiplier)
        return resolved

    @property
    def outputs(self) -> dict[str, str]:
        """Get the generated code outputs.

        :return: A copy of the filename to code mapping, so callers cannot mutate the internal state.
        """
        return dict(self._outputs)

    @property
    def test_outputs(self) -> dict[str, str]:
        """Get the generated placeholder test outputs.

        :return: A copy of the filename to code mapping, so callers cannot mutate the internal state.
        """
        return dict(self._test_outputs)

    def write(self, output_directory=None, test_output_directory=None):
        """Generate Chaste code and write to file.

        :param output_directory: The output directory for the model code.
            Defaults to the current directory.
        :param test_output_directory: The output directory for the generated
            placeholder tests. Defaults to ``output_directory``.
        """
        # Generate the code
        self._generate_outputs()

        # Write the model code to file (formatted with clang-format)
        self._renderer.write(self._outputs, output_directory)

        # Write the placeholder tests, defaulting to the model output directory
        if self._test_outputs:
            test_dir = test_output_directory if test_output_directory is not None else output_directory
            self._renderer.write(self._test_outputs, test_dir)

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

    def _generate_test_output(self) -> None:
        """Generate the placeholder test file for the model.

        Rendered into the separate test outputs so it can be written to its own directory.
        """
        code = self._renderer.render("test/test.hpp", self._template_vars)
        self._test_outputs[self._test_hpp_filename] = code

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

        # Generate a placeholder test skeleton for the model
        if self._generate_tests:
            self._generate_test_output()

    def _populate_template_vars(self) -> None:
        """Populate the template variables for generating C++ code."""
        template_vars: dict[str, "Any"] = dict(
            **self._builder.template_data(),
            ode_class_name=self._ode_class_name,
            ode_header_guard=generate_header_guard(self._ode_hpp_filename),
            ode_hpp_file=self._ode_hpp_filename,
            test_header_guard=generate_header_guard(self._test_hpp_filename),
            model_type=self._model_type,
            scale_time=self._time_unit.multiplier != 1.0,
            time_multiplier=self._time_unit.multiplier_literal,
            time_unit_display=self._time_unit.display,
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
