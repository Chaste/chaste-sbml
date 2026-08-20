"""Module for the ChasteSbmlModel class."""

import abc
import logging
import os
from dataclasses import dataclass
from typing import TYPE_CHECKING, Optional

from ._config import ModelType, TimeUnit
from ._model_builder import ModelBuilder
from ._names import NameManager, generate_header_guard, to_camel_case, to_cpp_name
from ._rendering import CodeRenderer
from ._sbml_reader import load_sbml_model

if TYPE_CHECKING:
    from typing import Any

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class _WrapperSpec:
    """Descriptor for a :class:`ModelType` that wraps the ODE system in a Chaste model class.

    Collects the per-type facts that would otherwise be spread across parallel ``if model_type ==``
    branches, so ``ChasteSbmlModel`` can drive naming, output selection and template variables from a
    single lookup, and the SRN and cell-cycle wrappers can share one pair of templates
    (``templates/wrapper``). ``ModelType.GENERIC`` has no wrapper (it generates only the ODE system)
    and so has no entry in :data:`_WRAPPER_SPECS`.

    :param class_suffix: Appended to the model name to form the wrapper class name (e.g. ``SrnModel``).
    :param base_class: The Chaste base class the wrapper derives from (e.g. ``AbstractSbmlSrnModel``).
    :param abstract_class: The abstract type returned by the ``Create*Model`` builder method
        (e.g. ``AbstractSrnModel``).
    :param create_method: The overridden builder method name (e.g. ``CreateSrnModel``).
    :param output_method: The overridden parameter-output method name (e.g. ``OutputSrnModelParameters``).
    :param model_noun: How the model is referred to in doc comments (e.g. ``SRN model``).
    """

    class_suffix: str
    base_class: str
    abstract_class: str
    create_method: str
    output_method: str
    model_noun: str


# Wrapper descriptor per model type. GENERIC is absent: it generates only the ODE system.
_WRAPPER_SPECS = {
    ModelType.SRN: _WrapperSpec(
        class_suffix="SrnModel",
        base_class="AbstractSbmlSrnModel",
        abstract_class="AbstractSrnModel",
        create_method="CreateSrnModel",
        output_method="OutputSrnModelParameters",
        model_noun="SRN model",
    ),
    ModelType.CELL_CYCLE: _WrapperSpec(
        class_suffix="CellCycleModel",
        base_class="AbstractSbmlCellCycleModel",
        abstract_class="AbstractCellCycleModel",
        create_method="CreateCellCycleModel",
        output_method="OutputCellCycleModelParameters",
        model_noun="cell-cycle model",
    ),
}


class ChasteSbmlModel:
    """Generates Chaste C++ code from an SBML model.

    Orchestrates the pipeline:

    - Load the SBML (``_sbml_reader``)
    - Resolve identifier conflicts and reserve synthetic names (``NameManager``)
    - Build the internal representation (``ModelBuilder``)
    - Render and write the C++ (``CodeRenderer``)
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

        # The SRN / cell-cycle wrapper class around the ODE system (none for a GENERIC model).
        self._wrapper_spec = _WRAPPER_SPECS.get(self._model_type)
        self._wrapper_class_name = ""
        self._wrapper_hpp_filename = ""
        self._wrapper_cpp_filename = ""
        if self._wrapper_spec is not None:
            self._wrapper_class_name = self._model_name + self._wrapper_spec.class_suffix
            self._wrapper_hpp_filename = f"{self._wrapper_class_name}.hpp"
            self._wrapper_cpp_filename = f"{self._wrapper_class_name}.cpp"

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
        """Resolve the model's native time unit into Chaste's hours.

        Precedence:

        - An explicit ``--timescale`` override wins (warning if it contradicts a declared unit)
        - Otherwise a declared unit is used
        - Otherwise the SBML default applies:

          - Level 2 predefines ``time`` as seconds
          - Level 3 (and any other unsupported level) leaves it undefined (no conversion).

        When no unit is declared and no override is given, a warning hints at the ``--timescale`` option.

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
        elif sbml_level == 3:
            # SBML Level 3 leaves an unset time unit undefined, so apply no conversion.
            resolved = TimeUnit.NONE
            logger.warning(
                "The model declares no time unit (SBML Level 3 leaves it undefined); no time "
                "conversion is applied. Pass --timescale ms|s|m|h to convert the model to Chaste "
                "hours.",
            )
        else:
            # Unsupported SBML level (already flagged when loading); do not guess a unit.
            resolved = TimeUnit.NONE
            logger.warning(
                "The model declares no time unit and its SBML level (%d) is unsupported; no time "
                "conversion is applied. Pass --timescale ms|s|m|h to convert the model to Chaste "
                "hours.",
                sbml_level,
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

        # Generate code for the SRN or cell-cycle wrapper, if this model type has one
        if self._wrapper_spec is not None:
            self._generate_output("wrapper/wrapper.hpp", self._wrapper_hpp_filename)
            self._generate_output("wrapper/wrapper.cpp", self._wrapper_cpp_filename)

        # Generate a placeholder test skeleton for the model
        if self._generate_tests:
            self._generate_test_output()

    def _populate_template_vars(self) -> None:
        """Populate the template variables for generating C++ code."""
        template_vars: dict[str, "Any"] = dict(
            **self._builder.template_data(),
            model_name=self._model_name,
            ode_class_name=self._ode_class_name,
            ode_header_guard=generate_header_guard(self._ode_hpp_filename),
            ode_hpp_file=self._ode_hpp_filename,
            test_header_guard=generate_header_guard(self._test_hpp_filename),
            model_type=self._model_type,
            scale_time=self._time_unit.multiplier != 1.0,
            time_multiplier=self._time_unit.multiplier_literal,
            time_unit_display=self._time_unit.display,
        )

        if self._wrapper_spec is not None:
            template_vars.update(
                dict(
                    wrapper_class_name=self._wrapper_class_name,
                    wrapper_header_guard=generate_header_guard(self._wrapper_hpp_filename),
                    wrapper_hpp_file=self._wrapper_hpp_filename,
                    wrapper_base_class=self._wrapper_spec.base_class,
                    wrapper_abstract_class=self._wrapper_spec.abstract_class,
                    wrapper_create_method=self._wrapper_spec.create_method,
                    wrapper_output_method=self._wrapper_spec.output_method,
                    wrapper_model_noun=self._wrapper_spec.model_noun,
                )
            )
        self._template_vars = template_vars
