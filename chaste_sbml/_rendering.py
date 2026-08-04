"""Rendering of generated Chaste C++ source from the packaged Jinja templates.

Holds the Jinja environment and the mechanics of turning a template plus a context into a
source string, and of writing the results to disk (formatted with clang-format). What data to
render and which templates to use stays with ChasteSbmlModel; this module is model-agnostic.
"""

import os
import pathlib
import shutil
import subprocess
import sys

from jinja2 import Environment, PackageLoader, StrictUndefined

from ._config import (
    AMOUNT_PREFIX,
    CONCENTRATION_PREFIX,
    PREFIX_SEP,
    ROOT_DIR,
    DerivedQuantityKind,
    EquationType,
    EventType,
    ModelType,
    VarType,
)


def _reject_none(value):
    """Jinja finalize hook: a ``None`` reaching the output is a bug, not the literal ``"None"``.

    StrictUndefined catches an *unprovided* variable, but a provided-yet-``None`` value (e.g. an
    unset ``initial_value`` rendered without its ``is not none`` guard) would otherwise silently
    emit the text ``None`` into the generated C++. Raise instead so it surfaces during generation.

    :param value: The value about to be rendered.
    :return: The value unchanged when it is not ``None``.
    """
    if value is None:
        raise ValueError("template rendered a None value (missing an 'is not none' guard?)")
    return value


class CodeRenderer:
    """Renders Chaste C++ source from the packaged Jinja templates and writes it to disk."""

    _env = Environment(
        loader=PackageLoader("chaste_sbml"),
        # These templates emit C++, never HTML/XML: autoescaping would turn operators like < > &
        # into &lt; &gt; &amp; and corrupt the output, so it is always off (never select_autoescape,
        # whose filename heuristic would silently enable it for e.g. an .xml template).
        autoescape=False,
        trim_blocks=True,
        lstrip_blocks=True,
        undefined=StrictUndefined,  # a template referencing an unprovided variable is an error
        finalize=_reject_none,  # a None reaching the output is a bug, not the literal "None"
    )
    _env.globals["AMOUNT_PREFIX"] = AMOUNT_PREFIX
    _env.globals["CONCENTRATION_PREFIX"] = CONCENTRATION_PREFIX
    _env.globals["DerivedQuantityKind"] = DerivedQuantityKind
    _env.globals["EquationType"] = EquationType
    _env.globals["EventType"] = EventType
    _env.globals["ModelType"] = ModelType
    _env.globals["PREFIX_SEP"] = PREFIX_SEP
    _env.globals["VarType"] = VarType

    def render(self, template_name: str, context: dict) -> str:
        """Render a template to a C++ source string.

        :param template_name: The template path within the package (e.g. ``ode/ode.cpp``).
        :param context: The template variables.
        :return: The rendered source.
        """
        return self._env.get_template(template_name).render(context)

    @staticmethod
    def write(outputs: dict, output_directory=None) -> None:
        """Write each output file and format it with clang-format.

        :param outputs: A mapping of filename to source code.
        :param output_directory: The directory to write into; defaults to the current directory.
        """
        root_dir = pathlib.Path(output_directory) if output_directory else pathlib.Path.cwd()

        for filename, code in outputs.items():
            with open(root_dir / filename, "w") as f:
                f.write(code)

        # Format with clang-format — search the current interpreter's bin/ first
        # so the clang-format installed as a Python dependency is always found.
        python_bin = os.path.dirname(sys.executable)
        search_path = os.pathsep.join([python_bin, os.environ.get("PATH", "")])
        clang_format = shutil.which("clang-format", path=search_path)
        if clang_format is not None:
            for filename in outputs:
                subprocess.run(
                    [clang_format, "-i", f"-style=file:{ROOT_DIR}/.clang-format", str(root_dir / filename)],
                    check=True,
                )
