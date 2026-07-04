"""Tests for console script."""

import logging
import subprocess

from chaste_sbml import __version__

logger = logging.getLogger(__name__)


def test_help():
    """The top-level help lists both subcommands."""
    output = subprocess.check_output(["chaste-sbml", "-h"]).decode("ascii")

    assert "generate" in output
    assert "copy-base-classes" in output


def test_version():
    """Test version message."""
    expected = f"chaste-sbml {__version__}"

    output = subprocess.check_output(["chaste-sbml", "--version"]).decode("ascii").strip()

    assert output == expected


def test_requires_subcommand():
    """Running with no subcommand is a usage error (exit code 2)."""
    result = subprocess.run(["chaste-sbml"], capture_output=True, text=True)

    assert result.returncode == 2


def test_generate_requires_sbml_file():
    """generate without an SBML file is a usage error (exit code 2)."""
    result = subprocess.run(["chaste-sbml", "generate"], capture_output=True, text=True)

    assert result.returncode == 2


def test_copy_base_classes_rejects_generation_options():
    """copy-base-classes takes only --output-dir; generation options are usage errors (exit 2)."""
    for extra in (["model.xml"], ["--model-type", "srn"]):
        result = subprocess.run(
            ["chaste-sbml", "copy-base-classes", *extra],
            capture_output=True,
            text=True,
        )

        assert result.returncode == 2, extra


def test_copy_base_classes_accepts_output_dir(tmp_path):
    """copy-base-classes with --output-dir copies the base classes."""
    result = subprocess.run(
        ["chaste-sbml", "copy-base-classes", "--output-dir", str(tmp_path)],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr
    assert (tmp_path / "AbstractSbmlOdeSystem.hpp").is_file()
