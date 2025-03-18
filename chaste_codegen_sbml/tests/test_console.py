import logging
import subprocess

from chaste_codegen_sbml import __version__

logger = logging.getLogger(__name__)


def test_help():
    """Test help message."""

    expected = """usage: chaste_codegen_sbml [-h] [--version]
                  [--output-dir OUTPUT_DIR]
                  [--srn] [--cell-cycle]
                  sbml_file"""
    expected = " ".join(expected.split())

    output = subprocess.check_output(["chaste_codegen_sbml", "-h"]).decode("ascii")
    output = " ".join(output.split())

    assert output.startswith(expected)


def test_version():
    """Test version message."""

    expected = f"chaste_codegen_sbml {__version__}"

    output = subprocess.check_output(["chaste_codegen_sbml", "--version"]).decode("ascii").strip()

    assert output == expected
