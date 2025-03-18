import logging
import subprocess

logger = logging.getLogger(__name__)


def test_console_help():
    """Test help message"""

    expected = """usage: chaste_codegen_sbml [-h]
                  [--output-dir OUTPUT_DIR]
                  [--srn] [--cell-cycle]
                  sbml_file"""
    expected = " ".join(expected.split())

    out = subprocess.check_output(["chaste_codegen_sbml", "-h"]).decode("ascii")
    out = " ".join(out.split())

    assert out.startswith(expected)
