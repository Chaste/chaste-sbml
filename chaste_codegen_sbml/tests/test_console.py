import logging
import subprocess

logger = logging.getLogger(__name__)


def test_console_help(capsys):
    """Test help message"""

    expected = """
usage: chaste_codegen_sbml [-h] [--output-dir OUTPUT_DIR] [--srn] [--cell-cycle] sbml_file

Generate C++ code from SBML models for the Chaste C++ library

positional arguments:
  sbml_file             The sbml file to convert to Chaste C++ code

options:
  -h, --help            show this help message and exit
  --output-dir OUTPUT_DIR
                        The directory to place output files in

ModelTypes:
  The types of models code can be generated for; defaults to "srn"

  --srn                 Generate SRN model
  --cell-cycle          Generate Cell Cycle model
    """
    expected = " ".join(expected.split())

    out = subprocess.check_output(["chaste_codegen_sbml", "-h"]).decode("ascii")
    out = " ".join(out.split())

    assert out == expected
