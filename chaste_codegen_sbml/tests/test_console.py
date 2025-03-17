import logging
import subprocess

logger = logging.getLogger(__name__)


def test_console_help(capsys):
    """Test help message"""
    testargs = ["chaste_codegen_sbml", "-h"]

    expected_output = """
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
    expected_output = " ".join(expected_output.split())

    output = subprocess.check_output(testargs).decode("ascii")
    output = " ".join(output.split())

    assert output == expected_output
