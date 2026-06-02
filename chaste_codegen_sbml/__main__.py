"""Entry point for the command line interface."""

import argparse

from chaste_codegen_sbml import ChasteSbmlModel
from chaste_codegen_sbml._config import ModelType

from ._version import __version__


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        prog="chaste_codegen_sbml",
        description="Convert SBML models to Chaste C++ code",
    )

    parser.add_argument(
        "sbml_file",
        help="The SBML file to convert",
    )

    parser.add_argument("--version", action="version", version="%(prog)s " + __version__)

    parser.add_argument("--output-dir", action="store", help="The directory to place output files in", default=None)

    parser.add_argument(
        "--model-type",
        help="The type of model to generate",
        choices=["generic", "srn", "cell-cycle"],
        default="generic",
        const="generic",
        nargs="?",
    )

    args = parser.parse_args()

    return args


def process_command_line(args: "argparse.Namespace"):
    """Run the command line interface.

    :args: The parsed command line arguments.
    """
    model_type = ModelType.GENERIC
    if args.model_type == "srn":
        model_type = ModelType.SRN
    elif args.model_type == "cell-cycle":
        model_type = ModelType.CELL_CYCLE

    chaste_model = ChasteSbmlModel(args.sbml_file, model_type=model_type)
    chaste_model.write(args.output_dir)


def main():
    """Run the command line interface."""
    args = parse_args()
    process_command_line(args)
