"""Entry point for the command line interface."""

import argparse

from chaste_codegen_sbml import ChasteSbmlCellCycleModel, ChasteSbmlSrnModel

from ._version import __version__


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        prog="chaste_codegen_sbml",
        description="Convert SBML models to C++ code for the Chaste library",
    )

    parser.add_argument(
        "sbml_file",
        help="The SBML file to convert",
    )

    parser.add_argument("--version", action="version", version="%(prog)s " + __version__)

    parser.add_argument(
        "--output-dir", action="store", help="The directory to place output files in", default=None
    )

    group = parser.add_argument_group(
        "ModelTypes", 'The types of models code can be generated for; defaults to "srn"'
    )
    group.add_argument("--srn", help="Generate SRN model", action="store_true")
    group.add_argument("--cell-cycle", help="Generate Cell Cycle model", action="store_true")

    args = parser.parse_args()

    return args


def process_command_line(args: "argparse.Namespace"):
    """Run the command line interface."""

    if args.cell_cycle:
        chaste_model = ChasteSbmlCellCycleModel(args.sbml_file)
    else:
        chaste_model = ChasteSbmlSrnModel(args.sbml_file)

    chaste_model.write(args.output_dir)


def main():
    """Main entrypoint."""
    args = parse_args()
    process_command_line(args)
