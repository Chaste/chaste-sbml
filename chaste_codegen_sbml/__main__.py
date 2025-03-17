"""Entry point for the command line interface."""

from argparse import ArgumentParser

from chaste_codegen_sbml import ChasteSbmlCellCycleModel, ChasteSbmlSrnModel


def parse_args():
    """Parse command line arguments."""
    parser = ArgumentParser(
        prog="chaste_codegen_sbml",
        description="Generate C++ code from SBML models for the Chaste C++ library",
    )

    parser.add_argument(
        "sbml_file",
        help="The sbml file to convert to Chaste C++ code",
    )

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


def generate_code(args):
    """Run the code generation."""

    if args.cell_cycle:
        chaste_model = ChasteSbmlCellCycleModel(args.sbml_file)
    else:
        chaste_model = ChasteSbmlSrnModel(args.sbml_file)

    chaste_model.write(args.output_dir)


def main():
    """Run the command line interface."""
    args = parse_args()
    generate_code(args)
