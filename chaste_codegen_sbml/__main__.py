"""Entry point for the command line interface."""

import os
from argparse import ArgumentParser

from . import ChasteSbmlSrnModel


def parse_args():
    """Parse command line arguments."""
    parser = ArgumentParser(
        prog="chaste_codegen_sbml",
        description="Generate C++ code from SBML models for the Chaste C++ library",
    )

    parser.add_argument(
        "sbml_file",
        metavar="sbml_file",
        help="The sbml file to convert to Chaste C++ code",
    )

    args = parser.parse_args()

    if not os.path.isfile(args.sbml_file):
        raise FileNotFoundError(f"Could not find sbml file {args.sbml_file}")

    return args


def generate_code(args):
    """Run the code generation."""
    chaste_model = ChasteSbmlSrnModel(args.sbml_file)
    chaste_model.write_chaste_code()


def main():
    """Run the command line interface."""
    args = parse_args()
    generate_code(args)
