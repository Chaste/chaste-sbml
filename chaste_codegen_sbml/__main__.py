import os
from argparse import ArgumentParser

from chaste_codegen_sbml.chaste_model import Generate


def parse_args():
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
    Generate(args.sbml_file)


def main():
    args = parse_args()
    generate_code(args)
