"""Entry point for the command line interface."""

import argparse

from chaste_sbml import ChasteSbmlModel, copy_base_classes
from chaste_sbml._config import ModelType

from ._version import __version__


def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        prog="chaste-sbml",
        description="Convert SBML models to Chaste C++ code",
    )
    parser.add_argument("--version", action="version", version="%(prog)s " + __version__)

    subparsers = parser.add_subparsers(dest="command", required=True, metavar="command")

    # generate: create Chaste C++ code from an SBML file.
    generate = subparsers.add_parser(
        "generate",
        help="Generate Chaste C++ code from an SBML file",
        description="Generate Chaste C++ code from an SBML file",
    )
    generate.add_argument("sbml_file", help="The SBML file to convert")
    generate.add_argument("--output-dir", default=None, help="The directory to place output files in")
    generate.add_argument(
        "--model-type",
        help="The type of model to generate",
        choices=["generic", "srn", "cell-cycle"],
        default="generic",
        const="generic",
        nargs="?",
    )

    # copy-base-classes: copy the C++ base classes the generated code depends on.
    copy_parser = subparsers.add_parser(
        "copy-base-classes",
        help="Copy the C++ base classes the generated code depends on",
        description="Copy the C++ base classes the generated code depends on",
    )
    copy_parser.add_argument("--output-dir", default=None, help="The directory to place the base classes in")

    return parser.parse_args()


def process_command_line(args: "argparse.Namespace"):
    """Run the command line interface.

    :args: The parsed command line arguments.
    """
    if args.command == "copy-base-classes":
        copy_base_classes(args.output_dir)
        return

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


if __name__ == "__main__":
    main()
