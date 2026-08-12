"""Entry point for the command line interface."""

import argparse

from chaste_sbml import ChasteSbmlModel, copy_base_classes
from chaste_sbml._config import ModelType, TimeUnit

from ._version import __version__


def parse_args() -> argparse.Namespace:
    """Parse command line arguments.

    Generating Chaste C++ code from an SBML file is the default action. Passing ``--copy-base-classes``
    switches to copying the C++ base classes the generated code depends on instead; in that mode no
    SBML file is taken and the only other option that applies is ``--output-dir``.
    """
    parser = argparse.ArgumentParser(
        prog="chaste-sbml",
        description="Convert SBML models to Chaste C++ code",
    )
    parser.add_argument("--version", action="version", version="%(prog)s " + __version__)

    parser.add_argument("sbml_file", nargs="?", help="The SBML file to convert")
    parser.add_argument(
        "--copy-base-classes",
        action="store_true",
        help="Copy the C++ base classes the generated code depends on, instead of generating code",
    )
    parser.add_argument("--output-dir", default=None, help="The directory to place output files in")
    parser.add_argument(
        "--model-type",
        help="The type of model to generate",
        choices=["generic", "srn", "cell-cycle"],
        default="generic",
        const="generic",
        nargs="?",
    )
    parser.add_argument(
        "--tests",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Generate placeholder test files (default: on). Use --no-tests to disable.",
    )
    parser.add_argument(
        "--timescale",
        choices=["ms", "s", "m", "h"],
        default=None,
        help="The model's native time unit (milliseconds/seconds/minutes/hours), used to convert "
        "derivatives to Chaste's hours. Overrides auto-detection from the SBML; omit to auto-detect.",
    )
    parser.add_argument(
        "--test-output-dir",
        default=None,
        help="The directory to place generated test files in (defaults to --output-dir)",
    )

    args = parser.parse_args()

    if args.copy_base_classes:
        if args.sbml_file is not None:
            parser.error("--copy-base-classes does not take an SBML file")
    elif args.sbml_file is None:
        parser.error("an SBML file is required (or pass --copy-base-classes to copy the base classes)")

    return args


def process_command_line(args: "argparse.Namespace"):
    """Run the command line interface.

    :args: The parsed command line arguments.
    """
    if args.copy_base_classes:
        copy_base_classes(args.output_dir)
        return

    model_type = ModelType.GENERIC
    if args.model_type == "srn":
        model_type = ModelType.SRN
    elif args.model_type == "cell-cycle":
        model_type = ModelType.CELL_CYCLE

    time_unit = TimeUnit.from_cli(args.timescale) if args.timescale else None

    chaste_model = ChasteSbmlModel(
        args.sbml_file, model_type=model_type, generate_tests=args.tests, time_unit=time_unit
    )
    chaste_model.write(args.output_dir, args.test_output_dir)


def main():
    """Run the command line interface."""
    args = parse_args()
    process_command_line(args)


if __name__ == "__main__":  # pragma: no cover
    main()
