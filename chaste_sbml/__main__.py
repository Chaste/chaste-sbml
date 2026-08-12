"""Entry point for the command line interface."""

import argparse

from chaste_sbml import ChasteSbmlModel, copy_base_classes
from chaste_sbml._config import ModelType, TimeUnit

from ._version import __version__

# Options that apply only when generating, mapped to the default each takes. In copy mode none may
# be supplied (only --output-dir applies), so they use ``argparse.SUPPRESS`` as their default: the
# attribute is absent unless the option was passed, which lets copy mode reject an explicitly-supplied
# option while generation fills in the real default below.
_GENERATION_DEFAULTS = {
    "model_type": "generic",
    "tests": True,
    "timescale": None,
    "test_output_dir": None,
}


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
        help="The type of model to generate (default: generic)",
        choices=["generic", "srn", "cell-cycle"],
        default=argparse.SUPPRESS,
        const="generic",
        nargs="?",
    )
    parser.add_argument(
        "--tests",
        action=argparse.BooleanOptionalAction,
        default=argparse.SUPPRESS,
        help="Generate placeholder test files (default: on). Use --no-tests to disable.",
    )
    parser.add_argument(
        "--timescale",
        choices=["ms", "s", "m", "h"],
        default=argparse.SUPPRESS,
        help="The model's native time unit (milliseconds/seconds/minutes/hours), used to convert "
        "derivatives to Chaste's hours. Overrides auto-detection from the SBML; omit to auto-detect.",
    )
    parser.add_argument(
        "--test-output-dir",
        default=argparse.SUPPRESS,
        help="The directory to place generated test files in (defaults to --output-dir)",
    )

    args = parser.parse_args()

    # Options are stored only when supplied (default=SUPPRESS), so a missing attribute means unset.
    supplied_generation_opts = [name for name in _GENERATION_DEFAULTS if hasattr(args, name)]

    if args.copy_base_classes:
        # Copy mode takes only --output-dir; reject the SBML file and any generation-only option
        # rather than silently ignoring them.
        rejected = ["an SBML file"] if args.sbml_file is not None else []
        rejected += [f"--{name.replace('_', '-')}" for name in supplied_generation_opts]
        if rejected:
            parser.error("--copy-base-classes takes only --output-dir, not: " + ", ".join(rejected))
    else:
        if args.sbml_file is None:
            parser.error("an SBML file is required (or pass --copy-base-classes to copy the base classes)")
        # Fill in the real default for every generation option that was not supplied.
        for name, default in _GENERATION_DEFAULTS.items():
            if not hasattr(args, name):
                setattr(args, name, default)

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
