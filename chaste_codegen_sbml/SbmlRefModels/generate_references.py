"""Regenerate the Chaste reference models from their SBML sources.

Each reference model lives in its own directory under the reference root and contains a
single SBML file (``<Name>.xml``) alongside the generated Chaste code. The model name is
``<Name>Sbml`` and the model type is inferred from the companion files already present:
an ``*SrnModel.*`` makes it an SRN model, a ``*CellCycleModel.*`` a cell-cycle model,
otherwise a generic ODE system.
"""

import argparse
import glob
import logging
import os
import sys

from chaste_codegen_sbml import ChasteSbmlModel
from chaste_codegen_sbml._config import ModelType

logger = logging.getLogger(__name__)


def _infer_model_type(model_dir: str) -> ModelType:
    """Infer a reference model's type from the generated files already present.

    :param model_dir: The reference model directory.
    :return: The inferred model type.
    """
    if glob.glob(os.path.join(model_dir, "*SrnModel.*")):
        return ModelType.SRN
    if glob.glob(os.path.join(model_dir, "*CellCycleModel.*")):
        return ModelType.CELL_CYCLE
    return ModelType.GENERIC


def generate_references(reference_dir: str) -> None:
    """Regenerate every reference model found under the reference directory.

    :param reference_dir: Directory containing one subdirectory per reference model.
    """
    for entry in sorted(os.listdir(reference_dir)):
        model_dir = os.path.join(reference_dir, entry)
        if not os.path.isdir(model_dir):
            continue

        sbml_files = glob.glob(os.path.join(model_dir, "*.xml"))
        if len(sbml_files) != 1:
            logger.warning(f"Skipping '{entry}': expected exactly one SBML file, found {len(sbml_files)}.")
            continue

        sbml_file = sbml_files[0]
        model_name = os.path.splitext(os.path.basename(sbml_file))[0] + "Sbml"
        model_type = _infer_model_type(model_dir)

        logger.info(f"Regenerating {model_name} ({model_type.name}) from {os.path.basename(sbml_file)}")
        model = ChasteSbmlModel(sbml_file, model_name=model_name, model_type=model_type)
        model.write(model_dir)


def parse_args() -> argparse.Namespace:
    """Parse command line arguments.

    :return: The parsed command line arguments.
    """
    parser = argparse.ArgumentParser(
        prog="generate_references",
        description="Regenerate the Chaste reference models from their SBML sources.",
    )
    parser.add_argument(
        "--reference-dir",
        type=str,
        default="src/reference",
        help="Directory containing one subdirectory per reference model.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    logging.basicConfig(level=logging.INFO)

    args.reference_dir = os.path.abspath(args.reference_dir)
    if not os.path.isdir(args.reference_dir):
        logger.error(f"No reference directory @ '{args.reference_dir}'")
        sys.exit(1)

    generate_references(args.reference_dir)


if __name__ == "__main__":
    main()
