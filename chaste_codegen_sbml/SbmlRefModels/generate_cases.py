"""Generate Chaste code from the cases in the SBML Test Suite."""

import argparse
import csv
import logging
import os
import sys
from enum import Enum
from typing import TYPE_CHECKING

from chaste_codegen_sbml import ChasteSbmlModel
from chaste_codegen_sbml._config import ROOT_DIR, ModelType
from chaste_codegen_sbml._utils import generate_header_guard

if TYPE_CHECKING:
    from typing import Any

logger = logging.getLogger(__name__)

sbml_versions = ["l2v5", "l3v2", "l3v1"]


class TestType(Enum):
    """Enumeration of test types for code generation."""

    SEMANTIC = 0
    STOCHASTIC = 1
    SYNTACTIC = 2


class ChasteSbmlTestSuiteModel(ChasteSbmlModel):
    """Extended ChasteSbmlModel with additional functionality for test suite generation."""

    def __init__(self, sbml_file: str, sbml_version: str, test_params: dict[str, "Any"], **kwargs) -> None:
        test_type = test_params["type"]
        if test_type == TestType.SEMANTIC:
            prefix = "Semantic"
        elif test_type == TestType.STOCHASTIC:
            prefix = "Stochastic"
        elif test_type == TestType.SYNTACTIC:
            prefix = "Syntactic"
        else:
            raise ValueError(f"Unknown test type: {test_type}")

        model_name = f"{prefix}{test_params['case']}{sbml_version.upper()}Sbml"

        super().__init__(sbml_file, model_name=model_name, model_type=ModelType.GENERIC, **kwargs)

        self._test_type = test_type
        self._test_hpp_filename = f"Test{self._model_name}.hpp"

        test_result_columns = ", ".join(f'"{col}"' for col in test_params["results"][0])
        test_result_data = ["{ " + ", ".join(row) + " }" for row in test_params["results"][1:]]
        test_result_data = ",\n".join(test_result_data)

        test_amounts = test_params["settings"]["amount"].split(",")
        test_amounts = ", ".join(f'"{amt.strip()}"' for amt in test_amounts if amt.strip())

        test_concentrations = test_params["settings"]["concentration"].split(",")
        test_concentrations = ", ".join(f'"{conc.strip()}"' for conc in test_concentrations if conc.strip())

        steady_state = test_params["settings"]["start"].strip() == ""

        self._template_vars.update(
            {
                "test_header_guard": generate_header_guard(self._test_hpp_filename),
                "test_result_columns": test_result_columns,
                "test_result_data": test_result_data,
                "test_amounts": test_amounts,
                "test_concentrations": test_concentrations,
                "test_settings": test_params["settings"],
                "test_steady_state": steady_state,
            }
        )

    def _generate_outputs(self) -> None:
        """Generate Chaste code for the model and tests."""
        # Generate the model code
        super()._generate_outputs()

        # Generate the test code
        template_path = "cases"
        if self._test_type == TestType.SEMANTIC:
            template_path += "/semantic.hpp"
        elif self._test_type == TestType.STOCHASTIC:
            template_path += "/stochastic.hpp"
        elif self._test_type == TestType.SYNTACTIC:
            template_path += "/syntactic.hpp"
        else:
            raise ValueError(f"Unknown test type: {self._test_type}")

        self._generate_output(template_path, self._test_hpp_filename)


def generate_semantic_cases(
    selection: list[int] | None = None,
    sbml_test_suite_dir: str = "sbml_test_suite",
    gen_src_dir: str = "src",
    gen_test_dir: str = "test",
    test_pack_file: str = "test/WeeklyTestPack.txt",
) -> None:
    """Generate Chaste code from sbml_test_suite semantic cases.

    :param selection: List of case numbers to generate. If None, all cases will be generated.
    :param sbml_test_suite_dir: Path to the SBML test suite repository.
    :param gen_src_dir: Path to the directory where generated source code will be placed.
    :param gen_test_dir: Path to the directory where generated tests will be placed.
    :param test_pack_file: Path to the file where generated test paths will be stored.
    """
    semantic_path = os.path.join(sbml_test_suite_dir, "cases", "semantic")
    cases = sorted(os.listdir(semantic_path))

    test_pack = []
    if selection:
        cases = [cases[i - 1] for i in selection]

    for case_ in cases:
        case_path = os.path.join(semantic_path, case_)
        if not os.path.isdir(case_path):
            # Skip non-directory entries (e.g. README files)
            continue

        test_params = {
            "type": TestType.SEMANTIC,
            "case": case_,
        }

        # Load settings from settings file
        settings_file = os.path.join(case_path, f"{case_}-settings.txt")
        if not os.path.isfile(settings_file):
            raise ValueError(f"Settings not found for semantic {case_}: '{settings_file}'")

        settings = {}
        with open(settings_file, "r") as f:
            for line in f:
                if ":" in line:
                    key, value = line.strip().split(":")
                    settings[key.strip()] = value.strip()

        test_params["settings"] = settings

        # Load results from results file
        results_file = os.path.join(case_path, f"{case_}-results.csv")
        if not os.path.isfile(results_file):
            logger.warning(f"Results not found for semantic {case_}: '{results_file}'")
            continue

        with open(results_file, "r") as csvfile:
            reader = csv.reader(csvfile)
            test_params["results"] = list(reader)

        # Make output directories
        model_dir = os.path.join(gen_src_dir, case_)
        if not os.path.exists(model_dir):
            os.makedirs(model_dir)

        test_dir = os.path.join(gen_test_dir, case_)
        if not os.path.exists(test_dir):
            os.makedirs(test_dir)

        # Generate code for models and tests
        for sbml_version in sbml_versions:
            sbml_file = os.path.join(case_path, f"{case_}-sbml-{sbml_version}.xml")
            if not os.path.isfile(sbml_file):
                continue

            logger.info(f"Converting semantic {case_} {sbml_version}")

            try:
                chaste_model = ChasteSbmlTestSuiteModel(sbml_file, sbml_version, test_params)
                chaste_model.write()
                for src_file, _ in chaste_model.outputs.items():
                    if src_file.startswith("Test"):
                        dst_dir = test_dir
                        test_pack.append(f"cases/semantic/{case_}/{src_file}")
                    else:
                        dst_dir = model_dir
                    dst_file = os.path.join(dst_dir, src_file)
                    os.replace(src_file, dst_file)
            except NotImplementedError as e:
                logger.warning(f"Skipping semantic {case_} {sbml_version}: {e}")
            break  # Only generate for the first available SBML version

    # Update WeeklyTestPack
    logger.info(f"Updating {test_pack_file}")
    with open(test_pack_file, "w") as f:
        f.write("\n".join(test_pack) + "\n")


def parse_args() -> argparse.Namespace:
    """Parse command line arguments.

    :return: The parsed command line arguments.
    """
    parser = argparse.ArgumentParser(
        prog="generate_cases",
        description="Generate Chaste code from the cases in the SBML Test Suite.",
    )

    parser.add_argument(
        "--sbml-test-suite-dir",
        type=str,
        default="sbml-test-suite",
        help="Directory containing the SBML test suite.",
    )

    parser.add_argument(
        "--gen-src-dir",
        type=str,
        default="src",
        help="Directory where generated code will be placed.",
    )

    parser.add_argument(
        "--gen-test-dir",
        type=str,
        default="test",
        help="Directory where generated tests will be placed.",
    )

    parser.add_argument(
        "--test-pack-file",
        type=str,
        default="test/WeeklyTestPack.txt",
        help="File where generated test paths will be stored.",
    )

    parser.add_argument(
        "--first-case",
        type=int,
        default=0,
        help="First test case to generate.",
    )

    parser.add_argument(
        "--last-case",
        type=int,
        default=100,
        help="Last test case to generate.",
    )

    return parser.parse_args()


def main() -> None:
    args = parse_args()

    logging.basicConfig(level=logging.INFO)

    args.sbml_test_suite_dir = os.path.abspath(args.sbml_test_suite_dir)
    if not os.path.isdir(args.sbml_test_suite_dir):
        logger.error(f"No SBML test suite @ '{args.sbml_test_suite_dir}'")
        sys.exit(1)

    args.test_pack_file = os.path.abspath(args.test_pack_file)
    if not os.path.isfile(args.test_pack_file):
        logger.error(f"No test pack file @ '{args.test_pack_file}'")
        sys.exit(1)

    # TODO: Only semantic cases are currently supported
    generate_semantic_cases(
        selection=list(range(args.first_case, args.last_case + 1)),
        sbml_test_suite_dir=args.sbml_test_suite_dir,
        gen_src_dir=args.gen_src_dir,
        gen_test_dir=args.gen_test_dir,
        test_pack_file=args.test_pack_file,
    )


if __name__ == "__main__":
    main()
