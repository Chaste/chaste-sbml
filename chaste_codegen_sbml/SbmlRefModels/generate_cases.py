"""Generate Chaste code from the cases in the SBML Test Suite."""

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

sbml_versions = ["l2v5"]


class TestType(Enum):
    """Enumeration of test types for code generation."""

    SEMANTIC = 0
    STOCHASTIC = 1
    SYNTACTIC = 2


class ChasteSbmlTestSuiteModel(ChasteSbmlModel):
    """Extended ChasteSbmlModel with additional functionality for test suite generation."""

    def __init__(
        self, sbml_file: str, sbml_version: str, test_params: dict[str, "Any"], **kwargs
    ) -> None:
        test_type = test_params["type"]
        match test_type:
            case TestType.SEMANTIC:
                prefix = "Semantic"
            case TestType.STOCHASTIC:
                prefix = "Stochastic"
            case TestType.SYNTACTIC:
                prefix = "Syntactic"
            case _:
                raise ValueError(f"Unknown test type: {test_type}")

        model_name = f"{prefix}{test_params['case']}{sbml_version.upper()}Sbml"

        super().__init__(sbml_file, model_name=model_name, model_type=ModelType.GENERIC, **kwargs)

        self._test_type = test_type
        self._test_settings = test_params["settings"]

        test_results = ["{ " + ", ".join(result) + " }" for result in test_params["results"][1:]]
        self._test_results = ",\n".join(test_results)

        self._test_hpp_filename = f"Test{self._model_name}.hpp"

        self._template_vars.update(
            {
                "test_header_guard": generate_header_guard(self._test_hpp_filename),
                "test_results": self._test_results,
                "test_settings": self._test_settings,
            }
        )

    def _generate_outputs(self) -> None:
        """Generate Chaste code for the model and tests."""
        # Generate the model code
        super()._generate_outputs()

        # Generate the test code
        template_path = "cases"
        match self._test_type:
            case TestType.SEMANTIC:
                template_path += "/semantic.hpp"
            case TestType.STOCHASTIC:
                template_path += "/stochastic.hpp"
            case TestType.SYNTACTIC:
                template_path += "/syntactic.hpp"
            case _:
                raise ValueError(f"Unknown test type: {self._test_type}")

        self._generate_output(template_path, self._test_hpp_filename)


def generate_semantic_cases(suite_path: str) -> None:
    """Generate Chaste code from sbml_test_suite semantic cases.

    :param suite_path: Path to the sbml_test_suite repository.
    """
    semantic_path = os.path.join(suite_path, "cases", "semantic")

    for case in sorted(os.listdir(semantic_path))[1538:]:
        case_path = os.path.join(semantic_path, case)
        if not os.path.isdir(case_path):
            continue

        test_params = {
            "type": TestType.SEMANTIC,
            "case": case,
        }

        # Load settings from settings file
        settings_file = os.path.join(case_path, f"{case}-settings.txt")
        if not os.path.isfile(settings_file):
            logger.warning(f"Settings not found for semantic {case}: '{settings_file}'")
            continue

        settings = {}
        with open(settings_file, "r") as f:
            for line in f:
                if ":" in line:
                    key, value = line.strip().split(":")
                    settings[key.strip()] = value.strip()

        test_params["settings"] = settings

        # Load results from results file
        results_file = os.path.join(case_path, f"{case}-results.csv")
        if not os.path.isfile(results_file):
            logger.warning(f"Results not found for semantic {case}: '{results_file}'")
            continue

        with open(results_file, "r") as csvfile:
            reader = csv.reader(csvfile)
            test_params["results"] = list(reader)

        # Make output directories
        model_dir = ROOT_DIR / "SbmlRefModels" / "src" / "cases" / "semantic" / case
        if not os.path.exists(model_dir):
            os.makedirs(model_dir)

        test_dir = ROOT_DIR / "SbmlRefModels" / "test" / "cases" / "semantic" / case
        if not os.path.exists(test_dir):
            os.makedirs(test_dir)

        # Generate code for models and tests
        for sbml_version in sbml_versions:
            sbml_file = os.path.join(case_path, f"{case}-sbml-{sbml_version}.xml")
            if not os.path.isfile(sbml_file):
                continue

            logger.info(f"Converting semantic {case} {sbml_version}")

            try:
                chaste_model = ChasteSbmlTestSuiteModel(sbml_file, sbml_version, test_params)
                chaste_model.write()
                for src_file, _ in chaste_model.outputs.items():
                    if src_file.startswith("Test"):
                        dst_dir = test_dir
                    else:
                        dst_dir = model_dir
                    dst_file = os.path.join(dst_dir, src_file)
                    os.replace(src_file, dst_file)
            except NotImplementedError as e:
                logger.warning(f"Skipping semantic {case} {sbml_version}: {e}")


def generate_stochastic_cases(suite_path: str) -> None:
    """Generate Chaste code from sbml_test_suite stochastic cases.

    :param suite_path: Path to the sbml_test_suite repository.
    """
    pass


def generate_syntactic_cases(suite_path: str) -> None:
    """Generate Chaste code from sbml_test_suite syntactic cases.

    :param suite_path: Path to the sbml_test_suite repository.
    """
    pass


def parse_args() -> str:
    """Parse command line arguments.

    :return: Path to the sbml_test_suite repository.
    """
    if len(sys.argv) != 2:
        print(f"Usage: python {sys.argv[0]} /path/to/sbml_test_suite")
        sys.exit(1)

    suite_path = sys.argv[1]
    if not os.path.isdir(suite_path):
        print(f"Error: Path '{suite_path}' does not exist.")
        sys.exit(1)

    return suite_path


def main() -> None:
    logging.basicConfig(level=logging.INFO)

    suite_path = parse_args()

    generate_semantic_cases(suite_path)
    # generate_stochastic_cases(suite_path)
    # generate_syntactic_cases(suite_path)


if __name__ == "__main__":
    main()
