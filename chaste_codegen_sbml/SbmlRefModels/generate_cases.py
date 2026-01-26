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

sbml_versions = ["l2v5", "l3v2", "l3v1"]


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

        if test_params["settings"]["start"].strip() == "":
            raise NotImplementedError("Steady state test with no start time.")

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
        test_concentrations = ", ".join(
            f'"{conc.strip()}"' for conc in test_concentrations if conc.strip()
        )

        self._template_vars.update(
            {
                "test_header_guard": generate_header_guard(self._test_hpp_filename),
                "test_result_columns": test_result_columns,
                "test_result_data": test_result_data,
                "test_amounts": test_amounts,
                "test_concentrations": test_concentrations,
                "test_settings": test_params["settings"],
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
    cases = sorted(os.listdir(semantic_path))

    test_pack = []
    selection = [1103, 1106, 1117, 1118, 1121, 1122, 1123, 1182, 1183, 1184, 1185, 1198]
    cases = [cases[i - 1] for i in selection]
    for case_ in cases:
        case_path = os.path.join(semantic_path, case_)
        if not os.path.isdir(case_path):
            continue

        test_params = {
            "type": TestType.SEMANTIC,
            "case": case_,
        }

        # Load settings from settings file
        settings_file = os.path.join(case_path, f"{case_}-settings.txt")
        if not os.path.isfile(settings_file):
            logger.warning(f"Settings not found for semantic {case_}: '{settings_file}'")
            continue

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
        model_dir = ROOT_DIR / "SbmlRefModels" / "src" / "cases" / "semantic" / case_
        if not os.path.exists(model_dir):
            os.makedirs(model_dir)

        test_dir = ROOT_DIR / "SbmlRefModels" / "test" / "cases" / "semantic" / case_
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
    with open(ROOT_DIR / "S"
    "bmlRefModels" / "test" / "WeeklyTestPack.txt", "w") as f:
        f.write("\n".join(test_pack) + "\n")


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
