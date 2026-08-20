"""Tests for placeholder test-file generation."""

import logging

import pytest

from chaste_sbml import ChasteSbmlModel
from chaste_sbml._config import ROOT_DIR, ModelType

logger = logging.getLogger(__name__)

REFERENCE_DIR = ROOT_DIR / "SbmlRefModels" / "src" / "reference"


def sbml_file(model_name: str) -> str:
    """Path to a reference model's SBML input.

    :param model_name: The reference model name e.g. "Goldbeter1991".
    :return: Path to the SBML file.
    """
    return str(REFERENCE_DIR / model_name / f"{model_name}.xml")


@pytest.mark.parametrize(
    ("model_name", "model_type", "extra_class", "extra_method"),
    [
        ("Chen2000", ModelType.GENERIC, None, None),
        ("Goldbeter1991", ModelType.SRN, "TestGoldbeter1991SbmlSrnModel", "TestModel"),
        ("Chen2000", ModelType.CELL_CYCLE, "TestChen2000SbmlCellCycleModel", "TestModel"),
    ],
)
def test_placeholder_test_generated(tmp_path, model_name, model_type, extra_class, extra_method):
    """A placeholder test file is generated with an OdeSystem suite, plus an SRN/CellCycle suite."""
    model = ChasteSbmlModel(sbml_file(model_name), model_type=model_type)
    model.write(tmp_path)

    test_filename = f"Test{model_name}Sbml.hpp"
    assert test_filename in model.test_outputs
    assert (tmp_path / test_filename).is_file()

    code = model.test_outputs[test_filename]

    # Every Chaste test must initialise PETSc; the placeholder is never run in parallel.
    assert '#include "FakePetscSetup.hpp"' in code

    # The OdeSystem suite is always present and uses the lightweight CxxTest base class.
    assert f"class Test{model_name}SbmlOdeSystem : public CxxTest::TestSuite" in code
    assert "void TestOdeSystem()" in code
    assert f"TS_ASSERT_THROWS_NOTHING({model_name}SbmlOdeSystem ode_system);" in code

    if extra_class is None:
        # Generic models get only the OdeSystem suite -- no cell-based dependency.
        assert "AbstractCellBasedTestSuite" not in code
        assert "SrnModel" not in code
        assert "CellCycleModel" not in code
    else:
        # SRN/CellCycle models get a separate cell-based suite for the extra role.
        assert f"class {extra_class} : public AbstractCellBasedTestSuite" in code
        assert f"void {extra_method}()" in code


def test_no_tests_disables_generation(tmp_path):
    """generate_tests=False produces no placeholder test, on disk or in memory."""
    model = ChasteSbmlModel(sbml_file("Goldbeter1991"), model_type=ModelType.SRN, generate_tests=False)
    model.write(tmp_path)

    assert model.test_outputs == {}
    written = {p.name for p in tmp_path.iterdir()}
    assert "TestGoldbeter1991Sbml.hpp" not in written
    # The model code is still generated.
    assert "Goldbeter1991SbmlOdeSystem.hpp" in written


def test_test_output_directory_routing(tmp_path):
    """The placeholder test goes to test_output_directory, model code to output_directory."""
    src_dir = tmp_path / "src"
    test_dir = tmp_path / "test"
    src_dir.mkdir()
    test_dir.mkdir()

    model = ChasteSbmlModel(sbml_file("Goldbeter1991"), model_type=ModelType.SRN)
    model.write(src_dir, test_dir)

    assert (test_dir / "TestGoldbeter1991Sbml.hpp").is_file()
    assert not (src_dir / "TestGoldbeter1991Sbml.hpp").exists()
    assert (src_dir / "Goldbeter1991SbmlOdeSystem.hpp").is_file()


def test_test_defaults_to_output_directory(tmp_path):
    """Without a test directory, the placeholder test lands next to the model code."""
    model = ChasteSbmlModel(sbml_file("Goldbeter1991"), model_type=ModelType.SRN)
    model.write(tmp_path)

    written = {p.name for p in tmp_path.iterdir()}
    assert "TestGoldbeter1991Sbml.hpp" in written
    assert "Goldbeter1991SbmlOdeSystem.hpp" in written
