"""End-to-end tests for the ChasteSbmlModel orchestrator."""

from chaste_sbml import ChasteSbmlModel
from chaste_sbml._config import ROOT_DIR, ModelType

GOLDBETER = ROOT_DIR / "SbmlRefModels" / "src" / "reference" / "Goldbeter1991" / "Goldbeter1991.xml"


def test_write_produces_output_files(tmp_path):
    """Writing an SRN model renders and writes its OdeSystem and SRN sources."""
    model = ChasteSbmlModel(str(GOLDBETER), model_type=ModelType.SRN)
    model.write(tmp_path)

    written = {p.name for p in tmp_path.iterdir()}
    assert "Goldbeter1991SbmlOdeSystem.hpp" in written
    assert "Goldbeter1991SbmlOdeSystem.cpp" in written
    assert "Goldbeter1991SbmlSrnModel.hpp" in written
    assert "Goldbeter1991SbmlSrnModel.cpp" in written


def test_outputs_property_exposes_generated_code(tmp_path):
    """The outputs property returns the generated filename -> code mapping after write."""
    model = ChasteSbmlModel(str(GOLDBETER), model_type=ModelType.SRN)
    model.write(tmp_path)

    outputs = model.outputs
    assert isinstance(outputs, dict)
    assert "Goldbeter1991SbmlOdeSystem.cpp" in outputs
    assert outputs["Goldbeter1991SbmlOdeSystem.cpp"].strip()  # non-empty source
