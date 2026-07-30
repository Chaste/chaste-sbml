"""In-process tests for the CLI entry point.

These call ``chaste_sbml.__main__.main()`` directly (with a patched ``sys.argv``) so the CLI code is
measured by in-process coverage, unlike the subprocess smoke tests in ``test_console.py``.
"""

import pytest

import chaste_sbml.__main__ as cli
from chaste_sbml._config import ROOT_DIR

GOLDBETER = ROOT_DIR / "SbmlRefModels" / "src" / "reference" / "Goldbeter1991" / "Goldbeter1991.xml"


def _run(monkeypatch, *argv):
    """Invoke the CLI in-process with the given arguments."""
    monkeypatch.setattr("sys.argv", ["chaste-sbml", *argv])
    cli.main()


def test_version_exits_zero(monkeypatch):
    """--version prints and exits successfully."""
    monkeypatch.setattr("sys.argv", ["chaste-sbml", "--version"])
    with pytest.raises(SystemExit) as exc:
        cli.main()
    assert exc.value.code == 0


@pytest.mark.parametrize("model_type", ["generic", "srn", "cell-cycle"])
def test_generate_each_model_type(monkeypatch, tmp_path, model_type):
    """generate handles every --model-type and emits the model plus placeholder test."""
    _run(monkeypatch, "generate", str(GOLDBETER), "--model-type", model_type, "--output-dir", str(tmp_path))

    assert (tmp_path / "Goldbeter1991SbmlOdeSystem.hpp").is_file()
    assert (tmp_path / "TestGoldbeter1991Sbml.hpp").is_file()


def test_generate_no_tests(monkeypatch, tmp_path):
    """generate --no-tests writes the model but no placeholder test."""
    _run(monkeypatch, "generate", str(GOLDBETER), "--output-dir", str(tmp_path), "--no-tests")

    assert (tmp_path / "Goldbeter1991SbmlOdeSystem.hpp").is_file()
    assert not (tmp_path / "TestGoldbeter1991Sbml.hpp").exists()


def test_generate_test_output_dir(monkeypatch, tmp_path):
    """generate --test-output-dir routes the placeholder test to its own directory."""
    src_dir = tmp_path / "src"
    test_dir = tmp_path / "test"
    src_dir.mkdir()
    test_dir.mkdir()

    _run(
        monkeypatch,
        "generate",
        str(GOLDBETER),
        "--output-dir",
        str(src_dir),
        "--test-output-dir",
        str(test_dir),
    )

    assert (src_dir / "Goldbeter1991SbmlOdeSystem.hpp").is_file()
    assert (test_dir / "TestGoldbeter1991Sbml.hpp").is_file()


def test_copy_base_classes(monkeypatch, tmp_path):
    """copy-base-classes copies the C++ base classes into --output-dir."""
    _run(monkeypatch, "copy-base-classes", "--output-dir", str(tmp_path))

    assert (tmp_path / "AbstractSbmlOdeSystem.hpp").is_file()
