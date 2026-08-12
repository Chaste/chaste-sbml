"""Tests for console script."""

import logging
import subprocess

from chaste_sbml import __version__
from chaste_sbml._config import ROOT_DIR

logger = logging.getLogger(__name__)

GOLDBETER = ROOT_DIR / "SbmlRefModels" / "src" / "reference" / "Goldbeter1991" / "Goldbeter1991.xml"


def test_help():
    """The help documents the default SBML file argument and the copy-base-classes switch."""
    output = subprocess.check_output(["chaste-sbml", "-h"]).decode("ascii")

    assert "sbml_file" in output
    assert "--copy-base-classes" in output


def test_version():
    """Test version message."""
    expected = f"chaste-sbml {__version__}"

    output = subprocess.check_output(["chaste-sbml", "--version"]).decode("ascii").strip()

    assert output == expected


def test_requires_sbml_file():
    """Running with no SBML file and no --copy-base-classes is a usage error (exit code 2)."""
    result = subprocess.run(["chaste-sbml"], capture_output=True, text=True)

    assert result.returncode == 2


def test_copy_base_classes_rejects_generation_options():
    """--copy-base-classes takes only --output-dir; an SBML file or generation option is a usage error."""
    for extra in (
        ["model.xml"],
        ["--model-type", "srn"],
        ["--no-tests"],
        ["--timescale", "s"],
        ["--test-output-dir", "test/"],
    ):
        result = subprocess.run(
            ["chaste-sbml", "--copy-base-classes", *extra],
            capture_output=True,
            text=True,
        )

        assert result.returncode == 2, extra


def test_copy_base_classes_accepts_output_dir(tmp_path):
    """--copy-base-classes with --output-dir copies the base classes."""
    result = subprocess.run(
        ["chaste-sbml", "--copy-base-classes", "--output-dir", str(tmp_path)],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr
    assert (tmp_path / "AbstractSbmlOdeSystem.hpp").is_file()


def test_generate_emits_placeholder_test(tmp_path):
    """Generation produces a placeholder test file alongside the model by default."""
    result = subprocess.run(
        ["chaste-sbml", str(GOLDBETER), "--output-dir", str(tmp_path)],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr
    assert (tmp_path / "Goldbeter1991SbmlOdeSystem.hpp").is_file()
    assert (tmp_path / "TestGoldbeter1991Sbml.hpp").is_file()


def test_generate_no_tests_skips_placeholder_test(tmp_path):
    """--no-tests writes the model but no placeholder test."""
    result = subprocess.run(
        ["chaste-sbml", str(GOLDBETER), "--output-dir", str(tmp_path), "--no-tests"],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr
    assert (tmp_path / "Goldbeter1991SbmlOdeSystem.hpp").is_file()
    assert not (tmp_path / "TestGoldbeter1991Sbml.hpp").exists()


def test_generate_test_output_dir_routes_test(tmp_path):
    """--test-output-dir places the placeholder test in its own directory."""
    src_dir = tmp_path / "src"
    test_dir = tmp_path / "test"
    src_dir.mkdir()
    test_dir.mkdir()

    result = subprocess.run(
        [
            "chaste-sbml",
            str(GOLDBETER),
            "--output-dir",
            str(src_dir),
            "--test-output-dir",
            str(test_dir),
        ],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr
    assert (src_dir / "Goldbeter1991SbmlOdeSystem.hpp").is_file()
    assert (test_dir / "TestGoldbeter1991Sbml.hpp").is_file()
    assert not (src_dir / "TestGoldbeter1991Sbml.hpp").exists()
