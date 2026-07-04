"""Tests for copying the C++ base classes."""

from chaste_sbml import copy_base_classes
from chaste_sbml._base_classes import BASE_CLASS_DIR, base_class_files


def test_base_class_files_found():
    """The package ships a non-empty set of base-class files, including fortests helpers."""
    names = {p.name for p in base_class_files()}
    assert "AbstractSbmlOdeSystem.hpp" in names
    assert "SbmlMath.cpp" in names
    assert "SbmlTestHelpers.hpp" in names


def test_copy_base_classes(tmp_path):
    """Copying writes every base-class file, preserving the fortests/ layout and contents."""
    written = copy_base_classes(tmp_path)

    assert written, "expected at least one file to be written"

    for source in base_class_files():
        destination = tmp_path / source.relative_to(BASE_CLASS_DIR)
        assert destination.is_file()
        assert destination.read_bytes() == source.read_bytes()

    # The fortests/ layout is preserved.
    assert (tmp_path / "fortests" / "SbmlTestHelpers.hpp").is_file()
