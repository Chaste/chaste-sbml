"""Tests for code generation."""

import logging
import re
from difflib import Differ

import pytest

from chaste_codegen_sbml import ChasteSbmlModel
from chaste_codegen_sbml._config import ROOT_DIR, ModelType

logger = logging.getLogger(__name__)


def load_source_lines(source_file: str) -> list[str]:
    """
    Load C++ source code from file and strip comments and whitespace.

    :param source_file: Path to source file.
    :return: Source code as string.
    """
    with open(source_file, "r") as sf:
        source = sf.read()

    # Strip comments
    source = re.sub(r"//.*", "", source)
    source = re.sub(r"/\*.*?\*/", " ", source, flags=re.DOTALL)

    source_lines = []
    for line in source.split("\n"):
        # Exclude blank lines
        line = line.strip()
        if line:
            # Merge whitespace sequences
            line = re.sub(r"\s+", " ", line)
            source_lines.append(line)

    return source_lines


def code_diff(file_a: str, file_b: str) -> str:
    """
    Compare two C++ source files and return the first line that differs.

    :param file_a: Path to first source file.
    :param file_b: Path to second source file.
    :return: First differing line.
    """
    lines_a = load_source_lines(file_a)
    lines_b = load_source_lines(file_b)

    if lines_a and not lines_b:
        return lines_a[0]

    if lines_b and not lines_a:
        return lines_b[0]

    differ = Differ()

    for line_a, line_b in zip(lines_a, lines_b):
        if line_a != line_b:
            diff = differ.compare([line_a], [line_b])
            return "\n".join(diff)

    return ""


@pytest.mark.parametrize(
    ("model_name", "model_type"),
    [
        ("Goldbeter1991", ModelType.SRN),
        ("Tan2014", ModelType.SRN),
        ("VanLeeuwen2007", ModelType.SRN),
        ("VanLeeuwen2007NonDim", ModelType.SRN),
        ("Chen2000", ModelType.CELL_CYCLE),
        ("Chen2004", ModelType.CELL_CYCLE),
        ("Gardner1998", ModelType.CELL_CYCLE),
        ("TysonNovak2001", ModelType.CELL_CYCLE),
    ],
)
def test_generation(tmp_path, model_name, model_type):
    """Check generated model against reference."""
    type_string = ""
    if model_type == ModelType.SRN:
        type_string = "Srn"
    elif model_type == ModelType.CELL_CYCLE:
        type_string = "CellCycle"

    ref_dir = ROOT_DIR / "SbmlRefModels" / "src" / "reference" / model_name
    gen_dir = tmp_path

    # import pathlib
    # gen_dir = pathlib.Path(".").absolute()

    gen_ode_hpp = gen_dir / f"{model_name}SbmlOdeSystem.hpp"
    gen_ode_cpp = gen_dir / f"{model_name}SbmlOdeSystem.cpp"

    gen_type_hpp = gen_dir / f"{model_name}Sbml{type_string}Model.hpp"
    gen_type_cpp = gen_dir / f"{model_name}Sbml{type_string}Model.cpp"

    ref_ode_hpp = ref_dir / f"{model_name}SbmlOdeSystem.hpp"
    ref_ode_cpp = ref_dir / f"{model_name}SbmlOdeSystem.cpp"

    ref_type_hpp = ref_dir / f"{model_name}Sbml{type_string}Model.hpp"
    ref_type_cpp = ref_dir / f"{model_name}Sbml{type_string}Model.cpp"

    sbml_file = ref_dir / f"{model_name}.xml"
    logger.info(f"Converting: {sbml_file}")
    chaste_model = ChasteSbmlModel(sbml_file, model_type=model_type)
    chaste_model.write(gen_dir)

    ode_hpp_diff = code_diff(ref_ode_hpp, gen_ode_hpp)
    assert ode_hpp_diff == "", ode_hpp_diff

    ode_cpp_diff = code_diff(ref_ode_cpp, gen_ode_cpp)
    assert ode_cpp_diff == "", ode_cpp_diff

    type_hpp_diff = code_diff(ref_type_hpp, gen_type_hpp)
    assert type_hpp_diff == "", type_hpp_diff

    type_cpp_diff = code_diff(ref_type_cpp, gen_type_cpp)
    assert type_cpp_diff == "", type_cpp_diff
