import logging
import re

import pytest

import chaste_codegen_sbml as cg
from chaste_codegen_sbml._config import ROOT_DIR

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
        # Simplify whitespace and exclude blank lines
        line = re.sub(r"\s+", " ", line)
        line = line.strip()
        if line:
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

    for line_a, line_b in zip(lines_a, lines_b):
        if line_a != line_b:
            return f'"{line_a}" != "{line_b}"'

    return ""


@pytest.mark.parametrize(
    (
        "filename",
        "model_name",
    ),
    [("Goldbeter1991SrnModel", "Goldbeter1991")],
)
def test_generation(tmp_path, filename, model_name):
    """
    Check generated model against reference.
    """
    ref_dir = ROOT_DIR / "data" / "reference_models" / filename
    ref_sbml = ref_dir / f"{filename}.xml"
    ref_cpp = ref_dir / f"{filename}.cpp"
    ref_hpp = ref_dir / f"{filename}.hpp"

    logger.info(f"Converting: {ref_sbml}")
    chaste_model = cg.ChasteSRNModel(ref_sbml, model_name)
    chaste_model.write_chaste_code(output_directory=tmp_path)

    gen_hpp = tmp_path / chaste_model.hpp_filename
    gen_cpp = tmp_path / chaste_model.cpp_filename

    hpp_diff = code_diff(ref_hpp, gen_hpp)
    assert hpp_diff == "", hpp_diff

    cpp_diff = code_diff(ref_cpp, gen_cpp)
    assert cpp_diff == "", cpp_diff
