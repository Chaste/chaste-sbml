import logging
import re

import pytest

import chaste_codegen_sbml as cg


logger = logging.getLogger(__name__)


def load_source(source_file):
    """
    Load C++ source code from file and strip comments and whitespace.

    :param source_file: Path to source file.
    :return: Source code as string.
    """
    with open(source_file, "r") as sf:
        source = "\n".join(line.strip() for line in sf)

    # Strip comments
    source = re.sub(r"//.*", "", source)
    source = re.sub(r"/\*.*?\*/", " ", source, flags=re.DOTALL)

    # Simplify whitespace
    source = re.sub(r"\b\s+|\s+\b", " ", source)
    source = re.sub(r"\B\s+|\s+\B", " ", source)
    print(source)
    return source


def compare_code(file_a, file_b):
    """
    Compare two C++ source files.

    :param file_a: Path to first source file.
    :param file_b: Path to second source file.
    :return: True if source code is identical, False otherwise
    """
    code_a = load_source(file_a)
    code_b = load_source(file_b)
    return code_a == code_b


@pytest.mark.parametrize("model", ["goldbeter_1991"])
def test_generation(tmp_path, model):
    """
    Check generated model against reference.
    """
    sbml = f"{model}.xml"
    cpp = f"{model}.cpp"
    hpp = f"{model}.hpp"

    ref_dir = cg.ROOT_DIR / "data" / "reference_models" / model
    ref_sbml = ref_dir / sbml
    ref_cpp = ref_dir / cpp
    ref_hpp = ref_dir / hpp

    gen_cpp = tmp_path / cpp
    gen_hpp = tmp_path / hpp

    logger.info(f"Converting: {ref_sbml}")
    cg.Generate(str(ref_sbml), output_directory=str(tmp_path))

    assert compare_code(ref_cpp, gen_cpp)
    assert compare_code(ref_hpp, gen_hpp)
