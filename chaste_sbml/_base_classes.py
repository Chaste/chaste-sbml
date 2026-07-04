"""Copy the C++ base classes that generated Chaste code depends on.

Generated model code ``#include``\\s and subclasses these files (e.g. ``AbstractSbmlOdeSystem``),
so a user project needs its own copy of them. They are shipped as package data, so a consumer
can obtain a matching set for whichever version of this package is installed, rather than
vendoring their own copy and letting it drift.
"""

import pathlib
import shutil

from ._config import ROOT_DIR

# The base classes live alongside the reference models and test cases under SbmlRefModels/src.
# Only the top-level files here and those in fortests/ are base classes; the model subdirectories
# (cases/, reference/) are deliberately excluded.
BASE_CLASS_DIR = ROOT_DIR / "SbmlRefModels" / "src"


def base_class_files() -> list[pathlib.Path]:
    """Return the C++ base-class files shipped with this package.

    :return: Absolute paths to the base-class headers and sources.
    """
    files = sorted(BASE_CLASS_DIR.glob("*.[hc]pp"))
    files += sorted((BASE_CLASS_DIR / "fortests").glob("*.[hc]pp"))
    return files


def copy_base_classes(output_directory=None) -> list[pathlib.Path]:
    """Copy the C++ base classes into ``output_directory`` (default: current directory).

    The fortests/ layout is preserved. Existing files are overwritten.

    :param output_directory: The directory to copy the base classes into.
    :return: The paths of the files written.
    """
    root_dir = pathlib.Path(output_directory) if output_directory else pathlib.Path.cwd()

    written = []
    for source in base_class_files():
        destination = root_dir / source.relative_to(BASE_CLASS_DIR)
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, destination)
        written.append(destination)

    return written
