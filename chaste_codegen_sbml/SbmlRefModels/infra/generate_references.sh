#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<EOF
Usage: $(basename "$0")

Regenerate the Chaste reference models from their SBML sources.

Environment:
  VENV_DIR  Python virtualenv directory (default: <repo>/venv)
EOF
}

if [[ $# -ne 0 ]]; then
  usage >&2
  exit 1
fi

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd -- "${script_dir}/../../.." && pwd)"
venv_dir="${VENV_DIR:-${repo_root}/venv}"
python_bin="${venv_dir}/bin/python"
generator="${repo_root}/chaste_codegen_sbml/SbmlRefModels/generate_references.py"
reference_dir="${repo_root}/chaste_codegen_sbml/SbmlRefModels/src/reference"

if [[ ! -x "${python_bin}" ]]; then
  echo "Error: Python interpreter not found at '${python_bin}'." >&2
  echo "Set VENV_DIR to override the default virtualenv location." >&2
  exit 1
fi

if [[ ! -f "${generator}" ]]; then
  echo "Error: generator script not found at '${generator}'." >&2
  exit 1
fi

if [[ ! -d "${reference_dir}" ]]; then
  echo "Error: reference directory not found at '${reference_dir}'." >&2
  exit 1
fi

cd "${repo_root}"
"${python_bin}" "${generator}" \
  --reference-dir "${reference_dir}"
