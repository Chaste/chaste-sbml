#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<EOF
Usage: $(basename "$0") FIRST_CASE LAST_CASE

Generate semantic test cases from a sibling SBML Test Suite checkout.

Environment:
  VENV_DIR             Python virtualenv directory (default: <repo>/venv)
  SBML_TEST_SUITE_DIR  Upstream sbml-test-suite checkout (default: <repo>/../sbml-test-suite)
EOF
}

if [[ $# -ne 2 ]]; then
  usage >&2
  exit 1
fi

first_case="${1}"
last_case="${2}"

if ! [[ "${first_case}" =~ ^[0-9]+$ && "${last_case}" =~ ^[0-9]+$ ]]; then
  echo "Error: FIRST_CASE and LAST_CASE must be integers." >&2
  exit 1
fi

if [[ "${first_case}" -gt "${last_case}" ]]; then
  echo "Error: FIRST_CASE must be less than or equal to LAST_CASE." >&2
  exit 1
fi

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd -- "${script_dir}/../../.." && pwd)"
venv_dir="${VENV_DIR:-${repo_root}/venv}"
python_bin="${venv_dir}/bin/python"
generator="${repo_root}/chaste_sbml/SbmlRefModels/generate_cases.py"
sbml_test_suite_dir="${SBML_TEST_SUITE_DIR:-${repo_root}/../sbml-test-suite}"
gen_src_dir="${repo_root}/chaste_sbml/SbmlRefModels/src/cases/semantic"
gen_test_dir="${repo_root}/chaste_sbml/SbmlRefModels/test/cases/semantic"
test_pack_file="${repo_root}/chaste_sbml/SbmlRefModels/test/WeeklyTestPack.txt"

if [[ ! -x "${python_bin}" ]]; then
  echo "Error: Python interpreter not found at '${python_bin}'." >&2
  echo "Set VENV_DIR to override the default virtualenv location." >&2
  exit 1
fi

if [[ ! -f "${generator}" ]]; then
  echo "Error: generator script not found at '${generator}'." >&2
  exit 1
fi

if [[ ! -d "${sbml_test_suite_dir}" ]]; then
  echo "Error: SBML Test Suite checkout not found at '${sbml_test_suite_dir}'." >&2
  echo "Set SBML_TEST_SUITE_DIR to override the default sibling checkout path." >&2
  exit 1
fi

mkdir -p "${gen_src_dir}" "${gen_test_dir}"

cd "${repo_root}"
"${python_bin}" "${generator}" \
  --sbml-test-suite-dir "${sbml_test_suite_dir}" \
  --gen-src-dir "${gen_src_dir}" \
  --gen-test-dir "${gen_test_dir}" \
  --test-pack-file "${test_pack_file}" \
  --first-case "${first_case}" \
  --last-case "${last_case}"
