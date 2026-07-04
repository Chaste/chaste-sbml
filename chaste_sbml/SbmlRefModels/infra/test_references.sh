#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<EOF
Usage: $(basename "$0") [PARTIAL_NAME]

Run tests on reference models.

Default behavior:
  With no argument, runs all tests labeled Continuous_project_SbmlRefModels.

Arguments:
  PARTIAL_NAME  Optional string to match against test names in addition to label filtering.

Environment:
  JOBS          Number of parallel ctest jobs (default: 16)
  NICENESS      niceness value passed to nice (default: 9)
EOF
}

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
build_dir="${script_dir}/build"
output_dir="${script_dir}/output"
jobs="${JOBS:-16}"
niceness="${NICENESS:-9}"
partial_name="${1:-}"
test_label="Continuous_project_SbmlRefModels"

if [[ "${partial_name}" == "-h" || "${partial_name}" == "--help" ]]; then
  usage
  exit 0
fi

if [[ $# -gt 1 ]]; then
  usage >&2
  exit 1
fi

if ! [[ "${jobs}" =~ ^[0-9]+$ ]] || [[ "${jobs}" -lt 1 ]]; then
  echo "Error: JOBS must be a positive integer (got '${jobs}')." >&2
  exit 1
fi

if ! [[ "${niceness}" =~ ^[0-9]+$ ]]; then
  echo "Error: NICENESS must be a non-negative integer (got '${niceness}')." >&2
  exit 1
fi

if ! command -v ctest >/dev/null 2>&1; then
  echo "Error: ctest is not available on PATH." >&2
  exit 1
fi

if [[ ! -f "${build_dir}/CMakeCache.txt" ]]; then
  echo "Error: build directory is not configured at '${build_dir}'. Run infra/configure.sh first." >&2
  exit 1
fi

mkdir -p "${output_dir}"
export CHASTE_TEST_OUTPUT="${output_dir}"
echo "CHASTE_TEST_OUTPUT=${CHASTE_TEST_OUTPUT}"

test_regex=""
if [[ -n "${partial_name}" ]]; then
  test_regex=".*${partial_name}.*"
  echo "Running tests with label '${test_label}' and name filter '${partial_name}'"
else
  echo "Running all tests with label '${test_label}'"
fi

cd "${build_dir}"
if [[ -n "${test_regex}" ]]; then
  nice -n "${niceness}" ctest -j"${jobs}" -V --output-on-failure -L "${test_label}" -R "${test_regex}"
else
  nice -n "${niceness}" ctest -j"${jobs}" -V --output-on-failure -L "${test_label}"
fi
