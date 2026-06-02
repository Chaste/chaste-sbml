#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<EOF
Usage: $(basename "$0") [FIRST_CASE LAST_CASE]

Run generated semantic tests from scripts/build.

Default behavior:
  With no argument, reads sbml_test_suite_status.csv and runs rows with status
  'pass', 'fail_known', or 'fail_unknown'.

Arguments:
  FIRST_CASE    Optional integer first semantic case number in a range.
  LAST_CASE     Optional integer last semantic case number in a range.

Environment:
  JOBS          Number of parallel ctest jobs (default: 16)
  NICENESS      niceness value passed to nice (default: 9)
EOF
}

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
build_dir="${script_dir}/build"
output_dir="${script_dir}/output"
status_csv="${script_dir}/../chaste_codegen_sbml/SbmlRefModels/test/data/sbml_test_suite_status.csv"
filter_script="${script_dir}/../chaste_codegen_sbml/SbmlRefModels/infra/filter_test_suite_cases.py"
jobs="${JOBS:-16}"
niceness="${NICENESS:-9}"
case_filter_start="${1:-}"
case_filter_end="${2:-}"

set_pass_test_regex_from_status_csv() {
  local filter_output
  local selected_count

  if [[ ! -x "${filter_script}" && ! -f "${filter_script}" ]]; then
    echo "Error: filter script not found at '${filter_script}'." >&2
    exit 1
  fi

  if ! filter_output="$(python3 "${filter_script}" \
    --status-csv "${status_csv}" \
    --first-case 0 \
    --last-case 99999)"; then
    exit 1
  fi

  selected_count="$(printf "%s\n" "${filter_output}" | grep '^selected_count=' | tail -n1 | cut -d'=' -f2-)"
  test_regex="$(printf "%s\n" "${filter_output}" | grep '^test_regex=' | tail -n1 | cut -d'=' -f2-)"

  if [[ -z "${selected_count}" || "${selected_count}" == "0" || -z "${test_regex}" || "${test_regex}" == "^$" ]]; then
    echo "Error: no passing tests found in '${status_csv}'." >&2
    exit 1
  fi

  echo "Running ${selected_count} selected tests (pass/fail_known/fail_unknown) from ${status_csv}"
}

set_pass_test_regex_from_status_csv_range() {
  local start_case="$1"
  local end_case="$2"
  local filter_output
  local selected_count

  if [[ ! -x "${filter_script}" && ! -f "${filter_script}" ]]; then
    echo "Error: filter script not found at '${filter_script}'." >&2
    exit 1
  fi

  if ! filter_output="$(python3 "${filter_script}" \
    --status-csv "${status_csv}" \
    --first-case "${start_case}" \
    --last-case "${end_case}")"; then
    exit 1
  fi

  selected_count="$(printf "%s\n" "${filter_output}" | grep '^selected_count=' | tail -n1 | cut -d'=' -f2-)"
  test_regex="$(printf "%s\n" "${filter_output}" | grep '^test_regex=' | tail -n1 | cut -d'=' -f2-)"

  if [[ -z "${selected_count}" || "${selected_count}" == "0" || -z "${test_regex}" || "${test_regex}" == "^$" ]]; then
    echo "Error: no passing tests found in '${status_csv}' for range ${start_case}-${end_case}." >&2
    exit 1
  fi

  echo "Running ${selected_count} selected tests (pass/fail_known/fail_unknown) from ${status_csv} for range ${start_case}-${end_case}"
}

if [[ "${case_filter_start}" == "-h" || "${case_filter_start}" == "--help" ]]; then
  usage
  exit 0
fi

if [[ $# -gt 2 ]]; then
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

test_regex="TestSemantic"
if [[ -n "${case_filter_start}" && -n "${case_filter_end}" ]]; then
  if [[ "${case_filter_start}" =~ ^[0-9]+$ && "${case_filter_end}" =~ ^[0-9]+$ ]]; then
    start_case="$(printf "%05d" "$((10#${case_filter_start}))")"
    end_case="$(printf "%05d" "$((10#${case_filter_end}))")"
    set_pass_test_regex_from_status_csv_range "${start_case}" "${end_case}"
  else
    echo "Error: when passing two arguments, both must be integer case numbers." >&2
    usage >&2
    exit 1
  fi
elif [[ -n "${case_filter_start}" || -n "${case_filter_end}" ]]; then
  echo "Error: pass either no arguments or a two-number semantic case range." >&2
  usage >&2
  exit 1
else
  set_pass_test_regex_from_status_csv
fi

cd "${build_dir}"
nice -n "${niceness}" ctest -j"${jobs}" -V -R "${test_regex}"
