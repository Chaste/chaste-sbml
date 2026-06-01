#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: filter_test_suite_cases.sh --status-csv PATH --first-case N --last-case N [--output-file PATH]

Reads the status CSV and prints selected_count and test_regex as key=value lines.
Rows are selected when status is pass, fail_known, or fail_unknown, and case is in [first-case, last-case].

When --output-file is provided, the same key=value lines are also appended to that file.
EOF
}

status_csv=""
first_case=""
last_case=""
output_file=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --status-csv)
      status_csv="${2:-}"
      shift 2
      ;;
    --first-case)
      first_case="${2:-}"
      shift 2
      ;;
    --last-case)
      last_case="${2:-}"
      shift 2
      ;;
    --output-file)
      output_file="${2:-}"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Error: unknown argument '$1'" >&2
      usage >&2
      exit 1
      ;;
  esac
done

if [[ -z "${status_csv}" || -z "${first_case}" || -z "${last_case}" ]]; then
  echo "Error: --status-csv, --first-case, and --last-case are required." >&2
  usage >&2
  exit 1
fi

if [[ ! -f "${status_csv}" ]]; then
  echo "Error: status file not found at '${status_csv}'" >&2
  exit 1
fi

if ! [[ "${first_case}" =~ ^[0-9]+$ && "${last_case}" =~ ^[0-9]+$ ]]; then
  echo "Error: --first-case and --last-case must be integers." >&2
  exit 1
fi

selected_cases="$({
  awk -F',' -v first="${first_case}" -v last="${last_case}" '
    {
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", $1)
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", $2)
    }

    /^[[:space:]]*$/ { next }
    tolower($1) == "test" { next }

    $1 ~ /^[0-9]+$/ && ($2 == "pass" || $2 == "fail_known" || $2 == "fail_unknown") {
      case_id = $1 + 0
      if (case_id >= first && case_id <= last) {
        printf "%05d\n", case_id
      }
    }
  ' "${status_csv}" | paste -sd'|' -
} || true)"

selected_count="0"
test_regex="^$"

if [[ -n "${selected_cases}" ]]; then
  selected_count="$(awk -F'|' '{print NF}' <<< "${selected_cases}")"
  test_regex="^TestSemantic(${selected_cases})L([0-9]|[0-9][0-9])V([0-9]|[0-9][0-9])Sbml$"
fi

echo "selected_count=${selected_count}"
echo "test_regex=${test_regex}"

if [[ -n "${output_file}" ]]; then
  echo "selected_count=${selected_count}" >> "${output_file}"
  echo "test_regex=${test_regex}" >> "${output_file}"
fi

if [[ "${selected_count}" == "0" ]]; then
  echo "No relevant semantic tests found in matrix range ${first_case}-${last_case}"
else
  echo "Selected ${selected_count} semantic tests for range ${first_case}-${last_case}"
fi
