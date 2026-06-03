#!/usr/bin/env bash
set -euo pipefail

usage() {
	echo "Usage: $(basename "$0")" >&2
}

if [[ $# -ne 0 ]]; then
	usage
	exit 1
fi

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd -- "${script_dir}/../../.." && pwd)"
src_cases_dir="${repo_root}/chaste_codegen_sbml/SbmlRefModels/src/cases/semantic"
test_cases_dir="${repo_root}/chaste_codegen_sbml/SbmlRefModels/test/cases/semantic"
weekly_test_pack="${repo_root}/chaste_codegen_sbml/SbmlRefModels/test/WeeklyTestPack.txt"

for dir in "${src_cases_dir}" "${test_cases_dir}"; do
	if [[ ! -d "${dir}" ]]; then
		echo "Error: expected directory '${dir}' was not found." >&2
		exit 1
	fi

	find "${dir}" -mindepth 1 -maxdepth 1 -type d -name '[0-9][0-9][0-9][0-9][0-9]' -exec rm -rf {} +
done

: > "${weekly_test_pack}"
echo "Reset generated semantic cases and cleared ${weekly_test_pack}"

