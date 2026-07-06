#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd -- "${script_dir}/.." && pwd)"
output_dir="${script_dir}/output"
build_dir="${script_dir}/build"
chaste_dir="${CHASTE_SOURCE_DIR:-${repo_root}/../Chaste}"

if [[ $# -gt 0 ]]; then
	echo "Usage: $(basename "$0")" >&2
	exit 1
fi

if ! command -v cmake >/dev/null 2>&1; then
	echo "Error: cmake is not available on PATH." >&2
	exit 1
fi

if [[ ! -d "${chaste_dir}" ]]; then
	echo "Error: Chaste source directory not found at '${chaste_dir}'." >&2
	echo "Set CHASTE_SOURCE_DIR to override the default checkout path." >&2
	exit 1
fi

mkdir -p "${output_dir}" "${build_dir}"
export CHASTE_TEST_OUTPUT="${output_dir}"

cd "${build_dir}"
nice -n 9 cmake -DChaste_UPDATE_PROVENANCE=OFF "${chaste_dir}"

