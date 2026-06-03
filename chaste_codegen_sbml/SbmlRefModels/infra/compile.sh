#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
build_dir="${script_dir}/build"
target="${CHASTE_BUILD_TARGET:-project_SbmlRefModels}"
jobs="${JOBS:-16}"

if [[ $# -gt 0 ]]; then
	echo "Usage: $(basename "$0")" >&2
	exit 1
fi

if ! [[ "${jobs}" =~ ^[0-9]+$ ]] || [[ "${jobs}" -lt 1 ]]; then
	echo "Error: JOBS must be a positive integer (got '${jobs}')." >&2
	exit 1
fi

if ! command -v cmake >/dev/null 2>&1; then
	echo "Error: cmake is not available on PATH." >&2
	exit 1
fi

if [[ ! -f "${build_dir}/CMakeCache.txt" ]]; then
	echo "Error: build directory is not configured at '${build_dir}'. Run infra/configure.sh first." >&2
	exit 1
fi

cd "${build_dir}"
nice -n 9 cmake --build . --target "${target}" --parallel "${jobs}"
