#!/usr/bin/env python3
"""Filter semantic SBML test cases based on status CSV and case range."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

VALID_STATUSES = {"pass", "fail_known", "fail_unknown"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="filter_test_suite_cases.py",
        description=(
            "Reads the status CSV and prints selected_count and test_regex as key=value lines. "
            "Rows are selected when status is pass, fail_known, or fail_unknown, and case is in range."
        ),
    )
    parser.add_argument("--status-csv", required=True, help="Path to sbml_test_suite_status.csv")
    parser.add_argument("--first-case", required=True, help="First case number (inclusive)")
    parser.add_argument("--last-case", required=True, help="Last case number (inclusive)")
    parser.add_argument(
        "--output-file",
        default="",
        help="Optional file path to append selected_count/test_regex key=value lines",
    )
    return parser.parse_args()


def select_cases(status_csv: Path, first_case: int, last_case: int) -> list[str]:
    selected_case_ids: list[str] = []
    with status_csv.open("r", newline="") as handle:
        reader = csv.reader(handle)
        next(reader)  # Skip header
        for row in reader:
            if not row:
                continue

            status_val = row[1].strip() if len(row) > 1 else ""
            if status_val not in VALID_STATUSES:
                continue

            case_val = row[0].strip()
            if not case_val.isdigit():
                continue

            case_id = int(case_val)
            if first_case <= case_id <= last_case:
                selected_case_ids.append(f"{case_id:05d}")
    return selected_case_ids


def create_test_regex(selected_case_ids: list[str]) -> str:
    if not selected_case_ids:
        return "^$"
    selected_cases = "|".join(selected_case_ids)
    return f"^TestSemantic({selected_cases})L([0-9]|[0-9][0-9])V([0-9]|[0-9][0-9])Sbml$"


def write_output(output_file: Path, key_value_lines: list[str]) -> None:
    with output_file.open("a") as handle:
        for line in key_value_lines:
            handle.write(f"{line}\n")


def main() -> int:
    args = parse_args()

    status_csv = Path(args.status_csv)
    if not status_csv.is_file():
        print(f"Error: status file not found at '{status_csv}'", file=sys.stderr)
        return 1

    if not re.fullmatch(r"\d+", args.first_case) or not re.fullmatch(r"\d+", args.last_case):
        print("Error: --first-case and --last-case must be integers.", file=sys.stderr)
        return 1

    first_case = int(args.first_case)
    last_case = int(args.last_case)

    selected_case_ids = select_cases(status_csv, first_case, last_case)
    selected_count = len(selected_case_ids)

    if selected_count == 0:
        print(f"No relevant semantic tests found in range {first_case:05d}-{last_case:05d}")
        return 1

    test_regex = create_test_regex(selected_case_ids)

    output_lines = [f"selected_count={selected_count}", f"test_regex={test_regex}"]
    for line in output_lines:
        print(line)

    if args.output_file:
        output_path = Path(args.output_file)
        write_output(output_path, output_lines)

    print(f"Selected {selected_count} semantic tests for range {first_case:05d}-{last_case:05d}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
