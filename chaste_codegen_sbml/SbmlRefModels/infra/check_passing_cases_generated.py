#!/usr/bin/env python3
"""Ensure every case marked 'pass' in the status CSV generated a runnable test.

A case can be marked 'pass' in the status CSV yet not be generatable by the
codegen (e.g. its selected SBML version uses an unsupported feature such as a
fast reaction). Such a case silently never runs in CI, hiding the
misclassification. This script cross-checks the 'pass' rows in a case range
against the generated test files and fails if any are missing.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

PASS_STATUS = "pass"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="check_passing_cases_generated.py",
        description="Fail if any 'pass'-marked case in the range has no generated test.",
    )
    parser.add_argument("--status-csv", required=True, help="Path to sbml_test_suite_status.csv")
    parser.add_argument("--gen-test-dir", required=True, help="Generated semantic tests directory")
    parser.add_argument("--first-case", required=True, help="First case number (inclusive)")
    parser.add_argument("--last-case", required=True, help="Last case number (inclusive)")
    return parser.parse_args()


def passing_cases(status_csv: Path, first_case: int, last_case: int) -> list[int]:
    cases: list[int] = []
    with status_csv.open("r", newline="") as handle:
        reader = csv.reader(handle)
        next(reader)  # Skip header
        for row in reader:
            if not row:
                continue
            status_val = row[1].strip() if len(row) > 1 else ""
            case_val = row[0].strip()
            if status_val != PASS_STATUS or not case_val.isdigit():
                continue
            case_id = int(case_val)
            if first_case <= case_id <= last_case:
                cases.append(case_id)
    return cases


def has_generated_test(gen_test_dir: Path, case_id: int) -> bool:
    case_dir = gen_test_dir / f"{case_id:05d}"
    if not case_dir.is_dir():
        return False
    return any(p.name.startswith("Test") and p.suffix == ".hpp" for p in case_dir.iterdir())


def main() -> int:
    args = parse_args()

    status_csv = Path(args.status_csv)
    if not status_csv.is_file():
        print(f"Error: status file not found at '{status_csv}'", file=sys.stderr)
        return 1

    gen_test_dir = Path(args.gen_test_dir)
    if not gen_test_dir.is_dir():
        print(f"Error: generated test directory not found at '{gen_test_dir}'", file=sys.stderr)
        return 1

    if not re.fullmatch(r"\d+", args.first_case) or not re.fullmatch(r"\d+", args.last_case):
        print("Error: --first-case and --last-case must be integers.", file=sys.stderr)
        return 1

    first_case = int(args.first_case)
    last_case = int(args.last_case)

    expected = passing_cases(status_csv, first_case, last_case)
    missing = [c for c in expected if not has_generated_test(gen_test_dir, c)]

    if missing:
        print(
            f"Error: {len(missing)} case(s) marked 'pass' in range "
            f"{first_case:05d}-{last_case:05d} generated no test:",
            file=sys.stderr,
        )
        for case_id in missing:
            print(f"  {case_id:05d}", file=sys.stderr)
        print(
            "Each must either generate a runnable test or be reclassified in the status CSV.",
            file=sys.stderr,
        )
        return 1

    print(f"All {len(expected)} passing cases in range {first_case:05d}-{last_case:05d} generated a test.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
