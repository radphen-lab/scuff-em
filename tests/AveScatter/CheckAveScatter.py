#!/usr/bin/env python3
"""Compare AVSCAT tables using the repository's KEY/DATA checklist format."""

import argparse
import math
from pathlib import Path
import sys


def checklist(path):
    keys, data = [], []
    absolute, relative = 0.0, 0.01
    for line in path.read_text().splitlines():
        tokens = line.split()
        if not tokens or tokens[0].startswith("#"):
            continue
        if tokens[0] == "KEY" and len(tokens) == 3:
            keys.append(int(tokens[2]) - 1)
        elif tokens[0] == "DATA" and len(tokens) == 3:
            data.append((tokens[1], int(tokens[2]) - 1))
        elif tokens[0] == "ABSTOL" and len(tokens) == 2:
            absolute = float(tokens[1])
        elif tokens[0] == "RELTOL" and len(tokens) == 2:
            relative = float(tokens[1])
        else:
            raise ValueError(f"{path}: invalid checklist line {line!r}")
    if not keys or not data or min(keys + [column for _, column in data]) < 0:
        raise ValueError(f"{path}: expected positive KEY and DATA columns")
    if not math.isfinite(absolute) or not math.isfinite(relative) or min(absolute, relative) < 0:
        raise ValueError(f"{path}: invalid comparison tolerances")
    return keys, data, absolute, relative


def records(path, keys, data):
    result = {}
    columns = max(keys + [column for _, column in data]) + 1
    for number, line in enumerate(path.read_text().splitlines(), 1):
        tokens = line.split()
        if not tokens or tokens[0].startswith("#"):
            continue
        if len(tokens) != columns:
            raise ValueError(f"{path}:{number}: expected {columns} columns")
        key = tuple(tokens[column] for column in keys)
        values = [float(tokens[column]) for _, column in data]
        if not all(math.isfinite(value) for value in values):
            raise ValueError(f"{path}:{number}: non-finite data")
        if key in result:
            raise ValueError(f"{path}:{number}: duplicate key {key}")
        result[key] = values
    if not result:
        raise ValueError(f"{path}: empty data file")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--checklist", required=True, type=Path)
    args = parser.parse_args()
    try:
        keys, data, absolute, relative = checklist(args.checklist)
        observed = records(args.data, keys, data)
        expected = records(args.reference, keys, data)
        if observed.keys() != expected.keys():
            raise ValueError(f"{args.data}: missing or extra transform/frequency/surface keys")
        passed = True
        for key, values in expected.items():
            for (name, _), actual, reference in zip(data, observed[key], values):
                if not math.isclose(actual, reference, rel_tol=relative, abs_tol=absolute):
                    print(f"FAIL {key} {name}: computed={actual:.12g} reference={reference:.12g}")
                    passed = False
        print(f"{'PASS' if passed else 'FAIL'} {args.data.name}, {len(expected)} rows")
        return 0 if passed else 1
    except (OSError, ValueError) as error:
        print(f"FAIL {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
