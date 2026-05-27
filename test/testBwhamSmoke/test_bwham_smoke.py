#!/usr/bin/env python3

import math
import subprocess
import sys
import tempfile
from pathlib import Path


def fail(message):
    print(message, file=sys.stderr)
    return 1


def main():
    if len(sys.argv) != 4:
        return fail("usage: test_bwham_smoke.py WHAM_EXE INPUT_FILE DATA_ROOT")

    wham_exe = Path(sys.argv[1])
    input_file = Path(sys.argv[2])
    data_root = Path(sys.argv[3])

    with tempfile.TemporaryDirectory() as tmpdir:
        result = subprocess.run(
            [str(wham_exe), str(input_file), "-abspath", str(data_root)],
            cwd=tmpdir,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        if result.returncode != 0:
            print(result.stdout, end="")
            print(result.stderr, end="", file=sys.stderr)
            return fail(f"BWHAM smoke run failed with exit code {result.returncode}")

        output = Path(tmpdir) / "p.out"
        if not output.exists():
            return fail("BWHAM smoke run did not write p.out")

        lines = output.read_text().splitlines()
        if len(lines) != 31:
            return fail(f"BWHAM smoke output expected 31 lines, got {len(lines)}")

        rows = []
        for line in lines[1:]:
            fields = line.split()
            if len(fields) != 3:
                return fail(f"BWHAM smoke output row has {len(fields)} fields: {line}")
            rows.append([float(field) for field in fields])

        if not math.isclose(rows[0][0], 0.5, abs_tol=1e-12):
            return fail(f"BWHAM smoke first bin center expected 0.5, got {rows[0][0]}")
        if not math.isclose(rows[-1][0], 29.5, abs_tol=1e-12):
            return fail(f"BWHAM smoke last bin center expected 29.5, got {rows[-1][0]}")

        for row in rows:
            if not all(math.isfinite(value) for value in row):
                return fail(f"BWHAM smoke output contains non-finite value: {row}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
