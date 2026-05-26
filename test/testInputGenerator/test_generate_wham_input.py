#!/usr/bin/env python3
import subprocess
import sys
from pathlib import Path


def main() -> int:
    repo = Path(__file__).resolve().parents[2]
    spec = repo / "test" / "testInputGenerator" / "uwham_1d.json"
    expected = repo / "test" / "testInputGenerator" / "uwham_1d.dat"
    generator = repo / "tools" / "generate_wham_input.py"

    result = subprocess.run(
        [sys.executable, str(generator), str(spec)],
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    if result.returncode != 0:
        sys.stderr.write(result.stderr)
        return result.returncode

    expected_text = expected.read_text()
    if result.stdout != expected_text:
        sys.stderr.write("Generated WHAM input did not match expected output.\n")
        sys.stderr.write("--- expected ---\n")
        sys.stderr.write(expected_text)
        sys.stderr.write("\n--- actual ---\n")
        sys.stderr.write(result.stdout)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
