# AGENTS.md

## Purpose

This repository contains the C++ WHAM implementation used to build the `Wham` executable and supporting libraries. This file is an operational guide for coding agents; use `README.md` for scientific background and validation context.

## Repository Layout

- `src/`: WHAM core implementation, bias models, calculation strategies, reweighting, and time-series handling.
- `tools/`: input parsing, command-line handling, filesystem helpers, and shared utilities.
- `parallel/`: OpenMP and MPI-related helpers.
- `test/`: CTest wiring, shell test runner, fixture inputs, and golden reference outputs.
- `scripts/`: example WHAM input files.
- `cmake/`: CMake helper modules.
- `Eigen/` and `LBFGS/`: vendored dependencies. Avoid editing these unless the task explicitly targets vendored code.

## Build

Use an out-of-source build:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

The main executable is built as `build/bin/Wham`.

## Official Tests

Run the registered CTest suite with:

```bash
ctest --test-dir build --output-on-failure
```

The current official CTest suite registers `testAdaptive` and `testLBFGS` with `OMP_NUM_THREADS=1`, `4`, and `8`.

Other test-like fixture directories are present under `test/`, including `testAC`, `testReweight`, `testSparseSampling`, and `testTS`, but they are not part of the official CTest command in this first-pass cleanup.

## Working Practices

- Inspect `git status --short` before editing.
- Treat existing uncommitted changes as intentional user work.
- Avoid editing generated build directories such as `build/`, `RELEASE/`, and `DEBUG/`.
- Keep reference-output updates separate from source changes unless the task explicitly updates tests.
- Avoid broad formatting churn.
- Prefer narrow edits that can be verified with the build and CTest commands above.

## Test Data Notes

The current official CTest suite passes `test/testdata/ModelPotential1d` as the data root used by `testAdaptive` and `testLBFGS`. Golden outputs live beside the corresponding test inputs, for example under `test/testAdaptive/` and `test/testLBFGS/`.

When changing behavior that affects numerical output, do not update reference files silently. Explain why the output changed and verify the change with the relevant test command.
