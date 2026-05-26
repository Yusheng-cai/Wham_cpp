# WHAM Agent and Test Cleanup Design

## Purpose

Make the repository easier for coding agents to work in without changing source behavior or expanding test coverage in the first pass.

The first pass should provide a concise operational guide and clarify the current official test path. It should not explain WHAM theory, refactor production code, or wire dormant fixture directories into CTest.

## Current Context

The repository is a C++ implementation of WHAM with a top-level executable, shared libraries, bundled dependencies, CMake build files, and shell-driven golden-output tests.

Relevant layout:

- `src/`: WHAM implementation and calculation strategies.
- `tools/`: parsing, utility, and helper code.
- `parallel/`: OpenMP and MPI-related helper code.
- `test/`: CTest definitions, test runner script, fixture inputs, and reference outputs.
- `scripts/`: example input files.
- `Eigen/` and `LBFGS/`: vendored dependencies.

The current CTest suite runs `testAdaptive` and `testLBFGS` with OpenMP thread counts `1`, `4`, and `8`. Other test-like fixture directories exist, but they are not part of the official first-pass test command.

## Design

Add a top-level `AGENTS.md` as a short operational guide for agents. It should answer the practical questions an agent needs before editing:

- What is this repository for?
- Which directories matter?
- How do I configure and build it?
- How do I run the current official tests?
- Which generated directories should I avoid editing?
- What repository-specific cautions should I follow?

The guide should remain concise. It should link or point to `README.md` for scientific context instead of duplicating WHAM background.

## Build and Test Guidance

Document the canonical first-pass build flow:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

Document that this test command currently covers only the registered CTest suite. Do not imply that all fixture directories under `test/` are active.

## Test Streamlining Scope

For this pass, preserve test behavior.

Allowed:

- Document the official test command in `AGENTS.md`.
- Clarify in `AGENTS.md` that extra fixture directories are present but not part of the official CTest suite.
- Optionally add a small comment near `MAIN_TESTS` in `test/CMakeLists.txt` explaining that it defines the current official suite.

Out of scope:

- Registering `testAC`, `testReweight`, `testSparseSampling`, or `testTS` with CTest.
- Updating golden reference outputs.
- Reworking `test/run_test.sh`.
- Changing production C++ behavior.
- Broad formatting or CMake modernization.

## Agent Cautions

The operational guide should tell agents to:

- Inspect `git status --short` before edits.
- Treat existing uncommitted changes as intentional user work.
- Avoid generated build directories such as `build/`, `RELEASE/`, and `DEBUG/`.
- Keep reference-output updates separate from code changes unless explicitly updating tests.
- Avoid broad formatting churn.
- Prefer narrow, locally verified changes.

## Verification

After implementation, verify with:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

If the existing dirty working tree prevents clean verification, report exactly what was run and what remains unverified.

## Open Decisions

No open decisions remain for the first pass. Future work can separately decide whether dormant fixture directories should become official CTest tests.

## Follow-Up Note

During implementation, the official CTest wiring was found to be stale. The narrow test cleanup also needs to use `BASH_PROGRAM`, point `testAdaptive` and `testLBFGS` at `test/testdata/ModelPotential1d`, and compare `kl.out` against `klref.out`.
