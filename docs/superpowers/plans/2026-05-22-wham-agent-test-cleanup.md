# WHAM Agent and Test Cleanup Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a concise operational guide for agents and clarify the current official CTest suite without changing source behavior or expanding test coverage.

**Architecture:** This is a documentation-first cleanup. The top-level `AGENTS.md` becomes the operational entry point for agents, while `test/CMakeLists.txt` keeps the current CTest behavior and receives only a clarifying comment near the `MAIN_TESTS` list.

**Tech Stack:** CMake, CTest, Bash test runner, C++14 project layout.

---

## File Structure

- Create: `AGENTS.md`
  - Responsibility: Short operational guide for agents working in this repository.
- Modify: `test/CMakeLists.txt`
  - Responsibility: Keep existing registered CTest tests unchanged, with one comment clarifying that `MAIN_TESTS` is the current official suite.
- Reference only: `docs/superpowers/specs/2026-05-22-wham-agent-test-cleanup-design.md`
  - Responsibility: Approved design spec for this implementation pass.

## Task 1: Add Top-Level Agent Guide

**Files:**
- Create: `AGENTS.md`

- [ ] **Step 1: Confirm the guide does not already exist**

Run:

```bash
test ! -e AGENTS.md
```

Expected: exit code `0`. If this fails because `AGENTS.md` already exists, stop and inspect the file before editing.

- [ ] **Step 2: Create `AGENTS.md`**

Create `AGENTS.md` with exactly this content:

```markdown
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
```

- [ ] **Step 3: Verify the guide content**

Run:

```bash
sed -n '1,220p' AGENTS.md
```

Expected: the output matches the content from Step 2, including the build and test commands.

## Task 2: Clarify Current CTest Scope

**Files:**
- Modify: `test/CMakeLists.txt`

- [ ] **Step 1: Inspect the current `MAIN_TESTS` block**

Run:

```bash
sed -n '1,80p' test/CMakeLists.txt
```

Expected: the file contains:

```cmake
set(OMP_THREADS 1 4 8)
set(MAIN_TESTS "")
list(APPEND MAIN_TESTS "testAdaptive")
list(APPEND MAIN_TESTS "testLBFGS")
```

- [ ] **Step 2: Add one clarifying comment**

Change this block:

```cmake
set(OMP_THREADS 1 4 8)
set(MAIN_TESTS "")
list(APPEND MAIN_TESTS "testAdaptive")
list(APPEND MAIN_TESTS "testLBFGS")
```

to this block:

```cmake
set(OMP_THREADS 1 4 8)

# Current official CTest suite. Other fixture directories under test/ are
# intentionally not registered until their expected behavior is reviewed.
set(MAIN_TESTS "")
list(APPEND MAIN_TESTS "testAdaptive")
list(APPEND MAIN_TESTS "testLBFGS")
```

Do not add or remove tests.

- [ ] **Step 3: Verify only the intended CMake section changed**

Run:

```bash
git diff -- test/CMakeLists.txt
```

Expected diff:

```diff
@@
 set(OMP_THREADS 1 4 8)
+
+# Current official CTest suite. Other fixture directories under test/ are
+# intentionally not registered until their expected behavior is reviewed.
 set(MAIN_TESTS "")
 list(APPEND MAIN_TESTS "testAdaptive")
 list(APPEND MAIN_TESTS "testLBFGS")
```

## Task 3: Verify Build and Tests

**Files:**
- Read: `AGENTS.md`
- Read: `test/CMakeLists.txt`
- Build output: `build/`

- [ ] **Step 1: Configure**

Run:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
```

Expected: exit code `0`. If configuration fails because of pre-existing local source changes, capture the failure and do not broaden this cleanup task.

- [ ] **Step 2: Build**

Run:

```bash
cmake --build build
```

Expected: exit code `0`. If the build fails in source files unrelated to `AGENTS.md` or the CMake comment, report that verification is blocked by the existing workspace state.

- [ ] **Step 3: Run CTest**

Run:

```bash
ctest --test-dir build --output-on-failure
```

Expected: CTest runs the registered suite and exits with code `0`. The expected registered tests are:

```text
testAdaptive::OMP_1
testLBFGS::OMP_1
testAdaptive::OMP_4
testLBFGS::OMP_4
testAdaptive::OMP_8
testLBFGS::OMP_8
```

If test names are listed in a different order, that is acceptable. Do not add dormant fixtures to make this step pass.

## Task 4: Review Final Diff

**Files:**
- Review: `AGENTS.md`
- Review: `test/CMakeLists.txt`

- [ ] **Step 1: Check repo status**

Run:

```bash
git status --short
```

Expected: this task's intended changes are:

```text
?? AGENTS.md
 M test/CMakeLists.txt
```

Other pre-existing changes may still be present. Do not revert or stage unrelated user work.

- [ ] **Step 2: Review the intended diff**

Run:

```bash
git diff -- AGENTS.md test/CMakeLists.txt
```

Expected: the diff contains the full new `AGENTS.md` and only the comment addition in `test/CMakeLists.txt`.

- [ ] **Step 3: Commit if the repository permits it**

Run:

```bash
git add AGENTS.md test/CMakeLists.txt
git commit -m "Add agent operations guide"
```

Expected: commit succeeds if `.git` is writable. If Git reports a read-only filesystem or cannot create `.git/index.lock`, leave the files uncommitted and report the exact error.

## Follow-Up Scope Added During Execution

The official CTest command was blocked by stale test wiring. The approved narrow expansion is:

- Use `${BASH_PROGRAM}` instead of `${BASH}` in `test/CMakeLists.txt`.
- Use `test/testdata/ModelPotential1d` as the current official data root for `testAdaptive` and `testLBFGS`.
- Compare `kl.out` to `klref.out` instead of stale `h.out`/`href.out` names.
- Fix the `test/run_test.sh` shebang.
- Remove stale generated outputs before each test run and stop if the `Wham` executable exits nonzero.
