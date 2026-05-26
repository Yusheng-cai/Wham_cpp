# WHAM Core Optimization Sprint Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Create a reliable fast baseline for WHAM math behavior, fix the most suspicious correctness issues, and gather a timing baseline before larger optimization work.

**Architecture:** This sprint deliberately avoids broad data-layout rewrites. It adds a small C++ math test executable, fixes isolated numerical bugs in existing functions, keeps the current official CTest harness intact, and records profiling commands for the next sprint.

**Tech Stack:** C++14, CMake, CTest, Bash, OpenMP-enabled WHAM library.

---

## File Structure

- Create: `test/testWhamTools/test_wham_tools.cpp`
  - Responsibility: Fast executable tests for pure `WhamTools` numerical helpers.
- Modify: `test/CMakeLists.txt`
  - Responsibility: Build and register the fast `testWhamTools` CTest target without changing the existing long WHAM fixture tests.
- Modify: `src/Wham.cpp`
  - Responsibility: Fix isolated BAR convergence behavior.
- Modify: `src/Uwham.cpp`
  - Responsibility: Fix KL-divergence accumulation.
- Create: `docs/superpowers/specs/2026-05-23-wham-core-optimization-notes.md`
  - Responsibility: Record baseline test/profiling results and remaining optimization candidates.

## Task 1: Add Fast WhamTools Tests

**Files:**
- Create: `test/testWhamTools/test_wham_tools.cpp`
- Modify: `test/CMakeLists.txt`

- [ ] **Step 1: Create the test directory**

Run:

```bash
mkdir -p test/testWhamTools
```

Expected: directory exists at `test/testWhamTools`.

- [ ] **Step 2: Add the test source**

Create `test/testWhamTools/test_wham_tools.cpp` with this content:

```cpp
#include "src/Wham.h"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace
{
bool near(double actual, double expected, double tolerance)
{
    return std::abs(actual - expected) <= tolerance;
}

int require_near(const std::string& name, double actual, double expected, double tolerance)
{
    if (!near(actual, expected, tolerance))
    {
        std::cerr << name << " expected " << expected << " but got " << actual << "\n";
        return 1;
    }

    return 0;
}
}

int main()
{
    int failures = 0;

    failures += require_near(
        "LogSumExp",
        WhamTools::LogSumExp(std::vector<double>{0.0, 0.0}, std::vector<double>{1.0, 1.0}),
        std::log(2.0),
        1e-12);

    failures += require_near(
        "EXP",
        WhamTools::EXP(std::vector<double>{0.0, 0.0}),
        0.0,
        1e-12);

    failures += require_near(
        "CalculateDeltaFBarIterative",
        WhamTools::CalculateDeltaFBarIterative(
            std::vector<double>{0.2, 0.4, 0.6},
            std::vector<double>{-0.1, 0.1, 0.3},
            500,
            1e-12),
        0.15,
        1e-8);

    return failures == 0 ? 0 : 1;
}
```

- [ ] **Step 3: Register the test target**

In `test/CMakeLists.txt`, add this block before `add_custom_target(build_test ...)`:

```cmake
add_executable(testWhamTools testWhamTools/test_wham_tools.cpp)
target_include_directories(testWhamTools PUBLIC ${CMAKE_SOURCE_DIR})
target_link_libraries(testWhamTools PUBLIC WHAMsrc WHAMtools WHAMparallel)
add_test(NAME testWhamTools COMMAND testWhamTools)
```

- [ ] **Step 4: Reconfigure**

Run:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
```

Expected: exit code `0`.

- [ ] **Step 5: Build the fast test target**

Run:

```bash
cmake --build build --target testWhamTools
```

Expected: exit code `0`.

- [ ] **Step 6: Verify the new test fails for the BAR bug**

Run:

```bash
ctest --test-dir build -R testWhamTools --output-on-failure
```

Expected: exit code nonzero. Output should include:

```text
CalculateDeltaFBarIterative expected 0.15
```

The current implementation returns after one iteration because it computes the relative change after overwriting the previous value.

## Task 2: Fix BAR Iterative Convergence

**Files:**
- Modify: `src/Wham.cpp`
- Test: `test/testWhamTools/test_wham_tools.cpp`

- [ ] **Step 1: Replace the BAR iteration update**

In `src/Wham.cpp`, replace `WhamTools::CalculateDeltaFBarIterative` with:

```cpp
WhamTools::Real WhamTools::CalculateDeltaFBarIterative(const std::vector<Real>& w_F, const std::vector<Real>& w_B, int max_iterations, Real tol)
{
    Real DeltaF = 0.0;

    for (int i=0;i<max_iterations;i++)
    {
        Real DeltaFold = DeltaF;
        DeltaF = DeltaFold - CalculateBAR(w_F, w_B, DeltaFold);

        Real scale = std::max(std::abs(DeltaFold), 1.0);
        Real relativeChange = std::abs(DeltaF - DeltaFold) / scale;

        if (relativeChange < tol)
        {
            break;
        }
    }

    return DeltaF;
}
```

- [ ] **Step 2: Run the fast test**

Run:

```bash
cmake --build build --target testWhamTools
ctest --test-dir build -R testWhamTools --output-on-failure
```

Expected: `testWhamTools` passes.

- [ ] **Step 3: Run one official adaptive test**

Run:

```bash
ctest --test-dir build -R 'testAdaptive::OMP_1' --output-on-failure
```

Expected: it may still fail against stale golden references, but it should run the executable and compare outputs rather than fail in the harness.

## Task 3: Fix KL Accumulation

**Files:**
- Modify: `src/Uwham.cpp`
- Test: official adaptive CTest fixture output.

- [ ] **Step 1: Change assignment to accumulation**

In `src/Uwham.cpp`, inside `Uwham::calculate()`, replace:

```cpp
KL_divergence_[i] = prob * (-ref_val + val);
```

with:

```cpp
KL_divergence_[i] += prob * (-ref_val + val);
```

- [ ] **Step 2: Build**

Run:

```bash
cmake --build build
```

Expected: exit code `0`.

- [ ] **Step 3: Run focused tests**

Run:

```bash
ctest --test-dir build -R 'testWhamTools|testAdaptive::OMP_1' --output-on-failure
```

Expected:

- `testWhamTools` passes.
- `testAdaptive::OMP_1` may still fail if the golden files represent older stochastic bootstrap output. Record whether `kl.out` now differs in direction or magnitude from `klref.out`.

## Task 4: Document Remaining Test Baseline State

**Files:**
- Create: `docs/superpowers/specs/2026-05-23-wham-core-optimization-notes.md`

- [ ] **Step 1: Run the current official tests once**

Run:

```bash
ctest --test-dir build --output-on-failure
```

Expected: collect the final state. If tests fail, do not update golden references in this sprint.

- [ ] **Step 2: Write the optimization notes**

Create `docs/superpowers/specs/2026-05-23-wham-core-optimization-notes.md` with this structure:

```markdown
# WHAM Core Optimization Notes

## Baseline Commands

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

## Fast Test Coverage

- `testWhamTools` covers `LogSumExp`, `EXP`, and `CalculateDeltaFBarIterative`.

## Current Official Test Status

- Record exact pass/fail status from the latest `ctest --test-dir build --output-on-failure` run.
- If `testAdaptive` fails only by golden-output differences, record which output file differs.
- If `testLBFGS` still aborts, record the exception text.

## Next Optimization Targets

- Extract `WhamTools` declarations from `src/Wham.h` into `src/WhamTools.h`.
- Move `WhamTools` definitions from `src/Wham.cpp` into `src/WhamTools.cpp`.
- Remove repeated temporary vector allocation in `WhamTools::calculatelnWi`, `WhamTools::Gradient`, and `Uwham::calculate`.
- Gate or remove routine `std::cout` progress output.
```

Replace the bullet text under “Current Official Test Status” with the actual result observed in Step 1.

## Task 5: Record Timing Baseline

**Files:**
- Modify: `docs/superpowers/specs/2026-05-23-wham-core-optimization-notes.md`

- [ ] **Step 1: Time one representative adaptive run**

Run:

```bash
/usr/bin/time -v build/bin/Wham test/testAdaptive/input.dat -abspath test/testdata/ModelPotential1d
```

Expected: command completes or fails with the same behavior as the corresponding CTest run. Capture elapsed time, max resident set size, and exit status.

- [ ] **Step 2: Append timing results**

Add this section to `docs/superpowers/specs/2026-05-23-wham-core-optimization-notes.md` using the actual values from Step 1:

```markdown
## Timing Baseline

- Command: `/usr/bin/time -v build/bin/Wham test/testAdaptive/input.dat -abspath test/testdata/ModelPotential1d`
- Exit status: record the integer exit status reported by `/usr/bin/time`.
- Elapsed wall time: record the `Elapsed (wall clock) time` line reported by `/usr/bin/time`.
- Maximum resident set size: record the `Maximum resident set size` line reported by `/usr/bin/time`.
```

## Task 6: Review Diff and Stop Before Larger Refactors

**Files:**
- Review: `test/testWhamTools/test_wham_tools.cpp`
- Review: `test/CMakeLists.txt`
- Review: `src/Wham.cpp`
- Review: `src/Uwham.cpp`
- Review: `docs/superpowers/specs/2026-05-23-wham-core-optimization-notes.md`

- [ ] **Step 1: Review intended diff**

Run:

```bash
git diff -- test/testWhamTools/test_wham_tools.cpp test/CMakeLists.txt src/Wham.cpp src/Uwham.cpp docs/superpowers/specs/2026-05-23-wham-core-optimization-notes.md
```

Expected: diff includes only the fast test target, the BAR convergence fix, the KL accumulation fix, and the notes file.

- [ ] **Step 2: Check status**

Run:

```bash
git status --short
```

Expected: unrelated pre-existing local changes may remain. Do not revert or stage unrelated user work.

- [ ] **Step 3: Stop and ask before extraction/data-layout work**

Do not extract `WhamTools` or change data layout in this sprint. Those are next-sprint changes after this correctness and timing baseline is documented.
