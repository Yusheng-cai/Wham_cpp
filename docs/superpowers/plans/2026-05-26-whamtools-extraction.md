# WhamTools Extraction Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Move the `WhamTools` numerical helper namespace out of `src/Wham.h` and `src/Wham.cpp` into focused `src/WhamTools.h` and `src/WhamTools.cpp` files without changing behavior.

**Architecture:** Keep `Wham` responsible for WHAM object state, input wiring, and output methods. Put free numerical helpers in `WhamTools.*`, and keep `Wham.h` including `WhamTools.h` initially so existing consumers do not need a broad include migration in this pass.

**Tech Stack:** C++14, CMake, CTest.

---

### Task 1: Add Red Include Test

**Files:**
- Modify: `test/testWhamTools/test_wham_tools.cpp`

- [x] Change the test include from `src/Wham.h` to `src/WhamTools.h`.
- [x] Run `cmake --build build --target testWhamTools`.
- [x] Expected: build fails because `src/WhamTools.h` does not exist yet.

### Task 2: Extract WhamTools

**Files:**
- Create: `src/WhamTools.h`
- Create: `src/WhamTools.cpp`
- Modify: `src/Wham.h`
- Modify: `src/Wham.cpp`
- Modify: `src/CMakeLists.txt`

- [x] Move all `namespace WhamTools` declarations from `src/Wham.h` to `src/WhamTools.h`.
- [x] Move all `WhamTools::...` definitions from `src/Wham.cpp` to `src/WhamTools.cpp`.
- [x] Include `WhamTools.h` from `Wham.h` to preserve existing downstream includes.
- [x] Add `WhamTools.cpp` and `WhamTools.h` to `src/CMakeLists.txt`.

### Task 3: Verify

**Files:**
- Test only.

- [x] Run `cmake --build build`.
- [x] Run `ctest --test-dir build -R testWhamTools --output-on-failure`.
- [x] Run `ctest --test-dir build --output-on-failure`.
- [x] Commit the extraction if all verification passes.
