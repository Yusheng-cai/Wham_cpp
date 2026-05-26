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

- Latest full command: `ctest --test-dir build --output-on-failure`.
- Result: 1 of 7 tests passed.
- Passing test: `testWhamTools`.
- Failing tests: `testAdaptive::OMP_1`, `testLBFGS::OMP_1`, `testAdaptive::OMP_4`, `testLBFGS::OMP_4`, `testAdaptive::OMP_8`, and `testLBFGS::OMP_8`.
- `testAdaptive` variants run the executable successfully and fail against `pji.out` golden references. A direct focused comparison after the KL accumulation fix showed `norm.out` matches, while `pji.out` and `kl.out` differ from the checked-in references.
- `testLBFGS::OMP_1` aborts with `the line search routine reached the maximum number of iterations` and exits through the runner with status `134`.
- `testLBFGS::OMP_4` and `testLBFGS::OMP_8` run to output comparison and fail against `pji.out` golden references.

## Timing Baseline

- Working directory: `build/test`
- Command: `/usr/bin/time -v ../bin/Wham ../../test/testAdaptive/input.dat -abspath ../../test/testdata/ModelPotential1d`
- Exit status: `0`
- Elapsed wall time: `30:32.51`
- Maximum resident set size: `70120 kB`
- Caveat: this timing streamed verbose solver output through the terminal, so wall time is inflated relative to the CTest adaptive runtime. Future profiling runs should redirect stdout and stderr from the start.

## Next Optimization Targets

- Extract `WhamTools` declarations from `src/Wham.h` into `src/WhamTools.h`.
- Move `WhamTools` definitions from `src/Wham.cpp` into `src/WhamTools.cpp`.
- Remove repeated temporary vector allocation in `WhamTools::calculatelnWi`, `WhamTools::Gradient`, and `Uwham::calculate`.
- Gate or remove routine `std::cout` progress output.
