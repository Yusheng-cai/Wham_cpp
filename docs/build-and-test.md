# Build and Test

This project uses CMake to build the `Wham` executable and CTest to run the
registered validation suite.

## Requirements

- CMake 3.18 or newer
- A C++14 compiler
- OpenMP
- FFTW3
- Bash
- Python 3

On Ubuntu:

```bash
sudo apt-get update
sudo apt-get install --yes cmake g++ libfftw3-dev python3
```

## Configure

Use an out-of-source build directory:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
```

The build writes executables to `build/bin/` and libraries to `build/lib/`.

## Build

```bash
cmake --build build --parallel
```

The main executable is:

```bash
build/bin/Wham
```

## Run Tests

Run the full registered CTest suite:

```bash
ctest --test-dir build --output-on-failure
```

The current suite includes:

- `testAdaptive::OMP_1`
- `testLBFGS::OMP_1`
- `testAdaptive::OMP_4`
- `testLBFGS::OMP_4`
- `testAdaptive::OMP_8`
- `testLBFGS::OMP_8`
- `testWhamTools`
- `testCoreUtilities`
- `testGenerateWhamInput`

The adaptive and L-BFGS executable tests compare generated outputs against
reference files under `test/testAdaptive/` and `test/testLBFGS/`. The
`OMP_*` suffix indicates the `OMP_NUM_THREADS` value used by the test runner.

## Run a Single Test

Use `-R` with a CTest name pattern:

```bash
ctest --test-dir build -R testCoreUtilities --output-on-failure
```

## Continuous Integration

GitHub Actions runs the same basic workflow on pushes and pull requests:

1. Install system dependencies.
2. Configure CMake in Release mode.
3. Build the project.
4. Run the registered CTest suite.

The workflow definition is `.github/workflows/build.yml`.

## Notes for Developers

- Keep builds out of the source tree. Use `build/` or another ignored build
  directory.
- Do not update reference outputs silently. If a numerical change is expected,
  document why the output changed and verify the relevant tests.
- Generated build directories such as `build/`, `RELEASE/`, and `DEBUG/`
  should not be edited or committed.
