# WHAM Input Generation

`Wham` reads native `.dat` files made from repeated `timeseries`, `bias`, and
`wham` blocks. For common umbrella-sampling workflows, writing those blocks by
hand is repetitive. The repository includes a small generator that converts a
compact JSON specification into the native input format.

## Quick Start

Generate a native input file from the bundled 1D UWHAM example:

```bash
python3 tools/generate_wham_input.py examples/input_specs/uwham_1d.json -o input.dat
```

Run WHAM with the generated file:

```bash
build/bin/Wham input.dat -abspath test/testdata/ModelPotential1d
```

The generator is dependency-free and uses Python's standard `json` module.

## Structured JSON Spec

The compact spec has two required top-level objects:

- `data`: time-series defaults and per-window bias parameters.
- `wham`: WHAM method, strategies, bins, and outputs.

Minimal 1D UWHAM example:

```json
{
  "data": {
    "columns": [2],
    "skipfrombeginning": -1000,
    "windows": [
      {
        "path": "US_-0.100.dat",
        "xstar": [-0.1],
        "kappa": [1000],
        "temperature": 300
      },
      {
        "path": "US_0.000.dat",
        "xstar": [0.0],
        "kappa": [1000],
        "temperature": 300
      }
    ]
  },
  "wham": {
    "name": "w",
    "type": "Uwham",
    "strategies": [
      {
        "type": "LBFGS",
        "name": "l",
        "max_iterations": 50,
        "printevery": -1
      },
      {
        "type": "adaptive",
        "name": "a",
        "printevery": -1
      }
    ],
    "bins": [
      {
        "dimension": 1,
        "range": [-1.5, 1.5],
        "numbins": 50
      }
    ],
    "outputs": {
      "pji": "pji.out",
      "normalization": "norm.out",
      "KL_divergence": "kl.out"
    }
  }
}
```

## Data Section

`data.columns` is the default list of 1-based columns to read from each
time-series file. Each item in `data.windows` must provide `path` and bias
parameters. The generator emits one native `timeseries` block and one native
`bias` block per window.

Common time-series keys:

- `path`
- `columns`
- `skipfrombeginning`
- `skip`
- `outputs`
- `outputNames`
- `verbose`

Common bias keys:

- `type`
- `dimension`
- `xstar`
- `kappa`
- `phi`
- `temperature`

If `dimension` is omitted, the generator infers it from `xstar`, `kappa`, `phi`,
or `data.columns`.

## WHAM Section

`wham.type` is required. For `Uwham`, use `wham.strategies` to emit modern
`Uwhamstrategy` blocks. If `wham.strategyNames` is omitted, the generator uses
the strategy names in list order.

`wham.bins` is a list so multidimensional inputs can emit one native `bins`
block per dimension.

`wham.outputs` can be written as an object mapping output names to files:

```json
"outputs": {
  "pji": "pji.out",
  "normalization": "norm.out"
}
```

The generator emits matching native `outputs` and `outputFile` vectors.

## Template

For users who prefer to edit native input directly, start from:

```text
examples/templates/uwham_1d.dat
```

Use the JSON generator for reproducible workflows and the native template for
quick manual edits.
