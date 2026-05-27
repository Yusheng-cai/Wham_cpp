# WHAM Methods

`Wham_cpp` implements weighted-histogram methods for estimating unbiased free
energies from biased simulation data. The executable is configured through
native `.dat` input files; see [input-format.md](input-format.md) for input
syntax and generation.

## Core Concepts

Each simulation window contributes:

- A `timeseries` block describing where to read sampled coordinates or order
  parameters.
- A `bias` block describing the biasing potential used for that window.

The `wham` block selects the estimator, optimization strategy, bins, and
outputs.

## UWHAM

`Uwham` is the unbinned WHAM path. It uses the raw time-series samples rather
than first reducing every window to a histogram. This is the main path used by
the current official validation tests.

Current UWHAM strategy configuration uses one or more `Uwhamstrategy` blocks:

```text
wham = {
    name = w
    type = Uwham
    Uwhamstrategy = {
        type = LBFGS
        name = l
        max_iterations = 50
        printevery = -1
    }
    Uwhamstrategy = {
        type = adaptive
        name = a
        printevery = -1
    }
    strategyNames = [ l a ]
}
```

`strategyNames` controls the order in which named strategies are applied.

### L-BFGS Strategy

The `LBFGS` strategy solves the UWHAM optimization problem using the vendored
L-BFGS implementation.

Common options:

- `name`: strategy name used by `strategyNames`.
- `max_iterations`: maximum optimizer iterations.
- `epsilon`: optimizer convergence tolerance.
- `epsilon_rel`: relative convergence tolerance.
- `printevery`: progress-print interval. Use `-1` to suppress routine output.

### Adaptive Strategy

The `adaptive` strategy performs iterative updates using the adaptive WHAM
implementation.

Common options:

- `name`: strategy name used by `strategyNames`.
- `tolerance`: convergence tolerance.
- `printevery`: progress-print interval. Use `-1` to suppress routine output.

## BWHAM

`Bwham` is the binned WHAM path. It first bins data according to the configured
`bins` block and then estimates the free energy from binned counts.

Example:

```text
wham = {
    type = Bwham
    strategy = LBFGS
    bins = {
        dimension = 1
        range = [ 0 30 ]
        numbins = 30
    }
    outputs = [ lnpl ]
    outputFile = [ p.out ]
}
```

`Bwham` currently uses the `LBFGS` calculation strategy.

## Bias Models

Bias blocks are factory-created by `type`. If `type` is omitted, the default is
`simplebias`.

### `simplebias`

`simplebias` combines a harmonic term and an optional linear term:

```text
bias = {
    type = simplebias
    dimension = 1
    xstar = [ 0.0 ]
    kappa = [ 1000 ]
    phi = [ 0.0 ]
    temperature = 300
}
```

Common keys:

- `dimension`: number of biased dimensions.
- `xstar`: bias center.
- `kappa`: harmonic force constant.
- `phi`: linear coefficient.
- `temperature`: temperature in Kelvin.

### `squaredbias`

`squaredbias` uses a quadratic form controlled by `phi`:

```text
bias = {
    type = squaredbias
    dimension = 1
    phi = [ 100 ]
}
```

## Binning

Each `bins` block describes one dimension:

```text
bins = {
    dimension = 1
    range = [ -1.5 1.5 ]
    numbins = 50
}
```

`range` is lower-inclusive and upper-exclusive. `dimension` is 1-based and
selects which coordinate dimension the bin applies to.

## Outputs

Outputs are selected by matching `outputs` and `outputFile` vectors:

```text
outputs = [ pji normalization KL_divergence ]
outputFile = [ pji.out norm.out kl.out ]
```

Common UWHAM outputs:

- `pji`
- `normalization`
- `lnwji`
- `derivative`
- `reweightFE`
- `KL_divergence`
- `FE_dim`
- `ErrorFE`

Common shared WHAM outputs:

- `histogram`
- `Autocorrelation`
- `forces`
- `Averages`
- `dataFE`

Common BWHAM output:

- `lnpl`

## Reweighting

Reweighting blocks act on a named WHAM calculation:

```text
Reweight = {
    wham = w
    type = UwhamReweight
    bias = {
        dimension = 1
        phi = [ 100 ]
    }
    outputs = [ ReweightAverages ]
    outputNames = [ ReweightAvg.out ]
}
```

The `wham` field must match a named `wham` block.
