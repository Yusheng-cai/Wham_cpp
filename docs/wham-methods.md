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

## Algorithm Overview

The implementation follows the same high-level flow for both UWHAM and BWHAM:

```text
time-series samples + bias definitions
        |
        v
evaluate each sample or bin under every bias potential
        |
        v
solve for WHAM normalization constants f_k
        |
        v
compute unbiased sample/bin probabilities
        |
        v
write requested free-energy, KL, derivative, or reweighting outputs
```

Notation used below:

- `k = 1, ..., K` indexes simulation windows.
- `i = 1, ..., M` indexes unbinned samples.
- `l = 1, ..., L` indexes histogram bins.
- `N_k` is the number of samples from window `k`.
- `N_tot = sum_k N_k`.
- `B_{k,i} = beta U_k(x_i)` is the reduced bias energy of sample `i` under
  window `k`.
- `B_{k,l}` is the reduced bias energy of bin `l` under window `k`.
- `f_k` is the reduced free-energy normalization constant for window `k`.
- `M_l` is the number of samples in bin `l`.

All logarithmic sums are evaluated with a log-sum-exp form for numerical
stability.

## UWHAM Equations

For unbinned WHAM, the code solves for `f_k` using the negative log-likelihood:

```text
A(f) = - sum_k N_k f_k
       + sum_i log[ sum_k (N_k / N_tot) exp(f_k - B_{k,i}) ]
```

The sample log weight is:

```text
ln w_i = -log[ sum_k N_k exp(f_k - B_{k,i}) ]
```

After solving for `f_k`, the implementation normalizes the sample weights:

```text
ln w_i <- ln w_i - log[ sum_i exp(ln w_i) ]
```

and shifts the `f_k` values by the same normalization so the output constants
are consistent with the normalized weights.

The gradient used by both solvers can be written as:

```text
p_{k,i} = exp(f_k - B_{k,i} + ln w_i)

grad_k = -(N_k - N_k sum_i p_{k,i})
```

At convergence, the self-consistency condition is:

```text
sum_i p_{k,i} = 1
```

for every window with nonzero samples.

### L-BFGS Solver

The `LBFGS` strategy minimizes `A(f)` directly using the vendored L-BFGS
implementation. The solver receives:

- objective value: `A(f)`
- gradient: `grad_k`
- options: `max_iterations`, `epsilon`, `epsilon_rel`, and `printevery`

Because adding a constant to all `f_k` values does not change the WHAM solution,
the implementation fixes the gauge by subtracting `f_0` during objective
evaluation and again normalizes the final constants after optimization.

### Adaptive Solver

The `adaptive` strategy alternates between two candidate updates and chooses the
one with the smaller gradient norm.

The self-consistent candidate is:

```text
f_k^SC = -log[ sum_i exp(ln w_i - B_{k,i}) ]
```

then shifted so `f_0^SC = 0`.

The Newton-Raphson candidate is:

```text
f^NR = f - H^+ grad
```

where `H^+ grad` is computed by solving a least-squares system using the WHAM
Hessian. In code this is:

```text
(H^T H) delta = H^T grad
f^NR = f - delta
```

The Hessian terms are computed from `p_{k,i}`:

```text
H_{a,a} = N_a sum_i p_{a,i} - N_a^2 sum_i p_{a,i}^2

H_{a,b} = -N_a N_b sum_i p_{a,i} p_{b,i},  a != b
```

At each iteration, the adaptive solver computes both candidates, evaluates the
gradient norm for each, keeps the candidate with the smaller norm, and stops
when the relative change in `f_k` is below `tolerance`.

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

For binned WHAM, the code solves for `f_k` with the binned objective:

```text
A_b(f) = - sum_k N_k f_k
         + sum_l [ -M_l log(M_l)
                   + M_l log( sum_k N_k exp(f_k - B_{k,l}) ) ]
```

Bins with `M_l = 0` are skipped in the objective. After solving for `f_k`, the
binned log probability is:

```text
ln p_l = log(M_l) - log[ sum_k N_k exp(f_k - B_{k,l}) ]
```

The binned gradient is:

```text
grad_k = N_k ( exp(f_k + log[sum_l exp(ln p_l - B_{k,l})]) - 1 )
```

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
