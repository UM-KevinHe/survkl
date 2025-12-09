# Calculate Survival Probabilities

Computes individual survival probabilities from a fitted linear
predictor `z%*%beta` using a stratified Breslow-type baseline hazard
estimate.

## Usage

``` r
cal_surv_prob(z, delta, time, beta, stratum)
```

## Arguments

- z:

  A numeric matrix (or data frame coercible to matrix) of covariates.
  Each row is an observation and each column a predictor.

- delta:

  A numeric vector of event indicators (1 = event, 0 = censored).

- time:

  A numeric vector of observed times (event or censoring).

- beta:

  A numeric vector of regression coefficients with length equal to the
  number of columns in `z`.

- stratum:

  An optional vector specifying the stratum for each observation. If
  missing, a single-stratum model is assumed.

## Value

A numeric matrix of survival probabilities with `nrow(z)` rows and
`length(time)` columns. Rows correspond to observations; columns are in
the internal sorted order of `(stratum, time)` (i.e., not collapsed to
unique event times). Entry `S[i, j]` is the estimated survival
probability for subject `i` evaluated at the `j`-th sorted time point.

## Details

Inputs are internally sorted by `stratum` and `time`. Within each
stratum, a baseline hazard increment is computed as `delta/S0`, where
`S0` is the risk set sum returned by `ddloglik_S0`. The stratified
baseline cumulative hazard `Lambda0` is then formed by a cumulative sum
within stratum, and individual survival curves are computed as
`S(t) = exp(-Lambda0(t) * exp(z %*% beta))`.
