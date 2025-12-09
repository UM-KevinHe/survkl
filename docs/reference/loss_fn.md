# Calculate the Log-Partial Likelihood for a Stratified Cox Model

Computes the stratified Cox partial log-likelihood for given covariates,
event indicators, times, and coefficients.

## Usage

``` r
loss_fn(z, delta, time, stratum, beta)
```

## Arguments

- z:

  A numeric matrix (or data frame coercible to matrix) of covariates.
  Each row is an observation and each column a predictor.

- delta:

  A numeric vector of event indicators (1 = event, 0 = censored).

- time:

  A numeric vector of observed times (event or censoring).

- stratum:

  An optional vector specifying the stratum for each observation
  (factor/character/numeric). If missing, a single-stratum model is
  assumed.

- beta:

  A numeric vector of regression coefficients with length equal to the
  number of columns in `z`.

## Value

A single numeric value giving the stratified Cox partial log-likelihood.

## Details

Inputs are internally sorted by `stratum` and `time`. The function
evaluates the stratified Cox partial log-likelihood using the supplied
`z`, `delta`, `beta`, and the stratum sizes.
