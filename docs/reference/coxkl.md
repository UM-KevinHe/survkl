# Cox Proportional Hazards Model with KL Divergence for Data Integration

Fits a Cox proportional hazards model that incorporates external
information via a Kullback–Leibler (KL) divergence penalty. External
information can be supplied either as external risk scores (`RS`) or as
external coefficients (`beta`). The tuning parameter(s) `etas` control
the strength of integration.

## Usage

``` r
coxkl(
  z,
  delta,
  time,
  stratum = NULL,
  RS = NULL,
  beta = NULL,
  etas,
  tol = 1e-04,
  Mstop = 100,
  backtrack = FALSE,
  message = FALSE,
  data_sorted = FALSE,
  beta_initial = NULL
)
```

## Arguments

- z:

  Numeric matrix of covariates with rows representing observations and
  columns representing predictor variables. All covariates must be
  numeric.

- delta:

  Numeric vector of event indicators (1 = event, 0 = censored).

- time:

  Numeric vector of observed event or censoring times. No sorting
  required.

- stratum:

  Optional numeric or factor vector defining strata.

- RS:

  Optional numeric vector or matrix of external risk scores. Length (or
  number of rows) must equal the number of observations. If not
  supplied, `beta` must be provided.

- beta:

  Optional numeric vector of external coefficients (e.g., from prior
  studies). Length must equal the number of columns in `z`. Use zeros to
  represent covariates without external information. If not supplied,
  `RS` must be provided.

- etas:

  Numeric vector of tuning parameters controlling the reliance on
  external information. Larger values place more weight on the external
  source.

- tol:

  Convergence tolerance for the optimization algorithm. Default is
  `1e-4`.

- Mstop:

  Maximum number of iterations for the optimization algorithm. Default
  is `100`.

- backtrack:

  Logical; if `TRUE`, backtracking line search is applied during
  optimization. Default is `FALSE`.

- message:

  Logical; if `TRUE`, progress messages are printed during model
  fitting. Default is `FALSE`.

- data_sorted:

  Logical; if `TRUE`, input data are assumed to be already sorted by
  stratum and time. Default is `FALSE`.

- beta_initial:

  Optional numeric vector of length `p` giving the starting value for
  the first `eta`. If `NULL`, a zero vector is used.

## Value

An object of class `"coxkl"` containing:

- `eta`: the fitted \\\eta\\ sequence.

- `beta`: estimated coefficient matrix (\\p \times \|\eta\|\\).

- `linear.predictors`: matrix of linear predictors.

- `likelihood`: vector of partial likelihoods.

- `data`: a list containing the input data used in fitting (`z`, `time`,
  `delta`, `stratum`, `data_sorted`).

## Details

If `beta` is supplied (length `ncol(z)`), external risk scores are
computed internally as `RS = z %*% beta`. If `RS` is supplied, it is
used directly. Data are optionally sorted by `stratum` (or a single
stratum if `NULL`) and increasing `time` when `data_sorted = FALSE`.
Estimation proceeds over the sorted data, and the returned
`linear.predictors` are mapped back to the original order. Optimization
uses warm starts across the (ascending) `etas` grid and supports
backtracking line search when `backtrack = TRUE`.

Internally, the routine computes a stratum-wise adjusted event indicator
(`delta_tilde`) and maximizes a KL-regularized partial likelihood. The
current implementation fixes `lambda = 0` in the low-level optimizer and
exposes `etas` as the primary tuning control.

## Examples

``` r
data(Exampledata_lowdim)
#> Warning: data set 'Exampledata_lowdim' not found

train_dat_lowdim <- ExampleData_lowdim$train
beta_external_good_lowdim <- ExampleData_lowdim$beta_external_good
eta_list <- generate_eta(method = "exponential", n = 10, max_eta = 5)

model <- coxkl(z = train_dat_lowdim$z,
               delta = train_dat_lowdim$status,
               time = train_dat_lowdim$time,
               stratum = train_dat_lowdim$stratum,
               RS = NULL,
               beta = beta_external_good_lowdim,
               etas = c(0:5))
```
