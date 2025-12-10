# Cross-Validation for CoxKL Ridge Model (eta tuning)

This function performs cross-validation on the Cox model with
Kullback–Leibler (KL) penalty and ridge (L2) regularization. It tunes
the parameter `eta` (external information weight) using user-specified
cross-validation criteria, while internally evaluating a `lambda` path
(provided or generated) and selecting the best `lambda` per `eta`

## Usage

``` r
cv.coxkl_ridge(
  z,
  delta,
  time,
  stratum = NULL,
  RS = NULL,
  beta = NULL,
  etas,
  lambda = NULL,
  nlambda = 100,
  lambda.min.ratio = ifelse(n_obs < n_vars, 0.01, 1e-04),
  nfolds = 5,
  cv.criteria = c("V&VH", "LinPred", "CIndex_pooled", "CIndex_foldaverage"),
  c_index_stratum = NULL,
  message = FALSE,
  seed = NULL,
  ...
)
```

## Arguments

- z:

  Numeric matrix of covariates with rows representing individuals and
  columns representing predictors.

- delta:

  Numeric vector of event indicators (1 = event, 0 = censored).

- time:

  Numeric vector of observed times (event or censoring).

- stratum:

  Optional factor or numeric vector indicating strata.

- RS:

  Optional numeric vector or matrix of external risk scores. If not
  provided, `beta` must be supplied.

- beta:

  Optional numeric vector of external coefficients (length equal to
  `ncol(z)`). If not provided, `RS` must be supplied.

- etas:

  Numeric vector of candidate `eta` values to be evaluated.

- lambda:

  Optional numeric scalar or vector of penalty parameters. If `NULL`, a
  sequence is generated automatically.

- nlambda:

  Integer number of lambda values to generate when `lambda` is `NULL`.
  Default `100`.

- lambda.min.ratio:

  Ratio of the smallest to the largest lambda when generating a sequence
  (when `lambda` is `NULL`). Default `0.01` when `n < p`, otherwise
  `1e-4`.

- nfolds:

  Integer; number of cross-validation folds. Default `5`.

- cv.criteria:

  Character string specifying the cross-validation criterion. Choices
  are:

  - `"V&VH"` (default): `"V&VH"` loss.

  - `"LinPred"`: loss based on cross-validated linear predictors.

  - `"CIndex_pooled"`: pool all held-out predictions and compute one
    overall C-index.

  - `"CIndex_foldaverage"`: average C-index across folds.

- c_index_stratum:

  Optional stratum vector. Used only when `cv.criteria` is
  `"CIndex_pooled"` or `"CIndex_foldaverage"` to compute a stratified
  C-index for a non-stratified fit; if supplied, it must be identical to
  `stratum`. Default `NULL`.

- message:

  Logical; whether to print progress messages. Default `FALSE`.

- seed:

  Optional integer random seed for fold assignment.

- ...:

  Additional arguments passed to
  [`coxkl_ridge`](https://um-kevinhe.github.io/survkl/reference/coxkl_ridge.md).

## Value

An object of class `"cv.coxkl_ridge"`:

- `integrated_stat.full_results`:

  Data frame with columns `eta`, `lambda`, and the aggregated CV score
  per `lambda`; for loss criteria an additional column
  `Loss = -2 * score`; for C-index criteria a column named
  `CIndex_pooled` or `CIndex_foldaverage`.

- `integrated_stat.best_per_eta`:

  Data frame with the best `lambda` (per `eta`) according to the chosen
  criterion.

- `external_stat`:

  Scalar baseline statistic computed from `RS` under the same
  `cv.criteria`.

- `criteria`:

  The evaluation criterion used.

- `nfolds`:

  Number of folds.

## Details

Data are sorted by `stratum` and `time`. External information must be
given via `RS` or `beta` (if `beta` has length `ncol(z)`, the function
computes `RS = z %*% beta`). For each candidate `eta`, a `lambda` path
is determined (generated if `lambda = NULL`, otherwise the supplied
`lambda` values are sorted decreasingly). Cross-validation folds are
created by `get_fold`. In each fold,
[`coxkl_ridge`](https://um-kevinhe.github.io/survkl/reference/coxkl_ridge.md)
is fit on the training split across the full `lambda` path with
`data_sorted = TRUE`, and the chosen criterion is evaluated on the test
split and aggregated:

- `"V&VH"`: sums `pl(full) - pl(train)` across folds (reported as loss
  via `Loss = -2 * score`).

- `"LinPred"`: aggregates test-fold linear predictors and evaluates
  partial log-likelihood on full data (reported as `Loss = -2 * score`).

- `"CIndex_pooled"`: pools comparable-pair numerators/denominators
  across folds to compute one C-index.

- `"CIndex_foldaverage"`: averages the per-fold stratified C-index.

The best `lambda` is chosen per `eta` (minimizing loss or maximizing
C-index). The function also computes an external baseline statistic from
`RS` under the same criterion.

## Examples

``` r
# \donttest{
data(ExampleData_highdim) 

train_dat_highdim <- ExampleData_highdim$train
beta_external_highdim <- ExampleData_highdim$beta_external

etas <- generate_eta(method = "exponential", n = 10, max_eta = 100)
etas <- sample(etas)
cv_res <- cv.coxkl_ridge(z = train_dat_highdim$z,
                         delta = train_dat_highdim$status,
                         time = train_dat_highdim$time,
                         stratum = NULL,
                         RS = NULL,
                         beta = beta_external_highdim,
                         etas = etas,
                         nfolds = 5, 
                         cv.criteria = "CIndex_pooled",
                         message = TRUE)
#> Warning: Stratum not provided. Treating all data as one stratum.
#> Cross-validation over eta sequence:
#>   |                                      |                              |   0%  |                                      |===                           |  10%  |                                      |======                        |  20%  |                                      |=========                     |  30%  |                                      |============                  |  40%  |                                      |===============               |  50%  |                                      |==================            |  60%  |                                      |=====================         |  70%  |                                      |========================      |  80%  |                                      |===========================   |  90%  |                                      |==============================| 100%
# }
```
