# Cross-Validated Selection of Integration Parameter （`eta`） for the Cox–KL Model

Performs K-fold cross-validation to select the integration parameter
`eta` for the Cox–KL model. Each fold fits the model on a training split
and evaluates on the held-out split using the specified performance
criterion.

## Usage

``` r
cv.coxkl(
  z,
  delta,
  time,
  stratum = NULL,
  RS = NULL,
  beta = NULL,
  etas = NULL,
  tol = 1e-04,
  Mstop = 100,
  backtrack = FALSE,
  nfolds = 5,
  criteria = c("V&VH", "LinPred", "CIndex_pooled", "CIndex_foldaverage"),
  c_index_stratum = NULL,
  message = FALSE,
  seed = NULL,
  ...
)
```

## Arguments

- z:

  Numeric matrix of covariates (rows = observations, columns =
  variables).

- delta:

  Numeric vector of event indicators (1 = event, 0 = censored).

- time:

  Numeric vector of observed event or censoring times.

- stratum:

  Optional numeric or factor vector defining strata. If `NULL`, all
  observations are treated as a single stratum.

- RS:

  Optional numeric vector or matrix of external risk scores. If omitted,
  `beta` must be supplied.

- beta:

  Optional numeric vector of external coefficients. If omitted, `RS`
  must be supplied.

- etas:

  Numeric vector of candidate tuning values to be cross-validated.
  (required). Values are internally sorted in ascending order.

- tol:

  Convergence tolerance for the optimizer used inside
  [`coxkl`](https://um-kevinhe.github.io/survkl/reference/coxkl.md).
  Default `1e-4`.

- Mstop:

  Maximum number of Newton iterations used inside
  [`coxkl`](https://um-kevinhe.github.io/survkl/reference/coxkl.md).
  Default `100`.

- backtrack:

  Logical; if `TRUE`, backtracking line search is applied during
  optimization. Default is `FALSE`.

- nfolds:

  Number of cross-validation folds. Default `5`.

- criteria:

  Character string specifying the performance criterion. Choices are
  `"V&VH"` (default), `"LinPred"`, `"CIndex_pooled"`, or
  `"CIndex_foldaverage"`.

- c_index_stratum:

  Optional stratum vector. Only required when `criteria` is set to
  `"CIndex_pooled"` or `"CIndex_foldaverage"`, and a stratified C-index
  is desired while the fitted model is non-stratified. Default `NULL`.

- message:

  Logical; if `TRUE`, prints progress messages and per-fold progress
  bars. Default `FALSE`.

- seed:

  Optional integer seed for reproducible fold assignment. Default
  `NULL`.

- ...:

  Additional arguments passed to
  [`coxkl`](https://um-kevinhe.github.io/survkl/reference/coxkl.md).

## Value

An object of class `"cv.coxkl"` with components:

- `internal_stat`:

  A data.frame with one row per `eta` containing `eta` and the
  cross-validated measure named according to `criteria` (one of
  `VVH_Loss`, `LinPred_Loss`, `CIndex_pooled`, `CIndex_foldaverage`).

- `external_stat`:

  Scalar baseline statistic computed from `RS` under the same
  `criteria`.

- `criteria`:

  The evaluation criterion used.

- `nfolds`:

  Number of folds.

## Details

External information is required: supply either `RS` or `beta` (if
`beta` is given, `RS` is computed as `z %*% beta`). Folds are created
with stratification by `stratum` and censoring status. Within each fold
and each candidate `eta`, the function fits `coxkl` on the training
split with warm-starts initialized to zero and evaluates on the test
split:

- `"V&VH"`: uses the difference of partial log-likelihoods between full
  and training fits; reported as \\-2\\ times the aggregated quantity.

- `"LinPred"`: aggregates the test-split linear predictors across folds
  and evaluates \\-2\\ times the partial log-likelihood on the full
  data.

- `"CIndex_pooled"`: pools pairwise comparable counts across folds
  (numerator/denominator).

- `"CIndex_foldaverage"`: averages the per-fold stratified C-index.

The function also computes an external baseline statistic from `RS`
using the same criterion for comparison.

## Examples

``` r
data(ExampleData_lowdim)

train_dat_lowdim <- ExampleData_lowdim$train
beta_external_good_lowdim <- ExampleData_lowdim$beta_external_good

etas <- generate_eta(method = "exponential", n = 10, max_eta = 5)

cv_res <- cv.coxkl(z = train_dat_lowdim$z,
                   delta = train_dat_lowdim$status,
                   time = train_dat_lowdim$time,
                   beta = beta_external_good_lowdim,
                   etas = etas)
#> Warning: Stratum not provided. Treating all data as one stratum.
                    
```
