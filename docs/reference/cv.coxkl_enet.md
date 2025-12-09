# Cross-Validation for CoxKL Model with elastic net & lasso penalty

This function performs cross-validation on the high-dimensional Cox
model with Kullback–Leibler (KL) penalty. It tunes the parameter `eta`
(external information weight) using user-specified cross-validation
criteria, while also evaluating a `lambda` path (either provided or
generated) and selecting the best `lambda` per `eta`.

## Usage

``` r
cv.coxkl_enet(
  z,
  delta,
  time,
  stratum = NULL,
  RS = NULL,
  beta = NULL,
  etas,
  alpha = 1,
  lambda = NULL,
  nlambda = 100,
  lambda.min.ratio = ifelse(n < p, 0.05, 0.001),
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

- alpha:

  Elastic-net mixing parameter in \\(0,1\]\\. Default = `1` (lasso
  penalty).

- lambda:

  Optional numeric scalar or vector of penalty parameters. If `NULL`, a
  decreasing path is generated using `nlambda` and `lambda.min.ratio`.

- nlambda:

  Integer number of lambda values to generate when `lambda` is `NULL`.
  Default `100`.

- lambda.min.ratio:

  Ratio of the smallest to the largest lambda when generating a sequence
  (when `lambda` is `NULL`). Default `0.05` when `n < p`, otherwise
  `1e-3`.

- nfolds:

  Integer; number of cross-validation folds. Default = `5`.

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
  C-index while the fitted model is non-stratified; if supplied, it must
  be identical to `stratum`. Default `NULL`.

- message:

  Logical; whether to print progress messages. Default = `FALSE`.

- seed:

  Optional integer random seed for fold assignment.

- ...:

  Additional arguments passed to
  [`coxkl_enet`](https://umkevinhe.github.io/survkl/reference/coxkl_enet.md).

## Value

An object of class `"cv.coxkl_enet"`:

- `integrated_stat.full_results`:

  Data frame with columns `eta`, `lambda`, and the aggregated CV score
  for each `lambda` under the chosen `cv.criteria`. For loss criteria,
  an additional column with the transformed loss (`Loss = -2 * score`);
  for C-index criteria, a column named `CIndex_pooled` or
  `CIndex_foldaverage`.

- `integrated_stat.best_per_eta`:

  Data frame with the best `lambda` (per `eta`) according to the chosen
  `cv.criteria` (minimizing loss or maximizing C-index).

- `integrated_stat.betahat_best`:

  Matrix of coefficient vectors (columns) corresponding to the best
  `lambda` for each `eta`.

- `external_stat`:

  Scalar baseline statistic computed from the external risk score `RS`
  under the same `cv.criteria`.

- `criteria`:

  The evaluation criterion used (as provided in `cv.criteria`).

- `alpha`:

  The elastic-net mixing parameter used.

- `nfolds`:

  Number of folds.

## Details

Data are sorted by `stratum` and `time`. External info must be from `RS`
or `beta` (if `beta` given with length `ncol(z)`, `RS = z %*% beta`);
`alpha` \\\in (0,1\]\\.

For each candidate `eta`, a decreasing `lambda` path is used (generated
from `nlambda`/`lambda.min.ratio` if `lambda = NULL`); CV folds are
created by `get_fold`. Each fold fits
[`coxkl_enet`](https://umkevinhe.github.io/survkl/reference/coxkl_enet.md)
on the training split (full `lambda` path) and evaluates the chosen
criterion on the test split.

Aggregation follows the code paths for `"V&VH"`, `"LinPred"`,
`"CIndex_pooled"`, or `"CIndex_foldaverage"`:

- `"V&VH"`: sums `pl(full) - pl(train)` across folds (reported as loss
  via `Loss = -2 * score`).

- `"LinPred"`: aggregates test-fold linear predictors and evaluates
  partial log-likelihood on full data (reported as `Loss = -2 * score`).

- `"CIndex_pooled"`: pools comparable-pair numerators/denominators
  across folds to compute one C-index.

- `"CIndex_foldaverage"`: averages the per-fold stratified C-index.

The best `lambda` is selected per `eta` (min loss / max C-index), and
the function returns full results, the per-`eta` optimum, corresponding
coefficients, and an external baseline from `RS`.

## Examples

``` r
data(ExampleData_highdim) 

train_dat_highdim <- ExampleData_highdim$train
beta_external_highdim <- ExampleData_highdim$beta_external

etas <- generate_eta(method = "exponential", n = 10, max_eta = 100)
etas <- sample(etas) 

cv_res <- cv.coxkl_enet(z = train_dat_highdim$z,
                        delta = train_dat_highdim$status,
                        time = train_dat_highdim$time,
                        stratum = NULL,
                        RS = NULL,
                        beta = beta_external_highdim,
                        etas = etas,
                        alpha = 1.0,
                        nfolds = 5, 
                        cv.criteria = "CIndex_pooled",
                        message = T)
#> Warning: Stratum not provided. Treating all data as one stratum.
#> Cross-validation over etas sequence:
#>   |                                      |                              |   0%  |                                      |===                           |  10%  |                                      |======                        |  20%  |                                      |=========                     |  30%  |                                      |============                  |  40%  |                                      |===============               |  50%  |                                      |==================            |  60%  |                                      |=====================         |  70%  |                                      |========================      |  80%  |                                      |===========================   |  90%  |                                      |==============================| 100%
   
```
