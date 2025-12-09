# Cox Proportional Hazards Model with KL Divergence for Data Integration and Lasso & Elastic Net Penalty

Fits a Cox proportional hazards model that incorporates external
information using Kullback–Leibler (KL) divergence, with an optional L1
(Lasso) or elastic net penalty on the coefficients. External information
can be supplied either as precomputed external risk scores (`RS`) or as
externally derived coefficients (`beta`). The integration strength is
controlled by the tuning parameter `eta`.

## Usage

``` r
coxkl_enet(
  z,
  delta,
  time,
  stratum = NULL,
  RS = NULL,
  beta = NULL,
  eta = NULL,
  alpha = NULL,
  lambda = NULL,
  nlambda = 100,
  lambda.min.ratio = ifelse(n < p, 0.05, 0.001),
  lambda.early.stop = FALSE,
  tol = 1e-04,
  Mstop = 1000,
  max.total.iter = (Mstop * nlambda),
  group = 1:ncol(z),
  group.multiplier = NULL,
  standardize = T,
  nvar.max = ncol(z),
  group.max = length(unique(group)),
  stop.loss.ratio = 0.001,
  actSet = TRUE,
  actIter = Mstop,
  actGroupNum = sum(unique(group) != 0),
  actSetRemove = F,
  returnX = FALSE,
  trace.lambda = FALSE,
  message = FALSE,
  data_sorted = FALSE,
  ...
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

- eta:

  Numeric tuning parameter controlling the reliance on external
  information. Larger values place more weight on the external source.

- alpha:

  Elastic-net mixing parameter in \\(0,1\]\\. When \\\alpha=1\\ the
  penalty is lasso.

- lambda:

  Optional nonnegative penalty parameter(s). If a numeric vector is
  supplied, the path is taken as-is. If `NULL`, a sequence is generated
  using `nlambda` and `lambda.min.ratio`.

- nlambda:

  Integer number of lambda values to generate when `lambda` is `NULL`.
  Default `100`.

- lambda.min.ratio:

  Ratio of the smallest to the largest lambda when generating a sequence
  (when `lambda` is `NULL`). Default `1e-3`.

- lambda.early.stop:

  Logical; if `TRUE`, stop traversing the lambda path early based on
  convergence or screening criteria. Default `FALSE`.

- tol:

  Convergence tolerance for the optimization algorithm. Default is
  `1e-3`.

- Mstop:

  Maximum number of iterations for the inner optimization at a given
  lambda. Default is `1000`.

- max.total.iter:

  Maximum total iterations across the entire lambda path. Default is
  `(Mstop * nlambda)`.

- group:

  Integer vector of group indices defining group membership of
  predictors for grouped penalties; use `0` to indicate unpenalized
  variables.

- group.multiplier:

  A vector of values representing multiplicative factors by which each
  covariate's penalty is to be multiplied. Default is a vector of 1's.

- standardize:

  Logical; if `TRUE`, columns of `z` are standardized prior to fitting,
  with coefficients re-scaled on output. Default `TRUE`.

- nvar.max:

  Integer cap on the number of active variables allowed during fitting.
  Default number of predictors.

- group.max:

  Integer cap on the number of active groups allowed during fitting.
  Default total number of groups.

- stop.loss.ratio:

  Relative improvement threshold for early stopping along the path;
  optimization may stop if objective gain falls below this value.
  Default `1e-3`.

- actSet:

  Logical; if `TRUE`, use an active-set strategy. Default `TRUE`.

- actIter:

  Maximum number of active-set refinement iterations per lambda. Default
  `Mstop`.

- actGroupNum:

  Maximum number of active groups allowed under the active-set scheme.

- actSetRemove:

  Logical; if `TRUE`, allow dropping variables/groups from the active
  set during iterations. Default `FALSE`.

- returnX:

  Logical; if `TRUE`, return standardized design and related internals
  in `result$returnX`. Default `FALSE`.

- trace.lambda:

  Logical; if `TRUE`, record path-wise traces across the lambda
  sequence. Default `FALSE`.

- message:

  Logical; if `TRUE`, progress messages are printed during model
  fitting. Default is `FALSE`.

- data_sorted:

  Logical; if `TRUE`, input is assumed already sorted by `stratum` then
  `time`. Default `FALSE`.

- ...:

  Additional arguments.

## Value

An object of class `"coxkl_enet"`, a list with components:

- `beta`:

  Coefficient estimates (vector or matrix across the path).

- `group`:

  A `factor` of the original group assignments.

- `lambda`:

  The lambda value(s) used or generated.

- `alpha`:

  The elastic-net mixing parameter used.

- `likelihood`:

  Vector of log-partial likelihoods for each lambda.

- `n`:

  Number of observations.

- `df`:

  Effective degrees of freedom (e.g., number of nonzero coefficients or
  group-adjusted count) along the path.

- `iter`:

  Number of iterations taken (per lambda and/or total).

- `W`:

  Exponentiated linear predictors on the original scale.

- `group.multiplier`:

  Group-specific penalty multipliers used.

- `returnX`:

  Only when `returnX = TRUE`: a list with elements `XX`
  (standardization/orthogonalization info from `std.Z`), `time`,
  `delta`, `stratum`, and `RS`.

## Details

Setting `lambda = 0` reduces to the unpenalized
[`coxkl`](https://umkevinhe.github.io/survkl/reference/coxkl.md) model.

When `lambda > 0`, the model fits a KL-regularized Cox objective with an
elastic-net penalty: \$\$\ell\_{\mathrm{KL}}(\beta;\eta) \\-\\
\lambda\Big\\ \alpha\\\beta\\\_1 \\+\\
(1-\alpha)\tfrac{1}{2}\\\beta\\\_2^2 \Big\\,\$\$ where \\\alpha=1\\
gives lasso and \\0\<\alpha\<1\\ gives elastic net. Grouped penalties
are supported via `group` (use `0` for unpenalized variables), with
optional per-group scaling through `group.multiplier`. If `lambda` is
`NULL`, a decreasing path of length `nlambda` is generated using
`lambda.min.ratio`; early stopping can prune the path
(`lambda.early.stop`, `stop.loss.ratio`). When `standardize = TRUE`,
predictors are standardized for fitting and coefficients are rescaled on
output. If `data_sorted = FALSE`, data are sorted by `stratum` then
`time` for optimization and predictions are returned in the original
order (reported via `W = exp(linear predictors)`). An active-set scheme
(`actSet`, `actIter`, `nvar.max`, `group.max`, `actGroupNum`,
`actSetRemove`) is used to accelerate the solution along the lambda
path.

## See also

[`coxkl`](https://umkevinhe.github.io/survkl/reference/coxkl.md)

## Examples

``` r
data(ExampleData_highdim) 

train_dat_highdim <- ExampleData_highdim$train
beta_external_highdim <- ExampleData_highdim$beta_external

model_enet <- coxkl_enet(z = train_dat_highdim$z,
                         delta = train_dat_highdim$status,
                         time = train_dat_highdim$time,
                         stratum = NULL,
                         RS = NULL,
                         beta = beta_external_highdim,
                         eta = 0,
                         alpha = 1.0,
                         message = TRUE)
#> External beta information is used.
#> Warning: Stratum information not provided. All data is assumed to originate from a single stratum!
```
