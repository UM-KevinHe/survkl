# Cox Proportional Hazards Model with Ridge Penalty and External Information

Fits a Cox proportional hazards model using a ridge-type penalty (L2) on
all covariates. The model can integrate external information either as
precomputed risk scores (`RS`) or externally supplied coefficients
(`beta`). A tuning parameter `eta` controls the relative weight of the
external information. If `lambda` is not provided, a lambda sequence is
automatically generated.

## Usage

``` r
coxkl_ridge(
  z,
  delta,
  time,
  stratum = NULL,
  RS = NULL,
  beta = NULL,
  eta = NULL,
  lambda = NULL,
  nlambda = 100,
  lambda.min.ratio = ifelse(n_obs < n_vars, 0.01, 1e-04),
  penalty.factor = 0.999,
  tol = 1e-04,
  Mstop = 50,
  backtrack = FALSE,
  message = FALSE,
  data_sorted = FALSE,
  beta_initial = NULL,
  ...
)
```

## Arguments

- z:

  Numeric matrix of covariates (observations in rows, predictors in
  columns).

- delta:

  Numeric vector of event indicators (1 = event, 0 = censored).

- time:

  Numeric vector of observed times.

- stratum:

  Optional numeric or factor vector specifying strata.

- RS:

  Optional numeric vector or matrix of external risk scores.

- beta:

  Optional numeric vector of externally derived coefficients.

- eta:

  Non-negative scalar controlling the strength of external information.

- lambda:

  Optional numeric scalar or vector of penalty parameters. If `NULL`, a
  sequence is generated automatically.

- nlambda:

  Number of lambda values to generate if `lambda` is `NULL`.

- lambda.min.ratio:

  Ratio defining the minimum lambda relative to `lambda.max`.

- penalty.factor:

  Numeric scalar in `[0, 1)`. Controls the overall strength of the
  penalty when generating the ridge regression lambda sequence. Smaller
  values correspond to stronger penalization. Only used when
  `lambda = NULL`.

- tol:

  Convergence tolerance for the iterative estimation algorithm.

- Mstop:

  Maximum number of iterations for estimation.

- backtrack:

  Logical; if `TRUE`, uses backtracking line search.

- message:

  Logical; if `TRUE`, progress messages are printed during model
  fitting. Default is `FALSE`.

- data_sorted:

  Logical; if `TRUE`, assumes input data is already sorted by strata and
  time.

- beta_initial:

  Optional; default NULL. When NULL, the algorithm initializes
  beta_initial to a zero vector as a warm start

- ...:

  Additional arguments.

## Value

An object of class `"coxkl_ridge"` containing:

- `lambda`: The lambda sequence used for estimation.

- `beta`: Matrix of estimated coefficients for each lambda.

- `linear.predictors`: Matrix of linear predictors.

- `likelihood`: Vector of log-partial likelihoods.

- `data`: A list containing the input data used in fitting (`z`, `time`,
  `delta`, `stratum`, `data_sorted`).

## Details

The estimator maximizes a KL-regularized Cox partial log-likelihood with
a ridge (L2) penalty on all coefficients. External information is
incorporated via a KL term weighted by `eta`: if `beta` is supplied
(length `ncol(z)`), external risk scores are computed internally as
`RS = z %*% beta`; otherwise `RS` must be provided. If `lambda` is
`NULL`, a decreasing lambda path of length `nlambda` is generated using
`lambda.min.ratio` (its overall scale is influenced by
`penalty.factor`). Optimization proceeds along the lambda path with warm
starts (re-using the previous solution as `beta_initial`); when
`beta_initial = NULL`, the first step uses zeros. If
`data_sorted = FALSE`, data are sorted by `stratum` and `time` for
fitting and the returned linear predictors are mapped back to the
original observation order. `tol`, `Mstop`, and `backtrack` control
convergence and line search.

## Examples

``` r
data(ExampleData_highdim) 

train_dat_highdim <- ExampleData_highdim$train
beta_external_highdim <- ExampleData_highdim$beta_external

model_ridge <- coxkl_ridge(z = train_dat_highdim$z,
                           delta = train_dat_highdim$status,
                           time = train_dat_highdim$time,
                           stratum = NULL,
                           RS = NULL,
                           beta = beta_external_highdim,
                           message = TRUE)
#> Warning: eta is not provided. Setting eta = 0 (no external information used).
#> External beta information is used.
#> Warning: Stratum information not provided. All data is assumed to originate from a single stratum!
#> Cross-validation over lambda sequence:
#>   |                                      |                              |   0%  |                                      |                              |   1%  |                                      |=                             |   2%  |                                      |=                             |   3%  |                                      |=                             |   4%  |                                      |==                            |   5%  |                                      |==                            |   6%  |                                      |==                            |   7%  |                                      |==                            |   8%  |                                      |===                           |   9%  |                                      |===                           |  10%  |                                      |===                           |  11%  |                                      |====                          |  12%  |                                      |====                          |  13%  |                                      |====                          |  14%  |                                      |====                          |  15%  |                                      |=====                         |  16%  |                                      |=====                         |  17%  |                                      |=====                         |  18%  |                                      |======                        |  19%  |                                      |======                        |  20%  |                                      |======                        |  21%  |                                      |=======                       |  22%  |                                      |=======                       |  23%  |                                      |=======                       |  24%  |                                      |========                      |  25%  |                                      |========                      |  26%  |                                      |========                      |  27%  |                                      |========                      |  28%  |                                      |=========                     |  29%  |                                      |=========                     |  30%  |                                      |=========                     |  31%  |                                      |==========                    |  32%  |                                      |==========                    |  33%  |                                      |==========                    |  34%  |                                      |==========                    |  35%  |                                      |===========                   |  36%  |                                      |===========                   |  37%  |                                      |===========                   |  38%  |                                      |============                  |  39%  |                                      |============                  |  40%  |                                      |============                  |  41%  |                                      |=============                 |  42%  |                                      |=============                 |  43%  |                                      |=============                 |  44%  |                                      |==============                |  45%  |                                      |==============                |  46%  |                                      |==============                |  47%  |                                      |==============                |  48%  |                                      |===============               |  49%  |                                      |===============               |  50%  |                                      |===============               |  51%  |                                      |================              |  52%  |                                      |================              |  53%  |                                      |================              |  54%  |                                      |================              |  55%  |                                      |=================             |  56%  |                                      |=================             |  57%  |                                      |=================             |  58%  |                                      |==================            |  59%  |                                      |==================            |  60%  |                                      |==================            |  61%  |                                      |===================           |  62%  |                                      |===================           |  63%  |                                      |===================           |  64%  |                                      |====================          |  65%  |                                      |====================          |  66%  |                                      |====================          |  67%  |                                      |====================          |  68%  |                                      |=====================         |  69%  |                                      |=====================         |  70%  |                                      |=====================         |  71%  |                                      |======================        |  72%  |                                      |======================        |  73%  |                                      |======================        |  74%  |                                      |======================        |  75%  |                                      |=======================       |  76%  |                                      |=======================       |  77%  |                                      |=======================       |  78%  |                                      |========================      |  79%  |                                      |========================      |  80%  |                                      |========================      |  81%  |                                      |=========================     |  82%  |                                      |=========================     |  83%  |                                      |=========================     |  84%  |                                      |==========================    |  85%  |                                      |==========================    |  86%  |                                      |==========================    |  87%  |                                      |==========================    |  88%  |                                      |===========================   |  89%  |                                      |===========================   |  90%  |                                      |===========================   |  91%  |                                      |============================  |  92%  |                                      |============================  |  93%  |                                      |============================  |  94%  |                                      |============================  |  95%  |                                      |============================= |  96%  |                                      |============================= |  97%  |                                      |============================= |  98%  |                                      |==============================|  99%  |                                      |==============================| 100%
```
