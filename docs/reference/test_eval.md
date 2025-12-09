# Evaluate model performance on test data

Evaluates model performance on a test dataset using either the
log-partial-likelihood loss or the concordance index (C-index).

This function accepts either:

- `test_z` and `betahat`, which will be multiplied to obtain risk
  scores; or

- `test_RS`, a pre-computed numeric vector of risk scores.

## Usage

``` r
test_eval(
  test_z = NULL,
  test_RS = NULL,
  test_delta,
  test_time,
  test_stratum = NULL,
  betahat = NULL,
  criteria = c("loss", "CIndex")
)
```

## Arguments

- test_z:

  Optional numeric matrix or data frame of covariates for the test
  dataset. Required if `test_RS` is not provided.

- test_RS:

  Optional numeric vector of pre-computed risk scores (e.g., linear
  predictors). If provided, `test_z` and `betahat` are ignored.

- test_delta:

  Numeric vector of event indicators (1 for event, 0 for censoring).

- test_time:

  Numeric vector of survival times for the test dataset.

- test_stratum:

  Optional vector indicating stratum membership for each test
  observation. If `NULL`, all observations are assumed to belong to a
  single stratum.

- betahat:

  Optional numeric vector of estimated regression coefficients. Required
  if `test_RS` is not provided.

- criteria:

  Character string specifying the evaluation criterion. Must be either
  `"loss"` for log-partial-likelihood loss or `"CIndex"` for concordance
  index.

## Value

A numeric value representing either:

- the negative twice log-partial-likelihood (`criteria = "loss"`);

- or the concordance index (`criteria = "CIndex"`).

## Details

Observations are automatically sorted by stratum and time to ensure
correct risk set ordering before evaluation.
