# Predict Linear Predictors from a coxkl Object

Computes linear predictors for new data based on a fitted `coxkl` model.
Users can specify one or more `eta` values; if not provided, predictions
are returned for all fitted `eta` values. Linear interpolation is
applied if an intermediate `eta` value is requested.

## Usage

``` r
# S3 method for class 'coxkl'
predict(object, newz, eta = NULL, ...)
```

## Arguments

- object:

  A fitted model object of class `"coxkl"`.

- newz:

  A numeric matrix or data frame of new covariates (must match the
  dimension of the training design matrix used to fit the model).

- eta:

  Optional numeric value(s) specifying which `eta` values to use for
  prediction. If `NULL`, predictions for all fitted `eta` values are
  returned.

- ...:

  Additional arguments (currently ignored).

## Value

A numeric matrix of linear predictors. Each column corresponds to one
value of `eta`, sorted in ascending order.
