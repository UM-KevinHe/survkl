# Predict Linear Predictors from a `coxkl` Object

Computes linear predictors for new data based on a fitted `coxkl` model.
If `eta` is supplied, predictions are returned for those `eta` values;
otherwise predictions are returned for all fitted `eta`s. Linear
interpolation is applied if an intermediate `eta` value is requested.

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

  Optional numeric vector of `eta` value(s) for which to predict. If
  `NULL`, predictions for all fitted `eta` values are returned.

- ...:

  Additional arguments.

## Value

A numeric matrix of linear predictors with one column per `eta` (sorted
ascending).

## Details

The linear predictors are computed as `as.matrix(newz) %*% beta`.

## See also

[`coef.coxkl`](https://um-kevinhe.github.io/survkl/reference/coef.coxkl.md)
