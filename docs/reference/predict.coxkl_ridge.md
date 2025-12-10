# Predict Linear Predictors from a coxkl_ridge Object

Computes linear predictors for new data using a fitted `coxkl_ridge`
model. If `lambda` is supplied, predictions are returned for those
`lambda` values; otherwise predictions are returned for all fitted
`lambda`s. When a requested `lambda` lies between fitted values,
coefficients are linearly interpolated.

## Usage

``` r
# S3 method for class 'coxkl_ridge'
predict(object, newz, lambda = NULL, ...)
```

## Arguments

- object:

  A fitted model object of class `"coxkl_ridge"`.

- newz:

  A numeric matrix or data frame of new covariates (same columns as in
  training data).

- lambda:

  Optional numeric value(s) specifying the regularization parameter(s)
  for which to predict. If `NULL`, predictions for all fitted `lambda`
  values are returned.

- ...:

  Additional arguments.

## Value

A numeric matrix of linear predictors. Each column corresponds to one
`lambda`, sorted in descending order.

## Details

The linear predictors are computed as `as.matrix(newz) %*% beta`.

## See also

[`coef.coxkl_ridge`](https://um-kevinhe.github.io/survkl/reference/coef.coxkl_ridge.md)
