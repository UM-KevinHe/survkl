# Predict Linear Predictors from a coxkl_enet Object

Computes linear predictors for new data from a fitted `coxkl_enet`
model. Users can specify one or more `lambda` values; if not provided,
predictions are returned for all fitted `lambda` values. Linear
interpolation is applied if an intermediate `lambda` value is requested.

## Usage

``` r
# S3 method for class 'coxkl_enet'
predict(object, newz, lambda = NULL, ...)
```

## Arguments

- object:

  A fitted model object of class `"coxkl_enet"`.

- newz:

  A numeric matrix or data frame of new covariates (same columns as in
  training data).

- lambda:

  Optional numeric value(s) specifying the regularization parameter(s)
  for which to predict. If `NULL`, predictions for all fitted `lambda`
  values are returned.

- ...:

  Additional arguments (currently ignored).

## Value

A numeric matrix of linear predictors. Each column corresponds to one
`lambda`, sorted in descending order.
