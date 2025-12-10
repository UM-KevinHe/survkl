# Extract Coefficients from a `coxkl` Object

Extracts the estimated regression coefficients (`beta`) from a fitted
`coxkl` object. Optionally, a value (or vector) of `eta` can be
supplied. If the requested `eta` values are not in the fitted sequence,
linear interpolation is performed between the nearest neighboring `eta`
values; out-of-range requests error.

## Usage

``` r
# S3 method for class 'coxkl'
coef(object, eta = NULL, ...)
```

## Arguments

- object:

  An object of class `"coxkl"`, typically the result of
  [`coxkl`](https://um-kevinhe.github.io/survkl/reference/coxkl.md).

- eta:

  Optional numeric value or vector specifying the \\\eta\\ values for
  which to extract (or interpolate) coefficients. If `NULL`, all
  estimated coefficients are returned.

- ...:

  Additional arguments (currently ignored).

## Value

A numeric matrix of regression coefficients. Each column corresponds to
one value of `eta`, sorted in ascending order.

## Examples

``` r
data(ExampleData_lowdim)

train_dat_lowdim <- ExampleData_lowdim$train
beta_external_good_lowdim <- ExampleData_lowdim$beta_external_good
eta_list <- generate_eta(method = "exponential", n = 5, max_eta = 5)

model <- coxkl(z = train_dat_lowdim$z,
               delta = train_dat_lowdim$status,
               time = train_dat_lowdim$time,
               stratum = train_dat_lowdim$stratum,
               beta = beta_external_good_lowdim,
               etas = eta_list)
coef(model)
#>                0       0.1092      0.4545     1.5466          5
#> [1,]  0.25982263  0.257655537  0.25342872  0.2484849  0.2452363
#> [2,] -0.49599768 -0.474154973 -0.43022578 -0.3771785 -0.3419523
#> [3,] -0.01972113  0.002327028  0.04893181  0.1101809  0.1549244
#> [4,] -0.84023506 -0.788594162 -0.67862767 -0.5311457 -0.4199306
#> [5,]  0.24808718  0.258493491  0.27915167  0.3031994  0.3181155
#> [6,] -0.78551866 -0.738669998 -0.63932919 -0.5063335 -0.4057148
```
