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
  [`coxkl`](https://umkevinhe.github.io/survkl/reference/coxkl.md).

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
data(Exampledata_lowdim)
#> Warning: data set 'Exampledata_lowdim' not found

train_dat_lowdim <- ExampleData_lowdim$train
beta_external_good_lowdim <- ExampleData_lowdim$beta_external_good

model <- coxkl(z = train_dat_lowdim$z,
               delta = train_dat_lowdim$status,
               time = train_dat_lowdim$time,
               stratum = train_dat_lowdim$stratum,
               RS = NULL,
               beta = beta_external_good_lowdim,
               etas = c(0:5))
coef(model)
#>                0           1          2          3          4          5
#> [1,]  0.25982263  0.25018262  0.2475909  0.2463857  0.2456895  0.2452363
#> [2,] -0.49599768 -0.39552949 -0.3674863 -0.3544124 -0.3468634 -0.3419523
#> [3,] -0.01972113  0.08827353  0.1221091  0.1386403  0.1484398  0.1549244
#> [4,] -0.84023506 -0.58441069 -0.5018446 -0.4608359 -0.4362843 -0.4199306
#> [5,]  0.24808718  0.29504012  0.3074159  0.3129777  0.3161117  0.3181155
#> [6,] -0.78551866 -0.55438847 -0.4798686 -0.4427800 -0.4205425 -0.4057148
```
