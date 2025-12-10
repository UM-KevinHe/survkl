# Extract Coefficients from a `coxkl_ridge` Object

Extracts the estimated regression coefficients (`beta`) from a fitted
`coxkl_ridge` object. Optionally, one or more `lambda` values can be
supplied. If requested `lambda` values are not in the fitted sequence,
linear interpolation is performed between nearest neighbors;
out-of-range requests error.

## Usage

``` r
# S3 method for class 'coxkl_ridge'
coef(object, lambda = NULL, ...)
```

## Arguments

- object:

  An object of class `"coxkl_ridge"`, typically the result of
  [`coxkl_ridge`](https://um-kevinhe.github.io/survkl/reference/coxkl_ridge.md).

- lambda:

  Optional numeric value or vector specifying the regularization
  parameter(s) for which to extract (or interpolate) coefficients. If
  `NULL`, all estimated coefficients are returned.

- ...:

  Additional arguments (currently ignored).

## Value

A numeric matrix of regression coefficients. Each column corresponds to
one value of `lambda`, sorted in *descending* order.

## Examples

``` r
data(ExampleData_highdim) 

train_dat_highdim <- ExampleData_highdim$train
beta_external_highdim <- ExampleData_highdim$beta_external

model_ridge <- coxkl_ridge(z = train_dat_highdim$z,
                           delta = train_dat_highdim$status,
                           time = train_dat_highdim$time,
                           beta = beta_external_highdim,
                           eta = 1)
#> Warning: Stratum information not provided. All data is assumed to originate from a single stratum!
coef(model_ridge)[1:5, 1:10]
#>          82.0989     74.8055     68.1599     62.1048     56.5876     51.5605
#> [1,]  0.03940019  0.04317345  0.04718185  0.05142086  0.05588316  0.06055855
#> [2,] -0.08703965 -0.09270340 -0.09863286 -0.10483116 -0.11129985 -0.11803870
#> [3,]  0.04654658  0.05031996  0.05434817  0.05864173  0.06321089  0.06806550
#> [4,] -0.05680339 -0.06147303 -0.06646655 -0.07179722 -0.07747730 -0.08351777
#> [5,]  0.10912724  0.11518227  0.12137484  0.12769148  0.13411779  0.14063857
#>            46.98     42.8064     39.0036     35.5387
#> [1,]  0.06543395  0.07049345  0.07571841  0.08108773
#> [2,] -0.12504531 -0.13231486 -0.13983987 -0.14760999
#> [3,]  0.07321496  0.07866810  0.08443310  0.09051738
#> [4,] -0.08992804 -0.09671563 -0.10388587 -0.11144168
#> [5,]  0.14723786  0.15389915  0.16060546  0.16733947
```
