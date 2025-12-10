# Extract Coefficients from a `coxkl_enet` Object

Extracts the estimated regression coefficients (`beta`) from a fitted
`coxkl_enet` object. Optionally, one or more `lambda` values can be
supplied. If requested `lambda` values are not in the fitted sequence,
linear interpolation is performed between nearest neighbors;
out-of-range requests error.

## Usage

``` r
# S3 method for class 'coxkl_enet'
coef(object, lambda = NULL, ...)
```

## Arguments

- object:

  An object of class `"coxkl_enet"`, typically the result of
  [`coxkl_enet`](https://um-kevinhe.github.io/survkl/reference/coxkl_enet.md).

- lambda:

  Optional numeric value or vector specifying the regularization
  parameter(s) for which to extract (or interpolate) coefficients. If
  `NULL`, all estimated coefficients are returned.

- ...:

  Additional arguments (currently ignored).

## Value

A numeric matrix of regression coefficients; each column corresponds to
one value of `lambda`, sorted in *descending* order.

## Examples

``` r
data(ExampleData_highdim)

train_dat_highdim <- ExampleData_highdim$train
beta_external_highdim <- ExampleData_highdim$beta_external


enet_model <- coxkl_enet(z = train_dat_highdim$z,
                         delta = train_dat_highdim$status,
                         time = train_dat_highdim$time,
                         beta = beta_external_highdim,
                         eta = 1,
                         alpha = 1.0)
coef(enet_model)[1:5, 1:10]
#>    0.0853     0.0796     0.0742     0.0692     0.0645     0.0602       0.0561
#> Z1      0 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000  0.000000000
#> Z2      0 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 -0.009001751
#> Z3      0 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000  0.000000000
#> Z4      0 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 -0.009194676
#> Z5      0 0.01945713 0.03849402 0.05620049 0.07282693 0.08844987  0.101159122
#>         0.0523      0.0488       0.0455
#> Z1  0.00000000  0.00000000  0.000000000
#> Z2 -0.02004127 -0.02981263 -0.038485677
#> Z3  0.00000000  0.00000000  0.005058644
#> Z4 -0.03565942 -0.06036793 -0.083632736
#> Z5  0.11223900  0.12286450  0.132785714
   
```
