# Generate a Sequence of Tuning Parameters (eta)

Produces a numeric vector of `eta` values to be used in Cox–KL model.

## Usage

``` r
generate_eta(method = "exponential", n = 10, max_eta = 5)
```

## Arguments

- method:

  Character string specifying the method to generate `eta`. Options are
  `"linear"` for a linearly spaced sequence, or `"exponential"` for an
  exponentially spaced sequence scaled to `max_eta`. Default is
  `"exponential"`.

- n:

  Integer, the number of `eta` values to generate. Default is 10.

- max_eta:

  Numeric, the maximum value of `eta` in the sequence. Default is 5.

## Value

Numeric vector of length `n` containing the generated `eta` values.

## Examples

``` r
# Generate 10 exponentially spaced eta values up to 5
generate_eta(method = "exponential", n = 10, max_eta = 5)
#>  [1] 0.00000000 0.03374245 0.09002825 0.18391863 0.34053721 0.60179276
#>  [7] 1.03759328 1.76455236 2.97719318 5.00000000

# Generate 5 linearly spaced eta values up to 3
generate_eta(method = "linear", n = 5, max_eta = 3)
#> Error in generate_eta(method = "linear", n = 5, max_eta = 3): object 'min_eta' not found
```
