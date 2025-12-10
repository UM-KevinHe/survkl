# Plot Cross-Validation Results vs Eta

Plots cross-validated performance across `eta` for `cv.coxkl`,
`cv.coxkl_ridge`, or `cv.coxkl_enet` results. The main CV curve is drawn
as a solid purple line; a green dotted horizontal reference line is
placed at the value corresponding to `eta = 0` (or the closest available
`eta`), with a solid green point marking that reference level.

## Usage

``` r
cv.plot(object, line_color = "#7570B3", baseline_color = "#1B9E77", ...)
```

## Arguments

- object:

  A fitted cross-validation result of class `"cv.coxkl"`,
  `"cv.coxkl_ridge"`, or `"cv.coxkl_enet"`.

- line_color:

  Color for the CV performance curve. Default `"#7570B3"`.

- baseline_color:

  Color for the horizontal reference line and point. Default
  `"#1B9E77"`.

- ...:

  Additional arguments (currently ignored).

## Value

A `ggplot` object showing cross-validation performance versus `eta`.

## Details

The function reads the performance metric from the object:

- For `"cv.coxkl"`: uses `object$internal_stat` (one row per `eta`).

- For `"cv.coxkl_ridge"` and `"cv.coxkl_enet"`: uses
  `object$integrated_stat.best_per_eta` (best `lambda` per `eta`).

The y-axis label is set to “Loss” if `criteria` in the object is “V&VH”
or “LinPred”; otherwise it is “C Index”. The horizontal reference
(“baseline”) is taken from the plotted series at `eta = 0` (or the
nearest `eta` present in the results).

## See also

[`cv.coxkl`](https://um-kevinhe.github.io/survkl/reference/cv.coxkl.md),
[`cv.coxkl_ridge`](https://um-kevinhe.github.io/survkl/reference/cv.coxkl_ridge.md),
[`cv.coxkl_enet`](https://um-kevinhe.github.io/survkl/reference/cv.coxkl_enet.md)

## Examples

``` r
data(Exampledata_lowdim)
#> Warning: data set 'Exampledata_lowdim' not found

train_dat_lowdim <- ExampleData_lowdim$train
beta_external_good_lowdim <- ExampleData_lowdim$beta_external_good

etas <- generate_eta(method = "exponential", n = 100, max_eta = 30)
cv_res <- cv.coxkl(z = train_dat_lowdim$z,
                   delta = train_dat_lowdim$status,
                   time = train_dat_lowdim$time,
                   stratrum = train_dat_lowdim$stratum,
                   beta = beta_external_good_lowdim,
                   etas = etas,
                   nfolds = 5,
                   criteria = c("V&VH"),
                   seed = 1)
#> Warning: Stratum not provided. Treating all data as one stratum.
cv.plot(cv_res)

```
