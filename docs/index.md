# survkl: An R Package for Transfer-Learning–Based Integrated Cox Models

**survkl** provides flexible and efficient tools for integrating
external risk scores into Cox proportional hazards models while
accounting for population heterogeneity. The package enables robust
estimation, improved predictive accuracy, and user-friendly workflows
for modern survival analysis.

## Key Features

- **Integration of External Risk Scores**  
  Seamlessly incorporate predictions from *any* external risk model
  (e.g., polygenic hazard scores, deep-learning–based risk scores).

- **Population Heterogeneity Adjustment**  
  Corrects for distributional differences between the external model’s
  training population and your study population.

- **Efficient Computation**  
  Designed for high-dimensional and large-scale survival datasets.

- **Improved Estimation and Prediction**  
  Demonstrated gains in estimation efficiency and predictive accuracy.

- **Built-In Cross-Validation**  
  Automated selection of tuning and penalization parameters.

## Installation

> Note: This package is under active development. Please report any
> issues you encounter.  
> Requires **R ≥ 4.0.0**.

Install from CRAN:

``` R
install.packages("survkl")
```

Or install the development version from GitHub:

``` R
require("devtools")
require("remotes")
remotes::install_github("UM-KevinHe/survkl")
```

## Detailed Tutorial

Full package documentation and parameter explanations:
[here](https://um-kevinhe.github.io/survkl/)

## Getting Help

If you encounter problems or bugs, please contact us:

- <ybshao@umich.edu>
- <lfluo@umich.edu>
- <xhliuu@umich.edu>
- <kevinhe@umich.edu>

------------------------------------------------------------------------

## References

1.  Wang, D., Ye, W., Zhu, J., Xu, G., Tang, W., Zawistowski, M.,
    Fritsche, L. G., & He, K. (2023). Incorporating external risk
    information with the Cox model under population heterogeneity:
    Applications to trans-ancestry polygenic hazard scores.
    *arXiv:2302.11123*.

2.  Luo, L., Taylor, J. M. G., Wang, D., & He, K. (2024). Flexible Deep
    Learning Techniques for Cox Models with Data Integration.
