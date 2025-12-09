# Package index

## Low-Dimensional Cox KL Integration

- [`coxkl()`](https://umkevinhe.github.io/survkl/reference/coxkl.md) :
  Cox Proportional Hazards Model with KL Divergence for Data Integration

- [`cv.coxkl()`](https://umkevinhe.github.io/survkl/reference/cv.coxkl.md)
  :

  Cross-Validated Selection of Integration Parameter （`eta`） for the
  Cox–KL Model

## High-Dimensional Cox KL Integration

- [`coxkl_ridge()`](https://umkevinhe.github.io/survkl/reference/coxkl_ridge.md)
  : Cox Proportional Hazards Model with Ridge Penalty and External
  Information
- [`coxkl_enet()`](https://umkevinhe.github.io/survkl/reference/coxkl_enet.md)
  : Cox Proportional Hazards Model with KL Divergence for Data
  Integration and Lasso & Elastic Net Penalty
- [`cv.coxkl_ridge()`](https://umkevinhe.github.io/survkl/reference/cv.coxkl_ridge.md)
  : Cross-Validation for CoxKL Ridge Model (eta tuning)
- [`cv.coxkl_enet()`](https://umkevinhe.github.io/survkl/reference/cv.coxkl_enet.md)
  : Cross-Validation for CoxKL Model with elastic net & lasso penalty

## Plotting Functions

- [`cv.plot()`](https://umkevinhe.github.io/survkl/reference/cv.plot.md)
  : Plot Cross-Validation Results vs Eta

- [`plot(`*`<coxkl>`*`)`](https://umkevinhe.github.io/survkl/reference/plot.coxkl.md)
  :

  Plot Model Performance vs Eta for `coxkl`

- [`plot(`*`<coxkl_ridge>`*`)`](https://umkevinhe.github.io/survkl/reference/plot.coxkl_ridge.md)
  :

  Plot Model Performance vs Lambda for `coxkl_ridge`

- [`plot(`*`<coxkl_enet>`*`)`](https://umkevinhe.github.io/survkl/reference/plot.coxkl_enet.md)
  :

  Plot Model Performance vs Lambda for `coxkl_enet`

## Coefficient & Prediction

- [`predict(`*`<coxkl>`*`)`](https://umkevinhe.github.io/survkl/reference/predict.coxkl.md)
  : Predict Linear Predictors from a coxkl Object

- [`predict(`*`<coxkl_ridge>`*`)`](https://umkevinhe.github.io/survkl/reference/predict.coxkl_ridge.md)
  : Predict Linear Predictors from a coxkl_ridge Object

- [`predict(`*`<coxkl_enet>`*`)`](https://umkevinhe.github.io/survkl/reference/predict.coxkl_enet.md)
  : Predict Linear Predictors from a coxkl_enet Object

- [`coef(`*`<coxkl>`*`)`](https://umkevinhe.github.io/survkl/reference/coef.coxkl.md)
  :

  Extract Coefficients from a `coxkl` Object

- [`coef(`*`<coxkl_ridge>`*`)`](https://umkevinhe.github.io/survkl/reference/coef.coxkl_ridge.md)
  :

  Extract Coefficients from a `coxkl_ridge` Object

- [`coef(`*`<coxkl_enet>`*`)`](https://umkevinhe.github.io/survkl/reference/coef.coxkl_enet.md)
  :

  Extract Coefficients from a `coxkl_enet` Object

## Utilities

- [`cal_surv_prob()`](https://umkevinhe.github.io/survkl/reference/cal_surv_prob.md)
  : Calculate Survival Probabilities
- [`loss_fn()`](https://umkevinhe.github.io/survkl/reference/loss_fn.md)
  : Calculate the Log-Partial Likelihood for a Stratified Cox Model
- [`generate_eta()`](https://umkevinhe.github.io/survkl/reference/generate_eta.md)
  : Generate a Sequence of Tuning Parameters (eta)
- [`test_eval()`](https://umkevinhe.github.io/survkl/reference/test_eval.md)
  : Evaluate model performance on test data

## Datasets

- [`ExampleData_lowdim`](https://umkevinhe.github.io/survkl/reference/ExampleData_lowdim.md)
  : Example low-dimensional survival data
- [`ExampleData_highdim`](https://umkevinhe.github.io/survkl/reference/ExampleData_highdim.md)
  : Example high-dimensional survival data
- [`support`](https://umkevinhe.github.io/survkl/reference/support.md) :
  Study to Understand Prognoses Preferences Outcomes and Risks of
  Treatment
