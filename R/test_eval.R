#' Evaluate model performance on test data
#'
#' @description
#' Evaluates model performance on a test dataset using either the
#' log-partial-likelihood loss or the concordance index (C-index).
#' 
#' This function accepts either:
#' \itemize{
#'   \item \code{test_z} and \code{betahat}, which will be multiplied to obtain risk scores; or
#'   \item \code{test_RS}, a pre-computed numeric vector of risk scores.
#' }
#'
#' @param test_z Optional numeric matrix or data frame of covariates for the test dataset.
#'   Required if \code{test_RS} is not provided.
#' @param test_RS Optional numeric vector of pre-computed risk scores (e.g., linear predictors).
#'   If provided, \code{test_z} and \code{betahat} are ignored.
#' @param test_delta Numeric vector of event indicators (1 for event, 0 for censoring).
#' @param test_time Numeric vector of survival times for the test dataset.
#' @param test_stratum Optional vector indicating stratum membership for each test observation.
#'   If \code{NULL}, all observations are assumed to belong to a single stratum.
#' @param betahat Optional numeric vector of estimated regression coefficients.
#'   Required if \code{test_RS} is not provided.
#' @param criteria Character string specifying the evaluation criterion.
#'   Must be either \code{"loss"} for log-partial-likelihood loss
#'   or \code{"CIndex"} for concordance index.
#'
#' @return
#' A numeric value representing either:
#' \itemize{
#'   \item the negative twice log-partial-likelihood (\code{criteria = "loss"});
#'   \item or the concordance index (\code{criteria = "CIndex"}).
#' }
#'
#' @details
#' Observations are automatically sorted by stratum and time to ensure correct
#' risk set ordering before evaluation.
#'
#' @export
test_eval <- function(test_z = NULL, test_RS = NULL, test_delta, test_time,
                      test_stratum = NULL, betahat = NULL,
                      criteria = c("loss", "CIndex")) {
  
  criteria <- match.arg(criteria)
  
  # ---- Validation ----
  if (is.null(test_RS)) {
    if (is.null(test_z) || is.null(betahat))
      stop("Either 'test_RS' must be provided, or both 'test_z' and 'betahat' must be specified.", call. = FALSE)
    test_RS <- as.vector(as.matrix(test_z) %*% as.matrix(betahat))
  }
  
  test_time <- as.numeric(test_time)
  test_delta <- as.numeric(test_delta)
  n <- length(test_time)
  if (is.null(test_stratum)) {
    test_stratum <- rep(1, n)
  } else {
    test_stratum <- match(test_stratum, unique(test_stratum))
  }
  
  # ---- Ensure proper ordering ----
  order_idx <- order(test_stratum, test_time)
  test_RS <- test_RS[order_idx]
  test_delta <- test_delta[order_idx]
  test_time <- test_time[order_idx]
  test_stratum <- test_stratum[order_idx]
  
  n.each_test_stratum <- as.numeric(table(test_stratum))
  
  # ---- Evaluation ----
  if (criteria == "loss") {
    test_loss <- -2 * pl_cal_theta(test_RS, test_delta, n.each_test_stratum)
    return(test_loss)
  } else if (criteria == "CIndex") {
    test_c_index <- c_stat_stratcox(test_time, test_RS, test_stratum, test_delta)$c_statistic
    return(test_c_index)
  } else {
    stop("'criteria' must be either 'loss' or 'CIndex'!", call. = FALSE)
  }
}
