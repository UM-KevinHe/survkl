#' Evaluate model performance on test data
#'
#' @description
#' This function evaluates the performance of a fitted model on a test dataset using either the
#' log-partial-likelihood loss or the concordance index (C-index).
#'
#' @param test_z A numeric matrix or data frame of covariates for the test dataset.
#' @param test_delta A numeric vector of event indicators (1 for event, 0 for censoring).
#' @param test_time A numeric vector of survival times for the test dataset.
#' @param test_stratum An optional vector indicating stratum membership for each test observation.
#'   If `NULL`, all observations are assumed to belong to a single stratum.
#' @param betahat A numeric vector of estimated regression coefficients from the training model.
#' @param criteria A character string specifying the evaluation criterion.
#'   Must be either `"loss"` for log-partial-likelihood loss
#'   or `"CIndex"` for concordance index.
#'
#' @return
#' A numeric value representing either:
#' \itemize{
#'   \item the negative twice log-partial-likelihood (`criteria = "loss"`);
#'   \item or the concordance index (`criteria = "CIndex"`).
#' }
#'
#' @export
test_eval <- function(test_z, test_delta, test_time, test_stratum = NULL,
                      betahat, criteria = c("loss", "CIndex")){
  if (is.null(test_stratum)) {
    test_stratum <- rep(1, nrow(test_z))
  } else {
    test_stratum <- match(test_stratum, unique(test_stratum))
  }
  test_RS <- as.matrix(test_z) %*% as.matrix(betahat)
  n.each_test_stratum <- as.numeric(table(test_stratum))
  
  if (criteria == "loss"){
    test_loss <- -2 * pl_cal_theta(as.vector(test_RS), test_delta, n.each_test_stratum)
    return(test_loss)
  } else if (criteria == "CIndex"){
    
    test_c_index <- c_stat_stratcox(test_time, test_RS, test_stratum, test_delta)$c_statistic
    return(test_c_index)
  } else {
    stop("'criteria' must be either 'loss' or 'CIndex'!", call. = FALSE)
  }
}