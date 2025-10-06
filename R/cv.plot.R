#' Plot Cross-Validation Results vs Eta
#'
#' @description
#' Plots cross-validation performance across eta values for
#' \code{cv.coxkl}, \code{cv.coxkl_ridge}, or \code{cv.coxkl_enet} objects.
#'
#' @param object A fitted CV object of class \code{"cv.coxkl"},
#'   \code{"cv.coxkl_ridge"}, or \code{"cv.coxkl_enet"}.
#' @param add_baseline Logical; if \code{TRUE} (default), adds a dashed line for external baseline.
#' @param point_color Color for points (default: "black").
#' @param line_color Color for the CV curve (default: "#0072B2").
#' @param baseline_color Color for the baseline line (default: "gray40").
#' @param ... Additional plotting arguments (ignored).
#'
#' @return A \code{ggplot} object showing CV performance across eta.
#' @export
cv.plot <- function(object,
                    add_baseline = TRUE,
                    point_color = "black",
                    line_color = "#0072B2",
                    baseline_color = "gray40",
                    ...) {
  
  # ---- Identify object type & pull data ----
  if (inherits(object, "cv.coxkl")) {
    df <- object$internal_stat
    criteria <- object$criteria
    external <- object$external_stat
  } else if (inherits(object, "cv.coxkl_ridge") || inherits(object, "cv.coxkl_enet")) {
    df <- object$integrated_stat.best_per_eta
    criteria <- object$criteria
    external <- object$external_stat
  } else {
    stop("Object must be a cv.coxkl, cv.coxkl_ridge, or cv.coxkl_enet.", call. = FALSE)
  }
  
  # ---- Decide label & candidate metric columns by criteria ----
  is_loss <- criteria %in% c("V&VH", "LinPred")
  ylab <- if (is_loss) "Loss" else "Concordance Index"
  
  loss_candidates   <- c("VVH_Loss", "LinPred_Loss", "Loss", "loss")
  cindex_candidates <- c("CIndex_pooled", "CIndex_foldaverage", "CIndex", "cindex")
  
  candidates <- if (is_loss) loss_candidates else cindex_candidates
  # pick the first existing candidate
  metric_col <- NULL
  for (nm in candidates) {
    if (nm %in% names(df)) { metric_col <- nm; break }
  }
  # Fallback (defensive): last non-eta/non-lambda numeric column
  if (is.null(metric_col)) {
    num_cols <- names(df)[vapply(df, is.numeric, logical(1))]
    num_cols <- setdiff(num_cols, c("eta", "lambda"))
    if (length(num_cols) == 0L) stop("Could not detect metric column in CV results.", call. = FALSE)
    metric_col <- num_cols[length(num_cols)]
    warning(sprintf("Metric column not found by criteria; using '%s' instead.", metric_col), call. = FALSE)
  }
  
  # ---- Clean & order ----
  if (!("eta" %in% names(df))) stop("CV results must contain 'eta' column.", call. = FALSE)
  df$eta <- as.numeric(df$eta)
  df <- df[order(df$eta), , drop = FALSE]
  
  # ---- Plot ----
  g <- ggplot2::ggplot(df, ggplot2::aes(x = eta, y = .data[[metric_col]])) +
    ggplot2::geom_line(linewidth = 1, color = line_color) +
    ggplot2::geom_point(size = 2, color = point_color) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::labs(x = expression(eta), y = ylab)
  
  # ---- External baseline ----
  if (add_baseline && !is.null(external) && is.finite(external)) {
    g <- g +
      ggplot2::geom_hline(yintercept = external, linetype = "dashed",
                          color = baseline_color, linewidth = 0.9) +
      ggplot2::annotate("text",
                        x = max(df$eta, na.rm = TRUE),
                        y = external,
                        label = "External baseline",
                        hjust = 1.1, vjust = -0.5,
                        color = baseline_color, size = 3.5)
  }
  
  g
}
