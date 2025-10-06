#' Plot Cross-Validation Results vs Eta (Biometrics-style)
#'
#' @description
#' Plots cross-validation performance across eta values for
#' \code{cv.coxkl}, \code{cv.coxkl_ridge}, or \code{cv.coxkl_enet} objects
#' in the Biometrics figure style, with a solid blue line for CV performance
#' and a dotted orange line for the external baseline.
#'
#' @param object A fitted cross-validation result of class \code{"cv.coxkl"},
#'   \code{"cv.coxkl_ridge"}, or \code{"cv.coxkl_enet"}.
#' @param line_color Color for the CV performance curve. Default is \code{"#7570B3"}.
#' @param baseline_color Color for the external baseline line. Default is \code{"#D95F02"}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A \code{ggplot} object showing cross-validation performance versus \code{eta}.
#' @export
cv.plot <- function(object,
                    line_color = "#7570B3",
                    baseline_color = "#D95F02",
                    ...) {
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
  
  is_loss <- criteria %in% c("V&VH", "LinPred")
  ylab <- if (is_loss) "Loss" else "Concordance Index"
  loss_candidates <- c("VVH_Loss", "LinPred_Loss", "Loss", "loss")
  cindex_candidates <- c("CIndex_pooled", "CIndex_foldaverage", "CIndex", "cindex")
  candidates <- if (is_loss) loss_candidates else cindex_candidates
  metric_col <- NULL
  for (nm in candidates) if (nm %in% names(df)) { metric_col <- nm; break }
  if (is.null(metric_col)) {
    num_cols <- names(df)[vapply(df, is.numeric, logical(1))]
    num_cols <- setdiff(num_cols, c("eta", "lambda"))
    if (length(num_cols) == 0L) stop("Could not detect metric column in CV results.", call. = FALSE)
    metric_col <- num_cols[length(num_cols)]
  }
  
  df$eta <- as.numeric(df$eta)
  df <- df[order(df$eta), , drop = FALSE]
  metric_vals <- df[[metric_col]]
  
  legend_df <- data.frame(
    x = rep(c(0, 1), 2),
    y = rep(1, 4),
    Method = factor(rep(c("CoxKL", "External"), each = 2),
                    levels = c("CoxKL", "External"))
  )
  
  g_main <- ggplot(df, aes(x = eta, y = metric_vals)) +
    geom_line(aes(color = "CoxKL", linetype = "CoxKL"), size = 1) +
    geom_point(color = line_color, size = 1.3) +
    geom_segment(aes(x = min(df$eta), xend = max(df$eta),
                     y = external, yend = external,
                     color = "External", linetype = "External"),
                 linewidth = 1) +
    scale_color_manual(values = c("CoxKL" = line_color, "External" = baseline_color)) +
    scale_linetype_manual(values = c("CoxKL" = "solid", "External" = "dotted")) +
    labs(x = expression(eta), y = ylab) +
    theme_minimal(base_size = 13) +
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = "black"),
          axis.ticks.length = unit(0.1, "cm"),
          axis.ticks = element_line(color = "black"),
          axis.text = element_text(size = 14),
          legend.position = "none")
  
  g_legend <- ggplot(legend_df, aes(x = x, y = y, color = Method, linetype = Method)) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = c("CoxKL" = line_color, "External" = baseline_color)) +
    scale_linetype_manual(values = c("CoxKL" = "solid", "External" = "dotted")) +
    theme_void(base_size = 13) +
    theme(legend.position = "top",
          legend.title = element_blank(),
          legend.text = element_text(size = 14),
          legend.key.width = unit(1.2, "lines"),
          legend.key.height = unit(0.6, "lines")) +
    guides(color = guide_legend(keywidth = 1.2, keyheight = 0.4, title = NULL),
           linetype = guide_legend(keywidth = 1.2, keyheight = 0.4, title = NULL))
  
  legend <- cowplot::get_legend(g_legend)
  cowplot::plot_grid(legend, g_main, ncol = 1, rel_heights = c(0.08, 1))
}
