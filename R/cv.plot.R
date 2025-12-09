#' Plot Cross-Validation Results vs Eta
#'
#' @description
#' Plots cross-validated performance across \code{eta} for
#' \code{cv.coxkl}, \code{cv.coxkl_ridge}, or \code{cv.coxkl_enet} results.
#' The main CV curve is drawn as a solid purple line; a green dotted horizontal
#' reference line is placed at the value corresponding to \code{eta = 0}
#' (or the closest available \code{eta}), with a solid green point marking that
#' reference level.
#'
#' @param object A fitted cross-validation result of class \code{"cv.coxkl"},
#'   \code{"cv.coxkl_ridge"}, or \code{"cv.coxkl_enet"}.
#' @param line_color Color for the CV performance curve. Default \code{"#7570B3"}.
#' @param baseline_color Color for the horizontal reference line and point.
#'   Default \code{"#1B9E77"}.
#' @param ... Additional arguments (currently ignored).
#' 
#' @details
#' The function reads the performance metric from the object:
#' \itemize{
#'   \item For \code{"cv.coxkl"}: uses \code{object$internal_stat} (one row per \code{eta}).
#'   \item For \code{"cv.coxkl_ridge"} and \code{"cv.coxkl_enet"}:
#'         uses \code{object$integrated_stat.best_per_eta} (best \code{lambda} per \code{eta}).
#' }
#' The y-axis label is set to \dQuote{Loss} if \code{criteria} in the object is
#' \dQuote{V&VH} or \dQuote{LinPred}; otherwise it is \dQuote{C Index}.
#' The horizontal reference (“baseline”) is taken from the plotted series at
#' \code{eta = 0} (or the nearest \code{eta} present in the results).
#' 
#' @return A \code{ggplot} object showing cross-validation performance versus \code{eta}.
#' 
#' @seealso \code{\link{cv.coxkl}}, \code{\link{cv.coxkl_ridge}}, \code{\link{cv.coxkl_enet}}
#' 
#' @examples
#' data(Exampledata_lowdim)
#' 
#' train_dat_lowdim <- ExampleData_lowdim$train
#' beta_external_good_lowdim <- ExampleData_lowdim$beta_external_good
#' 
#' etas <- generate_eta(method = "exponential", n = 10, max_eta = 5)
#' etas <- sample(etas)
#' 
#' cv_res <- cv.coxkl(z = train_dat_lowdim$z,
#'                    delta = train_dat_lowdim$status,
#'                    time = train_dat_lowdim$time,
#'                    stratum = NULL,
#'                    RS = NULL,
#'                    beta = beta_external_good_lowdim,
#'                    etas = etas,
#'                    nfolds = 5,
#'                    criteria = c("LinPred"),   #"V&VH", "LinPred", "CIndex_pooled", "CIndex_foldaverage"
#'                    message = TRUE)
#' 
#' 
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_segment labs theme_minimal theme
#' @importFrom ggplot2 element_blank element_line element_text
#' @importFrom grid unit
#' @importFrom cowplot plot_grid get_legend
#' @export
cv.plot <- function(object,
                    line_color = "#7570B3",      # CoxKL main curve (purple)
                    baseline_color = "#1B9E77",  # Internal baseline & dot (green)
                    ...) {
  if (inherits(object, "cv.coxkl")) {
    df <- object$internal_stat
    criteria <- object$criteria
  } else if (inherits(object, "cv.coxkl_ridge") || inherits(object, "cv.coxkl_enet")) {
    df <- object$integrated_stat.best_per_eta
    criteria <- object$criteria
  } else {
    stop("Object must be a cv.coxkl, cv.coxkl_ridge, or cv.coxkl_enet.", call. = FALSE)
  }
  
  is_loss <- criteria %in% c("V&VH", "LinPred")
  ylab <- if (is_loss) "Loss" else "C Index"
  
  loss_candidates   <- c("VVH_Loss", "LinPred_Loss", "Loss", "loss")
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
  df$metric <- as.numeric(df[[metric_col]])

  if (!any(df$eta == 0)) idx0 <- which.min(abs(df$eta - 0)) else idx0 <- which(df$eta == 0)[1]
  baseline_val <- df$metric[idx0]
  baseline_eta <- df$eta[idx0]
  
  xmin <- min(df$eta, na.rm = TRUE)
  xmax <- max(df$eta, na.rm = TRUE)
  
  g_main <- ggplot2::ggplot(df, ggplot2::aes(x = eta, y = metric, group = 1)) +
    ggplot2::geom_line(linewidth = 1, color = line_color) +
    ggplot2::geom_point(size = 1.3, color = line_color) +
    ggplot2::geom_segment(  # green dotted baseline across full x-range
      data = data.frame(xmin = xmin, xmax = xmax, y = baseline_val),
      ggplot2::aes(x = xmin, xend = xmax, y = y, yend = y),
      inherit.aes = FALSE, color = baseline_color, linetype = "dotted", linewidth = 1
    ) +
    ggplot2::geom_point(    # green solid dot at eta≈0 baseline
      data = data.frame(eta = baseline_eta, metric = baseline_val),
      ggplot2::aes(x = eta, y = metric),
      inherit.aes = FALSE, color = baseline_color, shape = 16, size = 2.4
    ) +
    ggplot2::labs(x = expression(eta), y = ylab) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   axis.line = ggplot2::element_line(color = "black"),
                   axis.ticks.length = grid::unit(0.1, "cm"),
                   axis.ticks = ggplot2::element_line(color = "black"),
                   axis.text = ggplot2::element_text(size = 14),
                   legend.position = "none")
  
  legend_df <- data.frame(
    x = rep(c(0, 1), 2),
    y = rep(1, 4),
    Method = factor(rep(c("CoxKL", "Internal"), each = 2),
                    levels = c("CoxKL", "Internal"))
  )
  g_legend <- ggplot2::ggplot(legend_df, ggplot2::aes(x = x, y = y, color = Method, linetype = Method)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(
      data = subset(legend_df, Method == "Internal"),
      ggplot2::aes(x = 0.5, y = 1, color = Method),
      inherit.aes = FALSE, shape = 16, size = 2.4
    ) +
    ggplot2::scale_color_manual(values = c("CoxKL" = line_color, "Internal" = baseline_color)) +
    ggplot2::scale_linetype_manual(values = c("CoxKL" = "solid", "Internal" = "dotted")) +
    ggplot2::theme_void(base_size = 13) +
    ggplot2::theme(legend.position = "top",
                   legend.title = ggplot2::element_blank(),
                   legend.text = ggplot2::element_text(size = 14),
                   legend.key.width = grid::unit(1.2, "lines"),
                   legend.key.height = grid::unit(0.6, "lines")) +
    ggplot2::guides(color = ggplot2::guide_legend(keywidth = 1.2, keyheight = 0.4, title = NULL),
                    linetype = ggplot2::guide_legend(keywidth = 1.2, keyheight = 0.4, title = NULL))
  
  cowplot::plot_grid(cowplot::get_legend(g_legend), g_main, ncol = 1, rel_heights = c(0.08, 1))
}