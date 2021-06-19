#' Generate a QC Histogram Plot for a single metric
#'
#' @param meta A data.frame containing metadata
#' @param column A character object specifying the metadata to display
#' @param name_x A character object specifying a name to display on the x-axis
#' @param log_x A logical indicating whether or not to log10-scale the x-axis. Default is TRUE.
#' @param fill A character object specifying the color to use for for the histogram. Default is "dodgerblue".
#' @param target A numeric value for a target line to display on the x-axis. Default is 2e4.
#' @param y_max A numeric value for the maximum value on the y-axis. Default is 2e3.
#'
#' @return a ggplot2 plot object
#' @export
qc_hist_plot <- function(meta,
                         column = "n_reads",
                         name_x = "N Reads per Cell",
                         log_x = TRUE,
                         fill = "dodgerblue",
                         target = 2e4,
                         y_max = 2e3) {
  
  
  assertthat::assert_that(sum(class(meta) %in% c("data.frame","data.table")) > 0)
  assertthat::assert_that(class(column) == "character")
  assertthat::assert_that(length(column) == 1)
  assertthat::assert_that(column %in% names(meta))
  assertthat::assert_that(class(name_x) == "character")
  assertthat::assert_that(length(name_x) == 1)
  assertthat::assert_that(class(log_x) == "logical")
  assertthat::assert_that(length(log_x) == 1)
  assertthat::assert_that(class(fill) == "character")
  assertthat::assert_that(length(fill) == 1)
  assertthat::assert_that(class(target) == "numeric")
  assertthat::assert_that(length(target) == 1)
  assertthat::assert_that(class(y_max) == "numeric")
  assertthat::assert_that(length(y_max) == 1)
  
  if(log_x) {
    binwidth <- 0.02
  } else {
    binwidth <- 1e3
  }
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_histogram(ggplot2::aes(x = meta[[column]]),
                            binwidth = binwidth,
                            fill = fill) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = stats::median(meta[[column]])),
                        linetype = "dashed",
                        color = "#000000") +
    ggplot2::geom_text(ggplot2::aes(x = stats::median(meta[[column]]) * .95,
                                    y = 1450,
                                    label = paste0("median: ", stats::median(meta[[column]]))),
                       color = "#000000",
                       hjust = 1,
                       vjust = 1) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = target),
                        linetype = "dashed",
                        color = "#808080") +
    ggplot2::geom_text(ggplot2::aes(x = target * 1.05,
                                    y = 1450,
                                    label = paste0(target / 1e3, "k")),
                       color = "#808080",
                       hjust = 0,
                       vjust = 1) +
    ggplot2::scale_color_identity() +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_y_continuous("N Cells", limits = c(0, y_max)) +
    ggplot2::theme_bw()
  
  if(log_x) {
    p <- p +
      ggplot2::scale_x_log10(paste0("log10(",name_x,")"),
                             limits = c(1e2, 2.5e5),
                             breaks = c(1e2, 5e2, 1e3, 5e3, 1e4, 5e4, 1e5, 2.5e5))
  } else {
    p <- p +
      ggplot2::scale_x_continuous(name_x,
                                  limits = c(1e2, 2.5e5),
                                  breaks = c(1e2, 5e2, 1e3, 5e3, 1e4, 5e4, 1e5, 2.5e5))
  }
  
  p
}

#' Generate a QC Histogram Plot for a single fractional metric
#'
#' @param meta A data.frame containing metadata
#' @param column A character object specifying the metadata to display
#' @param name_x A character object specifying a name to display on the x-axis
#' @param fill A character object specifying the color to use for for the histogram. Default is "dodgerblue".
#' @param target A numeric value for a target line to display on the x-axis. Default is 0.5,
#' @param y_max A numeric value for the maximum value on the y-axis. Default is 2e3.
#'
#' @return a ggplot2 plot object
#' @export
qc_frac_hist_plot <- function(meta,
                              column = "n_reads",
                              name_x = "N Reads per Cell",
                              fill = "dodgerblue",
                              target = 0.5,
                              y_max = 2e3) {
  
  
  assertthat::assert_that(sum(class(meta) %in% c("data.frame","data.table")) > 0)
  assertthat::assert_that(class(column) == "character")
  assertthat::assert_that(length(column) == 1)
  assertthat::assert_that(column %in% names(meta))
  assertthat::assert_that(class(name_x) == "character")
  assertthat::assert_that(length(name_x) == 1)
  assertthat::assert_that(class(fill) == "character")
  assertthat::assert_that(length(fill) == 1)
  assertthat::assert_that(class(target) == "numeric")
  assertthat::assert_that(length(target) == 1)
  assertthat::assert_that(class(y_max) == "numeric")
  assertthat::assert_that(length(y_max) == 1)
  
  binwidth <- 0.02
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_histogram(ggplot2::aes(x = meta[[column]]),
                            binwidth = binwidth,
                            fill = fill) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = stats::median(meta[[column]])),
                        linetype = "dashed",
                        color = "#000000") +
    ggplot2::geom_text(ggplot2::aes(x = stats::median(meta[[column]]) * .95,
                                    y = 1450,
                                    label = paste0("median: ", round(stats::median(meta[[column]]), 3))),
                       color = "#000000",
                       hjust = 1,
                       vjust = 1) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = target),
                        linetype = "dashed",
                        color = "#808080") +
    ggplot2::geom_text(ggplot2::aes(x = target * 1.05,
                                    y = 1450,
                                    label = target),
                       color = "#808080",
                       hjust = 0,
                       vjust = 1) +
    ggplot2::scale_color_identity() +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_y_continuous("N Cells", limits = c(0, y_max)) +
    ggplot2::theme_bw()
  
  
  p <- p +
    ggplot2::scale_x_continuous(name_x,
                                limits = c(-.02, 1.02),
                                breaks = seq(0, 1, by = 0.1))
  
  p
}

#' Generate a QC Scatter Plot for a pair of metrics
#'
#' @param meta A data.frame containing metadata
#' @param column_x A character object specifying the metadata to display on the x-axis
#' @param name_x A character object specifying a name to display on the x-axis
#' @param column_y A character object specifying the metadata to display on the y-axis
#' @param name_y A character object specifying a name to display on the y-axis
#' @param log_x A logical indicating whether or not to log10-scale the x-axis. Default is TRUE.
#' @param frac_x A logical indicating whether or not to scale the x-axis between 0 and 1. Default is FALSE
#' @param log_y A logical indicating whether or not to log10-scale the y-axis. Default is TRUE.
#' @param frac_y A logical indicating whether or not to scale the y-axis between 0 and 1. Default is FALSE
#' @param show_targets A logical indicating whether or not to plot lines displaying ratios of values. Default is TRUE.
#' @param color A character object specifying the color to use for for the points. Default is "dodgerblue".
#'
#' @return a ggplot2 plot object
#' @export
qc_scatter_plot <- function(meta,
                            column_x = "n_reads",
                            name_x = "N Reads per Cell",
                            column_y = "n_umis",
                            name_y = "N UMIs per Cell",
                            log_x = TRUE,
                            frac_x = FALSE,
                            log_y = TRUE,
                            frac_y = FALSE,
                            show_targets = TRUE,
                            color = "dodgerblue") {
  
  assertthat::assert_that(sum(class(meta) %in% c("data.frame","data.table")) > 0)
  assertthat::assert_that(class(column_x) == "character")
  assertthat::assert_that(length(column_x) == 1)
  assertthat::assert_that(column_x %in% names(meta))
  assertthat::assert_that(class(name_x) == "character")
  assertthat::assert_that(length(name_x) == 1)
  assertthat::assert_that(class(column_y) == "character")
  assertthat::assert_that(length(column_y) == 1)
  assertthat::assert_that(column_y %in% names(meta))
  assertthat::assert_that(class(name_y) == "character")
  assertthat::assert_that(length(name_y) == 1)
  assertthat::assert_that(class(log_x) == "logical")
  assertthat::assert_that(length(log_x) == 1)
  assertthat::assert_that(class(log_y) == "logical")
  assertthat::assert_that(length(log_y) == 1)
  assertthat::assert_that(class(show_targets) == "logical")
  assertthat::assert_that(length(show_targets) == 1)
  assertthat::assert_that(class(color) == "character")
  assertthat::assert_that(length(color) == 1)
  
  target_lines <- data.frame(x = c(1e2, 2e2, 4e2, 8e2),
                             xend = c(2.5e5, 2.5e5, 2.5e5, 2.5e5),
                             y = c(1e2, 1e2, 1e2, 1e2),
                             yend = c(2.5e5, 1.25e5, 6.25e4, 3.125e4),
                             group = c("1:1","1:2","1:4", "1:8"))
  
  p <- ggplot2::ggplot()
  
  if(show_targets) {
    p <- p +
      ggplot2::geom_segment(data = target_lines,
                            ggplot2::aes(x = x, xend = xend,
                                         y = y, yend = yend,
                                         group = group),
                            linetype = "dashed") +
      ggplot2::geom_text(data = target_lines,
                         ggplot2::aes(x = xend * 0.9,
                                      y = yend,
                                      label = group),
                         angle = 45,
                         hjust = 1,
                         vjust = 0,
                         size = 3)
  }
  
  p <- p +
    ggplot2::geom_point(ggplot2::aes(x = meta[[column_x]],
                                     y = meta[[column_y]]),
                        alpha = 0.2,
                        size = 0.2,
                        color = color) +
    ggplot2::scale_color_identity() +
    ggplot2::theme_bw()
  
  if(log_x) {
    p <- p +
      ggplot2::scale_x_log10(paste0("log10(",name_x,")"),
                             limits = c(1e2, 2.5e5),
                             breaks = c(1e2, 5e2, 1e3, 5e3, 1e4, 5e4, 1e5, 2.5e5),
                             labels = c("100", "500", "1k", "5k", "10k", "50k", "100k", "250k"))
  } else if(frac_x) {
    p <- p +
      ggplot2::scale_x_continuous(name_x,
                                  limits = c(0, 1),
                                  breaks = seq(0, 1, by = 0.1))
  } else {
    p <- p +
      ggplot2::scale_x_continuous(name_x)
  }
  
  if(log_y) {
    p <- p +
      ggplot2::scale_y_log10(paste0("log10(",name_y,")"),
                             limits = c(1e2, 2.5e5),
                             breaks = c(1e2, 5e2, 1e3, 5e3, 1e4, 5e4, 1e5, 2.5e5),
                             labels = c("100", "500", "1k", "5k", "10k", "50k", "100k", "250k"))
  } else if(frac_y) {
    p <- p +
      ggplot2::scale_y_continuous(name_y,
                                  limits = c(0, 1),
                                  breaks = seq(0, 1, by = 0.1))
  } else {
    p <- p +
      ggplot2::scale_y_continuous(name_y)
  }
  
  p
}

#' Generate a QC Violin Plot for a metric, grouped by a categorical metadata column.
#'
#' @param meta A data.frame containing metadata
#' @param category_x A character object specifying the metadata to use for grouping on the x-axis
#' @param name_x A character object specifying a name to display on the x-axis
#' @param column_y A character object specifying the metadata to display on the y-axis
#' @param name_y A character object specifying a name to display on the y-axis
#' @param log_y A logical indicating whether or not to log10-scale the y-axis. Default is TRUE.
#' @param fill A character object specifying the fill color to use for for the violins. Default is "skyblue".
#'
#' @return a ggplot2 plot object
#'
#' @importFrom rlang parse_expr !!
#'
#' @export
qc_violin_plot <- function(meta,
                           category_x = "well_id",
                           name_x = "Well ID",
                           column_y = "n_reads",
                           name_y = "N Reads per Cell",
                           log_y = TRUE,
                           fill = "skyblue") {
  
  assertthat::assert_that(sum(class(meta) %in% c("data.frame","data.table")) > 0)
  assertthat::assert_that(class(category_x) == "character")
  assertthat::assert_that(length(category_x) == 1)
  assertthat::assert_that(category_x %in% names(meta))
  assertthat::assert_that(class(name_x) == "character")
  assertthat::assert_that(length(name_x) == 1)
  assertthat::assert_that(class(column_y) == "character")
  assertthat::assert_that(length(column_y) == 1)
  assertthat::assert_that(column_y %in% names(meta))
  assertthat::assert_that(class(name_y) == "character")
  assertthat::assert_that(length(name_y) == 1)
  assertthat::assert_that(class(log_y) == "logical")
  assertthat::assert_that(length(log_y) == 1)
  assertthat::assert_that(class(fill) == "character")
  assertthat::assert_that(length(fill) == 1)
  
  tidy_x <- rlang::parse_expr(category_x)
  tidy_y <- rlang::parse_expr(column_y)
  
  meta <- as.data.table(meta)
  q_table <- meta[, .(q_25 = stats::quantile(get(column_y), 0.25),
                      q_50 = stats::quantile(get(column_y), 0.50),
                      q_75 = stats::quantile(get(column_y), 0.75)),
                  by = get(category_x)]
  names(q_table)[1] <- category_x
  
  global_median <- stats::median(meta[[column_y]])
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_violin(ggplot2::aes(x = as.factor(meta[[category_x]]),
                                      y = meta[[column_y]]),
                         fill = fill,
                         color = "#808080") +
    ggplot2::geom_errorbar(data = q_table,
                           ggplot2::aes(x = as.factor(!!tidy_x),
                                        ymin = q_25,
                                        ymax = q_75),
                           width = 0.25) +
    ggplot2::geom_point(data = q_table,
                        ggplot2::aes(x = as.factor(!!tidy_x),
                                     y = q_50)) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = global_median),
                        linetype = "dashed") +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_color_identity() +
    ggplot2::scale_x_discrete(name_x) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = element_text(angle = 90,
                                              hjust = 1,
                                              vjust = 0.3))
  
  
  if(log_y) {
    p <- p +
      ggplot2::scale_y_log10(paste0("log10(",name_y,")"),
                             limits = c(1e2, 2.5e5),
                             breaks = c(1e2, 5e2, 1e3, 5e3, 1e4, 5e4, 1e5, 2.5e5),
                             labels = c("100", "500", "1k", "5k", "10k", "50k", "100k", "250k"))
  } else {
    p <- p +
      ggplot2::scale_y_continuous(name_y)
  }
  
  p
  
}

#' Generate a QC Barplot for a metric at multiple cutoffs
#'
#' @param meta A data.frame containing metadata
#' @param column_x A character object specifying the metadata to use plotting
#' @param name_x A character object specifying a name to display on the x-axis
#' @param cutoffs A numeric vector specifying one or more cutoffs to use. Default is c(500, 750, 1000).
#' @param max_y A numeric value specifying the maximum value to use for the y-axis. Default is 3e4.
#' @param fill A character object specifying the fill color to use for for the violins. Default is "purple".
#'
#' @return a ggplot2 plot object
#' @export
qc_cutoff_barplot <- function(meta,
                              column_x = "n_umis",
                              name_x = "N UMIs",
                              cutoffs = c(500, 750, 1000),
                              max_y = 3e4,
                              fill = "purple") {
  
  assertthat::assert_that(sum(class(meta) %in% c("data.frame","data.table")) > 0)
  assertthat::assert_that(class(column_x) == "character")
  assertthat::assert_that(length(column_x) == 1)
  assertthat::assert_that(column_x %in% names(meta))
  assertthat::assert_that(class(name_x) == "character")
  assertthat::assert_that(length(name_x) == 1)
  assertthat::assert_that(class(cutoffs) == "numeric")
  assertthat::assert_that(length(cutoffs) > 0)
  assertthat::assert_that(class(fill) == "character")
  assertthat::assert_that(length(fill) == 1)
  
  cutoff_counts <- data.frame(cutoff = cutoffs,
                              n_cells = sapply(cutoffs,
                                               function(x) {
                                                 sum(meta[[column_x]] > x)
                                               }))
  
  cutoff_counts <- cutoff_counts[order(cutoff_counts$cutoff),]
  
  ggplot2::ggplot() +
    ggplot2::geom_bar(data = cutoff_counts,
                      ggplot2::aes(x = as.factor(cutoff),
                                   y = n_cells),
                      stat = "identity",
                      fill = fill) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_y_continuous("N Cells",
                                limits = c(0, max_y),
                                expand = c(0, 0)) +
    ggplot2::scale_x_discrete(paste0(name_x, " Cutoff")) +
    ggplot2::theme_bw()
  
}

#' Generate a stacked barplot for two categorical metrics
#'
#' The metric used for category_x will generate bars as columns on the x-axis.
#' The bar will be split vertically based on category_y.
#'
#' @param meta A data.frame containing metadata
#' @param category_x A character object specifying the metadata to use for grouping on the x-axis
#' @param name_x A character object specifying a name to display on the x-axis
#' @param category_x A character object specifying the metadata to use for splitting in the y-direction
#' @param category_name A character object specifying a name to display for the colors
#' @param colorset_y A colorset to use as fills for category_y. Currently supported: "rainbow" or "varibow". Default is "varibow"
#' @param name_y A character object specifying a name for the y-axis.
#' @param as_fraction A logical object specifying whether or not to display the stacked bars as fractions of the total count for each category_x. Default is FALSE.
#'
#' @return a ggplot2 plot object
#' @export
qc_stacked_barplot <- function(meta,
                               category_x = "batch_id",
                               name_x = "Batch ID",
                               category_y = "well_id",
                               category_name = "Well ID",
                               colorset_y = "varibow",
                               name_y = "N Cells",
                               as_fraction = FALSE) {
  
  assertthat::assert_that(sum(class(meta) %in% c("data.frame","data.table")) > 0)
  assertthat::assert_that(class(category_x) == "character")
  assertthat::assert_that(length(category_x) == 1)
  assertthat::assert_that(category_x %in% names(meta))
  assertthat::assert_that(class(name_x) == "character")
  assertthat::assert_that(length(name_x) == 1)
  assertthat::assert_that(class(category_y) == "character")
  assertthat::assert_that(length(category_y) == 1)
  assertthat::assert_that(category_y %in% names(meta))
  assertthat::assert_that(class(category_name) == "character")
  assertthat::assert_that(length(category_name) == 1)
  assertthat::assert_that(class(name_y) == "character")
  assertthat::assert_that(length(name_y) == 1)
  assertthat::assert_that(class(colorset_y) == "character")
  assertthat::assert_that(length(colorset_y) == 1)
  assertthat::assert_that(colorset_y %in% c("rainbow","varibow"))
  assertthat::assert_that(class(as_fraction) == "logical")
  assertthat::assert_that(length(as_fraction) == 1)
  
  
  meta <- as.data.table(meta)
  count_table <- meta[, .(n_cells = nrow(.SD)),
                      by = mget(c(category_x, category_y))]
  
  plot_xpos <- data.frame(unique(count_table[[category_x]]))
  names(plot_xpos) <- category_x
  plot_xpos <- plot_xpos[order(plot_xpos[[category_x]]),,drop = FALSE]
  plot_xpos$xpos <- 1:nrow(plot_xpos)
  
  count_table <- count_table[plot_xpos, on = category_x]
  
  plot_fills <- data.frame(unique(count_table[[category_y]]))
  names(plot_fills) <- category_y
  if(colorset_y == "rainbow") {
    set.seed(3030)
    plot_fills$fill <- sample(grDevices::rainbow(nrow(plot_fills)), nrow(plot_fills))
  } else if(colorset_y == "varibow") {
    set.seed(3030)
    plot_fills$fill <- sample(BarMixer::varibow(nrow(plot_fills)), nrow(plot_fills))
  }
  plot_fills <- plot_fills[order(plot_fills[[category_y]]),]
  
  count_table <- count_table[plot_fills, on = category_y]
  count_table <- count_table[order(get(category_y), decreasing = TRUE)]
  
  if(as_fraction) {
    count_table <- count_table[, ymax := cumsum(n_cells)/sum(n_cells), by = list(get(category_x))]
    count_table <- count_table[, ymin := shift(ymax, fill = 0, type = "lag"), by = list(get(category_x))]
  }
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(data = count_table,
                       ggplot2::aes(xmin = xpos - 0.4,
                                    xmax = xpos + 0.4,
                                    ymin = ymin,
                                    ymax = ymax,
                                    fill = fill)) +
    ggplot2::scale_fill_identity(category_name,
                                 breaks = plot_fills$fill,
                                 labels = plot_fills[[category_y]],
                                 guide = "legend") +
    ggplot2::scale_x_continuous(name_x,
                                breaks = plot_xpos$xpos,
                                labels = plot_xpos[[category_x]]) +
    ggplot2::scale_y_continuous(name_y) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                       hjust = 1,
                                                       vjust = 0.3))
  
  p
}


#' Generate a baseline-aligned barplot for two categorical metrics
#'
#' The metric used for category_x will generate bars as columns on the x-axis.
#' The bar will be split vertically based on category_y. Each group in category_y will be aligned to make them easier to compare.
#'
#' @param meta A data.frame containing metadata
#' @param category_x A character object specifying the metadata to use for grouping on the x-axis
#' @param name_x A character object specifying a name to display on the x-axis
#' @param category_x A character object specifying the metadata to use for splitting in the y-direction
#' @param category_name A character object specifying a name to display for the colors
#' @param colorset_y A colorset to use as fills for category_y. Currently supported: "rainbow" or "varibow". Default is "varibow"
#' @param name_y A character object specifying a name for the y-axis.
#' @param padding A numeric object specifying the fraction of the total vertical space to use for separating category_y groups. Default is 0.2.
#'
#' @return a ggplot2 plot object
#' @export
qc_aligned_barplot <- function(meta,
                               category_x = "batch_id",
                               name_x = "Batch ID",
                               category_y = "well_id",
                               category_name = "Well ID",
                               colorset_y = "varibow",
                               name_y = "N Cells",
                               padding = 0.2) {
  
  assertthat::assert_that(sum(class(meta) %in% c("data.frame","data.table")) > 0)
  assertthat::assert_that(class(category_x) == "character")
  assertthat::assert_that(length(category_x) == 1)
  assertthat::assert_that(category_x %in% names(meta))
  assertthat::assert_that(class(name_x) == "character")
  assertthat::assert_that(length(name_x) == 1)
  assertthat::assert_that(class(category_y) == "character")
  assertthat::assert_that(length(category_y) == 1)
  assertthat::assert_that(category_y %in% names(meta))
  assertthat::assert_that(class(category_name) == "character")
  assertthat::assert_that(length(category_name) == 1)
  assertthat::assert_that(class(name_y) == "character")
  assertthat::assert_that(length(name_y) == 1)
  assertthat::assert_that(class(colorset_y) == "character")
  assertthat::assert_that(length(colorset_y) == 1)
  assertthat::assert_that(colorset_y %in% c("rainbow","varibow"))
  assertthat::assert_that(class(padding) == "numeric")
  assertthat::assert_that(length(padding) == 1)
  assertthat::assert_that(padding < 1)
  
  
  tidy_x <- rlang::parse_expr(category_x)
  tidy_y <- rlang::parse_expr(category_y)
  
  meta <- as.data.table(meta)
  count_table <- meta[, .(n_cells = nrow(.SD)),
                      by = mget(c(category_x, category_y))]
  
  plot_xpos <- data.frame(unique(count_table[[category_x]]))
  names(plot_xpos) <- category_x
  plot_xpos <- plot_xpos[order(plot_xpos[[category_x]]),,drop = FALSE]
  plot_xpos$xpos <- 1:nrow(plot_xpos)
  
  count_table <- count_table[plot_xpos, on = category_x]
  
  plot_fills <- data.frame(unique(count_table[[category_y]]))
  names(plot_fills) <- category_y
  if(colorset_y == "rainbow") {
    set.seed(3030)
    plot_fills$fill <- sample(grDevices::rainbow(nrow(plot_fills)), nrow(plot_fills))
  } else if(colorset_y == "varibow") {
    set.seed(3030)
    plot_fills$fill <- sample(BarMixer::varibow(nrow(plot_fills)), nrow(plot_fills))
  }
  plot_fills <- plot_fills[order(plot_fills[[category_y]]),]
  count_table <- count_table[plot_fills, on = category_y]
  
  group_maxes <- count_table[, .(group_max = max(n_cells)), by = list(get(category_y))]
  names(group_maxes)[1] <- category_y
  group_maxes <- group_maxes[order(get(category_y), decreasing = TRUE)]
  group_maxes <- group_maxes[, cum_max := cumsum(group_max)]
  group_maxes <- group_maxes[, group_center := cum_max - group_max / 2]
  group_maxes <- group_maxes[, padded_center := group_center + (max(cum_max) * (padding/nrow(group_maxes))) * (1:nrow(group_maxes) - 1)]
  group_maxes <- group_maxes[, padded_base := padded_center - group_max/2]
  group_maxes <- group_maxes[, padded_top := padded_center + group_max/2]
  
  count_table <- count_table[group_maxes, on = category_y]
  
  count_table <- count_table[order(get(category_y), decreasing = TRUE)]
  count_table <- count_table[, ymax := cumsum(n_cells), by = list(get(category_x))]
  count_table <- count_table[, ymin := shift(ymax, fill = 0, type = "lag"), by = list(get(category_x))]
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(data = count_table,
                       ggplot2::aes(xmin = xpos - 0.4,
                                    xmax = xpos + 0.4,
                                    ymin = padded_base,
                                    ymax = padded_base + n_cells,
                                    fill = fill)) +
    ggplot2::geom_hline(data = count_table,
                        ggplot2::aes(yintercept = padded_base)) +
    ggplot2::geom_hline(data = count_table,
                        ggplot2::aes(yintercept = padded_top),
                        linetype = "dashed") +
    ggplot2::scale_fill_identity(category_name,
                                 breaks = plot_fills$fill,
                                 labels = plot_fills[[category_y]],
                                 guide = "legend") +
    ggplot2::scale_x_continuous(name_x,
                                breaks = plot_xpos$xpos,
                                labels = plot_xpos[[category_x]]) +
    ggplot2::scale_y_continuous(name_y,
                                breaks = c(group_maxes$padded_base,
                                           group_maxes$padded_top),
                                labels = c(rep("", nrow(group_maxes)),
                                           group_maxes$group_max),
                                expand = ggplot2::expand_scale(c(0, 0.02))) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor.y = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 90,
                                                       hjust = 1,
                                                       vjust = 0.3))
  
  p
}


#' Generate a rainbow palette with variation in saturation and value
#'
#' @param n_colors The number of colors to generate
#'
#' @return a character vector of hex color values of length n_colors.
#' @export
#'
varibow <- function(n_colors) {
  
  assertthat::assert_that(is.numeric(n_colors))
  assertthat::assert_that(n_colors %% 1 == 0)
  assertthat::assert_that(length(n_colors) == 1)
  
  sats <- rep_len(c(0.55, 0.7, 0.85, 1), length.out = n_colors)
  vals <- rep_len(c(1, 0.8, 0.6), length.out = n_colors)
  sub("FF$", "", grDevices::rainbow(n_colors,
                                    s = sats,
                                    v = vals))
}
