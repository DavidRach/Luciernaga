#' Internal function to render 2D density plots for expanded MFI data
#' 
#' @param data A data.frame containing a non-numeric identifier column and detector columns
#' @param theme Theme selection: 'viridis', 'aurora', or 'bigfoot'
#' @param bins Granularity of the 2D binning. Default is 100.
#' @param plotname Optional title for the plot
#' @param logScale Logical. If TRUE (default), applies a log10 transformation to the y-axis. 
#'   If FALSE, keeps the y-axis linear.
#' 
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_bin2d scale_y_continuous theme_bw theme scale_fill_gradientn ggtitle xlab ylab element_blank element_text element_rect
#' @importFrom scales breaks_log
#'
CartoonDensity2D <- function(data, theme = "viridis", bins = 100, plotname = NULL, logScale = TRUE) {
  
  # Separate non-numeric identity column from numeric intensity columns
  id_col <- names(data)[!sapply(data, is.numeric)][1]
  detector_data <- data |> dplyr::select(tidyselect::where(is.numeric))
  
  # Pivot to long format for ggplot2 2D binning
  dat_long <- tidyr::pivot_longer(
    data, 
    cols = dplyr::all_of(names(detector_data)),
    names_to = "Detector",
    values_to = "Intensity"
  )
  
  # Preserve original detector ordering
  dat_long$Detector <- factor(dat_long$Detector, levels = unique(dat_long$Detector))
  
  # Lower bound threshold for log transform compatibility
  if (logScale) {
    dat_long$Intensity <- ifelse(dat_long$Intensity <= 0, 1, dat_long$Intensity)
  }
  
  title_str <- if (!is.null(plotname)) plotname else "MFI Spectral Density"
  
  # Build base plot
  p <- ggplot2::ggplot(dat_long, ggplot2::aes(x = Detector, y = Intensity))
  
  # Conditionally apply y-axis scale transformation
  if (logScale) {
    p <- p + ggplot2::scale_y_continuous(trans = "log10", breaks = scales::breaks_log(7))
  } else {
    p <- p + ggplot2::scale_y_continuous()
  }

  # Add geoms, theme, and scale layers
  p <- p +
    ggplot2::geom_bin2d(bins = bins) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::ggtitle(title_str) +
    ggplot2::xlab("Detector") +
    ggplot2::ylab(if (logScale) "Log10 Intensity" else "Intensity") +
    ggplot2::scale_fill_gradientn(
      colours = c("white", "blue", "lightblue", "green", "yellow", "red"),
      values = c(0, 0.1, 0.2, 0.3, 0.4, 1)
    )

  return(p)
}