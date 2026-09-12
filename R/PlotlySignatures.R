#' Internal for expt_parse
#'
#' @param data A data.frame of Fluorophore and n-detector columns
#'
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr pull
#' @importFrom tidyselect where
#' @importFrom ggplot2 ggplot aes geom_line labs theme_bw theme element_text
#' @importFrom plotly ggplotly
#'
#' @return An interactive plotly object of the line signatures
#' @noRd
PlotlySignatures <- function(data, TheFactor="Fluorophore"){
  Tidyed <- data %>% pivot_longer(cols = where(is.numeric),
                                  names_to = "Detector",
                                  values_to = "Value")

  DetectorOrder <- Tidyed %>% pull(Detector) %>% unique()
  Tidyed$Detector <- factor(Tidyed$Detector, levels=DetectorOrder)

  plot <- ggplot(Tidyed, aes(x = Detector, y = Value, color = .data[[TheFactor]], group = .data[[TheFactor]])) +
    geom_line() + labs(x = "Detector", y = "Normalized MFI", color = "Fluorophore") +
    theme_bw() + theme(axis.text.x = element_text(size=5, angle = 70, hjust = 1))

  plot <- plotly::ggplotly(plot)
  return(plot)
}