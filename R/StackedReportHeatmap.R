
#' Internal for StackedReport
#'
#' @param data The data intermediate of Stacked Report
#' @param legend Default is "right", use "none"
#' @param transpose Default is FALSE, flips orientation
#'
#' @importFrom ggplot2 ggplot aes geom_tile geom_text coord_fixed 
#'  element_blank theme_bw scale_fill_gradient theme element_blank 
#'  element_text element_line
#'
#' @noRd
StackedReportHeatmap <- function(data, nameColumn, legend, transpose){
  data$Ratio <- round(data$Ratio, 2)

  if (transpose == FALSE){
  plot <- ggplot(data, aes(x=.data[[nameColumn]], y = Cluster, fill = Ratio)) +
    geom_tile() + geom_text(aes(label = Ratio)) + theme_bw() +
    scale_fill_gradient(name = "Ratio", low = "#FFFFFF", high = "#FF0000",
   limits = c(0, NA)) + theme(legend.position = legend,
   plot.title = element_text(hjust = 0.5), panel.grid.minor = element_line(
      linetype = "blank"), axis.title = element_text(size = 10), axis.title.y = element_blank(),
      axis.title.x = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(),
      axis.text.x = element_text(angle = 40, hjust = 1), legend.key.size = unit(0.4, "cm"))  +
    coord_fixed(ratio = 1.1)
   } else{
    plot <- ggplot(data, aes(y=.data[[nameColumn]], x = Cluster, fill = Ratio)) +
      geom_tile() + geom_text(aes(label = Ratio)) + theme_bw() +
      scale_fill_gradient(name = "Ratio", low = "#FFFFFF", high = "#FF0000",
     limits = c(0, NA)) + theme(legend.position = legend,
     plot.title = element_text(hjust = 0.5), panel.grid.minor = element_line(
        linetype = "blank"), axis.title = element_text(size = 10), axis.title.y = element_blank(),
        axis.title.x = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 40, hjust = 1), legend.key.size = unit(0.4, "cm"))  +
      coord_fixed(ratio = 1.1)
   }

  return(plot)
  }