#' Visualizes Cluster MFI over time, uses Luciernaga_QC output 
#' 
#' @param data The LuciernagaQC data.frame style output with raw detector values. 
#' @param columnname The columnname for Cluster, becomes the categorical defining lines and dots
#' @param datename The column name for Experiment/Dates, becomes the x-axis
#' @param Detector Default NULL, specify a detector, overrides the y-axis
#' 
#' @importFrom stringr str_detect
#' @importFrom dplyr pull
#' @importFrom purrr map_chr
#' @importFrom dplyr arrange
#' @importFrom dplyr desc
#' @importFrom dplyr slice
#' @importFrom dplyr select
#' @importFrom tidyselect all_of
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes 
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme
#' @importFrom viridis scale_color_viridis
#' @importFrom viridis scale_fill_viridis
#' 
#' @return A ggplot2 object
#' 
#' @export
#' 
#' @examples
#' library(Luciernaga) 
#' 
Luciernaga_BrigtnessOverTime <- function(data, columnname="Cluster",
 datename="Experiment", Detector=NULL, plotType="line"){

  PeakDetectors <- data |> pull(columnname)
  PeakDetectors <- as.character(PeakDetectors)
  PeakDetectors <- PeakDetectors |> strsplit("_") |> map_chr(1)
  Peaks <- as.data.frame(table(PeakDetectors))
  Peaks <- Peaks |> arrange(desc(Freq))

  if (is.null(Detector)){
  colnames(Peaks)[1] <- "Detector"
  Detector <- Peaks |> slice(1) |> pull(Detector)
  Detector <- paste0(Detector, "-A")
  } else {Detector <- Detector
    if(!str_detect(Detector, "-A")){
      Detector <- paste0(Detector, "-A")
    }


  }

  TheData <- data |> select(all_of(c(columnname, datename, Detector)))

  if (plotType == "line"){
  Plot <- ggplot(TheData, aes(x = .data[[datename]], y = .data[[Detector]],
   group = .data[[columnname]], color = .data[[columnname]],
   fill = .data[[columnname]])) + geom_line(linewidth = 1) +
   geom_point(shape = 21, size = 3, stroke = 0.5) +
   labs(title = NULL, x = NULL, y = paste0("MFI: ", Detector)) + 
   theme_bw() + theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
   scale_color_viridis(discrete = TRUE) + scale_fill_viridis(discrete = TRUE)
  } else {
  Plot <- ggplot(TheData, aes(x = .data[[datename]], y = .data[[Detector]],
      group = .data[[columnname]], color = .data[[columnname]],
      fill = .data[[columnname]])) +
      geom_boxplot(alpha = 0.3, width = 0.6, outlier.shape = NA) +
      geom_point(shape = 21, size = 3, stroke = 0.5) +
      labs(title = NULL, x = NULL, y = paste0("MFI: ", Detector)) + 
      theme_bw() + theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
      scale_color_viridis(discrete = TRUE) + scale_fill_viridis(discrete = TRUE)
  }

  return(Plot)
  }
