#' Takes a Lucierna SingleColorQC .fcs, splits into percentiles, and plots the results.
#'
#' @param x A GatingSet object
#' @param subset The desired gating hierarchy level to look at the data
#' @param sample.name Keyword for which sample name is stored
#' @param removestrings A list of values to remove from the name
#' @param stats Whether to use "mean" or "median"
#' @param returntype Whether to return a "raw" or "normalized" value lineplot.
#' @param probsratio Ratio increments to break quantiles into, default is set to 0.1.
#' @param Whether to return "plot"
#'
#' @importFrom flowCore keyword
#' @importFrom flowWorkspace gs_pop_get_data
#' @importFrom BiocGenerics nrow
#' @importFrom flowCore exprs
#' @importFrom dplyr arrange
#' @importFrom dplyr desc
#' @importFrom dplyr filter
#' @importFrom dplyr slice
#' @importFrom dplyr pull
#' @importFrom stats quantile
#' @importFrom dplyr mutate
#' @importFrom dplyr group_by
#' @importFrom tidyr nest
#' @importFrom dplyr select
#' @importFrom tidyr unnest
#' @importFrom dplyr ungroup
#' @importFrom tidyr pivot_longer
#' @importFrom tidyselect all_of
#' @importFrom ggplot2 ggplot
#'
#' @return Either ggplots or the summarized data.frame object preceeding
#' @export
#'
#' @examples NULL
Luciernaga_LinearSlices <- function(x, subset, sample.name, removestrings, stats,
                                    returntype, probsratio=0.1, output){

  name <- keyword(x, sample.name)
  name <- NameCleanUp(name, removestrings)

  cs <- gs_pop_get_data(x, subset)
  startingcells <- nrow(cs)[[1]]
  Data <- data.frame(exprs(cs[[1]]), check.names = FALSE)
  n <- Data[,-grep("Time|FS|SC|SS|Original|W$|H$", names(Data))]
  DetectorOrder <- colnames(n)

  # Triangulating on Detector
  n[n < 0] <- 0
  A <- do.call(pmax, n)
  Normalized <- n/A
  #Normalized <- round(Normalized, 1)
  na_counts <- colSums(is.na(Normalized))
  Normalized[is.na(Normalized)] <- 0
  Counts <- colSums(Normalized == 1)
  PeakDetectorCounts <- data.frame(Fluors = names(Counts),
                                   Counts = Counts) %>% arrange(desc(Counts))
  rownames(PeakDetectorCounts) <- NULL
  Detectors <- PeakDetectorCounts %>% filter(Counts > 0)

  if(nrow(Detectors) > 1){message("Luciernaga_LinearSlices is only meant to work
                          on LuciernagaQC output .fcs files, your file ", name,
                          " contained two peak detectors, only the first was selected")
                          Detectors <- Detectors %>% slice(1)
  }

  TheDetector <- Detectors %>% pull(Fluors)
  data <- n

  #Assigning by Percentiles
  #probsratio <- 0.1
  probslabel <- probsratio*100

  percentiles <- quantile(data[[TheDetector]], probs = seq(0, 1, by = probsratio),
                          na.rm = TRUE)

  ByPercentiles <- data %>% mutate(Percentiles = cut(.data[[TheDetector]],
                      breaks = c(-Inf, percentiles),
                      labels = seq(0, 100, by = probslabel), include.lowest = TRUE))

  if (returntype == "normalized"){
    ByPercentiles <- Normalized %>% mutate(Percentiles = ByPercentiles$Percentiles)
  }

  Samples <- ByPercentiles %>% group_by(Percentiles) %>%
    nest(data = where(is.numeric)) %>%
    mutate(Averaged = map(data, .f=AveragedSignature, stats=stats)) %>%
    select(Percentiles, Averaged) %>% unnest(Averaged) %>% ungroup()

  if (output == "data"){
    return(Samples)
  }

  if (returntype == "normalized"){Expression <- "Normalized Values"
  } else {Expression <- "Raw MFI"}

  if (output == "plot"){
  LineCols <- ncol(Samples)
  Melted <- Samples %>%
    pivot_longer(all_of(2:LineCols),names_to = "Detector", values_to = "value")

  Melted$Detector <- factor(Melted$Detector, levels = DetectorOrder)

  plot <- ggplot(Melted, aes(x = Detector, y = value, group = Percentiles,
          color = Percentiles)) + geom_line() + scale_color_hue(direction = 1) +
          labs(title = name, x = "Detectors", y = Expression) +
          theme_linedraw() + theme_bw() + theme(axis.title.x = element_text(
          face = "plain"), axis.title.y = element_text(face = "plain"),
          axis.text.x = element_text(size = 5,
          angle = 45, hjust = 1), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

  return(plot)
  }
}
