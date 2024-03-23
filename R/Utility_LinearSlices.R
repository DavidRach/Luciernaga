#' Takes a Lucierna SingleColorQC .fcs, splits into percentiles, and plots the results.
#'
#' @param x A GatingSet object
#' @param subset The desired gating hierarchy level to look at the data
#' @param sample.name Keyword for which sample name is stored
#' @param removestrings A list of values to remove from the name
#' @param TheReturn Whether to return a "raw" or "normalized" value lineplot.
#'
#' @importFrom BiocGenerics nrow
#' @importFrom flowCore exprs
#' @importFrom flowWorkspace keyword
#' @importFrom flowWorkspace gs_pop_get_data
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom dplyr arrange
#' @importFrom dplyr pull
#' @importFrom dplyr mutate
#' @importFrom dplyr summarise_all
#' @importFrom dplyr group_by
#' @importFrom dplyr ungroup
#' @importFrom tidyr nest
#' @importFrom tidyr unnest
#' @importFrom tidyr gather
#' @importFrom ggplot2 ggplot
#'
#'
#' @return A ggplot2 object of lineplots from the sliced .fcs file
#' @export
#'
#' @examples Not at this time

Utility_LinearSlices <- function(x, subset, sample.name, removestrings, TheReturn){
  name <- keyword(x, sample.name)
  name <- NameCleanUp(name, removestrings)

  cs <- gs_pop_get_data(x, subset)
  startingcells <- nrow(cs)[[1]]
  Data <- exprs(cs[[1]])
  Data <- data.frame(Data, check.names = FALSE)

  TheColumns <- Data[,-grep("Time|FS|SC|SS|Original|W$|H$", names(Data))]
  detector_order <- colnames(TheColumns)

  # Triangulating on Detector
  n <- TheColumns
  n[n < 0] <- 0
  A <- do.call(pmax, n)
  Normalized <- n/A
  Normalized <- round(Normalized, 1)
  na_counts <- colSums(is.na(Normalized))
  Normalized[is.na(Normalized)] <- 0
  Counts <- colSums(Normalized == 1)
  PeakDetectorCounts <- data.frame(Fluors = names(Counts), Counts = Counts)
  rownames(PeakDetectorCounts) <- NULL
  cutoff <- startingcells*0.0075
  Detectors <- PeakDetectorCounts %>% filter(Counts > cutoff) %>%
    arrange(desc(Counts))
  TheDetector <- Detectors[1,1]


  if (TheReturn == "raw"){
    #Assigning by Percentiles
    percentiles <- quantile(TheColumns[[TheDetector]], probs = seq(0, 1, by = 0.1), na.rm = TRUE)

    TheColumns1int <- TheColumns %>% mutate(Percentiles = cut(.data[[TheDetector]],
                                                              breaks = c(-Inf, percentiles), labels = seq(0, 100, by = 10),
                                                              include.lowest = TRUE))

    TheColumns1 <- TheColumns1int

  } else if (TheReturn == "normalized"){
    percentiles <- quantile(TheColumns[[TheDetector]], probs = seq(0, 1, by = 0.1), na.rm = TRUE)

    TheColumns1int <- TheColumns %>% mutate(Percentiles = cut(.data[[TheDetector]],
                                                              breaks = c(-Inf, percentiles), labels = seq(0, 100, by = 10),
                                                              include.lowest = TRUE))

    TheColumns1 <- Normalized %>% mutate(Percentiles = TheColumns1int$Percentiles)
  }

  TheValues <- table(TheColumns1$Percentiles) %>% data.frame() %>% pull(Var1)
  TheValues <- TheValues %>% as.character()

  if(stats == "mean"){Samples <- TheColumns1 %>% group_by(Percentiles) %>%
    nest(data = where(is.numeric)) %>% mutate(mean_data = map(data, ~ summarise_all(., ~ round(mean(., na.rm = TRUE),2)))) %>% select(Percentiles, mean_data) %>% unnest(mean_data) %>% ungroup()
  } else if (stats == "median"){Samples <- TheColumns1 %>% group_by(Percentiles) %>%  nest(data = where(is.numeric)) %>% mutate(median_data = map(data, ~ summarise_all(., ~ round(median(., na.rm = TRUE),2)))) %>% select(Percentiles, median_data) %>% unnest(median_data) %>% ungroup()} else(message("NA"))

  LineCols <- ncol(Samples)

  Melted <- gather(Samples, key = "Detector", value = "value", all_of(
    2:LineCols)) #Gather is my New Best Friend

  Melted$Detector <- factor(Melted$Detector, levels = detector_order)
  Melted$Percentiles <- factor(Melted$Percentiles)

  #Change this in case raw values provided instead of normalized?
  #Low <- 0
  #High <- 1.1

  Melted1 <- data.frame(Melted)

  plot <- ggplot(Melted1, aes(x = Detector, y = value, group = Percentiles,
                              color = Percentiles)) + geom_line() +
    scale_color_hue(direction = 1) +
    labs(title = TheDetector, x = "Detectors", y = "Normalized Values") +
    theme_linedraw() + theme_bw() +
    theme(axis.title.x = element_text(face = "plain"),
          axis.title.y = element_text(face = "plain",
                                      margin = margin(r = -120)),
          axis.text.x = element_text(size = 5, angle = 45, hjust = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
    )

  return(plot)
}
