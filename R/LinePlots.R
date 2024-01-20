#' Visualize normalized spectra of raw .fcs files to evaluate single color controls.
#'
#' @param thedata A data.frame with columns Fluorophore and Detector.
#' @param input The location where the .fcs files are stored.
#' @param stats Whether to use the median or mean measurement for MFI
#'
#' @importFrom tidyr nest
#'
#'
#' @return Visualized ggplots for each fluorophore.
#'
#' @export
#'
#' @examples NULL

LinePlots <- function(thedata, input, stats = NULL){
  data <- thedata
  data$Fluorophore <- gsub("-A$", "", data$Fluorophore)
  data$Fluorophore <- gsub(".", "", fixed = TRUE, data$Fluorophore)
  data$Fluorophore <- gsub("-", "", data$Fluorophore)
  data$Fluorophore <- gsub("_", "", data$Fluorophore)

  data$Detector <- gsub("-A$", "", data$Detector)
  data$Detector <- gsub(" ", "", data$Detector)

  InternalData <- data

  Variables <- InternalData %>% select(Fluorophore) %>% pull(.)

  inputfiles <- list.files(input, full.names = TRUE)

  #Present <- list()

  InternalList <- function(x, inputfiles){
    fcs_files <- inputfiles[str_detect(basename(inputfiles), x) &
                              str_detect(basename(inputfiles), ".fcs$")]
    if (length(fcs_files) > 0){return(x)}
  }

  Present <- map(Variables, InternalList, inputfiles = inputfiles)
  Present <- Filter(Negate(is.null), Present)
  Present <- unlist(Present)

  theplotlist <- list()

  #thex <- Present[2]

  InternalExprs <- function(thex, data, inputfiles){
    TheDetector <- data %>% filter(Fluorophore %in% thex) %>% pull(Detector)

    fcs_files <- inputfiles[str_detect(basename(inputfiles), thex) &
                              str_detect(basename(inputfiles), ".fcs$")]

    cs <- load_cytoset_from_fcs(fcs_files, truncate_max_range = FALSE, transform = FALSE)

    InternalExprs2 <- function(x, thex2){
      they <- x
      filename <- keyword(they, "FILENAME")
      filename <- sub(".*\\\\", "", filename)
      filename <- sub(paste0(".*", thex2), thex2, filename)
      filename <- gsub(".fcs$", "", filename)

      df <- exprs(they)
      TheDF <- data.frame(df, check.names = FALSE)
      TheDF <- TheDF[,-grep("Time|FS|SC|SS|Original|W$|H$", names(TheDF))]
      colnames(TheDF) <- gsub("-A$", "", colnames(TheDF))
      DFNames <- TheDF %>% mutate(Cluster = filename) %>% relocate(Cluster, .before = 1)
      return(DFNames)
    }

    TheDataFrames <- map(.x = cs, .f = InternalExprs2, thex2 = thex) %>% bind_rows()

    #Normalizing
    BackupNames <- colnames(TheDataFrames)
    BackupClusters <- TheDataFrames %>% select(Cluster)
    n <- TheDataFrames[,2:length(TheDataFrames)]

    #Normalizing By Peak Detector
    n[n < 0] <- 0
    A <- do.call(pmax, n)
    Normalized <- n/A
    Normalized <- round(Normalized, 1)
    detector_order <- colnames(Normalized)
    LinePlotData <- cbind(BackupClusters, Normalized)
    LinePlotData$Cluster <- factor(LinePlotData$Cluster)

    if(stats == "mean"){Samples <- LinePlotData %>% group_by(Cluster) %>%  nest(data = where(is.numeric)) %>%
      mutate(mean_data = map(data, ~ summarise_all(., ~ round(mean(., na.rm = TRUE),2)))) %>% select(Cluster, mean_data) %>% unnest(mean_data)
    } else if (stats == "median"){Samples <- LinePlotData %>% group_by(Cluster) %>%  nest(data = where(is.numeric)) %>%
      mutate(median_data = map(data, ~ summarise_all(., ~ round(median(., na.rm = TRUE),2)))) %>% select(Cluster, median_data) %>% unnest(median_data)
    } else(print("NA"))

    #Samples

    LineCols <- ncol(Samples)

    Melted <- gather(Samples, key = "Detector", value = "value", all_of(2:LineCols)) #Gather is my New Best Friend

    Melted$Detector <- factor(Melted$Detector, levels = detector_order)
    Melted$Cluster <- factor(Melted$Cluster)

    #Change this in case raw values provided instead of normalized?
    Low <- 0
    High <- 1.1

    Melted1 <- data.frame(Melted)

    plot <- ggplot(Melted1, aes(x = Detector, y = value, group = Cluster, color = Cluster)) + geom_line() + ylim(min = Low, max = High) +
      labs(title = "Fluorophores", x = "Detectors", y = "Normalized Values") + theme_bw() + scale_color_hue(direction = 1) + theme_linedraw() +
      theme(plot.title = element_text(size = 16L, face = "plain", hjust = 0.5), axis.title.y = element_text(size = 11L, face = "plain"),
            axis.title.x = element_text(size = 11L, face = "plain"), panel.grid.major = element_line(colour = "gray95", linetype = "twodash"),
            panel.grid.minor = element_line(colour = "gray95", linetype = "longdash"), panel.background = element_rect(fill = NA),
            plot.background = element_rect(colour = NA), legend.background = element_rect(fill = NA),
            axis.text.x = element_text(size = 5, angle = 45, hjust = 1))

    theplotlist[[thex]] <- plot
  }

  PlotTwist <- map(.x = Present, .f = InternalExprs, data = InternalData, inputfiles = inputfiles)

  return(PlotTwist)
}
