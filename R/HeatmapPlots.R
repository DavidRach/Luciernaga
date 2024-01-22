#' Returns a heatmap of subclusters proportion by Fluorophore.
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

HeatmapPlots <- function(thedata, input, stats = NULL){
  data <- thedata
  data$Fluorophore <- gsub("-A$", "", data$Fluorophore)
  #data$Fluorophore <- gsub(".", "", fixed = TRUE, data$Fluorophore)
  data$Fluorophore <- gsub("-", "", data$Fluorophore)
  data$Fluorophore <- gsub("_", "", data$Fluorophore)
  data$Fluorophore <- gsub(" ", "", data$Fluorophore)

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

    if (thex %in% c("PE", "APC")){thex <- paste0(thex, "_")}

    fcs_files <- inputfiles[str_detect(basename(inputfiles), thex) &
                              str_detect(basename(inputfiles), ".fcs$")]

    if (thex %in% c("PE_", "APC_")){thex <- gsub("_", "", thex)}

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

    #Removing the artificial negatives
    TheDataFrames <- TheDataFrames %>% mutate(Summed = rowSums(select_if(., is.numeric), na.rm = TRUE))
    TheDataFrames <- TheDataFrames %>% filter(!Summed == 0) %>% select(-Summed)

    TableData <- data.frame(table(TheDataFrames$Cluster))
    colnames(TableData)[1] <- "Cluster"

    BarChartData <- TableData %>% mutate(Ratio = round(Freq/sum(Freq), 2)) %>% mutate(sample = "")

    plot <- ggplot(BarChartData, aes(x= sample, y = Cluster, fill = Ratio)) + geom_tile() + geom_text(aes(label = Ratio)) +
      scale_fill_gradient(name = "Ratio", low = "#FFFFFF", high = "#FF0000", limits = c(0, NA)) + coord_fixed(ratio = 1.1) +
      theme_bw() + theme(axis.title.x = element_blank(),
                         axis.title.y = element_blank(),
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank()
      )

    theplotlist[[thex]] <- plot
  }

  PlotTwist <- map(.x = Present, .f = InternalExprs, data = InternalData, inputfiles = inputfiles)

  return(PlotTwist)
}
