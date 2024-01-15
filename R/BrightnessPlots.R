#' Visualize MFI of raw .fcs files to evaluate single color controls.
#'
#' @param thedata A data.frame with columns Fluorophore and Detector.
#' @param input The location where the .fcs files are stored.
#'
#' @return Visualized ggplots for each fluorophore. If multiple .fcs files of same fluorophore are
#'  present in the same folder, it overlays them.
#' @export
#'
#' @examples NULL

BrightnessPlots <- function(thedata, input){
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

  #thex <- Present[1]

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
      colnames(TheDF) <- gsub("-A$", "", colnames(TheDF))
      DFNames <- TheDF %>% select(all_of(TheDetector))
      DFNames <- DFNames %>% mutate(Cluster = filename)
      DFNames
    }

    TheDataFrames <- map(.x = cs, .f = InternalExprs2, thex2 = thex) %>% bind_rows()
    TheDataFrames$Cluster <- factor(TheDataFrames$Cluster)

    theXmin <- TheDataFrames[,1] %>% quantile(., 0.01)
    theXmax <- TheDataFrames[,1] %>% quantile(., 0.99)
    theXmin <- theXmin - abs((0.02*theXmin))
    theXmax <- theXmax + (0.02*theXmax)

    plot <- ggplot(TheDataFrames, aes(x =.data[[TheDetector]], fill = Cluster)) +
      geom_density(alpha = 0.5) + theme_bw() + coord_cartesian(xlim = c(theXmin, theXmax)) +
      labs(title = thex, x = TheDetector, y = "Frequency")

    theplotlist[[thex]] <- plot
  }

  PlotTwist <- map(.x = Present, .f = InternalExprs, data = InternalData, inputfiles = inputfiles)

  return(PlotTwist)
}
