#' Visualize cosine similarity of raw .fcs files to evaluate single color
#' controls.
#'
#' @param thedata A data.frame with columns Fluorophore and Detector.
#' @param input The location where the .fcs files are stored.
#' @param stats Whether to use the median or mean measurement for MFI
#'
#' @importFrom tidyr nest
#' @importFrom stringr str_detect
#' @importFrom flowCore exprs
#' @importFrom flowWorkspace keyword
#' @importFrom flowWorkspace load_cytoset_from_fcs
#' @importFrom dplyr select
#' @importFrom dplyr pull
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom dplyr relocate
#' @importFrom dplyr select_if
#' @importFrom dplyr bind_rows
#' @importFrom dplyr group_by
#' @importFrom dplyr ungroup
#' @importFrom dplyr summarise_all
#' @importFrom tidyr nest
#' @importFrom tidyr unnest
#' @importFrom purrr map
#' @importFrom lsa cosine
#' @importFrom reshape2 melt
#' @import ggplot2
#' @importFrom stats as.dist
#' @importFrom stats hclust
#'
#'
#'
#' @return Visualized ggplots for each fluorophore.
#'
#' @export
#'
#' @examples NULL

CosinePlots <- function(thedata, input, stats = NULL){
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

  #thex <- Present[20]

  InternalExprs <- function(thex, data, inputfiles){
    TheDetector <- data %>% filter(Fluorophore %in% thex) %>% pull(Detector)

    if (thex %in% c("PE", "APC")){thex <- paste0(thex, "_")}

    fcs_files <- inputfiles[str_detect(basename(inputfiles), thex) &
                              str_detect(basename(inputfiles), ".fcs$")]

    if (thex %in% c("PE_", "APC_")){thex <- gsub("_", "", thex)}

    cs <- load_cytoset_from_fcs(fcs_files, truncate_max_range = FALSE,
                                transformation = FALSE)

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
      DFNames <- TheDF %>% mutate(Cluster = filename) %>% relocate(Cluster,
          .before = 1)
      return(DFNames)
    }

    TheDataFrames <- map(.x = cs, .f = InternalExprs2, thex2 = thex) %>%
      bind_rows()

    #Removing the artificial negatives
    TheDataFrames <- TheDataFrames %>% mutate(Summed = rowSums(select_if(.,
              is.numeric), na.rm = TRUE))
    TheDataFrames <- TheDataFrames %>% filter(!Summed == 0) %>% select(-Summed)

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

    if(stats == "mean"){Samples <- LinePlotData %>% group_by(Cluster) %>%
      nest(data = where(is.numeric)) %>%
      mutate(mean_data = map(data, ~ summarise_all(., ~ round(mean(.,
          na.rm = TRUE),2)))) %>% select(Cluster, mean_data) %>%
      unnest(mean_data) %>% ungroup()
    } else if (stats == "median"){Samples <- LinePlotData %>%
      group_by(Cluster) %>%  nest(data = where(is.numeric)) %>%
      mutate(median_data = map(data, ~ summarise_all(., ~ round(median(.,
        na.rm = TRUE),2)))) %>% select(Cluster, median_data) %>%
      unnest(median_data) %>% ungroup()
    } else(print("NA"))

    if (!nrow(Samples) == 1){

      Names <- Samples$Cluster
      Numbers <- Samples %>% select_if(is.numeric)
      NumericsT <- t(Numbers)
      rownames(NumericsT) <- NULL
      colnames(NumericsT) <- Names

      NumericsT <- data.matrix(NumericsT)

      CosineMatrix <- cosine(NumericsT)
      CosineMatrix <- round(CosineMatrix, 2)

      # Reorder By Similarity
      reorder_cormat <- function(cormat){
        dd <- as.dist((1-cormat)/2)
        hc <- hclust(dd)
        cormat <-cormat[hc$order, hc$order]
      }

      # Reorder the correlation matrix
      cormat <- reorder_cormat(CosineMatrix)
      melted_cormat <- melt(cormat)

      #Generate a Heatmap
      plot <- ggplot(melted_cormat, aes(Var2, Var1, fill = value)) +
        geom_tile(color = "white") + geom_text(aes(Var2, Var1, label = value),
        color = "black", size = 2) +
        scale_fill_gradient2(low = "lightblue", high = "orange", mid = "white",
                             midpoint = 0.7, limit = c(0.4,1),
                             space = "Lab", name="Cosine\nSimilarity") +
        theme_bw() + theme(
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 25, hjust = 1)
        )

    theplotlist[[thex]] <- plot

    } else {
      LineCols <- ncol(Samples)

      Melted <- gather(Samples, key = "Detector", value = "value",
                       all_of(2:LineCols)) #Gather is my New Best Friend

      Melted$Detector <- factor(Melted$Detector, levels = detector_order)
      Melted$Cluster <- factor(Melted$Cluster)

      #Change this in case raw values provided instead of normalized?
      Low <- 0
      High <- 1.1

      Melted1 <- data.frame(Melted)

      plot <- ggplot(Melted1, aes(x = Detector, y = value, group = Cluster,
                                  color = Cluster)) + geom_line() +
        ylim(min = Low, max = High) + scale_color_hue(direction = 1) +
        labs(title = thex, x = "Detectors", y = "Normalized Values") +
        theme_linedraw() + theme_bw() +
        theme(axis.title.x = element_text(face = "plain"),
              axis.title.y = element_blank(),
              axis.text.x = element_text(size = 5, angle = 45, hjust = 1),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()
        )

      theplotlist[[thex]] <- plot
    }
  }

  PlotTwist <- map(.x = Present, .f = InternalExprs, data = InternalData,
                   inputfiles = inputfiles)

  return(PlotTwist)
}
