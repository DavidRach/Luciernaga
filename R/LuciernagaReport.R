#' Converts the Luciernaga outputs into .pdf plots.
#'
#' @param LuciernagaData A data.frame
#' @param FluorophoreColumnName  The name of the data.frame column containing your fluorophores.
#' @param ClusterColumnName The name of the data.frame column containing your cluster IDs.
#' @param outfolder The location that you want to save the .pdf output to.
#' @param filename The name you want to save your .pdf file as.
#' @param returntype Passed to Utility_Patchwork for "pdf" or "patchwork"
#'
#' @importFrom dplyr group_by
#' @importFrom dplyr summarize
#' @importFrom dplyr left_join
#' @importFrom dplyr relocate
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom dplyr rename
#' @importFrom dplyr across
#' @importFrom tidyselect everything
#' @importFrom dplyr bind_rows
#' @importFrom dplyr bind_cols
#' @importFrom dplyr ungroup
#' @importFrom purrr map
#' @importFrom purrr flatten
#'
#' @return A value to be determined later
#' @export
#'
#' @examples NULL

LuciernagaReport <- function(data, RetainedType, CellPopRatio, outfolder, filename,
                             LinePlots=TRUE, CosinePlots=TRUE,
                             StackedBarPlots=TRUE, HeatmapPlots=TRUE,
                             returntype = "patchwork"){

  ################################################
  # Filtered by CellPopRatio, and creating other #
  ################################################

  TheCounts <- data %>% group_by(Sample, Experiment, Condition) %>%
    summarize(TotalCells = sum(Count, na.rm = TRUE), .groups = 'drop')

  TheData <- data %>% left_join(TheCounts, by = c("Sample", "Experiment", "Condition"))

  TheData <- TheData %>% mutate(Ratio = round(Count/TotalCells, 3)) %>% relocate(Ratio, .after=Count) %>%
    select(-TotalCells)

  FilteredData <- TheData %>% filter(Ratio > CellPopRatio)

  OtherData <- FilteredData %>% group_by(Sample, Experiment, Condition) %>%
    summarize(LostRatio = 1 - sum(Ratio, na.rm = TRUE), .groups = 'drop')

  Other <- TheCounts %>% left_join(OtherData, by = c("Sample", "Experiment", "Condition")) %>%
    mutate(Count = round(TotalCells*LostRatio, 0)) %>% select(-TotalCells) %>%
    relocate(Count, .before=LostRatio) %>% rename(Ratio=LostRatio) %>% mutate(Cluster="Other")

  OtherN <- nrow(Other)
  FirstDetectorColumn <- which(grepl("\\d", colnames(data)))[1]
  LastDetectorColumn <- tail(which(grepl("\\d", colnames(data))), 1)

  Replacement <- data %>% select(FirstDetectorColumn:LastDetectorColumn) %>% head(OtherN) %>%
    mutate(across(everything(), ~0))

  Replacements <- bind_cols(Other, Replacement) %>% ungroup()
  Replaced <- bind_rows(FilteredData, Replacements)

  ##############
  # Lets Begin #
  ##############

  Items <- data.frame(table(data$Sample)) %>% pull(Var1) %>% as.character(.)
  #x <- Items[1]
  #data <- Replaced

  ThePlots <- map(.x=Items, .f=InternalReport, data=Replaced,
                  FirstDetectorColumn=FirstDetectorColumn,
                  LastDetectorColumn=LastDetectorColumn,
                  RetainedType=RetainedType, CellPopRatio=CellPopRatio,
                  LinePlots=LinePlots, CosinePlots=CosinePlots,
                  StackedBarPlots=StackedBarPlots, HeatmapPlots=HeatmapPlots)

  if (returntype == "pdf"){
  Utility_Patchwork(x=ThePlots, filename = filename, outfolder = outfolder,
                    thecolumns = 2, therows = 2, width = 9, height = 7, returntype = "pdf",
                    NotListofList = FALSE)
  } else {
  Hey <-Utility_Patchwork(x=ThePlots, filename = filename, outfolder = outfolder,
                      thecolumns = 2, therows = 2, width = 9, height = 7, returntype = "patchwork",
                      NotListofList = FALSE)
  return(Hey)
  }

}



#' Internal for LuciernagaReport
#'
#' @param x Passed Sample for filtering
#' @param data The data.frame
#' @param FirstDetectorColumn A passed parameter
#' @param LastDetectorColumn A passed parameter
#' @param RetainedType Whether "raw" or "normalized" values
#' @param CellPopRatio Mininum cutoff for cluster size
#' @param LinePlots Whether to return LinePlots
#' @param CosinePlots Whether to return CosinePlots
#' @param StackedBarPlots Whether to return StackedBarPlots
#' @param HeatmapPlots Whether to return Heatmap Plots
#'
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom tidyr gather
#' @importFrom tidyselect all_of
#' @importFrom ggplot2 ggplot
#' @importFrom tidyselect where
#' @importFrom lsa cosine
#' @importFrom reshape2 melt
#' @importFrom dplyr pull
#' @importFrom dplyr mutate
#' @importFrom viridis scale_fill_viridis
#' @importFrom figpatch fig
#'
#' @noRd

InternalReport <- function(x, data, FirstDetectorColumn, LastDetectorColumn,
                           RetainedType, CellPopRatio, LinePlots, CosinePlots,
                           StackedBarPlots, HeatmapPlots){
  First <- FirstDetectorColumn+1
  Last <- LastDetectorColumn+1

  subset <- data %>% filter(Sample %in% c(x))
  colnames(subset) <- Luciernaga:::NameCleanUp(colnames(subset), removestrings="-A")

  #ZeroBuggedRows <- subset %>% filter(rowSums(select(.,
  #    all_of(First:Last)), na.rm = TRUE) == 0) %>% nrow(.)

  #if (ZeroBuggedRows > 0) {subset <- subset %>% filter(rowSums(select(
  #  ., all_of(First:Last)), na.rm = TRUE) != 0)}

  if (LinePlots == TRUE){
  LinePlotData <- subset %>% filter(!Cluster %in% "Other") %>%
    select(Cluster, {{First}}:{{Last}})

  LineColN <- ncol(LinePlotData)
  DetectorOrder <- colnames(subset)[First:Last]

  Melted <- LinePlotData %>% gather(key = "Detector", value = "value", all_of(2:LineColN))
  Melted$Detector <- factor(Melted$Detector, levels = DetectorOrder)
  Melted$Cluster <- factor(Melted$Cluster)

  if (RetainedType == "raw"){message("FIX RAW LOW AND HIGH")
                             Low <- 0
                             High <- 1.1
                             Entry <- "Normalized Values"}

  if (RetainedType == "normalized"){Low <- 0
                                    High <- 1.1
                                    Entry <- "Normalized Values"}

  LinePlot <- ggplot(Melted, aes(x = Detector, y = value, group = Cluster,
    color = Cluster)) + geom_line() + ylim(min = Low, max = High) +
    labs(title = "Fluorophores", x = "Detectors", y = Entry) +
    theme_bw() + scale_color_hue(direction = 1) + theme_linedraw() +
    theme(plot.title = element_text(size = 16L, face = "plain", hjust = 0.5),
          axis.title.y = element_text(size = 11L, face = "plain"),
          axis.title.x = element_text(size = 11L, face = "plain"),
          panel.grid.major = element_line(colour = "gray95", linetype = "twodash"),
          panel.grid.minor = element_line(colour = "gray95",linetype = "longdash"),
          panel.background = element_rect(fill = NA), plot.background = element_rect(
          colour = NA), legend.background = element_rect(fill = NA),
          axis.text.x = element_text(size = 5, angle = 45, hjust = 1))
  }

  if (CosinePlots == TRUE){
    CosineData <- subset %>% filter(!Cluster %in% "Other") %>%
      select(Cluster, {{First}}:{{Last}})
    Names <- CosineData$Cluster
    Numbers <- CosineData %>% select(where(is.numeric))
    NumericsT <- t(Numbers)
    rownames(NumericsT) <- NULL
    colnames(NumericsT) <- Names
    NumericsT <- data.matrix(NumericsT)

  if (ncol(NumericsT) >= 2){
    CosineMatrix <- cosine(NumericsT)
    CosineMatrix <- round(CosineMatrix, 2)
    Reordered <- ReorderedCosine(CosineMatrix)
    MeltedCosine <- melt(Reordered)

    #Generate a Red to Blue Heatmap
    CosinePlot <- ggplot(MeltedCosine, aes(Var2, Var1, fill = value)) +
      geom_tile(color = "white") +
      scale_fill_gradient2(low = "lightblue", high = "orange", mid = "white",
                           midpoint = 0.7, limit = c(0.4,1), space = "Lab",
                           name="Cosine\nSimilarity") +
      theme_bw() + geom_text(aes(Var2, Var1, label = value), color = "black",
                             size = 2) + coord_fixed(ratio = 1.3) +
      theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
            panel.grid.major = element_blank(), panel.border = element_blank(),
            panel.background = element_blank(), axis.ticks = element_blank(),
            legend.position.inside = c(1.2, 0.5),
            legend.direction = "vertical", axis.text.x = element_text(
              angle = 45, vjust = 1, hjust = 1, size = 6),
            axis.text.y = element_text(size = 6),
            legend.key.size = unit(0.4, "cm"))

    CosineOrder <- data.frame(table(MeltedCosine$Var1)) %>% pull(Var1) %>%
      as.character(.)

  } else {image_path <- system.file("hex", "hex.png", package = "Luciernaga", mustWork = TRUE)
          CosinePlot <- fig(image_path)
          }

  }

  Bd <- subset %>% mutate(Ratio = round(Ratio, 2))

  if (exists("CosineOrder")) {Bd$Cluster <- factor(Bd$Cluster,
      levels = unique(Bd$Cluster)[order(match(unique(Bd$Cluster), CosineOrder))])
  }

  if (StackedBarPlots == TRUE){
  title <- as.character(x)

  StackedBarPlot <- ggplot(Bd, aes(x= Sample, y = Ratio,
    fill = Cluster)) + geom_col() + theme_bw() + scale_fill_viridis(
    discrete = TRUE, option = "inferno", direction = -1) + labs(title = title) +
    theme(plot.title = element_text(hjust = 0.5), panel.grid.minor = element_line(
    linetype = "blank"), axis.title = element_text(size = 10),
    axis.title.x = element_blank(), legend.key.size = unit(0.4, "cm")) +
    coord_fixed(ratio = 2)
  }

  if (HeatmapPlots == TRUE) {

  HeatmapPlot <- ggplot(Bd, aes(x= Sample, y = Cluster, fill = Ratio)) +
    geom_tile() + geom_text(aes(label = Ratio)) + theme_bw() +
    scale_fill_gradient(name = "Ratio", low = "#FFFFFF", high = "#FF0000",
    limits = c(0, NA)) + theme(plot.title = element_text(hjust = 0.5),
    panel.grid.minor = element_line(linetype = "blank"), axis.title =
    element_text(size = 10), axis.title.y = element_blank(), axis.title.x =
    element_blank(), axis.line = element_blank(), axis.ticks = element_blank(),
    legend.key.size = unit(0.4, "cm"))  + coord_fixed(ratio = 1.1)
  }

  ThePlots <- list()
  if (LinePlots == TRUE){ThePlots <- append(ThePlots, list(LinePlot))}

  if (CosinePlots == TRUE){ThePlots <- append(ThePlots, list(CosinePlot))}

  if (StackedBarPlots == TRUE){ThePlots <- append(ThePlots, list(StackedBarPlot))}

  if (HeatmapPlots == TRUE){ThePlots <- append(ThePlots, list(HeatmapPlot))}

  return(ThePlots)
}

#' Internal for LuciernagaReport
#'
#'
#' @noRd

ReorderedCosine <- function(CosineMatrix){
  Day <- as.dist((1-CosineMatrix)/2)
  Night <- hclust(Day)
  Twilight <- CosineMatrix[Night$order, Night$order]
}

