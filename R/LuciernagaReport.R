#' Converts the Luciernaga outputs into .pdf plots.
#'
#' @param LuciernagaData A data.frame
#' @param FluorophoreColumnName  The name of the data.frame column containing your fluorophores.
#' @param ClusterColumnName The name of the data.frame column containing your cluster IDs.
#' @param outfolder The location that you want to save the .pdf output to.
#' @param filename The name you want to save your .pdf file as.
#'
#' @importFrom dplyr pull
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#' @importFrom dplyr select_if
#' @importFrom tidyr gather
#' @import ggplot2
#' @importFrom lsa cosine
#' @importFrom reshape2 melt
#' @importFrom purrr flatten
#'
#'
#' @return A value to be determined later
#' @export
#'
#' @examples NULL
#'

LuciernagaReport <- function(data, FluorophoreColumnName, ClusterColumnName,
                             outfolder, filename, RetainedType, CellPopRatio){

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

  theplots <- map(.x=Items, .f=InternalReport, data=Replaced, FluorophoreColumnName=FluorophoreColumnName,
                  ClusterColumnName=ClusterColumnName, FirstDetectorColumn=FirstDetectorColumn,
                  LastDetectorColumn=LastDetectorColumn, RetainedType=RetainedType,
                  CellPopRatio=CellPopRatio)





   theflattenedplots <- purrr::flatten(theplots)

   Utility_Patchwork(x, filename=filename, outfolder=outfolder, thecolumns,
                     therows, width, weight, returntype)
}

InternalReport <- function(){
  First <- FirstDetectorColumn+1
  Last <- LastDetectorColumn+1

  subset <- data %>% filter(Sample %in% c(x))
  colnames(subset) <- Luciernaga:::NameCleanUp(colnames(subset), removestrings="-A")

  #ZeroBuggedRows <- subset %>% filter(rowSums(select(.,
  #    all_of(First:Last)), na.rm = TRUE) == 0) %>% nrow(.)

  #if (ZeroBuggedRows > 0) {subset <- subset %>% filter(rowSums(select(
  #  ., all_of(First:Last)), na.rm = TRUE) != 0)}

  LinePlotData <- subset %>% select(Cluster, {{First}}:{{Last}})

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

  Names <- LinePlotData$Cluster
  Numbers <- LinePlotData %>% select_if(is.numeric)
  NumericsT <- t(Numbers)
  rownames(NumericsT) <- NULL
  colnames(NumericsT) <- Names
  NumericsT <- data.matrix(NumericsT)

  if(!ncol(NumericsT) == 1){
    CosineMatrix <- lsa::cosine(NumericsT)
    CosineMatrix <- round(CosineMatrix, 2)

    # Reorder By Similarity
    reorder_cormat <- function(cormat){
      dd <- as.dist((1-cormat)/2)
      hc <- hclust(dd)
      cormat <-cormat[hc$order, hc$order]
    }

    # Reorder the correlation matrix
    cormat <- reorder_cormat(CosineMatrix)
    melted_cormat <- reshape2::melt(cormat)

    #Generate a Red to Blue Heatmap
    CosineHeatMap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value)) +
      geom_tile(color = "white") +
      scale_fill_gradient2(low = "lightblue", high = "orange", mid = "white",
                           midpoint = 0.7, limit = c(0.4,1), space = "Lab",
                           name="Cosine\nSimilarity") +
      theme_bw() + geom_text(aes(Var2, Var1, label = value), color = "black",
                             size = 2) + coord_fixed(ratio = 1.3) +
      theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
            panel.grid.major = element_blank(), panel.border = element_blank(),
            panel.background = element_blank(), axis.ticks = element_blank(),
            legend.position = c(1.2, 0.5),
            legend.direction = "vertical", axis.text.x = element_text(
              angle = 45, vjust = 1, hjust = 1, size = 6),
            axis.text.y = element_text(size = 6),
            legend.key.size = unit(0.4, "cm"))


    #print(CosineHeatMap)

    CosineOrder <- data.frame(table(melted_cormat$Var1)) %>% pull(Var1) %>%
      as.character(.)
  } else{CosineHeatMap <- NULL}

  #CosineOrder

  BarChartData <- LuciernagaSubset %>% mutate(Ratio = round(Ratio, 2))

  if(exists("CosineOrder")){BarChartData$Cluster <- factor(
    BarChartData$Cluster, levels = unique(BarChartData$Cluster)[
      order(match(unique(BarChartData$Cluster), CosineOrder))])}

  StackedBarChart <- ggplot(BarChartData, aes(x= sample, y = Ratio,
                                              fill = Cluster)) + geom_col() + theme_bw() + scale_fill_viridis(
                                                discrete = TRUE, option = "inferno", direction = -1) + labs(title = i) +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.minor = element_line(linetype = "blank"),
          axis.title = element_text(size = 10), axis.title.x = element_blank(),
          legend.key.size = unit(0.4, "cm")) + coord_fixed(ratio = 2)

  HeatMapChart <- ggplot(BarChartData, aes(x= sample, y = Cluster, fill = Ratio)) +
    geom_tile() + geom_text(aes(label = Ratio)) + theme_bw() +
    scale_fill_gradient(name = "Ratio", low = "#FFFFFF", high = "#FF0000",
                        limits = c(0, NA)) + theme(plot.title = element_text(hjust = 0.5),
                                                   panel.grid.minor = element_line(linetype = "blank"), axis.title =
                                                     element_text(size = 10), axis.title.y = element_blank(),
                                                   axis.title.x = element_blank(), axis.line = element_blank(),
                                                   axis.ticks = element_blank(), legend.key.size = unit(0.4, "cm"))  +
    coord_fixed(ratio = 1.1)

  theNegplots <- list(TheNormPlot, HeatMapChart, CosineHeatMap)
  #theNegplots <- list(TheNormPlot, StackedBarChart, HeatMapChart,
  #CosineHeatMap)

  theplots[[length(theplots) +1]] <- theNegplots
}

