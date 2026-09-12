
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
#' @importFrom dplyr filter select pull mutate
#' @importFrom tidyr gather
#' @importFrom tidyselect all_of where
#' @importFrom lsa cosine
#' @importFrom reshape2 melt
#' @importFrom viridis scale_fill_viridis
#' @importFrom figpatch fig
#' @importFrom ggplot2 scale_color_hue theme_linedraw element_rect
#'  geom_tile scale_fill_gradient2 unit geom_col scale_fill_gradient
#'  coord_fixed geom_text aes theme_bw theme element_blank element_text
#'  ylim ggplot
#' @importFrom stats as.dist hclust
#' 
#' @return An internal value
#'
#' @noRd
InternalReport <- function(x, data, FirstDetectorColumn,
   LastDetectorColumn, RetainedType, CellPopRatio, LinePlots,
   CosinePlots,StackedBarPlots, HeatmapPlots){
  
  First <- FirstDetectorColumn+1
  Last <- LastDetectorColumn+1

  subset <- data |> filter(Sample %in% c(x))
  colnames(subset) <- NameCleanUp(colnames(subset), removestrings="-A")

  #ZeroBuggedRows <- subset %>% filter(rowSums(select(.,
  #    all_of(First:Last)), na.rm = TRUE) == 0) %>% nrow(.)

  #if (ZeroBuggedRows > 0) {subset <- subset %>% filter(rowSums(select(
  #  ., all_of(First:Last)), na.rm = TRUE) != 0)}

  if (LinePlots == TRUE){
  LinePlotData <- subset |> filter(!Cluster %in% "Other") %>%
    select(Cluster, {{First}}:{{Last}})

  LineColN <- ncol(LinePlotData)
  DetectorOrder <- colnames(subset)[First:Last]

  Melted <- LinePlotData |>
    gather(key = "Detector", value = "value", all_of(2:LineColN))
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
          panel.grid.major = element_line(colour = "gray95",
           linetype = "twodash"),
          panel.grid.minor = element_line(colour = "gray95",
          linetype = "longdash"),
          panel.background = element_rect(fill = NA),
           plot.background = element_rect(
          colour = NA), legend.background = element_rect(fill = NA),
          axis.text.x = element_text(size = 5, angle = 45, hjust = 1))
  }

  if (CosinePlots == TRUE){
    CosineData <- subset |> filter(!Cluster %in% "Other") %>%
      select(Cluster, {{First}}:{{Last}})
    Names <- CosineData$Cluster
    Numbers <- CosineData |> select(where(is.numeric))
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

    CosineOrder <- data.frame(table(MeltedCosine$Var1)) |> pull(Var1) %>%
      as.character(.)

  } else {image_path <- system.file("hex", "hex.png", package = "Luciernaga",
                        mustWork = TRUE)
          CosinePlot <- fig(image_path)
          }

  }

  Bd <- subset |> mutate(Ratio = round(Ratio, 2))

  if (exists("CosineOrder")) {Bd$Cluster <- factor(Bd$Cluster,
      levels = unique(Bd$Cluster)[order(
        match(unique(Bd$Cluster), CosineOrder))])
  }

  if (StackedBarPlots == TRUE){
  title <- as.character(x)

  StackedBarPlot <- ggplot(Bd, aes(x= Sample, y = Ratio,
    fill = Cluster)) + geom_col() + theme_bw() + scale_fill_viridis(
    discrete = TRUE, option = "inferno", direction = -1) +
    labs(title = title) +
    theme(plot.title = element_text(hjust = 0.5),
    panel.grid.minor = element_line(
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
  if (LinePlots == TRUE){
    ThePlots <- append(ThePlots, list(LinePlot))}

  if (CosinePlots == TRUE){
    ThePlots <- append(ThePlots, list(CosinePlot))}

  if (StackedBarPlots == TRUE){
    ThePlots <- append(ThePlots, list(StackedBarPlot))}

  if (HeatmapPlots == TRUE){
    ThePlots <- append(ThePlots, list(HeatmapPlot))}

  return(ThePlots)
}