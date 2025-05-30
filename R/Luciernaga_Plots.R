#' Converts the Luciernaga outputs into .pdf plots
#'
#' @param data The data.frame output from LuciernagaQC
#' @param RetainedType Whether the data.frame contains "raw" or "normalized" values
#' @param CellPopRatio What mininum ratio needed to retain cluster.
#' @param outfolder The location that you want to save the .pdf output to.
#' @param filename The name you want to save your .pdf file as.
#' @param LinePlots Passed to Utility_Patchwork for "pdf" or "patchwork" or "plots"
#' @param CosinePlots Return this kind of plot, default is set to TRUE
#' @param StackedBarPlots Return this kind of plot, default is set to TRUE
#' @param HeatmapPlots Return this kind of plot, default is set to TRUE
#' @param returntype Return "pdf", "patchwork" or "plots"
#' @param reference path or data.frame containing Fluorophore column for ordering
#' @param thecolumns The number of columns per page
#' @param therows The number of rows per page
#' @param width Desired page width
#' @param height Desired page height
#'
#' @importFrom dplyr group_by
#' @importFrom dplyr summarize
#' @importFrom dplyr left_join
#' @importFrom dplyr relocate
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom dplyr rename
#' @importFrom dplyr across
#' @importFrom dplyr select
#' @importFrom tidyselect all_of
#' @importFrom tidyselect everything
#' @importFrom dplyr bind_rows
#' @importFrom dplyr bind_cols
#' @importFrom dplyr ungroup
#' @importFrom purrr map
#' @importFrom purrr flatten
#' @importFrom utils head
#' @importFrom utils tail
#'
#' @return A value to be determined later
#' @export
#'
#' @examples
#'
#' library(flowCore)
#' library(flowWorkspace)
#' library(openCyto)
#' library(data.table)
#' library(dplyr)
#' library(purrr)
#' library(stringr)
#'
#' File_Location <- system.file("extdata", package = "Luciernaga")
#' FCS_Files <- list.files(path = File_Location, pattern = ".fcs",
#'   full.names = TRUE)
#' CellSingleColorFiles <- FCS_Files[grep("Cells", FCS_Files)]
#' CellSingleColors <- CellSingleColorFiles[!str_detect("Unstained", CellSingleColorFiles)]
#' MyCytoSet <- load_cytoset_from_fcs(CellSingleColors[1:2],
#'   truncate_max_range = FALSE,transformation = FALSE)
#' MyGatingSet <- GatingSet(MyCytoSet)
#' MyGates <- fread(file.path(path = File_Location, pattern = 'Gates.csv'))
#' MyGatingTemplate <- gatingTemplate(MyGates)
#' gt_gating(MyGatingTemplate, MyGatingSet)
#' removestrings <-  c("DR_", "Cells", ".fcs", "-", " ")
#' StorageLocation <- file.path("C:", "Users", "JohnDoe", "Desktop")
#'
#' FileLocation <- system.file("extdata", package = "Luciernaga")
#' pattern = "AutofluorescentOverlaps.csv"
#' AFOverlap <- list.files(path=FileLocation, pattern=pattern, full.names = TRUE)
#'
#' SingleColor_Data <- map(.x=MyGatingSet[1:2], .f=Luciernaga_QC, subsets="lymphocytes",
#'  removestrings=removestrings, sample.name="GUID", unmixingcontroltype = "cells",
#'  Unstained = FALSE, ratiopopcutoff = 0.001, Verbose = FALSE, AFOverlap = AFOverlap,
#'  stats = "median", ExportType = "data", SignatureReturnNow = FALSE,
#'  outpath = TemporaryFolder, Increments=0.1, SecondaryPeaks=2,
#'  experiment = "FirstExperiment", condition = "ILTPanel", Subtraction = "Internal",
#'  CellAF=TheCellAF, SCData="subtracted", NegativeType="default") %>% bind_rows()
#'
#' pattern = "^Panel.csv"
#' CSV <- list.files(path=FileLocation, pattern=pattern, full.names=TRUE)
#' TheFluorophoreOrder <- read.csv(CSV, check.names = FALSE)
#'
#' ThePlots <- Luciernaga_Plots(data=SingleColor_Data, RetainedType="normalized",
#'  CellPopRatio=0.05, outfolder=NULL, filename="LuciernagaReport", returntype="plots",
#'  LinePlots=FALSE, CosinePlots=FALSE, StackedBarPlots = FALSE, HeatmapPlots = TRUE,
#'  reference = TheFluorophoreOrder)
#'
Luciernaga_Plots <- function(data, RetainedType, CellPopRatio, outfolder, filename,
                             LinePlots=TRUE, CosinePlots=TRUE,
                             StackedBarPlots=TRUE, HeatmapPlots=TRUE,
                             returntype = "patchwork", reference=NULL,
                             thecolumns=2, therows=2, width=9, height=7){

  if(!is.null(reference)){

  if (!is.data.frame(reference)){reference <- read.csv(reference, check.names=FALSE)}
  
  PreferredOrder <- reference %>% pull(Fluorophore)
  PreferredOrder <- gsub("-A", "", PreferredOrder)
  } else {PreferredOrder <- NULL}

  #################################################
  # Filtered by CellPopRatio, and creating other  #
  #################################################

  TheCounts <- data %>% group_by(Sample, Experiment, Condition) %>%
    summarize(TotalCells = sum(Count, na.rm = TRUE), .groups = 'drop')

  TheData <- data %>% left_join(TheCounts, by = c("Sample", "Experiment", "Condition"))

  TheData <- TheData %>% mutate(Ratio = round(Count/TotalCells, 3)) %>%
    relocate(Ratio, .after=Count) %>% select(-TotalCells)

  FilteredData <- TheData %>% filter(Ratio > CellPopRatio)

  OtherData <- FilteredData %>% group_by(Sample, Experiment, Condition) %>%
    summarize(LostRatio = 1 - sum(Ratio, na.rm = TRUE), .groups = 'drop')

  Other <- TheCounts %>% left_join(OtherData, by = c("Sample", "Experiment",
    "Condition")) %>% mutate(Count = round(TotalCells*LostRatio, 0)) %>%
    select(-TotalCells) %>% relocate(Count, .before=LostRatio) %>%
    rename(Ratio=LostRatio) %>% mutate(Cluster="Other")

  OtherN <- nrow(Other)
  FirstDetectorColumn <- which(grepl("\\d", colnames(data)))[1]
  LastDetectorColumn <- tail(which(grepl("\\d", colnames(data))), 1)

  Replacement <- data %>% select(all_of(FirstDetectorColumn:LastDetectorColumn)) %>%
    head(OtherN) %>% mutate(across(everything(), ~0))

  Replacements <- bind_cols(Other, Replacement) %>% ungroup()
  Replaced <- bind_rows(FilteredData, Replacements)

  ##############
  # Lets Begin #
  ##############

  Items <- data.frame(table(data$Sample)) %>% pull(Var1) %>% as.character(.)

  if (!is.null(PreferredOrder)){
  if (all(Items %in% PreferredOrder)){
    Items <- PreferredOrder
  } else {message("names not matching, no reorderring according to panel order")}
  }

  #x <- Items[1]
  #data <- Replaced

  ThePlots <- map(.x=Items, .f=Luciernaga:::InternalReport, data=Replaced,
                  FirstDetectorColumn=FirstDetectorColumn,
                  LastDetectorColumn=LastDetectorColumn,
                  RetainedType=RetainedType, CellPopRatio=CellPopRatio,
                  LinePlots=LinePlots, CosinePlots=CosinePlots,
                  StackedBarPlots=StackedBarPlots, HeatmapPlots=HeatmapPlots)

  if (returntype == "pdf"){
  Utility_Patchwork(x=ThePlots, filename = filename, outfolder = outfolder,
                    thecolumns = thecolumns, therows = therows, width = width,
                    height = height, returntype = "pdf", NotListofList = FALSE)
  }

  if (returntype == "patchwork"){
  Hey <-Utility_Patchwork(x=ThePlots, filename = filename, outfolder = outfolder,
                      thecolumns = thecolumns, therows = therows, width = width,
                      height = height, returntype = "patchwork",
                      NotListofList = FALSE)
  return(Hey)
  }

  if (returntype == "plots"){
    return(ThePlots)
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
#' @importFrom ggplot2 scale_color_hue
#' @importFrom ggplot2 theme_linedraw
#' @importFrom ggplot2 element_rect
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 scale_fill_gradient2
#' @importFrom ggplot2 unit
#' @importFrom ggplot2 geom_col
#' @importFrom ggplot2 scale_fill_gradient
#' @importFrom ggplot2 coord_fixed
#' @importFrom ggplot2 geom_text
#' @importFrom stats as.dist
#' @importFrom stats hclust
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 ylim
#'
#' @return An internal value
#'
#' @noRd
InternalReport <- function(x, data, FirstDetectorColumn, LastDetectorColumn,
                           RetainedType, CellPopRatio, LinePlots, CosinePlots,
                           StackedBarPlots, HeatmapPlots){
  First <- FirstDetectorColumn+1
  Last <- LastDetectorColumn+1

  subset <- data %>% filter(Sample %in% c(x))
  colnames(subset) <- NameCleanUp(colnames(subset), removestrings="-A")

  #ZeroBuggedRows <- subset %>% filter(rowSums(select(.,
  #    all_of(First:Last)), na.rm = TRUE) == 0) %>% nrow(.)

  #if (ZeroBuggedRows > 0) {subset <- subset %>% filter(rowSums(select(
  #  ., all_of(First:Last)), na.rm = TRUE) != 0)}

  if (LinePlots == TRUE){
  LinePlotData <- subset %>% filter(!Cluster %in% "Other") %>%
    select(Cluster, {{First}}:{{Last}})

  LineColN <- ncol(LinePlotData)
  DetectorOrder <- colnames(subset)[First:Last]

  Melted <- LinePlotData %>%
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

  } else {image_path <- system.file("hex", "hex.png", package = "Luciernaga",
                        mustWork = TRUE)
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
#' @importFrom stats as.dist
#' @importFrom stats hclust
#'
#' @return An internal value
#'
#' @noRd
ReorderedCosine <- function(CosineMatrix){
  Day <- as.dist((1-CosineMatrix)/2)
  Night <- hclust(Day)
  Twilight <- CosineMatrix[Night$order, Night$order]
}

