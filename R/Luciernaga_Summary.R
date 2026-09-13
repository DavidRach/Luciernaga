#' Returns a summary of signature metrics for a gated population
#' 
#' @param x Iterated in gate name
#' @param y Iterated in specimen index
#' @param gs The GatingSet 
#' @param externalAF_gs TBD
#' @param externalAF_gs_index TBD
#' @param externalAF_gs_gate TBD
#' @param GuessSimilar Default FALSE, will attempt to match to
#' similar fluorophores in the library
#' @param NumberDetectors Default 64 (5-laser Cytek Aurora)
#' @param excludeThese Excludes columns containing values, leaving just
#'  the detector
#' columns. The default is set to "FSC|SSC|Time|-H|-W" 
#' @param Unstained Default FALSE, set to TRUE if sample is unstained
#'  (and subtraction 
#' therefore is not needed)
#' @param inverse.transform Whether to inverse.transform, default is TRUE
#' @param AFOverlap See Luciernaga Vignette, default NULL falls back to the 
#' default shipped within Luciernaga extdata. 
#' 
#' @importFrom flowWorkspace gs_pop_get_parent gh_pop_get_indices
#'  gs_pop_get_data
#' @importFrom flowCore exprs
#' @importFrom stringr str_detect
#' @importFrom dplyr select filter mutate
#' @importFrom ggplot2 labs
#' 
#' @return A named list containing the various data.frames and plots for
#'  the respective gate
#' 
#' 
#' @export 
#' 
#' @examples A <- 2 + 2
#' 
Luciernaga_Summary <- function(x,
                                y,
                                gs,
                                externalAF_gs = NULL,
                                externalAF_gs_index = NULL,
                                externalAF_gs_gate = NULL,
                                GuessSimilar = FALSE,
                                NumberDetectors = 64,
                                excludeThese = "FSC|SSC|Time|-H|-W",
                                Unstained = FALSE,
                                inverse.transform = TRUE,
                                AFOverlap) {

  if (is.null(externalAF_gs)) {
    path <- gs_pop_get_parent(gs[y], x, inverse.transform = inverse.transform)
    parentgate <- basename(path)

    # Alternate assign to gh <- gs[[y]]
    parent_idx <- gh_pop_get_indices(gs[y], parentgate)
    # Alternate assign to gh <- gs[[y]]
    child_idx <- gh_pop_get_indices(gs[y], x)
    not_idx_in_parent <- !child_idx[parent_idx]

    parentCS <- gs_pop_get_data(gs[y], parentgate,
                                 inverse.transform = inverse.transform)
    TheseValues <- exprs(parentCS[[1]])
    TheseValues <- TheseValues[not_idx_in_parent, ]

  } else {
    parentCS <- gs_pop_get_data(externalAF_gs[externalAF_gs_index],
                                 externalAF_gs_gate,
                                 inverse.transform = inverse.transform)
    TheseValues <- exprs(parentCS[[1]])
  }

  parentData <- data.frame(TheseValues, check.names = FALSE)
  parentData <- parentData[!str_detect(colnames(parentData), excludeThese)]
  parentSignature <- AveragedSignature(parentData, stats = "median")

  # Stash <- parentSignature |> mutate(Fluorophore="Test") |>
  #   relocate(Fluorophore, .before=1)
  # VisualizeSignatures(Stash, 64, x="Test", columnname="Fluorophore")

  if (Unstained == TRUE) {

    parentSignature <- NULL

    QCdata <- Luciernaga_QC(subsets = x, x = gs[y], AFOverlap = AFOverlap,
                             experiment = "Test", condition = "Test",
                             CellAF = parentSignature, Unstained = TRUE,
                             inverse.transform = inverse.transform)

  } else {
    QCdata <- Luciernaga_QC(x = gs[y], subsets = x, AFOverlap = AFOverlap,
                             experiment = "Test", condition = "Test",
                             CellAF = parentSignature,
                             inverse.transform = inverse.transform)
  }

  # LuciernagaQC outputs

  QCdata <- QCdata |> select(-Experiment, -Condition)

  ThePlotName <- QCdata[, 1] |> unique() |> unname()

  QCPlot <- VisualizeSignatures(columnname = "Cluster",
                                 characterColumns = "Cluster", data = QCdata,
                                 Normalize = TRUE, plotname = ThePlotName)

  CosineData <- QCdata |> select(-Sample)
  # Cutoff <- sum(CosineData$Count)*0.01
  # CosineData <- CosineData |> filter(Count > Cutoff)

  AmalgamateData <- QC_Amalgamate(data = CosineData, samplecolumn = "Cluster",
                                   normalize = TRUE, countcolumn = "Count",
                                   returnType = "data",
                                   titlename = ThePlotName,
                                   linecolor = "blue", legend = FALSE)

  JustAmalgamate <- AmalgamateData |> filter(Cluster %in% "Average")
  colnames(JustAmalgamate)[1] <- "Fluorophore"

  if (GuessSimilar == TRUE) {

    SimilarData <- QC_WhatsThis(x = "Average", columnname = "Fluorophore",
                                 data = JustAmalgamate, NumberHits = 10,
                                 NumberDetectors = NumberDetectors,
                                 returnPlots = TRUE)

    TheSimilarData <- SimilarData[[1]]
    TheSimilarPlot <- SimilarData[[2]]

  } else {
    TheSimilarData <- NULL
    TheSimilarPlot <- NULL
  }

  NormalizedPlot <- QC_Amalgamate(data = CosineData, samplecolumn = "Cluster",
                                   normalize = TRUE, countcolumn = "Count",
                                   returnType = "plot",
                                   titlename = ThePlotName,
                                   linecolor = "blue", legend = FALSE)

  HeatmapData <- CosineData |>
    select(Cluster, Count) |>
    mutate(TotalCells = sum(Count, na.rm = TRUE)) |>
    mutate(Ratio = round(Count / TotalCells, 3)) |>
    mutate(Fluorophore = ThePlotName[[1]])

  HeatmapData$Fluorophore <- factor(HeatmapData$Fluorophore)

  HeatmapPlot <- StackedReportHeatmap(data = HeatmapData,
                                      nameColumn = "Fluorophore",
                                      legend = "right", transpose = FALSE)

  CosineData <- CosineData |> select(-Count)

  if (nrow(CosineData) > 1) {

    CosinePlot <- Luciernaga_Cosine(data = CosineData, returntype = "plot",
                                     rearrange = TRUE, limitlow = 0.90,
                                     limithigh = 1, colorlow = "navajowhite",
                                     colorhigh = "tan1", legend = TRUE) +
      labs(title = "ThePlotName")

  } else {
    message("Only one signature retrieved for ", x, " , Cosine not run")
    CosinePlot <- NULL
  }

  # subset <- GatesToAdd[14]
  LinearData <- Luciernaga_LinearSlices(x = gs[y], subset = x,
                                        sample.name = "TUBENAME",
                                        removestrings = ".fcs",
                                        stats = "median",
                                        returntype = "normalized",
                                        output = "data",
                                        parentSignature = parentSignature,
                                        inverse.transform = inverse.transform)

  LinearPlot <- Luciernaga_LinearSlices(x = gs[y], subset = x,
                                        sample.name = "TUBENAME",
                                        removestrings = ".fcs",
                                        stats = "median",
                                        returntype = "normalized",
                                        output = "plot",
                                        parentSignature = parentSignature,
                                        inverse.transform = inverse.transform)

  ReturnThisList <- list(
    "LuciernagaQC_Data" = QCdata,
    "LuciernagaQC_Signatures" = QCPlot,
    "AveragedSignature_Data" = AmalgamateData,
    "LuciernagaQC_AmalgamatedPlot" = NormalizedPlot,
    "GuessSimilar_Data" = TheSimilarData,
    "GuessSimilar_Plot" = TheSimilarPlot,
    "ProportionPlot" = HeatmapPlot,
    "CosineSimilarityPlot" = CosinePlot,
    "LinearSlices_Data" = LinearData,
    "LinearSlices_Plot" = LinearPlot)

  return(ReturnThisList)
}