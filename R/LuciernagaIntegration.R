#' Ease-of-life wrapper for Luciernaga_Summary, orchestrates passing of the
#' respective indices and gates to screen single-color unmixing controls
#' 
#' @param template Template designating single-color gates to screen
#' @param gs The GatingSet
#' @param AFOverlap See Luciernaga Vignette, default NULL falls back to the 
#' default shipped within Luciernaga extdata. 
#' @param externalAF_gs Default NULL, if you wish to provide an external AF for
#' subtraction, provide the gs name (and next two arguments). Leaving NULL will grab
#' the parent to provided gate on assumption both positive and negative cells will
#' be present in that gate. 
#' @param externalAF_gs_index The index value to select the unstained in combination
#'  with previous argument
#' @param externalAF_gs_gate The gate corresponding to the unstained population you
#' want to use for background subtraction
#' @param excludeThese Excludes columns containing values, leaving just the detector
#' columns. The default is set to "FSC|SSC|Time|-H|-W" 
#' @param inverse.transform Whether to inverse.transform, default is TRUE
#' @param GuessSimilar Default FALSE, will attempt to match to
#' similar fluorophores in the library
#' @param Unstained Default FALSE, set to TRUE if sample is unstained (and subtraction 
#' therefore is not needed)
#' 
#' @importFrom dplyr select filter pull
#' @importFrom utils read.csv
#' @importFrom flowWorkspace sampleNames
#' @importFrom purrr map2
#' 
#' @return A list of lists containing data.frame and plots for all the single-color
#' reference controls designated in the template. 
#' 
#' @export 
#' 
#' @examples A <- 2 + 2
#' 
LuciernagaIntegration <- function(template, gs, AFOverlap=NULL,
 externalAF_gs = NULL, externalAF_gs_index=NULL, externalAF_gs_gate=NULL,
 excludeThese="FSC|SSC|Time|-H|-W", inverse.transform=TRUE,
GuessSimilar=FALSE, Unstained=Unstained){

  if (is.null(AFOverlap)){
    FileLocation <- system.file("extdata", package = "Luciernaga")
    pattern = "AutofluorescentOverlaps.csv"
    AFOverlap <- list.files(path=FileLocation, pattern=pattern,
     full.names = TRUE)
    AFOverlap <- read.csv(AFOverlap, check.names=FALSE)
  } else {
    if(!is.data.frame(AFOverlap)){
      AFOverlap <- read.csv(AFOVerlap, check.names=FALSE)
      } else { # No Intervention needed
      }
  }

  These <- template |> select(name, Fluorophore, Detector)
  DetectorsPresent <- These |> filter(!is.na(Detector) & Detector != "")

  gs_index <- match(DetectorsPresent$name, sampleNames(gs))
  DetectorsPresent$gs_index <- gs_index

  GatesToAdd <- DetectorsPresent |> pull(Fluorophore)
  SpecimenIndeces <- DetectorsPresent |> pull(gs_index)

  # sampleNames(gs)
  # sampleNames(gs[SpecimenIndeces[14]])
  #x <- GatesToAdd[1]
  #y <- SpecimenIndeces[1]

  ListOfLists <- map2(.x=GatesToAdd, .y=SpecimenIndeces, .f=Luciernaga_Summary,
  gs=gs, externalAF_gs=externalAF_gs, externalAF_gs_index=externalAF_gs_index,
  externalAF_gs_gate=externalAF_gs_gate, GuessSimilar=GuessSimilar,
  NumberDetectors=NumberDetectors, excludeThese=excludeThese, Unstained=Unstained,
  inverse.transform=inverse.transform)

  return(ListOfLists)
}

#' Returns a summary of signature metrics for a gated population
#' 
#' @param x Iterated in gate name
#' @param y Iterated in specimen index
#' @param gs The GatingSet 
#' @param GuessSimilar Default FALSE, will attempt to match to
#' similar fluorophores in the library
#' @param NumberDetectors Default 64 (5-laser Cytek Aurora)
#' @param excludeThese Excludes columns containing values, leaving just the detector
#' columns. The default is set to "FSC|SSC|Time|-H|-W" 
#' @param Unstained Default FALSE, set to TRUE if sample is unstained (and subtraction 
#' therefore is not needed)
#' @param inverse.transform Whether to inverse.transform, default is TRUE
#' 
#' @importFrom  flowWorkspace gs_pop_get_parent
#' @importFrom flowCore exprs
#' @importFrom stringr str_detect
#' @importFrom dplyr select filter
#' 
#' @return A named list containing the various data.frames and plots for the respective gate
#' 
#' 
#' @export 
#' 
#' @examples A <- 2 + 2
#' 
Luciernaga_Summary <- function(x, y, gs,  
  externalAF_gs=NULL, externalAF_gs_index=NULL, externalAF_gs_gate=NULL,
  GuessSimilar=FALSE, NumberDetectors=64, excludeThese="FSC|SSC|Time|-H|-W", 
  Unstained=FALSE, inverse.transform=TRUE){

  if (is.null(externalAF_gs)){
    path <- gs_pop_get_parent(gs[y], x, inverse.transform=inverse.transform)
    parentgate <- basename(path)
    parentCS <- gs_pop_get_data(gs[y], parentgate)
  } else {
    parentCS <- gs_pop_get_data(externalAF_gs[externalAF_gs_index], externalAF_gs_gate)
  }

  parentData <- data.frame(exprs(parentCS[[1]]), check.names=FALSE)
  parentData <- parentData[!str_detect(colnames(parentData), excludeThese)]
  parentSignature <- AveragedSignature(parentData, stats="median")

  if (Unstained == TRUE){
    
    parentSignature <- NULL

    QCdata <- Luciernaga_QC(subsets=x, x=gs[y], AFOverlap=AFOverlap,
    experiment="Test", condition="Test", CellAF=parentSignature,
    Unstained = TRUE, inverse.transform=inverse.transform) 
  
  } else {
    QCdata <- Luciernaga_QC(x=gs[y], subsets=x, AFOverlap=AFOverlap,
    experiment="Test", condition="Test", CellAF=parentSignature,
    inverse.transform=inverse.transform) 
  }

  # LuciernagaQC outputs


  QCdata <- QCdata |> select(-Experiment, -Condition)
    
  ThePlotName <- QCdata[,1] |> unique() |> unname()

  QCPlot <- VisualizeSignatures(columnname="Cluster", characterColumns = "Cluster", 
    data=QCdata, Normalize=TRUE, plotname = ThePlotName)

  CosineData <- QCdata |> select(-Sample) 
  #Cutoff <- sum(CosineData$Count)*0.01
  #CosineData <- CosineData |> filter(Count > Cutoff)

  AmalgamateData <- QC_Amalgamate(data=CosineData, samplecolumn="Cluster",
   normalize=TRUE, countcolumn="Count", returnType="data",
    titlename=ThePlotName, linecolor="blue", legend=FALSE)

  JustAmalgamate <- AmalgamateData |> filter(Cluster %in% "Average")
  colnames(JustAmalgamate)[1] <- "Fluorophore"

  if(GuessSimilar == TRUE){

  SimilarData <- QC_WhatsThis(x="Average", columnname="Fluorophore",
   data=JustAmalgamate, NumberHits = 10, NumberDetectors = NumberDetectors,
   returnPlots = TRUE)

  TheSimilarData <- SimilarData[[1]]
  TheSimilarPlot <- SimilarData[[2]]
  
  } else {
    SimilarData <- NULL
    SimilarPlot <- NULL
  }

  NormalizedPlot <- QC_Amalgamate(data=CosineData, samplecolumn="Cluster",
   normalize=TRUE, countcolumn="Count", returnType="plot",
    titlename=ThePlotName, linecolor="blue", legend=FALSE)

  HeatmapData <- CosineData |> select(Cluster, Count) |>
     mutate(TotalCells = sum(Count, na.rm = TRUE)) |>
     mutate(Ratio = round(Count/TotalCells, 3)) |>
     mutate(Fluorophore=ThePlotName[[1]])

  HeatmapData$Fluorophore <-  factor(HeatmapData$Fluorophore) 

  HeatmapPlot <- StackedReportHeatmap(data=HeatmapData,
   nameColumn="Fluorophore", legend="right", transpose=FALSE)

  CosineData <- CosineData |> select(-Count) 

  CosinePlot <- Luciernaga_Cosine(data=CosineData,
  returntype="plot", rearrange=TRUE, limitlow=0.90, limithigh=1,
  colorlow="navajowhite", colorhigh="tan1", legend=TRUE) +
     labs(title="ThePlotName")

  # subset <- GatesToAdd[14]
  LinearData <- Luciernaga_LinearSlices(x=gs[y], subset=x,
   sample.name="TUBENAME", removestrings=".fcs", stats="median",
   returntype="normalized", output="data", parentSignature=parentSignature,
   inverse.transform=inverse.transform)

   LinearPlot <- Luciernaga_LinearSlices(x=gs[y], subset=x,
   sample.name="TUBENAME", removestrings=".fcs", stats="median",
   returntype="normalized", output="plot", parentSignature=parentSignature,
   inverse.transform=inverse.transform)

   ReturnThisList <- list(
    "LuciernagaQC_Data" = QCdata,
    "LuciernagaQC_Signatures"=QCPlot,
    "AveragedSignature_Data"=AmalgamateData,
    "LuciernagaQC_AmalgamatedPlot"=NormalizedPlot,
    "GuessSimilar_Data"=TheSimilarData,
    "GuessSimilar_Plot"=TheSimilarPlot,
    "ProportionPlot"=HeatmapPlot,
    "CosineSimilarityPlot"=CosinePlot,
    "LinearSlices_Data"=LinearData,
    "LinearSlices_Plot"=LinearPlot)

   return(ReturnThisList)
}