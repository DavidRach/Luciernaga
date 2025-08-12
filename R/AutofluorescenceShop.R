#' Handling function for AF processing and tag return for the AB LabNotebook
#' 
#' @param x An experiment folder, typically being iterated through
#' @param visualized A vector of experiment folders already processed
#' @param files A list of file.path to the Unstained raw .fcs files
#' @param experimentdesignation Default is "AB", used to identify
#' @param template A file.path to the openCyto gating template for our raw
#'  unstained .fcs files
#' @param GatePlots Default is TRUE, returns Utility_Gating plots for all
#' GatingSet files to verify gate placement worked as expected
#' @param TheN Selects the number of signature variants per peak detector
#' @param Display Default "selection" returns visual plots showing only TheN,
#'  alternatively "all" will show all signatures before filtering in the plots
#' 
#' @importFrom flowWorkspace load_cytoset_from_fcs GatingSet
#' @importFrom openCyto gatingTemplate gt_gating
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
#' @importFrom utils write.csv
#' @importFrom data.table fread
#' 
#' @return Selected .fcs tags to the Tags folder, visualized plot
#'  pdf and data .csv to the Autofluorescence folder. 
#' 
#' @export 
#' 
AutofluoresceShop <- function(x, visualized, files, experimentdesignation="AB",
template, GatePlots=TRUE, TheN=3, Display="selection"){

  Status <- x %in% visualized
  ExperimentName <- x
  if (Status == TRUE){return(Status)}

  internalfiles <- files[grep(ExperimentName, files)]

  # experimentdesignation <- "AB"
  Experiment <- sub(paste0("_", experimentdesignation, ".*"), "", ExperimentName)
  Experiment <- gsub("_", "-", Experiment)

  LabFiles <- list.files("LabNotebook", include.dirs=TRUE)
  if (!Experiment %in% LabFiles)(stop("LabNotebook for ", Experiment, " is not present"))
  
  Notebook <- file.path("LabNotebook", Experiment)
  NotebookFiles <- list.files(Notebook, include.dirs=TRUE)

  if (!"Autofluorescence" %in% NotebookFiles){
      Autofluorescence <- file.path(Notebook, "Autofluorescence")
      dir.create(Autofluorescence)
  } else {Autofluorescence <- file.path(Notebook, "Autofluorescence")}

  MyCytoSet <- load_cytoset_from_fcs(internalfiles,
   truncate_max_range = FALSE, transform = FALSE)
  MyGatingSet <- GatingSet(MyCytoSet)

  RawGates <- data.table::fread(template)
  RawGating <- gatingTemplate(RawGates)
  gt_gating(RawGating, MyGatingSet)

  if (GatePlots == TRUE){
  Plots <- purrr::map(.x=MyGatingSet, .f=Utility_GatingPlots,
   sample.name=c("GROUPNAME", "TUBENAME"),
   removestrings=c("Unmixed", "(", ")", ".fcs"),
   gtFile=RawGates, 
   outpath=NULL,
   returnType="patchwork",
   plotname=TRUE)

  fileName <- ExperimentName
  fileName <- paste(fileName, "AutofluorescenceGating", sep="_")

  Utility_Patchwork(x=Plots, filename=fileName, outfolder=Autofluorescence,
  thecolumns = 1, therows=1, returntype="pdf", NotListofList = FALSE,
  patches=TRUE)
  }

  FileLocation <- system.file("extdata", package = "Luciernaga")
  pattern = "AutofluorescentOverlaps.csv"
  AFOverlap <- list.files(path=FileLocation, pattern=pattern,
                          full.names = TRUE)

  Tags <- file.path(Autofluorescence, "Tags")
  if (!dir.exists(Tags)){dir.create(Tags)}

  ReturnedOutputs <- map(.x=MyGatingSet, .f=LuciernagaLocal,
   outpath=Tags, TheN=TheN, Display=Display,
   ExperimentName=ExperimentName, AFOverlap=AFOverlap)

  Dataset <- map(ReturnedOutputs, ~ .x$Data) |> bind_rows()
  ThePlots <- map(ReturnedOutputs, ~ .x$Plots)

  TheFileName <- paste(ExperimentName, "Signatures", sep="_")

  Utility_Patchwork(x=ThePlots, filename=TheFileName,
  outfolder=Autofluorescence, therows=3, thecolumns=1,
  NotListofList = FALSE)

  TheFileName <- paste(ExperimentName, "AFData", sep="_")
  TheFileName <- paste0(TheFileName, ".csv")
  StorageLocation <- file.path(Autofluorescence, TheFileName)

  write.csv(Dataset, StorageLocation, row.names=FALSE)
}


#' Internal for AutofluorescenceShop, parses individual .fcs file from
#' the GatingSet
#' 
#' @param x The iterated GatingSet object passed from AutofluorescenceShop
#' @param outpath Desired storage location for the AF Tags
#' @param ExperimentName The experiment name used by the respective folder
#' @param TheN Selects the number of signature variants per peak detector
#' @param Display Default "selection" returns visual plots showing only TheN,
#'  alternatively "all" will show all signatures before filtering in the plots
#' @param AFOverlap Luciernaga_QC default
#' 
#' @importFrom fs file_temp dir_ls file_copy dir_delete
#' @importFrom flowWorkspace keyword
#' @importFrom purrr map
#' @importFrom dplyr mutate group_by ungroup slice_head pull
#' arrange desc
#' 
#' @return A list containing Data and Plots, for subsequent splitting by 
#' the parent function
#' 
#' @noRd
LuciernagaLocal <- function(x, outpath, TheN, Display, ExperimentName,
AFOverlap){
  LuciernagaTemp <- file_temp("Luciernaga_Temp_")
  dir.create(LuciernagaTemp)
  # dir.exists(LuciernagaTemp)

  first <- keyword(x, "GROUPNAME")
  second <- keyword(x, "TUBENAME")
  LocalPlotName <- paste(first, second, sep="_")

  ReturnedFCS <- purrr::map(.x=x, .f=Luciernaga_QC, subsets="lymphocytes",
  removestrings=".fcs", 
  sample.name=c("GROUPNAME", "TUBENAME"), unmixingcontroltype = "cells",
  Unstained = TRUE, ratiopopcutoff = 0.01, Verbose = FALSE,
  AFOverlap = AFOverlap, stats = "median", ExportType = "fcs",
  SignatureReturnNow = FALSE, outpath = LuciernagaTemp, Increments=0.1,
  experiment=ExperimentName, condition="NA", minimalfcscutoff = 0.001, 
  NegativeType="artificial")

  TheTempFiles <- dir_ls(LuciernagaTemp, glob = "*.fcs")

  ReturnedFCS <- ReturnedFCS[[1]]

  DecisionData <- ReturnedFCS |>
     mutate(MainDetector = gsub("_.*", "", .data[["Cluster"]])) |>
     group_by(MainDetector) |> slice_head(n = TheN) |> ungroup()

  Clusters <- sub("_.*", "", DecisionData$Cluster)
  ClustersData <- data.frame(table(Clusters)) |>
     arrange(desc(Freq)) |> pull(Clusters)

  TheseSpecimens <- DecisionData |> mutate(Cluster = gsub("_", "", Cluster)) |>
    mutate(Cluster = gsub("-", "", Cluster)) |> 
    mutate(ID = paste(Sample, Cluster, sep="_")) |> pull(ID)

  TheseSpecimens <- paste0(TheseSpecimens, ".fcs")

  CopyThese <- TheTempFiles[basename(TheTempFiles) %in% TheseSpecimens]

  file_copy(CopyThese, outpath, overwrite=TRUE)

  if (Display=="selection"){
    DataToUse <- DecisionData
    TheLegend <- TRUE
    } else {
    DataToUse <- ReturnedFCS
    TheLegend <- FALSE}

  Plots <- purrr::map(.x=ClustersData, .f=SignatureVariants,
     data=DataToUse, returnType="plots", legend=TheLegend,
     plotname=LocalPlotName)

  dir_delete(LuciernagaTemp)
  #dir.exists(LuciernagaTemp)

  ReturnObjects <- list(Data=ReturnedFCS, Plots=Plots)
  return(ReturnObjects)
}

#' Internal for AutofluorescenceShop, filters down AF signatures to show
#' main variants as plots for later compiling via Patchwork
#' 
#' @param x The main detector being iterated on as grouping variable
#'  for the plot signatures
#' @param data The data return from Luciernaga_QC
#' @param returnType Whether to return plots, alternate is plotly
#' @param legend Whether to show legend or not
#' @param plotname Name to append to the plot for specimen identification
#' 
#' @importFrom dplyr select filter pull
#' @importFrom stringr str_detect
#' @importFrom plotly ggplotly
#' 
#' @return ggplot or plotly object
#' 
#' @noRd
SignatureVariants <- function(x, data, returnType, legend, plotname){
  x <- as.character(x)

  UnstainedSignature1 <- data |>
     dplyr::select(-Sample, -Experiment, -Condition, -Count)

  These <- UnstainedSignature1  |>
    dplyr::filter(stringr::str_detect(Cluster, paste0("^", x))) |>
    dplyr::pull(Cluster) |> unique()

  Plots <- Luciernaga::QC_ViewSignature(x=These,
     columnname="Cluster", data=UnstainedSignature1,
      Normalize=TRUE,TheFormat="wider", legend=legend,
  plotname=plotname)

  if (returnType == "plots"){
    return(Plots)
    } else {
      Plots <- plotly::ggplotly(Plots)
      return(Plots)
    }
}