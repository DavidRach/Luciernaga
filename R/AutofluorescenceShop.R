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
#' @param AFOverlap Default NULL, alternately a file.path to an AFOverlap .csv file
#' @param nameOverride Default FALSE, if true uses all fcs files provided regardless
#' if they match the experiment name string character provided in x
#' @param subsets Default lymphocytes, alternatively provide gate name for Luciernaga_QC to
#' retrieve autofluorescence signatures from. 
#' @param therows Default 3
#' @param thecolumns Default 1
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
template, GatePlots=TRUE, TheN=3, Display="selection", AFOverlap=NULL,
 nameOverride=FALSE, subsets="lymphocytes", therows=3, thecolumns=1){

  Status <- x %in% visualized
  ExperimentName <- x
  if (Status == TRUE){return(Status)}

  if (nameOverride==FALSE){
    internalfiles <- files[grep(ExperimentName, files)]
  } else {internalfiles <- files}

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

  if (is.null(AFOverlap)){
  FileLocation <- system.file("extdata", package = "Luciernaga")
  pattern = "AutofluorescentOverlaps.csv"
  AFOverlap <- list.files(path=FileLocation, pattern=pattern,
                          full.names = TRUE)
  } else {AFOverlap <- AFOverlap}

  Tags <- file.path(Autofluorescence, "Tags")
  if (!dir.exists(Tags)){dir.create(Tags)}

  ReturnedOutputs <- map(.x=MyGatingSet, .f=LuciernagaLocal,
   outpath=Tags, TheN=TheN, Display=Display,
   ExperimentName=ExperimentName, AFOverlap=AFOverlap,
  subsets=subsets)

  Dataset <- map(ReturnedOutputs, ~ .x$Data) |> bind_rows()
  ThePlots <- map(ReturnedOutputs, ~ .x$Plots)

  TheFileName <- paste(ExperimentName, "Signatures", sep="_")

  Utility_Patchwork(x=ThePlots, filename=TheFileName,
  outfolder=Autofluorescence, therows=therows, thecolumns=thecolumns,
  NotListofList = FALSE)

  TheFileName <- paste(ExperimentName, "AFData", sep="_")
  TheFileName <- paste0(TheFileName, ".csv")
  StorageLocation <- file.path(Autofluorescence, TheFileName)

  write.csv(Dataset, StorageLocation, row.names=FALSE)

  Today <- Sys.Date()
  DataRow <- data.frame(Experiment=ExperimentName,
     Date=Today)
  
  return(DataRow)
  }