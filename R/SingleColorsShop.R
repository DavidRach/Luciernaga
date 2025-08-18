AutofluoresceShop <- function(x, visualized, files, experimentdesignation="AB",
template, GatePlots=TRUE, TheN=3, Display="selection", AFOverlap=NULL, ExceptFor=NULL){

  Status <- x %in% visualized
  ExperimentName <- x
  if (Status == TRUE){return(Status)}

  internalfiles <- files[grep(ExperimentName, files)]
  internalfiles <- internalfiles[-grep("Beads", internalfiles)]
  internalfiles <- internalfiles[-grep("nstained", internalfiles)]

  if (!is.null(ExceptFor)){
    ExceptFor <- paste(ExceptFor, collapse="|")
    internalfiles <- internalfiles[-grep(ExceptFor, internalfiles)]
  }

  # experimentdesignation <- "AB"
  Experiment <- sub(paste0("_", experimentdesignation, ".*"), "", ExperimentName)
  Experiment <- gsub("_", "-", Experiment)

  LabFiles <- list.files("LabNotebook", include.dirs=TRUE)
  if (!Experiment %in% LabFiles)(stop("LabNotebook for ", Experiment, " is not present"))
  
  Notebook <- file.path("LabNotebook", Experiment)
  NotebookFiles <- list.files(Notebook, include.dirs=TRUE)

  if (!"SingleColors" %in% NotebookFiles){
      SingleColors <- file.path(Notebook, "SingleColors")
      dir.create(SingleColors)
  } else {SingleColors <- file.path(Notebook, "SingleColors")}

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
  fileName <- paste(fileName, "SingleColorGating", sep="_")

  Utility_Patchwork(x=Plots, filename=fileName, outfolder=SingleColors,
  thecolumns = 1, therows=1, returntype="pdf", NotListofList = FALSE,
  patches=TRUE)
  }

  if (is.null(AFOverlap)){
  FileLocation <- system.file("extdata", package = "Luciernaga")
  pattern = "AutofluorescentOverlaps.csv"
  AFOverlap <- list.files(path=FileLocation, pattern=pattern,
                          full.names = TRUE)
  } else {AFOverlap <- AFOverlap}

  Tags <- file.path(SingleColors, "Tags")
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

  Today <- Sys.Date()
  DataRow <- data.frame(Experiment=ExperimentName,
     Date=Today)
  
  return(DataRow)
  }