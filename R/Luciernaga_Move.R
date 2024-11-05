#' Transfer the .fcs files selected by LuciernagaTree to a new location
#'
#' @param x A vector containing fluorophore names found within data.
#' @param data The data.frame output of LuciernagaTree listing files to be moved.
#' @param input The path to the current storage location of the .fcs files.
#' @param output The path to the desired future storage location of the selected
#' files.
#'
#' @importFrom dplyr filter
#' @importFrom stringr str_detect
#' @importFrom dplyr pull
#'
#' @return No return, .fcs files are moved desired folder.
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
#' StorageLocation <- file.path(tempdir(), "LuciernagaTreeExample")
#' if (!dir.exists(StorageLocation)) {dir.create(StorageLocation)}
#'
#' File_Location <- system.file("extdata", package = "Luciernaga")
#' FCS_Files <- list.files(path = File_Location, pattern = ".fcs",
#'   full.names = TRUE)
#' CellSingleColorFiles <- FCS_Files[grep("Cells", FCS_Files)]
#' CellSingleColors <- CellSingleColorFiles[!str_detect("Unstained", CellSingleColorFiles)]
#' MyCytoSet <- load_cytoset_from_fcs(CellSingleColors,
#'   truncate_max_range = FALSE,transformation = FALSE)
#' MyGatingSet <- GatingSet(MyCytoSet)
#' MyGates <- fread(file.path(path = File_Location, pattern = 'Gates.csv'))
#' MyGatingTemplate <- gatingTemplate(MyGates)
#' gt_gating(MyGatingTemplate, MyGatingSet)
#'
#' UnstainedFCSFiles <- FCS_Files[grep("Unstained", FCS_Files)]
#' UnstainedCells <- UnstainedFCSFiles[-grep("Beads", UnstainedFCSFiles)]
#' MyUnstainedCytoSet <- load_cytoset_from_fcs(UnstainedCells[1],
#'   truncate_max_range = FALSE,transformation = FALSE)
#' MyUnstainedGatingSet <- GatingSet(MyUnstainedCytoSet)
#' gt_gating(MyGatingTemplate, MyUnstainedGatingSet)
#'
#' removestrings <-  c(".fcs", "(", ")", "Cells")
#' FileLocation <- system.file("extdata", package = "Luciernaga")
#' pattern = "AutofluorescentOverlaps.csv"
#' AFOverlap <- list.files(path=FileLocation, pattern=pattern, full.names = TRUE)
#'
#' TheCellAF  <- Luciernaga_QC(x=MyUnstainedGatingSet[[1]], subsets="nonDebris",
#'  removestrings=removestrings, sample.name="GUID", unmixingcontroltype = "cells",
#'  Unstained = TRUE, ratiopopcutoff = 0.001, Verbose = FALSE,
#'  AFOverlap = AFOverlap, stats = "median", ExportType = "fcs",
#'  Brightness=TRUE, SignatureReturnNow = TRUE, outpath = StorageLocation,
#'  Increments=0.1, SecondaryPeaks=2, experiment = "FirstExperiment",
#'  condition = "ILTPanel", Subtraction = "Internal", CellAF=TheCellAF,
#'  SCData="subtracted",NegativeType="default", minimalfcscutoff=0.01)
#'
#'  MyGatingSet <- subset(MyGatingSet, !str_detect(name, "CD14"))
#'
#' SingleColor_Data <- map(.x=MyGatingSet, .f=Luciernaga_QC,
#'  subsets="lymphocytes", removestrings=removestrings, sample.name="GUID",
#'  unmixingcontroltype = "cells", Unstained = FALSE, ratiopopcutoff = 0.001,
#'  Verbose = FALSE, AFOverlap = AFOverlap, stats = "median", ExportType = "fcs",
#'  Brightness=TRUE, SignatureReturnNow = FALSE,outpath = StorageLocation,
#'  Increments=0.1, SecondaryPeaks=2, experiment = "FirstExperiment",
#'  condition = "ILTPanel", Subtraction = "Internal", CellAF=TheCellAF,
#'   SCData="subtracted",NegativeType="default", minimalfcscutoff=0.02)
#'
#' Unstained_Data  <- map(.x=MyUnstainedGatingSet, .f=Luciernaga_QC,
#'  subsets="nonDebris", removestrings=removestrings, sample.name="GUID",
#'  unmixingcontroltype = "cells", Unstained = TRUE, ratiopopcutoff = 0.001,
#'  Verbose = FALSE, AFOverlap = AFOverlap, stats = "median", ExportType = "fcs",
#'  Brightness=TRUE, SignatureReturnNow = FALSE, outpath = StorageLocation,
#'   Increments=0.1, SecondaryPeaks=2, experiment = "FirstExperiment",
#'   condition = "ILTPanel", Subtraction = "Internal", CellAF=TheCellAF,
#'    SCData="subtracted",NegativeType="default", minimalfcscutoff=0.01)
#'
#' PanelPath <- file.path(File_Location, "UnmixingPanel.csv")
#' MoveThese <- Luciernaga_Tree(BrightnessFilePath = StorageLocation,
#'  PanelPath = PanelPath)
#'
#' SortedStorageLocation <- file.path(tempdir(), "LuciernagaSelected")
#' if (!dir.exists(SortedStorageLocation)) {dir.create(SortedStorageLocation)}
#'
#' UnmixingPanel <- read.csv(PanelPath, check.names=FALSE)
#' TheseFluorophores <- UnmixingPanel %>% pull(Fluorophore)
#'
#' walk(.x=TheseFluorophores, .f=Luciernaga_Move, data=MoveThese,
#'  input=StorageLocation, output=SortedStorageLocation)
#' MovedFiles <- list.files(SortedStorageLocation, pattern="fcs", full.names=TRUE)
#'
Luciernaga_Move <- function(x, data, input, output){
  OriginalX <- x
  x <- gsub("-A", "", x)
  ### x <- gsub(" ", "", x)

  Internal <- data %>% dplyr::filter(str_detect(sample, fixed(x, ignore_case = TRUE)))

  if (x %in% c("PE", "APC")){

    Internal <- Internal %>% dplyr::filter(!str_detect(sample, "PE-|APC-|Per"))

    } #ExceptionHandling


  Fluor <- Internal %>% pull(sample) %>% unique()
  Clusterlet <- Internal %>% pull(Cluster) %>% unique()

  internalstrings <- c("_","-")
  Clusterlet <- Luciernaga::NameCleanUp(Clusterlet, removestrings = internalstrings)

  inputfiles <- list.files(input, full.names = TRUE)

  files_to_move <- inputfiles[str_detect(basename(inputfiles), fixed(Fluor, ignore_case = TRUE)) &
                                str_detect(basename(inputfiles), paste0("_", Clusterlet, "\\."))]

  file.copy(files_to_move, output)

}

