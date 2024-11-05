
#' Select candidate Luciernaga output .fcs files for future use in unmixing.
#'
#' @param BrightnessFilePath A path to the folder the individual brightness files are in.
#' @param PanelPath A path to the .csv file containing the panel information. It should include only
#' the fluorophores captured by the BrightnessFiles
#'
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
#' @importFrom utils read.csv
#' @importFrom dplyr pull
#' @importFrom dplyr filter
#' @importFrom dplyr arrange
#' @importFrom dplyr mutate
#' @importFrom dplyr pull
#' @importFrom dplyr row_number
#' @importFrom dplyr relocate
#'
#' @return A data.frame listing the candidate .fcs files for future unmixing use.
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
#' MyCytoSet <- load_cytoset_from_fcs(CellSingleColors[1:2],
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
#' SingleColor_Data <- map(.x=MyGatingSet[1:2], .f=Luciernaga_QC,
#'  subsets="lymphocytes", removestrings=removestrings, sample.name="GUID",
#'  unmixingcontroltype = "cells", Unstained = FALSE, ratiopopcutoff = 0.001,
#'  Verbose = FALSE, AFOverlap = AFOverlap, stats = "median", ExportType = "fcs",
#'  Brightness=TRUE, SignatureReturnNow = FALSE,outpath = StorageLocation,
#'  Increments=0.1, SecondaryPeaks=2, experiment = "FirstExperiment",
#'  condition = "ILTPanel", Subtraction = "Internal", CellAF=TheCellAF,
#'   SCData="subtracted",NegativeType="default")
#'
#' Unstained_Data  <- map(.x=MyUnstainedGatingSet[1], .f=Luciernaga_QC,
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
Luciernaga_Tree <- function(BrightnessFilePath, PanelPath){

  TheCSVs <- list.files(BrightnessFilePath, pattern="RelativeBrightness", full.names = TRUE)
  TheData <- map(.x=TheCSVs, .f=CSVRead) %>% bind_rows()

  if(!is.data.frame(PanelPath)){Panel <- read.csv(PanelPath, check.names = FALSE)
  } else {Panel <- PanelPath}

  OriginalPanel <- Panel
  #Panel$Fluorophore <- gsub(" ", "", Panel$Fluorophore)

  TheFluorophores <- Panel %>% pull(Fluorophore)
  #TheFluorophores <- TheFluorophores[c(4, 16, 25)]
  TheFluorophores <- gsub("-A", "", TheFluorophores)

  #x <- TheFluorophores[19]

  NewData <- map(.x=TheFluorophores, .f=InternalTree, TheData=TheData) %>% bind_rows()

  return(NewData)
}

#' Internal for Luciernaga_Tree
#'
#' @param x The iterated Fluorophore to be filtered
#' @param TheData The Data for all fluorophores
#'
#' @importFrom dplyr filter
#' @importFrom stringr str_detect
#' @importFrom dplyr mutate
#' @importFrom dplyr relocate
#' @importFrom dplyr arrange
#' @importFrom dplyr desc
#' @importFrom dplyr pull
#' @importFrom dplyr row_number
#'
#' @return An internal value
#'
#' @noRd
InternalTree <- function(x, TheData){
  OriginalX <- x

  Internal <- TheData %>% dplyr::filter(str_detect(sample, fixed(x, ignore_case = TRUE)))

  if (x %in% c("PE", "APC")){

    Internal <- Internal %>% dplyr::filter(!str_detect(sample, "PE-|APC-|Per"))
  } #ExceptionHandling


  Total <- sum(Internal$Count, na.rm = TRUE)
  Internal <- Internal %>% mutate(Ratio = Count / Total) %>% relocate(Ratio, .after=Count)

  if(nrow(Internal) == 0){message("No Fluorophore for ", x, " was found")}

  if(nrow(Internal)>1){
    Internal <- Internal %>% arrange(desc(Detector1Raw)) #First Arrange
    MaxVal <- Internal %>% filter(row_number() == 1) %>% pull(Detector1Raw)
    Internal1 <- Internal %>% filter(Detector1Raw > MaxVal*0.8)

    if (nrow(Internal1)>1){Abundance <- Internal1 %>%
      filter(row_number() == 1) %>%  pull(Ratio)

    if(Abundance < 0.5){Internal2 <- Internal1 %>% arrange(desc(Ratio))
    #Second Arrange
    Proportion <- Internal2 %>% select(Ratio) %>% sum(.)
    Top <- Internal2 %>% filter(row_number() == 1) %>% pull(Ratio)

    if((Top/Proportion) < 0.5){MainComplexity <- Internal2 %>%
      filter(row_number() == 1) %>% pull(Brightness)
    Complexity <- Internal2 %>% arrange(desc(Brightness)) %>% pull(Brightness)
    ComplexityLowerBound <- Complexity[1]*0.6

    if(MainComplexity > ComplexityLowerBound){Internal3 <- Internal2 %>%
      arrange(Brightness)
    SubsetData <- Internal3 %>% filter(row_number() == 1)
    SubsetData <- SubsetData %>% mutate(Decision = "Contested") %>% relocate(
      Decision, .after = Cluster)
    return(SubsetData)

    } else {SubsetData <- Internal2 %>% filter(row_number() == 1)
    SubsetData <- SubsetData %>% mutate(Decision = "Fifth Level") %>% relocate(
      Decision, .after = Cluster)
    return(SubsetData)}
    } else {SubsetData <- Internal2 %>% filter(row_number() == 1)
    SubsetData <- SubsetData %>% mutate(Decision = "Fourth Level") %>% relocate(
      Decision, .after = Cluster)
    return(SubsetData)}
    } else {SubsetData <- Internal1 %>% filter(row_number() == 1)
    SubsetData <- SubsetData %>% mutate(Decision = "Third Level") %>% relocate(
      Decision, .after = Cluster)
    return(SubsetData)}
    } else {SubsetData <- Internal1 %>% filter(row_number() == 1)
    SubsetData <- SubsetData %>% mutate(Decision = "Second Level") %>% relocate(
      Decision, .after = Cluster)
    return(SubsetData)}
  } else {SubsetData <- Internal %>% filter(row_number() == 1)
  SubsetData <- SubsetData %>% mutate(Decision = "First Level") %>% relocate(
    Decision, .after = Cluster)
  return(SubsetData)}
}

#' Internal for Luciernaga_Tree
#'
#' @param x An iterated path to a .csv to be read.
#'
#' @param utils read.csv
#' @param dplyr mutate
#' @param dplyr relocate
#'
#' @return An internal value
#'
#' @noRd
CSVRead <- function(x){
  name <- basename(x)
  internalstrings <- c("RelativeBrightness", ".csv")
  name <- Luciernaga::NameCleanUp(name, removestrings=internalstrings)
  Data <- read.csv(x, check.names=FALSE)
  Data <- Data %>% mutate(sample = name) %>% relocate(sample, .before=Cluster)
}


