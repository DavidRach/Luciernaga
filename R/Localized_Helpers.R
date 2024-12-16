#' Scrambled Egg Unmixing of Single-Colors
#'
#' @param x A Gating Set object
#' @param sample.name The keyword containing the fcs file name
#' @param removestrings A list of values to remove from name
#' @param subset A gating hierarchy level to sort cells at, expression values retrieved
#' from these
#' @param multiplier A number to scale the OLS coefficients by
#' @param outpath The return folder for the .fcs files
#' @param returntype Whether to return "fcs" or "flowframe"
#' @param Verbose For troubleshooting name after removestrings
#' @param addon Additional addon to append to the new .fcs file name
#' @param ratiopopcutoff Desired cutoff for detector detection
#' @param NumberFluors Desired number additional fluors in matrix
#'
#' @importFrom flowCore keyword
#' @importFrom flowWorkspace gs_pop_get_data
#' @importFrom flowCore exprs
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom utils read.csv
#' @importFrom dplyr pull
#' @importFrom dplyr arrange
#' @importFrom stats lsfit
#' @importFrom flowCore write.FCS
#' @importFrom purrr map
#'
#' @return A new .fcs file with the new columns appended
#' @noRd
SC_Unmix <- function(x, sample.name, removestrings, subset,
                     multiplier, outpath, returntype, Verbose,
                     addon, ratiopopcutoff, NumberFluors=1){

  # Retrieving Single Color Sample data for Unmixing
  if (length(sample.name) == 2){
    first <- sample.name[[1]]
    second <- sample.name[[2]]
    first <- keyword(x, first)
    second <- keyword(x, second)
    name <- paste(first, second, sep="_")
  } else {name <- keyword(x, sample.name)}

  name <- NameCleanUp(name, removestrings=removestrings)
  if (Verbose == TRUE){message("After removestrings, name is ", name)}

  cs <- gs_pop_get_data(x, subset)
  Data <- exprs(cs[[1]])
  Data <- data.frame(Data, check.names = FALSE)

  OriginalColumnsVector <- colnames(Data)
  OriginalColumns <- colnames(Data)
  OriginalColumns <- data.frame(OriginalColumns, check.names = FALSE)
  OriginalColumnsIndex <- OriginalColumns %>% mutate(IndexLocation = 1:nrow(.))

  Backups <- Data %>% mutate(Backups = 1:nrow(Data)) %>% select(Backups)

  StashedDF <- Data[,grep("Time|FS|SC|SS|Original|W$|H$", names(Data))]
  StashedDF <- cbind(Backups, StashedDF)

  TheSampleData <- Data[,-grep("Time|FS|SC|SS|Original|W$|H$", names(Data))]
  BackupNames <- colnames(TheSampleData)
  startingcells <- nrow(TheSampleData)
  TheColsN <- ncol(TheSampleData)

  # Retrieving Single Color Average Brightness
  PeakReturn <- Luciernaga:::PeakDetectors(NewData=TheSampleData)

  WorkAround <- PeakReturn[[1]]
  PeakDetectorCounts <- PeakReturn[[2]]
  CellCutoff <- startingcells*ratiopopcutoff
  Detectors <- PeakDetectorCounts %>% dplyr::filter(Counts > CellCutoff)
  Retained <- Luciernaga:::RetainTheDetectors(AFOverlap=AFOverlap, Detectors=Detectors, name=name)

  if (!length(Retained)>1){
  RetainedA <- paste0(Retained, "-A")
  MaxValues <- WorkAround %>% filter(.data[[Retained]] == 1) %>%
    select(all_of(RetainedA)) %>% pull(.)
  } else {
    MaxValues <- map(.x=Retained, .f=RetainedParse, data=WorkAround) %>% unlist()
  }
  #hist(MaxValues)
  TargetValue <- quantile(MaxValues, 0.95)
  TargetValue <- TargetValue[[1]]

  # Returning SC Data from References

  controlData <- Luciernaga:::ReferenceScramble(name=name, NumberDetectors=TheColsN,
                                                NumberFluors=NumberFluors)
  NewNames <- c("Fluorophore", colnames(TheSampleData))
  colnames(controlData) <- NewNames

  TheNames <- controlData %>% select(!where(is.numeric))
  Numerics <- controlData %>% select(where(is.numeric))
  Numerics <- Numerics*TargetValue
  MFIedControlData <- cbind(TheNames, Numerics)

  NewNames <- MFIedControlData %>% select(Fluorophore) %>% pull()
  TheControlData <- MFIedControlData %>% select(!Fluorophore)

  LeastSquares <- lsfit(x = t(TheControlData), y = t(TheSampleData), intercept = FALSE)
  UnmixedData <- t(LeastSquares$coefficients)
  UnmixedData2 <- UnmixedData*multiplier

  colnames(UnmixedData2) <- NewNames
  TheData <- cbind(StashedDF, UnmixedData2)
  TheData <- TheData %>% select(-Backups)
  rownames(TheData) <- NULL

  Ligands <- NewNames

  new_fcs <- Luciernaga:::InternalUnmix(cs=cs, StashedDF=StashedDF, TheData=TheData,
                                        Ligands=Ligands)
  #View(new_fcs@description)

  if (!is.null(addon)){name <- paste0(name, addon)}

  CurrentTime <- format(Sys.time(), "%H:%M:%S")
  CurrentTime <- gsub(":", "", CurrentTime)

  AssembledName <- paste0(name, "_", CurrentTime, ".fcs")

  new_fcs@description$GUID <- AssembledName
  new_fcs@description$`$FIL` <- AssembledName

  if (is.null(outpath)) {outpath <- getwd()}

  fileSpot <- file.path(outpath, AssembledName)

  if (returntype == "fcs") {write.FCS(new_fcs, filename = fileSpot, delimiter="#")
  } else {return(new_fcs)}
}

#' Internal for SC_Unmix
#'
#' @param NewData The exprs for the .fcs file minus scatter params
#'
#' @importFrom dplyr arrange
#' @importFrom dplyr desc
#'
#' @return A data.frame of detectors and respective counts
#' @noRd
PeakDetectors <- function(NewData){
  NewData[NewData < 0] <- 0
  A <- do.call(pmax, NewData)
  Normalized <- NewData/A
  colnames(Normalized) <- gsub("-A", "", colnames(Normalized))

  ColsN <- ncol(NewData)
  ColsNormalized <- ncol(Normalized)
  StartNormalizedMergedCol <- ColsN + 1
  EndNormalizedMergedCol <- ColsN + ColsNormalized

  WorkAround <- cbind(NewData, Normalized)

  na_counts <- colSums(is.na(Normalized))
  Normalized[is.na(Normalized)] <- 0
  Counts <- colSums(Normalized == 1)
  PeakDetectorCounts <- data.frame(Fluors = names(Counts), Counts = Counts)
  rownames(PeakDetectorCounts) <- NULL
  PeakDetectorCounts <- PeakDetectorCounts %>% arrange(desc(Counts))
  PeakList <- list(WorkAround, PeakDetectorCounts)
  return(PeakList)
}

#' Internal for SC_Unmix
#'
#' @param AFOverlap The Overlap list to remove AF Detecotrs
#' @param Detectors The passed PeakDetectorList to decide on
#' @param name The passed name consisting ligand fluorophore
#'
#' @importFrom utils read.csv
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#' @importFrom stringr str_split
#' @importFrom stringr str_detect
#'
#' @return The retained detector(s)
#' @noRd
RetainTheDetectors <- function(AFOverlap, Detectors, name){
  if (is.data.frame(AFOverlap)){AFData <- AFOverlap
  } else {AFData <- read.csv(file=AFOverlap, check.names = FALSE)
  }

  AFChannels <- AFData %>% filter(Fluorophore %in% "Unstained") %>%
    pull(MainDetector) %>% str_split(",", simplify = TRUE)
  AFChannels <- AFChannels[1,]
  AFChannels <- gsub("-A", "", AFChannels)

  TheSCData <- AFData %>% dplyr::filter(Fluorophore != "Unstained")
  TheSCData$Fluorophore <- gsub("-A", "", TheSCData$Fluorophore)
  TroubleChannels <- TheSCData %>% pull(Fluorophore)

  results <- map(.x=TroubleChannels, .f=TroubleChannelExclusion,
                 TheSCData=TheSCData, MainDetector=MainDetector,
                 AFChannels=AFChannels) %>% set_names(TroubleChannels)

  OverlapFlag <- NULL
  matching_names <- names(results)[str_detect(name, names(results))]
  if (length(matching_names) > 0) {
    OverlapFlag <- "Yep"
    ExclusionList <- results[[matching_names[1]]]
    Retained <- Detectors %>% dplyr::filter(!Fluors %in% ExclusionList) %>% pull(Fluors)
  } else if (str_detect(name, "nstained")){Retained <- Detectors %>%
    pull(Fluors)
  } else {Retained <- Detectors %>% dplyr::filter(!Fluors %in% AFChannels) %>%
    pull(Fluors)}

  # else {Retained <- Detectors %>% pull(Fluors)}
  if (length(Retained) == 0) {stop("There were no Retained detectors in ", name)}
  return(Retained)
}

#' Internal for SC Unmix
#'
#' @param x Passed Retained Detector
#' @param data The Raw and Normed Data.frame
#'
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom dplyr pull
#' @importFrom tidyselect all_of
#'
#' @return A vector of values to be passed to quantile
#' @noRd
RetainedParse <- function(x, data){
  RetainedA <- paste0(x, "-A")
  MaxValues <- data %>% filter(.data[[x]] == 1) %>%
    select(all_of(RetainedA)) %>% pull(.)
}

#' Internal for SC_Unmix
#'
#' @param name The passed ligand fluorophore name
#' @param NumberDetectors The number of detectors for referencing
#' @param NumberFluors Desired number of additional fluorophores in the matrix
#'
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#' @importFrom dplyr group_by
#' @importFrom dplyr ungroup
#' @importFrom dplyr arrange
#' @importFrom dplyr slice
#' @importFrom dplyr slice_sample
#' @importFrom dplyr desc
#' @importFrom tidyr pivot_wider
#'
#' @return A data.frame of the normalized signatures for the randomly selected fluorophores
#'
#' @noRd
ReferenceScramble <- function(name, NumberDetectors, NumberFluors){
  Spaces <- str_count(name, " ")
  if (Spaces <= 1){Name <- strsplit(name, " ")[[1]]
  } else {Name <- c(sub(" .*", "", name), sub("^[^ ]+ ", "", name))}
  Name <- Name[2]
  TheFluor <- QC_ReferenceLibrary(Name, NumberDetectors=NumberDetectors)
  TheFluor <- TheFluor[[1]]
  Similar <- QC_SimilarFluorophores(TheFluorophore=TheFluor,
                                    NumberDetectors=NumberDetectors, NumberHits=20)
  TooSimilar <- Similar %>% filter(.data[[TheFluor]] > 0.97) %>% pull(Fluorophore)
  TooSimilar <- c(TheFluor, TooSimilar)

  ReferenceData <- InstrumentReferences(NumberDetectors=NumberDetectors)

  ThePrelimGroups <- ReferenceData %>% group_by(Fluorophore) %>%
    arrange(desc(AdjustedY)) %>% slice(1) %>% select(Fluorophore, Detector) %>%
    ungroup() %>% arrange(Detector)

  TheDetectorGroups <- ThePrelimGroups %>% filter(!Fluorophore %in% TooSimilar)

  TheDetectorList <- TheDetectorGroups %>% select(Detector) %>% unique() %>% pull()
  TheRandomDetectors <- sample(TheDetectorList, NumberFluors, replace=FALSE)
  TheDetectorGroups_subset <- TheDetectorGroups %>% filter(Detector %in% TheRandomDetectors)
  TheSampling <- TheDetectorGroups_subset %>% group_by(Detector) %>%
    slice_sample(n=1) %>% ungroup()
  #TheSampling

  TheseFluors <- TheSampling %>% pull(Fluorophore)
  TheseFluors <- c(TheFluor, TheseFluors)

  ThisData <- ReferenceData %>% filter(Fluorophore %in% TheseFluors)
  ThisData <- ThisData %>% select(-Instrument) %>%
    pivot_wider(., names_from="Detector", values_from="AdjustedY")
  return(ThisData)
}

#' Generates a GatingSet with a positive negative gate on the SC Fluor.
#'
#' @param x A flowframe object
#'
#' @importFrom flowCore exprs
#' @importFrom dplyr summarise
#' @importFrom dplyr across
#' @importFrom tidyselect where
#' @importFrom stats quantile
#' @importFrom flowWorkspace cytoset
#' @importFrom flowWorkspace flowFrame_to_cytoframe
#' @importFrom flowCore keyword
#' @importFrom flowWorkspace cs_add_cytoframe
#' @importFrom flowWorkspace GatingSet
#' @importFrom flowWorkspace flowjo_biexp_trans
#' @importFrom flowWorkspace transformerList
#' @importFrom flowWorkspace transform
#' @importFrom data.table fread
#' @importFrom openCyto gatingTemplate
#' @importFrom openCyto gt_gating
#'
#' @return A Gating Set object transformed and gated
#' @noRd
ToSmallGatingSet <- function(x){
  data <- exprs(x)
  data <- data.frame(data, check.names=FALSE)
  data <- data[,-grep("Time|FS|SC|SS|Original|W$|H$", names(data))]
  TheMedian <- data %>%
    summarise(across(where(is.numeric), \(x) quantile(x, probs = 0.95, na.rm = TRUE)))
  KeptMarkers <- colnames(TheMedian)
  Fluorophore <- names(TheMedian)[which.max(TheMedian)]

  MyCS <- cytoset()
  cf <- flowFrame_to_cytoframe(x)

  TheName <- keyword(cf, "$FIL")
  TheName <- unlist(TheName)
  cs_add_cytoframe(MyCS, TheName, cf)
  MyGS <- GatingSet(MyCS)

  MyBiexponentialTransform <- flowjo_biexp_trans(channelRange = 256,
                                                 maxValue = 4194303,
                                                 pos = 5.62, neg = 0,
                                                 widthBasis = -1000)

  TransformList <- transformerList(KeptMarkers, MyBiexponentialTransform)
  UnmixedGatingSet <- transform(MyGS, TransformList)
  # pData(UnmixedGatingSet)

  FileLocation <- system.file("extdata", package = "Luciernaga")
  UnmixedGates <- fread(file.path(path = FileLocation,
                                  pattern = 'GatesUnmixed.csv'))
  Example <- UnmixedGates[6]
  Example[1,1] <- "Positive"
  Example[1,2] <- "+"
  Example[1,3] <- "root"
  Example[1,4] <- Fluorophore
  Example <- rbind(Example, Example)
  Example[1,1] <- "Negative"
  Example[1,2] <- "-"

  UnmixedGating <- gatingTemplate(Example)
  gt_gating(UnmixedGating, UnmixedGatingSet)
  #plot(UnmixedGatingSet)
  #autoplot(UnmixedGatingSet, "Positive", bins=100)
  #autoplot(UnmixedGatingSet, "Negative", bins=100)
  return(UnmixedGatingSet)
}

#' Generates Approximate Staining Index and Associated Info
#'
#' @param x The GatingSet object with applied transformation and gates.
#' @param NumberDetectors The number detectors the original file was based on.
#'
#' @importFrom flowWorkspace gs_pop_get_data
#' @importFrom flowCore exprs
#' @importFrom dplyr summarise
#' @importFrom dplyr across
#' @importFrom tidyselect where
#' @importFrom stats quantile
#' @importFrom dplyr select
#' @importFrom tidyselect all_of
#' @importFrom dplyr pull
#' @importFrom dplyr filter
#' @importFrom tidyr pivot_wider
#'
#' @return A data.frame row containing the derrived information
#' @noRd
StainingIndexApproximation <- function(x, NumberDetectors){
  name <- keyword(x, "$FIL")

  Positive <- gs_pop_get_data(x, "Positive", inverse.transform = FALSE)
  Positive <- exprs(Positive[[1]])
  Positive <- data.frame(Positive, check.names=FALSE)
  Positive <- Positive[,-grep("Time|FS|SC|SS|Original|W$|H$", names(Positive))]
  TheMedian <- Positive %>%
    summarise(across(where(is.numeric), \(x) quantile(x, probs = 0.95, na.rm = TRUE)))
  KeptMarkers <- colnames(TheMedian)
  Fluorophore <- names(TheMedian)[which.max(TheMedian)]
  Positive <- Positive %>% select(all_of(Fluorophore))

  Negative <- gs_pop_get_data(x, "Negative", inverse.transform = FALSE)
  Negative <- exprs(Negative[[1]])
  Negative <- data.frame(Negative, check.names=FALSE)
  Negative <- Negative[,-grep("Time|FS|SC|SS|Original|W$|H$", names(Negative))]
  Negative <- Negative %>% select(all_of(Fluorophore))

  PositiveVals <- Positive %>% pull(1)
  MFI_Pos <- quantile(PositiveVals, probs=0.5, na.rm =TRUE)
  MFI_Pos <- MFI_Pos[[1]]
  NegativeVals <- Negative %>% pull(1)
  MFI_Neg <- quantile(NegativeVals, probs=0.5, na.rm =TRUE)
  MFI_Neg <- MFI_Neg[[1]]
  #hist(NegativeVals)
  MinNeg <- quantile(NegativeVals, probs=0.03, na.rm =TRUE)
  MinNeg <- MinNeg[[1]]
  MaxNeg <- quantile(NegativeVals, probs=0.97, na.rm =TRUE)
  MaxNeg <- MaxNeg[[1]]
  RSD <- (MaxNeg - MinNeg)/3.29
  StainingIndex <- (MFI_Pos - MFI_Neg)/(2*RSD)

  ReferenceData <- Luciernaga:::InstrumentReferences(NumberDetectors=NumberDetectors)
  Matrix <- ReferenceData %>% filter(Fluorophore %in% KeptMarkers) %>%
    select(-Instrument) %>% pivot_wider(., names_from="Detector", values_from="AdjustedY") %>%
    select(-Fluorophore) %>% as.matrix()
  kappa <- round(kappa(Matrix), 2)

  MatrixN <- length(KeptMarkers)
  MatrixMarkers <- paste(KeptMarkers, collapse = "_")
  FinalData <- cbind(name, Fluorophore,  MatrixN, MatrixMarkers, StainingIndex, kappa)
  colnames(FinalData)[1] <- "name"

  return(FinalData)
}

#' Wrapper Function for Scrambled Eggs Protocol
#'
#' @param NumberRepeats Desired number of bootstrap runs
#' @param NumberFluors Desired Number of Additional Fluorophores
#' @param NumberDetectors Detector configuration original FCS file
#' @param GS The GatingSet object corresponding to desired SC
#' @param sample.name The keyword designating single-color sample name
#' @param removestrings Values to be removed to leave just the sample name
#' @param subset The desired gating node
#' @param multiplier The multiplier for unmixing, default set to 50000
#' @param outpath Internal outpath argument
#' @param addon Internal unmix argument
#'
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
#'
#' @return A data.frame summarizing the Staining Index and Kappa for each bootstrap
#' @noRd
ScrambledEggs <- function(NumberRepeats, NumberFluors, NumberDetectors, GS, sample.name="TUBENAME",
                          removestrings, subset="lymphocytes", multiplier=50000,
                          outpath, addon=addon){
  Data <- list()
  for(i in seq_along(1:NumberRepeats)){
    Data[[i]] <- Luciernaga:::SC_Unmix(x=GS, sample.name=sample.name,
                                       removestrings=removestrings, subset=subset, multiplier=multiplier,
                                       outpath=outpath, returntype="data", Verbose=FALSE, addon="_Unmixed",
                                       ratiopopcutoff=0.01, NumberFluors=NumberFluors)
  }
  GatedData <- map(.x=Data, .f=Luciernaga:::ToSmallGatingSet)
  MyData <- map(.x=GatedData, .f=Luciernaga:::StainingIndexApproximation, NumberDetectors=NumberDetectors) %>% bind_rows()
  return(MyData)
}



