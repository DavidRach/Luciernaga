#' Main Luciernaga Function, normalized on individual cell level
#'
#' @param x A Gating Set object (ex. gs or gs[[1]])
#' @param subsets The Gating Hierarchy level you will be sampling at
#' @param sample.name Keyword variable which samples are stored (ex. "GUID")
#' @param removestrings A string of character values to remove from sample.name
#' @param Verbose Whether to print outputs as you go.
#' @param unmixingcontroltype Whether your inputs are "cells", "beads" or "both"
#' @param Unmixing Set to True when running unstained samples that don't have Unstained in the Name
#' @param ratiopopcutoff A numeric ratio for peak detector inclusion, default is set to 0.01 all startingcells
#' @param experiment Provide directly experiment name (ex. "JAN2024")
#' @param experiment.name Keyword variable which experiment information
#' is stored (ex. "TUBENAME")
#' @param condition Provide directly experiment name (ex. "JAN2024")
#' @param condition.name Keyword variable which experiment information
#' is stored (ex. "TUBENAME")
#' @param AFOverlap A data.frame or a filepath to the CSV containing the Autofluorescence
#' overlap of individual fluorophores for exclusion
#' @param stats Whether to take "mean" or "median"
#' @param desiredAF Main Autofluorescence Detector (ex. "V7-A")
#' @param externalAF A data.frame row containing external autofluorescence to subtract from single colors.
#' @param ExportType Whether to return "fcs", "data.frame" or "csv"
#' @param SignatureReturnNow Short circuits the function and returns signature for specified autofluorescence.
#' @param sourcelocation Location where .fcs creation file is stored
#' @param outpath  Location where created .fcs and .csv files are sent
#' @param artificial Whether an artificial 0 population should be added for
#'  a background autofluorescence stand in.
#' @param Brightness Whether sum of detectors should be returned.
#'
#' @importFrom flowWorkspace gs_pop_get_data
#' @importFrom BiocGenerics nrow
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom dplyr summarize_all
#' @importFrom dplyr pull
#' @importFrom dplyr arrange
#' @importFrom dplyr relocate
#' @importFrom dplyr left_join
#' @importFrom dplyr case_when
#' @importFrom dplyr rename
#' @importFrom flowWorkspace keyword
#' @importFrom stringr str_detect
#' @importFrom stringr str_split
#' @importFrom flowCore exprs
#' @importFrom purrr map
#' @importFrom purrr set_names
#'
#' @return Additional information to be added
#' @export
#'
#' @examples NULL

Utility_SingleColorQC <- function(x, subsets, sample.name, removestrings, Verbose = FALSE,
                                  unmixingcontroltype="both", Unstained=FALSE, ratiopopcutoff=0.01,
                                  experiment = NULL, experiment.name = NULL, condition = NULL, condition.name = NULL,
                                  AFOverlap, stats="median", desiredAF=NULL, externalAF = NULL,
                                  ExportType, SignatureReturnNow=FALSE, sourcelocation, outpath, artificial, Brightness=FALSE){

  name <- keyword(x, sample.name)
  Type <- Luciernaga:::Internal_Typing(name=name, unmixingcontroltype=unmixingcontroltype, Unstained=Unstained)
  #Type <- Internal_Typing(name=name, unmixingcontroltype=unmixingcontroltype, Unstained=Unstained)
  AggregateName <- Luciernaga:::NameForSample(x=x, sample.name=sample.name, removestrings=removestrings)
  #AggregateName <- NameForSample(x=x, sample.name=sample.name, removestrings=removestrings, ...)

  Experiment <- Luciernaga:::NameForSample(x=x, sample.name=sample.name, removestrings=removestrings, experiment.name = experiment.name,
                                           returnType = "experiment")
  #Experiment <- Luciernaga:::NameForSample(x=x, sample.name=sample.name, removestrings=removestrings, returnType = "experiment", ...)

  Condition <- Luciernaga:::NameForSample(x=x, sample.name=sample.name, removestrings=removestrings, condition.name = condition.name,
                                           returnType = "condition")
  #Condition <- Luciernaga:::NameForSample(x=x, sample.name=sample.name, removestrings=removestrings, returnType = "condition", ...)

  # Internal Name Cleanup, to make sure fluorophores match
  InternalCleanupList <- c(".fcs", "Cells", "Beads", " ", "_", "-", ".", "(", ")")
  name <- Luciernaga:::NameCleanUp(name, InternalCleanupList)
  #name <- NameCleanUp(name, InternalCleanupList)

  # Brought back until overlap is converted from name reliant to type reliant.
  if (Unstained == TRUE) {if(!str_detect(name, "stained")){name <- paste0(name, "_Unstained")}}

  ###############
  # Exprs Setup #
  ###############

  #Retrieving the exprs data for the target population
  ff <- gs_pop_get_data(x, subsets)
  startingcells <- nrow(ff)[[1]]
  DF <- as.data.frame(exprs(ff[[1]]), check.names=FALSE)

  #For Future Row Reordering
  Backups <- DF %>% mutate(Backups = 1:nrow(DF)) %>% select(Backups)

  #For Future Column Reordering
  OriginalColumns <- colnames(DF)
  OriginalColumns <- data.frame(OriginalColumns)
  OriginalColumnsIndex <- OriginalColumns %>% mutate(IndexLocation = 1:nrow(.))

  OriginalColumnsVector <- colnames(DF)

  # Storing Other Parameter (Time, FSC, SSC) for Later Return
  StashedDF <- DF[,grep("Time|FS|SC|SS|Original|W$|H$", names(DF))]
  StashedDF <- cbind(Backups, StashedDF)

  # Consolidating Desired Parameters
  n <- DF[,-grep("Time|FS|SC|SS|Original|W$|H$", names(DF))]

  # Next Step Will Remove Negative Values, to what extent is an ordinary .fcs file affected?
  if (Verbose == TRUE){
    TheTotal <- nrow(n) * ncol(n)
    BelowZero <- sum(apply(n, 2, function(x) x < 0))
    message(round(BelowZero/TheTotal,2), " of all events were negative and will be rounded to 0")
  }

  ##########################################
  # Generating Normalized by Peak Detector #
  ##########################################

  # Normalizing Individual Cells By Peak Detector
  n[n < 0] <- 0
  A <- do.call(pmax, n)
  Normalized <- n/A
  Normalized <- round(Normalized, 2) # Previously at 1
  colnames(Normalized) <- gsub("-A", "", colnames(Normalized))

  # Figuring out where Raw and Normalized Columns Start and End
  ColsN <- ncol(n)
  ColsNormalized <- ncol(Normalized)
  StartNormalizedMergedCol <- ColsN + 1
  EndNormalizedMergedCol <- ColsN + ColsNormalized

  # Bringing Together Raw And Normalized Data.Frames
  WorkAround <- cbind(n, Normalized)

  ####################################
  # Enumerating Peak Detector Counts #
  ####################################

  na_counts <- colSums(is.na(Normalized))
  Normalized[is.na(Normalized)] <- 0
  Counts <- colSums(Normalized == 1)
  PeakDetectorCounts <- data.frame(Fluors = names(Counts), Counts = Counts)
  rownames(PeakDetectorCounts) <- NULL
  PeakDetectorCounts <- PeakDetectorCounts %>% arrange(desc(Counts))

  if (Type == "Cells_Unstained") {CellCutoff <- startingcells*ratiopopcutoff
                                 Detectors <- PeakDetectorCounts %>% filter(Counts > CellCutoff)
  }

  if (Type == "Cells") {CellCutoff <- startingcells*ratiopopcutoff
                        Detectors <- PeakDetectorCounts %>% filter(Counts > CellCutoff)
  }

  if (Type == "Beads_Unstained") {BeadCutoff <- startingcells/ColsN #EqualDistributionAssumption
                                  Detectors <- PeakDetectorCounts %>% filter(Counts > BeadCutoff)
  }


  if (Type == "Beads") {BeadCutoff <- startingcells*ratiopopcutoff
                        Detectors <- PeakDetectorCounts %>% filter(Counts > BeadCutoff)
  }

  if (Verbose == TRUE){
  print(Detectors)
  }

  #####################################################
  # Handling Fluorophore and Autofluorescent Overlaps #
  #####################################################
  if (is.data.frame(AFOverlap)){AFData <- AFOverlap
  } else {AFData <- read.csv(file=AFOverlap, check.names = FALSE)
  }

  AFChannels <- AFData %>% filter(Fluorophore %in% "Unstained") %>%
    pull(MainDetector) %>% str_split(",", simplify = TRUE)
  AFChannels <- AFChannels[1,]
  AFChannels <- gsub("-A", "", AFChannels)

  SCData <- AFData %>% filter(Fluorophore != "Unstained")
  SCData$Fluorophore <- gsub("-A", "", SCData$Fluorophore)
  TroubleChannels <- SCData %>% pull(Fluorophore)

  results <- map(.x=TroubleChannels, .f=TroubleChannelExclusion, SCData=SCData, MainDetector=MainDetector, AFChannels=AFChannels) %>%
    set_names(TroubleChannels)

  matching_names <- names(results)[str_detect(name, names(results))]
  if (length(matching_names) > 0) {ExclusionList <- results[[matching_names[1]]]
  Retained <- Detectors %>% filter(!Fluors %in% ExclusionList) %>% pull(Fluors)
  } else if (str_detect(name, "nstained")){Retained <- Detectors %>%
    pull(Fluors)
  } else {Retained <- Detectors %>% filter(!Fluors %in% AFChannels) %>%
    pull(Fluors)}

  if (length(Retained) == 0) {stop("There were no Retained detectors in ", name)}

  if (Verbose == TRUE){print(Retained)}


  ##################################
  # Main Autofluorescence Detector #
  ##################################

  # Inverse of retained, then top for peak

  if (!str_detect(name, "nstained")){
    Intermediate <- Detectors %>% filter(!Fluors %in% Retained)
    if (nrow(Intermediate) >0){
      TheMainAF <- Intermediate %>% slice(1) %>% pull(Fluors)
    }

  # Adding way to switch to alternate AFs.
  if (is.null(desiredAF)){This <- WorkAround %>% filter(.data[[TheMainAF]] == 1) %>% select(all_of(1:ColsN))
  } else {desiredAF <- NameCleanUp(desiredAF, c("-A"))
          TheMainAF <- desiredAF
          This <- WorkAround %>% filter(.data[[desiredAF]] == 1) %>% select(all_of(1:ColsN))}

  if (is.null(externalAF)){Samples <- AveragedSignature(x=This, stats=stats)
  } else {if(is.data.frame(externalAF)){

    #Filter and designate the new MainAF from the example (without requiring matching name)
    Samples <- externalAF} else {stop("externalAF needs to be a single row of a data.frame")}}

  }

  if (str_detect(name, "nstained")){
    Intermediate <- Detectors
    if (nrow(Intermediate) >0){
      TheMainAF <- Intermediate %>% slice(1) %>% pull(Fluors)
    }

    # Adding way to switch to alternate AFs.
    if (is.null(desiredAF)){This <- WorkAround %>% filter(.data[[TheMainAF]] == 1) %>% select(all_of(1:ColsN))
    } else {desiredAF <- NameCleanUp(desiredAF, c("-A"))
    TheMainAF <- desiredAF
    This <- WorkAround %>% filter(.data[[desiredAF]] == 1) %>% select(all_of(1:ColsN))}

    if (is.null(externalAF)){Samples <- AveragedSignature(x=This, stats=stats)
    } else {if(is.data.frame(externalAF)){

      #Filter and designate the new MainAF from the example (without requiring matching name)
      Samples <- externalAF} else {stop("externalAF needs to be a single row of a data.frame")}}

  }

  if (SignatureReturnNow == TRUE){return(Samples)}

  ######################################
  # Handling the Double Peak Detectors #
  ######################################

  WorkAround1 <- WorkAround %>% mutate(Backups = Backups$Backups) %>%
    relocate(Backups, .before = 1) #This will change the start/end count

  StartN <- StartNormalizedMergedCol+1
  EndN <- EndNormalizedMergedCol+1
  DoublePeaks <- WorkAround1 %>% filter(rowSums(WorkAround1[, StartN:EndN] == 1.00) >= 2)
  TheDoublePeaks <- DoublePeaks %>% select(Backups) %>% pull()
  WorkAround2 <- WorkAround1 %>% filter(!Backups %in% TheDoublePeaks)

  ###############################
  # Back to Regular Programming #
  ###############################

  if (str_detect(name, "nstained")){
    # x <- Retained[1]
    RetainedDF <- map(.x= Retained, .f=Luciernaga:::UnstainedSignatures,
                      WorkAround1=WorkAround2, alternatename=AggregateName,
                      ColsN=ColsN, StartNormalizedMergedCol=StartNormalizedMergedCol,
                      EndNormalizedMergedCol=EndNormalizedMergedCol) %>% bind_rows()

  } else {

    RetainedDF <- map(.x= Retained, .f=SingleStainSignatures,
                      WorkAround1=WorkAround1, alternatename=AggregateName,
                      ColsN=ColsN, StartNormalizedMergedCol=StartNormalizedMergedCol,
                      EndNormalizedMergedCol=EndNormalizedMergedCol, Samples=Samples, name=name,
                      results=results, stats=stats, TheMainAF=TheMainAF,
                      AggregateName=AggregateName) %>% bind_rows()
  }

  Reintegrated <- left_join(RetainedDF, StashedDF, by = "Backups")

  BackupsCol <- "Backups"
  NormalizedColumns <- colnames(Normalized)
  ClusterCol <- "Cluster"
  RearrangedColumns <- c(BackupsCol, OriginalColumnsVector, NormalizedColumns,
                         ClusterCol)

  OriginalStart <- length(BackupsCol) + 1
  OriginalEnd <- length(BackupsCol) + length(OriginalColumnsVector)

  Reintegrated1 <- Reintegrated %>% relocate(all_of(RearrangedColumns))

  if (ExportType == "fcs"){source(sourcelocation, local = TRUE)
  }

  if (ExportType == "data.frame"){return(Reintegrated1)}
}



#' Internal for Utility_SingleColorQC
#'
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#' @noRd
TroubleChannelExclusion <- function(x, SCData, MainDetector, AFChannels){
  Internal <- SCData %>% filter(Fluorophore %in% x) %>% pull(MainDetector)
  Internal <- gsub("-A", "", Internal)
  Exclusion <- setdiff(AFChannels, Internal)
  return(Exclusion)
}

#' Internal for Utility_SingleColorQC
#'
#' @importFrom dplyr summarize_all
#' @noRd
AveragedSignature <- function(x, stats){
  Signature <- x %>% summarize_all(stats)
  return(Signature)
}

#' Internal for Utility_SingleColorQC
#'
#' @importFrom package function
#' @noRd
#'

DetectorPeakCounts <- function(x, StartN, EndN){
x <- x %>% select(-Backups)
Normalized <- x %>% select(all_of(
  StartN:EndN))
Counts <- colSums(Normalized == 1)
PeakDetectorCounts <- data.frame(Fluors = names(Counts), Counts = Counts)
rownames(PeakDetectorCounts) <- NULL
PeakDetectorCounts <- PeakDetectorCounts %>% arrange(desc(Counts))
return(PeakDetectorCounts)
}


