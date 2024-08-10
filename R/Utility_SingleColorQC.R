#' Main Luciernaga Function, normalized on individual cell level
#'
#' @param x A Gating Set object (ex. gs or gs[[1]])
#' @param subsets The Gating Hierarchy level you will be sampling at
#' @param sample.name Keyword variable which samples are stored (ex. "GUID")
#' @param group.name Keyword variable which groups are stored (ex. "GROUPNAME")
#' @param experiment Provide directly experiment name (ex. "JAN2024")
#' @param experiment.name Keyword variable which experiment information
#' is stored (ex. "TUBENAME")
#' @param stats Whether to take "mean" or "median"
#' @param Kept Whether "Raw" or "Normalized" values are retained in the
#' Luciernaga object.
#' @param external An external autofluorescence to subtract from single colors.
#' @param sourcelocation Location where .fcs creation file is stored
#' @param outpath  Location where created .fcs and .csv files are sent
#' @param artificial Whether an artificial 0 population should be added for
#'  a background autofluorescence stand in.
#' @param fcsexport Whether to export .fcs files, TRUE or FALSE
#' @param mainAF Main Autofluorescence Detector (ex. "V7-A")
#' @param AFOverlap Name of data.frame containing the Autofluorescence
#' overlap of individual fluorophores for exclusion
#' @param Beads  Whether the sample is Beads.
#' @param Brightness Whether sum of detectors should be returned.
#' @param Unstained Whether the sample is Unstained.
#'
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom dplyr summarize_all
#' @importFrom dplyr pull
#' @importFrom dplyr arrange
#' @importFrom dplyr relocate
#' @importFrom dplyr left_join
#' @importFrom dplyr case_when
#' @importFrom dplyr rename
#' @importFrom BiocGenerics nrow
#' @importFrom flowWorkspace keyword
#' @importFrom stringr str_detect
#' @importFrom stringr str_split
#' @importFrom flowWorkspace gs_pop_get_data
#' @importFrom flowCore exprs
#' @importFrom purrr map
#' @importFrom purrr set_names
#'
#'
#' @return Additional information to be added
#' @export
#'
#' @examples NULL

Utility_SingleColorQC <- function(x, subsets, sample.name, removestrings, unmixingcontroltype,
                                  experiment = NULL, experiment.name = NULL,
                                  mainAF, AFOverlap, stats, Unstained=FALSE, Beads=FALSE, Verbose = FALSE,
                                  external = NULL, fcsexport, sourcelocation, outpath, artificial, Brightness=FALSE, ...){

  name <- keyword(x, sample.name)
  Type <- Luciernaga:::Internal_Typing(name=name, unmixingcontroltype=unmixingcontroltype, Unstained=Unstained)
  #Type <- Internal_Typing(name=name, unmixingcontroltype=unmixingcontroltype, Unstained=Unstained)
  AlternateName <- Luciernaga:::NameForSample(x=x, sample.name=sample.name, removestrings=removestrings)
  #AlternateName <- NameForSample(x=x, sample.name=sample.name, removestrings=removestrings, ...)

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


  ###############
  # Exprs Setup #
  ###############

  #Retrieving the exprs data for the target population
  ff <- gs_pop_get_data(x, subsets)
  startingcells <- nrow(ff)[[1]]
  df <- exprs(ff[[1]])
  DF <- as.data.frame(df, check.names = FALSE)

  #For Future Column Reordering
  OriginalColumnsVector <- colnames(DF)
  OriginalColumns <- colnames(DF)
  OriginalColumns <- data.frame(OriginalColumns)
  OriginalColumnsIndex <- OriginalColumns %>% mutate(IndexLocation = 1:nrow(.))

  #For Future Row Reordering
  Backups <- DF %>% mutate(Backups = 1:nrow(DF)) %>% select(Backups)

  #Stashing Away Time FSC SSC For Later Use
  StashedDF <- DF[,grep("Time|FS|SC|SS|Original|W$|H$", names(DF))]
  StashedDF <- cbind(Backups, StashedDF)

  #Consolidating Columns Going Forward
  CleanedDF <- DF[,-grep("Time|FS|SC|SS|Original|W$|H$", names(DF))]
  BackupNames <- colnames(CleanedDF)
  n <- CleanedDF

  #################
  # Peak Detector #
  #################

  # Enumerating Negative Values Ratio
  if (Verbose == TRUE){
  TotalCells <- nrow(n) * ncol(n)
  BelowZero <- sum(apply(n, 2, function(x) x < 0))
  message(round(BelowZero/TotalCells,2), " of all events were negative and rounded to 0")
  }

  #Normalizing By Peak Detector
  n[n < 0] <- 0
  A <- do.call(pmax, n)
  Normalized <- n/A
  Normalized <- round(Normalized, 1)
  colnames(Normalized) <- gsub("-A", "", colnames(Normalized))

  #Towards Generalizing
  ColsN <- ncol(n)
  ColsNormalized <- ncol(Normalized)
  StartNormalizedMergedCol <- ColsN + 1
  EndNormalizedMergedCol <- ColsNormalized*2

  #Bringing Together Raw And Normalized and deriving Autofluorescence Negative
  WorkAround <- cbind(n, Normalized)

  ####################
  # Autofluorescence #
  ####################

  #Retrieving Main Auto fluorescent Channels signature
  if (Unstained == FALSE){MainAF <- mainAF
                          MainAF <- gsub("-A", "", MainAF)
  } else if (!grepl("nstained", name)){MainAF <- mainAF
                                      MainAF <- gsub("-A", "", MainAF)
  } else {na_counts <- colSums(is.na(Normalized))
    Normalized[is.na(Normalized)] <- 0
    Counts <- colSums(Normalized == 1)
    PeakDetectorCounts <- data.frame(Fluors = names(Counts), Counts = Counts)
    rownames(PeakDetectorCounts) <- NULL
    MainAF <- PeakDetectorCounts %>% arrange(desc(Counts)) %>% slice(1) %>% pull(Fluors)
  }

  # Grabbing Main AF Peak and the Associated Raw Values
  This <- WorkAround %>% filter(.data[[MainAF]] == 1) %>%
    select(all_of(1:ColsN))

  # Deriving Middle Autofluorescence Measurement
  if(stats == "mean"){Samples <- This %>% summarize_all(mean)
  } else if (stats == "median"){Samples <- This %>%
    summarize_all(median) #%>% select(-Backups)
  } else(stop("Please specify stats parameter mean or median"))

  ################
  # Fluorophores #
  ################

  #Deriving Peak Detector Counts and Detectors of Interest
  InitialRatio <- 0.0075  #Possible Parameter Add Here

  na_counts <- colSums(is.na(Normalized))
  Normalized[is.na(Normalized)] <- 0
  Counts <- colSums(Normalized == 1)
  PeakDetectorCounts <- data.frame(Fluors = names(Counts), Counts = Counts)
  rownames(PeakDetectorCounts) <- NULL
  cutoff <- startingcells*InitialRatio
  Detectors <- PeakDetectorCounts %>% filter(Counts > cutoff) %>%
    arrange(desc(Counts))

  ################################################################
  # Bringing in known AF detectors, and overlapping fluorophores #
  ################################################################

  AFData <- AFOverlap
  AFChannels <- AFData %>% filter(Fluorophore %in% "Unstained") %>%
    pull(MainDetector) %>% str_split(",", simplify = TRUE)
  AFChannels <- AFChannels[1,]
  AFChannels <- gsub("-A", "", AFChannels)

  SCData <- AFData %>% filter(Fluorophore != "Unstained")
  SCData$Fluorophore <- gsub("-A", "", SCData$Fluorophore)
  TroubleChannels <- SCData %>% pull(Fluorophore)

  # Handling the Exceptions
  results <- map(.x=TroubleChannels, .f=TroubleChannelExclusion, SCData=SCData, MainDetector=MainDetector, AFChannels=AFChannels) %>%
    set_names(TroubleChannels)

  matching_names <- names(results)[str_detect(name, names(results))]
  if (length(matching_names) > 0) {ExclusionList <- results[[matching_names[1]]]
  Retained <- Detectors %>% filter(!Fluors %in% ExclusionList) %>% pull(Fluors)
  } else if (str_detect(name, "nstained")){Retained <- Detectors %>%
    pull(Fluors)
  } else {Retained <- Detectors %>% filter(!Fluors %in% AFChannels) %>%
    pull(Fluors)}

  if (length(Retained) == 0) {
    stop("There were no Retained detectors in ", name)
  }

  if (!str_detect(name, "nstained")){
    Intermediate <- Detectors %>% filter(!Fluors %in% Retained)
    if (nrow(Intermediate) >0){
      TheMainAF <- Intermediate %>% slice(1) %>% pull(Fluors)
    } else {TheMainAF <- MainAF}
  }

  if(Beads == TRUE){Retained <- Retained[[1]]}

  #####################################
  # We resume our regular programming #
  #####################################

  WorkAround1 <- WorkAround %>% mutate(Backups = Backups$Backups) %>%
    relocate(Backups, .before = 1) #This will change the start/end count

  if (str_detect(name, "nstained")){

    RetainedDF <- map(.x= Retained, .f=UnstainedSignatures,
                      WorkAround1=WorkAround1, alternatename=AggregateName,
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

  if (fcsexport == TRUE){source(sourcelocation, local = TRUE)}

  return(Reintegrated1)
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




