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
ReferenceScramble <- function(name, NumberDetectors, NumberFluors=5){
  Name <- strsplit(name, " ")[[1]]
  Name <- Name[2]
  TheFluor <- QC_ReferenceLibrary(Name, NumberDetectors=NumberDetectors)
  TheFluor <- TheFluor[[1]]
  Similar <- QC_SimilarFluorophores(TheFluorophore=TheFluor,
                                    NumberDetectors=NumberDetectors, NumberHits=20)
  TooSimilar <- Similar %>% filter(.data[[TheFluor]] > 0.97) %>% pull(Fluorophore)
  TooSimilar <- c(TheFluor, TooSimilar)

  ReferenceData <- InstrumentReferences(NumberDetectors=ColsN)

  ThePrelimGroups <- ReferenceData %>% group_by(Fluorophore) %>%
    arrange(desc(AdjustedY)) %>% slice(1) %>% select(Fluorophore, Detector) %>%
    ungroup() %>% arrange(Detector)

  TheDetectorGroups <- ThePrelimGroups %>% filter(!Fluorophore %in% TooSimilar)

  TheDetectorList <- TheDetectorGroups %>% select(Detector) %>% unique() %>% pull()
  TheRandomDetectors <- sample(TheDetectorList, NumberFluors, replace=FALSE)
  TheDetectorGroups_subset <- TheDetectorGroups %>% filter(Detector %in% TheRandomDetectors)
  TheSampling <- TheDetectorGroups_subset %>% group_by(Detector) %>% slice_sample(1) %>% ungroup()
  #TheSampling

  TheseFluors <- TheSampling %>% pull(Fluorophore)
  TheseFluors <- c(TheFluor, TheseFluors)

  ThisData <- ReferenceData %>% filter(Fluorophore %in% TheseFluors)
  ThisData %>% select(-Instrument) %>%
    pivot_wider(., names_from="Detector", values_from="AdjustedY")
  return(ThisData)
}



