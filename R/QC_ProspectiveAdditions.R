#' From existing panel, figures out open detectors, and returns potential
#'  fluorophores that might fit between the existing ones
#'
#' @param path Filepath to the panel .csv
#' @param NumberDetectors Number of detectors for the instrument to pull references
#' @param TheCutoff Default is 0.9, cosine matrix value
#' @param returnAll Whether to return all variants, default is FALSE
#' @param returnCSV Whether to return as a csv to designated outpath and filename, a TRUE/FALSE
#' @param filename Desired name for the output .csv
#' @param outpath Path Location to store the the output .csv
#'
#' @importFrom utils read.csv
#' @importFrom dplyr pull
#' @importFrom dplyr filter
#' @importFrom dplyr group_by
#' @importFrom dplyr arrange
#' @importFrom dplyr desc
#' @importFrom dplyr slice
#' @importFrom dplyr select
#' @importFrom dplyr ungroup
#' @importFrom dplyr relocate
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
#' @importFrom utils write.csv
#'
#' @return A csv containing selected fluorophores
#' @export
#'
#' @examples
#' Folder_Location <- system.file("extdata", package = "Luciernaga")
#' ThePanelLocation <- list.files(Folder_Location, pattern="^Panel.csv",
#'  full.names=TRUE)
#' OutPath <- file.path("C:", "Users", "JohnDoe", "Desktop")
#'
#' ProspectiveAdditions <- QC_ProspectiveAdditions(path=ThePanelLocation, NumberDetectors=64,
#' TheCutoff=0.9, returnAll=FALSE, returnCSV=FALSE,
#' filename="ProspectiveAdditions", outpath=OutPath)
#'
QC_ProspectiveAdditions <- function(path, NumberDetectors, TheCutoff=0.9,
                                    returnAll=FALSE, filename, outpath,
                                    returnCSV){

  TheList <- read.csv(path, check.names=FALSE)
  ReferenceData <- InstrumentReferences(NumberDetectors=NumberDetectors)

  TheList <- TheList %>% pull(Fluorophore)
  #x <- TheList[1]

  TheReferenceList <- ReferenceData %>% dplyr::filter(Fluorophore %in% TheList)

  TheDetectorContained <- TheReferenceList %>% group_by(Fluorophore) %>%
    arrange(desc(AdjustedY)) %>% slice(1) %>% select(Fluorophore, Detector) %>%
    ungroup()

  TheOccupiedDetectors <- TheDetectorContained %>% pull(Detector) %>% sort()

  TheOtherFluors <-  ReferenceData %>% dplyr::filter(!Fluorophore %in% TheList)
  TheOtherDetectors <- TheOtherFluors %>% group_by(Fluorophore) %>%
    arrange(desc(AdjustedY)) %>% slice(1) %>% select(Fluorophore, Detector) %>%
    ungroup()
  TheOtherDetectors <- TheOtherDetectors %>% arrange(Detector)
  TheOtherDetectors <- TheOtherDetectors %>%
    dplyr::filter(!Detector %in% TheOccupiedDetectors)

  TheOtherDetectors <- TheOtherDetectors %>% relocate(Detector, .before=Fluorophore)

  PossibleDetectors <- TheOtherDetectors %>% pull(Detector) %>%
    unique() #%>% length()

  DataUnlocked <- map(.x=PossibleDetectors, .f=Comparison,
                      TheOtherDetectors=TheOtherDetectors,
                      TheList=TheList, ReferenceData=ReferenceData,
                      TheCutoff=TheCutoff) %>% bind_rows()

  if (returnAll == TRUE){PossibleLocations <- DataUnlocked %>%
    group_by(TheDetector) %>% arrange(RankValue) %>% ungroup()
  } else {PossibleLocations <- DataUnlocked %>% group_by(TheDetector) %>%
    arrange(RankValue) %>% slice(1) %>% ungroup()
  }

  TheOutput <- PossibleLocations %>% arrange(RankValue)

  if (returnCSV == TRUE) {
    TheFileName <- paste0(filename, ".csv")
    StoreHere <- file.path(outpath, TheFileName)
    write.csv(TheOutput, StoreHere, row.names=FALSE)
  }

  return(TheOutput)
}


#' Internal for QC_ProspectiveAdditions
#'
#' @param x Passed Argument
#' @param TheOtherDetectors Passed Argument
#' @param TheList Passed Argument
#' @param ReferenceData Passed Argument
#' @param TheCutoff Passed Argument
#'
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
#'
#' @return An internal value
#'
#' @noRd
Comparison <- function(x, TheOtherDetectors, TheList, ReferenceData, TheCutoff){

  TheDetector <- x

  TheComparisonList <- TheOtherDetectors %>% dplyr::filter(Detector %in% x) %>%
    pull(Fluorophore)

  TheIndividualDetector <- map(.x=TheComparisonList, .f=InternalComparison,
                               TheList=TheList, ReferenceData=ReferenceData,
                               TheCutoff=TheCutoff, TheDetector=TheDetector) %>%
    bind_rows()

  return(TheIndividualDetector)
}




#' Internal for QC_ProspectiveAdditions
#'
#' @param x Passed Argument
#' @param TheList Passed Argument
#' @param ReferenceData Passed Argument
#' @param TheCutoff Passed Argument
#' @param TheDetector Passed Argument
#'
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom tidyr pivot_wider
#' @importFrom dplyr pull
#' @importFrom lsa cosine
#' @importFrom tidyselect all_of
#' @importFrom dplyr sym
#'
#' @return An internal value
#'
#' @noRd
InternalComparison <- function(x, TheList, ReferenceData, TheCutoff, TheDetector){

  TheCandidate <- ReferenceData %>% dplyr::filter(Fluorophore %in% x)
  TheReferenceList <- ReferenceData %>% dplyr::filter(Fluorophore %in% TheList)
  TheData <- rbind(TheCandidate, TheReferenceList) %>% select(-Instrument)
  TheCosineData <- TheData %>%
    pivot_wider(names_from = Detector, values_from = AdjustedY)
  Names <- TheCosineData %>% pull(Fluorophore)
  TheValues <- TheCosineData %>% select(-Fluorophore)
  TheTransposed <- t(TheValues)
  colnames(TheTransposed) <- Names
  TheMatrix <- data.matrix(TheTransposed)
  CosineMatrix <- cosine(TheMatrix)
  CosineMatrix <- data.frame(CosineMatrix, check.names = FALSE)
  TheCandidateValues <- CosineMatrix %>% select(all_of(x))
  HighOverlaps <- TheCandidateValues %>% filter(!!sym(x) >= TheCutoff)
  HighOverlaps <- length(HighOverlaps)
  RankValue <- round(kappa(TheMatrix),2)
  Fluorophore <- x

  Prelim <- cbind(Fluorophore, TheDetector, HighOverlaps, RankValue)
  Data <- data.frame(Prelim, check.names = FALSE)
  return(Data)
}

#' Internal for QC_ProspectiveAdditions
#'
#' @param NumberDetectors A number corresponding to the number of detectors for
#' the spectral instrument
#'
#' @importFrom utils read.csv
#'
#' @return The stored reference data for that instrument contained within
#' Luciernaga
#' @noRd
InstrumentReferences <- function(NumberDetectors){
  if (NumberDetectors == 184){instrument <- "Sony_7L"
  FileLocation <- system.file("extdata", package = "Luciernaga")
  TheFile <- file.path(FileLocation, "SonyID700_7L_ReferenceLibrary.csv")
  ReferenceData <- read.csv(TheFile, check.names = FALSE)
  } else if (NumberDetectors == "182_DUV"){instrument <- "Sony_6L_DUV"
  FileLocation <- system.file("extdata", package = "Luciernaga")
  TheFile <- file.path(FileLocation, "SonyID700_6LDUV_ReferenceLibrary.csv")
  ReferenceData <- read.csv(TheFile, check.names = FALSE)
  } else if (NumberDetectors == 147){instrument <- "Sony_5L"
  FileLocation <- system.file("extdata", package = "Luciernaga")
  TheFile <- file.path(FileLocation, "SonyID700_5L_ReferenceLibrary.csv")
  ReferenceData <- read.csv(TheFile, check.names = FALSE)
  } else if (NumberDetectors == 112){instrument <- "Sony 4L"
  FileLocation <- system.file("extdata", package = "Luciernaga")
  TheFile <- file.path(FileLocation, "SonyID700_4L_ReferenceLibrary.csv")
  ReferenceData <- read.csv(TheFile, check.names = FALSE)
  } else if (NumberDetectors == 86){instrument <- "Sony 3L"
  FileLocation <- system.file("extdata", package = "Luciernaga")
  TheFile <- file.path(FileLocation, "SonyID700_3L_ReferenceLibrary.csv")
  ReferenceData <- read.csv(TheFile, check.names = FALSE)
  } else if (NumberDetectors == "52_7L"){instrument <- "ThermoFisher 7L_532-594"
  FileLocation <- system.file("extdata", package = "Luciernaga")
  TheFile <- file.path(FileLocation, "ThermoFisherBigfoot7L_532_594_ReferenceLibrary.csv")
  ReferenceData <- read.csv(TheFile, check.names = FALSE)
  } else if (NumberDetectors == 55){instrument <- "ThermoFisher 7L_488-561"
  FileLocation <- system.file("extdata", package = "Luciernaga")
  TheFile <- file.path(FileLocation, "ThermoFisherBigfoot7L_488_561_ReferenceLibrary.csv")
  ReferenceData <- read.csv(TheFile, check.names = FALSE)
  } else if (NumberDetectors == 51){instrument <- "ThermoFisher 6L_785"
  FileLocation <- system.file("extdata", package = "Luciernaga")
  TheFile <- file.path(FileLocation, "ThermoFisherBigfoot6L_785_ReferenceLibrary.csv")
  ReferenceData <- read.csv(TheFile, check.names = FALSE)
  } else if (NumberDetectors == "52_6L"){instrument <- "ThermoFisher 6L_445"
  FileLocation <- system.file("extdata", package = "Luciernaga")
  TheFile <- file.path(FileLocation, "ThermoFisherBigfoot6L_445_ReferenceLibrary.csv")
  ReferenceData <- read.csv(TheFile, check.names = FALSE)
  } else if (NumberDetectors == 78){instrument <- "BD S8"
  FileLocation <- system.file("extdata", package = "Luciernaga")
  TheFile <- file.path(FileLocation, "BDFACSDiscoverS8ReferenceLibrary.csv")
  ReferenceData <- read.csv(TheFile, check.names = FALSE)
  } else if (NumberDetectors == "48_S"){instrument <- "BD S6"
  FileLocation <- system.file("extdata", package = "Luciernaga")
  TheFile <- file.path(FileLocation, "BDFACSSymphonyS6ReferenceLibrary.csv")
  ReferenceData <- read.csv(TheFile, check.names = FALSE)
  } else if (NumberDetectors == "48_A"){instrument <- "BD A5"
  FileLocation <- system.file("extdata", package = "Luciernaga")
  TheFile <- file.path(FileLocation, "BDFACSSymphonyA5ReferenceLibrary.csv")
  ReferenceData <- read.csv(TheFile, check.names = FALSE)
  } else if (NumberDetectors == 64){instrument <- "Cytek 5L"
  FileLocation <- system.file("extdata", package = "Luciernaga")
  TheFile <- file.path(FileLocation, "CytekReferenceLibrary5L.csv")
  ReferenceData <- read.csv(TheFile, check.names = FALSE)
  } else if (NumberDetectors == 54){instrument <- "Cytek 4L UV"
  FileLocation <- system.file("extdata", package = "Luciernaga")
  TheFile <- file.path(FileLocation, "CytekReferenceLibrary4LUV.csv")
  ReferenceData <- read.csv(TheFile, check.names = FALSE)
  } else if (NumberDetectors == 48){instrument <- "Cytek 4L YG"
  FileLocation <- system.file("extdata", package = "Luciernaga")
  TheFile <- file.path(FileLocation, "CytekReferenceLibrary4LYG.csv")
  ReferenceData <- read.csv(TheFile, check.names = FALSE)
  } else if (NumberDetectors == 38){instrument <- "Cytek 3L"
  FileLocation <- system.file("extdata", package = "Luciernaga")
  TheFile <- file.path(FileLocation, "CytekReferenceLibrary3L.csv")
  ReferenceData <- read.csv(TheFile, check.names = FALSE)
  } else if (NumberDetectors == 30){instrument <- "Cytek 2L VB"
  FileLocation <- system.file("extdata", package = "Luciernaga")
  TheFile <- file.path(FileLocation, "CytekReferenceLibrary2LVB.csv")
  ReferenceData <- read.csv(TheFile, check.names = FALSE)
  } else if (NumberDetectors == 22){instrument <- "Cytek 2L BR"
  FileLocation <- system.file("extdata", package = "Luciernaga")
  TheFile <- file.path(FileLocation, "CytekReferenceLibrary2LBR.csv")
  ReferenceData <- read.csv(TheFile, check.names = FALSE)
  } else if (NumberDetectors == 14){instrument <- "Cytek 1L"
  FileLocation <- system.file("extdata", package = "Luciernaga")
  TheFile <- file.path(FileLocation, "CytekReferenceLibrary1L.csv")
  ReferenceData <- read.csv(TheFile, check.names = FALSE)
  } else {message("No References Found")}
}

#' Internal Helper for InstrumentComparison
#' Returns codes for all instruments
#' 
#' @return List of Instrument Codes
#' 
#' @noRd
GetInstruments <- function(){
  Available <- c( 64, 54, 48, 38, 30, 22, 14,
                  184, "182_DUV", 147, 112, 86,
                "52_7L", 55, 51, "52_6L",
              78, "48_A", "48_S")
  return(Available)
}

