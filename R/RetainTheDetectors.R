#' Internal for SC_Unmix
#'
#' @param AFOverlap The Overlap list to remove AF Detecotrs
#' @param Detectors The passed PeakDetectorList to decide on
#' @param name The passed name consisting ligand fluorophore
#'
#' @importFrom utils read.csv
#' @importFrom dplyr filter pull
#' @importFrom stringr str_split str_detect
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