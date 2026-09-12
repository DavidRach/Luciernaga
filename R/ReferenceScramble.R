#' Internal for SC_Unmix
#'
#' @param name The passed ligand fluorophore name
#' @param NumberDetectors The number of detectors for referencing
#' @param NumberFluors Desired number of additional fluorophores in the matrix
#'
#' @importFrom stringr str_count
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
  if (Spaces == 0){message("Provide argument like `CD4 BUV805` for correct formatting")
  } else if (Spaces <= 1){Name <- strsplit(name, " ")[[1]]
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