
#' Internal function for relative brightness
#'
#' @param x The passed cluster
#' @param data The passed data.frame
#'
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom tidyselect starts_with
#' @importFrom tidyselect ends_with
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr pull
#' @importFrom tidyselect all_of
#' @importFrom dplyr mutate
#' @importFrom dplyr across
#' @importFrom dplyr coalesce
#'
#' @return An internal value
#'
#' @noRd
FillIterate <- function(x, data){
  IndividualCluster <- data %>% dplyr::filter(Cluster %in% x)
  Detectors <- IndividualCluster %>% select(starts_with("Detector")) %>%
    select(-ends_with("Raw"), -ends_with("Value"))

  Detectors <- Detectors %>% pivot_longer(cols = everything(),
               names_to = "TheNames", values_to = "TheDetectors")

  Detectors <- Detectors %>% filter(!is.na(TheDetectors))

  TheNames <- Detectors %>% pull(TheNames)
  #i <- TheNames[1]

  for (i in TheNames){
    ThisDetector <- i
    ThisColumn <- Detectors %>% dplyr::filter(TheNames %in% i) %>% pull(TheDetectors)
    ThisValue <- IndividualCluster %>% select(all_of(ThisColumn)) %>% pull(.)

    ThisDetector <- paste0(ThisDetector, "Raw")

    IndividualCluster <- IndividualCluster %>% mutate(across(all_of(ThisDetector), ~ coalesce(.x, ThisValue)))
  }

  return(IndividualCluster)
}