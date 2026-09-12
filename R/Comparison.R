
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

  # x <- TheComparisonList[1]
  TheIndividualDetector <- map(.x=TheComparisonList, .f=InternalComparison,
                               TheList=TheList, ReferenceData=ReferenceData,
                               TheCutoff=TheCutoff, TheDetector=TheDetector) |>
    bind_rows()

  return(TheIndividualDetector)
}