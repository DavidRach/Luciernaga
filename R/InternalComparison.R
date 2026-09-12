
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

  TheCandidate <- ReferenceData |> dplyr::filter(Fluorophore %in% x)
  TheReferenceList <- ReferenceData |> dplyr::filter(Fluorophore %in% TheList)
  TheData <- rbind(TheCandidate, TheReferenceList) |> select(-Instrument)
  TheCosineData <- TheData |>
    pivot_wider(names_from = Detector, values_from = AdjustedY)
  Names <- TheCosineData |> pull(Fluorophore)
  TheValues <- TheCosineData |> select(-Fluorophore)
  TheTransposed <- t(TheValues)
  colnames(TheTransposed) <- Names
  TheMatrix <- data.matrix(TheTransposed)
  CosineMatrix <- cosine(TheMatrix)
  CosineMatrix <- data.frame(CosineMatrix, check.names = FALSE)
  TheCandidateValues <- CosineMatrix |> select(all_of(x))
  HighOverlaps <- TheCandidateValues %>% filter(!!sym(x) >= TheCutoff)
  HighOverlaps <- length(HighOverlaps)
  RankValue <- round(kappa(CosineMatrix, exact=TRUE),2) # I believe from Cosine, not the signature matrix
  Fluorophore <- x

  Prelim <- cbind(Fluorophore, TheDetector, HighOverlaps, RankValue)
  Data <- data.frame(Prelim, check.names = FALSE)
  return(Data)
}