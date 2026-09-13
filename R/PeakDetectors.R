#' Internal for SC_Unmix
#'
#' @param NewData The exprs for the .fcs file minus scatter params
#'
#' @importFrom dplyr arrange desc
#'
#' @return A data.frame of detectors and respective counts
#' @noRd
PeakDetectors <- function(NewData) {
  NewData[NewData < 0] <- 0
  A <- do.call(pmax, NewData)
  Normalized <- NewData / A
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
  PeakDetectorCounts <- PeakDetectorCounts |> arrange(desc(Counts))
  PeakList <- list(WorkAround, PeakDetectorCounts)
  return(PeakList)
}