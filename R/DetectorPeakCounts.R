#' Internal for LuciernagaQC
#'
#' @importFrom dplyr select arrange
#' @importFrom tidyselect all_of
#'
#' @return An internal value
#'
#' @noRd
DetectorPeakCounts <- function(x, StartN, EndN){

x <- x |> select(-Backups)
Normalized <- x |> select(all_of(
  StartN:EndN))
Counts <- colSums(Normalized == 1)
PeakDetectorCounts <- data.frame(Fluors = names(Counts), Counts = Counts)
rownames(PeakDetectorCounts) <- NULL
PeakDetectorCounts <- PeakDetectorCounts |>  arrange(desc(Counts))
return(PeakDetectorCounts)
}