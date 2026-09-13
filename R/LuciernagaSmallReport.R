#' Internal for LuciernagaQC
#'
#' @param x TBD
#' @param Data TBD
#' @param RetainedType TBD
#' @param ColsN TBD
#' @param StartNormalizedMergedCol TBD
#' @param EndNormalizedMergedCol TBD
#' @param stats TBD
#'
#' @importFrom dplyr filter select rename
#' @importFrom tidyselect all_of
#'
#' @return An internal value
#'
#' @noRd
LuciernagaSmallReport <- function(x,
                                   Data,
                                   RetainedType,
                                   ColsN,
                                   StartNormalizedMergedCol,
                                   EndNormalizedMergedCol,
                                   stats) {

  if (RetainedType == "raw") {
    Data <- Data |>
      filter(Cluster %in% x) |>
      select(all_of(1:ColsN))
    Averaged <- AveragedSignature(Data, stats)
  }
  if (RetainedType == "normalized") {
    # Data <- Data |> filter(Cluster %in% x) |>
    # select(all_of(StartNormalizedMergedCol:EndNormalizedMergedCol))
    Data <- Data |>
      filter(Cluster %in% x) |>
      select(all_of(1:ColsN))
    Averaged <- AveragedSignature(Data, stats, normalize = TRUE)
  }
  Summary <- cbind(x, Averaged) |> rename(Cluster = x)
  return(Summary)
}