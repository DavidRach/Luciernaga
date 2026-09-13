#' Internal for LJTracking Parse
#'
#' @param x The starting index
#' @param y The End Index
#' @param ReadInfo The ReadLines object
#'
#' @importFrom stringr str_detect
#' @importFrom dplyr select
#' @importFrom tidyselect all_of
#'
#' @return A data.frame chunk from the original
#'
#' @noRd
ParseThis <- function(x, y, ReadInfo) {
  DetectorSegment <- ReadInfo[x:y]
  DetectorLength <- length(DetectorSegment)
  header <- strsplit(DetectorSegment[1], ",")[[1]]
  data <- strsplit(DetectorSegment, ",")
  NotEmpty <- Filter(function(x) length(x) > 0 && !all(x == ""), data)
  TheData <- do.call(rbind, lapply(NotEmpty, function(x) as.data.frame(t(x),
    stringsAsFactors = FALSE)))
  colnames(TheData) <- header
  TheData <- TheData[-1, ]
  cols_to_select <- which(str_detect(colnames(TheData), "Out of Range Flag"))
  SelectedData <- TheData |> select(-all_of(cols_to_select))
  return(SelectedData)
}