
#' Internal for LuciernagaReportFromFCS
#'
#' @param x A passed single cytoset object
#' @param Fluorophore The detector
#'
#' @importFrom flowCore keyword exprs
#' @importFrom dplyr mutate relocate
#'
#' @return An internal value
#'
#' @noRd
FCSImportFile <- function(x, Fluorophore, sample.name = "FILENAME"){

  filename <- keyword(x, sample.name)
  filename <- sub(".*\\\\", "", filename)
  filename <- sub(paste0(".*", Fluorophore), Fluorophore, filename)
  filename <- gsub(".fcs$", "", filename)
  rownames(filename) <- NULL

  #df <- exprs(x[[1]])
  df <- exprs(x)
  TheDF <- data.frame(df, check.names = FALSE)
  TheDF <- TheDF[,-grep("Time|FS|SC|SS|Original|W$|H$", names(TheDF))]
  colnames(TheDF) <- gsub("-A$", "", colnames(TheDF))
  DFNames <- TheDF |> mutate(Cluster = filename) |>
    relocate(Cluster, .before = 1)
  return(DFNames)
}
