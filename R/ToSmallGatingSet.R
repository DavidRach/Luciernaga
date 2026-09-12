#' Generates a GatingSet with a positive negative gate on the SC Fluor.
#'
#' @param x A flowframe object
#'
#' @importFrom flowCore exprs keyword
#' @importFrom dplyr summarise across
#' @importFrom tidyselect where
#' @importFrom stats quantile
#' @importFrom flowWorkspace cytoset flowFrame_to_cytoframe cs_add_cytoframe
#'   GatingSet transform flowjo_biexp_trans transformerList
#' @importFrom data.table fread
#' @importFrom openCyto gatingTemplate gt_gating
#'
#' @return A Gating Set object transformed and gated
#' @noRd
ToSmallGatingSet <- function(x) {

  data <- exprs(x)
  data <- data.frame(data, check.names = FALSE)
  data <- data[, -grep("Time|FS|SC|SS|Original|W$|H$", names(data))]
  TheMedian <- data |>
    summarise(
      across(where(is.numeric), \(x) quantile(x, probs = 0.95, na.rm = TRUE))
    )
  KeptMarkers <- colnames(TheMedian)
  Fluorophore <- names(TheMedian)[which.max(TheMedian)]

  MyCS <- cytoset()
  cf <- flowFrame_to_cytoframe(x)

  TheName <- keyword(cf, "$FIL")
  TheName <- unlist(TheName)
  cs_add_cytoframe(MyCS, TheName, cf)
  MyGS <- GatingSet(MyCS)

  MyBiexponentialTransform <- flowjo_biexp_trans(
    channelRange = 256,
    maxValue = 4194303,
    pos = 5.62,
    neg = 0,
    widthBasis = -1000
  )

  TransformList <- transformerList(KeptMarkers, MyBiexponentialTransform)
  UnmixedGatingSet <- transform(MyGS, TransformList)

  FileLocation <- system.file("extdata", package = "Luciernaga")
  UnmixedGates <- fread(
    file.path(path = FileLocation, pattern = "GatesUnmixed.csv")
  )
  Example <- UnmixedGates[6]
  Example[1, 1] <- "Positive"
  Example[1, 2] <- "+"
  Example[1, 3] <- "root"
  Example[1, 4] <- Fluorophore
  Example <- rbind(Example, Example)
  Example[1, 1] <- "Negative"
  Example[1, 2] <- "-"

  UnmixedGating <- gatingTemplate(Example)
  gt_gating(UnmixedGating, UnmixedGatingSet)

  return(UnmixedGatingSet)
}