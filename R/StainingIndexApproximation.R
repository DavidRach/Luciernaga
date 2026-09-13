#' Generates Approximate Staining Index and Associated Info
#'
#' @param x The GatingSet object with applied transformation and gates.
#' @param NumberDetectors The number detectors the original file was based on.
#'
#' @importFrom flowWorkspace gs_pop_get_data
#' @importFrom flowCore exprs keyword
#' @importFrom dplyr summarise across select pull filter
#' @importFrom tidyselect where all_of
#' @importFrom stats quantile
#' @importFrom tidyr pivot_wider
#'
#' @return A data.frame row containing the derrived information
#' @noRd
StainingIndexApproximation <- function(x, NumberDetectors) {
  name <- keyword(x, "$FIL")

  Positive <- gs_pop_get_data(x, "Positive", inverse.transform = FALSE)
  Positive <- exprs(Positive[[1]])
  Positive <- data.frame(Positive, check.names = FALSE)
  Positive <- Positive[, -grep("Time|FS|SC|SS|Original|W$|H$",
                                names(Positive))]
  TheMedian <- Positive |>
    summarise(across(where(is.numeric),
                      \(x) quantile(x, probs = 0.95, na.rm = TRUE)))
  KeptMarkers <- colnames(TheMedian)
  Fluorophore <- names(TheMedian)[which.max(TheMedian)]
  Positive <- Positive |> select(all_of(Fluorophore))

  Negative <- gs_pop_get_data(x, "Negative", inverse.transform = FALSE)
  Negative <- exprs(Negative[[1]])
  Negative <- data.frame(Negative, check.names = FALSE)
  Negative <- Negative[, -grep("Time|FS|SC|SS|Original|W$|H$",
                                names(Negative))]
  Negative <- Negative |> select(all_of(Fluorophore))

  PositiveVals <- Positive |> pull(1)
  MFI_Pos <- quantile(PositiveVals, probs = 0.5, na.rm = TRUE)
  MFI_Pos <- MFI_Pos[[1]]
  NegativeVals <- Negative |> pull(1)
  MFI_Neg <- quantile(NegativeVals, probs = 0.5, na.rm = TRUE)
  MFI_Neg <- MFI_Neg[[1]]
  # hist(NegativeVals)
  MinNeg <- quantile(NegativeVals, probs = 0.03, na.rm = TRUE)
  MinNeg <- MinNeg[[1]]
  MaxNeg <- quantile(NegativeVals, probs = 0.97, na.rm = TRUE)
  MaxNeg <- MaxNeg[[1]]
  RSD <- (MaxNeg - MinNeg) / 3.29
  StainingIndex <- (MFI_Pos - MFI_Neg) / (2 * RSD)

  ReferenceData <- Luciernaga:::InstrumentReferences(
    NumberDetectors = NumberDetectors)
  Matrix <- ReferenceData |>
    filter(Fluorophore %in% KeptMarkers) |>
    select(-Instrument) |>
    pivot_wider(names_from = "Detector", values_from = "AdjustedY") |>
    select(-Fluorophore) |>
    as.matrix()
  kappa <- round(kappa(Matrix, exact = TRUE), 2)

  MatrixN <- length(KeptMarkers)
  MatrixMarkers <- paste(KeptMarkers, collapse = "_")
  FinalData <- cbind(name, Fluorophore, MatrixN, MatrixMarkers,
                      StainingIndex, kappa)
  colnames(FinalData)[1] <- "name"

  return(FinalData)
}