#' Dashboard Internal, summarizes passing status
#'
#' @param x The iterated parameter
#' @param data The data being compared against
#'
#' @importFrom dplyr filter slice pull
#'
#' @return The color-coded summary
#' @noRd
QCSummaryCheck <- function(x, data) {
  Subset <- data |> filter(Detector %in% x)

  GainValue <- Subset |> slice(1) |> pull(Gain)

  # When a FLAG is TRUE, then we have a problem. 
  if (any(Subset$Gain_Logical == TRUE)) {
    Followup <- Subset |> slice(1) |> pull(Gain_Logical)
    if (Followup == TRUE) {
      GainStatus <- "Red"
    } else {
      GainStatus <- "Yellow"
    }
  } else if (any(Subset$Gain_Logical == FALSE)) {
    GainStatus <- "Green"
  } else {
    GainStatus <- "Gray"
  }

  rCVValue <- Subset |> slice(1) |> pull(rCV) |> round(2)

  if (any(Subset$rCV_Logical == TRUE)) {
    Followup <- Subset |> slice(1) |> pull(rCV_Logical)
    if (Followup == TRUE) {
      rCVStatus <- "Red"
    } else {
      rCVStatus <- "Yellow"
    }
  } else if (any(Subset$rCV_Logical == FALSE)) {
    rCVStatus <- "Green"
  } else {
    rCVStatus <- "Gray"
  }

  Summary <- data.frame(Detector = x, GainValue = GainValue, Gain = GainStatus,
                         rCVValue = rCVValue, rCV = rCVStatus)
  return(Summary)
}