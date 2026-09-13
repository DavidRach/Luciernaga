#' Internal for Chorizo, returns Laser parameters for Description
#'
#' @param x The iterated laser
#' @param data The filtered laser list
#'
#' @importFrom dplyr filter pull
#'
#' @return The list of laser parameters for the Description list
#'
#' @noRd
LaserList <- function(x, data) {

  Subset <- data |> filter(LaserOrder %in% x)
  Number <- Subset |> pull(LaserNumber)

  LaserASF <- paste0("LASER", Number, "ASF")
  LaserDelay <- paste0("LASER", Number, "DELAY")
  LaserName <- paste0("LASER", Number, "NAME")

  LaserASFVal <- Subset |> pull(ASF)
  LaserDelayVal <- Subset |> pull(DELAY)
  LaserNameVal <- Subset |> pull(NAME)

  LaserList <- list(
    LaserASF = LaserASFVal,
    LaserDelay = LaserDelayVal,
    LaserName = LaserNameVal
  )

  # Optionally, set the names programmatically if dynamic
  names(LaserList) <- c(LaserASF, LaserDelay, LaserName)

  return(LaserList)
}