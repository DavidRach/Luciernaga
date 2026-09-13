#' Internal for Chorizo, generates a parameter data.frame from provided
#' data.
#'
#' @param x A data.frame with parameter column names in the correct order.
#'
#' @importFrom dplyr mutate pull case_when row_number
#'
#' @return A data.frame with the parameter data
#' @noRd
ParameterPrep <- function(x) {
  TheNames <- colnames(x)
  ParameterStandin <- data.frame(name = TheNames, desc = NA, range = 4194303,
    minRange = -111.00000, maxRange = 4194303)
  TheParams <- ParameterStandin |>
    mutate(TheParams = paste0("$P", row_number())) |> pull(TheParams)
  rownames(ParameterStandin) <- TheParams
  ParameterStandin$desc <- as.character(ParameterStandin$desc)

  ParameterStandin <- ParameterStandin |> mutate(
    range = case_when(name == "Time" ~ 532116, TRUE ~ range),
    minRange = case_when(name == "Time" ~ 0.00000, TRUE ~ minRange),
    minRange = case_when(name == "SSC-W" ~ 0.00000, TRUE ~ minRange),
    minRange = case_when(name == "SSC-H" ~ 0.00000, TRUE ~ minRange),
    minRange = case_when(name == "SSC-A" ~ 0.00000, TRUE ~ minRange),
    minRange = case_when(name == "FSC-W" ~ 0.00000, TRUE ~ minRange),
    minRange = case_when(name == "FSC-H" ~ 0.00000, TRUE ~ minRange),
    minRange = case_when(name == "FSC-A" ~ 0.00000, TRUE ~ minRange),
    minRange = case_when(name == "SSC-B-W" ~ 0.00000, TRUE ~ minRange),
    minRange = case_when(name == "SSC-B-H" ~ 0.00000, TRUE ~ minRange),
    minRange = case_when(name == "SSC-B-A" ~ 0.00000, TRUE ~ minRange),
    minRange = case_when(name == "UV8-A" ~ -91.25813, TRUE ~ minRange),
    minRange = case_when(name == "V6-A" ~ 0.00000, TRUE ~ minRange),
    minRange = case_when(name == "V7-A" ~ 0.00000, TRUE ~ minRange),
    minRange = case_when(name == "V8-A" ~ 0.00000, TRUE ~ minRange),
    minRange = case_when(name == "V9-A" ~ 0.00000, TRUE ~ minRange),
    minRange = case_when(name == "V10-A" ~ 0.00000, TRUE ~ minRange),
    minRange = case_when(name == "B3-A" ~ 0.00000, TRUE ~ minRange)
  )

  return(ParameterStandin)
}