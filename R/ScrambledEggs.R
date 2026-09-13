#' Wrapper Function for Scrambled Eggs Protocol
#'
#' @param NumberRepeats Desired number of bootstrap runs
#' @param NumberFluors Desired Number of Additional Fluorophores
#' @param NumberDetectors Detector configuration original FCS file
#' @param GS The GatingSet object corresponding to desired SC
#' @param sample.name The keyword designating single-color sample name
#' @param removestrings Values to be removed to leave just the sample name
#' @param subset The desired gating node
#' @param multiplier The multiplier for unmixing, default set to 50000
#' @param outpath Internal outpath argument
#' @param addon Internal unmix argument
#'
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
#'
#' @return A data.frame summarizing the Staining Index and Kappa for each bootstrap
#' @noRd
ScrambledEggs <- function(NumberRepeats,
                           NumberFluors,
                           NumberDetectors,
                           GS,
                           sample.name = "TUBENAME",
                           removestrings,
                           subset = "lymphocytes",
                           multiplier = 50000,
                           outpath,
                           addon = addon) {
  Data <- list()
  for (i in seq_along(1:NumberRepeats)) {
    Data[[i]] <- Luciernaga:::SC_Unmix(x = GS, sample.name = sample.name,
                                        removestrings = removestrings,
                                        subset = subset,
                                        multiplier = multiplier,
                                        outpath = outpath,
                                        returntype = "data", Verbose = FALSE,
                                        addon = "_Unmixed",
                                        ratiopopcutoff = 0.01,
                                        NumberFluors = NumberFluors)
  }
  GatedData <- map(.x = Data, .f = Luciernaga:::ToSmallGatingSet)
  MyData <- map(.x = GatedData, .f = Luciernaga:::StainingIndexApproximation,
                NumberDetectors = NumberDetectors) |>
    bind_rows()
  return(MyData)
}