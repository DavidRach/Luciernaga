#' Converts the extracted BeadData into the ArchivedData format
#'  needed to plot Gain/RCV fails
#' 
#' @param data A data.frame object of the bead data
#' @param manufacturer Options Cytek or other, basically whether
#'  to use the DailyQC CSV or the template
#' @param baselinecutoffs The DailyQC .csv or the template .csv
#'  or associated file.path
#' @param returnTemplate Returns DailyQC template that can be
#'  adjusted beyond Cytek settings
#' @param outpath Default NULL, specifies where to returnTemplate
#'  .csv to
#' @param gainmultiplier Gain times this value is the cutoff point
#'  at which Gain Fails
#' 
#' @importFrom stringr str_detect
#' @importFrom dplyr mutate bind_cols
#' @importFrom purrr map
#' 
#' @return An updated data.frame containing the necessary Flag
#'  columns for plotting
#' 
#' @export
HolisticToArchived <- function(data,
                                manufacturer = "Cytek",
                                baselinecutoffs,
                                returnTemplate = FALSE,
                                outpath = NULL,
                                gainmultiplier = 2) {

  Internal <- colnames(data)[str_detect(colnames(data), "-A")]
  RCVs <- Internal[str_detect(Internal, "rCV")]
  Gains <- Internal[str_detect(Internal, "Gain")]

  if (manufacturer == "Cytek") {
    if (returnTemplate == TRUE) {
      if (is.null(outpath)) {
        outpath <- getwd()
      }
      Cutoffs <- Luciernaga:::CytekDailyQC(x = baselinecutoffs,
                                            outpath = outpath,
                                            returnType = "csv")
    } else {
      Cutoffs <- Luciernaga:::CytekDailyQC(x = baselinecutoffs,
                                            returnType = "data")
    }
  } else {
    Cutoffs <- Luciernaga:::NotCytekDailyQC(x = baselinecutoffs)
  }

  Cutoffs1 <- Cutoffs |>
    mutate(GainBaseline = GainBaseline * gainmultiplier)

  TheRCVs <- map(.x = RCVs, .f = DerriveTheFlag, data = data,
                 cutoffs = Cutoffs1) |>
    bind_cols()
  TheGains <- map(.x = Gains, .f = DerriveTheFlag, data = data,
                  cutoffs = Cutoffs1) |>
    bind_cols()
  Assembled <- bind_cols(data, TheRCVs, TheGains)
  return(Assembled)
}