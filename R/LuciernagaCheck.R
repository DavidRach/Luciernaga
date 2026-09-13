#' Generates the plots for InteractiveLuciernaga
#' 
#' @param TheTarget The fluorophore being selected LuciernagaQC outputs
#' @param TheCleanStained The internal list of file.paths to the the
#'  single colors
#' @param returnType Default plotly
#' 
#' @importFrom stringr str_detect
#' @importFrom plotly ggplotly
#' 
#' @return A list of plotly objects for assembling
#' 
#' @noRd
#' 
#' @examples
#' A <- 2 + 2
#' 
LuciernagaCheck <- function(TheTarget, TheCleanStained, returnType = "plotly") {

  TheseFluorSigs <- TheCleanStained[str_detect(TheCleanStained,
                                                paste0(TheTarget, " "))]

  Signatures <- Luciernaga_FolderSignatures(FolderPath = TheseFluorSigs,
                                             sample.name = "GUID",
                                             fluorophore.name = "",
                                             StringRemoval = c(" (Cells)",
                                                                ".fcs"))

  MainDetector <- colSums(Signatures == 1, na.rm = TRUE)
  MainDetector <- names(which.max(MainDetector))

  LinePlot <- QC_ViewSignature(x = NULL, data = Signatures, Normalize = FALSE,
                                columnname = "Sample", legend = TRUE)

  BrightnessPlot <- Luciernaga_FolderBrightness(FolderPath = TheseFluorSigs,
                                                 sample.name = "GUID",
                                                 StringRemoval = c(" (Cells)",
                                                                    ".fcs"),
                                                 fluorophore.name = "",
                                                 returnType = "plot",
                                                 detector = MainDetector,
                                                 PanelCuts = c(0.1, 1))

  if (returnType == "plotly") {
    LinePlot <- ggplotly(LinePlot)
    BrightnessPlot <- ggplotly(BrightnessPlot)
  }

  TheList <- list(LinePlot, BrightnessPlot)

  return(TheList)
}