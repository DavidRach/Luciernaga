#' Takes a folder of Luciernaga_QC signatures and returns a
#'  data.frame of the signatures contained within, for use
#'  in downstream data analysis.
#'
#' @param FolderPath Location of the Luciernaga_QC .fcs outputs
#' @param sample.name The keyword where the identifying sample name
#' can be found
#' @param StringRemoval Default NULL, provide to remove items from
#'  sample.name
#' based on values found on the keyword
#' @param fluorophore.name Specify the name of the fluorophore,
#'  alternatively NULL
#' @param Verbose Default FALSE, provides info as it goes
#' @param stats Whether to use median or mean
#' @param PanelCuts Default NULL, provide a c(0.5,1) argument to specify
#' the brightness percentiles to retrieve signature from for the
#' individual files
#' @param normalize Default TRUE, whether to return normalized or
#' raw averaged MFI signatures
#' @param returnType Default is "Signatures"
#'
#' @importFrom flowWorkspace load_cytoset_from_fcs GatingSet
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
#'
#' @export
#'
#' @return A data.frame of the corresponding averaged signatures per file
#'
#' @examples
#' A <- 2 + 2
Luciernaga_FolderSignatures <- function(FolderPath, sample.name,
                                         StringRemoval = NULL,
                                         fluorophore.name, Verbose = FALSE,
                                         stats = "median", PanelCuts = NULL,
                                         normalize = TRUE,
                                         returnType = "Signatures") {
  if (length(FolderPath > 1)) {
    TheFCSFiles <- FolderPath
  } else {
    TheFCSFiles <- list.files(path = FolderPath, pattern = "fcs",
      full.names = TRUE)
  }

  Selected_CS <- load_cytoset_from_fcs(TheFCSFiles,
    truncate_max_range = FALSE, transformation = FALSE)
  Selected_GS <- GatingSet(Selected_CS)

  Returns <- map(.x = Selected_GS, .f = FolderSignatureIterator,
    sample.name = sample.name, StringRemoval = StringRemoval,
    fluorophore.name = fluorophore.name, Verbose = Verbose,
    stats = stats, PanelCuts = PanelCuts, normalize = normalize,
    returnType = returnType) |> bind_rows()

  return(Returns)
}