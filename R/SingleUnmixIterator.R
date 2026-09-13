#' Internal for SingleUnmix, unmixes with just single reference signature
#' 
#' @param x The iterated in fluorophore name, matching metadata in gs and fluorophore in matrix
#' @param matrix The signature matrix for the respective single-color unmixing controls
#' @param gs The GatingSet containing the raw single-color unmixing controls
#' @param subset The gate with the events you want to unmix. Default is root. 
#' @param outpath The location to store the unmixed .fcs files. 
#' @param inverse.transform Default equals TRUE
#' @param sample.name Keyword(s) to grab from the cytoset for renaming
#' 
#' @importFrom BiocGenerics subset
#' @importFrom dplyr filter select
#' @importFrom tidyselect where
#' @importFrom purrr walk
#' 
#' @return Unmixed single-color unmixing controls for use in stain reduction calculations. 
#' 
#' @noRd
SingleUnmixIterator <- function(x,
                                 gs,
                                 matrix,
                                 outpath,
                                 subset,
                                 inverse.transform,
                                 sample.name) {
  # x <- TheseFluorophores[1]
  InternalGS <- subset(gs, Fluorophore == x)
  Signature <- matrix |>
    dplyr::filter(Fluorophore %in% x) |>
    dplyr::select(Fluorophore, Antigen, where(is.numeric))
  Panel <- matrix |>
    dplyr::filter(Fluorophore %in% x) |>
    dplyr::select(Fluorophore, Antigen)
  colnames(Signature)[2] <- "Ligand"
  # colnames(Panel)[2] <- "Ligand"

  walk(.x = InternalGS, .f = Luciernaga_Unmix, controlData = Signature,
       sample.name = sample.name, addon = "_SingleUnmixed", subset = subset,
       removestrings = "fcs", outpath = outpath, PanelPath = Panel,
       Verbose = FALSE, inverse.transform = inverse.transform)

  # devtools::load_all("/home/david/Documents/Luciernaga")
}