#' Internal for SingleUnmix, drops the given single-color signature,
#'  unmixes full panel without it.
#'
#' @param x The iterated in fluorophore name, matching metadata in
#'  gs and fluorophore in matrix
#' @param gs The GatingSet containing the raw single-color unmixing
#'  controls
#' @param matrix The signature matrix for the respective single-color
#'  unmixing controls
#' @param outpath The location to store the unmixed .fcs files.
#' @param subset The gate with the events you want to unmix. Default
#'  is root.
#' @param inverse.transform Default equals TRUE
#' @param sample.name Keyword(s) to grab from the cytoset for renaming
#'
#' @importFrom BiocGenerics subset
#' @importFrom dplyr filter select
#' @importFrom tidyselect where
#' @importFrom purrr walk
#'
#' @return Unmixed single-color unmixing controls for use in stain
#'  reduction calculations.
#'
#' @noRd
FMOUnmixIterator <- function(x, gs, matrix, outpath, subset,
                              inverse.transform, sample.name) {
  # x <- TheseFluorophores[1]
  InternalGS <- subset(gs, Fluorophore != x)
  Signature <- matrix |> dplyr::filter(!Fluorophore %in% x) |>
    dplyr::select(Fluorophore, Antigen, where(is.numeric))
  Panel <- matrix |> dplyr::filter(!Fluorophore %in% x) |>
    dplyr::select(Fluorophore, Antigen)
  colnames(Signature)[2] <- "Ligand"
  #colnames(Panel)[2] <- "Ligand"

  TheName <- paste0("No", x)
  TheNameUnmixed <- paste0("_", TheName, "Unmixed")

  StashHere <- file.path(outpath, TheName)
  if (!dir.exists(StashHere)) {
    dir.create(StashHere)
  }

  walk(.x = InternalGS, .f = Luciernaga_Unmix, controlData = Signature,
    sample.name = sample.name,
    addon = TheNameUnmixed, subset = subset, removestrings = "fcs",
    outpath = StashHere, PanelPath = Panel,
    Verbose = FALSE, inverse.transform = inverse.transform)

  #devtools::load_all("/home/david/Documents/Luciernaga")
}