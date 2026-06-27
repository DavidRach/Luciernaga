#' Iterates through the InitialUnmixingSpanGates, allowing for their manual
#' adjustment without needing to search for the corresponding specimen by gate
#' combination
#' 
#' @param template The data.frame containing the name, Fluorophore and Detector columns
#' @param gs The GatingSet
#' @param AdjustAll Whether to adjust all specimens in the gating set to the updated gate,
#'  default is TRUE
#' 
#' @importFrom dplyr select filter pull
#' @importFrom flowWorkspace sampleNames
#' @importFrom purrr pmap
#' 
#' @return No object returned, but updates existing GatingSet gates
#' 
#' @export
#' 
#' @examples A <- 2+2
#' 
CheckUnmixingSpanGates <- function(template, gs, AdjustAll=TRUE){

  These <- template |> select(name, Fluorophore, Detector)
  DetectorsPresent <- These |> filter(!is.na(Detector) & Detector != "")

  gs_index <- match(DetectorsPresent$name, sampleNames(gs))
  DetectorsPresent$gs_index <- gs_index

  GatesToAdd <- DetectorsPresent |> pull(Fluorophore)
  SpecimenIndeces <- DetectorsPresent |> pull(gs_index)

  purrr::pmap(
    .l = list(gate = GatesToAdd, sample = SpecimenIndeces),
    .f = gs_apply_gate_check,
    gs = gs,
    AdjustAll=AdjustAll
  )
}