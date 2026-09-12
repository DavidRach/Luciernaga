#' Based on a template, creates the initial gates for use in retrieving the 
#' signature matrix. Expects columns name, Fluorophore and Detector
#' 
#' @param template The data.frame containing the name, Fluorophore and Detector columns
#' @param gs The GatingSet
#' @param subset The gate under which to create the new gates
#' @param minpercentile Default 0.90, sets the lower bound of the span gate being created
#' @param maxpercentile Default 0.99, sets the upper bound of the span gate being created
#' @param inverse.transform Whether to inverse a transformation, default set to FALSE
#' 
#' @importFrom dplyr select filter 
#' @importFrom purrr walk
#' @importFrom flowWorkspace sampleNames
#' 
#' 
#' @return Silent, creates gates in the GatingSet
#' 
#' @export
#' 
#' @examples A <- 2+2
#' 
InitialUnmixSpanGates <- function(template, gs, subset, minpercentile=0.5,
 maxpercentile=0.99, inverse.transform=FALSE){

    These <- template |> select(name, Fluorophore, Detector)
    DetectorsPresent <- These |> filter(!is.na(Detector) & Detector != "")

    gs_index <- match(DetectorsPresent$name, sampleNames(gs))
    DetectorsPresent$gs_index <- gs_index

    GatesToAdd <- DetectorsPresent |> pull(Fluorophore)

    purrr::walk(.x= GatesToAdd, .f=InitialSpans, gs=gs, subset=subset, data=DetectorsPresent,
    inverse.transform=inverse.transform, minpercentile=minpercentile,
    maxpercentile=maxpercentile, .progress = TRUE)
 }