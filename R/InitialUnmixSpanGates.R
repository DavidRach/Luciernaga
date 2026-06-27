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

#' Internal for InitialUnmixingGates, creates gates for designated fluorophores
#' at the corresponding dimension, using the assigned percentiles to draw
#' the initial span gates
#' 
#' @param x The iterated in Fluorophore
#' @param gs The GatingSet
#' @param subset The gate that will be the parent to the new one
#' @param data The template containing the relavent information
#' @param inverse.transform Whether to inverse.transform the exprs data
#' @param minpercentile The lower percentile bound to start the span gate at
#' @param maxpercentile The upper percentile bound to start the span gate at
#' 
#' @importFrom dplyr filter pull
#' @importFrom flowWorkspace gs_pop_get_data gs_get_pop_paths gs_pop_remove
#' @importFrom openCyto gs_add_gating_method
#' @importFrom stringr str_equal
#' @importFrom Biobase exprs
#' @importFrom stats quantile
#' 
#' @return Nothing, but creates the gate
#' 
#' @examples A <- 2+2
#' 
#' @noRd
#' 
InitialSpans <- function(x, gs, subset, data, inverse.transform,
 minpercentile, maxpercentile){

    Internal <- data |> filter(Fluorophore %in% x)
    filterId <- Internal |> pull(Fluorophore)
    dims <- Internal |> pull(Detector)
    theIndex <- Internal |> pull(gs_index)

    InternalData <- gs_pop_get_data(gs[theIndex], subset=subset, inverse.transform=inverse.transform)
    TheExprs <- data.frame(exprs(InternalData[[1]]), check.names=FALSE)
    ExprDim <- paste0(dims, "-A")
    TheValues <- TheExprs[[ExprDim]]
    MinBoundary <- quantile(TheValues, minpercentile) |> unname()
    MaxBoundary <- quantile(TheValues, maxpercentile) |> unname()
    Values <- paste(c(MinBoundary, MaxBoundary), collapse=",")

    ExistingGates <- gs_get_pop_paths(gs, path="auto")

    if(!any(str_equal(ExistingGates, filterId))){

        suppressMessages(
            gs_add_gating_method(gs, alias = filterId, pop = "+", parent = subset, 
                     dims = ExprDim, gating_method = "span_gate",
                     gating_args = Values)
        )

    } else {
        gs_pop_remove(gs, filterId, recompute = TRUE, recursive = TRUE) 

        suppressMessages(
        gs_add_gating_method(gs, alias = filterId, pop = "+", parent = subset, 
                     dims = ExprDim, gating_method = "span_gate",
                     gating_args = Values)
        )

    }
}