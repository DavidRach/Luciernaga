#' Similar to CytosetScreen, checks for mismatching cytoframes that
#'  throw inconvenient errors
#' 
#' @param gs A gating set object
#' @param gatingtemplate The gating template object
#' @param subsets Define what level to check
#' 
#' @return Purified Gating Set or a NULL Value
#' 
#' @noRd
GateCheck <- function(gs, gatingtemplate, subsets=NULL){
  
  Nodes <- gatingtemplate@nodes
  These <- Nodes
  #FinalNode <- Nodes[length(Nodes)]

  if (!is.null(subsets)){
    TheMatch <- paste0("/", subsets)
    Matching <- min(which(These %in% TheMatch))
    These <- These[1:Matching]
    if (length(These) == 1){CheckThis <- These
    } else {CheckThis <- paste(These, collapse="|")}
  } else {
    These <- These[-1]
    CheckThis <- paste(These, collapse="|")
  }

  ThisFluor <- gatingtemplate@edgeData@data[[CheckThis]]$gtMethod@dims

  TheIndex <- any(colnames(gs) %in% ThisFluor)
  Present <- gs[TheIndex]

  if (length(Present) == 0){gs <- NULL
  } else {gs <- Present}
  return(gs)
}