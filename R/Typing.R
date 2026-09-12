
#' Internal for LuciernagaQC, returns unmixing control type for downstream forking.
#'
#' @param name Passed from edits to sample.name
#' @param unmixingcontroltype Whether beads or cells
#' @param Unstained When set true pastes Unstained for downstream use
#'
#' @importFrom stringr str_detect
#'
#' @return An internal value
#'
#' @noRd
Typing <- function(name, unmixingcontroltype, Unstained = FALSE){
  if (unmixingcontroltype == "beads"){
    if (!str_detect(name, "ells")){Type <- "Beads"
    }
  }

  if (unmixingcontroltype == "cells"){
    if (!str_detect(name, "eads")) {Type <- "Cells"
    }
  }

  if (unmixingcontroltype == "both") {
    if (str_detect(name, "ells")){Type <- "Cells"
    } else if(str_detect(name, "eads")){Type <- "Beads"
    } else {Type <- "Unknown"}
  }

  if(!exists("Type")){
    if (str_detect(name, "ells")){Type <- "Cells"
    } else if(str_detect(name, "eads")){Type <- "Beads"
    } else {Type <- "Unknown"}
  }


  # Figuring out if Unstained
  if (Unstained == TRUE) {
    if(!str_detect(name, "stained")){Type <- paste0(Type, "_Unstained")}
  } else {
    if (str_detect(name, "stained")){Type <- paste0(Type, "_Unstained")}
  }
  return(Type)
}