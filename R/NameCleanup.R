#' Shorten a variable name
#'
#' @param name The variable that you wish to shorten.
#' @param removestrings A list of strings to be sequentially removed from the name.
#'  An example: removestrings = c("ReferenceGroup_", ".fcs")
#' @param substitutestrings A data.frame containing two columns, Original and
#' Substitute.
#'  When a row of Original is recognized, it's replaced with substitute value.
#'
#' @importFrom stringr str_replace_all
#' @importFrom stringr fixed
#'
#' @return The shortened variable
#'
#' @export
#'
#' @examples
#' name <- "DR BUV496 CD8 (Cells).fcs"
#' removestrings <- c("DR", "(Cells)", ".fcs", " ")
#' Cleaned_Name <- NameCleanUp(name, removestrings)
NameCleanUp <- function(name, removestrings, substitutestrings){
  for(i in removestrings){
    name <- str_replace_all(name, fixed(i), "")
  }

  return(name)
}

#' Generate a name from the sample for plot identification
#'
#' @param x A GatingSet object
#' @param sample.name The .fcs keyword that contains an unique name for that sample
#' @param removestrings A string of character values to remove from the sample.name
#' @param experiment.name The .fcs keyword that contains the name of the experiment
#' @param experiment Directly provide the name of the experiment (alternative to
#' experiment.name)
#' @param condition.name The .fcs keyword that contains the name of the condition.
#' @param condition Directly provide the name of the condition (alternative to
#' condition.name)
#'
#' @importFrom flowWorkspace keyword
#' @importFrom dplyr pull
#'
#' @return A name matching the provided information for use in plot labeling.
#'
#' @noRd
NameForSample <- function(x, sample.name, removestrings, experiment.name = NULL,
                          experiment = NULL, condition.name = NULL, condition = NULL,
                          returnType = "name"){

  if (length(sample.name) == 2){
    first <- sample.name[[1]]
    second <- sample.name[[2]]
    first <- keyword(x, first)
    second <- keyword(x, second)
    name <- paste(first, second, sep="_")
  } else {name <- keyword(x, sample.name)}

  name <- NameCleanUp(name = name, removestrings)

  if (exists("experiment")){
    if (!is.null(experiment)) {experiment <- experiment
    } else if (is.null(experiment)) {Internal <- "A"}
  }


  if (exists("experiment.name")){
    if(!is.null(experiment.name)){experiment <- keyword(x, experiment.name)
    experiment <- pull(experiment)
    } else if (is.null(experiment.name)) {Internal <- "A"}
  }

  if (exists("condition")){
    if (!is.null(condition)) {condition <- condition
    } else if (is.null(condition)) {Internal <- "A"}
  }

  if (exists("condition.name")){
    if(!is.null(condition.name)){condition <- keyword(x, condition.name)
    condition <- pull(condition)
    } else if (is.null(condition.name)) {Internal <- "A"}
  }

  if (returnType == "name"){
    if (is.null(experiment) && is.null(condition)){AggregateName <- name
    } else if (!is.null(experiment) && is.null(condition)){
      AggregateName <- paste(name, experiment, sep="_")
    } else if (is.null(experiment) && !is.null(condition)){
      AggregateName <- paste(name, condition, sep="_")
    } else if (!is.null(experiment) && !is.null(condition)){
      AggregateName <- paste(name, experiment, condition, sep="_")}
    return(AggregateName)
  } else if (returnType == "condition") {return(condition)
  } else if (returnType == "experiment") {return(experiment)}

}

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


