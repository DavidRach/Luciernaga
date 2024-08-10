#' Generate a name from the sample for plot identification
#'
#' @param x A GatingSet object
#' @param sample.name The .fcs keyword that contains an unique name for that sample
#' @param removestrings A string of character values to remove from the sample.name
#' @param experiment.name The .fcs keyword that contains the name of the experiment
#' @param experiment Directly provide the name of the experiment (alternative to experiment.name)
#' @param condition.name The .fcs keyword that contains the name of the condition.
#' @param condition Directly provide the name of the condition (alternative to condition.name)
#'
#' @importFrom flowWorkspace keyword
#' @importFrom dplyr pull
#'
#' @return A name matching the provided information for use in plot labeling.
#' @keywords internal
#'
#' @examples NULL

NameForSample <- function(x, sample.name, removestrings, experiment.name = NULL,
                          experiment = NULL, condition.name = NULL, condition = NULL, returnType = "name"){

  name <- keyword(x, sample.name)
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
    } else if (!is.null(experiment) && is.null(condition)){AggregateName <- paste(name, experiment, sep="_")
    } else if (is.null(experiment) && !is.null(condition)){AggregateName <- paste(name, condition, sep="_")
    } else if (!is.null(experiment) && !is.null(condition)){AggregateName <- paste(name, experiment, condition, sep="_")}
    return(AggregateName)
  } else if (returnType == "condition") {return(condition)
    } else if (returnType == "experiment") {return(experiment)}

}

