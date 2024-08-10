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
#'
#' @return A name matching the provided information for use in plot labeling.
#' @keywords internal
#'
#' @examples NULL

NameForSample <- function(x, sample.name, removestrings, experiment.name = NULL,
                          experiment = NULL, condition.name = NULL, condition = NULL){

  name <- keyword(x, sample.name)
  name <- NameCleanUp(name = name, removestrings)

  if (exists("experiment")){experiment <- experiment
  } else if (exists("experiment.name")){
     if(!is.null(experiment.name)){experiment <- keyword(x, experiment.name)}
   } else {experiment <- NULL}

  if (exists("condition")){condition <- condition
  } else if (exists("condition.name")){
    if(!is.null(condition.name)){condition <- keyword(x, condition.name)}
  } else {condition <- NULL}

  name
  experiment
  condition

  if (is.null(experiment) && is.null(condition)){AggregateName <- name}
  if (!is.null(experiment) && is.null(condition)){AggregateName <- paste(name, experiment, sep="_")}
  if (is.null(experiment) && !is.null(condition)){AggregateName <- paste(name, condition, sep="_")}
  if (!is.null(experiment) && !is.null(condition)){AggregateName <- paste(name, experiment, condition, sep="_")}

  return(AggregateName)
}

