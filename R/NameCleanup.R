#' Shorten a variable name
#'
#' @param name The variable that you wish to shorten.
#' @param removestrings A list of strings to be sequentially removed from the name.
#'  An example: removestrings = c("ReferenceGroup_", ".fcs")
#' @param substitutestrings A data.frame containing two columns, Original and Substitute.
#'  When a row of Original is recognized, it's replaced with substitute value.
#'
#' @importFrom stringr str_replace_all
#'
#' @return The shortened variable
#' @export
#'
#' @examples NULL

NameCleanUp <- function(name, removestrings, substitutestrings){
#InternalDevToDoList: The for-loop is iterating in place, don't believe its subject to replacing
# Rather than specifying everything, provide a data.frame for substitutions. Have it iterate afterwards.

  for(i in removestrings){
    name <- str_replace_all(name, i, "")
  }

  return(name)
}


