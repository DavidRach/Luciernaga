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