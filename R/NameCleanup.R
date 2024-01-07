#' Polish up single color names before unmixing
#'
#' @param name Original name
#' @param removestrings Strings to be removed from name (ex. c("DTR_ILT_", "cord"))
#'
#' @return NULL
#' @export
#'
#' @examples NULL
NameCleanUp <- function(name, removestrings){
  for(i in removestrings){
    name <- str_replace_all(name, i, "")
  }
  name <- gsub("_.*", "", name)
  name <- gsub("BUV", "_BUV", gsub("BV", "_BV", gsub("APC", "_APC", gsub("PE", "_PE", gsub(
    "Per", "_Per", gsub("FITC", "_FITC", gsub("Spark", "_Spark", gsub("Alexa", "_Alexa", gsub(
      "Pacific", "_Pacific", gsub("Zombie", "_Zombie", name))))))))))
  name <- gsub("Samples", "Unstained", gsub("FITCC", "C", name))
  return(name)
}
