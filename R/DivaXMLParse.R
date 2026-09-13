#' Takes a BD Diva .xml file and returns data.frame Gain and Laser settings
#' for all fcs files contained within the experiment.
#'
#' @param x A .xml file from BD Diva Software
#'
#' @importFrom xml2 read_xml xml_children xml_find_all xml_name
#' @importFrom purrr map compact
#' @importFrom dplyr bind_rows
#'
#' @return A data.frame row containing the parsed data
#'
#' @noRd
DivaXMLParse <- function(x) {
  Parsed <- read_xml(x)
  Landing <- xml_children(Parsed)
  Experiment <- Landing[xml_name(Landing) == "experiment"][[1]]
  Individual <- xml_find_all(Experiment, ".//specimen[@name]")
  Tubes <- xml_find_all(Individual, ".//tube[@name]")

  Data <- map(.x = Tubes, .f = TubeIterate)
  Data <- compact(Data)
  Data <- bind_rows(Data)

  return(Data)
}