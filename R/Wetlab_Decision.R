#' Takes Rest output, and referencing TheCondition csv makes decision for culture resuspensions
#'
#' @param data The Wetlab_Rest output
#' @param FinalConcentration_MillionperML Desired Final Concentration in Millions
#' @param MillionCellsPerTube Desired number of cells in Tube
#' @param TheConditions A filepath or data.frame to TheCondition csv
#' @param ReturnLeftover Whether to return leftover cells as own line
#'
#' @importFrom utils read.csv
#' @importFrom dplyr select filter mutate case_when bind_rows
#' @importFrom stringr str_detect
#' @importFrom purrr map
#'
#' @return A data.frame object
#' @export
#'
#' @examples NULL
Wetlab_Decision <- function(data, FinalConcentration_MillionperML, MillionCellsPerTube, TheConditions,
                            ReturnLeftover=TRUE){

  if(!is.data.frame(TheConditions)){TheConditions <- read.csv(TheConditions, check.names = FALSE)}

  data <- data %>% select(-c(CurrentConcentration, TotalVolume, IncreaseVolumeML))

  SingleColorStash <- data %>% filter(Specimen == FALSE)
  NotSingle <- data %>% filter(!Specimen == FALSE)

  Internal <- NotSingle %>% filter(!SpinDown %in% TRUE) %>% mutate(
    name = case_when(str_detect(name, "Spin") ~ gsub("Spin_", "", name), TRUE ~ name))

  Specimens <- Internal$name

  #Remove the single specimen select below when done
  TheReturn <- map(.x=Specimens, .f=DecisionInternal, data=Internal,
                   FinalConcentration_MillionperML=FinalConcentration_MillionperML,
                   MillionCellsPerTube=MillionCellsPerTube,
                   TheConditions=TheConditions, ReturnLeftover=ReturnLeftover) %>% bind_rows()

  return(TheReturn)
}








