#' Returns Detector Gains, Laser Delay and Scalings from individual .fcs files in a
#' CytoSet.
#'
#' @param x A CytoSet (when mapped) or an individual cytoset, example (x=MyCytoSet[[1]])
#' @param sample.name The keyword value that distinguishes individual .fcs files
#'
#' @importFrom flowCore keyword
#' @importFrom purrr map2
#' @importFrom dplyr bind_cols
#' @importFrom dplyr mutate
#' @importFrom dplyr across
#' @importFrom tidyselect everything
#' @importFrom lubridate dmy
#' @importFrom lubridate hms
#'
#' @return A dataframe row
#' @export
#'
#' @examples NULL
QC_Retrieval <- function(x, sample.name){
  KeywordsList <- keyword(x)
  KeywordsDF <- data.frame(KeywordsList, check.names = FALSE)
  TheColumnNames <- colnames(KeywordsDF)

  SAMPLE <- keyword(x)[[sample.name]]
  DATE <- keyword(x)$`$DATE`
  DATE <- dmy(DATE)
  TIME <- keyword(x)$`$BTIM`
  TIME <- hms(TIME)
  CYT <- keyword(x)$`$CYT`
  CYTSN <- keyword(x)$`$CYTSN`
  OP <- keyword(x)$`$OP`

  PN_Names1 <- TheColumnNames[grepl("^\\$P[0-9]{1}N$", TheColumnNames)]
  PN_Names2 <- TheColumnNames[grepl("^\\$P[0-9]{2}N$", TheColumnNames)]
  PV_Gains1 <- TheColumnNames[grepl("^\\$P[0-9]{1}V$", TheColumnNames)]
  PV_Gains2 <- TheColumnNames[grepl("^\\$P[0-9]{2}V$", TheColumnNames)]
  PN_Names <- c(PN_Names1, PN_Names2)
  PN_Names <- PN_Names[-1] # Remove Time
  PV_Gains <- c(PV_Gains1, PV_Gains2)

  ParameterRows <- map2(.x=PN_Names, .y=PV_Gains, .f=RetrievalMerge,
                        TheData=KeywordsDF) %>% bind_cols()

  Laser_Name <- TheColumnNames[grepl("^\\LASER[0-9]{1}NAME$", TheColumnNames)]
  Laser_Delay <- TheColumnNames[grepl("^\\LASER[0-9]{1}DELAY$", TheColumnNames)]
  Laser_ASF <- TheColumnNames[grepl("^\\LASER[0-9]{1}ASF$", TheColumnNames)]

  LaserDelayRows <- map2(.x=Laser_Name, .y=Laser_Delay, .f=RetrievalMerge,
                         TheData=KeywordsDF) %>% bind_cols()
  colnames(LaserDelayRows) <- paste0(colnames(LaserDelayRows), "_LaserDelay")
  LaserASFRows <- map2(.x=Laser_Name, .y=Laser_ASF, .f=RetrievalMerge,
                       TheData=KeywordsDF) %>% bind_cols()
  colnames(LaserASFRows) <- paste0(colnames(LaserASFRows), "_AreaScalingFactor")

  ParameterRows <- ParameterRows %>% mutate(across(everything(), as.numeric))
  colnames(ParameterRows) <- paste0(colnames(ParameterRows), "_Gain")
  LaserDelayRows <- LaserDelayRows %>% mutate(across(everything(), as.numeric))
  LaserASFRows <- LaserASFRows %>% mutate(across(everything(), as.numeric))

  RecoveredQC <- cbind(SAMPLE, DATE, TIME, CYT, CYTSN, OP, ParameterRows,
                       LaserDelayRows, LaserASFRows)

  return(RecoveredQC)
}

#' Internal for QC_Retrieval
#'
#' @param x Passed argument 1
#' @param y Passed argument 2
#' @param TheData The datset
#'
#' @importFrom dplyr select
#' @importFrom tidyselect all_of
#' @importFrom dplyr slice
#'
#' @keywords internal
RetrievalMerge <- function(x, y, TheData){
  Individual <- TheData %>% select(all_of(c(x, y))) %>% slice(1)
  Cell <- Individual %>% pivot_wider(names_from = 1, values_from = 2)
  return(Cell)
}
