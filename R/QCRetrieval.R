#' Returns the Detector and Laser settings from individual CytoSet files
#'
#' @param x A CytoSet object when mapped, or its frame (cs[1])
#' @param sample.name The keyword value that distinguishes individual .fcs files
#' @importFrom flowWorkspace keyword
#' @importFrom dplyr select
#' @importFrom dplyr bind_cols
#' @importFrom tidyr pivot_wider
#' @importFrom purrr map2
#'
#' @return A dataframe row
#' @export
#'
#' @examples Not at this time

QCRetrieval <- function(x, sample.name){
  KeywordsList <- keyword(x)
  KeywordsDF <- data.frame(KeywordsList, check.names = FALSE)
  TheColumnNames <- colnames(KeywordsDF)

  SAMPLE <- keyword(x)[[sample.name]]
  DATE <- keyword(x)$`$DATE`
  CYT <- keyword(x)$`$CYT`
  CYTSN <- keyword(x)$`$CYTSN`
  OP <- keyword(x)$`$OP`

  .Internal_Merge <- function(x, y, TheData){
    Individual <- TheData %>% select(all_of(c(x, y))) %>% slice(1)
    Cell <- Individual %>% pivot_wider(names_from = 1, values_from = 2)
    return(Cell)
  }

  PN_Names1 <- TheColumnNames[grepl("^\\$P[0-9]{1}N$", TheColumnNames)]
  PN_Names2 <- TheColumnNames[grepl("^\\$P[0-9]{2}N$", TheColumnNames)]
  PV_Gains1 <- TheColumnNames[grepl("^\\$P[0-9]{1}V$", TheColumnNames)]
  PV_Gains2 <- TheColumnNames[grepl("^\\$P[0-9]{2}V$", TheColumnNames)]
  PN_Names <- c(PN_Names1, PN_Names2)
  PN_Names <- PN_Names[-1] # Remove Time
  PV_Gains <- c(PV_Gains1, PV_Gains2)

  ParameterRows <- map2(.x=PN_Names, .y=PV_Gains, .f=.Internal_Merge, TheData=KeywordsDF) %>% bind_cols()

  Laser_Name <- TheColumnNames[grepl("^\\LASER[0-9]{1}NAME$", TheColumnNames)]
  Laser_Delay <- TheColumnNames[grepl("^\\LASER[0-9]{1}DELAY$", TheColumnNames)]
  Laser_ASF <- TheColumnNames[grepl("^\\LASER[0-9]{1}ASF$", TheColumnNames)]

  LaserDelayRows <- map2(.x=Laser_Name, .y=Laser_Delay, .f=.Internal_Merge, TheData=KeywordsDF) %>% bind_cols()
  colnames(LaserDelayRows) <- paste0(colnames(LaserDelayRows), "_LaserDelay")
  LaserASFRows <- map2(.x=Laser_Name, .y=Laser_ASF, .f=.Internal_Merge, TheData=KeywordsDF) %>% bind_cols()
  colnames(LaserASFRows) <- paste0(colnames(LaserASFRows), "_AreaScalingFactor")

  RecoveredQC <- cbind(SAMPLE, DATE, CYT, CYTSN, OP, ParameterRows, LaserDelayRows, LaserASFRows)

  return(RecoveredQC)
}
