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
#' @examples
#' library(flowWorkspace)
#'
#' File_Location <- system.file("extdata", package = "Luciernaga")
#' FCS_Pattern <- ".fcs$"
#' FCS_Files <- list.files(path = File_Location, pattern = FCS_Pattern,
#' full.names = TRUE, recursive = FALSE)
#' QCBeads <- FCS_Files[grep("Before", FCS_Files)]
#' MyCytoSet <- load_cytoset_from_fcs(QCBeads[1], truncate_max_range = FALSE,
#' transformation = FALSE)
#' MyGatingSet <- GatingSet(MyCytoSet)
#'
#' SingleSpecimen <- QC_Retrieval(x=MyGatingSet[[1]], sample.name="TUBENAME")
#'
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
  if(is.null(CYT)){CYTN <- "Unknown"}
  CYTSN <- keyword(x)$`$CYTSN`
  if(is.null(CYTSN)){CYTSN <- keyword(x)$`CYTNUM`}
  if(is.null(CYTSN)){CYTSN <- "Unknown"}
  OP <- keyword(x)$`$OP`
  if(is.null(OP)){OP <- "Unknown"}

  PN_Names1 <- TheColumnNames[grepl("^\\$P[0-9]{1}N$", TheColumnNames)]
  PN_Names2 <- TheColumnNames[grepl("^\\$P[0-9]{2}N$", TheColumnNames)]
  PN_Names3 <- TheColumnNames[grepl("^\\$P[0-9]{3}N$", TheColumnNames)]
  PV_Gains1 <- TheColumnNames[grepl("^\\$P[0-9]{1}V$", TheColumnNames)]
  PV_Gains2 <- TheColumnNames[grepl("^\\$P[0-9]{2}V$", TheColumnNames)]
  PV_Gains3 <- TheColumnNames[grepl("^\\$P[0-9]{3}V$", TheColumnNames)]

  if (length(PN_Names3) > 0){
    PN_Names <- c(PN_Names1, PN_Names2, PN_Names3)
  } else {  PN_Names <- c(PN_Names1, PN_Names2)}
  PN_Names <- PN_Names[-1] # Remove Time

  if (length(PV_Gains3) > 0){
    PV_Gains <- c(PV_Gains1, PV_Gains2, PV_Gains3)
  } else {PV_Gains <- c(PV_Gains1, PV_Gains2)}


  ParameterRows <- map2(.x=PN_Names, .y=PV_Gains, .f=Luciernaga:::RetrievalMerge,
                        TheData=KeywordsDF) %>% bind_cols()

  Laser_Name <- TheColumnNames[grepl("^\\LASER[0-9]{1}NAME$", TheColumnNames)]
  Laser_Delay <- TheColumnNames[grepl("^\\LASER[0-9]{1}DELAY$", TheColumnNames)]
  Laser_ASF <- TheColumnNames[grepl("^\\LASER[0-9]{1}ASF$", TheColumnNames)]

  LaserDelayRows <- map2(.x=Laser_Name, .y=Laser_Delay, .f=Luciernaga:::RetrievalMerge,
                         TheData=KeywordsDF) %>% bind_cols()
  colnames(LaserDelayRows) <- paste0(colnames(LaserDelayRows), "_LaserDelay")
  LaserASFRows <- map2(.x=Laser_Name, .y=Laser_ASF, .f=Luciernaga:::RetrievalMerge,
                       TheData=KeywordsDF) %>% bind_cols()
  colnames(LaserASFRows) <- paste0(colnames(LaserASFRows), "_AreaScalingFactor")

  ParameterRows <- ParameterRows %>% mutate(across(everything(), as.numeric))
  colnames(ParameterRows) <- paste0(colnames(ParameterRows), "_Gain")
  LaserDelayRows <- LaserDelayRows %>% mutate(across(everything(), as.numeric))
  if (ncol(LaserDelayRows) == 0) {LaserDelayRows <- NULL}
  LaserASFRows <- LaserASFRows %>% mutate(across(everything(), as.numeric))
  if (ncol(LaserASFRows) == 0) {LaserASFRows <- NULL}

  if(is.null(SAMPLE)){stop("sample.name keyword not recognized")}

  RecoveredQC <- cbind(SAMPLE, DATE, TIME, CYT, CYTSN, OP, ParameterRows)

  if(!is.null(LaserDelayRows)){RecoveredQC <- cbind(RecoveredQC, LaserDelayRows)}

  if(!is.null(LaserASFRows)){RecoveredQC <- cbind(RecoveredQC, LaserASFRows)}

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
#' @return An internal value
#'
#' @keywords internal
RetrievalMerge <- function(x, y, TheData){
  Individual <- TheData %>% select(all_of(c(x, y))) %>% slice(1)
  Cell <- Individual %>% pivot_wider(names_from = 1, values_from = 2)
  return(Cell)
}
