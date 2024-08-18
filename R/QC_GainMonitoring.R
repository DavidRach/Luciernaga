#' Takes Before After QC .fcs files of QC beads run on assay settings, returns the gains and average
#' measurement in a data.frame to allow for plotting to assess stability over time.
#'
#' @param x A CytoSet (when mapped) or an individual cytoset, example (x=MyCytoSet[[1]])
#' @param sample.name The keyword value that distinguishes individual .fcs files
#' @param stats Whether to summarize the MFIs by "mean" or "median"
#'
#' @importFrom flowCore exprs
#' @importFrom dplyr select
#' @importFrom stringr str_detect
#' @importFrom dplyr mutate
#' @importFrom dplyr relocate
#'
#' @return A data.frame row
#' @export
#'
#' @examples NULL
QC_GainMonitoring <- function(x, sample.name, stats){
  Guts <- Luciernaga::QC_Retrieval(x=x, sample.name=sample.name)
  Data <- data.frame(exprs(x), check.names=FALSE)
  Data <- Luciernaga:::AveragedSignature(Data, stats)
  Data <- Data %>% select(-Time)
  Bound <- cbind(Guts, Data)
  Bound[["SAMPLE"]] <- Luciernaga:::NameCleanUp(Bound[["SAMPLE"]], removestrings = ".fcs")
  Bound[["SAMPLE"]] <- Luciernaga:::NameCleanUp(Bound[["SAMPLE"]], removestrings = ".fcs")
  if (str_detect(Bound[["SAMPLE"]], "efore")){
    Bound <- Bound %>% mutate(Timepoint = "Before") %>% relocate(Timepoint, .after=TIME)
  } else if (str_detect(Bound[["SAMPLE"]], "fter")){
    Bound <- Bound %>% mutate(Timepoint = "After") %>% relocate(Timepoint, .after=TIME)
  } else {Bound <- Bound %>% mutate(Timepoint = "Unknown") %>% relocate(Timepoint, .after=TIME)}
  return(Bound)
}
