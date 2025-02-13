#' Takes Before After QC .fcs files of QC beads run on assay settings, returns the gains and average
#' measurement in a data.frame to allow for plotting to assess stability over time.
#'
#' @param x A CytoSet (when mapped) or an individual cytoset, example (x=MyCytoSet[[1]])
#' @param sample.name The keyword value that distinguishes individual .fcs files
#' @param stats Whether to summarize the MFIs by "mean" or "median"
#' @param subsets When provided with a GatingSet, selects this subset to provide the Cytoset
#' @param inverse.transform Default is FALSE. 
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
#' @examples
#' library(dplyr)
#' library(purrr)
#' library(flowWorkspace)
#'
#' File_Location <- system.file("extdata", package = "Luciernaga")
#' FCS_Pattern <- ".fcs$"
#' FCS_Files <- list.files(path = File_Location, pattern = FCS_Pattern,
#' full.names = TRUE, recursive = FALSE)
#' QCBeads <- FCS_Files[grep("Before|After", FCS_Files)]
#' BeforeAfter_CS <- load_cytoset_from_fcs(files=QCBeads,
#' transformation=FALSE, truncate_max_range = FALSE)
#' BeforeAfter <- map(.x=BeforeAfter_CS[1:2], .f=QC_GainMonitoring,
#'  sample.name = "TUBENAME", stats="median") %>% bind_rows()
QC_GainMonitoring <- function(x, sample.name, stats, subsets=NULL,
inverse.transform=FALSE){

  if (class(x) == "GatingHierarchy"){
    cs <- gs_pop_get_data(x, subsets, inverse.transform=inverse.transform)
    x <- cs[[1]]
  }

  Guts <- QC_Retrieval(x=x, sample.name=sample.name)
  Data <- data.frame(exprs(x), check.names=FALSE)
  Data <- AveragedSignature(Data, stats)
  Data <- Data %>% select(-Time)
  Bound <- cbind(Guts, Data)
  Bound[["SAMPLE"]] <- NameCleanUp(Bound[["SAMPLE"]], removestrings = ".fcs")
  Bound[["SAMPLE"]] <- NameCleanUp(Bound[["SAMPLE"]], removestrings = ".fcs")
  if (str_detect(Bound[["SAMPLE"]], "efore")){
    Bound <- Bound %>% mutate(Timepoint = "Before") %>% relocate(Timepoint, .after=TIME)
  } else if (str_detect(Bound[["SAMPLE"]], "fter")){
    Bound <- Bound %>% mutate(Timepoint = "After") %>% relocate(Timepoint, .after=TIME)
  } else {Bound <- Bound %>% mutate(Timepoint = "Unknown") %>% relocate(Timepoint, .after=TIME)}
  return(Bound)
}
