#' Calculate the RCV value for each marker for a specified GatingSet subset
#' 
#' @param x The Gating Set object being iterated
#' @param subset The Gating Set subset for population of interest
#' @param sample.name The keyword corresponding to the samples name
#' @param experiment.name Something
#' @param condition.name Something
#' 
#' @importFrom flowCore keyword
#' @importFrom flowWorkspace gs_pop_get_data
#' @importFrom Biobase exprs
#' @importFrom purrr map
#' @importFrom dplyr bind_cols
#' @importFrom dplyr mutate
#' @importFrom dplyr relocate
#' 
#' @return A data.frame row with the marker RCVs as a proportion. 
#' 
#' @export
RCVfromFCS <- function(x, subset, sample.name, experiment.name, condition.name){
  Name <- keyword(x, sample.name)
  Condition <- keyword(x, condition.name)
  Experiment <- keyword(x, experiment.name)
  Internal <- gs_pop_get_data(x, subset)
  Internal <- exprs(Internal[[1]])
  Internal <- data.frame(Internal, check.names=FALSE)
  These <- colnames(Internal)
  TheRCVs <- map(.x=These, .f=InternalRCV, data=Internal) |>
        bind_cols()
  TheRCVs <- TheRCVs |> mutate(Sample=Name, Condition=Condition,
     Experiment=Experiment) |> relocate(Sample, Condition, Experiment, .before=1)
  return(TheRCVs)
}

#' Internal for RCV from FCS
#' 
#' @param x The column name for an individual marker
#' @param data The data.frame of exprs values
#' 
#' @importFrom dplyr select
#' @importFrom stats mad
#' @importFrom stats median
#' 
#' @return An RCV for the individual marker
#' 
#' @noRd
InternalRCV <- function(x, data){
  Name <- x
  TheCol <- data |> select(x)
  TheCol <- TheCol[!is.na(TheCol)]
  if(length(TheCol) == 0) return(NULL)
  RCV <- mad(TheCol) / median(TheCol)
  RCV <- data.frame(RCV)
  colnames(RCV) <- Name
  return(RCV)
}
