#' Summarize a data.frame to desired stat
#'
#' @param x A data.frame containing double or numeric data.
#' @param stats Desired Stats "mean" or "median" to pass to summarize_all
#' @param normalize Default FALSE, TRUE peak detector normalizes.
#'
#' @importFrom dplyr summarize_all
#'
#' @return A data.frame row of summarized data
#' @export
#'
#' @examples
#'
#' library(flowCore)
#' library(flowWorkspace)
#' library(openCyto)
#' library(data.table)
#'
#' File_Location <- system.file("extdata", package = "Luciernaga")
#' FCS_Files <- list.files(path = File_Location, pattern = ".fcs",
#'   full.names = TRUE)
#' UnstainedFCSFiles <- FCS_Files[grep("Unstained", FCS_Files)]
#' UnstainedCells <- UnstainedFCSFiles[-grep("Beads", UnstainedFCSFiles)]
#' MyCytoSet <- load_cytoset_from_fcs(UnstainedCells[1],
#'   truncate_max_range = FALSE,transformation = FALSE)
#' MyGatingSet <- GatingSet(MyCytoSet)
#' MyGates <- fread(file.path(path = File_Location, pattern = 'Gates.csv'))
#' MyGatingTemplate <- gatingTemplate(MyGates)
#' gt_gating(MyGatingTemplate, MyGatingSet)
#' removestrings <-  c("DR_", "Cells", ".fcs", "-", " ")
#' StorageLocation <- file.path("C:", "Users", "JohnDoe", "Desktop")
#'
#' PopulationInterest <- gs_pop_get_data(MyGatingSet[1], subset="lymphocytes")
#' TheDataValues <- exprs(PopulationInterest[[1]])
#' TheDataValues <- data.frame(TheDataValues, check.names=FALSE)
#'
#' Signature <- AveragedSignature(TheDataValues, stats="median")
#'
AveragedSignature <- function(x, stats, normalize=FALSE){

  if (normalize == TRUE){
    x[x < 0] <- 0
    A <- do.call(pmax, x)
    x <- x/A
  }

  Signature <- x |> summarize_all(stats)

  if (normalize == TRUE){
  Signature <- round(Signature, 3)
  }

  return(Signature)
}