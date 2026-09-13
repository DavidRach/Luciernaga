#' Internal for AutofluorescenceShop, filters down AF signatures to show
#' main variants as plots for later compiling via Patchwork
#' 
#' @param x The main detector being iterated on as grouping variable
#'  for the plot signatures
#' @param data The data return from Luciernaga_QC
#' @param returnType Whether to return plots, alternate is plotly
#' @param legend Whether to show legend or not
#' @param plotname Name to append to the plot for specimen identification
#' 
#' @importFrom dplyr select filter pull
#' @importFrom stringr str_detect
#' @importFrom plotly ggplotly
#' 
#' @return ggplot or plotly object
#' 
#' @noRd
SignatureVariants <- function(x, data, returnType, legend, plotname) {
  x <- as.character(x)

  UnstainedSignature1 <- data |>
    select(-Sample, -Experiment, -Condition, -Count)

  These <- UnstainedSignature1 |>
    filter(str_detect(Cluster, paste0("^", x))) |>
    pull(Cluster) |>
    unique()

  Plots <- Luciernaga::QC_ViewSignature(x = These, columnname = "Cluster",
                                        data = UnstainedSignature1,
                                        Normalize = TRUE, TheFormat = "wider",
                                        legend = legend, plotname = plotname)

  if (returnType == "plots") {
    return(Plots)
  } else {
    Plots <- ggplotly(Plots)
    return(Plots)
  }
}