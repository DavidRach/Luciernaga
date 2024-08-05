#' Compare the plots for given parameters for all specimens in the gating set.
#'
#' @param x The desired x-axis parameter
#' @param y The desired y-axis parameter
#' @param GatingSet The Gating Set that contains the specimens
#' @param marginsubset The desired gate for the margins
#' @param gatesubset The desired gate of cells you want to plot
#' @param sample.name The keyword that designates different specimens
#' @param removestrings Values to remove from both the plot titles and the pdf
#' @param clearance A value of clearance multiplied to the margin
#' @param bins How many bins to class the cells into
#' @param gatelines Whether to plot the reference lines
#' @param reference A data.frame containing references
#' @param pdf Whether to also output a pdf
#' @param outpath The desired location to send the assembled pdf to
#'
#' @importFrom Biobase pData
#' @importFrom purrr map
#' @importFrom patchwork wrap_plots
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.off
#' @importFrom stats quantile
#' @importFrom flowWorkspace keyword
#' @importFrom flowWorkspace gs_pop_get_data
#' @importFrom flowCore exprs
#' @importFrom dplyr select
#' @importFrom tidyr all_of
#' @importFrom ggcyto as.ggplot
#' @importFrom ggcyto ggcyto
#' @importFrom ggplot2 geom_hex
#'
#' @return The ggplots for all the specimens, as well as the optional .pdf
#' @export
#'
#' @examples  NULL

Utility_UnityPlot <- function(x, y, GatingSet, marginsubset, gatesubset, sample.name, removestrings,
                              clearance, bins, gatelines, reference, pdf, outpath){
  TheX <- x
  TheY <- y
  PlotNumber <- length(pData(GatingSet)[["name"]])

  FileName <- NameCleanUp(TheX, removestrings)
  StorageLocation <- file.path(outpath, FileName)

 Plots <- map(GatingSet, .f=InternalUnity, TheX=TheX, TheY=y, marginsubset=marginsubset, gatesubset=gatesubset,
      sample.name=sample.name, removestrings=removestrings, clearance=clearance, bins=bins, gatelines=gatelines,
      reference=reference)

 if (pdf == TRUE){
   AssembledPlots <- Utility_Patchwork(x=Plots, filename=AggregateName, outfolder=outpath, returntype = "pdf")
 } else {
   AssembledPlots <- Utility_Patchwork(x=Plots, filename=AggregateName, outfolder=outpath, returntype = "patchwork")
 }

 return(AssembledPlots)
}
