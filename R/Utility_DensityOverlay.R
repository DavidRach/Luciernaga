#' Generates an overlay plot containing all objects in the GatingSet
#'
#' @param gs A GatingSet object
#' @param subset The desired Gating Hierarchy node (ex. "lymphocytes")
#' @param inverse.transform Default is FALSE, reverses transformation applied to
#' GatingSet as its converted to a CytoSet
#' @param TheX A desired marker to plot, leave NULL for all markers
#' @param TheFill A desired marker to color individual specimens by (as named in pData)
#' @param outpath Desired storage location
#' @param returntype Whether to return a "pdf", "patchwork" or "plots"
#' @param therows The desired number of rows per page
#' @param thecolumns The desired number of columns per page
#' @param width The desired page width
#' @param height The desired page height
#' @param filename The file name for the new .pdf
#'
#' @importFrom flowWorkspace gs_pop_get_data
#' @importFrom BiocGenerics colnames
#' @importFrom purrr map
#'
#' @return A pdf file or ggplot object
#' @export
#'
#' @examples NULL
Utility_DensityOverlay <- function(gs, subset, inverse.transform = FALSE, TheX=NULL, TheFill,
                                   returntype, outpath, filename, therows=3,
                                   thecolumns=3, width = 7, height=9){

  cs <- gs_pop_get_data(gs, subset, inverse.transform = inverse.transform)
  markers <- colnames(cs)
  DFNames <- markers[-grep("Time|FS|SC|SS|Original|-W$|-H$|AF", markers)]

  if(!is.null(TheX)){DFNames <- TheX}

  ThePlots <- map(x=DFNames, .f=InternalDensity, cs=cs, TheFill=TheFill)

  if (returntype == "pdf"){
    AssembledPlots <- Utility_Patchwork(x=ThePlots, filename=filename, outfolder=outpath,
    returntype = "pdf", therows=therows, thecolumns=thecolumns, width = width, height=height)
  }

  if (returntype == "patchwork"){
    AssembledPlots <- Utility_Patchwork(x=Plots, filename=filename, outfolder=outpath,
    returntype = "patchwork", therows=therows, thecolumns=thecolumns, width = width,
    height=height)
  }

  if (returntype == "plots"){
    AssembledPlots <- Utility_Patchwork(x=Plots, filename=filename, outfolder=outpath,
    returntype = "plots", therows=therows, thecolumns=thecolumns, width = width,
    height=height)
  }
}


#' Internal for Utility_DensityOverlay
#'
#' @param x A passed x-axis name
#' @param cs The passed CytoSet
#' @param TheFill The passed factor to fill plots by
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_density
#'
#' @keywords internal
InternalDensity <- function(x, cs, TheFill){
  Plot <- ggplot(cs, aes(x = .data[[x]], fill = .data[[TheFill]])) +
    geom_density(alpha = 0.2) + coord_cartesian(xlim = c(0, NA)) + theme_bw() + theme(
      legend.position = "none")
  return(Plot)
}
