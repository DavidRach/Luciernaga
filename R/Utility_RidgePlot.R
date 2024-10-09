
#' Draw Ridge Plots for comparison
#'
#' @param gs A Gating Set Object
#' @param subset The desired gate node, ie "lymphocytes"
#' @param TheX Optional, you can pass a list of column markers to plot.
#' @param TheY The desired Y-axis parameter (as named in pData)
#' @param TheFill The desired fill color (as named in pData)
#' @param inverse.transform Default is FALSE, reverts any transformation on
#' going from the GatingSet to the CytoSet object
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
#' @return A pdf file
#' @export
#'
#' @examples NULL
Utility_RidgePlots <- function(gs, subset, TheX=NULL, TheY, TheFill, inverse.transform = FALSE,
                               outpath, returntype, therows=2,
                               thecolumns=1, width=7, height=9, filename){

  cs <- gs_pop_get_data(gs, subset, inverse.transform = inverse.transform)
  markers <- colnames(cs)
  DFNames <- markers[-grep("Time|FS|SC|SS|Original|-W$|-H$|AF", markers)]

  if(!is.null(TheX)){DFNames <- TheX}

  ThePlots <- map(.x=DFNames, .f=RidgePlots, cs=cs, TheY=TheY, TheFill=TheFill)

  if (returntype == "pdf"){
    AssembledPlots <- Utility_Patchwork(x=ThePlots, filename=filename, outfolder=outpath,
       returntype = "pdf", therows=therows, thecolumns=thecolumns, width = width,
       height=height)

  } else if (returntype == "patchwork"){
    AssembledPlots <- Utility_Patchwork(x=ThePlots, filename=filename, outfolder=outpath,
        returntype = "patchwork", therows=therows, thecolumns=thecolumns, width = width,
        height=height)
    return(AssembledPlots)
  } else if (returntype == "plots"){return(ThePlots)}
}



#' Internal to Utility_RidgePlots
#'
#' @param x Passed Argument for X axis
#' @param cs The CytoSet with the data
#' @param TheY The passed argument for the Y axis
#' @param TheFill The passed argument for the Fill color by factor.
#'
#' @importFrom ggcyto ggcyto
#' @importFrom ggcyto axis_x_inverse_trans
#' @importFrom ggcyto as.ggplot
#' @importFrom ggridges geom_density_ridges
#' @importFrom ggplot2 facet_null
#'
#' @return An internal value
#'
#' @noRd
RidgePlots <- function(x, cs, TheY, TheFill){
  p <- ggcyto(cs, aes(x = .data[[x]]))
  p1 <- p + geom_density_ridges(aes(y = .data[[TheY]], fill = .data[[TheFill]], alpha = 0.2)) +
    facet_null() + theme_bw() + coord_cartesian(xlim = c(0, NA)) + theme(legend.position = "none") +
    axis_x_inverse_trans()
  p1 <- as.ggplot(p1)
  return(p1)
}
