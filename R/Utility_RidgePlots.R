
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
#' @examples
#'
#' library(BiocGenerics)
#' library(flowCore)
#' library(flowWorkspace)
#' library(openCyto)
#' library(data.table)
#'
#' File_Location <- system.file("extdata", package = "Luciernaga")
#' FCS_Files <- list.files(path = File_Location, pattern = ".fcs",
#'   full.names = TRUE)
#' Unmixed_FullStained <- FCS_Files[grep("Unmixed", FCS_Files)]
#' UnmixedFCSFiles <- Unmixed_FullStained[1:2]
#' UnmixedCytoSet <- load_cytoset_from_fcs(UnmixedFCSFiles[1:2],
#'   truncate_max_range = FALSE,transformation = FALSE)
#' UnmixedGatingSet <- GatingSet(UnmixedCytoSet)
#' Markers <- colnames(UnmixedCytoSet)
#' KeptMarkers <- Markers[-grep("Time|FS|SC|SS|Original|-W$|-H$|AF", Markers)]
#' MyBiexponentialTransform <- flowjo_biexp_trans(channelRange = 256,
#'   maxValue = 1000000,pos = 4.5, neg = 0, widthBasis = -1000)
#' TransformList <- transformerList(KeptMarkers, MyBiexponentialTransform)
#' UnmixedGatingSet <- flowWorkspace::transform(UnmixedGatingSet, TransformList)
#' FileLocation <- system.file("extdata", package = "Luciernaga")
#' UnmixedGates <- fread(file.path(path = FileLocation, pattern = 'GatesUnmixed.csv'))
#' UnmixedGating <- gatingTemplate(UnmixedGates)
#' gt_gating(UnmixedGating, UnmixedGatingSet)
#'
#' removestrings <-  c("DTR_", ".fcs")
#' StorageLocation <- file.path("C:", "Users", "JohnDoe", "Desktop")
#'
#' Condition <- data.frame(Condition=c("Ctrl", "Ctrl"))
#' pd <- pData(UnmixedGatingSet)
#' new_pd <- cbind(pd, Condition)
#' pData(UnmixedGatingSet) <- new_pd
#'
#' SinglePlot <- Utility_RidgePlots(gs=UnmixedGatingSet, subset="live",
#'   TheFill="Condition",TheX = "APC-Fire 810-A", TheY="name",
#'   returntype="plots", outpath=StorageLocation, filename="RidgePlot_Condition")
#'
#'
Utility_RidgePlots <- function(gs, subset, TheX=NULL, TheY, TheFill, inverse.transform = FALSE,
                               outpath, returntype, therows=2,
                               thecolumns=1, width=7, height=9, filename){

  if (class(gs) == "GatingHierarchy"){
    cs <- gs_pop_get_data(gs, subset, inverse.transform = inverse.transform)
  } else {cs <- gs}

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
