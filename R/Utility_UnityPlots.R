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
#' @param reference Reference for the gate cutoff lines
#' @param returntype Whether to return "pdf", "patchwork" or "plots"
#' @param outpath The desired location to send the assembled pdf to
#' @param filename Default NULL, provide name to set the filename.
#' @param cartesian Default TRUE, set to false to remove cartesian_coord centering
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
#'SingleUnityPlot <- Utility_UnityPlot(x="Spark Blue 550-A", y="BUV805-A",
#' GatingSet=UnmixedGatingSet, sample.name="GROUPNAME", bins=100, clearance=0.2,
#' removestrings=removestrings, marginsubset="lymphocytes", gatesubset="live",
#' gatelines=FALSE, reference=NULL, returntype="patchwork",outpath=StorageLocation)
#'
Utility_UnityPlot <- function(x, y, GatingSet, marginsubset, gatesubset,
  sample.name, removestrings, clearance, bins, gatelines, reference, returntype,
  outpath, filename=NULL, cartesian=TRUE){
  
  if (!is.null(reference)){
      if (is.data.frame(reference)){reference <- reference
      } else {reference <- read.csv(file=reference, check.names = FALSE)
      }
  }  

  TheX <- x
  TheY <- y
  FileName <- NameCleanUp(TheX, removestrings)

  #x <- GatingSet[[1]]
 Plots <- map(GatingSet, .f=Unity, TheX=TheX, TheY=TheY, marginsubset=marginsubset,
      gatesubset=gatesubset, sample.name=sample.name, removestrings=removestrings,
      clearance=clearance, bins=bins, gatelines=gatelines, reference=reference,
      cartesian=cartesian)

 if (!is.null(filename)){FileName <- filename} 
  
 if (returntype == "plots"){
   return(Plots)
 }

 if (returntype == "pdf"){
   AssembledPlots <- Utility_Patchwork(x=Plots, filename=FileName,
                                       outfolder=outpath, returntype = "pdf")
 }

 if (returntype == "patchwork"){
   AssembledPlots <- Utility_Patchwork(x=Plots, filename=FileName,
                                       outfolder=outpath, returntype = "patchwork")
 }
}

#' Internal to Utility_UnityPlots
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
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_hex
#' @importFrom ggplot2 coord_cartesian
#' @importFrom ggplot2 element_line
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_line
#' @importFrom ggplot2 geom_vline
#'
#' @return An internal value
#'
#' @noRd
Unity <- function(x, TheY, TheX, marginsubset, gatesubset, sample.name, removestrings,
                  clearance, bins, gatelines, reference, cartesian=TRUE){

  name <- keyword(x, sample.name)
  name <- NameCleanUp(name = name, removestrings)
  
  if (inherits(x, "flowFrame")){
    message("flowFrame detected, unable to use gatesubset, using root instead")
    df <- exprs(x)
  } else {
    ff <- gs_pop_get_data(x, marginsubset)
    df <- exprs(ff[[1]])
  }

  TheDF <- data.frame(df, check.names = FALSE)

  if (!TheX == TheY) {
    YExprsData <- TheDF |> select(all_of(TheY)) |> pull()
    theYmin <- YExprsData %>% quantile(., 0.001)
    theYmax <- YExprsData %>% quantile(., 0.999)
    theYmin <- theYmin - abs((clearance*theYmin))
    theYmax <- theYmax + (clearance*theYmax)

    XExprsData <- TheDF |> select(all_of(TheX)) |> pull()
    theXmin <- XExprsData %>% quantile(., 0.001)
    theXmax <- XExprsData %>% quantile(., 0.999)
    theXmin <- theXmin - abs((clearance*theXmin))
    theXmax <- theXmax + (clearance*theXmax)
  } else (stop("TheX and TheY have the same value"))

  if (inherits(x, "flowFrame")){
    ff1 <- x
  } else {ff1 <- gs_pop_get_data(x, gatesubset)}

  if (BiocGenerics::nrow(ff1) < 200) {

    if (cartesian == TRUE){
    Plot <- as.ggplot(ggcyto(ff1, aes(x = .data[[TheX]], y = .data[[TheY]]),
     subset = "root") + geom_point(size = 2, alpha = 0.8) + coord_cartesian(
       xlim = c(theXmin, theXmax), ylim = c(theYmin, theYmax), default = TRUE) +
       theme_bw() + labs(title = name) + theme(strip.background = element_blank(),
     strip.text.x = element_blank(), panel.grid.major = element_line(
     linetype = "blank"), panel.grid.minor = element_line(linetype = "blank"),
     axis.title = element_text(size = 10, face = "bold"),legend.position = "none"))
     } else {
      Plot <- as.ggplot(ggcyto(ff1, aes(x = .data[[TheX]], y = .data[[TheY]]),
     subset = "root") + geom_point(size = 2, alpha = 0.8) +
       theme_bw() + labs(title = name) + theme(strip.background = element_blank(),
     strip.text.x = element_blank(), panel.grid.major = element_line(
     linetype = "blank"), panel.grid.minor = element_line(linetype = "blank"),
     axis.title = element_text(size = 10, face = "bold"),legend.position = "none"))
     }
      
    } else {
    
    if (cartesian == TRUE){
    Plot <- as.ggplot(ggcyto(ff1, aes(x = .data[[TheX]], y = .data[[TheY]]),
     subset = "root") + geom_hex(bins=bins) + coord_cartesian(
       xlim = c(theXmin, theXmax), ylim = c(theYmin, theYmax), default = TRUE) +
       theme_bw() + labs(title = name) +
       theme(strip.background = element_blank(), strip.text.x = element_blank(),
     panel.grid.major = element_line(linetype = "blank"),
     panel.grid.minor = element_line(linetype = "blank"),
     axis.title = element_text(size = 10, face = "bold"),
     legend.position = "none"))
     } else {
      Plot <- as.ggplot(ggcyto(ff1, aes(x = .data[[TheX]], y = .data[[TheY]]),
     subset = "root") + geom_hex(bins=bins) +
       theme_bw() + labs(title = name) +
       theme(strip.background = element_blank(), strip.text.x = element_blank(),
     panel.grid.major = element_line(linetype = "blank"),
     panel.grid.minor = element_line(linetype = "blank"),
     axis.title = element_text(size = 10, face = "bold"),
     legend.position = "none"))
     }
    }

  if (gatelines == TRUE){Value <- reference %>% dplyr::filter(specimen %in% name) %>%
    select(all_of(TheX)) %>% pull(.)
  Plot <- Plot + geom_vline(xintercept = c(seq(0,200,25)), colour = "gray") +
    geom_vline(xintercept = c(seq(0,200,2)), colour = "white", alpha = 0.1) +
    geom_vline(xintercept = Value, colour = "red")}

  return(Plot)
}

