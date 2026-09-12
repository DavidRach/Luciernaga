#' Overlay plot two .fcs files from two different gs.
#'
#' @param x The first gs object
#' @param y The second gs object
#' @param sample.name The keyword under which sample names are stored
#' @param removestrings A list of values to remove from sample names
#' @param Override Exclude raw columns with -A$
#' @param marginsubset The gs subset that defines the plot margin
#' @param gatesubset The gs subset of interest
#' @param ycolumn The desired y-column for the comparisons
#' @param bins Desired number of hex bins
#' @param clearance A multiplication factor for margin wiggle room (0.2)
#' @param colorX Color for the x gs
#' @param colorY Color for the y gs
#' @param gatelines Whether to plot .csv specified gate lines
#' @param reference Reference for .csv specified gate lines
#' @param outpath Location which to store the output
#' @param pdf Whether to return as a pdf
#'
#' @importFrom flowWorkspace keyword gs_pop_get_data
#' @importFrom flowCore exprs
#' @importFrom patchwork wrap_plots plot_spacer
#' @importFrom purrr map
#' @importFrom dplyr select pull mutate
#' @importFrom ggplot2 ggplot
#'
#'
#' @return Either list ggplot objects or a pdf object
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
#' Plot <- Utility_ParallelNbyNPlots(x=UnmixedGatingSet[1], y = UnmixedGatingSet[2],
#'  sample.name = "GROUPNAME", removestrings = ".fcs", Override = FALSE,
#'  marginsubset = "lymphocytes", gatesubset = "live", ycolumn = "Spark Blue 550-A",
#'   bins = 120, clearance = 0.2, colorX = "lightblue", colorY = "orange", gatelines = FALSE,
#'   reference = NULL, outpath = StorageLocation,pdf = FALSE)
#'
Utility_ParallelNbyNPlots <- function(x, y, sample.name, removestrings,
    Override = FALSE, marginsubset, gatesubset, ycolumn, bins, clearance,
    colorX, colorY, gatelines, reference = NULL, outpath, pdf){

  nameX <- keyword(x, sample.name)
  AltNameX <- NameCleanUp(name = nameX, removestrings)

  nameY <- keyword(y, sample.name)
  AltNameY <- NameCleanUp(name = nameY, removestrings)

  if (AltNameX == AltNameY){message("X and Y names match, attaching numbers")
                           AltNameX <- paste0(AltNameX, "_1")
                           AltNameY <- paste0(AltNameY, "_2")}

  AlternateName <- paste(AltNameX, AltNameY, sep="_")
  AlternateName <- gsub("-", "", gsub(" ", "", AlternateName))

  #Would need to modify NameFor to handle multiple names if incorporating here.

  # Name and Location of Final PDF
  StorageLocation <- file.path(outpath, AlternateName)

  # Retrieving margin info for the x specimen
  xMargin <- gs_pop_get_data(x, marginsubset)
  xdf <- exprs(xMargin[[1]])
  TheXDF <- data.frame(xdf, check.names = FALSE)

  if(Override == TRUE){
    X_DFNames <- colnames(
      TheXDF[,-grep("Time|FS|SC|SS|Original|W$|H$|-A$", names(TheXDF))])
  } else {X_DFNames <- colnames(
    TheXDF[,-grep("Time|FS|SC|SS|Original|W$|H$", names(TheXDF))])
  }

  X_PlotNumber <- length(X_DFNames)

  # Retrieving margin info for the y specimen
  yMargin <- gs_pop_get_data(y, marginsubset)
  ydf <- exprs(yMargin[[1]])
  TheYDF <- data.frame(ydf, check.names = FALSE)

  if(Override == TRUE){
  Y_DFNames <- colnames(
    TheYDF[,-grep("Time|FS|SC|SS|Original|W$|H$|-A$", names(TheYDF))])
  } else {Y_DFNames <- colnames(
    TheYDF[,-grep("Time|FS|SC|SS|Original|W$|H$", names(TheYDF))])
  }

  Y_PlotNumber <- length(Y_DFNames)

  if (all(X_DFNames == Y_DFNames)){DFNames <- X_DFNames
                                   PlotNumber <- X_PlotNumber
  } else (stop("The two fcs files do not have matching column parameters"))

  TheDF <- rbind(TheXDF, TheYDF) # Merge the two margin data frames

  # Retrieving gating info for the x and y specimens
  x_ff <- gs_pop_get_data(x, gatesubset) #Sends cytoset forward
  y_ff <- gs_pop_get_data(y, gatesubset)

  if (ycolumn == "ALL"){

    Plots <- map(.x=DFNames, .f = ParallelUniversalIterator, x_ff=x_ff, y_ff=y_ff,
       TheDF=TheDF, yValue=ycolumn, columnlist=DFNames, gatelines=gatelines,
       reference=reference, clearance=clearance, bins=bins, AltNameX=AltNameX,
       AltNameY=AltNameY, colorX=colorX, colorY=colorY)
  } else {
    columnlist <- DFNames[DFNames != ycolumn] # Remove the universal Y value

    Plots <- map(.x = columnlist, .f = ParallelGating, x_ff=x_ff, y_ff=y_ff,
         TheDF=TheDF, yValue=ycolumn, columnlist=DFNames, gatelines=gatelines,
         reference=reference, clearance=clearance, bins=bins, AltNameX=AltNameX,
         AltNameY=AltNameY, colorX=colorX, colorY=colorY) #Name
    }

  if (pdf == TRUE){
    AssembledPlots <- Utility_Patchwork(x=Plots, filename=AggregateName,
                                        outfolder=outpath, returntype = "pdf")
  } else {
    AssembledPlots <- Utility_Patchwork(x=Plots, filename=AggregateName,
                                        outfolder=outpath, returntype = "patchwork")
  }

  return(AssembledPlots)
}

