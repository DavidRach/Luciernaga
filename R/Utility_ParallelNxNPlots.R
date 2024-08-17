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
#' @importFrom flowWorkspace keyword
#' @importFrom flowWorkspace gs_pop_get_data
#' @importFrom flowCore exprs
#' @importFrom patchwork wrap_plots
#' @importFrom patchwork plot_spacer
#' @importFrom purrr map
#' @importFrom dplyr select
#' @importFrom dplyr pull
#' @importFrom dplyr mutate
#' @importFrom ggplot2 ggplot
#'
#'
#' @return Either list ggplot objects or a pdf object
#' @export
#'
#' @examples NULL

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

#' Internal for Utility_ParallelNxNPlots
#'
#' @importFrom flowWorkspace keyword
#' @importFrom flowWorkspace gs_pop_get_data
#' @importFrom flowCore exprs
#' @importFrom patchwork wrap_plots
#' @importFrom patchwork plot_spacer
#' @importFrom purrr map
#' @importFrom dplyr select
#' @importFrom dplyr pull
#' @importFrom dplyr mutate
#' @importFrom ggplot2 ggplot
#'
#' @noRd
ParallelGating <- function(x, x_ff, y_ff, TheDF, yValue, columnlist, gatelines,
  reference, clearance, bins, AltNameX, AltNameY, colorX, colorY) {

  if (yValue == x){stop("x equals yValue and can't be plotted")}

  xValue <- x

  if (!grepl("FSC|SSC", yValue)) {
    ExprsData <- TheDF %>% select(all_of(yValue)) %>% pull()
    theYmin <- ExprsData %>% quantile(., 0.001)
    theYmax <- ExprsData %>% quantile(., 0.999)
    theYmin <- theYmin - abs((clearance*theYmin))
    theYmax <- theYmax + (clearance*theYmax)}

  if (!grepl("FSC|SSC", xValue)) {
    ExprsData <- TheDF %>% select(all_of(xValue)) %>% pull()
    theXmin <- ExprsData %>% quantile(., 0.001)
    theXmax <- ExprsData %>% quantile(., 0.999)
    theXmin <- theXmin - abs((clearance*theXmin))
    theXmax <- theXmax + (clearance*theXmax)}


  if (!exists("theYmax") || !exists("theXmax")){
    stop("Either theYmax or theXmax didn't exist, and since I didn't think
     it relavant to duplicate this code in the parallel NxN plot when coding,
          the function now crashed ")
  } else {

    x_ffXdata <- exprs(x_ff[[1]]) %>% data.frame(check.names = FALSE) %>%
      select(all_of(xValue))
    x_ffYdata <- exprs(x_ff[[1]]) %>% data.frame(check.names = FALSE) %>%
      select(all_of(yValue))
    Thex_ff <- cbind(x_ffXdata, x_ffYdata) %>% mutate(specimen = AltNameX)

    y_ffXdata <- exprs(y_ff[[1]]) %>% data.frame(check.names = FALSE) %>%
      select(all_of(xValue))
    y_ffYdata <- exprs(y_ff[[1]]) %>% data.frame(check.names = FALSE) %>%
      select(all_of(yValue))
    They_ff <- cbind(y_ffXdata, y_ffYdata) %>% mutate(specimen = AltNameY)

    TheData <- rbind(Thex_ff, They_ff)
    TheData$specimen <- as.factor(TheData$specimen)

    specimen_counts <- table(TheData$specimen)
    sorted_specimens <- names(sort(desc(specimen_counts)))
    TheData$specimen <- factor(TheData$specimen, levels = sorted_specimens)

    # Attempted Work Around for Specifying Colors while adjusting what
    #population is displayed forward.
    Xscheme <- cbind(AltNameX, colorX)
    Yscheme <- cbind(AltNameY, colorY)
    ColorFrame <- rbind(Xscheme, Yscheme)
    ColorFrame <- data.frame(ColorFrame, check.names = FALSE) %>%
      rename(specimen = AltNameX)
    ColorFrame$specimen <- factor(ColorFrame$specimen, levels = sorted_specimens)
    ColorFrame <- ColorFrame[order(ColorFrame$specimen), ]
    color1 <- ColorFrame[1,2]
    color2 <- ColorFrame[2,2]

    Plot <- ggplot(TheData, aes(x=.data[[xValue]], y = .data[[yValue]],
      fill = specimen)) + geom_hex(bins=bins, alpha = 0.5) + scale_fill_manual(
      values = c(color1, color2)) + coord_cartesian(xlim = c(
      theXmin, theXmax), ylim = c(theYmin, theYmax)) + theme_bw() + labs(
      title = NULL) + theme(strip.background = element_blank(),
      strip.text.x = element_blank(), panel.grid.major = element_line(
      linetype = "blank"), panel.grid.minor = element_line(linetype = "blank"),
      axis.title = element_text(size = 10, face = "bold"), legend.position = "none")

    #if (gatelines == TRUE){Value <- reference[reference$specimen == name, xValue]
    #Plot <- Plot + geom_vline(xintercept = c(seq(0,200,25)), colour = "gray") +
    #  geom_vline(xintercept = Value, colour = "red")}
  }

  return(Plot)
}
