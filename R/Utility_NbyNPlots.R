#' Generate NxN plots for one or all fluorophores.
#'
#' @param x A GatingSet object (ex. gs or gs[[1]])
#' @param sample.name Keyword for which samples names are stored (ex. "GUID")
#' @param removestrings A string of characters to remove from name (
#' ex. c("DR_ILT_2023_", "Cord"))
#' @param experiment Provide directly an experiment name (ex. "Jan2024")
#' @param experiment.name Keyword for which experiment name is stored (
#' ex. "GROUPNAME")
#' @param outpath Location for which to store the generated .pdf
#' @param marginsubset The Gating Hierarchy level from which values will be used to
#'  estimate the plot margins (ex. "nonDebris")
#' @param gatesubset The Gating Hierarchy level of the cells that you want to see
#'  plotted (ex. "lymph")
#' @param ycolumn The ycolumn that you want to see everything plotted by (ex. "APC-A")
#'  or ALL to see all comparisons
#' @param bins Bins for which the plotted cells will be divided into providing
#'  granularity
#' @param clearance The additional ratio added to the margins to avoid clipping
#'  main population but exclude outliers.
#' @param pdf Prints default NxN plot, TRUE or FALSE.
#' @param condition Provide a condition name
#' @param condition.name The keyword in the .fcs file storing the condition.name
#' @param gatelines Whether to add estimated gate cutoff lines
#' @param reference Reference for the gate cutoff lines
#' @param width Desired page width for a pdf, default is 9 inches.
#' @param height Desired page height for a pdf, default is 7 inches
#'
#' @importFrom flowWorkspace keyword
#' @importFrom flowWorkspace gs_pop_get_data
#' @importFrom flowCore exprs
#' @importFrom patchwork wrap_plots
#' @importFrom patchwork plot_spacer
#' @importFrom purrr map
#'
#' @return A value to be determined later
#' @export
#'
#' @examples NULL

Utility_NbyNPlots <- function(x, sample.name, removestrings, experiment = NULL,
  experiment.name = NULL, condition = NULL, condition.name = NULL, marginsubset,
  gatesubset, ycolumn, bins, clearance, gatelines, reference = NULL, outpath, pdf,
  width = 9, height = 7, ...) {
  #ycolumn <- ycolumn
  #x <- x

  AggregateName <- NameForSample(
    x=x, sample.name=sample.name, removestrings=removestrings, ...)

  StorageLocation <- file.path(outpath, AggregateName)

  mff <- gs_pop_get_data(x, marginsubset)
  df <- exprs(mff[[1]])
  TheDF <- data.frame(df, check.names = FALSE)
  DFNames <- colnames(TheDF[,-grep("Time|FS|SC|SS|Original|W$|H$", names(TheDF))])
  PlotNumber <- length(DFNames)

  ff <- gs_pop_get_data(x, gatesubset)

  if (ycolumn == "ALL"){

    Plots <- map(.x=DFNames, .f = UniversalIterator, x_ff=ff,
                 TheDF=TheDF, yValue=ycolumn, columnlist=DFNames, gatelines=gatelines,
                 reference=reference, clearance=clearance, bins=bins)

    Plots <- flatten(Plots)

  } else {
    columnlist <- DFNames[DFNames != ycolumn]
    Plots <- map(.x = columnlist, .f = GeneralGating, name = name, ff = ff,
      yValue = ycolumn, columnlist = DFNames, TheDF = TheDF, gatelines = gatelines,
      reference = reference, clearance, bins)
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

#' Internal for Utility_NbyNPlots
#' @importFrom flowWorkspace keyword
#' @importFrom flowWorkspace gs_pop_get_data
#' @importFrom flowCore exprs
#' @importFrom patchwork wrap_plots
#' @importFrom patchwork plot_spacer
#' @importFrom purrr map
#'
#' @noRd
UniversalIterator <- function(x, x_ff,
                              TheDF, yValue, columnlist, gatelines,
                              reference, clearance, bins, AltNameX,
                              AltNameY, colorX, colorY){
  ff <- x_ff
  DFNames <- columnlist
  columnlist <- columnlist[columnlist != x] # Remove the universal Y value

  Plots <- map(.x = columnlist, .f = GeneralGating, name = name, ff = ff,
    yValue = x, columnlist = DFNames, TheDF = TheDF, gatelines = gatelines,
    reference = reference, clearance=clearance, bins=bins)

  #Plots <- flatten(Plots)
  #Plots1 <- Plots
}

#' Internal for Utility_ParallelNbyNPlots
#' @importFrom flowWorkspace keyword
#' @importFrom flowWorkspace gs_pop_get_data
#' @importFrom flowCore exprs
#' @importFrom patchwork wrap_plots
#' @importFrom patchwork plot_spacer
#' @importFrom purrr map
#'
#' @noRd
ParallelUniversalIterator <- function(x, x_ff, y_ff,
                                      TheDF, yValue, columnlist, gatelines,
                                      reference, clearance, bins, AltNameX,
                                      AltNameY, colorX, colorY){
  DFNames <- columnlist
  columnlist <- columnlist[columnlist != x] # Remove the universal Y value

  Plots <- map(.x = columnlist, .f = ParallelGating, x_ff=x_ff, y_ff=y_ff,
               TheDF=TheDF, yValue=x, columnlist=DFNames, gatelines=gatelines,
               reference=reference, clearance=clearance, bins=bins, AltNameX=AltNameX,
               AltNameY=AltNameY, colorX=colorX, colorY=colorY) #Name
  Plots <- flatten(Plots)

  #Plots1 <- Plots
  #Plots <- flatten(Plots)
  return(Plots)
}

#' An internal function that generates ggplots for given channels
#'
#' @param x Passed channel
#' @param name The name that matches to the data.frame and gs object,
#' ex. altered.name
#' @param ff The dataframe data for that sample
#' @param yValue What wanted on the yaxis
#' @param bins Bins for which the plotted cells will be divided into
#'  providing granularity
#' @param clearance Space around the plot
#' @param columnlist list all channels with x removed
#' @param TheDF the external limits settings
#' @param gatelines whether to plot the ModernCutoffLines from reference
#'  dataframe
#' @param reference location ModernCutoff dataframe.
#'
#' @importFrom dplyr select
#' @importFrom dplyr pull
#' @importFrom flowCore keyword
#' @importFrom ggcyto ggcyto
#' @importFrom ggcyto as.ggplot
#' @import ggplot2
#'
#' @return A value to be determined later
#' @noRd
GeneralGating <- function(x, name, ff, yValue, clearance, bins,
                          columnlist, TheDF, gatelines, reference = NULL) {

  if (yValue == x){stop("x equals yValue and can't be plotted")}

  xValue <- x

  if (!grepl("FSC|SSC", yValue)) {ExprsData <- TheDF %>%
    select(all_of(yValue)) %>% pull()
  theYmin <- ExprsData %>% quantile(., 0.001)
  theYmax <- ExprsData %>% quantile(., 0.999)
  theYmin <- theYmin - abs((clearance*theYmin))
  theYmax <- theYmax + (clearance*theYmax)}

  if (!grepl("FSC|SSC", xValue)) {ExprsData <- TheDF %>%
    select(all_of(xValue)) %>% pull()
  theXmin <- ExprsData %>% quantile(., 0.001)
  theXmax <- ExprsData %>% quantile(., 0.999)
  theXmin <- theXmin - abs((clearance*theXmin))
  theXmax <- theXmax + (clearance*theXmax)}


  if (!exists("theYmax") || !exists("theXmax")){
    Plot <- as.ggplot(ggcyto(ff, aes(x = .data[[xValue]], y = .data[[yValue]]),
     subset = "root") + geom_hex(bins=bins) + theme_bw() + labs(title = NULL) +
     theme(strip.background = element_blank(), strip.text.x = element_blank(),
     panel.grid.major = element_line(linetype = "blank"),
     panel.grid.minor = element_line(linetype = "blank"),
     axis.title = element_text(size = 10, face = "bold"),
     legend.position = "none"))

    if (gatelines == TRUE){Value <- reference[reference$specimen == name, xValue]
    Plot <- Plot + geom_vline(xintercept = c(seq(0,200,25)), colour = "gray") +
      geom_vline(xintercept = Value, colour = "red")}

  } else {Plot <- as.ggplot(ggcyto(ff, aes(x = .data[[xValue]], y = .data[[yValue]]),
          subset = "root") + geom_hex(bins=bins) + coord_cartesian(
          xlim = c(theXmin, theXmax), ylim = c(theYmin, theYmax), default = TRUE) +
          theme_bw() + labs(title = NULL) + theme(strip.background = element_blank(),
          strip.text.x = element_blank(), panel.grid.major = element_line(
          linetype = "blank"), panel.grid.minor = element_line(linetype = "blank"),
          axis.title = element_text(size = 10, face = "bold"),
          legend.position = "none"))

  if (gatelines == TRUE){Value <- reference[reference$specimen == name, xValue]
  Plot <- Plot + geom_vline(xintercept = c(seq(0,200,25)), colour = "gray") +
    geom_vline(xintercept = Value, colour = "red")}
  }

  return(Plot)
}
