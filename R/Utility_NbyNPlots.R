#' Generate NxN plots for one or all fluorophores.
#'
#' @param x A GatingSet object (ex. gs or gs[[1]])
#' @param sample.name Keyword for which samples names are stored (ex. "GUID")
#' @param removestrings A string of characters to remove from name (ex. c("DR_ILT_2023_", "Cord"))
#' @param experiment Provide directly an experiment name (ex. "Jan2024")
#' @param experiment.name Keyword for which experiment name is stored (ex. "GROUPNAME")
#' @param outpath Location for which to store the generated .pdf
#' @param marginsubset The Gating Hierarchy level from which values will be used to estimate the plot margins (ex. "nonDebris")
#' @param gatesubset The Gating Hierarchy level of the cells that you want to see plotted (ex. "lymph")
#' @param ycolumn The ycolumn that you want to see everything plotted by (ex. "APC-A") or ALL to see all comparisons
#' @param bins Bins for which the plotted cells will be divided into providing granularity
#' @param clearance The additional ratio added to the margins to avoid clipping main population but exclude outliers.
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

Utility_NbyNPlots <- function(x, sample.name, removestrings, experiment = NULL, experiment.name = NULL,
                              condition = NULL, condition.name = NULL, marginsubset, gatesubset,
                              ycolumn, bins, clearance, gatelines, reference = NULL, outpath, pdf,
                              width = 9, height = 7) {
  #ycolumn <- ycolumn
  #x <- x

  name <- keyword(x, sample.name)
  name <- NameCleanUp(name = name, removestrings)

  if(!is.null(experiment)){experiment <- experiment
  } else if (exists("experiment.name")) {experiment <- keyword(x, experiment.name)
                                         experiment <- NameCleanUp(name=experiment, removestrings)}

  if(exists("condition")){condition <- condition
  } else if (exists("condition.name")){condition <- keyword(x, condition.name)
                                       condition <- NameCleanUp(name=condition, removestrings)
                                       } else {condition <- NULL}

  #We need to think about naming convention, especially for including condition.

  if (exists("experiment")){AggregateName <- paste(name, experiment, sep = "_")
  } else {AggregateName <- name}

  StorageLocation <- file.path(outpath, AggregateName)

  mff <- gs_pop_get_data(x, marginsubset)
  df <- exprs(mff[[1]])
  TheDF <- data.frame(df, check.names = FALSE)
  DFNames <- colnames(TheDF[,-grep("Time|FS|SC|SS|Original|W$|H$", names(TheDF))])
  PlotNumber <- length(DFNames)

  ff <- gs_pop_get_data(x, gatesubset)

  if (ycolumn == "ALL"){

    .UniversalIterator <- function(x, x_ff,
                                   TheDF, yValue, columnlist, gatelines,
                                   reference, clearance, bins, AltNameX,
                                   AltNameY, colorX, colorY){

      columnlist <- DFNames[DFNames != x] # Remove the universal Y value

      Plots <- map(.x = columnlist, .f = Utility_GeneralGating, name = name, ff = ff, yValue = x, columnlist = DFNames,
                   TheDF = TheDF, gatelines = gatelines, reference = reference, clearance=clearance, bins=bins)

      #Plots <- flatten(Plots)
      #Plots1 <- Plots
    }

    Plots <- map(.x=DFNames, .f = .UniversalIterator, x_ff=ff,
                 TheDF=TheDF, yValue=ycolumn, columnlist=DFNames, gatelines=gatelines,
                 reference=reference, clearance=clearance, bins=bins)

    Plots <- flatten(Plots)


  } else {
    columnlist <- DFNames[DFNames != ycolumn]
    Plots <- map(.x = columnlist, .f = Utility_GeneralGating, name = name, ff = ff, yValue = ycolumn, columnlist = DFNames,
                 TheDF = TheDF, gatelines = gatelines, reference = reference, clearance, bins)
  }


  if (pdf == TRUE){

    #Pass directly to Utility_Patchwork

    theList <- Plots
    theListLength <- length(Plots)

    thecolumns <- 4
    therows <- 3
    theoreticalitems <- therows*thecolumns

    DecimalLeftover <- (PlotNumber/theoreticalitems) %% 1
    AdditionalSpaces <- theoreticalitems*DecimalLeftover

    split_list <- function(input_list, chunk_size) {
      split(input_list, ceiling(seq_along(input_list) / chunk_size))
    }

    sublists <- split_list(theList, theoreticalitems)

    pdf(file = paste(StorageLocation, ".pdf", sep = "", collapse = NULL), width = width, height = height)

    for(i in sublists){p <- wrap_plots(i, ncol = thecolumns, nrow = therows, widths = 0.8, heights = 0.8)
    print(p)
    }

    dev.off()

  }

    return(Plots)
}









