Utility_ParallelNbyNPlots <- function(x, sample.name, removestrings, marginsubset,
    gatesubset, ycolumn, bins, clearance, gatelines, reference = NULL, outpath, fileName, pdf){

  # The Marker everything is plotted against
  ycolumn <- ycolumn

  # Retrieving both mapped items
  x <- x
  y <- y

  nameX <- keyword(x, sample.name)
  AltNameX <- NameCleanUp(name = nameX, removestrings)

  nameY <- keyword(y, sample.name)
  AltNameY <- NameCleanUp(name = nameY, removestrings)

  #if(!is.null(experiment)){experiment <- experiment
  #} else {experiment <- keyword(x, experiment.name)}

  #if(!is.null(condition)){condition <- condition
  #} else {condition <- keyword(x, condition.name)}

  #AggregateName <- paste(name, experiment, sep = "_") #Additional for condition (we need to think this through)
  #StorageLocation <- paste(outpath, AggregateName, sep = "/", collapse = NULL)

  # Name and Location of Final PDF
  StorageLocation <- file.path(outpath, fileName)

  # Retrieving margin info for the x specimen
  xMargin <- gs_pop_get_data(x, marginsubset)
  xdf <- exprs(xMargin[[1]])
  TheXDF <- data.frame(xdf, check.names = FALSE)
  X_DFNames <- colnames(TheXDF[,-grep("Time|FS|SC|SS|Original|W$|H$", names(TheXDF))])
  X_PlotNumber <- length(X_DFNames)

  # Retrieving margin info for the y specimen
  yMargin <- gs_pop_get_data(y, marginsubset)
  ydf <- exprs(yMargin[[1]])
  TheYDF <- data.frame(ydf, check.names = FALSE)
  Y_DFNames <- colnames(TheYDF[,-grep("Time|FS|SC|SS|Original|W$|H$", names(TheYDF))])
  Y_PlotNumber <- length(Y_DFNames)

  if (all(X_DFNames == Y_DFNames)){DFNames <- X_DFNames
  } else (stop("The two fcs files do not have matching column parameters"))

  TheDF <- rbind(TheXDF, TheYDF) # Merge the two margin data frames

  # Retrieving gating info for the x and y specimens
  x_ff <- gs_pop_get_data(x, gatesubset) #Sends cytoset forward
  y_ff <- gs_pop_get_data(y, gatesubset)

  if (ycolumn == "ALL"){message("Running all Y axis combinations has not yet been coded ")
  } else {
    columnlist <- DFNames[DFNames != ycolumn] # Remove the universal Y value

    Plots <- map(.x = columnlist, .f = .Internal_ParallelGating, x_ff=x_ff, y_ff=y_ff,
                 TheDF=TheDF, yValue=ycolumn, columnlist=DFNames, gatelines=gatelines,
                 reference=reference, clearance=clearance, bins=bins) #Name
    }

  if (pdf == TRUE){

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
    #length(sublists)

    #sublists[[length(sublists)]] <- c(sublists[[length(sublists)]], rep(plot_spacer(), AdditionalSpaces))

    pdf(file = paste(StorageLocation, ".pdf", sep = "", collapse = NULL), width = 9, height = 7) #Optional Adjustments for Second

    for(i in sublists){p <- wrap_plots(i, ncol = thecolumns, nrow = therows, widths = 0.8, heights = 0.8)
    print(p)
    }

    dev.off()

  }

  return(Plots)
}

.Internal_ParallelGating <- function(x, name, ff, yValue, clearance, bins,
                                  columnlist, TheDF, gatelines, reference = NULL) {

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
     stop("Either theYmax or theXmax didn't exist, and since I didn't think it relavant to
          duplicate this code in the parallel NxN plot when coding, the function now crashed ")
    } else {

      Plot <- as.ggplot(ggcyto(ff, aes(x = .data[[xValue]],
                                           y = .data[[yValue]]), subset = "root") + geom_hex(bins=bins) +
                              coord_cartesian(xlim = c(theXmin, theXmax), ylim = c(theYmin, theYmax),
                                              default = TRUE) + theme_bw() + labs(title = NULL) +
                              theme(strip.background = element_blank(),
                                    strip.text.x = element_blank(), panel.grid.major = element_line(
                                      linetype = "blank"),
                                    panel.grid.minor = element_line(linetype = "blank"),
                                    axis.title = element_text(size = 10, face = "bold"), legend.position = "none"))

  if (gatelines == TRUE){Value <- reference[reference$specimen == name, xValue]
  Plot <- Plot + geom_vline(xintercept = c(seq(0,200,25)), colour = "gray") +
    geom_vline(xintercept = Value, colour = "red")}
  }

  tryCatch({rm("theXmin", "theXmax", "theYmin", "theYmax")})

  return(Plot)
}
