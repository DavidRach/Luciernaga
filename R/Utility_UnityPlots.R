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

 Plots <- map(GatingSet, .f=.InternalUnity, TheX=TheX, TheY=y, marginsubset=marginsubset, gatesubset=gatesubset,
      sample.name=sample.name, removestrings=removestrings, clearance=clearance, bins=bins, gatelines=gatelines,
      reference=reference)

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

.InternalUnity <- function(x, TheY, TheX, marginsubset, gatesubset, sample.name, removestrings,
                           clearance, bins, gatelines, reference){

  name <- keyword(x, sample.name)
  name <- NameCleanUp(name = name, removestrings)

  ff <- gs_pop_get_data(x, marginsubset)
  df <- exprs(ff[[1]])
  TheDF <- data.frame(df, check.names = FALSE)

  if (!TheX == TheY) {
  YExprsData <- TheDF %>% select(all_of(TheY)) %>% pull()
  theYmin <- YExprsData %>% quantile(., 0.001)
  theYmax <- YExprsData %>% quantile(., 0.999)
  theYmin <- theYmin - abs((clearance*theYmin))
  theYmax <- theYmax + (clearance*theYmax)

  XExprsData <- TheDF %>% select(all_of(TheX)) %>% pull()
  theXmin <- XExprsData %>% quantile(., 0.001)
  theXmax <- XExprsData %>% quantile(., 0.999)
  theXmin <- theXmin - abs((clearance*theXmin))
  theXmax <- theXmax + (clearance*theXmax)
  } else (error("TheX and TheY have the same value"))

  ff1 <- gs_pop_get_data(x, gatesubset)

  if (BiocGenerics::nrow(ff1) < 200) {
    Plot <- as.ggplot(ggcyto(ff1, aes(x = .data[[TheX]], y = .data[[TheY]]),
                             subset = "root") + geom_point(size = 2, alpha = 0.8) + theme_bw() + labs(title = name) +
                        theme(strip.background = element_blank(), strip.text.x = element_blank(),
                              panel.grid.major = element_line(linetype = "blank"),
                              panel.grid.minor = element_line(linetype = "blank"),
                              axis.title = element_text(size = 10, face = "bold"),
                              legend.position = "none"))
  } else {
    Plot <- as.ggplot(ggcyto(ff1, aes(x = .data[[TheX]], y = .data[[TheY]]),
                             subset = "root") + geom_hex(bins=bins) + theme_bw() + labs(title = name) +
                        theme(strip.background = element_blank(), strip.text.x = element_blank(),
                              panel.grid.major = element_line(linetype = "blank"),
                              panel.grid.minor = element_line(linetype = "blank"),
                              axis.title = element_text(size = 10, face = "bold"),
                              legend.position = "none"))
  }



  if (gatelines == TRUE){Value <- reference[reference$specimen == name, TheX]
  Plot <- Plot + geom_vline(xintercept = c(seq(0,200,25)), colour = "gray") +
    geom_vline(xintercept = c(seq(0,200,5)), colour = "white", alpha = 0.1) +
    geom_vline(xintercept = Value, colour = "red")}

  return(Plot)
}
