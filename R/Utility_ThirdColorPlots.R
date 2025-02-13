#' Highlight the location of a particular cell population on a given bi-exponential
#' axis.
#'
#' @param x A GatingSet Object
#' @param subset Desired Gate of Interest
#' @param xaxis X-axis Marker
#' @param yaxis Y-axis Marker
#' @param zaxis The Marker you want visible
#' @param splitpoint Point that denotes positive and negative for that marker.
#' @param thecolor What color should positive cells be
#' @param sample.name The keyword that determines a specimens name
#' @param removestrings A list of string characters to remove from the name
#' @param tilesize Equivalent of bin, determines the height and width of the tile dots.
#' @param reference Provide Coereba reference when splitpoint equals crossreference
#' @param FactorNames Provide factor names when splitpoint equals categorical
#'
#' @return A ggplot object with z-axis cells highlighted from background cells
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
#' UnmixedFCSFiles <- Unmixed_FullStained[1]
#' UnmixedCytoSet <- load_cytoset_from_fcs(UnmixedFCSFiles[1],
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
#' SinglePlot <-  Utility_ThirdColorPlots(x=UnmixedGatingSet[1],
#'   subset = "live", xaxis="BUV496-A", yaxis = "Spark Blue 550-A",
#'   zaxis ="BUV805-A", splitpoint = "continuous", sample.name = "GROUPNAME",
#'   removestrings = removestrings,thecolor = "blue")
#'
Utility_ThirdColorPlots <- function(x, subset, xaxis, yaxis, zaxis, splitpoint = 100,
  thecolor="red", sample.name, removestrings, tilesize=0.7, FactorNames = NULL,
  reference = NULL){

  AggregateName <- NameForSample(x=x, sample.name=sample.name,
                                 removestrings=removestrings)


  # Retrieving margin info for the x specimen
  Margin <- gs_pop_get_data(x, subset)
  df <- exprs(Margin[[1]])
  TheDF <- data.frame(df, check.names = FALSE)

  if (is.data.frame(splitpoint)){message("Splitpoint is a Dataframe")
      TheDF <- TheDF %>% mutate(NewID = row_number())
      FilteringDF <- TheDF
      FilterArguments <- splitpoint %>% pull(Fluorophore)
      nrow(FilteringDF)
      #i <- FilterArguments[1]

      for (i in FilterArguments){
        zaxis <- splitpoint %>% filter(Fluorophore %in% i) %>% pull(Fluorophore)
        TheSplit <- splitpoint %>% filter(Fluorophore %in% i) %>% pull(Splitpoint)
        FilteringDF <- FilteringDF %>% dplyr::filter(.data[[zaxis]] >= TheSplit)
      }

      if (!nrow(FilteringDF)>0){stop("Nothing left after filtering arguments")}

      RetainedIDs <- FilteringDF %>% select(NewID) %>% pull()
      TheSubset <- FilteringDF %>% select(-NewID)
      TheBackground <- TheDF %>% filter(!NewID %in% RetainedIDs) %>% select(-NewID)

      Plot <- ggplot() +
        geom_tile(data = TheBackground, aes(x = .data[[xaxis]], y = .data[[yaxis]]),
        width = tilesize, height = tilesize, color = "lightgray", fill="lightgray") +
        geom_tile(data = TheSubset, aes(x = .data[[xaxis]], y = .data[[yaxis]]),
        width = tilesize, height = tilesize, color = thecolor, fill=thecolor) +
        scale_fill_gradient(low="lightgray", high=thecolor) + theme_bw() +
        labs(title = AggregateName) + theme(strip.background = element_blank(),
        strip.text.x = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.title = element_text(
        size = 10, face = "bold"), legend.position = "none")
    return(Plot)

  } else if (splitpoint == "crossreference"){

    message("Splitpoint is a crossreference")

  } else if (splitpoint == "continuous"){message("Splitpoint is a continuous")

    Plot <- ggplot(TheDF) + geom_tile(aes(x = .data[[xaxis]], y = .data[[yaxis]],
      fill = .data[[zaxis]]), width = tilesize, height = tilesize,
      color = "lightgray") + scale_fill_gradient(low = "lightgray", high = thecolor) +
      theme_bw() + labs(title = AggregateName) + theme(
      strip.background = element_blank(), strip.text.x = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      axis.title = element_text(size = 10, face = "bold"), legend.position = "none")
    return(Plot)
  } else if (splitpoint == "categorical"){

    TheSubset <- TheDF %>% dplyr::filter(.data[[zaxis]] %in% FactorNames)
    TheBackground <- TheDF %>% dplyr::filter(!.data[[zaxis]] %in% FactorNames)

    Plot <- ggplot() + geom_tile(data = TheBackground, aes(x = .data[[xaxis]],
      y = .data[[yaxis]],fill = .data[[zaxis]]), width = tilesize, height = tilesize,
      color = "lightgray") + geom_tile(data = TheSubset, aes(x = .data[[xaxis]],
      y = .data[[yaxis]], fill = .data[[zaxis]]), width = tilesize, height = tilesize,
      color = thecolor) + scale_fill_gradient(low="lightgray", high=thecolor) +
      theme_bw() + labs(title = AggregateName) + theme(
      strip.background = element_blank(), strip.text.x = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      axis.title = element_text(size = 10, face = "bold"), legend.position = "none")

    return(Plot)
  } else if (is.numeric(splitpoint)){splitpoint <- splitpoint
  } else if (is.character(splitpoint)){splitpoint <- as.numeric(splitpoint)
  } else {stop("splitpoint argument not recognized")}


  if (is.numeric(splitpoint)){
  TheSubset <- TheDF %>% dplyr::filter(.data[[zaxis]] >= splitpoint)
  TheBackground <- TheDF %>% dplyr::filter(.data[[zaxis]] < splitpoint)

  Plot <- ggplot() + geom_tile(data = TheBackground, aes(x = .data[[xaxis]],
    y = .data[[yaxis]], fill = .data[[zaxis]]), width = tilesize, height = tilesize,
    color = "lightgray") + geom_tile(data = TheSubset, aes(x = .data[[xaxis]],
    y = .data[[yaxis]], fill = .data[[zaxis]]), width = tilesize, height = tilesize,
    color = thecolor) + scale_fill_gradient(low="lightgray", high=thecolor) +
    theme_bw() + labs(title = AggregateName) + theme(
    strip.background = element_blank(), strip.text.x = element_blank(),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    axis.title = element_text(size = 10, face = "bold"), legend.position = "none")

  return(Plot)
  }
}

