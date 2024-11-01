#' Downsample a gs based on subset
#'
#' @param x A gating set object
#' @param sample.name Keyword specifying sample name
#' @param removestrings Value to be removed from sample name
#' @param subsets The gating hierarchy subset you want to include
#' @param subsample Total number of events to sub-sample from each specimen
#' @param inverse.transform Whether to reverse the GatingSet Transform on the data,
#' default is set to FALSE.
#' @param internal Whether to return as a data.frame (vs. a flow frame or .fcs file)
#' @param export Default is set to FALSE, when TRUE (and internal = FALSE) returns a
#' .fcs file to outpath.
#' @param outpath When export is true, the file.path to where you want the .fcs file
#' stored.
#' @param metadataCols column names from pData to append as metadata for the .fcs, default NULL
#'
#' @importFrom flowWorkspace keyword
#' @importFrom flowWorkspace gs_pop_get_data
#' @importFrom flowCore exprs
#' @importFrom dplyr slice_sample
#' @importFrom flowCore parameters
#' @importFrom flowCore write.FCS
#' @importFrom dplyr mutate
#' @importFrom methods new
#'
#' @return Either a data.frame, a flow.frame or an .fcs file depending on your
#' selected options
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
#' removestrings <- c("DTR_2023_ILT_15_Tetramers-","-Ctrl_Tetramer_Unmixed", ".fcs")
#' StorageLocation <- file.path("C:", "Users", "JohnDoe", "Desktop")
#'
#' CountData <- gs_pop_get_count_fast(UnmixedGatingSet)
#' Counts_Specimen <- CountData %>%
#'  filter(Population %in% "/singletsFSC/singletsSSC/singletsSSCB/nonDebris/lymphocytes/live") %>%
#'  select(name, Count)
#'
#' SingleSample <- Utility_Downsample(UnmixedGatingSet[1],
#'  sample.name = "GROUPNAME", removestrings=removestrings,
#'  subsets = "live", subsample = 2500, internal = FALSE, export = FALSE)
#'
Utility_Downsample <- function(x, sample.name, removestrings,
                                  subsets, subsample=NULL, inverse.transform,
                                  internal = FALSE, export = FALSE, outpath = NULL,
                                  metadataCols=NULL){

  if(length(sample.name) == 1){
    name <- keyword(x, sample.name)
    alternatename <- NameCleanUp(name, removestrings)
  } else {first <- sample.name[1]
          second <- sample.name[2]
          first <- keyword(x, first)
          first <- NameCleanUp(first, removestrings)
          second <- keyword(x, second)
          second <- NameCleanUp(second, removestrings)
          alternatename <- paste(first, second, sep="_")
  }

  #Retrieving the exprs data for my subset population of interest
  ff <- gs_pop_get_data(x, subsets, inverse.transform = inverse.transform)
  df <- exprs(ff[[1]])
  DF <- as.data.frame(df, check.names = FALSE)

  # If down-sampling is specified
  if (!is.null(subsample)) {DF <- slice_sample(DF, n = subsample,
                                              replace = FALSE)
  } else {DF <- DF}

  if (internal == FALSE){
    FCSSubset <- as.matrix(DF)
    FlowFrameTest <- ff[[1, returnType = "flowFrame"]]
    original_p <- parameters(FlowFrameTest)
    original_d <- keyword(FlowFrameTest)
    new_fcs <- new("flowFrame", exprs=FCSSubset, parameters=original_p,
                   description=original_d)

    if (export == TRUE) {
      if (!is.null(outpath)) {TheFileName <- paste0(alternatename, ".fcs")
                              fileSpot <- file.path(outpath, TheFileName)}
                              write.FCS(new_fcs, filename = fileSpot, delimiter="#")
    } else {return(new_fcs)}

  } else {

   if (is.null(metadataCols)){message(alternatename, " has been processed")
     DF <- DF %>% mutate(specimen = alternatename)
     return(DF)
  } else {
    Repeats <- nrow(DF)
    Metadata <- pData(x) %>% select(all_of(metadataCols)) %>%
      slice(rep(1:n(), length.out = Repeats))
    row.names(Metadata) <- NULL
    DF <- cbind(DF, Metadata)
    return(DF)
    }
  }
}




