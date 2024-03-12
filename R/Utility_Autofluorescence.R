#' Runs Luciernaga Comparison
#'
#' @param x A gating set object
#' @param sample.name Keyword for which sample name is stored
#' @param removestrings A list of things to remove from sample.name
#' @param subsets The subset of interest from gating hierarchy
#' @param columns A subset of columns to pass instead
#' @param subsample If downsampling is wanted.
#' @param outpath Location to store new .fcs files
#'
#' @importFrom flowWorkspace keyword
#' @importFrom flowWorkspace gs_pop_get_data
#' @importFrom flowCore exprs
#' @importFrom flowCore write.FCS
#' @importFrom dplyr slice_sample
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#'
#' @return Numbers
#' @export
#'
#' @examples Not applicable

Utility_Autofluorescence <- function(x, sample.name, removestrings, subsets, columns, subsample, outpath){
  # Retrieving the metadata # the abbreviated version
  name <- keyword(x, sample.name)
  alternatename <- NameCleanUp(name, removestrings)

  #Retrieving the exprs data for my subset population of interest
  ff <- gs_pop_get_data(x, subsets)
  newff <- realize_view(ff)

  startingcells <- nrow(ff)[[1]] # For Ratio Calculation, may not be needed here
  df <- exprs(ff[[1]])
  DF <- as.data.frame(df, check.names = FALSE)

  # If down-sampling is specified
  if(!is.null(subsample)){DF <- slice_sample(DF, n = subsample,
                                             replace = FALSE)
  } else{DF <- DF}

  # Saving Columns for future column reordering
  OriginalColumnsVector <- colnames(DF)
  OriginalColumns <- colnames(DF)
  OriginalColumns <- data.frame(OriginalColumns)
  OriginalColumnsIndex <- OriginalColumns %>% mutate(IndexLocation = 1:nrow(.))

  # Adding Backups for future row reordering
  Backups <- DF %>% mutate(Backups = 1:nrow(DF)) %>% select(Backups)

  #Stashing Away Time FSC SSC For Later Use
  StashedDF <- DF[,grep("Time|FS|SC|SS|Original|W$|H$", names(DF))]
  StashedDF <- cbind(Backups, StashedDF)

  #Consolidating Columns Going Forward
  CleanedDF <- DF[,-grep("Time|FS|SC|SS|Original|W$|H$", names(DF))]
  BackupNames <- colnames(CleanedDF)

  # If external columns interest specified
  if (!is.null(columns)){CleanedDF1 <- CleanedDF %>% select(all_of(columns))
  } else {CleanedDF1 <- CleanedDF}

  X <- CleanedDF1








  TheFileName <- paste0(alternatename, "_DR.fcs")

  fileSpot <- file.path(outpath, TheFileName)

  write.FCS(new_fcs, filename = fileSpot, delimiter="#")
}

