#' Runs openTSNE via Basilisk
#'
#' @param x A gating set object
#' @param sample.name Keyword for which sample name is stored
#' @param removestrings A list of things to remove from sample.name
#' @param subsets The subset of interest from gating hierarchy
#' @param columns A subset of columns to pass instead
#' @param notcolumns A subset of columns to remove
#' @param subsample If downsampling is wanted.
#' @param outpath Location to store new .fcs files
#'
#' @importFrom basilisk basiliskStart
#' @importFrom basilisk basiliskStop
#' @importFrom basilisk basiliskRun
#' @importFrom reticulate import
#' @importFrom reticulate py_to_r
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
#' @examples NULL

Utility_openTSNE <- function(x, sample.name, removestrings, subsets, columns, notcolumns, subsample, export, outpath){
  #Python Environment Startup # Early in Case of Failure
  proc <- basiliskStart(env1)
  on.exit(basiliskStop(proc))

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

  if (!is.null(notcolumns)){CleanedDF1 <- CleanedDF1 %>% select(-all_of(columns))
  } else {CleanedDF1 <- CleanedDF1}

  X <- CleanedDF1

  TheopenTSNE <- basiliskRun(proc, .Internal_openTSNE, X=X)

  new_fcs <- Utility_ColAppend(ff=newff, DF=DF, columnframe=TheopenTSNE)

  TheFileName <- paste0(alternatename, "_DR.fcs")

  fileSpot <- file.path(outpath, TheFileName)

  if(export == TRUE){write.FCS(new_fcs, filename = fileSpot, delimiter="#")
  }else{return(new_fcs)}
}

.Internal_openTSNE <-  function(X, ...){
  openTSNE <- import("openTSNE")
  np <- import("numpy")
  pd <- import("pandas", convert = TRUE)

  # Converting r.x to a pandas DataFrame
  df <- pd$DataFrame(X)
  X <- np$asarray(df)

  tsne_operator <- openTSNE$TSNE(perplexity = as.integer(30), metric="euclidean", verbose = TRUE)

  X_tsne <- tsne_operator$fit(X)

  colnames(X_tsne) <- c("openTSNE_1", "openTSNE_2")
  TheopenTSNE <- data.frame(X_tsne)
}
