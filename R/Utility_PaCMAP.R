#' Runs PaCMAP via Basilisk
#'
#' @param x A gating set object
#' @param sample.name Keyword for which sample name is stored
#' @param removestrings A list of things to remove from sample.name
#' @param subsets The subset of interest from gating hierarchy
#' @param columns A subset of columns to pass instead
#' @param subsample If downsampling is wanted.
#'
#' @importFrom basilisk basiliskStart
#' @importFrom basilisk basiliskStop
#' @importFrom basilisk basiliskRun
#' @importFrom reticulate import
#' @importFrom reticulate py_to_r
#' @importFrom flowWorkspace keyword
#' @importFrom flowWorkspace gs_pop_get_data
#' @importFrom flowCore exprs
#' @importFrom dplyr slice_sample
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#'
#' @return Numbers
#' @export
#'
#' @examples Not applicable

Utility_PaCMAP <- function(x, sample.name, removestrings, subsets, columns, subsample){
  #Python Environment Startup # Early in Case of Failure
  proc <- basiliskStart(env1)
  on.exit(basiliskStop(proc))

  # Retrieving the metadata # the abbreviated version
  name <- keyword(x, sample.name)
  alternatename <- NameCleanUp(name, removestrings)

  #Retrieving the exprs data for my subset population of interest
  ff <- gs_pop_get_data(x, subsets)
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

  MyDimVis <- basiliskRun(proc, .Internal_PaCMAP, X=X)

  return(MyDimVis)
}

.Internal_PaCMAP <-  function(X, ...){
  pacmap <- import("pacmap")
  np <- import("numpy")
  plt <- import("matplotlib.pyplot")
  pd <- import("pandas", convert = TRUE)

  # Converting r.x to a pandas DataFrame
  df <- pd$DataFrame(X)
  X <- np$asarray(df)

  # Initializing the pacmap instance
  embedding <- pacmap$PaCMAP(n_components = 2, n_neighbors=15, MN_ratio=0.5, FP_ratio=2.0)

  # Fit the data (The index of transformed data corresponds to the index of the original data)
  X_transformed <- embedding$fit_transform(X, init="pca")

  X_transformed <- py_to_r(X_transformed)
}
