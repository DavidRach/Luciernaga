#' Runs PHATE via Basilisk
#'
#' @param x A gating set object
#' @param sample.name Keyword for which sample name is stored
#' @param removestrings A list of things to remove from sample.name
#' @param subset The subset of interest from gating hierarchy
#' @param columns A subset of columns to pass instead
#' @param notcolumns A subset of columns to remove
#' @param subsample If downsampling is wanted.
#' @param outpath Location to store new .fcs files
#' @param export When set to TRUE returns fcs files to specified outpath.
#' @param ... Additional arguments passed to Phate for fine tuning
#'
#' @importFrom basilisk basiliskStart
#' @importFrom basilisk basiliskStop
#' @importFrom flowCore keyword
#' @importFrom flowWorkspace gs_pop_get_data
#' @importFrom flowWorkspace realize_view
#' @importFrom Biobase exprs
#' @importFrom dplyr slice_sample
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom tidyselect all_of
#' @importFrom basilisk basiliskRun
#' @importFrom flowCore write.FCS
#' @importFrom reticulate import
#' @importFrom reticulate py_to_r
#'
#' @return Numbers
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
#' Markers <- colnames(UnmixedCytoSet)
#' KeptMarkers <- Markers[-grep("Time|FS|SC|SS|Original|-W$|-H$|AF", Markers)]
#' SubsetMarkers <- c("BUV496-A", "BUV805-A", "Pacific Blue-A", "BV711-A",
#'   "BV786-A", "Spark Blue 550-A", "PE-A", "APC-Fire 750-A")
#'
Utility_Phate <- function(x, sample.name, removestrings, subset, columns=NULL,
                          notcolumns=NULL, subsample=NULL, export=FALSE, outpath=NULL,
                          ...){
  #Python Environment Startup # Early in Case of Failure
  proc <- basiliskStart(env1)
  on.exit(basiliskStop(proc))

  # Retrieving the metadata # the abbreviated version
  name <- keyword(x, sample.name)
  # alternatename <- NameCleanUp(name, removestrings)
  alternatename <- NameCleanUp(name, removestrings)

  #Retrieving the exprs data for my subset population of interest
  ff <- gs_pop_get_data(x, subset)
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

  if (!is.null(columns) && !is.null(notcolumns)) {
    stop("Columns and notcolumns are not currently combinable. Pick one")}

  # If external columns interest specified
  if (!is.null(columns)){CleanedDF1 <- CleanedDF %>% select(all_of(columns))
  } else {CleanedDF1 <- CleanedDF}

  if (!is.null(notcolumns)){CleanedDF1 <- CleanedDF1 %>% select(-all_of(columns))
  } else {CleanedDF1 <- CleanedDF1}

  X <- CleanedDF1

  ThePhate <- basiliskRun(proc, .Internal_Phate, X=X, ...)

  new_fcs <- Utility_ColAppend(ff=newff, DF=DF, columnframe=ThePhate, shift = TRUE)

  TheFileName <- paste0(alternatename, "_dimensionality.fcs")

  if (!is.null(outpath)) {fileSpot <- file.path(outpath, TheFileName)}

  if (export == TRUE) {write.FCS(new_fcs, filename = fileSpot, delimiter="#")
  } else {return(new_fcs)}
}

.Internal_Phate <-  function(X, ...){
  phate <- import("phate")
  np <- import("numpy")
  pd <- import("pandas", convert = TRUE)

  # Converting r.x to a pandas DataFrame
  df <- pd$DataFrame(X)
  X <- np$asarray(df)

  phate_operator <- phate$PHATE(n_components=as.integer(2),knn=as.integer(15), decay=as.integer(15), gamma=1)

  Y_phate <- phate_operator$fit_transform(X)

  colnames(Y_phate) <- c("Phate_1", "Phate_2")
  ThePhate <- data.frame(Y_phate)
}
