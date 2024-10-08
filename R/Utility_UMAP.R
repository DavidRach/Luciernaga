#' Runs UMAP (uwot implementation)
#'
#' @param x A gating set object
#' @param sample.name Keyword for which sample name is stored
#' @param removestrings A list of things to remove from sample.name
#' @param subset The subset of interest from gating hierarchy
#' @param columns Which columns to use. Not combinable with notcolumns, use one or other.
#' @param notcolumns Which columns not to use. Not combinable with columns,
#'  use one or other.
#' @param subsample If downsampling is wanted.
#' @param outpath Location to store new .fcs files
#' @param metric umap Argument, default is set to "euclidean"
#' @param n_neighbors umap Argument, default is set to 15
#' @param min_dist umap Argument, default is set to 0.5
#' @param export When set to TRUE returns fcs files to specified outpath.
#' @param ... Other arguments to pass to umap()
#'
#' @importFrom uwot umap
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

Utility_UMAP <- function(x, sample.name, removestrings, subset, columns=NULL,
  notcolumns=NULL, subsample=NULL, export=FALSE, outpath=NULL, metric = "euclidean",
  n_neighbors = 15, min_dist = 0.5, ...){
  # Retrieving the metadata # the abbreviated version
  name <- keyword(x, sample.name)
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

  TheUMAP <- umap(X, metric=metric, n_components=2, n_neighbors=n_neighbors,
                  min_dist=min_dist, ...)

  colnames(TheUMAP) <- c("UMAP_1", "UMAP_2")
  TheUMAP <- data.frame(TheUMAP)

  new_fcs <- Utility_ColAppend(ff=newff, DF=DF, columnframe=TheUMAP, shift = TRUE)

  TheFileName <- paste0(alternatename, "_Dimensionality.fcs")

  if (!is.null(outpath)) {fileSpot <- file.path(outpath, TheFileName)}

  if (export == TRUE) {write.FCS(new_fcs, filename = fileSpot, delimiter="#")
  } else {return(new_fcs)}
}


