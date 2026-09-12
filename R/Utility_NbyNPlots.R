#' Generate NxN plots for one or all fluorophores.
#'
#' @param x A GatingSet object (ex. gs or gs[[1]])
#' @param sample.name Keyword for which samples names are stored (ex. "GUID")
#' @param removestrings A string of characters to remove from name (
#' ex. c("DR_ILT_2023_", "Cord"))
#' @param experiment Provide directly an experiment name (ex. "Jan2024")
#' @param experiment.name Keyword for which experiment name is stored (
#' ex. "GROUPNAME")
#' @param outpath Location for which to store the generated .pdf
#' @param marginsubset The Gating Hierarchy level from which values will be used to
#'  estimate the plot margins (ex. "nonDebris")
#' @param gatesubset The Gating Hierarchy level of the cells that you want to see
#'  plotted (ex. "lymph")
#' @param ycolumn The ycolumn that you want to see everything plotted by (ex. "APC-A")
#'  or ALL to see all comparisons
#' @param bins Bins for which the plotted cells will be divided into providing
#'  granularity
#' @param clearance The additional ratio added to the margins to avoid clipping
#'  main population but exclude outliers.
#' @param returntype Whether to return as "pdf", "patchwork" or "plots"
#' @param condition Provide a condition name
#' @param condition.name The keyword in the .fcs file storing the condition.name
#' @param gatelines Whether to add estimated gate cutoff lines
#' @param reference Reference for the gate cutoff lines
#' @param width Desired page width for a pdf, default is 9 inches.
#' @param height Desired page height for a pdf, default is 7 inches
#' @param filename Default NULL, provide name to set the filename. 
#'
#' @importFrom flowWorkspace keyword gs_pop_get_data
#' @importFrom flowCore exprs
#' @importFrom patchwork wrap_plots plot_spacer
#' @importFrom purrr map flatten
#' @importFrom utils read.csv
#'
#' @return A one by N plot
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
#' IndividualNxN <- Utility_NbyNPlots(x=UnmixedGatingSet[[1]],
#'   sample.name = "GROUPNAME", removestrings = removestrings,
#'   marginsubset = "lymphocytes", gatesubset = "live",
#'   ycolumn = "Spark Blue 550-A", bins = 70, clearance = 0.2,
#'   gatelines = FALSE, reference = NULL, outpath = StorageLocation,
#'   returntype="patchwork")
#'
Utility_NbyNPlots <- function(x,
                               sample.name = c("GROUPNAME", "TUBENAME"),
                               removestrings = ".fcs",
                               experiment = NULL,
                               experiment.name = NULL,
                               condition = NULL,
                               condition.name = NULL,
                               marginsubset = "root",
                               gatesubset = "root",
                               ycolumn,
                               bins = 100,
                               clearance = 0.2,
                               gatelines = FALSE,
                               reference = NULL,
                               outpath = NULL,
                               returntype = "patchwork",
                               width = 9,
                               height = 7,
                               filename = NULL) {

  if (!is.null(reference)) {
    if (is.data.frame(reference)) {
      reference <- reference
    } else {
      reference <- read.csv(file = reference, check.names = FALSE)
    }
  }

  AggregateName <- NameForSample(x = x, sample.name = sample.name,
                                  removestrings = removestrings,
                                  experiment = experiment,
                                  experiment.name = experiment.name,
                                  condition = condition,
                                  condition.name = condition.name)

  if (length(sample.name) == 2) {
    first <- sample.name[[1]]
    second <- sample.name[[2]]
    first <- keyword(x, first)
    second <- keyword(x, second)
    name <- paste(first, second, sep = "_")
  } else {
    name <- keyword(x, sample.name)
  }

  StorageLocation <- file.path(outpath, AggregateName)

  if (inherits(x, "flowFrame")) {
    message("flowFrame detected, using entire population without gating")
    df <- exprs(x)
  } else {
    mff <- gs_pop_get_data(x, marginsubset)
    df <- exprs(mff[[1]])
  }

  TheDF <- data.frame(df, check.names = FALSE)
  DFNames <- colnames(
    TheDF[, -grep("Time|FS|SC|SS|Original|W$|H$", names(TheDF))])
  PlotNumber <- length(DFNames)

  if (inherits(x, "flowFrame")) {
    ff <- x
  } else {
    ff <- gs_pop_get_data(x, gatesubset)
  }

  if (ycolumn == "ALL") {

    Plots <- map(.x = DFNames, .f = UniversalIterator, x_ff = ff,
                 TheDF = TheDF, yValue = ycolumn, columnlist = DFNames,
                 gatelines = gatelines, reference = reference,
                 clearance = clearance, bins = bins)

    Plots <- flatten(Plots)

  } else {
    columnlist <- DFNames[DFNames != ycolumn]

    if (gatelines == TRUE) {
      columnlist <- intersect(columnlist, colnames(reference))
    }

    Plots <- map(.x = columnlist, .f = GeneralGating, name = name, ff = ff,
                 yValue = ycolumn, columnlist = DFNames, TheDF = TheDF,
                 gatelines = gatelines, reference = reference,
                 clearance = clearance, bins = bins)
  }

  if (!is.null(filename)) {
    AggregateName <- filename
  }

  if (returntype == "pdf") {
    AssembledPlots <- Utility_Patchwork(x = Plots, filename = AggregateName,
                                        outfolder = outpath,
                                        returntype = "pdf")
  }

  if (returntype == "patchwork") {
    AssembledPlots <- Utility_Patchwork(x = Plots, filename = AggregateName,
                                        outfolder = outpath,
                                        returntype = "patchwork")
  }

  if (returntype == "plots") {
    AssembledPlots <- Utility_Patchwork(x = Plots, filename = AggregateName,
                                        outfolder = outpath,
                                        returntype = "plots")
  }

  return(AssembledPlots)
}