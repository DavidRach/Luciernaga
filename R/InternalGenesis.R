
#' Internal for LuciernagaQC, creates .fcs files
#'
#' @param x Individual fluorescence cluster for filtering
#' @param Data The data.frame containing the many of the above
#' @param AggregateName Passed final name with modifications from name
#' @param outpath Location to export the fcs and .csv files to
#' @param OriginalStart Passed Argument indicating start column for Raw .fcs values
#' @param OriginalEnd Passed argument indicating end column for raw .fcs values
#' @param stats Whether "median" or "mean", default is "median"
#' @param NegativeType Whether to append a negative pop. Args are "artificial",
#' "internal" and "default"
#' @param TotalNegatives How many of the above rows to append, default is set to 500
#' @param Samples When Negative type = "Internal", the data.frame of averaged
#' fluorescence per detector
#' @param ExportType Passed from above, set to "fcs" for fcs.file return
#' @param parameters Passed parameters for .fcs creation
#' @param description Passed description for .fcs creation
#'
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom tidyselect all_of
#' @importFrom dplyr mutate
#' @importFrom dplyr across
#' @importFrom tidyselect one_of
#' @importFrom dplyr bind_cols
#' @importFrom dplyr relocate
#' @importFrom flowCore write.FCS
#'
#' @return An internal value
#'
#' @noRd
InternalGenesis <- function(x, Data, AggregateName, outpath=NULL, OriginalStart,
  OriginalEnd, stats="median", NegativeType="default", TotalNegatives = 500,
  Samples = NULL, ExportType, parameters, description){

  internalstrings <- c("-", "_")
  FCSname <- NameCleanUp(x, removestrings=internalstrings)
  FCSName <- paste(AggregateName, FCSname, sep = "_")
  #FCSName

  FCSSubset <- Data %>% filter(Cluster %in% x)

  # If Return Type No Add Ons
  RawFCSSubset <- FCSSubset %>% select(all_of(OriginalStart:OriginalEnd))

  HowBright <- AveragedSignature(RawFCSSubset, stats=stats)
  Count <- nrow(FCSSubset)
  HowBright <- cbind(x, Count, HowBright)
  colnames(HowBright)[1] <- "Cluster"
  #HowBright #Exported to bind_row with data.frame.
  #RawFCSSubset

  if (NegativeType == "artificial"){
    MeanFCS <- colMeans(RawFCSSubset)
    MeanFCS <- data.frame(t(MeanFCS), check.names = FALSE)
    mutateCols <- MeanFCS[,-grep("Time|FS|SC|SS|Original|W$|H$", names(MeanFCS))] %>%
      colnames(.)
    MeanFCS <- MeanFCS %>% mutate(across(all_of(mutateCols), ~ifelse(. >= 0, 0, .)))
    MeanFCS$Time <- round(MeanFCS$Time, 1)
    ArtificialNegative <- MeanFCS[rep(1, each = TotalNegatives),]
    rownames(ArtificialNegative) <- NULL
    FCSSubset <- rbind(RawFCSSubset, ArtificialNegative)
  }

  if (NegativeType == "samples"){
    if(!is.data.frame(Samples)){stop("Samples needs to be a single row of a data.frame
                                     for just the raw detectors")}
    SamplesCols <- colnames(Samples)
    MeanFCS <- colMeans(RawFCSSubset)
    MeanFCS <- data.frame(t(MeanFCS), check.names = FALSE)
    BackboneCols <- colnames(MeanFCS)
    Residual <- MeanFCS %>% select(-one_of(SamplesCols))
    Combined <- bind_cols(Residual, Samples)
    Combined <- Combined %>% relocate(all_of(BackboneCols))
    SampleNegative <- Combined[rep(1, each = TotalNegatives),]
    rownames(SampleNegative) <- NULL
    FCSSubset <- rbind(RawFCSSubset, SampleNegative)
  }

  if (NegativeType == "default"){
    FCSSubset <- RawFCSSubset
  }

  FCSSubset <- data.matrix(FCSSubset)
  new_fcs <- new("flowFrame", exprs=FCSSubset, parameters=parameters,
                 description=description)

  TheFileName <- paste(AggregateName, FCSname, sep="_")
  TheFileFCS <- paste0(TheFileName, ".fcs")
  if (is.null(outpath)) {outpath <- getwd()}
  fileSpot <- file.path(outpath, TheFileFCS)

  if (ExportType == "fcs") {write.FCS(new_fcs, filename = fileSpot, delimiter="#")}

  return(HowBright)
}