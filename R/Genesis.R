#' Internal for LuciernagaQC, creates .fcs files
#'
#' @param x The data.frame of Luciernaga data.
#' @param ff An individual cytoset object.
#' @param minimalfcscutoff A ratio indicating mininum of the total population needed to split off
#'  into own file, default is set to 0.05
#' @param AggregateName Passed final name with modifications from name
#' @param Brightness Whether to additionally return a brightness .csv to the outpath
#' @param outpath Location to export the fcs and .csv files to
#' @param OriginalStart Passed Argument indicating start column for Raw .fcs values
#' @param OrigingalEnd Passed argument indicating end column for raw .fcs values
#' @param stats Whether "median" or "mean", default is "median"
#' @param NegativeType Whether to append a negative pop. Args are "artificial", "internal" and "default"
#' @param TotalNegatives How many of the above rows to append, default is set to 500
#' @param Samples When Negative type = "Internal", the data.frame of averaged fluorescence per detector
#' @param ExportType Passed from above, set to "fcs" for fcs.file return
#'
#'
#' @importFrom flowCore parameters
#' @importFrom flowWorkspace keyword
#' @importFrom dplyr arrange
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
#' @importFrom utils write.csv
#'
#' @keywords internal

Genesis <- function(x, ff, minimalfcscutoff=0.05, AggregateName,
                    Brightness = FALSE, outpath=NULL, OriginalStart, OriginalEnd,
                    stats = "median", NegativeType="default", TotalNegatives=500,
                    Samples=NULL, ExportType){

  # Replicate the Original FCS Parameters
  FlowFrameTest <- ff[[1, returnType = "flowFrame"]]
  original_p <- parameters(FlowFrameTest)
  original_d <- keyword(FlowFrameTest)

  # Figure out what clusters to split from the file.
  x$Cluster <- factor(x$Cluster)
  ZZZ <- data.frame(table(x$Cluster))
  ZZZ <- ZZZ %>% arrange(desc(Freq))
  colnames(ZZZ)[1] <- "Cluster"
  colnames(ZZZ)[2] <- "Count"
  fcs_cutoff <- nrow(x)*minimalfcscutoff
  fcs_clusters <- ZZZ %>% filter(Count > fcs_cutoff) %>% pull(Cluster)

  Data <- x

  Brightness <- map(x=fcs_clusters, .f=InternalGenesis, Data=Data, AggregateName=AggregateName,
                    outpath=outpath, OriginalStart=OriginalStart, OriginalEnd=OriginalEnd,
                    stats=stats, NegativeType=NegativeType, TotalNegatives=TotalNegatives,
                    Samples=Samples, ExportType=ExportType) %>% bind_rows()

  if (Brightness == TRUE){
    RelativeBrightness <- Utility_RelativeBrightness(Brightness)
    CSVName <- paste0("RelativeBrightness", AggregateName, ".csv")
    CSVSpot <- file.path(outpath, CSVName)
    write.csv(RelativeBrightness, CSVSpot, row.names = FALSE)
  }
}


#' Internal for LuciernagaQC, creates .fcs files
#'
#' @param x Individual fluorescence cluster for filtering
#' @param Data The data.frame containing the many of the above
#' @param AggregateName Passed final name with modifications from name
#' @param outpath Location to export the fcs and .csv files to
#' @param OriginalStart Passed Argument indicating start column for Raw .fcs values
#' @param OriginalEnd Passed argument indicating end column for raw .fcs values
#' @param stats Whether "median" or "mean", default is "median"
#' @param NegativeType Whether to append a negative pop. Args are "artificial", "internal" and "default"
#' @param TotalNegatives How many of the above rows to append, default is set to 500
#' @param Samples When Negative type = "Internal", the data.frame of averaged fluorescence per detector
#' @param ExportType Passed from above, set to "fcs" for fcs.file return
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
#' @noRd

InternalGenesis <- function(x, Data, AggregateName, outpath=NULL, OriginalStart, OriginalEnd,
                            stats="median", NegativeType="default", TotalNegatives = 500,
                            Samples = NULL, ExportType){

  internalstrings <- c("-", "_")
  FCSname <- NameCleanUp(x, removestrings=internalstrings)
  FCSName <- paste(AggregateName, FCSname, sep = "_")
  #FCSName

  FCSSubset <- Data %>% filter(Cluster %in% x)

  # If Return Type No Add Ons
  RawFCSSubset <- FCSSubset %>% select(all_of(OriginalStart:OriginalEnd))

  HowBright <- AveragedSignature(RawFCSSubset, stats=stats)
  HowBright <- cbind(x, HowBright)
  colnames(HowBright)[1] <- "Cluster"
  #HowBright #Exported to bind_row with data.frame.
  #RawFCSSubset

  if (NegativeType == "artificial"){
    MeanFCS <- colMeans(RawFCSSubset)
    MeanFCS <- data.frame(t(MeanFCS), check.names = FALSE)
    mutateCols <- MeanFCS[,-grep("Time|FS|SC|SS|Original|W$|H$", names(MeanFCS))] %>% colnames(.)
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
  new_fcs <- new("flowFrame", exprs=FCSSubset, parameters=original_p, description=original_d)

  TheFileName <- paste(AggregateName, FCSname, sep="_")
  TheFileFCS <- paste0(TheFileName, ".fcs")
  if (is.null(outpath)) {outpath <- getwd()}
  fileSpot <- file.path(outpath, TheFileFCS)

  if (ExportType == "fcs") {write.FCS(new_fcs, filename = fileSpot, delimiter="#")
  }

  return(HowBright)
}



