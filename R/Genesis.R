Genesis(x=Reintegrated1, ff=ff, AggregateName=AggregateName,
        OriginalStart = OriginalStart, OriginalEnd = OriginalEnd,
        outpath=NULL)


#' Internal for Utility_SingleColorQC, creates .fcs files
#'
#' @param x The data.frame of Luciernaga data.
#' @param ff A cytoframe object.
#' @param outpath Location to export the fcs files to
#'
#' @importFrom flowWorkspace keyword
#'
#' @noRd

Genesis <- function(x, ff, minimalfcscutoff=0.05, AggregateName,
                    OriginalStart, OriginalEnd, outpath, stats){

  # Replicate the Original FCS Parameters
  FlowFrameTest <- ff[[1, returnType = "flowFrame"]]
  original_p <- parameters(FlowFrameTest)
  original_d <- keyword(FlowFrameTest)

  # Figure out what clusters to split from the file.
  x$Cluster <- factor(x$Cluster)
  ZZZ <- data.frame(table(x$Cluster))
  ZZZ <- ZZZ %>% dplyr::arrange(desc(Freq))
  colnames(ZZZ)[1] <- "Cluster"
  colnames(ZZZ)[2] <- "Count"
  fcs_cutoff <- nrow(x)*minimalfcscutoff
  fcs_clusters <- ZZZ %>% filter(Count > fcs_cutoff) %>% pull(Cluster)

  Data <- x

  Brightness <- map(x=fcs_clusters, .f=InternalGenesis, AggregateName) %>% bind_rows()
  #x=fcs_clusters[1]

  if (Brightness == TRUE){
    RelativeBrightness <- Utility_RelativeBrightness(Brightness)
    CSVName <- paste0("RelativeBrightness", AggregateName, ".csv")
    CSVSpot <- file.path(outpath, CSVName)
    write.csv(RelativeBrightness, CSVSpot, row.names = FALSE)
  }

}


#' Internal for Utility_SingleColorQC, creates .fcs files
#'
#' @param x
#' @param Data
#' @param AggregateName
#' @param OriginalStart
#' @param OriginalEnd
#' @param stats
#' @param Samples
#' @param NegativeType
#' @param TotalNegatives
#'
#' @noRd

InternalGenesis <- function(x, Data, AggregateName, OriginalStart, OriginalEnd, stats,
                            Samples, NegativeType, TotalNegatives = 500){

  internalstrings <- c("-", "_")
  FCSname <- Luciernaga:::NameCleanUp(x, removestrings=internalstrings)
  FCSName <- paste(AggregateName, FCSname, sep = "_")
  #FCSName

  FCSSubset <- Data %>% filter(Cluster %in% x)

  # If Return Type No Add Ons
  RawFCSSubset <- FCSSubset %>% select(all_of(OriginalStart:OriginalEnd))

  HowBright <- AveragedSignature(RawFCSSubset, stats=stats)
  HowBright <- cbind(x, HowBright)
  colnames(HowBright)[1] <- "Cluster"
  #HowBright #Exported to bind_row with data.frame.

  RawFCSSubset

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

  if (export == TRUE) {write.FCS(new_fcs, filename = fileSpot, delimiter="#")
  }

  return(HowBright)
}



