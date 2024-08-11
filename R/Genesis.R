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
#' @return
#' @export
#'
#' @examples
Genesis <- function(x, ff, minimalfcscutoff=0.05, AggregateName,
                    OriginalStart, OriginalEnd, outpath, stats){

  FlowFrameTest <- ff[[1, returnType = "flowFrame"]]
  original_p <- parameters(FlowFrameTest)
  original_d <- keyword(FlowFrameTest)

  x$Cluster <- factor(x$Cluster)
  ZZZ <- data.frame(table(x$Cluster))
  ZZZ <- ZZZ %>% dplyr::arrange(desc(Freq))
  colnames(ZZZ)[1] <- "Cluster"
  colnames(ZZZ)[2] <- "Count"
  fcs_cutoff <- nrow(x)*minimalfcscutoff
  fcs_clusters <- ZZZ %>% filter(Count > fcs_cutoff) %>% pull(Cluster)

  Data <- x

}


map(x=fcs_clusters, .f=InternalGenesis, AggregateName)
#x=fcs_clusters[1]

InternalGenesis <- function(x, Data, AggregateName, OriginalStart, OriginalEnd, stats){
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


}





RelativeBrightness <- data.frame()

#p <- fcs_clusters[[1]]
for(p in fcs_clusters){



  FCSSubset
  #HowMany <- nrow(FCSSubset)

  if (artificial == TRUE){
    MeanFCS <- colMeans(FCSSubset)
    MeanFCS <- data.frame(t(MeanFCS), check.names = FALSE)
    mutateCols <- MeanFCS[,-grep("Time|FS|SC|SS|Original|W$|H$", names(MeanFCS))] %>% colnames(.)
    MeanFCS <- MeanFCS %>% mutate(across(all_of(mutateCols), ~ifelse(. >= 0, 0, .)))
    MeanFCS$Time <- round(MeanFCS$Time, 1)
    ArtificialNegative <- MeanFCS[rep(1, each = 1000),]
    rownames(ArtificialNegative) <- NULL
    FCSSubset <- rbind(FCSSubset, ArtificialNegative)
  }

  FCSSubset <- data.matrix(FCSSubset)
  new_fcs <- new("flowFrame", exprs=FCSSubset, parameters=original_p, description=original_d)
  fileSpot <- paste(outpath, group, FCSName, ".fcs", sep ="")

  fileSpot <- gsub("Reference GroupDR", "", fileSpot)

  write.FCS(new_fcs, filename = fileSpot, delimiter="#")
}

if (Brightness == TRUE){
  RelativeBrightness1 <- Utility_RelativeBrightness(RelativeBrightness)
  csvSpot <- paste(outpath, "RelativeBrightness", alternate.name, ".csv", sep ="")
  write.csv(RelativeBrightness1, csvSpot, row.names = FALSE)
}
