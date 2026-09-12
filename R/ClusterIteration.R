
#' Internal for LuciernagaQC SingleStainSignatures
#'
#' @param x A cluster identity in the cluster column
#' @param data A data.frame
#' @param StartNormalizedMergedCol Indicated Start Normalized Columns
#' @param EndNormalizedMergedCol Indicated End Normalized Columns
#' @param ColsN Indicated end of Raw Value Columns
#' @param AggregateName The sample.name derrived name
#' @param Verbose Whether to return intermediate objects
#' @param LocalMaximaRatio Height of peaks to proceed
#' @param SecondaryPeaks Number of Secondary Peaks, default is set to 2.
#'
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom tidyselect all_of
#' @importFrom dplyr mutate
#' @importFrom dplyr relocate
#' @importFrom dplyr arrange
#' @importFrom dplyr left_join
#' @importFrom dplyr pull
#' @importFrom dplyr case_when
#' @importFrom dplyr near
#' @importFrom utils head
#'
#' @return An internal value
#'
#' @noRd
ClusterIteration <- function(x, data, TheDetector, StartNormalizedMergedCol,
                             EndNormalizedMergedCol, ColsN, AggregateName, Verbose,
                             LocalMaximaRatio = 0.15, SecondaryPeaks){



  subset <- data %>% filter(Cluster %in% x)
  StashedIDs <- subset %>% select(Backups)
  TheNormalized <- subset %>% select(-Backups) %>%
    select(all_of(StartNormalizedMergedCol:EndNormalizedMergedCol))
  MyRawData <- subset %>% select(-Backups) %>% select(all_of(1:ColsN))

  #Preparation for Local Maxima
  Conversion <- data.frame(t(TheNormalized), check.names = FALSE)
  Conversion <- cbind(Detectors = rownames(Conversion), Conversion)
  rownames(Conversion) <- NULL

  #Preparing Detector Stand Ins for left_join
  Decoys <- Conversion %>% select(Detectors)
  Decoys <- Decoys %>% mutate(TheDetector = 1:nrow(Decoys)) %>% relocate(
    TheDetector, .before = Detectors)

  #Deriving an average y-vector for local maxima
  Conversion <- Conversion %>% mutate(TheSums = rowSums(.[2:ncol(.)],
    na.rm = TRUE) /(ncol(Conversion) - 1)) %>% relocate(TheSums, .after = Detectors)
  Conversion$Detectors <- 1:nrow(Conversion)
  LocalX <- Conversion$Detectors
  LocalY <- Conversion$TheSums

  #I made it export, now just need to rebuild, then remove extra :
  alternatename <- AggregateName

  PointData <- LocalMaxima(theX = LocalX, theY = LocalY, therepeats = 3,
    w = 3, span = 0.11, alternatename = alternatename, Verbose = Verbose)

  colnames(PointData)[1] <- "TheDetector"
  colnames(PointData)[2] <- "TheHeight"

  Newest2 <- PointData %>% filter(TheHeight > LocalMaximaRatio) %>%
    arrange(desc(TheHeight))
  Assembled <- left_join(Newest2, Decoys, by = "TheDetector")
  if(nrow(Assembled) == 0){
    stop("Failed at Assembled, no local maxima greater than 0.15")}
  These <- Assembled %>% pull(Detectors)

  if (any(These %in% TheDetector)) {These <- These[These != TheDetector]}

  if(length(These) == 0){if (Verbose == TRUE) {message("Solitary Peak")}
  } else if (length(These) > SecondaryPeaks) {
    if (Verbose == TRUE) {
      message("More than ", SecondaryPeaks+1, " peaks. Abbreviated.")}
    These <- head(These, SecondaryPeaks)
  }

  DetectorName <- TheDetector
  MyData <- cbind(StashedIDs, MyRawData, TheNormalized)

  MyData <- LuciernagaClustering(MyData=MyData, These=These, DetectorName=DetectorName)

  return(MyData)
}