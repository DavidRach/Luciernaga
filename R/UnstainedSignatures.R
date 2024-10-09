#' Internal to LuciernagaQC
#'
#' @param x A passed detector used for filtering
#' @param WorkAround1 A passed data.frame
#' @param ColsN Indicated end of Raw Value Columns
#' @param StartNormalizedMergedCol Indicated Start Normalized Columns
#' @param EndNormalizedMergedCol Indicated End Normalized Columns
#' @param Increments A numeric to round the normalized bins by. Default is 0.1
#' @param Verbose Whether to return intermediate objects via print and plot for
#' progress monitoring
#' @param LocalMaximaRatio Height of peaks to proceed
#' @param SecondaryPeaks Number of Secondary Peaks, default is set to 2.
#'
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom tidyselect all_of
#' @importFrom dplyr mutate
#' @importFrom dplyr across
#' @importFrom tidyselect where
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
#' @keywords internal
UnstainedSignatures <- function(x, WorkAround1, alternatename, ColsN,
  StartNormalizedMergedCol, EndNormalizedMergedCol, Increments, Verbose = FALSE,
  LocalMaximaRatio=0.15, SecondaryPeaks){

  # Filter for Detector of Interest
  MySubset <- WorkAround1 %>% dplyr::filter(.data[[x]] == 1.000)
  #DetectorPeakCounts(x=MySubset, StartN=StartNormalizedMergedCol,
  #EndN=EndNormalizedMergedCol)

  # Segment into required components
  DetectorName <- x
  alternatename <- alternatename
  Verbose <- Verbose

  StashedIDs <- MySubset %>% select(Backups)
  MySubset <- MySubset %>% select(-Backups)
  MyRawData <- MySubset %>% select(all_of(1:ColsN))
  MyData <- MySubset %>% select(all_of(
    StartNormalizedMergedCol:EndNormalizedMergedCol)) %>%
    mutate(across(where(is.numeric), ~ ceiling(. / Increments) * Increments))

  #Preparation for Local Maxima
  Conversion <- data.frame(t(MyData), check.names = FALSE)
  Conversion <- cbind(Detectors = rownames(Conversion), Conversion)
  rownames(Conversion) <- NULL

  #Preparing Detector Stand Ins for left_join
  Decoys <- Conversion %>% select(Detectors)
  Decoys <- Decoys %>% mutate(TheDetector = 1:nrow(Decoys)) %>%
    relocate(TheDetector, .before = Detectors)

  #Deriving an average y-vector for local maxima
  Conversion <- Conversion %>% mutate(
    TheSums = rowSums(.[2:ncol(.)], na.rm = TRUE) /(ncol(Conversion) - 1)) %>%
    relocate(TheSums, .after = Detectors)

  Conversion$Detectors <- 1:nrow(Conversion)
  LocalX <- Conversion$Detectors
  LocalY <- Conversion$TheSums

  PointData <- LocalMaxima(theX = LocalX, theY = LocalY, therepeats = 3, w = 3,
               span = 0.11, alternatename = alternatename, Verbose=Verbose)

  colnames(PointData)[1] <- "TheDetector"
  colnames(PointData)[2] <- "TheHeight"

  Newest2 <- PointData %>% filter(TheHeight > LocalMaximaRatio) %>%
    arrange(desc(TheHeight))
  Assembled <- left_join(Newest2, Decoys, by = "TheDetector")
  if(nrow(Assembled) == 0){
    stop("Failed at Assembled, no local maxima greater than 0.15")}
  These <- Assembled %>% pull(Detectors)

  if (any(These %in% x)) {These <- These[These != x]}
  if (length(These) == 0) {
    if (Verbose == TRUE) {message("Solitary Peak")}
  } else if (length(These) > SecondaryPeaks) {
    if (Verbose == TRUE) {
      message("More than ", SecondaryPeaks+1, " peaks. Abbreviated.")}
    These <- head(These, SecondaryPeaks)
  }

  MyData <- cbind(StashedIDs, MyRawData, MyData)

  MyData <- LuciernagaClustering(MyData=MyData, These=These, DetectorName=DetectorName)

  return(MyData)
}





