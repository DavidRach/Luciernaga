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
#' @keywords internal
UnstainedSignatures <- function(x, WorkAround1, alternatename, ColsN,
  StartNormalizedMergedCol, EndNormalizedMergedCol, Increments, Verbose = FALSE,
  LocalMaximaRatio=0.15, SecondaryPeaks=2){

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
  MyData$Cluster <- paste(DetectorName, "10-", sep = "_")

  if (length(These) > 3){stop(
    "Only currently set up to handle up to 4 fluorescence peaks per fluorophore")
  } else if (length(These) == 3){second <- These[[1]]
  third <- These[[2]]
  fourth <- These[[3]]
  } else if (length(These) == 2){second <- These[[1]]
  third <- These[[2]]
  } else if (length(These) == 1){second <- These[[1]]
  } else if (length(These) == 0){
    if (Verbose == TRUE) {message("No second peak")}}

  if (length(These) >= 1){MyData <- MyData %>% mutate(Cluster = case_when(
    near(MyData[[second]], 0.0) ~ paste0(MyData$Cluster, second, "_00-"),
    near(MyData[[second]], 0.1) ~ paste0(MyData$Cluster, second, "_01-"),
    near(MyData[[second]], 0.2) ~ paste0(MyData$Cluster, second, "_02-"),
    near(MyData[[second]], 0.3) ~ paste0(MyData$Cluster, second, "_03-"),
    near(MyData[[second]], 0.4) ~ paste0(MyData$Cluster, second, "_04-"),
    near(MyData[[second]], 0.5) ~ paste0(MyData$Cluster, second, "_05-"),
    near(MyData[[second]], 0.6) ~ paste0(MyData$Cluster, second, "_06-"),
    near(MyData[[second]], 0.7) ~ paste0(MyData$Cluster, second, "_07-"),
    near(MyData[[second]], 0.8) ~ paste0(MyData$Cluster, second, "_08-"),
    near(MyData[[second]], 0.9) ~ paste0(MyData$Cluster, second, "_09-"),
    near(MyData[[second]], 1.0) ~ paste0(MyData$Cluster, second, "_10-")))
  }

  if(length(These) >= 2){MyData <- MyData %>% mutate(Cluster = case_when(
    near(MyData[[third]], 0.0) ~ paste0(MyData$Cluster, third, "_00"),
    near(MyData[[third]], 0.1) ~ paste0(MyData$Cluster, third, "_01"),
    near(MyData[[third]], 0.2) ~ paste0(MyData$Cluster, third, "_02"),
    near(MyData[[third]], 0.3) ~ paste0(MyData$Cluster, third, "_03"),
    near(MyData[[third]], 0.4) ~ paste0(MyData$Cluster, third, "_04"),
    near(MyData[[third]], 0.5) ~ paste0(MyData$Cluster, third, "_05"),
    near(MyData[[third]], 0.6) ~ paste0(MyData$Cluster, third, "_06"),
    near(MyData[[third]], 0.7) ~ paste0(MyData$Cluster, third, "_07"),
    near(MyData[[third]], 0.8) ~ paste0(MyData$Cluster, third, "_08"),
    near(MyData[[third]], 0.9) ~ paste0(MyData$Cluster, third, "_09"),
    near(MyData[[third]], 1.0) ~ paste0(MyData$Cluster, third, "_10")))
  }

  if(length(These) >= 3){MyData <- MyData %>% mutate(Cluster = case_when(
    near(MyData[[fourth]], 0.0) ~ paste0(MyData$Cluster, fourth, "_00"),
    near(MyData[[fourth]], 0.1) ~ paste0(MyData$Cluster, fourth, "_01"),
    near(MyData[[fourth]], 0.2) ~ paste0(MyData$Cluster, fourth, "_02"),
    near(MyData[[fourth]], 0.3) ~ paste0(MyData$Cluster, fourth, "_03"),
    near(MyData[[fourth]], 0.4) ~ paste0(MyData$Cluster, fourth, "_04"),
    near(MyData[[fourth]], 0.5) ~ paste0(MyData$Cluster, fourth, "_05"),
    near(MyData[[fourth]], 0.6) ~ paste0(MyData$Cluster, fourth, "_06"),
    near(MyData[[fourth]], 0.7) ~ paste0(MyData$Cluster, fourth, "_07"),
    near(MyData[[fourth]], 0.8) ~ paste0(MyData$Cluster, fourth, "_08"),
    near(MyData[[fourth]], 0.9) ~ paste0(MyData$Cluster, fourth, "_09"),
    near(MyData[[fourth]], 1.0) ~ paste0(MyData$Cluster, fourth, "_10")))
  }

  return(MyData)
}


