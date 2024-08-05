#' Internal to Utility_SingleColorQC
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom dplyr summarize_all
#' @importFrom dplyr pull
#' @importFrom dplyr arrange
#' @importFrom dplyr relocate
#' @importFrom dplyr left_join
#' @importFrom dplyr case_when
#' @importFrom dplyr rename
#' @importFrom BiocGenerics nrow
#' @importFrom flowWorkspace keyword
#' @importFrom stringr str_detect
#' @importFrom stringr str_split
#' @importFrom flowWorkspace gs_pop_get_data
#' @importFrom flowCore exprs
#' @importFrom purrr map
#' @importFrom purrr set_names
#' @noRd

ClusterIteration <- function(x, data, TheDetector, RatioCutoff, StartNormalizedMergedCol,
                             EndNormalizedMergedCol, ColsN, AggregateName){
  subset <- data %>% filter(Cluster %in% x)

  StashedIDs <- subset %>% select(Backups)
  TheNormalized <- subset %>% select(-Backups) %>% select(all_of(StartNormalizedMergedCol:EndNormalizedMergedCol))
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
                                                        na.rm = TRUE) /(ncol(Conversion) - 1)) %>% relocate(
                                                          TheSums, .after = Detectors)
  Conversion$Detectors <- 1:nrow(Conversion)
  LocalX <- Conversion$Detectors
  LocalY <- Conversion$TheSums

  #I made it export, now just need to rebuild, then remove extra :
  alternatename <- AggregateName

  PointData <- Luciernaga::Utility_LocalMaxima(theX = LocalX, theY = LocalY,
                                               therepeats = 3, w = 3, span = 0.11, alternatename = alternatename)

  colnames(PointData)[1] <- "TheDetector"
  colnames(PointData)[2] <- "TheHeight"

  LocalMaximaRatio <- 0.15 #Possible Relocate Out
  SecondaryPeaks <- 2

  Newest2 <- PointData %>% filter(TheHeight > LocalMaximaRatio) %>% arrange(desc(TheHeight))
  Assembled <- left_join(Newest2, Decoys, by = "TheDetector")
  if(nrow(Assembled) == 0){stop("Failed at Assembled, no local maxima greater than 0.15")}
  These <- Assembled %>% pull(Detectors)
  if (any(These %in% TheDetector)) {These <- These[These != TheDetector]}

  if(length(These) == 0){message("Solitary Peak")
  } else if (length(These) > SecondaryPeaks) {message("More than ", SecondaryPeaks+1, " peaks. Abbreviated.")
    These <- head(These, SecondaryPeaks)
  }

  MyData <- cbind(StashedIDs, MyRawData, TheNormalized)
  MyData$Cluster <- paste(TheDetector, "10-", sep = "_")

  if (length(These) == 3){second <- These[[1]]
  third <- These[[2]]
  fourth <- These[[3]]
  } else if (length(These) == 2){second <- These[[1]]
  third <- These[[2]]
  } else if (length(These) == 1){second <- These[[1]]
  } else if (length(These) == 0){message("No second peak")}

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
