#' Intended to assess changes at particular detectors within signatures, uses
#' LuciernagaQC internals with externally provided check at these points
#' 
#' @param x A GatingSet object being iterated on
#' @param subsets The node for the events of interest
#' @param CheckThese A vector of the detectors of interest
#' @param sample.name Keyword variable which samples are stored (ex. "GUID")
#' @param inverse.transform Default is FALSE
#' @param Increments Default is 0.1
#' @param stats Default is median
#' 
#' @importFrom flowCore keyword exprs
#' @importFrom flowWorkspace gs_pop_get_data
#' @importFrom BiocGenerics nrow
#' @importFrom dplyr mutate select across arrange desc pull
#' relocate bind_rows left_join
#' @importFrom tidyselect all_of where
#' @importFrom purrr map
#' 
#' @return A LuciernagaQC data style output for detectors of interest
#' 
#' @export
#' 
#' @examples
#' A <- 2+2
Luciernaga_SpotCheck <- function(x, subsets, sample.name, CheckThese,
   inverse.transform=FALSE, Increments=0.1, stats="median"){
  
  if (length(sample.name) == 2){
    first <- sample.name[[1]]
    second <- sample.name[[2]]
    first <- keyword(x, first)
    second <- keyword(x, second)
    name <- paste(first, second, sep="_")
  } else {name <- keyword(x, sample.name)}

  ff <- gs_pop_get_data(x, subsets, inverse.transform=inverse.transform)
  startingcells <- nrow(ff)[[1]]
  DF <- as.data.frame(exprs(ff[[1]]), check.names=FALSE)
  Backups <- DF |> mutate(Backups = 1:nrow(DF)) |> select(Backups)
  OriginalColumns <- colnames(DF)
  OriginalColumns <- data.frame(OriginalColumns)
  OriginalColumnsIndex <- OriginalColumns %>% mutate(IndexLocation = 1:nrow(.))
  OriginalColumnsVector <- colnames(DF)
  StashedDF <- DF[,grep("Time|FS|SC|SS|Original|W$|H$", names(DF))]
  StashedDF <- cbind(Backups, StashedDF)
  n <- DF[,-grep("Time|FS|SC|SS|Original|W$|H$", names(DF))]

  n[n < 0] <- 0
  A <- do.call(pmax, n)
  Normalized <- n/A
  colnames(Normalized) <- gsub("-A", "", colnames(Normalized))

  CheckThese <- gsub("-A", "", CheckThese)
  Normalized <- Normalized |> select(all_of(CheckThese))
  Normalized <- Normalized %>%
    mutate(across(where(is.numeric), ~ ceiling(. / Increments) * Increments))
  na_counts <- colSums(is.na(Normalized))
  Normalized[is.na(Normalized)] <- 0
  Counts <- colSums(Normalized == 1)
  PeakDetectorCounts <- data.frame(Fluors = names(Counts), Counts = Counts)
  rownames(PeakDetectorCounts) <- NULL
  PeakDetector <- PeakDetectorCounts |> arrange(desc(Counts)) |>
    slice(1) |> pull(Fluors)

  Normalized1 <- LuciernagaClustering(MyData=Normalized,
     These=CheckThese, DetectorName=PeakDetector, SkipMain=TRUE)
  
  ColsN <- ncol(n)
  ColsNormalized <- ncol(Normalized1)
  StartNormalizedMergedCol <- ColsN + 1
  EndNormalizedMergedCol <- ColsN + ColsNormalized

  WorkAround <- cbind(n, Normalized1)

  TheData <- data.frame(table(WorkAround$Cluster))
  TheData <- TheData |> arrange(desc(Freq))
  colnames(TheData)[1] <- "Cluster"
  colnames(TheData)[2] <- "Count"
  #Data

  name <- unname(name)[[1]]
  TheData <- TheData |> mutate(Sample=name)
  TheData <- TheData |> relocate(Sample, .before=Cluster)
  TheClusters <- TheData |> pull(Cluster)

  TheSummary <- map(.x=TheClusters, .f=LuciernagaSmallReport,
    Data=WorkAround, RetainedType="raw", ColsN=ColsN,
    StartNormalizedMergedCol=StartNormalizedMergedCol,
    EndNormalizedMergedCol=EndNormalizedMergedCol,
    stats=stats) |> bind_rows()

  FinalData <- left_join(TheData, TheSummary, by = "Cluster")
  return(FinalData)
}