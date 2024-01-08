#' Main Luciernaga Function, normalized on individual cell level
#'
#' @param x A Gating Set object (ex. gs or gs[[1]])
#' @param subsets The Gating Hierarchy level you will be sampling at
#' @param sample.name Keyword variable which samples are stored (ex. "GUID")
#' @param group.name Keyword variable which groups are stored (ex. "GROUPNAME")
#' @param experiment Provide directly experiment name (ex. "JAN2024")
#' @param experiment.name Keyword variable which experiment information is stored (ex. "TUBENAME")
#' @param stats Whether to take "mean" or "median"
#' @param Kept Whether "Raw" or "Normalized" values are retained in the Luciernaga object.
#' @param external An external autofluorescence to subtract from single colors.
#' @param sourcelocation Location where .fcs creation file is stored
#' @param outpath  Location where created .fcs and .csv files are sent
#' @param artificial Whether an artificial 0 population should be added for a background autofluorescence stand in.
#' @param fcsexport Whether to export .fcs files, TRUE or FALSE
#' @param mainAF Main Autofluorescence Detector (ex. "V7-A")
#' @param AFOverlap Name of data.frame containing the Autofluorescence overlap of individual fluorophores for exclusion
#' @param Beads  Whether the sample is Beads.
#' @param Brightness Whether sum of detectors should be returned.
#' @param Unstained Whether the sample is Unstained.
#'
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom flowCore keyword
#' @importFrom BiocGenerics nrow
#'
#' @return NULL
#' @export
#'
#' @examples NULL
Troubleshoot <- function(x, subsets, sample.name, group.name, experiment = NULL, experiment.name = NULL, stats, Kept, external, sourcelocation, outpath, artificial, fcsexport, mainAF, AFOverlap, Beads, Brightness, Unstained){

  gc()
  x <- x

  # Retrieving .fcs file name
  name <- keyword(x, sample.name)
  name <- gsub(".fcs", "", name)

  if (Unstained == TRUE) {name <- paste0(name, " Unstained")}

  #Retrieving group name #InfantID in this case
  group <- keyword(x, group.name)

  # Guessing Type
  if(str_detect(name, "(Cells)")){Type <- "Cells"} else if(str_detect(name, "(Beads)")){Type <- "Beads"} else {Type <- "NULL"}

  #Additional Name Cleanup #Using Cytek Reference Control Output Names
  name <- gsub(" (Cells)", "", fixed = TRUE, gsub(" (Beads)", "", fixed = TRUE,name))

  #Setting Up an Alternate Name #Removing All Separators
  alternate.name <- name #Before Removing Spaces et al
  alternate.name <- gsub(" ", "", gsub("(", "", fixed = TRUE, gsub(")", "", fixed = TRUE, gsub("_", "", fixed = TRUE, gsub("-", "", fixed = TRUE, alternate.name)))))

  #Retrieving Experiment Info #Switched to an exist statement.
  if(!is.null(experiment)){Experiment <- experiment
  } else {experiment <- keyword(x, experiment.name)
  experiment <- gsub("-", "", fixed = TRUE, experiment) #Some Cleanup
  Experiment <- experiment}
  #suppressWarnings(rm(experiment)) #Being Used Somewhere

  #Retrieving the exprs data for the target population
  ff <- gs_pop_get_data(x, subsets)
  startingcells <- nrow(ff)[[1]]
  df <- exprs(ff[[1]])
  DF <- as.data.frame(df, check.names = FALSE)

  #For Future Column Reordering
  OriginalColumnsVector <- colnames(DF)
  OriginalColumns <- colnames(DF)
  OriginalColumns <- data.frame(OriginalColumns)
  OriginalColumnsIndex <- OriginalColumns %>% mutate(IndexLocation = 1:nrow(.))

  #For Future Row Reordering
  Backups <- DF %>% mutate(Backups = 1:nrow(DF)) %>% select(Backups)

  #Stashing Away Time FSC SSC For Later Use
  StashedDF <- DF[,grep("Time|FS|SC|SS|Original|W$|H$", names(DF))]
  StashedDF <- cbind(Backups, StashedDF)

  #Consolidating Columns Going Forward
  CleanedDF <- DF[,-grep("Time|FS|SC|SS|Original|W$|H$", names(DF))]
  BackupNames <- colnames(CleanedDF)
  n <- CleanedDF
  #Normalizing By Peak Detector
  n[n < 0] <- 0
  A <- do.call(pmax, n)
  Normalized <- n/A
  Normalized <- round(Normalized, 1)
  colnames(Normalized) <- gsub("-A", "", colnames(Normalized))

  #Towards Generalizing
  ColsN <- ncol(n)
  ColsNormalized <- ncol(Normalized)
  StartNormalizedMergedCol <- ColsN + 1
  EndNormalizedMergedCol <- ColsNormalized*2

  #Bringing Together Raw And Normalized and deriving Autofluorescence Negative
  WorkAround <- cbind(n, Normalized)

  #Retrieving Main Auto fluorescent Channels signature
  MainAF <- mainAF
  MainAF <- gsub("-A", "", MainAF)

  This <- WorkAround %>% filter(.data[[MainAF]] == 1) %>% select(all_of(1:ColsN))

  if(stats == "mean"){Samples <- This %>% summarize_all(mean) #%>% select(-Backups)
  } else if (stats == "median"){Samples <- This %>% summarize_all(median) #%>% select(-Backups)
  } else(print("NA"))

  #Deriving Peak Detector Counts and Detectors of Interest
  na_counts <- colSums(is.na(Normalized))
  Normalized[is.na(Normalized)] <- 0
  Counts <- colSums(Normalized == 1)
  PeakDetectorCounts <- data.frame(Fluors = names(Counts), Counts = Counts)
  rownames(PeakDetectorCounts) <- NULL
  cutoff <- startingcells*0.0075

  return(startingcells)

}
