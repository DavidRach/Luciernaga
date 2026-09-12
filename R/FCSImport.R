
#' Internal for LuciernagaReportFromFCS
#'
#' @param x Passed Fluorophore Name
#' @param data Passed data.frame of Fluorophore with Detectors
#' @param inputfiles List of .fcs files from path
#' @param RetainedType Whether the data.frame contains "raw" or "normalized" values
#' @param TheSummary Whether to return individual cells or summarized by stats.
#' RetainedType Whether the data.frame contains "raw" or "normalized" values
#'
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#' @importFrom stringr str_detect
#' @importFrom flowWorkspace load_cytoset_from_fcs
#' @importFrom purrr map
#' @importFrom dplyr mutate
#' @importFrom dplyr across
#' @importFrom tidyselect where
#' @importFrom dplyr select
#' @importFrom dplyr group_by
#' @importFrom dplyr ungroup
#' @importFrom dplyr arrange
#' @importFrom dplyr left_join
#'
#' @return An internal value
#'
#' @noRd
FCSImport <- function(x, data, inputfiles, RetainedType, TheSummary, stats){

  # For each Fluorophore
  TheDetector <- data %>% dplyr::filter(Fluorophore %in% x) %>% pull(Detector)

  fcs_files <- inputfiles[str_detect(basename(inputfiles), x) &
                            str_detect(basename(inputfiles), ".fcs$")]

  if (x %in% c("PE", "APC")){x <- paste0(x, "-") #ExceptionHandling
  fcs_files <- fcs_files[!str_detect(basename(fcs_files), x)]
  x <- gsub("-", "", x)
  } #ExceptionHandling

  if (!length(fcs_files) == 0){

    cs <- load_cytoset_from_fcs(fcs_files, truncate_max_range = FALSE,
                                transformation = FALSE)
    thex <- x

    # Retrieve exprs data from each .fcs file and create cluster column
    # x <- cs[1]
    TheData <- map(.x = cs, .f = FCSImportFile, Fluorophore = thex) %>%
      bind_rows()

    # Removing Any Artificial Negatives Inserted By Luciernaga
    TheData <- TheData %>%
      mutate(Summed = rowSums(across(where(is.numeric)), na.rm = TRUE))

    TheData <- TheData %>% group_by(Summed) %>% dplyr::filter(n() <= 5) %>%
      ungroup() %>% dplyr::filter(!Summed == 0) %>% dplyr::select(-Summed)

    # Return Summarized Data
    if (RetainedType == "normalized"){
      Cluster <- TheData %>% dplyr::select(Cluster)
      DetectorData <- TheData %>% dplyr::select(-Cluster)
      DetectorData[DetectorData < 0] <- 0
      A <- do.call(pmax, DetectorData)
      Normalized <- DetectorData/A
      Normalized <- round(Normalized, 3)
      TheData <- cbind(Cluster, Normalized)
    }

    if (TheSummary == TRUE){
      TheTable <- data.frame(table(TheData$Cluster), check.names = FALSE)
      TheTable <- TheTable %>% dplyr::arrange(desc(Freq))
      colnames(TheTable)[1] <- "Cluster"
      colnames(TheTable)[2] <- "Count"

      TheClusters <- TheTable$Cluster

      SmallHelper <- function(x, data, stats){
      Cluster <- x
      Internal <- data %>% dplyr::filter(Cluster %in% x) %>% dplyr::select(where(is.numeric))
      Summarized <- AveragedSignature(x=Internal, stats=stats)
      Summarized <- cbind(Cluster, Summarized)
      }

      Summarized <- map(.x=TheClusters, .f=SmallHelper, data = TheData,
                        stats=stats) %>% bind_rows()

      ReturnFrame <- left_join(TheTable, Summarized, by="Cluster")
      return(ReturnFrame)
    } else {return(TheData)}
  }

}