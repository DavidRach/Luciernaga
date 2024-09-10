#' Querries reference signatures and returns most similar fluorophores by cosine similarity.
#'
#' @param x Name from Sample Column
#' @param data The dataframe object from QC_LibraryParse
#' @param NumberHits Number of most similar fluorophores by cosine.
#'
#' @importFrom dplyr filter
#' @importFrom dplyr slice
#' @importFrom dplyr mutate
#' @importFrom dplyr pull
#' @importFrom dplyr select
#' @importFrom dplyr group_by
#' @importFrom tidyr pivot_wider
#' @importFrom dplyr ungroup
#' @importFrom dplyr bind_rows
#' @importFrom dplyr select_if
#' @importFrom lsa cosine
#' @importFrom tibble rownames_to_column
#' @importFrom tidyselect starts_with
#' @importFrom dplyr arrange
#' @importFrom dplyr slice_head
#'
#' @returns A dataframe of similar fluorophores
#' @export
#'
#' @examples NULL

QC_WhatsThis <- function(x, data, NumberHits) {
  StartingData <- data %>% filter(Sample %in% x)
  if(nrow(StartingData) > 1){StartingData <- StartingData %>% slice(1)}

  Fluorophore <- StartingData %>% mutate(Fluorophore = paste0(
    "ID_", Fluorochrome, "_", Sample)) %>% pull()
  DetectorCols <- StartingData %>% select(where(is.numeric))
  WhoseThis <- cbind(Fluorophore, DetectorCols)

  TotalDetectors <- length(DetectorCols)

  if (TotalDetectors == 64){instrument <- "FiveLaser"
  FileLocation <- system.file("extdata", package = "Luciernaga")
  TheFile <- file.path(FileLocation, "CytekReferenceLibrary5L.csv")
  ReferenceData <- read.csv(TheFile, check.names = FALSE)
  } else if (TotalDetectors == 54){instrument <- "FourLaser"
  FileLocation <- system.file("extdata", package = "Luciernaga")
  TheFile <- file.path(FileLocation, "CytekReferenceLibrary4LUV.csv")
  ReferenceData <- read.csv(TheFile, check.names = FALSE)
  } else if (TotalDetectors == 38){instrument <- "ThreeLaser"
  FileLocation <- system.file("extdata", package = "Luciernaga")
  TheFile <- file.path(FileLocation, "CytekReferenceLibrary3L.csv")
  ReferenceData <- read.csv(TheFile, check.names = FALSE)
  } else {message("No References Found")}

  ReferenceData <- ReferenceData %>% select(-Instrument) %>%
    group_by(Fluorophore) %>% pivot_wider(
      names_from = Detector, values_from = AdjustedY) %>% ungroup()

  #if (!nrow(WhoseThis) == 1){}
  CombinedView <- bind_rows(WhoseThis, ReferenceData)
  Names <- CombinedView %>% select(Fluorophore) %>% pull()
  Numbers <- CombinedView %>% select_if(is.numeric)
  NumericsT <- t(Numbers)
  rownames(NumericsT) <- NULL
  colnames(NumericsT) <- Names

  CosineMatrix <- cosine(NumericsT)
  CosineMatrix <- round(CosineMatrix, 2)
  CosineFrame <- data.frame(CosineMatrix, check.names = FALSE)

  #CosineFrame <- CosineFrame[1,] #If Want To Work With Cols
  CosineFrame <- CosineFrame %>% select(starts_with("ID_"))
  TheData <- rownames_to_column(CosineFrame, var="Fluorophore")
  TheID <- TheData %>% select(starts_with("ID_")) %>% colnames()


  TheHits <- TheData %>% filter(!Fluorophore %in% TheID) %>%
    arrange(desc(.data[[TheID]])) %>% slice_head(n=NumberHits)
  return(TheHits)
}
