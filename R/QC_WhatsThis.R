#' Querries reference signatures and returns most similar fluorophores by cosine similarity.
#'
#' @param x Name in the Sample column you want to filter for
#' @param data A data.frame object from QC_LibraryParse containing Fluorophore name column and numeric detector columns.
#' @param NumberHits Number of most similar fluorophores by cosine.
#'
#' @importFrom dplyr filter
#' @importFrom dplyr slice
#' @importFrom dplyr mutate
#' @importFrom dplyr pull
#' @importFrom dplyr select
#' @importFrom tidyselect where
#' @importFrom dplyr group_by
#' @importFrom tidyr pivot_wider
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr ungroup
#' @importFrom dplyr bind_rows
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

QC_WhatsThis <- function(x, data, NumberHits, returnPlots) {
  StartingData <- data %>% filter(Sample %in% x)
  if(nrow(StartingData) > 1){message("Selecting the first row")
                             StartingData <- StartingData %>% slice(1)}

  CharacterLength <- StartingData %>% select(!where(is.numeric)) %>% length()

  if (CharacterLength == 0){stop("Please add a column Fluorophore with a name")}

  if (CharacterLength > 1){message("Combining character columns")
    Identity <- StartingData %>% select(!where(is.numeric)) %>%
      paste0(collapse = "_")
    Identity <- data.frame(Fluorophore=Identity)
    Identity <- Identity %>% mutate(Fluorophore=paste0("ID_", Fluorophore))
  } else {colnames(StartingData)[1] <- "Fluorophore"
          Identity <- StartingData %>% select(Fluorophore)
          Identity <- Identity %>% mutate(Fluorophore=paste0("ID_", Fluorophore))
  }

  DetectorCols <- StartingData %>% select(where(is.numeric))

  if(any(DetectorCols > 1)){

    message("Normalizing Data for Signature Comparison")

    n <- DetectorCols

    n[n < 0] <- 0
    A <- do.call(pmax, n)
    Normalized <- n/A

    DetectorCols <- Normalized
  }


  WhoseThis <- cbind(Identity, DetectorCols)

  TotalDetectors <- length(DetectorCols)

  ReferenceData <- InstrumentReferences(NumberDetectors=TotalDetectors)

  if (returnPlots == TRUE){ReferenceData1 <- ReferenceData %>% select(-Instrument)
                           WhoseThis1 <- WhoseThis %>%
                             pivot_longer(cols= where(is.numeric), names_to = "Detector",
                                          values_to = "AdjustedY")
                           WhoseThis1 <- WhoseThis1 %>% mutate(Detector = 1:n())
                           ReferenceData1 <- bind_rows(WhoseThis1, ReferenceData1)
  }

  ReferenceData <- ReferenceData %>% select(-Instrument) %>%
    group_by(Fluorophore) %>% pivot_wider(
      names_from = Detector, values_from = AdjustedY) %>% ungroup()

  colnames(ReferenceData) <- colnames(WhoseThis)

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

  TheHits <- TheData %>% dplyr::filter(!Fluorophore %in% TheID) %>%
    arrange(desc(.data[[TheID]])) %>% slice_head(n=NumberHits)

  if (returnPlots==TRUE){

    TheseFluorophores <- TheHits %>% pull(Fluorophore)
    #TheseFluorophores <- c(TheID, TheseFluorophores)
    TheFluorophore <- TheID

    ThePlot <- SimilarFluorPlots(TheseFluorophores=TheseFluorophores,
                                 TheFluorophore=TheFluorophore, data=ReferenceData1)
    ReturnThese <- list(TheHits, ThePlot)
    return(ReturnThese)
  } else {return(TheHits)}
}
