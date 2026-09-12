#' Queries fluorophore and returns similar fluorophores.
#'
#' @param TheFluorophore The name of the Fluorophore compare, see QC_ReferenceLibrary
#' @param NumberDetectors Number of detectors of the instrument
#' @param NumberHits Number of most similar fluorophores by cosine
#' @param returnSynonymns Returns only fluorophores > 0.98 cosine value, default FALSE
#' @param returnPlots Whether to also return signature plots, default is set FALSE
#' @param returnSynonyms Something
#' @param plotlinecolor Default NULL, otherwise if single line provide desired color
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
#' @importFrom tidyselect where
#' @importFrom lsa cosine
#' @importFrom tibble rownames_to_column
#' @importFrom tidyselect starts_with
#' @importFrom dplyr arrange
#' @importFrom dplyr slice_head
#'
#' @returns A dataframe of similar fluorophores
#' @export
#'
#' @examples
#' Results <- QC_SimilarFluorophores(TheFluorophore="Spark Blue 550",
#'  NumberDetectors=64, returnSynonymns=FALSE, NumberHits = 10, returnPlots=FALSE)

QC_SimilarFluorophores <- function(TheFluorophore, NumberDetectors,
   returnSynonyms=FALSE, NumberHits=10, returnPlots=FALSE, plotlinecolor=NULL) {

  ReferenceData <- Luciernaga:::InstrumentReferences(NumberDetectors=NumberDetectors)
  #nrow(ReferenceData)
  #ReferenceData %>% pull(Fluorophore) %>% unique()
  #ReferenceData1 <- ReferenceData |> unique()
  #nrow(ReferenceData1)
  #ReferenceData1 %>% pull(Fluorophore) %>% unique()

  if (returnPlots == TRUE){ReferenceData1 <- ReferenceData}

  ReferenceData <- ReferenceData |> select(-Instrument) |>
    group_by(Fluorophore) |> pivot_wider(
      names_from = Detector, values_from = AdjustedY) |> ungroup()
  
  #RowNAs <- ReferenceData[rowSums(is.na(ReferenceData)) > 0, ]
  #nrow(RowNAs)

  #CleanNAs <- ReferenceData[rowSums(is.na(ReferenceData)) == 0, ]
  #nrow(CleanNAs)
  #ncol(CleanNAs)-1
  #View(CleanNAs)

  TheAvailableFluors <- ReferenceData |> pull(Fluorophore)
  if (!TheFluorophore %in% TheAvailableFluors) {stop("Fluorophore not found")}

  CombinedView <- ReferenceData
  Names <- CombinedView |> pull(Fluorophore)
  #SanitizedNames <- Luciernaga::NameCleanUp(Names, removestrings = c(",", "-", " ", "."))

  Numbers <- CombinedView %>% select(where(is.numeric))
  #zero_columns <- colSums(Numbers) == 0
  #zero_rows <- rowSums(Numbers) == 0
  #print(which(zero_columns))
  #print(which(zero_rows))
  
  NumericsT <- t(Numbers)
  rownames(NumericsT) <- NULL
  colnames(NumericsT) <- Names
  #colnames(NumericsT) <- SanitizedNames

  #zero_columns <- colSums(NumericsT) == 0
  #zero_rows <- rowSums(NumericsT) == 0
  #print(which(zero_columns))
  #print(which(zero_rows))

  CosineMatrix <- lsa::cosine(NumericsT)
  CosineMatrix <- round(CosineMatrix, 2)
  
  CosineFrame <- data.frame(CosineMatrix, check.names = FALSE)

  CosineFrame <- CosineFrame |> select(all_of(TheFluorophore))
  TheData <- rownames_to_column(CosineFrame, var="Fluorophore")
  TheID <- TheData |> select(all_of(TheFluorophore)) |> colnames()

  if (returnSynonyms == FALSE){
  TheHits <- TheData |> filter(!Fluorophore %in% TheID) |>
    arrange(desc(.data[[TheID]])) |> slice_head(n=NumberHits)
  } else {
    TheHits <- TheData |> filter(!Fluorophore %in% TheID) |>
    arrange(desc(.data[[TheID]])) |> filter(.data[[TheID]] > 0.98)
  }

  if (returnPlots==TRUE){
    TheseFluorophores <- TheHits |> pull(Fluorophore)

    ThePlot <- SimilarFluorPlots(TheseFluorophores=TheseFluorophores,
                                 TheFluorophore=TheFluorophore, data=ReferenceData1,
                                 plotlinecolor=plotlinecolor)
    ReturnThese <- list(TheHits, ThePlot)
    return(ReturnThese)
  } else {return(TheHits)}
  }



