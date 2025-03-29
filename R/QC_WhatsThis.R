#' Querries reference signatures and returns most similar fluorophores by cosine similarity.
#'
#' @param x Name in the Sample column you want to filter for
#' @param columnname Designating column to filter the x argument from
#' @param data A data.frame object from QC_LibraryParse containing Fluorophore name column and numeric detector columns.
#' @param NumberHits Number of most similar fluorophores by cosine.
#' @param NumberDetectors Default NULL, estimated from number of numeric columns in passed data.
#' @param Normalize Default TRUE, needed for ReferenceLibrary comparison. 
#' @param returnPlots Whether to return signature plots, default is set to FALSE.
#' @param TheFormat Default wider for detectors in columns, specify longer if providing detectors as rows
#' @param detectorcolumn Default NULL, when TheFormat="longer" specify detector column name
#' @param valuecolumn Default NULL, when TheFormat="longer" specify value column name
#'
#' @importFrom dplyr filter
#' @importFrom dplyr slice
#' @importFrom dplyr select
#' @importFrom tidyselect where 
#' @importFrom dplyr rename
#' @importFrom dplyr mutate
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr bind_rows
#' @importFrom dplyr group_by
#' @importFrom tidyr pivot_wider
#' @importFrom dplyr ungroup
#' @importFrom dplyr pull
#' @importFrom lsa cosine
#' @importFrom tidyselect starts_with
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr arrange
#' @importFrom dplyr arrange
#' @importFrom dplyr desc
#' @importFrom dplyr slice_head
#'
#' @returns A dataframe of similar fluorophores or a ggplot2 object
#' @export
#'
#' @examples
#' library(dplyr)
#' Folder_Location <- system.file("extdata", package = "Luciernaga")
#' XML_Pattern <- ".XML$"
#' XML_Files <- list.files(path = Folder_Location, pattern = XML_Pattern,
#'                         full.names = TRUE, recursive = FALSE)
#' Data <- QC_LibraryParse(XML_Files[2], returntype="data", references=FALSE)
#' TheFluorophore <- Data |> pull(Sample)
#'
#' Results <- QC_WhatsThis(x=TheFluorophore, columnname="Sample", data=Data, NumberHits = 10, returnPlots=FALSE)

QC_WhatsThis <- function(x, columnname="Sample", data, NumberHits, NumberDetectors=NULL,
 Normalize=TRUE, returnPlots=FALSE, TheFormat="wider", detectorcolumn=NULL, valuecolumn=NULL) {

  StartingData <- data |> filter(.data[[columnname]] %in% x)

  if(nrow(StartingData) > 1){message("Selecting the first row for comparisons")
                             StartingData <- StartingData |> slice(1)}

  CharacterLength <- StartingData |> select(!where(is.numeric))|> length()

  if (CharacterLength == 0){
    stop("Please add a non-numeric column, and provide its columnname")}
  if (CharacterLength > 1){message("Combining character columns")
    Identity <- StartingData |> select(!where(is.numeric)) |>
      paste0(collapse = "_")
    Identity <- data.frame(Fluorophore=Identity)
    Identity <- Identity |> rename("Fluorophore"=1)
    Identity <- Identity |> mutate(Fluorophore=paste0("ID_", Fluorophore))
  } else {
    Identity <- StartingData |> select(!where(is.numeric)) |> rename("Fluorophore"=1)
    Identity <- Identity |> mutate(Fluorophore=paste0("ID_", Fluorophore))
  }

  DetectorCols <- StartingData |> select(where(is.numeric))

  if (Normalize == TRUE){
    if (any(DetectorCols > 1)){
      message("Normalizing Data for Signature Comparison")
      n <- DetectorCols
      n[n < 0] <- 0
      A <- do.call(pmax, n)
      Normalized <- n/A
      DetectorCols <- Normalized
    }
  }
  
  WhoseThis <- cbind(Identity, DetectorCols)

  if (is.null(NumberDetectors)){NumberDetectors <- length(DetectorCols)
  } else {NumberDetectors <- NumberDetectors}


  ReferenceData <- Luciernaga:::InstrumentReferences(NumberDetectors=NumberDetectors)

  # Longer-format
  if (returnPlots == TRUE){
    ReferenceData1 <- ReferenceData |> select(-Instrument)
    WhoseThis1 <- WhoseThis |> pivot_longer(
      cols= where(is.numeric), names_to = "Detector", values_to = "AdjustedY")
    ReferenceData1 <- bind_rows(WhoseThis1, ReferenceData1)
  }

  # Wider-format
  ReferenceData <- ReferenceData |> select(-Instrument) |>
    group_by(Fluorophore)|> pivot_wider(
      names_from = Detector, values_from = AdjustedY)|> ungroup()
  if (length(colnames(ReferenceData)) == length(colnames(WhoseThis))){
      colnames(ReferenceData) <- colnames(WhoseThis)
  } else {warning("ReferenceData and WhoseThis have differing dimensions")}
  CombinedView <- bind_rows(WhoseThis, ReferenceData)

  # Generating Cosine Matrix
  Names <- CombinedView |> select(Fluorophore) |> pull()
  Numbers <- CombinedView |> select(where(is.numeric)) 
  NumericsT <- t(Numbers)
  rownames(NumericsT) <- NULL
  colnames(NumericsT) <- Names
  CosineMatrix <- cosine(NumericsT)
  CosineMatrix <- round(CosineMatrix, 2)
  CosineFrame <- data.frame(CosineMatrix, check.names = FALSE)

  #CosineFrame <- CosineFrame[1,] #If Want To Work With Cols
  CosineFrame <- CosineFrame |> select(starts_with("ID_"))
  TheData <- rownames_to_column(CosineFrame, var="Fluorophore")
  TheID <- TheData |> select(starts_with("ID_")) |> colnames()
  TheHits <- TheData |> filter(!Fluorophore %in% TheID) |>
    arrange(desc(.data[[TheID]])) |> slice_head(n=NumberHits)

  if (returnPlots==TRUE){
    TheseFluorophores <- TheHits |> pull(Fluorophore)
    TheFluorophore <- TheID
    ThePlot <- Luciernaga:::SimilarFluorPlots(TheseFluorophores=TheseFluorophores,
                                 TheFluorophore=TheFluorophore, data=ReferenceData1)
    ReturnThese <- list(TheHits, ThePlot)
    return(ReturnThese)
  } else {return(TheHits)}
}
