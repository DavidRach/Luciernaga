#' Visualizes the Signature for given row in an averaged signature data.frame.
#'
#' @param x Name in the Sample column you want to filter for
#' @param data A data.frame object from QC_LibraryParse containing Fluorophore name column and numeric detector columns.
#' @param Normalize Whether to normalize the data based on peak detector value, default is TRUE
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

QC_ViewSignature <- function(x, data, Normalize = TRUE) {

  StartingData <- data %>% filter(Sample %in% x)

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

  TotalDetectors <- length(DetectorCols)

  WhoseThis1 <- WhoseThis %>%
    pivot_longer(cols= where(is.numeric), names_to = "Detector",
                 values_to = "AdjustedY")

  TheseFluorophores <- WhoseThis %>% pull(Fluorophore)

  ThePlot <- SimilarFluorPlots(TheseFluorophores=TheseFluorophores,
                                 TheFluorophore=NULL, data=WhoseThis1)

  return(ThePlot)
}
