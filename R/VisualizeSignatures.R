
#' Visualizes the provided fluorophore signatures from a template data.frame
#'
#' @param x Default NULL, otherwise name for a value in the columnname column
#' that you want to individually visualize
#' @param columnname Default is Flurophore, specifies column name
#'  from which x is filtered from 
#' @param data A template data.frame containing the fluorophore signature
#' @param Normalize Whether to normalize the data based on peak
#'  detector value, default is TRUE
#' @param legend Default TRUE, alternately removes plot legend
#' @param plotname Default NULL, alternately specifies the plot
#'  title
#' @param plotlinecolor Default NULL, alternatively provide color
#'  when only a single line
#'
#' @importFrom dplyr filter select rename mutate pull
#' @importFrom tidyselect everything where
#' @importFrom tidyr unite pivot_longer

#'
#' @returns Visualized plot
#' @export
#'
#' @examples A <- 2 + 2
#'
VisualizeSignatures <- function(x=NULL, columnname="Sample", data, 
 characterColumns=NULL, Normalize = TRUE,
 legend=TRUE, plotname=NULL, plotlinecolor=NULL) {

  if (is.null(x)){x <- data |> pull(columnname)}

  StartingData <- data |> filter(.data[[columnname]] %in% x)
  
  CharacterLength <- StartingData |> select(!where(is.numeric)) |> length()

  if (CharacterLength == 0){
    stop("Please add a non-numeric column, and provide its columnname")
  }
    
  if (CharacterLength > 1 && !is.null(characterColumns)){
    Identity <- StartingData |> select(all_of(characterColumns)) |>
      unite("combined", everything(), sep = "_") |> pull()
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
      # n[n < 0] <- 0
      A <- do.call(pmax, n)
      Normalized <- n/A
      DetectorCols <- Normalized
    }
  }
    
  WhoseThis <- cbind(Identity, DetectorCols)
  TotalDetectors <- length(DetectorCols)
  TheseFluorophores <- WhoseThis |> pull(Fluorophore)
    
  WhoseThis1 <- WhoseThis |>
    pivot_longer(cols= where(is.numeric), names_to = "Detector",
                 values_to = "AdjustedY") 

  ThePlot <- SimilarFluorPlots(TheseFluorophores=TheseFluorophores,
                                TheFluorophore=NULL, data=WhoseThis1,
                                legend=legend, plotname=plotname,
                                plotlinecolor=plotlinecolor)

  return(ThePlot)
  }