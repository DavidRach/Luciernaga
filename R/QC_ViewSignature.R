#' Visualizes the Signature for given row in an averaged signature
#'  data.frame.
#'
#' @param x Name in the Sample column you want to filter for
#' @param columnname Default is Sample, specifies column name
#'  from which x is filtered from 
#' @param data A data.frame object from QC_LibraryParse containing
#'  Fluorophore name column 
#' and numeric detector columns.
#' @param Normalize Whether to normalize the data based on peak
#'  detector value, default is TRUE
#' @param TheFormat Default wider for detectors in columns, specify
#'  longer if providing detectors as rows
#' @param detectorcolumn Default NULL, when TheFormat="longer" specify
#'  detector column name
#' @param valuecolumn Default NULL, when TheFormat="longer" specify
#'  value column name
#' @param legend Default TRUE, alternately removes plot legend
#' @param plotname Default NULL, alternately specifies the plot
#'  title
#' @param plotlinecolor Default NULL, alternatively provide color
#'  when only a single line
#'
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom tidyselect everything
#' @importFrom tidyr unite
#' @importFrom tidyselect where
#' @importFrom dplyr rename
#' @importFrom dplyr mutate
#' @importFrom dplyr pull
#' @importFrom tidyr pivot_longer
#'
#' @returns A dataframe of similar fluorophores
#' @export
#'
#' @examples
#'
#' library(flowCore)
#' library(flowWorkspace)
#' library(openCyto)
#' library(data.table)
#' library(dplyr)
#'
#' File_Location <- system.file("extdata", package = "Luciernaga")
#' FCS_Files <- list.files(path = File_Location, pattern = ".fcs",
#'   full.names = TRUE)
#' UnstainedFCSFiles <- FCS_Files[grep("Unstained", FCS_Files)]
#' UnstainedCells <- UnstainedFCSFiles[-grep("Beads", UnstainedFCSFiles)]
#' MyCytoSet <- load_cytoset_from_fcs(UnstainedCells[1],
#'   truncate_max_range = FALSE,transformation = FALSE)
#' MyGatingSet <- GatingSet(MyCytoSet)
#' MyGates <- fread(file.path(path = File_Location, pattern = 'Gates.csv'))
#' MyGatingTemplate <- gatingTemplate(MyGates)
#' gt_gating(MyGatingTemplate, MyGatingSet)
#' removestrings <-  c("DR_", "Cells", ".fcs", "-", " ")
#' StorageLocation <- file.path("C:", "Users", "JohnDoe", "Desktop")
#'
#' PopulationInterest <- gs_pop_get_data(MyGatingSet[1], subset="lymphocytes")
#' TheDataValues <- exprs(PopulationInterest[[1]])
#' TheDataValues <- data.frame(TheDataValues, check.names=FALSE)
#' Signature <- AveragedSignature(TheDataValues, stats="median")
#' TheData <- Signature[,-grep("Time|FS|SC|SS|Original|W$|H$", names(Signature))]
#' TheData <- TheData |> mutate(Sample="TestSignature") |>
#'  relocate(Sample, .before=1)
#'
#' Plot <- QC_ViewSignature(x="TestSignature", data=TheData, Normalize=TRUE)
#'
QC_ViewSignature <- function(x, columnname="Sample", data, Normalize = TRUE,
 TheFormat="wider", detectorcolumn=NULL, valuecolumn=NULL,
 legend=TRUE, plotname=NULL, plotlinecolor=NULL) {

  if (is.null(x)){x <- data |> dplyr::pull(columnname)
  }

  if (TheFormat=="wider"){
  StartingData <- data |> filter(.data[[columnname]] %in% x)

  CharacterLength <- StartingData |> select(!where(is.numeric)) |> length()
  if (CharacterLength == 0){
    stop("Please add a non-numeric column, and provide its columnname")}
  if (CharacterLength > 1){message("Combining character columns")
    Identity <- StartingData |> select(!where(is.numeric)) |>
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
      n[n < 0] <- 0
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
    
  } else {
    StartingData <- data |> filter(.data[[columnname]] %in% x)
    StartingData <- StartingData |> rename(Fluorophore=columnname)
    StartingData <- StartingData |> mutate(Fluorophore=paste0("ID_", Fluorophore))
    TheseFluorophores <- StartingData |> pull(Fluorophore) |> unique()
    StartingData <- StartingData |> rename(Detector=detectorcolumn)
    StartingData <- StartingData |> rename(AdjustedY=valuecolumn)
    WhoseThis1 <- StartingData
  
  } 

  ThePlot <- SimilarFluorPlots(TheseFluorophores=TheseFluorophores,
                                TheFluorophore=NULL, data=WhoseThis1,
                                legend=legend, plotname=plotname,
                                plotlinecolor=plotlinecolor)

  return(ThePlot)
  }
