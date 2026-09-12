#' Takes Concentration output, returns re-suspension amounts for Wetlab users.
#'
#' @param data The Wetlab_Concentration output (modified correct Total mL amounts)
#' @param DesiredConcentration_MillionperML Final desired concentration (ex. 1, 3, etc.)
#' @param MaxMLperTube The Final Volume of the rest tube
#' @param returntype Whether to return "data", "plot" or "both"
#' @param outpath File path of location to store the table
#' @param filename File name to store the table as, ex
#' @param ColorSelection Default NULL, else provide desired color palette. 
#' @param outputType Default png, alternate pdf
#'
#' @importFrom dplyr select mutate relocate bind_rows
#' @importFrom tidyselect any_of
#' @importFrom purrr map
#'
#' @return A data.frame of resuspension measurements to get to the desired parameters
#' @export
#'
#' @examples
#'
#' library(flowCore)
#' library(flowWorkspace)
#' library(CytoML)
#' library(dplyr)
#' library(purrr)
#'
#' File_Location <- system.file("extdata", package = "Luciernaga")
#' WSP_File <- list.files(File_Location, pattern=".wsp", full.names = TRUE)
#' ws <- open_flowjo_xml(WSP_File[1])
#' gs <- flowjo_to_gatingset(ws, name=1, path = File_Location)
#' nameKeyword <- c("GROUPNAME", "TUBENAME")
#'
#' TheData <- map(.x=gs, Wetlab_Concentration, subset = "CD45+",
#'   nameKeyword=nameKeyword, DilutionMultiplier=100, TotalVolume=1) %>%
#'    bind_rows()
#'
#' UpdatedData <- TheData %>% select(-TotalScientific, -TimeSeconds)
#'
#' Results <- Wetlab_Rest(data=UpdatedData, DesiredConcentration_MillionperML=3,
#'  MaxMLperTube=1, returntype="data", outpath=path)
#'
Wetlab_Rest <- function(data, DesiredConcentration_MillionperML, MaxMLperTube, returntype,
                        outpath=NULL, filename="CellResuspensions",
                        ColorSelection=NULL, outputType="png"){

  TheColNames <- colnames(data)
  RemoveThese <- c("TimeSeconds", "TotalScientific")
  if (any(RemoveThese %in% TheColNames)) {data <- data |> select(-any_of(RemoveThese))}

  data$Cells <- as.numeric(data$Cells)
  data$Volume <- as.numeric(data$Volume)
  data$ConcentrationScientific <- as.numeric(data$ConcentrationScientific)

  Updated <- data |> mutate(TotalCells=ConcentrationScientific*TotalVolume) |>
    relocate(TotalCells, .before=Instrument)
  Updated$TotalCells <- format(Updated$TotalCells, scientific = TRUE, digits = 2)
  Specimens <- Updated$name

  DesiredConcentration <- DesiredConcentration_MillionperML*1000000
  DesiredConcentration <- format(DesiredConcentration, scientific = TRUE, digits = 2)
  TubeMaxML <- MaxMLperTube

  Instructions <- map(.x=Specimens, .f=RestInternal, Updated=Updated,
                      DesiredConcentration_MillionperML=DesiredConcentration_MillionperML,
                      TubeMaxML=TubeMaxML, DesiredConcentration=DesiredConcentration) |> bind_rows()

  if (returntype == "data"){return(Instructions)
  } else if (returntype == "plot"){
    plot <- RestTable(data=Instructions, outpath=outpath, filename=filename,
       ColorSelection=ColorSelection, outputType=outputType)
  } else if (returntype == "both"){
    plot <- RestTable(data=Instructions, outpath=outpath, filename=filename,
    outputType=outputType)
    return(Instructions)
  }
}