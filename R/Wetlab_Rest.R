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


#' Internal for Wetlab Rest
#'
#' @param x Iterated specimen name
#' @param Updated The data.frame
#' @param DesiredConcentration_MillionperML Passed parameter
#' @param TubeMaxML Passed parameter
#' @param DesiredConcentration Redundant passed parameter
#'
#' @importFrom dplyr filter select pull mutate relocate case_when
#'
#' @return An internal value
#'
#' @noRd
RestInternal <- function(x, Updated, DesiredConcentration_MillionperML,
   TubeMaxML, DesiredConcentration){
  Internal <- Updated |>  filter(name %in% x)
  Name <- Internal |> select(name)
  Date <- Internal |> select(Date)
  TotalCells <- Internal |> pull(TotalCells) |> as.numeric()
  CurrentConcentration <- Internal |> select(ConcentrationScientific) |> pull()

  TotalVolume <- Internal |> pull(TotalVolume) |> as.numeric()
  A <- TotalCells/(DesiredConcentration_MillionperML*1000000)
  IncreaseVolumeML <- A-TotalVolume
  IncreaseVolumeML <- format(IncreaseVolumeML, digits=2)
  TotalTubes <- A/TubeMaxML
  TotalTubes <- format(TotalTubes, digits=2)
  TotalCells <- format(TotalCells, scientific = TRUE, digits=2)
  NeededVolume <- format(A, digits=2)

  Instructions <- cbind(Name, Date, CurrentConcentration, TotalVolume, TotalCells,
     DesiredConcentration, NeededVolume, IncreaseVolumeML, TubeMaxML, TotalTubes)
  Instructions <- Instructions |>
    mutate(SpinDown=ifelse(IncreaseVolumeML < 0, TRUE, FALSE)) |>
    relocate(SpinDown, .before=TubeMaxML)

  if(Instructions$SpinDown == TRUE){
    SpinProtocol <- Instructions
    SpinProtocol <- SpinProtocol |> mutate(TotalVolume = case_when(
      TotalVolume > 0 ~ NA_real_,TRUE ~ TotalVolume))
    SpinProtocol <- SpinProtocol |> mutate(CurrentConcentration = case_when(
      CurrentConcentration > 0 ~ NA_real_,TRUE ~ CurrentConcentration))
    TotalCells <- SpinProtocol |> pull(TotalCells) |> as.double()
    SpinProtocol$NeededVolume <- as.double(SpinProtocol$NeededVolume)
    SpinProtocol$IncreaseVolumeML <- as.double(SpinProtocol$IncreaseVolumeML)
    SpinProtocol$TubeMaxML <-  as.double(SpinProtocol$TubeMaxML)
    SpinProtocol$TotalTubes <-  as.double(SpinProtocol$TotalTubes)
    A <- TotalCells/(DesiredConcentration_MillionperML*1000000)
    SpinProtocol <- SpinProtocol |> mutate(NeededVolume = case_when(
      NeededVolume > 0 ~ A, TRUE ~ NeededVolume))
    SpinProtocol <- SpinProtocol |> mutate(IncreaseVolumeML = case_when(
      is.double(IncreaseVolumeML) ~ A, TRUE ~ IncreaseVolumeML))
    if (SpinProtocol$IncreaseVolumeML > 0){
      SpinProtocol <- SpinProtocol |> mutate(SpinDown = case_when(
        is.logical(SpinDown) ~ FALSE, TRUE ~ SpinDown))
    }
    TotalTubesTwo <- A/TubeMaxML
    SpinProtocol <- SpinProtocol |> mutate(TotalTubes = case_when(
      is.double(TotalTubes) ~ TotalTubesTwo, TRUE ~ TotalTubes))
    NewName <- SpinProtocol |> pull(name) |> paste0("Spin_", .)
    SpinProtocol <- SpinProtocol |> mutate(name = case_when(
      is.character(name) ~ NewName, TRUE ~ name))
    SpinProtocol$NeededVolume <- format(SpinProtocol$NeededVolume, digits=2)
    SpinProtocol$IncreaseVolumeML <- format(SpinProtocol$IncreaseVolumeML, digits=2)
    SpinProtocol$TotalTubes <- format(SpinProtocol$TotalTubes, digits=2)
    Instructions <- rbind(Instructions, SpinProtocol)
  }
  return(Instructions)
}


#' Internal for Wetlab_Rest
#'
#' @param data The resuspension data from WetlabRest
#' @param outpath The desired storage location for the .png file
#' @param filename The desired name of the .png file
#' @param ColorSelection Default NULL
#' @param vwidth Default 1200
#' @param outputType Default png, alternate pdf
#'
#' @importFrom gt cells_body gt tab_style cell_fill cell_text opt_table_font
#' cols_label cols_align opt_table_outline tab_options gtsave px html pct
#' @importFrom rlang sym
#' @importFrom dplyr arrange desc slice pull
#' @importFrom ghibli ghibli_palette
#'
#' @return An internal value
#'
#' @noRd
RestTable <- function(data, outpath=NULL, filename="CellResuspensions", ColorSelection,
  vwidth=1200, outputType){

  if (is.null(ColorSelection)){
  Palette <- ghibli_palette(name="PonyoLight", n=7, direction=1, type="discrete")
  ColorSelection <- Palette[5]
  }

  builder <- function(x, Limit){
    cells_body(columns = !!sym(x), rows = !!sym(x) >= Limit)
  }

  revbuilder <- function(x, Limit){
    cells_body(columns = !!sym(x), rows = !!sym(x) < Limit)
  }

  betweenbuilder <- function(x, Limit1, Limit2){
    cells_body(columns = !!sym(x), rows = !!sym(x) >= Limit1 & !!sym(x) < Limit2)
  }

  Date <- data.frame(table(data$Date))
  Date <- Date |> arrange(desc(Freq)) |> slice(1) |>
    pull(Var1) |> as.character()

  Table <- data |> gt() |> 
    tab_style(
      style = cell_fill(color = "#e8f5e9"),
      locations = cells_body(
        rows = str_detect(name, "_00_"),
        columns = c(TotalCells, TotalVolume, DesiredConcentration, TubeMaxML, 
                   NeededVolume, CurrentConcentration)
      )
    ) |>
    tab_style(
      style = cell_fill(color = "#e3f2fd"),
      locations = cells_body(
        rows = str_detect(name, "_09_"),
        columns = c(TotalCells, TotalVolume, DesiredConcentration, TubeMaxML, 
                   NeededVolume, CurrentConcentration) 
      )
    ) |>
    tab_style(style = list(
      cell_fill(color = ColorSelection),
      cell_text(weight = "bold")
    ), locations = cells_body(
      columns = c(IncreaseVolumeML, TotalTubes)
    )
    ) |> tab_style(style = cell_text(weight = "bold"), locations = cells_body(
      columns = c(name, IncreaseVolumeML, TotalTubes)
    )
    )

  Bolded <- Table |> opt_table_font(font = "Montserrat") |>
    cols_label(name ~ "Specimen") |>
    cols_label(TotalCells ~ html("Total<br> Cells<br>")) |>
    cols_label(TotalVolume ~ html("Total<br> Volume<br>")) |>
    cols_label(DesiredConcentration ~ html("Desired<br> Concentration<br>")) |>
    cols_label(TubeMaxML ~ html("Max mLs<br> per Tube<br>")) |>
    cols_label(IncreaseVolumeML ~ html("Increase<br> Volume mLs<br>")) |>
    cols_label(TotalTubes ~ html("Total<br> Tubes<br>")) |>
    cols_label(NeededVolume ~ html("Required<br> Volume mLs<br>")) |>
    cols_label(CurrentConcentration ~ html("Current<br> Concentration<br>")) |>
    cols_align(align = "center")

  FinalTable <- Bolded |> opt_table_outline(style="solid", width=px(2), color="black") |>
    tab_options(column_labels.border.top.style = "solid",
                column_labels.border.top.width = px(2),
                column_labels.border.top.color = "black",
                table.border.bottom.style = "double",
                table.border.bottom.width = px(2),
                table.border.bottom.color = "black",
                table.border.right.style = "solid",
                table.border.right.width = px(2),
                table.border.right.color = "black",
                table.width = pct(100))
  
  if (is.null(outpath)){outpath <- getwd()}
  
  
  if (outputType == "pdf") {
    if (!grepl("\\.pdf$", filename, ignore.case = TRUE)) {
      filename <- paste0(filename, ".pdf")
    }
  } else {
    if (!grepl("\\.png$", filename, ignore.case = TRUE)) {
      filename <- paste0(filename, ".png")
    }
  }
  
  TheFile <- file.path(outpath, filename)

  if (outputType == "png") {
    gtsave(FinalTable, TheFile, vwidth = vwidth)
  }

}


