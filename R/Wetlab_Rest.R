#' Takes Concentration output, returns re-suspension amounts for Wetlab users.
#'
#' @param x A list of specimens to be mapped
#' @param data The Wetlab_Concentration output (modified correct Total mL amounts)
#' @param DesiredConcentration_MillionperML Final desired concentration (ex. 1, 3, etc.)
#' @param MaxMLperTube The Final Volume of the rest tube
#' @param returntype Whether to return "data", "plot" or "both"
#' @param outpath File path of location to store the table
#' @param filename File name to store the table as
#'
#' @importFrom dplyr select
#' @importFrom tidyselect any_of
#' @importFrom dplyr mutate
#' @importFrom dplyr relocate
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#'
#' @return A data.frame of resuspension measurements to get to the desired parameters
#' @export
#'
#' @examples NULL
Wetlab_Rest <- function(x, data, DesiredConcentration_MillionperML, MaxMLperTube, returntype,
                        outpath=NULL, filename="CellResuspensions"){

  TheColNames <- colnames(data)
  RemoveThese <- c("TimeSeconds", "TotalScientific")
  if (any(RemoveThese %in% TheColNames)) {data <- data %>% select(-any_of(RemoveThese))}

  Updated <- data %>% mutate(TotalCells=ConcentrationScientific*TotalVolume) %>%
    relocate(TotalCells, .before=Instrument)
  Updated$TotalCells <- format(Updated$TotalCells, scientific = TRUE, digits = 2)

  DesiredConcentration <- DesiredConcentration_MillionperML*1000000
  DesiredConcentration <- format(DesiredConcentration, scientific = TRUE, digits = 2)
  TubeMaxML <- MaxMLperTube

  Internal <- Updated %>%  dplyr::filter(name %in% x)
  Name <- Internal %>% select(name)
  Date <- Internal %>% select(Date)
  TotalCells <- Internal %>% select(TotalCells) %>% pull(.) %>% as.numeric()
  CurrentConcentration <- Internal %>% select(ConcentrationScientific) %>% pull(.)

  TotalVolume <- Internal %>% select(TotalVolume) %>% pull(.) %>% as.numeric()
  A <- TotalCells/(DesiredConcentration_MillionperML*1000000)
  IncreaseVolumeML <- A-TotalVolume
  IncreaseVolumeML <- format(IncreaseVolumeML, digits=2)
  TotalTubes <- A/MaxMLperTube
  TotalTubes <- format(TotalTubes, digits=2)
  TotalCells <- format(TotalCells, scientific = TRUE, digits=2)
  NeededVolume <- format(A, digits=2)

  Instructions <- cbind(Name, Date, CurrentConcentration, TotalVolume, TotalCells, DesiredConcentration, NeededVolume, IncreaseVolumeML, TubeMaxML, TotalTubes)
  Instructions <- Instructions %>% mutate(SpinDown=ifelse(IncreaseVolumeML < 0, TRUE, FALSE)) %>%
    relocate(SpinDown, .before=TubeMaxML)

  if (returntype == "data"){return(Instructions)
  } else if (returntype == "plot"){
    plot <- RestTable(data=Instructions, outpath=outpath, filename=filename)
  } else if (returntype == "both"){
    plot <- RestTable(data=Instructions, outpath=outpath, filename=filename)
    return(Instructions)
  }
}



#' Internal for Wetlab_Rest
#'
#' @param data The resuspension data from WetlabRest
#' @param outpath The desired storage location for the .png file
#' @param filename The desired name of the .png file
#'
#' @importFrom gt cells_body
#' @importFrom dplyr arrange
#' @importFrom dplyr select
#' @importFrom dplyr desc
#' @importFrom dplyr pull
#' @importFrom gt gt
#' @importFrom gt tab_style
RestTable <- function(data, outpath=NULL, filename="CellResuspensions"){

  message("Make sure to library(gt)")

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
  Date <- Date %>% arrange(desc(Freq)) %>% slice(1) %>%
    pull(Var1) %>% as.character(.)

  Table <- data |> gt() |>
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
                table.border.bottom.color = "black")

  filename <- paste0("CellResuspensions", "_", Date, ".png")
  TheFile <- file.path(outpath, filename)

  gtsave(FinalTable, TheFile)
}


