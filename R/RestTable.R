
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
#'  cols_label cols_align opt_table_outline tab_options gtsave px html pct
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
