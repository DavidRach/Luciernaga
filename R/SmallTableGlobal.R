#' Dashboard Internal, fills color code for global view
#'
#' @param data Assembled data for cytometer pass fails past three months
#'
#' @importFrom gt gt data_color sub_values opt_table_font cols_align
#' @importFrom tidyselect everything
#' @importFrom scales col_factor
#' @importFrom dplyr case_when
#' 
#' @return A gt table object with color-coded, bolded, and aligned cells
#' 
#' @noRd
SmallTableGlobal <- function(data) {
  table <- data |>
    gt() |>
    data_color(
      columns = c(everything()),
      fn = function(x) {
        x <- as.character(x)

        case_when(
          x == "Green" ~ "#0B6623",
          x == "Orange" ~ "#BA8E23",
          x == "Yellow" ~ "#BA8E23",
          x == "Red" ~ "#C80815",
          is.na(x) ~ "#ECECEC",
          TRUE ~ "#FFFFFF"
        )
      }
    )

  Substituted <- table |>
    sub_values(values = c("Green"), replacement = "Pass") |>
    sub_values(values = c("Orange"), replacement = "Caution") |>
    sub_values(values = c("Yellow"), replacement = "Caution") |>
    sub_values(values = c("Red"), replacement = "Fail")

  Bolded <- Substituted |>
    opt_table_font(font = "Montserrat") |>
    cols_align(align = "center")

  Final <- Bolded

  return(Final)
}