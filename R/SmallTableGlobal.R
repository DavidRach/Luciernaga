#' Dashboard Internal, fills color code for global view
#'
#' @param data Assembled data for cytometer pass fails past three months
#'
#' @importFrom gt gt
#' @importFrom gt data_color
#' @importFrom tidyselect everything
#' @importFrom scales col_factor
#' @importFrom gt sub_values
#' @importFrom gt opt_table_font
#' @importFrom gt cols_align
#' @noRd
SmallTableGlobal <- function(data){
  table <- data %>%
    gt() %>%
    data_color(
      columns = c(everything()),
      fn=function(x) {
        x <- as.character(x)

        dplyr::case_when(
          x == "Green" ~ "#0B6623",
          x == "Orange" ~ "#BA8E23",
          x == "Yellow" ~ "#BA8E23",
          x == "Red" ~ "#C80815",
          is.na(x) ~ "#ECECEC",
          TRUE ~ "#FFFFFF"
        )
      }
    )

  Substituted <- table  |>
    sub_values(values= c("Green"), replacement = "Pass") |>
    sub_values(values= c("Orange"), replacement = "Caution") |>
    sub_values(values= c("Yellow"), replacement = "Caution") |>
    sub_values(values= c("Red"), replacement = "Fail")

  Bolded <- Substituted |>
    opt_table_font(font = "Montserrat") |>
    cols_align(align = "center")

  Final <- Bolded

  return(Final)
}