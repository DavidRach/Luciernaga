#' Dashboard Internal, returns gt table
#'
#' @param data The QC color status returns
#'
#' @importFrom gt gt data_color sub_values opt_table_font 
#'  cols_align tab_spanner cols_label
#' @importFrom scales col_factor
#' 
#' @noRd
SmallTable <- function(data){

  table <- data %>%
    gt() %>%
    data_color(
      columns = c(Gain, rCV),
      fn=function(x) {
        dplyr::case_when(
          x == "Green" ~ "#0B6623",
          x == "Orange" ~ "#FF6E00",
          x == "Yellow" ~ "#BA8E23",
          x == "Red" ~ "#C80815",
          x == "Gray" ~ "#D3D3D3",
          TRUE ~ NA_character_
        )
      }
  )

    Substituted <- table  |>
      sub_values(values= c("Green"), replacement = "Pass") |>
      sub_values(values= c("Orange"), replacement = "Warning") |>
      sub_values(values= c("Yellow"), replacement = "Caution") |>
      sub_values(values= c("Red"), replacement = "Fail") |>
      sub_values(values= c("Gray"), replacement = "")

    Bolded <- Substituted |>
      opt_table_font(font = "Montserrat") |>
      cols_align(align = "center")

    Final <- Bolded |> tab_spanner(
      label = "Gain ",
      columns = c(GainValue, Gain)
    ) |> tab_spanner(
      label = "%RCV ",
      columns = c(rCVValue, rCV)
    ) |> cols_label(
      GainValue = "Value",
      Gain = "Status"
    ) |> cols_label(
      rCVValue = "Value",
      rCV = "Status"
    )

  return(Final)
}