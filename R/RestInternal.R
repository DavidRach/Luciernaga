#' Internal for Wetlab Rest
#'
#' @param x Iterated specimen name
#' @param Updated The data.frame
#' @param DesiredConcentration_MillionperML Passed parameter
#' @param TubeMaxML Passed parameter
#' @param DesiredConcentration Redundant passed parameter
#'
#' @importFrom dplyr filter select pull mutate relocate case_when %>%
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
    NewName <- SpinProtocol |> pull(name) %>% paste0("Spin_", .)
    SpinProtocol <- SpinProtocol |> mutate(name = case_when(
      is.character(name) ~ NewName, TRUE ~ name))
    SpinProtocol$NeededVolume <- format(SpinProtocol$NeededVolume, digits=2)
    SpinProtocol$IncreaseVolumeML <- format(SpinProtocol$IncreaseVolumeML, digits=2)
    SpinProtocol$TotalTubes <- format(SpinProtocol$TotalTubes, digits=2)
    Instructions <- rbind(Instructions, SpinProtocol)
  }
  return(Instructions)
}