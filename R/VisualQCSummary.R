#' Dashboard Internal, processes to did parameter pass in past week
#'
#' @param x The data.frame output from LevyJennings or QCBeads Parse
#' @param detectorType Default is "-A"
#'
#' @importFrom dplyr filter select left_join pull bind_rows
#' @importFrom lubridate weeks
#' @importFrom tidyselect starts_with contains all_of
#' @importFrom tidyr pivot_longer
#' @importFrom purrr map
#' @importFrom stringr str_detect
#'
#' @return Data frame of passing status for respective parameters
#' @noRd
VisualQCSummary <- function(x, detectorType="-A"){

  WindowOfInterest <- Sys.time() - weeks(1)

  if (nrow(x) > 1){
  Data <- x |> filter(DateTime > WindowOfInterest)

  if (nrow(Data) == 0){Data <- x |> slice(1)}

  } else {Data <- x}

  if (any(stringr::str_detect(colnames(Data), "Flag"))){
    Flags <- Data |> select(starts_with("Flag"))
    colnames(Flags) <- gsub("Flag-", "", colnames(Flags))
    Gains <- Flags |> select(contains("Gain"))
    TheGains <- colnames(Gains)
    rCV <- Flags |> select(contains("rCV"))
    TherCV <- colnames(rCV)

    colnames(Gains) <- gsub("-Gain", "", fixed=TRUE, colnames(Gains))
    colnames(Gains) <- gsub("_Gain", "", fixed=TRUE, colnames(Gains))
    colnames(rCV) <- gsub("-% rCV", "", fixed=TRUE, colnames(rCV))
  } else {
    Gains <- Data |> select(contains("Gain"))
    TheGains <- colnames(Gains)
    rCV <- Data |> select(contains("rCV"))
    TherCV <- colnames(rCV) 

    colnames(Gains) <- gsub("-Gain", "", fixed=TRUE, colnames(Gains))
    colnames(Gains) <- gsub("_Gain", "", fixed=TRUE, colnames(Gains))
    Gains[] <- "NA"

    colnames(rCV) <- gsub("-% rCV", "", fixed=TRUE, colnames(rCV))
    rCV[] <- "NA"
  }

  TheGainData <- Data |> select(all_of(c("DateTime", TheGains)))
  colnames(TheGainData) <- gsub("_Gain", "", fixed=TRUE, colnames(TheGainData))
  colnames(TheGainData) <- gsub("-Gain", "", fixed=TRUE, colnames(TheGainData))

  TheGainData <- TheGainData |>
    pivot_longer(!DateTime, names_to = "Detector", values_to = "Gain")

  Gains <- Gains |> mutate(DateTime=Data$DateTime) |> relocate(DateTime, .before=1)
  Gains <- Gains |>
    pivot_longer(!DateTime, names_to = "Detector", values_to = "Gain_Logical")

  TherCVData <- Data |> select(all_of(c("DateTime", TherCV)))
  colnames(TherCVData) <- gsub("-% rCV", "", fixed=TRUE, colnames(TherCVData))
  TherCVData <- TherCVData |> pivot_longer(!DateTime, names_to = "Detector", values_to = "rCV")

  rCV <- rCV |> mutate(DateTime=Data$DateTime) |> relocate(DateTime, .before=1)
  rCV <- rCV |> pivot_longer(!DateTime, names_to = "Detector", values_to = "rCV_Logical")

  Tidy <- TheGainData |>
    left_join(Gains, by = c("Detector", "DateTime")) |>
    left_join(TherCVData, by = c("Detector", "DateTime")) |>
    left_join(rCV, by = c("Detector", "DateTime"))

  TheDetectors <- Tidy |> pull(Detector) |> unique()
  #TheDetectors <- TheDetectors[str_detect(TheDetectors, detectorType)]

  Summary <- map(.x=TheDetectors, .f=Luciernaga:::QCSummaryCheck, data=Tidy) |> bind_rows()
  return(Summary)
}