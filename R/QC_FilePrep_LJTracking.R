#' Converts the Cytek Aurora (TM)'s QC report into a data frame.
#'
#' @param x  Takes a Levy-Jennings QC tracking report saved as a .csv file, and
#' converts into a "tidyed" dataframe for plotting. Currently works on our 3L, 4L, 5L
#' Auroras. Please reach out if you find an issue, the .csv export varies a bit and
#' I want to continue to improve on the code to handle these odd exceptions.
#' @param DailyQC A single DailyQCReport .csv file, used to import baseline settings.
#'
#' @importFrom purrr map2
#' @importFrom tidyr pivot_longer
#' @importFrom tidyr pivot_wider
#' @importFrom tidyselect all_of
#' @importFrom dplyr "%>%"
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr left_join
#' @importFrom dplyr case_when
#' @importFrom dplyr bind_cols
#' @importFrom dplyr across
#' @importFrom tidyselect everything
#' @importFrom lubridate mdy_hms
#'
#' @return A dataframe.
#' @export
#'
#' @examples
#' File_Location <- system.file("extdata", package = "Luciernaga")
#' CSV_Pattern <- ".CSV$"
#' CSV_Files <- list.files(path=File_Location, pattern=CSV_Pattern,
#'                        full.names=TRUE)
#' TidyData <- QC_FilePrep(x=CSV_Files, DailyQC=DailyQC)

QC_FilePrep_LJTracking <- function(x, DailyQC){
  ReadInfo <- readLines(x)
  ReadInfo <- ReadInfo[ReadInfo != ""]

  GainIndex <- grep("^Gain", ReadInfo)
  rCVIndex <- grep("^% rCV", ReadInfo)
  LaserIndex <- grep("^Laser Delay", ReadInfo)
  AreaIndex <- grep("^Area Scaling Factor", ReadInfo)
  FSCIndex <- grep("^FSC Area Scaling Factor", ReadInfo)
  LaserPowerIndex <- grep("^Laser Power", ReadInfo)

  ThePositions <- c(GainIndex, rCVIndex, LaserIndex,
                    AreaIndex, FSCIndex, LaserPowerIndex)

  Final <- length(ReadInfo)
  StartPositions <- ThePositions + 1
  EndPositions <- ThePositions - 1
  EndPositions <- EndPositions[-1]
  EndPosition <- c(EndPositions, Final)

  if (StartPositions[length(StartPositions)] > Final){
    StartHere <- StartPositions[-length(StartPositions)]
    EndHere <- EndPosition[-length(EndPosition)]
  } else {
    StartHere <- StartPositions
    EndHere <- EndPosition
  }

  TheFrames <- map2(.x=StartHere[3], .y=EndHere[3], ReadInfo=ReadInfo, .f=ParseThis)
  #TheFrames <- map2(.x=StartHere, .y=EndHere, ReadInfo=ReadInfo, .f=ParseThis)

  Updated <- TheFrames

  Gains <- Updated[[1]]
  GainLength <- ncol(Gains)

  GainLonger <- Gains %>%
    pivot_longer(all_of(2:GainLength), names_to="DetectorMain", values_to="Value")

  Baseline <- DailyQCBaseline(DailyQC)
  Comparison <- cbind(GainLonger, Baseline)
  Comparison$Value <- as.numeric(Comparison$Value)
  Comparison$Comparison <- as.numeric(Comparison$Comparison)

  Comparison <- Comparison %>%
    mutate(GainFlag = Value > Comparison) %>%
    select(-Detector, -Comparison) %>% mutate(DetectorFlag=paste0("Flag-", DetectorMain))

  GainsMain <- Comparison %>% select(DateTime, DetectorMain, Value) %>%
    pivot_wider(names_from = DetectorMain, values_from = Value)

  GainsFlag <- Comparison %>% select(DateTime, DetectorFlag, GainFlag) %>%
    pivot_wider(names_from = DetectorFlag, values_from = GainFlag)

  TheGains <- left_join(GainsMain, GainsFlag, by="DateTime")

  RCVs <- Updated[[2]]
  RCVsLength <- ncol(RCVs)
  RCVsLonger <- RCVs %>%
    pivot_longer(all_of(2:RCVsLength), names_to = "RCVMain", values_to = "Value")

  RCVsLonger$Value <- as.numeric(RCVsLonger$Value)

  Comparison <- RCVsLonger %>% mutate(RCVFlag = Value > 6) %>%
    mutate(Detector=gsub("-% rCV", "", RCVMain)) %>%
    mutate(RCVFlag = case_when(
      Detector == "SSC" | Detector == "SSC-B" & Value < 8 ~ FALSE,
      TRUE ~ RCVFlag)) %>% select(-Detector) %>%
      mutate(FlagRCV=paste0("Flag-", RCVMain))

  RCVsMain <- Comparison %>% select(DateTime, RCVMain, Value) %>%
    pivot_wider(names_from=RCVMain, values_from = Value)

  RCVsFlag <- Comparison %>% select(DateTime, FlagRCV, RCVFlag) %>%
    pivot_wider(names_from=FlagRCV, values_from = RCVFlag)

  TheRCVs <- left_join(RCVsMain, RCVsFlag, by="DateTime")

  UpdatedItems <- length(Updated)

  Others <- Updated[3:UpdatedItems] %>%
    lapply(function(df) df %>% select(-DateTime)) %>%
    bind_cols()

  Others[] <- lapply(Others, as.numeric)

  OtherAdditional <- Others %>%
    mutate(across(everything(), ~ FALSE, .names = "Flag-{.col}"))

  TheDataset <- left_join(TheGains, TheRCVs, by="DateTime")
  TheDataset <- cbind(TheDataset, OtherAdditional)

  TheDataset$DateTime <- mdy_hms(TheDataset$DateTime)

  TheDataset$DateTime <-
  return(TheDataset)
}