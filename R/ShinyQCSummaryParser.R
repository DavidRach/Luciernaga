#' Dashboard Internal, processes to did parameter pass in past week
#'
#' @param x The iterated date
#' @param Intermediate The original instrument data
#'
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom lubridate weeks
#' @importFrom tidyselect starts_with
#' @importFrom tidyselect contains
#' @importFrom tidyselect all_of
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr left_join
#' @importFrom dplyr pull
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
#'
#' @return Data frame of passing status for respective parameters
#' @noRd
ShinyQCSummaryParser <- function(x, Intermediate){
  Date <- x

  x <- Intermediate %>% filter(Date %in% x)

  if (nrow(x) > 1){
    #Data <- x %>% filter(DateTime > WindowOfInterest)
    Data <- x %>% slice(1)
  } else {Data <- x}

  Flags <- Data %>% select(starts_with("Flag"))
  colnames(Flags) <- gsub("Flag-", "", colnames(Flags))
  Gains <- Flags %>% select(contains("Gain"))
  TheGains <- colnames(Gains)
  rCV <- Flags %>% select(contains("rCV"))
  TherCV <- colnames(rCV)

  TheGainData <- Data %>% select(all_of(c("DateTime", TheGains)))
  colnames(Gains) <- gsub("-Gain", "", fixed=TRUE, colnames(Gains))
  colnames(TheGainData) <- gsub("-Gain", "", fixed=TRUE, colnames(TheGainData))

  TheGainData <- TheGainData %>%
    pivot_longer(!DateTime, names_to = "Detector", values_to = "Gain")

  Gains <- Gains %>% mutate(DateTime=Data$DateTime) %>% relocate(DateTime, .before=1)

  Gains <- Gains %>%
    pivot_longer(!DateTime, names_to = "Detector", values_to = "Gain_Logical")

  TherCVData <- Data %>% select(all_of(c("DateTime", TherCV)))
  colnames(rCV) <- gsub("-% rCV", "", fixed=TRUE, colnames(rCV))
  colnames(TherCVData) <- gsub("-% rCV", "", fixed=TRUE, colnames(TherCVData))

  TherCVData <- TherCVData %>% pivot_longer(!DateTime, names_to = "Detector", values_to = "rCV")

  rCV <- rCV %>% mutate(DateTime=Data$DateTime) %>% relocate(DateTime, .before=1)

  rCV <- rCV %>% pivot_longer(!DateTime, names_to = "Detector", values_to = "rCV_Logical")

  Tidy <- TheGainData %>%
    left_join(Gains, by = c("Detector", "DateTime")) %>%
    left_join(TherCVData, by = c("Detector", "DateTime")) %>%
    left_join(rCV, by = c("Detector", "DateTime"))

  TheDetectors <- Tidy %>% pull(Detector) %>% unique()

  Summary <- map(.x=TheDetectors, .f=QCSummaryCheck, data=Tidy) %>% bind_rows()

  Summary <- Summary %>% mutate(Date=Date) %>% relocate(Date, .before=1)
  return(Summary)
}