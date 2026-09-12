#' Internal for QC_ChorusPDF, widens Main Data
#' 
#' @param MainData The passed MainData data.frame containing QC for
#'  the respective detectors
#' 
#' @importFrom tidyr pivot_wider
#' @importFrom dplyr rename mutate case_when select
#' 
#' @noRd
MainDataWiden <- function(MainData){
    MainData <- MainData |> rename(DetectorGain=`Detector Gain`,
     MFI=`MFI-A`, rCV=`%rCV`, LimitResolution= `Limit of Resolution`,
    SystemBackground=`System Background`)
    MainData <- MainData |> mutate(Detector = case_when(
        grepl("^LightLoss", Name) ~ Name,
        grepl("^SSC", Name) ~ Name,
        TRUE ~ sub(" .*", "", Name)))

    MainData$Detector <- gsub("(", "", fixed=TRUE, MainData$Detector)
    MainData$Detector <- gsub(")", "", fixed=TRUE, MainData$Detector)
    MainData$MFI <- gsub(",", "", fixed=TRUE, MainData$MFI)
    MainData$MFI <- as.numeric(MainData$MFI)

    Check <- MainData |> select(-Name, -Filter) |>
  pivot_wider(
    names_from = Detector,
    values_from = DetectorGain:SystemBackground,
    names_glue = "{Detector}_{.value}"
  )
    
    return(Check)
}