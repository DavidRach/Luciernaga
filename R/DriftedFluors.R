
#' Simulation function to approximate fluorophore drift
#' 
#' @param Residual The output of SignatureShifts
#' @param NumberDetectors Number of detectors corresponding to
#'  your cytometer
#' @param TheFluoruophore The fluorophore of interest
#' @param RestingMFI Value by which reference signature gets
#'  multiplied by for this simulation
#' @param legend Default FALSE, TRUE sets on right side plot
#' 
#' @importFrom dplyr filter mutate select left_join filter pull
#' @importFrom tidyselect all_of
#' @importFrom tidyr pivot_longer
#' @importFrom stringr str_detect
#' 
#' @return A signature plot of the reference control vs all. 
#' 
#' @noRd
DriftedFluors <- function(Residual, NumberDetectors=64,
     TheFluorophore, RestingMFI=100000, legend=FALSE){
 
 References <- Luciernaga:::InstrumentReferences(
     NumberDetectors=NumberDetectors)
 Internal <- References |>
      filter(Fluorophore %in% TheFluorophore) |>
      mutate(AdjustedY=AdjustedY*RestingMFI)
 
 FinalCol <- ncol(Residual)
 
 Internal <- Internal |> select(-Instrument, -Fluorophore)
 Residual <- Residual |> tidyr::pivot_longer(all_of(2:FinalCol),
  names_to="Detector", values_to="Adjustment")
 Residual$Detector <- gsub("-A", "", Residual$Detector)
 Merge <- left_join(Residual, Internal, by="Detector")
 Merge <- Merge |> mutate(AdjustedMFI=Adjustment*AdjustedY)
 Merge <- Merge |> filter(!str_detect(Detector, "SC"))
  
 Merge <- Merge |>
      mutate(Signature = AdjustedMFI/max(AdjustedMFI, na.rm = TRUE)) |>
     ungroup()
 
 Merge2 <- Merge |> select(-Adjustment, -AdjustedY, -AdjustedMFI)
  
 Merge2$TheSample <- sub(".*ter_", "", Merge2$TheSample )
 Merge2$TheSample <- sub(".*ore_", "", Merge2$TheSample )
 
 TheseDates <- Merge2 |> pull(TheSample) |> unique()
 
 Plot <- QC_ViewSignature(x=TheseDates,
      columnname="TheSample", data=Merge2,
  TheFormat="longer", detectorcolumn = "Detector",
     valuecolumn = "Signature",
  Normalize=FALSE, legend=legend)

 #plotly::ggplotly(Plot)
 
 return(Plot)
}