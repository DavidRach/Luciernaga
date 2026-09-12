
#' Takes the normalize false data output from NetSignatures and
#' returns the shift matrix
#' 
#' @param data The data output from NetSignatures
#' @param sampleName Column name of the sample, default is TheSample
#' 
#' @importFrom dplyr filter select across cur_column bind_cols mutate
#' @importFrom tidyselect where all_of everything 
#' 
#' @return A drift data.frame for use in simulating fluorophores
#' 
#' @noRd
SignatureShifts <- function(data, sampleName="TheSample"){
  TheNames <- colnames(data)
  TheNames <- TheNames[-1]
  AverageReference <- data |> filter(.data[[sampleName]] %in% "Average") |>
    select(where(is.numeric)) |> unlist()

  OtherReferences <- data |> filter(!.data[[sampleName]] %in% "Average")
  OtherMetadata <- OtherReferences |> select(all_of(sampleName))
  OtherNumeric <- OtherReferences |> select(!all_of(sampleName))

  SubtractedData <- OtherNumeric %>% mutate(across(
    .cols = everything(), .fns  = ~ .x - AverageReference[cur_column()]))

  DriftData <- SubtractedData %>% mutate(across(
      everything(),~ .x / AverageReference[cur_column()]))
  
  OnesData <- data.frame(matrix(1, 
    nrow = nrow(DriftData), 
    ncol = ncol(DriftData),
    dimnames = list(rownames(DriftData), colnames(DriftData))))
  
  Result <- OnesData + DriftData

  colnames(Result) <- TheNames
  
  DriftData <- bind_cols(OtherMetadata, Result)
  
  return(DriftData)
}