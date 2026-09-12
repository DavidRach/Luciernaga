 
#' Internal for Holistic to Archived, screens vs cutoff value, filling Flag column
#' 
#' @param x The column being iterated
#' @param data The data being selected on
#' @param cutoffs The processed Gain and RCV cutoff criteria
#' 
#' @importFrom stringr str_detect
#' @importFrom dplyr select filter case_when mutate pull
#' @importFrom tidyselect all_of
#' 
#' @return Individual flag columns
#' 
#' @noRd
DerriveTheFlag <- function(x, data, cutoffs){
  
    if (str_detect(x, "Gain")) {Type <- "Gain"
    } else if (str_detect(x, "rCV")) {Type <- "rCV"
    } else {Type <- "MFI"
    } 
  
    TheColumn <- gsub("-% rCV", "", gsub("_Gain", "", x))
    TheColumn <- gsub("-A", "", TheColumn)
    internaldata <- data |> dplyr::select(all_of(x))
    internalcutoffs <- cutoffs |> dplyr::filter(Detector %in% TheColumn)
  
    if (nrow(internalcutoffs) == 1){
      if (Type == "Gain"){
        Comparison <- internalcutoffs |> dplyr::pull(GainBaseline)
      } else if (Type == "rCV") {
        Comparison <- internalcutoffs |> dplyr::pull(RCVCutoff)
      } else {message("Skipped")
      }
    }
  
    Status <- internaldata |> mutate(Flag = .data[[x]] >= Comparison)
    Status <- Status |> dplyr::select(-all_of(x))
    NewName <- paste0("Flag-", x)
    colnames(Status)[1] <- NewName
    return(Status)
  }