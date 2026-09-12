#' Internal For Expt Parse
#'
#' @param x The tidyed data
#'
#' @return The tidyed data plus renamed columns for detectors
#' @noRd
ColumnNaming <- function(x){
    TotalDetectors <- ncol(x)-2

    The5L <- c("DateTime", "Fluorophore", "UV1", "UV2", "UV3", "UV4",
     "UV5", "UV6", "UV7", "UV8", "UV9", "UV10", "UV11", "UV12", "UV13",
     "UV14", "UV15", "UV16",
     "V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11",
     "V12", "V13", "V14", "V15", "V16",
     "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9", "B10", "B11",
     "B12", "B13", "B14",
     "YG1", "YG2", "YG3", "YG4", "YG5", "YG6", "YG7", "YG8", "YG9", "YG10",
     "R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8")
    The4LUV <- c("DateTime", "Fluorophore", "UV1", "UV2", "UV3", "UV4", "UV5", "UV6", "UV7", "UV8",
      "UV9", "UV10", "UV11", "UV12", "UV13", "UV14", "UV15", "UV16",
      "V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8",
      "V9", "V10", "V11","V12", "V13", "V14", "V15", "V16",
      "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8",
      "B9", "B10", "B11", "B12", "B13", "B14",
      "R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8")
    The4LYG <- c("DateTime", "Fluorophore", "V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8",
      "V9", "V10", "V11","V12", "V13", "V14", "V15", "V16",
      "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8",
      "B9", "B10", "B11", "B12", "B13", "B14",
      "YG1", "YG2", "YG3", "YG4", "YG5", "YG6", "YG7", "YG8", "YG9", "YG10",
      "R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8")
    The3L <- c("DateTime", "Fluorophore", "V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8",
      "V9", "V10", "V11","V12", "V13", "V14", "V15", "V16",
      "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8",
      "B9", "B10", "B11", "B12", "B13", "B14",
      "R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8")
    The2LVB <- c("DateTime", "Fluorophore", "V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8",
      "V9", "V10", "V11","V12", "V13", "V14", "V15", "V16",
      "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8",
      "B9", "B10", "B11", "B12", "B13", "B14")
    The2LBR <- c("DateTime", "Fluorophore", "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8",
      "B9", "B10", "B11", "B12", "B13", "B14",
      "R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8")
    The1L <- c("DateTime", "Fluorophore", "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8",
      "B9", "B10", "B11", "B12", "B13", "B14")

    if (TotalDetectors == 64){colnames(x) <- The5L
    } else if (TotalDetectors == 54){colnames(x) <- The4LUV
    } else if (TotalDetectors == 48){colnames(x) <- The4LYG
    } else if (TotalDetectors == 38){colnames(x) <- The3L
    } else if (TotalDetectors == 30){colnames(x) <- The2LVB
    } else if (TotalDetectors == 22){colnames(x) <- The2LBR
    } else if (TotalDetectors == 14){colnames(x) <- The1L
    } else {message("Number of Columns didn't match known Instrument")
    }

    return(x)
  }