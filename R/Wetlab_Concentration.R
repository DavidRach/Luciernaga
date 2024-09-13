#' Takes gated .fcs files and returns concentration and other info useful for Wetlab users.
#'
#' @param x A mapped gating set object
#' @param subset Desired population node to derrive counts from
#' @param nameKeyword Keyword containing samples name (ex. "GROUPNAME" or c("GROUPNAME", "TUBENAME"))
#' @param DilutionMultiplier The dilution multiplier for the sample
#' @param TotalVolume Volume the specimen was resuspended in.
#'
#' @importFrom flowWorkspace gs_pop_get_data
#' @importFrom flowCore keyword
#' @importFrom BiocGenerics nrow
#' @importFrom lubridate hms
#'
#' @return A data.frame of useful information
#' @export
#'
#' @examples NULL
Wetlab_Concentration <- function(x, subset, nameKeyword, DilutionMultiplier,
                                     TotalVolume){
  CS <- gs_pop_get_data(x, subset)

  if (length(nameKeyword)==2){
    first <- nameKeyword[[1]]
    second <- nameKeyword[[2]]
    first <- keyword(x, first)
    second <- keyword(x, second)
    name <- paste(first, second, sep="_")
  } else {name <- keyword(x, nameKeyword)}

  TotalEvents <- keyword(x)$`$TOT`
  Cells <- BiocGenerics::nrow(CS)[[1]]
  Volume <- as.numeric(keyword(x)$`$VOL`)
  Concentration <- (Cells*1000)/Volume
  Concentration <- Concentration*DilutionMultiplier
  ConcentrationScientific <- format(Concentration, scientific = TRUE,
                                    digits = 2)
  TheTotal <- TotalVolume*Concentration
  TotalScientific <- format(TheTotal, scientific=TRUE, digits = 2)

  StartTime <- lubridate::hms(keyword(x)$`$BTIM`)
  EndTime <- lubridate::hms(keyword(x)$`$ETIM`)
  TimeDifference <- EndTime - StartTime
  TimeSeconds <- as.numeric(TimeDifference, units = "secs")
  TimeSeconds <- round(TimeSeconds, 1)

  FlowRate <- keyword(x)$`$FLOWRATE`

  Instrument <- keyword(x)$`$CYTSN`

  fr <- CS[[1, returnType = "flowFrame"]]
  new_kw <- fr@description
  TypeParams <- new_kw[grepl("^\\$P\\d+TYPE\\d*", names(new_kw))]
  TypeParams <- TypeParams[grepl("^\\Raw_Fluorescence\\d*", TypeParams)]
  Detectors <- length(TypeParams)

  if (Detectors == 64) {Instrument <- "5L"
  } else if (Detectors == 54) {Instrument <- "4L_UV"
  } else if (Detectors == 48) {Instrument <- "4L_YG"
  } else if (Detectors == 38) {Instrument <- "3L"
  } else if (Detectors == 30) {Instrument <- "2L_V"
  } else if (Detectors == 22) {Instrument <- "2L_R"
  } else if (Detectors == 14) {Instrument <- "1L"}

  Date <- keyword(x)$`$DATE`

  TheRow <- cbind(name, Cells, Volume, ConcentrationScientific, TotalScientific,
                  TimeSeconds, FlowRate, Instrument, Date)
  TheRow <- data.frame(TheRow)
  return(TheRow)
}
