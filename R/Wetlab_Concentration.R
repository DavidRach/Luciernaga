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
#' @importFrom dplyr mutate
#' @importFrom dplyr relocate
#'
#' @return A data.frame of useful information
#' @export
#'
#' @examples
#'
#' library(flowCore)
#' library(flowWorkspace)
#' library(CytoML)
#' library(dplyr)
#' library(purrr)
#'
#' File_Location <- system.file("extdata", package = "Luciernaga")
#' WSP_File <- list.files(File_Location, pattern=".wsp", full.names = TRUE)
#' ws <- open_flowjo_xml(WSP_File[1])
#' gs <- flowjo_to_gatingset(ws, name=1, path = File_Location)
#' nameKeyword <- c("GROUPNAME", "TUBENAME")
#'
#' TheData <- map(.x=gs, Wetlab_Concentration, subset = "CD45+",
#'   nameKeyword=nameKeyword, DilutionMultiplier=100, TotalVolume=1) %>%
#'    bind_rows()
#'
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

  #InstrumentIdentifier <-keyword(x)[[1]][["$CYT"]]

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
                  TimeSeconds, Instrument, Date)
  TheRow <- data.frame(TheRow)

  TheRow <- TheRow %>% mutate(TotalVolume=TotalVolume) %>%
    relocate(TotalVolume, .after=ConcentrationScientific)
  return(TheRow)
}
