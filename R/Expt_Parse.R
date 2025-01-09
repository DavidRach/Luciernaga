#' Parse single color control signatures from .Expt file
#'
#' @param x File.path to the .Expt file
#' @param ColumnNames Default is "detector", else X1 numbers
#' @param returnType Either "data" or "plot"
#'
#' @importFrom xml2 read_xml
#' @importFrom xml2 xml_children
#' @importFrom xml2 xml_name
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
#'
#' @return A tidy-data frame of all fluorophores and
#' normalized MFI values by detector
#' @noRd
SpectroFloSignatureParser <- function(x, ColumnNames="detector", returnType="data"){
  Parsed <- read_xml(x)
  Landing <- xml_children(Parsed)
  Info <- Landing[xml_name(Landing) == "Info"][[1]]
  Info_child <- xml_children(Info)
  if (any(xml_name(Info_child) == "ExperimentDesc")){
  ExperimentDesc <- Info_child[xml_name(Info_child) == "ExperimentDesc"][[1]]
  Experiment_child <- xml_children(ExperimentDesc)
  RefSetUp <- Experiment_child[xml_name(Experiment_child) == "_RefSetupResult"][[1]]
  RefSetUp_child <- xml_children(RefSetUp)

  if (length(RefSetUp_child) != 0){
  SpillOverColumn <- RefSetUp_child[xml_name(RefSetUp_child) == "SpilloverColumnList"][[1]]
  Spill_child <- xml_children(SpillOverColumn) # Number Children
  Data <- map(.x=Spill_child, .f=NormalizedParser) %>% bind_rows()

  if (ColumnNames=="detector"){
    Data <- ColumnNaming(x=Data)
  }

  if (returnType == "data"){
    return(Data)
  } else if (returnType == "plot"){
    plot <- PlotlySignatures(data=Data)
    return(plot)
  }

  } else {message(".Expt lacked Unmixing parameters (was Raw), returning a NULL value instead of intended data.frame row")
    Value <- NULL
    return(Value)
  }

  } else {message("Old software version")}

}

#' Internal for expt_parse
#'
#' @param data A data.frame of Fluorophore and n-detector columns
#'
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr pull
#' @importFrom tidyselect where
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#' @importFrom plotly ggplotly
#'
#' @return An interactive plotly object of the line signatures
#' @noRd
PlotlySignatures <- function(data, TheFactor="Fluorophore"){
  Tidyed <- data %>% pivot_longer(cols = where(is.numeric),
                                  names_to = "Detector",
                                  values_to = "Value")

  DetectorOrder <- Tidyed %>% pull(Detector) %>% unique()
  Tidyed$Detector <- factor(Tidyed$Detector, levels=DetectorOrder)

  plot <- ggplot(Tidyed, aes(x = Detector, y = Value, color = .data[[TheFactor]], group = .data[[TheFactor]])) +
    geom_line() + labs(x = "Detector", y = "Normalized MFI", color = "Fluorophore") +
    theme_bw() + theme(axis.text.x = element_text(size=5, angle = 70, hjust = 1))

  plot <- plotly::ggplotly(plot)
  return(plot)
}

#' Internal For Expt Parse
#'
#' @param x Passed xml_node child to extract signature values from
#'
#' @importFrom xml2 xml_children
#' @importFrom xml2 xml_name
#' @importFrom xml2 xml_text
#' @importFrom xml2 xml_find_all
#' @importFrom xml2 xml_ns
#' @importFrom lubridate ymd
#' @importFrom lubridate hms
#'
#' @return An assembled row corresponding to the nodes fluorophore
#' @noRd
  NormalizedParser <- function(x){
    Parameters <- xml_children(x)
    Date <- Parameters[xml_name(Parameters) == "_DateTimeCreated"]
    DateTime <- xml_text(Date)
    TheDate <- sub("^(\\d{4}-\\d{2}-\\d{2})T.*", "\\1", DateTime)
    TheDate <- ymd(TheDate)
    TheTime <- sub("^.*T(\\d{2}:\\d{2}:\\d{2}).*", "\\1", DateTime)
    TheTime <- hms(TheTime)
    DateTime <- TheDate + TheTime
    DateTime <- data.frame(DateTime)

    DateValues <- xml_find_all(Date, ".//d4p1:float", ns = xml_ns(Date))

    RefControl <- Parameters[xml_name(Parameters) == "_RefControlDesc"]
    RefControl_child <- xml_children(RefControl)
    FluorophoreFloat <- RefControl_child[xml_name(RefControl_child) == "Fluorochrome"]
    Fluorophore <- xml_text(FluorophoreFloat)
    Fluorophore <- data.frame(Fluorophore)

    Param_child <- Parameters[xml_name(Parameters) == "_SpilloverVectorArea"]
    FloatingValues <- xml_find_all(Param_child, ".//d7p1:float", ns = xml_ns(Param_child))
    ValueVector <- as.numeric(xml2::xml_text(FloatingValues))
    Data <- data.frame(t(ValueVector))
    Data <- cbind(DateTime, Fluorophore, Data)
    return(Data)
  }

#' Internal For Expt Parse
#'
#' @param x The tidyed data
#'
#' @return The tidyed data plus renamed columns for detectors
#' @noRd
ColumnNaming <- function(x){
    TotalDetectors <- ncol(x)-2

    The5L <- c("DateTime", "Fluorophore", "UV1", "UV2", "UV3", "UV4", "UV5", "UV6", "UV7", "UV8",
               "UV9", "UV10", "UV11", "UV12", "UV13", "UV14", "UV15", "UV16",
               "V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8",
               "V9", "V10", "V11","V12", "V13", "V14", "V15", "V16",
               "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8",
               "B9", "B10", "B11", "B12", "B13", "B14",
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

