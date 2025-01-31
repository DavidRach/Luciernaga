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
  IsConventional <- Info_child[xml_name(Info_child) == "IsConventional"]

  if (length(IsConventional) > 0){
    Conventional <- xml_text(IsConventional)
    if (Conventional == "true"){
      message("Conventional experiment, skipping")
      Value <- NULL
      return(Value)
    }
  }

  # Two Older Versions
  if (any(xml_name(Landing) == "ExperimentData")){
    ExperimentData <- Landing[xml_name(Landing) == "ExperimentData"]
    Experiment_child <- xml_children(ExperimentData)
    if (length(Experiment_child) == 0){
      Info <- Landing[xml_name(Landing) == "Info"][[1]]
      Info_child <- xml_children(Info)
      if (any(xml_name(Info_child) == "ExperimentDesc")){
        ExperimentDesc <- Info_child[xml_name(Info_child) == "ExperimentDesc"][[1]]
        Experiment_child <- xml_children(ExperimentDesc)
      } else {message("Missed version for ", x)}
    }
  } else { # More Recent Version
    Info <- Landing[xml_name(Landing) == "Info"][[1]]
    Info_child <- xml_children(Info)
    if (any(xml_name(Info_child) == "ExperimentDesc")){
      ExperimentDesc <- Info_child[xml_name(Info_child) == "ExperimentDesc"][[1]]
      Experiment_child <- xml_children(ExperimentDesc)
    } else {message("Missed version for ", x)}
  }

  RefSetUp <- Experiment_child[xml_name(Experiment_child) == "_RefSetupResult"][[1]]
  RefSetUp_child <- xml_children(RefSetUp)

  if (length(RefSetUp_child) != 0){
  SpillOverColumn <- RefSetUp_child[xml_name(RefSetUp_child) == "SpilloverColumnList"][[1]]
  Spill_child <- xml_children(SpillOverColumn) # Number Children
  Data <- map(.x=Spill_child, .f=Luciernaga:::NormalizedParser) %>% bind_rows()

  if (ColumnNames=="detector"){
    Data <- Luciernaga:::ColumnNaming(x=Data)
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

  #} else {message("Old software version")}

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

    #DateValues <- xml_find_all(Date, ".//d4p1:float", ns = xml_ns(Date))

    RefControl <- Parameters[xml_name(Parameters) == "_RefControlDesc"]
    RefControl_child <- xml_children(RefControl)
    FluorophoreFloat <- RefControl_child[xml_name(RefControl_child) == "Fluorochrome"]
    Fluorophore <- xml_text(FluorophoreFloat)
    Fluorophore <- data.frame(Fluorophore)

    Param_child <- Parameters[xml_name(Parameters) == "_SpilloverVectorArea"]

    FloatingValues <- xml_find_all(Param_child, ".//d7p1:float", ns = xml_ns(Param_child))

    if (length(FloatingValues) == 0){
      FloatingValues <- xml_find_all(Param_child, ".//d6p1:float", ns = xml_ns(Param_child))
    }

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

#' Currently grabs Gain Information Held Parent Level of BD Diva .xml file
#' Poised to expand to grab from specimen_name tube_name for individualized tube level
#' 
#' @param x A .xml file from BD Diva Software
#' 
#' @importFrom xml2 read_xml
#' @importFrom xml2 xml_children
#' @importFrom xml2 xml_name
#' @importFrom xml2 xml_text
#' @importFrom xml2 xml_find_all
#' @importFrom lubridate ymd_hms
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
#' @importFrom tidyr pivot_wider
#' 
#' @return A data.frame row containing the parsed data
#' 
#' @noRd
DivaXMLParse <- function(x){

    Parsed <- read_xml(x)
    Landing <- xml_children(Parsed)
    Experiment <- Landing[xml_name(Landing) == "experiment"][[1]]
    Experiment_child <- xml_children(Experiment)
    Info <- Experiment_child[xml_name(Experiment_child ) == "date"][[1]]
    Date <- xml_text(Info)
    Date <- ymd_hms(Date)
    Info <- Experiment_child[xml_name(Experiment_child ) == "owner_name"][[1]]
    Owner <- xml_text(Info)
    
    # Failed to filter by xml_node, selecting by position, REVISE
    Keywords <- Experiment_child[xml_name(Experiment_child ) == "keywords"][[1]]
    Keywords_child <- xml_children(Keywords)
    CytometerConfig <- Keywords_child[[3]]
    CytometerConfig_child <- xml_children(CytometerConfig)
    Nozzle <- CytometerConfig_child[xml_name(CytometerConfig_child) == "value"]
    Nozzle <- xml_text(Nozzle)
    
    Instrument <- Experiment_child[xml_name(Experiment_child ) == "instrument_settings"][[1]]
    TheFluors <- xml_find_all(Instrument, ".//parameter[not(contains(@name, '-H') or contains(@name, '-W'))]")
    
    FluorophoreGains <- map(.x=TheFluors, .f=DivaParseInternal) |> 
       bind_rows()
    FluorophoreGains <- FluorophoreGains |> 
        pivot_wider(names_from = "Fluorophore", values_from = "Gain")
    
    Dataset <- cbind(Date, Owner, Nozzle, FluorophoreGains)
    
    # REVISE: Discovered .fcs files were indeed present in the folder, so no need to further process. 
    # Note, .xml file proceed to specimen_name, identify tubename, and process for all the above individually
    # as the current structure keeps QC for entire month in a single xml.
    
    return(Dataset)
    }
    
#' Internal for DivaParse, iterates on parameters for Fluorophore
#' 
#' @param x The iterated xml_node for each fluorophore being parsed
#' @importFrom xml2 xml_children
#' @importFrom xml2 xml_name
#' @importFrom xml2 xml_text
#' 
#' @return The iterated data.frame row of Fluorophore and Gain
#' @noRd
DivaParseInternal <- function(x){
    TheData <- xml_children(x)
    TheFluor <- TheData[xml_name(TheData) == "fl"]
    Fluorophore <- xml_text(TheFluor)
    Gain <- TheData[xml_name(TheData) == "voltage"]
    Gain <- xml_text(Gain)
    Data <- data.frame(cbind(Fluorophore, Gain))
    Data$Gain <- as.numeric(Data$Gain)
    return(Data)
    }

