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