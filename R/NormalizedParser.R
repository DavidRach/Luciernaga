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