#' Internal for DivaXMLParse, operates on iterated tubes to parse the
#' needed information
#' 
#' @param x An iterated xml_node corresponding to the tube
#' 
#' @importFrom xml2 xml_children xml_name xml_text xml_find_all
#'  xml_attr
#' @importFrom lubridate ymd_hms
#' @importFrom purrr map2
#' @importFrom dplyr bind_rows select
#' @importFrom tidyr pivot_wider
#' 
#' @return A data.frame row containing the parsed data
#' 
#' @noRd
TubeIterate <- function(x){

  Landing <- xml_children(x)
  Date <- Landing[xml_name(Landing) == "date"]
  Date <- xml_text(Date)
  DateTime <- ymd_hms(Date)

  FileName <- Landing[xml_name(Landing) == "data_filename"]
  FileName <- xml_text(FileName)

  Instrument <- Landing[xml_name(Landing) == "data_instrument_name"]
  Instrument <- xml_text(Instrument)

  InstrumentSN <- Landing[xml_name(Landing) == "data_instrument_serial_number"]
  InstrumentSN <- xml_text(InstrumentSN)

  FSC_Area_Scaling <- Landing[xml_name(Landing) == "fsc_area_scaling"]
  FSCAreaScaling <- xml_text(FSC_Area_Scaling)

  User <- Landing[xml_name(Landing) == "record_user"]
  User <- xml_text(User)

  Keywords <- Landing[xml_name(Landing) == "keywords"][[1]]
  Keywords_child <- xml_children(Keywords)
  CytometerConfig <- xml_find_all(Keywords, ".//keyword[@name='CYTOMETER CONFIG NAME']")
  CytometerConfig_child <- xml_children(CytometerConfig)
  Nozzle <- CytometerConfig_child[xml_name(CytometerConfig_child) == "value"]
  Nozzle <- xml_text(Nozzle)

  if(!any(xml_name(Landing) == "instrument_settings")){
    message("Instrument Settings Node Absent for ", FileName, " returning NULL")
    Dataset <- NULL
  } else{

  InstrumentSettings <- Landing[xml_name(Landing) == "instrument_settings"][[1]]
  TheFluors <- xml_find_all(InstrumentSettings, ".//parameter[@name]")
  Fluorophores <- xml_attr(TheFluors, "name")

  FluorophoreGains <- map2(.x=TheFluors, .y=Fluorophores, .f=DivaParseInternal) |> 
     bind_rows()

  FluorophoreGains <- FluorophoreGains |> 
      pivot_wider(names_from = "Fluorophore", values_from = "Gain")

  if(!any(xml_name(Landing) == "lasers")){
      message("Laser Settings Node Absent for ", FileName, " returning NULL")
      Dataset <- NULL
  } else{
  
  TheLasers <- Landing[xml_name(Landing) == "lasers"][[1]]
  TheLaser_child <- xml_children(TheLasers)
  Laser <- xml_attr(TheLaser_child, "name")

  Lasers <- map2(.x=TheLaser_child, .y=Laser, .f=DivaLaserParseInternal)|> bind_rows()
  Lasers$Delay <- as.numeric(Lasers$Delay)
  Lasers$AreaScaling <- as.numeric(Lasers$AreaScaling)

  LaserDelays <- Lasers |> select(Laser, Delay) |>
    pivot_wider(names_from = "Laser", values_from = "Delay")
  colnames(LaserDelays) <- paste0(colnames(LaserDelays), "-Laser Delay")

  LaserASF <- Lasers |> select(Laser, AreaScaling) |>
    pivot_wider(names_from = "Laser", values_from = "AreaScaling")
  colnames(LaserASF ) <- paste0(colnames(LaserASF ), "-Area Scaling Factor")

  Dataset <- cbind(DateTime, User, FileName, Instrument,
     InstrumentSN, Nozzle, FluorophoreGains, LaserDelays,
      LaserASF, FSCAreaScaling)
  }
    
  }
  return(Dataset)
  }