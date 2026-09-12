
#' Dashboard Internal, holistic passing summary
#'
#' @param x The Instrument Name
#' @param y The Instrument Data
#'
#' @return Global Passing Status
#' @noRd
ColorCodeStatus <- function(x, y){

  data <- y

  CytekMandate <- c("FSC", "SSC", "SSC-B")
  CytekData <- data %>%
    filter(Detector %in% CytekMandate)
  RCVdata <- data %>%
    filter(str_detect(Detector, "3")) %>%
    filter(!str_detect(Detector, "1"))
  RCVdata <- rbind(CytekData, RCVdata)

  if(nrow(data)== 0){ColorCode <- "Unknown"}

  if (nrow(data) > 0){

  if (any(data$Gain == "Red") || any(RCVdata$rCV == "Red")) {
    ColorCode <- "Red" # Overall QC Fail
  } else if (any(data$rCV == "Red")) {
    ColorCode <- "Orange" # Non-primary RCV Fail
  } else if (any(data$Gain == "Yellow") || any(data$rCV == "Yellow")) {
    ColorCode <- "Yellow"
  } else {ColorCode <- "Green"}

  }

  QCResults <- data.frame(Instrument = x, QCStatus=ColorCode)
  return(QCResults)
}