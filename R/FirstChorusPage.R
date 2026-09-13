#' Internal for QC_ChorusPDF, handles the first page
#'
#' @param x An iterated pdf page read in as text
#'
#' @importFrom stringr str_extract
#' @importFrom lubridate mdy_hm parse_date_time
#' @importFrom purrr map compact
#' @importFrom dplyr bind_cols
#'
#' @noRd
FirstChorusPage <- function(x) {
  lines <- strsplit(x, "\n")[[1]]
  lines <- lines[nzchar(trimws(lines))]

  DateTimeLine <- grep("- Daily Setup and QC", lines)

  if (length(DateTimeLine) == 1) {
    DateTimeLine <- lines[DateTimeLine]
    DateTime_str <- str_extract(DateTimeLine,
      "\\d{2}/\\d{2}/\\d{4}\\s+\\d{1,2}:\\d{2}\\s*(AM|PM)")
    PDFName <- gsub("/", "-", DateTime_str)
    PDFName <- gsub(":", "-", PDFName)
    PDFName <- paste0("Daily Setup and QC-", PDFName, ".pdf")
    DateTime <- lubridate::mdy_hm(DateTime_str)
    DateTime <- data.frame(DateTime = DateTime, PDFName = PDFName)
  }

  ImagingStatusLine <- grep("Imaging Status:", lines)

  if (length(ImagingStatusLine) == 1) {
    ImagingStatusLine <- lines[ImagingStatusLine]
    Status_DF <- TwoPartSplits(ImagingStatusLine)
  }

  ReportTypeLine <- grep("Report Type:", lines)
  if (length(ReportTypeLine) == 1) {
    ReportTypeLine <- lines[ReportTypeLine]
    Reports_DF <- TwoPartSplits(ReportTypeLine)
  }

  UserLine <- grep("User:", lines)
  if (length(UserLine) == 1) {
    UserLine <- lines[UserLine]
    User_DF <- TwoPartSplits(UserLine)
  }

  ReportTypeLine <- grep("Report Type:", lines)
  UserLine <- grep("User:", lines)
  ConfigurationLines <- (ReportTypeLine + 1):(UserLine - 1)

  if (length(ConfigurationLines) > 1) {
    TotalLength <- length(ConfigurationLines)
    ConfigurationLine <- lines[ConfigurationLines[1]]
    ConfigurationDateTime <- str_extract(ConfigurationLine,
      "\\d{2}/\\d{2}/\\d{4}\\s+\\d{1,2}:\\d{2}\\s*(AM|PM)")
    ConfigurationDateTimeValue <- lubridate::mdy_hm(ConfigurationDateTime)
    ConfigurationLine <- gsub(
      ConfigurationDateTime, "TIME    ", ConfigurationLine)
    ConfigurationLine <- gsub(
      "Configuration:", "    Configuration:", ConfigurationLine)
    ConfigurationLine <- gsub("/", "", ConfigurationLine)
    Hmm <- TwoPartSplits(ConfigurationLine)

    if (Hmm[1, 1] == "TIME") {
      Hmm[1, 1] <- as.character(ConfigurationDateTimeValue)
      Hmm[, 1] <- lubridate::parse_date_time(Hmm[, 1], orders = "ymd HMS")
    }

    Sequence <- 1:TotalLength
    Sequence <- Sequence[-1]
    ConfigurationLine <- lines[ConfigurationLines[Sequence]] # Scale to a Third?
    ConfigurationLine <- trimws(ConfigurationLine)

    UpdatedValue <- paste(Hmm$Configuration, ConfigurationLine, sep = " ")
    Hmm[1, "Configuration"] <- UpdatedValue
    SerialNumber_DF <- Hmm
  }

  NameLine <- grep("Name", lines)
  NameLine <- NameLine[-1]

  if (length(NameLine) == 1) {
    FinalLine <- grep("For Research Use Only", lines) - 1
    Data <- lines[NameLine:FinalLine]
    #Data <- gsub("%", "", Data)
    #x <- Data[1]

    Data <- purrr::map(.f = SplitLines, .x = Data)
    Data <- purrr::compact(Data)

    df <- do.call(rbind, Data)
    colnames(df) <- df[1, ]
    QCData <- as.data.frame(df[-1, ],
      stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    stop("Multiple name lines for page one")
    QCData <- NULL
  }

  Metadata <- bind_cols(DateTime, Status_DF, Reports_DF,
    User_DF, SerialNumber_DF)
  Cargo <- list(Metadata, QCData)

  return(Cargo)
}