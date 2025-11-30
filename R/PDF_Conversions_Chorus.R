
#' Internal for PDF Processing, handles Chorus 7 column QC params
#' 
#' @param x A read-in line from the pdf file
#' 
#' @noRd
SplitLines <- function(x){
    parts <- strsplit(trimws(x), "\\s{2,}")[[1]]
    if (!(length(parts) %in% c(3, 7))){parts <- NULL}
    return(parts)
}

#' Internal for QC_ChorusPDF, handles the first page
#' 
#' @param x An iterated pdf page read in as text
#' 
#' @importFrom stringr str_extract
#' @importFrom lubridate mdy_hm
#' @importFrom lubridate parse_date_time
#' @importFrom purrr map
#' @importFrom purrr compact
#' @importFrom dplyr bind_cols
#' 
#' @noRd
FirstChorusPage <- function(x){
lines <- strsplit(x, "\n")[[1]]
lines <- lines[nzchar(trimws(lines))]

DateTimeLine <- grep("- Daily Setup and QC", lines)

if (length(DateTimeLine) == 1){
   DateTimeLine <-lines[DateTimeLine]
   DateTime_str <- str_extract(DateTimeLine,
    "\\d{2}/\\d{2}/\\d{4}\\s+\\d{1,2}:\\d{2}\\s*(AM|PM)")
   PDFName <- gsub("/", "-", DateTime_str)
   PDFName <- gsub(":", "-", PDFName)
   PDFName <- paste0("Daily Setup and QC-", PDFName, ".pdf")
   DateTime <- lubridate::mdy_hm(DateTime_str)
   DateTime <- data.frame(DateTime=DateTime, PDFName=PDFName)
}

ImagingStatusLine <- grep("Imaging Status:", lines)

if (length(ImagingStatusLine) == 1){
    ImagingStatusLine <-lines[ImagingStatusLine]
    Status_DF <- TwoPartSplits(ImagingStatusLine)
}

ReportTypeLine <- grep("Report Type:", lines)
if (length(ReportTypeLine) == 1){
    ReportTypeLine <-lines[ReportTypeLine]
    Reports_DF <- TwoPartSplits(ReportTypeLine)
}

UserLine <- grep("User:", lines)
if (length(UserLine) == 1){
    UserLine <-lines[UserLine]
    User_DF <- TwoPartSplits(UserLine)
}

ReportTypeLine <- grep("Report Type:", lines)
UserLine <- grep("User:", lines)
ConfigurationLines <- (ReportTypeLine+1):(UserLine-1)

if (length(ConfigurationLines) > 1){
    TotalLength <- length(ConfigurationLines)
    ConfigurationLine <-lines[ConfigurationLines[1]]
    ConfigurationDateTime <- str_extract(ConfigurationLine,
        "\\d{2}/\\d{2}/\\d{4}\\s+\\d{1,2}:\\d{2}\\s*(AM|PM)")
    ConfigurationDateTimeValue <- lubridate::mdy_hm(ConfigurationDateTime)
    ConfigurationLine <- gsub(ConfigurationDateTime, "TIME    ", ConfigurationLine)
    ConfigurationLine <- gsub("Configuration:", "    Configuration:", ConfigurationLine)
    ConfigurationLine <- gsub("/", "", ConfigurationLine)
    Hmm <- TwoPartSplits(ConfigurationLine)

    if (Hmm[1,1] == "TIME"){
        Hmm[1,1] <- as.character(ConfigurationDateTimeValue)
        Hmm[,1] <- lubridate::parse_date_time(Hmm[,1], orders = "ymd HMS")
    }

    Sequence <- 1:TotalLength
    Sequence <- Sequence[-1]
    ConfigurationLine <- lines[ConfigurationLines[Sequence]] # Scale to a Third?
    ConfigurationLine <- trimws(ConfigurationLine)

    UpdatedValue <- paste(Hmm$Configuration, ConfigurationLine, sep=" ")
    Hmm[1, "Configuration"] <- UpdatedValue
    SerialNumber_DF <- Hmm
}

NameLine <- grep("Name", lines)
NameLine <- NameLine[-1]

if (length(NameLine) == 1){
    FinalLine <- grep("For Research Use Only", lines) - 1
    Data <- lines[NameLine:FinalLine]
    #Data <- gsub("%", "", Data)
    #x <- Data[1]

    Data <- purrr::map(.f=SplitLines, .x=Data)
    Data <- purrr::compact(Data)

    df <- do.call(rbind, Data)
    colnames(df) <- df[1,]
    QCData <- as.data.frame(df[-1,], stringsAsFactors = FALSE, check.names=FALSE)
} else {stop("Multiple name lines for page one")
        QCData <- NULL
}

Metadata <- bind_cols(DateTime, Status_DF, Reports_DF, User_DF, SerialNumber_DF)
Cargo <- list(Metadata, QCData)

return(Cargo)
}

#' Internal for QC_ChorusPDF, handles scrambled first line artefact
#' for the second pages and beyond, returning an unscrambled version
#' 
#' @param lines The cleaned read in lines from the QC pdf
#' 
#' @importFrom stringr str_extract
#' 
#' @noRd
FixFirstLine <- function(lines){

MetadataLine <- grep("Daily Setup and QC", lines)
MetadataLine <- lines[MetadataLine]

Date_Search <- stringr::str_extract(MetadataLine, "\\d{2}-\\d{2}-\\d{4}\\s+\\d{2}-\\d{2}\\s+[AP]M")
Parts <- strsplit(Date_Search, " ")[[1]]
Combined <- gsub("-", "/", Parts[1])
lines <- gsub(Combined, "", lines)
Combined <- gsub("-", ":", Parts[2])
lines <- gsub(Combined, "", lines)
lines <- gsub(Parts[3], "", lines)
lines <- gsub("Setup and QC", "", lines)
lines <- gsub("- Daily", "", lines)
x <- trimws(lines[2])

if (grepl("\\d", x)){#Is a detector page
DetectorPortion <- sub("\\s*\\(.*", "", x)
DetectorPortion <- gsub("\\s+", "", DetectorPortion)
Prefix <- gsub("\\d+", "", DetectorPortion)
Number <- as.integer(gsub("\\D+", "", DetectorPortion))
Number <- Number + 1
NextLine <- paste0(Prefix, Number)
NextLineIndex <- grep(NextLine, lines)

InitialDetector <- paste(lines[2:(NextLineIndex-1)], collapse = " ") 
TheFirstDetector <- sub("\\s*\\(\\s*", " (", InitialDetector)
Cargo <- c(TheFirstDetector, NextLineIndex)
} else {
    StartIndex <- grep(x, lines)
    InitialDetector <- paste(lines[StartIndex:(StartIndex+1)], collapse = " ") 
    NextLineIndex <- StartIndex+2
    Cargo <- c(InitialDetector, NextLineIndex)
}

return(Cargo)
}


#' Processes BD Chorus QC PDF files into .csv files
#' 
#' @param x A file.path to the desired QC pdf
#' @param returnPreference When both modes present, whether
#'  to return QC values from Imaging or High-Speed settings
#' @param returnType Default data, alternative is csv
#' @param outpath When return type is csv, file.path to store the csv
#' 
#' @importFrom pdftools pdf_text
#' @importFrom purrr map
#' @importFrom utils write.csv
#' 
#' @export
#' 
#' @examples A <- 2+2
#' 
QC_ChorusPDF <- function(x, returnPreference="Imaging", returnType="data", outpath=NULL){
    text <- pdftools::pdf_text(x)
    NumberPages <- length(text)
    FirstPageCargo <- FirstChorusPage(x=text[1])
    Metadata <- FirstPageCargo[[1]]
    FirstPage <- FirstPageCargo[[2]]

    Works <- purrr::map(.x=text[2:NumberPages], .f=AdditionalPageHandler)
    Dataset <- Consolidator(x=Works, Metadata=Metadata,
     FirstPage=FirstPage, returnPreference=returnPreference)
    
    if (returnType != "data"){
        if (is.null(outpath)){outpath <- getwd()}
        PDFName <- Dataset$PDFName |> unique()
        if (returnPreference != "Imaging"){
            AppendValue <- paste0(returnPreference, ".csv")
        } else {AppendValue <- ".csv"}
        PDFName <- gsub(".pdf", AppendValue, PDFName)
        StorageLocation <- file.path(outpath, PDFName)
        write.csv(Dataset, StorageLocation, row.names=FALSE)
    }

    return(Dataset)
}

#' Takes AdditionalPageHandler outputs and consolidates them
#' 
#' @param x The list of data.frames from AdditionalPageHandler
#' @param Metadata The Metadata data.frame from QC_Chorus
#' @param FirstPage The First page data.frame from QC_Chorus
#' @param returnPreference Whether to return Imaging or HighSpeed data
#' 
#' @importFrom dplyr bind_rows
#' @importFrom dplyr rename
#' @importFrom tidyr pivot_wider
#' 
#' @noRd
Consolidator <- function(x, Metadata, FirstPage, returnPreference){
    #Metadata
    MainColumns <- ncol(FirstPage)
    MainColumnsNames <-colnames(FirstPage)
    TheList <- x

    FlatList <- unlist(lapply(TheList, function(x) {
        if (is.data.frame(x)){list(x)
            } else if (is.list(x)) {x
            } else {NULL}}), recursive = FALSE)

    Ncols <- sapply(FlatList, ncol)
    RunID <- cumsum(c(TRUE, diff(Ncols) != 0))
    Grouped <- split(FlatList, RunID)
    CombinedList <- lapply(Grouped, function(grp) do.call(rbind, grp))

    FirstData <- CombinedList[[1]]
    colnames(FirstData) <- MainColumnsNames
    MainData <- bind_rows(FirstPage, FirstData)

    SecondData <- CombinedList[[2]]
    colnames(SecondData) <- as.character(unlist(SecondData[1, ]))
    SecondData <- SecondData[-1, , drop = FALSE]

    if (length(CombinedList) == 4){
        ThirdData <- CombinedList[[3]]
        colnames(ThirdData) <- as.character(unlist(ThirdData[1, ]))
        ThirdData <- ThirdData[-1, , drop = FALSE]

        FourthData <- CombinedList[[2]]
        colnames(FourthData) <- as.character(unlist(FourthData[1, ]))
        FourthData <- FourthData[-1, , drop = FALSE]
    }
    
    if (length(CombinedList) == 4 && returnPreference != "Imaging"){
    MainData <- ThirdData

    FourthDataWide <- FourthData |>
        rename(LaserDelay = `Laser Delay`, 
        LaserPowerWithinSpec = `Laser Power within Spec`) |>
        pivot_wider(names_from = Laser,
             values_from = c(LaserDelay, LaserPowerWithinSpec),
             names_glue = "{Laser}_{.value}")

    MetaFull <- cbind(Metadata, FourthDataWide)
    MainData <- MainDataWiden(MainData)
    MetaExpanded <- MetaFull
    #MetaExpanded <- MetaFull[rep(1, nrow(MainData)), , drop = FALSE]
    FinalData <- cbind(MetaExpanded, MainData)
    } else if (returnPreference== "Imaging"){

    SecondDataWide <- SecondData |>
        rename(LaserDelay = `Laser Delay`, 
        LaserPowerWithinSpec = `Laser Power within Spec`) |>
        pivot_wider(names_from = Laser,
             values_from = c(LaserDelay, LaserPowerWithinSpec),
             names_glue = "{Laser}_{.value}")
    
    MainData <- MainDataWiden(MainData)
    MetaExpanded <- MetaFull

    MetaFull <- cbind(Metadata, SecondDataWide)
    #MetaExpanded <- MetaFull[rep(1, nrow(MainData)), , drop = FALSE]
    FinalData <- cbind(MetaExpanded, MainData)
    }



   return(FinalData)
}

#' Internal for QC_ChorusPDF, widens Main Data
#' 
#' @param MainData The passed MainData data.frame containing QC for
#'  the respective detectors
#' 
#' @importFrom tidyr pivot_wider
#' @importFrom dplyr rename mutate case_when select
#' 
#' @noRd
MainDataWiden <- function(MainData){
    MainData <- MainData |> rename(DetectorGain=`Detector Gain`,
     MFI=`MFI-A`, rCV=`%rCV`, LimitResolution= `Limit of Resolution`,
    SystemBackground=`System Background`)
    MainData <- MainData |> mutate(Detector = case_when(
        grepl("^LightLoss", Name) ~ Name,
        grepl("^SSC", Name) ~ Name,
        TRUE ~ sub(" .*", "", Name)))

    MainData$Detector <- gsub("(", "", fixed=TRUE, MainData$Detector)
    MainData$Detector <- gsub(")", "", fixed=TRUE, MainData$Detector)
    MainData$MFI <- gsub(",", "", fixed=TRUE, MainData$MFI)
    MainData$MFI <- as.numeric(MainData$MFI)

    Check <- MainData |> select(-Name, -Filter) |>
  pivot_wider(
    names_from = Detector,
    values_from = DetectorGain:SystemBackground,
    names_glue = "{Detector}_{.value}"
  )
    
    return(Check)
}

#' Internal for QC_ChorusPDF, handles second plus pages
#' 
#' @param x The iterated in pdf page as text lines
#' 
#' @importFrom purrr map
#' @importFrom purrr compact
#' 
#' @noRd
AdditionalPageHandler <- function(x){
lines <- strsplit(x, "\n")[[1]]
lines <- lines[nzchar(trimws(lines))]
    
if (length(lines) == 3){PageData <- NULL
    return(PageData)
}

DetectorLines <- grep("            DE TE CTO R SE T TING S", lines)

if (length(DetectorLines) == 0){

LaserLines <- grep("TINGS", lines)

if (length(LaserLines) == 0){
    lines

    Cargo <- FixFirstLine(lines)
    NextLineIndex <- Cargo[[2]]
    TheFirstDetector <- Cargo[[1]]

    Final <- grep("For Research Use Only", lines)
    Segment <- lines[NextLineIndex:(Final-1)]
    Segment <- c(TheFirstDetector, Segment)

    Data <- purrr::map(.f=SplitLines, .x=Segment)
    Data <- purrr::compact(Data)
    df <- do.call(rbind, Data)
    PageData <- as.data.frame(df,
    stringsAsFactors = FALSE, check.names=FALSE)
} else {#Laser Case
    Cargo <- FixFirstLine(lines)
    NextLineIndex <- Cargo[[2]]
    TheFirstDetector <- Cargo[[1]]
    Segment <- lines[NextLineIndex:(LaserLines-1)]
    Segment <- c(TheFirstDetector, Segment)
    Data <- purrr::map(.f=SplitLines, .x=Segment)
    Data <- purrr::compact(Data)
    df <- do.call(rbind, Data)
    PageDataA <- as.data.frame(df,
    stringsAsFactors = FALSE, check.names=FALSE)

    Final <- grep("For Research Use Only", lines)
    Segment <- lines[(LaserLines+1):(Final-1)]
    Data <- purrr::map(.f=SplitLines, .x=Segment)
    Data <- purrr::compact(Data)
    df <- do.call(rbind, Data)
    PageDataB <- as.data.frame(df,
    stringsAsFactors = FALSE, check.names=FALSE)

    PageData <- list(PageDataA, PageDataB)   
}
} else {#High-speed page
    Final <- grep("For Research Use Only", lines)
    Segment <- lines[DetectorLines+2:Final]
    Data <- purrr::map(.f=SplitLines, .x=Segment)
    Data <- purrr::compact(Data)
    df <- do.call(rbind, Data)
    PageData <- as.data.frame(df,
    stringsAsFactors = FALSE, check.names=FALSE)
}

 return(PageData)
}
