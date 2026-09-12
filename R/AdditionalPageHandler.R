#' Internal for QC_ChorusPDF, handles second plus pages
#' 
#' @param x The iterated in pdf page as text lines
#' 
#' @importFrom purrr map compact 
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