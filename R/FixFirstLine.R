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