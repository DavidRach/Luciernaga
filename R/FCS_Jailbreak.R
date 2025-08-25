#' Helpful function for extracting .fcs files from SpectroFlo .Zip folders. 
#' 
#' @param x Provide either the Zip folder path, or a folder containing a single Zip folder within. 
#' I have not yet implemented a handling condition for when two zipped folders are present!
#' @param Type Default is "Raw", alternate is "Unmixed", provide "All" for both.
#' @param FileType Default is "Reference", alternates are "Unstained", and "Samples"
#' @param outpath The folder in which to store the unzipped .fcs files
#' 
#' @importFrom stringr str_detect
#' @importFrom purrr walk
#' @importFrom utils unzip
#' 
#' @return Files transferred to the desired folder
#' 
#' @export
#' 
#' @examples A <- 2 + 2
FCS_Jailbreak <- function(x, Type="Raw", FileType="Reference", outpath){
    if (grepl("\\.zip$", x)){message("Zip file detected")
    } else {
        ZipFile <- list.files(x, pattern=".zip", full.names=TRUE)
        if (length(ZipFile == 1)){message("Folder with Zip file detected")
            x <- ZipFile[1]
        } else {stop("No Zip File Found")}
    }

    ZippedFolder <- x

    ZipContents <- tryCatch({unzip(x, list = TRUE)$Name}, 
                          error=function(e){
                          return(NULL)
                          })

    ZipContents <- ZipContents[stringr::str_detect(ZipContents, ".fcs")]

    if (Type == "Raw"){ZipContents <- ZipContents[stringr::str_detect(ZipContents, "Raw")]}

    if (Type == "Unmixed"){ZipContents <- ZipContents[stringr::str_detect(ZipContents, "Unmixed")]}

    if (FileType == "Reference"){
        ZipContents <- ZipContents[stringr::str_detect(ZipContents, "Reference")]}

    if (FileType == "Unstained"){
        ZipContents <- ZipContents[stringr::str_detect(ZipContents, "nstained")]}

    if (FileType == "Samples"){
        ZipContents <- ZipContents[!stringr::str_detect(ZipContents, "Reference")]
        ZipContents <- ZipContents[!stringr::str_detect(ZipContents, "nstained")]
        }

    purrr::walk(.x=ZipContents, .f= ~ ZippedFileTransfer(x=.x, ZippedFolder=ZippedFolder,  outpath=outpath),
     .progress=TRUE)

    #message("Done")
}


#' Internal for FCS_Jailbreak, copies .fcs files to a temp folder
#'  in route to outpath
#' 
#' @param x The file.path to the file
#' @param ZippedFolder The original Zipped folder path
#' @param outpath The desired storage location
#' 
#' @importFrom utils unzip
#' 
#' @return Nothing, just transfers the file. 
#' 
#' @noRd
ZippedFileTransfer <- function(x, ZippedFolder, outpath){
  Tempd <- tempfile()
  dir.create(Tempd)
  #list.files(Tempd)
  unzip(ZippedFolder, files = x, exdir = Tempd)
  #list.files(Tempd)

  CleanName <- basename(x)
  CleanName <- gsub("[/\\\\]", ";", CleanName)
  CleanName <- gsub("Raw;", "", CleanName)
  CleanName <- gsub("Unmixed;", "", CleanName)
  CleanName <- gsub(";", "-", CleanName)

  ExtractedPath <- file.path(Tempd, x)
  OutwardPath <- file.path(outpath, CleanName)
  file.copy(ExtractedPath, OutwardPath, overwrite = TRUE)
  #list.files(Tempd)
  unlink(Tempd, recursive = TRUE, force = TRUE)
}