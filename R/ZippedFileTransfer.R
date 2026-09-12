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