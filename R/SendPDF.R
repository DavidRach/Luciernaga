#' Internal, sends patchwork objects to pdf
#' 
#' @param x The list of patchwork plots
#' @param outpath The file.path to desired storage location
#' @param filename The desired file name
#' @param width Default 7
#' @param height Default 9
#' 
#' @importFrom grDevices dev.off pdf
#' 
#' @return A pdf to the desired location
#' 
#' @noRd 
SendPDF <- function(x, outpath, filename, width = 7, height = 9) {
  TheName <- paste0(filename, ".pdf")
  StorageLocation <- file.path(outpath, TheName)

  pdf(file = StorageLocation, width = width, height = height)
  print(x)
  dev.off()
}