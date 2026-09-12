#' Internal for LuciernagaReportFromFCS
#'
#' @param x Passed Fluorophore Name
#' @param inputfiles List of .fcs files from path
#'
#' @importFrom stringr str_detect
#'
#' @return An internal value
#'
#' @noRd
FluorophoreFilePresent <- function(x, inputfiles){
  fcs_files <- inputfiles[str_detect(basename(inputfiles), x) &
                            str_detect(basename(inputfiles), ".fcs$")]
  if (x %in% c("PE", "APC")){
      x <- paste0(x, "-")
      fcs_files <- fcs_files[!str_detect(basename(fcs_files), x)]
      x <- gsub("-", "", x)
    } #ExceptionHandling

  if (length(fcs_files) > 0){return(x)}
}