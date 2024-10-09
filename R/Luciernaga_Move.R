#' Transfer the .fcs files selected by LuciernagaTree to a new location
#'
#' @param x A vector containing fluorophore names found within data.
#' @param data The data.frame output of LuciernagaTree listing files to be moved.
#' @param input The path to the current storage location of the .fcs files.
#' @param output The path to the desired future storage location of the selected
#' files.
#'
#' @importFrom dplyr filter
#' @importFrom stringr str_detect
#' @importFrom dplyr pull
#'
#' @return No return, .fcs files are moved desired folder.
#' @export
#'
#' @examples NULL
Luciernaga_Move <- function(x, data, input, output){
  OriginalX <- x
  x <- gsub("-A", "", x)
  x <- gsub(" ", "", x)

  Internal <- data %>% dplyr::filter(str_detect(sample, fixed(x, ignore_case = TRUE)))

  if (x %in% c("PE", "APC")){

    Internal <- Internal %>% dplyr::filter(!str_detect(sample, "PE-|APC-|Per"))

    } #ExceptionHandling


  Fluor <- Internal %>% pull(sample) %>% unique()
  Clusterlet <- Internal %>% pull(Cluster) %>% unique()

  internalstrings <- c("_","-")
  Clusterlet <- Luciernaga::NameCleanUp(Clusterlet, removestrings = internalstrings)

  inputfiles <- list.files(input, full.names = TRUE)

  files_to_move <- inputfiles[str_detect(basename(inputfiles), fixed(Fluor, ignore_case = TRUE)) &
                                str_detect(basename(inputfiles), paste0("_", Clusterlet, "\\."))]

  file.copy(files_to_move, output)

}

