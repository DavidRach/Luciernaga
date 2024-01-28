

#' Transfer the .fcs files selected by LuciernagaTree to a new location
#'
#' @param x A vector containing fluorophore names found within data.
#' @param data The data.frame output of LuciernagaTree listing files to be moved.
#' @param input The path to the current storage location of the .fcs files.
#' @param output The path to the desired future storage location of the selected
#' files.
#'
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#' @importFrom stringr str_detect
#'
#' @return No return, .fcs files are moved desired folder.
#' @export
#'
#' @examples NULL
LuciernagaMove <- function(x, data, input, output){
  Internal <- data %>% filter(sample %in% x)

  Fluor <- Internal %>% pull(sample)
  Clusterlet <- Internal %>% pull(Cluster)

  inputfiles <- list.files(input, full.names = TRUE)

  files_to_move <- inputfiles[str_detect(basename(inputfiles), Fluor) &
                                str_detect(basename(inputfiles), Clusterlet)]

  file.copy(files_to_move, output)

}
