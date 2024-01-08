

#' Move the LuciernaTree selected fcs files to a new location
#'
#' @param x A vector of fluorophores in the data
#' @param data The data.frame output of LuciernagaTree
#' @param input The location of the fcs stored files
#' @param output The location to move the selected files
#'
#' @return NULL
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
