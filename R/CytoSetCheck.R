#' Checks and filters .fcs files into compatible sublist to 
#' pass on to CytoSet without crashing out due to mismatch 
#' 
#' @param files A list containing fcs files file.paths
#' 
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
#' @importFrom dplyr mutate
#' @importFrom dplyr row_number
#' @importFrom dplyr relocate
#' @importFrom dplyr pull
#' 
#' @return A list containing lists of compatible fcs file paths
#' 
#' @export
CytosetScreen <- function(files){
  Objects <- map(.x=files, .f=CytoSetScreenInternal) |>
       bind_rows()

  Objects <- Objects %>% mutate(Iteration = row_number()) |>
       relocate(Iteration, .before=1)

  TheTable <- data.frame(table(Objects$ID))
  colnames(TheTable)[[1]] <- "Identity"
  TheIDs <- TheTable |> pull(Identity)

  ListOfList <- map(.x=TheIDs, .f=ListLocationFind, data=Objects, TheList=files)
  return(ListOfList)
}

#' Internal for CytoSetCheck, generates ID
#' 
#' @param x Iterated in path to individual fcs files
#' 
#' @importFrom flowWorkspace load_cytoset_from_fcs
#' 
#' @return A concatinated string corresponding column name order
#' 
#' @noRd
CytoSetScreenInternal <- function(x){
  TheCytoSet <- load_cytoset_from_fcs(files=x)
  Values <- colnames(TheCytoSet[[1]])
  ID <- paste(Values, collapse = " ")
  ID <- data.frame(ID)
  return(ID)
}

#' Internal for CytoSet check, actual filtering of subset list
#' 
#' @param x The iterated identity
#' @param data Intermediate data.frame containing RowNumber column
#' @param TheList The original list of fcs files to be sorted from
#' 
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#' 
#' @return A compatible subset list
#' 
#' @noRd
ListLocationFind <- function(x, data, TheList){
  TheRows <- data |> filter(ID %in% x) |> pull(Iteration)
  TheSubset <- TheList[TheRows]
  return(TheSubset)
}