#' Checks and filters .fcs files into compatible sublist to
#' pass on to CytoSet without crashing out due to mismatch
#'
#' @param files A list containing fcs files file.paths
#'
#' @importFrom purrr map
#' @importFrom dplyr bind_rows mutate row_number relocate pull
#'
#' @return A list containing lists of compatible fcs file paths
#'
#' @export
CytosetScreen <- function(files) {

  Objects <- map(.x=files, .f=CytoSetScreenInternal) |>
    bind_rows()

  Objects <- Objects |> mutate(Iteration = row_number()) |>
    relocate(Iteration, .before=1)

  TheTable <- data.frame(table(Objects$ID))
  colnames(TheTable)[[1]] <- "Identity"
  TheIDs <- TheTable |> pull(Identity)

  ListOfList <- map(.x=TheIDs, .f=ListLocationFind, data=Objects,
      TheList=files)

  return(ListOfList)
}