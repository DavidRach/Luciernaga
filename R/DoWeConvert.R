#' Internal for Utility_Concatinate
#'
#' @param data Selected metadata columns
#'
#' @noRd
DoWeConvert <- function(data) {

results <- data.frame(Column=character(), Numbers=logical(),
                        Letters=logical(), stringsAsFactors = FALSE)

for (i in names(data)) {
  if (is.character(data[[i]])) {
    NumbersPaa <- any(grepl("\\d", data[[i]]))
    LettersPaa <- any(grepl("[A-Za-z]", data[[i]]))

    results <- rbind(results, data.frame(Column=i, Numbers=NumbersPaa,
                                     Letters=LettersPaa,
                                     stringsAsFactors = FALSE))
  }
}

return(results)
}