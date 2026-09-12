#' Internal for CytekQCFilePrep
#'
#' @importFrom purrr map
#' @importFrom dplyr select rename
#' @importFrom tidyr starts_with
#' @importFrom lubridate mdy_hms mdy_hm
#'
#' @return An internal value
#'
#' @noRd
Internal_ChangeCalcs <- function(x, y, TheData){
  xx <- TheData %>% select(all_of(x))
  yy <- TheData %>% select(all_of(y))
  FinalValue <- xx[nrow(xx),]
  z <- xx[2:nrow(xx),]
  z <- c(z, FinalValue)
  Prelim <- xx %>% mutate(TheSubtract = z)
  Result <- Prelim %>% mutate(Crazy = .[[1]] - .[[2]])
  FinalResult <- Result %>% select(-2)
  FinalResult <- FinalResult %>% bind_cols(yy)
  FinalResult <- FinalResult %>% mutate(FlagCrazy = .[[3]])
  DiffColName <- paste0("Change_", colnames(FinalResult)[1])
  DiffFlagColName <- gsub("Flag-", "Flag-Change_", colnames(FinalResult)[3])
  colnames(FinalResult)[2] <- DiffColName
  colnames(FinalResult)[4] <- DiffFlagColName
  #FinalResult

  return(FinalResult)
}