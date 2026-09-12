
#' Internal for LuciernagaQC
#'
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom tidyselect all_of
#' @importFrom dplyr rename
#'
#' @return An internal value
#'
#' @noRd
LuciernagaSmallReport <- function( x, Data, RetainedType, ColsN,
    StartNormalizedMergedCol, EndNormalizedMergedCol, stats){

    if (RetainedType == "raw"){Data <- Data %>% filter(Cluster %in% x) %>%
      select(all_of(1:ColsN))
    Averaged <- AveragedSignature(Data, stats)
    }
    if (RetainedType == "normalized"){
      # Data <- Data %>% filter(Cluster %in% x) %>%
      # select(all_of(StartNormalizedMergedCol:EndNormalizedMergedCol))
      Data <- Data %>% filter(Cluster %in% x) %>%
        select(all_of(1:ColsN))
      Averaged <- AveragedSignature(Data, stats, normalize=TRUE)
      }
    Summary <- cbind(x, Averaged) %>% rename(Cluster = x)
    return(Summary)
}