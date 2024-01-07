#' Internal function for relative brightness
#'
#' @param x The data.frame object passed by relative brightness
#'
#' @return NULL
#' @export
#'
#' @examples NULL
Utility_Fill <- function(x){
  colnames(x) <- gsub("-A", "",  colnames(x))

  TheClusters <- x %>% select(Cluster) %>% pull()

  Generated <- data.frame()

  for(k in TheClusters){IndividualCluster <- x %>% filter(Cluster %in% k)
  IndividualCluster <- data.frame(IndividualCluster, check.names = FALSE)
  Detectors <- IndividualCluster %>% select(Detector1, Detector2, Detector3) %>% unlist() %>% as.vector()

  #Detectors <- Detectors[!is.na(Detectors)]

  D1 <- Detectors[[1]]
  D2 <- Detectors[[2]]
  D3 <- Detectors[[3]]

  D1Value <- IndividualCluster %>% pull(all_of(D1))
  D2Value <- IndividualCluster %>% pull(all_of(D2))

  if (!is.na(D3)){D3Value <- IndividualCluster %>% pull(all_of(D3))
  FilledCluster <- IndividualCluster %>% mutate(Detector1Raw = coalesce(Detector1Raw, D1Value), Detector2Raw = coalesce(Detector2Raw, D2Value),
                                                Detector3Raw = coalesce(Detector3Raw, D3Value))
  } else {FilledCluster <- IndividualCluster %>% mutate(Detector1Raw = coalesce(Detector1Raw, D1Value), Detector2Raw = coalesce(Detector2Raw, D2Value))}

  Generated <- rbind(Generated, FilledCluster)
  }

  return(Generated)

}
