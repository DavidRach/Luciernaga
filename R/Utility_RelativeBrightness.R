#' Helper function splits Cluster into individual components
#'
#' @param x A data.frame containing column Cluster
#'
#' @return NULL
#' @export
#'
#' @examples NULL
Utility_RelativeBrightness <- function(x){
  Regular <- x %>% filter(!str_detect(Cluster, "-$"))

  if(nrow(Regular) > 0){
    Regular <- Regular %>% mutate(regular_split = str_split(as.character(Cluster), "-")) %>% relocate(regular_split, .after = Cluster)
    Regular <- Regular %>%  unnest_wider(regular_split, names_sep = "_") %>% rename_with(~ c("Detector1", "Detector2", "Detector3"), starts_with("regular"))

    Regular <- Regular %>% mutate(Detector1 = str_split(as.character(Detector1), "_"))
    Regular <- Regular %>%  unnest_wider(Detector1, names_sep = "_") %>% rename_with(~ c("Detector1", "Detector1Value"), starts_with("Detector1"))

    Regular <- Regular %>% mutate(Detector2 = str_split(as.character(Detector2), "_"))
    Regular <- Regular %>%  unnest_wider(Detector2, names_sep = "_") %>% rename_with(~ c("Detector2", "Detector2Value"), starts_with("Detector2"))

    Regular <- Regular %>% mutate(Detector3 = str_split(as.character(Detector3), "_"))
    Regular <- Regular %>%  unnest_wider(Detector3, names_sep = "_") %>% rename_with(~ c("Detector3", "Detector3Value"), starts_with("Detector3"))
  }

  Exceptions <- x %>% filter(str_detect(Cluster, "-$"))

  if(nrow(Exceptions) > 0){

    Exceptions$Cluster <- gsub("-$", "", Exceptions$Cluster)
    Exceptions <- Exceptions %>% mutate(regular_split = str_split(as.character(Cluster), "-")) %>% relocate(regular_split, .after = Cluster)
    Exceptions <- Exceptions %>%  unnest_wider(regular_split, names_sep = "_") %>% rename_with(~ c("Detector1", "Detector2"), starts_with("regular"))

    Exceptions <- Exceptions %>% mutate(Detector1 = str_split(as.character(Detector1), "_"))
    Exceptions <- Exceptions %>%  unnest_wider(Detector1, names_sep = "_") %>% rename_with(~ c("Detector1", "Detector1Value"), starts_with("Detector1"))

    Exceptions <- Exceptions %>% mutate(Detector2 = str_split(as.character(Detector2), "_"))
    Exceptions <- Exceptions %>%  unnest_wider(Detector2, names_sep = "_") %>% rename_with(~ c("Detector2", "Detector2Value"), starts_with("Detector2"))

    Exceptions <- Exceptions %>% mutate(Detector3 = rep(NA, nrow(Exceptions))) %>% relocate(Detector3, .after = Detector2Value)
    Exceptions$Detector3 <- as.character(Exceptions$Detector3)

    Exceptions <- Exceptions %>% mutate(Detector3Value = rep(NA, nrow(Exceptions))) %>% relocate(Detector3Value, .after = Detector3)
    Exceptions$Detector3Value <- as.character(Exceptions$Detector3Value)

  }

  Combined <- rbind(Regular, Exceptions)

  Combined$Detector1Value <- as.numeric(Combined$Detector1Value)
  Combined$Detector2Value <- as.numeric(Combined$Detector2Value)
  Combined$Detector3Value <- as.numeric(Combined$Detector3Value)

  Combined <- Combined %>% mutate(Brightness = rowSums(select(., Detector1Value, Detector2Value, Detector3Value), na.rm = TRUE)) %>% relocate(Brightness, .after = Cluster)

  Combined <- Combined %>% mutate(Detector1Raw = rep(NA, nrow(Combined)), Detector2Raw = rep(NA, nrow(Combined)), Detector3Raw = rep(NA, nrow(Combined))) %>% relocate(Detector1Raw, .after = Detector1Value) %>% relocate(Detector2Raw, .after = Detector2Value) %>% relocate(Detector3Raw, .after = Detector3Value)

  Combined <- Utility_Fill(Combined)

  return(Combined)
}
