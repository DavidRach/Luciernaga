#' Helper function splits Cluster into individual components
#'
#' @param x A data.frame containing column Cluster
#'
#' @importFrom dplyr mutate relocate across rename_with relocate
#'  select pull bind_rows
#' @importFrom stringr str_split
#' @importFrom tidyr unnest_wider starts_with
#' @importFrom tidyselect ends_with
#' @importFrom purrr map 
#'
#' @return A value to be determined later
#'
#' @noRd
RelativeBrightness <- function(x){
  Regular <- x

  if(nrow(Regular) > 0){
    Regular <- Regular %>% mutate(regular_split = str_split(
      as.character(Cluster), "-")) %>% relocate(regular_split, .after = Cluster)
    Regular <- Regular %>%  unnest_wider(regular_split, names_sep = "_")
    Regular <- Regular %>% mutate(
      across(starts_with("regular_split"), ~ str_split(as.character(.), "_")))
    Regular <- Regular %>%  unnest_wider(starts_with("Regular"), names_sep = "_")

    Regular <- Regular %>% rename_with(~paste0("Detector", seq_along(.)), ends_with("_1"))
    Regular <- Regular %>%
      rename_with(~paste0("Detector", seq_along(.), "Value"), ends_with("_2"))
  } else {stop("No retained Clusters at this minimalfcsccutoff")}

  Combined <- Regular %>% mutate(across(ends_with("Value"), as.numeric))

  Combined <- Combined %>% mutate(
    Brightness = rowSums(select(., ends_with("Value")), na.rm = TRUE)) %>%
    relocate(Brightness, .after = Cluster)

  Values <- names(select(Combined, ends_with("Value")))

  Combined <- Combined %>% mutate(across(all_of(Values), ~ NA_real_, .names = "{.col}Raw")) %>%
    relocate(ends_with("Raw"), .after = all_of(Values))

  colnames(Combined) <- gsub("ValueRaw", "Raw", colnames(Combined))
  colnames(Combined) <- gsub("-A", "",  colnames(Combined))

  TheClusters <- Combined %>% select(Cluster) %>% pull()

  #x <- TheClusters[1]
  Generated <- map(.x=TheClusters, data=Combined, .f=FillIterate) %>% bind_rows()

  return(Generated)
}