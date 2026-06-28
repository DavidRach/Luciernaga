#' Working from a completed template, subtracts the respective negative for
#'  each fluorophore
#' 
#' @param data The updated template
#' @param Normalize Default is TRUE, will normalize the subtracted data so
#' that scale is from 0 to 100. 
#' @param keepNegatives Default is FALSE, when TRUE, returns the negatives used
#' to subtract back to their original rows
#' 
#' @importFrom dplyr select rename_with left_join mutate across cur_column
#' pull rows_update filter bind_cols
#' @importFrom tidyselect all_of ends_with
#' 
#' @return A data.frame with fluorophore signatures
#' 
#' @export
#' 
#' @examples A <- 2 + 2
#' 
BackgroundSubtraction <- function(data, Normalize=TRUE, keepNegatives=FALSE){
  StartingDetector <- which(colnames(data) %in% "Negative") + 1
  DetectorColumns <- colnames(data)[StartingDetector:ncol(data)]

NegativeTheoretical <- UpdatedTemplate |>
  select(Fluorophore, all_of(DetectorColumns)) |>
  rename_with(~ paste0(.x, "_neg"), all_of(DetectorColumns))

BackgroundSubtracted <- UpdatedTemplate |>
  left_join(NegativeTheoretical, by = c("Negative" = "Fluorophore")) |>
  mutate(across(
    all_of(DetectorColumns),
    ~ . - get(paste0(cur_column(), "_neg"))
  )) |> select(-ends_with("_neg"))

Negatives <- BackgroundSubtracted |> pull(Negative) |> unique()

if (keepNegatives==TRUE){
    Final <- BackgroundSubtracted |>
    rows_update(
      UpdatedTemplate |>
        filter(Fluorophore %in% Negatives) |>
        select(Fluorophore, all_of(DetectorColumns)),
      by = "Fluorophore"
    )

} else {
    Final <- BackgroundSubtracted |> filter(!Fluorophore %in% Negatives)
}

if (Normalize == TRUE){
    x <- Final |> select(all_of(DetectorColumns))
    # x[x < 0] <- 0
    A <- do.call(pmax, x)
    x <- x/A
    x <- round(x, 3)

    Final <- Final |>  select(-all_of(DetectorColumns)) |>
       bind_cols(x)
}

return(Final)

}