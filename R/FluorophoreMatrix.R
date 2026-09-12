#' Takes a .csv with Fluorophore and Antigen columns, and generates
#' a corresponding Cytek Aurora matrix for inter/intra comparisons
#' 
#' @param data A path to the .csv or a data.frame object
#' @param NumberDetectors The corresponding number of detectors for the Cytek Aurora
#' @param returnType Default is plot
#' 
#' @importFrom utils read.csv
#' @importFrom dplyr pull filter group_by arrange desc slice
#'  select ungroup mutate left_join across rename_with full_join
#' @importFrom stringr str_extract str_replace
#' @importFrom tidyselect everything
#' @importFrom tidyr replace_na
#' @importFrom purrr imap reduce
#' @importFrom stats setNames
#' @importFrom tidyselect all_of
#' @importFrom gt gt tab_style cell_text cells_body cols_width
#'  cols_label px
#' 
#' @return An assembled plot object for visualization
#' 
#' @export
#' 
#' @examples
#' A <- 2 + 2
FluorophoreMatrix <- function(data, NumberDetectors, returnType = "plot"){
  Vaiya <- InstrumentReferences(NumberDetectors = NumberDetectors)

  if (!is.data.frame(data)){Data <- read.csv(data, check.names=FALSE)
  } else {Data <- data}

  Fluorophores <- Data |> pull(Fluorophore)

  Subset <- Vaiya |> filter(Fluorophore %in% Fluorophores)

  Locations <- Subset |> group_by(Fluorophore) |>
    arrange(desc(AdjustedY)) |> slice(1) |>
    select(Fluorophore, Detector) |> ungroup()

  Order <- c("UV", "V", "B", "YG", "R")

  Sequence <- Locations |> mutate(
    prefix = str_extract(Detector, "^[A-Z]+"),
    num = as.numeric(str_extract(Detector, "\\d+")),
    group_order = match(prefix, Order)
  ) |> arrange(group_order, num)

  Lasers <- Sequence |> pull(prefix) |> unique()

  Dataset <- left_join(Sequence, Data, by="Fluorophore")

  if (returnType == "data"){return(Dataset)}

  TheList <- list()

  if(any(Lasers %in% "UV")){
    TheUV <- UVLaser()
    Selection <- Dataset |> filter(prefix %in% "UV") |>
      select(Fluorophore, Detector, Antigen)
    New <- left_join(TheUV, Selection, by="Detector")
    NewUV <- New %>% mutate(across(everything(), ~ replace_na(., "")))
    TheList <- append(TheList, list(NewUV))
  }
  if(any(Lasers %in% "V")){
    TheV <- VLaser()
    Selection <- Dataset |> filter(prefix %in% "V") |>
      select(Fluorophore, Detector, Antigen)
    New <- left_join(TheV, Selection, by="Detector")
    NewV <- New %>% mutate(across(everything(), ~ replace_na(., "")))
    TheList <- append(TheList, list(NewV))
  }
  if(any(Lasers %in% "B")){
    TheB <- BLaser()
    Selection <- Dataset |> filter(prefix %in% "B") |>
      select(Fluorophore, Detector, Antigen)
    New <- left_join(TheB, Selection, by="Detector")
    NewB <- New %>% mutate(across(everything(), ~ replace_na(., "")))
    TheList <- append(TheList, list(NewB))
  }
  if(any(Lasers %in% "YG")){
    TheYG <- YGLaser()
    Selection <- Dataset |> filter(prefix %in% "YG") |>
      select(Fluorophore, Detector, Antigen)
    New <- left_join(TheYG, Selection, by="Detector")
    NewYG <- New %>% mutate(across(everything(), ~ replace_na(., "")))
    TheList <- append(TheList, list(NewYG))
  }
  if(any(Lasers %in% "R")){
    TheR <- RLaser()
    Selection <- Dataset |> filter(prefix %in% "R") |>
      select(Fluorophore, Detector, Antigen)
    New <- left_join(TheR, Selection, by="Detector")
    NewR <- New %>% mutate(across(everything(), ~ replace_na(., "")))
    TheList <- append(TheList, list(NewR))
  }

  NamedList <- imap(TheList, function(df, name) {
    df %>% rename_with(~ paste0(., "_", name), .cols = -Wavelength)
  })

  Combined <- reduce(NamedList, full_join, by = "Wavelength") |>
  arrange(Wavelength)
  
  Combined1 <- Combined %>% mutate(across(everything(), ~ replace_na(., "")))

  Bolded <- names(Combined1)[grepl("^Fluorophore|^Antigen", names(Combined1))]
  colnames(Combined1) <- gsub("Detector", "D", colnames(Combined1))
  colnames(Combined1) <- gsub("Wavelength", "WV", colnames(Combined1))
  Narrow <- names(Combined1)[grepl("^(D|W)", names(Combined1), ignore.case = TRUE)]

  Cleaned <- names(Combined1) %>%
    setNames(., str_replace(., "_[1-5]$", ""))

  GTed <- Combined1 |> gt() |> tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(columns = all_of(Bolded))
  )
  
  #GTed <- GTed |> cols_width(all_of(Narrow) ~ px(30)) 

  Clean <- str_replace(names(Combined1), "_[0-9]+$", "")
  Labels <- setNames(Clean, names(Combined1))
  GTed <- GTed |> cols_label(.list = Labels)

  return(GTed)
}