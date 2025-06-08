#' Takes a .csv with Fluorophore and Antigen columns, and generates
#' a corresponding Cytek Aurora matrix for inter/intra comparisons
#' 
#' @param data A path to the .csv or a data.frame object
#' @param NumberDetectors The corresponding number of detectors for the Cytek Aurora
#' @param returnType Default is plot
#' 
#' @importFrom utils read.csv
#' @importFrom dplyr pull
#' @importFrom dplyr filter
#' @importFrom dplyr group_by
#' @importFrom dplyr arrange
#' @importFrom dplyr desc
#' @importFrom dplyr slice
#' @importFrom dplyr select
#' @importFrom dplyr ungroup
#' @importFrom dplyr mutate
#' @importFrom stringr str_extract
#' @importFrom dplyr left_join
#' @importFrom dplyr across
#' @importFrom tidyselect everything
#' @importFrom tidyr replace_na
#' @importFrom purrr imap
#' @importFrom dplyr rename_with
#' @importFrom purrr reduce
#' @importFrom dplyr full_join
#' @importFrom stats setNames
#' @importFrom stringr str_replace
#' @importFrom tidyselect all_of
#' @importFrom gt gt
#' @importFrom gt tab_style
#' @importFrom gt cell_text
#' @importFrom gt cells_body
#' @importFrom gt cols_width
#' @importFrom gt cols_label
#' @importFrom gt px
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

#' Generates a UV laser template
#' 
#' @return An UV laser data.frame
#' 
#' @noRd
UVLaser <- function(){
 Laser <- tibble(
  Detector=c("UV1", "UV2", "UV3", "UV4", "UV5", "UV6", "UV7", "UV8",
  "UV9", "UV10", "UV11", "UV12","UV13","UV14","UV15","UV16"),
  Wavelength=c(373, 388, 428, 443, 458, 473, 514, 
  542, 582, 613, 662, 695, 718, 750, 781, 812)
 )
 return(Laser)
}

#' Generates a V laser template
#' 
#' @return An V laser data.frame
#' 
#' @noRd
VLaser <- function(){
  Laser <- tibble(
   Detector=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8",
   "V9", "V10", "V11", "V12","V13","V14","V15","V16"),
   Wavelength=c(428, 443, 458, 473, 508, 525,
   542, 581, 598, 615, 662, 695, 718, 750, 781, 812)
  )
  return(Laser)
 }

#' Generates a B laser template
#' 
#' @return An B laser data.frame
#' 
#' @noRd
BLaser <- function(){
  Laser <- tibble(
   Detector=c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8",
   "B9", "B10", "B11", "B12","B13","B14"),
   Wavelength=c(508, 525, 542, 581, 598, 615, 662,
   679, 695, 718, 738, 760, 781, 812)
  )
  return(Laser)
 }

#' Generates a YG laser template
#' 
#' @return An YG laser data.frame
#' 
#' @noRd
YGLaser <- function(){
  Laser <- tibble(
   Detector=c("YG1", "YG2", "YG3", "YG4", "YG5", "YG6", "YG7", "YG8",
   "YG9", "YG10"),
   Wavelength=c(577, 598, 615, 662, 679, 695, 718, 750, 781, 812)
  )
  return(Laser)
 }

#' Generates a R laser template
#' 
#' @return An R laser data.frame
#' 
#' @noRd
RLaser <- function(){
  Laser <- tibble(
   Detector=c("R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8"),
   Wavelength=c(662, 679, 695, 718, 738, 760, 781, 812)
  )
  return(Laser)
 }
