#' Carries out the initial unmixing using just the single-color controls own reference signature. 
#' 
#' @param matrix The signature matrix for the respective single-color unmixing controls
#' @param matrix_columnname The column name used to match respective metadata in the GatingSet. 
#' Default is Fluorophore.
#' @param gs The GatingSet containing the raw single-color unmixing controls
#' @param subset The gate with the events you want to unmix. Default is root. 
#' @param outpath The location to store the unmixed .fcs files. 
#' @param inverse.transform Default equals TRUE
#' @param sample.name Keyword(s) to grab from the cytoset for renaming
#' 
#' @importFrom flowWorkspace pData
#' @importFrom dplyr select pull mutate row_number filter
#' @importFrom tidyselect all_of
#' @importFrom purrr walk
#' 
#' @return Unmixed single-color unmixing controls for use in stain reduction calculations. 
#' 
#' @export
#' 
#' @examples A <- 2 + 2
#' 
SingleUnmix <- function(matrix,
                         matrix_columnname = "Fluorophore",
                         gs,
                         subset = "root",
                         outpath,
                         inverse.transform = TRUE,
                         sample.name = "TUBENAME") {

  # Metadata Fluors
  pd <- pData(gs)
  Location <- pd |>
    select(all_of(matrix_columnname)) |>
    mutate(Row = row_number())
  TheseFluorophores <- Location |> pull(matrix_columnname)

  # Corrected Matrix Order
  RetainedMatrix <- matrix |> filter(Fluorophore %in% TheseFluorophores)
  new_order <- match(TheseFluorophores, RetainedMatrix$Fluorophore)
  RetainedMatrix <- RetainedMatrix[new_order, ]

  FolderName <- file.path(outpath, "SingleUnmix")
  if (!dir.exists(FolderName)) {
    dir.create(FolderName)
  }

  walk(.x = TheseFluorophores, .f = SingleUnmixIterator, gs = gs,
       matrix = RetainedMatrix, outpath = FolderName, subset = subset,
       inverse.transform = inverse.transform, sample.name = sample.name)

  FolderName2 <- file.path(outpath, "AllUnmix")
  if (!dir.exists(FolderName2)) {
    dir.create(FolderName2)
  }

  walk(.x = TheseFluorophores, .f = AllUnmixIterator, gs = gs,
       matrix = RetainedMatrix, outpath = FolderName2, subset = subset,
       inverse.transform = inverse.transform, sample.name = sample.name)

  FolderName3 <- file.path(outpath, "FMO")
  if (!dir.exists(FolderName3)) {
    dir.create(FolderName3)
  }

  walk(.x = TheseFluorophores, .f = FMOUnmixIterator, gs = gs,
       matrix = RetainedMatrix, outpath = FolderName3, subset = subset,
       inverse.transform = inverse.transform, sample.name = sample.name)

  message("Done!")
}