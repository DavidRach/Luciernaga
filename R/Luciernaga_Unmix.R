#' OLS unmixing of a Gating Set object
#'
#' @param x A Gating Set object
#' @param controlData The matrix of single color controls generated
#'  from Luciernaga
#' @param sample.name The keyword containing the fcs file name
#' @param addon Additional addon to append to the new .fcs file name
#' @param removestrings A list of values to remove from name
#' @param subset A gating hierarchy level to sort cells at, expression
#'  values retrieved
#' from these
#' @param outpath The return folder for the .fcs files
#' @param Verbose For troubleshooting name after removestrings
#' @param PanelPath Location to a panel.csv containing correct order of
#'  fluorophores
#' @param returnType Whether to return "fcs" or "flowframe"
#' @param inverse.transform Default is FALSE, set to TRUE if data is
#'  already transformed
#' and needs to be reversed. 
#'
#' @importFrom flowWorkspace gs_pop_get_data
#' @importFrom flowCore exprs keyword write.FCS
#' @importFrom dplyr mutate select pull arrange bind_cols
#' @importFrom tidyselect where
#' @importFrom utils read.csv
#' @importFrom stats lsfit
#'
#' @return A new .fcs file with the new columns appended
#' 
#' @export
#'
#' @examples A <- 2 + 2
#' 
Luciernaga_Unmix <- function(x,
                              controlData,
                              sample.name,
                              removestrings,
                              Verbose,
                              addon,
                              subset = "root",
                              outpath,
                              PanelPath,
                              returnType = "fcs",
                              inverse.transform = FALSE) {

  if (length(sample.name) == 2) {
    first <- sample.name[[1]]
    second <- sample.name[[2]]
    first <- keyword(x, first)
    second <- keyword(x, second)
    name <- paste(first, second, sep = "_")
  } else {
    name <- keyword(x, sample.name)
  }

  name <- NameCleanUp(name, removestrings = removestrings)
  if (Verbose == TRUE) {
    message("After removestrings, name is ", name)
  }

  cs <- gs_pop_get_data(x, subset, inverse.transform = inverse.transform)
  Data <- exprs(cs[[1]])
  Data <- data.frame(Data, check.names = FALSE)

  OriginalColumnsVector <- colnames(Data)
  OriginalColumns <- colnames(Data)
  OriginalColumns <- data.frame(OriginalColumns, check.names = FALSE)
  OriginalColumnsIndex <- OriginalColumns %>%
    mutate(IndexLocation = 1:nrow(.)) # TODO: `.` refers to lhs (magrittr-only)

  Backups <- Data |> mutate(Backups = 1:nrow(Data)) |> select(Backups)

  StashedDF <- Data[, grep("Time|FS|SC|SS|Original|W$|H$", names(Data))]
  StashedDF <- cbind(Backups, StashedDF)

  TheSampleData <- Data[, -grep("Time|FS|SC|SS|Original|W$|H$", names(Data))]
  BackupNames <- colnames(TheSampleData)

  if (!is.data.frame(PanelPath)) {
    Panel <- read.csv(PanelPath, check.names = FALSE)
  } else {
    Panel <- PanelPath
  }

  CorrectColumnOrder <- Panel |> pull(Fluorophore)
  CorrectColumnOrder <- gsub("-A$", "", CorrectColumnOrder)

  if (any(controlData |> select(where(is.numeric)) > 1)) {
    Metadata <- controlData |> select(!where(is.numeric))
    Numerics <- controlData |> select(where(is.numeric))
    n <- Numerics
    n[n < 0] <- 0
    A <- do.call(pmax, n)
    Normalized <- n / A
    controlData <- bind_cols(Metadata, Normalized)
  }

  controlData$Fluorophore <- gsub("-A$", "", controlData$Fluorophore)

  NewControlData <- controlData |>
    arrange(match(Fluorophore, CorrectColumnOrder))
  Newest <- NewControlData |> pull(Fluorophore)

  if (!identical(CorrectColumnOrder, Newest)) {
    message(Newest)
    stop("Column Reordering Failed, printed order output for troubleshooting:")
  }

  TheControlData <- NewControlData[names(TheSampleData)]
  NewNames <- NewControlData |> select(Fluorophore)
  NewNames$Fluorophore <- paste0(NewNames$Fluorophore, "-A")
  NewNames <- NewNames |> pull(Fluorophore)
  Ligands <- NewControlData |> pull(Ligand)

  LeastSquares <- lsfit(x = t(TheControlData), y = t(TheSampleData),
                         intercept = FALSE)
  UnmixedData <- t(LeastSquares$coefficients)
  UnmixedData2 <- UnmixedData

  colnames(UnmixedData2) <- NewNames
  TheData <- cbind(StashedDF, UnmixedData2)
  TheData <- TheData |> select(-Backups)
  rownames(TheData) <- NULL

  new_fcs <- InternalUnmix(cs = cs, StashedDF = StashedDF, TheData = TheData,
                            Ligands = Ligands)
  # View(new_fcs@description)

  if (!is.null(addon)) {
    name <- paste0(name, addon)
  }

  AssembledName <- paste0(name, ".fcs")

  new_fcs@description$GUID <- AssembledName
  new_fcs@description$`$FIL` <- AssembledName

  if (is.null(outpath)) {
    outpath <- getwd()
  }

  fileSpot <- file.path(outpath, AssembledName)

  if (returnType == "fcs") {
    write.FCS(new_fcs, filename = fileSpot, delimiter = "#")
  } else {
    return(new_fcs)
  }
}