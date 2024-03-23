#' OLS unmixing of a Gating Set object
#'
#' @param x A Gating Set object
#' @param controlData The matrix of single color controls generated from Luciernaga
#' @param sample.name The keyword containing the fcs file name
#' @param addon Additional addon to append to the new .fcs file name
#' @param removestrings A list of values to remove from name
#' @param subset A gating hierarchy level to sort cells at, expression values retrieved from these
#' @param multiplier A number to scale the OLS coefficients by
#' @param outpath The return folder for the .fcs files
#'
#' @importFrom flowCore exprs
#' @importFrom flowCore write.FCS
#' @importFrom flowWorkspace keyword
#' @importFrom flowWorkspace gs_pop_get_data
#' @importFrom flowWorkspace realize_view
#' @importFrom flowWorkspace cf_append_cols
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom stats lsfit
#'
#'
#' @return A new .fcs file with the new columns appended
#' @export
#'
#' @examples Not at this moment.

LuciernagaUnmix <- function(x, controlData, sample.name, addon, removestrings, subset, multiplier, outpath){

  name <- keyword(x, sample.name)
  name <- NameCleanUp(name, removestrings=removestrings)

  cs <- gs_pop_get_data(x, subset)
  Data <- exprs(cs[[1]])
  Data <- data.frame(Data, check.names = FALSE)

  #For Future Column Reordering
  OriginalColumnsVector <- colnames(Data)
  OriginalColumns <- colnames(Data)
  OriginalColumns <- data.frame(OriginalColumns, check.names = FALSE)
  OriginalColumnsIndex <- OriginalColumns %>% mutate(IndexLocation = 1:nrow(.))

  #For Future Row Reordering
  Backups <- Data %>% mutate(Backups = 1:nrow(Data)) %>% select(Backups)

  #Stashing Away Time FSC SSC For Later Use
  StashedDF <- Data[,grep("Time|FS|SC|SS|Original|W$|H$", names(Data))]
  StashedDF <- cbind(Backups, StashedDF)

  #Consolidating Columns Going Forward
  TheSampleData <- Data[,-grep("Time|FS|SC|SS|Original|W$|H$", names(Data))]
  BackupNames <- colnames(TheSampleData)

  # Ordering the Single Color Control Columns to Match the Samples
  TheControlData <- controlData[names(expresionData)]

  #OLS
  LeastSquares <- lsfit(x = t(TheControlData), y = t(TheSampleData), intercept = FALSE)
  UnmixedData <- t(LeastSquares$coefficients)
  UnmixedData2 <- UnmixedData*multiplier

  #StashedResults <- summary(UnmixedData)
  #StashedExpresion <- summary(TheSampleData)

  cf <- realize_view(cs[[1]])
  NewCF <- cf_append_cols(cf, UnmixedData2)
  #keyword(NewCF)

  if (!is.null(addon)){name <- paste0(name, addon)
  }

  AssembledName <- paste0(name, ".fcs")
  FileName <- file.path(outpath, AssembledName)

  write.FCS(NewCF, filename=FileName)
}
