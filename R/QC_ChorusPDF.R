#' Processes BD Chorus QC PDF files into .csv files
#'
#' @param x A file.path to the desired QC pdf
#' @param returnPreference When both modes present, whether
#'  to return QC values from Imaging or High-Speed settings
#' @param returnType Default data, alternative is csv
#' @param outpath When return type is csv, file.path to store the csv
#'
#' @importFrom pdftools pdf_text
#' @importFrom purrr map
#' @importFrom utils write.csv
#'
#' @export
#'
#' @examples A <- 2+2
#'
QC_ChorusPDF <- function(x, returnPreference = "Imaging",
                          returnType = "data", outpath = NULL) {
  text <- pdftools::pdf_text(x)
  NumberPages <- length(text)
  FirstPageCargo <- FirstChorusPage(x = text[1])
  Metadata <- FirstPageCargo[[1]]
  FirstPage <- FirstPageCargo[[2]]

  Works <- purrr::map(.x = text[2:NumberPages], .f = AdditionalPageHandler)
  Dataset <- Consolidator(x = Works, Metadata = Metadata,
    FirstPage = FirstPage, returnPreference = returnPreference)

  if (returnType != "data") {
    if (is.null(outpath)) {
      outpath <- getwd()
    }
    PDFName <- Dataset$PDFName |> unique()
    if (returnPreference != "Imaging") {
      AppendValue <- paste0(returnPreference, ".csv")
    } else {
      AppendValue <- ".csv"
    }
    PDFName <- gsub(".pdf", AppendValue, PDFName)
    StorageLocation <- file.path(outpath, PDFName)
    write.csv(Dataset, StorageLocation, row.names = FALSE)
  }

  return(Dataset)
}