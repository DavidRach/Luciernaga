#' Ease-of-life wrapper for Luciernaga_Summary, orchestrates passing of the
#' respective indices and gates to screen single-color unmixing controls
#' 
#' @param template Template designating single-color gates to screen
#' @param gs The GatingSet
#' @param AFOverlap See Luciernaga Vignette, default NULL falls back to the 
#' default shipped within Luciernaga extdata. 
#' @param externalAF_gs Default NULL, if you wish to provide an external AF for
#' subtraction, provide the gs name (and next two arguments). Leaving NULL will grab
#' the parent to provided gate on assumption both positive and negative cells will
#' be present in that gate. 
#' @param externalAF_gs_index The index value to select the unstained in combination
#'  with previous argument
#' @param externalAF_gs_gate The gate corresponding to the unstained population you
#' want to use for background subtraction
#' @param excludeThese Excludes columns containing values, leaving just the detector
#' columns. The default is set to "FSC|SSC|Time|-H|-W" 
#' @param inverse.transform Whether to inverse.transform, default is TRUE
#' @param GuessSimilar Default FALSE, will attempt to match to
#' similar fluorophores in the library
#' @param Unstained Default FALSE, set to TRUE if sample is unstained (and subtraction 
#' therefore is not needed)
#' @param NumberDetectors Default 64 (5-laser Cytek Aurora)
#' 
#' @importFrom dplyr select filter pull
#' @importFrom utils read.csv
#' @importFrom flowWorkspace sampleNames
#' @importFrom purrr map2
#' 
#' @return A list of lists containing data.frame and plots for all the single-color
#' reference controls designated in the template. 
#' 
#' @export 
#' 
#' @examples A <- 2 + 2
#' 
LuciernagaIntegration <- function(template, gs, AFOverlap=NULL,
 externalAF_gs = NULL, externalAF_gs_index=NULL, externalAF_gs_gate=NULL,
 excludeThese="FSC|SSC|Time|-H|-W", inverse.transform=TRUE,
 GuessSimilar=FALSE, Unstained=FALSE, NumberDetectors=64){

  if (is.null(AFOverlap)){
    FileLocation <- system.file("extdata", package = "Luciernaga")
    pattern = "AutofluorescentOverlaps.csv"
    AFOverlap <- list.files(path=FileLocation, pattern=pattern,
     full.names = TRUE)
    AFOverlap <- read.csv(AFOverlap, check.names=FALSE)
  } else {
    if(!is.data.frame(AFOverlap)){
      AFOverlap <- read.csv(AFOVerlap, check.names=FALSE)
      } else { # No Intervention needed
      }
  }

  These <- template |> dplyr::select(name, Fluorophore, Detector)
  DetectorsPresent <- These |> dplyr::filter(!is.na(Detector) & Detector != "")

  gs_index <- match(DetectorsPresent$name, sampleNames(gs))
  DetectorsPresent$gs_index <- gs_index

  GatesToAdd <- DetectorsPresent |> dplyr::pull(Fluorophore)
  SpecimenIndeces <- DetectorsPresent |> dplyr::pull(gs_index)

  # sampleNames(gs)
  # sampleNames(gs[SpecimenIndeces[14]])
  # x <- GatesToAdd[1]
  # y <- SpecimenIndeces[1]

  ListOfLists <- map2(.x=GatesToAdd, .y=SpecimenIndeces, .f=Luciernaga_Summary,
  gs=gs, externalAF_gs=externalAF_gs, externalAF_gs_index=externalAF_gs_index,
  externalAF_gs_gate=externalAF_gs_gate, GuessSimilar=GuessSimilar,
  NumberDetectors=NumberDetectors, excludeThese=excludeThese, Unstained=Unstained,
  inverse.transform=inverse.transform, AFOverlap=AFOverlap)

  return(ListOfLists)
}

