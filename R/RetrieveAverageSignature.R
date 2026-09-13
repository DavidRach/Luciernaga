#' Helper function, retrieves the average signature for each gate specified in
#' the template
#' 
#' @param template The template data.frame
#' @param gs A gating set
#' @param inverse.transform Whether to inverse.transform, default is set to TRUE
#' @param excludeThese Removes columns whose name contain these values via str_detect, 
#' default set for "FSC|SSC|Time|-H" which is typical for Cytek Aurora raw .fcs files
#' 
#' @importFrom dplyr select filter pull bind_rows left_join
#' @importFrom flowWorkspace sampleNames
#' @importFrom purrr map2
#' 
#' @return An updated template data.frame containing the average for each detector
#' 
#' @export
#' 
#' @examples A <- 2 + 2
#' 
RetrieveAverageSignature <- function(template,
                                      gs,
                                      inverse.transform = TRUE,
                                      excludeThese = "FSC|SSC|Time|-H") {

  These <- template |> select(name, Fluorophore, Detector)
  DetectorsPresent <- These |> filter(!is.na(Detector) & Detector != "")

  gs_index <- match(DetectorsPresent$name, sampleNames(gs))
  DetectorsPresent$gs_index <- gs_index

  GatesToAdd <- DetectorsPresent |> pull(Fluorophore)
  SpecimenIndeces <- DetectorsPresent |> pull(gs_index)

  Signatures <- map2(
    .x = SpecimenIndeces,
    .y = GatesToAdd,
    .f = \(x, y) GateExprsIterated(
      gs = gs,
      sample = x,
      gate = y,
      inverse.transform = inverse.transform,
      excludeThese = excludeThese
    )
  )

  TheSignatureMatrix <- Signatures |> bind_rows()

  UpdatedTemplate <- left_join(template, TheSignatureMatrix, by = "Fluorophore")
  return(UpdatedTemplate)
}