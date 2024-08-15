#' Visualize cosine similarity of raw .fcs files to evaluate single color
#' controls.
#'
#' @param path The location to the folder where the Luciernaga .fcs files are stored
#' @param reference A path to a .csv file or a dataframe containing Fluorophore and Detector column
#' information for the panel.
#' @param stats Whether to use the median or mean for fluorescent intensity.
#' @param LinePlots Return this kind of plot, default is set to TRUE
#' @param CosinePlots Return this kind of plot, default is set to TRUE
#' @param StackedBarPlots Return this kind of plot, default is set to TRUE
#' @param HeatmapPlots Return this kind of plot, default is set to TRUE
#' @param RetainedType Whether the data.frame contains "raw" or "normalized" values
#'
#' @importFrom dplyr select
#' @importFrom dplyr pull
#'
#' @importFrom tidyr nest
#' @importFrom stringr str_detect
#' @importFrom flowCore exprs
#' @importFrom flowWorkspace keyword
#' @importFrom flowWorkspace load_cytoset_from_fcs
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom dplyr relocate
#' @importFrom dplyr select_if
#' @importFrom dplyr bind_rows
#' @importFrom dplyr group_by
#' @importFrom dplyr ungroup
#' @importFrom dplyr summarise_all
#' @importFrom tidyr nest
#' @importFrom tidyr unnest
#' @importFrom purrr map
#' @importFrom lsa cosine
#' @importFrom reshape2 melt
#' @import ggplot2
#' @importFrom stats as.dist
#' @importFrom stats hclust
#'
#'
#'
#' @return Visualized ggplots for each fluorophore.
#'
#' @export
#'
#' @examples NULL

LuciernagaReportFromFCS <- function(path, reference, stats = "median",
                                    LinePlots = TRUE, CosinePlots = TRUE,
                                    StackedBarPlots = TRUE, Heatmaps = TRUE,
                                    RetainedType){

  if (!is.data.frame(reference)){CSV <- read.csv(reference, check.names = FALSE)
  } else {CSV <- reference}

  internalstrings <- c("-A", ".", "_", " ")
  CSV$Fluorophore <- Luciernaga:::NameCleanUp(name=CSV$Fluorophore, removestrings=internalstrings)
  CSV$Detector <- Luciernaga:::NameCleanUp(name=CSV$Detector, removestrings=internalstrings)
  Variables <- CSV %>% select(Fluorophore) %>% pull(.)
  fcsfiles <- list.files(path, pattern=".fcs", full.names = TRUE)
  #x <- Variables[2]
  #inputfiles <- fcsfiles

  TheseFluorophores <- map(.x=Variables, .f=FluorophoreFilePresent, inputfiles = fcsfiles)
  TheseFluorophores <- Filter(Negate(is.null), TheseFluorophores)
  TheseFluorophores <- unlist(TheseFluorophores)
  #x <- TheseFluorophores[1]
  #data <- CSV
  #inputfiles = fcsfiles

  PlotTwist <- map(.x = TheseFlurophores, .f = FCSImport, data = CSV,
                   inputfiles = fcsfiles, RetainedType=RetainedType, stats=stats)
  #PlotTwist <- Filter(Negate(is.null), PlotTwist)
  #PlotTwist <- unlist( PlotTwist)

  return(PlotTwist)
}


#' Internal for LuciernagaReportFromFCS
#'
#' @param x Passed Fluorophore Name
#' @param inputfiles List of .fcs files from path
#'
#' @importFrom stringr str_detect
#' @noRd

FluorophoreFilePresent <- function(x, inputfiles){
  fcs_files <- inputfiles[str_detect(basename(inputfiles), x) &
                            str_detect(basename(inputfiles), ".fcs$")]
  if (length(fcs_files) > 0){return(x)}
}


#' Internal for LuciernagaReportFromFCS
#'
#' @param x A passed single cytoset object
#' @param Fluorophore The detector
#'
#' @importFrom flowCore keyword
#' @importFrom flowCore exprs
#' @importFrom dplyr mutate
#' @importFrom dplyr relocate
#'
#' @noRd

FCSImportFile <- function(x, Fluorophore, sample.name = "FILENAME"){
  filename <- keyword(x, sample.name)
  filename <- sub(".*\\\\", "", filename)
  filename <- sub(paste0(".*", Fluorophore), Fluorophore, filename)
  filename <- gsub(".fcs$", "", filename)
  rownames(filename) <- NULL

  #df <- exprs(x[[1]])
  df <- exprs(x)
  TheDF <- data.frame(df, check.names = FALSE)
  TheDF <- TheDF[,-grep("Time|FS|SC|SS|Original|W$|H$", names(TheDF))]
  colnames(TheDF) <- gsub("-A$", "", colnames(TheDF))
  DFNames <- TheDF %>% mutate(Cluster = filename) %>% relocate(Cluster,
                                                               .before = 1)
  return(DFNames)
}

#' Internal for LuciernagaReportFromFCS
#'
#' @param x Passed Fluorophore Name
#' @param data Passed data.frame of Fluorophore with Detectors
#' @param inputfiles List of .fcs files from path
#'
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#' @importFrom stringr str_detect
#' @importFrom flowWorkspace load_cytoset_from_fcs
#' @importFrom purrr map
#' @noRd


FCSImport <- function(x, data, inputfiles, RetainedType, stats){

  # For each Fluorophore
  TheDetector <- data %>% filter(Fluorophore %in% x) %>% pull(Detector)

  if (x %in% c("PE", "APC")){x <- paste0(x, "_")} #ExceptionHandling
  fcs_files <- inputfiles[str_detect(basename(inputfiles), x) &
                            str_detect(basename(inputfiles), ".fcs$")]
  if (x %in% c("PE_", "APC_")){x <- gsub("_", "", x)} #ExceptionHandling

  if (!length(fcs_files) == 0){

    cs <- load_cytoset_from_fcs(fcs_files, truncate_max_range = FALSE,
                                transformation = FALSE)
    thex <- x

    # Retrieve exprs data from each .fcs file and create cluster column
    # x <- cs[1]
    TheData <- map(.x = cs, .f = FCSImportFile, Fluorophore = thex) %>%
      bind_rows()

    # Removing Any Artificial Negatives Inserted By Luciernaga
    TheData <- TheData %>% mutate(Summed = rowSums(across(where(is.numeric)), na.rm = TRUE))

    TheData <- TheData %>% group_by(Summed) %>% dplyr::filter(n() <= 5) %>%
      ungroup() %>% dplyr::filter(!Summed == 0) %>% select(-Summed)

    # Return Summarized Data
    if (RetainedType == "normalized"){


    }


  }

}







