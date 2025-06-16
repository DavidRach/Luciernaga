#' Processes experiment folder for Luciernaga Signatures, intended
#' for David's panel, not yet generalized
#' 
#' @param x File.path to a target folder
#' 
#' @importFrom stringr str_detect
#' @importFrom data.table fread
#' @importFrom openCyto gatingTemplate
#' @importFrom flowWorkspace load_cytoset_from_fcs
#' @importFrom flowWorkspace GatingSet
#' @importFrom openCyto gt_gating
#' 
#' @return Luciernaga_QC processed signatures to a new Luciernaga subfolder
#' 
#' @noRd
#' 
#' @examples 
#' A <- 2+2
HandlingFunction <- function(x){

  # Isolate unmixing control types
  files <- list.files(x, pattern=".fcs", full.names=TRUE)
  Unstained <- files[grep("Unstained", files)]
  BeadsUnstained <- Unstained[grep("Beads", Unstained )]
  CellsUnstained <- Unstained[!stringr::str_detect(
      Unstained, "Beads|Dead")]

  Stained <- files[!stringr::str_detect(files, "Unstained")]
  BeadsStained <- Stained[grep("Beads", Stained)]
  CellsStained <- Stained[!stringr::str_detect(Stained, "Beads")]

  # References
  FileLocation <- system.file("extdata", package = "Luciernaga")
  CellGates <- data.table::fread(file.path(
      path = FileLocation, pattern = 'Gates.csv'))
  CellTemplate <- gatingTemplate(CellGates)
  BeadGates <- data.table::fread(file.path(
      path = FileLocation, pattern = 'GatesBeads.csv'))
  BeadTemplate <- gatingTemplate(BeadGates)
  pattern = "AutofluorescentOverlaps.csv"
  AFOverlap <- list.files(path=FileLocation,
   pattern=pattern, full.names = TRUE)

  # Gating
  UnstainedCell_CS <- load_cytoset_from_fcs(CellsUnstained,
   truncate_max_range = FALSE, transform = FALSE)
  UnstainedCell_GS <- GatingSet(UnstainedCell_CS)
  gt_gating(CellTemplate, UnstainedCell_GS)
  UnstainedCellPlot <- Utility_GatingPlots(x=UnstainedCell_GS[[1]],
   sample.name = "GUID", removestrings = ".fcs",
   gtFile = CellGates, DesiredGates = NULL,
   outpath = NULL, returnType="patchwork")

  StainedCell_CS <- load_cytoset_from_fcs(CellsStained,
   truncate_max_range = FALSE, transform = FALSE)
  StainedCell_GS <- GatingSet(StainedCell_CS)
  gt_gating(CellTemplate, StainedCell_GS)
  StainedCellPlot <- Utility_GatingPlots(x=StainedCell_GS[[1]],
   sample.name = "GUID", removestrings = ".fcs",
   gtFile = CellGates, DesiredGates = NULL,
   outpath = NULL, returnType="patchwork")

  UnstainedBeads_CS <- load_cytoset_from_fcs(BeadsUnstained,
   truncate_max_range = FALSE, transform = FALSE)
  UnstainedBeads_GS <- GatingSet(UnstainedBeads_CS)
  gt_gating(BeadTemplate, UnstainedBeads_GS)
  UnstainedBeadPlot <- Utility_GatingPlots(x=UnstainedBeads_GS[[1]],
   sample.name = "GUID", removestrings = ".fcs",
   gtFile = BeadGates, DesiredGates = NULL,
   outpath = NULL, returnType="patchwork")

  StainedBeads_CS <- load_cytoset_from_fcs(BeadsStained,
   truncate_max_range = FALSE, transform = FALSE)
  StainedBeads_GS <- GatingSet(StainedBeads_CS)
  gt_gating(BeadTemplate, StainedBeads_GS)
  StainedBeadPlot <- Utility_GatingPlots(x=StainedBeads_GS[[1]],
   sample.name = "GUID", removestrings = ".fcs",
   gtFile = BeadGates, DesiredGates = NULL,
   outpath = NULL, returnType="patchwork")

  These <- c(UnstainedCellPlot, StainedCellPlot, UnstainedBeadPlot, StainedBeadPlot)

  OnLocation <- list.dirs(x, full.names=TRUE, recursive=TRUE)
  Presence <- list.files(OnLocation, pattern="Luciernaga", full.names=TRUE)
  if (length(Presence) == 0){dir.create(file.path(x, "Luciernaga"))}
  Presence <- list.files(OnLocation, pattern="Luciernaga", full.names=TRUE)

  Utility_Patchwork(These, filename="GatingReport", outfolder=Presence,
  thecolumns=1, therows=2)

  # Single-color Cells

  UnstainedCellAF <- Luciernaga_QC(x=UnstainedCell_GS[[1]],
      subsets="nonDebris", removestrings=".fcs", sample.name="GUID",unmixingcontroltype = "cells", Unstained = FALSE, 
      ratiopopcutoff = 0.001, Verbose = FALSE, AFOverlap = AFOverlap,
      stats = "median", ExportType = "fcs", SignatureReturnNow = TRUE,
      outpath = Presence, Increments=0.1, minimalfcscutoff=0.001,
      experiment="SingleColor",condition="Cells")

  Unstained <- purrr::map(.f=Luciernaga_QC, .x=UnstainedCell_GS,
      subsets="nonDebris", removestrings=".fcs", sample.name="GUID",unmixingcontroltype = "cells", Unstained = FALSE, 
      ratiopopcutoff = 0.001, Verbose = FALSE, AFOverlap = AFOverlap,
      stats = "median", ExportType = "fcs", SignatureReturnNow = FALSE,
      outpath = Presence, Increments=0.1, minimalfcscutoff=0.001,
      experiment="SingleColor",condition="Cells")

  Colors <- purrr::map(.f=Luciernaga_QC, .x=StainedCell_GS,
      subsets="nonDebris", removestrings=".fcs", sample.name="GUID",unmixingcontroltype = "cells", Unstained = FALSE, 
      ratiopopcutoff = 0.001, Verbose = FALSE, AFOverlap = AFOverlap,
      stats = "median", ExportType = "fcs", SignatureReturnNow = FALSE,
      outpath = Presence, Increments=0.1, minimalfcscutoff=0.001,
      experiment="SingleColor",condition="Cells", CellAF=UnstainedCellAF)
}


#' Internal, used as basic html framework InteractiveLuciernaga
#' 
#' @param plot_list A list of plot objects
#' 
#' @importFrom htmltools div
#' 
#' @noRd
grid_layout <- function(plot_list, ncol = 2) {
  rows <- split(plot_list, ceiling(seq_along(plot_list) / ncol))
  html <- lapply(rows, function(row_plots) {
    div(style = "display: flex; justify-content: center;", 
        lapply(row_plots, function(p) {
          div(style = "width: 50%; padding: 10px;", p)
        })
    )
  })
  return(html)
}


#' Generates the plots for InteractiveLuciernaga
#' 
#' @param TheTarget The fluorophore being selected LuciernagaQC outputs
#' @param TheCleanStained The internal list of file.paths to the the single colors
#' @param returnType Default plotly
#' 
#' @importFrom stringr str_detect
#' @importFrom plotly ggplotly
#' 
#' @return A list of plotly objects for assembling
#' 
#' @noRd
#' 
#' @examples
#' A <- 2+2
#' 
LuciernagaCheck <- function(TheTarget, TheCleanStained, returnType="plotly"){

  TheseFluorSigs <- TheCleanStained[stringr::str_detect(TheCleanStained, paste0(TheTarget, " "))]
  
  Signatures <- Luciernaga_FolderSignatures(FolderPath=TheseFluorSigs, sample.name="GUID", fluorophore.name="", 
      StringRemoval=c(" (Cells)", ".fcs"))
  
  MainDetector <- colSums(Signatures == 1, na.rm = TRUE)
  MainDetector <- names(which.max(MainDetector))
  
  LinePlot <- QC_ViewSignature(x=NULL, data=Signatures, Normalize = FALSE, columnname="Sample", legend=TRUE)
  
  BrightnessPlot <- Luciernaga_FolderBrightness(FolderPath=TheseFluorSigs, sample.name="GUID",
      StringRemoval=c(" (Cells)", ".fcs"), fluorophore.name="", returnType = "plot",
      detector = MainDetector, PanelCuts=c(0.1,1))
  
  if (returnType == "plotly"){
      LinePlot <- plotly::ggplotly(LinePlot)
      BrightnessPlot <- plotly::ggplotly(BrightnessPlot)  
  }
  
  TheList <- list(LinePlot, BrightnessPlot)
  
  return(TheList)
}


#' Makes the website
#' 
#' @param x The file path to the parent folder
#' @param outpath The desired storage location
#' @param filename The desired file name
#' 
#' @importFrom stringr str_detect
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom dplyr group_by
#' @importFrom dplyr slice
#' @importFrom dplyr desc
#' @importFrom dplyr arrange
#' @importFrom dplyr mutate
#' @importFrom dplyr pull
#' @importFrom stringr str_extract
#' @importFrom purrr flatten
#' @importFrom purrr map
#' @importFrom htmltools tagList
#' @importFrom htmltools h1
#' @importFrom htmltools save_html
#' 
#' @return A htmlwebpage to desired location
#' 
#' @noRd
#' 
#' @examples
#' A <- 2 + 2 
GenerateInteractives <- function(x, outpath=NULL, filename=NULL){
  TheHoller <- file.path(x, "Luciernaga")
  TheFiles <- list.files(TheHoller, full.names=TRUE)
  TheFiles <- TheFiles[!stringr::str_detect(TheFiles, "GatingReport.pdf")]
  
  TheUnstaineds <- TheFiles[grep("Unstained", TheFiles)]
  RemoveThese <- c("Unstained (Cells)_", ".fcs", "PMA")
  TheseTypes <- NameCleanUp(basename(TheUnstaineds), removestrings=RemoveThese)
  MainAFs <- sub("10.*", "", TheseTypes) |> unique()
  
  if ("V" %in% MainAFs){MainAFs[MainAFs == "V"] <- "V10"}
  if ("UV" %in% MainAFs){MainAFs[MainAFs == "UV"] <- "UV10"}
  if ("B" %in% MainAFs){MainAFs[MainAFs == "B"] <- "B10"} 
  if ("YG" %in% MainAFs){MainAFs[MainAFs == "YG"] <- "YG10"} 
  
  TheStained <- TheFiles[!stringr::str_detect(TheFiles, "Unstained")]
  WatchForThese <- paste0("_", MainAFs, collapse = "|")
  
  TheCleanStained <- TheStained[!stringr::str_detect(TheStained, WatchForThese)]
  TheseDudes <- sub(" \\(Cells\\).*", "", basename(TheCleanStained)) |> unique()
  TheseDudettes <- sub("^[^ ]+ ", "", TheseDudes)
  TheseDudettes <- TheseDudettes[!stringr::str_detect(TheseDudettes, "Unstim")] |> unique()
  
  TheRefs <- Luciernaga:::InstrumentReferences(NumberDetectors=64)
  
  Ranking <- TheRefs |> dplyr::filter(Fluorophore %in% TheseDudettes) |>
       group_by(Fluorophore) |> arrange(desc(AdjustedY)) |> slice(1) |>
       select(-Instrument, -AdjustedY)
  
  Order <- c("UV", "V", "B", "YG", "R")
  
  Sequence <- Ranking |> mutate(
      prefix = stringr::str_extract(Detector, "^[A-Z]+"),
      num = as.numeric(stringr::str_extract(Detector, "\\d+")),
      group_order = match(prefix, Order)
    ) |> arrange(group_order, num) |> pull(Fluorophore)
  
  ThePlots <- map(.x=Sequence, .f=LuciernagaCheck, TheCleanStained=TheCleanStained,
   returnType="plotly")
  
  ThePlots <- flatten(ThePlots)
  
  if (is.null(outpath)){outpath <- getwd()}
  if (is.null(filename)){filename <- basename(x)}
  
  html_content <- tagList(
    h1("Experiment ", filename),
    grid_layout(ThePlots, ncol = 2)
  )
  
  Assembled <- paste0(filename, ".html")
  StorageLocation <- file.path(outpath, StorageLocation)
  
  # Save to an HTML file
  save_html(html_content, file = StorageLocation )
}