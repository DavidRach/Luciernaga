#' Internal for AutofluorescenceShop, parses individual .fcs file from
#' the GatingSet
#' 
#' @param x The iterated GatingSet object passed from AutofluorescenceShop
#' @param outpath Desired storage location for the AF Tags
#' @param ExperimentName The experiment name used by the respective folder
#' @param TheN Selects the number of signature variants per peak detector
#' @param Display Default "selection" returns visual plots showing only TheN,
#'  alternatively "all" will show all signatures before filtering in the plots
#' @param AFOverlap Luciernaga_QC default
#' @param subsets The subset gate from which Luciernaga_QC should sample
#' 
#' @importFrom fs file_temp dir_ls file_copy dir_delete
#' @importFrom flowWorkspace keyword
#' @importFrom purrr map
#' @importFrom dplyr mutate group_by ungroup slice_head pull arrange desc
#' @importFrom rlang .data
#' 
#' @return A list containing Data and Plots, for subsequent splitting by 
#' the parent function
#' 
#' @noRd
LuciernagaLocal <- function(x,
                             outpath,
                             TheN,
                             Display,
                             ExperimentName,
                             AFOverlap,
                             subsets) {
  LuciernagaTemp <- file_temp("Luciernaga_Temp_")
  dir.create(LuciernagaTemp)
  # dir.exists(LuciernagaTemp)

  first <- keyword(x, "GROUPNAME")
  second <- keyword(x, "TUBENAME")
  LocalPlotName <- paste(first, second, sep = "_")

  ReturnedFCS <- map(.x = x, .f = Luciernaga_QC,
                      subsets = subsets,
                      removestrings = ".fcs",
                      sample.name = c("GROUPNAME", "TUBENAME"),
                      unmixingcontroltype = "cells",
                      Unstained = TRUE,
                      ratiopopcutoff = 0.01,
                      Verbose = FALSE,
                      AFOverlap = AFOverlap,
                      stats = "median",
                      ExportType = "fcs",
                      SignatureReturnNow = FALSE,
                      outpath = LuciernagaTemp,
                      Increments = 0.1,
                      experiment = ExperimentName,
                      condition = "NA",
                      minimalfcscutoff = 0.001,
                      NegativeType = "artificial")

  TheTempFiles <- dir_ls(LuciernagaTemp, glob = "*.fcs")

  ReturnedFCS <- ReturnedFCS[[1]]

  DecisionData <- ReturnedFCS |>
    mutate(MainDetector = gsub("_.*", "", .data[["Cluster"]])) |>
    group_by(MainDetector) |>
    slice_head(n = TheN) |>
    ungroup()

  Clusters <- sub("_.*", "", DecisionData$Cluster)
  ClustersData <- data.frame(table(Clusters)) |>
    arrange(desc(Freq)) |>
    pull(Clusters)

  TheseSpecimens <- DecisionData |>
    mutate(Cluster = gsub("_", "", Cluster)) |>
    mutate(Cluster = gsub("-", "", Cluster)) |>
    mutate(ID = paste(Sample, Cluster, sep = "_")) |>
    pull(ID)

  TheseSpecimens <- paste0(TheseSpecimens, ".fcs")

  CopyThese <- TheTempFiles[basename(TheTempFiles) %in% TheseSpecimens]

  file_copy(CopyThese, outpath, overwrite = TRUE)

  if (Display == "selection") {
    DataToUse <- DecisionData
    TheLegend <- TRUE
  } else {
    DataToUse <- ReturnedFCS
    TheLegend <- FALSE
  }

  Plots <- map(.x = ClustersData, .f = SignatureVariants, data = DataToUse,
               returnType = "plots", legend = TheLegend,
               plotname = LocalPlotName)

  dir_delete(LuciernagaTemp)
  # dir.exists(LuciernagaTemp)

  ReturnObjects <- list(Data = ReturnedFCS, Plots = Plots)
  return(ReturnObjects)
}