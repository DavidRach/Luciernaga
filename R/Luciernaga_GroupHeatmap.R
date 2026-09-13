#' A group version of the heatmap option from Luciernaga_Plots
#'
#' @param reports A data.frame or a list of Luciernaga_QC report
#'  data objects
#' @param nameColumn The name of the column that differentiates
#'  between the reports
#' @param cutoff Proportion of cells that at least 1 report needs
#'  to exceed for retention.
#' @param returntype Either "plot" or underlying "data"
#' @param legend Default is "right", use "none" to remove
#' @param transpose Default is FALSE, flips orientation
#'
#' @importFrom purrr map
#' @importFrom dplyr bind_rows group_by summarize left_join mutate
#'  relocate select pull rename filter
#' @importFrom rlang .data
#'
#' @return Either a plot or underlying data
#' @export
#'
#' @examples
#'
#' library(flowCore)
#' library(flowWorkspace)
#' library(openCyto)
#' library(data.table)
#' library(dplyr)
#' library(purrr)
#' File_Location <- system.file("extdata", package = "Luciernaga")
#' FCS_Files <- list.files(path = File_Location, pattern = ".fcs", full.names = TRUE)
#' UnstainedFCSFiles <- FCS_Files[grep("Unstained", FCS_Files)]
#' UnstainedCells <- UnstainedFCSFiles[-grep(
#'   "Beads", UnstainedFCSFiles)]
#' MyCytoSet <- load_cytoset_from_fcs(UnstainedCells[c(1,3,5)],
#'                                    truncate_max_range = FALSE,transformation = FALSE)
#' MyGatingSet <- GatingSet(MyCytoSet)
#' MyGates <- fread(file.path(path = File_Location, pattern = 'Gates.csv'))
#' MyGatingTemplate <- gatingTemplate(MyGates)
#' gt_gating(MyGatingTemplate, MyGatingSet)
#' removestrings <-  c("DR_", "Cells", ".fcs", "-", " ")
#' StorageLocation <- file.path("C:", "Users", "JohnDoe", "Desktop")
#' FileLocation <- system.file("extdata", package = "Luciernaga")
#' pattern = "AutofluorescentOverlaps.csv"
#' AFOverlap <- list.files(path=FileLocation, pattern=pattern, full.names = TRUE)
#'
#' reports <- map(.x=MyGatingSet[1:3], .f=Luciernaga_QC, subsets="lymphocytes",
#'                removestrings=removestrings, sample.name="GROUPNAME",
#'                unmixingcontroltype = "cells", Unstained = TRUE,
#'                ratiopopcutoff = 0.001, Verbose = FALSE, AFOverlap = AFOverlap,
#'                stats = "median", ExportType = "data", SignatureReturnNow = FALSE,
#'                outpath = TemporaryFolder, Increments=0.1, SecondaryPeaks=2,
#'                experiment = "Lymphocytes", condition = "Ctrl",
#'                Subtraction = "Internal", SCData="subtracted",
#'                NegativeType="default")
#'
#' plot <- Luciernaga_GroupHeatmap(reports=reports, nameColumn="Sample",
#'  cutoff=0.02, returntype = "plot")
#'
Luciernaga_GroupHeatmap <- function(reports, nameColumn, cutoff = 0.01,
                                     returntype = "plot", legend = "right",
                                     transpose = FALSE) {
  #nameColumn <- "Experiment"
  Columns <- c(nameColumn, "Cluster", "Count")

  if (!is.data.frame(reports)) {
    Processed <- map(.x = reports, .f = Luciernaga:::ReportProcess,
      columns = Columns) |> bind_rows()
  } else {
    Processed <- Luciernaga:::ReportProcess(x = reports, columns = Columns)
  }
  #nrow(Processed)

  Processed <- Processed |> unique()

  TheCounts <- Processed |> group_by(.data[[nameColumn]]) |>
    summarize(TotalCells = sum(Count, na.rm = TRUE), .groups = "drop")

  Processed <- Processed |> left_join(TheCounts, by = c(
    nameColumn))

  TheData <- Processed |> mutate(Ratio = round(Count / TotalCells, 3)) |>
    relocate(Ratio, .after = Count) |> select(-TotalCells)

  TheClusters <- TheData |> pull(Cluster) |> unique()

  Values <- TheData |> group_by(.data[[nameColumn]], Cluster) |>
    mutate(cutoff = Ratio > cutoff) #Set as cutoff value

  Clusters <- map(.x = TheClusters, .f = Luciernaga:::ClusterAbundance,
    data = Values)
  Clusters <- Filter(Negate(is.null), Clusters)
  Clusters <- unlist(Clusters)
  ExcludedClusters <- setdiff(TheClusters, Clusters)

  FilteredData <- TheData |> filter(Cluster %in% Clusters)

  OtherData <- FilteredData |> group_by(.data[[nameColumn]]) |>
    summarize(LostRatio = 1 - sum(Ratio, na.rm = TRUE), .groups = "drop")

  Other <- TheCounts |> left_join(OtherData, by = nameColumn) |>
    mutate(Count = round(TotalCells * LostRatio, 0)) |> select(-TotalCells) |>
    relocate(Count, .before = LostRatio) |> rename(Ratio = LostRatio) |>
    mutate(Cluster = "Other") |> relocate(Cluster, .before = "Count")

  UpdatedDataset <- bind_rows(FilteredData, Other)

  if (returntype == "plot") {
    plot <- Luciernaga:::StackedReportHeatmap(data = UpdatedDataset,
      nameColumn = nameColumn, legend = legend, transpose = transpose)

    return(plot)
  } else {
    return(UpdatedDataset)
  }
}