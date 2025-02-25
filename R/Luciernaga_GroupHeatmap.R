#' A group version of the heatmap option from Luciernaga_Plots
#'
#' @param reports A data.frame or a list of Luciernaga_QC report data objects
#' @param nameColumn The name of the column that differentiates between the reports
#' @param cutoff Proportion of cells that at least 1 report needs to exceed for retention.
#' @param returntype Either "plot" or underlying "data"
#'
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
#' @importFrom dplyr group_by
#' @importFrom dplyr summarize
#' @importFrom dplyr left_join
#' @importFrom dplyr mutate
#' @importFrom dplyr relocate
#' @importFrom dplyr select
#' @importFrom dplyr pull
#' @importFrom dplyr rename
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
Luciernaga_GroupHeatmap <- function(reports, nameColumn, cutoff=0.01, returntype="plot"){
  #nameColumn <- "Experiment"
  Columns <- c(nameColumn, "Cluster", "Count")

  if (!is.data.frame(reports)){
  Processed <- map(.x=reports, .f=Luciernaga:::ReportProcess,
                   columns=Columns) |> bind_rows()
  } else {Processed <- Luciernaga:::ReportProcess(x=reports, columns=Columns)}
  #nrow(Processed)
  
  Processed <- Processed |> unique()

  TheCounts <- Processed |> group_by(.data[[nameColumn]]) |>
    summarize(TotalCells = sum(Count, na.rm = TRUE), .groups = 'drop')

  Processed <- Processed |> left_join(TheCounts, by = c(
    nameColumn))

  TheData <- Processed |> mutate(Ratio = round(Count/TotalCells, 3)) |>
    relocate(Ratio, .after=Count) |> select(-TotalCells)

  TheClusters <- TheData |> pull(Cluster) |> unique()

  Values <- TheData |> group_by(.data[[nameColumn]], Cluster) |>
    mutate(cutoff = Ratio > cutoff) #Set as cutoff value

  Clusters <- map(.x=TheClusters, .f=Luciernaga:::ClusterAbundance, data=Values)
  Clusters <- Filter(Negate(is.null), Clusters)
  Clusters <- unlist(Clusters)
  ExcludedClusters <- setdiff(TheClusters, Clusters)

  FilteredData <- TheData |> filter(Cluster %in% Clusters)

  OtherData <- FilteredData |> group_by(.data[[nameColumn]]) |>
    summarize(LostRatio = 1 - sum(Ratio, na.rm = TRUE), .groups = 'drop')

  Other <- TheCounts |> left_join(OtherData, by = nameColumn) |>
  mutate(Count = round(TotalCells*LostRatio, 0)) |> select(-TotalCells) |>
  relocate(Count, .before=LostRatio) |> rename(Ratio=LostRatio) |>
    mutate(Cluster="Other") |> relocate(Cluster, .before="Count")

  UpdatedDataset <- bind_rows(FilteredData, Other)

  if (returntype == "plot"){

    plot <- Luciernaga:::StackedReportHeatmap(data=UpdatedDataset, nameColumn=nameColumn)

    return(plot)

  } else {return(UpdatedDataset)}
  }


#' Internal for StackedReport
#'
#' @param data The data intermediate of Stacked Report
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 scale_fill_gradient
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 element_line
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 coord_fixed
#'
#' @noRd
StackedReportHeatmap <- function(data, nameColumn){
  data$Ratio <- round(data$Ratio, 2)

  plot <- ggplot(data, aes(x=.data[[nameColumn]], y = Cluster, fill = Ratio)) +
    geom_tile() + geom_text(aes(label = Ratio)) + theme_bw() +
    scale_fill_gradient(name = "Ratio", low = "#FFFFFF", high = "#FF0000", limits = c(0, NA)) +
    theme(plot.title = element_text(hjust = 0.5), panel.grid.minor = element_line(
      linetype = "blank"), axis.title = element_text(size = 10), axis.title.y = element_blank(),
      axis.title.x = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(),
      axis.text.x = element_text(angle = 40, hjust = 1), legend.key.size = unit(0.4, "cm"))  +
    coord_fixed(ratio = 1.1)

  return(plot)
}

#' Internal for Stacked Reports
#'
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#'
#' @noRd
ClusterAbundance <- function(x, data){
  #x <- TheClusters[1]

  Subset <- data |> dplyr::filter(Cluster %in% x)
  TheValues <- Subset |> pull(cutoff) |> unique()

  if(length(TheValues) == 1 && TheValues== TRUE){Value <- x
  } else if(length(TheValues) == 2 && any(TheValues==TRUE)){Value <- x
  } else {Value <- NULL}

 return(Value)
}

#' Internal for Stacked Report
#' @importFrom dplyr select
#' @importFrom tidyselect all_of
#'
#' @noRd
ReportProcess <- function(x, columns){
  # x <- reports[[1]]
  data <- x %>% dplyr::select(all_of(columns))
}
