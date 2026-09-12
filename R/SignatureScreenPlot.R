#' Takes LuciernagaQC data output for dataset, returns amalgamated
#' 
#' @param data A LuciernagaQC dataset for all specimens
#' @param RetainThese A message of starting values is displayed, providing an
#' matching argument character string will filter for those to display
#' @param normalize Default FALSE
#' @param legend Default FALSE, setting TRUE will display legend on right
#' @param linecolor Default is red for the averaged signature line
#' @param name Default NULL, sets plot title
#' 
#' @importFrom dplyr group_by mutate filter ungroup arrange desc
#'  pull relocate slice select bind_rows bind_cols
#' @importFrom tidyselect where all_of
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_line scale_color_manual
#'  labs theme_linedraw theme_bw theme element_text element_blank
#' @importFrom stats setNames median
#' 
#' @return A ggplot object
#' 
#' @noRd
SignatureScreenPlot <- function(data, RetainThese=NULL, normalize=FALSE,
  legend=FALSE, linecolor="red", name=NULL){
WorkingData <- data |> group_by(Sample, Cluster) |>
   mutate(Ratio = Count / sum(Count)) |> filter(Ratio > 0.01) |>
   ungroup()

WorkingData$Cluster <- as.character(WorkingData$Cluster)

Filtering <- WorkingData |>
   mutate(ClusterStart = sub("_.*", "", Cluster))

Table <- data.frame(table(Filtering$ClusterStart))
Table <- Table |> arrange(desc(Freq))
These <- Table |> pull(Var1) |> paste(collapse = ", ")
message("Values are ", These)

MainValue <- Table |> slice(1) |> pull(Var1)
Values <- Filtering |> filter(ClusterStart %in% MainValue) |>
   select(-Count, -Ratio) |> select(where(is.numeric))
Average <- AveragedSignature(x=Values, stats=median)
Average <- Average |> mutate(TheSample = "Average") |>
   relocate(TheSample, .before=1)

if (!is.null(RetainThese)){
  Filtering1 <- Filtering |> filter(ClusterStart %in% RetainThese)
} else {Filtering1 <- Filtering}

Filtering1 <- Filtering1 |>
  mutate(TheSample=paste(Sample, Cluster, sep="_")) |>
  select(-Count, -Ratio) |> select(TheSample, where(is.numeric))

Dataset <- bind_rows(Filtering1, Average)

Numerics <- sapply(Dataset, is.numeric)
Dataset[Numerics] <- lapply(Dataset[Numerics], round, digits = 3)
   
if (normalize == TRUE){
  if (any(Dataset |> select(where(is.numeric)) > 1)){
      Metadata <- Dataset |> select(!where(is.numeric))
      Numerics <- Dataset |> select(where(is.numeric))
      n <- Numerics
      n[n < 0] <- 0
      A <- do.call(pmax, n)
      Normalized <- n/A
      Dataset <- bind_cols(Metadata, Normalized)
      }
  }
   
colnames(Dataset) <- gsub("Comp-", "", colnames(Dataset))
colnames(Dataset) <- gsub("-A", "", colnames(Dataset))
LineCols <- ncol(Dataset)
DetectorOrder <- colnames(Dataset)
DetectorOrder <- DetectorOrder[-1]
Melted <- Dataset |> pivot_longer(all_of(2:LineCols),
  names_to = "Detector", values_to = "value")

Melted$Detector <- factor(Melted$Detector, levels = DetectorOrder)

if (legend == TRUE){LegendPosition <- "right"
  } else {LegendPosition <- "none"}

Melted[["TheSample"]] <- factor(Melted[["TheSample"]],
    levels = c(setdiff(
      unique(Melted[["TheSample"]]), "Average"),"Average")
  )

if (normalize == TRUE){Expression <- "Normalized MFI"
  } else {Expression <- "MFI"}
   

plot <- ggplot(Melted, aes(x = Detector, y = value, group = TheSample,
             color = TheSample)) + geom_line(alpha=0.5, size=0.2) +
             scale_color_manual(values = setNames(
                  ifelse(unique(Melted[["TheSample"]]) == "Average", linecolor, "gray"),
                  unique(Melted[["TheSample"]]))) +
             labs(title = name, x = "Detectors", y = Expression) +
             theme_linedraw() + theme_bw() + theme(axis.title.x = element_text(
             face = "plain"), axis.title.y = element_text(face = "plain"),
             axis.text.x = element_text(size = 5,
             angle = 45, hjust = 1), panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(), legend.position=LegendPosition)

return(plot)

}




