#' Used to retrieve effective net signatures
#' 
#' @param data A LuciernagaQC dataset for all specimens
#' @param groupThese Default is c("Sample", "Condition")
#' @param stats Default is "median"
#' @param normalize Default FALSE
#' @param returnType Default is "data", alternate is "plot"
#' @param legend Default FALSE, setting TRUE will display legend on right
#' @param linecolor Default is red for the averaged signature line
#' @param name Default NULL, sets plot title
#' 
#' @importFrom dplyr group_by mutate filter ungroup arrange desc
#'  pull relocate slice select bind_rows bind_cols
#' @importFrom tidyselect where all_of
#' @importFrom tidyr pivot_longer uncount
#' @importFrom rlang syms
#' @importFrom ggplot2 ggplot aes geom_line scale_color_manual
#'  labs theme_linedraw theme_bw theme element_text element_blank
#' @importFrom stats setNames median
#' 
#' @return A ggplot object
#' 
#' @noRd
NetSignatures <- function(data, groupThese=c("Sample", 
  "Condition"), stats="median", normalize=FALSE, returnType="data",
  legend=FALSE, linecolor="red", name=NULL){

data <- data |> ungroup()
Internal <- data |> select(all_of(groupThese), where(is.numeric))
WorkingData <- Internal %>%
  tidyr::uncount(weights=Count, .remove=TRUE, .id=NULL)
WorkingData <- WorkingData |>
  mutate(TheSample = paste(!!!syms(groupThese), sep = "_")) |>
  relocate(TheSample, .before=1) |> select(-all_of(groupThese))

NetAverage <- WorkingData |> group_by(TheSample) |> 
  AveragedSignature(stats=stats)
  
Averaged <- WorkingData |> AveragedSignature(stats=stats)
Averaged[1,1] <- "Average"
  
Dataset <- bind_rows(NetAverage, Averaged)
  
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
  
if (returnType == "data"){return(Dataset)
} else {

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

}