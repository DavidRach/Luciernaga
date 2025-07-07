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
#' pull relocate slice select bind_rows bind_cols
#' @importFrom tidyselect where all_of
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_line scale_color_manual
#' labs theme_linedraw theme_bw theme element_text element_blank
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
#' pull relocate slice select bind_rows bind_cols
#' @importFrom tidyselect where all_of
#' @importFrom tidyr pivot_longer uncount
#' @importFrom rlang syms
#' @importFrom ggplot2 ggplot aes geom_line scale_color_manual
#' labs theme_linedraw theme_bw theme element_text element_blank
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

#' Takes the normalize false data output from NetSignatures and
#' returns the shift matrix
#' 
#' @param data The data output from NetSignatures
#' @param sampleName Column name of the sample, default is TheSample
#' 
#' @importFrom dplyr filter select across cur_column bind_cols mutate
#' @importFrom tidyselect where all_of everything 
#' 
#' @return A drift data.frame for use in simulating fluorophores
#' 
#' @noRd
SignatureShifts <- function(data, sampleName="TheSample"){
  TheNames <- colnames(data)
  TheNames <- TheNames[-1]
  AverageReference <- data |> filter(.data[[sampleName]] %in% "Average") |>
    select(where(is.numeric)) |> unlist()

  OtherReferences <- data |> filter(!.data[[sampleName]] %in% "Average")
  OtherMetadata <- OtherReferences |> select(all_of(sampleName))
  OtherNumeric <- OtherReferences |> select(!all_of(sampleName))

  SubtractedData <- OtherNumeric %>% mutate(across(
    .cols = everything(), .fns  = ~ .x - AverageReference[cur_column()]))

  DriftData <- SubtractedData %>% mutate(across(
      everything(),~ .x / AverageReference[cur_column()]))
  
  OnesData <- data.frame(matrix(1, 
    nrow = nrow(DriftData), 
    ncol = ncol(DriftData),
    dimnames = list(rownames(DriftData), colnames(DriftData))))
  
  Result <- OnesData + DriftData

  colnames(Result) <- TheNames
  
  DriftData <- bind_cols(OtherMetadata, Result)
  
  return(DriftData)
}

#' Simulation function to approximate fluorophore drift
#' 
#' @param Residual The output of SignatureShifts
#' @param NumberDetectors Number of detectors corresponding to your cytometer
#' @param TheFluoruophore The fluorophore of interest
#' @param RestingMFI Value by which reference signature gets multiplied by for
#' this simulation
#' @param legend Default FALSE, TRUE sets on right side plot
#' 
#' @importFrom dplyr filter mutate select left_join filter pull
#' @importFrom tidyselect all_of
#' @importFrom tidyr pivot_longer
#' @importFrom stringr str_detect
#' 
#' @return A signature plot of the reference control vs all. 
#' 
#' @noRd
DriftedFluors <- function(Residual, NumberDetectors=64, TheFluorophore,
  RestingMFI=100000, legend=FALSE){
 
 References <- Luciernaga:::InstrumentReferences(NumberDetectors=NumberDetectors)
 Internal <- References |> filter(Fluorophore %in% TheFluorophore) |>
      mutate(AdjustedY=AdjustedY*RestingMFI)
 
 FinalCol <- ncol(Residual)
 
 Internal <- Internal |> select(-Instrument, -Fluorophore)
 Residual <- Residual |> tidyr::pivot_longer(all_of(2:FinalCol),
  names_to="Detector", values_to="Adjustment")
 Residual$Detector <- gsub("-A", "", Residual$Detector)
 Merge <- left_join(Residual, Internal, by="Detector")
 Merge <- Merge |> mutate(AdjustedMFI=Adjustment*AdjustedY)
 Merge <- Merge |> filter(!str_detect(Detector, "SC"))
  
 Merge <- Merge |>
      mutate(Signature = AdjustedMFI/max(AdjustedMFI, na.rm = TRUE)) |>
     ungroup()
 
 Merge2 <- Merge |> select(-Adjustment, -AdjustedY, -AdjustedMFI)
  
 Merge2$TheSample <- sub(".*ter_", "", Merge2$TheSample )
 Merge2$TheSample <- sub(".*ore_", "", Merge2$TheSample )
 
 TheseDates <- Merge2 |> pull(TheSample) |> unique()
 
 Plot <- QC_ViewSignature(x=TheseDates, columnname="TheSample", data=Merge2,
  TheFormat="longer", detectorcolumn = "Detector", valuecolumn = "Signature",
  Normalize=FALSE, legend=legend)

 #plotly::ggplotly(Plot)
 
 return(Plot)
}
