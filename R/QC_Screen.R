#' Takes Luciernaga_QC data output, returns amalgamate plots, used
#' to screen QC beads for major issues that notably shift the average
#' 
#' @param data The LuciernagaQC data output
#' @param groupColumns Names of columns you want to combine to make the 
#' name being used for comparison, ex: c("Sample", "Condition")
#' @param clusterColumn Default is "Cluster"
#' @param filterDate Used to quickly return a screened plot, use the plot title
#' name in it's entirety and set returnType to "plot"
#' @param cutoff Default is 0.005
#' @param normalize Default is FALSE
#' @param returnType Default is pdf, alternatively plot
#' @param outpath The desired storage folder
#' @param filename The desired filename for the output to be saved as
#' @param therows For pdf output, number of rows per page, default 3
#' @param thecolumns For pdf output, number of columns per page, default 3
#' 
#' @importFrom tidyr unite
#' @importFrom dplyr mutate group_by filter ungroup pull
#' @importFrom tidyselect all_of 
#' @importFrom purrr map
#' 
#' @return Either individual plots or a pdf desired location
#' 
#' @export 
QC_Screen <- function(data, groupColumns, clusterColumn="Cluster",
 filterDate=NULL, cutoff=0.005, normalize=FALSE, returnType="pdf", 
 outpath, filename="BeadQCVisualized", therows=3, thecolumns=3){

    data[[clusterColumn]] <- as.character(data[[clusterColumn]])
    data <- data |> unite(TheSample,
     all_of(groupColumns), sep="_", remove=FALSE)
    Internal <- data |> mutate(
        Split_Cluster = strsplit(.data[[clusterColumn]], "-"))
    Filtered <- Internal |> group_by(TheSample) |>
         mutate(Ratio=Count/sum(Count)) |> filter(Ratio > cutoff) |>
         ungroup()

    TheseDates <- Filtered |> pull(TheSample) |> unique()
    # x <- TheseDates[1] 
    if (!is.null(filterDate)){
        TheseDates <- filterDate
    }
    
    plots <- map(.x=TheseDates, .f=QCPlot_Amalgamate, data=Filtered)

    if (returnType == "pdf"){
        Utility_Patchwork(x=plots, filename=filename, outfolder=outpath,
        therows=therows, thecolumns=thecolumns, returntype="pdf")
    } else {
        return(plots)
    }
}


#' Internal for QC_Screen, produced the amalgamate plots for the passed inputs
#' 
#' @param x The identifying name for which data will be filtered for
#' @param data The passed data 
#' @param countcolumn Default is Count
#' @param samplecolumn Default is TheSample
#' @param linecolor Default is "red", used to show Averaged Signature
#' @param returnType Default is "plot"
#' @param titlename Default is NULL, alternatively provided sets title
#' @param legend Default is FALSE, when TRUE sets legend.position to "right"
#' 
#' @importFrom dplyr filter pull select summarize across relocate bind_rows
#' bind_cols 
#' @importFrom tidyr uncount pivot_longer
#' @importFrom tidyselect where everything
#' @importFrom rlang sym
#' @importFrom stats median setNames
#' @importFrom ggplot2 ggplot aes geom_line scale_color_manual labs theme_linedraw
#' theme_bw theme element_text element_blank
#' 
#' @return An individual ggplot object
#' 
#' @noRd
QCPlot_Amalgamate <- function(x, data, countcolumn="Count",
 samplecolumn="TheSample", linecolor="red", returnType="plot",
 titlename=NULL, legend=FALSE){

    SmallPortion <- data |> filter(TheSample %in% x)
    if (is.null(titlename)){
        TheTitle <- SmallPortion |> pull(TheSample)}
    NumericPortion <- SmallPortion |>
        tidyr::uncount(weights = !!sym(countcolumn), .remove = TRUE) |>
        select(where(is.numeric)) |> select(-Ratio)
    AverageCalced <- NumericPortion |> summarize(across(everything(), median))  |>
     mutate(Cluster = "Average") |> relocate(Cluster, .before = 1)
    Rest <- SmallPortion |> select(Cluster, where(is.numeric)) |>
        tidyr::uncount(weights = !!sym(countcolumn), .remove = TRUE) |>
        select(-Ratio)
    ComparisonTable <- bind_rows(Rest, AverageCalced)

    Numerics <- sapply(ComparisonTable, is.numeric)
    ComparisonTable[Numerics] <- lapply(ComparisonTable[Numerics], round, digits = 3)
     
    if (normalize == TRUE){
        if (any(ComparisonTable |> select(where(is.numeric)) > 1)){
        Metadata <- ComparisonTable |> select(!where(is.numeric))
        Numerics <- ComparisonTable |> select(where(is.numeric))
        n <- Numerics
        n[n < 0] <- 0
        A <- do.call(pmax, n)
        Normalized <- n/A
        ComparisonTable <- bind_cols(Metadata, Normalized)
        }
    }
     
    colnames(ComparisonTable) <- gsub("Comp-", "", colnames(ComparisonTable))

    LineCols <- ncol(ComparisonTable)
    DetectorOrder <- colnames(ComparisonTable)
    DetectorOrder <- DetectorOrder[-1]
     
    Melted <- ComparisonTable |> 
          pivot_longer(all_of(2:LineCols),
          names_to = "Detector", values_to = "value")

    Melted$Detector <- factor(Melted$Detector, levels = DetectorOrder)
    
    if (!is.null(titlename)){name <- titlename
    } else {name <- TheTitle}

    if (legend == TRUE){LegendPosition <- "right"
    } else {LegendPosition <- "none"}

    Melted[["Cluster"]] <- factor(Melted[["Cluster"]],
          levels = c("Average", setdiff(unique(Melted[["Cluster"]]), "Average")))

    if (normalize == TRUE){Expression <- "Normalized MFI"
     } else {Expression <- "MFI"}
     

    plot <- ggplot(Melted, aes(x = Detector, y = value, group = Cluster,
               color = Cluster)) + geom_line(alpha=0.5, linewidth=0.5) +
               scale_color_manual(values = setNames(
                    ifelse(unique(Melted[["Cluster"]]) == "Average", linecolor, "gray"),
                    unique(Melted[["Cluster"]]))) +
               labs(title = name, x = "Detectors", y = Expression) +
               theme_linedraw() + theme_bw() + theme(axis.title.x = element_text(
               face = "plain"), axis.title.y = element_text(face = "plain"),
               axis.text.x = element_text(size = 5,
               angle = 45, hjust = 1), panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(), legend.position=LegendPosition)

    return(plot)
}

