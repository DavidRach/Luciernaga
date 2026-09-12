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
    
    plots <- map(.x=TheseDates, .f=QCPlot_Amalgamate, data=Filtered,
    normalize=normalize)

    if (returnType == "pdf"){
        Utility_Patchwork(x=plots, filename=filename, outfolder=outpath,
        therows=therows, thecolumns=thecolumns, returntype="pdf")
    } else {
        return(plots)
    }
}



