#' Internal for LuciernagaQC, creates .fcs files
#'
#' @param x The data.frame of Luciernaga data.
#' @param ff An individual cytoset object.
#' @param minimalfcscutoff A ratio indicating mininum of the total
#'  population needed  to split off into own file, default is set to 0.05
#' @param AggregateName Passed final name with modifications from name
#' @param Brightness Whether to additionally return a brightness .csv to
#'  the outpath
#' @param outpath Location to export the fcs and .csv files to
#' @param OriginalStart Passed Argument indicating start column for
#'  Raw .fcs values
#' @param OrigingalEnd Passed argument indicating end column for 
#' raw .fcs values
#' @param stats Whether "median" or "mean", default is "median"
#' @param NegativeType Whether to append a negative pop. Args are
#'  "artificial", "internal" and "default"
#' @param TotalNegatives How many of the above rows to append, default
#'  is set to 500
#' @param Samples When Negative type = "Internal", the data.frame of averaged
#' fluorescence per detector
#' @param ExportType Passed from above, set to "fcs" for fcs.file return
#'
#' @importFrom flowCore parameters
#' @importFrom flowWorkspace keyword
#' @importFrom dplyr arrange filter pull bind_rows
#' @importFrom purrr map
#' @importFrom utils write.csv
#'
#' @return An internal value
#'
#' @noRd
Genesis <- function(x, ff, minimalfcscutoff, AggregateName,
  Brightness, outpath=NULL, OriginalStart, OriginalEnd,
  stats = "median", NegativeType="default", TotalNegatives=500,
  Samples=NULL, ExportType, Consolidate){

  # Replicate the Original FCS Parameters
  FlowFrameTest <- ff[[1, returnType = "flowFrame"]]
  original_p <- parameters(FlowFrameTest)
  original_d <- keyword(FlowFrameTest)

  if(!is.null(Consolidate)){

    ConsolidatePaa <- function(x, data){
      Testing <- data |> dplyr::filter(str_detect(Cluster, x))
      x <- gsub("^", "", fixed=TRUE, x)
      x <- gsub("|", "and", fixed=TRUE, x)
      Testing$Cluster <- x
      Testing$Cluster <- factor(Testing$Cluster)
      return(Testing)
    }

    if (length(Consolidate) > 1){
      data <- x
       TheConsolidated <- map(
        .x=Consolidate, data=data, .f=ConsolidatePaa) |>
         bind_rows()
       x <- TheConsolidated
    } else {
      data <- x
      TheConsolidated <- ConsolidatePaa(x=Consolidate, data=data)
      x <- TheConsolidated 
    }
    x$Cluster <- factor(x$Cluster)
  } else {
    x$Cluster <- factor(x$Cluster)
  }

  # Figure out what clusters to split from the file.

  ZZZ <- data.frame(table(x$Cluster))
  ZZZ <- ZZZ |> arrange(desc(Freq))
  colnames(ZZZ)[1] <- "Cluster"
  colnames(ZZZ)[2] <- "Count"
  fcs_cutoff <- nrow(x)*minimalfcscutoff
  fcs_clusters <- ZZZ |> filter(Count > fcs_cutoff) |> pull(Cluster)

  Data <- x

  TheBrightness <- map(.x=fcs_clusters, .f=Luciernaga:::InternalGenesis,
     Data=Data,
    AggregateName=AggregateName, outpath=outpath, OriginalStart=OriginalStart,
    OriginalEnd=OriginalEnd, stats=stats, NegativeType=NegativeType,
    TotalNegatives=TotalNegatives, Samples=Samples, ExportType=ExportType,
    parameters=original_p, description=original_d) |> bind_rows()

  #message("TargetReached")

  if (Brightness == TRUE){
    RelativeBrightness <- RelativeBrightness(TheBrightness)
    CSVName <- paste0("RelativeBrightness", AggregateName, ".csv")
    CSVSpot <- file.path(outpath, CSVName)
    write.csv(RelativeBrightness, CSVSpot, row.names = FALSE)
  }
}