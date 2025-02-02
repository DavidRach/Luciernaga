#' Generates cosine comparison from a data.frame of fluorescent signatures
#'
#' @param data A data.frame with a single name column and rest numeric columns
#' @param returntype Default returns "plot", alternatively "matrix" for underlying data
#'
#' @importFrom dplyr select
#' @importFrom tidyselect where
#' @importFrom lsa cosine
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 scale_fill_gradient2
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 coord_fixed
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_text
#'
#' @return Either a ggplot or matrix object
#' @export
#'
#' @examples
#'
#' library(flowCore)
#' library(flowWorkspace)
#' library(openCyto)
#' library(data.table)
#' library(dplyr)
#'
#' File_Location <- system.file("extdata", package = "Luciernaga")
#' FCS_Files <- list.files(path = File_Location, pattern = ".fcs",
#'   full.names = TRUE)
#' UnstainedFCSFiles <- FCS_Files[grep("Unstained", FCS_Files)]
#' UnstainedCells <- UnstainedFCSFiles[-grep("Beads", UnstainedFCSFiles)]
#' MyCytoSet <- load_cytoset_from_fcs(UnstainedCells[1],
#'   truncate_max_range = FALSE,transformation = FALSE)
#' MyGatingSet <- GatingSet(MyCytoSet)
#' MyGates <- fread(file.path(path = File_Location, pattern = 'Gates.csv'))
#' MyGatingTemplate <- gatingTemplate(MyGates)
#' gt_gating(MyGatingTemplate, MyGatingSet)
#' removestrings <-  c("DR_", "Cells", ".fcs", "-", " ")
#' StorageLocation <- file.path("C:", "Users", "JohnDoe", "Desktop")
#'
#' PopulationInterest <- gs_pop_get_data(MyGatingSet[1], subset="lymphocytes")
#' TheDataValues <- exprs(PopulationInterest[[1]])
#' TheDataValues <- data.frame(TheDataValues, check.names=FALSE)
#' Signature <- AveragedSignature(TheDataValues, stats="median")
#' TheData1 <- Signature[,-grep("Time|FS|SC|SS|Original|W$|H$", names(Signature))]
#' TheData1 <- TheData1 %>% mutate(Sample="lymphocytes") %>%
#'  relocate(Sample, .before=1)
#'
#' PopulationInterest <- gs_pop_get_data(MyGatingSet[1], subset="nonDebris")
#' TheDataValues <- exprs(PopulationInterest[[1]])
#' TheDataValues <- data.frame(TheDataValues, check.names=FALSE)
#' Signature <- AveragedSignature(TheDataValues, stats="median")
#' TheData2 <- Signature[,-grep("Time|FS|SC|SS|Original|W$|H$", names(Signature))]
#' TheData2 <- TheData2 %>% mutate(Sample="nonDebris") %>%
#'  relocate(Sample, .before=1)
#'
#' FinalData <- rbind(TheData1, TheData2)
#'
#' Plot <- Luciernaga_Cosine(data=FinalData, returntype="plot")
#'
Luciernaga_Cosine <- function(data, returntype="plot", rearrange=TRUE){

    Names <- data %>% select(!where(is.numeric))
    if (ncol(Names) > 1){stop("Please use single column for names")}
    Names <- Names[[1]]

    Numbers <- data %>% select(where(is.numeric))
    NumericsT <- t(Numbers)
    rownames(NumericsT) <- NULL
    colnames(NumericsT) <- Names
    NumericsT <- data.matrix(NumericsT)

    data <- as.matrix(data)
    CosineMatrix <- cosine(NumericsT)
    CosineMatrix <- round(CosineMatrix, 2)
  
    if (rearrange==TRUE){
    Reordered <- ReorderedCosine(CosineMatrix)
    } else {Reordered <- CosineMatrix}
    MeltedCosine <- melt(Reordered)

    #Generate a Red to Blue Heatmap
    CosinePlot <- ggplot(MeltedCosine, aes(Var2, Var1, fill = value)) +
      geom_tile(color = "white") +
      scale_fill_gradient2(low = "lightblue", high = "orange", mid = "white",
                           midpoint = 0.7, limit = c(0.4,1), space = "Lab",
                           name="Cosine\nSimilarity") +
      theme_bw() + geom_text(aes(Var2, Var1, label = value), color = "black",
                             size = 2) + coord_fixed(ratio = 1.3) +
      theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
            panel.grid.major = element_blank(), panel.border = element_blank(),
            panel.background = element_blank(), axis.ticks = element_blank(),
            legend.position.inside = c(1.2, 0.5),
            legend.direction = "vertical", axis.text.x = element_text(
              angle = 45, vjust = 1, hjust = 1, size = 6),
            axis.text.y = element_text(size = 6),
            legend.key.size = unit(0.4, "cm"))

    if (returntype == "plot"){
      return(CosinePlot)
    } else {return(Reordered)}
}
