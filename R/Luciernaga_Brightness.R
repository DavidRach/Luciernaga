#' Visualize MFI of raw .fcs files to evaluate single color controls.
#'
#' @param x A passed Fluorophore name corresponding to those found in the data.
#' @param data A data.frame containing the Luciernaga_FCSToReport raw MFI values.
#' @param Downsample Whether to Downsample, default is set to TRUE to not overly
#' affect the y-axis
#' @param subsample If desired, a certain number of cells to take from each cluster.
#' Default is NULL, selecting the number of cells found in the smallest cluster.
#' @param reference A .csv path or data.frame containing Fluorophore and Detector info.
#' @param clearance A buffer area multiplier to xmin and xmax when Scaled equals False
#' @param Scaled Default is set to TRUE, returning logicle transformed MFI. FALSE
#' returns raw MFI.
#'
#' @importFrom dplyr filter
#' @importFrom dplyr arrange
#' @importFrom dplyr pull
#' @importFrom dplyr group_by
#' @importFrom dplyr slice_sample
#' @importFrom dplyr ungroup
#' @importFrom dplyr select
#' @importFrom ggplot2 ggplot
#' @importFrom ggcyto scale_x_logicle
#' @importFrom utils head
#' @importFrom utils read.csv
#'
#' @return ggplot objects for each fluorophore containing the various clusters
#' @export
#'
#' @examples NULL

Luciernaga_Brightness <- function(x, data, Downsample=TRUE, subsample = NULL,
                                  reference, clearance=0.02, Scaled = TRUE){
    TheFluorophore <- x


    TheData <- data %>% filter(Sample %in% x)
    TheTable <- data.frame(table(TheData$Cluster), check.names = FALSE)
    colnames(TheTable)[1] <- "Cluster"
    colnames(TheTable)[2] <- "Count"
    TheSlice <- TheTable %>% arrange(Count) %>% head(1) %>% pull(Count)

    if (Downsample == TRUE) {
      if (is.null(subsample)){
       TheData <- TheData %>% group_by(Cluster) %>%
         slice_sample(n=TheSlice, replace = FALSE) %>% ungroup()
      } else {TheData <- TheData %>% group_by(Cluster) %>%
        slice_sample(n=subsample, replace = FALSE) %>% ungroup()
      }
    }

    if (!is.data.frame(reference)){CSV <- read.csv(reference, check.names = FALSE)
      } else {CSV <- reference}

    internalstrings <- c(" ", ".", "_", "-A")
    CSV$Fluorophore <- NameCleanUp(CSV$Fluorophore, removestrings=internalstrings)
    CSV$Detector <- NameCleanUp(CSV$Detector, removestrings=internalstrings)
    TheDetector <- CSV %>% dplyr::filter(Fluorophore %in% TheFluorophore) %>%
      pull(Detector)

    Values <- TheData %>% select(all_of(TheDetector)) %>% as.matrix()
    theXmin <- Values %>% quantile(., 0.01)
    theXmax <- Values %>% quantile(., 0.99)
    theXmin <- theXmin - abs((clearance*theXmin))
    theXmax <- theXmax + (clearance*theXmax)

    if (Scaled == TRUE){
      plot <- ggplot(TheData, aes(x=.data[[TheDetector]], fill=Cluster)) +
        geom_density(alpha=0.5) + scale_x_logicle(w=1.5, t=4200000, m=5.62) +
        labs(title=TheFluorophore, x=TheDetector, y="Frequency") + theme_bw() +
        theme(axis.title.x=element_text(face="plain"), axis.title.y=element_text(
            face="plain"))
      } else {
        plot <- ggplot(TheData, aes(x=.data[[TheDetector]], fill=Cluster)) +
          geom_density(alpha=0.5) + xlim(theXmin, theXmax) +
          labs(title=TheFluorophore, x=TheDetector, y="Frequency") + theme_bw() +
          theme(axis.title.x=element_text(face="plain"), axis.title.y=element_text(
              face="plain"))
    }

 return(plot)
}

