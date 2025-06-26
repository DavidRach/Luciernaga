#' Visualize MFI of raw .fcs files to evaluate single color controls.
#'
#' @param fluorophore.name The name of the fluorophore to filter from data
#' @param data The data.frame derrived by Luciernaga_FolderBrightness, consisting
#' of fluorophore and cluster metadata, and the raw detector values
#' @param fluorophore.column The name of the column containing the fluorophore 
#' information
#' @param cluster.column The name of the column containing the cluster information
#' @param downsample Default is TRUE, to avoid having large populations take over
#' the entire y-axis.
#' @param subsample When downsample is true, number of cells from each cluster to
#' keep. The default NULL will select the number of cells found in the smallest cluster
#' @param Detector Default NULL, when reference is NULL sets the detector to plot
#' on the x-axis
#' @param reference A .csv path or data.frame containing Fluorophore
#' and Detector information from which to retrieve the x-axis detector, use same 
#' naming as with the passed fluorophore argument in x
#' @param clearance A buffer area multiplier to xmin and xmax when
#'  scaled equals false
#' @param Scaled Default is set to TRUE, returning log_10 transformed data, 
#' FALSE will return raw MFI.
#' @param maxtik Default is 1e6
#' @param legend Default is "right"
#'
#' @importFrom dplyr filter
#' @importFrom dplyr arrange
#' @importFrom utils head 
#' @importFrom dplyr pull
#' @importFrom dplyr group_by
#' @importFrom dplyr slice
#' @importFrom dplyr slice_sample
#' @importFrom dplyr ungroup
#' @importFrom utils read.csv
#' @importFrom dplyr pull
#' @importFrom dplyr select
#' @importFrom tidyselect where
#' @importFrom stats quantile
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_density
#' @importFrom ggplot2 scale_x_log10
#' @importFrom scales trans_format
#' @importFrom scales math_format
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 element_text
#'
#' @return ggplot objects for each fluorophore containing the various clusters
#' @export
#'
#' @examples NULL
Luciernaga_Brightness <- function(fluorophore.name, data,
  fluorophore.column, cluster.column, downsample=TRUE,
  subsample = NULL, detector=NULL, reference=NULL, 
  clearance=0.02, Scaled = TRUE, maxtik=1e6, legend="right"){
  
    TheFluorophore <- fluorophore.name
    TheData <- data |> filter(.data[[fluorophore.column]] %in% TheFluorophore)
    TheTable <- data.frame(table(TheData[[cluster.column]]), check.names = FALSE)
    colnames(TheTable)[1] <- "Cluster"
    colnames(TheTable)[2] <- "Count"
    TheSlice <- TheTable |> arrange(Count) |> slice(1) |> pull(Count)

    if (downsample == TRUE) {
      if (is.null(subsample)){
       TheData <- TheData |> group_by(.data[[cluster.column]]) |>
         slice_sample(n=TheSlice, replace = FALSE) |> ungroup()
      } else {TheData <- TheData |> group_by(.data[[cluster.column]]) |>
        slice_sample(n=subsample, replace = FALSE) |> ungroup()
      }
    }

    if (is.null(reference)){TheDetector <- detector
    } else {
      if(!is.data.frame(reference)){
        CSV <- read.csv(reference, check.names = FALSE)
      } else {CSV <- reference}
    internalstrings <- c(" ", ".", "_", "-A")
    CSV$Fluorophore <- NameCleanUp(CSV$Fluorophore, removestrings=internalstrings)
    CSV$Detector <- NameCleanUp(CSV$Detector, removestrings=internalstrings)
    TheDetector <- CSV |> dplyr::filter(Fluorophore %in% TheFluorophore) |>
      pull(Detector)
    }

    Values <- TheData |> select(TheDetector) |> as.matrix()
    theXmin <- Values %>% quantile(., 0.00)
    theXmax <- Values %>% quantile(., 1.00)
    theXmin <- theXmin - abs((clearance*theXmin))
    theXmax <- theXmax + (clearance*theXmax)

    if (Scaled == TRUE){
      if (theXmax > 1000000){custom_breaks <- c(1e3, 1e4, 1e5, 1e6, 1e7, 1e8)
      } else if (theXmax > 100000){custom_breaks <- c(1e3, 1e4, 1e5, 1e6, 1e7)
      } else if (theXmax > 10000){custom_breaks <- c(1e3, 1e4, 1e5, 1e6)
      } else {custom_breaks <- c(1e3, 1e4, 1e5)}

      if (theXmin < -100000){lower_breaks <- c(-1e6, -1e5, -1e4, -1e3, 0.1)
      } else if (theXmin < -10000){lower_breaks <- c(-1e5, -1e4, -1e3, 0.1)
      } else if (theXmin < -1000){lower_breaks <- c(-1e4, -1e3, 0.1)
      } else if (theXmin < 0){lower_breaks <- c(-1e3, 0.1)
      } else {lower_breaks <- c(1)}

      custom_breaks <- c(lower_breaks, custom_breaks)
      HighEnd <- tail(custom_breaks, 1)
      LowEnd <- head(custom_breaks, 1)     

      plot <- ggplot(TheData, aes(x=.data[[TheDetector]],
         fill=.data[[cluster.column]])) + geom_density(alpha=0.5) +
          scale_x_log10(limits = c(LowEnd, HighEnd), breaks = custom_breaks, 
        labels = scales:::trans_format("log10", scales::math_format(10^.x)))  +
        labs(title=TheFluorophore, x=TheDetector, y="Frequency") + theme_bw() +
        theme(axis.title.x=element_text(face="plain"), axis.title.y=element_text(
            face="plain"), legend.position=legend)
      } else {
        plot <- ggplot(TheData, aes(x=.data[[TheDetector]], fill=Cluster)) +
          geom_density(alpha=0.5) + xlim(theXmin, theXmax) +
          labs(title=TheFluorophore, x=TheDetector, y="Frequency") + theme_bw() +
          theme(axis.title.x=element_text(face="plain"), axis.title.y=element_text(
              face="plain"), legend.position=legend)
    }

 return(plot)
}


