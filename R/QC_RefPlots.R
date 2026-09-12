#' Internal for QC_UserLibraries
#'
#' @param x x
#' @param Data x
#' @param references x
#' @param refData x
#'
#' @importFrom dplyr filter select pull rename  mutate
#' @importFrom ggplot2 ggplot
#'
#' @return An internal value
#'
#' @noRd
QC_RefPlots <- function(x, Data, references=FALSE, refData=NULL){

  if (!references == TRUE){
    TheValue <- Data |> filter(TheSamples %in% x)
    TheFluorochrome <- TheValue |> select(Fluorochrome) |> unique() %>% pull()
    ThePlot <- ggplot(TheValue, aes(x=Detector, y=value, group=Sample)) +
      geom_line() + theme_bw() + labs(title=paste0(TheFluorochrome, " ", x),
                                      x=NULL, y="Normalized Value") +
      geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
      theme(plot.title = element_text(size = 8), axis.text.x = element_text(
        size = 6, angle = 45), panel.grid = element_blank(), axis.ticks.x = element_blank(),
      axis.title.y =  element_text(size=8)) + scale_x_discrete(breaks = unique(
        TheValue$Detector)[c(TRUE, rep(FALSE, 4))])
  } else {
    TheValue <- Data |> filter(TheSamples %in% x)
    TheFluorochrome <- TheValue |> select(Fluorochrome) |> unique() %>% pull()

    #ReferenceRetrieval
    TheFluorochrome1 <- TheFluorochrome #Name Backup
    TheFluorochrome <- gsub("AF", "Alexa Fluor", gsub("efl", "eFl", gsub(
      "Spk", "Spark", TheFluorochrome)))
    TheFluorochrome <- gsub(" ", "", gsub("-", "", gsub(".", "", fixed=TRUE,
                                                        TheFluorochrome)))
    ReferenceFluorList <- refData |> select(Fluorophore) |> unique() %>% pull()

    if (any(ReferenceFluorList == TheFluorochrome)) { #TheCleanedVersion
      ReferenceData1 <- refData |> filter(Fluorophore %in% TheFluorochrome) |>
        rename(value=TheValue)
      Iterations <- nrow(ReferenceData1)
      MyVector <- seq_len(Iterations)
      ReferenceData1$Detector <- factor(ReferenceData1$Detector, levels=MyVector)
      ReferenceData1 <- ReferenceData1 %>% mutate(TheSamples="Nope")

      ThePlot <- ggplot(TheValue, aes(x=Detector, y=value, group=TheSamples)) +
        geom_line() + theme_bw() + labs(title=paste0(x), x=NULL, y="Normalized Value") +
        geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
        theme(plot.title = element_text(size = 8),
              axis.text.x = element_text(size = 6, angle = 45),
              panel.grid = element_blank(),
              axis.ticks.x = element_blank(),
              axis.title.y =  element_text(size=8)) +
        scale_x_discrete(breaks = unique(TheValue$Detector)[c(TRUE, rep(FALSE, 4))])

      ThePlot <- ThePlot + geom_line(data = ReferenceData1, aes(
        x=Detector, y=value, group=TheSamples), color="red")
    } else {TheValue <- Data |> filter(TheSamples %in% x)
    TheFluorochrome <- TheValue |> select(Fluorochrome) |> unique() %>% pull()
    ThePlot <- ggplot(TheValue, aes(x=Detector, y=value, group=Sample)) + geom_line() +
      theme_bw() + labs(title=paste0(TheFluorochrome, " ", x), x=NULL, y="Normalized Value") +
      geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
      theme(plot.title = element_text(size = 8),
            axis.text.x = element_text(size = 6, angle = 45),
            panel.grid = element_blank(), axis.ticks.x = element_blank(),
            axis.title.y =  element_text(size=8)) +
      scale_x_discrete(breaks = unique(TheValue$Detector)[c(TRUE, rep(FALSE, 4))])
    }
  }
  return(ThePlot)
}