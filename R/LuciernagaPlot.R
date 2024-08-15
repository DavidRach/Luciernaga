#' Generate Luciernaga plot outputs from List of List
#'
#' @param ListOfList A list containing the list (BrightnessPlot, LinePlot,
#' HeatmapPlot, CosinePlot)
#' @param PlotType Whether "html" or "pdf"
#' @param ReturnFolder Location to store file
#' @param CurrentExperiment Name of Current Experiment
#'
#' @importFrom plotly ggplotly
#' @importFrom patchwork wrap_plots
#' @importFrom htmltools tagList
#' @importFrom htmltools save_html
#' @importFrom purrr map
#' @importFrom purrr flatten
#'
#' @return A file containing the Luciernaga plots in the specified format
#' @export
#'
#' @examples NULL
LuciernagaPlot <- function(ListOfList, PlotType, ReturnFolder, CurrentExperiment){

  #ItemSelect(Combination, n)

  indices <- length(ListOfList[[1]])

  Ultimate <- map(1:indices, ~ItemSelect(ListOfList, .))

  if (PlotType == "pdf"){
    theList <- Ultimate

    PDFPath <- paste0(ReturnFolder, CurrentExperiment, '_All', '.pdf')

    pdf(file = PDFPath, width = 8.5, height = 11) #Optional Adjustments for Second

    SubplotsPDF <- function(i, data){
      Components <- flatten(data[i])

      layout <- "
                AAAAAA
                BBBBBB
                CCC#DD
                "

      Patchworked <- wrap_plots(Components, design = layout, heights = c(1, 1, 2),
                                widths = c(1, 1, 1, 1, 1, 1))
      Patchworked <- Patchworked + plot_annotation(caption =
                                                     "Made with Luciernaga")
      print(Patchworked)
    }

    Rendered <- map(1:indices, ~SubplotsPDF(data = theList, .))

    dev.off()
  }

  if (PlotType == "html"){
    theList <- Ultimate #List of Lists

    #Subplots(data = Ultimate, i = 2)
    Rendered <- map(1:indices, ~Subplots(data = theList, .))

    HtmlPath <- paste0(ReturnFolder, CurrentExperiment, '_All', '.html')

    htmltools::save_html(html = Rendered, file = HtmlPath)
  }
}

ItemSelect <- function(ListOfList, n) {
  result <- lapply(ListOfList, function(innerList) innerList[[n]])
  #result <- list(result)
  return(result)
}

Subplots <- function(i, data) {
  Components <- flatten(data[i])

  plotlyobjs <- lapply(Components, ggplotly)

  subplot <- htmltools::tagList(
    lapply(plotlyobjs, function(x) {
      div(
        x,
        style = "float:left; width:50%;",
        tags$br()
      )
    }),
    tags$br(style = "clear:both;"), # Clear both to prevent float overlap
    tags$br(), tags$br(), tags$br() # Additional breaks as needed
  )

  return(subplot)
}
