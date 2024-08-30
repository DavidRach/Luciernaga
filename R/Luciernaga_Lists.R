#' Generate Luciernaga plot outputs from List of List
#'
#' @param ListOfList A list containing the returns from Luciernaga_Plots plot option.
#' @param PlotType Whether "html" or "pdf"
#' @param thecolumns The number of columns per page
#' @param therows The number of rows per page
#' @param width Desired page width
#' @param height Desired page height
#' @param ReturnFolder Location to store file
#' @param CurrentExperiment Name of Current Experiment
#'
#' @importFrom purrr map
#' @importFrom purrr transpose
#' @importFrom htmltools save_html
#'
#' @return A file containing the bound Luciernaga plots in specified format
#' @export
#'
#' @examples NULL
Luciernaga_Lists <- function(ListOfList, SecondaryList=NULL, PlotType,
                             thecolumns=2, therows=3, width=7, height=9,
                             ReturnFolder, CurrentExperiment){

  indices <- length(ListOfList[[1]])

  Ultimate <- map(1:indices, ~ItemSelect(ListOfList, .))

  if (!is.null(SecondaryList)){Ultimate <- c(Ultimate, list(SecondaryList))}

  indices <- length(Ultimate)

  Transposed <- transpose(Ultimate) #Test

  if (PlotType == "pdf"){

    Utility_Patchwork(x=Transposed, filename=CurrentExperiment, outfolder=ReturnFolder,
                      thecolumns=thecolumns, therows=therows, width=width, height=height,
                      returntype="pdf", NotListofList = FALSE)
  }

  if (PlotType == "html"){
    Rendered <- map(1:indices, ~Subplots(data = Transposed, .))
    TheFile <- paste0(CurrentExperiment, ".html")
    HtmlPath <- file.path(ReturnFolder, TheFile)
    save_html(html = Rendered, file = HtmlPath)
  }
}

#' Internal for Lucierna_Lists
#'
#' @param ListOfList A list of list
#' @param n Passed number of indices in the above
#'
#' @noRd
ItemSelect <- function(ListOfList, n) {
  result <- lapply(ListOfList, function(innerList) innerList[[n]])
  return(result)
}

#' Internal for Luciernaga_Lists
#'
#' @param i Passed Indicies
#' @param data The Transposed List of Lists
#'
#' @importFrom purrr flatten
#' @importFrom plotly ggplotly
#' @importFrom htmltools tagList
#' @importFrom htmltools div
#'
#' @keywords internal
Subplots <- function(i, data) {
  Components <- flatten(data[i])

  plotlyobjs <- lapply(Components, ggplotly)

  subplot <- tagList(
    lapply(plotlyobjs, function(x) {
      div(
        x,
        style = "float:left; width:50%;",
        tags$br()
      )
    }),
    tags$br(style = "clear:both;"),
    tags$br(), tags$br(), tags$br()
  )

  return(subplot)
}
