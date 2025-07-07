#' Check gate placement for individual .fcs files in a GatingSet
#'
#' @param x A GatingSet object
#' @param sample.name The .fcs keyword that contains an unique name for that sample
#' @param removestrings A string of character values to remove from the sample.name
#' @param subset The GatingSet subset that you want to visualize for data plotting,
#' "root" is the default.
#' @param bins Argument to geom_hex for number of bins to visualize the plotted
#' data density.
#' @param clearance A buffer area around the plot edge
#' @param gtFile The data.table imported .csv file containing the gating template.
#' @param DesiredGates A vector of gates that you want plotted, for example
#' Desired <- c("nonDebris, "lymphocytes")
#' @param returnType Whether to return "pdf", "patchwork" or "plots".
#' @param thecolumns Number of desired columns for the .pdf file
#' @param therows Number of desired rows for the .pdf file
#' @param width Desired page width
#' @param height Desired page height
#' @param outpath Location to store the generated .pdf file
#' @param filename Default NULL, overrides name
#' @param optionalX When gtFile is NULL, provides x-axis argument for the subset gated population
#' @param optionalY When gtFile is NULL, provides y-axis argument for the subset gated population
#' @param optionalGate Default NULL, if using optional arguments and correct X and Y, the gate arg
#' @param optionalName Default NULL, alternatively sets title for "plots"
#' 
#' @importFrom flowWorkspace keyword
#' @importFrom flowWorkspace gs_pop_get_data
#' @importFrom flowCore exprs
#' @importFrom dplyr pull
#' @importFrom purrr map
#' @importFrom stats alias
#'
#' @return Additional information to be added
#' @export
#'
#' @examples
#'
#' library(flowCore)
#' library(flowWorkspace)
#' library(openCyto)
#' library(data.table)
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
#' IndividualPlot <- Utility_GatingPlots(x=MyGatingSet[[1]],
#'   sample.name = "GUID",removestrings = removestrings,
#'   gtFile = MyGates, DesiredGates = NULL,
#'   outpath = StorageLocation, filename=NULL, 
#'   returnType = "patchwork")
#'
Utility_GatingPlots <- function(x, sample.name, removestrings,
  subset="root", gtFile=NULL, DesiredGates = NULL, outpath = NULL,
  filename=NULL, returnType, bins=270, therows=2, thecolumns=2,
  width=7, height=9, clearance=0.2, optionalX=NULL, optionalY=NULL,
  optionalGate=NULL, optionalName=NULL){

  # Setting up individual file name
  if (is.null(outpath)){outpath <- getwd()}

  AggregateName <- Luciernaga:::NameForSample(x=x, sample.name=sample.name,
                                 removestrings=removestrings)

  # Pulling Gating Information
  if(!is.null(gtFile)){
  TheXYZgates <- gtFile |> pull(alias)
  } else {
    message("No gating reference file provided, returning provided arguments")
    TheXYZgates <- NULL}

  # Desired Gate
  if(!is.null(TheXYZgates) && !is.null(DesiredGates)){
    TheXYZgates <- intersect(DesiredGates, TheXYZgates)}

  # Pulling Gating Set Data

  #Plot Generation
  if(!is.null(TheXYZgates)){

    if (!is.null(optionalName)){
      if (length(optionalName) == 2){
        first <- optionalName[[1]]
        second <- optionalName[[2]]
        first <- keyword(x, first)
        second <- keyword(x, second)
        name <- paste(first, second, sep="_")
      } else {name <- keyword(x, optionalName)}
    } else {name <- NULL}

    ff <- gs_pop_get_data(x, subset)
    df <- exprs(ff[[1]])
    TheDF <- data.frame(df, check.names = FALSE)
    x2 <- x

    CompiledPlots <- map(.x = TheXYZgates, .f = GatePlot, data=x2, TheDF = TheDF,
                        gtFile = gtFile, bins=bins, clearance=clearance, name=name)
    } else {
    CompiledPlots <- Utility_IterativeGating(x=x, subset=subset,
     xValue=optionalX, yValue=optionalY, gate=optionalGate, sample.name=sample.name,
    removestrings=removestrings, bins=bins)
  } 

  if (!is.null(filename)){AggregateName <- filename}

  if (returnType == "pdf"){
    AssembledPlots <- Utility_Patchwork(x=CompiledPlots, filename=AggregateName,
                                        outfolder=outpath, returntype = "pdf",
                                        therows=therows, thecolumns=thecolumns,
                                        width = width, height = height)
  } else if (returnType == "patchwork"){
    AssembledPlots <- Utility_Patchwork(x=CompiledPlots, filename=AggregateName,
                                        outfolder=outpath, returntype = "patchwork",
                                        therows=therows, thecolumns=thecolumns,
                                        width = width, height = height)
  } else if (returnType == "plots"){AssembledPlots <- CompiledPlots}

  return(AssembledPlots)
}

#' Generates called plots from Utility_GatingPlots
#'
#' @param x A specific gate, ex. "nonDebris"
#' @param data A GatingSet object
#' @param TheDF A data.frame object of the flow file's expr data
#' @param gtFile The data.table imported .csv file containing the gating template.
#' @param bins Argument to geom_hex for number of bins to visualize the plotted
#' data density.
#' @param clearance A buffer area around the plot edge
#' @param name Sets the title for the plot, default is NULL
#'
#' @importFrom dplyr filter
#' @importFrom stringr str_split
#' @importFrom dplyr select
#' @importFrom tidyselect all_of
#' @importFrom ggcyto ggcyto
#' @importFrom ggcyto as.ggplot
#' @importFrom ggcyto geom_gate
#' @importFrom ggplot2 ggplot
#'
#' @return A ggplot corresponding to the given inputs
#'
#' @noRd
GatePlot <- function(x, data, TheDF, gtFile, bins=270, clearance = 0.2,
  name){
    i <- x
    gtFile <- data.frame(gtFile, check.names = FALSE)
    RowData <- gtFile |> filter(alias %in% i)
    theSubset <- RowData |> pull(parent)
    theGate <- RowData |> pull(alias)
    theParameters <- RowData |> pull(dims) |> str_split(",", simplify = TRUE)

    theParameters <- gsub("^\\s+|\\s+$", "", theParameters)

    if(length(theParameters) == 2){xValue <- theParameters[[1]]
    yValue <- theParameters[[2]]
    } else if (length(theParameters) == 1){xValue <- theParameters[[1]]
    yValue <- "SSC-A" #or an alternate variable specify
    } else {message(
    "Plotting Parameters for Axis were not 1 or 2, please check the .csv file")
    }


  #Please Note, All the Below Are Raw Values With No Transforms Yet Applied.

  if (!grepl("FSC|SSC", xValue)) {ExprsData <- TheDF |>
    select(all_of(xValue)) |> pull()
  theXmin <- ExprsData %>% quantile(., 0.001)
  theXmax <- ExprsData %>% quantile(., 0.999)
  theXmin <- theXmin - abs((clearance*theXmin))
  theXmax <- theXmax + (clearance*theXmax)}
  if (!grepl("FSC|SSC", yValue)) {ExprsData <- TheDF |>
    select(all_of(yValue)) |> pull()
  theYmin <- ExprsData %>% quantile(., 0.001)
  theYmax <- ExprsData %>% quantile(., 0.999)
  theYmin <- theYmin - abs((clearance*theYmin))
  theYmax <- theYmax + (clearance*theYmax)}

  if (!exists("theYmax") || !exists("theXmax")){
    Plot <- as.ggplot(ggcyto(data, aes(x = .data[[xValue]], y = .data[[yValue]]),
       subset = theSubset) + geom_hex(bins=bins) + geom_gate(theGate) + theme_bw() +
       labs(title = name) + theme(strip.background = element_blank(),
       strip.text.x = element_blank(), panel.grid.major = element_line(
       linetype = "blank"), panel.grid.minor = element_line(linetype = "blank"),
       axis.title = element_text(size = 10, face = "bold"), legend.position = "none"))
  } else {
    Plot <- as.ggplot(ggcyto(data, aes(x = .data[[xValue]], y = .data[[yValue]]), subset = theSubset)) +
      geom_hex(bins=bins) +
      coord_cartesian(xlim = c(theXmin, theXmax), ylim = c(theYmin, theYmax), default = TRUE) +
      geom_gate(theGate) + theme_bw() + labs(title = name) +
      theme(strip.background = element_blank(),
            strip.text.x = element_blank(),
            panel.grid.major = element_line(linetype = "blank"),
            panel.grid.minor = element_line( linetype = "blank"),
            axis.title = element_text(size = 10),
            legend.position = "none")
  }
}
