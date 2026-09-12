#' Check gate placement for individual .fcs files in a GatingSet
#'
#' @param x A GatingSet object
#' @param sample.name The .fcs keyword that contains an unique name for sample
#' @param removestrings Character values to remove from the sample.name
#' @param subset The GatingSet subset that you want to visualize for data
#'   plotting, "root" is the default.
#' @param gtFile The data.table imported .csv file containing gating template.
#' @param DesiredGates A vector of gates that you want plotted, for example
#'   Desired <- c("nonDebris, "lymphocytes")
#' @param outpath Location to store the generated .pdf file
#' @param filename Default NULL, overrides name
#' @param returnType Whether to return "pdf", "patchwork" or "plots".
#' @param bins Argument to geom_hex for number of bins to visualize the plotted
#'   data density.
#' @param therows Number of desired rows for the .pdf file
#' @param thecolumns Number of desired columns for the .pdf file
#' @param width Desired page width
#' @param height Desired page height
#' @param clearance A buffer area around the plot edge
#' @param optionalX When gtFile is NULL, provides x-axis argument for the
#'   subset gated population
#' @param optionalY When gtFile is NULL, provides y-axis argument for the
#'   subset gated population
#' @param optionalGate Default NULL, if using optional arguments and correct X
#'   and Y, the gate arg
#' @param optionalName Default NULL, alternatively sets title for "plots"
#' @param plotname TBD
#'
#' @importFrom flowWorkspace keyword gs_pop_get_data
#' @importFrom flowCore exprs
#' @importFrom dplyr pull
#' @importFrom purrr map
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
#' FCS_Files <- list.files(
#'   path = File_Location, pattern = ".fcs",
#'   full.names = TRUE
#' )
#' UnstainedFCSFiles <- FCS_Files[grep("Unstained", FCS_Files)]
#' UnstainedCells <- UnstainedFCSFiles[-grep("Beads", UnstainedFCSFiles)]
#' MyCytoSet <- load_cytoset_from_fcs(
#'   UnstainedCells[1],
#'   truncate_max_range = FALSE, transformation = FALSE
#' )
#' MyGatingSet <- GatingSet(MyCytoSet)
#' MyGates <- fread(file.path(path = File_Location, pattern = "Gates.csv"))
#' MyGatingTemplate <- gatingTemplate(MyGates)
#' gt_gating(MyGatingTemplate, MyGatingSet)
#' removestrings <- c("DR_", "Cells", ".fcs", "-", " ")
#' StorageLocation <- file.path("C:", "Users", "JohnDoe", "Desktop")
#'
#' IndividualPlot <- Utility_GatingPlots(
#'   x = MyGatingSet[[1]],
#'   sample.name = "GUID", removestrings = removestrings,
#'   gtFile = MyGates, DesiredGates = NULL,
#'   outpath = StorageLocation, filename = NULL,
#'   returnType = "patchwork"
#' )
#'
Utility_GatingPlots <- function(x, sample.name, removestrings,
                                subset = "root", gtFile = NULL,
                                DesiredGates = NULL, outpath = NULL,
                                filename = NULL, returnType, bins = 270,
                                therows = 2, thecolumns = 2, width = 7,
                                height = 9, clearance = 0.2, optionalX = NULL,
                                optionalY = NULL, optionalGate = NULL,
                                optionalName = NULL, plotname = FALSE) {

  if (plotname == TRUE) {
    optionalName <- sample.name
  }

  # Setting up individual file name
  if (is.null(outpath)) {
    outpath <- getwd()
  }

  AggregateName <- Luciernaga:::NameForSample(
    x = x, sample.name = sample.name,
    removestrings = removestrings
  )

  # Retrieve Gates present in the GatingSet.
  if (!is.null(gtFile)) {
    TheXYZgates <- gtFile |> pull(alias)
    if ("*" %in% TheXYZgates) {
      gtFile <- TemplateConverter(gtFile)
      TheXYZgates <- gtFile |> pull(alias)
    }
  } else {
    message("No gating reference file provided, returning provided arguments")
    TheXYZgates <- NULL
  }

  # Desired Gate
  if (!is.null(TheXYZgates) && !is.null(DesiredGates)) {
    TheXYZgates <- intersect(DesiredGates, TheXYZgates)
  }

  if (!is.null(optionalName)) {
    if (length(optionalName) == 2) {
      first <- optionalName[[1]]
      second <- optionalName[[2]]
      first <- keyword(x, first)
      second <- keyword(x, second)
      name <- paste(first, second, sep = "_")
    } else {
      name <- keyword(x, optionalName)
    }
  } else {
    name <- NULL
  }

  # Pulling Gating Set Data
  # Plot Generation
  if (!is.null(TheXYZgates)) {
    ff <- gs_pop_get_data(x, subset)
    df <- exprs(ff[[1]])
    TheDF <- data.frame(df, check.names = FALSE)
    x2 <- x

    # x <- TheXYZgates[1]
    # data <- x2
    CompiledPlots <- map(
      .x = TheXYZgates, .f = GatePlot, data = x2,
      TheDF = TheDF, gtFile = gtFile, bins = bins,
      clearance = clearance, name = name
    )
  } else {
    CompiledPlots <- Utility_IterativeGating(
      x = x, subset = subset, xValue = optionalX,
      yValue = optionalY, gate = optionalGate,
      sample.name = sample.name,
      removestrings = removestrings, bins = bins,
      plotname = plotname
    )
  }

  if (!is.null(filename)) {
    AggregateName <- filename
  }

  if (returnType == "pdf") {
    AssembledPlots <- Utility_Patchwork(
      x = CompiledPlots, filename = AggregateName,
      outfolder = outpath, returntype = "pdf",
      therows = therows, thecolumns = thecolumns,
      width = width, height = height
    )
  } else if (returnType == "patchwork") {
    AssembledPlots <- Utility_Patchwork(
      x = CompiledPlots, filename = AggregateName,
      outfolder = outpath, returntype = "patchwork",
      therows = therows, thecolumns = thecolumns,
      width = width, height = height
    )
  } else if (returnType == "plots") {
    AssembledPlots <- CompiledPlots
  } else {
    AssembledPlots <- CompiledPlots
  }

  return(AssembledPlots)
}

#' Handles the plus-slash-minus verbiage to avoid crash out of the plots
#'
#' @param gtFile TBD
#'
#' @importFrom dplyr bind_rows slice add_row pull filter
#'
#' @noRd
TemplateConverter <- function(gtFile) {

  gtFile1 <- gtFile
  Values <- gtFile1 |> pull(alias)
  Indices <- which(Values == "*")

  Expansions <- gtFile1 |>
    dplyr::filter(alias %in% "*") |>
    pull(dims)

  for (i in seq_along(Indices)) {
    Values <- gtFile1 |> pull(alias)
    Current <- which(Values == "*")[1]

    Row <- gtFile1[Current, ]
    Differential <- Row[[1, 4]]
    Rows <- bind_rows(Row, Row)
    Rows[1, 1] <- paste0(Differential, "+")
    Rows[2, 1] <- paste0(Differential, "-")

    gtFile1 <- gtFile1 |>
      slice(-Current) |>
      add_row(Rows, .before = Current)
  }

  Parents <- gtFile1 |> pull(parent)

  Expanded <- c(paste0(Expansions, "+"), paste0(Expansions, "-"))
  Utilized <- Expanded[!Expanded %in% Parents]
  Removal <- Utilized[-length(Utilized)]
  gtFile1 <- gtFile1 |> filter(!alias %in% Removal)

  return(gtFile1)
}