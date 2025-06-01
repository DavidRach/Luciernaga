#' Adds openCyto boundary gates at designated location, selectively cleaning out margin
#' events that mess with various algorithms. Can return template for redo editing, as well
#' as returns plots for troubleshooting. 
#' 
#' @param gs The GatingSet object you want to clean margins on
#' @param desiredCols Provide a vector of fluorophore names to clean margins for
#'  (see colnames(gs)), default NULL cleans margins for all fluorophores
#' @param subset Subset at which to start margin cleanup, default is root. 
#' @param themin Default is NULL, provide a numeric value to establish a lower 
#' boundary to clean up small debris. 
#' @param themax Provide a numeric value to set the upper boundary to exclude
#'  margin events.
#' @param returnTemplate Default is FALSE, returns assembled openCyto gating 
#' template allowing for adjustments that can be brought in with importTemplate
#' and inpath
#' @param importTemplate Default is FALSE, when TRUE, retrieves previous 
#' returnTemplate .csv from the inpath location and uses for the openCyto gating
#' @param returnPlots Default is FALSE, when TRUE, returns a pdf of margin clean
#'  events to the outpath to verify didn't cut of the population of interest. 
#' @param Verbose Default is FALSE, will print to console the frequency of retained
#'  cells after margin cleanup for each specimen
#' @param inpath Default NULL, alternatively a file.path to the template.csv
#'  being imported for openCyto gating
#' @param outpath Default NULL, provide a file.path to desired location to store 
#' either the returnTemplate or the returnPlot objects
#' @param filename Default NULL, alternatively set a name for returnTemplate or
#' returnPlot objects
#' @param yaxis Sets yaxis fluorophore on returnPlots, the default NULL utilizes
#'  the first fluorophore in the panel
#' @param sample.name Used when returningPlots, default NULL uses TUBENAME as the keyword
#' value provided when retrieving sample name for plot titles.
#' @param inverse.transform Default is FALSE, retaining input GatingSet transformation setting
#' 
#' @importFrom BiocGenerics colnames
#' @importFrom flowWorkspace gs_get_pop_paths
#' @importFrom data.table fread
#' @importFrom purrr map2
#' @importFrom dplyr bind_rows
#' @importFrom utils write.csv
#' @importFrom openCyto gatingTemplate
#' @importFrom openCyto gt_gating
#' @importFrom flowWorkspace gs_pop_get_count_fast
#' @importFrom dplyr filter
#' @importFrom flowWorkspace gs_pop_get_data
#' @importFrom flowWorkspace GatingSet
#' 
#' @return A margin cleaned GatingSet, alternatively an openCyto template or visualized 
#' cleaned data
#' 
#' @export
MarginClean <- function(gs, desiredCols=NULL, subset="root", themin=NULL, themax,
  returnTemplate=FALSE, importTemplate=FALSE, returnPlots=FALSE, Verbose=FALSE, 
  inpath=NULL, outpath=NULL, filename=NULL, yaxis=NULL, sample.name=NULL,
  inverse.transform=FALSE){

  if (is.null(themin)){
    thegatingargs <- paste0("max=", themax)
  } else {thegatingargs <- paste0("min=", themin, ",max=", themax)}

  TheCols <- colnames(gs)

  if (is.null(desiredCols)){
        These <- TheCols[-grep("Time|FS|SC|SS|Original|W$|H$|AF", TheCols)]
      } else {These <- TheCols[TheCols %in% desiredCols]
      }
  
  ExistingGates <- gs_get_pop_paths(gs, path="auto")
  DesiredGateCheck <- subset %in% ExistingGates

  if (DesiredGateCheck==FALSE){
    stop("Make sure subset is an existing gate in your GatingSet")
  }

  if (!importTemplate == TRUE){
  FileLocation <- system.file("extdata", package = "Luciernaga")
  UnmixedGates <- fread(file.path(path = FileLocation, pattern = 'GatesUnmixed.csv'))
  Template <- UnmixedGates[1,]
  Template[,1] <- These[1]
  Template[,3] <- subset
  Template[,4] <- These[1]
  Template[,5] <- "boundary"
  Template[,6] <- thegatingargs

  thex <- These[-1]
  they <- These[-length(These)]

  Data <- map2(.x=thex, .y=they, .f=TemplateAssembly, template=Template) |> bind_rows()
  Data <- bind_rows(Template, Data)

  if (returnTemplate == TRUE){
    
    if (!is.null(filename)){TheFileName <- paste0(filename, "_template.csv")
      } else {TheFileName <- "MarginClean_template.csv"}

    if (!is.null(outpath)){outpath <- outpath
    } else {outpath <- getwd()}

    StorageLocation <- file.path(outpath, TheFileName)
    write.csv(Data, StorageLocation, row.names=FALSE)
    #return(Data)
    }
  } else {
  if (!is.null(inpath)){Data <- fread(inpath)
    } else {
      stop("Please provide a file.path to the template to the inpath argument")
      }
  }

  PullTheLever <- gatingTemplate(Data)
  gt_gating(PullTheLever, gs)

  ExistingGates <- gs_get_pop_paths(gs, path="auto")
  FirstGate <- ExistingGates[2]
  FinalGate <- ExistingGates[length(ExistingGates)]
  SortThese <- c(FirstGate, FinalGate)
  Dataset <- gs_pop_get_count_fast(gs, statistic="count")
  Dataset$Population <- basename(Dataset$Population)
  Dataset$Parent <- basename(Dataset$Parent)
  Dataset <- Dataset |> filter(Population %in% SortThese)
  if (Verbose==TRUE){print(Dataset)}

  if (returnPlots==TRUE){

  if (!is.null(filename)){TheFileName <- paste0(filename, "_plot.csv")
      } else {TheFileName <- "MarginClean_plot.csv"}

    if (!is.null(outpath)){outpath <- outpath
    } else {outpath <- getwd()}

    StorageLocation <- file.path(outpath, TheFileName)

    if (!is.null(yaxis)){TheY <- yaxis
    } else {TheY <- FirstGate}

    if (!is.null(sample.name)){thesample <- sample.name
    } else {thesample <- "TUBENAME"}

    Utility_NbyNPlots(x=gs[1], y=TheY, sample.name="TUBENAME", marginsubset="root",
     gatesubset=FinalGate, bins=100, clearance=0.1, gatelines=FALSE,  reference=NULL,
      outpath=outpath, returntype="pdf", removestrings=".fcs", filename=TheFileName,
      experiment="Test", condition="Test")
  }

  MarginClean <- gs_pop_get_data(gs, FinalGate, inverse.transform = inverse.transform)
  MarginClean_GS <- GatingSet(MarginClean)
  #plot(MarginClean_GS)
  #pData(MarginClean_GS)
  return(MarginClean_GS)
}

#' Internal for MarginClean, assembles openCyto template
#' 
#' @param x The fluorophore list for alias and dims
#' @param y The fluorophore list for parent gate
#' @param template The openCyto template for the first non-root gate
#' 
#' @return A data.table row with updated gate information
#' 
#' @noRd
TemplateAssembly <- function(x, y, template){
  template[1,1] <- x
  template[1,3] <- y
  template[1,4] <- x
  return(template)
}