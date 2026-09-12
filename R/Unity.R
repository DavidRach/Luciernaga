
#' Internal to Utility_UnityPlots
#'
#' @importFrom Biobase pData
#' @importFrom purrr map
#' @importFrom patchwork wrap_plots
#' @importFrom grDevices pdf dev.off
#' @importFrom stats quantile
#' @importFrom flowWorkspace keyword gs_pop_get_data
#' @importFrom flowCore exprs
#' @importFrom dplyr select
#' @importFrom tidyr all_of
#' @importFrom ggcyto as.ggplot ggcyto
#' @importFrom ggplot2 aes geom_point geom_hex coord_cartesian
#'  element_line theme_bw theme element_line geom_vline
#'
#' @return An internal value
#'
#' @noRd
Unity <- function(x, TheY, TheX, marginsubset, gatesubset, sample.name, removestrings,
                  clearance, bins, gatelines, reference, cartesian=TRUE){

  if (length(sample.name) == 2){
      first <- sample.name[[1]]
      second <- sample.name[[2]]
      first <- keyword(x, first)
      second <- keyword(x, second)
      name <- paste(first, second, sep="_")
    } else {name <- keyword(x, sample.name)}
  
  name <- NameCleanUp(name = name, removestrings)
  
  if (inherits(x, "flowFrame")){
    message("flowFrame detected, unable to use gatesubset, using root instead")
    df <- exprs(x)
  } else {
    ff <- gs_pop_get_data(x, marginsubset)
    df <- exprs(ff[[1]])
  }

  TheDF <- data.frame(df, check.names = FALSE)

  if (!TheX == TheY) {
    YExprsData <- TheDF |> select(all_of(TheY)) |> pull()
    theYmin <- YExprsData %>% quantile(., 0.001)
    theYmax <- YExprsData %>% quantile(., 0.999)
    theYmin <- theYmin - abs((clearance*theYmin))
    theYmax <- theYmax + (clearance*theYmax)

    XExprsData <- TheDF |> select(all_of(TheX)) |> pull()
    theXmin <- XExprsData %>% quantile(., 0.001)
    theXmax <- XExprsData %>% quantile(., 0.999)
    theXmin <- theXmin - abs((clearance*theXmin))
    theXmax <- theXmax + (clearance*theXmax)
  } else (stop("TheX and TheY have the same value"))

  if (inherits(x, "flowFrame")){
    ff1 <- x
  } else {ff1 <- gs_pop_get_data(x, gatesubset)}

  if (BiocGenerics::nrow(ff1) < 200) {

    if (cartesian == TRUE){
    Plot <- as.ggplot(ggcyto(ff1, aes(x = .data[[TheX]], y = .data[[TheY]]),
     subset = "root") + geom_point(size = 2, alpha = 0.8) + coord_cartesian(
       xlim = c(theXmin, theXmax), ylim = c(theYmin, theYmax), default = TRUE) +
       theme_bw() + labs(title = name) + theme(strip.background = element_blank(),
     strip.text.x = element_blank(), panel.grid.major = element_line(
     linetype = "blank"), panel.grid.minor = element_line(linetype = "blank"),
     axis.title = element_text(size = 10, face = "bold"),legend.position = "none"))
     } else {
      Plot <- as.ggplot(ggcyto(ff1, aes(x = .data[[TheX]], y = .data[[TheY]]),
     subset = "root") + geom_point(size = 2, alpha = 0.8) +
       theme_bw() + labs(title = name) + theme(strip.background = element_blank(),
     strip.text.x = element_blank(), panel.grid.major = element_line(
     linetype = "blank"), panel.grid.minor = element_line(linetype = "blank"),
     axis.title = element_text(size = 10, face = "bold"),legend.position = "none"))
     }
      
    } else {
    
    if (cartesian == TRUE){
    Plot <- as.ggplot(ggcyto(ff1, aes(x = .data[[TheX]], y = .data[[TheY]]),
     subset = "root") + geom_hex(bins=bins) + coord_cartesian(
       xlim = c(theXmin, theXmax), ylim = c(theYmin, theYmax), default = TRUE) +
       theme_bw() + labs(title = name) +
       theme(strip.background = element_blank(), strip.text.x = element_blank(),
     panel.grid.major = element_line(linetype = "blank"),
     panel.grid.minor = element_line(linetype = "blank"),
     axis.title = element_text(size = 10, face = "bold"),
     legend.position = "none"))
     } else {
      Plot <- as.ggplot(ggcyto(ff1, aes(x = .data[[TheX]], y = .data[[TheY]]),
     subset = "root") + geom_hex(bins=bins) +
       theme_bw() + labs(title = name) +
       theme(strip.background = element_blank(), strip.text.x = element_blank(),
     panel.grid.major = element_line(linetype = "blank"),
     panel.grid.minor = element_line(linetype = "blank"),
     axis.title = element_text(size = 10, face = "bold"),
     legend.position = "none"))
     }
    }

  if (gatelines == TRUE){
  Value <- reference |> dplyr::filter(specimen %in% name) |>
    select(all_of(TheX)) %>% pull(.)
 
  Plot <- Plot +
    geom_vline(xintercept = c(seq(0,10000,25)), colour = "gray", alpha=0.1) +
    geom_vline(xintercept = c(seq(0,10000,2)), colour = "white", alpha = 0) +
    geom_vline(xintercept = Value, colour = "red")
  }

  return(Plot)
}