#' Internal to Utility_UnityPlots
#'
#' @importFrom Biobase pData
#' @importFrom purrr map
#' @importFrom patchwork wrap_plots
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.off
#' @importFrom stats quantile
#' @importFrom flowWorkspace keyword
#' @importFrom flowWorkspace gs_pop_get_data
#' @importFrom flowCore exprs
#' @importFrom dplyr select
#' @importFrom tidyr all_of
#' @importFrom ggcyto as.ggplot
#' @importFrom ggcyto ggcyto
#' @importFrom ggplot2 geom_hex
#'
#' @keywords internal


InternalUnity <- function(x, TheY, TheX, marginsubset, gatesubset, sample.name, removestrings,
                          clearance, bins, gatelines, reference){

  name <- keyword(x, sample.name)
  name <- NameCleanUp(name = name, removestrings)

  ff <- gs_pop_get_data(x, marginsubset)
  df <- exprs(ff[[1]])
  TheDF <- data.frame(df, check.names = FALSE)

  if (!TheX == TheY) {
    YExprsData <- TheDF %>% select(all_of(TheY)) %>% pull()
    theYmin <- YExprsData %>% quantile(., 0.001)
    theYmax <- YExprsData %>% quantile(., 0.999)
    theYmin <- theYmin - abs((clearance*theYmin))
    theYmax <- theYmax + (clearance*theYmax)

    XExprsData <- TheDF %>% select(all_of(TheX)) %>% pull()
    theXmin <- XExprsData %>% quantile(., 0.001)
    theXmax <- XExprsData %>% quantile(., 0.999)
    theXmin <- theXmin - abs((clearance*theXmin))
    theXmax <- theXmax + (clearance*theXmax)
  } else (error("TheX and TheY have the same value"))

  ff1 <- gs_pop_get_data(x, gatesubset)

  if (BiocGenerics::nrow(ff1) < 200) {
    Plot <- as.ggplot(ggcyto(ff1, aes(x = .data[[TheX]], y = .data[[TheY]]),
                             subset = "root") + geom_point(size = 2, alpha = 0.8) + theme_bw() + labs(title = name) +
                        theme(strip.background = element_blank(), strip.text.x = element_blank(),
                              panel.grid.major = element_line(linetype = "blank"),
                              panel.grid.minor = element_line(linetype = "blank"),
                              axis.title = element_text(size = 10, face = "bold"),
                              legend.position = "none"))
  } else {
    Plot <- as.ggplot(ggcyto(ff1, aes(x = .data[[TheX]], y = .data[[TheY]]),
                             subset = "root") + geom_hex(bins=bins) + theme_bw() + labs(title = name) +
                        theme(strip.background = element_blank(), strip.text.x = element_blank(),
                              panel.grid.major = element_line(linetype = "blank"),
                              panel.grid.minor = element_line(linetype = "blank"),
                              axis.title = element_text(size = 10, face = "bold"),
                              legend.position = "none"))
  }



  if (gatelines == TRUE){Value <- reference[reference$specimen == name, TheX]
  Plot <- Plot + geom_vline(xintercept = c(seq(0,200,25)), colour = "gray") +
    geom_vline(xintercept = c(seq(0,200,5)), colour = "white", alpha = 0.1) +
    geom_vline(xintercept = Value, colour = "red")}

  return(Plot)
}
