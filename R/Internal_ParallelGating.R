#' Internal for Utility_ParallelNxNPlots
#'
#' @importFrom flowWorkspace keyword
#' @importFrom flowWorkspace gs_pop_get_data
#' @importFrom flowCore exprs
#' @importFrom patchwork wrap_plots
#' @importFrom patchwork plot_spacer
#' @importFrom purrr map
#' @importFrom dplyr select
#' @importFrom dplyr pull
#' @importFrom dplyr mutate
#' @importFrom ggplot2 ggplot
#'
#' @keywords internal

Internal_ParallelGating <- function(x, x_ff, y_ff, TheDF, yValue, columnlist, gatelines,
                                     reference, clearance, bins, AltNameX, AltNameY, colorX, colorY) {

  if (yValue == x){stop("x equals yValue and can't be plotted")}

  xValue <- x

  if (!grepl("FSC|SSC", yValue)) {
    ExprsData <- TheDF %>% select(all_of(yValue)) %>% pull()
    theYmin <- ExprsData %>% quantile(., 0.001)
    theYmax <- ExprsData %>% quantile(., 0.999)
    theYmin <- theYmin - abs((clearance*theYmin))
    theYmax <- theYmax + (clearance*theYmax)}

  if (!grepl("FSC|SSC", xValue)) {
    ExprsData <- TheDF %>% select(all_of(xValue)) %>% pull()
    theXmin <- ExprsData %>% quantile(., 0.001)
    theXmax <- ExprsData %>% quantile(., 0.999)
    theXmin <- theXmin - abs((clearance*theXmin))
    theXmax <- theXmax + (clearance*theXmax)}


  if (!exists("theYmax") || !exists("theXmax")){
    stop("Either theYmax or theXmax didn't exist, and since I didn't think
     it relavant to duplicate this code in the parallel NxN plot when coding,
          the function now crashed ")
  } else {

    x_ffXdata <- exprs(x_ff[[1]]) %>% data.frame(check.names = FALSE) %>% select(all_of(xValue))
    x_ffYdata <- exprs(x_ff[[1]]) %>% data.frame(check.names = FALSE) %>% select(all_of(yValue))
    Thex_ff <- cbind(x_ffXdata, x_ffYdata) %>% mutate(specimen = AltNameX)

    y_ffXdata <- exprs(y_ff[[1]]) %>% data.frame(check.names = FALSE) %>% select(all_of(xValue))
    y_ffYdata <- exprs(y_ff[[1]]) %>% data.frame(check.names = FALSE) %>% select(all_of(yValue))
    They_ff <- cbind(y_ffXdata, y_ffYdata) %>% mutate(specimen = AltNameY)

    TheData <- rbind(Thex_ff, They_ff)
    TheData$specimen <- as.factor(TheData$specimen)

    specimen_counts <- table(TheData$specimen)
    sorted_specimens <- names(sort(desc(specimen_counts)))
    TheData$specimen <- factor(TheData$specimen, levels = sorted_specimens)

    # Attempted Work Around for Specifying Colors while adjusting what population is displayed forward.
    Xscheme <- cbind(AltNameX, colorX)
    Yscheme <- cbind(AltNameY, colorY)
    ColorFrame <- rbind(Xscheme, Yscheme)
    ColorFrame <- data.frame(ColorFrame, check.names = FALSE) %>% rename(specimen = AltNameX)
    ColorFrame$specimen <- factor(ColorFrame$specimen, levels = sorted_specimens)
    ColorFrame <- ColorFrame[order(ColorFrame$specimen), ]
    color1 <- ColorFrame[1,2]
    color2 <- ColorFrame[2,2]

    Plot <- ggplot(TheData, aes(x=.data[[xValue]], y = .data[[yValue]], fill = specimen)) +
      geom_hex(bins=bins, alpha = 0.5) + scale_fill_manual(values = c(color1, color2)) +
      coord_cartesian(xlim = c(theXmin, theXmax), ylim = c(theYmin, theYmax)) + theme_bw() +
      labs(title = NULL) + theme(strip.background = element_blank(),
                                 strip.text.x = element_blank(), panel.grid.major = element_line(linetype = "blank"),
                                 panel.grid.minor = element_line(linetype = "blank"),
                                 axis.title = element_text(size = 10, face = "bold"),
                                 legend.position = "none")

    #if (gatelines == TRUE){Value <- reference[reference$specimen == name, xValue]
    #Plot <- Plot + geom_vline(xintercept = c(seq(0,200,25)), colour = "gray") +
    #  geom_vline(xintercept = Value, colour = "red")}
  }

  return(Plot)
}
