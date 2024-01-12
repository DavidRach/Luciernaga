#' Generate ggplot2s for different channels
#'
#' @param x Passed channel
#' @param ff The dataframe data for that sample
#' @param yValue What wanted on the yaxis
#' @param columnlist list all channels with x removed
#' @param TheDF the external limits settings
#' @param gatelines whether to plot the ModernCutoffLines from reference dataframe
#' @param reference location ModernCutoff dataframe.
#'
#' @importFrom flowCore keyword
#' @importFrom ggcyto ggcyto
#' @importFrom ggcyto as.ggplot
#'
#' @return NULL
#' @export
#'
#' @examples NULL
Utility_GeneralGating <- function(x, ff, yValue, columnlist, TheDF, gatelines, reference = NULL) {

  if (ycolumn == x){stop("Error: x equals ycolumn and can't be plotted")}

  xValue <- x

  if (!grepl("FSC|SSC", yValue)) {ExprsData <- TheDF %>% select(all_of(yValue)) %>% pull()
  theYmin <- ExprsData %>% quantile(., 0.001)
  theYmax <- ExprsData %>% quantile(., 0.999)
  theYmin <- theYmin - abs((clearance*theYmin))
  theYmax <- theYmax + (clearance*theYmax)}

  if (!grepl("FSC|SSC", xValue)) {ExprsData <- TheDF %>% select(all_of(xValue)) %>% pull()
  theXmin <- ExprsData %>% quantile(., 0.001)
  theXmax <- ExprsData %>% quantile(., 0.999)
  theXmin <- theXmin - abs((clearance*theXmin))
  theXmax <- theXmax + (clearance*theXmax)}


if (!exists("theYmax") || !exists("theXmax")){Plot <- as.ggplot(ggcyto(ff, aes(x = .data[[xValue]], y = .data[[yValue]]), subset = "root") + geom_hex(bins=bins) + theme_bw() + labs(title = NULL) +
                                                                  theme(strip.background = element_blank(), strip.text.x = element_blank(), panel.grid.major = element_line(linetype = "blank"),
                                                                        panel.grid.minor = element_line(linetype = "blank"), axis.title = element_text(size = 10, face = "bold"), legend.position = "none"))

  if (gatelines == TRUE){Value <- reference[reference$specimen == name, xValue]
    Plot <- Plot + geom_vline(xintercept = c(seq(0,200,25)), colour = "gray") + geom_vline(xintercept = Value, colour = "red")
  }

} else {Plot <- as.ggplot(ggcyto(ff, aes(x = .data[[xValue]], y = .data[[yValue]]), subset = "root") + geom_hex(bins=bins) + coord_cartesian(xlim = c(theXmin, theXmax), ylim = c(theYmin, theYmax), default = TRUE) +
                             theme_bw() + labs(title = NULL) + theme(strip.background = element_blank(), strip.text.x = element_blank(), panel.grid.major = element_line(linetype = "blank"),
                                                                                         panel.grid.minor = element_line(linetype = "blank"), axis.title = element_text(size = 10, face = "bold"), legend.position = "none"))

if (gatelines == TRUE){Value <- reference[reference$specimen == name, xValue]
Plot <- Plot + geom_vline(xintercept = c(seq(0,200,25)), colour = "gray") + geom_vline(xintercept = Value, colour = "red")}
}

tryCatch({rm("theXmin", "theXmax", "theYmin", "theYmax")})

return(Plot)
}
