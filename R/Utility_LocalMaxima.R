#' Detecting Fluorophore Peak Detectors by Local Maxima
#'
#' @param theX A vector of detectors from 1:n
#' @param theY The corresponding y values corresponding to the measurements of theX
#' @param therepeats Additional values to temporarily add to the edges to allow for peak detection
#' @param w The span around which rolling will happen
#' @param alternatename The cleaned up name passed to the plots (internal)
#' @param ... Additional arguments passed to zoo package
#'
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom stats loess
#' @importFrom zoo rollapply
#' @import ggplot2
#'
#' @return NULL
#'
#' @examples NULL
Utility_LocalMaxima <- function(theX, theY, therepeats, w, alternatename, ...){

  #Adding Margins
  repeats <- therepeats*2
  NewX <- length(theX) + repeats
  NewX <- 1:NewX
  NewYmin <- theY[[1]]*0.99
  LengthY <- length(theY)
  NewYmax <- theY[[LengthY]]*0.80
  replicatedYmin <- rep(NewYmin, each = therepeats)
  replicatedYmax <- rep(NewYmax, each = therepeats)
  NewY <- c(replicatedYmin, theY, replicatedYmax)

  x <- NewX
  y <- NewY

  ##Local Maxima Functions
  n <- length(y)
  #y.smooth <- loess(y ~ x)$fitted
  y.smooth <- loess(y ~ x, ...)$fitted
  y.max <- rollapply(zoo(y.smooth), 2*w+1, max, align="center")
  delta <- y.max - y.smooth[-c(1:w, n+1-1:w)]
  i.max <- which(delta <= 0) + w
  peaks <- list(x=x[i.max]-therepeats, i=i.max-therepeats, y.hat=y.smooth[(therepeats + 1):(length(y.smooth) - therepeats)])

  peak_points <- peaks$x

  MainData <- data.frame(x = theX, y = theY, yhat = peaks$y.hat)

  PointData <- MainData %>% filter(x %in% peak_points)

  Views <- ggplot(MainData, aes(x = x, y = y)) + geom_point(size = 2, color = "Gray") + geom_line(aes(y = yhat), linewidth = 1) + geom_point(data = PointData, aes(x, yhat), color = "Red", shape = 19, size = 2) + geom_segment(data = PointData, aes(x = x, xend = x, y = 0, yend = yhat), color = "Red", linewidth = 1, linetype = "dashed") + labs(title = alternatename) + theme_bw() +  theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(), axis.title.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(Views)

  PointData <- PointData %>% select(-y)

  return(PointData)
}
