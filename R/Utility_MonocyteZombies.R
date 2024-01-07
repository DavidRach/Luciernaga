#' Adds additional gates for monocytes, dead cells, or their respective unstained negatives
#'
#' @param x A Gating Set object (ex. gs or gs[[1]])
#' @param sample.name Keyword for which sample name is stored
#' @param bins Number of bins resulting ggplot2 object should be plotted with
#'
#' @return NULL
#' @export
#'
#' @examples NULL
Utility_MonocyteZombies <- function(x, sample.name, bins){
  x <- x
  name <- keyword(x, sample.name)

  #Unstained
  if (str_detect(name, "Zombie")){
    testthis <- gs_get_leaf_nodes(x)
    if(!str_detect(testthis, "Dead")){gs_add_gating_method(x, parent = "singletsSSCB", pop = "+", alias = "Dead", gating_method = "flowClust", gating_args = "K=2, target=c(5e5,3e5), quantile = 0.5", dims = "FSC-A,SSC-A")}
    cells <- "Dead"
    SingleColor <- as.ggplot(ggcyto(x, aes(x= "FSC-A", y = "SSC-A", ), subset = "singletsSSCB") + geom_hex(bins=bins) + geom_gate(cells) + theme_bw() + labs(title = NULL) + theme(strip.background = element_blank(), strip.text.x = element_blank(), panel.grid.major = element_line(linetype = "blank"), panel.grid.minor = element_line(linetype = "blank"), axis.title = element_text(size = 10, face = "bold"), legend.position = "none"))
    print(SingleColor)
  } else if (str_detect(name, "CD14")){
    testthis <- gs_get_leaf_nodes(x)
    if(!any(str_detect(testthis, "Monocytes"))){gs_add_gating_method(x, parent = "nonDebris", pop = "+", alias = "Monocytes", gating_method = "flowClust", gating_args = "K=2, target=c(2e6,2e6), quantile = 0.8", dims = "FSC-A,SSC-A")}
    cells <- "Monocytes"
    SingleColor1 <- as.ggplot(ggcyto(x, aes(x = "FSC-A", y = "SSC-A"), subset = "nonDebris") + geom_hex(bins=bins) + geom_gate(cells) + theme_bw() + labs(title = NULL) + theme(strip.background = element_blank(), strip.text.x = element_blank(), panel.grid.major = element_line(linetype = "blank"), panel.grid.minor = element_line(linetype = "blank"), axis.title = element_text(size = 10, face = "bold"), legend.position = "none"))
    print(SingleColor1)
  } else if (str_detect(name, "Unstained")){
    testthis <- gs_get_leaf_nodes(x)
    if(!str_detect(testthis, "Dead")){gs_add_gating_method(x, parent = "singletsSSCB", pop = "+", alias = "Dead", gating_method = "flowClust", gating_args = "K=2, target=c(5e5,3e5), quantile = 0.5", dims = "FSC-A,SSC-A")}
    cells <- "Dead"
    SingleColor2 <- as.ggplot(ggcyto(x, aes(x= "FSC-A", y = "SSC-A", ), subset = "singletsSSCB") + geom_hex(bins=bins) + geom_gate(cells) + theme_bw() + labs(title = NULL) + theme(strip.background = element_blank(), strip.text.x = element_blank(), panel.grid.major = element_line(linetype = "blank"), panel.grid.minor = element_line(linetype = "blank"), axis.title = element_text(size = 10, face = "bold"), legend.position = "none"))

    if(!any(str_detect(testthis, "Monocytes"))){gs_add_gating_method(x, parent = "nonDebris", pop = "+", alias = "Monocytes", gating_method = "flowClust", gating_args = "K=2, target=c(2e6,2e6), quantile = 0.8", dims = "FSC-A,SSC-A")}
    cells <- "Monocytes"
    SingleColor3 <- as.ggplot(ggcyto(x, aes(x = "FSC-A", y = "SSC-A"), subset = "nonDebris") + geom_hex(bins=bins) + geom_gate(cells) + theme_bw() + labs(title = NULL) + theme(strip.background = element_blank(), strip.text.x = element_blank(), panel.grid.major = element_line(linetype = "blank"), panel.grid.minor = element_line(linetype = "blank"), axis.title = element_text(size = 10, face = "bold"), legend.position = "none"))
    plot <- SingleColor2 + SingleColor3
    print(plot)
  }
}
