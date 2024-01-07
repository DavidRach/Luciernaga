#' Forming a pdf to plot Luciernaga returns
#'
#' @param x A list of plots
#' @param filename Name of the generated .pdf
#' @param outfolder Location to store the .pdf
#'
#' @return NULL
#' @export
#'
#' @examples NULL
Utility_SCPatchwork <- function(x, filename, outfolder){
  theList <- x
  #theList <- Filter(Negate(is.null), theList)
  theListLength <- length(theList)

  theoreticalitems <- 3
  #theoreticalitems <- 4
  #theoreticalpages <- theListLength/theoreticalitems
  #theoreticalpages <- round(theoreticalpages, 0)

  split_list <- function(input_list, chunk_size) {
    split(input_list, ceiling(seq_along(input_list) / chunk_size))
  }

  sublists <- split_list(theList, theoreticalitems)
  length(sublists)

  MergedName <- paste(outfolder, filename, sep = "/")

  pdf(file = paste(MergedName, ".pdf", sep = "", collapse = NULL), width = 11, height = 8.5) #Optional Adjustments for Second
  #pdf(file = paste(MergedName, ".pdf", sep = "", collapse = NULL), width = 8.5, height = 11) #Optional Adjustments for Second

  spacer <- plot_spacer()

  for(i in sublists){
    p1 <- i[[1]]
    p2 <- i[[2]]
    p3 <- i[[3]]
    #p4 <- i[[4]]

    arranged_plot <- p1 / (p2 + p3)
    #arranged_plot <- p1 / (p2 + p3) / p4

    p <- wrap_plots(list(arranged_plot), widths = 0.8, heights = 0.8)
    print(p)
  }

  dev.off()

}
