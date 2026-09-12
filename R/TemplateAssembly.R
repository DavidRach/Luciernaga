
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