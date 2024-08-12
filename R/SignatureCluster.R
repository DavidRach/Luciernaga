#' Internal for LuciernagaQC SingleStainSignatures
#'
#' @param Arg1 The first Detector Column
#' @param Arg2 The second Detector Column
#' @param data A data.frame
#'
#' @noRd

SignatureCluster <- function(Arg1, Arg2, data){
  data <- data %>% mutate(Cluster = case_when(
    near(data[[Arg1]], 0.0) ~ paste(data$Cluster, Arg1, "_00-", sep = "", collapse = NULL),
    near(data[[Arg1]], 0.1) ~ paste(data$Cluster, Arg1, "_01-", sep = "", collapse = NULL),
    near(data[[Arg1]], 0.2) ~ paste(data$Cluster, Arg1, "_02-", sep = "", collapse = NULL),
    near(data[[Arg1]], 0.3) ~ paste(data$Cluster, Arg1, "_03-", sep = "", collapse = NULL),
    near(data[[Arg1]], 0.4) ~ paste(data$Cluster, Arg1, "_04-", sep = "", collapse = NULL),
    near(data[[Arg1]], 0.5) ~ paste(data$Cluster, Arg1, "_05-", sep = "", collapse = NULL),
    near(data[[Arg1]], 0.6) ~ paste(data$Cluster, Arg1, "_06-", sep = "", collapse = NULL),
    near(data[[Arg1]], 0.7) ~ paste(data$Cluster, Arg1, "_07-", sep = "", collapse = NULL),
    near(data[[Arg1]], 0.8) ~ paste(data$Cluster, Arg1, "_08-", sep = "", collapse = NULL),
    near(data[[Arg1]], 0.9) ~ paste(data$Cluster, Arg1, "_09-", sep = "", collapse = NULL),
    near(data[[Arg1]], 1.0) ~ paste(data$Cluster, Arg1, "_10-", sep = "", collapse = NULL)))

  Second <- data %>% mutate(Cluster = case_when(
    near(data[[Arg2]], 0.0) ~ paste(data$Cluster, Arg2, "_00-", sep = "", collapse = NULL),
    near(data[[Arg2]], 0.1) ~ paste(data$Cluster, Arg2, "_01-", sep = "", collapse = NULL),
    near(data[[Arg2]], 0.2) ~ paste(data$Cluster, Arg2, "_02-", sep = "", collapse = NULL),
    near(data[[Arg2]], 0.3) ~ paste(data$Cluster, Arg2, "_03-", sep = "", collapse = NULL),
    near(data[[Arg2]], 0.4) ~ paste(data$Cluster, Arg2, "_04-", sep = "", collapse = NULL),
    near(data[[Arg2]], 0.5) ~ paste(data$Cluster, Arg2, "_05-", sep = "", collapse = NULL),
    near(data[[Arg2]], 0.6) ~ paste(data$Cluster, Arg2, "_06-", sep = "", collapse = NULL),
    near(data[[Arg2]], 0.7) ~ paste(data$Cluster, Arg2, "_07-", sep = "", collapse = NULL),
    near(data[[Arg2]], 0.8) ~ paste(data$Cluster, Arg2, "_08-", sep = "", collapse = NULL),
    near(data[[Arg2]], 0.9) ~ paste(data$Cluster, Arg2, "_09-", sep = "", collapse = NULL),
    near(data[[Arg2]], 1.0) ~ paste(data$Cluster, Arg2, "_10-", sep = "", collapse = NULL)))

  return(Second)
}
