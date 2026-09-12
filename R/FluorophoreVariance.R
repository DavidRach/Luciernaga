
#' Internal for SimulatedData, generates total variant events and multiplies by signature
#'
#' @param x Iterated in Fluorophore
#' @param LocalNumber Iterated in number of events (derrived from pop and abundance)
#' @param ScaledData The Scaled signature data
#' @param distribution The desired distribution for the population (1, 2, 3)
#'
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#' @importFrom stats rlnorm
#'
#' @return A matrix of detector data specific for that fluorophore by respective events
#'
#' @noRd
FluorophoreVariance <- function(x, LocalNumber, ScaledData, distribution){
LocalData <- ScaledData %>% filter(Fluorophore %in% x) %>% pull(AdjustedY)
LocalDist <- distribution %>% filter(Markers %in% x) %>% pull(Distribution)

if (LocalDist == 1){
mu <- 0 # Mean
sigma <- 0.05 #SD
Values <- rlnorm(LocalNumber, meanlog = mu, sdlog=sigma)
#mean(Values)
#hist(Values, main = "Log-normal Distribution", xlab = "Value", breaks = 50)
#LocalData

ResultMatrix <- sapply(Values, function(value) LocalData * value)
ResultMatrix <- t(ResultMatrix)

#Results <- data.frame(ResultMatrix)
#A <- do.call(pmax, Results)
#Normalized <- Results/A
#View(Normalized)

return(ResultMatrix)
}
}