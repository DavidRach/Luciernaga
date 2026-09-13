#' Internal for Simulated Data, generates scatter data and adds to data.frame
#'
#' @param x The existing data.frame of just detectors
#'
#' @importFrom stats rnorm
#' @importFrom stringr str_detect
#' @importFrom dplyr relocate
#'
#' @return The updated data.frame now including time and scatter parameters.
#' @noRd
ScatterParty <- function(x) {

  values <- 0:500000
  Time <- sample(values, size = nrow(x), replace = FALSE)

  SSCW <- data.frame("SSC-W" = rnorm(nrow(x), mean = 1000000, sd = 125000),
                      check.names = FALSE)
  SSCH <- data.frame("SSC-H" = rnorm(nrow(x), mean = 1000000, sd = 125000),
                      check.names = FALSE)
  SSCA <- data.frame("SSC-A" = rnorm(nrow(x), mean = 1000000, sd = 125000),
                      check.names = FALSE)

  FSCW <- data.frame("FSC-W" = rnorm(nrow(x), mean = 1000000, sd = 125000),
                      check.names = FALSE)
  FSCH <- data.frame("FSC-H" = rnorm(nrow(x), mean = 1000000, sd = 125000),
                      check.names = FALSE)
  FSCA <- data.frame("FSC-A" = rnorm(nrow(x), mean = 1000000, sd = 125000),
                      check.names = FALSE)

  SSCBW <- data.frame("SSC-B-W" = rnorm(nrow(x), mean = 1000000, sd = 125000),
                       check.names = FALSE)
  SSCBH <- data.frame("SSC-B-H" = rnorm(nrow(x), mean = 1000000, sd = 125000),
                       check.names = FALSE)
  SSCBA <- data.frame("SSC-B-A" = rnorm(nrow(x), mean = 1000000, sd = 125000),
                       check.names = FALSE)

  TheDetectors <- colnames(x)

  Data <- cbind(Time, x)

  if (any(str_detect(TheDetectors, "UV"))) {
    Data <- cbind(Data, SSCW, SSCH, SSCA) |>
      relocate("SSC-W", "SSC-H", "SSC-A", .before = "V1-A")
  } else {
    Data <- cbind(Data, SSCW, SSCH, SSCA) |>
      relocate("SSC-W", "SSC-H", "SSC-A", .after = "Time")
  }

  Data <- cbind(Data, FSCW, FSCH, FSCA, SSCBW, SSCBH, SSCBA) |>
    relocate("FSC-W", "FSC-H", "FSC-A", "SSC-B-W", "SSC-B-H", "SSC-B-A",
             .before = "B1-A")

  return(Data)
}