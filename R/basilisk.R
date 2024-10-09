#' Internal function for Basilisk environment setup
#'
#' @importFrom basilisk BasiliskEnvironment
#'
#' @return An internal value
#'
#' @noRd
env1 <- BasiliskEnvironment("env1", pkgname="Luciernaga",
                            packages=c("pandas==2.2.1", "numpy==1.26.4", "pacmap==0.7.2",
                                       "phate==1.0.11"))

