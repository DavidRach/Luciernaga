#' @importFrom basilisk BasiliskEnvironment
env1 <- BasiliskEnvironment("env1", pkgname="Luciernaga",
                            packages=c("pandas==2.2.1", "numpy==1.26.4", "pacmap==0.7.2"))

env2 <- BasiliskEnvironment("env2", pkgname="Luciernaga",
                            packages=c("pandas==2.2.1", "numpy==1.26.4", "phate==1.0.11"))

