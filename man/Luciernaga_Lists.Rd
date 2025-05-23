% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Luciernaga_Lists.R
\name{Luciernaga_Lists}
\alias{Luciernaga_Lists}
\title{Generate Luciernaga plot outputs from List of List}
\usage{
Luciernaga_Lists(
  ListOfList,
  SecondaryList = NULL,
  PlotType,
  thecolumns = 2,
  therows = 3,
  width = 7,
  height = 9,
  ReturnFolder,
  CurrentExperiment
)
}
\arguments{
\item{ListOfList}{A list containing the returns from Luciernaga_Plots plot option.}

\item{SecondaryList}{Default NULL, if given a list will incorporate the plot in.}

\item{PlotType}{Whether "html" or "pdf"}

\item{thecolumns}{The number of columns per page}

\item{therows}{The number of rows per page}

\item{width}{Desired page width}

\item{height}{Desired page height}

\item{ReturnFolder}{Location to store file}

\item{CurrentExperiment}{Name of Current Experiment}
}
\value{
A file containing the bound Luciernaga plots in specified format
}
\description{
Generate Luciernaga plot outputs from List of List
}
\examples{
NULL
}
