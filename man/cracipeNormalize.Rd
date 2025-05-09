% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/methods.R
\name{cracipeNormalize}
\alias{cracipeNormalize}
\alias{cracipeNormalize,cRacipeSE-method}
\title{Normalize the simulated concentrations}
\usage{
cracipeNormalize(.object)

\S4method{cracipeNormalize}{cRacipeSE}(.object)
}
\arguments{
\item{.object}{cRacipeSE object}
}
\value{
A cRacipeSE object
}
\description{
Log2 normalize the simulated concentrations and center them
about 0.
}
\section{Related Functions}{


\code{\link{cracipeSimulate}}, \code{\link{cracipePlotData}}
}

\examples{
wilhelm <- xml2::read_xml(system.file("extdata", "wilhelm.xml", package = "racipeCRN"))
cSet <- racipeCRN::cracipeSimulate(network = wilhelm, numModels = 20,
integrateStepSize = 0.01, numIterations = 30)
cSet <- cracipeNormalize(cSet)
}
