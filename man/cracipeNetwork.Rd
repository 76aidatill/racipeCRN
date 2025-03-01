% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/methods.R
\name{cracipeNetwork}
\alias{cracipeNetwork}
\alias{cracipeNetwork,cRacipeSE-method}
\title{Method to get the CRN Stoichiometry Matrix and Rate Vector}
\usage{
cracipeNetwork(.object)

\S4method{cracipeNetwork}{cRacipeSE}(.object)
}
\arguments{
\item{.object}{RacipeSE object}
}
\value{
list of dataframe and list
@examples
cSet <- cRacipeSE()
wilhelm <- xml2::read_xml(system.file("extdata", "wilhelm.xml", package = "racipeCRN"))
cracipeNetwork(rs) <- wilhelm
networkTopo <- cracipeNetwork(cSet)
rm(cSet, wilhelm,networkTopo)
}
\description{
This method will return the Stoichiometry matrix and the
Reactant coefficient matrix in species by reaction form.
}
