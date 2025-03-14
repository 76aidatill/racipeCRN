% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data-raw.R
\docType{data}
\name{wilhelm txt form}
\alias{wilhelm txt form}
\title{Small bistable Chemical Reaction Network (tpo)}
\format{
.tpo file with 4 lines representing 2 species with 4 reactions.
Please see README for details on structure of .tpo files.
}
\description{
This is a small CRN known to be capable of bistability. This
is useful for testing. For further details and the mathematical derivation,
see
Wilhelm T. (2009). The smallest chemical reaction system with bistability.
BMC systems biology, 3, 90. https://doi.org/10.1186/1752-0509-3-90
}
\examples{
system.file("extdata", "wilhelm.tpo", package = "racipeCRN")
}
