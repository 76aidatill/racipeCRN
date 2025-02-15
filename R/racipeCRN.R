# File for generating NAMESPACE and package info using roxygen

#' racipCRN: A package for the random circuit perturbation of chemical reaction
#' networks.
#'
#' racipeCRN is a tool for analyzing chemical reaction networks (CRNs) under
#' changing kinetic parameters. With just the CRN topology, racipeCRN
#' generates an ensemble of models for the CRN using uniform random sampling
#' of the rate constant parameters.
#'
#' @import Rcpp
#' @useDynLib racipeCRN, .registration = TRUE
#' @docType package
#' @name racipeCRN
NULL
