#' @docType data
#' @name wilhelm xml form
#' @title Small bistable Chemical Reaction Network (xml)
#' @examples
#' xml2::read_xml(system.file("extdata", "wilhelm.xml", package = "racipeCRN"))
#' @description This is a small CRN known to be capable of bistability.
#' Since it is capable of multistability under changing parameters and is
#' designed to prevent explosions, it is ideal for examples using the racipeCRN
#' simulation. For further details and the mathematical derivation, see
#' Wilhelm T. (2009). The smallest chemical reaction system with bistability.
#' BMC systems biology, 3, 90. https://doi.org/10.1186/1752-0509-3-90
#' @format xml_document in SBML level 3, version 2 format with 2 species and
#' 4 reactions
NULL

#' @docType data
#' @name wilhelm txt form
#' @title Small bistable Chemical Reaction Network (tpo)
#' @examples
#'  system.file("extdata", "wilhelm.tpo", package = "racipeCRN")
#' @description This is a small CRN known to be capable of bistability. This
#' is useful for testing. For further details and the mathematical derivation,
#' see
#' Wilhelm T. (2009). The smallest chemical reaction system with bistability.
#' BMC systems biology, 3, 90. https://doi.org/10.1186/1752-0509-3-90
#' @format .tpo file with 4 lines representing 2 species with 4 reactions.
#' Please see README for details on structure of .tpo files.
NULL
