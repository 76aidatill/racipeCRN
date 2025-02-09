
#' @export
#' @import SummarizedExperiment
#' @title Method to get the CRN Stoichiometry Matrix and Rate Vector
#' @description This method will return the Stoichiometry matrix in species by
#' reaction form and the rate vector as a list where the element for reaction
#' "R1" is a character vector of the reactants in reaction "R1", including
#' repeats.
#' @param .object RacipeSE object
#' @return list of dataframe and list
#'  @examples
#' cSet <- cRacipeSE()
#' wilhelm <- xml2::read_xml(system.file("extdata", "wilhelm.xml", package = "racipeCRN"))
#' cracipeNetwork(rs) <- wilhelm
#' networkTopo <- cracipeNetwork(cSet)
#' rm(cSet, wilhelm,networkTopo)
#'

setGeneric(name="cracipeNetwork",
           def=function(.object)
           {
             standardGeneric("cracipeNetwork")
           }
)

#' @export
#' @import xml2
#' @title Initialize the Network
#' @description
#' Initialize a cRacipeSE object using an SBML file containing a chemical
#' reaction network (CRN). All reactions are assumed to have mass-action
#' kinetics, and any defined rate law in the SBML file will be ignored. Please
#' see \link{makeNetwork} and \link{addReaction} for methods to make a
#' compatible SBML file
#' @param .object cRacipeSE object.
#' @param value xml_document or character. SBML file or .tpo text file containing
#' the network topology. If it is a character, it should be a path to the desired
#' txt file. See README for more information about the format of .tpo files.
#' @returns cRacipeSE
#'
#' @examples
#' cSet <- cRacipeSE()
#' wilhelm <- xml2::read_xml(system.file("extdata", "wilhelm.xml", package = "racipeCRN"))
#' cracipeNetwork(cSet) <- wilhelm
#' rm(cSet, wilhelm)
setGeneric("cracipeNetwork<-",
           def = function(.object, value)
           {
             standardGeneric("cracipeNetwork<-")
           }
)
