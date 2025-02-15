
#' @export
#' @import SummarizedExperiment
#' @title Method to get the CRN Stoichiometry Matrix and Rate Vector
#' @description This method will return the Stoichiometry matrix and the
#' Reactant coefficient matrix in species by reaction form.
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

#' @export
#' @import SummarizedExperiment
#' @importFrom S4Vectors metadata
#' @title A method to access the kinetic rate constants for each model in the
#' ensemble
#' @description Gets the kinetic rate constants for each model in the form of
#' a data.frame.
#' @param .object cRacipeSE object
#' @examples
#' wilhelm <- xml2::read_xml(system.file("extdata", "wilhelm.xml", package = "racipeCRN"))
#' cSet <- cracipeSimulate(wilhelm, integrate = FALSE, numModels=20)
#' parameters <- cracipeParams(cSet)
#' rm(parameters,cSet)
#' @return A  data.frame with numReactions columns and numModels rows.
#'
setGeneric("cracipeParams",
           def = function(.object)
           {
             standardGeneric("cracipeParams")
           }
)

#' @export
#' @import SummarizedExperiment
#' @title  A method to get the initial conditions used for simulations
#' @description The initial conditions of each of the models.
#' @param .object cRacipeSE object
#' @examples
#' wilhelm <- xml2::read_xml(system.file("extdata", "wilhelm.xml", package = "racipeCRN"))
#' cSet <- cracipeSimulate(wilhelm, integrate = FALSE, numModels=20)
#' ic <- cracipeParams(cSet)
#' rm(ic,cSet)
#' @return DataFrame with numModels columns and numSpecies rows.

setGeneric("cracipeIC",
           def = function(.object)
           {
             standardGeneric("cracipeIC")
           }
)
