
#' @export
#' @title Simulate a chemical reaction network
#' @import SummarizedExperiment
#' @importFrom utils read.table write.table data
#' @importFrom S4Vectors metadata
#' @description Simulate an ensemble of models for a chemical reaction network
#' using its topology as the only explicit input. The ensemble is randomly
#' generated.
#' @param network xml_document or cRacipeSE object or character containing the
#' file path to a text file with the network's topology. xml_document objects
#' must be in SBML format.
#' @param config (optional) List. Contains the simulation parameters and
#' hyperparameters. It should contain the vectors simParams, hyperParams, and
#' options.
#' @param numModels (optional) Integer. Default 2000. Number of random models
#' to be simulated.
#' @param numIC (optional) Integer. Default 1. Number of initial conditions to
#' simulate for each model.
#' @param outputPrecision (optional) integer. Default 12.
#' The decimal point precision of the output parameters and concentrations.
#' @param kMin (optional) numeric. Default 0. Minimum of uniform distribution
#' for sampling kinetic rate constants.
#' @param kMax (optional) numeric. Default 10. Maximum of uniform distribution
#' for sampling kinetic rate constants.
#' @param icMin (optional) numeric. Default 0. Minimum of uniform distribution
#' for sampling initial conditions.
#' @param icMax (optional) numeric. Default 100. Maximum of uniform distribution
#' for sampling initial conditions.
#' @param stepper (optional) character. Default "E". Defines the integration
#' method used for simulation. Use "E" for the 1st-Order explicit Euler method.
#' @param integrate (optional) logical. Default \code{TRUE}. Whether to
#' integrate the differential equations or not. If \code{FALSE}, the function
#' will only generate the parameters and initial conditions. This can be used
#' for setting the kinetic rate constants and initial conditions and then
#' adjusting them as needed.
#' @param genParams (optional) logical. Default \code{TRUE}. Whether or not to
#' generate kinetic parameters for each model. If set to \code{FALSE}, the
#' network should be a cRacipeSE object with preexisting parameters to use.
#' @param genIC (optional) logical. Default \code{TRUE}. Whether or not to
#' generate initial conditions for each model. If set to \code{FALSE}, the
#' network should be a cRacipeSE object with preexisting initial conditions to
#' use.
#'
#' @return cRacipeSE object
#' @examples
#' wilhelm <- system.file("extdata", "wilhelm.tpo", package = "racipeCRN")
#' cSet <- cracipeSimulate(wilhelm)
#'
#' @section Related Functions:
#'
#' \code{\link{cracipeNetwork}}
#'

cracipeSimulate <- function(network, config = config, numModels = 2000,
                            numIC = 1, outputPrecision = 12,
                            kMin = 0, kMax = 10, icMin = 0, icMax = 100,
                            stepper = "E",
                            integrate = TRUE, genParams = TRUE, genIC = TRUE){
  cSet <- cRacipeSE()
  metadataTmp <- metadata(cSet)
  configData <- NULL
  data("configData",envir = environment(), package = "racipeCRN")
  configuration <- configData

  if(methods::is(network,"cRacipeSE"))
  {
    cSet <- network
    metadataTmp <- metadata(cSet)
    configuration <- metadata(cSet)$config
  }
  if((methods::is(network,"character")) | (methods::is(network, "xml_document")))
  {
    cracipeNetwork(cSet) <- network
    metadataTmp$numReactions <- dim(as.matrix(rowData(cSet)))[2]
    metadataTmp$reactantMatrix <- metadata(cSet)$reactantMatrix
  }

  if(missing(network)){
    message("Please specify a network either as a file or as an xml_document")
    return()
  }
  if(!missing(config)){
    if(methods::is(config,"list")){
      configuration <- config
    }
    else{
      message("Incorrect config provided!")
      return()
    }
  }


  if(!missing(numModels)){
    configuration$simParams["numModels"] <- numModels
  }
  if(!missing(numIC)){
    configuration$simParams["numIC"] <- numIC
  }
  if(!missing(outputPrecision)){
    configuration$simParams["outputPrecision"] <- outputPrecision
  }
  if(!missing(kMin)){
    configuration$hyperParams["kMin"] <- kMin
  }
  if(!missing(kMax)){
    configuration$hyperParams["kMax"] <- kMax
  }
  if(!missing(icMin)){
    configuration$hyperParams["icMin"] <- icMin
  }
  if(!missing(icMax)){
    configuration$hyperParams["icMax"] <- icMax
  }
  if(!missing(integrate)){
    configuration$options["integrate"] <- integrate
  }
  if(!missing(genParams)){
    configuration$options["genParams"] <- genParams
  }
  if(!missing(genIC)){
    configuration$options["genIC"] <- genIC
  }

  if(!genParams){
    if(is.null(colData(cSet))){
      message("Please specify the parameters as genParams is FALSE")
      return(cSet)
    } else{
      params <- as.data.frame(cracipeParams(cSet))

      utils::write.table(params, file = outFileParams,
                         sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    }
  }

  if(!genIC){
    if(is.null(colData(cSet))){
      message("Please specify the initial conditions
              as genIC is FALSE")
      return(cSet)
    } else {
      ic <- as.data.frame(t(cracipeIC(cSet)))
      utils::write.table(ic, file = outFileIC,
                         sep = "\t", quote = FALSE, row.names = FALSE,
                         col.names = FALSE)
    }
  }

  configuration$stepper <- stepper
  stepperInt <- 1L

  numSpecies <- length(cSet)
  stoichMatrix <- as.matrix(rowData(cSet))
  reactantMatrix <- as.matrix(metadataTmp$reactantMatrix)
  speciesNames <- names(cSet)

  outFileSC <- tempfile(fileext = ".txt") #File for species concentration
  outFileParams <- tempfile(fileext = ".txt")
  outFileIC <- tempfile(fileext = ".txt")
  outFileConverge <- tempfile(fileext = ".txt")



  annotationTmp <- outFileSC
  message("Running the simulations")

  timeEvolutionTest <- simulateCRN(stoichMatrix, reactantMatrix, configuration,
                                   outFileSC, outFileParams, outFileIC,
                                   outFileConverge, stepperInt)

  if(configuration$options["integrate"]){
    print("INTEGRATION NOT CURRENTLY AVAILABLE")
    return(cSet)
  } else {
    parameters <- utils::read.table(outFileParams, header = FALSE)

    ic <- utils::read.table(outFileIC, header = FALSE)
    colnames(ic) <- speciesNames

    colData <- (cbind(parameters, ic))
    metadataTmp$config <- configuration

    cSet <- cRacipeSE(rowData = stoichMatrix, colData = colData,

                     metadata = metadataTmp)
    return(cSet)
  }
}
