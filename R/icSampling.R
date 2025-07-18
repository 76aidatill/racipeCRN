
#' @export
#' @title Sample from a specified Stoichiometric Compatibility class for a CRN
#' @import SummarizedExperiment
#' @import volesti
#' @importFrom utils read.table write.table data
#' @importFrom S4Vectors metadata
#' @description Given the topology of a CRN, this function samples points from
#' the phase space of the CRN which all belong to the same stoichiometric
#' compatibility class (SCC). The function uses the stoichiometry matrix to
#' distinguish between conservative CRNs, non-conservative CRNs with geometric
#' constraints, and non-conservative CRNs without geometric constraints.
#' For non-conservative CRNs, a bounding parameter is used so the sampling space
#' is finite. Sampling of polytopes is performed using random walks in volesti.
#'
#' @param network xml_document or cRacipeSE object or character containing the
#' file path to a text file with the network's topology. xml_document objects
#' must be in SBML format. If it is a cRacipeSE object, the output will be that
#' object with updated initial conditions. Else, the output will be a dataframe.
#' @param classPoint (optional) vector or matrix. Default NULL. This argument
#' can be used to define the stoichiometric compatibility class. If
#' \code{sameClass} is true, it should be a vector with numSpecies elements.
#' If \code{sameClass} is false, it should be a matrix with numSpecies rows and
#' numModels columns. Every element should be nonnegative. If not provided, the
#' classPoints will be randomly generated.
#' @param numModels (optional) Integer. Default 2000. Number of random models
#' to be simulated.
#' @param numIC (optional) Integer. Default 1. Number of initial conditions to
#' simulate for each model.
#' @param walkLength (optional) Integer. Default 1. Number of steps between
#' generated points in the random walk.
#' @param numBurns (optional) Integer. Default 0. Number of initial points from
#' the random walk to discard before sampling starts.
#' @param lambda (optional) numeric. Default 10. Bounding parameter used in for
#' sampling non-conservative networks to make the sampling space have finite
#' size
#' @param sameClass (optional) logical. Default \code{FALSE}.
#' Whether or not to sample from the same SCC for each model in the ensemble.
#'
#' @return cRacipeSE object or dataframe
#' @examples
#' wilhelm <- system.file("extdata", "wilhelm.tpo", package = "racipeCRN")
#' wilhelmICs <- cracipeSample(network = wilhelm, numModels = 10, lambda = 1)
#'
#' @section Related Functions:
#'
#' \code{\link{cracipeNetwork}}
#' \code{\link{cracipeSimulate}}


cracipeSample <- function(network, classPoint = NULL, numModels = 2000, numIC=1,
                          walkLength = 1, numBurns = 0, lambda = 10.0,
                          sameClass = FALSE){

  cSet <- cRacipeSE()
  metadataTmp <- metadata(cSet)

  if(methods::is(network,"cRacipeSE"))
  {
    cSet <- network
    metadataTmp <- metadata(cSet)
  }
  if((methods::is(network,"character")) | (methods::is(network, "xml_document")))
  {
    cracipeNetwork(cSet) <- network
    metadataTmp$reactantMatrix <- metadata(cSet)$reactantMatrix
  }

  if(missing(network)){
    message("Please specify a network either as a file or as an xml_document")
    return()
  }

  numSpecies <- length(cSet)

  #classPoint validation
  if(!is.null(classPoint)){
    if(any(is.na(classPoint))){
      message("NA values in classPoint, stopping function")
      return()
    }
    if(!(all(classPoint >= 0))){
      message("Negative values in classPoint, stopping function")
    }
    if(!sameClass){
      if(!all(dim(classPoint) == c(numSpecies, numModels))){
        message("Invalid Dimensions for classPoint")
        return()
      }
    }
  }

  stoichMatrix <- as.matrix(rowData(cSet))
  speciesNames <- names(cSet)

  ##Construct basis for stoichiometric subspace
  stoichDecomp <- qr(stoichMatrix)
  stoichSubspace <- qr.Q(stoichDecomp)
  stoichSubspace <- stoichSubspace[, 1:stoichDecomp$rank]


  ##Sample based on case
  if(stoichDecomp$rank == numSpecies){
    message("Sampling full-rank network")
    reactionSimplexV <- volesti::gen_cube(dimension = numSpecies,
                                          representation = 'V')
    uniformVSampling <- volesti::sample_points(P = reactionSimplexV,
                                               n=numSpecies*numIC,
                                               random_walk = list("walk_length" = walkLength, "nburns" = numBurns),
                                               distribution = list("density" = "uniform"))
    #stretch points
    uniformVSampling <- (uniformVSampling + 1)*lambda

  } else{
    message("Sampling network with conservation laws")
    stoichSubspace <- qr.Q(stoichDecomp)
    stoichSubspace <- stoichSubspace[, 1:stoichDecomp$rank]
    hspaceDim <- dim(stoichSubspace)[[2]]

    if(!sameClass){
      classPoints <- matrix(runif(n=numModels*numSpecies, min = 0, max = lambda),
                              nrow = numSpecies, ncol = numModels)
      uniformVSampling <- list()

      for(i in seq_len(numModels)){
        classPoint <- classPoints[,i]
        reactionSimplexH <- volesti::Hpolytope(A = rbind(-stoichSubspace, diag(hspaceDim), -diag(hspaceDim)),
                                               b = c(classPoint, rep(lambda, hspaceDim), rep(lambda, hspaceDim)))
        uniformHSampling <- volesti::sample_points(P = reactionSimplexH, n = numIC,
                                                   random_walk = list("walk_length" = walkLength, "nburns" = numBurns),
                                                   distribution = list("density" = "uniform"))

        uniformVSampling[[i]] <- apply(uniformHSampling, MARGIN = 2, function(x) {stoichSubspace %*% x+ classPoint})
      }

      uniformVSampling <- do.call(cbind, args = uniformVSampling)


    } else{
      if(is.null(classPoint)){
        classPoint <- runif(numSpecies, min = 0, max = lambda)
      }
      reactionSimplexH <- volesti::Hpolytope(A = rbind(-stoichSubspace, diag(hspaceDim), -diag(hspaceDim)),
                                    b = c(classPoint, rep(lambda, hspaceDim), rep(lambda, hspaceDim)))
      uniformHSampling <- volesti::sample_points(P = reactionSimplexH, n = numModels*numIC,
                                                 random_walk = list("walk_length" = walkLength, "nburns" = numBurns),
                                                 distribution = list("density" = "uniform"))
      uniformVSampling <- apply(uniformHSampling, MARGIN = 2, function(x) {stoichSubspace %*% x+ classPoint})
    }
  }
  uniformVSampling <- t(uniformVSampling)
  colnames(uniformVSampling) <- speciesNames
  return(uniformVSampling)


}
