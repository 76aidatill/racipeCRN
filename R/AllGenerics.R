
#' @export
#' @import SummarizedExperiment
#' @title Method to get the complete CRN in matrix form
#' @description This method will return the Stoichiometry matrix and the
#' Reactant coefficient matrix in species by reaction form. These two matrices
#' completely define the CRN with mass-action kinetics
#' @param .object cRacipeSE object
#' @return list of dataframes
#' @examples
#' cSet <- cRacipeSE()
#' wilhelm <- xml2::read_xml(system.file("extdata", "wilhelm.xml", package = "racipeCRN"))
#' cracipeNetwork(cSet) <- wilhelm
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
#' cSet <- cracipeSimulate(wilhelm, integrate = FALSE, numModels=10)
#' parameters <- cracipeParams(cSet)
#' rm(parameters,cSet)
#' @return A data.frame with numReactions columns and numModels rows.
#'
setGeneric("cracipeParams",
           def = function(.object)
           {
             standardGeneric("cracipeParams")
           }
)

#' @export
#' @title  A method to set the ensemble kinetic parameters
#' @description Set the parameters
#' @param .object cRacipeSE object
#' @param value DataFrame containing the parameters. Dimensions should be
#' (numModels) rows by (numReactions) columns.
#' @examples
#' wilhelm <- xml2::read_xml(system.file("extdata", "wilhelm.xml", package = "racipeCRN"))
#' cSet <- racipeCRN::cracipeSimulate(network = wilhelm, numModels = 10,
#' integrate = FALSE)
#' parameters <- cracipeParams(cSet)
#' cracipeParams(cSet) <- parameters
#' rm(parameters, cSet)
#' @return A cRacipeSE object
#'

setGeneric("cracipeParams<-",
           def = function(.object, value)
           {
             standardGeneric("cracipeParams<-")
           }
)


#' @export
#' @import SummarizedExperiment
#' @title  A method to get the initial conditions used for simulations
#' @description The initial conditions of each of the models.
#' @param .object cRacipeSE object
#' @examples
#' wilhelm <- xml2::read_xml(system.file("extdata", "wilhelm.xml", package = "racipeCRN"))
#' cSet <- cracipeSimulate(wilhelm, integrate = FALSE, numModels=10)
#' ic <- cracipeParams(cSet)
#' rm(ic,cSet)
#' @return DataFrame with numModels columns and numSpecies rows.

setGeneric("cracipeIC",
           def = function(.object)
           {
             standardGeneric("cracipeIC")
           }
)

#' @export
#' @import SummarizedExperiment
#' @title  A method to set the initial conditions
#' @description Set the initial conditions
#' @param .object cRacipeSE object
#' @param value DataFrame containing the initial conditions, dimensions should
#' be (# species) rows by (numModels*numIC) columns
#' @examples
#' wilhelm <- xml2::read_xml(system.file("extdata", "wilhelm.xml", package = "racipeCRN"))
#' cSet <- racipeCRN::cracipeSimulate(network = wilhelm, numModels = 10,
#' integrate=FALSE)
#' ics <- cracipeIC(cSet)
#' cracipeIC(cSet) <- ics
#' rm(cSet, ics)
#' @return A cRacipeSE object

setGeneric("cracipeIC<-",
           def = function(.object, value)
           {
             standardGeneric("cracipeIC<-")
           }
)

#' @export
#' @import SummarizedExperiment
#' @title  A method to access the simulation hyperparameters
#' @description The hyperparameters like number of models, range from which
#' parameters are to be sampled, simulation time etc.
#' @param .object cRacipeSE object
#' @examples
#' cSet <- cRacipeSE()
#' wilhelm <- xml2::read_xml(system.file("extdata", "wilhelm.xml", package = "racipeCRN"))
#' cracipeNetwork(cSet) <- wilhelm
#' cracipeConfig(cSet)
#' rm(cSet)
#' @return list
#'

setGeneric(name="cracipeConfig",
           def=function(.object)
           {
             standardGeneric("cracipeConfig")
           }
)

#' @export
#' @import SummarizedExperiment
#' @title  A method to get the convergence results
#' @description Gathers the convergence and speed of convergence for each model
#' and initial condition.
#' @param .object cRacipeSE object

#' @examples
#' wilhelm <- xml2::read_xml(system.file("extdata", "wilhelm.xml", package = "racipeCRN"))
#' cSet <- racipeCRN::cracipeSimulate(network = wilhelm, numModels = 100, numIC = 3)
#' cd <- cracipeConverge(cSet)
#' rm(cSet, cd)
#' @return DataFrame

setGeneric("cracipeConverge",
           def = function(.object)
           {
             standardGeneric("cracipeConverge")
           }
)


#' @export
#' @import SummarizedExperiment
#' @title  A method for grabbing unique states
#' @description This method selects the unique expression states for every model
#' in a racipeCRN object. Non-converged expressions are  filtered out using the
#' simulation convergence data. To minimize error from numerical integration,
#' concentration values are rounded before comparison.
#' @param .object cRacipeSE object generated by \code{\link{cracipeSimulate}}
#'  function.
#' @param roundDigits (optional) integer. Default \code{3}. Number of digits for
#' concentration rounding.
#' @examples
#' wilhelm <- xml2::read_xml(system.file("extdata", "wilhelm.xml", package = "racipeCRN"))
#' \dontrun{
#' cSet <- racipeCRN::cracipeSimulate(network = wilhelm, numModels = 20,
#' integrateStepSize = 0.001, numIterations = 30)
#' stateList <- cracipeUniqueStates(cSet)
#' combinedStates <- do.call(cbind, stateList)
#' }
#' @return \code{list} object. Element i of the list is a data frame containing
#' the unique expressions of model i in the input cRacipeSE object
#'@section Related Functions:
#'
#' \code{\link{cracipeSimulate}}
#'
setGeneric("cracipeUniqueStates",
           def = function(.object,
                          roundDigits = 3)
           {
             standardGeneric("cracipeUniqueStates")
           }
)

#' @export
#' @import SummarizedExperiment
#' @importFrom S4Vectors SimpleList
#' @importFrom stats sd
#' @title  Normalize the simulated concentrations
#' @description Log2 normalize the simulated concentrations and center them
#' about 0.
#' @param .object cRacipeSE object
#' @return A cRacipeSE object
#' @examples
#' wilhelm <- xml2::read_xml(system.file("extdata", "wilhelm.xml", package = "racipeCRN"))
#' cSet <- racipeCRN::cracipeSimulate(network = wilhelm, numModels = 20,
#' integrateStepSize = 0.01, numIterations = 30)
#' cSet <- cracipeNormalize(cSet)
#' @section Related Functions:
#'
#' \code{\link{cracipeSimulate}}, \code{\link{cracipePlotData}}
#'
setGeneric("cracipeNormalize",
           def = function(.object)
           {
             standardGeneric("cracipeNormalize")
           }
)

#' @export
#' @import SummarizedExperiment
#' @importFrom gplots heatmap.2
#' @importFrom graphics barplot hist image layout par
#' @importFrom stats as.dendrogram prcomp
#' @import ggplot2
#' @import gridExtra
#' @import umap
#' @import grDevices
#' @import RColorBrewer
#' @title Plot sRACIPE data
#' @description Plots heatmap, pca, umap of the data simulated using racipeCRN.
#' Only considers unique steady states for each model.
#' @param  .object List A list returned by \code{\link{cracipeSimulate}} function
#' @param plotToFile (optional) logical. Default \code{FALSE}. Whether to save
#' plots to a file.
#' @param nClusters (optional) Integer. Default 2. Expected number of clusters
#' in the simulated data. Hierarchical clustering will be used to cluster the
#' data and the the models will be colored in UMAP and PCA plots according to
#' these clustering results. The clusters can be also supplied using
#' \code{assignedClusters}.
#' @param heatmapPlot (optional) logical. Default \code{TRUE}. Whether to plot
#' hierarchichal clustering.
#' @param pcaPlot (optional) logical. Default \code{TRUE}. Whether to plot PCA
#' embedding.
#' @param umapPlot (optional) logical. Default \code{TRUE}. Whether to plot
#' UMAP embedding
#' @param clustMethod (optional) character. Default \code{"ward.D2"}. Clustering
#' method for heatmap. See \code{\link[gplots]{heatmap.2}}
#' @param col (optional) Color palette
#' @param distType (optional) Distance type.  Used only if specified
#' explicitly. Otherwise, 1-cor is used. See \code{\link[stats]{dist}},
#' \code{\link[stats]{hclust}}
#' @param assignedClusters vector integer or character. Default NULL.
#' Cluster assignment of models.
#' @param corMethod (optional) character. Default \code{"spearman"}. Correlation
#' method for distance function.
#' @param ... Other arguments
#' @examples
#' wilhelm <- xml2::read_xml(system.file("extdata", "wilhelm.xml", package = "racipeCRN"))
#' \dontrun{
#' cSet <- racipeCRN::cracipeSimulate(network = wilhelm, numModels = 20,
#' integrateStepSize = 0.001, numIterations = 30)
#' cSet <- cracipePlotData(cSet)
#' }
#' @return \code{cRacipeSE} object
#' @section Related Functions:
#'
#' \code{\link{cracipeSimulate}},  \code{\link{cracipeUniqueStates}}
setGeneric("cracipePlotData",
           def = function(.object, plotToFile = FALSE, nClusters = 2,
                          heatmapPlot = TRUE,
                          pcaPlot = TRUE, umapPlot = TRUE,
                          clustMethod = "ward.D2", col = col,
                          distType = "euclidean",
                          assignedClusters = NULL, corMethod = "spearman", ...)
           {
             standardGeneric("cracipePlotData")
           }
)

#' @export
#' @import SummarizedExperiment
#' @importFrom graphics barplot hist image layout par polygon
#' @import ggplot2
#' @title  A method to visualize convergence distributions
#' @description
#' Generates a time plot of the ratio of converged models as the simulation went
#' on. These visualizations can be used for evaluating the ideal simulation run
#' time, as defined by the number of iterations, number of steps per iteration,
#' and the step size. The final proportion of converged concentrations at the
#' end of the simulation is also computed and placed in the output metadata. If
#' the final proportion is greater than 99%, the earliest iteration with 99%
#' convergence is also reported to the metadata.
#' @param .object cRacipeSE object generated by \code{\link{cracipeSimulate}}
#'  function.
#' @param plotToFile (optional) logical. Default \code{FALSE}. Whether to save
#'  plots to a file.
#' @examples
#' wilhelm <- xml2::read_xml(system.file("extdata", "wilhelm.xml", package = "racipeCRN"))
#' \dontrun{
#' cSet <- racipeCRN::cracipeSimulate(network = wilhelm, numModels = 20,
#' integrateStepSize = 0.001, numIterations = 30)
#' cSet <- cracipeConvergeDist(cSet)
#' }
#' @return \code{cRacipeSE} object
#'@section Related Functions:
#'
#' \code{\link{cracipeSimulate}},  \code{\link{cracipePlotData}}
#'
setGeneric("cracipeConvergeDist",
           def = function(.object, plotToFile = FALSE)
           {
             standardGeneric("cracipeConvergeDist")
           }
)
