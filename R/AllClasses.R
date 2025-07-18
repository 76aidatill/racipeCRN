#' @export
#' @import methods
#' @import SummarizedExperiment
#' @title cRacipeSE
#' @description An S4 class for Random Circuit Perturbation (RACIPE)
#' simulations of chemical reaction networks.
#' Extends the \link{SummarizedExperiment} class.
#'
#'
.cRacipeSE <- setClass("cRacipeSE",
                      contains = "SummarizedExperiment")
#' @export
#' @title cRacipeSE constructor
#' @description Create an RacipeSE object. cRacipeSE is an S4 class for
#' Random Circuit Perturbation (RACIPE) simulations of chemical reaction
#' networks in which a large number of models with randomized parameters are
#' used for simulation of the circuit. Each model can be considered as a sample.
#' It extends the \link{SummarizedExperiment} class to store and access
#' the network, simulated expressions, parameters, intial conditions and
#' other meta information.
#' SummarizedExperiment slot assays is used for storing simulated
#' species expressions. The rows of these
#' matrix-like elements correspond to various species in the circuit and columns
#' correspond to models.
#'
#' @import SummarizedExperiment
#' @importFrom S4Vectors metadata DataFrame isEmpty SimpleList
#' @importFrom utils data
#' @param .object (optional) Another cRacipeSE object.
#' @param assays (optional) assay object for initialization
#' @param rowData (optional) rowData for initialization
#' @param colData (optional) colData for initialization
#' @param metadata (optional) metadata for initialization
#' @param ... Arguments passed to SummarizedExperiment
#' @return cRacipeSE object
#' @examples
#' cSet <- cRacipeSE()
#'

cRacipeSE <- function(.object = NULL, assays = SimpleList(),
                     rowData = NULL,
                     colData = DataFrame(),
                     metadata = list(), ...) {

  if(is(.object,"cRacipeSE")) {
    if(isEmpty(assays)) assays <- assays(.object)
    if(is.null(rowData)) rowData <- rowData(.object)
    if(isEmpty(colData)) colData <- colData(.object)
    if(length(metadata)==0) metadata <- metadata(.object)
  }
  objectTmp <- SummarizedExperiment(
    assays=assays, rowData=rowData, colData=colData, metadata=metadata)

  if(is.null(metadata(objectTmp)$config))
  {
    configData <- NULL
    data("configData",envir = environment(), package = "racipeCRN")
    metadata(objectTmp)$config <- configData
  }

  .cRacipeSE(objectTmp)
}
