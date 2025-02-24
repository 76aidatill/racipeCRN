#' @export
#' @title Generate parameter names for a network
#' @param network cRacipeSE object
#' @examples
#' cSet <- cRacipeSE()
#' wilhelm <- xml2::read_xml(system.file("extdata", "wilhelm.xml", package = "racipeCRN"))
#' cracipeNetwork(cSet) <- wilhelm
#' paramNames <- cracipeGenParamNames(cSet)
#'
#' @return character vector
cracipeGenParamNames <- function(network){
  if (!(methods::is(network, "cRacipeSE"))){
    stop("Function only accepts cRacipeSE objects")
  }
  stoichMatrix <- as.matrix(rowData(network))
  reacNames <- colnames(stoichMatrix)

  paramNames <- paste0("k_", reacNames)
  return(paramNames)
}
