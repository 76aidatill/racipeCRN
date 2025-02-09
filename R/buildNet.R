

#' @export
#' @import xml2
#' @title A method to quickly build an SBML frame
#' @description
#' A method to quickly make an SBML level 3, version 2, List object using xml2.
#' It builds the model and adds a singular compartment for the purpose of adding
#' species. It uses the minimum required parameters and it is designed solely for
#' building chemical reaction networks with mass action kinetics in a single
#' compartment.
#'
#' @param netName character. Name of the newly made chemical reaction network
#'
#' @returns xml_document
#'
#' @examples
#' crn <- makeNetwork("example")
#' rm(crn)
#' @section Related Functions:
#'
#' \code{\link{addReaction}}
#'
makeNetwork <- function(netName){
  # Make new xml file
  network <- xml2::xml_new_root("sbml", version = "2", level = "3", xmlns = "http://www.sbml.org/sbml/level3/version2/core")
  model <- xml2::xml_add_child(network, "model", id = netName)

  # Define compartment, which is necessary for defining species
  compartment_list <- xml2::xml_add_child(model, "listOfCompartments")
  xml2::xml_add_child(compartment_list, "compartment", id = "cell", constant = "true")

  return(network)
}



#' @export
#' @import xml2
#' @title For adding reactions to a pre-existing network
#' @description
#' This function adds reactions to a pre-existing chemical reaction network as
#' defined by the \code{makeNetwork} function. The reactions are added in a way
#' that keeps the model xml object completely compatible with System Biology
#' Markup Language (SBML) Level 3, Version 2. Kinetic laws are not included in
#' SBML file as in crnRACIPE, all reactions are assumed to follow mass-action
#' kinetics. The reaction must include at least one product or reactant.
#' @param network xml_document. An SBML file which defines the chemical reaction
#' network. It should have been originally created with \code{makeNetwork}.
#' @param reactants Named Integer vector. Default \code{NULL}.
#' Each element specifies the stoichiometric coefficient of a reactant in the
#' reaction, and the name of the element specifies the name of the reactant
#' species. If not included, the reaction is assumed to have no reactants.
#' @param products Named Integer vector. Default \code{NULL}. Each element
#' specifies the stoichiometric coefficient of a product in the reaction, and
#' the name of the element specifies the name of the product species. If not
#' included, the reaction is presumed to have no products.
#' @param reversible (optional) Logical. Default FALSE. Whether or not the input
#' reaction is reversible.
#' @param reactionName (optional) Character. Default NULL. The ID of the
#' reaction for the SBML file. If not provided, the SBML ID for the reaction
#' will be "Rn+1", where n is the number of reactions already in the network
#' before this was one added.
#'
#' @returns xml_document. The input network modified to contain the input
#' reaction.
#'
#' @examples
#' example <- makeNetwork("example")
#' example <- addReaction(network = example, reactants = c(A = 1, B = 1),
#'                        products = c(AB = 1), reversible = TRUE,
#'                        reactionName = "Fusion")
#' rm(example)
#' @section Related Functions:
#'
#' \code{\link{makeNetwork}}
#'
addReaction <- function(network, reactants = NULL, products = NULL,
                        reversible = FALSE, reactionName = NULL){
  # Check if listOfSpecies and listOfReactions exist, add if not
  model_node <- xml_find_first(network, "./model")
  if(is.na(xml_find_first(model_node, "./listOfSpecies"))){
    listOfSpecies <- xml_add_child(model_node, "listOfSpecies")
  } else{
    listOfSpecies <- xml_find_first(model_node, "./listOfSpecies")
  }
  if(is.na(xml_find_first(model_node, "./listOfReactions"))){
    listOfReactions <- xml_add_child(model_node, "listOfReactions")
  } else{
    listOfReactions <- xml_find_first(model_node, "./listOfReactions")
  }

  # Add new species to listOfSpecies
  reactantNames <- names(reactants)
  productNames <- names(products)
  if(is.null(reactantNames) && !(is.null(reactants))){
    message("Error: reactants are missing names")
    return()
  }
  if(is.null(productNames) && !(is.null(products))){
    message("Error: products are missing names")
    return()
  }
  reactionSpecies <- c(reactantNames, productNames)
  speciesNodes <- xml_find_all(listOfSpecies, "./species")
  speciesIDs <- xml_attr(speciesNodes, "id")
  for (species in reactionSpecies){
    if(!(species %in% speciesIDs)){
      xml_add_child(listOfSpecies, "species", id = species, compartment = "cell",
                    hasOnlySubstanceUnits = "false", boundaryCondition = "false", constant = "false")
    }
  }

  #Clean Reversible tag
  if(reversible){
    reversible <- "true"
  } else{
    reversible <- "false"
  }

  # Add Reaction to listOfReactions
  if(is.null(reactionName)){
    reactionName <- paste0("R", length(xml_find_all(listOfReactions, "./reaction"))+1)
  }
  reaction <- xml_add_child(listOfReactions, "reaction", id = reactionName, reversible = reversible)

  listOfReactants <- xml_add_child(reaction, "listOfReactants")
  for(reacNo in seq_len(length(reactants))){
    xml_add_child(listOfReactants, "speciesReference", species = reactantNames[reacNo],
                  stoichiometry = as.numeric(reactants[reacNo]), constant = "false")
  }
  listOfProducts <- xml_add_child(reaction, "listOfProducts")
  for(prodNo in seq_len(length(products))){
    xml_add_child(listOfProducts, "speciesReference", species = productNames[prodNo],
                  stoichiometry = as.numeric(products[prodNo]), constant = "false")
  }

  return(network)
}

