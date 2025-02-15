
#' @export
#' @rdname cracipeNetwork
#' @aliases cracipeNetwork
setMethod(f="cracipeNetwork",
          signature="cRacipeSE",
          definition=function(.object)
          {
            stoichMatrix <- as.matrix(rowData(.object))
            reactantMatrix <- metadata(.object)$reactantMatrix
            return(list(stoichMatrix, reactantMatrix))
          }
)


#' @rdname cracipeNetwork-set
#' @aliases cracipeNetwork-set
setMethod("cracipeNetwork<-", signature(.object= "cRacipeSE"),
          function(.object, value) {

            if(is(value, "xml_document")){
              #Strip namespace to simplify process
              xml_ns_strip(value)

              modelNode <- xml_find_first(value, "./model")
              netID <- xml_attr(modelNode, "id")
              if (is.na(netID)) {
                stop("Error: <model> node is missing the 'id' attribute.")
              }

              listOfSpecies <- xml_find_first(modelNode, "./listOfSpecies")
              listOfReactions <- xml_find_first(modelNode, "./listOfReactions")

              speciesNodes <- xml_find_all(listOfSpecies, "./species")
              speciesIDs <- xml_attr(speciesNodes, "id")
              nSpecies <- length(speciesIDs)

              reactionNodes <- xml_find_all(listOfReactions, "./reaction")
              reactionIDs <- xml_attr(reactionNodes, "id")
              numReactions <- length(reactionIDs)

              #Build Stoichiometry matrix and Reactant Matrix
              stoichMatrix <- matrix(data = 0L, nrow = nSpecies, ncol = numReactions)
              rownames(stoichMatrix) <- speciesIDs
              colnames(stoichMatrix) <- reactionIDs
              reactantMatrix <- stoichMatrix

              for(reactionNode in reactionNodes){
                reactionID <- xml_attr(reactionNode, "id")
                isReversible <- as.logical(xml_attr(reactionNode, "reversible"))
                listOfProducts <- xml_find_first(reactionNode, "./listOfProducts")
                productNodes <- xml_find_all(listOfProducts, "./speciesReference")

                if(isReversible){
                  revName <- paste0(reactionID, "Rev")
                  revReactant <- reactantMatrix[, reactionID, drop=FALSE]
                  colnames(revReactant) <- revName
                  numReactions <- numReactions + 1
                }

                for(speciesRef in productNodes){
                  speciesName <- xml_attr(speciesRef, "species")
                  speciesStoich <- xml_attr(speciesRef, "stoichiometry")
                  stoichMatrix[speciesName, reactionID] <- (stoichMatrix[speciesName, reactionID]
                                                            + as.numeric(speciesStoich))
                  if(isReversible){
                    revReactant[speciesName, revName] <- (revReactant[speciesName, revName] +
                                                            as.numeric(speciesStoich))
                  }
                }

                listOfReactants <- xml_find_first(reactionNode, "./listOfReactants")
                reactantNodes <- xml_find_all(listOfReactants, "./speciesReference")

                for(speciesRef in reactantNodes){
                  speciesName <- xml_attr(speciesRef, "species")
                  speciesStoich <- xml_attr(speciesRef, "stoichiometry")
                  stoichMatrix[speciesName, reactionID] <- (stoichMatrix[speciesName, reactionID]
                                                            - as.numeric(speciesStoich))
                  reactantMatrix[speciesName, reactionID] <- (reactantMatrix[speciesName, reactionID]
                                                            + as.numeric(speciesStoich))
                }

                if(isReversible){
                  revStoich <- -stoichMatrix[, reactionID, drop=FALSE]
                  colnames(revStoich) <- revName
                  stoichMatrix <- cbind(stoichMatrix, revStoich)
                  reactantMatrix <- cbind(reactantMatrix, revReactant)
                }
              }
            }else if(is(value, "character")){
              if(file.exists(value)){
                netID <- basename(value)
                netID <- tools::file_path_sans_ext(netID)

                con <- file(value, "r")
                numReactions <- length(utils::count.fields(value, sep = "\n"))

                reactionList <- list()
                reactNum <- 0

                #Read reactions
                while(TRUE) {
                  reaction <- readLines(con, n = 1, warn = FALSE)
                  if (length(reaction) == 0) break

                  reactNum <- reactNum + 1

                  reversible <- grepl("<->", reaction, fixed = TRUE)
                  if(reversible) {numReactions <- numReactions + 1}
                  delimiter <- if (reversible) "<->" else "->"
                  parts <- strsplit(reaction, delimiter)[[1]]
                  #handle case of no products
                  if(length(parts) == 1){
                    parts[[2]] <- ""
                  }
                  if(length(parts) != 2){
                    close(con = con)
                    stop(".tpo file incorrectly configured")
                  }
                  reactionInfo <- list(reactants = NULL, products = NULL,
                                          reversible = reversible)

                  for(type in 1:2){
                    side <- parts[[type]]
                    if (nchar(side) == 0){
                      reactionInfo[[type]] <- stats::setNames(integer(0),
                                                          character(0))
                      next
                    }
                    elements <- strsplit(trimws(side), " ")[[1]]
                    if (length(elements) %% 2 != 0){
                      close(con = con)
                      stop(".tpo file incorrectly configured")
                    }
                    numbers <- as.integer(elements[seq(1, length(elements), by = 2)])
                    species <- elements[seq(2, length(elements), by = 2)]
                    numbers <- stats::setNames(numbers, species)
                    reactionInfo[[type]] <- numbers
                  }

                  reactionList[[reactNum]] <- reactionInfo
                }

                close(con = con)

                #Generate reacttant matrix and stoichiometry matrix
                allSpecies <- unique(unlist(lapply(reactionList, function(r)
                  c(names(r$reactants), names(r$products)))))
                nSpecies <- length(allSpecies)

                # Generate reaction names
                reactionLabels <- character(numReactions)
                reacIdx <- 1
                for (i in seq_along(reactionList)) {
                  reactionLabel <- paste0("R", i)
                  reactionLabels[reacIdx] <- reactionLabel

                  if (reactionList[[i]]$reversible) {
                    reverseLabel <- paste0("R", i, "Rev")
                    reactionLabels[reacIdx + 1] <- reverseLabel

                    reacIdx <- reacIdx + 1
                  }
                  reacIdx <- reacIdx + 1
                }

                stoichMatrix <- matrix(0, nrow = length(allSpecies), ncol = numReactions,
                                          dimnames = list(allSpecies, reactionLabels))
                reactantMatrix <- stoichMatrix

                colIdx <- 1
                for (r in reactionList) {
                  reactionSpecies <- unique(c(names(r$reactants), names(r$products)))
                  for (species in reactionSpecies) {
                    prodCoef <- ifelse(species %in% names(r$products), r$products[[species]], 0)
                    reacCoef <- ifelse(species %in% names(r$reactants), r$reactants[[species]], 0)
                    stoichMatrix[species, colIdx] <- prodCoef - reacCoef
                    reactantMatrix[species, colIdx] <- reacCoef
                  }
                  if (r$reversible) {
                    colIdx <- colIdx + 1
                    for (species in reactionSpecies) {
                      prodCoef <- ifelse(species %in% names(r$reactants), r$reactants[[species]], 0)
                      reacCoef <- ifelse(species %in% names(r$products), r$products[[species]], 0)
                      stoichMatrix[species, colIdx] <- prodCoef - reacCoef
                      reactantMatrix[species, colIdx] <- reacCoef
                    }
                  }
                  colIdx <- colIdx + 1
                }

              } else{
                stop("File not found")
              }
            } else{
              stop("Incorrect network file format. It should be either an SBML
                   xml_document or a file path to a .tpo file")
            }

            configData <- NULL
            data("configData",envir = environment(), package = "racipeCRN")

            .object <- cRacipeSE(
              assays = SimpleList(matrix(NA, nrow = nSpecies,ncol = configData$simParams["numModels"])),
              rowData = DataFrame(stoichMatrix),
              colData = DataFrame(matrix(NA,nrow = configData$simParams["numModels"],ncol=0)),
              metadata = list(
                annotation = netID,
                numReactions = numReactions,
                config = configData,
                reactantMatrix = reactantMatrix)
            )
            message("network successfully loaded")
            return(.object)
          }
)

#' @rdname cracipeParams
#' @aliases cracipeParams
setMethod("cracipeParams", signature("cRacipeSE"), function(.object) {
  return(as(colData(.object)[,seq_len(metadata(.object)$numReactions),drop=FALSE], "data.frame"))
})

#' @rdname cracipeIC
#' @aliases cracipeIC
setMethod(f="cracipeIC",
          signature="cRacipeSE",
          definition=function(.object)
          {
            numReactions <- metadata(.object)$numReactions
            return(t(as.data.frame(colData(.object)[,(numReactions+1):
                                                      (numReactions + length(names(.object)))])))
          }
)
