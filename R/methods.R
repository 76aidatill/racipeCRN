#' @export
#' @import SummarizedExperiment
#' @importMethodsFrom BiocGenerics annotation
#' @title  A method to get the annotation
#' @param object cRacipeSE object
#' @param ... Additional arguments, for use in specific methods.
#' @examples
#'
#' cSet <- cRacipeSE()
#' wilhelm <- xml2::read_xml(system.file("extdata", "wilhelm.xml", package = "racipeCRN"))
#' cracipeNetwork(cSet) <- wilhelm
#' ann <- annotation(cSet)
#' rm(cSet, wilhelm,ann)
#' @return character
#'
setMethod("annotation", signature(object = "cRacipeSE"), function(object,...) {
  return(metadata(object)$annotation)
})

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

#' @rdname cracipeConfig
#' @aliases cracipeConfig
setMethod(f="cracipeConfig",
          signature="cRacipeSE",
          definition=function(.object)
          {
            return(metadata(.object)$config)
          }
)

#' @rdname cracipeConverge
#' @aliases cracipeConverge
setMethod(f="cracipeConverge",
          signature="cRacipeSE",
          definition=function(.object)
          {
            configuration <- cracipeConfig(.object)
            stoichMatrix <- rowData(.object)
            numReactions <- ncol(stoichMatrix)
            return(as.data.frame(colData(.object)[,(numReactions +
                                                      length(names(.object)) + 1)
                                                  :(dim(colData(.object))[2])]))
          }
)

#' @export
#' @rdname cracipeUniqueStates
#' @aliases cracipeUniqueStates
setMethod(f="cracipeUniqueStates",
          signature = "cRacipeSE",
          definition = function(.object, roundDigits = 3){
            speciesConc <- assay(.object)
            objMetadata <- metadata(.object)
            configuration <- cracipeConfig(.object)

            converge <- cracipeConverge(.object)
            numIC <- configuration$simParams["numIC"]
            numModels <- configuration$simParams["numModels"]

            uniqueConcList <- list() #stores unique expressions in a list

            for(modelCount in seq_len(numModels)){
              modelStates <- data.frame()
              #grab ICs and convergence data for each model
              startIdx <- (modelCount - 1)*numIC + 1
              endIdx <- modelCount*numIC

              finalModelExpressions <- speciesConc[, startIdx:endIdx]
              #filters out non-converging states
              ICconvergences <- converge[startIdx:endIdx, 1]
              convergedICs <- finalModelExpressions[, as.logical(ICconvergences)]

              if(!is.null(ncol(convergedICs))){
                #Taking unique states up until roundDigits
                uniqueIdx <- which(!(duplicated(round(convergedICs, digits = roundDigits), MARGIN = 2)))
                uniqueExprx <- convergedICs[,uniqueIdx]
                uniqueConcList[[modelCount]] <- as.data.frame(uniqueExprx)
              }else{
                uniqueConcList[[modelCount]] <- data.frame()
              }
            }
            return(uniqueConcList)
          }
)

#' @rdname cracipeNormalize
#' @aliases cracipeNormalize
setMethod(f="cracipeNormalize",
          signature="cRacipeSE",
          definition=function(.object)
          {
            metadataTmp <- metadata(.object)
            assayDataTmp <- assay(.object)
            speciesConcentration <- assayDataTmp
            speciesConcentration <- log2(1+speciesConcentration)
            means <- rowMeans(speciesConcentration, na.rm = TRUE)
            sds <-  apply(speciesConcentration, 1, function(x) sd(x, na.rm=TRUE))
            speciesConcentration <- sweep(speciesConcentration, 1, means, FUN = "-")
            speciesConcentration <- sweep(speciesConcentration, 1, sds, FUN = "/")
            assayDataTmp2 <- speciesConcentration


            metadataTmp$normalized <- TRUE
            assay(.object) <- assayDataTmp2
            metadata(.object) <- metadataTmp
            return(.object)
          }
)


#' @rdname cracipePlotData
#' @aliases cracipePlotData
setMethod(f="cracipePlotData",
          signature="cRacipeSE",
          definition=function(.object, plotToFile = TRUE, nClusters = 2,
                              heatmapPlot = TRUE,
                              pcaPlot = TRUE, umapPlot = TRUE,
                              clustMethod = "ward.D2", col = col,
                              distType = "euclidean",
                              assignedClusters = NULL,
                              corMethod = "spearman", ...)
          {


            if(missing(col)){
              col <-  grDevices::colorRampPalette(rev(
                RColorBrewer::brewer.pal(11, 'Spectral')))
            }
            col2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
                      "#D55E00", "#CC79A7")
            p <- list()
            ts <- list()
            tsCounter <- 1

            i=1;

            metadataTmp <- metadata(.object)

            #Get unique steady states
            uniqueStates <- cracipeUniqueStates(.object)
            assayDataTmp <- Filter(function(df) !(ncol(df) == 0), uniqueStates)
            assayDataTmp <- do.call(cbind, assayDataTmp)

            if(missing(distType)){
              corCol <- stats::cor((assayDataTmp), method = corMethod)
              distanceCol <- stats::as.dist((1 - corCol) / 2)
              corRow <- stats::cor(t(assayDataTmp), method = corMethod)
              distanceRow <- stats::as.dist((1 - corRow) / 2)
            }
            else{
              # distType = "manhattan"
              distanceCol <- stats::dist(t(assayDataTmp), method = distType)
              distanceRow <- stats::dist((assayDataTmp), method = distType)
            }
            clustersCol <- stats::hclust(distanceCol, method = clustMethod)
            ddCol <- as.dendrogram(clustersCol)


            clustersRow <- stats::hclust(distanceRow, method = clustMethod)
            ddRow <- stats::as.dendrogram(clustersRow)
            if(is.null(assignedClusters)){
              clustCut <- stats::cutree(clustersCol, nClusters)
              clustColors <- col2[clustCut]
              assignedClusters <- clustCut
            }

            if(!missing(assignedClusters)){
              clustNames <- unique(assignedClusters)
              nClusters <- length(clustNames)
              clustColors <- numeric(length(assignedClusters))
              for(tmp1 in seq_len(length(clustColors))){
                clustColors[tmp1] <- which(clustNames == assignedClusters[tmp1] )
              }
              clustColors <- col2[clustColors]
            }
            names(clustColors) <- assignedClusters

            if(plotToFile){
              fileName <- paste0(annotation(.object),"_heatmap.pdf")

              pdf(fileName, onefile = TRUE)
            }
            if(heatmapPlot) {
              gplots::heatmap.2((as.matrix(assayDataTmp)),
                                Colv = ddCol,
                                Rowv = ddRow,
                                trace = "none",
                                col = col,
                                ColSideColors = clustColors
              )
            }


            if(plotToFile){

              dev.off()
            }
            V1 <- NULL
            V2 <- NULL
            PC1 <- NULL
            PC2 <- NULL
            if(umapPlot){
              umapGE <- umap::umap(t(assayDataTmp))
              p[[i]] <-
                ggplot2::ggplot(data = as.data.frame(umapGE$layout)) +
                geom_point(aes(x = V1, y=V2), color = clustColors, shape = 1) +
                labs(x = "Umap1", y="Umap2") +
                theme(text = element_text(size=30),
                      panel.background = element_rect(fill = "white", color = "black"),
                      panel.grid.major = element_line(color="gray", size=0.25))
              #  panel.border = element_rect(color = "black"))
              i <- i+1
            }

            if(pcaPlot){

              pca1 = summary(prcomp(t(assayDataTmp), scale. = FALSE))
              p[[i]] <-
                ggplot2::ggplot(data = as.data.frame(pca1$x)) +
                geom_point(aes(x = PC1, y=PC2), color = clustColors, shape = 1) +
                labs(x = paste0("PC1(",100*pca1$importance[2,1],"%)"),
                     y=paste0("PC2(",100*pca1$importance[2,2],"%)")) +
                theme(text = element_text(size=30),
                      panel.background = element_rect(fill = "white", color = "black"),
                      panel.grid.major = element_line(color="gray", size=0.25))

              if(!is.null(metadataTmp$tsSimulations)){
                tsPca <- assayDataTmp[
                  2:(1+length(metadataTmp$tsSimulations))]
                for(j in seq_len(length(tsPca))){
                  tsPca[[j]] <-
                    t(scale(t(tsPca[[j]]), pca1$center, pca1$scale) %*%
                        pca1$rotation)

                  ts[[j]] <-
                    ggplot2::ggplot(data = as.data.frame(t(tsPca[[j]]))) +
                    geom_point(aes(x = PC1, y=PC2), shape = 1) +
                    labs(x = paste0("PC1(",100*pca1$importance[2,1],"%)"),
                         y=paste0("PC2(",100*pca1$importance[2,2],"%)"),
                         title = names(metadataTmp$tsSimulations)[j]) +
                    theme(text = element_text(size=30),
                          panel.background = element_rect(fill = "white",
                                                          color = "black"),
                          panel.grid.major = element_line(color="gray", size=0.25))
                }

              }
            }

            if(plotToFile){
              fileName <- paste0(annotation(.object),"_plots.pdf")
              pdf(fileName, onefile = TRUE)
            }

            for (i in seq(length(p))) {
              gridExtra::grid.arrange(p[[i]])
              # do.call("grid.arrange", p[[i]])
            }

            if(plotToFile){
              message("Plots saved as pdf files in the working directory.")
              dev.off()
            }


            if(umapPlot)metadataTmp$umap <- umapGE
            if(pcaPlot) metadataTmp$pca <- pca1
            metadataTmp$assignedClusters <- assignedClusters
            metadata(.object) <- metadataTmp
            return(.object)
          }
        )

#' @export
#' @rdname cracipeConvergeDist
#' @aliases cracipeConvergeDist
setMethod(f="cracipeConvergeDist",
          signature="cRacipeSE",
          definition=function(.object, plotToFile = FALSE)
          {
            metadataTmp <- metadata(.object)
            configuration <- cracipeConfig(.object)

            converging <- cracipeConverge(.object)
            numModels <- configuration$simParams["numModels"]
            numIC <- configuration$simParams["numIC"]
            numIterations <- configuration$simParams["numIterations"]
            numStepsIteration <- configuration$simParams["numStepsIteration"]
            integrateStepSize <- configuration$simParams["integrateStepSize"]
            numExprx <- numModels*numIC
            iterationLength <- numStepsIteration*integrateStepSize

            #Initialize proportions
            convergedProportions <- numeric(numIterations)

            if(plotToFile){
              fileName <- paste0(annotation(.object),"_ConvergDist.pdf")
              pdf(fileName) #Opens graphics object to store file in
            }

            #Getting rid of non-converged models
            convergedICs <- converging[converging[, 1] == 1, ]
            testScores <- convergedICs[,2]


            for (i in 1:numIterations){
              convergedProportions[i] <- sum(testScores <= i) / numExprx
            }

            #convert iterations to time
            xvals <- seq(1,numIterations)
            xvals <- xvals*iterationLength

            title = paste0("Time plot of stable ", annotation(.object), " concentrations")
            plot(xvals, convergedProportions, type="l", col="blue",
                 xlab="Time Units", ylab = "Fraction Converged",
                 main = title, ylim = c(0,1))
            graphics::polygon(c(xvals, rev(xvals)),
                              c(convergedProportions, rep(0, length(convergedProportions))),
                              col = rgb(0, 0, 1, 0.5), border = NA)

            if(plotToFile){
              message("Plot saved as pdf files in the working directory.")
              dev.off() #closes graphics object and send it to working directory
            }

            metadataTmp$stableProportion <- convergedProportions[numIterations]
            if(convergedProportions[numIterations] > 0.99){
              ninetyNineIndex <- which(convergedProportions > 0.99)[1]
              metadataTmp$ninetyNineConvergedNum <- ninetyNineIndex
            }
            metadata(.object) <- metadataTmp

            return(.object)
          }

)

