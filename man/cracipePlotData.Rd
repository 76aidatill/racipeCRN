% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/methods.R
\name{cracipePlotData}
\alias{cracipePlotData}
\alias{cracipePlotData,cRacipeSE-method}
\title{Plot sRACIPE data}
\usage{
cracipePlotData(
  .object,
  plotToFile = FALSE,
  nClusters = 2,
  heatmapPlot = TRUE,
  pcaPlot = TRUE,
  umapPlot = TRUE,
  clustMethod = "ward.D2",
  col = col,
  distType = "euclidean",
  assignedClusters = NULL,
  corMethod = "spearman",
  ...
)

\S4method{cracipePlotData}{cRacipeSE}(
  .object,
  plotToFile = TRUE,
  nClusters = 2,
  heatmapPlot = TRUE,
  pcaPlot = TRUE,
  umapPlot = TRUE,
  clustMethod = "ward.D2",
  col = col,
  distType = "euclidean",
  assignedClusters = NULL,
  corMethod = "spearman",
  ...
)
}
\arguments{
\item{.object}{List A list returned by \code{\link{cracipeSimulate}} function}

\item{plotToFile}{(optional) logical. Default \code{FALSE}. Whether to save
plots to a file.}

\item{nClusters}{(optional) Integer. Default 2. Expected number of clusters
in the simulated data. Hierarchical clustering will be used to cluster the
data and the the models will be colored in UMAP and PCA plots according to
these clustering results. The clusters can be also supplied using
\code{assignedClusters}.}

\item{heatmapPlot}{(optional) logical. Default \code{TRUE}. Whether to plot
hierarchichal clustering.}

\item{pcaPlot}{(optional) logical. Default \code{TRUE}. Whether to plot PCA
embedding.}

\item{umapPlot}{(optional) logical. Default \code{TRUE}. Whether to plot
UMAP embedding}

\item{clustMethod}{(optional) character. Default \code{"ward.D2"}. Clustering
method for heatmap. See \code{\link[gplots]{heatmap.2}}}

\item{col}{(optional) Color palette}

\item{distType}{(optional) Distance type.  Used only if specified
explicitly. Otherwise, 1-cor is used. See \code{\link[stats]{dist}},
\code{\link[stats]{hclust}}}

\item{assignedClusters}{vector integer or character. Default NULL.
Cluster assignment of models.}

\item{corMethod}{(optional) character. Default \code{"spearman"}. Correlation
method for distance function.}

\item{...}{Other arguments}
}
\value{
\code{cRacipeSE} object
}
\description{
Plots heatmap, pca, umap of the data simulated using racipeCRN.
Only considers unique steady states for each model.
}
\section{Related Functions}{


\code{\link{cracipeSimulate}},  \code{\link{cracipeUniqueStates}}
}

\examples{
wilhelm <- xml2::read_xml(system.file("extdata", "wilhelm.xml", package = "racipeCRN"))
\dontrun{
cSet <- racipeCRN::cracipeSimulate(network = wilhelm, numModels = 20,
integrateStepSize = 0.001, numIterations = 30)
cSet <- cracipePlotData(cSet)
}
}
