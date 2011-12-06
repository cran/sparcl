\name{HierarchicalSparseCluster.wrapper}
\alias{HierarchicalSparseCluster.wrapper}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A wrapper for the hierarchical sparse clustering algorithm}
\description{
A wrapper for HierarchicalSparseCluster which reads in the data in  GCT file
format, and then  automatically chooses the
optimal tuning parameter value using HierarchicalSparseCluster.permute
if not specified.
}
\usage{
HierarchicalSparseCluster.wrapper(file,  method=c("average", "complete", "single", "centroid"),
wbound=NULL, silent=FALSE, cluster.features=FALSE,
method.features=c("average", "complete",
"single","centroid"),output.cluster.files=TRUE,outputfile.prefix=NULL,maxnumgenes=5000,standardize.arrays=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{A GCT filename in the working directory containing the
    data to be clustered.}
  \item{method}{The type of linkage to use in the hierarchical
    clustering - "single", "complete", "average", or "centroid".}
  \item{wbound}{The L1 bound on w to use; this is the tuning parameter
    for sparse hierarchical clustering. If NULL, then it will be
    chosen via HierarchicalSparseCluster.permute.}
  \item{silent}{Print out progress?}
  \item{cluster.features}{Is a clustering for the features with non-zero
    weights also desired? Default is FALSE.}
  \item{method.features}{If cluster.features is TRUE, then the type of
    linkage used to cluster the features with non-zero weights: one of
    "single", "complete", "average", or "centroid".}
  \item{output.cluster.files}{Should files containing the clustering be output? Default is
    TRUE.}
  \item{outputfile.prefix}{The prefix for the output files. If NULL,
    then the prefix of the input file is used.}
  \item{maxnumgenes}{Limit the analysis to some number of genes with highest marginal
    variance, for computational reasons. This is recommended when the number of genes is very
    large. If NULL, then all genes are used.}
  \item{standardize.arrays}{Should the arrays first be standardized?
    Default is TRUE.}
}
\value{
  \item{hc}{The output of a call to "hclust", giving the results of
    hierarchical sparse clustering.}
  \item{ws}{The p-vector of feature weights.}
  \item{u}{The nxn dissimilarity matrix passed into hclust, of the form
    $(sum_j w_j d_ii'j)_ii'$.}
  \item{dists}{The (n*n)xp dissimilarity matrix for the data matrix
    x. This is useful if additional calls to HierarchicalSparseCluster
    will be made.}
}
\references{Witten and Tibshirani (2009) A framework for feature
  selection in clustering.}
\author{Daniela M. Witten and Robert Tibshirani}
\seealso{\link{HierarchicalSparseCluster.permute},\link{KMeansSparseCluster},\link{KMeansSparseCluster.permute}}

