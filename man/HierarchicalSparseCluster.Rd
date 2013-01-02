\name{HierarchicalSparseCluster}
\alias{HierarchicalSparseCluster}
\alias{print.HierarchicalSparseCluster}
\alias{plot.HierarchicalSparseCluster}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Hierarchical sparse clustering}
\description{
Performs sparse hierarchical  clustering. If $d_ii'j$ is the
dissimilarity between observations i and i' for feature j, seek a sparse
weight
vector w and then use $(sum_j (d_ii'j w_j))_ii'$ as a nxn dissimilarity
matrix for hierarchical clustering.
}
\details{
  We seek a p-vector of weights w (one per feature) and a nxn matrix U
  that optimize

  $maximize_U,w sum_j w_j sum_ii' d_ii'j U_ii'$ subject to $||w||_2 <= 1,
  ||w||_1 <= wbound, w_j >= 0, sum_ii' U_ii'^2 <= 1$.

  Here, $d_ii'j$ is the dissimilarity between observations i and i' with
  along feature j. The resulting matrix U is used as a dissimilarity
  matrix for hierarchical clustering. "wbound" is a tuning parameter for
  this method, which controls the L1 bound on w, and as a result the
  number of features with non-zero $w_j$ weights.
  The non-zero elements of w indicate features that are used in the
  sparse clustering.

  We optimize the above criterion with an iterative approach: hold U
  fixed and optimize with respect to w. Then, hold w fixed and optimize
  with respect to U.

  Note that the arguments described as "Not for use" are included for
  the sparcl package to function with GenePattern but should be ignored
  by the R user.
}
\usage{
HierarchicalSparseCluster(x=NULL, dists=NULL,method=c("average","complete", "single","centroid"),
wbound=NULL,niter=15,dissimilarity=c("squared.distance","absolute.value"), uorth=NULL,silent=FALSE,
cluster.features=FALSE,method.features=c("average", "complete",
"single","centroid"),output.cluster.files=FALSE,
outputfile.prefix="output",genenames=NULL,genedesc=NULL,standardize.arrays=FALSE)
\method{print}{HierarchicalSparseCluster}(x,...)
\method{plot}{HierarchicalSparseCluster}(x,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A nxp data matrix; n is the number of observations and p the
    number of features. If NULL, then specify dists instead.}
  \item{dists}{For advanced users, can be entered instead of x. If
    HierarchicalSparseCluster has already been run on this data, then
    the dists value of the previous output can be entered here. Under normal circumstances, leave this argument NULL and
    pass in x instead.}
  \item{method}{The type of linkage to use in the hierarchical
    clustering - "single", "complete", "centroid", or "average".}
  \item{wbound}{The L1 bound on w to use; this is the tuning parameter
    for sparse hierarchical clustering. Should be greater than 1.}
  \item{niter}{The number of iterations to perform in the sparse
    hierarchical clustering algorithm.}
  \item{dissimilarity}{The type of dissimilarity measure to use. One of
    "squared.distance" or "absolute.value". Only use this if x was
    passed in (rather than dists).}
  \item{uorth}{If complementary sparse clustering is desired, then this
    is the nxn dissimilarity matrix obtained in the original sparse
    clustering.}
  \item{standardize.arrays}{Should the arrays be standardized? Default
    is FALSE.}
  \item{silent}{Print out progress?}
  \item{cluster.features}{Not for use.}
  \item{method.features}{Not for use.}
  \item{output.cluster.files}{Not for use.}
  \item{outputfile.prefix}{Not for use.}
  \item{genenames}{Not for use.}
  \item{genedesc}{Not for use.}
  \item{\dots}{ not used. }
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
\examples{
  # Generate 2-class data
  set.seed(1)
  x <- matrix(rnorm(100*50),ncol=50)
  y <- c(rep(1,50),rep(2,50))
  x[y==1,1:25] <- x[y==1,1:25]+2
  # Do tuning parameter selection for sparse hierarchical clustering
  perm.out <- HierarchicalSparseCluster.permute(x, wbounds=c(1.5,2:6),
nperms=5)
  print(perm.out)
  plot(perm.out)
  # Perform sparse hierarchical clustering
  sparsehc <- HierarchicalSparseCluster(dists=perm.out$dists,
wbound=perm.out$bestw, method="complete")
  # faster than   sparsehc <- HierarchicalSparseCluster(x=x,wbound=perm.out$bestw, method="complete")
  par(mfrow=c(1,2))
  plot(sparsehc)
  plot(sparsehc$hc, labels=rep("", length(y)))
  print(sparsehc)
  # Plot using knowledge of class labels in order to compare true class
  #   labels to clustering obtained
  par(mfrow=c(1,1))
  ColorDendrogram(sparsehc$hc,y=y,main="My Simulated Data",branchlength=.007)
  # Now, what if we want to see if out data contains a *secondary*
  #   clustering after accounting for the first one obtained. We
  #   look for a complementary sparse clustering:
  sparsehc.comp <- HierarchicalSparseCluster(x,wbound=perm.out$bestw,
     method="complete",uorth=sparsehc$u)
  # Redo the analysis, but this time use "absolute value" dissimilarity:
  perm.out <- HierarchicalSparseCluster.permute(x, wbounds=c(1.5,2:6),
    nperms=5, dissimilarity="absolute.value")
  print(perm.out)
  plot(perm.out)
  # Perform sparse hierarchical clustering
  sparsehc <- HierarchicalSparseCluster(dists=perm.out$dists, wbound=perm.out$bestw, method="complete", dissimilarity="absolute.value")
  par(mfrow=c(1,2))
  plot(sparsehc)
}
