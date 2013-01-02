\name{KMeansSparseCluster}
\alias{KMeansSparseCluster}
\alias{print.KMeansSparseCluster}
\alias{plot.KMeansSparseCluster}
\title{Performs sparse k-means clustering}
\description{
This function performs sparse k-means clustering. You must specify a number of clusters K and an L1 bound on w, the feature weights.
}
\details{
We seek a p-vector of weights w (one per feature) and a set of clusters
C1,...,CK that optimize

$maximize_C1,...,CK,w sum_j w_j BCSS_j$ subject to $||w||_2 <= 1,
||w||_1 <= wbound, w_j >= 0$

where $BCSS_j$ is the between cluster sum of squares for feature j. An
iterative approach is taken: with w fixed, optimize with respect to
C1,...,CK, and with C1,...,CK fixed, optimize with respect to w. Here,
wbound is a tuning parameter which determines the L1 bound on w.

The non-zero elements of w indicate features that are used in the sparse clustering.
}
\usage{
KMeansSparseCluster(x, K=NULL, wbounds = NULL, nstart = 20, silent =
FALSE, maxiter=6, centers=NULL)
\method{plot}{KMeansSparseCluster}(x,...)
\method{print}{KMeansSparseCluster}(x,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An nxp data matrix. There are n observations and p features.}
  \item{K}{The number of clusters desired ("K" in K-means
    clustering). Must provide either K or centers.}
  \item{wbounds}{A single L1 bound on w (the feature weights), or a vector of L1 bounds on w. If wbound is small,
    then few features will have non-zero weights. If wbound is large then all
  features will have non-zero weights. Should be greater than 1.}
  \item{nstart}{The number of random starts for the k-means
    algorithm.}
  \item{silent}{Print out progress?}
  \item{maxiter}{The maximum number of iterations.}
  \item{centers}{Optional argument. If you want to run the k-means
    algorithm starting from a particular set of clusters, then you can
    enter the Kxp matrix of cluster centers here. Default use case
    involves taking centers=NULL and instead specifying K.}
  \item{\dots}{ not used. }
}
\value{
  If wbounds is a vector, then a list with elements as follows (one per
  element of wbounds). If
  wbounds is just a single value, then elements as follows:
  \item{ws}{The p-vector of feature weights.}
  \item{Cs}{The clustering obtained.}
}
\references{Witten and Tibshirani (2009) A framework for feature
  selection in clustering.}
\author{Daniela M. Witten and Robert Tibshirani}
\seealso{\link{KMeansSparseCluster.permute},\link{HierarchicalSparseCluster}}
\examples{
# generate data
set.seed(11)
x <- matrix(rnorm(50*70),ncol=70)
x[1:25,1:20] <- x[1:25,1:20]+1
x <- scale(x, TRUE, TRUE)
# choose tuning parameter
km.perm <- KMeansSparseCluster.permute(x,K=2,wbounds=seq(3,7,len=15),nperms=5)
print(km.perm)
plot(km.perm)
# run sparse k-means
km.out <- KMeansSparseCluster(x,K=2,wbounds=km.perm$bestw)
print(km.out)
plot(km.out)
# run sparse k-means for a range of tuning parameter values
km.out <- KMeansSparseCluster(x,K=2,wbounds=seq(1.3,4,len=8))
print(km.out)
plot(km.out)
# Run sparse k-means starting from a particular set of cluster centers
#in the k-means algorithm.
km.out <- KMeansSparseCluster(x,wbounds=2:7,centers=x[c(1,3,5),])
}



