\name{KMeansSparseCluster.permute}
\alias{KMeansSparseCluster.permute}
\title{Choose tuning parameter for sparse k-means clustering}
\description{
The tuning parameter controls the L1 bound on w, the feature weights. A
permutation approach is used to select the tuning parameter.
}
\usage{
KMeansSparseCluster.permute(x, K=NULL, nperms = 25, wbounds = NULL,
silent = FALSE, nvals = 10, centers=NULL)
}
\arguments{
  \item{x}{The nxp data matrix, n is the number of observations and p
    the number of features.}
  \item{K}{The number of clusters desired - that is, the "K" in K-means
    clustering. Must specify K or centers.}
  \item{nperms}{Number of permutations.}
  \item{wbounds}{The range of tuning parameters to consider. This is the
  L1 bound on w, the feature weights. If NULL, then a range of values
  will be chosen automatically. Should be greater than 1 if non-null.}
  \item{silent}{Print out progress?}
  \item{nvals}{If wbounds is NULL, then the number of candidate tuning
    parameter values to consider.}
  \item{centers}{Optional argument. If you want to run the k-means
    algorithm starting from a particular set of clusters, then you
    can enter the Kxp matrix of cluster centers here. Default use
    case involves taking centers=NULL and instead specifying K.} 
}
\value{ 
  \item{gaps}{The gap statistics obtained (one for each of the tuning
    parameters tried). If O(s) is the objective function evaluated at
    the tuning parameter s, and O*(s) is the same quantity but for the
    permuted data, then Gap(s)=log(O(s))-mean(log(O*(s))).}
  \item{sdgaps}{The standard deviation of log(O*(s)), for each value of the
    tuning parameter s.}
  \item{nnonzerows}{The number of features with non-zero weights, for
    each value of the tuning parameter.}
  \item{wbounds}{The tuning parameters considered.}
  \item{bestw}{The value of the tuning parameter corresponding to the
    highest gap statistic.}
}
\details{
  Sparse k-means clustering seeks a p-vector of weights w (one per feature) and a set of
  clusters  C1,...,CK that optimize
  $maximize_C1,...,CK,w sum_j w_j BCSS_j$ subject to $||w||_2 <= 1,
  ||w||_1 <= s, w_j >= 0$,
  where $BCSS_j$ is the between cluster sum of squares for feature
  j, and s is a value for the L1 bound on w. Let O(s) denote the
  objective function with tuning parameter s: i.e. $O(s)=sum_j w_j
  BCSS_j$.

  We permute the data as follows: within each feature, we permute the
  observations. Using the permuted data, we can run sparse K-means with
  tuning parameter s, yielding the objective function O*(s). If we do
  this repeatedly we can get a number of O*(s) values.

  Then, the Gap statistic is given by $Gap(s)=log(O(s))-mean(log(O*(s)))$. The
  optimal s is that which results in the highest Gap statistic. Or, we
  can choose the smallest s such that its Gap statistic is within
  $sd(log(O*(s)))$ of the largest Gap statistic.
}
\references{Witten and Tibshirani (2009) A framework for feature
  selection in clustering.}
\author{Daniela M. Witten and Robert Tibshirani}
\seealso{\link{KMeansSparseCluster}, \link{HierarchicalSparseCluster}, \link{HierarchicalSparseCluster.permute}}
\examples{
# generate data
set.seed(11)
x <- matrix(rnorm(50*300),ncol=300)
x[1:25,1:50] <- x[1:25,1:50]+1
x <- scale(x, TRUE, TRUE)
# choose tuning parameter
km.perm <-
KMeansSparseCluster.permute(x,K=2,wbounds=seq(3,9,len=15),nperms=5)
print(km.perm)
plot(km.perm)
# run sparse k-means
km.out <- KMeansSparseCluster(x,K=2,wbounds=km.perm$bestw)
print(km.out)
plot(km.out)
# run sparse k-means for a range of tuning parameter values
km.out <- KMeansSparseCluster(x,K=2,wbounds=2:7)
print(km.out)
plot(km.out)
# Repeat, but this time start with a particular choice of cluster
# centers.
# This will do 4-means clustering starting with this particular choice
# of cluster centers.
km.perm.out <- KMeansSparseCluster.permute(x,wbounds=2:7, centers=x[1:4,])
print(km.out)
plot(km.out)
}





