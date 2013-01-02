\name{HierarchicalSparseCluster.permute}
\alias{HierarchicalSparseCluster.permute}
\alias{print.HierarchicalSparseCluster.permute}
\alias{plot.HierarchicalSparseCluster.permute}
\title{Choose tuning parameter for sparse hierarchical clustering}
\description{
 The tuning parameter controls the L1 bound on w, the feature weights. A permutation approach is used to select the tuning parameter.
}
\usage{
HierarchicalSparseCluster.permute(x, nperms = 10, wbounds = NULL,
dissimilarity=c("squared.distance",
"absolute.value"),standardize.arrays=FALSE)
\method{plot}{HierarchicalSparseCluster.permute}(x,...) 
\method{print}{HierarchicalSparseCluster.permute}(x,...)
}
\arguments{
  \item{x}{A nxp data matrix, with n observations and p feaures.}
  \item{nperms}{The number of permutations to perform.}
  \item{wbounds}{The sequence of tuning parameters to consider. The
    tuning parameters are the L1 bound on w, the feature weights. If
    NULL, then a default sequence will be used. If non-null, should be
    greater than 1.}
  \item{dissimilarity}{How should dissimilarity be computed? Default is
    squared.distance.}
  \item{standardize.arrays}{Should the arrays first be standardized?
    Default is FALSE.}
  \item{\dots}{ not used. }
}
\details{
  Let $d_ii'j$ denote the dissimilarity between observations i and i'
  along feature j.
  
  Sparse hierarchical clustering seeks a p-vector of weights w (one per
  feature) and a nxn  matrix U that optimize
  $maximize_U,w sum_j w_j  sum_ii' d_ii'j U_ii'$ subject to $||w||_2 <= 1,
  ||w||_1 <= s, w_j >= 0, sum_ii' U_ii'^2 <= 1$,
  where  s is a value for the L1 bound on w. Let O(s) denote the
  objective function with tuning parameter s: i.e. $O(s)=sum_j w_j
  sum_ii' d_ii'j U_ii'$.

  We permute the data as follows: within each feature, we permute the
  observations. Using the permuted data, we can run sparse hierarchical clustering
  with tuning parameter s, yielding the objective function O*(s). If we do
  this repeatedly we can get a number of O*(s) values.

  Then, the Gap statistic is given by
  $Gap(s)=log(O(s))-mean(log(O*(s)))$. The
  optimal s is that which results in the highest Gap
  statistic. Or, we
  can choose the smallest s such that its Gap statistic is within
  $sd(log(O*(s)))$ of the largest Gap statistic.
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
\references{Witten and Tibshirani (2009) A framework for feature
  selection in clustering.}
\author{Daniela M. Witten and Robert Tibshirani}
\seealso{\link{HierarchicalSparseCluster}, \link{KMeansSparseCluster}, \link{KMeansSparseCluster.permute}}
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
  sparsehc <- HierarchicalSparseCluster(dists=perm.out$dists, wbound=perm.out$bestw, method="complete")
  par(mfrow=c(1,2))
  plot(sparsehc)
  plot(sparsehc$hc, labels=rep("", length(y)))
  print(sparsehc)
  # Plot using knowledge of class labels in order to compare true class
  #   labels to clustering obtained
  par(mfrow=c(1,1))
  ColorDendrogram(sparsehc$hc,y=y,main="My Simulated
Data",branchlength=.007)
}

