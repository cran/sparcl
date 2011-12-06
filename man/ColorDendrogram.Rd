\name{ColorDendrogram}
\alias{ColorDendrogram}
\title{Color the leaves in a hierarchical clustering dendrogram}
\description{
Pass in the output of "hclust" and a class label for each
observation. A colored dendrogram will result, with the leaf colors indicating
the classes.
}
\usage{
ColorDendrogram(hc, y, main = "", branchlength = 0.7, labels = NULL, xlab = NULL, sub = NULL, ylab = "", cex.main = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hc}{The output of running "hclust" on a nxn dissimilarity matrix}
  \item{y}{A vector of n class labels for the observations that were clustered
    using "hclust". If labels are numeric from 1 to K, then colors will
    be determine automatically. Otherwise the labels can take the form
    of colors (e.g. c("red", "red", "orange", "orange")).}
  \item{main}{The main title for the dendrogram.}
  \item{branchlength}{How long to make the colored part of the
    branches. Adjustment will be needed for each dissimilarity matrix}
  \item{labels}{The labels for the n observations.}
  \item{xlab}{X-axis label.}
  \item{sub}{Sub-x-axis label.}
  \item{ylab}{Y-axis label.}
  \item{cex.main}{The amount by which to enlarge the main title for the figure.}
}
\references{Witten and Tibshirani (2009) A framework for feature
  selection in clustering.}
\author{Daniela M. Witten and Robert Tibshirani}
\seealso{HierarchicalSparseCluster, HierarchicalSparseCluster.permute}
\examples{
# Generate 2-class data
set.seed(1)
x <- matrix(rnorm(100*20),ncol=20)
y <- c(rep(1,50),rep(2,50))
x[y==1,] <- x[y==1,]+2
# Perform hierarchical clustering
hc <- hclust(dist(x),method="complete")
# Plot
ColorDendrogram(hc,y=y,main="My Simulated Data",branchlength=3)
}
