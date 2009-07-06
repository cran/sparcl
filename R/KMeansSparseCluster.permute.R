`KMeansSparseCluster.permute` <-
function(x, K=2,  nperms=25, wbounds=NULL,silent=FALSE, nvals=10){
  if(is.null(wbounds)) wbounds <- exp(seq(log(1.2), log(sqrt(ncol(x))*.9), len=nvals))
  # was seq(1.2, sqrt(ncol(x))*.6, len=10)
  permx <- list()
  nnonzerows <- NULL
  for(i in 1:nperms){
    permx[[i]] <- matrix(NA, nrow=nrow(x), ncol=ncol(x))
    for(j in 1:ncol(x)) permx[[i]][,j] <- sample(x[,j])
  }
  tots <- NULL
  out <- KMeansSparseCluster(x, K, wbounds=wbounds, silent=silent)
  for(i in 1:length(out)){
    nnonzerows <- c(nnonzerows, sum(out[[i]]$ws!=0))
    bcss <- GetWCSS(x,out[[i]]$Cs)$bcss.perfeature 
    tots <- c(tots, sum(out[[i]]$ws*bcss))
  }
  permtots <- matrix(NA, nrow=length(wbounds), ncol=nperms)
  for(k in 1:nperms){
    if(!silent) cat("Permutation ", k, "of ", nperms, fill=TRUE)
    perm.out <- KMeansSparseCluster(permx[[k]], K, wbounds=wbounds, silent=silent)
    for(i in 1:length(perm.out)){
      perm.bcss <- GetWCSS(permx[[k]], perm.out[[i]]$Cs)$bcss.perfeature
      permtots[i,k] <- sum(perm.out[[i]]$ws*perm.bcss)
    }
  }
  gaps <- (log(tots)-apply(log(permtots),1,mean))
  out <- list(tots=tots, permtots=permtots, nnonzerows=nnonzerows, gaps=gaps, sdgaps=apply(log(permtots),1,sd), wbounds=wbounds, bestw=wbounds[which.max(gaps)])
  if(!silent) cat(fill=TRUE)
  class(out) <- "kmeanssparseperm"
  return(out)
}

print.kmeanssparseperm <- function(x,...){
  cat("Tuning parameter selection results for Sparse K-means Clustering:", fill=TRUE)
  mat <- round(cbind(x$wbounds, x$nnonzerows, x$gaps, x$sdgaps),4)
  dimnames(mat) <- list(1:length(x$wbounds), c("Wbound", "# Non-Zero W's", "Gap Statistic", "Standard Deviation"))
  print(mat, quote=FALSE)                        
  cat("Tuning parameter that leads to largest Gap statistic: ", x$bestw, fill=TRUE)
}

plot.kmeanssparseperm <- function(x,...){
  plot(x$nnonzerows, x$gaps, log="x", main="Gap Statistics", xlab="# Non-zero Wj's", ylab="")
  lines(x$nnonzerows, x$gaps)
}

