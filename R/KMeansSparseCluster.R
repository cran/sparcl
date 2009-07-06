`KMeansSparseCluster` <-
function(x, K, wbounds=NULL, nstart=20, silent=FALSE, maxiter=6){
  # The criterion is : minimize_{w, C} sum_j w_j (WCSS_j - TSS_j) s.t. ||w||_2=1, ||w||_1<=s, w_j>=0
  # x is the data, nxp
  # K is the number of clusters desired
  # wbounds is a vector of L1 constraints on w, of the form  sum(abs(w))<=wbounds[i]
  if(is.null(wbounds)) wbounds <- seq(1, sqrt(ncol(x)), len=20)
  if(min(wbounds)<1) stop("wbounds should be greater than or equal to 1")
  wbounds <- c(wbounds) # In case wbounds is a single number, turn it into a vector
  out <- list()
  Cs <- kmeans(x, centers=K, nstart=nstart)$cluster
  for(i in 1:length(wbounds)){
    if(length(wbounds)>1 && !silent) cat(i,fill=FALSE)
    ws <- rep(1/sqrt(ncol(x)), ncol(x)) # Start with equal weights on each feature
    ws.old <- rnorm(ncol(x))
    store.bcss.ws <- NULL
    niter <- 0
    while((sum(abs(ws-ws.old))/sum(abs(ws.old)))>1e-4 && niter<maxiter){
      if(!silent) cat(niter, fill=FALSE)
      niter <- niter+1
      ws.old <- ws
      if(niter>1) Cs <- UpdateCs(x, K, ws, Cs) # if niter=1, no need to update!!
      ws <- UpdateWs(x, Cs, wbounds[i])
      store.bcss.ws <- c(store.bcss.ws, sum(GetWCSS(x, Cs)$bcss.perfeature*ws))
    }
    out[[i]] <- list(ws=ws, Cs=Cs, wcss=GetWCSS(x, Cs, ws), crit=store.bcss.ws, wbound=wbounds[i])
  }
  if(!silent) cat(fill=TRUE)
  if(length(wbounds)==1){
    out <- out[[1]]
    class(out) <- "kmeanssparse"
    return(out)
  }
  class(out) <- "multikmeanssparse"
  return(out)
}

plot.multikmeanssparse <- function(x,...){
  N <- length(x)
  par(mfrow=c(ceiling(N/2),2))
  for(i in 1:N){
    plot(x[[i]]$ws, main=paste("Wbound is ", sep="", round(x[[i]]$wbound,3)), xlab="Feature Index", ylab="Wj")
  }
}

plot.kmeanssparse <- function(x,...){
  plot(x$ws, main=paste("Wbound is ", sep="", round(x$wbound,3)), xlab="Feature Index", ylab="Wj")
}

PrintIt <- function(x){
  cat("Number of non-zero weights: ", sum(x$ws!=0), fill=TRUE)
  cat("Sum of weights: ", sum(x$ws), fill=TRUE)
  cat("Clustering: ", x$Cs, fill=TRUE)
  cat(fill=TRUE)
}

print.kmeanssparse <- function(x,...){
  cat("Wbound is ", x$wbound, ":", fill=TRUE)
  PrintIt(x)
}

print.multikmeanssparse <- function(x,...){
  for(i in 1:length(x)){
    cat("Wbound is ", x[[i]]$wbound, ":", fill=TRUE)
    PrintIt(x[[i]])
  }  
}
