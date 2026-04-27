library(foreach)
library(doSNOW)
library(progress)

array_to_matrix <- function(A){
  p <- dim(A)[1]
  r <- dim(A)[2]
  N <- dim(A)[3]
  
  M <- matrix(NA, N, p*r)
  
  for(k in 1:N){
    M[k,] <- as.vector(A[,,k])
  }
  
  M
}

### Load here the .rData with matrix-variate results ###

### Prepare for parallel computing ###

nThreads <- 28
cluster <- makeCluster(nThreads, type = "SOCK")
registerDoSNOW(cluster)

pb2 <- progress_bar$new(
  format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
  total = nrow(scenarios),
  complete = "=",   # Completion bar character
  incomplete = "-", # Incomplete bar character
  current = ">",    # Current bar character
  width = 100)

progress <- function(n){
  pb2$tick()
} 
opts <- list(progress = progress)
comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

#### Simulation and Fitting ####

res_multi <- foreach(l = 1:nrow(scenarios), .combine = "comb", .multicombine = TRUE, 
                  .init = list(list(), list(), list(), list(), list(), list()), .packages = "CensMFM", .options.snow = opts) %dopar% {
  
  nsims <- 100

  X.or <- resF[[1]][[l]]
  X.cens <- resF[[2]][[l]]
  cc <- resF[[3]][[l]]
  LS <- resF[[4]][[l]]

  X.or.M <- X.cens.M <- cc.M <- LS.M <- LI.M <- res <- vector(mode = "list", length = nsims)
  
  for (j in 1:nsims) {
  
    X.or.M[[j]] <- array_to_matrix(X.or[[j]])
    X.cens.M[[j]] <- array_to_matrix(X.cens[[j]])
    cc.M[[j]] <- array_to_matrix(cc[[j]])
    LS.M[[j]] <- array_to_matrix(LS[[j]])
    LI.M[[j]] <- array(0, dim = dim(LS.M[[j]]))
    LI.M[[j]][cc.M[[j]] == 1] <- X.cens.M[[j]][cc.M[[j]] == 1]
    
    ### Fit multivariate model ###
      
    temp.fit <- vector(mode = "list", length = ntry)
    temp.bic <- numeric(ntry)
    
    for (m in 1:ntry) {
    
      temp.fit[[m]] <- tryCatch(CensMFM::fit.FMMSNC(cc = cc.M[[j]], LI = LI.M[[j]], LS = LS.M[[j]], y = X.or.M[[j]], g = k, get.init = TRUE,
                                                    criteria = TRUE, family = "Normal", error = 0.001, iter.max = 500,
                                                    uni.Gama = FALSE, cal.im = FALSE), error = function(e) {NA})
      
      temp.bic[m] <- tryCatch(temp.fit[[m]][["res"]][["bic"]], error=function(e){NA})
      
    }
    
    wn <- which.min(temp.bic)
    res[[j]]<- tryCatch(temp.fit[[wn]], error = function(e) {NA})
    
   }
    
  list(X.or.M, X.cens.M, cc.M, LS.M, LI.M, res)
  
}  







