library(foreach)
library(doSNOW)
library(progress)
# library(mclust) must be installed

#### Definition of the simulation objects ####

n   <- c(100, 200, 400)
bal <- c(1,2) # 1 balance, 2 unbalance
ov  <- c(1,2) # 1 close, 2 far
ty  <- c(1,2,3) # 1 only missing, 2 only int. censored, 3 both types
w   <- c(0.05,0.15,0.30) # proportion of missing/cens values
k <- 2
ntry <- 10

scenarios <-  expand.grid(list(n,bal,ov,ty,w))
colnames(scenarios) <- c("size","bal","overlap", "type","prop")
clo.scen <- which(scenarios$overlap==1)
far.scen <- which(scenarios$overlap==2)

### Prepare for parallel computing ###

nThreads <- 28 # Select here the number of cores for parallel computing
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

# load the "Parameters.RData" file here before running the simulation #

resF <- foreach(l = 1:nrow(scenarios), .combine = "comb", .multicombine = TRUE, 
               .init = list(list(), list(), list(), list(), list(), list(), list()), .options.snow = opts) %dopar% {
  
  nsims <- 100
  disc <- 0
  
  X.or <- X.cens <- cc <- LS <- res <- vector(mode = "list", length = nsims)
  class <- matrix(NA, scenarios[l,1], nsims)
  
  for (j in 1:nsims) {
  
    check1 <- check2 <- 0
    
    while (check2 < 1) {
      
      while (check1 < 1) {
        
        temp.par <- param[[scenarios[l,2]]][[1]][[scenarios[l,3]]]   
        
        dat <- r_MVN_MX(num = scenarios[l,1], pix=temp.par[[1]], M=temp.par[[2]], U=temp.par[[3]], V=temp.par[[4]])
        
        p <- dim(dat$Y)[1]
        r <- dim(dat$Y)[2]
        
        iniMV.true <- list()
        iniMV.true[["prior"]]<- temp.par[[1]]
        iniMV.true[["M"]]<- temp.par[[2]]
        iniMV.true[["U"]]<- temp.par[[3]]
        iniMV.true[["V"]]<- temp.par[[4]]
        iniMV.true[["class"]]<- dat$obs.class
        
        fitMV.true <- tryCatch(fit_MVN_MX(Y=dat$Y, k = k, init.par = iniMV.true) , error = function(e) {NA})  
        ari.true <- tryCatch(round(mclust::adjustedRandIndex(dat$obs.class,fitMV.true$class), digits = 2), error = function(e) {0}) 
        
        ### Check for ARI Interval ###
        
        if(l %in% clo.scen){
          if(ari.true >= 0.80 & ari.true <=0.85){
            check1 <- 1
          }else{
            disc <- disc +1
          }
          
        }else{
          if(ari.true >= 0.95 & ari.true <=1){
            check1 <- 1
          }else{
            disc <- disc +1
          }
        }
        
      }
      
      X.or[[j]] <- dat$Y
      class[,j] <- dat$obs.class

      groups <- unique(class[,j])
      perc <- scenarios[l,5]
      
      cc[[j]] <- LS[[j]] <- array(0, dim = dim(X.or[[j]]))
      X.cens[[j]] <- X.or[[j]]
      
      for(g in groups){
        
        ind_g <- which(class[,j] == g)
        Xg <- X.or[[j]][,,ind_g]
        
        tot_g <- length(Xg)
        n_tot <- round(perc * tot_g)
        
        idx_local <- sample(tot_g, n_tot)
        
        coord_local <- arrayInd(idx_local, dim(Xg))
        
        k_global <- ind_g[coord_local[,3]]
        
        coord_global <- cbind(coord_local[,1], coord_local[,2], k_global)
        
        if(scenarios[l,4]==1){   # only MISSING
          
          X.cens[[j]][coord_global] <- -Inf
          LS[[j]][coord_global] <- +Inf
          cc[[j]][coord_global] <- 1
          
        }
        
        if(scenarios[l,4]==2){   # only INTERVAL CENSORING
          
          sd_mat <- apply(Xg, c(1,2), sd)
          
          sd_vals <- sd_mat[cbind(coord_local[,1], coord_local[,2])]
          
          v <- X.or[[j]][coord_global]
          
          u1 <- runif(length(v), 0, sd_vals)
          u2 <- runif(length(v), 0, sd_vals)
          
          a <- v - u1
          b <- v + u2
          
          X.cens[[j]][coord_global] <- a
          LS[[j]][coord_global] <- b
          cc[[j]][coord_global] <- 1
          
        }
        
        if(scenarios[l,4]==3){   # 50% MISSING 50% INTERVAL
          
          n_miss <- floor(n_tot/2)
          n_int  <- n_tot - n_miss
          
          miss_idx <- 1:n_miss
          int_idx  <- (n_miss+1):n_tot
          
          coord_miss <- coord_global[miss_idx,,drop=FALSE]
          coord_int  <- coord_global[int_idx,,drop=FALSE]
          
          # ---- MISSING
          
          X.cens[[j]][coord_miss] <- -Inf
          LS[[j]][coord_miss] <- +Inf
          cc[[j]][coord_miss] <- 1
          
          # ---- INTERVAL
          
          sd_mat <- apply(Xg, c(1,2), sd)
          
          coord_local_int <- coord_local[int_idx,,drop=FALSE]
          
          sd_vals <- sd_mat[cbind(coord_local_int[,1], coord_local_int[,2])]
          
          v <- X.or[[j]][coord_int]
          
          u1 <- runif(length(v), 0, sd_vals)
          u2 <- runif(length(v), 0, sd_vals)
          
          a <- v - u1
          b <- v + u2
          
          X.cens[[j]][coord_int] <- a
          LS[[j]][coord_int] <- b
          cc[[j]][coord_int] <- 1
        }
        
      }
      
      ### Fit our model ###
      
      temp.fit <- vector(mode = "list", length = ntry)
      temp.bic <- numeric(ntry)
      
      for (m in 1:ntry) {
        
        init <- EM.MatrixNC.init.mix(dados = X.cens[[j]], cc = cc[[j]], LS = LS[[j]], g = k, seed = m)
        
        temp.fit[[m]] <- tryCatch(EM.MatrixNC.fit.mix(dados = X.cens[[j]], cc = cc[[j]], LS = LS[[j]], init = init, precision=0.001, MaxIter=500, plotll = F),error=function(e){NA})
        temp.bic[m] <- tryCatch(temp.fit[[m]][["BIC"]], error=function(e){NA})

      }
      
      wn <- which.min(temp.bic)
      res[[j]]<- tryCatch(temp.fit[[wn]], error = function(e) {NA})
      
      ### Check for fitting/computational issues ###
      
      if(all(!is.na(res[[j]]))){
        check2 <- 1
      }
      
    }

   }
    
  list(X.or, X.cens, cc, LS, res, class, disc)
  
}  







