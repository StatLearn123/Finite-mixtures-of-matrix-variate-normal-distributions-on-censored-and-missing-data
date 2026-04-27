library(progress)
library(foreach)
library(doSNOW)
library(CopulaCenR)
library(dplyr)

data("AREDS")

#### Data extraction and arrangement ####

p <- 2
r <- 2
n <- 629

X.cens <- array(NA, dim = c(p,r,n), dimnames = list(c("SX","DX"),
                                                    c("V1","V2"),
                                                    c(1:n)))
cc <- LS <- array(0, dim = c(p,r,n))

for (i in 1:n) {
  
  tempV1   <- filter(AREDS, id==i)
  
  if(!is.infinite(tempV1[1, 3]) & !is.infinite(tempV1[1, 4])){ # interval-censored
    X.cens[1,1,i] <- tempV1[1, 3] # lower bound 
    LS[1,1,i] <- tempV1[1, 4] # upper bound 
    cc[1,1,i] <- 1
  }
  if(!is.infinite(tempV1[2, 3]) & !is.infinite(tempV1[2, 4])){ # interval-censored
    X.cens[2,1,i] <- tempV1[2, 3] # lower bound 
    LS[2,1,i] <- tempV1[2, 4] # upper bound 
    cc[2,1,i] <- 1
  }
  
  if(!is.infinite(tempV1[1, 3]) &  is.infinite(tempV1[1, 4])){ # right-censored
    X.cens[1,1,i] <- tempV1[1, 3] # lower bound 
    LS[1,1,i] <- +Inf # upper bound 
    cc[1,1,i] <- 1
  }
  if(!is.infinite(tempV1[2, 3]) &  is.infinite(tempV1[2, 4])){ # right-censored
    X.cens[2,1,i] <- tempV1[2, 3] # lower bound 
    LS[2,1,i] <- +Inf # upper bound 
    cc[2,1,i] <- 1
  }
  
  X.cens[1,2,i] <- tempV1[1,6]
  X.cens[2,2,i] <- tempV1[2,6]
  
}

#### Fitting Part ####

ntry <- nThreads <- 10

### Prepare for parallel computing ###

cluster <- makeCluster(nThreads, type = "SOCK")
registerDoSNOW(cluster)

pb2 <- progress_bar$new(
  format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
  total = ntry,
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

resALL <- foreach(l = 1:ntry, .combine = "comb", .multicombine = TRUE, 
                  .init = list(list()), .options.snow = opts) %dopar% {
                    
    res.fit <- vector(mode = "list", length = 8)
    
    for (j in 1:8) {
      
      init <- EM.MatrixNC.init.mix(dados = X.cens, cc = cc, LS = LS, g = j, seed = l)
      res.fit[[j]]<- tryCatch(EM.MatrixNC.fit.mix(dados = X.cens, cc = cc, LS = LS, init = init, precision=0.001, MaxIter=500, plotll = F),error=function(e){NA})

    }
                    
    list(res.fit)
    
}

