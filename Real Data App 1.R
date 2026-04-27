library(progress)
library(foreach)
library(doSNOW)
library(baytrends)
library(dplyr)

data("dataCensored")
dataOR <- unSurvDF(dataCensored)

stat <- unique(dataOR$station)
vars.cens <- c("chla", "din", "nh4", "no23", "po4", "tdn", "tdp", "tn", "tp")
col.nam <- c("AP","B","BP","S")

stat.sel <- 5 # select the station
sel  <- c(2,5,8) # select the variables in vars.cens (must be sorted)
p <- length(sel) # number of variables 
col.sel <- c(1,2,3,4) # select the columns in col.nam (must be sorted)
r <- length(col.sel) # number of columns

vars.sel <- character(length = length(sel))
for (i in 1:length(sel)) {
  
  vars.sel[i] <- paste(vars.cens[sel[i]],"_lo",sep = "")
  
}

#### Data extraction and arrangement ####

filt1 <- filter(dataOR, station==stat[stat.sel])
times1 <- unique(filt1$date)
elim1 <- numeric(length(times1))

for (i in 1:length(times1)) {
  
  temp1 <- filter(filt1, date==times1[i])
  if(nrow(temp1)!=r) # check for each observation having each layer
    elim1[i] <- 1
  
  temp2 <- select(temp1, vars.sel)
  if(all(is.na(temp2))) # check for each observation having not all NA values
    elim1[i] <- 1
  
}

times1<- times1[-which(elim1==1)]
n <- length(times1)

X.cens <- array(NA, dim = c(p,r,n), dimnames = list(c(vars.cens[sel]),
                                                    c(col.nam[col.sel]),
                                                    c(1:n)))
cc <- LS <- array(0, dim = c(p,r,n))

for (i in 1:n) {
  
  temp3   <- filter(filt1, date==times1[i])
  temp4   <- select(temp3, vars.sel)
  
  X.cens[,,i] <- t(temp4)
  
  ## Handling Missing Values (NA) ##
  
  if(any(is.na(X.cens[,,i]))){
    res0 <- which(is.na(X.cens[,,i]))
    cc[,,i][res0] <- 1
    LS[,,i][res0] <- +Inf
    X.cens[,,i][res0] <- -Inf
  }
  
  ## Handling Censoring Values ##
  
  sel2 <- which(colnames(temp3) %in% vars.sel) # positions of the variables in temp3
  
  for (j in 1:p) {
    
    dif <- temp3[,sel2[j]]-temp3[,(sel2[j]+1)]
    if(any(dif!=0, na.rm = T)){
      
      res1 <- which(dif!=0)
      LS[j,res1,i] <- temp3[res1,(sel2[j]+1)]
      cc[j,res1,i] <- 1
      
    }
    
  }
  
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


