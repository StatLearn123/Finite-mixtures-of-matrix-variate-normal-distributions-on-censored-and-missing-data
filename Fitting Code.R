# library(matrixNormal) #must be installed
# library(mnormt) #must be installed   
# library(MomTrunc) #must be installed
# library(Matrix) #must be installed
# library(withr) #must be installed
# library(mclust) #must be installed
# library(tensor) #must be installed
# library(tidyr) #must be installed
# library(utils) #must be installed
# library(data.table) #must be installed
# library(LaplacesDemon) #must be installed

#### Utils ECM functions ####

logLikCens <- function(cc, LS, dados, muM, SigmaM, PsiM){
  
  p<-dim(dados)[1]
  q<-dim(dados)[2]
  n<-dim(dados)[3]
  ver <- matrix(0,n,1)
  
  Vari<- kronecker(PsiM,SigmaM)
  
  for (j in 1:n){
    
    cc1<-as.vector(matrixNormal::vec(cc[,,j]))
    LS1<-as.vector(matrixNormal::vec(LS[,,j]))
    y1<-as.vector(matrixNormal::vec(dados[,,j]))
    mu1<- as.vector(matrixNormal::vec(muM))
    
    if(sum(cc1)==0){
      
      ver[j]<-matrixNormal::dmatnorm(dados[,,j], muM, SigmaM, PsiM, tol = .Machine$double.eps^0.5,log = FALSE)
    }
    
    if(sum(cc1)>=1){
      if(sum(cc1)==p*q){
        ver[j]<- matrixNormal::pmatnorm(Lower = dados[,,j], Upper = LS[,,j], muM, SigmaM, PsiM, tol = .Machine$double.eps^0.5)
        
      }
      else{
        muc <- mu1[cc1==1]+Vari[cc1==1,cc1==0]%*%solve(Vari[cc1==0,cc1==0],tol = .Machine$double.eps^0.5)%*%(y1[cc1==0]-mu1[cc1==0])
        Sc  <- Vari[cc1==1,cc1==1]-Vari[cc1==1,cc1==0]%*%solve(Vari[cc1==0,cc1==0],tol = .Machine$double.eps^0.5)%*%Vari[cc1==0,cc1==1]
        Sc  <- (Sc+t(Sc))/2
        
        auxM<- Vari[cc1==0,cc1==0]
        auxM<-(auxM+t(auxM))/2
        
        auxcdf<-mnormt::sadmvn(lower=as.vector(y1[cc1==1]), upper=as.vector(LS1[cc1==1]), as.vector(muc), Sc, abseps=.Machine$double.eps^0.5)
        auxpdf<-mnormt::dmnorm(as.vector(y1[cc1==0]),as.vector(mu1[cc1==0]),auxM)
        ver[j]<- auxcdf*auxpdf
      }
    }
  }

  return(ver)
  
}

d.mixedMatNCens <- function(cc, LS, dados, pi1, muM, SigmaM, PsiM){

  g    <- length(pi1)
  dens <- 0
  for (j in 1:g) dens <- dens + pi1[j]*logLikCens(cc, LS, dados, muM[[j]], SigmaM[[j]], PsiM[[j]])
  return(dens)
}

somaL3<-function(L1,Sigma,Psi){
  n<-ncol(L1)
  p<-ncol(Sigma)
  q<-ncol(Psi)
  suma2<-matrix(0,q,q)
  suma3<-matrix(0,p,p)
  for(j in 1:n){
    aux<- matrix(L1[,j],p,q)
    suma2<-suma2+t(aux)%*%solve(Sigma)%*%aux
    suma3<-suma3+aux%*%solve(Psi)%*%t(aux)
  }
  return(list(sPsi=suma2,sSigma=suma3))
}

#### ECM algorithm for interval censored and missing data ####

EM.MatrixNC.init.mix<- function(dados, cc, LS, g, seed = NULL){
  
  p <- dim(dados)[1]
  repl_neg_inf <- function(arr) {
    for (i in 1:dim(arr)[1]) {    
      for (j in 1:dim(arr)[2]) {  
        values <- arr[i, j, ]
        
        mean_val <- mean(values[values != -Inf], na.rm = TRUE)
        
        arr[i, j, values == -Inf] <- mean_val
      }
    }
    return(arr)
  }
  
  dados1 <- dados
  dados2 <- repl_neg_inf(dados)
  
  if(is.null(seed)){
    seed <- p
  }else{
    seed <- seed
  }
  
  ini0 <- init_MVN_MX(dados2, g, nstartR = 25, maxit = 1, seed = seed)
  inis <- fit_MVN_MX(dados2, g, init.par = ini0, tol = 0.01)
  
  return(init = inis)
  
}

EM.MatrixNC.fit.mix<-function(dados, cc, LS, init = NULL, precision=0.0000001, MaxIter=50, plotll = TRUE){
  
  ## dados: if censoring it should contain the lower limit, if missing use -Inf (array)
  ## cc: indicator matrix 1: censoring 0: observed    (array)
  ## LS: superior Limit of the interval censoring, if missing then use +Inf   (array)
  ## init : initialization object 
  ## precision:  is the precision used in the convergence
  ## MaxIter: Maximum number of iterations
  ## plotll: log-lik plot
  
  ptm <- proc.time()
  
  p <- dim(dados)[1]
  q <- dim(dados)[2]
  n <- dim(dados)[3]
  
  dados1 <- dados
  
  pii<- init$prior
  mu<-init$M
  Sigma<-init$sigmaU
  Psi<-init$sigmaV
  g <- init$k
  
  mu<-lapply(seq(dim(mu)[3]), function(x) mu[ , , x])
  Sigma<-lapply(seq(dim(Sigma)[3]), function(x) Sigma[ , , x])
  Psi<-lapply(seq(dim(Psi)[3]), function(x) Psi[ , , x])
  
  loglik<-numeric()
  criterio<-1
  count<-0
  
  while(criterio > precision){
    
    cat("*")
    count <- count + 1
    tal   <- matrix(0, n, g)
    
    for (j in 1:g){
      
      # sort pii
      if( sqrt(sum((order(pii)-g:1)*(order(pii)-g:1))) != 0 ){
        musi  <- mu
        Sigsi <- Sigma
        Psisi <- Psi
        for(l in 1:g){
          mu[[l]]    <- musi[[(order(pii, decreasing = TRUE)[l])]]
          Sigma[[l]] <- Sigsi[[(order(pii, decreasing = TRUE)[l])]]
          Psi[[l]] <- Psisi[[(order(pii, decreasing = TRUE)[l])]]
        }
        pii <- pii[order(pii, decreasing = TRUE)]
      }
      
      suma0<-matrix(0,p,q)
      suma2<-matrix(0,q,q)
      suma3<-matrix(0,p,p)
      
      d1 <- logLikCens(cc, LS, dados, mu[[j]], Sigma[[j]], Psi[[j]])
      if(length(which(d1 == 0)) > 0) d1[which(d1 == 0)] <- .Machine$double.xmin
      
      d2 <- d.mixedMatNCens(cc, LS, dados, pii, mu, Sigma, Psi)
      
      if(length(which(d2 == 0)) > 0) d2[which(d2 == 0)] <- .Machine$double.xmin
      
      tal[,j] <- d1*pii[j] / d2
      
      pii[j] <- (1/n)*sum(tal[,j])
      
      mus    <- mu[[j]]
      Sigmas <- Sigma[[j]]
      Psis   <- Psi[[j]]
      Varis  <- kronecker(Psis,Sigmas)
      
      for (i in 1:n){
        cc1<-as.vector(matrixNormal::vec(cc[,,i]))
        LS1<-as.vector(matrixNormal::vec(LS[,,i]))
        y1<-as.vector(matrixNormal::vec(dados[,,i]))
        mu1<- as.vector(matrixNormal::vec(mus))
        
        if(sum(cc1)==0){
          omega<- matrix(y1,p,q)
          suma0<- suma0 + tal[i,j]*omega
        }
        
        if(sum(cc1) >= 1){
          
          if(sum(cc1) == p*q) {
            muc  <- mu1
            Sc   <- Varis
            aux  <- MomTrunc::meanvarTMD(y1,LS1,muc,Sc,dist="normal")
            tuy  <- aux$mean
            dados1[,,i]<- matrix(tuy,p,q)
            
          } else {
            muc <- mu1[cc1==1]+Varis[cc1==1,cc1==0]%*%solve(Varis[cc1==0,cc1==0])%*%(y1[cc1==0]-mu1[cc1==0])
            Sc  <- Varis[cc1==1,cc1==1]-Varis[cc1==1,cc1==0]%*%solve(Varis[cc1==0,cc1==0])%*%Varis[cc1==0,cc1==1]
            Sc  <- (Sc+t(Sc))/2
            aux <- MomTrunc::meanvarTMD(y1[cc1==1],LS1[cc1==1],muc,Sc,dist="normal")
            tuy <- matrix(y1,p*q,1)
            tuy[cc1==1] <- aux$mean
            dados1[,,i]<- matrix(tuy,p,q)
          }
          
          suma0<-suma0+tal[i,j]*matrix(tuy,p,q)
          
        }
        
      }
      
      mu[[j]]<- suma0/sum(tal[,j])
      mus    <- mu[[j]]
      
      for (i in 1:n){
        cc1<-as.vector(matrixNormal::vec(cc[,,i]))
        LS1<-as.vector(matrixNormal::vec(LS[,,i]))
        y1<-as.vector(matrixNormal::vec(dados[,,i]))
        mu1<- as.vector(matrixNormal::vec(mus))
        
        if(sum(cc1)==0){
          omega<- matrix(y1,p,q)
          suma2<-suma2+tal[i,j]*(t(omega-mu[[j]])%*%solve(Sigma[[j]])%*%(omega - mu[[j]]))
          
        }
        
        if(sum(cc1) >= 1){
          
          if(sum(cc1) == p*q) {
            muc  <- mu1
            Sc   <- Varis
            aux  <- MomTrunc::meanvarTMD(y1,LS1,muc,Sc,dist="normal")
            tuy  <- aux$mean
            tuyy <- aux$EYY
            dados1[,,j]<- matrix(tuy,p,q)
            
          } else {
            muc <- mu1[cc1==1]+Varis[cc1==1,cc1==0]%*%solve(Varis[cc1==0,cc1==0])%*%(y1[cc1==0]-mu1[cc1==0])
            Sc  <- Varis[cc1==1,cc1==1]-Varis[cc1==1,cc1==0]%*%solve(Varis[cc1==0,cc1==0])%*%Varis[cc1==0,cc1==1]
            Sc  <- (Sc+t(Sc))/2
            aux <- MomTrunc::meanvarTMD(y1[cc1==1],LS1[cc1==1],muc,Sc,dist="normal")
            tuy <- matrix(y1,p*q,1)
            tuy[cc1==1] <- aux$mean
            tuyy <- matrix(0,p*q,p*q)
            tuyy[cc1==1,cc1==1] <- aux$varcov
            tuyy <- tuyy+tuy%*%t(tuy)
            dados1[,,i]<- matrix(tuy,p,q)
          }
          
          omega1<- tuyy-(tuy)%*%t(mu1)-(mu1)%*%t(tuy)+mu1%*%t(mu1)
          omega<- as.matrix(Matrix::nearPD(omega1)$mat)
          omega<- (omega+t(omega))/2
          L1<- t(Matrix::chol(omega))
          aux<- somaL3(L1,Sigma[[j]],Psi[[j]])
          suma2<-suma2+tal[i,j]*aux$sPsi
          
        }
        
      }
      Psiaux<-suma2
      Psi[[j]] <- Psiaux/det(Psiaux)^(1/q)
      Psi[[j]] <- (Psi[[j]]+t(Psi[[j]]) )/2
      
      Varis<- kronecker(Psi[[j]],Sigma[[j]])
      
      for (i in 1:n){
        cc1<-as.vector(matrixNormal::vec(cc[,,i]))
        LS1<-as.vector(matrixNormal::vec(LS[,,i]))
        y1<-as.vector(matrixNormal::vec(dados[,,i]))
        mu1<- as.vector(matrixNormal::vec(mus))
        
        if(sum(cc1)==0){
          omega<- matrix(y1,p,q)
          suma3<-suma3+tal[i,j]*((omega- mu[[j]])%*%solve(Psi[[j]])%*%t(omega- mu[[j]]))
          
        }
        
        if(sum(cc1)>= 1){
          
          if(sum(cc1) == p*q) {
            muc  <- mu1
            Sc   <- Varis
            aux  <- MomTrunc::meanvarTMD(y1,LS1,muc,Sc,dist="normal")
            tuy  <- aux$mean
            tuyy <- aux$EYY
            dados1[,,i]<- matrix(tuy,p,q)
            
          } else {
            muc <- mu1[cc1==1]+Varis[cc1==1,cc1==0]%*%solve(Varis[cc1==0,cc1==0])%*%(y1[cc1==0]-mu1[cc1==0])
            Sc  <- Varis[cc1==1,cc1==1]-Varis[cc1==1,cc1==0]%*%solve(Varis[cc1==0,cc1==0])%*%Varis[cc1==0,cc1==1]
            Sc  <- (Sc+t(Sc))/2
            aux <- MomTrunc::meanvarTMD(y1[cc1==1],LS1[cc1==1],muc,Sc,dist="normal")
            tuy <- matrix(y1,p*q,1)
            tuy[cc1==1] <- aux$mean
            tuyy <- matrix(0,p*q,p*q)
            tuyy[cc1==1,cc1==1] <- aux$varcov
            tuyy <- tuyy+tuy%*%t(tuy)
            dados1[,,i]<- matrix(tuy,p,q)
          }
          
          omega1<- tuyy-(tuy)%*%t(mu1)-(mu1)%*%t(tuy)+mu1%*%t(mu1)
          omega<- as.matrix(Matrix::nearPD(omega1)$mat)
          omega<- (omega+t(omega))/2
          L1<- t(Matrix::chol(omega))
          aux<- somaL3(L1,Sigma[[j]],Psi[[j]])
          suma3<-suma3+tal[i,j]*aux$sSigma
          
        }
        
      }
      
      Sigma[[j]] <- suma3/(sum(tal[,j])*(q))
      Sigma[[j]] <- (Sigma[[j]]+t(Sigma[[j]]) )/2
      
      Varis  <- kronecker(Psi[[j]],Sigma[[j]])
      pii[g] <- 1 - (sum(pii) - pii[g])
      
    }
    
    loglik[count]<- sum(log(d.mixedMatNCens(cc, LS, dados, pii, mu, Sigma, Psi)))#logLikCens(cc, LS, dados, mu, Sigma, Psi)
    if (count>2){
      at <- (loglik[count]-loglik[count-1])/(loglik[count-1]-loglik[count-2])
      criterio <- (loglik[count]-loglik[count-1])/(1-at)
      if((loglik[count]-loglik[count-1])<0)
        warning("The complete-data loglikelihood values are not monotonic")
    }
    
    if (count==MaxIter){
      criterio <- 0.000000000000001
    }
    
  }
  
  if(plotll == TRUE){
    plot(1:count,loglik,type="l",xlab="Iterations",ylab="Complete-data log-likelihood values")
  }
  
  npar <- g*((p*q)+(p*(p+1)/2)+(q*(q+1)/2)-1) + (g-1)
  BIC  <- -2*loglik[count] + npar*log(n)    # to be minimized
  
  # Computing the Kronecker products among sow and column covariance matrices
  
  CovMatKro <- lapply(seq_along(Psi), function(i) kronecker(Psi[[i]], Sigma[[i]]))
  
  # Classification #
  
  if(g==1){
    classification<- rep(1,n)
  }else{
    colnames(tal)<-c(1:g)
    classification <- as.numeric(colnames(tal)[max.col(tal,ties.method="first")])
  }
  
  ptm2 <- proc.time() - ptm
  time <- ptm2[3]
  
  obj.out <- list(pii = pii, 
                  mu = mu, 
                  Sigma = Sigma, 
                  Psi = Psi, 
                  CovMatKro = CovMatKro,
                  dadosPred = dados1, 
                  loglik = loglik[count], 
                  class = classification,
                  BIC = BIC, 
                  iter = count,
                  time = time)
  
  class(obj.out) <- "EM.MatrixNC.mix"
  return(obj.out)
  
}

#### ECM algorithm for matrix-variate normal mixtures ####

init_MVN_MX <- function(Y, k, nstartR = 50, maxit = 5, seed = NULL) {
  # Dimensions
  
  num <- dim(Y)[3] # sample size
  p <- nrow(Y) # rows of Y ;
  r <- ncol(Y) # columns of Y;
  
  # Create some objects
  
  prior <- matrix(NA, nstartR, k)
  M <- array(NA, dim = c(p, r, k, nstartR))
  sigmaU <- array(NA, dim = c(p, p, k, nstartR))
  sigmaV <- array(NA, dim = c(r, r, k, nstartR))
  
  WR <- array(NA, dim = c(p, p, k))
  WC <- array(NA, dim = c(r, r, k))
  
  sel <- array(NA, dim = c(p, r, k))
  dens <- array(NA, c(num, k), dimnames = list(1:(num), paste("comp.", 1:k, sep = "")))
  llk <- rep(NA, nstartR)
  
  ## Random initialization ##
  
  eu <- matrix(0, nrow = num, ncol = k)
  classy <- matrix(0, nrow = num, ncol = nstartR)
  rand.start <- matrix(0, nstartR, k)
  
  if(is.null(seed)){
    seed <- p
  }else{
    seed <- seed
  }
  withr::with_seed(seed, for (i in 1:nstartR) {
    rand.start[i, ] <- sample(c(1:num), k)
  })

  for (l in 1:nstartR) {
    skip_to_next <- FALSE
    
    tryCatch(
      {
        ### part 0 ###
        
        sec <- rand.start[l, ]
        
        for (j in 1:k) {
          sel[, , j] <- Y[, , sec[j]]
        }
        
        for (j in 1:k) {
          for (i in 1:num) {
            eu[i, j] <- norm((Y[, , i] - sel[, , j]), type = "F")
          }
        }
        
        for (i in 1:num) {
          classy[i, l] <- which.min(eu[i, ])
        }
        
        z <- mclust::unmap(classy[, l])
        
        ### part 1 ###
        
        for (j in 1:k) {
          M[, , j, l] <- rowSums(Y * z[, j][slice.index(Y, 3)], dims = 2) / sum(z[, j])
          
          Mtemp <- array(M[, , j, l], dim = c(p, r, num))
          
          WR[, , j] <- tensor::tensor(aperm(tensor::tensor((Y - Mtemp) * z[, j][slice.index(Y - Mtemp, 3)], diag(r), 2, 1), c(1, 3, 2)), aperm((Y - Mtemp), c(2, 1, 3)), c(2, 3), c(1, 3))
          
          sigmaU[, , j, l] <- WR[, , j] / (r * sum(z[, j]))
          
          WC[, , j] <- tensor::tensor(aperm(tensor::tensor(aperm(Y - Mtemp, c(2, 1, 3)) * z[, j][slice.index(aperm(Y - Mtemp, c(2, 1, 3)), 3)], solve(sigmaU[, , j, l]), 2, 1), c(1, 3, 2)), (Y - Mtemp), c(2, 3), c(1, 3))
          
          sigmaV[, , j, l] <- WC[, , j] / (p * sum(z[, j]))
        }
        
        if (k == 1) {
          prior[l, ] <- 1
        } else {
          prior[l, ] <- colMeans(z)
        }
        
        m.iter <- 0
        loglik.new <- NULL
        ll <- NULL
        
        ### part 2 ###
        
        tryCatch(while (m.iter < maxit) {
          m.iter <- m.iter + 1
          
          ### E - STEP ###
          
          for (j in 1:k) {
            dens[, j] <- as.vector(dMVnorm(X = Y, M = M[, , j, l], U = sigmaU[, , j, l], V = sigmaV[, , j, l]))
          }
          
          numerator <- matrix(rep(prior[l, ], num), num, k, byrow = TRUE) * dens
          mixt.dens <- rowSums(numerator)
          
          post <- numerator / mixt.dens
          
          ### CM - STEPS ###
          
          for (j in 1:k) {
            M[, , j, l] <- rowSums(Y * post[, j][slice.index(Y, 3)], dims = 2) / sum(post[, j])
            
            Mtemp <- array(M[, , j, l], dim = c(p, r, num))
            
            WR[, , j] <- tensor::tensor(aperm(tensor::tensor((Y - Mtemp) * post[, j][slice.index(Y - Mtemp, 3)], solve(sigmaV[, , j, l]), 2, 1), c(1, 3, 2)), aperm((Y - Mtemp), c(2, 1, 3)), c(2, 3), c(1, 3))
            
            sigmaU[, , j, l] <- WR[, , j] / (r * sum(post[, j]))
            
            WC[, , j] <- tensor::tensor(aperm(tensor::tensor(aperm(Y - Mtemp, c(2, 1, 3)) * post[, j][slice.index(aperm(Y - Mtemp, c(2, 1, 3)), 3)], solve(sigmaU[, , j, l]), 2, 1), c(1, 3, 2)), (Y - Mtemp), c(2, 3), c(1, 3))
            
            sigmaV[, , j, l] <- WC[, , j] / (p * sum(post[, j]))
          }
          
          if (k == 1) {
            prior[l, ] <- 1
          } else {
            prior[l, ] <- colMeans(post)
          }
          
          for (j in 1:k) {
            dens[, j] <- as.vector(dMVnorm(X = Y, M = M[, , j, l], U = sigmaU[, , j, l], V = sigmaV[, , j, l]))
          }
          
          # mixture density
          
          numerator <- matrix(rep(prior[l, ], num), num, k, byrow = TRUE) * dens
          mixt.dens <- rowSums(numerator)
          loglik.new <- sum(log(mixt.dens))
          ll <- c(ll, loglik.new)
        }, error = function(e) {
          llk[l] <- loglik.new
        })
        
        llk[l] <- loglik.new
        
        if (any(prior[l, ] <= 0.05)) {
          llk[l] <- NA
        }
      },
      error = function(e) {
        skip_to_next <<- TRUE
      }
    )
    
    if (skip_to_next) {
      next
    }
  }
  
  df <- data.frame(llk = llk, pos = c(1:nstartR))
  df <- tidyr::drop_na(df)
  df <- df[!is.infinite(rowSums(df)), ]
  
  bestR <- utils::head(data.table::setorderv(df, cols = "llk", order = -1), n = 1)$pos
  
  return(list(prior = prior[bestR, ], M = M[, , , bestR], U = sigmaU[, , , bestR], V = sigmaV[, , , bestR]))
}

fit_MVN_MX <- function(Y, k, init.par = NULL, tol = 0.001, maxit = 500) {
  ptm <- proc.time()
  
  dec <- function(x) {
    if ((x %% 1) != 0) {
      nchar(strsplit(sub("0+$", "", as.character(x)), ".", fixed = TRUE)[[1]][[2]])
    } else {
      return(0)
    }
  }
  
  # Dimensions
  
  num <- dim(Y)[3] # sample size
  p <- nrow(Y) # rows of Y ;
  r <- ncol(Y) # columns of Y;
  
  ## Objects related to the two covariance matrices
  
  WR <- array(0, dim = c(p, p, k))
  WC <- array(0, dim = c(r, r, k))
  
  ## Other objects
  
  post <- dens <- array(0, c(num, k), dimnames = list(1:(num), paste("comp.", 1:k, sep = "")))
  check <- 0
  loglik.old <- -Inf
  loglik.new <- NULL
  ll <- NULL
  mark <- 0
  m.iter <- 0
  cnt <- 0
  
  ### Algorithm ###
  
  if (k == 1) {
    M <- array(init.par$M, dim = c(p, r, k))
    sigmaU <- array(init.par$U, dim = c(p, p, k))
    sigmaV <- array(init.par$V, dim = c(r, r, k))
  } else {
    M <- init.par$M
    sigmaU <- init.par$U
    sigmaV <- init.par$V
  }
  
  prior <- init.par$prior
  
  ### Estimation ###
  
  tryCatch(while (check < 1) {
    m.iter <- m.iter + 1
    
    ### E - STEP ###
    
    for (j in 1:k) {
      dens[, j] <- as.vector(dMVnorm(X = Y, M = M[, , j], U = sigmaU[, , j], V = sigmaV[, , j]))
    }
    
    numerator <- matrix(rep(prior, num), num, k, byrow = TRUE) * dens
    mixt.dens <- rowSums(numerator)
    
    post <- numerator / mixt.dens
    
    ### CM - STEPS ###
    
    for (j in 1:k) {
      M[, , j] <- rowSums(Y * post[, j][slice.index(Y, 3)], dims = 2) / sum(post[, j])
      
      Mtemp <- array(M[, , j], dim = c(p, r, num))
      
      WR[, , j] <- tensor::tensor(aperm(tensor::tensor((Y - Mtemp) * post[, j][slice.index(Y - Mtemp, 3)], solve(sigmaV[, , j]), 2, 1), c(1, 3, 2)), aperm((Y - Mtemp), c(2, 1, 3)), c(2, 3), c(1, 3))
      
      sigmaU[, , j] <- WR[, , j] / (r * sum(post[, j]))
      
      WC[, , j] <- tensor::tensor(aperm(tensor::tensor(aperm(Y - Mtemp, c(2, 1, 3)) * post[, j][slice.index(aperm(Y - Mtemp, c(2, 1, 3)), 3)], solve(sigmaU[, , j]), 2, 1), c(1, 3, 2)), (Y - Mtemp), c(2, 3), c(1, 3))
      
      sigmaV[, , j] <- WC[, , j] / (p * sum(post[, j]))
    }
    
    if (k == 1) {
      prior <- 1
    } else {
      prior <- colMeans(post)
    }
    
    for (j in 1:k) {
      dens[, j] <- as.vector(dMVnorm(X = Y, M = M[, , j], U = sigmaU[, , j], V = sigmaV[, , j]))
    }
    
    # mixture density
    
    numerator <- matrix(rep(prior, num), num, k, byrow = TRUE) * dens
    mixt.dens <- rowSums(numerator)
    loglik.new <- sum(log(mixt.dens))
    ll <- c(ll, loglik.new)
    
    if ((round(loglik.new, dec(tol)) - round(loglik.old, dec(tol))) < tol) {
      check <- 1
    }
    
    if (m.iter > 1) {
      if (round(loglik.new, dec(tol)) < round(loglik.old, dec(tol))) {
        mark <- 1
      }
    }
    
    if (m.iter == maxit) {
      check <- 1
    }
    
    loglik.old <- loglik.new
  }, error = function(e) {
    cnt <<- 1
  })
  
  #### Output ####
  
  if (cnt == 0) {
    if (k == 1) {
      classification <- rep(1, num)
    } else {
      colnames(post) <- c(1:k)
      classification <- as.numeric(colnames(post)[max.col(post, ties.method = "first")])
    }
    
    # Number of parameters
    
    meanpar <- (p * r) * k
    rowpar <- k * p * (p + 1) / 2
    colpar <- k * (r * ((r + 1) / 2) - 1)
    
    weights <- k - 1
    
    npar <- meanpar + rowpar + colpar + weights
    
    ptm2 <- proc.time() - ptm
    time <- ptm2[3]
    
    if (any(round(prior, digits = 2) <= 0.05)) {
      sp <- 1
    } else {
      sp <- 0
    }
    
    return(list(
      num = num, k = k, prior = prior, M = M, sigmaU = sigmaU, sigmaV = sigmaV,
      loglik = loglik.new, mark = mark, check = check, npar = npar, time = time,
      z = post, class = classification, sp = sp
    ))
  } else {
    return(NA)
  }
} # fit

dMVnorm <- function(X, M, U, V) {
  num <- dim(X)[3] # sample size
  p <- nrow(X) # rows of X
  r <- ncol(X) # columns of X
  
  if (is.na(num)) {
    X <- as.matrix(X)
    delta <- tr(solve(U) %*% (X - M) %*% solve(V) %*% t(X - M))
  } else {
    delta <- sapply(1:num, function(i) tr(solve(U) %*% (X[, , i] - M) %*% solve(V) %*% t(X[, , i] - M)))
  }
  
  pdf <- (2 * pi)^(-(p * r) / 2) * det(U)^(-r / 2) * det(V)^(-p / 2) * exp(-1 / 2 * delta)
  
  return(pdf)
}

tr <- function(x) {
  return(sum(diag(x)))
}

r_MVN_MX <- function (num, pix, M, U, V) {
  row.Y <- dim(U)[1]
  col.Y <- dim(V)[2]
  k <- dim(V)[3]
  Y <- array(NA, dim = c(row.Y, col.Y, num))
  temp <- sample(1:k, size = num, replace = TRUE, prob = pix)
  numE <- numeric(k)
  for (i in 1:k) {
    numE[i] <- length(temp[temp == i])
  }
  for (i in 1:num) {
    Y[, , i] <- LaplacesDemon::rmatrixnorm(M = M[, , temp[i]], 
                                           U = U[, , temp[i]], V = V[, , temp[i]])
  }
  return(list(Y = Y, obs.class = temp))
}


