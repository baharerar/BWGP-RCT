
###### ---------------- Supplemental functions to BWGP.RCT --------------- ######

# Author: Bahar Erar

#################################################################################

require(compiler)

##################################################################################################
#The welcome function that will appear every time that your run the program
welcome <- function()
{
  cat("\n");
  cat("---------------------------------------------------------------------\n");
  cat("                                                                   \n");
  cat(" BWGP.RCT                                                          \n");
  cat("---------------------------------------------------------------------\n");
  cat(" Bayesian Whole Genome Prediction in                               \n");
  cat(" Randomized Control Trials                                         \n");
  cat("                                                                   \n");
  cat(" Bahar Erar                                                        \n");
  cat(" Brown University, 2016                                            \n");
  cat("                                                                   \n");
  cat(" Current version date: 3/1/2016                                   \n");                                                                   #\n");
  cat("---------------------------------------------------------------------\n");
  cat("\n");
}
##################################################################################################

##################################################################################################
#The finish function that will appear every time the program finishes running
# finish <- function(tot.time, nIter, n, J)
finish <- function(tot.time, nIter)
{
  cat("\n");
  cat("---------------------------------------------------------------------\n")
  cat("                                                                     \n");
  cat(" BGWP.RCT complete.                                                  \n")
  cat("---------------------------------------------------------------------\n")
#   cat(" Data size:                                                          \n");
#   cat(c(paste(c("    N = ", "\n    P = "), 
#               c(n,  J), sep = "")), 
#       "\n")  
  cat("                                                                     \n");
  cat(c(paste(c(" Total iteration time elapsed = ", "s\n Average Time/Iter= "), 
              c(round(tot.time, 2),  round(tot.time/nIter, 3)), sep = "")), 
      "s\n")  
  cat("---------------------------------------------------------------------\n")
  cat("\n");
}
##################################################################################################


## Fixed Effects ##################################################################
#Function for initializing regression coefficients for Fixed effects.
set.list.Fixed <- function(X.fixed, saveAt, e, saveChain)
{
  
  FL <- list(X=X.fixed, p=ncol(X.fixed), colNames=colnames(X.fixed))
  
  #Sum of the squares for each term
  C <- solve(t(FL$X) %*% FL$X) # (W'W)^(-1)
  CtX <- C %*% t(FL$X)
    
  #Objects for saving posterior means from MCMC
  FL$b <- rep(0,FL$p)
  names(FL$b) <- FL$colNames
  FL$post_b <- rep(0,FL$p)
  names(FL$post_b) <- FL$colNames
  FL$post_b2 <- rep(0,FL$p)   # this is updated in the loop to get SD(post_b)
  names(FL$post_b2) <- FL$colNames
    
  # Initialize the coefficients to LS estimates
  FL$b <- CtX %*% e  # e was set to yStar initially
    
  e <- e - FL$X %*% FL$b # (e was set to yStar initially so this is actually yStar - FL$X %*% FL$b)
  
  if(saveChain){
    fname=paste(saveAt,"FIXED_b.dat",sep="")
    
    FL$NamefileOut <- fname
    unlink(fname) 
    FL$fileOut <- file(description=fname,open="w")
    
    write(FL$colNames, ncolumns  =  FL$p, file = FL$fileOut, append = TRUE)
  }

  
  FL$X <- as.vector(FL$X)  #store it as vector to save memory (stacked by columns)
  
  return(list(FL=FL,e=e))
}


### Marker component ########################################################################################################################################                 
#Function for initializing regression coefficients for marker effects. - Normal point & slab prior

set.list.Marker <- function(X, dfb, Sb, probLarge, saveAt, par.app, int.cond=NULL, gammaIndINT=NULL, saveChain)
{
  
  ML <- list(X=as.matrix(X), p=ncol(X), colNames=colnames(X), probLarge=probLarge)
  
  if (par.app == "common"){
    ML$txName <- "common"
  } else if (par.app == "group"){
    ML$txName <- ifelse(grepl("tx0", dimnames(X)[[2]], fixed=TRUE)[1], "tx0", "tx1")
  }else{
    ML$txName <- ifelse(grepl("tx1", dimnames(X)[[2]], fixed=TRUE)[1], "int", "main")
  }
  
  #Sum of the squares for each marker
  ML$x2 <- apply(ML$X, 2, function(x) sum(x^2))  
  
  if(!is.null(gammaIndINT)){
    tmp <- (gammaIndINT == 1)
    ML$gammaInd <- gammaIndINT
    ML$gammaInd[!tmp] <- rbinom(n = sum(!tmp), size = 1, prob = ML$probLarge)
    rm(tmp)
  }else{
    ML$gammaInd <- rbinom(n = ML$p, size = 1, prob = ML$probLarge)
  }
  
  
  # Number of markers with large effects (|gammaInd|):
  ML$countsLarge <- sum(ML$gammaInd)
  
  # Initial effect size variance
    ML$varB  <- dfb*Sb/(dfb + 2)
  # ML$varB  <- Sb/(dfb + 2)
  
  
  # Initial effects of markers with large effects (from prior):
  ML$b <- rep(0, ML$p)
  ML$b[ML$gammaInd==1] <- rnorm(ML$countsLarge, sd=sqrt(ML$varB))
  
  
  
  # Output files
  if(saveChain){
    fname <- paste(saveAt, ML$txName,"_MARKER_pars.dat",sep="")
    unlink(fname) 
    ML$fileOut <- file(description=fname,open="w")
    
    write(c("probLarge","varB"), ncolumns  =  2, file = ML$fileOut, append = TRUE)
  }

  #Objects for storing MCMC information 
  ML$post_varB <- 0
  ML$post_varB2 <- 0 # this is updated in the loop to get SD(post_varB)
  ML$post_gammaInd <- 0
  ML$post_probLarge <- 0
  ML$post_probLarge2 <- 0
  ML$post_b <- rep(0,ML$p)
  ML$post_b2 <- rep(0,ML$p) # this is updated in the loop to get SD(post_b)
  
  ML$X <- as.vector(ML$X)
  
  return(ML)  
}



### Marker component ########################################################################################################################################                 
#Function for initializing regression coefficients for marker effects. - Bayesian LASSO prior

set.list.Marker.BL <- function(X, lambda, saveAt, par.app, saveChain)
{
  
  ML <- list(X=as.matrix(X), p=ncol(X), colNames=colnames(X), lambda=lambda)
  
  ML$minAbsBeta <- 1e-8
  ML$lambda2=ML$lambda^2
  
  if (par.app == "common"){
    ML$txName <- "common"
  } else if (par.app == "group"){
    ML$txName <- ifelse(grepl("tx0", dimnames(X)[[2]], fixed=TRUE)[1], "tx0", "tx1")
  }else{
    ML$txName <- ifelse(grepl("tx1", dimnames(X)[[2]], fixed=TRUE)[1], "int", "main")
  }
  
  #Sum of the squares for each marker
  ML$x2 <- apply(ML$X, 2, function(x) sum(x^2))  
  
  # Initial effect size variances
  tau2 <- 1/ML$lambda2  # prior mean
  ML$tau2 <- rep(tau2, ML$p)
  ML$post_tau2 <- 0  
  
  # Initial effects:
  ML$b <- rep(0, ML$p)
  
  # Output files
  if(saveChain){
    fname <- paste(saveAt, ML$txName,"_lambda.dat",sep="")
    unlink(fname) 
    ML$fileOut <- file(description=fname,open="w")
  }
  
  #Objects for storing MCMC information 
  ML$post_lambda <- 0
  ML$post_lambda2 <- 0
  ML$post_b <- rep(0,ML$p)
  ML$post_b2 <- rep(0,ML$p) # this is updated in the loop to get SD(post_b)
  
  ML$X <- as.vector(ML$X)
  
  return(ML)  
}


### Random genetic component - NAIVE APPROACH ########################################################################################################################################                 
#Function for initializing the random genetic component.

set.list.G.naive <- function(Kg, txName, dfg, Sg, saveAt, tolD, saveChain)
{
  Kg <- as.matrix(Kg)
  GL <- list(nVarGcolumns = 1, txName = txName)
  
  tmp <- eigen(Kg)
  
  GL$V <- tmp$vectors
  GL$d <- tmp$values
    
  # Only those eigenvectors whose eigenvalues> tolD are kept. (in case K is not full rank)
  #Removing elements whose eigenvalues < tolD
  d.keep <- (GL$d > tolD)
  GL$numNonZeroD <- sum(d.keep)
  GL$d <- GL$d[d.keep]
  GL$V <- GL$V[, d.keep] # dim(V)=c(n,numNonZeroD)
  
  GL$Vt <- t(GL$V)
  
  GL$g <- rep(0, nrow(GL$V)) # gen. random effects
  GL$u <- rep(0, GL$numNonZeroD) # g=V %*% u where dim(V)=c(n,numNonZeroD)

  # Initial sigma_g^2
  GL$varG <- dfg*Sg/(dfg + 2)
  
  #Output files
  if(saveChain){
    fname <- paste(saveAt,"varG", GL$txName,".dat",sep="")
    unlink(fname) 
    GL$fileOut <- file(description=fname,open="w")
  }

  #Objects for storing information for MCMC iterations
  GL$post_varG <- 0
  GL$post_varG2 <- 0 # this is updated in the loop to get SD(post_varG)
  GL$post_g <- rep(0, nrow(GL$V))
  GL$post_g2 <- rep(0,nrow(GL$V))
  
  return(GL) 
}

### Random genetic component - IMPROVED APPROACH ########################################################################################################################################                 
#Function for initializing the random genetic component.

set.list.G <- function(K, tx, n, n_g, dfg, Sg, saveAt, tolD, saveChain)
{
  
  # Everything within g calculations will be done on objects ordered by tx
  # Outputted g will be reordered to match the original ordering
  ord.tx <- order(tx)
  ord.ord.tx <- order(ord.tx)
  
  K <- as.matrix(K)[ord.tx,ord.tx]
  GL <- list(ord.tx=ord.tx, ord.ord.tx=ord.ord.tx, nVarGcolumns=3)
  
  n1 <- n_g[1]
  n2 <- n_g[2]

  tmp <- eigen(K)
  
  V <- tmp$vectors
  GL$d <- tmp$values
  rm(tmp, K)
  
  # Only those eigenvectors whose eigenvalues> tolD are kept. (in case K is not full rank)
  #Removing elements whose eigenvalues < tolD
  d.keep <- (GL$d > tolD)
  GL$numNonZeroD <- sum(d.keep)
  GL$d <- GL$d[d.keep]
  V <- V[, d.keep] # dim(V)=c(n,numNonZeroD) 
    
  GL$V1 <- V[1:n1,]
  GL$V2 <- V[(n1+1):n,]
  
  
  # Set up things needed to create the sparse matrix within the updates

  A <- crossprod(GL$V1)
  B <- crossprod(GL$V2)
  
  zm1 <- matrix(0,nrow=GL$numNonZeroD,ncol=GL$numNonZeroD)
  zm2 <- matrix(0,nrow=GL$numNonZeroD,ncol=GL$numNonZeroD)
  IVC2 <- rbind(cbind(A,zm2),cbind(zm1,B)) 
  GL$Bbig <- as.vector(IVC2)  # this will be updated every iteration
  
  
  GL$block.diag.inds <- list(rep((1:GL$numNonZeroD),GL$numNonZeroD) + sort(rep((0:(GL$numNonZeroD-1))*2*GL$numNonZeroD,GL$numNonZeroD)),
                             rep((GL$numNonZeroD+2*GL$numNonZeroD^2)+(1:GL$numNonZeroD),GL$numNonZeroD) + sort(rep((0:(GL$numNonZeroD-1))*2*GL$numNonZeroD,GL$numNonZeroD)))
  
  GL$diag.inds <- ((1:(2*GL$numNonZeroD))*(2*GL$numNonZeroD)) - (((2*GL$numNonZeroD):1)-1)
  
  GL$offdiag.inds <- c(2*((1:GL$numNonZeroD)*GL$numNonZeroD) - ((GL$numNonZeroD:1)-1), 
                      2*GL$numNonZeroD*GL$numNonZeroD + ((1:(GL$numNonZeroD))*(2*GL$numNonZeroD)) - (((2*GL$numNonZeroD):(GL$numNonZeroD+1))-1))
  
  rm(V, A, B, IVC2, zm1, zm2)
  
  # Initial random effects
  # GL$gStar <- rep(0, 2*n) # gen. random effects with counter-factuals
  GL$u <- rep(0, 2*GL$numNonZeroD) # gStar= (I_2 %x% V) %*% u where dim(V)=c(n,numNonZeroD)
  GL$g <- rep(0, n) # gen. random effects, g= (C %*% gStar)[order(ord.tx)]
  
  # Initial Sigma_g^2
  SigmaG <- Sg/(dfg + 2 + 1) # prior mode of inverse wishart density
  GL$varG <- diag(SigmaG)
  GL$varG <- c(GL$varG, cov2cor(SigmaG)[1,2])
  names(GL$varG) <- c("varG1","varG2","rhoG")
  
  #Output files
  if(saveChain){
    fname <- paste(saveAt,"varG.dat",sep="")
    unlink(fname) 
    GL$fileOut <- file(description=fname,open="w")
    write(c("varG1","varG2","rhoG"), ncolumns  =  3, file = GL$fileOut, append = TRUE)
  }

  #Objects for storing information for MCMC iterations
  GL$post_varG1 <- 0
  GL$post_varG2 <- 0 
  GL$post_rhoG <- 0 
  GL$post_varG12 <- 0 # this is updated in the loop to get SD(post_varG1)
  GL$post_varG22 <- 0 # this is updated in the loop to get SD(post_varG2)
  GL$post_rhoG2 <- 0 # this is updated in the loop to get SD(post_rhoG)
  GL$post_g <- rep(0, n)
  GL$post_g2 <- rep(0, n)
  
  return(GL) 
}


### Sample Betas - FIXED ########################################################################################################################################                 
# Function to sample betas for the fixed part of the model.

sample_beta_FIXED <- function(b, X, n, p, e, varE) { 
  
  # Note: varE is a vector here and the cov. matrix of e is diag(varE).
  
  for(j in 1:p){
    
    xj <- X[(n*(j-1)+1):(j*n)]
    e <- e + xj * b[j]
    
    xjd <- xj/sqrt(varE)
    x2 <- sum(xjd * xjd)
    xe <- sum(e * xj / varE) # = t(e) %*% xj
    
    mu <- xe/x2
    sigma <- 1/sqrt(x2)
    
    b[j] <- rnorm(1, mean=mu, sd=sigma)
    
    e <- e - xj * b[j]
    
  }
  
  return(list(b=b,e=e))
}
 

### Sample Genetic Random effects ########################################################################################################################################                 
# Function to sample genetic random effects.

sample_G_naive <- function(g, varG, V, Vt, d, numNonZeroD, e, varE) { 

  e <- e + g # Update the residuals e =  y - W*alpha - X*beta (-g + g)

  # tcrossprod(Vt,t(e)) is faster than crossprod(e, V) and crossprod(V, e)
  
  Ve <- tcrossprod(Vt,t(e))/varE 
  varg.tmp <- varG * d
  C <- as.numeric(1/varg.tmp + 1/varE)
  sdU <- 1/sqrt(C)
  muU <- Ve/C
  u <- rnorm(n = numNonZeroD, mean = muU, sd = sdU) # length(u)=number of #nonzero diag elements in original D
  
  g <- as.vector(tcrossprod(V,t(u))) # length(g)= n (same as V %*% u but faster)

  
  # Update residuals
  e <- e - g
  
  return(list(g=g,u=u,e=e))
}


sample_G_C <- function(g, u, varG, varGinv, V1, V2, d, numNonZeroD, ord.tx, ord.ord.tx, Bbig, block.diag.inds, diag.inds, offdiag.inds, e, varE, n_g, n) { 
  
  # varE vector of size 2

  u.old <- u  
  
  n1 <- n_g[1]
  
  e <- e + g # Update the residuals e =  y - W*alpha - X*beta (-g + g)
  e <- e[ord.tx]
  
  IVe <- c(as.numeric(crossprod(V1,e[1:n1]/varE[1])),
           as.numeric(crossprod(V2,e[(n1+1):n]/varE[2]))) # this avoids tcrossprod(IVC,t(e))
  
  # varGinv <- get_Sigma_inv.vector(varG)
  
  Bbig[block.diag.inds[[1]]] <-  (1/varE[1])*Bbig[block.diag.inds[[1]]]
  Bbig[block.diag.inds[[2]]] <-  (1/varE[2])*Bbig[block.diag.inds[[2]]]
  
  Bbig[offdiag.inds] <- Bbig[offdiag.inds] + varGinv[3]*(1/d)
  Bbig[diag.inds] <-  Bbig[diag.inds] + c(varGinv[1]*(1/d),varGinv[2]*(1/d))
  Bdiag <- Bbig[diag.inds]
  
  u.samp <- .Call("sample_G", u.old, IVe, Bdiag, Bbig, 2*numNonZeroD, 2*numNonZeroD)
  u <- u.samp[[1]]
  
  u1 <- u[1:numNonZeroD]
  u2 <- u[(numNonZeroD+1):(2*numNonZeroD)]
  
  g <- c(as.numeric(tcrossprod(V1,t(u1))), 
         as.numeric(tcrossprod(V2,t(u2))))
  
  # faster than but same as:
  # gStar <- as.numeric(tcrossprod(IV,t(u))) 
  # g <- tcrossprod(C,t(gStar))
  # length(gStar)= 2n ordered by tx
  
  g <- g[ord.ord.tx] # length(g)= n
  
  # Update residuals (original order)
  e <- e[ord.ord.tx] - g
  
  return(list(g=g,u=u,e=e,gs=gs))
}


################################


#########################

get.inv.v <- function(varG) {
  # varG[1] = sigma_1^2, varG[2] = sigma_2^2, varG[3] = rho
  c <- 1/(varG[1]*varG[2]*(1-varG[3]^2))
  tmp <- c(varG[2],varG[1],-(sqrt(varG[1]*varG[2])*varG[3]))
  names(tmp) <- NULL
  return(c*tmp)
}

get_Sigma_inv.vector <- cmpfun(get.inv.v)

get.inv.m <- function(varG) {
  # varG[1] = sigma_1^2, varG[2] = sigma_2^2, varG[3] = rho
  c <- 1/(varG[1]*varG[2]*(1-varG[3]^2))
  tmp <- c(varG[2],varG[1],-(sqrt(varG[1]*varG[2])*varG[3]))
  c.tmp <- c*tmp
  varGinv.mat <- matrix(c(c.tmp[1],c.tmp[3],c.tmp[3],c.tmp[2]),nrow=2)
  return(varGinv.mat)
}

get_Sigma_inv.matrix <- cmpfun(get.inv.m)


get.Bbig_O <- function(varGinv, d, varE, Bdiag, IVC2offdiag_vec, i, j){
  
  part1 <- Bdiag
  part2 <- IVC2offdiag_vec
  part3 <- varE*varGinv[3]*(1/d)
  
  x <- c(part1,part2,part3)
  
  Bbig <- sparseMatrix(i = i, j = j, x = x, symmetric = TRUE, check=FALSE)
  return(Bbig)
}

get_Bbig <- cmpfun(get.Bbig_O)




##################################################################################################

#### Function to find beta parameters given initial probLarge and expected interval for probLarge



find.beta.pars <- function(probLarge0, r.int){
  r1 <- r.int[1]
  r2 <- r.int[2]
  
  probLarge0 <- probLarge0
  
  beta.par.optim.f <- function(a){
    b <- a*(1-probLarge0)/probLarge0
    abs((pbeta(r2, a, b) - pbeta(r1, a, b)) - 0.95)  # 95% Credible interval between r1 and r2
  }
  
  
  betaPARS <- numeric(2)
  opt.par  <- optimize(beta.par.optim.f, interval=c(0.01,15))
  betaPARS[1] <- opt.par$minimum
  betaPARS[2] <- betaPARS[1]*(1-probLarge0)/probLarge0
  
  return(betaPARS)
}



##################################################################################################


#### Function to get sets for cross validation

get.sets.CV <- function(n, fold, seeed){
  left <- fold - n %% fold
  set.seed(seeed) 
  if(left==fold){
    sets <- (rep(1:fold,ceiling(n/fold)))[order(runif(ceiling(n/fold)*fold))]
  } else{
    sets <- (rep(1:fold,ceiling(n/fold))[-(sample(1:fold,left,replace=F))])[order(runif(ceiling(n/fold)*fold)[-(1:left)])]
  }
return(sets)
}

##################################################################################################


#### Function to get training - test sets 

get.sets.test <- function(n, test.p, seeed){
  set.seed(seeed) 
  sets <- rep(c("train","test"),times=c(n-ceiling(n*test.p),ceiling(n*test.p)))[order(runif(n))]
  return(sets)
}

##################################################################################################
# Function to print out timing info
get.time.BGWP <- function(fm)
{
  cat("---------------------------------------------------------------------\n")
  cat(c(paste(c("  Total iteration time elapsed = ", "s\n  Average Time/Iter= "), 
              c(round(fm$Performance$totalIterTime, 2),  round(fm$Performance$AveIterTime, 3)), sep = "")), 
      "s\n")  
  cat("---------------------------------------------------------------------\n")
}
##################################################################################################

##################################################################################################
#rtrun draws from a truncated univariate normal distribution using the inverse CDF algorithm
#Arguments:
#mu: mean
#sigma: standard deviation
#a: lower bound
#b: upper bound
#NOTES: 1) This routine was taken from bayesm package, December 18, 2012
#       2) The inputs are not checked, 
#It is assumed that are ok.
rtrun <- function (mu, sigma, a, b) 
{
  FA = pnorm(((a - mu)/sigma))
  FB = pnorm(((b - mu)/sigma))
  return(mu + sigma * qnorm(runif(length(mu)) * (FB - FA) + FA))
}


#Extract the values of z such that y[i]=j
#z,y vectors, j integer
extract <- function(z,y,j) subset(as.data.frame(z,y),subset=(y==j))


#log-likelihood for ordinal data
#y: response vector
#predicted response vector, yHat=X%*%beta
#threshold
loglik_ordinal=function(y,yHat,threshold)
{
  sum=0
  n=length(y)
  for(i in 1:n)
  {
    sum=sum + log(pnorm(threshold[y[i] + 1]-yHat[i])-pnorm(threshold[y[i]]-yHat[i]))
  }
  return(sum)
}


##################################################################################################
# Function to calculate relative predictive gain (RPG)
get.RPG <- function(betahat, beta, X, y)
{
  n <- length(y)
  sj <- apply(X, 2, function(x) sum(x^2))/n
  
  mspe.hat <- sum(sj*((betahat-beta)^2))
  mspe.0 <- sum((mean(y)-y)^2)
  
  return((mspe.0-mspe.hat)/mspe.0)
}
##################################################################################################

### Prediction of future individuals (no variance yet) ########################################################################################################################################                 

# Note: Non-genetic variables not incorporated in this function yet.
# By default it assumes data.fixed=NULL in the original function call for BWGP.RCT.

# Note: Only implemented for par.app=="group" for now.

get.yHat.future <- function(X.f, K, out, tst, txgroup=NULL) 
{
  # X.f: Genotype data for the future individual
  # K: GRM of all subjects (both training and future)
  # out: BWGP.RCT output from training set
  # tst: vector of indices of future individuals
  
  if(!(out$par.app == "group")) stop("Not implemented.\n")
  
  if(is.null(txgroup)) txgroup <- out$txgroup-1
  
  n.f <- length(tst)
  
  tx.f <- rep(c(0,1),n.f)
  rows <- rep(1:nrow(X.f) , 2)
  X.f <- X.f[sort(rows),]
  
  w.f <- diag(2)
  rows <- rep(1:2 , n.f)
  w.f <- w.f[rows,]
  
  mu.f <- w.f %*% out$FIXED$b + ((tx.f==0)*1)*X.f %*% out$MARKER[[1]]$b + ((tx.f==1)*1)*X.f%*% out$MARKER[[2]]$b
  
  if (out$cov.app == "none"){ # i.e g=0
    yHatf <- mu.f
  } else {
    gHat <- out$G$g
    K00 <- K[-tst,-tst]
    Kf0 <- K[tst,-tst]
    if(n.f==1) Kf0 <- t(as.matrix(Kf0))
    
    if(out$cov.app == "adaptive"){
      Sigma00 <- rbind(cbind(out$G$varG1 * K00[txgroup==0,txgroup==0], sqrt(out$G$varG1*out$G$varG2)*out$G$rhoG* K00[txgroup==0,txgroup==1]),
                       cbind(sqrt(out$G$varG1*out$G$varG2)*out$G$rhoG* K00[txgroup==1,txgroup==0], out$G$varG2 * K00[txgroup==1,txgroup==1]))
      mug <- solve(Sigma00, gHat)
      
      if(n.f==1){
        Sigmaf0 <- rbind(c(out$G$varG1 * Kf0[,txgroup==0], sqrt(out$G$varG1*out$G$varG2)*out$G$rhoG* Kf0[,txgroup==1]),
                         c(sqrt(out$G$varG1*out$G$varG2)*out$G$rhoG* Kf0[,txgroup==0], out$G$varG2 * Kf0[,txgroup==1]))
      } else {
        Sigmaf0 <- rbind(cbind(out$G$varG1 * Kf0[,txgroup==0], sqrt(out$G$varG1*out$G$varG2)*out$G$rhoG* Kf0[,txgroup==1]),
                         cbind(sqrt(out$G$varG1*out$G$varG2)*out$G$rhoG* Kf0[,txgroup==0], out$G$varG2 * Kf0[,txgroup==1]))
      }
      yHatf <- mu.f + tcrossprod(Sigmaf0, t(mug))
    } else if(out$cov.app == "stratified"){
      mug1 <- solve(out$G$G1$varG * K00[txgroup==0,txgroup==0], gHat[txgroup==0])
      mug2 <- solve(out$G$G2$varG * K00[txgroup==1,txgroup==1], gHat[txgroup==1])
      
      yHatf.tmp <- c(out$G$G1$varG *tcrossprod(Kf0[,txgroup==0], t(mug1)),
                     out$G$G2$varG *tcrossprod(Kf0[,txgroup==1], t(mug2))) 
      yHatf <- mu.f + yHatf.tmp
    } 
  }
  
  ac <- matrix(c(-1,1),nrow=1)
  ac <- diag(n.f) %x% ac
  
  diffHat <- ac %*% yHatf
  return(list(diffHat=diffHat, yHatf=yHatf))  
}

##################################################################################################

