

###### ----------- Bayesian Whole Genome Prediction for RCTs ------------- ######

# Author: Bahar Erar

#################################################################################


BWGP.RCT <- function (y, tx, data.fixed=NULL, X.marker=NULL, K=NULL, J=NULL,
                      response.type="gaussian", a = NULL, b = NULL,
                      par.app=c("group","interaction","common","none"), 
                      fixed.prior=NULL, int.cond=NULL,
                      cov.app=c("stratified","adaptive","pooled","none"),
                      varyingE = TRUE,
                      nIter = 1500, burnIn = 500, thin = 5, saveAt = "", 
                      HPL, 
                      probLarge = NULL, probLargeTYPE="uniform", probLarge.betaPARS=NULL, probLargeLIMS=NULL,
                      lambda = NULL, lambdaTYPE="gamma",
                      tolDg = 1e-10, verbose = TRUE, saveChain=TRUE) 
{
  
  if (verbose) {
    welcome()
  }
  
  
  require(MCMCpack)
  if(cov.app == "adaptive"){
    require(Matrix)
  }
  dyn.load("get_samples.so")
  
  
  y0 <- y
  n <- length(y)
  
  if(is.null(names(y))){
    IDs <- 1:n
  }else{
    IDs <- names(y)    
  }
  
  
  if (!(response.type %in% c("gaussian", "binary", "ordinal"))){
    stop(" Only gaussian and ordinal responses are allowed\n")
  } 
  if (response.type == "binary") response.type <- "ordinal"
    
  if (response.type == "ordinal") {
    y <- factor(y, ordered = TRUE)
    lev <- levels(y)
    nclass <- length(lev)
    if (nclass == n) stop("The number of classes in y must be smaller than the number of observations.\n")
    if (nclass <= 1) stop("Data vector y should have at least 2 different values.\n")
    y <- as.integer(y)
    z <- y
  }
  
  Censored <- FALSE
  if (response.type == "gaussian") {
    if ((!is.null(a)) | (!is.null(b))) {
      Censored <- TRUE
      a <- as.vector(a)
      b <- as.vector(b)
      if ((length(a) != n) | (length(b) != n)) stop(" y, a and b must have the same dimension.\n")
    }
  }

  
  # Check that the given approaches are implemented
  if (!(par.app %in% c("group","interaction","common","none"))) stop(" Given parameterization approach (par.app) is not implemented.")
  if (par.app == "interaction") {
    if(is.null(int.cond) | (!(int.cond %in% c("independent","dependent")))) stop("Specify int.cond as 'independent' or 'dependent' if par.app='interaction'")
  }
  if (par.app %in% c("group","common")) int.cond <- "independent"
  int.cond <- (int.cond=="dependent")
  
  if (!(cov.app %in% c("stratified","adaptive","pooled","none"))) stop(" Given covariance modeling approach (cov.app) is not implemented.\n")
  if ((cov.app == "pooled") & (varyingE == TRUE)) stop(" Varying error variances under the pooled approach is not implemented.\n")
  
  
  # Check that prior model for fixed genetic effects is implemented
  if (!(par.app=="none")){
    if (is.null(fixed.prior)) fixed.prior <- "NormalMixture"
    if (!(fixed.prior %in% c("NormalMixture","BL"))) {
      fixed.prior <- "NormalMixture"
      warning("Given prior is not implemented for fixed genetic effects. Set to Normal point & slab mixture as default.\n", call.=FALSE)
    }
    if(fixed.prior == "BL") {
      if(int.cond) {
        int.cond <- FALSE
        warning("Interaction effects cannot be dependent under the BL prior. Set to independent.\n", call.=FALSE)
      }
    }
  }

    
  # Check that all design matrix and tx vector are of matching size with y

  if (!(length(tx) == n))   stop(" Length of treatment vector does not match the response vector.\n")
  
  if (!(cov.app=="none")){
    if (is.null(K)){
      stop(" Provide a genetic relationship matrix (K).\n")
    } else {
      if (!(sum(dim(K)== c(n,n))==2))   stop(" Matrix size of K does not match the response vector.\n")
      if(any(is.na(K))){ stop("NAs in K not allowed.")}
    }
  } else {
    if (par.app=="none") stop("Provide at least one of par.app or cov.app as not 'none'.\n")
  }
  
  if (!(par.app=="none")){
    if (is.null(X.marker)){
      stop(" Provide a genotype matrix.\n")
    } else {
      if (!(dim(X.marker)[1]  == n))   stop(" Matrix size (of genotypes) does not match the response vector.\n")
      if(any(is.na(X.marker))){ stop("NAs in the genotype matrix not allowed.")}
    }
  }

  # Check that tx only includes 0 and 1s (0: control, 1:tx)
  if (length(unique(tx)) < 2){
    stop(" Treatment vector should include two groups: 0: control, 1:tx.")
  } else if (!(length(unique(tx)) == 2)){
    stop(" Only 0 and 1s allowed in treatment vector: 0: control, 1:tx.")
  } else if (sum(!(sort(unique(tx)) == 0:1 ))>0L) {
    stop(" Only 0 and 1s allowed in treatment vector: 0: control, 1:tx.")
  }

  if (par.app %in% c("group","interaction")){
    if(fixed.prior == "NormalMixture") {
      if(length(probLarge)==1) probLarge <- rep(probLarge,2)
      if(length(HPL$Sb0)==1) HPL$Sb0 <- rep(HPL$Sb0,2)
      if(length(HPL$dfb0)==1) HPL$dfb0 <- rep(HPL$dfb0,2)
    }
  }
  

  # Convert tx into factor:
  tx <- factor(tx, levels=c("0","1"))
  n_g <- as.numeric(table(tx))
  tx <- as.integer(tx)   # tx is now 1 and 2
  
  if (verbose) {
    cat("Initializing...\n")
  }
  
  # The path that will be used if saveAt wasn't given
  if (saveAt == "") {
    saveAt <- paste(getwd(), "/", sep = "")
  }
  
  whichNa <- which(is.na(y))
  nNa <- length(whichNa)
  
  whichNa_g <- list(which(is.na(y[tx==1])), which(is.na(y[tx==2])))
  nNa_g <- c(length(whichNa_g[[1]]),length(whichNa_g[[2]]))
  
  post_logLik <- 0
  
  if (response.type == "ordinal") {
    countsZ <- table(z)
    threshold <- qnorm(p = c(0, cumsum(as.vector(countsZ)/n)))
    y <- rtrun(mu = 0, sigma = 1, a = threshold[z], b = threshold[(z + 1)])
    yInitGroupMeans <- rep(0,2)
    post_threshold <- 0
    post_threshold2 <- 0
    post_prob <- matrix(nrow = n, ncol = nclass, 0)
    post_prob2 <- post_prob
  } else {
    yInitGroupMeans <- c(mean(y[tx==1], na.rm = TRUE), mean(y[tx==2], na.rm = TRUE))
  }
  
  yStar <- y 
  
  yHat <- numeric(n)
  yHat[tx==1] <- yInitGroupMeans[1]
  yHat[tx==2] <- yInitGroupMeans[2]
  
  if (nNa > 0) {
    yStar[whichNa] <- yHat[whichNa] 
  }
  
  e <- yStar # this is updated in set.list.Fixed() as   e <- e - FL$X %*% FL$b
  
  post_yHat <- rep(0, n)
  post_yHat2 <- rep(0, n)
  
  
  # Set fixed model 
  tx.n <- tx
  tx <- factor(tx)
  
  if(is.null(data.fixed)){ # data.fixed is taken as intercept + tx or tx1 + tx2 based on parameterization
    data.fixed.model <- data.frame(tx)
    if (par.app == "group"){
      mf.fixed <- model.frame(formula="~ -1 + tx",data=data.fixed.model)  
    } else {
      mf.fixed <- model.frame(formula="~ tx",data=data.fixed.model)
    }
    rm(data.fixed.model) 
  } else if(is.list(data.fixed)){  # data.fixed is a list containing formula and data
    if(!(length(data.fixed)==2)) stop("If data.fixed is given as a list, it has to have 2 components: formula and data.\n")
    if(is.null(names(data.fixed))) names(data.fixed) <- c("formula","data")
    fixed.form <- data.fixed$formula
    fixed.data <- data.fixed$data
    if (!(dim(fixed.data)[1] == n) )  stop(" Number of rows in data.fixed does not match the number of phenotypes in y.\n")
    if(any(is.na(fixed.data))){ stop("NAs in design matrix of fixed effects not allowed.")}
    if(is.matrix(fixed.data)) fixed.data <- data.frame(fixed.data)
    mf.fixed <- model.frame(formula=fixed.form, data=fixed.data)  
    rm(fixed.form,fixed.data)
  } else if(is.matrix(data.fixed) | is.data.frame(data.fixed)){
    if (!(dim(data.fixed)[1] == n) )  stop(" Number of rows in data.fixed does not match the number of phenotypes in y.\n")
    if(any(is.na(data.fixed))){ stop("NAs in design matrix of fixed effects not allowed.")}
    if(is.matrix(data.fixed)) data.fixed <- data.frame(data.fixed)
    if (par.app == "group"){
      mf.fixed <- model.frame(formula=paste("~ -1 + tx + ",paste(names(data.fixed), collapse=" + "),sep=""), data=data.fixed)  
    } else {
      mf.fixed <- model.frame(formula=paste("~ tx + ",paste(names(data.fixed), collapse=" + "),sep=""), data=data.fixed)
    }
  } else {
    stop("data.fixed has to be either a data.frame, a matrix, a list containing a formula and data or NULL.\n")
  }
  
  # Initialize fixed effects
  X.fixed <- model.matrix(attr(mf.fixed, "terms"), data=mf.fixed)
  tmp <- set.list.Fixed(X.fixed = X.fixed, saveAt = saveAt, e = e, saveChain=saveChain) 
  FIXED.list <- tmp$FL
  e <- tmp$e # yStar - FL$X %*% FL$b
  
  tx <- tx.n
  rm(X.fixed, tmp, mf.fixed, tx.n)
  
  rm(data.fixed)

  
  if(response.type=="gaussian"){
    # Initial residual variance = mode of the prior scaled-inv chisq:
    #   varE = HPL$df0*HPL$S0/(HPL$df0 + 2)
    if(varyingE){
      varE <- c(as.numeric(var(e[tx==1], na.rm = TRUE)* HPL$S0/(HPL$df0 + 2)),
                as.numeric(var(e[tx==2], na.rm = TRUE)* HPL$S0/(HPL$df0 + 2)))
    } else {
      varE <- rep(as.numeric(var(e, na.rm = TRUE)* HPL$S0/(HPL$df0 + 2)), 2)
    }
  } else {
    varE <- rep(1,2)
  }
  
  sdE <- sqrt(varE)
  
  if(varyingE){
    post_varE <- rep(0,2)
    post_varE2 <- rep(0,2)
    
    if(saveChain){
      fname <- paste(saveAt, "varE.dat", sep = "")
      # Remove existing output files from previous runs
      unlink(fname) 
      fileOutVarE <- file(description = fname, open = "w")
      write(c("varE1","varE2"), ncolumns  =  2, file = fileOutVarE, append = TRUE)
    }

  } else {
    post_varE <- 0
    post_varE2 <- 0
    
    if(saveChain){
      fname <- paste(saveAt, "varE.dat", sep = "")
      # Remove existing output files from previous runs
      unlink(fname) 
      fileOutVarE <- file(description = fname, open = "w")
      write("varE", ncolumns  =  1, file = fileOutVarE, append = TRUE)
    }
  }
  
  
  # Set marker model 
  
  if (!(par.app == "none")){
    J <- dim(X.marker)[2]
    
    if (fixed.prior == "BL") {
      if(is.null(lambda)) {
        lambda <- 30
        warning("Initial value for probLarge not provided. Set to 30 arbitrarily.\n", call.=FALSE)
      } else {
        if(lambda < 0) {
          lambda <- 30
          warning("lambda cannot be negative. Set to 30 arbitrarily.\n", call.=FALSE)
        }
      }
      if(!(lambdaTYPE %in% c("fixed","gamma"))) {
        lambdaTYPE <- "gamma"
        warning("Only 'fixed' or 'gamma' options are implemented for lambdaTYPE. Set to 'gamma' by default.\n", call.=FALSE)
      }  
      if(lambdaTYPE == "gamma"){
        if(!("shape" %in% names(HPL))) {
          HPL$shape <- 1.1
          cat(paste("Shape parameter not given for BL prior. Set to ", HPL$shape, ".\n",sep=""))
        }
        if(!("rate" %in% names(HPL))) {
          HPL$rate <-(HPL$shape-1)/(lambda^2)
          cat(paste("Rate parameter not given for BL prior. Set to ", HPL$rate, ".\n",sep=""))
        }
      }
    } 
    

    if(is.null(dimnames(X.marker)[[2]])){
      dimnames(X.marker)[[2]] <- paste("M",1:J,sep="")
    }
    
    X.marker.save <- X.marker
    
    if (par.app == "common"){
      nmg <- 1
      # Initialize marker effects
      if(fixed.prior == "NormalMixture"){
        MARKER.list <- list(set.list.Marker(X = X.marker, HPL$dfb0, HPL$Sb0, probLarge, saveAt = saveAt, par.app = par.app, gammaIndINT=NULL, saveChain=saveChain))
      } else {
        MARKER.list <- list(set.list.Marker.BL(X = X.marker, lambda=lambda, saveAt = saveAt, par.app = par.app, saveChain=saveChain))
      }
      
    } else if(par.app == "group") {
      nmg <- 2
      X.marker.tmp <- matrix(NA,dim(X.marker)[1],2*J)
      X.marker.tmp[,1:J] <-((tx==1)*1)* X.marker
      X.marker.tmp[,(dim(X.marker)[2]+1):(2*J)] <- ((tx==2)*1)*X.marker
      dimnames(X.marker.tmp)[[1]] <- dimnames(X.marker)[[1]]
      dimnames(X.marker.tmp)[[2]] <- c(paste(dimnames(X.marker)[[2]], "tx0", sep=":"),
                                       paste(dimnames(X.marker)[[2]], "tx1", sep=":"))
      X.marker <- X.marker.tmp
      
      termInd <- grepl("tx1", dimnames(X.marker)[[2]], fixed=TRUE)
      M.names0 <- dimnames(X.marker)[[2]][!termInd]
      M.names1 <- dimnames(X.marker)[[2]][termInd]
      # Initialize marker effects
      if(fixed.prior == "NormalMixture"){
        MARKER.list <- list(set.list.Marker(X = X.marker[,M.names0], HPL$dfb0[1], HPL$Sb0[1], probLarge[1], saveAt = saveAt, par.app = par.app, gammaIndINT=NULL, saveChain=saveChain),
                            set.list.Marker(X = X.marker[,M.names1], HPL$dfb0[2], HPL$Sb0[2], probLarge[2], saveAt = saveAt, par.app = par.app, gammaIndINT=NULL, saveChain=saveChain))
      } else { # BL
        MARKER.list <- list(set.list.Marker.BL(X = X.marker[,M.names0], lambda=lambda, saveAt = saveAt, par.app = par.app, saveChain=saveChain),
                            set.list.Marker.BL(X = X.marker[,M.names1], lambda=lambda, saveAt = saveAt, par.app = par.app, saveChain=saveChain))
      }
      rm(X.marker.tmp)
      
    } else {
      nmg <- 2
      X.marker.tmp <- matrix(NA,dim(X.marker)[1],2*J)
      X.marker.tmp[,1:J] <- X.marker
      X.marker.tmp[,(J+1):(2*J)] <- ((tx==2)*1)*X.marker
      dimnames(X.marker.tmp)[[1]] <- dimnames(X.marker)[[1]]
      dimnames(X.marker.tmp)[[2]] <- c(dimnames(X.marker)[[2]],paste(dimnames(X.marker)[[2]], "tx1", sep=":"))
      X.marker <- X.marker.tmp
      
      termInd <- grepl("tx1", dimnames(X.marker)[[2]], fixed=TRUE)
      M.namesM <- dimnames(X.marker)[[2]][!termInd]
      M.namesI <- dimnames(X.marker)[[2]][termInd]
      
      # Initialize marker effects
      if(fixed.prior == "NormalMixture"){
        tmp1 <- set.list.Marker(X = X.marker[,M.namesI], dfb=HPL$dfb0[2], Sb=HPL$Sb0[2], probLarge[2], saveAt = saveAt, par.app = par.app, gammaIndINT=NULL, saveChain=saveChain)
        if (int.cond){
          tmp2 <- set.list.Marker(X = X.marker[,M.namesM], dfb=HPL$dfb0[1], Sb=HPL$Sb0[1], probLarge[1], saveAt = saveAt, par.app = par.app, gammaIndINT = tmp$gammaInd, saveChain=saveChain)
        } else {
          tmp2 <- set.list.Marker(X = X.marker[,M.namesM], dfb=HPL$dfb0[1], Sb=HPL$Sb0[1], probLarge[1], saveAt = saveAt, par.app = par.app, gammaIndINT = NULL, saveChain=saveChain)
        }
      } else { #BL
        tmp1 <- set.list.Marker.BL(X = X.marker[,M.namesI], lambda=lambda, saveAt = saveAt, par.app = par.app, saveChain=saveChain)
        tmp2 <- set.list.Marker.BL(X = X.marker[,M.namesM], lambda=lambda, saveAt = saveAt, par.app = par.app, saveChain=saveChain)
      }
      
      MARKER.list <- list(tmp1, tmp2)
      # Note: MARKER.list[[1]]=interaction effects and MARKER.list[[2]]=main effects
      
      rm(X.marker.tmp,tmp1, tmp2)
    }
    rm(X.marker)
    
    if (fixed.prior == "NormalMixture") {
      # Set prior on pi if probLargeTYPE="beta"
      if(probLargeTYPE=="beta") {
        if(is.null(probLarge.betaPARS)){
          if(is.null(probLargeLIMS)) probLargeLIMS <- c(1/J, (J-1)/J)
          if (verbose) cat("Beta prior hyperparameters not provided. Set to:\n")
          for (k in 1:nmg) {
            bpars <- find.beta.pars(MARKER.list[[k]]$probLarge, probLargeLIMS)
            MARKER.list[[k]]$betaPARS <- round(bpars,2)
            if (verbose) cat(paste("     Group ", k, ": (", round(bpars,2)[1], ", ", round(bpars,2)[2], ")\n", sep=""))
          }
        } else {
          for (k in 1:nmg) MARKER.list[[k]]$betaPARS <- probLarge.betaPARS
        }
      } else if (probLargeTYPE=="uniform") {
        for (k in 1:nmg) MARKER.list[[k]]$betaPARS <- c(1,1)
      } else if (probLargeTYPE=="loguniform") {
        if(is.null(probLargeLIMS)) probLargeLIMS <- c(1/J, (J-1)/J)
        for (k in 1:nmg){
          if(MARKER.list[[k]]$probLarge < probLargeLIMS[1]){
            MARKER.list[[k]]$probLarge <- round(probLargeLIMS[1] + (probLargeLIMS[2]-probLargeLIMS[1])/2,10)
            cat(paste("Given probLarge for group ", k," outside possible limits. Set to: ", MARKER.list[[k]]$probLarge,".\n", sep=""))
          } else if(MARKER.list[[k]]$probLarge > probLargeLIMS[2]){
            MARKER.list[[k]]$probLarge <- round(probLargeLIMS[1] + (probLargeLIMS[2]-probLargeLIMS[1])/2,10)
            cat(paste("Given probLarge for group ", k," outside possible limits. Set to: ", MARKER.list[[k]]$probLarge,".\n", sep=""))
          }
        }
      }
    }
    
  }

  
  # Initialize random genetic component, g
  if (!(cov.app == "none")){
    if (cov.app == "stratified"){
      if(!(length(HPL$Sg0)==1)) stop("The scale parameter (Sg) of the genetic variance hyperprior distribution has to be a scalar under the stratified approach.")
      G.list <- vector(mode="list",length=2)
      for(k in 1:2){
        Kg <- K[tx==k,tx==k]
        G.list[[k]] <- set.list.G.naive(Kg = Kg, txName=k,
                                        dfg=HPL$dfg0, Sg=HPL$Sg0, saveAt = saveAt, tolD=tolDg, saveChain=saveChain)
      }
    } else if(cov.app == "pooled") { 
      if(!(length(HPL$Sg0)==1)) stop("The scale parameter (Sg) of the genetic variance hyperprior distribution has to be a scalar under the pooled approach.")
      G.list <- set.list.G.naive(Kg = K, txName="pooled",
                                 dfg=HPL$dfg0, Sg=HPL$Sg0, saveAt = saveAt, tolD=tolDg, saveChain=saveChain)
      
    } else {
      if(length(HPL$Sg0)==1) stop("The scale parameter (Sg) of the genetic variance hyperprior distribution has to be a 2x2 matrix under the improved improach.")
      G.list <- set.list.G(K = K, tx=tx, n=n, n_g=n_g, HPL$dfg0, HPL$Sg0, saveAt = saveAt, tolD=tolDg, saveChain=saveChain)
      varGinv <- get_Sigma_inv.vector(G.list$varG)
    }
    
  }
  
  rm(K)
    
  
  # THE MAIN LOOP:
  
  if (verbose) cat("Starting MCMC...\n")
  
  
  nSums <- 0 # to count the number of samples used to obtain the posterior means (in the end, will be equal to (nIter - burnIn)/thin )
  tot.time <- 0

  time <- proc.time()[3]
  
  for (i in 1:nIter) {
    
    # Update the fixed effects
    varErep <- rep(varE[1], n)
    if(varyingE) varErep[tx==2] <- varE[2] 
    samp.F <- sample_beta_FIXED(b=FIXED.list$b, X=FIXED.list$X,
                                n=n, p=FIXED.list$p, e=e, varE=varErep)  
    FIXED.list$b <- samp.F$b
    e <- samp.F$e
    
    
    # Update the genetic random effects
    if (!(cov.app == "none")){
      if (cov.app == "stratified"){
        for(k in 1:2){
          # Update g
          samp.G <- sample_G_naive(g=G.list[[k]]$g, varG=G.list[[k]]$varG, V=G.list[[k]]$V, Vt=G.list[[k]]$Vt,
                                   d=G.list[[k]]$d, numNonZeroD=G.list[[k]]$numNonZeroD, 
                                   e=e[tx==k], varE=varE[k])  
          G.list[[k]]$g <- samp.G$g
          G.list[[k]]$u <- samp.G$u
          e[tx==k] <- samp.G$e
          
          # Update sigma2_g
          Sg <- as.numeric(crossprod(G.list[[k]]$u/sqrt(G.list[[k]]$d))) + HPL$Sg0
          dfg <- G.list[[k]]$numNonZeroD + HPL$dfg0
          # first term = t(u) %*% solve(diag(d)) %*% u but faster
          # first term = g'Kg = g'VD^(-1)V'g = u'D^(-1)u since g=Vu
          G.list[[k]]$varG <- Sg/rchisq(n = 1, df = dfg)
          rm(samp.G)
        }

      } else if (cov.app == "pooled"){ 
        varEcommon <- varE[1]
        # Update g
        samp.G <- sample_G_naive(g=G.list$g, varG=G.list$varG, V=G.list$V, Vt=G.list$Vt,
                                 d=G.list$d, numNonZeroD=G.list$numNonZeroD, 
                                 e=e, varE=varEcommon)  
        G.list$g <- samp.G$g
        G.list$u <- samp.G$u
        e <- samp.G$e
        
        # Update sigma2_g
        Sg <- as.numeric(crossprod(G.list$u/sqrt(G.list$d))) + HPL$Sg0
        dfg <- G.list$numNonZeroD + HPL$dfg0
        # first term = t(u) %*% solve(diag(d)) %*% u but faster
        # first term = g'Kg = g'VD^(-1)V'g = u'D^(-1)u since g=Vu
        G.list$varG <- Sg/rchisq(n = 1, df = dfg)
        rm(samp.G)
        
      } else {
        # Update g
        samp.G <- sample_G_C(g=G.list$g, u=G.list$u, varG=G.list$varG, varGinv=varGinv, V1=G.list$V1, V2=G.list$V2, d=G.list$d, G.list$numNonZeroD, ord.tx=G.list$ord.tx, ord.ord.tx=G.list$ord.ord.tx, Bbig=G.list$Bbig, block.diag.inds=G.list$block.diag.inds, diag.inds=G.list$diag.inds, offdiag.inds=G.list$offdiag.inds, e=e, varE=varE, n_g=n_g, n=n)        
        G.list$g <- samp.G$g
        G.list$u <- samp.G$u
        e <- samp.G$e
        
        
        # Update Sigma_g
        Gmat.tmp <- matrix(G.list$u,ncol=2)/sqrt(G.list$d)   # D^(-1) %*% U where U is the nx2 matrix of u's (which are ordered by tx) 
        Sg <- HPL$Sg0 + crossprod(Gmat.tmp)
        dfg <- G.list$numNonZeroD + HPL$dfg0
        tmp <- riwish(dfg, Sg)
        G.list$varG[1:3] <-  c(diag(tmp), cov2cor(tmp)[1,2])
        rm(tmp,Gmat.tmp,samp.G)
        varGinv <- get_Sigma_inv.vector(G.list$varG)

      }
    }

    if (!(par.app == "none")){
      # Update large marker effects
      
      if (int.cond){
        
        # Large interaction effects (beta_I)
        samp.M <- .Call("sample_beta_MARKER_DEP", MARKER.list[[1]]$b, MARKER.list[[1]]$gammaInd, MARKER.list[[1]]$probLarge, MARKER.list[[1]]$varB, MARKER.list[[1]]$X, MARKER.list[[1]]$x2, n, MARKER.list[[1]]$p, e, varE, rep(0,n), 0)
        names(samp.M) <- c("gammaInd", "e", "b")
        
        MARKER.list[[1]]$b <- samp.M$b
        MARKER.list[[1]]$gammaInd <- samp.M$gammaInd
        e <- samp.M$e 
        MARKER.list[[1]]$countsLarge <- sum(MARKER.list[[1]]$gammaInd)
        
        # Large main effects
        samp.M <- .Call("sample_beta_MARKER_DEP", MARKER.list[[2]]$b, MARKER.list[[2]]$gammaInd, MARKER.list[[2]]$probLarge, MARKER.list[[2]]$varB, MARKER.list[[2]]$X, MARKER.list[[2]]$x2, n, MARKER.list[[2]]$p, e, varE, MARKER.list[[1]]$gammaInd, 1)
        names(samp.M) <- c("gammaInd", "e", "b")
        
        MARKER.list[[2]]$b <- samp.M$b
        MARKER.list[[2]]$gammaInd <- samp.M$gammaInd
        e <- samp.M$e 
        MARKER.list[[2]]$countsLarge <- sum(MARKER.list[[2]]$gammaInd)
        
        rm(samp.M)
        
        for(k in 1:nmg){
          # Update pi_I
          if (probLargeTYPE=="beta" | probLargeTYPE=="uniform"){
            MARKER.list[[k]]$probLarge <- rbeta(n = 1, 
                                                shape1 = MARKER.list[[k]]$betaPARS[1] + MARKER.list[[k]]$countsLarge,
                                                shape2 = MARKER.list[[k]]$betaPARS[2] + MARKER.list[[k]]$p - MARKER.list[[k]]$countsLarge)
          } else if (probLargeTYPE=="loguniform"){
            proposed.pL <- rbeta(n = 1, 
                                 shape1 = MARKER.list[[k]]$countsLarge,
                                 shape2 = 1 + MARKER.list[[k]]$p - MARKER.list[[k]]$countsLarge)
            accept <- (probLargeLIMS[1] <= proposed.pL) & (probLargeLIMS[2] >= proposed.pL)
            if(accept){
              MARKER.list[[k]]$probLarge <- proposed.pL
            }
          }
          
          # Update sigma2_beta 
          dfb <- HPL$dfb0[k] + MARKER.list[[k]]$p
          bb <- sum(MARKER.list[[k]]$b^2)
          
          Sb <- HPL$Sb0[k]*HPL$dfb0[k] + bb
          MARKER.list[[k]]$varB <- Sb/rchisq(n = 1, df = dfb) # generate scaled inv chisq: (df * scale)/rchisq(n, df = df) where scale is Sb/dfb
        }
        
      } else {
        if(fixed.prior == "NormalMixture"){
          #For each marker effect group
          for(k in 1:nmg){
            samp.M <- .Call("sample_beta_MARKER", MARKER.list[[k]]$b, MARKER.list[[k]]$gammaInd, MARKER.list[[k]]$probLarge, MARKER.list[[k]]$varB, MARKER.list[[k]]$X, MARKER.list[[k]]$x2, n, MARKER.list[[k]]$p, e, varE[k])
            names(samp.M) <- c("gammaInd", "e", "b")
            
            MARKER.list[[k]]$b <- samp.M$b
            MARKER.list[[k]]$gammaInd <- samp.M$gammaInd
            e <- samp.M$e 
            MARKER.list[[k]]$countsLarge <- sum(MARKER.list[[k]]$gammaInd)
            
            rm(samp.M)
            
            # Update pi
            if (probLargeTYPE=="beta" | probLargeTYPE=="uniform"){
              MARKER.list[[k]]$probLarge <- rbeta(n = 1, 
                                                  shape1 = MARKER.list[[k]]$betaPARS[1] + MARKER.list[[k]]$countsLarge,
                                                  shape2 = MARKER.list[[k]]$betaPARS[2] + MARKER.list[[k]]$p - MARKER.list[[k]]$countsLarge)
            } else if (probLargeTYPE=="loguniform"){
              proposed.pL <- rbeta(n = 1, 
                                   shape1 = MARKER.list[[k]]$countsLarge,
                                   shape2 = 1 + MARKER.list[[k]]$p - MARKER.list[[k]]$countsLarge)
              accept <- (probLargeLIMS[1] <= proposed.pL) & (probLargeLIMS[2] >= proposed.pL)
              if(accept) {
                MARKER.list[[k]]$probLarge <- proposed.pL
              }
            }
            
            # Update sigma2_beta 
            dfb <- HPL$dfb0[k] + MARKER.list[[k]]$p
            bb <- sum(MARKER.list[[k]]$b^2)
            
            Sb <- HPL$Sb0[k]*HPL$dfb0[k] + bb
            MARKER.list[[k]]$varB <- Sb/rchisq(n = 1, df = dfb) # generate scaled inv chisq: (df * scale)/rchisq(n, df = df) where scale is Sb/dfb
          }
        } else { # BL
          for(k in 1:nmg){
            varBj <- MARKER.list[[k]]$tau2 * varE[k]
            samp.M <- .Call("sample_beta_MARKER_BL", MARKER.list[[k]]$b, varBj, MARKER.list[[k]]$X, MARKER.list[[k]]$x2, n, MARKER.list[[k]]$p, e, varE[k], MARKER.list[[k]]$minAbsBeta)
            names(samp.M) <- c("e", "b")
            
            MARKER.list[[k]]$b <- samp.M$b
            e <- samp.M$e 
            
            nu <- sqrt(varE[k]) * MARKER.list[[k]]$lambda/abs(MARKER.list[[k]]$b)
            try(tmp <- rinvGauss(n = MARKER.list[[k]]$p, nu = nu, 
                                 lambda = MARKER.list[[k]]$lambda2))
            if (!is.null(tmp) && !any(tmp < 0)) {
              if (!any(is.na(sqrt(tmp)))) {
                MARKER.list[[k]]$tau2 = 1/tmp
              } else {
                warning(paste("tau2 was not updated in iteration", i, "due to numeric problems with beta\n", sep = " "), immediate. = TRUE)
              }
            } else {
              warning(paste("tau2 was not updated  in iteration", i, "due to numeric problems with beta\n", sep = " "), immediate. = TRUE)
            }
            if (lambdaTYPE == "gamma") {
              rate <- sum(MARKER.list[[k]]$tau2)/2 + HPL$rate
              shape <- MARKER.list[[k]]$p + HPL$shape
              MARKER.list[[k]]$lambda2 <- rgamma(rate = rate, shape = shape, n = 1)
              if (!is.na(MARKER.list[[k]]$lambda2)) {
                MARKER.list[[k]]$lambda <- sqrt(MARKER.list[[k]]$lambda2)
              } else {
                warning(paste("lambda was not updated in iteration", i, "due to numeric problems with beta\n", sep = " "), immediate. = TRUE)
              }
            }
            
          }
        }
        
      }
    }
    
    
    yHat <- yStar - e
    
    
    # Update sigma2_e
    
    if (response.type == "gaussian") {
      if(varyingE){
        for(k in 1:2){
          Se <- sum(e[tx==k] * e[tx==k]) + HPL$S0
          dfe <- n_g[k] + HPL$df0 
          varE[k] <- Se/rchisq(n = 1, df = dfe)
          sdE <- sqrt(varE)
        }
      } else {
        Se <- sum(e * e) + HPL$S0
        dfe <- n + HPL$df0 
        tmp.e <- Se/rchisq(n = 1, df = dfe)
        varE <- rep(tmp.e,2)
        sdE <- sqrt(tmp.e)
      }
      
      if(nNa>0){
        if (Censored) {
          if(varyingE){
            for(k in 1:2){
              if(nNa_g[k]>0){
                yStar[tx==k][whichNa_g[[k]]] <- rtrun(mu = yHat[tx==k][whichNa_g[[k]]], 
                                                      a = a[tx==k][whichNa_g[[k]]], b = b[tx==k][whichNa_g[[k]]], sigma = sdE[k])
                e[tx==k][whichNa_g[[k]]] <- yStar[tx==k][whichNa_g[[k]]] - yHat[tx==k][whichNa_g[[k]]]
              }
            }
          } else {
            yStar[whichNa] <- rtrun(mu = yHat[whichNa], a = a[whichNa], b = b[whichNa], sigma = sdE)
            e[whichNa] <- yStar[whichNa] - yHat[whichNa]
          }
        } else {
          if(varyingE){
            for(k in 1:2){
              if(nNa_g[k]>0){
                yStar[tx==k][whichNa_g[[k]]] <- yHat[tx==k][whichNa_g[[k]]] + rnorm(n = nNa_g[k], sd = sdE[k])
                e[tx==k][whichNa_g[[k]]] <- yStar[tx==k][whichNa_g[[k]]] - yHat[tx==k][whichNa_g[[k]]]
              }
            }
          } else {
            yStar[whichNa] <- yHat[whichNa] + rnorm(n = nNa, sd = sdE)
            e[whichNa] <- yStar[whichNa] - yHat[whichNa]
          }
        }
      }

      
    } else {
      
      varE <- 1
      sdE <- 1
      if (nNa == 0) {
        yStar <- rtrun(mu = yHat, sigma = 1, a = threshold[z], b = threshold[(z + 1)])
      } else {
        yStar[-whichNa] <- rtrun(mu = yHat[-whichNa], 
                                sigma = 1, a = threshold[z[-whichNa]], b = threshold[(z[-whichNa] + 1)])
        yStar[whichNa] <- yHat[whichNa] + rnorm(n = nNa, sd = sdE)
      }
      if (nNa == 0) {
        for (m in 2:nclass) {
          lo <- max(max(extract(yStar, z, m - 1)), threshold[m - 1])
          hi <- min(min(extract(yStar, z, m)), threshold[m + 1])
          threshold[m] <- runif(1, lo, hi)
        }
      } else {
        for (m in 2:nclass) {
          tmpY <- yStar[-whichNa]
          tmpZ <- z[-whichNa]
          lo <- max(max(extract(tmpY, tmpZ, m - 1)), threshold[m - 1])
          hi <- min(min(extract(tmpY, tmpZ, m)), threshold[m + 1])
          threshold[m] <- runif(1, lo, hi)
        }
      }
      e <- yStar - yHat
      varE <- rep(1,2)
    }
  
    
    # Output the samples - ALL
    # Can be used to assess convergence and to estimate Monte Carlo error
    if(saveChain) {
      
      # Fixed effects
      write(FIXED.list$b, file = FIXED.list$fileOut, append = TRUE)
      
      # Pi and sigma2_b
      if (!(par.app == "none")){
        for(k in 1:nmg){
          if(fixed.prior=="NormalMixture"){
            write(c(MARKER.list[[k]]$probLarge, MARKER.list[[k]]$varB), ncolumns = 2, file = MARKER.list[[k]]$fileOut, append = TRUE)
          } else { # BL
            write(MARKER.list[[k]]$lambda, ncolumns = 1, file = MARKER.list[[k]]$fileOut, append = TRUE)
          }
          
        }
      }
      
      # sigma2_g (or elements of Sigma_g)
      if (!(cov.app == "none")){
        if(cov.app == "stratified"){
          for(k in 1:2){
            write(G.list[[k]]$varG, ncolumns = G.list[[k]]$nVarGcolumns, file = G.list[[k]]$fileOut, append = TRUE)
          }
        } else {
          write(G.list$varG, ncolumns = G.list$nVarGcolumns, file = G.list$fileOut, append = TRUE)
        }
      }
      
      # sigma2_e
      if(varyingE){
        write(varE, ncolumns = 2, file = fileOutVarE, append = TRUE)
      } else {
        write(varE[1], ncolumns = 1, file = fileOutVarE, append = TRUE)
      }
      
    }
    
    
    if (i %% thin == 0) { 
    
      # Note: The posterior mean is updated at every iteration i s.t. i %% thin == 0  & i > burnIn
      
      # Update the posterior means only if i %% thin == 0 & i>burnIn
      if (i > burnIn) { 
        nSums <- nSums + 1
        nSumsm1 <- (nSums - 1)
        
        FIXED.list$post_b <- (FIXED.list$post_b * nSumsm1 + FIXED.list$b)/nSums # Posterior mean 
        FIXED.list$post_b2 <- (FIXED.list$post_b2 * nSumsm1 + (FIXED.list$b^2))/nSums # Posterior mean of b^2 (to get SD(b) later)
        
        if (!(par.app == "none")){
          for(k in 1:nmg){
            MARKER.list[[k]]$post_b <- (MARKER.list[[k]]$post_b * nSumsm1 + MARKER.list[[k]]$b)/nSums
            MARKER.list[[k]]$post_b2 <- (MARKER.list[[k]]$post_b2 * nSumsm1 + (MARKER.list[[k]]$b^2))/nSums
            
            if(fixed.prior=="NormalMixture"){
              MARKER.list[[k]]$post_varB <- (MARKER.list[[k]]$post_varB * nSumsm1 + MARKER.list[[k]]$varB)/nSums
              MARKER.list[[k]]$post_varB2 <- (MARKER.list[[k]]$post_varB2 * nSumsm1 + (MARKER.list[[k]]$varB^2))/nSums
              MARKER.list[[k]]$post_gammaInd <- (MARKER.list[[k]]$post_gammaInd * nSumsm1 + MARKER.list[[k]]$gammaInd)/nSums
              MARKER.list[[k]]$post_probLarge <- (MARKER.list[[k]]$post_probLarge * nSumsm1 + MARKER.list[[k]]$probLarge)/nSums
              MARKER.list[[k]]$post_probLarge2 <- (MARKER.list[[k]]$post_probLarge2 * nSumsm1 + (MARKER.list[[k]]$probLarge)^2)/nSums
            } else { # BL
              MARKER.list[[k]]$post_lambda <- (MARKER.list[[k]]$post_lambda * nSumsm1 + MARKER.list[[k]]$lambda)/nSums
              MARKER.list[[k]]$post_lambda2 <- (MARKER.list[[k]]$post_lambda2 * nSumsm1 + (MARKER.list[[k]]$lambda^2))/nSums
              MARKER.list[[k]]$post_tau2 <- (MARKER.list[[k]]$post_tau2 * nSumsm1 + MARKER.list[[k]]$tau2)/nSums
            }
            
          }
        }

        if (!(cov.app == "none")){
          if(cov.app == "stratified"){
            for(k in 1:2){
              G.list[[k]]$post_g <- (G.list[[k]]$post_g * nSumsm1 + G.list[[k]]$g)/nSums
              G.list[[k]]$post_g2 <- (G.list[[k]]$post_g2 * nSumsm1 + (G.list[[k]]$g^2))/nSums
              G.list[[k]]$post_varG <- (G.list[[k]]$post_varG * nSumsm1 + G.list[[k]]$varG)/nSums
              G.list[[k]]$post_varG2 <- (G.list[[k]]$post_varG2 * nSumsm1 + (G.list[[k]]$varG^2))/nSums 
            }
          } else if(cov.app == "pooled"){
            G.list$post_g <- (G.list$post_g * nSumsm1 + G.list$g)/nSums
            G.list$post_g2 <- (G.list$post_g2 * nSumsm1 + (G.list$g^2))/nSums
            G.list$post_varG <- (G.list$post_varG * nSumsm1 + G.list$varG)/nSums
            G.list$post_varG2 <- (G.list$post_varG2 * nSumsm1 + (G.list$varG^2))/nSums 
          } else {
            G.list$post_g <- (G.list$post_g * nSumsm1 + G.list$g)/nSums
            G.list$post_g2 <- (G.list$post_g2 * nSumsm1 + (G.list$g^2))/nSums
            G.list$post_varG1 <- (G.list$post_varG1 * nSumsm1 + G.list$varG[1])/nSums
            G.list$post_varG2 <- (G.list$post_varG2 * nSumsm1 + G.list$varG[2])/nSums
            G.list$post_rhoG <- (G.list$post_rhoG * nSumsm1 + G.list$varG[3])/nSums 
            G.list$post_varG12 <- (G.list$post_varG12 * nSumsm1 + (G.list$varG[1]^2))/nSums 
            G.list$post_varG22 <- (G.list$post_varG22 * nSumsm1 + (G.list$varG[2]^2))/nSums  
            G.list$post_rhoG2 <- (G.list$post_rhoG2 * nSumsm1 + (G.list$varG[3]^2))/nSums  
          }

        }
        
        post_varE <- (post_varE * nSumsm1 + varE)/nSums
        post_varE2 <- (post_varE2 * nSumsm1 + (varE^2))/nSums
        
        post_yHat <-(post_yHat * nSumsm1 + yHat)/nSums
        post_yHat2 <- (post_yHat2 * nSumsm1 + (yHat^2))/nSums
        
        
        if (response.type == "ordinal") {
          post_threshold <- (post_threshold * nSumsm1 + threshold)/nSums
          post_threshold2 <- (post_threshold2 * nSumsm1 + (threshold^2))/nSums
          TMP <- matrix(nrow = n, ncol = nclass, 0)
          TMP[, 1] <- pnorm(threshold[2] - yHat)
          if (nclass > 2) {
            for (m in 2:(nclass - 1)) {
              TMP[, m] <- pnorm(threshold[(m + 1)] - yHat) - rowSums(as.matrix(TMP[, 1:(m - 1)]))
            }
          }
          TMP[, nclass] <- 1 - rowSums(TMP)
          post_prob <- (post_prob * nSumsm1 + TMP)/nSums
          post_prob2 <- (post_prob2 * nSumsm1 + (TMP^2))/nSums
          if (nNa == 0) {
            logLik <- loglik_ordinal(z, yHat, threshold)
          }
          else {
            logLik <- loglik_ordinal(z[-whichNa], yHat[-whichNa], threshold)
          }
        } else {
          
          # Calculate log-likelihood only with observed 
          logLik <- 0
          tmpE <- e
          tmpTX <- tx
          if (nNa > 0) {
            tmpE <- tmpE[-whichNa]
            tmpTX <- tmpTX[-whichNa]
          }
          for(k in 1:2){
            tmpEk <- tmpE[tmpTX==k]
            tmpSD <- sqrt(varE[k])
            logLik <- logLik + sum(dnorm(tmpEk, sd = tmpSD, log = TRUE))
          }
        }
        post_logLik <- (post_logLik * nSumsm1 + logLik)/nSums
        
      }
    }
    
    tmp <- proc.time()[3]
    iterTime <- tmp - time
    tot.time <- tot.time + iterTime
    if (verbose) {
      cat("---------------------------------------------------------------------\n")
      
      if(varyingE){
        cat(c(paste(c("  Iter=", ", Time/Iter=", "s, varE1=", ", varE2=", ", Total time elapsed="), 
                    c(i, round(iterTime, 5), round(varE[1], 3), round(varE[2], 3), round(tot.time, 1)), sep = "")), "s\n")
      } else {
        cat(c(paste(c("  Iter=", ", Time/Iter=", "s, varE=", ", Total time elapsed="), 
                    c(i, round(iterTime, 5), round(varE[1], 3), round(tot.time, 1)), sep = "")), "s\n")
      }
      
    }
    time <- tmp
  } # End of loop  
  
  if(saveChain) {
    close(fileOutVarE)
    close(FIXED.list$fileOut)
    
    if (!(cov.app == "none")){
      if(cov.app=="stratified"){
        close(G.list[[1]]$fileOut)
        close(G.list[[2]]$fileOut)
      } else {
        close(G.list$fileOut)
      }
    }

    if (!(par.app == "none")){
      for(k in 1:nmg){
        close(MARKER.list[[k]]$fileOut)
      }
    }

  }
  
  
  ### Output 
  
  if (!(par.app == "none")){
    if(fixed.prior=="NormalMixture"){
      HPL$betapars<- vector(mode="list", length=nmg)
      for(k in 1:nmg) HPL$betapars[[k]] <- MARKER.list[[k]]$betaPARS
    }
  }

  
  out <- list(y = y0, response.type = response.type, txgroup = tx, whichNa = whichNa, saveAt = saveAt, nIter = nIter, burnIn = burnIn, thin = thin,  
              cov.app = cov.app, par.app = par.app, fixed.prior = fixed.prior,
              prior.pars = HPL)
  
  names(out$y) <- IDs
  
  out$yHat <- post_yHat
  names(out$yHat) <- IDs
  out$SD.yHat <- sqrt(post_yHat2 - (post_yHat^2))
  
  out$residuals <- as.vector(e)
  names(out$residuals) <- IDs
  out$varE <- post_varE
  out$SD.varE <- sqrt(post_varE2 - post_varE^2)
  
  out$fit <- vector(mode="list", length=4)
  if (response.type == "gaussian") {
    
    logLik <- 0
    tmpE <- yStar - post_yHat
    tmpTX <- tx
    if (nNa > 0) {
      tmpE <- tmpE[-whichNa]
      tmpTX <- tmpTX[-whichNa]
    }
    for(k in 1:2){
      tmpEk <- tmpE[tmpTX==k]
      tmpSD <- sqrt(post_varE[k])
      logLik <- logLik + sum(dnorm(tmpEk, sd = tmpSD, log = TRUE))
    }
    out$fit$logLikAtPostMean <- logLik
    if (Censored) {
      cdfA <- pnorm(q = a[whichNa], sd = sqrt(post_varE), mean = post_yHat[whichNa])
      cdfB <- pnorm(q = b[whichNa], sd = sqrt(post_varE), mean = post_yHat[whichNa])
      out$fit$logLikAtPostMean <- out$fit$logLikAtPostMean + sum(log(cdfB - cdfA))
    }
  } else {
    out$fit$logLikAtPostMean <- loglik_ordinal(y, post_yHat, post_threshold)
    out$probs <- post_prob
    out$SD.probs <- sqrt(post_prob2 - post_prob^2)
    colnames(out$probs) <- lev
    colnames(out$SD.probs) <- lev
    out$threshold <- post_threshold[-c(1, nclass + 1)]
    out$SD.threshold <- sqrt(post_threshold2 - post_threshold^2)[-c(1, nclass + 1)]
    out$levels <- lev
    out$nlevels <- nclass
  }
  
  out$fit$postMeanLogLik <- post_logLik
  out$fit$pD <- -2 * (post_logLik - out$fit$logLikAtPostMean)
  out$fit$DIC <- out$fit$pD - 2 * post_logLik
  
  out$FIXED <- list(b = FIXED.list$post_b, SD.b = sqrt(FIXED.list$post_b2 - FIXED.list$post_b^2))
  
  if (!(par.app == "none")){
    out$MARKER <- vector(mode = "list", length = nmg)
    for(k in 1:nmg){
      if(fixed.prior=="NormalMixture"){
        out$MARKER[[k]] <- list(txName= MARKER.list[[k]]$txName, b = MARKER.list[[k]]$post_b, SD.b = sqrt(MARKER.list[[k]]$post_b2 - MARKER.list[[k]]$post_b^2),
                                varB = MARKER.list[[k]]$post_varB, SD.varB = sqrt(MARKER.list[[k]]$post_varB2 - MARKER.list[[k]]$post_varB^2),
                                PIP = MARKER.list[[k]]$post_gammaInd, 
                                probLarge = MARKER.list[[k]]$post_probLarge, SD.probLarge = sqrt(MARKER.list[[k]]$post_probLarge2 - MARKER.list[[k]]$post_probLarge^2)) 
      } else { # BL
        out$MARKER[[k]] <- list(txName= MARKER.list[[k]]$txName, b = MARKER.list[[k]]$post_b, SD.b = sqrt(MARKER.list[[k]]$post_b2 - MARKER.list[[k]]$post_b^2),
                                tau2 = MARKER.list[[k]]$post_tau2,
                                lambda = MARKER.list[[k]]$post_lambda, SD.lambda = sqrt(MARKER.list[[k]]$post_lambda2 - MARKER.list[[k]]$post_lambda^2))
      }
    }
  }

  if (!(cov.app == "none")){
    if(cov.app == "stratified"){
      out$G <- list(g = numeric(n), SD.g = numeric(n),
                    G1 = list(g = G.list[[1]]$post_g, SD.g = sqrt(G.list[[1]]$post_g2 - G.list[[1]]$post_g^2), 
                              varG = G.list[[1]]$post_varG, 
                              SD.varG = sqrt(G.list[[1]]$post_varG2 - G.list[[1]]$post_varG^2)),
                    G2 = list(g = G.list[[2]]$post_g, SD.g = sqrt(G.list[[2]]$post_g2 - G.list[[2]]$post_g^2), 
                              varG = G.list[[2]]$post_varG, 
                              SD.varG = sqrt(G.list[[2]]$post_varG2 - G.list[[2]]$post_varG^2)), tolDg = tolDg)
      out$G$g[tx==1] <- G.list[[1]]$post_g; 
      out$G$g[tx==2] <- G.list[[2]]$post_g;
      out$G$SD.g[tx==1] <- sqrt(G.list[[1]]$post_g2 - G.list[[1]]$post_g^2)
      out$G$SD.g[tx==2] <- sqrt(G.list[[2]]$post_g2 - G.list[[2]]$post_g^2)
    } else if(cov.app == "pooled"){ 
      out$G <- list(g = G.list$post_g, SD.g = sqrt(G.list$post_g2 - G.list$post_g^2),
                    varG = G.list$post_varG, SD.varG = sqrt(G.list$post_varG2 - G.list$post_varG^2),tolDg = tolDg)
    } else {
      out$G <- list(g = G.list$post_g, SD.g = sqrt(G.list$post_g2 - G.list$post_g^2),
                    varG1 = G.list$post_varG1, SD.varG1 = sqrt(G.list$post_varG12 - G.list$post_varG1^2),
                    varG2 = G.list$post_varG2, SD.varG2 = sqrt(G.list$post_varG22 - G.list$post_varG2^2),
                    rhoG = G.list$post_rhoG, SD.rhoG = sqrt(G.list$post_rhoG2 - G.list$post_rhoG^2), tolDg = tolDg)  
    }
  }

  
  if (verbose) {
    finish(tot.time, nIter)
  }
  
  out$Performance <- list(totalIterTime=tot.time, AveIterTime=tot.time/nIter)
  
  
  return(out)
}


