
source('BWGP.RCT.R')
source('Supplement.BWGP.RCT.R')


load("BWGP.RCT_demo_data.RData")

y <- OSD$y
tx <- OSD$tx
X <- OSD$X
K <- OSD$K


n <- dim(X)[1]
p <- dim(X)[2]
mx <- sum(colSums(X*X))/(n*p)


# Initial variance proportions
h <- 0.25
s.g <- 0.90

sg <- h*s.g
sb <- h - sg

Sb0star <- sb
Sg0star <- sg/dim(X.marker)[2]
rhog0star <- 0.9


data.fixed <- NULL

probLarge <- 0.001
probLargeLIMS <- c(1, 100)
probLargeLIMS <- probLargeLIMS/p
probLargeTYPE <- "beta"



#### Bayes C (Selection only)

prefix <- "DEMO_BayesC_"

dfb0 <- 4.2
Sb0 <- (Sb0star/(mx*probLarge*p))*(dfb0+2)/dfb0

HPL <- list(S0=1, df0=4.2, Sb0=Sb0, dfb0= dfb0)

fm <- BWGP.RCT(y = y, tx = tx, data.fixed=data.fixed, X.marker = X, 
               cov.app="none", par.app="group", 
               int.cond=NULL, varyingE = TRUE,
               nIter=nIter, burnIn=burnIn, thin=thin,
               saveAt = prefix, saveChain=TRUE,
               HPL=HPL, 
               probLarge=probLarge, probLargeTYPE=probLargeTYPE, probLargeLIMS=probLargeLIMS,
               verbose = TRUE) 
save(fm,file=paste(prefix,"fm.rda",sep=""))





#### Stratified BWGP (selection + g / Stratified)

  prefix <- "DEMO_Stratified_"
  
  tolDg <- 1e-6
  
  dfb0 <- 4.2
  Sb0 <- (Sb0star/(mx*probLarge*p))*(dfb0+2)/dfb0
  
  dfg0 <- 4.2
  Sg0 <- Sg0star*(dfg0+2)/dfg0
  
  HPL <- list(S0=1, df0=4.2, Sb0=Sb0, dfb0= dfb0, Sg0=Sg0, dfg0=dfg0)
  
  
  fm <- BWGP.RCT(y = y, tx = tx, data.fixed=data.fixed, X.marker = X.marker, K = K,
                 cov.app="stratified", par.app="group", 
                 int.cond=NULL, varyingE = TRUE,
                 nIter=nIter, burnIn=burnIn, thin=thin,
                 saveAt = saveAt, saveChain=TRUE,
                 HPL=HPL, 
                 probLarge=probLarge, probLargeTYPE=probLargeTYPE, probLargeLIMS=probLargeLIMS,
                 tolDg = tolDg, verbose =TRUE) 
  save(fm,file=paste(saveAt,"fm.rda",sep=""))
  
  
  
#### Adaptive BWGP (selection + g / adaptive)
  
  saveAt <- "DEMO_Adaptive_"
  
  tolDg <- 1e-6
  
  dfb0 <- 4.2
  Sb0 <- 10*(Sb0star/(mx*probLarge*p))*(dfb0+2)/dfb0
  
  dfg0 <- 4.2
  Sg0 <- Sg0star*(dfg0+2)/dfg0
  Sg0 <- matrix(c(Sg0,Sg0*rhog0star,Sg0*rhog0star,Sg0),nrow=2)
  
  HPL <- list(S0=1, df0=4.2, Sb0=Sb0, dfb0= dfb0, Sg0=Sg0, dfg0=dfg0)
  
  
  fm <- BWGP.RCT(y = y, tx = tx, data.fixed=data.fixed, X.marker = X.marker, K = K,
                 cov.app="adaptive", par.app="group", 
                 int.cond=NULL, varyingE = TRUE,
                 nIter=nIter, burnIn=burnIn, thin=thin,
                 saveAt = saveAt, saveChain=TRUE,
                 HPL=HPL, 
                 probLarge=probLarge, probLargeTYPE=probLargeTYPE, probLargeLIMS=probLargeLIMS,
                 tolDg = tolDg, verbose =TRUE, get.cfs = get.cfs) 
  save(fm,file=paste(saveAt,"fm.rda",sep=""))
  
  


