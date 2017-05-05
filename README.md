# BWGP-RCT
Bayesian Whole Genome Prediction for Randomized Clinical Trials

y
(numeric, n) the data-vector (NAs allowed).

data.fixed
(data.frame, n x d) data frame containing non-genetic covariates - should not include the treatment variable (NAs not allowed). Can be NULL. The fixed model will be constructed based on the specified parameterization approach (will include tx interactions). Interactions between variables (other than tx) should be included in data.fixed as a seperate variable to be included in the model. Entry order should match with tx and K.

tx
(integer/numeric, n) vector of treatment indicators (only 0 and 1s allowed: 0: control, 1:tx) - (NAs not allowed). Entry order should match with data.fixed and K.

X.marker
(matrix, n x J) the genotype matrix (coded as 0,1,2: count of minor alleles for each SNP) (NAs not allowed). Marker names are suggested to be stored in dimnames(X.marker)[[2]]. Entry order should match with data.fixed, tx and K.

K
(matrix, n x n) the genetic relationship matrix (GRM) (NAs not allowed). Entry order should match with data.fixed, tx and data.fixed.

par.app
(string) The parameterization approach. Can take values "group" for the group-specific approach; "interaction" for the interaction approach (specify int.cond if so); or "common" for common effects only. Note: The input is case sensitive.

int.cond
(string) Interaction approach assumption ("independent" or "dependent"). Dependent approach assumes a hierarchical relationship between main and interaction effects (main effects are included in the model if corresponding interaction is included).

fixed.prior
(string) The prior model for fixed genetic effects (two approaches implemented: "NormalMixture" and "BL" (Bayesian LASSO). 

cov.app
(string) The covariance modeling approach. Can take values c("stratified","adaptive","pooled","none"). Note: The input is case sensitive.

varyingE
(logical) if TRUE the error variance is allowed to vary between tx groups. Defaults to TRUE.

nIter,burnIn, thin
(integer) the number of iterations, burn-in and thinning.

saveAt
(string) this may include a path and a pre-fix that will be added to the name of the files that are saved as the program runs.

HPL
(list) A list containing hyperparameter values. Note that the list contents depend on parameter (par.app) and covariance (cov.app) approaches. Note: group 1=control, group2=tx.
       For par.app = "group" and cov.app = "adaptive":
       list(S0, df0, Sb=c(Sb1, Sb2), dfb=c(dfb1, dfb2), Sg=diag(Sg1,Sg2), dfg=4.2)

----------------------------------------------------------------------------------------------------
HPL list elements are described below:

S0, df0
(numeric) The scale parameter for the scaled inverse-chi squared prior assigned to the residual variance. The default value for the df parameter is 4.2.

Sb, dfb
(numeric) The scale parameter for the scaled inverse-chi squared prior assigned to the marker effect variances. The default value for the df parameter is 4.2.

Sg, dfg  (for naive/stratified approach)
(numeric) The scale parameter for the scaled inverse-chi squared prior assigned to the genetic variance component. The default value for the df parameter is 4.2.

Sg, dfg  (for adaptive approach)
(matrix, numeric) The scale parameter for the inverse-Wishart prior assigned to the genetic variance component. Recommended to be set to (dfg + k + 1)*Sigma_prior, where k=2 in this case (the # of tx groups) and Sigma_prior is an a priori value for the genetic covariance matrix. Can be set as a diagonal matrix.
The default value for the df parameter is 4.2.
----------------------------------------------------------------------------------------------------

prob.Large
(numeric, 0<prob.Large<1) Initial proportion markers with large effects when fixed.prior="NormalMixture", except when probLargeTYPE="fixed". When probLargeTYPE="fixed", its the given a priori fixed proportion of markers with large effects.

probLargeTYPE
(character) Can take values "fixed", "uniform" (default), "beta" (with parameters given in probLargePARS), or "loguniform" (development stage), representing the type of prior on the inclusion probabilities.

probLarge.betaPARS
(numeric, 2) Parameters of the prior on the inclusion probability when probLargeTYPE="beta". If NULL (default), probLargeLIMS is used to find the optimum parameters given prob.Large and probLargeLIMS (see documentation for details).

probLargeLIMS
(numeric, 2) Expected interval of the inclusion probability to be used to calculate probLarge.betaPARS when probLargeTYPE="beta". probLargeLIMS is used to find the optimum parameters given prob.Large and probLargeLIMS (see documentation for details). The default interval is c(1/dim(X.marker)[2], (dim(X.marker)[2]-1)/dim(X.marker)[2]).

lambda
(numeric, lambda>0) Initial value of the regularization parameter for the Bayesian LASSO prior when fixed.prior="BL".

lambdaTYPE
(character) Can take values "fixed" or "gamma" (default).

tolDg
(numeric) In generating the genetic random effects (g), only eigenvectors of K whose eigenvalues> tolD are kept.

verbose
(logical) if TRUE the iteration history is printed. Defaults to TRUE.

saveChain
(logical) if TRUE the MCMC samples are saved to file. Defaults to TRUE.
