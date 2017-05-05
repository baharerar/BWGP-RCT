#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rconfig.h>
#include <R_ext/Lapack.h>


/*
### Sample g  - C version ########################################################################################################################################                 
# Function to sample random effects sequentially.
*/

SEXP sample_G(SEXP u, SEXP IVe, SEXP Bdiag, SEXP X, SEXP n, SEXP p)
{
  int i,rows,cols;
  
  double m, SigmaU;
  double Bu;
  double *pX, *pBdiag, *pu, *pIVe;
  double *xi;
  int inc=1;
  SEXP list;
  
  cols=INTEGER_VALUE(p);
  rows=INTEGER_VALUE(n);
  
  PROTECT(X=AS_NUMERIC(X));
  pX=NUMERIC_POINTER(X);
  
  PROTECT(IVe=AS_NUMERIC(IVe));
  pIVe=NUMERIC_POINTER(IVe);
  
  PROTECT(u=AS_NUMERIC(u));
  pu=NUMERIC_POINTER(u);
  
  PROTECT(Bdiag=AS_NUMERIC(Bdiag));
  pBdiag=NUMERIC_POINTER(Bdiag);
  
  GetRNGstate();
  
  for(i=0; i<rows; i++)
  {
    xi=pX+i*cols; // Bbig[i,]
    pu[i]=0; // so that Bbig[i,-i] %*% u[-i] = Bbig %*% u where u[i]=0
    
    Bu=F77_NAME(ddot)(&cols,pu,&inc,xi,&inc);
    
    m=(pIVe[i]-Bu)/pBdiag[i];
    SigmaU = 1/pBdiag[i];
    
    pu[i] = m + sqrt(SigmaU)*norm_rand();
  }
  
  // Creating a list with 1 vector element
  PROTECT(list = allocVector(VECSXP,1));
  
  // attaching b vector to list
  SET_VECTOR_ELT(list, 0, u);
  
  PutRNGstate();
  
  UNPROTECT(5);
  
  return(list);
  
}


/*
### Sample Betas - MARKERS - C version ########################################################################################################################################                 
# Function to sample betas for the marker effects.
*/

SEXP sample_beta_MARKER(SEXP b, SEXP gammaInd, SEXP probLarge,SEXP varB, 
                        SEXP X, SEXP x2, SEXP n, SEXP p, SEXP e, SEXP varE)
{
  int j,rows,cols;
  double sigma2e, pvarB, probL, pGamma, tmp,betaj;
  
  double rhs,c,rest;
  double Xe;
  double *pX, *pe, *pb, *px2;
  double *xj;
  int inc=1;
  double c1;
  int *pgammaInd;
  int gammaInd0;
  SEXP list;
  
  cols=INTEGER_VALUE(p);
  rows=INTEGER_VALUE(n);
  sigma2e=NUMERIC_VALUE(varE);
  pvarB=NUMERIC_VALUE(varB);
  probL=NUMERIC_VALUE(probLarge);
  
  PROTECT(X=AS_NUMERIC(X));
  pX=NUMERIC_POINTER(X);
  
  PROTECT(x2=AS_NUMERIC(x2));
  px2=NUMERIC_POINTER(x2);
  
  PROTECT(gammaInd=AS_INTEGER(gammaInd));
  pgammaInd=INTEGER_POINTER(gammaInd);
  
  PROTECT(b=AS_NUMERIC(b));
  pb=NUMERIC_POINTER(b);
  
  PROTECT(e=AS_NUMERIC(e));
  pe=NUMERIC_POINTER(e);
  
  GetRNGstate();
  
  c1=0.5/sigma2e;
  
  for(j=0; j<cols; j++)
  {
    xj=pX+j*rows;
    Xe=F77_NAME(ddot)(&rows,pe,&inc,xj,&inc);
    
    if(pgammaInd[j])
    {
      rest = 2*Xe*pb[j] + pb[j]*pb[j]*px2[j];
    }else{
      rest = 2*Xe*pb[j] - pb[j]*pb[j]*px2[j];
    }
    
    tmp = exp(c1*rest);    
    pGamma=probL*tmp/(1-probL + probL*tmp);
    
    gammaInd0=pgammaInd[j];
    
    pgammaInd[j] = unif_rand()<pGamma ? 1 : 0;
    
    
    
    //Update residuals
    if(gammaInd0!=pgammaInd[j])
    {
      if(pgammaInd[j]>gammaInd0)
      {         
        betaj=-pb[j];
        F77_NAME(daxpy)(&rows, &betaj,xj,&inc, pe, &inc);
        Xe=F77_NAME(ddot)(&rows,pe,&inc,xj,&inc);
      }else{
        betaj=pb[j];
        F77_NAME(daxpy)(&rows, &betaj,xj,&inc, pe, &inc);
      }
    }
    
    //Sample the coefficients
    if(pgammaInd[j]==0)
    {
      //Sample from the prior
      pb[j]=sqrt(pvarB)*norm_rand();
      //pb[j]=0;
    }else{
      //Sampling from the conditional
      rhs=(px2[j]*pb[j] + Xe)/sigma2e;
      c=px2[j]/sigma2e + 1.0/pvarB;
      tmp=rhs/c + sqrt(1.0/c)*norm_rand();
      
      betaj=pb[j]-tmp;
      F77_NAME(daxpy)(&rows, &betaj,xj,&inc, pe, &inc);
      pb[j]=tmp;
    }
    
  }
  
  // Creating a list with 3 vector elements
  PROTECT(list = allocVector(VECSXP,3));
  
  // attaching b vector to list
  SET_VECTOR_ELT(list, 0,gammaInd);
  
  // attaching error vector to list
  SET_VECTOR_ELT(list, 1, e);
  
  // attaching b to the list
  SET_VECTOR_ELT(list,2,b);
  
  PutRNGstate();
  
  UNPROTECT(6);
  
  return(list);
  
}


/*
### Sample Betas - MARKERS - C version ########################################################################################################################################                 
# Function to sample betas for the marker effects when gammaInd_main depends on gammaInd_interaction.
*/

SEXP sample_beta_MARKER_DEP(SEXP b, SEXP gammaInd, SEXP probLarge,SEXP varB, 
                        SEXP X, SEXP x2, SEXP n, SEXP p, SEXP e, SEXP varE, SEXP gammaIndINT, SEXP isInt)
{
  int j,rows,cols;
  double sigma2e, pvarB, probL, pGamma, tmp,betaj;
  
  double rhs,c,rest;
  double Xe;
  double *pX, *pe, *pb, *px2;
  double *xj;
  int inc=1;
  double c1;
  int *pgammaInd;
  int *pgammaIndINT;
  int pisInt;
  int gammaInd0;
  SEXP list;
  
  cols=INTEGER_VALUE(p);
  rows=INTEGER_VALUE(n);
  sigma2e=NUMERIC_VALUE(varE);
  pvarB=NUMERIC_VALUE(varB);
  pisInt=INTEGER_VALUE(isInt);
  probL=NUMERIC_VALUE(probLarge);
  
  PROTECT(X=AS_NUMERIC(X));
  pX=NUMERIC_POINTER(X);
  
  PROTECT(x2=AS_NUMERIC(x2));
  px2=NUMERIC_POINTER(x2);
  
  PROTECT(gammaInd=AS_INTEGER(gammaInd));
  pgammaInd=INTEGER_POINTER(gammaInd);
  
  PROTECT(gammaIndINT=AS_INTEGER(gammaIndINT));
  pgammaIndINT=INTEGER_POINTER(gammaIndINT);
  
  PROTECT(b=AS_NUMERIC(b));
  pb=NUMERIC_POINTER(b);
  
  PROTECT(e=AS_NUMERIC(e));
  pe=NUMERIC_POINTER(e);
  
  GetRNGstate();
  
  c1=0.5/sigma2e;
  
  for(j=0; j<cols; j++)
  {
    xj=pX+j*rows;
    Xe=F77_NAME(ddot)(&rows,pe,&inc,xj,&inc);
    
    if(pgammaInd[j])
    {
      rest = 2*Xe*pb[j] + pb[j]*pb[j]*px2[j];
    }else{
      rest = 2*Xe*pb[j] - pb[j]*pb[j]*px2[j];
    }
    
    tmp = exp(c1*rest);    
    pGamma=probL*tmp/(1-probL + probL*tmp);
    
    gammaInd0=pgammaInd[j];
    
    if(pisInt)
    {
      if(pgammaIndINT[j] == 1)
      {
        pgammaInd[j] = 1;
      }else{
        pgammaInd[j] = unif_rand()<pGamma ? 1 : 0;
      }      
    } else {
      pgammaInd[j] = unif_rand()<pGamma ? 1 : 0;
    }
    
    
    //Update residuals
    if(gammaInd0!=pgammaInd[j])
    {
      if(pgammaInd[j]>gammaInd0)
      {         
        betaj=-pb[j];
        F77_NAME(daxpy)(&rows, &betaj,xj,&inc, pe, &inc);
        Xe=F77_NAME(ddot)(&rows,pe,&inc,xj,&inc);
      }else{
        betaj=pb[j];
        F77_NAME(daxpy)(&rows, &betaj,xj,&inc, pe, &inc);
      }
    }
    
    //Sample the coefficients
    if(pgammaInd[j]==0)
    {
      //Sample from the prior
      pb[j]=sqrt(pvarB)*norm_rand();
      //pb[j]=0;
    }else{
      //Sampling from the conditional
      rhs=(px2[j]*pb[j] + Xe)/sigma2e;
      c=px2[j]/sigma2e + 1.0/pvarB;
      tmp=rhs/c + sqrt(1.0/c)*norm_rand();
      
      betaj=pb[j]-tmp;
      F77_NAME(daxpy)(&rows, &betaj,xj,&inc, pe, &inc);
      pb[j]=tmp;
    }
    
  }
  
  // Creating a list with 3 vector elements
  PROTECT(list = allocVector(VECSXP,3));
  
  // attaching b vector to list
  SET_VECTOR_ELT(list, 0,gammaInd);
  
  // attaching error vector to list
  SET_VECTOR_ELT(list, 1, e);
  
  // attaching b to the list
  SET_VECTOR_ELT(list,2,b);
  
  PutRNGstate();
  
  UNPROTECT(7);
  
  return(list);
  
}



/*
### Estimate Beta_star_g's  ########################################################################################################################################                 
*/

SEXP est_beta_star(SEXP b1, SEXP b2, SEXP varGinv1, SEXP varGinv2, SEXP varGinv3,
                   SEXP X1, SEXP X2, SEXP x12, SEXP x22, SEXP n1, SEXP n2, SEXP p, SEXP e1, SEXP e2, SEXP varE1, SEXP varE2)
{
  int j, rows1, rows2, cols;
  double sigma2e1, sigma2e2, pvarGinv1, pvarGinv2, pvarGinv3;
  
  double c;
  
  double Omega1, Omega2, Omega3;
  double OmegaInv1, OmegaInv2, OmegaInv3;
  
  double xytilde1, xytilde2;
  double Xe1, Xe2;
  double *pX1, *pX2, *pe1, *pe2, *pb1, *pb2, *px12, *px22;
  double *xj1, *xj2;
  
  int inc=1;
  
  SEXP list;
  
  cols=INTEGER_VALUE(p);
  rows1=INTEGER_VALUE(n1);
  rows2=INTEGER_VALUE(n2);
  
  sigma2e1=NUMERIC_VALUE(varE1);
  sigma2e2=NUMERIC_VALUE(varE2);
  pvarGinv1=NUMERIC_VALUE(varGinv1);
  pvarGinv2=NUMERIC_VALUE(varGinv2);
  pvarGinv3=NUMERIC_VALUE(varGinv3);
  
  PROTECT(X1=AS_NUMERIC(X1));
  pX1=NUMERIC_POINTER(X1);
  PROTECT(X2=AS_NUMERIC(X2));
  pX2=NUMERIC_POINTER(X2);
  
  PROTECT(x12=AS_NUMERIC(x12));
  px12=NUMERIC_POINTER(x12);
  PROTECT(x22=AS_NUMERIC(x22));
  px22=NUMERIC_POINTER(x22);
  
  PROTECT(b1=AS_NUMERIC(b1));
  pb1=NUMERIC_POINTER(b1);
  PROTECT(b2=AS_NUMERIC(b2));
  pb2=NUMERIC_POINTER(b2);
  
  PROTECT(e1=AS_NUMERIC(e1));
  pe1=NUMERIC_POINTER(e1);
  PROTECT(e2=AS_NUMERIC(e2));
  pe2=NUMERIC_POINTER(e2);
  
  GetRNGstate();
  
  for(j=0; j<cols; j++)
  {
    xj1=pX1+j*rows1;
    Xe1=F77_NAME(ddot)(&rows1,pe1,&inc,xj1,&inc);
    xytilde1= (Xe1+pb1[j]*px12[j])/sigma2e1;
    
    xj2=pX2+j*rows2;
    Xe2=F77_NAME(ddot)(&rows2,pe2,&inc,xj2,&inc);
    xytilde2= (Xe2+pb2[j]*px22[j])/sigma2e2;
    
    OmegaInv1= (px12[j]/sigma2e1) + cols*pvarGinv1;
    OmegaInv2= (px22[j]/sigma2e2) + cols*pvarGinv2;
    OmegaInv3= cols*pvarGinv3;
    
    c = 1/(OmegaInv1*OmegaInv2*(1-OmegaInv3*OmegaInv3));
    Omega1 = c*OmegaInv2;
    Omega2 = c*OmegaInv1;
    Omega3 = -c*(sqrt(OmegaInv1*OmegaInv2)*OmegaInv3);
      
    pb1[j] = Omega1*xytilde1 + Omega3*xytilde2;
    pb2[j] = Omega3*xytilde1 + Omega2*xytilde2;
    
  }
  
  // Creating a list with 3 vector elements
  PROTECT(list = allocVector(VECSXP,2));
  
  SET_VECTOR_ELT(list, 0, b1);
  SET_VECTOR_ELT(list, 1, b2);
  
  PutRNGstate();
  
  UNPROTECT(9);
  
  return(list);
  
}

