#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include <cstdlib>
#include <stdio.h>
#include "matrix.h"
#include <string.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

/*#include "randomc.h"
#include "stocc.h"
#ifndef MULTIFILE_PROJECT
// If compiled as a single file then include these cpp files,
// If compiled as a project then compile and link in these cpp files.
// Include code for the chosen random number generator:
#include "mersenne.cpp"
#include "stoc1.cpp"
// define system specific user interface:
#include "userintf.cpp"
#endif
*/
#ifndef _NO_NAMESPACE
using namespace std;
using namespace math;
#define STD std
#else
#define STD
#endif

#ifndef _NO_TEMPLATE
typedef matrix<double> Matrix;
#else
typedef matrix Matrix;
#endif

using namespace std;

//////////////////////////////////////////////////////////////////////////////

Matrix Yij;
Matrix Xijkl;
Matrix hijkl;
Matrix pijkl;
double fijkl;
double fikl;
Matrix Sijkl;
Matrix Dijkl;
Matrix Wik;
Matrix Zk;
Matrix Dp;
Matrix SNP;              //score of the fixed effects and the NP mass points
Matrix FNP;              //information matrix for the fixed effects and the NP mass points
Matrix spnp;             //score of the probability vector ... needed if npk > 1
Matrix fpnp;             //information for the probabilities
Matrix fmix;             //mixed second derivatives
Matrix infonpml;         //observed information matrix based of the square of the score function
int flag = 0;
int flagcvm = 0;
int flaginfo = 0;

// A function that makes the entries of a matrix/array equal to zero

void zerof(Matrix &Z, int nr, int nc)
{
     for (register int i = 0; i < nr; i++)
         for (register int j = 0; j < nc; j++)
            Z(i,j) = 0.0;
}

// A function that finds the max of the absolute values of a vector s

double maxabs(Matrix &s, int L)
{
       double maxi = -1;
       double g;
       for (register int i = 0; i < L; i++)
           if ((g=abs(s(i,0))) > maxi) maxi = g;
       return maxi;
}

// (i,j)th response vector

void setYij(int i, int j, int q, int *CN, int *resp)
{
     Matrix tt2(q,1);
     for (register int u=0; u<q; u++) tt2(u,0)=0;
     if (resp[CN[i]+j]<=q) tt2(resp[CN[i]+j]-1,0)=1;
     Yij=tt2;
}

// Design matrix

void setXijkl(int i, int j, int *CN, int q, int np, int npk, int npl, Matrix &EP, double *model, int rslpind, double *rslp, int Ntot, double *npoind, int T)
{
     int ncol = q+np+(npk-1)*(1+rslpind);
     Matrix tt(q,ncol);
     zerof(tt,q,ncol);
     // diagonal qxq for the intercepts
     for (register int k=0; k<q; k++)
         tt(k,k) = 1;
     // count the non proportional odds
     int count=0;
     // npoind is a vector of length equal to the number of predictors (T) that has 1's for non % odds and 0's o/w
     for (register int u=0; u<T; u++){
        for (register int k=0; k<q; k++)
           tt(k,q+u+k*npoind[u]+count*(q-1)) = model[CN[i]+j+u*Ntot];
        count += npoind[u];
     }
     if (npk > 1)                                          // This part accomodates the NP random intercept and slope
     {
             if (npl < npk)
             {
                   for (register int l1 = 0; l1 < q; l1++){
                       tt(l1,q+np+npl-1) = 1.0;
                       if (rslpind>0)
                          for (register int u=0; u<rslpind; u++)
                             tt(l1,q+np+(npk-1)*(u+1)+npl-1) = rslp[CN[i]+j+u*Ntot];
                    }
             }
             if (npl == npk)
             {
                   for (register int l1 = 0; l1 < q; l1++){
                      for (register int jj = 0; jj < (npk-1); jj++){
                         tt(l1,q+np+jj) = -EP(jj,0)/EP(npk-1,0);
                         if (rslpind>0)
                            for (register int u=0; u<rslpind; u++)
                               tt(l1,q+np+(npk-1)*(u+1)+jj) = -rslp[CN[i]+j+u*Ntot]*EP(jj,0)/EP(npk-1,0);}}
             }
     }
     Xijkl=tt;
}

// q-variate linear predictor

void sethijkl(int i, int j, int *CN, int q, int np, int npk, int npl, Matrix &EP, double *model, Matrix &beta, int rslpind, double *rslp, int Ntot, double *npoind, int T)
{
       setXijkl(i, j, CN, q, np, npk, npl, EP, model, rslpind, rslp, Ntot, npoind, T);
       hijkl=Xijkl*beta;
}

// q-variate probability vector

void setpijkl(int i, int j, int *CN, int q, int np, int npk, int npl, Matrix &EP, double *model, Matrix &beta, int link, int rslpind, double *rslp, int Ntot, double *npoind, int T)
{
     sethijkl(i, j, CN, q, np, npk, npl, EP, model, beta, rslpind, rslp, Ntot, npoind, T);
     Matrix R(q,1);
     Matrix Rtemp(q,1);
     // clogit link
     if (link==1){
     for (register int k=0; k< q; k++)
         R(k,0) = exp(hijkl(k,0))/(1+exp(hijkl(k,0)));
     Rtemp(0,0) = R(0,0);
     for (register int k=1; k< q; k++) Rtemp(k,0) = R(k,0) - R(k-1,0);
     pijkl=Rtemp;}
     // blogit link: baseline category is the last category
     if (link==0){
     for (register int k=0; k< q; k++)
         R(k,0) = exp(hijkl(k,0));
     double partot = 1;
     for (register int k=0; k<q; k++) partot += R(k,0);
     for (register int k=0; k<q; k++) Rtemp(k,0)=R(k,0)/partot;
     pijkl=Rtemp;}
}

// Function Sum: Sums the first column of a matrix R

double sum1f(Matrix &R, int q)
{
     double sum2 = 0.0;
     for (int k=0; k < q; k++)
         sum2 += R(k,0);
     return sum2;
}

// Probability of observing what we observed for cluster i at time j

void setfijkl(int i, int j, int *CN, int q, int np, int npk, int npl, Matrix &EP, double *model, Matrix &beta, int link, int rslpind, double *rslp, int Ntot, int *resp, double *npoind, int T)
{
    setpijkl(i, j, CN, q, np, npk, npl, EP, model, beta, link, rslpind, rslp, Ntot, npoind, T);
    setYij(i, j, q, CN, resp);
    double t=1.0;
    for (int i=0; i<q; i++)
        t = t * pow(pijkl(i,0),Yij(i,0));
    t = t * pow(1-sum1f(pijkl,q), 1-sum1f(Yij,q));
    fijkl=t;
}

// Probability of observing what we observed for cluster i

void setfikl(int i, int *CN, int q, int np, int npk, int npl, Matrix &EP, double *model, Matrix &beta, int link, int rslpind, double *rslp, int Ntot, int *resp, double *npoind, int T)
{
    double t = 1.0;
    for (register int j=0; j<(CN[i+1]-CN[i]); j++)
    {
           setfijkl(i,j,CN,q,np,npk,npl,EP,model,beta,link,rslpind,rslp,Ntot,resp,npoind,T);
           t *= fijkl;
    }
    fikl=t;
}

// Sij = covariance matrix based on pij

void setSijkl(int i, int j, int *CN, int q, int np, int npk, int npl, Matrix &EP, double *model, Matrix &beta, int link, int rslpind, double *rslp, int Ntot, double *npoind, int T)
{
     setpijkl(i,j,CN,q,np,npk,npl,EP,model,beta,link,rslpind,rslp,Ntot,npoind,T);
     Matrix K(q,q);
     for (register int k=0; k< q; k++)
     for (register int l=0; l< q; l++)
         if (k == l) K(k,l) = pijkl(l,0); else K(k,l) = 0;
     Sijkl = K - pijkl * (~pijkl);
}

// Dij = der h(h) / der h, where h is the inverse link function

void setDijkl(int i, int j, int *CN, int q, int np, int npk, int npl, Matrix &EP, double *model, Matrix &beta, int link, int rslpind, double *rslp, int Ntot, double *npoind, int T)
{
     Matrix K2(q,q);
     // clogit link
     if (link==1){
     sethijkl(i,j,CN,q,np,npk,npl,EP,model,beta,rslpind,rslp,Ntot,npoind,T);
     zerof(K2,q,q);
     for (register int k=0; k<q; k++)
         K2(k,k) = exp(hijkl(k,0))/pow((1+exp(hijkl(k,0))),2);
     for (register int k=0; k<(q-1); k++)
         K2(k+1,k) = -exp(hijkl(k,0))/pow((1+exp(hijkl(k,0))),2);
     Dijkl = K2;}
     // blogit link
     if (link==0){
     setpijkl(i, j, CN, q, np, npk, npl, EP, model, beta, link,rslpind,rslp,Ntot,npoind,T);
     for (register int k=0; k< q; k++)
     for (register int l=0; l< q; l++)
         if (k == l) K2(k,l) = pijkl(k,0)*(1-pijkl(k,0)); else K2(k,l) = -pijkl(k,0)*pijkl(l,0);
     Dijkl = K2;}
}

///////////////////////// Sij and Dij together

void setSDijkl(int i, int j, int *CN, int q, int np, int npk, int npl, Matrix &EP, double *model, Matrix &beta, int link, int rslpind, double *rslp, int Ntot, double *npoind, int T)
{
  //  hijkl=Xijkl*beta;
     sethijkl(i, j, CN, q, np, npk, npl, EP, model, beta, rslpind, rslp, Ntot, npoind, T);
     ///////////////////////////Probabilities
     Matrix R(q,1);
     Matrix Rtemp(q,1);
     // clogit link
     if (link==1){
     for (register int k=0; k<q; k++)
         R(k,0) = exp(hijkl(k,0))/(1+exp(hijkl(k,0)));
     Rtemp(0,0) = R(0,0);
     for (register int k=1; k< q; k++) Rtemp(k,0) = R(k,0) - R(k-1,0);
     pijkl=Rtemp;}
     // blogit link
     if (link==0){
     for (register int k=0; k< q; k++)
         R(k,0) = exp(hijkl(k,0));
     double partot = 1;
     for (register int k=0; k<q; k++) partot += R(k,0);
     for (register int k=0; k<q; k++) Rtemp(k,0)=R(k,0)/partot;
     pijkl=Rtemp;}
     // Common Sijkl
     Matrix K(q,q);
     for (register int k=0; k< q; k++)
     for (register int l=0; l< q; l++)
         if (k == l) K(k,l) = pijkl(l,0); else K(k,l) = 0;
     Sijkl = K - pijkl * (~pijkl);
     /////////////////////////// Dij
      Matrix K2(q,q);
     // clogit link
     if (link==1){
     sethijkl(i,j,CN,q,np,npk,npl,EP,model,beta,rslpind,rslp,Ntot,npoind,T);
     zerof(K2,q,q);
     for (register int k=0; k<q; k++)
         K2(k,k) = exp(hijkl(k,0))/pow((1+exp(hijkl(k,0))),2);
     for (register int k=0; k<(q-1); k++)
         K2(k+1,k) = -exp(hijkl(k,0))/pow((1+exp(hijkl(k,0))),2);
     Dijkl = K2;}
     // blogit link
     if (link==0){
     setpijkl(i, j, CN, q, np, npk, npl, EP, model, beta, link, rslpind, rslp,Ntot,npoind,T);
     for (register int k=0; k< q; k++)
     for (register int l=0; l< q; l++)
         if (k == l) K2(k,l) = pijkl(k,0)*(1-pijkl(k,0)); else K2(k,l) = -pijkl(k,0)*pijkl(l,0);
     Dijkl = K2;}
}

// Posterior probabilities of group membership

void setWik(int i, int *CN, int q, int np, int npk, Matrix &EP, double *model, Matrix &beta, int link, int rslpind, double *rslp, int Ntot, int *resp, double *npoind, int T)
{
     Matrix Tstore(npk,1);
     Matrix Tstore2(npk,1);
     double tot = 0.0;
     for (register int l = 0; l < npk; l++)
     {
         setfikl(i,CN,q,np,npk,l+1,EP,model,beta,link,rslpind,rslp,Ntot,resp,npoind,T);
         Tstore(l,0) = fikl*EP(l,0);
         tot += Tstore(l,0);
     }
     for (register int l = 0; l < npk; l++)
         Tstore2(l,0) = Tstore(l,0)/tot;
     Wik = Tstore2;
}

//###############################Score and Fisher information matrix

//2 functions for calculating the score and info of the probabilities

void setZk(int i, int j, int *CN, int npk, int q, int np, Matrix &beta, int rslpind, double *rslp, int Ntot)
{
   Matrix Temp(q,npk-1);
   for (int c = 0; c < npk-1; c++){
   for (int r = 0; r < q; r++){
      Temp(r,c) = beta(q+np+c,0);
      if (rslpind>0){
         for (int u=0; u<rslpind; u++){
            Temp(r,c) += rslp[(CN[i]+j+u*Ntot)]*beta(q+np+(npk-1)*(u+1)+c,0);}}}}
   Zk = Temp;
}

void setDp(int npk, Matrix &EP)
{
   Matrix t(npk-1,npk-1);
   for (int r = 0; r < npk-1; r++)
   for (int c = 0; c < npk-1; c++)
      t(r,c) = -EP(r,0)/pow(EP(npk-1,0),2);
   for (int r = 0; r < npk-1; r++)
      t(r,r) = (sum1f(EP,npk-1)-1-EP(r,0))/pow(EP(npk-1,0),2);
   Dp = t;
}

//Updates the score vector and Fisher information matrix

void NPML(int npk, int m, int q, int np, int *CN, Matrix &beta, Matrix &EP, double *model, int link, int rslpind, double *rslp, double tol, int Ntot, int *resp, double *npoind, int T)
{
     flag = 0;            // if there is at least one eigenvalue less than tol, flag will become one for this iteration

     int ncol = q+np+(npk-1)*(1+rslpind);          // for beta*
     Matrix VT(ncol,q);
     Matrix s1a(ncol,1);
     Matrix s2a(ncol,1);
     Matrix s3a(ncol,1);
     Matrix s1b(ncol,ncol);
     Matrix s2b(ncol,ncol);
     Matrix s3b(ncol,ncol);

     zerof(s3a,ncol,1);
     zerof(s3b,ncol,ncol);

     Matrix VT2(q,q);                             // for p
     Matrix sp2(npk-1,1);
     Matrix sp22(npk-1,1);
     Matrix sp222(npk-1,1);
     Matrix fp2(npk-1,npk-1);
     Matrix fp22(npk-1,npk-1);
     Matrix fp222(npk-1,npk-1);

     zerof(sp22,npk-1,1);
     zerof(fp22,npk-1,npk-1);
     zerof(sp222,npk-1,1);
     zerof(fp222,npk-1,npk-1);

     Matrix VT3(q,ncol);                             // mixed derivatives
     Matrix fmix1(npk-1,ncol);
     Matrix fmix2(npk-1,ncol);
     zerof(fmix2,npk-1,ncol);

     // needed for calculating the g-inverse of Sijkl
     gsl_matrix *Ag,*Vg,*DSg,*R1g,*R2g;
     gsl_vector *Sg,*work;
     // the next 4 are for Ag=U S V^T
     Ag=gsl_matrix_alloc(q,q);
     Vg=gsl_matrix_alloc(q,q);
     Sg=gsl_vector_alloc(q);
     work=gsl_vector_alloc(q);
     // used for inverting the nonzero eigenvalues
     Matrix eigen(q,1);
     //to become the diagonal with the inverse of nonzero eigenvalues
     DSg=gsl_matrix_calloc(q,q);
     //to become Vg * DSg
     R1g=gsl_matrix_calloc(q,q);
     //to become Vg * DSg * Vg^T
     R2g=gsl_matrix_calloc(q,q);
     // to write R2g
     Matrix ginvSijkl(q,q);

     for (register int i = 0; i < m; i++)
     {
         zerof(s2a,ncol,1);
         zerof(s2b,ncol,ncol);
         zerof(sp2,npk-1,1);
         zerof(fp2,npk-1,npk-1);
         zerof(fmix1,npk-1,ncol);

         setWik(i, CN, q, np, npk, EP, model, beta, link, rslpind, rslp, Ntot, resp,npoind,T);

         for (int npl = 1; npl <= npk; npl++)
         {
             zerof(s1a,ncol,1);
             zerof(s1b,ncol,ncol);
             for (register int j = 0; j<(CN[i+1]-CN[i]); j++)
             {
                    //setXijkl(i, j, CN, q, np, npk, npl, EP, model, rslpind, rslp,Ntot);
                    setSDijkl(i, j, CN, q, np, npk, npl, EP, model, beta, link, rslpind, rslp,Ntot,npoind,T);
                    setYij(i, j, q, CN, resp);

                   // g-inverse of Sijkl

                  for (int i2=0; i2<q; i2++){                  // Ag becomes Sijkl
                    for (int j2=0; j2<q; j2++) {
                       gsl_matrix_set(Ag,i2,j2,Sijkl(i2,j2));
                    }
                  }
                  gsl_linalg_SV_decomp(Ag,Vg,Sg,work);          // Ag = U S V^T, where Ag becomes U
                  for (int i3=0; i3<q; i3++)                   // Calculate 1/eigen for eigen > 0
		             if (gsl_vector_get(Sg,i3) > tol) eigen(i3,0)=1/gsl_vector_get(Sg,i3); else {eigen(i3,0)=0; flag = 1;}
                  for (int i4=0; i4<q; i4++)                   // Dsg=1/eigen, for eigen > 0
                      gsl_matrix_set(DSg,i4,i4,eigen(i4,0));
                  // R1g = Vg * DSg
                  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, Vg, DSg, 0, R1g);
                  gsl_matrix_transpose(Ag);
                  // R2g = R1g * Ag^T = Vg * DSg * Vg^T thus obtaining the g-inverse
                  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, R1g, Ag, 0, R2g);

                  for (int i5=0; i5<q; i5++)
                     for (int j5=0; j5<q; j5++)
                        ginvSijkl(i5,j5) = gsl_matrix_get(R2g,i5,j5);

                    VT2 = (~Dijkl) * (ginvSijkl);

                    VT = (~Xijkl) * VT2;
                    VT3 = Dijkl * Xijkl;
                    s1a = s1a + VT * (Yij - pijkl);
                    s1b = s1b + VT * VT3;
                    if ((npk == npl) & (npk > 1))
                    {
                         if (rslpind>0) setZk(i, j, CN, npk, q, np, beta, rslpind, rslp, Ntot);
                         if ((rslpind==0) & (i==0) & (j==0)) setZk(i, j, CN, npk, q, np, beta, rslpind, rslp, Ntot);
                         sp2 = sp2 + ~Zk * VT2 * (Yij - pijkl);
                         fp2 = fp2 + ~Zk * VT2 * Dijkl * Zk;
                         fmix1 = fmix1 + ~Zk * (VT2 * VT3);
                    }
             }
             s2a = s2a + Wik(npl-1,0)*s1a;
             s2b = s2b + Wik(npl-1,0)*s1b;
          }
         s3a = s3a + s2a;
         s3b = s3b + s2b;
         if (npk > 1)
         {
                 sp22 = sp22 + Wik(npk-1,0) * sp2;
                 fp22 = fp22 + Wik(npk-1,0) * fp2;
                 for (int a=0; a< npk-1; a++)
                     sp222(a,0) = sp222(a,0) + Wik(a,0)/EP(a,0) - Wik(npk-1,0)/EP(npk-1,0);
                 for (int a=0; a < npk-1; a++)
                 for (int b=0; b < npk-1; b++)
                 {
                     if (a == b) fp222(a,b) = fp222(a,b) + Wik(a,0)/pow(EP(a,0),2) + Wik(npk-1,0)/pow(EP(npk-1,0),2);
                     else fp222(a,b) = fp222(a,b) + Wik(npk-1,0)/pow(EP(npk-1,0),2);
                 }
                 fmix2 = fmix2 + Wik(npk-1,0) * fmix1;
         }
     }
     SNP = s3a;
     FNP = s3b;
     if (npk > 1)
     {
         setDp(npk, EP);
         spnp = ~Dp * sp22 + sp222;
         fpnp = ~Dp * fp22 * Dp + fp222;
         fmix = ~Dp * fmix2;
     }
    // free the matrices needed for calculating the g-inverse
    gsl_matrix_free(Ag);
    gsl_matrix_free(Vg);
    gsl_matrix_free(DSg);
    gsl_matrix_free(R1g);
    gsl_matrix_free(R2g);
    gsl_vector_free(Sg);
    gsl_vector_free(work);
}

// A function that iterates until the NP model is fit

void upbp(int npk, Matrix &beta, Matrix &EP, int m, int q, int np, int *CN, double *model, double eps, Matrix &iter, int link, int maxit, int rslpind, double *rslp, double tol, int Ntot, int *resp, double *npoind, int T)
{
     int nfpar = q+np+(2+rslpind)*(npk-1);
     Matrix ua(nfpar,1);
     Matrix scrll(nfpar,1);
     Matrix FB(nfpar,nfpar);
     scrll(0,0)=1;
     // needed for calculating the g-inverse
     gsl_matrix *Ag,*Vg,*DSg,*R1g,*R2g;
     gsl_vector *Sg,*work;
     // the next 4 are for Ag=U S V^T
     Ag=gsl_matrix_alloc(nfpar,nfpar);
     Vg=gsl_matrix_alloc(nfpar,nfpar);
     Sg=gsl_vector_alloc(nfpar);
     work=gsl_vector_alloc(nfpar);
     // used for inverting the nonzero eigenvalues
     Matrix eigen(nfpar,1);
     //to become the diagonal with the inverse of nonzero eigenvalues
     DSg=gsl_matrix_calloc(nfpar,nfpar);
     //to become Vg * DSg
     R1g=gsl_matrix_calloc(nfpar,nfpar);
     //to become Vg * DSg * Vg^T
     R2g=gsl_matrix_calloc(nfpar,nfpar);
     // to write R2g
     Matrix FbinvA(nfpar,nfpar);

     while ((maxabs(scrll,nfpar) > eps) & (iter(0,0) < maxit))
     {
           NPML(npk,m,q,np,CN,beta,EP,model,link,rslpind,rslp,tol,Ntot,resp,npoind,T);
           iter(0,0)=iter(0,0)+1;
           if (flag==1) flagcvm=iter(0,0);

           for (register int i=0; i < q+np+(npk-1)*(1+rslpind); i++)
               scrll(i,0) = SNP(i,0);
           for (register int i=0; i < npk-1; i++)
               scrll(q+np+(npk-1)*(1+rslpind)+i,0) = spnp(i,0);

           for (register int i = 0; i < q+np+(npk-1)*(1+rslpind); i++)
           for (register int j = 0; j < q+np+(npk-1)*(1+rslpind); j++)
               FB(i,j) = FNP(i,j);
           for (register int i = 0; i < npk-1; i++)
           for (register int j = 0; j < npk-1; j++)
               FB(q+np+(npk-1)*(1+rslpind)+i,q+np+(npk-1)*(1+rslpind)+j) = fpnp(i,j);
           for (register int i = 0; i < npk-1; i++)
           for (register int j = 0; j < q+np+(npk-1)*(1+rslpind); j++)
           {
               FB(q+np+(npk-1)*(1+rslpind)+i,j) = fmix(i,j);
               FB(j,q+np+(npk-1)*(1+rslpind)+i) = fmix(i,j);
           }
           // g-inverse of FB(nfpar,nfpar);

           for (int i=0; i<nfpar; i++){                  // Ag becomes FB
               for (int j=0; j<nfpar; j++) {
                       gsl_matrix_set(Ag,i,j,FB(i,j));
               }
           }
           gsl_linalg_SV_decomp(Ag,Vg,Sg,work);          // Ag = U S V^T, where Ag becomes U
           for (int i=0; i<nfpar; i++)                   // Calculate 1/eigen for eigen > 0
	          if (gsl_vector_get(Sg,i) > tol) eigen(i,0)=1/gsl_vector_get(Sg,i); else {eigen(i,0)=0; flaginfo=iter(0,0);}
           for (int i=0; i<nfpar; i++)                   // Dsg=1/eigen, for eigen > 0
               gsl_matrix_set(DSg,i,i,eigen(i,0));
           // R1g = Vg * DSg
           gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, Vg, DSg, 0, R1g);
           gsl_matrix_transpose(Ag);
           // R2g = R1g * Ag^T = Vg * DSg * Vg^T thus obtaining the g-inverse
           gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, R1g, Ag, 0, R2g);

           for (int i=0; i<nfpar; i++)
               for (int j=0; j<nfpar; j++)
                   FbinvA(i,j) = gsl_matrix_get(R2g,i,j);

            ua = FbinvA * scrll;
           //ua = !FB * scrll; inverse replaced by g-inverse

           for (register int i=0; i < q+np+(npk-1)*(1+rslpind); i++)
               beta(i,0) = beta(i,0) + ua(i,0);
           for (register int i=0; i < npk-1; i++)
               EP(i,0) = EP(i,0) + ua(q+np+(npk-1)*(1+rslpind)+i,0);
           EP(npk-1,0) = 1 - sum1f(EP,npk-1);
     }
     infonpml=FbinvA;              // not the information but the covariance matrix
     // free the matrices needed for calculating the g-inverse
     gsl_matrix_free(Ag);
     gsl_matrix_free(Vg);
     gsl_matrix_free(DSg);
     gsl_matrix_free(R1g);
     gsl_matrix_free(R2g);
     gsl_vector_free(Sg);
     gsl_vector_free(work);
}

// Finding the SE for the NPML -obs mat

void NPMLSE(int npk, int m, int q, int np, int *CN, Matrix &beta, Matrix &EP, double *model, int link, int rslpind, double *rslp, double tol, int Ntot, int *resp, double *npoind, int T)
{
     int npar = q+np+(npk-1)*(1+rslpind);                    // fixed effects and mass points
     int npar2 = q+np+(2+rslpind)*(npk-1);                   // all parameters
     Matrix s1a(npar,1);
     Matrix s2a(npar,1);
     Matrix pse(npk-1,1);
     Matrix pse2(npk-1,1);
     Matrix VT3(q,q);
     Matrix sp2(npk-1,1);
     Matrix sec(npar2,1);
     Matrix sec2(npar2,npar2);
     zerof(sec2,npar2,npar2);

     // needed for calculating the g-inverse of Sijkl
     gsl_matrix *Ag,*Vg,*DSg,*R1g,*R2g;
     gsl_vector *Sg,*work;
     // the next 4 are for Ag=U S V^T
     Ag=gsl_matrix_alloc(q,q);
     Vg=gsl_matrix_alloc(q,q);
     Sg=gsl_vector_alloc(q);
     work=gsl_vector_alloc(q);
     // used for inverting the nonzero eigenvalues
     Matrix eigen(q,1);
     //to become the diagonal with the inverse of nonzero eigenvalues
     DSg=gsl_matrix_calloc(q,q);
     //to become Vg * DSg
     R1g=gsl_matrix_calloc(q,q);
     //to become Vg * DSg * Vg^T
     R2g=gsl_matrix_calloc(q,q);
     // to write R2g
     Matrix ginvSijkl(q,q);

     for (register int i = 0; i < m; i++)
     {
         zerof(pse,npk-1,1);
         zerof(s2a,npar,1);
         zerof(sp2,npk-1,1);
         setWik(i, CN, q, np, npk, EP, model, beta, link, rslpind, rslp, Ntot, resp, npoind, T);
         for (register int npl = 1; npl <= npk; npl++)
         {
             if (npl < npk) pse(npl-1,0) =  Wik(npl-1,0)/EP(npl-1,0) - Wik(npk-1,0)/EP(npk-1,0);
             zerof(s1a,npar,1);
             for (register int j = 0; j<(CN[i+1]-CN[i]); j++)
             {
                    //setXijkl(i, j, CN, m, q, np, npk, npl, EP, model);
                    setSDijkl(i, j, CN, q, np, npk, npl, EP, model, beta, link, rslpind, rslp, Ntot, npoind, T);
                    setYij(i, j, q, CN, resp);

                    // g-inverse of Sijkl

                    for (int i2=0; i2<q; i2++){                  // Ag becomes Sijkl
                       for (int j2=0; j2<q; j2++) {
                          gsl_matrix_set(Ag,i2,j2,Sijkl(i2,j2));
                       }
                    }
                    gsl_linalg_SV_decomp(Ag,Vg,Sg,work);          // Ag = U S V^T, where Ag becomes U
                    for (int i3=0; i3<q; i3++)                   // Calculate 1/eigen for eigen > 0
                        if (gsl_vector_get(Sg,i3) > tol) eigen(i3,0)=1/gsl_vector_get(Sg,i3); else eigen(i3,0)=0;
                    for (int i4=0; i4<q; i4++)                   // Dsg=1/eigen, for eigen > 0
                        gsl_matrix_set(DSg,i4,i4,eigen(i4,0));
                    // R1g = Vg * DSg
                    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, Vg, DSg, 0, R1g);
                    gsl_matrix_transpose(Ag);
                    // R2g = R1g * Ag^T = Vg * DSg * Vg^T thus obtaining the g-inverse
                    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, R1g, Ag, 0, R2g);

                    for (int i5=0; i5<q; i5++)
                       for (int j5=0; j5<q; j5++)
                          ginvSijkl(i5,j5) = gsl_matrix_get(R2g,i5,j5);

                    VT3 = (~Dijkl) * (ginvSijkl);

                    s1a = s1a + (~Xijkl) * VT3 * (Yij - pijkl);
                    if ((npk == npl) & (npk > 1))
                    {
                         if (rslpind>0) setZk(i, j, CN, npk, q, np, beta, rslpind, rslp, Ntot);
                         if ((rslpind==0) & (i==0) & (j==0)) setZk(i, j, CN, npk, q, np, beta, rslpind, rslp, Ntot);
                         sp2 = sp2 + ~Zk * VT3 * (Yij - pijkl);
                    }
             }
             s2a = s2a + Wik(npl-1,0)*s1a;
         }

         if (npk > 1)
         {
            setDp(npk, EP);
            pse2 = Wik(npk-1,0) * ~Dp * sp2 + pse;
         }
         for (register int i6=0; i6 < npar; i6++)
               sec(i6,0) = s2a(i6,0);
         for (register int i7=0; i7 < npk-1; i7++)
               sec(npar+i7,0) = pse2(i7,0);
         sec2 = sec2 + sec * (~sec) ;
     }

     //gsl_matrix_free(Ag);
     //gsl_matrix_free(Vg);
     //gsl_matrix_free(DSg);
     //gsl_matrix_free(R1g);
     //gsl_matrix_free(R2g);
     //gsl_vector_free(Sg);
     //gsl_vector_free(work);

     // calculating the g-inverse of sec2
     //gsl_matrix *Ag,*Vg,*DSg,*R1g,*R2g;
     //gsl_vector *Sg,*work;
     // the next 4 are for Ag=U S V^T
     Ag=gsl_matrix_alloc(npar2,npar2);
     Vg=gsl_matrix_alloc(npar2,npar2);
     Sg=gsl_vector_alloc(npar2);
     work=gsl_vector_alloc(npar2);
     // used for inverting the nonzero eigenvalues
     Matrix eigen2(npar2,1);
     //to become the diagonal with the inverse of nonzero eigenvalues
     DSg=gsl_matrix_calloc(npar2,npar2);
     //to become Vg * DSg
     R1g=gsl_matrix_calloc(npar2,npar2);
     //to become Vg * DSg * Vg^T
     R2g=gsl_matrix_calloc(npar2,npar2);
     // to write R2g
     Matrix ginvsec2(npar2,npar2);

     for (int i=0; i<npar2; i++){                  // Ag becomes sec2
         for (int j=0; j<npar2; j++) {
            gsl_matrix_set(Ag,i,j,sec2(i,j));
         }
     }
     gsl_linalg_SV_decomp(Ag,Vg,Sg,work);          // Ag = U S V^T, where Ag becomes U
     for (int i=0; i<npar2; i++)                   // Calculate 1/eigen for eigen > 0
        if (gsl_vector_get(Sg,i) > tol) eigen2(i,0)=1/gsl_vector_get(Sg,i); else eigen2(i,0)=0;
     for (int i=0; i<npar2; i++)                   // Dsg=1/eigen, for eigen > 0
         gsl_matrix_set(DSg,i,i,eigen2(i,0));
     // R1g = Vg * DSg
     gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, Vg, DSg, 0, R1g);
     gsl_matrix_transpose(Ag);
     // R2g = R1g * Ag^T = Vg * DSg * Vg^T thus obtaining the g-inverse
     gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, R1g, Ag, 0, R2g);

     for (int i=0; i<npar2; i++)
        for (int j=0; j<npar2; j++)
           ginvsec2(i,j) = gsl_matrix_get(R2g,i,j);

     infonpml = ginvsec2; // the covariance matrix (not the info matrix as the name might(!) imply)

     gsl_matrix_free(Ag);
     gsl_matrix_free(Vg);
     gsl_matrix_free(DSg);
     gsl_matrix_free(R1g);
     gsl_matrix_free(R2g);
     gsl_vector_free(Sg);
     gsl_vector_free(work);
}

// Log Likelihood of the model

double logL(int *CN, int m, int q, int np, int npk, Matrix &EP, double *model, Matrix &beta, int link, int rslpind, double *rslp,int Ntot, int *resp, double *npoind, int T)
{
    double LL = 0.0;
    double temp;
    for (register int i = 0; i < m; i++)
    {
        temp = 0.0;
        for (int npl = 1; npl <= npk; npl++)
        {
            setfikl(i,CN,q,np,npk,npl,EP,model,beta,link,rslpind,rslp,Ntot,resp,npoind,T);
            temp += fikl*EP(npl-1,0);
        }
        LL += log(temp);
    }
    return LL;
}

// Setting up for the delta method for var(bi): sets up (del G / del beta* )

void deltanp(int q, int np, int npk, Matrix &beta, Matrix &EP, Matrix &deltam, int rslpind)
{
     Matrix templ(q+np+(2+rslpind)*(npk-1),(1+rslpind));
     zerof(templ,q+np+(2+rslpind)*(npk-1),(1+rslpind));
     //sets up sum ri*pi, where ri is the estimated random effect
     Matrix ath((1+rslpind),1);
     zerof(ath,1+rslpind,1);
     for (int c = 0; c < (1+rslpind); c++)
        for (int a = 0; a < npk-1; a++)
           ath(c,0) += beta(q+np+(npk-1)*c+a,0)*EP(a,0);
     // derivatives wrt random effects
     for (int c = 0; c < (1+rslpind); c++)
        for (int a = 0; a < npk-1; a++)
           templ(q+np+(npk-1)*c+a,c) = 2*beta(q+np+(npk-1)*c+a,0)*EP(a,0)+2*ath(c,0)*EP(a,0)/EP(npk-1,0);
     //derivatives wrt probabilities
     for (int c = 0; c < (1+rslpind); c++)
        for (int a = 0; a < npk-1; a++)
           templ(q+np+(1+rslpind)*(npk-1)+a,c) = pow(beta(q+np+(npk-1)*c+a,0),2)+(2*beta(q+np+(npk-1)*c+a,0)*
                                                    ath(c,0)/EP(npk-1,0))+pow(ath(c,0),2)/pow(EP(npk-1,0),2);
     deltam = templ;
}

// Setting up for the delta method for the Kth mass point: sets up (del G(b*) / del beta* )

void deltanpkp(int q, int np, int npk, Matrix &beta, Matrix &EP, Matrix &deltammp, int rslpind)
{
     Matrix templ(q+np+(2+rslpind)*(npk-1),(1+rslpind));
     zerof(templ,q+np+(2+rslpind)*(npk-1),(1+rslpind));
    //sets up sum ri*pi, where ri is the estimated random effect
     Matrix ath((1+rslpind),1);
     zerof(ath,1+rslpind,1);
     for (int c = 0; c < (1+rslpind); c++)
        for (int a = 0; a < npk-1; a++)
           ath(c,0) += beta(q+np+(npk-1)*c+a,0)*EP(a,0);
     // derivatives wrt random effects
     for (int c = 0; c < (1+rslpind); c++)
        for (int a = 0; a < npk-1; a++)
           templ(q+np+(npk-1)*c+a,c) = -EP(a,0)/EP(npk-1,0);
     // derivatives wrt probabilities
     for (int c = 0; c < (1+rslpind); c++)
        for (int a = 0; a < npk-1; a++)
         templ(q+np+(1+rslpind)*(npk-1)+a,c) = -beta(q+np+(npk-1)*c+a,0)/EP(npk-1,0) - ath(c,0)/pow(EP(npk-1,0),2);
     deltammp = templ;
}

/////////////////////////######################################################################################################################

extern "C"{
void npmltd(int *y, int *q1, int *N1, int *m1, int *CN, int *npk1, int *np1, double *model,
            double *eps1,double *strvlint, double *strvlreg, double *strvlmp, double *strvlm,
            double *out, int *EBind1, double *outEB, int *link1, int *obskey1, int *maxit1,
            double *rslp, int *rslpind1, double *outFitted, double *outProb, double *tol1,
            double *npoind, int *T1, double *outparcvmat)
{
int q = q1[0];
int N = N1[0];
int m = m1[0];
int k = q+1;
int npk = npk1[0];
int np = np1[0];
int T = T1[0];
double eps = eps1[0];
int EBind = EBind1[0];
int linkchoice = link1[0];
int obskey = obskey1[0];
int maxit = maxit1[0];
int rslpind = rslpind1[0];
double tol = tol1[0];

// Starting values for probabilities and the beta* vector
Matrix beta(q+np+(npk-1)*(1+rslpind),1);
Matrix EP(npk,1);

for (register int i=0; i<q; i++) beta(i,0) = strvlint[i];
for (register int i=0; i<np; i++) beta(q+i,0) = strvlreg[i];
for (register int i=0; i<(npk-1)*(1+rslpind); i++) beta(q+np+i,0) = strvlmp[i];
for (register int i=0; i<npk; i++) EP(i,0) = strvlm[i];

// number of iterations
Matrix itere(1,1); itere(0,0)=0;

// upbp iterates until the NP model is fit
upbp(npk,beta,EP,m,q,np,CN,model,eps,itere,linkchoice,maxit,rslpind,rslp,tol,N,y,npoind,T);

// For calculating the value of the last mass point: sum_j=1^(npk-1)ai*pi
Matrix ath1(1+rslpind,1);
for (register int a = 0; a < rslpind+1 ; a++)
   ath1(a,0)=0;
for (register int c = 0; c < rslpind+1 ; c++)
   for (register int a = 0; a < npk-1 ; a++)
      ath1(c,0) += beta(q+np+(npk-1)*c+a,0)*EP(a,0);

// Output

// beta
for (register int i=0; i< q+np; i++) out[i] = beta(i,0);

// Mass Points
for (register int c = 0; c < rslpind+1 ; c++){
   for (register int i=0; i< npk; i++){
      if (i<npk-1) out[q+np+npk*c+i] = beta(q+np+(npk-1)*c+i,0);
      if (i==npk-1) out[q+np+npk*(c+1)-1] = -ath1(c,0)/EP(npk-1,0);}}

// Masses
for (register int i=0; i<npk; i++) out[q+np+(rslpind+1)*npk+i] = EP(i,0);

// Random effects covariance
Matrix rec((1+rslpind)*(2+rslpind)/2,1);

for (register int c = 0; c < rslpind+1 ; c++){
   rec(c,0)=EP(npk-1,0)*pow(ath1(c,0)/EP(npk-1,0),2);
   for (register int i=0; i<npk-1 ; i++){
      rec(c,0) += EP(i,0)*pow(beta(q+np+(npk-1)*c+i,0),2);}}

int bc=0;

for (register int c = 0; c < rslpind ; c++){
   for (register int r = (c+1); r < rslpind+1 ; r++){
      rec(rslpind+1+bc,0)=ath1(c,0)*ath1(r,0)/EP(npk-1,0);
      for (register int i=0; i<npk-1 ; i++)
         rec(rslpind+1+bc,0) += EP(i,0)*beta(q+np+(npk-1)*c+i,0)*beta(q+np+(npk-1)*r+i,0);
      bc += 1;}}

for (register int i=0; i<((1+rslpind)*(2+rslpind)/2); i++)
   out[q+np+npk+rslpind*npk+npk+i] = rec(i,0);

// Covariance Matrix based on the suqred derivative of the log-likelihood
Matrix CovNP;
NPMLSE(npk,m,q,np,CN,beta,EP,model,linkchoice,rslpind,rslp,tol,N,y,npoind,T);
CovNP = infonpml;  // the inversion is now done in the function

for (register int i=0; i<(q+np+(rslpind+2)*(npk-1)); i++)
   for (register int j=0; j<(q+np+(rslpind+2)*(npk-1)); j++)
      outparcvmat[i*(q+np+(rslpind+2)*(npk-1))+j]=CovNP(i,j);

//Call delta method for Standard error of k^th mass point
Matrix sekmp(1+rslpind,1+rslpind);
if (npk>1){
   Matrix deltammp(q+np+(2+rslpind)*(npk-1),1+rslpind);
   deltanpkp(q,np,npk,beta,EP,deltammp,rslpind); //updates matrix deltammp for all dimensions
   sekmp = (~deltammp) * CovNP * (deltammp);
}
if (npk==1) zerof(sekmp,1+rslpind,1+rslpind);

// Standard error of the last probability
Matrix selpk(1,1);
if (npk>1){
   Matrix selp(q+np+(2+rslpind)*(npk-1),1);
   zerof(selp,q+np+(2+rslpind)*(npk-1),1);
   for (int i=0; i < npk-1; i++) selp(q+np+(npk-1)*(1+rslpind)+i,0) = 1;
   selpk = (~selp) * CovNP * (selp);
}
if (npk==1) selpk(0,0)=0;

// Standard errors of the q+np coefficients
for (register int i=0; i<q+np; i++)
   out[q+np+npk+rslpind*npk+npk+((1+rslpind)*(2+rslpind)/2)+i]=sqrt(CovNP(i,i));

// Standard errors of Mass points
for (register int c = 0; c < rslpind+1 ; c++){
   for (register int i=0; i< npk; i++){
      if (i<npk-1) out[q+np+npk+rslpind*npk+npk+((1+rslpind)*(2+rslpind)/2)+q+np+npk*c+i] = sqrt(CovNP(q+np+(npk-1)*c+i,q+np+(npk-1)*c+i));
      if (i==npk-1) out[q+np+npk+rslpind*npk+npk+((1+rslpind)*(2+rslpind)/2)+q+np+npk*c+i] = sqrt(sekmp(c,c));}}

// Standard errors of Masses
for (register int i=0; i<npk-1; i++)
   out[q+np+npk+rslpind*npk+npk+((1+rslpind)*(2+rslpind)/2)+q+np+npk*(rslpind+1)+i] = sqrt(CovNP(q+np+(npk-1)*(1+rslpind)+i,q+np+(npk-1)*(1+rslpind)+i));
out[q+np+npk+rslpind*npk+npk+((1+rslpind)*(2+rslpind)/2)+q+np+npk*(rslpind+1)+npk-1]=sqrt(selpk(0,0));

//Call delta method for Standard errors of var(ai), var(bi), ...
Matrix deltam(q+np+(2+rslpind)*(npk-1),1+rslpind);
deltanp(q,np,npk,beta,EP,deltam,rslpind);    //updates matrix deltam for all random effects
Matrix sevar = (~deltam) * CovNP * (deltam);

for (int u=0; u<(1+rslpind); u++)
   out[q+np+npk+rslpind*npk+npk+((1+rslpind)*(2+rslpind)/2)+q+np+npk*(rslpind+1)+npk+u]=sqrt(sevar(u,u));

// log likelihood
out[q+np+npk+rslpind*npk+npk+((1+rslpind)*(2+rslpind)/2)+q+np+npk*(rslpind+1)+npk+rslpind+1]=logL(CN,m,q,np,npk,EP,model,beta,linkchoice,rslpind,rslp,N,y,npoind,T);

//Number of iterations
out[q+np+npk+rslpind*npk+npk+((1+rslpind)*(2+rslpind)/2)+q+np+npk*(rslpind+1)+npk+rslpind+2]=itere(0,0);

//Flag multinomial covariance matrix
out[q+np+npk+rslpind*npk+npk+((1+rslpind)*(2+rslpind)/2)+q+np+npk*(rslpind+1)+npk+rslpind+3]=flagcvm;

//Flag information matrix
out[q+np+npk+rslpind*npk+npk+((1+rslpind)*(2+rslpind)/2)+q+np+npk*(rslpind+1)+npk+rslpind+4]=flaginfo;

// Empirical Bayes: E(bi|Yi)
if ((EBind==1) & (npk>1)){
  double denEB;
  Matrix numEB(1+rslpind,1);
  for (register int i = 0; i < m; i++)
  {
     denEB=0.0;
     for (register int h = 0; h < (1+rslpind); h++)
        numEB(h,0) = 0.0;
     for (int j=0; j < npk; j++)
     {
        setfikl(i,CN,q,np,npk,j+1,EP,model,beta,linkchoice,rslpind,rslp,N,y,npoind,T);
        denEB += fikl*EP(j,0);
        for (register int h = 0; h < (1+rslpind); h++)
           numEB(h,0) += fikl*EP(j,0)*out[q+np+npk*h+j];
     }
    for (register int h = 0; h < (1+rslpind); h++) outEB[i+m*h] = numEB(h,0) / denEB;
  }
}

// Fitted values for the linear predictor at random effects equal to zero if EB=FALSE or at their EB estimate if EB=TRUE

Matrix betanew((q+np+(npk-1)*(1+rslpind)),1);
zerof(betanew,q+np+(npk-1)*(1+rslpind),1);
for (register int i = 0; i < (q+np); i++)
   betanew(i,0)=beta(i,0);

for (register int i = 0; i < m; i++){
   for (register int j = 0; j<(CN[i+1]-CN[i]); j++){
       if (npk>1){
          for (register int h = 0; h<(1+rslpind); h++)
             betanew(q+np+(npk-1)*h,0) = outEB[i+m*h];
       }
       sethijkl(i, j, CN, q, np, npk, 1, EP, model, betanew, rslpind, rslp,N,npoind,T);
       for (register int k=0; k < q; k++)
          outFitted[(CN[i]+j)*q+k]=hijkl(k,0);
    }
}

//Estimated category probabilities at random effects equal to zero if EB=FALSE or at their EB estimate if EB=TRUE

for (register int i = 0; i < m; i++){
   for (register int j = 0; j<(CN[i+1]-CN[i]); j++){
      //clogit link
      if (linkchoice==1){
         outProb[(CN[i]+j)*k] = exp(outFitted[(CN[i]+j)*q])/(1+exp(outFitted[(CN[i]+j)*q]));
         for (register int c=1; c < q; c++)
            outProb[(CN[i]+j)*k+c] = exp(outFitted[(CN[i]+j)*q+c])/(1+exp(outFitted[(CN[i]+j)*q+c]))-exp(outFitted[(CN[i]+j)*q+c-1])/(1+exp(outFitted[(CN[i]+j)*q+c-1]));
         outProb[(CN[i]+j)*k+q] = 1-exp(outFitted[(CN[i]+j)*q+q-1])/(1+exp(outFitted[(CN[i]+j)*q+q-1]));}
      //blogit link: baseline category is the last category
      if (linkchoice==0){
         double partot = 1;
         for (register int c=0; c<q; c++)
            partot += exp(outFitted[(CN[i]+j)*q+c]);
         for (register int c=0; c<q; c++)
            outProb[(CN[i]+j)*k+c]=exp(outFitted[(CN[i]+j)*q+c])/partot;
         outProb[(CN[i]+j)*k+q]=1/partot;}
   }
}

}}
