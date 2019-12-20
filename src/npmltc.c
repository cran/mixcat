#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <R.h>

//Print gsl matrix
void print_matrix(gsl_matrix *A)
{
    int i, j;
    for (i = 0; i < A->size1; i++) {
        for (j = 0; j < A->size2; j++)
            Rprintf("%g\t", gsl_matrix_get(A, i, j));
        Rprintf("\n");
    }
    Rprintf("\n");
}

//Print gsl matrix
void print_vector(gsl_vector *V)
{
    int i;
    for (i = 0; i < V->size; i++) 
        Rprintf("%g\t", V->data[i * V->stride]);
    Rprintf("\n");
}

// Finds the max of the absolute values of a vector s of dim L
double maxfabs(gsl_vector *s, int L)
{
    int i;
    double maxi = -1;
    double g;
    for (i = 0; i < L; i++)
        if ( (g = fabs(gsl_vector_get(s,i)) ) > maxi) maxi = g;
    return maxi;
}

// (i,j)th response vector
void setYij(int i, int j, int q, int *CN, int *resp, gsl_vector *Yij)
{
    gsl_vector_set_zero(Yij);
    if (resp[CN[i]+j] <= q) 
        gsl_vector_set(Yij,resp[CN[i]+j]-1,1); 
}

// Design matrix
void setXijkl(int i, int j, int *CN, int q, int np, int npk, int npl, gsl_vector *EP, double *model, int rslpind, 
              double *rslp, int Ntot, double *npoind, int T, gsl_matrix *Xijkl)
{
    int k, u, l1, jj;
    double temp; 
    gsl_matrix_set_zero(Xijkl);    
    // diagonal qxq for the intercepts
    for (k = 0; k < q; k++)
         gsl_matrix_set(Xijkl,k,k,1);        
    // count the non proportional odds
    int count = 0;
    // npoind is a vector of length equal to the number of predictors (T) that has 1's for non % odds and 0's o/w
    for (u = 0; u < T; u++){
        for (k = 0; k < q; k++)
            gsl_matrix_set(Xijkl,k,q+u+k*npoind[u]+count*(q-1),model[CN[i]+j+u*Ntot]);
        count += npoind[u];
    } 
    // This part accomodates the NP random intercept and slope
    if (npk > 1){
        if (npl < npk){
            for (l1 = 0; l1 < q; l1++){
                gsl_matrix_set(Xijkl,l1,q+np+npl-1,1.0);
                if (rslpind > 0)
                    for (u = 0; u < rslpind; u++)
                        gsl_matrix_set(Xijkl,l1,q+np+(npk-1)*(u+1)+npl-1,rslp[CN[i]+j+u*Ntot]);
            }
        }
        if (npl == npk){           
            for (l1 = 0; l1 < q; l1++){
                for (jj = 0; jj < (npk-1); jj++){                                         
                    temp = - gsl_vector_get(EP,jj) / gsl_vector_get(EP,npk-1);
                    gsl_matrix_set(Xijkl,l1,q+np+jj,temp);
                    if (rslpind > 0)
                        for (u = 0; u < rslpind; u++)                           
                            gsl_matrix_set(Xijkl,l1,q+np+(npk-1)*(u+1)+jj, rslp[CN[i]+j+u*Ntot]*temp);
                }
            }
        }
    }
}

// q-variate linear predictor
void sethijkl(int i, int j, int *CN, int q, int np, int npk, int npl, gsl_vector *EP, double *model, gsl_vector *beta, 
              int rslpind, double *rslp, int Ntot, double *npoind, int T, gsl_matrix *Xijkl, gsl_vector *hijkl)
{
    setXijkl(i, j, CN, q, np, npk, npl, EP, model, rslpind, rslp, Ntot, npoind, T, Xijkl);
    gsl_blas_dgemv(CblasNoTrans,1.0,Xijkl,beta,0.0,hijkl);
}

// q-variate probability vector
void setpijkl(int i, int j, int *CN, int q, int np, int npk, int npl, gsl_vector *EP, double *model, gsl_vector *beta, 
              int link, int rslpind, double *rslp, int Ntot, double *npoind, int T, 
              gsl_matrix *Xijkl, gsl_vector *hijkl, gsl_vector *pijkl)
{
    int k;
    double temp[q], temp2;
    sethijkl(i, j, CN, q, np, npk, npl, EP, model, beta, rslpind, rslp, Ntot, npoind, T, Xijkl, hijkl);
    // clogit link
    if (link == 1){
        for (k = 0; k < q; k++){
            temp2 = gsl_vector_get(hijkl,k);
            temp[k] = exp(temp2)/(1+exp(temp2));
		}
        gsl_vector_set(pijkl,0,temp[0]);
        for (k=1; k< q; k++) 
            gsl_vector_set(pijkl,k,temp[k]-temp[k-1]);       
    }
    // blogit link: baseline category is the last category q+1
    if (link == 0){
        for (k = 0; k < q; k++)
            temp[k] = exp(gsl_vector_get(hijkl,k));
        double partot = 1;
        for (k = 0; k < q; k++) partot += temp[k];
        for (k = 0; k < q; k++) gsl_vector_set(pijkl,k,temp[k]/partot);        
    }
}

// Sums the elements of vector R
double sum1f(gsl_vector *R, int q)
{
    int k;
    double sum2 = 0.0;
    for (k = 0; k < q; k++)
        sum2 += gsl_vector_get(R,k);
    return sum2;
}

// Probability of observing what we observed for cluster i at time j
double setfijkl(int i, int j, int *CN, int q, int np, int npk, int npl, gsl_vector *EP, double *model, gsl_vector *beta, 
                int link, int rslpind, double *rslp, int Ntot, int *resp, double *npoind, int T, 
                gsl_matrix *Xijkl, gsl_vector *hijkl, gsl_vector *pijkl, gsl_vector *Yij)
{
    int t;
    double fijkl = 1.0;
    setpijkl(i, j, CN, q, np, npk, npl, EP, model, beta, link, rslpind, rslp, Ntot, npoind, T, Xijkl, hijkl, pijkl);
    setYij(i, j, q, CN, resp, Yij);
    for (t = 0; t < q; t++)
        fijkl *= pow(gsl_vector_get(pijkl,t),gsl_vector_get(Yij,t));
    fijkl *= pow(1-sum1f(pijkl,q), 1-sum1f(Yij,q));
    return fijkl;
}

// Probability of observing what we observed for cluster i
double setfikl(int i, int *CN, int q, int np, int npk, int npl, gsl_vector *EP, double *model, gsl_vector *beta, 
               int link, int rslpind, double *rslp, int Ntot, int *resp, double *npoind, int T,
               gsl_matrix *Xijkl, gsl_vector *hijkl, gsl_vector *pijkl, gsl_vector *Yij)
{
    int j;
    double fijkl;
    double fikl = 1.0;    
    for (j = 0; j < (CN[i+1]-CN[i]); j++){
        fijkl = setfijkl(i,j,CN,q,np,npk,npl,EP,model,beta,link,rslpind,rslp,Ntot,resp,npoind,T,Xijkl,hijkl,pijkl,Yij);
        fikl *= fijkl;
    }
    return fikl;
}

// Sij = covariance matrix based on pij
void setSijkl(int i, int j, int *CN, int q, int np, int npk, int npl, gsl_vector *EP, double *model, 
              gsl_vector *beta, int link, int rslpind, double *rslp, int Ntot, double *npoind, int T,
              gsl_matrix *Xijkl, gsl_vector *hijkl, gsl_vector *pijkl, gsl_matrix *Sijkl)
{
    int k;
    setpijkl(i,j,CN,q,np,npk,npl,EP,model,beta,link,rslpind,rslp,Ntot,npoind,T,Xijkl,hijkl,pijkl);        
    gsl_matrix_set_zero(Sijkl);           
    for (k = 0; k < q; k++)
        gsl_matrix_set(Sijkl,k,k,gsl_vector_get(pijkl,k)); 
    gsl_blas_dger(-1.0, pijkl, pijkl, Sijkl);
}

// Dij = der h(h) / der h, where h is the inverse link function
void setDijkl(int i, int j, int *CN, int q, int np, int npk, int npl, gsl_vector *EP, double *model, 
              gsl_vector *beta, int link, int rslpind, double *rslp, int Ntot, double *npoind, int T,
              gsl_matrix *Xijkl, gsl_vector *hijkl, gsl_vector *pijkl, gsl_matrix *Dijkl)
{
    int k, l;
    double temp, temp2;
    // clogit link
    if (link == 1){		
        sethijkl(i,j,CN,q,np,npk,npl,EP,model,beta,rslpind,rslp,Ntot,npoind,T,Xijkl,hijkl);
        gsl_matrix_set_zero(Dijkl);
        for (k = 0; k < q; k++){
            temp = gsl_vector_get(hijkl,k);
            gsl_matrix_set(Dijkl,k,k,exp(temp)/pow((1+exp(temp)),2));
            if (k < (q-1)) gsl_matrix_set(Dijkl,k+1,k,-exp(temp)/pow((1+exp(temp)),2));//this line replaces the loop below 
		}
        //for (k = 0; k < (q-1); k++){
        //    temp = gsl_vector_get(hijkl,k);
        //    gsl_matrix_set(Dijkl,k+1,k,-exp(temp)/pow((1+exp(temp)),2));
		//}
    }
    // blogit link
    if (link == 0){		               
        setpijkl(i, j, CN, q, np, npk, npl, EP, model, beta, link,rslpind,rslp,Ntot,npoind,T,Xijkl,hijkl,pijkl);
        for (k = 0; k < q; k++){
            for (l = 0; l < q; l++){
                temp = gsl_vector_get(pijkl,k);
                temp2 = gsl_vector_get(pijkl,l);
                if (k == l) {gsl_matrix_set(Dijkl,k,l,temp*(1-temp));}else{gsl_matrix_set(Dijkl,k,l,-temp*temp2);}
			}
		}
    }
}

// Sij and Dij together
void setSDijkl(int i, int j, int *CN, int q, int np, int npk, int npl, gsl_vector *EP, double *model, 
               gsl_vector *beta, int link, int rslpind, double *rslp, int Ntot, double *npoind, int T,
               gsl_matrix *Xijkl, gsl_vector *hijkl, gsl_vector *pijkl, gsl_matrix *Sijkl, gsl_matrix *Dijkl)
{
    int k, l;
    double temp[q], temp2, temp3;
    sethijkl(i, j, CN, q, np, npk, npl, EP, model, beta, rslpind, rslp, Ntot, npoind, T, Xijkl, hijkl);
    ///////////////////////////Probabilities
    // clogit link
    if (link == 1){
        for (k = 0; k < q; k++){
            temp2 = gsl_vector_get(hijkl,k);
            temp[k] = exp(temp2)/(1+exp(temp2));
		}
        gsl_vector_set(pijkl,0,temp[0]);
        for (k=1; k< q; k++) 
            gsl_vector_set(pijkl,k,temp[k]-temp[k-1]);       
    }
    // blogit link: baseline category is the last category q+1
    if (link == 0){
        for (k = 0; k < q; k++)
            temp[k] = exp(gsl_vector_get(hijkl,k));
        double partot = 1;
        for (k = 0; k < q; k++) partot += temp[k];
        for (k = 0; k < q; k++) gsl_vector_set(pijkl,k,temp[k]/partot);        
    }
    /////////////////////////// Sij    
    gsl_matrix_set_zero(Sijkl);           
    for (k = 0; k < q; k++)
        gsl_matrix_set(Sijkl,k,k,gsl_vector_get(pijkl,k)); 
    gsl_blas_dger(-1.0, pijkl, pijkl, Sijkl);        
    /////////////////////////// Dij
    // clogit link
    if (link == 1){		
        sethijkl(i,j,CN,q,np,npk,npl,EP,model,beta,rslpind,rslp,Ntot,npoind,T,Xijkl,hijkl);
        gsl_matrix_set_zero(Dijkl);
        for (k = 0; k < q; k++){
            temp3 = gsl_vector_get(hijkl,k);
            gsl_matrix_set(Dijkl,k,k,exp(temp3)/pow((1+exp(temp3)),2));
            if (k < (q-1)) gsl_matrix_set(Dijkl,k+1,k,-exp(temp3)/pow((1+exp(temp3)),2));//this line replaces the loop below 
		}
        //for (k = 0; k < (q-1); k++){
        //    temp3 = gsl_vector_get(hijkl,k);
        //    gsl_matrix_set(Dijkl,k+1,k,-exp(temp3)/pow((1+exp(temp3)),2));
		//}
    }
    // blogit link
    if (link == 0){		               
        setpijkl(i, j, CN, q, np, npk, npl, EP, model, beta, link,rslpind,rslp,Ntot,npoind,T,Xijkl,hijkl,pijkl);
        for (k = 0; k < q; k++){
            for (l = 0; l < q; l++){
                temp3 = gsl_vector_get(pijkl,k);
                temp2 = gsl_vector_get(pijkl,l);
                if (k == l) {gsl_matrix_set(Dijkl,k,l,temp3*(1-temp3));}else{gsl_matrix_set(Dijkl,k,l,-temp3*temp2);}
			}
		}
    }
}

// Posterior probabilities of group membership
void setWik(int i, int *CN, int q, int np, int npk, gsl_vector *EP, double *model, gsl_vector *beta, 
            int link, int rslpind, double *rslp, int Ntot, int *resp, double *npoind, int T,
            gsl_matrix *Xijkl, gsl_vector *hijkl, gsl_vector *pijkl, gsl_vector *Yij, gsl_vector *Wik)
{
    int l;
    double fikl;
    double probs[npk];
    double tot = 0.0;
    for (l = 0; l < npk; l++){
        fikl = setfikl(i,CN,q,np,npk,l+1,EP,model,beta,link,rslpind,rslp,Ntot,resp,npoind,T,Xijkl,hijkl,pijkl,Yij);
        probs[l] = fikl * gsl_vector_get(EP,l);
        tot += probs[l];
    }
    for (l = 0; l < npk; l++)
        gsl_vector_set(Wik,l,probs[l]/tot);
}

//###############################Score and Fisher information matrix

//2 functions for calculating the score and info of the probabilities

void setZk(int i, int j, int *CN, int npk, int q, int np, gsl_vector *beta, int rslpind, double *rslp, int Ntot, gsl_matrix *Zk)
{
    int c, r, u;
    double temp;
    for (c = 0; c < npk-1; c++){
        for (r = 0; r < q; r++){
            temp = gsl_vector_get(beta,q+np+c);
            if (rslpind>0){				
                for (u = 0; u < rslpind; u++)
                    temp += rslp[(CN[i]+j+u*Ntot)]*gsl_vector_get(beta,q+np+(npk-1)*(u+1)+c);
            }
            gsl_matrix_set(Zk,r,c,temp);
        }
    }
}

void setDp(int npk, gsl_vector *EP, gsl_matrix *Dp)
{
    int c, r;
    for (r = 0; r < npk-1; r++)
        for (c = 0; c < npk-1; c++)
            gsl_matrix_set(Dp,r,c,-gsl_vector_get(EP,r)/pow(gsl_vector_get(EP,npk-1),2));
    for (r = 0; r < npk-1; r++)
        gsl_matrix_set(Dp,r,r,(sum1f(EP,npk-1)-1-gsl_vector_get(EP,r))/pow(gsl_vector_get(EP,npk-1),2));
}

//Generalized inverse (with flag)
void ginv(int p, double tol, gsl_matrix *A, int* flag){
    int i; 
    double temp, max;
    *flag = 0;
    gsl_matrix *D = gsl_matrix_calloc(p,p);
    gsl_matrix *M = gsl_matrix_alloc(p,p);
    gsl_matrix *N = gsl_matrix_alloc(p,p);
    gsl_vector *eval = gsl_vector_alloc(p);
    gsl_matrix *evec = gsl_matrix_alloc(p,p);
    gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(p);
    gsl_eigen_symmv(A,eval,evec,w);
    gsl_eigen_symmv_sort(eval,evec,GSL_EIGEN_SORT_VAL_DESC);
    max = 1; //gsl_vector_get(eval,0);
    for (i=0; i < p; i++){ 
        temp = gsl_vector_get(eval,i);        
	    if (temp > (tol * max)) gsl_matrix_set(D,i,i,1/temp);else{gsl_matrix_set(D,i,i,0.0); *flag = 1;}
    }
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,evec,D,0.0,M);  //D1 = V D
    gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,M,evec,0.0,N); //D2 = D1 A = ginv(A)
    gsl_matrix_memcpy(A,N);
    gsl_matrix_free(D); gsl_matrix_free(M); gsl_matrix_free(N);
    gsl_vector_free(eval); gsl_matrix_free(evec); gsl_eigen_symmv_free(w);
}

//Updates the score vector and Fisher information matrix
void NPML(int npk, int m, int q, int np, int *CN, gsl_vector *beta, gsl_vector *EP, double *model, 
          int link, int rslpind, double *rslp, double tol, int Ntot, int *resp, double *npoind, int T, int* flag,
          gsl_matrix *Xijkl, gsl_vector *hijkl, gsl_vector *pijkl, gsl_vector *Yij, gsl_vector *Wik,
          gsl_vector *SNP, gsl_matrix *FNP, gsl_vector *spnp, gsl_matrix *fpnp, gsl_matrix *fmix){
    int i, j, npl;
    int ncol = q+np+(npk-1)*(1+rslpind); // dim of beta
    
    gsl_matrix *VT = gsl_matrix_alloc(ncol,q); 
    gsl_vector *s1a = gsl_vector_alloc(ncol); 
    gsl_vector *s2a = gsl_vector_alloc(ncol);  
    gsl_matrix *s1b = gsl_matrix_alloc(ncol,ncol); 
    gsl_matrix *s2b = gsl_matrix_alloc(ncol,ncol); 
    gsl_matrix *VT2 = gsl_matrix_alloc(q,q); 
    gsl_vector *sp2 = gsl_vector_alloc(npk-1); 
    gsl_vector *sp22 = gsl_vector_calloc(npk-1); 
    gsl_matrix *fp2 = gsl_matrix_alloc(npk-1,npk-1); 
    gsl_matrix *fp22 = gsl_matrix_calloc(npk-1,npk-1); 
    gsl_matrix *fp222 = gsl_matrix_calloc(npk-1,npk-1); 
    gsl_matrix *VT3 = gsl_matrix_alloc(q,ncol); 
    gsl_matrix *fmix1 = gsl_matrix_alloc(npk-1,ncol); 
    gsl_matrix *fmix2 = gsl_matrix_calloc(npk-1,ncol); 
    gsl_matrix *ginvSijkl = gsl_matrix_alloc(q,q);
    gsl_matrix *Dijkl = gsl_matrix_alloc(q,q); 
    gsl_matrix *Dp = gsl_matrix_alloc(npk-1,npk-1);
    gsl_matrix *Zk = gsl_matrix_alloc(q,npk-1);                                              
    gsl_matrix *ZkVT2 = gsl_matrix_alloc(npk-1,q);
    gsl_matrix *ZVD = gsl_matrix_alloc(npk-1,q);
    gsl_vector *WikCopy = gsl_vector_alloc(npk); 
    gsl_matrix *mat = gsl_matrix_calloc(npk-1,npk-1);
    gsl_matrix *Dpfp22 = gsl_matrix_alloc(npk-1,npk-1);
    gsl_matrix *DfD = gsl_matrix_alloc(npk-1,npk-1);
    gsl_vector_view X, DiagM;
    
    gsl_vector_set_zero(SNP);      
    gsl_matrix_set_zero(FNP);
    gsl_vector_set_zero(spnp);
    gsl_matrix_set_zero(fpnp);          
    
    for (i = 0; i < m; i++){
        gsl_vector_set_zero(s2a); 
        gsl_matrix_set_zero(s2b); 
        gsl_vector_set_zero(sp2); 
        gsl_matrix_set_zero(fp2); 
        gsl_matrix_set_zero(fmix1);                 
        setWik(i, CN, q, np, npk, EP, model, beta, link, rslpind, rslp, Ntot, resp, npoind, T, Xijkl, hijkl, pijkl, Yij, Wik);
        for (npl = 1; npl <= npk; npl++){
            gsl_vector_set_zero(s1a); 
            gsl_matrix_set_zero(s1b); 
            for (j = 0; j < (CN[i+1] - CN[i]); j++){
                //setXijkl(i, j, CN, q, np, npk, npl, EP, model, rslpind, rslp,Ntot);
                setSDijkl(i, j, CN, q, np, npk, npl, EP, model, beta, link, rslpind, rslp,Ntot,npoind,T,
                          Xijkl, hijkl, pijkl, ginvSijkl, Dijkl);                
                setYij(i, j, q, CN, resp, Yij);                               
                ginv(q,tol,ginvSijkl,flag);                
                gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,Dijkl,ginvSijkl,0.0,VT2);
                gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,Xijkl,VT2,0.0,VT);
                gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,Dijkl,Xijkl,0.0,VT3);                  
                gsl_vector_sub(Yij,pijkl);
                gsl_blas_dgemv(CblasNoTrans,1.0,VT,Yij,1.0,s1a);
                gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,VT,VT3,1.0,s1b);                                 
                if ((npk == npl) & (npk > 1)){
                    if (rslpind > 0)                    setZk(i, j, CN, npk, q, np, beta, rslpind, rslp, Ntot, Zk);
                    if ((rslpind==0) & (i==0) & (j==0)) setZk(i, j, CN, npk, q, np, beta, rslpind, rslp, Ntot, Zk);                                            
                    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,Zk,VT2,0.0,ZkVT2);                                                  
                    gsl_blas_dgemv(CblasNoTrans,1.0,ZkVT2,Yij,1.0,sp2);                  
                    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,ZkVT2,Dijkl,0.0,ZVD);
                    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,ZVD,Zk,1.0,fp2);                                                            
                    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,ZkVT2,VT3,1.0,fmix1);
                }
            }
            gsl_vector_scale(s1a,gsl_vector_get(Wik,npl-1));
            gsl_vector_add(s2a,s1a);
            gsl_matrix_scale(s1b,gsl_vector_get(Wik,npl-1));
            gsl_matrix_add(s2b,s1b);
        }
        gsl_vector_add(SNP,s2a);
        gsl_matrix_add(FNP,s2b);
        if (npk > 1){			
			gsl_vector_scale(sp2,gsl_vector_get(Wik,npk-1));
            gsl_vector_add(sp22,sp2);
            gsl_matrix_scale(fp2,gsl_vector_get(Wik,npk-1));
            gsl_matrix_add(fp22,fp2);                                    
            gsl_vector_memcpy(WikCopy,Wik);    
            gsl_vector_div(WikCopy,EP);
            gsl_vector_add_constant(WikCopy,-gsl_vector_get(Wik,npk-1)/gsl_vector_get(EP,npk-1));    
            X = gsl_vector_subvector(WikCopy,0,npk-1);                        
            gsl_vector_add(spnp,&X.vector);                        
            gsl_matrix_set_all(fp222,gsl_vector_get(Wik,npk-1)/pow(gsl_vector_get(EP,npk-1),2));            
            gsl_vector_memcpy(WikCopy,Wik);                
            gsl_vector_div(WikCopy,EP);            
            gsl_vector_div(WikCopy,EP);            
            gsl_matrix_set_zero(mat);
            DiagM = gsl_matrix_diagonal(mat);                                     
            X = gsl_vector_subvector(WikCopy,0,npk-1);
            gsl_vector_memcpy(&DiagM.vector,&X.vector);
            gsl_matrix_add(fp222,mat);                                    
            gsl_matrix_add(fpnp,fp222);            
            gsl_matrix_scale(fmix1,gsl_vector_get(Wik,npk-1));            
            gsl_matrix_add(fmix2,fmix1);               
        }
    }    
    if (npk > 1){
        setDp(npk, EP, Dp);        
        gsl_blas_dgemv(CblasTrans,1.0,Dp,sp22,1.0,spnp);        
        gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,Dp,fp22,0.0,Dpfp22);
        gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,Dpfp22,Dp,0.0,DfD);
        gsl_matrix_add(fpnp,DfD);
        gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,Dp,fmix2,0.0,fmix);
    }
    gsl_matrix_free(VT); gsl_vector_free(s1a); gsl_vector_free(s2a); gsl_matrix_free(s1b); gsl_matrix_free(s2b); 
    gsl_matrix_free(VT2); gsl_vector_free(sp2); gsl_vector_free(sp22); gsl_matrix_free(fp2); gsl_matrix_free(fp22); 
    gsl_matrix_free(fp222); gsl_matrix_free(VT3); gsl_matrix_free(fmix1); gsl_matrix_free(fmix2); gsl_matrix_free(ginvSijkl); 
    gsl_matrix_free(Dijkl); gsl_matrix_free(Dp); gsl_matrix_free(Zk); gsl_matrix_free(ZkVT2); gsl_matrix_free(ZVD);
    gsl_vector_free(WikCopy); gsl_matrix_free(mat); gsl_matrix_free(Dpfp22); gsl_matrix_free(DfD);
}

// A function that iterates until the NP model is fit
void upbp(int npk, gsl_vector *beta, gsl_vector *EP, int m, int q, int np, int *CN, double *model, double eps, int* iter, 
          int link, int maxit, int rslpind, double *rslp, double tol, int Ntot, int *resp, double *npoind, int T,
          int* flagcvm, int* flaginfo, gsl_matrix *Xijkl, gsl_vector *hijkl, gsl_vector *pijkl, gsl_vector *Yij, gsl_vector *Wik,
          gsl_vector *SNP, gsl_matrix *FNP, gsl_vector *spnp, gsl_matrix *fpnp, gsl_matrix *fmix, gsl_matrix *infonpml)
{
    int flag;
    *flagcvm = 0;
    *flaginfo = 0;
    int nfpar = q+np+(2+rslpind)*(npk-1);
    int ncol = q+np+(npk-1)*(1+rslpind);
    gsl_vector *ua = gsl_vector_alloc(nfpar); 
    gsl_vector *scrll = gsl_vector_alloc(nfpar); 
    gsl_matrix *FB = gsl_matrix_alloc(nfpar,nfpar);           
    gsl_matrix *tfmix = gsl_matrix_calloc(ncol,npk-1);
    gsl_vector_view X, EPsub;
    gsl_matrix_view M;
               
    gsl_vector_set_all(scrll,1);        
    //gsl_vector_set(scrll,1,1);
    *iter = 0;
    int i;
    while ((maxfabs(scrll,nfpar) > eps) && (*iter < maxit)){         
        NPML(npk,m,q,np,CN,beta,EP,model,link,rslpind,rslp,tol,Ntot,resp,npoind,T,
             &flag,Xijkl,hijkl,pijkl,Yij,Wik,SNP,FNP,spnp,fpnp,fmix);
        *iter += 1;
        if (flag==1) *flagcvm = *iter;        
        X = gsl_vector_subvector(scrll,0,q+np+(npk-1)*(1+rslpind));
        gsl_vector_memcpy(&X.vector,SNP);            
        X = gsl_vector_subvector(scrll,q+np+(npk-1)*(1+rslpind),npk-1);
        gsl_vector_memcpy(&X.vector,spnp);                    
        M = gsl_matrix_submatrix(FB,0,0,q+np+(npk-1)*(1+rslpind),q+np+(npk-1)*(1+rslpind));
        gsl_matrix_memcpy(&M.matrix,FNP);                                   
        M = gsl_matrix_submatrix(FB,q+np+(npk-1)*(1+rslpind),q+np+(npk-1)*(1+rslpind),npk-1,npk-1);
        gsl_matrix_memcpy(&M.matrix,fpnp);        
        M = gsl_matrix_submatrix(FB,q+np+(npk-1)*(1+rslpind),0,npk-1,q+np+(npk-1)*(1+rslpind));
        gsl_matrix_memcpy(&M.matrix,fmix);
        M = gsl_matrix_submatrix(FB,0,q+np+(npk-1)*(1+rslpind),q+np+(npk-1)*(1+rslpind),npk-1);
        gsl_matrix_transpose_memcpy(tfmix,fmix);
        gsl_matrix_memcpy(&M.matrix,tfmix);                
        gsl_matrix_memcpy(infonpml,FB);
		ginv(nfpar,tol,infonpml,&flag);
        if (flag==1) *flaginfo = *iter;          
        gsl_blas_dgemv(CblasTrans,1.0,infonpml,scrll,0.0,ua);
        X = gsl_vector_subvector(ua,0,q+np+(npk-1)*(1+rslpind));        
        gsl_vector_add(beta,&X.vector);
        X = gsl_vector_subvector(ua,q+np+(npk-1)*(1+rslpind),npk-1);        
        EPsub = gsl_vector_subvector(EP,0,npk-1);
        gsl_vector_add(&EPsub.vector,&X.vector);
        gsl_vector_set(EP,npk-1,1-sum1f(EP,npk-1));
    }
    gsl_vector_free(ua); gsl_vector_free(scrll); gsl_matrix_free(FB); gsl_matrix_free(tfmix);
}
    
// Finding the SE for the NPML -obs mat
void NPMLSE(int npk, int m, int q, int np, int *CN, gsl_vector *beta, gsl_vector *EP, double *model, int link, 
            int rslpind, double *rslp, double tol, int Ntot, int *resp, double *npoind, int T, 
            gsl_matrix *Xijkl, gsl_vector *hijkl, gsl_vector *pijkl, gsl_vector *Yij, gsl_vector *Wik, gsl_matrix *infonpml)
{
    int i, npl, j, flag;
    int npar = q+np+(npk-1)*(1+rslpind);                    // fixed effects and mass points
    int npar2 = q+np+(2+rslpind)*(npk-1);                   // all parameters    
    double temp;
    gsl_vector *s1a = gsl_vector_alloc(npar);
    gsl_vector *s2a = gsl_vector_alloc(npar);
    gsl_vector *pse = gsl_vector_alloc(npk-1);
    gsl_matrix *VT3 = gsl_matrix_alloc(q,q); 
    gsl_vector *sp2 = gsl_vector_alloc(npk-1);
    gsl_vector *sec = gsl_vector_alloc(npar2);
    gsl_matrix *sec2 = gsl_matrix_calloc(npar2,npar2);
    gsl_matrix *ginvSijkl = gsl_matrix_alloc(q,q);
    gsl_matrix *Dijkl = gsl_matrix_alloc(q,q);
    gsl_matrix *VT = gsl_matrix_alloc(npar,q); 
    gsl_matrix *Zk = gsl_matrix_alloc(q,npk-1);
    gsl_matrix *ZkVT3 = gsl_matrix_alloc(npk-1,q);
    gsl_matrix *Dp = gsl_matrix_alloc(npk-1,npk-1);
    gsl_vector *Dpsp2 = gsl_vector_alloc(npk-1);
    gsl_vector_view X;
    
    for (i = 0; i < m; i++){
        gsl_vector_set_zero(pse);//is this needed?
        gsl_vector_set_zero(s2a);
        gsl_vector_set_zero(sp2);                
        setWik(i, CN, q, np, npk, EP, model, beta, link, rslpind, rslp, Ntot, resp, npoind, T, Xijkl, hijkl, pijkl, Yij, Wik);                         
        for (npl = 1; npl <= npk; npl++){			
            if (npl < npk){ 
                temp = gsl_vector_get(Wik,npl-1)/gsl_vector_get(EP,npl-1) - gsl_vector_get(Wik,npk-1)/gsl_vector_get(EP,npk-1);
                gsl_vector_set(pse,npl-1,temp);  
			}			
            gsl_vector_set_zero(s1a);
            for (j = 0; j< (CN[i+1]-CN[i]); j++){
                //setXijkl(i, j, CN, m, q, np, npk, npl, EP, model);
                setSDijkl(i, j, CN, q, np, npk, npl, EP, model, beta, link, rslpind, rslp,Ntot,npoind,T,
                          Xijkl, hijkl, pijkl, ginvSijkl, Dijkl);
                setYij(i, j, q, CN, resp, Yij);                    
                ginv(q,tol,ginvSijkl,&flag);                                                        
                gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,Dijkl,ginvSijkl,0.0,VT3);                
                gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,Xijkl,VT3,0.0,VT);                                 
                gsl_vector_sub(Yij,pijkl);
                gsl_blas_dgemv(CblasNoTrans,1.0,VT,Yij,1.0,s1a);
                if ((npk == npl) & (npk > 1)){
                    if (rslpind>0)                      setZk(i, j, CN, npk, q, np, beta, rslpind, rslp, Ntot, Zk);
                    if ((rslpind==0) & (i==0) & (j==0)) setZk(i, j, CN, npk, q, np, beta, rslpind, rslp, Ntot, Zk);                    
                    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,Zk,VT3,0.0,ZkVT3);                                                                      
                    gsl_blas_dgemv(CblasNoTrans,1.0,ZkVT3,Yij,1.0,sp2);                                      
                }
            }
            gsl_vector_scale(s1a,gsl_vector_get(Wik,npl-1));
            gsl_vector_add(s2a,s1a);
        }         
        if (npk > 1){
			setDp(npk, EP, Dp);                    
            gsl_blas_dgemv(CblasTrans,1.0,Dp,sp2,0.0,Dpsp2);        
            gsl_vector_scale(Dpsp2,gsl_vector_get(Wik,npk-1));
            gsl_vector_add(pse,Dpsp2);
        }        
        X = gsl_vector_subvector(sec,0,npar);
        gsl_vector_memcpy(&X.vector,s2a);               
        X = gsl_vector_subvector(sec,npar,npk-1);
        gsl_vector_memcpy(&X.vector,pse);        
        gsl_blas_dger(1.0, sec, sec, sec2);
              
    }
    
    gsl_matrix_memcpy(infonpml,sec2);
    ginv(npar2,tol,infonpml,&flag); // the covariance matrix 
    gsl_vector_free(s1a); gsl_vector_free(s2a); gsl_vector_free(pse); gsl_matrix_free(VT3); gsl_vector_free(sp2);
    gsl_vector_free(sec); gsl_matrix_free(sec2); gsl_matrix_free(ginvSijkl); gsl_matrix_free(Dijkl); 
    gsl_matrix_free(VT); gsl_matrix_free(Zk); gsl_matrix_free(ZkVT3); gsl_matrix_free(Dp); gsl_vector_free(Dpsp2);
}

// Log Likelihood of the model
double logL(int *CN, int m, int q, int np, int npk, gsl_vector *EP, double *model, gsl_vector *beta, int link, 
            int rslpind, double *rslp,int Ntot, int *resp, double *npoind, int T,
            gsl_matrix *Xijkl, gsl_vector *hijkl, gsl_vector *pijkl, gsl_vector *Yij)
{
    int i, npl;
    double LL = 0.0;
    double fikl;
    double temp;
    for (i = 0; i < m; i++){
        temp = 0.0;
        for (npl = 1; npl <= npk; npl++){
            fikl = setfikl(i,CN,q,np,npk,npl,EP,model,beta,link,rslpind,rslp,Ntot,resp,npoind,T,Xijkl,hijkl,pijkl,Yij);
            temp += fikl * gsl_vector_get(EP,npl-1);
        }
        LL += log(temp);
    }
    return LL;
}

// Setting up for the delta method for var(bi): sets up (del G / del beta* )
void deltanp(int q, int np, int npk, gsl_vector *beta, gsl_vector *EP, gsl_matrix *deltam, int rslpind)
{
    int a, c;    
    double temp, temp2;
    //sets up sum ri*pi, where ri is the estimated random effect
    double ath[1+rslpind];
    for (c = 0; c < (1+rslpind); c++){
        temp = 0;
        for (a = 0; a < npk-1; a++)
            temp += gsl_vector_get(beta,q+np+(npk-1)*c+a) * gsl_vector_get(EP,a); 
        ath[c] = temp; 
	}
    // derivatives wrt random effects
    temp = gsl_vector_get(EP,npk-1);
    for (c = 0; c < (1+rslpind); c++){
        for (a = 0; a < npk-1; a++){  
            temp2 = 2*gsl_vector_get(EP,a)*(gsl_vector_get(beta,q+np+(npk-1)*c+a) + ath[c]/temp);    
            gsl_matrix_set(deltam,q+np+(npk-1)*c+a,c,temp2);
		}
	}
    //derivatives wrt probabilities
    for (c = 0; c < (1+rslpind); c++){
        for (a = 0; a < npk-1; a++){
			temp2 = pow(gsl_vector_get(beta,q+np+(npk-1)*c+a),2) + 2*gsl_vector_get(beta,q+np+(npk-1)*c+a)*ath[c]/temp +
			            pow(ath[c],2)/pow(temp,2);
            gsl_matrix_set(deltam,q+np+(1+rslpind)*(npk-1)+a,c,temp2); 
	    }
    }
}

// Setting up for the delta method for the Kth mass point: sets up (del G(b*) / del beta* )
void deltanpkp(int q, int np, int npk, gsl_vector *beta, gsl_vector *EP, gsl_matrix *deltammp, int rslpind)
{
    int a, c;    
    double temp, temp2;           
    //sets up sum ri*pi, where ri is the estimated random effect  
    double ath[1+rslpind]; 
    for (c = 0; c < (1+rslpind); c++){
        temp = 0;
        for (a = 0; a < npk-1; a++)
            temp += gsl_vector_get(beta,q+np+(npk-1)*c+a) * gsl_vector_get(EP,a);
        ath[c] = temp;
    }  
    // derivatives wrt random effects
    temp = gsl_vector_get(EP,npk-1);
    for (c = 0; c < (1+rslpind); c++)        
        for (a = 0; a < npk-1; a++)
            gsl_matrix_set(deltammp,q+np+(npk-1)*c+a,c,-gsl_vector_get(EP,a)/temp);
    // derivatives wrt probabilities
    for (c = 0; c < (1+rslpind); c++){
        for (a = 0; a < npk-1; a++){
            temp2 = - gsl_vector_get(beta,q+np+(npk-1)*c+a) / temp - ath[c] / pow(temp,2);
            gsl_matrix_set(deltammp,q+np+(1+rslpind)*(npk-1)+a,c,temp2);
		}
	}
}

/////////////////////////###################################################################

void npmltd(int *y, int *q1, int *N1, int *m1, int *CN, int *npk1, int *np1, double *model,
            double *eps1, double *strvlint, double *strvlreg, double *strvlmp, double *strvlm,
            double *out, int *EBind1, double *outEB, int *link1, int *maxit1,
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
    int maxit = maxit1[0];
    int rslpind = rslpind1[0];
    double tol = tol1[0];

    int i, a, c, r, j, h, u; 
    int flagcvm, flaginfo, iter;
    int nfpar = q+np+(2+rslpind)*(npk-1);
    int ncol =  q+np+(1+rslpind)*(npk-1);
    double fikl;
    
    gsl_vector *beta = gsl_vector_alloc(ncol);     
    gsl_vector *EP = gsl_vector_alloc(npk);     
    gsl_vector *Yij = gsl_vector_alloc(q); 
    gsl_vector *hijkl = gsl_vector_alloc(q);
    gsl_vector *pijkl = gsl_vector_alloc(q);
    gsl_vector *Wik = gsl_vector_alloc(npk);
    gsl_vector *SNP = gsl_vector_alloc(ncol); //score of the fixed effects and the NP mass points
    gsl_vector *spnp = gsl_vector_calloc(npk-1); //score of the probability vector ... needed if npk > 1 
    gsl_vector *selp = gsl_vector_calloc(q+np+(2+rslpind)*(npk-1));
    gsl_vector *is = gsl_vector_alloc(q+np+(2+rslpind)*(npk-1));
    gsl_vector *betanew = gsl_vector_calloc(q+np+(npk-1)*(1+rslpind));    
    gsl_matrix *Xijkl = gsl_matrix_alloc(q,ncol);
    gsl_matrix *Sijkl = gsl_matrix_alloc(q,q);
    gsl_matrix *Dijkl = gsl_matrix_alloc(q,q);
    gsl_matrix *FNP = gsl_matrix_alloc(ncol,ncol); //information matrix for the fixed effects and the NP mass points 
    gsl_matrix *fpnp = gsl_matrix_calloc(npk-1,npk-1); //information for the probabilities
    gsl_matrix *fmix = gsl_matrix_calloc(npk-1,ncol); //mixed second derivatives
    gsl_matrix *infonpml = gsl_matrix_alloc(nfpar,nfpar);
    gsl_matrix *deltammp = gsl_matrix_calloc(q+np+(2+rslpind)*(npk-1),(1+rslpind));
    gsl_matrix *deltam = gsl_matrix_calloc(q+np+(2+rslpind)*(npk-1),(1+rslpind));
    gsl_matrix *sevar = gsl_matrix_calloc(1+rslpind,1+rslpind);
    gsl_matrix *sekmp = gsl_matrix_calloc(1+rslpind,1+rslpind);
    gsl_matrix *DI = gsl_matrix_alloc(1+rslpind,q+np+(2+rslpind)*(npk-1));

    // Starting values for beta and probabilities        
    for (i = 0; i < q; i++) gsl_vector_set(beta,i,strvlint[i]);
    for (i = 0; i < np; i++) gsl_vector_set(beta,q+i,strvlreg[i]);
    for (i = 0; i < ((npk-1)*(1+rslpind)); i++) gsl_vector_set(beta,q+np+i,strvlmp[i]); 
    for (i = 0; i < npk; i++) gsl_vector_set(EP,i,strvlm[i]);       

    // iterates until the NP model is fit
    upbp(npk,beta,EP,m,q,np,CN,model,eps,&iter,linkchoice,maxit,rslpind,rslp,tol,N,y,npoind,T,                                                                                     
         &flagcvm,&flaginfo,Xijkl,hijkl,pijkl,Yij,Wik,SNP,FNP,spnp,fpnp,fmix,infonpml);

    // For calculating the value of the last mass point: sum_j=1^(npk-1)ai*pi
    double ath1[1+rslpind];
    for (a = 0; a < rslpind+1 ; a++) ath1[a] = 0;
    for (c = 0; c < rslpind+1 ; c++)
        for (a = 0; a < npk-1 ; a++)
            ath1[c] += gsl_vector_get(beta,q+np+(npk-1)*c+a) * gsl_vector_get(EP,a);

    // Output

    // beta
    for (i = 0; i< q+np; i++) out[i] = gsl_vector_get(beta,i);

    // Mass Points
    for (c = 0; c < rslpind+1 ; c++){
        for (i = 0; i< npk; i++){
            if (i < npk-1) out[q+np+npk*c+i] = gsl_vector_get(beta,q+np+(npk-1)*c+i);
            if (i==npk-1) out[q+np+npk*(c+1)-1] = -ath1[c]/gsl_vector_get(EP,npk-1);
        }
    }

    // Masses
    for (i = 0; i < npk; i++) out[q+np+(rslpind+1)*npk+i] = gsl_vector_get(EP,i);

    // Random effects covariance
    double rec[(1+rslpind)*(2+rslpind)/2];    
    for (c = 0; c < ((1+rslpind)*(2+rslpind)/2); c++) rec[c] = 0;
    
    for (c = 0; c < rslpind+1 ; c++){
        rec[c] = gsl_vector_get(EP,npk-1) * pow(ath1[c]/gsl_vector_get(EP,npk-1),2);
        for (i = 0; i < npk-1 ; i++){
            rec[c] += gsl_vector_get(EP,i) * pow(gsl_vector_get(beta,q+np+(npk-1)*c+i),2);
        }
    }

    int bc = 0;

    for (c = 0; c < rslpind ; c++){
        for (r = (c+1); r < rslpind+1 ; r++){
            rec[rslpind+1+bc] = ath1[c] * ath1[r] / gsl_vector_get(EP,npk-1);
            for (i = 0; i < npk-1 ; i++)
                rec[rslpind+1+bc] += gsl_vector_get(EP,i) * gsl_vector_get(beta,q+np+(npk-1)*c+i) * gsl_vector_get(beta,q+np+(npk-1)*r+i);
            bc += 1;
        }
    }

    for (i = 0; i < ((1+rslpind)*(2+rslpind)/2); i++)
        out[q+np+npk+rslpind*npk+npk+i] = rec[i];

    // Covariance Matrix based on the suqred derivative of the log-likelihood
    NPMLSE(npk,m,q,np,CN,beta,EP,model,linkchoice,rslpind,rslp,tol,N,y,npoind,T,Xijkl,hijkl,pijkl,Yij,Wik,infonpml);

    for (i = 0; i < (q+np+(rslpind+2)*(npk-1)); i++)
        for (j = 0; j < (q+np+(rslpind+2)*(npk-1)); j++)
            outparcvmat[i*(q+np+(rslpind+2)*(npk-1))+j] = gsl_matrix_get(infonpml,i,j);

    //Call delta method for Standard error of k^th mass point    
    if (npk > 1){        
        deltanpkp(q,np,npk,beta,EP,deltammp,rslpind);
        gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,deltammp,infonpml,0.0,DI); 
        gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,DI,deltammp,0.0,sekmp);
    }

    // Standard error of the last probability
    double selpk = 0;
    
    if (npk > 1){        
        for (i = 0; i < npk-1; i++) gsl_vector_set(selp,q+np+(npk-1)*(1+rslpind)+i,1);
        gsl_blas_dgemv(CblasTrans,1.0,infonpml,selp,0.0,is);
        gsl_blas_ddot(selp,is,&selpk);
    }

    // Standard errors of the q+np coefficients
    for (i = 0; i < (q+np); i++)
        out[q+np+npk+rslpind*npk+npk+((1+rslpind)*(2+rslpind)/2)+i] = sqrt(gsl_matrix_get(infonpml,i,i));

    // Standard errors of Mass points
    for (c = 0; c < rslpind+1 ; c++){
        for (i = 0; i < npk; i++){
            if (i < npk-1) out[q+np+npk+rslpind*npk+npk+((1+rslpind)*(2+rslpind)/2)+q+np+npk*c+i] = 
                         sqrt(gsl_matrix_get(infonpml,q+np+(npk-1)*c+i,q+np+(npk-1)*c+i));
            if (i == npk-1) out[q+np+npk+rslpind*npk+npk+((1+rslpind)*(2+rslpind)/2)+q+np+npk*c+i] = 
                          sqrt(gsl_matrix_get(sekmp,c,c));
        }
    }

    // Standard errors of Masses
    for (i = 0; i < npk-1; i++)
        out[q+np+npk+rslpind*npk+npk+((1+rslpind)*(2+rslpind)/2)+q+np+npk*(rslpind+1)+i] = 
        sqrt(gsl_matrix_get(infonpml,q+np+(npk-1)*(1+rslpind)+i,q+np+(npk-1)*(1+rslpind)+i));
    
    out[q+np+npk+rslpind*npk+npk+((1+rslpind)*(2+rslpind)/2)+q+np+npk*(rslpind+1)+npk-1] = sqrt(selpk);

    //Call delta method for Standard errors of var(ai), var(bi), ...
    
    deltanp(q,np,npk,beta,EP,deltam,rslpind);    //updates matrix deltam for all random effects
    
    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,deltam,infonpml,0.0,DI); 
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,DI,deltam,0.0,sevar);

    for (u = 0; u < (1+rslpind); u++)
        out[q+np+npk+rslpind*npk+npk+((1+rslpind)*(2+rslpind)/2)+q+np+npk*(rslpind+1)+npk+u] = 
        sqrt(gsl_matrix_get(sevar,u,u));

    // log likelihood
    out[q+np+npk+rslpind*npk+npk+((1+rslpind)*(2+rslpind)/2)+q+np+npk*(rslpind+1)+npk+rslpind+1] = 
    logL(CN,m,q,np,npk,EP,model,beta,linkchoice,rslpind,rslp,N,y,npoind,T,Xijkl,hijkl,pijkl,Yij);

    //Number of iterations
    out[q+np+npk+rslpind*npk+npk+((1+rslpind)*(2+rslpind)/2)+q+np+npk*(rslpind+1)+npk+rslpind+2] = iter;

    //Flag multinomial covariance matrix
    out[q+np+npk+rslpind*npk+npk+((1+rslpind)*(2+rslpind)/2)+q+np+npk*(rslpind+1)+npk+rslpind+3] = flagcvm;

    //Flag information matrix
    out[q+np+npk+rslpind*npk+npk+((1+rslpind)*(2+rslpind)/2)+q+np+npk*(rslpind+1)+npk+rslpind+4] = flaginfo;

    // Empirical Bayes: E(bi|Yi)
    if ((EBind == 1) & (npk > 1)){
        double denEB;
        double numEB[1+rslpind];
        for (i = 0; i < m; i++){
            denEB = 0.0;
            for (h = 0; h < (1+rslpind); h++) numEB[h] = 0.0;
            for (j = 0; j < npk; j++){
                fikl = setfikl(i,CN,q,np,npk,j+1,EP,model,beta,linkchoice,rslpind,rslp,N,y,npoind,T,Xijkl,hijkl,pijkl,Yij);
                denEB += fikl*gsl_vector_get(EP,j);
                for (h = 0; h < (1+rslpind); h++)
                    numEB[h] += fikl*gsl_vector_get(EP,j)*out[q+np+npk*h+j];
            }
            for (h = 0; h < (1+rslpind); h++) outEB[i+m*h] = numEB[h] / denEB;
        }
    }

    // Fitted values for the linear predictor at random effects equal to zero if EB=FALSE or at their EB estimate if EB=TRUE    
    for (i = 0; i < (q+np); i++) 
        gsl_vector_set(betanew,i,gsl_vector_get(beta,i));
    for (i = 0; i < m; i++){
        for (j = 0; j < (CN[i+1]-CN[i]); j++){
            if (npk > 1){
                for (h = 0; h < (1+rslpind); h++)
                    gsl_vector_set(betanew,q+np+(npk-1)*h,outEB[i+m*h]);
            }            
            sethijkl(i, j, CN, q, np, npk, 1, EP, model, betanew, rslpind, rslp, N, npoind, T, Xijkl, hijkl);            
            for (u = 0; u < q; u++)
                outFitted[(CN[i]+j)*q+u] = gsl_vector_get(hijkl,u);
        }
    }
 
    //Estimated category probabilities at random effects equal to zero if EB = FALSE or at their EB estimate if EB = TRUE
    for (i = 0; i < m; i++){
        for (j = 0; j < (CN[i+1]-CN[i]); j++){
            //clogit link
            if (linkchoice == 1){
                outProb[(CN[i]+j)*k] = exp(outFitted[(CN[i]+j)*q])/(1+exp(outFitted[(CN[i]+j)*q]));
                for (c = 1; c < q; c++)
                    outProb[(CN[i]+j)*k+c] = exp(outFitted[(CN[i]+j)*q+c])/(1+exp(outFitted[(CN[i]+j)*q+c])) - exp(outFitted[(CN[i]+j)*q+c-1])/(1+exp(outFitted[(CN[i]+j)*q+c-1]));
                outProb[(CN[i]+j)*k+q] = 1-exp(outFitted[(CN[i]+j)*q+q-1])/(1+exp(outFitted[(CN[i]+j)*q+q-1]));
            }
            //blogit link: baseline category is the last category
            if (linkchoice == 0){
                double partot = 1;
                for (c = 0; c < q; c++)
                    partot += exp(outFitted[(CN[i]+j)*q+c]);
                for (c = 0; c < q; c++)
                    outProb[(CN[i]+j)*k+c] = exp(outFitted[(CN[i]+j)*q+c])/partot;
                outProb[(CN[i]+j)*k+q] = 1 / partot;
            }
        }
    }                           
    gsl_vector_free(beta); gsl_vector_free(EP); gsl_vector_free(Yij); gsl_vector_free(hijkl); 
    gsl_vector_free(pijkl); gsl_vector_free(Wik); gsl_vector_free(SNP); gsl_vector_free(spnp);
    gsl_vector_free(selp); gsl_vector_free(is); gsl_vector_free(betanew);
    gsl_matrix_free(Xijkl); gsl_matrix_free(Sijkl); gsl_matrix_free(Dijkl); gsl_matrix_free(FNP);
    gsl_matrix_free(fpnp); gsl_matrix_free(fmix); gsl_matrix_free(infonpml); gsl_matrix_free(deltammp);
    gsl_matrix_free(deltam); gsl_matrix_free(sevar); gsl_matrix_free(sekmp); gsl_matrix_free(DI);        
}
