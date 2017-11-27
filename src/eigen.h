#include <config.h>

#ifndef EIGEN_H
#define EIGEN_H

#include "utilities.h"
#include "free.h"

#ifdef RWRAPPER
#include <R.h>
#endif

int ludcmp_1D(phydbl *a, int n, phydbl *d);
void det_1D(phydbl *a, int n, phydbl *d);

int Eigen(int job, phydbl *A, int n, phydbl *rr, phydbl *ri,
          phydbl *vr, phydbl *vi, phydbl *w);
void balance(phydbl *mat, int n, int *low, int *hi, phydbl *scale);
void unbalance(int n, phydbl *vr, phydbl *vi, int low, int hi,
               phydbl *scale);
int realeig(int job, phydbl *mat, int n,int low, int hi, phydbl *valr,
            phydbl *vali, phydbl *vr, phydbl *vi);
void elemhess(int job, phydbl *mat, int n, int low, int hi, 
            phydbl *vr, phydbl *vi, int *work);

int ludcmp(phydbl **a, int n, phydbl *d);
void det(phydbl **a, int n, phydbl *d);

/* complex functions */

typedef struct { phydbl re, im; } complex;
#define csize(a) (FABS(a.re)+FABS(a.im))

complex compl (phydbl re,phydbl im);
complex _conj (complex a);
complex cplus (complex a, complex b);
complex cminus (complex a, complex b);
complex cby (complex a, complex b);
complex cdiv (complex a,complex b);
/* complex local_cexp (complex a); */
complex cfactor (complex x, phydbl a);
int cxtoy (complex *x, complex *y, int n);
int cmatby (complex *a, complex *b, complex *c, int n,int m,int k);
int cmatout (FILE * fout, complex *x, int n, int m);
int cmatinv( complex *x, int n, int m, phydbl *space);
phydbl *Cholesky_Decomp(phydbl *A,  int dim);


#endif

