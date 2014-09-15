/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include <config.h>

#ifndef RWRAPPER_H
#define RWRAPPER_H

#include "utilities.h"
#include "rates.h"
#include "eigen.h"
#include "io.h"
#include "stats.h"
#include "tiporder.h"

/* void RWRAPPER_Dmu2_Given_Mu1_And_Min_N(phydbl *mu1, phydbl *mu2, phydbl *dt1, phydbl *dt2, int *n_min, phydbl *a, phydbl *b, phydbl *lexp, phydbl *dens); */
/* void RWRAPPER_Dnorm(phydbl *x, phydbl *mean, phydbl *var, double *dens); */
/* void RWRAPPER_Dmu(phydbl *mu, phydbl *dt, phydbl *a, phydbl *b, phydbl *lexp, double *dens); */
/* void RWRAPPER_Dr_X_Dx(double *r, double *mu, double *y, double *dt, double *a, double *b, double *lexp, double *dens); */
/* void RWRAPPER_Dmu_Given_Y(double *mu, double *y, double *dt, double *a, double *b, double *lexp, double *dens); */
/* void RWRAPPER_Dmu(phydbl *mu, phydbl *dt, phydbl *a, phydbl *b, phydbl *lexp, double *dens); */
/* void RWRAPPER_Dmu2_Given_Mu1(double *mu1, double *mu2, double *dt1, double *dt2, double *a, double *b, double *lexp, double *dens); */
/* void RWRAPPER_Dgamma(phydbl *x, phydbl *shape, phydbl *scale, phydbl *dens); */
/* void RWRAPPER_Bivariate_Normal_Density(phydbl *x, phydbl *y, phydbl *mux, phydbl *muy, phydbl *sdx, phydbl *sdy, phydbl *rho, phydbl *dens); */
/* void RWRAPPER_Dmu2_And_Mu1_Given_Min_N(phydbl *mu1, phydbl *mu2, phydbl *dt1, phydbl *dt2, int *n_min, phydbl *a, phydbl *b, phydbl *lexp, phydbl *dens); */
/* void RWRAPPER_RATES_Dmu1_And_Mu2_One_Jump_Two_Intervals(phydbl *dt1, phydbl *dt2, phydbl *mu1, phydbl *mu2, phydbl *alpha, phydbl *beta, phydbl *dens); */
/* void RWRAPPER_Dmu_One(phydbl *mu, phydbl *dt, phydbl *a, phydbl *b, phydbl *lexp, double *dens); */
/* void RWRAPPER_Lk_Rates_Core(phydbl *mu1, phydbl *mu2, phydbl *dt1, phydbl *dt2, phydbl *a, phydbl *b, phydbl *lexp, phydbl *eps, int *approx, phydbl *dens); */
/* void RWRAPPER_Dmu2_And_Mu1(double *mu1, double *mu2, double *dt1, double *dt2, double *a, double *b, double *lexp, double *dens); */
/* void RWRAPPER_Rnorm_Trunc(phydbl *mean, phydbl *sd, phydbl *min, phydbl *max, phydbl *res); */
/* void  RWRAPPER_Cholesky_Decomp(double *A, int *dim); */
/* void RWRAPPER_Min_Number_Of_Tip_Permut(char **tree_file_name, char **coord_file_name, phydbl *res); */
void RWRAPPER_Log_Dnorm(phydbl *x, phydbl *mean, phydbl *sd,  phydbl *res);
void RWRAPPER_Integrated_Geometric_Brownian_Bridge_Mean(phydbl *T, phydbl *A, phydbl *B, phydbl *u, phydbl *mean);
void RWRAPPER_Integrated_Geometric_Brownian_Bridge_Var(phydbl *T, phydbl *A, phydbl *B, phydbl *u, phydbl *var);


#endif
