/*

PHYML :  a program that  computes maximum likelihood  phylogenies from
DNA or AA homologous sequences 

Copyright (C) Stephane Guindon. Oct 2003 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/

#include <config.h>

#ifndef MODELS_H
#define MODELS_H

#include "utilities.h"
#include "eigen.h"
#include "free.h"
#include "stats.h"
#include "mixt.h"

void PMat(phydbl l, t_mod *mod, int pos, phydbl *Pij, phydbl *tPij);
void  PMat_K80(phydbl l,phydbl kappa, int pos, phydbl *Pij);
void  PMat_TN93(phydbl l, t_mod *mod, int pos, phydbl *Pij);
void  PMat_Empirical(phydbl l, const t_mod *mod, const int pos, phydbl *Pij, phydbl *tPij);
void PMat_Zero_Br_Len(t_mod *mod, int pos, phydbl *Pij);
void PMat_Gamma(phydbl l, t_mod *mod, int pos, phydbl *Pij);
int GetDaa (phydbl *daa, phydbl *pi, char *file_name);
void Update_Qmat_GTR(phydbl *rr, phydbl *rr_val, int *rr_num, phydbl *pi, phydbl *qmat);
void Update_Qmat_HKY(phydbl kappa, phydbl *pi, phydbl *qmat);
phydbl Update_Qmat_Generic(phydbl *rr, phydbl *pi, int ns, phydbl *qmat);
void Translate_Custom_Mod_String(t_mod *mod);
int Set_Model_Parameters(t_mod *mod);
phydbl GTR_Dist(phydbl *F, phydbl alpha, eigen *eigen_struct);
phydbl General_Dist(phydbl *F, t_mod *mod, eigen *eigen_struct);
void Switch_From_Mod_To_M4mod(t_mod *mod);
void Switch_From_M4mod_To_Mod(t_mod *mod);
void PMat_JC69(phydbl l, int pos, phydbl *Pij, t_mod *mod);
phydbl Get_Lambda_F84(phydbl *pi, phydbl *kappa);
int Update_Eigen(t_mod *mod);
int Update_RAS(t_mod *mod);
int Update_Efrq(t_mod *mod);
void PMat_MGF_Gamma(phydbl mu, phydbl sigsq, const t_mod *mod, const int pos, phydbl *Pij, phydbl *tPij);
int Update_Boundaries(t_mod *mod);
void Update_Qmat_TN93(phydbl kappa, phydbl lambda, phydbl *pi, phydbl *qmat);

#endif
