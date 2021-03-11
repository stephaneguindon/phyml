/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include <config.h>

#ifndef SLFV_H
#define SLFV_H

#include "utilities.h"

phydbl SLFV_Prob_Two_Lineages_Coal(t_ldsk *l0, t_ldsk *l1, t_tree *tree);
phydbl SLFV_Prob_Two_Random_Lineages_Coal_One_Event(phydbl w, phydbl h, phydbl mu, phydbl rad);
phydbl SLFV_Coalescence_Rate(t_tree *tree);
phydbl SLFV_Path_Logdensity(t_ldsk *beg, t_ldsk *end, phydbl *sd, t_tree *tree);
void SLFV_Sample_Path(t_ldsk *young, t_ldsk *old, phydbl *sd, phydbl *global_hr, t_tree *tree);
t_ldsk *SLFV_Generate_Path(t_ldsk *beg, t_ldsk *end, int n_evt, phydbl *sd, t_tree *tree);
phydbl SLFV_Effective_Density(t_tree *tree);
phydbl SLFV_Rate_Per_Unit_Area(t_tree *tree);
phydbl SLFV_Sample_Rad_From_Prior(t_tree *tree);
phydbl SLFV_Update_Radius(t_tree *tree);
phydbl SLFV_Generation_Length(t_tree *tree);
phydbl SLFV_Neighborhood_Size(t_tree *tree);
phydbl SLFV_Update_Sigsq(t_tree *tree);
phydbl SLFV_Neighborhood_Size_Regression(t_tree *tree);
phydbl SLFV_Lk_Gaussian_Range(t_dsk *young, t_dsk *old, t_tree *tree);
phydbl SLFV_Lk_Gaussian(t_tree *tree);
phydbl SLFV_Lk_Gaussian_Core(t_dsk *disk, t_tree *tree);
t_sarea *SLFV_Simulate_Forward_Core(int n_sites, t_tree *tree);
phydbl SLFV_Simulate_Backward_Core(t_dsk *init_disk, int avoid_multiple_mergers, t_tree *tree);
t_tree *SLFV_Simulate(int n_otu, int n_sites, phydbl w, phydbl h, phydbl  lbda, phydbl rad, phydbl mu, int r_seed);
t_tree *SLFV_Simulate_Independent_Loci(int n_otu, int n_loci, phydbl w, phydbl h, int r_seed);
void SLFV_Integrated_Coal_Rate(t_ldsk *l0, t_ldsk *l1, phydbl T, t_tree *tree);
void SLFV_Sum_Coal_Rate(t_ldsk *l0, t_ldsk *l1, int n, t_tree *tree);
phydbl SLFV_Lk_Range(t_dsk *young, t_dsk *old, t_tree *tree);
void SLFV_Generate_Ldsk_New_Location(t_ldsk *l, t_ldsk *prev_l, phydbl rad, phydbl *hr, int dim_idx, t_tree *tree);

#endif





