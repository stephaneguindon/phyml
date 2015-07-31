/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include <config.h>

#ifndef MCMC_H
#define MCMC_H

#include "spr.h"
#include "utilities.h"
#include "lk.h"
#include "optimiz.h"
#include "bionj.h"
#include "models.h"
#include "free.h"
#include "help.h"
#include "simu.h"
#include "eigen.h"
#include "pars.h"
#include "alrt.h"
#include "times.h"
#include "m4.h"
#include "draw.h"
#include "rates.h"
#include "stats.h"
#include "phyrex.h"
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <signal.h>

void MCMC_Lexp(t_tree *tree);
void MCMC_Print_Param(t_mcmc *mcmc, t_tree *tree);
t_mcmc *MCMC_Make_MCMC_Struct();
void MCMC_Free_MCMC(t_mcmc *mcmc);
void MCMC(t_tree *tree);
void MCMC_Alpha(t_tree *tree);
void MCMC_Randomize_Branch_Lengths(t_tree *tree);
void MCMC_Randomize_Node_Times(t_tree *tree);
void MCMC_Randomize_Node_Times_Pre(t_node *a, t_node *d, t_tree *tree);
void MCMC_Randomize_Lexp(t_tree *tree);
void MCMC_Randomize_Jumps(t_tree *tree);
void MCMC_Randomize_Alpha(t_tree *tree);
void MCMC_One_Rate(t_node *a, t_node *d, int traversal, t_tree *tree);
void MCMC_No_Change(t_tree *tree);
void MCMC_Nu(t_tree *tree);
void MCMC_Randomize_Nu(t_tree *tree);
t_node *MCMC_Select_Random_Node_Pair(phydbl t_sup, t_tree *tree);
void MCMC_Modify_Rates(t_tree *tree);
void MCMC_Modify_Subtree_Rate(t_node *a, t_node *d, phydbl new_rate, t_tree *tree);
void MCMC_Randomize_Rates(t_tree *tree);
void MCMC_Stick_Rates(t_tree *tree);
void MCMC_Stick_Rates_Pre(t_node *a, t_node *d, t_tree *tree);
void MCMC_Times_Global(t_tree *tree);
void MCMC_Times_Local(t_tree *tree);
void MCMC_One_Time(t_node *a, t_node *d, int traversal, t_tree *tree);
void MCMC_Rates_Global(t_tree *tree);
void MCMC_Rates_Local(t_tree *tree);
void MCMC_Rates_Pre(t_node *a, t_node *d, t_tree *tree);
void MCMC_Mixing_Step(t_tree *tree);
void MCMC_Jumps_Local(t_tree *tree);
void MCMC_Jumps_Pre(t_node *a, t_node *d, int local, t_tree *tree);
void MCMC_Randomize_Clock_Rate(t_tree *tree);
void MCMC_Clock_Rate(t_tree *tree);
void MCMC_Time_Root(t_tree *tree);
void MCMC_Randomize_Node_Times_Bottom_Up(t_node *a, t_node *d, t_tree *tree);
void MCMC_Randomize_Node_Times_Top_Down(t_node *a, t_node *d, t_tree *tree);
void MCMC_Randomize_Rates_Pre(t_node *a, t_node *d, t_tree *tree);
void MCMC_Print_Means(t_mcmc *mcmc, t_tree *tree);
void MCMC_Print_Last(t_mcmc *mcmc, t_tree *tree);
void MCMC_Close_MCMC(t_mcmc *mcmc);
void MCMC_Rates_Global(t_tree *tree);
void MCMC_Omega(t_tree *tree);
void MCMC_Adjust_Tuning_Parameter(int move, t_mcmc *mcmc);
void MCMC_Copy_MCMC_Struct(t_mcmc *ori, t_mcmc *cpy, char *filename);
void MCMC_Randomize_Node_Times_Bottom_Up(t_node *a, t_node *d, t_tree *tree);
void MCMC_One_Length(t_edge *b,  t_tree *tree);
void MCMC_Br_Lens(t_tree *tree);
void MCMC_Br_Lens_Pre(t_node *a, t_node *d, t_edge *b, t_tree *tree);
void MCMC_Tree_Height(t_tree *tree);
void MCMC_Subtree_Height(t_tree *tree);
void MCMC_Swing(t_tree *tree);
void MCMC_Single_Param_Generic(phydbl *val, 
			       phydbl lim_inf, 
			       phydbl lim_sup, 
			       int move_num,
			       phydbl *lnPrior,
			       phydbl *lnLike,
			       phydbl (*prior_func)(t_edge *,t_tree *,supert_tree *), 
			       phydbl (*like_func)(t_edge *,t_tree *,supert_tree *),
			       int move_type, int _log,
			       t_edge *branch, t_tree *tree, supert_tree *stree);
void MCMC_Scale_Br_Lens(t_tree *tree);
void MCMC_Update_Mean_Br_Len(t_tree *tree);
void MCMC_Update_Cov_Br_Len(t_tree *tree);
void MCMC_Sim_Rate(t_node *a, t_node *d, t_tree *tree);
void Fill_All_Param(t_mcmc *mcmc, t_rate *rate, t_tree *tree);
int Get_Param_Num(t_mcmc *mcmc, phydbl *param);
void MCMC_Complete_MCMC(t_mcmc *mcmc, t_tree *tree);
void MCMC_Sample_Joint_Rates_Prior(t_tree *tree);
void MCMC_Sample_Joint_Rates_Posterior(t_tree *tree);
void MCMC_Pair_Rates_Constraint(t_node *a, t_node *d, int random, int traversal, t_tree *tree);
void MCMC_Times_And_Rates(t_node *a, t_node *d, int random, int traversal, t_tree *tree);
void MCMC_Tree_Rates(t_tree *tree);
void MCMC_Pause(t_mcmc *mcmc);
void MCMC_Print_Param_Stdin(t_mcmc *mcmc, t_tree *tree);
void MCMC_Subtree_Rates(t_tree *tree);
void MCMC_Get_Acc_Rates(t_mcmc *mcmc);
void MCMC_Update_Effective_Sample_Size(int move_num, t_mcmc *mcmc, t_tree *tree);
void MCMC_Initialize_Param_Val(t_mcmc *mcmc, t_tree *tree);
void MCMC_Terminate();
void MCMC_Copy_To_New_Param_Val(t_mcmc *mcmc, t_tree *tree);
void MCMC_Randomize_Node_Rates(t_tree *tree);
void MCMC_One_Node_Rate(t_node *a, t_node *d, int traversal, t_tree *tree);
void MCMC_Tree_Rates_Bis(t_tree *tree);
void MCMC_Slice_One_Rate(t_node *a, t_node *d, int traversal, t_tree *tree);
void MCMC_Updown_Nu_Cr(t_tree *tree);
void MCMC_All_Rates(t_tree *tree);
void MCMC_Alpha(t_tree *tree);
void MCMC_Kappa(t_tree *tree);
void MCMC_Rate_Across_Sites(t_tree *tree);
void MCMC_Free_Mixt_Rate(t_tree *tree);
void MCMC_Make_Move(phydbl *cur, phydbl *new, phydbl inf, phydbl sup, phydbl *loghr, phydbl tune, int move_type);
void MCMC_Randomize_Rate_Across_Sites(t_tree *tree);
void MCMC_Randomize_Kappa(t_tree *tree);
void MCMC_Updown_T_Cr(t_tree *tree);
void MCMC_Linreg_Par(t_tree *tree);
void MCMC_Covarion_Rates(t_tree *tree);
void MCMC_Covarion_Switch(t_tree *tree);
void MCMC_Randomize_Covarion_Rates(t_tree *tree);
void MCMC_Randomize_Covarion_Switch(t_tree *tree);
void MCMC_Read_Param_Vals(t_tree *tree);
void MCMC_Birth_Rate(t_tree *tree);
void MCMC_Randomize_Birth(t_tree *tree);
void MCMC_Clock_R(t_tree *mixt_tree);
void MCMC_Free_MCMC(t_mcmc *mcmc);
void MCMC_Updown_T_Br(t_tree *tree);
void MCMC_Root_Time(t_tree *tree);
void MCMC_Jump_Calibration(t_tree *tree);
void MCMC_GEO_Lbda(t_tree *mixt_tree);
void MCMC_GEO_Sigma(t_tree *mixt_tree);
void MCMC_GEO_Tau(t_tree *mixt_tree);
void MCMC_GEO_Loc(t_tree *tree);
void MCMC_GEO_Dum(t_tree *mixt_tree);
void MCMC_PHYREX_Lbda(t_tree *mixt_tree);
void MCMC_PHYREX_Mu(t_tree *mixt_tree);
void MCMC_PHYREX_Radius(t_tree *mixt_tree);
void MCMC_PHYREX_Triplet(t_tree *tree);
void MCMC_PHYREX_Move_Disk_Centre(t_tree *tree);
void MCMC_PHYREX_Move_Disk_Updown(t_tree *tree);
void MCMC_PHYREX_Swap_Disk(t_tree *tree);
void MCMC_PHYREX_Move_Ldsk(t_tree *tree);
void MCMC_PHYREX_Prune_Regraft(t_tree *tree);
void MCMC_PHYREX_Scale_Times(t_tree *tree);
void MCMC_PHYREX_Ldscape_Limits(t_tree *tree);
void MCMC_PHYREX_Insert_Disk(phydbl hr, int n_insert_disks, t_tree *tree);
void MCMC_PHYREX_Delete_Disk(phydbl hr, int n_delete_disks, t_tree *tree);
void MCMC_PHYREX_Indel_Disk(t_tree *tree);
void MCMC_PHYREX_Insert_Hit(phydbl hr, int n_insert_disks, t_tree *tree);
void MCMC_PHYREX_Delete_Hit(phydbl hr, int n_delete_disks, t_tree *tree);
void MCMC_PHYREX_Indel_Hit(t_tree *tree);
void MCMC_PHYREX_Simulate_Backward(t_tree *tree);
void MCMC_Update_Mode(int move_num, t_mcmc *mcmc, t_tree *tree);
void MCMC_PHYREX_Lineage_Traj(t_tree *tree);

#endif
