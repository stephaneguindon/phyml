/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include <config.h>

#ifndef RATES_H
#define RATES_H

#include "utilities.h"
#include "spr.h"
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
#include "time.h"
#include "m4.h"
#include "draw.h"
#include "mcmc.h"
#include "stats.h"

void RATES_Monte_Carlo_Mean_Rates(t_tree *tree);
void RATES_Monte_Carlo_Mean_Rates_Pre(t_node *a, t_node *d, t_edge *b, phydbl curr_rate, t_tree *tree);
void RATES_Print_Rates(t_tree *tree);
void RATES_Print_Rates_Pre(t_node *a, t_node *d, t_edge *b, t_tree *tree);
t_rate *RATES_Make_Rate_Struct(int n_otu);
void RATES_Init_Rate_Struct(t_rate *rates, t_rate *existing_rates, int n_otu);
void RATES_Classify_Branches(t_tree *tree);
void RATES_Adjust_Rates(t_tree *tree);
void RATES_Adjust_Rates_Local_Pre(t_node *a, t_node *d, t_edge *b, t_tree *tree);
void RATES_Adjust_Rates_Local(t_node *a, t_node *d, t_edge *b1, t_tree *tree);
void RATES_Record_T(t_tree *tree);
void RATES_Restore_T(t_tree *tree);
void RATES_Monte_Carlo_Mean_Rates_Core(phydbl t_lim_sup, phydbl t_lim_inf, phydbl *curr_rate, phydbl *mean_rate, phydbl lexp, phydbl alpha);
phydbl RATES_Lk_Rates(t_tree *tree);
void RATES_Lk_Rates_Pre(t_node *a, t_node *d, t_edge *b, t_tree *tree);
void RATES_Fill_Node_Rates_Pre(t_node *a, t_node *d, t_edge *b, phydbl *node_r, t_tree *tree);
void RATES_Fill_Node_Rates(phydbl *node_r, t_tree *tree);
void RATES_Optimize_Node_Times_Serie_Fixed_Br_Len(t_node *a, t_node *d, t_tree *tree);
void RATES_Optimize_Lexp(t_tree *tree);
void RATES_Round_Optimize(t_tree *tree);
void RATES_Optimize_Lexp(t_tree *tree);
void RATES_Optimize_Alpha(t_tree *tree);
phydbl RATES_Dmu(phydbl mu, int n_jumps, phydbl dt, phydbl a, phydbl b, phydbl lexp, int min_n, int jps_dens);
phydbl RATES_Dr_X_Dx(phydbl r, phydbl mu, phydbl y, phydbl dt, phydbl a, phydbl b, phydbl lexp);
phydbl RATES_Dmu_Given_Y_Trpzd(phydbl mu, phydbl y, phydbl dt, phydbl a, phydbl b, phydbl lexp, 
			       int nsteps, phydbl beg, phydbl end, phydbl prevs);
phydbl RATES_Dmu_Given_Y_Std(phydbl mu, phydbl y, phydbl dt, phydbl a, phydbl b, phydbl lexp);
phydbl RATES_Dmu_Given_Y_Romb(phydbl mu, phydbl y, phydbl dt, phydbl a, phydbl b, phydbl lexp);
phydbl RATES_Dmu_Given_Y(phydbl mu, phydbl y, phydbl dt, phydbl a, phydbl b, phydbl lexp);
phydbl RATES_Dy_Given_Mu(phydbl mu, phydbl y, phydbl dt, phydbl a, phydbl b, phydbl lexp);
phydbl RATES_Dmu2_Given_Y_X_Dy_Given_Mu1(phydbl mu1, phydbl mu2, phydbl y, phydbl dt1, phydbl dt2, phydbl a, phydbl b, phydbl lexp);
phydbl RATES_Dmu2_Given_Mu1_Trpz(phydbl mu1, phydbl mu2, phydbl dt1, phydbl dt2, phydbl a, phydbl b, phydbl lexp,
				 int nsteps, phydbl beg, phydbl end, phydbl prevs);
phydbl RATES_Dmu2_And_Mu1(phydbl mu1, phydbl mu2, phydbl dt1, phydbl dt2, phydbl a, phydbl b, phydbl lexp);
phydbl RATES_Dmu2_Given_Mu1(phydbl mu1, phydbl mu2, phydbl dt1, phydbl dt2, phydbl a, phydbl b, phydbl lexp);
phydbl RATES_Dmu2_Given_Mu1_Romb(phydbl mu1, phydbl mu2, phydbl dt1, phydbl dt2, phydbl a, phydbl b, phydbl lexp);
void RATES_Random_Branch_Lengths(t_tree *tree);
void RATES_Bracket_N_Jumps(int *up, int *down, phydbl param);
void RATES_Set_Node_Times(t_tree *tree);
void RATES_Set_Node_Times_Pre(t_node *a, t_node *d, t_tree *tree);
void RATES_Optimize_Node_Times(t_tree *tree);
phydbl RATES_Exp_Y(phydbl mu1, phydbl mu2, phydbl dt1, phydbl lexp);
phydbl RATES_Dmu2_Given_Mu1_Bis(phydbl mu1, phydbl mu2, phydbl dt1, phydbl dt2, phydbl alpha, phydbl beta, phydbl lexp);
void RATES_Replace_Br_Lengths_By_Rates_Pre(t_node *a, t_node *d, t_edge *b, t_tree *tree);
void RATES_Replace_Br_Lengths_By_Rates(t_tree *tree);

void RATES_Get_Mean_Rates(t_tree *tree);
void RATES_Get_Mean_Rates_Pre(t_node *a, t_node *d, t_edge *b, phydbl r_a, t_tree *tree);
void RATES_Expect_Number_Subst(phydbl t_beg, phydbl t_end, phydbl r_beg,  int *n_jumps, phydbl *mean_r, phydbl *r_end, t_rate *rates, t_tree *tree);
void RATES_Optimize_Clock_Rate(t_tree *tree);
phydbl RATES_Dmu1_Given_Lbda_And_Mu2(phydbl lbda, phydbl mu1, phydbl mu2, phydbl alpha, phydbl beta);
phydbl RATES_Dmu1_And_Mu2_One_Jump_Trpz(phydbl mu1, phydbl mu2, phydbl a, phydbl b,
					  int nsteps, phydbl beg, phydbl end, phydbl prevs);
phydbl RATES_Dmu1_And_Mu2_One_Jump_One_Interval(phydbl mu1, phydbl mu2, phydbl a, phydbl b);
phydbl RATES_Dmu1_And_Mu2_One_Jump_Two_Intervals(phydbl dt1, phydbl dt2, phydbl mu1, phydbl mu2, phydbl a, phydbl b);

phydbl RATES_Dmu1_And_Mu2_One_Jump_Old(phydbl mu1, phydbl mu2, phydbl a, phydbl b);

phydbl RATES_Dmu2_And_Min_N_Given_Mu1(phydbl mu1, phydbl mu2, phydbl dt1, phydbl dt2, int n_min, phydbl a, phydbl b, phydbl lexp);
phydbl RATES_Dmu2_And_Mu1_Given_N(phydbl mu1, phydbl mu2, phydbl dt1, phydbl dt2, int n, phydbl a, phydbl b, phydbl lexp);
phydbl RATES_Lk_Rates_Core(phydbl br_r_a, phydbl br_r_d, phydbl nd_r_a, phydbl nd_r_d, int n_a, int n_d, phydbl dt_a, phydbl dt_d, t_tree *tree);
void RATES_Init_Triplets(t_tree *tree);
phydbl RATES_Lk_Change_One_Time(t_node *n, phydbl new_t, t_tree *tree);
void RATES_Update_Triplet(t_node *n, t_tree *tree);
void RATES_Print_Triplets(t_tree *tree);
phydbl RATES_Lk_Change_One_Rate(t_node *d, phydbl new_rate, t_tree *tree);
phydbl RATES_Dmu2_And_Mu1_Given_Min_N(phydbl mu1, phydbl mu2, phydbl dt1, phydbl dt2, int n_min, phydbl a, phydbl b, phydbl lexp);
phydbl RATES_Dmu2_And_Mu1_Given_N_Normal(phydbl mu1, phydbl mu2, phydbl dt1, phydbl dt2, int n, phydbl a, phydbl b, phydbl lexp);
phydbl RATES_Coeff_Corr(phydbl alpha, phydbl beta, int n1, int n2);
phydbl RATES_Dmu2_And_Mu1_Given_N_Full(phydbl mu1, phydbl mu2, phydbl dt1, phydbl dt2, int n, phydbl a, phydbl b, phydbl lexp);
phydbl RATES_Dmu1_Given_V_And_N(phydbl mu1, phydbl v, int n, phydbl dt1, phydbl a, phydbl b);
phydbl RATES_Yule(t_tree *tree);
phydbl RATES_Check_Mean_Rates(t_tree *tree);
void RATES_Check_Mean_Rates_Pre(t_node *a, t_node *d, t_edge *b, phydbl *sum, t_tree *tree);
void RATES_Discretize_Rates(t_tree *tree);
void RATES_Discretize_Rates_Pre(t_node *a, t_node *d, t_edge *b, t_tree *tree);
phydbl RATES_Dmu_Given_V_And_MinN(phydbl mu, phydbl dt, phydbl v, int minn, phydbl a, phydbl b, phydbl lexp);
phydbl RATES_Dmu_One(phydbl mu, phydbl dt, phydbl a, phydbl b, phydbl lexp);
phydbl RATES_Compound_Core(phydbl mu1, phydbl mu2, int n1, int n2, phydbl dt1, phydbl dt2, phydbl alpha, phydbl beta, phydbl lexp, phydbl eps, int approx);
void RATES_Record_Rates(t_tree *tree);
void RATES_Reset_Rates(t_tree *tree);
void RATES_Record_Times(t_tree *tree);
void RATES_Reset_Times(t_tree *tree);
void RATES_Update_T_Rates_Pre(t_node *a, t_node *d, t_tree *tree);
void RATES_Update_T_Rates(t_tree *tree);
void RATES_Get_Rates_From_Bl(t_tree *tree);
phydbl RATES_Compound_Core_Joint(phydbl mu1, phydbl mu2, int n1, int n2, phydbl dt1, phydbl dt2, 
				 phydbl alpha, phydbl beta, phydbl lexp, phydbl eps, int approx);
phydbl RATES_Dmu_Joint(phydbl mu, int n, phydbl dt, phydbl a, phydbl b, phydbl lexp, int min_n);
phydbl RATES_Compound_Core_Marginal(phydbl mu1, phydbl mu2, phydbl dt1, phydbl dt2, phydbl alpha, 
				    phydbl beta, phydbl lexp, phydbl eps, int approx);
phydbl RATES_Lk_Jumps(t_tree *tree);
void RATES_Posterior_Rates(t_tree *tree);
void RATES_Posterior_One_Rate(t_node *a, t_node *d, int traversal, t_tree *tree);
void RATES_Free_Rates(t_rate *rates);
void RATES_Initialize_True_Rates(t_tree *tree);
void RATES_Posterior_Times(t_tree *tree);
void RATES_Posterior_One_Time(t_node *a, t_node *d, int traversal, t_tree *tree);
void RATES_Update_Cur_Bl(t_tree *tree);
void RATES_Update_Cur_Bl_Pre(t_node *a, t_node *d, t_edge *b, t_tree *tree);
void RATES_Get_Cov_Matrix_Rooted(phydbl *unroot_cov, t_tree *tree);
void RATES_Get_Cov_Matrix_Rooted_Pre(t_node *a, t_node *d, t_edge *b, phydbl *cov, t_tree *tree);
void RATES_Bl_To_Ml(t_tree *tree);
void RATES_Bl_To_Ml_Pre(t_node *a, t_node *d, t_edge *b, t_tree *tree);
int RATES_Check_Node_Times(t_tree *tree);
void RATES_Check_Node_Times_Pre(t_node *a, t_node *d, int *err, t_tree *tree);
void RATES_Covariance_Mu(t_tree *tree);
void RATES_Variance_Mu_Pre(t_node *a, t_node *d, t_tree *tree);
void RATES_Fill_Lca_Table(t_tree *tree);
void RATES_Posterior_Clock_Rate(t_tree *tree);
void RATES_Get_Conditional_Variances(t_tree *tree);
void RATES_Get_All_Reg_Coeff(t_tree *tree);
void RATES_Posterior_Time_Root(t_tree *tree);
void RATES_Get_Trip_Conditional_Variances(t_tree *tree);
void RATES_Get_All_Trip_Reg_Coeff(t_tree *tree);
void RATES_Check_Lk_Rates(t_tree *tree, int *err);
phydbl RATES_Expected_Tree_Length(t_tree *tree);
void RATES_Expected_Tree_Length_Pre(t_node *a, t_node *d, phydbl eranc, phydbl *mean, int *n, t_tree *tree);
void RATES_Normalise_Rates(t_tree *tree);
phydbl RATES_Check_Mean_Rates_True(t_tree *tree);
phydbl RATES_Find_Min_Dt_Bisec(phydbl r0, phydbl r1, phydbl t0, phydbl t1, phydbl nu, phydbl threshp, int inf);
phydbl RATES_Find_Max_Dt_Bisec(phydbl r0, phydbl r1, phydbl t0, phydbl t1, phydbl nu, phydbl threshp, int inf);
void RATES_Min_Max_Interval(phydbl u0, phydbl u1, phydbl u2, phydbl u3, phydbl t0, phydbl t2, phydbl t3, 
			    phydbl *t_min, phydbl *t_max, phydbl nu, phydbl p_thresh, t_tree *tree);


phydbl RATES_Get_Correction_Factor(phydbl mode, phydbl sd, int *err, t_tree *tree);
phydbl RATES_Average_Substitution_Rate(t_tree *tree);
void RATES_Update_Norm_Fact(t_tree *tree);
void RATES_Update_Mean_Br_Len(int iter, t_tree *tree);
void RATES_Update_Cov_Br_Len(int iter, t_tree *tree);
void RATES_Set_Mean_L(t_tree *tree);
void RATES_Fill_All_Param(t_rate *rate, t_tree *tree);
void RATES_Record_Rates(t_tree *tree);
void RATES_Reset_Rates(t_tree *tree);
phydbl RATES_Average_Rate(t_tree *tree);
void RATES_Set_Clock_And_Nu_Max(t_tree *tree);
void RATES_Write_Mean_R_On_Edge_Label(t_node *a, t_node *d, t_edge *b, t_tree *tree);
phydbl RATES_Lk_Linreg(t_tree *tree);
phydbl RATES_Get_Mean_Rate_In_Subtree(t_node *root, t_tree *tree);
void RATES_Get_Mean_Rate_In_Subtree_Pre(t_node *a, t_node *d, phydbl *sum, int *n, t_tree *tree);
char *RATES_Get_Model_Name(int model);
void RATES_Get_Survival_Ranks(t_tree *tree);
void RATES_Bl_To_Bl(t_tree *tree);
void RATES_Bl_To_Bl_Pre(t_node *a, t_node *d, t_edge *b, t_tree *tree);
void RATES_Set_Birth_Rate_Boundaries(t_tree *tree);

#endif
