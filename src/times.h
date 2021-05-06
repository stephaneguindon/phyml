/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include <config.h>

#ifndef TIMES_H
#define TIMES_H

#include "utilities.h"

int  TIMES_Main(int argc, char **argv);
void TIMES_Bl_From_T_Post(t_node *a, t_node *d, t_edge *b, t_tree *tree);
void TIMES_Bl_From_T(t_tree *tree);
void TIMES_Optimize_Node_Times_Serie(t_node *a, t_node *d, t_tree *tree);
void TIMES_Round_Optimize(t_tree *tree);
void TIMES_Print_Node_Times(t_node *a, t_node *d, t_tree *tree);
t_edge *TIMES_Find_Best_Root_Position(t_tree *tree);
void TIMES_Least_Square_Node_Times(t_edge *e_root, t_tree *tree);
void TIMES_Least_Square_Node_Times_Pre(t_node *a, t_node *d, phydbl *A, phydbl *b, int n, t_tree *tree);
void TIMES_Mult_Time_Stamps(t_tree *tree);
void TIMES_Div_Time_Stamps(t_tree *tree);
void TIMES_Optimize_Tree_Height(t_tree *tree);
void TIMES_Adjust_Node_Times(t_tree *tree);
void TIMES_Adjust_Node_Times_Pre(t_node *a, t_node *d, t_tree *tree);
void TIMES_Optimize_Root_Height(t_tree *tree);
void TIMES_Estimate_Branch_Rates(t_tree *tree);
t_edge *TIMES_Find_Best_Root_Position_Approx(t_tree *tree);
void TIMES_Estimate_Branch_Rate_Parameter(t_tree *tree);
phydbl TIMES_Classify_Branch_In_Rate_Class(t_tree *tree);
void TIMES_Compute_Rates_And_Times_Least_Square_Adjustments(t_tree *tree);
void TIMES_Compute_Rates_And_Times_Least_Square_Adjustments_Post(t_node *a, t_node *d, t_edge *b, t_tree *tree);
void TIMES_Classify_Branch_Rates(t_tree *tree);
int TIMES_Check_MC(t_tree *tree);
void TIMES_Set_All_Node_Priors(t_tree *tree);
void TIMES_Set_All_Node_Priors_Bottom_Up(t_node *a, t_node *d, t_tree *tree);
void TIMES_Set_All_Node_Priors_Top_Down(t_node *a, t_node *d, t_tree *tree);
void TIMES_Set_Floor(t_tree *tree);
void TIMES_Set_Floor_Post(t_node *a, t_node *d, t_tree *tree);
phydbl TIMES_Log_Conditional_Uniform_Density(t_tree *tree);
phydbl TIMES_Log_Yule(t_tree *tree);
void TIMES_Lk_Times_Trav(t_node *a, t_node *d, phydbl lim_inf, phydbl lim_sup, phydbl *logdens, t_tree *tree);
phydbl TIMES_Log_Number_Of_Ranked_Labelled_Histories(t_node *root, int per_slice, t_tree *tree);
void TIMES_Log_Number_Of_Ranked_Labelled_Histories_Post(t_node *a, t_node *d, int per_slice, phydbl *logn, t_tree *tree);
phydbl TIMES_Lk_Uniform_Core(t_tree *tree);
void TIMES_Get_Number_Of_Time_Slices(t_tree *tree);
void TIMES_Get_Number_Of_Time_Slices_Post(t_node *a, t_node *d, t_tree *tree);
void TIMES_Get_N_Slice_Spans(t_tree *tree);
void TIMES_Allocate_Vectors_Time_Slice_Combin(t_tree *tree);
void TIMES_Allocate_Vectors_Time_Slice_Combin_Post(t_node *a, t_node *d, t_tree *tree);
void TIMES_Update_Curr_Slice(t_tree *tree);
void TIMES_Lk_Uniform_Post(t_node *a, t_node *d, t_tree *tree);
void TIMES_Set_Root_Given_Tip_Dates(t_tree *tree);
void Get_Survival_Duration(t_tree *tree);
void Get_Survival_Duration_Post(t_node *a, t_node *d, t_tree *tree);
phydbl TIMES_Lk_Yule_Root_Marginal(t_tree *tree);
phydbl TIMES_Lk_Yule_Joint(t_tree *tree);
void TIMES_Update_Node_Ordering(t_tree *tree);
phydbl TIMES_Lk_Yule_Order(t_tree *tree);
void TIMES_Record_Prior_Times(t_tree *tree);
void TIMES_Reset_Prior_Times(t_tree *tree);
phydbl TIMES_Lk_Yule_Order_Root_Cond(t_tree *tree);
void TIMES_Connect_List_Of_Taxa(t_node **tax_list, int list_size, phydbl t_mrca, phydbl *times, int *nd_num, t_tree *mixt_tree);
void TIMES_Randomize_Tree_With_Time_Constraints(t_cal *cal_list, t_tree *tree);
phydbl TIMES_Lk_Birth_Death(int verbose, t_tree *tree);
int TIMES_Check_Node_Height_Ordering(t_tree *tree);
int TIMES_Check_Node_Height_Ordering_Post(t_node *a, t_node *d, t_tree *tree);
phydbl TIMES_Lk_Birth_Death_One_Node(phydbl t, phydbl min, phydbl max, t_tree *tree);
phydbl TIMES_Least_Square_Criterion(t_tree *tree);
void TIMES_Pre_Least_Square_Criterion(t_node *a, t_node *d, t_edge *b, phydbl *score, t_tree *tree);
void TIMES_Post_Randomize_Node_Ages(t_node *a, t_node *d, t_tree *tree);
void TIMES_Randomize_Node_Ages(t_tree *tree);
int TIMES_Calibrations_Apply_To_Tips_Only(t_tree *tree);
void TIMES_Randomize_Tip_Times_Given_Calibrations(t_tree *tree);
void TIMES_Time_To_Bl(t_tree *tree);
void TIMES_Time_To_Bl_Pre(t_node *a, t_node *d, t_edge *b, t_tree *tree);
void TIMES_Bl_To_Times_Post(t_node *a, t_node *d, t_edge *b, t_tree *tree);
void TIMES_Bl_To_Times(t_tree *tree);
phydbl TIMES_Tree_Length(t_tree *tree);
void TIMES_Copy_Time_Struct(t_time *from, t_time *to, int n_otu);
phydbl TIMES_Lk_Coalescent(t_tree *tree);
phydbl TIMES_Wrap_Lk_Coalescent(t_edge *b, t_tree *tree, supert_tree *stree);
void TIMES_Simulate_Coalescent(t_tree *tree);
phydbl TIMES_Lk_Coalescent_Range(t_dsk *young, t_dsk *old, t_tree *tree);
phydbl TIMES_Lk_SLFV(t_tree *tree);
phydbl TIMES_Lk_SLFV_Range(t_dsk *young, t_dsk *old, t_tree *tree);
phydbl TIMES_Lk_Range(t_dsk *young, t_dsk *old, t_tree *tree);
phydbl TIMES_Lk(t_tree *tree);
phydbl TIMES_Prior_Coalescent(t_tree *tree);
phydbl TIMES_Prior(t_tree *tree);
phydbl TIMES_Prior_SLFV(t_tree *tree);


#endif
