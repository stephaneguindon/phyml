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

int  TIMES_main(int argc, char **argv);
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
phydbl TIMES_Lk_Times(t_tree *tree);
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
void TIMES_Label_Edges_With_Calibration_Intervals(t_tree *tree);
void TIMES_Record_Prior_Times(t_tree *tree);
void TIMES_Reset_Prior_Times(t_tree *tree);
phydbl TIMES_Lk_Yule_Order_Root_Cond(t_tree *tree);

#endif
