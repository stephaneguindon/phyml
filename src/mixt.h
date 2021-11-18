/*

PHYML :  a program that  computes maximum likelihood  phylogenies from
DNA or AA homologous sequences 

Copyright (C) Stephane Guindon. Oct 2003 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/

#include <config.h>

#ifndef MIXT_H
#define MIXT_H

#include "utilities.h"

void MIXT_Connect_Edges_To_Next_Prev_Child_Parent(t_tree *tree);
void MIXT_Connect_Nodes_To_Next_Prev_Child_Parent(t_tree *tree);
void MIXT_Connect_Sprs_To_Next_Prev_Child_Parent(t_tree *tree);
void MIXT_Turn_Branches_OnOff(int onoff,t_tree *tree);
phydbl *MIXT_Get_Lengths_Of_This_Edge(t_edge *mixt_b, t_tree *tree);
void MIXT_Set_Lengths_Of_This_Edge(phydbl *lens,t_edge *mixt_b, t_tree *tree);
void MIXT_Post_Order_Lk(t_node *mixt_a,t_node *mixt_d,t_tree *mixt_tree);
void MIXT_Pre_Order_Lk(t_node *mixt_a,t_node *mixt_d,t_tree *mixt_tree);
phydbl MIXT_Lk(t_edge *mixt_b,t_tree *mixt_tree);
void MIXT_Update_Partial_Lk(t_tree *mixt_tree,t_edge *mixt_b,t_node *mixt_d);
void MIXT_Update_PMat_At_Given_Edge(t_edge *mixt_b,t_tree *mixt_tree);
int *MIXT_Get_Number_Of_Classes_In_All_Mixtures(t_tree *mixt_tree);
t_tree **MIXT_Record_All_Mixtures(t_tree *mixt_tree);
void MIXT_Break_All_Mixtures(int *c_max,t_tree *mixt_tree);
void MIXT_Reconnect_All_Mixtures(t_tree **tree_list,t_tree *mixt_tree);
int *MIXT_Record_Has_Invariants(t_tree *mixt_tree);
void MIXT_Reset_Has_Invariants(int *has_invariants,t_tree *mixt_tree);
void MIXT_Check_Invar_Setup(t_tree *mixt_tree);
void MIXT_Prune_Subtree(t_node *mixt_a,t_node *mixt_d,t_edge **mixt_target,t_edge **mixt_residual,t_tree *mixt_tree);
void MIXT_Graft_Subtree(t_edge *mixt_target, t_node *mixt_link, t_node *mixt_link_daughter, t_edge *mixt_residual, t_node *mixt_target_nd, t_tree *mixt_tree);
void MIXT_Br_Len_Opt(t_edge *mixt_b, t_tree *mixt_tree);
void MIXT_Check_Number_Of_Invar_Classes(t_tree *mixt_tree);
void MIXT_Prepare_Tree_For_Lk(t_tree *tree);
void MIXT_Check_Invar_Struct_In_Each_Partition_Elem(t_tree *mixt_tree);
void MIXT_Check_RAS_Struct_In_Each_Partition_Elem(t_tree *mixt_tree);
void MIXT_Br_Len_Involving_Invar(t_tree *mixt_tree);
void MIXT_Br_Len_Not_Involving_Invar(t_tree *mixt_tree);
phydbl MIXT_Unscale_Br_Len_Multiplier_Tree(t_tree *mixt_tree);
phydbl MIXT_Rescale_Br_Len_Multiplier_Tree(t_tree *mixt_tree);
void MIXT_Set_Alias_Subpatt(int onoff, t_tree *mixt_tree);
void MIXT_Check_Single_Edge_Lens(t_tree *mixt_tree);
void MIXT_Update_Eigen(t_mod *mixt_mod);
int MIXT_Pars(t_edge *mixt_b, t_tree *mixt_tree);
void MIXT_Set_Pars_Thresh(t_tree *mixt_tree);
void MIXT_Bootstrap(char *best_tree, xml_node *root);
void MIXT_Chain_All(t_tree *mixt_tree);
void MIXT_Chain_String(t_string *curr, t_string *next);
void MIXT_Chain_Scalar_Dbl(scalar_dbl *curr, scalar_dbl *next);
void MIXT_Chain_Rmat(t_rmat *curr, t_rmat *next);
void MIXT_Chain_Rmat(t_rmat *curr, t_rmat *next);
void MIXT_Chain_Efrq(t_efrq *curr, t_efrq *next);
void MIXT_Chain_RAS(t_ras *curr, t_ras *next);
void MIXT_Chain_Eigen(eigen *curr, eigen *next);
void MIXT_Chain_Vector_Dbl(vect_dbl *curr, vect_dbl *next);
void MIXT_Chain_Sprs(t_tree *tree);
void MIXT_Chain_Nodes(t_tree *tree);
void MIXT_Chain_Edges(t_tree *tree);
phydbl MIXT_Get_Mean_Edge_Len(t_edge *mixt_b, t_tree *tree);
phydbl MIXT_Get_Sum_Chained_Scalar_Dbl(scalar_dbl *s);
phydbl MIXT_Get_Sum_Of_Probas_Across_Mixtures(phydbl r_mat_weight_sum, phydbl e_frq_weight_sum, t_tree *mixt_tree);
phydbl MIXT_Rescale_Free_Rate_Tree(t_tree *mixt_tree);
void MIXT_Set_Br_Len_Var(t_edge *mixt_b, t_tree *mixt_tree);
void MIXT_Optimize_Br_Len_Multiplier(t_tree *mixt_tree);
void MIXT_Update_Br_Len_Multipliers(t_mod *mod);
void MIXT_Init_Model(t_mod *mod);
t_tree *MIXT_Starting_Tree(t_tree *mixt_tree);
void MIXT_Connect_Cseqs_To_Nodes(t_tree *mixt_tree);
void MIXT_Check_Edge_Lens_In_One_Elem(t_tree *mixt_tree);
void MIXT_Check_Edge_Lens_In_All_Elem(t_tree *mixt_tree);
void MIXT_Turn_Branches_OnOff_In_One_Elem(int onoff, t_tree *mixt_tree);
void MIXT_Turn_Branches_OnOff_In_All_Elem(int onoff, t_tree *mixt_tree);
void MIXT_Init_T_Beg(t_tree *mixt_tree);
void MIXT_Prepare_All(int num_rand_tree, t_tree *mixt_tree);
void MIXT_Init_T_End(t_tree *mixt_tree);
void MIXT_Add_Root(t_edge *mixt_b, t_tree *mixt_tree);
void MIXT_Check_Model_Validity(t_tree *mixt_tree);
void MIXT_Ancestral_Sequences_One_Node(t_node *mixt_d, t_tree *mixt_tree, int print);
void MIXT_Update_Partial_Pars(t_tree *mixt_tree, t_edge *mixt_b, t_node *mixt_d);
void MIXT_Chain_Cal(t_tree *mixt_tree);
void MIXT_Chain_Rates(t_rate *curr, t_rate *next);
void MIXT_RATES_Update_Edge_Lengths(t_tree *mixt_tree);
phydbl MIXT_dLk(phydbl *l, t_edge *mixt_b, t_tree *mixt_tree);
void MIXT_Update_Eigen_Lr(t_edge *mixt_b, t_tree *mixt_tree);
int MIXT_Part_Mixt_Size(t_tree *mixt_tree);
int MIXT_Mixt_Size(t_tree *mixt_tree);
void MIXT_Set_Use_Eigen_Lr(int yn, t_tree *mixt_tree);
void MIXT_Set_Update_Eigen_Lr(int yn, t_tree *mixt_tree);
void MIXT_Backup_Partial_Pars(t_node *mixt_d, t_edge *mixt_b, t_tree *mixt_tree);
void MIXT_Restore_Partial_Pars(t_node *mixt_d, t_edge *mixt_b, t_tree *mixt_tree);
void MIXT_Backup_Partial_Lk(t_node *mixt_d, t_edge *mixt_b, t_tree *mixt_tree);
void MIXT_Restore_Partial_Lk(t_node *mixt_d, t_edge *mixt_b, t_tree *mixt_tree);
void MIXT_Backup_Partial_Scale(t_node *mixt_d, t_edge *mixt_b, t_tree *mixt_tree);
void MIXT_Restore_Partial_Scale(t_node *mixt_d, t_edge *mixt_b, t_tree *mixt_tree);
void MIXT_Set_Both_Sides(int yesno, t_tree *mixt_tree);
void MIXT_Set_Br_Len(phydbl val, t_edge *mixt_b, t_tree *mixt_tree);
void MIXT_Multiply_Scalar_Dbl(scalar_dbl *this, phydbl scalar);
void MIXT_Sample_Ancestral_Seq(int mutmap, int fromprior, t_tree *mixt_tree);
void MIXT_Set_Update_Eigen(int yn, t_mod *mixt_mod);
void MIXT_Set_Ignore_Root(int yesno, t_tree *mixt_tree);
void MIXT_Make_Spr(t_tree *mixt_tree);
void MIXT_Make_Tree_For_Lk(t_tree *mixt_tree);
void MIXT_Make_Tree_For_Pars(t_tree *mixt_tree);
void MIXT_Connect_Tip_Disks(t_tree *mixt_tree);
void MIXT_Propagate_Tree_Update(t_tree *mixt_tree);
void MIXT_Set_Bl_From_Rt(int yn, t_tree *mixt_tree);
void MIXT_Copy_Tree(t_tree *ori, t_tree *cpy);
void MIXT_Init_NNI_Score(phydbl val, t_edge *mixt_b, t_tree *mixt_tree);
t_tree *MIXT_Duplicate_Tree(t_tree *ori);
void MIXT_Set_Model_Parameters(t_mod *mixt_mod);
void MIXT_Print_Site_Lk(t_tree *mixt_tree, FILE *fp);
void MIXT_Exchange_Nodes(t_node *a, t_node *d, t_node *w, t_node *v, t_tree *mixt_tree);
void MIXT_Chain_Models(t_tree *mixt_tree);
void MIXT_Repeat_Task(void (*Task_Function)(),t_tree *mixt_tree);
void MIXT_Free_Tree(t_tree *mixt_tree);

#endif
