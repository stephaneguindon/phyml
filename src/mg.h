#include <config.h>

#ifndef PART_H
#define PART_H

#include "utilities.h"

void Menu_Supertree(option *input);
void PART_Print_Nodes(t_node *a, t_node *d, supert_tree *st);
supert_tree *PART_Make_Supert_tree_Light(option *input);
void PART_Make_Supert_tree_Full(supert_tree *st, option *input, calign **data);
void PART_Get_List_Of_Reachable_Tips(t_node *a, t_node *d, calign **data, supert_tree *st);
void PART_Get_List_Of_Reachable_Tips_Pre(t_node *a, t_node *d, supert_tree *st);
void PART_Get_List_Of_Reachable_Tips_Post(t_node *a, t_node *d, supert_tree *st);
void PART_Prune_St_Topo(t_tree *tree, calign *data, supert_tree *st);
void PART_Match_St_Nodes_In_Gt_Recurr(t_node *a_gt, t_node *d_gt, t_node *a_st, t_node *d_st, t_tree *gt, supert_tree *st);
void PART_Match_St_Nodes_In_Gt(t_tree *tree, supert_tree *st);
void PART_Match_St_Edges_In_Gt(t_tree *gt, supert_tree *st);
void PART_Match_St_Edges_In_Gt_Recurr(t_node *a, t_node *d, t_node *a_st, t_node *d_st, t_tree *gt, supert_tree *st);
void PART_Simu(supert_tree *tr);
int PART_Mov_Backward_Topo_Bl(supert_tree *st, phydbl lk_old, t_edge **tested_b, int n_tested);
int PART_Get_Species_Found_In_St(supert_tree *st, calign *data);
void PART_Map_St_Nodes_In_Gt_Pre(t_node *a_st, t_node *d_st, t_tree *gt, supert_tree *st);
void PART_Map_St_Nodes_In_Gt_Post(t_node *a_st, t_node *d_st, t_tree *gt, supert_tree *st);
void PART_Map_St_Nodes_In_Gt(t_tree *gt, supert_tree *st);
void PART_Map_St_Nodes_In_Gt_One_Edge(t_node *a_st, t_node *d_st, t_edge *b_st, t_tree *gt, supert_tree *st);
void PART_Map_St_Edges_In_Gt(t_tree *gt, supert_tree *st);
phydbl PART_Lk(supert_tree *st);
int PART_Pars(supert_tree *st);
int PART_Spr(phydbl init_lnL, supert_tree *st);
void PART_Speed_Spr(supert_tree *st);
int Map_Spr_Move(t_edge *st_pruned, t_edge *st_target, t_node *st_link, t_tree *gt, supert_tree *st);
void PART_Test_All_Spr_Targets(t_edge *pruned, t_node *n_link, supert_tree *st);
void PART_Test_One_Spr_Target_Recur(t_node *a, t_node *d, t_edge *target, t_edge *pruned, t_node *n_link, supert_tree *st);
void PART_Test_One_Spr_Target(t_edge *st_p, t_edge *st_t, t_node *n_link, supert_tree *st);
int PART_Test_List_Of_Regraft_Pos(t_spr **st_spr_list, int list_size, supert_tree *st);
int PART_Try_One_Spr_Move(t_spr *st_move, supert_tree *st);
void PART_Map_Gt_Edges_In_St(t_tree *gt, supert_tree *st);
void PART_NNI(t_edge *st_b, supert_tree *st);
void PART_Swap(t_node *st_a, t_node *st_b, t_node *st_c, t_node *st_d, supert_tree *st);
void PART_Set_Bl(phydbl **bl, supert_tree *st);
void PART_Restore_Br_Len(supert_tree *st);
void PART_Record_Br_Len(supert_tree *st);
phydbl PART_Lk_At_Given_Edge(t_edge *st_b, supert_tree *st);
phydbl PART_Update_Lk_At_Given_Edge(t_edge *st_b, supert_tree *st);
void PART_Fill_Model_Partitions_Table(supert_tree *st);
phydbl PART_Br_Len_Brent(t_edge *st_b, int quickdirty, supert_tree *tree);
void PART_Initialise_Bl_Partition(supert_tree *st);
void PART_Update_P_Lk(t_edge *st_b, t_node *st_n, supert_tree *st);
void PART_Optimize_Br_Len_Serie(t_node *st_a, t_node *st_d, t_edge *st_b, supert_tree *st);
void PART_Update_Bl_Swaped(t_edge **st_b, int n, supert_tree *st);
void PART_Update_Bl(phydbl fact, supert_tree *st);
void PART_Make_N_Swap(t_edge **st_b, int beg, int end, supert_tree *st);
void PART_Do_Mapping(supert_tree *st);
void PART_Update_PMat(t_edge *st_b, supert_tree *st);
void PART_Print_Bl(supert_tree *st);
void PART_Check_Extra_Taxa(supert_tree *st);
int PART_main(int argc, char **argv);
#endif

