/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include <config.h>

#ifndef PHYREX_H
#define PHYREX_H

#include "utilities.h"

void PHYREX_XML(char *xml_filename);
int PHYREX_Main(int argc, char **argv);
int PHYREX_Main_Simulate(int argc, char **argv);
phydbl PHYREX_Lk(t_tree *tree);
phydbl PHYREX_Lk_RW(t_tree *tree);
phydbl PHYREX_Wrap_Lk(t_edge *b, t_tree *tree, supert_tree *stree);
phydbl *PHYREX_MCMC(t_tree *tree);
int PHYREX_Is_In_Disk(t_geo_coord *coord, t_dsk *disk, t_phyrex_mod *mmod);
void PHYREX_Remove_Disk(t_dsk *disk);
void PHYREX_Insert_Disk(t_dsk *ins, t_tree *tree);
t_ldsk *PHYREX_Prev_Coal_Lindisk(t_ldsk *t);
t_ldsk *PHYREX_Next_Coal_Lindisk(t_ldsk *t);
int PHYREX_Get_Next_Direction(t_ldsk *young, t_ldsk *old);
void PHYREX_Update_Lindisk_List(t_tree *tree);
void PHYREX_Update_Lindisk_List_Pre(t_dsk *disk, t_tree *tree);
/* void PHYREX_Update_Lindisk_List(phydbl time, t_ldsk **list, int *pos, t_dsk *disk); */
/* void PHYREX_Update_Lindisk_List_Pre(t_ldsk *ldsk, phydbl time, t_ldsk **list, int *pos); */
void PHYREX_Connect_Ldsk_Given_Disk(t_dsk **disk, int n_disk, t_ldsk *y_ldsk, t_ldsk *o_ldsk, int dir_o_y);
void PHYREX_Print_Struct(char sign, t_tree *tree);
phydbl PHYREX_Uniform_Path_Density(t_ldsk *y_ldsk, t_ldsk *o_ldsk, t_tree *tree);
int PHYREX_Check_Struct(t_tree *tree, int exit);
void PHYREX_Store_Geo_Coord(t_geo_coord *t);
void PHYREX_Restore_Geo_Coord(t_geo_coord *t);
int PHYREX_Total_Number_Of_Intervals(t_tree *tree);
int PHYREX_Total_Number_Of_Coal_Disks(t_tree *tree);
int PHYREX_Total_Number_Of_Hit_Disks(t_tree *tree);
int PHYREX_Total_Number_Of_Floating_Disks(t_tree *tree);
phydbl PHYREX_Log_Dunif_Rectangle_Overlap(t_ldsk *ldsk, t_dsk *disk, t_phyrex_mod *mmod);
phydbl PHYREX_Runif_Rectangle_Overlap(t_ldsk *ldsk, t_dsk *disk, t_phyrex_mod *mod);
int PHYREX_One_New_Traj(t_ldsk *y_ldsk, t_ldsk *o_ldsk, int dir_o_y, t_dsk *xtra_dsk, int n_cur_disk, t_tree *tree);
phydbl PHYREX_Wrap_Prior_Radius(t_edge *e, t_tree *tree, supert_tree *st);
phydbl PHYREX_LnPrior_Radius(t_tree *tree);
void PHYREX_Initial_Ldsk_Pos(t_tree *tree);
phydbl PHYREX_Min_Radius(t_tree *tree);
void PHYREX_Get_Min_Max_Disk_Given_Ldsk(t_dsk *disk, phydbl **min, phydbl **max, t_tree *tree);
void PHYREX_Get_Min_Max_Ldsk_Given_Disk(t_ldsk *ldsk, phydbl **min, phydbl **max, t_tree *tree);
void PHYREX_One_New_Traj_Given_Disk(t_ldsk *y_ldsk, t_ldsk *o_ldsk, t_tree *tree);
void PHYREX_Update_Disk_Ldsk_Subtree(t_ldsk *root_ldsk, t_tree *tree);
void PHYREX_Update_Disk_Ldsk_Subtree_Pre(t_ldsk *old_ldsk, t_ldsk *young_ldsk, t_ldsk *root_ldsk, t_tree *tree);
void PHYREX_Restore_Disk_Ldsk_Subtree(t_ldsk *root_ldsk, t_tree *tree);
void PHYREX_Restore_Disk_Ldsk_Subtree_Pre(t_ldsk *old_ldsk, t_ldsk *young_ldsk, t_tree *tree);
void PHYREX_Proposal_Disk_Ldsk_Subtree(t_ldsk *root_ldsk, phydbl *logdens, t_tree *tree);
void PHYREX_Proposal_Disk_Ldsk_Subtree_Pre(t_ldsk *old_ldsk, t_ldsk *young_ldsk, t_ldsk *root_ldsk, phydbl *logdens, t_tree *tree);
phydbl PHYREX_LnPrior_Lbda(t_tree *tree);
phydbl PHYREX_LnPrior_Mu(t_tree *tree);
void PHYREX_Ldsk_To_Tree(t_tree *tree);
void PHYREX_Ldsk_To_Tree_Post(t_node *a, t_ldsk *ldsk, int *available, t_tree *tree);
phydbl PHYREX_Rnorm_Trunc(t_ldsk *ldsk, t_dsk *disk, t_phyrex_mod *mod);
void PHYREX_Remove_Lindisk_Next(t_ldsk *ldsk, t_ldsk *rm);
phydbl *PHYREX_Mean_Pairwise_Distance_Between_Lineage_Locations(t_tree *tree);
phydbl PHYREX_Random_Select_Time_Between_Jumps(t_tree *tree);
int PHYREX_Is_In_Ldscape(t_ldsk *ldsk, t_phyrex_mod *mmod);
void PHYREX_Update_Lindisk_List_Core(t_dsk *disk, t_tree *tree);
phydbl PHYREX_Mean_Time_Between_Events(t_tree *tree);
void PHYREX_All_Pairs_Coal_Times_Dist(t_tree *tree);
void PHYREX_Rand_Pairs_Coal_Times_Dist(t_tree *tree);
void PHYREX_Read_Tip_Coordinates(t_tree *tree);
phydbl PHYREX_LnPrior_Sigsq(t_tree *tree);
phydbl PHYREX_Tree_Height(t_tree *tree);
int PHYREX_Random_Insert_Ldsk_In_Next_List(t_ldsk *ins, t_ldsk *where);
void PHYREX_Insert_Ldsk_In_Next_List(t_ldsk *ins, int pos, t_ldsk *where);
t_ldsk *PHYREX_Remove_Path(t_ldsk *beg, t_ldsk *end, int *pos_end, t_tree *tree);
void PHYREX_Insert_Path(t_ldsk *beg, t_ldsk *end, t_ldsk *path, int pos, t_tree *tree);
phydbl PHYREX_Time_Tree_Length(t_tree *tree);
void PHYREX_Time_Tree_Length_Pre(t_ldsk *a, t_ldsk *d, phydbl *len, t_tree *tree);
int PHYREX_Is_On_Path(t_ldsk *target, t_ldsk *beg, t_ldsk *end);
int PHYREX_Path_Len(t_ldsk *beg, t_ldsk *end);
phydbl PHYREX_Lk_Core(t_dsk *disk, t_tree *tree);
phydbl PHYREX_Lk_RW_Core(t_dsk *disk, t_tree *tree);
void PHYREX_Print_Disk_Lk(t_tree *tree);
t_ldsk *PHYREX_Find_Lca_Pair_Of_Ldsk(t_ldsk *n1, t_ldsk *n2, t_tree *tree);
void PHYREX_Get_List_Of_Ancestors(t_ldsk *start, t_ldsk ***list, int *len, t_tree *tree);
phydbl PHYREX_Dist_To_Lca(t_ldsk *d, t_ldsk *lca);
phydbl PHYREX_Dist_Between_Two_Ldsk(t_ldsk *n1,  t_ldsk *n2, t_tree *tree);
phydbl PHYREX_Lk_Range(t_dsk *young, t_dsk *old, t_tree *tree);
int PHYREX_Number_Of_Sampled_Demes(t_tree *tree);
void PHYREX_Tree_To_Ldsk(t_tree *tree);
void PHYREX_Tree_To_Ldsk_Post(t_node *a, t_node *d, t_dsk *a_disk, t_tree *tree);
void PHYREX_Make_And_Connect_Tip_Disks(t_tree *tree);
void PHYREX_Simulate_Disk_And_Node_Times(t_tree *tree);
int PHYREX_Number_Of_Outgoing_Ldsks(t_dsk *disk);
t_ldsk *PHYREX_Random_Select_Outgoing_Ldsk(t_dsk *disk);
void PHYREX_Strip_And_Reconnect_Tree(t_tree *tree);
int PHYREX_Scale_All(phydbl scale, t_dsk *start_disk, t_tree *tree);
void  PHYREX_Update_Lindisk_List_Range(t_dsk *young, t_dsk *old, t_tree *tree);
phydbl PHYREX_Lk_Core_Range(t_dsk *young, t_dsk *old, t_tree *tree);
int PHYREX_Number_Of_Intervals_Range(t_dsk *young, t_dsk *old, t_tree *tree);
void PHYREX_Oldest_Sampled_Disk(t_tree *tree);
phydbl PHYREX_Time_Of_Descendants(t_ldsk *ldsk, t_tree *tree);
phydbl PHYREX_Time_Of_Prev_Sampled_Disk(t_dsk *disk, t_tree *tree);
phydbl PHYREX_Time_Of_Next_Sampled_Disk(t_dsk *disk, t_tree *tree);
phydbl PHYREX_Lk_Time_Component(t_tree *tree);
phydbl PHYREX_Lk_Space_Component(t_tree *tree);
t_ldsk *PHYREX_Find_Ldsk_From_Id(char *id, t_ldsk *root);
t_geo_coord *PHYREX_Mean_Next_Loc(t_ldsk *ldsk, t_tree *tree);
phydbl PHYREX_Root_To_Tip_Realized_Sigsq(t_tree *tree);
phydbl PHYREX_Tip_To_Root_Realized_Sigsq(t_tree *tree);
phydbl PHYREX_Realized_Dispersal_Dist(t_tree *tree);
phydbl PHYREX_Tip_To_Root_Realized_Bis_Sigsq(t_tree *tree);
void PHYREX_Label_Nodes_With_Locations(t_tree *tree);
void PHYREX_Label_Edges(t_tree *tree);
phydbl PHYREX_Tip_To_Root_Realized_Ter_Sigsq(t_tree *tree);
void PHYREX_Remove_All_Disks_Except_Coal_And_Tips(t_tree *tree);
phydbl PHYREX_Path_Logdensity(t_ldsk *young, t_ldsk *old, phydbl *sd, t_tree *tree);
void PHYREX_Sample_Path(t_ldsk *young, t_ldsk *old, phydbl *sd, phydbl *global_hr, t_tree *tree);
t_ldsk *PHYREX_Generate_Path(t_ldsk *young, t_ldsk *old, phydbl n_evt, phydbl *sd, t_tree *tree);
void PHYREX_Simulate_Backward_Core(t_dsk *disk,int avoid_multiple_mergers, t_tree *tree);
t_tree *PHYREX_Simulate(int n_otu, int n_sites, phydbl w, phydbl h, phydbl  lbda, phydbl rad, phydbl mu, int r_seed, int modid);
phydbl PHYREX_Update_Sigsq(t_tree *tree);
void PHYREX_Get_Baselines(t_tree *tree);
void PHYREX_Get_Baselines_Post(t_ldsk *ldsk, t_tree *tree);
int PHYREX_Total_Number_Of_Single_Hit_Disks(t_tree *tree);
t_edge *PHYREX_Edge_Between_Two_Ldsks(t_ldsk *a, t_ldsk *d, t_tree *tree);
void PHYREX_Update_Node_Times_Given_Disks(t_tree *tree);
void PHYREX_Update_Ldsk_Rates_Given_Edges(t_tree *tree);
void PHYREX_Update_Edge_Rates_Given_Ldsks(t_tree *tree);
void PHYREX_Update_Ldsk_Rates_Given_One_Edge(t_node *d, t_tree *tree);
void PHYREX_Update_Edge_Sigsq_Given_Ldsks(t_tree *tree);
void PHYREX_Update_Ldsk_Sigsq_Given_One_Edge(t_node *d, t_tree *tree);
void PHYREX_Update_Ldsk_Sigsq_Given_Edges(t_tree *tree);
void PHYREX_Duplicate_Ldsk_Struct(t_tree *from, t_tree *where);
phydbl PHYREX_Get_Posterior(t_tree *tree);
phydbl PHYREX_Realized_Dispersal_Dist_Alt(t_tree *tree);
void PHYREX_Evolve_All(t_tree *tree);
void PHYREX_Insert_Disk_At(t_dsk *ins, t_dsk *disk);
void PHYREX_Move_Disk_Updown(t_dsk *this, phydbl target_time, t_tree *tree);
void PHYREX_Restore_Disk_Times(t_tree *tree);
void PHYREX_Record_Disk_Times(t_tree *tree);
t_dsk *PHYREX_Next_Floating_Disk(t_dsk *disk);
void PHYREX_Swap_Coords(t_ldsk *a, t_ldsk *b, t_tree *tree);
void PHYREX_Print_Disk(t_tree *tree);
void PHYREX_Check_Disk_Times(t_tree *tree);
void PHYREX_Exchange_Ldsk(t_ldsk *a, t_ldsk *d, t_ldsk *w, t_ldsk *v, int aw, int dv, t_tree *tree);

#endif
