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

int PHYREX_Main(int argc, char **argv);
int PHYREX_Main_Simulate(int argc, char **argv);
int PHYREX_Main_Estimate(int argc, char **argv);
t_tree *PHYREX_Simulate(int n_otu, int n_sites, phydbl width, phydbl height, int r_seed);
phydbl PHYREX_Lk(t_tree *tree);
phydbl PHYREX_Wrap_Lk(t_edge *b, t_tree *tree, supert_tree *stree);
phydbl *PHYREX_MCMC(t_tree *tree);
int PHYREX_Is_In_Disk(t_geo_coord *coord, t_dsk *disk, t_phyrex_mod *mmod);
void PHYREX_New_Traj(t_dsk *start, t_dsk *end, t_tree *tree);
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
void PHYREX_Check_Struct(t_tree *tree);
void PHYREX_Store_Geo_Coord(t_geo_coord *t);
void PHYREX_Restore_Geo_Coord(t_geo_coord *t);
int PHYREX_Total_Number_Of_Intervals(t_tree *tree);
int PHYREX_Total_Number_Of_Coal_Disks(t_tree *tree);
int PHYREX_Total_Number_Of_Hit_Disks(t_tree *tree);
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
phydbl PHYREX_Simulate_Backward_Core(int new_loc, t_dsk *init_disk, t_tree *tree);
phydbl *PHYREX_Mean_Pairwise_Distance_Between_Lineage_Locations(t_tree *tree);
phydbl PHYREX_Random_Select_Time_Between_Jumps(t_tree *tree);
phydbl PHYREX_Simulate_Forward_Core(int n_sites, t_tree *tree);
int PHYREX_Is_In_Ldscape(t_ldsk *ldsk, t_phyrex_mod *mmod);
void PHYREX_Update_Lindisk_List_Core(t_dsk *disk, t_tree *tree);
phydbl PHYREX_Mean_Time_Between_Events(t_tree *tree);
void PHYREX_All_Pairs_Coal_Times_Dist(t_tree *tree);
void PHYREX_Rand_Pairs_Coal_Times_Dist(t_tree *tree);
phydbl PHYREX_Neighborhood_Size_Regression(t_tree *tree);
phydbl PHYREX_Neighborhood_Size(t_tree *tree);
phydbl PHYREX_Update_Radius(t_tree *tree);
phydbl PHYREX_Update_Sigsq(t_tree *tree);
void PHYREX_Read_Tip_Coordinates(t_ldsk **ldsk_a, t_tree *tree);
phydbl PHYREX_Sample_Rad_From_Prior(t_tree *tree);
void MCMC_PHYREX_Sigsq(t_tree *tree);
phydbl PHYREX_LnPrior_Sigsq(t_tree *tree);
phydbl PHYREX_Rate_Per_Unit_Area(t_tree *tree);
phydbl PHYREX_Tree_Height(t_tree *tree);
int PHYREX_Random_Insert_Ldsk_In_Next_List(t_ldsk *ins, t_ldsk *where);
void PHYREX_Insert_Ldsk_In_Next_List(t_ldsk *ins, int pos, t_ldsk *where);
t_ldsk *PHYREX_Remove_Path(t_ldsk *beg, t_ldsk *end, t_tree *tree);
void PHYREX_Insert_Path(t_ldsk *beg, t_ldsk *end, t_ldsk *path, t_tree *tree);
t_ldsk *PHYREX_Generate_Path(t_ldsk *beg, t_ldsk *end, phydbl n_evt, t_tree *tree);
phydbl PHYREX_Path_Logdensity(t_ldsk *beg, t_ldsk *end, phydbl cur_n_evt, t_tree *tree);
phydbl PHYREX_Time_Tree_Length(t_tree *tree);
void PHYREX_Time_Tree_Length_Pre(t_ldsk *a, t_ldsk *d, phydbl *len, t_tree *tree);
int PHYREX_Is_On_Path(t_ldsk *target, t_ldsk *beg, t_ldsk *end);
int PHYREX_Path_Len(t_ldsk *beg, t_ldsk *end);


#endif
