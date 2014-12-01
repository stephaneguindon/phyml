/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include <config.h>

#ifndef MIGREP_H
#define MIGREP_H

#include <gtk/gtk.h>
#include <glib.h>
#include <pthread.h>
#include "utilities.h"

t_tree *MIGREP_Simulate_Backward(int n_otu, phydbl width, phydbl height);
int MIGREP_Main(int argc, char **argv);
phydbl MIGREP_Lk(t_tree *tree);
phydbl MIGREP_Wrap_Lk(t_edge *b, t_tree *tree, supert_tree *stree);
void MIGREP_MCMC(t_tree *tree);
int MIGREP_Is_In_Disk(t_geo_coord *coord, t_dsk *disk, t_migrep_mod *mmod);
void MIGREP_New_Traj(t_dsk *start, t_dsk *end, t_tree *tree);
int MIGREP_Draw(GtkWidget *widget, cairo_t *cr, gpointer *data);
void MIGREP_Remove_Disk(t_dsk *disk);
void MIGREP_Insert_Disk(t_dsk *disk);
t_ldsk *MIGREP_Prev_Coal_Lindisk(t_ldsk *t);
t_ldsk *MIGREP_Next_Coal_Lindisk(t_ldsk *t);
int MIGREP_Get_Next_Direction(t_ldsk *young, t_ldsk *old);
void MIGREP_Update_Lindisk_List(t_tree *tree);
void MIGREP_Update_Lindisk_List_Pre(t_dsk *disk);
/* void MIGREP_Update_Lindisk_List(phydbl time, t_ldsk **list, int *pos, t_dsk *disk); */
/* void MIGREP_Update_Lindisk_List_Pre(t_ldsk *ldsk, phydbl time, t_ldsk **list, int *pos); */
void MIGREP_Connect_Ldsk_Given_Disk(t_dsk **disk, int n_disk, t_ldsk *y_ldsk, t_ldsk *o_ldsk, int dir_o_y);
void MIGREP_Print_Struct(char sign, t_tree *tree);
phydbl MIGREP_Uniform_Path_Density(t_ldsk *y_ldsk, t_ldsk *o_ldsk, t_tree *tree);
void MIGREP_Check_Struct(t_tree *tree);
void MIGREP_Store_Geo_Coord(t_geo_coord *t);
void MIGREP_Restore_Geo_Coord(t_geo_coord *t);
int MIGREP_Total_Number_Of_Intervals(t_tree *tree);
int MIGREP_Total_Number_Of_Coal_Disks(t_tree *tree);
int MIGREP_Total_Number_Of_Hit_Disks(t_tree *tree);
phydbl MIGREP_Log_Dunif_Rectangle_Overlap(t_ldsk *ldsk, t_dsk *disk, t_migrep_mod *mmod);
phydbl MIGREP_Runif_Rectangle_Overlap(t_ldsk *ldsk, t_dsk *disk, t_migrep_mod *mod);
int MIGREP_One_New_Traj(t_ldsk *y_ldsk, t_ldsk *o_ldsk, int dir_o_y, t_dsk *xtra_dsk, int n_cur_disk, t_tree *tree);
phydbl MIGREP_Wrap_Prior_Radius(t_edge *e, t_tree *tree, supert_tree *st);
phydbl MIGREP_LnPrior_Radius(t_tree *tree);
void MIGREP_Initial_Ldsk_Pos(t_tree *tree);
phydbl MIGREP_Min_Radius(t_tree *tree);
void MIGREP_Get_Min_Max_Disk_Given_Ldsk(t_dsk *disk, phydbl **min, phydbl **max, t_tree *tree);
void MIGREP_Get_Min_Max_Ldsk_Given_Disk(t_ldsk *ldsk, phydbl **min, phydbl **max, t_tree *tree);
void MIGREP_One_New_Traj_Given_Disk(t_ldsk *y_ldsk, t_ldsk *o_ldsk, t_tree *tree);
void MIGREP_Update_Disk_Ldsk_Subtree(t_ldsk *root_ldsk, t_tree *tree);
void MIGREP_Update_Disk_Ldsk_Subtree_Pre(t_ldsk *old_ldsk, t_ldsk *young_ldsk, t_ldsk *root_ldsk, t_tree *tree);
void MIGREP_Restore_Disk_Ldsk_Subtree(t_ldsk *root_ldsk, t_tree *tree);
void MIGREP_Restore_Disk_Ldsk_Subtree_Pre(t_ldsk *old_ldsk, t_ldsk *young_ldsk, t_tree *tree);
void MIGREP_Proposal_Disk_Ldsk_Subtree(t_ldsk *root_ldsk, phydbl *logdens, t_tree *tree);
void MIGREP_Proposal_Disk_Ldsk_Subtree_Pre(t_ldsk *old_ldsk, t_ldsk *young_ldsk, t_ldsk *root_ldsk, phydbl *logdens, t_tree *tree);
phydbl MIGREP_LnPrior_Lbda(t_tree *tree);
phydbl MIGREP_LnPrior_Mu(t_tree *tree);
void MIGREP_Ldsk_To_Tree(t_tree *tree);
void MIGREP_Ldsk_To_Tree_Post(t_node *a, t_ldsk *ldsk, int *available, t_tree *tree);
phydbl MIGREP_Rnorm_Trunc(t_ldsk *ldsk, t_dsk *disk, t_migrep_mod *mod);
void MIGREP_Remove_Lindisk_Next(t_ldsk *ldsk, t_ldsk *rm);
void MIGREP_Simulate_Backward_Core(t_tree *tree);

#endif
