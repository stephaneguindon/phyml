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
int MIGREP_Is_In_Disk(t_geo_coord *coord, t_dsk *disk);
void MIGREP_New_Traj(t_dsk *start, t_dsk *end, t_tree *tree);
int MIGREP_Draw(GtkWidget *widget, cairo_t *cr, gpointer *data);
void MIGREP_Remove_Disk(t_dsk *disk);
void MIGREP_Insert_Disk(t_dsk *disk);
t_ldsk *MIGREP_Prev_Coal_Lindisk(t_ldsk *t);
t_ldsk *MIGREP_Next_Coal_Lindisk(t_ldsk *t);
int MIGREP_Get_Next_Direction(t_ldsk *young, t_ldsk *old);
void MIGREP_One_New_Traj_Given_Disk(t_ldsk *s_ldsk, t_ldsk *e_ldsk);
void MIGREP_Update_Lindisk_List(phydbl time, t_ldsk **list, int *pos, t_dsk *disk);
void MIGREP_Update_Lindisk_List_Pre(t_ldsk *ldsk, phydbl time, t_ldsk **list, int *pos);
void MIGREP_Connect_Ldsk_Given_Disk(t_dsk **disk, int n_disk, t_ldsk *y_ldsk, t_ldsk *o_ldsk, int dir_o_y);
void MIGREP_Print_Struct(char sign, t_tree *tree);
phydbl MIGREP_Uniform_Path_Density(t_ldsk *y_ldsk, t_ldsk *o_ldsk);
void MIGREP_Check_Struct(t_tree *tree);
t_geo_coord *MIGREP_Copy_Geo_Coord(t_geo_coord *orig);
int MIGREP_Total_Number_Of_Intervals(t_tree *tree);
int MIGREP_Total_Number_Of_Coal_Disks(t_tree *tree);
int MIGREP_Total_Number_Of_Hit_Disks(t_tree *tree);
phydbl MIGREP_Log_Uniform_Rectangle_Overlap(t_dsk *disk, t_migrep_mod *mmod);
void MIGREP_Runif_Rectangle_Overlap(t_ldsk *ldsk, t_dsk *disk, t_migrep_mod *mod);
void MIGREP_One_New_Traj(t_ldsk *y_ldsk, t_ldsk *o_ldsk, int dir_o_y, t_dsk *xtra_dsk, int *n_add_disk, t_tree *tree);

#endif
