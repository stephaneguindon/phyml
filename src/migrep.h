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
phydbl MIGREP_Lk(t_disk_evt *devt, t_migrep_mod *mmod);
phydbl MIGREP_Wrap_Lk(t_edge *b, t_tree *tree, supert_tree *stree);
void MIGREP_MCMC(t_tree *tree);
int MIGREP_Is_In_Disk(t_geo_coord *coord, t_disk_evt *devt);
void MIGREP_New_Traj(t_disk_evt *start, t_disk_evt *end, t_tree *tree);
void MIGREP_Copy_Coord(t_geo_coord *ori, t_geo_coord *cpy);
int MIGREP_Draw(GtkWidget *widget, cairo_t *cr, gpointer *data);
void MIGREP_Remove_Devt(t_disk_evt *devt);
void MIGREP_Insert_Devt(t_disk_evt *devt);
t_lindisk_nd *MIGREP_Prev_Coal_Lindisk(t_lindisk_nd *t);
t_lindisk_nd *MIGREP_Next_Coal_Lindisk(t_lindisk_nd *t);
int MIGREP_Get_Next_Direction(t_lindisk_nd *young, t_lindisk_nd *old);
void MIGREP_One_New_Traj(t_lindisk_nd *s_ldsk, t_lindisk_nd *e_ldsk);

#endif
