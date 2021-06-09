/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include <config.h>

#ifndef RRW_H
#define RRW_H

#include "utilities.h"

phydbl RRW_Lk(t_tree *tree);
phydbl RRW_Forward_Lk_Path(t_ldsk *a, t_ldsk *d, t_tree *tree);
void RRW_Rescale_Times(int prod, t_tree *tree);
void RRW_Rescale_Times_Pre(t_node *a, t_node *d, phydbl cur_ta, int prod, t_tree *tree);
phydbl RRW_Independent_Contrasts(t_tree *tree);
phydbl RRW_Lk_Range(t_dsk *young, t_dsk *old, t_tree *tree);
phydbl RRW_Forward_Lk_Range(t_dsk *young, t_dsk *old, t_tree *tree);
phydbl RRW_Lk_Core(t_dsk *disk, t_tree *tree);
void RRW_Generate_Ldsk_New_Location(t_ldsk *l, phydbl rad,  int dim_idx, int child_only, t_tree *tree);
phydbl RRW_Density_Ldsk_Location(t_ldsk *l, phydbl rad, int dim_idx, int child_only, t_tree *tree);
short int RRW_Is_Rw(t_phyrex_mod *mod);
phydbl RRW_Prior_Sigsq_Scale(t_tree *tree);
void RRW_Update_Normalization_Factor(t_tree *tree);
phydbl RRW_Mean_Displacement_Rate(t_tree *tree);
phydbl RRW_Prior_Sigsq(t_tree *tree);
phydbl RRW_Prior(t_tree *tree);
void RRW_Generate(t_tree *tree);


#endif
