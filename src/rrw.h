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

phydbl RRW_Forward_Lk(t_tree *tree);
void RRW_Forward_Lk_Pre(t_node *a, t_node *d, phydbl *lnP, t_tree *tree);
phydbl RRW_Lk(t_tree *tree);
phydbl RRW_Prior_Sigsq_Scale(t_tree *tree);
phydbl RRW_Forward_Lk_Path(t_node *a, t_node *d, t_tree *tree);
void RRW_Rescale_Times(int prod, t_tree *tree);
void RRW_Rescale_Times_Pre(t_node *a, t_node *d, phydbl cur_ta, int prod, t_tree *tree);
phydbl RRW_Independent_Contrasts(t_tree *tree);


#endif
