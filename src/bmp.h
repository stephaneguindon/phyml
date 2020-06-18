/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include <config.h>

#ifndef BMP_H
#define BMP_H

#include "utilities.h"

int BMP_Main(int argc, char **argv);
phydbl BMP_Independent_Contrasts(t_tree *tree);
void BMP_Independent_Contrasts_Post(t_node *a, t_node *d, phydbl *lnP, t_tree *tree);
phydbl BMP_Forward_Lk(t_tree *tree);
void BMP_Forward_Lk_Pre(t_node *a, t_node *d, phydbl *lnP, t_tree *tree);
phydbl BMP_Lk(t_tree *tree);

#endif
