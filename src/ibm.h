/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include <config.h>

#ifndef IBM_H
#define IBM_H

#include "utilities.h"

phydbl IBM_Lk_Locations(t_tree *tree);
void IBM_Lk_Locations_Post(t_node *a, t_node *d, phydbl sigsq, int dim, t_tree *tree, short int print);
phydbl IBM_Lk_Velocities(t_tree *tree);
phydbl IBM_Lk(t_tree *tree);
void IBM_Sample_Locations(t_tree *tree);
phydbl IBM_Velocities_Forward_Lk_Path(t_ldsk *a, t_ldsk *d, t_tree *tree);
phydbl IBM_Lk_Velocities_Core(t_dsk *disk, t_tree *tree);
phydbl IBM_Lk_Velocities(t_tree *tree);
phydbl IBM_Locations_Conditional(short int sample, int *node_order, t_tree *tree);
phydbl IBM_Locations_Forward_Lk_Path(t_ldsk *a, t_ldsk *d, t_tree *tree);
phydbl IBM_Lk_Locations_Core(t_dsk *disk, t_tree *tree);
phydbl IBM_Sample_Velocities_And_Locations(int *node_order, t_tree *tree);
phydbl IBM_Velocities_Conditional(short int sample, int *node_order, t_tree *tree);
phydbl IBM_Velocity_One_Node(t_node *n, short int sample, t_tree *tree);
phydbl IBM_Location_One_Node(t_node *n, short int sample, t_tree *tree);
short int IBM_Is_Ibm(t_phyrex_mod *mod);

#endif
