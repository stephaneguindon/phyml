/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include <config.h>

#ifndef VELOC_H
#define VELOC_H

#include "utilities.h"

phydbl VELOC_Wrap_Lk(t_edge *b, t_tree *tree, supert_tree *stree);
short int VELOC_Is_Integrated_Velocity(t_phyrex_mod *mod);
phydbl VELOC_Lk(t_tree *tree);

phydbl VELOC_Augmented_Lk_Locations(t_tree *tree);
phydbl VELOC_Augmented_Lk_Locations_Core(t_dsk *disk, t_tree *tree);
phydbl VELOC_Locations_Forward_Lk_Path(t_ldsk *a, t_ldsk *d, t_tree *tree);

phydbl VELOC_Augmented_Lk_Velocities(t_tree *tree);
phydbl VELOC_Augmented_Lk_Velocities_Core(t_dsk *disk, t_tree *tree);
phydbl VELOC_Velocities_Forward_Lk_Path(t_ldsk *a, t_ldsk *d, t_tree *tree);

phydbl VELOC_Integrated_Lk_Location(t_tree *tree);
void VELOC_Integrated_Lk_Location_Post(t_node *a, t_node *d, short int dim, t_tree *tree, short int print);
void VELOC_Integrated_Lk_Location_Pre(t_node *a, t_node *d, short int dim, t_tree *tree);

phydbl VELOC_Velocity_Variance_Along_Edge(t_node *d, short int dim, t_tree *tree);
phydbl VELOC_Location_Variance_Along_Edge(t_node *d, short int dim, t_tree *tree);
phydbl VELOC_Velocity_Mean_Along_Edge(t_node *d, short int dim, t_tree *tree);
phydbl VELOC_Location_Mean_Along_Edge(t_node *d, short int dim, t_tree *tree);

void VELOC_Sample_Node_Locations(t_tree *tree);

void VELOC_Integrated_Lk_Velocity_Post(t_node *a, t_node *d, short int dim, t_tree *tree, short int print);
void VELOC_Integrated_Lk_Velocity_Pre(t_node *a, t_node *d, short int dim, t_tree *tree);
void VELOC_Veloc_Gibbs_Mean_Var(t_node *n, phydbl *mean, phydbl *var, short int dim, t_tree *tree);
phydbl VELOC_Mean_Speed(t_tree *tree);
phydbl VELOC_Veloc_To_Speed(t_geo_veloc *v, t_tree *tree);
phydbl VELOC_Mean_Velocity(short int dim, t_tree *tree);
phydbl PHYREX_Degrees_To_Km(phydbl deg, t_tree *tree);


#endif
