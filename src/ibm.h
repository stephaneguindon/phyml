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

short int IBM_Is_Ibm(t_phyrex_mod *mod);
phydbl IBM_Sample_Velocities_And_Locations(int *node_order, int n_nodes, t_tree *tree);
phydbl IBM_Velocities_Conditional(short int sample, int *node_order, short int dim, t_tree *tree);
phydbl IBM_Velocity_One_Node(t_node *n, short int sample, short int dim, t_tree *tree);
phydbl IBM_Locations_Conditional(short int sample, int *node_order, short int dim, t_tree *tree);
phydbl IBM_Location_One_Node(t_node *n, short int sample, short int dim, t_tree *tree);
void IBM_Augmented_Lk_Locations_Post(t_node *a, t_node *d, phydbl sigsq, int dim, t_tree *tree, short int print);
void IBM_Generate_Velocities_Then_Locations(t_tree *tree);
void IBM_Generate_Velocities(t_tree *tree);
void IBM_Generate_Velocities_Pre(t_node *a, t_node *d, t_tree *tree);
void IBM_Generate_Locations_Given_Velocities(t_tree *tree);
void IBM_Generate_Locations_Given_Velocities_Pre(t_node *a, t_node *d, t_tree *tree);
phydbl IBM_Velocity_Variance_Along_Edge(t_node *d, short int dim, t_tree *tree);
phydbl IBM_Velocity_Mean_Along_Edge(t_node *d, short int dim, t_tree *tree);
phydbl IBM_Location_Variance_Along_Edge(t_node *d, short int dim, t_tree *tree);
phydbl IBM_Location_Mean_Along_Edge(t_node *d, short int dim, t_tree *tree);
phydbl IBM_Prior(t_tree *tree);


void IBM_Integrated_Location_Down(phydbl dt1, phydbl dt2,
                                  phydbl av1, phydbl bv1, phydbl v1mu, phydbl v1var, phydbl dv1var,
                                  phydbl av2, phydbl bv2, phydbl v2mu, phydbl v2var, phydbl dv2var,
                                  phydbl v1logrem, phydbl v2logrem,
                                  phydbl *mean, phydbl *var, phydbl *logrem);

void IBM_Integrated_Location_Up(phydbl dt1, phydbl dt2,
                                phydbl av1, phydbl bv1, phydbl v1mu, phydbl v1var, phydbl av1var,
                                phydbl av2, phydbl bv2, phydbl v2mu, phydbl v2var, phydbl av2var,
                                phydbl v1logrem, phydbl v2logrem,
                                phydbl *mean, phydbl *var, phydbl *logrem,
                                short int is_root);
phydbl IBM_Mean_Speed(t_tree *tree);

#endif
