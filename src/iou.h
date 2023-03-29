/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include <config.h>

#ifndef IOU_H
#define IOU_H

#include "utilities.h"

short int IOU_Is_Iou(t_phyrex_mod *mod);

phydbl IOU_Velocity_Variance_Along_Edge(t_node *d, short int dim, t_tree *tree);
phydbl IOU_Velocity_Mean_Along_Edge(t_node *d, short int dim, t_tree *tree);

phydbl IOU_Location_Variance_Along_Edge(t_node *d, short int dim, t_tree *tree);
phydbl IOU_Location_Mean_Along_Edge(t_node *d, short int dim, t_tree *tree);

phydbl IOU_Prior(t_tree *tree);
phydbl IOU_Prior_Theta(t_tree *tree);
phydbl IOU_Prior_Mu(t_tree *tree);

void IOU_Integrated_Location_Down(phydbl dt1, phydbl dt2,
                                  phydbl av1, phydbl bv1, phydbl v1mu, phydbl v1var, phydbl dv1var,
                                  phydbl av2, phydbl bv2, phydbl v2mu, phydbl v2var, phydbl dv2var,
                                  phydbl v1logrem, phydbl v2logrem,
                                  phydbl *mean, phydbl *var, phydbl *logrem);

void IOU_Integrated_Location_Up(phydbl dt1, phydbl dt2,
                                phydbl av1, phydbl bv1, phydbl v1mu, phydbl v1var, phydbl av1var,
                                phydbl av2, phydbl bv2, phydbl v2mu, phydbl v2var, phydbl av2var,
                                phydbl v1logrem, phydbl v2logrem,
                                phydbl *mean, phydbl *var, phydbl *logrem,
                                short int a_is_root);



#endif
