/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include <config.h>

#ifndef IWN_H
#define IWN_H

#include "utilities.h"

short int IWN_Is_Iwn(t_phyrex_mod *mod);

int IWN_Get_Gt(phydbl t, phydbl omega);

phydbl IWN_Velocity_Variance_Along_Edge(t_node *d, short int dim, t_tree *tree);
phydbl IWN_Velocity_Mean_Along_Edge(t_node *d, short int dim, t_tree *tree);

phydbl IWN_Location_Variance_Along_Edge(t_node *d, short int dim, t_tree *tree);
phydbl IWN_Location_Mean_Along_Edge(t_node *d, short int dim, t_tree *tree);

phydbl IWN_Prior(t_tree *tree);
phydbl IWN_Prior_Omega(t_tree *tree);

void IWN_Integrated_Location_Down(phydbl dt1, phydbl dt2,
                                  phydbl av1, phydbl bv1, phydbl v1mu, phydbl v1var, phydbl dv1var,
                                  phydbl av2, phydbl bv2, phydbl v2mu, phydbl v2var, phydbl dv2var,
                                  phydbl v1logrem, phydbl v2logrem,
                                  phydbl veloc_d, phydbl veloc_v1, phydbl veloc_v2,
                                  phydbl omega,
                                  phydbl *mean, phydbl *var, phydbl *logrem);

void IWN_Integrated_Location_Up(phydbl dt1, phydbl dt2,
                                phydbl av1, phydbl bv1, phydbl v1mu, phydbl v1var, phydbl av1var,
                                phydbl av2, phydbl bv2, phydbl v2mu, phydbl v2var, phydbl av2var,
                                phydbl v1logrem, phydbl v2logrem,
                                phydbl veloc_a, phydbl veloc_v1, phydbl veloc_v2,
                                phydbl omega,
                                phydbl *mean, phydbl *var, phydbl *logrem,
                                short int a_is_root);



#endif
