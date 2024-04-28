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

void IOU_Integrated_Location_Down(phydbl son_a, phydbl son_b, phydbl son_mu_down, phydbl son_var_down, phydbl son_var,
                                  phydbl bro_a, phydbl bro_b, phydbl bro_mu_down, phydbl bro_var_down, phydbl bro_var,
                                  phydbl son_logrem, phydbl bro_logrem,
                                  phydbl *mean, phydbl *var, phydbl *logrem);

void IOU_Integrated_Location_Up(phydbl dad_mu_up, phydbl dad_var_up, phydbl dad_logrem_up,
                                phydbl son_a, phydbl son_b, phydbl son_var,
                                phydbl bro_a, phydbl bro_b, phydbl bro_mu_down, phydbl bro_var_down, phydbl bro_var, phydbl bro_logrem_down,
                                

                                phydbl *mean, phydbl *var, phydbl *logrem);


#endif
