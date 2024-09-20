/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include <config.h>

#ifndef RW_H
#define RW_H

#include "utilities.h"

phydbl RW_Independent_Contrasts(t_tree *tree);
void RW_Independent_Contrasts_Post(t_node *a, t_node *d, phydbl sd, phydbl *lnP, t_tree *tree);
phydbl RW_Lk(t_tree *tree);
phydbl RW_Lk_Range(t_dsk *young, t_dsk *old, t_tree *tree);
phydbl RW_Prior(t_tree *tree);
phydbl RW_Prior_Sigsq(t_tree *tree);
void RW_Integrated_Lk_Down(phydbl son_a, phydbl son_b, phydbl son_mu_down, phydbl son_var_down, phydbl son_var,
                           phydbl bro_a, phydbl bro_b, phydbl bro_mu_down, phydbl bro_var_down, phydbl bro_var,
                           phydbl son_logrem, phydbl bro_logrem,
                           phydbl *mean, phydbl *var, phydbl *logrem);
void RW_Integrated_Lk_Up(phydbl dad_mu_up, phydbl dad_var_up, phydbl dad_logrem_up,
                         phydbl son_a, phydbl son_b, phydbl son_var,
                         phydbl bro_a, phydbl bro_b, phydbl bro_mu_down, phydbl bro_var_down, phydbl bro_var, phydbl bro_logrem_down,
                         phydbl *mean, phydbl *var, phydbl *logrem);
phydbl RW_Prior_Observational_Model(t_tree *tree);

#endif
