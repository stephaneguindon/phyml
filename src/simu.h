/*

PHYML :  a program that  computes maximum likelihood  phylogenies from
DNA or AA homologous sequences 

Copyright (C) Stephane Guindon. Oct 2003 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/

#include <config.h>

#ifndef CURR_H
#define CURR_H

#include "utilities.h"

void Simu_Loop(t_tree *tree);
int Simu(t_tree *tree,int n_step_max);
void Select_Edges_To_Swap(t_tree *tree,t_edge **sorted_b,int *n_neg);
void Update_Bl(t_tree *tree,phydbl fact);
void Make_N_Swap(t_tree *tree,t_edge **b,int beg,int end);
int Make_Best_Swap(t_tree *tree);
int Mov_Backward_Topo_Bl(t_tree *tree,phydbl lk_old,t_edge **tested_b,int n_tested);
void Unswap_N_Branch(t_tree *tree,t_edge **b,int beg,int end);
void Swap_N_Branch(t_tree *tree,t_edge **b,int beg,int end);
void Check_NNI_Scores_Around(t_node *a, t_node *d, t_edge *b, phydbl *best_score, t_tree *tree);
int Mov_Backward_Topo_Pars(t_tree *tree, int pars_old, t_edge **tested_b, int n_tested);
void Simu_Pars(t_tree *tree, int n_step_max);

#endif
