/*

PHYML :  a program that  computes maximum likelihood  phylogenies from
DNA or AA homologous sequences 

Copyright (C) Stephane Guindon. Oct 2003 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/

#include <config.h>

#ifndef EVOLVE_H
#define EVOLVE_H

#include "utilities.h"

int EVOLVE_Main(int argc, char **argv);
void EVOLVE_Coalescent(t_tree *tree);
void EVOLVE_Seq(calign *data, t_mod *mod, FILE *fp, t_tree *tree);
int EVOLVE_Pick_State(int n,phydbl *prob);
void EVOLVE_Seq_Recur(t_node *a,t_node *d,t_edge *b,int a_state,int r_class,int site_num,calign *gen_data,t_mod *mod,t_tree *tree);
phydbl *EVOLVE_Site_Lk(int site_idx, int rate_class, calign *data, t_tree *tree);

#endif
