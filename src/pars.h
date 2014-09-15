/*

PHYML :  a program that  computes maximum likelihood  phylogenies from
DNA or AA homologous sequences 

Copyright (C) Stephane Guindon. Oct 2003 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/

#include <config.h>

#ifndef PARS_H
#define PARS_H

#include "utilities.h"
#include "lk.h"
#include "optimiz.h"
#include "models.h"
#include "free.h"

int  Pars(t_edge *b, t_tree *tree);
void Post_Order_Pars(t_node *a, t_node *d, t_tree *tree);
void Pre_Order_Pars(t_node *a, t_node *d, t_tree *tree);
void Get_Partial_Pars(t_tree *tree, t_edge *b_fcus, t_node *a, t_node *d);
void Site_Pars(t_tree *tree);
void Init_Ui_Tips(t_tree *tree);
void Update_P_Pars(t_tree *tree, t_edge *b_fcus, t_node *n);
int Pars_At_Given_Edge(t_edge *b, t_tree *tree);
void Get_All_Partial_Pars(t_tree *tree, t_edge *b_fcus, t_node *a, t_node *d);
int Update_Pars_At_Given_Edge(t_edge *b_fcus, t_tree *tree);
void Init_P_Pars_Tips(t_tree *tree);
void Get_Step_Mat(t_tree *tree);
int Pars_Core(t_edge *b, t_tree *tree);
int One_Pars_Step(t_edge *b,t_tree *tree);

#endif
