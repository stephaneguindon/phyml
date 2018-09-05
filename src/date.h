/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include <config.h>

#ifndef DATE_H
#define DATE_H

#include "utilities.h"

int DATE_Main(int argc, char **argv);
void DATE_XML(char *xml_filename);
void DATE_Update_Secondary_Cal(t_tree *tree);
void DATE_Update_Secondary_Cal_Post(t_node *a, t_node *d, t_tree *tree);
void DATE_Update_Secondary_Cal_Pre(t_node *a, t_node *d, t_tree *tree);
phydbl *DATE_Splitted_Calibration(t_tree *tree);
void DATE_Assign_Primary_Calibration(t_tree *tree);
void DATE_Update_T_Prior_MinMax(t_tree *tree);
phydbl DATE_J(phydbl birth_r, phydbl death_r, phydbl t_min, phydbl t_pls);
int DATE_Is_Split_Accessible(t_node *d, int which, phydbl *splitted_cal, t_tree *tree);
phydbl *DATE_Splitted_Calibration(t_tree *tree);
phydbl DATE_J_Sum_Product(t_tree *tree);
int DATE_J_Sum_Product_Pre(t_node *d, int split_idx_d, int split_idx_a, phydbl prod, int fact, phydbl *total, phydbl *splitted_cal, int rk, t_tree *tree);
void DATE_Chain_Cal(t_tree *mixt_tree);
int DATE_Check_Calibration_Constraints(t_tree *tree);
int DATE_Check_Time_Constraints(t_tree *tree);
phydbl *DATE_MCMC(t_tree *tree);
void DATE_List_Of_Nodes_Younger_Than(t_node *a, t_node *d, phydbl lim, t_ll **list, t_tree *tree);
void DATE_List_Of_Nodes_And_Ancestors_Younger_Than(t_node *a, t_node *d, phydbl lim, t_ll **list, t_tree *tree);
t_ll *DATE_List_Of_Regraft_Nodes(t_node *prune, t_node *prune_daughter, phydbl *t_min, phydbl *t_max, int verbose, t_tree *tree);
phydbl DATE_Lk_Calib(t_tree *tree);




#endif
