/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include <config.h>

#ifndef TIPORDER_H
#define TIPORDER_H

#include "times.h"
#include "spr.h"
#include "utilities.h"
#include "lk.h"
#include "optimiz.h"
#include "bionj.h"
#include "models.h"
#include "free.h"
#include "help.h"
#include "simu.h"
#include "eigen.h"
#include "pars.h"
#include "alrt.h"
#include "m4.h"
#include "draw.h"
#include "rates.h"
#include "mcmc.h"
#include "io.h"
#include "stats.h"

void TIPO_Get_Tips_Y_Rank(t_tree *tree);
void TIPO_Get_Tips_Y_Rank_Pre(t_node *a, t_node *d, phydbl *curr_rank, t_tree *tree);
void TIPO_Get_All_Y_Rank(t_tree *tree);
void TIPO_Get_All_Y_Rank_Pre(t_node *a, t_node *d, t_tree *tree);
void TIPO_Swap_One_Node(t_node *d, t_tree *tree);
void TIPO_Minimize_Tip_Order_Score(int n_trees, t_tree **list_tree, t_tree *ref_tree);
void TIPO_Print_Tip_Ordered(t_tree *ref_tree);
void TIPO_Print_Tip_Ordered_Pre(t_node *a, t_node *d, t_tree *ref_tree);
phydbl TIPO_Untangle_Tree(t_tree *tree);
void TIPO_Untangle_Node(t_node *a, t_node *d, t_node **node_table, int *conflict, t_tree *tree);
int TIPO_Untangle_Tree_List(int n_trees, t_tree **list_tree, t_tree *ref_tree);
int TIPO_Check_Tip_Ranks(t_tree *tree);
void TIPO_Read_Taxa_Coordinates(FILE *fp_coord, t_tree *tree);
void TIPO_Get_Tips_Y_Rank_From_Zscores(t_tree *tree);
void TIPO_Init_Tip_Num(t_tree *tree);
void TIPO_Read_Taxa_Zscores(FILE *fp_coord, t_tree *tree);
void  TIPO_Sort_Translation_Table(t_tree *tree);
void TIPO_Randomize_Tip_Y_Ranks(t_tree *tree);
phydbl TIPO_Read_One_Taxon_Zscore(FILE *fp, char *seqname_qry, int col, t_tree *tree);
void TIPO_Normalize_Zscores(t_tree *tree);
void TIPO_Get_Min_Number_Of_Tip_Permut(t_tree *tree);
phydbl TIPO_Lk(t_tree *tree);
phydbl TIPO_Lk_Post(t_node *a, t_node *d, t_tree *tree);
phydbl TIPO_Lk_Core(t_node *a, t_node *d, t_tree *tree);


#endif
