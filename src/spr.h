/*
** spr.h: Header file for the SPR routines.
**
** Wim Hordijk   Last modified: 28 August 2006
*/

#include <config.h>

#ifndef _SPR_H_
#define _SPR_H_

#include "utilities.h"

#define ALL   1
#define BEST  2
#define ONE   3

/*
** _move_: Structure for holding the relevant information for candidate SPR moves.
*/

typedef struct
{
  t_node   *v_prune, *u_prune, *v_n, *v_nx1, *u_n, **path;
  t_edge   *e_prune, *e_regraft;
  phydbl  l_connect, l_est[3], delta_lk, d_L, d_up_v, d_un_v;
  int     dist, rgrft_rank, optim_rank, globl_rank;
} _move_;



void Init_SPR          (t_tree *tree);
void Clean_SPR         (t_tree *tree);
void Optim_SPR         (t_tree *tree, int max_size, int method);
int  Perform_SPR_Moves (t_tree *tree, int max_size);
int  Perform_Best_SPR  (t_tree *tree, int max_size);
int  Perform_One_SPR   (t_tree *tree, int max_size);

void Prune            (t_edge *e, t_node *v, t_edge **e_connect, t_edge **e_avail,
		       t_tree *tree);
void Regraft          (t_edge *e, t_node *v, t_edge *avail, t_tree *tree);
void PostOrder_v      (t_tree *tree, t_node *v, t_edge *e);
void PostOrder_w      (t_tree *tree, t_node *v, t_edge *v_e, t_node *w, t_edge *e);





void Speed_Spr(t_tree *tree, phydbl prop_spr, int max_cycles, phydbl delta_lnL);
void Global_Spr_Search(t_tree *tree);
void Make_Spr_List(t_tree *tree);
void Init_One_Spr(t_spr *a_spr);
t_spr *Make_One_Spr(t_tree *tree);
int Spr(phydbl init_lnL, phydbl prop_spr, t_tree *tree);
int Spr_Recur(t_node *a, t_node *d, t_tree *tree);
int Test_All_Spr_Targets(t_edge *pulled, t_node *link, t_tree *tree);
void Randomize_Spr_List(t_tree *tree);
void Test_One_Spr_Target_Recur(t_node *a, t_node *d, t_edge *pulled, t_node *link, t_edge *residual, t_edge *init_target, int *best_found, t_spr *prev_move, t_tree *tree);
t_spr *Test_One_Spr_Target(t_edge *b_target, t_edge *b_arrow, t_node *n_link, t_edge *b_residual, t_edge *init_target, t_node *polarity, t_tree *tree);
void Apply_Spr_Moves_One_By_One(t_tree *tree);
int Try_One_Spr_Move_Triple(t_spr *move, t_tree *tree);
int Try_One_Spr_Move_Full(t_spr *move, short int apply_move, t_tree *tree);
void Make_Best_Spr(t_tree *tree);
void Random_Spr(int n_moves, t_tree *tree);
unsigned int Include_One_Spr_To_List_Of_Spr(t_spr **list, int list_size, t_spr *move, t_tree *tree);
void Reset_Spr_List(t_spr **list, int size);
int Evaluate_List_Of_Regraft_Pos_Triple(t_spr **spr_list, int list_size, t_tree *tree);
void Best_Spr(t_tree *tree);
int Check_Spr_Move_Validity(t_spr *this_spr_move, t_tree *tree);
void Spr_Subtree(t_edge *b, t_node *link, t_tree *tree);
void Spr_Pars(int threshold, int n_round_max, t_tree *tree);
void Spr_Shuffle(t_tree *tree);
void Sort_Spr_List_Depth(t_tree *tree);
void Sort_Spr_List_LnL(t_spr **list, int list_size, t_tree *tree);
void Spr_Random_Explore(t_tree *tree, phydbl anneal_temp, phydbl prop_spr, int do_rnd, int max_cycles);
void Sort_Spr_List_Pars(t_tree *tree);
void Spr_List_Of_Trees(t_tree *tree);
void Spr_Pre_Order(t_node *a, t_node *d, t_edge *b, t_tree *tree);



#endif  /* _SPR_H_ */


/*
** EOF: spr.h
*/
