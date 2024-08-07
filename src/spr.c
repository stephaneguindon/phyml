/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

/*
** spr.c: Routines for performing SPR moves on the tree.
**
** Wim Hordijk   Last modified: 28 August 2006
** Stephane Guindon 2007
*/

#include "spr.h"

#ifdef BEAGLE
#include "beagle_utils.h"
#endif

/*********************************************************/
/* Sort list of SPR move by putting the shallowest moves first */
void Sort_Spr_List_Depth(t_tree *tree)
{
  int i,j;
  t_spr *buff;
  short int swap;

  do
    {
      swap = NO;
  
      for(i=0;i<tree->n_moves-1;i++)
        {
          for(j=i+1;j<tree->n_moves;j++)
            {
              if(tree->spr_list_one_edge[j]->depth_path < tree->spr_list_one_edge[i]->depth_path)
                {
                  buff                       = tree->spr_list_one_edge[i];
                  tree->spr_list_one_edge[i] = tree->spr_list_one_edge[j];
                  tree->spr_list_one_edge[j] = buff;
                  swap = YES;
                }
            }
        }
    }
  while(swap == YES);
}

/*********************************************************/
/*********************************************************/
/* Sort list of SPR move by putting the more likely moves first */
void Sort_Spr_List_LnL(t_spr **list, int list_size, t_tree *tree)
{
  int i,j;
  t_spr *buff;
  short int swap;

  do
    {
      swap = NO;
      for(i=0;i<list_size-1;i++)
        {
          for(j=i+1;j<list_size;j++)
            {
              if(list[j]->lnL > list[i]->lnL)
                {
                  buff    = list[i];
                  list[i] = list[j];
                  list[j] = buff;
                  swap = YES;
                }
            }
        }
    }
  while(swap == YES);
}


/*********************************************************/
/*********************************************************/
/* Sort list of SPR move by putting the more parsimonious moves first */
void Sort_Spr_List_Pars(t_tree *tree)
{
  int i,j;
  t_spr *buff;
  short int swap;

  do
    {
      swap = NO;
      
      for(i=0;i<tree->size_spr_list_one_edge-1;i++)
        {
          for(j=i+1;j<tree->size_spr_list_one_edge;j++)
            {
              if(tree->spr_list_one_edge[j]->pars < tree->spr_list_one_edge[i]->pars)
                {
                  buff              = tree->spr_list_one_edge[i];
                  tree->spr_list_one_edge[i] = tree->spr_list_one_edge[j];
                  tree->spr_list_one_edge[j] = buff;
                  swap = YES;
                }
            }
        }
    }
  while(swap == YES);
}


/*********************************************************/
/*********************************************************/
/*********************************************************/

void Randomize_Spr_List(t_tree *tree)
{
  int i,j;
  t_spr *buff;

  for(i=0;i<tree->n_moves;i++)
    {
      j = Rand_Int(0,tree->n_moves-1);
      buff                       = tree->spr_list_one_edge[i];
      tree->spr_list_one_edge[i] = tree->spr_list_one_edge[j];
      tree->spr_list_one_edge[j] = buff;
    }
}

/*********************************************************/

int Spr(phydbl init_lnL, phydbl prop_spr, t_tree *tree)
{
  int i,br;
  int *br_idx;
  t_edge *b;

  tree->mod->s_opt->n_improvements = 0;
  tree->mod->s_opt->max_spr_depth  = 0;
  tree->mod->s_opt->max_rank_triple_move = 0;
  tree->mod->s_opt->max_delta_lnL_spr_current = 0.0;
  
  Reset_Spr_List(tree->spr_list_all_edge,tree->size_spr_list_all_edge);

  br_idx = Permutate(2*tree->n_otu-3);

  Set_Both_Sides(YES,tree);
  Lk(NULL,tree);
  tree->best_lnL = tree->c_lnL;
  
  for(i=0;i<MAX(1,(int)((2*tree->n_otu-3)*prop_spr));++i)
    {
      br = br_idx[i];
      
      if(!(br%10)) if(tree->io->print_json_trace == YES) JSON_Tree_Io(tree,tree->io->fp_out_json_trace);

      b = tree->a_edges[br];
      
      if(b->l->v > tree->mod->s_opt->l_min_spr)
        {
          Spr_Subtree(b,b->left,tree);
          Spr_Subtree(b,b->rght,tree);
        }
    }
  
  Free(br_idx);
    
  if(tree->mod->s_opt->n_improvements == 0 && tree->mod->s_opt->spr_lnL == YES)
    {
      Optimize_Br_Len_Serie(2,tree);
      tree->best_lnL = tree->c_lnL;

      Sort_Spr_List_LnL(tree->spr_list_all_edge,tree->size_spr_list_all_edge,tree);
      for(i=0;i<MIN(20,2*tree->n_otu-3);++i)
        {
          Try_One_Spr_Move_Full(tree->spr_list_all_edge[i],NO,tree);
          if(tree->spr_list_all_edge[i]->lnL > tree->best_lnL + tree->mod->s_opt->min_diff_lk_move) break;
        }
      Sort_Spr_List_LnL(tree->spr_list_all_edge,tree->size_spr_list_all_edge,tree);
      if(tree->spr_list_all_edge[0]->lnL > tree->best_lnL + tree->mod->s_opt->min_diff_lk_move)
        Try_One_Spr_Move_Full(tree->spr_list_all_edge[0],YES,tree);
    }
  
  return tree->mod->s_opt->n_improvements;
}

/*********************************************************/

void Spr_Pre_Order(t_node *a, t_node *d, t_edge *b, t_tree *tree)
{
  if(d->tax) return;
  else
    {
      unsigned int i;
      

      /* printf("\n. a: %d d: %d score: %d d1: %d d2: %d ",a->num,d->num,tree->c_pars); */

      /* for(i=0;i<3;++i) */
      /*   { */
      /*     if(d->v[i] != a) */
      /*       { */
      /*         Spr_Subtree(d->b[i],d->v[i],tree); */
      /*       } */
      /*   } */

      Spr_Subtree(b,a,tree);
      
      for(i=0;i<3;++i)
        {
          if(d->v[i] != a)
            {
              Spr_Pre_Order(d,d->v[i],d->b[i],tree);
            }
        }

    }
}

/*********************************************************/

void Spr_Subtree(t_edge *b, t_node *link, t_tree *tree)
{
  int i;
  int n_moves_pars, n_moves, min_pars, best_move_idx;
  t_spr *best_pars_move;
  t_edge *init_target, *dummy, *residual;
    
  if(link->v[0] == NULL || link->v[1] == NULL || link->v[2] == NULL) return;

  Reset_Spr_List(tree->spr_list_one_edge,tree->size_spr_list_one_edge);
 
  best_move_idx = -1;
  tree->n_moves = 0;

  MIXT_Set_Pars_Thresh(tree);

  if((link != b->left) && (link != b->rght))
    {
      PhyML_Fprintf(stderr,"\n== Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
      Exit("\n");
    }
  else
    {
      tree->mod->s_opt->worst_lnL_spr = BIG;

      if(!link->tax) Test_All_Spr_Targets(b,link,tree);

      if(tree->n_moves)
        {
          n_moves_pars = MIN(tree->mod->s_opt->min_n_triple_moves,tree->n_moves);
          n_moves      = MIN(tree->mod->s_opt->min_n_triple_moves,tree->n_moves);

          if(tree->mod->s_opt->spr_lnL == NO) n_moves = n_moves_pars;
          n_moves = MAX(1,n_moves);
          
          if(tree->mod->s_opt->spr_pars == YES)
            {
              min_pars = 1E+8;
              best_pars_move = NULL;
              
              for(i=0;i<n_moves;++i)
                {              
                  if(tree->spr_list_one_edge[i]->pars < min_pars)
                    {
                      best_pars_move = tree->spr_list_one_edge[i];
                      min_pars = tree->spr_list_one_edge[i]->pars;
                    }
                }
              
              assert(best_pars_move);
              
              if(best_pars_move->pars < tree->best_pars)
                {
                  Prune_Subtree(best_pars_move->n_link,best_pars_move->n_opp_to_link,&init_target,&residual,tree);
                  Graft_Subtree(best_pars_move->b_target,best_pars_move->n_link,NULL,residual,NULL,tree);                  
                  if(!Check_Topo_Constraints(tree,tree->io->cstr_tree))
                    {
                      Prune_Subtree(best_pars_move->n_link,best_pars_move->n_opp_to_link,&dummy,&residual,tree);
                      Graft_Subtree(init_target,best_pars_move->n_link,NULL,residual,NULL,tree);                  
                      Set_Both_Sides(YES,tree);
                      Pars(NULL,tree);
                    }
                  else
                    {
                      if(best_pars_move->depth_path > tree->mod->s_opt->max_spr_depth) tree->mod->s_opt->max_spr_depth = best_pars_move->depth_path;
                      Set_Both_Sides(YES,tree);
                      Pars(NULL,tree);
                      tree->best_pars = tree->c_pars;
                      if(tree->best_pars != best_pars_move->pars)
                        {
                          PhyML_Fprintf(stderr,"\n== best_pars = %d move_pars = %d",tree->best_pars,best_pars_move->pars);
                          PhyML_Fprintf(stderr,"\n== Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
                          Exit("\n");
                        }
                      tree->mod->s_opt->n_improvements++;
                    }
                  Reset_Spr_List(tree->spr_list_one_edge,tree->size_spr_list_one_edge);
                }
              else
                {
                  Set_Both_Sides(YES,tree);
                  Pars(NULL,tree);
                }
            }
          else
            {
              int apply_move = NO;
              phydbl accept_prob,u;

              if(tree->mod->s_opt->spr_lnL == YES) 
                {
                  Sort_Spr_List_LnL(tree->spr_list_one_edge,tree->size_spr_list_one_edge,tree);

                  if(tree->spr_list_one_edge[0]->lnL > tree->best_lnL)
                    {
                      best_move_idx = 0;
                    }
                  else if(tree->mod->s_opt->eval_list_regraft == YES)
                    {
                      best_move_idx = Evaluate_List_Of_Regraft_Pos_Triple(tree->spr_list_one_edge,n_moves,tree);
                    }
                  else
                    {
                      best_move_idx = -1;
                    }
                }
              else
                {
                  best_move_idx = Evaluate_List_Of_Regraft_Pos_Triple(tree->spr_list_one_edge,n_moves,tree);
                }

              
              if(best_move_idx > -1)
                {
                  if(Are_Equal(tree->annealing_temp,0.0,1.E-3) == NO)
                    {
                      accept_prob = exp((tree->spr_list_one_edge[best_move_idx]->lnL - tree->best_lnL)/tree->annealing_temp);
                      u = Uni();
                      if(!(u > accept_prob)) apply_move = YES;
                    }
                  else
                    {
                      if(tree->spr_list_one_edge[best_move_idx]->lnL > tree->best_lnL + tree->mod->s_opt->min_diff_lk_move)
                        apply_move = YES;
                    }
                }

              if((best_move_idx > -1) && (apply_move == YES))
                {
                  Try_One_Spr_Move_Triple(tree->spr_list_one_edge[best_move_idx],tree);
                }
              else
                {
                  Pars(NULL,tree);
                }
            }
        }
      Reset_Spr_List(tree->spr_list_one_edge,tree->size_spr_list_one_edge);
    }
}

/*********************************************************/

int Test_All_Spr_Targets(t_edge *b_pulled, t_node *n_link, t_tree *tree)
{
  t_node *n_opp_to_link,*n_v1,*n_v2;
  t_edge *b_target,*b_residual;
  int i,dir1,dir2;
  scalar_dbl *init_l_v1, *init_l_v2, *init_l_pulled;
  scalar_dbl *init_v_v1, *init_v_v2, *init_v_pulled;
  int best_found;
  phydbl init_lnL;

  if(tree->mixt_tree != NULL)
    {
      PhyML_Fprintf(stderr,"\n== Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
      Exit("\n");
    }

  init_lnL = tree->c_lnL;
  b_target = b_residual = NULL;
  n_opp_to_link  = (n_link == b_pulled->rght)?(b_pulled->left):(b_pulled->rght);

  init_l_pulled = Duplicate_Scalar_Dbl(b_pulled->l);
  init_v_pulled = Duplicate_Scalar_Dbl(b_pulled->l_var);

  dir1 = dir2 = -1;
  for(i=0;i<3;i++)
    if(n_link->v[i] != n_opp_to_link)
      {
        if(dir1<0) dir1 = i;
        else       dir2 = i;
      }

  assert(dir1 > -1);
  assert(dir2 > -1);

  if(n_link->v[dir1]->num < n_link->v[dir2]->num)
    {
      n_v1      = n_link->v[dir1];
      n_v2      = n_link->v[dir2];
      init_l_v1 = Duplicate_Scalar_Dbl(n_link->b[dir1]->l);
      init_l_v2 = Duplicate_Scalar_Dbl(n_link->b[dir2]->l);
      init_v_v1 = Duplicate_Scalar_Dbl(n_link->b[dir1]->l_var);
      init_v_v2 = Duplicate_Scalar_Dbl(n_link->b[dir2]->l_var);
    }
  else
    {
      n_v1      = n_link->v[dir2];
      n_v2      = n_link->v[dir1];
      init_l_v1 = Duplicate_Scalar_Dbl(n_link->b[dir2]->l);
      init_l_v2 = Duplicate_Scalar_Dbl(n_link->b[dir1]->l);
      init_v_v1 = Duplicate_Scalar_Dbl(n_link->b[dir2]->l_var);
      init_v_v2 = Duplicate_Scalar_Dbl(n_link->b[dir1]->l_var);
    }

  if(!(n_v1->tax && n_v2->tax)) /*! Pruning is meaningless otherwise */
    {      
      Prune_Subtree(n_link,n_opp_to_link,&b_target,&b_residual,tree);

      if(tree->mod->s_opt->spr_lnL == YES) Update_PMat_At_Given_Edge(b_target,tree);
      
      for(i=0;i<tree->size_spr_list_one_edge;++i) tree->spr_list_one_edge[i]->path_prev = NULL;

      tree->edge_list = NULL;
      tree->node_list = NULL;
      best_found = NO;
      tree->depth_curr_path = 0;
      tree->curr_path[0] = b_target->left;
      Test_One_Spr_Target_Recur(b_target->rght,
                                b_target->left,
                                b_pulled,n_link,b_residual,b_target,&best_found,NULL,tree);

      if(best_found == NO || tree->perform_spr_right_away == NO)
        {
          tree->depth_curr_path = 0;
          tree->curr_path[0] = b_target->rght;
          Test_One_Spr_Target_Recur(b_target->left,
                                    b_target->rght,
                                    b_pulled,n_link,b_residual,b_target,&best_found,NULL,tree);
        }

      Graft_Subtree(b_target,n_link,NULL,b_residual,NULL,tree);

      if((n_link->v[dir1] != n_v1) || (n_link->v[dir2] != n_v2)) PhyML_Printf("\n== Warning: -- SWITCH NEEDED -- ! \n");

      Copy_Scalar_Dbl(init_l_v1,n_link->b[dir1]->l);
      Copy_Scalar_Dbl(init_v_v1,n_link->b[dir1]->l_var);

      Copy_Scalar_Dbl(init_l_v2,n_link->b[dir2]->l);
      Copy_Scalar_Dbl(init_v_v2,n_link->b[dir2]->l_var);

      Copy_Scalar_Dbl(init_l_pulled,b_pulled->l);
      Copy_Scalar_Dbl(init_v_pulled,b_pulled->l_var);

      if(tree->mod->s_opt->spr_pars == NO)      
        {
          Update_PMat_At_Given_Edge(n_link->b[dir1],tree);
          Update_PMat_At_Given_Edge(n_link->b[dir2],tree);
          Update_PMat_At_Given_Edge(b_pulled,tree);
        }

      if(tree->mod->s_opt->spr_pars == NO)
      	{
          Update_Partial_Lk(tree,b_pulled,  n_link);
          Update_Partial_Lk(tree,b_target,  n_link);
          Update_Partial_Lk(tree,b_residual,n_link);
        }
      else
      	{
          Update_Partial_Pars(tree,b_pulled,  n_link);
          Update_Partial_Pars(tree,b_target,  n_link);
          Update_Partial_Pars(tree,b_residual,n_link);
        }

      t_ll *e_ll = tree->edge_list->head;
      t_ll *n_ll = tree->node_list->head;
      t_edge *e;
      t_node *n;
      do
        {
          assert(e_ll);
          assert(n_ll);

          e = (t_edge *)e_ll->v;
          n = (t_node *)n_ll->v;

          /* printf("\n. update on edge %d node %d",e->num,n->num); fflush(NULL); */
          if(tree->mod->s_opt->spr_lnL)
            Update_Partial_Lk(tree,e,n);
          else
            Update_Partial_Pars(tree,e,n);
            
          e_ll = e_ll->next;
          n_ll = n_ll->next;
        }
      while(e_ll != NULL);


      
      Free_Linked_List(tree->edge_list);
      Free_Linked_List(tree->node_list);
    }

  tree->c_lnL = init_lnL;

  Free_Scalar_Dbl(init_l_v1);
  Free_Scalar_Dbl(init_l_v2);
  Free_Scalar_Dbl(init_l_pulled);

  Free_Scalar_Dbl(init_v_v1);
  Free_Scalar_Dbl(init_v_v2);
  Free_Scalar_Dbl(init_v_pulled);

  return 0;
}

/*********************************************************/

void Test_One_Spr_Target_Recur(t_node *a, t_node *d, t_edge *pulled, t_node *link, t_edge *residual, t_edge *init_target, int *best_found, t_spr *prev_move, t_tree *tree)
{
  unsigned int i;
  t_spr *move,*next_move;
  
  move = next_move = NULL;
  
  if(*best_found == YES && tree->perform_spr_right_away == YES) return;
  
  if(d->tax) return;
  else
    {      
      for(i=0;i<3;++i)
        {
          if(d->v[i] != a)
            {
              
              if(tree->mod->s_opt->spr_pars == NO)
                Update_Partial_Lk(tree,d->b[i],d);
              else
                Update_Partial_Pars(tree,d->b[i],d);
              
              /* PhyML_Printf("\n push edge %d node %d",d->b[i]->num,d->num); */
              
              Push_Bottom_Linked_List(d->b[i],&tree->edge_list,NO);
              Push_Bottom_Linked_List(d,&tree->node_list,NO);

              tree->depth_curr_path++;

              tree->curr_path[tree->depth_curr_path] = d->v[i];


              if((tree->depth_curr_path <= tree->mod->s_opt->max_depth_path) &&
                 (tree->depth_curr_path >= tree->mod->s_opt->min_depth_path))
                {
                  
                  move = Test_One_Spr_Target(d->b[i],pulled,link,residual,init_target,d,tree);

                  move->path_prev = prev_move;


                  if((tree->mod->s_opt->spr_pars == NO  && move->lnL > tree->best_lnL + tree->mod->s_opt->min_diff_lk_move) ||
                     (tree->mod->s_opt->spr_pars == YES && move->pars < tree->best_pars))
                    {
                      *best_found = YES;
                    }
                }

              bool go_to_next =
                (tree->depth_curr_path < tree->mod->s_opt->max_depth_path &&
                ((tree->mod->s_opt->spr_pars == NO  && move && move->lnL > tree->best_lnL - tree->mod->s_opt->max_delta_lnL_spr) ||
                 tree->mod->s_opt->spr_pars == YES));

              /* bool go_to_next = tree->depth_curr_path < tree->mod->s_opt->max_depth_path; */
                           
              if(go_to_next == YES) Test_One_Spr_Target_Recur(d,d->v[i],pulled,link,residual,init_target,best_found,move,tree);

              tree->depth_curr_path--;
            }
        }
    }
}

/*********************************************************/

t_spr *Test_One_Spr_Target(t_edge *b_target, t_edge *b_arrow, t_node *n_link, t_edge *b_residual, t_edge *init_target, t_node *polarity, t_tree *tree)
{
  scalar_dbl *init_target_l, *init_arrow_l, *init_residual_l;
  scalar_dbl *init_target_v, *init_arrow_v, *init_residual_v;
  int i,dir_v0,dir_v1,dir_v2;
  scalar_dbl *l0,*l1,*l2;
  scalar_dbl *v0,*v1,*v2;
  t_node *n1,*n2;
  phydbl init_lnL;
  int init_pars;
  t_spr *move;
  unsigned int rk;
  
  if(tree->mixt_tree != NULL)
    {
      PhyML_Fprintf(stderr,"\n== Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
      Exit("\n");
    }

  tree->n_moves++;
  
  init_lnL  = tree->c_lnL;
  init_pars = tree->c_pars;
  
  move = tree->spr_list_one_edge[tree->size_spr_list_one_edge];

  if(move->init_target_l == NULL)
    {
      move->init_target_l = Duplicate_Scalar_Dbl(init_target->l);
      move->init_target_v = Duplicate_Scalar_Dbl(init_target->l_var);
    }
  else
    {
      Copy_Scalar_Dbl(init_target->l,    move->init_target_l);
      Copy_Scalar_Dbl(init_target->l_var,move->init_target_v);
    }


  /* PhyML_Printf("\n. TEST_ONE_SPR_TARGET"); */
  // Save edge lengths so that they can be recovered in the end
  init_target_l   = Duplicate_Scalar_Dbl(b_target->l);
  init_target_v   = Duplicate_Scalar_Dbl(b_target->l_var);

  init_arrow_l    = Duplicate_Scalar_Dbl(b_arrow->l);
  init_arrow_v    = Duplicate_Scalar_Dbl(b_arrow->l_var);

  init_residual_l = Duplicate_Scalar_Dbl(b_residual->l);
  init_residual_v = Duplicate_Scalar_Dbl(b_residual->l_var);

  Graft_Subtree(b_target,n_link,NULL,b_residual,NULL,tree);

  if(tree->mod->s_opt->spr_lnL == YES)
    {
      Update_PMat_At_Given_Edge(b_target,tree);
      Update_PMat_At_Given_Edge(b_residual,tree);
      Update_Partial_Lk(tree,b_arrow,n_link);
      Lk(b_arrow,tree);
    }
  else
    {
      Update_Partial_Pars(tree,b_arrow,n_link);
      Pars(b_arrow,tree);
    }

  n1 = (b_residual->left == n_link)?(b_residual->rght):(b_residual->left);
  n2 = (b_target->left   == n_link)?(b_target->rght):(b_target->left);
  dir_v1 = dir_v2 = dir_v0 = -1;
  for(i=0;i<3;++i)
    {
      if(n_link->v[i]      == n1) dir_v1 = i;
      else if(n_link->v[i] == n2) dir_v2 = i;
      else                        dir_v0 = i;
    }

  l0 = Duplicate_Scalar_Dbl(n_link->b[dir_v0]->l);
  v0 = Duplicate_Scalar_Dbl(n_link->b[dir_v0]->l_var);

  if(n_link->v[dir_v1]->num > n_link->v[dir_v2]->num)
    {
      l1 = Duplicate_Scalar_Dbl(n_link->b[dir_v2]->l);
      v1 = Duplicate_Scalar_Dbl(n_link->b[dir_v2]->l_var);
      l2 = Duplicate_Scalar_Dbl(n_link->b[dir_v1]->l);
      v2 = Duplicate_Scalar_Dbl(n_link->b[dir_v1]->l_var);
    }
  else
    {
      l1 = Duplicate_Scalar_Dbl(n_link->b[dir_v1]->l);
      v1 = Duplicate_Scalar_Dbl(n_link->b[dir_v1]->l_var);
      l2 = Duplicate_Scalar_Dbl(n_link->b[dir_v2]->l);
      v2 = Duplicate_Scalar_Dbl(n_link->b[dir_v2]->l_var);
    }

  for(i=0;i<tree->depth_curr_path+1;++i) move->path[i] = tree->curr_path[i];

  if(move->l0 != NULL)
    {
      Free_Scalar_Dbl(move->l0);
      Free_Scalar_Dbl(move->v0);
    }

  if(move->l1 != NULL)
    {
      Free_Scalar_Dbl(move->l1);
      Free_Scalar_Dbl(move->v1);
    }

  if(move->l2 != NULL)
    {
      Free_Scalar_Dbl(move->l2);
      Free_Scalar_Dbl(move->v2);
    }

  move->l0 = l0;
  move->v0 = v0;

  move->l1 = l1;
  move->v1 = v1;

  move->l2 = l2;
  move->v2 = v2;

  move->depth_path    = tree->depth_curr_path;
  move->pars          = tree->c_pars;
  move->lnL           = tree->c_lnL;
  move->b_target      = b_target;
  move->n_link        = n_link;
  move->b_opp_to_link = b_arrow;
  move->b_init_target = init_target;
  move->dist          = b_target->topo_dist_btw_edges;
  move->n_opp_to_link = (n_link==b_arrow->left)?(b_arrow->rght):(b_arrow->left);

  /* PhyML_Printf("\n. MOVE TEST_ONE_SPR_TARGET : %f",move->lnL); */

  /* for(int i=0;i<tree->size_spr_list_one_edge;++i) */
  /*   { */
  /*     PhyML_Printf("\n. MOVE : %p LNL: %f",tree->spr_list_one_edge[i],tree->spr_list_one_edge[i]->lnL); */
  /*   } */
  rk = Include_One_Spr_To_List_Of_Spr(tree->spr_list_one_edge,tree->size_spr_list_one_edge,move,tree);
  /* PhyML_Printf("\n. RK TEST_ONE_SPR_TARGET : %d",rk); */
  Include_One_Spr_To_List_Of_Spr(tree->spr_list_all_edge,tree->size_spr_list_all_edge,move,tree);
  
  Prune_Subtree(n_link,
                (n_link==b_arrow->left)?(b_arrow->rght):(b_arrow->left),
                &b_target,
                &b_residual,
                tree);

  
  Copy_Scalar_Dbl(init_target_l,b_target->l);
  Copy_Scalar_Dbl(init_target_v,b_target->l_var);

  Copy_Scalar_Dbl(init_arrow_l,b_arrow->l);
  Copy_Scalar_Dbl(init_arrow_v,b_arrow->l_var);

  Copy_Scalar_Dbl(init_residual_l,b_residual->l);
  Copy_Scalar_Dbl(init_residual_v,b_residual->l_var);

  if(tree->mod->s_opt->spr_lnL == YES) Update_PMat_At_Given_Edge(b_target,tree);
  
  tree->c_lnL   = init_lnL;
  tree->c_pars  = init_pars;

  Free_Scalar_Dbl(init_target_l);
  Free_Scalar_Dbl(init_arrow_l);
  Free_Scalar_Dbl(init_residual_l);
  Free_Scalar_Dbl(init_target_v);
  Free_Scalar_Dbl(init_arrow_v);
  Free_Scalar_Dbl(init_residual_v);

  return tree->spr_list_one_edge[rk];
}

/*********************************************************/

void Global_Spr_Search(t_tree *tree)
{
  unsigned int i,iter;
  phydbl best_lnL;
  t_tree *best_tree;
  time_t t_cur;
  phydbl max_delta_lnL_spr,tune_l_mult;
  
  unsigned int no_improv  = 0;
  unsigned int last_best_found = 0;
  unsigned int freq = 1;
  const unsigned int round_freq = 10;
  
  best_lnL      = UNLIKELY;
  tree->verbose = (tree->verbose == VL0) ? VL0 : VL1;
  best_tree     = Duplicate_Tree(tree);
    
  if(tree->io->print_json_trace == YES) JSON_Tree_Io(tree,tree->io->fp_out_json_trace);

  if(tree->verbose > VL0 && tree->io->quiet == NO) PhyML_Printf("\n\n. Score of initial tree: %.2f",tree->c_lnL);
                 
  tree->mod->s_opt->min_diff_lk_move  = 1.E-1;
  tree->mod->s_opt->min_diff_lk_local = 1.E-1;
  Round_Optimize(tree,1000);
  best_lnL = tree->c_lnL;
  Copy_Tree(tree,best_tree);
  
  if(tree->verbose > VL0 && tree->io->quiet == NO) PhyML_Printf("\n\n. Starting first round of SPRs...\n");

  tree->mod->s_opt->max_depth_path            = MAX(40,1+(int)(tree->n_otu/4));
  tree->mod->s_opt->max_delta_lnL_spr         = 1000;
  tree->mod->s_opt->l_min_spr                 = 0.0;
  tree->mod->s_opt->spr_lnL                   = YES;
  tree->mod->s_opt->spr_pars                  = NO;
  tree->mod->s_opt->min_diff_lk_move          = 1.E-0;
  tree->mod->s_opt->min_diff_lk_local         = 1.E-0;
  tree->perform_spr_right_away                = YES;
  tree->mod->s_opt->eval_list_regraft         = NO;
  tree->mod->s_opt->max_delta_lnL_spr_current = 0.0;
  tree->mod->s_opt->min_n_triple_moves        = 1;
  max_delta_lnL_spr                           = 0.0;
  tune_l_mult                                 = 0.01;
  best_lnL                                    = tree->c_lnL;
  
  iter = 0;
  do
    {
      if(!(iter%freq))
        {
          for(int i=0;i<2*tree->n_otu-3;++i) MIXT_Multiply_Scalar_Dbl(tree->a_edges[i]->l,Rgamma((phydbl)(tune_l_mult*iter+1),(phydbl)(1./(tune_l_mult*iter+1))));
          for(int i=0;i<2*tree->n_otu-3;++i) Set_Scalar_Dbl_Min_Thresh(tree->mod->l_min,tree->a_edges[i]->l);
          for(int i=0;i<2*tree->n_otu-3;++i) Set_Scalar_Dbl_Max_Thresh(tree->mod->l_max,tree->a_edges[i]->l);
          freq--;
          if(freq < 1) freq=1;
        }
      
      if(iter> 10) tree->mod->s_opt->l_min_spr = 1.E-3;

      Spr(tree->c_lnL,1.0,tree);
      Optimize_Br_Len_Serie(2,tree);
      
      if(!(iter%round_freq) && iter > 0) Round_Optimize(tree,1000);

      if(tree->verbose > VL0 && tree->io->quiet == NO)
        {
          time(&t_cur);
          PhyML_Printf("\n\t%8ds | %3d | lnL=%12.1f | depth=%5d/%5d | improvements=%4d | delta_lnL=%7.1f/%7.1f %c",
                       (int)(t_cur-tree->t_beg),
                       iter+1,
                       tree->c_lnL,
                       tree->mod->s_opt->max_spr_depth,
                       tree->mod->s_opt->max_depth_path,
                       tree->mod->s_opt->n_improvements,
                       tree->mod->s_opt->max_delta_lnL_spr_current,
                       tree->mod->s_opt->max_delta_lnL_spr,
                       (tree->numerical_warning == YES) ? '!' : ' ');
        }

      if(tree->mod->s_opt->n_improvements > (int)(tree->n_otu/5))
        {
          tune_l_mult *= 1.2;
        }
      
      tree->mod->s_opt->max_depth_path = MAX(5,MAX(tree->mod->s_opt->max_spr_depth+6,(int)(0.6*tree->mod->s_opt->max_depth_path)));
       
      if((iter%4) > 0 || iter == 0)
        {
          if(tree->mod->s_opt->max_delta_lnL_spr_current > max_delta_lnL_spr) max_delta_lnL_spr = tree->mod->s_opt->max_delta_lnL_spr_current;
        }
      else if(iter > 0)
        {
          tree->mod->s_opt->max_delta_lnL_spr = MAX(50.,2.*max_delta_lnL_spr);
          max_delta_lnL_spr = 0.0;
        }
      
      
      if(tree->c_lnL > best_lnL)
        {
          no_improv = 0;
          best_lnL = tree->c_lnL;
          Copy_Tree(tree,best_tree);
          if(tree->verbose > VL0 && tree->io->quiet == NO) PhyML_Printf(" +");
          if(tree->io->print_json_trace == YES) JSON_Tree_Io(tree,tree->io->fp_out_json_trace);
        }
            
      no_improv++;
      iter++;
    }
  while(tree->mod->s_opt->n_improvements > 15 && no_improv < 10);
        

  if(tree->verbose > VL0 && tree->io->quiet == NO) PhyML_Printf("\n\n. Second round of optimization...\n");

  
  tree->mod->s_opt->max_depth_path            = MAX(10,(int)(1.5*tree->mod->s_opt->max_depth_path));
  tree->mod->s_opt->l_min_spr                 = 1.E-3;
  tree->mod->s_opt->spr_lnL                   = YES;
  tree->mod->s_opt->spr_pars                  = NO;
  tree->mod->s_opt->min_diff_lk_move          = 1.E-2;
  tree->mod->s_opt->min_diff_lk_local         = 1.E-2;
  tree->mod->s_opt->eval_list_regraft         = NO;
  tree->mod->s_opt->max_delta_lnL_spr_current = 0.0;
  tree->mod->s_opt->min_n_triple_moves        = 1;
  tree->mod->s_opt->deepest_path              = 0;
  max_delta_lnL_spr                           = 0.0;
  no_improv                                   = 0;
  
  do
    {
      if(!(iter%freq))
        {
          for(int i=0;i<2*tree->n_otu-3;++i) MIXT_Multiply_Scalar_Dbl(tree->a_edges[i]->l,Rgamma((phydbl)(tune_l_mult*iter+1),(phydbl)(1./(tune_l_mult*iter+1))));
          for(int i=0;i<2*tree->n_otu-3;++i) Set_Scalar_Dbl_Min_Thresh(tree->mod->l_min,tree->a_edges[i]->l);
          for(int i=0;i<2*tree->n_otu-3;++i) Set_Scalar_Dbl_Max_Thresh(tree->mod->l_max,tree->a_edges[i]->l);
          freq--;
          if(freq < 1) freq=1;
        }
      
      Spr(tree->c_lnL,1.0,tree);
      Optimize_Br_Len_Serie(2,tree);
      
      if(!(iter%round_freq) && iter > 0) Round_Optimize(tree,1000);
      
      if(tree->verbose > VL0 && tree->io->quiet == NO)
        {
          time(&t_cur);
          PhyML_Printf("\n\t%8ds | %3d | lnL=%12.1f | depth=%5d/%5d | improvements=%4d | delta_lnL=%7.1f/%7.1f %c",
                       (int)(t_cur-tree->t_beg),
                       iter+1,
                       tree->c_lnL,
                       tree->mod->s_opt->max_spr_depth,
                       tree->mod->s_opt->max_depth_path,
                       tree->mod->s_opt->n_improvements,
                       tree->mod->s_opt->max_delta_lnL_spr_current,
                       tree->mod->s_opt->max_delta_lnL_spr,
                       (tree->numerical_warning == YES) ? '!' : ' ');
        }

      tree->mod->s_opt->max_depth_path = MAX(5,MAX(tree->mod->s_opt->max_spr_depth+4,(int)(0.8*tree->mod->s_opt->max_depth_path)));
      tree->mod->s_opt->max_depth_path = MIN(20,tree->mod->s_opt->max_depth_path);
      
      if((iter%4) > 0 || iter == 0)
        {
          if(tree->mod->s_opt->max_delta_lnL_spr_current > max_delta_lnL_spr) max_delta_lnL_spr = tree->mod->s_opt->max_delta_lnL_spr_current;
        }
      else if(iter > 0)
        {
          tree->mod->s_opt->max_delta_lnL_spr = MAX(20.,2.*max_delta_lnL_spr);
          max_delta_lnL_spr = 0.0;
        }
      
      
      if(tree->c_lnL > best_lnL)
        {
          no_improv = 0;
          best_lnL = tree->c_lnL;
          Copy_Tree(tree,best_tree);
          if(tree->verbose > VL0 && tree->io->quiet == NO) PhyML_Printf(" +");
          if(tree->io->print_json_trace == YES) JSON_Tree_Io(tree,tree->io->fp_out_json_trace);
        }
      
      no_improv++;
      iter++;
    }
  while(tree->mod->s_opt->n_improvements > 5 && no_improv < 10);

  if(tree->verbose > VL0 && tree->io->quiet == NO) PhyML_Printf("\n\n. Third round of optimization...\n");
  last_best_found = 0;

  tree->mod->s_opt->max_depth_path            = MAX(10,tree->mod->s_opt->max_depth_path);
  tree->mod->s_opt->max_delta_lnL_spr         = MAX(100.,tree->mod->s_opt->max_delta_lnL_spr);
  tree->mod->s_opt->spr_lnL                   = YES;
  tree->mod->s_opt->spr_pars                  = NO;
  tree->mod->s_opt->l_min_spr                 = 1.E-4;
  tree->mod->s_opt->min_diff_lk_move          = 1.E-2;
  tree->mod->s_opt->min_diff_lk_local         = 1.E-2;
  tree->mod->s_opt->eval_list_regraft         = YES;
  tree->mod->s_opt->max_delta_lnL_spr_current = 0.0;
  tree->mod->s_opt->min_n_triple_moves        = 5;
  tree->mod->s_opt->max_rank_triple_move      = 0;
  tree->mod->s_opt->max_no_better_tree_found  = 10;
  
  do
    {
      Spr(tree->c_lnL,1.0,tree);
      Optimize_Br_Len_Serie(2,tree);
            
      if(!(iter%round_freq)) Round_Optimize(tree,1000);

      if(tree->verbose > VL0 && tree->io->quiet == NO)
        {
          time(&t_cur);
          PhyML_Printf("\n\t%8ds | %3d | lnL=%12.1f | depth=%5d/%5d | improvements=%4d | delta_lnL=%7.1f/%7.1f | triple moves=%4d %c",
                       (int)(t_cur-tree->t_beg),
                       iter+1,
                       tree->c_lnL,
                       tree->mod->s_opt->max_spr_depth,
                       tree->mod->s_opt->max_depth_path,
                       tree->mod->s_opt->n_improvements,
                       tree->mod->s_opt->max_delta_lnL_spr_current,
                       tree->mod->s_opt->max_delta_lnL_spr,
                       tree->mod->s_opt->min_n_triple_moves,
                       (tree->numerical_warning == YES) ? '!' : ' ');
        }
      
      tree->mod->s_opt->max_depth_path     = MIN(30,MAX(5,MAX(2*tree->mod->s_opt->max_spr_depth,(int)(0.8*tree->mod->s_opt->max_depth_path))));
      tree->mod->s_opt->max_delta_lnL_spr  = MAX(100.,0.8*tree->mod->s_opt->max_delta_lnL_spr_current);
      tree->mod->s_opt->min_n_triple_moves = MAX(5,2*tree->mod->s_opt->max_rank_triple_move);
      
      if(tree->c_lnL > best_lnL)
        {
          last_best_found = 0;;
          best_lnL = tree->c_lnL;
          Copy_Tree(tree,best_tree);
          if(tree->verbose > VL0 && tree->io->quiet == NO) PhyML_Printf(" +");
          if(tree->io->print_json_trace == YES) JSON_Tree_Io(tree,tree->io->fp_out_json_trace);
        }
      
      last_best_found++;
      iter++;
    }
  while(tree->mod->s_opt->n_improvements > 0 && last_best_found <= tree->mod->s_opt->max_no_better_tree_found);

  Copy_Tree(best_tree,tree);
  
  if(tree->verbose > VL0 && tree->io->quiet == NO) PhyML_Printf("\n\n. Final optimisation steps...\n");
  
  i = tree->verbose;
  tree->verbose = VL0;
  tree->mod->s_opt->min_diff_lk_move  = 1.E-03;
  tree->mod->s_opt->min_diff_lk_local = 1.E-03;
  do
    {
      tree->mod->s_opt->fast_nni = NO;
      Round_Optimize(tree,10);
      if(!Check_NNI_Five_Branches(tree)) break;
    }
  while(1);
  tree->verbose = i;
  
  Free_Tree(best_tree);

  return;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Speed_Spr(t_tree *tree, phydbl prop_spr, int max_cycles, phydbl delta_lnL)
{
  int step,old_pars;
  phydbl old_lnL;

  if(tree->lock_topo == YES)
    {
      PhyML_Fprintf(stderr,"\n== The tree topology is locked.");
      PhyML_Fprintf(stderr,"\n== Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
      Exit("\n");
    }


  Set_Both_Sides(NO,tree);
  Pars(NULL,tree);
  if(tree->mod->s_opt->spr_pars == NO) Lk(NULL,tree);
  Record_Br_Len(tree);
 
  tree->mod->s_opt->deepest_path  = 0;
  tree->best_pars                 = tree->c_pars;
  tree->best_lnL                  = tree->c_lnL;
  old_lnL                         = tree->c_lnL;
  old_pars                        = tree->c_pars;
  step                            = 0;
  do
    {
      ++step;

      old_lnL  = tree->c_lnL;
      old_pars = tree->c_pars;

      Set_Both_Sides(YES,tree);
      Pars(NULL,tree);
      if(tree->mod->s_opt->spr_pars == NO) Lk(NULL,tree);
      Spr(UNLIKELY,prop_spr,tree);

      // Set maximum depth for future spr rounds to deepest spr found so far
      tree->mod->s_opt->max_depth_path = tree->mod->s_opt->max_spr_depth;

      if(tree->mod->s_opt->spr_pars == NO)
        {
          if(tree->mod->s_opt->n_improvements > 0)
            {
              /* Optimise branch lengths */
              Optimize_Br_Len_Serie(2,tree);
              /* Print log-likelihood and parsimony scores */
              if(tree->verbose > VL2 && tree->io->quiet == NO) Print_Lk(tree,"[Branch lengths     ]");
            }
        }
      else
        {
          if(tree->verbose > VL2 && tree->io->quiet == NO) Print_Pars(tree);
        }

      Pars(NULL,tree);
      
      if(tree->io->print_trace)
        {
          char *s = Write_Tree(tree);
          PhyML_Fprintf(tree->io->fp_out_trace,"[%f]%s\n",tree->c_lnL,s); fflush(tree->io->fp_out_trace);
          if((tree->io->print_site_lnl) && (!tree->mod->s_opt->spr_pars))
            {
              Print_Site_Lk(tree,tree->io->fp_out_lk);
              fflush(tree->io->fp_out_lk);
            }
          Free(s);
        }

      if(tree->io->print_json_trace == YES) JSON_Tree_Io(tree,tree->io->fp_out_json_trace); 
      
      /* Record the current best log-likelihood and parsimony */
      if(tree->c_lnL  > tree->best_lnL)  tree->best_lnL  = tree->c_lnL; 
      if(tree->c_pars < tree->best_pars) tree->best_pars = tree->c_pars;

      if(tree->mod->s_opt->spr_pars == NO)
        {
          if(tree->c_lnL < old_lnL-tree->mod->s_opt->min_diff_lk_local)
            {
              PhyML_Fprintf(stderr,"\n== old_lnL = %f c_lnL = %f",old_lnL,tree->c_lnL);
              PhyML_Fprintf(stderr,"\n== Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
              Exit("");
            }
        }
      else
        {
          if(tree->c_pars > old_pars)
            {
              PhyML_Fprintf(stderr,"\n== old_pars = %d c_pars = %d",old_pars,tree->c_pars);
              PhyML_Fprintf(stderr,"\n== Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
              Exit("");
            }
        }

      /* Record the current best branch lengths  */
      Record_Br_Len(tree);

      /* Exit if no improvements after complete optimization */
      if(step+1 > max_cycles) break;
      if((tree->mod->s_opt->spr_pars == NO)  && (fabs(old_lnL-tree->c_lnL)   < delta_lnL)) break;
      if((tree->mod->s_opt->spr_pars == YES) && (fabs((phydbl)(old_pars-tree->c_pars)) < 1)) break;
      if(!tree->mod->s_opt->n_improvements) break;
    }
  while(1);
}

/*********************************************************/

int Evaluate_List_Of_Regraft_Pos_Triple(t_spr **spr_list, int list_size, t_tree *tree)
{
  t_spr *move;
  t_edge *init_target, *b_residual;
  int i,j,best_move,better_found;
  int dir_v0, dir_v1, dir_v2;
  scalar_dbl *recorded_l,*recorded_v;
  phydbl best_lnL,init_lnL;
  int recorded;

  if(tree->mixt_tree != NULL)
    {
      PhyML_Fprintf(stderr,"\n== Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
      Exit("\n");
    }

  best_lnL = UNLIKELY;
  init_target = b_residual = NULL;
  best_move = -1;
  init_lnL = tree->c_lnL;
  recorded_v = recorded_l = NULL;
  better_found = NO;

  if(list_size == 0)
    {
      PhyML_Fprintf(stderr,"\n== List size is 0 !");
      PhyML_Fprintf(stderr,"\n== Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
      Exit("\n");
    }

  
  recorded = NO;
  for(i=0;i<list_size;i++)
    {
      move = spr_list[i];

      if(!move)
        {
          PhyML_Fprintf(stderr,"\n== move is NULL\n");
          PhyML_Fprintf(stderr,"\n== Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
          Exit("\n");
        }
      
      if(move->b_target)
        {
          /* Record t_edge lengths */
          Record_Br_Len(tree);
          
          /* Prune subtree */
          Prune_Subtree(move->n_link,
                        move->n_opp_to_link,
                        &init_target,
                        &b_residual,
                        tree);
          
          if(recorded == NO)
            {
              /*! Rough optimisation of the branch length at prune site
               *  We only need to perform this optimisation for the first
               *  element of spr_list because the pruned subtree is the
               *  same across all the elements of spr_list. It would not
               *  be true in the general case
               */
              recorded = YES;
              
              Br_Len_Opt(&(init_target->l->v),init_target,tree);

              /*! Record branch length at prune site */
              if(recorded_l == NULL)
                {
                  recorded_l = Duplicate_Scalar_Dbl(init_target->l);
                  recorded_v = Duplicate_Scalar_Dbl(init_target->l_var);
                }
              else
                {
                  Copy_Scalar_Dbl(init_target->l,recorded_l);
                  Copy_Scalar_Dbl(init_target->l_var,recorded_v);
                }

              Copy_Scalar_Dbl(recorded_l,move->init_target_l);
              Copy_Scalar_Dbl(recorded_v,move->init_target_v);
            }
          else
            {
              Copy_Scalar_Dbl(recorded_l,move->b_init_target->l);
              Copy_Scalar_Dbl(recorded_v,move->b_init_target->l_var);

              Copy_Scalar_Dbl(recorded_l,move->init_target_l);
              Copy_Scalar_Dbl(recorded_v,move->init_target_v);
            }
          
          /* Update the change proba matrix at prune position */
          Update_PMat_At_Given_Edge(init_target,tree);
          
          /* Update conditional likelihoods along the path from the prune to
             the regraft position */
          MIXT_Set_Alias_Subpatt(YES,tree);
          Update_Partial_Lk_Along_A_Path(move->path,move->depth_path+1,tree);
          MIXT_Set_Alias_Subpatt(NO,tree);

          
          /* Regraft subtree */
          Graft_Subtree(move->b_target,move->n_link,NULL,b_residual,NULL,tree);
                    
          MIXT_Set_Alias_Subpatt(YES,tree);
          move->lnL = Triple_Dist(move->n_link,tree);
          MIXT_Set_Alias_Subpatt(NO,tree);

          /* PhyML_Printf("\n. %d/%d init_lnL: %12G move->lnL= %12G best_lnL=%12G absolute_best=%12G", */
          /*              i, */
          /*              list_size, */
          /*              init_lnL, */
          /*              move->lnL, */
          /*              best_lnL, */
          /*              tree->best_lnL); */
          
          /* Record updated branch lengths for this move */
          dir_v1 = dir_v2 = dir_v0 = -1;
          for(j=0;j<3;j++)
            {
              if(move->n_link->v[j] == move->n_opp_to_link) dir_v0 = j;
              else if(dir_v1 < 0)                           dir_v1 = j;
              else                                          dir_v2 = j;
            }

          Copy_Scalar_Dbl(move->n_link->b[dir_v0]->l,    move->l0);
          Copy_Scalar_Dbl(move->n_link->b[dir_v0]->l_var,move->v0);

          if(move->n_link->v[dir_v1]->num > move->n_link->v[dir_v2]->num)
            {
              Copy_Scalar_Dbl(move->n_link->b[dir_v2]->l,    move->l1);
              Copy_Scalar_Dbl(move->n_link->b[dir_v2]->l_var,move->v1);

              Copy_Scalar_Dbl(move->n_link->b[dir_v1]->l,    move->l2);
              Copy_Scalar_Dbl(move->n_link->b[dir_v1]->l_var,move->v2);
            }
          else
            {
              Copy_Scalar_Dbl(move->n_link->b[dir_v1]->l,    move->l1);
              Copy_Scalar_Dbl(move->n_link->b[dir_v1]->l_var,move->v1);

              Copy_Scalar_Dbl(move->n_link->b[dir_v2]->l,    move->l2);
              Copy_Scalar_Dbl(move->n_link->b[dir_v2]->l_var,move->v2);
            }
          
          if(move->lnL > best_lnL)
            {
              best_lnL  = move->lnL;
              best_move = i;
            }
          
          /* Regraft the subtree at its original position */
          Prune_Subtree(move->n_link,
                        move->n_opp_to_link,
                        &move->b_target,
                        &b_residual,
                        tree);
          
          Graft_Subtree(init_target,
                        move->n_link,
                        NULL,
                        b_residual,
                        NULL,
                        tree);
          
          /* Restore branch lengths */
          Restore_Br_Len(tree);
          
          /* Update relevant change proba matrices */
          Update_PMat_At_Given_Edge(move->b_target,tree);
                    
          tree->c_lnL = init_lnL;
        }
      
      /* PhyML_Printf("\n. [ %4d/%4d ] %f %f %s", */
      /*              i,list_size,tree->best_lnL,move->lnL, */
      /*              (move->lnL > tree->best_lnL + tree->mod->s_opt->min_diff_lk_move) ? "**" : ""); */


      /* Bail out as soon as you've found a true improvement */
      if(move->lnL > tree->best_lnL + tree->mod->s_opt->min_diff_lk_move) 
        {
          better_found = YES;
          if(i > tree->mod->s_opt->max_rank_triple_move) tree->mod->s_opt->max_rank_triple_move = i; 
          break;      
        }
    }
  
  /*   PhyML_Printf("\n. max_improv = %f",max_improv); */
  
  
  if(better_found == NO)
    {
      MIXT_Set_Alias_Subpatt(YES,tree);
      for(i=0;i<list_size;i++)
        {
          move = spr_list[i];
          if(move->b_target)
            {
              for(j=0;j<3;++j) Update_PMat_At_Given_Edge(move->n_link->b[j],tree);
              for(j=0;j<3;++j) Update_Partial_Lk(tree,move->n_link->b[j],move->n_link);
              
              /* TO DO : we don't need to update all these partial likelihoods here.
                 Would need to record only those that were along the paths examined
                 above */
              
              for(j=0;j<3;++j)
                if(move->n_link->v[j] != move->n_opp_to_link)
                  Pre_Order_Lk(move->n_link,move->n_link->v[j],tree);
              
              break;
            }
        }
      MIXT_Set_Alias_Subpatt(NO,tree);
    }

#ifdef DEBUG
  if(best_move < 0 && list_size > 0)
    {
      PhyML_Printf("\n\n== Best_move < 0 !");
      PhyML_Printf("\n== List size = %d",list_size);
      PhyML_Printf("\n== Best lnL = %f",best_lnL);
      for(i=0;i<list_size;i++)
        {
          move = spr_list[i];
          PhyML_Printf("\n== move %p %p lnL: %f",move,move->b_target,move->lnL);
        }
      
      PhyML_Fprintf(stderr,"\n== Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
      Exit("\n");
    }
#endif
  
  Free_Scalar_Dbl(recorded_l);
  Free_Scalar_Dbl(recorded_v);
  
  return best_move;
}

/*********************************************************/

int Try_One_Spr_Move_Triple(t_spr *move, t_tree *tree)
{
  t_edge *init_target, *b_residual;
  int j;
  int dir_v0, dir_v1, dir_v2;
  int accept;
  
  assert(move);
  if(move->n_link == NULL) return -1;
  
  if(tree->mixt_tree != NULL)
    {
      PhyML_Fprintf(stderr,"\n== Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
      Exit("\n");
    }

  Record_Br_Len(tree);

  Prune_Subtree(move->n_link,
                move->n_opp_to_link,
                &init_target,
                &b_residual,
                tree);
    
  Copy_Scalar_Dbl(move->init_target_l,init_target->l);
  Copy_Scalar_Dbl(move->init_target_v,init_target->l_var);

  Graft_Subtree(move->b_target,move->n_link,NULL,b_residual,NULL,tree);

  dir_v1 = dir_v2 = dir_v0 = -1;
  for(j=0;j<3;j++)
    {
      if(move->n_link->v[j] == move->n_opp_to_link) dir_v0 = j;
      else if(dir_v1 < 0)                           dir_v1 = j;
      else                                          dir_v2 = j;
    }

  Copy_Scalar_Dbl(move->l0,move->n_link->b[dir_v0]->l);
  Copy_Scalar_Dbl(move->v0,move->n_link->b[dir_v0]->l_var);

  if(move->n_link->v[dir_v1]->num > move->n_link->v[dir_v2]->num)
    {
      Copy_Scalar_Dbl(move->l1,move->n_link->b[dir_v2]->l);
      Copy_Scalar_Dbl(move->v1,move->n_link->b[dir_v2]->l_var);

      Copy_Scalar_Dbl(move->l2,move->n_link->b[dir_v1]->l);
      Copy_Scalar_Dbl(move->v2,move->n_link->b[dir_v1]->l_var);
    }
  else
    {
      Copy_Scalar_Dbl(move->l1,move->n_link->b[dir_v1]->l);
      Copy_Scalar_Dbl(move->v1,move->n_link->b[dir_v1]->l_var);

      Copy_Scalar_Dbl(move->l2,move->n_link->b[dir_v2]->l);
      Copy_Scalar_Dbl(move->v2,move->n_link->b[dir_v2]->l_var);
    }

  accept = YES;
  if(!Check_Topo_Constraints(tree,tree->io->cstr_tree)) accept = NO;

  
  if(accept == YES) /* Apply the move */
    {
      time(&(tree->t_current));
      Pars(NULL,tree);

      Update_PMat_At_Given_Edge(init_target,tree);
      for(int i=0;i<3;++i) Update_PMat_At_Given_Edge(move->n_link->b[i],tree);
      Post_Order_Lk(move->n_opp_to_link,move->n_link,tree);
      Pre_Order_Lk(move->n_opp_to_link,move->n_link,tree);
      Pre_Order_Lk(move->n_link,move->n_opp_to_link,tree);      
      Lk(move->b_opp_to_link,tree);

      if(fabs(tree->c_lnL - move->lnL) > tree->mod->s_opt->min_diff_lk_move)
        {
          PhyML_Fprintf(stderr,"\n== c_lnL = %f move_lnL = %f", tree->c_lnL,move->lnL);
          PhyML_Fprintf(stderr,"\n== %d l0=%G l1=%G l2=%G v0=%G v1=%G v2=%G",move->n_link->num,move->l0->v,move->l1->v,move->l2->v,move->v0->v,move->v1->v,move->v2->v);
          PhyML_Fprintf(stderr,"\n== Gamma MGF? %d",tree->io->mod->gamma_mgf_bl);
          PhyML_Fprintf(stderr,"\n== Err. in file %s at line %d.\n",__FILE__,__LINE__);
          Check_Lk_At_Given_Edge(YES,tree);
          Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
        }

      if(tree->verbose > VL2 && tree->io->quiet == NO)
        {
          Print_Lk_And_Pars(tree);
          PhyML_Printf(" [depth=%5d]",move->depth_path); fflush(NULL);
        }


      tree->mod->s_opt->n_improvements++;

      t_spr *dum_move = move;
      phydbl delta = 0.0;
      while(dum_move)
        {
          delta = move->lnL - dum_move->lnL;
          if(delta > tree->mod->s_opt->max_delta_lnL_spr_current)
            tree->mod->s_opt->max_delta_lnL_spr_current = delta;
          dum_move = dum_move->path_prev;
        }

      
      if(tree->c_lnL > tree->best_lnL) tree->best_lnL = tree->c_lnL;
      Record_Br_Len(tree);

      if(move->depth_path > tree->mod->s_opt->deepest_path)
        tree->mod->s_opt->deepest_path = move->depth_path;

      if(move->depth_path > tree->mod->s_opt->max_spr_depth) tree->mod->s_opt->max_spr_depth = move->depth_path;
        
      return 1;
    }
  else // Go back to original topology
    {
      Prune_Subtree(move->n_link,
                    move->n_opp_to_link,
                    &move->b_target,
                    &b_residual,
                    tree);
      
      Graft_Subtree(init_target,
                    move->n_link,
                    NULL,
                    b_residual,
                    NULL,
                    tree);
  
      Restore_Br_Len(tree);

      return 0;
    }
}

/*********************************************************/

int Try_One_Spr_Move_Full(t_spr *move, short int apply_move, t_tree *tree)
{
  t_edge *init_target, *b_residual;
  phydbl init_lnL;
  
  assert(move);
  if(move->n_link == NULL) return -1;
  init_lnL = tree->c_lnL;
  
  Record_Br_Len(tree);

  Prune_Subtree(move->n_link,
                move->n_opp_to_link,
                &init_target,
                &b_residual,
                tree);
  Graft_Subtree(move->b_target,move->n_link,NULL,b_residual,NULL,tree);

  Optimize_Br_Len_Serie(2,tree);

  move->lnL = tree->c_lnL;

  
  if(tree->c_lnL > tree->best_lnL + tree->mod->s_opt->min_diff_lk_move && apply_move == YES)
    {
      tree->best_lnL = tree->c_lnL;
      tree->mod->s_opt->n_improvements++;
      return 1;
    }
  else
    {
      Prune_Subtree(move->n_link,
            move->n_opp_to_link,
            &move->b_target,
            &b_residual,
            tree);

      Graft_Subtree(init_target,
                    move->n_link,
                    NULL,
                    b_residual,
                    NULL,
                    tree);

      Restore_Br_Len(tree);
      tree->c_lnL = init_lnL;
      return 0;
    }

  return -1;
}

/*********************************************************/

unsigned int Include_One_Spr_To_List_Of_Spr(t_spr **list, int list_size, t_spr *move, t_tree *tree)
{
  unsigned int i, rk;
  t_spr *buff_spr,*orig_move, *orig_move_list, *move_list;
  t_tree *orig_tree;

  if(tree->mixt_tree != NULL)
    {
      PhyML_Fprintf(stderr,"\n== Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
      Exit("\n");
    }

  rk = 0;
  
  /* printf("\n. INCLUDE: %d %G %G",list_size,move->lnL,list[list_size-1]->lnL); */
  
  if(((tree->mod->s_opt->spr_lnL == YES) && (move->lnL  > list[list_size-1]->lnL)) ||
     ((tree->mod->s_opt->spr_lnL == NO) && (move->pars <= list[list_size-1]->pars)))
    {
      move_list = list[list_size-1];

      /* printf("\n. Include move with lnL: %G to list %p",move->lnL,(void *)list); */
      
      move_list->depth_path    = move->depth_path;
      move_list->pars          = move->pars;
      move_list->lnL           = move->lnL;
      move_list->dist          = move->dist;

      if(move_list->l0 == NULL)
        {
          move_list->l0 = Duplicate_Scalar_Dbl(move->l0);
          move_list->v0 = Duplicate_Scalar_Dbl(move->v0);
        }
      else
        {
          Copy_Scalar_Dbl(move->l0,move_list->l0);
          Copy_Scalar_Dbl(move->v0,move_list->v0);
        }

      if(move_list->l1 == NULL)
        {
          move_list->l1 = Duplicate_Scalar_Dbl(move->l1);
          move_list->v1 = Duplicate_Scalar_Dbl(move->v1);
        }
      else
        {
          Copy_Scalar_Dbl(move->l1,move_list->l1);
          Copy_Scalar_Dbl(move->v1,move_list->v1);
        }

      if(move_list->l2 == NULL)
        {
          move_list->l2 = Duplicate_Scalar_Dbl(move->l2);
          move_list->v2 = Duplicate_Scalar_Dbl(move->v2);
        }
      else
        {
          Copy_Scalar_Dbl(move->l2,move_list->l2);
          Copy_Scalar_Dbl(move->v2,move_list->v2);
        }


      if(move_list->init_target_l == NULL)
        {
          move_list->init_target_l = Duplicate_Scalar_Dbl(move->init_target_l);     
          move_list->init_target_v = Duplicate_Scalar_Dbl(move->init_target_v);     
        }
      else
        { 
          Copy_Scalar_Dbl(move->init_target_l,move_list->init_target_l);     
          Copy_Scalar_Dbl(move->init_target_v,move_list->init_target_v);     
        }

      orig_move      = move;
      orig_move_list = move_list;
      orig_tree      = tree;
      do
        {
          move_list->b_target      = move->b_target;
          move_list->n_link        = move->n_link;
          move_list->n_opp_to_link = move->n_opp_to_link;
          move_list->b_opp_to_link = move->b_opp_to_link;
          move_list->b_init_target = move->b_init_target;

          move      = move->next;
          move_list = move_list->next;
          tree      = tree->next;
        }
      while(tree);

      move      = orig_move;
      move_list = orig_move_list;
      tree      = orig_tree;
      
      for(i=0;i<move_list->depth_path+1;++i) move_list->path[i] = move->path[i];
      
      /* PhyML_Printf("\n. LIST_SIZE: %d",list_size); */
      for(i=list_size-1;i>0;i--)
        {
          if(((tree->mod->s_opt->spr_lnL == YES) && (list[i]->lnL > list[i-1]->lnL)) ||
             ((tree->mod->s_opt->spr_lnL == NO) && (list[i]->pars <=  list[i-1]->pars)))
            {  
              buff_spr  = list[i-1];
              list[i-1] = list[i];
              list[i]   = buff_spr;
            }
          else
            {
              rk = i;
              break;
            }
        }
    }
  return rk;
}

/*********************************************************/

void Random_Spr(int n_moves, t_tree *tree)
{
  int i;
  int br_pulled, br_target;
  t_spr *spr_struct;
  t_edge *target, *residual;

  spr_struct = Make_One_Spr(tree);
  Init_One_Spr(spr_struct);
  target = residual = NULL;

  for(i=0;i<n_moves;i++)
    {
      /* br_pulled = (int)((phydbl)rand()/RAND_MAX * (2*tree->n_otu-3-1)); */
      br_pulled = Rand_Int(0,2*tree->n_otu-3-1);
      
      do
        {
          /* br_target = (int)((phydbl)rand()/RAND_MAX * (2*tree->n_otu-3-1)); */
          br_target = Rand_Int(0,2*tree->n_otu-3-1);
        }
      while(br_target == br_pulled);

      spr_struct->n_link        = tree->a_edges[br_pulled]->left;
      spr_struct->n_opp_to_link = tree->a_edges[br_pulled]->rght;
      spr_struct->b_opp_to_link = tree->a_edges[br_pulled];
      spr_struct->b_target      = tree->a_edges[br_target];
      spr_struct->b_init_target = NULL;
      
      if(!Check_Spr_Move_Validity(spr_struct,tree))
        {
          spr_struct->n_link        = tree->a_edges[br_pulled]->rght;
          spr_struct->n_opp_to_link = tree->a_edges[br_pulled]->left;
        }
      
#ifdef DEBUG
      if(!Check_Spr_Move_Validity(spr_struct,tree))
        {
          Warn_And_Exit("\n== Could not find a valid move...\n");
        }
#endif

      if(spr_struct->n_link == spr_struct->b_target->left ||
         spr_struct->n_link == spr_struct->b_target->rght)
        {
          n_moves++;
          continue;
        }

      Prune_Subtree(spr_struct->n_link,
                    spr_struct->n_opp_to_link,
                    &target,
                    &residual,
                    tree);

      Graft_Subtree(spr_struct->b_target,
                    spr_struct->n_link,                    
                    NULL,
                    residual,NULL,tree);
    }
  Free(spr_struct);
}

/*********************************************************/

void Reset_Spr_List(t_spr **list, int size)
{
  int i;

  for(i=0;i<size;++i)
    {
      list[i]->depth_path     = 0;
      list[i]->pars           = MAX_PARS;
      list[i]->lnL            = UNLIKELY;
      list[i]->n_link         = NULL;
      list[i]->n_opp_to_link  = NULL;
      list[i]->b_target       = NULL;
  }
}

/*********************************************************/

int Check_Spr_Move_Validity(t_spr *this_spr_move, t_tree *tree)
{
  int match;

  match = 0;
  Found_In_Subtree(this_spr_move->n_link,
           this_spr_move->n_opp_to_link,
           this_spr_move->b_target->left,
           &match,
           tree);

  if(match) return 0;
  else      return 1;
}

/*********************************************************/

void Spr_Pars(int threshold, int n_round_max, t_tree *tree)
{  
  int curr_pars,round;

  if(tree->verbose > VL2 && tree->io->quiet == NO) PhyML_Printf("\n. Minimizing parsimony...\n");

  tree->best_pars                  = 1E+8;
  tree->best_lnL                   = UNLIKELY;
  tree->mod->s_opt->spr_lnL        = NO;
  tree->mod->s_opt->spr_pars       = YES;
  curr_pars                        = tree->c_pars;
  tree->mod->s_opt->max_depth_path = tree->n_otu;
  round                            = 0;
  do
    {
      curr_pars = tree->c_pars;
      Speed_Spr(tree,1.0,1,0.0);      
    }
  while(tree->mod->s_opt->n_improvements && fabs((phydbl)(tree->c_pars - curr_pars)) > threshold && round++ < n_round_max);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Spr_Shuffle(t_tree *mixt_tree)
{
  phydbl lk_old;
  int *orig_catg,n,n_iter;
  t_tree *tree,**tree_list;

  if(mixt_tree->verbose > VL0) PhyML_Printf("\n\n. Refining the tree...\n");

  /*! Get the number of classes in each mixture */
  orig_catg = MIXT_Get_Number_Of_Classes_In_All_Mixtures(mixt_tree);


  /*! Set the number of rate classes to (at most) 2.
    ! Propagate this to every mixture tree in the analysis
  */
  tree = mixt_tree;
  n = 0;
  do
    {
#ifdef BEAGLE
      tree->b_inst = create_beagle_instance(tree, tree->io->quiet, tree->io);
      //Instead of capping the rate categories at 2, just
      //give the other categories 0 weight
      if(orig_catg[n] > 2) //should we even bother?
      {
          double cat_wghts[orig_catg[n]];
          //Give the first two categories equal weights
          cat_wghts[0] = 0.5;
          cat_wghts[1] = 0.5;
          int i;
          for(i=2;i<orig_catg[n];++i){
            cat_wghts[i] = 0.0;
          }
          int ret = beagleSetCategoryWeights(tree->b_inst,0,cat_wghts);
          if(ret<0) {fprintf(stderr, "beagleSetCategoryWeights() on instance %i failed:%i\n\n",tree->b_inst,ret);Exit(""); }
          tree->mod->optimizing_topology = true;
      }
#endif
      tree->mod->ras->n_catg = MIN(2,orig_catg[n]);
      if(tree->mod->ras->invar == YES) tree->mod->ras->n_catg--;
      tree = tree->next_mixt;
      n++;
    }
  while(tree);


  /*! Make sure the number of trees in each mixture is at most 2
   */
  tree_list = MIXT_Record_All_Mixtures(mixt_tree);
  MIXT_Break_All_Mixtures(orig_catg,mixt_tree);

  Set_Both_Sides(YES,mixt_tree);
  Lk(NULL,mixt_tree);


  mixt_tree->best_pars                     = 1E+8;
  mixt_tree->mod->s_opt->spr_lnL           = NO;
  mixt_tree->mod->s_opt->spr_pars          = NO;
  mixt_tree->mod->s_opt->quickdirty        = NO;
  mixt_tree->best_lnL                      = mixt_tree->c_lnL;
  mixt_tree->mod->s_opt->max_depth_path    = mixt_tree->n_otu;
  mixt_tree->mod->s_opt->min_diff_lk_move  = 0.1;
  mixt_tree->annealing_temp                = 0.0;

  n_iter = 0;
  do
    {
      Set_Both_Sides(YES,mixt_tree);
      Lk(NULL,mixt_tree);
      Pars(NULL,mixt_tree);
      Record_Br_Len(mixt_tree);

      mixt_tree->best_pars = mixt_tree->c_pars;
      mixt_tree->best_lnL  = mixt_tree->c_lnL;

      lk_old = mixt_tree->c_lnL;
      Spr(UNLIKELY,1.0,mixt_tree);

      mixt_tree->annealing_temp -= 2.;

      if(mixt_tree->annealing_temp < 0.0) mixt_tree->annealing_temp = 0.0;

      if(mixt_tree->mod->s_opt->n_improvements < 5      || 
         mixt_tree->mod->s_opt->max_spr_depth  < 2      || 
         FABS(lk_old-mixt_tree->c_lnL) < 5. ||
         ++n_iter > 10) break;
      
    }
  while(1);

  mixt_tree->annealing_temp = 0.0;

  if(mixt_tree->verbose > VL0 && mixt_tree->io->quiet == NO)
    {
      PhyML_Printf("\n\n. End of refining stage...\n");
    }

  /*! Go back to the original data structure, with potentially more
    ! than 2 trees per mixture
   */
  MIXT_Reconnect_All_Mixtures(tree_list,mixt_tree);
  Free(tree_list);

  /*! Set the number of rate classes to their original values
   */
  tree = mixt_tree;
  n = 0;
  do
    {
      tree->mod->ras->n_catg = orig_catg[n];
#ifdef BEAGLE
      tree->mod->optimizing_topology = false;
      //Reset the rate categories to their original weights
      if(orig_catg[n] > 2){
          update_beagle_ras(tree->mod);
      }
#endif
      if(tree->mod->ras->invar == YES) tree->mod->ras->n_catg--;
      tree = tree->next_mixt;
      n++;
    }
  while(tree);

  Free(orig_catg);

  /*! Only the first two trees for each mixture have been modified so
    ! far -> need to update the other trees by copying the modified trees
    ! onto them.
   */
  tree = mixt_tree;
  do
    {
      if(tree != mixt_tree) Copy_Tree(mixt_tree,tree);
      tree = tree->next;
    }
  while(tree);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Spr_Random_Explore(t_tree *tree, phydbl anneal_temp, phydbl prop_spr, int do_rnd, int max_cycles)
{
  int step,i,n_targets,n_rand;
  t_tree *best_tree;
  scalar_dbl **best_bl;
  t_node *rnd_node;
  t_edge *b_target,*b_residual,**target_list,*rnd_edge;
  phydbl true_best_lnL;

  if(tree->lock_topo == YES)
    {
      PhyML_Fprintf(stderr,"\n== The tree topology is locked.");
      PhyML_Fprintf(stderr,"\n== Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
      Exit("\n");
    }

  Set_Both_Sides(NO,tree);
  Pars(NULL,tree);
  Lk(NULL,tree);

  tree->mod->s_opt->max_depth_path    = (int)(tree->n_otu/3);
  tree->mod->s_opt->min_diff_lk_move  = 0.1;
  tree->mod->s_opt->spr_lnL           = NO;
  tree->mod->s_opt->spr_pars          = NO;
  tree->mod->s_opt->deepest_path      = 0;
  tree->best_pars                     = tree->c_pars;
  step                                = 0;
  true_best_lnL                       = tree->c_lnL;
  best_tree                           = Make_Tree_From_Scratch(tree->n_otu,tree->data);
  best_bl                             = Copy_Br_Len(tree);
  target_list                         = (t_edge **)mCalloc(2*tree->n_otu-3,sizeof(t_edge *));
  n_targets                           = 0;
  tree->annealing_temp                = anneal_temp;
  Copy_Tree(tree,best_tree);

  do
    {

      if(do_rnd == YES)
        {
          n_rand = 0;
          do
            {
              rnd_node = tree->a_nodes[Rand_Int(tree->n_otu,2*tree->n_otu-3)];
              assert(rnd_node != tree->n_root && rnd_node->tax == NO);
              
              rnd_edge = rnd_node->b[Rand_Int(0,2)];
              
              Prune_Subtree(rnd_node,
                            rnd_node == rnd_edge->left ? rnd_edge->rght : rnd_edge->left,
                            &b_target,
                            &b_residual,
                            tree);
              
              n_targets = 0;
              for(i=0;i<3;i++)
                if(b_target->left->v[i] != b_target->rght)
                  Get_List_Of_Adjacent_Targets(b_target->left,b_target->left->v[i],NULL,&target_list,&n_targets,0,5);
              
              for(i=0;i<3;i++)
                if(b_target->rght->v[i] != b_target->left)
                  Get_List_Of_Adjacent_Targets(b_target->rght,b_target->rght->v[i],NULL,&target_list,&n_targets,0,5);
              
              if(n_targets > 0) b_target = target_list[Rand_Int(0,n_targets-1)];
              
              assert(b_target != NULL);
              
              Graft_Subtree(b_target,rnd_node,NULL,b_residual,NULL,tree);
              
              n_rand++;
            }
          while(n_rand != 1);
        }
          
      Set_Both_Sides(YES,tree);
      Lk(NULL,tree);
      Pars(NULL,tree);
      
      Print_Lk_And_Pars(tree);
            
      if(tree->annealing_temp < 0.0) tree->annealing_temp = 0.0;
      if(prop_spr > 1.0)             prop_spr = 1.0;
      
      tree->best_lnL       = tree->c_lnL;
      tree->best_pars      = tree->c_pars;
      Spr(UNLIKELY,prop_spr,tree);
      
      tree->annealing_temp -= 0.5;
      prop_spr+=0.2;

      Optimiz_All_Free_Param(tree,(tree->io->quiet == YES)?(0):(tree->verbose > VL0));
      Optimize_Br_Len_Serie(2,tree);

      if(tree->io->print_trace)
        {
          char *s = Write_Tree(tree);
          PhyML_Fprintf(tree->io->fp_out_trace,"[%f]%s\n",tree->c_lnL,s); fflush(tree->io->fp_out_trace);
          if((tree->io->print_site_lnl) && (!tree->mod->s_opt->spr_pars))
            {
              Print_Site_Lk(tree,tree->io->fp_out_lk);
              fflush(tree->io->fp_out_lk);
            }
          Free(s);
        }

      if(tree->io->print_json_trace == YES) JSON_Tree_Io(tree,tree->io->fp_out_json_trace); 
      
      /* Record the current best log-likelihood and parsimony */
      if(tree->c_lnL > true_best_lnL)
        {
          true_best_lnL = tree->c_lnL;
          For(i,2*tree->n_otu-1) Free_Scalar_Dbl(best_bl[i]);
          Free(best_bl);
          best_bl = Copy_Br_Len(tree);
          Copy_Tree(tree,best_tree); /* Record tree topology, branch lengths and model parameters */
        }
      Transfer_Br_Len_To_Tree(best_bl,tree);
      Copy_Tree(best_tree,tree);

    }
  while(++step <= max_cycles && tree->mod->s_opt->n_improvements > 0 && tree->mod->s_opt->max_spr_depth  > 1);

  Free(target_list);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

  
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

// ** EOF: spr.c
