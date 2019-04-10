/*

PHYML :  a program that  computes maximum likelihood  phylogenies from
DNA or AA homologous sequences

Copyright (C) Stephane Guindon. Oct 2003 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/

#include "simu.h"
#ifdef BEAGLE
#include "beagle_utils.h"
#endif


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Simu_Loop(t_tree *tree)
{
  return;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int Simu(t_tree *tree, int n_step_max, phydbl delta_lnL, phydbl init_T, phydbl delta_T, int min_n_edges_traversed)
{
  phydbl old_loglk,delta;
  unsigned int n_round;
  time_t t_cur;
  
  tree->c_lnL = UNLIKELY;
  delta = UNLIKELY;
  old_loglk = UNLIKELY;
  n_round = 0;
  tree->annealing_temp = init_T;
  do
    {
      for(int i=0;i<2*tree->n_otu-3;++i) tree->a_edges[i]->l->v *= Rgamma((phydbl)(0.2*n_round+1),(phydbl)(1./(0.2*n_round+1)));
      old_loglk = tree->c_lnL;
      Set_Both_Sides(NO,tree);
      tree->tip_root = Rand_Int(0,tree->n_otu-1);
      Lk(NULL,tree);
      tree->n_edges_traversed = 0;
      tree->fully_nni_opt = YES;
      NNI_Traversal(tree->a_nodes[tree->tip_root],
                    tree->a_nodes[tree->tip_root]->v[0],
                    NULL,
                    tree->a_nodes[tree->tip_root]->b[0],
                    YES,
                    tree);
      delta = tree->c_lnL - old_loglk;
      tree->annealing_temp -= delta_T;
      if(tree->annealing_temp < 0.0) tree->annealing_temp = 0.0;
      n_round++;
      time(&t_cur);
 
      PhyML_Printf("\n. %5ds lnL: %12G T: %12G %4d/%4d",
                   (int)(t_cur-tree->t_beg),
                   tree->c_lnL,
                   tree->annealing_temp,
                   tree->n_edges_traversed,
                   tree->n_otu);

      if((n_round >= n_step_max || tree->fully_nni_opt == YES) && Are_Equal(tree->annealing_temp,0.0,1.E-3) && delta < delta_lnL) break;
    }
  while(1);

  return 1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Simu_Pars(t_tree *tree, int n_step_max)
{
  phydbl old_pars,n_iter,lambda;
  int i,n_neg,n_tested,n_without_swap,n_tot_swap,step;
  t_edge **sorted_b,**tested_b;
  int each;

  sorted_b = (t_edge **)mCalloc(tree->n_otu-3,sizeof(t_edge *));
  tested_b = (t_edge **)mCalloc(tree->n_otu-3,sizeof(t_edge *));

  old_pars            = 0;
  tree->c_pars        = 0;
  n_iter              = 1.0;
  n_tested            = 0;
  n_without_swap      = 0;
  step                = 0;
  each                = 4;
  lambda              = 0.5;
  n_tot_swap          = 0;

  Update_Dirs(tree);

  if((tree->verbose > VL0) && (tree->io->quiet == NO)) PhyML_Printf("\n\n. Starting simultaneous NNI moves (parsimony criterion)...\n");

  do
    {
      ++step;
      
      if(step > n_step_max) break;
      
      each--;
      
      Set_Both_Sides(YES,tree);
      Pars(NULL,tree);
      
      if((tree->verbose > VL0) && (tree->io->quiet == NO))
        {
          Print_Pars(tree);
          if(step > 1) (n_tested > 1)?(printf("[%4d NNIs]",n_tested)):(printf("[%4d NNI ]",n_tested));
        }
      
      
      if(FABS(old_pars - tree->c_pars) < SMALL) break;
      
      if((tree->c_pars > old_pars) && (step > 1))
        {
          if(tree->verbose > VL0 && tree->io->quiet == NO)
            PhyML_Printf("\n\n. Moving backward (topology) \n");
          if(!Mov_Backward_Topo_Pars(tree,old_pars,tested_b,n_tested))
            Generic_Exit(__FILE__,__LINE__,__FUNCTION__);    
          if(!tree->n_swap) n_neg = 0;
          
          Set_Both_Sides(YES,tree);
          Pars(NULL,tree);
        }
      else
        {          
          old_pars = tree->c_pars;
          
          n_neg = 0;
          For(i,2*tree->n_otu-3)
            if((!tree->a_edges[i]->left->tax) &&
               (!tree->a_edges[i]->rght->tax))
              NNI_Pars(tree,tree->a_edges[i],NO);
          
          Select_Edges_To_Swap(tree,sorted_b,&n_neg);
          Sort_Edges_NNI_Score(tree,sorted_b,n_neg);
          
          n_tested = 0;
          For(i,(int)CEIL((phydbl)n_neg*(lambda)))
            tested_b[n_tested++] = sorted_b[i];
          
          Make_N_Swap(tree,tested_b,0,n_tested);
          
          n_tot_swap += n_tested;
          
          if(n_tested > 0) n_without_swap = 0;
          else             n_without_swap++;
        }
      n_iter+=1.0;
    }
  while(1);
  
  Free(sorted_b);
  Free(tested_b);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Select_Edges_To_Swap(t_tree *tree, t_edge **sorted_b, int *n_neg)
{
  int i;
  t_edge *b;
  /* phydbl best_score; */

  *n_neg = 0;

  For(i,2*tree->n_otu-3)
    {
      b = tree->a_edges[i];
      /* best_score = b->nni->score; */
      
      if((!b->left->tax) && (!b->rght->tax) && (b->nni->score < -tree->mod->s_opt->min_diff_lk_move))
        {
          /* // Evaluate NNIs on edges at distance 1 */
          /* Check_NNI_Scores_Around(b->left,b->rght,b,&best_score,tree); */
          /* Check_NNI_Scores_Around(b->rght,b->left,b,&best_score,tree);           */
          
          /* // Evaluate NNIs on edges at distance 2 */
          /* Check_NNI_Scores_Around(b->left,b->left->v[b->l_v1],b,&best_score,tree); */
          /* Check_NNI_Scores_Around(b->left,b->left->v[b->l_v2],b,&best_score,tree); */

          /* // Evaluate NNIs on edges at distance 2 */
          /* Check_NNI_Scores_Around(b->rght,b->rght->v[b->r_v1],b,&best_score,tree); */
          /* Check_NNI_Scores_Around(b->rght,b->rght->v[b->r_v2],b,&best_score,tree); */

          /* if(best_score < b->nni->score) continue; */

          sorted_b[*n_neg] = b;
          (*n_neg)++;
        }


      /* if((!b->left->tax) && (!b->rght->tax) && (b->nni->score < -tree->mod->s_opt->min_diff_lk_move)) */
      /*   { */
      /*     sorted_b[*n_neg] = b; */
      /*     (*n_neg)++; */
      /*   } */

    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Update_Bl(t_tree *tree, phydbl fact)
{
  int i;
  scalar_dbl *l,*l_old,*l0;

  For(i,2*tree->n_otu-3)
    {
      l     = tree->a_edges[i]->l;
      l_old = tree->a_edges[i]->l_old;
      l0    = tree->a_edges[i]->nni->l0;

      do
        {
          l->v  = l_old->v + (l0->v - l_old->v)*fact;
          l     = l->next;
          l_old = l_old->next;
          l0    = l0->next;
        }
      while(l);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Make_N_Swap(t_tree *tree,t_edge **b, int beg, int end)
{
  int i;
  /* t_edge *orig; */
  t_node *n1,*n2,*n3,*n4;

  n1 = n2 = n3 = n4 = NULL;

  tree->n_swap = 0;
  for(i=beg;i<end;i++)
    {
      n1 = n2 = n3 = n4 = NULL;

      n1 = b[i]->nni->swap_node_v1;
      n2 = b[i]->nni->swap_node_v2;
      n3 = b[i]->nni->swap_node_v3;
      n4 = b[i]->nni->swap_node_v4;

      /* if(b[i]->nni->best_conf == 1) */
      /*   { */
      /*     n1 = b[i]->left->v[b[i]->l_v2]; */
      /*     n2 = b[i]->left; */
      /*     n3 = b[i]->rght; */
      /*     n4 = b[i]->rght->v[b[i]->r_v1]; */
      /*   } */
      /* else if(b[i]->nni->best_conf == 2) */
      /*   { */
      /*     n1 = b[i]->left->v[b[i]->l_v2]; */
      /*     n2 = b[i]->left; */
      /*     n3 = b[i]->rght; */
      /*     n4 = b[i]->rght->v[b[i]->r_v2]; */
      /*   } */

      if(b[i]->nni->best_conf == 1)
        {
          if(n1 != b[i]->left->v[b[i]->l_v2] ||
             /* n2 != b[i]->left || */
             /* n3 != b[i]->rght || */
             n4 != b[i]->rght->v[b[i]->r_v1]) continue;
        }
      else if(b[i]->nni->best_conf == 2)
        {
          if(n1 != b[i]->left->v[b[i]->l_v2] ||
             /* n2 != b[i]->left || */
             /* n3 != b[i]->rght || */
             n4 != b[i]->rght->v[b[i]->r_v2]) continue;
        }
            
      Swap(n1,n2,n3,n4,tree);
      
      if(!Check_Topo_Constraints(tree,tree->io->cstr_tree))
        {
          /* Undo this swap as it violates one of the topological constraints
             defined in the input constraint tree
          */
          Swap(n4,n2,n3,n1,tree);
        }
      
      if(tree->n_root)
        {
          tree->n_root->v[2] = tree->e_root->left;
          tree->n_root->v[1] = tree->e_root->rght;
        }
      
      Copy_Scalar_Dbl(b[i]->nni->best_l,b[i]->l);
      Copy_Scalar_Dbl(b[i]->nni->best_v,b[i]->l_var);

      /* orig = b[i]; */
      /* do */
      /*   { */
      /*     b[i]->l->v = b[i]->nni->best_l; */
      /*     b[i] = b[i]->next; */
      /*   } */
      /* while(b[i]); */
      /* b[i] = orig; */
      
      tree->n_swap++;
    }


  /* PhyML_Printf("\n. End Actually performing swaps\n"); */

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int Make_Best_Swap(t_tree *tree)
{
  int i,j,return_value;
  t_edge *b,**sorted_b;
  /* t_edge *orig; */
  t_node *n1,*n2,*n3,*n4;


  sorted_b = (t_edge **)mCalloc(tree->n_otu-3,sizeof(t_edge *));
  
  j=0;
  For(i,2*tree->n_otu-3) 
    if((!tree->a_edges[i]->left->tax) &&
       (!tree->a_edges[i]->rght->tax))
      sorted_b[j++] = tree->a_edges[i];
  
  Sort_Edges_NNI_Score(tree,sorted_b,tree->n_otu-3);
  
  if(sorted_b[0]->nni->score < -0.0)
    {
      b = sorted_b[0];
      return_value = 1;
      
      n1 = n2 = n3 = n4 = NULL;

      /* n1 = b->nni->swap_node_v1; */
      /* n2 = b->nni->swap_node_v2; */
      /* n3 = b->nni->swap_node_v3; */
      /* n4 = b->nni->swap_node_v4; */

      if(b->nni->best_conf == 1)
        {
          n1 = b->left->v[b->l_v2];
          n2 = b->left;
          n3 = b->rght;
          n4 = b->rght->v[b->r_v1];
        }
      else if(b->nni->best_conf == 2)
        {
          n1 = b->left->v[b->l_v2];
          n2 = b->left;
          n3 = b->rght;
          n4 = b->rght->v[b->r_v2];
        }

      Swap(n1,n2,n3,n4,tree);
      
      if(!Check_Topo_Constraints(tree,tree->io->cstr_tree))
        {
          /* Undo this swap as it violates one of the topological constraints
             defined in the input constraint tree
          */
          Swap(n4,n2,n3,n1,tree);
        }
      
      /* b->l->v = b->nni->best_l; */
      
      Copy_Scalar_Dbl(b->nni->best_l,b->l);
      Copy_Scalar_Dbl(b->nni->best_v,b->l_var);

      /* orig = b; */
      /* do */
      /*   { */
      /*     b->l->v = b->nni->best_l; */
      /*     b = b->next; */
      /*   } */
      /* while(b); */
      /* b = orig; */
    }
  else return_value = 0;
  
  Free(sorted_b);

  return return_value;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int Mov_Backward_Topo_Bl(t_tree *tree, phydbl lk_old, t_edge **tested_b, int n_tested)
{
  scalar_dbl **l_init,**v_init;
  int i,step,beg,end;
  t_edge *b,*orig;

  l_init = (scalar_dbl **)mCalloc(2*tree->n_otu-3,sizeof(scalar_dbl *));
  v_init = (scalar_dbl **)mCalloc(2*tree->n_otu-3,sizeof(scalar_dbl *));

  For(i,2*tree->n_otu-3) 
    {
      l_init[i] = Duplicate_Scalar_Dbl(tree->a_edges[i]->l);
      v_init[i] = Duplicate_Scalar_Dbl(tree->a_edges[i]->l_var);
    }

  step = 2;
  do
    {
      For(i,2*tree->n_otu-3)
        {
          b = tree->a_edges[i];

          orig = b;
          do
            {
              b->l->v = b->l_old->v + (1./step) * (l_init[i]->v - b->l_old->v);
              if(b->next == NULL) break;
              b = b->next;
              l_init[i] = l_init[i]->next;
            }
          while(b);
          b = orig;
        }
      
      beg = (int)FLOOR((phydbl)n_tested/(step-1));
      end = 0;
      Unswap_N_Branch(tree,tested_b,beg,end);
      beg = 0;
      end = (int)FLOOR((phydbl)n_tested/step);
      Swap_N_Branch(tree,tested_b,beg,end);
      
      if(!end) tree->n_swap = 0;
      
      Set_Both_Sides(NO,tree);
      Lk(NULL,tree);

      step++;
      
    }while((tree->c_lnL < lk_old) && (step < 1000));
  
  
  if(step == 1000)
    {
      if(tree->n_swap)  Exit("\n== Err. in Mov_Backward_Topo_Bl (n_swap > 0)\n");
      
      Restore_Br_Len(tree);
      
      Set_Both_Sides(NO,tree);
      Lk(NULL,tree);
    }

  For(i,2*tree->n_otu-3) 
    {
      while(l_init[i]->prev) l_init[i] = l_init[i]->prev; 
      Free_Scalar_Dbl(l_init[i]);
    }
  Free(l_init);

  For(i,2*tree->n_otu-3) 
    {
      while(v_init[i]->prev) v_init[i] = v_init[i]->prev;
      Free_Scalar_Dbl(v_init[i]);
    }
  Free(v_init);

  tree->n_swap = 0;
  For(i,2*tree->n_otu-3)
    {
      if(tree->a_edges[i]->nni->score < 0.0) tree->n_swap++;
      tree->a_edges[i]->nni->score = +1.0;
    }


  if(tree->c_lnL > lk_old)                                return  1;
  else if((tree->c_lnL > lk_old-tree->mod->s_opt->min_diff_lk_local) &&
          (tree->c_lnL < lk_old+tree->mod->s_opt->min_diff_lk_local)) return -1;
  else                                                    return  0;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int Mov_Backward_Topo_Pars(t_tree *tree, int pars_old, t_edge **tested_b, int n_tested)
{
  int i,step,beg,end;

  step = 2;
  do
    {
      beg = (int)FLOOR((phydbl)n_tested/(step-1));
      end = 0;
      Unswap_N_Branch(tree,tested_b,beg,end);
      beg = 0;
      end = (int)FLOOR((phydbl)n_tested/step);
      Swap_N_Branch(tree,tested_b,beg,end);

      if(!end) tree->n_swap = 0;

      Set_Both_Sides(NO,tree);
      Pars(NULL,tree);

      step++;

    }
  while((tree->c_pars > pars_old) && (step < 1000));


  if(step == 1000)
    {
      if(tree->n_swap)  Exit("\n. Err. in Mov_Backward_Topo_Bl (n_swap > 0)\n");

      Set_Both_Sides(NO,tree);
      Pars(NULL,tree);
    }

  tree->n_swap = 0;
  For(i,2*tree->n_otu-3)
    {
      if(tree->a_edges[i]->nni->score < 0.0) tree->n_swap++;
      tree->a_edges[i]->nni->score = +1.0;
    }


  if(tree->c_pars < pars_old)       return  1;
  else if(tree->c_pars == pars_old) return -1;
  else                              return  0;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Unswap_N_Branch(t_tree *tree, t_edge **b, int beg, int end)
{
  int i;
  /* t_edge *orig; */
  t_node *n1,*n2,*n3,*n4;

  n1 = n2 = n3 = n4 = NULL;
          
  if(end>beg)
    {
      for(i=beg;i<end;i++)
        {
          n1 = n2 = n3 = n4 = NULL;
          
          n1 = b[i]->nni->swap_node_v4;
          n2 = b[i]->nni->swap_node_v2;
          n3 = b[i]->nni->swap_node_v3;
          n4 = b[i]->nni->swap_node_v1;

          /* if(b[i]->nni->best_conf == 1) */
          /*   { */
          /*     n1 = b[i]->left->v[b[i]->l_v2]; */
          /*     n2 = b[i]->left; */
          /*     n3 = b[i]->rght; */
          /*     n4 = b[i]->rght->v[b[i]->r_v1]; */
          /*   } */
          /* else if(b[i]->nni->best_conf == 2) */
          /*   { */
          /*     n1 = b[i]->left->v[b[i]->l_v2]; */
          /*     n2 = b[i]->left; */
          /*     n3 = b[i]->rght; */
          /*     n4 = b[i]->rght->v[b[i]->r_v2]; */
          /*   } */

          Swap(n1,n2,n3,n4,tree);
          
          if(!Check_Topo_Constraints(tree,tree->io->cstr_tree))
            {
              /* Undo this swap as it violates one of the topological constraints
                 defined in the input constraint tree
              */
              Swap(n4,n2,n3,n1,tree);
            }
          
          
          /* 	  (b[i]->nni->best_conf == 1)? */
          /* 	    (Swap(b[i]->left->v[b[i]->l_v2],b[i]->left,b[i]->rght,b[i]->rght->v[b[i]->r_v1],tree)): */
          /* 	    (Swap(b[i]->left->v[b[i]->l_v2],b[i]->left,b[i]->rght,b[i]->rght->v[b[i]->r_v2],tree)); */
          
          /* b[i]->l->v = b[i]->l_old->v; */
          
          Copy_Scalar_Dbl(b[i]->l_old,b[i]->l);
          Copy_Scalar_Dbl(b[i]->l_var_old,b[i]->l_var);
          
          /* orig = b[i]; */
          /* do */
          /*   { */
          /*     b[i]->l->v = b[i]->l_old->v; */
          /*     b[i] = b[i]->next; */
          /*   } */
          /* while(b[i]); */
          /* b[i] = orig; */

        }
    }
  else
    {
      for(i=beg-1;i>=end;i--)
        {
          n1 = n2 = n3 = n4 = NULL;
          
          n1 = b[i]->nni->swap_node_v4;
          n2 = b[i]->nni->swap_node_v2;
          n3 = b[i]->nni->swap_node_v3;
          n4 = b[i]->nni->swap_node_v1;

          /* if(b[i]->nni->best_conf == 1) */
          /*   { */
          /*     n1 = b[i]->left->v[b[i]->l_v2]; */
          /*     n2 = b[i]->left; */
          /*     n3 = b[i]->rght; */
          /*     n4 = b[i]->rght->v[b[i]->r_v1]; */
          /*   } */
          /* else if(b[i]->nni->best_conf == 2) */
          /*   { */
          /*     n1 = b[i]->left->v[b[i]->l_v2]; */
          /*     n2 = b[i]->left; */
          /*     n3 = b[i]->rght; */
          /*     n4 = b[i]->rght->v[b[i]->r_v2]; */
          /*   } */

          Swap(n1,n2,n3,n4,tree);
          
          if(!Check_Topo_Constraints(tree,tree->io->cstr_tree))
            {
              /* Undo this swap as it violates one of the topological constraints
                 defined in the input constraint tree
              */
              Swap(n4,n2,n3,n1,tree);
            }
          
          
          /* b[i]->l->v = b[i]->l_old->v; */
          Copy_Scalar_Dbl(b[i]->l_old,b[i]->l);
          Copy_Scalar_Dbl(b[i]->l_var_old,b[i]->l_var);
          

          /* orig = b[i]; */
          /* do */
          /*   { */
          /*     b[i]->l->v = b[i]->l_old->v; */
          /*     b[i] = b[i]->next; */
          /*   } */
          /* while(b[i]); */
          /* b[i] = orig; */
          
        }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Swap_N_Branch(t_tree *tree,t_edge **b, int beg, int end)
{
  int i;
  /* t_edge *orig; */
  t_node *n1,*n2,*n3,*n4;

  n1 = n2 = n3 = n4 = NULL;

  if(end>beg)
    {
      for(i=beg;i<end;i++)
        {
          n1 = n2 = n3 = n4 = NULL;
          
          n1 = b[i]->nni->swap_node_v1;
          n2 = b[i]->nni->swap_node_v2;
          n3 = b[i]->nni->swap_node_v3;
          n4 = b[i]->nni->swap_node_v4;


          /* if(b[i]->nni->best_conf == 1) */
          /*   { */
          /*     n1 = b[i]->left->v[b[i]->l_v2]; */
          /*     n2 = b[i]->left; */
          /*     n3 = b[i]->rght; */
          /*     n4 = b[i]->rght->v[b[i]->r_v1]; */
          /*   } */
          /* else if(b[i]->nni->best_conf == 2) */
          /*   { */
          /*     n1 = b[i]->left->v[b[i]->l_v2]; */
          /*     n2 = b[i]->left; */
          /*     n3 = b[i]->rght; */
          /*     n4 = b[i]->rght->v[b[i]->r_v2]; */
          /*   } */
          
          Swap(n1,n2,n3,n4,tree);
          
          if(!Check_Topo_Constraints(tree,tree->io->cstr_tree))
            {
              /* Undo this swap as it violates one of the topological constraints
                 defined in the input constraint tree
              */
              Swap(n4,n2,n3,n1,tree);
            }
          
          /* b[i]->l->v = b[i]->nni->best_l; */

          Copy_Scalar_Dbl(b[i]->nni->best_l,b[i]->l);
          Copy_Scalar_Dbl(b[i]->nni->best_v,b[i]->l_var);

          /* orig = b[i]; */
          /* do */
          /*   { */
          /*     b[i]->l->v = b[i]->nni->best_l; */
          /*     b[i] = b[i]->next; */
          /*   } */
          /* while(b[i]); */
          /* b[i] = orig; */
        }
    }
  else
    {
      for(i=beg-1;i>=end;i--)
        {

          n1 = n2 = n3 = n4 = NULL;

          n1 = b[i]->nni->swap_node_v1;
          n2 = b[i]->nni->swap_node_v2;
          n3 = b[i]->nni->swap_node_v3;
          n4 = b[i]->nni->swap_node_v4;


          /* if(b[i]->nni->best_conf == 1) */
          /*   { */
          /*     n1 = b[i]->left->v[b[i]->l_v2]; */
          /*     n2 = b[i]->left; */
          /*     n3 = b[i]->rght; */
          /*     n4 = b[i]->rght->v[b[i]->r_v1]; */
          /*   } */
          /* else if(b[i]->nni->best_conf == 2) */
          /*   { */
          /*     n1 = b[i]->left->v[b[i]->l_v2]; */
          /*     n2 = b[i]->left; */
          /*     n3 = b[i]->rght; */
          /*     n4 = b[i]->rght->v[b[i]->r_v2]; */
          /*   } */
          
          Swap(n1,n2,n3,n4,tree);
          
          if(!Check_Topo_Constraints(tree,tree->io->cstr_tree))
            {
              /* Undo this swap as it violates one of the topological constraints
                 defined in the input constraint tree
              */
              Swap(n4,n2,n3,n1,tree);
            }
          
          /* b[i]->l->v = b[i]->nni->best_l; */
          
          Copy_Scalar_Dbl(b[i]->nni->best_l,b[i]->l);
          Copy_Scalar_Dbl(b[i]->nni->best_v,b[i]->l_var);

          /* orig = b[i]; */
          /* do */
          /*   { */
          /*     b[i]->l->v = b[i]->nni->best_l; */
          /*     b[i] = b[i]->next; */
          /*   } */
          /* while(b[i]); */
          /* b[i] = orig; */
        }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Check_NNI_Scores_Around(t_node *a, t_node *d, t_edge *b, phydbl *best_score, t_tree *tree)
{
  int i;

  if(d->tax) return;
  
  for(i=0;i<3;i++)
    {
      if((d->v[i] != a) && (!d->v[i]->tax))
        {
          if((d->b[i]->nni->score > *best_score-1.E-10) &&
             (d->b[i]->nni->score < *best_score+1.E-10)) /* ties */
            {
              d->b[i]->nni->score = *best_score+1.;
            }

          if(d->b[i]->nni->score < *best_score)
            {
              *best_score = d->b[i]->nni->score;
            }
        }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
/*
        v
        |
        |
        a
       / \
      d   u
     / \
    v1 v2   
*/
void NNI_Traversal(t_node *a, t_node *d, t_node *v, t_edge *b, int opt_edges, t_tree *tree)
{
  int i,dir1, dir2;

  /* printf("\n. a: %d d: %d b->is_alive ? %d -- [%d - %d]", a->num,d->num,b->is_alive,a->tax,d->tax); */
  
  if(d->tax == YES)
    {
      if(opt_edges == YES) Br_Len_Opt(&(b->l->v),b,tree);
      return;
    }
  else if(a->tax == YES)
    {
      if(opt_edges == YES && a->tax == YES)  Br_Len_Opt(&(b->l->v),b,tree);
      for(i=0;i<3;++i)
        if(d->v[i] != a)
          {
            Update_Partial_Lk(tree,d->b[i],d);
            NNI_Traversal(d,d->v[i],a,d->b[i],opt_edges,tree);
          }
      Update_Partial_Lk(tree,b,d);
    }
  else
    {
      tree->n_edges_traversed++;
      NNI_Core(a,d,v,b,opt_edges,tree);

      dir1 = dir2 = -1;
      for(i=0;i<3;++i)
        if(d->v[i] != a)
          {
            if(dir1 < 0) dir1 = i;
            else dir2 = i;
          }

      Update_Partial_Lk(tree,d->b[dir1],d);
      NNI_Traversal(d,d->v[dir1],a,d->b[dir1],opt_edges,tree);
      Update_Partial_Lk(tree,d->b[dir2],d);
      NNI_Traversal(d,d->v[dir2],a,d->b[dir2],opt_edges,tree);

      Update_Partial_Lk(tree,b,d);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void NNI_Core(t_node *a, t_node *d, t_node *v, t_edge *b, int opt_edges, t_tree *tree)
{
  phydbl lk0,lk1,lk2;
  phydbl rnd;
  t_node *v1,*v2,*u,*dum;
  scalar_dbl *l0,*l1,*l2;
  phydbl p_accept;
  int i;
  
  
  l0 = l1 = l2 = NULL;
  
  lk0 = UNLIKELY;
  lk1 = UNLIKELY;
  lk2 = UNLIKELY;
  
  v1 = v2 = NULL;
  for(i=0;i<3;++i)
    if(d->v[i] != a)
      {
        if(v1 == NULL) v1 = d->v[i];
        else           v2 = d->v[i];
      }
  assert(v1 != NULL);
  assert(v2 != NULL);
  
  dum = NULL;
  rnd = Uni();
  if(rnd < .5)
    {
      dum = v1;
      v1  = v2;
      v2  = dum;
    }
  
  u = NULL;
  for(i=0;i<3;++i) if(a->v[i] != d && a->v[i] != v) { u = a->v[i]; break; }
  
  if(opt_edges == YES) Br_Len_Opt(&(b->l->v),b,tree);
  lk0 = Lk(b,tree);
  l0 = Duplicate_Scalar_Dbl(b->l);
  
  /* Swap_Partial_Lk_Extra(b,a,0,tree); */
  /* Swap_Partial_Lk_Extra(b,d,1,tree); */
  
  // First NNI
  Swap(v1,d,a,u,tree);
  // Update partial likelihood looking up
  Update_Partial_Lk(tree,b,a);
  // Update partial likelihood looking down
  Update_Partial_Lk(tree,b,d);
  // Evaluate likelihood
  if(opt_edges == YES) Br_Len_Opt(&(b->l->v),b,tree);
  lk1 = Lk(b,tree);
  /* if(lk1 > lk0 + tree->mod->s_opt->min_diff_lk_move) */
  /*   { */
  /*     tree->fully_nni_opt = NO; */
  /*     return; */
  /*   } */
  l1 = Duplicate_Scalar_Dbl(b->l);
  // Unswap
  Swap(u,d,a,v1,tree);              
  

  
  // Second NNI
  Swap(v2,d,a,u,tree);
  // Update partial likelihood looking up
  Update_Partial_Lk(tree,b,a);
  // Update partial likelihood looking down
  Update_Partial_Lk(tree,b,d);
  // Evaluate likelihood
  if(opt_edges == YES) Br_Len_Opt(&(b->l->v),b,tree);
  lk2 = Lk(b,tree);
  /* if(lk2 > lk0 + tree->mod->s_opt->min_diff_lk_move) */
  /*   { */
  /*     tree->fully_nni_opt = NO; */
  /*     Free_Scalar_Dbl(l1); */
  /*     return; */
  /*   } */
  l2 = Duplicate_Scalar_Dbl(b->l);
  // Unswap
  Swap(u,d,a,v2,tree);

  /* Swap_Partial_Lk_Extra(b,a,0,tree); */
  /* Swap_Partial_Lk_Extra(b,d,1,tree); */
    
  /* if((u->tax == YES && !strcmp(u->name,"tax57")) || */
  /*    (v1->tax == YES && !strcmp(v1->name,"tax57")) || */
  /*    (v2->tax == YES && !strcmp(v2->name,"tax57"))) */
  /*   printf("\n. lk0: %G lk1: %G lk2: %G l0: %G l1: %G l2: %G",lk0,lk1,lk2,l0->v,l1->v,l2->v); */

  /* printf("\n. a: %d d: %d -- lk0: %f lk1: %f lk2: %f p: %G %G %G %G",a->num,d->num,lk0,lk1,lk2,p_accept,l0->v,l1->v,l2->v); */

  p_accept = exp((lk1-lk0)/(tree->annealing_temp+1.E-6));
  if(Are_Equal(lk1,lk0,tree->mod->s_opt->min_diff_lk_local) && Are_Equal(tree->annealing_temp,0.0,1.E-3)) p_accept = .0;
  rnd = Uni();
  
  
  if(rnd < p_accept && lk2 < lk1)
    {
      Swap(v1,d,a,u,tree);
      Copy_Scalar_Dbl(l1,b->l);
      tree->c_lnL = lk1;
      Update_Partial_Lk(tree,b,a);
      Update_Partial_Lk(tree,b,d);
      tree->fully_nni_opt = NO;
    }
  else
    {
      p_accept = exp((lk2-lk0)/(tree->annealing_temp+1.E-6));
      if(Are_Equal(lk2,lk0,tree->mod->s_opt->min_diff_lk_local) && Are_Equal(tree->annealing_temp,0.0,1.E-3)) p_accept = .0;
      rnd = Uni();
      if(rnd < p_accept)
        {
          Swap(v2,d,a,u,tree);
          Copy_Scalar_Dbl(l2,b->l);
          tree->c_lnL = lk2;
          Update_Partial_Lk(tree,b,a);
          Update_Partial_Lk(tree,b,d);
          tree->fully_nni_opt = NO;
        }
      else
        {
          Update_Partial_Lk(tree,b,a);
          Update_Partial_Lk(tree,b,d);
          Copy_Scalar_Dbl(l0,b->l);
          tree->c_lnL = lk0; 
        }
    }
  
  Update_PMat_At_Given_Edge(b,tree);
  
  if(l0) Free_Scalar_Dbl(l0);
  if(l1) Free_Scalar_Dbl(l1);
  if(l2) Free_Scalar_Dbl(l2);
}
