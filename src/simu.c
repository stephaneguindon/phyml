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
  Spr_List_Of_Trees(tree);
  return;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int Simu(t_tree *tree, int n_step_max)
{
  phydbl old_loglk,n_iter,lambda;
  int i,n_neg,n_tested,n_without_swap,n_tot_swap,step,it_lim_without_swap;
  t_edge **sorted_b,**tested_b;

  sorted_b = (t_edge **)mCalloc(tree->n_otu-3,sizeof(t_edge *));
  tested_b = (t_edge **)mCalloc(tree->n_otu-3,sizeof(t_edge *));

  old_loglk           = UNLIKELY;
  tree->c_lnL         = UNLIKELY;
  n_iter              = 1.0;
  it_lim_without_swap = (tree->mod->ras->invar)?(1):(1);
  n_tested            = 0;
  n_without_swap      = 0;
  step                = 0;
  lambda              = .75;
  n_tot_swap          = 0;

  Update_Dirs(tree);

  if(tree->lock_topo)
    {
      PhyML_Printf("\n== The tree topology is locked.");
      PhyML_Printf("\n== Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
      Warn_And_Exit("");
    }

  do
    {
      ++step;

      if(n_tested || step == 1) MIXT_Set_Alias_Subpatt(YES,tree);
      old_loglk = tree->c_lnL;
      Set_Both_Sides(YES,tree);
      Lk(NULL,tree);
      MIXT_Set_Alias_Subpatt(NO,tree);

      
      if(tree->c_lnL > old_loglk - 5.0)
        {
          MIXT_Set_Alias_Subpatt(YES,tree);
          Optimize_Br_Len_Serie(tree);
          Set_Both_Sides(YES,tree);
          Lk(NULL,tree);
          MIXT_Set_Alias_Subpatt(NO,tree);
        }

      if(tree->c_lnL < old_loglk)
        {
          if(tree->verbose > VL2 && tree->io->quiet == NO) PhyML_Printf("\n\n. Moving backward\n");
          if(!Mov_Backward_Topo_Bl(tree,old_loglk,tested_b,n_tested))
            {
              PhyML_Printf("\n== tree->c_lnL: %f old_loglk: %f",tree->c_lnL,old_loglk);
              Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
            }
          if(!tree->n_swap) n_neg = 0;
          Record_Br_Len(tree);
          Set_Both_Sides(YES,tree);
          Lk(NULL,tree);
        }

      if(step > n_step_max) break;

      if(tree->io->print_trace)
        {
          char *s = Write_Tree(tree,NO);
          PhyML_Fprintf(tree->io->fp_out_trace,"[%f]%s\n",tree->c_lnL,s); fflush(tree->io->fp_out_trace);
          if(tree->io->print_site_lnl) Print_Site_Lk(tree,tree->io->fp_out_lk); fflush(tree->io->fp_out_lk);
          Free(s);
        }

      if(tree->io->print_json_trace == YES) JSON_Tree_Io(tree,tree->io->fp_out_json_trace); 

      if(tree->verbose > VL2 && tree->io->quiet == NO) Print_Lk(tree,"[Topology           ]");

      if((FABS(old_loglk-tree->c_lnL) < tree->mod->s_opt->min_diff_lk_global) || (n_without_swap > it_lim_without_swap)) break;

      Fix_All(tree);
      n_neg = 0;
      For(i,2*tree->n_otu-3)
        if((!tree->a_edges[i]->left->tax) &&
           (!tree->a_edges[i]->rght->tax))
          {
            NNI(tree,tree->a_edges[i],NO);
          }
      Select_Edges_To_Swap(tree,sorted_b,&n_neg);
      Sort_Edges_NNI_Score(tree,sorted_b,n_neg);
      Optimiz_Ext_Br(tree);
      Update_Bl(tree,lambda);

      n_tested = 0;
      For(i,(int)CEIL((phydbl)n_neg*(lambda)))
        tested_b[n_tested++] = sorted_b[i];

      Make_N_Swap(tree,tested_b,0,n_tested);

      if((tree->verbose > VL2) && (tree->io->quiet == NO)) PhyML_Printf("[# nnis=%3d]",n_tested);

      n_tot_swap += n_tested;

      if(n_tested > 0) n_without_swap = 0;
      else             n_without_swap++;

      n_iter+=1.0;
    }
  while(1);

  Free(sorted_b);
  Free(tested_b);

  return n_tot_swap;
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
  phydbl best_score;

  *n_neg = 0;

  For(i,2*tree->n_otu-3)
    {
      b = tree->a_edges[i];
      best_score = b->nni->score;
      
      if((!b->left->tax) && (!b->rght->tax) && (b->nni->score < -tree->mod->s_opt->min_diff_lk_move))
        {
          Check_NNI_Scores_Around(b->left,b->rght,b,&best_score,tree);
          Check_NNI_Scores_Around(b->rght,b->left,b,&best_score,tree);
          if(best_score < b->nni->score) continue;
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

      if(b[i]->nni->best_conf == 1)
        {
          n1 = b[i]->left->v[b[i]->l_v2];
          n2 = b[i]->left;
          n3 = b[i]->rght;
          n4 = b[i]->rght->v[b[i]->r_v1];
        }
      else if(b[i]->nni->best_conf == 2)
        {
          n1 = b[i]->left->v[b[i]->l_v2];
          n2 = b[i]->left;
          n3 = b[i]->rght;
          n4 = b[i]->rght->v[b[i]->r_v2];
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
          
          if(b[i]->nni->best_conf == 1)
            {
              n1 = b[i]->left->v[b[i]->l_v2];
              n2 = b[i]->left;
              n3 = b[i]->rght;
              n4 = b[i]->rght->v[b[i]->r_v1];
            }
          else if(b[i]->nni->best_conf == 2)
            {
              n1 = b[i]->left->v[b[i]->l_v2];
              n2 = b[i]->left;
              n3 = b[i]->rght;
              n4 = b[i]->rght->v[b[i]->r_v2];
            }

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
          
          if(b[i]->nni->best_conf == 1)
            {
              n1 = b[i]->left->v[b[i]->l_v2];
              n2 = b[i]->left;
              n3 = b[i]->rght;
              n4 = b[i]->rght->v[b[i]->r_v1];
            }
          else if(b[i]->nni->best_conf == 2)
            {
              n1 = b[i]->left->v[b[i]->l_v2];
              n2 = b[i]->left;
              n3 = b[i]->rght;
              n4 = b[i]->rght->v[b[i]->r_v2];
            }

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
          
          if(b[i]->nni->best_conf == 1)
            {
              n1 = b[i]->left->v[b[i]->l_v2];
              n2 = b[i]->left;
              n3 = b[i]->rght;
              n4 = b[i]->rght->v[b[i]->r_v1];
            }
          else if(b[i]->nni->best_conf == 2)
            {
              n1 = b[i]->left->v[b[i]->l_v2];
              n2 = b[i]->left;
              n3 = b[i]->rght;
              n4 = b[i]->rght->v[b[i]->r_v2];
            }
          
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
          
          if(b[i]->nni->best_conf == 1)
            {
              n1 = b[i]->left->v[b[i]->l_v2];
              n2 = b[i]->left;
              n3 = b[i]->rght;
              n4 = b[i]->rght->v[b[i]->r_v1];
            }
          else if(b[i]->nni->best_conf == 2)
            {
              n1 = b[i]->left->v[b[i]->l_v2];
              n2 = b[i]->left;
              n3 = b[i]->rght;
              n4 = b[i]->rght->v[b[i]->r_v2];
            }
          
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
  For(i,3)
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

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

