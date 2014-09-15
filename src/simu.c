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

void Simu_Loop(t_tree *mixt_tree)
{
  phydbl lk_old;

  SPR_Shuffle(mixt_tree);

  Set_Both_Sides(YES,mixt_tree);
  Lk(NULL,mixt_tree);
  mixt_tree->best_lnL = mixt_tree->c_lnL;
  PhyML_Printf("\n. Current value of log-likelihood: %f",mixt_tree->c_lnL);

  int n_tot_moves = 0;
  do
    {
      lk_old = mixt_tree->c_lnL;
      Optimiz_All_Free_Param(mixt_tree,(mixt_tree->io->quiet)?(0):(mixt_tree->mod->s_opt->print));
      n_tot_moves = Simu(mixt_tree,10);
      if(!n_tot_moves) break;
    }
  while(mixt_tree->c_lnL > lk_old + mixt_tree->mod->s_opt->min_diff_lk_local);

  do
    {
      Round_Optimize(mixt_tree,mixt_tree->data,ROUND_MAX);
      if(!Check_NNI_Five_Branches(mixt_tree)) break;
    }
  while(1);

  if((mixt_tree->mod->s_opt->print) &&
     (!mixt_tree->io->quiet)) PhyML_Printf("\n");
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
  /* it_lim_without_swap = (tree->mod->ras->invar)?(5):(2); */
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

      if(tree->c_lnL < old_loglk)
        {
          if((tree->mod->s_opt->print) && (!tree->io->quiet)) PhyML_Printf("\n\n. Moving backward\n");
          if(!Mov_Backward_Topo_Bl(tree,old_loglk,tested_b,n_tested))
            Exit("\n== Err. mov_back failed\n");
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

      if((tree->mod->s_opt->print) && (!tree->io->quiet)) Print_Lk(tree,"[Topology           ]");

/*       if(((tree->c_lnL > old_loglk) && (FABS(old_loglk-tree->c_lnL) < tree->mod->s_opt->min_diff_lk_local)) || (n_without_swap > it_lim_without_swap)) break; */
      if((FABS(old_loglk-tree->c_lnL) < tree->mod->s_opt->min_diff_lk_global) || (n_without_swap > it_lim_without_swap)) break;

      Fill_Dir_Table(tree);
      Fix_All(tree);
      n_neg = 0;
      For(i,2*tree->n_otu-3)
        if((!tree->a_edges[i]->left->tax) &&
           (!tree->a_edges[i]->rght->tax))
          NNI(tree,tree->a_edges[i],0);

      Select_Edges_To_Swap(tree,sorted_b,&n_neg);
      Sort_Edges_NNI_Score(tree,sorted_b,n_neg);
      Optimiz_Ext_Br(tree);
      Update_Bl(tree,lambda);

      n_tested = 0;
      For(i,(int)CEIL((phydbl)n_neg*(lambda)))
        tested_b[n_tested++] = sorted_b[i];

      Make_N_Swap(tree,tested_b,0,n_tested);

      if((tree->mod->s_opt->print) && (!tree->io->quiet)) PhyML_Printf("[# nnis=%3d]",n_tested);

      n_tot_swap += n_tested;

      if(n_tested > 0) n_without_swap = 0;
      else             n_without_swap++;

      n_iter+=1.0;
    }
  while(1);

/*   Round_Optimize(tree,tree->data); */

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
  lambda              = .75;
  n_tot_swap          = 0;

  Update_Dirs(tree);

  if((tree->mod->s_opt->print) && (!tree->io->quiet)) PhyML_Printf("\n\n. Starting simultaneous NNI moves (parsimony criterion)...\n");

  do
    {
      ++step;
      
      if(step > n_step_max) break;
      
      each--;
      
      Set_Both_Sides(YES,tree);
      Pars(NULL,tree);
      
      if((tree->mod->s_opt->print) && (!tree->io->quiet))
        {
          Print_Pars(tree);
          if(step > 1) (n_tested > 1)?(printf("[%4d NNIs]",n_tested)):(printf("[%4d NNI ]",n_tested));
        }
      
      
      if(FABS(old_pars - tree->c_pars) < SMALL) break;
      
      if((tree->c_pars > old_pars) && (step > 1))
        {
          if((tree->mod->s_opt->print) && (!tree->io->quiet))
            PhyML_Printf("\n\n. Moving backward (topoLlogy) \n");
          if(!Mov_Backward_Topo_Pars(tree,old_pars,tested_b,n_tested))
            Exit("\n. Err: mov_back failed\n");
          if(!tree->n_swap) n_neg = 0;
          
          
          Set_Both_Sides(YES,tree);
          Pars(NULL,tree);
        }
      else
        {
          
          old_pars = tree->c_pars;
          Fill_Dir_Table(tree);
          
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
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Update_Bl(t_tree *tree, phydbl fact)
{
  int i;
  t_edge *b,*orig;

  For(i,2*tree->n_otu-3)
    {
      b = tree->a_edges[i];
      /* b->l->v = b->l_old->v + (b->nni->l0 - b->l_old->v)*fact;       */

      orig = b;
      do
    {
      b->l->v = b->l_old->v + (b->nni->l0 - b->l_old->v)*fact;
      if(b->next) b = b->next;
      else         b = b->next;
    }
      while(b);
      b = orig;

    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Make_N_Swap(t_tree *tree,t_edge **b, int beg, int end)
{
  int i;
  int dim;
  t_edge *orig;

  dim = 2*tree->n_otu-2;

  /* PhyML_Printf("\n. Beg Actually performing swaps\n"); */
  tree->n_swap = 0;
  for(i=beg;i<end;i++)
    {
      /* we use t_dir here to take into account previous modifications of the topology */
      /* printf("\n. Swap on edge %d [%d %d %d %d]",b[i]->num, */
      /* 	     b[i]->nni->swap_node_v2->v[tree->t_dir[b[i]->nni->swap_node_v2->num*dim+b[i]->nni->swap_node_v1->num]]->num, */
      /* 	     b[i]->nni->swap_node_v2->num, */
      /* 	     b[i]->nni->swap_node_v3->num, */
      /* 	     b[i]->nni->swap_node_v3->v[tree->t_dir[b[i]->nni->swap_node_v3->num*dim+b[i]->nni->swap_node_v4->num]]->num); */

      Swap(b[i]->nni->swap_node_v2->v[tree->t_dir[b[i]->nni->swap_node_v2->num*dim+b[i]->nni->swap_node_v1->num]],
       b[i]->nni->swap_node_v2,
       b[i]->nni->swap_node_v3,
       b[i]->nni->swap_node_v3->v[tree->t_dir[b[i]->nni->swap_node_v3->num*dim+b[i]->nni->swap_node_v4->num]],
       tree);

      if(!Check_Topo_Constraints(tree,tree->io->cstr_tree))
    {
      /* Undo this swap as it violates one of the topological constraints
         defined in the input constraint tree
      */
      Swap(b[i]->nni->swap_node_v2->v[tree->t_dir[b[i]->nni->swap_node_v2->num*dim+b[i]->nni->swap_node_v1->num]],
           b[i]->nni->swap_node_v2,
           b[i]->nni->swap_node_v3,
           b[i]->nni->swap_node_v3->v[tree->t_dir[b[i]->nni->swap_node_v3->num*dim+b[i]->nni->swap_node_v4->num]],
           tree);
    }

      if(tree->n_root)
    {
      tree->n_root->v[2] = tree->e_root->left;
      tree->n_root->v[1] = tree->e_root->rght;
    }

      orig = b[i];
      do
    {
      b[i]->l->v = b[i]->nni->best_l;
      if(b[i]->next) b[i] = b[i]->next;
      else            b[i] = b[i]->next;
    }
      while(b[i]);
      b[i] = orig;

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
  int dim;
  t_edge *orig;

  dim = 2*tree->n_otu-2;

  sorted_b = (t_edge **)mCalloc(tree->n_otu-3,sizeof(t_edge *));

  j=0;
  For(i,2*tree->n_otu-3) if((!tree->a_edges[i]->left->tax) &&
                (!tree->a_edges[i]->rght->tax))
                              sorted_b[j++] = tree->a_edges[i];

  Sort_Edges_NNI_Score(tree,sorted_b,tree->n_otu-3);

  if(sorted_b[0]->nni->score < -0.0)
    {
      b = sorted_b[0];
      return_value = 1;

      Swap(b->nni->swap_node_v2->v[tree->t_dir[b->nni->swap_node_v2->num*dim+b->nni->swap_node_v1->num]],
       b->nni->swap_node_v2,
       b->nni->swap_node_v3,
       b->nni->swap_node_v3->v[tree->t_dir[b->nni->swap_node_v3->num*dim+b->nni->swap_node_v4->num]],
       tree);

      if(!Check_Topo_Constraints(tree,tree->io->cstr_tree))
    {
      /* Undo this swap as it violates one of the topological constraints
         defined in the input constraint tree
      */
      Swap(b->nni->swap_node_v2->v[tree->t_dir[b->nni->swap_node_v2->num*dim+b->nni->swap_node_v1->num]],
           b->nni->swap_node_v2,
           b->nni->swap_node_v3,
           b->nni->swap_node_v3->v[tree->t_dir[b->nni->swap_node_v3->num*dim+b->nni->swap_node_v4->num]],
           tree);
    }

      /* b->l->v = b->nni->best_l; */

      orig = b;
      do
    {
      b->l->v = b->nni->best_l;
      if(b->next) b = b->next;
      else         b = b->next;
    }
      while(b);
      b = orig;


/*       (b->nni->best_conf == 1)? */
/* 	(Swap(b->left->v[b->l_v2],b->left,b->rght,b->rght->v[b->r_v1],tree)): */
/* 	(Swap(b->left->v[b->l_v2],b->left,b->rght,b->rght->v[b->r_v2],tree)); */

/*       b->l->v =  */
/* 	(b->nni->best_conf == 1)? */
/* 	(b->nni->l1): */
/* 	(b->nni->l2); */


    }
  else return_value = 0;

  Free(sorted_b);

  return return_value;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


int Mov_Backward_Topo_Bl(t_tree *tree, phydbl lk_old, t_edge **tested_b, int n_tested)
{
  phydbl **l_init;
  int i,j,step,beg,end;
  t_edge *b,*orig;

  l_init = (phydbl **)mCalloc(2*tree->n_otu-3,sizeof(phydbl *));

  For(i,2*tree->n_otu-3) l_init[i] = MIXT_Get_Lengths_Of_This_Edge(tree->a_edges[i],tree);

  step = 2;
  do
    {
      For(i,2*tree->n_otu-3)
    {
      b = tree->a_edges[i];

      /* b->l->v = b->l_old->v + (1./step) * (l_init[i] - b->l_old->v); */

      j = 0;
      orig = b;
      do
        {
          b->l->v = b->l_old->v + (1./step) * (l_init[i][j] - b->l_old->v);
          if(b->next) b = b->next;
          else         b = b->next;
          j++;
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

      For(i,2*tree->n_otu-3)
    {
      b = tree->a_edges[i];

      orig = b;
      do
        {
          b->l->v = b->l_old->v;
          if(b->next) b = b->next;
          else         b = b->next;
        }
      while(b);
      b = orig;
    }


      Set_Both_Sides(NO,tree);
      Lk(NULL,tree);
    }

  For(i,2*tree->n_otu-3) Free(l_init[i]);
  Free(l_init);

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

    }while((tree->c_pars > pars_old) && (step < 1000));


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
  int dim;
  t_edge *orig;

  dim = 2*tree->n_otu-2;

  if(end>beg)
    {
      for(i=beg;i<end;i++)
    {

/* 	  PhyML_Printf("MOV BACK UNSWAP Edge %d Swap nodes %d(%d) %d %d %d(%d)\n", */
/* 		 b[i]->num, */
/* 		 b[i]->nni->swap_node_v2->v[tree->t_dir[b[i]->nni->swap_node_v2->num][b[i]->nni->swap_node_v1->num]]->num, */
/* 		 b[i]->nni->swap_node_v4->num, */
/* 		 b[i]->nni->swap_node_v2->num, */
/* 		 b[i]->nni->swap_node_v3->num, */
/* 		 b[i]->nni->swap_node_v3->v[tree->t_dir[b[i]->nni->swap_node_v3->num][b[i]->nni->swap_node_v4->num]]->num, */
/* 		 b[i]->nni->swap_node_v1->num */
/* 		 ); */

      Swap(b[i]->nni->swap_node_v2->v[tree->t_dir[b[i]->nni->swap_node_v2->num*dim+b[i]->nni->swap_node_v1->num]],
           b[i]->nni->swap_node_v2,
           b[i]->nni->swap_node_v3,
           b[i]->nni->swap_node_v3->v[tree->t_dir[b[i]->nni->swap_node_v3->num*dim+b[i]->nni->swap_node_v4->num]],
           tree);

      if(!Check_Topo_Constraints(tree,tree->io->cstr_tree))
        {
          /* Undo this swap as it violates one of the topological constraints
         defined in the input constraint tree
          */
          Swap(b[i]->nni->swap_node_v2->v[tree->t_dir[b[i]->nni->swap_node_v2->num*dim+b[i]->nni->swap_node_v1->num]],
           b[i]->nni->swap_node_v2,
           b[i]->nni->swap_node_v3,
           b[i]->nni->swap_node_v3->v[tree->t_dir[b[i]->nni->swap_node_v3->num*dim+b[i]->nni->swap_node_v4->num]],
           tree);
        }


/* 	  (b[i]->nni->best_conf == 1)? */
/* 	    (Swap(b[i]->left->v[b[i]->l_v2],b[i]->left,b[i]->rght,b[i]->rght->v[b[i]->r_v1],tree)): */
/* 	    (Swap(b[i]->left->v[b[i]->l_v2],b[i]->left,b[i]->rght,b[i]->rght->v[b[i]->r_v2],tree)); */

      /* b[i]->l->v = b[i]->l_old->v; */


      orig = b[i];
      do
        {
          b[i]->l->v = b[i]->l_old->v;
          if(b[i]->next) b[i] = b[i]->next;
          else            b[i] = b[i]->next;
        }
      while(b[i]);
      b[i] = orig;

    }
    }
  else
    {
      for(i=beg-1;i>=end;i--)
    {
      Swap(b[i]->nni->swap_node_v2->v[tree->t_dir[b[i]->nni->swap_node_v2->num*dim+b[i]->nni->swap_node_v1->num]],
           b[i]->nni->swap_node_v2,
           b[i]->nni->swap_node_v3,
           b[i]->nni->swap_node_v3->v[tree->t_dir[b[i]->nni->swap_node_v3->num*dim+b[i]->nni->swap_node_v4->num]],
           tree);

      if(!Check_Topo_Constraints(tree,tree->io->cstr_tree))
        {
          /* Undo this swap as it violates one of the topological constraints
         defined in the input constraint tree
          */
          Swap(b[i]->nni->swap_node_v2->v[tree->t_dir[b[i]->nni->swap_node_v2->num*dim+b[i]->nni->swap_node_v1->num]],
           b[i]->nni->swap_node_v2,
           b[i]->nni->swap_node_v3,
           b[i]->nni->swap_node_v3->v[tree->t_dir[b[i]->nni->swap_node_v3->num*dim+b[i]->nni->swap_node_v4->num]],
           tree);
        }


      /* b[i]->l->v = b[i]->l_old->v; */

      orig = b[i];
      do
        {
          b[i]->l->v = b[i]->l_old->v;
          if(b[i]->next) b[i] = b[i]->next;
          else            b[i] = b[i]->next;
        }
      while(b[i]);
      b[i] = orig;

    }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Swap_N_Branch(t_tree *tree,t_edge **b, int beg, int end)
{
  int i;
  int dim;
  t_edge *orig;

  dim = 2*tree->n_otu-2;

  if(end>beg)
    {
      for(i=beg;i<end;i++)
    {
      Swap(b[i]->nni->swap_node_v2->v[tree->t_dir[b[i]->nni->swap_node_v2->num*dim+b[i]->nni->swap_node_v1->num]],
           b[i]->nni->swap_node_v2,
           b[i]->nni->swap_node_v3,
           b[i]->nni->swap_node_v3->v[tree->t_dir[b[i]->nni->swap_node_v3->num*dim+b[i]->nni->swap_node_v4->num]],
           tree);

      if(!Check_Topo_Constraints(tree,tree->io->cstr_tree))
        {
          /* Undo this swap as it violates one of the topological constraints
         defined in the input constraint tree
          */
          Swap(b[i]->nni->swap_node_v2->v[tree->t_dir[b[i]->nni->swap_node_v2->num*dim+b[i]->nni->swap_node_v1->num]],
           b[i]->nni->swap_node_v2,
           b[i]->nni->swap_node_v3,
           b[i]->nni->swap_node_v3->v[tree->t_dir[b[i]->nni->swap_node_v3->num*dim+b[i]->nni->swap_node_v4->num]],
           tree);
        }

      /* b[i]->l->v = b[i]->nni->best_l; */

      orig = b[i];
      do
        {
          b[i]->l->v = b[i]->nni->best_l;
          if(b[i]->next) b[i] = b[i]->next;
          else            b[i] = b[i]->next;
        }
      while(b[i]);
      b[i] = orig;

    }
    }
  else
    {
      for(i=beg-1;i>=end;i--)
    {
      Swap(b[i]->nni->swap_node_v2->v[tree->t_dir[b[i]->nni->swap_node_v2->num*dim+b[i]->nni->swap_node_v1->num]],
           b[i]->nni->swap_node_v2,
           b[i]->nni->swap_node_v3,
           b[i]->nni->swap_node_v3->v[tree->t_dir[b[i]->nni->swap_node_v3->num*dim+b[i]->nni->swap_node_v4->num]],
           tree);

      if(!Check_Topo_Constraints(tree,tree->io->cstr_tree))
        {
          /* Undo this swap as it violates one of the topological constraints
         defined in the input constraint tree
          */
          Swap(b[i]->nni->swap_node_v2->v[tree->t_dir[b[i]->nni->swap_node_v2->num*dim+b[i]->nni->swap_node_v1->num]],
           b[i]->nni->swap_node_v2,
           b[i]->nni->swap_node_v3,
           b[i]->nni->swap_node_v3->v[tree->t_dir[b[i]->nni->swap_node_v3->num*dim+b[i]->nni->swap_node_v4->num]],
           tree);
        }

      /* b[i]->l->v = b[i]->nni->best_l; */

      orig = b[i];
      do
        {
          b[i]->l->v = b[i]->nni->best_l;
          if(b[i]->next) b[i] = b[i]->next;
          else            b[i] = b[i]->next;
        }
      while(b[i]);
      b[i] = orig;

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

