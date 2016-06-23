/*

PhyML:  a program that  computes maximum likelihood phyLOGenies from
DNA or AA homoLOGous sequences.

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

/*
** BIG: Some big number.
*/

/* #define BIG  1e05 */

/*
** Global vars.
**
**   - cur_lk:        The current likelihood of the tree.
**   - subtree_dist:  The average subtree distances matrix.
**   - seq_dist:      The sequence distance matrix.
**   - optim_cand:    Array for holding candidate moves for local and global branch
**                    length optimization.
**   - rgrft_cand:    Array for holding candidate regraft positions.
**   - v_tmp:         The central t_node of the temporary regraft structure for
**                    estimating changes in likelihood.
**   - path:          The path through the tree during the recursive tree length
**                    calculation.
**   - sum_scale_tmp  Array for temporarily storing scaling factors.
**   - p_lk_tmp:      Temporary partial likelihood storage.
**   - e_brent:       A temporary t_edge to use for estimating distances using Brent.

**   - tree->mod->s_opt->wim_n_rgrft:       Number of promising regraft positions to consider when
                      performing all improving SPR moves.
**   - tree->mod->s_opt->wim_n_optim:       Number of candidate moves on which to perform local branch
                      length optimization.
**   - tree->mod->s_opt->wim_max_dist:      Maximum regraft distance to consider.
**   - tree->mod->s_opt->wim_n_globl:       Number of candidates moves on which to perform global branch
                      length optimization.
**   - tree->mod->s_opt->wim_n_best:        Number of promising regraft positions to consider when
                      performing only the best SPR move.

**   - nr_d_l:        Total number of change in tree length calculations done.
**   - nr_d_lk:       Total number of change in likelihood calculations done.
**   - nr_loc:        Total number of local branch length optimizations done.
**   - nr_glb:        Total number of global branch length optimizations done.
*/

phydbl   cur_lk, **subtree_dist, *sum_scale_tmp, *p_lk_tmp;
matrix  *seq_dist;
_move_ **optim_cand, **rgrft_cand;
t_node    *v_tmp=NULL, **path;
t_edge    *e_brent=NULL;
int      nr_d_L, nr_d_lk, nr_loc, nr_glb;

/*
** Init_SPR: Initialize the SPR algorithm: allocate memory and set variables.
**
** Parameters:
**   - tree: The current tree to use for initialization.
*/

void Init_SPR (t_tree *tree)
{
  int   i, nr_nodes, nr_edges;
  t_node *u_0, *u_1, *u_2;

  /*
  ** Get the SPR parameter values.
  */
  nr_edges = 2*tree->n_otu-3;

  if(tree->mod->s_opt->wim_n_rgrft  < 0) tree->mod->s_opt->wim_n_rgrft = 1 + nr_edges / 5;
  if(tree->mod->s_opt->wim_n_globl  < 0) tree->mod->s_opt->wim_n_globl = 1 + nr_edges / 10;
  if(tree->mod->s_opt->wim_max_dist < 0) tree->mod->s_opt->wim_max_dist = 1 + nr_edges / 10;
  if(tree->mod->s_opt->wim_n_optim  < 0) tree->mod->s_opt->wim_n_optim = 100;
  if(tree->mod->s_opt->wim_n_best   < 0) tree->mod->s_opt->wim_n_best = tree->mod->s_opt->wim_n_rgrft; /* can't
                                                    *  be
                                                    *  anything else
                                                        */


  /*
  ** If it doesn't exist yet, create the temporary regraft structure:
  ** a central t_node with three edges and tip nodes adjacent to it.
  */
  if (v_tmp == NULL)
  {

    v_tmp=Make_Node_Light(0);
    v_tmp->tax = 0;
    u_0=Make_Node_Light(1);
    u_0->tax = 1;
    u_1=Make_Node_Light(2);
    u_1->tax = 1;
    u_2=Make_Node_Light(3);
    u_2->tax = 1;

    v_tmp->v[0] = u_0;
    v_tmp->v[1] = u_1;
    v_tmp->v[2] = u_2;
    u_0->v[0]   = v_tmp;
    u_1->v[0]   = v_tmp;
    u_2->v[0]   = v_tmp;

    t_edge *edge_0 = Make_Edge_Light (v_tmp, u_0, 0);
    Make_Edge_Lk (edge_0, tree);
    t_edge *edge_1 = Make_Edge_Light (v_tmp, u_1,1);
    Make_Edge_Lk (edge_1, tree);
    t_edge *edge_2 = Make_Edge_Light (v_tmp, u_2,2);
    Make_Edge_Lk (edge_2, tree);


/*     For(i,tree->data->crunch_len) */
/*       { */
/* 	For(j,tree->mod->ras->n_catg) */
/* 	  { */
/* 	    Free(edge_0->p_lk_rght[i][j]); */
/* 	  } */
/* 	Free(edge_0->p_lk_rght[i]); */
/*       } */
    Free(edge_0->p_lk_rght);

    if(!edge_0->rght->tax) Free(edge_0->sum_scale_rght);


/*     For(i,tree->data->crunch_len) */
/*       { */
/* 	For(j,tree->mod->ras->n_catg) */
/* 	  { */
/* 	    Free(edge_1->p_lk_rght[i][j]); */
/* 	  } */
/* 	Free(edge_1->p_lk_rght[i]); */
/*       } */
    Free(edge_1->p_lk_rght);

    if(!edge_1->rght->tax) Free(edge_1->sum_scale_rght);


/*     For(i,tree->data->crunch_len) */
/*       { */
/* 	For(j,tree->mod->ras->n_catg) */
/* 	  { */
/* 	    Free(edge_2->p_lk_rght[i][j]); */
/* 	  } */
/* 	Free(edge_2->p_lk_rght[i]); */
/*       } */
    Free(edge_2->p_lk_rght);

    if(!edge_2->rght->tax) Free(edge_2->sum_scale_rght);
  }

  /*
  ** If it doesn't exist yet, create the temporary edge.
  */
  if (e_brent == NULL)
  {

    u_1=Make_Node_Light(4);
    u_1->tax = 1;
    u_2=Make_Node_Light(5);
    u_2->tax = 1;
    u_1->v[0] = u_2;
    u_2->v[0] = u_1;
    t_edge *edge_4 = Make_Edge_Light (u_1, u_2, 3);
    Make_Edge_Lk (edge_4, tree);
    e_brent = u_1->b[0];


/*     For(i,tree->data->crunch_len) */
/*       { */
/* 	For(j,tree->mod->ras->n_catg) */
/* 	  { */
/* 	    Free(edge_4->p_lk_rght[i][j]); */
/* 	  } */
/* 	Free(edge_4->p_lk_rght[i]); */
/*       } */
    Free(edge_4->p_lk_rght);

    if(!edge_4->rght->tax) Free(edge_4->sum_scale_rght);


/*     For(i,tree->data->crunch_len) */
/*       { */
/* 	For(j,tree->mod->ras->n_catg) */
/* 	  { */
/* 	    Free(edge_4->p_lk_left[i][j]); */
/* 	  } */
/* 	Free(edge_4->p_lk_left[i]); */
/*       } */
    Free(edge_4->p_lk_left);

    if(!edge_4->left->tax) Free(edge_4->sum_scale_left);
  }

  /*
  ** Allocate memory for temporarily storing partial likelihoods and
  ** scaling factors.
  */
  p_lk_tmp = (phydbl *)mCalloc (tree->n_pattern*tree->mod->ras->n_catg*tree->mod->ns, sizeof (phydbl));
/*   p_lk_tmp = (phydbl ***)mCalloc (tree->n_pattern, sizeof (phydbl **)); */
/*   for (i = 0; i < tree->n_pattern; i++) */
/*   { */
/*     p_lk_tmp[i] = (phydbl **)mCalloc (tree->mod->ras->n_catg, sizeof (phydbl *)); */
/*     for (j = 0; j < tree->mod->ras->n_catg; j++) */
/*     { */
/*       p_lk_tmp[i][j] = (phydbl *)mCalloc (tree->mod->ns, sizeof (phydbl)); */
/*     } */
/*   } */
  sum_scale_tmp = (phydbl *)mCalloc (tree->n_pattern, sizeof (phydbl));

  /*
  ** Allocate memory for storing the average subtree distances.
  */
  nr_nodes = 2*tree->n_otu-2;
  subtree_dist = (phydbl **)malloc (nr_nodes * sizeof (phydbl *));
  for (i = 0; i < nr_nodes; i++)
    {
      subtree_dist[i] = (phydbl *)malloc (nr_nodes * sizeof (phydbl));
    }

  /*
  ** Allocate memory for storing the candidate regraft positions and
  ** t_edge length optimization moves.
  */
  rgrft_cand = (_move_ **)malloc (MAX(tree->mod->s_opt->wim_n_rgrft,tree->mod->s_opt->wim_n_best) * sizeof (_move_ *));
  for (i = 0; i < MAX(tree->mod->s_opt->wim_n_rgrft,tree->mod->s_opt->wim_n_best); i++)
    {
      rgrft_cand[i] = (_move_ *)malloc (sizeof (_move_));
      rgrft_cand[i]->path = (t_node **)malloc ((tree->mod->s_opt->wim_max_dist+2) * sizeof (t_node *));
    }
  optim_cand = (_move_ **)malloc (tree->mod->s_opt->wim_n_optim * sizeof (_move_ *));
  for (i = 0; i < tree->mod->s_opt->wim_n_optim; i++)
    {
      optim_cand[i] = (_move_ *)malloc (sizeof (_move_));
      optim_cand[i]->path = (t_node **)malloc ((tree->mod->s_opt->wim_max_dist+2) * sizeof (t_node *));
    }
  path = (t_node **)malloc ((tree->mod->s_opt->wim_max_dist+2) * sizeof (t_node *));

  if(!tree->mat)
    {
      seq_dist = ML_Dist (tree->data, tree->mod);
      tree->mat = seq_dist;
    }
  else
    seq_dist = tree->mat;

  /*
  ** Set variables.
  */
  nr_d_L = 0;
  nr_d_lk = 0;
  nr_loc = 0;
  nr_glb = 0;
}


/*
** Clean_SPR: Free up the used memory.
**
** Parameters:
**   - tree: The current tree.
*/

void Clean_SPR (t_tree *tree)
{
  int i;

  /*
  ** Clean up the temporary regraft structure.
  */
  Free_Node (v_tmp->v[0]);
  Free_Node (v_tmp->v[1]);
  Free_Node (v_tmp->v[2]);
  v_tmp->b[0]->p_lk_rght = NULL;
  Free_Edge_Lk (v_tmp->b[0]);
  Free_Edge (v_tmp->b[0]);
  v_tmp->b[1]->p_lk_rght = NULL;
  Free_Edge_Lk(v_tmp->b[1]);
  Free_Edge(v_tmp->b[1]);
  v_tmp->b[2]->p_lk_rght = NULL;
  Free_Edge_Lk (v_tmp->b[2]);
  Free_Edge (v_tmp->b[2]);
  Free_Node (v_tmp);
  v_tmp = NULL;

  /*
  ** Clean up the temporary edge.
  */
  Free_Node (e_brent->left);
  Free_Node (e_brent->rght);
  e_brent->p_lk_left = NULL;
  e_brent->p_lk_rght = NULL;
  Free_Edge_Lk (e_brent);
  Free_Edge (e_brent);
  e_brent = NULL;

  /*
  ** Free the temporary partial likelihood and scaling memory.
  */
/*   for (i = 0; i < tree->n_pattern; i++) */
/*   { */
/*     for (j = 0; j < tree->mod->ras->n_catg; j++) */
/*     { */
/*       free (p_lk_tmp[i][j]); */
/*     } */
/*     free (p_lk_tmp[i]); */
/*   } */
  free (p_lk_tmp);
  free (sum_scale_tmp);

  /*
  ** Free the subtree distance matrix.
  */
  for (i = 0; i < 2*tree->n_otu - 2; i++)
  {
    free (subtree_dist[i]);
  }
  free (subtree_dist);

  /*
  ** Free the arrays for storing the candidate regrafting positions and
  ** t_edge length optimization moves.
  */
  for (i = 0; i < MAX(tree->mod->s_opt->wim_n_rgrft,tree->mod->s_opt->wim_n_best); i++)
  {
    free (rgrft_cand[i]->path);
    free (rgrft_cand[i]);
  }
  free (rgrft_cand);
  for (i = 0; i < tree->mod->s_opt->wim_n_optim; i++)
  {
    free (optim_cand[i]->path);
    free (optim_cand[i]);
  }
  free (optim_cand);
  free (path);


  /*
  ** Print some statistics (for "research" purposes only).
  */
/*   PhyML_Printf ("nr_d_L:  %d\n", nr_d_L); */
/*   PhyML_Printf ("nr_d_lk: %d\n", nr_d_lk); */
/*   PhyML_Printf ("nr_loc:  %d\n", nr_loc); */
/*   PhyML_Printf ("nr_glb:  %d\n", nr_glb); */
}


/*
** Optim_SPR: Optimize the tree using SPR moves.
**
** Parameters:
**   - tree:     The tree to optimize.
**   - max_size: The maximum size (= number of taxa) of the subtrees to be
**               pruned. If m=0 or m>ntax, all possible prunings will be
**               considered.
**   - method:   The optimization method to use ("ALL" or "BEST").
*/

void Optim_SPR (t_tree *tree, int max_size, int method)
{
  int   nr_moves, improvement;

  if(tree->mod->s_opt->print) PhyML_Printf("\n\n. Starting SPR moves...\n");

  /*
  ** Calculate the current likelihood value.
  */
  Set_Both_Sides(YES,tree);
  cur_lk = Lk(NULL,tree);
  time(&(tree->t_current));
  if(tree->mod->s_opt->print) Print_Lk(tree,"topology");

  /*
  ** Optimize all t_edge lengths and calculate the new likelihood value.
  */
/*   PhyML_Printf("\n. Optimizing t_edge lengths."); */
  Optimize_Br_Len_Serie (tree);
  Set_Both_Sides(YES,tree);
  cur_lk = Lk(NULL,tree);
  time(&(tree->t_current));
  if(tree->mod->s_opt->print) Print_Lk(tree,"topology");

  /*
  ** While improvements were found, perform another round of SPR moves.
  */
  nr_moves = 0;
  improvement = 1;
  while (improvement)
    {
      /*
      ** Perform one round of SPR moves.
      */
      if (method == ALL)
    {
      improvement = Perform_SPR_Moves (tree, max_size);
    }
      else if (method == BEST)
    {
      improvement = Perform_Best_SPR (tree, max_size);
    }
      else if (method == ONE)
    {
      improvement = Perform_One_SPR (tree, max_size);
    }
      else
    {
      PhyML_Printf ("\n. Unknown SPR optimization method, bailing out...\n");
      exit (1);
    }

      /*     If an improvement was found, update statistics. */
      if(improvement)
    {
    nr_moves++;
    if((nr_moves == 1) || (nr_moves % 4 == 1))
      {
        /*
        ** Optimize model parameters.
        */
        Optimiz_All_Free_Param (tree,(tree->io->quiet)?(0):(tree->mod->s_opt->print));
        Set_Both_Sides(YES,tree);
        Lk(NULL,tree);
      }
      }

      /* Beg SG 28 May 2007 */
      if(method == BEST || method == ONE) break;
      /* Beg SG 28 May 2007 */
  }

  if(tree->mod->s_opt->print) PhyML_Printf ("\n\n. Number of SPR moves: %d\n", nr_moves);

  /*
  ** Perform a last round of optimization steps (for t_edge lengths).
  */
  Round_Optimize(tree,tree->data,ROUND_MAX);
  Check_NNI_Five_Branches(tree);
}


/*
** Perform_SPR_Moves: Perform a round of SPR moves on the tree. Prune each subtree in
**                    turn and calculate the change in tree length for each candidate
**                    regraft position. Estimate change in likelihood for the most
**                    promising moves, and perform all moves that result in an
**                    improvement. If no improvements were found at all, try local edge
**                    length optimization. If still no improvement, try global edge
**                    length optimization.
**
** Parameters:
**   - tree:     The tree to perform the SPR moves on.
**   - max_size: The maximum size (= number of taxa) of the subtrees to be
**               pruned. If m=0 or m>ntax, all possible prunings will be
**               considered.
**
** Returns:
**   If the current tree could be improved: 1.
**   Otherwise:                             0.
*/

int Perform_SPR_Moves (t_tree *tree, int max_size)
{
  int   nr_edges, i, j, candidate, improvement;
  t_node *root, *v_prune;
  t_edge *e_prune;

  /*
  ** Calculate the average subtree distances.
  */
  root = tree->a_nodes[0];
  PostOrder_v (tree, root->v[2], root->b[2]);

  /*
  ** Initialize the array of optimization candidates.
  */
  for (i = 0; i < tree->mod->s_opt->wim_n_optim; i++)
  {
    optim_cand[i]->delta_lk = -1.0*BIG;
    optim_cand[i]->d_L = -1.0*BIG;
  }

  /*
  ** Try all possible SPR moves and perform the ones that give an improvement.
  */
  nr_edges = 2*tree->n_otu - 3;
  cur_lk = tree->c_lnL;
  improvement = 0;


/*   PhyML_Printf("\n >>>>>>>>>>>>>>>>>>"); */
/*   PhyML_Printf("\n. cur_lk = %f %f",cur_lk,Lk(NULL,tree)); */

/*   PhyML_Printf ("\n. Trying SPR moves"); */
/*   PhyML_Printf ("\n.  - calculating tree distances and estimating likelihoods"); */

  for(i = 0; i < nr_edges; i++)
    {
      /*
      ** Get the next prune edge.
      */
      e_prune = tree->a_edges[i];
      /*
      ** Try right subtree if appropriate.
      */
      if (!e_prune->left->tax)
    {
      /*
      ** Clear the regraft candidate list.
      */
      for (j = 0; j < tree->mod->s_opt->wim_n_rgrft; j++)
        {
          rgrft_cand[j]->d_L = -1.0*BIG;
        }
      v_prune = e_prune->left;
/* 	  if ((max_size == 0) || (e_prune->num_tax_rght <= max_size)) */
        {
          /*
          ** Calculate changes in tree length, and estimate changes in likelihood for
          ** the most promising candidates. Perform moves that give an improvement.
          */
          Calc_Tree_Length (e_prune, v_prune, tree);
          if ((candidate = Est_Lk_Change (e_prune, v_prune, tree)) >= 0)
        {
          improvement = 1;
          Make_Move (rgrft_cand[candidate],0,tree);
/* 		  PhyML_Printf("\n. Make simple move"); */
/* 		  PhyML_Printf("\n. lk after simple move = %f",Lk(tree)); */
        }
        }
    }
      /*
      ** Try left subtree if appropriate.
      */
      if (!e_prune->rght->tax)
    {
      /*
      ** Clear the regraft candidate list.
      */
      for (j = 0; j < tree->mod->s_opt->wim_n_rgrft; j++)
        {
          rgrft_cand[j]->d_L = -1.0*BIG;
        }
      v_prune = e_prune->rght;
/* 	  if ((max_size == 0) || (e_prune->num_tax_left <= max_size)) */
        {
          /*
          ** Calculate changes in tree length, and estimate changes in likelihood for
          ** the most promising candidates. Perform moves that give an improvement.
          */
          Calc_Tree_Length (e_prune, v_prune, tree);
          if ((candidate = Est_Lk_Change (e_prune, v_prune, tree)) >= 0)
        {
          improvement = 1;
          Make_Move (rgrft_cand[candidate],0,tree);
/* 		  PhyML_Printf("\n. Make simple move"); */
/* 		  PhyML_Printf("\n. lk after simple move = %f",Lk(tree)); */
        }
        }
    }
    }

  /*
  ** If there was no improvement at all, try local t_edge length optimization at the
  ** regraft position.
  */

/*   PhyML_Printf("\n. before local = %f %f",tree->c_lnL,Lk(tree)); */

  if (!improvement)
    {
      /*       PhyML_Printf ("\n.  - performing local t_edge length optimizations"); */
      if ((candidate = Find_Optim_Local (tree)) >= 0)
    {
/*  	  PhyML_Printf("\n. make local move"); */
      improvement = 1;
      Make_Move (optim_cand[candidate],1,tree);
/* 	  PhyML_Printf("\n. lk after local move = %f",Lk(tree)); */
    }
    }


  /*
  ** If there was still no improvement, try global t_edge length optimization.
  */
/*   PhyML_Printf("\n. before global = %f %f",tree->c_lnL,Lk(tree)); */

  if (!improvement)
    {
/*       PhyML_Printf ("\n.  - performing global t_edge length optimization"); */
      if ((candidate = Find_Optim_Globl (tree)) >= 0)
    {
/*  	  PhyML_Printf("\n. make global move"); */
      improvement = 1;
      Make_Move (optim_cand[candidate],2,tree);
/* 	  PhyML_Printf("\n. lk after global move = %f",Lk(tree)); */
    }
    }

/*   PhyML_Printf("\n. after all = %f %f",tree->c_lnL,Lk(tree)); */

  /*
  ** Optimize all t_edge lengths again to make sure we got an updated
  ** likelihood value.
  */
  Set_Both_Sides(YES,tree);
  cur_lk = Lk(NULL,tree);
  root = tree->a_nodes[0];
  Optimize_Br_Len_Serie (tree);
  Set_Both_Sides(YES,tree);
  cur_lk = Lk(NULL,tree);
  time(&(tree->t_current));
  if(tree->mod->s_opt->print) Print_Lk(tree,"topoLOGy");

  /*
  ** Return the result.
  */
  return (improvement);
}


/*
** Perform_Best_SPR: Perform the best SPR move on the tree. Prune each subtree in
**                   turn and calculate the change in tree length for each candidate
**                   regraft position. Estimate change in likelihood for the most
**                   promising regraft positions, and store the best one. Then choose
**                   the best candidate over all moves. If no improving move can be
**                   found, try local t_edge length optimization, and if necessary
**                   global t_edge length optimization.
**
** Parameters:
**   - tree:     The tree to perform the SPR moves on.
**   - max_size: The maximum size (= number of taxa) of the subtrees to be
**               pruned. If m=0 or m>ntax, all possible prunings will be
**               considered.
**
** Returns:
**   If an improving move could be performed: 1.
**   Otherwise:                               0.
*/

int Perform_Best_SPR (t_tree *tree, int max_size)
{
  int   nr_edges, i, j, candidate, improvement;
  t_node *root, *v_prune;
  t_edge *e_prune;

  /*
  ** Calculate the average subtree distances.
  */
  root = tree->a_nodes[0];
  PostOrder_v (tree, root->v[2], root->b[2]);

  /*
  ** Initialize the array of optimization candidates.
  */
  for (i = 0; i < tree->mod->s_opt->wim_n_optim; i++)
  {
    optim_cand[i]->delta_lk = -1.0*BIG;
    optim_cand[i]->d_L = -1.0*BIG;
  }

  /*
  ** Try all possible SPR moves and perform the best one.
  */
  nr_edges = 2*tree->n_otu - 3;
  cur_lk = tree->c_lnL;
  improvement = 0;
/*   PhyML_Printf ("\n. Trying SPR moves"); */
/*   PhyML_Printf ("\n.  -calculating tree distances and estimating likelihoods"); */
  for (i = 0; i < nr_edges; i++)
  {
    /*
    ** Get the next prune edge.
    */
    e_prune = tree->a_edges[i];
    /*
    ** Try right subtree if appropriate.
    */
    if (!e_prune->left->tax)
    {
      /*
      ** Clear the regraft candidate list.
      */
      for (j = 0; j < tree->mod->s_opt->wim_n_best; j++)
      {
    rgrft_cand[j]->d_L = -1.0*BIG;
      }
      v_prune = e_prune->left;
/*       if ((max_size == 0) || (e_prune->num_tax_rght <= max_size)) */
      {
    /*
    ** Calculate changes in tree length, and estimate changes in likelihood for
    ** the most promising candidates. Store the best one in the optimization list.
    */
    Calc_Tree_Length (e_prune, v_prune, tree);
    candidate = Best_Lk_Change (e_prune, v_prune, tree);
      }
    }
    /*
    ** Try left subtree if appropriate.
    */
    if (!e_prune->rght->tax)
    {
      /*
      ** Clear the regraft candidate list.
      */
      for (j = 0; j < tree->mod->s_opt->wim_n_rgrft; j++)
      {
    rgrft_cand[j]->d_L = -1.0*BIG;
      }
      v_prune = e_prune->rght;
/*       if ((max_size == 0) || (e_prune->num_tax_left <= max_size)) */
      {
    /*
    ** Calculate changes in tree length, and estimate changes in likelihood for
    ** the most promising candidates. Perform moves that give an improvement.
    */
    Calc_Tree_Length (e_prune, v_prune, tree);
    candidate = Best_Lk_Change (e_prune, v_prune, tree);
      }
    }
  }

  /* If the best candidate has a positive estimated change in
  ** likelihood, perform that move.
  */
  if (optim_cand[0]->delta_lk > 1.0/BIG)
  {
    improvement = 1;
    Make_Move (optim_cand[0],0,tree);
  }

  /*
  ** If there was no improvement at all, try local t_edge length optimization at the
  ** regraft position.
  */
  if (!improvement)
  {
/*     PhyML_Printf ("\n.  - performing local t_edge length optimizations"); */
    if ((candidate = Find_Optim_Local (tree)) >= 0)
    {
      improvement = 1;
      Make_Move (optim_cand[candidate],1,tree);
    }
  }

  /*
  ** If there was still no improvement, try global t_edge length optimization.
  */
  if (!improvement)
  {
/*     PhyML_Printf ("\n.  - performing global t_edge length optimization"); */
    if ((candidate = Find_Optim_Globl (tree)) >= 0)
    {
      improvement = 1;
      Make_Move (optim_cand[candidate],2,tree);
    }
  }

  /*
  ** Optimize all t_edge lengths again to make sure we got an updated
  ** likelihood value.
  */
  Set_Both_Sides(YES,tree);
  cur_lk = Lk(NULL,tree);
  root = tree->a_nodes[0];
  Optimize_Br_Len_Serie (tree);
  Set_Both_Sides(YES,tree);
  cur_lk = Lk(NULL,tree);
  time(&(tree->t_current));
  if(tree->mod->s_opt->print) Print_Lk(tree,"topology");

  /*
  ** Return the result.
  */
  return (improvement);
}


/*
** Perform_One_Moves: Perform a round of SPR moves on the tree. Prune each subtree in
**                    turn and calculate the change in tree length for each candidate
**                    regraft position. Estimate change in likelihood for the most
**                    promising moves, and perform the first move that results in an
**                    improvement. If no improvements were found at all, try local edge
**                    length optimization. If still no improvement, try global edge
**                    length optimization.
**
** Parameters:
**   - tree:     The tree to perform the SPR moves on.
**   - max_size: The maximum size (= number of taxa) of the subtrees to be
**               pruned. If m=0 or m>ntax, all possible prunings will be
**               considered.
**
** Returns:
**   If the current tree could be improved: 1.
**   Otherwise:                             0.
*/

int Perform_One_SPR(t_tree *tree, int max_size)
{
  int   nr_edges, i, j, candidate, improvement;
  t_node *root, *v_prune;
  t_edge *e_prune;

  /*
  ** Calculate the average subtree distances.
  */

  root = tree->a_nodes[0];
  PostOrder_v (tree, root->v[2], root->b[2]);

  /*
  ** Initialize the array of optimization candidates.
  */
  for (i = 0; i < tree->mod->s_opt->wim_n_optim; i++)
  {
    optim_cand[i]->delta_lk = -1.0*BIG;
    optim_cand[i]->d_L = -1.0*BIG;
  }

  /*
  ** Try all possible SPR moves and perform the ones that give an improvement.
  */
  nr_edges = 2*tree->n_otu - 3;
  cur_lk = tree->c_lnL;
  improvement = 0;

/*   PhyML_Printf("\n >>>>>>>>>>>>>>>>>>"); */
/*   PhyML_Printf("\n. cur_lk = %f %f",cur_lk,Lk(tree)); */

/*   PhyML_Printf ("\n. Trying SPR moves"); */
/*   PhyML_Printf ("\n.  - calculating tree distances and estimating likelihoods"); */

  for(i = 0; i < nr_edges; i++)
    {
      /*
      ** Get the next prune edge.
      */
      e_prune = tree->a_edges[i];
      /*
      ** Try right subtree if appropriate.
      */
      if (!e_prune->left->tax)
    {
      /*
      ** Clear the regraft candidate list.
      */
      for (j = 0; j < tree->mod->s_opt->wim_n_rgrft; j++)
        {
          rgrft_cand[j]->d_L = -1.0*BIG;
        }
      v_prune = e_prune->left;
/* 	  if ((max_size == 0) || (e_prune->num_tax_rght <= max_size)) */
        {
          /*
          ** Calculate changes in tree length, and estimate changes in likelihood for
          ** the most promising candidates. Perform moves that give an improvement.
          */
          Calc_Tree_Length (e_prune, v_prune, tree);
          if ((candidate = Est_Lk_Change (e_prune, v_prune, tree)) >= 0)
        {
          improvement = 1;
          Make_Move (rgrft_cand[candidate],0,tree);
        }
        }
    }
      /*
      ** Try left subtree if appropriate.
      */
      if (!e_prune->rght->tax && !improvement)
    {
      /*
      ** Clear the regraft candidate list.
      */
      for (j = 0; j < tree->mod->s_opt->wim_n_rgrft; j++)
        {
          rgrft_cand[j]->d_L = -1.0*BIG;
        }
      v_prune = e_prune->rght;
/* 	  if ((max_size == 0) || (e_prune->num_tax_left <= max_size)) */
        {
          /*
          ** Calculate changes in tree length, and estimate changes in likelihood for
          ** the most promising candidates. Perform moves that give an improvement.
          */
          Calc_Tree_Length (e_prune, v_prune, tree);
          if ((candidate = Est_Lk_Change (e_prune, v_prune, tree)) >= 0)
        {
          improvement = 1;
          Make_Move (rgrft_cand[candidate],0,tree);
/* 		  PhyML_Printf("\n. Make simple move"); */
/* 		  PhyML_Printf("\n. lk after simple move = %f",Lk(tree)); */
        }
        }
    }
      if(improvement) break;
    }

  /*
  ** If there was no improvement at all, try local t_edge length optimization at the
  ** regraft position.
  */

/*   PhyML_Printf("\n. before local = %f %f",tree->c_lnL,Lk(tree)); */

  if (!improvement)
    {
      /*       PhyML_Printf ("\n.  - performing local t_edge length optimizations"); */
      if ((candidate = Find_Optim_Local (tree)) >= 0)
    {
/*  	  PhyML_Printf("\n. make local move"); */
      improvement = 1;
      Make_Move (optim_cand[candidate],1,tree);
/* 	  PhyML_Printf("\n. lk after local move = %f",Lk(tree)); */
    }
    }


  /*
  ** If there was still no improvement, try global t_edge length optimization.
  */
/*   PhyML_Printf("\n. before global = %f %f",tree->c_lnL,Lk(tree)); */

  if (!improvement)
    {
/*       PhyML_Printf ("\n.  - performing global t_edge length optimization"); */
      if ((candidate = Find_Optim_Globl (tree)) >= 0)
    {
/*  	  PhyML_Printf("\n. make global move"); */
      improvement = 1;
      Make_Move (optim_cand[candidate],2,tree);
/* 	  PhyML_Printf("\n. lk after global move = %f",Lk(tree)); */
    }
    }

/*   PhyML_Printf("\n. after all = %f %f",tree->c_lnL,Lk(tree)); */

  /*
  ** Optimize all t_edge lengths again to make sure we got an updated
  ** likelihood value.
  */
  Set_Both_Sides(YES,tree);
  cur_lk = Lk(NULL,tree);
  root = tree->a_nodes[0];
  Optimize_Br_Len_Serie (tree);
  Set_Both_Sides(YES,tree);
  cur_lk = Lk(NULL,tree);
  time(&(tree->t_current));
  if(tree->mod->s_opt->print) Print_Lk(tree,"topology");

  /*
  ** Return the result.
  */
  return (improvement);
}


/*
**  Calc_Tree_Length: Calculate the change in tree length, given a pruned subtree,
**                    for each possible regraft position.
**
** Parameters:
**   - e_prune: The t_edge at which the subtree is pruned.
**   - v_prune: The root of the pruned subtree.
*/

void Calc_Tree_Length (t_edge *e_prune, t_node *v_prune, t_tree *tree)
{
  int     i, d0, d1, d2;
  phydbl  d_uu;
  t_node   *u_prune, *u1, *u2;

  /*
  ** Get the directions from t_node v_prune.
  */
  d0 = -1;
  u_prune = NULL;
  for (i = 0; i < 3; i++)
  {
    if (v_prune->b[i] == e_prune)
    {
      d0 = i;
      u_prune = v_prune->v[i];
      break;
    }
  }
  d1 = (d0 + 1) % 3;
  d2 = 3 - d0 - d1;

  /*
  ** Get the relevant average subtree distance within the pruned subtree.
  */
  if (!u_prune->tax)
  {
    u1 = u2 = NULL;
    for (i = 0; i < 3; i++)
    {
      if (u_prune->b[i] != e_prune)
      {
    if (u1 == NULL)
    {
      u1 = u_prune->v[i];
    }
    else
    {
      u2 = u_prune->v[i];
    }
      }
    }
    d_uu = subtree_dist[u1->num][u2->num];
  }
  else
  {
    d_uu = 0.0;
  }

  /*
  ** Recursively calculate the change in tree length for each
  ** possible regraft position.
  **
  ** First recurse into direction d1.
  */
  if (!v_prune->v[d1]->tax)
  {
    u1 = u2 = NULL;
    for (i = 0; i < 3; i++)
    {
      if (v_prune->v[d1]->b[i] != v_prune->b[d1])
      {
    if (u1 == NULL)
    {
      u1 = v_prune->v[d1]->v[i];
    }
    else
    {
      u2 = v_prune->v[d1]->v[i];
    }
      }
    }
    Tree_Length(v_prune, u_prune, v_prune->v[d1], v_prune->v[d2], u1, v_prune->v[d2],
        u2, subtree_dist[u_prune->num][v_prune->v[d2]->num], d_uu, 0.0, 1, tree);
    Tree_Length(v_prune, u_prune, v_prune->v[d1], v_prune->v[d2], u2, v_prune->v[d2],
        u1, subtree_dist[u_prune->num][v_prune->v[d2]->num], d_uu, 0.0, 1, tree);
  }
  /*
  ** Next recurse into direction d2.
  */
  if (!v_prune->v[d2]->tax)
  {
    u1 = u2 = NULL;
    for (i = 0; i < 3; i++)
    {
      if (v_prune->v[d2]->b[i] != v_prune->b[d2])
      {
    if (u1 == NULL)
    {
      u1 = v_prune->v[d2]->v[i];
    }
    else
    {
      u2 = v_prune->v[d2]->v[i];
    }
      }
    }
    Tree_Length(v_prune, u_prune, v_prune->v[d2], v_prune->v[d1], u1, v_prune->v[d1],
        u2, subtree_dist[u_prune->num][v_prune->v[d1]->num], d_uu, 0.0, 1, tree);
    Tree_Length(v_prune, u_prune, v_prune->v[d2], v_prune->v[d1], u2, v_prune->v[d1],
        u1, subtree_dist[u_prune->num][v_prune->v[d1]->num], d_uu, 0.0, 1, tree);
  }
}


/*
** Tree_Length: Recursively calculate the change in tree length for a given pruned
**              subtree and regraft position.
**
** Parameters:
**   - v_prune: The root of the pruned subtree.
**   - u_prune: The t_node adjacent to v_p along the pruned edge.
**   - v_n:     The t_node adjacent to the regraft t_edge in the "backward" direction.
**   - v_n_1:   The previous v_n.
**   - v_nx1:   The t_node adjacent to the regrafting t_edge in the "forward" direction.
**   - v_0:     The other t_node originally adjacent to v_p;
**   - u_n:     The other t_node adjecent to v_n (besides v_n_1 and v_nx1);
**   - d_uv_1:  The distance between u_p and v_n_1;
**   - d_uu:    The subtree distance between descendants of u_prune.
**   - d_L_1:   The previous change in tree length.
**   - n:       The current distance from the prune position.
*/

void Tree_Length (t_node *v_prune, t_node *u_prune, t_node *v_n, t_node *v_n_1, t_node *v_nx1,
          t_node *v_0, t_node *u_n, phydbl d_up_v_1, phydbl d_L_1, phydbl d_uu,
          int n, t_tree *tree)
{
  int     i, j;
  phydbl  d_un_v, d_up_v, d_L;
  t_node   *u1, *u2;
  t_edge   *e_prune, *e_regraft;
  _move_ *tmp_cand;

  /*
  ** Update the path and number of calculations.
  */
  path[n] = v_n;
  nr_d_L++;
  e_prune = NULL;
  e_regraft = NULL;

  /*
  ** Calculate the change in tree length for the current pruned subtree and regraft
  ** position.
  */
  if (n == 1)
  {
    d_un_v = subtree_dist[u_n->num][v_0->num];
  }
  else
  {
    d_un_v = subtree_dist[u_n->num][v_n_1->num] -
        (pow (0.5, n) * subtree_dist[u_n->num][u_prune->num]) +
        (pow (0.5, n) * subtree_dist[u_n->num][v_0->num]);
  }
  d_up_v = 0.5 * (d_up_v_1 + subtree_dist[u_prune->num][u_n->num]);
  /*
  ** Alternative method for calculating d_up_v. Just kept it around for reference...
  **
  d_up_v = subtree_dist[u_prune->num][u_n->num] - (0.5 * d_uu);
  if (!u_n->tax)
  {
    u1 = u2 = NULL;
    for (i = 0; i < 3; i++)
    {
      if (u_n->v[i] != v_n)
      {
    if (u1 == NULL)
    {
      u1 = u_n->v[i];
    }
    else
    {
      u2 = u_n->v[i];
    }
      }
    }
    d_up_v -= 0.5 * subtree_dist[u1->num][u2->num];
  }
  for (i = 0; i < 3; i++)
  {
    if (u_n->v[i] == v_n)
    {
      d_up_v -= u_n->b[i]->l->v;
      break;
    }
  }
  */
  d_L = d_L_1 + 0.25*((d_up_v_1 + subtree_dist[u_n->num][v_nx1->num]) -
              (d_un_v + subtree_dist[u_prune->num][v_nx1->num]));

  /*
  ** If the change is within the tree->mod->s_opt->wim_n_rgrft best ones so far, save it.
  */
  if (d_L > rgrft_cand[tree->mod->s_opt->wim_n_rgrft-1]->d_L)
  {
    for (i = 0; i < 3; i++)
    {
      if (v_prune->v[i] == u_prune)
      {
    e_prune = v_prune->b[i];
      }
      if (v_n->v[i] == v_nx1)
      {
    e_regraft = v_n->b[i];
      }
    }
    i = tree->mod->s_opt->wim_n_rgrft-1;
    rgrft_cand[i]->v_prune = v_prune;
    rgrft_cand[i]->u_prune = u_prune;
    rgrft_cand[i]->v_n = v_n;
    rgrft_cand[i]->v_nx1 = v_nx1;
    rgrft_cand[i]->u_n = u_n;
    rgrft_cand[i]->e_prune = e_prune;
    rgrft_cand[i]->e_regraft = e_regraft;
    rgrft_cand[i]->d_L = d_L;
    rgrft_cand[i]->d_up_v = d_up_v;
    rgrft_cand[i]->d_un_v = d_un_v;
    rgrft_cand[i]->dist = n;
    for (j = 1; j <= n; j++)
      {
    rgrft_cand[i]->path[j] = path[j];
      }

    rgrft_cand[i]->path[n+1] = v_nx1;
    /*
    ** Move the candidate to the appropriate position in the list, so the list
    ** remains sorted in decreasing d_L value.
    */
    while ((i > 0) && (rgrft_cand[i]->d_L > rgrft_cand[i-1]->d_L))
    {
      tmp_cand = rgrft_cand[i];
      rgrft_cand[i] = rgrft_cand[i-1];
      rgrft_cand[i-1] = tmp_cand;
      i--;
    }
  }

  /*
  ** Recurse.
  */
  if (n < tree->mod->s_opt->wim_max_dist)
  {
    if (!v_nx1->tax)
    {
      u1 = u2 = NULL;
      for (i = 0; i < 3; i++)
      {
    if (v_nx1->v[i] != v_n)
    {
      if (u1 == NULL)
      {
        u1 = v_nx1->v[i];
      }
      else
      {
        u2 = v_nx1->v[i];
      }
    }
      }
      Tree_Length (v_prune, u_prune, v_nx1, v_n, u1, v_0, u2, d_up_v, d_uu, d_L, n+1, tree);
      Tree_Length (v_prune, u_prune, v_nx1, v_n, u2, v_0, u1, d_up_v, d_uu, d_L, n+1, tree);
    }
  }
}


/*
** Est_Lk_Change: Estimate the changes in likelihood for the most promising candidate
**                regraft positions given a pruned subtree.
**
** Parameters:
**   - e_prune: The t_edge at which the subtree was pruned.
**   - v_prune: The root of the pruned subtree.
**   - tree:    The tree on which to do the calculations.
**
** Returns:
**   If an improvement as found: The candidate which gives the improvement (which
**                               will be the first one found).
**   Otherwise:                  -1.
*/

int Est_Lk_Change (t_edge *e_prune, t_node *v_prune, t_tree *tree)
{
  int     i, j, cand, best_cand, d0, d1, d2, n, pat, cat, ste;
  phydbl  d_uu, best_d_lk, l_connect, l_01, l_02, l_12, l_est[3], new_lk,
          l_simple[3], l_dist[3];
  phydbl *p_lk1_tmp, *p_lk2_tmp, *p_lk;
  int *p_sum;
  t_node   *u_prune, *v_n, *v_nx1, *u1, *u2;
  t_edge   *e_regraft, *e_tmp;
  _move_ *tmp_cand;
  int dim1, dim2;


  dim1 = tree->mod->ns * tree->mod->ras->n_catg;
  dim2 = tree->mod->ras->n_catg;

  /*
  ** Get the directions from t_node v_prune.
  */
  d0 = -1;
  u_prune = NULL;
  for (i = 0; i < 3; i++)
  {
    if (v_prune->b[i] == e_prune)
    {
      d0 = i;
      u_prune = v_prune->v[i];
      break;
    }
  }
  d1 = (d0 + 1) % 3;
  d2 = 3 - d0 - d1;

  /*
  ** Copy the relevant partial likelihoods to the temporary regraft structure.
  ** We can point to the original matrices, cos they won't be changed anyway.
  */
  if (v_prune == e_prune->left)
  {
    v_tmp->b[0]->p_lk_rght = e_prune->p_lk_rght;
    v_tmp->b[0]->sum_scale_rght = e_prune->sum_scale_rght;
  }
  else
  {
    v_tmp->b[0]->p_lk_rght = e_prune->p_lk_left;
    v_tmp->b[0]->sum_scale_rght = e_prune->sum_scale_left;
  }
  v_tmp->num = v_prune->num;
  v_tmp->v[0]->num = u_prune->num;
  v_tmp->b[0]->num = e_prune->num;

  /*
  ** Estimate the length of the t_edge that will connect the two "detached" nodes
  ** after pruning. (The average of the sum of the lengths of the original two
  ** edges and the average subtree distance based estimate.)
  */
  l_connect = subtree_dist[v_prune->v[d1]->num][v_prune->v[d2]->num];
  if (!v_prune->v[d1]->tax)
  {
    u1 = u2 = NULL;
    for (i = 0; i < 3; i++)
    {
      if (v_prune->v[d1]->b[i] != v_prune->b[d1])
      {
    if (u1 == NULL)
    {
      u1 = v_prune->v[d1]->v[i];
    }
    else
    {
      u2 = v_prune->v[d1]->v[i];
    }
      }
    }
    l_connect -= 0.5 * subtree_dist[u1->num][u2->num];
  }
  if (!v_prune->v[d2]->tax)
  {
    u1 = u2 = NULL;
    for (i = 0; i < 3; i++)
    {
      if (v_prune->v[d2]->b[i] != v_prune->b[d2])
      {
    if (u1 == NULL)
    {
      u1 = v_prune->v[d2]->v[i];
    }
    else
    {
      u2 = v_prune->v[d2]->v[i];
    }
      }
    }
    l_connect -= 0.5 * subtree_dist[u1->num][u2->num];
  }
  l_connect += (v_prune->b[d1]->l->v + v_prune->b[d2]->l->v);
  l_connect /= 2.0;

  /*
  ** Temporarily swap the relevant partial likelihoods at the prune site.
  **
  ** Direction d1.
  */
  if (v_prune == v_prune->b[d1]->left)
  {
    p_lk1_tmp = v_prune->b[d1]->p_lk_left;
    if (v_prune == v_prune->b[d2]->left)
    {
      v_prune->b[d1]->p_lk_left = v_prune->b[d2]->p_lk_rght;
    }
    else
    {
      v_prune->b[d1]->p_lk_left = v_prune->b[d2]->p_lk_left;
    }
  }
  else
  {
    p_lk1_tmp = v_prune->b[d1]->p_lk_rght;
    if (v_prune == v_prune->b[d2]->left)
    {
      v_prune->b[d1]->p_lk_rght = v_prune->b[d2]->p_lk_rght;
    }
    else
    {
      v_prune->b[d1]->p_lk_rght = v_prune->b[d2]->p_lk_left;
    }
  }
  /*
  ** Direction d2.
  */
  if (v_prune == v_prune->b[d2]->left)
  {
    p_lk2_tmp = v_prune->b[d2]->p_lk_left;
    if (v_prune == v_prune->b[d1]->left)
    {
      v_prune->b[d2]->p_lk_left = v_prune->b[d1]->p_lk_rght;
    }
    else
    {
      v_prune->b[d2]->p_lk_left = v_prune->b[d1]->p_lk_left;
    }
  }
  else
  {
    p_lk2_tmp = v_prune->b[d2]->p_lk_rght;
    if (v_prune == v_prune->b[d1]->left)
    {
      v_prune->b[d2]->p_lk_rght = v_prune->b[d1]->p_lk_rght;
    }
    else
    {
      v_prune->b[d2]->p_lk_rght = v_prune->b[d1]->p_lk_left;
    }
  }

  /*
  ** Temporarily set the t_edge lengths and update transition prob's at the
  ** prune site.
  */
  v_prune->b[d1]->l_old->v = v_prune->b[d1]->l->v;
  v_prune->b[d2]->l_old->v = v_prune->b[d2]->l->v;
  v_prune->b[d1]->l->v = l_connect;
  v_prune->b[d2]->l->v = l_connect;
  Update_PMat_At_Given_Edge (v_prune->b[d1], tree);
  Update_PMat_At_Given_Edge (v_prune->b[d2], tree);

  /*
  ** Get the relevant average subtree distance within the pruned subtree.
  */
  if (!u_prune->tax)
  {
    u1 = u2 = NULL;
    for (i = 0; i < 3; i++)
    {
      if (u_prune->b[i] != e_prune)
      {
    if (u1 == NULL)
    {
      u1 = u_prune->v[i];
    }
    else
    {
      u2 = u_prune->v[i];
    }
      }
    }
    d_uu = subtree_dist[u1->num][u2->num];
  }
  else
  {
    d_uu = 0.0;
  }

  /*
  ** Try each candidate SPR and estimate the change in likelihood.
  */
  best_d_lk = 1.0/BIG;
  best_cand = -1;
  for (cand = 0; cand < tree->mod->s_opt->wim_n_rgrft; cand++)
  {
    /*
    ** If there are no more candidates, bail out...
    */
    if (FABS(rgrft_cand[cand]->d_L - 1.0*BIG) < SMALL)
    {
      break;
    }
    else
    {
      nr_d_lk++;
    }

    /*
    ** Get the relevant nodes and edges.
    */
    v_n = rgrft_cand[cand]->v_n;
    v_nx1 = rgrft_cand[cand]->v_nx1;
    e_regraft = rgrft_cand[cand]->e_regraft;

    /*
    ** Update the relevant partial likelihoods along the path between the prune
    ** and regraft positions (temporarily save the first one).
    */
    n = rgrft_cand[cand]->dist;
    e_tmp = NULL;
    p_lk = NULL;
    p_sum = NULL;
    for (i = 1; i <= n; i++)
    {
      /*
      ** Get the next t_edge along the path.
      */
      for (j = 0; j < 3; j++)
      {
    if (rgrft_cand[cand]->path[i]->v[j] == rgrft_cand[cand]->path[i+1])
    {
      e_tmp = rgrft_cand[cand]->path[i]->b[j];
      break;
    }
      }
      if (i == 1)
      {
    /*
    ** Save the first partial likelihood along the path.
    */
    if (rgrft_cand[cand]->path[i] == e_tmp->left)
      {
        p_lk  = e_tmp->p_lk_left;
        p_sum = e_tmp->sum_scale_left;
      }
    else
      {
        p_lk  = e_tmp->p_lk_rght;
        p_sum = e_tmp->sum_scale_rght;
      }

    for (pat = 0; pat < tree->n_pattern; pat++)
      {
        sum_scale_tmp[pat] = p_sum[pat];
        for (cat = 0; cat < tree->mod->ras->n_catg; cat++)
          {
        for (ste = 0; ste < tree->mod->ns; ste++)
          {
            p_lk_tmp[pat*dim1+cat*dim2+ste] = p_lk[pat*dim1+cat*dim2+ste];
          }
          }
      }
      }
      Update_P_Lk (tree, e_tmp, rgrft_cand[cand]->path[i]);
    }
    if (v_n == e_regraft->left)
    {
      v_tmp->b[1]->p_lk_rght = e_regraft->p_lk_left;
      v_tmp->b[2]->p_lk_rght = e_regraft->p_lk_rght;
      v_tmp->b[1]->sum_scale_rght = e_regraft->sum_scale_left;
      v_tmp->b[2]->sum_scale_rght = e_regraft->sum_scale_rght;
    }
    else
    {
      v_tmp->b[1]->p_lk_rght = e_regraft->p_lk_rght;
      v_tmp->b[2]->p_lk_rght = e_regraft->p_lk_left;
      v_tmp->b[1]->sum_scale_rght = e_regraft->sum_scale_rght;
      v_tmp->b[2]->sum_scale_rght = e_regraft->sum_scale_left;
    }

    /*
    ** Estimate t_edge lengths of the three relevant regraft edges based on
    ** average subtree distances.
    **
    ** l_01
    */
    /*
    ** Alternative method of estimating l_01. Kept it around for reference...
    **
    l_01 = subtree_dist[u_prune->num][u_n->num] - (0.5 * d_uu);
    if (!u_n->tax)
    {
      u1 = u2 = NULL;
      for (i = 0; i < 3; i++)
      {
    if (u_n->v[i] != v_n)
    {
      if (u1 == NULL)
      {
        u1 = u_n->v[i];
      }
      else
      {
        u2 = u_n->v[i];
      }
    }
      }
      l_01 -= 0.5 * subtree_dist[u1->num][u2->num];
    }
    for (i = 0; i < 3; i++)
    {
      if (u_n->v[i] == v_n)
      {
    l_01 -= u_n->b[i]->l->v;
    break;
      }
    }
    */
    l_01 = rgrft_cand[cand]->d_up_v - (0.5 * rgrft_cand[cand]->d_un_v) -
           (0.5 * d_uu);
    /*
    ** l_02
    */
    l_02 = subtree_dist[u_prune->num][v_nx1->num] - (0.5 * d_uu);
    if (!v_nx1->tax)
    {
      u1 = u2 = NULL;
      for (i = 0; i < 3; i++)
      {
    if (v_nx1->v[i] != v_n)
    {
      if (u1 == NULL)
      {
        u1 = v_nx1->v[i];
      }
      else
      {
        u2 = v_nx1->v[i];
      }
    }
      }
      l_02 -= (0.5 * subtree_dist[u1->num][u2->num]);
    }
    /*
    ** l_12
    */
    l_12 = e_regraft->l->v;
    /*
    ** Simple estimates.
    */
    l_simple[0] = l_02 - (0.5*e_regraft->l->v);
    l_simple[1] = 0.5 * e_regraft->l->v;
    l_simple[2] = 0.5 * e_regraft->l->v;
    /*
    ** Average subtree distance based estimates.
    */
    l_dist[0] = 0.5 * ( l_01 + l_02 - l_12);
    l_dist[1] = 0.5 * ( l_01 - l_02 + l_12);
    l_dist[2] = 0.5 * (-l_01 + l_02 + l_12);
    /*
    ** Take the average of the two estimates.
    */
    l_est[0] = (l_simple[0] + l_dist[0]) / 2.0;
    l_est[1] = (l_simple[1] + l_dist[1]) / 2.0;
    l_est[2] = (l_simple[2] + l_dist[2]) / 2.0;

    /*
    ** Set the t_edge lengths and update the relevant transition prob's and
    ** partial likelihoods in the temporary regraft structure.
    */

    for (i = 0; i < 3; i++)
    {
      v_tmp->b[i]->l->v = l_est[i]; /* TO DO */
      Update_PMat_At_Given_Edge (v_tmp->b[i], tree);
    }

    /* Beg SG 18 May 2007 */
    if(tree->mod->s_opt->wim_inside_opt)
      {
    Triple_Dist(v_tmp,tree,0);
    For(i,3) l_est[i] = v_tmp->b[i]->l->v;
      }
    /* End SG 18 May 2007 */


    /*
    ** Calculate the change in likelihood locally. Save it and the estimated edge
    ** lengths in the current candidate in the list.
    */
    Update_P_Lk (tree, v_tmp->b[0], v_tmp);
    new_lk = Lk(v_tmp->b[0],tree);
/*     PhyML_Printf("\n. new_lk = %f",new_lk); */

    rgrft_cand[cand]->delta_lk = new_lk - cur_lk;
    rgrft_cand[cand]->rgrft_rank = cand;
    rgrft_cand[cand]->optim_rank = -1;
    rgrft_cand[cand]->globl_rank = -1;
    rgrft_cand[cand]->l_connect = l_connect;
    for (i = 0; i < 3; i++)
    {
      rgrft_cand[cand]->l_est[i] = v_tmp->b[i]->l->v;
    }
    if (rgrft_cand[cand]->delta_lk > best_d_lk)
    {
      best_d_lk = rgrft_cand[cand]->delta_lk;
      best_cand = cand;
    }

    /*
    ** If the change is within the tree->mod->s_opt->wim_n_optim best ones, save it in the list of
    ** optimization candidates.
    */
    if (rgrft_cand[cand]->delta_lk > optim_cand[tree->mod->s_opt->wim_n_optim-1]->delta_lk)
    {
      i = tree->mod->s_opt->wim_n_optim-1;
      optim_cand[i]->v_prune = rgrft_cand[cand]->v_prune;
      optim_cand[i]->u_prune = rgrft_cand[cand]->u_prune;
      optim_cand[i]->v_n = rgrft_cand[cand]->v_n;
      optim_cand[i]->v_nx1 = rgrft_cand[cand]->v_nx1;
      optim_cand[i]->u_n = rgrft_cand[cand]->u_n;
      optim_cand[i]->e_prune = rgrft_cand[cand]->e_prune;
      optim_cand[i]->e_regraft = rgrft_cand[cand]->e_regraft;
      optim_cand[i]->d_L = rgrft_cand[cand]->d_L;
      optim_cand[i]->dist = rgrft_cand[cand]->dist;
      optim_cand[i]->rgrft_rank = rgrft_cand[cand]->rgrft_rank;
      optim_cand[i]->optim_rank = rgrft_cand[cand]->optim_rank;
      optim_cand[i]->globl_rank = rgrft_cand[cand]->globl_rank;
      optim_cand[i]->l_connect = rgrft_cand[cand]->l_connect;
      for (j = 0; j < 3; j++)
      {
    optim_cand[i]->l_est[j] = rgrft_cand[cand]->l_est[j];
      }
      optim_cand[i]->delta_lk = rgrft_cand[cand]->delta_lk;
      /*
      ** Move the candidate to the appropriate position in the list, so the list
      ** remains sorted in decreasing delta_Lk value.
      */
      while ((i > 0) && (optim_cand[i]->delta_lk > optim_cand[i-1]->delta_lk))
      {
    tmp_cand = optim_cand[i];
    optim_cand[i] = optim_cand[i-1];
    optim_cand[i-1] = tmp_cand;
    i--;
      }
    }

    /*
    ** Reset the partial likelihoods along the path.
    */
    for (pat = 0; pat < tree->n_pattern; pat++)
    {
      p_sum[pat] = sum_scale_tmp[pat];
      for (cat = 0; cat < tree->mod->ras->n_catg; cat++)
      {
    for (ste = 0; ste < tree->mod->ns; ste++)
    {
      p_lk[pat*dim1+cat*dim2+ste] = p_lk_tmp[pat*dim1+cat*dim2+ste];
    }
      }
    }
    n = rgrft_cand[cand]->dist;
    for (i = 2; i <= n; i++)
    {
      for (j = 0; j < 3; j++)
      {
    if (rgrft_cand[cand]->path[i]->v[j] == rgrft_cand[cand]->path[i+1])
    {
      e_tmp = rgrft_cand[cand]->path[i]->b[j];
      break;
    }
      }
      Update_P_Lk (tree, e_tmp, rgrft_cand[cand]->path[i]);
    }

    /*
    ** If an improvement was found, forget the other candidates...
    */
    if (best_cand >= 0)
    {
      break;
    }
  }

  /*
  ** Swap back the relevant partial likelihoods at the prune site.
  */
  if (v_prune == v_prune->b[d1]->left)
  {
    v_prune->b[d1]->p_lk_left = p_lk1_tmp;
  }
  else
  {
    v_prune->b[d1]->p_lk_rght = p_lk1_tmp;
  }
  if (v_prune == v_prune->b[d2]->left)
  {
    v_prune->b[d2]->p_lk_left = p_lk2_tmp;
  }
  else
  {
    v_prune->b[d2]->p_lk_rght = p_lk2_tmp;
  }

  /*
  ** Reset the relevant t_edge lengths and transition prob's at the prune site.
  */
  v_prune->b[d1]->l->v = v_prune->b[d1]->l_old->v;
  v_prune->b[d2]->l->v = v_prune->b[d2]->l_old->v;
  Update_PMat_At_Given_Edge (v_prune->b[d1], tree);
  Update_PMat_At_Given_Edge (v_prune->b[d2], tree);

  /*
  ** Return the best candidate.
  */
  return (best_cand);
}


/*
** Best_Lk_Change: Estimate the changes in likelihood for the most promising candidate
**                 regraft positions given a pruned subtree and save the best one.
**
** Parameters:
**   - e_prune: The t_edge at which the subtree was pruned.
**   - v_prune: The root of the pruned subtree.
**   - tree:    The tree on which to do the calculations.
**
** Returns:
**   The candidate which gives the best (possibly negative) improvement.
*/

int Best_Lk_Change (t_edge *e_prune, t_node *v_prune, t_tree *tree)
{
  int     i, j, cand, best_cand, d0, d1, d2, n, pat, cat, ste;
  phydbl  d_uu, best_d_lk, l_connect, l_01, l_02, l_12, l_est[3], new_lk, l_simple[3], l_dist[3];
  phydbl *p_lk1_tmp, *p_lk2_tmp, *p_lk;
  int *p_sum;
  t_node   *u_prune, *v_n, *v_nx1, *u1, *u2;
  t_edge   *e_regraft, *e_tmp;
  _move_ *tmp_cand;
  int dim1, dim2;

  dim1 = tree->mod->ns * tree->mod->ras->n_catg;
  dim2 = tree->mod->ns ;

  /*
  ** Get the directions from t_node v_prune.
  */
  d0 = -1;
  u_prune = NULL;
  for (i = 0; i < 3; i++)
  {
    if (v_prune->b[i] == e_prune)
    {
      d0 = i;
      u_prune = v_prune->v[i];
      break;
    }
  }
  d1 = (d0 + 1) % 3;
  d2 = 3 - d0 - d1;

  /*
  ** Copy the relevant partial likelihoods to the temporary regraft structure.
  ** We can point to the original matrices, cos they won't be changed anyway.
  */
  if (v_prune == e_prune->left)
  {
    v_tmp->b[0]->p_lk_rght = e_prune->p_lk_rght;
    v_tmp->b[0]->sum_scale_rght = e_prune->sum_scale_rght;
  }
  else
  {
    v_tmp->b[0]->p_lk_rght = e_prune->p_lk_left;
    v_tmp->b[0]->sum_scale_rght = e_prune->sum_scale_left;
  }
  v_tmp->num = v_prune->num;
  v_tmp->v[0]->num = u_prune->num;
  v_tmp->b[0]->num = e_prune->num;

  /*
  ** Estimate the length of the t_edge that will connect the two "detached" nodes
  ** after pruning. (The average of the sum of the lengths of the original two
  ** edges and the average subtree distance based estimate.)
  */
  l_connect = subtree_dist[v_prune->v[d1]->num][v_prune->v[d2]->num];
  if (!v_prune->v[d1]->tax)
  {
    u1 = u2 = NULL;
    for (i = 0; i < 3; i++)
    {
      if (v_prune->v[d1]->b[i] != v_prune->b[d1])
      {
    if (u1 == NULL)
    {
      u1 = v_prune->v[d1]->v[i];
    }
    else
    {
      u2 = v_prune->v[d1]->v[i];
    }
      }
    }
    l_connect -= 0.5 * subtree_dist[u1->num][u2->num];
  }
  if (!v_prune->v[d2]->tax)
  {
    u1 = u2 = NULL;
    for (i = 0; i < 3; i++)
    {
      if (v_prune->v[d2]->b[i] != v_prune->b[d2])
      {
    if (u1 == NULL)
    {
      u1 = v_prune->v[d2]->v[i];
    }
    else
    {
      u2 = v_prune->v[d2]->v[i];
    }
      }
    }
    l_connect -= 0.5 * subtree_dist[u1->num][u2->num];
  }
  l_connect += (v_prune->b[d1]->l->v + v_prune->b[d2]->l->v);
  l_connect /= 2.0;

  /*
  ** Temporarily swap the relevant partial likelihoods at the prune site.
  **
  ** Direction d1.
  */
  if (v_prune == v_prune->b[d1]->left)
  {
    p_lk1_tmp = v_prune->b[d1]->p_lk_left;
    if (v_prune == v_prune->b[d2]->left)
    {
      v_prune->b[d1]->p_lk_left = v_prune->b[d2]->p_lk_rght;
    }
    else
    {
      v_prune->b[d1]->p_lk_left = v_prune->b[d2]->p_lk_left;
    }
  }
  else
  {
    p_lk1_tmp = v_prune->b[d1]->p_lk_rght;
    if (v_prune == v_prune->b[d2]->left)
    {
      v_prune->b[d1]->p_lk_rght = v_prune->b[d2]->p_lk_rght;
    }
    else
    {
      v_prune->b[d1]->p_lk_rght = v_prune->b[d2]->p_lk_left;
    }
  }
  /*
  ** Direction d2.
  */
  if (v_prune == v_prune->b[d2]->left)
  {
    p_lk2_tmp = v_prune->b[d2]->p_lk_left;
    if (v_prune == v_prune->b[d1]->left)
    {
      v_prune->b[d2]->p_lk_left = v_prune->b[d1]->p_lk_rght;
    }
    else
    {
      v_prune->b[d2]->p_lk_left = v_prune->b[d1]->p_lk_left;
    }
  }
  else
  {
    p_lk2_tmp = v_prune->b[d2]->p_lk_rght;
    if (v_prune == v_prune->b[d1]->left)
    {
      v_prune->b[d2]->p_lk_rght = v_prune->b[d1]->p_lk_rght;
    }
    else
    {
      v_prune->b[d2]->p_lk_rght = v_prune->b[d1]->p_lk_left;
    }
  }

  /*
  ** Temporarily set the t_edge lengths and update transition prob's at the
  ** prune site.
  */
  v_prune->b[d1]->l_old->v = v_prune->b[d1]->l->v;
  v_prune->b[d2]->l_old->v = v_prune->b[d2]->l->v;
  v_prune->b[d1]->l->v = l_connect;
  v_prune->b[d2]->l->v = l_connect;
  Update_PMat_At_Given_Edge (v_prune->b[d1], tree);
  Update_PMat_At_Given_Edge (v_prune->b[d2], tree);

  /*
  ** Get the relevant average subtree distance within the pruned subtree.
  */
  if (!u_prune->tax)
  {
    u1 = u2 = NULL;
    for (i = 0; i < 3; i++)
    {
      if (u_prune->b[i] != e_prune)
      {
    if (u1 == NULL)
    {
      u1 = u_prune->v[i];
    }
    else
    {
      u2 = u_prune->v[i];
    }
      }
    }
    d_uu = subtree_dist[u1->num][u2->num];
  }
  else
  {
    d_uu = 0.0;
  }

  /*
  ** Try the best candidate SPRs and estimate the change in likelihood.
  */
  best_d_lk = -1.0*BIG;
  best_cand = 0;
  for (cand = 0; cand < tree->mod->s_opt->wim_n_best; cand++)
  {
    /*
    ** If there are no more candidates, bail out...
    */
    if (FABS(rgrft_cand[cand]->d_L - 1.0*BIG) < SMALL)
    {
      break;
    }
    else
    {
      nr_d_lk++;
    }

    /*
    ** Get the relevant nodes and edges.
    */
    v_n = rgrft_cand[cand]->v_n;
    v_nx1 = rgrft_cand[cand]->v_nx1;
    e_regraft = rgrft_cand[cand]->e_regraft;

    /*
    ** Update the relevant partial likelihoods along the path between the prune
    ** and regraft positions (temporarily save the first one).
    */
    n = rgrft_cand[cand]->dist;
    e_tmp = NULL;
    p_lk = NULL;
    p_sum = NULL;
    for (i = 1; i <= n; i++)
    {
      /*
      ** Get the next t_edge along the path.
      */
      for (j = 0; j < 3; j++)
      {
    if (rgrft_cand[cand]->path[i]->v[j] == rgrft_cand[cand]->path[i+1])
    {
      e_tmp = rgrft_cand[cand]->path[i]->b[j];
      break;
    }
      }
      if (i == 1)
      {
    /*
    ** Save the first partial likelihood along the path.
    */
    if (rgrft_cand[cand]->path[i] == e_tmp->left)
    {
      p_lk = e_tmp->p_lk_left;
      p_sum = e_tmp->sum_scale_left;
    }
    else
    {
      p_lk = e_tmp->p_lk_rght;
      p_sum = e_tmp->sum_scale_rght;
    }
    for (pat = 0; pat < tree->n_pattern; pat++)
    {
      sum_scale_tmp[pat] = p_sum[pat];
      for (cat = 0; cat < tree->mod->ras->n_catg; cat++)
      {
        for (ste = 0; ste < tree->mod->ns; ste++)
        {
          p_lk_tmp[pat*dim1+cat*dim2+ste] = p_lk[pat*dim1+cat*dim2+ste];
        }
      }
    }
      }
      Update_P_Lk (tree, e_tmp, rgrft_cand[cand]->path[i]);
    }
    if (v_n == e_regraft->left)
    {
      v_tmp->b[1]->p_lk_rght = e_regraft->p_lk_left;
      v_tmp->b[2]->p_lk_rght = e_regraft->p_lk_rght;
      v_tmp->b[1]->sum_scale_rght = e_regraft->sum_scale_left;
      v_tmp->b[2]->sum_scale_rght = e_regraft->sum_scale_rght;
    }
    else
    {
      v_tmp->b[1]->p_lk_rght = e_regraft->p_lk_rght;
      v_tmp->b[2]->p_lk_rght = e_regraft->p_lk_left;
      v_tmp->b[1]->sum_scale_rght = e_regraft->sum_scale_rght;
      v_tmp->b[2]->sum_scale_rght = e_regraft->sum_scale_left;
    }

    /*
    ** Estimate t_edge lengths of the three relevant regraft edges based on
    ** average subtree distances.
    **
    ** l_01
    */
    /*
    ** Alternative method of estimating l_01. Kept it around for reference...
    **
    l_01 = subtree_dist[u_prune->num][u_n->num] - (0.5 * d_uu);
    if (!u_n->tax)
    {
      u1 = u2 = NULL;
      for (i = 0; i < 3; i++)
      {
    if (u_n->v[i] != v_n)
    {
      if (u1 == NULL)
      {
        u1 = u_n->v[i];
      }
      else
      {
        u2 = u_n->v[i];
      }
    }
      }
      l_01 -= 0.5 * subtree_dist[u1->num][u2->num];
    }
    for (i = 0; i < 3; i++)
    {
      if (u_n->v[i] == v_n)
      {
    l_01 -= u_n->b[i]->l->v;
    break;
      }
    }
    */
    l_01 = rgrft_cand[cand]->d_up_v - (0.5 * rgrft_cand[cand]->d_un_v) -
           (0.5 * d_uu);
    /*
    ** l_02
    */
    l_02 = subtree_dist[u_prune->num][v_nx1->num] - (0.5 * d_uu);
    if (!v_nx1->tax)
    {
      u1 = u2 = NULL;
      for (i = 0; i < 3; i++)
      {
    if (v_nx1->v[i] != v_n)
    {
      if (u1 == NULL)
      {
        u1 = v_nx1->v[i];
      }
      else
      {
        u2 = v_nx1->v[i];
      }
    }
      }
      l_02 -= (0.5 * subtree_dist[u1->num][u2->num]);
    }
    /*
    ** l_12
    */
    l_12 = e_regraft->l->v;
    /*
    ** Simple estimates.
    */
    l_simple[0] = l_02 - (0.5*e_regraft->l->v);
    l_simple[1] = 0.5 * e_regraft->l->v;
    l_simple[2] = 0.5 * e_regraft->l->v;
    /*
    ** Average subtree distance based estimates.
    */
    l_dist[0] = 0.5 * ( l_01 + l_02 - l_12);
    l_dist[1] = 0.5 * ( l_01 - l_02 + l_12);
    l_dist[2] = 0.5 * (-l_01 + l_02 + l_12);
    /*
    ** Take the average of the two estimates.
    */
    l_est[0] = (l_simple[0] + l_dist[0]) / 2.0;
    l_est[1] = (l_simple[1] + l_dist[1]) / 2.0;
    l_est[2] = (l_simple[2] + l_dist[2]) / 2.0;

    /*
    ** Set the t_edge lengths and update the relevant transition prob's and
    ** partial likelihoods in the temporary regraft structure.
    */
    for (i = 0; i < 3; i++)
    {
      v_tmp->b[i]->l->v = l_est[i];
      Update_PMat_At_Given_Edge (v_tmp->b[i], tree);
    }
    Update_P_Lk (tree, v_tmp->b[0], v_tmp);

    /*
    ** Calculate the change in likelihood locally. Save it and the estimated edge
    ** lengths in the current candidate in the list.
    */
    new_lk = Lk(v_tmp->b[0],tree);
    rgrft_cand[cand]->delta_lk = new_lk - cur_lk;
    rgrft_cand[cand]->rgrft_rank = cand;
    rgrft_cand[cand]->optim_rank = -1;
    rgrft_cand[cand]->globl_rank = -1;
    rgrft_cand[cand]->l_connect = l_connect;
    for (i = 0; i < 3; i++)
    {
      rgrft_cand[cand]->l_est[i] = v_tmp->b[i]->l->v;
    }
    if (rgrft_cand[cand]->delta_lk > best_d_lk)
    {
      best_d_lk = rgrft_cand[cand]->delta_lk;
      best_cand = cand;
    }

    /*
    ** Reset the partial likelihoods along the path.
    */
    for (pat = 0; pat < tree->n_pattern; pat++)
    {
      p_sum[pat] = sum_scale_tmp[pat];
      for (cat = 0; cat < tree->mod->ras->n_catg; cat++)
      {
    for (ste = 0; ste < tree->mod->ns; ste++)
    {
      p_lk[pat*dim1+cat*dim2+ste] = p_lk_tmp[pat*dim1+cat*dim2+ste];
    }
      }
    }
    n = rgrft_cand[cand]->dist;
    for (i = 2; i <= n; i++)
    {
      for (j = 0; j < 3; j++)
      {
    if (rgrft_cand[cand]->path[i]->v[j] == rgrft_cand[cand]->path[i+1])
    {
      e_tmp = rgrft_cand[cand]->path[i]->b[j];
      break;
    }
      }
      Update_P_Lk (tree, e_tmp, rgrft_cand[cand]->path[i]);
    }
  }

  /*
  ** If the best candidate is within the tree->mod->s_opt->wim_n_optim best ones, save it in the list of
  ** optimization candidates.
  */
  if (rgrft_cand[best_cand]->delta_lk > optim_cand[tree->mod->s_opt->wim_n_optim-1]->delta_lk)
    {
      i = tree->mod->s_opt->wim_n_optim-1;
      optim_cand[i]->v_prune = rgrft_cand[best_cand]->v_prune;
      optim_cand[i]->u_prune = rgrft_cand[best_cand]->u_prune;
      optim_cand[i]->v_n = rgrft_cand[best_cand]->v_n;
      optim_cand[i]->v_nx1 = rgrft_cand[best_cand]->v_nx1;
      optim_cand[i]->u_n = rgrft_cand[best_cand]->u_n;
      optim_cand[i]->e_prune = rgrft_cand[best_cand]->e_prune;
      optim_cand[i]->e_regraft = rgrft_cand[best_cand]->e_regraft;
      optim_cand[i]->d_L = rgrft_cand[best_cand]->d_L;
      optim_cand[i]->dist = rgrft_cand[best_cand]->dist;
      optim_cand[i]->rgrft_rank = rgrft_cand[best_cand]->rgrft_rank;
      optim_cand[i]->optim_rank = rgrft_cand[best_cand]->optim_rank;
      optim_cand[i]->globl_rank = rgrft_cand[best_cand]->globl_rank;
      optim_cand[i]->l_connect = rgrft_cand[best_cand]->l_connect;

      for (j = 0; j < 3; j++)
    {
      optim_cand[i]->l_est[j] = rgrft_cand[best_cand]->l_est[j];
    }
      optim_cand[i]->delta_lk = rgrft_cand[best_cand]->delta_lk;
      /*
      ** Move the candidate to the appropriate position in the list, so the list
      ** remains sorted in decreasing delta_Lk value.
      */
      while ((i > 0) && (optim_cand[i]->delta_lk > optim_cand[i-1]->delta_lk))
    {
      tmp_cand = optim_cand[i];
      optim_cand[i] = optim_cand[i-1];
      optim_cand[i-1] = tmp_cand;
      i--;
    }
    }

  /*
  ** Swap back the relevant partial likelihoods at the prune site.
  */
  if (v_prune == v_prune->b[d1]->left)
    {
      v_prune->b[d1]->p_lk_left = p_lk1_tmp;
    }
  else
    {
      v_prune->b[d1]->p_lk_rght = p_lk1_tmp;
    }
  if (v_prune == v_prune->b[d2]->left)
    {
      v_prune->b[d2]->p_lk_left = p_lk2_tmp;
    }
  else
    {
      v_prune->b[d2]->p_lk_rght = p_lk2_tmp;
    }

  /*
  ** Reset the relevant t_edge lengths and transition prob's at the prune site.
  */
  v_prune->b[d1]->l->v = v_prune->b[d1]->l_old->v;
  v_prune->b[d2]->l->v = v_prune->b[d2]->l_old->v;
  Update_PMat_At_Given_Edge (v_prune->b[d1], tree);
  Update_PMat_At_Given_Edge (v_prune->b[d2], tree);

  /*
  ** Return the best candidate.
  */
  return (best_cand);
}


/*
** Make_Move: Perform an actual SPR move and calculate the new likelihood.
**
** Parameters:
**   - candidate: The candidate move to perform.
**   - tree:      The tree on which to perform the move.
**
*/

void Make_Move (_move_ *move, int type, t_tree *tree)
{
  int     i;
  t_node   *v_prune, *u_prune, *v_n, *root;
  t_edge   *e_prune, *e_regraft, *e_connect, *e_avail;
  phydbl  new_lk;

  /*
  ** Get the relevant nodes and edges.
  */
  v_prune = move->v_prune;
  u_prune = move->u_prune;
  v_n = move->v_n;
  e_prune = move->e_prune;
  e_regraft = move->e_regraft;
  /*   PhyML_Printf ("  making move: %d -> %d (%f)\n", e_prune->num, e_regraft->num, move->delta_lk); */
  /*
  ** Perform the move and set t_edge lengths.
  */
  Prune (e_prune, v_prune, &(e_connect), &(e_avail), tree);
  Regraft (e_regraft, v_prune, e_avail, tree);
  e_connect->l->v = move->l_connect;

  for (i = 0; i < 3; i++)
    {
      if (v_prune->v[i] == u_prune)
    {
      v_prune->b[i]->l->v = move->l_est[0];
    }
      else if (v_prune->v[i] == v_n)
    {
      v_prune->b[i]->l->v = move->l_est[1];
    }
      else
    {
      v_prune->b[i]->l->v = move->l_est[2];
    }
    }


  if(type > 0) /* local or global move */
    {
      Restore_Br_Len(tree);
    }

  /*
  ** Calculate the new likelihood.
  */
  Set_Both_Sides(YES,tree);
  new_lk = Lk(NULL,tree);

  if(tree->c_lnL < cur_lk-tree->mod->s_opt->min_diff_lk_local)
    {
      PhyML_Printf("\n== tree->c_lnL = %f cur_lk = %f",tree->c_lnL,cur_lk);
      PhyML_Printf("\n== Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
      Warn_And_Exit("");
    }
  cur_lk = new_lk;

  /*
  ** Recalculate the average distances between all (non-overlapping) subtrees.
  */
  root = tree->a_nodes[0];
  PostOrder_v (tree, root->v[2], root->b[2]);
}


/*
** Find_Optim_Local: Perform local t_edge length optimization on the candidates in the
**                   optimization list, and return the first one that gives an
**                   improvement in likelihood.
**
** Parameters:
**   - tree: The tree on which to check the moves.
**
** Returns:
**   If an improvement was found: The candidate that gives the improvement.
**   Otherwise:                   -1.
*/

int Find_Optim_Local (t_tree *tree)
{
  int     best_cand, cand, i;
  t_node   *v_prune, *u_prune, *v_n;
  t_edge   *e_prune, *e_regraft, *e_connect, *e_avail;
  phydbl  max_change, new_lk;
  _move_ *move, *tmp_cand;

  /*
  ** Try all candidate moves starting from the first one.
  */
  best_cand = -1;
  max_change = 1.0/BIG;
  for(cand = 0; cand < tree->mod->s_opt->wim_n_optim; cand++)
    {
      move = optim_cand[cand];
      if(move->delta_lk > -1.0*BIG)
    {
      /*
      ** Get the relevant nodes and edges.
      */
      v_prune   = move->v_prune;
      u_prune   = move->u_prune;
      v_n       = move->v_n;
      e_prune   = move->e_prune;
      e_regraft = move->e_regraft;

      /*
      ** Perform the move and set t_edge lengths.
      */
      Prune (e_prune, v_prune, &(e_connect), &(e_avail), tree);
      Regraft (e_regraft, v_prune, e_avail, tree);
      e_connect->l_old->v = e_connect->l->v;
      e_connect->l->v = move->l_connect;

      for (i = 0; i < 3; i++)
        {
          v_prune->b[i]->l_old->v = v_prune->b[i]->l->v;

          if (v_prune->v[i] == u_prune)
        {
          v_prune->b[i]->l->v = move->l_est[0];
        }
          else if (v_prune->v[i] == v_n)
        {
          v_prune->b[i]->l->v = move->l_est[1];
        }
          else
        {
          v_prune->b[i]->l->v = move->l_est[2];
        }
        }

      Set_Both_Sides(YES,tree);
      Lk(NULL,tree);  // Not sure anymore whether this is required...

      /*
      ** Use Brent optimization on the relevant edges at the regraft position
      ** and calculate the new likelihood value.
      */
      Br_Len_Brent (0.05,20.,v_prune->b[0], tree);
      Br_Len_Brent (0.05,20.,v_prune->b[1], tree);
      Br_Len_Brent (0.05,20.,v_prune->b[2], tree);

/* 	  Update_PMat_At_Given_Edge (v_prune->b[0], tree); */
/* 	  Update_PMat_At_Given_Edge (v_prune->b[1], tree); */
/* 	  Update_PMat_At_Given_Edge (v_prune->b[2], tree); */

      Update_P_Lk (tree, v_prune->b[0], v_prune);
      new_lk = Lk(v_prune->b[0],tree);

/* 	  PhyML_Printf("\n. local new_lk = %f",new_lk); */

      /*
      ** Save the change in likelihood and move the current candidate to the
      ** appropriate place in the list.
      */
      move->delta_lk = new_lk - cur_lk;
      move->optim_rank = cand;
      i = cand;
      while ((i > 0) && (optim_cand[i]->delta_lk > optim_cand[i-1]->delta_lk))
        {
          tmp_cand = optim_cand[i];
          optim_cand[i] = optim_cand[i-1];
          optim_cand[i-1] = tmp_cand;
          i--;
        }
      if (move->delta_lk > max_change)
        {
          best_cand = i;
          max_change = move->delta_lk;
          Record_Br_Len(tree);
        }


      /*
      ** Undo the move again.
      */
      Prune (e_prune, v_prune, &(e_regraft), &(e_avail), tree);
      Regraft (e_connect, v_prune, e_avail, tree);
      e_regraft->l->v = e_regraft->l_old->v;
      for (i = 0; i < 3; i++)
        {
          v_prune->b[i]->l->v = v_prune->b[i]->l_old->v;
        }
      Set_Both_Sides(YES,tree);
      Lk(NULL,tree);
      nr_loc++;
/* 	  PhyML_Printf("\n. local back to = %f",tree->c_lnL); */
    }
      else
    {
      break;
    }

      /*
      ** If an improvement was found, forget the other candidates...
      */
      if (best_cand >= 0)
    {
      break;
    }
    }

  /*
  ** Return the best candidate.
  */
  return (best_cand);
}


/*
** Find_Optim_Globl: Perform global t_edge length optimization on the candidates in the
**                   optimization list, and return the first one that gives an
**                   improvement in likelihood.
**
** Parameters:
**   - tree: The tree on which to check the moves.
**
** Returns:
**   If an improvement is found: The candidate that gives the improvement.
**   Otherwise:                  -1.
*/

int Find_Optim_Globl (t_tree *tree)
{
  int     best_cand, cand, i;
  t_node   *v_prune, *u_prune, *v_n;
  t_edge   *e_prune, *e_regraft, *e_connect, *e_avail;
  phydbl  max_change, new_lk;
  _move_ *move;

  /*
  ** Try all moves.
  */
  best_cand = -1;
  max_change = 1.0/BIG;
  for (cand = 0; cand < tree->mod->s_opt->wim_n_globl; cand++)
  {
    move = optim_cand[cand];
    if (move->delta_lk > -1.0*BIG)
    {
      Record_Br_Len(tree);

      /*
      ** Get the relevant nodes and edges.
      */
      v_prune = move->v_prune;
      u_prune = move->u_prune;
      v_n = move->v_n;
      e_prune = move->e_prune;
      e_regraft = move->e_regraft;

      /*
      ** Perform the move and optimize all t_edge lengths.
      */
      Prune (e_prune, v_prune, &(e_connect), &(e_avail), tree);
      Regraft (e_regraft, v_prune, e_avail, tree);
      e_connect->l_old->v = e_connect->l->v;
      e_connect->l->v = move->l_connect;

      for (i = 0; i < 3; i++)
    {
      v_prune->b[i]->l_old->v = v_prune->b[i]->l->v;
      if (v_prune->v[i] == u_prune)
        {
          v_prune->b[i]->l->v = move->l_est[0];
        }
      else if (v_prune->v[i] == v_n)
        {
          v_prune->b[i]->l->v = move->l_est[1];
        }
      else
        {
          v_prune->b[i]->l->v = move->l_est[2];
        }
    }

      Set_Both_Sides(YES,tree);
      Lk(NULL,tree);
      Optimize_Br_Len_Serie (tree);
      Set_Both_Sides(YES,tree);
      new_lk = Lk (NULL,tree);

/*       PhyML_Printf("\n. global new_lk = %f\n",tree->c_lnL); */

      /*
      ** Save the change in likelihood and undo the move.
      */

      move->delta_lk = new_lk - cur_lk;
      move->globl_rank = cand;
      if (move->delta_lk > max_change)
    {
      best_cand = cand;
      max_change = move->delta_lk;
      Record_Br_Len(tree);
    }

      Prune (e_prune, v_prune, &(e_regraft), &(e_avail), tree);
      Regraft (e_connect, v_prune, e_avail, tree);
      e_regraft->l->v = e_regraft->l_old->v;
      for (i = 0; i < 3; i++)
    {
      v_prune->b[i]->l->v = v_prune->b[i]->l_old->v;
    }
      Set_Both_Sides(YES,tree);
      Restore_Br_Len(tree);
      Lk(NULL,tree);
      nr_glb++;
/*       PhyML_Printf("\n. global back to = %f",tree->c_lnL); */

    }
    else break;

    /*
    ** If an improvement was found, forget the other candidates...
    */
    if (best_cand >= 0) break;
  }

  /*
  ** Return the best candidate.
  */
  return (best_cand);
}


/*
** Prune: Prune the subtree at a certain t_edge and node. Note that edge
**        lengths are not set and partial likelihoods are not updated.
**
** Parameters:
**   - e:         The t_edge at which to prune the subtree.
**   - v:         The t_node adjacent to t_edge e which forms the root of the subtree.
**   - e_connect: An t_edge pointer which will point to the t_edge that was left
**                after pruning.
**   - e_avail:   The t_edge that is "left over" and should be used in the
**                regrafting step.
**   - tree:      The tree on which the pruning is done.
**
**
**
**          \ /
**           o
**           |
**           | e    --> subtree to be pruned
**           |
**           o v
**          / \
**       e1/   \e2
**        /     \
**       o       o
**      u1       u2   --> such that u1->num < u2->num
*/

void Prune (t_edge *e, t_node *v, t_edge **e_connect, t_edge **e_avail, t_tree *tree)
{
  int     dir0, dir1, dir2, v0, v1, v2, tmp_dir, i, j, k;
  t_node   *u1, *u2, *tmp_node;
  t_edge   *e1, *e2;
  int *sum_scale_f;
  phydbl *p_lk;
  int dim1, dim2;


  dim1 = tree->mod->ns * tree->mod->ras->n_catg;
  dim2 = tree->mod->ns;

  /*
  ** Get the relevant directions, nodes and edges.
  ** Make sure that t_node u1 is the t_node with the smaller number.
  */
  dir0 = -1;
  for (i = 0; i < 3; i++)
  {
    if (v->b[i] == e)
    {
      dir0 = i;
      break;
    }
  }
  dir1 = (dir0 + 1) % 3;
  dir2 = 3 - dir0 - dir1;
  u1 = v->v[dir1];
  u2 = v->v[dir2];
  if (u1->num > u2->num)
  {
    tmp_node = u1;
    u1 = u2;
    u2 = tmp_node;
    tmp_dir = dir1;
    dir1 = dir2;
    dir2 = tmp_dir;
  }
  e1 = v->b[dir1];
  e2 = v->b[dir2];

  /*
  ** Detach t_node v from the tree.
  */
  v->v[dir1] = NULL;
  v->v[dir2] = NULL;
  v->b[dir1] = NULL;
  v->b[dir2] = NULL;

  /*
  ** Connect nodes u1 and u2 via t_edge e1 and copy relevant partial likelihoods.
  */
  if (u2 == e2->left)
  {
    v0 = e2->l_r;
    v1 = e2->l_v1;
    v2 = e2->l_v2;
    sum_scale_f = e2->sum_scale_left;
    p_lk = e2->p_lk_left;
  }
  else
  {
    v0 = e2->r_l;
    v1 = e2->r_v1;
    v2 = e2->r_v2;
    sum_scale_f = e2->sum_scale_rght;
    p_lk = e2->p_lk_rght;
  }
  if (u1 == e1->left)
  {
    e1->rght = u2;
    e1->r_l = v0;
    e1->r_v1 = v1;
    e1->r_v2 = v2;
    for (i = 0; i < tree->n_pattern; i++)
    {
      e1->sum_scale_rght[i] = sum_scale_f[i];
      for (j = 0; j < tree->mod->ras->n_catg; j++)
      {
    for (k = 0; k < tree->mod->ns; k++)
    {
      e1->p_lk_rght[i*dim1+j*dim2+k] = p_lk[i*dim1+j*dim2+k];
    }
      }
    }
  }
  else
  {
    e1->left = u2;
    e1->l_r = v0;
    e1->l_v1 = v1;
    e1->l_v2 = v2;
    for (i = 0; i < tree->n_pattern; i++)
    {
      e1->sum_scale_left[i] = sum_scale_f[i];
      for (j = 0; j < tree->mod->ras->n_catg; j++)
      {
    for (k = 0; k < tree->mod->ns; k++)
    {
      e1->p_lk_left[i*dim1+j*dim2+k] = p_lk[i*dim1+j*dim2+k];
    }
      }
    }
  }
  for (i = 0; i < 3; i++)
  {
    if (u1->v[i] == v)
    {
      u1->v[i] = u2;
    }
    if (u2->v[i] == v)
    {
      u2->v[i] = u1;
      u2->b[i] = e1;
      u2->l[i] = e1->l->v;
    }
  }

  /*
  ** Make sure that a possible tip t_node is still on the right side.
  */
  if (e1->left->tax)
  {
    /*
    ** Swap left and right.
    */
    tmp_node = e1->left;
    e1->left = e1->rght;
    e1->rght = tmp_node;
    tmp_dir = e1->l_r;
    e1->l_r = e1->r_l;
    e1->r_l = tmp_dir;
    tmp_dir = e1->l_v1;
    e1->l_v1 = e1->r_v1;
    e1->r_v1 = tmp_dir;
    tmp_dir = e1->l_v2;
    e1->l_v2 = e1->r_v2;
    e1->r_v2 = tmp_dir;
    p_lk = e1->p_lk_left;
    e1->p_lk_left = e1->p_lk_rght;
    e1->p_lk_rght = p_lk;
    sum_scale_f = e1->sum_scale_left;
    e1->sum_scale_left = e1->sum_scale_rght;
    e1->sum_scale_rght = sum_scale_f;
  }

  /*
  ** Set the connecting and available edges.
  */
  *(e_connect) = e1;
  *(e_avail) = e2;
}


/*
** Regraft: Regraft a subtree at a certain edge. Note that t_edge lengths
**          are not set and partial likelihoods are not updated.
**
** Parameters:
**   - e:     The t_edge to regraft the subtree on.
**   - v:     The root of the subtree to regraft.
**   - avail: A previously deleted t_edge now available for insertion again.
**   - tree:  The tree on which the regrafting is done.
**
**
**      \ /
**       o
**       |
**       |   --> subtree to regraft
**       |
**       o v
**
**   o--------o
**  u1   e    u2   --> such that u1->num < u2->num
*/

void Regraft (t_edge *e, t_node *v, t_edge *avail, t_tree *tree)
{
  int     dir0, dir1, dir2, i, j, k;
  int *sum_scale_f;
  phydbl *p_lk;
  t_node   *u1, *u2;
  int dim1, dim2;

  dim1 = tree->mod->ns * tree->mod->ras->n_catg;
  dim2 = tree->mod->ns ;

  /*
  ** Get the relevant directions and nodes.
  */
  dir0 = -1;
  for (i = 0; i < 3; i++)
  {
    if (v->b[i] != NULL)
    {
      dir0 = i;
      break;
    }
  }
  dir1 = (dir0 + 1) % 3;
  dir2 = 3 - dir0 - dir1;
  if (e->left->num < e->rght->num)
  {
    u1 = e->left;
    u2 = e->rght;
    sum_scale_f = e->sum_scale_rght;
    p_lk = e->p_lk_rght;
  }
  else
  {
    u1 = e->rght;
    u2 = e->left;
    sum_scale_f = e->sum_scale_left;
    p_lk = e->p_lk_left;
  }

  /*
  ** Connect nodes v and u2 via the available t_edge 'avail' and copy the
  ** relevant partial likelihood.
  ** (We want to do this first, cos we need some of the values of edge
  **  e before changing them below).
  */
  avail->left = v;
  avail->rght = u2;
  avail->l_r = dir2;
  avail->l_v1 = dir0;
  avail->l_v2 = dir1;
  v->v[dir2] = u2;
  v->b[dir2] = avail;
  for (i = 0; i < 3; i++)
  {
    if (e == u2->b[i])
    {
      u2->v[i] = v;
      u2->b[i] = avail;
      avail->r_l = i;
      avail->r_v1 = (i + 1) % 3;
      avail->r_v2 = 3 - i - avail->r_v1;
      break;
    }
  }
  for (i = 0; i < tree->n_pattern; i++)
  {
    avail->sum_scale_rght[i] = sum_scale_f[i];
    for (j = 0; j < tree->mod->ras->n_catg; j++)
    {
      for (k = 0; k < tree->mod->ns; k++)
      {
    avail->p_lk_rght[i*dim1+j*dim2+k] = p_lk[i*dim1+j*dim2+k];
      }
    }
  }

  /*
  ** Connect nodes v and u1 via t_edge e.
  */
  if (u1 == e->left)
  {
    e->rght = v;
    e->r_l = dir1;
    e->r_v1 = dir0;
    e->r_v2 = dir2;
  }
  else
  {
    e->left = v;
    e->l_r = dir1;
    e->l_v1 = dir0;
    e->l_v2 = dir2;
  }
  v->v[dir1] = u1;
  v->b[dir1] = e;
  for (i = 0; i < 3; i++)
  {
    if (u1->v[i] == u2)
    {
      u1->v[i] = v;
      break;
    }
  }
}


/*
** PostOrder_v: Recursively visit all nodes v in postorder to calculate
**              the average distance between subtrees.
**
** Parameters:
**   - tree: The tree for which to calculate the average distances.
**   - v:    The current node.
**   - e:    The t_edge we came from.
*/

void PostOrder_v (t_tree *tree, t_node *v, t_edge *e)
{
  int   i;
  t_node *w;

  /*
  ** If v is not a taxon, recurse.
  */
  if (!v->tax)
  {
    for (i = 0; i < 3; i++)
    {
      if (v->b[i] != e)
      {
    PostOrder_v (tree, v->v[i], v->b[i]);
      }
    }
  }

  /*
  ** Recurse over all nodes w not in the current subtree and calculate
  ** the average distance between the current subtree and all others.
  */
  if (v == e->left)
  {
    w = e->rght;
  }
  else
  {
    w = e->left;
  }
  PostOrder_w (tree, v, e, w, e);
}


/*
** PostOrder_w: Recursively visit all nodes w not in the subtree of v in
**              postorder and calculate the average distance between the
**              subtree of v and all others.
**
** Parameters:
**   - tree: The tree for which to calculate the average distances.
**   - v:    The root t_node of the first subtree.
**   - v_e:  The t_edge adjacent to the first subtree.
**   - w:    The current node.
**   - e:    The t_edge we came from.
*/

void PostOrder_w (t_tree *tree, t_node *v, t_edge *v_e, t_node *w, t_edge *e)
{
  int   i;
  t_node *w1, *w2, *v1, *v2;

  /*
  ** If w is not a taxon, recurse.
  */
  if (!w->tax)
  {
    for (i = 0; i < 3; i++)
    {
      if (w->b[i] != e)
      {
    PostOrder_w (tree, v, v_e, w->v[i], w->b[i]);
      }
    }
  }

  /*
  ** Calculate the average distance between the subtrees defined by
  ** nodes v and w.
  */
  if (v->tax && w->tax)
  {
    subtree_dist[v->num][w->num] = seq_dist->dist[v->num][w->num];
  }
  else if (!v->tax)
  {
    v1 = v2 = NULL;
    for (i = 0; i < 3; i++)
    {
      if (v->b[i] != v_e)
      {
    if (v1 == NULL)
    {
      v1 = v->v[i];
    }
    else
    {
      v2 = v->v[i];
    }
      }
    }
    subtree_dist[v->num][w->num] = 0.5*(subtree_dist[v1->num][w->num] +
                    subtree_dist[v2->num][w->num]);
  }
  else
  {
    w1 = w2 = NULL;
    for (i = 0; i < 3; i++)
    {
      if (w->b[i] != e)
      {
    if (w1 == NULL)
    {
      w1 = w->v[i];
    }
    else
    {
      w2 = w->v[i];
    }
      }
    }
    subtree_dist[v->num][w->num] = 0.5*(subtree_dist[v->num][w1->num] +
                    subtree_dist[v->num][w2->num]);
  }
  subtree_dist[w->num][v->num] = subtree_dist[v->num][w->num];
}


/*********************************************************/
/*********************************************************/
/* Sort list of SPR move by putting the shallowest moves first */
void Sort_Spr_List_Depth(t_tree *tree)
{
  int i,j;
  t_spr *buff;

  For(i,tree->n_moves-1)
    {
      for(j=i+1;j<tree->n_moves;j++)
        {
          if(tree->spr_list[j]->depth_path < tree->spr_list[i]->depth_path)
            {
              buff              = tree->spr_list[i];
              tree->spr_list[i] = tree->spr_list[j];
              tree->spr_list[j] = buff;
            }
        }
    }
}

/*********************************************************/
/*********************************************************/
/* Sort list of SPR move by putting the more likely moves first */
void Sort_Spr_List_LnL(t_tree *tree)
{
  int i,j;
  t_spr *buff;

  For(i,tree->size_spr_list-1)
    {
      for(j=i+1;j<tree->size_spr_list;j++)
        {
          if(tree->spr_list[j]->lnL > tree->spr_list[i]->lnL)
            {
              buff              = tree->spr_list[i];
              tree->spr_list[i] = tree->spr_list[j];
              tree->spr_list[j] = buff;
            }
        }
    }
}


/*********************************************************/
/*********************************************************/
/* Sort list of SPR move by putting the more parsimonious moves first */
void Sort_Spr_List_Pars(t_tree *tree)
{
  int i,j;
  t_spr *buff;

  For(i,tree->size_spr_list-1)
    {
      for(j=i+1;j<tree->size_spr_list;j++)
        {
          if(tree->spr_list[j]->pars < tree->spr_list[i]->pars)
            {
              buff              = tree->spr_list[i];
              tree->spr_list[i] = tree->spr_list[j];
              tree->spr_list[j] = buff;
            }
        }
    }
}


/*********************************************************/
/*********************************************************/
/*********************************************************/

void Randomize_Spr_List(t_tree *tree)
{
  int i,j;
  t_spr *buff;

  For(i,tree->n_moves)
    {
      j = Rand_Int(0,tree->n_moves-1);
      buff              = tree->spr_list[i];
      tree->spr_list[i] = tree->spr_list[j];
      tree->spr_list[j] = buff;
    }
}

/*********************************************************/

int Spr(phydbl init_lnL, phydbl prop_spr, t_tree *tree)
{
  int i,br;
  int *br_idx;
  t_edge *b;

  tree->n_improvements = 0;
  tree->max_spr_depth  = 0;

  Set_Both_Sides(YES,tree);

  Reset_Spr_List(tree);

  br_idx = Permutate(2*tree->n_otu-3);

  For(i,MAX(1,(int)((2*tree->n_otu-3)*prop_spr)))
    {
      br = br_idx[i];

      if(!(br%10)) if(tree->io->print_json_trace == YES) JSON_Tree_Io(tree,tree->io->fp_out_json_trace); 

      b = tree->a_edges[br];
      Spr_Subtree(b,b->left,tree);
      Spr_Subtree(b,b->rght,tree);
    }

  Free(br_idx);

  return tree->n_improvements;
}

/*********************************************************/

void Spr_Subtree(t_edge *b, t_node *link, t_tree *tree)
{
  int i;
  int n_moves_pars, n_moves, min_pars, best_move_idx;
  t_spr *best_pars_move;
  t_edge *target, *residual;
  
  best_move_idx = -1;
  tree->n_moves = 0;

  MIXT_Set_Pars_Thresh(tree);

  if((link != b->left) && (link != b->rght))
    {
      PhyML_Printf("\n== Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
      Exit("\n");
    }
  else
    {
      if(!link->tax) Test_All_Spr_Targets(b,link,tree);

      if(tree->n_moves)
        {
          n_moves_pars = MIN(5,tree->n_moves);
          n_moves      = MIN(5,tree->n_moves);

          if(tree->mod->s_opt->spr_lnL == NO)       n_moves = n_moves_pars;
          if(tree->io->fp_in_constraint_tree == NO) n_moves = MAX(1,n_moves);
          
          if(tree->mod->s_opt->spr_pars == YES)
            {
              if(tree->io->fp_in_constraint_tree)
                {
                  PhyML_Printf("\n== Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
                  Exit("\n");
                }

              min_pars = 1E+8;
              best_pars_move = NULL;

              For(i,n_moves)
                if(tree->spr_list[i]->pars < min_pars)
                  {
                    best_pars_move = tree->spr_list[i];
                    min_pars = tree->spr_list[i]->pars;
                  }

              if(best_pars_move->pars < tree->best_pars)
                {
                  Prune_Subtree(best_pars_move->n_link,best_pars_move->n_opp_to_link,&target,&residual,tree);
                  Graft_Subtree(best_pars_move->b_target,best_pars_move->n_link,residual,tree);
                  Set_Both_Sides(YES,tree);
                  Pars(NULL,tree);
                  tree->best_pars = tree->c_pars;
                  if(tree->best_pars != best_pars_move->pars)
                    {
                      PhyML_Printf("\n== best_pars = %d move_pars = %d",tree->best_pars,best_pars_move->pars);
                      PhyML_Printf("\n== Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
                      Exit("\n");
                    }
                  tree->n_improvements++;
                }
              else
                Pars(NULL,tree);
            }
          else
            {
              int apply_move = NO;
              phydbl accept_prob,u;

              if(tree->mod->s_opt->spr_lnL == YES) 
                {
                  Sort_Spr_List_LnL(tree);
                  if(tree->spr_list[0]->lnL > tree->best_lnL)
                    {
                      /* n_moves = 1; */
                      /* best_move_idx = Evaluate_List_Of_Regraft_Pos_Triple(tree->spr_list,n_moves,tree); */
                      best_move_idx = 0;
                    }
                  else
                    {
                      best_move_idx = Evaluate_List_Of_Regraft_Pos_Triple(tree->spr_list,n_moves,tree);
                    }
                }
              else
                {
                  best_move_idx = Evaluate_List_Of_Regraft_Pos_Triple(tree->spr_list,n_moves,tree);
                }

              /* For(i,n_moves) printf("\n. %d %f %d",i,tree->spr_list[i]->lnL,tree->spr_list[i]->pars); */
              /* fflush(NULL); */

              if(best_move_idx > -1)
                {
                  if(Are_Equal(tree->annealing_temp,0.0,1.E-3) == NO)
                    {
                      accept_prob = exp((tree->spr_list[best_move_idx]->lnL - tree->best_lnL)/tree->annealing_temp);
                      u = Uni();
                      if(!(u > accept_prob)) apply_move = YES;
                    }
                  else
                    {
                      if(tree->spr_list[best_move_idx]->lnL > tree->best_lnL + tree->mod->s_opt->min_diff_lk_move)
                        apply_move = YES;
                    }
                }
                
              if((best_move_idx > -1) && (apply_move == YES))
                {
                  Try_One_Spr_Move_Triple(tree->spr_list[best_move_idx],tree);
                }
              else
                {
                  Pars(NULL,tree);
                }
            }
        }
      Reset_Spr_List(tree);
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
      PhyML_Printf("\n== Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
      Exit("\n");
    }

  init_lnL = tree->c_lnL;
  b_target = b_residual = NULL;
  n_opp_to_link  = (n_link == b_pulled->rght)?(b_pulled->left):(b_pulled->rght);

  init_l_pulled = Duplicate_Scalar_Dbl(b_pulled->l);
  init_v_pulled = Duplicate_Scalar_Dbl(b_pulled->l_var);

  dir1 = dir2 = -1;
  For(i,3)
    if(n_link->v[i] != n_opp_to_link)
      {
        if(dir1<0) dir1 = i;
        else       dir2 = i;
      }

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

      if(tree->mod->s_opt->spr_lnL)
        {
          /* Fast_Br_Len(b_target,tree,YES); */
          Update_PMat_At_Given_Edge(b_target,tree);
        }

      best_found = NO;
      tree->depth_curr_path = 0;
      tree->curr_path[0] = b_target->left;
      Test_One_Spr_Target_Recur(b_target->rght,
                                b_target->left,
                                b_pulled,n_link,b_residual,b_target,&best_found,tree);
      
      if(best_found == NO || tree->mod->s_opt->spr_lnL == NO)
        {
          tree->depth_curr_path = 0;
          tree->curr_path[0] = b_target->rght;
          Test_One_Spr_Target_Recur(b_target->left,
                                    b_target->rght,
                                    b_pulled,n_link,b_residual,b_target,&best_found,tree);
        }

      Graft_Subtree(b_target,n_link,b_residual,tree);

      if((n_link->v[dir1] != n_v1) || (n_link->v[dir2] != n_v2)) PhyML_Printf("\n== Warning: -- SWITCH NEEDED -- ! \n");

      Copy_Scalar_Dbl(init_l_v1,n_link->b[dir1]->l);
      Copy_Scalar_Dbl(init_v_v1,n_link->b[dir1]->l_var);

      Copy_Scalar_Dbl(init_l_v2,n_link->b[dir2]->l);
      Copy_Scalar_Dbl(init_v_v2,n_link->b[dir2]->l_var);

      Copy_Scalar_Dbl(init_l_pulled,b_pulled->l);
      Copy_Scalar_Dbl(init_v_pulled,b_pulled->l_var);

      // really useful if spr_lnL == No ?
      Update_PMat_At_Given_Edge(n_link->b[dir1],tree);
      Update_PMat_At_Given_Edge(n_link->b[dir2],tree);
      Update_PMat_At_Given_Edge(b_pulled,tree);

      if(tree->mod->s_opt->spr_lnL == YES)
      	{
          MIXT_Set_Alias_Subpatt(YES,tree);
          Update_P_Lk(tree,b_pulled,  n_link);
          Update_P_Lk(tree,b_target,  n_link);
          Update_P_Lk(tree,b_residual,n_link);
          MIXT_Set_Alias_Subpatt(NO,tree);
        }
      else
      	{
          Update_P_Pars(tree,b_pulled,  n_link);
          Update_P_Pars(tree,b_target,  n_link);
          Update_P_Pars(tree,b_residual,n_link);
        }

      For(i,3)
        if(n_link->v[i] != n_opp_to_link)
          {
            if(tree->mod->s_opt->spr_lnL == YES)
              {
                MIXT_Set_Alias_Subpatt(YES,tree);
                Pre_Order_Lk(n_link,n_link->v[i],tree);
                MIXT_Set_Alias_Subpatt(NO,tree);
              }
            else
              Pre_Order_Pars(n_link,n_link->v[i],tree);
          }
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

void Test_One_Spr_Target_Recur(t_node *a, t_node *d, t_edge *pulled, t_node *link, t_edge *residual, t_edge *init_target, int *best_found, t_tree *tree)
{
  int i;
  
  if(*best_found == YES) return;
  
  if(d->tax) return;
  else
    {
      phydbl move_lnL;
      
      For(i,3)
        {
          if(d->v[i] != a)
            {
              if(tree->mod->s_opt->spr_lnL == YES)
                {
                  MIXT_Set_Alias_Subpatt(YES,tree);
                  Update_P_Lk(tree,d->b[i],d);
                  MIXT_Set_Alias_Subpatt(NO,tree);
                }
              else
                Update_P_Pars(tree,d->b[i],d);
              
              tree->depth_curr_path++;
              tree->curr_path[tree->depth_curr_path] = d->v[i];

              if((tree->depth_curr_path <= tree->mod->s_opt->max_depth_path) &&
                 (tree->depth_curr_path >= tree->mod->s_opt->min_depth_path))
                {
                  move_lnL = Test_One_Spr_Target(d->b[i],pulled,link,residual,init_target,tree);
                  if(move_lnL > tree->best_lnL + tree->mod->s_opt->min_diff_lk_move) 
                    {
                      *best_found = YES;
                      return;
                    }
                }
                            
              if(tree->depth_curr_path < tree->mod->s_opt->max_depth_path)
                Test_One_Spr_Target_Recur(d,d->v[i],pulled,link,residual,init_target,best_found,tree);
              
              tree->depth_curr_path--;
            }
        }
    }
}

/*********************************************************/

phydbl Test_One_Spr_Target(t_edge *b_target, t_edge *b_arrow, t_node *n_link, t_edge *b_residual, t_edge *init_target, t_tree *tree)
{
  scalar_dbl *init_target_l, *init_arrow_l, *init_residual_l;
  scalar_dbl *init_target_v, *init_arrow_v, *init_residual_v;
  int i,dir_v0,dir_v1,dir_v2;
  scalar_dbl *l0,*l1,*l2;
  scalar_dbl *v0,*v1,*v2;
  t_node *n1,*n2;
  phydbl init_lnL, move_lnL;
  int init_pars;
  t_spr *move;

  if(tree->mixt_tree != NULL)
    {
      PhyML_Printf("\n== Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
      Exit("\n");
    }

  tree->n_moves++;

  move_lnL  = UNLIKELY;
  init_lnL  = tree->c_lnL;
  init_pars = tree->c_pars;

  move = tree->spr_list[tree->size_spr_list];

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

  Graft_Subtree(b_target,n_link,b_residual,tree);

  // Save edge lengths so that they can be recovered in the end
  init_target_l   = Duplicate_Scalar_Dbl(b_target->l);
  init_target_v   = Duplicate_Scalar_Dbl(b_target->l_var);

  init_arrow_l    = Duplicate_Scalar_Dbl(b_arrow->l);
  init_arrow_v    = Duplicate_Scalar_Dbl(b_arrow->l_var);

  init_residual_l = Duplicate_Scalar_Dbl(b_residual->l);
  init_residual_v = Duplicate_Scalar_Dbl(b_residual->l_var);

  if(tree->mod->s_opt->spr_lnL == YES)
    {
      MIXT_Set_Alias_Subpatt(YES,tree);
      Update_PMat_At_Given_Edge(b_target,tree);
      Update_PMat_At_Given_Edge(b_arrow,tree);
      Update_P_Lk(tree,b_residual,n_link);
      move_lnL = Lk(b_residual,tree);
      MIXT_Set_Alias_Subpatt(NO,tree);

      /* if(FABS(move_lnL - tree->best_lnL < 10.)) */
      /*   { */
      /*     MIXT_Set_Alias_Subpatt(YES,tree); */
      /*     move->lnL = Triple_Dist(n_link,tree,YES); */
      /*     MIXT_Set_Alias_Subpatt(NO,tree); */
      /*   } */
    }
  else
    {
      Update_P_Pars(tree,b_residual,n_link);
      Pars(b_residual,tree);
    }

  n1 = (b_residual->left == n_link)?(b_residual->rght):(b_residual->left);
  n2 = (b_target->left   == n_link)?(b_target->rght):(b_target->left);
  dir_v1 = dir_v2 = dir_v0 = -1;
  For(i,3)
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

  For(i,tree->depth_curr_path+1) move->path[i] = tree->curr_path[i];

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
  
  Include_One_Spr_To_List_Of_Spr(move,tree);

  Copy_Scalar_Dbl(init_target_l,b_target->l);
  Copy_Scalar_Dbl(init_target_v,b_target->l_var);

  Copy_Scalar_Dbl(init_arrow_l,b_arrow->l);
  Copy_Scalar_Dbl(init_arrow_v,b_arrow->l_var);

  Copy_Scalar_Dbl(init_residual_l,b_residual->l);
  Copy_Scalar_Dbl(init_residual_v,b_residual->l_var);

  Prune_Subtree(n_link,
                (n_link==b_arrow->left)?(b_arrow->rght):(b_arrow->left),
                &b_target,
                &b_residual,
                tree);

  if(tree->mod->s_opt->spr_lnL == YES) Update_PMat_At_Given_Edge(b_target,tree);

  tree->c_lnL   = init_lnL;
  tree->c_pars  = init_pars;


  Free_Scalar_Dbl(init_target_l);
  Free_Scalar_Dbl(init_arrow_l);
  Free_Scalar_Dbl(init_residual_l);
  Free_Scalar_Dbl(init_target_v);
  Free_Scalar_Dbl(init_arrow_v);
  Free_Scalar_Dbl(init_residual_v);

  return move_lnL;
}

/*********************************************************/

void Speed_Spr_Loop(t_tree *tree)
{
  phydbl lk_old,delta_lnL;
  int i;


  Spr_List_Of_Trees(tree);
  return;

  tree->best_pars                  = 1E+8;
  tree->mod->s_opt->spr_lnL        = NO;
  tree->mod->s_opt->spr_pars       = NO;
  tree->mod->s_opt->quickdirty     = NO;

  if((tree->mod->s_opt->print) && (!tree->io->quiet)) PhyML_Printf("\n\n. Maximizing likelihood (using SPR moves)...\n");

  tree->mod->s_opt->max_depth_path = tree->n_otu;
  Spr_Pars(0,10,tree);
  Set_Both_Sides(NO,tree);
  Lk(NULL,tree);

  Round_Optimize(tree,tree->data,1);

  tree->best_pars = tree->c_pars;
  tree->best_lnL  = tree->c_lnL;

  /*****************************/
  if(tree->mod->s_opt->print == YES && tree->io->quiet == NO) PhyML_Printf("\n\n. First round of SPR moves...\n");
  lk_old = tree->c_lnL;
  tree->mod->s_opt->max_depth_path    = tree->n_otu;
  tree->mod->s_opt->max_delta_lnL_spr = (tree->io->datatype == NT)?(-1.):(-1.);
  tree->mod->s_opt->spr_lnL           = NO;
  tree->mod->s_opt->spr_pars          = NO;
  tree->mod->s_opt->min_diff_lk_move  = 0.1;
  delta_lnL                           = 5.0;
  Speed_Spr(tree,0.5,tree->n_otu,delta_lnL);
  Optimiz_All_Free_Param(tree,(tree->io->quiet)?(0):(tree->mod->s_opt->print));

  /*****************************/
  if(tree->mod->s_opt->print == YES && tree->io->quiet == NO) PhyML_Printf("\n\n. Second round of SPR moves...\n");
  lk_old = tree->c_lnL;
  tree->mod->s_opt->max_depth_path    = tree->n_otu;
  tree->mod->s_opt->max_delta_lnL_spr = (tree->io->datatype == NT)?(2.):(0.);
  tree->mod->s_opt->spr_lnL           = YES;
  tree->mod->s_opt->spr_pars          = NO;
  tree->mod->s_opt->min_diff_lk_move  = 0.01;
  delta_lnL                           = 1.0;
  Speed_Spr(tree,1.0,20,delta_lnL);
  Optimiz_All_Free_Param(tree,(tree->io->quiet)?(0):(tree->mod->s_opt->print));
  /*****************************/

  /*****************************/
  if(tree->mod->s_opt->print == YES && tree->io->quiet == NO) PhyML_Printf("\n\n. NNI moves...\n");
  tree->mod->s_opt->min_diff_lk_move  = 0.001;
  lk_old = UNLIKELY;
  do
    {
      lk_old = tree->c_lnL;
      if(!Simu(tree,5)) break;
    }
  while(FABS(lk_old - tree->c_lnL) > tree->mod->s_opt->min_diff_lk_local);
  /*****************************/

  For(i,2*tree->n_otu-3) if(tree->a_edges[i]->l->v < 1.E-3) tree->a_edges[i]->l->v = 1.E-3;
  Round_Optimize(tree,tree->data,ROUND_MAX);

  /* /\*****************************\/ */
  /* do */
  /*   { */
  /*     Round_Optimize(tree,tree->data,ROUND_MAX); */
  /*     if(!Check_NNI_Five_Branches(tree)) break; */
  /*   } */
  /* while(1); */
  /* /\*****************************\/ */

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Speed_Spr(t_tree *tree, phydbl prop_spr, int max_cycles, phydbl delta_lnL)
{
  int step,old_pars;
  phydbl old_lnL;

  if(tree->lock_topo == YES)
    {
      PhyML_Printf("\n== The tree topology is locked.");
      PhyML_Printf("\n== Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
      Exit("\n");
    }


  Set_Both_Sides(YES,tree);
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

      Spr(UNLIKELY,prop_spr,tree);

      // Set maximum depth for future spr rounds to deepest spr found so far
      tree->mod->s_opt->max_depth_path = tree->max_spr_depth;

      if(tree->mod->s_opt->spr_pars == NO)
        {
          if(tree->n_improvements > 0)
            {
              /* Optimise branch lengths */
              Optimize_Br_Len_Serie(tree);
              /* Update partial likelihoods */
              Set_Both_Sides(YES,tree);
              Lk(NULL,tree);
              /* Print log-likelihood and parsimony scores */
              if((tree->mod->s_opt->print) && (!tree->io->quiet)) Print_Lk(tree,"[Branch lengths     ]");
            }
        }
      else
        {
          if((tree->mod->s_opt->print) && (!tree->io->quiet)) Print_Pars(tree);
        }

      Pars(NULL,tree);

      if(tree->io->print_trace)
        {
          char *s = Write_Tree(tree,NO);
          PhyML_Fprintf(tree->io->fp_out_trace,"[%f]%s\n",tree->c_lnL,s); fflush(tree->io->fp_out_trace);
          if((tree->io->print_site_lnl) && (!tree->mod->s_opt->spr_pars)) Print_Site_Lk(tree,tree->io->fp_out_lk); fflush(tree->io->fp_out_lk);
          Free(s);
        }

      if(tree->io->print_json_trace == YES) JSON_Tree_Io(tree,tree->io->fp_out_json_trace); 
      
      /* Record the current best log-likelihood and parsimony */
      tree->best_lnL  = tree->c_lnL;
      tree->best_pars = tree->c_pars;

      if(tree->mod->s_opt->spr_pars == NO)
        {
          if(tree->c_lnL < old_lnL-tree->mod->s_opt->min_diff_lk_local)
            {
              PhyML_Printf("\n== old_lnL = %f c_lnL = %f",old_lnL,tree->c_lnL);
              PhyML_Printf("\n== Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
              Exit("");
            }
        }
      else
        {
          if(tree->c_pars > old_pars)
            {
              PhyML_Printf("\n== old_pars = %d c_pars = %d",old_pars,tree->c_pars);
              PhyML_Printf("\n== Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
              Exit("");
            }
        }

      /* Record the current best branch lengths  */
      Record_Br_Len(tree);

      /* Exit if no improvements after complete optimization */
      if(step+1 > max_cycles) break;
      if((tree->mod->s_opt->spr_pars == NO)  && (FABS(old_lnL-tree->c_lnL)   < delta_lnL)) break;
      if((tree->mod->s_opt->spr_pars == YES) && (FABS(old_pars-tree->c_pars) < 1)) break;
      if(!tree->n_improvements) break;
    }
  while(1);
}

/*********************************************************/

int Evaluate_List_Of_Regraft_Pos_Triple(t_spr **spr_list, int list_size, t_tree *tree)
{
  t_spr *move;
  t_edge *init_target, *b_residual;
  int i,j,best_move;
  int dir_v0, dir_v1, dir_v2;
  scalar_dbl *recorded_l,*recorded_v;
  phydbl best_lnL,init_lnL;
  int recorded;

  if(tree->mixt_tree != NULL)
    {
      PhyML_Printf("\n== Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
      Exit("\n");
    }

  best_lnL = UNLIKELY;
  init_target = b_residual = NULL;
  best_move = -1;
  init_lnL = tree->c_lnL;
  recorded_v = recorded_l = NULL;

  if(!list_size && !tree->io->fp_in_constraint_tree)
    {
      PhyML_Printf("\n== List size is 0 !");
      PhyML_Printf("\n== Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
      Exit("\n");
    }

  recorded = NO;
  For(i,list_size)
    {
      move = spr_list[i];

      if(!move)
        {
          PhyML_Printf("\n== move is NULL\n");
          PhyML_Printf("\n== Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
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
              
              Fast_Br_Len(init_target,tree,NO);

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
          Update_P_Lk_Along_A_Path(move->path,move->depth_path+1,tree);
          MIXT_Set_Alias_Subpatt(NO,tree);
          
          /* Regraft subtree */
          Graft_Subtree(move->b_target,move->n_link,b_residual,tree);
                    
          MIXT_Set_Alias_Subpatt(YES,tree);
          move->lnL = Triple_Dist(move->n_link,tree,YES);

          MIXT_Set_Alias_Subpatt(NO,tree);

          if((move->lnL < best_lnL) && (move->lnL > best_lnL - tree->mod->s_opt->max_delta_lnL_spr))
            {
              /* Estimate the three t_edge lengths at the regraft site */
              MIXT_Set_Alias_Subpatt(YES,tree);
              move->lnL = Triple_Dist(move->n_link,tree,NO);
              MIXT_Set_Alias_Subpatt(NO,tree);
            }

          /* printf("\n. %d/%d move->lnL= %f best_lnL=%f absolute_best=%f",i,list_size,move->lnL,best_lnL,tree->best_lnL); */
          
          /* Record updated branch lengths for this move */
          dir_v1 = dir_v2 = dir_v0 = -1;
          For(j,3)
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
          
          if(move->lnL > best_lnL + tree->mod->s_opt->min_diff_lk_move)
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
                        b_residual,
                        tree);
          
          /* Restore branch lengths */
          Restore_Br_Len(tree);
          
          /* Update relevant change proba matrices */
          /* 	  Update_PMat_At_Given_Edge(move->n_link->b[0],tree); */
          /* 	  Update_PMat_At_Given_Edge(move->n_link->b[1],tree); */
          /* 	  Update_PMat_At_Given_Edge(move->n_link->b[2],tree); */
          Update_PMat_At_Given_Edge(move->b_target,tree);
          
          /* 	  Update_P_Lk(tree,move->n_link->b[0],move->n_link); */
          /* 	  Update_P_Lk(tree,move->n_link->b[1],move->n_link); */
          /* 	  Update_P_Lk(tree,move->n_link->b[2],move->n_link); */
          
          /* Update conditional likelihoods along the path from the prune to
             the regraft position */
          /* 	  Update_P_Lk_Along_A_Path(move->path,move->depth_path+1,tree); */
          
          tree->c_lnL = init_lnL;
        }
      

      /* Bail out as soon as you've found a true improvement */
      /* if(move->lnL > tree->best_lnL + 1.0) break; */
      if(move->lnL > tree->best_lnL + tree->mod->s_opt->min_diff_lk_move) break;
      
    }
  
  /* PhyML_Printf("\n. [ %4d/%4d ] %f",i,list_size,tree->best_lnL); */
  /*   PhyML_Printf("\n. max_improv = %f",max_improv); */
  
  
  MIXT_Set_Alias_Subpatt(YES,tree);
  For(i,list_size)
    {
      move = spr_list[i];
      if(move->b_target)
        {
          For(j,3) Update_PMat_At_Given_Edge(move->n_link->b[j],tree);
          For(j,3) Update_P_Lk(tree,move->n_link->b[j],move->n_link);
          
          /* TO DO : we don't need to update all these partial likelihoods here.
             Would need to record only those that were along the paths examined
             above */
          
          For(j,3)
            if(move->n_link->v[j] != move->n_opp_to_link)
              Pre_Order_Lk(move->n_link,move->n_link->v[j],tree);
          
          break;
        }
    }  
  MIXT_Set_Alias_Subpatt(NO,tree);
  
#ifdef DEBUG
  if(best_move < 0 && list_size > 0)
    {
      PhyML_Printf("\n\n. Best_move < 0 !");
      
      PhyML_Printf("\n. List size = %d",list_size);
      For(i,list_size)
        {
          move = spr_list[i];
          PhyML_Printf("\n. %p %p",move,move->b_target);
        }
      
      PhyML_Printf("\n== Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
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

  if(tree->mixt_tree != NULL)
    {
      PhyML_Printf("\n== Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
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

  Graft_Subtree(move->b_target,move->n_link,b_residual,tree);

  dir_v1 = dir_v2 = dir_v0 = -1;
  For(j,3)
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
      Set_Both_Sides(YES,tree);
      MIXT_Set_Alias_Subpatt(YES,tree);
      Lk(NULL,tree);
      MIXT_Set_Alias_Subpatt(NO,tree);

      if(FABS(tree->c_lnL - move->lnL) > tree->mod->s_opt->min_diff_lk_move)
        {
          PhyML_Printf("\n== c_lnL = %f move_lnL = %f", tree->c_lnL,move->lnL);
          PhyML_Printf("\n== %d l0=%G l1=%G l2=%G v0=%G v1=%G v2=%G",move->n_link->num,move->l0->v,move->l1->v,move->l2->v,move->v0->v,move->v1->v,move->v2->v);
          PhyML_Printf("\n== Gamma MGF? %d",tree->io->mod->gamma_mgf_bl);
          PhyML_Printf("\n== Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
          Check_Lk_At_Given_Edge(YES,tree);
          /* { */
          /*   int i; */
          /*   For(i,2*tree->n_otu-1) */
          /*     printf("\n. Edge %d %G %G", */
          /*            i, */
          /*            tree->a_edges[i]->l->v, */
          /*            tree->a_edges[i]->l_var->v); */
          /* } */
          Exit("\n");
        }

      if((tree->mod->s_opt->print) && (!tree->io->quiet))
        {
          Print_Lk_And_Pars(tree);
          PhyML_Printf(" [depth=%5d]",move->depth_path); fflush(NULL);
        }

      if(move->depth_path > tree->max_spr_depth) tree->max_spr_depth = move->depth_path;

      tree->n_improvements++;
      if(tree->c_lnL > tree->best_lnL) tree->best_lnL = tree->c_lnL;
      Record_Br_Len(tree);

      if(move->depth_path > tree->mod->s_opt->deepest_path)
        tree->mod->s_opt->deepest_path = move->depth_path;

      return 1;
    }

  Prune_Subtree(move->n_link,
                move->n_opp_to_link,
                &move->b_target,
                &b_residual,
                tree);

  Graft_Subtree(init_target,
                move->n_link,
                b_residual,
                tree);
  
  Restore_Br_Len(tree);
  
  Set_Both_Sides(YES,tree);
  MIXT_Set_Alias_Subpatt(YES,tree);
  Lk(NULL,tree);
  MIXT_Set_Alias_Subpatt(NO,tree);
  Pars(NULL,tree);
  return 0;
}

/*********************************************************/

int Try_One_Spr_Move_Full(t_spr *move, t_tree *tree)
{
  t_edge *init_target, *b_residual;

  Record_Br_Len(tree);

  Prune_Subtree(move->n_link,
        move->n_opp_to_link,
        &init_target,
        &b_residual,
        tree);

  Graft_Subtree(move->b_target,move->n_link,b_residual,tree);

  MIXT_Set_Alias_Subpatt(YES,tree);
  Set_Both_Sides(YES,tree);
  Lk(NULL,tree);
  MIXT_Set_Alias_Subpatt(NO,tree);
  Optimize_Br_Len_Serie(tree);

  Set_Both_Sides(YES,tree);
  Lk(NULL,tree);

  if(tree->c_lnL > tree->best_lnL + tree->mod->s_opt->min_diff_lk_move)
    {
      Pars(NULL,tree);
      if((tree->mod->s_opt->print) && (!tree->io->quiet)) Print_Lk(tree,"[Topology           ]");
      tree->n_improvements++;
      tree->best_lnL = tree->c_lnL;
      Record_Br_Len(tree);
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
            b_residual,
            tree);

      Restore_Br_Len(tree);
      Set_Both_Sides(YES,tree);

      MIXT_Set_Alias_Subpatt(YES,tree);
      Lk(NULL,tree);
      MIXT_Set_Alias_Subpatt(NO,tree);
      Pars(NULL,tree);
      return 0;
    }

  return -1;
}

/*********************************************************/

void Include_One_Spr_To_List_Of_Spr(t_spr *move, t_tree *tree)
{
  int i;
  t_spr *buff_spr,*orig_move, *orig_move_list, *move_list;
  t_tree *orig_tree;

  if(tree->mixt_tree != NULL)
    {
      PhyML_Printf("\n== Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
      Exit("\n");
    }

  if(((tree->mod->s_opt->spr_lnL == YES) && (move->lnL  > tree->spr_list[tree->size_spr_list-1]->lnL)) ||
     ((tree->mod->s_opt->spr_lnL == NO) && (move->pars <= tree->spr_list[tree->size_spr_list-1]->pars)))
    {
      move_list = tree->spr_list[tree->size_spr_list-1];

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

          if(tree->next)
            {
              move      = move->next;
              move_list = move_list->next;
              tree      = tree->next;
            }
          else
            {
              move      = move->next;
              move_list = move_list->next;
              tree      = tree->next;
            }
        }
      while(tree);
      move      = orig_move;
      move_list = orig_move_list;
      tree      = orig_tree;
      
      For(i,move_list->depth_path+1) move_list->path[i] = move->path[i];
      
      for(i=tree->size_spr_list-1;i>0;i--)
        {
          if(((tree->mod->s_opt->spr_lnL == YES) && (tree->spr_list[i]->lnL > tree->spr_list[i-1]->lnL)) ||
             ((tree->mod->s_opt->spr_lnL == NO) && (tree->spr_list[i]->pars <=  tree->spr_list[i-1]->pars)))
            {
              
              orig_tree = tree;
              do
                {
                  buff_spr            = tree->spr_list[i-1];
                  tree->spr_list[i-1] = tree->spr_list[i];
                  tree->spr_list[i]   = buff_spr;
                  
                  if(tree->next) tree = tree->next;
                  else           tree = tree->next;
                }
              while(tree);
              tree = orig_tree;
              
            }
          else  break;
        }
    }
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

  For(i,n_moves)
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
                    residual,tree);
    }
  Free(spr_struct);
}

/*********************************************************/

void Reset_Spr_List(t_tree *tree)
{
  int i;

  For(i,tree->size_spr_list)
    {
      tree->spr_list[i]->depth_path     = 0;
      tree->spr_list[i]->pars           = MAX_PARS;
      tree->spr_list[i]->lnL            = UNLIKELY;
      tree->spr_list[i]->n_link         = NULL;
      tree->spr_list[i]->n_opp_to_link  = NULL;
      tree->spr_list[i]->b_target       = NULL;
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
  int curr_pars;

  PhyML_Printf("\n. Minimizing parsimony...\n");

  tree->best_pars                  = 1E+8;
  tree->best_lnL                   = UNLIKELY;
  tree->mod->s_opt->spr_lnL        = NO;
  tree->mod->s_opt->spr_pars       = YES;
  curr_pars                        = tree->c_pars;
  tree->mod->s_opt->max_depth_path = tree->n_otu;
  do
    {
      curr_pars = tree->c_pars;
      Speed_Spr(tree,1.0,1,0.0);
    }
  while(tree->n_improvements && FABS(tree->c_pars - curr_pars) > threshold);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Spr_Shuffle(t_tree *mixt_tree)
{
  phydbl lk_old;
  int *orig_catg,n,n_iter;
  t_tree *tree,**tree_list;

  if(mixt_tree->mod->s_opt->print) PhyML_Printf("\n\n. Refining the tree...\n");

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
  mixt_tree->mod->s_opt->max_delta_lnL_spr = 0.;
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

      if(mixt_tree->n_improvements < 5      || 
         mixt_tree->max_spr_depth  < 2      || 
         FABS(lk_old-mixt_tree->c_lnL) < 5. ||
         ++n_iter > 10) break;
      
    }
  while(1);

  mixt_tree->annealing_temp = 0.0;

  if(mixt_tree->mod->s_opt->print && (!mixt_tree->io->quiet))
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
  int step,i,n_targets,n_rand,no_improvement;
  t_tree *best_tree;
  scalar_dbl **best_bl;
  t_node *rnd_node;
  t_edge *b_target,*b_residual,**target_list,*rnd_edge;
  phydbl true_best_lnL;

  if(tree->lock_topo == YES)
    {
      PhyML_Printf("\n== The tree topology is locked.");
      PhyML_Printf("\n== Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
      Exit("\n");
    }

  Set_Both_Sides(NO,tree);
  Pars(NULL,tree);
  Lk(NULL,tree);

  tree->mod->s_opt->max_depth_path    = (int)(tree->n_otu/3);
  tree->mod->s_opt->max_delta_lnL_spr = (tree->io->datatype == NT)?(5.):(0.);
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
  no_improvement                      = 0;
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
              For(i,3)
                if(b_target->left->v[i] != b_target->rght)
                  Get_List_Of_Adjacent_Targets(b_target->left,b_target->left->v[i],NULL,&target_list,&n_targets,0,5);
              
              For(i,3)
                if(b_target->rght->v[i] != b_target->left)
                  Get_List_Of_Adjacent_Targets(b_target->rght,b_target->rght->v[i],NULL,&target_list,&n_targets,0,5);
              
              if(n_targets > 0) b_target = target_list[Rand_Int(0,n_targets-1)];
              
              assert(b_target != NULL);
              
              Graft_Subtree(b_target,rnd_node,b_residual,tree);
              
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

      Set_Both_Sides(YES,tree);
      Optimiz_All_Free_Param(tree,(tree->io->quiet)?(0):(tree->mod->s_opt->print));
      Optimize_Br_Len_Serie(tree);

      if(tree->io->print_trace)
        {
          char *s = Write_Tree(tree,NO);
          PhyML_Fprintf(tree->io->fp_out_trace,"[%f]%s\n",tree->c_lnL,s); fflush(tree->io->fp_out_trace);
          if((tree->io->print_site_lnl) && (!tree->mod->s_opt->spr_pars)) Print_Site_Lk(tree,tree->io->fp_out_lk); fflush(tree->io->fp_out_lk);
          Free(s);
        }

      if(tree->io->print_json_trace == YES) JSON_Tree_Io(tree,tree->io->fp_out_json_trace); 
      
      /* Record the current best log-likelihood and parsimony */
      if(tree->c_lnL > true_best_lnL)
        {
          no_improvement = 0;
          true_best_lnL = tree->c_lnL;
          For(i,2*tree->n_otu-1) Free_Scalar_Dbl(best_bl[i]);
          Free(best_bl);
          best_bl = Copy_Br_Len(tree);
          Copy_Tree(tree,best_tree); /* Record tree topology, branch lengths and model parameters */
        }
      else
        {
          no_improvement++;
        }
      
      Transfer_Br_Len_To_Tree(best_bl,tree);
      Copy_Tree(best_tree,tree);

    }
  while(++step <= max_cycles && tree->n_improvements > 0 && tree->max_spr_depth  > 1);

  Free(target_list);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Spr_List_Of_Trees(t_tree *tree)
{
  int i,j,list_size,max_list_size,*rk;
  t_tree **tree_list;
  phydbl delta_lnL,*lnL_list,anneal_temp,best_lnL,worst_lnL;
  
  max_list_size                    = tree->n_otu;
  tree->mod->s_opt->max_depth_path = tree->n_otu;
  anneal_temp                      = 0.0;

  tree_list = (t_tree **)mCalloc(max_list_size,sizeof(t_tree *));
  lnL_list  = (phydbl *)mCalloc(max_list_size,sizeof(phydbl));

  For(i,max_list_size) lnL_list[i] = UNLIKELY;

  Spr_Pars(0,100,tree);
  Round_Optimize(tree,tree->data,5);
  tree_list[0] = Make_Tree_From_Scratch(tree->n_otu,tree->data);
  Copy_Tree(tree,tree_list[0]);
  lnL_list[0] = tree->c_lnL;


  if(tree->mod->s_opt->print == YES && tree->io->quiet == NO) PhyML_Printf("\n\n. First round of SPR moves...\n");
  list_size = 1;
  best_lnL = UNLIKELY;
  worst_lnL = -UNLIKELY;
  do
    {
      printf("\n");
      Randomize_Tree(tree,1);

      Set_Both_Sides(NO,tree);
      Lk(NULL,tree);
      Optimize_Br_Len_Serie (tree);

      Set_Both_Sides(YES,tree);
      Lk(NULL,tree);

      tree->annealing_temp                = anneal_temp;
      tree->mod->s_opt->max_depth_path    = 10;
      tree->mod->s_opt->max_delta_lnL_spr = (tree->io->datatype == NT)?(-1.0):(-1.);
      tree->mod->s_opt->spr_lnL           = NO;
      tree->mod->s_opt->spr_pars          = NO;
      tree->mod->s_opt->min_diff_lk_move  = 0.1;
      tree->best_lnL                      = tree->c_lnL;

      Spr(UNLIKELY,0.5,tree);
      Set_Both_Sides(NO,tree);
      Lk(NULL,tree);
      Optimize_Br_Len_Serie (tree);

      printf("\n>> lnL: %f",tree->c_lnL);

      tree_list[list_size] = Make_Tree_From_Scratch(tree->n_otu,tree->data);
      Copy_Tree(tree,tree_list[list_size]);
      lnL_list[list_size] = tree->c_lnL;
      if(tree->c_lnL > best_lnL) best_lnL   = tree->c_lnL;
      if(tree->c_lnL < worst_lnL) worst_lnL = tree->c_lnL;
    }
  while(list_size++ < 1 + (int)tree->n_otu/10);

  delta_lnL = (best_lnL - worst_lnL) / 4.;
  
  rk = Ranks(lnL_list,max_list_size);

  printf("\n. delta_lnL = %f",delta_lnL);

  do
    {
      best_lnL = UNLIKELY;
      worst_lnL = -UNLIKELY;
      For(i,list_size)
        {
          printf("\n"); printf("\n. current delta: %f [%f]",best_lnL-worst_lnL,delta_lnL);
          Copy_Tree(tree_list[rk[i]],tree);
          Set_Both_Sides(YES,tree);
          Lk(NULL,tree);
          tree->annealing_temp                = anneal_temp;
          tree->mod->s_opt->max_depth_path    = 10;
          tree->mod->s_opt->max_delta_lnL_spr = (tree->io->datatype == NT)?(-1.0):(-1.);
          tree->mod->s_opt->spr_lnL           = YES;
          tree->mod->s_opt->spr_pars          = NO;
          tree->mod->s_opt->min_diff_lk_move  = 0.01;
          tree->best_lnL                      = tree->c_lnL;
          Spr(UNLIKELY,1.0,tree);
          For(j,2*tree->n_otu-3) if(tree->a_edges[j]->l->v < 1.E-3) tree->a_edges[j]->l->v = 1.E-3;
          Set_Both_Sides(NO,tree);
          Lk(NULL,tree);
          Optimize_Br_Len_Serie (tree);
          /* tree->mod->s_opt->fast_nni = YES; */
          /* Simu(tree,10); */
          printf("\n>> lnL: %f",tree->c_lnL);
          Copy_Tree(tree,tree_list[rk[i]]);
          lnL_list[rk[i]] = tree->c_lnL;          
          if(tree->c_lnL > best_lnL) best_lnL   = tree->c_lnL;
          if(tree->c_lnL < worst_lnL) worst_lnL = tree->c_lnL;
          if(tree->c_lnL < best_lnL - delta_lnL) { list_size = i; break; }
        }

      delta_lnL = (best_lnL - worst_lnL) / 2.;
      anneal_temp *= 0.5;

      rk = Ranks(lnL_list,max_list_size);
      delta_lnL *= 0.7;
      printf("\n. delta_lnL: %f best_lnL: %f list_size: %d anneal: %f",delta_lnL,best_lnL,list_size,anneal_temp); fflush(NULL);

      Copy_Tree(tree_list[rk[0]],tree);
      Round_Optimize(tree,tree->data,5);

    }
  while(delta_lnL > 1.0);

  Copy_Tree(tree_list[rk[0]],tree);


  /* if(tree->mod->s_opt->print == YES && tree->io->quiet == NO) PhyML_Printf("\n\n. First round of SPR moves...\n"); */
  /* list_size = 0; */
  /* do */
  /*   { */
  /*     printf("\n"); */
  /*     tree->annealing_temp                = anneal_temp; */
  /*     tree->mod->s_opt->max_depth_path    = 10; */
  /*     tree->mod->s_opt->max_delta_lnL_spr = (tree->io->datatype == NT)?(-1.):(-1.); */
  /*     tree->mod->s_opt->spr_lnL           = NO; */
  /*     tree->mod->s_opt->spr_pars          = NO; */
  /*     tree->mod->s_opt->min_diff_lk_move  = 0.1; */
  /*     Set_Both_Sides(YES,tree); */
  /*     Lk(NULL,tree); */
  /*     if(!Spr(UNLIKELY,0.5,tree)) break; */
  /*     Optimize_Br_Len_Serie (tree); */
  /*     tree_list[list_size] = Make_Tree_From_Scratch(tree->n_otu,tree->data); */
  /*     Copy_Tree(tree,tree_list[list_size]); */
  /*     lnL_list[list_size] = tree->c_lnL; */
  /*     anneal_temp -= 0.2; */
  /*     if(anneal_temp < 0.0) anneal_temp = 0.0; */
  /*     list_size++; */
  /*   } */
  /* while(1); */



  /* rk = Ranks(lnL_list,max_list_size); */

  /* if(tree->mod->s_opt->print == YES && tree->io->quiet == NO) PhyML_Printf("\n\n. Second round of SPR moves...\n"); */
  /* list_size = (int)list_size/2; */
  /* For(i,list_size) */
  /*   { */
  /*     printf("\n"); */
  /*     printf("\n. tree %p -> %f",tree_list[rk[i]],lnL_list[rk[i]]); */
  /*     Copy_Tree(tree_list[rk[i]],tree); */
  /*     tree->annealing_temp                = 0.; */
  /*     tree->mod->s_opt->max_depth_path    = tree->n_otu; */
  /*     tree->mod->s_opt->max_delta_lnL_spr = (tree->io->datatype == NT)?(-1.):(-1.); */
  /*     tree->mod->s_opt->spr_lnL           = NO; */
  /*     tree->mod->s_opt->spr_pars          = NO; */
  /*     tree->mod->s_opt->min_diff_lk_move  = 0.1; */
  /*     delta_lnL                           = 2.0; */
  /*     Speed_Spr(tree,1.0,tree->n_otu,delta_lnL); */
  /*     Copy_Tree(tree,tree_list[rk[i]]); */
  /*     lnL_list[rk[i]] = tree->c_lnL; */
  /*   } */


  /* rk = Ranks(lnL_list,max_list_size); */

  /* if(tree->mod->s_opt->print == YES && tree->io->quiet == NO) PhyML_Printf("\n\n. Third round of SPR moves...\n"); */
  /* list_size = 1; */
  /* For(i,list_size) */
  /*   { */
  /*     Copy_Tree(tree_list[rk[i]],tree); */
  /*     printf("\n. tree %p -> %f",tree_list[rk[i]],lnL_list[rk[i]]); */
  /*     tree->mod->s_opt->max_depth_path    = 5; */
  /*     tree->mod->s_opt->max_delta_lnL_spr = (tree->io->datatype == NT)?(2.):(0.); */
  /*     tree->mod->s_opt->spr_lnL           = YES; */
  /*     tree->mod->s_opt->spr_pars          = NO; */
  /*     tree->mod->s_opt->min_diff_lk_move  = 0.01; */
  /*     delta_lnL                           = 1.0; */
  /*     Speed_Spr(tree,1.0,20,delta_lnL); */
  /*     Copy_Tree(tree,tree_list[rk[i]]); */
  /*     lnL_list[rk[i]] = tree->c_lnL; */
  /*   } */


  /* rk = Ranks(lnL_list,max_list_size); */

  /* list_size = 1; */
  /* For(i,list_size) */
  /*   {       */
  /*     Copy_Tree(tree_list[rk[i]],tree); */
  /*     if(tree->mod->s_opt->print == YES && tree->io->quiet == NO) PhyML_Printf("\n\n. Last optimization step...\n"); */
  /*     For(i,2*tree->n_otu-3) if(tree->a_edges[i]->l->v < 1.E-3) tree->a_edges[i]->l->v = 1.E-3; */
  /*     Round_Optimize(tree,tree->data,ROUND_MAX); */
  /*     Copy_Tree(tree,tree_list[rk[i]]); */
  /*     lnL_list[rk[i]] = tree->c_lnL; */
  /*   } */

  Free(tree_list);
  Free(lnL_list);
  Free(rk);
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

/*
** EOF: spr.c
*/
