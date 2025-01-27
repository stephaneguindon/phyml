/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include "cv.h"
// Routine for cross validation. 
phydbl *CV_Tip_Cv(t_tree *tree)
{
  t_edge    *b;
  int        obs_d_state, ambiguity_check;
  char       obs_state;
  int        ns, nsncatg, n_otu;
  phydbl    *p_lk_left, *Pij, *state_probs_all_sites, *weights;
  phydbl     sum;
  short int *truth;

  if (tree->is_mixt_tree == YES)
  {
    MIXT_Cv(tree);
    return (NULL);
  }

  Set_Both_Sides(YES, tree);
  Lk(NULL, tree);

  b         = NULL;
  p_lk_left = NULL;

  state_probs_all_sites = (phydbl *)mCalloc(
      tree->data->n_pattern * tree->n_otu * tree->mod->ns, sizeof(phydbl));
  truth   = (short int *)mCalloc(tree->data->n_pattern * tree->n_otu * tree->mod->ns,
                                 sizeof(phydbl));
  weights = (phydbl *)mCalloc(tree->data->n_pattern * tree->n_otu, sizeof(phydbl));

  ns      = tree->mod->ns;
  nsncatg = ns * tree->mod->ras->n_catg;
  n_otu   = tree->n_otu;

  for (int tax_id = 0; tax_id < tree->n_otu; ++tax_id)
  {
    b = tree->a_nodes[tax_id]->b[0];

    p_lk_left = b->p_lk_left;
    Pij       = b->Pij_rr;

    for (int site = 0; site < tree->data->n_pattern; ++site)
    {
      if (tree->data->wght[site] > SMALL)
      {
        ambiguity_check = tree->a_nodes[tax_id]->c_seq->is_ambigu[site];
        if (ambiguity_check == NO)
        {
          obs_state   = tree->a_nodes[tax_id]->c_seq->state[site];
          obs_d_state = tree->a_nodes[tax_id]->c_seq->d_state[site];

          tree->a_nodes[tax_id]->c_seq->is_ambigu[site] = YES;
          tree->a_nodes[tax_id]->c_seq->state[site]     = '?';
          tree->a_nodes[tax_id]->c_seq->d_state[site]   = -1;

          Init_Partial_Lk_Tips_Double_One_Character(tax_id, site, tree);

          Br_Len_Opt(&(b->l->v), b, tree);
          // Optimize_Br_Len_Serie(1000,tree);
          // Optimiz_All_Free_Param(tree,NO);

          for (int tip_state = 0; tip_state < ns; ++tip_state)
          {
            state_probs_all_sites[site * ns * n_otu + tax_id * ns + tip_state] =
                0.0;
            for (int int_state = 0; int_state < ns; ++int_state)
            {
              for (int catg = 0; catg < tree->mod->ras->n_catg; ++catg)
              {
                state_probs_all_sites[site * ns * n_otu + tax_id * ns +
                                      tip_state] +=
                    tree->mod->ras->gamma_r_proba->v[catg] *
                    tree->mod->e_frq->pi->v[tip_state] *
                    Pij[catg * ns * ns + tip_state * ns + int_state] *
                    p_lk_left[catg * ns + int_state];
              }
            }
          }

          sum = 0.0;
          for (int state = 0; state < ns; ++state)
            sum +=
                state_probs_all_sites[site * ns * n_otu + tax_id * ns + state];

          for (int state = 0; state < ns; ++state)
            state_probs_all_sites[site * ns * n_otu + tax_id * ns + state] /=
                sum;

          for (int state = 0; state < ns; ++state)
            truth[site * ns * n_otu + tax_id * ns + state] = 0;

          truth[site * ns * n_otu + tax_id * ns + obs_d_state] = 1;

          PhyML_Printf("\n###%s,%s,%d,%g", tree->mod->modelname->s,
                       tree->a_nodes[tax_id]->name, site,
                       log(state_probs_all_sites[site * ns * n_otu +
                                                 tax_id * ns + obs_d_state]));

          tree->a_nodes[tax_id]->c_seq->state[site]     = obs_state;
          tree->a_nodes[tax_id]->c_seq->d_state[site]   = obs_d_state;
          tree->a_nodes[tax_id]->c_seq->is_ambigu[site] = NO;

          Init_Partial_Lk_Tips_Double_One_Character(tax_id, site, tree);
        }
      }
      p_lk_left += nsncatg;
    }
  }
  
  for (int tax_id = 0; tax_id < tree->n_otu; ++tax_id)
    for (int site = 0; site < tree->data->n_pattern; ++site)
      weights[site * tree->n_otu + tax_id] = tree->data->wght[site];

  // for (int tax_id = 0; tax_id < tree->n_otu; ++tax_id)
  //   for (int site = 0; site < tree->n_pattern; ++site)
  //     for (int state = 0; state < tree->mod->ns; ++state)
  //     {
  //       PhyML_Printf("\n. CV tax: %s site: %d state: %d prob: %f truth: %d weight: %f",
  //                    tree->a_nodes[tax_id]->name, site, state,
  //                    state_probs_all_sites[site * tree->n_otu * tree->mod->ns +
  //                                          tax_id * tree->mod->ns + state],
  //                    truth[site * tree->n_otu * tree->mod->ns +
  //                          tax_id * tree->mod->ns + state],
  //                    weights[site * tree->n_otu + tax_id]);
  //     }

  ROC(state_probs_all_sites, truth, ns, tree->n_otu * tree->data->n_pattern, weights,
      tree->mod->modelname->s);

  Free(weights);
  
  return (state_probs_all_sites);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void CV_Hide_Characters_At_Random(calign *data, phydbl mask_prob)
{
  phydbl u;

  for (int site = 0; site < data->n_pattern; ++site)
  {
    for (int tax_id = 0; tax_id < data->n_otu; ++tax_id)
    {
      u = Uni();

      if (!(u > mask_prob))
      {

        // PhyML_Printf("\n Hidding '%c' for taxon %s at site %d",
        //              data->c_seq[tax_id]->state[site],
        //              data->c_seq[tax_id]->name, site);

        data->c_seq[tax_id]->state[site]     = '?';
        data->c_seq[tax_id]->d_state[site]   = -1;
        data->c_seq[tax_id]->is_ambigu[site] = YES;
      }
    }
  }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void CV_State_Probs_At_Hidden_Positions(phydbl **state_probs, short int **truth, phydbl **weights, int *n_prob_vectors, calign *masked_data, calign *orig_data, t_tree *tree)
{
  t_tree *orig_tree;
  phydbl  r_mat_weight_sum, e_frq_weight_sum, sum_probas;
  int     ns;
  phydbl *Pij, *p_lk_left;

  orig_tree = tree;

  r_mat_weight_sum = -1.;
  e_frq_weight_sum = -1.;
  sum_probas       = -1.;
  
  ns = tree->mod->ns;

  if (tree->is_mixt_tree == YES)
  {
    r_mat_weight_sum =
        MIXT_Get_Sum_Chained_Scalar_Dbl(tree->next->mod->r_mat_weight);
    e_frq_weight_sum =
        MIXT_Get_Sum_Chained_Scalar_Dbl(tree->next->mod->e_frq_weight);
    sum_probas = MIXT_Get_Sum_Of_Probas_Across_Mixtures(
        r_mat_weight_sum, e_frq_weight_sum, tree);
  }

  for (int tax_id = 0; tax_id < tree->n_otu; ++tax_id) // for each row in the alignment
  {
    for (int site = 0; site < tree->data->n_pattern; ++site) // for each column
    {
      if(masked_data->c_seq[tax_id]->state[site] != orig_data->c_seq[tax_id]->state[site]) // position was masked
      {

        // PhyML_Printf("\n. tax_id: %d site: %d", tax_id, site);
        
        assert(!strcmp(tree->a_nodes[tax_id]->name,
                       masked_data->c_seq[tax_id]->name));
        
        // PhyML_Printf("\n. state_probs: %p %p", state_probs, state_probs ? state_probs[0] : NULL);
        
        if (*n_prob_vectors == 0)
        {
          (*state_probs) = (phydbl *)mCalloc(ns, sizeof(phydbl));
          (*truth)       = (short int *)mCalloc(ns, sizeof(short int));
          (*weights)     = (phydbl *)mCalloc(1, sizeof(phydbl));
        }
        else
        {
          (*state_probs) = (phydbl *)mRealloc(
              *state_probs, (*n_prob_vectors + 1) * ns, sizeof(phydbl));
          (*truth) = (short int *)mRealloc(*truth, (*n_prob_vectors + 1) * ns,
                                           sizeof(short int));
          (*weights) =
              (phydbl *)mRealloc(*weights, *n_prob_vectors + 1, sizeof(phydbl));
        }
        
        // PhyML_Printf("\n. state_probs: %p", state_probs);
        
        for (int tip_state = 0; tip_state < ns; ++tip_state)
          (*truth)[*n_prob_vectors * ns + tip_state] = 0;

        (*truth)[*n_prob_vectors * ns + orig_data->c_seq[tax_id]->d_state[site]] = 1;

        (*weights)[*n_prob_vectors] = orig_data->wght[site];

        if (tree->is_mixt_tree == YES)
        {
          for (int tip_state = 0; tip_state < ns; ++tip_state)
            (*state_probs)[*n_prob_vectors * ns + tip_state] = 0.0;

          tree = tree->next;
          do
          {
            Pij       = tree->a_nodes[tax_id]->b[0]->Pij_rr;
            p_lk_left = tree->a_nodes[tax_id]->b[0]->p_lk_left;

            for (int tip_state = 0; tip_state < ns; ++tip_state)
            {
              for (int int_state = 0; int_state < ns; ++int_state)
              {
                (*state_probs)[*n_prob_vectors * ns + tip_state] +=
                    tree->mod->r_mat_weight->v / r_mat_weight_sum *
                    tree->mod->e_frq_weight->v / e_frq_weight_sum / sum_probas *
                    tree->mod->e_frq->pi->v[tip_state] *
                    Pij[tip_state * ns + int_state] * p_lk_left[int_state];
              }
            }
            tree = tree->next;
          } while (tree && tree->is_mixt_tree == NO);

          tree = orig_tree;
          }
        else
        {
          Pij       = tree->a_nodes[tax_id]->b[0]->Pij_rr;
          p_lk_left = tree->a_nodes[tax_id]->b[0]->p_lk_left;

          for (int tip_state = 0; tip_state < ns; ++tip_state)
            (*state_probs)[*n_prob_vectors * ns + tip_state] = 0.0;

          for (int tip_state = 0; tip_state < ns; ++tip_state)
          {
            for (int int_state = 0; int_state < ns; ++int_state)
            {
              for (int catg = 0; catg < tree->mod->ras->n_catg; ++catg)
              {
                (*state_probs)[*n_prob_vectors + tip_state] +=
                    tree->mod->ras->gamma_r_proba->v[catg] *
                    tree->mod->e_frq->pi->v[tip_state] *
                    Pij[catg * ns * ns + tip_state * ns + int_state] *
                    p_lk_left[catg * ns + int_state];
              }
            }
          }
        }
        *n_prob_vectors += 1;        
      }
    }
  }
}