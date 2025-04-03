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
      tree->mod->modelname->s,NULL);

  Free(weights);
  
  return (state_probs_all_sites);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// Mask positions uniformly at random (more than one position may
// be masked at any given site)
void CV_Hide_Align_At_Random_Pos(calign *data, phydbl mask_prob)
{
  phydbl u;

  data->n_masked = 0;

  for (int site = 0; site < data->n_pattern; ++site)
  {
    for (int tax_id = 0; tax_id < data->n_otu; ++tax_id)
    {
      u = Uni();

      if (!(u > mask_prob))
      {
        data->c_seq[tax_id]->state[site]     = '?';
        data->c_seq[tax_id]->d_state[site]   = -1;
        data->c_seq[tax_id]->is_ambigu[site] = YES;

        if (data->n_masked == 0)
          data->masked_pos = (int *)mCalloc(data->n_pattern * 100, sizeof(int));
        // data->masked_pos = (int *)mCalloc(1, sizeof(int));
        // else
        //   data->masked_pos = (int *)mRealloc(data->masked_pos,
        //                                      data->n_masked + 1, sizeof(int));

        data->masked_pos[data->n_masked] = tax_id * data->n_pattern + site;
        data->n_masked++;

        // PhyML_Printf("\n. Mask @ %d %d", tax_id, site);
      }
    }
  }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// Mask exactly one position per site 
void CV_Hide_Align_At_Random_One_Per_Site(calign *data)
{
  int tax_id;

  data->n_masked = 0;

  for (int site = 0; site < data->n_pattern; ++site)
  {
    tax_id = Rand_Int(0, data->n_otu - 1);

    data->c_seq[tax_id]->state[site]     = '?';
    data->c_seq[tax_id]->d_state[site]   = -1;
    data->c_seq[tax_id]->is_ambigu[site] = YES;

    if (data->n_masked == 0)
      data->masked_pos = (int *)mCalloc(data->n_pattern * 100, sizeof(int));
    // data->masked_pos = (int *)mCalloc(1, sizeof(int));
    // else
    //   data->masked_pos =
    //       (int *)mRealloc(data->masked_pos, data->n_masked + 1, sizeof(int));

    data->masked_pos[data->n_masked] = tax_id * data->n_pattern + site;
    data->n_masked++;

    // PhyML_Printf("\n. Mask @ %d %d", tax_id, site);
  }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void CV_Hide_Align_At_Random_Col(calign *data, phydbl mask_prob)
{
  phydbl u;

  data->n_masked = 0;
  
  for (int site = 0; site < data->n_pattern; ++site)
  {
    u = Uni();

    if (!(u > mask_prob))
    {
      for (int tax_id = 0; tax_id < data->n_otu; ++tax_id)
      {
        data->c_seq[tax_id]->state[site]     = '?';
        data->c_seq[tax_id]->d_state[site]   = -1;
        data->c_seq[tax_id]->is_ambigu[site] = YES;
      }

      assert(data->n_masked >= 0);

      if (data->n_masked == 0)
      {
        // PhyML_Printf("\n. First alloc @ %p", data);
        // data->masked_pos = (int *)mCalloc(1, sizeof(int));
        data->masked_pos = (int *)mCalloc(data->n_pattern * 100, sizeof(int));
      }
      // else
      // {
        // PhyML_Printf("\n. Realloc @ %p", data);
        // data->masked_pos =
        //     (int *)mRealloc(data->masked_pos, data->n_masked + 1, sizeof(int));
      // }
      
      data->masked_pos[data->n_masked] = site;
      data->n_masked++;
      // PhyML_Printf("\n. Hiding colmun at position %4d", site);
    }
  }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void CV_Hide_Align_At_Given_Pos(calign *data, int tax_id, int site)
{
  data->c_seq[tax_id]->state[site]     = '?';
  data->c_seq[tax_id]->d_state[site]   = -1;
  data->c_seq[tax_id]->is_ambigu[site] = YES;

  assert(data->n_masked >= 0);
  
  if (data->n_masked == 0)
    // data->masked_pos = (int *)mCalloc(1, sizeof(int));
  data->masked_pos = (int *)mCalloc(data->n_pattern * 100, sizeof(int));
  // else data->masked_pos =
  //     (int *)mRealloc(data->masked_pos, data->n_masked + 1, sizeof(int));

  data->masked_pos[data->n_masked] = tax_id * data->n_pattern + site;
  data->n_masked++;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void CV_State_Probs_At_Hidden_Positions(phydbl **state_probs, short int **truth,
                                        phydbl **site_loglk, phydbl **weights,
                                        int *n_prob_vectors, t_tree *tree)
{

  int tax_id, site;
  
  for (int m = 0; m < tree->data->n_masked; ++m)
  {
    tax_id = floor((phydbl)tree->data->masked_pos[m] / tree->data->n_pattern);
    site   = tree->data->masked_pos[m] - tax_id * tree->data->n_pattern;

    CV_State_Probs_Core(state_probs, truth, site_loglk, weights, n_prob_vectors,
                        tax_id, site, tree->data->c_seq[tax_id]->d_state[site],
                        tree->data->wght[site], tree);
  }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void CV_State_Probs_Core(phydbl **state_probs, short int **truth,
                         phydbl **site_loglk, phydbl **weights,
                         int *n_prob_vectors, int tax_id, int site,
                         int true_d_state, phydbl patt_weight, t_tree *tree)
{
  int     ns;
  phydbl *Pij, *p_lk_left, sum;
  
  
  ns = tree->mod->ns;
  
  assert(*n_prob_vectors >= 0);

  if (*n_prob_vectors == 0)
  {
    // (*state_probs) = (phydbl *)mCalloc(ns,
    //                                    sizeof(phydbl));
    // (*site_loglk)  = (phydbl *)mCalloc(1,
    //                                    sizeof(phydbl));
    // (*truth)   = (short int *)mCalloc(ns,
    //                                   sizeof(short int));
    // (*weights) = (phydbl *)mCalloc(1,
    //                                sizeof(phydbl));
    (*state_probs) = (phydbl *)mCalloc(ns * tree->data->n_pattern * tree->n_otu,
                                       sizeof(phydbl));
    (*site_loglk)  = (phydbl *)mCalloc(1 * tree->data->n_pattern * tree->n_otu,
                                       sizeof(phydbl));
    (*truth) = (short int *)mCalloc(ns * tree->data->n_pattern * tree->n_otu,
                                    sizeof(short int));
    (*weights)     = (phydbl *)mCalloc(1 * tree->data->n_pattern * tree->n_otu,
                                       sizeof(phydbl));
  }
  else
   {
    (*state_probs) = (phydbl *)mRealloc(
        *state_probs, (*n_prob_vectors + 1) * ns, sizeof(phydbl));
    // (*site_loglk) =
    //     (phydbl *)mRealloc(*site_loglk, *n_prob_vectors + 1, sizeof(phydbl));
    // (*truth) = (short int *)mRealloc(*truth, (*n_prob_vectors + 1) * ns,
    //                                  sizeof(short int));
    // (*weights) =
    //     (phydbl *)mRealloc(*weights, *n_prob_vectors + 1, sizeof(phydbl));
  }
  // PhyML_Printf("\n. state_probs: %p", state_probs);

  for (int tip_state = 0; tip_state < ns; ++tip_state)
    (*truth)[*n_prob_vectors * ns + tip_state] = 0;

  (*truth)[*n_prob_vectors * ns + true_d_state] = 1;

  (*weights)[*n_prob_vectors] = patt_weight;

  (*site_loglk)[*n_prob_vectors] = tree->c_lnL_sorted[site];

  if (tree->is_mixt_tree == YES)
  {
    phydbl r_mat_weight_sum =
        MIXT_Get_Sum_Chained_Scalar_Dbl(tree->next->mod->r_mat_weight);

    phydbl e_frq_weight_sum =
        MIXT_Get_Sum_Chained_Scalar_Dbl(tree->next->mod->e_frq_weight);

    phydbl sum_probas = MIXT_Get_Sum_Of_Probas_Across_Mixtures(
        r_mat_weight_sum, e_frq_weight_sum, tree);

    for (int tip_state = 0; tip_state < ns; ++tip_state)
      (*state_probs)[*n_prob_vectors * ns + tip_state] = 0.0;

    tree = tree->next;
    do
    {
      Pij       = tree->a_nodes[tax_id]->b[0]->Pij_rr;
      p_lk_left = tree->a_nodes[tax_id]->b[0]->p_lk_left + site * ns;

      // PhyML_Printf("\n. ncatg: %d truth: %d p_lk_left: %f %f %f %f Pij: %f %f
      // %f %f",
      //              tree->mod->ras->n_catg,
      //              orig_data->c_seq[tax_id]->d_state[site],
      //              p_lk_left[0],
      //              p_lk_left[1], p_lk_left[2], p_lk_left[3],
      //              Pij[orig_data->c_seq[tax_id]->d_state[site] * ns + 0],
      //              Pij[orig_data->c_seq[tax_id]->d_state[site] * ns + 1],
      //              Pij[orig_data->c_seq[tax_id]->d_state[site] * ns + 2],
      //              Pij[orig_data->c_seq[tax_id]->d_state[site] * ns + 3]);

      for (int tip_state = 0; tip_state < ns; ++tip_state)
      {
        for (int int_state = 0; int_state < ns; ++int_state)
        {
          (*state_probs)[*n_prob_vectors * ns + tip_state] +=
              tree->mixt_tree->mod->ras->gamma_r_proba
                  ->v[tree->mod->ras->parent_class_number] *
              tree->mod->r_mat_weight->v / r_mat_weight_sum *
              tree->mod->e_frq_weight->v / e_frq_weight_sum / sum_probas *
              tree->mod->e_frq->pi->v[tip_state] *
              Pij[tip_state * ns + int_state] * p_lk_left[int_state];
        }
      }

      tree = tree->next;
    } while (tree && tree->is_mixt_tree == NO);
    // We only consider one partition element here.
    // Calls to this function are required for every
    // partition element.
  }
  else
  {
    Pij       = tree->a_nodes[tax_id]->b[0]->Pij_rr;
    p_lk_left = tree->a_nodes[tax_id]->b[0]->p_lk_left +
                site * ns * tree->mod->ras->n_catg;

    for (int tip_state = 0; tip_state < ns; ++tip_state)
      (*state_probs)[*n_prob_vectors * ns + tip_state] = 0.0;

    for (int tip_state = 0; tip_state < ns; ++tip_state)
    {
      for (int int_state = 0; int_state < ns; ++int_state)
      {
        for (int catg = 0; catg < tree->mod->ras->n_catg; ++catg)
        {
          (*state_probs)[*n_prob_vectors * ns + tip_state] +=
              tree->mod->ras->gamma_r_proba->v[catg] *
              tree->mod->e_frq->pi->v[tip_state] *
              Pij[catg * ns * ns + tip_state * ns + int_state] *
              p_lk_left[catg * ns + int_state];
        }
      }
    }
  }

  sum = 0.0;
  for (int tip_state = 0; tip_state < ns; ++tip_state)
    sum += (*state_probs)[*n_prob_vectors * ns + tip_state];

  for (int tip_state = 0; tip_state < ns; ++tip_state)
    (*state_probs)[*n_prob_vectors * ns + tip_state] /= sum;

  // PhyML_Printf("\n. state_probs: %f %f %f %f truth: %d",
  //              (*state_probs)[*n_prob_vectors * ns + 0],
  //              (*state_probs)[*n_prob_vectors * ns + 1],
  //              (*state_probs)[*n_prob_vectors * ns + 2],
  //              (*state_probs)[*n_prob_vectors * ns + 3],
  //              orig_data->c_seq[tax_id]->d_state[site]);

  *n_prob_vectors += 1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void CV_Score_At_Hidden_Cols(phydbl **site_loglk, phydbl **weights,
                             int *n_prob_vectors, t_tree *tree)
{

  int site;

  for (int m = 0; m < tree->data->n_masked; ++m)
  {
    site = tree->data->masked_pos[m];

    // PhyML_Printf("\n. CV prob @ site %d", site);

assert(*n_prob_vectors >= 0);

    if (*n_prob_vectors == 0)
    {
      (*site_loglk) = (phydbl *)mCalloc(1, sizeof(phydbl));
      (*weights)    = (phydbl *)mCalloc(1, sizeof(phydbl));
    }
    else
    {
      (*site_loglk) =
          (phydbl *)mRealloc(*site_loglk, *n_prob_vectors + 1, sizeof(phydbl));
      (*weights) =
          (phydbl *)mRealloc(*weights, *n_prob_vectors + 1, sizeof(phydbl));
    }

    (*site_loglk)[*n_prob_vectors] = tree->c_lnL_sorted[site];
    (*weights)[*n_prob_vectors]    = tree->data->wght[site];

    // PhyML_Printf(" lnL: %15g w: %12f #: %d site: %d", (*site_loglk)[*n_prob_vectors],
    //              (*weights)[*n_prob_vectors],*n_prob_vectors,site);

    (*n_prob_vectors)++;
  }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
