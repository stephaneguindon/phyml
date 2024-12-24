/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include "cv.h"

void CV_Tip_Cv(t_tree *tree)
{
    t_edge *b;
    int obs_d_state, ambiguity_check;
    char obs_state;
    int ns, nsncatg, n_otu;
    phydbl *p_lk_left, *Pij, *state_probs_all_sites;
    phydbl sum;
    short int *truth;
    

    Set_Both_Sides(YES, tree);
    Lk(NULL, tree);


    b = NULL;
    p_lk_left = NULL;
    state_probs_all_sites = (phydbl *)mCalloc(tree->n_pattern * tree->n_otu * tree->mod->ns, sizeof(phydbl));
    truth = (short int *)mCalloc(tree->n_pattern * tree->n_otu * tree->mod->ns, sizeof(phydbl));
    
    ns = tree->mod->ns;
    nsncatg = ns * tree->mod->ras->n_catg;
    n_otu = tree->n_otu;


    for (int tax_id = 0; tax_id < tree->n_otu; ++tax_id)
    {
        b = tree->a_nodes[tax_id]->b[0];

        p_lk_left = b->p_lk_left;
        Pij = b->Pij_rr;

        for (int site = 0; site < tree->n_pattern; ++site)
        {
            if (tree->data->wght[site] > SMALL)
            {
                ambiguity_check = tree->a_nodes[tax_id]->c_seq->is_ambigu[site];
                if (ambiguity_check == NO)
                {
                    obs_state = tree->a_nodes[tax_id]->c_seq->state[site];
                    obs_d_state = tree->a_nodes[tax_id]->c_seq->d_state[site];

                    tree->a_nodes[tax_id]->c_seq->is_ambigu[site] = YES;
                    tree->a_nodes[tax_id]->c_seq->state[site] = '?';
                    tree->a_nodes[tax_id]->c_seq->d_state[site] = -1;

                    Init_Partial_Lk_Tips_Double_One_Character(tax_id, site, tree);

                    Br_Len_Opt(&(b->l->v), b, tree);
                    // Optimize_Br_Len_Serie(1000,tree);
                    // Optimiz_All_Free_Param(tree,NO);

                    for (int tip_state = 0; tip_state < ns; ++tip_state)
                    {
                      state_probs_all_sites[site * ns * n_otu + tax_id * ns + tip_state] = 0.0;
                      for (int int_state = 0; int_state < ns; ++int_state)
                      {
                        for (int catg = 0; catg < tree->mod->ras->n_catg; ++catg)
                        {
                          state_probs_all_sites[site * ns * n_otu + tax_id * ns + tip_state] +=
                              tree->mod->ras->gamma_r_proba->v[catg] *
                              tree->mod->e_frq->pi->v[tip_state] *
                              Pij[catg * ns * ns + tip_state * ns + int_state] *
                              p_lk_left[catg * ns + int_state];                        
                        }
                      }
                    }

                    sum = 0.0;
                    for (int state = 0; state < ns; ++state)
                      sum += state_probs_all_sites[site * ns * n_otu + tax_id * ns + state];

                    for (int state = 0; state < ns; ++state)
                      state_probs_all_sites[site * ns * n_otu + tax_id * ns + state] /= sum;


                    for (int state = 0; state < ns; ++state)
                      truth[site * ns * n_otu + tax_id * ns + state] = 0;

                    truth[site * ns * n_otu + tax_id * ns + obs_d_state] = 1;

                    PhyML_Printf("\n###%s,%s,%d,%g", tree->mod->modelname->s, tree->a_nodes[tax_id]->name, site, log(state_probs_all_sites[site * ns * n_otu + tax_id * ns + obs_d_state]));

                    tree->a_nodes[tax_id]->c_seq->state[site] = obs_state;
                    tree->a_nodes[tax_id]->c_seq->d_state[site] = obs_d_state;
                    tree->a_nodes[tax_id]->c_seq->is_ambigu[site] = NO;

                    Init_Partial_Lk_Tips_Double_One_Character(tax_id, site, tree);
                }
            }
            p_lk_left += nsncatg;
        }
    }

    ROC(state_probs_all_sites, truth, ns, tree->n_otu * tree->n_pattern, tree->mod->modelname->s);
    
}
