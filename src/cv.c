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
    int obs_d_state,ambiguity_check;
    char obs_state;
    int ns,nsncatg;
    phydbl *p_lk_left,*Pij,*state_prob;
    phydbl sum;
    int argmax_state_prob;
    int n_errors,n_characters,n_correct;
    phydbl prob_error,prob_truth;

    Set_Both_Sides(YES, tree);
    Lk(NULL, tree);

    
    b = NULL;
    p_lk_left = NULL;
    state_prob = (phydbl *)mCalloc(tree->mod->ns,sizeof(phydbl));

    ns = tree->mod->ns;
    nsncatg = ns * tree->mod->ras->n_catg;

    n_errors = 0;
    n_correct = 0;
    n_characters = 0;
    prob_error = 0.0;
    prob_truth = 0.0;

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

                    Br_Len_Opt(&(b->l->v),b,tree);

                    for (int tip_state = 0; tip_state < ns; ++tip_state)
                    {
                        state_prob[tip_state] = 0.0;
                        for (int int_state = 0; int_state < ns; ++int_state)
                        {
                            for (int catg = 0; catg < tree->mod->ras->n_catg; ++catg)
                            {
                                state_prob[tip_state] +=
                                    tree->mod->ras->gamma_r_proba->v[catg] *
                                    tree->mod->e_frq->pi->v[tip_state] *
                                    Pij[catg * ns * ns + tip_state * ns + int_state] *
                                    p_lk_left[catg * ns + int_state];
                            }
                        }
                    }

                    sum = 0.0;
                    argmax_state_prob = 0;
                    for (int state = 0; state < ns; ++state) 
                    {
                        sum += state_prob[state];
                        if (state_prob[state] > state_prob[argmax_state_prob])
                            argmax_state_prob = state;
                    }

                    if (argmax_state_prob != obs_d_state)
                    {
                        for (int state = 0; state < ns; ++state)
                        {
                            PhyML_Printf("\n Tax.: %s site: %d current: %c obs.state: %c obs. dstate: %d predicted: %c prob: %f%c%c",
                                         tree->a_nodes[tax_id]->name,
                                         site + 1, 
                                         D_State_To_Character(state, tree),
                                         obs_state, obs_d_state, 
                                         D_State_To_Character(argmax_state_prob, tree),
                                         state_prob[state] / sum,
                                         (state == argmax_state_prob) ? '*' : ' ',
                                         (state == obs_d_state) ? '!' : ' ');
                        }
                        prob_error += state_prob[argmax_state_prob]/sum;                        
                        n_errors++;
                    }
                    else
                    {
                        prob_truth += state_prob[argmax_state_prob]/sum;
                        n_correct++;
                    }

                    n_characters++;
                    tree->a_nodes[tax_id]->c_seq->state[site] = obs_state;
                    tree->a_nodes[tax_id]->c_seq->d_state[site] = obs_d_state;
                    tree->a_nodes[tax_id]->c_seq->is_ambigu[site] = NO;

                    Init_Partial_Lk_Tips_Double_One_Character(tax_id, site, tree);
                }
            }
            p_lk_left += nsncatg;
        }
    }
    
    PhyML_Printf("\n. Avg. prob of error: %G",prob_error/n_errors);
    PhyML_Printf("\n. Avg. prob of correct: %G",prob_truth/n_correct);
    PhyML_Printf("\n. Error frq.: %G",(phydbl)n_errors/n_characters);

    Free(state_prob);
}
