/*

PHYML :  a program that  computes maximum likelihood  phylogenies from
DNA or AA homologous sequences 

Copyright (C) Stephane Guindon. Oct 2003 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/

#include <config.h>

#ifndef CV_H
#define CV_H

#include "utilities.h"

phydbl *CV_Tip_Cv(t_tree *tree);
void CV_State_Probs_At_Hidden_Positions(phydbl **state_probs, short int **truth,
                                        phydbl **site_loglk, phydbl **weights,
                                        int *n_prob_vectors, t_tree *tree);
void CV_Hide_Align_At_Random_Pos(calign *data, phydbl mask_prob);
void CV_Hide_Align_At_Random_Col(calign *data, phydbl mask_prob);
void CV_Hide_Align_At_Given_Pos(calign *data, int tax_id, int site);
void CV_Hide_Align_At_Random_One_Per_Site(calign *data);
void CV_State_Probs_Core(phydbl **state_probs, short int **truth,
                         phydbl **site_loglk, phydbl **weights,
                         int *n_prob_vectors, int tax_id, int site,
                         int true_d_state, phydbl patt_weight, t_tree *tree);
void CV_Score_At_Hidden_Cols(phydbl **site_loglk, phydbl **weights,
                             int *n_prob_vectors, t_tree *tree);

#endif
