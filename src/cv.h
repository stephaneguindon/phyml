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
void    CV_State_Probs_At_Hidden_Positions(phydbl *state_probs, short int *truth,
                                           phydbl *weights, int *n_prob_vectors,
                                           calign *masked_data, calign *orig_data,
                                           t_tree *tree);
void    CV_Hide_Characters_At_Random(calign *data, phydbl mask_prob);

#endif
