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

    b = NULL;

    for(int i=0; i<tree->n_otu; ++i)
    {
        b = tree->a_nodes[i]->b[0];

        for(int j=0; j<tree->n_pattern; ++j)
        {
            
        }
    }



}