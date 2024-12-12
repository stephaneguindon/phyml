/*

PHYML :  a program that  computes maximum likelihood  phylogenies from
DNA or AA homologous sequences 

Copyright (C) Stephane Guindon. Oct 2003 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/

#include <config.h>

#ifndef EVOLVE_H
#define EVOLVE_H

#include "utilities.h"

int EVOLVE_Main(int argc, char **argv);
void EVOLVE_Coalescent(t_tree *tree);

#endif
