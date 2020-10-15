/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include <config.h>

#ifndef LOCATION_H
#define LOCATION_H

#include "utilities.h"

phydbl LOCATION_Lk(t_tree *tree);
phydbl LOCATION_Lk_Range(t_dsk *young, t_dsk *old, t_tree *tree);
phydbl LOCATION_Lk_Core(t_dsk *disk, t_tree *tree);


#endif
