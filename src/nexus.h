/*

PHYML :  a program that  computes maximum likelihood  phylogenies from
DNA or AA homologous sequences 

Copyright (C) Stephane Guindon. Oct 2003 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/

#include <config.h>

#ifndef NEXUS_H
#define NEXUS_H

#include "utilities.h"

void Find_Nexus_Com(char *token,nexcom **found_com,nexparm **default_parm,nexcom **com_list);
void Find_Nexus_Parm(char *token,nexparm **found_parm,nexcom *curr_com);
int Read_Nexus_Taxa(char *token,nexparm *curr_parm,option *io);
int Read_Nexus_Translate(char *token,nexparm *curr_parm,option *io);
int Read_Nexus_Matrix(char *token,nexparm *curr_parm,option *io);
int Read_Nexus_Tree(char *token,nexparm *curr_parm,option *io);
int Read_Nexus_Begin(char *token,nexparm *curr_parm,option *io);
int Read_Nexus_Dimensions(char *token,nexparm *curr_parm,option *io);
int Read_Nexus_Format(char *token,nexparm *curr_parm,option *io);
int Read_Nexus_Eliminate(char *token,nexparm *curr_parm,option *io);
int Read_Nexus_Taxlabel(char *token,nexparm *curr_parm,option *io);
int Read_Nexus_Charstatelabels(char *token,nexparm *curr_parm,option *io);
int Read_Nexus_Charlabels(char *token,nexparm *curr_parm,option *io);
int Read_Nexus_Statelabels(char *token,nexparm *curr_parm,option *io);

#endif
