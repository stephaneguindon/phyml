/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include <config.h>

#ifndef NJ_H
#define NJ_H

#include "utilities.h"
#include "optimiz.h"
#include "free.h"

void   Bionj(matrix *mat);
void   Finish(matrix *mat);
void   Compute_Sx(matrix *mat);
phydbl Sum_S(matrix *mat, int i);
phydbl Dist(matrix *mat, int x, int y);
phydbl Q_Agglo(matrix *mat, int x, int y);
phydbl Variance(matrix *mat, int x, int y);
phydbl Br_Length(matrix *mat, int x, int y);
void   Update_Dist(matrix *mat, int x, int y);
phydbl Lamda(matrix *mat, int x, int y, phydbl vxy);
void   Best_Pair(matrix *mat, int *x, int *y, phydbl *score);
phydbl Var_Red(matrix *mat, int x, int y, int i, phydbl lamda, phydbl vxy);
void   Update_Tree(matrix *mat, int x, int y, phydbl lx, phydbl ly, phydbl score);
void   Update_Mat(matrix *mat, int x, int y, 
		  phydbl lx, phydbl ly, phydbl vxy, phydbl lamda);
phydbl Dist_Red(matrix *mat, int x, phydbl lx, int y, 
		phydbl ly, int i, phydbl lamda);
int    Bionj_Br_Length_Post(t_node *a, t_node *d, matrix *mat);
void   Bionj_Br_Length(matrix *mat);

#endif
