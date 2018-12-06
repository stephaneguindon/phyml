/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include <config.h>

#ifndef ANCESTRAL_H
#define ANCESTRAL_H

#include "utilities.h"


void Sample_Ancestral_Seq(int fullmutmap, int fromprior, t_tree *tree);
void Sample_Ancestral_Seq_Pre(t_node *a, t_node *d, t_edge *b,
                              int site, int r_cat,
                              int *muttype, phydbl *muttime, int *muttax, int *n_mut,
                              int fullmutmap, int fromprior, t_tree *tree);
int Sample_Ancestral_Seq_Core(t_node *a, t_node *d, t_edge *b, int r_cat, int site, t_tree *tree);
void Map_Mutations(t_node *a, t_node *d, int sa, int sd, t_edge *b, int site, int rcat, int *muttype, phydbl *muttime, int *muttax, int *n_mut, t_tree *tree);
void Ancestral_Sequences(t_tree *tree, int print);
void Ancestral_Sequences_One_Node(t_node *d, t_tree *tree, int print);
int MPEE_Score(const phydbl *alpha, int *idx, const phydbl *p, const int ns);
int MPEE_Infer(const phydbl *p, const int ns);
void Reverse_Muttime(phydbl *muttime, int n_mut, t_tree *tree);

#endif
