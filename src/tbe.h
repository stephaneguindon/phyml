/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#ifndef TBE_H
#define TBE_H

#include "utilities.h"

/* Functions to compute Transfer distances, used in transfert bootstrap */
void Update_All_IC_Ref_Tree(t_tree* ref_tree, t_tree* boot_tree,
			    short unsigned** i_matrix, short unsigned** c_matrix, int* cluster_sizes);
void Update_IC_Ref_Tree(t_tree *ref_tree, t_node * orig, t_node* target, t_edge *my_br, t_tree *boot_tree,
			short unsigned** i_matrix, short unsigned** c_matrix, int* cluster_sizes);
void Update_IC_Boot_Tree(t_tree* ref_tree, t_tree* boot_tree, t_node* orig, t_node* target,
			 t_edge *my_br, short unsigned** i_matrix, short unsigned** c_matrix,
			 short unsigned** hamming, short unsigned* min_dist,
			 short unsigned* min_dist_edge, int* cluster_sizes);
void Update_All_IC_Boot_Tree(t_tree* ref_tree, t_tree* boot_tree,
			     short unsigned** i_matrix, short unsigned** c_matrix,
			     short unsigned** hamming, short unsigned* min_dist,
			     short unsigned* min_dist_edge, int* cluster_sizes);
#endif
