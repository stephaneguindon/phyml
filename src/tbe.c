/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include "tbe.h"
#include "utilities.h"

/* UNION AND INTERSECT CALCULATIONS (FOR THE TRANSFER METHOD) */
void Update_IC_Ref_Tree(t_tree *ref_tree, t_node * orig, t_node* target, t_edge *my_br, t_tree *boot_tree,
			short unsigned** i_matrix, short unsigned** c_matrix, int* cluster_sizes){
  /* this function does the post-order traversal (recursive from the pseudoroot to the leaves, updating knowledge for the subtrees)
     of the reference tree, examining only leaves (terminal edges) of the bootstrap tree.
     It sends a probe from the orig node to the target node (nodes in ref_tree), calculating I_ij and C_ij
     (see Brehelin, Gascuel, Martin 2008). */
  int j, k;
  int edge_id; /* its id */
  int next_edge_id;
  t_node *tip;
  int boot_edge_id;
  edge_id = my_br->num; /* all this is in ref_tree */
  assert(target==my_br->rght); /* the descendant should always be the right side of the edge ? */
  
  if(target->tax) {
    cluster_sizes[edge_id] = 1;
    for (j=0; j < 2*boot_tree->n_otu-3; j++) { /* for all the terminal edges of boot_tree */
      boot_edge_id = boot_tree->a_edges[j]->num;
      // If not tip, continue
      tip = NULL;
      if(boot_tree->a_edges[j]->rght->tax)
	tip = boot_tree->a_edges[j]->rght;
      if(boot_tree->a_edges[j]->left->tax)
	tip = boot_tree->a_edges[j]->left;
      if(tip == NULL) continue;
      /* we only want to scan terminal edges of boot_tree, where the right son is a leaf */
      /* else we update all the I_ij and C_ij with i = edge_id */
      if (strcmp(target->name,tip->name)) {
	/* here the taxa are different */
	i_matrix[edge_id][boot_edge_id] = 0;
	c_matrix[edge_id][boot_edge_id] = 1;
      } else {
	/* same taxa here in T_ref and T_boot */
	i_matrix[edge_id][boot_edge_id] = 1;
	c_matrix[edge_id][boot_edge_id] = 0;
      }
    } /* end for on all edges of T_boot, for my_br being terminal */
  } else {
    cluster_sizes[edge_id] = 0;
    /* now the case where my_br is not a terminal edge */
    /* first initialise (zero) the cells we are going to update */
    for (j=0; j < 2*boot_tree->n_otu-3; j++){
      /* We initialize the i and c matrices for the edge edge_id with :
       * 0 for i : because afterwards we do i[edge_id] = i[edge_id] || i[next_edge_id]
       * 1 for c : because afterwards we do c[edge_id] = c[edge_id] && c[next_edge_id]
       */
      if(boot_tree->a_edges[j]->rght->tax || boot_tree->a_edges[j]->left->tax){
	boot_edge_id = boot_tree->a_edges[j]->num;
 	i_matrix[edge_id][boot_edge_id] = 0;
	c_matrix[edge_id][boot_edge_id] = 1;
      }
    }
    
    for (k = 0; k < 3; k++) {
      if(target->v[k]!=orig){
	Update_IC_Ref_Tree(ref_tree, target, target->v[k], target->b[k], boot_tree, i_matrix, c_matrix, cluster_sizes);
	next_edge_id = target->b[k]->num;
	cluster_sizes[edge_id] += cluster_sizes[next_edge_id];
	for (j=0; j < 2*boot_tree->n_otu-3; j++) { /* for all the terminal edges of boot_tree */
	  boot_edge_id = boot_tree->a_edges[j]->num;
	  // If not a tip, continue
	  if(!(boot_tree->a_edges[j]->rght->tax) && !(boot_tree->a_edges[j]->left->tax)) continue;
	  i_matrix[edge_id][boot_edge_id] = i_matrix[edge_id][boot_edge_id] || i_matrix[next_edge_id][boot_edge_id];
	  /* above is an OR between two integers, result is 0 or 1 */
	  c_matrix[edge_id][boot_edge_id] = c_matrix[edge_id][boot_edge_id] && c_matrix[next_edge_id][boot_edge_id];
	  /* above is an AND between two integers, result is 0 or 1 */
	} /* end for j */
      }
    } /* end for on all edges of T_boot, for my_br being internal */
  } /* ending the case where my_br is an internal edge */
} /* end update_i_c_post_order_ref_tree */


void Update_All_IC_Ref_Tree(t_tree* ref_tree, t_tree* boot_tree,
			    short unsigned** i_matrix, short unsigned** c_matrix, int* cluster_sizes) {
  /* this function is the first step of the union and intersection calculations */
  int i;
  t_node *root;
  
  if(!ref_tree->n_root){
    i = 0;
    while((!ref_tree->a_nodes[ref_tree->n_otu+i]->v[0]) ||
	  (!ref_tree->a_nodes[ref_tree->n_otu+i]->v[1]) ||
	  (!ref_tree->a_nodes[ref_tree->n_otu+i]->v[2])) i++;
    root=ref_tree->a_nodes[ref_tree->n_otu+i];
    Update_IC_Ref_Tree(ref_tree, root, root->v[0], root->b[0], boot_tree, i_matrix, c_matrix, cluster_sizes);
    Update_IC_Ref_Tree(ref_tree, root, root->v[1], root->b[1], boot_tree, i_matrix, c_matrix, cluster_sizes);
    Update_IC_Ref_Tree(ref_tree, root, root->v[2], root->b[2], boot_tree, i_matrix, c_matrix, cluster_sizes);
  } else {
    root=ref_tree->n_root;
   Update_IC_Ref_Tree(ref_tree, root, root->v[0], root->b[0], boot_tree, i_matrix, c_matrix, cluster_sizes);
   Update_IC_Ref_Tree(ref_tree, root, root->v[1], root->b[1], boot_tree, i_matrix, c_matrix, cluster_sizes);
  }
} /* end update_all_i_c_post_order_ref_tree */

void Update_IC_Boot_Tree(t_tree* ref_tree, t_tree* boot_tree, t_node* orig, t_node* target,
			 t_edge *my_br, short unsigned** i_matrix, short unsigned** c_matrix,
			 short unsigned** hamming, short unsigned* min_dist,
			 short unsigned* min_dist_edge, int* cluster_sizes) {
  /* here we implement the second part of the Brehelin/Gascuel/Martin algorithm:
     post-order traversal of the bootstrap tree, and numerical recurrence. 
     in this function, orig and target are nodes of boot_tree (aka T_boot).
     min_dist is an array whose size is equal to the number of edges in T_ref.
     It gives for each edge of T_ref its min distance to a split in T_boot. */
  int i, k;
  int boot_edge_id /* its id */, next_boot_edge_id /* id of descending branches. */;
  int N = ref_tree->n_otu;
  /* we first have to determine which is the direction of the edge (orig -> target and target -> orig) */
  boot_edge_id = my_br->num; /* here this is an edge_id corresponding to T_boot */

  if(!target->tax){

    /* because nothing to do in the case where target is a leaf: intersection and union already ok. */
    /* otherwise, keep on posttraversing in all other directions */
    
    /* first initialise (zero) the cells we are going to update */
    for (i=0; i < 2*ref_tree->n_otu-3; i++) i_matrix[i][boot_edge_id] = c_matrix[i][boot_edge_id] = 0;
    
    for(k=0;k<3;k++) {
      if(target->v[k] != orig){
	next_boot_edge_id = target->b[k]->num;
	Update_IC_Boot_Tree(ref_tree, boot_tree, target, target->v[k], target->b[k],
					i_matrix, c_matrix, hamming, min_dist, min_dist_edge, cluster_sizes);
	for (i=0; i < 2*ref_tree->n_otu-3; i++) { /* for all the edges of ref_tree */ 
	  i_matrix[i][boot_edge_id] += i_matrix[i][next_boot_edge_id];
	  c_matrix[i][boot_edge_id] += c_matrix[i][next_boot_edge_id];
	} /* end for i */
      }
    } 
  } /* end if target is not a leaf: the following loop is performed in all cases */

  for (i=0; i< 2*ref_tree->n_otu-3; i++) { /* for all the edges of ref_tree */ 
    /* at this point we can calculate in all cases (internal branch or not) the Hamming distance at [i][boot_edge_id], */
    /* card of union minus card of intersection */
    hamming[i][boot_edge_id] = cluster_sizes[i] /* #taxa in the cluster i of T_ref */
      + c_matrix[i][boot_edge_id] /* #taxa in cluster edge_id of T_boot BUT NOT in cluster i of T_ref */
      - i_matrix[i][boot_edge_id]; /* #taxa in the intersection of the two clusters */
    
    /* Let's immediately calculate the right ditance, taking into account the fact that the true disance is min (dist, N-dist) */
    if (hamming[i][boot_edge_id] > N/2 /* floor value */) hamming[i][boot_edge_id] = N - hamming[i][boot_edge_id];
	  
    /*   and update the min of all Hamming (TRANSFER) distances hamming[i][j] over all j */
    if (hamming[i][boot_edge_id] < min_dist[i]){
      min_dist[i] = hamming[i][boot_edge_id];
      min_dist_edge[i] = boot_edge_id;
    }
  } /* end for on all edges of T_ref */
} /* end update_i_c_post_order_boot_tree */


void Update_All_IC_Boot_Tree(t_tree* ref_tree, t_tree* boot_tree, short unsigned** i_matrix,
			     short unsigned** c_matrix, short unsigned** hamming,
			     short unsigned* min_dist, short unsigned* min_dist_edge, int* cluster_sizes) {
	/* this function is the second step of the union and intersection calculations */
  int i;
  t_node *root;
  if(!boot_tree->n_root){
    i = 0;
    while((!boot_tree->a_nodes[boot_tree->n_otu+i]->v[0]) ||
	  (!boot_tree->a_nodes[boot_tree->n_otu+i]->v[1]) ||
	  (!boot_tree->a_nodes[boot_tree->n_otu+i]->v[2])) i++;
    root=boot_tree->a_nodes[boot_tree->n_otu+i];
    Update_IC_Boot_Tree(ref_tree, boot_tree, root, root->v[0], root->b[0], i_matrix, c_matrix, hamming, min_dist, min_dist_edge, cluster_sizes);
    Update_IC_Boot_Tree(ref_tree, boot_tree, root, root->v[1], root->b[1], i_matrix, c_matrix, hamming, min_dist, min_dist_edge, cluster_sizes);
    Update_IC_Boot_Tree(ref_tree, boot_tree, root, root->v[2], root->b[2], i_matrix, c_matrix, hamming, min_dist, min_dist_edge, cluster_sizes);
  } else {
    root=boot_tree->n_root;
    Update_IC_Boot_Tree(ref_tree, boot_tree, root, root->v[0], root->b[0], i_matrix, c_matrix, hamming, min_dist, min_dist_edge, cluster_sizes);
    Update_IC_Boot_Tree(ref_tree, boot_tree, root, root->v[1], root->b[2], i_matrix, c_matrix, hamming, min_dist, min_dist_edge, cluster_sizes);
  }
  
  /* and then some checks to make sure everything went ok */
  for(i=0; i<2*ref_tree->n_otu-3; i++) {
    assert(min_dist[ref_tree->a_edges[i]->num] >= 0);
    if(ref_tree->a_edges[i]->rght->tax || ref_tree->a_edges[i]->left->tax){
      assert(min_dist[ref_tree->a_edges[i]->num] == 0); /* any terminal edge should have an exact match in any bootstrap tree */
    }
  }
} /* end update_all_i_c_post_order_boot_tree */
