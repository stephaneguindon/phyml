/*
 * author: Imran Fanaswala
 */

#ifndef BEAGLE_UTILS_H
#define BEAGLE_UTILS_H

#include "assert.h"
#include "libhmsbeagle/beagle.h"
#include "utilities.h"
#define UNINITIALIZED -42

/*
BEAGLE is a library that abstracts the core computations common to most
phylogenetic packages. The following paragraphs serve as a "design document" for
the implementation details of how BEAGLE has been integrated with PhyML. The
reader is expected to be familiar with phylogenetics, skimmed the BEAGLE API,
and seen a couple of its examples.

PARTIAL LIKELIHOODS:
- For each internal edge, PhyML maintains a vector of partials that lie to the
"right" and "left", to faciliate both post-order and pre-order tree traversal.
Therefore each t_edge has a `p_lk_right` and `p_lk_left` field
- PhyML represents a tip as a vector of partials (where one bit is "1", and
the rest are "0"). Therefore, a t_edge connected to a tip has `p_lk_tip` field.
Also by convention, the tip always lies to the right of the edge its connected
to. In other words, if b->right->tax is True, then b->p_lk_tip holds the tip
partials.
- Therefore for BEAGLE, we allocate N+(2*num_edges) partial buffers.
- Recall that BEAGLE requires each partial buffer to be uniquely indexed. Thus
each edge also has an additional `p_lk_rght_idx` and `p_lk_left_idx` field.
These fields are used to create a BeagleOperation struct during the call to
beagleUpdatePartials().
- In PhyML, Update_P_Lk() updates the partials for a specific edge lieing along
a specific "direction". The "direction" depends on which tree-traversal is being
employed. In other words, the tree-traversal dictates whether we are computing
"up"/"right" or "down"/"left" partials relative to the t_node*. Specifically,
Set_All_P_Lk() handles the nitty gritty of determining whether we are computing
"right" or "left" partials. Therefore this function has been modified to also
set the appropriate BEAGLE indices.

P-MATS and EDGELIKELIHOODS:
- (1) Observe that PhyML does eigen decomposition for only GTR and AA models and
thus the P-matrix is computed using Beagle. This implies the usage of
SetEigenDecomposition() and UpdateTransitionMatricies() calls.
- (2) On the other hand, in case of simpler models (JC69, F81, etc), the
P-Matrix is explicitly set in BEAGLE thus implying the usage of
SetTransitionMatrix().
- The choice between (1) and (2) is made in Update_PMat_At_Given_Edge(). BTW,
upon reading Update_PMat_At_Given_Edge(), it may appear that the P-Matrix is
being computed in PhyML *and* in BEAGLE; afterall because PMat() function in all
cases. However, when BEAGLE is enabled, observe that in PMat() essentially does
nothing.
- Lk_Core() computes the edgelikelihoods, and calc_edgelks_beagle() is the
corresponding BEAGLE function that does the same.


*/



#define TODO_BEAGLE "TODO. This codepath has not been implemented in PhyML-X, please post your usecase on the PhyML discussion list"

int  create_beagle_instance(t_tree* tree, int quiet, option* io);
int  finalize_beagle_instance(t_tree* tree);
void update_beagle_partials(t_tree* tree, t_edge* b, t_node* d);
void update_beagle_ras(t_mod* mod);
void update_beagle_efrqs(t_mod* mod);
void update_beagle_eigen(t_mod* mod);
void calc_edgelks_beagle(t_edge* b, t_tree* tree);
double* int_to_double(const int* src, int num_elems);
double* short_to_double(const short* src, int num_elems);
double* float_to_double(const phydbl *src, int num_elems);



#endif // BEAGLE_UTILS_H
