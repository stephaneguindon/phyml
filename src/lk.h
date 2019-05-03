/*

PHYML :  a program that  computes maximum likelihood  phylogenies from
DNA or AA homologous sequences

Copyright (C) Stephane Guindon. Oct 2003 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/

#include <config.h>

#ifndef LK_H
#define LK_H

#include "utilities.h"
#include "optimiz.h"
#include "models.h"
#include "free.h"
#include "times.h"
#include "mixt.h"
#include "avx.h"
#include "sse.h"

void Update_All_Partial_Lk(t_tree *tree);
void Init_Tips_At_One_Site_Nucleotides_Float(char state, int pos, phydbl *p_lk);
void Init_Tips_At_One_Site_AA_Float(char aa, int pos, phydbl *p_lk);
void Get_All_Partial_Lk(t_tree *tree,t_edge *b_fcus,t_node *a,t_node *d);
void Get_All_Partial_Lk_Scale(t_tree *tree,t_edge *b_fcus,t_node *a,t_node *d);
void Post_Order_Lk(t_node *pere, t_node *fils, t_tree *tree);
void Pre_Order_Lk(t_node *pere, t_node *fils, t_tree *tree);
phydbl Lk(t_edge *b, t_tree *tree);
void Site_Lk(t_tree *tree);
phydbl Return_Abs_Lk(t_tree *tree);
matrix *ML_Dist(calign *data, t_mod *mod);
phydbl Lk_Given_Two_Seq(calign *data, int numseq1, int numseq2, phydbl dist, t_mod *mod, phydbl *loglk);
void Unconstraint_Lk(t_tree *tree);
void Update_Partial_Lk(t_tree *tree,t_edge *b_fcus,t_node *n);
void Update_Partial_Lk_Generic(t_tree *tree,t_edge *b_fcus,t_node *n);
/* void Default_Update_Partial_Lk(t_tree *tree,t_edge *b_fcus,t_node *n); */
void Init_Partial_Lk_Tips_Double(t_tree *tree);
void Init_Partial_Lk_Tips_Int(t_tree *tree);
void Init_Partial_Lk_At_One_Node(t_node *a, t_tree *tree);
void Update_PMat_At_Given_Edge(t_edge *b_fcus, t_tree *tree);
void Sort_Sites_Based_On_Lk(t_tree *tree);
void Get_Partial_Lk_Scale(t_tree *tree, t_edge *b_fcus, t_node *a, t_node *d);
void Get_Partial_Lk(t_tree *tree, t_edge *b_fcus, t_node *a, t_node *d);
void Init_Tips_At_One_Site_Nucleotides_Int(char state, int pos, short int *p_pars);
void Init_Tips_At_One_Site_AA_Int(char aa, int pos, short int *p_pars);
void Update_Partial_Lk_Along_A_Path(t_node **path, int path_length, t_tree *tree);
phydbl Lk_Dist(phydbl *F, phydbl dist, t_mod *mod);
phydbl Update_Lk_At_Given_Edge(t_edge *b_fcus, t_tree *tree);
void Update_Partial_Lk_Greedy(t_tree *tree, t_edge *b_fcus, t_node *n);
void Get_All_Partial_Lk_Scale_Greedy(t_tree *tree, t_edge *b_fcus, t_node *a, t_node *d);
phydbl Lk_Triplet(t_node *a, t_node *d, t_tree *tree);
phydbl *Post_Prob_Rates_At_Given_Edge(t_edge *b, phydbl *post_prob, t_tree *tree);
phydbl Lk_With_MAP_Branch_Rates(t_tree *tree);
void Init_Tips_At_One_Site_Generic_Int(char *state, int ns, int state_len, int pos, short int *p_pars);
void Init_Tips_At_One_Site_Generic_Float(char *state, int ns, int state_len, int pos, phydbl *p_lk);
void Backup_Partial_Lk(t_node *d, t_edge *b, t_tree *tree);
void Restore_Partial_Lk(t_node *d, t_edge *b, t_tree *tree);
void Backup_Partial_Scale(t_node *d, t_edge *b, t_tree *tree);
void Restore_Partial_Scale(t_node *d, t_edge *b, t_tree *tree);
void Init_Partial_Lk_Loc(t_tree *tree);
phydbl Lk_Normal_Approx(t_tree *tree);
phydbl Wrap_Lk(t_edge *b, t_tree *tree, supert_tree *stree);
phydbl Wrap_Lk_At_Given_Edge(t_edge *b, t_tree *tree, supert_tree *stree);
phydbl Wrap_Part_Lk_At_Given_Edge(t_edge *b, t_tree *tree, supert_tree *stree);
phydbl Wrap_Part_Lk(t_edge *b, t_tree *tree, supert_tree *stree);
phydbl Wrap_Geo_Lk(t_edge *b, t_tree *tree, supert_tree *stree);
phydbl Wrap_Diff_Lk_Norm_At_Given_Edge(t_edge *b, t_tree *tree, supert_tree *stree);
phydbl Wrap_Lk_Rates(t_edge *b, t_tree *tree, supert_tree *stree);
phydbl Wrap_Lk_Linreg(t_edge *b, t_tree *tree, supert_tree *stree);
phydbl Wrap_Lk_Times(t_edge *b, t_tree *tree, supert_tree *stree);
phydbl Lk_LastFirst(t_tree *tree);
int Check_Lk_At_Given_Edge(int verbose, t_tree *tree);
void Stepwise_Add_Lk(t_tree *tree);
void Update_Eigen_Lr(t_edge *b, t_tree *tree);
phydbl dLk(phydbl *l, t_edge *b, t_tree *tree);
phydbl Lk_Core(int state, int ambiguity_check,phydbl *p_lk_left, phydbl *p_lk_rght,phydbl *Pij_rr,phydbl *tPij_rr,t_edge *b,t_tree *tree);
/* void Lk_dLk_Core_Eigen_Lr(phydbl *expl, phydbl *dot_prod, t_edge *b, phydbl *lk, phydbl *dlk, t_tree *tree); */
phydbl Lk_Core_Eigen_Lr(phydbl *expl, phydbl *dot_prod, t_edge *b, t_tree *tree);
phydbl Invariant_Lk(int fact_sum_scale, int site, int *num_prec_issue, t_tree *tree);
void Composite_Lk(t_tree *tree);
int Sample_Ancestral_Seq_Core(t_node *a, t_node *d, t_edge *b, int r_cat, int site, t_tree *tree);
void Default_Update_Partial_Lk(t_tree *tree, t_edge *b, t_node *d);
void Lk_dLk_Core_Eigen_Lr(phydbl *expl, phydbl *dot_prod, t_edge *b, phydbl *lk, phydbl *dlk, t_tree *tree);
void Lk_dLk_Core_One_Class_Eigen_Lr(phydbl *dot_prod, phydbl *expl, unsigned int ns, phydbl *lk, phydbl *dlk);
phydbl Lk_Core_One_Class_Eigen_Lr(phydbl *dot_prod, phydbl *expl, int ns);
void Core_Default_Update_Partial_Lk(const t_node *n_v1, const t_node *n_v2,
                                    phydbl *plk0, const phydbl *plk1, const phydbl *plk2,
                                    const phydbl *Pij1, const phydbl *Pij2,
                                    int *sum_scale0, const int *sum_scale1, const int *sum_scale2,
                                    const int ns, const int ncatg, const int npatterns, const int apply_scaling,
                                    const phydbl *wght);






/* phydbl Lk_Core_One_Class_Eigen_Lr(phydbl *dot_prod, phydbl *expl, int ns); */
/* void Lk_dLk_Core_One_Class_Eigen_Lr(phydbl *dot_prod, phydbl *expl, unsigned int ns, phydbl *lk, phydbl *dlk); */
phydbl Lk_Core_One_Class_No_Eigen_Lr(phydbl *p_lk_left, phydbl *p_lk_rght, phydbl *Pij,phydbl *pi, int ns, int ambiguity_check, int state);

void Pull_Scaling_Factors(int site,t_edge *b,t_tree *tree);
void Set_All_Partial_Lk(t_node **n_v1, t_node **n_v2,phydbl **p_lk, int **sum_scale, int **p_lk_loc,phydbl **Pij1, phydbl **tPij1, phydbl **p_lk_v1, int **sum_scale_v1,phydbl **Pij2, phydbl **tPij2, phydbl **p_lk_v2, int **sum_scale_v2,t_node *d, t_edge *b, t_tree *tree
#ifdef BEAGLE
                                , int *dest_p_idx, int *child1_p_idx, int* child2_p_idx, int* Pij1_idx, int* Pij2_idx
#endif
                                );
void Set_Partial_Lk_One_Side(phydbl **Pij, phydbl **tPij, phydbl **p_lk,  int **sum_scale, t_node *d, t_edge *b, t_tree *tree
#ifdef BEAGLE
                                     , int* child_p_idx, int* Pij_idx
#endif
                                     );
void Rate_Correction(int exponent, phydbl *site_lk_cat);
void Switch_Partial_Lk_Pre(t_node *a, t_node *d, t_edge *b, short int yesno, t_tree *tree);
void Switch_Partial_Lk_Post(t_node *a, t_node *d, t_edge *b, short int yesno, t_tree *tree);

void Partial_Lk_Inin(const phydbl *Pij1, const phydbl *plk1, const phydbl *Pij2, const phydbl *plk2, const int ns, phydbl *plk0);
void Partial_Lk_Exex(const phydbl *Pij1, const int state1, const phydbl *Pij2, const int state2, const int ns, phydbl *plk0);
void Partial_Lk_Exin(const phydbl *Pij1, const int state1, const phydbl *Pij2, const phydbl *plk2, const int ns, phydbl *plk0);


#endif






