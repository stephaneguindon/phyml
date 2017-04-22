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
void Default_Update_Partial_Lk(t_tree *tree,t_edge *b_fcus,t_node *n);
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
void Print_Lk_Given_Edge_Recurr(t_node *a, t_node *d, t_edge *b, t_tree *tree);
phydbl *Post_Prob_Rates_At_Given_Edge(t_edge *b, phydbl *post_prob, t_tree *tree);
phydbl Lk_With_MAP_Branch_Rates(t_tree *tree);
void Init_Tips_At_One_Site_Generic_Int(char *state, int ns, int state_len, int pos, short int *p_pars);
void Init_Tips_At_One_Site_Generic_Float(char *state, int ns, int state_len, int pos, phydbl *p_lk);
void Alias_Subpatt(t_tree *tree);
void Alias_One_Subpatt(t_node *a, t_node *d, t_tree *tree);
void Alias_Subpatt_Post(t_node *a, t_node *d, t_tree *tree);
void Alias_Subpatt_Pre(t_node *a, t_node *d, t_tree *tree);
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
void Sample_Ancestral_Seq(int mutmap, int fromprior, t_tree *tree);
void Map_Mutations(t_node *a, t_node *d, int sa, int sd, t_edge *b, int site, int rate_cat, int *muttype, phydbl *muttime, int *n_mut, t_tree *tree);
void Sample_Ancestral_Seq_Pre(t_node *a, t_node *d, t_edge *b,int site, int rate_cat,int *muttype, phydbl *muttime, int *n_mut,int mutmap, int fromprior, t_tree *tree);
phydbl Wrap_Lk_Times(t_edge *b, t_tree *tree, supert_tree *stree);
phydbl Lk_LastFirst(t_tree *tree);
int Check_Lk_At_Given_Edge(int verbose, t_tree *tree);
void Ancestral_Sequences_One_Node(t_node *mixt_d, t_tree *mixt_tree, int print);
void Ancestral_Sequences(t_tree *tree, int print);
void Stepwise_Add_Lk(t_tree *tree);
void Update_Eigen_Lr(t_edge *b, t_tree *tree);
phydbl dLk(phydbl *l, t_edge *b, t_tree *tree);
phydbl Lk_Core(int state, int ambiguity_check, short int derivative,
                             phydbl *p_lk_left, phydbl *p_lk_rght,
                             phydbl *Pij_rr,
                             t_edge *b,
                             t_tree *tree);
phydbl Lk_Core_Eigen_Lr(phydbl *expl, phydbl *dot_prod, short int derivative, t_edge *b, t_tree *tree);
phydbl Invariant_Lk(int fact_sum_scale, int site, int *num_prec_issue, t_tree *tree);


#if defined(__AVX__)
phydbl AVX_Lk_Core(int state, int ambiguity_check, t_edge *b, t_tree *tree);
phydbl AVX_Lk_Core_Nucl(int state, int ambiguity_check, t_edge *b, t_tree *tree);
phydbl AVX_Lk_Core_AA(int state, int ambiguity_check, t_edge *b, t_tree *tree);
#elif defined(__SSE3__)
phydbl SSE_Lk_Core(int state, int ambiguity_check, t_edge *b, t_tree *tree);
phydbl SSE_Lk_Core_Nucl(int state, int ambiguity_check, t_edge *b, t_tree *tree);
phydbl SSE_Lk_Core_AA(int state, int ambiguity_check, t_edge *b, t_tree *tree);
#endif


#endif






