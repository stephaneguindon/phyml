/*

PHYML :  a program that  computes maximum likelihood  phylogenies from
DNA or AA homologous sequences

Copyright (C) Stephane Guindon. Oct 2003 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/

#include <config.h>

#ifndef SSE_H
#define SSE_H

#include "utilities.h"
#include "optimiz.h"
#include "models.h"
#include "free.h"
#include "times.h"
#include "mixt.h"

#if (defined(__SSE__) || defined(__SSE2__) || defined(__SSE3__))

void SSE_Update_Partial_Lk(t_tree *tree,t_edge *b_fcus,t_node *n);
void SSE_Update_Eigen_Lr(t_edge *b, t_tree *tree);
phydbl SSE_Lk_Core_One_Class_Eigen_Lr(phydbl *dot_prod, phydbl *expl, int ns);
phydbl SSE_Lk_Core_One_Class_No_Eigen_Lr(phydbl *p_lk_left, phydbl *p_lk_rght, phydbl *Pij, phydbl *tPij, phydbl *pi, int ns, int ambiguity_check, int state);
phydbl SSE_Lk_Core(int state, int ambiguity_check, t_edge *b, t_tree *tree);
phydbl SSE_Lk_Core_Nucl(int state, int ambiguity_check, t_edge *b, t_tree *tree);
phydbl SSE_Lk_Core_AA(int state, int ambiguity_check, t_edge *b, t_tree *tree);
void SSE_Lk_dLk_Core_One_Class_Eigen_Lr(phydbl *dot_prod, phydbl *expl, unsigned int ns, phydbl *lk, phydbl *dlk);

void SSE_Partial_Lk_Exex(const __m128d *_tPij1, const int state1, const __m128d *_tPij2, const int state2, const int ns, __m128d *plk0);
void SSE_Partial_Lk_Exin(const __m128d *_tPij1, const int state1, const __m128d *_tPij2, const phydbl *_plk2, __m128d *_pmat2plk2, const int ns, __m128d *_plk0);
void SSE_Partial_Lk_Inin(const __m128d *_tPij1, const phydbl *plk1, __m128d *_pmat1plk1, const __m128d *_tPij2, const phydbl *plk2, __m128d *_pmat2plk2, const int ns, __m128d *_plk0);
void SSE_Matrix_Vect_Prod(const __m128d *_m_transpose, const phydbl *_v, const int ns, __m128d *_u);

#endif
#endif
