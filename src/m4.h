/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include <config.h>

#ifndef M4_H
#define M4_H

#include "spr.h"
#include "utilities.h"
#include "lk.h"
#include "optimiz.h"
#include "bionj.h"
#include "models.h"
#include "free.h"
#include "help.h"
#include "simu.h"
#include "eigen.h"
#include "pars.h"
#include "alrt.h"
#include "time.h"
#include "draw.h"
#ifdef PART
#include "mg.h"
#endif

int M4_main(int argc, char **argv);
void M4_Make_Complete(int n_h, int n_o, m4 *m4mod);
m4 *M4_Make_Light();
void M4_Free_M4_Model(m4 *m4mod);
void M4_Init_Qmat(m4 *m4mod, calign *data, t_mod *mod);
void M4_Update_Qmat(m4 *m4mod, t_mod *mod);
void M4_Init_Model(m4 *m4mod, calign *data, t_mod *mod);
void M4_Init_P_Lk_Tips_Double(t_tree *tree);
void M4_Init_P_Lk_Tips_Int(t_tree *tree);
void M4_Post_Prob_H_Class_Edge_Site(t_edge *b, phydbl ****integral, phydbl *postprob, t_tree *tree);
phydbl ****M4_Integral_Term_On_One_Edge(t_edge *b, t_tree *tree);
phydbl ***M4_Compute_Proba_Hidden_States_On_Edges(t_tree *tree);
void M4_Free_Integral_Term_On_One_Edge(phydbl ****integral, t_tree *tree);
void M4_Compute_Posterior_Mean_Rates(phydbl ***post_probs, t_tree *tree);
void M4_Scale_Br_Len(t_tree *tree);
void M4_Detect_Site_Switches_Experiment(t_tree *tree);
m4 *M4_Copy_M4_Model(t_mod *ori_mod, m4 *ori_m4mod);
void M4_Posterior_Prediction_Experiment(t_tree *tree);
void M4_Set_M4mod_Default(m4 *m4mod);
phydbl **M4_Site_Branch_Classification(phydbl ***post_probs, t_tree *tree);
void M4_Site_Branch_Classification_Experiment(t_tree *tree);

#endif
