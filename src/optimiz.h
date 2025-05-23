/*

PHYML :  a program that  computes maximum likelihood  phylogenies from
DNA or AA homologous sequences

Copyright (C) Stephane Guindon. Oct 2003 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/

#include <config.h>

#ifndef OPTIMIZ_H
#define OPTIMIZ_H

#include "free.h"
#include "lk.h"
#include "mg.h"
#include "models.h"
#include "tiporder.h"
#include "utilities.h"

void   Optimize_Br_Len_Serie(int n_max_iter, t_tree *tree);
void   Optimiz_Ext_Br(t_tree *tree);
void   Optimize_Param_Parall(t_tree *tree);
phydbl Optimize_Branch_Quad(t_tree *tree, calign *cdata, t_edge *b_fcus);
void   Optimize_After_Hide(t_tree *tree, calign *cdata, t_node *h);
void   Round_Optimize(t_tree *tree, int n_round_max);
int    Dist_Seq_Brak(phydbl *ax, phydbl *bx, phydbl *cx, phydbl *fa, phydbl *fb,
                     phydbl *fc, calign *data, int num1, int num2, t_mod *mod);
phydbl Dist_Seq_Brent(phydbl ax, phydbl bx, phydbl cx, phydbl tol, phydbl *xmin,
                      calign *data, int num1, int num2, t_mod *mod);
phydbl Kappa_Golden(phydbl ax, phydbl bx, phydbl cx, phydbl tol, phydbl *xmin,
                    t_tree *tree, calign *cdata);
phydbl Lambda_Golden(phydbl ax, phydbl bx, phydbl cx, phydbl tol, phydbl *xmin,
                     t_tree *tree, calign *cdata);
phydbl Alpha_Golden_Br_Opt(phydbl ax, phydbl bx, phydbl cx, phydbl tol,
                           phydbl *xmin, t_tree *tree, calign *cdata, int n_opt,
                           phydbl *init_l);
phydbl Alpha_Golden(phydbl ax, phydbl bx, phydbl cx, phydbl tol, phydbl *xmin,
                    t_tree *tree, calign *cdata);
phydbl Br_Len_Golden(phydbl ax, phydbl bx, phydbl cx, phydbl tol, phydbl *xmin,
                     t_edge *b_fcus, t_tree *tree);
phydbl Br_Len_Opt(phydbl *l, t_edge *b, t_tree *tree);
int    Br_Len_Brak(phydbl *ax, phydbl *bx, phydbl *cx, phydbl *fa, phydbl *fb,
                   phydbl *fc, t_edge *b_fcus, t_tree *tree);
phydbl Optimize_Path_Length(t_mod *mod, calign *cdata, t_edge *a, int lra,
                            t_edge *b, int lrb, phydbl i_len);
void   Optimize_Param_Serie(t_node *a, t_node *d, t_edge *b_fcus, t_tree *tree,
                            calign *cdata, int n_passes);
phydbl Optimize_Dist(t_mod *mod, phydbl init, calign *twoseqs);
phydbl Pinvar_Golden(phydbl ax, phydbl bx, phydbl cx, phydbl tol, phydbl *xmin,
                     t_tree *tree, calign *cdata, int n_iter_max);
void   Optimize_Pinvar(t_tree *tree);
int    Lambda_Brak(phydbl *ax, phydbl *bx, phydbl *cx, phydbl *fa, phydbl *fb,
                   phydbl *fc, t_tree *tree);
int    Kappa_Brak(phydbl *ax, phydbl *bx, phydbl *cx, phydbl *fa, phydbl *fb,
                  phydbl *fc, t_tree *tree);
int    Alpha_Brak(phydbl *ax, phydbl *bx, phydbl *cx, phydbl *fa, phydbl *fb,
                  phydbl *fc, t_tree *tree);
int    Pinvar_Brak(phydbl *ax, phydbl *bx, phydbl *cx, phydbl *fa, phydbl *fb,
                   phydbl *fc, t_tree *tree);
void   Optimiz_All_Free_Param(t_tree *tree, int verbose);
void   Optimiz_RRparam_GTR(t_tree *tree, int num_param);
phydbl RRparam_GTR_Golden(phydbl ax, phydbl bx, phydbl cx, phydbl tol,
                          phydbl *xmin, t_tree *tree, calign *cdata,
                          phydbl *param, int n_iter_max);

int    Powell_GTR_Param(t_tree *tree, phydbl *p, int n, phydbl ftol);
phydbl Linmin_GTR_Param(t_tree *tree, phydbl *p, phydbl *xi, int n);
phydbl F1dim(t_tree *tree, phydbl x, phydbl *p, phydbl *xi, phydbl n);
int    Mnbrak_1dim(phydbl *ax, phydbl *bx, phydbl *cx, phydbl *fa, phydbl *fb,
                   phydbl *fc, t_tree *tree, phydbl *p, phydbl *xi, phydbl n);
phydbl Brent_1dim(phydbl ax, phydbl bx, phydbl cx, phydbl tol, phydbl *xmin,
                  t_tree *tree, phydbl *p, phydbl *xi, phydbl n);

int  Min_With_Derivatives(t_tree *tree, phydbl *p, int n, phydbl ftol,
                          phydbl step_size, phydbl (*func)(void),
                          void (*dfunc)(void), phydbl (*linmin)(void));
void BFGS(t_tree *tree, phydbl *p, int n, phydbl gtol, phydbl difff,
          phydbl step_size, short int logt, short int expt,
          short int is_positive, phydbl (*func)(t_tree *tree),
          int (*dfunc)(t_tree *tree, phydbl *param, int n_param,
                       phydbl stepsize, short int logt, short int expt,
                       phydbl (*func)(t_tree *tree), phydbl *derivatives,
                       int is_positive),
          int (*lnsrch)(t_tree *tree, int n, phydbl *xold, phydbl fold,
                        phydbl *g, phydbl *p, phydbl *x, phydbl *f,
                        phydbl stpmax, int *check, short int logt,
                        short int expt, int is_positive),
          int *failed);

void BFGS_Nonaligned(
    t_tree *tree, phydbl **p, int n, phydbl gtol, phydbl difff,
    phydbl step_size, short int logt, short int expt, short int is_positive,
    phydbl (*func)(t_tree *tree),
    int (*dfunc_nonaligned)(t_tree *tree, phydbl **param, int n_param,
                            phydbl stepsize, short int logt, short int expt,
                            phydbl (*func)(t_tree *tree), phydbl *derivatives,
                            int is_positive),
    int (*lnsrch_nonaligned)(t_tree *tree, int n, phydbl **xold, phydbl fold,
                             phydbl *g, phydbl *p, phydbl *x, phydbl *f,
                             phydbl stpmax, int *check, short int logt,
                             short int expt, int is_positive),
    int *failed);

void Optimize_Single_Param_Generic(t_tree *tree, phydbl *param, phydbl lim_inf,
                                   phydbl lim_sup, phydbl tol, int n_max_iter,
                                   int quickdirty);
int  Generic_Brak_Lk(phydbl *param, phydbl *ax, phydbl *bx, phydbl *cx,
                     phydbl *fa, phydbl *fb, phydbl *fc, phydbl min, phydbl max,
                     phydbl (*obj_func)(t_edge *, t_tree *, supert_tree *),
                     t_edge *branch, t_tree *tree, supert_tree *stree);
int  Generic_Brak(phydbl *param, phydbl *ax, phydbl *bx, phydbl *cx, phydbl *fa,
                  phydbl *fb, phydbl *fc, phydbl lim_inf, phydbl lim_sup,
                  t_tree *tree);
void Optimize_Br_Len_Serie_Post(t_node *a, t_node *d, t_edge *b_fcus,
                                t_tree *tree);
void Optimize_Global_Rate(t_tree *tree);
phydbl Fast_Br_Len(t_edge *b, t_tree *tree, int approx);

void   EM_Dist(t_mod *mod, calign *data);
phydbl Dist_F_Brent(phydbl ax, phydbl bx, phydbl cx, phydbl tol, int n_iter_max,
                    phydbl *param, phydbl *F, t_mod *mod);
int    Dist_F_Brak(phydbl *ax, phydbl *bx, phydbl *cx, phydbl *F, phydbl *param,
                   t_mod *mod);
void   Opt_Dist_F(phydbl *dist, phydbl *F, t_mod *mod);
phydbl Missing_Dist_Brent(phydbl ax, phydbl bx, phydbl cx, phydbl tol,
                          int n_iter_max, int x, int y, matrix *mat);
int    Missing_Dist_Brak(phydbl *ax, phydbl *bx, phydbl *cx, int x, int y,
                         matrix *mat);
void   Opt_Missing_Dist(int x, int y, matrix *mat);
int    Optimiz_Alpha_And_Pinv(t_tree *tree, int verbose);
phydbl Node_Time_Brent(phydbl ax, phydbl bx, phydbl cx, phydbl tol, t_node *anc,
                       t_node *des, t_tree *tree, int n_iter_max);
phydbl Time_Stamps_Mult_Brent(phydbl ax, phydbl bx, phydbl cx, phydbl tol,
                              t_tree *tree, int n_iter_max);
phydbl Branch_Rate_Shape_Brent(phydbl ax, phydbl bx, phydbl cx, phydbl tol,
                               phydbl *xmin, t_tree *tree, int n_iter_max);
phydbl Node_Time_Brent_Fixed_Br_Len(phydbl ax, phydbl bx, phydbl cx, phydbl tol,
                                    t_node *n, t_tree *tree, int n_iter_max);

int  Generic_Brent_Lk(phydbl *param, phydbl ax, phydbl cx, phydbl tol,
                      int n_iter_max, int quickdirty,
                      phydbl (*obj_func)(t_edge *, t_tree *, supert_tree *),
                      t_edge *branch, t_tree *tree, supert_tree *stree,
                      short int logt, short int expt);
void Round_Optimize_Node_Heights(t_tree *tree);
void Opt_Node_Heights_Recurr_Pre(t_node *a, t_node *d, t_tree *tree);
void Opt_Node_Heights_Recurr(t_tree *tree);

int Lnsrch(t_tree *tree, int n, phydbl *xold, phydbl fold, phydbl *g, phydbl *p,
           phydbl *x, phydbl *f, phydbl stpmax, int *check, short int logt,
           short int expt, int is_positive);

int Lnsrch_Nonaligned(t_tree *tree, int n, phydbl **xold, phydbl fold,
                      phydbl *g, phydbl *p, phydbl *x, phydbl *f, phydbl stpmax,
                      int *check, short int logt, short int expt,
                      int is_positive);

void   Optimize_RR_Params(t_tree *mixt_tree, int verbose);
void   Optimize_TsTv(t_tree *mixt_tree, int verbose);
void   Optimize_Lambda(t_tree *mixt_tree, int verbose);
void   Optimize_Alpha(t_tree *mixt_tree, int verbose);
void   Optimize_Pinv(t_tree *mixt_tree, int verbose);
void   Optimize_State_Freqs(t_tree *mixt_tree, int verbose);
void   Optimize_Rmat_Weights(t_tree *mixt_tree, int verbose);
void   Optimize_Efrq_Weights(t_tree *mixt_tree, int verbose);
void   Optimize_Free_Rate(t_tree *mixt_tree, int verbose);
void   Optimize_Free_Rate_Weights(t_tree *tree, int fast, int verbose);
void   Optimize_Free_Rate_Rr(t_tree *tree, int fast, int verbose);
void   Optimize_Br_Len_Multiplier(t_tree *tree, int verbose);
phydbl Generic_Brent(phydbl *param, phydbl ax, phydbl cx, phydbl tol,
                     int     n_iter_max, phydbl (*obj_func)(t_tree *),
                     t_tree *tree);
void   Least_Square_Node_Ages(t_tree *tree);
void   Optimize_Lvar(t_tree *mixt_tree, int verbose);
void   Optimize_M4mod(t_tree *tree, int verbose);

#endif
