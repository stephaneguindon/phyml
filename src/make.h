/*

PHYML :  a program that  computes maximum likelihood  phylogenies from
DNA or AA homologous sequences 

Copyright (C) Stephane Guindon. Oct 2003 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/

#include <config.h>

#ifndef MAKE_H
#define MAKE_H

#include "utilities.h"

void Make_All_Edges_Lk(t_node *a,t_node *d,t_tree *tree);
void Make_New_Edge_Label(t_edge *b);
t_edge *Make_Edge_Light(t_node *a,t_node *d,int num);
void Make_Edge_Pars(t_edge *b,t_tree *tree);
void Make_Edge_Pars_Left(t_edge *b,t_tree *tree);
void Make_Edge_Pars_Rght(t_edge *b,t_tree *tree);
void Make_Edge_Lk(t_edge *b,t_tree *tree);
void Make_Edge_Lk_Left(t_edge *b,t_tree *tree);
void Make_Edge_Lk_Rght(t_edge *b,t_tree *tree);
void Make_Edge_NNI(t_edge *b);
nni *Make_NNI();
t_node *Make_Node_Light(int num);
void Make_Node_Lk(t_node *n);
nexcom **Make_Nexus_Com();
nexparm *Make_Nexus_Parm();
matrix *Make_Mat(int n_otu);
t_tree *Make_Tree_From_Scratch(int n_otu,calign *data);
t_tree *Make_Tree(int n_otu);
void Make_Tree_Path(t_tree *tree);
void Make_All_Tree_Nodes(t_tree *tree);
void Make_All_Tree_Edges(t_tree *tree);
calign *Make_Cseq(int n_otu,int crunch_len,int state_len,int init_len,char **sp_names);
t_treelist *Make_Treelist(int list_size);
t_opt *Make_Optimiz();
void Make_Custom_Model(t_mod *mod);
t_mod *Make_Model_Basic();
void Make_Model_Complete(t_mod *mod);
t_efrq *Make_Efrq(int ns);
t_rmat *Make_Rmat(int ns);
option *Make_Input();
eigen *Make_Eigen_Struct(int ns);
triplet *Make_Triplet_Struct(t_mod *mod);
void Make_Short_L(t_tree *tree);
void Make_RAS_Complete(t_ras *ras);
t_ras *Make_RAS_Basic();
void Make_Best_Spr(t_tree *tree);
void Make_Spr_List(t_tree *tree);
t_spr *Make_One_Spr(t_tree *tree);
void Make_Tree_4_Pars(t_tree *tree, calign *cdata, int n_site);
t_string *Make_String(int len);
t_mcmc *MCMC_Make_MCMC_Struct();
void Make_Tree_4_Lk(t_tree *tree,calign *cdata,int n_site);
t_rate *RATES_Make_Rate_Struct(int n_otu);
t_cal *Make_Calib();
void Make_Efrq_Weight(t_tree *mixt_tree);
void Make_Rmat_Weight(t_tree *mixt_tree);
t_geo *GEO_Make_Geo_Basic();
void GEO_Make_Geo_Complete(int ldscape_sz,int n_dim,int n_tax,t_geo *t);
t_geo_coord *GEO_Make_Geo_Coord(int n_dim);
t_phyrex_mod *PHYREX_Make_Migrep_Model();
t_dsk *PHYREX_Make_Disk_Event(int n_dim, int n_otu);
t_ldsk *PHYREX_Make_Lindisk_Node(int n_dim);
void PHYREX_Make_Lindisk_Next(t_ldsk *t);
t_poly *Make_Poly(int n);

#endif
