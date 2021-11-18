/*

PHYML :  a program that  computes maximum likelihood  phylogenies from
DNA or AA homologous sequences 

Copyright (C) Stephane Guindon. Oct 2003 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/

#include <config.h>

#ifndef FREE_H
#define FREE_H

#include "utilities.h"

void Free_All_Nodes_Light(t_tree *tree);
void Free_All_Edges_Light(t_tree *tree);
void Free_Mat(matrix *mat);
void Free_Partial_Lk(phydbl *p_lk,int len,int n_catg);
void Free_Tree(t_tree *tree);
void Free_Bip(t_tree *tree);
void Free_Edge(t_edge *b);
void Free_Node(t_node *n);
void Free_Calign(calign *data);
void Free_Seq(align **d,int n_otu);
void Free_All(align **d,calign *cdata,t_tree *tree);
void Free_SubTree(t_edge *b_fcus,t_node *a,t_node *d,t_tree *tree);
void Free_Tree_Ins_Tar(t_tree *tree);
void Free_Tree_Pars(t_tree *tree);
void Free_Edge_Pars(t_edge *b);
void Free_Edge_Pars_Left(t_edge *b);
void Free_Edge_Pars_Rght(t_edge *b);
void Free_Tree_Lk(t_tree *tree);
void Free_Node_Lk(t_node *n);
void Free_Edge_Lk(t_edge *b);
void Free_Model_Complete(t_mod *mod);
void Free_Model_Basic(t_mod *mod);
void Free_Custom_Model(t_mod *mod);
void Free_Efrq(t_efrq *e_frq);
void Free_Rmat(t_rmat *r_mat);
void Free_Model(t_mod *mod);
void Free(void *p);
void Free_Input(option *io);
void Free_Tree_List(t_treelist *list);
void Free_St(supert_tree *st);
void Free_Eigen(eigen *eigen_struct);
void Free_One_Spr(t_spr *this_spr);
void Free_Spr_List_One_Edge(t_tree *tree);
void Free_Spr_List_All_Edge(t_tree *mixt_tree);
void Free_Actual_CSeq(calign *data);
void Free_Prefix_Tree(pnode *n,int size);
void Free_Pnode(pnode *n);
void Free_Optimiz(t_opt *s_opt);
void Free_Nexus(option *io);
void Free_Nexus_Com(nexcom **com);
void Free_Nexus_Parm(nexparm *parm);
void Free_RAS(t_ras *ras);
void XML_Free_XML_Tree(xml_node *node);
void XML_Free_XML_Node(xml_node *node);
void XML_Free_XML_Attr(xml_attr *attr);
void XML_Free_XML_Ds(t_ds *ds);
void Free_String(t_string *ts);
void Free_Vect_Dbl(vect_dbl *v);
void Free_Vect_Int(vect_int *v);
void Free_Scalar_Dbl(scalar_dbl *v);
void Free_Scalar_Int(scalar_int *v);
void Free_Edge_Core(t_edge *b);
void Free_Rates(t_rate *rates);
void Free_Calib(t_cal *cal);
void Free_Edge_Lk_Left(t_edge *b);
void Free_Edge_Lk_Rght(t_edge *b);
void Free_Geo(t_geo *t);
void Free_Geo_Coord(t_geo_coord *t);
void PHYREX_Free_Disk(t_dsk *t);
void PHYREX_Free_Ldisk(t_ldsk *t);
void PHYREX_Free_Ldsk_Struct(t_tree *tree);
void Free_Poly(t_poly *p);
void Free_Mmod(t_phyrex_mod *mmod);
void Free_Efrq_Weights(t_mod *mixt_mod);
void Free_Rmat_Weights(t_mod *mixt_mod);
void JSON_Free_KeyVal(json_kv *kv);
void JSON_Free_Object(json_o *o);
void JSON_Free_Array(json_a *a);
void Free_Edge_Loc_Rght(t_edge *b);
void Free_Edge_Loc_Left(t_edge *b);
void Free_Edge_Loc(t_edge *b);
void Free_NNI(t_nni *t);
void Free_Linked_List(t_ll *t);
void Free_TBE_Matrices(int n_otu,  short unsigned*** i_matrix, short unsigned*** c_matrix,
		       short unsigned*** hamming, short unsigned** min_dist,
		       short unsigned**  min_dist_edge, int** cluster_sizes);
void Free_Extra_Edge_Lk(t_tree *tree);
void Free_Edge_Length(t_edge *b);
void Free_Clade(t_clad *this);
void Free_Label(t_label *lab);
void Free_Contrasts(t_ctrst *ctrst);
void TIMES_Free_Times(t_time *times);
void Free_Best_Spr(t_tree *tree);
void Free_All_Edges_Light(t_tree *tree);
void Free_All_Edges_Lens(t_tree *tree);

#endif
