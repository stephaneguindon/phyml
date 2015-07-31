/*

PHYML :  a program that  computes maximum likelihood  phylogenies from
DNA or AA homologous sequences 

Copyright (C) Stephane Guindon. Oct 2003 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/

#include <config.h>

#ifndef INIT_H
#define INIT_H

#include "utilities.h"

void Init_Eigen_Struct(eigen *this);
void Init_Scalar_Dbl(scalar_dbl *p);
void Init_Scalar_Int(scalar_int *p);
void Init_Vect_Dbl(int len,vect_dbl *p);
void Init_Vect_Int(int len,vect_int *p);
void Init_Tree(t_tree *tree,int n_otu);
void Init_Edge_Light(t_edge *b,int num);
void Init_Node_Light(t_node *n,int num);
void Init_NNI(nni *a_nni);
void Init_Nexus_Format(nexcom **com);
void Init_Mat(matrix *mat,calign *data);
void Set_Defaults_Input(option *io);
void Set_Defaults_Model(t_mod *mod);
void Set_Defaults_Optimiz(t_opt *s_opt);
void XML_Init_Node(xml_node *prev,xml_node *new_node,char *name);
void Init_One_Spr(t_spr *a_spr);
void Init_Model(calign *data,t_mod *mod,option *io);
int Init_Qmat_Dayhoff(phydbl *daa,phydbl *pi);
int Init_Qmat_DCMut(phydbl *daa,phydbl *pi);
int Init_Qmat_MtArt(phydbl *daa,phydbl *pi);
int Init_Qmat_HIVb(phydbl *daa,phydbl *pi);
int Init_Qmat_HIVw(phydbl *daa,phydbl *pi);
int Init_Qmat_JTT(phydbl *daa,phydbl *pi);
int Init_Qmat_MtREV(phydbl *daa,phydbl *pi);
int Init_Qmat_LG(phydbl *daa,phydbl *pi);
int Init_Qmat_WAG(phydbl *daa,phydbl *pi);
int Init_Qmat_RtREV(phydbl *daa,phydbl *pi);
int Init_Qmat_CpREV(phydbl *daa,phydbl *pi);
int Init_Qmat_VT(phydbl *daa,phydbl *pi);
int Init_Qmat_Blosum62(phydbl *daa,phydbl *pi);
int Init_Qmat_MtMam(phydbl *daa,phydbl *pi);
int Init_Qmat_AB(phydbl *daa, phydbl *pi);
void XML_Init_Attribute(xml_attr *attr);
void Init_String(t_string *ts);
void Init_Triplet_Struct(triplet *triplet);
void Init_Efrq(t_efrq *f);
void M4_Init_Model(m4 *m4mod, calign *data, t_mod *mod);
void RATES_Init_Rate_Struct(t_rate *rates, t_rate *existing_rates, int n_otu);
void Init_Rmat(t_rmat *rmat);
void Init_MGF_Bl(t_tree *tree);
int Init_Qmat_FLU(phydbl *daa, phydbl *pi);
void Set_Defaults_Ras(t_ras *ras);
void GEO_Init_Coord(t_geo_coord *t, int n_dim);
void PHYREX_Init_Disk_Event(t_dsk *t, int n_dim, t_phyrex_mod *mod);
void PHYREX_Init_Lindisk_Node(t_ldsk *t, t_dsk *devt, int n_dim);
void PHYREX_Init_Migrep_Mod(t_phyrex_mod *t, int n_dim, phydbl max_lat, phydbl max_lon);
void MCMC_Init_MCMC_Struct(char *filename, option *io, t_mcmc *mcmc);

#endif
