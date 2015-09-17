/*

PHYML :  a program that  computes maximum likelihood  phylogenies from
DNA or AA homologous sequences

Copyright (C) Stephane Guindon. Oct 2003 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/

#include <config.h>

#ifndef IO_H
#define IO_H

#include "utilities.h"

t_tree *Read_Tree(char **s_tree);
void R_rtree(char *s_tree_a,char *s_tree_d,t_node *a,t_tree *tree,int *n_int,int *n_ext);
void Read_Branch_Label(char *s_d,char *s_a,t_edge *b);
void Read_Branch_Length(char *s_d,char *s_a,t_tree *tree);
void Read_Node_Name(t_node *d,char *s_tree_d,t_tree *tree);
void Clean_Multifurcation(char **subtrees,int current_deg,int end_deg);
char **Sub_Trees(char *tree,int *degree);
int Next_Par(char *s,int pos);
void Print_Tree(FILE *fp,t_tree *tree);
char *Write_Tree(t_tree *tree,int custom);
void R_wtree(t_node *pere,t_node *fils,int *available,char **s_tree,t_tree *tree);
void R_wtree_Custom(t_node *pere,t_node *fils,int *available,char **s_tree,int *pos,t_tree *tree);
void Detect_Align_File_Format(option *io);
void Detect_Tree_File_Format(option *io);
align **Get_Seq(option *io);
void Get_Nexus_Data(FILE *fp,option *io);
int Get_Token(FILE *fp,char *token);
align **Get_Seq_Phylip(option *io);
void Read_Ntax_Len_Phylip(FILE *fp,int *n_otu,int *n_tax);
align **Read_Seq_Sequential(option *io);
align **Read_Seq_Interleaved(option *io);
int Read_One_Line_Seq(align ***data,int num_otu,FILE *in);
t_tree *Read_Tree_File(option *io);
char *Return_Tree_String_Phylip(FILE *fp_input_tree);
t_tree *Read_Tree_File_Phylip(FILE *fp_input_tree);
void Print_Site_Lk(t_tree *tree,FILE *fp);
void Print_Seq(FILE *fp, align **data, int n_otu);
void Print_CSeq(FILE *fp,int compressed,calign *cdata);
void Print_CSeq_Select(FILE *fp,int compressed,calign *cdata,t_tree *tree);
void Print_Dist(matrix *mat);
void Print_Node(t_node *a,t_node *d,t_tree *tree);
void Print_Model(t_mod *mod);
void Print_Mat(matrix *mat);
FILE *Openfile(char *filename,int mode);
void Print_Fp_Out(FILE *fp_out,time_t t_beg,time_t t_end,t_tree *tree,option *io,int n_data_set,int num_tree, int add_citation);
void Print_Fp_Out_Lines(FILE *fp_out,time_t t_beg,time_t t_end,t_tree *tree,option *io,int n_data_set);
void Print_Freq(t_tree *tree);
void Print_Settings(option *io);
void Print_Banner(FILE *fp);
void Print_Banner_Small(FILE *fp);
void Print_Data_Set_Number(option *io,FILE *fp);
void Print_Lk(t_tree *tree,char *string);
void Print_Pars(t_tree *tree);
void Print_Lk_And_Pars(t_tree *tree);
void Read_Qmat(phydbl *daa,phydbl *pi,FILE *fp);
void Print_Qmat_AA(phydbl *daa,phydbl *pi);
void Print_Square_Matrix_Generic(int n,phydbl *mat);
void Print_Diversity(FILE *fp,t_tree *tree);
void Print_Diversity_Pre(t_node *a,t_node *d,t_edge *b,FILE *fp,t_tree *tree);
t_tree *Read_User_Tree(calign *cdata,t_mod *mod,option *io);
void Print_Time_Info(time_t t_beg,time_t t_end);
void PhyML_Printf(char *format,...);
void PhyML_Fprintf(FILE *fp,char *format,...);
void Read_Clade_Priors(char *file_name,t_tree *tree);
option *Get_Input(int argc,char **argv);
void Print_Data_Structure(int final, FILE *fp, t_tree *root);
int Set_Whichmodel(int select);
void Print_Site(calign *cdata, int num, int n_otu, char *sep, int stepsize, FILE *fp);
option *PhyML_XML(char *xml_filename);
void Check_Taxa_Sets(t_tree *mixt_tree);
void Make_Ratematrice_From_XML_Node(xml_node *instance, option *io, t_mod *mod);
void Make_Efrq_From_XML_Node(xml_node *instance, option *io, t_mod *mod);
void Make_Topology_From_XML_Node(xml_node *instance, option *io, t_mod *mod);
void Make_RAS_From_XML_Node(xml_node *parent, t_mod *mod);
void Post_Process_Data(option *io);
int *Return_Int(int in);
void Print_All_Edge_PMats(t_tree* tree);
void Print_All_Edge_Likelihoods(t_tree* tree);
void Print_Edge_Likelihoods(t_tree* tree, t_edge* b, bool scientific);
void Print_Edge_PMats(t_tree* tree, t_edge* b);
void Print_Tip_Partials(t_tree* tree, t_node* d);
void Dump_Arr_D(phydbl* arr, int num);
void Dump_Arr_S(short int* arr, int num);
void Dump_Arr_I(int* arr, int num);
void Print_Tree_Structure(t_tree* tree);
void Print_Node_Brief(t_node *a, t_node *d, t_tree *tree, FILE *fp);
void Generic_Exit(const char *file, int line, const char *function);
#endif
