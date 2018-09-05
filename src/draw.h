#include <config.h>

#ifndef DRAW_H
#define DRAW_H

#include "utilities.h"

void DR_Dist_To_Root_Pre(t_node *a, t_node *d, t_edge *b, t_tree *tree);
void DR_Dist_To_Root(t_node *n_root, t_tree *tree);
void DR_Get_X_Coord_Pre(t_node *a, t_node *d, t_edge *b, tdraw *w, int fixed_tips, t_tree *tree);
void DR_Get_X_Coord(int fixed_tips, tdraw *w, t_tree *tree);
tdraw *DR_Make_Tdraw_Struct(t_tree *tree);
void DR_Init_Tdraw_Struct(tdraw *d);
void DR_Get_Tree_Box_Width(tdraw *w, t_tree *tree);
void DR_Get_Y_Coord_Post(t_node *a, t_node *d, t_edge *b, int *next_y_slot, int fixed_tips, tdraw *w, t_tree *tree);
void DR_Get_Y_Coord(int fixed_tips, tdraw *w, t_tree *tree);
void DR_Get_Tree_Coord(t_tree *tree);
phydbl DR_Get_Max_Dist_To_Root(t_tree *tree);
void DR_Print_Tree_Postscript(int tree_num, int render_name,FILE *fp, t_tree *tree);
void DR_Print_Tree_Postscript_Pre(t_node *a, t_node *d, t_edge *b, int render_name, FILE *fp, tdraw *w, t_tree *tree);
void DR_Print_Postscript_EOF(FILE *fp);
void DR_Print_Postscript_Header(int n_pages, FILE *fp);
void DR_Get_Tree_Coord_Scaled(tdraw *w, t_tree *tree);
void DR_Get_Cdf_Mat(t_tree *tree);
void DR_Draw_Tree(char *file_name, t_tree *tree);



#endif
