/*

PhyML:  a program that  computes maximum likelihood phyLOGenies from
DNA or AA homoLOGous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/


/* Routines for molecular clock trees and molecular dating */


#include "tiporder.h"

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


/* int TIPO_main(int argc, char **argv) */
/* { */
/*   t_tree **list_tree,*ref_tree,*tree; */
/*   FILE *fp_ref_tree,*fp_list_tree,*fp_coord,*ps_tree; */
/*   int i,j,k; */
/*   int n_trees; */
/*   option *ref_io,*list_io; */
/*   char **name_table; */
/*   int r_seed; */

/*   r_seed = time(NULL); */
/*   srand(r_seed); */


/*   ref_io  = (option *)Make_Input(); */
/*   list_io = (option *)Make_Input(); */

/*   fp_ref_tree  = (FILE *)fopen(argv[1],"r"); */
/*   fp_list_tree = (FILE *)fopen(argv[2],"r"); */
/*   fp_coord     = (FILE *)fopen(argv[3],"r"); */

/*   if(!fp_ref_tree)  */
/*     { */
/*       PhyML_Printf("\n. Can't find %s",argv[1]); */
/*       PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__); */
/*       Warn_And_Exit(""); */
/*     } */

/*   if(!fp_list_tree)  */
/*     { */
/*       PhyML_Printf("\n. Can't find %s",argv[2]); */
/*       PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__); */
/*       Warn_And_Exit(""); */
/*     } */

/*   ref_io->fp_in_tree = fp_ref_tree; */
/*   Read_Tree_File(ref_io); */
/*   fclose(ref_io->fp_in_tree); */
/*   ref_tree = ref_io->treelist->tree[0]; */
/*   ref_tree->io = ref_io; */

  
/*   list_io->fp_in_tree = fp_list_tree; */
/*   Read_Tree_File(list_io); */
/*   fclose(list_io->fp_in_tree); */
/*   list_tree = list_io->treelist->tree; */
/*   n_trees = list_io->treelist->list_size; */
/*   PhyML_Printf("\n. Read %d trees\n",n_trees); */

/*   For(i,n_trees) list_tree[i]->io = list_io; */

/*   name_table = (char **)mCalloc(ref_tree->n_otu,sizeof(char **)); */
/*   For(i,ref_tree->n_otu) name_table[i] = (char *)mCalloc(T_MAX_NAME,sizeof(char)); */


/*   /\* Sort translation table such that tree->a_nodes[i]->name == tree->io->short_tax_name[i] for all i *\/ */
/*   TIPO_Sort_Translation_Table(ref_tree); */

/*   ref_tree->io->z_scores = (phydbl *)mCalloc(ref_tree->n_otu,sizeof(phydbl)); */

/* /\*   TIPO_Read_Taxa_Zscores(fp_coord,ref_tree); *\/ */

/*   For(i,ref_tree->n_otu) ref_tree->io->z_scores[i] = TIPO_Read_One_Taxon_Zscore(fp_coord,ref_tree->a_nodes[i]->name,ref_tree); */

/*   TIPO_Normalize_Zscores(ref_tree); */

/*   /\* Find matching tips *\/ */
/*   For(i,n_trees) */
/*     { */
/*       TIPO_Sort_Translation_Table(list_tree[i]); */

/*       For(j,ref_tree->n_otu)  */
/* 	{ */
/* 	  For(k,ref_tree->n_otu)  */
/* 	    { */
/* 	      if(!strcmp(ref_io->long_tax_names[j],list_io->long_tax_names[k])) */
/* 		{ */
/* 		  list_tree[i]->a_nodes[k]->ext_node = ref_tree->a_nodes[j]; */
/* 		  break; */
/* 		} */
/* 	    } */
/* 	  if(k == ref_tree->n_otu) */
/* 	    { */
/* 	      PhyML_Printf("\n. Could not find matching tips for \"%s\" (tree %d)",ref_tree->a_nodes[j]->name,i); */
/* 	      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__); */
/* 	      Warn_And_Exit(""); */
/* 	    } */
/* 	} */
/*     } */

/*   PhyML_Printf("\n. Getting ancestors"); fflush(NULL); */
/*   Update_Ancestors(ref_tree->n_root,ref_tree->n_root->v[2],ref_tree); */
/*   Update_Ancestors(ref_tree->n_root,ref_tree->n_root->v[1],ref_tree); */

/*   For(i,n_trees)    */
/*     { */
/*       Update_Ancestors(list_tree[i]->n_root,list_tree[i]->n_root->v[2],list_tree[i]); */
/*       Update_Ancestors(list_tree[i]->n_root,list_tree[i]->n_root->v[1],list_tree[i]); */
/*       list_tree[i]->n_root->anc = NULL; */
/*     } */

/*   PhyML_Printf("\n. Getting bipartitions"); fflush(NULL); */

/*   Free_Bip(ref_tree); */
/*   Alloc_Bip(ref_tree); */
/*   Get_Bip(ref_tree->a_nodes[0], */
/* 	  ref_tree->a_nodes[0]->v[0], */
/* 	  ref_tree); */

/*   For(i,n_trees)  */
/*     { */
/*       if(!(i%10)) */
/*       PhyML_Printf("\n. Getting bipartition for tree %d",i); */
/*       Free_Bip(list_tree[i]); */
/*       Alloc_Bip(list_tree[i]); */
/*       Get_Bip(list_tree[i]->a_nodes[0], */
/* 	      list_tree[i]->a_nodes[0]->v[0], */
/* 	      list_tree[i]); */
/*     } */


/*   PhyML_Printf("\n. Getting tip ranks"); fflush(NULL); */
/* /\*   TIPO_Get_Tips_Y_Rank(ref_tree); *\/ */

/*   TIPO_Get_Tips_Y_Rank_From_Zscores(ref_tree); */

/* /\*   PhyML_Printf("\n. Minimizing"); fflush(NULL); *\/ */
/* /\*   TIPO_Minimize_Tip_Order_Score(n_trees,list_tree,ref_tree); *\/ */

/*   PhyML_Printf("\n. N_OTU = %d",ref_tree->n_otu); */
/*   TIPO_Untangle_Tree(ref_tree); */
/*   PhyML_Printf("\n ** ORIGINAL %f",ref_tree->tip_order_score); */

/* /\*   i = 0; *\/ */
/* /\*   do *\/ */
/* /\*     { *\/ */
/* /\* /\\*       TIPO_Get_Tips_Y_Rank_From_Zscores(ref_tree); *\\/ *\/ */
/* /\*       TIPO_Randomize_Tip_Y_Ranks(ref_tree); *\/ */
/* /\*       TIPO_Untangle_Tree(ref_tree); *\/ */
/* /\*       PhyML_Printf("\n ** %f",ref_tree->tip_order_score); *\/ */
/* /\*       i++; *\/ */
/* /\*     }while(i < 5000); *\/ */

/*   Test_Node_Table_Consistency(ref_tree); */

/*   ref_tree->ps_tree = DR_Make_Tdraw_Struct(ref_tree); */
/*   DR_Get_Tree_Coord(ref_tree); */
/*   For(j,ref_tree->n_otu)  */
/*     { */
/*       ref_tree->ps_tree->ycoord[j] =  */
/* 	(ref_tree->a_nodes[j]->y_rank/ref_tree->n_otu)* */
/* 	ref_tree->ps_tree->page_height; */
/*     } */

/*   ps_tree  = (FILE *)fopen("order_tree.ps","w"); */


/*   list_io->z_scores = (phydbl *)mCalloc(ref_tree->n_otu,sizeof(phydbl)); */

/*   DR_Print_Postscript_Header(1,ps_tree); */
/*   For(i,n_trees) */
/*     { */
/*       tree = list_tree[i]; */

/*       Test_Node_Table_Consistency(tree); */
/*       tree->rates = RATES_Make_Rate_Struct(tree->n_otu); */
/*       RATES_Init_Rate_Struct(tree->rates,tree->n_otu); */
/*       TIMES_Least_Square_Node_Times(tree->e_root,tree); */
/*       TIMES_Adjust_Node_Times(tree); */
/*       RATES_Update_Cur_Bl(tree); */

/*       tree->ps_tree = DR_Make_Tdraw_Struct(tree); */
/*       DR_Init_Tdraw_Struct(tree->ps_tree); */
/*       DR_Get_Tree_Box_Width(tree->ps_tree,tree); */
/*       Dist_To_Root(tree->n_root,tree); */
/*       tree->ps_tree->max_dist_to_root = DR_Get_Max_Dist_To_Root(tree); */
 
/*       For(j,ref_tree->n_otu) tree->io->z_scores[j] = ref_tree->io->z_scores[tree->a_nodes[j]->ext_node->num]; */
/*       TIPO_Get_Tips_Y_Rank_From_Zscores(tree); */
/*       TIPO_Untangle_Tree(tree); */
/*       For(j,ref_tree->n_otu) tree->ps_tree->ycoord[j] =  (tree->a_nodes[j]->y_rank/tree->n_otu)*tree->ps_tree->page_height; */

/* /\*       For(j,ref_tree->n_otu) tree->ps_tree->ycoord[j] = ref_tree->ps_tree->ycoord[tree->a_nodes[j]->ext_node->num]; *\/ */
/*       For(j,ref_tree->n_otu) list_io->z_scores[j] = ref_io->z_scores[tree->a_nodes[j]->ext_node->num]; */

/*       DR_Get_Y_Coord(YES,tree->ps_tree,tree); */
/*       DR_Get_X_Coord( NO,tree->ps_tree,tree); */

/*       if(!i) DR_Print_Tree_Postscript(1,YES,ps_tree,tree); */
/*       else   DR_Print_Tree_Postscript(1, NO,ps_tree,tree); */
/*     } */
/*   DR_Print_Postscript_EOF(ps_tree); */
/*   fclose(ps_tree); */


/*   ps_tree  = (FILE *)fopen("ref_tree.ps","w"); */

/*   DR_Print_Postscript_Header(1,ps_tree); */
/*   tree = ref_tree; */
/*   tree->rates = RATES_Make_Rate_Struct(tree->n_otu); */
/*   RATES_Init_Rate_Struct(tree->rates,tree->n_otu); */
/*   TIMES_Least_Square_Node_Times(tree->e_root,tree); */
/*   TIMES_Adjust_Node_Times(tree); */
/*   RATES_Update_Cur_Bl(tree); */
/*   DR_Init_Tdraw_Struct(tree->ps_tree); */
/*   DR_Get_Tree_Box_Width(tree->ps_tree,tree); */
/*   Dist_To_Root(tree->n_root,tree); */
/*   tree->ps_tree->max_dist_to_root = DR_Get_Max_Dist_To_Root(tree); */
  
/*   TIPO_Get_Tips_Y_Rank_From_Zscores(tree); */
/* /\*   TIPO_Untangle_Tree(tree); *\/ */
/*   For(j,tree->n_otu) tree->ps_tree->ycoord[j] =  (tree->a_nodes[j]->y_rank/tree->n_otu)*tree->ps_tree->page_height; */

/*   DR_Get_Y_Coord(YES,tree->ps_tree,tree); */
/*   DR_Get_X_Coord( NO,tree->ps_tree,tree); */
/*   DR_Print_Tree_Postscript(1,YES,ps_tree,tree); */
/*   DR_Print_Postscript_EOF(ps_tree); */
/*   fclose(ps_tree); */


/*   PhyML_Printf("\n"); */

/*   fclose(fp_ref_tree); */
/*   fclose(fp_list_tree); */
/*   fclose(fp_coord); */
/* } */ 

int TIPO_main(int argc, char **argv)
{
  t_tree *tree;
  option *io;
  FILE *fp_tree_file, *fp_coord_file;
  int i;

/*   Rprintf("%s\n",tree_file_name[0]); */
/*   Rprintf("%s\n",coord_file_name[0]); */


  srand(time(NULL)); rand();


  fp_tree_file  = (FILE *)fopen(argv[1],"r");
  fp_coord_file = (FILE *)fopen(argv[2],"r");

  io = (option *)Make_Input();
  io->fp_in_tree = fp_tree_file;
  Read_Tree_File(io);
  tree = io->treelist->tree[0];
  tree->io = io;

  tree->io->z_scores = (phydbl *)mCalloc(tree->n_otu,sizeof(phydbl));

  For(i,tree->n_otu) tree->io->z_scores[i] = TIPO_Read_One_Taxon_Zscore(fp_coord_file,tree->a_nodes[i]->name,1,tree);
  /* TIPO_Normalize_Zscores(tree); */
  Free_Bip(tree);
  Alloc_Bip(tree);
  Get_Bip(tree->a_nodes[0],tree->a_nodes[0]->v[0],tree);
  TIPO_Get_Tips_Y_Rank_From_Zscores(tree);
  /* TIPO_Get_Tips_Y_Rank(tree); */
  

  tree->geo_mig_sd = 0.1;
  Generic_Brent_Lk(&(tree->geo_mig_sd),
  		   1.E-3,1.E+2,1.E-6,
  		   100,NO,
  		   &Wrap_Geo_Lk,
  		   NULL,tree,NULL,NO);

  PhyML_Printf("\n. sd=%f",tree->geo_mig_sd);


  fclose(fp_tree_file);
  fclose(fp_coord_file);

  return 0;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* Z_scores have already been recorder here */
void TIPO_Get_Min_Number_Of_Tip_Permut(t_tree *tree)
{
  Update_Ancestors(tree->n_root,tree->n_root->v[2],tree);
  Update_Ancestors(tree->n_root,tree->n_root->v[1],tree);

  Free_Bip(tree);
  Alloc_Bip(tree);
  Get_Bip(tree->a_nodes[0],tree->a_nodes[0]->v[0],tree);

  TIPO_Get_Tips_Y_Rank_From_Zscores(tree);

  TIPO_Untangle_Tree(tree);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void TIPO_Get_Tips_Y_Rank(t_tree *tree)
{
  phydbl curr_rank;

  curr_rank = .0;
  TIPO_Get_Tips_Y_Rank_Pre(tree->n_root,tree->n_root->v[2],&curr_rank,tree);
  TIPO_Get_Tips_Y_Rank_Pre(tree->n_root,tree->n_root->v[1],&curr_rank,tree);
  
  if(curr_rank != tree->n_otu)
    {
      PhyML_Printf("\n. tree->n_otu = %d curr_rank = %d",tree->n_otu,curr_rank);
      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void TIPO_Get_Tips_Y_Rank_Pre(t_node *a, t_node *d, phydbl *curr_rank, t_tree *tree)
{
  if(d->tax) 
    {
      d->y_rank = *curr_rank;
      *curr_rank += 1.;
      return;
    }
  else
    {
      int i;
      For(i,3)
	{
	  if(d->v[i] != a && d->b[i] != tree->e_root)
	    {
	      TIPO_Get_Tips_Y_Rank_Pre(d,d->v[i],curr_rank,tree);
	    }
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void TIPO_Get_All_Y_Rank(t_tree *tree)
{
  tree->sum_y_dist_sq = .0;
  tree->sum_y_dist    = .0;
  TIPO_Get_All_Y_Rank_Pre(tree->n_root,tree->n_root->v[2],tree);
  TIPO_Get_All_Y_Rank_Pre(tree->n_root,tree->n_root->v[1],tree);
  tree->n_root->y_rank = (tree->n_root->v[2]->y_rank+tree->n_root->v[1]->y_rank)/2.;
  tree->n_root->y_rank_min = MIN(tree->n_root->v[2]->y_rank_min,tree->n_root->v[1]->y_rank_min);
  tree->n_root->y_rank_max = MAX(tree->n_root->v[2]->y_rank_max,tree->n_root->v[1]->y_rank_max);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void TIPO_Get_All_Y_Rank_Pre(t_node *a, t_node *d, t_tree *tree)
{
  if(d->tax) 
    {
      d->y_rank_min = d->y_rank;
      d->y_rank_max = d->y_rank;
      return;
    }
  else
    {
      int i;
      int dir1,dir2;
      phydbl v1,v2;

      For(i,3)
	{
	  if(d->v[i] != a && d->b[i] != tree->e_root)
	    {
	      TIPO_Get_All_Y_Rank_Pre(d,d->v[i],tree);
	    }
	}

      dir1 = dir2 = -1;
      For(i,3)
	{
	  if(d->v[i] != a && d->b[i] != tree->e_root)
	    {
	      if(dir1 < 0) dir1 = i;
	      else         dir2 = i;
	    }
	}

      v1 = d->v[dir1]->y_rank;
      v2 = d->v[dir2]->y_rank;

      d->y_rank            = (v1+v2)/2.;
      tree->sum_y_dist_sq += POW(v1-v2,2);
      tree->sum_y_dist    += FABS(v1-v2);
      d->y_rank_min        = MIN(d->v[dir1]->y_rank_min,d->v[dir2]->y_rank_min);
      d->y_rank_max        = MAX(d->v[dir1]->y_rank_max,d->v[dir2]->y_rank_max);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void TIPO_Swap_One_Node(t_node *d, t_tree *tree)
{
  if(d->tax) return;
  else
    {
      int i;
      int dir1, dir2;
      t_node *tmp_n;
      t_edge *tmp_e;

      if(d != tree->n_root)
	{
	  dir1 = dir2 = -1;
	  For(i,3)
	    {
	      if((d->v[i] != d->anc) && (d->b[i] != tree->e_root))
		{
		  if(dir1 < 0) dir1 = i;
		  else         dir2 = i;
		}
	    }
	}
      else
	{
	  dir1 = 0;
	  dir2 = 1;
	}

      tmp_n      = d->v[dir2];
      d->v[dir2] = d->v[dir1];
      d->v[dir1] = tmp_n;

      tmp_e      = d->b[dir2];
      d->b[dir2] = d->b[dir1];
      d->b[dir1] = tmp_e;
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void TIPO_Minimize_Tip_Order_Score(int n_trees, t_tree **list_tree, t_tree *ref_tree)
{
  int i,j;
  phydbl score,min_score,old_min_score;
  phydbl diff,eps;
  t_node **node_table;
  int swapped;
  t_node *tmp;

  eps = 1.E-3;
  old_min_score = min_score = score = INT_MAX;
  /*   TIPO_Print_Tip_Ordered(ref_tree); */

  do
    {
      for(i=ref_tree->n_otu;i<2*ref_tree->n_otu-1;i++)
	{	  
	  TIPO_Swap_One_Node(ref_tree->a_nodes[i],ref_tree);
	  TIPO_Get_Tips_Y_Rank(ref_tree);
	  score = (phydbl)TIPO_Untangle_Tree_List(n_trees,list_tree,ref_tree);
	  if(score == -1) 
	    {
	      return;
	    }

	  if(score < min_score) 
	    {
	      min_score = score;
	      PhyML_Printf("\n- Score = %f",score);
	    }
	  else 
	    {
	      TIPO_Swap_One_Node(ref_tree->a_nodes[i],ref_tree);    
	      TIPO_Get_Tips_Y_Rank(ref_tree);
/* 	      PhyML_Printf("\n+ Score = %f",score); */
	    }
	}
      diff = fabs(old_min_score - min_score);
      old_min_score = min_score;
    }while(diff > eps);

  PhyML_Printf("\n");

  node_table = (t_node **)mCalloc(ref_tree->n_otu,sizeof(t_node *));


  For(i,ref_tree->n_otu)
    {
      For(j,ref_tree->n_otu)
	{
	  if(!strcmp(ref_tree->io->short_tax_names[i],ref_tree->a_nodes[j]->name))
	    {
	      Free(ref_tree->a_nodes[j]->name);
	      ref_tree->a_nodes[j]->name = (char *)mCalloc((int)strlen(ref_tree->io->long_tax_names[i])+1,sizeof(char));
	      strcpy(ref_tree->a_nodes[j]->name,ref_tree->io->long_tax_names[i]);
	      break;
	    }
	}
    }

  For(i,ref_tree->n_otu) node_table[i] = ref_tree->a_nodes[i];

/*       bubble sort of conflict nodes according to their y_rank */
  do
    {
      swapped = NO;
      For(i,ref_tree->n_otu-1)
	{
	  if(node_table[i]->y_rank > node_table[i+1]->y_rank)
	    {
	      swapped = YES;
	      tmp             = node_table[i];
	      node_table[i]   = node_table[i+1];
	      node_table[i+1] = tmp;
	    }
	}
    }while(swapped == YES);

  For(i,ref_tree->n_otu)
    {
      PhyML_Printf("\n%s",node_table[i]->name,node_table[i]->y_rank);
    }



  Free(node_table);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void TIPO_Print_Tip_Ordered(t_tree *tree)
{
  TIPO_Print_Tip_Ordered_Pre(tree->n_root,tree->n_root->v[2],tree);
  TIPO_Print_Tip_Ordered_Pre(tree->n_root,tree->n_root->v[1],tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void TIPO_Print_Tip_Ordered_Pre(t_node *a, t_node *d, t_tree *tree)
{
  
  if(d->tax)
    {
      PhyML_Printf("\n. %f \"%s\"",d->y_rank,d->name);
    }
  else
    {
      int i,dir1,dir2;

      dir1 = dir2 = -1;
      For(i,3)
	{
	  if((d->v[i] != a) && (d->b[i] != tree->e_root))
	    {
	      if(dir1 < 0) dir1 = i;
	      else         dir2 = i;
	    }
	}
      if(d->v[dir1]->y_rank < d->v[dir2]->y_rank)
	{
	  TIPO_Print_Tip_Ordered_Pre(d,d->v[dir1],tree);
	  TIPO_Print_Tip_Ordered_Pre(d,d->v[dir2],tree);
	}
      else
	{
	  TIPO_Print_Tip_Ordered_Pre(d,d->v[dir2],tree);
	  TIPO_Print_Tip_Ordered_Pre(d,d->v[dir1],tree);
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


int TIPO_Untangle_Tree_List(int n_trees, t_tree **list_tree, t_tree *ref_tree)
{
  int i,j;
  int tree_score,score;

  score = 0;
  For(i,n_trees) 
    {
/*       PhyML_Printf("\n. Untangling tree %3d",i); */
      For(j,ref_tree->n_otu) list_tree[i]->a_nodes[j]->y_rank = list_tree[i]->a_nodes[j]->ext_node->y_rank;
      tree_score = TIPO_Untangle_Tree(list_tree[i]);
/*       PhyML_Printf(" score = %3d",tree_score); */
      score += tree_score;
      if(tree_score < 0) 
	{
	  return -1;
	}
    }

  return score;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl TIPO_Untangle_Tree(t_tree *tree)
{
  int conflict;
  int n_trials;
  t_node **node_table;
  int i,swapped;
  t_node *tmp;

  node_table = (t_node **)mCalloc(tree->n_otu,sizeof(t_node *));

  For(i,tree->n_otu) node_table[i] = tree->a_nodes[i];
  For(i,tree->n_otu) tree->a_nodes[i]->y_rank_ori = tree->a_nodes[i]->y_rank;


/* bubble sort of nodes according to their y_rank */
  do
    {
      swapped = NO;
      For(i,tree->n_otu-1)
	{
	  if(node_table[i]->y_rank > node_table[i+1]->y_rank)
	    {
	      swapped = YES;
	      tmp             = node_table[i];
	      node_table[i]   = node_table[i+1];
	      node_table[i+1] = tmp;
	    }
	}
    }
  while(swapped == YES);
  

  /* Work out the y_rank values for every internal node given the external node ranks */
  TIPO_Get_All_Y_Rank(tree);

  tree->tip_order_score = .0;
      
  n_trials = 0;
  do
    {
      conflict= NO;
      /* Recusrssive untangling of the tree */ 
      TIPO_Untangle_Node(tree->n_root,tree->n_root->v[2],node_table,&conflict,tree);
      TIPO_Untangle_Node(tree->n_root,tree->n_root->v[1],node_table,&conflict,tree);
      n_trials++;


      if(n_trials > 2) /* We should have been able to untangle the tree after just one tree traversal */
	{
	  int i;
	  FILE *ps_tree;

	  PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
	  ps_tree = (FILE *)fopen("failed_tree.ps","w");

	  Test_Node_Table_Consistency(tree);
	  tree->rates = RATES_Make_Rate_Struct(tree->n_otu);
	  RATES_Init_Rate_Struct(tree->rates,tree->io->rates,tree->n_otu);
	  TIMES_Least_Square_Node_Times(tree->e_root,tree);
	  TIMES_Adjust_Node_Times(tree);
	  RATES_Update_Cur_Bl(tree);

	  DR_Print_Postscript_Header(1,ps_tree);
	  tree->ps_tree = DR_Make_Tdraw_Struct(tree);
	  DR_Init_Tdraw_Struct(tree->ps_tree);
	  DR_Get_Tree_Box_Width(tree->ps_tree,tree);
	  Dist_To_Root(tree->n_root,tree);
	  tree->ps_tree->max_dist_to_root = DR_Get_Max_Dist_To_Root(tree);
	  For(i,tree->n_otu) tree->ps_tree->ycoord[i] = tree->a_nodes[i]->y_rank * (int)(tree->ps_tree->page_height / (tree->n_otu));
	  DR_Get_X_Coord(NO,tree->ps_tree,tree);
	  DR_Get_Y_Coord(YES,tree->ps_tree,tree);
	  DR_Print_Tree_Postscript(1,NO,ps_tree,tree);
	  DR_Print_Postscript_EOF(ps_tree);
	  fclose(ps_tree);
	  Warn_And_Exit("");	  
	}
    }while(conflict == YES);  

  Free(node_table);
  return tree->tip_order_score;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void TIPO_Untangle_Node(t_node *a, t_node *d, t_node **node_table, int *conflict, t_tree *tree)
{

  if(d->tax) return;
  else
    {
      int    i,j;
      int    d_a;
      int    beg,end;
      phydbl min,max;
      t_node *lca;
      t_node **conflict_tips, **anc_conflict;
      int    n_conflicts;
      phydbl eps,tmp_rank;
      t_node *tmp_node;
      int    n_moved;
      int    n_anc_conflicts;

      anc_conflict = NULL;
      d_a = -1;

      /* It is a post order traversal */
      For(i,3)
	{
	  if((d->v[i] != d->anc) && (d->b[i] != tree->e_root))
	    {
	      TIPO_Untangle_Node(d,d->v[i],node_table,conflict,tree);
	    }
	}
      
      
      lca = NULL;
      eps = 1.E-10;

      /* Find direction fron node d ((d)escendant) to a ((a)ncestor) */
      For(i,3)
	{
	  if((d->v[i] == d->anc) || (d->b[i] == tree->e_root))
	    {
	      d_a = i;
	      break;
	    }
	}

      
      /* y_rank_min is the minimum rank among all the ranks of the tips that 
	 can be reached when going from a to d */
      min = d->y_rank_min;
      max = d->y_rank_max;

      /* Get the list of tip nodes which ranks are between d->y_rank_min and d->y_rank_max */
      n_conflicts = 0;
      For(i,tree->n_otu)
	{
	  if((node_table[i]->y_rank > min - eps) && (node_table[i]->y_rank < max + eps))
	    {
	      n_conflicts++;	      
	    }
	}
      
      conflict_tips = NULL;
      For(i,tree->n_otu)
	{
	  if(node_table[i]->y_rank > min - eps)
	    {
	      conflict_tips = node_table+i;
	      break;
	    }
	}     


      beg = 0;
      end = n_conflicts;
      n_moved = 0;
      n_anc_conflicts = 0;
      do
	{
	  for(i=beg;i<end;i++)
	    {
	      For(j,d->bip_size[d_a]) if(conflict_tips[i] == d->bip_node[d_a][j]) break;
	      
	      if(j == d->bip_size[d_a]) 	      
		/* conflict_tips[i] does not belong to the list of descendant of node d. It is 
		   therefore responsible for a conflict */
		{
		  *conflict = YES;

/* 		  printf("\n. Moving %d with rank %f",conflict_tips[i]->num,conflict_tips[i]->y_rank); */
		  
		  n_moved++;
		 
		  /* Move from conflict_tips[i] towards the root as long as the rank of the node lca is between min and max */
		  lca = conflict_tips[i];
		  while(lca->y_rank_min > min && lca->y_rank_max < max) lca = lca->anc;


		  if(lca->y_rank_min > max && lca->y_rank_max > max)
		    {
		      PhyML_Printf("\n. lca = %d (%d)",lca->num,lca==tree->n_root);
		      PhyML_Printf("\n. Err in file %s at line %d",__FILE__,__LINE__);
		      Warn_And_Exit("");
		    }

		  if(lca->y_rank_min < min && lca->y_rank_max < min)
		    {
		      PhyML_Printf("\n. lca = %d (root ? %d) (tip ? %d)",lca->num,lca==tree->n_root,lca->tax);
		      PhyML_Printf("\n. lca->y_rank_min = %f lca->y_rank_max = %f",lca->y_rank_min,lca->y_rank_max);
		      PhyML_Printf("\n. min=%f max=%f",min,max);
		      PhyML_Printf("\n. Err in file %s at line %d",__FILE__,__LINE__);
		      Warn_And_Exit("");
		    }
		  if(lca->tax)
		    {
		      PhyML_Printf("\n. lca  (%d) cannot be a tip.",lca->num);
		      PhyML_Printf("\n. lca->anc->y_rank=%f",lca->anc->y_rank);
		      PhyML_Printf("\n. lca->y_rank=%f",lca->y_rank);
		      PhyML_Printf("\n. lca->y_rank_min=%f lca=>y_rank_max=%f",lca->y_rank_min,lca->y_rank_max);
		      PhyML_Printf("\n. min=%f max=%f",min,max);
		      PhyML_Printf("\n. %p %p",lca,conflict_tips[i]);
		      PhyML_Printf("\n. Err in file %s at line %d",__FILE__,__LINE__);
		      Warn_And_Exit("");
		    }
		  
		  /* Have you found lca previously ? */
		  For(j,n_anc_conflicts) if(anc_conflict[j] == lca) break;
		  if(j == n_anc_conflicts) /* if no, then update the tree score and the list of ancestral nodes at the origin of conflicts */
		    {
		      /* tree->tip_order_score+=1.; */
		      /* tree->tip_order_score+=(lca->y_rank_max-lca->y_rank_min); */
		      n_anc_conflicts++;
		      anc_conflict = (t_node **)realloc(anc_conflict,n_anc_conflicts*sizeof(t_node *));
		      anc_conflict[n_anc_conflicts-1] = lca;
		    }
		      
		  /* 		  PhyML_Printf("\n. Detected conflict for ``%s'' (rank:%f min=%f max=%f lca=%f)", */
		  /* 			       conflict_tips[i]->name, */
		  /* 			       conflict_tips[i]->y_rank, */
		  /* 			       min,max,lca->y_rank); */
		  
		  /* Solve the conflict by shifting tip nodes to the left or to the right */
		  if(lca->y_rank > d->y_rank)
		    {
		      end--;
/* 		      max-=1.; */
		      max = conflict_tips[end-1]->y_rank;
/* 		      printf("\n. max=%f %f",max,conflict_tips[end-1]->y_rank); */

		      for(j=i;j<n_conflicts-1;j++)			
			{
/* 			  PhyML_Printf("\n+ Moved (%d,%d) from (%f,%f)", */
/* 				       conflict_tips[j]->num,conflict_tips[j+1]->num, */
/* 				       conflict_tips[j]->y_rank,conflict_tips[j+1]->y_rank); */

			  tmp_rank                   = conflict_tips[j]->y_rank;
			  conflict_tips[j]->y_rank   = conflict_tips[j+1]->y_rank;
			  conflict_tips[j+1]->y_rank = tmp_rank;

			  tmp_node           = conflict_tips[j];
			  conflict_tips[j]   = conflict_tips[j+1];
			  conflict_tips[j+1] = tmp_node;

			  tree->tip_order_score += fabs(conflict_tips[j]->y_rank - conflict_tips[j+1]->y_rank)/n_conflicts;

/* 			  PhyML_Printf(" to (%d,%d) (%f,%f)", */
/* 				       conflict_tips[j]->num,conflict_tips[j+1]->num, */
/* 				       conflict_tips[j]->y_rank,conflict_tips[j+1]->y_rank); */
			}
		    }
		  else
		    {
		      beg++;
/* 		      min+=1.; */
		      min = conflict_tips[beg]->y_rank;
/* 		      printf("\n. min=%f %f",min,conflict_tips[beg]->y_rank); */

		      for(j=i;j>0;j--)
			{
/* 			  PhyML_Printf("\n- Moved (%d,%d) from (%f,%f)", */
/* 				       conflict_tips[j]->num,conflict_tips[j-1]->num, */
/* 				       conflict_tips[j]->y_rank,conflict_tips[j-1]->y_rank); */

		
			  tmp_rank                   = conflict_tips[j]->y_rank;
			  conflict_tips[j]->y_rank   = conflict_tips[j-1]->y_rank;
			  conflict_tips[j-1]->y_rank = tmp_rank;
			  
			  tmp_node           = conflict_tips[j];
			  conflict_tips[j]   = conflict_tips[j-1];
			  conflict_tips[j-1] = tmp_node;		     

			  tree->tip_order_score += fabs(conflict_tips[j]->y_rank - conflict_tips[j+1]->y_rank)/n_conflicts;

/* 			  PhyML_Printf(" to (%d,%d) (%f,%f)", */
/* 				       conflict_tips[j]->num,conflict_tips[j-1]->num, */
/* 				       conflict_tips[j]->y_rank,conflict_tips[j-1]->y_rank); */
			}
		    }

/* 		  printf("\n"); */
/* 		  For(j,n_conflicts) printf("%.0f ",conflict_tips[j]->y_rank); */
/* 		  printf("\n. min=%f max=%f",min,max); */
/* 		  printf("\n. Node %d has now rank %f",c_node->num,c_node->y_rank); */


		  /* Update internal nodes y ranks */
		  TIPO_Get_All_Y_Rank(tree);

		  break;
		}
	    }
	}while(n_moved + d->bip_size[d_a] != n_conflicts);


      For(i,tree->n_otu)
	{
	  if((tree->a_nodes[i]->y_rank > min - eps) && (tree->a_nodes[i]->y_rank < max + eps))
	    {
	      For(j,d->bip_size[d_a])
		{
		  if(tree->a_nodes[i] == d->bip_node[d_a][j]) break;
		}
	      if(j == d->bip_size[d_a])
		{
		  printf("\n. Conflict remaining for node %d (%d)",d->num,a->num);
		  PhyML_Printf("\n. Err in file %s at line %d",__FILE__,__LINE__);
		  Warn_And_Exit("");
		}
	    }
	}

      if(anc_conflict) Free(anc_conflict);

      return;
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


int TIPO_Check_Tip_Ranks(t_tree *tree)
{
  int i,j;
  phydbl eps;

  eps = 1.E-6;

  For(i,tree->n_otu-1)
    {
      for(j=i+1;j<tree->n_otu;j++)
	{
	  if(fabs(tree->a_nodes[i]->y_rank - tree->a_nodes[j]->y_rank) < eps)
	    {
	      return 0;
	    }
	}
    }
  return 1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void TIPO_Read_Taxa_Zscores(FILE *fp_coord, t_tree *tree)
{
  char *name,*line;
  phydbl z;
  int i;

  name = (char *)mCalloc(T_MAX_NAME,sizeof(char));
  line = (char *)mCalloc(T_MAX_LINE,sizeof(char));

  if(!fgets(line,T_MAX_LINE,fp_coord))
    {
      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
  Free(line);


  do
    {
      if(fscanf(fp_coord,"%s\t%lf\n",name,&z) == EOF) break;
      PhyML_Printf("\n. Read %s. Z-score: %f",name,z);

      For(i,tree->n_otu) if(!strcmp(tree->io->long_tax_names[i],name)) break;
      
      if(i == tree->n_otu)
	{
	  PhyML_Printf("\n. Could not find taxon '%s' in coordinate file.",name);
	  PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}

      tree->io->z_scores[i] = z;
      
    }while(1);

  Free(name);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void TIPO_Read_Taxa_Coordinates(FILE *fp_coord, t_tree *tree)
{
  char *name,*line;
  phydbl lon, lat;
  int i;

  name = (char *)mCalloc(T_MAX_NAME,sizeof(char));
  line = (char *)mCalloc(T_MAX_LINE,sizeof(char));

  if(!fgets(line,T_MAX_LINE,fp_coord))
    {
      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
  Free(line);

  tree->io->lat = (phydbl *)mCalloc(tree->n_otu,sizeof(phydbl));
  tree->io->lon = (phydbl *)mCalloc(tree->n_otu,sizeof(phydbl));

  do
    {
      if(fscanf(fp_coord,"%s\t%lf\t%lf\n",name,&lat,&lon) == EOF) break;
      PhyML_Printf("\n. Read %s %f %f",name,lat,lon);

      For(i,tree->n_otu) if(!strcmp(tree->io->long_tax_names[i],name)) break;
      
      if(i == tree->n_otu)
	{
	  PhyML_Printf("\n. Could not find taxon '%s' in coordinate file.",name);
	  PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}

      tree->io->lat[i] = lat;
      tree->io->lon[i] = lon;
      
    }while(1);
	
  Free(name);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void TIPO_Get_Tips_Y_Rank_From_Zscores(t_tree *tree)
{
  int i;

  For(i,tree->n_otu) tree->a_nodes[i]->y_rank = .0;

  /* Randomization in order to avoid ties */
  For(i,tree->n_otu) tree->io->z_scores[i] += Rnorm(0.0,0.001);

  For(i,tree->n_otu) tree->a_nodes[i]->y_rank = tree->io->z_scores[i];

/*   For(i,tree->n_otu) tree->a_nodes[i]->y_rank = .0; */

/*   For(i,tree->n_otu-1) */
/*     { */
/*       for(j=i+1;j<tree->n_otu;j++) */
/* 	{ */
/* 	  if(tree->io->z_scores[i] > tree->io->z_scores[j]) */
/* 	    { */
/* 	      tree->a_nodes[i]->y_rank += 1.0; */
/* 	    } */
/* 	  else */
/* 	  if(tree->io->z_scores[i] < tree->io->z_scores[j]) */
/* 	    { */
/* 	      tree->a_nodes[j]->y_rank += 1.0; */
/* 	    } */
/* 	  else */
/* 	    { */
/* 	      PhyML_Printf("\n. Ties not allowed.\n"); */
/* 	      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__); */
/* 	      Warn_And_Exit(""); */
/* 	    } */
/* 	} */
/*     } */

/*   For(i,tree->n_otu) printf("- %f\n",tree->a_nodes[i]->y_rank); */

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* Sort translation table such that tree->a_nodes[i]->name == tree->io->short_tax_name[i] for all i */
void  TIPO_Sort_Translation_Table(t_tree *tree)
{
  int i,j;
  char *s;


  Test_Node_Table_Consistency(tree);

  For(i,tree->n_otu-1)
    {
      for(j=i+1;j<tree->n_otu;j++)
	{
	  if(!strcmp(tree->a_nodes[i]->name,tree->io->short_tax_names[j]))
	    {
	      s = tree->io->short_tax_names[i];
	      tree->io->short_tax_names[i] = tree->io->short_tax_names[j];
	      tree->io->short_tax_names[j] = s;

	      s = tree->io->long_tax_names[i];
	      tree->io->long_tax_names[i] = tree->io->long_tax_names[j];
	      tree->io->long_tax_names[j] = s;
	      
	      break;
	    }
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void TIPO_Randomize_Tip_Y_Ranks(t_tree *tree)
{

  int i;
  phydbl rk_tmp;
  int rnd_node_num;
  
  For(i,tree->n_otu) tree->a_nodes[i]->y_rank_ori = tree->a_nodes[i]->y_rank;

  For(i,tree->n_otu)
    {
      rnd_node_num = Rand_Int(0,tree->n_otu-1);

      rk_tmp                            = tree->a_nodes[rnd_node_num]->y_rank;
      tree->a_nodes[rnd_node_num]->y_rank = tree->a_nodes[i]->y_rank;
      tree->a_nodes[i]->y_rank            = rk_tmp;
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl TIPO_Read_One_Taxon_Zscore(FILE *fp_coord, char *seqname_qry, int col, t_tree *tree)
{
  char *seqname, *place;
  phydbl lat;

  seqname = (char *)mCalloc(T_MAX_NAME,sizeof(char));
  place   = (char *)mCalloc(T_MAX_NAME,sizeof(char));

  rewind(fp_coord);

  /* skip first line */
/*   if(!fgets(line,T_MAX_LINE,fp_coord)) */
/*     { */
/*       PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__); */
/*       Warn_And_Exit(""); */
/*     } */
/*   Free(line); */

  do
    {
/*       if(fscanf(fp_coord,"%s\t%s\t%lf\t%lf\t%d\n",seqname,place,&lat,&lon,&year) == EOF) */
/*       if(fscanf(fp_coord,"%s\t%s\t%lf\t%lf\n",seqname,place,&lat,&lon) == EOF) */
/*       if(fscanf(fp_coord,"%s\t%lf\t%lf\n",seqname,&lat,&lon) == EOF) */
      if(fscanf(fp_coord,"%s %lf\n",seqname,&lat) == EOF)
	{
	  PhyML_Printf("\n. Could not find sequence '%s' in coordinate file",seqname_qry);
	  PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}


      if(!strcmp(seqname,seqname_qry)) break;
      
    }while(1);
	
/*   PhyML_Printf("\n. Found %20s %s @ %10.2f %10.2f. Recording %10.2f",seqname,place,lat,lon,(col==1)?lat:lon); */

  Free(seqname);
  Free(place);
  
/*   if(col == 1)      return lat; */
/*   else if(col == 2) return lon; */
/*   else              return -1.; */

  return lat;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void TIPO_Normalize_Zscores(t_tree *tree)
{
  int i;
  phydbl min_z,max_z;
  phydbl eps;

  eps = 1.E-10;

  min_z = FLT_MAX;
  For(i,tree->n_otu)
    {
      if(tree->io->z_scores[i] < min_z)
	{
	  min_z = tree->io->z_scores[i];
	}
    }
  
  max_z = -FLT_MAX;
  For(i,tree->n_otu)
    {
      if(tree->io->z_scores[i] > max_z)
	{
	  max_z = tree->io->z_scores[i];
	}
    }

  For(i,tree->n_otu) tree->io->z_scores[i] = (tree->io->z_scores[i] - min_z)/(max_z-min_z+eps); 

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl TIPO_Lk(t_tree *tree)
{
  tree->geo_lnL = 0.0;
  TIPO_Lk_Post(tree->n_root,tree->n_root->v[2],tree);
  TIPO_Lk_Post(tree->n_root,tree->n_root->v[1],tree);
  TIPO_Lk_Core(NULL,tree->n_root,tree);
  return(tree->geo_lnL);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl TIPO_Lk_Post(t_node *a, t_node *d, t_tree *tree)
{
  if(!d->tax)
    {
      int i;

      For(i,3)
	{
	  if(d->v[i] != a && d->b[i] != tree->e_root)
	    {
	      TIPO_Lk_Post(d,d->v[i],tree);
	    }
	}
      TIPO_Lk_Core(a,d,tree);
    }

  return .0;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl TIPO_Lk_Core(t_node *a, t_node *d, t_tree *tree)
{

  int i,j;
  int d_v1,d_v2,v1_d,v2_d;
  t_node *v1, *v2;
  phydbl dist,dens,min_dist;
    

  if(d->tax)
    {
      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  if(d == tree->n_root)
    {
      d_v1 = 0;
      d_v2 = 1;
    }
  else
    {
      d_v1 = d_v2 = -1;
      For(i,3)
	{
	  if(d->v[i] != a && d->b[i] != tree->e_root)
	    {
	      if(d_v1 < 0) d_v1 = i;
	      else         d_v2 = i;
	    }
	}
    }

  v1 = d->v[d_v1];
  v2 = d->v[d_v2];

  v1_d = v2_d = -1;
  if(d == tree->n_root)
    {
      For(i,3)
	{
	  if(v1->b[i] == tree->e_root) v1_d = i;
	  if(v2->b[i] == tree->e_root) v2_d = i;
	}
    }
  else
    {
      For(i,3)
	{
	  if(v1->v[i] == d) v1_d = i;
	  if(v2->v[i] == d) v2_d = i;
	}
    }
  

  dens = 0.0;
  min_dist = FLT_MAX;
  For(i,v1->bip_size[v1_d])
    {
      For(j,v2->bip_size[v2_d])
	{
	  dist = fabs(v1->bip_node[v1_d][i]->y_rank - 
		      v2->bip_node[v2_d][j]->y_rank);

	  if(dist < min_dist) min_dist = dist;

	  dens += Dnorm(dist,0.0,tree->geo_mig_sd);
	  /* printf("\n. dist=%f dens=%f %f %f", */
	  /* 	 dist,Dnorm(dist,0.0,tree->geo_mig_sd), */
	  /* 	 v1->bip_node[v1_d][i]->y_rank, */
	  /* 	 v2->bip_node[v2_d][j]->y_rank); */
	}
    }

  /* printf("\n. min_dist=%f dens=%f", */
  /* 	 min_dist,Dnorm(dist,0.0,tree->geo_mig_sd)); */

  /* dens = Dnorm(min_dist,0.0,tree->geo_mig_sd); */
  tree->geo_lnL += LOG(dens);

  return tree->geo_lnL;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

