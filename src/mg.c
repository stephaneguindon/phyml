/*

PHYML :  a program that  computes maximum likelihood  phylogenies from
DNA or AA homologous sequences 

Copyright (C) Stephane Guindon. Oct 2003 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/

#include <config.h>

#include "mg.h"
#include "free.h"
#include "help.h"
#include "utilities.h"
#include "optimiz.h"
#include "models.h"
#include "simu.h"
#include "lk.h"
#include "pars.h"
#include "interface.h"

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


int PART_main(int argc, char **argv)
{
  option *io;
  char *s_tree;
  FILE *fp_phyml_tree,*fp_phyml_stats,*fp_phyml_lk;
  int part;
  time_t t_beg,t_end;
  div_t hour,min;
  int r_seed;
  int i;

  fflush(NULL);
  io = (option *)Get_Input(argc,argv);
  r_seed = (io->r_seed < 0)?(time(NULL)):(io->r_seed);
  srand(r_seed);
  fp_phyml_stats = Openfile(io->out_stats_file,io->out_stats_file_open_mode);
  PhyML_Fprintf(fp_phyml_stats,"\n- PHYML %s -\n\n", VERSION);
  fp_phyml_tree = Openfile(io->out_tree_file,io->out_tree_file_open_mode);
  fp_phyml_lk = fopen(io->out_lk_file,"w");

  time(&t_beg);

  if(io->multigene)
    {
      align ***data;
      calign **cdata;
      t_mod **mod;
      matrix **mat;
      t_treelist *treelist;
      supert_tree *st;

      data  = (align ***)  mCalloc(io->n_part,sizeof(align **));
      cdata = (calign **)mCalloc(io->n_part,sizeof(calign *));
      mod   = (t_mod **) mCalloc(io->n_part,sizeof(t_mod *));
      mat   = (matrix **)mCalloc(io->n_part,sizeof(matrix *));

      /* Read the sequences (for each partition) */
      For(part,io->n_part)
	{
	  Make_Model_Complete(io->st->optionlist[part]->mod); /* Complete model for each data part */
	  data[part] = Get_Seq(io->st->optionlist[part]);
	  Make_Model_Complete(io->st->optionlist[part]->mod);
	  Set_Model_Name(io->st->optionlist[part]->mod);
	  mod[part]  = io->st->optionlist[part]->mod;

	  PhyML_Printf("\n. Data part [#%d]\n",part+1);
	  PhyML_Printf("\n. Compressing sequences...\n");
	  cdata[part] = Compact_Data(data[part],io->st->optionlist[part]);
	  fclose(io->st->optionlist[part]->fp_in_align);
	  Free_Seq(data[part],cdata[part]->n_otu);
	  Init_Model(cdata[part],mod[part],io->st->optionlist[part]);
	  Check_Ambiguities(cdata[part],
			    io->st->optionlist[part]->mod->io->datatype,
			    io->st->optionlist[part]->mod->io->state_len);
	}

      PART_Make_Supert_tree_Full(io->st,io,cdata);
      st = io->st;
      treelist = st->treelist;
      Fill_Dir_Table(st->tree);
      Update_Dirs(st->tree);

      For(part,io->n_part)
	{
	  st->curr_cdata = cdata[part];
	  if(!PART_Get_Species_Found_In_St(st,cdata[part])) break;
	  treelist->tree[part] = Make_Tree_From_Scratch(st->tree->n_otu,NULL);
	  Copy_Tree_Topology_With_Labels(st->tree,treelist->tree[part]);
 	  treelist->tree[part]->num_curr_branch_available = 0;
	  Connect_Edges_To_Nodes_Recur(treelist->tree[part]->a_nodes[0],
				       treelist->tree[part]->a_nodes[0]->v[0],
				       treelist->tree[part]);
	  PART_Prune_St_Topo(treelist->tree[part],cdata[part],st);

	  if(treelist->tree[part]->n_otu != cdata[part]->n_otu)
	    {
	      PhyML_Printf("\n. Problem with sequence file '%s'\n",io->st->optionlist[part]->in_align_file);
	      PhyML_Printf("\n. # taxa found in supertree restricted to '%s' taxa = %d\n",
		     io->st->optionlist[part]->in_align_file,
		     treelist->tree[part]->n_otu);
	      PhyML_Printf("\n. # sequences in '%s' = %d\n",
		     io->st->optionlist[part]->in_align_file,
		     cdata[part]->n_otu);
	      Exit("\n");
	    }

	  treelist->tree[part]->dp         = part;
	  treelist->tree[part]->n_otu      = cdata[part]->n_otu;
	  treelist->tree[part]->mod        = mod[part];
	  treelist->tree[part]->io         = io->st->optionlist[part];
	  treelist->tree[part]->data       = cdata[part];
	  treelist->tree[part]->n_pattern  = treelist->tree[part]->data->crunch_len/
	                                     treelist->tree[part]->io->state_len;

          Set_Both_Sides(YES,treelist->tree[part]);
	  Connect_CSeqs_To_Nodes(cdata[part],treelist->tree[part]);
	  Fill_Dir_Table(treelist->tree[part]);
	  Update_Dirs(treelist->tree[part]);
	  Make_Tree_4_Lk(treelist->tree[part],cdata[part],cdata[part]->init_len);
	  Make_Tree_4_Pars(treelist->tree[part],cdata[part],cdata[part]->init_len);
	  treelist->tree[part]->triplet_struct = Make_Triplet_Struct(treelist->tree[part]->mod);
          Init_Triplet_Struct(treelist->tree[part]->triplet_struct);
	}

      if(part != io->n_part)
	{
	  PhyML_Printf("\n. Sequence data part found in '%s' has one or more taxa not found in the '%s' tree file\n",
		 io->st->optionlist[part]->in_align_file,
		 io->in_tree_file);
	  Exit("\n");
	}

      PART_Check_Extra_Taxa(st);
	
      st->tree->c_lnL = .0;
      For(part,io->n_part)
	{
	  Lk(NULL,treelist->tree[part]);
/* 	  Get_List_Of_Reachable_Tips(treelist->tree[part]); */
	  PART_Match_St_Nodes_In_Gt(treelist->tree[part],st);
	  PART_Match_St_Edges_In_Gt(treelist->tree[part],st);
	  PART_Map_St_Nodes_In_Gt(treelist->tree[part],st);
	  PART_Map_St_Edges_In_Gt(treelist->tree[part],st);
	  PART_Map_Gt_Edges_In_St(treelist->tree[part],st);
	  st->tree->c_lnL += treelist->tree[part]->c_lnL;
	}

      PART_Initialise_Bl_Partition(st);
      PART_Set_Bl(st->bl,st);

      time(&(st->tree->t_beg));
      time(&(st->tree->t_current));
      PhyML_Printf("\n. (%5d sec) [00] [%10.2f] [%5d]\n",
	     (int)(st->tree->t_current-st->tree->t_beg),
	     PART_Lk(st),
	     PART_Pars(st));


      int n_iter=0;
      do
	{
	  PART_Optimize_Br_Len_Serie(st->tree->a_nodes[0],
				   st->tree->a_nodes[0]->v[0],
				   st->tree->a_nodes[0]->b[0],
				   st);

          Set_Both_Sides(YES,st->tree);
	  PART_Lk(st);
	  PhyML_Printf("\n. %f",st->tree->c_lnL);
/* 	  For(part,st->n_part) PhyML_Printf("\n. %s",Write_Tree(st->treelist->tree[part],NO)); */
	  n_iter++;
	}while(n_iter < 5);
/*       Exit("\n"); */


/*       PART_Lk(st); */
/*       PhyML_Printf("\n> %f",st->tree->c_lnL); */
/*       For(i,2*st->tree->n_otu-3) */
/* 	{ */
/* 	  if((!st->tree->a_edges[i]->left->tax) && (!st->tree->a_edges[i]->rght->tax)) */
/* 	    { */
/* 	      PART_NNI(st->tree->a_edges[i],st); */
/* 	    } */
/* 	} */

      
/*       if(io->mod->s_opt->topo_search == NNI_MOVE) */
      PART_Simu(st);
/*       else */
/*       PART_Speed_Spr(st); */


      time(&t_end);

      hour = div(t_end-t_beg,3600);
      min  = div(t_end-t_beg,60  );

      min.quot -= hour.quot*60;


      PhyML_Fprintf(fp_phyml_stats,"\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n");
      PhyML_Fprintf(fp_phyml_stats,"\n. Number of partitions = %d\n\n",st->n_part);
      PhyML_Fprintf(fp_phyml_stats,"\n. Full data set -- lnL = %f\n\n",st->tree->c_lnL);
      PhyML_Fprintf(fp_phyml_stats,"\n. Tree search algorithm : %s\n\n",(io->mod->s_opt->topo_search == NNI_MOVE)?("NNIs"):("SPRs"));
      PhyML_Fprintf(fp_phyml_stats,"\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n\n");
      For(part,io->n_part)
	{
	  Print_Fp_Out(fp_phyml_stats,t_beg,t_end,st->treelist->tree[part],
		       io->st->optionlist[part],1,1,YES);
	}


      PhyML_Printf("\n\n. Time used %dh%dm%ds\n", hour.quot,min.quot,(int)(t_end-t_beg)%60);
      PhyML_Printf("\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n");

      For(i,2*st->tree->n_otu-3) st->tree->a_edges[i]->l->v = 0.1;
      s_tree = Write_Tree(st->tree,NO);
      PhyML_Fprintf(fp_phyml_tree,"Supertree\n");
      PhyML_Fprintf(fp_phyml_tree,"%s\n",s_tree);
      Free(s_tree);
      For(part,st->n_part)
	{
	  PhyML_Fprintf(fp_phyml_tree,"Gene tree number %d\n",part+1);
	  s_tree = Write_Tree(st->treelist->tree[part],NO);
	  PhyML_Fprintf(fp_phyml_tree,"%s\n",s_tree);
	  Free(s_tree);
	}


      For(part,st->n_part)
	{
	  if(io->mod->s_opt->topo_search == SPR_MOVE) Free_Spr_List(treelist->tree[part]);
	  Free_Tree_Lk(treelist->tree[part]);
	  Free_Tree_Pars(treelist->tree[part]);
	  Free_Triplet(treelist->tree[part]->triplet_struct);
	  Free_Tree(treelist->tree[part]);
	  Free_Cseq(cdata[part]);
	  Free_Model(mod[part]);
	  Free_Input(io->st->optionlist[part]);

	}

      if(io->mod->s_opt->topo_search == SPR_MOVE) Free_Spr_List(st->tree);
      Free(mat);
      Free(mod);
      Free(data);
      Free(cdata);
      Free(treelist->tree);
      Free(treelist);
      Free_St(st);
    }


  if(io->fp_in_align ) fclose(io->fp_in_align);
  if(io->fp_in_tree) fclose(io->fp_in_tree);


  Free_Model(io->mod);
  Free_Input(io);

  fclose(fp_phyml_lk);
  fclose(fp_phyml_tree);
  fclose(fp_phyml_stats);


  return 0;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void PART_Print_Nodes(t_node *a, t_node *d, supert_tree *st)
{
  int i;
  PhyML_Printf(">>>>>>>>>>>>>>>>>>>>\n");
  PhyML_Printf("num t_node = %d\n",d->num);
  if(d->tax) PhyML_Printf("name='%s'\n",d->name);
  else
    {
      PhyML_Printf("n_of_reachable_tips : \n");
      For(i,3)
	{
/* 	  PhyML_Printf("dir%d=%d; ",i,st->n_of_reachable_tips[st->num_data_of_interest][d->num][i]); */
/* 	  For(j,st->n_of_reachable_tips[st->num_data_of_interest][d->num][i]) */
/* 	    { */
/* 	      PhyML_Printf("%s ", */
/* 		     st->list_of_reachable_tips[st->num_data_of_interest][d->num][i][j]->name); */
/* 	    } */
	  
	  PhyML_Printf("\n");
	}
    }
  PhyML_Printf("<<<<<<<<<<<<<<<<<<<\n\n");
  if(d->tax) return;
  else
    {
      For(i,3) if(d->v[i] != a) PART_Print_Nodes(d,d->v[i],st);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


supert_tree *PART_Make_Supert_tree_Light(option *input)
{
  supert_tree *st;

  st = (supert_tree *)mCalloc(1,sizeof(supert_tree));
  st->optionlist = (option **)mCalloc(input->n_part,sizeof(option *));
  st->bl_partition = (int *)mCalloc(input->n_part,sizeof(int ));
  return st;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void PART_Make_Supert_tree_Full(supert_tree *st, option *io, calign **data)
{
  int i,j,k;

  if(io->in_tree == 2)
    {
      PhyML_Printf("\n. Reading user tree...\n");
      rewind(io->fp_in_tree);
      
      st->tree = Read_Tree_File_Phylip(io->fp_in_tree);
      
      if(!st->tree->has_branch_lengths)
	{
	  PhyML_Printf("\n. Branch lengths are all set to 0.1...\n");
	  For(i,2*st->tree->n_otu-3) st->tree->a_edges[i]->l->v = 0.1;
	}
    }
  else
    {
      Warn_And_Exit("\n. A user-defined input tree is needed\n");
    }

  st->tree->io      = io;
  st->treelist      = (t_treelist *)Make_Treelist(io->n_part);
  st->n_part          = io->n_part;
  st->tree->mod     = io->mod;
  st->lock_br_len   = 0;

  st->map_st_node_in_gt = (t_node *****)mCalloc(st->n_part,sizeof(t_node ****));
  For(i,st->n_part) 
    {
      st->map_st_node_in_gt[i] = (t_node ****)mCalloc(2*st->tree->n_otu-2,sizeof(t_node ***));
      For(j,2*st->tree->n_otu-2) 
	{
	  st->map_st_node_in_gt[i][j] = (t_node ***)mCalloc(3,sizeof(t_node **));
	  For(k,3) st->map_st_node_in_gt[i][j][k] = (t_node **)mCalloc(2,sizeof(t_node *));
	}
    }

  st->map_st_edge_in_gt = (t_edge ***)mCalloc(st->n_part,sizeof(t_edge **));
  For(i,st->n_part) st->map_st_edge_in_gt[i] = (t_edge **)mCalloc(2*st->tree->n_otu-3,sizeof(t_edge *));

  st->map_gt_edge_in_st = (t_edge ****)mCalloc(st->n_part,sizeof(t_edge ***));
  For(i,st->n_part)
    {
      st->map_gt_edge_in_st[i] = (t_edge ***)mCalloc(2*st->tree->n_otu-3,sizeof(t_edge **));
      For(j,2*st->tree->n_otu-3) st->map_gt_edge_in_st[i][j] = (t_edge **)mCalloc(2*st->tree->n_otu-3,sizeof(t_edge *));
    }

  st->size_map_gt_edge_in_st = (int **)mCalloc(st->n_part,sizeof(int *));
  For(i,st->n_part) st->size_map_gt_edge_in_st[i] = (int *)mCalloc(2*st->tree->n_otu-3,sizeof(int));


  st->match_st_edge_in_gt = (t_edge ***)mCalloc(st->n_part,sizeof(t_edge **));
  For(i,st->n_part) st->match_st_edge_in_gt[i] = (t_edge **)mCalloc(2*st->tree->n_otu-3,sizeof(t_edge *));

  st->match_gt_edge_in_st = (t_edge ***)mCalloc(st->n_part,sizeof(t_edge **));
  For(i,st->n_part) st->match_gt_edge_in_st[i] = (t_edge **)mCalloc(2*st->tree->n_otu-3,sizeof(t_edge *));

  st->bl = (phydbl **)mCalloc(st->n_part,sizeof(phydbl *));
  For(i,st->n_part) st->bl[i] = (phydbl *)mCalloc(2*st->tree->n_otu-3,sizeof(phydbl));

  st->bl_cpy = (phydbl **)mCalloc(st->n_part,sizeof(phydbl *));
  For(i,st->n_part) st->bl_cpy[i] = (phydbl *)mCalloc(2*st->tree->n_otu-3,sizeof(phydbl));

  st->bl0 = (phydbl **)mCalloc(st->n_part,sizeof(phydbl *));
  For(i,st->n_part) st->bl0[i] = (phydbl *)mCalloc(2*st->tree->n_otu-3,sizeof(phydbl));

  st->bl1 = (phydbl **)mCalloc(st->n_part,sizeof(phydbl *));
  For(i,st->n_part) st->bl1[i] = (phydbl *)mCalloc(2*st->tree->n_otu-3,sizeof(phydbl));

  st->bl2 = (phydbl **)mCalloc(st->n_part,sizeof(phydbl *));
  For(i,st->n_part) st->bl2[i] = (phydbl *)mCalloc(2*st->tree->n_otu-3,sizeof(phydbl));

  st->s_mod = (t_mod **)mCalloc(st->n_part,sizeof(t_mod *));

  For(i,2*st->tree->n_otu-3) Make_Edge_NNI(st->tree->a_edges[i]);

  st->match_st_node_in_gt = (t_node ***)mCalloc(io->n_part,sizeof(t_node **));
  For(i,io->n_part) st->match_st_node_in_gt[i] = (t_node **)mCalloc(2*st->tree->n_otu-2,sizeof(t_node *));
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void PART_Prune_St_Topo(t_tree *tree, calign *data, supert_tree *st)
{
  int i,j,not_found;
  int curr_ext_node, curr_int_node, curr_br, n_pruned_nodes;;
  t_node **pruned_nodes;
  t_edge **residual_edges;

  pruned_nodes   = (t_node **)mCalloc(st->tree->n_otu,sizeof(t_node *));
  residual_edges = (t_edge **)mCalloc(st->tree->n_otu,sizeof(t_edge *));

  n_pruned_nodes = 0;
  For(i,st->tree->n_otu)
    {
      For(j,data->n_otu)
	{
	  if(!strcmp(data->c_seq[j]->name,st->tree->a_nodes[i]->name))
	    break;
	}

      not_found = 1;
      if(j == data->n_otu)
	{
	  For(j,tree->n_otu)
	    {
	      if(!strcmp(tree->a_nodes[j]->name,st->tree->a_nodes[i]->name))
		{
		  Prune_Subtree(tree->a_nodes[j]->v[0],
				tree->a_nodes[j],
				NULL,&(residual_edges[n_pruned_nodes]),
				tree);

		  pruned_nodes[n_pruned_nodes] = tree->a_nodes[j];
		  n_pruned_nodes++;
		  not_found = 0;
		  break;
		}	      
	    }


	  if(not_found)	    
	    {
	      PhyML_Printf("\n. Taxon '%s'",st->tree->a_nodes[i]->name);
	      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
	      Warn_And_Exit("");
	    }
	}
    }

  Free(tree->t_dir);

  tree->n_otu -= n_pruned_nodes;

  curr_ext_node = 0;
  curr_int_node = tree->n_otu;  
  curr_br = 0;
  For(i,st->tree->n_otu)
    {
      For(j,n_pruned_nodes)
	{
	  if(!strcmp(pruned_nodes[j]->name,st->tree->a_nodes[i]->name))
	    break;
	}
      if(j == n_pruned_nodes) /* That t_node still belongs to the tree */
	{
	  Reassign_Node_Nums(tree->a_nodes[i],tree->a_nodes[i]->v[0], 
			     &curr_ext_node, &curr_int_node,tree);
	  break;
	}
    }
  
  Reassign_Edge_Nums(tree->a_nodes[0],tree->a_nodes[0]->v[0],&curr_br,tree);

  tree->t_dir = (short int *)mCalloc((2*tree->n_otu-2)*(2*tree->n_otu-2),sizeof(short int));

  For(i,n_pruned_nodes) 
    {
      Free_Edge(residual_edges[i]);
      Free_Edge(pruned_nodes[i]->b[0]);
      Free_Node(pruned_nodes[i]->v[0]);
      Free_Node(pruned_nodes[i]);
    }

  Free(pruned_nodes);
  Free(residual_edges);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void PART_Match_St_Nodes_In_Gt(t_tree *gt, supert_tree *st)
{
  int i,j;

  For(i,2*st->tree->n_otu-2) st->match_st_node_in_gt[gt->dp][i] = NULL; /* don't forget that step ! */

  /* Map tips */
  For(i,st->tree->n_otu)
    {
      For(j,gt->n_otu)
	{
	  if(!strcmp(st->tree->a_nodes[i]->name,gt->a_nodes[j]->name))
	    {
	      st->match_st_node_in_gt[gt->dp][st->tree->a_nodes[i]->num] = gt->a_nodes[j];
	      break;
	    }
	}
    }

#ifdef DEBUG
  /* Checking that the results are correct so far */
  int n_matches;
  n_matches = 0;
  For(i,2*st->tree->n_otu-2)
    if(st->match_st_node_in_gt[gt->dp][i])
      n_matches++;

  if(n_matches != gt->n_otu)
    {
      PhyML_Printf("\n");
      PhyML_Printf("\n. n_matches = %d 2*gt->n_otu-2 = %d\n",n_matches,2*gt->n_otu-2);
      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
#endif


  /* Map internal nodes */
  For(i,st->tree->n_otu)
    {
      if(st->match_st_node_in_gt[gt->dp][st->tree->a_nodes[i]->num])
	{
	  PART_Match_St_Nodes_In_Gt_Recurr(st->match_st_node_in_gt[gt->dp][st->tree->a_nodes[i]->num],
					 st->match_st_node_in_gt[gt->dp][st->tree->a_nodes[i]->num]->v[0],
					 st->tree->a_nodes[i],
					 st->tree->a_nodes[i]->v[0],
					 gt,
					 st);
	  break;
	}
    }
  


#ifdef DEBUG
  /* Checking that the results are correct */
  n_matches = 0;
  For(i,2*st->tree->n_otu-2) 
    if(st->match_st_node_in_gt[gt->dp][st->tree->a_nodes[i]->num])
	n_matches++;

  if(n_matches != 2*gt->n_otu-2)
    {
      int j;
      PhyML_Printf("\n");
      PhyML_Printf("\n. n_matches = %d 2*gt->n_otu-2 = %d\n",n_matches,2*gt->n_otu-2);
      For(j,2*gt->n_otu-2)
	{
	  For(i,2*st->tree->n_otu-2) 
	    if(st->match_st_node_in_gt[gt->dp][i] == gt->a_nodes[j])
	      break;

 	  if(i == 2*st->tree->n_otu-2)
	    {
	      PhyML_Printf("\n. Gt %3d t_node %3d (%3d %3d %3d) (%s %s %s) (%f %f %f) does not match\n",
		     gt->dp,
		     gt->a_nodes[j]->num,
		     gt->a_nodes[j]->v[0] ? gt->a_nodes[j]->v[0]->num : -1,
		     gt->a_nodes[j]->v[1] ? gt->a_nodes[j]->v[1]->num : -1,
		     gt->a_nodes[j]->v[2] ? gt->a_nodes[j]->v[2]->num : -1,
		     gt->a_nodes[j]->v[0]->tax ? gt->a_nodes[j]->v[0]->name : NULL,
		     gt->a_nodes[j]->v[1]->tax ? gt->a_nodes[j]->v[1]->name : NULL,
		     gt->a_nodes[j]->v[2]->tax ? gt->a_nodes[j]->v[2]->name : NULL,
		     gt->a_nodes[j]->v[0] ? gt->a_nodes[j]->b[0]->l->v : -1.,
		     gt->a_nodes[j]->v[1] ? gt->a_nodes[j]->b[1]->l->v : -1.,
		     gt->a_nodes[j]->v[2] ? gt->a_nodes[j]->b[2]->l->v : -1.);
	    }
	}

      PhyML_Printf("oooooooo\n");
      Print_Node(st->tree->a_nodes[0],
		 st->tree->a_nodes[0]->v[0],
		 st->tree);
      PhyML_Printf(">>>>>>>\n");
      For(i,st->n_part)
	{
	  Print_Node(st->treelist->tree[i]->a_nodes[0],
		     st->treelist->tree[i]->a_nodes[0]->v[0],
		     st->treelist->tree[i]);
	  PhyML_Printf("<<<<<<<\n");
	}
      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
#endif
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void PART_Match_St_Nodes_In_Gt_Recurr(t_node *a_gt, t_node *d_gt, t_node *a_st, t_node *d_st, t_tree *gt, supert_tree *st)
{
  int i,j,k;
  int *score_d_st;


  if((d_gt->tax) || (d_st->tax)) return;
  else
    {
      score_d_st = (int *)mCalloc(3,sizeof(int));
      
      /* Might be wrong. Check function Match_Nodes_In_Small_Tree */
      For(i,3)
	{
	  For(j,3)
	    {
	      For(k,d_st->bip_size[j])
		{
		  if(!strcmp(d_gt->bip_node[i][0]->name,d_st->bip_node[j][k]->name))
		    {
		      score_d_st[j] += 1;
		      break;
		    }
		}
	    }
	}


      if((score_d_st[0] == 2) && (score_d_st[1] == 2) && (score_d_st[2] == 2))
	{
	  st->match_st_node_in_gt[gt->dp][d_st->num] = d_gt;

	  For(i,3)
	    {
	      if(d_gt->v[i] != a_gt)
		{
		  For(j,3)
		    {		      

		      if(score_d_st[j] != 3)
			{
			  PART_Match_St_Nodes_In_Gt_Recurr(d_gt,d_gt->v[i],d_st,d_st->v[j],gt,st);
			  break;
			}

/* 		      For(k,d_st->n_of_reachable_tips[j]) */
/* 			if(!strcmp(d_gt->list_of_reachable_tips[i][0]->name, */
/* 				   d_st->list_of_reachable_tips[j][k]->name)) */
/* 			  { */
/* 			    PART_Match_St_Nodes_In_Gt_Recurr(d_gt,d_gt->v[i],d_st,d_st->v[j],gt,st); */
/* 			    break; */
/* 			  } */
/* 		      if(k != d_st->n_of_reachable_tips[j]) break; */
		    }
		}
	    }
	}
      else
	{
	  For(i,3)
	    if(d_st->v[i] != a_st)
	      PART_Match_St_Nodes_In_Gt_Recurr(a_gt,d_gt,d_st,d_st->v[i],gt,st);
	}
      Free(score_d_st);	
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void PART_Match_St_Edges_In_Gt(t_tree *gt, supert_tree *st)
{
  int i;

  For(i,2*st->tree->n_otu-3) 
    {
      st->match_st_edge_in_gt[gt->dp][i] = NULL; 
      st->match_gt_edge_in_st[gt->dp][i] = NULL;
    }

  For(i,st->tree->n_otu) 
    if(st->match_st_node_in_gt[gt->dp][i])
      {
	PART_Match_St_Edges_In_Gt_Recurr(st->match_st_node_in_gt[gt->dp][i],
				       st->match_st_node_in_gt[gt->dp][i]->v[0],
				       st->tree->a_nodes[i],
				       st->tree->a_nodes[i]->v[0],
				       gt,st);
	break;
      }

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void PART_Match_St_Edges_In_Gt_Recurr(t_node *a_gt, t_node *d_gt, t_node *a_st, t_node *d_st, t_tree *gt, supert_tree *st)
{
  t_edge *b_gt, *b_st;
  int i,j;

  b_gt = b_st = NULL;

  if((st->match_st_node_in_gt[gt->dp][a_st->num] == a_gt) &&
     (st->match_st_node_in_gt[gt->dp][d_st->num] == d_gt))
    {
      For(i,3) if((a_st->v[i]) && (a_st->v[i] == d_st)) {b_st = a_st->b[i]; break;}
      For(i,3) if((a_gt->v[i]) && (a_gt->v[i] == d_gt)) {b_gt = a_gt->b[i]; break;}

      st->match_st_edge_in_gt[gt->dp][b_st->num] = b_gt;
      st->match_gt_edge_in_st[gt->dp][b_gt->num] = b_st;
    }


  if(!d_gt)
    {
      PhyML_Printf("\n");
      PhyML_Printf("\n. a_gt->num = %d\n",a_gt->num);
      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  if(d_gt->tax || d_st->tax) return;
  else
    {
      if(st->match_st_node_in_gt[gt->dp][d_st->num] == d_gt)
	{
	  For(i,3)
	    {
	      if(d_gt->v[i] != a_gt)
		{
		  For(j,3)
		    {
/* 		      For(k,d_st->n_of_reachable_tips[j]) */
/* 			if(!strcmp(d_gt->list_of_reachable_tips[i][0]->name,d_st->list_of_reachable_tips[j][k]->name)) */
			  {
			    PART_Match_St_Edges_In_Gt_Recurr(d_gt,d_gt->v[i],d_st,d_st->v[j],gt,st);
			    break;
			  }
/* 		      if(k != d_st->n_of_reachable_tips[j]) break; */
		    }
		}
	    }
	}
      else
	{
	  For(i,3)
	    if(d_st->v[i] != a_st)
	      PART_Match_St_Edges_In_Gt_Recurr(a_gt,d_gt,d_st,d_st->v[i],gt,st);
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void PART_Simu(supert_tree *st)
{
  int i,j,step,n_without_swap,it_lim_without_swap;
  t_edge **sorted_b,*st_b,**tested_b;
  int n_neg,n_tested,each;
  phydbl lambda,old_loglk;

  sorted_b = (t_edge **)mCalloc(st->tree->n_otu-3,sizeof(t_edge *));
  tested_b = (t_edge **)mCalloc(st->tree->n_otu-3,sizeof(t_edge *));

  For(i,st->n_part) Update_Dirs(st->treelist->tree[i]);
  Update_Dirs(st->tree);

  each                = 4;
  step                = 0;
  lambda              = .75;
  n_tested            = 0;
  n_without_swap      = 0;
  old_loglk           = UNLIKELY; 
  it_lim_without_swap = 2;
  st_b                = NULL; 
  do
    {
      For(i,st->n_part) Check_Dirs(st->treelist->tree[i]);

      ++step;     
      each--;

/*       PART_Print_Bl(st); */

      /* Compute the likelihood of the supertreee */
      st->tree->c_lnL  = PART_Lk(st);
      st->tree->c_pars = PART_Pars(st);
/*       For(i,st->n_part) PhyML_Printf("\n. %s",Write_Tree(st->treelist->tree[i],NO)); */
/*       PhyML_Printf("\n"); */

      time(&(st->tree->t_current));
      PhyML_Printf("\n. (%5d sec) [tot lnL=%15.5f] [# swaps=%3d]",
	     (int)(st->tree->t_current-st->tree->t_beg),
	     st->tree->c_lnL,n_tested);
/*       For(i,st->n_part) PhyML_Printf("\n[gt %3d lnL=%15.5f]",i,st->treelist->tree[i]->c_lnL); */
      
      if((FABS(old_loglk-st->tree->c_lnL) < st->tree->mod->s_opt->min_diff_lk_global) || 
	 (n_without_swap > it_lim_without_swap)) break;

      if(st->tree->c_lnL < old_loglk)
	{
	  PhyML_Printf("\n. Moving backward (topology + branch lengths) \n");
	  
	  if(!PART_Mov_Backward_Topo_Bl(st,old_loglk,tested_b,n_tested))
	    Warn_And_Exit("\n. Err: mov_back failed\n");

	  if(!st->tree->n_swap) n_neg = 0;
	    
	  PART_Record_Br_Len(st);
	  For(i,st->n_part) Optimiz_All_Free_Param(st->treelist->tree[i],0);
	}
      else 
	{
	  if(!each)
	    {
	      each = 4;
	      /* Markov model parameters are free to vary across data partitions */
	      For(i,st->n_part) Optimiz_All_Free_Param(st->treelist->tree[i],0);	      
	      For(i,st->n_part) Set_Both_Sides(YES,st->treelist->tree[i]);
	      st->tree->c_lnL  = PART_Lk(st);
	      st->tree->c_pars = PART_Pars(st);
	    }
	  
	  old_loglk = st->tree->c_lnL;
	  

	  For(i,2*st->tree->n_otu-3) Init_NNI(st->tree->a_edges[i]->nni);

	  /* Test NNIs */
	  For(i,2*st->tree->n_otu-3)
	    {
	      st_b = st->tree->a_edges[i];
	      if((!st_b->left->tax) && (!st_b->rght->tax)) PART_NNI(st_b,st);
	    }
	  
	  /* Optimise external branch lengths */
	  For(i,2*st->tree->n_otu-3)
	    {
	      st_b = st->tree->a_edges[i];
	      if((st_b->left->tax) || (st_b->rght->tax))
		{
		  PART_Record_Br_Len(st);
		  PART_Br_Len_Brent(st_b,0,st);
		  For(j,st->n_part) st->bl0[st->bl_partition[j]][st_b->num] = st->bl[st->bl_partition[j]][st_b->num];
		  st_b->nni->score     = .0;
		  st_b->nni->best_conf =  0;
		  PART_Restore_Br_Len(st);
		  PART_Lk_At_Given_Edge(st_b,st);
		}
	    }

	  /* Select and sort swaps */
	  n_neg = 0;
	  Select_Edges_To_Swap(st->tree,sorted_b,&n_neg);
	  Sort_Edges_NNI_Score(st->tree,sorted_b,n_neg);	  

	  n_tested = 0;
	  For(i,(int)CEIL((phydbl)n_neg*(lambda)))
	    tested_b[n_tested++] = sorted_b[i];

	  if(n_tested > 0) n_without_swap = 0;
	  else             n_without_swap++;

	  PART_Record_Br_Len(st);
	  
	  /* Apply swaps */
	  PART_Make_N_Swap(tested_b,0,n_tested,st);

	  /* Update branch lengths (all edges first and then swaped edges) */
	  PART_Update_Bl(lambda,st);
	  PART_Update_Bl_Swaped(tested_b,n_tested,st);

	    
	}
    }
  while(1);
  
  PhyML_Printf("\n\n. End of PART_Simu \n");
  Free(sorted_b);
  Free(tested_b);

  if((n_without_swap > it_lim_without_swap))
    {
      PhyML_Printf("\n. Last optimization step...\n");
      For(i,st->n_part) Round_Optimize(st->treelist->tree[i],st->treelist->tree[i]->data,ROUND_MAX);
      st->tree->c_lnL = PART_Lk(st);
      PART_Simu(st);
    }
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


int PART_Mov_Backward_Topo_Bl(supert_tree *st, phydbl lk_old, t_edge **tested_b, int n_tested)
{
  int i,j,step,beg,end;
  t_edge *st_b;
  phydbl **l_init;
  int dim;

  l_init = (phydbl **)mCalloc(st->n_part,sizeof(phydbl *));
  For(i,st->n_part) l_init[i] = (phydbl *)mCalloc(2*st->tree->n_otu-3,sizeof(phydbl ));

  For(i,2*st->tree->n_otu-3) For(j,st->n_part) l_init[st->bl_partition[j]][i] = st->bl[st->bl_partition[j]][i];

  step = 2;
  do
    {
      For(i,2*st->tree->n_otu-3)
	{
	  For(j,st->n_part)
	    {
	      st->bl[st->bl_partition[j]][i] = st->bl_cpy[st->bl_partition[j]][i] + 
		(1./step) * (l_init[st->bl_partition[j]][i] - st->bl_cpy[st->bl_partition[j]][i]);
/* 	      st->bl[st->bl_partition[j]][i] = st->bl_cpy[st->bl_partition[j]][i]; */
	    }
	}
      
      beg = (int)FLOOR((phydbl)n_tested/(step-1));
      end = 0;
      st_b = NULL;
      dim = 2*st->tree->n_otu-2;

      for(i=beg-1;i>=end;i--)
	{
	  st_b = tested_b[i];

	  PART_Swap(st_b->nni->swap_node_v2->v[st->tree->t_dir[st_b->nni->swap_node_v2->num*dim+st_b->nni->swap_node_v1->num]],
		  st_b->nni->swap_node_v2,
		  st_b->nni->swap_node_v3,
		  st_b->nni->swap_node_v3->v[st->tree->t_dir[st_b->nni->swap_node_v3->num*dim+st_b->nni->swap_node_v4->num]],
		  st);

	  Swap(st_b->nni->swap_node_v2->v[st->tree->t_dir[st_b->nni->swap_node_v2->num*dim+st_b->nni->swap_node_v1->num]],
	       st_b->nni->swap_node_v2,
	       st_b->nni->swap_node_v3,
	       st_b->nni->swap_node_v3->v[st->tree->t_dir[st_b->nni->swap_node_v3->num*dim+st_b->nni->swap_node_v4->num]],
	       st->tree);

	  PART_Do_Mapping(st);

	}

      beg = 0;
      end = (int)FLOOR((phydbl)n_tested/step);
      st_b = NULL;
      dim = 2*st->tree->n_otu-2;

      for(i=beg;i<end;i++)
	{
	  st_b = tested_b[i];
	  
	  PART_Swap(st_b->nni->swap_node_v2->v[st->tree->t_dir[st_b->nni->swap_node_v2->num*dim+st_b->nni->swap_node_v1->num]],
		  st_b->nni->swap_node_v2,
		  st_b->nni->swap_node_v3,
		  st_b->nni->swap_node_v3->v[st->tree->t_dir[st_b->nni->swap_node_v3->num*dim+st_b->nni->swap_node_v4->num]],
		  st);

	  Swap(st_b->nni->swap_node_v2->v[st->tree->t_dir[st_b->nni->swap_node_v2->num*dim+st_b->nni->swap_node_v1->num]],
	       st_b->nni->swap_node_v2,
	       st_b->nni->swap_node_v3,
	       st_b->nni->swap_node_v3->v[st->tree->t_dir[st_b->nni->swap_node_v3->num*dim+st_b->nni->swap_node_v4->num]],
	       st->tree);

	  PART_Do_Mapping(st);

	}
     
      if(!end) st->tree->n_swap = 0;      

      PART_Lk(st);

      PhyML_Printf("\n. lnL = %15.5f",st->tree->c_lnL);
      step++;
    }
  while((st->tree->c_lnL < lk_old) && (step < 100));

  if(step == 100)
    {
      For(i,2*st->tree->n_otu-3) For(j,st->n_part)
	st->bl[st->bl_partition[j]][i] = st->bl_cpy[st->bl_partition[j]][i];
    }

  st->tree->n_swap = 0;
  For(i,2*st->tree->n_otu-3) 
    {
      if(st->tree->a_edges[i]->nni->score < 0.0) st->tree->n_swap++;
      st->tree->a_edges[i]->nni->score = +1.0;
    }

  PART_Lk(st);

  if(st->tree->c_lnL > lk_old)                                                    return  1;
  else if(FABS(st->tree->c_lnL-lk_old) < st->tree->mod->s_opt->min_diff_lk_local) return -1;
  else                                                                            return  0;

  For(i,st->n_part) Free(l_init[i]);
  Free(l_init);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void PART_Check_Extra_Taxa(supert_tree *st)
{
  int i,j,k;
  int sum;
  int *st_taxa;

  st_taxa = (int *)mCalloc(st->tree->n_otu,sizeof(int));

  For(i,st->tree->n_otu)
    {
      For(j,st->n_part)
	{
	  For(k,st->treelist->tree[j]->n_otu) 
	    if(!strcmp(st->treelist->tree[j]->a_nodes[k]->name,st->tree->a_nodes[i]->name)) break;
	  if(k != st->treelist->tree[j]->n_otu) { st_taxa[i] = 1; break; }
	}
    }

  sum = 0;
  For(i,st->tree->n_otu) if(st_taxa[i]) sum++;
  if(sum != st->tree->n_otu) 
    {
      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
  Free(st_taxa);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


int PART_Get_Species_Found_In_St(supert_tree *st, calign *data)
{
  int i,j;
  
  For(i,data->n_otu)
    {
      For(j,st->tree->n_otu)
	{
	  if(!strcmp(data->c_seq[i]->name,st->tree->a_nodes[j]->name))
	    {
	      break;
	    }
	}
      if(j == st->tree->n_otu) return 0;
    }
  return 1;

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void PART_Map_St_Nodes_In_Gt(t_tree *gt, supert_tree *st)
{
  int i;

  For(i,2*st->tree->n_otu-2)
    {
      st->map_st_node_in_gt[gt->dp][i][0][0] = NULL;
      st->map_st_node_in_gt[gt->dp][i][1][0] = NULL;
      st->map_st_node_in_gt[gt->dp][i][2][0] = NULL;

      st->map_st_node_in_gt[gt->dp][i][0][1] = NULL;
      st->map_st_node_in_gt[gt->dp][i][1][1] = NULL;
      st->map_st_node_in_gt[gt->dp][i][2][1] = NULL;
    }

  
  /* Root */
  PART_Map_St_Nodes_In_Gt_One_Edge(st->tree->a_nodes[0]->v[0],
				 st->tree->a_nodes[0],
				 st->tree->a_nodes[0]->b[0],
				 gt,st);

  /* Internal nodes */
  PART_Map_St_Nodes_In_Gt_Post(st->tree->a_nodes[0],st->tree->a_nodes[0]->v[0],gt,st);
  PART_Map_St_Nodes_In_Gt_Pre (st->tree->a_nodes[0],st->tree->a_nodes[0]->v[0],gt,st);
  
  /* Root */
  PART_Map_St_Nodes_In_Gt_One_Edge(st->tree->a_nodes[0],
				 st->tree->a_nodes[0]->v[0],
				 st->tree->a_nodes[0]->b[0],
				 gt,st);
  
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void PART_Map_St_Nodes_In_Gt_Post(t_node *a_st, t_node *d_st, t_tree *gt, supert_tree *st)
{
  int i;

  if(d_st->tax) return;
  else
    {
      For(i,3)
	if(d_st->v[i] != a_st)
	  PART_Map_St_Nodes_In_Gt_Post(d_st,d_st->v[i],gt,st);
      
      For(i,3)
	if(d_st->v[i] != a_st)
	  {
	    PART_Map_St_Nodes_In_Gt_One_Edge(d_st,d_st->v[i],d_st->b[i],gt,st);
	  }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void PART_Map_St_Nodes_In_Gt_Pre(t_node *a_st, t_node *d_st, t_tree *gt, supert_tree *st)
{
  int i;

  if(d_st->tax) return;
  else
    {
      For(i,3)
	{
	  if(d_st->v[i] != a_st)
	    {	      
	      PART_Map_St_Nodes_In_Gt_One_Edge(d_st->v[i],d_st,d_st->b[i],gt,st);
	      PART_Map_St_Nodes_In_Gt_Pre(d_st,d_st->v[i],gt,st);
	    }
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void PART_Map_St_Nodes_In_Gt_One_Edge(t_node *a_st, t_node *d_st, t_edge *b_st, t_tree *gt, supert_tree *st)
{
  if(d_st->tax)
    {
#ifdef DEBUG
      if(b_st->rght != d_st)
	{
	  PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}
#endif
      
      st->map_st_node_in_gt[gt->dp][d_st->num][0][0] = st->match_st_node_in_gt[gt->dp][d_st->num];
    }
  else
    {
      t_node **list_of_nodes_d, **list_of_nodes_v1, **list_of_nodes_v2;
      int dir1, dir2;
      int i;

      list_of_nodes_d  = NULL;
      list_of_nodes_v1 = NULL;
      list_of_nodes_v2 = NULL;

      dir1 = dir2 = -1;
      For(i,3) 
	{
	  if(d_st->v[i] != a_st) (dir1 < 0)?(dir1 = i):(dir2 = i);
	  else list_of_nodes_d = st->map_st_node_in_gt[gt->dp][d_st->num][i];
	}

      For(i,3) 
	if((d_st->v[dir1]->v[i]) && (d_st->v[dir1]->v[i] == d_st)) 
	  {
	    list_of_nodes_v1 = st->map_st_node_in_gt[gt->dp][d_st->v[dir1]->num][i];
	    break;
	  }

      For(i,3) 
	if((d_st->v[dir2]->v[i]) && (d_st->v[dir2]->v[i] == d_st)) 
	  {
	    list_of_nodes_v2 = st->map_st_node_in_gt[gt->dp][d_st->v[dir2]->num][i];
	    break;
	  }

      /* d_st matches one t_node in gt */
      if(st->match_st_node_in_gt[gt->dp][d_st->num])
	{
	  list_of_nodes_d[0] = st->match_st_node_in_gt[gt->dp][d_st->num];
	  list_of_nodes_d[1] = NULL;
	}
      else
	{
	  /* list_of_nodes = union of  list_of_nodes_v1  &  list_of_nodes_v2 */
	  
	  if(!list_of_nodes_v1[0])
	    {
	      list_of_nodes_d[0] = list_of_nodes_v2[0];
	      list_of_nodes_d[1] = list_of_nodes_v2[1];
	    }
	  else if(!list_of_nodes_v2[0])
	    {
	      list_of_nodes_d[0] = list_of_nodes_v1[0];
	      list_of_nodes_d[1] = list_of_nodes_v1[1];
	    }
	  else
	    {
	      list_of_nodes_d[0] = list_of_nodes_v1[0];
	      list_of_nodes_d[1] = list_of_nodes_v2[0];

	      if(list_of_nodes_v1[1] || list_of_nodes_v2[1])
		{
		  Print_Node(st->tree->a_nodes[0],
			     st->tree->a_nodes[0]->v[0],
			     st->tree);
		  
		  PhyML_Printf("\n\n--------------------------\n\n");
		  Print_Node(gt->a_nodes[0],
			     gt->a_nodes[0]->v[0],
			     gt);

		  PhyML_Printf("\n\n--------------------------\n\n");
		  
		  PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
		  Warn_And_Exit("");
		}
	    }
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void PART_Map_St_Edges_In_Gt(t_tree *gt, supert_tree *st)
{
  int i,j;
  t_edge *st_b;
  t_node *gt_a, *gt_d;
  
  gt_a = NULL;
  gt_d = NULL;

  For(i,2*st->tree->n_otu-3) st->map_st_edge_in_gt[gt->dp][i] = NULL;

  For(i,2*st->tree->n_otu-3)
    {
      st_b = st->tree->a_edges[i];

      if(!st->map_st_node_in_gt[gt->dp][st_b->left->num][st_b->l_r][0])
	{
	  gt_a = st->map_st_node_in_gt[gt->dp][st_b->rght->num][st_b->r_l][0];
	  gt_d = st->map_st_node_in_gt[gt->dp][st_b->rght->num][st_b->r_l][1];

	  For(j,3)
	    {
	      if((gt_a->v[j]) && (gt_a->v[j] == gt_d))
		{
		  st->map_st_edge_in_gt[gt->dp][st_b->num] = gt_a->b[j];
		  break;
		}
	    }
	}
      else if(!st->map_st_node_in_gt[gt->dp][st_b->rght->num][st_b->r_l][0])
	{
	  gt_a = st->map_st_node_in_gt[gt->dp][st_b->left->num][st_b->l_r][0];
	  gt_d = st->map_st_node_in_gt[gt->dp][st_b->left->num][st_b->l_r][1];

	  For(j,3)
	    {
	      if((gt_a->v[j]) && (gt_a->v[j] == gt_d))
		{
		  st->map_st_edge_in_gt[gt->dp][st_b->num] = gt_a->b[j];
		  break;
		}
	    }
	}
      else
	{	  
	  gt_a = st->map_st_node_in_gt[gt->dp][st_b->left->num][st_b->l_r][0];
	  gt_d = st->map_st_node_in_gt[gt->dp][st_b->rght->num][st_b->r_l][0];

	  #ifdef DEBUG
	  if((!gt_a) || (!gt_d))
	    {
	      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
	      Warn_And_Exit("");
	    }
	  #endif

	  For(j,3)
	    {
	      if((gt_a->v[j]) && (gt_a->v[j] == gt_d))
		{
		  st->map_st_edge_in_gt[gt->dp][st_b->num] = gt_a->b[j];
		  break;
		}
	    }
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void PART_Map_Gt_Edges_In_St(t_tree *gt, supert_tree *st)
{
  int i;
  t_edge *st_b, *gt_b;
  
  For(i,2*st->tree->n_otu-3) st->size_map_gt_edge_in_st[gt->dp][i] = 0;

  st_b = NULL;
  gt_b = NULL;
  For(i,2*st->tree->n_otu-3)
    {
      st_b = st->tree->a_edges[i];
      gt_b = st->map_st_edge_in_gt[gt->dp][st_b->num];

      if(gt_b)
	{
	  st->map_gt_edge_in_st[gt->dp][gt_b->num][st->size_map_gt_edge_in_st[gt->dp][gt_b->num]] = st_b;
	  st->size_map_gt_edge_in_st[gt->dp][gt_b->num]++;
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


int PART_Pars(supert_tree *st)
{
  int i;
  
  st->tree->c_pars = 0;
  For(i,st->n_part) 
    {
      Set_Both_Sides(YES,st->treelist->tree[i]);	  
      Pars(NULL,st->treelist->tree[i]);
      st->tree->c_pars += st->treelist->tree[i]->c_pars;
    }

  return st->tree->c_pars;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


int PART_Spr(phydbl init_lnL, supert_tree *st)
{
  int gt;
  int i;
  t_edge *pruned;
  int best_move;
  t_node *gt_a, *gt_d;

  st->tree->n_root     = st->tree->a_nodes[0];
  pruned               = NULL;
  gt_a                 = NULL;
  gt_d                 = NULL;
  
  Set_Both_Sides(YES,st->tree);

  For(i,2*st->tree->n_otu-3)
    {
      pruned            = st->tree->a_edges[i];
      st->tree->n_moves = 0;

      Reset_Spr_List(st->tree);
      For(gt,st->n_part) Reset_Spr_List(st->treelist->tree[gt]);
      
      if(!pruned->rght->tax)
	{
	  For(gt,st->n_part)
	    {
	      /* Check constraints at prune site on gt tree */
 	      gt_a = st->map_st_node_in_gt[gt][pruned->rght->num][pruned->r_l][0];
 	      gt_d = st->map_st_node_in_gt[gt][pruned->left->num][pruned->l_r][0];

	      if((gt_a) && (gt_d) && (!st->map_st_edge_in_gt[gt][pruned->num]->rght->tax))
		{
		  Test_All_Spr_Targets(st->map_st_edge_in_gt[gt][pruned->num],
				       st->map_st_edge_in_gt[gt][pruned->num]->rght,
				       st->treelist->tree[gt]);
		}
	    }
	}

      if(!pruned->left->tax)
	{
	  For(gt,st->n_part)
	    {
	      /* Check constraints at prune site on gt tree */	      
 	      gt_a = st->map_st_node_in_gt[gt][pruned->rght->num][pruned->r_l][0];
 	      gt_d = st->map_st_node_in_gt[gt][pruned->left->num][pruned->l_r][0];

	      if((gt_a) && (gt_d) && (!st->map_st_edge_in_gt[gt][pruned->num]->left->tax))
		{
		  Test_All_Spr_Targets(st->map_st_edge_in_gt[gt][pruned->num],
				       st->map_st_edge_in_gt[gt][pruned->num]->left,
				       st->treelist->tree[gt]);
		}
	    }
	}

      
      if(!pruned->left->tax)
	{
	  PART_Test_All_Spr_Targets(st->tree->a_edges[i],
				  st->tree->a_edges[i]->left,
				  st);      
	}
      
      if(!pruned->rght->tax)
	{
	  PART_Test_All_Spr_Targets(st->tree->a_edges[i],
				  st->tree->a_edges[i]->rght,
				  st);      
	}


      if(st->tree->n_moves)
	{
	  best_move = PART_Test_List_Of_Regraft_Pos(st->tree->spr_list,
						  (int)CEIL(0.1*(st->tree->n_moves)),
						  st);	  
	  
	  if(st->tree->spr_list[best_move]->lnL > init_lnL)
	    {
	      PART_Try_One_Spr_Move(st->tree->spr_list[best_move],st);
	    }
	  else
	    {
              Set_Both_Sides(YES,st->tree);
	      st->tree->c_lnL  = PART_Lk(st);
	      st->tree->c_pars = PART_Pars(st);	      
	    }
	}
    }
  return 1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void PART_Speed_Spr(supert_tree *st)
{
  int step;
  int gt;
  phydbl old_lnL;

  Make_Spr_List(st->tree);
  For(gt,st->n_part) Make_Spr_List(st->treelist->tree[gt]);

  Set_Both_Sides(YES,st->tree); 
  For(gt,st->n_part) 
    {
      Set_Both_Sides(YES,st->treelist->tree[gt]);
      Record_Br_Len(st->treelist->tree[gt]);
    }
  
  st->tree->c_pars = PART_Pars(st);
  st->tree->c_lnL  = PART_Lk(st);
  

  st->tree->best_lnL = st->tree->c_lnL;
  old_lnL            = st->tree->c_lnL;
  step               = 0;
  do
    {
      ++step;

      PhyML_Printf("\n. Starting a SPR cycle... \n");

      old_lnL = st->tree->c_lnL;

      st->tree->n_improvements         = 0;
      st->tree->perform_spr_right_away = 1;
      PART_Spr(UNLIKELY,st);

 
      time(&(st->tree->t_current));      
      PhyML_Printf("\n. (%5d sec) [00] [%10.2f] [%5d]\n",
	     (int)(st->tree->t_current-st->tree->t_beg),
	     PART_Lk(st),PART_Pars(st));

      /* Optimise parameters of the Markov t_mod */
      For(gt,st->n_part) Optimiz_All_Free_Param(st->treelist->tree[gt],
					      st->treelist->tree[gt]->mod->s_opt->print);

      time(&(st->tree->t_current));      
      PhyML_Printf("\n. (%5d sec) [ 0] [%10.2f] [%5d]\n",
	     (int)(st->tree->t_current-st->tree->t_beg),
	     PART_Lk(st),PART_Pars(st));

      /* Optimise branch lengths */
      For(gt,st->n_part)
	{
	  Optimize_Br_Len_Serie(st->treelist->tree[gt]);
	}


      /* Update partial likelihoods & parsimony */
      Set_Both_Sides(YES,st->tree); 
      st->tree->c_pars = PART_Pars(st);
      st->tree->c_lnL  = PART_Lk(st);
      
      
      time(&(st->tree->t_current));      
      PhyML_Printf("\n. (%5d sec) [**] [%10.2f] [%5d]\n",
	     (int)(st->tree->t_current-st->tree->t_beg),
	     st->tree->c_lnL,st->tree->c_pars);

      /* Record the current best log-likleihood  */
      st->tree->best_lnL = st->tree->c_lnL;

      if(st->tree->c_lnL < old_lnL)
	{
	  PhyML_Printf("\n. old_lnL = %f c_lnL = %f\n",old_lnL,st->tree->c_lnL); 
	  PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}

      /* Record the current best branch lengths  */
      For(gt,st->n_part) Record_Br_Len(st->treelist->tree[gt]);

      /* Exit if no improvements after complete optimization */
      if((!st->tree->n_improvements) || 
	 (FABS(old_lnL-st->tree->c_lnL) < st->tree->mod->s_opt->min_diff_lk_global)) break;
            
    }while(1);
  
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////



void PART_Test_All_Spr_Targets(t_edge *pruned, t_node *n_link, supert_tree *st)
{
  int i,j;

  For(i,3)
    {
      if(n_link->b[i] != pruned)
	{
	  For(j,3)
	    {
	      if((n_link->v[i]->v[j]) && (n_link->v[i]->v[j] != n_link))
		{
		  PART_Test_One_Spr_Target_Recur(n_link->v[i],n_link->v[i]->v[j],n_link->v[i]->b[j],pruned,n_link,st);
		}
	    }
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void PART_Test_One_Spr_Target_Recur(t_node *a, t_node *d, t_edge *target, t_edge *pruned, t_node *n_link, supert_tree *st)
{

  PART_Test_One_Spr_Target(pruned,target,n_link,st);

  if(d->tax) return;
  else
    {
      int i;

      For(i,3)
	if(d->v[i] != a)
	  {
	    PART_Test_One_Spr_Target_Recur(d,d->v[i],d->b[i],pruned,n_link,st);
	  }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void PART_Test_One_Spr_Target(t_edge *st_p, t_edge *st_t, t_node *n_link, supert_tree *st) 
{
  int gt, move;

  st->tree->n_moves++;
  st->tree->spr_list[st->tree->size_spr_list]->b_target      = st_t;
  st->tree->spr_list[st->tree->size_spr_list]->n_link        = n_link;
  st->tree->spr_list[st->tree->size_spr_list]->n_opp_to_link = (n_link == st_p->left)?(st_p->rght):(st_p->left);
  st->tree->spr_list[st->tree->size_spr_list]->b_opp_to_link = st_p;
  st->tree->spr_list[st->tree->size_spr_list]->pars          = 0;

  For(gt,st->n_part)
    {
      move = Map_Spr_Move(st_p,st_t,n_link,st->treelist->tree[gt],st);
      
      if(move > -1)
	st->tree->spr_list[st->tree->size_spr_list]->pars += st->treelist->tree[gt]->spr_list[move]->pars;
      else if(move == -1 || move == -2)
	st->tree->spr_list[st->tree->size_spr_list]->pars += st->treelist->tree[gt]->c_pars;
    }

  Include_One_Spr_To_List_Of_Spr(st->tree->spr_list[st->tree->size_spr_list],st->tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


int Map_Spr_Move(t_edge *st_pruned, t_edge *st_target, t_node *st_link, t_tree *gt, supert_tree *st)
{
  int i;
  t_edge *gt_pruned, *gt_target;
  t_node *gt_link, *gt_a, *gt_d;

  gt_pruned = NULL;
  gt_target = NULL;
  gt_link   = NULL;
  gt_a      = NULL;
  gt_d      = NULL;

  /* Check contraints at prune and regraft sites on the gt tree */

  /* st_pruned is not on a path that matches a branch in gt */
  gt_a = st->map_st_node_in_gt[gt->dp][st_pruned->left->num][st_pruned->l_r][0];
  gt_d = st->map_st_node_in_gt[gt->dp][st_pruned->rght->num][st_pruned->r_l][0];

  if((!gt_a) || (!gt_d)) return -1;
  else
    {
      /* which gt nodes matches st_link ? */
      gt_link = (st_pruned->left == st_link)?(gt_a):(gt_d);
      
      if(gt_link->tax) return -1;
      else
	{
	  gt_pruned = st->map_st_edge_in_gt[gt->dp][st_pruned->num];
	  gt_target = st->map_st_edge_in_gt[gt->dp][st_target->num];
	   
	  if((gt_pruned->left == gt_target->left) ||
	     (gt_pruned->left == gt_target->rght) ||
	     (gt_pruned->rght == gt_target->left) ||
	     (gt_pruned->rght == gt_target->rght)) return -1;
	  else
	    {
	      For(i,gt->size_spr_list)
		{
		  if((gt_pruned == gt->spr_list[i]->b_opp_to_link) && (gt_target == gt->spr_list[i]->b_target))
		    return i;
		}
	    }
	}
    }
  return -2;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


int PART_Test_List_Of_Regraft_Pos(t_spr **st_spr_list, int list_size, supert_tree *st)
{

  int i,j,best_move;
  t_spr *move;
  t_edge *init_target, *b_residual;
  phydbl best_lnL, init_lnL;
  int dir_v0, dir_v1, dir_v2;
  int gt;
  int move_num;
  

  best_lnL = UNLIKELY;
  init_target = b_residual = NULL;
  best_move = -1;

#ifdef DEBUG
  if(!list_size)
    {
      PhyML_Printf("\n\n. List size is 0 !");
      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
      Warn_And_Exit(""); 
    }
#endif
  
  init_lnL = UNLIKELY;

  For(i,list_size)
    {
      st->tree->spr_list[i]->lnL = .0;

      For(gt,st->n_part) 
	{
	  move_num = Map_Spr_Move(st->tree->spr_list[i]->b_opp_to_link,
				  st->tree->spr_list[i]->b_target,
				  st->tree->spr_list[i]->n_link,
				  st->treelist->tree[gt],st);

	  if(move_num > -1)
	    {
	      move = st->treelist->tree[gt]->spr_list[move_num];

	      if(move->b_target)
		{
		  init_lnL = st->treelist->tree[gt]->c_lnL;

		  /* Record t_edge lengths */
		  Record_Br_Len(st->treelist->tree[gt]);
		  
		  /* Prune subtree */
		  Prune_Subtree(move->n_link,
				move->n_opp_to_link,			    
				&init_target,
				&b_residual,
				st->treelist->tree[gt]);

		  /* Rough optimisation of the branch length */
		  Fast_Br_Len(init_target,st->treelist->tree[gt],0);
		  
		  /* Update the change proba matrix at prune position */
		  Update_PMat_At_Given_Edge(init_target,st->treelist->tree[gt]);
	      
		  /* Update partial likelihood along the path from the prune to
		     the regraft position */
		  Update_P_Lk_Along_A_Path(move->path,move->depth_path,st->treelist->tree[gt]);

		  /* Regraft subtree */
		  Graft_Subtree(move->b_target,move->n_link,b_residual,st->treelist->tree[gt]);
	      
		  /* Estimate the three t_edge lengths at the regraft site */
		  Triple_Dist(move->n_link,st->treelist->tree[gt],-1);
	      
		  /* Update the transition proba matrices along edges defining 
		     the regraft site */
		  For(j,3)
		    if(move->n_link->v[j] != move->n_opp_to_link)
		      Update_PMat_At_Given_Edge(move->n_link->b[j],st->treelist->tree[gt]);
	      
		  /* Compute the likelihood */
		  Update_P_Lk(st->treelist->tree[gt],
			      move->b_opp_to_link,
			      move->n_link);
		  
		  move->lnL = Lk(move->b_opp_to_link,st->treelist->tree[gt]);

		  
		  st->tree->spr_list[i]->lnL += move->lnL;

		  /* Record branch lengths */
		  dir_v1 = dir_v2 = dir_v0 = -1;
		  For(j,3)
		    {
		      if(move->n_link->v[j] == move->n_opp_to_link) dir_v0 = j;
		      else if(dir_v1 < 0)                           dir_v1 = j;
		      else                                          dir_v2 = j;
		    }
		  
		  move->l0 = move->n_link->b[dir_v0]->l->v;
		  
		  if(move->n_link->v[dir_v1]->num > move->n_link->v[dir_v2]->num)
		    {
		      move->l1 = move->n_link->b[dir_v2]->l->v;
		      move->l2 = move->n_link->b[dir_v1]->l->v;
		    }
		  else
		    {
		      move->l1 = move->n_link->b[dir_v1]->l->v;
		      move->l2 = move->n_link->b[dir_v2]->l->v;
		    }
		  	  
		  /* Regraft the subtree at its original position */
		  Prune_Subtree(move->n_link,
				move->n_opp_to_link,
				&move->b_target,
				&b_residual,
				st->treelist->tree[gt]);
		  
		  Graft_Subtree(init_target,
				move->n_link,
				b_residual,
				st->treelist->tree[gt]);
		  
		  /* Restore branch lengths */
		  Restore_Br_Len(st->treelist->tree[gt]);
	      
		  /* Update relevant change proba matrices */
		  Update_PMat_At_Given_Edge(move->b_target,st->treelist->tree[gt]);
		  For(j,3) Update_PMat_At_Given_Edge(move->n_link->b[j],st->treelist->tree[gt]);
		  
		  /* Update relevant partial likelihoods */
		  For(j,3) Update_P_Lk(st->treelist->tree[gt],move->n_link->b[j],move->n_link);
		  
		  st->treelist->tree[gt]->c_lnL = init_lnL;
		}
	    }
	  else
	    {
	      st->tree->spr_list[i]->lnL += st->treelist->tree[gt]->c_lnL;
	    }
	}

      if(st->tree->spr_list[i]->lnL > best_lnL)
	{
	  best_lnL  = st->tree->spr_list[i]->lnL;
	  best_move = i;
	}
    }

  return best_move;  

}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


int PART_Try_One_Spr_Move(t_spr *st_move, supert_tree *st)
{
  int j;
  t_spr **gt_move;
  t_edge **init_target, **b_residual;
  int dir_v0, dir_v1, dir_v2;
  int gt;
  int gt_move_num;
  int n_moves;

  
  init_target = (t_edge **)mCalloc(st->n_part,sizeof(t_edge *));
  b_residual  = (t_edge **)mCalloc(st->n_part,sizeof(t_edge *));
  gt_move     = (t_spr **)mCalloc(st->n_part,sizeof(t_spr *));


  n_moves = 0;
  For(gt,st->n_part) 
    {
      gt_move_num = Map_Spr_Move(st_move->b_opp_to_link,
				 st_move->b_target,
				 st_move->n_link,
				 st->treelist->tree[gt],st);
      
      if(gt_move_num > -1)
	{
	  n_moves++;

	  gt_move[gt] = st->treelist->tree[gt]->spr_list[gt_move_num];
	  
	  if(gt_move[gt]->b_target)
	    {
	      /* Record t_edge lengths */
	      Record_Br_Len(st->treelist->tree[gt]);

	      /* Prune subtree */
	      Prune_Subtree(gt_move[gt]->n_link,
			    gt_move[gt]->n_opp_to_link,
			    &(init_target[gt]),
			    &(b_residual[gt]),
			    st->treelist->tree[gt]);
	      
	      /* Rough optimisation of the branch length */
	      Fast_Br_Len(init_target[gt],st->treelist->tree[gt],0);
	      
	      /* Update the change proba matrix at prune position */
	      Update_PMat_At_Given_Edge(init_target[gt],st->treelist->tree[gt]); /* TO DO : NECESSARY ?? */
	      
	      /* Update partial likelihood along the path from the prune to
		 the regraft position */
	      Update_P_Lk_Along_A_Path(gt_move[gt]->path,gt_move[gt]->depth_path,st->treelist->tree[gt]); /* TO DO : NECESSARY ?? */
	      
	      /* Regraft subtree */
	      Graft_Subtree(gt_move[gt]->b_target,gt_move[gt]->n_link,b_residual[gt],st->treelist->tree[gt]);
	      
	      dir_v1 = dir_v2 = dir_v0 = -1;
	      For(j,3)
		{
		  if(gt_move[gt]->n_link->v[j] == gt_move[gt]->n_opp_to_link) dir_v0 = j;
		  else if(dir_v1 < 0)                                         dir_v1 = j;
		  else                                                        dir_v2 = j;
		}
	      
	      gt_move[gt]->n_link->b[dir_v0]->l->v = gt_move[gt]->l0;
		  
	      if(gt_move[gt]->n_link->v[dir_v1]->num > gt_move[gt]->n_link->v[dir_v2]->num)
		{
		  gt_move[gt]->n_link->b[dir_v2]->l->v = gt_move[gt]->l1;
		  gt_move[gt]->n_link->b[dir_v1]->l->v = gt_move[gt]->l2;
		}
	      else
		{
		  gt_move[gt]->n_link->b[dir_v1]->l->v = gt_move[gt]->l1;
		  gt_move[gt]->n_link->b[dir_v2]->l->v = gt_move[gt]->l2;
		}
	    }
	}
    }
  

  if(n_moves)
    {
      if(st_move->lnL > st->tree->best_lnL)
	{
	  t_edge *st_target, *st_residual;

	  /* Apply the move on the super-tree */
	  Prune_Subtree(st_move->n_link,
			st_move->n_opp_to_link,
			&st_target,
			&st_residual,
			st->tree);
	  
	  Graft_Subtree(st_move->b_target,
			st_move->n_link,
			st_residual,
			st->tree);
	  
	  	  
	  /* Map gt and st nodes and edges */
	  PART_Do_Mapping(st);

	  time(&(st->tree->t_current));

	  Set_Both_Sides(YES,st->tree);	  
	  st->tree->c_lnL      = PART_Lk(st);
	  st->tree->c_pars     = PART_Pars(st);
	  
	  
	  if(FABS(st->tree->c_lnL - st_move->lnL) > st->tree->mod->s_opt->min_diff_lk_local)
	    {
	      PhyML_Printf("\n. st->tree->c_lnL = %f st_move->lnL = %f\n",
		     st->tree->c_lnL,st_move->lnL);

	      For(gt,st->n_part)
		{
		  PhyML_Printf("\n. truth -> %f ; move -> %f",
			       Lk(NULL,st->treelist->tree[gt]),
			       gt_move[gt] ? gt_move[gt]->lnL : -1.);
		}
	    }
	  
	  PhyML_Printf("\n. (%5d sec) [+ ] [%10.2f] [%5d] -- ",
		 (int)(st->tree->t_current - st->tree->t_beg),
		 st->tree->c_lnL,st->tree->c_pars);	  
	  
	  For(gt,st->n_part)
	    PhyML_Printf("[%10.2f] ",st->treelist->tree[gt]->c_lnL);
	  
	  
	  st->tree->n_improvements++;
	  st->tree->best_lnL = st->tree->c_lnL;
	  For(gt,st->n_part) Record_Br_Len(st->treelist->tree[gt]);
	  
	  Free(init_target);
	  Free(b_residual);
	  Free(gt_move);
	  
	  return 1;
	}
/*       else */
/* 	{ */
/* 	  For(gt,st->n_part)  */
/* 	    { */
/* 	      if(gt_move[gt]) */
/* 		{ */
/* 		  Lk(st->treelist->tree[gt]); */
/* 		  Fast_Br_Len_Recur(st->treelist->tree[gt]->a_nodes[0], */
/* 				    st->treelist->tree[gt]->a_nodes[0]->v[0], */
/* 				    st->treelist->tree[gt]->a_nodes[0]->b[0], */
/* 				    st->treelist->tree[gt]); */
/* 		} */
/* 	    } */
	  
/* 	  time(&(st->tree->t_current)); */
/* 	  st->tree->both_sides = 1; */
/* 	  st->tree->c_lnL      = PART_Lk(st); */
	  
/* 	  if(st->tree->c_lnL > st->tree->best_lnL) */
/* 	    { */
/* 	      t_edge *st_target, *st_residual; */
	      
/* 	      /\* Apply the move on the super-tree *\/ */
/* 	      Prune_Subtree(st_move->n_link, */
/* 			    st_move->n_opp_to_link,			     */
/* 			    &st_target, */
/* 			    &st_residual, */
/* 			    st->tree); */
	      
/* 	      Graft_Subtree(st_move->b_target, */
/* 			    st_move->n_link, */
/* 			    st_residual, */
/* 			    st->tree); */
	      
	      
/* 	      /\* Map gt and st nodes and edges *\/ */
/* 	      PART_Do_Mapping(st); */


/* 	      st->tree->c_pars = PART_Pars(st); */
/* 	      PhyML_Printf("\n. (%5d sec) [++] [%10.2f] [%5d] -- ", */
/* 		     (int)(st->tree->t_current-st->tree->t_beg), */
/* 		     st->tree->c_lnL, */
/* 		     st->tree->c_pars); */
/* 	      For(gt,st->n_part) */
/* 		PhyML_Printf("[%10.2f] ",st->treelist->tree[gt]->c_lnL); */
	      
/* 	      st->tree->n_improvements++; */
/* 	      st->tree->best_lnL = st->tree->c_lnL; */
/* 	      For(gt,st->n_part) Record_Br_Len(st->treelist->tree[gt]); */

/* 	      Free(init_target); */
/* 	      Free(b_residual); */
/* 	      Free(gt_move); */

/* 	      return 1; */
/* 	    } */
/* 	} */
    }
  
  For(gt,st->n_part) 
    {
      if(gt_move[gt])
	{	  
	  /* Regraft the subtree at its original position */
	  Prune_Subtree(gt_move[gt]->n_link,
			gt_move[gt]->n_opp_to_link,
			&(gt_move[gt]->b_target),
			&(b_residual[gt]),
			st->treelist->tree[gt]);

	  Graft_Subtree(init_target[gt],
			gt_move[gt]->n_link,
			b_residual[gt],
			st->treelist->tree[gt]);	  

	  /* Restore branch lengths */
	  Restore_Br_Len(st->treelist->tree[gt]);
	}
    }
  
  Set_Both_Sides(YES,st->tree);
  st->tree->c_lnL      = PART_Lk(st);
  st->tree->c_pars     = PART_Pars(st);

  time(&(st->tree->t_current));
  
  PhyML_Printf("\n. (%5d sec) [--] [%10.2f] [%5d] -- ",
	 (int)(st->tree->t_current - st->tree->t_beg),
	 st->tree->c_lnL,st->tree->c_pars);	  
  
  For(gt,st->n_part) PhyML_Printf("[%10.2f] ",st->treelist->tree[gt]->c_lnL);

  Free(init_target);
  Free(b_residual);
  Free(gt_move);

  return 0;

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void PART_NNI(t_edge *st_b, supert_tree *st)
{  
  t_node *v1, *v2, *v3, *v4;
  phydbl lk0_opt, lk1_opt, lk2_opt;
  int i,j;
  phydbl *init_bl;
  t_edge **map_edge_bef_swap, **map_edge_aft_swap;


  init_bl = (phydbl *)mCalloc(st->n_bl_part,sizeof(phydbl));
  map_edge_bef_swap = (t_edge **)mCalloc(st->n_part,sizeof(t_edge *));
  map_edge_aft_swap = (t_edge **)mCalloc(st->n_part,sizeof(t_edge *));


  v1 = st_b->left->v[st_b->l_v1];
  v2 = st_b->left->v[st_b->l_v2];
  v3 = st_b->rght->v[st_b->r_v1];
  v4 = st_b->rght->v[st_b->r_v2];

  lk0_opt  = lk1_opt  = lk2_opt  = UNLIKELY;

  if(v1->num < v2->num)
    {
      Check_Dirs(st->tree);
      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  if(v3->num < v4->num)
    {
      Check_Dirs(st->tree);
      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  
/*   PhyML_Printf("oooooooo\n"); */
/*   Print_Node(st->tree->a_nodes[0], */
/* 	     st->tree->a_nodes[0]->v[0], */
/* 	     st->tree); */
/*   PhyML_Printf(">>>>>>>\n"); */
/*   For(i,st->n_part) */
/*     { */
/*       Print_Node(st->treelist->tree[i]->a_nodes[0], */
/* 		 st->treelist->tree[i]->a_nodes[0]->v[0], */
/* 		 st->treelist->tree[i]); */
/*       PhyML_Printf("<<<<<<<\n"); */
/*     } */

  
  PART_Record_Br_Len(st);

  For(i,st->n_part) map_edge_bef_swap[i] = NULL;
  For(i,st->n_part) if(st->map_st_edge_in_gt[i][st_b->num]) map_edge_bef_swap[i] = st->map_st_edge_in_gt[i][st_b->num];

  /* First alternative topological configuration */
  /* Swap */
  PART_Swap(v2,st_b->left,st_b->rght,v3,st);
  Swap(v2,st_b->left,st_b->rght,v3,st->tree);
  PART_Do_Mapping(st);
  PART_Set_Bl(st->bl,st);
  For(i,st->n_part) map_edge_aft_swap[i] = NULL;
  For(i,st->n_part) if(st->map_st_edge_in_gt[i][st_b->num]) map_edge_aft_swap[i] = st->map_st_edge_in_gt[i][st_b->num];
  For(i,st->n_part) if(map_edge_bef_swap[i]) Update_PMat_At_Given_Edge(map_edge_bef_swap[i],st->treelist->tree[i]);
  For(i,st->n_part) if(map_edge_aft_swap[i]) Update_PMat_At_Given_Edge(map_edge_aft_swap[i],st->treelist->tree[i]);
  For(i,st->n_part) if(map_edge_bef_swap[i] && map_edge_aft_swap[i]) 
    {
      For(j,3) if((map_edge_aft_swap[i]->left->v[j]) && 
		  (map_edge_aft_swap[i]->left->b[j] == map_edge_bef_swap[i]) &&
		  (map_edge_aft_swap[i] != map_edge_bef_swap[i])) 
	Update_P_Lk(st->treelist->tree[i],
		    map_edge_aft_swap[i],
		    map_edge_aft_swap[i]->left);
      For(j,3) if((map_edge_aft_swap[i]->rght->v[j]) && 
		  (map_edge_aft_swap[i]->rght->b[j] == map_edge_bef_swap[i]) &&
		  (map_edge_aft_swap[i] != map_edge_bef_swap[i])) 
	Update_P_Lk(st->treelist->tree[i],
		    map_edge_aft_swap[i],
		    map_edge_aft_swap[i]->rght);
    }
  PART_Update_Lk_At_Given_Edge(st_b,st);
  lk1_opt  = PART_Br_Len_Brent(st_b,0,st);
  For(i,st->n_part) st->bl1[st->bl_partition[i]][st_b->num] = st->bl[st->bl_partition[i]][st_b->num];
  /* Unswap */
  PART_Swap(v3,st_b->left,st_b->rght,v2,st);
  Swap(v3,st_b->left,st_b->rght,v2,st->tree);
  PART_Do_Mapping(st);
  PART_Restore_Br_Len(st);
  PART_Set_Bl(st->bl,st);
  For(i,st->n_part) if(map_edge_bef_swap[i]) Update_PMat_At_Given_Edge(map_edge_bef_swap[i],st->treelist->tree[i]);
  For(i,st->n_part) if(map_edge_aft_swap[i]) Update_PMat_At_Given_Edge(map_edge_aft_swap[i],st->treelist->tree[i]);
  For(i,st->n_part) if(map_edge_bef_swap[i] && map_edge_aft_swap[i]) 
    {
      For(j,3) if((map_edge_aft_swap[i]->left->v[j]) && 
		  (map_edge_aft_swap[i]->left->b[j] == map_edge_bef_swap[i]) &&
		  (map_edge_aft_swap[i] != map_edge_bef_swap[i])) 
	Update_P_Lk(st->treelist->tree[i],
		    map_edge_aft_swap[i],
		    map_edge_aft_swap[i]->left);
      For(j,3) if((map_edge_aft_swap[i]->rght->v[j]) && 
		  (map_edge_aft_swap[i]->rght->b[j] == map_edge_bef_swap[i]) &&
		  (map_edge_aft_swap[i] != map_edge_bef_swap[i])) 
	Update_P_Lk(st->treelist->tree[i],
		    map_edge_aft_swap[i],
		    map_edge_aft_swap[i]->rght);
    }




  /* Second alternative topological configuration */
  /* Swap */
  PART_Swap(v2,st_b->left,st_b->rght,v4,st);
  Swap(v2,st_b->left,st_b->rght,v4,st->tree);
  PART_Do_Mapping(st);
  PART_Set_Bl(st->bl,st);
  For(i,st->n_part) map_edge_aft_swap[i] = NULL;
  For(i,st->n_part) if(st->map_st_edge_in_gt[i][st_b->num]) map_edge_aft_swap[i] = st->map_st_edge_in_gt[i][st_b->num];
  For(i,st->n_part) if(map_edge_bef_swap[i]) Update_PMat_At_Given_Edge(map_edge_bef_swap[i],st->treelist->tree[i]);
  For(i,st->n_part) if(map_edge_aft_swap[i]) Update_PMat_At_Given_Edge(map_edge_aft_swap[i],st->treelist->tree[i]);
  For(i,st->n_part) if(map_edge_bef_swap[i] && map_edge_aft_swap[i]) 
    {
      For(j,3) if((map_edge_aft_swap[i]->left->v[j]) && 
		  (map_edge_aft_swap[i]->left->b[j] == map_edge_bef_swap[i]) &&
		  (map_edge_aft_swap[i] != map_edge_bef_swap[i])) 
	Update_P_Lk(st->treelist->tree[i],
		    map_edge_aft_swap[i],
		    map_edge_aft_swap[i]->left);
      For(j,3) if((map_edge_aft_swap[i]->rght->v[j]) && 
		  (map_edge_aft_swap[i]->rght->b[j] == map_edge_bef_swap[i]) &&
		  (map_edge_aft_swap[i] != map_edge_bef_swap[i])) 
	Update_P_Lk(st->treelist->tree[i],
		    map_edge_aft_swap[i],
		    map_edge_aft_swap[i]->rght);
    }

  PART_Update_Lk_At_Given_Edge(st_b,st);
  lk2_opt  = PART_Br_Len_Brent(st_b,0,st);
  For(i,st->n_part) st->bl2[st->bl_partition[i]][st_b->num] = st->bl[st->bl_partition[i]][st_b->num];
  /*   PhyML_Printf("\n. lk2_init = %f lk2_opt = %f",lk2_init,lk2_opt); */
  /* Unswap */
  PART_Swap(v4,st_b->left,st_b->rght,v2,st);
  Swap(v4,st_b->left,st_b->rght,v2,st->tree);
  PART_Do_Mapping(st);
  PART_Restore_Br_Len(st);
  PART_Set_Bl(st->bl,st);
  For(i,st->n_part) if(map_edge_bef_swap[i]) Update_PMat_At_Given_Edge(map_edge_bef_swap[i],st->treelist->tree[i]);
  For(i,st->n_part) if(map_edge_aft_swap[i]) Update_PMat_At_Given_Edge(map_edge_aft_swap[i],st->treelist->tree[i]);
  For(i,st->n_part) if(map_edge_bef_swap[i] && map_edge_aft_swap[i]) 
    {
      For(j,3) if((map_edge_aft_swap[i]->left->v[j]) && 
		  (map_edge_aft_swap[i]->left->b[j] == map_edge_bef_swap[i]) &&
		  (map_edge_aft_swap[i] != map_edge_bef_swap[i])) 
	Update_P_Lk(st->treelist->tree[i],
		    map_edge_aft_swap[i],
		    map_edge_aft_swap[i]->left);
      For(j,3) if((map_edge_aft_swap[i]->rght->v[j]) && 
		  (map_edge_aft_swap[i]->rght->b[j] == map_edge_bef_swap[i]) &&
		  (map_edge_aft_swap[i] != map_edge_bef_swap[i])) 
	Update_P_Lk(st->treelist->tree[i],
		    map_edge_aft_swap[i],
		    map_edge_aft_swap[i]->rght);
    }


  /* Back to the initial topological configuration 
   * and branch lengths.
   */
  PART_Do_Mapping(st);
  PART_Set_Bl(st->bl,st);
  PART_Restore_Br_Len(st);
  For(i,st->n_part) map_edge_aft_swap[i] = NULL;
  For(i,st->n_part) if(st->map_st_edge_in_gt[i][st_b->num]) map_edge_aft_swap[i] = st->map_st_edge_in_gt[i][st_b->num];
  For(i,st->n_part) if(map_edge_bef_swap[i]) Update_PMat_At_Given_Edge(map_edge_bef_swap[i],st->treelist->tree[i]);
  For(i,st->n_part) if(map_edge_aft_swap[i]) Update_PMat_At_Given_Edge(map_edge_aft_swap[i],st->treelist->tree[i]);
  PART_Update_Lk_At_Given_Edge(st_b,st);
  lk0_opt  = PART_Br_Len_Brent(st_b,0,st);
  For(i,st->n_part) st->bl0[st->bl_partition[i]][st_b->num] = st->bl[st->bl_partition[i]][st_b->num];

  PART_Restore_Br_Len(st);
  PART_Set_Bl(st->bl,st);
  For(i,st->n_part) if(map_edge_bef_swap[i]) Update_PMat_At_Given_Edge(map_edge_bef_swap[i],st->treelist->tree[i]);
  For(i,st->n_part) if(map_edge_aft_swap[i]) Update_PMat_At_Given_Edge(map_edge_aft_swap[i],st->treelist->tree[i]);
  PART_Update_Lk_At_Given_Edge(st_b,st);



/*   For(i,2*st->tree->n_otu-3) */
/*     PhyML_Printf("\n. 3 Edge %3d --> lnL=%f",i,PART_Lk_At_Given_Edge(st->tree->a_edges[i],st)); */



  st_b->nni->lk0 = lk0_opt;
  st_b->nni->lk1 = lk1_opt;
  st_b->nni->lk2 = lk2_opt;

  st_b->nni->score = lk0_opt - MAX(lk1_opt,lk2_opt);

  if((st_b->nni->score <  st->tree->mod->s_opt->min_diff_lk_local) &&
     (st_b->nni->score > -st->tree->mod->s_opt->min_diff_lk_local))
    {
      st_b->nni->score = .0;
      st_b->nni->lk1   = st_b->nni->lk0;
      st_b->nni->lk2   = st_b->nni->lk0;
     }

  PART_Restore_Br_Len(st);
  PART_Update_Lk_At_Given_Edge(st_b,st); /* to replace by PART_Update_PMat_At_Given_Edge(st_b,st); */
/*   PhyML_Printf("\n. lk_end = %f",st->tree->c_lnL); */
/*   For(i,2*st->tree->n_otu-3) PhyML_Printf("\n. %f",PART_Lk_At_Given_Edge(st->tree->a_edges[i],st)); */
/*   PhyML_Printf("\n. lk_end = %f",PART_Lk(st)); */
/*   PhyML_Printf("\n"); */

/*   PhyML_Printf("\n. Edge %3d, score = %20f",st_b->num,st_b->nni->score); */

  if(st_b->num == 90)
    PhyML_Printf("\n. v1=%d v2=%d v3=%d v4=%d left-%d right-%d",
	   v1->num,
	   v2->num,
	   v3->num,
	   v4->num,
	   st_b->left->num,
	   st_b->rght->num);



  if(lk0_opt > MAX(lk1_opt,lk2_opt))
    {
      st_b->nni->best_conf = 0;
      st_b->nni->swap_node_v1 = NULL;
      st_b->nni->swap_node_v2 = NULL;
      st_b->nni->swap_node_v3 = NULL;
      st_b->nni->swap_node_v4 = NULL;
    }
  else if(lk1_opt > MAX(lk0_opt,lk2_opt))
    {
      st_b->nni->best_conf    = 1;
      st_b->nni->swap_node_v1 = v2;
      st_b->nni->swap_node_v2 = st_b->left;
      st_b->nni->swap_node_v3 = st_b->rght;
      st_b->nni->swap_node_v4 = v3;
    }
  else if(lk2_opt > MAX(lk0_opt,lk1_opt))
    {
      st_b->nni->best_conf    = 2;
      st_b->nni->swap_node_v1 = v2;
      st_b->nni->swap_node_v2 = st_b->left;
      st_b->nni->swap_node_v3 = st_b->rght;
      st_b->nni->swap_node_v4 = v4;
    }
  else
    {
      st_b->nni->score        = +1.0;
      st_b->nni->best_conf    = 0;
      st_b->nni->swap_node_v1 = NULL;
      st_b->nni->swap_node_v2 = NULL;
      st_b->nni->swap_node_v3 = NULL;
      st_b->nni->swap_node_v4 = NULL;
    }
  
  Free(init_bl);
  Free(map_edge_aft_swap);
  Free(map_edge_bef_swap);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void PART_Swap(t_node *st_a, t_node *st_b, t_node *st_c, t_node *st_d, supert_tree *st)
{
  int i,j;
  t_node *gt_a, *gt_b, *gt_c, *gt_d;
  int ab, ba, cd, dc, bc;

  ab = ba = cd = dc = bc = -1;

  For(i,3) if(st_a->v[i] == st_b) { ab = i; break; }
  For(i,3) if(st_b->v[i] == st_a) { ba = i; break; }
  For(i,3) if(st_c->v[i] == st_d) { cd = i; break; }
  For(i,3) if(st_d->v[i] == st_c) { dc = i; break; }
  For(i,3) if(st_b->v[i] == st_c) { bc = i; break; }

  if(ab < 0 || ba < 0 || cd < 0 || dc < 0)
    {
      PhyML_Printf("\n. Nodes %d %d %d %d\n",st_a->num,st_b->num,st_c->num,st_d->num);
      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  gt_a = gt_b = gt_c = gt_d = NULL;
  
  For(i,st->n_part)
    {
      gt_b = st->match_st_node_in_gt[i][st_b->num];
      gt_c = st->match_st_node_in_gt[i][st_c->num];
      
      if(gt_b && gt_c) /* The st t_edge with st_b and st_c at its extremities
			* matches an t_edge in gt 
		        */
	{
#ifdef DEBUG
	  For(j,3) if((gt_b->v[j]) && (gt_b->v[j] == gt_c)) break;
	  if(j == 3)
	    {
	      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
	      Warn_And_Exit("");
	    }
#endif
	  gt_a = st->map_st_node_in_gt[i][st_a->num][ab][0];
	  gt_d = st->map_st_node_in_gt[i][st_d->num][dc][0];
	  Swap(gt_a,gt_b,gt_c,gt_d,st->treelist->tree[i]);
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void PART_Set_Bl(phydbl **bl, supert_tree *st)
{
  int i,j;
  t_edge *gt_b;

  gt_b = NULL;
						
  /* Set all the actual branch lengths to 0.0 
   */
  For(i,st->n_part)
    {
      For(j,2*st->treelist->tree[i]->n_otu-3)
	{
	  gt_b = st->treelist->tree[i]->a_edges[j];
	  gt_b->l->v = .0;
	}
    }

  /* Update every branch length 
   */  
  For(i,2*st->tree->n_otu-3)
    {
      For(j,st->n_part)
	{
	  gt_b = st->map_st_edge_in_gt[j][i];	
	  
	  /* Need to make sure that st->tree->a_edges[i] is on an existing path in gt */
	  if((st->map_st_node_in_gt[j][st->tree->a_edges[i]->left->num][st->tree->a_edges[i]->l_r][0]) &&
	     (st->map_st_node_in_gt[j][st->tree->a_edges[i]->rght->num][st->tree->a_edges[i]->r_l][0]))
	    {
	      gt_b->l->v += bl[st->bl_partition[j]][i];
	    }
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void PART_Record_Br_Len(supert_tree *st)
{
  int i,j;
  For(i,st->n_part) For(j,2*st->tree->n_otu-3) st->bl_cpy[i][j] = st->bl[i][j];
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void PART_Restore_Br_Len(supert_tree *st)
{
  int i,j;
  For(i,st->n_part) For(j,2*st->tree->n_otu-3) st->bl[i][j] = st->bl_cpy[i][j];
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl PART_Lk(supert_tree *st)
{
  int i;

  PART_Do_Mapping(st);
  PART_Set_Bl(st->bl,st);  

  st->tree->c_lnL = .0;
  For(i,st->n_part) 
    {
      Set_Both_Sides(YES,st->treelist->tree[i]);	  
      Lk(NULL,st->treelist->tree[i]);
/*       PhyML_Printf("\n. Tree %3d lnL = %f",i+1,st->treelist->tree[i]->c_lnL); */
      st->tree->c_lnL += st->treelist->tree[i]->c_lnL;
    }
  return st->tree->c_lnL;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl PART_Lk_At_Given_Edge(t_edge *st_b, supert_tree *st)
{
  int i;
  t_edge *gt_b;
  phydbl lnL;

  PART_Set_Bl(st->bl,st);

  gt_b = NULL;
  st->tree->c_lnL = .0;
  lnL = .0;
  For(i,st->n_part)
    {      
      gt_b = st->map_st_edge_in_gt[i][st_b->num];
      lnL = Lk(gt_b,st->treelist->tree[i]);
      st->tree->c_lnL += lnL;
/*       PhyML_Printf("\n. gt %d st t_edge %d gt t_edge %d lnL=%f l=%f ",i,st_b->num,gt_b->num,lnL,gt_b->l->v); */
    }
  return st->tree->c_lnL;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl PART_Update_Lk_At_Given_Edge(t_edge *st_b, supert_tree *st)
{
  int i;
  t_edge *gt_b;

  PART_Set_Bl(st->bl,st);
  
  gt_b = NULL;
  st->tree->c_lnL = .0;
  For(i,st->n_part)
    {
      gt_b = st->map_st_edge_in_gt[i][st_b->num];
      if(gt_b) st->tree->c_lnL += Update_Lk_At_Given_Edge(gt_b,st->treelist->tree[i]);
      else     st->tree->c_lnL += st->treelist->tree[i]->c_lnL;
    }
  return st->tree->c_lnL;
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void PART_Fill_Model_Partitions_Table(supert_tree *st)
{
  int i,j;
  char *c;
  char *abc;
  int lig, col;
  int n_groups;
  int *encountered_vals;

  c = (char *)mCalloc(10,sizeof(char));
  abc = (char *)mCalloc(20,sizeof(char));
  encountered_vals = (int *)mCalloc(st->n_part,sizeof(int));

  strcpy(abc,"ABCDEFGHIJKLMNOP\0");

  PhyML_Printf("\n\n\n");
  lig = col = 0;
  while(1)
    {
      PhyML_Printf("\n\n");
      For(i,st->n_part)
	PhyML_Printf(". Data set %3d : %s\n",i+1,st->optionlist[i]->in_align_file);

      PhyML_Printf("\n. Data set             ");
      For(i,st->n_part) PhyML_Printf("%3d ",i+1);
      PhyML_Printf("\n. -A- t_edge lengths     ");
      For(i,st->n_part) PhyML_Printf("%3d ",st->bl_partition[i]);

      if(lig == 1) break;

      PhyML_Printf("\n. (%c-%2d)> ",abc[lig],col+1);
      Getstring_Stdin(c);
      
      switch(lig)
	{
	case 0 :
	  {
	    st->bl_partition[col] = atoi(c);
	    break;
	  }
	default :
	  {
	    break;
	  }
	}

      col++;

      if(col == st->n_part)
	{
	  col = 0;
	  lig++;
	}
    }
    
  n_groups = 0;
  For(i,st->n_part) 
    {
      For(j,n_groups)
	if(encountered_vals[j] == st->bl_partition[i])
	  break;

      if(j == n_groups) 
	{
	  encountered_vals[n_groups] = st->bl_partition[i];
	  n_groups++;
	}
      st->bl_partition[i] = j;
    }
  
  st->n_bl_part = n_groups;

  Free(encountered_vals);
  Free(c);
  Free(abc);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl PART_Br_Len_Brent(t_edge *st_b, int quickdirty, supert_tree *st)
{
  phydbl ax,  cx;
  int part;
  phydbl cur_l;


  For(part,st->n_bl_part)
    {
      cur_l = st->bl[part][st_b->num];
      
      ax = 10.*cur_l;
      cx = st->tree->mod->l_min;

      Generic_Brent_Lk(&(st->bl[part][st_b->num]),
		       ax,cx,
		       st->tree->mod->s_opt->min_diff_lk_local,
		       st->tree->mod->s_opt->brent_it_max,
		       st->tree->mod->s_opt->quickdirty,
		       Wrap_Part_Lk_At_Given_Edge,st_b,NULL,st,NO);
      
    }
  return st->tree->c_lnL;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void PART_Initialise_Bl_Partition(supert_tree *st)
{
  int i,j;
  
  For(i,st->n_bl_part)
    {
      For(j,2*st->tree->n_otu-3)
	{
	  st->bl[i][j] = .1;
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void PART_Optimize_Br_Len_Serie(t_node *st_a, t_node *st_d, t_edge *st_b, supert_tree *st)
{
  phydbl lk_init;
  int i;

  lk_init = st->tree->c_lnL;
  
  PART_Br_Len_Brent(st_b,0,st);

  if(st->tree->c_lnL < lk_init - st->tree->mod->s_opt->min_diff_lk_local)
    { 
      PhyML_Printf("\n== %f -- %f",lk_init,st->tree->c_lnL);
      PhyML_Printf("\n== Err. in file %s at line %d\n\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
    
  if(st_d->tax) return;
  else For(i,3) if(st_d->v[i] != st_a)
    {
      PART_Update_P_Lk(st_d->b[i],st_d,st);
      PART_Optimize_Br_Len_Serie(st_d,st_d->v[i],st_d->b[i],st);
    }

  For(i,3) if((st_d->v[i] == st_a) && (!st_d->v[i]->tax)) PART_Update_P_Lk(st_d->b[i],st_d,st);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void PART_Update_P_Lk(t_edge *st_b, t_node *st_n, supert_tree *st)
{
  int i,dir;

  dir = -1;
  For(i,3) if((st_n->b[i]) && (st_n->b[i] == st_b)) {dir = i; break;}
  
  if(dir < 0)
    {
      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  For(i,st->n_part)
    {
      if((st->map_st_node_in_gt[i][st_n->num][dir][0]) && (!st->map_st_node_in_gt[i][st_n->num][dir][0]->tax))
	{	  
	  Update_P_Lk(st->treelist->tree[i],
		      st->map_st_edge_in_gt[i][st_b->num],
		      st->map_st_node_in_gt[i][st_n->num][dir][0]);
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void PART_Make_N_Swap(t_edge **st_b, int beg, int end, supert_tree *st)
{
  int i;
  int dim;

  dim = 2*st->tree->n_otu-2;

  st->tree->n_swap = 0;
  for(i=beg;i<end;i++)
    {
      if(st_b[i]->left->tax || st_b[i]->rght->tax)
	{
	  PhyML_Printf("\n. Edge %d is external.",st_b[i]->num);
	  PhyML_Printf("\n. Err in file %s at line %d\n\n.",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}
      
      PART_Swap(st_b[i]->nni->swap_node_v2->v[st->tree->t_dir[st_b[i]->nni->swap_node_v2->num*dim+st_b[i]->nni->swap_node_v1->num]],
	      st_b[i]->nni->swap_node_v2,
	      st_b[i]->nni->swap_node_v3,
	      st_b[i]->nni->swap_node_v3->v[st->tree->t_dir[st_b[i]->nni->swap_node_v3->num*dim+st_b[i]->nni->swap_node_v4->num]],
	      st);

      Swap(st_b[i]->nni->swap_node_v2->v[st->tree->t_dir[st_b[i]->nni->swap_node_v2->num*dim+st_b[i]->nni->swap_node_v1->num]],
	   st_b[i]->nni->swap_node_v2,
	   st_b[i]->nni->swap_node_v3,
	   st_b[i]->nni->swap_node_v3->v[st->tree->t_dir[st_b[i]->nni->swap_node_v3->num*dim+st_b[i]->nni->swap_node_v4->num]],
	   st->tree);

      PART_Do_Mapping(st);

      st->tree->n_swap++;
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void PART_Update_Bl(phydbl fact, supert_tree *st)
{
  int i,j;
  
  For(i,2*st->tree->n_otu-3)
    {
      For(j,st->n_part)
	st->bl[st->bl_partition[j]][i] = 
	st->bl_cpy[st->bl_partition[j]][i] + 
	(st->bl0[st->bl_partition[j]][i] - st->bl_cpy[st->bl_partition[j]][i]) * fact;
    }
  PART_Set_Bl(st->bl,st);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void PART_Update_Bl_Swaped(t_edge **st_b, int n, supert_tree *st)
{
  int i,j;
  
  For(i,n)
    {
      For(j,st->n_part)
	{
	  st->bl[st->bl_partition[j]][st_b[i]->num] = 
	    (st_b[i]->nni->best_conf == 1)?
	    (st->bl1[st->bl_partition[j]][st_b[i]->num]):
	    (st->bl2[st->bl_partition[j]][st_b[i]->num]);
	}
    }
  PART_Set_Bl(st->bl,st);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void PART_Do_Mapping(supert_tree *st)
{
  int k;

  Fill_Dir_Table(st->tree);
  For(k,st->n_part)
    {
      Fill_Dir_Table(st->treelist->tree[k]);
      PART_Match_St_Nodes_In_Gt(st->treelist->tree[k],st);	      
      PART_Match_St_Edges_In_Gt(st->treelist->tree[k],st);	      
      PART_Map_St_Nodes_In_Gt(st->treelist->tree[k],st);
      PART_Map_St_Edges_In_Gt(st->treelist->tree[k],st);
      PART_Map_Gt_Edges_In_St(st->treelist->tree[k],st);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void PART_Print_Bl(supert_tree *st)
{
  int i,j;
  
  For(j,2*st->tree->n_otu-3)
    { 
      PhyML_Printf("\n. t_edge %4d ",j);
      For(i,st->n_bl_part)
	{
	  PhyML_Printf("%f ",st->bl[i][j]);
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

