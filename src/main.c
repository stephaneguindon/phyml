/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

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
#include "mixt.h"
#include "invitee.h"

#ifdef MPI
#include "mpi_boot.h"
#endif

#ifdef BEAGLE
#include "beagle_utils.h"
#endif



#if (defined PHYML || EVOLVE)

int main(int argc, char **argv)
{
  calign *cdata;
  option *io;
  t_tree *tree;
  int num_data_set;
  int num_tree,num_rand_tree;
  t_mod *mod;
  time_t t_beg,t_end;
  phydbl best_lnL;
  int r_seed;
  char *most_likely_tree=NULL;
  int orig_random_input_tree;

  
#ifdef MPI
  int rc;
  rc = MPI_Init(&argc,&argv);
  if (rc != MPI_SUCCESS) 
    {
      PhyML_Fprintf(stderr,"\n. Err. starting MPI program. Terminating.\n");
      MPI_Abort(MPI_COMM_WORLD, rc);
    }
  if(MPI_Comm_size(MPI_COMM_WORLD,&Global_numTask) != MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD, rc);
  if(MPI_Comm_rank(MPI_COMM_WORLD,&Global_myRank) != MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD, rc);
  PhyML_Fprintf(stdout,"\n\n. Running the analysis on %d CPU%s.",
                Global_numTask,
                Global_numTask>1?"s.":".");
#endif

#ifdef QUIET
  setvbuf(stdout,NULL,_IOFBF,2048);
#endif

  tree     = NULL;
  mod      = NULL;
  best_lnL = UNLIKELY;


  io = (option *)Get_Input(argc,argv);
  if(!io) return(0);
  else if(io->use_xml == YES)
    {
      Free(io);
      return(0);
    }
  
#ifdef EVOLVE
  io->colalias = NO;
#endif

  r_seed = (io->r_seed < 0)?(time(NULL)):(io->r_seed);
  srand(r_seed);
  io->r_seed = r_seed;

  if(io->in_tree == 2) Test_Multiple_Data_Set_Format(io);
  else io->n_trees = 1;

  if(io->n_trees == 0 && io->in_tree == 2)
    {
      PhyML_Printf("\n== Err.: the input tree file does not provide a tree in valid format.");
      Exit("\n");
    }

  if((io->n_data_sets > 1) && (io->n_trees > 1))
    {
      io->n_data_sets = MIN(io->n_trees,io->n_data_sets);
      io->n_trees     = MIN(io->n_trees,io->n_data_sets);
    }

  for(num_data_set=0;num_data_set<io->n_data_sets;num_data_set++)
    {
      best_lnL = UNLIKELY;
      Get_Seq(io);
      Make_Model_Complete(io->mod);
      Set_Model_Name(io->mod);
      Print_Settings(io);
      mod = io->mod;
      orig_random_input_tree = io->mod->s_opt->random_input_tree;
      

      if(io->data)
        {
          if(io->n_data_sets > 1) PhyML_Printf("\n. Data set [#%d]\n",num_data_set+1);
          cdata = Compact_Data(io->data,io);
          
          Free_Seq(io->data,cdata->n_otu);

          for(num_tree=(io->n_trees == 1)?(0):(num_data_set);num_tree < io->n_trees;num_tree++)
            {
              if(io->mod->s_opt->random_input_tree == NO) io->mod->s_opt->n_rand_starts = 1;

              if(orig_random_input_tree == YES && io->n_trees > 1)
                {
                  PhyML_Printf("\n== Cannot combine random starting trees with multiple input trees.");
                  Exit("\n");
                }

              for(num_rand_tree=0;num_rand_tree<io->mod->s_opt->n_rand_starts;num_rand_tree++)
                {
                  if((io->mod->s_opt->random_input_tree) && (io->mod->s_opt->topo_search != NNI_MOVE))
                    if(!io->quiet) PhyML_Printf("\n\n. [Random start %3d/%3d]",num_rand_tree+1,io->mod->s_opt->n_rand_starts);

                  Init_Model(cdata,mod,io);
                  Set_Model_Parameters(mod);
                  
#ifdef M4
                  if(io->mod->use_m4mod) M4_Init_Model(mod->m4mod,cdata,mod);
#endif
                  
                  // Make the initial tree
                  switch(io->in_tree)
                    {
                    case 0 : case 1 :
                      {
                        tree = Dist_And_BioNJ(cdata,mod,io);
                        if(io->print_mat_and_exit == YES)
                          {
                            Print_Mat(ML_Dist(cdata,mod));
                            exit(-1);
                          }
                        break;
                      }
                    case 2 :
                      {
                        tree = Read_User_Tree(cdata,mod,io);
                        break;
                      }
                    }

                  if(io->mod->s_opt->opt_topo == YES) Remove_Duplicates(cdata,io,tree);

                  if(io->fp_in_constraint_tree != NULL)
                    {
                      char *s;

                      PhyML_Printf("\n. Reading constraint tree file...");
                      
                      io->cstr_tree = Read_Tree_File_Phylip(io->fp_in_constraint_tree);

                      if(io->cstr_tree->n_root != NULL)
                        {
                          PhyML_Printf("\n== The constraint tree file must be unrooted");
                          Exit("\n");
                        }

                      s = Add_Taxa_To_Constraint_Tree(io->fp_in_constraint_tree,cdata);
                      fflush(NULL);
                      Free_Tree(tree);
                      tree = Read_Tree(&s);
                      io->in_tree = 2;
                      Free(s);
                      Check_Constraint_Tree_Taxa_Names(io->cstr_tree,cdata);
                      Alloc_Bip(io->cstr_tree);
                      Get_Bip(io->cstr_tree->a_nodes[0],
                              io->cstr_tree->a_nodes[0]->v[0],
                              io->cstr_tree);
                      if(tree->has_branch_lengths == NO) Add_BioNJ_Branch_Lengths(tree,cdata,mod,NULL);
                    }

                  
                  if(!tree) continue;
                  
                  time(&t_beg);
                  time(&(tree->t_beg));

                  tree->mod          = mod;
                  tree->io           = io;
                  tree->data         = cdata;
                  tree->n_pattern    = tree->data->crunch_len;
                  tree->n_root       = NULL;
                  tree->e_root       = NULL;
                  tree->n_tot_bl_opt = 0;

                  Set_Both_Sides(YES,tree);

                  if((!num_data_set) && (!num_tree) && (!num_rand_tree)) Check_Memory_Amount(tree);

                  if(io->cstr_tree && !Check_Topo_Constraints(tree,io->cstr_tree))
                    {
                      PhyML_Printf("\n\n== The initial tree does not satisfy the topological constraint.");
                      PhyML_Printf("\n== Please use the user input tree option with an adequate tree topology.");
                      Exit("\n");
                    }

                  Connect_CSeqs_To_Nodes(tree->data,tree->io,tree);
                  Make_Tree_For_Pars(tree);
                  Make_Tree_For_Lk(tree);
                  Make_Spr(tree);
                  Br_Len_Not_Involving_Invar(tree);
                  Unscale_Br_Len_Multiplier_Tree(tree);


#ifdef BEAGLE
                  if(mod->bootstrap == YES)
                    {
                      PhyML_Printf("\n== PhyML-BEAGLE does not support bootstrap analysis yet... ");
                      Exit("\n");
                    }
                  if(mod->ras->invar == YES)
                    {
                      PhyML_Printf("\n== PhyML-BEAGLE does not support invariant site models yet... ");
                      Exit("\n");
                    }
#endif

#ifdef PHYML
                  if(tree->io->print_json_trace == YES) JSON_Tree_Io(tree,tree->io->fp_out_json_trace); 

                  
                  Set_Update_Eigen(YES,tree->mod);                 
                  Lk(NULL,tree);
                  Set_Update_Eigen(NO,tree->mod);
                  
                  if(tree->mod->s_opt->opt_topo)
		    {
                      Global_Spr_Search(tree);                      
                      if(tree->n_root) Add_Root(tree->a_edges[0],tree);
                    }
                  else
                    {
#ifdef BEAGLE
                      tree->b_inst = create_beagle_instance(tree, io->quiet, io);
#endif
                      if(tree->mod->s_opt->opt_subst_param || tree->mod->s_opt->opt_bl) Round_Optimize(tree,ROUND_MAX);                      
                    }

                  
                  if(tree->mod->gamma_mgf_bl) Best_Root_Position_IL_Model(tree);

                  Set_Both_Sides(YES,tree);
                  Lk(NULL,tree);
                  Pars(NULL,tree);
                  Get_Tree_Size(tree);
                  PhyML_Printf("\n\n. Log likelihood of the current tree: %.*f.",DECIMAL_DIG,tree->c_lnL);

                  if(tree->io->ancestral == YES) Ancestral_Sequences(tree,YES);
                  
                  Check_Br_Lens(tree);
                  Br_Len_Involving_Invar(tree);
                  Rescale_Br_Len_Multiplier_Tree(tree);
                  
                  
#elif defined EVOLVE
                  Evolve(tree->data,tree->mod,tree);
                  Exit("\n. Exiting 'evolve'\n");
#endif

                  if(!tree->n_root) Get_Best_Root_Position(tree);

                  /* Print the tree estimated using the current random (or BioNJ) starting tree */
                  /* if(io->mod->s_opt->n_rand_starts > 1) */
                  if(orig_random_input_tree == YES)
                    {
                      Print_Tree(io->fp_out_trees,tree);
                      fflush(NULL);
                    }

                  /* Record the most likely tree in a string of characters */
                  if(tree->c_lnL > best_lnL)
                    {
                      best_lnL = tree->c_lnL;
                      if(most_likely_tree) Free(most_likely_tree);
                      most_likely_tree = Write_Tree(tree);

                      time(&t_end);
                      
                      Print_Fp_Out(io->fp_out_stats,t_beg,t_end,tree,
                                   io,num_data_set+1,
                                   (orig_random_input_tree == YES)?(num_rand_tree):(num_tree),
                                   (num_rand_tree == io->mod->s_opt->n_rand_starts-1)?(YES):(NO), io->precision);

                      if(tree->io->print_site_lnl) Print_Site_Lk(tree,io->fp_out_lk);
                    }



                  /* Start from BioNJ tree */
                  if((num_rand_tree == io->mod->s_opt->n_rand_starts-1) && (tree->mod->s_opt->random_input_tree))
                    {
                      /* Do one more iteration in the loop, but don't randomize the tree */
                      tree->mod->s_opt->n_rand_starts++;
                      tree->mod->s_opt->random_input_tree = NO;
                    }
#ifdef BEAGLE
                  finalize_beagle_instance(tree);
#endif
                  Free_Best_Spr(tree);
                  Free_Spr_List_One_Edge(tree);
                  Free_Spr_List_All_Edge(tree);
                  Free_Tree_Pars(tree);
                  Free_Tree_Lk(tree);
                  Free_Tree(tree);
                } //Tree done

              if(io->n_data_sets == 1) rewind(io->fp_out_tree);
              if(most_likely_tree) PhyML_Fprintf(io->fp_out_tree,"%s\n",most_likely_tree);


              /* Launch bootstrap analysis */
              if(io->do_boot || io->do_tbe)
                {
                  if(!io->quiet) PhyML_Printf("\n\n. Launch bootstrap analysis on the most likely tree...");

#ifdef MPI
                  MPI_Bcast (most_likely_tree, strlen(most_likely_tree)+1, MPI_CHAR, 0, MPI_COMM_WORLD);
                  if(!io->quiet) PhyML_Printf("\n\n. The bootstrap analysis will use %d CPU%c.",Global_numTask,Global_numTask>1?'s':'\0');
#endif
                  
                  most_likely_tree = Bootstrap_From_String(most_likely_tree,cdata,mod,io);

                  PhyML_Printf("\n\n. Completed the bootstrap analysis succesfully."); fflush(NULL);
                }
              else
                if(io->ratio_test != NO)
                  {
                    /* Launch aLRT */
                    most_likely_tree = aLRT_From_String(most_likely_tree,cdata,mod,io);
                  }


              /* Print the most likely tree in the output file */
              if(!io->quiet) PhyML_Printf("\n\n. Printing the most likely tree in file '%s'.", Basename(io->out_tree_file));
              if(io->n_data_sets == 1) rewind(io->fp_out_tree);
              
              t_tree *dum;
              dum = Read_Tree(&most_likely_tree);
              dum->data = cdata;
              dum->mod  = mod;
              dum->io   = io;
              Connect_CSeqs_To_Nodes(cdata,io,dum);
              Insert_Duplicates(dum);
              Free(most_likely_tree);
              most_likely_tree = Write_Tree(dum);
              Free_Tree(dum);
              
              PhyML_Fprintf(io->fp_out_tree,"%s\n",most_likely_tree);

              if(io->n_trees > 1 && io->n_data_sets > 1) break;
            }
            Free_Calign(cdata);
        }
    else
        {
          PhyML_Printf("\n== No data was found.\n");
          PhyML_Printf("\n== Err. in file %s at line %d\n",__FILE__,__LINE__);
          Exit("\n");
        }
      Free_Model_Complete(mod);
    }

  if(most_likely_tree) Free(most_likely_tree);

  if(mod->s_opt->n_rand_starts > 1) PhyML_Printf("\n. Best log likelihood: %f\n",best_lnL);

  Free_Optimiz(mod->s_opt);
  Free_Model_Basic(mod);

  if(io->fp_in_constraint_tree) fclose(io->fp_in_constraint_tree);
  if(io->fp_in_align)           fclose(io->fp_in_align);
  if(io->fp_in_tree)            fclose(io->fp_in_tree);
  if(io->fp_out_lk)             fclose(io->fp_out_lk);
  if(io->fp_out_tree)           fclose(io->fp_out_tree);
  if(io->fp_out_trees)          fclose(io->fp_out_trees);
  if(io->fp_out_stats)          fclose(io->fp_out_stats);
  if(io->fp_out_trace)          fclose(io->fp_out_trace);
  if(io->fp_out_json_trace)     fclose(io->fp_out_json_trace);

  if(io->fp_in_constraint_tree != NULL) Free_Tree(io->cstr_tree);
  Free_Input(io);

  time(&t_end);
  Print_Time_Info(t_beg,t_end);

#ifdef MPI
  MPI_Finalize();
#endif

  return 0;
}

#elif(M4)
#include "m4.h"
int main(int argc, char **argv)
{
  M4_main(argc, argv);
  return 1;
}

#elif(PART)
#include "mg.h"
int main(int argc, char **argv)
{
  PART_main(argc, argv);
  return 1;
}

/* #elif(PHYTIME) */
/* #include "times.h" */
/* int main(int argc, char **argv) */
/* { */
/*   TIMES_main(argc, argv); */
/*   return 1; */
/* } */

#elif(PHYCONT)
#include "continuous.h"
int main(int argc, char **argv)
{
  CONT_main(argc, argv);
  return 1;
}

#elif(RF)
int main(int argc, char **argv)
{
  option *io;
  int r_seed;


  io = (option *)Get_Input(argc,argv);
  if(!io) return(0);

  r_seed = (io->r_seed < 0)?(time(NULL)):(io->r_seed);
  srand(r_seed);
  io->r_seed = r_seed;
  
  Get_Seq(io);
  Shuffle_Sites(io->mod->ras->pinvar->v,io->data,io->n_otu);
  Print_Seq(stdout,io->data,io->n_otu);
  
  /* t_tree *tree1, *tree2; */
  /* FILE *fp_tree1, *fp_tree2; */
  /* int i,j; */

  /* fp_tree1 = (FILE *)fopen(argv[1],"r"); */
  /* fp_tree2 = (FILE *)fopen(argv[2],"r"); */

  /* tree1 = Read_Tree_File_Phylip(fp_tree1); */
  /* tree2 = Read_Tree_File_Phylip(fp_tree2); */


  /* Prune_Tree(tree1,tree2); */
  /* Prune_Tree(tree2,tree1); */
  

  /* PhyML_Printf("%s\n%s", */
  /*              Write_Tree(tree1), */
  /*              Write_Tree(tree2)); */
               
  
  /* Match_Nodes_In_Small_Tree(tree1,tree2); */

  /* For(i,2*tree1->n_otu-2) */
  /*   { */
  /*     printf("\n. Node %d in tree1 matches node %d in tree2",i,(tree1->noeud[i]->match_node)?(tree1->noeud[i]->match_node->num):(-1)); */
  /*   } */




/*   t_tree *tree1, *tree2; */
/*   FILE *fp_tree1, *fp_tree2; */
/*   int i,j,rf,n_edges,n_common,bip_size; */
/*   phydbl thresh; */
/*   t_edge *b; */


/*   fp_tree1 = (FILE *)fopen(argv[1],"r"); */
/*   fp_tree2 = (FILE *)fopen(argv[2],"r"); */
/*   thresh = (phydbl)atof(argv[3]); */

/*   tree1 = Read_Tree_File(fp_tree1); */
/*   tree2 = Read_Tree_File(fp_tree2); */

/*   Get_Rid_Of_Prefix('_',tree1); */

/* /\*   Find_Common_Tips(tree1,tree2); *\/ */

/*   Alloc_Bip(tree1); */
/*   Alloc_Bip(tree2); */

/*   Get_Bip(tree1->noeud[0],tree1->noeud[0]->v[0],tree1); */
/*   Get_Bip(tree2->noeud[0],tree2->noeud[0]->v[0],tree2); */

/* /\*   PhyML_Printf("\n. rf=%f\n",Compare_Bip_On_Existing_Edges(thresh,tree1,tree2)); *\/ */
/*   For(i,2*tree1->n_otu-3) tree1->a_edges[i]->bip_score = 0; */
/*   For(i,2*tree2->n_otu-3) tree2->a_edges[i]->bip_score = 0; */

/*   rf = 0; */
/*   n_edges = 0; */

/*   /\* First tree *\/ */
/*   For(i,2*tree1->n_otu-3)  */
/*     { */
/*       /\* Consider the branch only if the corresponding bipartition has size > 1 *\/ */
/*       b = tree1->a_edges[i]; */
/*       bip_size = MIN(b->left->bip_size[b->l_r],b->rght->bip_size[b->r_l]); */

/*       if(bip_size > 1) */
/* 	{ */
/* 	  /\* with non-zero length *\/ */
/* 	  if(tree1->a_edges[i]->l > thresh)   */
/* 	    { */
/* 	      n_edges++; */
/* 	      /\* This t_edge is not found in tree2 *\/ */
/* 	      if(!tree1->a_edges[i]->bip_score) rf++; ; */
/* 	    } */
/* 	} */
/*     } */


/*   /\* Second tree *\/ */
/*   For(i,2*tree2->n_otu-3)  */
/*     { */
/*       b = tree2->a_edges[i]; */
/*       bip_size = MIN(b->left->bip_size[b->l_r],b->rght->bip_size[b->r_l]); */

/*       if(bip_size > 1) */
/* 	{ */
/* 	  if(tree2->a_edges[i]->l > thresh)   */
/* 	    { */
/* 	      n_edges++; */
/* 	      /\* This t_edge is not found in tree1 *\/ */
/* 	      if(!tree2->a_edges[i]->bip_score) rf++; ; */
/* 	    } */
/* 	} */
/*     } */

/*   if(!n_edges) */
/*     { */
/*       Exit("\n. No comparable internal edges were found.\n"); */
/*     } */
/*   else */
/*     { */
/*       PhyML_Printf("\n. Robinson and Foulds distance: %f.",(double)rf/(n_edges)); */
/* /\*       PhyML_Printf("\n. %d internal edges were processed (%d in the first tree, %d in the second).\n",n_edges,n_edges_t1,n_edges-n_edges_t1); *\/ */
/*       PhyML_Printf("\n"); */
/*     } */

  return 1;
}

#elif(TIPORDER)
#include "tiporder.h"
int main(int argc, char **argv)
{
  TIPO_main(argc, argv);
  return 1;
}

#elif(TEST)
#include "xml.h"
int main(int argc, char **argv)
{
  option *io;
  int i;
  int year;
  
  io = (option *)Get_Input(argc,argv);
  if(!io) return(0);

  Get_Seq(io);

  for(i=0;i<io->n_otu;i++)
    {
      sscanf(io->data[i]->name,"%d",&year);
      PhyML_Printf("\n%s\t%d",io->data[i]->name,year);
    }

  /* for(i=0;i<io->n_otu;i++) */
  /*   { */
  /*     sscanf(io->data[i]->name,"%d",&year); */
  /*     PhyML_Printf("\n<clade id=\"clad%d\">",i+1); */
  /*     PhyML_Printf("\n\t<taxon value=\"%s\"/>",io->data[i]->name); */
  /*     PhyML_Printf("\n</clade>"); */
  /*     PhyML_Printf("\n<calibration id=\"cal%d\">",i+1); */
  /*     PhyML_Printf("\n\t<lower>%d</lower>",year); */
  /*     PhyML_Printf("\n\t<upper>%d</upper>",year); */
  /*     PhyML_Printf("\n\t<appliesto clade.id=\"clad%d\"/>",i+1); */
  /*     PhyML_Printf("\n</calibration>"); */
  /*   } */

  /* FILE *fp; */
  /* char *name,*date; */
  /* int i; */
  
  /* name = (char *)mCalloc(T_MAX_NAME,sizeof(char)); */
  /* date = (char *)mCalloc(T_MAX_NAME,sizeof(char)); */

  /* i = 0; */
  /* fp = Openfile(argv[1],READ); */
  /* do */
  /*   { */
  /*     if(fscanf(fp,"%s",name) == EOF) break; */
  /*     if(fscanf(fp,"%s",date) == EOF) break; */
  /*     PhyML_Printf("\n<clade id=\"clad%d\">",i+1); */
  /*     PhyML_Printf("\n\t<taxon value=\"%s\"/>",name); */
  /*     PhyML_Printf("\n</clade>"); */
  /*     PhyML_Printf("\n<calibration id=\"cal%d\">",i+1); */
  /*     PhyML_Printf("\n\t<lower>%s</lower>",date); */
  /*     PhyML_Printf("\n\t<upper>%s</upper>",date); */
  /*     PhyML_Printf("\n\t<appliesto clade.id=\"clad%d\"/>",i+1); */
  /*     PhyML_Printf("\n</calibration>"); */
  /*     ++i; */
  /*   } */
  /* while(1); */
  
}

#elif(INVITEE)
#include "invitee.h"
int main(int argc, char **argv)
{
  /*My_Function(argc, argv);*/
  /* PhyTime_XML(argc, argv); */
  Get_Input(argc,argv);
  return 1;
}

#elif(GEO)
#include "geo.h"
int main(int argc, char **argv)
{
  GEO_Main(argc,argv);
  return 1;
}

#elif(defined PHYREX || PHYREXSIM)
#include "phyrex.h"
int main(int argc, char **argv)
{
  PHYREX_Main(argc,argv);
  return 1;
}

#elif(PHYTIME)
#include "date.h"
int main(int argc, char **argv)
{
  DATE_Main(argc,argv);
  return 1;
}

#elif(CHECKPOINT)
#include "checkpoint.h"
int main(int argc, char **argv)
{
  CHECK_Main(argc,argv);
  return 1;
}

#endif

