

// PhyML:  a program that  computes maximum likelihood phylogenies from
// DNA or AA homologous sequences.

// Copyright (C) Stephane Guindon. Oct 2003 onward.

// All parts of the source except where indicated are distributed under
// the GNU public licence. See http://www.opensource.org for details.

#include "alrt.h"
#include "bionj.h"
#include "eigen.h"
#include "free.h"
#include "help.h"
#include "invitee.h"
#include "lk.h"
#include "mixt.h"
#include "models.h"
#include "optimiz.h"
#include "pars.h"
#include "simu.h"
#include "spr.h"
#include "utilities.h"

#ifdef MPI
#include "mpi_boot.h"
#endif

#ifdef BEAGLE
#include "beagle_utils.h"
#endif

#if (defined PHYML)
int main(int argc, char **argv)
{
  calign *cdata;
  option *io;
  t_tree *tree;
  int     num_data_set;
  int     num_tree, num_rand_tree;
  t_mod  *mod;
  time_t  t_beg, t_end;
  phydbl  best_lnL;
  int     r_seed;
  char   *most_likely_tree = NULL;
  int     orig_random_input_tree;

#ifdef MPI
  int rc;
  rc = MPI_Init(&argc, &argv);
  if (rc != MPI_SUCCESS)
  {
    PhyML_Fprintf(stderr, "\n. Err. starting MPI program. Terminating.\n");
    MPI_Abort(MPI_COMM_WORLD, rc);
  }
  if (MPI_Comm_size(MPI_COMM_WORLD, &Global_numTask) != MPI_SUCCESS)
    MPI_Abort(MPI_COMM_WORLD, rc);
  if (MPI_Comm_rank(MPI_COMM_WORLD, &Global_myRank) != MPI_SUCCESS)
    MPI_Abort(MPI_COMM_WORLD, rc);
  PhyML_Fprintf(stdout, "\n\n. Running the analysis on %d CPU%s.",
                Global_numTask, Global_numTask > 1 ? "s." : ".");
#endif

#ifdef QUIET
  setvbuf(stdout, NULL, _IOFBF, 2048);
#endif

  tree     = NULL;
  mod      = NULL;
  best_lnL = UNLIKELY;

  io = (option *)Get_Input(argc, argv);
  if (!io)
    return (0);
  else if (io->use_xml == YES)
  {
    Free(io);
    return (0);
  }

  r_seed = (io->r_seed < 0) ? (time(NULL)) : (io->r_seed);
#ifdef MPI
  srand(r_seed + Global_myRank);
#else
  srand(r_seed);
#endif
  io->r_seed = r_seed;

  if (io->in_tree == 2)
    Test_Multiple_Data_Set_Format(io);
  else
    io->n_trees = 1;

  if (io->n_trees == 0 && io->in_tree == 2)
  {
    PhyML_Printf("\n== Err.: the input tree file does not provide a tree in "
                 "valid format.");
    Exit("\n");
  }

  if ((io->n_data_sets > 1) && (io->n_trees > 1))
  {
    io->n_data_sets = MIN(io->n_trees, io->n_data_sets);
    io->n_trees     = MIN(io->n_trees, io->n_data_sets);
  }

  for (num_data_set = 0; num_data_set < io->n_data_sets; num_data_set++)
  {
    best_lnL = UNLIKELY;
    Get_Seq(io);
    Make_Model_Complete(io->mod);
    Set_Model_Name(io->mod);
    Print_Settings(io);
    mod                    = io->mod;
    orig_random_input_tree = io->mod->s_opt->random_input_tree;

    if (io->data)
    {
      if (io->n_data_sets > 1)
        PhyML_Printf("\n. Data set [#%d]\n", num_data_set + 1);
      cdata = Compact_Data(io->data, io);

      Free_Seq(io->data, cdata->n_otu);

      for (num_tree = (io->n_trees == 1) ? (0) : (num_data_set);
           num_tree < io->n_trees; num_tree++)
      {
        if (io->mod->s_opt->random_input_tree == NO)
          io->mod->s_opt->n_rand_starts = 1;

        if (orig_random_input_tree == YES && io->n_trees > 1)
        {
          PhyML_Printf("\n== Cannot combine random starting trees with "
                       "multiple input trees.");
          Exit("\n");
        }

        for (num_rand_tree = 0; num_rand_tree < io->mod->s_opt->n_rand_starts;
             num_rand_tree++)
        {
          if ((io->mod->s_opt->random_input_tree) &&
              (io->mod->s_opt->topo_search != NNI_MOVE))
            if (!io->quiet)
              PhyML_Printf("\n\n. [Random start %3d/%3d]", num_rand_tree + 1,
                           io->mod->s_opt->n_rand_starts);

          Init_Model(cdata, mod, io);
          Set_Model_Parameters(mod);

#ifdef M4
          if (io->mod->use_m4mod) M4_Init_Model(mod->m4mod, cdata, mod);
#endif

          // Make the initial tree
          switch (io->in_tree)
          {
          case 0:
          case 1:
          {
            tree = Dist_And_BioNJ(cdata, mod, io);
            if (io->print_mat_and_exit == YES)
            {
              Print_Mat(ML_Dist(cdata, mod));
              exit(-1);
            }
            break;
          }
          case 2:
          {
            tree = Read_User_Tree(cdata, mod, io);
            break;
          }
          }

          if (io->mod->s_opt->opt_topo == YES)
            Remove_Duplicates(cdata, io, tree);

          if (io->fp_in_constraint_tree != NULL)
          {
            char *s;

            PhyML_Printf("\n. Reading constraint tree file...");

            io->cstr_tree = Read_Tree_File_Phylip(io->fp_in_constraint_tree);

            if (io->cstr_tree->n_root != NULL)
            {
              PhyML_Printf("\n== The constraint tree file must be unrooted");
              Exit("\n");
            }

            s = Add_Taxa_To_Constraint_Tree(io->fp_in_constraint_tree, cdata);
            fflush(NULL);
            Free_Tree(tree);
            tree        = Read_Tree(&s);
            io->in_tree = 2;
            Free(s);
            Check_Constraint_Tree_Taxa_Names(io->cstr_tree, cdata);
            Alloc_Bip(io->cstr_tree);
            Get_Bip(io->cstr_tree->a_nodes[0], io->cstr_tree->a_nodes[0]->v[0],
                    io->cstr_tree);
            if (tree->has_branch_lengths == NO)
              Add_BioNJ_Branch_Lengths(tree, cdata, mod, NULL);
          }

          if (!tree) continue;

          time(&t_beg);
          time(&(tree->t_beg));

          tree->mod          = mod;
          tree->io           = io;
          tree->data         = cdata;
          tree->n_root       = NULL;
          tree->e_root       = NULL;
          tree->n_tot_bl_opt = 0;

          Set_Both_Sides(YES, tree);

          if ((!num_data_set) && (!num_tree) && (!num_rand_tree))
            Check_Memory_Amount(tree);

          if (io->cstr_tree && !Check_Topo_Constraints(tree, io->cstr_tree))
          {
            PhyML_Printf("\n\n== The initial tree does not satisfy the "
                         "topological constraint.");
            PhyML_Printf("\n== Please use the user input tree option with an "
                         "adequate tree topology.");
            Exit("\n");
          }

          Connect_CSeqs_To_Nodes(tree->data, tree->io, tree);
          Make_Tree_For_Pars(tree);
          Make_Tree_For_Lk(tree);
          Make_Spr(tree);
          Br_Len_Not_Involving_Invar(tree);
          Unscale_Br_Len_Multiplier_Tree(tree);

#ifdef BEAGLE
          if (mod->bootstrap == YES)
          {
            PhyML_Printf("\n== PhyML-BEAGLE does not support bootstrap "
                         "analysis yet... ");
            Exit("\n");
          }
          if (mod->ras->invar == YES)
          {
            PhyML_Printf("\n== PhyML-BEAGLE does not support invariant site "
                         "models yet... ");
            Exit("\n");
          }
#endif

          if (tree->io->print_json_trace == YES)
            JSON_Tree_Io(tree, tree->io->fp_out_json_trace);

          Set_Update_Eigen(YES, tree->mod);
          Lk(NULL, tree);
          Set_Update_Eigen(NO, tree->mod);

          PhyML_Printf("\n. Init log-likelihood: %f", tree->c_lnL);

          if (tree->mod->s_opt->opt_topo)
          {
            Global_Spr_Search(tree);
            if (tree->n_root) Add_Root(tree->a_edges[0], tree);
          }
          else
          {
#ifdef BEAGLE
            tree->b_inst = create_beagle_instance(tree, io->quiet, io);
#endif
            if (tree->mod->s_opt->opt_subst_param ||
                tree->mod->s_opt->opt_bl_one_by_one)
              Round_Optimize(tree, ROUND_MAX);
          }

        
          /* if(tree->mod->gamma_mgf_bl)
           * Best_Root_Position_IL_Model(tree); */

          Set_Both_Sides(YES, tree);
          Lk(NULL, tree);
          Pars(NULL, tree);
          Get_Tree_Size(tree);
          PhyML_Printf("\n\n. Log likelihood of the current tree: %.*f.",
                       DECIMAL_DIG, tree->c_lnL);

          if (tree->io->ancestral == YES) Ancestral_Sequences(tree, YES);

          Check_Br_Lens(tree);
          Br_Len_Involving_Invar(tree);
          Rescale_Br_Len_Multiplier_Tree(tree);

          if (!tree->n_root) Get_Best_Root_Position(tree);

          /* Print the tree estimated using the current random (or
           * BioNJ) starting tree */
          /* if(io->mod->s_opt->n_rand_starts > 1) */
          if (orig_random_input_tree == YES)
          {
            Print_Tree(io->fp_out_trees, tree);
            fflush(NULL);
          }

          /* Record the most likely tree in a string of characters */
          if (tree->c_lnL > best_lnL)
          {
            best_lnL = tree->c_lnL;
            if (most_likely_tree) Free(most_likely_tree);
            most_likely_tree = Write_Tree(tree);

            time(&t_end);

            Print_Fp_Out(
                io->fp_out_stats, t_beg, t_end, tree, io, num_data_set + 1,
                (orig_random_input_tree == YES) ? (num_rand_tree) : (num_tree),
                (num_rand_tree == io->mod->s_opt->n_rand_starts - 1) ? (YES)
                                                                     : (NO),
                io->precision);

            if (tree->io->print_site_lnl) Print_Site_Lk(tree, io->fp_out_lk);
          }

          /* Start from BioNJ tree */
          if ((num_rand_tree == io->mod->s_opt->n_rand_starts - 1) &&
              (tree->mod->s_opt->random_input_tree))
          {
            /* Do one more iteration in the loop, but don't randomize the tree
             */
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
        } // Tree done

        if (io->n_data_sets == 1) rewind(io->fp_out_tree);
        if (most_likely_tree)
          PhyML_Fprintf(io->fp_out_tree, "%s\n", most_likely_tree);

        /* Launch bootstrap analysis */
        if (io->do_boot || io->do_tbe)
        {
          if (!io->quiet)
            PhyML_Printf(
                "\n\n. Launch bootstrap analysis on the most likely tree...");

#ifdef MPI
          MPI_Bcast(most_likely_tree, strlen(most_likely_tree) + 1, MPI_CHAR, 0,
                    MPI_COMM_WORLD);
          if (!io->quiet)
            PhyML_Printf("\n\n. The bootstrap analysis will use %d CPU%c.",
                         Global_numTask, Global_numTask > 1 ? 's' : '\0');
#endif

          most_likely_tree =
              Bootstrap_From_String(most_likely_tree, cdata, mod, io);

          PhyML_Printf("\n\n. Completed the bootstrap analysis succesfully.");
          fflush(NULL);
        }
        else if (io->ratio_test != NO)
        {
          /* Launch aLRT */
          most_likely_tree = aLRT_From_String(most_likely_tree, cdata, mod, io);
        }

        /* Print the most likely tree in the output file */
        if (!io->quiet)
          PhyML_Printf("\n\n. Printing the most likely tree in file '%s'.",
                       Basename(io->out_tree_file));
        if (io->n_data_sets == 1) rewind(io->fp_out_tree);

        t_tree *dum;
        dum       = Read_Tree(&most_likely_tree);
        dum->data = cdata;
        dum->mod  = mod;
        dum->io   = io;
        Connect_CSeqs_To_Nodes(cdata, io, dum);
        Insert_Duplicates(dum);
        Free(most_likely_tree);
        most_likely_tree = Write_Tree(dum);
        Free_Tree(dum);

        PhyML_Fprintf(io->fp_out_tree, "%s\n", most_likely_tree);

        if (io->n_trees > 1 && io->n_data_sets > 1) break;
      }
      Free_Calign(cdata);
    }
    else
    {
      PhyML_Printf("\n== No data was found.\n");
      PhyML_Printf("\n== Err. in file %s at line %d\n", __FILE__, __LINE__);
      Exit("\n");
    }
    Free_Model_Complete(mod);
  }

  if (most_likely_tree) Free(most_likely_tree);

  if (mod->s_opt->n_rand_starts > 1)
    PhyML_Printf("\n. Best log likelihood: %f\n", best_lnL);

  Free_Optimiz(mod->s_opt);
  Free_Model_Basic(mod);

  if (io->fp_in_constraint_tree) fclose(io->fp_in_constraint_tree);
  if (io->fp_in_align) fclose(io->fp_in_align);
  if (io->fp_in_tree) fclose(io->fp_in_tree);
  if (io->fp_out_lk) fclose(io->fp_out_lk);
  if (io->fp_out_tree) fclose(io->fp_out_tree);
  if (io->fp_out_trees) fclose(io->fp_out_trees);
  if (io->fp_out_stats) fclose(io->fp_out_stats);
  if (io->fp_out_trace) fclose(io->fp_out_trace);
  if (io->fp_out_json_trace) fclose(io->fp_out_json_trace);

  if (io->fp_in_constraint_tree != NULL) Free_Tree(io->cstr_tree);
  Free_Input(io);

  time(&t_end);
  Print_Time_Info(t_beg, t_end);

#ifdef MPI
  MPI_Finalize();
#endif

  return 0;
}

#elif (M4)
#include "m4.h"
int main(int argc, char **argv)
{
  M4_main(argc, argv);
  return 1;
}

#elif (PART)
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

#elif (PHYCONT)
#include "continuous.h"
int main(int argc, char **argv)
{
  CONT_main(argc, argv);
  return 1;
}

#elif (RF)
int main(int argc, char **argv)
{
  // option *io;
  // int     r_seed;

  // io = (option *)Get_Input(argc, argv);
  // if (!io) return (0);

  // r_seed = (io->r_seed < 0) ? (time(NULL)) : (io->r_seed);
  // srand(r_seed);
  // io->r_seed = r_seed;

  // Get_Seq(io);
  // Shuffle_Sites(io->mod->ras->pinvar->v, io->data, io->n_otu);
  // Print_Seq(stdout, io->data, io->n_otu);

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
  /*     printf("\n. Node %d in tree1 matches node %d in
   * tree2",i,(tree1->noeud[i]->match_node)?(tree1->noeud[i]->match_node->num):(-1));
   */
  /*   } */

  t_tree *tree1, *tree2;
  FILE *fp_tree1, *fp_tree2;
  int i,j,rf,n_edges,n_common,bip_size;
  phydbl thresh;
  t_edge *b; 

  fp_tree1 = (FILE *)fopen(argv[1],"r");
  fp_tree2 = (FILE *)fopen(argv[2],"r");

  tree1 = Read_Tree_File_Phylip(fp_tree1); 
  tree2 = Read_Tree_File_Phylip(fp_tree2); 

  // Get_Rid_Of_Prefix('_',tree1);
  // Find_Common_Tips(tree1,tree2); 

  Alloc_Bip(tree1); 
  Alloc_Bip(tree2); 

  Get_Bip(tree1->a_nodes[0],tree1->a_nodes[0]->v[0],tree1); 
  Get_Bip(tree2->a_nodes[0],tree2->a_nodes[0]->v[0],tree2);

  PhyML_Printf("\n. rf=%f\n",
               .5 * (Compare_Bip(tree1, tree2, NO) + Compare_Bip(tree2, tree1, NO)));

  // for (i = 0; i < 2 * tree1->n_otu - 3; ++i) tree1->a_edges[i]->bip_score = 0;
  // for (i = 0; i < 2 * tree2->n_otu - 3; ++i) tree2->a_edges[i]->bip_score = 0;

  /*   rf = 0; */
  /*   n_edges = 0; */

  /*   /\* First tree *\/ */
  /*   For(i,2*tree1->n_otu-3)  */
  /*     { */
  /*       /\* Consider the branch only if the corresponding bipartition has
   * size > 1 *\/ */
  /*       b = tree1->a_edges[i]; */
  /*       bip_size = MIN(b->left->bip_size[b->l_r],b->rght->bip_size[b->r_l]);
   */

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
  /*       bip_size = MIN(b->left->bip_size[b->l_r],b->rght->bip_size[b->r_l]);
   */

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
  /*       PhyML_Printf("\n. Robinson and Foulds distance:
   * %f.",(double)rf/(n_edges)); */
  /* /\*       PhyML_Printf("\n. %d internal edges were processed (%d in the
   * first tree, %d in the second).\n",n_edges,n_edges_t1,n_edges-n_edges_t1);
   * *\/ */
  /*       PhyML_Printf("\n"); */
  /*     } */

  return 1;
}

#elif (TIPORDER)
#include "tiporder.h"
int main(int argc, char **argv)
{
  TIPO_main(argc, argv);
  return 1;
}

#elif (TEST)
#include "xml.h"
int main(int argc, char **argv)
{

  //   /* Prediction using linear extrapolation of velocities estimated at the
  //   tips. Takes as input */
  //   /* the XML file used for running the MCMC analysis (so as to get the
  //   coordinates + sampling dates), */
  //   /* the list of sampled trees, the time interval for which prediction will
  //   be made and the burnin */
  //   /* proportion */

  //    option *io;
  //    t_tree *tree;
  //    xml_node *root;
  //    phydbl date_old,date_recent,burnin;
  //    short int first;
  //    int stepsize,sampsize;

  //    date_recent = atof(argv[4]);
  //    date_old = atof(argv[5]);
  //    burnin = atof(argv[6]);
  //    sampsize = 1000;
  //    first = YES;

  //    PhyML_Printf("\n. Time interval considered:
  //    [%f,%f]",date_old,date_recent); assert(date_old < date_recent);

  //    io = (option *)Get_Input(argc,argv);
  //    if(!io) return(0);

  //    Read_User_Tree(NULL,NULL,io);

  //    root = XML_Load_File(io->fp_in_xml);

  //    stepsize = io->treelist->list_size*(1. - burnin) / sampsize;

  //    for(int i=0;i<io->treelist->list_size;i++)
  //      {
  //        if(i > (int)(burnin*io->treelist->list_size))
  //          {
  //            PhyML_Fprintf(stderr,"\n. Processing tree %d",i+1);
  //            tree = io->treelist->tree[i];

  //            tree->times = TIMES_Make_Time_Struct(tree->n_otu);
  //            TIMES_Init_Time_Struct(tree->times,NULL,tree->n_otu);

  //            tree->rates = RATES_Make_Rate_Struct(tree->n_otu);
  //            RATES_Init_Rate_Struct(tree->rates,NULL,tree->n_otu);

  //            tree->mmod = PHYREX_Make_Migrep_Model(tree->n_otu,2);
  //            tree->mmod->n_dim = 2;
  //            PHYREX_Set_Default_Migrep_Mod(tree->n_otu,tree->mmod);

  //            XML_Read_Calibration(root,tree);
  //            MIXT_Chain_Cal(tree);

  //            TIMES_Randomize_Tip_Times_Given_Calibrations(tree); // Topology
  //            is unchanged tree->rates->clock_r = 1.0;
  //            TIMES_Bl_To_Times(1,tree);
  //            Update_Ancestors(tree->n_root,tree->n_root->v[2],tree->n_root->b[2],tree);
  //            Update_Ancestors(tree->n_root,tree->n_root->v[1],tree->n_root->b[1],tree);

  //            PHYREX_Make_And_Connect_Tip_Disks(tree);
  //            PHYREX_Tree_To_Ldsk(tree);

  //            Node_Labels_To_Velocities(tree);
  //            Node_Labels_To_Locations(tree);

  //            if(first == YES) PhyML_Printf("\n. XXX Tree\t Time\t Tax\t
  //            Longitude\t Latitude\t NextLongitude\t NextLatitude"); first =
  //            NO;

  //            for(int j=0;j<2*tree->n_otu-1;++j)
  //              {
  //                if(tree->a_nodes[j]->tax == YES)
  //                  {
  //                    phydbl survival_duration,survival_rate,t,delta_t;
  //                    int k,K;

  //                    survival_rate = 1.0;
  //                    K = 4;
  //                    delta_t = (date_recent - date_old) / K;

  //                    survival_duration = Rexp(survival_rate);

  //                    if(tree->times->nd_t[j] + survival_duration > date_old)
  //                      {
  //                        t = date_old;
  //                        k = 0;
  //                        do
  //                          {
  //                            PhyML_Printf("\n. XXX %d\t %f\t %d\t %f\t %f\t
  //                            %f\t %f",
  //                                         i+1,
  //                                         tree->times->nd_t[j],
  //                                         tree->a_nodes[j]->tax,
  //                                         tree->a_nodes[j]->ldsk->coord->lonlat[0],
  //                                         tree->a_nodes[j]->ldsk->coord->lonlat[1],
  //                                         tree->a_nodes[j]->ldsk->coord->lonlat[0]
  //                                         + (k*delta_t) *
  //                                         tree->a_nodes[j]->ldsk->veloc->deriv[0],
  //                                         tree->a_nodes[j]->ldsk->coord->lonlat[1]
  //                                         + (k*delta_t) *
  //                                         tree->a_nodes[j]->ldsk->veloc->deriv[1]);
  //                            k++;
  //                            t += delta_t;
  //                          }
  //                        while(t < MIN(date_recent,tree->times->nd_t[j] +
  //                        survival_duration));
  //                      }
  //                  }
  //              }

  //            PHYREX_Free_Ldsk_Struct(tree);
  //            RATES_Free_Rates(tree->rates);
  //            TIMES_Free_Times(tree->times);
  //            Free_Mmod(tree->mmod);
  //            Free_Tree(tree);

  //          }
  //      }

  //    Exit("\n");

  // /* Select samples in the 801 sequence WNV data set within a user defined
  // date range ./test --xml=../WNV_RRW_tree_1.xml 2000 1990 1.0 */ option *io;
  // xml_node
  // *beast_root,*beast_taxa,*beast_taxon,*beast_date,*beast_loc,*beast_seq,**valid_taxa;
  // xml_node *phyrex_root,*nd,*ndnd;
  // int i,n_selected,n_selected_max,r_seed;
  // FILE *fp_xml,*fp_coord,*fp_seq;
  // char *filename,*dum,*xml_filename;
  // phydbl date_recent, date_f;
  // int *permut;
  // phydbl *select_proba,sum;

  // valid_taxa = (xml_node **)mCalloc(800,sizeof(xml_node *));
  // select_proba = (phydbl *)mCalloc(800,sizeof(phydbl));

  // filename = (char *)mCalloc(100,sizeof(char));
  // xml_filename = (char *)mCalloc(100,sizeof(char));
  // r_seed = time(NULL);
  // PhyML_Printf("\n seed: %d",r_seed);
  // srand(r_seed);

  // io = (option *)Get_Input(argc,argv);
  // if(!io) return(0);

  // date_recent = atof(argv[2]);
  // n_selected_max = atoi(argv[3]);

  // dum = (char *)mCalloc(100,sizeof(char));
  // sprintf(dum,"%s%d%s","wnv_config_",r_seed,".xml");
  // strcpy(xml_filename,dum);
  // Free(dum);

  // PhyML_Printf("\n. Upper limit of time: %f",date_recent);

  // strcpy(filename,"coord.txt");
  // fp_coord = Openfile(filename,WRITE);

  // PhyML_Fprintf(fp_coord,"traits lat long\n");
  // PhyML_Fprintf(fp_coord,"|NorthEast| 1000 1000\n");
  // PhyML_Fprintf(fp_coord,"|SouthWest| -1000 -1000");

  // strcpy(filename,"seq.txt");
  // fp_seq = Openfile(filename,WRITE);

  // dum = (char *)mCalloc(100,sizeof(char));
  // sprintf(dum,"%s%d","wnv_predict_",r_seed);
  // phyrex_root = Generate_PhyREX_XMLObj("ibm",dum,"seq.txt","coord.txt");
  // Free(dum);

  // beast_root = XML_Load_File(io->fp_in_xml);

  // beast_taxa = XML_Search_Node_Name("taxa",NO,beast_root);

  // i = 1;
  // n_selected = 0;
  // beast_taxon = beast_taxa->child;
  // do
  //   {
  //     beast_date = XML_Search_Node_Name("date",NO,beast_taxon);
  //     assert(beast_date);

  //     date_f = atof(XML_Get_Attribute_Value(beast_date,"value"));

  //     if(date_f < date_recent)
  //       {
  //         valid_taxa[n_selected] = beast_taxon;
  //         select_proba[n_selected] = exp(-fabs(date_f-date_recent));
  //         n_selected++;
  //       }

  //     beast_taxon = beast_taxon->next;
  //   }
  // while(beast_taxon);

  // PhyML_Fprintf(fp_seq,"%d 10302\n",MIN(n_selected,n_selected_max));
  // permut = Permutate(MIN(n_selected,n_selected_max));

  // sum = 0.0;
  // for(i=0;i<n_selected;++i) sum += select_proba[i];
  // for(i=0;i<n_selected;++i) select_proba[i] /= sum;

  // for(i=0;i<n_selected;++i) PhyML_Printf("\n. Target %s (%d) prob:
  // %f",XML_Get_Attribute_Value(valid_taxa[i],"id"),i,select_proba[i]);

  // for(i=0;i<MIN(n_selected,n_selected_max);++i)
  //   {
  //     permut[i] = Sample_i_With_Proba_pi(select_proba,n_selected);

  //     beast_date = XML_Search_Node_Name("date",NO,valid_taxa[permut[i]]);
  //     assert(beast_date);
  //     date_f = atof(XML_Get_Attribute_Value(beast_date,"value"));

  //     nd = XML_Add_Node(phyrex_root,"clade");
  //     dum = (char *)mCalloc(100,sizeof(char));
  //     sprintf(dum,"%s%d","clad",i+1);
  //     nd->attr = XML_Make_Attribute(NULL,"id",dum);
  //     ndnd = XML_Add_Node(nd,"taxon");
  //     ndnd->attr =
  //     XML_Make_Attribute(NULL,"value",XML_Get_Attribute_Value(valid_taxa[permut[i]],"id"));

  //     nd = XML_Add_Node(phyrex_root,"calibration");
  //     dum = (char *)mCalloc(100,sizeof(char));
  //     sprintf(dum,"%s%d","cal",i+1);
  //     nd->attr = XML_Make_Attribute(NULL,"id",dum);
  //     Free(dum);

  //     ndnd = XML_Add_Node(nd,"lower");
  //     XML_Set_Node_Value(ndnd,XML_Get_Attribute_Value(beast_date,"value"));

  //     ndnd = XML_Add_Node(nd,"upper");
  //     XML_Set_Node_Value(ndnd,XML_Get_Attribute_Value(beast_date,"value"));

  //     ndnd = XML_Add_Node(nd,"appliesto");
  //     dum = (char *)mCalloc(100,sizeof(char));
  //     sprintf(dum,"%s%d","clad",i+1);
  //     ndnd->attr = XML_Make_Attribute(NULL,"clade.id",dum);
  //     Free(dum);

  //     PhyML_Fprintf(fp_coord,"\n");
  //     PhyML_Fprintf(fp_coord,"
  //     %s",XML_Get_Attribute_Value(valid_taxa[permut[i]],"id"));

  //     beast_loc =
  //     XML_Search_Node_Generic("attr","name","latitude",NO,valid_taxa[permut[i]]);
  //     assert(beast_loc);
  //     PhyML_Fprintf(fp_coord,"\t %s",beast_loc->value);

  //     beast_loc =
  //     XML_Search_Node_Generic("attr","name","longitude",NO,valid_taxa[permut[i]]);
  //     assert(beast_loc);
  //     PhyML_Fprintf(fp_coord,"\t %s",beast_loc->value);

  //     beast_seq =
  //     XML_Search_Node_Attribute_Value("idref",XML_Get_Attribute_Value(valid_taxa[permut[i]],"id"),NO,beast_root);
  //     assert(beast_seq);
  //     PhyML_Printf("\n. Found match with %s (%d) prob:
  //     %f",XML_Get_Attribute_Value(valid_taxa[permut[i]],"id"),permut[i],select_proba[permut[i]]);
  //     PhyML_Fprintf(fp_seq,"%s",XML_Get_Attribute_Value(valid_taxa[permut[i]],"id"));
  //     PhyML_Fprintf(fp_seq,"\t\t%s\n",beast_seq->value);

  //     select_proba[permut[i]] = 0.0;
  //     sum = 0.0;
  //     for(int i=0;i<n_selected;++i) sum += select_proba[i];
  //     for(int i=0;i<n_selected;++i) select_proba[i] /= sum;

  //   }

  // PhyML_Printf("\n. Number of selected sequences:
  // %d",MIN(n_selected,n_selected_max));

  // fp_xml = Openfile(xml_filename,WRITE);
  // XML_Write_XML_Graph(fp_xml,phyrex_root);
  // fclose(fp_xml);

  // fclose(fp_coord);
  // fclose(fp_seq);

  // Philippe Lemey data preprocessisng steps
  option *io;
  int     i;
  char   *year, *lat, *lon;
  t_tree *tree;

  io = (option *)Get_Input(argc, argv);
  if (!io) return (0);

  tree = Read_User_Tree(NULL, NULL, io);

  PhyML_Printf("\n. n_otus: %d", tree->n_otu);
  for (i = 0; i < tree->n_otu; ++i)
  {
    PhyML_Printf("\n<clade id=\"clad%d\">", i + 1);
    PhyML_Printf("\n\t<taxon value=\"%s\"/>", tree->a_nodes[i]->name);
    PhyML_Printf("\n</clade>");
    PhyML_Printf("\n<calibration id=\"cal%d\">", i + 1);
    PhyML_Printf("\n\t<lower>-%f</lower>",
                 atof(strrchr(tree->a_nodes[i]->name, '|') + 1));
    PhyML_Printf("\n\t<upper>-%f</upper>",
                 atof(strrchr(tree->a_nodes[i]->name, '|') + 1));
    PhyML_Printf("\n\t<appliesto clade.id=\"clad%d\"/>", i + 1);
    PhyML_Printf("\n</calibration>");
  }

  PhyML_Printf("\n\n%d\t 4", tree->n_otu);
  for (i = 0; i < tree->n_otu; ++i)
  {
    PhyML_Printf("\n%s\tATGC", tree->a_nodes[i]->name);
  }
  Exit("\n");

  //    Get_Seq(io);

  //    for(i=0;i<io->n_otu;i++)
  //       {
  //         sscanf(io->data[i]->name,"%d",&year);
  //         PhyML_Printf("\n%s\t%d",io->data[i]->name,year);
  //       }

  for (i = 0; i < io->n_otu; i++)
  {
    //    sscanf(io->data[i]->name,"%[^_]_%[^_]_%lf",a,b,&year);
    //    PhyML_Printf("\n. a: %s b: %s",a,b);
    strtok(io->data[i]->name, "_");
    char *a = strtok(NULL, "_");
    year    = strtok(NULL, "_");
    lat     = strtok(NULL, "_");
    lon     = strtok(NULL, "_");

    PhyML_Printf("\n\t <taxon id=\"%s_%s_%s\">", io->data[i]->name, a, year);
    PhyML_Printf(
        "\n\t\t <date value=\"%s\" direction=\"forwards\" units=\"years\"/>",
        year);
    PhyML_Printf("\n\t\t <attr name=\"lat\">");
    PhyML_Printf("\n\t\t\t %s", lat);
    PhyML_Printf("\n\t\t </attr>");
    PhyML_Printf("\n\t\t <attr name=\"lon\">");
    PhyML_Printf("\n\t\t\t %s", lon);
    PhyML_Printf("\n\t\t </attr>");
    PhyML_Printf("\n\t\t <attr name=\"location\">");
    PhyML_Printf("\n\t\t\t %s %s", lat, lon);
    PhyML_Printf("\n\t\t </attr>");
    PhyML_Printf("\n\t </taxon>");

    // PhyML_Printf("\n<clade id=\"clad%d\">",i+1);
    // PhyML_Printf("\n\t<taxon value=\"%s\"/>",io->data[i]->name);
    // PhyML_Printf("\n</clade>");
    // PhyML_Printf("\n<calibration id=\"cal%d\">",i+1);
    // PhyML_Printf("\n\t<lower>%d</lower>",year);
    // PhyML_Printf("\n\t<upper>%d</upper>",year);
    // PhyML_Printf("\n\t<appliesto clade.id=\"clad%d\"/>",i+1);
    // PhyML_Printf("\n</calibration>");
  }

  //  FILE *fp;
  //  char *name, *date;
  //  int i;

  //  name = (char *)mCalloc(T_MAX_NAME, sizeof(char));
  //  date = (char *)mCalloc(T_MAX_NAME, sizeof(char));

  //  i = 0;
  //  fp = Openfile(argv[1], READ);
  //  do
  //  {
  //      if (fscanf(fp, "%s", name) == EOF)
  //          break;
  //      if (fscanf(fp, "%s", date) == EOF)
  //          break;
  //      PhyML_Printf("\n<clade id=\"clad%d\">", i + 1);
  //      PhyML_Printf("\n\t<taxon value=\"%s\"/>", name);
  //      PhyML_Printf("\n</clade>");
  //      PhyML_Printf("\n<calibration id=\"cal%d\">", i + 1);
  //      PhyML_Printf("\n\t<lower>%s</lower>", date);
  //      PhyML_Printf("\n\t<upper>%s</upper>", date);
  //      PhyML_Printf("\n\t<appliesto clade.id=\"clad%d\"/>", i + 1);
  //      PhyML_Printf("\n</calibration>");
  //      ++i;
  //  } while (1);
}

#elif (INVITEE)
#include "invitee.h"
int main(int argc, char **argv)
{
  // My_Function(argc, argv);
  // PhyTime_XML(argc, argv);
  // Get_Input(argc, argv);
  // return 1;
}

#elif (GEO)
#include "geo.h"
int main(int argc, char **argv)
{
  GEO_Main(argc, argv);
  return 1;
}

#elif (defined PHYREX || PHYREXSIM)
#include "phyrex.h"
int main(int argc, char **argv)
{
  PHYREX_Main(argc, argv);
  return 1;
}

#elif (PHYTIME)
#include "date.h"
int main(int argc, char **argv)
{
  DATE_Main(argc, argv);
  return 1;
}

#elif (CHECKPOINT)
#include "checkpoint.h"
int main(int argc, char **argv)
{
  CHECK_Main(argc, argv);
  return 1;
}

#elif (EVOLVE)
#include "evolve.h"
int main(int argc, char **argv)
{
  EVOLVE_Main(argc, argv);
  return 1;
}

#endif
