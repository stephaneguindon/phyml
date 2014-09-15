/*

PhyML:  a program that  computes maximum likelihood phyLOGenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include "mpi_boot.h"

/* #ifdef MPI */

/*********************************************************/

void Bootstrap_MPI(t_tree *tree)
{
  int *site_num, n_site;
  int replicate,j,k;
  int position, init_len, nbRep;

  calign *boot_data;
  t_tree *boot_tree;
  t_mod *boot_mod;
  matrix *boot_mat;
  char *s;

  MPI_Status Stat;
  int randomRecv, bootRecv, nbElem, i;
  int *score_par, *score_tot;
  char *bootStr, *t;
  
  randomRecv = nbElem = bootRecv = 0;

  tree->print_boot_val       = 1;
  tree->print_alrt_val       = 0;
  boot_tree                  = NULL;

  site_num = (int *)mCalloc(tree->data->init_len,sizeof(int));
  
  Free_Bip(tree);
  Alloc_Bip(tree);
  Get_Bip(tree->a_nodes[0],tree->a_nodes[0]->v[0],tree);

  n_site = 0;
  For(j,tree->data->crunch_len) For(k,tree->data->wght[j])
    {
      site_num[n_site] = j;
      n_site++;
    }

  boot_data = Copy_Cseq(tree->data,tree->io);

  if (Global_myRank == 0)
    PhyML_Printf("\n. Non parametric bootstrap analysis \n");
  
  //number of bootstraps for each process
  if (tree->mod->bootstrap%Global_numTask != 0) 
    {
      nbRep = (tree->mod->bootstrap / Global_numTask) + 1;
      tree->mod->bootstrap = nbRep * Global_numTask;
      if (Global_myRank == 0) {
        PhyML_Printf("\n. The number of replicates is not a multiple of %d CPUs.\n", Global_numTask);
        PhyML_Printf("\n. Will run %d replicates analysis.\n", tree->mod->bootstrap);
      }
    }
  else
    nbRep = tree->mod->bootstrap/Global_numTask;
  
  //Bip score
  if (Global_myRank == 0) 
    {
      score_tot = (int *)mCalloc((2*tree->n_otu - 3),sizeof(int));
      For(i,2*tree->n_otu-3)
        score_tot[i] = 0;
    }
  else
    score_tot = NULL;

  score_par = (int *)mCalloc((2*tree->n_otu - 3),sizeof(int));
  For(i,2*tree->n_otu-3)
    score_par[i] = 0;

  if (Global_myRank == 0)
    PhyML_Printf("\n  [");

  For(replicate, nbRep)
    {
      For(j,boot_data->crunch_len) boot_data->wght[j] = 0;

      // Send random data to other process
      if (Global_myRank == 0) {
        // Compute number of data to send
        if (tree->mod->bootstrap - randomRecv > Global_numTask)
          nbElem = Global_numTask;
        else
          nbElem = tree->mod->bootstrap - randomRecv;

        For(i,nbElem) {
          For(j,boot_data->crunch_len) boot_data->wght[j] = 0;
          init_len = 0;
          // Create random data
          For(j,boot_data->init_len)
            {
              position = Rand_Int(0,(int)(tree->data->init_len-1.0));
              boot_data->wght[site_num[position]] += 1;
              init_len++;
            }
            

          if (init_len != tree->data->init_len) {
            MPI_Finalize();
            Warn_And_Exit("\n== Pb. when copying sequences...\n");
          }
          // Send random data to other process, not to current process
          if (i < nbElem-1) {
            MPI_Ssend (boot_data->wght, boot_data->crunch_len, MPI_INT, i+1, Global_myRank, MPI_COMM_WORLD);
#ifdef MPI_DEBUG
fprintf (stderr, "\ntask %d, sending random to %d done\n", Global_myRank, i+1);
fflush(stderr);
#endif
          }
          randomRecv++;
        }
      }
      else {
        MPI_Recv (boot_data->wght, boot_data->crunch_len, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &Stat);
#ifdef MPI_DEBUG
fprintf (stderr, "\ntask %d, receiving random from task %d done\n", Global_myRank, Stat.MPI_SOURCE);
fflush(stderr);
#endif
      }

      init_len = 0;
      For(j,boot_data->crunch_len) init_len += boot_data->wght[j];

      if(init_len != tree->data->init_len) {
        MPI_Finalize();
        Warn_And_Exit("\n== Pb when copying sequences\n");
      }

      (tree->mod->io->datatype == NT)?
        (Get_Base_Freqs(boot_data)):
        (Get_AA_Freqs(boot_data));

      if(tree->io->random_boot_seq_order) Randomize_Sequence_Order(boot_data);

      Set_D_States(boot_data,tree->io->datatype,tree->io->state_len);

      boot_mod        = Copy_Model(tree->mod);
      boot_mod->s_opt = tree->mod->s_opt; /* WARNING: re-using the same address here instead of creating a copying
					     requires to leave the value of s_opt unchanged during the boostrap. */
      boot_mod->io    = tree->io; /* WARNING: re-using the same address here instead of creating a copying
				     requires to leave the value of io unchanged during the boostrap. */
      Init_Model(boot_data,boot_mod,tree->io);

      if(tree->io->mod->use_m4mod) M4_Init_Model(boot_mod->m4mod,boot_data,boot_mod);

      if(tree->io->in_tree == 2)
        {
	  switch(tree->io->tree_file_format)
	    {
	    case PHYLIP: 
	      {
		rewind(tree->io->fp_in_tree);
		boot_tree = Read_Tree_File_Phylip(tree->io->fp_in_tree);
		break;
	      }
	    case NEXUS:
	      {
		PhyML_Printf("\n. Unfortunately, PhyML cannot read NEXUS files and perform a bootstrap analysis."); 
		PhyML_Printf("\n. Please use the PHYLIP format.."); 
		PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		Warn_And_Exit("");
		break;
	      }
	    default:
	      {
		PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		Warn_And_Exit("");
		break;
	      }
	    }
	}
      else
        {
          boot_mat = ML_Dist(boot_data,boot_mod);
          boot_mat->tree = Make_Tree_From_Scratch(boot_data->n_otu,boot_data);
          Fill_Missing_Dist(boot_mat);
          Bionj(boot_mat);
          boot_tree = boot_mat->tree;
          boot_tree->mat = boot_mat;
        }

      boot_tree->mod                = boot_mod;
      boot_tree->io                 = tree->io;
      boot_tree->data               = boot_data;
      boot_tree->mod->s_opt->print  = 0;
      boot_tree->n_pattern          = boot_tree->data->crunch_len;
      boot_tree->io->print_site_lnl = 0;
      boot_tree->io->print_trace    = 0;

      Set_Both_Sides(YES,boot_tree);
      
      if((boot_tree->mod->s_opt->random_input_tree) && (boot_tree->mod->s_opt->topo_search == SPR_MOVE)) Random_Tree(boot_tree);

      Connect_CSeqs_To_Nodes(boot_data,boot_tree);
      Share_Lk_Struct(tree,boot_tree);
      Share_Spr_Struct(tree,boot_tree);
      Share_Pars_Struct(tree,boot_tree);
      Fill_Dir_Table(boot_tree);
      Update_Dirs(boot_tree);

      if(tree->mod->s_opt->greedy) Init_P_Lk_Tips_Double(boot_tree);
      else                         Init_P_Lk_Tips_Int(boot_tree);
      Init_Ui_Tips(boot_tree);
      Init_P_Pars_Tips(boot_tree);
      Br_Len_Not_Involving_Invar(boot_tree);
      
      if(boot_tree->mod->s_opt->opt_topo)
        {
          if(boot_tree->mod->s_opt->topo_search == NNI_MOVE) 
            {
              Simu_Loop(boot_tree);
            }
          else if((boot_tree->mod->s_opt->topo_search == SPR_MOVE) ||
                  (boot_tree->mod->s_opt->topo_search == BEST_OF_NNI_AND_SPR))
            {
              Speed_Spr_Loop(boot_tree);
            }
        }
      else
        {
          if(boot_tree->mod->s_opt->opt_subst_param || boot_tree->mod->s_opt->opt_bl)
            Round_Optimize(boot_tree,boot_tree->data,ROUND_MAX);
          else
            Lk(NULL,boot_tree);
        }

      Free_Bip(boot_tree);

      Alloc_Bip(boot_tree);

      Match_Tip_Numbers(tree,boot_tree);

      Get_Bip(boot_tree->a_nodes[0],
              boot_tree->a_nodes[0]->v[0],
              boot_tree);

      Compare_Bip(tree,boot_tree,NO);
      
      Br_Len_Involving_Invar(boot_tree);

      if(tree->io->print_boot_trees)
        {
          s = Write_Tree(boot_tree,NO);
          t=(char *)mCalloc(T_MAX_LINE,sizeof(char));
          Print_Fp_Out_Lines_MPI(boot_tree, tree->io, replicate+1, t);
          
          // Get bootstrap trees from other process and write to boot file
          if (Global_myRank == 0) {
            fprintf(tree->io->fp_out_boot_tree,"%s\n",s);
            fprintf(tree->io->fp_out_boot_stats,"%s\n",t);
            bootRecv++;
            PhyML_Printf(".");
	    if(!((bootRecv)%tree->io->boot_prog_every)) 
	      {
		PhyML_Printf("] %4d/%4d\n  ",bootRecv,tree->mod->bootstrap);
		if(bootRecv != tree->mod->bootstrap) PhyML_Printf("[");
	      }
	    
            // Compute number of bootstraps to receive
            if (tree->mod->bootstrap - bootRecv > Global_numTask)
              nbElem = Global_numTask;
            else
              nbElem = tree->mod->bootstrap - bootRecv + 1;
              
            bootStr=(char *)mCalloc(T_MAX_LINE,sizeof(char));
            for (i=1; i<nbElem; i++) 
	      {
		MPI_Recv (bootStr, T_MAX_LINE, MPI_CHAR, i, MPI_ANY_TAG, MPI_COMM_WORLD, &Stat);
#ifdef MPI_DEBUG
		PhyML_Fprintf (stderr, "\ntask %d, receiving bootstrap from task %d tag %d done\n", Global_myRank, Stat.MPI_SOURCE, Stat.MPI_TAG);
		fflush(stderr);
#endif
		if (Stat.MPI_TAG == BootTreeTag)
		  fprintf(tree->io->fp_out_boot_tree,"%s\n", bootStr);
		if (Stat.MPI_TAG == BootStatTag)
		  fprintf(tree->io->fp_out_boot_stats,"%s\n", bootStr);
		
		MPI_Recv (bootStr, T_MAX_LINE, MPI_CHAR, i, MPI_ANY_TAG, MPI_COMM_WORLD, &Stat);
#ifdef MPI_DEBUG
		fprintf (stderr, "\ntask %d, receiving bootstrap from task %d tag %d done\n", Global_myRank, Stat.MPI_SOURCE, Stat.MPI_TAG);
		fflush(stderr);
#endif
		if (Stat.MPI_TAG == BootTreeTag)
		  fprintf(tree->io->fp_out_boot_tree,"%s\n", bootStr);
		if (Stat.MPI_TAG == BootStatTag)
		  fprintf(tree->io->fp_out_boot_stats,"%s\n", bootStr);
		
		bootRecv++;
		PhyML_Printf(".");
		if(!((bootRecv)%tree->io->boot_prog_every)) {
		  PhyML_Printf("] %4d/%4d\n  ",bootRecv,tree->mod->bootstrap);
		  if(bootRecv != tree->mod->bootstrap) PhyML_Printf("[");
		}
	      }
            Free(bootStr);
          }
          else {
	    MPI_Ssend (s, T_MAX_LINE, MPI_CHAR, 0, BootTreeTag, MPI_COMM_WORLD);
	    MPI_Ssend (t, T_MAX_LINE, MPI_CHAR, 0, BootStatTag, MPI_COMM_WORLD);
#ifdef MPI_DEBUG
fprintf (stderr, "\ntask %d, sending bootstraps done\n", Global_myRank);
fflush(stderr);
#endif
          }
          Free(t);
          Free(s);
        }

#ifndef QUIET
fflush(stdout);
#endif

      /* if(boot_tree->mat) Free_Mat(boot_tree->mat); */
      Free_Tree(boot_tree);
      Free_Model(boot_mod);

      //Each process computes the Bip score sum for all its bootstrap trees
      For(i,2*tree->n_otu-3) score_par[i] = tree->a_edges[i]->bip_score;

      //Each process sends its Bip score sum. The sums are summed.
      MPI_Reduce(score_par, score_tot, 2*tree->n_otu - 3, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    }

  
  if (Global_myRank == 0) 
    {
      For(i,2*tree->n_otu-3)
        tree->a_edges[i]->bip_score = score_tot[i];
      Free (score_tot);
    }
  Free (score_par);

  if (Global_myRank == 0)
    if(((bootRecv)%tree->io->boot_prog_every)) PhyML_Printf("] %4d/%4d\n ",bootRecv,tree->mod->bootstrap);

PhyML_Printf("\n\n. Exiting bootstrap function normally."); fflush(NULL);
  tree->lock_topo = 1; /* TopoLOGy should not be modified afterwards */

  if(tree->io->print_boot_trees)
    {
      fclose(tree->io->fp_out_boot_tree);
      fclose(tree->io->fp_out_boot_stats);
    }

  Free_Cseq(boot_data);
  Free(site_num);
}

/*********************************************************/

void Print_Fp_Out_Lines_MPI(t_tree *tree, option *io, int n_data_set, char *bootStr)
{
  char *s, *tmp;

  // Build a string to be sent to the writing process
  s = (char *)mCalloc(T_MAX_LINE,sizeof(char));
  tmp=(char *)mCalloc(T_MAX_LINE,sizeof(char));
  
  if (Global_myRank == 0 && n_data_set == 1) {
    snprintf(tmp, T_MAX_LINE, ". Sequence file : [%s]\n\n", Basename(io->in_align_file)); strncat (s, tmp, T_MAX_LINE);
    
    (tree->mod->io->datatype == NT)?
      (snprintf(tmp, T_MAX_LINE, ". Model of nucleotides substitution : %s\n\n", io->mod->modelname->s)):
      (snprintf(tmp, T_MAX_LINE, ". Model of amino acids substitution : %s\n\n", io->mod->modelname->s));
    strncat (s, tmp, T_MAX_LINE);
    
    switch(io->in_tree)
      {
      case 0: { snprintf(tmp, T_MAX_LINE, ". Initial tree : [BioNJ]\n\n");               break; }
      case 1: { snprintf(tmp, T_MAX_LINE, ". Initial tree : [parsimony]\n\n");           break; }
      case 2: { snprintf(tmp, T_MAX_LINE, ". Initial tree : [%s]\n\n",io->in_tree_file); break; }
      }
    strncat (s, tmp, T_MAX_LINE);

    strncat (s, "\n", T_MAX_LINE);
    
    /*headline 1*/
    strncat (s, ". Data\t", T_MAX_LINE);

    strncat (s, "Nb of \t", T_MAX_LINE);

    strncat (s, "Likelihood\t", T_MAX_LINE);

    strncat (s, "Discrete   \t", T_MAX_LINE);

    if(tree->mod->ras->n_catg > 1)
      strncat (s, "Number of \tGamma shape\t", T_MAX_LINE);

    strncat (s, "Proportion of\t", T_MAX_LINE);

    if(tree->mod->whichmodel <= 6)
      strncat (s, "Transition/ \t", T_MAX_LINE);
    
    strncat (s, "Nucleotides frequencies               \t", T_MAX_LINE);

    if((tree->mod->whichmodel == GTR) ||
       (tree->mod->whichmodel == CUSTOM))
      strncat (s, "Instantaneous rate matrix              \t", T_MAX_LINE);

    strncat (s, "\n", T_MAX_LINE);

    /*headline 2*/
    strncat (s, "  set\t", T_MAX_LINE);

    strncat (s, "taxa\t", T_MAX_LINE);

    strncat (s, "LOGlk     \t", T_MAX_LINE);

    strncat (s, "gamma model\t", T_MAX_LINE);

    if(tree->mod->ras->n_catg > 1)
      strncat (s, "categories\tparameter  \t", T_MAX_LINE);

    strncat (s, "invariant    \t", T_MAX_LINE);

    if(tree->mod->whichmodel <= 6)
      strncat (s, "transversion\t", T_MAX_LINE);

    strncat (s, "f(A)      f(C)      f(G)      f(T)    \t", T_MAX_LINE);

    if((tree->mod->whichmodel == GTR) ||
       (tree->mod->whichmodel == CUSTOM))
      strncat (s, "[A---------C---------G---------T------]\t", T_MAX_LINE);

    strncat (s, "\n", T_MAX_LINE);

    /*headline 3*/
    if(tree->mod->whichmodel == TN93) {
      strncat (s, "    \t      \t          \t           \t", T_MAX_LINE);
      if(tree->mod->ras->n_catg > 1)
        strncat (s, "         \t         \t", T_MAX_LINE);
      strncat (s, "             \t", T_MAX_LINE);
      strncat (s, "purines pyrimid.\t", T_MAX_LINE);
      strncat (s, "\n", T_MAX_LINE);
    }

    strncat (s, "\n", T_MAX_LINE);
  }

  /*line items*/

  snprintf(tmp, T_MAX_LINE, "  #%d\t", (((n_data_set-1)*Global_numTask)+Global_myRank+1));
  strncat (s, tmp, T_MAX_LINE);
  
  snprintf(tmp, T_MAX_LINE, "%d   \t",tree->n_otu); strncat (s, tmp, T_MAX_LINE);
  
  snprintf(tmp, T_MAX_LINE, "%.5f\t",tree->c_lnL); strncat (s, tmp, T_MAX_LINE);
  
  snprintf(tmp, T_MAX_LINE, "%s        \t",
          (tree->mod->ras->n_catg>1)?("Yes"):("No ")); strncat (s, tmp, T_MAX_LINE);
  
  if(tree->mod->ras->n_catg > 1)
    {
      snprintf(tmp, T_MAX_LINE, "%d        \t",tree->mod->ras->n_catg); strncat (s, tmp, T_MAX_LINE);
      snprintf(tmp, T_MAX_LINE, "%.3f    \t",tree->mod->ras->alpha->v); strncat (s, tmp, T_MAX_LINE);
    }
  
  snprintf(tmp, T_MAX_LINE, "%.3f    \t",tree->mod->ras->pinvar->v); strncat (s, tmp, T_MAX_LINE);
  
  if(tree->mod->whichmodel <= 5)
    {
      snprintf(tmp, T_MAX_LINE, "%.3f     \t",tree->mod->kappa->v); strncat (s, tmp, T_MAX_LINE);
    }
  else if(tree->mod->whichmodel == TN93)
    {
      snprintf(tmp, T_MAX_LINE, "%.3f   ",
              tree->mod->kappa->v*2.*tree->mod->lambda->v/(1.+tree->mod->lambda->v));
      strncat (s, tmp, T_MAX_LINE);
      snprintf(tmp, T_MAX_LINE, "%.3f\t",
              tree->mod->kappa->v*2./(1.+tree->mod->lambda->v));
      strncat (s, tmp, T_MAX_LINE);
    }

  if(tree->mod->io->datatype == NT)
    {
      snprintf(tmp, T_MAX_LINE, "%8.5f  ",tree->mod->e_frq->pi->v[0]); strncat (s, tmp, T_MAX_LINE);
      snprintf(tmp, T_MAX_LINE, "%8.5f  ",tree->mod->e_frq->pi->v[1]); strncat (s, tmp, T_MAX_LINE);
      snprintf(tmp, T_MAX_LINE, "%8.5f  ",tree->mod->e_frq->pi->v[2]); strncat (s, tmp, T_MAX_LINE);
      snprintf(tmp, T_MAX_LINE, "%8.5f\t",tree->mod->e_frq->pi->v[3]); strncat (s, tmp, T_MAX_LINE);
    }
    
  if((tree->mod->whichmodel == GTR) || (tree->mod->whichmodel == CUSTOM))
    {
      int i,j;

      For(i,4)
        {
          if (i!=0) {
            /*format*/
            snprintf(tmp, T_MAX_LINE, "      \t     \t          \t           \t");
            strncat (s, tmp, T_MAX_LINE);
            if(tree->mod->ras->n_catg > 1) {
              snprintf(tmp, T_MAX_LINE, "          \t           \t");
              strncat (s, tmp, T_MAX_LINE);
            }
            snprintf(tmp, T_MAX_LINE, "             \t                                      \t");
            strncat (s, tmp, T_MAX_LINE);
          }
          For(j,4) {
            snprintf(tmp, T_MAX_LINE, "%8.5f  ",tree->mod->r_mat->qmat->v[i*4+j]);
            strncat (s, tmp, T_MAX_LINE);
          }
          if (i<3) {
            snprintf(tmp, T_MAX_LINE, "\n");
            strncat (s, tmp, T_MAX_LINE);
          }
        }
    }
    
    Free (tmp);
    strncpy (bootStr, s, T_MAX_LINE);
    Free (s);

  return;
}
/* #endif */
