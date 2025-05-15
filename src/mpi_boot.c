/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include "mpi_boot.h"

int Global_numTask;
int Global_myRank;


/* #ifdef MPI */

/*********************************************************/

/* 
   if tbe_bootstrap == 0 => Classical FBP (Felsenstein bootstrap proportions) 
   else => TBE (Transfer bootstrap expectation)
*/
void Bootstrap_MPI(t_tree *tree)
{
	int *site_num, n_site;
	int replicate, j, k;
	int position, init_len;
	calign *boot_data;
	t_tree *boot_tree;
	t_mod *boot_mod;
	matrix *boot_mat;
	char *s;
	MPI_Status Stat;
	int iterations_per_process, bootRecv, nbElem, i;
	int *score_par, *score_tot;
	char *bootStr, *t;
	
	if(tree->is_mixt_tree == YES)
	{
		PhyML_Printf("\n== Bootstrap option not yet available for partition/mixture analysis.");
		assert(FALSE);
	}
	
	nbElem = bootRecv = 0;
	tree->io->print_support_val = YES;
	boot_tree = NULL;
	site_num = (int *)mCalloc(tree->data->init_len,sizeof(int));
	
	Free_Bip(tree);
	Alloc_Bip(tree);
	Get_Bip(tree->a_nodes[0],tree->a_nodes[0]->v[0],tree);
	
	n_site = 0;
	for(j=0;j<tree->data->n_pattern;j++)
	{
		for(k=0;k<tree->data->wght[j];++k)
		{
			site_num[n_site] = j;
			n_site++;
		}
	}
	
	boot_data = Copy_Cseq(tree->data,tree->io, NULL);
	
	if(Global_numTask <= 1)
	{
		PhyML_Printf("\n\n. The number of CPUs used in the MPI bootstrap analysis should be");
		PhyML_Printf("\n. strictly greater than 1 (it is equal to %d here).",Global_numTask);
		assert(FALSE);
	}
	if (Global_myRank == 0)
	{
		PhyML_Printf("\n\n. Non parametric bootstrap analysis \n\n");
		PhyML_Printf("\n  [");
	}
	
	iterations_per_process = tree->io->n_boot_replicates / Global_numTask;
	// When number of bootstrap is not multiple of the number of tasks, provide one more iteration
	if (tree->io->n_boot_replicates % Global_numTask != 0) iterations_per_process++;
	
	// Bip score
	if (Global_myRank == 0) 
	{
		score_tot = (int *)mCalloc((2*tree->n_otu - 3),sizeof(int));
		for(i=0;i<2*tree->n_otu-3;++i) score_tot[i] = 0;
	}
	else
		score_tot = NULL;
	
	score_par = (int *)mCalloc((2*tree->n_otu - 3),sizeof(int));
	for(i=0;i<2*tree->n_otu-3;++i) score_par[i] = 0;
	
	for (replicate = 0; replicate < iterations_per_process; replicate++)
	{
		// Avoid computing too many replicates
		if ( (replicate*Global_numTask) + Global_myRank < tree->io->n_boot_replicates)
		{
#ifdef MPI_DEBUG
			fprintf (stderr, "\nThread %d, computing bootstrap replicate %d\n", Global_myRank, (replicate*Global_numTask) + Global_myRank + 1);
			fflush(stderr);
#endif
			// Compute one bootstrap replicate	
			for(j=0;j<boot_data->n_pattern;j++) boot_data->wght[j] = 0;

			init_len = 0;
			for(j=0;j<boot_data->init_len;j++)
			{
				position = Rand_Int(0,(int)(tree->data->init_len-1.0));
				boot_data->wght[site_num[position]] += 1;
				init_len++;
			}
			if (init_len != tree->data->init_len) 
			{
				PhyML_Printf("\n== thread: %d || init: %d %d here. Pb when copying sequences\n",Global_myRank,init_len,tree->data->init_len);
				fflush(NULL);
				Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
			}

			init_len = 0;
			for(j=0;j<boot_data->n_pattern;j++)
				init_len += boot_data->wght[j];
			
			if (init_len != tree->data->init_len)
			{
				printf("\n== thread: %d init: %d %d. Pb when copying sequences\n",Global_myRank,init_len,tree->data->init_len);
				Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
			}
			
			if(tree->io->datatype == NT)      Get_Base_Freqs(boot_data);
			else if(tree->io->datatype == AA) Get_AA_Freqs(boot_data);
			
			if(tree->io->random_boot_seq_order) Randomize_Sequence_Order(boot_data);
			
			Set_D_States(boot_data,tree->io->datatype,tree->io->state_len);
			
			boot_mod        = Copy_Model(tree->mod);
			boot_mod->s_opt = tree->mod->s_opt; /* WARNING: re-using the same address here instead of creating a copying
												requires to leave the value of s_opt unchanged during the boostrap. */
			boot_mod->io    = tree->io; /* WARNING: re-using the same address here instead of creating a copying
										requires to leave the value of io unchanged during the boostrap. */
			
			Init_Model(boot_data,boot_mod,tree->io);
			Set_Model_Parameters(boot_mod);
			
			if(tree->io->in_tree == 2)
			{
				rewind(tree->io->fp_in_tree);
				boot_tree = Read_Tree_File_Phylip(tree->io->fp_in_tree);
				Remove_Duplicates_From_Tree(boot_data,boot_tree);
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
			
			boot_tree->mod                  = boot_mod;
			boot_tree->io                   = tree->io;
			boot_tree->data                 = boot_data;
			boot_tree->verbose              = VL0;
			boot_tree->io->print_site_lnl   = NO;
			boot_tree->io->print_trace      = NO;
			boot_tree->io->print_json_trace = NO;
			boot_tree->n_root               = NULL;
			boot_tree->e_root               = NULL;
			boot_tree->l_ev                 = tree->l_ev;
			boot_tree->p_lk_left_pi         = tree->p_lk_left_pi;
			
#if (defined(__AVX__) || defined(__AVX2__) || defined(__SSE__) || defined(__SSE2__) || defined(__SSE3__))
			boot_tree->_tPij1               = tree->_tPij1;
			boot_tree->_tPij2               = tree->_tPij2;
			boot_tree->_pmat1plk1           = tree->_pmat1plk1;
			boot_tree->_pmat2plk2           = tree->_pmat2plk2;
			boot_tree->_plk0                = tree->_plk0;
			boot_tree->_l_ev                = tree->_l_ev;
			boot_tree->_r_ev                = tree->_r_ev;
			boot_tree->_prod_left           = tree->_prod_left;
			boot_tree->_prod_rght           = tree->_prod_rght;
#endif

			Set_Both_Sides(YES,boot_tree);
			
			if((boot_tree->mod->s_opt->random_input_tree) && (boot_tree->mod->s_opt->topo_search == SPR_MOVE)) Random_Tree(boot_tree);
			
			Connect_CSeqs_To_Nodes(boot_data,boot_tree->io,boot_tree);
			Share_Lk_Struct(tree,boot_tree);
			Share_Spr_Struct(tree,boot_tree);
			Share_Pars_Struct(tree,boot_tree);
			Fill_Dir_Table(boot_tree);
			Update_Dirs(boot_tree);
			
			Init_Partial_Lk_Tips_Double(boot_tree);
			Init_Ui_Tips(boot_tree);
			Init_Partial_Pars_Tips(boot_tree);
			Br_Len_Not_Involving_Invar(boot_tree);
			
			Set_Update_Eigen(YES,boot_tree->mod);
			Lk(NULL,boot_tree);
			Set_Update_Eigen(NO,boot_tree->mod);
			
			if (boot_tree->mod->s_opt->opt_topo)
			{
				Global_Spr_Search(boot_tree);
			}
			else
			{
				if (boot_tree->mod->s_opt->opt_subst_param || boot_tree->a_edges[0]->l->optimize == YES)
					Round_Optimize(boot_tree,ROUND_MAX);
				else
					Lk(NULL,boot_tree);
			}
			
			Free_Bip(boot_tree);
			Alloc_Bip(boot_tree);
			Match_Tip_Numbers(tree,boot_tree);
			Get_Bip(boot_tree->a_nodes[0],
				boot_tree->a_nodes[0]->v[0],
				boot_tree);
			
			if(tree->io->do_boot)     Compare_Bip(tree,boot_tree,NO, TREE_COMP_RF_PLAIN);
			else if(tree->io->do_tbe) Compare_Bip_Distance(tree, boot_tree);
			else assert(FALSE);
			
			Check_Br_Lens(boot_tree);
			Br_Len_Involving_Invar(boot_tree);
			
			if(tree->io->print_boot_trees)
			{
				s = Write_Tree(boot_tree);
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
						PhyML_Printf("] %4d/%4d\n  ",bootRecv,tree->io->n_boot_replicates);
						if (bootRecv != tree->io->n_boot_replicates) PhyML_Printf("[");
					}
					// Compute number of bootstraps to receive
					nbElem = Global_numTask -1;
					if (bootRecv + nbElem > tree->io->n_boot_replicates) nbElem = tree->io->n_boot_replicates - bootRecv;
					
					bootStr=(char *)mCalloc(T_MAX_LINE,sizeof(char));
					i=0;
					while (i<nbElem)
					{
						i++;
						MPI_Recv (bootStr, T_MAX_LINE, MPI_CHAR, i, MPI_ANY_TAG, MPI_COMM_WORLD, &Stat);
#ifdef MPI_DEBUG
						fprintf (stderr, "\nThread %d, bootstrap data received from thread %d tag %d\n", Global_myRank, Stat.MPI_SOURCE, Stat.MPI_TAG);
						fflush(stderr);
#endif
						if (Stat.MPI_TAG == BootTreeTag)
							fprintf(tree->io->fp_out_boot_tree,"%s\n", bootStr);
						if (Stat.MPI_TAG == BootStatTag)
							fprintf(tree->io->fp_out_boot_stats,"%s\n", bootStr);
						
						MPI_Recv (bootStr, T_MAX_LINE, MPI_CHAR, i, MPI_ANY_TAG, MPI_COMM_WORLD, &Stat);
#ifdef MPI_DEBUG
						fprintf (stderr, "\nThread %d, bootstrap data received from thread %d tag %d\n", Global_myRank, Stat.MPI_SOURCE, Stat.MPI_TAG);
						fflush(stderr);
#endif
						if (Stat.MPI_TAG == BootTreeTag)
							fprintf(tree->io->fp_out_boot_tree,"%s\n", bootStr);
						if (Stat.MPI_TAG == BootStatTag)
							fprintf(tree->io->fp_out_boot_stats,"%s\n", bootStr);
						
						bootRecv++;
						PhyML_Printf(".");
						if (!((bootRecv)%tree->io->boot_prog_every)) {
							PhyML_Printf("] %4d/%4d\n  ",bootRecv,tree->io->n_boot_replicates);
							if (bootRecv != tree->io->n_boot_replicates) PhyML_Printf("[");
						}
					}
					Free (bootStr);
				}
				else
				{
					MPI_Ssend (s, T_MAX_LINE, MPI_CHAR, 0, BootTreeTag, MPI_COMM_WORLD);
					MPI_Ssend (t, T_MAX_LINE, MPI_CHAR, 0, BootStatTag, MPI_COMM_WORLD);
#ifdef MPI_DEBUG
					fprintf (stderr, "\nThread %d, bootstrap data sent\n", Global_myRank);
					fflush(stderr);
#endif
				}
				Free(t);
				Free(s);
			}
			
#ifndef QUIET
fflush(stdout);
#endif

			Free_Tree(boot_tree);
			Free_Model(boot_mod);

			//Each process computes the Bip score sum for all its bootstrap trees
			for(i=0;i<2*tree->n_otu-3;++i) score_par[i] = tree->a_edges[i]->bip_score;

			//Each process sends its Bip score sum. The sums are summed.
			MPI_Reduce(score_par, score_tot, 2*tree->n_otu - 3, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		}
		else
		{
			// For the remaining process that did not have a last replicate to compute
			for(i=0;i<2*tree->n_otu-3;++i) score_par[i] = 0;
			
			MPI_Reduce(score_par, score_tot, 2*tree->n_otu - 3, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		}
	}
	
	if (Global_myRank == 0)
	{
		for(i=0;i<2*tree->n_otu-3;++i) tree->a_edges[i]->bip_score = score_tot[i];
		Free (score_tot);
	}
	Free (score_par);
	
	if(Global_myRank == 0)
		if(((bootRecv)%tree->io->boot_prog_every))
			PhyML_Printf("] %4d/%4d\n ",bootRecv,tree->io->n_boot_replicates);
	
	PhyML_Printf("\n\n. Exiting bootstrap function normally."); fflush(NULL);
	tree->lock_topo = 1; /* Topology should not be modified afterwards */
	
	if(tree->io->print_boot_trees)
	{
		fclose(tree->io->fp_out_boot_tree);
		fclose(tree->io->fp_out_boot_stats);
	}
	
	Free_Calign(boot_data);
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

    strncat (s, "loglk     \t", T_MAX_LINE);

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

      for(i=0;i<4;i++)
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
          for(j=0;j<4;j++) {
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
