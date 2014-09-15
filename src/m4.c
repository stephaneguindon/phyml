/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homoLOGous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/


/* Routines for Markov-Modulated Markov Models (M4) */

#include "m4.h"

int M4_main(int argc, char **argv)
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

  
#ifdef MPI
  int rc;
  rc = MPI_Init(&argc,&argv);
  if (rc != MPI_SUCCESS) {
    PhyML_Printf("\n== Err. starting MPI program. Terminating.\n");
    MPI_Abort(MPI_COMM_WORLD, rc);
  }
  MPI_Comm_size(MPI_COMM_WORLD,&Global_numTask);
  MPI_Comm_rank(MPI_COMM_WORLD,&Global_myRank);
#endif

#ifdef QUIET
  setvbuf(stdout,NULL,_IOFBF,2048);
#endif


  tree             = NULL;
  mod              = NULL;
  best_lnL         = UNLIKELY;
      
  io = (option *)Get_Input(argc,argv);
  r_seed = (io->r_seed < 0)?(time(NULL)):(io->r_seed);
  srand(r_seed);
  io->r_seed = r_seed;


  if(io->in_tree == 2) Test_Multiple_Data_Set_Format(io);
  else io->n_trees = 1;


  if((io->n_data_sets > 1) && (io->n_trees > 1))
    {
      io->n_data_sets = MIN(io->n_trees,io->n_data_sets);
      io->n_trees     = MIN(io->n_trees,io->n_data_sets);
    }

  For(num_data_set,io->n_data_sets)
    {
      best_lnL = UNLIKELY;
      Get_Seq(io);
      Make_Model_Complete(io->mod);
      Set_Model_Name(io->mod);
      Print_Settings(io);
      mod = io->mod;
      
      if(io->data)
	{
	  if(io->n_data_sets > 1) PhyML_Printf("\n. Data set [#%d]\n",num_data_set+1);
	  cdata = Compact_Data(io->data,io);

	  Free_Seq(io->data,cdata->n_otu);
	  
	  if(cdata) Check_Ambiguities(cdata,io->datatype,io->state_len);
	  else
	    {
	      PhyML_Printf("\n== Err. in file %s at line %d\n",__FILE__,__LINE__);
	      Warn_And_Exit("");
	    }

	  for(num_tree=(io->n_trees == 1)?(0):(num_data_set);num_tree < io->n_trees;num_tree++)
	    {
	      if(!io->mod->s_opt->random_input_tree) io->mod->s_opt->n_rand_starts = 1;

	      For(num_rand_tree,io->mod->s_opt->n_rand_starts)
		{
		  if((io->mod->s_opt->random_input_tree) && (io->mod->s_opt->topo_search != NNI_MOVE))
		    if(!io->quiet) PhyML_Printf("\n. [Random start %3d/%3d]\n",num_rand_tree+1,io->mod->s_opt->n_rand_starts);

		  Init_Model(cdata,mod,io);

		  if(io->mod->use_m4mod) M4_Init_Model(mod->m4mod,cdata,mod);

		  switch(io->in_tree)
		    {
		    case 0 : case 1 : { tree = Dist_And_BioNJ(cdata,mod,io); break; }
		    case 2 :          { tree = Read_User_Tree(cdata,mod,io); break; }
		    }


		  if(!tree) continue;

		  time(&t_beg);
		  time(&(tree->t_beg));
      

		  tree->mod         = mod;
		  tree->io          = io;
		  tree->data        = cdata;
		  tree->n_pattern   = tree->data->crunch_len;

                  Set_Both_Sides(YES,tree);
    
		  if(mod->s_opt->random_input_tree) Random_Tree(tree);

		  if((!num_data_set) && (!num_tree) && (!num_rand_tree)) Check_Memory_Amount(tree);

		  Prepare_Tree_For_Lk(tree);

		  if(io->in_tree == 1) Spr_Pars(tree);
		 
		  if(io->do_alias_subpatt)
		    {
		      MIXT_Set_Alias_Subpatt(YES,tree);
		      Lk(NULL,tree);
		      MIXT_Set_Alias_Subpatt(NO,tree);
		    }

		  if(tree->mod->s_opt->opt_topo)
		    {
		      if(tree->mod->s_opt->topo_search      == NNI_MOVE) Simu_Loop(tree);
		      else if(tree->mod->s_opt->topo_search == SPR_MOVE) Speed_Spr_Loop(tree);
		      else                                               Best_Of_NNI_And_SPR(tree);
		    }
		  else
		    {
		      if(tree->mod->s_opt->opt_subst_param || 
			 tree->mod->s_opt->opt_bl)                       Round_Optimize(tree,tree->data,ROUND_MAX);
		      else                                               Lk(NULL,tree);
		    }
		  
		  
                  Set_Both_Sides(YES,tree);     
		  Lk(NULL,tree);
		  Pars(NULL,tree);
		  Get_Tree_Size(tree);
		  PhyML_Printf("\n. Log likelihood of the current tree: %f.\n",tree->c_lnL);

		  Exit("\n");

		  /* */
		  M4_Compute_Proba_Hidden_States_On_Edges(tree);
		  /* */

		  Get_Best_Root_Position(tree);

		  /* Print the tree estimated using the current random (or BioNJ) starting tree */
		  if(io->mod->s_opt->n_rand_starts > 1)
		    {
		      Br_Len_Involving_Invar(tree);
		      Print_Tree(io->fp_out_trees,tree);
		      fflush(NULL);
		    }

		  /* Record the most likely tree in a string of characters */
		  if(tree->c_lnL > best_lnL)
		    {
		      best_lnL = tree->c_lnL;
		      Br_Len_Involving_Invar(tree);
		      if(most_likely_tree) Free(most_likely_tree);
		      most_likely_tree = Write_Tree(tree,NO);
		      Get_Tree_Size(tree);
		    }

/* 		  JF(tree); */

		  time(&t_end);
		  Print_Fp_Out(io->fp_out_stats,t_beg,t_end,tree,
			       io,num_data_set+1,
			       (tree->mod->s_opt->n_rand_starts > 1)?
			       (num_rand_tree):(num_tree),YES);
		  
		  if(tree->io->print_site_lnl) Print_Site_Lk(tree,io->fp_out_lk);

		  /* Start from BioNJ tree */
		  if((num_rand_tree == io->mod->s_opt->n_rand_starts-1) && (tree->mod->s_opt->random_input_tree))
		    {
		      /* Do one more iteration in the loop, but don't randomize the tree */
		      num_rand_tree--;
		      tree->mod->s_opt->random_input_tree = 0;
		    }
		  
		  Free_Spr_List(tree);
		  Free_One_Spr(tree->best_spr);
		  if(tree->mat) Free_Mat(tree->mat);
		  Free_Triplet(tree->triplet_struct);
		  Free_Tree_Pars(tree);
		  Free_Tree_Lk(tree);
		  Free_Tree(tree);
		}


	      /* Launch bootstrap analysis */
	      if(mod->bootstrap) 
		{
		  if(!io->quiet) PhyML_Printf("\n. Launch bootstrap analysis on the most likely tree...\n");

                  #ifdef MPI
		  MPI_Bcast (most_likely_tree, strlen(most_likely_tree)+1, MPI_CHAR, 0, MPI_COMM_WORLD);
		  if(!io->quiet)  PhyML_Printf("\n. The bootstrap analysis will use %d CPUs.\n",Global_numTask);
		  #endif

		  most_likely_tree = Bootstrap_From_String(most_likely_tree,cdata,mod,io);
		}
	      else if(io->ratio_test) 
		{
		  /* Launch aLRT */
		  if(!io->quiet) PhyML_Printf("\n. Compute aLRT branch supports on the most likely tree...\n");
		  most_likely_tree = aLRT_From_String(most_likely_tree,cdata,mod,io);
		}

	      /* Print the most likely tree in the output file */
	      if(!io->quiet) PhyML_Printf("\n. Printing the most likely tree in file '%s'...\n", Basename(io->out_tree_file));
	      if(io->n_data_sets == 1) rewind(io->fp_out_tree);
	      PhyML_Fprintf(io->fp_out_tree,"%s\n",most_likely_tree);
	      

	      if(io->n_trees > 1 && io->n_data_sets > 1) break;
	    }
	  Free_Cseq(cdata);
	}
      else
	{
	  PhyML_Printf("\n. No data was found.\n");
	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}
      Free_Model_Complete(mod);
    }
  
  if(most_likely_tree) Free(most_likely_tree);

  if(mod->s_opt->n_rand_starts > 1) PhyML_Printf("\n. Best log likelihood: %f\n",best_lnL);

  Free_Optimiz(mod->s_opt);
  Free_Model_Basic(mod);

  if(io->fp_in_align)  fclose(io->fp_in_align);
  if(io->fp_in_tree)   fclose(io->fp_in_tree);
  if(io->fp_out_lk)    fclose(io->fp_out_lk);
  if(io->fp_out_tree)  fclose(io->fp_out_tree);
  if(io->fp_out_trees) fclose(io->fp_out_trees);
  if(io->fp_out_stats) fclose(io->fp_out_stats);

  Free_Input(io);

  time(&t_end);
  Print_Time_Info(t_beg,t_end);

#ifdef MPI
  MPI_Finalize();
#endif

  return 0;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


/* Allocate memory */
m4 *M4_Make_Light()
{
  m4 *m4mod;

  m4mod = (m4 *)mCalloc(1,sizeof(m4));
  M4_Set_M4mod_Default(m4mod);
  return m4mod;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void M4_Set_M4mod_Default(m4 *m4mod)
{
  m4mod->use_cov_alpha = 1;
  m4mod->use_cov_alpha = 0;
  m4mod->n_h           = 3;
  m4mod->n_o           = 4;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* Allocate memory */
void M4_Make_Complete(int n_h, int n_o, m4 *m4mod)
{
  int i;

  m4mod->n_h = n_h;
  m4mod->n_o = n_o;
  m4mod->n_o = n_o;
  m4mod->o_rr = (phydbl *)mCalloc(n_o*n_o,sizeof(phydbl));
  m4mod->o_fq = (phydbl *)mCalloc(n_o,sizeof(phydbl));
  m4mod->o_mats = (phydbl **)mCalloc(n_h,sizeof(phydbl *));
  For(i,n_h) m4mod->o_mats[i] = (phydbl *)mCalloc(n_o*n_o,sizeof(phydbl));
  m4mod->h_mat = (phydbl *)mCalloc(n_h*n_h,sizeof(phydbl));
  m4mod->h_rr = (phydbl *)mCalloc(n_h*n_h,sizeof(phydbl));
  m4mod->h_fq = (phydbl *)mCalloc(n_h,sizeof(phydbl));
  m4mod->multipl = (phydbl *)mCalloc(n_h,sizeof(phydbl));
  m4mod->multipl_unscaled = (phydbl *)mCalloc(n_h,sizeof(phydbl));
  m4mod->h_fq_unscaled = (phydbl *)mCalloc(n_h,sizeof(phydbl));
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


/* Fill in the (big) rate matrix of the M4 t_mod */ 
void M4_Update_Qmat(m4 *m4mod, t_mod *mod)
{
  int i,j;
  int n_s, n_o, n_h;
  phydbl mr, sum;

  /* The number of states in M4 models is the product 
     of the number of hidden states (or classes) by the
     number of observable states 
   */
  n_s = mod->ns;
  n_o = m4mod->n_o;
  n_h = m4mod->n_h;
  
  /* Set the relative substitution rates */
  if(mod->m4mod->use_cov_alpha)
    {
      DiscreteGamma(m4mod->h_fq,m4mod->multipl,m4mod->alpha,m4mod->alpha,m4mod->n_h,mod->ras->gamma_median);
    }
  else if(mod->m4mod->use_cov_free)
    {
      sum = .0;
      For(i,mod->m4mod->n_h) sum += FABS(mod->m4mod->h_fq_unscaled[i]);
      For(i,mod->m4mod->n_h) mod->m4mod->h_fq[i] = FABS(mod->m4mod->h_fq_unscaled[i])/sum;
      
      do
	{
	  sum = .0;
	  For(i,mod->m4mod->n_h)
	    {
	      if(mod->m4mod->h_fq[i] < 0.01) mod->m4mod->h_fq[i]=0.01;
	      if(mod->m4mod->h_fq[i] > 0.99) mod->m4mod->h_fq[i]=0.99;
	      sum += mod->m4mod->h_fq[i];
	    }

	  For(i,mod->m4mod->n_h) mod->m4mod->h_fq[i]/=sum;
	}
      while((sum > 1.01) || (sum < 0.99));


      /* Make sure the multipliers are centered on 1.0 */
      sum = .0;
      For(i,mod->m4mod->n_h) sum += FABS(mod->m4mod->multipl_unscaled[i]) * mod->m4mod->h_fq[i];
      For(i,mod->m4mod->n_h) mod->m4mod->multipl[i] = mod->m4mod->multipl_unscaled[i] / sum;
      
      /* printf("\n. WARNING\n"); */
      /* mod->m4mod->h_fq[0] = 1./3; */
      /* mod->m4mod->h_fq[1] = 1./3; */
      /* mod->m4mod->h_fq[2] = 1./3; */

      /* mod->m4mod->multipl[0] = 1.0; */
      /* mod->m4mod->multipl[1] = 1.0; */
      /* mod->m4mod->multipl[2] = 1.0; */

      sum = 0;
      For(i,mod->m4mod->n_h) sum += mod->m4mod->multipl[i] * mod->m4mod->h_fq[i];
      if(sum < 0.99 || sum > 1.01)
	{
	  PhyML_Printf("\n. sum = %f",sum);
	  PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
	  Warn_And_Exit("\n");
	}

      /* PhyML_Printf("\n__ "); */
      /* For(i,mod->m4mod->n_h) PhyML_Printf("\n.%f %f %f", */
      /* 				    mod->m4mod->h_fq[i], */
      /* 				    mod->m4mod->h_fq_unscaled[i], */
      /* 				    mod->m4mod->multipl[i]); */
    }


  /* PhyML_Printf("\n."); */
  /* PhyML_Printf("\n. M4 model parameters"); */
  /* m4mod->delta=.0; */
  /* PhyML_Printf("\n. Delta = %f",m4mod->delta); */
  /* For(i,mod->m4mod->n_h) PhyML_Printf("\n. multipl %d = %f",i,m4mod->multipl[i]); */
  /* For(i,mod->m4mod->n_h) PhyML_Printf("\n. fq %d = %f",i,m4mod->h_fq[i]); */


  /* Set up the stationary frequency vector */
  For(i,n_s) mod->e_frq->pi->v[i] = m4mod->o_fq[i%n_o] * m4mod->h_fq[i/n_o];


  if(mod->whichmodel != CUSTOM &&
     mod->whichmodel != GTR    &&
     mod->io->datatype == NT)    
    {
      phydbl kappa1,kappa2;

      if((mod->whichmodel != F84) && (mod->whichmodel != TN93)) mod->lambda->v = 1.; 
      else if(mod->whichmodel == F84)
	{
	  mod->lambda->v = Get_Lambda_F84(mod->e_frq->pi->v,&(mod->kappa->v));
	}

      kappa2 = mod->kappa->v*2./(1.+mod->lambda->v);
      kappa1 = kappa2 * mod->lambda->v;

      /* A <-> C */ m4mod->o_rr[0] = 1.0;
      /* A <-> G */ m4mod->o_rr[1] = kappa2;
      /* A <-> T */ m4mod->o_rr[2] = 1.0;
      /* C <-> G */ m4mod->o_rr[3] = 1.0;
      /* C <-> T */ m4mod->o_rr[4] = kappa1;
    }

  /* Fill in the matrices of nucleotide or amino-acid substitution rates here */
  Update_Qmat_Generic(m4mod->o_rr, m4mod->o_fq, m4mod->n_o, m4mod->o_mats[0]);

  /* Print_Square_Matrix_Generic(n_o,m4mod->o_mats[0]); */

  /* Multiply each of these matrices by a relative substitution rate */
  for(i=1;i<m4mod->n_h;i++) For(j,n_o*n_o) m4mod->o_mats[i][j] = m4mod->o_mats[0][j]*m4mod->multipl[i];
  For(j,n_o*n_o) m4mod->o_mats[0][j] *= m4mod->multipl[0];

  For(i,n_s*n_s) mod->r_mat->qmat->v[i] = .0;

  /* Diagonal blocks (i.e, nucleotide substitutions), symmetric */
  For(i,n_s)
    {
      for(j=i+1;j<n_s;j++)
	{
	  if((int)(j/n_o) == (int)(i/n_o))
	    {
	      mod->r_mat->qmat->v[i*n_s+j] = m4mod->o_mats[(int)(i/n_o)][(i%n_o)*n_o+j%n_o];
	      mod->r_mat->qmat->v[j*n_s+i] = mod->r_mat->qmat->v[i*n_s+j] * m4mod->o_fq[i%n_o] / m4mod->o_fq[j%n_o];
	    }
	}
    }

  /* Work out scaling factor such that the  expected number of observed state substitution
     along a branch of length 1 is 1.*/
  mr = .0;
  For(i,n_s)
    {
      sum = .0;
      For(j,n_s) sum += mod->r_mat->qmat->v[i*n_s+j];
      mr += sum * m4mod->o_fq[i%n_o] * m4mod->h_fq[(int)(i/n_o)]; 
    }
  
  /* Scale the diagonal blocks */
  For(i,n_s*n_s) mod->r_mat->qmat->v[i] /= mr;
  
  /* We are done with the diagonal blocks. Let's fill the non-diagonal ones now. */

  /* Fill the matrix of substitution rate across classes (switches) here */
  Update_Qmat_Generic(m4mod->h_rr, m4mod->h_fq, m4mod->n_h, m4mod->h_mat);

/*   Print_Square_Matrix_Generic(m4mod->n_h,m4mod->h_mat); */

  /* Multiply this matrix by the switching rate */
  For(i,n_h*n_h) m4mod->h_mat[i] *= m4mod->delta;

  /* Fill the non diagonal blocks */
  For(i,n_s)
    {
      for(j=i+1;j<n_s;j++)
	{
	  if((int)(j/n_o) != (int)(i/n_o))
	    {
	      if(i%n_o == j%n_o)
		{
		  mod->r_mat->qmat->v[i*n_s+j] = m4mod->h_mat[(int)(i/n_o)*n_h+(int)(j/n_o)];
		  mod->r_mat->qmat->v[j*n_s+i] = mod->r_mat->qmat->v[i*n_s+j] * m4mod->h_fq[(int)(i/n_o)] / m4mod->h_fq[(int)(j/n_o)]; 
		}
	    }
	}
    }

  /* Note: class equilibrium frequencies are already built in the h_mat matrix.
     No need to 'add' these frequencies later on. */

  /* We are done with the non diagonal blocks */


  /* Diagonal cells */
  For(i,n_s)
    {
      sum = .0;
      For(j,n_s)
	{
	  if(j != i)
	    sum += mod->r_mat->qmat->v[i*n_s+j];
	}
      mod->r_mat->qmat->v[i*n_s+i] = -sum;
    }

  For(i,n_s*n_s) mod->eigen->q[i] = mod->r_mat->qmat->v[i];
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void M4_Init_P_Lk_Tips_Double(t_tree *tree)
{
  int curr_site,i,j,k,l,dim1,dim2,dim3;
  
  dim1 = tree->mod->ras->n_catg * tree->mod->m4mod->n_h * tree->mod->m4mod->n_o;
  dim2 = tree->mod->m4mod->n_h * tree->mod->m4mod->n_o;
  dim3 = tree->mod->m4mod->n_o;


  Fors(curr_site,tree->data->crunch_len,tree->mod->io->state_len)
    {
      For(i,tree->n_otu)
	{
	  for(j=1;j<tree->mod->m4mod->n_h;j++)
	    {
	      For(k,tree->mod->m4mod->n_o)
		{
		  tree->a_nodes[i]->b[0]->p_lk_rght[curr_site*dim1 + 0*dim2 + j*dim3+k] = 
		    tree->a_nodes[i]->b[0]->p_lk_rght[curr_site*dim1 + 0*dim2 + 0*dim3+k];
		  
		  printf("\n() i=%d plk=%f",
			 curr_site*dim1 + 0*dim2 + j*dim3+k,
			 tree->a_nodes[i]->b[0]->p_lk_rght[curr_site*dim1 + 0*dim2 + j*dim3+k]);

		  /* tree->a_nodes[i]->b[0]->p_lk_rght[curr_site][0][j*tree->mod->m4mod->n_o+k] =  */
		  /* tree->a_nodes[i]->b[0]->p_lk_rght[curr_site][0][k]; */
		}


	      For(k,tree->mod->m4mod->n_o)
		for(l=1;l<tree->mod->ras->n_catg;l++)
		  tree->a_nodes[i]->b[0]->p_lk_rght[curr_site*dim1 + l*dim2 + j*dim3+k] = 
		  tree->a_nodes[i]->b[0]->p_lk_rght[curr_site*dim1 + 0*dim2 + j*dim3+k];
		  /* tree->a_nodes[i]->b[0]->p_lk_rght[curr_site][l][j*tree->mod->m4mod->n_o+k] =  */
		  /* tree->a_nodes[i]->b[0]->p_lk_rght[curr_site][0][j*tree->mod->m4mod->n_o+k]; */
	    }
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void M4_Init_P_Lk_Tips_Int(t_tree *tree)
{
  int curr_site,i,j,k,dim2,dim3;

  dim2 = tree->mod->m4mod->n_h * tree->mod->m4mod->n_o;
  dim3 = tree->mod->m4mod->n_o;

  Fors(curr_site,tree->data->crunch_len,tree->mod->io->state_len)
    {
      For(i,tree->n_otu)
	{
	  for(j=1;j<tree->mod->m4mod->n_h;j++)
	    {
	      For(k,tree->mod->m4mod->n_o)
		{
		  tree->a_nodes[i]->b[0]->p_lk_tip_r[curr_site*dim2 + j*dim3+k] = 
		    tree->a_nodes[i]->b[0]->p_lk_tip_r[curr_site*dim2 + 0*dim3+k];
		  /* tree->a_nodes[i]->b[0]->p_lk_tip_r[curr_site][j*tree->mod->m4mod->n_o+k] =  */
		  /* tree->a_nodes[i]->b[0]->p_lk_tip_r[curr_site][k]; */
		}
	    }
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl ****M4_Integral_Term_On_One_Edge(t_edge *b, t_tree *tree)
{
  phydbl ****integral,*P1,*P2;  
  int ns;
  int g,i,j,k,l;
  int step;
  

  ns = tree->mod->ns;


  P1 = (phydbl *)mCalloc(tree->mod->ras->n_catg*ns*ns,sizeof(phydbl));
  P2 = (phydbl *)mCalloc(tree->mod->ras->n_catg*ns*ns,sizeof(phydbl));


  integral = (phydbl ****)mCalloc(tree->mod->ras->n_catg,sizeof(phydbl ***));
  For(g,tree->mod->ras->n_catg)
    {
      integral[g] = (phydbl ***)mCalloc(ns,sizeof(phydbl **));
      For(j,ns)
	{
	  integral[g][j] = (phydbl **)mCalloc(ns,sizeof(phydbl *));
	  For(k,ns) integral[g][j][k] = (phydbl *)mCalloc(ns,sizeof(phydbl));
	}
    }

  /* Integral calculation */
  step = 100;

  PhyML_Printf("\n. [");
  For(i,step)
    {
      For(g,tree->mod->ras->n_catg)
	{
	  PMat(((phydbl)(i+0.5)/step)*b->l->v*tree->mod->ras->gamma_rr->v[g],tree->mod,g*ns*ns,P1);
	  PMat(((phydbl)(step-i-0.5)/step)*b->l->v*tree->mod->ras->gamma_rr->v[g],tree->mod,g*ns*ns,P2);

	  For(j,ns)
	    {
	      For(k,ns)
		{
		  For(l,ns)
		    {
		      /* integral[g][j][k][l] += P1[g][j][k] * P2[g][j][l]  / ((phydbl)(step)); */
		      integral[g][j][k][l] += P1[g*ns*ns + j*ns+k] * P2[g*ns*ns + j*ns+l] / ((phydbl)(step));
		    }
		}
	    }      
	}
      PhyML_Printf("."); fflush(NULL);
    }
  PhyML_Printf("]\n");

  Free(P1);
  Free(P2);

  return integral;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void M4_Post_Prob_H_Class_Edge_Site(t_edge *b, phydbl ****integral, phydbl *postprob, t_tree *tree)
{
  /* Calculation of the expected frequencies of each hidden
     class at a given site. */

  phydbl site_lk;
  int g,i,j,k,l;
  int n_h;
  phydbl sum;
  int dim1,dim2;

  dim1 = tree->mod->ras->n_catg * tree->mod->ns;
  dim2 = tree->mod->ns;

  n_h = tree->mod->m4mod->n_h; /* number of classes, i.e., number of hidden states */

  site_lk = (phydbl)EXP(tree->cur_site_lk[tree->curr_site]);

  if(b->rght->tax)
    {
      sum = .0;
      For(i,n_h)
	{
	  postprob[i] = .0;
	  For(j,tree->mod->m4mod->n_o)
	    {
	      For(g,tree->mod->ras->n_catg)
		{
		  For(k,tree->mod->ns)
		    {
		      For(l,tree->mod->ns)
			{
			  postprob[i] +=

			    (1./site_lk) *
			    tree->mod->ras->gamma_r_proba->v[g] *
			    tree->mod->m4mod->h_fq[i] *
			    tree->mod->m4mod->o_fq[j] *
			    b->p_lk_left[tree->curr_site*dim1 + g*dim2 + k] *
			    b->p_lk_tip_r[tree->curr_site*dim2 + l] *
			    integral[g][i*tree->mod->m4mod->n_o+j][k][l];

			    /* (1./site_lk) * */
			    /* tree->mod->ras->gamma_r_proba[g] * */
			    /* tree->mod->m4mod->h_fq[i] * */
			    /* tree->mod->m4mod->o_fq[j] * */
			    /* b->p_lk_left[tree->curr_site][g][k] * */
			    /* b->p_lk_tip_r[tree->curr_site][l] * */
			    /* /\* 			b->p_lk_rght[tree->curr_site][0][l] * *\/ */
			    /* integral[g][i*tree->mod->m4mod->n_o+j][k][l]; */
			}
		    }
		}
	    }
	  sum += postprob[i];
	}

      /* TO DO */
      For(i,n_h) postprob[i] *= EXP(b->sum_scale_left[tree->curr_site]); 

    }
  else
    {
      sum = .0;
      For(i,n_h)
	{
	  postprob[i] = .0;
	  For(j,tree->mod->m4mod->n_o)
	    {
	      For(g,tree->mod->ras->n_catg)
		{
		  For(k,tree->mod->ns)
		    {
		      For(l,tree->mod->ns)
			{
			  postprob[i] +=

			    (1./site_lk) *
			    tree->mod->ras->gamma_r_proba->v[g] *
			    tree->mod->m4mod->h_fq[i] *
			    tree->mod->m4mod->o_fq[j] *
			    b->p_lk_left[tree->curr_site*dim1 + g*dim2 + k] *
			    b->p_lk_rght[tree->curr_site*dim1 + g*dim2 + l] *
			    integral[g][i*tree->mod->m4mod->n_o+j][k][l];

			    /* (1./site_lk) * */
			    /* tree->mod->ras->gamma_r_proba[g] * */
			    /* tree->mod->m4mod->h_fq[i] * */
			    /* tree->mod->m4mod->o_fq[j] * */
			    /* b->p_lk_left[tree->curr_site][g][k] * */
			    /* b->p_lk_rght[tree->curr_site][g][l] * */
			    /* integral[g][i*tree->mod->m4mod->n_o+j][k][l]; */
			}
		    }
		}
	    }
	  sum += postprob[i];
	}

      /* TO DO */
      For(i,n_h) postprob[i] *= EXP(b->sum_scale_left[tree->curr_site] + b->sum_scale_rght[tree->curr_site]); 

    }

  For(i,n_h) 
    if((postprob[i] < -1.E-5) || (postprob[i] > 1.0+1.E-5))
      {
	PhyML_Printf("\n. Cat : %d Prob : %f\n",i,postprob[i]);
	PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
	Warn_And_Exit("\n");
      }

  sum = 0.0;
  For(i,n_h) sum += postprob[i];

  if((sum > 1.0+1.E-2) || (sum < 1.0-1.E-2))
    {
      PhyML_Printf("\n. Sum = %f\n",sum);
      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
      Exit("\n");
    }

  return;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl ***M4_Compute_Proba_Hidden_States_On_Edges(t_tree *tree)
{
  int i;
  phydbl ***post_probs;
  phydbl ****integral;


  post_probs = (phydbl ***)mCalloc(2*tree->n_otu-3,sizeof(phydbl **));

  For(i,2*tree->n_otu-3)
    {
      post_probs[i] = (phydbl **)mCalloc(tree->n_pattern,sizeof(phydbl *));
      For(tree->curr_site,tree->n_pattern) 
	post_probs[i][tree->curr_site] = (phydbl *)mCalloc(tree->mod->m4mod->n_h,sizeof(phydbl));
    }


  /* Compute posterior probabilities of each hidden class (usually a rate class) 
     on each edge, at each site.
  */
  For(i,2*tree->n_otu-3) 
    {
      PhyML_Printf("\n. Edge %4d/%4d",i+1,2*tree->n_otu-3);

      integral = M4_Integral_Term_On_One_Edge(tree->a_edges[i],tree);

      For(tree->curr_site,tree->n_pattern)
	M4_Post_Prob_H_Class_Edge_Site(tree->a_edges[i],
				       integral,
				       post_probs[i][tree->curr_site],
				       tree);
      
      M4_Free_Integral_Term_On_One_Edge(integral,tree);
    }
  return post_probs;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


/* Estimate the (posterior) mean relative rate of substitution on each branch
   at each site. The posterior mean rates averaged over sites is also estimated
   for each edge. The corresponding trees are printed in a postscript file. Tree 0
   is the tree with posterior mean rates averaged over the sites. The following trees
   have posterior mean rates computed for each site.
*/
void M4_Compute_Posterior_Mean_Rates(phydbl ***post_probs, t_tree *tree)
{
  char *s;
  int i;
  phydbl **mean_post_probs;
  phydbl *mrr;
  phydbl sum;
  int patt,br,rcat;
  phydbl *mean_br_len;
  int best_r,len_var;
  phydbl max_prob;

  mean_br_len = (phydbl *)mCalloc(2*tree->n_otu-3,sizeof(phydbl));
  mean_post_probs = (phydbl **)mCalloc(2*tree->n_otu-3,sizeof(phydbl *));
  For(i,2*tree->n_otu-3) mean_post_probs[i] = (phydbl *)mCalloc(tree->mod->m4mod->n_h,sizeof(phydbl ));
  mrr = (phydbl *)mCalloc(2*tree->n_otu-3,sizeof(phydbl));

  Record_Br_Len(tree);
  M4_Scale_Br_Len(tree);

  /* Compute the posterior mean relative rate on each branch averaged over the 
     whole set of patterns (sites) */
  len_var = 0;
  For(patt,tree->n_pattern) 
    {
      if(!Is_Invar(patt,1,NT,tree->data))
	{
	  For(br,2*tree->n_otu-3)
	    {
	      max_prob = -1.;
	      best_r = -1;
	      For(rcat,tree->mod->m4mod->n_h)
		{
		  if(post_probs[br][patt][rcat] > max_prob)
		    {
		      max_prob = post_probs[br][patt][rcat];
		      best_r = rcat;
		    }
		}

/* /\* 	      Add weight on each category, weight is proportional to the corresponding posterior probability *\/ */
/* 	      For(rcat,tree->mod->m4mod->n_h) */
/* 		{ */
/* 		  mean_post_probs[br][rcat] += post_probs[br][patt][rcat] * tree->data->wght[patt]; */
/* 		} */

	      /* Add weight on the most probable rate category only */
	      mean_post_probs[br][best_r] += tree->data->wght[patt];
	    }
	  len_var += tree->data->wght[patt];
	}
    }

  For(br,2*tree->n_otu-3) 
    {
      For(rcat,tree->mod->m4mod->n_h)
	{
	  mean_post_probs[br][rcat] /= (phydbl)len_var;
	}
    }

  /* Compute the posterior mean relative rate and scale
     each branch length using this factor */
  For(br,2*tree->n_otu-3)
    {
      For(rcat,tree->mod->m4mod->n_h)
	{
	  mrr[br] += mean_post_probs[br][rcat] * tree->mod->m4mod->multipl[rcat];
	}
      tree->a_edges[br]->l->v *= mrr[br];
    }

  PhyML_Fprintf(tree->io->fp_out_stats,"\n. Mean posterior probabilities & rates\n");
  For(rcat,tree->mod->m4mod->n_h) PhyML_Fprintf(tree->io->fp_out_stats,"%2.4f ",tree->mod->m4mod->multipl[rcat]);
  PhyML_Fprintf(tree->io->fp_out_stats,"\n");
  For(br,2*tree->n_otu-3) 
    {
      For(rcat,tree->mod->m4mod->n_h)
	{
	  PhyML_Fprintf(tree->io->fp_out_stats,"%2.4f ",mean_post_probs[br][rcat]);
	}
/*       PhyML_Fprintf(tree->io->fp_out_stats," -- %f -> %f x %f = %f",mrr[br],tree->a_edges[br]->l->v,mrr[br],tree->a_edges[br]->l->v*mrr[br]); */

      PhyML_Fprintf(tree->io->fp_out_stats," mrr=%f ",mrr[br]);

      if(mrr[br] > 1.) PhyML_Fprintf(tree->io->fp_out_stats,"FAST ");
      else             PhyML_Fprintf(tree->io->fp_out_stats,"SLOW ");
      
      PhyML_Fprintf(tree->io->fp_out_stats,"%s",tree->a_edges[br]->labels[0]);

      PhyML_Fprintf(tree->io->fp_out_stats,"\n");
    }

  /* Print the tree */
  PhyML_Fprintf(tree->io->fp_out_tree,"Constrained tree with corrected branch lengths = ");
  s = Write_Tree(tree,NO);
  PhyML_Fprintf(tree->io->fp_out_tree,"%s\n",s);
  Free(s);
  tree->ps_tree = DR_Make_Tdraw_Struct(tree);
  DR_Print_Postscript_Header(tree->n_pattern,tree->io->fp_out_ps);
  tree->ps_page_number = 0;
  DR_Print_Tree_Postscript(tree->ps_page_number++,YES,tree->io->fp_out_ps,tree);

  /* Go back to the initial scaled branch lengths */
  For(br,2*tree->n_otu-3) tree->a_edges[br]->l->v /= mrr[br];

  /* Compute the posterior mean relative rate at each site, for each branch
     and each rate category. Scale branch lengths using these factors and
     print each tree (i.e., on tree per site pattern) */
  For(patt,tree->n_pattern) 
    {
      For(br,2*tree->n_otu-3) 
	{
	  mrr[br] = .0;
	  max_prob = -1.;
	  best_r = -1;
	  For(rcat,tree->mod->m4mod->n_h) /* For each rate class */
	    {
	      mrr[br] += post_probs[br][patt][rcat] * tree->mod->m4mod->multipl[rcat];
	      if(post_probs[br][patt][rcat] > max_prob)
		{
		  max_prob = post_probs[br][patt][rcat];
		  best_r = rcat;
		}
	    }
/* 	  mrr[br] = tree->mod->m4mod->multipl[best_r]; /\* Use the most probable relative rate instead of mean *\/ */
	  tree->a_edges[br]->l->v *= mrr[br];
	}

      For(br,2*tree->n_otu-3) mean_br_len[br] += tree->a_edges[br]->l->v * tree->data->wght[patt];

      PhyML_Fprintf(tree->io->fp_out_stats,"\n. Posterior probabilities site %4d (weight=%d, is_inv=%d)\n",
	     patt,
	     tree->data->wght[patt],
	     Is_Invar(patt,1,NT,tree->data));

      For(rcat,tree->mod->m4mod->n_h) PhyML_Fprintf(tree->io->fp_out_stats,"%2.4f ",tree->mod->m4mod->multipl[rcat]);
      PhyML_Fprintf(tree->io->fp_out_stats,"\n");
      For(br,2*tree->n_otu-3)
	{
	  PhyML_Fprintf(tree->io->fp_out_stats,"Edge %3d ",br);
	  max_prob = -1.0;
	  best_r = -1;
	  For(rcat,tree->mod->m4mod->n_h)
	    {
	      if(post_probs[br][patt][rcat] > max_prob)
		{
		  max_prob = post_probs[br][patt][rcat];
		  best_r = rcat;
		}
	    }

	  For(rcat,tree->mod->m4mod->n_h)
	    {
	      PhyML_Fprintf(tree->io->fp_out_stats,"%2.4f",post_probs[br][patt][rcat]);
	      if(rcat == best_r) PhyML_Fprintf(tree->io->fp_out_stats,"* ");
	      else               PhyML_Fprintf(tree->io->fp_out_stats,"  ");
	    }

/* 	  PhyML_Fprintf(tree->io->fp_out_stats," -- %f -> %f x %f = %f",mrr[br],tree->a_edges[br]->l->v,mrr[br],tree->a_edges[br]->l->v*mrr[br]); */
	  
	  if(mrr[br] > 1.01)      PhyML_Fprintf(tree->io->fp_out_stats," %s ","FAST");
	  else if(mrr[br] < 0.99) PhyML_Fprintf(tree->io->fp_out_stats," %s ","SLOW");
	  else 	                  PhyML_Fprintf(tree->io->fp_out_stats," %s ","MEDIUM");
	  PhyML_Fprintf(tree->io->fp_out_stats,"%s ",tree->a_edges[br]->labels[0]);
	  PhyML_Fprintf(tree->io->fp_out_stats,"\n");
	}

      PhyML_Fprintf(tree->io->fp_out_tree,"tree %d = ",patt+1);
      s = Write_Tree(tree,NO);
      PhyML_Fprintf(tree->io->fp_out_tree,"%s\n",s);
      Free(s);
      DR_Print_Tree_Postscript(tree->ps_page_number++,YES,tree->io->fp_out_ps,tree);

      /* Go back to the initial scaled branch lengths */
      For(br,2*tree->n_otu-3) tree->a_edges[br]->l->v /= mrr[br];

      For(br,2*tree->n_otu-3) 
	{
	  sum = .0;
	  For(rcat,tree->mod->m4mod->n_h)
	    {
	      sum += post_probs[br][patt][rcat];
	    }
	  
	  if((sum < 0.99) || (sum > 1.01))
	    {
	      PhyML_Fprintf(tree->io->fp_out_stats,"\n. sum = %f\n",sum);
	      PhyML_Fprintf(tree->io->fp_out_stats,"\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
	      Warn_And_Exit("\n");
	    }
	}
    }
  
  /* Mean branch lengths */
  For(br,2*tree->n_otu-3)
    {
      mean_br_len[br] /= (phydbl)tree->data->init_len;
      tree->a_edges[br]->l->v = mean_br_len[br];
    }
  PhyML_Fprintf(tree->io->fp_out_tree,"Mean branch lengths=");
  s = Write_Tree(tree,NO);
  PhyML_Fprintf(tree->io->fp_out_tree,"%s\n",s);
  Free(s);
/*   DR_Print_Tree_Postscript(tree->ps_page_number++,tree->io->fp_out_ps,tree); */

  Restore_Br_Len(tree);

  DR_Print_Postscript_EOF(tree->io->fp_out_ps);

  For(br,2*tree->n_otu-3)
    {
      For(tree->curr_site,tree->n_pattern)
	Free(post_probs[br][tree->curr_site]);
      Free(post_probs[br]);
    }
  Free(post_probs);
  For(i,2*tree->n_otu-3) Free(mean_post_probs[i]);
  Free(mean_post_probs);
  Free(mrr);
  Free(mean_br_len);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


/* Classifiy each branch, at each site, among one of the rate classes */
phydbl **M4_Site_Branch_Classification(phydbl ***post_probs, t_tree *tree)
{
  int patt, br, rcat, i;
  phydbl **best_probs;
  phydbl post_prob_fast, post_prob_slow;

  best_probs = (phydbl **)mCalloc(tree->n_pattern,sizeof(phydbl *));
  For(i,tree->n_pattern) best_probs[i] = (phydbl *)mCalloc(2*tree->n_otu-3,sizeof(phydbl));

  tree->write_labels = YES;

  For(patt,tree->n_pattern) 
    {
      For(br,2*tree->n_otu-3) 
	{
	  post_prob_fast = .0;
	  post_prob_slow = .0;

	  For(rcat,tree->mod->m4mod->n_h) /* For each rate class */
	    {	      
	      if(tree->mod->m4mod->multipl[rcat] > 1.0) 
		post_prob_fast += post_probs[br][patt][rcat];
	      else
		post_prob_slow += post_probs[br][patt][rcat];
	    }

	  best_probs[patt][br] = (post_prob_fast > post_prob_slow)?(post_prob_fast):(post_prob_slow);

	  if(!(tree->a_edges[br]->n_labels%BLOCK_LABELS)) Make_New_Edge_Label(tree->a_edges[br]);

/* 	  if((post_prob_fast > post_prob_slow) && (best_probs[patt][br] > 0.95)) */
/* 	    strcpy(tree->a_edges[br]->labels[tree->a_edges[br]->n_labels],"FASTER"); */
/* 	  else if((post_prob_fast < post_prob_slow) && (best_probs[patt][br] > 0.95)) */
/* 	    strcpy(tree->a_edges[br]->labels[tree->a_edges[br]->n_labels],"SLOWER"); */
/* 	  else */
/* 	    strcpy(tree->a_edges[br]->labels[tree->a_edges[br]->n_labels],"UNKNOWN"); */

	  if(post_prob_fast > post_prob_slow)
	    strcpy(tree->a_edges[br]->labels[tree->a_edges[br]->n_labels],"FASTER");
	  else 
	    strcpy(tree->a_edges[br]->labels[tree->a_edges[br]->n_labels],"SLOWER");

	  tree->a_edges[br]->n_labels++;
	}
    }
  return best_probs;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void M4_Site_Branch_Classification_Experiment(t_tree *tree)
{
  calign *cpy_data;
  short int **true_rclass, **est_rclass;
  phydbl **best_probs;
  int i,j;
  phydbl correct_class, mis_class, unknown;
  
  true_rclass = (short int **)mCalloc(tree->data->init_len, sizeof(short int *));
  est_rclass  = (short int **)mCalloc(tree->data->init_len, sizeof(short int *));
 
  For(i,tree->data->init_len)
    {
      true_rclass[i] = (short int *)mCalloc(2*tree->n_otu-3,sizeof(short int));
      est_rclass[i]  = (short int *)mCalloc(2*tree->n_otu-3,sizeof(short int));
    }

  if(tree->io->datatype != NT && tree->io->datatype != AA)
    {
      PhyML_Printf("\n. Not implemented yet.");
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
  
  cpy_data = Copy_Cseq(tree->data,tree->io);

  /* Generate a simulated data set under H0, with the right sequence length. */
  PhyML_Printf("\n. Evolving sequences (delta=%f, alpha=%f) ...\n",tree->mod->m4mod->delta,tree->mod->m4mod->alpha);
  Evolve(cpy_data,tree->mod,tree);

  For(i,cpy_data->init_len)
    {
      For(j,2*tree->n_otu-3)
	{
	  if(!strcmp(tree->a_edges[j]->labels[i],"FASTER"))
	    {
	      true_rclass[i][j] = 1;
	    }
	  else if(!strcmp(tree->a_edges[j]->labels[i],"SLOWER"))
	    {
	      true_rclass[i][j] = 0;
	    }
	  else
	    {
	      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
	      Warn_And_Exit("\n");
	    }
	}
    }
  
  For(j,2*tree->n_otu-3) 
    {
      Free_Edge_Labels(tree->a_edges[j]);
      tree->a_edges[j]->n_labels = 0;
    }

  /* Generate the memory needed for likelihood calculation because
     we will need bigger arrays
  */
  Free_Tree_Lk(tree);
  Free_Tree_Pars(tree);

  tree->data      = cpy_data;
  tree->n_pattern = tree->data->crunch_len;

  /* Allocate memory and initialize likelihood structure with
     simulated sequences (i.e., columns are not compressed)
  */
  Make_Tree_4_Pars(tree,cpy_data,cpy_data->init_len);
  Make_Tree_4_Lk(tree,cpy_data,cpy_data->init_len);

  /* Estimate model parameters */
  PhyML_Printf("\n. Estimating model parameters...\n");
  tree->mod->s_opt->opt_cov_alpha = 1;
  tree->mod->s_opt->opt_cov_delta = 1;
  Round_Optimize(tree,tree->data,ROUND_MAX);

  tree->both_sides = 1;
  Lk(NULL,tree);

  /* Classify branches */
  best_probs = M4_Site_Branch_Classification(M4_Compute_Proba_Hidden_States_On_Edges(tree),tree);

  For(i,tree->data->init_len)
    {
      For(j,2*tree->n_otu-3)
	{
	  if(!strcmp(tree->a_edges[j]->labels[i],"FASTER"))
	    {
	      est_rclass[i][j] = 1;
	    }
	  else if(!strcmp(tree->a_edges[j]->labels[i],"SLOWER"))
	    {
	      est_rclass[i][j] = 0;
	    }
	  else if(!strcmp(tree->a_edges[j]->labels[i],"UNKNOWN"))
	    {
	      est_rclass[i][j] = -1;
	    }
	  else
	    {
	      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
	      Warn_And_Exit("\n");
	    }
	}
    }

  unknown       = .0;
  correct_class = .0;
  mis_class     = .0;
  For(i,tree->data->init_len)
    {
      For(j,2*tree->n_otu-3)
	{
/* 	  PhyML_Printf("\n. Edge %3d %4d %4d - %f",j,true_rclass[i][j],est_rclass[i][j],best_probs[i][j]); */
	  if(est_rclass[i][j] == -1)
	    {
	      unknown += 1.;
	    }
	  else if(est_rclass[i][j] == true_rclass[i][j])
	    {
	      correct_class += 1.;
	    }
	  else if(est_rclass[i][j] != true_rclass[i][j])
	    {
	      mis_class += 1.;
	    }
	  else
	    {
	      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
	      Warn_And_Exit("\n");
	    }
	}
/*       PhyML_Printf("\n"); */
    }

  correct_class /= ((tree->data->init_len * (2*tree->n_otu-3)) - unknown);
  mis_class     /= ((tree->data->init_len * (2*tree->n_otu-3)) - unknown);
  unknown       /= (tree->data->init_len  * (2*tree->n_otu-3));
  
  PhyML_Printf("\n. correct_class = %f mis_class = %f unknown = %f",
	 correct_class,
	 mis_class,
	 unknown);


  For(i,tree->data->init_len)
    {
      Free(true_rclass[i]);
      Free(est_rclass[i]);
      Free(best_probs[i]);
    }
  Free(true_rclass);
  Free(est_rclass);
  Free(best_probs);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


  /* Scale branch lengths such that they express expected number
     of nucleotide or amino-acid substitutions */

void M4_Scale_Br_Len(t_tree *tree)
{
  phydbl scale_fact,mrs;
  int i,j;

  /* (1) Work out the relative mean rate of switches */
  mrs = .0;
  For(i,tree->mod->m4mod->n_h)
    {
      For(j,tree->mod->m4mod->n_h)
	{
	  if(j != i)
	    mrs += tree->mod->m4mod->h_fq[i] * tree->mod->m4mod->h_mat[i*tree->mod->m4mod->n_h+j];
	}
    }

  /* (2) scale_fact = (1 + delta x mrs) */
  scale_fact = 1.0 + tree->mod->m4mod->delta * mrs;

  /* (3) Scale branch lengths */
  For(i,2*tree->n_otu-3) tree->a_edges[i]->l->v /= scale_fact;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void M4_Free_Integral_Term_On_One_Edge(phydbl ****integral, t_tree *tree)
{
  int g,i,j;

  For(g,tree->mod->ras->n_catg)
    {
      For(i,tree->mod->m4mod->n_h)
	{
	  For(j,tree->mod->m4mod->n_h)
	    {
	      Free(integral[g][i][j]);
	    }
	  Free(integral[g][i]);
	}      
      Free(integral[g]);
    }
  Free(integral);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void M4_Detect_Site_Switches_Experiment(t_tree *tree)
{
  t_mod *nocov_mod,*cov_mod,*ori_mod;
  calign *ori_data,*cpy_data;
  int i,n_iter;
  phydbl *nocov_bl,*cov_bl;
  phydbl *site_lnl_nocov, *site_lnl_cov;

  nocov_bl       = (phydbl *)mCalloc(2*tree->n_otu-3,sizeof(phydbl));
  cov_bl         = (phydbl *)mCalloc(2*tree->n_otu-3,sizeof(phydbl));
  site_lnl_nocov = (phydbl *)mCalloc(tree->data->init_len,sizeof(phydbl));
  site_lnl_cov   = (phydbl *)mCalloc(tree->data->init_len,sizeof(phydbl));

  ori_data = tree->data;
  ori_mod  = tree->mod;

  if(tree->io->datatype != NT && tree->io->datatype != AA)
    {
      PhyML_Printf("\n== Not implemented yet.");
      PhyML_Printf("\n== Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  cpy_data = Copy_Cseq(tree->data,tree->io);

  PhyML_Printf("\n. Estimate model parameters under non-switching substitution model.\n");
  Switch_From_M4mod_To_Mod(tree->mod);
  Simu_Loop(tree);
  nocov_mod = (t_mod *)Copy_Model(tree->mod); /* Record model parameters */
  For(i,2*tree->n_otu-3) nocov_bl[i] = tree->a_edges[i]->l->v; /* Record branch lengths */
  For(i,tree->data->crunch_len) site_lnl_nocov[i] = tree->cur_site_lk[i];
  Print_Lk(tree,"[LnL under non-switching substitution model]");
  
  PhyML_Printf("\n. Estimate model parameters under switching substitution model.\n");
  Switch_From_Mod_To_M4mod(tree->mod);
  Simu_Loop(tree);
  cov_mod = (t_mod *)Copy_Model(tree->mod); /* Record model parameters */
  For(i,2*tree->n_otu-3) cov_bl[i] = tree->a_edges[i]->l->v; /* Record branch lengths */
  For(i,tree->data->crunch_len) site_lnl_cov[i] = tree->cur_site_lk[i];
  Print_Lk(tree,"[LnL under switching substitution model]");
  

  PhyML_Printf("\n");
  For(i,tree->data->crunch_len) PhyML_Printf("TRUTH %f %f\n",site_lnl_nocov[i],site_lnl_cov[i]);

  /* Generate a simulated data set under H0, with the right sequence length. */
  tree->mod = nocov_mod;
  Evolve(cpy_data, nocov_mod, tree);

  /* Generate the memory needed for likelihood calculation because
     we will need bigger arrays 
  */
  tree->mod = cov_mod;
  Free_Tree_Lk(tree);
  Free_Tree_Pars(tree);

  tree->data      = cpy_data;
  tree->n_pattern = tree->data->crunch_len;
  tree->mod       = cov_mod;

  /* Allocate memory and initialize likelihood structure with
     simulated sequences (i.e., columns are not compressed)
  */
  Make_Tree_4_Pars(tree,cpy_data,cpy_data->init_len);
  Make_Tree_4_Lk(tree,cpy_data,cpy_data->init_len);

 
  n_iter = 0;
  do
    {
      /* Get the transition proba right to generate sequences */
      tree->mod = nocov_mod;
      For(i,2*tree->n_otu-3) tree->a_edges[i]->l->v = nocov_bl[i];
      For(i,2*tree->n_otu-3) Update_PMat_At_Given_Edge(tree->a_edges[i],tree);
      
      /* Generate sequences */
      Evolve(cpy_data, nocov_mod, tree);
      tree->data = cpy_data;

      if(tree->mod->s_opt->greedy) Init_P_Lk_Tips_Double(tree);
      else Init_P_Lk_Tips_Int(tree);
      
      tree->mod = nocov_mod;
      For(i,2*tree->n_otu-3) tree->a_edges[i]->l->v = nocov_bl[i];
      Lk(NULL,tree);
      For(i,tree->data->crunch_len) site_lnl_nocov[i] = tree->cur_site_lk[i];
      Print_Lk(tree,"[CPY LnL under non-switching substitution model]");

      tree->mod = cov_mod;
      For(i,2*tree->n_otu-3) tree->a_edges[i]->l->v = cov_bl[i];
      Lk(NULL,tree);
      For(i,tree->data->crunch_len) site_lnl_cov[i] = tree->cur_site_lk[i];
      Print_Lk(tree,"[CPY LnL under switching substitution model]");

      PhyML_Printf("\n");
      For(i,tree->data->crunch_len) PhyML_Printf("SYNTH %f %f\n",site_lnl_nocov[i],site_lnl_cov[i]);
    }
  while(++n_iter < 200);

  Free_Tree_Lk(tree);
  Free_Tree_Pars(tree);

  /* Back to the original data set to check that everything is ok */
  tree->data      = ori_data;
  tree->n_pattern = tree->data->crunch_len;

  Make_Tree_4_Pars(tree,ori_data,ori_data->init_len);
  Make_Tree_4_Lk(tree,ori_data,ori_data->init_len);

  tree->mod = nocov_mod;
  For(i,2*tree->n_otu-3) tree->a_edges[i]->l->v = nocov_bl[i];  
  Lk(NULL,tree);
  Print_Lk(tree,"[FINAL LnL under non-switching substitution model]");
  
  tree->mod = cov_mod;
  For(i,2*tree->n_otu-3) tree->a_edges[i]->l->v = cov_bl[i];  
  Lk(NULL,tree);
  Print_Lk(tree,"[FINAL LnL under switching substitution model]");

  tree->mod = ori_mod;

  Free_Model(cov_mod);
  Free_Model(nocov_mod);
  Free(site_lnl_cov);
  Free(site_lnl_nocov);

  Free_Cseq(cpy_data);
  Free(nocov_bl);
  Free(cov_bl);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void M4_Posterior_Prediction_Experiment(t_tree *tree)
{
  calign *ori_data,*cpy_data;
  int i,n_iter,n_simul;
  FILE *fp_nocov,*fp_cov,*fp_obs;
  char *s;
  t_edge *best_edge;

  s = (char *)mCalloc(100,sizeof(char));
  
  best_edge = NULL;

  strcpy(s,tree->io->in_align_file);
  fp_nocov = Openfile(strcat(s,"_nocov"),1);
  strcpy(s,tree->io->in_align_file);
  fp_cov   = Openfile(strcat(s,"_cov"),1);
  strcpy(s,tree->io->in_align_file);
  fp_obs = Openfile(strcat(s,"_obs"),1);
  
  Free(s);

  Print_Diversity_Header(fp_nocov, tree);
  Print_Diversity_Header(fp_cov, tree);
  Print_Diversity_Header(fp_obs, tree);

  ori_data = tree->data;

  if(tree->io->datatype != NT && tree->io->datatype != AA)
   {
      PhyML_Printf("\n. Not implemented yet.");
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  cpy_data = Copy_Cseq(tree->data,tree->io);

  /* Generate a simulated data set under H0, with the right sequence length. */
  Set_Model_Parameters(tree->mod);
  For(i,2*tree->n_otu-3) Update_PMat_At_Given_Edge(tree->a_edges[i],tree);
  Evolve(cpy_data,tree->mod,tree);

  /* Generate the memory needed for likelihood calculation because
     we will need bigger arrays
  */
  Free_Tree_Lk(tree);
  Free_Tree_Pars(tree);

  tree->data      = cpy_data;
  tree->n_pattern = tree->data->crunch_len;

  /* Allocate memory and initialize likelihood structure with
     simulated sequences (i.e., columns are not compressed)
  */
  Make_Tree_4_Pars(tree,cpy_data,cpy_data->init_len);
  Make_Tree_4_Lk(tree,cpy_data,cpy_data->init_len);

  /* Go back to the original data set */
  tree->data      = ori_data;
  tree->n_pattern = ori_data->crunch_len;
  
  if(tree->mod->s_opt->greedy) Init_P_Lk_Tips_Double(tree);
  else Init_P_Lk_Tips_Int(tree);

  PhyML_Printf("\n. Estimate model parameters under non-switching substitution model.\n");
  Switch_From_M4mod_To_Mod(tree->mod);

  tree->bl_from_node_stamps = 1;
  /* best_edge = TIMES_Find_Best_Root_Position(tree); */
  PhyML_Printf("\n. Put root on t_edge %3d",i);
  TIMES_Least_Square_Node_Times(best_edge,tree);
  TIMES_Adjust_Node_Times(tree);
  /* TIMES_Round_Optimize(tree); */

/*   Round_Optimize(tree,tree->data); */
/*   Simu_Loop(tree); */
  Print_Lk(tree,"[LnL under non-switching substitution model]");
  Print_Tree(tree->io->fp_out_tree,tree);
  
  /* Print observed diversities */
  Init_Ui_Tips(tree);
  Site_Diversity(tree);
  Print_Diversity(fp_obs,tree);

  n_simul = 100;
  n_iter = 0;
  do
    {
      Evolve(cpy_data,tree->mod,tree);
      tree->data      = cpy_data;
      tree->n_pattern = cpy_data->init_len;

      if(tree->mod->s_opt->greedy) Init_P_Lk_Tips_Double(tree);
      else Init_P_Lk_Tips_Int(tree);

      Lk(NULL,tree);

      Init_Ui_Tips(tree);
      Site_Diversity(tree);
      Print_Diversity(fp_nocov,tree);

      Print_Lk(tree,"[CPY under non-switching substitution model]");
    }while(++n_iter < n_simul);


  /* Go back to the original data set */
  tree->data      = ori_data;
  tree->n_pattern = ori_data->crunch_len;
  
  if(tree->mod->s_opt->greedy) Init_P_Lk_Tips_Double(tree);
  else Init_P_Lk_Tips_Int(tree);

  PhyML_Printf("\n. Estimate model parameters under switching substitution model.\n");
  Switch_From_Mod_To_M4mod(tree->mod);
  /* TIME_Round_Optimize(tree); */
/*   Round_Optimize(tree,tree->data); */
  /*   Simu_Loop(tree); */
  Print_Lk(tree,"[LnL under switching substitution model]");
  Print_Tree(tree->io->fp_out_tree,tree);

  n_iter = 0;
  do
    {
      Evolve(cpy_data,tree->mod,tree);
      tree->data      = cpy_data;
      tree->n_pattern = cpy_data->init_len;
      if(tree->mod->s_opt->greedy) Init_P_Lk_Tips_Double(tree);
      else Init_P_Lk_Tips_Int(tree);

      Lk(NULL,tree);

      Init_Ui_Tips(tree);
      Site_Diversity(tree);
      Print_Diversity(fp_cov,tree);

      Print_Lk(tree,"[LnL under switching substitution model]");
    }while(++n_iter < n_simul);

  fclose(fp_obs);
  fclose(fp_nocov);
  fclose(fp_cov);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

m4 *M4_Copy_M4_Model(t_mod *ori_mod, m4 *ori_m4mod)
{
  int i,j,n_h, n_o;
  m4 *cpy_m4mod;
  
  if(ori_mod->io->datatype != NT && ori_mod->io->datatype != AA)
    {
      PhyML_Printf("\n== Not implemented yet.");
      PhyML_Printf("\n== Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");
    }


  cpy_m4mod = (m4 *)M4_Make_Light();
  cpy_m4mod->n_o = ori_m4mod->n_o;
  cpy_m4mod->n_h = ori_m4mod->n_h;

  if(ori_mod->use_m4mod)
    {
      M4_Make_Complete(cpy_m4mod->n_h,
		       cpy_m4mod->n_o,
		       cpy_m4mod);
      
      n_h = cpy_m4mod->n_h;
      n_o = cpy_m4mod->n_o;
      
      cpy_m4mod->n_h = ori_m4mod->n_h;
      cpy_m4mod->n_o = ori_m4mod->n_o;
      For(i,n_h) For(j,n_o*n_o) cpy_m4mod->o_mats[i][j] = ori_m4mod->o_mats[i][j];
      For(i,n_h) cpy_m4mod->multipl[i] = ori_m4mod->multipl[i];
      For(i,n_h) cpy_m4mod->multipl_unscaled[i] = ori_m4mod->multipl_unscaled[i];  
      For(i,n_o*n_o) cpy_m4mod->o_rr[i] = ori_m4mod->o_rr[i];
      For(i,n_h*n_h) cpy_m4mod->h_rr[i] = ori_m4mod->h_rr[i];
      For(i,n_h*n_h) cpy_m4mod->h_mat[i] = ori_m4mod->h_mat[i];
      For(i,n_o) cpy_m4mod->o_fq[i] = ori_m4mod->o_fq[i];
      For(i,n_h) cpy_m4mod->h_fq[i] = ori_m4mod->h_fq[i];
      For(i,n_h) cpy_m4mod->h_fq_unscaled[i] = ori_m4mod->h_fq_unscaled[i];
      cpy_m4mod->delta = ori_m4mod->delta;
      cpy_m4mod->alpha = ori_m4mod->alpha;
    }

  return cpy_m4mod;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

