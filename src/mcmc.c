/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homoLOGous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include "mcmc.h"

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MCMC(t_tree *tree)
{
  int move;
  phydbl u;
  int first,secod;
  int i;

  RATES_Set_Clock_And_Nu_Max(tree);
  RATES_Set_Birth_Rate_Boundaries(tree);
 
  if(tree->mcmc->randomize == YES)
    {
      MCMC_Randomize_Birth(tree);
      MCMC_Randomize_Nu(tree);
      MCMC_Randomize_Node_Times(tree); 
      MCMC_Sim_Rate(tree->n_root,tree->n_root->v[2],tree);
      MCMC_Sim_Rate(tree->n_root,tree->n_root->v[1],tree);
      MCMC_Randomize_Node_Rates(tree);
      MCMC_Randomize_Clock_Rate(tree); /* Clock Rate must be the last parameter to be randomized */
      MCMC_Randomize_Rate_Across_Sites(tree);
      MCMC_Randomize_Kappa(tree);
      MCMC_Randomize_Covarion_Rates(tree);
      MCMC_Randomize_Covarion_Switch(tree);
    }
  else
    {
      MCMC_Read_Param_Vals(tree);
    }
/* #ifdef INVITEE */
/*   tree -> rates -> node_height_dens_log_norm_const_update = Norm_Constant_Prior_Times(tree); */
/* #endif */
  Switch_Eigen(YES,tree->mod);

  MCMC_Initialize_Param_Val(tree->mcmc,tree);
 
  Update_Ancestors(tree->n_root,tree->n_root->v[2],tree);
  Update_Ancestors(tree->n_root,tree->n_root->v[1],tree);

  For(i,2*tree->n_otu-2) tree->rates->br_do_updt[i] = YES;  
  RATES_Update_Cur_Bl(tree);
  RATES_Lk_Rates(tree);	

#ifdef INVITEE
  tree->rates->update_time_norm_const = YES;
#endif
  TIMES_Lk_Times(tree); 
#ifdef INVITEE
  tree->rates->update_time_norm_const = NO;
#endif
  /* printf("\n"); */
  /* printf("\n. current calibration: %d ", tree->rates->cur_comb_numb); */
  /* printf("\n"); */
  /* for(i = tree -> n_otu; i < 2 * tree -> n_otu -1; i++) printf("\nJUMP0 Node number:[%d] Lower bound:[%f] Upper bound:[%f] Node time:[%f].", i, */
  /*                                                              tree -> rates -> t_prior_min[i], */
  /*                                                              tree -> rates -> t_prior_max[i], */
  /*                                                              tree -> rates -> nd_t[i]); */
  /* Exit("\n"); */

  Set_Both_Sides(NO,tree);  
  if(tree->mcmc->use_data) Lk(NULL,tree);
  else tree->c_lnL = 0.0;
  Switch_Eigen(NO,tree->mod);
  MCMC_Print_Param(tree->mcmc,tree);
  

  //////////////////
  if(tree->io->mutmap == YES)
    {
      int j;
      char *s,*t;
      FILE *fp;
      
      Make_MutMap(tree);
 
      For(j,tree->n_otu)
	{
	  s = (char *)mCalloc(T_MAX_NAME,sizeof(char));
	  strcpy(s,tree->a_nodes[j]->name);
	  tree->a_nodes[j]->name = s;
	}
      
      s = (char *)mCalloc(T_MAX_NAME,sizeof(char));
      t = (char *)mCalloc(T_MAX_NAME,sizeof(char));
      
      tree->write_tax_names = YES;
      For(i,tree->n_pattern)
	{
	  For(j,tree->n_otu)
	    {
	      strcpy(t,tree->a_nodes[j]->name);
	      s[0]=tree->data->c_seq[j]->state[i];
	      s[1]='\0';
	      strcat(s,"--");
	      sprintf(s+strlen(s),"%.0f",tree->rates->nd_t[j]);
	      /* strcat(s,tree->a_nodes[j]->name); */
	      strcpy(tree->a_nodes[j]->name,s);
	    }
	  
	  strcpy(s,"rosettatree.");
	  sprintf(s+strlen(s),"%d",i);
	  fp = fopen(s,"w");
	  s = Write_Tree(tree,NO);
	  PhyML_Fprintf(fp,"%s",s);
	  fclose(fp);
	  
	  For(j,tree->n_otu) strcpy(tree->a_nodes[j]->name,t);
	}
      Free(s);
      Free(t);
    }


  first = 2;
  secod = 1;
  do
    {
      /* if(tree->mcmc->ess[tree->mcmc->num_move_tree_height] > 100 &&  */
      /* 	 tree->mcmc->ess[tree->mcmc->num_move_nu] > 100          && */
      /* 	 tree->mcmc->ess[tree->mcmc->num_move_clock_r] > 100     && */
      /* 	 tree->mcmc->run > 1000) */
      /* 	{ */
      /* 	  FILE *fp; */
      /* 	  char *s; */

      /* 	  s = (char *)mCalloc(100,sizeof(char)); */

      /* 	  sprintf(s,"simul_par.%d",getpid()); */
      /* 	  fclose(tree->mcmc->out_fp_stats); */
      /* 	  tree->mcmc->out_fp_stats = fopen(s,"w"); */
      /* 	  tree->mcmc->run = 0; */
      /* 	  tree->mcmc->nd_t_digits = 4; */
      /* 	  MCMC_Print_Param(tree->mcmc,tree); */

      /* 	  RATES_Update_Cur_Bl(tree); */
      /* 	  printf("\n. %s",Write_Tree(tree,NO)); */
      /* 	  Evolve(tree->data,tree->mod,tree); */

      /* 	  sprintf(s,"simul_seq.%d",getpid()); */
      /* 	  fp = fopen(s,"w"); */
      /* 	  Print_CSeq(fp,NO,tree->data); */
      /* 	  fflush(NULL); */
      /* 	  fclose(fp); */
      /* 	  Free(s); */

      /* 	  Exit("\n"); */
      /* 	} */


      /* if(tree->mcmc->ess[tree->mcmc->num_move_tree_height] > 100 && */
      /* 	 tree->mcmc->ess[tree->mcmc->num_move_nu] > 100          && */
      /* 	 tree->mcmc->ess[tree->mcmc->num_move_clock_r] > 100     && */
      /* 	 tree->mcmc->run > 1000) */
      /* 	{ */
      /* 	  FILE *fp; */
      /* 	  char *s,*t; */

      /* 	  s = (char *)mCalloc(100,sizeof(char)); */
	  
      /* 	  t = strrchr(tree->io->in_align_file,'.'); */
      /* 	  sprintf(s,"res%s",t); */
      /* 	  fp = fopen(s,"w"); */
      /* 	  fclose(tree->mcmc->out_fp_stats); */
      /* 	  tree->mcmc->out_fp_stats = fopen(s,"w"); */
      /* 	  tree->mcmc->run = 0; */
      /* 	  MCMC_Print_Param(tree->mcmc,tree); */
      /* 	  fclose(fp); */
      /* 	  Free(s); */
      /* 	  Exit("\n"); */
      /* 	} */

      u = Uni();

      For(move,tree->mcmc->n_moves) if(tree->mcmc->move_weight[move] > u) break;
      
      if(u < .5) { first = 2; secod = 1; }
      else       { first = 1; secod = 2; }
 



      /* Clock rate */
      if(!strcmp(tree->mcmc->move_name[move],"clock"))
      	{
      	  For(i,2*tree->n_otu-2) tree->rates->br_do_updt[i] = YES;	  
          MCMC_Clock_R(tree);  
      	}

      /* Nu */
      else if(!strcmp(tree->mcmc->move_name[move],"nu"))
      	{
      	  For(i,2*tree->n_otu-2) tree->rates->br_do_updt[i] = YES;
	  MCMC_Nu(tree);
      	}

      /* Tree height */
      else if(!strcmp(tree->mcmc->move_name[move],"tree_height"))
      	{  
      	  MCMC_Tree_Height(tree);
      	}

      /* Subtree height */
      else if(!strcmp(tree->mcmc->move_name[move],"subtree_height"))
      	{ 
      	  MCMC_Subtree_Height(tree);  
      	}

      /* Subtree rates */
      else if(!strcmp(tree->mcmc->move_name[move],"subtree_rates"))
      	{
      	  MCMC_Subtree_Rates(tree);
      	}

      /* Birth rate */
      else if(!strcmp(tree->mcmc->move_name[move],"birth_rate"))
      	{
      	  MCMC_Birth_Rate(tree);
      	}

      /* Swing rates */
      else if(!strcmp(tree->mcmc->move_name[move],"tree_rates"))
      	{
      	  MCMC_Tree_Rates(tree);
      	}

      else if(!strcmp(tree->mcmc->move_name[move],"updown_nu_cr"))
      	{
      	  MCMC_Updown_Nu_Cr(tree);
      	}

      else if(!strcmp(tree->mcmc->move_name[move],"updown_t_cr"))
      	{
      	  MCMC_Updown_T_Cr(tree);
      	}

      else if(!strcmp(tree->mcmc->move_name[move],"updown_t_br"))
      	{
      	  MCMC_Updown_T_Br(tree);
      	}

      /* Ts/tv ratio */
      else if(!strcmp(tree->mcmc->move_name[move],"kappa"))
      	{
      	  MCMC_Kappa(tree);
	}

      /* Gamma shape parameter */
      else if(!strcmp(tree->mcmc->move_name[move],"ras"))
      	{
      	  MCMC_Rate_Across_Sites(tree);
	}

      /* Covarion change calibration interval */
      else if(!strcmp(tree->mcmc->move_name[move],"jump_calibration"))
      	{
      	  MCMC_Jump_Calibration(tree);
	}

      /* Covarion model parameters */
      else if(!strcmp(tree->mcmc->move_name[move],"cov_rates"))
      	{
      	  MCMC_Covarion_Rates(tree);
	}

      /* Covarion model parameters */
      else if(!strcmp(tree->mcmc->move_name[move],"cov_switch"))
      	{
      	  MCMC_Covarion_Switch(tree);
	}


      /* Times */
      else if(!strcmp(tree->mcmc->move_name[move],"time"))
      	{
          Set_Both_Sides(YES,tree);     
      	  if(tree->mcmc->use_data == YES) Lk(NULL,tree);
          Set_Both_Sides(NO,tree);     

      	  if(tree->mcmc->is == NO || tree->rates->model_log_rates == YES)
      	    {
              MCMC_Root_Time(tree);
	      MCMC_One_Time(tree->n_root,tree->n_root->v[first],YES,tree);
	      MCMC_One_Time(tree->n_root,tree->n_root->v[secod],YES,tree);
      	    }
      	  else
      	    {
	      /* MCMC_One_Time(tree->n_root,tree->n_root->v[first],YES,tree); */
	      /* MCMC_One_Time(tree->n_root,tree->n_root->v[secod],YES,tree); */
      	      RATES_Posterior_One_Time(tree->n_root,tree->n_root->v[first],YES,tree);
      	      RATES_Posterior_One_Time(tree->n_root,tree->n_root->v[secod],YES,tree);
      	    }
      	}
      
      /* Node Rates */

      else if(!strcmp(tree->mcmc->move_name[move],"nd_rate"))
      	{
      	  MCMC_One_Node_Rate(tree->n_root,tree->n_root->v[first],YES,tree);
      	  MCMC_One_Node_Rate(tree->n_root,tree->n_root->v[secod],YES,tree);
      	}

      /* Edge Rates */
      else if(!strcmp(tree->mcmc->move_name[move],"br_rate"))
      	{

      	  Set_Both_Sides(YES,tree);
      	  if(tree->mcmc->use_data == YES) Lk(NULL,tree);
      	  Set_Both_Sides(NO,tree);
      	  
      	  if(tree->mcmc->is == NO)
      	    {
      	      /* MCMC_Slice_One_Rate(tree->n_root,tree->n_root->v[first],YES,tree); */
      	      /* MCMC_Slice_One_Rate(tree->n_root,tree->n_root->v[secod],YES,tree); */

	      MCMC_One_Rate(tree->n_root,tree->n_root->v[first],YES,tree);
	      MCMC_One_Rate(tree->n_root,tree->n_root->v[secod],YES,tree);
      	    }
      	  else
      	    {
      	      RATES_Posterior_One_Rate(tree->n_root,tree->n_root->v[first],YES,tree);
      	      RATES_Posterior_One_Rate(tree->n_root,tree->n_root->v[secod],YES,tree);
      	    }

	  /* MCMC_Sim_Rate(tree->n_root,tree->n_root->v[2],tree); */
	  /* MCMC_Sim_Rate(tree->n_root,tree->n_root->v[1],tree); */
      	  /* if(tree->mcmc->use_data == YES) Lk(NULL,tree); */
	  /* RATES_Lk_Rates(tree); */          
      	}


      /* printf("\n. move: '%s' lnL: %f",tree->mcmc->move_name[move],tree->rates->c_lnL_times); */
      /* int i; */
      /* for(i = tree -> n_otu; i < 2 * tree -> n_otu -1; i++) printf("\nLOOP Node number:[%d] Lower bound:[%f] Upper bound:[%f] Node time:[%f].", i, */
      /*                                                              tree -> rates -> t_prior_min[i], */
      /*                                                              tree -> rates -> t_prior_max[i], */
      /*                                                              tree -> rates -> nd_t[i]); */

      tree->mcmc->run++;
      MCMC_Get_Acc_Rates(tree->mcmc);

      MCMC_Print_Param(tree->mcmc,tree);
      MCMC_Print_Param_Stdin(tree->mcmc,tree);

      if(tree->io->mutmap == YES)
	{
	  if(!(tree->mcmc->run%tree->mcmc->sample_interval)) 
	    {
	      Sample_Ancestral_Seq(YES,!tree->mcmc->use_data,tree);
	      
	      phydbl sum = 0.0;
	      int edge,site,mut;
	      char *s;
	      FILE *fp;
	      
	      s = (char *)mCalloc(T_MAX_NAME,sizeof(char));
	      
	      strcpy(s,tree->mcmc->io->in_align_file);
	      strcat(s,"_");
	      strcat(s,tree->mcmc->out_filename);
	      strcat(s,".mutmap");
	      fp = fopen(s,"w");
	      
	      Free(s);
	      
	      For(i,(2*tree->n_otu-3)*(tree->n_pattern)*6) sum += tree->mutmap[i];
	      PhyML_Fprintf(fp,"edge\t site\t mut\t count");
	      For(i,(2*tree->n_otu-3)*(tree->n_pattern)*6) 
		{
		  Get_Mutmap_Coord(i,&edge,&site,&mut,tree);
		  PhyML_Fprintf(fp,"\n%4d\t %4d\t %4d\t %10f",edge,site,mut,(phydbl)tree->mutmap[i]/sum);
		}
	      
	      fclose(fp);	      
	    }
	}

      (void)signal(SIGINT,MCMC_Terminate);
    }
  while(tree->mcmc->run < tree->mcmc->chain_len);


}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MCMC_Single_Param_Generic(phydbl *val, 
			       phydbl lim_inf, 
			       phydbl lim_sup, 
			       int move_num,
			       phydbl *lnPrior,
			       phydbl *lnLike,
			       phydbl (*prior_func)(t_edge *,t_tree *,supert_tree *), 
			       phydbl (*like_func)(t_edge *,t_tree *,supert_tree *),
			       int move_type,
			       int _log, /* _log == YES: the model describes the distribution of log(val) but the move applies to val. Need a correction factor */
			       t_edge *branch, t_tree *tree, supert_tree *stree)
{
  phydbl cur_val,new_val,new_lnLike,new_lnPrior,cur_lnLike,cur_lnPrior;
  phydbl u,alpha,ratio;
  phydbl K;
  phydbl new_lnval, cur_lnval;
 
  Record_Br_Len(tree);  

  cur_val       = *val;
  new_val       = -1.0;
  ratio         =  0.0;
  K             = tree->mcmc->tune_move[move_num];
  cur_lnval     = LOG(*val);
  new_lnval     = cur_lnval;

  if(lnLike)
    {
      cur_lnLike = *lnLike;
      new_lnLike = *lnLike;
    }
  else
    {
      cur_lnLike = 0.0;
      new_lnLike = 0.0;
    }
  
  if(lnPrior)
    {
      cur_lnPrior = *lnPrior;
      new_lnPrior = *lnPrior;
    }
  else
    {
      cur_lnPrior = 0.0;
      new_lnPrior = 0.0;
    }

  MCMC_Make_Move(&cur_val,&new_val,lim_inf,lim_sup,&ratio,K,move_type);

  if(new_val < lim_sup && new_val > lim_inf)
    {
      *val = new_val;
  
      if(tree->rates) RATES_Update_Cur_Bl(tree);
            
      if(_log == YES) ratio += (cur_lnval - new_lnval);
      
      if(prior_func) /* Prior ratio */
	{ 
	  new_lnPrior = (*prior_func)(branch,tree,stree); 
	  ratio += (new_lnPrior - cur_lnPrior); 
	}
      
      /* if(like_func && tree->mcmc->use_data == YES)  /\* Likelihood ratio *\/ */
      if(like_func)  /* Likelihood ratio */
	{ 
	  new_lnLike  = (*like_func)(branch,tree,stree);  
	  ratio += (new_lnLike - cur_lnLike);  
	}
      
      /* printf("\n. %s cur_val: %f new_val:%f cur_lnL: %f new_lnL: %f cur_lnPrior: %f new_lnPrior: %f ratio: %f", */
      /*        tree->mcmc->move_name[move_num], */
      /*        cur_val, */
      /*        new_val, */
      /*        cur_lnLike, */
      /*        new_lnLike, */
      /*        cur_lnPrior, */
      /*        new_lnPrior, */
      /*        ratio); */

      ratio = EXP(ratio);
      alpha = MIN(1.,ratio);
      
      if(new_lnLike < UNLIKELY + 0.1) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);

      /* Always accept move */
      if(tree->mcmc->always_yes == YES && new_lnLike > UNLIKELY) alpha = 1.0;

      u = Uni();
      if(u > alpha) /* Reject */
	{
	  *val    = cur_val;
	  new_val = cur_val;
	  if(lnPrior) *lnPrior = cur_lnPrior;
	  if(lnLike)  *lnLike  = cur_lnLike;
	  Restore_Br_Len(tree);
	  if(tree->mod && tree->mod->update_eigen) Update_Eigen(tree->mod);
	}
      else /* Accept */
	{
	  tree->mcmc->acc_move[move_num]++;
	  if(lnPrior) *lnPrior = new_lnPrior;
	  if(lnLike)  *lnLike  = new_lnLike;
	}
    }
      
  tree->mcmc->run_move[move_num]++;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
/* Formula from ``Markov chain Monte Carlo in practice: a roundtable discussion'', 
   Kass, Robert E and Carlin, Bradley P and Gelman, Andrew and Neal, Radford M
   The American Statistician, 1998
*/
void MCMC_Update_Effective_Sample_Size(int move_num, t_mcmc *mcmc, t_tree *tree)
{
  int i,N,lag;
  phydbl rho,mean,var,old_rho,act,ess; 
  int burnin;

  N = mcmc->sample_num+1;
  
  burnin = (int)(0.1*N);
  if(burnin < 1) return;

  N -= burnin;

  mean = Weighted_Mean(mcmc->sampled_val+move_num*mcmc->sample_size+burnin,NULL,N);
  var  = Variance(mcmc->sampled_val+move_num*mcmc->sample_size+burnin,N);
  
  /* if(move_num == tree->mcmc->num_move_phyrex_sigsq)  */
  /*   printf("\n. %d %f %f\n",N,mean,var); */


  act = -1.0;
  old_rho = 1.0;
  For(lag,MIN(N,mcmc->max_lag))
    {      
      rho = 0.0;
      For(i,N-lag) rho += 
        (mcmc->sampled_val[move_num*mcmc->sample_size+burnin+i]     - mean) *
        (mcmc->sampled_val[move_num*mcmc->sample_size+burnin+i+lag] - mean) ;
      
      rho /= (N - lag)*var;

      if(old_rho + rho < 0.0) 
        {
          break; /* Geyer (1992) stopping criterion */
        }

      old_rho = rho;

      act += 2.*rho;
    }

  if(act > 0.0) ess = N/act;
  else          ess = 0.0;

  mcmc->ess[move_num] = ess;

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MCMC_Update_Mode(int move_num, t_mcmc *mcmc, t_tree *tree)
{
  int i,j,N,best_bin;
  int burnin,breaks,*bin_score,best_score;
  phydbl min,max;

  breaks = 100;
  bin_score = (int *)mCalloc(breaks,sizeof(int));

  N = mcmc->sample_num+1;
  
  burnin = (int)(0.1*N);
  if(burnin < 1) return;

  N -= burnin;
  
  min = +INFINITY;
  For(i,N) 
    if(mcmc->sampled_val[move_num*mcmc->sample_size+burnin+i] < min)
      min = mcmc->sampled_val[move_num*mcmc->sample_size+burnin+i];

  max = -INFINITY;
  For(i,N) 
    if(mcmc->sampled_val[move_num*mcmc->sample_size+burnin+i] > max)
      max = mcmc->sampled_val[move_num*mcmc->sample_size+burnin+i];


  For(i,N) 
    {
      for(j=1;j<breaks;j++)
        if(min+j*(max-min)/breaks > mcmc->sampled_val[move_num*mcmc->sample_size+burnin+i])
          {
            bin_score[j-1]++;
            break;
          }
    }

  best_score = 0;
  best_bin = 0;
  For(j,breaks) 
    if(bin_score[j] > best_score) 
      {
        best_score = bin_score[j];
        best_bin = j;
      }
        
  mcmc->mode[move_num] = min + best_bin*(max-min)/breaks;

  Free(bin_score);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MCMC_Clock_R(t_tree *mixt_tree)
{
  t_tree *tree;

  tree = mixt_tree;
  do
    {
      MCMC_Single_Param_Generic(&(tree->rates->clock_r),
                                mixt_tree->rates->min_clock,
                                mixt_tree->rates->max_clock,
                                mixt_tree->mcmc->num_move_clock_r,
                                NULL,&(mixt_tree->c_lnL),
                                NULL,Wrap_Lk,
                                mixt_tree->mcmc->move_type[mixt_tree->mcmc->num_move_clock_r],
                                NO,NULL,mixt_tree,NULL);

      tree = tree->next;
    }
  while(tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

#ifdef GEO
// Sample dispersal parameter from a phylogeo model
void MCMC_GEO_Sigma(t_tree *mixt_tree)
{
  t_tree *tree;

  tree = mixt_tree;
  do
    {
      MCMC_Single_Param_Generic(&(tree->geo->sigma),
                                mixt_tree->geo->min_sigma,
                                mixt_tree->geo->max_sigma,
                                mixt_tree->mcmc->num_move_geo_sigma,
                                NULL,&(mixt_tree->geo->c_lnL),
                                NULL,GEO_Wrap_Lk,
                                mixt_tree->mcmc->move_type[mixt_tree->mcmc->num_move_geo_sigma],
                                NO,NULL,mixt_tree,NULL);

      GEO_Update_Fmat(tree->geo);
      tree = tree->next;
    }
  while(tree);
}
#endif

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

#ifdef GEO
// Sample competition parameter from a phylogeo model
void MCMC_GEO_Lbda(t_tree *mixt_tree)
{
  t_tree *tree;

  tree = mixt_tree;
  do
    {

      MCMC_Single_Param_Generic(&(tree->geo->lbda),
                                mixt_tree->geo->min_lbda,
                                mixt_tree->geo->max_lbda,
                                mixt_tree->mcmc->num_move_geo_lambda,
                                NULL,&(mixt_tree->geo->c_lnL),
                                NULL,GEO_Wrap_Lk,
                                mixt_tree->mcmc->move_type[mixt_tree->mcmc->num_move_geo_lambda],
                                NO,NULL,mixt_tree,NULL);

      tree = tree->next;
    }
  while(tree);
}
#endif

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

#ifdef GEO
// Sample global migration rate parameter from a phylogeo model
void MCMC_GEO_Tau(t_tree *mixt_tree)
{
  t_tree *tree;

  tree = mixt_tree;
  do
    {
      MCMC_Single_Param_Generic(&(tree->geo->tau),
                                mixt_tree->geo->min_tau,
                                mixt_tree->geo->max_tau,
                                mixt_tree->mcmc->num_move_geo_tau,
                                NULL,&(mixt_tree->geo->c_lnL),
                                NULL,GEO_Wrap_Lk,
                                mixt_tree->mcmc->move_type[mixt_tree->mcmc->num_move_geo_tau],
                                NO,NULL,mixt_tree,NULL);

      tree = tree->next;
    }
  while(tree);
}
#endif

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

#ifdef GEO
void MCMC_GEO_Dum(t_tree *mixt_tree)
{
  t_tree *tree;

  tree = mixt_tree;
  do
    {
      MCMC_Single_Param_Generic(&(tree->geo->dum),
                                mixt_tree->geo->min_dum,
                                mixt_tree->geo->max_dum,
                                mixt_tree->mcmc->num_move_geo_dum,
                                NULL,&(mixt_tree->geo->c_lnL),
                                NULL,GEO_Wrap_Lk,
                                mixt_tree->mcmc->move_type[mixt_tree->mcmc->num_move_geo_dum],
                                NO,NULL,mixt_tree,NULL);

      tree = tree->next;
    }
  while(tree);
}
#endif

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

#ifdef GEO
void MCMC_GEO_Loc(t_tree *tree)
{
  int target;
  int *rec_loc; // recorded locations
  int i;
  phydbl cur_lnL, new_lnL;
  phydbl u, ratio, alpha;
  phydbl sum;
  phydbl *probs;

  cur_lnL = tree->geo->c_lnL;

  rec_loc = (int *)mCalloc(2*tree->n_otu-1,sizeof(int));

  For(i,2*tree->n_otu-1) rec_loc[i] = tree->geo->idx_loc[i];

  // Choose an internal node (including the root) at random 
  target = Rand_Int(tree->n_otu,2*tree->n_otu-2); 

  target = 2*tree->n_otu-2;
  
  // Root node is special. Select new location uniformly at random
  if(tree->a_nodes[target] == tree->n_root) 
    {
      probs = (phydbl *)mCalloc(tree->geo->ldscape_sz,sizeof(phydbl));

      sum = 0.0;
      For(i,tree->geo->ldscape_sz) sum += tree->geo->idx_loc_beneath[tree->n_root->num * tree->geo->ldscape_sz + i];
      For(i,tree->geo->ldscape_sz) probs[i] = tree->geo->idx_loc_beneath[tree->n_root->num * tree->geo->ldscape_sz + i]/sum;
      
      tree->geo->idx_loc[tree->n_root->num] = Sample_i_With_Proba_pi(probs,tree->geo->ldscape_sz);      

      Free(probs);
    }


  // Randomize the locations below the selected node
  GEO_Randomize_Locations(tree->a_nodes[target], 
                          tree->geo,
                          tree);
  
  new_lnL = GEO_Lk(tree->geo,tree);

  ratio = (new_lnL - cur_lnL);        
  ratio = EXP(ratio);
  alpha = MIN(1.,ratio);      
  u = Uni();

  if(u > alpha) /* Reject */
    {
      For(i,2*tree->n_otu-1) tree->geo->idx_loc[i] = rec_loc[i];
      tree->geo->c_lnL = GEO_Lk(tree->geo,tree); // TO DO: you only need to update the occupation vector here...
    }

  Free(rec_loc);
}
#endif

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MCMC_Sample_Joint_Rates_Prior(t_tree *tree)
{
  int i,dim;
  phydbl T;
  phydbl *r,*t,*lambda;
  phydbl *min_r,*max_r;
  phydbl k;

  dim    = 2*tree->n_otu-2;
  lambda = tree->rates->_2n_vect1;
  min_r  = tree->rates->_2n_vect2;
  max_r  = tree->rates->_2n_vect3;
  r      = tree->rates->br_r;
  t      = tree->rates->nd_t;

  For(i,dim) tree->rates->mean_r[i] = 1.0;

  RATES_Fill_Lca_Table(tree);
  RATES_Covariance_Mu(tree);

  T = .0;
  For(i,dim) T += (t[tree->a_nodes[i]->num] - t[tree->a_nodes[i]->anc->num]);
  For(i,dim) lambda[i] = (t[tree->a_nodes[i]->num] - t[tree->a_nodes[i]->anc->num])/T;
  For(i,dim) r[i] = 1.0;
  For(i,dim) min_r[i] = tree->rates->min_rate;
  For(i,dim) max_r[i] = tree->rates->max_rate;

  k = 1.; /* We want \sum_i lambda[i] r[i] = 1 */

  Rnorm_Multid_Trunc_Constraint(tree->rates->mean_r, 
				tree->rates->cov_r, 
				min_r,max_r, 
				lambda,
				k, 
				r,
				dim);
}
  
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MCMC_One_Rate(t_node *a, t_node *d, int traversal, t_tree *tree)
{
  t_edge *b;
  int i;
  phydbl u;
  phydbl new_lnL_data, cur_lnL_data, new_lnL_rate, cur_lnL_rate;
  phydbl ratio, alpha;
  phydbl new_mu, cur_mu;
  phydbl r_min, r_max;
  t_edge *b1,*b2,*b3;
  t_node *v2,*v3;
  int move_num;
  phydbl K;

  if(tree->rates->model == STRICTCLOCK) return;

  b = NULL;
  if(a == tree->n_root) b = tree->e_root;
  else
    For(i,3) if(d->v[i] == a) { b = d->b[i]; break; }
   
  cur_mu       = tree->rates->br_r[d->num];
  cur_lnL_data = tree->c_lnL;
  new_lnL_data = tree->c_lnL;
  cur_lnL_rate = tree->rates->c_lnL_rates;
  new_lnL_rate = tree->rates->c_lnL_rates;
  r_min        = tree->rates->min_rate;
  r_max        = tree->rates->max_rate;
  ratio        = 0.0;
  move_num     = d->num+tree->mcmc->num_move_br_r;
  K            = tree->mcmc->tune_move[move_num];

  Record_Br_Len(tree);
  
  u = Uni();
  
  MCMC_Make_Move(&cur_mu,&new_mu,r_min,r_max,&ratio,K,tree->mcmc->move_type[tree->mcmc->num_move_br_r+d->num]);

  /* phydbl dt,sd,mean; */
  /* int err; */
  /* dt     = tree->rates->nd_t[d->num] - tree->rates->nd_t[a->num]; */
  /* sd     = SQRT(tree->rates->nu * dt); */
  /* mean   = tree->rates->br_r[a->num]; */
  /* new_mu = Rnorm_Trunc(mean,sd,r_min,r_max,&err); */
  /* ratio  = Log_Dnorm_Trunc(cur_mu,mean,sd,r_min,r_max,&err) - Log_Dnorm_Trunc(new_mu,mean,sd,r_min,r_max,&err); */

  if(new_mu > r_min && new_mu < r_max)
    {      
      tree->rates->br_r[d->num] = new_mu;
            
      v2 = v3 = NULL;
      For(i,3)
	if((d->v[i] != a) && (d->b[i] != tree->e_root))
	  {
	    if(!v2) { v2 = d->v[i]; }
	    else    { v3 = d->v[i]; }
	  }
      
      
      b1 = NULL;
      if(a == tree->n_root) b1 = tree->e_root;
      else For(i,3) if(d->v[i] == a) { b1 = d->b[i]; break; }
      
      b2 = b3 = NULL;
      if(!d->tax)
	{
	  For(i,3)
	    if((d->v[i] != a) && (d->b[i] != tree->e_root))
	      {
		if(!b2) { b2 = d->b[i]; }
		else    { b3 = d->b[i]; }
	      }
	}
      
      tree->rates->br_do_updt[d->num] = YES;
      if(!d->tax)
      	{
      	  tree->rates->br_do_updt[v2->num] = YES;
      	  tree->rates->br_do_updt[v3->num] = YES;
      	}
      
      RATES_Update_Cur_Bl(tree);
      
      /* printf("\n. r0=%f r1=%f cr=%f mean=%f var=%f nu=%f dt=%f", */
      /* 	 r0,r1,tree->rates->clock_r,b1->gamma_prior_mean,b1->gamma_prior_var,nu,t1-t0); */
            
      if(tree->mcmc->use_data)
	{
	  if(tree->io->lk_approx == EXACT)
	    {
	      Update_PMat_At_Given_Edge(b1,tree);
	      if(!d->tax)
		{
		  Update_PMat_At_Given_Edge(b2,tree);
		  Update_PMat_At_Given_Edge(b3,tree);
		}
	      Update_P_Lk(tree,b1,d);
	    }
	  new_lnL_data = Lk(b1,tree);
	  
	  /* tree->both_sides = NO; */
	  /* new_lnL_data = Lk(tree); */
	}
      
      new_lnL_rate = RATES_Lk_Rates(tree);
      
      /* Likelihood ratio */
      if(tree->mcmc->use_data) ratio += (new_lnL_data - cur_lnL_data);

      /* Prior ratio */
      ratio += (new_lnL_rate - cur_lnL_rate);
            
      ratio = EXP(ratio);
      alpha = MIN(1.,ratio);
      
      u = Uni();
      
      if(u > alpha) /* Reject */
	{
	  tree->rates->br_r[d->num] = cur_mu;	  
	  tree->c_lnL               = cur_lnL_data;
	  tree->rates->c_lnL_rates  = cur_lnL_rate;
	  
	  Restore_Br_Len(tree);
	  
	  if(tree->mcmc->use_data && tree->io->lk_approx == EXACT)
	    {
	      Update_PMat_At_Given_Edge(b1,tree);
	      if(!d->tax)
		{
		  Update_PMat_At_Given_Edge(b2,tree);
		  Update_PMat_At_Given_Edge(b3,tree);
		}
	      Update_P_Lk(tree,b1,d);
	    }
	  
	  /* tree->both_sides = YES; */
	  /* new_lnL_data = Lk(tree); */
	  /* tree->both_sides = NO; */
	  
	}
      else
	{
	  tree->mcmc->acc_move[tree->mcmc->num_move_br_r+d->num]++;
	}
    }
  tree->mcmc->run_move[tree->mcmc->num_move_br_r+d->num]++;
  
  if(traversal == YES)
    {
      if(d->tax == YES) return;
      else
	{
	  For(i,3)
	    if(d->v[i] != a && d->b[i] != tree->e_root)
	      {
		if(tree->io->lk_approx == EXACT && tree->mcmc->use_data) Update_P_Lk(tree,d->b[i],d);
		/* if(tree->io->lk_approx == EXACT && tree->mcmc->use_data) {tree->both_sides = YES; Lk(tree); } */
		MCMC_One_Rate(d,d->v[i],YES,tree);
	      }
	}
      if(tree->io->lk_approx == EXACT && tree->mcmc->use_data) Update_P_Lk(tree,b,d);
      /* if(tree->io->lk_approx == EXACT && tree->mcmc->use_data) {tree->both_sides = YES; Lk(tree); } */
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MCMC_One_Node_Rate(t_node *a, t_node *d, int traversal, t_tree *tree)
{
  t_edge *b;
  int i;

  b = NULL;
  if(a == tree->n_root) b = tree->e_root;
  else
    For(i,3) if(d->v[i] == a) { b = d->b[i]; break; }

  /* Only the LOG_RANDWALK move seems to work here. Change with caution then. */
  tree->rates->br_do_updt[d->num] = YES;
  MCMC_Single_Param_Generic(&(tree->rates->nd_r[d->num]),
			    tree->rates->min_rate,
			    tree->rates->max_rate,
			    tree->mcmc->num_move_nd_r+d->num,
			    &(tree->rates->c_lnL_rates),NULL,
			    Wrap_Lk_Rates,NULL,
			    tree->mcmc->move_type[tree->mcmc->num_move_nd_r+d->num],
			    NO,NULL,tree,NULL);


  Update_PMat_At_Given_Edge(b,tree);

  if(traversal == YES)
    {
      if(d->tax == YES) return;
      else
	{
	  For(i,3)
	    if(d->v[i] != a && d->b[i] != tree->e_root)
	      {
		MCMC_One_Node_Rate(d,d->v[i],YES,tree);
	      }
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MCMC_One_Time(t_node *a, t_node *d, int traversal, t_tree *tree)
{
  phydbl u;
  phydbl t_min,t_max;
  phydbl t1_cur, t1_new;
  phydbl cur_lnL_data, new_lnL_data;
  phydbl cur_lnL_rate, new_lnL_rate;
  phydbl cur_lnL_time, new_lnL_time;
  phydbl ratio,alpha;
  t_edge *b1,*b2,*b3;
  int    i;
  phydbl t0,t2,t3;
  t_node *v2,*v3;
  phydbl K;
  int move_num;
  
  if(d->tax) return; /* Won't change time at tip */

  /* if(FABS(tree->rates->t_prior_min[d->num] - tree->rates->t_prior_max[d->num]) < 1.E-10) return; */

  Record_Br_Len(tree);
  RATES_Record_Rates(tree);
  RATES_Record_Times(tree);
  
  move_num       = d->num-tree->n_otu+tree->mcmc->num_move_nd_t;
  K              = tree->mcmc->tune_move[move_num];
  cur_lnL_data   = tree->c_lnL;
  cur_lnL_rate   = tree->rates->c_lnL_rates;
  t1_cur         = tree->rates->nd_t[d->num];
  new_lnL_data   = cur_lnL_data;
  new_lnL_rate   = cur_lnL_rate;
  ratio          = 0.0;
  cur_lnL_time   = tree->rates->c_lnL_times;
  new_lnL_time   = cur_lnL_time;

  v2 = v3 = NULL;
  For(i,3)
    if((d->v[i] != a) && (d->b[i] != tree->e_root))
      {
	if(!v2) { v2 = d->v[i]; }
	else    { v3 = d->v[i]; }
      }


  b1 = NULL;
  if(a == tree->n_root) b1 = tree->e_root;
  else For(i,3) if(d->v[i] == a) { b1 = d->b[i]; break; }

  b2 = b3 = NULL;
  For(i,3)
    if((d->v[i] != a) && (d->b[i] != tree->e_root))
      {
	if(!b2) { b2 = d->b[i]; }
	else    { b3 = d->b[i]; }
      }
 
  t0 = tree->rates->nd_t[a->num];
  t2 = tree->rates->nd_t[v2->num];
  t3 = tree->rates->nd_t[v3->num];


  /* t_min = MAX(t0,tree->rates->t_prior_min[d->num]);        */
  /* t_max = MIN(MIN(t2,t3),tree->rates->t_prior_max[d->num]);*/

  t_min = t0;
  t_max = MIN(t2,t3);

  t_min += tree->rates->min_dt;
  t_max -= tree->rates->min_dt;

  if(t_min > t_max) 
    {
      PhyML_Printf("\n== t_min = %f t_max = %f",t_min,t_max);
      PhyML_Printf("\n== prior_min = %f prior_max = %f",tree->rates->t_prior_min[d->num],tree->rates->t_prior_max[d->num]);
      PhyML_Printf("\n== Err. in file %s at line %d\n",__FILE__,__LINE__);
      /* Exit("\n"); */
    }
 
  MCMC_Make_Move(&t1_cur,&t1_new,t_min,t_max,&ratio,K,tree->mcmc->move_type[move_num]);
  

  if(t1_new > t_min && t1_new < t_max) 
    {
      tree->rates->nd_t[d->num] = t1_new;

      new_lnL_time = TIMES_Lk_Times(tree);
      ratio += (new_lnL_time - cur_lnL_time);

      if(isinf(new_lnL_time) == NO) // Proposed value of t is inside its boundary
        {
          /* Update branch lengths */
          tree->rates->br_do_updt[d->num]  = YES;
          tree->rates->br_do_updt[v2->num] = YES;
          tree->rates->br_do_updt[v3->num] = YES;
          RATES_Update_Cur_Bl(tree);
          
          new_lnL_rate = RATES_Lk_Rates(tree);
          ratio += (new_lnL_rate - cur_lnL_rate);

          if(tree->mcmc->use_data)
            {
              if(tree->io->lk_approx == EXACT)
                {
                  Update_PMat_At_Given_Edge(b1,tree);
                  Update_PMat_At_Given_Edge(b2,tree);
                  Update_PMat_At_Given_Edge(b3,tree);
                  Update_P_Lk(tree,b1,d);
                }
              new_lnL_data = Lk(b1,tree);
              
              /* /\* !!!!!!!!!!!!!!!!!!!!1 *\/ */
              /* if(FABS(Lk(tree) - new_lnL_data) > 1.E-5) */
              /*   { */
              /*     PhyML_Printf("\n. b1->l->v=%f b2->l->v=%f b3->l->v=%f",b1->l->v,b2->l->v,b3->l->v); */
              /*     PhyML_Printf("\n. a->num=%d d->num=%d root=%d (%d %d)",a->num,d->num,a==tree->n_root,tree->n_root->v[2]->num,tree->n_root->v[1]->num); */
              /*     PhyML_Printf("\n. t_min=%f t_max=%f",t_min,t_max); */
              /*     PhyML_Printf("\n. t1_new=%f t1_cur=%f",t1_new,t1_cur); */
              /*     PhyML_Printf("\n. %f %f",cur_lnL_data,tree->c_lnL); */
              /*     PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__); */
              /*     Warn_And_Exit(""); */
              /*   } */
            }
          
          if(tree->mcmc->use_data) ratio += (new_lnL_data - cur_lnL_data);
        }
      
      /* if(d->num == 7) */
      /*   { */
      /*     printf("\n. nd_t: %f %f new_lnL_rate: %f cur_lnL_rate: %f new_lnL_time: %f cur_lnL_time: %f ratio: %f", */
      /*            t1_cur, */
      /*            tree->rates->nd_t[d->num], */
      /*            new_lnL_rate, */
      /*            cur_lnL_rate, */
      /*            new_lnL_time, */
      /*            cur_lnL_time, */
      /*            ratio); */
      /*   } */

          
      ratio = EXP(ratio);
      alpha = MIN(1.,ratio);
      u = Uni();
	           
      if(u > alpha) /* Reject */
	{
          //if(d -> num == 7) PhyML_Printf("\n. t_cur = %f t_rej = %f", t1_cur, t1_new);
	  tree->rates->nd_t[d->num] = t1_cur;
	  tree->c_lnL              = cur_lnL_data;
	  tree->rates->c_lnL_rates = cur_lnL_rate;
	  tree->rates->c_lnL_times = cur_lnL_time;
          
          if(isinf(new_lnL_time) == NO)
            {
              Restore_Br_Len(tree);
              RATES_Reset_Rates(tree);
              RATES_Reset_Times(tree);
              
              if(tree->io->lk_approx == EXACT && tree->mcmc->use_data) 
                {
                  Update_PMat_At_Given_Edge(b1,tree);
                  Update_PMat_At_Given_Edge(b2,tree);
                  Update_PMat_At_Given_Edge(b3,tree);
                  Update_P_Lk(tree,b1,d);
                }
            }
	}
      else
	{
	  /* printf("\n. A new_lnL_data = %f cur_lnL_data = %f t_new=%f t_cur=%f tmin=%f tmax=%f", */
	  /* 	 new_lnL_data,cur_lnL_data,t1_new,t1_cur,t_min,t_max); */
	  tree->mcmc->acc_move[move_num]++;
	}
      if(t1_new < t0)
	{
	  t1_new = t0+1.E-4;
	  PhyML_Printf("\n");
	  PhyML_Printf("\n== a is root -> %s",(a == tree->n_root)?("YES"):("NO"));
	  PhyML_Printf("\n== t0 = %f t1_new = %f",t0,t1_new);
	  PhyML_Printf("\n== t_min=%f t_max=%f",t_min,t_max);
	  PhyML_Printf("\n== (t1-t0)=%f (t2-t1)=%f",t1_cur-t0,t2-t1_cur);
	  PhyML_Printf("\n== Err. in file %s at line %d\n",__FILE__,__LINE__);
	  /*       Exit("\n"); */
	}
      if(t1_new > MIN(t2,t3))
	{
	  PhyML_Printf("\n");
	  PhyML_Printf("\n== a is root -> %s",(a == tree->n_root)?("YES"):("NO"));
	  PhyML_Printf("\n== t0 = %f t1_new = %f t1 = %f t2 = %f t3 = %f MIN(t2,t3)=%f",t0,t1_new,t1_cur,t2,t3,MIN(t2,t3));
	  PhyML_Printf("\n== t_min=%f t_max=%f",t_min,t_max);
	  PhyML_Printf("\n== (t1-t0)=%f (t2-t1)=%f",t1_cur-t0,t2-t1_cur);
	  PhyML_Printf("\n== Err. in file %s at line %d\n",__FILE__,__LINE__);
	  /*       Exit("\n"); */
	}
      
      if(isnan(t1_new))
	{
	  PhyML_Printf("\n== run=%d",tree->mcmc->run);
	  PhyML_Printf("\n== Err in file %s at line %d\n",__FILE__,__LINE__);
	  /*       Exit("\n"); */
	}
    }
  
  tree->mcmc->run_move[move_num]++;

  if(traversal == YES)
    {
      if(d->tax == YES) return;
      else
	For(i,3)
	  if(d->v[i] != a && d->b[i] != tree->e_root)
	    {
	      if(tree->io->lk_approx == EXACT && tree->mcmc->use_data) Update_P_Lk(tree,d->b[i],d);
	      /* if(tree->io->lk_approx == EXACT && tree->mcmc->use_data) {tree->both_sides = YES; Lk(tree); } */
	      MCMC_One_Time(d,d->v[i],YES,tree);
	    }
      if(tree->io->lk_approx == EXACT && tree->mcmc->use_data) Update_P_Lk(tree,b1,d);
      /* if(tree->io->lk_approx == EXACT && tree->mcmc->use_data) {tree->both_sides = YES; Lk(tree); } */
    }	    
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MCMC_Jump_Calibration(t_tree *tree)
{
#ifdef INVITEE
  phydbl u;
  phydbl *calib_prior_cumprob, *calib_prior_prob, *uniform_prob; 
  phydbl cur_lnL_data, new_lnL_data;
  phydbl cur_lnL_rate, new_lnL_rate;
  phydbl cur_lnL_time, new_lnL_time;
  phydbl cur_lnL_K;
  phydbl times_log_hastings_ratio;
  int new_calib_comb_num, cur_calib_comb_num;
  phydbl ratio, alpha;
  int i, result;
  int move_num;
  int tot_num_of_calib_comb;
  int *buff_calib_num;
  int n_pos_probs;

  tot_num_of_calib_comb = Number_Of_Comb(tree -> rates -> calib);

  tree->rates->update_time_norm_const = YES;

  if(tot_num_of_calib_comb > 1)
    {
      calib_prior_prob = tree -> rates -> times_partial_proba;
      calib_prior_cumprob = (phydbl *)mCalloc(tot_num_of_calib_comb + 1, sizeof(phydbl)); 
      uniform_prob =  (phydbl *)mCalloc(tot_num_of_calib_comb + 1, sizeof(phydbl));

      ///////////////////////////////////////////////////////////////////////////////////////////
      /* for(i = tree -> n_otu; i < 2 * tree -> n_otu - 1; i++) printf("\n. '%f' '%f' \n", tree -> rates -> t_prior_min[i], tree -> rates -> t_prior_max[i]); */
      //for(i = tree -> n_otu; i < 2 * tree -> n_otu - 1; i++) printf("\n. '%f' \n", tree -> rates -> nd_t[i]);
      ///////////////////////////////////////////////////////////////////////////////////////////
      Record_Br_Len(tree);
      RATES_Record_Rates(tree);
      RATES_Record_Times(tree);
      TIMES_Record_Prior_Times(tree);
         
      move_num = tree -> mcmc -> num_move_jump_calibration;

      cur_lnL_data = tree -> c_lnL;
      cur_lnL_rate = tree -> rates -> c_lnL_rates;

      new_lnL_data = cur_lnL_data;
      new_lnL_rate = cur_lnL_rate;

      ratio = 0.0;
      cur_lnL_time = tree -> rates -> c_lnL_times;
      new_lnL_time = cur_lnL_time;

      cur_calib_comb_num = tree -> rates -> cur_comb_numb;
      new_calib_comb_num = cur_calib_comb_num;

      cur_lnL_K = tree->rates->log_K_cur;

      /* new_calib_comb_num = Sample_i_With_Proba_pi(calib_prior_prob,tot_num_of_calib_comb); */

      u = Uni();
      if(u < 0.5)
        {
          For(i,tot_num_of_calib_comb) uniform_prob[i] = calib_prior_prob[i] > .0 ? 1. : 0.;
          new_calib_comb_num = Sample_i_With_Proba_pi(uniform_prob,tot_num_of_calib_comb);
        }
      else
        {
          buff_calib_num = (int *)mCalloc(tot_num_of_calib_comb,sizeof(int));
          
          n_pos_probs = 0;
          For(i,tot_num_of_calib_comb) 
            if(calib_prior_prob[i] > 1.E-10)
              {
                buff_calib_num[n_pos_probs] = i;
                n_pos_probs++;
              }
          
          int curr_pos=0;
          while(buff_calib_num[curr_pos] != cur_calib_comb_num) { curr_pos++; }
          u = Uni();
          if(u < 0.5) new_calib_comb_num = buff_calib_num[Modulo(curr_pos-1,n_pos_probs)];
          else        new_calib_comb_num = buff_calib_num[Modulo(curr_pos+1,n_pos_probs)];
          
          Free(buff_calib_num);
        }

      /* printf("\n"); */
      /* printf("\n. current calibration: %d ", cur_calib_comb_num); */
      /* printf("\n"); */
      /* for(i = tree -> n_otu; i < 2 * tree -> n_otu -1; i++) printf("\nJUMP0 Node number:[%d] Lower bound:[%f] Upper bound:[%f] Node time:[%f].", i, */
      /*                                                              tree -> rates -> t_prior_min[i], */
      /*                                                              tree -> rates -> t_prior_max[i], */
      /*                                                              tree -> rates -> nd_t[i]); */

      /* for(i = tree -> n_otu; i < 2 * tree -> n_otu -1; i++) printf("\n. Node number:[%d] Lower bound:[%f] Upper bound:[%f] Node time:[%f]. \n", i, tree -> rates -> t_prior_min[i], tree -> rates -> t_prior_max[i], tree -> rates -> nd_t[i]); */ /* Exit("\n"); */
      if(!Are_Equal(cur_calib_comb_num, new_calib_comb_num, 1.E-10))
        {
          /* printf("\n"); */
          /* printf("\n. current calibration: %d == new_calibration: %d", cur_calib_comb_num, new_calib_comb_num); */
          /* printf("\n"); */
          /* Print_Node(tree->n_root,tree->n_root->v[1],tree); */
          /* Print_Node(tree->n_root,tree->n_root->v[2],tree); */
          /* for(i = tree -> n_otu; i < 2 * tree -> n_otu -1; i++) printf("\nAV Node number:[%d] Lower bound:[%f] Upper bound:[%f] Node time:[%f].", i, */
          /*                                                              tree -> rates -> t_prior_min[i], */
          /*                                                              tree -> rates -> t_prior_max[i], */
          /*                                                              tree -> rates -> nd_t[i]); */

          Set_Current_Calibration(new_calib_comb_num, tree);
          TIMES_Set_All_Node_Priors(tree);
          
          result = TRUE;
          
          /* Check_Node_Time(tree -> n_root, tree -> n_root -> v[1], &result, tree); */
          /* Check_Node_Time(tree -> n_root, tree -> n_root -> v[2], &result, tree); */
          
          /* printf("\n. New calibration %d", new_calib_comb_num); */
          /* printf("\n"); */
          /* Print_Node(tree->n_root,tree->n_root->v[1],tree); */
          /* Print_Node(tree->n_root,tree->n_root->v[2],tree); */
          /* for(i = tree -> n_otu; i < 2 * tree -> n_otu -1; i++) printf("\nAP Node number:[%d] Lower bound:[%f] Upper bound:[%f] Node time:[%f].", i, */
          /*                                                              tree -> rates -> t_prior_min[i], */
          /*                                                              tree -> rates -> t_prior_max[i], */
          /*                                                              tree -> rates -> nd_t[i]); */


          times_log_hastings_ratio = 0.0;  
          Jump_Calibration_Move_Pre(tree->n_root,tree->n_root->v[1],tree->rates->nd_t[tree->n_root->num],&times_log_hastings_ratio,tree);
          Jump_Calibration_Move_Pre(tree->n_root,tree->n_root->v[2],tree->rates->nd_t[tree->n_root->num],&times_log_hastings_ratio,tree);

          /* Update_Current_Times_Down_Tree(tree -> n_root, tree -> n_root -> v[1], tree); */
          /* Update_Current_Times_Down_Tree(tree -> n_root, tree -> n_root -> v[2], tree); */
          /* new_lnL_proposal_density = 0.0; */
          /* Multiple_Time_Proposal_Density(tree -> n_root, tree -> n_root -> v[1], &new_lnL_proposal_density, tree); */
          /* Multiple_Time_Proposal_Density(tree -> n_root, tree -> n_root -> v[2], &new_lnL_proposal_density, tree); */
          
          result = TRUE;
          
          Check_Node_Time(tree -> n_root, tree -> n_root -> v[1], &result, tree);
          Check_Node_Time(tree -> n_root, tree -> n_root -> v[2], &result, tree);
          
         
          if(result != TRUE)
            {
              PhyML_Printf("\n. ...................... OLD CALIBRATION.....................................\n");
              for(i = tree -> n_otu; i < 2 * tree -> n_otu -1; i++) printf("\n. Node number:[%d] Lower bound:[%f] Upper bound:[%f] Node time:[%f]. \n", i, tree -> rates -> t_prior_min_ori[i], tree -> rates -> t_prior_max_ori[i], tree -> rates -> buff_t[i]);
              PhyML_Printf("\n. ...........................................................................\n");
              PhyML_Printf("\n. ................. NEW PROPOSED CALIBRATION ................................\n");
              for(i = tree -> n_otu; i < 2 * tree -> n_otu -1; i++) printf("\n. Node number:[%d] Lower bound:[%f] Upper bound:[%f] Node time:[%f]. \n", i, tree -> rates -> t_prior_min[i], tree -> rates -> t_prior_max[i], tree -> rates -> nd_t[i]);
              PhyML_Printf("\n. ...........................................................................\n");
              PhyML_Printf("\n== There is a problem with calibration information.\n");
              Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
            }
          
          For(i,2*tree->n_otu-2) tree->rates->br_do_updt[i] = YES;
          RATES_Update_Cur_Bl(tree);
          if(tree->mcmc->use_data) new_lnL_data = Lk(NULL,tree);
          
          new_lnL_rate = RATES_Lk_Rates(tree);
          new_lnL_time = TIMES_Lk_Times(tree);
          
          /* printf("\n. JUMP cur_lnL_time: %f new_lnL_time: %f",cur_lnL_time,new_lnL_time); */
          /* for(i = tree -> n_otu; i < 2 * tree -> n_otu -1; i++) printf("\nJUMP1 Node number:[%d] Lower bound:[%f] Upper bound:[%f] Node time:[%f].", i, */
          /*                                                          tree -> rates -> t_prior_min[i], */
          /*                                                          tree -> rates -> t_prior_max[i], */
          /*                                                          tree -> rates -> nd_t[i]); */

          /* Likelihood ratio */
          if(tree->mcmc->use_data) ratio += (new_lnL_data - cur_lnL_data);
          
          /* Prior ratio */
          ratio += (new_lnL_rate - cur_lnL_rate);
          ratio += (new_lnL_time - cur_lnL_time);
          ratio += (LOG(calib_prior_prob[new_calib_comb_num]) - LOG(calib_prior_prob[cur_calib_comb_num]));
          ratio += times_log_hastings_ratio;
          
          ratio = EXP(ratio);
          alpha = MIN(1.,ratio);
          u = Uni();
          
          if(u > alpha)
            {
              RATES_Reset_Times(tree);
              Restore_Br_Len(tree);
              TIMES_Reset_Prior_Times(tree);
              tree -> c_lnL                         = cur_lnL_data;
              tree -> rates -> c_lnL_rates          = cur_lnL_rate;
              tree -> rates -> c_lnL_times          = cur_lnL_time;
              tree -> rates -> log_K_cur            = cur_lnL_K;
              tree -> rates -> cur_comb_numb        = cur_calib_comb_num;
              /* for(i = tree -> n_otu; i < 2 * tree -> n_otu -1; i++) printf("\nJUMP3 Node number:[%d] Lower bound:[%f,%f] Upper bound:[%f,%f] Node time:[%f].", i, */
              /*                                                              tree->rates->t_prior_min[i],tree->rates->t_prior_min_ori[i], */
              /*                                                              tree->rates->t_prior_max[i],tree->rates->t_prior_max_ori[i], */
              /*                                                              tree->rates->nd_t[i]); */
              Set_Current_Calibration(tree -> rates -> cur_comb_numb, tree);
              /* for(i = tree -> n_otu; i < 2 * tree -> n_otu -1; i++) printf("\nJUMP4 Node number:[%d] Lower bound:[%f,%f] Upper bound:[%f,%f] Node time:[%f].", i, */
              /*                                                              tree->rates->t_prior_min[i],tree->rates->t_prior_min_ori[i], */
              /*                                                              tree->rates->t_prior_max[i],tree->rates->t_prior_max_ori[i], */
              /*                                                              tree->rates->nd_t[i]); */
              TIMES_Set_All_Node_Priors(tree);
              /* printf("\n. ......................REJECTED.....................................\n"); */
              tree -> rates -> numb_calib_chosen[cur_calib_comb_num]++;
            }
          else
            {
              tree -> rates -> cur_comb_numb = new_calib_comb_num;
              tree->mcmc->acc_move[move_num]++;
              /* printf("\n. ......................ACCEPTED.....................................\n"); */
              tree -> rates -> numb_calib_chosen[new_calib_comb_num]++;
            }      

          /* for(i = tree -> n_otu; i < 2 * tree -> n_otu -1; i++) printf("\nJUMP2 Node number:[%d] Lower bound:[%f,%f] Upper bound:[%f,%f] Node time:[%f].", i, */
          /*                                                              tree->rates->t_prior_min[i],tree->rates->t_prior_min_ori[i], */
          /*                                                              tree->rates->t_prior_max[i],tree->rates->t_prior_max_ori[i], */
          /*                                                              tree->rates->nd_t[i]); */

          Check_Node_Time(tree -> n_root, tree -> n_root -> v[1], &result, tree);
          Check_Node_Time(tree -> n_root, tree -> n_root -> v[2], &result, tree);
          
         
          if(result != TRUE) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);

          tree->mcmc->run_move[move_num]++;
          Free(calib_prior_cumprob);
          Free(uniform_prob);
        }                                             
    }
  
  tree->rates->update_time_norm_const = NO;

#endif
}
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////

void MCMC_Root_Time(t_tree *tree)
{
  phydbl u;
  phydbl t_min,t_max;
  phydbl t1_cur, t1_new;
  phydbl cur_lnL_data, new_lnL_data;
  phydbl cur_lnL_rate, new_lnL_rate;
  phydbl cur_lnL_time, new_lnL_time;
  phydbl ratio,alpha;
  t_edge *b1;
  phydbl t0,t2,t3;
  t_node *v2,*v3;
  phydbl K;
  int move_num;
  t_node *root;

  root = tree->n_root;

  if(FABS(tree->rates->t_prior_min[root->num] - tree->rates->t_prior_max[root->num]) < 1.E-10) return;

  Record_Br_Len(tree);
  RATES_Record_Rates(tree);
  RATES_Record_Times(tree);
  
  move_num       = root->num-tree->n_otu+tree->mcmc->num_move_nd_t;
  K              = tree->mcmc->tune_move[move_num];
  cur_lnL_data   = tree->c_lnL;
  cur_lnL_rate   = tree->rates->c_lnL_rates;
  t1_cur         = tree->rates->nd_t[root->num];
  new_lnL_data   = cur_lnL_data;
  new_lnL_rate   = cur_lnL_rate;
  ratio          = 0.0;
  cur_lnL_time   = tree->rates->c_lnL_times;
  new_lnL_time   = cur_lnL_time;


  v2 = root->v[2];
  v3 = root->v[1];

  b1 = tree->e_root;
  
  t0 = tree->rates->t_prior_min[root->num];
  t2 = tree->rates->nd_t[v2->num];
  t3 = tree->rates->nd_t[v3->num];
  /* printf("\n %d t0= %f %d t2= %f %d t3= %f \n", root -> num, t0, v2 -> num, t2, v3 -> num, t3);  */

  t_min = t0;
  /* t_max = MIN(MIN(t2,t3),tree->rates->t_prior_max[root->num]); */
  t_max = MIN(t2,t3);

  t_min += tree->rates->min_dt;
  t_max -= tree->rates->min_dt;

  if(t_min > t_max) 
    {
      PhyML_Printf("\n== t0 = %f t2 = %f t3 = %f",t0,t2,t3);
      PhyML_Printf("\n== t_min = %f t_max = %f",t_min,t_max);
      PhyML_Printf("\n== prior_min = %f prior_max = %f",tree->rates->t_prior_min[root->num],tree->rates->t_prior_max[root->num]);
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
    }

  MCMC_Make_Move(&t1_cur,&t1_new,t_min,t_max,&ratio,K,tree->mcmc->move_type[move_num]);
 
  if(t1_new > t_min && t1_new < t_max) 
    {
      tree->rates->nd_t[root->num] = t1_new;

      /* Update branch lengths */
      RATES_Update_Cur_Bl(tree);

      if(tree->mcmc->use_data) new_lnL_data = Lk(b1,tree);

      new_lnL_rate = RATES_Lk_Rates(tree);
      new_lnL_time = TIMES_Lk_Times(tree); 

      if(tree->mcmc->use_data) ratio += (new_lnL_data - cur_lnL_data);
      ratio += (new_lnL_rate - cur_lnL_rate);
      ratio += (new_lnL_time - cur_lnL_time);

      ratio = EXP(ratio);
      alpha = MIN(1.,ratio);
      u = Uni();
         
      if(u > alpha) /* Reject */
	{
	  tree->rates->nd_t[root->num] = t1_cur;
	  tree->c_lnL              = cur_lnL_data;
	  tree->rates->c_lnL_rates = cur_lnL_rate;
	  tree->rates->c_lnL_times = cur_lnL_time;
	  Restore_Br_Len(tree);
	  RATES_Reset_Rates(tree);
	  RATES_Reset_Times(tree);
	}
      else
	{
	  tree->mcmc->acc_move[move_num]++;
	}
      
      if(t1_new < t0)
	{
	  t1_new = t0+1.E-4;
	  PhyML_Printf("\n");
	  PhyML_Printf("\n== t0 = %f t1_new = %f",t0,t1_new);
	  PhyML_Printf("\n== t_min=%f t_max=%f",t_min,t_max);
	  PhyML_Printf("\n== (t1-t0)=%f (t2-t1)=%f",t1_cur-t0,t2-t1_cur);
	  PhyML_Printf("\n== Err. in file %s at line %d\n",__FILE__,__LINE__);
	  /*       Exit("\n"); */
	}
      if(t1_new > MIN(t2,t3))
	{
	  PhyML_Printf("\n");
	  PhyML_Printf("\n== t0 = %f t1_new = %f t1 = %f t2 = %f t3 = %f MIN(t2,t3)=%f",t0,t1_new,t1_cur,t2,t3,MIN(t2,t3));
	  PhyML_Printf("\n== t_min=%f t_max=%f",t_min,t_max);
	  PhyML_Printf("\n== (t1-t0)=%f (t2-t1)=%f",t1_cur-t0,t2-t1_cur);
	  PhyML_Printf("\n== Err. in file %s at line %d\n",__FILE__,__LINE__);
	  /*       Exit("\n"); */
	}
      
      if(isnan(t1_new))
	{
	  PhyML_Printf("\n== run=%d",tree->mcmc->run);
	  PhyML_Printf("\n== Err/ in file %s at line %d\n",__FILE__,__LINE__);
	  /*       Exit("\n"); */
	}
    }
  
  tree->mcmc->run_move[move_num]++;

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MCMC_Tree_Height(t_tree *tree)
{
  int i;
  phydbl K,mult,u,alpha,ratio;
  phydbl cur_lnL_data,new_lnL_data;
  phydbl cur_lnL_rate,new_lnL_rate;
  phydbl cur_lnL_time,new_lnL_time;
  phydbl floor;
  int n_nodes;
  
  if(FABS(tree->rates->t_prior_max[tree->n_root->num] - tree->rates->t_prior_min[tree->n_root->num]) < 1.E-10) return;

  RATES_Record_Times(tree);
  Record_Br_Len(tree);

  K            = tree->mcmc->tune_move[tree->mcmc->num_move_tree_height];
  cur_lnL_data = tree->c_lnL;
  new_lnL_data = tree->c_lnL;
  ratio        = 0.0;
  cur_lnL_rate = tree->rates->c_lnL_rates;
  new_lnL_rate = tree->rates->c_lnL_rates;
  cur_lnL_time = tree->rates->c_lnL_times;
  
  u = Uni();
  mult = EXP(K*(u-0.5));
 
 /* WARNING: It must not be floor = tree->rates->t_prior_max[tree->n_root->num]; 
     floor is the maximum value a node height can take when one ignores the 
     calibration nodes, i.e., floor is set by the height of the tips
  */
  /* floor = 0.0; */
  floor = tree->rates->t_floor[tree->n_root->num];
  /* floor = tree->rates->t_prior_max[tree->n_root->num]; */

  Scale_Subtree_Height(tree->n_root,mult,floor,&n_nodes,tree);

  /* For(i,2*tree->n_otu-1) */
  /*   { */
  /*     if(tree->rates->nd_t[i] > tree->rates->t_prior_max[i] || */
  /* 	 tree->rates->nd_t[i] < tree->rates->t_prior_min[i]) */
  /* 	{ */
  /* 	  RATES_Reset_Times(tree); */
  /* 	  Restore_Br_Len(tree); */
  /* 	  tree->mcmc->run_move[tree->mcmc->num_move_tree_height]++; */
  /*         return; */
  /* 	} */
  /*   } */
  
  For(i,2*tree->n_otu-2) tree->rates->br_do_updt[i] = YES;
  RATES_Update_Cur_Bl(tree);
  if(tree->mcmc->use_data) new_lnL_data = Lk(NULL,tree);

  new_lnL_rate = RATES_Lk_Rates(tree);
  new_lnL_time = TIMES_Lk_Times(tree);
 

  /* The Hastings ratio is actually mult^(n) when changing the absolute
     node heights. When considering the relative heights, this ratio combined
     to the Jacobian for the change of variable ends up to being equal to mult. 
  */
  /* ratio += LOG(mult); */
  ratio += (phydbl)(n_nodes)*LOG(mult);

  /* Likelihood ratio */
  if(tree->mcmc->use_data) ratio += (new_lnL_data - cur_lnL_data);

  /* Prior ratio */
  ratio += (new_lnL_rate - cur_lnL_rate);
  ratio += (new_lnL_time - cur_lnL_time);

  /* phydbl diff = (new_lnL_time - cur_lnL_time)+(phydbl)(n_nodes-1.)*LOG(mult); */
  /* printf("\n. diff_lk = %12f -(n-1)*log(mult) = %12f n_nodes=%d [%12f] log(mult)=%f", */
  /* 	 (new_lnL_time - cur_lnL_time), */
  /* 	 -(phydbl)(n_nodes-1.)*LOG(mult), */
  /* 	 n_nodes, */
  /* 	 diff, */
  /* 	 LOG(mult)); */

  /* !!!!!!!!!!!! */
  /* ratio += LOG(Dexp(FABS(new_height-floor),1./10.) / Dexp(FABS(cur_height-floor),1./10.)); */
  
  ratio = EXP(ratio);
  alpha = MIN(1.,ratio);
  u = Uni();

  /* printf("\n. cur_lnL_time: %f new_lnL_time: %f",cur_lnL_time,new_lnL_time); */
  
  if(u > alpha)
    {
      RATES_Reset_Times(tree);
      Restore_Br_Len(tree);
      tree->c_lnL = cur_lnL_data;
      tree->rates->c_lnL_rates = cur_lnL_rate;
      tree->rates->c_lnL_times = cur_lnL_time;
    }
  else
    {
      tree->mcmc->acc_move[tree->mcmc->num_move_tree_height]++;
      tree->mcmc->acc_move[tree->mcmc->num_move_nd_t+tree->n_root->num-tree->n_otu]++;
      /* printf("\n. ACCEPT\n"); */
    }

  tree->mcmc->run_move[tree->mcmc->num_move_tree_height]++;
  tree->mcmc->run_move[tree->mcmc->num_move_nd_t+tree->n_root->num-tree->n_otu]++;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MCMC_Updown_T_Cr(t_tree *tree)
{
  /*! TO DO: make sure to change the values of clock_r across 
    the different data partitions
  */
  
  int i;
  phydbl K,mult,u,alpha,ratio;
  phydbl cur_lnL_data,new_lnL_data;
  phydbl cur_lnL_rate,new_lnL_rate;
  phydbl cur_lnL_time,new_lnL_time;
  phydbl floor;
  int n_nodes;

  /*! Check that sequences are isochronous. */
  For(i,tree->n_otu-1) if(!Are_Equal(tree->rates->nd_t[i+1],tree->rates->nd_t[i],1.E-10)) return; 

  if(FABS(tree->rates->t_prior_max[tree->n_root->num] - tree->rates->t_prior_min[tree->n_root->num]) < 1.E-10) return;

  RATES_Record_Times(tree);
  Record_Br_Len(tree);

  K            = tree->mcmc->tune_move[tree->mcmc->num_move_updown_t_cr];
  cur_lnL_data = tree->c_lnL;
  new_lnL_data = tree->c_lnL;
  ratio        = 0.0;
  cur_lnL_rate = tree->rates->c_lnL_rates;
  new_lnL_rate = tree->rates->c_lnL_rates;
  cur_lnL_time = tree->rates->c_lnL_times;

  u = Uni();
  mult = EXP(K*(u-0.5));

  floor = 0.0;

  Scale_Subtree_Height(tree->n_root,mult,floor,&n_nodes,tree);

  For(i,2*tree->n_otu-1)
    {
      if(tree->rates->nd_t[i] > tree->rates->t_prior_max[i] ||
  	 tree->rates->nd_t[i] < tree->rates->t_prior_min[i])
  	{
  	  RATES_Reset_Times(tree);
	  Restore_Br_Len(tree);
	  tree->mcmc->run_move[tree->mcmc->num_move_updown_t_cr]++;
  	  return;
  	}
    }

  if(RATES_Check_Node_Times(tree)) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);

  tree->rates->clock_r /= mult;
  if(tree->rates->clock_r < tree->rates->min_clock || tree->rates->clock_r > tree->rates->max_clock)
    {
      tree->rates->clock_r *= mult;
      RATES_Reset_Times(tree);
      Restore_Br_Len(tree);
      tree->mcmc->run_move[tree->mcmc->num_move_updown_t_cr]++;
      return;
    }
  
  For(i,2*tree->n_otu-2) tree->rates->br_do_updt[i] = YES;
  RATES_Update_Cur_Bl(tree);
  if(tree->mcmc->use_data) new_lnL_data = Lk(NULL,tree);

  new_lnL_rate = RATES_Lk_Rates(tree);
  new_lnL_time = TIMES_Lk_Times(tree);

  /* The Hastings ratio is actually mult^(n) when changing the absolute
     node heights. When considering the relative heights, this ratio combined
     to the Jacobian for the change of variable ends up to being equal to mult. 
  */
  ratio += (n_nodes - 1)*LOG(mult);
  /* ratio += -LOG(mult) + LOG(Dgamma(1./mult,1./K,K)/Dgamma(mult,1./K,K)); */

  /* Likelihood ratio */
  if(tree->mcmc->use_data) ratio += (new_lnL_data - cur_lnL_data);

  /* Prior ratio */
  ratio += (new_lnL_rate - cur_lnL_rate);
  ratio += (new_lnL_time - cur_lnL_time);

  /* !!!!!!!!!!!!1 */
  /* ratio += LOG(Dexp(FABS(new_height-floor),1./10.) / Dexp(FABS(cur_height-floor),1./10.)); */
  
  ratio = EXP(ratio);
  alpha = MIN(1.,ratio);
  u = Uni();
  
  /* printf("\n. t_old = %f t_new = %f cr_old = %f cr_new = %f", */
  /*        tree->rates->nd_t[tree->n_root->num]/mult, */
  /*        tree->rates->nd_t[tree->n_root->num], */
  /*        tree->rates->clock_r*mult, */
  /*        tree->rates->clock_r); */

  if(u > alpha)
    {
      RATES_Reset_Times(tree);
      tree->rates->clock_r *= mult;
      Restore_Br_Len(tree);
      tree->c_lnL = cur_lnL_data;
      tree->rates->c_lnL_rates = cur_lnL_rate;
      tree->rates->c_lnL_times = cur_lnL_time;
    }
  else
    {
      tree->mcmc->acc_move[tree->mcmc->num_move_updown_t_cr]++;
    }

  tree->mcmc->run_move[tree->mcmc->num_move_updown_t_cr]++;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MCMC_Updown_T_Br(t_tree *tree)
{
  
  int i;
  phydbl K,mult,u,alpha,ratio;
  phydbl cur_lnL_data,new_lnL_data;
  phydbl cur_lnL_rate,new_lnL_rate;
  phydbl cur_lnL_time,new_lnL_time;
  phydbl floor;
  int n_nodes;

  /*! Check that sequences are isochronous. */
  For(i,tree->n_otu-1) if(!Are_Equal(tree->rates->nd_t[i+1],tree->rates->nd_t[i],1.E-10)) return; 

  if(FABS(tree->rates->t_prior_max[tree->n_root->num] - tree->rates->t_prior_min[tree->n_root->num]) < 1.E-10) return;

  RATES_Record_Times(tree);
  Record_Br_Len(tree);

  K            = tree->mcmc->tune_move[tree->mcmc->num_move_updown_t_br];
  cur_lnL_data = tree->c_lnL;
  new_lnL_data = tree->c_lnL;
  ratio        = 0.0;
  cur_lnL_rate = tree->rates->c_lnL_rates;
  new_lnL_rate = tree->rates->c_lnL_rates;
  cur_lnL_time = tree->rates->c_lnL_times;

  u = Uni();
  mult = EXP(K*(u-0.5));


  floor = 0.0;

  Scale_Subtree_Height(tree->n_root,mult,floor,&n_nodes,tree);

  For(i,2*tree->n_otu-1)
    {
      if(tree->rates->nd_t[i] > tree->rates->t_prior_max[i] ||
  	 tree->rates->nd_t[i] < tree->rates->t_prior_min[i])
  	{
  	  RATES_Reset_Times(tree);
	  Restore_Br_Len(tree);
	  tree->mcmc->run_move[tree->mcmc->num_move_updown_t_br]++;
  	  return;
  	}
    }

  if(RATES_Check_Node_Times(tree)) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);

  tree->rates->birth_rate /= mult;

  if(tree->rates->birth_rate < tree->rates->birth_rate_min || tree->rates->birth_rate > tree->rates->birth_rate_max)
    {
      tree->rates->birth_rate *= mult;
      RATES_Reset_Times(tree);
      Restore_Br_Len(tree);
      tree->mcmc->run_move[tree->mcmc->num_move_updown_t_br]++;
      return;
    }
  
  For(i,2*tree->n_otu-2) tree->rates->br_do_updt[i] = YES;
  RATES_Update_Cur_Bl(tree);
  if(tree->mcmc->use_data) new_lnL_data = Lk(NULL,tree);

  new_lnL_rate = RATES_Lk_Rates(tree);
  new_lnL_time = TIMES_Lk_Times(tree);

  /* The Hastings ratio is actually mult^(n) when changing the absolute
     node heights. When considering the relative heights, this ratio combined
     to the Jacobian for the change of variable ends up to being equal to mult. 
  */
  ratio += (n_nodes - 1)*LOG(mult);
  /* ratio += -LOG(mult) + LOG(Dgamma(1./mult,1./K,K)/Dgamma(mult,1./K,K)); */

  /* Likelihood ratio */
  if(tree->mcmc->use_data) ratio += (new_lnL_data - cur_lnL_data);

  /* Prior ratio */
  ratio += (new_lnL_rate - cur_lnL_rate);
  ratio += (new_lnL_time - cur_lnL_time);

  /* !!!!!!!!!!!!1 */
  /* ratio += LOG(Dexp(FABS(new_height-floor),1./10.) / Dexp(FABS(cur_height-floor),1./10.)); */
  
  ratio = EXP(ratio);
  alpha = MIN(1.,ratio);
  u = Uni();
  
  /* printf("\n. t_old = %f t_new = %f br_old = %f br_new = %f mult = %f K=%f", */
  /*        tree->rates->nd_t[tree->n_root->num]/mult, */
  /*        tree->rates->nd_t[tree->n_root->num], */
  /*        tree->rates->birth_rate*mult, */
  /*        tree->rates->birth_rate,mult,K); */

  if(u > alpha)
    {
      RATES_Reset_Times(tree);
      tree->rates->birth_rate *= mult;
      Restore_Br_Len(tree);
      tree->c_lnL = cur_lnL_data;
      tree->rates->c_lnL_rates = cur_lnL_rate;
      tree->rates->c_lnL_times = cur_lnL_time;
    }
  else
    {
      tree->mcmc->acc_move[tree->mcmc->num_move_updown_t_br]++;
    }

  tree->mcmc->run_move[tree->mcmc->num_move_updown_t_br]++;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MCMC_Subtree_Height(t_tree *tree)
{
  int i;
  phydbl K,mult,u,alpha,ratio;
  phydbl cur_lnL_data,new_lnL_data;
  phydbl cur_lnL_rate,new_lnL_rate;
  phydbl cur_lnL_time,new_lnL_time;
  phydbl floor;
  int target;
  int n_nodes;

  RATES_Record_Times(tree);
  Record_Br_Len(tree);

  K = tree->mcmc->tune_move[tree->mcmc->num_move_subtree_height];
  cur_lnL_data = tree->c_lnL;
  new_lnL_data = tree->c_lnL;
  cur_lnL_rate = tree->rates->c_lnL_rates;
  new_lnL_rate = tree->rates->c_lnL_rates;
  ratio        = 0.0;
  cur_lnL_time = tree->rates->c_lnL_times;

  u = Uni();
  mult = EXP(K*(u-0.5));
  /* mult = Rgamma(1./K,K); */

  target = Rand_Int(tree->n_otu,2*tree->n_otu-3);

  floor = tree->rates->t_floor[target];

  if(tree->a_nodes[target] == tree->n_root) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);

  if(!Scale_Subtree_Height(tree->a_nodes[target],mult,floor,&n_nodes,tree))
    {
      RATES_Reset_Times(tree);
      Restore_Br_Len(tree);
      tree->mcmc->run_move[tree->mcmc->num_move_subtree_height]++;
      return;
    }

  
  For(i,2*tree->n_otu-1)
    {
      if(tree->rates->nd_t[i] > tree->rates->t_prior_max[i] ||
  	 tree->rates->nd_t[i] < tree->rates->t_prior_min[i])
  	{
  	  RATES_Reset_Times(tree);
	  Restore_Br_Len(tree);
	  tree->mcmc->run_move[tree->mcmc->num_move_subtree_height]++;
  	  return;
  	}
    }

     
  For(i,2*tree->n_otu-2) tree->rates->br_do_updt[i] = YES;
  RATES_Update_Cur_Bl(tree);
  if(tree->mcmc->use_data) new_lnL_data = Lk(NULL,tree);

  new_lnL_rate = RATES_Lk_Rates(tree);
  new_lnL_time = TIMES_Lk_Times(tree);

  /* The Hastings ratio here is mult^(n_nodes) and the ratio of the prior joint densities
     of the modified node heigths given the unchanged one is 1. This is different from the 
     case where all the nodes, including the root node, are scaled. 
  */
  ratio += (phydbl)(n_nodes)*LOG(mult);

  /* Likelihood ratio */
  if(tree->mcmc->use_data) ratio += (new_lnL_data - cur_lnL_data);

  /* Prior ratio */
  ratio += (new_lnL_rate - cur_lnL_rate);
  ratio += (new_lnL_time - cur_lnL_time);


  ratio = EXP(ratio);
  alpha = MIN(1.,ratio);
  u = Uni();
  
  if(u > alpha)
    {
      RATES_Reset_Times(tree);
      Restore_Br_Len(tree);
      tree->c_lnL = cur_lnL_data;
      tree->rates->c_lnL_rates = cur_lnL_rate;
      tree->rates->c_lnL_times = cur_lnL_time;
    }
  else
    {
      tree->mcmc->acc_move[tree->mcmc->num_move_subtree_height]++;
    }

  tree->mcmc->run_move[tree->mcmc->num_move_subtree_height]++;

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Tree_Rates(t_tree *tree)
{
  phydbl K,mult,u,alpha,ratio;
  phydbl cur_lnL_rate,new_lnL_rate;
  phydbl cur_lnL_data,new_lnL_data;
  int n_nodes;
  phydbl init_clock;
  
  if(tree->rates->model == STRICTCLOCK) return;

  RATES_Record_Rates(tree);
  Record_Br_Len(tree);

  K             = tree->mcmc->tune_move[tree->mcmc->num_move_tree_rates];
  cur_lnL_data  = tree->c_lnL;
  new_lnL_data  = tree->c_lnL;
  cur_lnL_rate  = tree->rates->c_lnL_rates;
  new_lnL_rate  = tree->rates->c_lnL_rates;
  init_clock    = tree->rates->clock_r;
  ratio         = 0.0;

  u = Uni();
  mult = EXP(K*(u-0.5));
  /* mult = Rgamma(1./K,K); */
    
  /* Multiply branch rates (or add to log of rates) */
  if(!Scale_Subtree_Rates(tree->n_root,mult,&n_nodes,tree))
    {
      RATES_Reset_Rates(tree);
      Restore_Br_Len(tree);
      tree->mcmc->run_move[tree->mcmc->num_move_tree_rates]++;
      return;
    }

  if(n_nodes != 2*tree->n_otu-2) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);    

  /* Divide clock_r */
  tree->rates->clock_r /= mult;
  if(tree->rates->clock_r < tree->rates->min_clock || tree->rates->clock_r > tree->rates->max_clock)
    {
      tree->rates->clock_r = init_clock;
      RATES_Reset_Rates(tree);
      Restore_Br_Len(tree);
      tree->mcmc->run_move[tree->mcmc->num_move_tree_rates]++;
      return;
    }

/*   if(tree->rates->model == GUINDON) */
/*     { */
      int i;
      For(i,2*tree->n_otu-2) tree->rates->br_do_updt[i] = YES;
      RATES_Update_Cur_Bl(tree);
      if(tree->mcmc->use_data) new_lnL_data = Lk(NULL,tree);
/*     } */

  new_lnL_rate = RATES_Lk_Rates(tree);


  /* Proposal ratio: 2n-2=> number of multiplications, 1=>number of divisions */
  ratio += (+(2*tree->n_otu-2)-1)*LOG(mult);
  /* ratio += (+(2*tree->n_otu-2)-1-2)*LOG(mult) + LOG(Dgamma(1./mult,1./K,K)/Dgamma(mult,1./K,K)); */

  /* If modelling log of rates instead of rates */
  if(tree->rates->model_log_rates == YES) ratio -= (2*tree->n_otu-2)*LOG(mult);

  /* Prior density ratio */
  ratio += (new_lnL_rate - cur_lnL_rate);

  /* Likelihood density ratio */
  ratio += (new_lnL_data - cur_lnL_data);

  ratio = EXP(ratio);
  alpha = MIN(1.,ratio);
  u = Uni();
  
/*   printf("\n. cur_lnL=%f new_lnL=%f ratio=%f mult=%f %f [%f %f]", */
/*   	 cur_lnL_rate,new_lnL_rate,ratio,mult,(+(2*tree->n_otu-2)-1-2)*LOG(mult) + LOG(Dgamma(1./mult,1./K,K)/Dgamma(mult,1./K,K)),new_lnL_data,cur_lnL_data); */

  if(u > alpha)
    {
/*       PhyML_Printf("\n. Reject mult=%f",mult); */
      tree->rates->clock_r = init_clock;
      RATES_Reset_Rates(tree);
      Restore_Br_Len(tree);
      tree->rates->c_lnL_rates = cur_lnL_rate;
      tree->c_lnL        = cur_lnL_data;
    }
  else
    {
/*       PhyML_Printf("\n. Accept mult=%f",mult); */
      tree->mcmc->acc_move[tree->mcmc->num_move_tree_rates]++;
    }

  tree->mcmc->run_move[tree->mcmc->num_move_tree_rates]++;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Subtree_Rates(t_tree *tree)
{  phydbl K,mult,u,alpha,ratio;
  phydbl cur_lnL_rate,new_lnL_rate;
  phydbl cur_lnL_data,new_lnL_data;
  int target;
  int n_nodes;
  
  if(tree->rates->model == STRICTCLOCK) return;

  RATES_Record_Rates(tree);
  Record_Br_Len(tree);

  K             = tree->mcmc->tune_move[tree->mcmc->num_move_subtree_rates];
  cur_lnL_rate  = tree->rates->c_lnL_rates;
  new_lnL_rate  = cur_lnL_rate;
  cur_lnL_data  = tree->c_lnL;
  new_lnL_data  = cur_lnL_data;
  ratio         = 0.0;

  u = Uni();
  mult = EXP(K*(u-0.5));
  /* mult = Rgamma(1./K,K); */

  target = Rand_Int(tree->n_otu,2*tree->n_otu-3);

  /* Multiply branch rates */
  if(!Scale_Subtree_Rates(tree->a_nodes[target],mult,&n_nodes,tree))
    {
      RATES_Reset_Rates(tree);
      Restore_Br_Len(tree);
      tree->mcmc->run_move[tree->mcmc->num_move_subtree_rates]++;
      return;
    }

  new_lnL_rate = RATES_Lk_Rates(tree);
  if(tree->mcmc->use_data) new_lnL_data = Lk(NULL,tree);

  /* Proposal ratio: 2n-2=> number of multiplications, 1=>number of divisions */
  ratio += (+n_nodes)*LOG(mult);
  /* ratio += (n_nodes-2)*LOG(mult) + LOG(Dgamma(1./mult,1./K,K)/Dgamma(mult,1./K,K)); */
  
  /* If modelling log of rates instead of rates */
  if(tree->rates->model_log_rates == YES) ratio -= (n_nodes)*LOG(mult);

  /* Prior density ratio */
  ratio += (new_lnL_rate - cur_lnL_rate);

  /* Likelihood density ratio */
  ratio += (new_lnL_data - cur_lnL_data);


  ratio = EXP(ratio);
  alpha = MIN(1.,ratio);
  u = Uni();
  
  if(u > alpha)
    {
      RATES_Reset_Rates(tree);
      Restore_Br_Len(tree);
      tree->rates->c_lnL_rates = cur_lnL_rate;
      tree->c_lnL        = cur_lnL_data;
    }
  else
    {
      tree->mcmc->acc_move[tree->mcmc->num_move_subtree_rates]++;
    }

  tree->mcmc->run_move[tree->mcmc->num_move_subtree_rates]++;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Swing(t_tree *tree)
{
  int i;
  phydbl K,mult,u,alpha,ratio;
  phydbl cur_lnL_data,new_lnL_data;
  phydbl cur_lnL_rate,new_lnL_rate;

  if(tree->rates->model == STRICTCLOCK) return;

  RATES_Record_Times(tree);
  RATES_Record_Rates(tree);
  Record_Br_Len(tree);

  K             = 3.;
  cur_lnL_data  = tree->c_lnL;
  new_lnL_data  = cur_lnL_data;
  cur_lnL_rate  = tree->rates->c_lnL_rates;
  new_lnL_rate  = cur_lnL_rate;
  ratio         = 0.0;

  u = Uni();
  /* mult = EXP(K*(u-0.5)); */
  mult = u*(K - 1./K) + 1./K;


  For(i,2*tree->n_otu-1)
    {
      if(tree->a_nodes[i]->tax == NO) 
	{
	  tree->rates->nd_t[i] *= mult;
	}

      if(tree->rates->nd_t[i] > tree->rates->t_prior_max[i] ||
         tree->rates->nd_t[i] < tree->rates->t_prior_min[i])
        {
          RATES_Reset_Times(tree);
          Restore_Br_Len(tree);
          return;
        }
    }

  For(i,2*tree->n_otu-2)
    {
      tree->rates->br_r[i] /= mult;

      if(tree->rates->br_r[i] > tree->rates->max_rate ||
         tree->rates->br_r[i] < tree->rates->min_rate)
        {
          RATES_Reset_Times(tree);
          RATES_Reset_Rates(tree);
          Restore_Br_Len(tree);
          return;
        }
    }

  For(i,2*tree->n_otu-2) tree->rates->br_do_updt[i] = YES;
  RATES_Update_Cur_Bl(tree);
  if(tree->mcmc->use_data) new_lnL_data = Lk(NULL,tree);

  new_lnL_rate = RATES_Lk_Rates(tree);

  ratio += (-(tree->n_otu-1.)-2.)*LOG(mult);
  ratio += (new_lnL_rate - cur_lnL_rate);
  if(tree->mcmc->use_data) ratio += (new_lnL_data - cur_lnL_data);

  ratio = EXP(ratio);
  alpha = MIN(1.,ratio);
  u = Uni();

  if(u > alpha)
    {
      RATES_Reset_Times(tree);
      RATES_Reset_Rates(tree);
      Restore_Br_Len(tree);
      tree->c_lnL = cur_lnL_data;
      tree->rates->c_lnL_rates = cur_lnL_rate;
      /* printf("\n. Reject %8f",mult); */
    }
  else
    {
      /* printf("\n. Accept %8f",mult); */
    }


}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Updown_Nu_Cr(t_tree *tree)
{
  phydbl K,mult,u,alpha,ratio;
  phydbl cur_lnL_rate,new_lnL_rate;
  phydbl cur_lnL_data,new_lnL_data;
  int i;
  
  RATES_Record_Rates(tree);
  Record_Br_Len(tree);

  K             = tree->mcmc->tune_move[tree->mcmc->num_move_updown_nu_cr];
  cur_lnL_data  = tree->c_lnL;
  new_lnL_data  = tree->c_lnL;
  cur_lnL_rate  = tree->rates->c_lnL_rates;
  new_lnL_rate  = tree->rates->c_lnL_rates;
  
  u = Uni();
  mult = EXP(K*(u-0.5));

  /* Multiply branch rates */
  /* if(!Scale_Subtree_Rates(tree->n_root,mult,&n_nodes,tree)) */
  /*   { */
  /*     RATES_Reset_Rates(tree); */
  /*     Restore_Br_Len(tree); */
  /*     tree->mcmc->run_move[tree->mcmc->num_move_updown_nu_cr]++; */
  /*     return; */
  /*   } */

  tree->rates->clock_r /= mult;
  if(tree->rates->clock_r < tree->rates->min_clock || tree->rates->clock_r > tree->rates->max_clock)
    {
      tree->rates->clock_r *= mult;
      RATES_Reset_Rates(tree);
      Restore_Br_Len(tree);
      tree->mcmc->run_move[tree->mcmc->num_move_updown_nu_cr]++;
      return;
    }

  tree->rates->nu *= mult;
  if(tree->rates->nu < tree->rates->min_nu || tree->rates->nu > tree->rates->max_nu)
    {
      tree->rates->nu      /= mult;
      tree->rates->clock_r *= mult;
      RATES_Reset_Rates(tree);
      Restore_Br_Len(tree);
      tree->mcmc->run_move[tree->mcmc->num_move_updown_nu_cr]++;
      return;
    }
  
  For(i,2*tree->n_otu-2) tree->rates->br_do_updt[i] = YES;
  RATES_Update_Cur_Bl(tree);
  if(tree->mcmc->use_data) new_lnL_data = Lk(NULL,tree);

  new_lnL_rate = RATES_Lk_Rates(tree);

  ratio = 0.0;
  /* Proposal ratio: 2n-2=> number of multiplications, 1=>number of divisions */
  /* ratio += n_nodes*LOG(mult); /\* (1-1)*LOG(mult); *\/ */
  ratio += 0.0*LOG(mult); /* (1-1)*LOG(mult); */
  /* Prior density ratio */
  ratio += (new_lnL_rate - cur_lnL_rate);
  /* Likelihood density ratio */
  ratio += (new_lnL_data - cur_lnL_data);

  ratio = EXP(ratio);
  alpha = MIN(1.,ratio);
  u = Uni();
  
  if(u > alpha)
    {
      tree->rates->clock_r *= mult;
      tree->rates->nu /= mult;
      RATES_Reset_Rates(tree);
      Restore_Br_Len(tree);
      tree->rates->c_lnL_rates = cur_lnL_rate;
      tree->c_lnL        = cur_lnL_data;
    }
  else
    {
      tree->mcmc->acc_move[tree->mcmc->num_move_updown_nu_cr]++;
    }

  tree->mcmc->run_move[tree->mcmc->num_move_updown_nu_cr]++;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Print_Param_Stdin(t_mcmc *mcmc, t_tree *tree)
{
  time_t cur_time;
  phydbl min;
  int i;
  time(&cur_time);
  
  min = MDBL_MAX;
  For(i,tree->n_otu-1)
    {
      /*       printf("\n. %d %f %f %f %f",i, */
      /* 	     tree->mcmc->new_param_val[tree->mcmc->num_move_nd_t+i], */
      /* 	     tree->mcmc->old_param_val[tree->mcmc->num_move_nd_t+i], */
      /* 	     tree->mcmc->ess[tree->mcmc->num_move_nd_t+i], */
      /* 	     tree->mcmc->sum_val[tree->mcmc->num_move_nd_t+i]); */

      if(tree->mcmc->ess[tree->mcmc->num_move_nd_t+i] < min)
	min = tree->mcmc->ess[tree->mcmc->num_move_nd_t+i];
    }


  if(mcmc->run == 1)
    {
      PhyML_Printf("\n\n");
      PhyML_Printf("%9s","Run");
      PhyML_Printf("  %5s","Time");
      PhyML_Printf("  %10s","Likelihood");
      PhyML_Printf("  %10s","Prior");
      PhyML_Printf("  %19s","SubstRate[ ESS ]");
      PhyML_Printf("  %17s","TreeHeight[ ESS ]");    
      if(tree->rates->model == THORNE || tree->rates->model == GUINDON) PhyML_Printf("  %16s","AutoCor[ ESS ]");    
      else PhyML_Printf("  %16s","RateVar[ ESS ]");
      PhyML_Printf("  %15s","BirthR[ ESS ]");
      PhyML_Printf("  %8s","MinESS");    
    }

  if((cur_time - mcmc->t_last_print) >  mcmc->print_every)
    {
      mcmc->t_last_print = cur_time;
      PhyML_Printf("\n");
      PhyML_Printf("%9d",tree->mcmc->run);
      PhyML_Printf("  %5d",(int)(cur_time-mcmc->t_beg));
      PhyML_Printf("  %10.2f",tree->c_lnL);
      PhyML_Printf("  %10.2f",(tree->rates ? tree->rates->c_lnL_rates+tree->rates->c_lnL_times : +1));
      PhyML_Printf("  %12.6f[%5.0f]",RATES_Average_Substitution_Rate(tree),tree->mcmc->ess[tree->mcmc->num_move_clock_r]);
      /* PhyML_Printf("\t%12.6f[%5.0f]",tree->rates->clock_r,tree->mcmc->ess[tree->mcmc->num_move_clock_r]); */
      PhyML_Printf("  %10.1f[%5.0f]",
                   (tree->rates ? tree->rates->nd_t[tree->n_root->num] : -1.),
                   tree->mcmc->ess[tree->mcmc->num_move_nd_t+tree->n_root->num-tree->n_otu]);
      PhyML_Printf("  %9f[%5.0f]",
                   (tree->rates ? tree->rates->nu : -1.),
                   tree->mcmc->ess[tree->mcmc->num_move_nu]);
      PhyML_Printf("  %8f[%5.0f]",
                   (tree->rates ? tree->rates->birth_rate : -1.),
                   tree->mcmc->ess[tree->mcmc->num_move_birth_rate]);
      PhyML_Printf("  %8.0f",min);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Print_Param(t_mcmc *mcmc, t_tree *tree)
{
  int i;
  FILE *fp;
  char *s;
  int orig_approx;
  phydbl orig_lnL;
  char *s_tree;
  
  if(tree->mcmc->run > mcmc->chain_len) return;

  s = (char *)mCalloc(100,sizeof(char));

  fp = mcmc->out_fp_stats;

/*   if(tree->mcmc->run == 0) */
/*     { */
/*       PhyML_Fprintf(stdout," ["); */
/*       fflush(NULL); */
/*     } */
  
/*   if(!(mcmc->run%(mcmc->chain_len/10))) */
/*     { */
/*       PhyML_Fprintf(stdout,"."); */
/*       fflush(NULL); */
/*     } */

/*   if(tree->mcmc->run == mcmc->chain_len) */
/*     { */
/*       PhyML_Fprintf(stdout,"]"); */
/*       fflush(NULL); */
/*     } */


/*   MCMC_Print_Means(mcmc,tree); */
/*   MCMC_Print_Last(mcmc,tree); */



  if(!(mcmc->run%mcmc->sample_interval)) 
    {
      MCMC_Copy_To_New_Param_Val(tree->mcmc,tree);

      For(i,tree->mcmc->n_moves)
	{
	  if((tree->mcmc->acc_rate[i] > .1) && (tree->mcmc->start_ess[i] == NO)) tree->mcmc->start_ess[i] = YES;
	  if(tree->mcmc->start_ess[i] == YES) MCMC_Update_Effective_Sample_Size(i,tree->mcmc,tree);
	  if(tree->mcmc->run > (int) tree->mcmc->chain_len * 0.1) tree->mcmc->adjust_tuning[i] = NO;
	}

      if(tree->mcmc->run == 0)
	{
	  time(&(mcmc->t_beg));
	  time(&(mcmc->t_last_print));

	  PhyML_Fprintf(fp,"# Random seed: %d",tree->io->r_seed);
	  PhyML_Fprintf(fp,"\n");
	  PhyML_Fprintf(fp,"Run\t");
/* 	  PhyML_Fprintf(fp,"Time\t"); */
	  /* PhyML_Fprintf(fp,"MeanRate\t"); */
/* 	  PhyML_Fprintf(fp,"NormFact\t"); */

	  /* For(i,mcmc->n_moves) */
	  /*   { */
	  /*     strcpy(s,"Acc."); */
	  /*     PhyML_Fprintf(fp,"%s%d\t",strcat(s,mcmc->move_name[i]),i); */
	  /*   } */

	  /* For(i,mcmc->n_moves) */
	  /*   { */
	  /*     strcpy(s,"Tune."); */
	  /*     PhyML_Fprintf(fp,"%s%d\t",strcat(s,mcmc->move_name[i]),i); */
	  /*   } */

	  /* For(i,mcmc->n_moves) */
	  /*   { */
	  /*     strcpy(s,"Run."); */
	  /*     PhyML_Fprintf(fp,"%s\t",strcat(s,mcmc->move_name[i])); */
	  /*   } */

	  PhyML_Fprintf(fp,"LnLike[Exact]\t");
	  PhyML_Fprintf(fp,"LnLike[Approx]\t");
	  PhyML_Fprintf(fp,"LnRates\t");
	  PhyML_Fprintf(fp,"LnTimes\t");
	  PhyML_Fprintf(fp,"LnPosterior\t");
	  PhyML_Fprintf(fp,"ClockRate\t");
	  PhyML_Fprintf(fp,"EvolRate\t");
	  PhyML_Fprintf(fp,"Nu\t");
	  PhyML_Fprintf(fp,"BirthRate\t");
	  PhyML_Fprintf(fp,"TsTv\t");


	  if(tree->mod->ras->n_catg > 1)
	    {
	      if(tree->mod->ras->free_mixt_rates == NO) PhyML_Fprintf(fp,"Alpha\t");
	      else
		{
		  For(i,tree->mod->ras->n_catg) PhyML_Fprintf(fp,"p%d\t",i);
		  For(i,tree->mod->ras->n_catg) PhyML_Fprintf(fp,"r%d\t",i);
		}
	    }


	  if(tree->mod->m4mod->n_h > 1 && tree->mod->use_m4mod == YES)
	    {
	      For(i,tree->mod->m4mod->n_h) PhyML_Fprintf(fp,"cov_p%d\t",i);
	      For(i,tree->mod->m4mod->n_h) PhyML_Fprintf(fp,"cov_r%d\t",i);
	      PhyML_Fprintf(fp,"cov_switch\t");
	    }


	  if(fp != stdout)
	    {
	      for(i=tree->n_otu;i<2*tree->n_otu-1;i++)
	  	{
	  	  if(tree->a_nodes[i] == tree->n_root->v[2])
		    PhyML_Fprintf(fp,"T%d%s\t",i,"[LeftRoot]");
		  else if(tree->a_nodes[i] == tree->n_root->v[1])
		    PhyML_Fprintf(fp,"T%d%s\t",i,"[RightRoot]");
		  else if(tree->a_nodes[i] == tree->n_root)
		    PhyML_Fprintf(fp,"T%d%s\t",i,"[Root]");
		  else
		    PhyML_Fprintf(fp,"T%d[%d]\t",i,tree->a_nodes[i]->anc->num);
	  	}
	    }


/* 	  if(fp != stdout) */
/* 	    { */
/* 	      For(i,2*tree->n_otu-1) */
/* 	  	{ */
/* 	  	  if(tree->a_nodes[i] == tree->n_root->v[2]) */
/* 	  	    PhyML_Fprintf(fp,"R%d[LeftRoot]\t",i); */
/* 	  	  else if(tree->a_nodes[i] == tree->n_root->v[1]) */
/* 	  	    PhyML_Fprintf(fp,"R%d[RightRoot]\t",i); */
/* 	  	  else if(tree->a_nodes[i] != tree->n_root) */
/* 	  	    PhyML_Fprintf(fp," R%d[%d]\t",i,tree->a_nodes[i]->anc->num); */
/* 		  else */
/* 	  	    PhyML_Fprintf(fp," R%d[Root]\t",i); */

/* /\* 		  PhyML_Fprintf(fp," R%d[%f]\t",i,tree->rates->mean_l[i]); *\/ */
/* 	  	} */
/* 	    } */


	  if(fp != stdout)
	    {
	      For(i,2*tree->n_otu-2)
	  	{
	  	  if(tree->a_nodes[i] == tree->n_root->v[2])
	  	    PhyML_Fprintf(fp,"B%d[LeftRoot]\t",i);
	  	  else if(tree->a_nodes[i] == tree->n_root->v[1])
	  	    PhyML_Fprintf(fp,"B%d[RightRoot]\t",i);
	  	  else
	  	    PhyML_Fprintf(fp," B%d[%d]\t",i,tree->a_nodes[i]->anc->num);

/* 		  PhyML_Fprintf(fp," R%d[%f]\t",i,tree->rates->mean_l[i]); */
	  	}
	    }
	  

      /* 	  if(fp != stdout) */
      /* 	    { */
      /* 	      For(i,2*tree->n_otu-2) */
      /* 	  	{ */
      /* 	  	  if(tree->a_nodes[i] == tree->n_root->v[2]) */
      /* 	  	    PhyML_Fprintf(fp,"G%d[LeftRoot]\t",i); */
      /* 	  	  else if(tree->a_nodes[i] == tree->n_root->v[1]) */
      /* 	  	    PhyML_Fprintf(fp,"G%d[RightRoot]\t",i); */
      /* 	  	  else */
      /* 	  	    PhyML_Fprintf(fp," G%d[%f]\t",i,tree->rates->ml_l[i]); */


      /* /\* 		  PhyML_Fprintf(fp," R%d[%f]\t",i,tree->rates->mean_l[i]); *\/ */
      /* 	  	} */
      /* 	    } */


	  /* if(fp != stdout) */
	  /*   { */
	  /*     For(i,2*tree->n_otu-3) */
	  /*       { */
	  /*         if(tree->a_edges[i] == tree->e_root) */
	  /*           PhyML_Fprintf(fp,"*L[%f]%d\t",i,tree->rates->u_ml_l[i]); */
	  /*         else */
	  /*           PhyML_Fprintf(fp," L[%f]%d\t",i,tree->rates->u_ml_l[i]); */
	  /*       } */
	  /*   } */


	  PhyML_Fprintf(mcmc->out_fp_trees,"#NEXUS\n");
	  PhyML_Fprintf(mcmc->out_fp_trees,"BEGIN TREES;\n");
	  PhyML_Fprintf(mcmc->out_fp_trees,"\tTRANSLATE\n");
	  For(i,tree->n_otu-1) PhyML_Fprintf(mcmc->out_fp_trees,"\t%3d\t%s,\n",tree->a_nodes[i]->num+1,tree->a_nodes[i]->name);
	  PhyML_Fprintf(mcmc->out_fp_trees,"\t%3d\t%s;\n",tree->a_nodes[i]->num+1,tree->a_nodes[i]->name);
	  tree->write_tax_names = NO;
	}

      PhyML_Fprintf(fp,"\n");

      PhyML_Fprintf(fp,"%6d\t",tree->mcmc->run);

/*       time(&mcmc->t_cur); */
/*       PhyML_Fprintf(fp,"%6d\t",(int)(mcmc->t_cur-mcmc->t_beg)); */
      
/*       RATES_Update_Cur_Bl(tree); */
/*       PhyML_Fprintf(fp,"%f\t",RATES_Check_Mean_Rates(tree)); */

/*       PhyML_Fprintf(fp,"%f\t",tree->rates->norm_fact); */
      /* For(i,tree->mcmc->n_moves) PhyML_Fprintf(fp,"%f\t",tree->mcmc->acc_rate[i]); */
      /* For(i,tree->mcmc->n_moves) PhyML_Fprintf(fp,"%f\t",(phydbl)(tree->mcmc->tune_move[i])); */
/*       For(i,tree->mcmc->n_moves) PhyML_Fprintf(fp,"%d\t",(int)(tree->mcmc->run_move[i])); */

      orig_approx = tree->io->lk_approx;
      orig_lnL = tree->c_lnL;
      tree->io->lk_approx = EXACT;
      if(tree->mcmc->use_data)  Lk(NULL,tree);  else tree->c_lnL = 0.0;
      PhyML_Fprintf(fp,"%.1f\t",tree->c_lnL);
      tree->io->lk_approx = NORMAL;
      tree->c_lnL = 0.0;
      if(tree->mcmc->use_data)  Lk(NULL,tree);  else tree->c_lnL = 0.0;
      PhyML_Fprintf(fp,"%.1f\t",tree->c_lnL);
      tree->io->lk_approx = orig_approx;
      tree->c_lnL = orig_lnL;

/*       PhyML_Fprintf(fp,"0\t0\t"); */

      PhyML_Fprintf(fp,"%G\t",tree->rates->c_lnL_rates);
      PhyML_Fprintf(fp,"%G\t",tree->rates->c_lnL_times);
      PhyML_Fprintf(fp,"%G\t",tree->c_lnL+tree->rates->c_lnL_rates+tree->rates->c_lnL_times);
      PhyML_Fprintf(fp,"%G\t",tree->rates->clock_r);
      PhyML_Fprintf(fp,"%G\t",RATES_Average_Substitution_Rate(tree));
      PhyML_Fprintf(fp,"%G\t",tree->rates->nu);
      PhyML_Fprintf(fp,"%G\t",tree->rates->birth_rate);
      PhyML_Fprintf(fp,"%G\t",tree->mod->kappa->v);
      
      if(tree->mod->ras->n_catg > 1)
	{
	  if(tree->mod->ras->free_mixt_rates == NO)
	    PhyML_Fprintf(fp,"%G\t",tree->mod->ras->alpha->v);
	  else
	    {
	      For(i,tree->mod->ras->n_catg) PhyML_Fprintf(fp,"%G\t",tree->mod->ras->gamma_r_proba->v[i]);
	      For(i,tree->mod->ras->n_catg) PhyML_Fprintf(fp,"%G\t",tree->mod->ras->gamma_rr->v[i]);
	      /* For(i,tree->mod->ras->n_catg) PhyML_Fprintf(fp,"%G\t",tree->mod->ras->gamma_r_proba_unscaled[i]); */
	      /* For(i,tree->mod->ras->n_catg) PhyML_Fprintf(fp,"%G\t",tree->mod->ras->gamma_rr_unscaled[i]); */
	    }
	}

      if(tree->mod->m4mod->n_h > 1 && tree->mod->use_m4mod == YES)
	{
	  For(i,tree->mod->m4mod->n_h) PhyML_Fprintf(fp,"%G\t",tree->mod->m4mod->h_fq[i]);
	  For(i,tree->mod->m4mod->n_h) PhyML_Fprintf(fp,"%G\t",tree->mod->m4mod->multipl[i]);
	  PhyML_Fprintf(fp,"%G\t",tree->mod->m4mod->delta);
	}


      char *format = (char *)mCalloc(100,sizeof(char));
      sprintf(format,"%%.%df\t",tree->mcmc->nd_t_digits);
      for(i=tree->n_otu;i<2*tree->n_otu-1;i++) PhyML_Fprintf(fp,format,tree->rates->nd_t[i]);
      Free(format);

      /* for(i=0;i<2*tree->n_otu-1;i++) PhyML_Fprintf(fp,"%.4f\t",LOG(tree->rates->nd_r[i])); */
      
      // Average rate along edges: length divided by elapsed time
      For(i,2*tree->n_otu-2)
	PhyML_Fprintf(fp,"%.4f\t",
		      tree->rates->cur_l[i]/(tree->rates->nd_t[tree->a_nodes[i]->num] - tree->rates->nd_t[tree->a_nodes[i]->anc->num]));

      /* fp_pred = fopen("predict.txt","a");       */
      /* for(i=0;i<2*tree->n_otu-2;i++)  */
      /* 	PhyML_Fprintf(fp_pred,"B%d\t%12f\t%12f\t%4d\n",i,EXP(tree->rates->br_r[i]),tree->rates->nd_t[i],tree->rates->has_survived[i]); */
      /* fclose(fp_pred); */


      /* phydbl p,sd,mean; */
      
      /* for(i=0;i<2*tree->n_otu-2;i++) */
      /* 	{ */
      /* 	  sd = tree->rates->nu * (tree->rates->nd_t[i] - tree->rates->nd_t[tree->a_nodes[i]->anc->num]); */
      /* 	  mean = tree->rates->br_r[tree->a_nodes[i]->anc->num] - .5*sd*sd; */
      /* 	  p = Pnorm(tree->rates->br_r[i],mean,sd); */
      /* 	  PhyML_Fprintf(fp,"%f\t",p); */
      /* 	  tree->rates->mean_r[i] += p; */
      /* 	} */

      /* for(i=0;i<2*tree->n_otu-2;i++) PhyML_Fprintf(fp,"%.4f\t",tree->rates->cur_gamma_prior_mean[i]); */
      /* if(fp != stdout) for(i=tree->n_otu;i<2*tree->n_otu-1;i++) PhyML_Fprintf(fp,"%G\t",tree->rates->t_prior[i]); */
/*       For(i,2*tree->n_otu-3) PhyML_Fprintf(fp,"%f\t",EXP(tree->a_edges[i]->l->v)); */


      /* RATES_Update_Cur_Bl(tree); */
      /* For(i,2*tree->n_otu-3) PhyML_Fprintf(fp,"%f\t",tree->a_edges[i]->l->v); */


      For(i,2*tree->n_otu-2) tree->rates->mean_r[i] = EXP(tree->rates->br_r[i]);
      For(i,2*tree->n_otu-1) tree->rates->mean_t[i] = tree->rates->nd_t[i];
      
      /* Time_To_Branch(tree); */
      /* char *s; */
      /* s = (char *)mCalloc(T_MAX_NAME,sizeof(char)); */
      /* strcpy(s,mcmc->io->in_align_file); */
      /* strcat(s,"_"); */
      /* strcat(s,mcmc->out_filename); */
      /* strcat(s,".ps"); */
      /* DR_Draw_Tree(s,tree); */
      /* Free(s);	 */

      /* FILE *fp; */
      /* int j; */
      /* t_node *d, *v1, *v2; */
      /* int n1, n2; */
      /* phydbl r1, r2; */

      /* s = (char *)mCalloc(T_MAX_NAME,sizeof(char)); */
      /* strcpy(s,mcmc->io->in_align_file); */
      /* strcat(s,"_"); */
      /* strcat(s,mcmc->out_filename); */
      /* strcat(s,".corr"); */
      /* fp = fopen(s,"w"); */

      /* n1 = n2 = 0; */
      /* r1 = r2 = 0.f; */
      /* For(i,2*tree->n_otu-3) */
      /* 	{ */
      /* 	  if(tree->a_nodes[i]->tax == NO) */
      /* 	    { */
      /* 	      d = tree->a_nodes[i]; */
      /* 	      v1 = v2 = NULL; */
      /* 	      For(j,3) */
      /* 		{ */
      /* 		  if(d->v[j] != d->anc && d->b[j] != tree->e_root) */
      /* 		    { */
      /* 		      if(!v1) v1 = d->v[j]; */
      /* 		      else    v2 = d->v[j]; */
      /* 		    } */
      /* 		} */

      /* 	      For(j,3) */
      /* 		{ */
      /* 		  if(v1->v[j] && v1->v[j] == d) */
      /* 		    { */
      /* 		      n1 = v1->bip_size[j]; */
      /* 		      /\* r1 = RATES_Get_Mean_Rate_In_Subtree(v1,tree); *\/ */
      /* 		      r1 = EXP(tree->rates->br_r[v1->num]); */
      /* 		      break; */
      /* 		    } */
      /* 		} */


      /* 	      For(j,3) */
      /* 		{ */
      /* 		  if(v2->v[j] && v2->v[j] == d) */
      /* 		    { */
      /* 		      n2 = v2->bip_size[j]; */
      /* 		      /\* r2 = RATES_Get_Mean_Rate_In_Subtree(v2,tree); *\/ */
      /* 		      r2 = EXP(tree->rates->br_r[v2->num]); */
      /* 		      break; */
      /* 		    } */
      /* 		} */

      /* 	      fprintf(fp,"\n%4d %4d %15f %15f",n1,n2,r1,r2); */
      /* 	    } */
      /* 	} */

      /* fclose(fp); */
      
      // TREES
      Time_To_Branch(tree);
      tree->bl_ndigits = 3;
      s_tree = Write_Tree(tree,NO);
      tree->bl_ndigits = 7;
      PhyML_Fprintf(mcmc->out_fp_trees,"TREE %8d [%f] = [&R] %s\n",mcmc->run,tree->c_lnL,s_tree);
      Free(s_tree);
    }

  if(tree->mcmc->run == mcmc->chain_len) PhyML_Fprintf(mcmc->out_fp_trees,"END;\n");

  fflush(NULL);   

  Free(s);
  
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Print_Means(t_mcmc *mcmc, t_tree *tree)
{

  if(!(mcmc->run%mcmc->sample_interval)) 
    {
      int i;
      char *s;

      s = (char *)mCalloc(T_MAX_FILE,sizeof(char));

      strcpy(s,tree->mcmc->out_filename);
      strcat(s,".means");
      
      fclose(mcmc->out_fp_means);

      mcmc->out_fp_means = fopen(s,"w");
      
      PhyML_Fprintf(mcmc->out_fp_means,"#");
      for(i=tree->n_otu;i<2*tree->n_otu-1;i++) PhyML_Fprintf(mcmc->out_fp_means,"T%d\t",i);	  

      PhyML_Fprintf(mcmc->out_fp_means,"\n");      

      for(i=tree->n_otu;i<2*tree->n_otu-1;i++) tree->rates->t_mean[i] *= (phydbl)(mcmc->run / mcmc->sample_interval);

      for(i=tree->n_otu;i<2*tree->n_otu-1;i++)
	{
	  tree->rates->t_mean[i] += tree->rates->nd_t[i];
	  tree->rates->t_mean[i] /= (phydbl)(mcmc->run / mcmc->sample_interval + 1);

/* 	  PhyML_Fprintf(tree->mcmc->out_fp_means,"%d\t",tree->mcmc->run / tree->mcmc->sample_interval);	   */
	  PhyML_Fprintf(tree->mcmc->out_fp_means,"%.1f\t",tree->rates->t_mean[i]);
	}

      PhyML_Fprintf(tree->mcmc->out_fp_means,"\n");
      fflush(NULL);

      Free(s);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Print_Last(t_mcmc *mcmc, t_tree *tree)
{

  if(!(mcmc->run%mcmc->sample_interval)) 
    {
      int i;
      char *s;

      s = (char *)mCalloc(T_MAX_FILE,sizeof(char));

      strcpy(s,tree->mcmc->out_filename);
      strcat(s,".lasts");
      
      fclose(mcmc->out_fp_last);

      mcmc->out_fp_last = fopen(s,"w");

/*       rewind(mcmc->out_fp_last); */

      PhyML_Fprintf(mcmc->out_fp_last,"#");
      PhyML_Fprintf(tree->mcmc->out_fp_last,"Time\t");

      for(i=tree->n_otu;i<2*tree->n_otu-1;i++)
	PhyML_Fprintf(tree->mcmc->out_fp_last,"T%d\t",i);

      PhyML_Fprintf(tree->mcmc->out_fp_last,"\n");

      if(mcmc->run)
	{
	  time(&(mcmc->t_cur));
	  PhyML_Fprintf(tree->mcmc->out_fp_last,"%d\t",(int)(mcmc->t_cur-mcmc->t_beg));
/* 	  PhyML_Fprintf(tree->mcmc->out_fp_last,"%d\t",(int)(mcmc->t_beg)); */
	}

      for(i=tree->n_otu;i<2*tree->n_otu-1;i++) PhyML_Fprintf(tree->mcmc->out_fp_last,"%.1f\t",tree->rates->nd_t[i]);

      PhyML_Fprintf(tree->mcmc->out_fp_last,"\n");
      fflush(NULL);

      Free(s);
  }
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MCMC_Pause(t_mcmc *mcmc)
{
  char choice;
  char *s;
  int len;

  s = (char *)mCalloc(100,sizeof(char));

  
  if(!(mcmc->run%mcmc->chain_len) && (mcmc->is_burnin == NO))
    {
      PhyML_Printf("\n. Do you wish to stop the analysis [N/y] ");
      if(!scanf("%c",&choice)) Exit("\n");
      if(choice == '\n') choice = 'N';
      else getchar(); /* \n */
      
      Uppercase(&choice);
	  
      switch(choice)
	{
	case 'N': 
	  {	    
	    len = 1E+4;
	    PhyML_Printf("\n. How many extra generations is required [default: 1E+4] ");
	    Getstring_Stdin(s);
	    if(s[0] == '\0') len = 1E+4; 
	    else len = (int)atof(s); 

	    if(len < 0)
	      {
		PhyML_Printf("\n. The value entered must be an integer greater than 0.\n");
		Exit("\n");
	      }	    
	    mcmc->chain_len += len;
	    break;
	  }
	      
	case 'Y': 
	  {
	    PhyML_Printf("\n. Ok. Done.\n");
	    Exit("\n");
	    break;
	  }
	  
	default: 
	  {
	    PhyML_Printf("\n. Please enter 'Y' or 'N'.\n");
	    Exit("\n");
	  }
	}
    }

  Free(s);

}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Terminate()
{
  char choice;
  PhyML_Printf("\n\n. Do you really want to terminate [Y/n]: ");
  if(!scanf("%c",&choice)) Exit("\n");
  if(choice == '\n') choice = 'Y';
  else getchar(); /* \n */      
  Uppercase(&choice);	  
  if(choice == 'Y') raise(SIGTERM);
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Copy_MCMC_Struct(t_mcmc *ori, t_mcmc *cpy, char *filename)
{
  int pid;
  int i;
  
  cpy->use_data           = ori->use_data        ;
  cpy->sample_interval    = ori->sample_interval ;
  cpy->chain_len          = ori->chain_len       ;
  cpy->randomize          = ori->randomize       ;
  cpy->norm_freq          = ori->norm_freq       ;
  cpy->n_moves            = ori->n_moves         ;
  cpy->max_tune           = ori->max_tune        ;
  cpy->min_tune           = ori->min_tune        ;
  cpy->print_every        = ori->print_every     ;
  cpy->is_burnin          = ori->is_burnin       ;
  cpy->is                 = ori->is              ;
  cpy->io                 = ori->io              ;
  cpy->in_fp_par          = ori->in_fp_par       ;
  cpy->nd_t_digits        = ori->nd_t_digits     ;
  cpy->max_lag            = ori->max_lag         ;

  For(i,cpy->n_moves) 
    {
      cpy->start_ess[i]          = ori->start_ess[i];
      cpy->ess_run[i]            = ori->ess_run[i];
      cpy->ess[i]                = ori->ess[i];
      cpy->move_weight[i]        = ori->move_weight[i];
      cpy->run_move[i]           = ori->run_move[i];
      cpy->acc_move[i]           = ori->acc_move[i];
      cpy->prev_run_move[i]      = ori->prev_run_move[i];
      cpy->prev_acc_move[i]      = ori->prev_acc_move[i];
      cpy->acc_rate[i]           = ori->acc_rate[i];
      cpy->tune_move[i]          = ori->tune_move[i];
      strcpy(cpy->move_name[i],ori->move_name[i]);
      cpy->adjust_tuning[i]      = ori->adjust_tuning[i];
    }


  
  if(filename) 
    {
      char *s;

      s = (char *)mCalloc(T_MAX_NAME,sizeof(char));

      strcpy(cpy->out_filename,filename);
      pid = getpid();
      sprintf(cpy->out_filename+strlen(cpy->out_filename),"_%d",pid);

      strcpy(s,cpy->io->in_align_file);
      strcat(s,"_");
      strcat(s,cpy->out_filename);
      strcat(s,"_stats");
      cpy->out_fp_stats = fopen(s,"w");

      strcpy(s,cpy->io->in_align_file);
      strcat(s,"_");
      strcat(s,cpy->out_filename);
      strcat(s,"_trees");
      cpy->out_fp_trees = fopen(s,"w");

      strcpy(s,cpy->io->in_align_file);
      strcat(s,"_");
      strcat(s,cpy->out_filename);
      strcat(s,"_constree");
      cpy->out_fp_constree = fopen(s,"w");
 
      Free(s);
    }
  else 
    {
      cpy->out_fp_stats = stderr;
      cpy->out_fp_trees = stderr;
    }
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Close_MCMC(t_mcmc *mcmc)
{
  fclose(mcmc->out_fp_trees);
  fclose(mcmc->out_fp_stats);
  fclose(mcmc->out_fp_constree);
  /* fclose(mcmc->out_fp_means); */
  /* fclose(mcmc->out_fp_last); */
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Randomize_Kappa(t_tree *tree)
{
  tree->mod->kappa->v = Uni()*5.;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Randomize_Rate_Across_Sites(t_tree *tree)
{
  if(tree->mod->ras->n_catg == 1) return;
  
  if(tree->mod->ras->free_mixt_rates == YES)
    {
      int i;
      For(i,tree->mod->ras->n_catg-1) tree->mod->ras->gamma_r_proba_unscaled->v[i] = Uni()*100.;
      tree->mod->ras->gamma_r_proba_unscaled->v[tree->mod->ras->n_catg-1] = 100.;
      For(i,tree->mod->ras->n_catg) tree->mod->ras->gamma_rr_unscaled->v[i] = (phydbl)i+1.; /* Do not randomize those as their ordering matter */
    }
  else
    {
      tree->mod->ras->alpha->v = Uni()*5.;
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Randomize_Covarion_Rates(t_tree *tree)
{
  if(tree->mod->use_m4mod == NO) return;
 
  if(tree->mod->m4mod->n_h == 1) return;
  
  int i;

  For(i,tree->mod->m4mod->n_h)
    {
      tree->mod->m4mod->multipl_unscaled[i] = (phydbl)i+1.;
      tree->mod->m4mod->h_fq_unscaled[i] = Uni()*(100.-0.01) + 0.01;
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Randomize_Covarion_Switch(t_tree *tree)
{
  if(tree->mod->use_m4mod == NO) return;

  tree->mod->m4mod->delta = Uni()*(10.-0.01)+0.01;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Randomize_Branch_Lengths(t_tree *tree)
{
  int i;

  if(tree->mod->log_l == NO)
    For(i,2*tree->n_otu-3) tree->a_edges[i]->l->v = Rexp(10.);
  else
    For(i,2*tree->n_otu-3) tree->a_edges[i]->l->v = -4* Uni();
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Randomize_Node_Rates(t_tree *tree)
{

  int i,err;
  phydbl mean_r, var_r;
  phydbl min_r, max_r;

  mean_r = 1.0;
  var_r  = 0.5;
  min_r  = tree->rates->min_rate;
  max_r  = tree->rates->max_rate;

  For(i,2*tree->n_otu-2)
    if(tree->a_nodes[i] != tree->n_root)
      tree->rates->nd_r[i] = Rnorm_Trunc(mean_r,SQRT(var_r),min_r,max_r,&err);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Randomize_Rates(t_tree *tree)
{

  /* Should be called once t_node times have been determined */

  int i;

  For(i,2*tree->n_otu-2) tree->rates->br_r[i] = 1.0;

/*   For(i,2*tree->n_otu-2) */
/*     { */
/*       u = Uni(); */
/*       u = u * (r_max-r_min) + r_min; */
/*       tree->rates->br_r[i] = u; */

/*       if(tree->rates->br_r[i] < tree->rates->min_rate) tree->rates->br_r[i] = tree->rates->min_rate;  */
/*       if(tree->rates->br_r[i] > tree->rates->max_rate) tree->rates->br_r[i] = tree->rates->max_rate;  */
/*     } */

  MCMC_Randomize_Rates_Pre(tree->n_root,tree->n_root->v[2],tree);
  MCMC_Randomize_Rates_Pre(tree->n_root,tree->n_root->v[1],tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Randomize_Rates_Pre(t_node *a, t_node *d, t_tree *tree)
{
  int i;
  phydbl mean_r, var_r;
  phydbl min_r, max_r;
  int err;

/*   mean_r = tree->rates->br_r[a->num]; */
/*   var_r  = tree->rates->nu * (tree->rates->nd_t[d->num] - tree->rates->nd_t[a->num]); */
  mean_r = 1.0;
  var_r  = 0.5;
  min_r  = tree->rates->min_rate;
  max_r  = tree->rates->max_rate;
  
  tree->rates->br_r[d->num] = Rnorm_Trunc(mean_r,SQRT(var_r),min_r,max_r,&err);

  if(d->tax) return;
  else
    {
      For(i,3)
	if(d->v[i] != a && d->b[i] != tree->e_root)
	  MCMC_Randomize_Rates_Pre(d,d->v[i],tree);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MCMC_Randomize_Birth(t_tree *tree)
{
  phydbl min_b,max_b;
  phydbl u;

  min_b = tree->rates->birth_rate_min;
  max_b = MIN(0.5,tree->rates->birth_rate_max);
  
  u = Uni();
  tree->rates->birth_rate = (max_b - min_b) * u + min_b;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MCMC_Randomize_Nu(t_tree *tree)
{
  phydbl min_nu,max_nu;
  phydbl u;

  /* It is preferable to start with small values of nu 
     as if is difficult for the MCMC sampler to sample
     equal rates on edge (i.e., molecular clock) since
     such combination of rate lies on the boundary of 
     the space of all edge rate combination. We give
     here a bit of help to the sampler by considering 
     starting points close to the molecular clock
     constraint.
  */
  min_nu = tree->rates->min_nu;
  max_nu = tree->rates->max_nu/10.;
  
  u = Uni();
  tree->rates->nu = (max_nu - min_nu) * u + min_nu;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Randomize_Clock_Rate(t_tree *tree)
{
  phydbl u;
  u = Uni();
  tree->rates->clock_r = u * (tree->rates->max_clock - tree->rates->min_clock) + tree->rates->min_clock;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Randomize_Alpha(t_tree *tree)
{
  phydbl u;

  u = Uni();
  tree->rates->alpha = u*6.0+1.0;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Randomize_Node_Times(t_tree *tree)
{
  phydbl t_sup, t_inf;
  phydbl u;
  int iter;
  int i;
  phydbl dt,min_dt;
  int min_node;

  t_inf = tree->rates->t_prior_min[tree->n_root->num];
  t_sup = tree->rates->t_prior_max[tree->n_root->num];

  u = Uni();
  u *= (t_sup - t_inf);
  u += t_inf;
  
  tree->rates->nd_t[tree->n_root->num] = u;

  MCMC_Randomize_Node_Times_Top_Down(tree->n_root,tree->n_root->v[2],tree);
  MCMC_Randomize_Node_Times_Top_Down(tree->n_root,tree->n_root->v[1],tree);


  min_node = -1;
  iter = 0;
  do
    {
      min_dt = MDBL_MAX;
      For(i,2*tree->n_otu-2) 
	{
	  dt = tree->rates->nd_t[i] - tree->rates->nd_t[tree->a_nodes[i]->anc->num];
	  if(dt < min_dt)
	    {
	      min_dt = dt;
	      min_node = i;
	    }
	}

      if(min_dt > 0.01 * FABS(tree->rates->nd_t[tree->n_root->num])/(phydbl)(tree->n_otu-1)) break;

      RATES_Record_Times(tree);
      For(i,2*tree->n_otu-1) 
	{
	  if(tree->a_nodes[i]->tax == NO) 
	    tree->rates->nd_t[i] -= 0.1*FABS(tree->rates->nd_t[tree->n_root->num])/(phydbl)(tree->n_otu-1);

	  if(tree->rates->nd_t[i] < tree->rates->t_prior_min[i] || 
	     tree->rates->nd_t[i] > tree->rates->t_prior_max[i])
	    {
	      RATES_Reset_Times(tree);
	      break;
	    }
	}

      MCMC_Randomize_Node_Times_Bottom_Up(tree->n_root,tree->n_root->v[2],tree);
      MCMC_Randomize_Node_Times_Bottom_Up(tree->n_root,tree->n_root->v[1],tree);
      

      iter++;
    }
  while(iter < 1000);

  if(iter == 1000)
    {      
      PhyML_Printf("\n== min_dt = %f",min_dt);
      PhyML_Printf("\n== min->t=%f min->anc->t=%f",tree->rates->nd_t[min_node],tree->rates->nd_t[tree->a_nodes[min_node]->anc->num]);
      PhyML_Printf("\n== d up=%f down=%f",tree->rates->t_prior_min[min_node],tree->rates->t_prior_max[min_node]);
      PhyML_Printf("\n== a up=%f down=%f",tree->rates->t_prior_min[tree->a_nodes[min_node]->anc->num],tree->rates->t_prior_max[tree->a_nodes[min_node]->anc->num]);
      PhyML_Printf("\n== up=%f down=%f",tree->rates->t_prior_min[min_node],tree->rates->t_floor[tree->a_nodes[min_node]->anc->num]);
      PhyML_Printf("\n== min_node = %d",min_node);
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
    }


/*   PhyML_Printf("\n. Needed %d iterations to randomize node heights.",iter); */
/*   TIMES_Print_Node_Times(tree->n_root,tree->n_root->v[2],tree); */
/*   TIMES_Print_Node_Times(tree->n_root,tree->n_root->v[1],tree); */

  if(RATES_Check_Node_Times(tree)) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);  
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Randomize_Node_Times_Bottom_Up(t_node *a, t_node *d, t_tree *tree)
{
  if(d->tax) return;
  else
    {
      int i;      
      phydbl u;
      phydbl t_inf, t_sup;
      t_node *v1, *v2;


      For(i,3)
	{
	  if((d->v[i] != a) && (d->b[i] != tree->e_root))
	    {
	      MCMC_Randomize_Node_Times_Bottom_Up(d,d->v[i],tree);
	    }
	}

      v1 = v2 = NULL;
      For(i,3)
      	{
      	  if(d->v[i] != a && d->b[i] != tree->e_root)
      	    {
      	      if(!v1) v1 = d->v[i];
      	      else    v2 = d->v[i];
      	    }
      	}

      t_sup = MIN(tree->rates->nd_t[v1->num],tree->rates->nd_t[v2->num]);
      t_inf = tree->rates->nd_t[a->num];

      u = Uni();
      u *= (t_sup - t_inf);
      u += t_inf;
      
 
      if(u > tree->rates->t_prior_min[d->num] && u < tree->rates->t_prior_max[d->num])
	tree->rates->nd_t[d->num] = u;
    }  
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Randomize_Node_Times_Top_Down(t_node *a, t_node *d, t_tree *tree)
{
  if(d->tax) return;
  else
    {
      int i;      
      phydbl u;
      phydbl t_inf, t_sup;

      t_inf = MAX(tree->rates->nd_t[a->num],tree->rates->t_prior_min[d->num]);
      t_sup = tree->rates->t_prior_max[d->num];

      u = Uni();
      u *= (t_sup - t_inf);
      u += t_inf;
      
      tree->rates->nd_t[d->num] = u;

      For(i,3)
	{
	  if((d->v[i] != a) && (d->b[i] != tree->e_root))
	    {
	      MCMC_Randomize_Node_Times_Top_Down(d,d->v[i],tree);
	    }
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MCMC_Get_Acc_Rates(t_mcmc *mcmc)
{
  int i;
  phydbl eps;
  int lag;


  if(mcmc->run < (int)(0.01*mcmc->chain_len)) lag = 100;
  else lag = 100;

  eps = 1.E-6;

  For(i,mcmc->n_moves)
    {
      if(mcmc->run_move[i] - mcmc->prev_run_move[i] > lag)
	{
	  mcmc->acc_rate[i] = 
	    (phydbl)(mcmc->acc_move[i] - mcmc->prev_acc_move[i] + eps) / 
	    (phydbl)(mcmc->run_move[i] - mcmc->prev_run_move[i] + eps) ;


	  mcmc->prev_run_move[i] = mcmc->run_move[i];
	  mcmc->prev_acc_move[i] = mcmc->acc_move[i];
	  
	  MCMC_Adjust_Tuning_Parameter(i,mcmc);
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Adjust_Tuning_Parameter(int move, t_mcmc *mcmc)
{
  if(mcmc->adjust_tuning[move])
    {
      phydbl scale;
      phydbl rate;
      phydbl rate_inf,rate_sup;
      
      if(mcmc->run < (int)(0.01*mcmc->chain_len)) scale = 1.5;
      else scale = 1.2;

      if(!strcmp(mcmc->move_name[move],"tree_height"))
	{
	  rate_inf = 0.1;
	  rate_sup = 0.1;
	}
      else if(!strcmp(mcmc->move_name[move],"subtree_height"))
	{
	  rate_inf = 0.2;
	  rate_sup = 0.2;
	}
      else if(!strcmp(mcmc->move_name[move],"updown_t_cr"))
	{
	  rate_inf = 0.1;
	  rate_sup = 0.1;
	}
      else if(!strcmp(mcmc->move_name[move],"clock"))
	{
	  rate_inf = 0.1;
	  rate_sup = 0.1;
	}

      /* if(!strcmp(mcmc->move_name[move],"tree_rates")) */
      /* 	{ */
      /* 	  rate_inf = 0.05; */
      /* 	  rate_sup = 0.05; */
      /* 	} */
      else
	{
	  rate_inf = 0.234; // Gareth Robert's magic number !
	  rate_sup = 0.234;
	}


      rate = mcmc->acc_rate[move];
      
      if(rate < rate_inf)      
	{
	  mcmc->tune_move[move] /= scale;
	}

      else if(rate > rate_sup) 
	{
	  mcmc->tune_move[move] *= scale;
	}

      if(mcmc->tune_move[move] > mcmc->max_tune) mcmc->tune_move[move] = mcmc->max_tune;
      if(mcmc->tune_move[move] < mcmc->min_tune) mcmc->tune_move[move] = mcmc->min_tune;
	  
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_One_Length(t_edge *b, t_tree *tree)
{
  phydbl u;
  phydbl new_lnL_data, cur_lnL_data;
  phydbl ratio, alpha;
  phydbl new_l, cur_l;
  phydbl K,mult;


  cur_l        = b->l->v;
  cur_lnL_data = tree->c_lnL;
  K            = 0.1;
  
  u = Uni();
  mult = EXP(K*(u-0.5));
  /* mult = u*(K-1./K)+1./K; */
  new_l = cur_l * mult;

  if(new_l < tree->mod->l_min || new_l > tree->mod->l_max) return;

  b->l->v = new_l;
  new_lnL_data = Lk(b,tree);
/*   tree->both_sides = NO; */
/*   new_lnL_data = Lk(tree); */


  ratio =
    (new_lnL_data - cur_lnL_data) +
    (LOG(mult));


  ratio = EXP(ratio);
  alpha = MIN(1.,ratio);
  
  u = Uni();
  
  if(u > alpha) /* Reject */
    {
      b->l->v = cur_l;
      Update_PMat_At_Given_Edge(b,tree);
      tree->c_lnL = cur_lnL_data;
    }

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Scale_Br_Lens(t_tree *tree)
{
  phydbl u;
  phydbl new_lnL_data, cur_lnL_data;
  phydbl ratio, alpha;
  phydbl K,mult;
  int i;

  Record_Br_Len(tree);

  cur_lnL_data = tree->c_lnL;
  K            = 1.2;
  
  u = Uni();
  mult = u*(K-1./K)+1./K;

  For(i,2*tree->n_otu-3) 
    {
      tree->a_edges[i]->l->v *= mult;
      if(tree->a_edges[i]->l->v < tree->mod->l_min || 
	 tree->a_edges[i]->l->v > tree->mod->l_max) return;
    }

  Set_Both_Sides(NO,tree);
  new_lnL_data = Lk(NULL,tree);

  ratio =
    (new_lnL_data - cur_lnL_data) +
    (2*tree->n_otu-5) * (LOG(mult));

  ratio = EXP(ratio);
  alpha = MIN(1.,ratio);
  
  u = Uni();
  
  if(u > alpha) /* Reject */
    {
      Restore_Br_Len(tree);
      tree->c_lnL = cur_lnL_data;
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Br_Lens(t_tree *tree)
{
  MCMC_Br_Lens_Pre(tree->a_nodes[0],
  		   tree->a_nodes[0]->v[0],
  		   tree->a_nodes[0]->b[0],tree);

  /* int i; */
  /* For(i,2*tree->n_otu-3) */
  /*   { */
  /*     MCMC_One_Length(tree->a_edges[Rand_Int(0,2*tree->n_otu-4)],acc,run,tree); */
  /*   } */
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Br_Lens_Pre(t_node *a, t_node *d, t_edge *b, t_tree *tree)
{
  int i;

  if(a == tree->n_root || d == tree->n_root) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
    
  MCMC_One_Length(b,tree);
  if(d->tax) return;
  else 
    {
      For(i,3) 
	if(d->v[i] != a)
	  {
	    Update_P_Lk(tree,d->b[i],d);
	    MCMC_Br_Lens_Pre(d,d->v[i],d->b[i],tree);
	  }
      Update_P_Lk(tree,b,d);
    }  
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MCMC_Kappa(t_tree *tree)
{
  int change,i;

  change = NO;

  Switch_Eigen(YES,tree->mod);
		
  if(tree->io->lk_approx == NORMAL)
    {
      tree->io->lk_approx = EXACT;
      if(tree->mcmc->use_data == YES) Lk(NULL,tree);
      change = YES;
    }
  
  For(i,2*tree->n_otu-2) tree->rates->br_do_updt[i] = NO;
  MCMC_Single_Param_Generic(&(tree->mod->kappa->v),0.,100.,tree->mcmc->num_move_kappa,
			    NULL,&(tree->c_lnL),
			    NULL,Wrap_Lk,tree->mcmc->move_type[tree->mcmc->num_move_kappa],NO,NULL,tree,NULL);
	  
  if(change == YES)
    {
      tree->io->lk_approx = NORMAL;
      if(tree->mcmc->use_data == YES) Lk(NULL,tree);
    }
  
  Switch_Eigen(NO,tree->mod);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MCMC_Rate_Across_Sites(t_tree *tree)
{
  if(tree->mod->ras->n_catg == 1) return;

  if(tree->mod->ras->free_mixt_rates == YES)
    {
      MCMC_Free_Mixt_Rate(tree);
    }
  else
    {
      MCMC_Alpha(tree);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Alpha(t_tree *tree)
{
  int i;
  
  For(i,2*tree->n_otu-2) tree->rates->br_do_updt[i] = NO;
  MCMC_Single_Param_Generic(&(tree->mod->ras->alpha->v),0.,100.,tree->mcmc->num_move_ras,
			    NULL,&(tree->c_lnL),
			    NULL,Wrap_Lk,tree->mcmc->move_type[tree->mcmc->num_move_ras],NO,NULL,tree,NULL);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Free_Mixt_Rate(t_tree *tree)
{
  phydbl num,denom;
  phydbl Jnow,Jthen;
  phydbl *z,*y;
  phydbl y_cur,z_cur;
  phydbl low_bound,up_bound;
  int c,i;
  int c2updt; 
  phydbl u;
  phydbl hr;
  int n_moves;
  phydbl cur_lnL_data, new_lnL_data;
  phydbl ratio,alpha;
  
  cur_lnL_data = tree->c_lnL;
  new_lnL_data = tree->c_lnL;

  c = tree->mod->ras->n_catg;

  z = tree->mod->ras->gamma_rr_unscaled->v;
  y = tree->mod->ras->gamma_r_proba_unscaled->v;

  num = z[c-1]*(y[c-1]-y[c-2]);
  denom = z[0]*y[0];
  for(i=1;i<c;i++) denom += z[i]*(y[i]-y[i-1]);
  denom = POW(denom,c);
  Jthen = num/denom;

  n_moves = 0;
  do
    {
      n_moves++;

      // Update frequencies

      // Choose the class freq to update at random.
      c2updt = Rand_Int(0,c-2); 

      // Proposal is uniform. Determine upper and lower bounds.
      u = Uni();
      low_bound = (c2updt==0)?(.0):(y[c2updt-1]);
      up_bound = (c2updt==c-1)?(100.0):(y[c2updt+1]);
      y_cur = y[c2updt];
      y[c2updt] = low_bound + u*(up_bound - low_bound);
      
      // Calculate the Jacobian for the change of variable from unscaled 
      // frequencies to the frequencies themselves.
      num = z[c-1]*(y[c-1]-y[c-2]);
      denom = z[0]*y[0];
      for(i=1;i<c;i++) denom += z[i]*(y[i]-y[i-1]);
      denom = POW(denom,c);
      Jnow = num/denom;
 
      hr = Jnow/Jthen;

      new_lnL_data = Lk(NULL,tree);

      /* Metropolis-Hastings step */
      ratio = 0.;
      if(tree->mcmc->use_data == YES) ratio += (new_lnL_data - cur_lnL_data);
      ratio += LOG(hr);
      ratio = EXP(ratio);
      alpha = MIN(1.,ratio);
      
      /* printf("\n. class=%d new_val=%f cur_val=%f ratio=%f hr=%f y=%f denom=%f",c2updt,y[c2updt],y_cur,ratio,Jthen/Jnow,y[c-1],denom); */

      u = Uni();
      if(u > alpha) /* Reject */
 	{
	  y[c2updt] = y_cur;
	  tree->c_lnL = cur_lnL_data;
	}
      else /* Accept */
	{
	  cur_lnL_data = new_lnL_data;
	  /* Update the Jacobian */
	  num = z[c-1]*(y[c-1]-y[c-2]);
	  denom = z[0]*y[0];
	  for(i=1;i<c;i++) denom += z[i]*(y[i]-y[i-1]);
	  denom = POW(denom,c);
	  Jthen = num/denom;
	}
      



      // Update rates

      // Choose the class freq to update at random.
      c2updt = Rand_Int(0,tree->mod->ras->n_catg-1);

      // Proposal move.
      u = Uni();

      /* K = tree->mcmc->tune_move[tree->mcmc->num_move_ras+c2updt+c]; */
      /* z_cur = z[c2updt]; */
      /* mult = EXP(K*(u-0.5)); */
      /* z[c2updt] *= mult; */

      u = Uni();
      low_bound = (c2updt==0)?(.0):(z[c2updt-1]);
      up_bound = (c2updt==c-1)?(100.0):(z[c2updt+1]);
      z_cur = z[c2updt];
      z[c2updt] = low_bound + u*(up_bound - low_bound);

      
      // Calculate the Jacobian for the change of variable from unscaled 
      // frequencies to the frequencies themselves.
      num = z[c-1]*(y[c-1]-y[c-2]);
      denom = z[0]*y[0];
      for(i=1;i<c;i++) denom += z[i]*(y[i]-y[i-1]);
      denom = POW(denom,c);
      Jnow = num/denom;

      hr = Jnow/Jthen;
      
      new_lnL_data = Lk(NULL,tree);

      // Metropolis-Hastings step
      ratio = 0.;
      if(tree->mcmc->use_data == YES) ratio += (new_lnL_data - cur_lnL_data);
      ratio += LOG(hr);
      /* ratio += LOG(mult); */

      ratio = EXP(ratio);
      alpha = MIN(1.,ratio);

      u = Uni();
      if(u > alpha)
	{
	  z[c2updt] = z_cur;
	  tree->c_lnL = cur_lnL_data;
	}
      else
	{
	  cur_lnL_data = new_lnL_data;

	  // Update the Jacobian
	  num = z[c-1]*(y[c-1]-y[c-2]);
	  denom = z[0]*y[0];
	  for(i=1;i<c;i++) denom += z[i]*(y[i]-y[i-1]);
	  denom = POW(denom,c);
	  Jthen = num/denom;
	}

    }while(n_moves != c);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Covarion_Rates(t_tree *tree)
{
  int i, class;
  phydbl u;
  phydbl min,max;

  if(tree->mod->use_m4mod == NO) return;

  Switch_Eigen(YES,tree->mod);

  For(i,2*tree->n_otu-2) tree->rates->br_do_updt[i] = YES;

  class = Rand_Int(0,tree->mod->m4mod->n_h-1);


  min = 0.01;
  max = 100.;
  u = Uni();
  if(u < .5)
    {
      if(!class)
	{
	  min = 0.01;
	  max = tree->mod->m4mod->multipl_unscaled[1];
	}
      else if(class == tree->mod->m4mod->n_h-1)
	{
	  min = tree->mod->m4mod->multipl_unscaled[tree->mod->m4mod->n_h-2];
	  max = +100.;
	}
      else
	{
	  min = MIN(tree->mod->m4mod->multipl_unscaled[class-1],tree->mod->m4mod->multipl_unscaled[class+1]);
	  max = MAX(tree->mod->m4mod->multipl_unscaled[class-1],tree->mod->m4mod->multipl_unscaled[class+1]);
	}

      MCMC_Single_Param_Generic(&(tree->mod->m4mod->multipl_unscaled[class]),min,max,tree->mcmc->num_move_cov_rates+class+tree->mod->m4mod->n_h,
				NULL,&(tree->c_lnL),
				NULL,Wrap_Lk,tree->mcmc->move_type[tree->mcmc->num_move_cov_rates+class+tree->mod->m4mod->n_h],NO,NULL,tree,NULL);
    }
  else
    {
      MCMC_Single_Param_Generic(&(tree->mod->m4mod->h_fq_unscaled[class]),0.01,+100.,tree->mcmc->num_move_cov_rates+class,
				NULL,&(tree->c_lnL),
				NULL,Wrap_Lk,tree->mcmc->move_type[tree->mcmc->num_move_cov_rates+class],NO,NULL,tree,NULL);
    }

  Switch_Eigen(NO,tree->mod);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Covarion_Switch(t_tree *tree)
{
  if(tree->mod->use_m4mod == NO) return;

  Switch_Eigen(YES,tree->mod);
  MCMC_Single_Param_Generic(&(tree->mod->m4mod->delta),0.01,+100.,tree->mcmc->num_move_cov_switch,
			    NULL,&(tree->c_lnL),
			    NULL,Wrap_Lk,tree->mcmc->move_type[tree->mcmc->num_move_cov_switch],NO,NULL,tree,NULL); 
  Switch_Eigen(NO,tree->mod);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MCMC_Birth_Rate(t_tree *tree)
{
 
  /* printf("\n. BIRTH cur_lnL_time: %f ",tree->rates->c_lnL_times); */

#ifdef INVITEE
  tree->rates->update_time_norm_const = YES;
#endif

  MCMC_Single_Param_Generic(&(tree->rates->birth_rate),
			    tree->rates->birth_rate_min,
			    tree->rates->birth_rate_max,
			    tree->mcmc->num_move_birth_rate,
			    &(tree->rates->c_lnL_times),NULL,
			    Wrap_Lk_Times,NULL,tree->mcmc->move_type[tree->mcmc->num_move_birth_rate],NO,NULL,tree,NULL); 

#ifdef INVITEE
  TIMES_Lk_Times(tree); 
  tree->rates->update_time_norm_const = NO;
#endif
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MCMC_Nu(t_tree *tree)
{
  phydbl cur_nu,new_nu,cur_lnL_rate,new_lnL_rate;
  phydbl u,alpha,ratio;
  phydbl min_nu,max_nu;
  phydbl K;
  phydbl cur_lnL_data, new_lnL_data;
  int i;

  Record_Br_Len(tree);

  cur_nu        = -1.0;
  new_nu        = -1.0;
  ratio         =  0.0;

  K = tree->mcmc->tune_move[tree->mcmc->num_move_nu];

  cur_lnL_rate = tree->rates->c_lnL_rates;
  new_lnL_rate = tree->rates->c_lnL_rates;

  cur_lnL_data = tree->c_lnL;
  new_lnL_data = tree->c_lnL;
  
  cur_nu       = tree->rates->nu;

  min_nu = tree->rates->min_nu;
  max_nu = tree->rates->max_nu;

  MCMC_Make_Move(&cur_nu,&new_nu,min_nu,max_nu,&ratio,K,tree->mcmc->move_type[tree->mcmc->num_move_nu]);
  
  if(new_nu < max_nu && new_nu > min_nu) 
    {
      tree->rates->nu = new_nu;
      
      new_lnL_rate = RATES_Lk_Rates(tree);      

      if((tree->rates->model == GUINDON) && (tree->mcmc->use_data))
	{
	  For(i,2*tree->n_otu-2) tree->rates->br_do_updt[i] = YES;
	  new_lnL_data = Lk(NULL,tree);
	}
      
      ratio +=
      	(new_lnL_rate - cur_lnL_rate);

      ratio +=
	(new_lnL_data - cur_lnL_data);


      /* !!!!!!!!!!!!!!!! */
      /* Modelling exp(nu) and making move on nu */
      /* ratio += (new_nu - cur_nu); */
	 
      /* Exponential prior on nu */
      /* ratio += LOG(Dexp(new_nu,10.) / Dexp(cur_nu,10.)); */

 
      ratio = EXP(ratio);
      alpha = MIN(1.,ratio);
      
      u = Uni();
      if(u > alpha) /* Reject */
	{
	  tree->rates->nu    = cur_nu;
	  tree->rates->c_lnL_rates = cur_lnL_rate;
	  tree->c_lnL        = cur_lnL_data;
	  Restore_Br_Len(tree);
	}
      else
	{
	  tree->mcmc->acc_move[tree->mcmc->num_move_nu]++;
	}
    }
  tree->mcmc->run_move[tree->mcmc->num_move_nu]++;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_All_Rates(t_tree *tree)
{
  phydbl cur_lnL_data, new_lnL_data, cur_lnL_rate;
  phydbl u, ratio, alpha;

  new_lnL_data = tree->c_lnL;
  cur_lnL_data = tree->c_lnL;
  cur_lnL_rate = tree->rates->c_lnL_rates;
  ratio        = 0.0;
  
  Record_Br_Len(tree);
  RATES_Record_Rates(tree);
  
  MCMC_Sim_Rate(tree->n_root,tree->n_root->v[2],tree);
  MCMC_Sim_Rate(tree->n_root,tree->n_root->v[1],tree);

  new_lnL_data = Lk(NULL,tree);
  
  ratio += (new_lnL_data - cur_lnL_data);
  ratio = EXP(ratio);

  alpha = MIN(1.,ratio);
  u = Uni();

  if(u > alpha) /* Reject */
    {
      Restore_Br_Len(tree);
      RATES_Reset_Rates(tree);
      tree->rates->c_lnL_rates = cur_lnL_rate;
    }
  else
    {
      tree->rates->c_lnL_rates = RATES_Lk_Rates(tree);;
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


/* Only works when simulating from prior */
void MCMC_Sim_Rate(t_node *a, t_node *d, t_tree *tree)
{
  int err;
  phydbl mean,sd,br_r_a,dt_d;

  br_r_a = tree->rates->br_r[a->num];
  dt_d   = tree->rates->nd_t[d->num] - tree->rates->nd_t[a->num];
  sd     = SQRT(dt_d*tree->rates->nu);
  
  if(tree->rates->model_log_rates == YES)
    {
      mean = br_r_a - .5*sd*sd;
    }
  else
    {
      mean = br_r_a;
    }
      
  if(tree->rates->model == STRICTCLOCK) tree->rates->br_r[d->num] = 1.0;
  else
    tree->rates->br_r[d->num] = Rnorm_Trunc(mean,sd,tree->rates->min_rate,tree->rates->max_rate,&err);

  if(d->tax) return;
  else
    {
      int i;

      For(i,3)
	if(d->v[i] != a && d->b[i] != tree->e_root)
	  MCMC_Sim_Rate(d,d->v[i],tree);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MCMC_Complete_MCMC(t_mcmc *mcmc, t_tree *tree)
{
  int i;
  phydbl sum;

  mcmc->io      = tree->io;
  mcmc->n_moves = 0;

  mcmc->num_move_nd_r = mcmc->n_moves;
  mcmc->n_moves += 2*tree->n_otu-1;

  mcmc->num_move_br_r = mcmc->n_moves;
  mcmc->n_moves += 2*tree->n_otu-2;

  mcmc->num_move_nd_t = mcmc->n_moves;
  mcmc->n_moves += tree->n_otu-1;

  mcmc->num_move_nu                      = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_clock_r                 = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_tree_height             = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_subtree_height          = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_kappa                   = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_tree_rates              = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_subtree_rates           = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_updown_nu_cr            = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_ras                     = mcmc->n_moves; mcmc->n_moves += (tree->mod ? 2*tree->mod->ras->n_catg : 1);
  mcmc->num_move_updown_t_cr             = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_cov_rates               = mcmc->n_moves; mcmc->n_moves += (tree->mod->m4mod ? 2*tree->mod->m4mod->n_h : 1);
  mcmc->num_move_cov_switch              = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_birth_rate              = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_updown_t_br             = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_jump_calibration        = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_geo_lambda              = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_geo_sigma               = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_geo_tau                 = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_geo_updown_tau_lbda     = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_geo_updown_lbda_sigma   = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_geo_dum                 = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_phyrex_lbda             = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_phyrex_mu               = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_phyrex_rad              = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_phyrex_indel_disk       = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_phyrex_move_disk_ct     = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_phyrex_move_disk_ud     = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_phyrex_swap_disk        = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_phyrex_indel_hit        = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_phyrex_move_ldsk        = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_phyrex_spr              = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_phyrex_scale_times      = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_phyrex_ldscape_lim      = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_phyrex_sigsq            = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_phyrex_sim              = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_phyrex_traj             = mcmc->n_moves; mcmc->n_moves += 1;

  mcmc->run_move       = (int *)mCalloc(mcmc->n_moves,sizeof(int));
  mcmc->acc_move       = (int *)mCalloc(mcmc->n_moves,sizeof(int));
  mcmc->prev_run_move  = (int *)mCalloc(mcmc->n_moves,sizeof(int));
  mcmc->prev_acc_move  = (int *)mCalloc(mcmc->n_moves,sizeof(int));
  mcmc->acc_rate       = (phydbl *)mCalloc(mcmc->n_moves,sizeof(phydbl));
  mcmc->move_weight    = (phydbl *)mCalloc(mcmc->n_moves,sizeof(phydbl));
  mcmc->move_type      = (int *)mCalloc(mcmc->n_moves,sizeof(int));
  
  /* TO DO: instead of n_moves here we should have something like n_param */
  mcmc->ess            = (phydbl *)mCalloc(mcmc->n_moves,sizeof(phydbl));
  mcmc->ess_run        = (int *)mCalloc(mcmc->n_moves,sizeof(int));
  mcmc->start_ess      = (int *)mCalloc(mcmc->n_moves,sizeof(int));
  mcmc->adjust_tuning  = (int *)mCalloc(mcmc->n_moves,sizeof(int));
  mcmc->tune_move      = (phydbl *)mCalloc(mcmc->n_moves,sizeof(phydbl));
  mcmc->sampled_val    = (phydbl *)mCalloc((int)mcmc->n_moves*(mcmc->chain_len/mcmc->sample_interval + 1),sizeof(phydbl));
  mcmc->mode           = (phydbl *)mCalloc((int)mcmc->n_moves,sizeof(phydbl));
  mcmc->move_name      = (char **)mCalloc(mcmc->n_moves,sizeof(char *));
  For(i,mcmc->n_moves) mcmc->move_name[i] = (char *)mCalloc(T_MAX_MCMC_MOVE_NAME,sizeof(char));

  For(i,mcmc->n_moves) mcmc->adjust_tuning[i] = YES;

  for(i=mcmc->num_move_br_r;i<mcmc->num_move_br_r+2*tree->n_otu-2;i++) strcpy(mcmc->move_name[i],"br_rate");
  for(i=mcmc->num_move_nd_r;i<mcmc->num_move_nd_r+2*tree->n_otu-1;i++) strcpy(mcmc->move_name[i],"nd_rate");
  for(i=mcmc->num_move_nd_t;i<mcmc->num_move_nd_t+tree->n_otu-1;i++)   strcpy(mcmc->move_name[i],"time");
  strcpy(mcmc->move_name[mcmc->num_move_nu],"nu");
  strcpy(mcmc->move_name[mcmc->num_move_clock_r],"clock");
  strcpy(mcmc->move_name[mcmc->num_move_tree_height],"tree_height");
  strcpy(mcmc->move_name[mcmc->num_move_subtree_height],"subtree_height");
  strcpy(mcmc->move_name[mcmc->num_move_kappa],"kappa");
  strcpy(mcmc->move_name[mcmc->num_move_tree_rates],"tree_rates");
  strcpy(mcmc->move_name[mcmc->num_move_subtree_rates],"subtree_rates");
  strcpy(mcmc->move_name[mcmc->num_move_updown_nu_cr],"updown_nu_cr");
  for(i=mcmc->num_move_ras;i<mcmc->num_move_ras+(tree->mod ? 2*tree->mod->ras->n_catg : 1);i++) 
    strcpy(mcmc->move_name[i],"ras");  
  strcpy(mcmc->move_name[mcmc->num_move_updown_t_cr],"updown_t_cr");
  for(i=mcmc->num_move_cov_rates;i<mcmc->num_move_cov_rates+(tree->mod->m4mod ? 2*tree->mod->m4mod->n_h : 1);i++) 
    strcpy(mcmc->move_name[i],"cov_rates");  
  strcpy(mcmc->move_name[mcmc->num_move_cov_switch],"cov_switch");
  strcpy(mcmc->move_name[mcmc->num_move_birth_rate],"birth_rate");
  strcpy(mcmc->move_name[mcmc->num_move_updown_t_br],"updown_t_br");
  strcpy(mcmc->move_name[mcmc->num_move_jump_calibration],"jump_calibration");
  strcpy(mcmc->move_name[mcmc->num_move_geo_lambda],"geo_lambda");
  strcpy(mcmc->move_name[mcmc->num_move_geo_sigma],"geo_sigma");
  strcpy(mcmc->move_name[mcmc->num_move_geo_tau],"geo_tau");
  strcpy(mcmc->move_name[mcmc->num_move_geo_updown_tau_lbda],"geo_updown_tau_lbda");
  strcpy(mcmc->move_name[mcmc->num_move_geo_updown_lbda_sigma],"geo_updown_lbda_sigma");
  strcpy(mcmc->move_name[mcmc->num_move_geo_dum],"geo_dum");
  strcpy(mcmc->move_name[mcmc->num_move_phyrex_lbda],"phyrex_lbda");
  strcpy(mcmc->move_name[mcmc->num_move_phyrex_mu],"phyrex_mu");
  strcpy(mcmc->move_name[mcmc->num_move_phyrex_rad],"phyrex_rad");
  strcpy(mcmc->move_name[mcmc->num_move_phyrex_indel_disk],"phyrex_indel_disk");
  strcpy(mcmc->move_name[mcmc->num_move_phyrex_move_disk_ct],"phyrex_move_disk_ct");
  strcpy(mcmc->move_name[mcmc->num_move_phyrex_move_disk_ud],"phyrex_move_disk_ud");
  strcpy(mcmc->move_name[mcmc->num_move_phyrex_swap_disk],"phyrex_swap_disk");
  strcpy(mcmc->move_name[mcmc->num_move_phyrex_indel_hit],"phyrex_indel_hit");
  strcpy(mcmc->move_name[mcmc->num_move_phyrex_move_ldsk],"phyrex_move_ldsk");
  strcpy(mcmc->move_name[mcmc->num_move_phyrex_spr],"phyrex_spr");
  strcpy(mcmc->move_name[mcmc->num_move_phyrex_scale_times],"phyrex_scale_times");
  strcpy(mcmc->move_name[mcmc->num_move_phyrex_ldscape_lim],"phyrex_ldscape_lim");
  strcpy(mcmc->move_name[mcmc->num_move_phyrex_sigsq],"phyrex_sigsq");
  strcpy(mcmc->move_name[mcmc->num_move_phyrex_sim],"phyrex_sim");
  strcpy(mcmc->move_name[mcmc->num_move_phyrex_traj],"phyrex_traj");
  
  if(tree->rates && tree->rates->model_log_rates == YES)
    for(i=mcmc->num_move_br_r;i<mcmc->num_move_br_r+2*tree->n_otu-2;i++) mcmc->move_type[i] = MCMC_MOVE_RANDWALK_NORMAL;
  else
    for(i=mcmc->num_move_br_r;i<mcmc->num_move_br_r+2*tree->n_otu-2;i++) mcmc->move_type[i] = MCMC_MOVE_SCALE_THORNE;

  for(i=mcmc->num_move_nd_r;i<mcmc->num_move_nd_r+2*tree->n_otu-1;i++) mcmc->move_type[i] = MCMC_MOVE_SCALE_THORNE;
  /* for(i=mcmc->num_move_nd_t;i<mcmc->num_move_nd_t+tree->n_otu-1;i++)   mcmc->move_type[i] = MCMC_MOVE_RANDWALK_UNIFORM; */
  for(i=mcmc->num_move_nd_t;i<mcmc->num_move_nd_t+tree->n_otu-1;i++)   mcmc->move_type[i] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_nu] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_clock_r] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_tree_height] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_subtree_height] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_kappa] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_tree_rates] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_subtree_rates] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_updown_nu_cr] = MCMC_MOVE_RANDWALK_NORMAL;
  for(i=mcmc->num_move_ras;i<mcmc->num_move_ras+(tree->mod ? 2*tree->mod->ras->n_catg : 1);i++) mcmc->move_type[i] = MCMC_MOVE_RANDWALK_NORMAL;  
  mcmc->move_type[mcmc->num_move_updown_t_cr] = MCMC_MOVE_SCALE_THORNE;
  for(i=mcmc->num_move_cov_rates;i<mcmc->num_move_cov_rates+(tree->mod->m4mod ? 2*tree->mod->m4mod->n_h : 1);i++) mcmc->move_type[i] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_cov_switch] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_birth_rate] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_updown_t_br] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_jump_calibration] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_geo_lambda] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_geo_sigma] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_geo_tau] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_geo_updown_tau_lbda] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_geo_updown_lbda_sigma] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_geo_dum] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_phyrex_lbda] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_phyrex_mu] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_phyrex_rad] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_phyrex_scale_times] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_phyrex_sigsq] = MCMC_MOVE_SCALE_THORNE;

  /* We start with small tuning parameter values in order to have inflated ESS
     for clock_r */
  For(i,mcmc->n_moves)
    {
      switch(mcmc->move_type[i])
	{
	case MCMC_MOVE_RANDWALK_NORMAL:
	  {
	    /* mcmc->tune_move[i] = 1.E-1; */
	    mcmc->tune_move[i] = 100.;
	    break;
	  }
	case MCMC_MOVE_RANDWALK_UNIFORM:
	  {
	    /* mcmc->tune_move[i] = 2.0; */
	    mcmc->tune_move[i] = 20.0;
	    break;
	  }
	case MCMC_MOVE_SCALE_GAMMA:
	  {
	    mcmc->tune_move[i] = 1.0;
	    /* mcmc->tune_move[i] = 10.0; */
	    break;
	  }
	case MCMC_MOVE_SCALE_THORNE:
	  {
	    mcmc->tune_move[i] = 1.0;
	    break;
	  }
        default :
          {
            mcmc->tune_move[i] = 1.0;
	    break;
          }
	}
    }
  
  for(i=mcmc->num_move_br_r;i<mcmc->num_move_br_r+2*tree->n_otu-2;i++) mcmc->move_weight[i] = (phydbl)(1./(2.*tree->n_otu-2)); /* Rates */
  for(i=mcmc->num_move_nd_r;i<mcmc->num_move_nd_r+2*tree->n_otu-1;i++) mcmc->move_weight[i] = 0.0; /* Node rates */
  for(i=mcmc->num_move_nd_t;i<mcmc->num_move_nd_t+tree->n_otu-1;i++)   mcmc->move_weight[i] = (phydbl)(1./(tree->n_otu-1));  /* Times */
  mcmc->move_weight[mcmc->num_move_clock_r]          = 1.0;
  mcmc->move_weight[mcmc->num_move_tree_height]      = 2.0;
  mcmc->move_weight[mcmc->num_move_subtree_height]   = 0.0;
  mcmc->move_weight[mcmc->num_move_nu]               = 2.0;
  mcmc->move_weight[mcmc->num_move_kappa]            = 0.5;
  mcmc->move_weight[mcmc->num_move_tree_rates]       = 1.0;
  mcmc->move_weight[mcmc->num_move_subtree_rates]    = 0.5;
  mcmc->move_weight[mcmc->num_move_updown_nu_cr]     = 0.0;
  for(i=mcmc->num_move_ras;i<mcmc->num_move_ras+(tree->mod ? 2*tree->mod->ras->n_catg : 1);i++) mcmc->move_weight[i] = 0.5*(1./(tree->mod ? (phydbl)tree->mod->ras->n_catg : 1));
  mcmc->move_weight[mcmc->num_move_updown_t_cr]      = 0.0; /* Does not seem to work well (does not give uniform prior on root height
  							      when sampling from prior) */
  for(i=mcmc->num_move_cov_rates;i<mcmc->num_move_cov_rates+(tree->mod->m4mod ? 2*tree->mod->m4mod->n_h : 1);i++) mcmc->move_weight[i] = 0.5*(1./(tree->mod->m4mod ? (phydbl)tree->mod->m4mod->n_h : 1));
  mcmc->move_weight[mcmc->num_move_cov_switch]            = 1.0;
  mcmc->move_weight[mcmc->num_move_birth_rate]            = 2.0;
  mcmc->move_weight[mcmc->num_move_updown_t_br]           = 1.0;
#if defined (INVITEE)
  mcmc->move_weight[mcmc->num_move_jump_calibration]      = 0.1;
#else
  mcmc->move_weight[mcmc->num_move_jump_calibration]      = 0.0;
#endif
  mcmc->move_weight[mcmc->num_move_geo_lambda]            = 1.0;
  mcmc->move_weight[mcmc->num_move_geo_sigma]             = 1.0;
  mcmc->move_weight[mcmc->num_move_geo_tau]               = 1.0;
  mcmc->move_weight[mcmc->num_move_geo_updown_tau_lbda]   = 1.0;
  mcmc->move_weight[mcmc->num_move_geo_updown_lbda_sigma] = 1.0;
  mcmc->move_weight[mcmc->num_move_geo_dum]               = 1.0;

# if defined (PHYREX)
  mcmc->move_weight[mcmc->num_move_phyrex_lbda]                  =  2.0;
  mcmc->move_weight[mcmc->num_move_phyrex_mu]                    =  2.0;
  mcmc->move_weight[mcmc->num_move_phyrex_rad]                   =  2.0;
  mcmc->move_weight[mcmc->num_move_phyrex_sigsq]                 =  0.0;
  mcmc->move_weight[mcmc->num_move_phyrex_indel_disk]            = 10.0;
  mcmc->move_weight[mcmc->num_move_phyrex_move_disk_ct]          =  1.0;
  mcmc->move_weight[mcmc->num_move_phyrex_move_disk_ud]          =  5.0;
  mcmc->move_weight[mcmc->num_move_phyrex_swap_disk]             =  1.0;
  mcmc->move_weight[mcmc->num_move_phyrex_indel_hit]             =  4.0;
  mcmc->move_weight[mcmc->num_move_phyrex_move_ldsk]             =  1.0;
  mcmc->move_weight[mcmc->num_move_phyrex_spr]                   =  1.0;
  mcmc->move_weight[mcmc->num_move_phyrex_scale_times]           =  2.0;
  mcmc->move_weight[mcmc->num_move_phyrex_ldscape_lim]           =  0.0;
  mcmc->move_weight[mcmc->num_move_phyrex_sim]                   =  0.1;
  mcmc->move_weight[mcmc->num_move_phyrex_traj]                  =  1.0;
# else
  mcmc->move_weight[mcmc->num_move_phyrex_lbda]                  = 0.0;
  mcmc->move_weight[mcmc->num_move_phyrex_mu]                    = 0.0;
  mcmc->move_weight[mcmc->num_move_phyrex_rad]                   = 0.0;
  mcmc->move_weight[mcmc->num_move_phyrex_sigsq]                 = 0.0;
  mcmc->move_weight[mcmc->num_move_phyrex_indel_disk]            = 0.0;
  mcmc->move_weight[mcmc->num_move_phyrex_move_disk_ct]          = 0.0;
  mcmc->move_weight[mcmc->num_move_phyrex_move_disk_ud]          = 0.0;
  mcmc->move_weight[mcmc->num_move_phyrex_swap_disk]             = 0.0;
  mcmc->move_weight[mcmc->num_move_phyrex_indel_hit]             = 0.0;
  mcmc->move_weight[mcmc->num_move_phyrex_move_ldsk]             = 0.0;
  mcmc->move_weight[mcmc->num_move_phyrex_spr]                   = 0.0;
  mcmc->move_weight[mcmc->num_move_phyrex_scale_times]           = 0.0;
  mcmc->move_weight[mcmc->num_move_phyrex_ldscape_lim]           = 0.0;
  mcmc->move_weight[mcmc->num_move_phyrex_sim]                   = 0.0;
  mcmc->move_weight[mcmc->num_move_phyrex_traj]                  = 0.0;
#endif

  /* for(i=mcmc->num_move_br_r;i<mcmc->num_move_br_r+2*tree->n_otu-2;i++) mcmc->move_weight[i] = 0.0; /\* Rates *\/ */
  /* for(i=mcmc->num_move_nd_r;i<mcmc->num_move_nd_r+2*tree->n_otu-1;i++) mcmc->move_weight[i] = 0.0; /\* Node rates *\/ */
  /* for(i=mcmc->num_move_nd_t;i<mcmc->num_move_nd_t+tree->n_otu-1;i++)   mcmc->move_weight[i] = 0.0;  /\* Times *\/ */
/*   mcmc->move_weight[mcmc->num_move_clock_r]          = 0.0; */
/*   mcmc->move_weight[mcmc->num_move_tree_height]      = 0.0; */
/*   mcmc->move_weight[mcmc->num_move_subtree_height]   = 0.0; */
/*   mcmc->move_weight[mcmc->num_move_tree_height]      = 0.0; */
/*   mcmc->move_weight[mcmc->num_move_subtree_height]   = 0.0; */
/*   mcmc->move_weight[mcmc->num_move_nu]               = 0.0; */
/*   mcmc->move_weight[mcmc->num_move_kappa]            = 0.0; */
/*   mcmc->move_weight[mcmc->num_move_tree_rates]       = 0.0; */
/*   mcmc->move_weight[mcmc->num_move_subtree_rates]    = 0.0; */
/*   mcmc->move_weight[mcmc->num_move_updown_nu_cr]     = 0.0; */
/*   for(i=mcmc->num_move_ras;i<mcmc->num_move_ras+(tree->mod ? 2*tree->mod->ras->n_catg : 1);i++) mcmc->move_weight[i] = 0.0; */
/*   mcmc->move_weight[mcmc->num_move_updown_t_cr]      = 0.0; /\* Does not seem to work well (does not give uniform prior on root height */
/*   							      when sampling from prior) *\/ */
/*   for(i=mcmc->num_move_cov_rates;i<mcmc->num_move_cov_rates+(tree->mod ? 2*tree->mod->m4mod->n_h : 1);i++) mcmc->move_weight[i] = 0.5*(1./(tree->mod ? (phydbl)tree->mod->m4mod->n_h : 1)); */
/*   mcmc->move_weight[mcmc->num_move_cov_switch]       = 0.0; */
/*   mcmc->move_weight[mcmc->num_move_birth_rate]       = 0.0; */
/*   mcmc->move_weight[mcmc->num_move_updown_t_br]      = 0.0; */
/* #if defined (INVITEE) */
/*   mcmc->move_weight[mcmc->num_move_jump_calibration] = 1.0; */
/* #else */
/*   mcmc->move_weight[mcmc->num_move_jump_calibration] = 0.0; */
/* #endif */
/*   mcmc->move_weight[mcmc->num_move_geo_lambda]       = 0.0; */
/*   mcmc->move_weight[mcmc->num_move_geo_sigma]        = 0.0; */
/*   mcmc->move_weight[mcmc->num_move_geo_tau]          = 0.0; */
/*   mcmc->move_weight[mcmc->num_move_geo_updown_tau_lbda]   = 0.0; */
/*   mcmc->move_weight[mcmc->num_move_geo_updown_lbda_sigma] = 0.0; */


  sum = 0.0;
  For(i,mcmc->n_moves) sum += mcmc->move_weight[i];
  For(i,mcmc->n_moves) mcmc->move_weight[i] /= sum;
  for(i=1;i<mcmc->n_moves;i++) mcmc->move_weight[i] += mcmc->move_weight[i-1];
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Initialize_Param_Val(t_mcmc *mcmc, t_tree *tree)
{
  /* int i; */

  /* mcmc->lag_param_val[mcmc->num_move_nu]          = tree->rates->nu; */
  /* mcmc->lag_param_val[mcmc->num_move_clock_r]     = tree->rates->clock_r; */
  /* mcmc->lag_param_val[mcmc->num_move_tree_height] = tree->rates->nd_t[tree->n_root->num]; */
  /* mcmc->lag_param_val[mcmc->num_move_kappa]       = tree->mod->kappa->v; */

  /* For(i,2*tree->n_otu-2) */
  /*   mcmc->lag_param_val[mcmc->num_move_br_r+i] = tree->rates->br_r[i]; */
  
  /* For(i,tree->n_otu-1) */
  /*   mcmc->lag_param_val[mcmc->num_move_nd_t+i] = tree->rates->nd_t[tree->n_otu+i]; */

  /* For(i,2*tree->n_otu-1) */
  /*   mcmc->lag_param_val[mcmc->num_move_nd_r+i] = tree->rates->nd_r[i]; */

  /* For(i,mcmc->n_moves) tree->mcmc->new_param_val[i] = tree->mcmc->lag_param_val[i]; */
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Copy_To_New_Param_Val(t_mcmc *mcmc, t_tree *tree)
{
  mcmc->sampled_val[mcmc->num_move_nu*mcmc->sample_size+mcmc->sample_num]          = tree->rates->nu;
  mcmc->sampled_val[mcmc->num_move_clock_r*mcmc->sample_size+mcmc->sample_num]     = tree->rates->clock_r;
  mcmc->sampled_val[mcmc->num_move_tree_height*mcmc->sample_size+mcmc->sample_num] = tree->rates->nd_t[tree->n_root->num];
  mcmc->sampled_val[mcmc->num_move_kappa*mcmc->sample_size+mcmc->sample_num]       = tree->mod ? tree->mod->kappa->v : -1.;
  mcmc->sampled_val[mcmc->num_move_birth_rate*mcmc->sample_size+mcmc->sample_num]  = tree->rates->birth_rate;

  /* For(i,2*tree->n_otu-2) */
  /*   mcmc->sampled_val[(mcmc->num_move_br_r+i)*mcmc->sample_size+mcmc->sample_num] = tree->rates->br_r[i]; */
  
  /* For(i,tree->n_otu-1) */
  /*   mcmc->sampled_val[(mcmc->num_move_nd_t+i)*mcmc->sample_size+mcmc->sample_num] = tree->rates->nd_t[tree->n_otu+i]; */

  /* For(i,2*tree->n_otu-1) */
  /*   mcmc->sampled_val[(mcmc->num_move_nd_r+i)*mcmc->sample_size+mcmc->sample_num] = tree->rates->nd_r[i]; */

  mcmc->sampled_val[mcmc->num_move_geo_tau*mcmc->sample_size+mcmc->sample_num]      = tree->geo ? tree->geo->tau   : -1.;
  mcmc->sampled_val[mcmc->num_move_geo_lambda*mcmc->sample_size+mcmc->sample_num]   = tree->geo ? tree->geo->lbda  : -1.;
  mcmc->sampled_val[mcmc->num_move_geo_sigma*mcmc->sample_size+mcmc->sample_num]    = tree->geo ? tree->geo->sigma : -1.;
  mcmc->sampled_val[mcmc->num_move_geo_dum*mcmc->sample_size+mcmc->sample_num]      = tree->geo ? tree->geo->dum   : -1.;  
  #ifdef PHYREX
  mcmc->sampled_val[mcmc->num_move_phyrex_lbda*mcmc->sample_size+mcmc->sample_num]  = tree->mmod ? tree->mmod->lbda          : -1.;
  mcmc->sampled_val[mcmc->num_move_phyrex_mu*mcmc->sample_size+mcmc->sample_num]    = tree->mmod ? tree->mmod->mu            : -1.;  
  mcmc->sampled_val[mcmc->num_move_phyrex_sigsq*mcmc->sample_size+mcmc->sample_num] = tree->mmod ? PHYREX_Update_Sigsq(tree) : -1.;
  mcmc->sampled_val[mcmc->num_move_phyrex_rad*mcmc->sample_size+mcmc->sample_num]   = tree->mmod ? tree->mmod->rad           : -1.;
  #endif
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MCMC_Slice_One_Rate(t_node *a, t_node *d, int traversal, t_tree *tree)
{
  phydbl L,R; /* Left and Right limits of the slice */
  phydbl w; /* window width */
  phydbl u;
  phydbl x0,x1;
  phydbl logy;
  t_edge *b;
  int i;



  b = NULL;
  if(a == tree->n_root) b = tree->e_root;
  else For(i,3) if(d->v[i] == a) { b = d->b[i]; break; }
  
  w = 0.05;
  /* w = 10.; */

  x0 = tree->rates->br_r[d->num];
  logy = tree->c_lnL+tree->rates->c_lnL_rates - Rexp(1.);
  
  u = Uni();

  L = x0 - w*u;  
  R = L + w;
  
  do
    {
      tree->rates->br_r[d->num] = L;
      tree->rates->br_do_updt[d->num] = YES;
      RATES_Update_Cur_Bl(tree);
      Lk(b,tree);
      RATES_Lk_Rates(tree);
      if(L < tree->rates->min_rate) { L = tree->rates->min_rate - w; break;}
      L = L - w;
    }
  while(tree->c_lnL + tree->rates->c_lnL_rates > logy);
  L = L + w;


  do
    {
      tree->rates->br_r[d->num] = R;
      tree->rates->br_do_updt[d->num] = YES;
      RATES_Update_Cur_Bl(tree);
      Lk(b,tree);
      RATES_Lk_Rates(tree);
      if(R > tree->rates->max_rate) { R = tree->rates->max_rate + w; break;}
      R = R + w;
    }
  while(tree->c_lnL + tree->rates->c_lnL_rates > logy);
  R = R - w;


  for(;;)
    {
      u = Uni();
      x1 = L + u*(R-L);

      tree->rates->br_r[d->num] = x1;
      tree->rates->br_do_updt[d->num] = YES;
      RATES_Update_Cur_Bl(tree);
      Lk(b,tree);
      RATES_Lk_Rates(tree);
      
      if(tree->c_lnL + tree->rates->c_lnL_rates > logy) break;
      
      if(x1 < x0) L = x1;
      else        R = x1;
    }


  if(traversal == YES)
    {
      if(d->tax == YES) return;
      else
	{
	  For(i,3)
	    if(d->v[i] != a && d->b[i] != tree->e_root)
	      {
		if(tree->io->lk_approx == EXACT && tree->mcmc->use_data) Update_P_Lk(tree,d->b[i],d);
		/* if(tree->io->lk_approx == EXACT && tree->mcmc->use_data) { tree->both_sides = YES; Lk(tree); } */
		MCMC_Slice_One_Rate(d,d->v[i],YES,tree);
	      }
	}
      if(tree->io->lk_approx == EXACT && tree->mcmc->use_data) Update_P_Lk(tree,b,d);
      /* if(tree->io->lk_approx == EXACT && tree->mcmc->use_data) { tree->both_sides = YES; Lk(tree); } */
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MCMC_Make_Move(phydbl *cur, phydbl *new, phydbl inf, phydbl sup, phydbl *loghr, phydbl tune, int move_type)
{
  phydbl u;


      
  switch(move_type)
    {

    case MCMC_MOVE_RANDWALK_NORMAL:
      {
	*new   = *cur + Rnorm(0.0,tune);
	/* Do not use reflection here */
	*loghr = 0.0;

	break;
      }

    case MCMC_MOVE_RANDWALK_UNIFORM:
      {
	u = Uni();  
	/* *new   = u * (2.*tune) + (*cur) - tune; */
	/* *new   = Reflect(*new,inf,sup); */
	*new   = u*(sup-inf)+inf;  
	*loghr = 0.0;
	break;
      }

    case MCMC_MOVE_SCALE_THORNE:
      {
	u = Uni();
	*new   = (*cur) * EXP(tune*(u-.5));
	*loghr = LOG((*new)/(*cur));
	break;
      }

    case MCMC_MOVE_SCALE_GAMMA:
      {
	phydbl r;
	*new = (*cur) * Rgamma(1./tune,tune);
	r = (*new)/(*cur);
	*loghr = -LOG(r) + LOG(Dgamma(1./r,1./tune,tune)/Dgamma(r,1./tune,tune));
	break;
      }

    default : 
      {
	PhyML_Printf("\n. Move not implemented");
	Exit("");
	break;
      }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MCMC_Read_Param_Vals(t_tree *tree)
{
  char *token;
  int sizemax;
  FILE *in_fp;
  phydbl val;
  int i,v;
  
  in_fp = tree->mcmc->in_fp_par;
  
  sizemax = T_MAX_LINE;

  token = (char *)mCalloc(sizemax,sizeof(char));
  
  if(fgets(token,sizemax,in_fp) == NULL) // Skip first line
    {
      PhyML_Printf("\n== Wrong file format.");
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);    
    }
    
  if(fgets(token,sizemax,in_fp) == NULL) // Skip second
    {
      PhyML_Printf("\n== Wrong file format.");
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
    }

  v=fscanf(in_fp,"%lf\t",&val); // Run
  /* PhyML_Printf("\n. Run = %d",(int)val); */
  v=fscanf(in_fp,"%lf\t",&val); // LnLike[Exact]
  /* PhyML_Printf("\n. LnLike = %f",val); */
  v=fscanf(in_fp,"%lf\t",&val); // LnLike[Approx]
  /* PhyML_Printf("\n. LnLike = %f",val); */
  v=fscanf(in_fp,"%lf\t",&val); // LnPriorRate
  /* PhyML_Printf("\n. LnPrior = %f",val); */
  v=fscanf(in_fp,"%lf\t",&val); // LnPriorTime
  /* PhyML_Printf("\n. LnPrior = %f",val); */
  v=fscanf(in_fp,"%lf\t",&val); // LnPosterior
  /* PhyML_Printf("\n. LnPost = %f",val); */

  v=fscanf(in_fp,"%lf\t",&val); // ClockRate
  tree->rates->clock_r = val;
  /* PhyML_Printf("\n. Clock = %f",val); */

  v=fscanf(in_fp,"%lf\t",&val); // EvolRate

  v=fscanf(in_fp,"%lf\t",&val); // Nu
  tree->rates->nu = val;
  /* PhyML_Printf("\n. Nu = %f",val); */

  v=fscanf(in_fp,"%lf\t",&val); // Birth rate
  tree->rates->birth_rate = val;

  v=fscanf(in_fp,"%lf\t",&val); // TsTv
  tree->mod->kappa->v = val;
  /* PhyML_Printf("\n. TsTv = %f",val); */

  For(i,tree->n_otu-1)
    {
      v=fscanf(in_fp,"%lf\t",&val); // Node heights
      tree->rates->nd_t[i+tree->n_otu] = val;
    }

  For(i,2*tree->n_otu-2)
    {
      v=fscanf(in_fp,"%lf\t",&val); // Edge average rates
      tree->rates->br_r[i] = LOG(val);
    }
  
  v++;

  Free(token);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

#ifdef PHYREX
void MCMC_PHYREX_Lbda(t_tree *tree)
{
  MCMC_Single_Param_Generic(&(tree->mmod->lbda),
                            tree->mmod->min_lbda,
                            tree->mmod->max_lbda,
                            tree->mcmc->num_move_phyrex_lbda,
                            NULL,&(tree->mmod->c_lnL),
                            NULL,PHYREX_Wrap_Lk,
                            tree->mcmc->move_type[tree->mcmc->num_move_phyrex_lbda],
                            NO,NULL,tree,NULL);
  
}
#endif

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

#ifdef PHYREX
void MCMC_PHYREX_Mu(t_tree *tree)
{
  MCMC_Single_Param_Generic(&(tree->mmod->mu),
                            tree->mmod->min_mu,
                            tree->mmod->max_mu,
                            tree->mcmc->num_move_phyrex_mu,
                            NULL,&(tree->mmod->c_lnL),
                            NULL,PHYREX_Wrap_Lk,
                            tree->mcmc->move_type[tree->mcmc->num_move_phyrex_mu],
                            NO,NULL,tree,NULL);
  
}
#endif

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

#ifdef PHYREX
void MCMC_PHYREX_Radius(t_tree *tree)
{
  MCMC_Single_Param_Generic(&(tree->mmod->rad),
                            tree->mmod->min_rad,
                            tree->mmod->max_rad,
                            tree->mcmc->num_move_phyrex_rad,
                            NULL,&(tree->mmod->c_lnL),
                            NULL,PHYREX_Wrap_Lk,
                            tree->mcmc->move_type[tree->mcmc->num_move_phyrex_rad],
                            NO,NULL,tree,NULL);

}
#endif 

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

#ifdef PHYREX
void MCMC_PHYREX_Sigsq(t_tree *tree)
{
  tree->mmod->update_rad = YES;
  MCMC_Single_Param_Generic(&(tree->mmod->sigsq),
                            tree->mmod->min_sigsq,
                            tree->mmod->max_sigsq,
                            tree->mcmc->num_move_phyrex_sigsq,
                            NULL,&(tree->mmod->c_lnL),
                            NULL,PHYREX_Wrap_Lk,
                            tree->mcmc->move_type[tree->mcmc->num_move_phyrex_sigsq],
                            NO,NULL,tree,NULL);
  tree->mmod->rad = PHYREX_Update_Radius(tree);
  tree->mmod->update_rad = NO;
}
#endif 

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
#ifdef PHYREX
void MCMC_PHYREX_Indel_Disk(t_tree *tree)
{
  int n_disks_cur, n_disks_new;
  t_dsk  *disk;
  phydbl hr,K;

  K  = 1.0;
  hr = 0.0;

  disk = tree->disk->prev;
  n_disks_cur = 0;
  do
    {
      if(!disk->ldsk) n_disks_cur++;
      disk = disk->prev;
    }
  while(disk && disk->prev);

  /* n_disks_new = (int)Rpois(K*n_disks_cur); */
  
  /* hr += Dpois(n_disks_cur,K*n_disks_new,YES); */
  /* hr -= Dpois(n_disks_new,K*n_disks_cur,YES); */

  phydbl e,v;

  e = (phydbl)n_disks_cur;
  v = 1.+e*5.;

  n_disks_new = (int)Rnbinom(e*e/(v-e),e/v);
  
  hr += Dnbinom(n_disks_cur,e*e/(v-e),e/v,YES);
  hr -= Dnbinom(n_disks_new,e*e/(v-e),e/v,YES);

  if(n_disks_new < n_disks_cur) MCMC_PHYREX_Delete_Disk(hr, n_disks_cur - n_disks_new ,tree);
  else                          MCMC_PHYREX_Insert_Disk(hr, n_disks_new - n_disks_cur,tree);  
}
#endif

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
/* Insert or delete a disk that does not affect any ldisk */
#ifdef PHYREX
void MCMC_PHYREX_Delete_Disk(phydbl hr, int n_delete_disks, t_tree *tree)
{
  phydbl u,alpha,ratio;
  phydbl cur_glnL, new_glnL;
  phydbl T;
  t_dsk  *disk,**target_disk,**valid_disks;
  int i,j,block,n_valid_disks,*permut;
  phydbl new_lbda, cur_lbda;
  phydbl new_rad, cur_rad;
  phydbl mean,sd;
  int n_tot_disks_cur,n_tot_disks_new;
  int err;

  if(n_delete_disks == 0) return;

  valid_disks     = NULL;
  disk            = NULL;
  new_glnL        = tree->mmod->c_lnL;
  cur_glnL        = tree->mmod->c_lnL;
  ratio           = 0.0;
  block           = 100;
  cur_lbda        = tree->mmod->lbda;
  cur_rad         = tree->mmod->rad;
  n_tot_disks_cur = PHYREX_Total_Number_Of_Intervals(tree);


  if(tree->disk->next) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  disk = tree->disk->prev;

  if(!disk->prev) return;

  n_valid_disks = 0;
  do
    {
      if(!disk->ldsk)
        {
          if(!n_valid_disks) valid_disks = (t_dsk **)mCalloc(block,sizeof(t_dsk *));
          else if(!(n_valid_disks%block)) valid_disks = (t_dsk **)mRealloc(valid_disks,n_valid_disks+block,sizeof(t_dsk *));
          valid_disks[n_valid_disks] = disk;
          n_valid_disks++;
        }
      disk = disk->prev;
    }
  while(disk && disk->prev);

  if(!n_valid_disks) return;

  if(n_valid_disks - n_delete_disks < 0) return; 


  target_disk = (t_dsk **)mCalloc(n_delete_disks,sizeof(t_dsk *));

  permut = Permutate(n_valid_disks);

  For(j,n_delete_disks)
    {
      /* Uniform selection of a disk where no coalescent nor 'hit' occurred */
      target_disk[j] = valid_disks[permut[j]];
  
      PHYREX_Remove_Disk(target_disk[j]);
  
      For(i,tree->mmod->n_dim) hr += LOG(1./tree->mmod->lim->lonlat[i]);
    }

  T = PHYREX_Tree_Height(tree);

  hr += LnChoose(n_valid_disks,n_delete_disks);
  hr -= n_delete_disks * LOG(-T);
  hr += LnFact(n_delete_disks);

  Free(valid_disks);

  /* Adjust value of lambda */
  n_tot_disks_new = n_tot_disks_cur - n_delete_disks;

  mean = n_tot_disks_new/FABS(T);
  sd   = 0.1*n_tot_disks_new/FABS(T);
  new_lbda = Rnorm_Trunc(mean,sd,
                         tree->mmod->min_lbda,
                         mean + 2.*sd,
                         &err);
  
  mean = n_tot_disks_cur/FABS(T);
  sd   = 0.1*n_tot_disks_cur/FABS(T);
  hr += Log_Dnorm_Trunc(cur_lbda,
                        mean,sd,
                        tree->mmod->min_lbda,
                        mean + 2.*sd,
                        &err);

  mean = n_tot_disks_new/FABS(T);
  sd   = 0.1*n_tot_disks_new/FABS(T);
  hr -= Log_Dnorm_Trunc(new_lbda,
                        mean,sd,
                        tree->mmod->min_lbda,
                        mean + 2.*sd,
                        &err);
  
  tree->mmod->lbda = new_lbda;



  /* Adjust value of radius */
  mean = 1.1*cur_rad;
  sd   = 0.1*mean;
  new_rad = Rnorm_Trunc(mean,sd,
                        tree->mmod->min_rad,
                        tree->mmod->max_rad,
                        &err);

  
  mean = 1.1*new_rad;
  sd   = 0.1*mean;
  hr += Log_Dnorm_Trunc(cur_rad,
                        mean,sd,
                        tree->mmod->min_rad,
                        tree->mmod->max_rad,
                        &err);

  mean = 1.1*cur_rad;
  sd   = 0.1*mean;
  hr -= Log_Dnorm_Trunc(new_rad,
                        mean,sd,
                        tree->mmod->min_rad,
                        tree->mmod->max_rad,
                        &err);

  tree->mmod->rad = new_rad;


  new_glnL = PHYREX_Lk(tree);
  ratio += (new_glnL - cur_glnL);
  ratio += hr;
  
  if(new_glnL < UNLIKELY + 0.1) 
    {
      PhyML_Printf("\n== %f %f ",cur_glnL,new_glnL);
      PhyML_Printf("\n== %s ",Write_Tree(tree,NO));
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
    }

  ratio = EXP(ratio);
  alpha = MIN(1.,ratio);

  /* Always accept move */
  if(tree->mcmc->always_yes == YES && new_glnL > UNLIKELY) alpha = 1.0;

  u = Uni();
  
  /* printf("\n- Delete new_glnL: %f [%f] hr: %f u:%f alpha: %f",new_glnL,cur_glnL,hr,u,alpha); */

  if(u > alpha) /* Reject */
    {
      tree->mmod->lbda = cur_lbda;
      tree->mmod->rad  = cur_rad;

      /* printf("\nD Reject %f",target_disk[0]->time); */
      For(j,n_delete_disks) PHYREX_Insert_Disk(target_disk[j],tree);

      tree->mmod->c_lnL = cur_glnL;

      new_glnL = PHYREX_Lk(tree);
      if(Are_Equal(new_glnL,cur_glnL,1.E-3) == NO)
        {
          PhyML_Printf("\n. new_glnL: %f cur_glnL: %f",new_glnL,cur_glnL);
          Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
        }
    }
  else
    {
      /* printf("\nD Accept %f",target_disk[0]->time); */
      For(j,n_delete_disks) Free_Disk(target_disk[j]);
      tree->mcmc->acc_move[tree->mcmc->num_move_phyrex_indel_disk]++;
    }

  tree->mcmc->run_move[tree->mcmc->num_move_phyrex_indel_disk]++;

  Free(target_disk);
  Free(permut);
}
#endif

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

#ifdef PHYREX
void MCMC_PHYREX_Insert_Disk(phydbl hr, int n_insert_disks, t_tree *tree)
{
  t_dsk  *disk,**new_disk,**target_disk;
  phydbl T,t;
  phydbl cur_glnL, new_glnL;
  phydbl u,alpha,ratio;
  int i,j,n_valid_disks;
  phydbl cur_lbda, new_lbda;
  phydbl cur_rad, new_rad;
  int n_tot_disks_cur,n_tot_disks_new;
  int err;
  phydbl mean,sd;

  if(n_insert_disks == 0) return; 

  disk            = NULL;
  new_glnL        = tree->mmod->c_lnL;
  cur_glnL        = tree->mmod->c_lnL;
  cur_lbda        = tree->mmod->lbda;
  cur_rad         = tree->mmod->rad;
  n_tot_disks_cur = PHYREX_Total_Number_Of_Intervals(tree);


  if(tree->disk->next) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);

  T = PHYREX_Tree_Height(tree);

  disk = tree->disk->prev;
  n_valid_disks = 0;
  do
    {
      if(!disk->ldsk) n_valid_disks++;
      disk = disk->prev;
    }
  while(disk && disk->prev);

  target_disk = (t_dsk **)mCalloc(n_insert_disks,sizeof(t_dsk *));
  new_disk    = (t_dsk **)mCalloc(n_insert_disks,sizeof(t_dsk *));

  For(j,n_insert_disks)
    {
      t = Uni()*T;
      disk = tree->disk;
      while(disk && disk->time > t) disk = disk->prev;
      target_disk[j] = disk->next;
      
      if(disk->next == NULL) 
        {
          target_disk[j] = tree->disk;
          PhyML_Printf("\n. run: %d",tree->mcmc->run);
          PhyML_Printf("\n. t: %f T: %f n_valid: %d n_insert: %d",t,T,n_valid_disks,n_insert_disks);
          PhyML_Printf("\n. c_lnL: %f",tree->mmod->c_lnL);
          PHYREX_Print_Struct('=',tree);
          Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
        }

      new_disk[j] = PHYREX_Make_Disk_Event(tree->mmod->n_dim,tree->n_otu);
      PHYREX_Init_Disk_Event(new_disk[j],tree->mmod->n_dim,tree->mmod);
      new_disk[j]->time = t;
        
      PHYREX_Insert_Disk(new_disk[j],tree);
      
      For(i,tree->mmod->n_dim) new_disk[j]->centr->lonlat[i] = Uni()*tree->mmod->lim->lonlat[i];

      For(i,tree->mmod->n_dim) hr -= LOG(1./tree->mmod->lim->lonlat[i]);
    }

  hr -= LnChoose(n_valid_disks+n_insert_disks,n_insert_disks);
  hr += n_insert_disks * LOG(-T);
  hr -= LnFact(n_insert_disks);

  /* Adjust value of lambda */
  n_tot_disks_new = n_tot_disks_cur + n_insert_disks;
  mean = n_tot_disks_new/FABS(T);
  sd   = 0.1*n_tot_disks_new/FABS(T);
  new_lbda = Rnorm_Trunc(mean,sd,
                         tree->mmod->min_lbda,
                         mean + 2.*sd,
                         &err);

  
  mean = n_tot_disks_cur/FABS(T);
  sd   = 0.1*n_tot_disks_cur/FABS(T);
  hr += Log_Dnorm_Trunc(cur_lbda,
                        mean,sd,
                        tree->mmod->min_lbda,
                        mean + 2.*sd,
                        &err);

  mean = n_tot_disks_new/FABS(T);
  sd   = 0.1*n_tot_disks_new/FABS(T);
  hr -= Log_Dnorm_Trunc(new_lbda,
                        mean,sd,
                        tree->mmod->min_lbda,
                        mean + 2.*sd,
                        &err);

  tree->mmod->lbda = new_lbda;



  /* Adjust value of radius */
  mean = 0.9*cur_rad;
  sd   = 0.1*mean;
  new_rad = Rnorm_Trunc(mean,sd,
                        tree->mmod->min_rad,
                        tree->mmod->max_rad,
                        &err);

  
  mean = 0.9*new_rad;
  sd   = 0.1*mean;
  hr += Log_Dnorm_Trunc(cur_rad,
                        mean,sd,
                        tree->mmod->min_rad,
                        tree->mmod->max_rad,
                        &err);

  mean = 0.9*cur_rad;
  sd   = 0.1*mean;
  hr -= Log_Dnorm_Trunc(new_rad,
                        mean,sd,
                        tree->mmod->min_rad,
                        tree->mmod->max_rad,
                        &err);

  tree->mmod->rad = new_rad;





  new_glnL = PHYREX_Lk(tree);
  ratio = (new_glnL - cur_glnL);
  ratio += hr;
  
  if(new_glnL < UNLIKELY + 0.1) 
    {
      PhyML_Printf("\n== %f %f",cur_glnL,new_glnL);
      PhyML_Printf("\n== %s",Write_Tree(tree,NO));
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
    }

  ratio = EXP(ratio);
  alpha = MIN(1.,ratio);
  
  /* Always accept move */
  if(tree->mcmc->always_yes == YES && new_glnL > UNLIKELY) alpha = 1.0;

  u = Uni();
  
  /* printf("\n+ Insert new_glnL: %f [%f] hr: %f u: %f alpha: %f",new_glnL,cur_glnL,hr,u,alpha); */

  if(u > alpha) /* Reject */
    {
      tree->mmod->lbda = cur_lbda;
      tree->mmod->rad  = cur_rad;
      
      /* printf("\nI Reject"); */

      For(j,n_insert_disks)
        {
          PHYREX_Remove_Disk(new_disk[j]);      
          Free_Disk(new_disk[j]);
        }

      tree->mmod->c_lnL = cur_glnL;

      new_glnL = PHYREX_Lk(tree);
      if(Are_Equal(new_glnL,cur_glnL,1.E-3) == NO)
        {
          PhyML_Printf("\n. new_glnL: %f cur_glnL: %f",new_glnL,cur_glnL);
          Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
        }
    }
  else
    {
      /* printf("\nI Accept"); */
      /* if(indel > 0) printf("\n. Accept"); */
      tree->mcmc->acc_move[tree->mcmc->num_move_phyrex_indel_disk]++;
    }

  tree->mcmc->run_move[tree->mcmc->num_move_phyrex_indel_disk]++;
  Free(target_disk);
  Free(new_disk);
}
#endif

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

#ifdef PHYREX
void MCMC_PHYREX_Move_Disk_Centre(t_tree *tree)
{
  phydbl u,alpha,ratio;
  phydbl cur_glnL, new_glnL, hr;
  t_dsk  *disk,**target_disk,**all_disks;
  int i,j,block,n_all_disks,n_move_disks,*permut;
  phydbl max, min;
  int err;

  disk          = NULL;
  new_glnL      = tree->mmod->c_lnL;
  cur_glnL      = tree->mmod->c_lnL;
  hr            = 0.0;
  ratio         = 0.0;
  block         = 100;
  all_disks     = NULL;

  if(tree->disk->next) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  disk = tree->disk->prev;
  n_all_disks = 0;
  do
    {
      if(!n_all_disks) all_disks = (t_dsk **)mCalloc(block,sizeof(t_dsk *));
      else if(!(n_all_disks%block)) all_disks = (t_dsk **)mRealloc(all_disks,n_all_disks+block,sizeof(t_dsk *));
      all_disks[n_all_disks] = disk;
      n_all_disks++;
      disk = disk->prev;
    }
  while(disk);

  if(!n_all_disks) return;

  target_disk = (t_dsk **)mCalloc(n_all_disks,sizeof(t_dsk *));
  
  n_move_disks = Rand_Int(1,1+(int)(n_all_disks/20));
  /* n_move_disks = n_all_disks; */

  permut = Permutate(n_all_disks);

  For(i,n_move_disks)
    {
      target_disk[i] = all_disks[permut[i]];
  
      PHYREX_Store_Geo_Coord(target_disk[i]->centr);
      
      if(target_disk[i]->ldsk != NULL)
        {

           For(j,tree->mmod->n_dim) hr += Log_Dnorm_Trunc(target_disk[i]->centr->lonlat[j],
                                                          target_disk[i]->ldsk->coord->lonlat[j],
                                                          1.0*tree->mmod->rad,
                                                          0.0,
                                                          tree->mmod->lim->lonlat[j],&err);
          
          For(j,tree->mmod->n_dim)
            target_disk[i]->centr->lonlat[j] =
            Rnorm_Trunc(target_disk[i]->ldsk->coord->lonlat[j],
                        1.0*tree->mmod->rad,
                        0.0,
                        tree->mmod->lim->lonlat[j],&err);
          
          For(j,tree->mmod->n_dim) hr -= Log_Dnorm_Trunc(target_disk[i]->centr->lonlat[j],
                                                         target_disk[i]->ldsk->coord->lonlat[j],
                                                         1.0*tree->mmod->rad,
                                                         0.0,
                                                         tree->mmod->lim->lonlat[j],&err);

        }
      else
        {
          For(j,tree->mmod->n_dim)
            {
              max = tree->mmod->lim->lonlat[j];
              min = 0.0;
              target_disk[i]->centr->lonlat[j] = Uni()*(max - min) + min;
            }
        }
    }

  Free(permut);

  new_glnL = PHYREX_Lk(tree);
  ratio += (new_glnL - cur_glnL);
  ratio += hr;
  
  if(new_glnL < UNLIKELY + 0.1) 
    {
      PhyML_Printf("\n== %f %f",cur_glnL,new_glnL);
      PhyML_Printf("\n== %s",Write_Tree(tree,NO));
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
    }

  ratio = EXP(ratio);
  alpha = MIN(1.,ratio);
  
  /* Always accept move */
  if(tree->mcmc->always_yes == YES && new_glnL > UNLIKELY) alpha = 1.0;

  u = Uni();
  
  /* printf("\n- Delete new_glnL: %f [%f] hr: %f u:%f alpha: %f",new_glnL,cur_glnL,hr,u,alpha); */

  if(u > alpha) /* Reject */
    {
      /* printf("- Reject"); */
      
      For(i,n_move_disks) PHYREX_Restore_Geo_Coord(target_disk[i]->centr);
      
      tree->mmod->c_lnL = cur_glnL;

      new_glnL = PHYREX_Lk(tree);
      if(Are_Equal(new_glnL,cur_glnL,1.E-3) == NO)
        {
          PhyML_Printf("\n. new_glnL: %f cur_glnL: %f",new_glnL,cur_glnL);
          Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
        }
    }
  else
    {
      tree->mcmc->acc_move[tree->mcmc->num_move_phyrex_move_disk_ct]++;
    }

  tree->mcmc->run_move[tree->mcmc->num_move_phyrex_move_disk_ct]++;

  Free(all_disks);
  Free(target_disk);
}
#endif

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
/* Move ldsk on landscape */ 
#ifdef PHYREX
void MCMC_PHYREX_Move_Ldsk(t_tree *tree)
{
  phydbl u,alpha,ratio;
  phydbl cur_glnL, new_glnL, hr;
  t_dsk  *disk,**target_disk,**all_disks;
  int i,j,block,n_all_disks,n_move_ldsk,*permut;
  phydbl max,min;

  disk          = NULL;
  new_glnL      = tree->mmod->c_lnL;
  cur_glnL      = tree->mmod->c_lnL;
  hr            = 0.0;
  ratio         = 0.0;
  block         = 100;
  all_disks     = NULL;

  if(tree->disk->next) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  disk = tree->disk->prev;
  n_all_disks = 0;
  do
    {
      if(disk->ldsk != NULL)
        {
          if(!n_all_disks) all_disks = (t_dsk **)mCalloc(block,sizeof(t_dsk *));
          else if(!(n_all_disks%block)) all_disks = (t_dsk **)mRealloc(all_disks,n_all_disks+block,sizeof(t_dsk *));
          all_disks[n_all_disks] = disk;
          n_all_disks++;
        }
      disk = disk->prev;
    }
  while(disk);

  if(!n_all_disks) return;
  
  n_move_ldsk = Rand_Int(1,1+(int)(n_all_disks/20));
  /* n_move_ldsk = n_all_disks; */
  
  target_disk = (t_dsk **)mCalloc(n_all_disks,sizeof(t_dsk *));

  permut = Permutate(n_all_disks);

  For(i,n_move_ldsk)
    {
      target_disk[i] = all_disks[permut[i]];
  
      PHYREX_Store_Geo_Coord(target_disk[i]->ldsk->coord);
      PHYREX_Store_Geo_Coord(target_disk[i]->centr);

      max = min = -1.;
      if(tree->mmod->name == PHYREX_UNIFORM)
        {
          For(j,tree->mmod->n_dim)
            {
              max = target_disk[i]->centr->lonlat[j] + tree->mmod->rad;
              min = target_disk[i]->centr->lonlat[j] - tree->mmod->rad;
              /* need to add condition here that new position also has to be into disk above */
              max = MIN(tree->mmod->lim->lonlat[j],max);
              min = MAX(0.0,min);
              target_disk[i]->ldsk->coord->lonlat[j] = Uni()*(max - min) + min;
            }
        }
      else if(tree->mmod->name == PHYREX_NORMAL)
        {
          int err;

          For(j,tree->mmod->n_dim) hr += Log_Dnorm_Trunc(target_disk[i]->ldsk->coord->lonlat[j],
                                                         target_disk[i]->centr->lonlat[j],
                                                         1.0*tree->mmod->rad,
                                                         0.0,
                                                         tree->mmod->lim->lonlat[j],&err);
          

          For(j,tree->mmod->n_dim)
            target_disk[i]->ldsk->coord->lonlat[j] =
            Rnorm_Trunc(target_disk[i]->centr->lonlat[j],
                        1.0*tree->mmod->rad,
                        0.0,
                        tree->mmod->lim->lonlat[j],&err);

          For(j,tree->mmod->n_dim) hr -= Log_Dnorm_Trunc(target_disk[i]->ldsk->coord->lonlat[j],
                                                         target_disk[i]->centr->lonlat[j],
                                                         1.0*tree->mmod->rad,
                                                         0.0,
                                                         tree->mmod->lim->lonlat[j],&err);

        }
    }

  Free(permut);

  new_glnL = PHYREX_Lk(tree);
  ratio += (new_glnL - cur_glnL);
  ratio += hr;
  
  if(new_glnL < UNLIKELY + 0.1) 
    {
      PhyML_Printf("\n== %f %f",cur_glnL,new_glnL);
      PhyML_Printf("\n== %s",Write_Tree(tree,NO));
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
    }

  ratio = EXP(ratio);
  alpha = MIN(1.,ratio);
  
  /* Always accept move */
  if(tree->mcmc->always_yes == YES && new_glnL > UNLIKELY) alpha = 1.0;

  u = Uni();
  
  /* printf("\n- Move_ldsk %15f",new_glnL-cur_glnL); */

  if(u > alpha) /* Reject */
    {
      /* printf("- Reject"); */
      
      For(i,n_move_ldsk) 
        {
          PHYREX_Restore_Geo_Coord(target_disk[i]->ldsk->coord);
          PHYREX_Restore_Geo_Coord(target_disk[i]->centr);
        }

      tree->mmod->c_lnL = cur_glnL;

      new_glnL = PHYREX_Lk(tree);
      if(Are_Equal(new_glnL,cur_glnL,1.E-3) == NO)
        {
          PhyML_Printf("\n. new_glnL: %f cur_glnL: %f",new_glnL,cur_glnL);
          Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
        }
    }
  else
    {
      tree->mcmc->acc_move[tree->mcmc->num_move_phyrex_move_ldsk]++;
    }

  tree->mcmc->run_move[tree->mcmc->num_move_phyrex_move_ldsk]++;
   
  Free(all_disks);
  Free(target_disk);
}
#endif

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
/* Update time of disk */
#ifdef PHYREX
void MCMC_PHYREX_Move_Disk_Updown(t_tree *tree)
{
  phydbl u,alpha,ratio;
  phydbl cur_glnL, new_glnL, hr;
  phydbl cur_alnL, new_alnL;
  phydbl *ori_time, new_time;
  phydbl max, min;
  t_dsk  *disk,**target_disk,**all_disks;
  int i,block,n_all_disks,n_move_disks,*permut;

  disk          = NULL;
  new_alnL      = tree->c_lnL;
  cur_alnL      = tree->c_lnL;
  new_glnL      = tree->mmod->c_lnL;
  cur_glnL      = tree->mmod->c_lnL;
  hr            = 0.0;
  ratio         = 0.0;
  block         = 100;
  all_disks     = NULL;


  if(tree->disk->next) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  disk = tree->disk->prev;
  n_all_disks = 0;
  do
    {
      if(1)
      /* if(disk->ldsk && disk->ldsk->n_next > 1) /\* Moving disks other than coalescent is pointless *\/  */
        {
          if(!n_all_disks) all_disks = (t_dsk **)mCalloc(block,sizeof(t_dsk *));
          else if(!(n_all_disks%block)) all_disks = (t_dsk **)mRealloc(all_disks,n_all_disks+block,sizeof(t_dsk *));
          all_disks[n_all_disks] = disk;
          n_all_disks++;
        }
      disk = disk->prev;
    }
  while(disk);

  if(!n_all_disks) return;
  
  /* n_move_disks = Rand_Int(1,1+(int)(n_all_disks/5)); */
  n_move_disks = n_all_disks;
  /* n_move_disks = 1; */
  
  target_disk = (t_dsk **)mCalloc(n_all_disks,sizeof(t_dsk *));
  ori_time    = (phydbl *)mCalloc(n_all_disks,sizeof(phydbl));

  permut = Permutate(n_all_disks);

  For(i,n_move_disks)
    {
      target_disk[i] = all_disks[permut[i]];
      
      ori_time[i] = target_disk[i]->time;

      if(target_disk[i]->prev)
        {
          max = target_disk[i]->next->time;
          min = target_disk[i]->prev->time;
          new_time = Uni()*(max - min) + min;
        }
      else
        {
          
          /* phydbl scale,K; */

          /* K = 1.; */
          /* scale = exp(K*(Uni()-0.5)); */
          
          /* max = target_disk[i]->next->time; */
          /* new_time = scale * ori_time[i] + max*(1.-scale); */
          /* hr += LOG(K*((new_time - max)/(ori_time[i] - max)) - max); */
          /* hr -= LOG(K*((ori_time[i] - max)/(new_time - max)) - max); */


          phydbl new_plusmin, cur_plusmin;

          max = target_disk[i]->next->time;

          cur_plusmin = FABS(ori_time[i] - max);
          new_plusmin = Rexp(1./cur_plusmin);

          new_time = max - new_plusmin;

          hr += LOG(Dexp(cur_plusmin,1./new_plusmin));
          hr -= LOG(Dexp(new_plusmin,1./cur_plusmin));

        }
            
      target_disk[i]->time = new_time;
    }
  
  Free(permut);

  new_glnL = PHYREX_Lk(tree);
  if(tree->mcmc->use_data == YES) new_alnL = Lk(NULL,tree);

  ratio += (new_alnL - cur_alnL);
  ratio += (new_glnL - cur_glnL);
  ratio += hr;
  
  if(new_glnL < UNLIKELY + 0.1) 
    {
      PhyML_Printf("\n== %f %f",cur_glnL,new_glnL);
      PhyML_Printf("\n== %s",Write_Tree(tree,NO));
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
    }

  ratio = EXP(ratio);
  alpha = MIN(1.,ratio);
  
  /* Always accept move */
  if(tree->mcmc->always_yes == YES && new_glnL > UNLIKELY) alpha = 1.0;

  u = Uni();
  
  /* printf("\n- Move disk new_glnL: %f [%f] hr: %f u:%f alpha: %f",new_glnL,cur_glnL,hr,u,alpha); */

  if(u > alpha) /* Reject */
    {
      /* printf("- Reject"); */
      
      For(i,n_move_disks) target_disk[i]->time = ori_time[i]; 

      tree->mmod->c_lnL = cur_glnL;
      tree->c_lnL       = cur_alnL;

      /* if(tree->mcmc->use_data == YES) Lk(NULL,tree); */

      new_glnL = PHYREX_Lk(tree);

      if(Are_Equal(new_glnL,cur_glnL,1.E-3) == NO)
        {
          PhyML_Printf("\n. new_glnL: %f cur_glnL: %f",new_glnL,cur_glnL);
          Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
        }
    }
  else
    {
      /* printf("- Accept"); */
    }  

  Free(target_disk);
  Free(all_disks);
  Free(ori_time);

}
#endif

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
#ifdef PHYREX
void MCMC_PHYREX_Scale_Times(t_tree *tree)
{
  phydbl u,alpha,ratio;
  phydbl cur_glnL, new_glnL, hr;
  phydbl cur_alnL, new_alnL;
  phydbl scale_fact_times;
  int n_disks;
  t_dsk  *disk;
  phydbl K;
  phydbl cur_lbda, new_lbda;

  disk            = NULL;
  new_alnL        = tree->c_lnL;
  cur_alnL        = tree->c_lnL;
  new_glnL        = tree->mmod->c_lnL;
  cur_glnL        = tree->mmod->c_lnL;
  hr              = 0.0;
  ratio           = 0.0;
  K               = tree->mcmc->tune_move[tree->mcmc->num_move_phyrex_scale_times];
  cur_lbda        = tree->mmod->lbda;

  u = Uni();
  scale_fact_times = EXP(K*(u-.5));
  
  if(tree->disk->next) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);

  n_disks = 0;
  disk = tree->disk->prev;
  do
    {
      disk->time = disk->time * scale_fact_times;
      n_disks++;
      disk = disk->prev;
    }
  while(disk);

  
  /* The Hastings ratio involves (n_disk-2) when considering a uniform distrib
     for the multiplier, which is not the case here.
  */
  hr += (n_disks)*LOG(scale_fact_times);

  /* Adjust the value of lambda */
  new_lbda = cur_lbda * (1./scale_fact_times);
  hr += LOG(1./scale_fact_times);
  tree->mmod->lbda = new_lbda;

  new_glnL = PHYREX_Lk(tree);
  if(tree->mcmc->use_data == YES) new_alnL = Lk(NULL,tree);

  ratio += (new_alnL - cur_alnL);
  ratio += (new_glnL - cur_glnL);
  ratio += hr;
  
  if(new_glnL < UNLIKELY + 0.1) 
    {
      PhyML_Printf("\n== %f %f",cur_glnL,new_glnL);
      PhyML_Printf("\n== %s",Write_Tree(tree,NO));
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
    }

  ratio = EXP(ratio);
  alpha = MIN(1.,ratio);
  
  /* Always accept move */
  if(tree->mcmc->always_yes == YES && new_glnL > UNLIKELY) alpha = 1.0;

  u = Uni();
  
  /* PhyML_Printf("\n. Scale times hr: %f new_glnL: %f cur_glnL: %f ratio: %f", */
  /*              hr, */
  /*              new_glnL,cur_glnL, */
  /*              ratio) */

  if(u > alpha) /* Reject */
    {
      tree->mmod->lbda = cur_lbda;

      /* printf("- Reject"); */      
      disk = tree->disk->prev;
      do
        {
          disk->time /= scale_fact_times;
          disk = disk->prev;
        }
      while(disk);

      tree->mmod->c_lnL = cur_glnL;
      tree->c_lnL       = cur_alnL;

      /* if(tree->mcmc->use_data == YES) Lk(NULL,tree); /\* Not necessary. Remove once tested *\/ */

      new_glnL = PHYREX_Lk(tree); /* Same here */
      
      if(Are_Equal(new_glnL,cur_glnL,1.E-3) == NO)
        {
          PhyML_Printf("\n. new_glnL: %f cur_glnL: %f",new_glnL,cur_glnL);
          Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
        }
    }
  else
    {
      /* printf("- Accept"); */
      tree->mcmc->acc_move[tree->mcmc->num_move_phyrex_scale_times]++;
    }  
  tree->mcmc->run_move[tree->mcmc->num_move_phyrex_scale_times]++;
}
#endif

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

#ifdef PHYREX
void MCMC_PHYREX_Swap_Disk(t_tree *tree)
{
  phydbl u,alpha,ratio;
  phydbl cur_glnL, new_glnL;
  phydbl cur_alnL, new_alnL;
  phydbl hr,t,t_min,t_max,ori_time;
  t_dsk  *disk,*target_disk,*ori_disk_old,*ori_disk_young,**valid_disks;
  int i,j,block,n_valid_disks;

  For(j,tree->n_otu)
    {
      disk           = NULL;
      target_disk    = NULL;
      ori_disk_old   = NULL;
      ori_disk_young = NULL;
      valid_disks    = NULL;
      new_glnL       = tree->mmod->c_lnL;
      cur_glnL       = tree->mmod->c_lnL;
      new_alnL       = tree->c_lnL;
      cur_alnL       = tree->c_lnL;
      hr             = 0.0;
      ratio          = 0.0;
      block          = 100;
      
      if(tree->disk->next) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
      disk = tree->disk->prev;
      if(!disk->prev) return;
      n_valid_disks = 0;
      do
        {
          /* Record disk with a lineage displacement or coalescent that is not the root disk */
          if(disk && disk->prev && disk->ldsk && disk->ldsk->n_next >= 1)
            {
              if(!n_valid_disks) valid_disks = (t_dsk **)mCalloc(block,sizeof(t_dsk *));
              else if(!(n_valid_disks%block)) valid_disks = (t_dsk **)mRealloc(valid_disks,n_valid_disks+block,sizeof(t_dsk *));
              valid_disks[n_valid_disks] = disk;
              n_valid_disks++;
            }
          disk = disk->prev;
        }
      while(disk && disk->prev);
      
      if(n_valid_disks < 2) return;
      
      /* Uniform selection of a disk where either a coalescent or a 'hit/displacement' occurred */
      i = Rand_Int(0,n_valid_disks-1);
      target_disk = valid_disks[i];
      Free(valid_disks);
      
      ori_disk_old   = target_disk->prev;
      ori_disk_young = target_disk->next;
      ori_time       = target_disk->time;
      
      t_min = target_disk->ldsk->prev->disk->time;
      t_max = +INFINITY;
      For(i,target_disk->ldsk->n_next)
        if(target_disk->ldsk->next[i]->disk->time < t_max) 
          t_max = target_disk->ldsk->next[i]->disk->time;
      
      t_max = t_max - 1.E-10;
      t_min = t_min + 1.E-10;
      
      if(t_max < t_min) return;
      
      t = Uni()*(t_max - t_min) + t_min;
      target_disk->time = t;
      
      PHYREX_Remove_Disk(target_disk);
      
      PHYREX_Insert_Disk(target_disk,tree);
      
      new_glnL = PHYREX_Lk(tree);
      if(tree->mcmc->use_data == YES) new_alnL = Lk(NULL,tree);
      
      ratio += (new_alnL - cur_alnL);
      ratio += (new_glnL - cur_glnL);
      ratio += hr;
      
      if(new_glnL < UNLIKELY + 0.1) 
        {
          PhyML_Printf("\n== %f %f",cur_glnL,new_glnL);
          PhyML_Printf("\n== %s",Write_Tree(tree,NO));
          Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
        }
      
      ratio = EXP(ratio);
      alpha = MIN(1.,ratio);
      
      /* Always accept move */
      if(tree->mcmc->always_yes == YES && new_glnL > UNLIKELY) alpha = 1.0;
      
      u = Uni();
      
      /* printf("\n- Swap new_glnL: %f [%f] hr: %f u:%f alpha: %f",new_glnL,cur_glnL,hr,u,alpha); */
      
      if(u > alpha) /* Reject */
        {
          /* printf("\n- Reject %f %f",target_disk->time,ori_time); */
          
          PHYREX_Remove_Disk(target_disk);
          
          target_disk->time = ori_time;
          target_disk->prev = ori_disk_old;
          target_disk->next = ori_disk_young;
          
          PHYREX_Insert_Disk(target_disk,tree);
          
          tree->mmod->c_lnL = cur_glnL;
          tree->c_lnL       = cur_alnL;
          
          /* if(tree->mcmc->use_data == YES) Lk(NULL,tree); /\* Remove once checked *\/ */
          
          new_glnL = PHYREX_Lk(tree); /* Same */
          if(Are_Equal(new_glnL,cur_glnL,1.E-3) == NO)
            {
              PhyML_Printf("\n. new_glnL: %f cur_glnL: %f",new_glnL,cur_glnL);
              Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
            }
        }
      else
        {
          /* printf("\n- Accept %f %f",target_disk->time,ori_time); */
        }
    }
}
#endif
  
/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
#ifdef PHYREX
void MCMC_PHYREX_Indel_Hit(t_tree *tree)
{
  int n_disks_cur, n_disks_new;
  t_dsk  *disk;
  phydbl hr,K;
  
  K  = 1.0;
  hr = 0.0;

  disk = tree->disk->prev;
  n_disks_cur = 0;
  do
    {
      if(disk->ldsk != NULL && disk->ldsk->n_next == 1) n_disks_cur++;
      disk = disk->prev;
    }
  while(disk && disk->prev);
  

  /* n_disks_new = (int)Rpois(K*n_disks_cur); */
  
  /* hr += Dpois(n_disks_cur,K*n_disks_new,YES); */
  /* hr -= Dpois(n_disks_new,K*n_disks_cur,YES); */


  phydbl e,v;

  e = (phydbl)n_disks_cur;
  v = 1.+e*5.;

  n_disks_new = (int)Rnbinom(e*e/(v-e),e/v);
  
  hr += Dnbinom(n_disks_cur,e*e/(v-e),e/v,YES);
  hr -= Dnbinom(n_disks_new,e*e/(v-e),e/v,YES);


  if(n_disks_new < n_disks_cur) MCMC_PHYREX_Delete_Hit(hr, n_disks_cur - n_disks_new ,tree);
  else                          MCMC_PHYREX_Insert_Hit(hr, n_disks_new - n_disks_cur,tree);  
}
#endif


/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
/* Insert a disk with a new lineage displacement */
#ifdef PHYREX
void MCMC_PHYREX_Insert_Hit(phydbl hr, int n_insert_disks, t_tree *tree)
{
  t_dsk  *disk,**new_disk,*young_disk;
  t_ldsk **young_ldsk, **old_ldsk, **new_ldsk;
  phydbl T,t;
  phydbl cur_glnL, new_glnL;
  phydbl u,alpha,ratio;
  int i,j,*dir_old_young,n_valid_disks,err;
  /* phydbl cur_lbda, new_lbda; */
  /* int n_tot_disks_new, n_tot_disks_cur; */

  disk            = NULL;
  new_glnL        = tree->mmod->c_lnL;
  cur_glnL        = tree->mmod->c_lnL;
  /* cur_lbda        = tree->mmod->lbda; */
  /* n_tot_disks_cur = PHYREX_Total_Number_Of_Intervals(tree); */

  if(tree->disk->next) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);

  disk = tree->disk->prev;
  n_valid_disks = 0;
  do
    {
      if(disk->ldsk != NULL && disk->ldsk->n_next == 1) n_valid_disks++;
      disk = disk->prev;
    }
  while(disk && disk->prev);
  
  if(n_insert_disks == 0)  { return; }

  new_disk      = (t_dsk **)mCalloc(n_insert_disks,sizeof(t_dsk *));
  new_ldsk      = (t_ldsk **)mCalloc(n_insert_disks,sizeof(t_ldsk *));
  old_ldsk      = (t_ldsk **)mCalloc(n_insert_disks,sizeof(t_ldsk *));
  young_ldsk    = (t_ldsk **)mCalloc(n_insert_disks,sizeof(t_ldsk *));
  dir_old_young = (int *)mCalloc(n_insert_disks,sizeof(int));

  T = PHYREX_Tree_Height(tree);

  For(j,n_insert_disks)
    {  
      /* Time of insertion of new disk */
      t = Uni()*T;
      disk = tree->disk;
      while(disk && disk->time > t) disk = disk->prev;
      
      /* Disks located just prior and after inserted disk */
      young_disk = disk->next;

      /* Make and initialize new disk */ 
      new_disk[j] = PHYREX_Make_Disk_Event(tree->mmod->n_dim,tree->n_otu);
      PHYREX_Init_Disk_Event(new_disk[j],tree->mmod->n_dim,tree->mmod);
        
      /* Time of the new disk */
      new_disk[j]->time = t;
      
      /* Insert it */
      PHYREX_Insert_Disk(new_disk[j],tree);

      if(!young_disk->n_ldsk_a) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);

      /* Which lineage is going to be hit ? */
      hr += LOG(young_disk->n_ldsk_a);
      
      young_ldsk[j] = young_disk->ldsk_a[Rand_Int(0,young_disk->n_ldsk_a-1)];
      old_ldsk[j]   = young_ldsk[j]->prev;
  
      if(old_ldsk[j]->disk->time > young_ldsk[j]->disk->time) 
        {
          PhyML_Printf("\n. young_ldsk: %f",young_ldsk[j]->disk->time);
          PhyML_Printf("\n. old_ldsk: %f",old_ldsk[j]->disk->time);
          Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
        }
      
      /* Direction from old to young ldsk */
      dir_old_young[j] = PHYREX_Get_Next_Direction(young_ldsk[j],old_ldsk[j]);
  
      if(dir_old_young[j] == -1) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);

      /* Make and initialize new ldsk */
      new_ldsk[j] = PHYREX_Make_Lindisk_Node(tree->mmod->n_dim);
      PHYREX_Init_Lindisk_Node(new_ldsk[j],new_disk[j],tree->mmod->n_dim);
      PHYREX_Make_Lindisk_Next(new_ldsk[j]);
      new_disk[j]->ldsk = new_ldsk[j];

      /* Connect it */
      new_ldsk[j]->prev                   = old_ldsk[j];
      new_ldsk[j]->next[0]                = young_ldsk[j];
      young_ldsk[j]->prev                 = new_ldsk[j];  
      old_ldsk[j]->next[dir_old_young[j]] = new_ldsk[j];

      /* PHYREX_Update_Lindisk_List_Core(young_disk); */
      /* PHYREX_Update_Lindisk_List_Core(new_disk[j]); */
      /* PHYREX_Update_Lindisk_List_Core(old_disk); */

      PHYREX_Update_Lindisk_List(tree);
      
      /* Sample position of the displaced ldsk */
      For(i,tree->mmod->n_dim)
        {
          new_ldsk[j]->coord->lonlat[i] = Rnorm_Trunc(0.5*(young_ldsk[j]->coord->lonlat[i]+old_ldsk[j]->coord->lonlat[i]),
                                                      SQRT(2.0*POW(tree->mmod->rad,2)),
                                                      0.0,tree->mmod->lim->lonlat[i],&err);
          hr -= Log_Dnorm_Trunc(new_ldsk[j]->coord->lonlat[i],
                                0.5*(young_ldsk[j]->coord->lonlat[i]+old_ldsk[j]->coord->lonlat[i]),
                                SQRT(2.0*POW(tree->mmod->rad,2)),
                                0.0,tree->mmod->lim->lonlat[i],&err);
        }

      /* Sample position of the center of new_disk */
      For(i,tree->mmod->n_dim)
        {
          new_disk[j]->centr->lonlat[i] = Rnorm_Trunc(new_ldsk[j]->coord->lonlat[i],
                                                      1.0*tree->mmod->rad,
                                                      0.0,tree->mmod->lim->lonlat[i],&err);
          hr -= Log_Dnorm_Trunc(new_disk[j]->centr->lonlat[i],
                                new_ldsk[j]->coord->lonlat[i],
                                1.0*tree->mmod->rad,0.0,tree->mmod->lim->lonlat[i],&err);

        }

    }

  T = PHYREX_Tree_Height(tree);

  hr -= LnChoose(n_valid_disks+n_insert_disks,n_insert_disks);
  hr += n_insert_disks * LOG(-T);  
  hr -= LnFact(n_insert_disks);

  new_glnL = PHYREX_Lk(tree);
  
  if(new_glnL < UNLIKELY + 0.1) 
    {
      PhyML_Printf("\n== %f %f",cur_glnL,new_glnL);
      PhyML_Printf("\n== %s",Write_Tree(tree,NO));
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
    }

  ratio = (new_glnL - cur_glnL);
  ratio += hr;  

  ratio = EXP(ratio);
  alpha = MIN(1.,ratio);

  if(tree->mcmc->always_yes == YES && new_glnL > UNLIKELY) alpha = 1.0;

  u = Uni();
  
  /* printf("\n- Insert hit %15f %15f %5d",new_glnL - cur_glnL, alpha, n_insert_disks); */

  if(u > alpha) /* Reject */
    {
      /* printf("+ Reject"); */

      /* tree->mmod->lbda = cur_lbda; */
      
      for(j=n_insert_disks-1;j>=0;j--)
        {
          PHYREX_Remove_Disk(new_disk[j]);      
          Free_Disk(new_disk[j]);
          Free_Ldisk(new_ldsk[j]);

          old_ldsk[j]->next[dir_old_young[j]] = young_ldsk[j];
          young_ldsk[j]->prev = old_ldsk[j];
        }
      
      tree->mmod->c_lnL = cur_glnL;

      new_glnL = PHYREX_Lk(tree);
      if(Are_Equal(new_glnL,cur_glnL,1.E-3) == NO)
        {
          PhyML_Printf("\n. new_glnL: %f cur_glnL: %f",new_glnL,cur_glnL);
          Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
        }
    }
  else
    {
      /* printf("+ Accept"); */
      /* if(indel > 0) printf("\n. Accept"); */
      tree->mcmc->acc_move[tree->mcmc->num_move_phyrex_indel_hit]++;
    }

  tree->mcmc->run_move[tree->mcmc->num_move_phyrex_indel_hit]++;

  Free(new_disk);
  Free(new_ldsk);
  Free(old_ldsk);
  Free(young_ldsk);
  Free(dir_old_young);


}
#endif

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
/* Remove one or more disks with a lineage displacement */
#ifdef PHYREX
void MCMC_PHYREX_Delete_Hit(phydbl hr, int n_delete_disks, t_tree *tree)
{
  phydbl u,alpha,ratio;
  phydbl cur_glnL,new_glnL,T;
  t_dsk  *disk,**target_disk,**valid_disks;
  t_ldsk **target_ldsk,**old_ldsk,**young_ldsk;
  int i,j,block,n_valid_disks,*dir_old_young,*permut,err;
  /* phydbl cur_lbda, new_lbda; */
  /* int n_tot_disks_new, n_tot_disks_cur; */

  disk            = NULL;
  valid_disks     = NULL;
  new_glnL        = tree->mmod->c_lnL;
  cur_glnL        = tree->mmod->c_lnL;  
  ratio           = 0.0;
  block           = 100;
  /* cur_lbda        = tree->mmod->lbda; */
  /* n_tot_disks_cur = PHYREX_Total_Number_Of_Intervals(tree); */

  if(tree->disk->next) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  disk = tree->disk->prev;
  if(!disk->prev) return;

  n_valid_disks = 0;
  do
    {
      /* Include only disks with displacement that are not coalescent events */
      if(disk->ldsk != NULL && disk->ldsk->n_next == 1) 
        {
          if(!n_valid_disks) valid_disks = (t_dsk **)mCalloc(block,sizeof(t_dsk *));
          else if(!(n_valid_disks%block)) valid_disks = (t_dsk **)mRealloc(valid_disks,n_valid_disks+block,sizeof(t_dsk *));
          valid_disks[n_valid_disks] = disk;
          n_valid_disks++;
        }
      disk = disk->prev;
    }
  while(disk && disk->prev);

  if(!n_valid_disks) return;
  
  if(n_valid_disks - n_delete_disks < 0) { return; }

  target_disk   = (t_dsk **)mCalloc(n_delete_disks,sizeof(t_dsk *));
  target_ldsk   = (t_ldsk **)mCalloc(n_delete_disks,sizeof(t_ldsk *));
  old_ldsk      = (t_ldsk **)mCalloc(n_delete_disks,sizeof(t_ldsk *));
  young_ldsk    = (t_ldsk **)mCalloc(n_delete_disks,sizeof(t_ldsk *));
  dir_old_young = (int *)mCalloc(n_delete_disks,sizeof(int));

  permut = Permutate(n_valid_disks);

  For(j,n_delete_disks)
    {
      target_disk[j] = valid_disks[permut[j]];  
      target_ldsk[j] = target_disk[j]->ldsk;
  
      if(!target_ldsk[j])             Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
      if(target_ldsk[j]->n_next != 1) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);      
      
      old_ldsk[j]   = target_ldsk[j]->prev;
      young_ldsk[j] = target_ldsk[j]->next[0];
            
      dir_old_young[j] = PHYREX_Get_Next_Direction(young_ldsk[j],old_ldsk[j]);
      
      if(dir_old_young[j] == -1) 
        {
          PhyML_Printf("\n. j=%d n_delete_disks=%d",j,n_delete_disks);
          Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
        }
      
      /* Part of the Hastings ratio corresponding to the probability of selecting */
      /* one of target_disk->n_ldsk_a to be hit */ 
      hr -= LOG(target_disk[j]->next->n_ldsk_a);
            
      /* Density for position of the displaced ldsk */
      For(i,tree->mmod->n_dim)
        {
          hr += Log_Dnorm_Trunc(target_ldsk[j]->coord->lonlat[i],
                                0.5*(young_ldsk[j]->coord->lonlat[i]+old_ldsk[j]->coord->lonlat[i]),
                                SQRT(2.0*POW(tree->mmod->rad,2)),
                                0.0,
                                tree->mmod->lim->lonlat[i],&err);
        }

      /* Density for position of the center of target_disk */
      For(i,tree->mmod->n_dim)
        {
          hr += Log_Dnorm_Trunc(target_disk[j]->centr->lonlat[i],
                                target_ldsk[j]->coord->lonlat[i],
                                1.0*tree->mmod->rad,
                                0.0,tree->mmod->lim->lonlat[i],&err);
        }

      /* New connections between old_ldsk and young_ldsk */
      old_ldsk[j]->next[dir_old_young[j]] = young_ldsk[j];
      young_ldsk[j]->prev                 = old_ldsk[j];

      /* Remove target disk */
      PHYREX_Remove_Disk(target_disk[j]);
    }
  
  T = PHYREX_Tree_Height(tree);
  
  hr += LnChoose(n_valid_disks,n_delete_disks);
  hr -= n_delete_disks * LOG(-T);
  hr += LnFact(n_delete_disks);

  Free(valid_disks);

  new_glnL = PHYREX_Lk(tree);
  ratio += (new_glnL - cur_glnL);
  ratio += hr;

  if(new_glnL < UNLIKELY + 0.1) 
    {
      PhyML_Printf("\n== %f %f",cur_glnL,new_glnL);
      PhyML_Printf("\n== %s",Write_Tree(tree,NO));
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
    }

  ratio = EXP(ratio);
  alpha = MIN(1.,ratio);
  
  /* Always accept move */
  if(tree->mcmc->always_yes == YES && new_glnL > UNLIKELY) alpha = 1.0;

  u = Uni();
  
  /* printf("\n- Delete hit %15f %15f %5d",new_glnL - cur_glnL, alpha, n_delete_disks); */

  if(u > alpha) /* Reject */
    {
      /* tree->mmod->lbda = cur_lbda; */

      /* printf("- Reject"); */
      for(j=n_delete_disks-1;j>=0;j--) 
        {
          PHYREX_Insert_Disk(target_disk[j],tree);      
          old_ldsk[j]->next[dir_old_young[j]] = target_ldsk[j];
          young_ldsk[j]->prev                 = target_ldsk[j];
        }
      
      tree->mmod->c_lnL = cur_glnL;

      new_glnL = PHYREX_Lk(tree);
      if(Are_Equal(new_glnL,cur_glnL,1.E-3) == NO)
        {
          PhyML_Printf("\n. new_glnL: %f cur_glnL: %f",new_glnL,cur_glnL);
          Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
        }
    }
  else
    {
      /* printf(" ***"); */

      For(j,n_delete_disks) Free_Disk(target_disk[j]);
      For(j,n_delete_disks) Free_Ldisk(target_ldsk[j]);
      tree->mcmc->acc_move[tree->mcmc->num_move_phyrex_indel_hit]++;
    }
  
  tree->mcmc->run_move[tree->mcmc->num_move_phyrex_indel_hit]++;

  Free(target_disk);
  Free(target_ldsk);
  Free(old_ldsk);
  Free(young_ldsk);
  Free(dir_old_young);
  Free(permut);
}
#endif

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

#ifdef PHYREX
void MCMC_PHYREX_Prune_Regraft(t_tree *tree)
{
  phydbl u,alpha,ratio,hr;
  phydbl cur_glnL, new_glnL;
  phydbl cur_alnL, new_alnL;
  t_dsk  *disk,*prune_disk,*regraft_disk,**valid_disks;
  t_ldsk *prune_ldsk,*regraft_ldsk,*prune_daughter_ldsk,*cur_path,*new_path,*ldsk,*ldsk_dum;
  int i,block,n_valid_disks,prune_next_num;
  phydbl rate,dt,sizeT;
  int cur_path_len;
  int n_hits,n_iter;

  n_iter = MAX(1,(int)(tree->n_otu/5));

  while(n_iter--)
    {

      valid_disks = NULL;
      disk        = NULL;
      new_glnL    = tree->mmod->c_lnL;
      cur_glnL    = tree->mmod->c_lnL;
      new_alnL    = tree->c_lnL;
      cur_alnL    = tree->c_lnL;
      hr          = 0.0;
      ratio       = 0.0;
      block       = 100;
      
      if(tree->disk->next) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
      
      /* Get a ldsk from which you can prune a lineage */
      disk = tree->disk->prev;
      n_valid_disks = 0;
      do
        {
          /* Include only disks with displacement that are coalescent events */
          if(disk->ldsk != NULL && disk->ldsk->n_next > 1)
            {
              /* Root has degree 2: don't pull any lineage */
              if(!(disk->prev == NULL && disk->ldsk->n_next == 2))
                {
                  if(!n_valid_disks) valid_disks = (t_dsk **)mCalloc(block,sizeof(t_dsk *));
                  else if(!(n_valid_disks%block)) valid_disks = (t_dsk **)mRealloc(valid_disks,n_valid_disks+block,sizeof(t_dsk *));
                  valid_disks[n_valid_disks] = disk;
                  n_valid_disks++;
                }
            }
          disk = disk->prev;
        }
      while(disk);
      
      if(!n_valid_disks) return;
      
      hr -= LOG(1./n_valid_disks);
      
      /* Uniform selection of a disk where a coalescent occurred */
      i = Rand_Int(0,n_valid_disks-1);
      prune_disk = valid_disks[i];
      Free(valid_disks);
      
      prune_ldsk = prune_disk->ldsk;
      /* Which daughter lineage are we pruning? */
      prune_next_num = Rand_Int(0,prune_ldsk->n_next-1);
      prune_daughter_ldsk = prune_ldsk->next[prune_next_num];
      while(prune_daughter_ldsk->n_next < 2 &&
            prune_daughter_ldsk->disk->next) prune_daughter_ldsk = prune_daughter_ldsk->next[0];
      
      /* prune_daughter_ldsk has to be the next coalescent node under prune_ldsk for this move to work */
      
      /* Proba of pruning this particular one */
      hr -= LOG(1./prune_ldsk->n_next);
      
      /* /\* Proba of regrafting at that particular place in prune_ldsk->next *\/ */
      /* hr += LOG(1./prune_ldsk->n_next); */
      
      /* Get a ldsk to reattach the pruned lineage to */
      disk = tree->disk->prev;
      n_valid_disks = 0;
      do
        {
          /* Include only disks with displacement that are coalescent that are not younger
             than prune_daughter_ldsk */
          if((disk->ldsk != NULL) && 
             (disk->ldsk->n_next > 0) && 
             (disk->time < prune_daughter_ldsk->disk->time) && 
             (PHYREX_Is_On_Path(disk->ldsk,prune_daughter_ldsk,prune_ldsk) == NO))
            {
              if(!n_valid_disks) valid_disks = (t_dsk **)mCalloc(block,sizeof(t_dsk *));
              else if(!(n_valid_disks%block)) valid_disks = (t_dsk **)mRealloc(valid_disks,n_valid_disks+block,sizeof(t_dsk *));
              valid_disks[n_valid_disks] = disk;
              n_valid_disks++;
            }
          disk = disk->prev;
        }
      while(disk && disk->prev);
      
      if(!n_valid_disks) return;
      
      /* Uniform selection of a disk among the list of valid ones */
      i = Rand_Int(0,n_valid_disks-1);
      regraft_disk = valid_disks[i];
      Free(valid_disks);
      
      hr -= LOG(1./n_valid_disks);
      
      regraft_ldsk = regraft_disk->ldsk;
      
      /* /\* Proba of regafting at that particular place in regraft_ldsk->next *\/ */
      /* hr -= LOG(1./(phydbl)(regraft_ldsk->n_next+1)); */
      
      /* Proba of pruning that particular lineage */
      hr += LOG(1./(phydbl)(regraft_ldsk->n_next+1));
      
      
      /* /\* Prune and regraft *\/ */
      /* n_hits = PHYREX_Total_Number_Of_Hit_Disks(tree); */
      /* sizeT  = PHYREX_Time_Tree_Length(tree); */
      /* rate = (phydbl)n_hits/sizeT; */
      /* dt = FABS(prune_daughter_ldsk->disk->time - regraft_ldsk->disk->time); */
      /* new_path = PHYREX_Generate_Path(prune_daughter_ldsk,regraft_ldsk,rate*dt,tree); */
      
      /* cur_path_len = PHYREX_Path_Len(prune_daughter_ldsk,prune_ldsk)-2; */
      /* new_path_len = PHYREX_Path_Len(new_path,NULL)-1; */
      
      /* rate = (phydbl)(n_hits - cur_path_len + new_path_len)/sizeT; */
      /* dt = FABS(prune_daughter_ldsk->disk->time - prune_ldsk->disk->time); */
      /* hr += PHYREX_Path_Logdensity(prune_daughter_ldsk,prune_ldsk,rate*dt,tree); */
      
      /* cur_path = PHYREX_Remove_Path(prune_daughter_ldsk,prune_ldsk,tree); */
      /* PHYREX_Insert_Path(prune_daughter_ldsk,regraft_ldsk,new_path,tree); */
      
      /* rate = (phydbl)n_hits/sizeT; */
      /* dt = FABS(prune_daughter_ldsk->disk->time - regraft_ldsk->disk->time); */
      /* hr -= PHYREX_Path_Logdensity(prune_daughter_ldsk,regraft_ldsk,rate*dt,tree); */
      
      
      /* Prune and regraft */
      n_hits = PHYREX_Total_Number_Of_Hit_Disks(tree) - PHYREX_Total_Number_Of_Coal_Disks(tree);
      sizeT  = PHYREX_Time_Tree_Length(tree);
      dt = FABS(prune_daughter_ldsk->disk->time - prune_ldsk->disk->time);
      cur_path_len = PHYREX_Path_Len(prune_daughter_ldsk,prune_ldsk)-2;
      rate = (phydbl)(n_hits - cur_path_len)/(sizeT - dt);
      
      hr += PHYREX_Path_Logdensity(prune_daughter_ldsk,prune_ldsk,rate*dt,tree);
      
      dt = FABS(prune_daughter_ldsk->disk->time - regraft_ldsk->disk->time);
      new_path = PHYREX_Generate_Path(prune_daughter_ldsk,regraft_ldsk,rate*dt,tree);
      cur_path = PHYREX_Remove_Path(prune_daughter_ldsk,prune_ldsk,tree);
      PHYREX_Insert_Path(prune_daughter_ldsk,regraft_ldsk,new_path,tree);
      
      hr -= PHYREX_Path_Logdensity(prune_daughter_ldsk,regraft_ldsk,rate*dt,tree);
      
      
      
      
      /* prune_daughter_ldsk->prev = regraft_ldsk; */
      /* PHYREX_Remove_Lindisk_Next(prune_ldsk,prune_daughter_ldsk); */
      /* PHYREX_Random_Insert_Ldsk_In_Next_List(prune_daughter_ldsk,regraft_ldsk); */
      
      
      
      PHYREX_Ldsk_To_Tree(tree);
      Update_Ancestors(tree->n_root,tree->n_root->v[2],tree);
      Update_Ancestors(tree->n_root,tree->n_root->v[1],tree);
      RATES_Fill_Lca_Table(tree);
      RATES_Update_Cur_Bl(tree);
      new_glnL = PHYREX_Lk(tree);
      
      /* Get a ldsk from which you can prune a lineage */
      disk = tree->disk->prev;
      n_valid_disks = 0;
      do
        {
          /* Include only disks with displacement that are coalescent events */
          if(disk->ldsk != NULL && disk->ldsk->n_next > 1) 
            {
              /* Root has degree 2: don't pull any lineage */
              if(!(disk->prev == NULL && disk->ldsk->n_next == 2)) n_valid_disks++;
            }
          disk = disk->prev;
        }
      while(disk);
      hr += LOG(1./n_valid_disks);
      
      
      disk = tree->disk->prev;
      n_valid_disks = 0;
      do
        {
          if((disk->ldsk != NULL) && 
             (disk->ldsk->n_next > 0) && 
             (disk->time < prune_daughter_ldsk->disk->time) && 
             (PHYREX_Is_On_Path(disk->ldsk,prune_daughter_ldsk,regraft_ldsk) == NO))
            n_valid_disks++;
          disk = disk->prev;
        }
      while(disk && disk->prev);
      hr += LOG(1./n_valid_disks);
      
      if(tree->mcmc->use_data == YES) new_alnL = Lk(NULL,tree);
      
      ratio += (new_alnL - cur_alnL);
      ratio += (new_glnL - cur_glnL);
      ratio += hr;
      
      if(new_glnL < UNLIKELY + 0.1) 
        {
          PhyML_Printf("\n== %f %f",cur_glnL,new_glnL);
          PhyML_Printf("\n== %s",Write_Tree(tree,NO));
          Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
        }
      
      ratio = EXP(ratio);
      alpha = MIN(1.,ratio);
      
      /* Always accept move */
      if(tree->mcmc->always_yes == YES && new_glnL > UNLIKELY) alpha = 1.0;
      
      /* printf("\n. %3d %3d alpha:%12f %3d",prune_ldsk->n_next+1,regraft_ldsk->n_next-1,alpha,new_n_coal); */
      
      u = Uni();
      
      if(u > alpha) /* Reject */
        {
          new_path = PHYREX_Remove_Path(prune_daughter_ldsk,regraft_ldsk,tree);
          
          ldsk = new_path;
          while(ldsk)
            {
              Free_Disk(ldsk->disk);
              ldsk_dum = ldsk;
              ldsk = ldsk->prev;
              Free_Ldisk(ldsk_dum);
            }
          
          PHYREX_Insert_Path(prune_daughter_ldsk,prune_ldsk,cur_path,tree);
          
          /* prune_daughter_ldsk->prev = prune_ldsk; */
          /* PHYREX_Remove_Lindisk_Next(regraft_ldsk,prune_daughter_ldsk); */
          /* PHYREX_Insert_Ldsk_In_Next_List(prune_daughter_ldsk,prune_next_num,prune_ldsk); */
          
          PHYREX_Ldsk_To_Tree(tree);
          Update_Ancestors(tree->n_root,tree->n_root->v[2],tree);
          Update_Ancestors(tree->n_root,tree->n_root->v[1],tree);
          RATES_Fill_Lca_Table(tree);
          RATES_Update_Cur_Bl(tree);
          
          tree->c_lnL       = cur_alnL;
          tree->mmod->c_lnL = cur_glnL;
          
          /* if(tree->mcmc->use_data == YES) Lk(NULL,tree); */
          
          new_glnL = PHYREX_Lk(tree);
          if(Are_Equal(new_glnL,cur_glnL,1.E-3) == NO)
            {
              PhyML_Printf("\n. new_glnL: %f cur_glnL: %f",new_glnL,cur_glnL);
              Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
            }
        }
      else
        {
          ldsk = cur_path;
          while(ldsk)
            {
              Free_Disk(ldsk->disk);
              ldsk_dum = ldsk;
              ldsk = ldsk->prev;
              Free_Ldisk(ldsk_dum);
            }
          
          tree->mcmc->acc_move[tree->mcmc->num_move_phyrex_spr]++;
        }
      
      tree->mcmc->run_move[tree->mcmc->num_move_phyrex_spr]++;
      
    }
}
#endif
  
/*////////////////////////////////////////////////////////////l
////////////////////////////////////////////////////////////*/


/* #ifdef PHYREX */
/* void MCMC_PHYREX_Simulate_Backward(t_tree *tree) */
/* { */
/*   phydbl u,alpha,ratio,hr,T,t; */
/*   phydbl cur_glnL, new_glnL; */
/*   phydbl cur_alnL, new_alnL; */
/*   t_dsk  *disk,*bkp_disk,*target_disk; */
/*   t_ldsk **bkp_ldsk; */
/*   int i; */


/*   disk        = NULL; */
/*   bkp_disk    = NULL; */
/*   target_disk = NULL; */
/*   bkp_ldsk    = NULL; */

/*   new_glnL    = tree->mmod->c_lnL; */
/*   cur_glnL    = tree->mmod->c_lnL; */
/*   new_alnL    = tree->c_lnL; */
/*   cur_alnL    = tree->c_lnL; */
/*   hr          = 0.0; */
/*   ratio       = 0.0; */
/*   T           = 0.0; */
/*   t           = 0.0; */

/*   disk = tree->disk; */
/*   while(disk->prev) disk = disk->prev; */
/*   do */
/*     { */
/*       disk = disk->next; */
/*       if(disk->ldsk && disk->ldsk->n_next > 1) break; */
/*     } */
/*   while(1); */
/*   target_disk = disk; */

/*   /\* Chop off the tree at time t and simulate upwards from here *\/ */
/*   /\* T =  PHYREX_Tree_Height(tree); *\/ */
/*   /\* t = Uni()*T; *\/ */

/*   /\* disk = tree->disk; *\/ */
/*   /\* while(disk && disk->time > t) disk = disk->prev; *\/ */
/*   /\* target_disk = disk->next; *\/ */

/*   /\* hr -= LOG(FABS((target_disk->prev->time - target_disk->time)/T)); *\/ */

/*   bkp_disk = target_disk->prev; */

/*   bkp_ldsk = (t_ldsk **)mCalloc(target_disk->n_ldsk_a,sizeof(t_ldsk *)); */
/*   For(i,target_disk->n_ldsk_a) bkp_ldsk[i] = target_disk->ldsk_a[i]->prev; */
  
/*   PHYREX_Simulate_Backward_Core(NO,target_disk,tree); */

/*   /\* T =  PHYREX_Tree_Height(tree);   *\/ */
/*   /\* hr += LOG(FABS((target_disk->prev->time - target_disk->time)/T)); *\/ */
    
/*   if(tree->mcmc->use_data == YES) new_alnL = Lk(NULL,tree); */
/*   ratio += (new_alnL - cur_alnL); */
/*   ratio += hr; */
   
/*   ratio = EXP(ratio); */
/*   alpha = MIN(1.,ratio); */
    
/*  /\* Always accept move *\/ */
/*   if(tree->mcmc->always_yes == YES) alpha = 1.0; */
    
/*   u = Uni(); */
        
/*   if(u > alpha) /\* Reject *\/ */
/*     { */

/*       disk = target_disk->prev; */
/*       while(disk->prev) */
/*         { */
/*           disk = disk->prev; */
/*           if(disk->next->ldsk != NULL) Free_Ldisk(disk->next->ldsk); */
/*           Free_Disk(disk->next); */
/*         } */

/*       /\* Root *\/ */
/*       Free_Ldisk(disk->ldsk); */
/*       Free_Disk(disk); */

/*       target_disk->prev = bkp_disk; */
/*       For(i,target_disk->n_ldsk_a) target_disk->ldsk_a[i]->prev = bkp_ldsk[i]; */
      
/*       tree->c_lnL       = cur_alnL; */
/*       tree->mmod->c_lnL = cur_glnL; */

/*       /\* if(tree->mcmc->use_data == YES) Lk(NULL,tree); *\/ */

/*       new_glnL = PHYREX_Lk(tree); */

/*       if(Are_Equal(new_glnL,cur_glnL,1.E-3) == NO) */
/*         { */
/*           PhyML_Printf("\n. new_glnL: %f cur_glnL: %f",new_glnL,cur_glnL); */
/*           Generic_Exit(__FILE__,__LINE__,__FUNCTION__); */
/*         } */
/*     } */
/*   else */
/*     { */
/*       /\* Accept *\/ */
/*       disk = bkp_disk; */
/*       while(disk->prev) */
/*         { */
/*           disk = disk->prev; */
/*           if(disk->next->ldsk != NULL) Free_Ldisk(disk->next->ldsk); */
/*           Free_Disk(disk->next); */
/*         } */
      
/*       /\* Root *\/ */
/*       Free_Ldisk(disk->ldsk); */
/*       Free_Disk(disk); */

/*       /\* Likelihood needs to be updated *\/ */
/*       PHYREX_Lk(tree); */

/*       tree->mcmc->acc_move[tree->mcmc->num_move_phyrex_sim]++; */
/*     } */
  
/*   tree->mcmc->run_move[tree->mcmc->num_move_phyrex_sim]++; */

/*   Free(bkp_ldsk); */
/* } */
/* #endif */



#ifdef PHYREX
void MCMC_PHYREX_Simulate_Backward(t_tree *tree)
{
  phydbl u,alpha,ratio,hr,T,t;
  phydbl cur_glnL, new_glnL;
  phydbl cur_alnL, new_alnL;
  phydbl cur_lbda, new_lbda;
  phydbl cur_rad, new_rad;
  phydbl cur_mu, new_mu;
  t_dsk  *disk,*bkp_disk,*target_disk;
  t_ldsk **bkp_ldsk;
  int i;


  disk        = NULL;
  bkp_disk    = NULL;
  target_disk = NULL;
  bkp_ldsk    = NULL;

  new_glnL    = tree->mmod->c_lnL;
  cur_glnL    = tree->mmod->c_lnL;
  new_alnL    = tree->c_lnL;
  cur_alnL    = tree->c_lnL;
  hr          = 0.0;
  ratio       = 0.0;
  T           = 0.0;
  t           = 0.0;

  hr -= PHYREX_LnPrior_Radius(tree);
  hr -= PHYREX_LnPrior_Mu(tree);
  hr -= PHYREX_LnPrior_Lbda(tree);

  cur_lbda    = tree->mmod->lbda;
  cur_rad     = tree->mmod->rad;
  cur_mu      = tree->mmod->mu;

  new_lbda    = cur_lbda;
  new_rad     = cur_rad;
  new_mu      = cur_mu;


  /* u = Uni(); */
  /* if(u < 0.33)                  new_lbda    = cur_lbda + Rnorm(0.0,0.1); */
  /* else if(u < 0.66 && u > 0.33) new_mu      = cur_mu   + Rnorm(0.0,0.1); */
  /* else                          new_rad     = cur_rad  + Rnorm(0.0,0.1); */


  /* u = Uni(); */
  /* if(u < 0.33)                   */
  /*   { */
  /*     new_lbda = cur_lbda * EXP(2.0*(Uni()-.5)); */
  /*     hr += LOG(new_lbda/cur_lbda); */
  /*   } */
  /* else if(u < 0.66 && u > 0.33)  */
  /*   { */
  /*     new_mu = cur_mu * EXP(0.5*(Uni()-.5)); */
  /*     hr += LOG(new_mu/cur_mu); */
  /*   } */
  /* else                           */
  /*   { */
  /*     new_rad = cur_rad * EXP(2.0*(Uni()-.5)); */
  /*     hr += LOG(new_rad/cur_rad); */
  /*   } */


  /* new_lbda = cur_lbda * EXP(0.2*(Uni()-.5)); */
  /* hr += LOG(new_lbda/cur_lbda); */
  
  /* new_mu = cur_mu * EXP(0.05*(Uni()-.5)); */
  /* hr += LOG(new_mu/cur_mu); */

  /* new_rad = cur_rad * EXP(0.2*(Uni()-.5)); */
  /* hr += LOG(new_rad/cur_rad); */


  if(new_lbda > tree->mmod->max_lbda || new_lbda < tree->mmod->min_lbda) return;
  if(new_rad  > tree->mmod->max_rad  || new_rad  < tree->mmod->min_rad)  return;
  if(new_mu   > tree->mmod->max_mu   || new_mu   < tree->mmod->min_mu)   return;

  tree->mmod->lbda = new_lbda;
  tree->mmod->mu   = new_mu;
  tree->mmod->rad  = new_rad;

  hr += PHYREX_LnPrior_Radius(tree);
  hr += PHYREX_LnPrior_Mu(tree);
  hr += PHYREX_LnPrior_Lbda(tree);

  /* disk = tree->disk; */
  /* while(disk->prev) disk = disk->prev; */
  /* do */
  /*   { */
  /*     disk = disk->next; */
  /*     if(disk->ldsk && disk->ldsk->n_next > 1) break; */
  /*   } */
  /* while(1); */
  /* target_disk = disk; */

  /* Chop off the tree at time t and simulate upwards from here */
  T =  PHYREX_Tree_Height(tree);
  t = Uni()*T;

  disk = tree->disk;
  while(disk && disk->time > t) disk = disk->prev;
  target_disk = disk->next;

  hr -= LOG(FABS((target_disk->prev->time - target_disk->time)/T));

  bkp_disk = target_disk->prev;

  bkp_ldsk = (t_ldsk **)mCalloc(target_disk->n_ldsk_a,sizeof(t_ldsk *));
  For(i,target_disk->n_ldsk_a) bkp_ldsk[i] = target_disk->ldsk_a[i]->prev;
  
  PHYREX_Simulate_Backward_Core(NO,target_disk,tree);

  T =  PHYREX_Tree_Height(tree);
  hr += LOG(FABS((target_disk->prev->time - target_disk->time)/T));
    
  if(tree->mcmc->use_data == YES) new_alnL = Lk(NULL,tree);
  ratio += (new_alnL - cur_alnL);
  ratio += hr;
   
  ratio = EXP(ratio);
  alpha = MIN(1.,ratio);
    
 /* Always accept move */
  if(tree->mcmc->always_yes == YES) alpha = 1.0;
    
  u = Uni();
        
  if(u > alpha) /* Reject */
    {
      tree->mmod->lbda = cur_lbda;
      tree->mmod->mu   = cur_mu;
      tree->mmod->rad  = cur_rad;

      disk = target_disk->prev;
      while(disk->prev)
        {
          disk = disk->prev;
          if(disk->next->ldsk != NULL) Free_Ldisk(disk->next->ldsk);
          Free_Disk(disk->next);
        }

      /* Root */
      Free_Ldisk(disk->ldsk);
      Free_Disk(disk);

      target_disk->prev = bkp_disk;
      For(i,target_disk->n_ldsk_a) target_disk->ldsk_a[i]->prev = bkp_ldsk[i];
      
      tree->c_lnL       = cur_alnL;
      tree->mmod->c_lnL = cur_glnL;

      /* if(tree->mcmc->use_data == YES) Lk(NULL,tree); */

      new_glnL = PHYREX_Lk(tree);

      if(Are_Equal(new_glnL,cur_glnL,1.E-3) == NO)
        {
          PhyML_Printf("\n. new_glnL: %f cur_glnL: %f",new_glnL,cur_glnL);
          Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
        }
    }
  else
    {
      /* Accept */
      disk = bkp_disk;
      while(disk->prev)
        {
          disk = disk->prev;
          if(disk->next->ldsk != NULL) Free_Ldisk(disk->next->ldsk);
          Free_Disk(disk->next);
        }
      
      /* Root */
      Free_Ldisk(disk->ldsk);
      Free_Disk(disk);

      /* Likelihood needs to be updated */
      PHYREX_Lk(tree);

      tree->mcmc->acc_move[tree->mcmc->num_move_phyrex_sim]++;
    }
  
  tree->mcmc->run_move[tree->mcmc->num_move_phyrex_sim]++;

  Free(bkp_ldsk);
}
#endif


/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

#ifdef PHYREX
void MCMC_PHYREX_Lineage_Traj(t_tree *tree)
{
  phydbl u,alpha,ratio,hr;
  phydbl cur_glnL, new_glnL;
  t_dsk  *disk,**valid_disks;
  t_ldsk *start_ldsk,*end_ldsk,*cur_path,*new_path,*ldsk,*ldsk_dum;
  int i,block,n_valid_disks;
  phydbl rate,dt,sizeT;
  int cur_path_len;
  int n_hits,n_iter;

  /* n_hits = PHYREX_Total_Number_Of_Hit_Disks(tree); */
  /* n_iter = MAX(tree->n_otu,(int)(n_hits/5)); */
  n_iter = tree->n_otu;

  while(n_iter--)
    {
      
      valid_disks = NULL;
      disk        = NULL;
      new_glnL    = tree->mmod->c_lnL;
      cur_glnL    = tree->mmod->c_lnL;
      hr          = 0.0;
      ratio       = 0.0;
      block       = 100;
      
      if(tree->disk->next) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
      
      /* Get a ldsk from which you can prune a lineage */
      disk = tree->disk->prev;
      n_valid_disks = 0;
      do
        {
          if(disk->ldsk && disk->ldsk->n_next > 1)
            /* if(disk->ldsk) */
            {
              if(!n_valid_disks) valid_disks = (t_dsk **)mCalloc(block,sizeof(t_dsk *));
              else if(!(n_valid_disks%block)) valid_disks = (t_dsk **)mRealloc(valid_disks,n_valid_disks+block,sizeof(t_dsk *));
              valid_disks[n_valid_disks] = disk;
              n_valid_disks++;
            }
          disk = disk->prev;
        }
      while(disk);
      
      if(!n_valid_disks) return;
      
      /* Uniform selection of a target disk */
      i = Rand_Int(0,n_valid_disks-1);
      end_ldsk = valid_disks[i]->ldsk;
      Free(valid_disks);
      
      i = Rand_Int(0,end_ldsk->n_next-1);
      start_ldsk = end_ldsk->next[i];
      while(start_ldsk->n_next < 2 && start_ldsk->disk->next) start_ldsk = start_ldsk->next[0];
      
      if(end_ldsk == NULL)   Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
      if(start_ldsk == NULL) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
      
      /* n_hits = PHYREX_Total_Number_Of_Hit_Disks(tree) - PHYREX_Total_Number_Of_Coal_Disks(tree); */
      /* sizeT  = PHYREX_Time_Tree_Length(tree); */
      /* dt     = FABS(start_ldsk->disk->time - end_ldsk->disk->time); */
      
      /* rate = (phydbl)n_hits/sizeT;  */
      /* new_path = PHYREX_Generate_Path(start_ldsk,end_ldsk,rate*dt,tree); */
      
      /* cur_path_len = PHYREX_Path_Len(start_ldsk,end_ldsk)-2; */
      /* new_path_len = PHYREX_Path_Len(new_path,NULL)-1; */
      
      /* rate = (phydbl)(n_hits - cur_path_len + new_path_len)/sizeT; */
      /* hr += PHYREX_Path_Logdensity(start_ldsk,end_ldsk,rate*dt,tree); */
      
      /* cur_path = PHYREX_Remove_Path(start_ldsk,end_ldsk,tree); */
      /* PHYREX_Insert_Path(start_ldsk,end_ldsk,new_path,tree); */
      
      /* rate = (phydbl)n_hits/sizeT; */
      /* hr -= PHYREX_Path_Logdensity(start_ldsk,end_ldsk,rate*dt,tree); */
      
      
      n_hits       = PHYREX_Total_Number_Of_Hit_Disks(tree) - PHYREX_Total_Number_Of_Coal_Disks(tree);
      sizeT        = PHYREX_Time_Tree_Length(tree);
      dt           = FABS(start_ldsk->disk->time - end_ldsk->disk->time);
      cur_path_len = PHYREX_Path_Len(start_ldsk,end_ldsk)-2;
      rate         = (phydbl)(n_hits - cur_path_len)/(sizeT - dt);
      
      hr += PHYREX_Path_Logdensity(start_ldsk,end_ldsk,rate*dt,tree);
      
      new_path = PHYREX_Generate_Path(start_ldsk,end_ldsk,rate*dt,tree);
      cur_path = PHYREX_Remove_Path(start_ldsk,end_ldsk,tree);
      PHYREX_Insert_Path(start_ldsk,end_ldsk,new_path,tree);
      
      hr -= PHYREX_Path_Logdensity(start_ldsk,end_ldsk,rate*dt,tree);
      
      new_glnL = PHYREX_Lk(tree);
      
      ratio += (new_glnL - cur_glnL);
      ratio += hr;
      
      if(new_glnL < UNLIKELY + 0.1) 
        {
          PhyML_Printf("\n== %f %f",cur_glnL,new_glnL);
          PhyML_Printf("\n== %s",Write_Tree(tree,NO));
          Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
        }
      
      ratio = EXP(ratio);
      alpha = MIN(1.,ratio);
      
      /* PhyML_Printf("\n= Traj. new_lnL:%12f cur_lnL:%12f hr:%12f alpha:%12f",new_glnL,cur_glnL,hr,alpha); */
      
      if(tree->mcmc->always_yes == YES && new_glnL > UNLIKELY) alpha = 1.0;
      
      u = Uni();
      
      if(u > alpha) /* Reject */
        {
          new_path = PHYREX_Remove_Path(start_ldsk,end_ldsk,tree);
          
          ldsk = new_path;
          while(ldsk) 
            {
              Free_Disk(ldsk->disk);
              ldsk_dum = ldsk;
              ldsk = ldsk->prev;
              Free_Ldisk(ldsk_dum);
            }
          
          PHYREX_Insert_Path(start_ldsk,end_ldsk,cur_path,tree);
          
          tree->mmod->c_lnL = cur_glnL;
          
          new_glnL = PHYREX_Lk(tree);
          
          if(Are_Equal(new_glnL,cur_glnL,1.E-3) == NO)
            {
              PhyML_Printf("\n. new_glnL: %f cur_glnL: %f",new_glnL,cur_glnL);
              Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
            }
        }
      else
        {
          ldsk = cur_path;
          while(ldsk) 
            {
              Free_Disk(ldsk->disk);
              ldsk_dum = ldsk;
              ldsk = ldsk->prev;
              Free_Ldisk(ldsk_dum);
            }      
          
          tree->mcmc->acc_move[tree->mcmc->num_move_phyrex_traj]++;
        }
      
      tree->mcmc->run_move[tree->mcmc->num_move_phyrex_traj]++;
    }
}
#endif
  
/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

