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
  Lk(NULL,tree);
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
      	  Lk(NULL,tree);
          Set_Both_Sides(NO,tree);     

      	  if(tree->mcmc->is == NO || tree->rates->model_log_rates == YES)
      	    {
              MCMC_Root_Time(tree);
	      MCMC_Time_Recur(tree->n_root,tree->n_root->v[first],YES,tree);
	      MCMC_Time_Recur(tree->n_root,tree->n_root->v[secod],YES,tree);
      	    }
      	  else
      	    {
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
      	  Lk(NULL,tree);
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
      	  /* if(tree->eval_alnL == YES) Lk(NULL,tree); */
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
	      Sample_Ancestral_Seq(YES,!tree->eval_alnL,tree);
	      
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
      
      /* if(like_func && tree->eval_alnL == YES)  /\* Likelihood ratio *\/ */
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
  t_edge *b2,*b3;
  t_node *v2,*v3;
  int move_num;
  phydbl K;

  if(tree->rates->model == STRICTCLOCK) return;

  b = NULL;
  if(a == tree->n_root) b = tree->e_root;
  else For(i,3) if(d->v[i] == a) { b = d->b[i]; break; }
   
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

  
  u = Uni();
  
  MCMC_Make_Move(&cur_mu,&new_mu,r_min,r_max,&ratio,K,tree->mcmc->move_type[tree->mcmc->num_move_br_r+d->num]);

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
      
      if(tree->eval_alnL == YES && tree->io->lk_approx == EXACT)
        {
          Update_PMat_At_Given_Edge(b,tree);
          if(!d->tax)
            {
              Update_PMat_At_Given_Edge(b2,tree);
              Update_PMat_At_Given_Edge(b3,tree);
              Update_P_Lk(tree,b,d);
            }
        }
      new_lnL_data = Lk(b,tree);
      
      /* tree->both_sides = NO; */
      /* new_lnL_data = Lk(tree); */
    
      
      new_lnL_rate = RATES_Lk_Rates(tree);
      
      ratio += (new_lnL_data - cur_lnL_data);
      ratio += (new_lnL_rate - cur_lnL_rate);
            
      ratio = EXP(ratio);
      alpha = MIN(1.,ratio);
      
      u = Uni();
      
      if(u > alpha) /* Reject */
	{
	  tree->rates->br_r[d->num] = cur_mu;	  
	  tree->c_lnL               = cur_lnL_data;
	  tree->rates->c_lnL_rates  = cur_lnL_rate;
	  
          RATES_Update_Cur_Bl(tree);
	  
	  if(tree->io->lk_approx == EXACT)
	    {
	      Update_PMat_At_Given_Edge(b,tree);
	      if(!d->tax)
		{
		  Update_PMat_At_Given_Edge(b2,tree);
		  Update_PMat_At_Given_Edge(b3,tree);
                  Update_P_Lk(tree,b,d);
		}
	    }
	  
          /* /\* !!!!!!!!!!!!!! *\/ */
          /* printf("\n. Rej"); */
          /* new_lnL_data = Lk(NULL,tree); /\* Not necessary. Remove once tested *\/ */
          /* if(Are_Equal(new_lnL_data,cur_lnL_data,1.E-3) == NO) */
          /*   { */
          /*     PhyML_Printf("\n== a: %d d: %d v2: %d v3: %d",a->num,d->num,v2?v2->num:-1,v3?v3->num:-1); */
          /*     PhyML_Printf("\n== new_alnL: %f cur_alnL: %f",new_lnL_data,cur_lnL_data); */
          /*     Generic_Exit(__FILE__,__LINE__,__FUNCTION__); */
          /*   } */
        }
      else
	{
          tree->mcmc->acc_move[tree->mcmc->num_move_br_r+d->num]++;
          /* /\* !!!!!!!!!!!!!! *\/ */
          /* printf("\n. Acc"); */
          /* cur_lnL_data = Lk(NULL,tree); /\* Not necessary. Remove once tested *\/ */
          /* if(Are_Equal(new_lnL_data,cur_lnL_data,1.E-3) == NO) */
          /*   { */
          /*     PhyML_Printf("\n== a: %d d: %d v2: %d v3: %d",a->num,d->num,v2?v2->num:-1,v3?v3->num:-1); */
          /*     PhyML_Printf("\n== new_alnL: %f cur_alnL: %f",new_lnL_data,cur_lnL_data); */
          /*     Generic_Exit(__FILE__,__LINE__,__FUNCTION__); */
          /*   } */
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
                Update_P_Lk(tree,d->b[i],d);
                MCMC_One_Rate(d,d->v[i],YES,tree);
              }
          Update_P_Lk(tree,b,d);
        }
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

void MCMC_Time_All(t_tree *tree)
{
  // Down partials may not be up to date.
  Set_Both_Sides(YES,tree);
  Lk(NULL,tree);
  Set_Both_Sides(NO,tree);
  MCMC_Root_Time(tree);
  MCMC_Time_Recur(tree->n_root,tree->n_root->v[1],YES,tree);
  MCMC_Time_Recur(tree->n_root,tree->n_root->v[2],YES,tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MCMC_Time_Recur(t_node *a, t_node *d, int traversal, t_tree *tree)
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
  int move_num;
  
  if(d->tax) return; /* Won't change time at tip */
  
  move_num     = tree->mcmc->num_move_times;
  cur_lnL_data = tree->c_lnL;
  cur_lnL_rate = tree->rates->c_lnL_rates;
  t1_cur       = tree->rates->nd_t[d->num];
  new_lnL_data = cur_lnL_data;
  new_lnL_rate = cur_lnL_rate;
  ratio        = 0.0;
  cur_lnL_time = tree->rates->c_lnL_times;
  new_lnL_time = cur_lnL_time;

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

  t_min = MAX(t0,tree->rates->t_prior_min[d->num]);
  t_max = MIN(tree->rates->t_prior_max[d->num],MIN(t2,t3));

  t_min += tree->rates->min_dt;
  t_max -= tree->rates->min_dt;
       
  /* MCMC_Make_Move(&t1_cur,&t1_new,t_min,t_max,&ratio,K,tree->mcmc->move_type[move_num]); */
  t1_new = Uni()*(t_max - t_min) + t_min;

  if(t1_new > t_min && t1_new < t_max) 
    {
      RATES_Record_Times(tree);

      tree->rates->nd_t[d->num] = t1_new;

      new_lnL_time = TIMES_Lk_Times(tree);

      ratio += (new_lnL_time - cur_lnL_time);

      if(new_lnL_time > UNLIKELY)
        {
          RATES_Update_Cur_Bl(tree);
          
          new_lnL_rate = RATES_Lk_Rates(tree);
          ratio += (new_lnL_rate - cur_lnL_rate);
          
          if(tree->eval_alnL == YES && tree->io->lk_approx == EXACT)
            {
              Update_PMat_At_Given_Edge(b1,tree);
              Update_PMat_At_Given_Edge(b2,tree);
              Update_PMat_At_Given_Edge(b3,tree);
              Update_P_Lk(tree,b1,d);
            }
          new_lnL_data = Lk(b1,tree);
          
          ratio += (new_lnL_data - cur_lnL_data);
        }
      
      /* printf("\n. One_Time cur_t: %f new_t: %f ratio: %f alnL:%f->%f tlnL: %f->%f rlnL: %f->%f", */
      /*        t1_cur,t1_new, */
      /*        ratio, */
      /*        cur_lnL_data,new_lnL_data, */
      /*        cur_lnL_time,new_lnL_time, */
      /*        cur_lnL_rate,new_lnL_rate); */

      ratio = EXP(ratio);
      alpha = MIN(1.,ratio);
      u = Uni();
	           
      if(u > alpha) /* Reject */
	{
          /* printf("\n. rej"); */
          RATES_Reset_Times(tree);
          RATES_Update_Cur_Bl(tree);
              
          if(tree->eval_alnL == YES && tree->io->lk_approx == EXACT) 
            {
              Update_PMat_At_Given_Edge(b1,tree);
              Update_PMat_At_Given_Edge(b2,tree);
              Update_PMat_At_Given_Edge(b3,tree);
              Update_P_Lk(tree,b1,d);
            }

          if(isinf(FABS(new_lnL_time)) == YES || isnan(new_lnL_time) == YES)
            {
              Print_Node(tree->n_root,tree->n_root->v[1],tree);
              Print_Node(tree->n_root,tree->n_root->v[2],tree);              
              assert(FALSE);
            }

	  tree->c_lnL               = cur_lnL_data;
	  tree->rates->c_lnL_rates  = cur_lnL_rate;
          tree->rates->c_lnL_times  = cur_lnL_time;
          DATE_Update_T_Prior_MinMax(tree);


          if(Are_Equal(tree->rates->c_lnL_times,cur_lnL_time,1.E-3) == NO)
            {
              PhyML_Printf("\n\n");
              PhyML_Printf("\n. moved node %d from %f to %f\n",d->num,t1_cur,t1_new);
              Print_Node(tree->n_root,tree->n_root->v[1],tree);
              Print_Node(tree->n_root,tree->n_root->v[2],tree);              
              PhyML_Printf("\n== new_glnL: %f cur_glnL: %f",tree->rates->c_lnL_times,cur_lnL_time);
              Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
            }
          
          /* !!!!!!!!!!!!!! */
          /* new_lnL_data = Lk(NULL,tree); /\* Not necessary. Remove once tested *\/ */
          /* if(Are_Equal(new_lnL_data,cur_lnL_data,1.E-3) == NO) */
          /*   { */
          /*     PhyML_Printf("\n== a: %d d: %d v2: %d v3: %d",a->num,d->num,v2->num,v3->num); */
          /*     PhyML_Printf("\n== t1_cur: %f t1_new: %f",t1_cur,t1_new); */
          /*     PhyML_Printf("\n== new_alnL: %f cur_alnL: %f",new_lnL_data,cur_lnL_data); */
          /*     Generic_Exit(__FILE__,__LINE__,__FUNCTION__); */
          /*   } */
	}
      else
	{
          /* !!!!!!!!!!!!!! */
          /* cur_lnL_data = Lk(NULL,tree); /\* Not necessary. Remove once tested *\/ */
          /* if(Are_Equal(new_lnL_data,cur_lnL_data,1.E-3) == NO) */
          /*   { */
          /*     PhyML_Printf("\n== a: %d d: %d v2: %d v3: %d",a->num,d->num,v2->num,v3->num); */
          /*     PhyML_Printf("\n== t1_cur: %f t1_new: %f",t1_cur,t1_new); */
          /*     PhyML_Printf("\n== new_alnL: %f cur_alnL: %f",new_lnL_data,cur_lnL_data); */
          /*     Generic_Exit(__FILE__,__LINE__,__FUNCTION__); */
          /*   } */
	  tree->mcmc->acc_move[move_num]++;
	}

      /* printf("\n. %f",new_lnL_data); fflush(NULL); */

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
	  PhyML_Printf("\n== Err. in file %s at line %d\n",__FILE__,__LINE__);
	}
    }
  
  tree->mcmc->run_move[move_num]++;

  if(traversal == YES)
    {
      if(d->tax == YES) return;
      else
	{
          For(i,3)
            if(d->v[i] != a && d->b[i] != tree->e_root)
              {
                Update_P_Lk(tree,d->b[i],d);
                MCMC_Time_Recur(d,d->v[i],YES,tree);
              }
          Update_P_Lk(tree,b1,d);
        }
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
          new_lnL_data = Lk(NULL,tree);
          
          new_lnL_rate = RATES_Lk_Rates(tree);
          new_lnL_time = TIMES_Lk_Times(tree);
          
          /* printf("\n. JUMP cur_lnL_time: %f new_lnL_time: %f",cur_lnL_time,new_lnL_time); */
          /* for(i = tree -> n_otu; i < 2 * tree -> n_otu -1; i++) printf("\nJUMP1 Node number:[%d] Lower bound:[%f] Upper bound:[%f] Node time:[%f].", i, */
          /*                                                          tree -> rates -> t_prior_min[i], */
          /*                                                          tree -> rates -> t_prior_max[i], */
          /*                                                          tree -> rates -> nd_t[i]); */

          /* Likelihood ratio */
          ratio += (new_lnL_data - cur_lnL_data);
          
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

  /* printf("\n IN ROOT: %f",tree->rates->c_lnL_times); */

  root = tree->n_root;

  if(FABS(tree->rates->t_prior_min[root->num] - tree->rates->t_prior_max[root->num]) < 1.E-10) return;
  
  move_num       = tree->mcmc->num_move_root_time;
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

  t_min = MAX(t0,2.*MIN(t2,t3));
  t_max = MIN(tree->rates->t_prior_max[root->num],MIN(t2,t3));

  t_min += tree->rates->min_dt;
  t_max -= tree->rates->min_dt;

  if(t_min > t_max) 
    {
      PhyML_Printf("\n== glnL:%f",TIMES_Lk_Times(tree));
      PhyML_Printf("\n== t:%f",tree->rates->nd_t[tree->n_root->num]);
      PhyML_Printf("\n== t0 = %f t2 = %f t3 = %f",t0,t2,t3);
      PhyML_Printf("\n== t_min = %f t_max = %f",t_min,t_max);
      PhyML_Printf("\n== prior_min = %f prior_max = %f",tree->rates->t_prior_min[root->num],tree->rates->t_prior_max[root->num]);
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
    }


  MCMC_Make_Move(&t1_cur,&t1_new,t_min,t_max,&ratio,K,tree->mcmc->move_type[move_num]);
  /* t1_new = Uni()*(t_max-t_min) + t_min; */

  /* PhyML_Printf("\n. K: %f t1: %f->%f acc.num: %d run: %d acc.rate: %f", */
  /*              K,t1_cur,t1_new, */
  /*              tree->mcmc->acc_move[move_num], */
  /*              tree->mcmc->run_move[move_num], */
  /*              tree->mcmc->run_move[move_num] > 0 ? (phydbl)tree->mcmc->acc_move[move_num]/tree->mcmc->run_move[move_num] : -1.); */

  if(t1_new > t_min && t1_new < t_max) 
    {
      RATES_Record_Times(tree);

      tree->rates->nd_t[root->num] = t1_new;

      RATES_Update_Cur_Bl(tree);

      /* !!!!!! What if model==GUINDON? */
      new_lnL_data = Lk(b1,tree);

      new_lnL_rate = RATES_Lk_Rates(tree);
      new_lnL_time = TIMES_Lk_Times(tree); 

      ratio += (new_lnL_data - cur_lnL_data);
      ratio += (new_lnL_rate - cur_lnL_rate);
      ratio += (new_lnL_time - cur_lnL_time);

      ratio = EXP(ratio);
      alpha = MIN(1.,ratio);
      u = Uni();
         
      /* printf("\n. alnL:%f->%f tlnL: %f->%f rlnL: %f->%f ratio: %f", */
      /*        cur_lnL_data,new_lnL_data, */
      /*        cur_lnL_time,new_lnL_time, */
      /*        cur_lnL_rate,new_lnL_rate, */
      /*        ratio); */

      if(u > alpha) /* Reject */
	{
	  RATES_Reset_Times(tree);
          RATES_Update_Cur_Bl(tree);
          if(tree->eval_alnL == YES) Update_PMat_At_Given_Edge(b1,tree);

	  tree->c_lnL              = cur_lnL_data;
	  tree->rates->c_lnL_rates = cur_lnL_rate;
	  tree->rates->c_lnL_times = TIMES_Lk_Times(tree); // required as some t_prior_min/max have been modified 
          
          if(Are_Equal(tree->rates->c_lnL_times,cur_lnL_time,1.E-3) == NO)
            {
              PhyML_Printf("\n\n");
              Print_Node(tree->n_root,tree->n_root->v[1],tree);
              Print_Node(tree->n_root,tree->n_root->v[2],tree);
              PhyML_Printf("\n== new_glnL: %f cur_glnL: %f",tree->rates->c_lnL_times,cur_lnL_time);
              Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
            }
	}
      else
	{
	  tree->mcmc->acc_move[move_num]+=1;
	}
      
      if(t1_new < t0)
	{
	  t1_new = t0+1.E-4;
	  PhyML_Printf("\n");
	  PhyML_Printf("\n== t0 = %f t1_new = %f",t0,t1_new);
	  PhyML_Printf("\n== t_min=%f t_max=%f",t_min,t_max);
	  PhyML_Printf("\n== (t1-t0)=%f (t2-t1)=%f",t1_cur-t0,t2-t1_cur);
	  PhyML_Printf("\n== Err. in file %s at line %d\n",__FILE__,__LINE__);
          /* Exit("\n"); */
	}
      if(t1_new > MIN(t2,t3))
	{
	  PhyML_Printf("\n");
	  PhyML_Printf("\n== t0 = %f t1_new = %f t1 = %f t2 = %f t3 = %f MIN(t2,t3)=%f",t0,t1_new,t1_cur,t2,t3,MIN(t2,t3));
	  PhyML_Printf("\n== t_min=%f t_max=%f",t_min,t_max);
	  PhyML_Printf("\n== (t1-t0)=%f (t2-t1)=%f",t1_cur-t0,t2-t1_cur);
	  PhyML_Printf("\n== Err. in file %s at line %d\n",__FILE__,__LINE__);
          /* Exit("\n"); */
	}
      
      if(isnan(t1_new))
	{
	  PhyML_Printf("\n== run=%d",tree->mcmc->run);
	  PhyML_Printf("\n== Err. in file %s at line %d\n",__FILE__,__LINE__);
          /* Exit("\n"); */
	}
    }
    
  tree->mcmc->run_move[move_num]+=1;

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

  K            = tree->mcmc->tune_move[tree->mcmc->num_move_tree_height];
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
  
  For(i,2*tree->n_otu-2) tree->rates->br_do_updt[i] = YES;
  RATES_Update_Cur_Bl(tree);
  new_lnL_data = Lk(NULL,tree);

  new_lnL_rate = RATES_Lk_Rates(tree);
  new_lnL_time = TIMES_Lk_Times(tree);
 
  ratio += (phydbl)(n_nodes)*LOG(mult);

  /* Likelihood ratio */
  ratio += (new_lnL_data - cur_lnL_data);

  /* Prior ratio */
  ratio += (new_lnL_rate - cur_lnL_rate);
  ratio += (new_lnL_time - cur_lnL_time);

  /* printf("\n. data: %f %f times: %f %f rates: %f %f", */
  /*        new_lnL_data,cur_lnL_data, */
  /*        new_lnL_time,cur_lnL_time, */
  /*        new_lnL_rate,cur_lnL_rate); */


  ratio = EXP(ratio);
  alpha = MIN(1.,ratio);
  u = Uni();
  
  if(u > alpha)
    {
      RATES_Reset_Times(tree);
      tree->rates->c_lnL_times = TIMES_Lk_Times(tree); // Required in order to set t_prior_min/max to their original values
      tree->c_lnL              = cur_lnL_data;
      tree->rates->c_lnL_rates = cur_lnL_rate;

      if(Are_Equal(tree->rates->c_lnL_times,cur_lnL_time,1.E-3) == NO)
        {
          PhyML_Printf("\n== new_glnL: %f cur_glnL: %f",tree->rates->c_lnL_times,cur_lnL_time);
          Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
        }
    }
  else
    {
      tree->mcmc->acc_move[tree->mcmc->num_move_tree_height]++;
    }

  tree->mcmc->run_move[tree->mcmc->num_move_tree_height]++;
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

  K = 0.1;

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
  new_lnL_data = Lk(NULL,tree);

  new_lnL_rate = RATES_Lk_Rates(tree);
  new_lnL_time = TIMES_Lk_Times(tree);

  /* The Hastings ratio is actually mult^(n) when changing the absolute
     node heights. When considering the relative heights, this ratio combined
     to the Jacobian for the change of variable ends up to being equal to mult. 
  */
  ratio += (n_nodes - 1)*LOG(mult);
  /* ratio += -LOG(mult) + LOG(Dgamma(1./mult,1./K,K)/Dgamma(mult,1./K,K)); */

  /* Likelihood ratio */
  ratio += (new_lnL_data - cur_lnL_data);

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
  
  For(i,2*tree->n_otu-2) tree->rates->br_do_updt[i] = YES;
  RATES_Update_Cur_Bl(tree);
  new_lnL_data = Lk(NULL,tree);

  new_lnL_rate = RATES_Lk_Rates(tree);
  new_lnL_time = TIMES_Lk_Times(tree);

  /* The Hastings ratio is actually mult^(n) when changing the absolute
     node heights. When considering the relative heights, this ratio combined
     to the Jacobian for the change of variable ends up to being equal to mult. 
  */
  ratio += (n_nodes - 1)*LOG(mult);
  /* ratio += -LOG(mult) + LOG(Dgamma(1./mult,1./K,K)/Dgamma(mult,1./K,K)); */

  /* Likelihood ratio */
  ratio += (new_lnL_data - cur_lnL_data);

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
  new_lnL_data = Lk(NULL,tree);

  new_lnL_rate = RATES_Lk_Rates(tree);
  new_lnL_time = TIMES_Lk_Times(tree);

  /* The Hastings ratio here is mult^(n_nodes) and the ratio of the prior joint densities
     of the modified node heigths given the unchanged one is 1. This is different from the 
     case where all the nodes, including the root node, are scaled. 
  */
  ratio += (phydbl)(n_nodes)*LOG(mult);

  /* Likelihood ratio */
  ratio += (new_lnL_data - cur_lnL_data);

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
      new_lnL_data = Lk(NULL,tree);
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
  new_lnL_data = Lk(NULL,tree);

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
  new_lnL_data = Lk(NULL,tree);

  new_lnL_rate = RATES_Lk_Rates(tree);

  ratio += (-(tree->n_otu-1.)-2.)*LOG(mult);
  ratio += (new_lnL_rate - cur_lnL_rate);
  ratio += (new_lnL_data - cur_lnL_data);

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
  new_lnL_data = Lk(NULL,tree);

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

void MCMC_Time_Slice(t_tree *tree)
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
  
  RATES_Record_Times(tree);
  
  K            = tree->mcmc->tune_move[tree->mcmc->num_move_time_slice];
  cur_lnL_data = tree->c_lnL;
  new_lnL_data = tree->c_lnL;
  ratio        = 0.0;
  cur_lnL_rate = tree->rates->c_lnL_rates;
  new_lnL_rate = tree->rates->c_lnL_rates;
  cur_lnL_time = tree->rates->c_lnL_times;
  
  u = Uni();
  mult = EXP(K*(u-0.5));

  floor = Uni()*tree->rates->nd_t[tree->n_root->num];

  Scale_Subtree_Height(tree->n_root,mult,floor,&n_nodes,tree);

  if(TIMES_Check_Node_Height_Ordering(tree) != NO)
    {
      For(i,2*tree->n_otu-2) tree->rates->br_do_updt[i] = YES;
      RATES_Update_Cur_Bl(tree);
      new_lnL_data = Lk(NULL,tree);
      
      new_lnL_rate = RATES_Lk_Rates(tree);
      new_lnL_time = TIMES_Lk_Times(tree);
      
      ratio += (phydbl)(n_nodes)*LOG(mult);

      /* Likelihood ratio */
      ratio += (new_lnL_data - cur_lnL_data);
      
      /* Prior ratio */
      ratio += (new_lnL_rate - cur_lnL_rate);
      ratio += (new_lnL_time - cur_lnL_time);
      
      /* printf("\n. data: %f<-%f times: %f<-%f rates: %f<-%f acc.: [%d/%d] k: %f t: %f x: %f", */
      /*        new_lnL_data,cur_lnL_data, */
      /*        new_lnL_time,cur_lnL_time, */
      /*        new_lnL_rate,cur_lnL_rate, */
      /*        tree->mcmc->acc_move[tree->mcmc->num_move_time_slice], */
      /*        tree->mcmc->run_move[tree->mcmc->num_move_time_slice], */
      /*        K,floor,mult); */
      
      ratio = EXP(ratio);
    }
  else
    {
      ratio = 0.0;
    }

  alpha = MIN(1.,ratio);
  u = Uni();
  
  if(u > alpha)
    {
      RATES_Reset_Times(tree);
      tree->rates->c_lnL_times = TIMES_Lk_Times(tree); // Required in order to set t_prior_min/max to their original values
      tree->c_lnL              = cur_lnL_data;
      tree->rates->c_lnL_rates = cur_lnL_rate;

      if(Are_Equal(tree->rates->c_lnL_times,cur_lnL_time,1.E-3) == NO)
        {
          PhyML_Printf("\n== new_glnL: %f cur_glnL: %f",tree->rates->c_lnL_times,cur_lnL_time);
          Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
        }
    }
  else
    {
      tree->mcmc->acc_move[tree->mcmc->num_move_time_slice]++;
    }

  tree->mcmc->run_move[tree->mcmc->num_move_time_slice]++;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MCMC_Print_Param_Stdin(t_mcmc *mcmc, t_tree *tree)
{
  time_t cur_time;

  time(&cur_time);
  


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
      PhyML_Printf("  %9f[%5.0f]",
                   (tree->rates ? tree->rates->nu : -1.),
                   tree->mcmc->ess[tree->mcmc->num_move_nu]);
      PhyML_Printf("  %8f[%5.0f]",
                   (tree->rates ? tree->rates->birth_rate : -1.),
                   tree->mcmc->ess[tree->mcmc->num_move_birth_rate]);
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
      Lk(NULL,tree);
      PhyML_Fprintf(fp,"%.1f\t",tree->c_lnL);
      tree->io->lk_approx = NORMAL;
      tree->c_lnL = 0.0;
      Lk(NULL,tree);
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

void MCMC_Randomize_Death(t_tree *tree)
{
  phydbl min_d,max_d;
  phydbl u;

  min_d = tree->rates->death_rate_min;
  max_d = MIN(MIN(0.5,tree->rates->death_rate_max),tree->rates->birth_rate);
  
  u = Uni();
  tree->rates->death_rate = (max_d - min_d) * u + min_d;
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

  printf("\n. ROOT: %f %f %d",t_inf,t_sup,tree->n_root->num);
  printf("\n. ROOT: %f",u);

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
	  rate_inf = 0.234;
	  rate_sup = 0.234;
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
	  rate_inf = 0.234;
	  rate_sup = 0.234;
	}

      /* if(!strcmp(mcmc->move_name[move],"tree_rates")) */
      /* 	{ */
      /* 	  rate_inf = 0.05; */
      /* 	  rate_sup = 0.05; */
      /* 	} */
      else if(!strcmp(mcmc->move_name[move],"phyrex_lbda"))
	{
	  rate_inf = 0.234;
	  rate_sup = 0.234;
	}
      else if(!strcmp(mcmc->move_name[move],"phyrex_mu"))
	{
	  rate_inf = 0.234;
	  rate_sup = 0.234;
	}
      else if(!strcmp(mcmc->move_name[move],"phyrex_rad"))
	{
	  rate_inf = 0.234;
	  rate_sup = 0.234;
	}
      else if(!strcmp(mcmc->move_name[move],"phyrex_ldsk_and_disk"))
	{
	  rate_inf = 0.234;
	  rate_sup = 0.234;
	}
      else if(!strcmp(mcmc->move_name[move],"phyrex_ldsk_multi"))
	{
	  rate_inf = 0.234;
	  rate_sup = 0.234;
	}
      else if(!strcmp(mcmc->move_name[move],"phyrex_disk_multi"))
	{
	  rate_inf = 0.234;
	  rate_sup = 0.234;
	}
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
      Lk(NULL,tree);
      change = YES;
    }
  
  For(i,2*tree->n_otu-2) tree->rates->br_do_updt[i] = NO;
  MCMC_Single_Param_Generic(&(tree->mod->kappa->v),0.,100.,tree->mcmc->num_move_kappa,
			    NULL,&(tree->c_lnL),
			    NULL,Wrap_Lk,tree->mcmc->move_type[tree->mcmc->num_move_kappa],NO,NULL,tree,NULL);
	  
  if(change == YES)
    {
      tree->io->lk_approx = NORMAL;
      Lk(NULL,tree);
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
      ratio += (new_lnL_data - cur_lnL_data);
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
      ratio += (new_lnL_data - cur_lnL_data);
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
  MCMC_Single_Param_Generic(&(tree->rates->birth_rate),
        		    tree->rates->birth_rate_min,
        		    tree->rates->birth_rate_max,
        		    tree->mcmc->num_move_birth_rate,
        		    &(tree->rates->c_lnL_times),NULL,
        		    Wrap_Lk_Times,NULL,tree->mcmc->move_type[tree->mcmc->num_move_birth_rate],NO,NULL,tree,NULL);
 

  /* phydbl cur_birth_rate,new_birth_rate; */
  /* phydbl cur_lnL_time,new_lnL_time; */
  /* phydbl cur_lnL_time_ghost,new_lnL_time_ghost; */
  /* phydbl u,alpha,ratio; */
  /* phydbl birth_rate_min,birth_rate_max; */
  /* phydbl K; */
  /* int i,n_mcmc_steps,move; */

  /* cur_birth_rate = -1.0; */
  /* new_birth_rate = -1.0; */
  /* ratio          =  0.0; */
  /* n_mcmc_steps   =  100; /\* !!!!!!!!!!! *\/ */
  /* move           = -1; */

  /* K = tree->mcmc->tune_move[tree->mcmc->num_move_birth_rate]; */

  /* cur_lnL_time = tree->rates->c_lnL_times; */
  /* new_lnL_time = tree->rates->c_lnL_times; */

  /* cur_lnL_time_ghost = tree->ghost_tree->rates->c_lnL_times; */
  /* new_lnL_time_ghost = tree->ghost_tree->rates->c_lnL_times; */
  
  /* cur_birth_rate = tree->rates->birth_rate; */

  /* birth_rate_min = tree->rates->birth_rate_min; */
  /* birth_rate_max = tree->rates->birth_rate_max; */

  /* MCMC_Make_Move(&cur_birth_rate,&new_birth_rate,birth_rate_min,birth_rate_max,&ratio,K,tree->mcmc->move_type[tree->mcmc->num_move_birth_rate]); */
  
  /* if(new_birth_rate < birth_rate_max && new_birth_rate > birth_rate_min && new_birth_rate > tree->rates->death_rate) */
  /*   { */
  /*     tree->rates->birth_rate = new_birth_rate; */
      
  /*     new_lnL_time = TIMES_Lk_Times(tree); */
  /*     ratio += (new_lnL_time - cur_lnL_time); */
      
  /*     tree->ghost_tree->rates->birth_rate = cur_birth_rate; */
  /*     cur_lnL_time_ghost = TIMES_Lk_Times(tree->ghost_tree); */

  /*     tree->ghost_tree->eval_alnL         = NO; */
  /*     tree->ghost_tree->eval_rlnL         = NO; */
  /*     tree->ghost_tree->eval_glnL         = YES; */
  /*     tree->ghost_tree->rates->birth_rate = new_birth_rate; */
  /*     TIMES_Lk_Times(tree->ghost_tree); */

  /*     i = 0; */
  /*     do */
  /*       { */
  /*         u = Uni(); */
  /*         For(move,tree->mcmc->n_moves) if(tree->mcmc->move_weight[move] > u-1.E-10) break; */
          
  /*         /\* PhyML_Printf("\n<< Move '%s' %f",tree->mcmc->move_name[move],tree->ghost_tree->rates->c_lnL_times); *\/ */

  /*         if(!strcmp(tree->mcmc->move_name[move],"tree_height")) { MCMC_Tree_Height(tree->ghost_tree); i++; } */
  /*         if(!strcmp(tree->mcmc->move_name[move],"times"))       { MCMC_Time_All(tree->ghost_tree); i++; } */
  /*         if(!strcmp(tree->mcmc->move_name[move],"spr"))         { MCMC_Prune_Regraft(tree->ghost_tree); i++; } */
  /*         if(!strcmp(tree->mcmc->move_name[move],"spr_local"))   { MCMC_Prune_Regraft_Local(tree->ghost_tree); i++; } */
  /*       } */
  /*     while(i < n_mcmc_steps); */

  /*     new_lnL_time_ghost = TIMES_Lk_Times(tree->ghost_tree); */
 
  /*     ratio += (cur_lnL_time_ghost - new_lnL_time_ghost); */

  /*     ratio = EXP(ratio); */
  /*     alpha = MIN(1.,ratio); */
      
  /*     u = Uni(); */
  /*     if(u > alpha) /\* Reject *\/ */
  /*       { */
  /*         tree->rates->birth_rate  = cur_birth_rate; */
  /*         tree->ghost_tree->rates->birth_rate  = cur_birth_rate; */
  /*         tree->rates->c_lnL_times = cur_lnL_time; */
  /*       } */
  /*     else */
  /*       { */
  /*         tree->mcmc->acc_move[tree->mcmc->num_move_birth_rate]++; */
  /*       } */
  /*   } */
  /* tree->mcmc->run_move[tree->mcmc->num_move_birth_rate]++; */
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MCMC_Death_Rate(t_tree *tree)
{
  MCMC_Single_Param_Generic(&(tree->rates->death_rate),
        		    0.0, // instead of tree->rates->death_rate_min as death rate can be equal to 0 (Yule model)
        		    tree->rates->death_rate_max,
        		    tree->mcmc->num_move_death_rate,
        		    &(tree->rates->c_lnL_times),NULL,
        		    Wrap_Lk_Times,NULL,tree->mcmc->move_type[tree->mcmc->num_move_death_rate],NO,NULL,tree,NULL);

 
  /* phydbl cur_death_rate,new_death_rate; */
  /* phydbl cur_lnL_time,new_lnL_time; */
  /* phydbl cur_lnL_time_ghost,new_lnL_time_ghost; */
  /* phydbl u,alpha,ratio; */
  /* phydbl death_rate_min,death_rate_max; */
  /* phydbl K; */
  /* int i,n_mcmc_steps,move; */

  /* cur_death_rate = -1.0; */
  /* new_death_rate = -1.0; */
  /* ratio          =  0.0; */
  /* n_mcmc_steps   =  100; /\* !!!!!!!!!!! *\/ */
  /* move           = -1; */

  /* K = tree->mcmc->tune_move[tree->mcmc->num_move_death_rate]; */

  /* cur_lnL_time = tree->rates->c_lnL_times; */
  /* new_lnL_time = tree->rates->c_lnL_times; */

  /* cur_lnL_time_ghost = tree->ghost_tree->rates->c_lnL_times; */
  /* new_lnL_time_ghost = tree->ghost_tree->rates->c_lnL_times; */
  
  /* cur_death_rate = tree->rates->death_rate; */

  /* death_rate_min = tree->rates->death_rate_min; */
  /* death_rate_max = tree->rates->death_rate_max; */

  /* MCMC_Make_Move(&cur_death_rate,&new_death_rate,death_rate_min,death_rate_max,&ratio,K,tree->mcmc->move_type[tree->mcmc->num_move_death_rate]); */
  
  /* if(new_death_rate < death_rate_max && new_death_rate > death_rate_min && new_death_rate < tree->rates->birth_rate) */
  /*   { */
  /*     tree->rates->death_rate = new_death_rate; */
      
  /*     new_lnL_time = TIMES_Lk_Times(tree); */
  /*     ratio += (new_lnL_time - cur_lnL_time); */
      
  /*     tree->ghost_tree->rates->death_rate = cur_death_rate; */
  /*     cur_lnL_time_ghost = TIMES_Lk_Times(tree->ghost_tree); */

  /*     tree->ghost_tree->eval_alnL         = NO; */
  /*     tree->ghost_tree->eval_rlnL         = NO; */
  /*     tree->ghost_tree->eval_glnL         = YES; */
  /*     tree->ghost_tree->rates->death_rate = new_death_rate; */
  /*     TIMES_Lk_Times(tree->ghost_tree); */

  /*     i = 0; */
  /*     do */
  /*       { */
  /*         u = Uni(); */
  /*         For(move,tree->mcmc->n_moves) if(tree->mcmc->move_weight[move] > u-1.E-10) break; */
          
  /*         /\* PhyML_Printf("\n>> Move '%s' %f %f b:%f d:%f-%f", *\/ */
  /*         /\*              tree->mcmc->move_name[move], *\/ */
  /*         /\*              tree->ghost_tree->rates->c_lnL_times, *\/ */
  /*         /\*              tree->rates->c_lnL_times, *\/ */
  /*         /\*              tree->rates->birth_rate, *\/ */
  /*         /\*              cur_death_rate, *\/ */
  /*         /\*              new_death_rate); *\/ */

  /*         if(!strcmp(tree->mcmc->move_name[move],"tree_height")) { MCMC_Tree_Height(tree->ghost_tree); i++; } */
  /*         if(!strcmp(tree->mcmc->move_name[move],"times"))       { MCMC_Time_All(tree->ghost_tree); i++; } */
  /*         if(!strcmp(tree->mcmc->move_name[move],"spr"))         { MCMC_Prune_Regraft(tree->ghost_tree); i++; } */
  /*         if(!strcmp(tree->mcmc->move_name[move],"spr_local"))   { MCMC_Prune_Regraft_Local(tree->ghost_tree); i++; } */
  /*       } */
  /*     while(i < n_mcmc_steps); */

  /*     new_lnL_time_ghost = TIMES_Lk_Times(tree->ghost_tree); */
 
  /*     ratio += (cur_lnL_time_ghost - new_lnL_time_ghost); */

  /*     ratio = EXP(ratio); */
  /*     alpha = MIN(1.,ratio); */
      
  /*     u = Uni(); */
  /*     if(u > alpha) /\* Reject *\/ */
  /*       { */
  /*         tree->rates->death_rate  = cur_death_rate; */
  /*         tree->ghost_tree->rates->death_rate  = cur_death_rate; */
  /*         tree->rates->c_lnL_times = cur_lnL_time; */
  /*       } */
  /*     else */
  /*       { */
  /*         tree->mcmc->acc_move[tree->mcmc->num_move_death_rate]++; */
  /*       } */
  /*   } */
  /* tree->mcmc->run_move[tree->mcmc->num_move_death_rate]++; */
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

      if(tree->rates->model == GUINDON)
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

void MCMC_Prune_Regraft(t_tree *tree)
{
  phydbl u,alpha,ratio,hr;
  phydbl t_min,t_max;
  phydbl cur_glnL, new_glnL;
  phydbl cur_alnL, new_alnL;
  phydbl regraft_t_min,regraft_t_max;
  int i,prune_idx,n_iter,n_regraft_nd,regraft_idx;
  phydbl *times;
  int rnd_dir,dir_v1,dir_v2;
  t_node *prune,*prune_daughter,*regraft_nd,*init_regraft_nd,*effective_prune;
  t_ll *regraft_nd_list;
  t_edge *target, *ori_target, *residual,*regraft_edge;

  n_iter = MAX(1,(int)(tree->n_otu/5));
  regraft_nd_list = NULL;
  times = tree->rates->nd_t;

  while(n_iter--)
    {
      tree->mcmc->run_move[tree->mcmc->num_move_spr]++;
      
      RATES_Record_Times(tree);
  
      new_glnL       = tree->rates->c_lnL_times;
      cur_glnL       = tree->rates->c_lnL_times;
      new_alnL       = tree->c_lnL;
      cur_alnL       = tree->c_lnL;
      ratio          = 0.0;
      hr             = 0.0;
      regraft_nd     = NULL;
      regraft_edge   = NULL;


      // Select prune node (any internal node)
      prune_idx = Rand_Int(tree->n_otu,2*tree->n_otu-2);

      prune = tree->a_nodes[prune_idx];

      assert(prune && prune->tax == NO);

      // Select a daughter of prune node
      dir_v1 = dir_v2 = -1;
      For(i,3) 
        if(prune->v[i] != prune->anc && prune->b[i] != tree->e_root)
          {
            if(dir_v1 < 0) dir_v1 = i;
            else           dir_v2 = i;
          }
      
      u = Uni();
      if(u < 0.5) rnd_dir = dir_v1;
      else        rnd_dir = dir_v2;
          
      prune_daughter  = prune->v[rnd_dir];
      init_regraft_nd = prune->v[rnd_dir == dir_v1 ? dir_v2 : dir_v1];      

      assert(prune_daughter->anc == prune);

      if(prune == tree->n_root)
        {
          init_regraft_nd = tree->n_root;

          if(prune_daughter == prune->v[dir_v1] && prune->v[dir_v2]->tax == YES)
            prune_daughter = prune->v[dir_v2];

          if(prune_daughter == prune->v[dir_v2] && prune->v[dir_v1]->tax == YES)
            prune_daughter = prune->v[dir_v1];

          if(prune_daughter == tree->n_root->v[1])
            effective_prune = tree->n_root->v[2];
          else
            effective_prune = tree->n_root->v[1];
        }
      else
        {
          effective_prune = prune;
        }


      /* printf("\n\n. prune: %d [%d-%d-%d] [%s-%s-%s] [%f-%f;%f] prune_daughter: %d [%f-%f;%f] prune->anc: %d [%f-%f;%f] effective: %d", */
      /*        prune->num, */
      /*        prune->v[0] ? prune->v[0]->num : -1, */
      /*        prune->v[1]->num, */
      /*        prune->v[2]->num, */
      /*        prune->v[0] ? prune->v[0]->tax ? prune->v[0]->name : "XXX" : "XXX", */
      /*        prune->v[1]->tax ? prune->v[1]->name : "XXX", */
      /*        prune->v[2]->tax ? prune->v[2]->name : "XXX", */
      /*        tree->rates->nd_t[prune->num], */
      /*        tree->rates->t_prior_min[prune->num], */
      /*        tree->rates->t_prior_max[prune->num], */
      /*        prune_daughter->num, */
      /*        tree->rates->nd_t[prune_daughter->num], */
      /*        tree->rates->t_prior_min[prune_daughter->num], */
      /*        tree->rates->t_prior_max[prune_daughter->num], */
      /*        prune->anc ? prune->anc->num : -1, */
      /*        prune->anc ? tree->rates->nd_t[prune->anc->num] : -1, */
      /*        prune->anc ? tree->rates->t_prior_min[prune->anc->num] : -1, */
      /*        prune->anc ? tree->rates->t_prior_max[prune->anc->num] : -1, */
      /*        effective_prune->num); fflush(NULL); */

      // Time and corresponding (partial) Hastings ratio
      dir_v1 = dir_v2 = -1;
      For(i,3) 
        if(prune->v[i] != prune->anc && prune->b[i] != tree->e_root)
          {
            if(dir_v1 < 0) dir_v1 = i;
            else           dir_v2 = i;
          }
      
      t_max = MIN(times[prune->v[dir_v1]->num],times[prune->v[dir_v2]->num]);
      t_min = prune->anc ? times[prune->anc->num] : 1.2*t_max;
      
      hr += LOG(1./(t_max - t_min));

      // Get the list of potential regraft nodes (oldest node on regraft edge)
      regraft_nd_list = DATE_List_Of_Regraft_Nodes(effective_prune,prune_daughter,&regraft_t_min,&regraft_t_max,tree);

      // Add root node as one can regraft above it
      /* Push_Bottom_Linked_List(tree->n_root,&regraft_nd_list); */

      // Number of regraft nodes
      n_regraft_nd = Linked_List_Len(regraft_nd_list);

      if(!(n_regraft_nd > 0)) // Should be at least 1, since original graft site is in the list
        {
          printf("\n\n. prune: %d [%d-%d-%d] [%s-%s-%s] [%f-%f;%f] prune_daughter: %d [%f-%f;%f] prune->anc: %d [%f-%f;%f] effective: %d glnL: %f",
                 prune->num,
                 prune->v[0] ? prune->v[0]->num : -1,
                 prune->v[1]->num,
                 prune->v[2]->num,
                 prune->v[0] ? prune->v[0]->tax ? prune->v[0]->name : "XXX" : "XXX",
                 prune->v[1]->tax ? prune->v[1]->name : "XXX",
                 prune->v[2]->tax ? prune->v[2]->name : "XXX",
                 tree->rates->nd_t[prune->num],
                 tree->rates->t_prior_min[prune->num],
                 tree->rates->t_prior_max[prune->num],
                 prune_daughter->num,
                 tree->rates->nd_t[prune_daughter->num],
                 tree->rates->t_prior_min[prune_daughter->num],
                 tree->rates->t_prior_max[prune_daughter->num],
                 prune->anc ? prune->anc->num : -1,
                 prune->anc ? tree->rates->nd_t[prune->anc->num] : -1,
                 prune->anc ? tree->rates->t_prior_min[prune->anc->num] : -1,
                 prune->anc ? tree->rates->t_prior_max[prune->anc->num] : -1,
                 effective_prune->num,
                 tree->rates->c_lnL_times); fflush(NULL);
          assert(FALSE);
        }
      
      // Prune
      target = residual = NULL;
      Prune_Subtree(prune,prune_daughter,&target,&residual,tree);
      ori_target = target;
      

      /* printf("\n. after prune"); */


      // Randomly select one (uniform)
      regraft_idx = Rand_Int(0,n_regraft_nd-1);
      
      regraft_nd = Linked_List_Elem(regraft_idx,regraft_nd_list);
      
      Free_Linked_List(regraft_nd_list);

      /* printf("\n. regraft: %d [%d-%d-%d] regraft->anc: %d [%f] len: %d residual_nd: %d", */
      /*        regraft_nd ? regraft_nd->num : -1, */
      /*        regraft_nd->v[0] ? regraft_nd->v[0]->num : -1, */
      /*        regraft_nd->v[1] ? regraft_nd->v[1]->num : -1, */
      /*        regraft_nd->v[2] ? regraft_nd->v[2]->num : -1, */
      /*        regraft_nd->anc ? regraft_nd->anc->num : -1, */
      /*        regraft_nd->anc ? tree->rates->nd_t[regraft_nd->anc->num] : +1., */
      /*        n_regraft_nd,residual->num); fflush(NULL); */

      /* printf("\n> root: %d %d %d edge: %d [%d %d]", */
      /*        tree->n_root->num, */
      /*        tree->n_root->v[1]->num, */
      /*        tree->n_root->v[2]->num, */
      /*        tree->e_root->num, */
      /*        tree->e_root->left->num, */
      /*        tree->e_root->rght->num); */
      /* fflush(NULL); */

      
      // Regraft edge is the one sitting above regraft_nd
      if(regraft_nd == tree->n_root->v[1] ||
         regraft_nd == tree->n_root->v[2] ||
         regraft_nd == tree->n_root) regraft_edge = tree->e_root;
      else
        {
          For(i,3) if(regraft_nd->v[i] == regraft_nd->anc) break;
          assert(i!=3);
          regraft_edge = regraft_nd->b[i];
        }

      /* printf("\n. effective_prune: %d",effective_prune->num); */
      /* fflush(NULL); */

      assert(regraft_edge);
      

      // Regraft
      Graft_Subtree(regraft_edge,
                    effective_prune,
                    residual,
                    regraft_nd,tree);

      /* printf("\n>> root: %d %d %d edge: %d [%d %d]", */
      /*        tree->n_root->num, */
      /*        tree->n_root->v[1]->num, */
      /*        tree->n_root->v[2]->num, */
      /*        tree->e_root->num, */
      /*        tree->e_root->left->num, */
      /*        tree->e_root->rght->num); */
      /* fflush(NULL); */


      // Time and corresponding (partial) Hastings ratio
      dir_v1 = dir_v2 = -1;
      For(i,3) 
        if(effective_prune->v[i] != effective_prune->anc && effective_prune->b[i] != tree->e_root)
          {
            if(dir_v1 < 0) dir_v1 = i;
            else           dir_v2 = i;
          }
      
      t_max = MIN(times[effective_prune->v[dir_v1]->num],times[effective_prune->v[dir_v2]->num]);
      t_min = effective_prune->anc ? times[effective_prune->anc->num] : 2.*t_max;
      /* printf("\n. mx: %f min:%f",t_max,t_min); */
      t_min = MAX(t_min,regraft_t_min);

      times[effective_prune->num] = Uni() * (t_max - t_min) + t_min;

      /* printf("\nxx> new prune[%d]->t:%.3f prune_anc[%d]->t:%.3f regraft->t:%.3f regraft_anc->t:%.3f [effective:%d] t_prior_min/max: [prune:[%.3f %.3f] regraft:[%.3f %.3f]] reg_t_min:%.3f reg_t_max:%.3f", */
      /*        effective_prune->num, */
      /*        times[effective_prune->num], */
      /*        effective_prune->anc->num, */
      /*        times[effective_prune->anc->num], */
      /*        times[regraft_nd->num], */
      /*        regraft_nd->anc ? times[regraft_nd->anc->num] : +1., */
      /*        effective_prune->num, */
      /*        tree->rates->t_prior_min[effective_prune->num], */
      /*        tree->rates->t_prior_max[effective_prune->num], */
      /*        tree->rates->t_prior_min[regraft_nd->num], */
      /*        tree->rates->t_prior_max[regraft_nd->num], */
      /*        regraft_t_min,regraft_t_max); */
            
      hr -= LOG(1./(t_max - t_min));

      if(!TIMES_Check_Node_Height_Ordering(tree))
        {
          printf("\n. root: %d %d %d",tree->n_root->num,tree->n_root->v[1]->num,tree->n_root->v[2]->num);
          Exit("\n");
        }

      DATE_Assign_Primary_Calibration(tree);
      new_glnL = TIMES_Lk_Times(tree); 
      ratio += (new_glnL - cur_glnL);

      if(new_glnL > UNLIKELY)
        {
          Set_Both_Sides(NO,tree);
          new_alnL = Lk(NULL,tree);
          ratio += (new_alnL - cur_alnL);
        }
      

      ratio = EXP(ratio);
      alpha = MIN(1.,ratio);

      /* Always accept move */
      if(tree->mcmc->always_yes == YES && new_glnL > UNLIKELY) alpha = 1.0;
      
      u = Uni();
            
      if(u > alpha)
        {
          // Reject
          Prune_Subtree(effective_prune,prune_daughter,&target,&residual,tree);
          /* printf("\nx root: %d %d %d edge: %d [%d %d]", */
          /*        tree->n_root->num, */
          /*        tree->n_root->v[1]->num, */
          /*        tree->n_root->v[2]->num, */
          /*        tree->e_root->num, */
          /*        tree->e_root->left->num, */
          /*        tree->e_root->rght->num); */

          Graft_Subtree(ori_target,effective_prune,residual,init_regraft_nd,tree);
          /* printf("\n< %d %d %d - %d %d",init_regraft_nd->num,ori_target->left->num,ori_target->rght->num,residual->left->num,residual->rght->num); */

          /* printf("\nxx root: %d %d %d edge: %d [%d %d]", */
          /*        tree->n_root->num, */
          /*        tree->n_root->v[1]->num, */
          /*        tree->n_root->v[2]->num, */
          /*        tree->e_root->num, */
          /*        tree->e_root->left->num, */
          /*        tree->e_root->rght->num); */
          /* fflush(NULL); */

          RATES_Reset_Times(tree);
          DATE_Assign_Primary_Calibration(tree);
          TIMES_Lk_Times(tree); 
          tree->c_lnL = cur_alnL;
          tree->rates->c_lnL_times = cur_glnL;

          /* /\* !!!!!!!!!!!!!! *\/ */
          /* new_alnL = Lk(NULL,tree); /\* Not necessary. Remove once tested *\/ */
          /* if(Are_Equal(new_alnL,cur_alnL,1.E-3) == NO) */
          /*   { */
          /*     PhyML_Printf("\n== new_alnL: %f cur_alnL: %f",new_alnL,cur_alnL); */
          /*     Generic_Exit(__FILE__,__LINE__,__FUNCTION__); */
          /*   } */

          /* /\* !!!!!!!!!!!!!! *\/ */
          /* new_glnL = TIMES_Lk_Times(tree); /\* Not necessary. Remove once tested *\/ */
          /* if(Are_Equal(new_glnL,cur_glnL,1.E-3) == NO) */
          /*   { */
          /*     PhyML_Printf("\n== new_glnL: %f cur_glnL: %f",new_glnL,cur_glnL); */
          /*     Generic_Exit(__FILE__,__LINE__,__FUNCTION__); */
          /*   } */
        }
      else
        {
          tree->mcmc->acc_move[tree->mcmc->num_move_spr]++;
        }      
    }
}
  
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MCMC_Prune_Regraft_Local(t_tree *tree)
{
  phydbl u,alpha,ratio,hr;
  phydbl t_min,t_max;
  phydbl cur_glnL, new_glnL;
  phydbl cur_alnL, new_alnL;
  phydbl regraft_t_min,regraft_t_max;
  int i,prune_idx,n_iter,n_regraft_nd,regraft_idx;
  phydbl *times;
  int rnd_dir,dir_v1,dir_v2;
  t_node *prune,*prune_daughter,*regraft_nd,*init_regraft_nd;
  t_ll *regraft_nd_list;
  t_edge *target, *ori_target, *residual,*regraft_edge;
  t_node *a, *b, *c, *d, *e;

  n_iter = MAX(1,(int)(tree->n_otu/5));
  regraft_nd_list = NULL;
  times = tree->rates->nd_t;

  while(n_iter--)
    {
      tree->mcmc->run_move[tree->mcmc->num_move_spr_local]++;
      
      RATES_Record_Times(tree);
  
      new_glnL       = tree->rates->c_lnL_times;
      cur_glnL       = tree->rates->c_lnL_times;
      new_alnL       = tree->c_lnL;
      cur_alnL       = tree->c_lnL;
      ratio          = 0.0;
      hr             = 0.0;
      regraft_nd     = NULL;

      // Select prune node (any internal node except the root)
      i = 0;
      do
        {
          prune_idx = Rand_Int(tree->n_otu,2*tree->n_otu-2);
          i++;
          assert(i<1000);
        }
      while(tree->a_nodes[prune_idx] == tree->n_root);

      prune = tree->a_nodes[prune_idx];

      assert(prune && prune->tax == NO);

      // Select a daughter of prune node
      dir_v1 = dir_v2 = -1;
      For(i,3) 
        if(prune->v[i] != prune->anc && prune->b[i] != tree->e_root)
          {
            if(dir_v1 < 0) dir_v1 = i;
            else           dir_v2 = i;
          }
      
      u = Uni();
      if(u < 0.5) rnd_dir = dir_v1;
      else        rnd_dir = dir_v2;
      
      prune_daughter  = prune->v[rnd_dir];
      init_regraft_nd = prune->v[rnd_dir == dir_v1 ? dir_v2 : dir_v1];

      // Time and corresponding (partial) Hastings ratio
      dir_v1 = dir_v2 = -1;
      For(i,3) 
        if(prune->v[i] != prune->anc && prune->b[i] != tree->e_root)
          {
            if(dir_v1 < 0) dir_v1 = i;
            else           dir_v2 = i;
          }
      
      t_max = MIN(times[prune->v[dir_v1]->num],times[prune->v[dir_v2]->num]);
      t_min = times[prune->anc->num];
      
      hr += LOG(1./(t_max - t_min));

      /*
          |
          |
          d
         / \
        /   \c
       / \
    e /   \prune_daughter
     /\
    /  \
    a   b
        */

      /* printf("\n> prune: %d prune->anc: %d [%d %d %d] daughter: %d root: %d root->v1: %d root->v2: %d", */
      /*        prune->num, */
      /*        prune->anc->num, */
      /*        prune->v[0]->num, */
      /*        prune->v[1]->num, */
      /*        prune->v[2]->num, */
      /*        prune_daughter->num, */
      /*        tree->n_root->num, */
      /*        tree->n_root->v[1]->num, */
      /*        tree->n_root->v[2]->num); */
      /* fflush(NULL); */

      // Get the list of potential regraft nodes (oldest node on regraft edge)
      regraft_nd_list = DATE_List_Of_Regraft_Nodes(prune,prune_daughter,&regraft_t_min,&regraft_t_max,tree);

      assert(regraft_nd_list);
      
      /* { */
      /*   t_ll *ll; */
      /*   t_node *n; */
      /*   ll = regraft_nd_list->head; */
      /*   do */
      /*     { */
      /*       n = ll->v; */
      /*       printf("\n> LIST %d",n->num); */
      /*       ll = ll->next;             */
      /*     } */
      /*   while(ll); */
      /* } */
      
      a = b = c = d = e = NULL;
      
      For(i,3)
        {
          if(prune->v[i] != prune->anc && prune->b[i] != tree->e_root && prune->v[i] != prune_daughter)
            {
              e = prune->v[i];
              break;
            }
        }
      
      assert(e);

      if(e->tax == NO)
        {
          For(i,3)
            {
              if(e->v[i] != e->anc && e->b[i] != tree->e_root)
                {
                  if(a == NULL) a = e->v[i];
                  else          b = e->v[i];
                }
              }
        }
      
      d = prune->anc;
      
      For(i,3)
        {
          if(d->v[i] != d->anc && d->b[i] != tree->e_root && d->v[i] != prune)
            {
              c = d->v[i];
              break;
            }
        }
      
      /* printf("\n> a: %d b: %d c: %d d; %d e: %d", */
      /*        a?a->num:-1, */
      /*        b?b->num:-1, */
      /*        c?c->num:-1, */
      /*        d?d->num:-1, */
      /*        e?e->num:-1); */

      For(i,2*tree->n_otu-1)
        {
          if(tree->a_nodes[i] != a &&
             tree->a_nodes[i] != b &&
             tree->a_nodes[i] != c &&
             tree->a_nodes[i] != d &&
             tree->a_nodes[i] != e)
            Remove_From_Linked_List(NULL,tree->a_nodes[i],&regraft_nd_list);
        }

      // Cannot regraft above root node
      Remove_From_Linked_List(NULL,tree->n_root,&regraft_nd_list);
      
      // Number of regraft nodes
      n_regraft_nd = Linked_List_Len(regraft_nd_list);
      
      assert(n_regraft_nd > 0); // Should be at least 1, since original graft site is in the list
      
      hr -= 1./n_regraft_nd;


      
      // Randomly select one (uniform)
      regraft_idx = Rand_Int(0,n_regraft_nd-1);      
      regraft_nd = Linked_List_Elem(regraft_idx,regraft_nd_list);
      Free_Linked_List(regraft_nd_list);
            
      // Prune
      target = residual = NULL;
      Prune_Subtree(prune,prune_daughter,&target,&residual,tree);
      ori_target = target;
            
      // Regraft edge is the one sitting above regraft_nd
      if(regraft_nd == tree->n_root->v[1] ||
         regraft_nd == tree->n_root->v[2]) regraft_edge = tree->e_root;
      else
        {
          For(i,3) if(regraft_nd->v[i] == regraft_nd->anc) break;
          regraft_edge = regraft_nd->b[i];
        }

      assert(regraft_edge);

      // Regraft
      Graft_Subtree(regraft_edge,prune,residual,regraft_nd,tree);

      /* printf("\n. REGRAFT: %d",regraft_nd->num); */
      
      /* printf("\n< prune: %d prune->anc: %d [%d %d %d] daughter: %d root: %d root->v1: %d root->v2: %d", */
      /*        prune->num, */
      /*        prune->anc->num, */
      /*        prune->v[0]->num, */
      /*        prune->v[1]->num, */
      /*        prune->v[2]->num, */
      /*        prune_daughter->num, */
      /*        tree->n_root->num, */
      /*        tree->n_root->v[1]->num, */
      /*        tree->n_root->v[2]->num); */
      /* fflush(NULL); */

      
      // Probability of inverse move
      // Get the list of potential regraft nodes (oldest node on regraft edge)
      regraft_nd_list = DATE_List_Of_Regraft_Nodes(prune,prune_daughter,&regraft_t_min,&regraft_t_max,tree);

      /* { */
      /*   t_ll *ll; */
      /*   t_node *n; */
      /*   ll = regraft_nd_list->head; */
      /*   do */
      /*     { */
      /*       n = ll->v; */
      /*       printf("\n< LIST %d",n->num); */
      /*       ll = ll->next;             */
      /*     } */
      /*   while(ll); */
      /* } */

      a = b = c = d = e = NULL;
      
      For(i,3)
        {
          if(prune->v[i] != prune->anc && prune->b[i] != tree->e_root && prune->v[i] != prune_daughter)
            {
              e = prune->v[i];
              break;
            }
        }
      
      assert(e);

      if(e->tax == NO)
        {
          For(i,3)
            {
              if(e->v[i] != e->anc && e->b[i] != tree->e_root)
                {
                  if(a == NULL) a = e->v[i];
                  else          b = e->v[i];
                }
              }
        }
      
      d = prune->anc;
      
      For(i,3)
        {
          if(d->v[i] != d->anc && d->b[i] != tree->e_root && d->v[i] != prune)
            {
              c = d->v[i];
              break;
            }
        }
      
      /* printf("\n< a: %d b: %d c: %d d; %d e: %d", */
      /*        a?a->num:-1, */
      /*        b?b->num:-1, */
      /*        c?c->num:-1, */
      /*        d?d->num:-1, */
      /*        e?e->num:-1); */

      For(i,2*tree->n_otu-1)
        {
          if(tree->a_nodes[i] != a &&
             tree->a_nodes[i] != b &&
             tree->a_nodes[i] != c &&
             tree->a_nodes[i] != d &&
             tree->a_nodes[i] != e)
            Remove_From_Linked_List(NULL,tree->a_nodes[i],&regraft_nd_list);
        }
          
      // Cannot regraft above root node
      Remove_From_Linked_List(NULL,tree->n_root,&regraft_nd_list);

      // Number of regraft nodes
      n_regraft_nd = Linked_List_Len(regraft_nd_list);
      
      assert(n_regraft_nd > 0); // Should be at least 1, since original graft site is in the list
      
      hr += 1./n_regraft_nd;

      Free_Linked_List(regraft_nd_list);




      // Time of regraft node and corresponding (partial) Hastings ratio
      dir_v1 = dir_v2 = -1;
      For(i,3) 
        if(prune->v[i] != prune->anc && prune->b[i] != tree->e_root)
          {
            if(dir_v1 < 0) dir_v1 = i;
            else           dir_v2 = i;
          }
      
      t_max = MIN(times[prune->v[dir_v1]->num],times[prune->v[dir_v2]->num]);
      t_min = times[prune->anc->num];

      // !!!!!!!!!!
      t_min = MAX(t_min,regraft_t_min);


      
      times[prune->num] = Uni() * (t_max - t_min) + t_min;
            
      hr -= LOG(1./(t_max - t_min));

      if(!TIMES_Check_Node_Height_Ordering(tree))
        {
          printf("\n. root: %d %d %d",tree->n_root->num,tree->n_root->v[1]->num,tree->n_root->v[2]->num);
          Exit("\n");
        }

      DATE_Assign_Primary_Calibration(tree);
      new_glnL = TIMES_Lk_Times(tree); 
      ratio += (new_glnL - cur_glnL);

      if(new_glnL > UNLIKELY)
        {
          Set_Both_Sides(NO,tree);
          new_alnL = Lk(NULL,tree);
          ratio += (new_alnL - cur_alnL);
        }
      
      ratio = EXP(ratio);
      alpha = MIN(1.,ratio);

      /* Always accept move */
      if(tree->mcmc->always_yes == YES && new_glnL > UNLIKELY) alpha = 1.0;
      
      u = Uni();
            
      if(u > alpha)
        {
          // Reject
          Prune_Subtree(prune,prune_daughter,&target,&residual,tree);
          /* printf("\n< %d %d %d - %d %d",init_regraft_nd->num,ori_target->left->num,ori_target->rght->num,residual->left->num,residual->rght->num); */
          /* printf("\nx root: %d %d %d",tree->n_root->num,tree->n_root->v[1]->num,tree->n_root->v[2]->num);           */
          fflush(NULL);
          Graft_Subtree(ori_target,prune,residual,init_regraft_nd,tree);

          RATES_Reset_Times(tree);
          DATE_Assign_Primary_Calibration(tree);
          TIMES_Lk_Times(tree); 
          tree->c_lnL = cur_alnL;
          tree->rates->c_lnL_times = cur_glnL;

          /* DATE_Assign_Primary_Calibration(tree); */
          /* DATE_Update_T_Prior_MinMax(tree); */

          /* /\* !!!!!!!!!!!!!! *\/ */
          /* new_alnL = Lk(NULL,tree); /\* Not necessary. Remove once tested *\/ */
          /* if(Are_Equal(new_alnL,cur_alnL,1.E-3) == NO) */
          /*   { */
          /*     PhyML_Printf("\n== new_alnL: %f cur_alnL: %f",new_alnL,cur_alnL); */
          /*     Generic_Exit(__FILE__,__LINE__,__FUNCTION__); */
          /*   } */

          /* /\* !!!!!!!!!!!!!! *\/ */
          /* new_glnL = TIMES_Lk_Times(tree); /\* Not necessary. Remove once tested *\/ */
          /* if(Are_Equal(new_glnL,cur_glnL,1.E-3) == NO) */
          /*   { */
          /*     PhyML_Printf("\n== new_glnL: %f cur_glnL: %f",new_glnL,cur_glnL); */
          /*     Generic_Exit(__FILE__,__LINE__,__FUNCTION__); */
          /*   } */
        }
      else
        {
          tree->mcmc->acc_move[tree->mcmc->num_move_spr_local]++;
        }      
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

  mcmc->num_move_times                    = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_root_time                = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_nu                       = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_clock_r                  = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_tree_height              = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_time_slice               = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_subtree_height           = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_kappa                    = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_tree_rates               = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_subtree_rates            = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_updown_nu_cr             = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_ras                      = mcmc->n_moves; mcmc->n_moves += (tree->mod ? 2*tree->mod->ras->n_catg : 1);
  mcmc->num_move_updown_t_cr              = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_cov_rates                = mcmc->n_moves; mcmc->n_moves += (tree->mod && tree->mod->m4mod ? 2*tree->mod->m4mod->n_h : 1);
  mcmc->num_move_cov_switch               = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_birth_rate               = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_spr                      = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_spr_local                = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_death_rate               = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_updown_t_br              = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_jump_calibration         = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_geo_lambda               = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_geo_sigma                = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_geo_tau                  = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_geo_updown_tau_lbda      = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_geo_updown_lbda_sigma    = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_geo_dum                  = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_phyrex_lbda              = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_phyrex_mu                = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_phyrex_rad               = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_phyrex_indel_disk        = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_phyrex_move_disk_ud      = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_phyrex_swap_disk         = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_phyrex_indel_hit         = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_phyrex_spr               = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_phyrex_scale_times       = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_phyrex_ldscape_lim       = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_phyrex_sigsq             = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_phyrex_sim               = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_phyrex_traj              = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_phyrex_lbda_times        = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_phyrex_indel_disk_serial = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_phyrex_sim_plus          = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_phyrex_indel_hit_serial  = mcmc->n_moves; mcmc->n_moves += 1;
  
  mcmc->num_move_phyrex_ldsk_given_disk   = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_phyrex_disk_given_ldsk   = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_phyrex_ldsk_and_disk     = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_phyrex_ldsk_multi        = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_phyrex_disk_multi        = mcmc->n_moves; mcmc->n_moves += 1;

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
  strcpy(mcmc->move_name[mcmc->num_move_times],"times");
  strcpy(mcmc->move_name[mcmc->num_move_root_time],"root_time");
  strcpy(mcmc->move_name[mcmc->num_move_nu],"nu");
  strcpy(mcmc->move_name[mcmc->num_move_clock_r],"clock");
  strcpy(mcmc->move_name[mcmc->num_move_tree_height],"tree_height");
  strcpy(mcmc->move_name[mcmc->num_move_time_slice],"time_slice");
  strcpy(mcmc->move_name[mcmc->num_move_subtree_height],"subtree_height");
  strcpy(mcmc->move_name[mcmc->num_move_kappa],"kappa");
  strcpy(mcmc->move_name[mcmc->num_move_spr],"spr");
  strcpy(mcmc->move_name[mcmc->num_move_spr_local],"spr_local");
  strcpy(mcmc->move_name[mcmc->num_move_tree_rates],"tree_rates");
  strcpy(mcmc->move_name[mcmc->num_move_subtree_rates],"subtree_rates");
  strcpy(mcmc->move_name[mcmc->num_move_updown_nu_cr],"updown_nu_cr");
  for(i=mcmc->num_move_ras;i<mcmc->num_move_ras+(tree->mod ? 2*tree->mod->ras->n_catg : 1);i++) 
    strcpy(mcmc->move_name[i],"ras");  
  strcpy(mcmc->move_name[mcmc->num_move_updown_t_cr],"updown_t_cr");
  for(i=mcmc->num_move_cov_rates;i<mcmc->num_move_cov_rates+(tree->mod && tree->mod->m4mod ? 2*tree->mod->m4mod->n_h : 1);i++) 
    strcpy(mcmc->move_name[i],"cov_rates");  
  strcpy(mcmc->move_name[mcmc->num_move_cov_switch],"cov_switch");
  strcpy(mcmc->move_name[mcmc->num_move_birth_rate],"birth_rate");
  strcpy(mcmc->move_name[mcmc->num_move_death_rate],"death_rate");
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
  strcpy(mcmc->move_name[mcmc->num_move_phyrex_move_disk_ud],"phyrex_move_disk_ud");
  strcpy(mcmc->move_name[mcmc->num_move_phyrex_swap_disk],"phyrex_swap_disk");
  strcpy(mcmc->move_name[mcmc->num_move_phyrex_indel_hit],"phyrex_indel_hit");
  strcpy(mcmc->move_name[mcmc->num_move_phyrex_spr],"phyrex_spr");
  strcpy(mcmc->move_name[mcmc->num_move_phyrex_scale_times],"phyrex_scale_times");
  strcpy(mcmc->move_name[mcmc->num_move_phyrex_ldscape_lim],"phyrex_ldscape_lim");
  strcpy(mcmc->move_name[mcmc->num_move_phyrex_sigsq],"phyrex_sigsq");
  strcpy(mcmc->move_name[mcmc->num_move_phyrex_sim],"phyrex_sim");
  strcpy(mcmc->move_name[mcmc->num_move_phyrex_traj],"phyrex_traj");
  strcpy(mcmc->move_name[mcmc->num_move_phyrex_lbda_times],"phyrex_lbda_times");
  strcpy(mcmc->move_name[mcmc->num_move_phyrex_indel_disk_serial],"phyrex_indel_disk_serial");
  strcpy(mcmc->move_name[mcmc->num_move_phyrex_sim_plus],"phyrex_sim_plus");
  strcpy(mcmc->move_name[mcmc->num_move_phyrex_indel_hit_serial],"phyrex_indel_hit_serial");
  strcpy(mcmc->move_name[mcmc->num_move_phyrex_ldsk_multi],"phyrex_ldsk_multi");
  strcpy(mcmc->move_name[mcmc->num_move_phyrex_disk_multi],"phyrex_disk_multi");
  strcpy(mcmc->move_name[mcmc->num_move_phyrex_ldsk_and_disk],"phyrex_ldsk_and_disk");
  strcpy(mcmc->move_name[mcmc->num_move_phyrex_disk_given_ldsk],"phyrex_disk_given_ldsk");
  strcpy(mcmc->move_name[mcmc->num_move_phyrex_ldsk_given_disk],"phyrex_ldsk_given_disk");


  For(i,mcmc->n_moves) mcmc->move_type[i] = -1;
  
  if(tree->rates && tree->rates->model_log_rates == YES)
    for(i=mcmc->num_move_br_r;i<mcmc->num_move_br_r+2*tree->n_otu-2;i++) mcmc->move_type[i] = MCMC_MOVE_RANDWALK_NORMAL;
  else
    for(i=mcmc->num_move_br_r;i<mcmc->num_move_br_r+2*tree->n_otu-2;i++) mcmc->move_type[i] = MCMC_MOVE_SCALE_THORNE;

  for(i=mcmc->num_move_nd_r;i<mcmc->num_move_nd_r+2*tree->n_otu-1;i++) mcmc->move_type[i] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_times] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_root_time] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_nu] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_clock_r] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_tree_height] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_time_slice] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_subtree_height] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_kappa] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_tree_rates] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_subtree_rates] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_updown_nu_cr] = MCMC_MOVE_RANDWALK_NORMAL;
  for(i=mcmc->num_move_ras;i<mcmc->num_move_ras+(tree->mod ? 2*tree->mod->ras->n_catg : 1);i++) mcmc->move_type[i] = MCMC_MOVE_RANDWALK_NORMAL;  
  mcmc->move_type[mcmc->num_move_updown_t_cr] = MCMC_MOVE_SCALE_THORNE;
  for(i=mcmc->num_move_cov_rates;i<mcmc->num_move_cov_rates+(tree->mod && tree->mod->m4mod ? 2*tree->mod->m4mod->n_h : 1);i++) mcmc->move_type[i] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_cov_switch] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_birth_rate] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_death_rate] = MCMC_MOVE_SCALE_THORNE;
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
  
  for(i=mcmc->num_move_br_r;i<mcmc->num_move_br_r+2*tree->n_otu-2;i++) mcmc->move_weight[i] = (phydbl)(1./tree->n_otu); /* Rates */
  for(i=mcmc->num_move_nd_r;i<mcmc->num_move_nd_r+2*tree->n_otu-1;i++) mcmc->move_weight[i] = 0.0; /* Node rates */
  mcmc->move_weight[mcmc->num_move_times]            = 1.0;
  mcmc->move_weight[mcmc->num_move_root_time]        = 0.0;
  mcmc->move_weight[mcmc->num_move_clock_r]          = 1.0;
  mcmc->move_weight[mcmc->num_move_tree_height]      = 1.0;
  mcmc->move_weight[mcmc->num_move_time_slice]       = 1.0;
  mcmc->move_weight[mcmc->num_move_subtree_height]   = 1.0;
  mcmc->move_weight[mcmc->num_move_nu]               = 1.0;
  mcmc->move_weight[mcmc->num_move_kappa]            = 0.5;
  mcmc->move_weight[mcmc->num_move_spr]              = 3.0;
  mcmc->move_weight[mcmc->num_move_spr_local]        = 2.0;
  mcmc->move_weight[mcmc->num_move_tree_rates]       = 1.0;
  mcmc->move_weight[mcmc->num_move_subtree_rates]    = 0.5;
  mcmc->move_weight[mcmc->num_move_updown_nu_cr]     = 0.0;
  for(i=mcmc->num_move_ras;i<mcmc->num_move_ras+(tree->mod ? 2*tree->mod->ras->n_catg : 1);i++) mcmc->move_weight[i] = 0.5*(1./(tree->mod ? (phydbl)tree->mod->ras->n_catg : 1));
  mcmc->move_weight[mcmc->num_move_updown_t_cr]      = 1.0; /* Does not seem to work well (does not give uniform prior on root height
  							      when sampling from prior) */
  for(i=mcmc->num_move_cov_rates;i<mcmc->num_move_cov_rates+(tree->mod && tree->mod->m4mod ? 2*tree->mod->m4mod->n_h : 1);i++) mcmc->move_weight[i] = 0.5*(1./(tree->mod && tree->mod->m4mod ? (phydbl)tree->mod->m4mod->n_h : 1));
  mcmc->move_weight[mcmc->num_move_cov_switch]            = 1.0;
  mcmc->move_weight[mcmc->num_move_birth_rate]            = 2.0;
  mcmc->move_weight[mcmc->num_move_death_rate]            = 2.0;
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

  mcmc->move_weight[mcmc->num_move_phyrex_lbda]                  = 5.0;
  mcmc->move_weight[mcmc->num_move_phyrex_mu]                    = 8.0;
  mcmc->move_weight[mcmc->num_move_phyrex_rad]                   = 5.0;
  mcmc->move_weight[mcmc->num_move_phyrex_sigsq]                 = 0.0;
  mcmc->move_weight[mcmc->num_move_phyrex_indel_disk]            = 2.0;
  mcmc->move_weight[mcmc->num_move_phyrex_indel_hit]             = 2.0;
  mcmc->move_weight[mcmc->num_move_phyrex_move_disk_ud]          = 5.0;
  mcmc->move_weight[mcmc->num_move_phyrex_swap_disk]             = 1.0;
  mcmc->move_weight[mcmc->num_move_phyrex_spr]                   = 1.0;
  mcmc->move_weight[mcmc->num_move_phyrex_scale_times]           = 2.0;
  mcmc->move_weight[mcmc->num_move_phyrex_ldscape_lim]           = 0.0;
  mcmc->move_weight[mcmc->num_move_phyrex_sim]                   = 2.0;
  mcmc->move_weight[mcmc->num_move_phyrex_traj]                  = 3.0;
  mcmc->move_weight[mcmc->num_move_phyrex_lbda_times]            = 1.0;
  mcmc->move_weight[mcmc->num_move_phyrex_sim_plus]              = 0.0;
  mcmc->move_weight[mcmc->num_move_phyrex_indel_disk_serial]     = 2.0;
  mcmc->move_weight[mcmc->num_move_phyrex_indel_hit_serial]      = 2.0;

  mcmc->move_weight[mcmc->num_move_phyrex_ldsk_given_disk]       = 1.0;
  mcmc->move_weight[mcmc->num_move_phyrex_disk_given_ldsk]       = 1.0;
  mcmc->move_weight[mcmc->num_move_phyrex_ldsk_and_disk]         = 1.0;
  mcmc->move_weight[mcmc->num_move_phyrex_ldsk_multi]            = 2.0;
  mcmc->move_weight[mcmc->num_move_phyrex_disk_multi]            = 2.0;

# else

  mcmc->move_weight[mcmc->num_move_phyrex_lbda]                  = 0.0;
  mcmc->move_weight[mcmc->num_move_phyrex_mu]                    = 0.0;
  mcmc->move_weight[mcmc->num_move_phyrex_rad]                   = 0.0;
  mcmc->move_weight[mcmc->num_move_phyrex_sigsq]                 = 0.0;
  mcmc->move_weight[mcmc->num_move_phyrex_indel_disk]            = 0.0;
  mcmc->move_weight[mcmc->num_move_phyrex_indel_hit]             = 0.0;
  mcmc->move_weight[mcmc->num_move_phyrex_move_disk_ud]          = 0.0;
  mcmc->move_weight[mcmc->num_move_phyrex_swap_disk]             = 0.0;
  mcmc->move_weight[mcmc->num_move_phyrex_spr]                   = 0.0;
  mcmc->move_weight[mcmc->num_move_phyrex_scale_times]           = 0.0;
  mcmc->move_weight[mcmc->num_move_phyrex_ldscape_lim]           = 0.0;
  mcmc->move_weight[mcmc->num_move_phyrex_sim]                   = 0.0;
  mcmc->move_weight[mcmc->num_move_phyrex_traj]                  = 0.0;
  mcmc->move_weight[mcmc->num_move_phyrex_lbda_times]            = 0.0;
  mcmc->move_weight[mcmc->num_move_phyrex_sim_plus]              = 0.0;
  mcmc->move_weight[mcmc->num_move_phyrex_indel_disk_serial]     = 0.0;
  mcmc->move_weight[mcmc->num_move_phyrex_indel_hit_serial]      = 0.0;

  mcmc->move_weight[mcmc->num_move_phyrex_ldsk_given_disk]       = 0.0;
  mcmc->move_weight[mcmc->num_move_phyrex_disk_given_ldsk]       = 0.0;
  mcmc->move_weight[mcmc->num_move_phyrex_ldsk_and_disk]         = 0.0;
  mcmc->move_weight[mcmc->num_move_phyrex_ldsk_multi]            = 0.0;
  mcmc->move_weight[mcmc->num_move_phyrex_disk_multi]            = 0.0;
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
  mcmc->sampled_val[mcmc->num_move_death_rate*mcmc->sample_size+mcmc->sample_num]  = tree->rates->death_rate;

  /* For(i,2*tree->n_otu-2) */
  /*   mcmc->sampled_val[(mcmc->num_move_br_r+i)*mcmc->sample_size+mcmc->sample_num] = tree->rates->br_r[i]; */
  

  /* For(i,2*tree->n_otu-1) */
  /*   mcmc->sampled_val[(mcmc->num_move_nd_r+i)*mcmc->sample_size+mcmc->sample_num] = tree->rates->nd_r[i]; */

  mcmc->sampled_val[mcmc->num_move_geo_tau*mcmc->sample_size+mcmc->sample_num]      = tree->geo ? tree->geo->tau   : -1.;
  mcmc->sampled_val[mcmc->num_move_geo_lambda*mcmc->sample_size+mcmc->sample_num]   = tree->geo ? tree->geo->lbda  : -1.;
  mcmc->sampled_val[mcmc->num_move_geo_sigma*mcmc->sample_size+mcmc->sample_num]    = tree->geo ? tree->geo->sigma : -1.;
  mcmc->sampled_val[mcmc->num_move_geo_dum*mcmc->sample_size+mcmc->sample_num]      = tree->geo ? tree->geo->dum   : -1.;  
  #ifdef PHYREX
  mcmc->sampled_val[mcmc->num_move_phyrex_lbda*mcmc->sample_size+mcmc->sample_num]  = tree->mmod ? tree->mmod->lbda               : -1.;
  mcmc->sampled_val[mcmc->num_move_phyrex_mu*mcmc->sample_size+mcmc->sample_num]    = tree->mmod ? PHYREX_Neighborhood_Size(tree) : -1.;  
  mcmc->sampled_val[mcmc->num_move_phyrex_sigsq*mcmc->sample_size+mcmc->sample_num] = tree->mmod ? PHYREX_Update_Sigsq(tree)      : -1.;
  mcmc->sampled_val[mcmc->num_move_phyrex_rad*mcmc->sample_size+mcmc->sample_num]   = tree->mmod ? tree->mmod->rad                : -1.;
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
		if(tree->io->lk_approx == EXACT) Update_P_Lk(tree,d->b[i],d);
		MCMC_Slice_One_Rate(d,d->v[i],YES,tree);
	      }
	}
      if(tree->io->lk_approx == EXACT) Update_P_Lk(tree,b,d);
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
  /* phydbl u,alpha,ratio; */
  /* phydbl cur_glnL, new_glnL, hr; */
  /* phydbl ori_mu; */
  /* phydbl ori_beta,new_beta; */
  /* phydbl K; */

  /* tree->mcmc->run_move[tree->mcmc->num_move_phyrex_mu]++; */

  /* new_glnL       = UNLIKELY; */
  /* cur_glnL       = tree->mmod->c_lnL; */
  /* hr             = 0.0; */
  /* ratio          = 0.0; */
  /* ori_mu         = tree->mmod->mu; */
  /* ori_beta       = (1./tree->mmod->mu)-1.; */
  /* K              = tree->mcmc->tune_move[tree->mcmc->num_move_phyrex_mu]; */

  /* u = Uni(); */
  /* new_beta = ori_beta * EXP(K*(u-.5)); */
  /* hr += LOG(new_beta/ori_beta); */

  /* tree->mmod->mu = 1./(new_beta+1.); */

  /* new_glnL = PHYREX_Lk(tree); */
  
  /* ratio += (new_glnL - cur_glnL); */
  /* ratio += hr; */
  
  /* ratio = EXP(ratio); */
  /* alpha = MIN(1.,ratio); */
  
  /* /\* Always accept move *\/ */
  /* if(tree->mcmc->always_yes == YES && new_glnL > UNLIKELY) alpha = 1.0; */

  /* u = Uni(); */
  

  /* if(u > alpha) /\* Reject *\/ */
  /*   { */
  /*     tree->mmod->mu = ori_mu; */

  /*     tree->mmod->c_lnL = cur_glnL; */

  /*     /\* new_glnL = MIGREP_Lk(tree); *\/ */
  /*     /\* if(Are_Equal(new_glnL,cur_glnL,1.E-3) == NO) *\/ */
  /*     /\*   { *\/ */
  /*     /\*     PhyML_Printf("\n. new_glnL: %f cur_glnL: %f",new_glnL,cur_glnL); *\/ */
  /*     /\*     Generic_Exit(__FILE__,__LINE__,__FUNCTION__); *\/ */
  /*     /\*   } *\/ */
  /*   } */
  /* else */
  /*   { */
  /*     tree->mcmc->acc_move[tree->mcmc->num_move_phyrex_mu]++; */
  /*     /\* printf("- Accept"); *\/ */
  /*   } */

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


  /* phydbl u,alpha,ratio; */
  /* phydbl cur_glnL, new_glnL, hr; */
  /* phydbl cur_rad, new_rad; */
  /* phydbl rate; */

  /* tree->mcmc->run_move[tree->mcmc->num_move_phyrex_rad]++; */

  /* new_glnL = tree->mmod->c_lnL; */
  /* cur_glnL = tree->mmod->c_lnL; */
  /* hr       = 0.0; */
  /* ratio    = 0.0; */
  /* cur_rad  = tree->mmod->rad; */
  /* new_rad  = tree->mmod->rad; */
  /* rate     = 1.0; */

  /* new_rad = Rexp_Trunc(rate, */
  /*                      tree->mmod->min_rad, */
  /*                      tree->mmod->max_rad); */

  /* tree->mmod->rad = new_rad; */

  /* hr += rate*(new_rad - cur_rad); */

  /* new_glnL = PHYREX_Lk(tree); */


  /* ratio += (new_glnL - cur_glnL); */
  /* ratio += hr; */
  
  /* ratio = EXP(ratio); */
  /* alpha = MIN(1.,ratio); */
  
  /* /\* Always accept move *\/ */
  /* if(tree->mcmc->always_yes == YES && new_glnL > UNLIKELY) alpha = 1.0; */

  /* u = Uni(); */
  
  /* if(u > alpha) /\* Reject *\/ */
  /*   { */
  /*     tree->mmod->rad = cur_rad; */
  /*     PHYREX_Lk(tree); */
  /*   } */
  /* else */
  /*   { */
  /*     tree->mcmc->acc_move[tree->mcmc->num_move_phyrex_rad]++; */
  /*   } */
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
  int cur_n_disk, new_n_disk;
  int cur_n_int;
  int cur_n_hit;
  phydbl hr;
  phydbl cur_lbda,new_lbda;
  phydbl cur_mu,new_mu;
  phydbl cur_rad,new_rad;
  phydbl T;
  
  T = PHYREX_Tree_Height(tree);
  T = FABS(T);

  hr = 0.0;

  new_lbda = tree->mmod->lbda;
  cur_lbda = tree->mmod->lbda;
  new_mu = tree->mmod->mu;
  cur_mu = tree->mmod->mu;
  new_rad = tree->mmod->rad;
  cur_rad = tree->mmod->rad;

  /* new_lbda = cur_lbda * EXP(0.5*(Uni()-.5)); */
  /* hr += LOG(new_lbda/cur_lbda); */

  /* new_mu = cur_mu * EXP(0.5*(Uni()-.5)); */
  /* hr += LOG(new_mu/cur_mu); */

  /* new_rad = cur_rad * EXP(0.5*(Uni()-.5)); */
  /* hr += LOG(new_rad/cur_rad); */

  if(new_rad > tree->mmod->max_rad || new_rad < tree->mmod->min_rad)     return; 
  if(new_mu > tree->mmod->max_mu || new_mu < tree->mmod->min_mu)         return; 
  if(new_lbda > tree->mmod->max_lbda || new_lbda < tree->mmod->min_lbda) return; 

  tree->mmod->lbda = new_lbda;
  tree->mmod->mu   = new_mu;
  tree->mmod->rad  = new_rad;
  
  cur_n_hit = PHYREX_Total_Number_Of_Hit_Disks(tree);
  cur_n_int = PHYREX_Total_Number_Of_Intervals(tree);

  cur_n_disk = cur_n_int - cur_n_hit;

  if(cur_n_hit == 0)
    {
      new_n_disk = (int)Rpois(cur_n_disk+1);
      
      hr += Dpois(cur_n_disk,new_n_disk+1,YES);
      hr -= Dpois(new_n_disk,cur_n_disk+1,YES);
    }
  else
    {
      new_n_disk = (int)Rpois(T*new_lbda*(phydbl)cur_n_disk/(cur_n_hit+cur_n_disk));
      
      hr += Dpois(cur_n_disk,T*cur_lbda*(phydbl)new_n_disk/(cur_n_hit+new_n_disk),YES);
      hr -= Dpois(new_n_disk,T*new_lbda*(phydbl)cur_n_disk/(cur_n_hit+cur_n_disk),YES);
    }

  if(new_n_disk < cur_n_disk) MCMC_PHYREX_Delete_Disk(hr, cur_n_disk - new_n_disk , cur_lbda, cur_mu, cur_rad, tree);
  else                        MCMC_PHYREX_Insert_Disk(hr, new_n_disk - cur_n_disk,  cur_lbda, cur_mu, cur_rad, tree);  
}
#endif

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
/* Insert or delete a disk that does not affect any ldisk */
#ifdef PHYREX
void MCMC_PHYREX_Delete_Disk(phydbl hr, int n_delete_disks, phydbl cur_lbda, phydbl cur_mu, phydbl cur_rad, t_tree *tree)
{
  phydbl u,alpha,ratio;
  phydbl cur_glnL, new_glnL;
  phydbl T;
  t_dsk  *disk,**target_disk,**valid_disks;
  int i,j,block,n_valid_disks,*permut;

  tree->mcmc->run_move[tree->mcmc->num_move_phyrex_indel_disk]++;

  if(n_delete_disks == 0) 
    {
      tree->mmod->lbda = cur_lbda;
      tree->mmod->rad = cur_rad;
      tree->mmod->mu = cur_mu;
      return; 
    }

  valid_disks     = NULL;
  disk            = NULL;
  new_glnL        = tree->mmod->c_lnL;
  cur_glnL        = tree->mmod->c_lnL;
  ratio           = 0.0;
  block           = 100;

  disk = tree->disk->prev;

  if(!disk->prev) 
    {
      tree->mmod->lbda = cur_lbda;
      tree->mmod->rad = cur_rad;
      tree->mmod->mu = cur_mu;
      return; 
    }

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

  if(!n_valid_disks)
    {
      tree->mmod->lbda = cur_lbda;
      tree->mmod->rad = cur_rad;
      tree->mmod->mu = cur_mu;
      return; 
    }


  if(n_valid_disks - n_delete_disks < 0) 
    {
      tree->mmod->lbda = cur_lbda;
      tree->mmod->rad = cur_rad;
      tree->mmod->mu = cur_mu;
      return; 
    }

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

  new_glnL = PHYREX_Lk(tree);
  ratio += (new_glnL - cur_glnL);
  ratio += hr;
  
  ratio = EXP(ratio);
  alpha = MIN(1.,ratio);

  /* Always accept move */
  if(tree->mcmc->always_yes == YES && new_glnL > UNLIKELY) alpha = 1.0;

  u = Uni();
  
  /* printf("\n- Delete new_glnL: %f [%f] hr: %f u:%f alpha: %f",new_glnL,cur_glnL,hr,u,alpha); */

  if(u > alpha) /* Reject */
    {
      tree->mmod->lbda = cur_lbda;
      tree->mmod->rad = cur_rad;
      tree->mmod->mu = cur_mu;

      For(j,n_delete_disks) PHYREX_Insert_Disk(target_disk[j],tree);

      if(tree->mmod->safe_phyrex == YES)
        {
          new_glnL = PHYREX_Lk(tree);
          if(Are_Equal(new_glnL,cur_glnL,1.E-3) == NO)
            {
              PhyML_Printf("\n== new_glnL: %f cur_glnL: %f",new_glnL,cur_glnL);
              Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
            }
        }
      else
        {
          PHYREX_Update_Lindisk_List(tree);
          tree->mmod->c_lnL = cur_glnL;
        }
    }
  else
    {
      For(j,n_delete_disks) Free_Disk(target_disk[j]);
      tree->mcmc->acc_move[tree->mcmc->num_move_phyrex_indel_disk]++;
    }

  Free(valid_disks);
  Free(target_disk);
  Free(permut);
}
#endif

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

#ifdef PHYREX
void MCMC_PHYREX_Insert_Disk(phydbl hr, int n_insert_disks, phydbl cur_lbda, phydbl cur_mu, phydbl cur_rad, t_tree *tree)
{
  t_dsk  *disk,**new_disk,**target_disk;
  phydbl T,t;
  phydbl cur_glnL, new_glnL;
  phydbl u,alpha,ratio;
  int i,j,n_valid_disks;

  tree->mcmc->run_move[tree->mcmc->num_move_phyrex_indel_disk]++;

  if(n_insert_disks == 0) 
    {
      tree->mmod->lbda = cur_lbda;
      tree->mmod->rad = cur_rad;
      tree->mmod->mu = cur_mu;
      return; 
    }

  disk            = NULL;
  new_glnL        = tree->mmod->c_lnL;
  cur_glnL        = tree->mmod->c_lnL;

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
      disk = tree->disk->prev;
      while(disk && disk->time > t) disk = disk->prev;
      assert(disk->next);
      target_disk[j] = disk->next;
      
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

  new_glnL = PHYREX_Lk(tree);
  ratio = (new_glnL - cur_glnL);
  ratio += hr;
  
  ratio = EXP(ratio);
  alpha = MIN(1.,ratio);
  
  /* Always accept move */
  if(tree->mcmc->always_yes == YES && new_glnL > UNLIKELY) alpha = 1.0;

  u = Uni();
  
  /* printf("\n+ Insert new_glnL: %f [%f] hr: %f u: %f alpha: %f",new_glnL,cur_glnL,hr,u,alpha); */

  if(u > alpha) /* Reject */
    {
      tree->mmod->lbda = cur_lbda;
      tree->mmod->rad = cur_rad;
      tree->mmod->mu = cur_mu;
      
      /* printf("\nI Reject"); */

      For(j,n_insert_disks)
        {
          PHYREX_Remove_Disk(new_disk[j]);
          Free_Disk(new_disk[j]);
        }

      if(tree->mmod->safe_phyrex == YES)
        {
          new_glnL = PHYREX_Lk(tree);
          if(Are_Equal(new_glnL,cur_glnL,1.E-3) == NO)
            {
              PhyML_Printf("\n== new_glnL: %f cur_glnL: %f",new_glnL,cur_glnL);
              Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
            }
        }
      else
        {
          PHYREX_Update_Lindisk_List(tree);
          tree->mmod->c_lnL = cur_glnL;
        }

      /* printf("\nI Reject %f",tree->mmod->lbda); */
    }
  else
    {
      /* printf("\nI Accept %f",tree->mmod->lbda); */
      tree->mcmc->acc_move[tree->mcmc->num_move_phyrex_indel_disk]++;
    }

  Free(target_disk);
  Free(new_disk);
}
#endif

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

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
  phydbl max, min, log_lbda;
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
  log_lbda      = tree->mmod->lbda;

  if(tree->disk->next) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  disk = tree->disk->prev;
  n_all_disks = 0;
  do
    {
      /* if(1) */
      if(disk->ldsk && disk->ldsk->n_next > 1) /* Moving disks other than coalescent is pointless */
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
  
  n_move_disks = Rand_Int(1,(int)(n_all_disks));
  
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
          phydbl new_plusmin, cur_plusmin;

          max = target_disk[i]->next->time;

          cur_plusmin = FABS(ori_time[i] - max);
          new_plusmin = Rexp(1./cur_plusmin);

          new_time = max - new_plusmin;

          hr += LOG(Dexp(cur_plusmin,1./new_plusmin));
          hr -= LOG(Dexp(new_plusmin,1./cur_plusmin));

          new_glnL -= (log_lbda - FABS(ori_time[i]-max)*tree->mmod->lbda);
          new_glnL += (log_lbda - FABS(new_time   -max)*tree->mmod->lbda);
        }
            
      target_disk[i]->time = new_time;
    }
  

  tree->mmod->c_lnL = new_glnL;

  new_alnL = Lk(NULL,tree);

  ratio += (new_alnL - cur_alnL);
  ratio += (new_glnL - cur_glnL);
  ratio += hr;
  
  ratio = EXP(ratio);
  alpha = MIN(1.,ratio);
  
  /* Always accept move */
  if(tree->mcmc->always_yes == YES && new_glnL > UNLIKELY) alpha = 1.0;

  u = Uni();
  
  /* printf("\n- Move disk new_glnL: %f [%f] hr: %f u:%f alpha: %f",new_alnL,cur_alnL,hr,u,alpha); */

  if(u > alpha) /* Reject */
    {
      /* printf("- Reject"); */
      
      For(i,n_move_disks) target_disk[i]->time = ori_time[i]; 

      if(tree->mmod->safe_phyrex == YES)
        {
          new_alnL = Lk(NULL,tree);
              
          if(Are_Equal(new_alnL,cur_alnL,1.E-3) == NO)
            {
              PhyML_Printf("\n== new_alnL: %f cur_alnL: %f",new_alnL,cur_alnL);
              Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
            }
            
          
          new_glnL = PHYREX_Lk(tree);
          
          if(Are_Equal(new_glnL,cur_glnL,1.E-3) == NO)
            {
              PhyML_Printf("\n== new_glnL: %f cur_glnL: %f",new_glnL,cur_glnL);
              Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
            }
        }
      else
        {
          tree->mmod->c_lnL = cur_glnL;
          tree->c_lnL       = cur_alnL;
        }
    }
  else
    {
      /* printf("- Accept"); */
    }  

  Free(target_disk);
  Free(all_disks);
  Free(ori_time);
  Free(permut);

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
  cur_alnL        = tree->c_lnL;
  new_alnL        = tree->c_lnL;
  new_glnL        = tree->mmod->c_lnL;
  cur_glnL        = tree->mmod->c_lnL;
  hr              = 0.0;
  ratio           = 0.0;
  K               = tree->mcmc->tune_move[tree->mcmc->num_move_phyrex_scale_times];
  cur_lbda        = tree->mmod->lbda;
  new_lbda        = tree->mmod->lbda;

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
  new_alnL = Lk(NULL,tree);

  ratio += (new_alnL - cur_alnL);
  ratio += (new_glnL - cur_glnL);
  ratio += hr;
  
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

      /* printf("- Reject");       */
      disk = tree->disk->prev;
      do
        {
          disk->time /= scale_fact_times;
          disk = disk->prev;
        }
      while(disk);


      if(tree->mmod->safe_phyrex == YES)
        {
          new_alnL = Lk(NULL,tree); /* Not necessary. Remove once tested */
          
          if(Are_Equal(new_alnL,cur_alnL,1.E-3) == NO)
            {
              PhyML_Printf("\n== new_alnL: %f cur_alnL: %f lbda: %f",new_alnL,cur_alnL,tree->mmod->lbda);
              Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
            }
            
          
          new_glnL = PHYREX_Lk(tree); /* Same here */
          
          if(Are_Equal(new_glnL,cur_glnL,1.E-3) == NO)
            {
              PhyML_Printf("\n== new_glnL: %f cur_glnL: %f lbda: %f",new_glnL,cur_glnL,tree->mmod->lbda);
              Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
            }
        }
      else
        {
          tree->mmod->c_lnL = cur_glnL;
          tree->c_lnL       = cur_alnL;
        }
    }
  else
    {
      /* printf("- Accept -> %f",tree->c_lnL); */
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

  /* For(j,1+(int)(tree->n_otu/5)) */
  For(j,1+tree->n_otu)
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
      new_alnL = Lk(NULL,tree);
      
      ratio += (new_alnL - cur_alnL);
      ratio += (new_glnL - cur_glnL);
      ratio += hr;
            
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
          
          if(tree->mmod->safe_phyrex == YES)
            {
              new_alnL = Lk(NULL,tree); /* Remove once checked */
              
              if(Are_Equal(new_alnL,cur_alnL,1.E-3) == NO)
                {
                  PhyML_Printf("\n== new_alnL: %f cur_alnL: %f",new_alnL,cur_alnL);
                  Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
                }
                
              
              new_glnL = PHYREX_Lk(tree); /* Same */
              
              if(Are_Equal(new_glnL,cur_glnL,1.E-3) == NO)
                {
                  PhyML_Printf("\n== new_glnL: %f cur_glnL: %f",new_glnL,cur_glnL);
                  Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
                }
            }
          else
            {
              PHYREX_Update_Lindisk_List(tree);
              tree->mmod->c_lnL = cur_glnL;
              tree->c_lnL       = cur_alnL;          
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
  phydbl hr;
  phydbl cur_rad, new_rad;
  phydbl cur_mu, new_mu;

  hr      = 0.0;
  new_rad = tree->mmod->rad;
  cur_rad = tree->mmod->rad;
  new_mu  = tree->mmod->mu;
  cur_mu  = tree->mmod->mu;

  /* new_rad = cur_rad * EXP(0.5*(Uni()-.5)); */
  /* hr += LOG(new_rad/cur_rad); */

  /* new_mu = cur_mu * EXP(0.5*(Uni()-.5)); */
  /* hr += LOG(new_mu/cur_mu); */

  if(new_rad > tree->mmod->max_rad || new_rad < tree->mmod->min_rad) return; 
  if(new_mu > tree->mmod->max_mu || new_mu < tree->mmod->min_mu)     return; 

  tree->mmod->rad = new_rad;
  tree->mmod->mu  = new_mu;

  disk = tree->disk->prev;
  n_disks_cur = 0;
  do
    {
      if(disk->ldsk != NULL && disk->ldsk->n_next == 1) n_disks_cur++;
      disk = disk->prev;
    }
  while(disk && disk->prev);
  
  n_disks_new = (int)Rpois(n_disks_cur+1);
  
  hr += Dpois(n_disks_cur,n_disks_new+1,YES);
  hr -= Dpois(n_disks_new,n_disks_cur+1,YES);

  /* K = 0.2; */
  /* n_disks_new = (int)Rgamma((phydbl)n_disks_cur/K,K); */

  /* hr += LOG(Pgamma(n_disks_cur+1,(phydbl)n_disks_new/K,K) - Pgamma(n_disks_cur,(phydbl)n_disks_new/K,K)); */
  /* hr -= LOG(Pgamma(n_disks_new+1,(phydbl)n_disks_cur/K,K) - Pgamma(n_disks_new,(phydbl)n_disks_cur/K,K)); */

  if(n_disks_new < n_disks_cur) MCMC_PHYREX_Delete_Hit(hr, n_disks_cur - n_disks_new, cur_rad, cur_mu, tree);
  else                          MCMC_PHYREX_Insert_Hit(hr, n_disks_new - n_disks_cur, cur_rad, cur_mu, tree);  
}
#endif


/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
/* Insert a disk with a new lineage displacement */
#ifdef PHYREX
void MCMC_PHYREX_Insert_Hit(phydbl hr, int n_insert_disks, phydbl cur_rad, phydbl cur_mu, t_tree *tree)
{
  t_dsk  *disk,**new_disk,*young_disk;
  t_ldsk **young_ldsk, **old_ldsk, **new_ldsk;
  phydbl T,t;
  phydbl cur_glnL, new_glnL;
  phydbl u,alpha,ratio;
  int i,j,*dir_old_young,n_valid_disks,err;

  tree->mcmc->run_move[tree->mcmc->num_move_phyrex_indel_hit]++;

  disk            = NULL;
  new_glnL        = tree->mmod->c_lnL;
  cur_glnL        = tree->mmod->c_lnL;

  disk = tree->disk->prev;
  n_valid_disks = 0;
  do
    {
      if(disk->ldsk != NULL && disk->ldsk->n_next == 1) n_valid_disks++;
      disk = disk->prev;
    }
  while(disk && disk->prev);
  
  if(n_insert_disks == 0) 
    {
      tree->mmod->rad = cur_rad;
      tree->mmod->mu = cur_mu;
      return; 
    }

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

      assert(young_disk->n_ldsk_a);

      /* Which lineage is going to be hit ? */
      hr += LOG(young_disk->n_ldsk_a);
      
      young_ldsk[j] = young_disk->ldsk_a[Rand_Int(0,young_disk->n_ldsk_a-1)];
      old_ldsk[j]   = young_ldsk[j]->prev;
  
      if(old_ldsk[j]->disk->time > young_ldsk[j]->disk->time) 
        {
          PhyML_Printf("\n== young_ldsk: %f",young_ldsk[j]->disk->time);
          PhyML_Printf("\n== old_ldsk: %f",old_ldsk[j]->disk->time);
          assert(FALSE);
        }
      
      /* Direction from old to young ldsk */
      dir_old_young[j] = PHYREX_Get_Next_Direction(young_ldsk[j],old_ldsk[j]);
      assert(dir_old_young[j] != -1);

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

      /* Really necessary given that lk function is called afterwards? */
      PHYREX_Update_Lindisk_List(tree);
      
      /* Sample position of the displaced ldsk */
      For(i,tree->mmod->n_dim)
        {
          new_ldsk[j]->coord->lonlat[i] = Rnorm_Trunc(0.5*(young_ldsk[j]->coord->lonlat[i]+old_ldsk[j]->coord->lonlat[i]),
                                                      1.0*tree->mmod->rad,
                                                      0.0,tree->mmod->lim->lonlat[i],&err);
          hr -= Log_Dnorm_Trunc(new_ldsk[j]->coord->lonlat[i],
                                0.5*(young_ldsk[j]->coord->lonlat[i]+old_ldsk[j]->coord->lonlat[i]),
                                1.0*tree->mmod->rad,
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
                                1.0*tree->mmod->rad,
                                0.0,tree->mmod->lim->lonlat[i],&err);

        }

    }

  T = PHYREX_Tree_Height(tree);

  hr -= LnChoose(n_valid_disks+n_insert_disks,n_insert_disks);
  hr += n_insert_disks * LOG(-T);  
  hr -= LnFact(n_insert_disks);

  new_glnL = PHYREX_Lk(tree);
  
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

      tree->mmod->rad = cur_rad;
      tree->mmod->mu  = cur_mu;
      
      for(j=n_insert_disks-1;j>=0;j--)
        {
          PHYREX_Remove_Disk(new_disk[j]);      
          Free_Disk(new_disk[j]);
          Free_Ldisk(new_ldsk[j]);

          old_ldsk[j]->next[dir_old_young[j]] = young_ldsk[j];
          young_ldsk[j]->prev = old_ldsk[j];
        }
      

      if(tree->mmod->safe_phyrex == YES)
        {
          new_glnL = PHYREX_Lk(tree);
          if(Are_Equal(new_glnL,cur_glnL,1.E-3) == NO)
            {
              PhyML_Printf("\n== new_glnL: %f cur_glnL: %f",new_glnL,cur_glnL);
              Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
            }
        }
      else
        {
          PHYREX_Update_Lindisk_List(tree);
          tree->mmod->c_lnL = cur_glnL;
        }

    }
  else
    {
      /* printf("+ Accept"); */
      /* if(indel > 0) printf("\n. Accept"); */
      tree->mcmc->acc_move[tree->mcmc->num_move_phyrex_indel_hit]++;
    }


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
void MCMC_PHYREX_Delete_Hit(phydbl hr, int n_delete_disks, phydbl cur_rad, phydbl cur_mu, t_tree *tree)
{
  phydbl u,alpha,ratio;
  phydbl cur_glnL,new_glnL,T;
  t_dsk  *disk,**target_disk,**valid_disks;
  t_ldsk **target_ldsk,**old_ldsk,**young_ldsk;
  int i,j,block,n_valid_disks,*dir_old_young,*permut,err;

  tree->mcmc->run_move[tree->mcmc->num_move_phyrex_indel_hit]++;

  disk            = NULL;
  valid_disks     = NULL;
  new_glnL        = tree->mmod->c_lnL;
  cur_glnL        = tree->mmod->c_lnL;  
  ratio           = 0.0;
  block           = 100;

  if(tree->disk->next) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  disk = tree->disk->prev;
  if(!disk->prev) 
    {
      tree->mmod->rad = cur_rad;
      tree->mmod->mu = cur_mu;
      return;
    }

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

  if(!n_valid_disks) 
    {
      tree->mmod->rad = cur_rad;
      tree->mmod->mu = cur_mu;
      return;
    }
  
  if(n_valid_disks - n_delete_disks < 0) 
    { 
      tree->mmod->rad = cur_rad;
      tree->mmod->mu = cur_mu;
      return; 
    }

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
  
      assert(target_ldsk[j] != NULL);
      assert(target_ldsk[j]->n_next == 1);

      old_ldsk[j]   = target_ldsk[j]->prev;
      young_ldsk[j] = target_ldsk[j]->next[0];
            
      dir_old_young[j] = PHYREX_Get_Next_Direction(young_ldsk[j],old_ldsk[j]);
      assert(dir_old_young[j] != -1);
      
      /* Part of the Hastings ratio corresponding to the probability of selecting */
      /* one of target_disk->n_ldsk_a to be hit */ 
      hr -= LOG(target_disk[j]->next->n_ldsk_a);
            
      /* Density for position of the displaced ldsk */
      For(i,tree->mmod->n_dim)
        {
          hr += Log_Dnorm_Trunc(target_ldsk[j]->coord->lonlat[i],
                                0.5*(young_ldsk[j]->coord->lonlat[i]+old_ldsk[j]->coord->lonlat[i]),
                                1.0*tree->mmod->rad,
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

  ratio = EXP(ratio);
  alpha = MIN(1.,ratio);
  
  /* Always accept move */
  if(tree->mcmc->always_yes == YES && new_glnL > UNLIKELY) alpha = 1.0;

  u = Uni();
  
  /* printf("\n- Delete hit %15f %15f %5d",new_glnL - cur_glnL, alpha, n_delete_disks); */

  if(u > alpha) /* Reject */
    {
      tree->mmod->rad = cur_rad;
      tree->mmod->mu = cur_mu;
      
      /* printf("- Reject"); */
      for(j=n_delete_disks-1;j>=0;j--) 
        {
          PHYREX_Insert_Disk(target_disk[j],tree);      
          old_ldsk[j]->next[dir_old_young[j]] = target_ldsk[j];
          young_ldsk[j]->prev                 = target_ldsk[j];
        }
      
      if(tree->mmod->safe_phyrex == YES)
        {
          new_glnL = PHYREX_Lk(tree);
          if(Are_Equal(new_glnL,cur_glnL,1.E-3) == NO)
            {
              PhyML_Printf("\n== new_glnL: %f cur_glnL: %f",new_glnL,cur_glnL);
              Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
            }
        }
      else
        {
          PHYREX_Update_Lindisk_List(tree);
          tree->mmod->c_lnL = cur_glnL;
        }
    }
  else
    {
      /* printf(" ***"); */

      For(j,n_delete_disks) Free_Disk(target_disk[j]);
      For(j,n_delete_disks) Free_Ldisk(target_ldsk[j]);
      tree->mcmc->acc_move[tree->mcmc->num_move_phyrex_indel_hit]++;
    }
  

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
  t_dsk  *disk,*prune_disk,*regraft_disk,**valid_disks,*dum_dum_disk;
  t_ldsk *prune_ldsk,*regraft_ldsk,*prune_daughter_ldsk,*cur_path,*new_path,*ldsk,*ldsk_dum;
  phydbl *prob_disks;
  int i,block,n_valid_disks,prune_next_num,num_prune_disk,n_prune_disks;
  phydbl rate,dt,sizeT,sum;
  phydbl max_dist, param_exp;
  int cur_path_len;
  int n_hits,n_iter;
  int cur_pos,new_pos;

  n_iter   = MAX(1,(int)(tree->n_otu/5));

  while(n_iter--)
    {
      tree->mcmc->run_move[tree->mcmc->num_move_phyrex_spr]++;
  
      valid_disks    = NULL;
      disk           = NULL;
      prob_disks     = NULL;
      new_glnL       = tree->mmod->c_lnL;
      cur_glnL       = tree->mmod->c_lnL;
      new_alnL       = tree->c_lnL;
      cur_alnL       = tree->c_lnL;
      hr             = 0.0;
      ratio          = 0.0;
      block          = 100;
      cur_pos        = -1;
      new_pos        = -1;
      num_prune_disk = -1;
      sizeT          = PHYREX_Time_Tree_Length(tree);

      if(tree->disk->next) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
      
      /* Get a ldsk from which you can prune a lineage */
      disk = tree->disk->prev;
      n_prune_disks = 0;
      do
        {
          /* Include only disks with displacement that are coalescent events */
          if(disk->ldsk != NULL && disk->ldsk->n_next > 1)
            {
              /* Root has degree 2: don't pull any lineage */
              if(!(disk->prev == NULL && disk->ldsk->n_next == 2))
                {
                  if(!n_prune_disks) valid_disks = (t_dsk **)mCalloc(block,sizeof(t_dsk *));
                  else if(!(n_prune_disks%block)) valid_disks = (t_dsk **)mRealloc(valid_disks,n_prune_disks+block,sizeof(t_dsk *));
                  valid_disks[n_prune_disks] = disk;
                  n_prune_disks++;
                }
            }
          disk = disk->prev;
        }
      while(disk);
      
      if(!n_prune_disks) return;
            
      /* Uniform selection of a disk where a coalescent occurred */
      i = Rand_Int(0,n_prune_disks-1);
      prune_disk = valid_disks[i];
      Free(valid_disks);
      
      hr -= LOG(1./n_prune_disks);

      prune_ldsk = prune_disk->ldsk;
      /* Which daughter lineage are we pruning? */
      prune_next_num = Rand_Int(0,prune_ldsk->n_next-1);
      prune_daughter_ldsk = prune_ldsk->next[prune_next_num];
      while(prune_daughter_ldsk->n_next < 2 && /* prune_daughter is a coalescent node */
            prune_daughter_ldsk->disk->next) prune_daughter_ldsk = prune_daughter_ldsk->next[0];
      
      /* prune_daughter_ldsk has to be the next coalescent node under prune_ldsk for this move to work */
      
      /* Proba of pruning this particular one */
      hr -= LOG(1./prune_ldsk->n_next);
            
      /* Get a ldsk to reattach the pruned lineage to */
      disk = tree->disk->prev;
      n_valid_disks = 0;
      do
        {
          /* Include only disks with displacement that are not younger
             than prune_daughter_ldsk */
          if((disk->ldsk != NULL) && 
             (disk->ldsk->n_next > 0) && 
             (disk->time < prune_daughter_ldsk->disk->time) && 
             (PHYREX_Is_On_Path(disk->ldsk,prune_daughter_ldsk,prune_ldsk) == NO))
            {
              if(!n_valid_disks) 
                {
                  valid_disks = (t_dsk **)mCalloc(block,sizeof(t_dsk *));
                  prob_disks  = (phydbl *)mCalloc(block,sizeof(phydbl));
                }
              else if(!(n_valid_disks%block)) 
                {
                  valid_disks = (t_dsk **)mRealloc(valid_disks,n_valid_disks+block,sizeof(t_dsk *));
                  prob_disks  = (phydbl *)mRealloc(prob_disks,n_valid_disks+block,sizeof(phydbl));
                }

              valid_disks[n_valid_disks] = disk;

              prob_disks[n_valid_disks] = PHYREX_Dist_Between_Two_Ldsk(disk->ldsk,
                                                                       prune_daughter_ldsk,
                                                                       tree);

              if(disk == prune_disk) 
                {
                  num_prune_disk = n_valid_disks;
                }

              n_valid_disks++;
            }
          disk = disk->prev;
        }
      while(disk);
      
      if(!n_valid_disks) return;

      if(n_valid_disks > 1)
        {
          max_dist = -INFINITY;
          For(i,n_valid_disks) if(prob_disks[i] > max_dist) max_dist = prob_disks[i];
          
          sum = 0.0;
          For(i,n_valid_disks) prob_disks[i] /= max_dist;

          param_exp = 1.0;
          For(i,n_valid_disks) prob_disks[i] = Dexp(prob_disks[i],param_exp);

          sum = 0.0;
          For(i,n_valid_disks) sum += prob_disks[i];
          For(i,n_valid_disks) prob_disks[i] /= sum;
          
          i = Sample_i_With_Proba_pi(prob_disks,n_valid_disks);
          
          /* Prob of selecting this node as regraft site */
          hr -= LOG(prob_disks[i]);

          assert(num_prune_disk > -1);

          /* Prob of selecting this node as regraft site (opposite move) */
          hr += LOG(prob_disks[num_prune_disk]);
        }
      else
        {
          i = 0;
        }

      Free(prob_disks);

      /* Uniform selection of a disk among the list of valid ones */
      /* i = Rand_Int(0,n_valid_disks-1); */      

      regraft_disk = valid_disks[i];
      Free(valid_disks);
            
      regraft_ldsk = regraft_disk->ldsk;

      dum_dum_disk = (regraft_disk->time < prune_disk->time ? regraft_disk : prune_disk); 
            
      new_glnL = cur_glnL;
      new_glnL -= PHYREX_Lk_Range(prune_daughter_ldsk->disk->prev,dum_dum_disk,tree);

      /* Proba of pruning that particular lineage in the opposite move 
       */
      hr += LOG(1./(phydbl)(regraft_ldsk->n_next+1));
      
      /* Prob of selecting that prune ldsk in the opposite move. Number
         of prune sites increases by one in case actual regraft ldsk
         (i.e., regraft ldsk in the proposed move) has n_next = 1 
      */
      if(regraft_ldsk->n_next == 1) n_prune_disks += 1;
      if(prune_ldsk->n_next == 2) n_prune_disks   -= 1;
      hr += LOG(1./n_prune_disks);


      /* Prune and regraft */
      n_hits = PHYREX_Total_Number_Of_Hit_Disks(tree) - PHYREX_Total_Number_Of_Coal_Disks(tree);
      dt = FABS(prune_daughter_ldsk->disk->time - prune_ldsk->disk->time);
      cur_path_len = PHYREX_Path_Len(prune_daughter_ldsk,prune_ldsk)-2;
      rate = (phydbl)(n_hits - cur_path_len)/(sizeT - dt);
      
      hr += PHYREX_Path_Logdensity(prune_daughter_ldsk,prune_ldsk,rate*dt,1.0*tree->mmod->rad,tree);
      
      dt = FABS(prune_daughter_ldsk->disk->time - regraft_ldsk->disk->time);
      new_path = PHYREX_Generate_Path(prune_daughter_ldsk,regraft_ldsk,rate*dt,1.0*tree->mmod->rad,tree);
      cur_path = PHYREX_Remove_Path(prune_daughter_ldsk,prune_ldsk,&cur_pos,tree);
      PHYREX_Insert_Path(prune_daughter_ldsk,regraft_ldsk,new_path,regraft_ldsk->n_next,tree);
      
      hr -= PHYREX_Path_Logdensity(prune_daughter_ldsk,regraft_ldsk,rate*dt,1.0*tree->mmod->rad,tree);
      
      /* prune_daughter_ldsk->prev = regraft_ldsk; */
      /* PHYREX_Remove_Lindisk_Next(prune_ldsk,prune_daughter_ldsk); */
      /* PHYREX_Random_Insert_Ldsk_In_Next_List(prune_daughter_ldsk,regraft_ldsk); */
      
      PHYREX_Ldsk_To_Tree(tree);
      Update_Ancestors(tree->n_root,tree->n_root->v[2],tree);
      Update_Ancestors(tree->n_root,tree->n_root->v[1],tree);
      RATES_Fill_Lca_Table(tree);
      RATES_Update_Cur_Bl(tree);


      /* new_glnL = PHYREX_Lk(tree); */

      new_glnL += PHYREX_Lk_Range(prune_daughter_ldsk->disk->prev,dum_dum_disk,tree);
      tree->mmod->c_lnL = new_glnL;

      new_alnL = Lk(NULL,tree);
      
      ratio += (new_alnL - cur_alnL);
      ratio += (new_glnL - cur_glnL);
      ratio += hr;
            
      ratio = EXP(ratio);
      alpha = MIN(1.,ratio);
      
      /* Always accept move */
      if(tree->mcmc->always_yes == YES && new_glnL > UNLIKELY) alpha = 1.0;
      
      /* printf("\n. %3d %3d alpha:%12f %3d",prune_ldsk->n_next+1,regraft_ldsk->n_next-1,alpha,new_n_coal); */
      
      u = Uni();
      
      if(u > alpha) /* Reject */
        {
          new_path = PHYREX_Remove_Path(prune_daughter_ldsk,regraft_ldsk,&new_pos,tree);
          
          ldsk = new_path;
          while(ldsk)
            {
              Free_Disk(ldsk->disk);
              ldsk_dum = ldsk;
              ldsk = ldsk->prev;
              Free_Ldisk(ldsk_dum);
            }
          
          PHYREX_Insert_Path(prune_daughter_ldsk,prune_ldsk,cur_path,cur_pos,tree);
          
          /* prune_daughter_ldsk->prev = prune_ldsk; */
          /* PHYREX_Remove_Lindisk_Next(regraft_ldsk,prune_daughter_ldsk); */
          /* PHYREX_Insert_Ldsk_In_Next_List(prune_daughter_ldsk,prune_next_num,prune_ldsk); */
          
          PHYREX_Ldsk_To_Tree(tree);
          Update_Ancestors(tree->n_root,tree->n_root->v[2],tree);
          Update_Ancestors(tree->n_root,tree->n_root->v[1],tree);
          RATES_Fill_Lca_Table(tree);
          RATES_Update_Cur_Bl(tree);
          

          if(tree->mmod->safe_phyrex == YES)
            {
              new_alnL = Lk(NULL,tree);
              if(Are_Equal(new_alnL,cur_alnL,1.E-3) == NO)
                {
                  PhyML_Printf("\n== new_alnL: %f cur_alnL: %f",new_alnL,cur_alnL);
                  Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
                }
                
              
              new_glnL = PHYREX_Lk(tree);
              if(Are_Equal(new_glnL,cur_glnL,1.E-3) == NO)
                {
                  PhyML_Printf("\n== new_glnL: %f cur_glnL: %f",new_glnL,cur_glnL);
                  Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
                }
            }
          else
            {
              PHYREX_Update_Lindisk_List(tree);
              tree->c_lnL       = cur_alnL;
              tree->mmod->c_lnL = cur_glnL;          
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
    }
}
#endif
  
/*////////////////////////////////////////////////////////////l
////////////////////////////////////////////////////////////*/

#ifdef PHYREX
void MCMC_PHYREX_Simulate_Backward(t_tree *tree)
{
  phydbl u,alpha,ratio,hr,T,t;
  phydbl cur_alnL, new_alnL;
  phydbl cur_glnL, new_glnL;
  t_dsk  *disk,*bkp_disk,*target_disk;
  t_ldsk **bkp_ldsk;
  int i;

  tree->mcmc->run_move[tree->mcmc->num_move_phyrex_sim]++;

  disk        = NULL;
  bkp_disk    = NULL;
  target_disk = NULL;
  bkp_ldsk    = NULL;

  new_alnL      = tree->c_lnL;
  cur_alnL      = tree->c_lnL;
  cur_glnL      = tree->mmod->c_lnL;
  new_glnL      = tree->mmod->c_lnL;
  hr            = 0.0;
  ratio         = 0.0;
  T             = 0.0;
  t             = 0.0;

  /* Chop off the tree at time t and simulate upwards from here */
  T = PHYREX_Tree_Height(tree);
  t = Uni()*T;
  
  disk = tree->disk->prev;
  while(disk && disk->time > t) disk = disk->prev;
  target_disk = disk->next;

  hr -= LOG(FABS((target_disk->prev->time - target_disk->time)/T));

  bkp_disk = target_disk->prev;

  bkp_ldsk = (t_ldsk **)mCalloc(target_disk->n_ldsk_a,sizeof(t_ldsk *));
  For(i,target_disk->n_ldsk_a) bkp_ldsk[i] = target_disk->ldsk_a[i]->prev;

  PHYREX_Simulate_Backward_Core(NO,target_disk,tree);
  
  T =  PHYREX_Tree_Height(tree);
  hr += LOG(FABS((target_disk->prev->time - target_disk->time)/T));

  new_alnL = Lk(NULL,tree);
  ratio += (new_alnL - cur_alnL);
  ratio += hr;
  
  ratio = EXP(ratio);
  alpha = MIN(1.,ratio);
    
 /* Always accept move */
  if(tree->mcmc->always_yes == YES) alpha = 1.0;
    
  u = Uni();
        
  if(u > alpha) /* Reject */
    {
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
      
      if(tree->mmod->safe_phyrex == YES)
        {
          new_alnL = Lk(NULL,tree);
          if(Are_Equal(new_alnL,cur_alnL,1.E-3) == NO)
            {
              PhyML_Printf("\n== new_alnL: %f cur_alnL: %f",new_alnL,cur_alnL);
              Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
            }
        
          
          new_glnL = PHYREX_Lk(tree);
          
          if(Are_Equal(new_glnL,cur_glnL,1.E-3) == NO)
            {
              PhyML_Printf("\n. new_glnL: %f cur_glnL: %f",new_glnL,cur_glnL);
              Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
            }
        }
      else
        {
          PHYREX_Lk(tree);
          tree->c_lnL       = cur_alnL;
          tree->mmod->c_lnL = cur_glnL;
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
  

  Free(bkp_ldsk);
}
#endif

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/


#ifdef PHYREX
void MCMC_PHYREX_Simulate_Backward_Plus(t_tree *tree)
{
  phydbl u,alpha,ratio,hr,T,t;
  phydbl cur_alnL, new_alnL;
  phydbl cur_glnL_do, new_glnL_do;
  phydbl cur_glnL, new_glnL;
  phydbl cur_lbda, new_lbda;
  phydbl cur_rad, new_rad;
  phydbl cur_mu, new_mu;
  t_dsk  *disk,*bkp_disk,*target_disk;
  t_ldsk **bkp_ldsk;
  int i;

  tree->mcmc->run_move[tree->mcmc->num_move_phyrex_sim_plus]++;

  disk          = NULL;
  bkp_disk      = NULL;
  target_disk   = NULL;
  bkp_ldsk      = NULL;
  new_alnL      = tree->c_lnL;
  cur_alnL      = tree->c_lnL;
  cur_glnL      = tree->mmod->c_lnL;
  new_glnL      = tree->mmod->c_lnL;
  cur_glnL_do   = UNLIKELY;
  new_glnL_do   = UNLIKELY;
  hr            = 0.0;
  ratio         = 0.0;
  T             = 0.0;
  t             = 0.0;
  cur_lbda      = tree->mmod->lbda;
  cur_rad       = tree->mmod->rad;
  cur_mu        = tree->mmod->mu;
  new_lbda      = cur_lbda;
  new_rad       = cur_rad;
  new_mu        = cur_mu;

  
  new_lbda = cur_lbda * EXP(0.2*(Uni()-.5));
  hr += LOG(new_lbda/cur_lbda);

  new_mu = cur_mu * EXP(0.2*(Uni()-.5));
  hr += LOG(new_mu/cur_mu);

  new_rad = cur_rad * EXP(0.2*(Uni()-.5));
  hr += LOG(new_rad/cur_rad);

  if(new_lbda < tree->mmod->min_lbda) return;
  if(new_lbda > tree->mmod->max_lbda) return;
  if(new_mu   < tree->mmod->min_mu)   return;
  if(new_mu   > tree->mmod->max_mu)   return;
  if(new_rad  < tree->mmod->min_rad)  return;
  if(new_rad  > tree->mmod->max_rad)  return;


  hr += PHYREX_LnPrior_Lbda(tree);
  hr += PHYREX_LnPrior_Mu(tree);
  hr += PHYREX_LnPrior_Radius(tree);

  /* Chop off the tree at time t and simulate upwards from here */
  T = PHYREX_Tree_Height(tree);
  t = Uni()*T;
  
  disk = tree->disk->prev;
  while(disk && disk->time > t) disk = disk->prev;
  target_disk = disk->next;

  hr -= LOG(FABS((target_disk->prev->time - target_disk->time)/T));

  bkp_disk = target_disk->prev;

  bkp_ldsk = (t_ldsk **)mCalloc(target_disk->n_ldsk_a,sizeof(t_ldsk *));
  For(i,target_disk->n_ldsk_a) bkp_ldsk[i] = target_disk->ldsk_a[i]->prev;

  cur_glnL_do = PHYREX_Lk_Range(tree->disk->prev,target_disk,tree);

  tree->mmod->lbda = new_lbda;
  tree->mmod->mu   = new_mu;
  tree->mmod->rad  = new_rad;

  new_glnL_do = PHYREX_Lk_Range(tree->disk->prev,target_disk,tree);

  hr -= PHYREX_LnPrior_Lbda(tree);
  hr -= PHYREX_LnPrior_Mu(tree);
  hr -= PHYREX_LnPrior_Radius(tree);

  PHYREX_Simulate_Backward_Core(NO,target_disk,tree);
  
  T =  PHYREX_Tree_Height(tree);
  hr += LOG(FABS((target_disk->prev->time - target_disk->time)/T));

  new_alnL = Lk(NULL,tree);
  ratio += (new_alnL - cur_alnL);
  ratio += (new_glnL_do - cur_glnL_do);
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
      
      if(tree->mmod->safe_phyrex == YES)
        {
          new_alnL = Lk(NULL,tree);
          if(Are_Equal(new_alnL,cur_alnL,1.E-3) == NO)
            {
              PhyML_Printf("\n== new_alnL: %f cur_alnL: %f",new_alnL,cur_alnL);
              Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
            }
        
          
          new_glnL = PHYREX_Lk(tree);
          
          if(Are_Equal(new_glnL,cur_glnL,1.E-3) == NO)
            {
              PhyML_Printf("\n. new_glnL: %f cur_glnL: %f",new_glnL,cur_glnL);
              Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
            }
        }
      else
        {
          PHYREX_Lk(tree);
          tree->c_lnL = cur_alnL;
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

      tree->mcmc->acc_move[tree->mcmc->num_move_phyrex_sim_plus]++;
    }
  

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
  int j,block,n_valid_disks;
  phydbl rate,dt,sizeT;
  int n_hits,n_iter,cur_path_len;
  int pos,*permut;

  n_iter      = 0;
  valid_disks = NULL;
  disk        = NULL;
  new_glnL    = tree->mmod->c_lnL;
  cur_glnL    = tree->mmod->c_lnL;
  hr          = 0.0;
  ratio       = 0.0;
  block       = 100;
  pos         = -1;
  
      
  /* Get a ldsk from which you can prune a lineage */
  disk = tree->disk->prev;
  n_valid_disks = 0;
  do
    {
      if(disk->ldsk && disk->ldsk->prev && disk->ldsk->n_next > 1)
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
      
  permut = Permutate(n_valid_disks);
  n_iter = Rand_Int(1,1+(int)(n_valid_disks/5));
  
  For(j,n_iter)
    {
      tree->mcmc->run_move[tree->mcmc->num_move_phyrex_traj]++;

      new_glnL    = tree->mmod->c_lnL;
      cur_glnL    = tree->mmod->c_lnL;
      hr          = 0.0;
      ratio       = 0.0;

      start_ldsk = valid_disks[permut[j]]->ldsk;
      assert(start_ldsk != NULL);

      end_ldsk = start_ldsk->prev;
      while(end_ldsk->n_next < 2) end_ldsk = end_ldsk->prev;

      assert(end_ldsk != NULL); 

      new_glnL = cur_glnL;
      new_glnL -= PHYREX_Lk_Range(start_ldsk->disk->prev,end_ldsk->disk,tree);
      
      n_hits       = PHYREX_Total_Number_Of_Hit_Disks(tree) - PHYREX_Total_Number_Of_Coal_Disks(tree);
      sizeT        = PHYREX_Time_Tree_Length(tree);
      dt           = FABS(start_ldsk->disk->time - end_ldsk->disk->time);
      cur_path_len = PHYREX_Path_Len(start_ldsk,end_ldsk)-2;
      rate         = (phydbl)(n_hits - cur_path_len)/(sizeT - dt);
      
      hr += PHYREX_Path_Logdensity(start_ldsk,end_ldsk,rate*dt,1.*tree->mmod->rad,tree);
      
      new_path = PHYREX_Generate_Path(start_ldsk,end_ldsk,rate*dt,1.*tree->mmod->rad,tree);
      cur_path = PHYREX_Remove_Path(start_ldsk,end_ldsk,&pos,tree);
      PHYREX_Insert_Path(start_ldsk,end_ldsk,new_path,pos,tree);
      
      hr -= PHYREX_Path_Logdensity(start_ldsk,end_ldsk,rate*dt,1.*tree->mmod->rad,tree);
      
      /* new_glnL = PHYREX_Lk(tree); */

      new_glnL += PHYREX_Lk_Range(start_ldsk->disk->prev,end_ldsk->disk,tree);
      tree->mmod->c_lnL = new_glnL;

      ratio += (new_glnL - cur_glnL);
      ratio += hr;
            
      ratio = EXP(ratio);
      alpha = MIN(1.,ratio);
      
      if(tree->mcmc->always_yes == YES && new_glnL > UNLIKELY) alpha = 1.0;
            
      /* PhyML_Printf("\n= Traj. new_alnL:%12f cur_alnL:%12f hr:%12f alpha:%12f",new_alnL,cur_alnL,hr,alpha); */

      u = Uni();
      
      if(u > alpha) /* Reject */
        {
          new_path = PHYREX_Remove_Path(start_ldsk,end_ldsk,&pos,tree);
          
          ldsk = new_path;
          while(ldsk) 
            {
              Free_Disk(ldsk->disk);
              ldsk_dum = ldsk;
              ldsk = ldsk->prev;
              Free_Ldisk(ldsk_dum);
            }
          
          PHYREX_Insert_Path(start_ldsk,end_ldsk,cur_path,pos,tree);
          
          if(tree->mmod->safe_phyrex == YES)
            {
              new_glnL = PHYREX_Lk(tree);
              
              if(Are_Equal(new_glnL,cur_glnL,1.E-3) == NO)
                {
                  PhyML_Printf("\n== new_glnL: %f cur_glnL: %f",new_glnL,cur_glnL);
                  Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
                }
            }
          else
            {
              PHYREX_Update_Lindisk_List(tree);
              tree->mmod->c_lnL = cur_glnL;
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
    }
  
  Free(valid_disks);
  Free(permut);
}
#endif

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

#ifdef PHYREX
void MCMC_PHYREX_Lbda_Times(t_tree *tree)
{
  phydbl u,alpha,ratio,hr;
  phydbl cur_glnL, new_glnL;
  phydbl cur_alnL, new_alnL;
  phydbl cur_lbda, new_lbda;
  phydbl *cur_times,dt;
  t_dsk  *disk;
  int n_inter,i;

  tree->mcmc->run_move[tree->mcmc->num_move_phyrex_lbda_times]++;

  dt          = 0.0;
  hr          = 0.0;
  ratio       = 0.0;
  new_glnL    = tree->mmod->c_lnL;
  cur_glnL    = tree->mmod->c_lnL;
  new_alnL    = tree->c_lnL;
  cur_alnL    = tree->c_lnL;
  cur_lbda    = tree->mmod->lbda;
  new_lbda    = cur_lbda;

  n_inter = PHYREX_Total_Number_Of_Intervals(tree);

  cur_times = (phydbl *)mCalloc(n_inter,sizeof(phydbl));

  disk = tree->disk;
  while(disk->prev != NULL) disk = disk->prev;
  hr += (n_inter)*LOG(tree->mmod->lbda) + tree->mmod->lbda*disk->time;

  new_glnL  = cur_glnL;
  new_glnL -= (n_inter)*LOG(tree->mmod->lbda) + tree->mmod->lbda*disk->time;
  new_lbda = cur_lbda * EXP(0.3*(Uni()-.5));
  hr += LOG(new_lbda/cur_lbda);

  tree->mmod->lbda = new_lbda;

  i = 0;
  disk = tree->disk->prev;  
  do
    {
      cur_times[i++] = disk->time;
      dt = Rexp(tree->mmod->lbda);
      disk->time = disk->next->time - dt;
      disk = disk->prev;
    }
  while(disk != NULL);

  disk = tree->disk;
  while(disk->prev != NULL) disk = disk->prev;
  hr -= n_inter*LOG(tree->mmod->lbda) + tree->mmod->lbda*disk->time;
  
  new_glnL += n_inter*LOG(tree->mmod->lbda) + tree->mmod->lbda*disk->time;;
  tree->mmod->c_lnL = new_glnL;

  new_alnL = Lk(NULL,tree);
    
  ratio += (new_glnL - cur_glnL);
  ratio += (new_alnL - cur_alnL);
  ratio += hr;
            
  ratio = EXP(ratio);
  alpha = MIN(1.,ratio);

  u = Uni();

  if(u > alpha) /* Reject */
    {
      tree->mmod->lbda = cur_lbda;

      i = 0;
      disk = tree->disk->prev;
      do
        {
          disk->time = cur_times[i++];
          disk = disk->prev;
        }
      while(disk != NULL);
      
      if(tree->mmod->safe_phyrex == YES)
        {
          new_alnL = Lk(NULL,tree);
          
          if(Are_Equal(new_alnL,cur_alnL,1.E-3) == NO)
            {
              PhyML_Printf("\n== new_alnL: %f cur_alnL: %f",new_alnL,cur_alnL);
              Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
            }
            
          
          new_glnL = PHYREX_Lk(tree);
          
          if(Are_Equal(new_glnL,cur_glnL,1.E-3) == NO)
            {
              PhyML_Printf("\n== new_glnL: %f cur_glnL: %f",new_glnL,cur_glnL);
              Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
            }
        }
      else
        {
          tree->c_lnL       = cur_alnL;
          tree->mmod->c_lnL = cur_glnL;
        }
    }
  else
    {
      tree->mcmc->acc_move[tree->mcmc->num_move_phyrex_lbda_times]++;
    }

  Free(cur_times);
}
#endif

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

#ifdef PHYREX
void MCMC_PHYREX_Disk_Multi(t_tree *tree)
{
  phydbl u,alpha,ratio;
  phydbl cur_glnL, new_glnL, hr;
  t_dsk  *disk,**target_disk,**all_disks;
  int i,j,block,n_all_disks,n_move_disks,*permut;
  int err;

  disk          = NULL;
  new_glnL      = tree->mmod->c_lnL;
  cur_glnL      = tree->mmod->c_lnL;
  hr            = 0.0;
  ratio         = 0.0;
  block         = 100;
  all_disks     = NULL;

  tree->mcmc->run_move[tree->mcmc->num_move_phyrex_disk_multi]++;

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
  
  n_move_disks = Rand_Int(1,1+(int)(n_all_disks/5));
  /* n_move_disks = n_all_disks; */

  permut = Permutate(n_all_disks);

  For(i,n_move_disks)
    {
      target_disk[i] = all_disks[permut[i]];
  
      assert(target_disk[i]);

      PHYREX_Store_Geo_Coord(target_disk[i]->centr);
      
      if(target_disk[i]->ldsk != NULL)
        {
          For(j,tree->mmod->n_dim)
            {
              target_disk[i]->centr->lonlat[j] = 
                Rnorm_Trunc(target_disk[i]->centr->lonlat[j],
                            1.*tree->mmod->rad,
                            0.0,
                            tree->mmod->lim->lonlat[j],&err);
              
              hr -= Log_Dnorm_Trunc(target_disk[i]->centr->lonlat[j],
                                    target_disk[i]->centr->cpy->lonlat[j],
                                    1.*tree->mmod->rad,                                                         
                                    0.0,
                                    tree->mmod->lim->lonlat[j],&err);
              
              hr += Log_Dnorm_Trunc(target_disk[i]->centr->cpy->lonlat[j],
                                    target_disk[i]->centr->lonlat[j],
                                    1.*tree->mmod->rad,
                                    0.0,
                                    tree->mmod->lim->lonlat[j],&err);
            }
        }
      else
        {
          For(j,tree->mmod->n_dim)
            {
              target_disk[i]->centr->lonlat[j] = 
                Rnorm_Trunc(target_disk[i]->centr->lonlat[j],
                            1.*tree->mmod->rad,
                            0.0,
                            tree->mmod->lim->lonlat[j],&err);
              
              hr -= Log_Dnorm_Trunc(target_disk[i]->centr->lonlat[j],
                                    target_disk[i]->centr->cpy->lonlat[j],
                                    1.*tree->mmod->rad,                                                         
                                    0.0,
                                    tree->mmod->lim->lonlat[j],&err);
              
              hr += Log_Dnorm_Trunc(target_disk[i]->centr->cpy->lonlat[j],
                                    target_disk[i]->centr->lonlat[j],
                                    1.*tree->mmod->rad,
                                    0.0,
                                    tree->mmod->lim->lonlat[j],&err);
            }
        }
    }
  
  Free(permut);

  new_glnL = PHYREX_Lk(tree);
  tree->mmod->c_lnL = new_glnL;

  ratio += (new_glnL - cur_glnL);
  ratio += hr;
  
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
      
      if(tree->mmod->safe_phyrex == YES)
        {
          new_glnL = PHYREX_Lk(tree);
          if(Are_Equal(new_glnL,cur_glnL,1.E-3) == NO)
            {
              PhyML_Printf("\n. new_glnL: %f cur_glnL: %f",new_glnL,cur_glnL);
              Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
            }
        }
      else
        {
          tree->mmod->c_lnL = cur_glnL;
        }
    }
  else
    {
      tree->mcmc->acc_move[tree->mcmc->num_move_phyrex_disk_multi]++;
    }


  Free(all_disks);
  Free(target_disk);
}
#endif

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
/* Move ldsk on landscape */ 
#ifdef PHYREX
void MCMC_PHYREX_Ldsk_Multi(t_tree *tree)
{
  phydbl u,alpha,ratio;
  phydbl cur_glnL, new_glnL, hr, c, o, f;
  t_dsk  *disk,**target_disk,**all_disks;
  int i,j,block,n_all_disks,n_move_ldsk,*permut;
  int err;

  disk          = NULL;
  new_glnL      = tree->mmod->c_lnL;
  cur_glnL      = tree->mmod->c_lnL;
  hr            = 0.0;
  ratio         = 0.0;
  block         = 100;
  all_disks     = NULL;
  o             = -1.;
  c             = -1.;
  f             = -1.;

  tree->mcmc->run_move[tree->mcmc->num_move_phyrex_ldsk_multi]++;

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
  
  n_move_ldsk = Rand_Int(1,1+(int)(n_all_disks/5));
  /* n_move_ldsk = n_all_disks; */
  
  target_disk = (t_dsk **)mCalloc(n_all_disks,sizeof(t_dsk *));

  permut = Permutate(n_all_disks);

  For(i,n_move_ldsk)
    {
      target_disk[i] = all_disks[permut[i]];
      
      PHYREX_Store_Geo_Coord(target_disk[i]->ldsk->coord);

      For(j,tree->mmod->n_dim)
        {
          /* c: center; o: pos of direct ldsk ancestor */
          c = target_disk[i]->centr->lonlat[j];
          o = target_disk[i]->ldsk->prev ? target_disk[i]->ldsk->prev->coord->lonlat[j] : target_disk[i]->centr->lonlat[j];
          f = target_disk[i]->ldsk->prev ? (2./3.) : (1.0);
          
          target_disk[i]->ldsk->coord->lonlat[j] =
            Rnorm_Trunc((o+2.*c)/3.,
                        1.*SQRT(f*tree->mmod->rad*tree->mmod->rad),
                        0.0,
                        tree->mmod->lim->lonlat[j],&err);
          
          
          For(j,tree->mmod->n_dim) hr -= Log_Dnorm_Trunc(target_disk[i]->ldsk->coord->lonlat[j],
                                                         (o+2.*c)/3.,
                                                         1.*SQRT(f*tree->mmod->rad*tree->mmod->rad),
                                                         0.0,
                                                         tree->mmod->lim->lonlat[j],&err);
          
          For(j,tree->mmod->n_dim) hr += Log_Dnorm_Trunc(target_disk[i]->ldsk->coord->cpy->lonlat[j],
                                                         (o+2.*c)/3.,
                                                         1.*SQRT(f*tree->mmod->rad*tree->mmod->rad),
                                                         0.0,
                                                         tree->mmod->lim->lonlat[j],&err);
        }
    }

  Free(permut);

  new_glnL = PHYREX_Lk(tree);

  ratio += (new_glnL - cur_glnL);
  ratio += hr;
  
  ratio = EXP(ratio);
  alpha = MIN(1.,ratio);
  
  /* Always accept move */
  if(tree->mcmc->always_yes == YES && new_glnL > UNLIKELY) alpha = 1.0;

  u = Uni();
  
  /* printf("\n- Move_ldsk %15f",new_glnL-cur_glnL); */

  if(u > alpha) /* Reject */
    {
      /* printf("- Reject"); */
      
      For(i,n_move_ldsk) PHYREX_Restore_Geo_Coord(target_disk[i]->ldsk->coord);

      if(tree->mmod->safe_phyrex == YES)
        {
          new_glnL = PHYREX_Lk(tree);
          if(Are_Equal(new_glnL,cur_glnL,1.E-3) == NO)
            {
              PhyML_Printf("\n. new_glnL: %f cur_glnL: %f",new_glnL,cur_glnL);
              Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
            }
        }
      else
        {
          tree->mmod->c_lnL = cur_glnL;
        }
    }
  else
    {
      tree->mcmc->acc_move[tree->mcmc->num_move_phyrex_ldsk_multi]++;
    }

   
  Free(all_disks);
  Free(target_disk);
}
#endif

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

#ifdef PHYREX
void MCMC_PHYREX_Ldsk_And_Disk(t_tree *tree)
{
  phydbl u,alpha,ratio;
  phydbl cur_glnL, new_glnL, hr;
  t_dsk  *disk,**target_disk,**all_disks;
  int i,j,block,n_all_disks,n_move_ldsk,*permut;
  int err;
  phydbl cur_rad, new_rad;

  disk          = NULL;
  new_glnL      = tree->mmod->c_lnL;
  cur_glnL      = tree->mmod->c_lnL;
  hr            = 0.0;
  ratio         = 0.0;
  block         = 100;
  all_disks     = NULL;
  cur_rad       = tree->mmod->rad;
  new_rad       = tree->mmod->rad;

  tree->mcmc->run_move[tree->mcmc->num_move_phyrex_ldsk_and_disk]++;

  new_rad = cur_rad * EXP(0.2*(Uni()-.5));
  hr += LOG(new_rad/cur_rad);

  if(new_rad > tree->mmod->max_rad) return;
  if(new_rad < tree->mmod->min_rad) return;

  tree->mmod->rad = new_rad;

  if(tree->disk->next) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  disk = tree->disk->prev;
  n_all_disks = 0;
  do
    {
      if(disk->ldsk != NULL && disk->ldsk->n_next >= 1)
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
  
  n_move_ldsk = Rand_Int(1,1+(int)(n_all_disks/5));
  /* n_move_ldsk = n_all_disks; */
  
  target_disk = (t_dsk **)mCalloc(n_all_disks,sizeof(t_dsk *));

  permut = Permutate(n_all_disks);

  For(i,n_move_ldsk)
    {
      target_disk[i] = all_disks[permut[i]];
      
      PHYREX_Store_Geo_Coord(target_disk[i]->ldsk->coord);
            
      For(j,tree->mmod->n_dim)
        target_disk[i]->ldsk->coord->lonlat[j] =
        Rnorm_Trunc(target_disk[i]->centr->lonlat[j],
                    1.*new_rad,
                    0.0,
                    tree->mmod->lim->lonlat[j],&err);
      

      For(j,tree->mmod->n_dim) hr -= Log_Dnorm_Trunc(target_disk[i]->ldsk->coord->lonlat[j],
                                                     target_disk[i]->centr->lonlat[j],
                                                     1.*new_rad,
                                                     0.0,
                                                     tree->mmod->lim->lonlat[j],&err);
      
      For(j,tree->mmod->n_dim) hr += Log_Dnorm_Trunc(target_disk[i]->ldsk->coord->cpy->lonlat[j],
                                                     target_disk[i]->centr->lonlat[j],
                                                     1.*cur_rad,
                                                     0.0,
                                                     tree->mmod->lim->lonlat[j],&err);
    }


  Free(permut);

  new_glnL = PHYREX_Lk(tree);

  ratio += (new_glnL - cur_glnL);
  ratio += hr;
  
  ratio = EXP(ratio);
  alpha = MIN(1.,ratio);
  
  /* Always accept move */
  if(tree->mcmc->always_yes == YES && new_glnL > UNLIKELY) alpha = 1.0;

  u = Uni();
  
  /* printf("\n- Move_ldsk %15f",new_glnL-cur_glnL); */

  if(u > alpha) /* Reject */
    {
      /* printf("- Reject"); */
      
      tree->mmod->rad = cur_rad;

      For(i,n_move_ldsk) PHYREX_Restore_Geo_Coord(target_disk[i]->ldsk->coord);
      /* For(i,n_move_ldsk) PHYREX_Restore_Geo_Coord(target_disk[i]->centr); */

      if(tree->mmod->safe_phyrex == YES)
        {
          new_glnL = PHYREX_Lk(tree);
          if(Are_Equal(new_glnL,cur_glnL,1.E-3) == NO)
            {
              PhyML_Printf("\n. new_glnL: %f cur_glnL: %f",new_glnL,cur_glnL);
              Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
            }
        }
      else
        {
          tree->mmod->c_lnL = cur_glnL;
        }
    }
  else
    {
      tree->mcmc->acc_move[tree->mcmc->num_move_phyrex_ldsk_and_disk]++;
    }
   
  Free(all_disks);
  Free(target_disk);
}
#endif

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

#ifdef PHYREX
void MCMC_PHYREX_Ldsk_Given_Disk(t_tree *tree)
{
  phydbl u,alpha,ratio,hr,o,c,f;
  phydbl cur_glnL, new_glnL;
  t_dsk  *disk,**all_disks;
  int i,j,err,n_all_disks,block,n_move_ldsk,*permut;

  block       = 100;
  all_disks   = NULL;
  n_all_disks = 0;
  o           = -1.;
  c           = -1.;
  f           = -1.;

  disk = tree->disk->prev;
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
  
  /* n_move_ldsk = Rand_Int(1,1+(int)(n_all_disks/5)); */
  n_move_ldsk = n_all_disks;
  
  permut = Permutate(n_all_disks);

  For(i,n_move_ldsk)
    {
      tree->mcmc->run_move[tree->mcmc->num_move_phyrex_ldsk_given_disk]++;
      
      disk = all_disks[permut[i]];
            
      hr          = 0.0;
      ratio       = 0.0;
      cur_glnL    = tree->mmod->c_lnL;
      
      new_glnL = cur_glnL;
      new_glnL -= PHYREX_Lk_Range(disk,disk->ldsk->prev ? disk->ldsk->prev->disk : NULL,tree);
      
      PHYREX_Store_Geo_Coord(disk->ldsk->coord);
                  
      For(j,tree->mmod->n_dim)
        {
          /* c: center; o: pos of direct ldsk ancestor */
          c = disk->centr->lonlat[j];
          o = disk->ldsk->prev ? disk->ldsk->prev->coord->lonlat[j] : disk->centr->lonlat[j];
          f = disk->ldsk->prev ? (2./3.) : (1.0);

          disk->ldsk->coord->lonlat[j] =
            Rnorm_Trunc((o+2.*c)/3.,
                        1.*SQRT(f*tree->mmod->rad*tree->mmod->rad),
                        0.0,
                        tree->mmod->lim->lonlat[j],&err);
        

          hr -= Log_Dnorm_Trunc(disk->ldsk->coord->lonlat[j],
                                (o+2.*c)/3.,
                                1.*SQRT(f*tree->mmod->rad*tree->mmod->rad),
                                0.0,
                                tree->mmod->lim->lonlat[j],&err);
          
          hr += Log_Dnorm_Trunc(disk->ldsk->coord->cpy->lonlat[j],
                                (o+2.*c)/3.,
                                1.*SQRT(f*tree->mmod->rad*tree->mmod->rad),
                                0.0,
                                tree->mmod->lim->lonlat[j],&err);
        }
      
      new_glnL += PHYREX_Lk_Range(disk,disk->ldsk->prev ? disk->ldsk->prev->disk : NULL,tree);
      tree->mmod->c_lnL = new_glnL;
      
      ratio += (new_glnL - cur_glnL);
      ratio += hr;
      
      ratio = EXP(ratio);
      alpha = MIN(1.,ratio);
      
      u = Uni();
      
      if(u > alpha) /* Reject */
        {
          PHYREX_Restore_Geo_Coord(disk->ldsk->coord);
          
          if(tree->mmod->safe_phyrex == YES)
            {
              new_glnL = PHYREX_Lk(tree);
              
              if(Are_Equal(new_glnL,cur_glnL,1.E-3) == NO)
                {
                  PhyML_Printf("\n. disk->ldsk->prev: %p",disk->ldsk->prev);
                  PhyML_Printf("\n== new_glnL: %f cur_glnL: %f",new_glnL,cur_glnL);
                  Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
                }
            }
          else
            {
              tree->mmod->c_lnL = cur_glnL;
            }
        }
      else
        {
          tree->mcmc->acc_move[tree->mcmc->num_move_phyrex_ldsk_given_disk]++;
        }
    }

  Free(permut);
  Free(all_disks);
}
#endif

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

#ifdef PHYREX
void MCMC_PHYREX_Disk_Given_Ldsk(t_tree *tree)
{
  phydbl u,alpha,ratio,hr;
  phydbl cur_glnL, new_glnL;
  t_dsk  *disk,**all_disks;
  int i,j,n_all_disks,block,n_move_ldsk,*permut;

  block       = 100;
  all_disks   = NULL;
  n_all_disks = 0;

  disk = tree->disk->prev;
  do
    {
      if(TRUE)
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
  
  /* n_move_ldsk = Rand_Int(1,1+(int)(n_all_disks/10)); */
  n_move_ldsk = n_all_disks;
  
  permut = Permutate(n_all_disks);

  For(i,n_move_ldsk)
    {
      tree->mcmc->run_move[tree->mcmc->num_move_phyrex_disk_given_ldsk]++;
      
      disk = all_disks[permut[i]];
      
      hr          = 0.0;
      ratio       = 0.0;
      cur_glnL    = tree->mmod->c_lnL;
      
      new_glnL = cur_glnL;

      if(disk->ldsk != NULL)
        new_glnL -= PHYREX_Lk_Range(disk,disk->ldsk->prev ? disk->ldsk->prev->disk : NULL,tree);
      else
        new_glnL -= PHYREX_Lk_Range(disk,disk->prev,tree);

      PHYREX_Store_Geo_Coord(disk->centr);

      For(j,tree->mmod->n_dim) disk->centr->lonlat[j] = Uni()*tree->mmod->lim->lonlat[j];

      if(disk->ldsk != NULL)
        new_glnL += PHYREX_Lk_Range(disk,disk->ldsk->prev ? disk->ldsk->prev->disk : NULL,tree);
      else
        new_glnL += PHYREX_Lk_Range(disk,disk->prev,tree);
      
      tree->mmod->c_lnL = new_glnL;
      
      ratio += (new_glnL - cur_glnL);
      ratio += hr;
      
      ratio = EXP(ratio);
      alpha = MIN(1.,ratio);
      
      u = Uni();
      
      if(u > alpha) /* Reject */
        {
          PHYREX_Restore_Geo_Coord(disk->centr);
          
          if(tree->mmod->safe_phyrex == YES)
            {
              new_glnL = PHYREX_Lk(tree);
              
              if(Are_Equal(new_glnL,cur_glnL,1.E-3) == NO)
                {
                  PhyML_Printf("\n== disk->ldsk->prev: %p",disk->ldsk->prev);
                  PhyML_Printf("\n== new_glnL: %f cur_glnL: %f",new_glnL,cur_glnL);
                  Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
                }
            }
          else
            {
              tree->mmod->c_lnL = cur_glnL;
            }
        }
      else
        {
          tree->mcmc->acc_move[tree->mcmc->num_move_phyrex_disk_given_ldsk]++;
        }
    }

  Free(permut);
  Free(all_disks);
}
#endif

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

#ifdef PHYREX
void MCMC_PHYREX_Indel_Hit_Serial(t_tree *tree)
{
  t_dsk *disk,*new_disk,*target_disk,*young_disk;
  t_ldsk *young_ldsk, *old_ldsk, *new_ldsk;
  int i,j,n_trials,dir_old_young,err;
  phydbl ratio, u, alpha, hr, type;
  phydbl cur_glnL, new_glnL;
  phydbl log_one_on_T;
  phydbl T,t,pindel;

  cur_glnL     = tree->mmod->c_lnL;
  new_glnL     = tree->mmod->c_lnL;
  hr           = 0.0;
  ratio        = 0.0;
  type         = -1.0;
  n_trials     = (int)(PHYREX_Total_Number_Of_Hit_Disks(tree)/10);
  T            = PHYREX_Tree_Height(tree);
  log_one_on_T = -LOG(FABS(T));
  pindel       = 0.5;

  disk = tree->disk->prev;

  For(i,n_trials)
    {
      tree->mcmc->run_move[tree->mcmc->num_move_phyrex_indel_hit_serial]++;

      t = Uni()*T;
      disk = tree->disk->prev;
      while(disk && disk->time > t) disk = disk->prev;
                                                
      assert(disk->next);
      
      cur_glnL = tree->mmod->c_lnL;
      new_glnL = tree->mmod->c_lnL;
      hr       = 0.0;
      ratio    = 0.0;

      type = Uni();
      
      if(type < pindel) /* Insert */
        {
          hr += LOG(FABS((t - disk->next->time)/T));
          hr -= log_one_on_T;
          
          young_disk = disk->next;
          assert(young_disk->n_ldsk_a);
          
          /* Which lineage is going to be hit ? */
          hr -= LOG(1./young_disk->n_ldsk_a);

          new_disk = PHYREX_Make_Disk_Event(tree->mmod->n_dim,tree->n_otu);
          PHYREX_Init_Disk_Event(new_disk,tree->mmod->n_dim,tree->mmod);
          new_disk->time = t;
          
          young_ldsk = young_disk->ldsk_a[Rand_Int(0,young_disk->n_ldsk_a-1)];
          old_ldsk   = young_ldsk->prev;

          if(old_ldsk->disk->time > young_ldsk->disk->time) 
            {
              PhyML_Printf("\n== young_ldsk: %f",young_ldsk->disk->time);
              PhyML_Printf("\n== old_ldsk: %f",old_ldsk->disk->time);
              assert(FALSE);
            }
          
          new_glnL -= PHYREX_Lk_Range(young_disk->prev,old_ldsk->disk,tree);
          /* new_glnL -= PHYREX_Lk_Range(tree->disk->prev,NULL,tree); */

          /* Direction from old to young ldsk */
          dir_old_young = PHYREX_Get_Next_Direction(young_ldsk,old_ldsk);
          assert(dir_old_young != -1);
          
          /* Make and initialize new ldsk */
          new_ldsk = PHYREX_Make_Lindisk_Node(tree->mmod->n_dim);
          PHYREX_Init_Lindisk_Node(new_ldsk,new_disk,tree->mmod->n_dim);
          PHYREX_Make_Lindisk_Next(new_ldsk);
          new_disk->ldsk = new_ldsk;
          
          /* Connect it */
          new_ldsk->prev                = old_ldsk;
          new_ldsk->next[0]             = young_ldsk;
          young_ldsk->prev              = new_ldsk;  
          old_ldsk->next[dir_old_young] = new_ldsk;
          
          /* Insert disk */
          PHYREX_Insert_Disk(new_disk,tree);

          For(j,tree->mmod->n_dim)
            {
              new_ldsk->coord->lonlat[j] = Rnorm_Trunc(0.5*(young_ldsk->coord->lonlat[j]+old_ldsk->coord->lonlat[j]),
                                                       1.*tree->mmod->rad,
                                                       0.0,tree->mmod->lim->lonlat[j],&err);
              
              hr -= Log_Dnorm_Trunc(new_ldsk->coord->lonlat[j],
                                    0.5*(young_ldsk->coord->lonlat[j]+old_ldsk->coord->lonlat[j]),
                                    1.*tree->mmod->rad,
                                    0.0,tree->mmod->lim->lonlat[j],&err);
            }
          
          /* Sample position of the center of new_disk */
          For(j,tree->mmod->n_dim)
            {
              new_disk->centr->lonlat[j] = Rnorm_Trunc(new_ldsk->coord->lonlat[j],
                                                       1.*tree->mmod->rad,
                                                       0.0,tree->mmod->lim->lonlat[j],&err);

              hr -= Log_Dnorm_Trunc(new_disk->centr->lonlat[j],
                                    new_ldsk->coord->lonlat[j],
                                    1.*tree->mmod->rad,
                                    0.0,tree->mmod->lim->lonlat[j],&err);
            }
          

          new_glnL += PHYREX_Lk_Range(young_disk->prev,old_ldsk->disk,tree);
          tree->mmod->c_lnL = new_glnL;

          ratio  = (new_glnL - cur_glnL);
          ratio += hr;
          
          ratio = EXP(ratio);
          alpha = MIN(1.,ratio);
          
          /* Always accept move */
          if(tree->mcmc->always_yes == YES && new_glnL > UNLIKELY) alpha = 1.0;
          
          u = Uni();
          
          if(u > alpha) /* Reject */
            {
              PHYREX_Remove_Disk(new_disk);
              Free_Disk(new_disk);             
              Free_Ldisk(new_ldsk);
              
              old_ldsk->next[dir_old_young] = young_ldsk;
              young_ldsk->prev = old_ldsk;

              if(tree->mmod->safe_phyrex == YES)
                {
                  new_glnL = PHYREX_Lk(tree);
                  if(Are_Equal(new_glnL,cur_glnL,1.E-3) == NO)
                    {
                      PhyML_Printf("\n== new_glnL: %f cur_glnL: %f",new_glnL,cur_glnL);
                      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
                    }
                }
              else
                {
                  PHYREX_Lk_Range(young_disk->prev,old_ldsk->disk,tree);
                  tree->mmod->c_lnL = cur_glnL;
                }
            }
          else
            {
              tree->mcmc->acc_move[tree->mcmc->num_move_phyrex_indel_hit_serial]++;
            }
        }

      else /* Remove disk */

        {
          if(disk->ldsk == NULL || (disk->ldsk != NULL && disk->ldsk->n_next > 1)) continue;

          target_disk = disk;

          old_ldsk   = target_disk->ldsk->prev;
          young_ldsk = target_disk->ldsk->next[0];

          dir_old_young = PHYREX_Get_Next_Direction(young_ldsk,old_ldsk);
          assert(dir_old_young != -1);

          /* Part of the Hastings ratio corresponding to the probability of selecting */
          /* one of disk->n_ldsk_a to be hit */ 
          hr += LOG(1./target_disk->next->n_ldsk_a);
          
          /* Density for position of the displaced ldsk */
          For(j,tree->mmod->n_dim)
            {
              hr += Log_Dnorm_Trunc(target_disk->ldsk->coord->lonlat[j],
                                    0.5*(young_ldsk->coord->lonlat[j]+old_ldsk->coord->lonlat[j]),
                                    1.0*tree->mmod->rad,
                                    0.0,
                                    tree->mmod->lim->lonlat[j],&err);
            }
          
          /* Density for position of the center of target_disk */
          For(j,tree->mmod->n_dim)
            {
              hr += Log_Dnorm_Trunc(target_disk->centr->lonlat[j],
                                    target_disk->ldsk->coord->lonlat[j],
                                    1.0*tree->mmod->rad,
                                    0.0,tree->mmod->lim->lonlat[j],&err);
            }

          new_glnL -= PHYREX_Lk_Range(target_disk,target_disk->ldsk->prev->disk,tree);
          /* new_glnL -= PHYREX_Lk_Range(tree->disk->prev,NULL,tree); */

          /* New connections between old_ldsk and young_ldsk */
          old_ldsk->next[dir_old_young] = young_ldsk;
          young_ldsk->prev              = old_ldsk;

          hr -= LOG(FABS((target_disk->time - target_disk->next->time)/T));
          hr += log_one_on_T;
          
          PHYREX_Remove_Disk(target_disk);
          new_glnL += PHYREX_Lk_Range(target_disk->prev,target_disk->ldsk->prev->disk,tree);
          tree->mmod->c_lnL = new_glnL;
                              
          ratio  = (new_glnL - cur_glnL);
          ratio += hr;
          
          ratio = EXP(ratio);
          alpha = MIN(1.,ratio);
          
          /* Always accept move */
          if(tree->mcmc->always_yes == YES && new_glnL > UNLIKELY) alpha = 1.0;
          
          u = Uni();
          
          if(u > alpha) /* Reject */
            {
              PHYREX_Insert_Disk(target_disk,tree);
              old_ldsk->next[dir_old_young] = target_disk->ldsk;
              young_ldsk->prev              = target_disk->ldsk;
              
              if(tree->mmod->safe_phyrex == YES)
                {
                  new_glnL = PHYREX_Lk(tree);
                  if(Are_Equal(new_glnL,cur_glnL,1.E-3) == NO)
                    {
                      PhyML_Printf("\n== new_glnL: %f cur_glnL: %f",new_glnL,cur_glnL);
                      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
                    }
                }
              else
                {
                  PHYREX_Lk_Range(target_disk->prev,target_disk->ldsk->prev->disk,tree);
                  tree->mmod->c_lnL = cur_glnL;
                }
            }
          else
            {
              Free_Ldisk(disk->ldsk);
              Free_Disk(disk);
              tree->mcmc->acc_move[tree->mcmc->num_move_phyrex_indel_hit_serial]++;
            }
        }
    }
  
}
#endif

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

#ifdef PHYREX
void MCMC_PHYREX_Indel_Disk_Serial(t_tree *tree)
{
  t_dsk *disk,*new_disk;
  int i,j,n_trials;
  phydbl ratio, u, alpha, hr, type;
  phydbl cur_glnL, new_glnL;
  phydbl log_lk_centr,log_one_on_T;
  phydbl T,t,pindel;

  cur_glnL     = tree->mmod->c_lnL;
  new_glnL     = tree->mmod->c_lnL;
  hr           = 0.0;
  ratio        = 0.0;
  type         = -1.0;
  n_trials     = (int)(PHYREX_Total_Number_Of_Intervals(tree)/2);
  T            = PHYREX_Tree_Height(tree);
  log_one_on_T = -LOG(FABS(T));
  pindel       = 0.5;

  log_lk_centr = 0.0;
  For(j,tree->mmod->n_dim) log_lk_centr += LOG(1./tree->mmod->lim->lonlat[j]);

  disk = tree->disk->prev;

  For(i,n_trials)
  /* do */
    {
      t = Uni()*T;
      disk = tree->disk->prev;
      while(disk && disk->time > t) disk = disk->prev;
                                                
      assert(disk->next);

      /* prev_disk = disk->prev; */

      cur_glnL = tree->mmod->c_lnL;
      new_glnL = tree->mmod->c_lnL;
      hr       = 0.0;
      ratio    = 0.0;

      type = Uni();
      
      if(type < pindel) /* Insert */
        {
          hr += LOG(FABS((t - disk->next->time)/T));
          hr -= log_one_on_T;
          hr -= log_lk_centr;
          
          new_glnL -= PHYREX_Lk_Range(disk,disk,tree);
          new_disk = PHYREX_Make_Disk_Event(tree->mmod->n_dim,tree->n_otu);
          PHYREX_Init_Disk_Event(new_disk,tree->mmod->n_dim,tree->mmod);
          new_disk->time = t;
          PHYREX_Insert_Disk(new_disk,tree);

          For(j,tree->mmod->n_dim) new_disk->centr->lonlat[j] = Uni()*tree->mmod->lim->lonlat[j];

          new_glnL += PHYREX_Lk_Range(new_disk,disk,tree);
          tree->mmod->c_lnL = new_glnL;
          
          /* new_glnL = PHYREX_Lk(tree); */
                    
          ratio  = (new_glnL - cur_glnL);
          ratio += hr;
          
          ratio = EXP(ratio);
          alpha = MIN(1.,ratio);
          
          /* Always accept move */
          if(tree->mcmc->always_yes == YES && new_glnL > UNLIKELY) alpha = 1.0;
          
          u = Uni();
          
          if(u > alpha) /* Reject */
            {
              PHYREX_Remove_Disk(new_disk);
              Free_Disk(new_disk);
              
              if(tree->mmod->safe_phyrex == YES)
                {
                  new_glnL = PHYREX_Lk(tree);
                  if(Are_Equal(new_glnL,cur_glnL,1.E-3) == NO)
                    {
                      PhyML_Printf("\n== new_glnL: %f cur_glnL: %f",new_glnL,cur_glnL);
                      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
                    }
                }
              else
                {
                  tree->mmod->c_lnL = cur_glnL;
                }
            }
        }
      else /* Remove disk */
        {
          if(disk->ldsk != NULL) continue; /* Reject */

          hr -= LOG(FABS((disk->time - disk->next->time)/T));
          hr += log_one_on_T;
          hr += log_lk_centr;

          new_glnL -= PHYREX_Lk_Range(disk,disk->prev,tree);
          PHYREX_Remove_Disk(disk);
          new_glnL += PHYREX_Lk_Range(disk->prev,disk->prev,tree);
          tree->mmod->c_lnL = new_glnL;
                    

          ratio  = (new_glnL - cur_glnL);
          ratio += hr;
          
          ratio = EXP(ratio);
          alpha = MIN(1.,ratio);
          
          /* Always accept move */
          if(tree->mcmc->always_yes == YES && new_glnL > UNLIKELY) alpha = 1.0;
          
          u = Uni();
          
          if(u > alpha) /* Reject */
            {
              PHYREX_Insert_Disk(disk,tree);
              
              if(tree->mmod->safe_phyrex == YES)
                {
                  new_glnL = PHYREX_Lk(tree);
                  if(Are_Equal(new_glnL,cur_glnL,1.E-3) == NO)
                    {
                      PhyML_Printf("\n== new_glnL: %f cur_glnL: %f",new_glnL,cur_glnL);
                      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
                    }
                }
              else
                {
                  tree->mmod->c_lnL = cur_glnL;
                }
            }
          else
            {
              Free_Disk(disk);
            }
        }
    
      /* disk = prev_disk; */
    }
  
  /* while(disk->prev); */
}
#endif

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
