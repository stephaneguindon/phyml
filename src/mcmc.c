/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include "mcmc.h"


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
  cur_lnval     = log(*val);
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

      ratio = exp(ratio);
      alpha = MIN(1.,ratio);
      
      /* Always accept move */
      if(tree->mcmc->always_yes == YES && new_lnLike > UNLIKELY) alpha = 1.0;

      u = Uni();

      assert(isnan(u) == NO && isinf(fabs(u)) == NO);
      
      if(u > alpha) /* Reject */
	{
	  *val    = cur_val;
	  new_val = cur_val;
	  if(lnPrior) *lnPrior = cur_lnPrior;
	  if(lnLike)  *lnLike  = cur_lnLike;
	  Restore_Br_Len(tree);
	  if(tree->mod && tree->mod->update_eigen)
            {
              if(!Update_Eigen(tree->mod))
                {
                  PhyML_Fprintf(stderr,"\n. Problem in move %s",tree->mcmc->move_name[tree->mcmc->move_idx]);
                  Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
                }
            }
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
      for(i=0;i<N-lag;i++) rho += 
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
  for(i=0;i<N;i++) 
    if(mcmc->sampled_val[move_num*mcmc->sample_size+burnin+i] < min)
      min = mcmc->sampled_val[move_num*mcmc->sample_size+burnin+i];

  max = -INFINITY;
  for(i=0;i<N;i++) 
    if(mcmc->sampled_val[move_num*mcmc->sample_size+burnin+i] > max)
      max = mcmc->sampled_val[move_num*mcmc->sample_size+burnin+i];


  for(i=0;i<N;i++) 
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
  for(j=0;j<breaks;j++) 
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

void MCMC_Clock_R(t_tree *tree)
{
  /* t_tree *tree; */

  /* tree = mixt_tree; */
  /* do */
  /*   { */
  /*     MCMC_Single_Param_Generic(&(tree->rates->clock_r), */
  /*                               mixt_tree->rates->min_clock, */
  /*                               mixt_tree->rates->max_clock, */
  /*                               mixt_tree->mcmc->num_move_clock_r, */
  /*                               NULL,&(mixt_tree->c_lnL), */
  /*                               NULL,Wrap_Lk, */
  /*                               mixt_tree->mcmc->move_type[mixt_tree->mcmc->num_move_clock_r], */
  /*                               NO,NULL,mixt_tree,NULL); */

  /*     tree = tree->next; */
  /*   } */
  /* while(tree); */

  phydbl new_lnL_data, cur_lnL_data;
  phydbl new_lnL_rate, cur_lnL_rate;
  phydbl u, ratio, alpha;
  phydbl new_clock_r, cur_clock_r;
  phydbl min, max;
  int move_num;
  phydbl K;

  if(tree->eval_alnL == NO) return;
  if(tree->rates->clock_r_fixed == YES) return;
  
  new_clock_r  = tree->rates->clock_r;
  cur_clock_r  = tree->rates->clock_r;
  min          = tree->rates->min_clock;
  max          = tree->rates->max_clock;
  ratio        = 0.0;
  move_num     = tree->mcmc->num_move_clock_r;
  K            = tree->mcmc->tune_move[move_num];
  cur_lnL_data = tree->c_lnL;
  new_lnL_data = tree->c_lnL;
  cur_lnL_rate = tree->rates->c_lnL_rates;
  new_lnL_rate = tree->rates->c_lnL_rates;
  
  MCMC_Make_Move(&cur_clock_r,&new_clock_r,min,max,&ratio,K,tree->mcmc->move_type[move_num]);
  
  if(new_clock_r > min && new_clock_r < max)
    {      
      tree->rates->clock_r = new_clock_r;

      RATES_Update_Cur_Bl(tree);
      
      if(tree->eval_alnL == YES) new_lnL_data = Lk(NULL,tree);
      if(tree->eval_rlnL == YES) new_lnL_rate = RATES_Lk_Rates(tree);
      
      ratio += (new_lnL_data - cur_lnL_data);
      ratio += (new_lnL_rate - cur_lnL_rate);
            
      ratio = exp(ratio);
      alpha = MIN(1.,ratio);
      
      u = Uni();
      
      assert(isnan(u) == NO && isinf(fabs(u)) == NO);
      
      if(u > alpha) /* Reject */
	{
          /* PhyML_Printf("\n. reject"); */
	  tree->rates->clock_r      = cur_clock_r;
	  tree->c_lnL               = cur_lnL_data;
	  tree->rates->c_lnL_rates  = cur_lnL_rate;
          RATES_Update_Cur_Bl(tree);
        }
      else
	{
          /* PhyML_Printf("\n. accept"); */
          tree->mcmc->acc_move[move_num]++;
	}
    }
  tree->mcmc->run_move[move_num]++;
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
      for(i=0;i<tree->geo->ldscape_sz;i++) sum += tree->geo->idx_loc_beneath[tree->n_root->num * tree->geo->ldscape_sz + i];
      for(i=0;i<tree->geo->ldscape_sz;i++) probs[i] = tree->geo->idx_loc_beneath[tree->n_root->num * tree->geo->ldscape_sz + i]/sum;
      
      tree->geo->idx_loc[tree->n_root->num] = Sample_i_With_Proba_pi(probs,tree->geo->ldscape_sz);      

      Free(probs);
    }


  // Randomize the locations below the selected node
  GEO_Randomize_Locations(tree->a_nodes[target], 
                          tree->geo,
                          tree);
  
  new_lnL = GEO_Lk(tree->geo,tree);

  ratio = (new_lnL - cur_lnL);        
  ratio = exp(ratio);
  alpha = MIN(1.,ratio);      
  u = Uni();

  assert(isnan(u) == NO && isinf(fabs(u)) == NO);
  
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

  for(i=0;i<dim;i++) tree->rates->mean_r[i] = 1.0;

  RATES_Fill_Lca_Table(tree);
  RATES_Covariance_Mu(tree);

  T = .0;
  for(i=0;i<dim;i++) T += (t[tree->a_nodes[i]->num] - t[tree->a_nodes[i]->anc->num]);
  for(i=0;i<dim;i++) lambda[i] = (t[tree->a_nodes[i]->num] - t[tree->a_nodes[i]->anc->num])/T;
  for(i=0;i<dim;i++) r[i] = 1.0;
  for(i=0;i<dim;i++) min_r[i] = tree->rates->min_rate;
  for(i=0;i<dim;i++) max_r[i] = tree->rates->max_rate;

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

void MCMC_Rates_All(t_tree *tree)
{
  Set_Both_Sides(NO,tree);
  if(tree->eval_alnL == YES) Lk(NULL,tree);
  MCMC_One_Rate(tree->n_root,tree->n_root->v[1],YES,tree);
  Update_Partial_Lk(tree,tree->e_root,tree->n_root->v[1]);
  MCMC_One_Rate(tree->n_root,tree->n_root->v[2],YES,tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MCMC_One_Rate(t_node *a, t_node *d, int traversal, t_tree *tree)
{
  t_edge *b;
  int i;
  phydbl u;
  phydbl new_lnL_seq, cur_lnL_seq, new_lnL_rate, cur_lnL_rate;
  phydbl ratio, alpha;
  phydbl new_mu, cur_mu;
  phydbl r_min, r_max;
  int move_num,err;
  phydbl K;
  
  if(tree->rates->model == STRICTCLOCK) return;
  if(tree->rates->model == GUINDON) return;

  b = NULL;
  if(a == tree->n_root) b = tree->e_root;
  else for(i=0;i<3;i++) if(d->v[i] == a) { b = d->b[i]; break; }
   
  cur_mu       = tree->rates->br_r[d->num];
  r_min        = tree->rates->min_rate;
  r_max        = tree->rates->max_rate;
  ratio        = 0.0;
  move_num     = tree->mcmc->num_move_br_r;
  K            = tree->mcmc->tune_move[move_num];
  cur_lnL_seq  = tree->c_lnL;
  new_lnL_seq  = tree->c_lnL;
  cur_lnL_rate = tree->rates->c_lnL_rates;
  new_lnL_rate = tree->rates->c_lnL_rates;
  
  /* MCMC_Make_Move(&cur_mu,&new_mu,r_min,r_max,&ratio,K,tree->mcmc->move_type[move_num]); */

  /* new_mu = Rnorm_Trunc(1.0,3.*tree->rates->nu,r_min,r_max,&err); */
  /* ratio -= Log_Dnorm_Trunc(new_mu,1.0,3.*tree->rates->nu,r_min,r_max,&err); */
  /* ratio += Log_Dnorm_Trunc(cur_mu,1.0,3.*tree->rates->nu,r_min,r_max,&err); */

  err = NO;
  new_mu = Rnorm_Trunc(cur_mu,K*tree->rates->nu,r_min,r_max,&err);

  if(err == YES) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
    
  ratio -= Log_Dnorm_Trunc(new_mu,cur_mu,K*tree->rates->nu,r_min,r_max,&err);
  ratio += Log_Dnorm_Trunc(cur_mu,new_mu,K*tree->rates->nu,r_min,r_max,&err);
  
  if(new_mu > r_min && new_mu < r_max)
    {      
      tree->rates->br_r[d->num] = new_mu;      
      tree->rates->br_do_updt[d->num] = YES;
      
      if(tree->eval_alnL == YES) new_lnL_seq  = Lk(b,tree);
      if(tree->eval_rlnL == YES) new_lnL_rate = RATES_Lk_Rates(tree);
      
      ratio += (new_lnL_seq - cur_lnL_seq);
      ratio += (new_lnL_rate - cur_lnL_rate);
            
      ratio = exp(ratio);
      alpha = MIN(1.,ratio);
      
      u = Uni();
      
      /* PhyML_Printf("\n. %f->%f %f->%f %f->%f ratio:%f", */
      /*              cur_mu,new_mu, */
      /*              cur_lnL_seq,new_lnL_seq, */
      /*              cur_lnL_rate,new_lnL_rate, */
      /*              ratio); */

      assert(isnan(u) == NO && isinf(fabs(u)) == NO);
      
      if(u > alpha) /* Reject */
	{
	  tree->rates->br_r[d->num] = cur_mu;	  
	  tree->c_lnL               = cur_lnL_seq;
	  tree->rates->c_lnL_rates  = cur_lnL_rate;

          RATES_Update_Cur_Bl(tree);
	  Update_PMat_At_Given_Edge(b,tree);          
        }
      else
	{
          tree->mcmc->acc_move[move_num]++;
	}
    }

  tree->mcmc->run_move[move_num]++;
  
  if(traversal == YES)
    {
      if(d->tax == YES) return;
      else
	{
	  for(i=0;i<3;++i)
	    if(d->v[i] != a && d->b[i] != tree->e_root)
              {
                Update_Partial_Lk(tree,d->b[i],d);
                MCMC_One_Rate(d,d->v[i],YES,tree);
              }
          Update_Partial_Lk(tree,b,d);
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
    for(i=0;i<3;i++) if(d->v[i] == a) { b = d->b[i]; break; }

  /* Only the log_RANDWALK move seems to work here. Change with caution then. */
  tree->rates->br_do_updt[d->num] = YES;
  MCMC_Single_Param_Generic(&(tree->rates->nd_r[d->num]),
			    tree->rates->min_rate,
			    tree->rates->max_rate,
			    tree->mcmc->num_move_nd_r,
			    &(tree->rates->c_lnL_rates),NULL,
			    Wrap_Lk_Rates,NULL,
			    tree->mcmc->move_type[tree->mcmc->num_move_nd_r],
			    NO,NULL,tree,NULL);


  Update_PMat_At_Given_Edge(b,tree);

  if(traversal == YES)
    {
      if(d->tax == YES) return;
      else
	{
	  for(i=0;i<3;i++)
	    if(d->v[i] != a && d->b[i] != tree->e_root)
	      {
		MCMC_One_Node_Rate(d,d->v[i],YES,tree);
	      }
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MCMC_Times_And_Rates_All(t_tree *tree)
{
  MCMC_Times_And_Rates_Root(tree);
  MCMC_Times_And_Rates_Recur(tree->n_root,tree->n_root->v[1],YES,tree);
  MCMC_Times_And_Rates_Recur(tree->n_root,tree->n_root->v[2],YES,tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MCMC_Times_And_Rates_Root(t_tree *tree)
{
  phydbl u;
  phydbl t_min,t_max;
  phydbl r_min,r_max;
  phydbl t1_cur, t1_new;
  phydbl cur_lnL_rate, new_lnL_rate;
  phydbl cur_lnL_time, new_lnL_time;
  phydbl cur_lnL_seq, new_lnL_seq;
  phydbl ratio,alpha;
  phydbl t0,t2,t3;
  t_node *v2,*v3;
  phydbl K;
  int move_num;
  t_node *root;
  phydbl r2_cur,r3_cur;
  phydbl r2_new,r3_new;

  root = tree->n_root;

  if(FABS(tree->rates->t_prior_min[root->num] - tree->rates->t_prior_max[root->num]) < 1.E-10) return;
  
  move_num     = tree->mcmc->num_move_times_and_rates_root;
  K            = tree->mcmc->tune_move[move_num];
  t1_cur       = tree->rates->nd_t[root->num];
  ratio        = 0.0;
  cur_lnL_rate = tree->rates->c_lnL_rates;
  new_lnL_rate = tree->rates->c_lnL_rates;
  cur_lnL_time = tree->rates->c_lnL_times;
  new_lnL_time = tree->rates->c_lnL_times;
  cur_lnL_seq  = tree->c_lnL;
  new_lnL_seq  = tree->c_lnL;
  r_min        = tree->rates->min_rate;
  r_max        = tree->rates->max_rate;

  v2 = root->v[2];
  v3 = root->v[1];
  
  t0 = tree->rates->t_prior_min[root->num];
  t2 = tree->rates->nd_t[v2->num];
  t3 = tree->rates->nd_t[v3->num];

  t_min = -INFINITY;
  t_max = MIN(t2,t3);

  t_min += tree->rates->min_dt;
  t_max -= tree->rates->min_dt;

  u = Uni();
  t1_new = t_max - (t_max - t1_cur)*exp(K*(u-.5));
  ratio += log((t1_new - t_max) / (t1_cur - t_max));

  r2_cur = r3_cur = -1.0;
  
  if(tree->rates->model == THORNE ||
     tree->rates->model == LOGNORMAL ||
     tree->rates->model == STRICTCLOCK)
    {
      r2_cur = tree->rates->br_r[v2->num];
      r3_cur = tree->rates->br_r[v3->num];
    }
  else if(tree->rates->model == GUINDON)
    {
      r2_cur = tree->rates->nd_r[v2->num];
      r3_cur = tree->rates->nd_r[v3->num];
    }
  else assert(FALSE);

  
  r2_new = r2_cur * (t2 - t1_cur) / (t2 - t1_new);
  r3_new = r3_cur * (t3 - t1_cur) / (t3 - t1_new);


  if(t_min > t_max) 
    {
      PhyML_Fprintf(stderr,"\n. glnL:%f",TIMES_Lk_Times(NO,tree));
      PhyML_Fprintf(stderr,"\n. t:%f",tree->rates->nd_t[tree->n_root->num]);
      PhyML_Fprintf(stderr,"\n. t0 = %f t2 = %f t3 = %f",t0,t2,t3);
      PhyML_Fprintf(stderr,"\n. t_min = %f t_max = %f",t_min,t_max);
      PhyML_Fprintf(stderr,"\n. prior_min = %f prior_max = %f",tree->rates->t_prior_min[root->num],tree->rates->t_prior_max[root->num]);
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
    }



  if(t1_new > t_min &&
     t1_new < t_max &&
     r2_new > r_min &&
     r2_new < r_max &&
     r3_new > r_min &&
     r3_new < r_max) 
    {
      RATES_Record_Times(tree);

      tree->rates->nd_t[root->num] = t1_new;
      
      if(tree->rates->model == THORNE ||
         tree->rates->model == LOGNORMAL ||
         tree->rates->model == STRICTCLOCK)
        {
          tree->rates->br_r[v2->num]   = r2_new;
          tree->rates->br_r[v3->num]   = r3_new;
        }
      else if(tree->rates->model == GUINDON)
        {
          tree->rates->nd_r[v2->num]   = r2_new;
          tree->rates->nd_r[v3->num]   = r3_new;
        }
      else assert(FALSE);
      
      RATES_Update_Cur_Bl(tree);
         
      if(tree->eval_glnL == YES) new_lnL_time = TIMES_Lk_Times(NO,tree); 
      if(tree->eval_rlnL == YES) new_lnL_rate = RATES_Lk_Rates(tree);

      if(tree->rates->model == GUINDON)
        {
          RATES_Update_Cur_Bl(tree);
          if(tree->eval_alnL == YES) new_lnL_seq = Lk(NULL,tree);        
        }
      else new_lnL_seq = cur_lnL_seq;
      
      ratio += (new_lnL_rate - cur_lnL_rate);
      ratio += (new_lnL_time - cur_lnL_time);
      ratio += (new_lnL_seq  - cur_lnL_seq);
      
      ratio = exp(ratio);
      alpha = MIN(1.,ratio);
      u = Uni();
         
      /* printf("\n. alnL:%f->%f tlnL: %f->%f rlnL: %f->%f ratio: %f", */
      /*        cur_lnL_seq,new_lnL_seq, */
      /*        cur_lnL_time,new_lnL_time, */
      /*        cur_lnL_rate,new_lnL_rate, */
      /*        ratio); */

      assert(isnan(u) == NO && isinf(fabs(u)) == NO);

      if(u > alpha) /* Reject */
	{
          tree->rates->nd_t[root->num] = t1_cur;

          if(tree->rates->model == THORNE ||
             tree->rates->model == LOGNORMAL ||
             tree->rates->model == STRICTCLOCK)
            {
              tree->rates->br_r[v2->num] = r2_cur;
              tree->rates->br_r[v3->num] = r3_cur;
            }
          else if(tree->rates->model == GUINDON)
            {
              tree->rates->nd_r[v2->num] = r2_cur;
              tree->rates->nd_r[v3->num] = r3_cur;
            }
          else assert(FALSE);
          
          RATES_Update_Cur_Bl(tree);

          tree->rates->c_lnL_rates = cur_lnL_rate;
	  tree->rates->c_lnL_times = cur_lnL_time;
	  tree->c_lnL              = cur_lnL_seq;
          
	}
      else
	{
	  tree->mcmc->acc_move[move_num]+=1;
	}

      // Ignore boundaries when updating tuning parameter 
      tree->mcmc->run_move[move_num]+=1;

    }   
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MCMC_Times_And_Rates_Recur(t_node *a, t_node *d, int traversal, t_tree *tree)
{
  phydbl u;
  phydbl t_min,t_max;
  phydbl r_min,r_max;
  phydbl t1_cur, t1_new;
  phydbl cur_lnL_rate, new_lnL_rate;
  phydbl cur_lnL_time, new_lnL_time;
  phydbl cur_lnL_seq, new_lnL_seq;
  phydbl ratio,alpha;
  int    i;
  phydbl t0,t2,t3;
  t_node *v2,*v3;
  int move_num;
  phydbl r1_cur,r2_cur,r3_cur;
  phydbl r1_new,r2_new,r3_new;
  
  if(d->tax) return; /* Won't change time at tip */
  
  move_num     = tree->mcmc->num_move_times_and_rates;
  t1_cur       = tree->rates->nd_t[d->num];
  ratio        = 0.0;
  cur_lnL_rate = tree->rates->c_lnL_rates;
  new_lnL_rate = tree->rates->c_lnL_rates;
  cur_lnL_time = tree->rates->c_lnL_times;
  new_lnL_time = tree->rates->c_lnL_times;
  cur_lnL_seq  = tree->c_lnL;
  new_lnL_seq  = tree->c_lnL;
  r_min        = tree->rates->min_rate;
  r_max        = tree->rates->max_rate;
  r2_cur       = -1.0;
  r1_cur       = -1.0;
  r3_cur       = -1.0;
  
  
  v2 = v3 = NULL;
  for(i=0;i<3;++i)
    if((d->v[i] != a) && (d->b[i] != tree->e_root))
      {
	if(!v2) { v2 = d->v[i]; }
	else    { v3 = d->v[i]; }
      }
  
  if(tree->rates->model == THORNE ||
     tree->rates->model == LOGNORMAL ||
     tree->rates->model == STRICTCLOCK)
    {
      r1_cur = tree->rates->br_r[d->num];
      r2_cur = tree->rates->br_r[v2->num];
      r3_cur = tree->rates->br_r[v3->num];
    }
  else if(tree->rates->model == GUINDON)
    {
      r1_cur = tree->rates->nd_r[d->num];
      r2_cur = tree->rates->nd_r[v2->num];
      r3_cur = tree->rates->nd_r[v3->num];
    }
  else assert(FALSE);

  t0 = tree->rates->nd_t[a->num];
  t2 = tree->rates->nd_t[v2->num];
  t3 = tree->rates->nd_t[v3->num];

  t_min = MAX(t0,tree->rates->t_prior_min[d->num]);
  t_max = MIN(tree->rates->t_prior_max[d->num],MIN(t2,t3));

  t_min += tree->rates->min_dt;
  t_max -= tree->rates->min_dt;

  t1_new = Uni()*(t_max - t_min) + t_min;

  r1_new = r1_cur * (t1_cur - t0) / (t1_new - t0);
  r2_new = r2_cur * (t2 - t1_cur) / (t2 - t1_new);
  r3_new = r3_cur * (t3 - t1_cur) / (t3 - t1_new);

  ratio +=
    log((t1_cur - t0) / (t1_new - t0) *
        (t2 - t1_cur) / (t2 - t1_new) *
        (t3 - t1_cur) / (t3 - t1_new));
  
  if(t1_new > t_min && t1_new < t_max &&
     r1_new > r_min && r1_new < r_max &&
     r2_new > r_min && r2_new < r_max &&
     r3_new > r_min && r3_new < r_max) 
    {
      RATES_Record_Times(tree);

      tree->rates->nd_t[d->num] = t1_new;

      if(tree->rates->model == THORNE ||
         tree->rates->model == LOGNORMAL ||
         tree->rates->model == STRICTCLOCK)
        {
          tree->rates->br_r[d->num]  = r1_new;
          tree->rates->br_r[v2->num] = r2_new;
          tree->rates->br_r[v3->num] = r3_new;
        }
      else if(tree->rates->model == GUINDON)
        {
          tree->rates->nd_r[d->num]  = r1_new;
          tree->rates->nd_r[v2->num] = r2_new;
          tree->rates->nd_r[v3->num] = r3_new;          
        }
      else assert(FALSE);

      if(tree->eval_glnL == YES) new_lnL_time = TIMES_Lk_Times(NO,tree);
      if(tree->eval_rlnL == YES) new_lnL_rate = RATES_Lk_Rates(tree);

      if(tree->rates->model == GUINDON)
        {
          RATES_Update_Cur_Bl(tree);
          if(tree->eval_alnL == YES) new_lnL_seq = Lk(NULL,tree);
        }
      else new_lnL_seq = cur_lnL_seq;

      ratio += (new_lnL_time - cur_lnL_time);
      ratio += (new_lnL_rate - cur_lnL_rate);
      ratio += (new_lnL_seq - cur_lnL_seq);

      
      ratio = exp(ratio);
      alpha = MIN(1.,ratio);
      u = Uni();

      assert(isnan(u) == NO && isinf(fabs(u)) == NO);
      
      if(u > alpha) /* Reject */
	{
          /* printf("\n. rej"); */
          RATES_Reset_Times(tree);
              
          if(tree->rates->model == THORNE ||
             tree->rates->model == LOGNORMAL ||
             tree->rates->model == STRICTCLOCK)
            {
              tree->rates->br_r[d->num]  = r1_cur;
              tree->rates->br_r[v2->num] = r2_cur;
              tree->rates->br_r[v3->num] = r3_cur;
            }
          else if(tree->rates->model == GUINDON)
            {
              tree->rates->nd_r[d->num]  = r1_cur;
              tree->rates->nd_r[v2->num] = r2_cur;
              tree->rates->nd_r[v3->num] = r3_cur;
            }
          else assert(FALSE);
          
          RATES_Update_Cur_Bl(tree);

          tree->rates->c_lnL_rates = cur_lnL_rate;
          tree->rates->c_lnL_times = cur_lnL_time;          
	  tree->c_lnL              = cur_lnL_seq;
	}
      else
	{
          /* printf("\n. acc"); */
	  tree->mcmc->acc_move[move_num]++;
	}
    }
  
  tree->mcmc->run_move[move_num]++;

  if(d->tax == YES) return;
  else
    {
      for(i=0;i<3;i++)
        if(d->v[i] != a && d->b[i] != tree->e_root)
          {
            MCMC_Times_And_Rates_Recur(d,d->v[i],YES,tree);
          }
    }	    
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MCMC_Times_All(t_tree *tree)
{
  // Down partials may not be up to date.
  Set_Both_Sides(YES,tree);
  if(tree->eval_alnL == YES) Lk(NULL,tree);
  Set_Both_Sides(NO,tree);
  MCMC_Root_Time(tree);
  MCMC_Time_Recur(tree->n_root,tree->n_root->v[1],YES,tree);
  if(tree->eval_alnL == YES) Update_Partial_Lk(tree,tree->e_root,tree->n_root->v[1]);
  MCMC_Time_Recur(tree->n_root,tree->n_root->v[2],YES,tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MCMC_Time_Recur(t_node *a, t_node *d, int traversal, t_tree *tree)
{
  phydbl u;
  phydbl t_min,t_max;
  phydbl t1_cur, t1_new;
  phydbl cur_lnL_seq, new_lnL_seq;
  phydbl cur_lnL_rate, new_lnL_rate;
  phydbl cur_lnL_time, new_lnL_time;
  /* phydbl K; */
  phydbl ratio,alpha;
  t_edge *b1,*b2,*b3;
  int    i;
  phydbl t0,t2,t3;
  t_node *v2,*v3;
  int move_num;
  
  if(d->tax) return; /* Won't change time at tip */
  
  move_num     = tree->mcmc->num_move_times;
  t1_cur       = tree->rates->nd_t[d->num];
  ratio        = 0.0;
  cur_lnL_seq = tree->c_lnL;
  new_lnL_seq = tree->c_lnL;
  cur_lnL_rate = tree->rates->c_lnL_rates;
  new_lnL_rate = tree->rates->c_lnL_rates;
  cur_lnL_time = tree->rates->c_lnL_times;
  new_lnL_time = tree->rates->c_lnL_times;

  
  v2 = v3 = NULL;
  for(i=0;i<3;i++)
    if((d->v[i] != a) && (d->b[i] != tree->e_root))
      {
	if(!v2) { v2 = d->v[i]; }
	else    { v3 = d->v[i]; }
      }


  b1 = NULL;
  if(a == tree->n_root) b1 = tree->e_root;
  else for(i=0;i<3;i++) if(d->v[i] == a) { b1 = d->b[i]; break; }

  b2 = b3 = NULL;
  for(i=0;i<3;i++)
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

  /* K = tree->mcmc->tune_move[tree->mcmc->num_move_times] * (t_max - t_min);   */
  /* MCMC_Make_Move(&t1_cur,&t1_new,t_min,t_max,&ratio,K,tree->mcmc->move_type[move_num]); */
  t1_new = Uni()*(t_max - t_min) + t_min;

  if(t1_new > t_min && t1_new < t_max) 
    {
      RATES_Record_Times(tree);

      tree->rates->nd_t[d->num] = t1_new;

      if(tree->eval_glnL == YES) new_lnL_time = TIMES_Lk_Times(NO,tree);
      ratio += (new_lnL_time - cur_lnL_time);

      if(new_lnL_time > UNLIKELY)
        {
          RATES_Update_Cur_Bl(tree);
          
          if(tree->eval_rlnL == YES) new_lnL_rate = RATES_Lk_Rates(tree);
          ratio += (new_lnL_rate - cur_lnL_rate);
          
          if(tree->eval_alnL == YES && tree->io->lk_approx == EXACT)
            {
              Update_PMat_At_Given_Edge(b1,tree);
              Update_PMat_At_Given_Edge(b2,tree);
              Update_PMat_At_Given_Edge(b3,tree);
              Update_Partial_Lk(tree,b1,d);
            }
          if(tree->eval_alnL == YES) new_lnL_seq = Lk(b1,tree);
          
          ratio += (new_lnL_seq - cur_lnL_seq);
        }
      
      /* printf("\n. One_Time cur_t: %f new_t: %f ratio: %f alnL:%f->%f tlnL: %f->%f rlnL: %f->%f", */
      /*        t1_cur,t1_new, */
      /*        ratio, */
      /*        cur_lnL_seq,new_lnL_seq, */
      /*        cur_lnL_time,new_lnL_time, */
      /*        cur_lnL_rate,new_lnL_rate); */

      ratio = exp(ratio);
      alpha = MIN(1.,ratio);
      u = Uni();

      assert(isnan(u) == NO && isinf(fabs(u)) == NO);

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
              Update_Partial_Lk(tree,b1,d);
            }

          if(isinf(FABS(new_lnL_time)) == YES || isnan(new_lnL_time) == YES)
            {
              Print_Node(tree->n_root,tree->n_root->v[1],tree);
              Print_Node(tree->n_root,tree->n_root->v[2],tree);              
              assert(FALSE);
            }

	  tree->c_lnL               = cur_lnL_seq;
	  tree->rates->c_lnL_rates  = cur_lnL_rate;
          tree->rates->c_lnL_times  = cur_lnL_time;

          if(Are_Equal(tree->rates->c_lnL_times,cur_lnL_time,1.E-3) == NO)
            {
              PhyML_Fprintf(stderr,"\n\n");
              PhyML_Fprintf(stderr,"\n. moved node %d from %f to %f\n",d->num,t1_cur,t1_new);
              Print_Node(tree->n_root,tree->n_root->v[1],tree);
              Print_Node(tree->n_root,tree->n_root->v[2],tree);              
              PhyML_Fprintf(stderr,"\n. new_glnL: %f cur_glnL: %f",tree->rates->c_lnL_times,cur_lnL_time);
              Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
            }
          
          /* new_lnL_seq = Lk(NULL,tree); /\* Not necessary. Remove once tested *\/ */
          /* if(Are_Equal(new_lnL_seq,cur_lnL_seq,1.E-3) == NO) */
          /*   { */
          /*     PhyML_Printf("\n. a: %d d: %d v2: %d v3: %d",a->num,d->num,v2->num,v3->num); */
          /*     PhyML_Printf("\n. t1_cur: %f t1_new: %f",t1_cur,t1_new); */
          /*     PhyML_Printf("\n. new_alnL: %f cur_alnL: %f",new_lnL_seq,cur_lnL_seq); */
          /*     Generic_Exit(__FILE__,__LINE__,__FUNCTION__); */
          /*   } */
	}
      else
	{
          /* cur_lnL_seq = Lk(NULL,tree); /\* Not necessary. Remove once tested *\/ */
          /* if(Are_Equal(new_lnL_seq,cur_lnL_seq,1.E-3) == NO) */
          /*   { */
          /*     PhyML_Printf("\n. a: %d d: %d v2: %d v3: %d",a->num,d->num,v2->num,v3->num); */
          /*     PhyML_Printf("\n. t1_cur: %f t1_new: %f",t1_cur,t1_new); */
          /*     PhyML_Printf("\n. new_alnL: %f cur_alnL: %f",new_lnL_seq,cur_lnL_seq); */
          /*     Generic_Exit(__FILE__,__LINE__,__FUNCTION__); */
          /*   } */
	  tree->mcmc->acc_move[move_num]++;
	}

      /* printf("\n. %f",new_lnL_seq); fflush(NULL); */

      if(t1_new < t0)
	{
	  t1_new = t0+1.E-4;
	  PhyML_Fprintf(stderr,"\n");
	  PhyML_Fprintf(stderr,"\n. a is root -> %s",(a == tree->n_root)?("YES"):("NO"));
	  PhyML_Fprintf(stderr,"\n. t0 = %f t1_new = %f",t0,t1_new);
	  PhyML_Fprintf(stderr,"\n. t_min=%f t_max=%f",t_min,t_max);
	  PhyML_Fprintf(stderr,"\n. (t1-t0)=%f (t2-t1)=%f",t1_cur-t0,t2-t1_cur);
	  PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d\n",__FILE__,__LINE__);
	  /*       Exit("\n"); */
	}
      if(t1_new > MIN(t2,t3))
	{
	  PhyML_Fprintf(stderr,"\n");
	  PhyML_Fprintf(stderr,"\n. a is root -> %s",(a == tree->n_root)?("YES"):("NO"));
	  PhyML_Fprintf(stderr,"\n. t0 = %f t1_new = %f t1 = %f t2 = %f t3 = %f MIN(t2,t3)=%f",t0,t1_new,t1_cur,t2,t3,MIN(t2,t3));
	  PhyML_Fprintf(stderr,"\n. t_min=%f t_max=%f",t_min,t_max);
	  PhyML_Fprintf(stderr,"\n. (t1-t0)=%f (t2-t1)=%f",t1_cur-t0,t2-t1_cur);
	  PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d\n",__FILE__,__LINE__);
	  /*       Exit("\n"); */
	}
      
      if(isnan(t1_new))
	{
	  PhyML_Fprintf(stderr,"\n. run=%d",tree->mcmc->run);
	  PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d\n",__FILE__,__LINE__);
	}
    }
  
  tree->mcmc->run_move[move_num]++;

  if(traversal == YES)
    {
      if(d->tax == YES) return;
      else
	{
          for(i=0;i<3;i++)
            if(d->v[i] != a && d->b[i] != tree->e_root)
              {
                if(tree->eval_alnL == YES) Update_Partial_Lk(tree,d->b[i],d);
                MCMC_Time_Recur(d,d->v[i],YES,tree);
              }
          if(tree->eval_alnL == YES) Update_Partial_Lk(tree,b1,d);
        }
    }	    
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MCMC_Root_Time(t_tree *tree)
{
  phydbl u;
  phydbl t_min,t_max;
  phydbl t1_cur, t1_new;
  phydbl cur_lnL_seq, new_lnL_seq;
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
  
  move_num     = tree->mcmc->num_move_root_time;
  K            = tree->mcmc->tune_move[move_num];
  t1_cur       = tree->rates->nd_t[root->num];
  ratio        = 0.0;
  cur_lnL_seq = tree->c_lnL;
  new_lnL_seq = tree->c_lnL;
  cur_lnL_rate = tree->rates->c_lnL_rates;
  new_lnL_rate = tree->rates->c_lnL_rates;
  cur_lnL_time = tree->rates->c_lnL_times;
  new_lnL_time = tree->rates->c_lnL_times;

  v2 = root->v[2];
  v3 = root->v[1];

  b1 = tree->e_root;
  
  t0 = tree->rates->t_prior_min[root->num];
  t2 = tree->rates->nd_t[v2->num];
  t3 = tree->rates->nd_t[v3->num];

  t_min = -INFINITY;
  t_max = MIN(t2,t3);

  t_min += tree->rates->min_dt;
  t_max -= tree->rates->min_dt;

  if(t_min > t_max) 
    {
      PhyML_Fprintf(stderr,"\n. glnL:%f",TIMES_Lk_Times(NO,tree));
      PhyML_Fprintf(stderr,"\n. t:%f",tree->rates->nd_t[tree->n_root->num]);
      PhyML_Fprintf(stderr,"\n. t0 = %f t2 = %f t3 = %f",t0,t2,t3);
      PhyML_Fprintf(stderr,"\n. t_min = %f t_max = %f",t_min,t_max);
      PhyML_Fprintf(stderr,"\n. prior_min = %f prior_max = %f",tree->rates->t_prior_min[root->num],tree->rates->t_prior_max[root->num]);
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
    }

  u = Uni();
  t1_new = t_max - (t_max - t1_cur)*exp(K*(u-.5));
  ratio += log((t1_new - t_max) / (t1_cur - t_max));
  /* MCMC_Make_Move(&t1_cur,&t1_new,t_min,t_max,&ratio,K,tree->mcmc->move_type[move_num]); */
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

      if(tree->eval_glnL == YES) new_lnL_time = TIMES_Lk_Times(NO,tree); 

      if(new_lnL_time > UNLIKELY)
        {
          if(tree->eval_alnL == YES) new_lnL_seq = Lk(b1,tree);
          if(tree->eval_rlnL == YES) new_lnL_rate = RATES_Lk_Rates(tree);
        }
      
      ratio += (new_lnL_seq - cur_lnL_seq);
      ratio += (new_lnL_rate - cur_lnL_rate);
      ratio += (new_lnL_time - cur_lnL_time);

      ratio = exp(ratio);
      alpha = MIN(1.,ratio);
      u = Uni();
         
      /* printf("\n. alnL:%f->%f tlnL: %f->%f rlnL: %f->%f ratio: %f", */
      /*        cur_lnL_seq,new_lnL_seq, */
      /*        cur_lnL_time,new_lnL_time, */
      /*        cur_lnL_rate,new_lnL_rate, */
      /*        ratio); */

      assert(isnan(u) == NO && isinf(fabs(u)) == NO);

      if(u > alpha) /* Reject */
	{
	  RATES_Reset_Times(tree);
          RATES_Update_Cur_Bl(tree);
          if(tree->eval_alnL == YES) Update_PMat_At_Given_Edge(b1,tree);

	  tree->c_lnL              = cur_lnL_seq;
	  tree->rates->c_lnL_rates = cur_lnL_rate;
	  tree->rates->c_lnL_times = cur_lnL_time;
          
          /* if(Are_Equal(tree->rates->c_lnL_times,cur_lnL_time,1.E-3) == NO) */
          /*   { */
          /*     PhyML_Printf("\n\n"); */
          /*     Print_Node(tree->n_root,tree->n_root->v[1],tree); */
          /*     Print_Node(tree->n_root,tree->n_root->v[2],tree); */
          /*     PhyML_Printf("\n. new_glnL: %f cur_glnL: %f",tree->rates->c_lnL_times,cur_lnL_time); */
          /*     Generic_Exit(__FILE__,__LINE__,__FUNCTION__); */
          /*   } */
	}
      else
	{
	  tree->mcmc->acc_move[move_num]++;
	}

      // Ignore boundaries when updating tuning parameter
      tree->mcmc->run_move[move_num]+=1;

    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MCMC_Tree_Height(t_tree *tree)
{
  phydbl K,mult,u,alpha,ratio;
  phydbl cur_lnL_seq,new_lnL_seq;
  phydbl cur_lnL_rate,new_lnL_rate;
  phydbl cur_lnL_time,new_lnL_time;
  phydbl floor;
  int n_nodes;
      
  if(FABS(tree->rates->t_prior_max[tree->n_root->num] - tree->rates->t_prior_min[tree->n_root->num]) < 1.E-10) return;

  RATES_Record_Times(tree);

  K            = tree->mcmc->tune_move[tree->mcmc->num_move_tree_height];
  ratio        = 0.0;
  cur_lnL_seq = tree->c_lnL;
  new_lnL_seq = tree->c_lnL;
  cur_lnL_rate = tree->rates->c_lnL_rates;
  new_lnL_rate = tree->rates->c_lnL_rates;
  cur_lnL_time = tree->rates->c_lnL_times;
  new_lnL_time = tree->rates->c_lnL_times;
  
  u = Uni();
  mult = exp(K*(u-0.5));

  floor = 0.0;
  Scale_Subtree_Height(tree->n_root,mult,floor,&n_nodes,tree);
  
  RATES_Update_Cur_Bl(tree);

  if(tree->eval_glnL == YES) new_lnL_time = TIMES_Lk_Times(NO,tree);

  if(new_lnL_time > UNLIKELY)
    {
      if(tree->eval_alnL == YES) new_lnL_seq = Lk(NULL,tree);
      if(tree->eval_rlnL == YES) new_lnL_rate = RATES_Lk_Rates(tree);
    }
  
  ratio += (phydbl)(n_nodes)*log(mult);

  ratio += (new_lnL_seq - cur_lnL_seq);
  ratio += (new_lnL_rate - cur_lnL_rate);
  ratio += (new_lnL_time - cur_lnL_time);

  /* printf("\n. seq: %f %f times: %f %f rates: %f %f", */
  /*        new_lnL_seq,cur_lnL_seq, */
  /*        new_lnL_time,cur_lnL_time, */
  /*        new_lnL_rate,cur_lnL_rate); */


  ratio = exp(ratio);
  alpha = MIN(1.,ratio);
  u = Uni();
  
  assert(isnan(u) == NO && isinf(fabs(u)) == NO);

  if(u > alpha)
    {
      RATES_Reset_Times(tree);
      RATES_Update_Cur_Bl(tree);
      tree->rates->c_lnL_times = TIMES_Lk_Times(NO,tree); // Required in order to set t_prior_min/max to their original values
      tree->c_lnL              = cur_lnL_seq;
      tree->rates->c_lnL_rates = cur_lnL_rate;

      if(Are_Equal(tree->rates->c_lnL_times,cur_lnL_time,1.E-3) == NO)
        {
          PhyML_Fprintf(stderr,"\n. new_glnL: %f cur_glnL: %f",tree->rates->c_lnL_times,cur_lnL_time);
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
    the different seq partitions
  */
  
  int i;
  phydbl K,mult,u,alpha,ratio;
  phydbl cur_lnL_rate,new_lnL_rate;
  phydbl cur_lnL_time,new_lnL_time;
  phydbl cur_lnL_seq,new_lnL_seq;
  phydbl floor;
  int n_nodes;

  /*! Check that sequences are isochronous. */
  for(i=0;i<tree->n_otu-1;i++) if(!Are_Equal(tree->rates->nd_t[i+1],tree->rates->nd_t[i],1.E-10)) return; 

  if(FABS(tree->rates->t_prior_max[tree->n_root->num] - tree->rates->t_prior_min[tree->n_root->num]) < 1.E-10) return;

  RATES_Record_Times(tree);

  K            = tree->mcmc->tune_move[tree->mcmc->num_move_updown_t_cr];
  ratio        = 0.0;
  cur_lnL_rate = tree->rates->c_lnL_rates;
  new_lnL_rate = tree->rates->c_lnL_rates;
  cur_lnL_time = tree->rates->c_lnL_times;
  new_lnL_time = tree->rates->c_lnL_times;
  cur_lnL_seq  = tree->c_lnL;
  new_lnL_seq  = tree->c_lnL;


  u = Uni();
  mult = exp(K*(u-0.5));

  floor = 0.0;
  Scale_Subtree_Height(tree->n_root,mult,floor,&n_nodes,tree);

  if(tree->rates->clock_r_fixed == NO)
    {
      tree->rates->clock_r /= mult;
      if(tree->rates->clock_r < tree->rates->min_clock || tree->rates->clock_r > tree->rates->max_clock)
        {
          tree->rates->clock_r *= mult;
          RATES_Reset_Times(tree);
          tree->mcmc->run_move[tree->mcmc->num_move_updown_t_cr]++;
          return;
        }
    }
  
  if(tree->eval_rlnL == YES) new_lnL_rate = RATES_Lk_Rates(tree);
  if(tree->eval_glnL == YES) new_lnL_time = TIMES_Lk_Times(NO,tree);

  /* The Hastings ratio is actually mult^(n) when changing the absolute
     node heights. When considering the relative heights, this ratio combined
     to the Jacobian for the change of variable ends up to being equal to mult. 
  */
  ratio += (n_nodes - 1)*log(mult);

  ratio += (new_lnL_rate - cur_lnL_rate);
  ratio += (new_lnL_time - cur_lnL_time);

  
  if(tree->rates->model == GUINDON)
    {
      RATES_Update_Cur_Bl(tree);
      if(tree->eval_alnL == YES) new_lnL_seq = Lk(NULL,tree);
    }
  else new_lnL_seq = cur_lnL_seq;

  ratio += (new_lnL_seq - cur_lnL_seq);
  
  ratio = exp(ratio);
  alpha = MIN(1.,ratio);
  u = Uni();
  
  /* printf("\n. t_old = %f t_new = %f cr_old = %f cr_new = %f", */
  /*        tree->rates->nd_t[tree->n_root->num]/mult, */
  /*        tree->rates->nd_t[tree->n_root->num], */
  /*        tree->rates->clock_r*mult, */
  /*        tree->rates->clock_r); */

  assert(isnan(u) == NO && isinf(fabs(u)) == NO);

  if(u > alpha)
    {
      /* printf(" reject"); */
      RATES_Reset_Times(tree);
      if(tree->rates->clock_r_fixed == NO) tree->rates->clock_r *= mult;
      tree->rates->c_lnL_rates = cur_lnL_rate;
      tree->rates->c_lnL_times = cur_lnL_time;
      if(tree->rates->model == GUINDON && tree->eval_alnL == YES)
        {
          RATES_Update_Cur_Bl(tree);
          tree->c_lnL = cur_lnL_seq;
        }
    }
  else
    {
      /* printf("\ accept"); */
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
  phydbl cur_lnL_seq,new_lnL_seq;
  phydbl cur_lnL_rate,new_lnL_rate;
  phydbl cur_lnL_time,new_lnL_time;
  phydbl floor;
  int n_nodes, target;

  /*! Check that sequences are isochronous. */
  for(i=0;i<tree->n_otu-1;i++) if(!Are_Equal(tree->rates->nd_t[i+1],tree->rates->nd_t[i],1.E-10)) return; 

  if(FABS(tree->rates->t_prior_max[tree->n_root->num] - tree->rates->t_prior_min[tree->n_root->num]) < 1.E-10) return;

  RATES_Record_Times(tree);

  K            = tree->mcmc->tune_move[tree->mcmc->num_move_updown_t_br];
  ratio        = 0.0;
  cur_lnL_seq = tree->c_lnL;
  new_lnL_seq = tree->c_lnL;
  cur_lnL_rate = tree->rates->c_lnL_rates;
  new_lnL_rate = tree->rates->c_lnL_rates;
  cur_lnL_time = tree->rates->c_lnL_times;
  new_lnL_time = tree->rates->c_lnL_times;

  u = Uni();
  mult = exp(K*(u-0.5));

  target = Rand_Int(tree->n_otu,2*tree->n_otu-3);
  floor = 0.0;
  Scale_Subtree_Height(tree->a_nodes[target],mult,floor,&n_nodes,tree);

  for(i=0;i<2*tree->n_otu-1;++i)
    {
      if(tree->rates->nd_t[i] > tree->rates->t_prior_max[i] ||
  	 tree->rates->nd_t[i] < tree->rates->t_prior_min[i])
  	{
  	  RATES_Reset_Times(tree);
	  tree->mcmc->run_move[tree->mcmc->num_move_updown_t_br]++;
  	  return;
  	}
    }

  if(RATES_Check_Node_Times(tree)) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  
  if(tree->eval_alnL == YES) new_lnL_seq = Lk(NULL,tree);
  if(tree->eval_rlnL == YES) new_lnL_rate = RATES_Lk_Rates(tree);
  if(tree->eval_glnL == YES) new_lnL_time = TIMES_Lk_Times(NO,tree);

  ratio += n_nodes*log(mult);

  /* Likelihood ratio */
  ratio += (new_lnL_seq - cur_lnL_seq);
  ratio += (new_lnL_rate - cur_lnL_rate);
  ratio += (new_lnL_time - cur_lnL_time);
  
  ratio = exp(ratio);
  alpha = MIN(1.,ratio);
  u = Uni();
  
  assert(isnan(u) == NO && isinf(fabs(u)) == NO);

  if(u > alpha)
    {
      RATES_Reset_Times(tree);
      tree->c_lnL              = cur_lnL_seq;
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
  phydbl K,mult,u,alpha,ratio;
  phydbl cur_lnL_seq,new_lnL_seq;
  phydbl cur_lnL_rate,new_lnL_rate;
  phydbl cur_lnL_time,new_lnL_time;
  phydbl floor;
  int target;
  int n_nodes;

  RATES_Record_Times(tree);

  K = tree->mcmc->tune_move[tree->mcmc->num_move_subtree_height];
  ratio        = 0.0;
  cur_lnL_seq = tree->c_lnL;
  new_lnL_seq = tree->c_lnL;
  cur_lnL_rate = tree->rates->c_lnL_rates;
  new_lnL_rate = tree->rates->c_lnL_rates;
  cur_lnL_time = tree->rates->c_lnL_times;
  new_lnL_time = tree->rates->c_lnL_times;

  u = Uni();
  mult = exp(K*(u-0.5));

  target = Rand_Int(tree->n_otu,2*tree->n_otu-3);
  floor = 0.0;
  if(tree->a_nodes[target] == tree->n_root) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);

  if(!Scale_Subtree_Height(tree->a_nodes[target],mult,floor,&n_nodes,tree))
    {
      RATES_Reset_Times(tree);
      tree->mcmc->run_move[tree->mcmc->num_move_subtree_height]++;
      return;
    }

  if(RATES_Check_Node_Times(tree)) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);

  if(tree->eval_glnL == YES) new_lnL_time = TIMES_Lk_Times(NO,tree);

  if(new_lnL_time > UNLIKELY)
    {
      if(tree->eval_alnL == YES) new_lnL_seq = Lk(NULL,tree);
      if(tree->eval_rlnL == YES) new_lnL_rate = RATES_Lk_Rates(tree);
    }
  
  ratio += (phydbl)(n_nodes)*log(mult);
  ratio += (new_lnL_seq - cur_lnL_seq);
  ratio += (new_lnL_rate - cur_lnL_rate);
  ratio += (new_lnL_time - cur_lnL_time);

  ratio = exp(ratio);
  alpha = MIN(1.,ratio);
  u = Uni();
  
  assert(isnan(u) == NO && isinf(fabs(u)) == NO);

  if(u > alpha)
    {
      RATES_Reset_Times(tree);
      RATES_Update_Cur_Bl(tree);
      tree->c_lnL              = cur_lnL_seq;
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
  phydbl cur_lnL_seq,new_lnL_seq;
  int n_nodes;
  phydbl init_clock;

  if(tree->eval_alnL == NO) return;
  if(tree->rates->model == STRICTCLOCK) return;

  RATES_Record_Rates(tree);

  tree->mcmc->run_move[tree->mcmc->num_move_tree_rates]++;

  K             = tree->mcmc->tune_move[tree->mcmc->num_move_tree_rates];
  cur_lnL_rate  = tree->rates->c_lnL_rates;
  new_lnL_rate  = UNLIKELY;
  init_clock    = tree->rates->clock_r;
  ratio         = 0.0;
  cur_lnL_seq   = tree->c_lnL;
  new_lnL_seq   = UNLIKELY;
    
  u = Uni();
  mult = exp(K*(u-0.5));
    
  /* Multiply branch rates (or add to log of rates) */
  if(!Scale_Subtree_Rates(tree->n_root,mult,&n_nodes,tree))
    {
      RATES_Reset_Rates(tree);
      return;
    }

  if(n_nodes != 2*tree->n_otu-2) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);    

  /* Divide clock_r */
  if(tree->rates->clock_r_fixed == NO)
    {
      tree->rates->clock_r /= mult;
      if(tree->rates->clock_r < tree->rates->min_clock || tree->rates->clock_r > tree->rates->max_clock)
        {
          tree->rates->clock_r = init_clock;
          RATES_Reset_Rates(tree);
          return;
        }
    }
  
  if((tree->rates->model == GUINDON) ||
     (tree->rates->is_asynchronous == YES) ||
     (tree->rates->clock_r_fixed == YES))
    {
      RATES_Update_Cur_Bl(tree);
      if(tree->eval_alnL == YES) new_lnL_seq = Lk(NULL,tree);
    }
  else
    new_lnL_seq = cur_lnL_seq;
  
  if(tree->eval_rlnL == YES) new_lnL_rate = RATES_Lk_Rates(tree);

  ratio += (n_nodes-1)*log(mult);
  ratio += (new_lnL_rate - cur_lnL_rate);
  ratio += (new_lnL_seq - cur_lnL_seq);

  ratio = exp(ratio);
  alpha = MIN(1.,ratio);
  u = Uni();
  
  assert(isnan(u) == NO && isinf(fabs(u)) == NO);

  if(u > alpha)
    {
      tree->rates->clock_r = init_clock;
      RATES_Reset_Rates(tree);
      tree->rates->c_lnL_rates = cur_lnL_rate;
      if(tree->rates->model == GUINDON)
        {
          RATES_Update_Cur_Bl(tree);
          tree->c_lnL = cur_lnL_seq;
        }
    }
  else
    {
      tree->mcmc->acc_move[tree->mcmc->num_move_tree_rates]++;
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MCMC_Subtree_Rates(t_tree *tree)
{  phydbl K,mult,u,alpha,ratio;
  phydbl cur_lnL_rate,new_lnL_rate;
  phydbl cur_lnL_seq,new_lnL_seq;
  int target;
  int n_nodes;
  
  if(tree->rates->model == STRICTCLOCK) return;

  RATES_Record_Rates(tree);
  Record_Br_Len(tree);

  K             = tree->mcmc->tune_move[tree->mcmc->num_move_subtree_rates];
  cur_lnL_rate  = tree->rates->c_lnL_rates;
  new_lnL_rate  = cur_lnL_rate;
  cur_lnL_seq   = tree->c_lnL;
  new_lnL_seq   = cur_lnL_seq;
  ratio         = 0.0;

  u = Uni();
  mult = exp(K*(u-0.5));
  /* mult = Rgamma(1./K,K); */

  target = Rand_Int(tree->n_otu,2*tree->n_otu-3);

  /* Multiply branch rates */
  if(!Scale_Subtree_Rates(tree->a_nodes[target],mult,&n_nodes,tree))
    {
      RATES_Reset_Rates(tree);
      tree->mcmc->run_move[tree->mcmc->num_move_subtree_rates]++;
      return;
    }

  if(tree->eval_rlnL == YES) new_lnL_rate = RATES_Lk_Rates(tree);
  if(tree->eval_alnL == YES) new_lnL_seq  = Lk(NULL,tree);

  /* Proposal ratio: 2n-2=> number of multiplications, 1=>number of divisions */
  ratio += (+n_nodes)*log(mult);
  /* ratio += (n_nodes-2)*log(mult) + log(Dgamma(1./mult,1./K,K)/Dgamma(mult,1./K,K)); */
  

  ratio += (new_lnL_rate - cur_lnL_rate);
  ratio += (new_lnL_seq - cur_lnL_seq);


  ratio = exp(ratio);
  alpha = MIN(1.,ratio);
  u = Uni();
  
  assert(isnan(u) == NO && isinf(fabs(u)) == NO);

  if(u > alpha)
    {
      RATES_Reset_Rates(tree);
      Restore_Br_Len(tree);
      tree->rates->c_lnL_rates = cur_lnL_rate;
      tree->c_lnL = cur_lnL_seq;
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
  phydbl cur_lnL_seq,new_lnL_seq;
  phydbl cur_lnL_rate,new_lnL_rate;

  if(tree->rates->model == STRICTCLOCK) return;

  RATES_Record_Times(tree);
  RATES_Record_Rates(tree);
  Record_Br_Len(tree);

  K             = 3.;
  cur_lnL_seq  = tree->c_lnL;
  new_lnL_seq  = cur_lnL_seq;
  cur_lnL_rate  = tree->rates->c_lnL_rates;
  new_lnL_rate  = cur_lnL_rate;
  ratio         = 0.0;

  u = Uni();
  /* mult = exp(K*(u-0.5)); */
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

  RATES_Update_Cur_Bl(tree);
  if(tree->eval_alnL == YES) new_lnL_seq = Lk(NULL,tree);
  if(tree->eval_rlnL == YES) new_lnL_rate = RATES_Lk_Rates(tree);

  ratio += (-(tree->n_otu-1.)-2.)*log(mult);
  ratio += (new_lnL_rate - cur_lnL_rate);
  ratio += (new_lnL_seq - cur_lnL_seq);

  ratio = exp(ratio);
  alpha = MIN(1.,ratio);
  u = Uni();

  assert(isnan(u) == NO && isinf(fabs(u)) == NO);

  if(u > alpha)
    {
      RATES_Reset_Times(tree);
      RATES_Reset_Rates(tree);
      Restore_Br_Len(tree);
      tree->c_lnL = cur_lnL_seq;
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
  phydbl cur_lnL_seq,new_lnL_seq;
  
  RATES_Record_Rates(tree);
  Record_Br_Len(tree);

  K            = tree->mcmc->tune_move[tree->mcmc->num_move_updown_nu_cr];
  cur_lnL_seq  = tree->c_lnL;
  new_lnL_seq  = tree->c_lnL;
  cur_lnL_rate = tree->rates->c_lnL_rates;
  new_lnL_rate = tree->rates->c_lnL_rates;
  
  u = Uni();
  mult = exp(K*(u-0.5));

  /* Multiply branch rates */
  /* if(!Scale_Subtree_Rates(tree->n_root,mult,&n_nodes,tree)) */
  /*   { */
  /*     RATES_Reset_Rates(tree); */
  /*     Restore_Br_Len(tree); */
  /*     tree->mcmc->run_move[tree->mcmc->num_move_updown_nu_cr]++; */
  /*     return; */
  /*   } */


  if(tree->rates->clock_r_fixed == NO)
    {
      tree->rates->clock_r /= mult;
      if(tree->rates->clock_r < tree->rates->min_clock || tree->rates->clock_r > tree->rates->max_clock)
        {
          tree->rates->clock_r *= mult;
          RATES_Reset_Rates(tree);
          Restore_Br_Len(tree);
          tree->mcmc->run_move[tree->mcmc->num_move_updown_nu_cr]++;
          return;
        }
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
  
  RATES_Update_Cur_Bl(tree);
  if(tree->eval_alnL == YES) new_lnL_seq = Lk(NULL,tree);
  if(tree->eval_rlnL == YES) new_lnL_rate = RATES_Lk_Rates(tree);

  ratio = 0.0;
  /* Proposal ratio: 2n-2=> number of multiplications, 1=>number of divisions */
  /* ratio += n_nodes*log(mult); /\* (1-1)*log(mult); *\/ */
  ratio += 0.0*log(mult); /* (1-1)*log(mult); */
  /* Prior density ratio */
  ratio += (new_lnL_rate - cur_lnL_rate);
  /* Likelihood density ratio */
  ratio += (new_lnL_seq - cur_lnL_seq);

  ratio = exp(ratio);
  alpha = MIN(1.,ratio);
  u = Uni();
  
  assert(isnan(u) == NO && isinf(fabs(u)) == NO);

  if(u > alpha)
    {
      if(tree->rates->clock_r_fixed == NO) tree->rates->clock_r *= mult;
      tree->rates->nu /= mult;
      RATES_Reset_Rates(tree);
      Restore_Br_Len(tree);
      tree->rates->c_lnL_rates = cur_lnL_rate;
      tree->c_lnL        = cur_lnL_seq;
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
  phydbl cur_lnL_seq,new_lnL_seq;
  phydbl cur_lnL_rate,new_lnL_rate;
  phydbl cur_lnL_time,new_lnL_time;
  phydbl floor;
  int n_nodes;
  
  /*! Check that sequences are isochronous. */
  for(i=0;i<tree->n_otu-1;i++) if(!Are_Equal(tree->rates->nd_t[i+1],tree->rates->nd_t[i],1.E-10)) return; 
  
  RATES_Record_Times(tree);
  
  K            = tree->mcmc->tune_move[tree->mcmc->num_move_time_slice];
  cur_lnL_seq = tree->c_lnL;
  new_lnL_seq = UNLIKELY;
  cur_lnL_rate = tree->rates->c_lnL_rates;
  new_lnL_rate = UNLIKELY;
  cur_lnL_time = tree->rates->c_lnL_times;
  new_lnL_time = UNLIKELY;
  ratio        = 0.0;
  
  u = Uni();
  mult = exp(K*(u-0.5));

  floor = Uni()*tree->rates->nd_t[tree->n_root->num];

  Scale_Subtree_Height(tree->n_root,mult,floor,&n_nodes,tree);

  if(TIMES_Check_Node_Height_Ordering(tree) != NO)
    {
      RATES_Update_Cur_Bl(tree);
      
      new_lnL_time = TIMES_Lk_Times(NO,tree);

      if(new_lnL_time > UNLIKELY)
        {
          new_lnL_seq = Lk(NULL,tree);      
          new_lnL_rate = RATES_Lk_Rates(tree);
        }
      
      ratio += (phydbl)(n_nodes)*log(mult);

      /* Likelihood ratio */
      ratio += (new_lnL_seq - cur_lnL_seq);
      
      /* Prior ratio */
      ratio += (new_lnL_rate - cur_lnL_rate);
      ratio += (new_lnL_time - cur_lnL_time);
      
      /* printf("\n. seq: %f<-%f times: %f<-%f rates: %f<-%f acc.: [%d/%d] k: %f t: %f x: %f", */
      /*        new_lnL_seq,cur_lnL_seq, */
      /*        new_lnL_time,cur_lnL_time, */
      /*        new_lnL_rate,cur_lnL_rate, */
      /*        tree->mcmc->acc_move[tree->mcmc->num_move_time_slice], */
      /*        tree->mcmc->run_move[tree->mcmc->num_move_time_slice], */
      /*        K,floor,mult); */
      
      ratio = exp(ratio);
    }
  else
    {
      ratio = 0.0;
    }

  alpha = MIN(1.,ratio);
  u = Uni();
  
  assert(isnan(u) == NO && isinf(fabs(u)) == NO);

  if(u > alpha)
    {
      RATES_Reset_Times(tree);
      RATES_Update_Cur_Bl(tree);
      tree->rates->c_lnL_times = TIMES_Lk_Times(NO,tree); // Required in order to set t_prior_min/max to their original values
      tree->c_lnL              = cur_lnL_seq;
      tree->rates->c_lnL_rates = cur_lnL_rate;

      if(Are_Equal(tree->rates->c_lnL_times,cur_lnL_time,1.E-3) == NO)
        {
          PhyML_Fprintf(stderr,"\n. new_glnL: %f cur_glnL: %f",tree->rates->c_lnL_times,cur_lnL_time);
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

      for(i=0;i<tree->mcmc->n_moves;i++)
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

	  /* for(i=0;i<mcmc->n_moves;i++) */
	  /*   { */
	  /*     strcpy(s,"Acc."); */
	  /*     PhyML_Fprintf(fp,"%s%d\t",strcat(s,mcmc->move_name[i]),i); */
	  /*   } */

	  /* for(i=0;i<mcmc->n_moves;i++) */
	  /*   { */
	  /*     strcpy(s,"Tune."); */
	  /*     PhyML_Fprintf(fp,"%s%d\t",strcat(s,mcmc->move_name[i]),i); */
	  /*   } */

	  /* for(i=0;i<mcmc->n_moves;i++) */
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
		  for(i=0;i<tree->mod->ras->n_catg;i++) PhyML_Fprintf(fp,"p%d\t",i);
		  for(i=0;i<tree->mod->ras->n_catg;i++) PhyML_Fprintf(fp,"r%d\t",i);
		}
	    }


	  if(tree->mod->m4mod->n_h > 1 && tree->mod->use_m4mod == YES)
	    {
	      for(i=0;i<tree->mod->m4mod->n_h;i++) PhyML_Fprintf(fp,"cov_p%d\t",i);
	      for(i=0;i<tree->mod->m4mod->n_h;i++) PhyML_Fprintf(fp,"cov_r%d\t",i);
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
	  for(i=0;i<tree->n_otu-1;i++) PhyML_Fprintf(mcmc->out_fp_trees,"\t%3d\t%s,\n",tree->a_nodes[i]->num+1,tree->a_nodes[i]->name);
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
      /* for(i=0;i<tree->mcmc->n_moves;i++) PhyML_Fprintf(fp,"%f\t",tree->mcmc->acc_rate[i]); */
      /* for(i=0;i<tree->mcmc->n_moves;i++) PhyML_Fprintf(fp,"%f\t",(phydbl)(tree->mcmc->tune_move[i])); */
/*       for(i=0;i<tree->mcmc->n_moves;i++) PhyML_Fprintf(fp,"%d\t",(int)(tree->mcmc->run_move[i])); */

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
	      for(i=0;i<tree->mod->ras->n_catg;i++) PhyML_Fprintf(fp,"%G\t",tree->mod->ras->gamma_r_proba->v[i]);
	      for(i=0;i<tree->mod->ras->n_catg;i++) PhyML_Fprintf(fp,"%G\t",tree->mod->ras->gamma_rr->v[i]);
	      /* for(i=0;i<tree->mod->ras->n_catg;i++) PhyML_Fprintf(fp,"%G\t",tree->mod->ras->gamma_r_proba_unscaled[i]); */
	      /* for(i=0;i<tree->mod->ras->n_catg;i++) PhyML_Fprintf(fp,"%G\t",tree->mod->ras->gamma_rr_unscaled[i]); */
	    }
	}

      if(tree->mod->m4mod->n_h > 1 && tree->mod->use_m4mod == YES)
	{
	  for(i=0;i<tree->mod->m4mod->n_h;i++) PhyML_Fprintf(fp,"%G\t",tree->mod->m4mod->h_fq[i]);
	  for(i=0;i<tree->mod->m4mod->n_h;i++) PhyML_Fprintf(fp,"%G\t",tree->mod->m4mod->multipl[i]);
	  PhyML_Fprintf(fp,"%G\t",tree->mod->m4mod->delta);
	}


      char *format = (char *)mCalloc(100,sizeof(char));
      sprintf(format,"%%.%df\t",tree->mcmc->nd_t_digits);
      for(i=tree->n_otu;i<2*tree->n_otu-1;i++) PhyML_Fprintf(fp,format,tree->rates->nd_t[i]);
      Free(format);

      /* for(i=0;i<2*tree->n_otu-1;i++) PhyML_Fprintf(fp,"%.4f\t",log(tree->rates->nd_r[i])); */
      
      // Average rate along edges: length divided by elapsed time
      For(i,2*tree->n_otu-2)
	PhyML_Fprintf(fp,"%.4f\t",
		      tree->rates->cur_l[i]/(tree->rates->nd_t[tree->a_nodes[i]->num] - tree->rates->nd_t[tree->a_nodes[i]->anc->num]));

      /* fp_pred = fopen("predict.txt","a");       */
      /* for(i=0;i<2*tree->n_otu-2;i++)  */
      /* 	PhyML_Fprintf(fp_pred,"B%d\t%12f\t%12f\t%4d\n",i,exp(tree->rates->br_r[i]),tree->rates->nd_t[i],tree->rates->has_survived[i]); */
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
/*       For(i,2*tree->n_otu-3) PhyML_Fprintf(fp,"%f\t",exp(tree->a_edges[i]->l->v)); */


      /* RATES_Update_Cur_Bl(tree); */
      /* For(i,2*tree->n_otu-3) PhyML_Fprintf(fp,"%f\t",tree->a_edges[i]->l->v); */


      For(i,2*tree->n_otu-2) tree->rates->mean_r[i] = exp(tree->rates->br_r[i]);
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
      /* 	      for(j=0;j<3;j++) */
      /* 		{ */
      /* 		  if(d->v[j] != d->anc && d->b[j] != tree->e_root) */
      /* 		    { */
      /* 		      if(!v1) v1 = d->v[j]; */
      /* 		      else    v2 = d->v[j]; */
      /* 		    } */
      /* 		} */

      /* 	      for(j=0;j<3;j++) */
      /* 		{ */
      /* 		  if(v1->v[j] && v1->v[j] == d) */
      /* 		    { */
      /* 		      n1 = v1->bip_size[j]; */
      /* 		      /\* r1 = RATES_Get_Mean_Rate_In_Subtree(v1,tree); *\/ */
      /* 		      r1 = exp(tree->rates->br_r[v1->num]); */
      /* 		      break; */
      /* 		    } */
      /* 		} */


      /* 	      for(j=0;j<3;j++) */
      /* 		{ */
      /* 		  if(v2->v[j] && v2->v[j] == d) */
      /* 		    { */
      /* 		      n2 = v2->bip_size[j]; */
      /* 		      /\* r2 = RATES_Get_Mean_Rate_In_Subtree(v2,tree); *\/ */
      /* 		      r2 = exp(tree->rates->br_r[v2->num]); */
      /* 		      break; */
      /* 		    } */
      /* 		} */

      /* 	      fprintf(fp,"\n%4d %4d %15f %15f",n1,n2,r1,r2); */
      /* 	    } */
      /* 	} */

      /* fclose(fp); */
      
      // TREES
      TIMES_Time_To_Bl(tree);
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
		PhyML_Fprintf(stderr,"\n. The value entered must be an integer greater than 0.\n");
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
  cpy->move_idx           = ori->move_idx        ;
  
  for(i=0;i<cpy->n_moves;i++) 
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
      for(i=0;i<tree->mod->ras->n_catg-1;i++) tree->mod->ras->gamma_r_proba_unscaled->v[i] = Uni();
      tree->mod->ras->gamma_r_proba_unscaled->v[tree->mod->ras->n_catg-1] = 1.;
      for(i=0;i<tree->mod->ras->n_catg-1;i++) tree->mod->ras->gamma_rr_unscaled->v[i] = (phydbl)i+0.1; /* Do not randomize those as their ordering matter */
      tree->mod->ras->gamma_rr_unscaled->v[tree->mod->ras->n_catg-1] = (phydbl)tree->mod->ras->n_catg;
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

  for(i=0;i<tree->mod->m4mod->n_h;i++)
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

  for(i=0;i<2*tree->n_otu-2;++i)
    if(tree->a_nodes[i] != tree->n_root)
      {
        tree->rates->nd_r[i] = Rnorm_Trunc(mean_r,SQRT(var_r),min_r,max_r,&err);
        if(err == YES) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
      }
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
  phydbl mean_r, sd_r;
  phydbl min_r, max_r;
  int err;

  mean_r = 1.0;
  sd_r   = 1.0;
  min_r  = tree->rates->min_rate;
  max_r  = tree->rates->max_rate;
  err    = -1;
    
  tree->rates->br_r[d->num] = Rnorm_Trunc(mean_r,sd_r,min_r,max_r,&err);
  if(err == YES)
    {
      PhyML_Fprintf(stderr,"\n. mean_r: %f sd_r: %f min_r: %f max_r: %f",
                    mean_r,
                    sd_r,
                    min_r,
                    max_r);      
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
    }
  if(d->tax) return;
  else
    {
      for(i=0;i<3;++i)
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
  if(tree->rates->clock_r_fixed == NO) tree->rates->clock_r = u * (1.0 - tree->rates->min_clock) + tree->rates->min_clock;
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
      PhyML_Fprintf(stderr,"\n. min_dt = %f",min_dt);
      PhyML_Fprintf(stderr,"\n. min->t=%f min->anc->t=%f",tree->rates->nd_t[min_node],tree->rates->nd_t[tree->a_nodes[min_node]->anc->num]);
      PhyML_Fprintf(stderr,"\n. d up=%f down=%f",tree->rates->t_prior_min[min_node],tree->rates->t_prior_max[min_node]);
      PhyML_Fprintf(stderr,"\n. a up=%f down=%f",tree->rates->t_prior_min[tree->a_nodes[min_node]->anc->num],tree->rates->t_prior_max[tree->a_nodes[min_node]->anc->num]);
      PhyML_Fprintf(stderr,"\n. up=%f down=%f",tree->rates->t_prior_min[min_node],tree->rates->t_floor[tree->a_nodes[min_node]->anc->num]);
      PhyML_Fprintf(stderr,"\n. min_node = %d",min_node);
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


      for(i=0;i<3;i++)
	{
	  if((d->v[i] != a) && (d->b[i] != tree->e_root))
	    {
	      MCMC_Randomize_Node_Times_Bottom_Up(d,d->v[i],tree);
	    }
	}

      v1 = v2 = NULL;
      for(i=0;i<3;i++)
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

      for(i=0;i<3;i++)
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


  lag = 100;

  eps = 1.E-6;

  for(i=0;i<mcmc->n_moves;i++)
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
  if(mcmc->adjust_tuning[move] == YES)
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
      else if(!strcmp(mcmc->move_name[move],"phyrex_indel_disk"))
	{
	  rate_inf = 0.234;
	  rate_sup = 0.234;
	}
      else if(!strcmp(mcmc->move_name[move],"phyrex_indel_hit"))
	{
	  rate_inf = 0.234;
	  rate_sup = 0.234;
	}
      else if(!strcmp(mcmc->move_name[move],"phyrex_scale_times"))
	{
	  rate_inf = 0.234;
	  rate_sup = 0.234;
	}
      else if(!strcmp(mcmc->move_name[move],"phyrex_ldsk_given_disk"))
	{
	  rate_inf = 0.234;
	  rate_sup = 0.234;
	}
      else if(!strcmp(mcmc->move_name[move],"phyrex_disk_given_ldsk"))
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
  phydbl new_lnL_seq, cur_lnL_seq;
  phydbl ratio, alpha;
  phydbl new_l, cur_l;
  phydbl K,mult;


  cur_l       = b->l->v;
  cur_lnL_seq = tree->c_lnL;
  new_lnL_seq = tree->c_lnL;
  K           = 0.1;
  
  u = Uni();
  mult = exp(K*(u-0.5));
  /* mult = u*(K-1./K)+1./K; */
  new_l = cur_l * mult;

  if(new_l < tree->mod->l_min || new_l > tree->mod->l_max) return;

  b->l->v = new_l;
  if(tree->eval_alnL == YES) new_lnL_seq = Lk(b,tree);

  ratio =
    (new_lnL_seq - cur_lnL_seq) +
    (log(mult));


  ratio = exp(ratio);
  alpha = MIN(1.,ratio);
  
  u = Uni();
  
  assert(isnan(u) == NO && isinf(fabs(u)) == NO);

  if(u > alpha) /* Reject */
    {
      b->l->v = cur_l;
      Update_PMat_At_Given_Edge(b,tree);
      tree->c_lnL = cur_lnL_seq;
    }

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MCMC_Scale_Br_Lens(t_tree *tree)
{
  phydbl u;
  phydbl new_lnL_seq, cur_lnL_seq;
  phydbl ratio, alpha;
  phydbl K,mult;
  int i;

  Record_Br_Len(tree);

  cur_lnL_seq = tree->c_lnL;
  new_lnL_seq = tree->c_lnL;
  K            = 1.2;
  
  u = Uni();
  mult = u*(K-1./K)+1./K;

  for(i=0;i<2*tree->n_otu-3;++i) 
    {
      tree->a_edges[i]->l->v *= mult;
      if(tree->a_edges[i]->l->v < tree->mod->l_min || 
	 tree->a_edges[i]->l->v > tree->mod->l_max) return;
    }

  Set_Both_Sides(NO,tree);
  if(tree->eval_alnL == YES) new_lnL_seq = Lk(NULL,tree);

  ratio =
    (new_lnL_seq - cur_lnL_seq) +
    (2*tree->n_otu-5) * (log(mult));

  ratio = exp(ratio);
  alpha = MIN(1.,ratio);
  
  u = Uni();
  
  assert(isnan(u) == NO && isinf(fabs(u)) == NO);

  if(u > alpha) /* Reject */
    {
      Restore_Br_Len(tree);
      tree->c_lnL = cur_lnL_seq;
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
      for(i=0;i<3;i++) 
	if(d->v[i] != a)
	  {
	    Update_Partial_Lk(tree,d->b[i],d);
	    MCMC_Br_Lens_Pre(d,d->v[i],d->b[i],tree);
	  }
      Update_Partial_Lk(tree,b,d);
    }  
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// Potentially one kappa parameter for each tree in a mixture -> update each of them
void MCMC_Kappa(t_tree *mixt_tree)
{
  t_tree *tree;
  phydbl cur_kappa,new_kappa;
  phydbl u,alpha,ratio;
  phydbl min_kappa,max_kappa;
  phydbl K;
  phydbl cur_lnL_seq, new_lnL_seq;

  if(mixt_tree->eval_alnL == NO) return;
  
  Set_Update_Eigen(YES,mixt_tree->mod);

  tree = mixt_tree;
    
  do
    {
      if(tree->is_mixt_tree == YES) tree = tree->next;
      
      if(!(tree->mod->whichmodel == HKY85 || tree->mod->whichmodel == K80 || tree->mod->whichmodel == TN93)) tree = tree->next;

      if(tree == NULL) return;
      
      cur_kappa     = -1.0;
      new_kappa     = -1.0;
      ratio         =  0.0;
      
      K = mixt_tree->mcmc->tune_move[mixt_tree->mcmc->num_move_kappa];

      cur_lnL_seq = mixt_tree->c_lnL;
      new_lnL_seq = mixt_tree->c_lnL;
      cur_kappa   = tree->mod->kappa->v;
      
      min_kappa = 0.1;
      max_kappa = 10.;

      MCMC_Make_Move(&cur_kappa,&new_kappa,min_kappa,max_kappa,&ratio,K,mixt_tree->mcmc->move_type[mixt_tree->mcmc->num_move_kappa]);
      
      if(new_kappa < max_kappa && new_kappa > min_kappa) 
        {
          tree->mod->kappa->v = new_kappa;
          
          if(mixt_tree->eval_alnL == YES) new_lnL_seq = Lk(NULL,mixt_tree);
          
          ratio += (new_lnL_seq - cur_lnL_seq);
                    
          ratio = exp(ratio);
          alpha = MIN(1.,ratio);
          
          u = Uni();

          assert(isnan(u) == NO && isinf(fabs(u)) == NO);

          if(u > alpha) /* Reject */
            {
              tree->mod->kappa->v = cur_kappa;
              mixt_tree->c_lnL    = cur_lnL_seq;
              if(!Set_Model_Parameters(mixt_tree->mod))
                {
                  PhyML_Fprintf(stderr,"\n. Problem in move %s",tree->mcmc->move_name[tree->mcmc->move_idx]);
                  Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
                }
            }
          else
            {
              mixt_tree->mcmc->acc_move[mixt_tree->mcmc->num_move_kappa]++;
            }
          
          /* PhyML_Printf("\n. MCMC cur_k: %f new_k: %f cur: %f new: %f k: %f -> %f",cur_kappa,new_kappa,cur_lnL_seq,new_lnL_seq,tree->mod->kappa->v,mixt_tree->c_lnL); */

        }
      mixt_tree->mcmc->run_move[mixt_tree->mcmc->num_move_kappa]++;
      
      tree = tree->next;
    }
  while(tree != NULL);
  
  Set_Update_Eigen(NO,mixt_tree->mod);
  
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// Potentially one set of relative rate parameters for each tree in a mixture -> update each of them

void MCMC_RR(t_tree *mixt_tree)
{
  t_tree *tree;
  phydbl cur_rr,new_rr;
  phydbl u,alpha,ratio;
  int n_r_mat,*permut;  
  phydbl K;
  phydbl cur_lnL_seq, new_lnL_seq;
  t_rmat **r_mat;
  int i;

  if(mixt_tree->eval_alnL == NO) return;

  Set_Update_Eigen(YES,mixt_tree->mod);

  tree    = mixt_tree;
  n_r_mat = 0;
  r_mat   = NULL;
  permut  = NULL;
  
  do
    {
      if(tree->is_mixt_tree == YES) tree = tree->next;

      if(!(tree->mod->whichmodel == GTR || tree->mod->whichmodel == CUSTOM)) tree = tree->next;

      if(tree == NULL) return;
      
      for(i=0;i<n_r_mat;i++) if(tree->mod->r_mat == r_mat[i]) break;

      if(i == n_r_mat &&
         (tree->mod->whichmodel == GTR || tree->mod->whichmodel == CUSTOM) &&
         tree->mod->r_mat->n_diff_rr > 1)
        {          
          permut = Permutate(tree->mod->r_mat->n_diff_rr);

          for(i=0;i<tree->mod->r_mat->n_diff_rr;++i)
            {
              cur_rr      = -1.0;
              new_rr      = -1.0;
              ratio       =  0.0;          
              K           = 1.0;
              cur_lnL_seq = mixt_tree->c_lnL;
              new_lnL_seq = UNLIKELY;
      
              mixt_tree->mcmc->run_move[mixt_tree->mcmc->num_move_rr]++;

              cur_rr = tree->mod->r_mat->rr_val->v[permut[i]];
      
              MCMC_Make_Move(&cur_rr,&new_rr,UNSCALED_RR_MIN,UNSCALED_RR_MAX,&ratio,K,mixt_tree->mcmc->move_type[mixt_tree->mcmc->num_move_rr]);
              
              if(new_rr < UNSCALED_RR_MAX && new_rr > UNSCALED_RR_MIN)
                {
                  tree->mod->r_mat->rr_val->v[permut[i]] = new_rr;

                  new_lnL_seq = Lk(NULL,mixt_tree);
          
                  ratio += (new_lnL_seq - cur_lnL_seq);
                  ratio += log(new_rr) - log(cur_rr); /* Because we are updating log(rr) instead of rr */
                  
          
                  ratio = exp(ratio);
                  alpha = MIN(1.,ratio);
          
                  u = Uni();

                  assert(isnan(u) == NO && isinf(fabs(u)) == NO);

                  if(u > alpha) /* Reject */
                    {
                      tree->mod->r_mat->rr_val->v[permut[i]] = cur_rr;
                      mixt_tree->c_lnL = cur_lnL_seq;
                      if(!Set_Model_Parameters(mixt_tree->mod))
                        {
                          PhyML_Fprintf(stderr,"\n. Problem in move %s",tree->mcmc->move_name[tree->mcmc->move_idx]);
                          Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
                        }
                    }
                  else
                    {
                      mixt_tree->mcmc->acc_move[mixt_tree->mcmc->num_move_rr]++;
                    }
                  
                  /* PhyML_Printf("\n. MCMC cur_rr: %f new_rr: %f cur: %f new: %f",cur_rr,new_rr,cur_lnL_seq,new_lnL_seq); */
                }
            }
          Free(permut);
        }
      tree = tree->next;
    }
  while(tree != NULL);
  
  Set_Update_Eigen(NO,mixt_tree->mod);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MCMC_Rate_Across_Sites(t_tree *tree)
{
  int i;

  if(tree->mod->ras->n_catg == 1) return;

  for(i=0;i<tree->mod->ras->n_catg;i++)
    {
      if(tree->mod->ras->free_mixt_rates == YES) MCMC_Free_Mixt_Rate(tree);
      else                                       MCMC_Alpha(tree);
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

void MCMC_Free_Mixt_Rate(t_tree *mixt_tree)
{
  t_tree *tree;
  phydbl *u,*v;
  phydbl *f;
  phydbl cur_val;
  phydbl cur_A, new_A;
  phydbl low_bound,up_bound;
  unsigned int k,i,idx; 
  phydbl hr,z;
  int n_moves;
  phydbl cur_lnL_seq, new_lnL_seq;
  phydbl ratio,alpha;

  tree = mixt_tree;

  do
    {
      
      tree->mod->ras->sort_rate_classes = YES;
      tree->mod->ras->normalise_rr      = YES;

      cur_lnL_seq = mixt_tree->c_lnL;
      new_lnL_seq = mixt_tree->c_lnL;

      k = tree->mod->ras->n_catg;

      f = tree->mod->ras->gamma_r_proba->v;

      u = tree->mod->ras->gamma_r_proba_unscaled->v;
      v = tree->mod->ras->gamma_rr_unscaled->v;


      // Setting the values of the unscaled rr freq as is done below
      // does not lead to changing the likelihood.
      u[0] = f[0];
      for(i=1;i<k;++i) u[i] = u[i-1] + f[i];      

            
      n_moves = 0;
      do
        {
          n_moves++;
          
          z = Uni();
          
          if(z < 0.5)
            {
              // Update frequencies

              cur_A = 0.0;
              for(i=0;i<k;++i) cur_A += f[i]*v[i];
              
              hr = +2.*log(cur_A);
              for(i=0;i<k-1;++i) hr -= log(cur_A - f[i]*v[i]);
              
              idx = Rand_Int(0,k-2); 
              
              // Proposal is uniform. Determine upper and lower bounds.
              z = Uni();
              low_bound = (idx==0)?(.0):(u[idx-1]);
              up_bound  = u[idx+1];
              cur_val = u[idx];
              u[idx] = low_bound + z*(up_bound - low_bound);

              Update_RAS(tree->mod);
              for(i=0;i<k;++i)
                {
                  if(u[i] < GAMMA_RR_UNSCALED_MIN || u[i] > GAMMA_RR_UNSCALED_MAX) break;
                  if(v[i] < GAMMA_R_PROBA_UNSCALED_MIN || v[i] > GAMMA_R_PROBA_UNSCALED_MAX) break;
                }
              if(i != k) new_lnL_seq = -INFINITY;
              else new_lnL_seq = Lk(NULL,mixt_tree);
              
              new_A = 0.0;
              for(i=0;i<k;++i) new_A += f[i]*v[i];
              
              hr -= 2.*log(new_A);
              for(i=0;i<k-1;++i) hr += log(new_A - f[i]*v[i]);
                            
              /* Metropolis-Hastings step */
              ratio = 0.;
              ratio += (new_lnL_seq - cur_lnL_seq);
              ratio += hr;
              ratio = exp(ratio);
              alpha = MIN(1.,ratio);
              
              /* printf("\nf class=%d new_val=%f cur_val=%f cur: %f -> new: %f",idx,u[idx],cur_val,cur_lnL_seq,new_lnL_seq); */
              
              z = Uni();
              if(z > alpha) /* Reject */
                {
                  u[idx] = cur_val;

                  if(!Set_Model_Parameters(mixt_tree->mod))
                    {
                      PhyML_Fprintf(stderr,"\n. Problem in move %s",tree->mcmc->move_name[tree->mcmc->move_idx]);
                      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
                    }

                  /* phydbl cur_lk = cur_lnL_seq; */
                  /* phydbl new_lk = Lk(NULL,mixt_tree); */
                  /* if(Are_Equal(cur_lk,new_lk,1.E-5) == NO) */
                  /*   { */
                  /*     PhyML_Printf("\n. new: %f cur: %f",new_lk,cur_lk); */
                  /*     assert(FALSE); */
                  /*   } */

                  mixt_tree->c_lnL = cur_lnL_seq;
                }
              else /* Accept */
                {
                  cur_lnL_seq = new_lnL_seq;
                }
              
            }
          else
            {
              // Update rates
                            
              cur_A = 0.0;
              for(i=0;i<k;++i) cur_A += f[i]*v[i];
              
              hr = +2.*log(cur_A);
              for(i=0;i<k-1;++i) hr -= log(cur_A - f[i]*v[i]);
              
              idx = Rand_Int(0,k-2);
              
              // Proposal is uniform. Determine upper and lower bounds.
              z = Uni();
              low_bound = (idx==0)?(.0):(v[idx-1]);
              up_bound = v[idx+1];
              cur_val = v[idx];
              v[idx] = low_bound + z*(up_bound - low_bound);
              
              Update_RAS(tree->mod);
              for(i=0;i<k;++i)
                {
                  if(u[i] < GAMMA_RR_UNSCALED_MIN || u[i] > GAMMA_RR_UNSCALED_MAX) break;
                  if(v[i] < GAMMA_R_PROBA_UNSCALED_MIN || v[i] > GAMMA_R_PROBA_UNSCALED_MAX) break;
                }
              if(i != k) new_lnL_seq = -INFINITY;
              else new_lnL_seq = Lk(NULL,mixt_tree);

              new_A = 0.0;
              for(i=0;i<k;++i) new_A += f[i]*v[i];
              
              hr -= 2.*log(new_A);
              for(i=0;i<k-1;++i) hr += log(new_A - f[i]*v[i]);

              
              /* Metropolis-Hastings step */
              ratio = 0.;
              ratio += (new_lnL_seq - cur_lnL_seq);
              ratio += hr;
              ratio = exp(ratio);
              alpha = MIN(1.,ratio);
              
              /* printf("\nr class=%d new_val=%f cur_val=%f cur: %f -> new: %f",idx,v[idx],cur_val,cur_lnL_seq,new_lnL_seq); */
              

              z = Uni();
              if(z > alpha) /* Reject */
                {
                  v[idx] = cur_val;
                  
                  if(!Set_Model_Parameters(mixt_tree->mod))
                    {
                      PhyML_Fprintf(stderr,"\n. Problem in move %s",tree->mcmc->move_name[tree->mcmc->move_idx]);
                      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
                    }


                  /* phydbl cur_lk = cur_lnL_seq; */
                  /* phydbl new_lk = Lk(NULL,mixt_tree); */
                  /* if(Are_Equal(cur_lk,new_lk,1.E-5) == NO) */
                  /*   { */
                  /*     PhyML_Printf("\n. new: %f cur: %f",new_lk,cur_lk); */
                  /*     assert(FALSE); */
                  /*   } */

                  mixt_tree->c_lnL = cur_lnL_seq;
                }
              else /* Accept */
                {
                  cur_lnL_seq = new_lnL_seq;
                }
              
            }
        }
      while(n_moves != k);
      
      tree = tree->next_mixt;
    }
  while(tree);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MCMC_Covarion_Rates(t_tree *tree)
{
  int i, class;
  phydbl u;
  phydbl min,max;

  if(tree->mod->use_m4mod == NO) return;

  Set_Update_Eigen(YES,tree->mod);

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

  Set_Update_Eigen(NO,tree->mod);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Covarion_Switch(t_tree *tree)
{
  if(tree->mod->use_m4mod == NO) return;

  Set_Update_Eigen(YES,tree->mod);
  MCMC_Single_Param_Generic(&(tree->mod->m4mod->delta),0.01,+100.,tree->mcmc->num_move_cov_switch,
			    NULL,&(tree->c_lnL),
			    NULL,Wrap_Lk,tree->mcmc->move_type[tree->mcmc->num_move_cov_switch],NO,NULL,tree,NULL); 
  Set_Update_Eigen(NO,tree->mod);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MCMC_Birth_Rate(t_tree *tree)
{
  /* MCMC_Single_Param_Generic(&(tree->rates->birth_rate), */
  /*       		    tree->rates->birth_rate_min, */
  /*       		    tree->rates->birth_rate_max, */
  /*       		    tree->mcmc->num_move_birth_rate, */
  /*       		    &(tree->rates->c_lnL_times),NULL, */
  /*       		    Wrap_Lk_Times,NULL,tree->mcmc->move_type[tree->mcmc->num_move_birth_rate],NO,NULL,tree,NULL); */
 
  phydbl cur_birth_rate,new_birth_rate;
  phydbl cur_lnL_time,new_lnL_time;
  phydbl cur_lnL_time_ghost,new_lnL_time_ghost;
  /* phydbl cur_lnL_time_pivot,new_lnL_time_pivot; */
  phydbl u,alpha,ratio;
  phydbl birth_rate_min,birth_rate_max;
  phydbl K;
  int i,n_mcmc_steps,move;

  cur_birth_rate = -1.0;
  new_birth_rate = -1.0;
  ratio          =  0.0;
  n_mcmc_steps   =  tree->mcmc->run_move[tree->mcmc->num_move_birth_rate] == 1 ? 1000 : 100;
  move           = -1;

  K = tree->mcmc->tune_move[tree->mcmc->num_move_birth_rate];

  cur_lnL_time = tree->rates->c_lnL_times;
  new_lnL_time = tree->rates->c_lnL_times;

  cur_lnL_time_ghost = UNLIKELY;
  new_lnL_time_ghost = UNLIKELY;

  /* cur_lnL_time_pivot = UNLIKELY; */
  /* new_lnL_time_pivot = UNLIKELY; */

  cur_birth_rate = tree->rates->birth_rate;

  birth_rate_min = MAX(tree->rates->birth_rate_min,tree->rates->death_rate);
  birth_rate_max = tree->rates->birth_rate_max;
  
  MCMC_Make_Move(&cur_birth_rate,&new_birth_rate,birth_rate_min,birth_rate_max,&ratio,K,tree->mcmc->move_type[tree->mcmc->num_move_birth_rate]);
  /* new_birth_rate = Uni()*(birth_rate_max - birth_rate_min) + birth_rate_min; */
  
  if(new_birth_rate < birth_rate_max && new_birth_rate > birth_rate_min && new_birth_rate > tree->rates->death_rate)
    {
      tree->rates->birth_rate = new_birth_rate;
      new_lnL_time = TIMES_Lk_Times(NO,tree);
      ratio += (new_lnL_time - cur_lnL_time);
      
      /* printf("\n.  b :%4d: %12G -> %12G ratio : %12G [real: %12G %12G]", */
      /*        tree->mcmc->run_move[tree->mcmc->num_move_birth_rate], */
      /*        cur_birth_rate, */
      /*        new_birth_rate, */
      /*        ratio, */
      /*        cur_lnL_time, */
      /*        new_lnL_time); */


      if(tree->mcmc->run_move[tree->mcmc->num_move_birth_rate] == 500)
        {
          tree->rates->birth_rate_pivot = tree->rates->birth_rate;
          tree->rates->death_rate_pivot = tree->rates->death_rate;
        }
      
      if(tree->mcmc->run_move[tree->mcmc->num_move_birth_rate] >= 0)
        {
          /* if(tree->mcmc->run_move[tree->mcmc->num_move_birth_rate] >= 500) */
          /*   { */
          /*     tree->extra_tree->rates->birth_rate = tree->rates->birth_rate_pivot; */
          /*     tree->extra_tree->rates->death_rate = tree->rates->death_rate_pivot; */
          /*     cur_lnL_time_pivot = TIMES_Lk_Times(NO,tree->extra_tree); */
          /*   } */
          
          /* tree->extra_tree->rates->birth_rate = cur_birth_rate; */
          /* tree->extra_tree->rates->death_rate = tree->rates->death_rate; */
          /* cur_lnL_time_ghost = TIMES_Lk_Times(NO,tree->extra_tree); */
          
          /* Copy_Tree(tree->extra_tree,tree->extra_tree->extra_tree); */
          /* RATES_Copy_Rate_Struct(tree->extra_tree->rates,tree->extra_tree->extra_tree->rates,tree->n_otu); */
          /* DATE_Assign_Primary_Calibration(tree->extra_tree->extra_tree); */
          


          tree->extra_tree->eval_alnL          = NO;
          tree->extra_tree->eval_rlnL          = NO;
          tree->extra_tree->eval_glnL          = YES;
          tree->extra_tree->rates->birth_rate  = new_birth_rate;
          tree->extra_tree->rates->c_lnL_rates = UNLIKELY;
          tree->extra_tree->c_lnL              = UNLIKELY;
          
          TIMES_Randomize_Tree_With_Time_Constraints(tree->extra_tree->rates->a_cal[0],tree->extra_tree);
          tree->extra_tree->rates->birth_rate = new_birth_rate;
          tree->extra_tree->rates->death_rate = tree->rates->death_rate;
          TIMES_Lk_Times(NO,tree->extra_tree);
          
          if(!(tree->extra_tree->rates->c_lnL_times > UNLIKELY))
            {
              PhyML_Fprintf(stderr,"\n. glnL=%f",tree->extra_tree->rates->c_lnL_times);
              PhyML_Fprintf(stderr,"\n. birth=%G death=%G [%G]",new_birth_rate,tree->rates->death_rate,tree->extra_tree->rates->death_rate);
              TIMES_Lk_Times(YES,tree->extra_tree);
              assert(FALSE);
            }
          
          i = 0;
          do
            {
              u = Uni();
              for(move=0;move<tree->mcmc->n_moves;move++) if(tree->mcmc->move_weight[move] > u-1.E-10) break;
              

              /* if(!(i%10)) PhyML_Printf("\n<< %5d Move '%20s' %12f",i,tree->mcmc->move_name[move],tree->extra_tree->rates->c_lnL_times); */
              
              if(!strcmp(tree->mcmc->move_name[move],"tree_height")) { MCMC_Tree_Height(tree->extra_tree); i++; }
              if(!strcmp(tree->mcmc->move_name[move],"times"))       { MCMC_Times_All(tree->extra_tree); i++; }
              if(!strcmp(tree->mcmc->move_name[move],"spr"))         { MCMC_Prune_Regraft(tree->extra_tree); i++; }
              if(!strcmp(tree->mcmc->move_name[move],"spr_local"))   { MCMC_Prune_Regraft_Local(tree->extra_tree); i++; }
              
              if(!(tree->extra_tree->rates->c_lnL_times > UNLIKELY))
                {
                  PhyML_Fprintf(stderr,"\n. move: %s",tree->mcmc->move_name[move]);
                  PhyML_Fprintf(stderr,"\n. glnL=%f",tree->extra_tree->rates->c_lnL_times);
                  TIMES_Lk_Times(YES,tree->extra_tree);
                  assert(FALSE);
                }

              /* PhyML_Printf("\n.> %4d %15f",i,tree->extra_tree->rates->c_lnL_times); */
            }
          while(i < n_mcmc_steps);
          
          tree->extra_tree->rates->birth_rate = cur_birth_rate;
          tree->extra_tree->rates->death_rate = tree->rates->death_rate;
          cur_lnL_time_ghost = TIMES_Lk_Times(NO,tree->extra_tree);

          tree->extra_tree->rates->birth_rate = new_birth_rate;
          tree->extra_tree->rates->death_rate = tree->rates->death_rate;
          new_lnL_time_ghost = TIMES_Lk_Times(NO,tree->extra_tree);

          ratio += (cur_lnL_time_ghost - new_lnL_time_ghost);

          /* if(tree->mcmc->run_move[tree->mcmc->num_move_birth_rate] >= 500) */
          /*   { */
          /*     tree->extra_tree->rates->birth_rate = tree->rates->birth_rate_pivot; */
          /*     tree->extra_tree->rates->death_rate = tree->rates->death_rate_pivot; */
          /*     new_lnL_time_pivot = TIMES_Lk_Times(NO,tree->extra_tree); */
          /*     ratio += (new_lnL_time_pivot - cur_lnL_time_pivot); */
          /*   } */
        }
      
      ratio = exp(ratio);
      alpha = MIN(1.,ratio);
      
      /* printf("\n.  b :%4d: %12G -> %12G ratio : %12G [ghost: %12G %12G -- real: %12G %12G]", */
      /*        tree->mcmc->run_move[tree->mcmc->num_move_birth_rate], */
      /*        cur_birth_rate, */
      /*        new_birth_rate, */
      /*        ratio, */
      /*        cur_lnL_time_ghost, */
      /*        new_lnL_time_ghost, */
      /*        cur_lnL_time, */
      /*        new_lnL_time); */
            
      u = Uni();

      assert(isnan(u) == NO && isinf(fabs(u)) == NO);

      if(u > alpha) /* Reject */
        {
          /* PhyML_Printf("  reject"); */
          tree->rates->birth_rate  = cur_birth_rate;
          tree->rates->c_lnL_times = cur_lnL_time;

          /* Copy_Tree(tree->extra_tree->extra_tree,tree->extra_tree); */
          /* RATES_Copy_Rate_Struct(tree->extra_tree->extra_tree->rates,tree->extra_tree->rates,tree->n_otu); */
          /* DATE_Assign_Primary_Calibration(tree->extra_tree); */
        }
      else
        {
          /* PhyML_Printf("  accept"); */
          tree->mcmc->acc_move[tree->mcmc->num_move_birth_rate]++;
        }
    }
  tree->mcmc->run_move[tree->mcmc->num_move_birth_rate]++;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MCMC_Death_Rate(t_tree *tree)
{
  
  /* MCMC_Single_Param_Generic(&(tree->rates->death_rate), */
  /*       		    0.0, // instead of tree->rates->death_rate_min as death rate can be equal to 0 (Yule model) */
  /*       		    tree->rates->death_rate_max, */
  /*       		    tree->mcmc->num_move_death_rate, */
  /*       		    &(tree->rates->c_lnL_times),NULL, */
  /*       		    Wrap_Lk_Times,NULL,tree->mcmc->move_type[tree->mcmc->num_move_death_rate],NO,NULL,tree,NULL); */
 
  phydbl cur_death_rate,new_death_rate;
  phydbl cur_lnL_time,new_lnL_time;
  phydbl cur_lnL_time_ghost,new_lnL_time_ghost;
  /* phydbl cur_lnL_time_pivot,new_lnL_time_pivot; */
  phydbl u,alpha,ratio;
  phydbl death_rate_min,death_rate_max;
  phydbl K;
  int i,n_mcmc_steps,move;

  cur_death_rate = -1.0;
  new_death_rate = -1.0;
  ratio          =  0.0;
  n_mcmc_steps   =  tree->mcmc->run_move[tree->mcmc->num_move_death_rate] == 1 ? 1000 : 100;
  move           = -1;

  K = tree->mcmc->tune_move[tree->mcmc->num_move_death_rate];

  cur_lnL_time = tree->rates->c_lnL_times;
  new_lnL_time = tree->rates->c_lnL_times;

  cur_lnL_time_ghost = UNLIKELY;
  new_lnL_time_ghost = UNLIKELY;

  /* cur_lnL_time_pivot = UNLIKELY; */
  /* new_lnL_time_pivot = UNLIKELY; */
  
  cur_death_rate = tree->rates->death_rate;

  death_rate_min = tree->rates->death_rate_min;
  death_rate_max = MIN(tree->rates->death_rate_max,tree->rates->birth_rate);

  MCMC_Make_Move(&cur_death_rate,&new_death_rate,death_rate_min,death_rate_max,&ratio,K,tree->mcmc->move_type[tree->mcmc->num_move_death_rate]);
  /* new_death_rate = Uni()*(death_rate_max - death_rate_min) + death_rate_min; */
  
  if(new_death_rate < death_rate_max && new_death_rate > death_rate_min && new_death_rate < tree->rates->birth_rate)
    {
      tree->rates->death_rate = new_death_rate;
      new_lnL_time = TIMES_Lk_Times(NO,tree);
      ratio += (new_lnL_time - cur_lnL_time);

      if(tree->mcmc->run_move[tree->mcmc->num_move_birth_rate] >= 0)
        {
          /* if(tree->mcmc->run_move[tree->mcmc->num_move_birth_rate] >= 500) */
          /*   { */
          /*     tree->extra_tree->rates->birth_rate = tree->rates->birth_rate_pivot; */
          /*     tree->extra_tree->rates->death_rate = tree->rates->death_rate_pivot; */
          /*     cur_lnL_time_pivot = TIMES_Lk_Times(NO,tree->extra_tree); */
          /*   } */
          
          /* tree->extra_tree->rates->death_rate = cur_death_rate; */
          /* tree->extra_tree->rates->birth_rate = tree->rates->birth_rate; */
          /* cur_lnL_time_ghost = TIMES_Lk_Times(NO,tree->extra_tree); */
          
          /* Copy_Tree(tree->extra_tree,tree->extra_tree->extra_tree); */
          /* RATES_Copy_Rate_Struct(tree->extra_tree->rates,tree->extra_tree->extra_tree->rates,tree->n_otu); */
          /* DATE_Assign_Primary_Calibration(tree->extra_tree->extra_tree); */
          


          tree->extra_tree->eval_alnL          = NO;
          tree->extra_tree->eval_rlnL          = NO;
          tree->extra_tree->eval_glnL          = YES;
          tree->extra_tree->rates->c_lnL_rates = UNLIKELY;
          tree->extra_tree->c_lnL              = UNLIKELY;

          TIMES_Randomize_Tree_With_Time_Constraints(tree->extra_tree->rates->a_cal[0],tree->extra_tree);
          tree->extra_tree->rates->birth_rate = tree->rates->birth_rate;
          tree->extra_tree->rates->death_rate = new_death_rate;
          TIMES_Lk_Times(NO,tree->extra_tree);
          
          if(!(tree->extra_tree->rates->c_lnL_times > UNLIKELY))
            {
              PhyML_Fprintf(stderr,"\n. glnL=%f",tree->extra_tree->rates->c_lnL_times);
              PhyML_Fprintf(stderr,"\n. death=%G birth=%G [%G]",new_death_rate,tree->rates->birth_rate,tree->extra_tree->rates->birth_rate);
              PhyML_Fprintf(stderr,"\n");
              for(int i=0;i<tree->extra_tree->rates->n_cal;++i)
                {
                  t_cal *cal = tree->extra_tree->rates->a_cal[i];
                  PhyML_Fprintf(stderr,"\n. calibration %s applies to clade %s\t",cal->id,cal->clade_list[cal->current_clade_idx]->id);
                  for(int j=0;j<cal->clade_list_size;++j)
                    {
                      t_clad *clade = cal->clade_list[j];
                      PhyML_Fprintf(stderr,"time:%G\t",tree->extra_tree->rates->nd_t[clade->target_nd->num]);
                    }
                }

              PhyML_Fprintf(stderr,"\n");
              TIMES_Lk_Times(YES,tree->extra_tree);
              assert(FALSE);
            }
                    
          i = 0;
          do
            {
              u = Uni();
              for(move=0;move<tree->mcmc->n_moves;move++) if(tree->mcmc->move_weight[move] > u-1.E-10) break;
              
              /* if(!(i%10)) PhyML_Printf("\n>> %5d Move '%20s' %12f",i,tree->mcmc->move_name[move],tree->extra_tree->rates->c_lnL_times); */
              
              if(!strcmp(tree->mcmc->move_name[move],"tree_height")) { MCMC_Tree_Height(tree->extra_tree); i++; }
              if(!strcmp(tree->mcmc->move_name[move],"times"))       { MCMC_Times_All(tree->extra_tree); i++; }
              if(!strcmp(tree->mcmc->move_name[move],"spr"))         { MCMC_Prune_Regraft(tree->extra_tree); i++; }
              if(!strcmp(tree->mcmc->move_name[move],"spr_local"))   { MCMC_Prune_Regraft_Local(tree->extra_tree); i++; }
              
              if(!(tree->extra_tree->rates->c_lnL_times > UNLIKELY))
                {
                  PhyML_Fprintf(stderr,"\n. move: %s",tree->mcmc->move_name[move]);
                  PhyML_Fprintf(stderr,"\n. glnL=%f",tree->extra_tree->rates->c_lnL_times);
                  TIMES_Lk_Times(YES,tree->extra_tree);
                  Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
                }
            }
          while(i < n_mcmc_steps);
          
          tree->extra_tree->rates->death_rate = cur_death_rate;
          tree->extra_tree->rates->birth_rate = tree->rates->birth_rate;
          cur_lnL_time_ghost = TIMES_Lk_Times(NO,tree->extra_tree);

          tree->extra_tree->rates->birth_rate = tree->rates->birth_rate;
          tree->extra_tree->rates->death_rate = new_death_rate;
          new_lnL_time_ghost = TIMES_Lk_Times(NO,tree->extra_tree);

          ratio += (cur_lnL_time_ghost - new_lnL_time_ghost);
          
          /* if(tree->mcmc->run_move[tree->mcmc->num_move_birth_rate] >= 500) */
          /*   { */
          /*     tree->extra_tree->rates->birth_rate = tree->rates->birth_rate_pivot; */
          /*     tree->extra_tree->rates->death_rate = tree->rates->death_rate_pivot; */
          /*     new_lnL_time_pivot = TIMES_Lk_Times(NO,tree->extra_tree); */
          /*     ratio += (new_lnL_time_pivot - cur_lnL_time_pivot); */
          /*   } */
        }
      
      ratio = exp(ratio);
      alpha = MIN(1.,ratio);
      
      /* printf("\n.  d :%4d: %12G -> %12G ratio : %12G [ghost: %12G %12G -- real: %12G %12G]", */
      /*        tree->mcmc->run_move[tree->mcmc->num_move_death_rate], */
      /*        cur_death_rate, */
      /*        new_death_rate, */
      /*        ratio, */
      /*        cur_lnL_time_ghost, */
      /*        new_lnL_time_ghost, */
      /*        cur_lnL_time, */
      /*        new_lnL_time); */


      u = Uni();

      assert(isnan(u) == NO && isinf(fabs(u)) == NO);

      if(u > alpha) /* Reject */
        {
          /* PhyML_Printf("  reject"); */
          tree->rates->death_rate  = cur_death_rate;
          tree->rates->c_lnL_times = cur_lnL_time;

          /* Copy_Tree(tree->extra_tree->extra_tree,tree->extra_tree); */
          /* RATES_Copy_Rate_Struct(tree->extra_tree->extra_tree->rates,tree->extra_tree->rates,tree->n_otu); */
          /* DATE_Assign_Primary_Calibration(tree->extra_tree); */
        }
      else
        {
          /* PhyML_Printf("  accept"); */
          tree->mcmc->acc_move[tree->mcmc->num_move_death_rate]++;
        }
    }
  tree->mcmc->run_move[tree->mcmc->num_move_death_rate]++;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MCMC_Birth_Death_Updown(t_tree *tree)
{ 
  phydbl cur_death_rate,new_death_rate;
  phydbl cur_birth_rate,new_birth_rate;
  phydbl cur_lnL_time,new_lnL_time;
  phydbl cur_lnL_time_ghost,new_lnL_time_ghost;
  /* phydbl cur_lnL_time_pivot,new_lnL_time_pivot; */
  phydbl u,alpha,ratio;
  phydbl death_rate_min,death_rate_max;
  phydbl birth_rate_min,birth_rate_max;
  phydbl K,scale;
  int i,n_mcmc_steps,move;

  cur_death_rate = -1.0;
  new_death_rate = -1.0;
  cur_birth_rate = -1.0;
  new_birth_rate = -1.0;
  ratio          =  0.0;
  n_mcmc_steps   =  tree->mcmc->run_move[tree->mcmc->num_move_birth_death_updown] == 1 ? 1000 : 100;
  move           = -1;

  K = tree->mcmc->tune_move[tree->mcmc->num_move_birth_death_updown];

  cur_lnL_time = tree->rates->c_lnL_times;
  new_lnL_time = tree->rates->c_lnL_times;

  cur_lnL_time_ghost = UNLIKELY;
  new_lnL_time_ghost = UNLIKELY;

  /* cur_lnL_time_pivot = UNLIKELY; */
  /* new_lnL_time_pivot = UNLIKELY; */
  
  cur_death_rate = tree->rates->death_rate;
  cur_birth_rate = tree->rates->birth_rate;

  death_rate_min = tree->rates->death_rate_min;
  death_rate_max = tree->rates->death_rate_max;

  birth_rate_min = tree->rates->birth_rate_min;
  birth_rate_max = tree->rates->birth_rate_max;

  scale = exp(K*(Uni()-.5));
  new_birth_rate = cur_birth_rate * scale;
  new_death_rate = cur_death_rate * scale;
  ratio += 2.*log(scale);
  
  if(new_death_rate < death_rate_max && new_death_rate > death_rate_min &&
     new_birth_rate < birth_rate_max && new_birth_rate > birth_rate_min &&
     new_death_rate < new_birth_rate)
    {
      tree->rates->death_rate = new_death_rate;
      tree->rates->birth_rate = new_birth_rate;      
      new_lnL_time = TIMES_Lk_Times(NO,tree);
      ratio += (new_lnL_time - cur_lnL_time);      

      if(tree->mcmc->run_move[tree->mcmc->num_move_birth_rate] >= 0)
        {
          /* if(tree->mcmc->run_move[tree->mcmc->num_move_birth_rate] >= 500) */
          /*   { */
          /*     tree->extra_tree->rates->death_rate = tree->rates->death_rate_pivot; */
          /*     tree->extra_tree->rates->birth_rate = tree->rates->birth_rate_pivot; */
          /*     cur_lnL_time_pivot = TIMES_Lk_Times(NO,tree->extra_tree); */
          /*   } */
          
          /* tree->extra_tree->rates->death_rate = cur_death_rate; */
          /* tree->extra_tree->rates->birth_rate = cur_birth_rate; */
          /* cur_lnL_time_ghost = TIMES_Lk_Times(NO,tree->extra_tree); */
          
          /* Copy_Tree(tree->extra_tree,tree->extra_tree->extra_tree); */
          /* RATES_Copy_Rate_Struct(tree->extra_tree->rates,tree->extra_tree->extra_tree->rates,tree->n_otu); */
          /* DATE_Assign_Primary_Calibration(tree->extra_tree->extra_tree); */
          
          tree->extra_tree->eval_alnL          = NO;
          tree->extra_tree->eval_rlnL          = NO;
          tree->extra_tree->eval_glnL          = YES;
          tree->extra_tree->rates->c_lnL_rates = UNLIKELY;
          tree->extra_tree->c_lnL              = UNLIKELY;
          
          TIMES_Randomize_Tree_With_Time_Constraints(tree->extra_tree->rates->a_cal[0],tree->extra_tree);
          tree->extra_tree->rates->death_rate = new_death_rate;
          tree->extra_tree->rates->birth_rate = new_birth_rate;
          TIMES_Lk_Times(NO,tree->extra_tree);
          
          if(!(tree->extra_tree->rates->c_lnL_times > UNLIKELY))
            {
              PhyML_Fprintf(stderr,"\n. glnL=%f",tree->extra_tree->rates->c_lnL_times);
              PhyML_Fprintf(stderr,"\n. death=%G birth=%G [%G]",new_death_rate,tree->rates->birth_rate,tree->extra_tree->rates->birth_rate);
              TIMES_Lk_Times(YES,tree->extra_tree);
              assert(FALSE);
            }
          
          
          i = 0;
          do
            {
              u = Uni();
              for(move=0;move<tree->mcmc->n_moves;move++) if(tree->mcmc->move_weight[move] > u-1.E-10) break;
              
              /* PhyML_Printf("\n>> Move '%s' %f %f b:%f d:%f-%f", */
              /*              tree->mcmc->move_name[move], */
              /*              tree->extra_tree->rates->c_lnL_times, */
              /*              tree->rates->c_lnL_times, */
              /*              tree->rates->birth_rate, */
              /*              cur_death_rate, */
              /*              new_death_rate); */
              
              if(!strcmp(tree->mcmc->move_name[move],"tree_height")) MCMC_Tree_Height(tree->extra_tree);
              if(!strcmp(tree->mcmc->move_name[move],"times"))       MCMC_Times_All(tree->extra_tree);
              if(!strcmp(tree->mcmc->move_name[move],"spr"))         MCMC_Prune_Regraft(tree->extra_tree);
              if(!strcmp(tree->mcmc->move_name[move],"spr_local"))   MCMC_Prune_Regraft_Local(tree->extra_tree);
              
              if(!(tree->extra_tree->rates->c_lnL_times > UNLIKELY))
                {
                  PhyML_Fprintf(stderr,"\n. move: %s",tree->mcmc->move_name[move]);
                  PhyML_Fprintf(stderr,"\n. glnL=%f",tree->extra_tree->rates->c_lnL_times);
                  TIMES_Lk_Times(YES,tree->extra_tree);
                  Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
                }
              i++;
            }
          while(i < n_mcmc_steps);
          
          tree->extra_tree->rates->death_rate = cur_death_rate;
          tree->extra_tree->rates->birth_rate = cur_birth_rate;
          cur_lnL_time_ghost = TIMES_Lk_Times(NO,tree->extra_tree);

          tree->extra_tree->rates->death_rate = new_death_rate;
          tree->extra_tree->rates->birth_rate = new_birth_rate;
          new_lnL_time_ghost = TIMES_Lk_Times(NO,tree->extra_tree);

          ratio += (cur_lnL_time_ghost - new_lnL_time_ghost);
          
          /* if(tree->mcmc->run_move[tree->mcmc->num_move_birth_rate] == 500) */
          /*   { */
          /*     tree->extra_tree->rates->death_rate = tree->rates->death_rate_pivot; */
          /*     tree->extra_tree->rates->birth_rate = tree->rates->birth_rate_pivot; */
          /*     new_lnL_time_pivot = TIMES_Lk_Times(NO,tree->extra_tree); */
          /*     ratio += (new_lnL_time_pivot - cur_lnL_time_pivot); */
          /*   } */
        }
      
      ratio = exp(ratio);
      alpha = MIN(1.,ratio);
      
      /* printf("\n. bd :%4d: %12G -> %12G ratio : %12G [ghost: %12G %12G -- real: %12G %12G]", */
      /*        tree->mcmc->run_move[tree->mcmc->num_move_birth_death_updown], */
      /*        cur_death_rate, */
      /*        new_death_rate, */
      /*        ratio, */
      /*        cur_lnL_time_ghost, */
      /*        new_lnL_time_ghost, */
      /*        cur_lnL_time, */
      /*        new_lnL_time); */
      
      
      u = Uni();

      assert(isnan(u) == NO && isinf(fabs(u)) == NO);

      if(u > alpha) /* Reject */
        {
          /* PhyML_Printf("  reject"); */
          tree->rates->death_rate  = cur_death_rate;
          tree->rates->birth_rate  = cur_birth_rate;
          tree->rates->c_lnL_times = cur_lnL_time;

          /* Copy_Tree(tree->extra_tree->extra_tree,tree->extra_tree); */
          /* RATES_Copy_Rate_Struct(tree->extra_tree->extra_tree->rates,tree->extra_tree->rates,tree->n_otu); */
          /* DATE_Assign_Primary_Calibration(tree->extra_tree); */
        }
      else
        {
          /* PhyML_Printf("  accept"); */
          tree->mcmc->acc_move[tree->mcmc->num_move_birth_death_updown]++;
        }
    }
  tree->mcmc->run_move[tree->mcmc->num_move_birth_death_updown]++;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MCMC_Nu(t_tree *tree)
{
  phydbl cur_nu,new_nu,cur_lnL_rate,new_lnL_rate;
  phydbl u,alpha,ratio;
  phydbl min_nu,max_nu;
  phydbl K;
  phydbl cur_lnL_seq, new_lnL_seq;

  cur_nu        = -1.0;
  new_nu        = -1.0;
  ratio         =  0.0;

  K = tree->mcmc->tune_move[tree->mcmc->num_move_nu];

  cur_lnL_rate = tree->rates->c_lnL_rates;
  new_lnL_rate = tree->rates->c_lnL_rates;

  cur_lnL_seq = tree->c_lnL;
  new_lnL_seq = tree->c_lnL;
  
  cur_nu      = tree->rates->nu;

  min_nu = tree->rates->min_nu;
  max_nu = tree->rates->max_nu;

  MCMC_Make_Move(&cur_nu,&new_nu,min_nu,max_nu,&ratio,K,tree->mcmc->move_type[tree->mcmc->num_move_nu]);
  
  if(new_nu < max_nu && new_nu > min_nu) 
    {
      tree->rates->nu = new_nu;
      
      if(tree->eval_rlnL == YES) new_lnL_rate = RATES_Lk_Rates(tree);      

      if(tree->rates->model == GUINDON)
	{
          RATES_Update_Cur_Bl(tree);
	  if(tree->eval_alnL == YES) new_lnL_seq = Lk(NULL,tree);
	}
      else
        new_lnL_seq = cur_lnL_seq;
      
      ratio += (new_lnL_rate - cur_lnL_rate);
      ratio += (new_lnL_seq - cur_lnL_seq);
 
      ratio = exp(ratio);
      alpha = MIN(1.,ratio);
      
      u = Uni();

      assert(isnan(u) == NO && isinf(fabs(u)) == NO);

      if(u > alpha) /* Reject */
	{
	  tree->rates->nu          = cur_nu;
	  tree->rates->c_lnL_rates = cur_lnL_rate;
	  tree->c_lnL              = cur_lnL_seq;
          if(tree->rates->model == GUINDON && tree->eval_alnL == YES) RATES_Update_Cur_Bl(tree);
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

void MCMC_Clade_Change(t_tree *tree)
{
  phydbl u,alpha,ratio;
  phydbl cur_lnL_calib, new_lnL_calib;
  phydbl cur_lnL_time, new_lnL_time;
  int target_cal_idx,current_clade_idx,new_clade_idx;
  t_cal *cal;
  
  ratio         = 0.0;
  cal           = NULL;
  cur_lnL_calib = DATE_Lk_Calib(tree);

  new_lnL_time = tree->rates->c_lnL_times;
  cur_lnL_time = tree->rates->c_lnL_times;
      

  
  // Choose a calibration uniformly at random
  target_cal_idx = Rand_Int(0,tree->rates->n_cal-1);
  cal = tree->rates->a_cal[target_cal_idx];

  /* printf("\n. CURRENT cal: %s clade: %s target: %d",cal->id,cal->clade_list[cal->current_clade_idx]->id,cal->clade_list[cal->current_clade_idx]->target_nd->num); */

  // Choose a new clade uniformly at random
  current_clade_idx = cal->current_clade_idx;
  new_clade_idx = Rand_Int(0,cal->clade_list_size-1);
  cal->current_clade_idx = new_clade_idx;

  new_lnL_time = TIMES_Lk_Times(NO,tree);
  ratio += new_lnL_time - cur_lnL_time;
  
  /* printf("\n. NEW cal: %s clade: %s target: %d",cal->id,cal->clade_list[cal->current_clade_idx]->id,cal->clade_list[cal->current_clade_idx]->target_nd->num); */

  new_lnL_calib = DATE_Lk_Calib(tree);
  ratio += new_lnL_calib - cur_lnL_calib;
      
  ratio = exp(ratio);
  alpha = MIN(1.,ratio);
  
  u = Uni();

  assert(isnan(u) == NO && isinf(fabs(u)) == NO);

  if(u > alpha) /* Reject */
    {
      /* printf("\n. reject\n"); */
      cal->current_clade_idx = current_clade_idx;
      tree->rates->c_lnL_times = cur_lnL_time;
      TIMES_Lk_Times(NO,tree); // Required in order to update node targeted by selected calibration
    }
  else
    {
      /* printf("\n. accept\n"); */
      tree->mcmc->acc_move[tree->mcmc->num_move_clade_change]++;
    }
  tree->mcmc->run_move[tree->mcmc->num_move_clade_change]++;
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MCMC_All_Rates(t_tree *tree)
{
  phydbl cur_lnL_seq, new_lnL_seq, cur_lnL_rate;
  phydbl u, ratio, alpha;

  new_lnL_seq = tree->c_lnL;
  cur_lnL_seq = tree->c_lnL;
  cur_lnL_rate = tree->rates->c_lnL_rates;
  ratio        = 0.0;
  
  Record_Br_Len(tree);
  RATES_Record_Rates(tree);
  
  MCMC_Sim_Rate(tree->n_root,tree->n_root->v[2],tree);
  MCMC_Sim_Rate(tree->n_root,tree->n_root->v[1],tree);

  if(tree->eval_alnL == YES) new_lnL_seq = Lk(NULL,tree);
  
  ratio += (new_lnL_seq - cur_lnL_seq);
  ratio = exp(ratio);

  alpha = MIN(1.,ratio);
  u = Uni();

  assert(isnan(u) == NO && isinf(fabs(u)) == NO);

  if(u > alpha) /* Reject */
    {
      Restore_Br_Len(tree);
      RATES_Reset_Rates(tree);
      tree->rates->c_lnL_rates = cur_lnL_rate;
      tree->c_lnL = cur_lnL_seq;
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
  
  mean = br_r_a;
      
  if(tree->rates->model == STRICTCLOCK) tree->rates->br_r[d->num] = 1.0;
  else
    {
      tree->rates->br_r[d->num] = Rnorm_Trunc(mean,sd,tree->rates->min_rate,tree->rates->max_rate,&err);
      if(err == YES) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
    }
        
  if(d->tax) return;
  else
    {
      int i;

      for(i=0;i<3;i++)
	if(d->v[i] != a && d->b[i] != tree->e_root)
	  MCMC_Sim_Rate(d,d->v[i],tree);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MCMC_Prune_Regraft(t_tree *tree)
{
  phydbl u,alpha,ratio;
  phydbl t_min,t_max;
  phydbl cur_lnL_seq,new_lnL_seq;
  phydbl cur_lnL_rate,new_lnL_rate;
  phydbl cur_lnL_time,new_lnL_time;
  phydbl new_t;
  int i,prune_idx,n_iter,n_regraft_nd,regraft_idx,dir_prune;
  phydbl *times;
  int rnd_dir,dir_v1,dir_v2;
  t_node *prune,*prune_daughter,*new_regraft_nd,*cur_regraft_nd;
  t_ll *regraft_nd_list;
  t_edge *target, *ori_target, *residual,*regraft_edge;
  phydbl regraft_t_min,regraft_t_max;
  /* phydbl *prob_idx; */
  /* phydbl r,sum,prob_select_fwd,prob_select_bwd; */
  
  n_iter = MAX(1,(int)(tree->n_otu/5));
  /* n_iter = tree->n_otu-2; */
  times = tree->rates->nd_t;

  /* prob_idx = (phydbl *)mCalloc(tree->n_otu-1,sizeof(phydbl)); */

  // Shallowest node has probability proportional to r to be
  // to be selected as prune node while deepest a proba prop to 1.
  /* r = 0.5; */
  /* for(i=0;i<tree->n_otu-1;i++) prob_idx[i] = i*(r-1.0)/(tree->n_otu-2) + 1.0; */

  /* sum = .0; */
  /* for(i=0;i<tree->n_otu-1;i++) sum += prob_idx[i]; */
  /* for(i=0;i<tree->n_otu-1;i++) prob_idx[i] /= sum; */

  
  while(n_iter--)
    {      

      TIMES_Update_Node_Ordering(tree);

      /* if(tree->eval_alnL == YES) Lk(NULL,tree); */
      
      tree->mcmc->run_move[tree->mcmc->num_move_spr]++;
      
      RATES_Record_Times(tree);
      
      cur_lnL_seq  = tree->c_lnL;
      new_lnL_seq  = tree->c_lnL;
      cur_lnL_rate = tree->rates->c_lnL_rates;
      new_lnL_rate = tree->rates->c_lnL_rates;
      cur_lnL_time = tree->rates->c_lnL_times;
      new_lnL_time = tree->rates->c_lnL_times;

      ratio          = 0.0;
      regraft_edge   = NULL;
      new_regraft_nd = NULL;
      cur_regraft_nd = NULL;
      new_t          = 0.0;

      // Select prune node (any internal node)
      prune_idx = Rand_Int(tree->n_otu,2*tree->n_otu-2);

      /* prune_idx = Sample_i_With_Proba_pi(prob_idx,tree->n_otu-1); */
      /* prob_select_fwd = prob_idx[prune_idx]; */
      /* ratio -= log(prob_select_fwd); */
      /* prune_idx = tree->rates->t_rank[prune_idx]; */
      
      prune = tree->a_nodes[prune_idx];

      
      assert(prune && prune->tax == NO);

      // Select a daughter of prune node
      dir_v1 = dir_v2 = -1;
      for(i=0;i<3;i++) 
        if(prune->v[i] != prune->anc && prune->b[i] != tree->e_root)
          {
            if(dir_v1 < 0) dir_v1 = i;
            else           dir_v2 = i;
          }
      
      u = Uni();
      if(u < 0.5) rnd_dir = dir_v1;
      else        rnd_dir = dir_v2;
          
      prune_daughter = prune->v[rnd_dir];
      cur_regraft_nd = prune->v[rnd_dir == dir_v1 ? dir_v2 : dir_v1];            
      
      if(prune == tree->n_root)
        {
          if(prune_daughter == prune->v[dir_v1] && prune->v[dir_v2]->tax == YES)
            {
              prune_daughter = prune->v[dir_v2];
              cur_regraft_nd = prune->v[dir_v1];
            }

          if(prune_daughter == prune->v[dir_v2] && prune->v[dir_v1]->tax == YES)
            {
              prune_daughter = prune->v[dir_v1];
              cur_regraft_nd = prune->v[dir_v2];
            }
        }
           
      if(prune_daughter->anc != prune) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);

      dir_prune = -1;
      for(i=0;i<3;i++)
        {
          if(prune_daughter->v[i] == prune || prune_daughter->b[i] == tree->e_root)
            {
              dir_prune = i;
              break;
            }
        }
      assert(dir_prune > -1);
            

      // Get the list of potential regraft nodes (oldest node on regraft edge)
      regraft_nd_list = DATE_List_Of_Regraft_Nodes(prune_daughter->v[dir_prune],prune_daughter,&regraft_t_min,&regraft_t_max,NO,tree);
      assert(regraft_nd_list);
      if(prune == tree->n_root) Push_Bottom_Linked_List(prune,&regraft_nd_list,YES);

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
                 prune_daughter->v[dir_prune]->num,
                 tree->rates->c_lnL_times); fflush(NULL);

          regraft_nd_list = DATE_List_Of_Regraft_Nodes(prune_daughter->v[dir_prune],prune_daughter,&regraft_t_min,&regraft_t_max,YES,tree);

          assert(FALSE);
        }
      
      // Randomly select one (uniform)
      regraft_idx = Rand_Int(0,n_regraft_nd-1);      
      new_regraft_nd = Linked_List_Elem(regraft_idx,regraft_nd_list);      
      Free_Linked_List(regraft_nd_list);

            
      // Time of regraft node and corresponding (partial) Hastings ratio
      t_max = MIN(times[prune_daughter->num],times[cur_regraft_nd->num]);      
      if(prune == tree->n_root) t_min = 10.0*t_max;
      else t_min = times[prune->anc->num];
      t_min = MAX(t_min,regraft_t_min);
      ratio += log(1./(t_max - t_min));

      t_max = MIN(times[prune_daughter->num],times[new_regraft_nd->num]);      
      if(new_regraft_nd == tree->n_root) t_min = 10.0*t_max;
      else t_min = times[new_regraft_nd->anc->num];
      t_min = MAX(t_min,regraft_t_min);
      ratio -= log(1./(t_max - t_min));
      
      new_t = Uni()*(t_max-t_min) + t_min;

      // Age of root node changes when pruned subtree is on one side of that node
      // Change here, not after the prune and regraft move
      /* if(prune == tree->n_root) */
      /*   { */
      /*     if(prune_daughter == tree->n_root->v[1]) */
      /*       times[tree->n_root->num] = times[tree->n_root->v[2]->num]; */
      /*     else if(prune_daughter == tree->n_root->v[2]) */
      /*       times[tree->n_root->num] = times[tree->n_root->v[1]->num]; */
      /*     else assert(false); */
      /*   } */
      
          
      // New age
      if(prune == tree->n_root || new_regraft_nd == tree->n_root)
        {
          if(prune == tree->n_root)
            {
              if(prune_daughter == tree->n_root->v[1]) times[tree->n_root->num] = times[tree->n_root->v[2]->num];
              else                                     times[tree->n_root->num] = times[tree->n_root->v[1]->num];
              times[prune_daughter->v[dir_prune]->num] = new_t;
            }
          if(new_regraft_nd == tree->n_root)
            {
              times[prune_daughter->v[dir_prune]->num] = times[tree->n_root->num];
              times[tree->n_root->num] = new_t;
            }
        }
      else
        {
          times[prune->num] = new_t;
        }


      // Prune
      target = residual = NULL;
      Prune_Subtree(prune_daughter->v[dir_prune],
                    prune_daughter,
                    &target,&residual,tree);
      ori_target = target;
      

      // Regraft edge is the one sitting above regraft_nd
      if(new_regraft_nd == tree->n_root->v[1] ||
         new_regraft_nd == tree->n_root->v[2] ||
         new_regraft_nd == tree->n_root) regraft_edge = tree->e_root;
      else
        {
          for(i=0;i<3;i++) if(new_regraft_nd->v[i] == new_regraft_nd->anc) break;
          assert(i!=3);
          regraft_edge = new_regraft_nd->b[i];
        }

      assert(regraft_edge);      
      assert(residual->left != residual->rght);
      assert(regraft_edge->left != prune_daughter->v[dir_prune]);
      assert(regraft_edge->rght != prune_daughter->v[dir_prune]);
      

      // Regraft
      Graft_Subtree(regraft_edge,
                    prune_daughter->v[dir_prune],
                    prune_daughter,
                    residual,
                    new_regraft_nd,tree);
      

      /* TIMES_Update_Node_Ordering(tree); */
      /* For(i,2*tree->n_otu-2) if(prune->num == tree->rates->t_rank[i]) break; */
      /* prob_select_bwd = prob_idx[i]; */
      /* ratio += log(prob_select_bwd); */

      if(!TIMES_Check_Node_Height_Ordering(tree))
        {
          printf("\n. prune[%d]->t:%.3f daughter[%d]->t:%.3f prune_anc[%d]->t:%.3f regraft[%d]->t:%.3f regraft_anc[%d]->t:%.3f [effective:%d] t_prior_min/max: [prune:[%.3f %.3f] regraft:[%.3f %.3f]] ",
                 prune->num,
                 times[prune->num],
                 prune_daughter->num,
                 times[prune_daughter->num],
                 prune->anc ? prune->anc->num : -1,
                 prune->anc ? times[prune->anc->num] : -1.,
                 new_regraft_nd->num,
                 times[new_regraft_nd->num],
                 new_regraft_nd->anc ? new_regraft_nd->anc->num : -1,
                 new_regraft_nd->anc ? times[new_regraft_nd->anc->num] : +1.,
                 prune->num,
                 tree->rates->t_prior_min[prune->num],
                 tree->rates->t_prior_max[prune->num],
                 tree->rates->t_prior_min[new_regraft_nd->num],
                 tree->rates->t_prior_max[new_regraft_nd->num]);
          PhyML_Fprintf(stderr,"\n. root: %d %d %d",tree->n_root->num,tree->n_root->v[1]->num,tree->n_root->v[2]->num);
          Generic_Exit(__FILE__,__LINE__,__FUNCTION__);    
        }

      DATE_Assign_Primary_Calibration(tree);
      RATES_Update_Cur_Bl(tree);

      if(tree->eval_glnL == YES) new_lnL_time = TIMES_Lk_Times(NO,tree); 

      if(new_lnL_time > UNLIKELY)
        {
          Set_Both_Sides(NO,tree);
          if(tree->eval_alnL == YES) new_lnL_seq = Lk(NULL,tree);
          if(tree->eval_rlnL == YES) new_lnL_rate = RATES_Lk_Rates(tree);
        }
      
      ratio += (new_lnL_seq - cur_lnL_seq);
      ratio += (new_lnL_rate - cur_lnL_rate);
      ratio += (new_lnL_time - cur_lnL_time);
      
      ratio = exp(ratio);
      alpha = MIN(1.,ratio);

      /* Always accept move */
      if(tree->mcmc->always_yes == YES && new_lnL_time > UNLIKELY) alpha = 1.0;
      
      u = Uni();

      assert(isnan(u) == NO && isinf(fabs(u)) == NO);

      if(u > alpha)
        {
          // Reject
          Prune_Subtree(prune_daughter->v[dir_prune],
                        prune_daughter,
                        &target,&residual,tree);
          assert(residual->left != residual->rght);
          assert(ori_target->left != prune_daughter->v[dir_prune]);
          assert(ori_target->rght != prune_daughter->v[dir_prune]);
          Graft_Subtree(ori_target,
                        prune_daughter->v[dir_prune],
                        prune_daughter,residual,prune == tree->n_root ? tree->n_root : cur_regraft_nd,tree);

          RATES_Reset_Times(tree);
          RATES_Update_Cur_Bl(tree);
          DATE_Assign_Primary_Calibration(tree);
          
          new_lnL_time = TIMES_Lk_Times(NO,tree); 
          if(Are_Equal(new_lnL_time,cur_lnL_time,1.E-5) == NO)
            {
              PhyML_Printf("\n. new_lnL_time: %f cur_lnL_time: %f",new_lnL_time,cur_lnL_time);
              assert(FALSE);
            }

          /* if(tree->eval_alnL == YES) */
          /*   { */
          /*     new_lnL_seq = Lk(NULL,tree); */
          /*     if(Are_Equal(new_lnL_seq,cur_lnL_seq,1.E-5) == NO) */
          /*       { */
          /*         PhyML_Printf("\n. new: %f cur: %f",new_lnL_seq,cur_lnL_seq); */
          /*         assert(FALSE); */
          /*       } */
          /*   } */
              
          if(!(tree->rates->c_lnL_times > UNLIKELY))
            {
              printf("\n. time prune: %f",times[prune->num]);
              printf("\n. time prune_daughter: %f",times[prune_daughter->num]);
              printf("\n. prune: %d prune_daughter: %d prune_daughter->v[dir_prune]: %d cur_regraft_nd: %d new_regraft_nd: %d",
                     prune->num,
                     prune_daughter->num,
                     prune_daughter->v[dir_prune]->num,
                     cur_regraft_nd->num,
                     new_regraft_nd->num);
              TIMES_Lk_Times(YES,tree); 
              fflush(NULL);
            }
          assert(tree->rates->c_lnL_times > UNLIKELY);
          
          tree->c_lnL              = cur_lnL_seq;
          tree->rates->c_lnL_times = cur_lnL_time;
          tree->rates->c_lnL_rates = cur_lnL_rate;
        }
      else
        {
          tree->mcmc->acc_move[tree->mcmc->num_move_spr]++;
        }      
    }
  /* Free(prob_idx); */
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MCMC_Prune_Regraft_Weighted(t_tree *tree)
{
  phydbl u,alpha,ratio;
  phydbl cur_lnL_seq,new_lnL_seq;
  phydbl cur_lnL_rate,new_lnL_rate;
  phydbl cur_lnL_time,new_lnL_time;  
  phydbl new_t;
  int i,prune_idx,n_iter,dir_prune;
  phydbl *times;
  int rnd_dir,dir_v1,dir_v2;
  t_node *prune,*prune_daughter,*new_regraft_nd,*cur_regraft_nd;
  t_edge *target, *ori_target, *residual,*regraft_edge;
  phydbl radius;
  
  n_iter = MAX(1,(int)(tree->n_otu/5));
  times  = tree->rates->nd_t;
  
  while(n_iter--)
    {      
      tree->mcmc->run_move[tree->mcmc->num_move_spr_weighted]++;
      
      RATES_Record_Times(tree);
      
      cur_lnL_seq  = tree->c_lnL;
      new_lnL_seq  = tree->c_lnL;
      cur_lnL_rate = tree->rates->c_lnL_rates;
      new_lnL_rate = tree->rates->c_lnL_rates;
      cur_lnL_time = tree->rates->c_lnL_times;
      new_lnL_time = tree->rates->c_lnL_times;

      ratio          = 0.0;
      regraft_edge   = NULL;
      new_regraft_nd = NULL;
      cur_regraft_nd = NULL;
      new_t          = 0.0;

      // Select prune node (any internal node)
      prune_idx = Rand_Int(tree->n_otu,2*tree->n_otu-2);      
      prune = tree->a_nodes[prune_idx];

      assert(prune && prune->tax == NO);

      // Select a daughter of prune node
      dir_v1 = dir_v2 = -1;
      for(i=0;i<3;i++) 
        if(prune->v[i] != prune->anc && prune->b[i] != tree->e_root)
          {
            if(dir_v1 < 0) dir_v1 = i;
            else           dir_v2 = i;
          }
      
      u = Uni();
      if(u < 0.5) rnd_dir = dir_v1;
      else        rnd_dir = dir_v2;
          
      prune_daughter = prune->v[rnd_dir];
      cur_regraft_nd = prune->v[rnd_dir == dir_v1 ? dir_v2 : dir_v1];            

      if(prune == tree->n_root)
        {
          if(prune_daughter == prune->v[dir_v1] && prune->v[dir_v2]->tax == YES)
            {
              prune_daughter = prune->v[dir_v2];
              cur_regraft_nd = prune->v[dir_v1];
            }

          if(prune_daughter == prune->v[dir_v2] && prune->v[dir_v1]->tax == YES)
            {
              prune_daughter = prune->v[dir_v1];
              cur_regraft_nd = prune->v[dir_v2];
            }
        }
           
      assert(prune_daughter->anc == prune);

      dir_prune = -1;
      for(i=0;i<3;i++)
        {
          if(prune_daughter->v[i] == prune || prune_daughter->b[i] == tree->e_root)
            {
              dir_prune = i;
              break;
            }
        }
      assert(dir_prune > -1);


      radius = Rnorm(0.0,0.05);
      radius = fabs(radius);
      radius += tree->rates->cur_l[prune_daughter->num]; // valid too when prune == tree->n_root

      /* printf("\n. Init rad: %G l: %G %G root: %d",radius,tree->rates->cur_l[prune_daughter->num],prune_daughter->b[dir_prune]->l->v,prune_daughter->b[dir_prune] == tree->e_root); */
      // Get the list of potential regraft nodes (oldest node on regraft edge)
      Random_Walk_Along_Tree_On_Radius(prune_daughter,
                                       prune_daughter->v[dir_prune],                                       
                                       prune_daughter->b[dir_prune],
                                       &radius,
                                       &regraft_edge,
                                       &new_regraft_nd,
                                       &new_t,
                                       tree);
      
      
      /* printf("\n. new_regraft_edge: %d, new_regraft_nd: %d time: %f cur_regraft_nd: %d prune: %d prune_daughter: %d", */
      /*        regraft_edge ? regraft_edge->num : -1, */
      /*        new_regraft_nd ? new_regraft_nd->num : -1, */
      /*        new_t, */
      /*        cur_regraft_nd->num, */
      /*        prune->num, */
             /* prune_daughter->num); fflush(NULL); */

      if(new_regraft_nd == NULL) continue;
      if(new_regraft_nd == prune) continue;
      if(new_regraft_nd == cur_regraft_nd) continue;
      if(new_t > times[prune_daughter->num]) continue;
      assert(new_regraft_nd != prune_daughter);
      
      // New age
      if(prune == tree->n_root || new_regraft_nd == tree->n_root)
        {
          if(prune == tree->n_root)
            {
              if(prune_daughter == tree->n_root->v[1]) times[tree->n_root->num] = times[tree->n_root->v[2]->num];
              else                                     times[tree->n_root->num] = times[tree->n_root->v[1]->num];
              times[prune_daughter->v[dir_prune]->num] = new_t;
            }
          if(new_regraft_nd == tree->n_root)
            {
              times[prune_daughter->v[dir_prune]->num] = times[tree->n_root->num];
              times[tree->n_root->num] = new_t;
            }
        }
      else
        {
          times[prune->num] = new_t;
        }


      // Prune
      target = residual = NULL;
      Prune_Subtree(prune_daughter->v[dir_prune],
                    prune_daughter,
                    &target,&residual,tree);
      ori_target = target;
      

      // Regraft edge is the one sitting above new_regraft_nd
      if(new_regraft_nd == tree->n_root->v[1] ||
         new_regraft_nd == tree->n_root->v[2] ||
         new_regraft_nd == tree->n_root) regraft_edge = tree->e_root;
      else
        {
          for(i=0;i<3;++i) if(new_regraft_nd->v[i] == new_regraft_nd->anc) break;
          assert(i!=3);
          regraft_edge = new_regraft_nd->b[i];
        }

      assert(regraft_edge);
      assert(residual->left != residual->rght);
      assert(regraft_edge->left != prune_daughter->v[dir_prune]);
      assert(regraft_edge->rght != prune_daughter->v[dir_prune]);

      // Regraft
      Graft_Subtree(regraft_edge,
                    prune_daughter->v[dir_prune],
                    prune_daughter,
                    residual,
                    new_regraft_nd,tree);

      
      if(!TIMES_Check_Node_Height_Ordering(tree))
        {
          PhyML_Fprintf(stderr,"\n. prune[%d]->t:%.3f daughter[%d]->t:%.3f prune_anc[%d]->t:%.3f regraft[%d]->t:%.3f regraft_anc[%d]->t:%.3f [effective:%d] t_prior_min/max: [prune:[%.3f %.3f] regraft:[%.3f %.3f]] ",
                 prune->num,
                 times[prune->num],
                 prune_daughter->num,
                 times[prune_daughter->num],
                 prune->anc ? prune->anc->num : -1,
                 prune->anc ? times[prune->anc->num] : -1.,
                 new_regraft_nd->num,
                 times[new_regraft_nd->num],
                 new_regraft_nd->anc ? new_regraft_nd->anc->num : -1,
                 new_regraft_nd->anc ? times[new_regraft_nd->anc->num] : +1.,
                 prune->num,
                 tree->rates->t_prior_min[prune->num],
                 tree->rates->t_prior_max[prune->num],
                 tree->rates->t_prior_min[new_regraft_nd->num],
                 tree->rates->t_prior_max[new_regraft_nd->num]);
          PhyML_Fprintf(stderr,"\n. root: %d %d %d",tree->n_root->num,tree->n_root->v[1]->num,tree->n_root->v[2]->num);
          assert(FALSE);
        }

      RATES_Update_Cur_Bl(tree);
      DATE_Assign_Primary_Calibration(tree);
      if(tree->eval_glnL == YES) new_lnL_time = TIMES_Lk_Times(NO,tree); 

      if(new_lnL_time > UNLIKELY)
        {
          Set_Both_Sides(NO,tree);
          if(tree->eval_alnL == YES) new_lnL_seq = Lk(NULL,tree);
          if(tree->eval_rlnL == YES) new_lnL_rate = RATES_Lk_Rates(tree);
        }
      
      ratio += (new_lnL_seq - cur_lnL_seq);
      ratio += (new_lnL_rate - cur_lnL_rate);
      ratio += (new_lnL_time - cur_lnL_time);
      
      ratio = exp(ratio);
      alpha = MIN(1.,ratio);

      /* Always accept move */
      if(tree->mcmc->always_yes == YES && new_lnL_time > UNLIKELY) alpha = 1.0;
      
      u = Uni();

      assert(isnan(u) == NO && isinf(fabs(u)) == NO);

      if(u > alpha)
        {
          // Reject
          Prune_Subtree(prune_daughter->v[dir_prune],
                        prune_daughter,
                        &target,&residual,tree);
          assert(residual->left != residual->rght);
          assert(ori_target->left != prune_daughter->v[dir_prune]);
          assert(ori_target->rght != prune_daughter->v[dir_prune]);
          Graft_Subtree(ori_target,
                        prune_daughter->v[dir_prune],
                        prune_daughter,residual,prune == tree->n_root ? tree->n_root : cur_regraft_nd,tree);

          RATES_Reset_Times(tree);
          RATES_Update_Cur_Bl(tree);
          DATE_Assign_Primary_Calibration(tree);
          new_lnL_time = TIMES_Lk_Times(NO,tree); 
          if(Are_Equal(new_lnL_time,cur_lnL_time,1.E-5) == NO)
            {
              PhyML_Printf("\n. new_lnL_time: %f cur_lnL_time: %f",new_lnL_time,cur_lnL_time);
              assert(FALSE);
            }

          /* new_lnL_seq = Lk(NULL,tree); */
          /* if(Are_Equal(new_lnL_seq,cur_lnL_seq,1.E-5) == NO) */
          /*   { */
          /*     PhyML_Printf("\n. new: %f cur: %f",new_lnL_seq,cur_lnL_seq); */
          /*     assert(FALSE); */
          /*   } */
          
          if(!(tree->rates->c_lnL_times > UNLIKELY))
            {
              printf("\n. time prune: %f",times[prune->num]);
              printf("\n. time prune_daughter: %f",times[prune_daughter->num]);
              printf("\n. prune: %d prune_daughter: %d prune_daughter->v[dir_prune]: %d cur_regraft_nd: %d new_regraft_nd: %d",
                     prune->num,
                     prune_daughter->num,
                     prune_daughter->v[dir_prune]->num,
                     cur_regraft_nd->num,
                     new_regraft_nd->num);
              TIMES_Lk_Times(YES,tree); 
              fflush(NULL);
            }
          assert(tree->rates->c_lnL_times > UNLIKELY);
          
          tree->c_lnL              = cur_lnL_seq;
          tree->rates->c_lnL_times = cur_lnL_time;
          tree->rates->c_lnL_rates = cur_lnL_rate;
        }
      else
        {
          tree->mcmc->acc_move[tree->mcmc->num_move_spr_weighted]++;
        }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MCMC_Prune_Regraft_Local(t_tree *tree)
{
  phydbl u,alpha,ratio;
  phydbl t_min,t_max;
  phydbl cur_lnL_seq,new_lnL_seq;
  phydbl cur_lnL_rate,new_lnL_rate;
  phydbl cur_lnL_time,new_lnL_time;
  phydbl new_t;
  int i,prune_idx,n_iter,n_regraft_nd,regraft_idx,dir_prune;
  phydbl *times;
  int rnd_dir,dir_v1,dir_v2;
  t_node *prune,*prune_daughter,*new_regraft_nd,*cur_regraft_nd;
  t_ll *regraft_nd_list;
  t_edge *target, *ori_target, *residual,*regraft_edge;
  t_node *a, *b, *c, *d, *e;
  
  n_iter = MAX(1,(int)(tree->n_otu/5));
  times = tree->rates->nd_t;

  while(n_iter--)
    {
      tree->mcmc->run_move[tree->mcmc->num_move_spr_local]++;
      
      RATES_Record_Times(tree);
      
      cur_lnL_seq  = tree->c_lnL;
      new_lnL_seq  = tree->c_lnL;
      cur_lnL_rate = tree->rates->c_lnL_rates;
      new_lnL_rate = tree->rates->c_lnL_rates;
      cur_lnL_time = tree->rates->c_lnL_times;
      new_lnL_time = tree->rates->c_lnL_times;

      ratio          = 0.0;
      regraft_edge   = NULL;
      new_regraft_nd = NULL;
      cur_regraft_nd = NULL;
      new_t          = 0.0;

      // Select prune node (any internal node)
      prune_idx = Rand_Int(tree->n_otu,2*tree->n_otu-2);
      /* prune_idx = tree->n_otu+n_iter; */

      prune = tree->a_nodes[prune_idx];

      assert(prune && prune->tax == NO);

      // Select a daughter of prune node
      dir_v1 = dir_v2 = -1;
      for(i=0;i<3;i++) 
        if(prune->v[i] != prune->anc && prune->b[i] != tree->e_root)
          {
            if(dir_v1 < 0) dir_v1 = i;
            else           dir_v2 = i;
          }
      
      u = Uni();
      if(u < 0.5) rnd_dir = dir_v1;
      else        rnd_dir = dir_v2;
          
      prune_daughter = prune->v[rnd_dir];
      cur_regraft_nd = prune->v[rnd_dir == dir_v1 ? dir_v2 : dir_v1];            

      if(prune == tree->n_root)
        {
          if(prune_daughter == prune->v[dir_v1] && prune->v[dir_v2]->tax == YES)
            {
              prune_daughter = prune->v[dir_v2];
              cur_regraft_nd = prune->v[dir_v1];
            }

          if(prune_daughter == prune->v[dir_v2] && prune->v[dir_v1]->tax == YES)
            {
              prune_daughter = prune->v[dir_v1];
              cur_regraft_nd = prune->v[dir_v2];
            }
        }
           
      assert(prune_daughter->anc == prune);

      dir_prune = -1;
      for(i=0;i<3;i++)
        {
          if(prune_daughter->v[i] == prune || prune_daughter->b[i] == tree->e_root)
            {
              dir_prune = i;
              break;
            }
        }
      assert(dir_prune > -1);
            
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


      // Probability of forward move
      regraft_nd_list = NULL;
      a = b = c = d = e = NULL;
      e = cur_regraft_nd;

      if(e->tax == NO && times[e->num] < times[prune_daughter->num])
        {
          for(i=0;i<3;i++)
            {
              if(e->v[i] != e->anc && e->b[i] != tree->e_root)
                {
                  if(a == NULL) a = e->v[i];
                  else          b = e->v[i];
                }
            }
          Push_Bottom_Linked_List(a,&regraft_nd_list,YES);
          Push_Bottom_Linked_List(b,&regraft_nd_list,YES);
        }

      if(prune_daughter->anc != tree->n_root)
        {          
          d = prune_daughter->anc->anc;      

          for(i=0;i<3;i++)
            {
              if(d->v[i] != d->anc && d->b[i] != tree->e_root && d->v[i] != prune)
                {
                  c = d->v[i];
                  break;
                }
            }
          Push_Bottom_Linked_List(c,&regraft_nd_list,YES);
          Push_Bottom_Linked_List(d,&regraft_nd_list,YES);
        }
      
      n_regraft_nd = Linked_List_Len(regraft_nd_list);
      
      if(n_regraft_nd == 0)
        {
          Free_Linked_List(regraft_nd_list);
          continue;
        }

      assert(n_regraft_nd > 0);
      ratio -= log(1./n_regraft_nd);

      // Randomly select one node among potential regraft sites (uniform)
      regraft_idx = Rand_Int(0,n_regraft_nd-1);
      new_regraft_nd = Linked_List_Elem(regraft_idx,regraft_nd_list);      
      Free_Linked_List(regraft_nd_list);

            
      // Time of regraft node and corresponding (partial) Hastings ratio
      t_max = MIN(times[prune_daughter->num],times[cur_regraft_nd->num]);      
      if(prune == tree->n_root) t_min = 10.0*t_max;
      else t_min = times[prune->anc->num];
      ratio += log(1./(t_max - t_min));


      t_max = MIN(times[prune_daughter->num],times[new_regraft_nd->num]);      
      if(new_regraft_nd == tree->n_root) t_min = 10.0*t_max;
      else t_min = times[new_regraft_nd->anc->num];
      ratio -= log(1./(t_max - t_min));
      
      new_t = Uni()*(t_max-t_min) + t_min;

            
      // Age of root node changes when pruned subtree is on one side of that node
      // Change here, not after the prune and regraft move
      /* if(prune == tree->n_root) */
      /*   { */
      /*     if(prune_daughter == tree->n_root->v[1]) */
      /*       times[tree->n_root->num] = times[tree->n_root->v[2]->num]; */
      /*     else if(prune_daughter == tree->n_root->v[2]) */
      /*       times[tree->n_root->num] = times[tree->n_root->v[1]->num]; */
      /*     else assert(false); */
      /*   } */

      /* printf("\n-> prune: %d prune_daughter: %d prune_daughter->v[dir_prune]: %d cur_regraft_nd: %d new_regraft_nd: %d", */
      /*        prune->num, */
      /*        prune_daughter->num, */
      /*        prune_daughter->v[dir_prune]->num, */
      /*        cur_regraft_nd->num, */
      /*        new_regraft_nd->num); */

      /* printf("\n. a:%d b:%d c:%d d:%d e:%d", */
      /*        a ? a->num : -1, */
      /*        b ? b->num : -1, */
      /*        c ? c->num : -1, */
      /*        d ? d->num : -1, */
      /*        e ? e->num : -1); */
      /* fflush(NULL); */
             
      
      // New age
     if(prune == tree->n_root || new_regraft_nd == tree->n_root)
        {
          if(prune == tree->n_root)
            {
              if(prune_daughter == tree->n_root->v[1]) times[tree->n_root->num] = times[tree->n_root->v[2]->num];
              else                                     times[tree->n_root->num] = times[tree->n_root->v[1]->num];
              times[prune_daughter->v[dir_prune]->num] = new_t;
            }
          if(new_regraft_nd == tree->n_root)
            {
              times[prune_daughter->v[dir_prune]->num] = times[tree->n_root->num];
              times[tree->n_root->num] = new_t;
            }
        }
      else
        {
          times[prune->num] = new_t;
        }



      // Prune
      target = residual = NULL;
      Prune_Subtree(prune_daughter->v[dir_prune],
                    prune_daughter,
                    &target,&residual,tree);
      ori_target = target;
      
      // Regraft edge is the one sitting above regraft_nd
      if(new_regraft_nd == tree->n_root->v[1] ||
         new_regraft_nd == tree->n_root->v[2] ||
         new_regraft_nd == tree->n_root) regraft_edge = tree->e_root;
      else
        {
          for(i=0;i<3;i++) if(new_regraft_nd->v[i] == new_regraft_nd->anc) break;
          assert(i!=3);
          regraft_edge = new_regraft_nd->b[i];
        }

      assert(regraft_edge);      
      assert(residual->left != residual->rght);
      assert(regraft_edge->left != prune_daughter->v[dir_prune]);
      assert(regraft_edge->rght != prune_daughter->v[dir_prune]);
      
      // Regraft
      Graft_Subtree(regraft_edge,
                    prune_daughter->v[dir_prune],
                    prune_daughter,
                    residual,
                    new_regraft_nd,tree);
                         
      
      a = b = c = d = e = NULL;
      regraft_nd_list = NULL;
      e = new_regraft_nd;
      if(new_regraft_nd == tree->n_root)
        {
          if(prune_daughter == tree->n_root->v[1])
            e = tree->n_root->v[2];
          else if(prune_daughter == tree->n_root->v[2])
            e = tree->n_root->v[1];
          else assert(false);
        }
      
      if(e->tax == NO && times[e->num] < times[prune_daughter->num])
        {
          for(i=0;i<3;i++)
            {
              if(e->v[i] != e->anc && e->b[i] != tree->e_root)
                {
                  if(a == NULL) a = e->v[i];
                  else          b = e->v[i];
                }
            }
          Push_Bottom_Linked_List(a,&regraft_nd_list,YES);
          Push_Bottom_Linked_List(b,&regraft_nd_list,YES);
        }

      // prune is different from prune_daughter->anc here if prune == tree->n_root
      if(prune_daughter->anc != tree->n_root)
        {          
          d = prune_daughter->anc->anc;      

          for(i=0;i<3;i++)
            {
              if(d->v[i] != d->anc && d->b[i] != tree->e_root && d->v[i] != prune)
                {
                  c = d->v[i];
                  break;
                }
            }
          Push_Bottom_Linked_List(c,&regraft_nd_list,YES);
          Push_Bottom_Linked_List(d,&regraft_nd_list,YES);
        }

      /* printf("\n+ a:%d b:%d c:%d d:%d e:%d", */
      /*        a ? a->num : -1, */
      /*        b ? b->num : -1, */
      /*        c ? c->num : -1, */
      /*        d ? d->num : -1, */
      /*        e ? e->num : -1); */
      /* fflush(NULL); */

      // Number of regraft nodes
      n_regraft_nd = Linked_List_Len(regraft_nd_list);
      assert(n_regraft_nd > 0);
      ratio += log(1./n_regraft_nd);
      Free_Linked_List(regraft_nd_list);
      
      // New age
      /* if(new_regraft_nd == tree->n_root) */
      /*   { */
      /*     if(prune_daughter == tree->n_root->v[1]) */
      /*       times[tree->n_root->v[2]->num] = times[tree->n_root->num]; */
      /*     else if(prune_daughter == tree->n_root->v[2]) */
      /*       times[tree->n_root->v[1]->num] = times[tree->n_root->num]; */
      /*     else assert(false); */

      /*     times[tree->n_root->num] = new_t; */
      /*   } */
      /* else if(prune == tree->n_root) */
      /*   { */
      /*     times[prune_daughter->v[dir_prune]->num] = new_t; */
      /*   } */
      /* else */
      /*   { */
      /*     times[prune->num] = new_t; */
      /*   } */





      if(!TIMES_Check_Node_Height_Ordering(tree))
        {
          printf("\n. prune[%d]->t:%.3f daughter[%d]->t:%.3f prune_anc[%d]->t:%.3f regraft[%d]->t:%.3f regraft_anc[%d]->t:%.3f [effective:%d] t_prior_min/max: [prune:[%.3f %.3f] regraft:[%.3f %.3f]] ",
                 prune->num,
                 times[prune->num],
                 prune_daughter->num,
                 times[prune_daughter->num],
                 prune->anc ? prune->anc->num : -1,
                 prune->anc ? times[prune->anc->num] : -1.,
                 new_regraft_nd->num,
                 times[new_regraft_nd->num],
                 new_regraft_nd->anc ? new_regraft_nd->anc->num : -1,
                 new_regraft_nd->anc ? times[new_regraft_nd->anc->num] : +1.,
                 prune->num,
                 tree->rates->t_prior_min[prune->num],
                 tree->rates->t_prior_max[prune->num],
                 tree->rates->t_prior_min[new_regraft_nd->num],
                 tree->rates->t_prior_max[new_regraft_nd->num]);
          PhyML_Fprintf(stderr,"\n. root: %d %d %d",tree->n_root->num,tree->n_root->v[1]->num,tree->n_root->v[2]->num);
          assert(FALSE);                
        }

      RATES_Update_Cur_Bl(tree);
      DATE_Assign_Primary_Calibration(tree);

      if(tree->eval_glnL == YES) new_lnL_time = TIMES_Lk_Times(NO,tree); 

      if(new_lnL_time > UNLIKELY)
        {
          Set_Both_Sides(NO,tree);
          if(tree->eval_alnL == YES) new_lnL_seq = Lk(NULL,tree);
          if(tree->eval_rlnL == YES) new_lnL_rate = RATES_Lk_Rates(tree);
        }
      
      ratio += (new_lnL_seq - cur_lnL_seq);
      ratio += (new_lnL_rate - cur_lnL_rate);
      ratio += (new_lnL_time - cur_lnL_time);
      
      ratio = exp(ratio);
      alpha = MIN(1.,ratio);

      /* Always accept move */
      if(tree->mcmc->always_yes == YES && new_lnL_time > UNLIKELY) alpha = 1.0;
      
      u = Uni();

      assert(isnan(u) == NO && isinf(fabs(u)) == NO);

      if(u > alpha)
        {
          // Reject
          Prune_Subtree(prune_daughter->v[dir_prune],
                        prune_daughter,
                        &target,&residual,tree);
          assert(residual->left != residual->rght);
          assert(ori_target->left != prune_daughter->v[dir_prune]);
          assert(ori_target->rght != prune_daughter->v[dir_prune]);
          Graft_Subtree(ori_target,
                        prune_daughter->v[dir_prune],
                        prune_daughter,residual,prune == tree->n_root ? tree->n_root : cur_regraft_nd,tree);

          RATES_Reset_Times(tree);
          RATES_Update_Cur_Bl(tree);
          DATE_Assign_Primary_Calibration(tree);
          new_lnL_time = TIMES_Lk_Times(NO,tree); 
          if(Are_Equal(new_lnL_time,cur_lnL_time,1.E-5) == NO)
            {
              PhyML_Printf("\n. new_lnL_time: %f cur_lnL_time: %f",new_lnL_time,cur_lnL_time);
              assert(FALSE);
            }
          
          /* new_lnL_seq = Lk(NULL,tree); */
          /* if(Are_Equal(new_lnL_seq,cur_lnL_seq,1.E-5) == NO) */
          /*   { */
          /*     PhyML_Printf("\n. new: %f cur: %f",new_lnL_seq,cur_lnL_seq); */
          /*     assert(FALSE); */
          /*   } */
          
          
          if(!(tree->rates->c_lnL_times > UNLIKELY))
            {
              printf("\n. time prune: %f",times[prune->num]);
              printf("\n. time prune_daughter: %f",times[prune_daughter->num]);
              printf("\n. prune: %d prune_daughter: %d prune_daughter->v[dir_prune]: %d cur_regraft_nd: %d new_regraft_nd: %d",
                     prune->num,
                     prune_daughter->num,
                     prune_daughter->v[dir_prune]->num,
                     cur_regraft_nd->num,
                     new_regraft_nd->num);
              TIMES_Lk_Times(YES,tree); 
              fflush(NULL);
            }
          assert(tree->rates->c_lnL_times > UNLIKELY);
          
          tree->c_lnL              = cur_lnL_seq;
          tree->rates->c_lnL_times = cur_lnL_time;
          tree->rates->c_lnL_rates = cur_lnL_rate;
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
  mcmc->move_idx = -1;

  mcmc->num_move_nd_r                     = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_br_r                     = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_times                    = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_times_and_rates          = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_times_and_rates_root     = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_root_time                = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_nu                       = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_clock_r                  = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_tree_height              = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_time_slice               = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_subtree_height           = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_kappa                    = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_rr                       = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_tree_rates               = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_subtree_rates            = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_updown_nu_cr             = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_ras                      = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_updown_t_cr              = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_cov_rates                = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_cov_switch               = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_birth_rate               = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_death_rate               = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_birth_death_updown       = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_spr                      = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_spr_weighted             = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_spr_local                = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_spr_root                 = mcmc->n_moves; mcmc->n_moves += 1;
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
  mcmc->num_move_phyrex_spr_local         = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_phyrex_scale_times       = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_phyrex_ldscape_lim       = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_phyrex_sigsq             = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_phyrex_sim               = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_phyrex_traj              = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_phyrex_indel_disk_serial = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_phyrex_sim_plus          = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_phyrex_indel_hit_serial  = mcmc->n_moves; mcmc->n_moves += 1;  
  mcmc->num_move_phyrex_ldsk_given_disk   = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_phyrex_disk_given_ldsk   = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_phyrex_ldsk_and_disk     = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_phyrex_ldsk_multi        = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_phyrex_disk_multi        = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_phyrex_add_remove_jump   = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_clade_change             = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_phyrex_ldsk_tip_to_root  = mcmc->n_moves; mcmc->n_moves += 1;
    
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
  for(i=0;i<mcmc->n_moves;i++) mcmc->move_name[i] = (char *)mCalloc(T_MAX_MCMC_MOVE_NAME,sizeof(char));

  for(i=0;i<mcmc->n_moves;i++) mcmc->adjust_tuning[i] = YES;

  strcpy(mcmc->move_name[mcmc->num_move_br_r],"br_rate");
  strcpy(mcmc->move_name[mcmc->num_move_nd_r],"nd_rate");
  strcpy(mcmc->move_name[mcmc->num_move_times],"times");
  strcpy(mcmc->move_name[mcmc->num_move_times_and_rates],"times_and_rates");
  strcpy(mcmc->move_name[mcmc->num_move_times_and_rates_root],"times_and_rates_root");
  strcpy(mcmc->move_name[mcmc->num_move_root_time],"root_time");
  strcpy(mcmc->move_name[mcmc->num_move_nu],"nu");
  strcpy(mcmc->move_name[mcmc->num_move_clock_r],"clock");
  strcpy(mcmc->move_name[mcmc->num_move_tree_height],"tree_height");
  strcpy(mcmc->move_name[mcmc->num_move_time_slice],"time_slice");
  strcpy(mcmc->move_name[mcmc->num_move_subtree_height],"subtree_height");
  strcpy(mcmc->move_name[mcmc->num_move_kappa],"kappa");
  strcpy(mcmc->move_name[mcmc->num_move_rr],"rr");
  strcpy(mcmc->move_name[mcmc->num_move_spr],"spr");
  strcpy(mcmc->move_name[mcmc->num_move_spr_weighted],"spr_weighted");
  strcpy(mcmc->move_name[mcmc->num_move_spr_local],"spr_local");
  strcpy(mcmc->move_name[mcmc->num_move_spr_root],"spr_root");
  strcpy(mcmc->move_name[mcmc->num_move_tree_rates],"tree_rates");
  strcpy(mcmc->move_name[mcmc->num_move_subtree_rates],"subtree_rates");
  strcpy(mcmc->move_name[mcmc->num_move_updown_nu_cr],"updown_nu_cr");
  strcpy(mcmc->move_name[mcmc->num_move_ras],"ras");  
  strcpy(mcmc->move_name[mcmc->num_move_updown_t_cr],"updown_t_cr");
  strcpy(mcmc->move_name[mcmc->num_move_cov_rates],"cov_rates");  
  strcpy(mcmc->move_name[mcmc->num_move_cov_switch],"cov_switch");
  strcpy(mcmc->move_name[mcmc->num_move_birth_rate],"birth_rate");
  strcpy(mcmc->move_name[mcmc->num_move_death_rate],"death_rate");
  strcpy(mcmc->move_name[mcmc->num_move_birth_death_updown],"birth_death_updown");
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
  strcpy(mcmc->move_name[mcmc->num_move_phyrex_spr_local],"phyrex_spr_local");
  strcpy(mcmc->move_name[mcmc->num_move_phyrex_scale_times],"phyrex_scale_times");
  strcpy(mcmc->move_name[mcmc->num_move_phyrex_ldscape_lim],"phyrex_ldscape_lim");
  strcpy(mcmc->move_name[mcmc->num_move_phyrex_sigsq],"phyrex_sigsq");
  strcpy(mcmc->move_name[mcmc->num_move_phyrex_sim],"phyrex_sim");
  strcpy(mcmc->move_name[mcmc->num_move_phyrex_traj],"phyrex_traj");
  strcpy(mcmc->move_name[mcmc->num_move_phyrex_indel_disk_serial],"phyrex_indel_disk_serial");
  strcpy(mcmc->move_name[mcmc->num_move_phyrex_sim_plus],"phyrex_sim_plus");
  strcpy(mcmc->move_name[mcmc->num_move_phyrex_indel_hit_serial],"phyrex_indel_hit_serial");
  strcpy(mcmc->move_name[mcmc->num_move_phyrex_ldsk_multi],"phyrex_ldsk_multi");
  strcpy(mcmc->move_name[mcmc->num_move_phyrex_disk_multi],"phyrex_disk_multi");
  strcpy(mcmc->move_name[mcmc->num_move_phyrex_ldsk_and_disk],"phyrex_ldsk_and_disk");
  strcpy(mcmc->move_name[mcmc->num_move_phyrex_disk_given_ldsk],"phyrex_disk_given_ldsk");
  strcpy(mcmc->move_name[mcmc->num_move_phyrex_ldsk_given_disk],"phyrex_ldsk_given_disk");
  strcpy(mcmc->move_name[mcmc->num_move_phyrex_add_remove_jump],"phyrex_add_remove_jump");
  strcpy(mcmc->move_name[mcmc->num_move_clade_change],"clade_change");
  strcpy(mcmc->move_name[mcmc->num_move_phyrex_ldsk_tip_to_root],"phyrex_ldsk_tip_to_root");


  for(i=0;i<mcmc->n_moves;i++) mcmc->move_type[i] = -1;
  mcmc->move_type[mcmc->num_move_nd_r] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_br_r] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_times] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_root_time] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_nu] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_clock_r] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_tree_height] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_time_slice] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_subtree_height] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_kappa] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_rr] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_tree_rates] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_subtree_rates] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_updown_nu_cr] = MCMC_MOVE_RANDWALK_NORMAL;
  mcmc->move_type[mcmc->num_move_ras] = MCMC_MOVE_RANDWALK_NORMAL;  
  mcmc->move_type[mcmc->num_move_updown_t_cr] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_cov_rates] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_cov_switch] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_birth_rate] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_death_rate] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_birth_death_updown] = MCMC_MOVE_SCALE_THORNE;
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
  mcmc->move_type[mcmc->num_move_phyrex_indel_hit] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_phyrex_indel_disk] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_phyrex_ldsk_tip_to_root] = MCMC_MOVE_SCALE_THORNE;

  /* We start with small tuning parameter values in order to have inflated ESS
     for clock_r */
  for(i=0;i<mcmc->n_moves;i++)
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
  
  mcmc->move_weight[mcmc->num_move_br_r]                  = 1.0; 
  mcmc->move_weight[mcmc->num_move_nd_r]                  = 0.0;
  mcmc->move_weight[mcmc->num_move_times]                 = 1.0;
  mcmc->move_weight[mcmc->num_move_times_and_rates]       = 5.0;
  mcmc->move_weight[mcmc->num_move_root_time]             = 0.0;
  mcmc->move_weight[mcmc->num_move_clock_r]               = 1.0;
  mcmc->move_weight[mcmc->num_move_tree_height]           = 1.0;
  mcmc->move_weight[mcmc->num_move_time_slice]            = 0.0;
  mcmc->move_weight[mcmc->num_move_subtree_height]        = 0.0;
  mcmc->move_weight[mcmc->num_move_nu]                    = 1.0;
  mcmc->move_weight[mcmc->num_move_kappa]                 = 0.5;
  mcmc->move_weight[mcmc->num_move_rr]                    = 0.5;
  mcmc->move_weight[mcmc->num_move_spr]                   = 3.0;
  mcmc->move_weight[mcmc->num_move_spr_weighted]          = 5.0;
  mcmc->move_weight[mcmc->num_move_spr_local]             = 3.0;
  mcmc->move_weight[mcmc->num_move_spr_root]              = 0.0;
  mcmc->move_weight[mcmc->num_move_tree_rates]            = 1.0;
  mcmc->move_weight[mcmc->num_move_subtree_rates]         = 0.0;
  mcmc->move_weight[mcmc->num_move_updown_nu_cr]          = 0.0;
  mcmc->move_weight[mcmc->num_move_ras]                   = 1.0;
  mcmc->move_weight[mcmc->num_move_updown_t_cr]           = 1.0; /* Does not seem to work well (does not give uniform prior on root height
                                                                    when sampling from prior) */
  mcmc->move_weight[mcmc->num_move_cov_rates]             = 0.0;    
  mcmc->move_weight[mcmc->num_move_cov_switch]            = 0.0;
  mcmc->move_weight[mcmc->num_move_birth_rate]            = 1.0;
  mcmc->move_weight[mcmc->num_move_death_rate]            = 1.0;
  mcmc->move_weight[mcmc->num_move_birth_death_updown]    = 1.0;
  mcmc->move_weight[mcmc->num_move_updown_t_br]           = 0.0;
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
  mcmc->move_weight[mcmc->num_move_clade_change]          = 1.0;

# if defined (PHYREX)

  mcmc->move_weight[mcmc->num_move_phyrex_lbda]                  = 1.0;
  mcmc->move_weight[mcmc->num_move_phyrex_mu]                    = 1.0;
  mcmc->move_weight[mcmc->num_move_phyrex_rad]                   = 1.0;
  mcmc->move_weight[mcmc->num_move_phyrex_sigsq]                 = 0.0;
  mcmc->move_weight[mcmc->num_move_phyrex_indel_disk]            = 1.0;
  mcmc->move_weight[mcmc->num_move_phyrex_indel_hit]             = 1.0;
  mcmc->move_weight[mcmc->num_move_phyrex_move_disk_ud]          = 2.0;
  mcmc->move_weight[mcmc->num_move_phyrex_swap_disk]             = 2.0;
  mcmc->move_weight[mcmc->num_move_phyrex_spr]                   = 0.5;
  mcmc->move_weight[mcmc->num_move_phyrex_spr_local]             = 1.0;
  mcmc->move_weight[mcmc->num_move_phyrex_scale_times]           = 2.0;
  mcmc->move_weight[mcmc->num_move_phyrex_ldscape_lim]           = 0.0;
  mcmc->move_weight[mcmc->num_move_phyrex_sim]                   = 0.0;
  mcmc->move_weight[mcmc->num_move_phyrex_traj]                  = 2.0;
  mcmc->move_weight[mcmc->num_move_phyrex_sim_plus]              = 0.0;
  mcmc->move_weight[mcmc->num_move_phyrex_indel_disk_serial]     = 2.0;
  mcmc->move_weight[mcmc->num_move_phyrex_indel_hit_serial]      = 2.0;

  mcmc->move_weight[mcmc->num_move_phyrex_ldsk_given_disk]       = 1.0;
  mcmc->move_weight[mcmc->num_move_phyrex_disk_given_ldsk]       = 1.0;
  mcmc->move_weight[mcmc->num_move_phyrex_ldsk_and_disk]         = 0.5;
  mcmc->move_weight[mcmc->num_move_phyrex_ldsk_multi]            = 1.0;
  mcmc->move_weight[mcmc->num_move_phyrex_disk_multi]            = 1.0;
  mcmc->move_weight[mcmc->num_move_phyrex_add_remove_jump]       = 3.0;
  mcmc->move_weight[mcmc->num_move_phyrex_ldsk_tip_to_root]      = 1.0;

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
  mcmc->move_weight[mcmc->num_move_phyrex_spr_local]             = 0.0;
  mcmc->move_weight[mcmc->num_move_phyrex_scale_times]           = 0.0;
  mcmc->move_weight[mcmc->num_move_phyrex_ldscape_lim]           = 0.0;
  mcmc->move_weight[mcmc->num_move_phyrex_sim]                   = 0.0;
  mcmc->move_weight[mcmc->num_move_phyrex_traj]                  = 0.0;
  mcmc->move_weight[mcmc->num_move_phyrex_sim_plus]              = 0.0;
  mcmc->move_weight[mcmc->num_move_phyrex_indel_disk_serial]     = 0.0;
  mcmc->move_weight[mcmc->num_move_phyrex_indel_hit_serial]      = 0.0;

  mcmc->move_weight[mcmc->num_move_phyrex_ldsk_given_disk]       = 0.0;
  mcmc->move_weight[mcmc->num_move_phyrex_disk_given_ldsk]       = 0.0;
  mcmc->move_weight[mcmc->num_move_phyrex_ldsk_and_disk]         = 0.0;
  mcmc->move_weight[mcmc->num_move_phyrex_ldsk_multi]            = 0.0;
  mcmc->move_weight[mcmc->num_move_phyrex_disk_multi]            = 0.0;
  mcmc->move_weight[mcmc->num_move_phyrex_add_remove_jump]       = 0.0;
  mcmc->move_weight[mcmc->num_move_phyrex_ldsk_tip_to_root]      = 0.0;
#endif

  sum = 0.0;
  for(i=0;i<mcmc->n_moves;i++) sum += mcmc->move_weight[i];
  for(i=0;i<mcmc->n_moves;i++) mcmc->move_weight[i] /= sum;
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
  
  /* for(i=0;i<tree->n_otu-1;i++) */
  /*   mcmc->lag_param_val[mcmc->num_move_nd_t+i] = tree->rates->nd_t[tree->n_otu+i]; */

  /* For(i,2*tree->n_otu-1) */
  /*   mcmc->lag_param_val[mcmc->num_move_nd_r+i] = tree->rates->nd_r[i]; */

  /* for(i=0;i<mcmc->n_moves;i++) tree->mcmc->new_param_val[i] = tree->mcmc->lag_param_val[i]; */
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
  else for(i=0;i<3;i++) if(d->v[i] == a) { b = d->b[i]; break; }
  
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
	  for(i=0;i<3;i++)
	    if(d->v[i] != a && d->b[i] != tree->e_root)
	      {
		if(tree->io->lk_approx == EXACT) Update_Partial_Lk(tree,d->b[i],d);
		MCMC_Slice_One_Rate(d,d->v[i],YES,tree);
	      }
	}
      if(tree->io->lk_approx == EXACT) Update_Partial_Lk(tree,b,d);
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
	*new   = (*cur) * exp(tune*(u-.5));
	*loghr = log((*new)/(*cur));
	break;
      }

    case MCMC_MOVE_SCALE_GAMMA:
      {
	phydbl r;
	*new = (*cur) * Rgamma(1./tune,tune);
	r = (*new)/(*cur);
	*loghr = -log(r) + log(Dgamma(1./r,1./tune,tune)/Dgamma(r,1./tune,tune));
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
      PhyML_Fprintf(stderr,"\n. Wrong file format.");
      assert(FALSE);
    }
    
  if(fgets(token,sizemax,in_fp) == NULL) // Skip second
    {
      PhyML_Fprintf(stderr,"\n. Wrong file format.");
      assert(FALSE);
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

  for(i=0;i<tree->n_otu-1;i++)
    {
      v=fscanf(in_fp,"%lf\t",&val); // Node heights
      tree->rates->nd_t[i+tree->n_otu] = val;
    }

  For(i,2*tree->n_otu-2)
    {
      v=fscanf(in_fp,"%lf\t",&val); // Edge average rates
      tree->rates->br_r[i] = log(val);
    }
  
  v++;

  Free(token);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

#ifdef PHYREX
void MCMC_PHYREX_Lbda(t_tree *tree)
{
  /* MCMC_Single_Param_Generic(&(tree->mmod->lbda), */
  /*                           tree->mmod->min_lbda, */
  /*                           tree->mmod->max_lbda, */
  /*                           tree->mcmc->num_move_phyrex_lbda, */
  /*                           NULL,&(tree->mmod->c_lnL), */
  /*                           NULL,(tree->eval_glnL == YES) ? PHYREX_Wrap_Lk : NULL, */
  /*                           tree->mcmc->move_type[tree->mcmc->num_move_phyrex_lbda], */
  /*                           NO,NULL,tree,NULL); */


  phydbl u,alpha,ratio;
  phydbl cur_glnL, new_glnL, hr;
  phydbl cur_lbda, new_lbda;
  phydbl K;

  tree->mcmc->run_move[tree->mcmc->num_move_phyrex_lbda]++;

  new_glnL  = tree->mmod->c_lnL;
  cur_glnL  = tree->mmod->c_lnL;
  hr        = 0.0;
  ratio     = 0.0;
  cur_lbda  = tree->mmod->lbda;
  new_lbda  = tree->mmod->lbda;
  K         = tree->mcmc->tune_move[tree->mcmc->num_move_phyrex_lbda];

  MCMC_Make_Move(&cur_lbda,&new_lbda,tree->mmod->min_lbda,tree->mmod->max_lbda,&ratio,K,tree->mcmc->move_type[tree->mcmc->num_move_phyrex_lbda]);
  
  if(new_lbda < tree->mmod->max_lbda && new_lbda > tree->mmod->min_lbda)
    {
      tree->mmod->lbda = new_lbda;
      
      if(tree->eval_glnL == YES) new_glnL = PHYREX_Lk(tree);

      ratio += (new_glnL - cur_glnL);
      ratio += hr;
  
      ratio = exp(ratio);
      alpha = MIN(1.,ratio);


      /* PhyML_Printf("\n. Lbda: %12f new_glnL: %12f cur_glnL: %12f",tree->mmod->lbda,new_glnL,cur_glnL); */
      
      /* Always accept move */
      if(tree->mcmc->always_yes == YES && new_glnL > UNLIKELY) alpha = 1.0;
      
      u = Uni();
      
      if(u > alpha) /* Reject */
        {
          tree->mmod->lbda = cur_lbda;
          tree->mmod->c_lnL = cur_glnL;
        }
      else
        {
          tree->mcmc->acc_move[tree->mcmc->num_move_phyrex_lbda]++;
        }
    }
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
                            NULL,(tree->eval_glnL == YES) ? PHYREX_Wrap_Lk : NULL,
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
                            NULL,(tree->eval_glnL == YES) ? PHYREX_Wrap_Lk : NULL,
                            tree->mcmc->move_type[tree->mcmc->num_move_phyrex_rad],
                            NO,NULL,tree,NULL);

  /* phydbl u,alpha,ratio; */
  /* phydbl cur_glnL, new_glnL, hr; */
  /* phydbl cur_rad, new_rad; */
  /* phydbl K; */

  /* tree->mcmc->run_move[tree->mcmc->num_move_phyrex_rad]++; */

  /* new_glnL  = tree->mmod->c_lnL; */
  /* cur_glnL  = tree->mmod->c_lnL; */
  /* hr        = 0.0; */
  /* ratio     = 0.0; */
  /* cur_rad  = tree->mmod->rad; */
  /* new_rad  = tree->mmod->rad; */
  /* K         = tree->mcmc->tune_move[tree->mcmc->num_move_phyrex_rad]; */

  /* MCMC_Make_Move(&cur_rad,&new_rad,tree->mmod->min_rad,tree->mmod->max_rad,&ratio,K,tree->mcmc->move_type[tree->mcmc->num_move_phyrex_rad]); */
  
  /* if(new_rad < tree->mmod->max_rad && new_rad > tree->mmod->min_rad) */
  /*   { */
  /*     tree->mmod->rad = new_rad; */
      
  /*     if(tree->eval_glnL == YES) new_glnL = PHYREX_Lk(tree); */

  /*     ratio += (new_glnL - cur_glnL); */
  /*     ratio += hr; */
  
  /*     ratio = exp(ratio); */
  /*     alpha = MIN(1.,ratio); */


  /*     PhyML_Printf("\n. new_rad: %12f cur_rad: %12f new_glnL: %12f cur_glnL: %12f",new_rad,cur_rad,new_glnL,cur_glnL); */
      
  /*     /\* Always accept move *\/ */
  /*     if(tree->mcmc->always_yes == YES && new_glnL > UNLIKELY) alpha = 1.0; */
      
  /*     u = Uni(); */
      
  /*     if(u > alpha) /\* Reject *\/ */
  /*       { */
  /*         tree->mmod->rad = cur_rad; */
  /*         tree->mmod->c_lnL = cur_glnL; */
  /*       } */
  /*     else */
  /*       { */
  /*         tree->mcmc->acc_move[tree->mcmc->num_move_phyrex_rad]++; */
  /*       } */
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
                            NULL,(tree->eval_glnL == YES) ? PHYREX_Wrap_Lk : NULL,
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
  phydbl hr;
  phydbl cur_lbda,new_lbda;
  phydbl cur_mu,new_mu;
  phydbl cur_rad,new_rad;
  t_dsk *disk;

  hr       = 0.0;  
  new_lbda = tree->mmod->lbda;
  cur_lbda = tree->mmod->lbda;
  new_mu   = tree->mmod->mu;
  cur_mu   = tree->mmod->mu;
  new_rad  = tree->mmod->rad;
  cur_rad  = tree->mmod->rad;

  new_lbda = cur_lbda * exp(1.0*(Uni()-.5));
  hr += log(new_lbda/cur_lbda);

  new_mu = cur_mu * exp(0.5*(Uni()-.5));
  hr += log(new_mu/cur_mu);

  new_rad = cur_rad * exp(0.5*(Uni()-.5));
  hr += log(new_rad/cur_rad);

  if(new_rad > tree->mmod->max_rad || new_rad < tree->mmod->min_rad)     return; 
  if(new_mu > tree->mmod->max_mu || new_mu < tree->mmod->min_mu)         return; 
  if(new_lbda > tree->mmod->max_lbda || new_lbda < tree->mmod->min_lbda) return; 

  tree->mmod->lbda = new_lbda;
  tree->mmod->mu   = new_mu;
  tree->mmod->rad  = new_rad;

  disk = tree->young_disk->prev;
  n_disks_cur = 0;
  do
    {
      if(disk->ldsk == NULL && disk->age_fixed == NO) n_disks_cur++;
      disk = disk->prev;
    }
  while(disk);
  
  n_disks_new = (int)Rpois(n_disks_cur+SMALL);
  hr += Dpois(n_disks_cur,n_disks_new+SMALL,YES);
  hr -= Dpois(n_disks_new,n_disks_cur+SMALL,YES);

  if(n_disks_new < n_disks_cur)
    {
      MCMC_PHYREX_Delete_Disk(hr, n_disks_cur - n_disks_new , cur_lbda, cur_mu, cur_rad, tree);
    }
  else
    {
      MCMC_PHYREX_Insert_Disk(hr, n_disks_new - n_disks_cur,  cur_lbda, cur_mu, cur_rad, tree);
    }

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

  disk = tree->young_disk->prev;

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
      if(!disk->ldsk && disk->age_fixed == NO)
        {
          if(!n_valid_disks) valid_disks = (t_dsk **)mCalloc(block,sizeof(t_dsk *));
          else if(!(n_valid_disks%block)) valid_disks = (t_dsk **)mRealloc(valid_disks,n_valid_disks+block,sizeof(t_dsk *));
          valid_disks[n_valid_disks] = disk;
          n_valid_disks++;
        }
      disk = disk->prev;
    }
  while(disk && disk->prev);


  if(!n_valid_disks || (n_valid_disks - n_delete_disks < 0))
    {
      tree->mmod->lbda = cur_lbda;
      tree->mmod->rad = cur_rad;
      tree->mmod->mu = cur_mu;
      return; 
    }

  target_disk = (t_dsk **)mCalloc(n_delete_disks,sizeof(t_dsk *));

  permut = Permutate(n_valid_disks);

  for(j=0;j<n_delete_disks;j++)
    {
      /* Uniform selection of a disk where no coalescent nor 'hit' occurred */
      target_disk[j] = valid_disks[permut[j]];
      PHYREX_Remove_Disk(target_disk[j]);  
      for(i=0;i<tree->mmod->n_dim;i++) hr += log(1./(tree->mmod->lim_up->lonlat[i]-tree->mmod->lim_do->lonlat[i]));
    }

  T = PHYREX_Tree_Height(tree);
  T = fabs(T-tree->young_disk->time);

  
  // Pr(n -> n-k) i.e, denominator
  hr += LnChoose(n_valid_disks,n_delete_disks);
  
  // Pr(n-k -> n) i.e., numerator
  hr += n_delete_disks * log(1./T);
  hr += LnFact(n_delete_disks);

  if(tree->eval_glnL == YES) new_glnL = PHYREX_Lk(tree);

  ratio += (new_glnL - cur_glnL);
  ratio += hr;
  
  ratio = exp(ratio);
  alpha = MIN(1.,ratio);

  /* Always accept move */
  if(tree->mcmc->always_yes == YES && new_glnL > UNLIKELY) alpha = 1.0;

  u = Uni();
  
  /* PhyML_Printf("\n- Delete new_glnL: %f [%f] hr: %f u:%f alpha: %f n_delete: %d [%d %d %f]", */
  /*              new_glnL,cur_glnL,hr,u,alpha,n_delete_disks, */
  /*              tree->mcmc->acc_move[tree->mcmc->num_move_phyrex_indel_disk], */
  /*              tree->mcmc->run_move[tree->mcmc->num_move_phyrex_indel_disk], */
  /*              tree->mcmc->acc_rate[tree->mcmc->num_move_phyrex_indel_disk]); */


  assert(isnan(u) == NO && isinf(fabs(u)) == NO);

  if(u > alpha) /* Reject */
    {
      tree->mmod->lbda = cur_lbda;
      tree->mmod->rad = cur_rad;
      tree->mmod->mu = cur_mu;

      for(j=0;j<n_delete_disks;j++) PHYREX_Insert_Disk(target_disk[j],tree);      
      PHYREX_Update_Lindisk_List(tree);
      tree->mmod->c_lnL = cur_glnL;
    }
  else
    {
      for(j=0;j<n_delete_disks;j++) Free_Disk(target_disk[j]);
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

  disk     = NULL;
  new_glnL = tree->mmod->c_lnL;
  cur_glnL = tree->mmod->c_lnL;
  ratio    = 0.0;

  if(tree->young_disk->next) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);

  T = PHYREX_Tree_Height(tree);
  
  disk = tree->young_disk->prev;
  n_valid_disks = 0;
  do
    {
      if(!disk->ldsk) n_valid_disks++;
      disk = disk->prev;
    }
  while(disk && disk->prev);

  target_disk = (t_dsk **)mCalloc(n_insert_disks,sizeof(t_dsk *));
  new_disk    = (t_dsk **)mCalloc(n_insert_disks,sizeof(t_dsk *));


  for(j=0;j<n_insert_disks;j++)
    {
      t = Uni()*(tree->young_disk->time - T) + T;
      disk = tree->young_disk->prev;
      while(disk && disk->time > t) disk = disk->prev;
      assert(disk->next);
      target_disk[j] = disk->next;
      
      new_disk[j] = PHYREX_Make_Disk_Event(tree->mmod->n_dim,tree->n_otu);
      PHYREX_Init_Disk_Event(new_disk[j],tree->mmod->n_dim,tree->mmod);
      new_disk[j]->time = t;
        
      PHYREX_Insert_Disk(new_disk[j],tree);
      
      for(i=0;i<tree->mmod->n_dim;i++) new_disk[j]->centr->lonlat[i] = Uni()*(tree->mmod->lim_up->lonlat[i]-tree->mmod->lim_do->lonlat[i])+tree->mmod->lim_do->lonlat[i];

      for(i=0;i<tree->mmod->n_dim;i++) hr -= log(1./(tree->mmod->lim_up->lonlat[i]-tree->mmod->lim_do->lonlat[i]));
    }

  /* // Pr(n+k -> n) i.e, numerator */
  hr -= LnChoose(n_valid_disks+n_insert_disks,n_insert_disks);
  
  /* // Pr(n -> n+k) i.e., denominator */
  T = fabs(T-tree->young_disk->time);
  hr -= n_insert_disks * log(1./T);
  hr -= LnFact(n_insert_disks);

  
  if(tree->eval_glnL == YES) new_glnL = PHYREX_Lk(tree);  

  ratio += (new_glnL - cur_glnL);
  ratio += hr;
  
  ratio = exp(ratio);
  alpha = MIN(1.,ratio);
  
  /* Always accept move */
  if(tree->mcmc->always_yes == YES && new_glnL > UNLIKELY) alpha = 1.0;

  u = Uni();
  
  /* PhyML_Printf("\n+ Insert new_glnL: %f [%f] hr: %f u: %f alpha: %f n_insert: %d [%d %d %d %d %f]", */
  /*              new_glnL,cur_glnL,hr,u,alpha,n_insert_disks, */
  /*              tree->mcmc->acc_move[tree->mcmc->num_move_phyrex_indel_disk], */
  /*              tree->mcmc->run_move[tree->mcmc->num_move_phyrex_indel_disk], */
  /*              tree->mcmc->prev_acc_move[tree->mcmc->num_move_phyrex_indel_disk], */
  /*              tree->mcmc->prev_run_move[tree->mcmc->num_move_phyrex_indel_disk], */
  /*              tree->mcmc->acc_rate[tree->mcmc->num_move_phyrex_indel_disk]); */

  assert(isnan(u) == NO && isinf(fabs(u)) == NO);

  if(u > alpha) /* Reject */
    {
      tree->mmod->lbda = cur_lbda;
      tree->mmod->rad = cur_rad;
      tree->mmod->mu = cur_mu;
      
      for(j=0;j<n_insert_disks;j++)
        {
          PHYREX_Remove_Disk(new_disk[j]);
          Free_Disk(new_disk[j]);
        }

      PHYREX_Update_Lindisk_List(tree);
      tree->mmod->c_lnL = cur_glnL;
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

/* Update time of disk */
#ifdef PHYREX

void MCMC_PHYREX_Move_Disk_Updown(t_tree *tree)
{
  phydbl u,alpha,ratio;
  phydbl cur_glnL, new_glnL, hr;
  phydbl cur_alnL, new_alnL;
  phydbl cur_rlnL, new_rlnL;
  phydbl *ori_time,new_time,cur_time;
  phydbl max, min;
  t_dsk  *disk,**target_disk,**all_disks;
  int i,block,n_all_disks,n_move_disks,*permut,update_alnL;

  disk          = NULL;
  block         = 100;
  all_disks     = NULL;
  
  disk = tree->young_disk->prev;
  n_all_disks = 0;
  do
    {
      if((disk->age_fixed) == NO && !(disk->ldsk && disk->ldsk->n_next > 1 && tree->mod->s_opt->opt_bl == NO))
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
  
  n_move_disks = (int)(1+0.2*n_all_disks);
  
  target_disk = (t_dsk **)mCalloc(n_all_disks,sizeof(t_dsk *));
  ori_time    = (phydbl *)mCalloc(n_all_disks,sizeof(phydbl));

  permut = Permutate(n_all_disks);
  
  for(i=0;i<n_move_disks;i++)
    {
      tree->mcmc->run_move[tree->mcmc->num_move_phyrex_move_disk_ud]++;

      cur_glnL = tree->mmod->c_lnL;
      new_glnL = cur_glnL;
      new_alnL = tree->c_lnL;
      cur_alnL = tree->c_lnL;
      new_rlnL = tree->rates->c_lnL_rates;
      cur_rlnL = tree->rates->c_lnL_rates;
      hr       = 0.0;
      ratio    = 0.0;

      target_disk[i] = all_disks[permut[i]];


      update_alnL = NO;
      if(target_disk[i]->ldsk && target_disk[i]->ldsk->n_next > 1) update_alnL = YES;
      
      ori_time[i] = target_disk[i]->time;
      cur_time = target_disk[i]->time;
      
      
      if(target_disk[i]->prev)
        {
          max = target_disk[i]->next->time;
          min = target_disk[i]->prev->time;
          new_time = Uni()*(max - min) + min;
          /* PhyML_Printf("\n. disk: %s max: %f min: %f curr_time: %f new_time: %f",target_disk[i]->id,max,min,target_disk[i]->time,new_time); */
        }
      else
        {
          /* phydbl new_plusmin, cur_plusmin; */
          phydbl K;

          K = 0.1;
                    
          max = target_disk[i]->next->time;

          /* cur_plusmin = fabs(ori_time[i] - max); */
          /* new_plusmin = Rexp(1./cur_plusmin); */

          /* new_time = max - new_plusmin; */
          new_time = max - Rexp(K*tree->mmod->lbda);
                    
          /* hr += log(Dexp(cur_plusmin,1./new_plusmin)); */
          /* hr -= log(Dexp(new_plusmin,1./cur_plusmin)); */
          
          hr -= -(max-new_time)*K*tree->mmod->lbda;
          hr += -(max-cur_time)*K*tree->mmod->lbda;
          /* PhyML_Printf("\n* disk: %s max: %f cur_time: %f new_time: %f",target_disk[i]->id,max,target_disk[i]->time,new_time); */


          new_glnL += -tree->mmod->lbda * fabs(new_time);
          new_glnL -= -tree->mmod->lbda * fabs(cur_time);
          
          tree->mmod->c_lnL = new_glnL;
        }
      

      target_disk[i]->time = new_time;
    
      if(tree->eval_alnL == YES && update_alnL == YES) new_alnL = Lk(NULL,tree);
      if(tree->eval_rlnL == YES && update_alnL == YES) new_rlnL = RATES_Lk_Rates(tree);
      
      ratio += (new_alnL - cur_alnL);
      ratio += (new_glnL - cur_glnL);
      ratio += (new_rlnL - cur_rlnL);
      ratio += hr;
      
      ratio = exp(ratio);
      alpha = MIN(1.,ratio);
      
      /* Always accept move */
      if(tree->mcmc->always_yes == YES && new_glnL > UNLIKELY) alpha = 1.0;
      
      u = Uni();
      
      /* if(target_disk[i]->prev == NULL) PhyML_Printf("\n- Move disk new_glnL: %12f [%12f] hr: %12f u:%6f alpha: %12f new_time: %12f cur_time: %12f max: %12f",new_glnL,cur_glnL,hr,u,alpha,new_time,cur_time,max); */
      
      assert(isnan(u) == NO && isinf(fabs(u)) == NO);
      
      if(u > alpha) /* Reject */
        {
          /* printf("- Reject"); */
          
          target_disk[i]->time = ori_time[i];
                    
          tree->mmod->c_lnL        = cur_glnL;
          tree->c_lnL              = cur_alnL;
          tree->rates->c_lnL_rates = cur_rlnL;
        }
      else
        {
          tree->mcmc->acc_move[tree->mcmc->num_move_phyrex_move_disk_ud]++;
          /* printf("- Accept"); */
        }
    }
       
  Free(target_disk);
  Free(all_disks);
  Free(ori_time);
  Free(permut);
}
#endif

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

#ifdef PHYREX
void MCMC_PHYREX_Scale_Times(t_tree *tree)
{
  phydbl u,alpha,ratio;
  phydbl cur_height,new_height;
  phydbl cur_glnL, new_glnL, hr;
  phydbl cur_alnL, new_alnL;
  phydbl cur_rlnL, new_rlnL;
  phydbl scale_fact_times;
  int n_disks;
  t_dsk  *start_disk;
  phydbl K;
  
  if(tree->mod->s_opt->opt_bl == NO) return;

  cur_alnL   = tree->c_lnL;
  new_alnL   = tree->c_lnL;
  new_glnL   = tree->mmod->c_lnL;
  cur_glnL   = tree->mmod->c_lnL;
  new_rlnL   = tree->rates->c_lnL_rates;
  cur_rlnL   = tree->rates->c_lnL_rates;
  hr         = 0.0;
  ratio      = 0.0;
  cur_height = fabs(PHYREX_Tree_Height(tree));
  K          = tree->mcmc->tune_move[tree->mcmc->num_move_phyrex_scale_times];
  
  tree->mcmc->run_move[tree->mcmc->num_move_phyrex_scale_times]++;
  
  u = Uni();
  K = 2.;
  scale_fact_times = exp(K*(u-.5));

  start_disk = tree->old_samp_disk;
  assert(start_disk);
  
  n_disks = PHYREX_Scale_All(scale_fact_times,start_disk,tree);
  if(!PHYREX_Check_Struct(tree))
    {
      PHYREX_Scale_All(1./scale_fact_times,start_disk,tree);
      return;
    }

    
  /* The Hastings ratio involves (n_disk-2) when considering a uniform distrib
     for the multiplier, which is not the case here.
  */
  hr += (n_disks)*log(scale_fact_times);
  
  if(tree->eval_glnL == YES)
    {
      new_height = cur_height * scale_fact_times + fabs(start_disk->time) * (1. - scale_fact_times);
      
      new_glnL -= -tree->mmod->lbda * cur_height;
      new_glnL += -tree->mmod->lbda * new_height;
      tree->mmod->c_lnL = new_glnL;

      /* PHYREX_Lk(tree); */
      /* new_glnL = tree->mmod->c_lnL; */
    }
  
  if(tree->eval_alnL == YES) new_alnL = Lk(NULL,tree);
  if(tree->eval_rlnL == YES) new_rlnL = RATES_Lk_Rates(tree);
  
  ratio += (new_alnL - cur_alnL);
  ratio += (new_glnL - cur_glnL);
  ratio += (new_rlnL - cur_rlnL);

  ratio += hr;
  
  ratio = exp(ratio);
  alpha = MIN(1.,ratio);
  
  /* Always accept move */
  if(tree->mcmc->always_yes == YES && new_glnL > UNLIKELY) alpha = 1.0;

  u = Uni();
  
  /* PhyML_Printf("\n. Scale times hr: %f new_glnL: %f cur_glnL: %f new_alnL: %f cur_alnL: %f new_rlnL: %f cur_rlnL: %f ratio: %f scale: %f", */
  /*              hr, */
  /*              new_glnL,cur_glnL, */
  /*              new_alnL,cur_alnL, */
  /*              new_rlnL,cur_rlnL, */
  /*              ratio,scale_fact_times); */

  assert(isnan(u) == NO && isinf(fabs(u)) == NO);

  
  if(u > alpha) /* Reject */
    {
      PHYREX_Scale_All(1./scale_fact_times,start_disk,tree);
      PHYREX_Update_Lindisk_List(tree);
      tree->mmod->c_lnL        = cur_glnL;
      tree->c_lnL              = cur_alnL;        
      tree->rates->c_lnL_rates = cur_rlnL;
    }
  else
    {
      tree->mcmc->acc_move[tree->mcmc->num_move_phyrex_scale_times]++;
    }  
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
  phydbl cur_rlnL, new_rlnL;
  phydbl hr,t,t_min,t_max,ori_time;
  t_dsk  *disk,*target_disk,*ori_disk_old,*ori_disk_young,**valid_disks;
  t_ldsk *ldsk_next;
  int i,j,block,n_valid_disks;

  valid_disks = NULL;
  block       = 100;

  disk = tree->young_disk->prev;
  if(!disk->prev) return;

  n_valid_disks = 0;
  do
    {
      /* Record disk with a lineage displacement or coalescent that is not the root disk and not a disk with sampled taxa on it */
      if(disk && disk->prev && disk->ldsk && disk->ldsk->n_next >= 1 && disk->age_fixed == NO)
        {
          if(!(disk->ldsk->n_next > 1 && tree->mod->s_opt->opt_bl == NO))
            {
              if(!n_valid_disks) valid_disks = (t_dsk **)mCalloc(block,sizeof(t_dsk *));
              else if(!(n_valid_disks%block)) valid_disks = (t_dsk **)mRealloc(valid_disks,n_valid_disks+block,sizeof(t_dsk *));
              valid_disks[n_valid_disks] = disk;
              n_valid_disks++;
            }
        }
      disk = disk->prev;
    }
  while(disk && disk->prev);
  if(n_valid_disks < 2) return;
  
  for(j=0;j<MIN(10,(int)(1.+0.1*tree->n_otu));++j)
    {
      disk           = NULL;
      target_disk    = NULL;
      ori_disk_old   = NULL;
      ori_disk_young = NULL;
      ldsk_next      = NULL;
      new_glnL       = tree->mmod->c_lnL;
      cur_glnL       = tree->mmod->c_lnL;
      new_alnL       = tree->c_lnL;
      cur_alnL       = tree->c_lnL;
      new_rlnL       = tree->rates->c_lnL_rates;
      cur_rlnL       = tree->rates->c_lnL_rates;
      hr             = 0.0;
      ratio          = 0.0;
      
          
      /* Uniform selection of a valid disk */
      i = Rand_Int(0,n_valid_disks-1);
      target_disk = valid_disks[i];
      
      ori_disk_old   = target_disk->prev;
      ori_disk_young = target_disk->next;
      ori_time       = target_disk->time;
      
      t_min = target_disk->ldsk->prev->disk->time;
      t_max = +INFINITY;
      for(i=0;i<target_disk->ldsk->n_next;++i)
        if(target_disk->ldsk->next[i]->disk->time < t_max) 
          {
            t_max = target_disk->ldsk->next[i]->disk->time;
            ldsk_next = target_disk->ldsk->next[i];
          }
      
      t_max = t_max - SMALL;
      t_min = t_min + SMALL;

      if(tree->eval_glnL == YES)
        {
          new_glnL -= PHYREX_Lk_Range(ldsk_next->disk,target_disk->ldsk->prev->disk,tree);
        }

      if(!(t_max > t_min))
        {
          PhyML_Printf("\n. t_max: %g t_min: %g",t_max,t_min);
          assert(t_max > t_min);
        }
      t = Uni()*(t_max - t_min) + t_min;
      target_disk->time = t;
      
      PHYREX_Remove_Disk(target_disk);      
      PHYREX_Insert_Disk(target_disk,tree);
      
      if(tree->eval_glnL == YES)
        {
          new_glnL += PHYREX_Lk_Range(ldsk_next->disk,target_disk->ldsk->prev->disk,tree);
          tree->mmod->c_lnL = new_glnL;
          /* new_glnL = PHYREX_Lk(tree); */
        }

      if(tree->eval_alnL == YES) new_alnL = Lk(NULL,tree);
      if(tree->eval_rlnL == YES) new_rlnL = RATES_Lk_Rates(tree);
      
      ratio += (new_alnL - cur_alnL);
      ratio += (new_glnL - cur_glnL);
      ratio += (new_rlnL - cur_rlnL);

      ratio += hr;
            
      ratio = exp(ratio);
      alpha = MIN(1.,ratio);
      
      /* Always accept move */
      if(tree->mcmc->always_yes == YES && new_glnL > UNLIKELY) alpha = 1.0;
      
      u = Uni();
      
      /* printf("\n- Swap new_glnL: %f [%f] new_alnL: %f [%f] hr: %f u:%f alpha: %f",new_glnL,cur_glnL,new_alnL,cur_alnL,hr,u,alpha); */
      
      assert(isnan(u) == NO && isinf(fabs(u)) == NO);

      if(u > alpha) /* Reject */
        {
          /* printf("\n- Reject %f %f",target_disk->time,ori_time); */
          
          PHYREX_Remove_Disk(target_disk);
          
          target_disk->time = ori_time;
          target_disk->prev = ori_disk_old;
          target_disk->next = ori_disk_young;

          PHYREX_Insert_Disk(target_disk,tree);          
          PHYREX_Update_Lindisk_List_Range(ldsk_next->disk,target_disk->ldsk->prev->disk,tree);

          tree->mmod->c_lnL        = cur_glnL;
          tree->c_lnL              = cur_alnL;        
          tree->rates->c_lnL_rates = cur_rlnL;        
        }
      else
        {
          /* printf("\n- Accept %f %f",target_disk->time,ori_time); */
        }
    }

  Free(valid_disks);
  
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
  
  /* new_rad = cur_rad * exp(0.5*(Uni()-.5)); */
  /* hr += log(new_rad/cur_rad); */

  /* new_mu = cur_mu * exp(0.5*(Uni()-.5)); */
  /* hr += log(new_mu/cur_mu); */

  if(new_rad > tree->mmod->max_rad || new_rad < tree->mmod->min_rad) return; 
  if(new_mu > tree->mmod->max_mu || new_mu < tree->mmod->min_mu)     return; 

  tree->mmod->rad = new_rad;
  tree->mmod->mu  = new_mu;


  disk = tree->young_disk->prev;
  n_disks_cur = 0;
  do
    {
      if(disk->ldsk != NULL && disk->ldsk->n_next == 1) n_disks_cur++;
      disk = disk->prev;
    }
  while(disk);
  
  n_disks_new = (int)Rpois(n_disks_cur+0.01);
  hr += Dpois(n_disks_cur,n_disks_new+0.01,YES);
  hr -= Dpois(n_disks_new,n_disks_cur+0.01,YES);

  if(n_disks_new < n_disks_cur)
    {
      MCMC_PHYREX_Delete_Hit(hr, n_disks_cur - n_disks_new, cur_rad, cur_mu, tree);
    }
  else
    {
      MCMC_PHYREX_Insert_Hit(hr, n_disks_new - n_disks_cur, cur_rad, cur_mu, tree);
    }
  
}
#endif


/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
/* Insert a disk with a new lineage displacement */
#ifdef PHYREX
void MCMC_PHYREX_Insert_Hit(phydbl hr, int n_insert_disks, phydbl cur_rad, phydbl cur_mu, t_tree *tree)
{
  t_dsk  *disk,**new_disk,*young_disk,*old_disk;
  t_ldsk **young_ldsk, **old_ldsk, **new_ldsk;
  phydbl T,t;
  phydbl cur_glnL, new_glnL;
  phydbl u,alpha,ratio;
  int i,j,*dir_old_young,n_valid_disks,err;

  tree->mcmc->run_move[tree->mcmc->num_move_phyrex_indel_hit]++;

  if(n_insert_disks == 0) 
    {
      tree->mmod->rad = cur_rad;
      tree->mmod->mu = cur_mu;
      return; 
    }


  disk      = NULL;
  new_glnL  = tree->mmod->c_lnL;
  cur_glnL  = tree->mmod->c_lnL;
  ratio     = 0.0;
  
  disk = tree->young_disk->prev;
  n_valid_disks = 0;
  do
    {
      if(disk->ldsk != NULL && disk->ldsk->n_next == 1) n_valid_disks++;
      disk = disk->prev;
    }
  while(disk && disk->prev);
  

  new_disk      = (t_dsk **)mCalloc(n_insert_disks,sizeof(t_dsk *));
  new_ldsk      = (t_ldsk **)mCalloc(n_insert_disks,sizeof(t_ldsk *));
  old_ldsk      = (t_ldsk **)mCalloc(n_insert_disks,sizeof(t_ldsk *));
  young_ldsk    = (t_ldsk **)mCalloc(n_insert_disks,sizeof(t_ldsk *));
  dir_old_young = (int *)mCalloc(n_insert_disks,sizeof(int));

  T = PHYREX_Tree_Height(tree);

  for(j=0;j<n_insert_disks;j++)
    {  
      /* Time of insertion of new disk */
      t = Uni()*(tree->young_disk->time-T) + T;
      disk = tree->young_disk;
      while(disk && disk->time > t) disk = disk->prev;
      
      /* Disks located just prior and after inserted disk */
      young_disk = disk->next;
      assert(young_disk->n_ldsk_a);

      old_disk = disk;
      assert(old_disk->n_ldsk_a);
      
      /* Make and initialize new disk */ 
      new_disk[j] = PHYREX_Make_Disk_Event(tree->mmod->n_dim,tree->n_otu);
      PHYREX_Init_Disk_Event(new_disk[j],tree->mmod->n_dim,tree->mmod);
        
      /* Time of the new disk */
      new_disk[j]->time = t;
      
      /* Insert it */
      PHYREX_Insert_Disk(new_disk[j],tree);

      /* Which lineage is going to be hit ? */    
      hr -= log(1./PHYREX_Number_Of_Outgoing_Ldsks(young_disk));
      
      young_ldsk[j] = PHYREX_Random_Select_Outgoing_Ldsk(young_disk);
      old_ldsk[j]   = young_ldsk[j]->prev;
      
      if(old_ldsk[j]->disk->time > young_ldsk[j]->disk->time ||
         young_ldsk[j]->disk->time < new_disk[j]->time ||
         old_ldsk[j]->disk->time > new_disk[j]->time) 
        {
          PhyML_Fprintf(stderr,"\n. young_disk->time: %f",young_disk->time);
          PhyML_Fprintf(stderr,"\n. old_disk->time: %f",old_disk->time);
          PhyML_Fprintf(stderr,"\n. young_ldsk->disk->time: %f",young_ldsk[j]->disk->time);
          PhyML_Fprintf(stderr,"\n. old_ldsk->disk->time: %f",old_ldsk[j]->disk->time);
          PhyML_Fprintf(stderr,"\n. new_disk->disk->time: %f",new_disk[j]->time);
          Generic_Exit(__FILE__,__LINE__,__FUNCTION__);    
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

      /* NOT OPTIMAL */
      /* Update ldsk_a arrays in the time interval affected by the insertion */
      PHYREX_Update_Lindisk_List_Range(new_ldsk[j]->disk,old_ldsk[j]->disk,tree);
      
      /* Sample position of the displaced ldsk */
      for(i=0;i<tree->mmod->n_dim;i++)
        {
          new_disk[j]->centr->lonlat[i] = Uni() * (tree->mmod->lim_up->lonlat[i]-tree->mmod->lim_do->lonlat[i])+tree->mmod->lim_do->lonlat[i];
          hr -= log(1./(tree->mmod->lim_up->lonlat[i]-tree->mmod->lim_do->lonlat[i]));
          
          err = NO;
          new_ldsk[j]->coord->lonlat[i] = Rnorm_Trunc(new_disk[j]->centr->lonlat[i],
                                                      1.0*tree->mmod->rad,
                                                      tree->mmod->lim_do->lonlat[i],
                                                      tree->mmod->lim_up->lonlat[i],&err);
          if(err == YES) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
          
          hr -= Log_Dnorm_Trunc(new_ldsk[j]->coord->lonlat[i],
                                new_disk[j]->centr->lonlat[i],
                                1.0*tree->mmod->rad,
                                tree->mmod->lim_do->lonlat[i],
                                tree->mmod->lim_up->lonlat[i],&err);
        }
    }

  T = PHYREX_Tree_Height(tree);
  T = fabs(T-tree->young_disk->time);
  
  hr -= LnChoose(n_valid_disks+n_insert_disks,n_insert_disks);
  hr -= n_insert_disks * log(1./T);
  hr -= LnFact(n_insert_disks);

  if(tree->eval_glnL == YES) new_glnL = PHYREX_Lk(tree);

  ratio += (new_glnL - cur_glnL);
  ratio += hr;  

  ratio = exp(ratio);
  alpha = MIN(1.,ratio);

  if(tree->mcmc->always_yes == YES && new_glnL > UNLIKELY) alpha = 1.0;

  u = Uni();
  
  /* printf("\n- Insert hit %15f %15f %5d",new_glnL - cur_glnL, alpha, n_insert_disks); */

  assert(isnan(u) == NO && isinf(fabs(u)) == NO);

  if(u > alpha) /* Reject */
    {
      /* printf("+ Reject"); */

      tree->mmod->rad = cur_rad;
      tree->mmod->mu  = cur_mu;
      
      for(j=n_insert_disks-1;j>=0;j--)
        {
          old_ldsk[j]->next[dir_old_young[j]] = young_ldsk[j];
          young_ldsk[j]->prev = old_ldsk[j];

          PHYREX_Remove_Disk(new_disk[j]);      
        }
      

      PHYREX_Update_Lindisk_List(tree);
      tree->mmod->c_lnL = cur_glnL;
        
      for(j=n_insert_disks-1;j>=0;j--)
        {
          Free_Disk(new_disk[j]);
          Free_Ldisk(new_ldsk[j]);
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

  if(n_delete_disks == 0) 
    {
      tree->mmod->rad = cur_rad;
      tree->mmod->mu = cur_mu;
      return; 
    }

  disk            = NULL;
  valid_disks     = NULL;
  new_glnL        = tree->mmod->c_lnL;
  cur_glnL        = tree->mmod->c_lnL;  
  ratio           = 0.0;
  block           = 100;

  if(tree->young_disk->next) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);

  disk = tree->young_disk->prev;
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

  if(!n_valid_disks || (n_valid_disks - n_delete_disks < 0)) 
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

  for(j=0;j<n_delete_disks;j++)
    {
      target_disk[j] = valid_disks[permut[j]];  
      target_ldsk[j] = target_disk[j]->ldsk;

      assert(target_disk[j]->age_fixed == NO);      
      assert(target_ldsk[j] != NULL);
      assert(target_ldsk[j]->n_next == 1);

      old_ldsk[j]   = target_ldsk[j]->prev;
      young_ldsk[j] = target_ldsk[j]->next[0];
            
      dir_old_young[j] = PHYREX_Get_Next_Direction(young_ldsk[j],old_ldsk[j]);
      assert(dir_old_young[j] != -1);


      /* Part of the Hastings ratio corresponding to the probability of selecting */
      /* target_disk->ldsk->next[0] to be hit (reverse move) */ 
      /* hr += log(1./target_disk[j]->next->n_ldsk_a); */
      
      /* Density for position of the displaced ldsk */
      for(i=0;i<tree->mmod->n_dim;i++)
        {
          hr += log(1./(tree->mmod->lim_up->lonlat[i]-tree->mmod->lim_do->lonlat[i]));
          
          hr += Log_Dnorm_Trunc(target_ldsk[j]->coord->lonlat[i],
                                target_disk[j]->centr->lonlat[i],
                                1.0*tree->mmod->rad,
                                tree->mmod->lim_do->lonlat[i],
                                tree->mmod->lim_up->lonlat[i],&err);
        }


      /* New connections between old_ldsk and young_ldsk */
      old_ldsk[j]->next[dir_old_young[j]] = young_ldsk[j];
      young_ldsk[j]->prev                 = old_ldsk[j];
      
      /* Remove target disk */
      PHYREX_Remove_Disk(target_disk[j]);

      /* Update ldsk_a arrays in the time interval affected by the deletion */
      PHYREX_Update_Lindisk_List_Range(young_ldsk[j]->disk,old_ldsk[j]->disk,tree);
    }
  
  T = PHYREX_Tree_Height(tree);
  T = fabs(T-tree->young_disk->time);
  
  hr += LnChoose(n_valid_disks,n_delete_disks);
  hr += n_delete_disks * log(1./T);
  hr += LnFact(n_delete_disks);
  
  Free(valid_disks);

  if(tree->eval_glnL == YES) new_glnL = PHYREX_Lk(tree);

  ratio += (new_glnL - cur_glnL);
  ratio += hr;

  ratio = exp(ratio);
  alpha = MIN(1.,ratio);
  
  /* Always accept move */
  if(tree->mcmc->always_yes == YES && new_glnL > UNLIKELY) alpha = 1.0;

  u = Uni();
  
  /* printf("\n- Delete hit %15f %15f %5d",new_glnL - cur_glnL, alpha, n_delete_disks); */

  assert(isnan(u) == NO && isinf(fabs(u)) == NO);

  if(u > alpha) /* Reject */
    {
      tree->mmod->rad = cur_rad;
      tree->mmod->mu = cur_mu;
      
      for(j=n_delete_disks-1;j>=0;j--) 
        {
          PHYREX_Insert_Disk(target_disk[j],tree);      
          old_ldsk[j]->next[dir_old_young[j]] = target_ldsk[j];
          young_ldsk[j]->prev                 = target_ldsk[j];
        }
      
      PHYREX_Update_Lindisk_List(tree);

      tree->mmod->c_lnL = cur_glnL;
    }
  else
    {
      for(j=0;j<n_delete_disks;j++) Free_Disk(target_disk[j]);
      for(j=0;j<n_delete_disks;j++) Free_Ldisk(target_ldsk[j]);
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
  phydbl cur_rlnL, new_rlnL;
  t_dsk  *disk,*regraft_disk,**valid_disks,*oldest_disk;
  t_ldsk *prune_ldsk,*regraft_ldsk,*prune_ldsk_daughter,*cur_path,*new_path,*ldsk,*ldsk_dum,**valid_ldsks;
  phydbl *prob_disks_prune,*prob_disks_regraft;
  int i,block,n_valid_disks,num_regraft_disk,num_prune_disk,n_valid_ldsks;
  phydbl max_dist, param_exp;
  int cur_path_len;
  int n_iter;
  int cur_pos,new_pos;
  phydbl sum;

  if(tree->mod->s_opt->opt_topo == NO) return;

  n_iter = 1+(int)(tree->n_otu/10);

  while(n_iter--)
    {      
      tree->mcmc->run_move[tree->mcmc->num_move_phyrex_spr]++;
  
      valid_ldsks        = NULL;
      valid_disks        = NULL;
      disk               = NULL;
      prob_disks_regraft = NULL;
      prob_disks_prune   = NULL;
      new_glnL           = tree->mmod->c_lnL;
      cur_glnL           = tree->mmod->c_lnL;
      new_alnL           = tree->c_lnL;
      cur_alnL           = tree->c_lnL;
      new_rlnL           = tree->rates->c_lnL_rates;
      cur_rlnL           = tree->rates->c_lnL_rates;
      hr                 = 0.0;
      ratio              = 0.0;
      block              = 100;
      cur_pos            = -1;
      new_pos            = -1;
      param_exp          = 1.0;
      num_regraft_disk   = -1;
      num_prune_disk     = -1;
      
      if(tree->young_disk->next) assert(FALSE);
            
      /* Get a ldsk from which you can prune a lineage */
      disk = tree->young_disk->prev;
      n_valid_ldsks = 0;
      do
        {
          /* Include only disks with displacement that are coalescent events */
          if(disk->ldsk != NULL && disk->ldsk->n_next > 1)
            {
              /* Root has degree 2: don't pull any lineage */
              if(!(disk->prev == NULL && disk->ldsk->n_next == 2))
                {
                  if(!n_valid_ldsks) valid_ldsks = (t_ldsk **)mCalloc(disk->ldsk->n_next,sizeof(t_ldsk *));
                  else valid_ldsks = (t_ldsk **)mRealloc(valid_ldsks,n_valid_ldsks+disk->ldsk->n_next,sizeof(t_ldsk *));
                  for(i=0;i<disk->ldsk->n_next;++i)
                    {
                      valid_ldsks[n_valid_ldsks] = disk->ldsk->next[i];
                      n_valid_ldsks++;
                    }
                }
            }
          disk = disk->prev;
        }
      while(disk);
      
      if(!n_valid_ldsks) return;
      
      prune_ldsk_daughter = valid_ldsks[Rand_Int(0,n_valid_ldsks-1)];
      prune_ldsk = prune_ldsk_daughter->prev;
      hr -= log(1./(phydbl)(n_valid_ldsks));
      Free(valid_ldsks);

      /* prune_ldsk_daughter is the next coalescent event or tip node */
      while(prune_ldsk_daughter->n_next == 1) prune_ldsk_daughter = prune_ldsk_daughter->next[0];

     
      /* Get a ldsk to reattach the pruned lineage to */
      disk = tree->young_disk->prev;
      n_valid_disks = 0;
      do
        {
          /* Include only disks with displacement that are not younger
             than prune_ldsk_daughter and different from prune_ldsk */
          if((disk->ldsk != NULL) &&
             (disk->ldsk->n_next >= 1) &&
             (disk->time < prune_ldsk_daughter->disk->time) && 
             (PHYREX_Is_On_Path(disk->ldsk,prune_ldsk_daughter,prune_ldsk) == NO))
            {
              if(!n_valid_disks) 
                {
                  valid_disks        = (t_dsk **)mCalloc(block,sizeof(t_dsk *));
                  prob_disks_prune   = (phydbl *)mCalloc(block,sizeof(phydbl));
                }
              else if(!(n_valid_disks%block)) 
                {
                  valid_disks        = (t_dsk **)mRealloc(valid_disks,n_valid_disks+block,sizeof(t_dsk *));
                  prob_disks_prune   = (phydbl *)mRealloc(prob_disks_prune,n_valid_disks+block,sizeof(phydbl));
                }

              valid_disks[n_valid_disks] = disk;

              prob_disks_prune[n_valid_disks] = PHYREX_Dist_Between_Two_Ldsk(disk->ldsk,prune_ldsk,tree);
              
              n_valid_disks++;
            }
          disk = disk->prev;
        }
      while(disk);
      
      assert(n_valid_disks);

      if(n_valid_disks > 1)
        {
          max_dist = -INFINITY;
          for(i=0;i<n_valid_disks;i++) if(prob_disks_prune[i] > max_dist) max_dist = prob_disks_prune[i];
          
          sum = 0.0;
          for(i=0;i<n_valid_disks;i++) prob_disks_prune[i] /= max_dist;
          
          for(i=0;i<n_valid_disks;i++) prob_disks_prune[i] = Dexp(prob_disks_prune[i],param_exp);
          
          sum = 0.0;
          for(i=0;i<n_valid_disks;i++) sum += prob_disks_prune[i];
          for(i=0;i<n_valid_disks;i++) prob_disks_prune[i] /= sum;
          
          num_regraft_disk = Sample_i_With_Proba_pi(prob_disks_prune,n_valid_disks);
          regraft_disk = valid_disks[num_regraft_disk];
          regraft_ldsk = regraft_disk->ldsk;
          
          /* Prob of selecting this node as regraft site */
          hr -= log(prob_disks_prune[num_regraft_disk]);
        }
      
      Free(prob_disks_prune);
      Free(valid_disks);
            
      if(n_valid_disks == 1) return;
      
      oldest_disk = (regraft_ldsk->disk->time < prune_ldsk->disk->time) ? (regraft_ldsk->disk) : (prune_ldsk->disk); 
            
      new_glnL = cur_glnL;
      if(tree->eval_glnL == YES)
        {
          new_glnL -= PHYREX_Lk_Range(prune_ldsk_daughter->disk,oldest_disk,tree);
        }
      
      /* Prune and regraft */
      cur_path_len = PHYREX_Path_Len(prune_ldsk_daughter,prune_ldsk)-2;
      
      hr += PHYREX_Path_Logdensity(prune_ldsk_daughter,prune_ldsk,1.0*tree->mmod->rad,tree);
      
      new_path = PHYREX_Generate_Path(prune_ldsk_daughter,regraft_ldsk,cur_path_len,1.0*tree->mmod->rad,tree);
      cur_path = PHYREX_Remove_Path(prune_ldsk_daughter,prune_ldsk,&cur_pos,tree);
      new_pos = Rand_Int(0,regraft_ldsk->n_next);
      hr -= log(1./(phydbl)(regraft_ldsk->n_next+1));
      hr += log(1./(phydbl)(prune_ldsk->n_next+1));
      PHYREX_Insert_Path(prune_ldsk_daughter,regraft_ldsk,new_path,new_pos,tree);

      hr -= PHYREX_Path_Logdensity(prune_ldsk_daughter,regraft_ldsk,1.0*tree->mmod->rad,tree);
            

      PHYREX_Ldsk_To_Tree(tree);
      Update_Ancestors(tree->n_root,tree->n_root->v[2],tree);
      Update_Ancestors(tree->n_root,tree->n_root->v[1],tree);
      RATES_Update_Cur_Bl(tree);

      assert(isinf(hr) == NO);

      PHYREX_Update_Lindisk_List(tree);
      
      /* Number of ldsk which can be pruned (reverse move) */
      if(regraft_disk->prev == NULL && regraft_ldsk->n_next == 3) // regraft_ldsk is root & went from 2 to 3 outgoing lineages
        n_valid_ldsks+=3;
      if(prune_ldsk->prev == NULL && prune_ldsk->n_next == 2) // prune_ldsk is root & went from 3 to 2 outgoing lineages
        n_valid_ldsks-=3;
      if(regraft_ldsk->n_next == 2) n_valid_ldsks += 2;
      if(prune_ldsk->n_next == 1) n_valid_ldsks -= 2;
        
      assert(n_valid_ldsks > 0);
      hr += log(1./(phydbl)(n_valid_ldsks));




      /* Get a ldsk to reattach the pruned lineage to (reverse move) */
      disk = tree->young_disk->prev;
      n_valid_disks = 0;
      do
        {
          /* Include only disks with displacement that are not younger
             than prune_ldsk_daughter and different from prune_ldsk */
          if((disk->ldsk != NULL) &&
             (disk->ldsk->n_next >= 1) &&
             (disk->time < prune_ldsk_daughter->disk->time) && 
             (PHYREX_Is_On_Path(disk->ldsk,prune_ldsk_daughter,regraft_ldsk) == NO))
            {
              if(!n_valid_disks) 
                {
                  valid_disks        = (t_dsk **)mCalloc(block,sizeof(t_dsk *));
                  prob_disks_regraft = (phydbl *)mCalloc(block,sizeof(phydbl));
                }
              else if(!(n_valid_disks%block)) 
                {
                  valid_disks        = (t_dsk **)mRealloc(valid_disks,n_valid_disks+block,sizeof(t_dsk *));
                  prob_disks_regraft = (phydbl *)mRealloc(prob_disks_regraft,n_valid_disks+block,sizeof(phydbl));
                }

              valid_disks[n_valid_disks] = disk;

              prob_disks_regraft[n_valid_disks] = PHYREX_Dist_Between_Two_Ldsk(disk->ldsk,regraft_ldsk,tree);

              if(disk == prune_ldsk->disk) num_prune_disk = n_valid_disks;
              
              n_valid_disks++;
            }
          disk = disk->prev;
        }
      while(disk);
      
      max_dist = -INFINITY;
      for(i=0;i<n_valid_disks;i++) if(prob_disks_regraft[i] > max_dist) max_dist = prob_disks_regraft[i];
      
      sum = 0.0;
      for(i=0;i<n_valid_disks;i++) prob_disks_regraft[i] /= max_dist;
      
      for(i=0;i<n_valid_disks;i++) prob_disks_regraft[i] = Dexp(prob_disks_regraft[i],param_exp);
      
      sum = 0.0;
      for(i=0;i<n_valid_disks;i++) sum += prob_disks_regraft[i];
      for(i=0;i<n_valid_disks;i++) prob_disks_regraft[i] /= sum;
      
      
      assert(num_prune_disk > -1);
      /* Prob of selecting this node as regraft site (reverse move) */
      hr += log(prob_disks_regraft[num_prune_disk]);
    
      Free(prob_disks_regraft);
      Free(valid_disks);
      
      if(tree->eval_glnL == YES)
        {
          new_glnL += PHYREX_Lk_Range(prune_ldsk_daughter->disk,oldest_disk,tree);
          tree->mmod->c_lnL = new_glnL;
          /* new_glnL = PHYREX_Lk(tree); */
        }
      
      if(tree->eval_alnL == YES) new_alnL = Lk(NULL,tree);
      if(tree->eval_rlnL == YES) new_rlnL = RATES_Lk_Rates(tree);
      
      ratio += (new_alnL - cur_alnL);
      ratio += (new_glnL - cur_glnL);
      ratio += (new_rlnL - cur_rlnL);

      ratio += hr;
            
      ratio = exp(ratio);
      alpha = MIN(1.,ratio);
      
      /* PhyML_Printf("\nYYY %12f %12f %4d %4d %12f %12f", */
      /*              prune_ldsk_daughter->disk->time - prune_ldsk->disk->time, */
      /*              prune_ldsk_daughter->disk->time - regraft_ldsk->disk->time, */
      /*              cur_path_len, */
      /*              new_path_len, */
      /*              PHYREX_Dist_Between_Two_Ldsk(regraft_ldsk,prune_ldsk,tree), */
      /*              alpha); */


      /* Always accept move */
      if(tree->mcmc->always_yes == YES && new_glnL > UNLIKELY) alpha = 1.0;
      
      /* printf("\n. %3d %3d alpha:%12f %3d",prune_ldsk->n_next+1,regraft_ldsk->n_next-1,alpha,new_n_coal); */
      /* PhyML_Printf("\n. new_glnL: %f cur_glnL: %f hr: %f alpha: %f",new_glnL,cur_glnL,hr,alpha); */
      
      u = Uni();
      
      assert(isnan(u) == NO && isinf(fabs(u)) == NO);

      if(u > alpha) /* Reject */
        {
          new_path = PHYREX_Remove_Path(prune_ldsk_daughter,regraft_ldsk,&new_pos,tree);                    
          PHYREX_Insert_Path(prune_ldsk_daughter,prune_ldsk,cur_path,cur_pos,tree);
                    
          PHYREX_Ldsk_To_Tree(tree);
          Update_Ancestors(tree->n_root,tree->n_root->v[2],tree);
          Update_Ancestors(tree->n_root,tree->n_root->v[1],tree);
          RATES_Fill_Lca_Table(tree);
          RATES_Update_Cur_Bl(tree);
          

          PHYREX_Update_Lindisk_List(tree);
          tree->c_lnL              = cur_alnL;
          tree->mmod->c_lnL        = cur_glnL;          
          tree->rates->c_lnL_rates = cur_rlnL;

          ldsk = new_path;
          while(ldsk)
            {
              Free_Disk(ldsk->disk);
              ldsk_dum = ldsk;
              ldsk = ldsk->prev;
              Free_Ldisk(ldsk_dum);
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
  
/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

#ifdef PHYREX
void MCMC_PHYREX_Prune_Regraft_Local(t_tree *tree)
{
  phydbl u,alpha,ratio,hr;
  phydbl cur_glnL, new_glnL;
  phydbl cur_alnL, new_alnL;
  phydbl cur_rlnL, new_rlnL;
  t_dsk  *disk,*regraft_disk,*oldest_disk,**valid_disks;
  t_ldsk *prune_ldsk,*prune_ldsk_daughter,*regraft_ldsk,*cur_path,**valid_ldsks;
  int i,block,n_valid_disks,n_valid_ldsks;
  int n_iter;
  int cur_pos,new_pos;

  if(tree->mod->s_opt->opt_topo == NO) return;
  
  n_iter = 1+(int)(tree->n_otu/10);

  while(n_iter--)
    {      
      tree->mcmc->run_move[tree->mcmc->num_move_phyrex_spr_local]++;
  
      valid_ldsks    = NULL;
      valid_disks    = NULL;
      disk           = NULL;
      new_glnL       = tree->mmod->c_lnL;
      cur_glnL       = tree->mmod->c_lnL;
      new_alnL       = tree->c_lnL;
      cur_alnL       = tree->c_lnL;
      new_rlnL       = tree->rates->c_lnL_rates;
      cur_rlnL       = tree->rates->c_lnL_rates;
      hr             = 0.0;
      ratio          = 0.0;
      block          = 100;
      cur_pos        = -1;
      new_pos        = -1;

      if(tree->young_disk->next) assert(FALSE);
      
      /* Get a ldsk from which you can prune a lineage */
      disk = tree->young_disk->prev;
      n_valid_ldsks = 0;
      do
        {
          /* Include only disks with displacement that are coalescent events */
          if(disk->ldsk != NULL && disk->ldsk->n_next > 1)
            {
              /* Root has degree 2: don't pull any lineage */
              if(!(disk->prev == NULL && disk->ldsk->n_next == 2))
                {
                  if(!n_valid_ldsks) valid_ldsks = (t_ldsk **)mCalloc(disk->ldsk->n_next,sizeof(t_ldsk *));
                  else valid_ldsks = (t_ldsk **)mRealloc(valid_ldsks,n_valid_ldsks+disk->ldsk->n_next,sizeof(t_ldsk *));
                  for(i=0;i<disk->ldsk->n_next;++i)
                    {
                      valid_ldsks[n_valid_ldsks] = disk->ldsk->next[i];
                      n_valid_ldsks++;
                    }
                }
            }
          disk = disk->prev;
        }
      while(disk);
      
      if(!n_valid_ldsks) return;
      
      prune_ldsk_daughter = valid_ldsks[Rand_Int(0,n_valid_ldsks-1)];
      prune_ldsk = prune_ldsk_daughter->prev;
      hr -= log(1./(phydbl)(n_valid_ldsks));
      Free(valid_ldsks);


      /* Get a ldsk to reattach the pruned lineage to */
      n_valid_disks = 0;
      if(prune_ldsk->prev != NULL)
        {
          valid_disks = (t_dsk **)mCalloc(block,sizeof(t_dsk *));
          valid_disks[0] = prune_ldsk->prev->disk;
          n_valid_disks++;
        }

      for(i=0;i<prune_ldsk->n_next;++i)
        {
          if(prune_ldsk->next[i] != prune_ldsk_daughter &&
             prune_ldsk->next[i]->n_next >= 1 &&
             prune_ldsk->next[i]->disk->time < prune_ldsk_daughter->disk->time)
            {
              if(n_valid_disks == 0)
                {
                  valid_disks = (t_dsk **)mCalloc(block,sizeof(t_dsk *));
                }
              else
                {
                  if(!(n_valid_disks%block)) 
                    {
                      valid_disks = (t_dsk **)mRealloc(valid_disks,n_valid_disks+block,sizeof(t_dsk *));
                    }
                }
              
              valid_disks[n_valid_disks] = prune_ldsk->next[i]->disk;

              n_valid_disks++;
            }
        }
      
      if(!n_valid_disks) return;

      regraft_disk = valid_disks[Rand_Int(0,n_valid_disks-1)];            
      regraft_ldsk = regraft_disk->ldsk;
      hr -= log(1./(phydbl)(n_valid_disks));
      Free(valid_disks);
      
      oldest_disk = (regraft_disk->time < prune_ldsk->disk->time) ? (regraft_disk) : (prune_ldsk->disk);
      
      new_glnL = cur_glnL;
      if(tree->eval_glnL == YES)
        {
          new_glnL -= PHYREX_Lk_Range(prune_ldsk_daughter->disk,oldest_disk,tree);
        }
      
      /* Prune and regraft */
      cur_path = PHYREX_Remove_Path(prune_ldsk_daughter,prune_ldsk,&cur_pos,tree);
      new_pos = Rand_Int(0,regraft_ldsk->n_next);
      hr -= log(1./(phydbl)(regraft_ldsk->n_next+1));
      hr += log(1./(phydbl)(prune_ldsk->n_next+1));
      PHYREX_Insert_Path(prune_ldsk_daughter,regraft_ldsk,cur_path,new_pos,tree);

      PHYREX_Ldsk_To_Tree(tree);
      Update_Ancestors(tree->n_root,tree->n_root->v[2],tree);
      Update_Ancestors(tree->n_root,tree->n_root->v[1],tree);
      RATES_Update_Cur_Bl(tree);

      assert(isinf(hr) == NO);

      PHYREX_Update_Lindisk_List_Range(prune_ldsk_daughter->disk,oldest_disk,tree);


      /* Number of ldsk which one can be pruned (reverse move) */
      if(regraft_disk->prev == NULL && regraft_ldsk->n_next == 3) // regraft_ldsk is root & went from 2 to 3 outgoing lineages
        n_valid_ldsks+=3;
      if(prune_ldsk->prev == NULL && prune_ldsk->n_next == 2) // prune_ldsk is root & went from 3 to 2 outgoing lineages
        n_valid_ldsks-=3;
      if(regraft_ldsk->n_next == 2) n_valid_ldsks += 2;
      if(prune_ldsk->n_next == 1) n_valid_ldsks -= 2;
        
            
      assert(n_valid_ldsks > 0);
      hr += log(1./(phydbl)(n_valid_ldsks));
      
      
      /* Get a ldsk to reattach the pruned lineage to (reverse move) */
      n_valid_disks = 0;
      if(regraft_ldsk->prev != NULL) n_valid_disks++;

      for(i=0;i<regraft_ldsk->n_next;++i)
        {
          if(regraft_ldsk->next[i] != prune_ldsk_daughter &&
             regraft_ldsk->next[i]->n_next >= 1 &&
             regraft_ldsk->next[i]->disk->time < prune_ldsk_daughter->disk->time)
            {
              n_valid_disks++;
            }
        }

      assert(n_valid_disks > 0);
      hr += log(1./(phydbl)(n_valid_disks));

      if(tree->eval_glnL == YES)
        {
          new_glnL += PHYREX_Lk_Range(prune_ldsk_daughter->disk,oldest_disk,tree);
          tree->mmod->c_lnL = new_glnL;
          /* new_glnL = PHYREX_Lk(tree); */
        }
      
      if(tree->eval_alnL == YES) new_alnL = Lk(NULL,tree);
      if(tree->eval_rlnL == YES) new_rlnL = RATES_Lk_Rates(tree);
      
      ratio += (new_alnL - cur_alnL);
      ratio += (new_glnL - cur_glnL);
      ratio += (new_rlnL - cur_rlnL);

      ratio += hr;
            
      ratio = exp(ratio);
      alpha = MIN(1.,ratio);
      
      /* PhyML_Printf("\nZZZ %12f %12f %12f %12f", */
      /*              prune_ldsk_daughter->disk->time - prune_ldsk->disk->time, */
      /*              prune_ldsk_daughter->disk->time - regraft_ldsk->disk->time, */
      /*              PHYREX_Dist_Between_Two_Ldsk(regraft_ldsk,prune_ldsk,tree), */
      /*              alpha); */


      /* Always accept move */
      if(tree->mcmc->always_yes == YES && new_glnL > UNLIKELY) alpha = 1.0;
      
      /* printf("\n. %3d %3d alpha:%12f %3d",prune_ldsk->n_next+1,regraft_ldsk->n_next-1,alpha,new_n_coal); */
      /* PhyML_Printf("\n. new_glnL: %f cur_glnL: %f hr: %f alpha: %f",new_glnL,cur_glnL,hr,alpha); */
      
      u = Uni();
      
      assert(isnan(u) == NO && isinf(fabs(u)) == NO);

      if(u > alpha) /* Reject */
        {
          cur_path = PHYREX_Remove_Path(prune_ldsk_daughter,regraft_ldsk,&new_pos,tree);                    
          PHYREX_Insert_Path(prune_ldsk_daughter,prune_ldsk,cur_path,cur_pos,tree);
                    
          PHYREX_Ldsk_To_Tree(tree);
          Update_Ancestors(tree->n_root,tree->n_root->v[2],tree);
          Update_Ancestors(tree->n_root,tree->n_root->v[1],tree);
          RATES_Fill_Lca_Table(tree);
          RATES_Update_Cur_Bl(tree);
          
          PHYREX_Update_Lindisk_List_Range(prune_ldsk_daughter->disk,oldest_disk,tree);
           
          tree->c_lnL              = cur_alnL;
          tree->mmod->c_lnL        = cur_glnL;          
          tree->rates->c_lnL_rates = cur_rlnL;
        }
      else
        {
          tree->mcmc->acc_move[tree->mcmc->num_move_phyrex_spr_local]++;
        }
    }
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
  t_ldsk *old_ldsk,*young_ldsk,*cur_path,*new_path,*ldsk,*ldsk_dum;
  int j,block,n_valid_disks;
  int n_iter,cur_path_len,new_path_len;
  int pos,*permut,dir_next;
  phydbl area;
  
  n_iter      = 0;
  valid_disks = NULL;
  disk        = NULL;
  new_glnL    = tree->mmod->c_lnL;
  cur_glnL    = tree->mmod->c_lnL;
  hr          = 0.0;
  ratio       = 0.0;
  block       = 100;
  pos         = -1;
  dir_next    = -1;
      
  area = 1.0;
  for(j=0;j<tree->mmod->n_dim;++j) area *= (tree->mmod->lim_up->lonlat[j] - tree->mmod->lim_do->lonlat[j]);

  disk = tree->young_disk->prev;
  n_valid_disks = 0;
  do
    {
      if(disk->ldsk && disk->ldsk->n_next >= 2)
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
  n_iter = 1+(int)(0.1*n_valid_disks);
  
  for(j=0;j<n_iter;j++)
    {
      tree->mcmc->run_move[tree->mcmc->num_move_phyrex_traj]++;

      new_glnL    = tree->mmod->c_lnL;
      cur_glnL    = tree->mmod->c_lnL;
      hr          = 0.0;
      ratio       = 0.0;
      dir_next    = -1;
      
      old_ldsk = valid_disks[permut[Rand_Int(0,n_valid_disks-1)]]->ldsk;
      assert(old_ldsk != NULL);

      dir_next = Rand_Int(0,old_ldsk->n_next-1);
      young_ldsk = old_ldsk->next[dir_next];
      while(young_ldsk->n_next == 1) young_ldsk = young_ldsk->next[0];

      
      /* dt = fabs(young_ldsk->disk->time - old_ldsk->disk->time); */
      
      assert(young_ldsk != NULL); 

      new_glnL = cur_glnL;
      if(tree->eval_glnL == YES)
        {
          new_glnL -= PHYREX_Lk_Range(young_ldsk->disk,old_ldsk->disk,tree);
        }
           
      /* phydbl dist = Euclidean_Dist(young_ldsk->coord,old_ldsk->coord); */
      cur_path_len = PHYREX_Path_Len(young_ldsk,old_ldsk)-2;
      
      /* new_path_len = Rpois(dt * 2. * PI * tree->mmod->lbda * tree->mmod->mu * pow(tree->mmod->rad,2) / area); */
      /* new_path_len = Rpois(cur_path_len+0.01); // To be avoided as create non-irreducible markov chain when new_path_len = 0 (since Pr(cur_path_len|new_path_len)=0) */
      /* new_path_len = Rpois((int)(pow(dist,2)/(4.*pow(tree->mmod->rad,2)))+0.01); */
      
      /* hr -= Dpois(new_path_len,dt * 2.* PI * tree->mmod->lbda * tree->mmod->mu * pow(tree->mmod->rad,2) / area,YES); */
      /* hr += Dpois(cur_path_len,dt * 2.* PI * tree->mmod->lbda * tree->mmod->mu * pow(tree->mmod->rad,2) / area,YES); */
      /* hr -= Dpois(new_path_len,cur_path_len+0.01,YES); */
      /* hr += Dpois(cur_path_len,new_path_len+0.01,YES); */
      /* hr -= Dpois(new_path_len,(int)(pow(dist,2)/(4.*pow(tree->mmod->rad,2)))+0.01,YES); */
      /* hr += Dpois(cur_path_len,(int)(pow(dist,2)/(4.*pow(tree->mmod->rad,2)))+0.01,YES); */
      
      /* hr -= (new_path_len) * log(1./dt); */
      /* hr += (cur_path_len) * log(1./dt); */
      
      /* hr -= LnFact(new_path_len); */
      /* hr += LnFact(cur_path_len); */

      new_path_len = cur_path_len;
            
      hr += PHYREX_Path_Logdensity(young_ldsk,old_ldsk,1.0*tree->mmod->rad,tree);
      
      new_path = PHYREX_Generate_Path(young_ldsk,old_ldsk,new_path_len,1.0*tree->mmod->rad,tree);
      cur_path = PHYREX_Remove_Path(young_ldsk,old_ldsk,&pos,tree);
      PHYREX_Insert_Path(young_ldsk,old_ldsk,new_path,pos,tree);
      
      hr -= PHYREX_Path_Logdensity(young_ldsk,old_ldsk,1.0*tree->mmod->rad,tree);
      
      if(tree->eval_glnL == YES)
        {
          new_glnL += PHYREX_Lk_Range(young_ldsk->disk,old_ldsk->disk,tree);
          tree->mmod->c_lnL = new_glnL;
          /* new_glnL = PHYREX_Lk(tree); */
        }

      ratio += (new_glnL - cur_glnL);
      ratio += hr;
            
      ratio = exp(ratio);
      alpha = MIN(1.,ratio);
      
      /* PhyML_Printf("\nXXX %12f %4d %4d %12f %12f", */
      /*              young_ldsk->disk->time - old_ldsk->disk->time, */
      /*              cur_path_len,new_path_len,alpha,Euclidean_Dist(young_ldsk->coord,old_ldsk->coord)); */

      if(tree->mcmc->always_yes == YES && new_glnL > UNLIKELY) alpha = 1.0;


      u = Uni();
      
      if(u > alpha) /* Reject */
        {
          new_path = PHYREX_Remove_Path(young_ldsk,old_ldsk,&pos,tree);                    
          PHYREX_Insert_Path(young_ldsk,old_ldsk,cur_path,pos,tree);
          
          PHYREX_Update_Lindisk_List(tree);
          tree->mmod->c_lnL = cur_glnL;
          
          ldsk = new_path;
          while(ldsk) 
            {
              Free_Disk(ldsk->disk);
              ldsk_dum = ldsk;
              ldsk = ldsk->prev;
              Free_Ldisk(ldsk_dum);
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

  disk = tree->young_disk->prev;
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
  
  /* n_move_disks = Rand_Int(1,1+(int)(n_all_disks/10)); */
  /* n_move_disks = MIN(10,(int)(1.+0.1*n_all_disks)); */
  n_move_disks = (int)(1+n_all_disks/20);

  permut = Permutate(n_all_disks);

  for(i=0;i<n_move_disks;i++)
    {
      target_disk[i] = all_disks[permut[i]];
  
      assert(target_disk[i]);

      PHYREX_Store_Geo_Coord(target_disk[i]->centr);
      
      if(target_disk[i]->ldsk != NULL)
        {
          for(j=0;j<tree->mmod->n_dim;j++)
            {
              err = NO;
              target_disk[i]->centr->lonlat[j] = 
                Rnorm_Trunc(target_disk[i]->centr->lonlat[j],
                            1.*tree->mmod->rad,
                            tree->mmod->lim_do->lonlat[j],
                            tree->mmod->lim_up->lonlat[j],&err);

              if(err == YES) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);

              hr -= Log_Dnorm_Trunc(target_disk[i]->centr->lonlat[j],
                                    target_disk[i]->centr->cpy->lonlat[j],
                                    1.*tree->mmod->rad,                                                         
                                    tree->mmod->lim_do->lonlat[j],
                                    tree->mmod->lim_up->lonlat[j],&err);
              
              hr += Log_Dnorm_Trunc(target_disk[i]->centr->cpy->lonlat[j],
                                    target_disk[i]->centr->lonlat[j],
                                    1.*tree->mmod->rad,
                                    tree->mmod->lim_do->lonlat[j],
                                    tree->mmod->lim_up->lonlat[j],&err);
            }
        }
      else
        {
          for(j=0;j<tree->mmod->n_dim;j++)
            {
              err = NO;
              target_disk[i]->centr->lonlat[j] = 
                Rnorm_Trunc(target_disk[i]->centr->lonlat[j],
                            1.*tree->mmod->rad,
                            tree->mmod->lim_do->lonlat[j],
                            tree->mmod->lim_up->lonlat[j],&err);

              if(err == YES) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);

              hr -= Log_Dnorm_Trunc(target_disk[i]->centr->lonlat[j],
                                    target_disk[i]->centr->cpy->lonlat[j],
                                    1.*tree->mmod->rad,                                                         
                                    tree->mmod->lim_do->lonlat[j],
                                    tree->mmod->lim_up->lonlat[j],&err);
              
              hr += Log_Dnorm_Trunc(target_disk[i]->centr->cpy->lonlat[j],
                                    target_disk[i]->centr->lonlat[j],
                                    1.*tree->mmod->rad,
                                    tree->mmod->lim_do->lonlat[j],
                                    tree->mmod->lim_up->lonlat[j],&err);
            }
        }
    }
  
  Free(permut);

  if(tree->eval_glnL == YES)
    {
      new_glnL = PHYREX_Lk(tree);
    }
  
  ratio += (new_glnL - cur_glnL);
  ratio += hr;
  
  ratio = exp(ratio);
  alpha = MIN(1.,ratio);
  
  /* Always accept move */
  if(tree->mcmc->always_yes == YES && new_glnL > UNLIKELY) alpha = 1.0;

  u = Uni();
  
  /* printf("\n- Delete new_glnL: %f [%f] hr: %f u:%f alpha: %f",new_glnL,cur_glnL,hr,u,alpha); */

  assert(isnan(u) == NO && isinf(fabs(u)) == NO);

  if(u > alpha) /* Reject */
    {
      /* printf("- Reject"); */
      
      for(i=0;i<n_move_disks;i++) PHYREX_Restore_Geo_Coord(target_disk[i]->centr);
  
      tree->mmod->c_lnL = cur_glnL;
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
  phydbl cur_glnL,new_glnL,hr,c;
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
  c             = -1.;
  
  tree->mcmc->run_move[tree->mcmc->num_move_phyrex_ldsk_multi]++;

  if(tree->young_disk->next) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  disk = tree->young_disk->prev;
  n_all_disks = 0;
  do
    {
      if(disk->ldsk != NULL && (disk->ldsk->nd != NULL && disk->ldsk->nd->tax == NO))
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
  
  n_move_ldsk = (int)(1+n_all_disks/20);
  
  target_disk = (t_dsk **)mCalloc(n_all_disks,sizeof(t_dsk *));

  permut = Permutate(n_all_disks);

  for(i=0;i<n_move_ldsk;i++)
    {
      target_disk[i] = all_disks[permut[i]];
      
      PHYREX_Store_Geo_Coord(target_disk[i]->ldsk->coord);

      for(j=0;j<tree->mmod->n_dim;j++)
        {
          /* c: center; o: pos of direct ldsk ancestor */
          /* c = target_disk[i]->centr->lonlat[j]; */
          c = target_disk[i]->ldsk->coord->lonlat[j];

          err = NO;
          target_disk[i]->ldsk->coord->lonlat[j] =
            Rnorm_Trunc(c,
                        1.*tree->mmod->rad,
                        tree->mmod->lim_do->lonlat[j],
                        tree->mmod->lim_up->lonlat[j],&err);
          
          if(err == YES) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
          
          hr -= Log_Dnorm_Trunc(target_disk[i]->ldsk->coord->lonlat[j],
                                c,
                                1.*tree->mmod->rad,
                                tree->mmod->lim_do->lonlat[j],
                                tree->mmod->lim_up->lonlat[j],&err);
          
          /* hr += Log_Dnorm_Trunc(target_disk[i]->ldsk->coord->cpy->lonlat[j], */
          /*                       c, */
          /*                       1.*tree->mmod->rad, */
          /*                       tree->mmod->lim_do->lonlat[j], */
          /*                       tree->mmod->lim_up->lonlat[j],&err); */
          hr += Log_Dnorm_Trunc(c,
                                target_disk[i]->ldsk->coord->lonlat[j],
                                1.*tree->mmod->rad,
                                tree->mmod->lim_do->lonlat[j],
                                tree->mmod->lim_up->lonlat[j],&err);
        }
    }

  Free(permut);

  if(tree->eval_glnL == YES) new_glnL = PHYREX_Lk(tree);

  ratio += (new_glnL - cur_glnL);
  ratio += hr;
  
  ratio = exp(ratio);
  alpha = MIN(1.,ratio);
  
  /* Always accept move */
  if(tree->mcmc->always_yes == YES && new_glnL > UNLIKELY) alpha = 1.0;

  u = Uni();
  
  /* printf("\n- Move_ldsk %15f",new_glnL-cur_glnL); */

  assert(isnan(u) == NO && isinf(fabs(u)) == NO);

  if(u > alpha) /* Reject */
    {
      /* printf("- Reject"); */
      
      for(i=0;i<n_move_ldsk;i++) PHYREX_Restore_Geo_Coord(target_disk[i]->ldsk->coord);

      tree->mmod->c_lnL = cur_glnL;
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
  phydbl cur_rad,new_rad;
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
  cur_rad       = tree->mmod->rad;
  new_rad       = tree->mmod->rad;
  
  tree->mcmc->run_move[tree->mcmc->num_move_phyrex_ldsk_and_disk]++;

  if(tree->young_disk->next) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  disk = tree->young_disk->prev;
  n_all_disks = 0;
  do
    {
      if(disk->ldsk != NULL && disk->ldsk->n_next >= 1 && (disk->ldsk->nd != NULL && disk->ldsk->nd->tax == NO))
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
  
  n_move_ldsk = (int)(1+n_all_disks/20);
  /* n_move_ldsk = n_all_disks; */
  
  target_disk = (t_dsk **)mCalloc(n_all_disks,sizeof(t_dsk *));

  permut = Permutate(n_all_disks);


  new_rad = cur_rad * exp(0.5*(Uni()-.5));
  hr += log(new_rad/cur_rad);

  
  for(i=0;i<n_move_ldsk;i++)
    {
      target_disk[i] = all_disks[permut[i]];
      
      PHYREX_Store_Geo_Coord(target_disk[i]->ldsk->coord);
      PHYREX_Store_Geo_Coord(target_disk[i]->centr);
            
      for(j=0;j<tree->mmod->n_dim;j++)
        {
          err = NO;
          
          target_disk[i]->centr->lonlat[j] =
            Rnorm_Trunc(target_disk[i]->ldsk->coord->lonlat[j],
                        1.*new_rad,
                        tree->mmod->lim_do->lonlat[j],
                        tree->mmod->lim_up->lonlat[j],&err);

          target_disk[i]->ldsk->coord->lonlat[j] =
            Rnorm_Trunc(target_disk[i]->centr->lonlat[j],
                        1.*new_rad,
                        tree->mmod->lim_do->lonlat[j],
                        tree->mmod->lim_up->lonlat[j],&err);
          

          hr -= Log_Dnorm_Trunc(target_disk[i]->centr->lonlat[j],
                                target_disk[i]->ldsk->coord->cpy->lonlat[j],
                                1.*new_rad,
                                tree->mmod->lim_do->lonlat[j],
                                tree->mmod->lim_up->lonlat[j],&err);

          hr += Log_Dnorm_Trunc(target_disk[i]->centr->cpy->lonlat[j],
                                target_disk[i]->ldsk->coord->lonlat[j],
                                1.*new_rad,
                                tree->mmod->lim_do->lonlat[j],
                                tree->mmod->lim_up->lonlat[j],&err);

          hr -= Log_Dnorm_Trunc(target_disk[i]->ldsk->coord->lonlat[j],
                                target_disk[i]->centr->lonlat[j],
                                1.*new_rad,
                                tree->mmod->lim_do->lonlat[j],
                                tree->mmod->lim_up->lonlat[j],&err);
          
          hr += Log_Dnorm_Trunc(target_disk[i]->ldsk->coord->cpy->lonlat[j],
                                target_disk[i]->centr->cpy->lonlat[j],
                                1.*new_rad,
                                tree->mmod->lim_do->lonlat[j],
                                tree->mmod->lim_up->lonlat[j],&err);

        }
    }
      

  Free(permut);

  if(tree->eval_glnL == YES) new_glnL = PHYREX_Lk(tree);

  ratio += (new_glnL - cur_glnL);
  ratio += hr;
  
  ratio = exp(ratio);
  alpha = MIN(1.,ratio);
  
  /* Always accept move */
  if(tree->mcmc->always_yes == YES && new_glnL > UNLIKELY) alpha = 1.0;

  u = Uni();
  
  /* printf("\n- Move_ldsk %15f",new_glnL-cur_glnL); */

  assert(isnan(u) == NO && isinf(fabs(u)) == NO);

  if(u > alpha) /* Reject */
    {
      /* printf("- Reject"); */
      for(i=0;i<n_move_ldsk;i++) PHYREX_Restore_Geo_Coord(target_disk[i]->ldsk->coord);
      for(i=0;i<n_move_ldsk;i++) PHYREX_Restore_Geo_Coord(target_disk[i]->centr);

      tree->mmod->rad = cur_rad;
      tree->mmod->c_lnL = cur_glnL;
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
  phydbl u,alpha,ratio,hr,c;
  phydbl cur_glnL, new_glnL;
  t_dsk  *disk,**all_disks;
  int i,j,err,n_all_disks,block,n_move_ldsk,*permut;

  block       = 100;
  all_disks   = NULL;
  n_all_disks = 0;
  c           = -1.;
  
  disk = tree->young_disk->prev;
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
  /* n_move_ldsk = MIN(10,(int)(1.+0.1*n_all_disks)); */
  n_move_ldsk = 1+(int)(n_all_disks/10);
  
  permut = Permutate(n_all_disks);
  cur_glnL = tree->mmod->c_lnL;
  new_glnL = tree->mmod->c_lnL;
  for(i=0;i<n_move_ldsk;i++)
    {
      tree->mcmc->run_move[tree->mcmc->num_move_phyrex_ldsk_given_disk]++;
      
      disk = all_disks[permut[i]];
            
      hr       = 0.0;
      ratio    = 0.0;
      new_glnL = cur_glnL;   
      
      PHYREX_Store_Geo_Coord(disk->ldsk->coord);
      
      if(tree->eval_glnL == YES)
        {
          if(disk->ldsk->prev != NULL) new_glnL -= PHYREX_Lk_Range(disk,disk->ldsk->prev->disk,tree);
          else new_glnL -= PHYREX_Lk_Core(disk,tree);
        }
      
      for(j=0;j<tree->mmod->n_dim;j++)
        {
          c = disk->centr->lonlat[j];
          
          err = NO;
          disk->ldsk->coord->lonlat[j] =
            Rnorm_Trunc(c,
                        1.*tree->mmod->rad,
                        tree->mmod->lim_do->lonlat[j],
                        tree->mmod->lim_up->lonlat[j],&err);
          
          if(err == YES) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);

          hr -= Log_Dnorm_Trunc(disk->ldsk->coord->lonlat[j],
                                c,
                                1.*tree->mmod->rad,
                                tree->mmod->lim_do->lonlat[j],
                                tree->mmod->lim_up->lonlat[j],&err);
          
          hr += Log_Dnorm_Trunc(disk->ldsk->coord->cpy->lonlat[j],
                                c,
                                1.*tree->mmod->rad,
                                tree->mmod->lim_do->lonlat[j],
                                tree->mmod->lim_up->lonlat[j],&err);

        }
      
     if(tree->eval_glnL == YES)
       {
          if(disk->ldsk->prev != NULL) new_glnL += PHYREX_Lk_Range(disk,disk->ldsk->prev->disk,tree);
          else new_glnL += PHYREX_Lk_Core(disk,tree);
          tree->mmod->c_lnL = new_glnL;
         /* new_glnL = PHYREX_Lk(tree); */
       }
     
      ratio += (new_glnL - cur_glnL);
      ratio += hr;
      
      ratio = exp(ratio);
      alpha = MIN(1.,ratio);
      
      u = Uni();
      
      assert(isnan(u) == NO && isinf(fabs(u)) == NO);

      if(u > alpha) /* Reject */
        {
          PHYREX_Restore_Geo_Coord(disk->ldsk->coord);          
          tree->mmod->c_lnL = cur_glnL;            
        }
      else
        {
          cur_glnL = new_glnL;
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


  disk = tree->young_disk->prev;
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
  
  /* n_move_ldsk = Rand_Int(1,1+(int)(n_all_disks/10)); */
  /* n_move_ldsk = MIN(10,(int)(1+0.1*n_all_disks)); */
  n_move_ldsk = (int)(1+n_all_disks/10);
  
  permut = Permutate(n_all_disks);


  for(i=0;i<n_move_ldsk;i++)
    {
      tree->mcmc->run_move[tree->mcmc->num_move_phyrex_disk_given_ldsk]++;
      
      disk = all_disks[permut[i]];
      
      hr       = 0.0;
      ratio    = 0.0;
      cur_glnL = tree->mmod->c_lnL;      
      new_glnL = tree->mmod->c_lnL;
      
      if(tree->eval_glnL == YES) new_glnL -= PHYREX_Lk_Core(disk,tree);
 
      PHYREX_Store_Geo_Coord(disk->centr);

      for(j=0;j<tree->mmod->n_dim;j++)
        {
          /* if(disk->ldsk == NULL) */
          disk->centr->lonlat[j] = Uni()*(tree->mmod->lim_up->lonlat[j]-tree->mmod->lim_do->lonlat[j])+tree->mmod->lim_do->lonlat[j];
            
          /* else */
          /*   { */
          /*     /\* c: center; o: pos of direct ldsk ancestor *\/ */
          /*     c = disk->ldsk->coord->lonlat[j]; */
              
          /*     err = NO; */
          /*     disk->centr->lonlat[j] = */
          /*       Rnorm_Trunc(c, */
          /*                   2.*tree->mmod->rad, */
          /*                   0.0, */
          /*                   tree->mmod->lim->lonlat[j],&err); */
              
          /*     if(err == YES) Generic_Exit(__FILE__,__LINE__,__FUNCTION__); */
              
          /*     hr -= Log_Dnorm_Trunc(disk->centr->lonlat[j], */
          /*                           c, */
          /*                           2.*tree->mmod->rad, */
          /*                           0.0, */
          /*                           tree->mmod->lim->lonlat[j],&err); */
              
          /*     hr += Log_Dnorm_Trunc(disk->centr->cpy->lonlat[j], */
          /*                           c, */
          /*                           2.*tree->mmod->rad, */
          /*                           0.0, */
          /*                           tree->mmod->lim->lonlat[j],&err); */
          /*   } */
        }
      
      if(tree->eval_glnL == YES)
        {
          new_glnL += PHYREX_Lk_Core(disk,tree);
          tree->mmod->c_lnL = new_glnL;
        }
      
      ratio += (new_glnL - cur_glnL);
      ratio += hr;
      
      ratio = exp(ratio);
      alpha = MIN(1.,ratio);
      
      u = Uni();
      
      assert(isnan(u) == NO && isinf(fabs(u)) == NO);

      if(u > alpha) /* Reject */
        {
          PHYREX_Restore_Geo_Coord(disk->centr);          
          tree->mmod->c_lnL = cur_glnL;
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
  t_dsk *disk,*new_disk,*target_disk,*young_disk,*old_disk,**valid_disks;
  t_ldsk *young_ldsk, *old_ldsk, *new_ldsk;
  int i,j,n_trials,dir_old_young,err,block;
  phydbl ratio, u, alpha, hr, type;
  phydbl cur_glnL, new_glnL;
  phydbl T,t,pindel;
  int n_valid_disks;
  
  cur_glnL     = tree->mmod->c_lnL;
  new_glnL     = tree->mmod->c_lnL;
  hr           = 0.0;
  ratio        = 0.0;
  type         = -1.0;
  n_trials     = 1 + (int)(0.1*PHYREX_Total_Number_Of_Intervals(tree));
  T            = PHYREX_Tree_Height(tree);
  pindel       = 0.5;
  block        = 50;
  
  for(i=0;i<n_trials;i++)
    {
      tree->mcmc->run_move[tree->mcmc->num_move_phyrex_indel_hit_serial]++;
      
      cur_glnL = tree->mmod->c_lnL;
      new_glnL = tree->mmod->c_lnL;
      hr       = 0.0;
      ratio    = 0.0;

      type = Uni();
      
      if(type < pindel) /* Insert */
        {          
          t = Uni()*(tree->young_disk->time-T) + T;
          disk = tree->young_disk->prev;
          while(disk && disk->time > t) disk = disk->prev;
          
          assert(disk->next);

          hr -= log(1./(tree->young_disk->time-T));
          
          young_disk = disk->next;
          assert(young_disk->n_ldsk_a);
          
          old_disk = disk;
          assert(old_disk->n_ldsk_a);
                      
          new_disk = PHYREX_Make_Disk_Event(tree->mmod->n_dim,tree->n_otu);
          PHYREX_Init_Disk_Event(new_disk,tree->mmod->n_dim,tree->mmod);
          new_disk->time = t;

          /* Which lineage is going to be hit ? */
          hr -= log(1./PHYREX_Number_Of_Outgoing_Ldsks(young_disk));
          
          /* young_ldsk = old_disk->ldsk_a[hit_ldsk_idx]; */
          young_ldsk = PHYREX_Random_Select_Outgoing_Ldsk(young_disk);
          old_ldsk   = young_ldsk->prev;
          
          assert(young_disk != old_ldsk->disk);
                    
          if(tree->eval_glnL == YES) new_glnL -= PHYREX_Lk_Range(young_ldsk->disk,old_ldsk->disk,tree);
          
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
          
          
          for(j=0;j<tree->mmod->n_dim;j++)
            {
              new_disk->centr->lonlat[j] = Uni()*(tree->mmod->lim_up->lonlat[j]-tree->mmod->lim_do->lonlat[j])+tree->mmod->lim_do->lonlat[j];
              hr -= log(1./(tree->mmod->lim_up->lonlat[j]-tree->mmod->lim_do->lonlat[j]));
              
              err = NO;
              new_ldsk->coord->lonlat[j] = Rnorm_Trunc(new_disk->centr->lonlat[j],
                                                       1.0*tree->mmod->rad,
                                                       tree->mmod->lim_do->lonlat[j],
                                                       tree->mmod->lim_up->lonlat[j],&err);

              if(err == YES) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);

              hr -= Log_Dnorm_Trunc(new_ldsk->coord->lonlat[j],
                                    new_disk->centr->lonlat[j],
                                    1.0*tree->mmod->rad,
                                    tree->mmod->lim_do->lonlat[j],
                                    tree->mmod->lim_up->lonlat[j],&err);
            }
          
          

          if(tree->eval_glnL == YES)
            {
              new_glnL += PHYREX_Lk_Range(young_ldsk->disk,old_ldsk->disk,tree);
              tree->mmod->c_lnL = new_glnL;
              /* new_glnL = PHYREX_Lk(tree); */
            }
          
          
          n_valid_disks = 0;
          disk = tree->young_disk;
          while(disk->prev)
            {
              if(disk->ldsk && disk->ldsk->n_next == 1) n_valid_disks++;
              disk = disk->prev;
            }

          assert(n_valid_disks);
          hr += log(1./n_valid_disks);
          
          ratio  = (new_glnL - cur_glnL);
          ratio += hr;
          
          ratio = exp(ratio);
          alpha = MIN(1.,ratio);
          
          /* PhyML_Printf("\n+ new_glnL: %f cur_glnL: %f hr: %f alpha: %f",new_glnL,cur_glnL,hr,alpha); */

          /* Always accept move */
          if(tree->mcmc->always_yes == YES && new_glnL > UNLIKELY) alpha = 1.0;
          
          u = Uni();

          assert(isnan(u) == NO && isinf(fabs(u)) == NO);
          
          if(u > alpha) /* Reject */
            {              
              old_ldsk->next[dir_old_young] = young_ldsk;
              young_ldsk->prev = old_ldsk;

              PHYREX_Remove_Disk(new_disk);
              
              PHYREX_Lk_Range(young_ldsk->disk,old_ldsk->disk,tree);
              tree->mmod->c_lnL = cur_glnL;
              /* PHYREX_Lk(tree); */
              
              Free_Disk(new_disk);
              Free_Ldisk(new_ldsk);

            }
          else
            {
              tree->mcmc->acc_move[tree->mcmc->num_move_phyrex_indel_hit_serial]++;
            }
        }

      else /* Remove hit */

        {
          disk = tree->young_disk->prev;
          valid_disks = NULL;
          n_valid_disks = 0;
          do
            {
              if(disk->ldsk && disk->ldsk->n_next == 1)
                {
                  if(!n_valid_disks) valid_disks = (t_dsk **)mCalloc(block,sizeof(t_dsk *));
                  else if(!(n_valid_disks%block)) valid_disks = (t_dsk **)mRealloc(valid_disks,n_valid_disks+block,sizeof(t_dsk *));
                  valid_disks[n_valid_disks] = disk;
                  n_valid_disks++;
                }
              disk = disk->prev;
            }
          while(disk);
          
          if(n_valid_disks == 0) continue;
          
          hr -= log(1./n_valid_disks);
          
          target_disk = valid_disks[Rand_Int(0,n_valid_disks-1)];
          Free(valid_disks);

          assert(target_disk->age_fixed == NO);
          
          old_ldsk   = target_disk->ldsk->prev;
          young_ldsk = target_disk->ldsk->next[0];

          dir_old_young = PHYREX_Get_Next_Direction(target_disk->ldsk,old_ldsk);
          assert(dir_old_young != -1);

          /* Part of the Hastings ratio corresponding to the probability of selecting */
          /* one of target_disk->n_ldsk_a to be hit (reverse move)  */ 
          hr += log(1./target_disk->next->n_ldsk_a);
          

          /* Density for position of the displaced ldsk */
          for(j=0;j<tree->mmod->n_dim;j++)
            {
              hr += log(1./(tree->mmod->lim_up->lonlat[j]-tree->mmod->lim_do->lonlat[j]));

              hr += Log_Dnorm_Trunc(target_disk->ldsk->coord->lonlat[j],
                                    target_disk->centr->lonlat[j],
                                    1.0*tree->mmod->rad,
                                    tree->mmod->lim_do->lonlat[j],
                                    tree->mmod->lim_up->lonlat[j],&err);
            }
          

          assert(target_disk->next);
          
          if(tree->eval_glnL == YES) new_glnL -= PHYREX_Lk_Range(target_disk->next,target_disk->ldsk->prev->disk,tree);
          
          /* New connections between old_ldsk and young_ldsk */
          old_ldsk->next[dir_old_young] = young_ldsk;
          young_ldsk->prev              = old_ldsk;

          hr += log(1./(tree->young_disk->time-T));
          
          PHYREX_Remove_Disk(target_disk);

          if(tree->eval_glnL == YES)
            {
              new_glnL += PHYREX_Lk_Range(target_disk->next,target_disk->ldsk->prev->disk,tree);
              tree->mmod->c_lnL = new_glnL;
              /* new_glnL = PHYREX_Lk(tree); */
            }
          

          ratio  = (new_glnL - cur_glnL);
          ratio += hr;
         
          ratio = exp(ratio);
          alpha = MIN(1.,ratio);

          /* PhyML_Printf("\n- new_glnL: %f cur_glnL: %f hr: %f alpha: %f",new_glnL,cur_glnL,hr,alpha); */
          
          /* Always accept move */
          if(tree->mcmc->always_yes == YES && new_glnL > UNLIKELY) alpha = 1.0;
          
          u = Uni();
          
          assert(isnan(u) == NO && isinf(fabs(u)) == NO);

          if(u > alpha) /* Reject */
            {
              PHYREX_Insert_Disk(target_disk,tree);
              old_ldsk->next[dir_old_young] = target_disk->ldsk;
              young_ldsk->prev              = target_disk->ldsk;
                            
              PHYREX_Lk_Range(target_disk->next,target_disk->ldsk->prev->disk,tree);
              tree->mmod->c_lnL = cur_glnL;
              /* PHYREX_Lk(tree); */
            }
          else
            {
              Free_Ldisk(target_disk->ldsk);
              Free_Disk(target_disk);
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
  t_dsk *disk,*new_disk,*target_disk,**valid_disks,*young_disk,*old_disk;
  int i,j,n_trials,n_valid_disks,block;
  phydbl ratio, u, alpha, hr, type;
  phydbl cur_glnL, new_glnL;
  phydbl log_lk_centr;
  phydbl T,t,pindel;

  cur_glnL     = tree->mmod->c_lnL;
  new_glnL     = tree->mmod->c_lnL;
  hr           = 0.0;
  ratio        = 0.0;
  type         = -1.0;
  n_trials     = (int)(1.+0.1*PHYREX_Total_Number_Of_Intervals(tree));
  T            = PHYREX_Tree_Height(tree);
  pindel       = 0.5;
  disk         = NULL;
  new_disk     = NULL;
  target_disk  = NULL;
  young_disk   = NULL;
  old_disk     = NULL;
  block        = 50;
  
  log_lk_centr = 0.0;
  for(j=0;j<tree->mmod->n_dim;j++) log_lk_centr += log(1./(tree->mmod->lim_up->lonlat[j]-tree->mmod->lim_do->lonlat[j]));

  disk = tree->young_disk->prev;

  for(i=0;i<n_trials;i++)
    {
      tree->mcmc->run_move[tree->mcmc->num_move_phyrex_indel_disk_serial]++;

      cur_glnL = tree->mmod->c_lnL;
      new_glnL = tree->mmod->c_lnL;
      hr       = 0.0;
      ratio    = 0.0;

      type = Uni();
      
      if(type < pindel) /* Insert */
        {
          t = Uni()*(tree->young_disk->time-T) + T;
          disk = tree->young_disk->prev;
          while(disk && disk->time > t) disk = disk->prev;

          young_disk = disk->next;
          old_disk   = disk;
          
          assert(disk->next);
                
          hr -= log(1./(tree->young_disk->time-T));
          hr -= log_lk_centr;
          
          if(tree->eval_glnL == YES) new_glnL -= PHYREX_Lk_Range(young_disk,old_disk,tree);
          
          new_disk = PHYREX_Make_Disk_Event(tree->mmod->n_dim,tree->n_otu);
          PHYREX_Init_Disk_Event(new_disk,tree->mmod->n_dim,tree->mmod);
          new_disk->time = t;
          PHYREX_Insert_Disk(new_disk,tree);
          
          for(j=0;j<tree->mmod->n_dim;j++) new_disk->centr->lonlat[j] = Uni()*(tree->mmod->lim_up->lonlat[j]-tree->mmod->lim_do->lonlat[j])+tree->mmod->lim_do->lonlat[j];

          n_valid_disks = 0;
          disk = tree->young_disk->prev;
          do
            {
              if(!disk->ldsk && disk->age_fixed == NO) n_valid_disks++;
              disk = disk->prev;
            }
          while(disk);
          
          assert(n_valid_disks);
          hr += log(1./n_valid_disks);
          
          if(tree->eval_glnL == YES)
            {
              new_glnL += PHYREX_Lk_Range(young_disk,old_disk,tree);
              tree->mmod->c_lnL = new_glnL;
              /* new_glnL = PHYREX_Lk(tree); */
            }
          
          ratio  = (new_glnL - cur_glnL);
          ratio += hr;
          
          ratio = exp(ratio);
          alpha = MIN(1.,ratio);
          
          /* Always accept move */
          if(tree->mcmc->always_yes == YES && new_glnL > UNLIKELY) alpha = 1.0;
          
          u = Uni();
          
          assert(isnan(u) == NO && isinf(fabs(u)) == NO);

          if(u > alpha) /* Reject */
            {
              PHYREX_Remove_Disk(new_disk);
              Free_Disk(new_disk);
              tree->mmod->c_lnL = cur_glnL;
            }
          else
            {
              tree->mcmc->acc_move[tree->mcmc->num_move_phyrex_indel_disk_serial]++;
            }
        }
      else /* Remove disk */
        {
          disk = tree->young_disk->prev;
          valid_disks = NULL;
          n_valid_disks = 0;
          do
            {
              if(!disk->ldsk && disk->age_fixed == NO)
                {
                  if(!n_valid_disks) valid_disks = (t_dsk **)mCalloc(block,sizeof(t_dsk *));
                  else if(!(n_valid_disks%block)) valid_disks = (t_dsk **)mRealloc(valid_disks,n_valid_disks+block,sizeof(t_dsk *));
                  valid_disks[n_valid_disks] = disk;
                  n_valid_disks++;
                }
              disk = disk->prev;
            }
          while(disk);

          if(n_valid_disks == 0) continue;
          
          assert(n_valid_disks);
          hr -= log(1./n_valid_disks);
          
          target_disk = valid_disks[Rand_Int(0,n_valid_disks-1)];
          Free(valid_disks);
          
          hr += log(1./(tree->young_disk->time-T));
          hr += log_lk_centr;

          assert(target_disk->next->prev);

          if(tree->eval_glnL == YES)
            {
              new_glnL -= PHYREX_Lk_Range(target_disk->next,target_disk->prev,tree);
            }
          
          PHYREX_Remove_Disk(target_disk);
          assert(target_disk->next->prev);

          if(tree->eval_glnL == YES)
            {
              new_glnL += PHYREX_Lk_Range(target_disk->next,target_disk->prev,tree);
              tree->mmod->c_lnL = new_glnL;
              /* new_glnL = PHYREX_Lk(tree); */
            }
          
          ratio  = (new_glnL - cur_glnL);
          ratio += hr;
          
          ratio = exp(ratio);
          alpha = MIN(1.,ratio);
          

          /* Always accept move */
          if(tree->mcmc->always_yes == YES && new_glnL > UNLIKELY) alpha = 1.0;
          
          u = Uni();
          
          assert(isnan(u) == NO && isinf(fabs(u)) == NO);

          if(u > alpha) /* Reject */
            {
              PHYREX_Insert_Disk(target_disk,tree);              
              tree->mmod->c_lnL = cur_glnL;
            }
          else
            {
              tree->mcmc->acc_move[tree->mcmc->num_move_phyrex_indel_disk_serial]++;
              Free_Disk(target_disk);
            }
        }
    }
}
#endif

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
// Add or remove hit but leave corresponding disk unchanged
#ifdef PHYREX
void MCMC_PHYREX_Add_Remove_Jump(t_tree *tree)
{
  t_dsk *disk,*target_disk,**valid_disks;
  int i,j,n_trials,n_valid_disks,block,err;
  phydbl ratio, u, alpha, hr;
  phydbl cur_glnL, new_glnL;
  int target_ldsk_idx,dir_next;
  t_ldsk *new_ldsk,*target_ldsk;
  
  n_trials = (int)(1.+0.1*PHYREX_Total_Number_Of_Intervals(tree));
  block    = 50;

  for(i=0;i<n_trials;++i)
    {      
      tree->mcmc->run_move[tree->mcmc->num_move_phyrex_add_remove_jump]++;


      cur_glnL     = tree->mmod->c_lnL;
      new_glnL     = tree->mmod->c_lnL;
      hr           = 0.0;
      ratio        = 0.0;
      disk         = NULL;
      target_disk  = NULL;
      target_ldsk  = NULL;
      new_ldsk     = NULL;
      
      disk = tree->young_disk->prev;
      valid_disks = NULL;
      n_valid_disks = 0;
      do
        {
          if(disk->prev != NULL &&
             disk->age_fixed == NO &&
             ((disk->ldsk != NULL && disk->ldsk->n_next == 1) || disk->ldsk == NULL))
            {
              if(!n_valid_disks) valid_disks = (t_dsk **)mCalloc(block,sizeof(t_dsk *));
              else if(!(n_valid_disks%block)) valid_disks = (t_dsk **)mRealloc(valid_disks,n_valid_disks+block,sizeof(t_dsk *));
              valid_disks[n_valid_disks] = disk;
              n_valid_disks++;
            }
          disk = disk->prev;
        }
      while(disk);
      
      if(n_valid_disks == 0) return;
      
      target_disk = valid_disks[Rand_Int(0,n_valid_disks-1)];
      Free(valid_disks);

      
      if(target_disk->ldsk == NULL) // Add jump
        {
          target_ldsk_idx = Rand_Int(0,target_disk->n_ldsk_a-1);
          target_ldsk = target_disk->ldsk_a[target_ldsk_idx];
          
          hr -= log(1./(phydbl)target_disk->n_ldsk_a);
          
          if(tree->eval_glnL == YES)
            {
              assert(target_disk->next);
              assert(target_disk->prev);
              new_glnL -= PHYREX_Lk_Range(target_disk,target_ldsk->prev->disk,tree);
            }
          
      
          dir_next = PHYREX_Get_Next_Direction(target_ldsk,target_ldsk->prev);
          assert(dir_next != -1);

          /* Make and initialize new ldsk */
          new_ldsk = PHYREX_Make_Lindisk_Node(tree->mmod->n_dim);
          PHYREX_Init_Lindisk_Node(new_ldsk,target_disk,tree->mmod->n_dim);
          PHYREX_Make_Lindisk_Next(new_ldsk);
          target_disk->ldsk = new_ldsk;
          
          

          new_ldsk->prev = target_ldsk->prev;;
          new_ldsk->next[0] = target_ldsk;
          new_ldsk->prev->next[dir_next] = new_ldsk;          
          target_ldsk->prev = new_ldsk;
          
          
          for(j=0;j<tree->mmod->n_dim;j++)
            {
              err = NO;
              new_ldsk->coord->lonlat[j] = Rnorm_Trunc(target_disk->centr->lonlat[j],
                                                       1.*tree->mmod->rad,
                                                       tree->mmod->lim_do->lonlat[j],
                                                       tree->mmod->lim_up->lonlat[j],&err);
              
              if(err == YES) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
              
              hr -= Log_Dnorm_Trunc(new_ldsk->coord->lonlat[j],
                                    target_disk->centr->lonlat[j],
                                    1.*tree->mmod->rad,
                                    tree->mmod->lim_do->lonlat[j],
                                    tree->mmod->lim_up->lonlat[j],&err);
            }
          
          if(tree->eval_glnL == YES)
            {
              new_glnL += PHYREX_Lk_Range(target_disk,new_ldsk->prev->disk,tree);
              tree->mmod->c_lnL = new_glnL;
              /* new_glnL = PHYREX_Lk(tree); */
            }
          
          ratio  = (new_glnL - cur_glnL);
          ratio += hr;
          
          ratio = exp(ratio);
          alpha = MIN(1.,ratio);
          
          /* Always accept move */
          if(tree->mcmc->always_yes == YES && new_glnL > UNLIKELY) alpha = 1.0;
          
          u = Uni();
          
          assert(isnan(u) == NO && isinf(fabs(u)) == NO);
          
          if(u > alpha) /* Reject */
            {              
              target_ldsk->prev = new_ldsk->prev;
              new_ldsk->prev->next[dir_next] = target_ldsk;
              target_disk->ldsk = NULL;
              
              PHYREX_Lk_Range(target_disk,target_ldsk->prev->disk,tree);
              tree->mmod->c_lnL = cur_glnL;
              /* PHYREX_Lk(tree); */
              
              Free_Ldisk(new_ldsk);
            }
          else
            {
              tree->mcmc->acc_move[tree->mcmc->num_move_phyrex_add_remove_jump]++;
            }      
        }
      else // Remove jump
        {
          target_ldsk = target_disk->ldsk;
          
          hr += log(1./(phydbl)target_disk->n_ldsk_a);
          
          if(tree->eval_glnL == YES)
            {
              assert(target_disk->next);
              assert(target_disk->prev);
              new_glnL -= PHYREX_Lk_Range(target_disk,target_ldsk->prev->disk,tree);
            }
          
          for(j=0;j<tree->mmod->n_dim;j++)
            {
              hr += Log_Dnorm_Trunc(target_ldsk->coord->lonlat[j],
                                    target_disk->centr->lonlat[j],
                                    1.*tree->mmod->rad,
                                    tree->mmod->lim_do->lonlat[j],
                                    tree->mmod->lim_up->lonlat[j],&err);
            }
          
          dir_next = PHYREX_Get_Next_Direction(target_ldsk,target_ldsk->prev);
          assert(dir_next != -1);
          
          target_ldsk->prev->next[dir_next] = target_ldsk->next[0];
          target_ldsk->next[0]->prev = target_ldsk->prev;
          target_disk->ldsk = NULL;
          
          
          if(tree->eval_glnL == YES)
            {
              new_glnL += PHYREX_Lk_Range(target_disk,target_ldsk->prev->disk,tree);
              tree->mmod->c_lnL = new_glnL;
              /* new_glnL = PHYREX_Lk(tree); */
            }
          
          ratio  = (new_glnL - cur_glnL);
          ratio += hr;
          
          ratio = exp(ratio);
          alpha = MIN(1.,ratio);
          
          /* Always accept move */
          if(tree->mcmc->always_yes == YES && new_glnL > UNLIKELY) alpha = 1.0;
          
          u = Uni();
          
          assert(isnan(u) == NO && isinf(fabs(u)) == NO);
          
          if(u > alpha) /* Reject */
            {              
              target_ldsk->prev->next[dir_next] = target_ldsk;
              target_ldsk->next[0]->prev = target_ldsk;              
              target_disk->ldsk = target_ldsk;
              
              PHYREX_Lk_Range(target_disk,target_ldsk->prev->disk,tree);
              tree->mmod->c_lnL = cur_glnL;
              /* PHYREX_Lk(tree); */
            }
          else
            {              
              Free_Ldisk(target_ldsk);
              tree->mcmc->acc_move[tree->mcmc->num_move_phyrex_add_remove_jump]++;
            }
        }
    }
}
#endif

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

#ifdef PHYREX
void MCMC_PHYREX_Simulate_Backward(t_tree *tree)
{
  phydbl u,alpha,ratio,hr,t;
  phydbl Told,Tyoung;
  phydbl cur_alnL, new_alnL;
  phydbl cur_glnL, new_glnL;
  phydbl cur_rlnL, new_rlnL;
  t_dsk  *disk,*bkp_disk,*target_disk;
  t_ldsk **bkp_ldsk;
  int i;
  t_node *n;

  tree->mcmc->run_move[tree->mcmc->num_move_phyrex_sim]++;

  disk        = NULL;
  bkp_disk    = NULL;
  target_disk = NULL;
  bkp_ldsk    = NULL;

  new_alnL      = tree->c_lnL;
  cur_alnL      = tree->c_lnL;
  cur_glnL      = tree->mmod->c_lnL;
  new_glnL      = tree->mmod->c_lnL;
  new_rlnL       = tree->rates->c_lnL_rates;
  cur_rlnL       = tree->rates->c_lnL_rates;
  hr            = 0.0;
  ratio         = 0.0;
  t             = 0.0;


  Get_Node_Ranks_From_Tip_Times(tree);
  
  // Get to the oldest sampled disk
  n = tree->a_nodes[0];
  while(n->rk_prev && n->rk_prev->tax == YES) n = n->rk_prev;
  Tyoung = n->ldsk->disk->time;
  

  
  /* Chop off the tree at time t and simulate upwards from here */
  Told = PHYREX_Tree_Height(tree);

  t = Uni()*(Tyoung-Told) + Told;
  
  disk = tree->young_disk->prev;
  while(disk && disk->time > t) disk = disk->prev;
  target_disk = disk->next;

  hr -= log(fabs((target_disk->prev->time - target_disk->time)/(Tyoung-Told)));

  bkp_disk = target_disk->prev;

  bkp_ldsk = (t_ldsk **)mCalloc(target_disk->n_ldsk_a,sizeof(t_ldsk *));
  for(i=0;i<target_disk->n_ldsk_a;i++) bkp_ldsk[i] = target_disk->ldsk_a[i]->prev;

  
  PHYREX_Simulate_Backward_Core(target_disk,NO,tree);
  PHYREX_Ldsk_To_Tree(tree);

  
  Told = PHYREX_Tree_Height(tree);
  hr += log(fabs((target_disk->prev->time - target_disk->time)/(Tyoung-Told)));


  if(tree->eval_alnL == YES) new_alnL = Lk(NULL,tree);
  if(tree->eval_rlnL == YES) new_rlnL = RATES_Lk_Rates(tree);
      
  ratio += (new_alnL - cur_alnL);
  ratio += (new_rlnL - cur_rlnL);
  ratio += hr;
  
  ratio = exp(ratio);
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
      for(i=0;i<target_disk->n_ldsk_a;i++) target_disk->ldsk_a[i]->prev = bkp_ldsk[i];
      
      if(tree->mmod->safe_phyrex == YES)
        {
          new_alnL = Lk(NULL,tree);
          if(Are_Equal(new_alnL,cur_alnL,1.E-3) == NO)
            {
              PhyML_Fprintf(stderr,"\n. new_alnL: %f cur_alnL: %f",new_alnL,cur_alnL);
              Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
            }
        
          
          new_glnL = PHYREX_Lk(tree);
          
          if(Are_Equal(new_glnL,cur_glnL,1.E-3) == NO)
            {
              PhyML_Fprintf(stderr,"\n. new_glnL: %f cur_glnL: %f",new_glnL,cur_glnL);
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
void MCMC_PHYREX_Ldsk_Tip_To_Root(t_tree *tree)
{
  phydbl u,alpha,ratio,hr,c;
  phydbl cur_glnL, new_glnL;
  t_dsk  *disk,**all_disks;
  int i,j,err,n_all_disks,block,n_move_ldsk,*permut;
  t_geo_coord *mean;
  phydbl K;
  
  block       = 100;
  all_disks   = NULL;
  n_all_disks = 0;
  c           = -1.;
  mean        = NULL;
  K           = tree->mcmc->tune_move[tree->mcmc->num_move_phyrex_ldsk_tip_to_root];
  
  disk = tree->young_disk->prev;
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
  
  n_move_ldsk = (int)(n_all_disks);
  
  permut = Permutate(n_all_disks);
  cur_glnL = tree->mmod->c_lnL;
  new_glnL = tree->mmod->c_lnL;
  for(i=0;i<n_move_ldsk;i++)
    {
      tree->mcmc->run_move[tree->mcmc->num_move_phyrex_ldsk_tip_to_root]++;
      
      disk = all_disks[permut[i]];
            
      hr       = 0.0;
      ratio    = 0.0;
      new_glnL = cur_glnL;   
      
      PHYREX_Store_Geo_Coord(disk->ldsk->coord);
      
      if(tree->eval_glnL == YES)
        {
          if(disk->ldsk->prev != NULL) new_glnL -= PHYREX_Lk_Range(disk,disk->ldsk->prev->disk,tree);
          else new_glnL -= PHYREX_Lk_Core(disk,tree);
        }
      
      mean = PHYREX_Mean_Next_Loc(disk->ldsk,tree);

      for(j=0;j<tree->mmod->n_dim;j++)
        {
          if(disk->ldsk->prev != NULL)
            c = 0.5*(mean->lonlat[j] + disk->ldsk->prev->coord->lonlat[j]);
          else
            c = mean->lonlat[j];
          
          err = NO;
          disk->ldsk->coord->lonlat[j] =
            Rnorm_Trunc(c,
                        (1./K)*(tree->mmod->lim_up->lonlat[j] - tree->mmod->lim_do->lonlat[j]),
                        tree->mmod->lim_do->lonlat[j],
                        tree->mmod->lim_up->lonlat[j],&err);
          
          if(err == YES) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);

          hr -= Log_Dnorm_Trunc(disk->ldsk->coord->lonlat[j],
                                c,
                                (1./K)*(tree->mmod->lim_up->lonlat[j] - tree->mmod->lim_do->lonlat[j]),
                                tree->mmod->lim_do->lonlat[j],
                                tree->mmod->lim_up->lonlat[j],&err);
          
          hr += Log_Dnorm_Trunc(disk->ldsk->coord->cpy->lonlat[j],
                                c,
                                (1./K)*(tree->mmod->lim_up->lonlat[j] - tree->mmod->lim_do->lonlat[j]),
                                tree->mmod->lim_do->lonlat[j],
                                tree->mmod->lim_up->lonlat[j],&err);

        }

      Free_Geo_Coord(mean);
      
     if(tree->eval_glnL == YES)
       {
          if(disk->ldsk->prev != NULL) new_glnL += PHYREX_Lk_Range(disk,disk->ldsk->prev->disk,tree);
          else new_glnL += PHYREX_Lk_Core(disk,tree);
          tree->mmod->c_lnL = new_glnL;
       }
     
      ratio += (new_glnL - cur_glnL);
      ratio += hr;
      
      ratio = exp(ratio);
      alpha = MIN(1.,ratio);
      
      u = Uni();
      
      assert(isnan(u) == NO && isinf(fabs(u)) == NO);

      if(u > alpha) /* Reject */
        {
          PHYREX_Restore_Geo_Coord(disk->ldsk->coord);          
          tree->mmod->c_lnL = cur_glnL;            
        }
      else
        {
          cur_glnL = new_glnL;
          tree->mcmc->acc_move[tree->mcmc->num_move_phyrex_ldsk_tip_to_root]++;
        }
    }

  Free(permut);
  Free(all_disks);


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
