/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

/* Routines for molecular clock trees and molecular dating */


#include "rates.h"

#ifdef RWRAPPER
#include <R.h>
#endif

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl RATES_Lk(t_tree *tree)
{
  int err;
  
  if(tree->eval_rlnL == NO) return UNLIKELY;

  tree->rates->c_lnL  = .0;

  RATES_Lk_Pre(tree->n_root,tree->n_root->v[2],NULL,tree);
  RATES_Lk_Pre(tree->n_root,tree->n_root->v[1],NULL,tree);
  
  err = NO;

  switch(tree->rates->model_id)
    {
    case LOGNORMAL :
      {
        tree->rates->c_lnL += Log_Dnorm(log(tree->rates->br_r[tree->n_root->num]),-tree->rates->nu*tree->rates->nu/2.,tree->rates->nu,&err);
        tree->rates->c_lnL -= log(tree->rates->br_r[tree->n_root->num]);        
        break;
      }
    case STRICTCLOCK :
      {
        tree->rates->c_lnL += 0.0;
	break;
      }      
    default : 
      {
	PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d\n",__FILE__,__LINE__);
	Warn_And_Exit("");
      }
    }

  if(tree->rates->clock_r_has_prior == YES)
    {
      tree->rates->c_lnL += RATES_Clock_R_Prior(tree);      
    }
  
  if(isnan(tree->rates->c_lnL) || err == YES) assert(false);
  
  return tree->rates->c_lnL;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl RATES_Clock_R_Prior(t_tree *tree)
{
  phydbl mean,sd,lnP;
  int err;

  if(tree->rates->clock_r_has_prior == NO) return(UNLIKELY);
  
  err = NO;
  lnP = 0.0;
  
  mean = log(tree->rates->clock_r_prior_mean);
  sd = sqrt(tree->rates->clock_r_prior_var);
  
  lnP += Log_Dnorm(log(tree->rates->clock_r),mean-sd*sd/2.,sd,&err);
  lnP -= log(tree->rates->clock_r);

  /* PhyML_Printf("\n. prior mean: %f var: %f || mean: %f sd: %f [%f;%f] lnL: %f", */
  /*              tree->rates->clock_r_prior_mean, */
  /*              tree->rates->clock_r_prior_var, */
  /*              mean, */
  /*              sd, */
  /*              tree->rates->clock_r, */
  /*              log(tree->rates->clock_r), */
  /*              lnP); */

  return(lnP);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl RATES_Prior(t_tree *tree)
{
  tree->rates->c_lnP = 0.0;
  return(tree->rates->c_lnP);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Lk_Pre(t_node *a, t_node *d, t_edge *b, t_tree *tree)
{
  int i;
  phydbl log_dens,mu_a,mu_d,r_a,r_d,dt_a,dt_d;
  int n_a,n_d;

  log_dens = -1.;
  
  if(d->anc != a)
    {
      PhyML_Fprintf(stderr,"\n. d=%d d->anc=%d a=%d root=%d",d->num,d->anc->num,a->num,tree->n_root->num);
      PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d\n",__FILE__,__LINE__);
      assert(FALSE);
    }

  dt_a = -1.;
  if(a != tree->n_root) dt_a = tree->times->nd_t[a->num] - tree->times->nd_t[a->anc->num];
  
  mu_a = tree->rates->br_r[a->num];
  r_a  = tree->rates->nd_r[a->num];
  n_a  = tree->times->n_jps[a->num];

  dt_d = FABS(tree->times->nd_t[d->num] - tree->times->nd_t[a->num]);
  mu_d = tree->rates->br_r[d->num];
  r_d  = tree->rates->nd_r[d->num];
  n_d  = tree->times->n_jps[d->num];

  log_dens = RATES_Lk_Core(mu_a,mu_d,r_a,r_d,n_a,n_d,dt_a,dt_d,tree);
  tree->rates->c_lnL += log_dens;

  if(isnan(tree->rates->c_lnL))
    {
      PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d\n",__FILE__,__LINE__);
      MCMC_Print_Param(tree->mcmc,tree);
      Exit("\n");
    }

  tree->rates->triplet[a->num] += log_dens;

  if(d->tax) return;
  else
    {
      for(i=0;i<3;i++)
	{
	  if((d->v[i] != a) && (d->b[i] != tree->e_root))
	    {
	      RATES_Lk_Pre(d,d->v[i],d->b[i],tree);
	    }
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl RATES_Lk_Core(phydbl br_r_a, phydbl br_r_d, phydbl nd_r_a, phydbl nd_r_d, int n_a, int n_d, phydbl dt_a, phydbl dt_d, t_tree *tree)
{
  phydbl log_dens,mean,sd,min_r, max_r,cr;
  

  log_dens  = UNLIKELY;
  mean = sd = -1.;
  min_r     = tree->rates->min_rate;
  max_r     = tree->rates->max_rate;
  cr        = tree->rates->clock_r;
  
  if(br_r_d > tree->rates->max_rate) return UNLIKELY;
  if(br_r_d < tree->rates->min_rate) return UNLIKELY;
  
  /* dt_d = MAX(0.5,dt_d); // We give only one decimal when printing out node heights. It is therefore a fair approximation */

  switch(tree->rates->model_id)
    {
    case THORNE :
      {
        int err;
        phydbl log_br_r_d = log(br_r_d);
        phydbl log_br_r_a = log(br_r_a);

        /* log_dens = Log_Dnorm_Trunc(log_br_r_d,log_br_r_a,sqrt(tree->rates->nu*dt_d),log(tree->rates->min_rate),log(tree->rates->max_rate),&err); */
        log_dens = Log_Dnorm(log_br_r_d,log_br_r_a,sqrt(tree->rates->nu*dt_d),&err);
        log_dens -= log_br_r_d;

        break;
      }
    case GUINDON :
      {
	int err;	

	min_r = tree->rates->min_rate;
	max_r = tree->rates->max_rate;

	nd_r_d = log(nd_r_d*cr);
	nd_r_a = log(nd_r_a*cr);
        min_r  = log(min_r*cr);
	max_r  = log(max_r*cr);

	sd   = sqrt(tree->rates->nu*dt_d);
	mean = nd_r_a - .5*sd*sd;

        // log(density(log(Rate)))
 	log_dens = Log_Dnorm_Trunc(nd_r_d,mean,sd,min_r,max_r,&err);
        // log(density(Rate))
        log_dens -= log(exp(nd_r_d)/cr);
        
	if(err)
	  {
	    PhyML_Fprintf(stderr,"\n. Run: %d",tree->mcmc->run);
	    PhyML_Fprintf(stderr,"\n. br_r_d=%f mean=%f sd=%f min_r=%f max_r=%f dt_d=%f",br_r_d,mean,sd,min_r,max_r,dt_d);
	    PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d\n",__FILE__,__LINE__);
	    Exit("\n");
	  }
	break;
      }
    case LOGNORMAL :
      {
        int err;
        phydbl log_br_r_d = log(br_r_d);
        log_dens = Log_Dnorm(log_br_r_d,-tree->rates->nu*tree->rates->nu/2.,tree->rates->nu,&err);
        log_dens -= log_br_r_d;
        break;
      }
    case STRICTCLOCK :
      {
	log_dens = 0.0;
	break;
      }

    default : 
      {
	PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d\n",__FILE__,__LINE__);
	Warn_And_Exit("");
      }
    }

  if(isnan(log_dens))
    {
      PhyML_Fprintf(stderr,"\n. Run=%4d br_r_d=%f br_r_a=%f dt_d=%f dt_a=%f nu=%f log_dens=%G sd=%f mean=%f\n",
                    tree->mcmc->run,
                    br_r_d,br_r_a,dt_d,dt_a,tree->rates->nu,log_dens,
                    sd,mean);
      assert(false);
    }

  return log_dens;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl RATES_Lk_Change_One_Rate(t_node *d, phydbl new_rate, t_tree *tree)
{
  tree->rates->br_r[d->num] = new_rate;
  RATES_Update_Triplet(d,tree);
  RATES_Update_Triplet(d->anc,tree);
  return(tree->rates->c_lnL);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl RATES_Lk_Change_One_Time(t_node *n, phydbl new_t, t_tree *tree)
{  
  if(n == tree->n_root)
    {
      PhyML_Fprintf(stderr,"\n. Moving the time of the root t_node is not permitted.");
      PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
  else
    {
      int i;
      
      tree->times->nd_t[n->num] = new_t;

      RATES_Update_Triplet(n,tree);
      
      for(i=0;i<3;i++)
	{
	  if(n->b[i] != tree->e_root) RATES_Update_Triplet(n->v[i],tree);
	  else RATES_Update_Triplet(tree->n_root,tree);
	}
    }
  return(tree->rates->c_lnL);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Update_Triplet(t_node *n, t_tree *tree)
{
  phydbl curr_triplet,new_triplet;
  phydbl dt0,dt1,dt2;
  phydbl mu1_mu0,mu2_mu0;
  phydbl mu0,mu1,mu2;
  phydbl r0,r1,r2;
  int n0,n1,n2;
  int i;
  t_node *v1,*v2;

  if(n->tax) return;

  curr_triplet = tree->rates->triplet[n->num];

  dt0 = dt1 = dt2 = -100.0;
  r0 = r1 = r2 = 0.0;

  if(n == tree->n_root)
    {
      phydbl log_dens;
      
      log_dens = 0.0;

      dt0 = tree->times->nd_t[tree->n_root->v[2]->num] - tree->times->nd_t[tree->n_root->num];
      dt1 = tree->times->nd_t[tree->n_root->v[1]->num] - tree->times->nd_t[tree->n_root->num];
      
      mu0 = tree->rates->br_r[tree->n_root->v[2]->num];
      mu1 = tree->rates->br_r[tree->n_root->v[1]->num];

      r0 = tree->rates->nd_r[tree->n_root->v[2]->num];
      r1 = tree->rates->nd_r[tree->n_root->v[1]->num];
      
      n0  = tree->times->n_jps[tree->n_root->v[2]->num];
      n1  = tree->times->n_jps[tree->n_root->v[1]->num];

      switch(tree->rates->model_id)
	{
	case COMPOUND_COR : case COMPOUND_NOCOR : 
	  {
	    log_dens  = RATES_Dmu(mu0,n0,dt0,tree->rates->nu,1./tree->rates->nu,tree->rates->lexp,0,1);
	    log_dens *= RATES_Dmu(mu1,n1,dt1,tree->rates->nu,1./tree->rates->nu,tree->rates->lexp,0,1);
	    log_dens  = log(log_dens);
	    break;
	  }
	case EXPONENTIAL : 
	  {
	    log_dens = Dexp(mu0,tree->rates->lexp) * Dexp(mu1,tree->rates->lexp);
	    log_dens = log(log_dens);
	    break;
	  }
	case LOGNORMAL :
	  {
	    log_dens = Dgamma(mu0,tree->rates->nu,1./tree->rates->nu) * Dgamma(mu1,tree->rates->nu,1./tree->rates->nu);
	    log_dens = log(log_dens);
	    break;
	  }
	case THORNE :
	  {
	    int err;
	    phydbl mean0,sd0;
	    phydbl mean1,sd1;

	    
	    sd0 = SQRT(tree->rates->nu*dt0);
	    mean0 = 1.0;

	    sd1 = SQRT(tree->rates->nu*dt1);
	    mean1 = 1.0;

	    log_dens = 
	      Log_Dnorm_Trunc(mu0,mean0,sd0,tree->rates->min_rate,tree->rates->max_rate,&err) + 
	      Log_Dnorm_Trunc(mu1,mean1,sd1,tree->rates->min_rate,tree->rates->max_rate,&err); 

	    break;
	  }
	case GUINDON :
	  {
	    PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d\n",__FILE__,__LINE__);
	    Exit("\n. Not implemented yet.\n");
	    break;
	  }
	default :
	  {
	    Exit("\n. Model not implemented yet.\n");
	    break;
	  }
	}
      new_triplet = log_dens;

      if(isnan(log_dens) || isinf(FABS(log_dens)))
	{
	  PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d\n",__FILE__,__LINE__);
	  MCMC_Print_Param(tree->mcmc,tree);
	  Exit("\n");
	}
    }
  else
    {
      mu0 = mu1 = mu2 = -1.;
      n0 = n1 = n2 = -1;

      mu0 = tree->rates->br_r[n->num];
      dt0 = FABS(tree->times->nd_t[n->num] - tree->times->nd_t[n->anc->num]);
      n0  = tree->times->n_jps[n->num];
      r0 = tree->rates->nd_r[n->num];

      v1 = v2 = NULL;
      for(i=0;i<3;i++)
	{
	  if((n->v[i] != n->anc) && (n->b[i] != tree->e_root))
	    {
	      if(!v1)
		{
		  v1  = n->v[i]; 
		  mu1 = tree->rates->br_r[v1->num];
		  dt1 = FABS(tree->times->nd_t[v1->num] - tree->times->nd_t[n->num]);
		  n1  = tree->times->n_jps[v1->num];
		  r1  = tree->rates->nd_r[v1->num];
		}
	      else
		{
		  v2  = n->v[i]; 
		  mu2 = tree->rates->br_r[v2->num];
		  dt2 = FABS(tree->times->nd_t[v2->num] - tree->times->nd_t[n->num]);
		  n2  = tree->times->n_jps[v2->num];
		  r2  = tree->rates->nd_r[v2->num];
		}
	    }
	}
 
      mu1_mu0 = RATES_Lk_Core(mu0,mu1,r0,r1,n0,n1,dt0,dt1,tree);
      mu2_mu0 = RATES_Lk_Core(mu0,mu2,r0,r1,n0,n2,dt0,dt2,tree);
      
      new_triplet = mu1_mu0 + mu2_mu0;
    }

  tree->rates->c_lnL = tree->rates->c_lnL + new_triplet - curr_triplet;
  tree->rates->triplet[n->num] = new_triplet;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Print_Triplets(t_tree *tree)
{
  int i;
  for(i=0;i<2*tree->n_otu-1;++i) PhyML_Printf("\n. Node %3d t=%f",i,tree->rates->triplet[i]);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Print_Rates(t_tree *tree)
{
  RATES_Print_Rates_Pre(tree->n_root,tree->n_root->v[2],NULL,tree);
  RATES_Print_Rates_Pre(tree->n_root,tree->n_root->v[1],NULL,tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Copy_Rate_Struct(t_rate *from, t_rate *to, int n_otu)
{
  int i;

  to->lexp = from->lexp;
  to->alpha = from->alpha;
  to->less_likely = from->less_likely;
  
  to->nu = from->nu;
  to->min_nu = from->min_nu;
  to->max_nu = from->max_nu;

  to->clock_r_has_prior = from->clock_r_has_prior;
  to->clock_r_prior_mean = from->clock_r_prior_mean;  
  to->clock_r_prior_var = from->clock_r_prior_var;  
  
  to->min_rate = from->min_rate;
  to->max_rate = from->max_rate;

  to->c_lnL1 = from->c_lnL1;
  to->c_lnL2 = from->c_lnL2;

  to->c_lnL = from->c_lnL;
  
  to->clock_r = from->clock_r;
  to->min_clock = from->min_clock;
  to->max_clock = from->max_clock;

  to->lbda_nu = from->lbda_nu;
  to->min_dt = from->min_dt;
  to->step_rate = from->step_rate;
  to->true_tree_size = from->true_tree_size;
  to->p_max = from->p_max;
  to->norm_fact = from->norm_fact;

  to->adjust_rates = from->adjust_rates;
  to->use_rates = from->use_rates;
  to->bl_from_rt = from->bl_from_rt;
  to->approx = from->approx;
  to->model_id = from->model_id;
  to->is_allocated = from->is_allocated;
  to->met_within_gibbs = from->met_within_gibbs;
  
  to->update_mean_l = from->update_mean_l;
  to->update_cov_l = from->update_cov_l;
  
  to->br_r_recorded = from->br_r_recorded;
  
  to->log_K_cur = from->log_K_cur;
  to->cur_comb_numb = from->cur_comb_numb;

  to->norm_fact = from->norm_fact;
  
  for(i=0;i<2*n_otu-1;++i) to->nd_r[i] = from->nd_r[i];
  for(i=0;i<2*n_otu-1;++i) to->br_r[i] = from->br_r[i];
  for(i=0;i<2*n_otu-1;++i) to->buff_br_r[i] = from->buff_br_r[i];
  for(i=0;i<2*n_otu-1;++i) to->buff_nd_r[i] = from->buff_nd_r[i];
  for(i=0;i<2*n_otu-1;++i) to->true_r[i] = from->true_r[i];
  for(i=0;i<2*n_otu-2;++i) to->dens[i] = from->dens[i];
  for(i=0;i<2*n_otu-1;++i) to->triplet[i] = from->triplet[i];
  for(i=0;i<(2*n_otu-2)*(2*n_otu-2);++i) to->cov_l[i] = from->cov_l[i];
  for(i=0;i<(2*n_otu-2)*(2*n_otu-2);++i) to->invcov[i] = from->invcov[i];
  for(i=0;i<2*n_otu-2;++i) to->mean_l[i] = from->mean_l[i];
  for(i=0;i<2*n_otu-2;++i) to->ml_l[i] = from->ml_l[i];
  for(i=0;i<2*n_otu-2;++i) to->cur_l[i] = from->cur_l[i];
  for(i=0;i<2*n_otu-3;++i) to->u_ml_l[i] = from->u_ml_l[i];
  for(i=0;i<2*n_otu-3;++i) to->u_cur_l[i] = from->u_cur_l[i];
  for(i=0;i<(2*n_otu-2)*(2*n_otu-2);++i) to->cov_r[i] = from->cov_r[i];
  for(i=0;i<2*n_otu-2;++i) to->cond_var[i] = from->cond_var[i];
  for(i=0;i<2*n_otu-2;++i) to->mean_r[i] = from->mean_r[i];
  for(i=0;i<(2*n_otu-1)*(2*n_otu-1);++i) to->lca[i] = from->lca[i];
  for(i=0;i<(2*n_otu-2)*(2*n_otu-2);++i) to->reg_coeff[i] = from->reg_coeff[i];
  for(i=0;i<(2*n_otu-2)*(6*n_otu-9);++i) to->trip_reg_coeff[i] = from->trip_reg_coeff[i];
  for(i=0;i<(2*n_otu-2)*9;++i) to->trip_cond_cov[i] = from->trip_cond_cov[i];
  for(i=0;i<2*n_otu;++i) to->_2n_vect1[i] = from->_2n_vect1[i];
  for(i=0;i<2*n_otu;++i) to->_2n_vect2[i] = from->_2n_vect2[i];
  for(i=0;i<2*n_otu;++i) to->_2n_vect3[i] = from->_2n_vect3[i];
  for(i=0;i<2*n_otu;++i) to->_2n_vect4[i] = from->_2n_vect4[i];
  for(i=0;i<2*n_otu;++i) to->_2n_vect5[i] = from->_2n_vect5[i];
  for(i=0;i<4*n_otu*n_otu;++i) to->_2n2n_vect1[i] = from->_2n2n_vect1[i];
  for(i=0;i<4*n_otu*n_otu;++i) to->_2n2n_vect2[i] = from->_2n2n_vect2[i];
  for(i=0;i<2*n_otu-1;++i) to->br_do_updt[i] = from->br_do_updt[i];
  for(i=0;i<2*n_otu-1;++i) to->cur_gamma_prior_mean[i] = from->cur_gamma_prior_mean[i];
  for(i=0;i<2*n_otu-1;++i) to->cur_gamma_prior_var[i] = from->cur_gamma_prior_var[i];
  for(i=0;i<2*n_otu-1;++i) to->n_tips_below[i] = from->n_tips_below[i];
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Duplicate_Calib_Struct(t_tree *from, t_tree *to)
{
  int i,j;

  to->times->n_cal = from->times->n_cal;
  for(i=0;i<from->times->n_cal;i++)
    {
      to->times->a_cal[i] = Duplicate_Calib(from->times->a_cal[i]);
      for(j=0;j<from->times->a_cal[i]->clade_list_size;j++)
        Init_Target_Tip(to->times->a_cal[i]->clade_list[j],to);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Print_Rates_Pre(t_node *a, t_node *d, t_edge *b, t_tree *tree)
{  
  if((d == tree->n_root->v[2] && d->tax) || (d == tree->n_root->v[1] && d->tax))
    PhyML_Printf("\n. a=%3d ++d=%3d rate=%12f t_left=%12f t_rght=%12f ml=%12f l=%12f %12f",
	   a->num,d->num,
	   tree->rates->br_r[d->num],
	   tree->times->nd_t[a->num],tree->times->nd_t[d->num],
	   tree->rates->ml_l[d->num],
	   tree->rates->cur_l[d->num],
	   (tree->times->nd_t[d->num]-tree->times->nd_t[a->num])*tree->rates->clock_r*tree->rates->br_r[d->num]);
  
  else if((d == tree->n_root->v[2]) || (d == tree->n_root->v[1]))
    PhyML_Printf("\n. a=%3d __d=%3d rate=%12f t_left=%12f t_rght=%12f ml=%12f l=%12f %12f",
	   a->num,d->num,
	   tree->rates->br_r[d->num],
	   tree->times->nd_t[a->num],tree->times->nd_t[d->num],
	   tree->rates->ml_l[d->num],
	   tree->rates->cur_l[d->num],
	   (tree->times->nd_t[d->num]-tree->times->nd_t[a->num])*tree->rates->clock_r*tree->rates->br_r[d->num]);
  else 
    PhyML_Printf("\n. a=%3d   d=%3d rate=%12f t_left=%12f t_rght=%12f ml=%12f l=%12f %12f",
	   a->num,d->num,
	   tree->rates->br_r[d->num],
	   tree->times->nd_t[a->num],tree->times->nd_t[d->num],
	   tree->rates->ml_l[d->num],
	   tree->rates->cur_l[d->num],
	   (tree->times->nd_t[d->num]-tree->times->nd_t[a->num])*tree->rates->clock_r*tree->rates->br_r[d->num]);
  
  if(d->tax) return;
  else
    {
      int i;

      for(i=0;i<3;i++) 
	{
	  if((d->v[i] != a) && (d->b[i] != tree->e_root))
	    {
	      RATES_Print_Rates_Pre(d,d->v[i],d->b[i],tree);
	    }
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl RATES_Average_Rate(t_tree *tree)
{
  int i;
  phydbl sum;
  sum = 0.0;
  for(i=0;i<2*tree->n_otu-2;++i) sum += tree->rates->br_r[i];
  return sum/(2*tree->n_otu-2);
}
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl RATES_Average_Substitution_Rate(t_tree *tree)
{
  phydbl l,dt;
  int i;

  l = 0.0;
  dt = 0.0;
  for(i=0;i<2*tree->n_otu-1;++i)
    {
      if(tree->a_nodes[i] != tree->n_root)
        {
          dt += FABS(tree->times->nd_t[tree->a_nodes[i]->num] - tree->times->nd_t[tree->a_nodes[i]->anc->num]); 
          l  += tree->rates->cur_l[tree->a_nodes[i]->num];
        }
    }
      
  return(l/dt);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl RATES_Check_Mean_Rates_True(t_tree *tree)
{
  phydbl sum;
  int i;
  
  sum = 0.0;
  for(i=0;i<2*tree->n_otu-2;++i) sum += tree->rates->true_r[i];
  return(sum/(phydbl)(2*tree->n_otu-2));
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int RATES_Check_Node_Times(t_tree *tree)
{
  int err;
  err = NO;
  RATES_Check_Node_Times_Pre(tree->n_root,tree->n_root->v[2],&err,tree);
  RATES_Check_Node_Times_Pre(tree->n_root,tree->n_root->v[1],&err,tree);
  return err;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Check_Node_Times_Pre(t_node *a, t_node *d, int *err, t_tree *tree)
{
  if((tree->times->nd_t[d->num] < tree->times->nd_t[a->num]) || (FABS(tree->times->nd_t[d->num] - tree->times->nd_t[a->num]) < 1.E-20))
    {
      PhyML_Printf("\n. a->t=%f d->t=%f",tree->times->nd_t[a->num],tree->times->nd_t[d->num]);
      PhyML_Printf("\n. a->t_prior_min=%f a->t_prior_max=%f",tree->times->t_prior_min[a->num],tree->times->t_prior_max[a->num]);
      PhyML_Printf("\n. d->t_prior_min=%f d->t_prior_max=%f",tree->times->t_prior_min[d->num],tree->times->t_prior_max[d->num]);
      *err = YES;
    }
  if(d->tax) return;
  else
    {
      int i;

      for(i=0;i<3;i++)
	if((d->v[i] != a) && (d->b[i] != tree->e_root))	  
	  RATES_Check_Node_Times_Pre(d,d->v[i],err,tree);
    }
}
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Bracket_N_Jumps(int *up, int *down, phydbl param)
{
  phydbl cdf,eps,a,b,c;
  int step;

  step = 10;
  eps = 1.E-10;
  cdf = 0.0;
  c = 1;
  
  while(cdf < 1.-eps)
    {
      c = (int)FLOOR(c * step);
      cdf = Ppois(c,param);      
    }
  
  a = 0.0;
  b = (c-a)/2.;
  step = 0;
  do
    {
      step++;
      cdf = Ppois(b,param);
      if(cdf < eps) a = b;
      else 
	{
	  break;
	}
      b = (c-a)/2.;
    }
  while(step < 1000);
  
  if(step == 1000)
    {
      PhyML_Fprintf(stderr,"\n. a=%f b=%f c=%f param=%f",a,b,c,param);
      PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
  *up = c;
  *down = a;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* 
   mu   : average rate of the time period dt
   dt   : time period to be considered
   a    : rate at a given time point is gamma distributed. a is the shape parameter
   b    : rate at a given time point is gamma distributed. b is the scale parameter
   lexp : the number of rate switches is Poisson distributed with parameter lexp * dt
*/ 
/* compute f(mu;dt,a,b,lexp), the probability density of mu. We need to integrate over the
   possible number of jumps (n) during the time interval dt */
phydbl RATES_Dmu(phydbl mu, int n_jumps, phydbl dt, phydbl a, phydbl b, phydbl lexp, int min_n, int jps_dens)
{
  if(n_jumps < 0) /* Marginal, i.e., the number of jumps is not fixed */
    {
      phydbl var,cumpoissprob,dens,mean,poissprob,ab2,gammadens,lexpdt;
      int n,up,down;
      
      var          = 0.0;
      cumpoissprob = 0.0;
      dens         = 0.0;
      n            = 0;
      mean         = a*b;
      ab2          = a*b*b;
      lexpdt       = lexp*dt;  
      
      RATES_Bracket_N_Jumps(&up,&down,lexpdt);
      For(n,MAX(down,min_n)-1) cumpoissprob += Dpois(n,lexpdt,NO);
      
      for(n=MAX(down,min_n);n<up+1;n++)
	{
	  poissprob    = Dpois(n,lexpdt,NO); /* probability of having n jumps */      
	  var          = (2./(n+2.))*ab2; /* var(mu|n) = var(mu|n=0) * 2 / (n+2) */
	  gammadens    = Dgamma_Moments(mu,mean,var);
	  dens         += poissprob * gammadens;
	  cumpoissprob += poissprob;
	  if(cumpoissprob > 1.-1.E-04) break;
	}
      
      if(dens < 1.E-70) dens = 1.E-70;

      return(dens);      
    }
  else /* Joint, i.e., return P(mu | dt, n_jumps) */
    {
      phydbl mean, var, density;


      mean = 1.0;
      var = (2./(n_jumps+2.))*a*b*b;

      if(jps_dens)
	density = Dgamma_Moments(mu,mean,var) * Dpois(n_jumps,dt*lexp,NO);
      else
	density = Dgamma_Moments(mu,mean,var);
      
      if(density < 1.E-70) density = 1.E-70;

      return density;
    }
}
  
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* Given the times of nodes a (ta) and d (td), the shape of the gamma distribution of instantaneous
   rates, the parameter of the exponential distribution of waiting times between rate jumps and the 
   instantaneous rate at t_node a, this function works out an expected number of (amino-acids or 
   nucleotide) substitutions per site.
*/
void RATES_Expect_Number_Subst(phydbl t_beg, phydbl t_end, phydbl r_beg,  int *n_jumps, phydbl *mean_r, phydbl *r_end, t_rate *rates, t_tree *tree)
{
  phydbl curr_r, curr_t, next_t;

  switch(rates->model_id)
    {
    case COMPOUND_COR:case COMPOUND_NOCOR:
      {
	/* Compound Poisson */
	if(rates->model_id == COMPOUND_COR)
	  {
	    curr_r  = r_beg;
	    *mean_r = r_beg;
	  }
	else
	  {
	    curr_r  = Rgamma(rates->nu,1./rates->nu);;
	    *mean_r = curr_r;
	  }

	curr_t = t_beg + Rexp(rates->lexp); /* Exponentially distributed waiting times */
	next_t = curr_t;
	
	*n_jumps = 0;
	while(curr_t < t_end)
	  {
	    curr_r = Rgamma(rates->nu,1./rates->nu); /* Gamma distributed random instantaneous rate */
	    
	    (*n_jumps)++;
	    
	    next_t = curr_t + Rexp(rates->lexp);
	    
	    if(next_t < t_end)
	      {
		*mean_r = (1./(next_t - t_beg)) * (*mean_r * (curr_t - t_beg) + curr_r * (next_t - curr_t));
	      }
	    else
	      {
		*mean_r = (1./(t_end - t_beg)) * (*mean_r * (curr_t - t_beg) + curr_r * (t_end - curr_t));
	      }
	    curr_t = next_t;
	  }
	
	/*   PhyML_Printf("\n. [%3d %f %f]",*n_jumps,*mean_r,r_beg); */
	
	if(*mean_r < rates->min_rate) *mean_r = rates->min_rate;
	if(*mean_r > rates->max_rate) *mean_r = rates->max_rate;

	*r_end = curr_r;
	break;
      }
    case EXPONENTIAL:
      {
	*mean_r = Rexp(rates->nu);

	if(*mean_r < rates->min_rate) *mean_r = rates->min_rate;
	if(*mean_r > rates->max_rate) *mean_r = rates->max_rate;

	*r_end  = *mean_r;
	break;
      }
    case LOGNORMAL:
      {
	*mean_r = Rgamma(rates->nu,1./rates->nu);

	if(*mean_r < rates->min_rate) *mean_r = rates->min_rate;
	if(*mean_r > rates->max_rate) *mean_r = rates->max_rate;

	*r_end  = *mean_r;
	break;
      }
    case THORNE:
      {
	phydbl sd,mean;
	int err;
	
	sd = SQRT(rates->nu*FABS(t_beg-t_end));
	mean = r_beg;

	*mean_r = Rnorm_Trunc(mean,sd,rates->min_rate,rates->max_rate,&err);

	if(err) PhyML_Printf("\n. %s %d %d",__FILE__,__LINE__,tree->mcmc->run);
	*r_end  = *mean_r;
	break;
      }
    default:
      {
	PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d\n",__FILE__,__LINE__);
	Exit("\n. Model not implemented yet.\n");
	break;
      }
    }
}
  
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Get_Mean_Rates_Pre(t_node *a, t_node *d, t_edge *b, phydbl r_a, t_tree *tree)
{
  phydbl a_t,d_t;
  phydbl mean_r;
  int n_jumps;
  phydbl r_d;

  a_t = tree->times->nd_t[a->num];
  d_t = tree->times->nd_t[d->num];

  mean_r = -1.;
  n_jumps = -1;
  r_d = -1.;
  RATES_Expect_Number_Subst(a_t,d_t,r_a,&n_jumps,&mean_r,&r_d,tree->rates,tree);
  
  tree->rates->br_r[d->num]   = mean_r;
  tree->rates->true_r[d->num] = mean_r;
  tree->times->t_jps[d->num]  = n_jumps;


  /* Move to the next branches */
  if(d->tax) return;
  else
    {
      int i;
      
      for(i=0;i<3;i++)
	{
	  if((d->v[i] != a) && (d->b[i] != tree->e_root))
	    {
	      RATES_Get_Mean_Rates_Pre(d,d->v[i],d->b[i],r_d,tree);
	    }
	}
    }

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Random_Branch_Lengths(t_tree *tree)
{
  phydbl r0;
   
  r0 = 1.0;

  tree->rates->br_r[tree->n_root->num] = r0;

  RATES_Get_Mean_Rates_Pre(tree->n_root,tree->n_root->v[2],NULL,r0,tree);
  RATES_Get_Mean_Rates_Pre(tree->n_root,tree->n_root->v[1],NULL,r0,tree);

  RATES_Check_Rates(tree);

  RATES_Update_Edge_Lengths(tree);
  RATES_Initialize_True_Rates(tree);

  tree->n_root_pos = 
    tree->rates->cur_l[tree->n_root->v[2]->num] /
    (tree->rates->cur_l[tree->n_root->v[2]->num] + tree->rates->cur_l[tree->n_root->v[1]->num]);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* Scale relative rates (on edges) so that they average to one */
void RATES_Update_Normalization_Factor(t_tree *tree)
{
  phydbl dt,rdt,T,RT;
  int i;

  rdt = 0.0;
  dt = 0.0;
  T = 0.0;
  RT = 0.0;
  for(i=0;i<2*tree->n_otu-2;++i)
    {
      assert(tree->a_nodes[i] != tree->n_root);
      dt = fabs(tree->times->nd_t[i] - tree->times->nd_t[tree->a_nodes[i]->anc->num]);
      rdt = dt*tree->rates->br_r[i];
      T+=dt;
      RT+=rdt;
    }
  assert(RT!=0.0);
  tree->rates->norm_fact = T/RT;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* Test whether rates are normalized */
void RATES_Check_Rates(t_tree *tree)
{
  phydbl scale,dt,T,eps;
  int i;

  eps = 1.E-6;
  
  scale = 0.0;
  dt = 0.0;
  T = 0.0;
  for(i=0;i<2*tree->n_otu-2;++i)
    {
      assert(tree->a_nodes[i] != tree->n_root);
      dt = fabs(tree->times->nd_t[i] - tree->times->nd_t[tree->a_nodes[i]->anc->num]);
      T+=dt;
      scale += tree->rates->br_r[i] * dt;
    }

  scale /= T;
  
  if(scale > 1.+eps || scale < 1.-eps)
    {
      PhyML_Fprintf(stderr,"\n. Relative rates are not normalised!...");
      assert(false);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Init_Triplets(t_tree *tree)
{
  int i;
  for(i=0;i<2*tree->n_otu-1;++i) tree->rates->triplet[i] = 0.0;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Initialize_True_Rates(t_tree *tree)
{
  int i;
  for(i=0;i<2*tree->n_otu-2;++i) tree->rates->true_r[i] = tree->rates->br_r[i];
}
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Get_Rates_From_Bl(t_tree *tree)
{
  phydbl dt,cr;
  t_node *left, *rght;
  int i;

  dt = -1.0;
  cr = tree->rates->clock_r;

  if(tree->n_root)
    {
      dt = FABS(tree->times->nd_t[tree->n_root->num] - tree->times->nd_t[tree->n_root->v[2]->num]);
      tree->rates->br_r[tree->n_root->v[2]->num] = 0.5 * tree->e_root->l->v / (dt*cr);
      dt = FABS(tree->times->nd_t[tree->n_root->num] - tree->times->nd_t[tree->n_root->v[1]->num]);
      tree->rates->br_r[tree->n_root->v[1]->num] = 0.5 * tree->e_root->l->v / (dt*cr);
    }
  

  for(i=0;i<2*tree->n_otu-3;++i)
    {
      if(tree->a_edges[i] != tree->e_root)
	{
	  left = tree->a_edges[i]->left;
	  rght = tree->a_edges[i]->rght;
	  dt = FABS(tree->times->nd_t[left->num] - tree->times->nd_t[rght->num]);	  
	  
	  if(left->anc == rght) tree->rates->br_r[left->num] = tree->a_edges[i]->l->v / (dt*cr);
	  else                  tree->rates->br_r[rght->num] = tree->a_edges[i]->l->v / (dt*cr);
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl RATES_Lk_Jumps(t_tree *tree)
{
  int i,n_jps;
  phydbl dens,dt,lexp;
  t_node *n;

  n = NULL;
  lexp = tree->rates->lexp;
  n_jps = 0;
  dt = 0.0;
  dens = 0.0;

  for(i=0;i<2*tree->n_otu-2;++i)
    {
      n = tree->a_nodes[i];
      dt = FABS(tree->times->nd_t[n->num]-tree->times->nd_t[n->anc->num]);
      n_jps = tree->times->n_jps[n->num];
      dens += Dpois(n_jps,lexp*dt,YES);
    }

  tree->times->c_lnL_jps = dens;
  
  return dens;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Update_Edge_Lengths(t_tree *tree)
{
  if(tree->is_mixt_tree == YES)
    {
      MIXT_RATES_Update_Edge_Lengths(tree);
      return;
    }
  
  RATES_Update_Normalization_Factor(tree);

  RATES_Update_Edge_Lengths_Pre(tree->n_root,tree->n_root->v[1],tree->n_root->b[1],tree);
  RATES_Update_Edge_Lengths_Pre(tree->n_root,tree->n_root->v[2],tree->n_root->b[2],tree);
  
  RATES_Update_One_Edge_Length(tree->e_root,tree);
  
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Update_Edge_Lengths_Pre(t_node *a, t_node *d, t_edge *b, t_tree *tree)
{
  RATES_Update_One_Edge_Length(b,tree);
  
  if(d->tax) return;
  else
    {
      int i;
      for(i=0;i<3;++i) 
	if((d->v[i] != a) && (d->b[i] != tree->e_root)) 
	  RATES_Update_Edge_Lengths_Pre(d,d->v[i],d->b[i],tree);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Update_One_Edge_Length(t_edge *b, t_tree *tree)
{

  if(b == tree->e_root || b == tree->n_root->b[1] || b == tree->n_root->b[2])
    {
      if(b == tree->e_root)
        {
          RATES_Update_One_Edge_Length_Core(tree->n_root->b[1],tree);
          RATES_Update_One_Edge_Length_Core(tree->n_root->b[2],tree);
        }
      else if(b == tree->n_root->b[1])
        {
          RATES_Update_One_Edge_Length_Core(tree->n_root->b[1],tree);
        }
      else if(b == tree->n_root->b[2])
        {
          RATES_Update_One_Edge_Length_Core(tree->n_root->b[2],tree);
        }
      
      if(tree->mod && tree->mod->log_l == YES)
        {
          tree->e_root->l->v = 
            exp(tree->n_root->b[1]->l->v) +
            exp(tree->n_root->b[2]->l->v) ;
          tree->e_root->l->v = log(tree->e_root->l->v);
        }
      else
        {
          tree->e_root->l->v = tree->n_root->b[1]->l->v + tree->n_root->b[2]->l->v;
        }
      
      tree->rates->u_cur_l[tree->e_root->num] = tree->e_root->l->v;
      tree->n_root_pos = tree->n_root->b[2]->l->v / tree->e_root->l->v;
  
      if(tree->rates->model_id == GUINDON)
        {
          phydbl t0,t1,t2;
          t_node *n0, *n1;
          
          n0 = tree->n_root->v[2];
          n1 = tree->n_root->v[1];
          t1 = tree->times->nd_t[tree->n_root->v[2]->num];
          t2 = tree->times->nd_t[tree->n_root->v[1]->num];
          t0 = tree->times->nd_t[tree->n_root->num];
          
          tree->e_root->l->v = 
            (t1-t0)/(t1+t2-2.*t0)*tree->rates->cur_gamma_prior_mean[n0->num] +
            (t2-t0)/(t1+t2-2.*t0)*tree->rates->cur_gamma_prior_mean[n1->num];
          
          tree->e_root->l_var->v = 
            pow((t1-t0)/(t1+t2-2.*t0),2)*tree->rates->cur_gamma_prior_var[n0->num] +
            pow((t2-t0)/(t1+t2-2.*t0),2)*tree->rates->cur_gamma_prior_var[n1->num];
        }
    }
  else
    {
      RATES_Update_One_Edge_Length_Core(b,tree);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Update_One_Edge_Length_Core(t_edge *b, t_tree *tree)
{      
  phydbl dt,rr,ra,rd,ta,td,nu,cr,Z;
  t_node *a, *d;

  if(b->left->anc == b->rght)
    {
      d = b->left;
      a = b->rght;
    }
  else
    {
      d = b->rght;
      a = b->left;
    }

  if(d->anc != a)
    {
      PhyML_Printf("\n. a->num: %d d->num: %d d->anc: %d b->left: %d b->rght: %d b->num: %d",
                   a->num,
                   d->num,
                   d->anc->num,
                   b->left->num,
                   b->rght->num,
                   b->num);
    }
  
  assert(a);
  assert(d);
  assert(d->anc == a);
  
  ra = rd = -1.;
  
  if(tree->rates->model_id == LOGNORMAL ||
     tree->rates->model_id == THORNE ||
     tree->rates->model_id == STRICTCLOCK)
    {
      rd = tree->rates->br_r[d->num];
      ra = tree->rates->br_r[a->num];
    }
  else if(tree->rates->model_id == GUINDON)
    {
      rd = tree->rates->nd_r[d->num];
      ra = tree->rates->nd_r[a->num];
    }
  else assert(FALSE);
  
  dt = fabs(tree->times->nd_t[d->num] - tree->times->nd_t[a->num]);
  cr = tree->rates->clock_r;
  td = tree->times->nd_t[d->num];
  ta = tree->times->nd_t[a->num];
  nu = tree->rates->nu;
  Z = tree->rates->norm_fact;
  rr = -1.0;
  
  if(tree->rates->model_id == LOGNORMAL)
    {
      tree->rates->cur_l[d->num] = dt*rd*cr*Z;
    }
  
  if(tree->rates->model_id == THORNE)
    {
      rr = (ra+rd)/2.;          
      tree->rates->cur_l[d->num] = dt*rr*cr*Z;
    }
  
  if(tree->rates->model_id == GUINDON)
    {
      phydbl m,v;
      
      Integrated_Geometric_Brownian_Bridge_Moments(dt,ra,rd,nu,&m,&v);
      
      if(isnan(m) || isnan(v) || m < 0.0 || v < 0.0)
        {
          PhyML_Fprintf(stderr,"\n. dt: %G ra: %G rd: %G nu: %G m: %G v: %G a is root ? %d d is root ? %d",
                        dt,ra,rd,nu,m,v,
                        a==tree->n_root,
                        d==tree->n_root);
        }
      
      m *= dt*cr; // the actual rate average is m * dt. We multiply by dt in order to derive the value for the branch length
      v *= (dt*dt)*(cr*cr);
      
      tree->rates->cur_gamma_prior_mean[d->num] = m;
      tree->rates->cur_gamma_prior_var[d->num]  = v;
      
      tree->rates->cur_l[d->num] = tree->rates->cur_gamma_prior_mean[d->num]; // Required for having proper branch lengths in Write_Tree function
    }
  
  if(tree->rates->model_id == STRICTCLOCK)
    {
      tree->rates->cur_l[d->num] = dt*cr;          
    }
  
  if(tree->mod && tree->mod->log_l == YES) tree->rates->cur_l[d->num] = log(tree->rates->cur_l[d->num]);
  
  if(b)
    {
      b->l->v                      = tree->rates->cur_l[d->num];
      tree->rates->u_cur_l[b->num] = tree->rates->cur_l[d->num];
      b->l_var->v                  = tree->rates->cur_gamma_prior_var[d->num];
    }
  
  if(b && (isnan(b->l->v) || isnan(b->l_var->v)))
    {
      PhyML_Fprintf(stderr,"\n. dt=%G rr=%G cr=%G ra=%G rd=%G nu=%G b->l_var=%f b->l=%f Z=%f",dt,rr,tree->rates->clock_r,ra,rd,nu,b->l_var->v,b->l->v,Z);	  
      PhyML_Fprintf(stderr,"\n. ta=%G td=%G ra*cr=%G rd*cr=%G sd=%G",ta,td,ra,rd,SQRT(dt*nu));
      /* assert(FALSE); */
    }
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Bl_To_Bl(t_tree *tree)
{
  RATES_Bl_To_Bl_Pre(tree->n_root,tree->n_root->v[2],NULL,tree);
  RATES_Bl_To_Bl_Pre(tree->n_root,tree->n_root->v[1],NULL,tree);
  /* tree->rates->cur_l[tree->n_root->v[2]->num] = tree->a_edges[tree->e_root->num]->l->v * tree->n_root_pos; */
  /* tree->rates->cur_l[tree->n_root->v[1]->num] = tree->a_edges[tree->e_root->num]->l->v * (1. - tree->n_root_pos); */
  tree->rates->cur_l[tree->n_root->v[2]->num] = tree->a_edges[tree->e_root->num]->l->v * 0.5;
  tree->rates->cur_l[tree->n_root->v[1]->num] = tree->a_edges[tree->e_root->num]->l->v * (1. - 0.5);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Bl_To_Bl_Pre(t_node *a, t_node *d, t_edge *b, t_tree *tree)
{

  if(b) 
    {
      tree->rates->cur_l[d->num] = b->l->v;
    }

  if(d->tax) return;
  else
    {
      int i;

      for(i=0;i<3;i++) 
	if((d->v[i] != a) && (d->b[i] != tree->e_root)) 
	  RATES_Bl_To_Bl_Pre(d,d->v[i],d->b[i],tree);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Bl_To_Ml(t_tree *tree)
{
  RATES_Bl_To_Ml_Pre(tree->n_root,tree->n_root->v[2],NULL,tree);
  RATES_Bl_To_Ml_Pre(tree->n_root,tree->n_root->v[1],NULL,tree);
  tree->rates->u_ml_l[tree->e_root->num] = tree->a_edges[tree->e_root->num]->l->v;
  tree->rates->ml_l[tree->n_root->v[2]->num] = tree->rates->u_ml_l[tree->e_root->num] * tree->n_root_pos;
  tree->rates->ml_l[tree->n_root->v[1]->num] = tree->rates->u_ml_l[tree->e_root->num] * (1. - tree->n_root_pos);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Bl_To_Ml_Pre(t_node *a, t_node *d, t_edge *b, t_tree *tree)
{

  if(b) 
    {
      tree->rates->u_ml_l[b->num] = b->l->v;
      tree->rates->ml_l[d->num]   = b->l->v;
    }

  if(d->tax) return;
  else
    {
      int i;

      for(i=0;i<3;i++) 
	if((d->v[i] != a) && (d->b[i] != tree->e_root)) 
	  RATES_Bl_To_Ml_Pre(d,d->v[i],d->b[i],tree);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void RATES_Get_Cov_Matrix_Rooted(phydbl *unroot_cov, t_tree *tree)
{
  int i,dim;

  dim = 2*tree->n_otu-3;

  RATES_Get_Cov_Matrix_Rooted_Pre(tree->n_root,tree->n_root->v[2],NULL,unroot_cov,tree);
  RATES_Get_Cov_Matrix_Rooted_Pre(tree->n_root,tree->n_root->v[1],NULL,unroot_cov,tree);

  for(i=0;i<dim+1;++i) tree->rates->cov_l[i*(dim+1)+tree->n_root->v[2]->num] /= 2.;
  for(i=0;i<dim+1;++i) tree->rates->cov_l[i*(dim+1)+tree->n_root->v[1]->num] /= 2.;
  for(i=0;i<dim+1;++i) tree->rates->cov_l[tree->n_root->v[2]->num*(dim+1)+i] /= 2.;
  for(i=0;i<dim+1;++i) tree->rates->cov_l[tree->n_root->v[1]->num*(dim+1)+i] /= 2.;

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Get_Cov_Matrix_Rooted_Pre(t_node *a, t_node *d, t_edge *b, phydbl *cov, t_tree *tree)
{
  int i, dim;
  t_node *n;

  dim = 2*tree->n_otu-3;
  n   = NULL;

  for(i=0;i<dim;i++) 
    { 
      if(tree->a_edges[i] != tree->e_root)
	{
	  n = 
	    (tree->a_edges[i]->left->anc == tree->a_edges[i]->rght)?
	    (tree->a_edges[i]->left):
	    (tree->a_edges[i]->rght);

	  if(b)
	    {
	      tree->rates->cov_l[d->num*(dim+1) + n->num] = cov[b->num*dim + i];
	    }
	  else
	    {
	      tree->rates->cov_l[d->num*(dim+1) + n->num] = cov[tree->e_root->num*dim + i];
	    }
	}
      else
	{
	  n = tree->e_root->left;
	  if(b)
	    tree->rates->cov_l[d->num*(dim+1) + n->num] = cov[b->num*dim + i];
	  else
	    tree->rates->cov_l[d->num*(dim+1) + n->num] = cov[tree->e_root->num*dim + i];

	  n = tree->e_root->rght;
	  if(b)
	    tree->rates->cov_l[d->num*(dim+1) + n->num] = cov[b->num*dim + i];
	  else
	    tree->rates->cov_l[d->num*(dim+1) + n->num] = cov[tree->e_root->num*dim + i];
	}
    }


  if(d->tax) return;
  else
    {
      for(i=0;i<3;i++)
	if((d->v[i] != a) && (d->b[i] != tree->e_root)) 
	  RATES_Get_Cov_Matrix_Rooted_Pre(d,d->v[i],d->b[i],cov,tree);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Covariance_Mu(t_tree *tree)
{
  int i,j;
  phydbl dt,var;
  int dim;
  int lca_num;

  dim = 2*tree->n_otu-2;

  for(i=0;i<dim*dim;++i) tree->rates->cov_r[i] = 0.0;

  dt =  tree->times->nd_t[tree->n_root->v[2]->num] - tree->times->nd_t[tree->n_root->num];
  var = dt * tree->rates->nu;
  tree->rates->cov_r[tree->n_root->v[2]->num*dim+tree->n_root->v[2]->num] = var;


  dt = tree->times->nd_t[tree->n_root->v[1]->num] - tree->times->nd_t[tree->n_root->num];
  var = dt * tree->rates->nu;
  tree->rates->cov_r[tree->n_root->v[1]->num*dim+tree->n_root->v[1]->num] = var;

  RATES_Variance_Mu_Pre(tree->n_root,tree->n_root->v[2],tree);
  RATES_Variance_Mu_Pre(tree->n_root,tree->n_root->v[1],tree);

  for(i=0;i<dim;i++)
    {
      for(j=i+1;j<dim;j++)
	{
	  lca_num = tree->rates->lca[i*(dim+1)+j]->num;
	  if(lca_num < dim) 
	    {
	      tree->rates->cov_r[i*dim+j] = tree->rates->cov_r[lca_num*dim+lca_num];	    
	      tree->rates->cov_r[j*dim+i] = tree->rates->cov_r[i*dim+j];
	    }
	  else if(lca_num == dim)
	    {
	      tree->rates->cov_r[i*dim+j] = 0.0;	    
	      tree->rates->cov_r[j*dim+i] = 0.0;
	    }
	  else
	    {
	      PhyML_Fprintf(stderr,"\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	      Exit("\n");
	    }
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Variance_Mu_Pre(t_node *a, t_node *d, t_tree *tree)
{
  int dim;
  phydbl var0;
  phydbl dt1,var1;
  phydbl dt2,var2;
  int i;
  int dir1, dir2;

  dim = 2*tree->n_otu-2;

  if(d->tax) return;

  dir1 = dir2 = -1;
  for(i=0;i<3;i++)
    {
      if((d->v[i] != a) && (d->b[i] != tree->e_root))
	{
	  if(dir1 < 0) dir1 = i;
	  else         dir2 = i;
	}
    }


  var0 = tree->rates->cov_r[d->num*dim+d->num];

  dt1  = tree->times->nd_t[d->v[dir1]->num] - tree->times->nd_t[d->num];
  var1 = tree->rates->nu*dt1;

  dt2  = tree->times->nd_t[d->v[dir2]->num] - tree->times->nd_t[d->num];
  var2 = tree->rates->nu*dt2;

  tree->rates->cov_r[d->v[dir1]->num*dim+d->v[dir1]->num] = var0+var1;
  tree->rates->cov_r[d->v[dir2]->num*dim+d->v[dir2]->num] = var0+var2;


  for(i=0;i<3;i++)
    {
      if((d->v[i] != a) && (d->b[i] != tree->e_root))
	{
	  RATES_Variance_Mu_Pre(d,d->v[i],tree);
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Fill_Lca_Table(t_tree *tree)
{
  int i,j;
  int dim;
  
  dim = 2*tree->n_otu-1;

  for(i=0;i<dim;i++)
    {
      for(j=i;j<dim;j++)
	{
	  tree->rates->lca[i*dim+j] = Find_Lca_Pair_Of_Nodes(tree->a_nodes[i],tree->a_nodes[j],tree);
	  tree->rates->lca[j*dim+i] = tree->rates->lca[i*dim+j];
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* Get V(L_{i} | L_{-i}) for all i */
void RATES_Get_Conditional_Variances(t_tree *tree)
{
  int i,j;
  short int *is_1;
  phydbl *a;
  int dim;
  t_edge *b;
  phydbl *cond_mu, *cond_cov;

  dim      = 2*tree->n_otu-3;
  a        = tree->rates->_2n_vect1;
  is_1     = tree->rates->_2n_vect5;
  b        = NULL;
  cond_mu  = tree->rates->_2n_vect2;
  cond_cov = tree->rates->_2n2n_vect1;

  for(i=0;i<dim;i++) a[i] = tree->rates->mean_l[i] * (Uni() * 0.2 + 0.9);

  for(i=0;i<dim;i++)
    {
      b = tree->a_edges[i];

      for(j=0;j<dim;j++) is_1[j] = 0;
      is_1[b->num]       = 1;

      For(j,dim*dim) cond_cov[j] = 0.0;
      for(j=0;j<dim;j++)     cond_mu[j]  = 0.0;
      Normal_Conditional(tree->rates->mean_l,tree->rates->cov_l,a,dim,is_1,1,cond_mu,cond_cov);
      
      tree->rates->cond_var[b->num] = cond_cov[b->num*dim+b->num];
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void RATES_Get_All_Reg_Coeff(t_tree *tree)
{
  int i,j;
  short int *is_1;
  phydbl *a;
  int dim;
  t_edge *b;

  dim  = 2*tree->n_otu-3;
  a    = tree->rates->_2n_vect1;
  is_1 = tree->rates->_2n_vect5;
  b    = NULL;

  for(i=0;i<dim;i++) a[i] = tree->rates->mean_l[i] * (Uni() * 0.2 + 0.9);

  for(i=0;i<dim;i++)
    {
      b = tree->a_edges[i];

      for(j=0;j<dim;j++) is_1[j] = 0;
      is_1[b->num]       = 1;
      
      Get_Reg_Coeff(tree->rates->mean_l,tree->rates->cov_l,a,dim,is_1,1,tree->rates->reg_coeff+b->num*dim);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


/* Get V(L_{i} | L_{-i}) for all i */
void RATES_Get_Trip_Conditional_Variances(t_tree *tree)
{
  int i,j;
  short int *is_1;
  phydbl *a;
  phydbl *cond_mu, *cond_cov;
  t_node *n;
  int n_otu;

  a        = tree->rates->_2n_vect1;
  is_1     = tree->rates->_2n_vect5;
  cond_mu  = tree->rates->_2n_vect2;
  cond_cov = tree->rates->_2n2n_vect1;
  n        = NULL;
  n_otu    = tree->n_otu;

  For(i,2*n_otu-3) a[i] = tree->rates->mean_l[i] * (Uni() * 0.2 + 0.9);

  for(i=0;i<2*n_otu-2;++i)
    {
      n = tree->a_nodes[i];
      if(!n->tax)
	{
	  For(j,2*n_otu-3) is_1[j] = 0;
	  is_1[n->b[0]->num] = 1;
	  is_1[n->b[1]->num] = 1;
	  is_1[n->b[2]->num] = 1;
	  
	  Normal_Conditional_Unsorted(tree->rates->mean_l,tree->rates->cov_l,a,2*n_otu-3,is_1,3,cond_mu,cond_cov);
	  
	  for(j=0;j<9;j++) tree->rates->trip_cond_cov[n->num*9+j] = cond_cov[j];
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void RATES_Get_All_Trip_Reg_Coeff(t_tree *tree)
{
  int i,j;
  short int *is_1;
  phydbl *a;
  t_node *n;
  int n_otu;

  a     = tree->rates->_2n_vect1;
  is_1  = tree->rates->_2n_vect5;
  n     = NULL;
  n_otu = tree->n_otu;

  for(i=0;i<2*n_otu-3;++i) a[i] = tree->rates->mean_l[i] * (Uni() * 0.2 + 0.9);

  for(i=0;i<2*n_otu-2;++i)
    {
      n = tree->a_nodes[i];
      if(!n->tax)
	{	  
	  For(j,2*n_otu-3) is_1[j] = 0;
	  is_1[n->b[0]->num] = 1;
	  is_1[n->b[1]->num] = 1;
	  is_1[n->b[2]->num] = 1;
	  
	  Get_Reg_Coeff(tree->rates->mean_l,tree->rates->cov_l,a,2*n_otu-3,is_1,3,tree->rates->trip_reg_coeff+n->num*(6*n_otu-9));
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Check_Lk_Rates(t_tree *tree, int *err)
{
  int i;
  phydbl u,u_anc;
  phydbl t,t_anc;
  
  *err = 0;

  For(i,2*tree->n_otu-2)
    {
      u     = tree->rates->br_r[i];
      u_anc = tree->rates->br_r[tree->a_nodes[i]->anc->num];
      t     = tree->times->nd_t[i];
      t_anc = tree->times->nd_t[tree->a_nodes[i]->anc->num];

      if(t_anc > t)
	{
	  PhyML_Printf("\n. %d %d u=%f u_anc=%f t=%f t_anc=%f",i,tree->a_nodes[i]->anc->num,u,u_anc,t,t_anc);
	  PhyML_Printf("\n. %d %d %d",tree->n_root->num,tree->n_root->v[2]->num,tree->n_root->v[1]->num);
	  *err = 1;
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl RATES_Realized_Substitution_Rate(t_tree *tree)
{
  RATES_Update_Edge_Lengths(tree);
  return(Tree_Length(tree)/TIMES_Tree_Length(tree));  
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl RATES_Expected_Tree_Length(t_tree *tree)
{
  int n;
  phydbl mean;

  n    = 0;
  mean = 0.0;
  RATES_Expected_Tree_Length_Pre(tree->n_root,tree->n_root->v[2],1.0,&mean,&n,tree);
  RATES_Expected_Tree_Length_Pre(tree->n_root,tree->n_root->v[1],1.0,&mean,&n,tree);

  if(n != 2*tree->n_otu-2)
    {
      PhyML_Fprintf(stderr,"\n. n=%d 2n-2=%d",n,2*tree->n_otu-2);
      PhyML_Fprintf(stderr,"\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");
    }

  return mean;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Expected_Tree_Length_Pre(t_node *a, t_node *d, phydbl eranc, phydbl *mean, int *n, t_tree *tree)
{
  phydbl erdes;
  int i;
  phydbl loc_mean;
  int loc_n;
  
  
/*  erdes = u_anc - */
/*     sd*(Dnorm((tree->rates->min_rate-u_anc)/sd,.0,1.) - Dnorm((tree->rates->max_rate-u_anc)/sd,.0,1.))/ */
/*        (Pnorm((tree->rates->max_rate-u_anc)/sd,.0,1.) - Pnorm((tree->rates->min_rate-u_anc)/sd,.0,1.)); */

  /* erdes = Norm_Trunc_Mean(eranc,sd,tree->rates->min_rate,tree->rates->max_rate); */
  erdes = 1.0;

  loc_mean = *mean;
  loc_n = *n;

  loc_mean *= loc_n;
  loc_mean += erdes;
  loc_mean /= (loc_n + 1);

  *mean = loc_mean;
  *n    = loc_n + 1;


  if(d->tax) return;
  else
    for(i=0;i<3;i++)
      if(d->v[i] != a && d->b[i] != tree->e_root)
	RATES_Expected_Tree_Length_Pre(d,d->v[i],erdes,mean,n,tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl Sample_Average_Rate(t_node *a, t_node *d, t_tree *tree)
{
  return(-1.);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Update_Mean_Br_Len(int iter, t_tree *tree)
{
  int i,dim;
  phydbl *mean;

  if(tree->rates->update_mean_l == NO) return;

  dim = 2*tree->n_otu-3;
  mean = tree->rates->mean_l;

  for(i=0;i<dim;i++)
    {      
     mean[i] *= (phydbl)iter;
     mean[i] += tree->a_edges[i]->l->v;
     mean[i] /= (phydbl)(iter+1);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void RATES_Update_Cov_Br_Len(int iter, t_tree *tree)
{
  int i,j,dim;
  phydbl *mean,*cov;

  if(tree->rates->update_cov_l == NO) return;

  dim = 2*tree->n_otu-3;
  mean = tree->rates->mean_l;
  cov  = tree->rates->cov_l;

  for(i=0;i<dim;i++)
    {     
      for(j=0;j<dim;j++)
	{
	  cov[i*dim+j] += mean[i]*mean[j];
	  cov[i*dim+j] *= (phydbl)tree->mcmc->run;
	  cov[i*dim+j] += tree->a_edges[i]->l->v * tree->a_edges[j]->l->v;
	  cov[i*dim+j] /= (phydbl)(tree->mcmc->run+1);
	  cov[i*dim+j] -= mean[i]*mean[j];

	  if(i == j && cov[i*dim+j] < MIN_VAR_BL) cov[i*dim+j] = MIN_VAR_BL;
	}      
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Set_Mean_L(t_tree *tree)
{
  int i;
  for(i=0;i<2*tree->n_otu-3;++i) 
    {
      tree->rates->mean_l[i] = tree->a_edges[i]->l->v;
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Record_Times(t_tree *mixt_tree)
{
  int i;
  t_tree *tree;

  tree = mixt_tree;
  do
    {
      if(tree->times->nd_t_recorded == YES)
        {
          PhyML_Fprintf(stderr,"\n. Overwriting recorded times is forbidden.\n");
          PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d\n",__FILE__,__LINE__);
          Exit("\n");
        }

      For(i,2*tree->n_otu-1) tree->times->buff_t[i] = tree->times->nd_t[i];
      tree = tree->next;
    }
  while(tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Reset_Times(t_tree *mixt_tree)
{
  int i;
  t_tree *tree;

  tree = mixt_tree;
  do
    {
      tree->times->nd_t_recorded = NO;
      for(i=0;i<2*tree->n_otu-1;++i) tree->times->nd_t[i] = tree->times->buff_t[i];
      tree = tree->next;
    }
  while(tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Record_Rates(t_tree *tree)
{
  int i;

  if(tree->rates->br_r_recorded == YES)
    {
      PhyML_Fprintf(stderr,"\n. Overwriting recorded rates is forbidden.\n");
      PhyML_Fprintf(stderr,"\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");
    }

  for(i=0;i<2*tree->n_otu-2;++i) tree->rates->buff_br_r[i] = tree->rates->br_r[i];
  for(i=0;i<2*tree->n_otu-1;++i) tree->rates->buff_nd_r[i] = tree->rates->nd_r[i];
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Reset_Rates(t_tree *tree)
{
  int i;
  tree->rates->br_r_recorded = NO;
  for(i=0;i<2*tree->n_otu-2;++i) tree->rates->br_r[i] = tree->rates->buff_br_r[i];
  for(i=0;i<2*tree->n_otu-1;++i) tree->rates->nd_r[i] = tree->rates->buff_nd_r[i];
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Set_Clock_And_Nu_Max(t_tree *tree)
{
  phydbl dt,nu;
  phydbl min_t;
  int i;
  phydbl step;
  phydbl l_max;
  phydbl max_clock;
  phydbl r_max;
  phydbl tune;
  phydbl pa,pb;

  if(tree->rates->model_id == THORNE || 
     tree->rates->model_id == GUINDON)
    {
      tune  = 1.05;
      r_max = tree->rates->max_rate;     
      l_max = tree->mod->l_max;
      
      min_t = .0;
      for(i=0;i<2*tree->n_otu-1;++i) if(tree->times->t_prior_min[i] < min_t) min_t = tree->times->t_prior_min[i];
      
      dt = FABS(min_t);
      max_clock = l_max / dt; 
      
      nu   = 1.E-10;
      step = 1.E-1;
      do
	{
	  do
	    {
	      nu += step;
	      pa = Dnorm(0.0,  0.0,SQRT(nu*dt)); 
	      pb = Dnorm(r_max,0.0,SQRT(nu*dt));
	    }while(pa/pb > tune);
	  nu -= step;
	  step /= 10.;
	}while(step > 1.E-10);
      
      tree->rates->max_nu    = nu;
      /* tree->rates->max_nu    = 1.0; */
      tree->rates->max_clock = max_clock;
      
      PhyML_Printf("\n. Clock rate parameter upper bound set to %f expected subst./site/time unit",tree->rates->max_clock);
      PhyML_Printf("\n. Autocorrelation parameter upper bound set to %f",tree->rates->max_nu);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Set_Birth_Rate_Boundaries(t_tree *tree)
{
  phydbl lbda;
  phydbl p_above_min,p_below_max;
  phydbl min,max;
  int assign = YES;

  min = -0.01*tree->times->t_prior_max[tree->n_root->num];
  max = -100.*tree->times->t_prior_min[tree->n_root->num];

  for(lbda = 0.0001; lbda < 10; lbda+=0.0001)
    {
      p_above_min = 1. - POW(1.-exp(-lbda*min),tree->n_otu);
      p_below_max = POW(1.-exp(-lbda*max),tree->n_otu);
 
      if(p_above_min < 1.E-10) 
	{ 
	  tree->times->birth_rate_max = lbda;
	  break;
	}
      if(p_below_max > 1.E-10 && assign==YES)
	{
	  assign = NO;
	  tree->times->birth_rate_min = lbda;
	}
    }
  
  /* tree->times->birth_rate_min = 1.E-6; */
  /* tree->times->birth_rate_max = 1.; */
  PhyML_Printf("\n. Birth rate lower bound set to %f.",tree->times->birth_rate_min);
  PhyML_Printf("\n. Birth rate upper bound set to %f.",tree->times->birth_rate_max);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl RATES_Get_Mean_Rate_In_Subtree(t_node *root, t_tree *tree)
{
  phydbl sum;
  int n;

  sum = 0.0;
  n   = 0;

  if(root->tax == NO)
    {
      if(root == tree->n_root)
	{
	  RATES_Get_Mean_Rate_In_Subtree_Pre(root,root->v[2],&sum,&n,tree);
	  RATES_Get_Mean_Rate_In_Subtree_Pre(root,root->v[1],&sum,&n,tree);
	}
      else
	{
	  int i;
	  for(i=0;i<3;i++)
	    {
	      if(root->v[i] != root->anc && root->b[i] != tree->e_root)
		{
		  RATES_Get_Mean_Rate_In_Subtree_Pre(root,root->v[i],&sum,&n,tree);
		}
	    }
	}
      return sum/(phydbl)n;
    }
  else
    {
      return 0.0;
    }  
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Get_Mean_Rate_In_Subtree_Pre(t_node *a, t_node *d, phydbl *sum, int *n, t_tree *tree)
{
  /* (*sum) += exp(tree->rates->nd_r[d->num]); */
  
  if(tree->rates->model_id == LOGNORMAL ||
     tree->rates->model_id == THORNE ||
     tree->rates->model_id == STRICTCLOCK)
    {
      (*sum) += tree->rates->br_r[d->num];
    }
  else if(tree->rates->model_id == GUINDON)
    {
      (*sum) += tree->rates->nd_r[d->num];
    }

  else assert(FALSE);

  (*n) += 1;

  if(d->tax == YES)  return;
  else
    {
      int i;
      for(i=0;i<3;++i)
	{
	  if(d->v[i] != a && d->b[i] != tree->e_root)
	    {
	      RATES_Get_Mean_Rate_In_Subtree_Pre(d,d->v[i],sum,n,tree);
	    }
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

char *RATES_Get_Model_Name(int model)
{
  char *s;

  s = (char *)mCalloc(T_MAX_NAME,sizeof(char));

  switch(model)
    {
    case GUINDON     : {strcpy(s,"integrated"); break;}
    case THORNE      : {strcpy(s,"autocorrelated"); break;}
    case LOGNORMAL   : {strcpy(s,"uncorrelated"); break;}
    case STRICTCLOCK : {strcpy(s,"strict clock"); break;}
    default : 
      {
	PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d\n",__FILE__,__LINE__);
	Exit("\n");
	break;
     }
    }
  return s;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int RATES_Check_Edge_Length_Consistency(t_tree *mixt_tree)
{
  int i;
  t_tree *tree;

  tree = mixt_tree;
  
  do
    {
      if(tree->is_mixt_tree == YES) tree = tree->next;
      
      for(i=0;i<2*tree->n_otu-3;++i)
        {
          if(tree->a_edges[i]->left->anc == tree->a_edges[i]->rght)
            {
              if(Are_Equal(tree->rates->cur_l[tree->a_edges[i]->left->num],
                           tree->a_edges[i]->l->v,
                           1.E-5) == NO)
                {
                  PhyML_Fprintf(stderr,"\n. cur_l: %G l: %G is_root: %d",
                                tree->rates->cur_l[tree->a_edges[i]->left->num],
                                tree->a_edges[i]->l->v,
                                tree->a_edges[i] == tree->e_root);
                  return 0;
                }
            }
          
          if(tree->a_edges[i]->rght->anc == tree->a_edges[i]->left)
            {
              if(Are_Equal(tree->rates->cur_l[tree->a_edges[i]->rght->num],
                           tree->a_edges[i]->l->v,
                           1.E-5) == NO)
                {
                  PhyML_Fprintf(stderr,"\n. cur_l: %G l: %G is_root: %d",
                                tree->rates->cur_l[tree->a_edges[i]->rght->num],
                                tree->a_edges[i]->l->v,
                                tree->a_edges[i] == tree->e_root);
                  return 0;
                }
            }
        }

      tree = tree->next;
    }
  while(tree);

  return 1;
  
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////



