/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include "iou.h"

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

short int IOU_Is_Iou(t_phyrex_mod *mod)
{
  if(mod->model_id == IOU) return(YES);
  return(NO);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
/* Mean velocity at the end of an edge, given the velocity at the start of it */
phydbl IOU_Velocity_Mean_Along_Edge(t_node *d, short int dim, t_tree *tree)
{
  phydbl y0,t,ta,mu;
  
  assert(d != tree->n_root);

  y0 = d->anc->ldsk->veloc->deriv[dim];
  t  = fabs(tree->times->nd_t[d->num] - tree->times->nd_t[d->anc->num]);
  ta = tree->mmod->ou_theta;
  mu = tree->mmod->ou_mu[dim];
  
  return(y0*exp(-ta*t)+mu*(1.-exp(-ta*t)));
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
/* Variance of velocity at the end of an edge, given the velocity at the start of it */
phydbl IOU_Velocity_Variance_Along_Edge(t_node *d, short int dim, t_tree *tree)
{
  phydbl var,t,sg2,ta;
  
  assert(d != tree->n_root);

  /* t = fabs(tree->times->nd_t[d->num] - tree->times->nd_t[d->anc->num]); */
  /* sg = tree->mmod->sigsq[dim]; */
  /* ta = tree->mmod->ou_theta; */
  
  /* var = sg/(2.*ta)*(1.-exp(-2.*ta*t)); */

  t = fabs(tree->times->nd_t[d->num] - tree->times->nd_t[d->anc->num]);
  sg2 = tree->mmod->sigsq[dim] * tree->mmod->sigsq_scale[d->num] * tree->mmod->sigsq_scale_norm_fact;
  ta = tree->mmod->ou_theta;
  
  var = sg2 /(2.*ta) * (1. - exp(-2.*ta*t));
  
  if(var < 0.0) var = 0.0;

  if(isinf(var) || isnan(var))
    {
      PhyML_Printf("\n. ta: %f t: %f sg2: %f sinh: %f",ta,t,sg2,sinh(ta*t));
    }
  
  assert(isnan(var) == NO);
  assert(isinf(var) == NO);
  
  return(var);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
/* Mean of location at the end of the edge, given velocities
   at both extremities */
phydbl IOU_Location_Mean_Along_Edge(t_node *d, short int dim, t_tree *tree)
{
  phydbl t,x0,y0,yt,mean,ta,mu;
  
  
  assert(d != tree->n_root);
  
  t = fabs(tree->times->nd_t[d->num] - tree->times->nd_t[d->anc->num]);

  if(t < SMALL) return(d->anc->ldsk->coord->lonlat[dim]);

  x0 = d->anc->ldsk->coord->lonlat[dim];
  y0 = d->anc->ldsk->veloc->deriv[dim];
  yt = d->ldsk->veloc->deriv[dim];

  ta = tree->mmod->ou_theta;
  mu = tree->mmod->ou_mu[dim];
  
  mean =
    /* x0 + (cosh(ta*t)-1.)/(sinh(ta*t)*ta)*(yt+y0) + */
    /* mu*t - (mu/ta)*(1.-exp(-ta*t))*(1.+(cosh(ta*t)-1.)/sinh(ta*t)); */
    x0 + (1./tanh(ta*t) - 1./sinh(ta*t))*(yt+y0)/ta +
    mu*t - (mu/ta)*(1.-exp(-ta*t))*(1.+1./tanh(ta*t)-1./sinh(ta*t));

  /* PhyML_Printf("\n. t: %f x0: %f y0: %f yt: %f ta: %f mu: %f cosh: %f sinh: %f", */
  /*              t,x0,y0,yt,ta,mu,cosh(ta*t),sinh(ta*t)); */
  
  assert(isnan(mean) == NO);
  
  return(mean);
    
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
/* Variance of location at the end of the edge, given velocities
   at both extremities */
phydbl IOU_Location_Variance_Along_Edge(t_node *d, short int dim, t_tree *tree)
{
  phydbl var,t,ta,sg;
  
  assert(d != tree->n_root);

  
  t  = fabs(tree->times->nd_t[d->num] - tree->times->nd_t[d->anc->num]);
  ta = tree->mmod->ou_theta;
  sg = tree->mmod->sigsq[dim];
  
  if(t < SMALL) return(0.0);

  /* var = sg/pow(ta,3)*(ta*t-2.*(cosh(ta*t)-1.)/sinh(ta*t)); */
  var = sg/pow(ta,3)*(ta*t-2.*(1./tanh(ta*t) - 1./sinh(ta*t)));

  if(var < 0.0) var = 0.0;
  
  assert(isnan(var) == NO);
  
  return(var);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl IOU_Prior(t_tree *tree)
{
  tree->mmod->c_lnP =
    IOU_Prior_Theta(tree) +
    IOU_Prior_Mu(tree) +
    RW_Prior_Sigsq(tree);

  return(tree->mmod->c_lnP);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl IOU_Prior_Theta(t_tree *tree)
{
  phydbl lbda,lnP;

  lbda = 1.0;
  lnP = log(lbda)-lbda*tree->mmod->ou_theta;

  return(lnP);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl IOU_Prior_Mu(t_tree *tree)
{
  return(0.0);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void IOU_Integrated_Location_Down(phydbl son_a, phydbl son_b, phydbl son_mu_down, phydbl son_var_down, phydbl son_var,
                                  phydbl bro_a, phydbl bro_b, phydbl bro_mu_down, phydbl bro_var_down, phydbl bro_var,
                                  phydbl son_logrem, phydbl bro_logrem,
                                  phydbl *mean, phydbl *var, phydbl *logrem)
{
  RW_Integrated_Lk_Down(son_a, son_b, son_mu_down, son_var_down, son_var,
                        bro_a, bro_b, bro_mu_down, bro_var_down, bro_var,
                        son_logrem, bro_logrem,
                        mean, var, logrem);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void IOU_Integrated_Location_Up(phydbl dad_mu_up, phydbl dad_var_up, phydbl dad_logrem_up,
                                phydbl son_a, phydbl son_b, phydbl son_var,
                                phydbl bro_a, phydbl bro_b, phydbl bro_mu_down, phydbl bro_var_down, phydbl bro_var, phydbl bro_logrem_down,
                                phydbl *mean, phydbl *var, phydbl *logrem)
{
  RW_Integrated_Lk_Up(dad_mu_up, dad_var_up, dad_logrem_up,
                      son_a, son_b, son_var,
                      bro_a, bro_b, bro_mu_down, bro_var_down, bro_var, bro_logrem_down,
                      mean, var, logrem);
}
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
