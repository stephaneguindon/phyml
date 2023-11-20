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

void IOU_Integrated_Location_Down(phydbl dt1, phydbl dt2,
                                  phydbl av1, phydbl bv1, phydbl v1mu, phydbl v1var, phydbl dv1var,
                                  phydbl av2, phydbl bv2, phydbl v2mu, phydbl v2var, phydbl dv2var,
                                  phydbl v1logrem, phydbl v2logrem,
                                  phydbl *mean, phydbl *var, phydbl *logrem)
{

  phydbl m,v,logr,epst;
  int err;

  err = 0;
  m = epst = -1.;
  v = 0.0;
  logr = 0.0;
  
  if(dv1var + v1var > SMALL && dv2var + v2var > SMALL) // Standard case
    {
      v = pow(av1,2)/(v1var + dv1var) + pow(av2,2)/(v2var + dv2var);
      v = 1/v;

      m = (av1*(v1mu-bv1)/(v1var + dv1var) + av2*(v2mu-bv2)/(v2var + dv2var)) * v;

      logr  = v1logrem + v2logrem;
      logr -= log(fabs(av2*av1));
      logr += Log_Dnorm((v1mu-bv1)/av1,(v2mu-bv2)/av2,sqrt((v1var+dv1var)/pow(av1,2)+(v2var+dv2var)/pow(av2,2)),&err);
    }
  else if(dv1var + v1var > SMALL) // Null variance along d - v2
    {
      m = (v2mu-bv2)/av2;
    }
  else if(dv2var + v2var > SMALL) // Null variance along d - v1
    {
      m = (v1mu-bv1)/av1;
    }
  else // Null variance along d - v1 and d - v2
    {
      m = (v1mu-bv1)/av1;
    }

  *mean = m;
  *var  = v;
  *logrem = logr;
  
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void IOU_Integrated_Location_Up(phydbl dt1, phydbl dt2,
                                phydbl av1, phydbl bv1, phydbl v1mu, phydbl v1var, phydbl av1var,
                                phydbl av2, phydbl bv2, phydbl v2mu, phydbl v2var, phydbl av2var,
                                phydbl v1logrem, phydbl v2logrem,
                                phydbl *mean, phydbl *var, phydbl *logrem,
                                short int a_is_root)
{
  phydbl m,v,logr;
  int err;
  
  v = logr = 0.0;
  err = 0;

  
  if(a_is_root == NO)
    {
      if(pow(av1,2)*v1var+av1var > SMALL && av2var + v2var > SMALL) // Standard case
        {
          v     = pow(av2,2)/(v2var + av2var) + 1./(pow(av1,2)*v1var+av1var);
          v     = 1./v;
          
          m     = (av2*(v2mu-bv2)/(v2var + av2var) + (av1*v1mu+bv1)/(pow(av1,2)*v1var+av1var)) * v;          
          
          logr  = v1logrem + v2logrem;
          logr -= log(fabs(av2));
          logr += Log_Dnorm((v2mu-bv2)/av2,av1*v1mu+bv1,sqrt((v2var+av2var)/pow(av2,2)+pow(av1,2)*v1var+av1var),&err);
        }
      else if(pow(av1,2)*v1var+av1var > SMALL) // Null variance along d - v2
        {
          m = (v2mu-bv2)/av2;
        }
      else if(av2var + v2var > SMALL) // Null variance along d - v1
        {
          m = (v1mu-bv1)/av1;
        }
      else
        {
          m = (v1mu-bv1)/av1;
        }
    }
  else
    {
      if(v2var + av2var > SMALL)
        {
          m    = (v2mu-bv2)/av2;
          v    = (v2var + av2var)/pow(av2,2);
          logr = v2logrem;
          logr -= log(fabs(av2));
        }
      else
        {
          m    = (v2mu-bv2)/av2;
          v    = 0.0;
          logr = 0.0;
        }
    }


  if(isnan(m))
    {
      PhyML_Printf("\n. a_is_root ? %d",a_is_root);
      PhyML_Printf("\n. dt1: %f dt2: %f av1: %f av2: %f bv1: %f bv2: %f v1mu: %f v1var: %f av1var: %f v2mu: %f v2var: %f av2var: %f",
                   dt1, dt2,
                   av1, bv1, v1mu, v1var, av1var,
                   av2, bv2, v2mu, v2var, av2var);
                                
    }
  
  assert((isinf(m) || isinf(v)) == NO);
  assert((isnan(m) || isnan(v)) == NO);

  *mean = m;
  *var  = v;
  *logrem = logr;

}
  
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
