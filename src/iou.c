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
  phydbl var,t,sg,ta;
  
  assert(d != tree->n_root);

  t = fabs(tree->times->nd_t[d->num] - tree->times->nd_t[d->anc->num]);
  sg = tree->mmod->sigsq[dim];
  ta = tree->mmod->ou_theta;
  
  var = sg/(2.*ta)*(1.-exp(-2.*ta*t));
     
  if(var < 0.0) var = 0.0;

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
    x0 + (cosh(ta*t)-1.)/sinh(ta*t)/ta*(yt+y0) +
    mu*t - (mu/ta)*(1.-exp(-ta*t))*(1.+(cosh(ta*t)-1.)/sinh(ta*t));
  
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
  
  var = sg/pow(ta,3)*(ta*t-2.*(cosh(ta*t)-1.)/sinh(ta*t));

  if(var < 0.0) var = 0.0;

  return(var);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl IOU_Prior(t_tree *tree)
{
  return(IOU_Prior_Theta(tree) +
         IOU_Prior_Mu(tree) +
         PHYREX_LnPrior_Sigsq(tree));
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
  
  *mean = m;
  *var  = v;
  *logrem = logr;

}
  
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
