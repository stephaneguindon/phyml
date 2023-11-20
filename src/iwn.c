/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include "iwn.h"

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

short int IWN_Is_Iwn(t_phyrex_mod *mod)
{
  if(mod->model_id == IWNc  ||
     mod->model_id == RIWNc ||
     mod->model_id == IWNu  ||
     mod->model_id == RIWNu) return(YES);
  return(NO);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl IWN_Velocity_Mean_Along_Edge(t_node *d, short int dim, t_tree *tree)
{
  assert(d != tree->n_root);
  return(d->anc->ldsk->veloc->deriv[dim]);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl IWN_Velocity_Variance_Along_Edge(t_node *d, short int dim, t_tree *tree)
{
  phydbl logvar,dt;
  
  assert(d != tree->n_root);

  dt = fabs(tree->times->nd_t[d->num] - tree->times->nd_t[d->anc->num]);

  /* PhyML_Printf("\n. omega: %12f dt: %12f omega > dt: %3d", */
  /*              tree->mmod->omega, */
  /*              dt, */
  /*              !(tree->mmod->omega < dt)); */
  
  if(!(tree->mmod->omega < dt))
    {
      return(0.0);
    }
  else
    {
      logvar =
        log(tree->mmod->sigsq[dim]) +
        log(tree->mmod->sigsq_scale[d->num]) +
        log(tree->mmod->sigsq_scale_norm_fact);
      
      return(exp(logvar));
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl IWN_Location_Mean_Along_Edge(t_node *d, short int dim, t_tree *tree)
{
  phydbl dt, loc_beg, veloc_beg, veloc_end, epst;
  int gt;
  
  assert(d != tree->n_root);
  
  dt = fabs(tree->times->nd_t[d->num] - tree->times->nd_t[d->anc->num]);

  loc_beg = d->anc->ldsk->coord->lonlat[dim];
  veloc_beg = d->anc->ldsk->veloc->deriv[dim];
  veloc_end = d->ldsk->veloc->deriv[dim];

  if(!(tree->mmod->omega < dt))
    {
      return(loc_beg + veloc_beg*dt);
    }
  else
    {
      gt = IWN_Get_Gt(dt,tree->mmod->omega);
      epst = dt-gt*tree->mmod->omega;
      return(loc_beg + veloc_beg*tree->mmod->omega + .5*dt*(veloc_beg + veloc_end) + veloc_end*epst);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
/* Variance of location at the end of the edge, given velocities
   at both extremities */
phydbl IWN_Location_Variance_Along_Edge(t_node *d, short int dim, t_tree *tree)
{
  phydbl logvar,dt;
  int gt;

  assert(d != tree->n_root);
  
  dt = fabs(tree->times->nd_t[d->num] - tree->times->nd_t[d->anc->num]);

  if(!(tree->mmod->omega < dt))
    {
      return(0.0);
    }
  else
    {
      gt = IWN_Get_Gt(dt,tree->mmod->omega);
      
      logvar =
        log(tree->mmod->sigsq[dim]) + 
        log(tree->mmod->sigsq_scale[d->num]) +
        log(tree->mmod->sigsq_scale_norm_fact) +
        2.*log(tree->mmod->omega) +
        log(gt-1.)+
        log((gt-1.)/4. + 1.);
      
      return(exp(logvar));
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int IWN_Get_Gt(phydbl t, phydbl omega)
{
  /* int gt; */
  
  /* gt = 0; */
  /* while(!(gt*omega > t)) { gt++; } */

  /* return(gt-1); */

  return((int)(t/omega));
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl IWN_Prior(t_tree *tree)
{
  tree->mmod->c_lnP =
    IWN_Prior_Omega(tree) +
    RW_Prior_Sigsq(tree);

  return(tree->mmod->c_lnP);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl IWN_Prior_Omega(t_tree *tree)
{
  phydbl lbda,lnP;

  lbda = 1.0;
  lnP = log(lbda)-lbda*tree->mmod->omega;

  return(lnP);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void IWN_Integrated_Location_Down(phydbl dt1, phydbl dt2,
                                  phydbl av1, phydbl bv1, phydbl v1mu, phydbl v1var, phydbl dv1var,
                                  phydbl av2, phydbl bv2, phydbl v2mu, phydbl v2var, phydbl dv2var,
                                  phydbl v1logrem, phydbl v2logrem,
                                  phydbl veloc_d, phydbl veloc_v1, phydbl veloc_v2,
                                  phydbl omega,
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
      assert(dt2 < 2.*omega);
      epst = dt2 - omega;
      if(epst > .0) m = v2mu + omega * veloc_v2 + epst * veloc_d;
      else m = v2mu + omega * veloc_v2;
    }
  else if(dv2var + v2var > SMALL) // Null variance along d - v1
    {
      assert(dt1 < 2.*omega);
      epst = dt1 - omega;
      if(epst > .0) m = v1mu + omega * veloc_v1 + epst * veloc_d;
      else m = v1mu + omega * veloc_v1;
    }
  else // Null variance along d - v1 and d - v2
    {
      phydbl m1;
      
      epst = dt1 - omega;
      if(epst > .0) m1 = v1mu + omega * veloc_v1 + epst * veloc_d;
      else m1 = v1mu + omega * veloc_v1;

      /* PhyML_Printf("\n. dt1: %12f dt2: %12f v1mu: %12f v2mu: %12f veloc_v1: %12f veloc_v2: %12f veloc_d: %12f epst1: %12f epst2: %12f", */
      /*              dt1,dt2, */
      /*              v1mu,v2mu, */
      /*              veloc_v1,veloc_v2,veloc_d, */
      /*              dt1 - omega, */
      /*              dt2 - omega); */
                   
      /* PhyML_Printf("\n. m1=%f m2=%f",m1,m2); */

      /* assert(fabs(m1 - m2) < 1.E-10); */

      m = m1;
    }

  *mean = m;
  *var  = v;
  *logrem = logr;
  
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void IWN_Integrated_Location_Up(phydbl dt1, phydbl dt2,
                                phydbl av1, phydbl bv1, phydbl v1mu, phydbl v1var, phydbl av1var,
                                phydbl av2, phydbl bv2, phydbl v2mu, phydbl v2var, phydbl av2var,
                                phydbl v1logrem, phydbl v2logrem,
                                phydbl veloc_a, phydbl veloc_v1, phydbl veloc_v2,
                                phydbl omega,
                                phydbl *mean, phydbl *var, phydbl *logrem,
                                short int a_is_root)
{
  phydbl m,v,logr,epst;
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
          assert(dt2 < 2.*omega);
          epst = dt2 - omega;
          if(epst > .0) m = v2mu + omega * veloc_v2 + epst * veloc_a;
          else m = v2mu + omega * veloc_v2;
        }
      else if(av2var + v2var > SMALL) // Null variance along d - v1
        {
          assert(dt1 < 2.*omega);
          epst = dt1 - omega;
          if(epst > .0) m = v1mu + omega * veloc_a + epst * veloc_v1;
          else m = v1mu + omega * veloc_v1;
        }
      else
        {
          phydbl m1;
          
          epst = dt1 - omega;
          if(epst > .0) m1 = v1mu + omega * veloc_a + epst * veloc_v1;
          else m1 = v1mu + omega * veloc_v1;
                    
          m = m1;
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
          assert(dt2 < 2.*omega);
          epst = dt2 - omega;
          if(epst > .0) m = v2mu + omega * veloc_v2 + epst * veloc_a;
          else m = v2mu + omega * veloc_v2;
        }
    }
  
  *mean = m;
  *var  = v;
  *logrem = logr;

}
  
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
