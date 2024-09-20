/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include "rw.h"
#include "rrw.h"

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl RW_Lk(t_tree *tree)
{
  RRW_Lk(tree);
  return(tree->mmod->c_lnL);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl RW_Prior(t_tree *tree)
{
  tree->mmod->c_lnP = RW_Prior_Sigsq(tree);
  tree->mmod->c_lnP += RW_Prior_Observational_Model(tree);
  return(tree->mmod->c_lnP);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl RW_Prior_Observational_Model(t_tree *tree)
{
    phydbl lnP;
    int i;

    lnP = 0.0;

    if (tree->contmod->obs_model == NO)
        return (0.0);

    for (i = 0; i < tree->mmod->n_dim; ++i)
    {
        lnP += Log_Dexp(tree->contmod->obs_var[i],
                        1. / tree->contmod->obs_var_prior_mean);
    }
    return (lnP);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl RW_Prior_Sigsq(t_tree *tree)
{
  phydbl lnP;
  int i,err;
  
  lnP = 0.0;
    
  if(tree->mmod->rw_prior_distrib == EXPONENTIAL_PRIOR)
    {
      for(i=0;i<tree->mmod->n_dim;++i)
        {
          lnP += Log_Dexp(tree->mmod->sigsq[i],
                          1./tree->mmod->rw_prior_mean);
        }
    }
  else if(tree->mmod->rw_prior_distrib == NORMAL_PRIOR)
    {
      for(i=0;i<tree->mmod->n_dim;++i)
        {
          lnP += Log_Dnorm(tree->mmod->sigsq[i],
                           tree->mmod->rw_prior_mean,
                           tree->mmod->rw_prior_sd,
                           &err);
        }
    }
  else if(tree->mmod->rw_prior_distrib == FLAT_PRIOR)
    {
      lnP += 0.0;
    }
  else assert(false);

  
  
  return(lnP);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl RW_Lk_Range(t_dsk *young, t_dsk *old, t_tree *tree)
{
  phydbl lnP;
  lnP = RRW_Lk_Range(young,old,tree);
  return(lnP);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
/* Returns probability density of observed lineage location as
   defined in Felsenstein's (1973) independent contrasts article.
*/
phydbl RW_Independent_Contrasts(t_tree *tree)
{
  phydbl lnP;
  t_node *n1,*n2;
  phydbl v1,v2,v1pv2,u,s;
  
  lnP = 0.0;

  for(int i=0;i<tree->mmod->n_dim;++i)
    {
      RW_Init_Contrasts(i,tree);
      RW_Independent_Contrasts_Post(tree->n_root,tree->n_root->v[1],tree->mmod->sigsq[i],&lnP,tree);
      RW_Independent_Contrasts_Post(tree->n_root,tree->n_root->v[2],tree->mmod->sigsq[i],&lnP,tree);
 
      n1 = tree->n_root->v[1];
      n2 = tree->n_root->v[2];
      
      // v1 and v2 in Equation (7) of Felsenstein's article     
      v1 = fabs(tree->ctrst->tprime[n1->num]-tree->times->nd_t[tree->n_root->num]);
      v2 = fabs(tree->ctrst->tprime[n2->num]-tree->times->nd_t[tree->n_root->num]);

      v1pv2 = v1+v2;
      if(v1pv2 < SMALL)
        {
          v1 = v2 = 1.;
          v1pv2 = 2.;
        }
            
      // x6
      tree->ctrst->x[tree->n_root->num] =
        (v2/(v1pv2))*tree->ctrst->x[n1->num] +
        (v1/(v1pv2))*tree->ctrst->x[n2->num] ;

      // u1
      u = tree->ctrst->x[n1->num] - tree->ctrst->x[n2->num];

      // t6'
      tree->ctrst->tprime[tree->n_root->num] = tree->times->nd_t[tree->n_root->num] + (v1*v2)/(v1pv2);
      
      s = tree->mmod->sigsq[i];
      lnP -= log(sqrt(2.*PI*s*(v1pv2)));
      lnP -= (1./(2.*s*(v1pv2)))*u*u;

      if(lnP < UNLIKELY)
        {
          PhyML_Printf("\n. v1: %f v2: %f",v1,v2);
          PhyML_Printf("\n. ctrst1: %f ctrst2: %f",tree->ctrst->x[n1->num],tree->ctrst->x[n2->num]);
          PhyML_Printf("\n. tprime: %f",tree->ctrst->tprime[tree->n_root->num]);
          PhyML_Printf("\n. t: %f",tree->times->nd_t[tree->n_root->num]);
          assert(false);
        } 
    }
  
  return(lnP);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RW_Independent_Contrasts_Post(t_node *a, t_node *d, phydbl sd, phydbl *lnP, t_tree *tree)
{
  if(d->tax == YES) return;
  else
    {
      t_node *n1,*n2;
      phydbl v1,v2,v1pv2,u;
      
      n1 = n2 = NULL;
      for(int i=0;i<3;++i)
        {
          if(d->v[i] != a && d->b[i] != tree->e_root)
            {
              if(!n1) n1 = d->v[i];
              else    n2 = d->v[i];
              RW_Independent_Contrasts_Post(d,d->v[i],sd,lnP,tree);
            }          
        }
      
      // v1 and v2 in Equation (7) of Felsenstein's article     
      v1 = fabs(tree->ctrst->tprime[n1->num]-tree->times->nd_t[d->num]);
      v2 = fabs(tree->ctrst->tprime[n2->num]-tree->times->nd_t[d->num]);

      v1pv2 = v1+v2;
      if(v1pv2 < SMALL)
        {
          v1 = v2 = 1.;
          v1pv2 = 2.;
        }

      // x6
      tree->ctrst->x[d->num] =
        (v2/(v1pv2))*tree->ctrst->x[n1->num] +
        (v1/(v1pv2))*tree->ctrst->x[n2->num] ;

      // u1
      u = tree->ctrst->x[n1->num] - tree->ctrst->x[n2->num];

      // t6'
      tree->ctrst->tprime[d->num] = tree->times->nd_t[d->num] + (v1*v2)/(v1pv2);
      
      *lnP -= log(sqrt(2.*PI*sd*(v1pv2)));
      *lnP -= (1./(2.*sd*(v1pv2)))*u*u;

      if(*lnP < UNLIKELY)
        {
          PhyML_Printf("\n. v1: %f v2: %f",v1,v2);
          PhyML_Printf("\n. ctrst1: %f ctrst2: %f",tree->ctrst->x[n1->num],tree->ctrst->x[n2->num]);
          PhyML_Printf("\n. tprime: %f",tree->ctrst->tprime[d->num]);
          PhyML_Printf("\n. t: %f",tree->times->nd_t[d->num]);
          assert(false);
        }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
void RW_Integrated_Lk_Down(phydbl son_a, phydbl son_b, phydbl son_mu_down, phydbl son_var_down, phydbl son_var,
                           phydbl bro_a, phydbl bro_b, phydbl bro_mu_down, phydbl bro_var_down, phydbl bro_var,
                           phydbl son_logrem, phydbl bro_logrem,
                           phydbl *mean, phydbl *var, phydbl *logrem)
{
  phydbl m,v,logr;
  int err;

  err = 0;
  v = 0.0;

  logr  = son_logrem + bro_logrem;
  logr -= log(fabs(bro_a*son_a));
  logr += Log_Dnorm((son_mu_down-son_b)/son_a,(bro_mu_down-bro_b)/bro_a,sqrt((son_var_down+son_var)/pow(son_a,2)+(bro_var_down+bro_var)/pow(bro_a,2)),&err);

  if((son_var + son_var_down > SMALL) && (bro_var + bro_var_down > SMALL)) // Standard case
    {
      v = pow(son_a,2)/(son_var_down + son_var) + pow(bro_a,2)/(bro_var_down + bro_var);
      v = 1/v;

      m = (son_a*(son_mu_down-son_b)/(son_var_down + son_var) + bro_a*(bro_mu_down-bro_b)/(bro_var_down + bro_var)) * v;
    }
  else if(son_var + son_var_down > SMALL) // Null variance along d - v2
    {
      m = (bro_mu_down-bro_b)/bro_a;
    }
  else if(bro_var + bro_var_down > SMALL) // Null variance along d - v1
    {
      m = (son_mu_down-son_b)/son_a; 
   }
  else
    {
      m = (son_mu_down-son_b)/son_a;
    }

  *mean = m;
  *var  = v;
  *logrem = logr;
  // PhyML_Printf("\n. m: %g v: %g logr: %g son_var_down: %g son_var: %g bro_var_down: %g bro_var: %g",m,v,logr,son_var_down,son_var,bro_var_down,bro_var);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RW_Integrated_Lk_Up(phydbl dad_mu_up, phydbl dad_var_up, phydbl dad_logrem_up,
                         phydbl son_a, phydbl son_b, phydbl son_var,
                         phydbl bro_a, phydbl bro_b, phydbl bro_mu_down, phydbl bro_var_down, phydbl bro_var, phydbl bro_logrem_down,
                         phydbl *mean, phydbl *var, phydbl *logrem)
{
  phydbl m, v, logr;
  int err;

  v = logr = 0.0;
  err = 0;

  logr = dad_logrem_up + bro_logrem_down;
  logr -= log(fabs(bro_a));
  logr += Log_Dnorm((bro_mu_down - bro_b) / bro_a, dad_mu_up, sqrt((bro_var_down + bro_var) / (bro_a * bro_a) + dad_var_up), &err);

  if (err == YES)
  {
    PhyML_Printf("\n. bro_mu_down: %f bro_b: %f bro_a: %f dad_mu_up: %f  bro_var_down: %f bro_var: %f dad_var_up: %f",
                 bro_mu_down,
                 bro_b,
                 bro_a,
                 dad_mu_up,
                 bro_var_down,
                 bro_var,
                 dad_var_up);
    assert(false);
  }

  if ((bro_var_down + bro_var) > SMALL && dad_var_up > SMALL) // Standard case
  {
    v += dad_var_up * (bro_var_down + bro_var);
    v /= bro_a * bro_a * dad_var_up + bro_var_down + bro_var;

    m = v * (bro_a * (bro_mu_down - bro_b) / (bro_var_down + bro_var) + dad_mu_up / dad_var_up);
    m = son_a * m + son_b;

    v = son_a * son_a * v + son_var;
  }
  else if (dad_var_up > SMALL) // Null variance along bro
  {
    v = son_var;
    m = son_a * bro_mu_down / bro_a + son_b;
  }
  else if ((bro_var_down + bro_var) > SMALL) // Null variance along dad
  {
    v = son_var;
    m = son_a * dad_mu_up + son_b;
  }
  else
  {
    v = son_var;
    m = son_a * bro_mu_down / bro_a + son_b;
    /* m = son_b; */
  }

  *mean = m;
  *var = v;
  *logrem = logr;

  /* !!!!!!!!!!!!!!!!!!!! */
  /* assert(v>0.0); */
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
