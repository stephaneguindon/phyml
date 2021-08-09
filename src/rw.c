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
  return(RW_Prior_Sigsq(tree));
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl RW_Prior_Sigsq(t_tree *tree)
{
  phydbl lbda,lnP;

  lnP = 0.0;
  lbda = 1.0 / tree->mmod->disp_prior_mean;
  for(int i=0;i<tree->mmod->n_dim;++i) lnP += log(lbda)-lbda*tree->mmod->sigsq[i];
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
  phydbl v1,v2,u,s;
  
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

      // x6
      tree->ctrst->x[tree->n_root->num] =
        (v2/(v2+v1))*tree->ctrst->x[n1->num] +
        (v1/(v2+v1))*tree->ctrst->x[n2->num] ;

      // u1
      u = tree->ctrst->x[n1->num] - tree->ctrst->x[n2->num];

      // t6'
      tree->ctrst->tprime[tree->n_root->num] = tree->times->nd_t[tree->n_root->num] + (v1*v2)/(v1+v2);
      
      s = tree->mmod->sigsq[i];
      lnP -= log(sqrt(2.*PI*s*(v1+v2)));
      lnP -= (1./(2.*s*(v1+v2)))*u*u;
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
      phydbl v1,v2,u;
      
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

      // x6
      tree->ctrst->x[d->num] =
        (v2/(v2+v1))*tree->ctrst->x[n1->num] +
        (v1/(v2+v1))*tree->ctrst->x[n2->num] ;

      // u1
      u = tree->ctrst->x[n1->num] - tree->ctrst->x[n2->num];

      // t6'
      tree->ctrst->tprime[d->num] = tree->times->nd_t[d->num] + (v1*v2)/(v1+v2);
      
      *lnP -= log(sqrt(2.*PI*sd*(v1+v2)));
      *lnP -= (1./(2.*sd*(v1+v2)))*u*u;
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
