/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include "bmp.h"

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl BMP_Lk(t_tree *tree)
{
  assert(tree->mmod->id == BMP);

  tree->mmod->c_lnL = 0.0;
  tree->mmod->c_lnL += BMP_Forward_Lk(tree);
  tree->mmod->c_lnL -= BMP_Independent_Contrasts(tree);

  
  return(tree->mmod->c_lnL);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
/* Returns probability density of observed lineage location as
   defined in Felsenstein's (1973) independent contrasts article.
*/
phydbl BMP_Independent_Contrasts(t_tree *tree)
{
  phydbl lnP;
  
  lnP = 0.0;

  for(int i=0;i<tree->mmod->n_dim;++i)
    {
      BMP_Init_Contrasts(i,tree);
      BMP_Independent_Contrasts_Post(tree->n_root,tree->n_root->v[1],&lnP,tree);
      BMP_Independent_Contrasts_Post(tree->n_root,tree->n_root->v[2],&lnP,tree);
    }
  
  return(lnP);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void BMP_Independent_Contrasts_Post(t_node *a, t_node *d, phydbl *lnP, t_tree *tree)
{
  if(d->tax == YES) return;
  else
    {
      t_node *n1,*n2;
      phydbl v1,v2,u,s;
      
      n1 = n2 = NULL;
      for(int i=0;i<3;++i)
        {
          if(d->v[i] != a && d->b[i] != tree->e_root)
            {
              if(!n1) n1 = d->v[i];
              else    n2 = d->v[i];
              BMP_Independent_Contrasts_Post(d,d->v[i],lnP,tree);
            }          
        }
      
      // v1 and v2 in Equation (7) of Felsenstein's article     
      v1 = fabs(tree->ctrst->tprime[n1->num]-tree->rates->nd_t[d->num]);
      v2 = fabs(tree->ctrst->tprime[n2->num]-tree->rates->nd_t[d->num]);

      // x6
      tree->ctrst->x[d->num] =
        (v2/(v2+v1))*tree->ctrst->x[n1->num] +
        (v1/(v2+v1))*tree->ctrst->x[n2->num] ;

      // u1
      u = tree->ctrst->x[n1->num] - tree->ctrst->x[n2->num];

      // t6'
      tree->ctrst->tprime[d->num] = tree->rates->nd_t[d->num] + (v1*v2)/(v1+v2);
      
      s = tree->mmod->sigsq;
      *lnP -= log(sqrt(2.*PI*s*(v1+v2)));
      *lnP -= (1./(2.*s))*u*u;
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl BMP_Forward_Lk(t_tree *tree)
{
  phydbl lnP = 0.0;
  BMP_Forward_Lk_Pre(tree->n_root,tree->n_root->v[1],&lnP,tree);
  BMP_Forward_Lk_Pre(tree->n_root,tree->n_root->v[2],&lnP,tree);
  for(int i=0;i<tree->mmod->n_dim;++i) lnP -= log(tree->mmod->lim_up->lonlat[i]-tree->mmod->lim_do->lonlat[i]);
  return(lnP);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void BMP_Forward_Lk_Pre(t_node *a, t_node *d, phydbl *lnP, t_tree *tree)
{
  int i,err;
  phydbl ld,la,sd;
  
  for(i=0;i<tree->mmod->n_dim;++i)
    {
      ld = d->ldsk->coord->lonlat[i];
      la = a->ldsk->coord->lonlat[i];
      sd = sqrt(tree->mmod->sigsq*fabs(tree->rates->nd_t[a->num]-tree->rates->nd_t[d->num]));
      
      *lnP += Log_Dnorm_Trunc(ld,
                              la,
                              sd,
                              tree->mmod->lim_do->lonlat[i],
                              tree->mmod->lim_up->lonlat[i],&err);

    }
  
  if(d->tax) return;
  else
    {
      for(i=0;i<3;++i)
        {
          if(d->v[i] != a && d->b[i] != tree->e_root)
            {
              BMP_Forward_Lk_Pre(d,d->v[i],lnP,tree);
            }
        }
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
