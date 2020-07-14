/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include "rrw.h"

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl RRW_Lk(t_tree *tree)
{
  phydbl d_fwd,d_ic,d_coal,d_sigsq_scale;
  phydbl t_dum;
  int idx_dum;
      
  assert(tree->mmod->id == RRW);

  #ifdef PHYREX
  if(PHYREX_Total_Number_Of_Intervals(tree) > tree->mmod->max_num_of_intervals) return UNLIKELY;
  #endif
  
  idx_dum = Rand_Int(0,2*tree->n_otu-3);
  t_dum = tree->times->nd_t[idx_dum];
  d_fwd = RRW_Forward_Lk(tree);
  d_coal = TIMES_Lk_Coalescent(tree);
  RRW_Rescale_Times(YES,tree);
  d_ic = RW_Independent_Contrasts(tree);
  RRW_Rescale_Times(NO,tree);
  d_sigsq_scale = RRW_Prior_Sigsq_Scale(tree);
  
  // Make sure node times are set back to their original values
  assert(fabs(t_dum - tree->times->nd_t[idx_dum]) < 1.E-4);

  /* tree->mmod->c_lnL = d_fwd - d_ic + d_coal + d_sigsq_scale; */
  tree->mmod->c_lnL = d_fwd + d_coal + d_sigsq_scale;
    
  return(tree->mmod->c_lnL);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl RRW_Prior_Sigsq_Scale(t_tree *tree)
{
  phydbl lnP;

  lnP = 0.0;
  
  for(int i=0;i<2*tree->n_otu-2;++i)
    {
      lnP += log(Dgamma(tree->mmod->sigsq_scale[i],
                        tree->mmod->nu/2.,
                        2./tree->mmod->nu));
    }
  
  return(lnP);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RRW_Rescale_Times(int prod, t_tree *tree)
{
  RRW_Rescale_Times_Pre(tree->n_root,tree->n_root->v[1],tree->times->nd_t[tree->n_root->num],prod,tree);
  RRW_Rescale_Times_Pre(tree->n_root,tree->n_root->v[2],tree->times->nd_t[tree->n_root->num],prod,tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RRW_Rescale_Times_Pre(t_node *a, t_node *d, phydbl cur_ta, int prod, t_tree *tree)
{
  phydbl ta,td;
  int i;
  
  td = tree->times->nd_t[d->num];
  ta = tree->times->nd_t[a->num];

  if(prod == YES)
    tree->times->nd_t[d->num] = ta + (td-cur_ta) * tree->mmod->sigsq_scale[d->num];
  else
    tree->times->nd_t[d->num] = ta + (td-cur_ta) / tree->mmod->sigsq_scale[d->num];
  
  if(d->tax == YES) return;
  else for(i=0;i<3;++i) if(d->v[i] != a && d->b[i] != tree->e_root) RRW_Rescale_Times_Pre(d,d->v[i],td,prod,tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl RRW_Forward_Lk(t_tree *tree)
{
  phydbl lnP,la,ld,sd;
  int i,err;

  lnP = 0.0;
  RRW_Forward_Lk_Pre(tree->n_root,tree->n_root->v[1],&lnP,tree);
  RRW_Forward_Lk_Pre(tree->n_root,tree->n_root->v[2],&lnP,tree);
  

  for(i=0;i<tree->mmod->n_dim;++i)
    {
      ld = tree->n_root->v[1]->ldsk->coord->lonlat[i];
      la = tree->n_root->ldsk->coord->lonlat[i];
      sd = sqrt(tree->mmod->sigsq *
                tree->mmod->sigsq_scale[tree->n_root->v[1]->num] *
                fabs(tree->times->nd_t[tree->n_root->v[1]->num]-tree->times->nd_t[tree->n_root->num]));
      
      lnP += Log_Dnorm(ld,la,sd,&err);

      
      ld = tree->n_root->v[2]->ldsk->coord->lonlat[i];
      la = tree->n_root->ldsk->coord->lonlat[i];
      sd = sqrt(tree->mmod->sigsq*
                tree->mmod->sigsq_scale[tree->n_root->v[2]->num] *
                fabs(tree->times->nd_t[tree->n_root->v[2]->num]-tree->times->nd_t[tree->n_root->num]));
      
      lnP += Log_Dnorm(ld,la,sd,&err);
    }
  
  for(i=0;i<tree->mmod->n_dim;++i) lnP -= log(tree->mmod->lim_up->lonlat[i]-tree->mmod->lim_do->lonlat[i]);
  return(lnP);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RRW_Forward_Lk_Pre(t_node *a, t_node *d, phydbl *lnP, t_tree *tree)
{
  int i,err;
  phydbl ld,la,sd;
  
  sd = sqrt(tree->mmod->sigsq *
            tree->mmod->sigsq_scale[d->num] *
            fabs(tree->times->nd_t[d->num]-tree->times->nd_t[a->num]));

  for(i=0;i<tree->mmod->n_dim;++i)
    {
      ld = d->ldsk->coord->lonlat[i];
      la = a->ldsk->coord->lonlat[i];
      
      *lnP += Log_Dnorm(ld,la,sd,&err);
    }
  
  if(d->tax) return;
  else
    {
      for(i=0;i<3;++i)
        {
          if(d->v[i] != a && d->b[i] != tree->e_root)
            {
              RRW_Forward_Lk_Pre(d,d->v[i],lnP,tree);
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
