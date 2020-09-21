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
  phydbl d_fwd,d_coal,d_sigsq_scale;
  phydbl t_dum;
  int idx_dum;
  
  assert(tree->mmod->model_id == RRW);


  d_fwd         = 0.0;
  d_coal        = 0.0;
  d_sigsq_scale = 0.0;
  
  idx_dum = Rand_Int(0,2*tree->n_otu-3);
  t_dum = tree->times->nd_t[idx_dum];

  d_fwd = RRW_Forward_Lk_Range(tree->young_disk,NULL,tree);

#ifdef PHYREX
  if(PHYREX_Total_Number_Of_Intervals(tree) > tree->mmod->max_num_of_intervals)
    {
      tree->mmod->c_lnL = UNLIKELY;
      return(UNLIKELY);
    }
#endif

  d_coal = TIMES_Lk_Coalescent(tree);

  d_sigsq_scale = RRW_Prior_Sigsq_Scale(tree);

  
  // Make sure node times are set back to their original values
  assert(fabs(t_dum - tree->times->nd_t[idx_dum]) < 1.E-4);

  tree->mmod->c_lnL = d_fwd + d_coal + d_sigsq_scale;
  
  return(tree->mmod->c_lnL);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl RRW_Lk_Range(t_dsk *young, t_dsk *old, t_tree *tree)
{
  phydbl lnP;
  
  lnP = 0.0;
  
  /* lnP += RRW_Forward_Lk_Range(young,old,tree); */
  /* lnP += TIMES_Lk_Coalescent_Range(young,old,tree); */

  /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
  lnP += RRW_Forward_Lk_Range(young,NULL,tree);
  lnP += TIMES_Lk_Coalescent_Range(young,old,tree);


  /* PhyML_Printf("\n. RANGE = %f",TIMES_Lk_Coalescent_Range(young,old,tree)); */
  /* PhyML_Printf("\n. FULL = %f",TIMES_Lk_Coalescent(tree)); */
  /* Exit("\n"); */

  
  return(lnP);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl RRW_Lk_Core(t_dsk *disk, t_tree *tree)
{
  phydbl lnP;
  int i;

  assert(disk);
  
  if(disk == tree->young_disk) return 0.0;
  if(disk->age_fixed == YES) return 0.0;

  lnP = 0.0;
  
  if(disk->ldsk != NULL)
    {
      for(i=0;i<disk->ldsk->n_next;++i)
        {
          lnP += RRW_Forward_Lk_Path(disk->ldsk,disk->ldsk->next[i],tree);
        }
    }

  return(lnP);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl RRW_Independent_Contrasts(t_tree *tree)
{
  phydbl lnL;
  
  RRW_Rescale_Times(YES,tree);
  lnL = RW_Independent_Contrasts(tree);
  RRW_Rescale_Times(NO,tree);

  return(lnL);
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

phydbl RRW_Forward_Lk_Range(t_dsk *young, t_dsk *old, t_tree *tree)
{
  phydbl lnP;
  int i;
  t_dsk *disk;
  
  lnP = 0.0;
  disk = young;
  
  do
    {
      lnP += RRW_Lk_Core(disk,tree);
      
      if(old && disk == old) break;
      
      disk = disk->prev;
    }  
  while(disk);

  return(lnP);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl RRW_Forward_Lk_Path(t_ldsk *a, t_ldsk *d, t_tree *tree)
{
  t_ldsk *ldsk;
  phydbl lnP,sd,ld,la,eps;
  int i,err;

  lnP = 0.0;
  eps = 1.E-5;
  ldsk = d;

  assert(!(a->disk->time > d->disk->time));
  assert(a!=d);
  
  do
    {
      assert(ldsk->prev);

      sd = log(tree->mmod->sigsq) + log(fabs(ldsk->disk->time-ldsk->prev->disk->time));
      sd = exp(sd);
      
      for(i=0;i<tree->mmod->n_dim;++i)
        {
          ld = ldsk->coord->lonlat[i];
          la = ldsk->prev->coord->lonlat[i];

          if(ld > tree->mmod->lim_up->lonlat[i] || ld < tree->mmod->lim_do->lonlat[i])
            {
              return UNLIKELY;
            }
          if(la > tree->mmod->lim_up->lonlat[i] || la < tree->mmod->lim_do->lonlat[i])
            {
              return UNLIKELY;
            }
          
          assert(!Are_Equal(ld,0.0,1.E-5));
          assert(!Are_Equal(la,0.0,1.E-5));

          lnP += Log_Dnorm(ld,la,sd,&err);
        }

      ldsk = ldsk->prev;
      assert(ldsk);
    }
  while(ldsk != a);

  return(lnP);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
