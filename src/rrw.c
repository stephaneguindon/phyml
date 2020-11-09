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
  phydbl d_fwd;
  
  assert(tree->mmod->model_id == RRW);
  
  d_fwd = RRW_Forward_Lk_Range(tree->young_disk,NULL,tree);
  
#ifdef PHYREX
  if(PHYREX_Total_Number_Of_Intervals(tree) > tree->mmod->max_num_of_intervals)
    {
      tree->mmod->c_lnL = UNLIKELY;
      return(UNLIKELY);
    }
#endif

  tree->mmod->c_lnL = d_fwd;

  return(tree->mmod->c_lnL);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl RRW_Lk_Range(t_dsk *young, t_dsk *old, t_tree *tree)
{
  phydbl lnP;
  lnP = RRW_Forward_Lk_Range(young,old,tree);
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
  phydbl lnP,disk_lnP;
  t_dsk *disk;

  disk_lnP = 0.0;
  lnP = 0.0;
  disk = young;
  
  do
    {
      disk_lnP = RRW_Lk_Core(disk,tree);
      lnP += disk_lnP;
      
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
  phydbl lnP,sd,ld,la,disk_lnP;
  int i,err;

  lnP = 0.0;
  ldsk = d;

  assert(a!=d);
  
  do
    {
      assert(ldsk->prev);

      disk_lnP = 0.0;
      for(i=0;i<tree->mmod->n_dim;++i)
        {
          sd = log(tree->mmod->sigsq[i]) + log(fabs(ldsk->disk->time-ldsk->prev->disk->time));
          sd = sqrt(exp(sd));
          
          ld = ldsk->coord->lonlat[i];
          la = ldsk->prev->coord->lonlat[i];

          /* if(ld > tree->mmod->lim_up->lonlat[i] || */
          /*    ld < tree->mmod->lim_do->lonlat[i] || */
          /*    la > tree->mmod->lim_up->lonlat[i] || */
          /*    la < tree->mmod->lim_do->lonlat[i]) return UNLIKELY; */
          
          /* assert(!Are_Equal(ld,0.0,1.E-5)); */
          /* assert(!Are_Equal(la,0.0,1.E-5)); */

          disk_lnP += Log_Dnorm(ld,la,sd,&err);

          if(isinf(lnP)) return(UNLIKELY);

          if(isnan(lnP))
            {
              PhyML_Printf("\n. la=%f ld=%f sd=%f dt=[%f,%f] sigsq=%f",la,ld,sd,ldsk->disk->time,ldsk->prev->disk->time,tree->mmod->sigsq);
              assert(FALSE);
            }
        }

      lnP += disk_lnP;
            
      ldsk = ldsk->prev;
      assert(ldsk);
    }
  while(ldsk != a);

  return(lnP);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RRW_Generate_Ldsk_New_Location(t_ldsk *l, t_ldsk *prev_l, phydbl rad, phydbl *hr, int dim_idx, t_tree *tree)
{
  phydbl dta,dtd,sum,sd,c,new_loc,w;
  int err,i;

  err = NO;              
  
  if(l->prev != NULL)
    {
      dta = fabs(l->prev->disk->time - l->disk->time);
      
      if(dta<SMALL)
        {
          l->coord->lonlat[dim_idx] = l->prev->coord->lonlat[dim_idx];
          return;
        }
      
      sum = 0.0;
      for(i=0;i<l->n_next;++i)
        {
          dtd = fabs(l->disk->time - l->next[i]->disk->time);
          w = dta/(dta+dtd);
          sum += w;
        }

      c = 0.0;
      sd = 0.0;
      for(i=0;i<l->n_next;++i)
        {
          dtd = fabs(l->disk->time - l->next[i]->disk->time);
          w = dta/(dta+dtd)/sum;
          c += w*((dtd/(dta+dtd))*l->prev->coord->lonlat[dim_idx] + (dta/(dta+dtd))*l->next[i]->coord->lonlat[dim_idx]);
          sd += w*sqrt((dta*dtd)/(dta+dtd)*rad);
        }
      
      assert((isnan(c) && isnan(sd)) == false);

      new_loc = Rnorm(c,sd);
      
      /* Forward move */
      (*hr) -= Log_Dnorm(new_loc,c,sd,&err);
      
      if(prev_l != NULL && prev_l != l)
        {
          /* Reverse move */
          dta = fabs(prev_l->prev->disk->time - prev_l->disk->time);

          sum = 0.0;
          for(i=0;i<prev_l->n_next;++i)
            {
              dtd = fabs(prev_l->disk->time - prev_l->next[i]->disk->time);
              w = dta/(dta+dtd);
              sum += w;
            }
          
          c = 0.0;
          sd = 0.0;
          for(i=0;i<prev_l->n_next;++i)
            {
              dtd = fabs(prev_l->disk->time - prev_l->next[i]->disk->time);
              w = dta/(dta+dtd)/sum;
              c += w*((dtd/(dta+dtd))*prev_l->prev->coord->lonlat[dim_idx] + (dta/(dta+dtd))*prev_l->next[i]->coord->lonlat[dim_idx]);
              sd += w*sqrt((dta*dtd)/(dta+dtd)*rad);
            }
        }
          
      assert((isnan(c) && isnan(sd)) == false);

      (*hr) += Log_Dnorm(prev_l->coord->lonlat[dim_idx],c,sd,&err);
      l->coord->lonlat[dim_idx] = new_loc;
    }
  else
    {
      sd = sqrt(rad);
      c = l->coord->lonlat[dim_idx];
  
      new_loc = Rnorm(c,sd);
      
      /* Forward move */
      (*hr) -= Log_Dnorm(new_loc,c,sd,&err);
      
      if(prev_l != NULL)
        {
          c = new_loc;
          
          /* Reverse move */
          (*hr) += Log_Dnorm(prev_l->coord->lonlat[dim_idx],c,sd,&err); 
        }
          
      l->coord->lonlat[dim_idx] = new_loc;          
    }
  assert(err == NO);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
