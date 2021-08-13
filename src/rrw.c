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
/* Relaxed random walk log-likelihood + log-prior densities for 
   relative diffusion rates on every edge */
phydbl RRW_Lk(t_tree *tree)
{
  phydbl d_fwd;
  
  d_fwd = UNLIKELY;
  
  assert(RRW_Is_Rw(tree->mmod) == YES);

  RRW_Update_Normalization_Factor(tree);

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

phydbl RRW_Prior(t_tree *tree)
{
  return(RW_Prior_Sigsq(tree) + RRW_Prior_Sigsq_Scale(tree));
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
/* Relaxed random walk log-likelihood + log-prior densities for 
   relative diffusion rates on a range of disks */
phydbl RRW_Lk_Range(t_dsk *young, t_dsk *old, t_tree *tree)
{
  RRW_Update_Normalization_Factor(tree);
  return(RRW_Forward_Lk_Range(young,old,tree));
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

phydbl RRW_Prior_Sigsq_Scale(t_tree *tree)
{
  phydbl lnP,sd;
  int err;
  
  if(tree->mmod->sigsq_scale == NULL) return(-1.);
  
  lnP = 0.0;
  err = NO;
  sd  = 2.0; // Value of hyper-prior governing the variance of relative dispersal rates across lineages
  
  if(tree->mmod->model_id == RW) return(-1.0);
                                              
  for(int i=0;i<2*tree->n_otu-2;++i)
    {
      if(tree->mmod->model_id == RRW_GAMMA)
        {
          lnP += log(Dgamma(tree->mmod->sigsq_scale[i],1./sd,sd));
        }
      else if(tree->mmod->model_id == RRW_LOGNORMAL)
        {
          lnP += Log_Dnorm(log(tree->mmod->sigsq_scale[i]),-sd*sd/2.,sd,&err);
          lnP -= log(tree->mmod->sigsq_scale[i]);
        }
    }

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
          /* PhyML_Printf("\n. RRW_LK_CORE disk: %p time: %12f lnP: %20f [%12f %12f]", */
          /*              disk,disk->time,lnP, */
          /*              disk->time,disk->ldsk->next[i]->disk->time); */
        }
    }

  return(lnP);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl RRW_Forward_Lk_Path(t_ldsk *a, t_ldsk *d, t_tree *tree)
{
  t_ldsk *ldsk;
  phydbl lnP,ld,la,disk_lnP,sd,dt;
  int i,err;
  t_node *nd_d;

  assert(a != NULL);
  assert(d != NULL);
  
  lnP = 0.0;
  
  ldsk = d;
  while(ldsk->n_next == 1) ldsk = ldsk->next[0];
  nd_d = ldsk->nd;
  
  ldsk = d;
  assert(a!=d);
  
  do
    {
      assert(ldsk->prev);

      disk_lnP = 0.0;
      
      for(i=0;i<tree->mmod->n_dim;++i)
        {
          dt = fabs(ldsk->disk->time-ldsk->prev->disk->time);

          sd =
            log(tree->mmod->sigsq[i]) +
            log(tree->mmod->sigsq_scale[nd_d->num]) +
            log(tree->mmod->rrw_norm_fact) +
            log(dt);

          sd = sqrt(exp(sd));

          /* sd = tree->mmod->sigsq[i]*tree->mmod->sigsq_scale[nd_d->num]*tree->mmod->rrw_norm_fact*dt; */
          /* sd = sqrt(sd); */
          
          ld = ldsk->coord->lonlat[i];
          la = ldsk->prev->coord->lonlat[i];

          disk_lnP += Log_Dnorm(ld,la,sd,&err);
          
          if(isinf(lnP) || isnan(lnP)) return(UNLIKELY);
          
          /* if(ld < tree->mmod->lim_do->lonlat[i] || ld > tree->mmod->lim_up->lonlat[i] || */
          /*    la < tree->mmod->lim_do->lonlat[i] || la > tree->mmod->lim_up->lonlat[i]) return(UNLIKELY); */

          if(isnan(lnP))
            {
              PhyML_Printf("\n. la=%f ld=%f sd=%f dt=[%f,%f] sigsq=%f",la,ld,sd,ldsk->disk->time,ldsk->prev->disk->time,tree->mmod->sigsq);
              assert(FALSE);
            }
          
          /* PhyML_Printf("\n. RRW_PATH time: %12f (%12f) disk_lnP: %12f sd: %12f [%12f %12f] [%12f %12f %12f %12f] %12f", */
          /*              ldsk->disk->time,ldsk->prev->disk->time, */
          /*              disk_lnP, */
          /*              sd, */
          /*              ld,la, */
          /*              tree->mmod->sigsq[i], */
          /*              tree->mmod->sigsq_scale[nd_d->num], */
          /*              tree->mmod->rrw_norm_fact, */
          /*              ldsk->disk->time-ldsk->prev->disk->time,Log_Dnorm(ld,la,sd,&err)); */
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

phydbl RRW_Density_Ldsk_Location(t_ldsk *l, phydbl rad, int dim_idx, int child_only, t_tree *tree)
{
  phydbl dta,dtd,sd,c,sumd,sum;
  int err,i;

  err = NO;  
  
  if(l->prev != NULL && l->next != NULL)
    {
      dta = fabs(l->prev->disk->time - l->disk->time);
      
      if(child_only == YES) dta = TWO_TO_THE_LARGE;

      if(dta<SMALL) return(1.0);

      sumd = 0.0;
      for(i=0;i<l->n_next;++i)
        {
          dtd = fabs(l->disk->time - l->next[i]->disk->time);
          sumd += dtd;
        }

      sum = 0.0;
      for(i=0;i<l->n_next;++i)
        {
          dtd = fabs(l->disk->time - l->next[i]->disk->time);
          sum += exp(-dtd/sumd);
        }

      
      c = 0.0;
      sd = 0.0;
      for(i=0;i<l->n_next;++i)
        {
          dtd = fabs(l->disk->time - l->next[i]->disk->time);
          /* c += (1./(phydbl)l->n_next)*((dtd/(dta+dtd))*l->prev->coord->lonlat[dim_idx] + (dta/(dta+dtd))*l->next[i]->coord->lonlat[dim_idx]); */
          /* sd += (1./(phydbl)l->n_next)*sqrt((dta*dtd)/(dta+dtd)*rad); */
          c += (exp(-dtd/sumd)/sum)*((dtd/(dta+dtd))*l->prev->coord->lonlat[dim_idx] + (dta/(dta+dtd))*l->next[i]->coord->lonlat[dim_idx]);
          sd += (exp(-dtd/sumd)/sum)*sqrt((dta*dtd)/(dta+dtd)*rad);
        }
      
      assert((isnan(c) && isnan(sd)) == false);

      /* PhyML_Printf("\n. loc: %f c: %f sd: %f dens=%f",l->coord->lonlat[dim_idx],c,sd,Log_Dnorm(l->coord->lonlat[dim_idx],c,sd,&err)); */

      return(Log_Dnorm(l->coord->lonlat[dim_idx],c,sd,&err));
    }
  else if(l->prev == NULL && l->next != NULL)
    {
      sd = sqrt(rad);
      c = l->coord->lonlat[dim_idx];
      return(Log_Dnorm(l->coord->lonlat[dim_idx],c,sd,&err));
    }
  else if(l->prev != NULL && l->next == NULL)
    {
      sd = sqrt(rad);
      c = l->coord->lonlat[dim_idx];
      return(Log_Dnorm(l->coord->lonlat[dim_idx],c,sd,&err));
    }
  else assert(false);
  
  assert(err == NO);
  return(-1.);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RRW_Generate_Ldsk_New_Location(t_ldsk *l, phydbl rad, int dim_idx, int child_only, t_tree *tree)
{
  phydbl dta,dtd,sd,c,new_loc,sum,sumd;
  int err,i;

  err = NO;              
  
  if(l->prev != NULL && l->next != NULL)
    {
      dta = fabs(l->prev->disk->time - l->disk->time);

      if(child_only == YES) dta = TWO_TO_THE_LARGE;
      
      if(dta<SMALL)
        {
          l->coord->lonlat[dim_idx] = l->prev->coord->lonlat[dim_idx];
          return;
        }
      
      sumd = 0.0;
      for(i=0;i<l->n_next;++i)
        {
          dtd = fabs(l->disk->time - l->next[i]->disk->time);
          sumd += dtd;
        }

      sum = 0.0;
      for(i=0;i<l->n_next;++i)
        {
          dtd = fabs(l->disk->time - l->next[i]->disk->time);
          sum += exp(-dtd/sumd);
        }

      c = 0.0;
      sd = 0.0;
      for(i=0;i<l->n_next;++i)
        {
          dtd = fabs(l->disk->time - l->next[i]->disk->time);
          c += (exp(-dtd/sumd)/sum)*((dtd/(dta+dtd))*l->prev->coord->lonlat[dim_idx] + (dta/(dta+dtd))*l->next[i]->coord->lonlat[dim_idx]);
          sd += (exp(-dtd/sumd)/sum)*sqrt((dta*dtd)/(dta+dtd)*rad);
          /* c += (1./(phydbl)l->n_next)*((dtd/(dta+dtd))*l->prev->coord->lonlat[dim_idx] + (dta/(dta+dtd))*l->next[i]->coord->lonlat[dim_idx]); */
          /* sd += (1./(phydbl)l->n_next)*sqrt((dta*dtd)/(dta+dtd)*rad); */
        }
      
      assert((isnan(c) && isnan(sd)) == false);

      new_loc = Rnorm(c,sd);
      
      l->coord->lonlat[dim_idx] = new_loc;
    }
  else if(l->prev == NULL && l->next != NULL)
    {
      sd = sqrt(rad);
      c = l->coord->lonlat[dim_idx];
      new_loc = Rnorm(c,sd);      
      l->coord->lonlat[dim_idx] = new_loc;          
    }
  else if(l->prev != NULL && l->next == NULL)
    {
      sd = sqrt(rad);
      c = l->coord->lonlat[dim_idx];
      new_loc = Rnorm(c,sd);      
      l->coord->lonlat[dim_idx] = new_loc;
    }
  else assert(false);
  assert(err == NO);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

short int RRW_Is_Rw(t_phyrex_mod *mod)
{
  if(mod->model_id == RW ||
     mod->model_id == RRW_GAMMA ||
     mod->model_id == RRW_LOGNORMAL) return(YES);
  return(NO);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RRW_Update_Normalization_Factor(t_tree *tree)
{
  phydbl dt,rdt,T,RT;
  int i;

  rdt = 0.0;
  dt  = 0.0;
  T   = 0.0;
  RT  = 0.0;
  
  for(i=0;i<2*tree->n_otu-2;++i)
    {
      assert(tree->a_nodes[i] != tree->n_root);
      dt = fabs(tree->times->nd_t[i] - tree->times->nd_t[tree->a_nodes[i]->anc->num]);
      rdt = dt*tree->mmod->sigsq_scale[tree->a_nodes[i]->num];
      T+=dt;
      RT+=rdt;
    }
  tree->mmod->rrw_norm_fact = T/RT;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl RRW_Mean_Displacement_Rate(t_tree *tree)
{
  phydbl dt,disp;

  disp = 0.0;
  dt = 0.0;
  for(int i=0;i<2*tree->n_otu-2;++i)
    {
      assert(tree->a_nodes[i] != tree->n_root);
      dt += fabs(tree->times->nd_t[i] - tree->times->nd_t[tree->a_nodes[i]->anc->num]);
      disp +=
        fabs(tree->times->nd_t[i] - tree->times->nd_t[tree->a_nodes[i]->anc->num]) *
        tree->mmod->sigsq_scale[tree->a_nodes[i]->num]*
        (tree->mmod->sigsq[0]+tree->mmod->sigsq[1])*
        tree->mmod->rrw_norm_fact;
    }
  return(disp/dt);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RRW_Generate(t_tree *tree)
{
  t_ldsk *ldsk;
  t_dsk *disk;
  int i,j;

  ldsk = tree->n_root->ldsk;
  disk = tree->n_root->ldsk->disk;
  
  for(j=0;j<tree->mmod->n_dim;++j) ldsk->coord->lonlat[j] = 0.5*(tree->mmod->lim_up->lonlat[j] + tree->mmod->lim_do->lonlat[j]);;

  /* disk = disk->next; */
  do
    {
      if(disk->ldsk != NULL)
        {
          ldsk = disk->ldsk;
          
          for(i=0;i<ldsk->n_next;++i)
            {
              for(j=0;j<tree->mmod->n_dim;++j)
                {
                  ldsk->next[i]->coord->lonlat[j] = Rnorm(ldsk->coord->lonlat[j],
                                                          sqrt((ldsk->next[i]->disk->time - ldsk->disk->time)*tree->mmod->sigsq[j]));
                }
            }
        }
      
      disk = disk->next;
    }
  while(disk);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
