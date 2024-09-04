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
  
  if(tree->mmod->integrateAncestralLocations == YES)
    {
      tree->mmod->c_lnL = RRW_Integrated_Lk(tree);
    }
  else
    {
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
    }
  
  return(tree->mmod->c_lnL);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl RRW_Prior(t_tree *tree)
{
  tree->mmod->c_lnP = RW_Prior_Sigsq(tree) + RRW_Prior_Sigsq_Scale(tree);
  return(tree->mmod->c_lnP);
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
  sd  = tree->mmod->rrw_prior_sd; // Value of hyper-prior governing the variance of relative dispersal rates across lineages
  
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
            log(tree->mmod->sigsq_scale_norm_fact) +
            log(dt);

          sd = sqrt(exp(sd));

          
          ld = ldsk->coord->lonlat[i];
          la = ldsk->prev->coord->lonlat[i];

          disk_lnP += Log_Dnorm(ld,la,sd,&err);
          
          if(isinf(lnP) || isnan(lnP)) return(UNLIKELY);
          
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
          /*              tree->mmod->sigsq_scale_norm_fact, */
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

  RRW_Update_Normalization_Factor(tree);
  TIMES_Record_Times(tree);
  RRW_Rescale_Times(YES,tree);
  lnL = RW_Independent_Contrasts(tree);
  TIMES_Reset_Times(tree);

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

  assert(!(cur_ta > td));

  if(prod == YES)
    tree->times->nd_t[d->num] = ta + (td-cur_ta) * (tree->mmod->sigsq_scale[d->num] * tree->mmod->sigsq_scale_norm_fact);
  else
    tree->times->nd_t[d->num] = ta + (td-cur_ta) / (tree->mmod->sigsq_scale[d->num] * tree->mmod->sigsq_scale_norm_fact);

  assert(!(tree->times->nd_t[d->num] < ta));
  
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
  /* !!!!!!!!!!!!!!!!!!!!! */
  tree->mmod->sigsq_scale_norm_fact = 1.0;
  
  /* phydbl dt,rdt,T,RT; */
  /* int i; */

  /* rdt = 0.0; */
  /* dt  = 0.0; */
  /* T   = 0.0; */
  /* RT  = 0.0; */
  
  /* for(i=0;i<2*tree->n_otu-2;++i) */
  /*   { */
  /*     assert(tree->a_nodes[i] != tree->n_root); */
  /*     dt = fabs(tree->times->nd_t[i] - tree->times->nd_t[tree->a_nodes[i]->anc->num]); */
  /*     rdt = dt*tree->mmod->sigsq_scale[tree->a_nodes[i]->num]; */
  /*     T+=dt; */
  /*     RT+=rdt; */
  /*   } */
  /* tree->mmod->sigsq_scale_norm_fact = T/RT; */
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
        tree->mmod->sigsq_scale_norm_fact;
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

void RRW_Sample_Arealin_Plot(t_tree *tree)
{
  int i,j;
  t_ldsk *ldsk;
  t_dsk *disk;
  t_geo_coord **coord_array;


  coord_array = (t_geo_coord **)mCalloc(tree->n_otu,sizeof(t_geo_coord *));
  for(j=0;j<tree->n_otu;++j) 
    {
      coord_array[j] = GEO_Make_Geo_Coord(tree->mmod->n_dim);
      GEO_Init_Coord(coord_array[j],tree->mmod->n_dim);      
    }
  
  RRW_Update_Normalization_Factor(tree);
  
  Init_Contmod_Locations(tree);

  for(i=0;i<tree->mmod->n_dim;++i)
    {
      RRW_Integrated_Lk_Location_Post(NULL,tree->n_root,tree);
      RRW_Integrated_Lk_Location_Pre(tree->n_root,tree->n_root->v[1],tree);
      RRW_Integrated_Lk_Location_Pre(tree->n_root,tree->n_root->v[2],tree);
    }
  
  /* Sample location on every ldsk */
  disk = tree->young_disk;      
  do
    {
      for(j=0;j<disk->n_ldsk_a;++j)
        {
          if(disk->ldsk_a[j] == disk->ldsk)
            {
              assert(disk->ldsk->nd != NULL && disk->ldsk->nd->tax == NO);
              
              for(i=0;i<tree->mmod->n_dim;++i)
                {
                  coord_array[j]->lonlat[i] =
                    Sample_Ancestral_Trait_Contmod(disk->ldsk->nd->anc,
                                                   disk->ldsk->nd,
                                                   fabs(tree->times->nd_t[disk->ldsk->nd->anc->num]-tree->times->nd_t[disk->ldsk->nd->num]),
                                                   0.0,
                                                   log(tree->mmod->sigsq[i]) +
                                                   log(tree->mmod->sigsq_scale[disk->ldsk->nd->num]) + 
                                                   log(tree->mmod->sigsq_scale_norm_fact),                                                   
                                                   NO,
                                                   tree);
                }
            }
          else if(disk->age_fixed == NO || (disk->age_fixed == YES && disk->ldsk_a[j]->disk != disk))
            {
              t_ldsk *ldsk_a,*ldsk_d;
              
              ldsk_a = NULL;
              ldsk_d = NULL;
              
              ldsk = disk->ldsk_a[j];
              do ldsk = ldsk->prev; while(ldsk->nd == NULL);
              ldsk_a = ldsk;
              
              ldsk_d = disk->ldsk_a[j];
              
              assert(ldsk_d->nd != NULL);
              
              for(i=0;i<tree->mmod->n_dim;++i)
                {
                  coord_array[j]->lonlat[i] =
                    Sample_Ancestral_Trait_Contmod(ldsk_a->nd,
                                                   ldsk_d->nd,
                                                   fabs(disk->time-tree->times->nd_t[ldsk_a->nd->num]),                                        
                                                   fabs(disk->time-tree->times->nd_t[ldsk_d->nd->num]),
                                                   log(tree->mmod->sigsq[i]) +
                                                   log(tree->mmod->sigsq_scale[ldsk_d->nd->num]) + 
                                                   log(tree->mmod->sigsq_scale_norm_fact),                                                   
                                                   NO,
                                                   tree);
                }
            }
          else 
            {
              assert(disk->age_fixed == YES);
              for(i=0;i<tree->mmod->n_dim;++i) coord_array[j]->lonlat[i] = disk->ldsk_a[j]->coord->lonlat[i];                    
            }
        }
      
      /* Get barycenter, radius and then area for every disk */
      phydbl max_dist,dist;
      t_geo_coord *barycentr;
      
      /* PhyML_Printf("\n<> "); */

      barycentr = GEO_Make_Geo_Coord(tree->mmod->n_dim);
      GEO_Init_Coord(barycentr,tree->mmod->n_dim);
          
      for(i=0;i<tree->mmod->n_dim;++i) 
        {
          barycentr->lonlat[i] = 0.0;
          for(j=0;j<disk->n_ldsk_a;++j) barycentr->lonlat[i] += coord_array[j]->lonlat[i];
          barycentr->lonlat[i] /= disk->n_ldsk_a;
        }
      
      max_dist = -1.;
      dist = -1.;
      for(j=0;j<disk->n_ldsk_a;++j)
        {
          dist = Euclidean_Distance(barycentr,coord_array[j]);
          if(dist > max_dist) max_dist = dist;
        }
      
            
      assert(!(max_dist < 0.0));
            
      /* PhyML_Printf("(%f;%f;%d)",area,disk->time,disk->n_ldsk_a); */
      
      Free_Geo_Coord(barycentr);
      
      disk = disk->prev;
    }
  while(disk);
  
  for(j=0;j<tree->n_otu;++j) Free_Geo_Coord(coord_array[j]);
  Free(coord_array);  
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RRW_Sample_Node_Locations_Joint(t_tree *tree)
{
  int i,j,start;
  phydbl root_var, root_mean, var, mean;
  
  root_var = tree->mmod->rw_root_var;
  root_mean = -BIG;
  
  RRW_Update_Normalization_Factor(tree);

  Init_Contmod_Locations(tree);

  RRW_Integrated_Lk_Location_Post(NULL, tree->n_root, tree);

  for(i=0;i<tree->mmod->n_dim;++i)
    {
      start = Contmod_Start(LOCATION,i,tree);
      
      var = 1./tree->contmod->var_down[start+tree->n_root->num] + 1./root_var;
      var = 1./var;
      
      root_mean = LOCATION_Mean_Lonlat(i,tree);

      mean = (tree->contmod->mu_down[start+tree->n_root->num]/tree->contmod->var_down[start+tree->n_root->num] +
              root_mean / root_var)*var;

      tree->n_root->ldsk->coord->lonlat[i] =
        Rnorm(tree->contmod->mu_down[start+tree->n_root->num],
              sqrt(tree->contmod->var_down[start+tree->n_root->num]));

      assert(isnan(tree->contmod->var_down[start + tree->n_root->num]) == NO);
    }

  RRW_Sample_Node_Locations_Joint_Post(tree->n_root, tree->n_root->v[1], tree);
  RRW_Sample_Node_Locations_Joint_Post(tree->n_root, tree->n_root->v[2], tree);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RRW_Sample_Node_Locations_Joint_Post(t_node *a, t_node *d, t_tree *tree)
{

  if(d->tax == TRUE) return;
  else
  {
    int start, i;
    phydbl son_a, son_b, son_var, son_var_down, son_mu, son_mu_down, var, mean;

    for (int i = 0; i < 3; ++i)
      if (d->v[i] != a && !(a == tree->n_root && d->b[i] == tree->e_root))
        RRW_Sample_Node_Locations_Joint_Post(d, d->v[i], tree);

    for (i = 0; i < tree->mmod->n_dim; ++i)
    {
      start = Contmod_Start(LOCATION, i, tree);

      son_a = 1.0;
      son_b = 0.0;

      son_var = RRW_Location_Variance_Along_Edge(d, i, tree);
      son_var_down = tree->contmod->var_down[start + d->num];

      son_mu = son_a * a->ldsk->coord->lonlat[i] + son_b;
      son_mu_down = tree->contmod->mu_down[start + d->num];

      if (son_var > SMALL && son_var_down > SMALL)
      {
        var = 1. / son_var + 1. / son_var_down;
        var = 1. / var;

        mean = (son_mu_down / son_var_down + son_mu / son_var) * var;
      }
      else if (son_var_down > SMALL)
      {
        var = 0.0;
        mean = son_mu;
      }
      else if (son_var > SMALL)
      {
        var = 0.0;
        mean = son_mu_down;
      }
      else
      {
        var = 0.0;
        mean = son_mu;
      }

      assert(!(var < 0.0));
      d->ldsk->coord->lonlat[i] = Rnorm(mean, sqrt(var));
    }
  }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RRW_Sample_Node_Locations_Marginal(t_tree *tree)
{
  int i,j,start;
  t_node *n;
  phydbl au,bu,varu,var,mean;
  phydbl root_var, root_mean;
  
  n = NULL;
  au = bu = varu = var = mean = -1;

  root_var = tree->mmod->rw_root_var;
  root_mean = -BIG;
  
  RRW_Update_Normalization_Factor(tree);

  Init_Contmod_Locations(tree);

  RRW_Integrated_Lk_Location_Post(NULL, tree->n_root, tree);
  RRW_Integrated_Lk_Location_Pre(tree->n_root, tree->n_root->v[1], tree);
  RRW_Integrated_Lk_Location_Pre(tree->n_root, tree->n_root->v[2], tree);      

  for(i=0;i<tree->mmod->n_dim;++i)
    {
      start = Contmod_Start(LOCATION,i,tree);

      for(j=0;j<2*tree->n_otu-2;++j)
        {
          if(tree->a_nodes[j]->tax == NO && tree->a_nodes[j] != tree->n_root)
            {
              n = tree->a_nodes[j];
              au = 1.0;
              bu = 0.0;
              varu = RRW_Location_Variance_Along_Edge(n,i,tree); 
        
              if(tree->contmod->var_down[start+n->num] > SMALL && pow(au,2)*tree->contmod->var_up[start+n->num]+varu > SMALL)
                {
                  var = 1./tree->contmod->var_down[start+n->num] + 1. / (pow(au,2)*tree->contmod->var_up[start+n->num]+varu);
                  var = 1./var;
                  mean = (tree->contmod->mu_down[start+n->num]/tree->contmod->var_down[start+n->num] + (au*tree->contmod->mu_up[start+n->num]+bu)/(pow(au,2)*tree->contmod->var_up[start+n->num]+varu))*var;
                }
              else if(tree->contmod->var_down[start+n->num] > SMALL)
                {
                  var = 0.0;
                  mean = au*tree->contmod->mu_up[start+n->num]+bu;
                }
              else if(pow(au,2)*tree->contmod->var_up[start+n->num]+varu > SMALL)
                {
                  var = 0.0;
                  mean = tree->contmod->mu_down[start+n->num];
                }
              else
                {
                  var = 0.0;
                  mean = tree->contmod->mu_down[start+n->num];
                }
              
              /* Below is only valid if prior distribution of location at that node is flat */
              tree->a_nodes[j]->ldsk->coord->lonlat[i] = Rnorm(mean,sqrt(var));
            }
        }

      assert(isnan(tree->contmod->var_down[start+tree->n_root->num]) == NO);
      
      var = 1./tree->contmod->var_down[start+tree->n_root->num] + 1./root_var;
      var = 1./var;
      
      root_mean = LOCATION_Mean_Lonlat(i,tree);

      mean = (tree->contmod->mu_down[start+tree->n_root->num]/tree->contmod->var_down[start+tree->n_root->num] +
              root_mean / root_var)*var;

      tree->n_root->ldsk->coord->lonlat[i] =
        Rnorm(tree->contmod->mu_down[start+tree->n_root->num],
              sqrt(tree->contmod->var_down[start+tree->n_root->num]));
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl RRW_Integrated_Lk(t_tree *tree)
{
  phydbl lnL,root_mean,root_var;
  int i,err,start;

  RRW_Update_Normalization_Factor(tree);

  lnL = 0.0;

  Init_Contmod_Locations(tree);

  RRW_Integrated_Lk_Location_Post(NULL, tree->n_root, tree);

  if (tree->contmod->both_sides[LOCATION] == YES)
  {
    RRW_Integrated_Lk_Location_Pre(tree->n_root, tree->n_root->v[1], tree);
    RRW_Integrated_Lk_Location_Pre(tree->n_root, tree->n_root->v[2], tree);
  }

  for (i = 0; i < tree->mmod->n_dim; ++i)
  {
    root_var = tree->mmod->rw_root_var;
    root_mean = LOCATION_Mean_Lonlat(i,tree);
    // root_mean = tree->mmod->rw_root_mean;

    start = Contmod_Start(LOCATION, i, tree);

    tree->contmod->lnL[LOCATION * tree->mmod->n_dim + i] = 0.0;
    tree->contmod->lnL[LOCATION * tree->mmod->n_dim + i] += tree->contmod->logrem_down[start + tree->n_root->num];
    tree->contmod->lnL[LOCATION * tree->mmod->n_dim + i] += Log_Dnorm(tree->contmod->mu_down[start + tree->n_root->num],
                                                                      root_mean,
                                                                      sqrt(root_var + tree->contmod->var_down[start + tree->n_root->num]),
                                                                      &err);

    lnL += tree->contmod->lnL[LOCATION * tree->mmod->n_dim + i];
  }
  return(lnL);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl RRW_Location_Variance_Along_Edge(t_node *d, short int dim, t_tree *tree)
{
  phydbl logvar;

  assert(d != tree->n_root);
  
  logvar =
    log(tree->mmod->sigsq[dim]) +
    log(tree->mmod->sigsq_scale[d->num]) + 
    log(tree->mmod->sigsq_scale_norm_fact) +
    log(fabs(tree->times->nd_t[d->num]-tree->times->nd_t[d->anc->num]));
  
  return(exp(logvar));
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// See Pybus et al. 10.1073/pnas.1206598109 + see my technical notes
void RRW_Integrated_Lk_Location_Post(t_node *a, t_node *d, t_tree *tree)
{
  if(d->tax == TRUE)return;
  else
    {
      int i;

      for (i = 0; i < 3; ++i)
      {
        if (d->v[i] != a && !(a == tree->n_root && d->b[i] == tree->e_root))
        {
          RRW_Integrated_Lk_Location_Post(d, d->v[i], tree);
        }
      }
      RRW_Update_Lk_Location_Down(a, d, tree);
    }
    return;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RRW_Integrated_Lk_Location_Pre(t_node *a, t_node *d, t_tree *tree)
{
  int i;

  if (d == tree->n_root)
  {
    RRW_Integrated_Lk_Location_Pre(tree->n_root, tree->n_root->v[1], tree);
    RRW_Integrated_Lk_Location_Pre(tree->n_root, tree->n_root->v[2], tree);
    return;
  }

  RRW_Update_Lk_Location_Up(a, d, tree);

  if (d->tax == TRUE)
    return;
  else
  {
    for (i = 0; i < 3; ++i)
    {
      if (d->v[i] != a && !(a == tree->n_root && d->b[i] == tree->e_root))
      {
        RRW_Integrated_Lk_Location_Pre(d, d->v[i], tree);
      }
    }
  }
  return;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RRW_Update_Lk_Location_Down(t_node *a, t_node *d, t_tree *tree)
{
  int i,start;
  t_node *son, *bro, *dad;
  phydbl son_mu_down,bro_mu_down;
  phydbl son_var_down,bro_var_down;
  phydbl son_var,bro_var;
  phydbl son_logrem_down,bro_logrem_down;
  phydbl son_a, bro_a;
  phydbl son_b, bro_b;
  
  if(d->tax) return;

  dad = d;
  son = bro = NULL;
  for(i=0;i<3;++i)
    {
      if(d->v[i] != a && !(a == tree->n_root && d->b[i] == tree->e_root))
        {
          if(son == NULL) son = d->v[i];
          else bro = d->v[i];
        }
    }
  
  for(i=0;i<tree->mmod->n_dim;++i)
    {
      start = Contmod_Start(LOCATION, i, tree);

      son_mu_down = tree->contmod->mu_down[start + son->num];
      bro_mu_down = tree->contmod->mu_down[start + bro->num];

      son_var_down = tree->contmod->var_down[start + son->num];
      bro_var_down = tree->contmod->var_down[start + bro->num];

      son_logrem_down = tree->contmod->logrem_down[start + son->num];
      bro_logrem_down = tree->contmod->logrem_down[start + bro->num];

      son_var = RRW_Location_Variance_Along_Edge(son, i, tree);
      bro_var = RRW_Location_Variance_Along_Edge(bro, i, tree);

      son_a = 1.0;
      bro_a = 1.0;

      son_b = 0.0;
      bro_b = 0.0;

      RW_Integrated_Lk_Down(son_a, son_b, son_mu_down, son_var_down, son_var,
                            bro_a, bro_b, bro_mu_down, bro_var_down, bro_var,
                            son_logrem_down, bro_logrem_down,
                            tree->contmod->mu_down + start + dad->num,
                            tree->contmod->var_down + start + dad->num,
                            tree->contmod->logrem_down + start + dad->num);

      assert(isnan(tree->contmod->mu_down[start + dad->num]) == NO);
      assert(isnan(tree->contmod->var_down[start + dad->num]) == NO && !(tree->contmod->var_down[start + dad->num] < 0.0));
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RRW_Update_Lk_Location_Up(t_node *a, t_node *d, t_tree *tree)
{
  int i,start;
  t_node *dad, *son, *bro;
  phydbl dad_mu_up, dad_var_up, dad_logrem_up;
  phydbl son_a, son_b, son_var;
  phydbl bro_a, bro_b, bro_mu_down, bro_var_down, bro_var, bro_logrem_down;
    
  dad = a;
  son = d;
  bro = NULL;

  for(i=0;i<3;++i)
    if(a->v[i] != d && a->v[i] != a->anc && a->b[i] != tree->e_root)
      {
        bro = a->v[i];
        break;
      }

  assert(bro->anc == dad);
  assert(bro != son);

  for(i=0;i<tree->mmod->n_dim;++i)
    {
      start = Contmod_Start(LOCATION,i,tree);

      bro_mu_down     = tree->contmod->mu_down[start+bro->num];
      bro_var_down    = tree->contmod->var_down[start+bro->num];
      bro_logrem_down = tree->contmod->logrem_down[start+bro->num];

      bro_var         = RRW_Location_Variance_Along_Edge(bro,i,tree);    
      bro_a           = 1.0;
      bro_b           = 0.0;

      son_var         = RRW_Location_Variance_Along_Edge(son,i,tree);
      son_a           = 1.0;
      son_b           = 0.0;
      
      if(dad != tree->n_root)
        {
          dad_mu_up     = tree->contmod->mu_up[start+dad->num];
          dad_var_up    = tree->contmod->var_up[start+dad->num];
          dad_logrem_up = tree->contmod->logrem_up[start+dad->num];
        }
      else
        {
          dad_mu_up     = LOCATION_Mean_Lonlat(i,tree);
          // dad_mu_up     = tree->mmod->rw_root_mean;
          dad_var_up    = tree->mmod->rw_root_var;
          dad_logrem_up = 0.0;
        }

        RW_Integrated_Lk_Up(dad_mu_up, dad_var_up, dad_logrem_up,
                            son_a, son_b, son_var,
                            bro_a, bro_b, bro_mu_down, bro_var_down, bro_var, bro_logrem_down,
                            tree->contmod->mu_up + start + son->num,
                            tree->contmod->var_up + start + son->num,
                            tree->contmod->logrem_up + start + son->num);

        assert(isnan(tree->contmod->mu_up[start + d->num]) == NO);
        assert(isnan(tree->contmod->var_up[start + d->num]) == NO && !(tree->contmod->var_down[start + d->num] < 0.0));

      
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
/* Compute the velocity at each tip under the RRW model. For each dimension, compute the slope of the line through the tip and the ancestor. Divide by the time difference
   to get the velocity */
void RRW_Tip_Velocities(t_tree *tree)
{
  int i,j;
  short int dist_type;
  t_node *n;
  phydbl dt;

  dist_type = tree->mmod->dist_type;

  RRW_Sample_Node_Locations_Joint(tree);

  for (i = 0; i < tree->n_otu; ++i)
  {
    n = tree->a_nodes[i];
    dt = tree->times->nd_t[n->num] - tree->times->nd_t[n->anc->num];

    for (j = 0; j < tree->mmod->n_dim; ++j)
    {
      n->ldsk->veloc->deriv[j] = n->ldsk->coord->lonlat[j] - n->anc->ldsk->coord->lonlat[j];
      n->ldsk->veloc->deriv[j] /= dt;
    }
  }
}
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
