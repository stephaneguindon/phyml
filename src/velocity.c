/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/



#include "velocity.h"

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl VELOC_Wrap_Lk(t_edge *b, t_tree *tree, supert_tree *stree)
{
  return(VELOC_Lk(tree));  
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

short int VELOC_Is_Integrated_Velocity(t_phyrex_mod *mod)
{
  return(IWN_Is_Iwn(mod) || IBM_Is_Ibm(mod) || IOU_Is_Iou(mod));
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl VELOC_Lk(t_tree *tree)
{
  phydbl lnL_loc, lnL_veloc;

  RRW_Update_Normalization_Factor(tree);

  lnL_loc   = UNLIKELY;
  lnL_veloc = UNLIKELY;
  
  lnL_loc = (tree->mmod->integrateAncestralLocations == NO) ? (VELOC_Augmented_Lk_Locations(tree)) : (VELOC_Integrated_Lk_Location(tree));
  lnL_veloc = VELOC_Augmented_Lk_Velocities(tree);

  
  if(tree->mmod->print_lk == YES)
    {
      PhyML_Printf("\n. Sigsq: %f %f",tree->mmod->sigsq[0],tree->mmod->sigsq[1]);
      PhyML_Printf("\n. IBM locations: %f",lnL_loc);
      PhyML_Printf("\n. IBM velocities: %f",lnL_veloc);
      PhyML_Printf("\n. Sum: %f",lnL_loc + lnL_veloc);
    }

  return(lnL_loc + lnL_veloc);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* Returns log[Pr(locations at all nodes | velocities at all nodes, tree)] */
phydbl VELOC_Augmented_Lk_Locations(t_tree *tree)
{
  phydbl lnP,disk_lnP;
  t_dsk *disk;
  
  disk_lnP = 0.0;
  lnP = 0.0;
  disk = tree->young_disk;
  
  do
    {
      disk_lnP = VELOC_Augmented_Lk_Locations_Core(disk,tree);
      lnP += disk_lnP;
      disk = disk->prev;
    }
  while(disk);
  
  return(lnP);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl VELOC_Augmented_Lk_Locations_Core(t_dsk *disk, t_tree *tree)
{
  phydbl lnP;
  int i,err;
  phydbl root_mean,root_var;

  root_mean = 0.0;
  root_var = 100.;
  
  assert(disk);
  
  if(disk == tree->young_disk) return 0.0;
  if(disk->age_fixed == YES) return 0.0;

  lnP = 0.0;

  if(disk->ldsk != NULL)
    {
      for(i=0;i<disk->ldsk->n_next;++i)
        {          
          lnP += VELOC_Locations_Forward_Lk_Path(disk->ldsk,disk->ldsk->next[i],tree);
        }

      if(disk->ldsk->prev == NULL) /* root node */
        {
          for(i=0;i<tree->mmod->n_dim;++i)
            {
              root_mean = LOCATION_Mean_Lonlat(i,tree);
              lnP += Log_Dnorm(disk->ldsk->coord->lonlat[i],root_mean,sqrt(root_var),&err); /* Marginal density of location at root */
            }
        }
    }

  return(lnP);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl VELOC_Locations_Forward_Lk_Path(t_ldsk *a, t_ldsk *d, t_tree *tree)
{
  t_ldsk *ldsk;
  phydbl lnP,ld,disk_lnP,var,dt,mean;
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
  mean = var = -1.;
  dt = -1.;
  
  do
    {
      assert(ldsk->prev);

      disk_lnP = 0.0;
      
      for(i=0;i<tree->mmod->n_dim;++i)
        {
          /* dt = fabs(ldsk->disk->time-ldsk->prev->disk->time); */
          
          /* mean = ldsk->prev->coord->lonlat[i] + dt*ldsk->prev->veloc->deriv[i]; */

          /* sd = */
          /*   log(tree->mmod->sigsq[i]) + */
          /*   log(tree->mmod->sigsq_scale[nd_d->num]) + */
          /*   log(tree->mmod->sigsq_scale_norm_fact) + */
          /*   3.*log(dt)- */
          /*   LOG3; */


          mean = VELOC_Location_Mean_Along_Edge(nd_d,i,tree);          
          var  = VELOC_Location_Variance_Along_Edge(nd_d,i,tree);
          
          ld = ldsk->coord->lonlat[i];

          disk_lnP += Log_Dnorm(ld,mean,sqrt(var),&err);

          if(tree->mmod->print_lk == YES)
            PhyML_Printf("\n. LOC %2d Time: %10f nd_d: %d sigsq: (%f,%f,%f) a->coord: %10f a->veloc: %10f mean: %10f var: %10f dt: %10f [%10f %10f] x: %10f lk: %10f",
                         i,
                         ldsk->disk->time,
                         nd_d->num,
                         tree->mmod->sigsq[i],
                         tree->mmod->sigsq_scale[nd_d->num],
                         tree->mmod->sigsq_scale_norm_fact,
                         ldsk->prev->coord->lonlat[i],
                         ldsk->prev->veloc->deriv[i],
                         mean,
                         var,
                         dt,
                         ldsk->disk->time,ldsk->prev->disk->time,
                         ld,
                         disk_lnP);
          
          
          if(isinf(lnP) || isnan(lnP)) return(UNLIKELY);
          
          if(isnan(lnP))
            {
              PhyML_Printf("\n. ld=%f var=%f dt=[%f,%f] sigsq=%f",ld,var,ldsk->disk->time,ldsk->prev->disk->time,tree->mmod->sigsq);
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

/* Returns log[Pr(velocities at all nodes | tree)] */
phydbl VELOC_Augmented_Lk_Velocities(t_tree *tree)
{
  phydbl lnP,disk_lnP;
  t_dsk *disk;

  disk_lnP = 0.0;
  lnP = 0.0;
  disk = tree->young_disk;
  
  do
    {
      disk_lnP = VELOC_Augmented_Lk_Velocities_Core(disk,tree);
      lnP += disk_lnP;      
      disk = disk->prev;
    }  
  while(disk);
  
  return(lnP);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl VELOC_Augmented_Lk_Velocities_Core(t_dsk *disk, t_tree *tree)
{
  int i,err;
  phydbl lnP,root_mean,root_var;

  root_mean = 0.0;
  root_var = 0.0;

  assert(disk);
  
  if(disk == tree->young_disk) return 0.0;
  if(disk->age_fixed == YES) return 0.0;

  lnP = 0.0;

  if(disk->ldsk != NULL)
    {
      for(i=0;i<disk->ldsk->n_next;++i)
        {          
          lnP += VELOC_Velocities_Forward_Lk_Path(disk->ldsk,disk->ldsk->next[i],tree);
        }

      if(disk->ldsk->prev == NULL) /* root node */
        {
          for(i=0;i<tree->mmod->n_dim;++i)
            {
              root_mean = 0.0;
              root_var  = 10.;
              lnP += Log_Dnorm(disk->ldsk->veloc->deriv[i],root_mean,sqrt(root_var),&err); /* Marginal density of velocity at root */
            }
        }
    }

  return(lnP);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl VELOC_Velocities_Forward_Lk_Path(t_ldsk *a, t_ldsk *d, t_tree *tree)
{
  t_ldsk *ldsk;
  phydbl lnP,ld,disk_lnP,var,dt,mean;
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
  mean = var = -1.;
  dt = -1.;
  
  do
    {
      assert(ldsk->prev);

      disk_lnP = 0.0;
      
      for(i=0;i<tree->mmod->n_dim;++i)
        {
          dt = fabs(ldsk->disk->time-ldsk->prev->disk->time);

          /* mean = 0.5*(3.*(ldsk->coord->lonlat[i] - ldsk->prev->coord->lonlat[i])/dt - ldsk->prev->veloc->deriv[i]); */

          /* sd = */
          /*   log(tree->mmod->sigsq[i]) + */
          /*   log(tree->mmod->sigsq_scale[nd_d->num]) + */
          /*   log(tree->mmod->sigsq_scale_norm_fact) + */
          /*   log(dt)- */
          /*   LOG4; */

          mean = VELOC_Velocity_Mean_Along_Edge(nd_d,i,tree);          
          var = VELOC_Velocity_Variance_Along_Edge(nd_d,i,tree);
          
          ld = ldsk->veloc->deriv[i];          
          
          disk_lnP += Log_Dnorm(ld,mean,sqrt(var),&err);
          
          if(tree->mmod->print_lk == YES)
            PhyML_Printf("\n. VEL %2d Time: %10f d->coord: %10f a->coord: %10f a->veloc: %10f d->veloc: %10f mean: %10f var: %10f dt: %10f [%10f; %10f] x: %10f lk: %10f",
                         i,
                         ldsk->disk->time,
                         ldsk->coord->lonlat[i],
                         ldsk->prev->coord->lonlat[i],
                         ldsk->prev->veloc->deriv[i],
                         ldsk->veloc->deriv[i],
                         mean,
                         var,
                         dt,
                         ldsk->disk->time,ldsk->prev->disk->time,
                         ld,
                         disk_lnP);

          if(isinf(lnP) || isnan(lnP)) return(UNLIKELY);
          
          if(isnan(lnP))
            {
              PhyML_Printf("\n. ld=%f var=%f dt=[%f,%f] sigsq=%f",ld,var,ldsk->disk->time,ldsk->prev->disk->time,tree->mmod->sigsq);
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

phydbl VELOC_Integrated_Lk_Location(t_tree *tree)
{
  phydbl lnL,root_mean,root_var;
  int i,err;
  
  RRW_Update_Normalization_Factor(tree);
  
  lnL = 0.0;
  /* root_var = 1.0; */
  root_var = 100.;
  root_mean = 0.0;
  
  for(i=0;i<tree->mmod->n_dim;++i)
    {
      VELOC_Init_Contmod_Locations(i,tree);
      VELOC_Integrated_Lk_Location_Post(NULL,tree->n_root,i,tree,NO);

      /* root_mean = 0.0; */
      root_mean = LOCATION_Mean_Lonlat(i,tree);
      
      lnL += tree->contmod->logrem_down[tree->n_root->num];
      lnL += Log_Dnorm(tree->contmod->mu_down[tree->n_root->num],root_mean,sqrt(root_var+tree->contmod->var_down[tree->n_root->num]),&err);
      /* PhyML_Printf("\n. root mean: %f var: %f", */
      /*              root_mean, */
      /*              root_var+tree->contmod->var_down[tree->n_root->num]); */
    }
  return(lnL);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// See Pybus et al. 10.1073/pnas.1206598109 + see my technical notes
void VELOC_Integrated_Lk_Location_Post(t_node *a, t_node *d, short int dim, t_tree *tree, short int print)
{
  if(d->tax == TRUE)
    {
      return;
    }
  else
    {
      int i;
      t_node *v1, *v2;
      phydbl v1mu,v2mu;
      phydbl v1var,v2var;
      phydbl dv1var,dv2var;
      phydbl v1logrem,v2logrem;
      phydbl dtv1, dtv2;
      phydbl av1, av2;
      phydbl bv1, bv2;
      
      for(i=0;i<3;++i)
        {
          if(d->v[i] != a && !(a == tree->n_root && d->b[i] == tree->e_root))
            {
              VELOC_Integrated_Lk_Location_Post(d,d->v[i],dim,tree,print);
            }
        }

      v1 = v2 = NULL;
      for(i=0;i<3;++i)
        {
          if(d->v[i] != a && !(a == tree->n_root && d->b[i] == tree->e_root))
            {
              if(v1 == NULL) v1 = d->v[i];
              else v2 = d->v[i];
            }
        }
      
      v1mu = tree->contmod->mu_down[v1->num];
      v2mu = tree->contmod->mu_down[v2->num];

      v1var = tree->contmod->var_down[v1->num];
      v2var = tree->contmod->var_down[v2->num];

      v1logrem = tree->contmod->logrem_down[v1->num];
      v2logrem = tree->contmod->logrem_down[v2->num];

      dv1var = VELOC_Location_Variance_Along_Edge(v1,dim,tree);
      dv2var = VELOC_Location_Variance_Along_Edge(v2,dim,tree);

      dtv1 = fabs(tree->times->nd_t[v1->num] - tree->times->nd_t[d->num]);
      dtv2 = fabs(tree->times->nd_t[v2->num] - tree->times->nd_t[d->num]);
      
      av1 = 1.0;
      av2 = 1.0;

      /* Below is only correct for IBM... !!!!!!!!!!!!!!!!!!! */
      /* bv1 = (v1->ldsk->veloc->deriv[dim] + d->ldsk->veloc->deriv[dim])*.5*dtv1; */
      /* bv2 = (v2->ldsk->veloc->deriv[dim] + d->ldsk->veloc->deriv[dim])*.5*dtv2; */

      bv1 = VELOC_Location_Mean_Along_Edge(v1,dim,tree) - d->ldsk->coord->lonlat[dim];
      bv2 = VELOC_Location_Mean_Along_Edge(v2,dim,tree) - d->ldsk->coord->lonlat[dim];

      if(d == tree->n_root && print == YES)
        {
          PhyML_Printf("\n. v1mu=%f v2mu=%f v1var=%f dv1var=%f v2var=%f dv2var=%f t=%f t1=%f t2=%f",
                       v1mu,
                       v2mu,
                       v1var,dv1var,
                       v2var,dv2var,
                       tree->times->nd_t[d->num],
                       tree->times->nd_t[v1->num],
                       tree->times->nd_t[v2->num]);
        }
      
                         
      if(IBM_Is_Ibm(tree->mmod) == YES)
        {
          IBM_Integrated_Location_Down(dtv1,dtv2,
                                       av1,bv1,v1mu,v1var,dv1var,
                                       av2,bv2,v2mu,v2var,dv2var,
                                       v1logrem,v2logrem,
                                       tree->contmod->mu_down+d->num,
                                       tree->contmod->var_down+d->num,
                                       tree->contmod->logrem_down+d->num);
        }
      else if(IWN_Is_Iwn(tree->mmod) == YES)
        {
          IWN_Integrated_Location_Down(dtv1,dtv2,
                                       av1,bv1,v1mu,v1var,dv1var,
                                       av2,bv2,v2mu,v2var,dv2var,
                                       v1logrem,v2logrem,
                                       d->ldsk->veloc->deriv[dim],v1->ldsk->veloc->deriv[dim],v2->ldsk->veloc->deriv[dim],
                                       tree->mmod->omega,
                                       tree->contmod->mu_down+d->num,
                                       tree->contmod->var_down+d->num,
                                       tree->contmod->logrem_down+d->num);
        }
      else if(IOU_Is_Iou(tree->mmod) == YES)
        {
          IOU_Integrated_Location_Down(dtv1,dtv2,
                                       av1,bv1,v1mu,v1var,dv1var,
                                       av2,bv2,v2mu,v2var,dv2var,
                                       v1logrem,v2logrem,
                                       tree->contmod->mu_down+d->num,
                                       tree->contmod->var_down+d->num,
                                       tree->contmod->logrem_down+d->num);
        }

      /* if(v1->tax && v2->tax) */
        /* { */
        /*   PhyML_Printf("\n. LOCATION  %c %d (%d,%d) loc: v1:%f v2:%f -- mean: %f sd: %f dt1: %f dt2: %f sigsq: %f derivatives: %f %f %f av1: %f bv1: %f v1mu: %f av2: %f bv2: %f v2mu: %f ", */
        /*                d == tree->n_root ? '*' : ' ', */
        /*                d->num, */
        /*                v1->num,v2->num, */
        /*                v1->ldsk->coord->lonlat[dim],v2->ldsk->coord->lonlat[dim], */
        /*                tree->contmod->mu_down[d->num],sqrt(tree->contmod->var_down[d->num]), */
        /*                dtv1,dtv2, */
        /*                tree->mmod->sigsq[dim], */
        /*                v1->ldsk->veloc->deriv[dim], */
        /*                v2->ldsk->veloc->deriv[dim], */
        /*                d->ldsk->veloc->deriv[dim], */
        /*                av1,bv1,v1mu, */
        /*                av2,bv2,v2mu); */
        /* } */

    }
  
  return;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void VELOC_Integrated_Lk_Location_Pre(t_node *a, t_node *d, short int dim, t_tree *tree)
{
  int i;
  t_node *v1, *v2;
  phydbl v1mu,v2mu;
  phydbl v1var,v2var;
  phydbl av1var,av2var;
  phydbl v1logrem,v2logrem;
  phydbl dtv1, dtv2;
  phydbl av1, av2;
  phydbl bv1, bv2;
  
  
  v1 = NULL;
  if(a != tree->n_root)
    {
      v1 = a->anc;
      assert(v1);
    }
  
  v2 = NULL;
  for(i=0;i<3;++i)
    {
      if(a->v[i] != d && a->v[i] != v1 && !(v1 == tree->n_root && a->b[i] == tree->e_root))
        {
          v2 = a->v[i];
          break;
        }
    }

  av1 = bv1 = -1.;
  v1mu = v1var = v1logrem = -1.;
  dtv1 = -1.;
  av1var = -1.;
  
  assert(v2->anc == a);
  assert(v2->anc->anc == v1);
    
  v2mu     = tree->contmod->mu_down[v2->num];
  v2var    = tree->contmod->var_down[v2->num];
  v2logrem = tree->contmod->logrem_down[v2->num];
  dtv2     = fabs(tree->times->nd_t[v2->num] - tree->times->nd_t[a->num]);      
  av2var   = VELOC_Location_Variance_Along_Edge(v2,dim,tree);    

  av2      = 1.0;
  /* Below is only correct for IBM... !!!!!!!!!!!!!!!!!!! */
  /* bv2      = (v2->ldsk->veloc->deriv[dim] + a->ldsk->veloc->deriv[dim])*.5*dtv2; */
  bv2      = VELOC_Location_Mean_Along_Edge(v2,dim,tree) - a->ldsk->coord->lonlat[dim];
      
  if(v1 != NULL)
    {
      v1mu     = tree->contmod->mu_up[a->num];
      v1var    = tree->contmod->var_up[a->num];
      v1logrem = tree->contmod->logrem_up[a->num];
      dtv1     = fabs(tree->times->nd_t[a->num] - tree->times->nd_t[v1->num]);          
      av1var   = VELOC_Location_Variance_Along_Edge(a,dim,tree);

      av1      = 1.0;
      /* Below is only correct for IBM... !!!!!!!!!!!!!!!!!!! */
      /* bv1      = (a->ldsk->veloc->deriv[dim] + v1->ldsk->veloc->deriv[dim])*.5*dtv1; */
      bv1      = VELOC_Location_Mean_Along_Edge(a,dim,tree) - v1->ldsk->coord->lonlat[dim];
  

      if(IBM_Is_Ibm(tree->mmod) == YES)
        {
          IBM_Integrated_Location_Up(dtv1,dtv2,
                                     av1,bv1,v1mu,v1var,av1var,
                                     av2,bv2,v2mu,v2var,av2var,
                                     v1logrem,v2logrem,
                                     tree->contmod->mu_up+d->num,
                                     tree->contmod->var_up+d->num,
                                     tree->contmod->logrem_up+d->num,
                                     NO);
        }
      else if(IWN_Is_Iwn(tree->mmod) == YES)
        {
          IWN_Integrated_Location_Up(dtv1,dtv2,
                                     av1,bv1,v1mu,v1var,av1var,
                                     av2,bv2,v2mu,v2var,av2var,
                                     v1logrem,v2logrem,
                                     a->ldsk->veloc->deriv[dim],v1->ldsk->veloc->deriv[dim],v2->ldsk->veloc->deriv[dim],
                                     tree->mmod->omega,
                                     tree->contmod->mu_up+d->num,
                                     tree->contmod->var_up+d->num,
                                     tree->contmod->logrem_up+d->num,
                                     NO);
        }
      else if(IOU_Is_Iou(tree->mmod) == YES)
        {
          IOU_Integrated_Location_Up(dtv1,dtv2,
                                     av1,bv1,v1mu,v1var,av1var,
                                     av2,bv2,v2mu,v2var,av2var,
                                     v1logrem,v2logrem,
                                     tree->contmod->mu_up+d->num,
                                     tree->contmod->var_up+d->num,
                                     tree->contmod->logrem_up+d->num,
                                     NO);
        }

    }
  else
    {
      if(IBM_Is_Ibm(tree->mmod) == YES)
        {
          IBM_Integrated_Location_Up(dtv1,dtv2,
                                     av1,bv1,v1mu,v1var,av1var,
                                     av2,bv2,v2mu,v2var,av2var,
                                     v1logrem,v2logrem,
                                     tree->contmod->mu_up+d->num,
                                     tree->contmod->var_up+d->num,
                                     tree->contmod->logrem_up+d->num,
                                     YES);
        }
      else if(IWN_Is_Iwn(tree->mmod) == YES)
        {
          IWN_Integrated_Location_Up(dtv1,dtv2,
                                     av1,bv1,v1mu,v1var,av1var,
                                     av2,bv2,v2mu,v2var,av2var,
                                     v1logrem,v2logrem,
                                     a->ldsk->veloc->deriv[dim],-INFINITY,v2->ldsk->veloc->deriv[dim],
                                     tree->mmod->omega,
                                     tree->contmod->mu_up+d->num,
                                     tree->contmod->var_up+d->num,
                                     tree->contmod->logrem_up+d->num,
                                     YES);
        }
      else if(IOU_Is_Iou(tree->mmod) == YES)
        {
          IOU_Integrated_Location_Up(dtv1,dtv2,
                                     av1,bv1,v1mu,v1var,av1var,
                                     av2,bv2,v2mu,v2var,av2var,
                                     v1logrem,v2logrem,
                                     tree->contmod->mu_up+d->num,
                                     tree->contmod->var_up+d->num,
                                     tree->contmod->logrem_up+d->num,
                                     YES);
        }
    }
  
  if(d->tax == TRUE) return;
  else
    {
      for(i=0;i<3;++i)
        {
          if(d->v[i] != a && !(a == tree->n_root && d->b[i] == tree->e_root))
            {
              VELOC_Integrated_Lk_Location_Pre(d,d->v[i],dim,tree);
            }
        }
    }
  return;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl VELOC_Velocity_Variance_Along_Edge(t_node *d, short int dim, t_tree *tree)
{
  if(IBM_Is_Ibm(tree->mmod) == YES)
    {
      return(IBM_Velocity_Variance_Along_Edge(d,dim,tree));
    }
  else if(IWN_Is_Iwn(tree->mmod) == YES)
    {
      return(IWN_Velocity_Variance_Along_Edge(d,dim,tree));
    }
  else if(IOU_Is_Iou(tree->mmod) == YES)
    {
      return(IOU_Velocity_Variance_Along_Edge(d,dim,tree));
    }
  else assert(false);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl VELOC_Location_Variance_Along_Edge(t_node *d, short int dim, t_tree *tree)
{
  if(IBM_Is_Ibm(tree->mmod) == YES)
    {
      return(IBM_Location_Variance_Along_Edge(d,dim,tree));
    }
  else if(IWN_Is_Iwn(tree->mmod) == YES)
    {
      return(IWN_Location_Variance_Along_Edge(d,dim,tree));
    }
  else if(IOU_Is_Iou(tree->mmod) == YES)
    {
      return(IOU_Location_Variance_Along_Edge(d,dim,tree));
    }
  else assert(false);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl VELOC_Velocity_Mean_Along_Edge(t_node *d, short int dim, t_tree *tree)
{
  if(IBM_Is_Ibm(tree->mmod) == YES)
    {
      return(IBM_Velocity_Mean_Along_Edge(d,dim,tree));
    }
  else if(IWN_Is_Iwn(tree->mmod) == YES)
    {
      return(IWN_Velocity_Mean_Along_Edge(d,dim,tree));
    }
  else if(IOU_Is_Iou(tree->mmod) == YES)
    {
      return(IOU_Velocity_Mean_Along_Edge(d,dim,tree));
    }
  else assert(false);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl VELOC_Location_Mean_Along_Edge(t_node *d, short int dim, t_tree *tree)
{
  if(IBM_Is_Ibm(tree->mmod) == YES)
    {
      return(IBM_Location_Mean_Along_Edge(d,dim,tree));
    }
  else if(IWN_Is_Iwn(tree->mmod) == YES)
    {
      return(IWN_Location_Mean_Along_Edge(d,dim,tree));
    }
  else if(IOU_Is_Iou(tree->mmod) == YES)
    {
      return(IOU_Location_Mean_Along_Edge(d,dim,tree));
    }
  else assert(false);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void VELOC_Sample_Node_Locations(t_tree *tree)
{
  int i,j;
  t_node *n;
  phydbl au,bu,varu,var,mean;

  n = NULL;
  au = bu = varu = var = mean = -1;
  
  RRW_Update_Normalization_Factor(tree);

  for(i=0;i<tree->mmod->n_dim;++i)
    {
      VELOC_Init_Contmod_Locations(i,tree);
      VELOC_Integrated_Lk_Location_Post(NULL,tree->n_root,i,tree,NO);
      VELOC_Integrated_Lk_Location_Pre(tree->n_root,tree->n_root->v[1],i,tree);
      VELOC_Integrated_Lk_Location_Pre(tree->n_root,tree->n_root->v[2],i,tree);

      for(j=0;j<2*tree->n_otu-2;++j)
        {
          if(tree->a_nodes[j]->tax == NO && tree->a_nodes[j] != tree->n_root)
            {
              n = tree->a_nodes[j];

              au = 1.0;
              /* Below is only correct for IBM model !!!!!!!!!!!!!!!!!!!!!! */
              /* bu = (n->ldsk->veloc->deriv[i] + n->anc->ldsk->veloc->deriv[i])*.5*dt; */
              bu   = VELOC_Location_Mean_Along_Edge(n,i,tree) - n->anc->ldsk->coord->lonlat[i];
              varu = VELOC_Location_Variance_Along_Edge(n,i,tree); 


              if(tree->contmod->var_down[n->num] > SMALL && pow(au,2)*tree->contmod->var_up[n->num]+varu > SMALL)
                {
                  var = 1./tree->contmod->var_down[n->num] + 1. / (pow(au,2)*tree->contmod->var_up[n->num]+varu);
                  var = 1./var;
                  mean = (tree->contmod->mu_down[n->num]/tree->contmod->var_down[n->num] + (au*tree->contmod->mu_up[n->num]+bu)/(pow(au,2)*tree->contmod->var_up[n->num]+varu))*var;
                }
              else if(tree->contmod->var_down[n->num] > SMALL)
                {
                  var = 1./tree->contmod->var_down[n->num];
                  var = 1./var;
                  mean = tree->contmod->mu_down[n->num];
                }
              else if(pow(au,2)*tree->contmod->var_up[n->num]+varu > SMALL)
                {
                  var = 1. / (pow(au,2)*tree->contmod->var_up[n->num]+varu);
                  var = 1./var;
                  mean = au*tree->contmod->mu_up[n->num]+bu;
                }
              else
                {
                  var = 0.0;
                  mean = tree->contmod->mu_down[n->num];
                }
              
              /* Below is only valid if prior distribution of location at that node is flat */
              tree->a_nodes[j]->ldsk->coord->lonlat[i] = Rnorm(mean,sqrt(var));
            }
        }

      assert(isnan(tree->contmod->var_down[tree->n_root->num]) == NO);
      
      /* PhyML_Printf("\n\n. root mean: %f var: %f", */
      /*              tree->contmod->mu_down[tree->n_root->num], */
      /*              tree->contmod->var_down[tree->n_root->num]); */
      
      /* Below is only valid if prior distribution at root is flat */
      tree->n_root->ldsk->coord->lonlat[i] =
        Rnorm(tree->contmod->mu_down[tree->n_root->num],
              sqrt(tree->contmod->var_down[tree->n_root->num]));
    }
}



//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////




/* DEPRECATED CODE BELOW */



//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////



// See Pybus et al. 10.1073/pnas.1206598109 + see my technical notes
void VELOC_Integrated_Lk_Velocity_Post(t_node *a, t_node *d, short int dim, t_tree *tree, short int print)
{
  if(d->tax == TRUE)
    {
      return;
    }
  else
    {
      int i,err;
      t_node *v1, *v2;
      phydbl v1mu,v2mu;
      phydbl v1var,v2var;
      phydbl dv1var,dv2var;
      phydbl v1logrem,v2logrem;
      phydbl dtv1, dtv2;
      phydbl av1, av2;
      phydbl bv1, bv2;
      
      for(i=0;i<3;++i)
        {
          if(d->v[i] != a && !(a == tree->n_root && d->b[i] == tree->e_root))
            {
              VELOC_Integrated_Lk_Velocity_Post(d,d->v[i],dim,tree,print);
            }
        }

      v1 = v2 = NULL;
      for(i=0;i<3;++i)
        {
          if(d->v[i] != a && !(a == tree->n_root && d->b[i] == tree->e_root))
            {
              if(v1 == NULL) v1 = d->v[i];
              else v2 = d->v[i];
            }
        }

      assert(v1->anc == d);
      assert(v2->anc == d);
      
      err = -1;
      
      v1mu = tree->contmod->mu_down[v1->num];
      v2mu = tree->contmod->mu_down[v2->num];

      v1var = tree->contmod->var_down[v1->num];
      v2var = tree->contmod->var_down[v2->num];

      v1logrem = tree->contmod->logrem_down[v1->num];
      v2logrem = tree->contmod->logrem_down[v2->num];

      dv1var = VELOC_Velocity_Variance_Along_Edge(v1,dim,tree);
      dv2var = VELOC_Velocity_Variance_Along_Edge(v2,dim,tree);

      dtv1 = fabs(tree->times->nd_t[v1->num] - tree->times->nd_t[d->num]);
      dtv2 = fabs(tree->times->nd_t[v2->num] - tree->times->nd_t[d->num]);

      av1 = -.5;
      av2 = -.5;

      bv1 = 1.5*(v1->ldsk->coord->lonlat[dim] - d->ldsk->coord->lonlat[dim])/dtv1;
      bv2 = 1.5*(v2->ldsk->coord->lonlat[dim] - d->ldsk->coord->lonlat[dim])/dtv2;

      if(d == tree->n_root && print == YES)
        {
          PhyML_Printf("\n. v1mu=%f v2mu=%f v1var=%f dv1var=%f v2var=%f dv2var=%f t=%f t1=%f t2=%f",
                       v1mu,
                       v2mu,
                       v1var,dv1var,
                       v2var,dv2var,
                       tree->times->nd_t[d->num],
                       tree->times->nd_t[v1->num],
                       tree->times->nd_t[v2->num]);
        }
      
      tree->contmod->var_down[d->num] = pow(av1,2)/(v1var + dv1var) + pow(av2,2)/(v2var + dv2var);
      tree->contmod->var_down[d->num] = 1./tree->contmod->var_down[d->num];

      tree->contmod->mu_down[d->num] = (av1*(v1mu-bv1)/(v1var + dv1var) + av2*(v2mu-bv2)/(v2var + dv2var)) * tree->contmod->var_down[d->num];
      
      tree->contmod->logrem_down[d->num]  = v1logrem + v2logrem;
      tree->contmod->logrem_down[d->num] -= log(fabs(av2*av1));
      tree->contmod->logrem_down[d->num] += Log_Dnorm((v1mu-bv1)/av1,(v2mu-bv2)/av2,sqrt((v1var+dv1var)/pow(av1,2)+(v2var+dv2var)/pow(av2,2)),&err);

      /* if(v1->tax && v2->tax) */
        /* { */
        /*   PhyML_Printf("\n. VELOC %c (%3d,%3d) veloc: %f %f -- mean: %f var: %f logrem: %f dt1: %f dt2: %f sigsq: %f locations: %f %f %f v1mu:%f v2mu:%f bv1:%f bv2:%f av1:%f av2:%f", */
        /*                d == tree->n_root ? '*':' ', */
        /*                v1->num,v2->num, */
        /*                v1->ldsk->veloc->deriv[dim],v2->ldsk->veloc->deriv[dim], */
        /*                tree->contmod->mu_down[d->num], */
        /*                tree->contmod->var_down[d->num], */
        /*                tree->contmod->logrem_down[d->num], */
        /*                dtv1,dtv2, */
        /*                tree->mmod->sigsq[dim], */
        /*                v1->ldsk->coord->lonlat[dim], */
        /*                v2->ldsk->coord->lonlat[dim], */
        /*                d->ldsk->coord->lonlat[dim], */
        /*                v1mu,v2mu, */
        /*                bv1,bv2, */
        /*                av1,av2); */
        /* } */

    }
  
  return;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void VELOC_Integrated_Lk_Velocity_Pre(t_node *a, t_node *d, short int dim, t_tree *tree)
{
  int i,err;
  t_node *v1, *v2;
  phydbl v1mu,v2mu;
  phydbl v1var,v2var;
  phydbl av1var,av2var;
  phydbl v1logrem,v2logrem;
  phydbl av1, av2;
  phydbl bv1, bv2;
  phydbl dtv1, dtv2;
  
  
  v1 = NULL;
  if(a != tree->n_root)
    {
      v1 = a->anc;
      assert(v1);
    }
  
  v2 = NULL;
  for(i=0;i<3;++i)
    {
      if(a->v[i] != d && a->v[i] != v1 && !(v1 == tree->n_root && a->b[i] == tree->e_root))
        {
          v2 = a->v[i];
          break;
        }
    }

  assert(v2->anc == a);
  
  err = -1;
  
  v2mu     = tree->contmod->mu_down[v2->num];
  v2var    = tree->contmod->var_down[v2->num];
  v2logrem = tree->contmod->logrem_down[v2->num];
  dtv2     = fabs(tree->times->nd_t[v2->num] - tree->times->nd_t[a->num]);      
  av2var   = VELOC_Velocity_Variance_Along_Edge(v2,dim,tree);
  av2      = -.5;
  bv2      = 1.5*(v2->ldsk->coord->lonlat[dim] - a->ldsk->coord->lonlat[dim])/dtv2;
  
  
  if(v1 != NULL)
    {
      v1mu     = tree->contmod->mu_up[a->num];
      v1var    = tree->contmod->var_up[a->num];
      v1logrem = tree->contmod->logrem_up[a->num];
      dtv1     = fabs(tree->times->nd_t[a->num] - tree->times->nd_t[v1->num]);
      av1var   = VELOC_Velocity_Variance_Along_Edge(a,dim,tree);
      av1      = -.5;
      bv1      = 1.5*(a->ldsk->coord->lonlat[dim] - v1->ldsk->coord->lonlat[dim])/dtv1;

      
      tree->contmod->var_up[d->num] = pow(av2,2)/(v2var + av2var) + 1./(pow(av1,2)*v1var+av1var);
      tree->contmod->var_up[d->num] = 1./tree->contmod->var_up[d->num];

      tree->contmod->mu_up[d->num] = (fabs(av2)*(v2mu-bv2)/(v2var + av2var) + (fabs(av1)*v1mu+bv1)/(av1*av1*v1var+av1var)) * tree->contmod->var_up[d->num];
      
      tree->contmod->logrem_up[d->num]  = v1logrem + v2logrem;
      tree->contmod->logrem_up[d->num] -= log(fabs(av2));
      tree->contmod->logrem_up[d->num] += Log_Dnorm((v2mu-bv2)/av2,av1*v1mu+bv1,sqrt((v2var+av2var)/pow(av2,2)+pow(av1,2)*v1var+av1var),&err);
    }
  else
    {
      tree->contmod->mu_up[d->num]     = (v2mu-bv2)/av2;
      tree->contmod->var_up[d->num]    = (v2var + av2var)/pow(av2,2);
      tree->contmod->logrem_up[d->num] = v2logrem;
    }
  
  if(d->tax == TRUE) return;
  else
    {
      for(i=0;i<3;++i)
        {
          if(d->v[i] != a && !(a == tree->n_root && d->b[i] == tree->e_root))
            {
              VELOC_Integrated_Lk_Velocity_Pre(d,d->v[i],dim,tree);
            }
        }
    }
  return;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void VELOC_Veloc_Gibbs_Mean_Var(t_node *n, phydbl *mean, phydbl *var, short int dim, t_tree *tree)
{

  if(n->tax == YES)
    {
      phydbl dt;
      
      (*mean) = n->anc->ldsk->veloc->deriv[dim];

      dt = fabs(tree->times->nd_t[n->num] - tree->times->nd_t[n->anc->num]);
      
      (*var) =
        log(tree->mmod->sigsq[dim]) +
        log(tree->mmod->sigsq_scale[n->num]) +
        log(tree->mmod->sigsq_scale_norm_fact) +
        log(dt);
      
      (*var) = exp(*var);      
    }
  else if(n == tree->n_root)
    {
      phydbl var1,var2;
      phydbl dt1,dt2;
      phydbl y1, y2;

      y1 = tree->n_root->v[1]->ldsk->veloc->deriv[dim];
      y2 = tree->n_root->v[2]->ldsk->veloc->deriv[dim];

      dt1 = fabs(tree->times->nd_t[tree->n_root->v[1]->num]-tree->times->nd_t[tree->n_root->num]);
      dt2 = fabs(tree->times->nd_t[tree->n_root->v[2]->num]-tree->times->nd_t[tree->n_root->num]);
      
      var1 =
        log(tree->mmod->sigsq[dim]) +
        log(tree->mmod->sigsq_scale[n->v[1]->num]) +
        log(tree->mmod->sigsq_scale_norm_fact) +
        log(dt1);
      var1 = exp(var1);

      var2 =
        log(tree->mmod->sigsq[dim]) +
        log(tree->mmod->sigsq_scale[n->v[2]->num]) +
        log(tree->mmod->sigsq_scale_norm_fact) +
        log(dt2);
      var2 = exp(var2);


      (*var) = 1./var1 + 1./var2;
      (*var) = 1./(*var);

      (*mean) = (*var) * (y1/var1 + y2/var2);      

      /* PhyML_Printf("\n. y1: %f y2: %f -> mean: %f var: %f", */
      /*              y1,y2, */
      /*              *mean,*var); */

    }
  else
    {
      phydbl var1,var2,var3;
      phydbl dt1,dt2,dt3;
      phydbl y1, y2, y3;
      t_node *v1, *v2, *v3;
      int i;
      
      v3 = n->anc;

      v1 = v2 = NULL;
      for(i=0;i<3;++i)
        {
          if(n->v[i] != n->anc && !(n->anc == tree->n_root && n->b[i] == tree->e_root))
            {
              if(v1 == NULL) v1 = n->v[i];
              else v2 = n->v[i];
            }
        }
      
      assert(v1->anc == n);
      assert(v2->anc == n);
      
      y1 = v1->ldsk->veloc->deriv[dim];
      y2 = v2->ldsk->veloc->deriv[dim];
      y3 = v3->ldsk->veloc->deriv[dim];
      

      dt1 = fabs(tree->times->nd_t[v1->num]-tree->times->nd_t[n->num]);
      dt2 = fabs(tree->times->nd_t[v2->num]-tree->times->nd_t[n->num]);
      dt3 = fabs(tree->times->nd_t[n->num]-tree->times->nd_t[v3->num]);

      var1 =
        log(tree->mmod->sigsq[dim]) +
        log(tree->mmod->sigsq_scale[v1->num]) +
        log(tree->mmod->sigsq_scale_norm_fact) +
        log(dt1);
      var1 = exp(var1);

      var2 =
        log(tree->mmod->sigsq[dim]) +
        log(tree->mmod->sigsq_scale[v2->num]) +
        log(tree->mmod->sigsq_scale_norm_fact) +
        log(dt2);
      var2 = exp(var2);

      var3 =
        log(tree->mmod->sigsq[dim]) +
        log(tree->mmod->sigsq_scale[n->num]) +
        log(tree->mmod->sigsq_scale_norm_fact) +
        log(dt3);
      var3 = exp(var3);


      (*var) = 1./var1 + 1./var2 + 1./var3;
      (*var) = 1./(*var);
      
      (*mean) = (*var) * (y1/var1 + y2/var2 + y3/var3);      


      /* PhyML_Printf("\n. y1: %f y2: %f y3: %f -> mean: %f var: %f", */
      /*              y1,y2,y3, */
      /*              *mean,*var); */
    }
}

