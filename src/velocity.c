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
  return(VELOC_Lk(NULL,tree));  
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

short int VELOC_Is_Integrated_Velocity(t_phyrex_mod *mod)
{
  return(IWN_Is_Iwn(mod) || IBM_Is_Ibm(mod) || IOU_Is_Iou(mod));
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl VELOC_Lk(t_node *z, t_tree *tree)
{
  phydbl lnL_loc, lnL_veloc;

  RRW_Update_Normalization_Factor(tree);


  lnL_loc   = UNLIKELY;
  lnL_veloc = UNLIKELY;

  lnL_loc = (tree->mmod->integrateAncestralLocations == NO) ? (VELOC_Augmented_Lk_Locations(z,tree)) : (VELOC_Integrated_Lk_Location(z,tree));
  lnL_veloc = VELOC_Augmented_Lk_Velocity(z,tree);
    
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
phydbl VELOC_Augmented_Lk_Locations(t_node *z, t_tree *tree)
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
  root_var  = tree->mmod->rw_root_var[LOCATION];
  
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
phydbl VELOC_Augmented_Lk_Velocity(t_node *z, t_tree *tree)
{
  phydbl lnL;
  phydbl root_mean,root_var;
  int i,err,start;

  lnL = UNLIKELY;

  root_var  = tree->mmod->rw_root_var[VELOCITY];
  root_mean = tree->mmod->rw_root_mean[VELOCITY];

  // PhyML_Printf("\n root mean: %f root_var: %f",root_mean,root_var);

  if(z == NULL)
    {
      Init_Contmod_Velocities(tree);

      RRW_Update_Normalization_Factor(tree);
      
      lnL = 0.0;

      VELOC_Augmented_Lk_Velocity_Post(NULL,tree->n_root,tree);

      if(tree->contmod->both_sides[VELOCITY] == YES)
        {
          VELOC_Augmented_Lk_Velocity_Pre(tree->n_root,tree->n_root->v[1],tree);
          VELOC_Augmented_Lk_Velocity_Pre(tree->n_root,tree->n_root->v[2],tree);
        }
          
      for(i=0;i<tree->mmod->n_dim;++i)
        {
          start = Contmod_Start(VELOCITY,i,tree);

          tree->contmod->lnL[VELOCITY*tree->mmod->n_dim+i] = tree->contmod->lnL_down[start + tree->n_root->num];
          tree->contmod->lnL[VELOCITY*tree->mmod->n_dim+i] += Log_Dnorm(tree->n_root->ldsk->veloc->deriv[i],
                                                                        root_mean,
                                                                        sqrt(root_var),
                                                                        &err);

          lnL += tree->contmod->lnL[VELOCITY*tree->mmod->n_dim+i];
        }
    }
  else
    {
      lnL = 0.0;
      
      for(i=0;i<tree->mmod->n_dim;++i)
        {
          start = Contmod_Start(VELOCITY,i,tree);

          tree->contmod->lnL[VELOCITY * tree->mmod->n_dim + i] = tree->contmod->lnL_down[start + z->num];
          tree->contmod->lnL[VELOCITY * tree->mmod->n_dim + i] += tree->contmod->lnL_up[start + z->num];
          if (z == tree->n_root)
            tree->contmod->lnL[VELOCITY * tree->mmod->n_dim + i] += Log_Dnorm(z->ldsk->veloc->deriv[i],
                                                                              root_mean,
                                                                              sqrt(root_var),
                                                                              &err);
          else
            tree->contmod->lnL[VELOCITY * tree->mmod->n_dim + i] += Log_Dnorm(z->ldsk->veloc->deriv[i],
                                                                              VELOC_Velocity_Mean_Along_Edge(z, i, tree),
                                                                              sqrt(VELOC_Velocity_Variance_Along_Edge(z, i, tree)), &err);

          // if(z != tree->n_root) tree->contmod->lnL[VELOCITY*tree->mmod->n_dim+i] += tree->contmod->lnL_up[start + z->num];
          // tree->contmod->lnL[VELOCITY*tree->mmod->n_dim+i] += Log_Dnorm(z->ldsk->veloc->deriv[i],
          //                                                               root_mean,
          //                                                               sqrt(root_var),
          //                                                               &err);

          lnL += tree->contmod->lnL[VELOCITY*tree->mmod->n_dim+i];
        }
    }
              
  return(lnL);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
void VELOC_Augmented_Lk_Velocity_Post(t_node *a, t_node *d, t_tree *tree)
{
  int i;
  
  if(d->tax == TRUE) return;
  else
    {
      
      for(i=0;i<3;++i)
        {
          if(d->v[i] != a && !(a == tree->n_root && d->b[i] == tree->e_root))
            {
              VELOC_Augmented_Lk_Velocity_Post(d,d->v[i],tree);
            }
        }
      
      VELOC_Update_Lk_Velocity_Down(a,d,tree);
    }
  return;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void VELOC_Augmented_Lk_Velocity_Pre(t_node *a, t_node *d, t_tree *tree)
{
  int i;

  if(d == tree->n_root)
    {
      VELOC_Augmented_Lk_Velocity_Pre(tree->n_root,tree->n_root->v[1],tree);
      VELOC_Augmented_Lk_Velocity_Pre(tree->n_root,tree->n_root->v[2],tree);
      return;
    }
  
  VELOC_Update_Lk_Velocity_Up(a,d,tree);
  
  if(d->tax == TRUE) return;
  else
    {
      for(i=0;i<3;++i)
        {
          if(d->v[i] != a && !(a == tree->n_root && d->b[i] == tree->e_root))
            {
              VELOC_Augmented_Lk_Velocity_Pre(d,d->v[i],tree);
            }
        }
    }
  return;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void VELOC_Update_Lk_Velocity_Up(t_node *a, t_node *d, t_tree *tree)
{
  int i,start,err;
  t_node *dad, *son, *bro;
  

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
  assert(son->anc == dad);
  assert(bro != son);
  
  for(i=0;i<tree->mmod->n_dim;++i)
    {
      start = Contmod_Start(VELOCITY,i,tree);
  
      tree->contmod->lnL_up[start+son->num] =
        tree->contmod->lnL_up[start+dad->num] +
        tree->contmod->lnL_down[start+bro->num] ;
  
      tree->contmod->lnL_up[start+son->num] +=
        Log_Dnorm(bro->ldsk->veloc->deriv[i],VELOC_Velocity_Mean_Along_Edge(bro,i,tree),sqrt(VELOC_Velocity_Variance_Along_Edge(bro,i,tree)),&err);

      if (dad != tree->n_root)
        tree->contmod->lnL_up[start + son->num] +=
            Log_Dnorm(dad->ldsk->veloc->deriv[i], VELOC_Velocity_Mean_Along_Edge(dad, i, tree), sqrt(VELOC_Velocity_Variance_Along_Edge(dad, i, tree)), &err);
      else
        tree->contmod->lnL_up[start + son->num] += Log_Dnorm(dad->ldsk->veloc->deriv[i],
                                                             tree->mmod->rw_root_mean[VELOCITY],
                                                             sqrt(tree->mmod->rw_root_var[VELOCITY]),
                                                             &err);

      // tree->contmod->lnL_up[start+son->num] +=
      //   -Log_Dnorm(son->ldsk->veloc->deriv[i],tree->mmod->rw_root_mean[VELOCITY],sqrt(tree->mmod->rw_root_var[VELOCITY]),&err) +
      //    Log_Dnorm(dad->ldsk->veloc->deriv[i],tree->mmod->rw_root_mean[VELOCITY],sqrt(tree->mmod->rw_root_var[VELOCITY]),&err);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void VELOC_Update_Lk_Velocity_Down(t_node *a, t_node *d, t_tree *tree)
{
  int i,start,err;
  t_node *son, *bro, *dad;
  

  
  if(d->tax == YES) return;
  
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
      start = Contmod_Start(VELOCITY,i,tree);
      
      tree->contmod->lnL_down[start+dad->num] =
        tree->contmod->lnL_down[start+son->num] +
        tree->contmod->lnL_down[start+bro->num] ;
      
      tree->contmod->lnL_down[start+dad->num] +=
        Log_Dnorm(son->ldsk->veloc->deriv[i],VELOC_Velocity_Mean_Along_Edge(son,i,tree),sqrt(VELOC_Velocity_Variance_Along_Edge(son,i,tree)),&err) +
        Log_Dnorm(bro->ldsk->veloc->deriv[i],VELOC_Velocity_Mean_Along_Edge(bro,i,tree),sqrt(VELOC_Velocity_Variance_Along_Edge(bro,i,tree)),&err);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl VELOC_Integrated_Lk_Location(t_node *z, t_tree *tree)
{
  phydbl lnL;

  lnL = UNLIKELY;
  
  if(z == NULL)
    {
      phydbl root_mean,root_var;
      int i,err,start;

      RRW_Update_Normalization_Factor(tree);

      Init_Contmod_Locations(tree);
      
      lnL = 0.0;
            
      VELOC_Integrated_Lk_Location_Post(NULL,tree->n_root,tree);

      if(tree->contmod->both_sides[LOCATION] == YES)
        {
          VELOC_Integrated_Lk_Location_Pre(tree->n_root,tree->n_root->v[1],tree);
          VELOC_Integrated_Lk_Location_Pre(tree->n_root,tree->n_root->v[2],tree);
        }
          

      for(i=0;i<tree->mmod->n_dim;++i)
        {
          root_var  = tree->mmod->rw_root_var[LOCATION];
          root_mean = LOCATION_Mean_Lonlat(i,tree);
          
          start = Contmod_Start(LOCATION,i,tree);


          tree->contmod->lnL[LOCATION*tree->mmod->n_dim+i] = 0.0;
          tree->contmod->lnL[LOCATION*tree->mmod->n_dim+i] += tree->contmod->logrem_down[start + tree->n_root->num];
          tree->contmod->lnL[LOCATION*tree->mmod->n_dim+i] += Log_Dnorm(tree->contmod->mu_down[start + tree->n_root->num],
                                                                        root_mean,
                                                                        sqrt(root_var+tree->contmod->var_down[start + tree->n_root->num]),
                                                                        &err);

          lnL += tree->contmod->lnL[LOCATION*tree->mmod->n_dim+i];
        }
    }
  else
    {
      lnL = VELOC_Integrated_Lk_Location_Node(z,tree);
    }
              
  return(lnL);
    
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl VELOC_Integrated_Lk_Location_Node(t_node *z, t_tree *tree)
{
  phydbl lnL;
  int i,start,err;
  phydbl z_mean_down, z_mean_up;
  phydbl z_var_down, z_var_up;
  phydbl z_logrem_down, z_logrem_up;
  
  lnL = 0.0;
  
  for(i=0;i<tree->mmod->n_dim;++i)
    {
      start = Contmod_Start(LOCATION,i,tree);

      z_mean_down   = tree->contmod->mu_down[start+z->num];
      z_var_down    = tree->contmod->var_down[start+z->num];
      z_logrem_down = tree->contmod->logrem_down[start+z->num];

      if(z == tree->n_root)
        {
          z_mean_up   = LOCATION_Mean_Lonlat(i,tree);
          z_var_up    = tree->mmod->rw_root_var[LOCATION];
          z_logrem_up = 0.0;
        }
      else
        {
          z_mean_up   = tree->contmod->mu_up[start+z->num];
          z_var_up    = tree->contmod->var_up[start+z->num];
          z_logrem_up = tree->contmod->logrem_up[start+z->num];
        }
      
      tree->contmod->lnL[LOCATION*tree->mmod->n_dim+i] = 0.0;
      tree->contmod->lnL[LOCATION*tree->mmod->n_dim+i] += z_logrem_down + z_logrem_up;
      tree->contmod->lnL[LOCATION*tree->mmod->n_dim+i] += Log_Dnorm(z_mean_down,z_mean_up,sqrt(z_var_up+z_var_down),&err);

      lnL += tree->contmod->lnL[LOCATION*tree->mmod->n_dim+i];
    }
    
  return(lnL);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// See Pybus et al. 10.1073/pnas.1206598109 + see my technical notes

void VELOC_Integrated_Lk_Location_Post(t_node *a, t_node *d, t_tree *tree)
{
  if(d->tax == TRUE) return;
  else
    {
      int i;
      
      for(i=0;i<3;++i)
        {
          if(d->v[i] != a && !(a == tree->n_root && d->b[i] == tree->e_root))
            {
              VELOC_Integrated_Lk_Location_Post(d,d->v[i],tree);
            }
        }
      
      VELOC_Update_Lk_Location_Down(a,d,tree);
      
    } 
  return;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void VELOC_Integrated_Lk_Location_Pre(t_node *a, t_node *d, t_tree *tree)
{
  int i;

  if(d == tree->n_root)
    {
      VELOC_Integrated_Lk_Location_Pre(tree->n_root,tree->n_root->v[1],tree);
      VELOC_Integrated_Lk_Location_Pre(tree->n_root,tree->n_root->v[2],tree);
      return;
    }

  VELOC_Update_Lk_Location_Up(a,d,tree);
  
  if(d->tax == TRUE) return;
  else
    {
      for(i=0;i<3;++i)
        {
          if(d->v[i] != a && !(a == tree->n_root && d->b[i] == tree->e_root))
            {
              VELOC_Integrated_Lk_Location_Pre(d,d->v[i],tree);
            }
        }
    }
  return;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void VELOC_Update_Lk_Location_Up(t_node *a, t_node *d, t_tree *tree)
{
  int i,start;
  t_node *dad, *son, *bro;
  phydbl dad_mu_up, dad_var_up, dad_logrem_up;
  phydbl son_a, son_b, son_var;
  phydbl bro_a, bro_b, bro_mu_down, bro_var_down, bro_var, bro_logrem_down;
  phydbl dt_son, dt_bro;

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

  dt_son = tree->times->nd_t[son->num] - tree->times->nd_t[dad->num];
  dt_bro = tree->times->nd_t[bro->num] - tree->times->nd_t[dad->num];

  for(i=0;i<tree->mmod->n_dim;++i)
    {
      start = Contmod_Start(LOCATION,i,tree);

      bro_mu_down     = tree->contmod->mu_down[start+bro->num];
      bro_var_down    = tree->contmod->var_down[start+bro->num];
      bro_logrem_down = tree->contmod->logrem_down[start+bro->num];

      bro_var         = VELOC_Location_Variance_Along_Edge(bro,i,tree);
      son_var         = VELOC_Location_Variance_Along_Edge(son,i,tree);

      son_a           = 1.0;
      bro_a           = 1.0;

      // son_b           = VELOC_Location_Mean_Along_Edge(son,i,tree) - son_a*dad->ldsk->coord->lonlat[i];
      // bro_b           = VELOC_Location_Mean_Along_Edge(bro,i,tree) - bro_a*dad->ldsk->coord->lonlat[i];

      if (IBM_Is_Ibm(tree->mmod) == YES)
      {
          son_b = dt_son * dad->ldsk->veloc->deriv[i];
          bro_b = dt_bro * dad->ldsk->veloc->deriv[i];
      }
      else if (IOU_Is_Iou(tree->mmod) == YES)
      {
          son_b = (dad->ldsk->veloc->deriv[i] - tree->mmod->ou_mu[i]) * (1. - exp(-tree->mmod->ou_theta * dt_son)) / tree->mmod->ou_theta + tree->mmod->ou_mu[i] * dt_son;
          bro_b = (dad->ldsk->veloc->deriv[i] - tree->mmod->ou_mu[i]) * (1. - exp(-tree->mmod->ou_theta * dt_bro)) / tree->mmod->ou_theta + tree->mmod->ou_mu[i] * dt_bro;
      }

      if(dad != tree->n_root)
        {
          dad_mu_up     = tree->contmod->mu_up[start+dad->num];
          dad_var_up    = tree->contmod->var_up[start+dad->num];
          dad_logrem_up = tree->contmod->logrem_up[start+dad->num];
        }
      else
        {
          dad_mu_up     = LOCATION_Mean_Lonlat(i,tree);
          dad_var_up    = tree->mmod->rw_root_var[LOCATION];
          dad_logrem_up = 0.0;
        }
      
      
      if(IBM_Is_Ibm(tree->mmod) == YES)
        {
          IBM_Integrated_Location_Up(dad_mu_up,dad_var_up,dad_logrem_up,
                                     son_a,son_b,son_var,
                                     bro_a,bro_b,bro_mu_down,bro_var_down,bro_var,bro_logrem_down,
                                     tree->contmod->mu_up+start+son->num,
                                     tree->contmod->var_up+start+son->num,
                                     tree->contmod->logrem_up+start+son->num);
        }
      else if(IWN_Is_Iwn(tree->mmod) == YES)
        {
          PhyML_Printf("\n. Not implemented yet");
          assert(false);
        }
      else if(IOU_Is_Iou(tree->mmod) == YES)
        {
          IOU_Integrated_Location_Up(dad_mu_up,dad_var_up,dad_logrem_up,
                                     son_a,son_b,son_var,
                                     bro_a,bro_b,bro_mu_down,bro_var_down,bro_var,bro_logrem_down,
                                     tree->contmod->mu_up+start+son->num,
                                     tree->contmod->var_up+start+son->num,
                                     tree->contmod->logrem_up+start+son->num);
        }

      assert(isnan(tree->contmod->logrem_up[start+d->num]) == NO);
      assert(isnan(tree->contmod->mu_up[start+d->num]) == NO);
      assert(isnan(tree->contmod->var_up[start+d->num]) == NO && !(tree->contmod->var_down[start+d->num] < 0.0));

      
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void VELOC_Update_Lk_Location_Down(t_node *a, t_node *d, t_tree *tree)
{
  int i,start;
  t_node *son, *bro, *dad;
  phydbl son_mu_down,bro_mu_down;
  phydbl son_var_down,bro_var_down;
  phydbl son_var,bro_var;
  phydbl son_logrem_down,bro_logrem_down;
  phydbl son_a, bro_a;
  phydbl son_b, bro_b;
  phydbl dt_son, dt_bro;

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
  
  dt_son = tree->times->nd_t[son->num] - tree->times->nd_t[dad->num];
  dt_bro = tree->times->nd_t[bro->num] - tree->times->nd_t[dad->num];

  for(i=0;i<tree->mmod->n_dim;++i)
    {
      start = Contmod_Start(LOCATION,i,tree);
      
      son_mu_down = tree->contmod->mu_down[start + son->num];
      bro_mu_down = tree->contmod->mu_down[start + bro->num];
      
      son_var_down = tree->contmod->var_down[start + son->num];
      bro_var_down = tree->contmod->var_down[start + bro->num];
      
      son_logrem_down = tree->contmod->logrem_down[start + son->num];
      bro_logrem_down = tree->contmod->logrem_down[start + bro->num];
      
      son_var = VELOC_Location_Variance_Along_Edge(son,i,tree);
      bro_var = VELOC_Location_Variance_Along_Edge(bro,i,tree);
      
      son_a = 1.0;
      bro_a = 1.0;
      
      // son_b = VELOC_Location_Mean_Along_Edge(son,i,tree) - son_a*dad->ldsk->coord->lonlat[i];
      // bro_b = VELOC_Location_Mean_Along_Edge(bro,i,tree) - bro_a*dad->ldsk->coord->lonlat[i];
      if (IBM_Is_Ibm(tree->mmod) == YES)
      {
          son_b = dt_son * dad->ldsk->veloc->deriv[i]; 
          bro_b = dt_bro * dad->ldsk->veloc->deriv[i]; 
      }
      else if (IOU_Is_Iou(tree->mmod) == YES)
      {
          son_b = (dad->ldsk->veloc->deriv[i] - tree->mmod->ou_mu[i]) * (1. - exp(-tree->mmod->ou_theta*dt_son))/tree->mmod->ou_theta + tree->mmod->ou_mu[i]*dt_son;
          bro_b = (dad->ldsk->veloc->deriv[i] - tree->mmod->ou_mu[i]) * (1. - exp(-tree->mmod->ou_theta*dt_bro))/tree->mmod->ou_theta + tree->mmod->ou_mu[i]*dt_bro;
      }

      if(IBM_Is_Ibm(tree->mmod) == YES)
        {
          IBM_Integrated_Location_Down(son_a,son_b,son_mu_down,son_var_down,son_var,
                                       bro_a,bro_b,bro_mu_down,bro_var_down,bro_var,
                                       son_logrem_down,bro_logrem_down,
                                       tree->contmod->mu_down+start+dad->num,
                                       tree->contmod->var_down+start+dad->num,
                                       tree->contmod->logrem_down+start+dad->num);
        }
      else if(IWN_Is_Iwn(tree->mmod) == YES)
        {
          IWN_Integrated_Location_Down(dt_son,dt_bro,
                                       son_a,son_b,son_mu_down,son_var_down,son_var,
                                       bro_a,bro_b,bro_mu_down,bro_var_down,bro_var,
                                       son_logrem_down,bro_logrem_down,
                                       d->ldsk->veloc->deriv[i],son->ldsk->veloc->deriv[i],bro->ldsk->veloc->deriv[i],
                                       tree->mmod->omega,
                                       tree->contmod->mu_down+start+dad->num,
                                       tree->contmod->var_down+start+dad->num,
                                       tree->contmod->logrem_down+start+dad->num);
        }
      else if(IOU_Is_Iou(tree->mmod) == YES)
        {
          IOU_Integrated_Location_Down(son_a,son_b,son_mu_down,son_var_down,son_var,
                                       bro_a,bro_b,bro_mu_down,bro_var_down,bro_var,
                                       son_logrem_down,bro_logrem_down,
                                       tree->contmod->mu_down+start+dad->num,
                                       tree->contmod->var_down+start+dad->num,
                                       tree->contmod->logrem_down+start+dad->num);
        }
      
      // PhyML_Printf("\n. node %d mu_down: %g var_down: %g logrem_down: %g", dad->num, tree->contmod->mu_down[start+dad->num], tree->contmod->var_down[start+dad->num], tree->contmod->logrem_down[start+dad->num]);

      assert(isnan(tree->contmod->logrem_down[start+dad->num]) == NO);
      assert(isnan(tree->contmod->mu_down[start+dad->num]) == NO);
      assert(isnan(tree->contmod->var_down[start+dad->num]) == NO && !(tree->contmod->var_down[start+dad->num] < 0.0));      
    }  
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

void VELOC_Sample_Node_Locations_Marginal(t_tree *tree)
{
  int i;

  if (tree->mmod->do_location_sampling == NO) return;

  RRW_Update_Normalization_Factor(tree);

  VELOC_Integrated_Lk_Location_Post(NULL, tree->n_root, tree);
  VELOC_Integrated_Lk_Location_Pre(tree->n_root, tree->n_root->v[1], tree);
  VELOC_Integrated_Lk_Location_Pre(tree->n_root, tree->n_root->v[2], tree);

  for (i = tree->n_otu; i < 2 * tree->n_otu - 1; ++i)
      VELOC_Sample_One_Node_Location(tree->a_nodes[i], tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void VELOC_Sample_One_Node_Location(t_node *z, t_tree *tree)
{
  int i,start;
  phydbl var_down,var_up;
  phydbl mu_down,mu_up;
  phydbl var,mean;
  
  for(i=0;i<tree->mmod->n_dim;++i)
    {
      start = Contmod_Start(LOCATION,i,tree);
      var = mean = -1.;
      
      mu_down  = tree->contmod->mu_down[start+z->num];
      var_down = tree->contmod->var_down[start+z->num];

      if(z != tree->n_root)
        {
          mu_up  = tree->contmod->mu_up[start+z->num];
          var_up = tree->contmod->var_up[start+z->num];
        }
      else
        {
          mu_up  = LOCATION_Mean_Lonlat(i,tree); 
          var_up = tree->mmod->rw_root_var[LOCATION];
        }

      if(var_up > SMALL && var_down > SMALL)
        {
          var = 1./var_up + 1./var_down;
          var = 1./var;

          mean = var * (mu_down / var_down + mu_up / var_up);
        }
      else if(var_up > SMALL)
        {
          var  = 0.0;
          mean = mu_down;
        }
      else if(var_down > SMALL)
        {
          var  = 0.0;
          mean = mu_up;          
        }
      else
        {
          var  = 0.0;
          mean = mu_down;          
        }
      
      assert(isnan(mean) == NO);
      assert(isnan(var) == NO);

      z->ldsk->coord->lonlat[i] = Rnorm(mean,sqrt(var));

    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void VELOC_Sample_Node_Locations_Joint(t_tree *tree)
{
  int i,start;
  phydbl root_var, root_mean, var, mean;
    
  if (tree->mmod->do_location_sampling == NO) return;

  root_var = tree->mmod->rw_root_var[LOCATION];
  root_mean = -BIG;
  
  RRW_Update_Normalization_Factor(tree);

  Init_Contmod_Locations(tree);

  VELOC_Integrated_Lk_Location_Post(NULL, tree->n_root, tree);

  for(i=0;i<tree->mmod->n_dim;++i)
    {
      start = Contmod_Start(LOCATION,i,tree);
      
      var = 1./tree->contmod->var_down[start+tree->n_root->num] + 1./root_var;
      var = 1./var;
      
      root_mean = LOCATION_Mean_Lonlat(i,tree);

      mean = (tree->contmod->mu_down[start+tree->n_root->num]/tree->contmod->var_down[start+tree->n_root->num] +
              root_mean / root_var)*var;

      tree->n_root->ldsk->coord->lonlat[i] = Rnorm(mean,sqrt(var));

      assert(isnan(tree->contmod->var_down[start + tree->n_root->num]) == NO);
    }

  VELOC_Sample_Node_Locations_Joint_Post(tree->n_root, tree->n_root->v[1], tree);
  VELOC_Sample_Node_Locations_Joint_Post(tree->n_root, tree->n_root->v[2], tree);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void VELOC_Sample_Node_Locations_Joint_Post(t_node *a, t_node *d, t_tree *tree)
{

  if(d->tax == TRUE) return;
  else
  {
    int start, i;
    phydbl son_a, son_b, son_var, son_var_down, son_mu, son_mu_down, var, mean;

    for (int i = 0; i < 3; ++i)
      if (d->v[i] != a && !(a == tree->n_root && d->b[i] == tree->e_root))
        VELOC_Sample_Node_Locations_Joint_Post(d, d->v[i], tree);

    for (i = 0; i < tree->mmod->n_dim; ++i)
    {
      start = Contmod_Start(LOCATION, i, tree);

      son_a = 1.0;
      son_b = 0.0;

      son_var = VELOC_Location_Variance_Along_Edge(d, i, tree);
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
phydbl VELOC_Mean_Velocity(short int dim, t_tree *tree)
{
  if(VELOC_Is_Integrated_Velocity(tree->mmod) == NO) return(-1);
  else
    {
      int i;
      phydbl mean;

      mean = 0.0;
      for(i=0;i<2*tree->n_otu-1;++i) mean += tree->a_nodes[i]->ldsk->veloc->deriv[dim];
      return(mean / (phydbl)(2*tree->n_otu-1));
    }
  return(-1.);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* Take the average of speed taken at each node of the treez */
phydbl VELOC_Mean_Speed(t_tree *tree)
{
  if(VELOC_Is_Integrated_Velocity(tree->mmod) == NO) return(-1);
  else
    {
      int i;
      phydbl mean_speed;

      mean_speed = 0.0;
      for(i=0;i<2*tree->n_otu-1;++i) mean_speed += VELOC_Veloc_To_Speed(tree->a_nodes[i]->ldsk->veloc,tree);
      return(mean_speed / (phydbl)(2*tree->n_otu-1));
    }
  return(-1.);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl VELOC_Veloc_To_Speed(t_geo_veloc *v, t_tree *tree)
{
  t_geo_coord *p1,*p2;
  int i;
  phydbl speed;
  
  p1 = GEO_Make_Geo_Coord(tree->mmod->n_dim);
  GEO_Init_Coord(p1,tree->mmod->n_dim);
  for(i=0;i<tree->mmod->n_dim;++i) p1->lonlat[i] = 0.0;
  
  p2 = GEO_Make_Geo_Coord(tree->mmod->n_dim);
  GEO_Init_Coord(p2,tree->mmod->n_dim);
  for(i=0;i<tree->mmod->n_dim;++i) p2->lonlat[i] = v->deriv[i];

  speed = -1.0;
  
  switch(tree->mmod->dist_type)
    {
    case HAVERSINE :
      {
        speed = Haversine_Distance(p1,p2);
        break;
      }
    case EUCLIDEAN :
      {
        speed = Euclidean_Distance(p1,p2);
        break;
      }
    case MANHATTAN :
      {
        speed = Manhattan_Distance(p1,p2);
        break;
      }
    default : assert(false);
    }

  Free_Geo_Coord(p1);
  Free_Geo_Coord(p2);

  return(speed);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl PHYREX_Degrees_To_Km(phydbl deg, t_tree *tree)
{
  t_geo_coord *p1,*p2;
  phydbl km;
  
  p1 = GEO_Make_Geo_Coord(1);
  GEO_Init_Coord(p1,1);
  p1->lonlat[0] = 0.0;

  p2 = GEO_Make_Geo_Coord(1);
  GEO_Init_Coord(p2,1);
  p2->lonlat[0] = fabs(deg);
  
  km = -1.0;
  
  switch(tree->mmod->dist_type)
    {
    case HAVERSINE :
      {
        km = Haversine_Distance(p1,p2);
        break;
      }
    case EUCLIDEAN :
      {
        km = Euclidean_Distance(p1,p2);
        break;
      }
    case MANHATTAN :
      {
        km = Manhattan_Distance(p1,p2);
        break;
      }
    default : assert(false);
    }

  Free_Geo_Coord(p1);
  Free_Geo_Coord(p2);

  return(km);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
/* WARNING : Code below only works for IBM model */
void VELOC_Simulate_Velocities(t_tree *tree)
{
  phydbl root_mean,root_var;
 
  root_mean = 0.0;
  root_var = 100;

  for(int i=0;i<tree->mmod->n_dim;++i) tree->n_root->ldsk->veloc->deriv[i] = Rnorm(root_mean,sqrt(root_var));

  VELOC_Simulate_Velocities_Pre(tree->n_root->ldsk,tree->n_root->ldsk->next[0],tree);
  VELOC_Simulate_Velocities_Pre(tree->n_root->ldsk,tree->n_root->ldsk->next[1],tree);
      
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void VELOC_Simulate_Velocities_Pre(t_ldsk *a, t_ldsk *d, t_tree *tree)
{
  phydbl mean,var;
  int i;
  
  for(i=0;i<tree->mmod->n_dim;++i)
    {
      mean = a->veloc->deriv[i];

      var =
        fabs(a->disk->time - d->disk->time) *
        tree->mmod->sigsq[i];
      
      d->veloc->deriv[i] = Rnorm(mean,sqrt(var));
    }

  if(d->nd->tax) return;
  else
    {
      for(i=0;i<d->n_next;++i) VELOC_Simulate_Velocities_Pre(d,d->next[i],tree);
    }
}
