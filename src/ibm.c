/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include "ibm.h"

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl IBM_Lk(t_tree *tree)
{
  phydbl lnL_loc, lnL_veloc;

  RRW_Update_Normalization_Factor(tree);

  lnL_loc   = UNLIKELY;
  lnL_veloc = UNLIKELY;
  
  lnL_loc = (tree->mmod->integrateAncestralLocations == NO) ? (IBM_Lk_Locations(tree)) : (IBM_Integrated_Lk_Location(tree));
  lnL_veloc = IBM_Lk_Velocities(tree);
  
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

short int IBM_Is_Ibm(t_phyrex_mod *mod)
{
  if(mod->model_id == IBM ||
     mod->model_id == RIBM) return(YES);
  return(NO);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl IBM_Sample_Velocities_And_Locations(int *node_order, int n_nodes, t_tree *tree)
{
  phydbl log_hr;
  int i,j,idx;
  
  RRW_Update_Normalization_Factor(tree);

  log_hr = 0.0;

  for(j=0;j<n_nodes;++j)
    {
      idx = (node_order != NULL) ? node_order[j] : j;

      for(i=0;i<tree->mmod->n_dim;++i)
        {
          log_hr += IBM_Velocity_One_Node(tree->a_nodes[idx],YES,i,tree);
          log_hr += IBM_Location_One_Node(tree->a_nodes[idx],YES,i,tree);
        }
    }
  
  return(log_hr);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl IBM_Velocities_Conditional(short int sample, int *node_order, short int dim, t_tree *tree)
{
  int j,idx;
  phydbl log_hr;
    
  log_hr = 0.0;
  idx = -1;
  
  for(j=0;j<2*tree->n_otu-1;++j)
    {
      idx = (node_order != NULL) ? node_order[j] : j;
      log_hr += IBM_Velocity_One_Node(tree->a_nodes[idx],sample,dim,tree);
    }

  return(log_hr);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl IBM_Velocity_One_Node(t_node *n, short int sample, short int dim, t_tree *tree)
{
  /*          w
              |
              |
              z
             / \
            /   \
           u     v
     
     The x's are the locations. The y's the velocities
  */ 

  t_node *z,*u,*v;
  int k;
  phydbl mean,var,dt,zi;
  phydbl xu,xv,xz,xw;
  phydbl yu,yv,yw;
  phydbl tu,tv,tz;
  phydbl muu,muv,muz;
  phydbl sigsqu,sigsqv,sigsqz;
  phydbl log_hr;
  short int diru;
  int err;

  err = NO;
  mean = var = -1.;
  
  xu = xv = xz = xw = -1.;
  yu = yv = yw = -1.;
  tu = tv = tz = -1.;
  muu = muv = muz = -1.;
  sigsqu = sigsqv = sigsqz = -1.;
  diru = -1;
  log_hr = 0.0;
  z = n;
  u = v = NULL;
  
  /* for(dim=0;dim<tree->mmod->n_dim;++dim) */
    {
      if(n->tax == NO && n != tree->n_root) /* Internal node that is not root node */
        {
          xz = n->ldsk->coord->lonlat[dim];
      
          diru = -1.;
          for(k=0;k<3;++k)
            {
              if(n->v[k] != n->anc && n->b[k] != tree->e_root)
                {
                  if(diru < 0)
                    {
                      diru = k;
                      u  = n->v[k];
                      xu = n->v[k]->ldsk->coord->lonlat[dim];
                      yu = n->v[k]->ldsk->veloc->deriv[dim];
                      tu = fabs(tree->times->nd_t[n->v[k]->num] - tree->times->nd_t[n->num]);
                    }
                  else
                    {
                      v  = n->v[k];
                      xv = n->v[k]->ldsk->coord->lonlat[dim];
                      yv = n->v[k]->ldsk->veloc->deriv[dim];
                      tv = fabs(tree->times->nd_t[n->v[k]->num] - tree->times->nd_t[n->num]);
                    }
                }
              else if(n->v[k] == n->anc)
                {
                  xw = n->v[k]->ldsk->coord->lonlat[dim];
                  yw = n->v[k]->ldsk->veloc->deriv[dim];
                  tz = fabs(tree->times->nd_t[n->v[k]->num] - tree->times->nd_t[n->num]);
                }
              else if(n->b[k] == tree->e_root)
                {
                  xw = tree->n_root->ldsk->coord->lonlat[dim];
                  yw = tree->n_root->ldsk->veloc->deriv[dim];
                  tz = fabs(tree->times->nd_t[tree->n_root->num] - tree->times->nd_t[n->num]);
                }
            }
          
          assert(tu > 0.);
          assert(tv > 0.);
          assert(tz > 0.);
          
          muu = 3.*(xu-xz)/(2.*tu) - yu/2.;
          
          sigsqu =
            log(tree->mmod->sigsq[dim]) +
            log(tree->mmod->sigsq_scale[u->num]) +
            log(tree->mmod->sigsq_scale_norm_fact) +
            log(tu)-
            LOG4;
          sigsqu = exp(sigsqu);
          
          
          muv = 3.*(xv-xz)/(2.*tv) - yv/2.;
          
          sigsqv =
            log(tree->mmod->sigsq[dim]) +
            log(tree->mmod->sigsq_scale[v->num]) +
            log(tree->mmod->sigsq_scale_norm_fact) +
            log(tv)-
            LOG4;
          sigsqv = exp(sigsqv);

          muz = 3.*(xz-xw)/(2.*tz) - yw/2.;
          
          sigsqz =
            log(tree->mmod->sigsq[dim]) +
            log(tree->mmod->sigsq_scale[z->num]) +
            log(tree->mmod->sigsq_scale_norm_fact) +
            log(tz)-
            LOG4;
          sigsqz = exp(sigsqz);
          
          
          assert(sigsqu > 0.0);
          assert(sigsqv > 0.0);
          assert(sigsqz > 0.0);
          
          var = 1./sigsqu + 1./sigsqv + 1./sigsqz;
          var = 1./var;
          
          mean = (muu/sigsqu + muv/sigsqv + muz/sigsqz)*var;
          
          /* PhyML_Printf("\n. VIT Node %d xu: %f yu: %f tu: %f xv: %f yv: %f tv: %f -- muu: %f sigsqu: %f | muv: %f sigsqv: %f | mean: %f var: %f [sigsq: %f]", */
          /*              n->num, */
          /*              xu,yu,tu, */
          /*              xv,yv,tv, */
          /*              muu,sigsqu, */
          /*              muv,sigsqv, */
          /*              mean,var, */
          /*              tree->mmod->sigsq[dim]); */
          
          assert(isinf(mean) == NO && isnan(mean) == NO);
          assert(isinf(var) == NO && isnan(var) == NO);
          
        }
      else if(n->tax == YES) /* Tip node */
        {
          dt = fabs(tree->times->nd_t[n->anc->num]-tree->times->nd_t[n->num]);
          zi = n->ldsk->coord->lonlat[dim] - n->anc->ldsk->coord->lonlat[dim];
          
          assert(dt > 0.0);
          
          mean = 0.5*(3.*zi/dt - n->anc->ldsk->veloc->deriv[dim]);
          
          var =
            log(dt) +
            log(tree->mmod->sigsq[dim]) +
            log(tree->mmod->sigsq_scale[n->num]) +
            log(tree->mmod->sigsq_scale_norm_fact) -
            LOG4;
          
          var = exp(var);
        }
      else if(n == tree->n_root) /* Root node */
        {
          xz = tree->n_root->ldsk->coord->lonlat[dim];
          
          u  = tree->n_root->v[1];
          xu = tree->n_root->v[1]->ldsk->coord->lonlat[dim];
          yu = tree->n_root->v[1]->ldsk->veloc->deriv[dim];
          tu = fabs(tree->times->nd_t[tree->n_root->v[1]->num] - tree->times->nd_t[tree->n_root->num]);
          
          v  = tree->n_root->v[2];
          xv = tree->n_root->v[2]->ldsk->coord->lonlat[dim];
          yv = tree->n_root->v[2]->ldsk->veloc->deriv[dim];
          tv = fabs(tree->times->nd_t[tree->n_root->v[2]->num] - tree->times->nd_t[tree->n_root->num]);
          

          
          assert(tu > 0.0);
          muu = 3.*(xu-xz)/(2.*tu) - yu/2.;
          sigsqu =
            log(tree->mmod->sigsq[dim]) +
            log(tree->mmod->sigsq_scale[u->num]) +
            log(tree->mmod->sigsq_scale_norm_fact) +
            log(tu)-
            LOG4;
          sigsqu = exp(sigsqu);
          
          assert(tv > 0.0);
          muv = 3.*(xv-xz)/(2.*tv) - yv/2.;
          sigsqv =
            log(tree->mmod->sigsq[dim]) +
            log(tree->mmod->sigsq_scale[v->num]) +
            log(tree->mmod->sigsq_scale_norm_fact) +
            log(tv)-
            LOG4;
          sigsqv = exp(sigsqv);
          
          
          assert(sigsqu > 0.0);
          assert(sigsqv > 0.0);
          
          var = 1./sigsqu + 1./sigsqv;
          var = 1./var;
          
          
          mean = (muu/sigsqu + muv/sigsqv)*var;
        }
      
      
      if(tree->mmod->print_lk == YES)
        PhyML_Printf("\n. VIT %d Root xz: %f xu: %f yu: %f tu: %f xv: %f yv: %f tv: %f -- muu: %f sigsqu: %f | muv: %f sigsqv: %f | mean: %f var: %f",
                     dim,
                     xz,
                     xu,yu,tu,
                     xv,yv,tv,
                     muu,sigsqu,
                     muv,sigsqv,
                     mean,var);
      
      
      assert(isinf(mean) == NO && isnan(mean) == NO);
      assert(isinf(var) == NO && isnan(var) == NO);


      /* !!!!!!!!!!!!!!!!!!!!1 */
      
      /* log_hr += Log_Dnorm_Trunc(n->ldsk->veloc->deriv[dim],mean,sqrt(var),tree->mmod->min_veloc,tree->mmod->max_veloc,&err); */
      
      /* if(sample == YES) */
      /*   { */
      /*     n->ldsk->veloc->deriv[dim] = Rnorm_Trunc(mean,sqrt(var),tree->mmod->min_veloc,tree->mmod->max_veloc,&err); */
      /*     if(err == YES) assert(false); */
      /*   } */
      
      /* log_hr -= Log_Dnorm_Trunc(n->ldsk->veloc->deriv[dim],mean,sqrt(var),tree->mmod->min_veloc,tree->mmod->max_veloc,&err); */


      log_hr += Log_Dnorm(n->ldsk->veloc->deriv[dim],mean,sqrt(var),&err);
      
      if(sample == YES)
        {
          n->ldsk->veloc->deriv[dim] = Rnorm(mean,sqrt(var));

          if(n->ldsk->veloc->deriv[dim] > tree->mmod->max_veloc)
            {
              n->ldsk->veloc->deriv[dim] = tree->mmod->max_veloc;
            }
          
          if(n->ldsk->veloc->deriv[dim] < tree->mmod->min_veloc)
            {
              n->ldsk->veloc->deriv[dim] = tree->mmod->min_veloc;
            }          
        }
      
      log_hr -= Log_Dnorm(n->ldsk->veloc->deriv[dim],mean,sqrt(var),&err);
    }

  return(log_hr);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl IBM_Locations_Conditional(short int sample, int *node_order, short int dim, t_tree *tree)
{
  int j,idx;
  phydbl log_hr;
  
  idx = -1;  
  log_hr = 0.0;

  
  for(j=0;j<2*tree->n_otu-1;++j)
    {
      idx = (node_order != NULL) ? node_order[j] : j;
      log_hr += IBM_Location_One_Node(tree->a_nodes[idx],sample,dim,tree);
    }

  return(log_hr);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl IBM_Location_One_Node(t_node *n, short int sample, short int dim, t_tree *tree)
{
  t_node *z,*u,*v;
  int k;
  phydbl mean,var;
  phydbl xu,xv,xw;
  phydbl yu,yv,yz,yw;
  phydbl tu,tv,tz;
  phydbl muu,muv,muz;
  phydbl sigsqu,sigsqv,sigsqz;
  phydbl log_hr;
  short int diru,dirv;
  int err;
  
  /*          w
              |
              |
              z
             / \
            /   \
           u     v
     
     The x's are the locations. The y's the velocities
  */ 
  
  err = NO;
  mean = var = -1.;
  
  xu = xv = xw = -1.;
  yu = yv = yz = yw = -1.;
  tu = tv = tz = -1.;
  muu = muv = muz = -1.;
  sigsqu = sigsqv = sigsqz = -1.;
  diru = -1;
  
  log_hr = 0.0;
  z = n;
  u = v = NULL;
  
  /* Sample ancestral locations given locations at tips and velocities everywhere */
  /* for(dim=0;dim<tree->mmod->n_dim;++dim) */
    {  
      if(n->tax == NO && n != tree->n_root)
        {
          yz = n->ldsk->veloc->deriv[dim];
          
          diru = dirv = -1.;
          for(k=0;k<3;++k)
            {
              if(n->v[k] != n->anc && n->b[k] != tree->e_root)
                {
                  if(diru < 0)
                    {
                      diru = k;
                      u  = n->v[k];
                      xu = n->v[k]->ldsk->coord->lonlat[dim];
                      yu = n->v[k]->ldsk->veloc->deriv[dim];
                      tu = fabs(tree->times->nd_t[n->v[k]->num] - tree->times->nd_t[n->num]);
                    }
                  else
                    {
                      dirv = k;
                      v  = n->v[k];
                      xv = n->v[k]->ldsk->coord->lonlat[dim];
                      yv = n->v[k]->ldsk->veloc->deriv[dim];
                      tv = fabs(tree->times->nd_t[n->v[k]->num] - tree->times->nd_t[n->num]);
                    }
                }
              else if(n->v[k] == n->anc)
                {
                  xw = n->v[k]->ldsk->coord->lonlat[dim];
                  yw = n->v[k]->ldsk->veloc->deriv[dim];
                  tz = fabs(tree->times->nd_t[n->v[k]->num] - tree->times->nd_t[n->num]);
                }
              else if(n->b[k] == tree->e_root)
                {
                  xw = tree->n_root->ldsk->coord->lonlat[dim];
                  yw = tree->n_root->ldsk->veloc->deriv[dim];
                  tz = fabs(tree->times->nd_t[tree->n_root->num] - tree->times->nd_t[n->num]);
                }
            }

          assert(tu > 0.);
          assert(tv > 0.);
          assert(tz > 0.);
          
          muu = xu - (yu+yz)/2.*tu;
          
          sigsqu =
            log(tree->mmod->sigsq[dim]) +
            log(tree->mmod->sigsq_scale[u->num]) +
            log(tree->mmod->sigsq_scale_norm_fact) +
            3.*log(tu)-
            LOG12;
          sigsqu = exp(sigsqu);
          
          
          muv = xv - (yv+yz)/2.*tv;
          
          sigsqv =
            log(tree->mmod->sigsq[dim]) +
            log(tree->mmod->sigsq_scale[v->num]) +
            log(tree->mmod->sigsq_scale_norm_fact) +
            3.*log(tv)-
            LOG12;
          sigsqv = exp(sigsqv);
          
          
          muz = xw + (yw+yz)/2.*tz;
          
          sigsqz =
            log(tree->mmod->sigsq[dim]) +
            log(tree->mmod->sigsq_scale[z->num]) +
            log(tree->mmod->sigsq_scale_norm_fact) +
            3.*log(tz)-
            LOG12;
          sigsqz = exp(sigsqz);
          
          
          assert(sigsqu > 0.0);
          assert(sigsqv > 0.0);
          assert(sigsqz > 0.0);
          
          var = 1./sigsqu + 1./sigsqv + 1./sigsqz;
          var = 1./var;
          
          mean = (muu/sigsqu + muv/sigsqv + muz/sigsqz)*var;
          
          if(tree->mmod->print_lk == YES)
            {
              PhyML_Printf("\n. LOC Node %d xu: %f yu: %f tu: %f xv: %f yv: %f tv: %f xw: %f yw: %f tz: %f -- muu: %f sigsqu: %f | muv: %f sigsqv: %f | muz: %f sigsqz: %f mean: %f var: %f [sigsq: %f] --> %f",
                           n->num,
                           xu,yu,tu,
                           xv,yv,tv,
                           xw,yw,tz,
                           muu,sigsqu,
                           muv,sigsqv,
                           muz,sigsqz,
                           mean,var,
                           tree->mmod->sigsq[dim],tree->mmod->c_lnL);
            }
          
          /*     PhyML_Printf(" %f --> %f",dum,n->ldsk->coord->lonlat[dim]); */
          
        }
      else if(n == tree->n_root) /* Root node */
        {
          yz = tree->n_root->ldsk->veloc->deriv[dim];

          u  = tree->n_root->v[1];
          xu = tree->n_root->v[1]->ldsk->coord->lonlat[dim];
          yu = tree->n_root->v[1]->ldsk->veloc->deriv[dim];
          tu = fabs(tree->times->nd_t[tree->n_root->v[1]->num] - tree->times->nd_t[tree->n_root->num]); 
          
          v  = tree->n_root->v[2];
          xv = tree->n_root->v[2]->ldsk->coord->lonlat[dim];
          yv = tree->n_root->v[2]->ldsk->veloc->deriv[dim];
          tv = fabs(tree->times->nd_t[tree->n_root->v[2]->num] - tree->times->nd_t[tree->n_root->num]); 
          
          
          assert(tu > 0.0);
          
          muu = xu - (yu+yz)/2.*tu;
          
          sigsqu =
            log(tree->mmod->sigsq[dim]) +
            log(tree->mmod->sigsq_scale[u->num]) +
            log(tree->mmod->sigsq_scale_norm_fact) +
            3.*log(tu)-
            LOG12;
          sigsqu = exp(sigsqu);
          
          
          muv = xv - (yv+yz)/2.*tv;
          
          sigsqv = 
            log(tree->mmod->sigsq[dim]) +
            log(tree->mmod->sigsq_scale[v->num]) +
            log(tree->mmod->sigsq_scale_norm_fact) +
            3.*log(tv)-
            LOG12;
          sigsqv = exp(sigsqv);
          
          assert(sigsqu > 0.0);
          assert(sigsqv > 0.0);
          
          var = 1./sigsqu + 1./sigsqv;
          var = 1./var;
          
          mean = (muu/sigsqu + muv/sigsqv)*var;
        }
      
      assert(isinf(mean) == NO && isnan(mean) == NO);
      assert(isinf(var) == NO && isnan(var) == NO);
      
      if(n->tax == NO)
        {
          log_hr += Log_Dnorm(n->ldsk->coord->lonlat[dim],mean,sqrt(var),&err);
          
          if(sample == YES) n->ldsk->coord->lonlat[dim] = Rnorm(mean,sqrt(var));
          
          log_hr -= Log_Dnorm(n->ldsk->coord->lonlat[dim],mean,sqrt(var),&err);

          RRW_Update_Normalization_Factor(tree);
        }
      
      /* PhyML_Printf("\n. LOC yz: %f Root xu: %f yu: %f tu: %f xv: %f yv: %f tv: %f -- muu: %f sigsqu: %f | muv: %f sigsqv: %f | mean: %f var: %f --> %f", */
      /*              yz, */
      /*              xu,yu,tu, */
      /*              xv,yv,tv, */
      /*              muu,sigsqu, */
      /*              muv,sigsqv, */
      /*              mean,var, */
      /*              LOCATION_Lk(tree)); */
      /* PhyML_Printf(" --> %f",tree->n_root->ldsk->coord->lonlat[dim]); */
    }

  return(log_hr);

}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
  
/* Returns log[Pr(locations at tips | velocities at all nodes, tree)] */
phydbl IBM_Lk_Locations(t_tree *tree)
{  
  phydbl lnP,disk_lnP;
  t_dsk *disk;
  
  disk_lnP = 0.0;
  lnP = 0.0;
  disk = tree->young_disk;
  
  do
    {
      disk_lnP = IBM_Lk_Locations_Core(disk,tree);
      lnP += disk_lnP;
      disk = disk->prev;
    }
  while(disk);
  
  return(lnP);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl IBM_Lk_Velocities(t_tree *tree)
{
  phydbl lnP,disk_lnP;
  t_dsk *disk;

  disk_lnP = 0.0;
  lnP = 0.0;
  disk = tree->young_disk;
  
  do
    {
      disk_lnP = IBM_Lk_Velocities_Core(disk,tree);
      lnP += disk_lnP;      
      disk = disk->prev;
    }  
  while(disk);
  
  return(lnP);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void IBM_Lk_Locations_Post(t_node *a, t_node *d, phydbl sigsq, int dim, t_tree *tree, short int print)
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
      
      for(i=0;i<3;++i)
        {
          if(d->v[i] != a && !(a == tree->n_root && d->b[i] == tree->e_root))
            {
              IBM_Lk_Locations_Post(d,d->v[i],sigsq,dim,tree,print);
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
      
      dv1var =
        log(sigsq) +
        log(tree->mmod->sigsq_scale[v1->num]) + 
        -log(12.) + 
        3.*log(fabs(tree->times->nd_t[v1->num]-tree->times->nd_t[d->num]));

      dv1var = exp(dv1var);
      

      dv2var =
        log(sigsq) +
        log(tree->mmod->sigsq_scale[v2->num]) + 
        -LOG12 +
        3.*log(fabs(tree->times->nd_t[v2->num]-tree->times->nd_t[d->num]));

      dv2var = exp(dv2var);
      
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
      
      tree->contmod->mu_down[d->num] = (v1mu*(v2var+dv2var) + v2mu*(v1var+dv1var))/(v2var+dv2var+v1var+dv1var);
      
      tree->contmod->var_down[d->num] = (v2var+dv2var)*(v1var+dv1var)/(v2var+dv2var+v1var+dv1var);

      tree->contmod->logrem_down[d->num] = v1logrem + v2logrem;
      tree->contmod->logrem_down[d->num] -= .5*log(2.*PI*(v2var+dv2var+v1var+dv1var));
      tree->contmod->logrem_down[d->num] -= .5*pow(v1mu-v2mu,2)/(v2var+dv2var+v1var+dv1var);      
    }
  
  return;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl IBM_Lk_Locations_Core(t_dsk *disk, t_tree *tree)
{
  phydbl lnP;
  int i,err;
  phydbl root_mean,root_var;

  root_mean = 0.0;
  root_var = 1000.0;
  
  assert(disk);
  
  if(disk == tree->young_disk) return 0.0;
  if(disk->age_fixed == YES) return 0.0;

  lnP = 0.0;

  if(disk->ldsk != NULL)
    {
      for(i=0;i<disk->ldsk->n_next;++i)
        {          
          lnP += IBM_Locations_Forward_Lk_Path(disk->ldsk,disk->ldsk->next[i],tree);
        }

      if(disk->ldsk->prev == NULL) /* root node */
        {
          for(i=0;i<tree->mmod->n_dim;++i)
            {
              root_mean = LOCATION_Mean_Lonlat(i,tree);
              lnP += Log_Dnorm(disk->ldsk->veloc->deriv[i],root_mean,sqrt(root_var),&err); /* Marginal density of velocity at root */
            }
        }
    }

  return(lnP);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl IBM_Lk_Velocities_Core(t_dsk *disk, t_tree *tree)
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
          lnP += IBM_Velocities_Forward_Lk_Path(disk->ldsk,disk->ldsk->next[i],tree);
        }
    }

  return(lnP);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl IBM_Velocities_Forward_Lk_Path(t_ldsk *a, t_ldsk *d, t_tree *tree)
{
  t_ldsk *ldsk;
  phydbl lnP,ld,disk_lnP,sd,dt,mean;
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
  mean = sd = -1.;
  
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

          mean = ldsk->prev->veloc->deriv[i];

          sd =
            log(tree->mmod->sigsq[i]) +
            log(tree->mmod->sigsq_scale[nd_d->num]) +
            log(tree->mmod->sigsq_scale_norm_fact) +
            log(dt);

          sd = sqrt(exp(sd));
          
          ld = ldsk->veloc->deriv[i];
          
          disk_lnP += Log_Dnorm(ld,mean,sd,&err);
          
          if(tree->mmod->print_lk == YES)
            PhyML_Printf("\n. VEL %2d Time: %10f d->coord: %10f a->coord: %10f a->veloc: %10f mean: %10f sd: %10f dt: %10f [%10f; %10f] x: %10f lk: %10f",
                         i,
                         ldsk->disk->time,
                         ldsk->coord->lonlat[i],
                         ldsk->prev->coord->lonlat[i],
                         ldsk->prev->veloc->deriv[i],
                         mean,
                         sd,
                         dt,
                         ldsk->disk->time,ldsk->prev->disk->time,
                         ld,
                         disk_lnP);

          if(isinf(lnP) || isnan(lnP)) return(UNLIKELY);
          
          if(isnan(lnP))
            {
              PhyML_Printf("\n. ld=%f sd=%f dt=[%f,%f] sigsq=%f",ld,sd,ldsk->disk->time,ldsk->prev->disk->time,tree->mmod->sigsq);
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

phydbl IBM_Locations_Forward_Lk_Path(t_ldsk *a, t_ldsk *d, t_tree *tree)
{
  t_ldsk *ldsk;
  phydbl lnP,ld,disk_lnP,sd,dt,mean;
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
  mean = sd = -1.;
  
  do
    {
      assert(ldsk->prev);

      disk_lnP = 0.0;
      
      for(i=0;i<tree->mmod->n_dim;++i)
        {
          dt = fabs(ldsk->disk->time-ldsk->prev->disk->time);
          
          /* mean = ldsk->prev->coord->lonlat[i] + dt*ldsk->prev->veloc->deriv[i]; */

          /* sd = */
          /*   log(tree->mmod->sigsq[i]) + */
          /*   log(tree->mmod->sigsq_scale[nd_d->num]) + */
          /*   log(tree->mmod->sigsq_scale_norm_fact) + */
          /*   3.*log(dt)- */
          /*   LOG3; */

          mean = ldsk->prev->coord->lonlat[i] + .5*(ldsk->prev->veloc->deriv[i] + ldsk->veloc->deriv[i])*dt;

          sd =
            log(tree->mmod->sigsq[i]) +
            log(tree->mmod->sigsq_scale[nd_d->num]) +
            log(tree->mmod->sigsq_scale_norm_fact) +
            3.*log(dt)-
            LOG12;
          
          sd = sqrt(exp(sd));

          
          ld = ldsk->coord->lonlat[i];

          disk_lnP += Log_Dnorm(ld,mean,sd,&err);

          if(tree->mmod->print_lk == YES)
            PhyML_Printf("\n. LOC %2d Time: %10f nd_d: %d sigsq: (%f,%f,%f) a->coord: %10f a->veloc: %10f mean: %10f sd: %10f dt: %10f [%10f %10f] x: %10f lk: %10f",
                         i,
                         ldsk->disk->time,
                         nd_d->num,
                         tree->mmod->sigsq[i],
                         tree->mmod->sigsq_scale[nd_d->num],
                         tree->mmod->sigsq_scale_norm_fact,
                         ldsk->prev->coord->lonlat[i],
                         ldsk->prev->veloc->deriv[i],
                         mean,
                         sd,
                         dt,
                         ldsk->disk->time,ldsk->prev->disk->time,
                         ld,
                         disk_lnP);
          
          
          if(isinf(lnP) || isnan(lnP)) return(UNLIKELY);
          
          if(isnan(lnP))
            {
              PhyML_Printf("\n. ld=%f sd=%f dt=[%f,%f] sigsq=%f",ld,sd,ldsk->disk->time,ldsk->prev->disk->time,tree->mmod->sigsq);
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

void IBM_Generate_Velocities_Then_Locations(t_tree *tree)
{
  RRW_Update_Normalization_Factor(tree);
  IBM_Generate_Velocities(tree);
  IBM_Generate_Locations_Given_Velocities(tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* Generate velocities under a Brownian process. Locations are 
   not known */
void IBM_Generate_Velocities(t_tree *tree)
{
  IBM_Generate_Velocities_Pre(tree->n_root,tree->n_root->v[1],tree);
  IBM_Generate_Velocities_Pre(tree->n_root,tree->n_root->v[2],tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void IBM_Generate_Velocities_Pre(t_node *a, t_node *d, t_tree *tree)
{
  int i;
  phydbl mean,sd,dt;

  mean = sd = dt = -1.;
  
  if(a == tree->n_root)
    {
      for(i=0;i<tree->mmod->n_dim;++i)
        {
          a->ldsk->veloc->deriv[i] = Rnorm(0.0,1.);
        }
    }
  
  for(i=0;i<tree->mmod->n_dim;++i)
    {
      dt = fabs(tree->times->nd_t[d->num] - tree->times->nd_t[a->num]);
      
      mean = a->ldsk->veloc->deriv[i];

      sd   =
        log(tree->mmod->sigsq[i]) +
        log(tree->mmod->sigsq_scale[i]) +
        log(tree->mmod->sigsq_scale_norm_fact) +
        log(dt);
      
      sd = exp(sd);
      sd = sqrt(sd);
        
      d->ldsk->veloc->deriv[i] = Rnorm(mean,sd);
    }
  
  if(d->tax == YES) return;
  else
    {
      for(i = 0; i<3; ++i)
        {
          if(d->v[i] != a && !(a == tree->n_root && d->b[i] == tree->e_root))
            {
              IBM_Generate_Velocities_Pre(d,d->v[i],tree);
            }
        }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void IBM_Generate_Locations_Given_Velocities(t_tree *tree)
{
  IBM_Generate_Locations_Given_Velocities_Pre(tree->n_root,tree->n_root->v[1],tree);
  IBM_Generate_Locations_Given_Velocities_Pre(tree->n_root,tree->n_root->v[2],tree);
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void IBM_Generate_Locations_Given_Velocities_Pre(t_node *a, t_node *d, t_tree *tree)
{
  int i;
  phydbl mean,sd,dt;

  mean = sd = dt = -1.;
  
  if(a == tree->n_root)
    {
      for(i=0;i<tree->mmod->n_dim;++i)
        {
          a->ldsk->coord->lonlat[i] = Rnorm(0.0,1.0);
        }
    }
  
  for(i=0;i<tree->mmod->n_dim;++i)
    {
      dt = fabs(tree->times->nd_t[d->num] - tree->times->nd_t[a->num]);
      
      mean = a->ldsk->coord->lonlat[i] + dt * (a->ldsk->veloc->deriv[i]+d->ldsk->veloc->deriv[i])/2.;

      sd   =
        log(tree->mmod->sigsq[i]) +
        log(tree->mmod->sigsq_scale[i]) +
        log(tree->mmod->sigsq_scale_norm_fact) +
        3.*log(dt) -
        LOG12;
      
      sd = exp(sd);
      sd = sqrt(sd);
        
      d->ldsk->coord->lonlat[i] = Rnorm(mean,sd);
    }
  
  if(d->tax == YES) return;
  else
    {
      for(i = 0; i<3; ++i)
        {
          if(d->v[i] != a && !(a == tree->n_root && d->b[i] == tree->e_root))
            {
              IBM_Generate_Locations_Given_Velocities_Pre(d,d->v[i],tree);
            }
        }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void IBM_Sample_Node_Locations(t_tree *tree)
{
  int i,j;
  t_node *n;
  phydbl dt,au,bu,varu,var,mean;

  n = NULL;
  dt = au = bu = varu = var = mean = -1;
  
  RRW_Update_Normalization_Factor(tree);

  for(i=0;i<tree->mmod->n_dim;++i)
    {
      IBM_Init_Contmod_Locations(i,tree);
      IBM_Integrated_Lk_Location_Post(NULL,tree->n_root,i,tree,NO);
      IBM_Integrated_Lk_Location_Pre(tree->n_root,tree->n_root->v[1],i,tree);
      IBM_Integrated_Lk_Location_Pre(tree->n_root,tree->n_root->v[2],i,tree);

      for(j=0;j<2*tree->n_otu-2;++j)
        {
          if(tree->a_nodes[j]->tax == NO && tree->a_nodes[j] != tree->n_root)
            {
              n = tree->a_nodes[j];
              dt = fabs(tree->times->nd_t[n->num]-tree->times->nd_t[n->anc->num]);
              au = 1.0;
              bu = (n->ldsk->veloc->deriv[i] + n->anc->ldsk->veloc->deriv[i])*.5*dt;
              varu = IBM_Location_Variance_Along_Edge(n,i,tree); 
              var = 1./tree->contmod->var_down[n->num] + 1. / (pow(au,2)*tree->contmod->var_up[n->num]+varu);
              var = 1./var;
              mean = (tree->contmod->mu_down[n->num]/tree->contmod->var_down[n->num] + (au*tree->contmod->mu_up[n->num]+bu)/(pow(au,2)*tree->contmod->var_up[n->num]+varu))*var;
              /* PhyML_Printf("\n. Location node %d%c mean: %f var: %f",j,(tree->a_nodes[j] == tree->n_root->v[1] || tree->a_nodes[j] == tree->n_root->v[2])?'*':' ',mean,var); */
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

phydbl IBM_Integrated_Lk_Location(t_tree *tree)
{
  phydbl lnL,root_mean,root_var;
  int i,err;
  
  RRW_Update_Normalization_Factor(tree);
  
  lnL = 0.0;
  root_var = 1000.0;
  root_mean = 0.0;
  
  for(i=0;i<tree->mmod->n_dim;++i)
    {
      IBM_Init_Contmod_Locations(i,tree);
      IBM_Integrated_Lk_Location_Post(NULL,tree->n_root,i,tree,NO);
      root_mean = LOCATION_Mean_Lonlat(i,tree);
      lnL += tree->contmod->logrem_down[tree->n_root->num];
      lnL += Log_Dnorm(tree->contmod->mu_down[tree->n_root->num],root_mean,sqrt(root_var+tree->contmod->var_down[tree->n_root->num]),&err);
    }
  return(lnL);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl IBM_Location_Variance_Along_Edge(t_node *d, short int dim, t_tree *tree)
{
  phydbl logvar;

  assert(d != tree->n_root);
  
  logvar =
    log(tree->mmod->sigsq[dim]) +
    log(tree->mmod->sigsq_scale[d->num]) +
    log(tree->mmod->sigsq_scale_norm_fact) +
    3.*log(fabs(tree->times->nd_t[d->num]-tree->times->nd_t[d->anc->num])) -
    LOG12;


  return(exp(logvar));
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// See Pybus et al. 10.1073/pnas.1206598109 + see my technical notes
void IBM_Integrated_Lk_Location_Post(t_node *a, t_node *d, short int dim, t_tree *tree, short int print)
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
              IBM_Integrated_Lk_Location_Post(d,d->v[i],dim,tree,print);
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

      err = -1;
      
      v1mu = tree->contmod->mu_down[v1->num];
      v2mu = tree->contmod->mu_down[v2->num];

      v1var = tree->contmod->var_down[v1->num];
      v2var = tree->contmod->var_down[v2->num];

      v1logrem = tree->contmod->logrem_down[v1->num];
      v2logrem = tree->contmod->logrem_down[v2->num];

      dv1var = IBM_Location_Variance_Along_Edge(v1,dim,tree);
      dv2var = IBM_Location_Variance_Along_Edge(v2,dim,tree);

      dtv1 = fabs(tree->times->nd_t[v1->num] - tree->times->nd_t[d->num]);
      dtv2 = fabs(tree->times->nd_t[v2->num] - tree->times->nd_t[d->num]);
      
      av1 = 1.0;
      av2 = 1.0;

      bv1 = (v1->ldsk->veloc->deriv[dim] + d->ldsk->veloc->deriv[dim])*.5*dtv1;
      bv2 = (v2->ldsk->veloc->deriv[dim] + d->ldsk->veloc->deriv[dim])*.5*dtv2;

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

void IBM_Integrated_Lk_Location_Pre(t_node *a, t_node *d, short int dim, t_tree *tree)
{
  int i,err;
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

  err = -1;
  
  assert(v2->anc == a);
  assert(v2->anc->anc == v1);
    
  v2mu     = tree->contmod->mu_down[v2->num];
  v2var    = tree->contmod->var_down[v2->num];
  v2logrem = tree->contmod->logrem_down[v2->num];
  dtv2     = fabs(tree->times->nd_t[v2->num] - tree->times->nd_t[a->num]);      
  av2var   = IBM_Location_Variance_Along_Edge(v2,dim,tree);    

  av2      = 1.0;
  bv2      = (v2->ldsk->veloc->deriv[dim] + a->ldsk->veloc->deriv[dim])*.5*dtv2;
      
  if(v1 != NULL)
    {
      v1mu     = tree->contmod->mu_up[a->num];
      v1var    = tree->contmod->var_up[a->num];
      v1logrem = tree->contmod->logrem_up[a->num];
      dtv1     = fabs(tree->times->nd_t[a->num] - tree->times->nd_t[v1->num]);          
      av1var   = IBM_Location_Variance_Along_Edge(a,dim,tree);

      av1      = 1.0;
      bv1      = (a->ldsk->veloc->deriv[dim] + v1->ldsk->veloc->deriv[dim])*.5*dtv1;
  
      tree->contmod->var_up[d->num] = pow(av2,2)/(v2var + av2var) + 1./(pow(av1,2)*v1var+av1var);
      tree->contmod->var_up[d->num] = 1./tree->contmod->var_up[d->num];

      tree->contmod->mu_up[d->num] = (av2*(v2mu-bv2)/(v2var + av2var) + (av1*v1mu+bv1)/(pow(av1,2)*v1var+av1var)) * tree->contmod->var_up[d->num];
      
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
              IBM_Integrated_Lk_Location_Pre(d,d->v[i],dim,tree);
            }
        }
    }
  return;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void IBM_Sample_Node_Velocities(t_tree *tree)
{
  int i,j;
  t_node *n;
  phydbl dt,au,bu,varu,var,mean;

  n = NULL;
  dt = au = bu = varu = var = mean = -1;

  RRW_Update_Normalization_Factor(tree);

  for(i=0;i<tree->mmod->n_dim;++i)
    {
      IBM_Init_Contmod_Velocities(i,tree);
      IBM_Integrated_Lk_Velocity_Post(NULL,tree->n_root,i,tree,NO);
      IBM_Integrated_Lk_Velocity_Pre(tree->n_root,tree->n_root->v[1],i,tree);
      IBM_Integrated_Lk_Velocity_Pre(tree->n_root,tree->n_root->v[2],i,tree);

      for(j=0;j<2*tree->n_otu-2;++j)
        {
          if(tree->a_nodes[j]->tax == NO && tree->a_nodes[j] != tree->n_root)
            {
              n = tree->a_nodes[j];
              dt = fabs(tree->times->nd_t[n->num]-tree->times->nd_t[n->anc->num]);
              au = -.5;
              bu = 1.5*(n->ldsk->coord->lonlat[i] + n->anc->ldsk->coord->lonlat[i])/dt;
              varu = IBM_Velocity_Variance_Along_Edge(n,i,tree); 
              var = 1./tree->contmod->var_down[n->num] + 1./(pow(au,2)*tree->contmod->var_up[n->num]+varu);
              var = 1./var;
              mean = (tree->contmod->mu_down[n->num]/tree->contmod->var_down[n->num] + (au*tree->contmod->mu_up[n->num]+bu)/(pow(au,2)*tree->contmod->var_up[n->num]+varu))*var;
              tree->a_nodes[j]->ldsk->veloc->deriv[i] = Rnorm(mean,sqrt(var));
            }
        }

      assert(isnan(tree->contmod->var_down[tree->n_root->num]) == NO);
      
      tree->n_root->ldsk->veloc->deriv[i] =
        Rnorm(tree->contmod->mu_down[tree->n_root->num],
              sqrt(tree->contmod->var_down[tree->n_root->num]));
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl IBM_Integrated_Lk_Velocity(t_tree *tree)
{
  phydbl lnL;
  int i;
  
  RRW_Update_Normalization_Factor(tree);

  lnL = 0.0;

  for(i=0;i<tree->mmod->n_dim;++i)
    {
      IBM_Init_Contmod_Velocities(i,tree);
      IBM_Integrated_Lk_Velocity_Post(NULL,tree->n_root,i,tree,NO);
      lnL += tree->contmod->logrem_down[tree->n_root->num];
    }
  return(lnL);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl IBM_Velocity_Variance_Along_Edge(t_node *d, short int dim, t_tree *tree)
{
  phydbl logvar;

  assert(d != tree->n_root);

  logvar =
    log(tree->mmod->sigsq[dim]) +
    log(tree->mmod->sigsq_scale[d->num]) + 
    log(tree->mmod->sigsq_scale_norm_fact) +
    log(fabs(tree->times->nd_t[d->num]-tree->times->nd_t[d->anc->num])) -
    LOG4;
  
  return(exp(logvar));
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// See Pybus et al. 10.1073/pnas.1206598109 + see my technical notes
void IBM_Integrated_Lk_Velocity_Post(t_node *a, t_node *d, short int dim, t_tree *tree, short int print)
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
              IBM_Integrated_Lk_Velocity_Post(d,d->v[i],dim,tree,print);
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

      dv1var = IBM_Velocity_Variance_Along_Edge(v1,dim,tree);
      dv2var = IBM_Velocity_Variance_Along_Edge(v2,dim,tree);

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

void IBM_Integrated_Lk_Velocity_Pre(t_node *a, t_node *d, short int dim, t_tree *tree)
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
  av2var   = IBM_Velocity_Variance_Along_Edge(v2,dim,tree);
  av2      = -.5;
  bv2      = 1.5*(v2->ldsk->coord->lonlat[dim] - a->ldsk->coord->lonlat[dim])/dtv2;
  
  
  if(v1 != NULL)
    {
      v1mu     = tree->contmod->mu_up[a->num];
      v1var    = tree->contmod->var_up[a->num];
      v1logrem = tree->contmod->logrem_up[a->num];
      dtv1     = fabs(tree->times->nd_t[a->num] - tree->times->nd_t[v1->num]);
      av1var   = IBM_Velocity_Variance_Along_Edge(a,dim,tree);
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
              IBM_Integrated_Lk_Velocity_Pre(d,d->v[i],dim,tree);
            }
        }
    }
  return;
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

