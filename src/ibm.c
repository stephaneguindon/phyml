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

short int IBM_Is_Ibm(t_phyrex_mod *mod)
{
  if(mod->model_id == IBM ||
     mod->model_id == RIBM) return(YES);
  return(NO);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl IBM_Prior(t_tree *tree)
{
  tree->mmod->c_lnP = RW_Prior(tree);
  return(tree->mmod->c_lnP);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl IBM_Velocity_Mean_Along_Edge(t_node *d, short int dim, t_tree *tree)
{
  assert(d != tree->n_root);
  return(d->anc->ldsk->veloc->deriv[dim]);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl IBM_Velocity_Variance_Along_Edge(t_node *d, short int dim, t_tree *tree)
{
  phydbl dt,logvar,eps;

  eps = 0.0;
  
  assert(d != tree->n_root);
  
  dt = fabs(tree->times->nd_t[d->num] - tree->times->nd_t[d->anc->num]+eps);
  
  logvar =
    log(tree->mmod->sigsq[dim]) +
    log(tree->mmod->sigsq_scale[d->num]) +
    log(tree->mmod->sigsq_scale_norm_fact) +
    log(dt);

  return(exp(logvar));
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl IBM_Location_Mean_Along_Edge(t_node *d, short int dim, t_tree *tree)
{
  phydbl dt, loc_beg, veloc_beg, veloc_end, eps;
  
  assert(d != tree->n_root);

  eps = 0.0;
  
  dt = fabs(tree->times->nd_t[d->num] - tree->times->nd_t[d->anc->num] + eps);

  loc_beg = d->anc->ldsk->coord->lonlat[dim];
  veloc_beg = d->anc->ldsk->veloc->deriv[dim];
  veloc_end = d->ldsk->veloc->deriv[dim];
  
  return(loc_beg + .5*dt*(veloc_beg + veloc_end));
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
/* Variance of location at the end of the edge, given velocities
   at both extremities */

phydbl IBM_Location_Variance_Along_Edge(t_node *d, short int dim, t_tree *tree)
{
  phydbl logvar,eps;

  eps = 0.0;
  
  assert(d != tree->n_root);
  
  logvar =
    log(tree->mmod->sigsq[dim]) +
    log(tree->mmod->sigsq_scale[d->num]) +
    log(tree->mmod->sigsq_scale_norm_fact) +
    3.*log(fabs(tree->times->nd_t[d->num]-tree->times->nd_t[d->anc->num]+eps)) -
    LOG12;

  return(exp(logvar));
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void IBM_Integrated_Location_Down(phydbl son_a, phydbl son_b, phydbl son_mu_down, phydbl son_var_down, phydbl son_var,
                                  phydbl bro_a, phydbl bro_b, phydbl bro_mu_down, phydbl bro_var_down, phydbl bro_var,
                                  phydbl son_logrem, phydbl bro_logrem,
                                  phydbl *mean, phydbl *var, phydbl *logrem)
{
  RW_Integrated_Lk_Down(son_a, son_b, son_mu_down, son_var_down, son_var,
                        bro_a, bro_b, bro_mu_down, bro_var_down, bro_var,
                        son_logrem, bro_logrem,
                        mean, var, logrem);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void IBM_Integrated_Location_Up(phydbl dad_mu_up, phydbl dad_var_up, phydbl dad_logrem_up,
                                phydbl son_a, phydbl son_b, phydbl son_var,
                                phydbl bro_a, phydbl bro_b, phydbl bro_mu_down, phydbl bro_var_down, phydbl bro_var, phydbl bro_logrem_down,
                                phydbl *mean, phydbl *var, phydbl *logrem)
{  
  RW_Integrated_Lk_Up(dad_mu_up, dad_var_up, dad_logrem_up, 
                      son_a, son_b, son_var, 
                      bro_a, bro_b, bro_mu_down, bro_var_down, bro_var, bro_logrem_down, 
                      mean, var, logrem);
  
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////



/* DEPRECATED CODE BELOW */



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

void IBM_Augmented_Lk_Locations_Post(t_node *a, t_node *d, phydbl sigsq, int dim, t_tree *tree, short int print)
{
  if(d->tax == TRUE)
    {
      return;
    }
  else
    {
      int i,start;
      t_node *v1, *v2;
      phydbl v1mu,v2mu;
      phydbl v1var,v2var;
      phydbl dv1var,dv2var;
      phydbl v1logrem,v2logrem;
      
      for(i=0;i<3;++i)
        {
          if(d->v[i] != a && !(a == tree->n_root && d->b[i] == tree->e_root))
            {
              IBM_Augmented_Lk_Locations_Post(d,d->v[i],sigsq,dim,tree,print);
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

      start = Contmod_Start(LOCATION,dim,tree);

      v1mu = tree->contmod->mu_down[start+v1->num];      
      v2mu = tree->contmod->mu_down[start+v2->num];

      v1var = tree->contmod->var_down[start+v1->num];
      v2var = tree->contmod->var_down[start+v2->num];

      v1logrem = tree->contmod->logrem_down[start+v1->num];
      v2logrem = tree->contmod->logrem_down[start+v2->num];
      
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
      
      tree->contmod->mu_down[start+d->num] = (v1mu*(v2var+dv2var) + v2mu*(v1var+dv1var))/(v2var+dv2var+v1var+dv1var);
      
      tree->contmod->var_down[start+d->num] = (v2var+dv2var)*(v1var+dv1var)/(v2var+dv2var+v1var+dv1var);

      tree->contmod->logrem_down[start+d->num] = v1logrem + v2logrem;
      tree->contmod->logrem_down[start+d->num] -= .5*log(2.*PI*(v2var+dv2var+v1var+dv1var));
      tree->contmod->logrem_down[start+d->num] -= .5*pow(v1mu-v2mu,2)/(v2var+dv2var+v1var+dv1var);      
    }
  
  return;
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
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
