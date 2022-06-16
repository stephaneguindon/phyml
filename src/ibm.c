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

  lnL_loc   = UNLIKELY;
  lnL_veloc = UNLIKELY;

  lnL_loc = IBM_Lk_Locations(tree);
  lnL_veloc = IBM_Lk_Velocities(tree);

  /* PhyML_Printf("\n. Sigsq: %f %f",tree->mmod->sigsq[0],tree->mmod->sigsq[1]); */
  /* PhyML_Printf("\n. IBM locations: %f",lnL_loc); */
  /* PhyML_Printf("\n. IBM velocities: %f",lnL_veloc); */
  /* PhyML_Printf("\n. Sum: %f",lnL_loc + lnL_veloc); */
  return(lnL_loc + lnL_veloc);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void IBM_Sample_Velocities_Conditional(t_tree *tree)
{
  int i,j,k,print;
  phydbl mean,var,dt,zi;
  phydbl xu,xv,xz,xw;
  phydbl yu,yv,yw;
  phydbl tu,tv,tz;
  phydbl muu,muv,muz;
  phydbl sigsqu,sigsqv,sigsqz;
  short int diru;

  /*          w
              |
              |
              z
             / \
            /   \
           u     v
     
     The x's are the locations. The y's the velocities
  */ 
  
  print = NO;

  xu = xv = xz = xw = -1.;
  yu = yv = yw = -1.;
  tu = tv = tz = -1.;
  muu = muv = muz = -1.;
  sigsqu = sigsqv = sigsqz = -1.;
  diru = -1;
  
  tree->mmod->sigsq_scale_norm_fact = 1.0;

  
  /* Sample ancestral velocities given locations and velocities at tips */
  for(i=0;i<tree->mmod->n_dim;++i)
    {
      for(j=0;j<2*tree->n_otu-2;++j)
        {
          if(tree->a_nodes[j]->tax == NO && tree->a_nodes[j] != tree->n_root)
            {
              xz = tree->a_nodes[j]->ldsk->coord->lonlat[i];
              
              diru = -1.;
              for(k=0;k<3;++k)
                {
                  if(tree->a_nodes[j]->v[k] != tree->a_nodes[j]->anc && tree->a_nodes[j]->b[k] != tree->e_root)
                    {
                      if(diru < 0)
                        {
                          diru = k;
                          xu = tree->a_nodes[j]->v[k]->ldsk->coord->lonlat[i];
                          yu = tree->a_nodes[j]->v[k]->ldsk->veloc->deriv[i];
                          tu = fabs(tree->times->nd_t[tree->a_nodes[j]->v[k]->num] - tree->times->nd_t[tree->a_nodes[j]->num]);
                        }
                      else
                        {
                          xv = tree->a_nodes[j]->v[k]->ldsk->coord->lonlat[i];
                          yv = tree->a_nodes[j]->v[k]->ldsk->veloc->deriv[i];
                          tv = fabs(tree->times->nd_t[tree->a_nodes[j]->v[k]->num] - tree->times->nd_t[tree->a_nodes[j]->num]);
                        }
                    }
                  else
                    {
                      xw = tree->a_nodes[j]->v[k]->ldsk->coord->lonlat[i];
                      yw = tree->a_nodes[j]->v[k]->ldsk->veloc->deriv[i];
                      tz = fabs(tree->times->nd_t[tree->a_nodes[j]->v[k]->num] - tree->times->nd_t[tree->a_nodes[j]->num]);
                    }
                }

              assert(tu > 0.);
              assert(tv > 0.);
              assert(tz > 0.);
              
              muu = 3.*(xu-xz)/(2.*tu) - yu/2.;

              sigsqu =
                log(tree->mmod->sigsq[i]) +
                log(tree->mmod->sigsq_scale_norm_fact) +
                log(tu)-
                LOG4;
              sigsqu = exp(sigsqu);


              muv = 3.*(xv-xz)/(2.*tv) - yv/2.;

              sigsqv =
                log(tree->mmod->sigsq[i]) +
                log(tree->mmod->sigsq_scale_norm_fact) +
                log(tv)-
                LOG4;
              sigsqv = exp(sigsqv);



              muz = 3.*(xz-xw)/(2.*tz) - yw/2.;
              
              sigsqz =
                log(tree->mmod->sigsq[i]) +
                log(tree->mmod->sigsq_scale_norm_fact) +
                3.*log(tz) -
                LOG4;
              sigsqz = exp(sigsqz);

              
              assert(sigsqu > 0.0);
              assert(sigsqv > 0.0);
              assert(sigsqz > 0.0);
              
              var = 1./sigsqu + 1./sigsqv + 1./sigsqz;
              var = 1./var;
              
              mean = (muu/sigsqu + muv/sigsqv + muz/sigsqz)*var;
              
              /* PhyML_Printf("\n. VIT Node %d xu: %f yu: %f tu: %f xv: %f yv: %f tv: %f -- muu: %f sigsqu: %f | muv: %f sigsqv: %f | mean: %f var: %f [sigsq: %f] --> %f", */
              /*              j, */
              /*              xu,yu,tu, */
              /*              xv,yv,tv, */
              /*              muu,sigsqu, */
              /*              muv,sigsqv, */
              /*              mean,var, */
              /*              tree->mmod->sigsq[i], */
              /*              LOCATION_Lk(tree)); */

              assert(isinf(mean) == NO && isnan(mean) == NO);
              assert(isinf(var) == NO && isnan(var) == NO);
              
              tree->a_nodes[j]->ldsk->veloc->deriv[i] = Rnorm(mean,sqrt(var));

              if(tree->a_nodes[j]->ldsk->veloc->deriv[i] < tree->mmod->min_veloc)
                tree->a_nodes[j]->ldsk->veloc->deriv[i] = tree->mmod->min_veloc;

              if(tree->a_nodes[j]->ldsk->veloc->deriv[i] > tree->mmod->max_veloc)
                tree->a_nodes[j]->ldsk->veloc->deriv[i] = tree->mmod->max_veloc;
            }
        }
      
      assert(isnan(tree->contmod->var_down[tree->n_root->num]) == NO);


      /* Root node */

      xz = tree->n_root->ldsk->coord->lonlat[i];

      xu = tree->n_root->v[1]->ldsk->coord->lonlat[i];
      yu = tree->n_root->v[1]->ldsk->veloc->deriv[i];
      tu = fabs(tree->times->nd_t[tree->n_root->v[1]->num] - tree->times->nd_t[tree->n_root->num]);

      xv = tree->n_root->v[2]->ldsk->coord->lonlat[i];
      yv = tree->n_root->v[2]->ldsk->veloc->deriv[i];
      tv = fabs(tree->times->nd_t[tree->n_root->v[2]->num] - tree->times->nd_t[tree->n_root->num]);

      
      assert(tu > 0.0);
      muu = 3.*(xu-xz)/(2.*tu) - yu/2.;
      sigsqu =
        log(tree->mmod->sigsq[i]) +
        log(tree->mmod->sigsq_scale_norm_fact) +
        log(tu)-
        LOG4;
      sigsqu = exp(sigsqu);
      
      assert(tv > 0.0);
      muv = 3.*(xv-xz)/(2.*tv) - yv/2.;
      sigsqv =
        log(tree->mmod->sigsq[i]) +
        log(tree->mmod->sigsq_scale_norm_fact) +
        log(tv)-
        LOG4;
      sigsqv = exp(sigsqv);


      assert(sigsqu > 0.0);
      assert(sigsqv > 0.0);
      
      var = 1./sigsqu + 1./sigsqv;
      var = 1./var;

      
      mean = (muu/sigsqu + muv/sigsqv)*var;

      /* PhyML_Printf("\n. VIT %d Root xz: %f xu: %f yu: %f tu: %f xv: %f yv: %f tv: %f -- muu: %f sigsqu: %f | muv: %f sigsqv: %f | mean: %f var: %f --> %f", */
      /*              i, */
      /*              xz, */
      /*              xu,yu,tu, */
      /*              xv,yv,tv, */
      /*              muu,sigsqu, */
      /*              muv,sigsqv, */
      /*              mean,var,LOCATION_Lk(tree)); */

      
      assert(isinf(mean) == NO && isnan(mean) == NO);
      assert(isinf(var) == NO && isnan(var) == NO);
      
      tree->n_root->ldsk->veloc->deriv[i] = Rnorm(mean,sqrt(var));


      if(tree->n_root->ldsk->veloc->deriv[i] < tree->mmod->min_veloc)
        tree->n_root->ldsk->veloc->deriv[i] = tree->mmod->min_veloc;
      
      if(tree->n_root->ldsk->veloc->deriv[i] > tree->mmod->max_veloc)
        tree->n_root->ldsk->veloc->deriv[i] = tree->mmod->max_veloc;

      /* PhyML_Printf(" --> %f",tree->n_root->ldsk->veloc->deriv[i]); */
    }


  /* Sample velocities at tips given ancestral velocities */
  for(i=0;i<tree->mmod->n_dim;++i)
    {
      for(j=0;j<tree->n_otu;++j)
        {
          dt = fabs(tree->times->nd_t[tree->a_nodes[j]->v[0]->num]-tree->times->nd_t[tree->a_nodes[j]->num]);
          zi = tree->a_nodes[j]->ldsk->coord->lonlat[i] - tree->a_nodes[j]->v[0]->ldsk->coord->lonlat[i];
          
          assert(dt > 0.0);
          
          mean = 0.5*(3.*zi/dt - tree->a_nodes[j]->v[0]->ldsk->veloc->deriv[i]);

          var =
            log(dt) +
            log(tree->mmod->sigsq[i]) +
            log(tree->mmod->sigsq_scale_norm_fact) -
            LOG4;
          
          var = exp(var);

          assert(isinf(mean) == NO && isnan(mean) == NO);
          assert(isinf(var) == NO && isnan(var) == NO);
          
          tree->a_nodes[j]->ldsk->veloc->deriv[i] = Rnorm(mean,sqrt(var));

          if(tree->a_nodes[j]->ldsk->veloc->deriv[i] < tree->mmod->min_veloc)
            tree->a_nodes[j]->ldsk->veloc->deriv[i] = tree->mmod->min_veloc;
          
          if(tree->a_nodes[j]->ldsk->veloc->deriv[i] > tree->mmod->max_veloc)
            tree->a_nodes[j]->ldsk->veloc->deriv[i] = tree->mmod->max_veloc;
        }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void IBM_Sample_Locations(t_tree *tree)
{
  int i,j;
  
  /* Sample ancestral locations given locations at tips and velocities at all nodes */
  for(i=0;i<tree->mmod->n_dim;++i)
    {
      RRW_Init_Contmod_Locations(i,tree);
      Lk_Contmod_Post(NULL,tree->n_root,tree->mmod->sigsq[i],tree,NO);
      Lk_Contmod_Pre(tree->n_root,tree->n_root->v[1],tree->mmod->sigsq[i],tree);
      Lk_Contmod_Pre(tree->n_root,tree->n_root->v[2],tree->mmod->sigsq[i],tree);      
      
      for(j=0;j<2*tree->n_otu-2;++j)
        {
          if(tree->a_nodes[j]->tax == NO && tree->a_nodes[j] != tree->n_root)
            {
              tree->a_nodes[j]->ldsk->coord->lonlat[i] =
                Sample_Ancestral_Trait_Contmod(tree->a_nodes[j]->anc,
                                               tree->a_nodes[j],
                                               fabs(tree->times->nd_t[tree->a_nodes[j]->anc->num]-tree->times->nd_t[tree->a_nodes[j]->num]),
                                               0.0,
                                               -log(12.)+
                                               log(tree->mmod->sigsq[i]) +
                                               log(tree->mmod->sigsq_scale[tree->a_nodes[j]->num]) +
                                               log(pow(fabs(tree->times->nd_t[tree->a_nodes[j]->anc->num]-tree->times->nd_t[tree->a_nodes[j]->num]),2)),
                                               NO,
                                               tree);
            }
        }
      
      assert(isnan(tree->contmod->var_down[tree->n_root->num]) == NO);
      
      tree->n_root->ldsk->coord->lonlat[i] =
        Rnorm(tree->contmod->mu_down[tree->n_root->num],
              sqrt(tree->contmod->var_down[tree->n_root->num]));
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void IBM_Sample_Locations_Conditional(t_tree *tree)
{
  int i,j,k,print;
  phydbl mean,var;
  phydbl xu,xv,xw;
  phydbl yu,yv,yz,yw;
  phydbl tu,tv,tz;
  phydbl muu,muv,muz;
  phydbl sigsqu,sigsqv,sigsqz;
  short int diru,dirv;
  
  
  /*          w
              |
              |
              z
             / \
            /   \
           u     v
     
     The x's are the locations. The y's the velocities
  */ 
  
  print = NO;

  xu = xv = xw = -1.;
  yu = yv = yz = yw = -1.;
  tu = tv = tz = -1.;
  muu = muv = muz = -1.;
  sigsqu = sigsqv = sigsqz = -1.;
  diru = -1;
  
  tree->mmod->sigsq_scale_norm_fact = 1.0;
  
  /* Sample ancestral locations given locations at tips and velocities everywhere */
  for(i=0;i<tree->mmod->n_dim;++i)
    {
      for(j=0;j<2*tree->n_otu-2;++j)
        {
          if(tree->a_nodes[j]->tax == NO && tree->a_nodes[j] != tree->n_root)
            {
              yz = tree->a_nodes[j]->ldsk->veloc->deriv[i];
              
              diru = dirv = -1.;
              for(k=0;k<3;++k)
                {
                  if(tree->a_nodes[j]->v[k] != tree->a_nodes[j]->anc && tree->a_nodes[j]->b[k] != tree->e_root)
                    {
                      if(diru < 0)
                        {
                          diru = k;
                          xu = tree->a_nodes[j]->v[k]->ldsk->coord->lonlat[i];
                          yu = tree->a_nodes[j]->v[k]->ldsk->veloc->deriv[i];
                          tu = fabs(tree->times->nd_t[tree->a_nodes[j]->v[k]->num] - tree->times->nd_t[tree->a_nodes[j]->num]);
                        }
                      else
                        {
                          dirv = k;
                          xv = tree->a_nodes[j]->v[k]->ldsk->coord->lonlat[i];
                          yv = tree->a_nodes[j]->v[k]->ldsk->veloc->deriv[i];
                          tv = fabs(tree->times->nd_t[tree->a_nodes[j]->v[k]->num] - tree->times->nd_t[tree->a_nodes[j]->num]);
                        }
                    }
                  else
                    {
                      xw = tree->a_nodes[j]->v[k]->ldsk->coord->lonlat[i];
                      yw = tree->a_nodes[j]->v[k]->ldsk->veloc->deriv[i];
                      tz = fabs(tree->times->nd_t[tree->a_nodes[j]->v[k]->num] - tree->times->nd_t[tree->a_nodes[j]->num]);
                    }
                }
              
              assert(tu > 0.);
              assert(tv > 0.);
              assert(tz > 0.);
              
              muu = xu - (yu+yz)/2.*tu;
              
              sigsqu =
                log(tree->mmod->sigsq[i]) +
                log(tree->mmod->sigsq_scale_norm_fact) +
                3.*log(tu)-
                LOG12;
              sigsqu = exp(sigsqu);
              
              
              muv = xv - (yv+yz)/2.*tv;
              
              sigsqv =
                log(tree->mmod->sigsq[i]) +
                log(tree->mmod->sigsq_scale_norm_fact) +
                3.*log(tv)-
                LOG12;
              sigsqv = exp(sigsqv);
              
              
              
              muz = xw + (yw+yz)/2.*tz;;
              
              sigsqz =
                log(tree->mmod->sigsq[i]) +
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
              
              
              
              assert(isinf(mean) == NO && isnan(mean) == NO);
              assert(isinf(var) == NO && isnan(var) == NO);
              
              phydbl dum = tree->a_nodes[j]->ldsk->coord->lonlat[i];
              tree->a_nodes[j]->ldsk->coord->lonlat[i] = Rnorm(mean,sqrt(var));
              
              /* PhyML_Printf("\n. LOC Node %d xu: %f yu: %f tu: %f xv: %f yv: %f tv: %f xw: %f yw: %f tz: %f -- muu: %f sigsqu: %f | muv: %f sigsqv: %f | muz: %f sigsqz: %f mean: %f var: %f [sigsq: %f] --> %f", */
              /*                  j, */
              /*                  xu,yu,tu, */
              /*                  xv,yv,tv, */
              /*                  xw,yw,tz, */
              /*                  muu,sigsqu, */
              /*                  muv,sigsqv, */
              /*                  muz,sigsqz, */
              /*                  mean,var, */
              /*                  tree->mmod->sigsq[i],tree->mmod->c_lnL); */
              /*     PhyML_Printf(" %f --> %f",dum,tree->a_nodes[j]->ldsk->coord->lonlat[i]); */
            
            }
        }
      

      /* Root node */

      yz = tree->n_root->ldsk->veloc->deriv[i];

      xu = tree->n_root->v[1]->ldsk->coord->lonlat[i];
      yu = tree->n_root->v[1]->ldsk->veloc->deriv[i];
      tu = fabs(tree->times->nd_t[tree->n_root->v[1]->num] - tree->times->nd_t[tree->n_root->num]); 

      xv = tree->n_root->v[2]->ldsk->coord->lonlat[i];
      yv = tree->n_root->v[2]->ldsk->veloc->deriv[i];
      tv = fabs(tree->times->nd_t[tree->n_root->v[2]->num] - tree->times->nd_t[tree->n_root->num]); 

      
      assert(tu > 0.0);

      muu = xu - (yu+yz)/2.*tu;

      sigsqu =
        log(tree->mmod->sigsq[i]) +
        log(tree->mmod->sigsq_scale_norm_fact) +
        3.*log(tu)-
        LOG12;
      sigsqu = exp(sigsqu);
      
      
      muv = xv - (yv+yz)/2.*tv;
      
      sigsqv = 
        log(tree->mmod->sigsq[i]) +
        log(tree->mmod->sigsq_scale_norm_fact) +
        3.*log(tv)-
        LOG12;
      sigsqv = exp(sigsqv);
      
      assert(sigsqu > 0.0);
      assert(sigsqv > 0.0);
      
      var = 1./sigsqu + 1./sigsqv;
      var = 1./var;

      mean = (muu/sigsqu + muv/sigsqv)*var;
      
      assert(isinf(mean) == NO && isnan(mean) == NO);
      assert(isinf(var) == NO && isnan(var) == NO);
      
      tree->n_root->ldsk->coord->lonlat[i] = Rnorm(mean,sqrt(var));

      /* PhyML_Printf("\n. LOC yz: %f Root xu: %f yu: %f tu: %f xv: %f yv: %f tv: %f -- muu: %f sigsqu: %f | muv: %f sigsqv: %f | mean: %f var: %f --> %f", */
      /*              yz, */
      /*              xu,yu,tu, */
      /*              xv,yv,tv, */
      /*              muu,sigsqu, */
      /*              muv,sigsqv, */
      /*              mean,var, */
      /*              LOCATION_Lk(tree)); */
      /* PhyML_Printf(" --> %f",tree->n_root->ldsk->coord->lonlat[i]); */
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
  
/* Returns log[Pr(locations at tips | velocities at all nodes, tree)] */
phydbl IBM_Lk_Locations(t_tree *tree)
{  
  /* phydbl lnL; */
  /* int i; */
  
  /* tree->mmod->sigsq_scale_norm_fact = 1.0; */

  /* lnL = 0.0; */

  /* for(i=0;i<tree->mmod->n_dim;++i) */
  /*   { */
  /*     IBM_Init_Contmod_Locations(i,tree); */
  /*     IBM_Lk_Locations_Post(NULL,tree->n_root,tree->mmod->sigsq[i],i,tree,NO); */
  /*     lnL += tree->contmod->logrem_down[tree->n_root->num]; */
  /*   } */
  /* return(lnL); */
  
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
  int i;

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

          mean = 0.5*(3.*(ldsk->coord->lonlat[i] - ldsk->prev->coord->lonlat[i])/dt - ldsk->prev->veloc->deriv[i]); 

          sd =
            log(tree->mmod->sigsq[i]) +
            log(tree->mmod->sigsq_scale[nd_d->num]) +
            log(dt)-
            LOG4;

          sd = sqrt(exp(sd));

          
          ld = ldsk->veloc->deriv[i];

          disk_lnP += Log_Dnorm(ld,mean,sd,&err);
                    
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
          
          /* mean = ldsk->prev->coord->lonlat[i] + 0.5*dt*(ldsk->veloc->deriv[i] + ldsk->prev->veloc->deriv[i]); */
          mean = ldsk->prev->coord->lonlat[i] + dt*ldsk->prev->veloc->deriv[i];

          /* sd = */
          /*   log(tree->mmod->sigsq[i]) + */
          /*   log(tree->mmod->sigsq_scale[nd_d->num]) + */
          /*   3.*log(dt)- */
          /*   LOG12; */
          sd =
            log(tree->mmod->sigsq[i]) +
            log(tree->mmod->sigsq_scale[nd_d->num]) +
            3.*log(dt)-
            LOG3;

          sd = sqrt(exp(sd));

          
          ld = ldsk->coord->lonlat[i];

          disk_lnP += Log_Dnorm(ld,mean,sd,&err);
                    
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

