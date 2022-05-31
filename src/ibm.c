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
  
/* Returns log[Pr(locations at tips | velocities at all nodes, tree)] */
phydbl IBM_Lk_Locations(t_tree *tree)
{  
  phydbl lnL;
  int i;
  
  lnL = 0.0;

  for(i=0;i<tree->mmod->n_dim;++i)
    {
      RRW_Init_Contmod_Locations(i,tree);
      IBM_Lk_Locations_Post(NULL,tree->n_root,tree->mmod->sigsq[i],i,tree,NO);
      lnL += tree->contmod->logrem_down[tree->n_root->num];
    }
  return(lnL);
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
          if(d->v[i] != a && d->b[i] != tree->e_root)
            {
              IBM_Lk_Locations_Post(d,d->v[i],sigsq,dim,tree,print);
            }
        }

      v1 = v2 = NULL;
      for(i=0;i<3;++i)
        {
          if(d->v[i] != a && d->b[i] != tree->e_root)
            {
              if(v1 == NULL) v1 = d->v[i];
              else v2 = d->v[i];
            }
        }

      v1mu =
        tree->contmod->mu_down[v1->num] -
        (v1->ldsk->veloc->deriv[dim] + d->ldsk->veloc->deriv[dim])/2. *
        fabs(tree->times->nd_t[v1->num]-tree->times->nd_t[d->num]);
      
      v2mu =
        tree->contmod->mu_down[v2->num] -
        (v2->ldsk->veloc->deriv[dim] + d->ldsk->veloc->deriv[dim])/2. *
        fabs(tree->times->nd_t[v2->num]-tree->times->nd_t[d->num]);

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
        -log(12.) +
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

void IBM_Sample_Velocities(t_tree *tree)
{
  int i,j,print;
  phydbl mean,sd,dt,zi;
  
  
  print = NO;
  
  /* Sample ancestral velocities given locations and velocities at tips */
  for(i=0;i<tree->mmod->n_dim;++i)
    {
      IBM_Init_Contmod_Velocities(i,tree);
      tree->mmod->sigsq_scale_norm_fact = 1.0;
      Lk_Contmod_Post(NULL,tree->n_root,tree->mmod->sigsq[i],tree,NO);
      Lk_Contmod_Pre(tree->n_root,tree->n_root->v[1],tree->mmod->sigsq[i],tree);
      Lk_Contmod_Pre(tree->n_root,tree->n_root->v[2],tree->mmod->sigsq[i],tree);      
      
      for(j=0;j<2*tree->n_otu-2;++j)
        {
          if(tree->a_nodes[j]->tax == NO && tree->a_nodes[j] != tree->n_root)
            {
              tree->a_nodes[j]->ldsk->veloc->deriv[i] =
                Sample_Ancestral_Trait_Contmod(tree->a_nodes[j]->anc,
                                               tree->a_nodes[j],
                                               fabs(tree->times->nd_t[tree->a_nodes[j]->anc->num]-tree->times->nd_t[tree->a_nodes[j]->num]),
                                               0.0,
                                               log(tree->mmod->sigsq[i]) +
                                               log(tree->mmod->sigsq_scale[tree->a_nodes[j]->num]) +
                                               log(tree->mmod->sigsq_scale_norm_fact),
                                               NO,
                                               tree);
            }
        }
      
      assert(isnan(tree->contmod->var_down[tree->n_root->num]) == NO);
      
      tree->n_root->ldsk->veloc->deriv[i] =
        Rnorm(tree->contmod->mu_down[tree->n_root->num],
              sqrt(tree->contmod->var_down[tree->n_root->num]));
    }


  /* Sample velocities at tips given ancestral velocities */
  for(i=0;i<tree->mmod->n_dim;++i)
    {      
      for(j=0;j<2*tree->n_otu-2;++j)
        {
          dt = fabs(tree->times->nd_t[tree->a_nodes[j]->v[0]->num]-tree->times->nd_t[tree->a_nodes[j]->num]);
          zi = tree->a_nodes[j]->ldsk->coord->lonlat[i] - tree->a_nodes[j]->v[0]->ldsk->coord->lonlat[i];
          
          mean = 0.5*(3.*zi/dt + tree->a_nodes[j]->v[0]->ldsk->veloc->deriv[i]);

          sd =
            -log(4.) +
            log(dt) + 
            log(tree->mmod->sigsq[i]) +
            log(tree->mmod->sigsq_scale_norm_fact) +
            log(tree->mmod->sigsq_scale[tree->a_nodes[j]->num]);
          
          sd = exp(sd);
          sd = sqrt(sd);
          
          tree->a_nodes[j]->ldsk->veloc->deriv[i] = Rnorm(mean,sd);
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
            log(dt);

          sd = sqrt(exp(sd));

          
          ld = ldsk->veloc->deriv[i];
          la = ldsk->prev->veloc->deriv[i];

          disk_lnP += Log_Dnorm(ld,la,sd,&err);
          
          if(isinf(lnP) || isnan(lnP)) return(UNLIKELY);
          
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

