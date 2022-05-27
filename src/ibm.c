/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include "ibm.h"


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

phydbl IBM_Lk_Veloc(t_tree *tree)
{
  phydbl lnL;
  int i;
  
  lnL = 0.0;

  for(i=0;i<tree->mmod->n_dim;++i)
    {
      RRW_Init_Contmod_Velocities(i,tree);
      Lk_Contmod_Post(NULL,tree->n_root,tree->mmod->sigsq[i],tree,NO);
      lnL += tree->contmod->logrem_down[tree->n_root->num];
    }
  return(lnL);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

