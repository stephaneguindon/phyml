/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/


/* Routines for molecular clock trees and molecular dating */


#include "times.h"

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void TIMES_Least_Square_Node_Times(t_edge *e_root, t_tree *tree)
{

  /* Solve A.x = b, where x is the vector of times estimated
     under the least square criterion and b is the vector
     of edge lengths

     A is a n x n matrix, with n being the number of
     nodes in a rooted tree (i.e. 2*n_otu-1).

   */

  phydbl *A, *b, *x;
  int n;
  int i,j;
  t_node *root;

  n = 2*tree->n_otu-1;
  
  A = (phydbl *)mCalloc(n*n,sizeof(phydbl));
  b = (phydbl *)mCalloc(n,  sizeof(phydbl));
  x = (phydbl *)mCalloc(n,  sizeof(phydbl));
    
  if(!tree->n_root && e_root) Add_Root(e_root,tree);
  else if(!e_root)            Add_Root(tree->a_edges[0],tree);
  
  root = tree->n_root;

  TIMES_Least_Square_Node_Times_Pre(root,root->v[1],A,b,n,tree);
  TIMES_Least_Square_Node_Times_Pre(root,root->v[2],A,b,n,tree);
  
  b[root->num] = tree->e_root->l->v/2.;
  
  A[root->num * n + root->num]       = 1.0;
  A[root->num * n + root->v[2]->num] = -.5;
  A[root->num * n + root->v[1]->num] = -.5;
    
  if(!Matinv(A, n, n,YES))
    {
      PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d (function '%s').\n",__FILE__,__LINE__,__FUNCTION__);
      Exit("\n");      
    }

  for(i=0;i<n;i++) x[i] = .0;
  for(i=0;i<n;i++) for(j=0;j<n;j++) x[i] += A[i*n+j] * b[j];

  for(i=0;i<n-1;i++) tree->times->nd_t[tree->a_nodes[i]->num] = -x[i];
  tree->times->nd_t[root->num] = -x[n-1];
  tree->n_root->b[2]->l->v = tree->times->nd_t[root->v[2]->num] - tree->times->nd_t[root->num];
  tree->n_root->b[1]->l->v = tree->times->nd_t[root->v[1]->num] - tree->times->nd_t[root->num];

  Free(A);
  Free(b);
  Free(x);

  return;
  
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void TIMES_Least_Square_Node_Times_Pre(t_node *a, t_node *d, phydbl *A, phydbl *b, int n, t_tree *tree)
{
  if(d->tax)
    {
      A[d->num * n + d->num] = 1.;
      
      /* Set the time stamp at tip nodes to 0.0 */
/*       PhyML_Printf("\n. Tip t_node date set to 0"); */
      b[d->num] = 0.0;
      return;
    }
  else
    {
      int i;
       
      for(i=0;i<3;++i)
	if((d->v[i] != a) && (d->b[i] != tree->e_root))
	  TIMES_Least_Square_Node_Times_Pre(d,d->v[i],A,b,n,tree);
      
      A[d->num * n + d->num] = 1.;
      b[d->num] = .0;
      for(i=0;i<3;++i)
	{
	  A[d->num * n + d->v[i]->num] = -1./3.;
	  if(d->v[i] != a) b[d->num] += d->b[i]->l->v;
	  else             b[d->num] -= d->b[i]->l->v;
	}
      b[d->num] /= 3.;
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


/* Adjust t_node times in order to have correct time stamp ranking with
 respect to the tree topology */

void TIMES_Adjust_Node_Times(t_tree *tree)
{
  TIMES_Adjust_Node_Times_Pre(tree->n_root->v[2],tree->n_root->v[1],tree);
  TIMES_Adjust_Node_Times_Pre(tree->n_root->v[1],tree->n_root->v[2],tree);

  if(tree->times->nd_t[tree->n_root->num] > MIN(tree->times->nd_t[tree->n_root->v[2]->num],
						tree->times->nd_t[tree->n_root->v[1]->num]))
    {
      tree->times->nd_t[tree->n_root->num] = MIN(tree->times->nd_t[tree->n_root->v[2]->num],
						 tree->times->nd_t[tree->n_root->v[1]->num]);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void TIMES_Adjust_Node_Times_Pre(t_node *a, t_node *d, t_tree *tree)
{
  if(d->tax) return;
  else
    {
      int i;
      phydbl min_height;

      for(i=0;i<3;i++)
	if((d->v[i] != a) && (d->b[i] != tree->e_root))
	  {
	    TIMES_Adjust_Node_Times_Pre(d,d->v[i],tree);
	  }

      min_height = 0.0;
      for(i=0;i<3;i++)
	{
	  if((d->v[i] != a) && (d->b[i] != tree->e_root))
	    {
	      if(tree->times->nd_t[d->v[i]->num] < min_height)
		{
		  min_height = tree->times->nd_t[d->v[i]->num];
		}
	    }
	}

      if(tree->times->nd_t[d->num] > min_height) tree->times->nd_t[d->num] = min_height;

      if(tree->times->nd_t[d->num] < -100.) tree->times->nd_t[d->num] = -100.;

    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


  /* Multiply each time stamp at each internal 
     t_node by  'tree->time_stamp_mult'.
   */

void TIMES_Mult_Time_Stamps(t_tree *tree)
{
  int i;
  For(i,2*tree->n_otu-2) tree->times->nd_t[tree->a_nodes[i]->num] *= FABS(tree->mod->s_opt->tree_size_mult);
  tree->times->nd_t[tree->n_root->num] *= FABS(tree->mod->s_opt->tree_size_mult);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void TIMES_Print_Node_Times(t_node *a, t_node *d, t_tree *tree)
{
  t_edge *b;
  int i;
  
  b = NULL;
  for(i=0;i<3;i++) if((d->v[i]) && (d->v[i] == a)) {b = d->b[i]; break;}

  PhyML_Printf("\n. (%3d %3d) a->t = %12f d->t = %12f (#=%12f) b->l->v = %12f [%12f;%12f]",
	       a->num,d->num,
	       tree->times->nd_t[a->num],
	       tree->times->nd_t[d->num],
	       tree->times->nd_t[a->num]-tree->times->nd_t[d->num],
	       (b)?(b->l->v):(-1.0),
	       tree->times->t_prior_min[d->num],
	       tree->times->t_prior_max[d->num]);
  if(d->tax) return;
  else
    {
      int i;
      for(i=0;i<3;i++)
	if((d->v[i] != a) && (d->b[i] != tree->e_root))
	  TIMES_Print_Node_Times(d,d->v[i],tree);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void TIMES_Set_All_Node_Priors(t_tree *tree)
{
  int i;
  phydbl min_prior;

  /* Set all t_prior_max values */
  TIMES_Set_All_Node_Priors_Bottom_Up(tree->n_root,tree->n_root->v[2],tree);
  TIMES_Set_All_Node_Priors_Bottom_Up(tree->n_root,tree->n_root->v[1],tree);

  tree->times->t_prior_max[tree->n_root->num] = 
    MIN(tree->times->t_prior_max[tree->n_root->num],
	MIN(tree->times->t_prior_max[tree->n_root->v[2]->num],
	    tree->times->t_prior_max[tree->n_root->v[1]->num]));


  /* Set all t_prior_min values */
  if(!tree->times->t_has_prior[tree->n_root->num])
    {
      min_prior = 1.E+10;
      for(i=0;i<2*tree->n_otu-2;++i)
	{
	  if(tree->times->t_has_prior[i])
	    {
	      if(tree->times->t_prior_min[i] < min_prior)
		min_prior = tree->times->t_prior_min[i];
	    }
	}
      tree->times->t_prior_min[tree->n_root->num] = 2.0 * min_prior;
      /* tree->times->t_prior_min[tree->n_root->num] = 10. * min_prior; */
    }
  
  if(tree->times->t_prior_min[tree->n_root->num] > 0.0)
    {
      PhyML_Fprintf(stderr,"\n. Failed to set the lower bound for the root node.");
      PhyML_Fprintf(stderr,"\n. Make sure at least one of the calibration interval");
      PhyML_Fprintf(stderr,"\n. provides a lower bound.");
      Exit("\n");
    }


  TIMES_Set_All_Node_Priors_Top_Down(tree->n_root,tree->n_root->v[2],tree);
  TIMES_Set_All_Node_Priors_Top_Down(tree->n_root,tree->n_root->v[1],tree);

  Get_Node_Ranks(tree);
  TIMES_Set_Floor(tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void TIMES_Set_All_Node_Priors_Bottom_Up(t_node *a, t_node *d, t_tree *tree)
{
  int i;
  phydbl t_sup;
  
  if(d->tax) return;
  else 
    {
      t_node *v1, *v2; /* the two sons of d */

      for(i=0;i<3;i++)
	{
	  if((d->v[i] != a) && (d->b[i] != tree->e_root))
	    {
	      TIMES_Set_All_Node_Priors_Bottom_Up(d,d->v[i],tree);	      
	    }
	}
      
      v1 = v2 = NULL;
      for(i=0;i<3;i++) if((d->v[i] != a) && (d->b[i] != tree->e_root)) 
	{
	  if(!v1) v1 = d->v[i]; 
	  else    v2 = d->v[i];
	}
      
      if(tree->times->t_has_prior[d->num] == YES)
	{
	  t_sup = MIN(tree->times->t_prior_max[d->num],
		      MIN(tree->times->t_prior_max[v1->num],
			  tree->times->t_prior_max[v2->num]));

	  tree->times->t_prior_max[d->num] = t_sup;

	  if(tree->times->t_prior_max[d->num] < tree->times->t_prior_min[d->num])
	    {
	      PhyML_Fprintf(stderr,"\n. prior_min=%f prior_max=%f",tree->times->t_prior_min[d->num],tree->times->t_prior_max[d->num]);
	      PhyML_Fprintf(stderr,"\n. Inconsistency in the prior settings detected at node %d",d->num);
	      PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d (function %s)\n\n",__FILE__,__LINE__,__FUNCTION__);
	      Warn_And_Exit("\n");
	    }
	}
      else
	{
	  tree->times->t_prior_max[d->num] = 
	    MIN(tree->times->t_prior_max[v1->num],
		tree->times->t_prior_max[v2->num]);
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void TIMES_Set_All_Node_Priors_Top_Down(t_node *a, t_node *d, t_tree *tree)
{
  if(d->tax) return;
  else
    {
      int i;      
      
      if(tree->times->t_has_prior[d->num] == YES)
	{
	  tree->times->t_prior_min[d->num] = MAX(tree->times->t_prior_min[d->num],tree->times->t_prior_min[a->num]);
	  
	  if(tree->times->t_prior_max[d->num] < tree->times->t_prior_min[d->num])
	    {
	      PhyML_Fprintf(stderr,"\n. prior_min=%f prior_max=%f",tree->times->t_prior_min[d->num],tree->times->t_prior_max[d->num]);
	      PhyML_Fprintf(stderr,"\n. Inconsistency in the prior settings detected at t_node %d",d->num);
	      PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d (function %s)\n\n",__FILE__,__LINE__,__FUNCTION__);
	      Warn_And_Exit("\n");
	    }
	}
      else
	{
	  tree->times->t_prior_min[d->num] = tree->times->t_prior_min[a->num];
	}
            
      for(i=0;i<3;i++)
	{
	  if((d->v[i] != a) && (d->b[i] != tree->e_root))
	    {
	      TIMES_Set_All_Node_Priors_Top_Down(d,d->v[i],tree);
	    }
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void TIMES_Set_Floor(t_tree *tree)
{
  TIMES_Set_Floor_Post(tree->n_root,tree->n_root->v[2],tree);
  TIMES_Set_Floor_Post(tree->n_root,tree->n_root->v[1],tree);
  tree->times->t_floor[tree->n_root->num] = MIN(tree->times->t_floor[tree->n_root->v[2]->num],
						tree->times->t_floor[tree->n_root->v[1]->num]);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void TIMES_Set_Floor_Post(t_node *a, t_node *d, t_tree *tree)
{
  if(d->tax)
    {
      tree->times->t_floor[d->num] = tree->times->nd_t[d->num];
      d->rank_max = d->rank;
      return;
    }
  else
    {
      int i;
      t_node *v1,*v2;

      v1 = v2 = NULL;
      for(i=0;i<3;i++)
	{
	  if(d->v[i] != a && d->b[i] != tree->e_root)
	    {
	      TIMES_Set_Floor_Post(d,d->v[i],tree);
	      if(!v1) v1 = d->v[i];
	      else    v2 = d->v[i];
	    }
	}
      tree->times->t_floor[d->num] = MIN(tree->times->t_floor[v1->num],
					 tree->times->t_floor[v2->num]);

      if(tree->times->t_floor[v1->num] < tree->times->t_floor[v2->num])
	{
	  d->rank_max = v1->rank_max;
	}
      else if(tree->times->t_floor[v2->num] < tree->times->t_floor[v1->num])
	{
	  d->rank_max = v2->rank_max;
	}
      else
	{
	  d->rank_max = MAX(v1->rank_max,v2->rank_max);
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* Does it work for serial samples? */
phydbl TIMES_Log_Conditional_Uniform_Density(t_tree *tree)
{
  phydbl min,max;
  phydbl dens;
  int i;

  min = tree->times->nd_t[tree->n_root->num];

  dens = 0.0;
  For(i,2*tree->n_otu-1)
    {
      if((tree->a_nodes[i]->tax == NO) && (tree->a_nodes[i] != tree->n_root))
	{
	  max = tree->times->t_floor[i];

	  dens += log(Dorder_Unif(tree->times->nd_t[i],
				  tree->a_nodes[i]->rank-1,
				  tree->a_nodes[i]->rank_max-2,
				  min,max));
	}
    }
  return dens;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// Returns the marginal density of tree height assuming the
// Yule model of speciation. 
phydbl TIMES_Lk_Yule_Root_Marginal(t_tree *tree)
{
  int n;
  int j;
  t_node *nd;
  phydbl *t,*ts;
  phydbl lbda;
  phydbl T;

  lbda = tree->times->birth_rate;
  t    = tree->times->nd_t;
  ts   = tree->times->time_slice_lims;
  T    = ts[0] - t[tree->n_root->num];

  n = 0;
  nd = NULL;
  For(j,2*tree->n_otu-2) 
    {
      nd = tree->a_nodes[j];

      if((t[nd->num] > ts[0] && t[nd->anc->num] < ts[0]) || // lineage that is crossing ts[0]
	 (nd->tax == YES && Are_Equal(t[nd->num],ts[0],1.E-6) == YES)) // tip that is lying on ts[0]
	n++;
    }

  return LnGamma(n+1) + log(lbda) - 2.*lbda*T + (n-2.)*log(1. - exp(-lbda*T));
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// Returns the joint density of internal node heights assuming
// the Yule model of speciation.
phydbl TIMES_Lk_Yule_Joint(t_tree *tree)
{
  int i,j;
  phydbl loglk;
  phydbl *t;
  phydbl dt;
  int n; // number of lineages at a given time point
  phydbl lbda;
  t_node *nd;
  phydbl *ts;
  int *tr;
  phydbl top_t;
  short int *interrupted;
  phydbl sumdt;

  interrupted = (short int *)mCalloc(tree->n_otu,sizeof(short int));

  t = tree->times->nd_t;
  ts = tree->times->time_slice_lims;
  tr = tree->times->t_rank;
  lbda = tree->times->birth_rate;

  TIMES_Update_Node_Ordering(tree);

  for(j=0;j<tree->n_otu;j++) interrupted[j] = NO;

  loglk = .0;

  sumdt = .0;
  n = 1;
  For(i,2*tree->n_otu-2) // t[tr[0]] is the time of the oldest node, t[tr[1]], the second oldest and so on...
    {

      for(j=0;j<tree->n_otu;j++)
	if((t[j] < t[tr[i]]) && (interrupted[j] == NO)) 
	  {
	    interrupted[j] = YES;
	    n--; // How many lineages have stopped above t[tr[i]]?
	  }
      
      top_t = t[tr[i+1]];
      dt = top_t - t[tr[i]];
      sumdt += dt;

      /* printf("\n. %d node up=%d [%f] node do=%d [%f] dt=%f",i,tr[i],t[tr[i]],tr[i+1],t[tr[i+1]],dt); */

      if(n<1)
	{
	  PhyML_Fprintf(stderr,"\n. i=%d tr[i]=%f",i,t[tr[i]]);
	  PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d\n",__FILE__,__LINE__);
	  Exit("\n");
	}

      if(dt > 1.E-10) loglk += log((n+1)*lbda) - (n+1)*lbda*dt;
      n++;      
    }

  /* printf("\n. sumdt = %f th=%f",sumdt,tree->times->nd_t[tree->n_root->num]); */
  /* printf("\n0 loglk = %f",loglk); */

  for(i=0;i<tree->times->n_time_slices-1;i++)
    {
      n = 0;
      dt = 0.;
      For(j,2*tree->n_otu-2)
  	{
  	  nd = tree->a_nodes[j];
  	  if(t[nd->num] > ts[i] && t[nd->anc->num] < ts[i]) // How many lineages are crossing this time slice limit?
  	    {
  	      n++;
  	      if(t[nd->num] < dt) dt = t[nd->num]; // take the oldest node younger than the time slice
  	    }
  	}
      dt -= ts[i];
      loglk += log(n*lbda) - n*lbda*dt;
    }

  /* printf("\n1 loglk = %f",loglk); */

  Free(interrupted);

  return loglk;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// Returns the conditional density of internal node heights 
// given the tree height under the Yule model. Uses the order
// statistics 'simplification' as described in Yang and Rannala, 2005. 
phydbl TIMES_Lk_Yule_Order(t_tree *tree)
{
  int j;
  phydbl *t,*tf;
  t_node *n;
  phydbl loglk;
  phydbl loglbda;
  phydbl lbda;
  phydbl *tp_min,*tp_max;
  phydbl lower_bound,upper_bound;
  /* phydbl root_height; */

  tp_min = tree->times->t_prior_min;
  tp_max = tree->times->t_prior_max;
  tf = tree->times->t_floor;
  t  = tree->times->nd_t;
  n = NULL;
  loglbda = log(tree->times->birth_rate);
  lbda = tree->times->birth_rate;
  lower_bound = -1.;
  upper_bound = -1.;
  /* root_height = FABS(tree->times->nd_t[tree->n_root->num]); */

  /*! Adapted from  Equation (6) in T. Stadler's Systematic Biology, 2012 paper with
      sampling fraction set to 1 and death rate set to 0. Dropped the 1/(n-1) scaling 
      factor. */

  /* loglk = 0.0; */
  /* For(j,2*tree->n_otu-2) */
  /*   { */
  /*     n = tree->a_nodes[j]; */
  /*     lower_bound = MAX(FABS(tf[j]),FABS(tp_max[j])); */
  /*     upper_bound = MIN(FABS(t[tree->n_root->num]),FABS(tp_min[j])); */

  /*     if(n->tax == NO) */
  /*       { */
  /*         loglk  += (loglbda - lbda * FABS(t[j])); */
  /*         /\* loglk -= log(exp(-lbda*lower_bound) - exp(-lbda*upper_bound)); // incorporate calibration boundaries here. *\/ */
  /*       } */
  /*   } */

  
  /*! Adapted from  Equation (7) in T. Stadler's Systematic Biology, 2012 paper with
      sampling fraction set to 1 and death rate set to 0. */

  // Check that each node is within its calibration-derived interval
  For(j,2*tree->n_otu-1) if(t[j] < tp_min[j] || t[j] > tp_max[j]) return(-INFINITY);

  loglk = 0.0;
  For(j,2*tree->n_otu-2)
    {
      n = tree->a_nodes[j];
      lower_bound = MAX(FABS(tf[j]),FABS(tp_max[j]));
      upper_bound = FABS(tp_min[j]);
      
      if(n->tax == NO)
        {
          loglk  += (loglbda - lbda * FABS(t[j]));
          loglk -= log(exp(-lbda*lower_bound) - exp(-lbda*upper_bound)); // incorporate calibration boundaries here.    
        }
      
      if(isinf(loglk) || isnan(loglk))
        {
          /* PhyML_Printf("\n. Lower bound: %f",lower_bound); */
          /* PhyML_Printf("\n. Upper bound: %f",upper_bound); */
          /* PhyML_Printf("\n. tf: %f tp_max: %f tp_min: %f ",tf[j],tp_max[j],tp_min[j]); */
          /* PhyML_Printf("\n. exp1: %f",exp(-lbda*lower_bound)); */
          /* PhyML_Printf("\n. exp2: %f",exp(-lbda*upper_bound)); */
          /* PhyML_Printf("\n. diff: %f",exp(-lbda*lower_bound) - exp(-lbda*upper_bound)); */
          /* Exit("\n"); */
          return(-INFINITY);
        }      
    }

  lower_bound = MAX(FABS(tf[tree->n_root->num]),FABS(tp_max[tree->n_root->num]));
  upper_bound = FABS(tp_min[tree->n_root->num]);
  loglk += log(2) + loglbda - 2.*lbda * FABS(t[tree->n_root->num]);
  loglk -= log(exp(-2.*lbda*lower_bound) - exp(-2.*lbda*upper_bound));

  if(isinf(loglk) || isnan(loglk))
    {
      /* PhyML_Printf("\n. * Lower bound: %f",lower_bound); */
      /* PhyML_Printf("\n. * Upper bound: %f",upper_bound); */
      /* PhyML_Printf("\n. * tf: %f tp_max: %f tp_min: %f",tf[tree->n_root->num],tp_max[tree->n_root->num],tp_min[tree->n_root->num]); */
      /* PhyML_Printf("\n. * exp1: %f",exp(-2.*lbda*lower_bound)); */
      /* PhyML_Printf("\n. * exp2: %f",exp(-2.*lbda*upper_bound)); */
      /* PhyML_Printf("\n. * diff: %f",exp(-2.*lbda*lower_bound) - exp(-2.*lbda*upper_bound)); */
      /* Exit("\n"); */
      return(-INFINITY);
    }


  return(loglk);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void TIMES_Lk_Times_Trav(t_node *a, t_node *d, phydbl lim_inf, phydbl lim_sup, phydbl *logdens, t_tree *tree)
{
  int i;
  
  if(!d->tax)
    {
      /* if(tree->times->nd_t[d->num] > lim_sup) */
      /* 	{ */
      /* 	  lim_inf = lim_sup; */
      /* 	  lim_sup = 0.0; */
      /* 	  For(i,2*tree->n_otu-2) */
      /* 	    if((tree->times->t_floor[i] < lim_sup) && (tree->times->t_floor[i] > tree->times->nd_t[d->num])) */
      /* 	      lim_sup = tree->times->t_floor[i]; */
      /* 	} */
      
      /* if(tree->times->nd_t[d->num] < lim_inf || tree->times->nd_t[d->num] > lim_sup) */
      /* 	{ */
      /* 	  PhyML_Printf("\n. nd_t = %f lim_inf = %f lim_sup = %f",tree->times->nd_t[d->num],lim_inf,lim_sup); */
      /* 	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__); */
      /* 	  Exit("\n"); */
      /* 	} */
  
      lim_inf = tree->times->nd_t[tree->n_root->num];
      lim_sup = tree->times->t_floor[d->num];
      
      *logdens = *logdens + log(lim_sup - lim_inf);   
    }
  
  if(d->tax == YES) return;
  else
    {      
      for(i=0;i<3;i++)
	{
	  if(d->v[i] != a && d->b[i] != tree->e_root)
	    {
	      TIMES_Lk_Times_Trav(d,d->v[i],lim_inf,lim_sup,logdens,tree);
	    }
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl TIMES_Log_Number_Of_Ranked_Labelled_Histories(t_node *root, int per_slice, t_tree *tree)
{
  int i;
  phydbl logn;
  t_node *v1,*v2;
  int n1,n2;
  
  TIMES_Update_Curr_Slice(tree);

  logn = .0;
  v1 = v2 = NULL;
  if(root == tree->n_root)
    {
      TIMES_Log_Number_Of_Ranked_Labelled_Histories_Post(root,root->v[2],per_slice,&logn,tree);
      TIMES_Log_Number_Of_Ranked_Labelled_Histories_Post(root,root->v[1],per_slice,&logn,tree);
      v1 = root->v[2];
      v2 = root->v[1];
    }
  else
    {
      for(i=0;i<3;i++)
	{
	  if(root->v[i] != root->anc && root->b[i] != tree->e_root)
	    {
	      TIMES_Log_Number_Of_Ranked_Labelled_Histories_Post(root,root->v[i],per_slice,&logn,tree);
	      if(!v1) v1 = root->v[i];
	      else    v2 = root->v[i];
	    }
	}
    }

 
  if(per_slice == NO)
    {
      n1 = tree->rates->n_tips_below[v1->num];
      n2 = tree->rates->n_tips_below[v2->num];
    }
  else
    {
      if(tree->times->curr_slice[v1->num] == tree->times->curr_slice[root->num])
  	n1 = tree->rates->n_tips_below[v1->num];
      else
  	n1 = 1;
      
      if(tree->times->curr_slice[v2->num] == tree->times->curr_slice[root->num])
  	n2 = tree->rates->n_tips_below[v2->num];
      else
  	n2 = 1;
    }

  tree->rates->n_tips_below[root->num] = n1+n2;

  logn += Factln(n1+n2-2) - Factln(n1-1) - Factln(n2-1);

  return(logn);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void TIMES_Log_Number_Of_Ranked_Labelled_Histories_Post(t_node *a, t_node *d, int per_slice, phydbl *logn, t_tree *tree)
{
  if(d->tax == YES) 
    {
      tree->rates->n_tips_below[d->num] = 1;
      return;
    }
  else
    {
      int i,n1,n2;
      t_node *v1, *v2;

      for(i=0;i<3;i++)
	{
	  if(d->v[i] != a && d->b[i] != tree->e_root)
	    {
	      TIMES_Log_Number_Of_Ranked_Labelled_Histories_Post(d,d->v[i],per_slice,logn,tree);
	    }
	}

      v1 = v2 = NULL;
      for(i=0;i<3;i++)
	{
	  if(d->v[i] != a && d->b[i] != tree->e_root)
	    {
	      if(v1 == NULL) {v1 = d->v[i];}
	      else           {v2 = d->v[i];}
	    }
	}


      if(per_slice == NO)
	{
	  n1 = tree->rates->n_tips_below[v1->num];
	  n2 = tree->rates->n_tips_below[v2->num];
	}
      else
	{
	  if(tree->times->curr_slice[v1->num] == tree->times->curr_slice[d->num])
	    n1 = tree->rates->n_tips_below[v1->num];
	  else
	    n1 = 1;

	  if(tree->times->curr_slice[v2->num] == tree->times->curr_slice[d->num])
	    n2 = tree->rates->n_tips_below[v2->num];
	  else
	    n2 = 1;
	}

      tree->rates->n_tips_below[d->num] = n1+n2;

      (*logn) += Factln(n1+n2-2) - Factln(n1-1) - Factln(n2-1);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void TIMES_Update_Curr_Slice(t_tree *tree)
{
  int i,j;

  for(i=0;i<2*tree->n_otu-1;++i)
    {
      for(j=0;j<tree->times->n_time_slices;j++)
	{
	  if(!(tree->times->nd_t[i] > tree->times->time_slice_lims[j])) break;
	}
      tree->times->curr_slice[i] = j;

      /* PhyML_Printf("\n. Node %3d [%12f] is in slice %3d.",i,tree->times->nd_t[i],j); */
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl TIMES_Wrap_Lk_Coalescent(t_edge *b, t_tree *tree, supert_tree *stree)
{
  return TIMES_Lk_Coalescent(tree);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl TIMES_Lk_Coalescent(t_tree *tree)
{  
  t_node *n;
  int n_lineages;
  phydbl lnP,Ne,dt,exp_g;
  phydbl T;
  
  Ne = tree->times->scaled_pop_size;
  exp_g = tree->times->neff_growth;
  
  if(Ne > tree->times->scaled_pop_size_max || Ne < tree->times->scaled_pop_size_min)
    {
      tree->times->c_lnL = UNLIKELY;
      return(tree->times->c_lnL);
    }
  
  Get_Node_Ranks_From_Times(tree);

  n = tree->n_root;
  while(n->rk_next) n = n->rk_next;
  T = tree->times->nd_t[n->num];
  
  lnP = 0.0;
  n_lineages = 1;
  n = tree->n_root;
  while(n->rk_next)
    {
      dt = fabs(tree->times->nd_t[n->num] - tree->times->nd_t[n->rk_next->num]);
      
      assert(!(tree->times->nd_t[n->num] > tree->times->nd_t[n->rk_next->num]));
      
      if(n->tax == YES) n_lineages--;
      else n_lineages++;


      if(tree->times->coalescent_model_id == EXPCOALESCENT)
        {
          dt =
            1./exp_g *
            (exp(exp_g * fabs(tree->times->nd_t[n->num]-T)) -
             exp(exp_g * fabs(tree->times->nd_t[n->rk_next->num]-T)));
        }
      else if(tree->times->coalescent_model_id == POWLAW)
        {
          dt =
            1./(1.+exp_g) *
            (pow(fabs(tree->times->nd_t[n->num]-T)+1,exp_g+1.) -
             pow(fabs(tree->times->nd_t[n->rk_next->num]-T)+1,exp_g+1.));
        }
      
      lnP += -n_lineages * (n_lineages-1.) / (2.*Ne) * dt;
      
      if(tree->times->coalescent_model_id == EXPCOALESCENT && n->tax == NO)
        lnP += exp_g * fabs(tree->times->nd_t[n->num]-T);
      else if(tree->times->coalescent_model_id == POWLAW && n->tax == NO)
        lnP += exp_g * log(fabs(tree->times->nd_t[n->num]-T) + 1.);
      
      if(n->tax == NO) lnP -= log(Ne);
      
      // Multifurcations not allowed under standard Kingman coalescent
      if(fabs(tree->times->nd_t[n->num] - tree->times->nd_t[n->rk_next->num]) < SMALL &&
         n->tax == NO && n->rk_next->tax == NO)
        {
          tree->times->c_lnL = UNLIKELY;
          return(UNLIKELY);
        }
            
      if((isnan(lnP) || isinf(lnP)) == TRUE)
        {
          tree->times->c_lnL = UNLIKELY;
          return(UNLIKELY);
        }
        
      n = n->rk_next;
    }
  
  tree->times->c_lnL = lnP;

  return(lnP);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl TIMES_Lk_Coalescent_Range(t_dsk *young, t_dsk *old, t_tree *tree)
{
  int n_lineages;
  phydbl disk_lnP,lnP,Ne,dt,exp_g,T;
  t_dsk *disk;

  Ne = tree->times->scaled_pop_size;  
  exp_g = tree->times->neff_growth;
  T = tree->young_disk->time;
  
  if(Ne > tree->times->scaled_pop_size_max || Ne < tree->times->scaled_pop_size_min) return(UNLIKELY);

  
  
  lnP = 0.0;
  disk_lnP = 0.0;
  disk = young->prev;
  do
    {
      n_lineages = disk->next->n_ldsk_a;
      n_lineages -= (disk->next->ldsk ? (disk->next->ldsk->n_next -1) : 0);

      dt = 0.0;
      
      if(tree->times->coalescent_model_id == STRICTCOALESCENT)
        dt = fabs(disk->time - disk->next->time);
      else if(tree->times->coalescent_model_id == EXPCOALESCENT)
        dt =
          1./exp_g *
          (exp(exp_g * fabs(disk->time-T)) -
           exp(exp_g * fabs(disk->next->time-T)));


      disk_lnP = -n_lineages * (n_lineages-1.) / (2.*Ne) * dt;

      if(tree->times->coalescent_model_id == EXPCOALESCENT)
        disk_lnP += exp_g * fabs(disk->time-T);
      
      lnP += disk_lnP;

      // Multifurcations not allowed under standard Kingman coalescent
      if(fabs(disk->time - disk->next->time) < SMALL &&
         disk->ldsk && disk->ldsk->n_next > 1 &&
         disk->next->ldsk && disk->next->ldsk->n_next > 1)
        {
          tree->times->c_lnL = UNLIKELY;
          return(UNLIKELY);
        }

      assert((isnan(lnP) || isinf(lnP)) == FALSE);
      disk = disk->prev;
    }
  while(disk != NULL && disk != old->prev);
  
  lnP -= (tree->n_otu-1) * log(Ne);

  return(lnP);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl TIMES_Prior_Coalescent(t_tree *tree)
{
  phydbl lbda;
  lbda = 1./tree->times->neff_prior_mean;
  return(log(lbda)-lbda*tree->times->scaled_pop_size);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl TIMES_Prior(t_tree *tree)
{
  tree->times->c_lnP = UNLIKELY;
    
  switch(tree->times->model_id)
    {
    case COALESCENT :
      {
        tree->times->c_lnP = TIMES_Prior_Coalescent(tree);
        break;
      }
    case SLFV_GAUSSIAN : case SLFV_UNIFORM :
      {
        tree->times->c_lnP = TIMES_Prior_SLFV(tree);
        break;
      }
    case BIRTHDEATH : 
      {
        tree->times->c_lnP = UNLIKELY;
        break;
      }     
    default : break;
    }

  return(tree->times->c_lnP);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void TIMES_Simulate_Coalescent(t_tree *tree)
{
  t_node *n;
  phydbl t,dt,Ne,*times,u,coal_rate,prob_coal;
  int n_lineages,*permut,available_idx,i;
  t_ll *lineages;
  t_node *v1, *v2, *a;
  
  
  Ne = tree->times->scaled_pop_size;
  times = tree->times->nd_t;
  lineages = NULL;
  
  Get_Node_Ranks_From_Tip_Times(tree);


  // Start from the most recent tip (or set of tips)
  n = tree->a_nodes[tree->n_otu-1];  
  while(n->rk_next) n = n->rk_next;
  
  n_lineages = 0;
  available_idx = tree->n_otu;
  dt = 0.0;
  t = tree->times->nd_t[n->num];
  do
    {
      if(n->tax == YES)
        {
          Push_Bottom_Linked_List(n,&lineages,NO);
          n_lineages++;
        }
      
      coal_rate = n_lineages * (n_lineages - 1.) / (2.*Ne);

      // Proba next coalescent event is in [T[n->num],T[n->rk_prev->num] 
      prob_coal = 1. - exp(-coal_rate * fabs(times[n->num] - (n->rk_prev ? times[n->rk_prev->num] : -INFINITY)));
            
      u = Uni();
      
      if(!(u > prob_coal))
        {
          dt = Rexp_Trunc(coal_rate,
                          0.0,
                          times[n->num] - (n->rk_prev ? times[n->rk_prev->num] : -INFINITY));
                          
          t = t - dt;
                    
          permut = Permutate(n_lineages);

          a = tree->a_nodes[available_idx];
          v1 = Linked_List_Elem(permut[0],lineages);
          v2 = Linked_List_Elem(permut[1],lineages);
          available_idx++;

          a->v[0]  = NULL;
          a->v[1]  = v1;
          a->v[2]  = v2;
          v1->v[0] = a;
          v2->v[0] = a;

          times[a->num] = t;
          
          Remove_From_Linked_List(NULL,v1,&lineages);
          Remove_From_Linked_List(NULL,v2,&lineages);
          Push_Bottom_Linked_List(a,&lineages,NO);

          a->rk_prev = n->rk_prev;
          n->rk_prev = a;
          a->rk_next = n;
          
          n_lineages--;
        }
      
      n = n->rk_prev;
      t = tree->times->nd_t[n->num];
    }
  while(!(n_lineages == 1 && n->rk_prev == NULL));

  for(i=0;i<3;++i) 
    if(tree->n_root->v[1]->v[i] == tree->n_root) 
      { 
        tree->n_root->v[1]->v[i] = tree->n_root->v[2]; 
        break; 
      }
  
  for(i=0;i<3;++i) 
    if(tree->n_root->v[2]->v[i] == tree->n_root) 
      { 
        tree->n_root->v[2]->v[i] = tree->n_root->v[1]; 
        break; 
      }
  
  Connect_Edges_To_Nodes_Serial(tree);
  
  tree->e_root = NULL;
  for(i=0;i<2*tree->n_otu-3;++i)
    {
      if((tree->a_edges[i]->left == tree->n_root->v[1] && tree->a_edges[i]->rght == tree->n_root->v[2]) ||
         (tree->a_edges[i]->left == tree->n_root->v[2] && tree->a_edges[i]->rght == tree->n_root->v[1]))
        {        
          tree->e_root = tree->a_edges[i];
          break;
        }
    }
  assert(!(tree->e_root == NULL));
  
  tree->n_root->b[1] = tree->a_edges[2*tree->n_otu-3];
  tree->n_root->b[2] = tree->a_edges[2*tree->n_otu-2];
  
  tree->n_root->b[1]->left = tree->n_root;
  tree->n_root->b[1]->rght = tree->n_root->v[1];
  
  tree->n_root->b[2]->left = tree->n_root;
  tree->n_root->b[2]->rght = tree->n_root->v[2];
  
  Update_Ancestors(tree->n_root,tree->n_root->v[2],tree->n_root->b[2],tree);
  Update_Ancestors(tree->n_root,tree->n_root->v[1],tree->n_root->b[1],tree);
  
  MIXT_Propagate_Tree_Update(tree);  
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl TIMES_Lk_Uniform_Core(t_tree *tree)
{  
  phydbl logn;

  logn = TIMES_Log_Number_Of_Ranked_Labelled_Histories(tree->n_root,YES,tree);

  tree->times->c_lnL = 0.0;
  TIMES_Lk_Uniform_Post(tree->n_root,tree->n_root->v[2],tree);
  TIMES_Lk_Uniform_Post(tree->n_root,tree->n_root->v[1],tree);

  /* printf("\n. ^ %f %f %f", */
  /* 	 (phydbl)(tree->rates->n_tips_below[tree->n_root->num]-2.), */
  /* 	 log(tree->times->time_slice_lims[tree->times->curr_slice[tree->n_root->num]] - */
  /* 	     tree->times->nd_t[tree->n_root->num]), */
  /* 	 (phydbl)(tree->rates->n_tips_below[tree->n_root->num]-2.) * */
  /* 	 log(tree->times->time_slice_lims[tree->times->curr_slice[tree->n_root->num]] - */
  /* 	     tree->times->nd_t[tree->n_root->num])); */

  tree->times->c_lnL +=
    Factln(tree->rates->n_tips_below[tree->n_root->num]-2.) -
    (phydbl)(tree->rates->n_tips_below[tree->n_root->num]-2.) *
    log(tree->times->time_slice_lims[tree->times->curr_slice[tree->n_root->num]] -
  	tree->times->nd_t[tree->n_root->num]);
  
  tree->times->c_lnL -= logn;
  
  return(tree->times->c_lnL);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void TIMES_Get_Number_Of_Time_Slices(t_tree *tree)
{
  int i;

  tree->times->n_time_slices=0;
  TIMES_Get_Number_Of_Time_Slices_Post(tree->n_root,tree->n_root->v[2],tree);
  TIMES_Get_Number_Of_Time_Slices_Post(tree->n_root,tree->n_root->v[1],tree);
  Qksort(tree->times->time_slice_lims,NULL,0,tree->times->n_time_slices-1);

  if(tree->times->n_time_slices > 1)
    {
      PhyML_Printf("\n");
      PhyML_Printf("\n. Sequences were collected at %d different time points.",tree->times->n_time_slices);
      for(i=0;i<tree->times->n_time_slices;i++) printf("\n+ [%3d] time point @ %12f ",i+1,tree->times->time_slice_lims[i]);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void TIMES_Get_Number_Of_Time_Slices_Post(t_node *a, t_node *d, t_tree *tree)
{
  int i;

  if(d->tax == YES)
    {
      for(i=0;i<tree->times->n_time_slices;i++) 
	if(Are_Equal(tree->times->t_floor[d->num],tree->times->time_slice_lims[i],1.E-6)) break;

      if(i == tree->times->n_time_slices) 
	{
	  tree->times->time_slice_lims[i] = tree->times->t_floor[d->num];
	  tree->times->n_time_slices++;
	}
      return;
    }
  else
    {
      for(i=0;i<3;i++)
	if(d->v[i] != a && d->b[i] != tree->e_root)
	  TIMES_Get_Number_Of_Time_Slices_Post(d,d->v[i],tree);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void TIMES_Get_N_Slice_Spans(t_tree *tree)
{
  int i,j;

  For(i,2*tree->n_otu-2)
    {
      if(tree->a_nodes[i]->tax == NO)
	{
	  for(j=0;j<tree->times->n_time_slices;j++)
	    {
	      if(Are_Equal(tree->times->t_floor[i],tree->times->time_slice_lims[j],1.E-6))
		{
		  tree->times->n_time_slice_spans[i] = j+1;
		  /* PhyML_Printf("\n. Node %3d spans %3d slices [%12f].", */
		  /* 	       i+1, */
		  /* 	       tree->rates->n_slice_spans[i], */
		  /* 	       tree->times->t_floor[i]); */
		  break;
		}
	    }
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void TIMES_Lk_Uniform_Post(t_node *a, t_node *d, t_tree *tree)
{
  if(d->tax == YES) return;
  else
    {
      int i;

      for(i=0;i<3;i++)
	{
	  if(d->v[i] != a && d->b[i] != tree->e_root)
	    {
	      TIMES_Lk_Uniform_Post(d,d->v[i],tree);
	    }
	}
      
      if(tree->times->curr_slice[a->num] != tree->times->curr_slice[d->num])
	{
	  tree->times->c_lnL += 
	    Factln(tree->rates->n_tips_below[d->num]-1.) - 
	    (phydbl)(tree->rates->n_tips_below[d->num]-1.) *
	    log(tree->times->time_slice_lims[tree->times->curr_slice[d->num]] -
		tree->times->nd_t[d->num]);
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* Set the root position so that most of the taxa in the outgroup 
   correspond to the most ancient time point.
*/
void TIMES_Set_Root_Given_Tip_Dates(t_tree *tree)
{
  int i,j;
  t_node *left,*rght;
  int n_left_in, n_left_out;
  int n_rght_in, n_rght_out;
  t_edge *b,*best;
  phydbl eps,score,max_score;
  
  Free_Bip(tree);
  Alloc_Bip(tree);
  Get_Bip(tree->a_nodes[0],tree->a_nodes[0]->v[0],tree);
  
  left = rght = NULL;
  b = best = NULL;
  n_left_in = n_left_out = -1;
  n_rght_in = n_rght_out = -1;
  eps = 1.E-6;
  score = max_score = -1.;

  For(i,2*tree->n_otu-3)
    {
      left = tree->a_edges[i]->left;
      rght = tree->a_edges[i]->rght;
      b    = tree->a_edges[i];

      n_left_in = 0;
      For(j,left->bip_size[b->l_r]) 
	if(FABS(tree->times->nd_t[left->bip_node[b->l_r][j]->num] - tree->times->time_slice_lims[0]) < eps)
	  n_left_in++;
      
      n_left_out = left->bip_size[b->l_r]-n_left_in;
      
      n_rght_in = 0;
      For(j,rght->bip_size[b->r_l]) 
	if(FABS(tree->times->nd_t[rght->bip_node[b->r_l][j]->num] - tree->times->time_slice_lims[0]) < eps)
	  n_rght_in++;

      n_rght_out = rght->bip_size[b->r_l]-n_rght_in;


      /* score = POW((phydbl)(n_left_in)/(phydbl)(n_left_in+n_left_out)- */
      /* 		  (phydbl)(n_rght_in)/(phydbl)(n_rght_in+n_rght_out),2); */
      /* score = (phydbl)(n_left_in * n_rght_out + eps)/(n_left_out * n_rght_in + eps); */
      /* score = (phydbl)(n_left_in * n_rght_out + eps); */
      score = FABS((phydbl)((n_left_in+1.) * (n_rght_out+1.)) - (phydbl)((n_left_out+1.) * (n_rght_in+1.)));
      
      if(score > max_score)
	{
	  max_score = score;
	  best = b;
	}
    }
  
  Add_Root(best,tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* Update the ranking of node heights. Use bubble sort algorithm */
/* t_rank[i] is the node number that has rank i */

void TIMES_Update_Node_Ordering(t_tree *tree)
{
  int buff;
  int i;
  phydbl *t;
  int swap = NO;

  for(i=0;i<2*tree->n_otu-1;++i) tree->times->t_rank[i] = i;

  t = tree->times->nd_t;

  do
    {
      swap = NO;
      for(i=0;i<2*tree->n_otu-2;++i)
	{
	  if(t[tree->times->t_rank[i]] > t[tree->times->t_rank[i+1]]) // Sort in ascending order
	    {
	      swap = YES;
	      buff                     = tree->times->t_rank[i];
	      tree->times->t_rank[i]   = tree->times->t_rank[i+1];
	      tree->times->t_rank[i+1] = buff;
	    }	    
	}
    }
  while(swap == YES);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void TIMES_Set_Calibration(t_tree *tree)
{
  t_cal *cal;
  int i;

  For(i,2*tree->n_otu-1)
    {
      tree->times->t_has_prior[i] = NO;
      tree->times->t_prior_min[i] = BIG;
      tree->times->t_prior_max[i] = BIG; 
   }

  cal = tree->times->a_cal[0];
  while(cal)
    {
      /* if(cal->is_active == YES) */
      /*   { */
          /* tree->times->t_has_prior[cal->node_num] = YES; */
          /* tree->times->t_prior_min[cal->node_num] = cal->lower; */
          /* tree->times->t_prior_max[cal->node_num] = cal->upper;           */
        /* } */
      cal = cal->next;
    }

  TIMES_Set_All_Node_Priors(tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void TIMES_Record_Prior_Times(t_tree *tree)
{
  int i;
  for(i=0;i<2*tree->n_otu-1;++i) 
    {
      tree->times->t_prior_min_ori[i] = tree->times->t_prior_min[i];
      tree->times->t_prior_max_ori[i] = tree->times->t_prior_max[i];
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void TIMES_Reset_Prior_Times(t_tree *tree)
{
  int i;
  For(i,2*tree->n_otu-1) 
    {
      tree->times->t_prior_min[i] = tree->times->t_prior_min_ori[i];
      tree->times->t_prior_max[i] = tree->times->t_prior_max_ori[i];
     }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// Returns the conditional density of internal node heights 
// given the tree height under the Yule model. Uses the order
// statistics 'simplification' as described in Yang and Rannala, 2005. 
phydbl TIMES_Lk_Yule_Order_Root_Cond(t_tree *tree)
{
  int j;
  phydbl *t,*tf;
  t_node *n;
  phydbl loglk;
  phydbl loglbda;
  phydbl lbda;
  phydbl *tp_min,*tp_max;
  phydbl lower_bound,upper_bound;
  phydbl root_height;

  tp_min = tree->times->t_prior_min;
  tp_max = tree->times->t_prior_max;
  tf = tree->times->t_floor;
  t  = tree->times->nd_t;
  n = NULL;
  loglbda = log(tree->times->birth_rate);
  lbda = tree->times->birth_rate;
  lower_bound = -1.;
  upper_bound = -1.;
  root_height = FABS(tree->times->nd_t[tree->n_root->num]);

  /*! Adapted from  Equation (6) in T. Stadler's Systematic Biology, 2012 paper with
      sampling fraction set to 1 and death rate set to 0. Dropped the 1/(n-1) scaling 
      factor. */

  /* loglk = 0.0; */
  /* For(j,2*tree->n_otu-2) */
  /*   { */
  /*     n = tree->a_nodes[j]; */
  /*     lower_bound = MAX(FABS(tf[j]),FABS(tp_max[j])); */
  /*     upper_bound = MIN(FABS(t[tree->n_root->num]),FABS(tp_min[j])); */

  /*     if(n->tax == NO) */
  /*       { */
  /*         loglk  += (loglbda - lbda * FABS(t[j])); */
  /*         /\* loglk -= log(exp(-lbda*lower_bound) - exp(-lbda*upper_bound)); // incorporate calibration boundaries here. *\/ */
  /*       } */
  /*   } */

  
  /*! Adapted from  Equation (7) in T. Stadler's Systematic Biology, 2012 paper with
      sampling fraction set to 1 and death rate set to 0. */

  // Check that each node is within its calibration-derived interval
  For(j,2*tree->n_otu-1) if(t[j] < tp_min[j] || t[j] > tp_max[j]) return(-INFINITY);

  loglk = 0.0;
  For(j,2*tree->n_otu-2)
    {
      n = tree->a_nodes[j];
      lower_bound = MAX(FABS(tf[j]),FABS(tp_max[j]));
      upper_bound = MIN(FABS(tp_min[j]),root_height);
      
      if(n->tax == NO)
        {
          loglk  += (loglbda - lbda * FABS(t[j]));
          loglk -= log(exp(-lbda*lower_bound) - exp(-lbda*upper_bound)); // incorporate calibration boundaries here.    
        }

        if(isinf(loglk) || isnan(loglk))
          {
            PhyML_Fprintf(stderr,"\n. Lower bound: %f",lower_bound);
            PhyML_Fprintf(stderr,"\n. Upper bound: %f",upper_bound);
            PhyML_Fprintf(stderr,"\n. tf: %f tp_max: %f tp_min: %f root: %f",tf[j],tp_max[j],tp_min[j],root_height);
            Exit("\n");
          }

    }

  /* lower_bound = MAX(FABS(tf[tree->n_root->num]),FABS(tp_max[tree->n_root->num])); */
  /* upper_bound = FABS(tp_min[tree->n_root->num]); */
  /* loglk += log(2) + loglbda - 2.*lbda * FABS(t[tree->n_root->num]); */
  /* loglk -= log(exp(-2.*lbda*lower_bound) - exp(-2.*lbda*upper_bound)); */


  return(loglk);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// Log of prob density of internal node ages conditional on tree height, under
// the birth death process with incomplete sampling.
phydbl TIMES_Lk_Birth_Death(int verbose, t_tree *tree)
{
  int i,n;
  phydbl lnL;
  phydbl t,b,d,bmd,logbmd,expmbmd,logb;
  phydbl bmin,bmax;
  phydbl dmin,dmax;
  phydbl lognut1,logp_1t,troot,pt,nut1;
  
  lnL     = 0.0;
  b       = tree->times->birth_rate;
  d       = tree->times->death_rate;
  bmin    = tree->times->birth_rate_min;
  bmax    = tree->times->birth_rate_max;
  dmin    = tree->times->death_rate_min;
  dmax    = tree->times->death_rate_min;
  bmd     = b-d;
  logbmd  = -1.;
  expmbmd = -1.;
  logb    = -1.;
  t       = 0.0;
  n       = tree->n_otu;
  troot   = fabs(tree->times->nd_t[tree->n_root->num]);
  logb    = log(b);
  
  if(b < d)
    {
      tree->times->c_lnL = UNLIKELY;
      if(verbose) PhyML_Printf("\n. b < d");
      return UNLIKELY;
    }
  
  // Verify that calibration constraints are satisfied
  for(i=0;i<2*tree->n_otu-1;++i)
    if(tree->a_nodes[i]->tax == NO)
      {
        if(tree->times->nd_t[i] < tree->times->t_prior_min[i] ||
           tree->times->nd_t[i] > tree->times->t_prior_max[i]) 
          {
            tree->times->c_lnL = UNLIKELY;
            if(verbose)
              {
                PhyML_Printf("\n. Time outside calibration range : %G [%G,%G]",tree->times->nd_t[i],tree->times->t_prior_min[i],tree->times->t_prior_max[i]);
                PhyML_Printf("\n. Clade incriminated: ");
                if(tree->a_nodes[i] == tree->n_root)
                  {
                    List_Taxa_In_Clade(tree->n_root,tree->n_root->v[1],tree);
                    List_Taxa_In_Clade(tree->n_root,tree->n_root->v[2],tree);
                  }
                else
                  {
                    assert(tree->a_nodes[i]->anc);
                    List_Taxa_In_Clade(tree->a_nodes[i]->anc,tree->a_nodes[i],tree);
                  }
                PhyML_Printf("\n");
              }
            return UNLIKELY;
          }
      }
  
  if(b > bmin && d > dmin && Are_Equal(bmd,0.0,bmin/10.) == NO)
    {
      logbmd = log(bmd);
      expmbmd = exp(d-b);
      
      pt = bmd/(b-d*pow(expmbmd,troot));
      nut1 = 1. - pt*pow(expmbmd,troot);
      lognut1 = log(nut1);
      /* printf("\n. lognut1: %G pt: %G b: %G d: %G exp: %G",nut1,pt,b,d,expmbmd); fflush(NULL); */
      
      for(i=0;i<2*tree->n_otu-1;++i)
        if(tree->a_nodes[i]->tax == NO && tree->a_nodes[i] != tree->n_root)
          {
            t = fabs(tree->times->nd_t[i]);            

            /* // Equation 3.19 in Tanja Stadler's PhD thesis */
            /* lnL += 2*logbmd - bmd*t - 2.*log(b-d*pow(expmbmd,t));             */

            // Equation 6 in Yang and Rannala, 1997 with rho=1
            logp_1t = 2.*logbmd - 2.*log(b-d*exp((d-b)*t)) + (d-b)*t;
            lnL += logb + logp_1t - lognut1; 

            
            if(!(lnL > UNLIKELY))
              {
                PhyML_Printf("\n. logb: %G pt: %G nut1: %G lognut1: %G t: %G logp_1t: %G\n",logb,pt,nut1,lognut1,t,logp_1t);
                tree->times->c_lnL = UNLIKELY;
                return UNLIKELY;
              }
          }

      lnL += LnGamma(n-1);

      
      /* t = FABS(tree->times->nd_t[tree->n_root->num]); */
      /* lnL += -bmd*t - log(b-d*pow(expmbmd,t)); */
      /* lnL += log(bmd) + (tree->n_otu-1)*log(b) + LnGamma(tree->n_otu+1); */

      // Divide joint density p(x(1),...,x(n-1)) by p(x(1)) in order to get
      // the conditional density p(x(2),...,x(n-1)|x(1)). The log of the marginal p(x(1))
      // (Theorem 3.4.11 in Tanja's PhD thesis) is given below and subtracted to lnL.
      /* k = 1; */
      /* s = FABS(tree->times->nd_t[tree->n_root->num]); */
      /* phydbl p_x1 = (k+1)*Choose(n,k+1)*pow(b,n-k)*pow(bmd,k+2)*exp(-bmd*(k+1)*s)*pow(1.-exp(-bmd*s),n-k-1)/pow(b-d*exp(-bmd*s),n+1); */
      /* lnL -= log(p_x1); */
    }
  else if(b > bmin && d < dmin) // Yule process 
    {
      lognut1 = log(1.-exp(-b*troot));
      
      for(i=0;i<2*tree->n_otu-1;++i)
        if(tree->a_nodes[i]->tax == NO && tree->a_nodes[i] != tree->n_root)
          {
            t = fabs(tree->times->nd_t[i]);            
            /* // Equation 3.19 in Tanja Stadler's PhD thesis (Yule case) */
            /* lnL -= b*t; */

            // Equation 6 in Yang and Rannala, 1997 with rho=1
            logp_1t = - b*t;
            lnL += logb + logp_1t - lognut1; 

            if(!(lnL > UNLIKELY))              
              {
                PhyML_Printf("\n. lognut1: %G t: %G logp_1t: %G",lognut1,t,logp_1t);
                tree->times->c_lnL = UNLIKELY;
                return UNLIKELY;
              }
          }

      lnL += LnGamma(n-1);

      /* t = fabs(tree->times->nd_t[tree->n_root->num]); */
      /* lnL -= b*t;       */
      /* lnL += (tree->n_otu-1)*log(b) + LnGamma(tree->n_otu+1); */


      // Divide joint density p(x(1),...,x(n-1)) by p(x(1)) in order to get
      // the conditional density p(x(2),...,x(n-1)|x(1)). The log of the marginal p(x(1))
      // (Theorem 3.4.11, remark 3.4.12 in Tanja's PhD thesis) is given below and subtracted to lnL.
      /* k = 1; */
      /* s = FABS(tree->times->nd_t[tree->n_root->num]); */
      /* phydbl p_x1 = (k+1)*Choose(n,k+1)*b*pow(exp(b*s)-1.,n-k-1)/exp(b*s*n); */
      /* lnL -= log(p_x1); */
    }
  else if(b < bmin && d > dmin) 
    {
      PhyML_Printf("\n. b: %G bmin: %G d: %G dmin: %G",b,bmin,d,dmin);
      tree->times->c_lnL = UNLIKELY;
      return UNLIKELY;
    }
  else if(Are_Equal(bmd,0.0,bmin/10.) == YES) // Critical birth-death process
    {
      logb = log(b);

      for(i=0;i<2*tree->n_otu-1;++i)
        if(tree->a_nodes[i]->tax == NO && tree->a_nodes[i] != tree->n_root)
          {
            t = fabs(tree->times->nd_t[i]);            

            /* // Equation 3.19 in Tanja Stadler's PhD thesis (Critical case) */
            /* lnL += logb - 2.*log(1.+b*t); */

            
            // Equation 7 in Yang and Rannala, 1997 with rho=1
            lnL += log((1.+d)/pow(1.+d*t,2)); 

            if(!(lnL > UNLIKELY))              
              {
                PhyML_Printf("\n. logb: %G t: %G",logb,t);
                tree->times->c_lnL = UNLIKELY;
                return UNLIKELY;
              }
          }

      lnL += LnGamma(n-1);
      
      /* t = fabs(tree->times->nd_t[tree->n_root->num]); */
      /* lnL -= log(b*t); */
      /* lnL += LnGamma(tree->n_otu+1); */
      
      // Divide joint density p(x(1),...,x(n-1)) by p(x(1)) in order to get
      // the conditional density p(x(2),...,x(n-1)|x(1)). The log of the marginal p(x(1))
      // (Theorem 3.4.11, remark 3.4.12 in Tanja's PhD thesis) is given below and subtracted to lnL.
      /* k = 1; */
      /* s = fabs(tree->times->nd_t[tree->n_root->num]); */
      /* phydbl p_x1 = (k+1)*Choose(n,k+1)*pow(b,n-k)*pow(s,n-k-1)/pow(1.+b*s,n+1); */
      /* lnL -= log(p_x1); */
    }
  else if(b < bmin && d < dmin) // Birth and death rates are below their limits
    {
      PhyML_Fprintf(stderr,"\n. b: %G bmin: %G d: %G dmin: %G",b,bmin,d,dmin);
      tree->times->c_lnL = UNLIKELY;
      return -INFINITY;
    }
  else if(b > bmax && d > dmax)
    {
      PhyML_Fprintf(stderr,"\n. b: %G bmax: %G d: %G dmax: %G",b,bmax,d,dmax);
      tree->times->c_lnL = UNLIKELY;
      return -INFINITY;
    }
  else
    {
      assert(FALSE);
    }

  if(isnan(lnL) || isinf(fabs(lnL)) || !(lnL > UNLIKELY))
    {
      PhyML_Fprintf(stderr,"\n. lnL times: %f",lnL);
      tree->times->c_lnL = UNLIKELY;
      return UNLIKELY;
    }

  tree->times->c_lnL = lnL;

  return(lnL);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

// Generate a subtree including all taxa in tax_list. The age of the root of that
// subtree is t_mrca. All nodes in the subtree are thus younger than that.
void TIMES_Connect_List_Of_Taxa(t_node **tax_list, int list_size, phydbl t_mrca, phydbl *times, int *nd_num, t_tree *mixt_tree)
{
  phydbl t_upper_bound, t_lower_bound,up,lo;
  int i,j,n_anc,*permut,rand_idx;
  t_node *n,**anc,*new_mrca;

  // Calibration of a tip node
  if(list_size == 1)
    {
      times[tax_list[0]->num] = t_mrca;
      return;
    }
  else
    {
      t_lower_bound = t_mrca;
      t_upper_bound = +INFINITY;
      n             = NULL;
      anc           = NULL;
      new_mrca      = NULL;
      permut        = NULL;
      
      // Find the upper bound for all the new node ages that 
      // will be created in this function
      for(i=0;i<list_size;i++)
        {
          n = tax_list[i];
          while(n->v[0] != NULL) n = n->v[0];
          if(times[n->num] < t_upper_bound) t_upper_bound = times[n->num];
        }
      
      if(t_upper_bound < t_lower_bound)
        {
          PhyML_Printf("\n. upper: %f lower: %f t_mrca: %f\n",t_upper_bound,t_lower_bound,t_mrca);
          assert(FALSE);
        }
      
      // Get the list of current mrcas to all taxa in tax_list. There should be
      // at least one of these
      n_anc = 0;
      for(i=0;i<list_size;i++)
        {
          n = tax_list[i];
          while(n->v[0] != NULL) n = n->v[0];
          for(j=0;j<n_anc;j++) if(anc[j] == n) break;
          if(j == n_anc)
            {
              if(n_anc == 0) anc = (t_node **)mCalloc(1,sizeof(t_node *));
              else           anc = (t_node **)mRealloc(anc,n_anc+1,sizeof(t_node *));
              anc[n_anc] = n;
              n_anc++;
            }
        }
      
      /* printf("\n. n_anc: %d",n_anc); */
      /* for(i=0;i<n_anc;i++) PhyML_Printf("\n. anc: %d",anc[i]->num); */
      
      if(n_anc == 1) // All the nodes in tax_list are already connected. Bail out.
        {
          Free(anc);
          return;
        }
      
      // Connect randomly and set ages
      permut = Permutate(n_anc);
      /* permut = (int *)mCalloc(n_anc,sizeof(int)); */
      /* for(i=0;i<n_anc;i++) permut[i] = i; */
      i = 0;
      do
        {
          new_mrca               = mixt_tree->a_nodes[*nd_num];
          anc[permut[i]]->v[0]   = new_mrca;
          anc[permut[i+1]]->v[0] = new_mrca;
          new_mrca->v[1]         = anc[permut[i]];
          new_mrca->v[2]         = anc[permut[i+1]];
          new_mrca->v[0]         = NULL;
          up                     = MIN(times[new_mrca->v[1]->num],times[new_mrca->v[2]->num]);
          lo                     = up - (up - t_mrca)/5.;
          times[new_mrca->num]   = Uni()*(up - lo) + lo;
          
          /* times[new_mrca->num]   = t_upper_bound - ((phydbl)(i+1.)/n_anc)*(t_upper_bound - t_lower_bound); */
                    
          /* printf("\n. new_mrca->num: %d time: %f [%f %f] t_mrca: %f %d connect to %d %d %d [%f %f]", */
          /*        new_mrca->num, */
          /*        times[new_mrca->num], */
          /*        t_lower_bound, */
          /*        t_upper_bound, */
          /*        t_mrca, */
          /*        new_mrca->num, */
          /*        new_mrca->v[0] ? new_mrca->v[0]->num : -1, */
          /*        new_mrca->v[1] ? new_mrca->v[1]->num : -1, */
          /*        new_mrca->v[2] ? new_mrca->v[2]->num : -1, */
          /*        times[new_mrca->v[1]->num], */
          /*        times[new_mrca->v[2]->num] */
          /*        ); */
          /* fflush(NULL); */

          anc[permut[i+1]] = new_mrca;

          rand_idx = Rand_Int(i+1,n_anc-1);
          n                     = anc[permut[i+1]];
          anc[permut[i+1]]      = anc[permut[rand_idx]];
          anc[permut[rand_idx]] = n;

          
          i++;
          (*nd_num) += 1;
          if(n_anc == i+1) { times[new_mrca->num] = t_mrca; break; }
        }
      while(1);
      
      Free(permut);
      Free(anc);
    }
}
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// Generate a  random rooted tree with node ages fullfiling the
// time constraints defined in the list of calibration cal_list 
void TIMES_Randomize_Tree_With_Time_Constraints(t_cal *cal_list, t_tree *mixt_tree)
{
  t_node **tips,**nd_list;
  phydbl *times,*cal_times,time_oldest_cal;
  int i,j,nd_num,*cal_ordering,n_cal,swap,list_size,tmp,orig_is_mixt_tree,repeat,n_max_repeats,tip_num,*no_cal_tip_num,n_no_cal_tip_num,*permut,*tip_is_in_cal_clad;
  t_cal *cal;
  t_clad *clade;

  assert(mixt_tree->rates);

  if(mixt_tree->mod->s_opt->opt_topo == NO)
    {
      PhyML_Fprintf(stderr,"\n. Fixing the tree topology is only allowed when calibrating tip nodes only.");
      PhyML_Fprintf(stderr,"\n. You are most likely calibrating here the MRCA of at least one clade with more than two tips.");
      PhyML_Fprintf(stderr,"\n. It is difficult to set the age of that clade within the limit of the calibration constraints,");
      PhyML_Fprintf(stderr,"\n. and fix the tree topology at the same time. Please contact me (guindon@lirmm.fr) for a more");
      PhyML_Fprintf(stderr,"\n. detailed diagnostic.");
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
    }
  
  tips               = (t_node **)mCalloc(mixt_tree->n_otu,sizeof(t_node *));
  nd_list            = (t_node **)mCalloc(2*mixt_tree->n_otu-1,sizeof(t_node *));
  tip_is_in_cal_clad = (int *)mCalloc(mixt_tree->n_otu,sizeof(t_node *));
                                      
  times                   = mixt_tree->times->nd_t;
  orig_is_mixt_tree       = mixt_tree->is_mixt_tree;
  mixt_tree->is_mixt_tree = NO;
  n_max_repeats           = 1000;
  cal                     = NULL;
  clade                   = NULL;
  

  // Is tip in calibrated clade or not?
  for(i=0;i<mixt_tree->n_otu;++i)
    {
      cal = cal_list;
      clade = NULL;
      while(cal != NULL)
        {
          while(cal && cal->clade_list == NULL) cal = cal->next;
          if(cal == NULL) break;
          clade = cal->clade_list[cal->current_clade_idx];
          for(j=0;j<clade->n_tax;++j) if(clade->tip_list[j] == mixt_tree->a_nodes[i]) break;
          if(j != clade->n_tax) break;
          cal = cal->next;
        }
      
      if(cal != NULL) tip_is_in_cal_clad[i] = YES;
      else tip_is_in_cal_clad[i] = NO;
    }

  // List node indices that are not in any calibration set
  no_cal_tip_num = NULL;
  n_no_cal_tip_num = 0;
  for(i=0;i<mixt_tree->n_otu;++i)
    {
      cal = cal_list;
      while(cal != NULL)
        {
          while(cal && cal->clade_list == NULL) cal = cal->next;
          if(cal == NULL) break;
          clade = cal->clade_list[cal->current_clade_idx];
          for(j=0;j<clade->n_tax;++j) if(clade->tip_list[j] == mixt_tree->a_nodes[i]) break;
          if(j != clade->n_tax) break;
          cal = cal->next;
        }
      
      if(cal == NULL)
        {
          if(n_no_cal_tip_num == 0) no_cal_tip_num = (int *)mCalloc(1,sizeof(int));
          else no_cal_tip_num = (int *)mRealloc(no_cal_tip_num,n_no_cal_tip_num+1,sizeof(int));
          no_cal_tip_num[n_no_cal_tip_num] = i;
          n_no_cal_tip_num++;
        }
    }

  
  for(repeat=0;repeat<n_max_repeats;++repeat)
    {
      nd_num = mixt_tree->n_otu;
  
      /* PhyML_Printf("\n\n. Repeat %d",repeat); */

      for(i=0;i<mixt_tree->n_otu;++i) tips[i] = mixt_tree->a_nodes[i];
      for(i=0;i<mixt_tree->n_otu;++i) mixt_tree->a_nodes[i]->v[0] = NULL; 


      // Set a time for each calibration 
      cal_times = (phydbl *)mCalloc(1,sizeof(phydbl));
      time_oldest_cal = 0.0;
      n_cal = 0;
      cal = cal_list;
      while(cal != NULL)
        {
          if(cal->is_primary == YES)
            {
              if(n_cal > 0) cal_times = (phydbl *)mRealloc(cal_times,n_cal+1,sizeof(phydbl));
              cal_times[n_cal] = Uni()*(cal->upper - cal->lower) + cal->lower;
              if(cal_times[n_cal] < time_oldest_cal) time_oldest_cal = cal_times[n_cal];
              n_cal++;
            }
          cal = cal->next;
        }

      cal_ordering = (int *)mCalloc(n_cal,sizeof(int));
      for(i=0;i<n_cal;++i) cal_ordering[i] = i;
      
      // Sort calibration times from youngest to oldest
      do
        {
          swap = NO;
          for(i=0;i<n_cal-1;i++) 
            {
              if(cal_times[cal_ordering[i]] < cal_times[cal_ordering[i+1]])
                {
                  tmp               = cal_ordering[i+1];
                  cal_ordering[i+1] = cal_ordering[i];
                  cal_ordering[i]   = tmp;
                  swap = YES;
                }
            }
        }
      while(swap == YES);
      
      for(i=0;i<n_cal-1;i++) assert(!(cal_times[cal_ordering[i]] < cal_times[cal_ordering[i+1]]));
              
      
      // Connect taxa that appear in every primary calibrations
      for(i=0;i<n_cal;++i)
        {
          cal = cal_list;
          j = 0;
          while(j != cal_ordering[i]) 
            { 
              if(cal->is_primary == YES) j++; 
              cal = cal->next; 
              assert(cal); 
            }
          
          list_size = 0;

          /* printf("\n. n_cal: %d n_no_cal_tip_num: %d [%d %d] n_otu: %d", */
          /*        n_cal, */
          /*        n_no_cal_tip_num, */
          /*        n_no_cal_tip_num-1, */
          /*        (int)(2*n_no_cal_tip_num/(n_cal)), */
          /*        mixt_tree->n_otu); fflush(NULL); */
          
          // Add some taxa that are not in any calibration set
          if(n_no_cal_tip_num > 0)
            {
              permut = Permutate(n_no_cal_tip_num);
              j = 0;
              do
                {
                  tip_num = no_cal_tip_num[permut[j]];
                  if(tips[tip_num]->v[0] == NULL)
                    {
                      nd_list[list_size] = tips[tip_num];
                      /* PhyML_Printf("\n# %s",tips[tip_num]->name); */
                      list_size++;
                      assert(list_size <= mixt_tree->n_otu);
                    }
                }
              while(j++ < MIN(n_no_cal_tip_num-1,(int)(2*n_no_cal_tip_num/(n_cal))));
              Free(permut);
            }
          

          // Add all the taxa that are in the calibration set of cal
          // This should be done here so that the last node in nd_list
          // belongs to the taxa in the calibration

          if(cal->clade_list != NULL)
            {
              clade = cal->clade_list[cal->current_clade_idx];
              for(j=0;j<clade->n_tax;j++)
                {
                  nd_list[list_size] = tips[clade->tip_list[j]->num];
                  list_size++;
                  assert(list_size <= mixt_tree->n_otu);
                }

              TIMES_Connect_List_Of_Taxa(nd_list, 
                                         list_size,
                                         cal_times[cal_ordering[i]], 
                                         times, 
                                         &nd_num,
                                         mixt_tree);
            }
        }
      
      Free(cal_times);
      Free(cal_ordering);
      
      
      // Connect all remaining taxa
      for(i=0;i<mixt_tree->n_otu;i++) nd_list[i] = NULL;
      list_size = 0;
      for(i=0;i<mixt_tree->n_otu;++i)
        {
          if(tip_is_in_cal_clad[tips[i]->num] == NO) // Tip is not in a calibrated clade
            {
              nd_list[list_size] = tips[i];
              list_size++;
              assert(list_size <= mixt_tree->n_otu);
            }
        }
      
      cal = cal_list;
      do
        {
          while(cal && cal->clade_list == NULL) cal = cal->next;
          if(cal == NULL) break;
          clade = cal->clade_list[cal->current_clade_idx];
          nd_list[list_size] = clade->tip_list[0];
          list_size++;
          assert(list_size <= 2*mixt_tree->n_otu-1);
          cal = cal->next;
        }
      while(cal);

      /* for(i=0;i<list_size;i++) printf("\n# To connect: %d lower: %f",nd_list[i]->num,time_oldest_cal); */

      
      TIMES_Connect_List_Of_Taxa(nd_list,
                                 list_size,
                                 10.*time_oldest_cal-1.E-3, // 1.E-3 required in case time_oldest_cal = 0.
                                 times,
                                 &nd_num,
                                 mixt_tree);

      // Adding root node 
      mixt_tree->n_root = mixt_tree->a_nodes[2*mixt_tree->n_otu-2];
      mixt_tree->n_root->v[1]->v[0] = mixt_tree->n_root->v[2];
      mixt_tree->n_root->v[2]->v[0] = mixt_tree->n_root->v[1];    
      mixt_tree->n_root->anc = NULL;
      
      Connect_Edges_To_Nodes_Serial(mixt_tree);
      

      // Adding root edge
      for(i=0;i<2*mixt_tree->n_otu-3;++i)
        {
          if(((mixt_tree->a_edges[i]->left == mixt_tree->n_root->v[1]) || (mixt_tree->a_edges[i]->rght == mixt_tree->n_root->v[1])) &&
             ((mixt_tree->a_edges[i]->left == mixt_tree->n_root->v[2]) || (mixt_tree->a_edges[i]->rght == mixt_tree->n_root->v[2])))
            {
              mixt_tree->e_root = mixt_tree->a_edges[i]; 
              break;
            }
        }
      
      Update_Ancestors(mixt_tree->n_root,mixt_tree->n_root->v[2],mixt_tree->n_root->b[2],mixt_tree);
      Update_Ancestors(mixt_tree->n_root,mixt_tree->n_root->v[1],mixt_tree->n_root->b[1],mixt_tree);
      DATE_Assign_Primary_Calibration(mixt_tree);
      DATE_Update_T_Prior_MinMax(mixt_tree);
      

      /* for(int i=0;i<mixt_tree->times->n_cal;i++) */
      /*   { */
      /*     int idx; */
      /*     cal   = mixt_tree->times->a_cal[i]; */

      /*     if(cal->clade_list != NULL) */
      /*       { */
      /*         clade = cal->clade_list[cal->current_clade_idx]; */
      /*         idx = Find_Clade(clade->tax_list,clade->n_tax,mixt_tree); */
      /*         PhyML_Printf("\n. Calibration %3d | Node number to which calibration applies is [%d]",i,idx); */
      /*         PhyML_Printf("\n. Calibration %3d | Lower bound set to: %15f time units.",i,mixt_tree->times->a_cal[i]->lower); */
      /*         PhyML_Printf("\n. Calibration %3d | Upper bound set to: %15f time units.",i,mixt_tree->times->a_cal[i]->upper); */
      /*         PhyML_Printf("\n. Calibration %3d | t_prior_min: %G t_prior_max: %G",i,mixt_tree->times->t_prior_min[idx],mixt_tree->times->t_prior_max[idx]); */
      /*         PhyML_Printf("\n. Calibration %3d | Time set to %G",i,mixt_tree->times->nd_t[idx]); */
      /*       } */
      /*   } */
      
      if(!DATE_Check_Calibration_Constraints(mixt_tree))
        {
          /* PhyML_Printf("\n. Could not generate tree (DATE_Check_Calibration_Constraints)\n\n"); */
        }
      else if(!DATE_Check_Time_Constraints(mixt_tree))
        {
          /* PhyML_Printf("\n. Could not generate tree (DATE_Check_Time_Constraints)\n\n"); */
        }
      else break; // Tree successfully generated
    }

    
  if(repeat == n_max_repeats)
    {
      PhyML_Fprintf(stderr,"\n\n");
      PhyML_Fprintf(stderr,"\n. A random tree satisfying the calibration constraints provided");
      PhyML_Fprintf(stderr,"\n. could not be generated. It probably means that there are some");
      PhyML_Fprintf(stderr,"\n. inconsistencies in the calibration data. For instance, the calibration");
      PhyML_Fprintf(stderr,"\n. time interval for the MRCA of a clade with taxa {X,Y} (noted as [a,b])");
      PhyML_Fprintf(stderr,"\n. cannot be strictly older than the interval corresponding to taxa ");
      PhyML_Fprintf(stderr,"\n. {X,Z,Y} (noted as [c,d]), i.e., b cannot be smaller (older) than c. ");
      PhyML_Fprintf(stderr,"\n. Also, please remember that the present time corresponds to a time");
      PhyML_Fprintf(stderr,"\n. value equal to zero and past events have negative time values.");

      PhyML_Fprintf(stderr,"\n\n");
      if(!DATE_Check_Calibration_Constraints(mixt_tree))
        {
          PhyML_Fprintf(stderr,"\n. Could not generate tree (DATE_Check_Calibration_Constraints)\n\n");
        }
      else if(!DATE_Check_Time_Constraints(mixt_tree))
        {
          PhyML_Fprintf(stderr,"\n. Could not generate tree (DATE_Check_Time_Constraints)\n\n");
        }

      Exit("\n");
    }

  /* assert(i != 2*mixt_tree->n_otu-3); */

  mixt_tree->is_mixt_tree = orig_is_mixt_tree;
    
  Free(tips);
  Free(nd_list);
  Free(no_cal_tip_num);
  Free(tip_is_in_cal_clad);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int TIMES_Check_Node_Height_Ordering(t_tree *tree)
{
  if(!TIMES_Check_Node_Height_Ordering_Post(tree->n_root,tree->n_root->v[1],tree)) return NO;
  if(!TIMES_Check_Node_Height_Ordering_Post(tree->n_root,tree->n_root->v[2],tree)) return NO;
  return YES;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int TIMES_Check_Node_Height_Ordering_Post(t_node *a, t_node *d, t_tree *tree)
{

  if(d->anc != a)
    {
      PhyML_Printf("\n. d=%d d->anc=%d a=%d root=%d",d->num,d->anc->num,a->num,tree->n_root->num);
      return NO;
    }
  if(tree->times->nd_t[d->num] < tree->times->nd_t[a->num])
    {
      PhyML_Printf("\n. a->t = %f [num:%d] d->t %f [num:%d]",
                   tree->times->nd_t[a->num],
                   a->num,
                   tree->times->nd_t[d->num],
                   d->num);
      return NO;
    }
  if(d->tax == YES) return YES;
  else
    {
      int i;
      for(i=0;i<3;i++)
        {
          if(d->v[i] != a && d->b[i] != tree->e_root)
            if(!TIMES_Check_Node_Height_Ordering_Post(d,d->v[i],tree))
              return NO;
        }
    }
  return YES;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* 
   2n-2 edge lengths (in a rooted tree), n-1 internal node
   ages. Return \sum_i (l_i - delta t_i)^2 where delta t_i is the
   calendar time elapsed along the edge of length l_i  
*/
phydbl TIMES_Least_Square_Criterion(t_tree *tree)
{
  phydbl score;
  score = 0.0;
  assert(tree->n_root);
  assert(tree->rates);
  TIMES_Pre_Least_Square_Criterion(tree->n_root,tree->n_root->v[1],tree->n_root->b[1],&score,tree);
  TIMES_Pre_Least_Square_Criterion(tree->n_root,tree->n_root->v[2],tree->n_root->b[2],&score,tree);
  return(score);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void TIMES_Pre_Least_Square_Criterion(t_node *a, t_node *d, t_edge *b, phydbl *score, t_tree *tree)
{
  int i;

  (*score) += pow(b->l->v - fabs(tree->times->nd_t[a->num] + tree->times->nd_t[d->num])*tree->rates->clock_r,2);

  if(d->tax) return;
  
  for(i=0;i<3;++i)
    if(d->v[i] != a && d->b[i] != tree->e_root)
      TIMES_Pre_Least_Square_Criterion(d,d->v[i],d->b[i],score,tree);  
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void TIMES_Randomize_Node_Ages(t_tree *tree)
{
  assert(tree->n_root);
  assert(tree->rates);
  TIMES_Post_Randomize_Node_Ages(tree->n_root,tree->n_root->v[1],tree);
  TIMES_Post_Randomize_Node_Ages(tree->n_root,tree->n_root->v[2],tree);
  tree->times->nd_t[tree->n_root->num] =
    MIN(tree->times->nd_t[tree->n_root->v[1]->num],
        tree->times->nd_t[tree->n_root->v[2]->num])
    - Rexp(1.);
  /* PhyML_Printf("\n. RAND TIMES node %3d %15f",tree->n_root->num,tree->times->nd_t[tree->n_root->num]); */
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void TIMES_Post_Randomize_Node_Ages(t_node *a, t_node *d, t_tree *tree)
{
  if(d->tax == YES) return;
  else
    {
      int i,dir1,dir2;
      
      for(i=0;i<3;++i)
        if(d->v[i] != a && d->b[i] != tree->e_root)
          TIMES_Post_Randomize_Node_Ages(d,d->v[i],tree);

      dir1 = dir2 = -1;
      for(i=0;i<3;++i)
        if(d->v[i] != a && d->b[i] != tree->e_root)
          {
            if(dir1 < 0) dir1 = i;
            else dir2 = i;
          }
      
      tree->times->nd_t[d->num] =
        MIN(tree->times->nd_t[d->v[dir1]->num],
            tree->times->nd_t[d->v[dir2]->num])
        - Rexp(1.);
      /* PhyML_Printf("\n. RAND TIMES node %3d %15f",d->num,tree->times->nd_t[d->num]); */
    }
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int TIMES_Calibrations_Apply_To_Tips_Only(t_tree *tree)
{
  t_cal *cal;
  t_clad *clade;

  cal = tree->times->a_cal[0];
  assert(cal);
  clade = NULL;
  
  do
    {
      while(cal && cal->clade_list == NULL) cal = cal->next;
      if(cal == NULL) break;
      clade = cal->clade_list[cal->current_clade_idx];
      if(clade && clade->n_tax > 1) break;
      cal = cal->next;
    }
  while(cal);

  if(cal == NULL) return YES;

  return NO;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void TIMES_Randomize_Tip_Times_Given_Calibrations(t_tree *tree)
{
  
  t_cal *cal;
  t_clad *clade;

  cal = tree->times->a_cal[0];
  assert(cal);
  
  do
    {
      clade = cal->clade_list[cal->current_clade_idx];

      if(clade->n_tax == 1)
        {
          assert(clade->target_nd->tax == YES);
          tree->times->nd_t[clade->target_nd->num] = Uni()*(cal->upper - cal->lower) + cal->lower;
        }
      
      cal = cal->next;
    }
  while(cal);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void TIMES_Time_To_Bl(t_tree *tree)
{
  TIMES_Time_To_Bl_Pre(tree->n_root,tree->n_root->v[1],tree->n_root->b[1],tree);
  TIMES_Time_To_Bl_Pre(tree->n_root,tree->n_root->v[2],tree->n_root->b[2],tree);
  tree->n_root->b[1]->l->v = tree->times->nd_t[tree->n_root->v[1]->num] - tree->times->nd_t[tree->n_root->num];
  tree->n_root->b[2]->l->v = tree->times->nd_t[tree->n_root->v[2]->num] - tree->times->nd_t[tree->n_root->num];
  tree->e_root->l->v = tree->n_root->b[1]->l->v + tree->n_root->b[2]->l->v;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void TIMES_Time_To_Bl_Pre(t_node *a, t_node *d, t_edge *b, t_tree *tree)
{
  int i;

  b->l->v = tree->times->nd_t[d->num] - tree->times->nd_t[a->num];
  
  if(d->tax) return;
  else
    {
      for(i=0;i<3;i++)
        if((d->v[i] != a) && (d->b[i] != tree->e_root))
          TIMES_Time_To_Bl_Pre(d,d->v[i],d->b[i],tree);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl TIMES_Tree_Length(t_tree *tree)
{
  phydbl sum;

  assert(tree->rates);
  assert(tree->e_root);
  
  sum = 0.0;
  for(int i=0;i<2*tree->n_otu-1;++i)
    if(tree->a_edges[i] != tree->e_root)
      sum += fabs(tree->times->nd_t[tree->a_edges[i]->left->num]-
                  tree->times->nd_t[tree->a_edges[i]->rght->num]);

  return(sum);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void TIMES_Bl_To_Times(t_tree *tree)
{
  t_node *v1,*v2;
  int dir1,dir2;
  phydbl t1,t2;
  
  TIMES_Bl_To_Times_Post(tree->n_root,tree->n_root->v[1],tree->n_root->b[1],tree);
  TIMES_Bl_To_Times_Post(tree->n_root,tree->n_root->v[2],tree->n_root->b[2],tree);

  dir1 = 1;
  dir2 = 2;
  
  v1 = tree->n_root->v[dir1];
  v2 = tree->n_root->v[dir2];
  
  t1 = tree->times->nd_t[v1->num] - MIXT_Get_Mean_Edge_Len(tree->n_root->b[1],tree) / (tree->rates->clock_r * tree->rates->br_r[v1->num]);
  t2 = tree->times->nd_t[v2->num] - MIXT_Get_Mean_Edge_Len(tree->n_root->b[2],tree) / (tree->rates->clock_r * tree->rates->br_r[v2->num]);
  
  if(Are_Equal(t1,t2,1.E-6) == NO)
    {
      PhyML_Fprintf(stderr,"\n. It looks as if the edge lengths suplied do not define an ultrametric tree.");
      PhyML_Fprintf(stderr,"\n. Please amend these lengths so as it becomes straightforward to transform your tree");
      PhyML_Fprintf(stderr,"\n. into a time-tree. Please contact me (guindon@lirmm.fr) for more information.");
      PhyML_Fprintf(stderr,"\n. l1: %f l2: %f",MIXT_Get_Mean_Edge_Len(tree->n_root->b[1],tree),MIXT_Get_Mean_Edge_Len(tree->n_root->b[2],tree));
      PhyML_Fprintf(stderr,"\n. t1: %f t2: %f",tree->times->nd_t[v1->num],tree->times->nd_t[v2->num]);
      PhyML_Fprintf(stderr,"\n. rr1: %f rr2: %f",tree->rates->br_r[v1->num],tree->rates->br_r[v2->num]);
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
    }
      
  tree->times->nd_t[tree->n_root->num] = t1;
  
  
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void TIMES_Bl_To_Times_Post(t_node *a, t_node *d, t_edge *b, t_tree *tree)
{
  if(d->tax == YES) return;
  else
    {
      t_node *v1,*v2;
      int dir1,dir2;
      phydbl t1,t2;
      
      dir1 = dir2 = -1;
      for(int i=0;i<3;i++)
        if((d->v[i] != a) && (d->b[i] != tree->e_root))
          {
            TIMES_Bl_To_Times_Post(d,d->v[i],d->b[i],tree);
            if(dir1 < 0) dir1 = i;
            else dir2 = i;
          }
      
      v1 = d->v[dir1];
      v2 = d->v[dir2];
      
      t1 = tree->times->nd_t[v1->num] - MIXT_Get_Mean_Edge_Len(d->b[dir1],tree) / (tree->rates->clock_r * tree->rates->br_r[v1->num]);
      t2 = tree->times->nd_t[v2->num] - MIXT_Get_Mean_Edge_Len(d->b[dir2],tree) / (tree->rates->clock_r * tree->rates->br_r[v2->num]);

      if(Are_Equal(t1,t2,1.E-6) == NO)
        {
          PhyML_Fprintf(stderr,"\n. It looks at if the edge lengths suplied do not define an ultrametric tree.");
          PhyML_Fprintf(stderr,"\n. Please amend these lengths so as it becomes straightforward to transform your tree");
          PhyML_Fprintf(stderr,"\n. into a time-tree.");
          PhyML_Fprintf(stderr,"\n. l1: %f l2: %f",MIXT_Get_Mean_Edge_Len(d->b[dir1],tree),MIXT_Get_Mean_Edge_Len(d->b[dir2],tree));
          PhyML_Fprintf(stderr,"\n. t1: %f t2: %f",tree->times->nd_t[v1->num],tree->times->nd_t[v2->num]);
          PhyML_Fprintf(stderr,"\n. rr1: %f rr2: %f",tree->rates->br_r[v1->num],tree->rates->br_r[v2->num]);
          PhyML_Fprintf(stderr,"\n. est: %f %f diff: %G",t1,t2,t1-t2);
          Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
        }
      
      tree->times->nd_t[d->num] = t1;
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void TIMES_Copy_Time_Struct(t_time *from, t_time *to, int n_otu)
{
  int i;

  to->birth_rate = from->birth_rate;
  to->birth_rate_min = from->birth_rate_min;
  to->birth_rate_max = from->birth_rate_max;
  to->birth_rate_pivot = from->birth_rate_pivot;

  to->death_rate = from->death_rate;
  to->death_rate_min = from->death_rate_min;
  to->death_rate_max = from->death_rate_max;
  to->death_rate_pivot = from->death_rate_pivot;
  
  to->c_lnL = from->c_lnL;
  to->c_lnL_jps = from->c_lnL_jps;
  
  to->nd_t_recorded = from->nd_t_recorded;
  to->update_time_norm_const = from->update_time_norm_const;

  to->scaled_pop_size = from->scaled_pop_size;
  to->scaled_pop_size_min = from->scaled_pop_size_min;
  to->scaled_pop_size_max = from->scaled_pop_size_max;

  to->neff_growth = from->neff_growth;
  to->neff_growth_min = from->neff_growth_min;
  to->neff_growth_max = from->neff_growth_max;

  to->model_id = from->model_id;
  to->coalescent_model_id = from->coalescent_model_id;

  to->update_time_norm_const = from->update_time_norm_const;

  to->augmented_coalescent = from->augmented_coalescent;

  to->neff_prior_mean = from->neff_prior_mean;
  
  for(i=0;i<2*n_otu-1;++i) to->nd_t[i] = from->nd_t[i];
  for(i=0;i<2*n_otu-1;++i) to->buff_t[i] = from->buff_t[i];
  for(i=0;i<2*n_otu-1;++i) to->true_t[i] = from->true_t[i];
  for(i=0;i<2*n_otu-1;++i) to->t_mean[i] = from->t_mean[i];
  for(i=0;i<2*n_otu-1;++i) to->t_prior[i] = from->t_prior[i];
  for(i=0;i<2*n_otu-1;++i) to->t_prior_min[i] = from->t_prior_min[i];
  for(i=0;i<2*n_otu-1;++i) to->t_prior_max[i] = from->t_prior_max[i];
  for(i=0;i<2*n_otu-1;++i) to->t_floor[i] = from->t_floor[i];
  for(i=0;i<2*n_otu-1;++i) to->t_rank[i] = from->t_rank[i];
  for(i=0;i<2*n_otu-1;++i) to->t_has_prior[i] = from->t_has_prior[i];
  for(i=0;i<2*n_otu-1;++i) to->n_jps[i] = from->n_jps[i];
  for(i=0;i<2*n_otu-2;++i) to->t_jps[i] = from->t_jps[i];
  for(i=0;i<2*n_otu-1;++i) to->mean_t[i] = from->mean_t[i];
  for(i=0;i<2*n_otu-1;++i) to->time_slice_lims[i] = from->time_slice_lims[i];
  for(i=0;i<2*n_otu-1;++i) to->n_time_slice_spans[i] = from->n_time_slice_spans[i];
  for(i=0;i<2*n_otu-1;++i) to->curr_slice[i] = from->curr_slice[i];
  for(i=0;i<2*n_otu-1;++i) to->has_survived[i] = from->has_survived[i];
  for(i=0;i<2*n_otu-1;++i) to->calib_prob[i] = from->calib_prob[i];
  for(i=0;i<2*n_otu-1;++i) to->t_prior_min_ori[i] = from->t_prior_max_ori[i];
  for(i=0;i<n_otu*n_otu;++i) to->times_partial_proba[i] = from->times_partial_proba[i];
  for(i=0;i<n_otu*n_otu;++i) to->numb_calib_chosen[i] = from->numb_calib_chosen[i];
  for(i=0;i<2*n_otu-1;++i) to->mean_t[i] = from->mean_t[i];
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl TIMES_Lk(t_tree *tree)
{
  switch(tree->times->model_id)
    {
    case COALESCENT :
      {
        tree->times->c_lnL = TIMES_Lk_Coalescent(tree);
        break;
      }
    case SLFV_GAUSSIAN : case SLFV_UNIFORM :
      {
        tree->times->c_lnL = TIMES_Lk_SLFV(tree);        
        break;
      }
    case BIRTHDEATH : 
      {
        DATE_Assign_Primary_Calibration(tree);
        DATE_Update_T_Prior_MinMax(tree);
        tree->times->c_lnL =  TIMES_Lk_Birth_Death(NO,tree);
        break;
      }     
    default : break;
    }

  return(tree->times->c_lnL);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl TIMES_Lk_Range(t_dsk *young, t_dsk *old, t_tree *tree)
{
  switch(tree->times->model_id)
    {
    case COALESCENT :
      {
        return(TIMES_Lk_Coalescent_Range(young,old,tree));
        break;
      }
    case SLFV_GAUSSIAN : case SLFV_UNIFORM :
      {
        return(TIMES_Lk_SLFV_Range(young,old,tree));
        break;
      }
    default : break;
    }
  return(-1.);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl TIMES_Prior_SLFV(t_tree *tree)
{
  return(UNLIKELY);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl TIMES_Lk_SLFV(t_tree *tree)
{
  phydbl lnL;
  t_dsk *disk;
  phydbl dt;
  int n_evt;

  n_evt = 0;
  dt    = 0.0;
  lnL   = 0.0;
  disk  = tree->young_disk->prev;
  do
    {
      if(disk->age_fixed == NO)
        {
          dt += fabs(disk->next->time - disk->time);
          n_evt++;
        }
      
      disk = disk->prev;
    }
  while(disk);

  lnL += n_evt * log(tree->mmod->lbda) - tree->mmod->lbda * dt;

  tree->times->c_lnL = lnL;
  
  return(tree->times->c_lnL);  
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl TIMES_Lk_SLFV_Range(t_dsk *young, t_dsk *old, t_tree *tree)
{
  t_dsk *disk;
  int n_evt;
  phydbl lnP,dt;

  assert(young);

  lnP   = 0.0;    
  dt    = 0.0;
  n_evt = 0;
  disk = young->prev;
  do
    {
      if(disk->age_fixed == NO)
        {
          dt += fabs(disk->next->time - disk->time);
          n_evt++;
        }

      assert(disk);
      if(disk->time > disk->next->time) return UNLIKELY;
      if(disk == old) break;
      disk = disk->prev;
    }
  while(disk);

  lnP = n_evt * log(tree->mmod->lbda) - tree->mmod->lbda * dt;
 
  return(lnP);
}
  
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

  
