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

  /* Solve A.x = b, where x are the t_node time estimated
     under the least square criterion.

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

  TIMES_Least_Square_Node_Times_Pre(root,root->v[2],A,b,n,tree);
  TIMES_Least_Square_Node_Times_Pre(root,root->v[1],A,b,n,tree);
  
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

  for(i=0;i<n-1;i++) { tree->rates->nd_t[tree->a_nodes[i]->num] = -x[i]; }
  tree->rates->nd_t[root->num] = -x[n-1];
  tree->n_root->l[2] = tree->rates->nd_t[root->v[2]->num] - tree->rates->nd_t[root->num];
  tree->n_root->l[1] = tree->rates->nd_t[root->v[1]->num] - tree->rates->nd_t[root->num];
  ////////////////////////////////////////
  return;
  ////////////////////////////////////////

  /* Rescale the t_node times such that the time at the root
     is -100. This constraint implies that the clock rate
     is fixed to the actual tree length divided by the tree
     length measured in term of differences of t_node times */

  phydbl scale_f,time_tree_length,tree_length;

  scale_f = -100./tree->rates->nd_t[root->num];
  For(i,2*tree->n_otu-1) tree->rates->nd_t[i] *= scale_f;
  For(i,2*tree->n_otu-1) if(tree->rates->nd_t[i] > .0) tree->rates->nd_t[i] = .0;

  time_tree_length = 0.0;
  For(i,2*tree->n_otu-3)
    if(tree->a_edges[i] != tree->e_root)
      time_tree_length +=
	FABS(tree->rates->nd_t[tree->a_edges[i]->left->num] -
	     tree->rates->nd_t[tree->a_edges[i]->rght->num]);
  time_tree_length += FABS(tree->rates->nd_t[root->num] - tree->rates->nd_t[root->v[2]->num]);
  time_tree_length += FABS(tree->rates->nd_t[root->num] - tree->rates->nd_t[root->v[1]->num]);
  
  tree_length = 0.0;
  For(i,2*tree->n_otu-3) tree_length += tree->a_edges[i]->l->v;

  tree->rates->clock_r = tree_length / time_tree_length;

  Free(A);
  Free(b);
  Free(x);

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
      
 
      for(i=0;i<3;i++)
	if((d->v[i] != a) && (d->b[i] != tree->e_root))
	  TIMES_Least_Square_Node_Times_Pre(d,d->v[i],A,b,n,tree);
      
      A[d->num * n + d->num] = 1.;
      b[d->num] = .0;
      for(i=0;i<3;i++)
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

  if(tree->rates->nd_t[tree->n_root->num] > MIN(tree->rates->nd_t[tree->n_root->v[2]->num],
						tree->rates->nd_t[tree->n_root->v[1]->num]))
    {
      tree->rates->nd_t[tree->n_root->num] = MIN(tree->rates->nd_t[tree->n_root->v[2]->num],
						 tree->rates->nd_t[tree->n_root->v[1]->num]);
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
	      if(tree->rates->nd_t[d->v[i]->num] < min_height)
		{
		  min_height = tree->rates->nd_t[d->v[i]->num];
		}
	    }
	}

      if(tree->rates->nd_t[d->num] > min_height) tree->rates->nd_t[d->num] = min_height;

      if(tree->rates->nd_t[d->num] < -100.) tree->rates->nd_t[d->num] = -100.;

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
  For(i,2*tree->n_otu-2) tree->rates->nd_t[tree->a_nodes[i]->num] *= FABS(tree->mod->s_opt->tree_size_mult);
  tree->rates->nd_t[tree->n_root->num] *= FABS(tree->mod->s_opt->tree_size_mult);
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
	       tree->rates->nd_t[a->num],
	       tree->rates->nd_t[d->num],
	       tree->rates->nd_t[a->num]-tree->rates->nd_t[d->num],
	       (b)?(b->l->v):(-1.0),
	       tree->rates->t_prior_min[d->num],
	       tree->rates->t_prior_max[d->num]);
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

  tree->rates->t_prior_max[tree->n_root->num] = 
    MIN(tree->rates->t_prior_max[tree->n_root->num],
	MIN(tree->rates->t_prior_max[tree->n_root->v[2]->num],
	    tree->rates->t_prior_max[tree->n_root->v[1]->num]));


  /* Set all t_prior_min values */
  if(!tree->rates->t_has_prior[tree->n_root->num])
    {
      min_prior = 1.E+10;
      For(i,2*tree->n_otu-2)
	{
	  if(tree->rates->t_has_prior[i])
	    {
	      if(tree->rates->t_prior_min[i] < min_prior)
		min_prior = tree->rates->t_prior_min[i];
	    }
	}
      tree->rates->t_prior_min[tree->n_root->num] = 2.0 * min_prior;
      /* tree->rates->t_prior_min[tree->n_root->num] = 10. * min_prior; */
    }
  
  if(tree->rates->t_prior_min[tree->n_root->num] > 0.0)
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
      
      if(tree->rates->t_has_prior[d->num] == YES)
	{
	  t_sup = MIN(tree->rates->t_prior_max[d->num],
		      MIN(tree->rates->t_prior_max[v1->num],
			  tree->rates->t_prior_max[v2->num]));

	  tree->rates->t_prior_max[d->num] = t_sup;

	  if(tree->rates->t_prior_max[d->num] < tree->rates->t_prior_min[d->num])
	    {
	      PhyML_Fprintf(stderr,"\n. prior_min=%f prior_max=%f",tree->rates->t_prior_min[d->num],tree->rates->t_prior_max[d->num]);
	      PhyML_Fprintf(stderr,"\n. Inconsistency in the prior settings detected at node %d",d->num);
	      PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d (function %s)\n\n",__FILE__,__LINE__,__FUNCTION__);
	      Warn_And_Exit("\n");
	    }
	}
      else
	{
	  tree->rates->t_prior_max[d->num] = 
	    MIN(tree->rates->t_prior_max[v1->num],
		tree->rates->t_prior_max[v2->num]);
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
      
      if(tree->rates->t_has_prior[d->num] == YES)
	{
	  tree->rates->t_prior_min[d->num] = MAX(tree->rates->t_prior_min[d->num],tree->rates->t_prior_min[a->num]);
	  
	  if(tree->rates->t_prior_max[d->num] < tree->rates->t_prior_min[d->num])
	    {
	      PhyML_Fprintf(stderr,"\n. prior_min=%f prior_max=%f",tree->rates->t_prior_min[d->num],tree->rates->t_prior_max[d->num]);
	      PhyML_Fprintf(stderr,"\n. Inconsistency in the prior settings detected at t_node %d",d->num);
	      PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d (function %s)\n\n",__FILE__,__LINE__,__FUNCTION__);
	      Warn_And_Exit("\n");
	    }
	}
      else
	{
	  tree->rates->t_prior_min[d->num] = tree->rates->t_prior_min[a->num];
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
  tree->rates->t_floor[tree->n_root->num] = MIN(tree->rates->t_floor[tree->n_root->v[2]->num],
						tree->rates->t_floor[tree->n_root->v[1]->num]);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void TIMES_Set_Floor_Post(t_node *a, t_node *d, t_tree *tree)
{
  if(d->tax)
    {
      tree->rates->t_floor[d->num] = tree->rates->nd_t[d->num];
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
      tree->rates->t_floor[d->num] = MIN(tree->rates->t_floor[v1->num],
					 tree->rates->t_floor[v2->num]);

      if(tree->rates->t_floor[v1->num] < tree->rates->t_floor[v2->num])
	{
	  d->rank_max = v1->rank_max;
	}
      else if(tree->rates->t_floor[v2->num] < tree->rates->t_floor[v1->num])
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

  min = tree->rates->nd_t[tree->n_root->num];

  dens = 0.0;
  For(i,2*tree->n_otu-1)
    {
      if((tree->a_nodes[i]->tax == NO) && (tree->a_nodes[i] != tree->n_root))
	{
	  max = tree->rates->t_floor[i];

	  dens += log(Dorder_Unif(tree->rates->nd_t[i],
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

  lbda = tree->rates->birth_rate;
  t    = tree->rates->nd_t;
  ts   = tree->rates->time_slice_lims;
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

  t = tree->rates->nd_t;
  ts = tree->rates->time_slice_lims;
  tr = tree->rates->t_rank;
  lbda = tree->rates->birth_rate;

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

  /* printf("\n. sumdt = %f th=%f",sumdt,tree->rates->nd_t[tree->n_root->num]); */
  /* printf("\n0 loglk = %f",loglk); */

  for(i=0;i<tree->rates->n_time_slices-1;i++)
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

  tp_min = tree->rates->t_prior_min;
  tp_max = tree->rates->t_prior_max;
  tf = tree->rates->t_floor;
  t  = tree->rates->nd_t;
  n = NULL;
  loglbda = log(tree->rates->birth_rate);
  lbda = tree->rates->birth_rate;
  lower_bound = -1.;
  upper_bound = -1.;
  /* root_height = FABS(tree->rates->nd_t[tree->n_root->num]); */

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

phydbl TIMES_Lk_Times(int verbose, t_tree *tree)
{
  DATE_Assign_Primary_Calibration(tree);
  DATE_Update_T_Prior_MinMax(tree);
  tree->rates->c_lnL_times =  TIMES_Lk_Birth_Death(verbose,tree);
  /* if(tree->rates->c_lnL_times > UNLIKELY) tree->rates->c_lnL_times  = -1.; */
  return(tree->rates->c_lnL_times);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void TIMES_Lk_Times_Trav(t_node *a, t_node *d, phydbl lim_inf, phydbl lim_sup, phydbl *logdens, t_tree *tree)
{
  int i;
  
  if(!d->tax)
    {
      /* if(tree->rates->nd_t[d->num] > lim_sup) */
      /* 	{ */
      /* 	  lim_inf = lim_sup; */
      /* 	  lim_sup = 0.0; */
      /* 	  For(i,2*tree->n_otu-2) */
      /* 	    if((tree->rates->t_floor[i] < lim_sup) && (tree->rates->t_floor[i] > tree->rates->nd_t[d->num])) */
      /* 	      lim_sup = tree->rates->t_floor[i]; */
      /* 	} */
      
      /* if(tree->rates->nd_t[d->num] < lim_inf || tree->rates->nd_t[d->num] > lim_sup) */
      /* 	{ */
      /* 	  PhyML_Printf("\n. nd_t = %f lim_inf = %f lim_sup = %f",tree->rates->nd_t[d->num],lim_inf,lim_sup); */
      /* 	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__); */
      /* 	  Exit("\n"); */
      /* 	} */
  
      lim_inf = tree->rates->nd_t[tree->n_root->num];
      lim_sup = tree->rates->t_floor[d->num];
      
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
      if(tree->rates->curr_slice[v1->num] == tree->rates->curr_slice[root->num])
  	n1 = tree->rates->n_tips_below[v1->num];
      else
  	n1 = 1;
      
      if(tree->rates->curr_slice[v2->num] == tree->rates->curr_slice[root->num])
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
	  if(tree->rates->curr_slice[v1->num] == tree->rates->curr_slice[d->num])
	    n1 = tree->rates->n_tips_below[v1->num];
	  else
	    n1 = 1;

	  if(tree->rates->curr_slice[v2->num] == tree->rates->curr_slice[d->num])
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

  For(i,2*tree->n_otu-1)
    {
      for(j=0;j<tree->rates->n_time_slices;j++)
	{
	  if(!(tree->rates->nd_t[i] > tree->rates->time_slice_lims[j])) break;
	}
      tree->rates->curr_slice[i] = j;

      /* PhyML_Printf("\n. Node %3d [%12f] is in slice %3d.",i,tree->rates->nd_t[i],j); */
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl TIMES_Lk_Uniform_Core(t_tree *tree)
{  
  phydbl logn;

  logn = TIMES_Log_Number_Of_Ranked_Labelled_Histories(tree->n_root,YES,tree);

  tree->rates->c_lnL_times = 0.0;
  TIMES_Lk_Uniform_Post(tree->n_root,tree->n_root->v[2],tree);
  TIMES_Lk_Uniform_Post(tree->n_root,tree->n_root->v[1],tree);

  /* printf("\n. ^ %f %f %f", */
  /* 	 (phydbl)(tree->rates->n_tips_below[tree->n_root->num]-2.), */
  /* 	 log(tree->rates->time_slice_lims[tree->rates->curr_slice[tree->n_root->num]] - */
  /* 	     tree->rates->nd_t[tree->n_root->num]), */
  /* 	 (phydbl)(tree->rates->n_tips_below[tree->n_root->num]-2.) * */
  /* 	 log(tree->rates->time_slice_lims[tree->rates->curr_slice[tree->n_root->num]] - */
  /* 	     tree->rates->nd_t[tree->n_root->num])); */

  tree->rates->c_lnL_times +=
    Factln(tree->rates->n_tips_below[tree->n_root->num]-2.) -
    (phydbl)(tree->rates->n_tips_below[tree->n_root->num]-2.) *
    log(tree->rates->time_slice_lims[tree->rates->curr_slice[tree->n_root->num]] -
  	tree->rates->nd_t[tree->n_root->num]);
  
  tree->rates->c_lnL_times -= logn;
  
  return(tree->rates->c_lnL_times);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void TIMES_Get_Number_Of_Time_Slices(t_tree *tree)
{
  int i;

  tree->rates->n_time_slices=0;
  TIMES_Get_Number_Of_Time_Slices_Post(tree->n_root,tree->n_root->v[2],tree);
  TIMES_Get_Number_Of_Time_Slices_Post(tree->n_root,tree->n_root->v[1],tree);
  Qksort(tree->rates->time_slice_lims,NULL,0,tree->rates->n_time_slices-1);

  if(tree->rates->n_time_slices > 1)
    {
      PhyML_Printf("\n");
      PhyML_Printf("\n. Sequences were collected at %d different time points.",tree->rates->n_time_slices);
      for(i=0;i<tree->rates->n_time_slices;i++) printf("\n+ [%3d] time point @ %12f ",i+1,tree->rates->time_slice_lims[i]);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void TIMES_Get_Number_Of_Time_Slices_Post(t_node *a, t_node *d, t_tree *tree)
{
  int i;

  if(d->tax == YES)
    {
      for(i=0;i<tree->rates->n_time_slices;i++) 
	if(Are_Equal(tree->rates->t_floor[d->num],tree->rates->time_slice_lims[i],1.E-6)) break;

      if(i == tree->rates->n_time_slices) 
	{
	  tree->rates->time_slice_lims[i] = tree->rates->t_floor[d->num];
	  tree->rates->n_time_slices++;
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
	  for(j=0;j<tree->rates->n_time_slices;j++)
	    {
	      if(Are_Equal(tree->rates->t_floor[i],tree->rates->time_slice_lims[j],1.E-6))
		{
		  tree->rates->n_time_slice_spans[i] = j+1;
		  /* PhyML_Printf("\n. Node %3d spans %3d slices [%12f].", */
		  /* 	       i+1, */
		  /* 	       tree->rates->n_slice_spans[i], */
		  /* 	       tree->rates->t_floor[i]); */
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
      
      if(tree->rates->curr_slice[a->num] != tree->rates->curr_slice[d->num])
	{
	  tree->rates->c_lnL_times += 
	    Factln(tree->rates->n_tips_below[d->num]-1.) - 
	    (phydbl)(tree->rates->n_tips_below[d->num]-1.) *
	    log(tree->rates->time_slice_lims[tree->rates->curr_slice[d->num]] -
		tree->rates->nd_t[d->num]);
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
	if(FABS(tree->rates->nd_t[left->bip_node[b->l_r][j]->num] - tree->rates->time_slice_lims[0]) < eps)
	  n_left_in++;
      
      n_left_out = left->bip_size[b->l_r]-n_left_in;
      
      n_rght_in = 0;
      For(j,rght->bip_size[b->r_l]) 
	if(FABS(tree->rates->nd_t[rght->bip_node[b->r_l][j]->num] - tree->rates->time_slice_lims[0]) < eps)
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


void Get_Survival_Duration(t_tree *tree)
{
  Get_Survival_Duration_Post(tree->n_root,tree->n_root->v[2],tree);
  Get_Survival_Duration_Post(tree->n_root,tree->n_root->v[1],tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Get_Survival_Duration_Post(t_node *a, t_node *d, t_tree *tree)
{
  if(d->tax)
    {
      tree->rates->survival_dur[d->num] = tree->rates->nd_t[d->num];
      return;
    }
  else
    {
      int i;
      t_node *v1, *v2;

      for(i=0;i<3;i++)
	if(d->v[i] != a && d->b[i] != tree->e_root)
	  Get_Survival_Duration_Post(d,d->v[i],tree);
      
      v1 = v2 = NULL;
      for(i=0;i<3;i++)
	{
	  if(d->v[i] != a && d->b[i] != tree->e_root)
	    {
	      if(!v1) v1 = d->v[i];
	      else    v2 = d->v[i];
	    }
	}

      tree->rates->survival_dur[d->num] = MAX(tree->rates->survival_dur[v1->num],
					      tree->rates->survival_dur[v2->num]);
    }
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

  For(i,2*tree->n_otu-1) tree->rates->t_rank[i] = i;

  t = tree->rates->nd_t;

  do
    {
      swap = NO;
      For(i,2*tree->n_otu-2)
	{
	  if(t[tree->rates->t_rank[i]] > t[tree->rates->t_rank[i+1]]) // Sort in ascending order
	    {
	      swap = YES;
	      buff                     = tree->rates->t_rank[i];
	      tree->rates->t_rank[i]   = tree->rates->t_rank[i+1];
	      tree->rates->t_rank[i+1] = buff;
	    }	    
	}
    }
  while(swap == YES);

  /* For(i,2*tree->n_otu-1) PhyML_Printf("\n. node %3d time: %12f", */
  /*                                     tree->rates->t_rank[i], */
  /*                                     tree->rates->nd_t[tree->rates->t_rank[i]]); */
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void TIMES_Label_Edges_With_Calibration_Intervals(t_tree *tree)
{
  char *s;
  int i;

  s = (char *)mCalloc(T_MAX_LINE,sizeof(char));
  
  tree->write_labels = YES;

  For(i,2*tree->n_otu-2)
    {
      if(tree->a_nodes[i]->tax == NO)
        {
          if(tree->rates->t_has_prior[i] == YES && tree->a_nodes[i]->b[0] != tree->e_root)
            {
              tree->a_nodes[i]->b[0]->n_labels = 1;
              Make_New_Edge_Label(tree->a_nodes[i]->b[0]);
              sprintf(s,"'>%f<%f'",FABS(tree->rates->t_prior_max[i])/100.,FABS(tree->rates->t_prior_min[i])/100.);
              tree->a_nodes[i]->b[0]->labels[0] = (char *)mCalloc(strlen(s)+1,sizeof(char));
              strcpy(tree->a_nodes[i]->b[0]->labels[0],s);
            }
        }
    }
  
  Free(s);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void TIMES_Set_Calibration(t_tree *tree)
{
  t_cal *cal;
  int i;

  For(i,2*tree->n_otu-1)
    {
      tree->rates->t_has_prior[i] = NO;
      tree->rates->t_prior_min[i] = BIG;
      tree->rates->t_prior_max[i] = BIG; 
   }

  cal = tree->rates->a_cal[0];
  while(cal)
    {
      /* if(cal->is_active == YES) */
      /*   { */
          /* tree->rates->t_has_prior[cal->node_num] = YES; */
          /* tree->rates->t_prior_min[cal->node_num] = cal->lower; */
          /* tree->rates->t_prior_max[cal->node_num] = cal->upper;           */
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
  For(i,2*tree->n_otu-1) 
    {
      tree->rates->t_prior_min_ori[i] = tree->rates->t_prior_min[i];
      tree->rates->t_prior_max_ori[i] = tree->rates->t_prior_max[i];
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void TIMES_Reset_Prior_Times(t_tree *tree)
{
  int i;
  For(i,2*tree->n_otu-1) 
    {
      tree->rates->t_prior_min[i] = tree->rates->t_prior_min_ori[i];
      tree->rates->t_prior_max[i] = tree->rates->t_prior_max_ori[i];
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

  tp_min = tree->rates->t_prior_min;
  tp_max = tree->rates->t_prior_max;
  tf = tree->rates->t_floor;
  t  = tree->rates->nd_t;
  n = NULL;
  loglbda = log(tree->rates->birth_rate);
  lbda = tree->rates->birth_rate;
  lower_bound = -1.;
  upper_bound = -1.;
  root_height = FABS(tree->rates->nd_t[tree->n_root->num]);

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
  int i;
  phydbl lnL;
  phydbl t,b,d,bmd,logbmd,expmbmd,logb;
  phydbl bmin,bmax;
  phydbl dmin,dmax;

  lnL     = 0.0;
  b       = tree->rates->birth_rate;
  d       = tree->rates->death_rate;
  bmin    = tree->rates->birth_rate_min;
  bmax    = tree->rates->birth_rate_max;
  dmin    = tree->rates->death_rate_min;
  dmax    = tree->rates->death_rate_max;
  bmd     = b-d;
  logbmd  = -1.;
  expmbmd = -1.;
  logb    = -1.;
  t       = 0.0;

  if(b < d)
    {
      if(verbose) printf("\n. b: %G d: %G",b,d);
      return UNLIKELY;
    }
  
  For(i,2*tree->n_otu-1)
    if(tree->a_nodes[i]->tax == NO)
      {
        if(tree->rates->nd_t[i] < tree->rates->t_prior_min[i] ||
           tree->rates->nd_t[i] > tree->rates->t_prior_max[i]) 
          {
            if(verbose) printf("\n. node %d @ time %f min: %f max: %f",i,tree->rates->nd_t[i],tree->rates->t_prior_min[i],tree->rates->t_prior_max[i]);
            return UNLIKELY;
          }
      }
  
  if(b > bmin && d > dmin && Are_Equal(bmd,0.0,bmin/10.) == NO)
    {
      logbmd = log(bmd);
      expmbmd = exp(d-b);

      For(i,2*tree->n_otu-1)
        if(tree->a_nodes[i]->tax == NO && tree->a_nodes[i] != tree->n_root)
          {
            t = FABS(tree->rates->nd_t[i]);            
            // Equation 3.19 in Tanja Stadler's PhD thesis
            lnL += 2*logbmd - bmd*t - 2.*log(b-d*POW(expmbmd,t));            
          }

      t = FABS(tree->rates->nd_t[tree->n_root->num]);
      lnL += -bmd*t - log(b-d*POW(expmbmd,t));

      lnL += log(bmd) + (tree->n_otu-1)*log(b) + LnGamma(tree->n_otu+1);
    }
  else if(b > bmin && d < dmin) // Yule process 
    {
      For(i,2*tree->n_otu-1)
        if(tree->a_nodes[i]->tax == NO && tree->a_nodes[i] != tree->n_root)
          {
            t = FABS(tree->rates->nd_t[i]);            
            // Equation 3.19 in Tanja Stadler's PhD thesis (Yule case)
            lnL -= b*t;
          }

      t = FABS(tree->rates->nd_t[tree->n_root->num]);
      lnL -= b*t;
      
      lnL += (tree->n_otu-1)*log(b) + LnGamma(tree->n_otu+1);
    }
  else if(b < bmin && d > dmin) 
    {
      if(verbose) printf("\n. b: %G bmin: %G d: %G dmin: %G",b,bmin,d,dmin);
      return UNLIKELY;
    }
  else if(Are_Equal(bmd,0.0,bmin/10.) == YES) // Critical birth-death process
    {
      logb = log(b);

      For(i,2*tree->n_otu-1)
        if(tree->a_nodes[i]->tax == NO && tree->a_nodes[i] != tree->n_root)
          {
            t = FABS(tree->rates->nd_t[i]);            
            // Equation 3.19 in Tanja Stadler's PhD thesis (Critical case)
            lnL += logb - 2.*log(1.+b*t);
          }

      t = FABS(tree->rates->nd_t[tree->n_root->num]);
      lnL -= log(b*t);

      lnL += LnGamma(tree->n_otu+1);
      
    }
  else if(b < bmin && d < dmin) // Birth and death rates are below their limits
    {
      if(verbose) printf("\n. b: %G bmin: %G d: %G dmin: %G",b,bmin,d,dmin);
      return -INFINITY;
    }
  else if(b > bmax && d > dmax)
    {
      if(verbose) printf("\n. b: %G bmax: %G d: %G dmax: %G",b,bmax,d,dmax);
      return -INFINITY;
    }
  else
    {
      assert(FALSE);
    }

  if(isnan(lnL) || isinf(FABS(lnL)))
    {
      if(verbose) printf("\n. lnL: %f",lnL);
      tree->rates->c_lnL_times = UNLIKELY;
      return UNLIKELY;
    }

  tree->rates->c_lnL_times = lnL;


  return(lnL);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

// Generate a subtree including all taxa in tax_list. The age of the root of that
// subtree is t_mrca. All nodes in the subtree are thus younger than that.
void TIMES_Connect_List_Of_Taxa(t_node **tax_list, int list_size, phydbl t_mrca, phydbl *times, int *nd_num, t_tree *mixt_tree)
{
  phydbl t_upper_bound, t_lower_bound,up,lo;
  int i,j,n_anc,*permut;
  t_node *n,**anc,*new_mrca;

  t_lower_bound = t_mrca;
  t_upper_bound = 0.0;
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
  
  /* printf("\n. upper: %f lower: %f t_mrca: %f\n",t_upper_bound,t_lower_bound,t_mrca); */
  assert(t_upper_bound > t_lower_bound);
  
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
  /* permut = Permutate(n_anc); */
  permut = (int *)mCalloc(n_anc,sizeof(int));
  for(i=0;i<n_anc;i++) permut[i] = i;
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
      lo                     = up - (up - t_mrca)/10.;
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
      i++;
      (*nd_num) += 1;
      if(n_anc == i+1) { times[new_mrca->num] = t_mrca; break; }
    }
  while(1);
  
  Free(permut);
  Free(anc);
}
  
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// Generate a  random rooted tree with node ages fullfiling the
// time constraints defined in the list of calibration cal_list 
void TIMES_Randomize_Tree_With_Time_Constraints(t_cal *cal_list, t_tree *mixt_tree)
{
  t_node **tips,**nd_list;
  phydbl *times,*cal_times,time_oldest_cal;
  int i,j,nd_num,*cal_ordering,n_cal,swap,list_size,tmp,orig_is_mixt_tree,repeat,n_max_repeats,tip_num,*no_cal_tip_num,n_no_cal_tip_num,*permut;
  t_cal *cal;
  t_clad *clade;
  

  
  assert(mixt_tree->rates);
  
  tips    = (t_node **)mCalloc(mixt_tree->n_otu,sizeof(t_node *));
  nd_list = (t_node **)mCalloc(mixt_tree->n_otu,sizeof(t_node *));
  
  times                   = mixt_tree->rates->nd_t;
  orig_is_mixt_tree       = mixt_tree->is_mixt_tree;
  mixt_tree->is_mixt_tree = NO;
  n_max_repeats           = 1000;
  cal                     = NULL;
  clade                   = NULL;
  
  // List node indices that are not in any calibration set
  no_cal_tip_num = NULL;
  n_no_cal_tip_num = 0;
  for(i=0;i<mixt_tree->n_otu;++i)
    {
      cal = cal_list;
      while(cal != NULL)
        {
          clade = cal->clade_list[cal->current_clade_idx];
          for(j=0;j<clade->n_tax;++j)
            {
              /* printf("\n>>> %s %s %d [%d %d] [%p %p]", */
              /*        clade->tip_list[j]->name, */
              /*        mixt_tree->a_nodes[i]->name, */
              /*        clade->tip_list[j] == mixt_tree->a_nodes[i], */
              /*        clade->tip_list[j]->num, */
              /*        mixt_tree->a_nodes[i]->num, */
              /*        clade->tip_list[j], */
              /*        mixt_tree->a_nodes[i]); */              
              if(clade->tip_list[j] == mixt_tree->a_nodes[i]) break;
            }
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
  
  
  for(repeat=0;repeat<n_max_repeats;repeat++)
    {
      nd_num = mixt_tree->n_otu;
  
      /* PhyML_Printf("\n\n. Repeat %d",repeat); */

      for(i=0;i<mixt_tree->n_otu;i++) tips[i] = mixt_tree->a_nodes[i];
      for(i=0;i<mixt_tree->n_otu;i++) mixt_tree->a_nodes[i]->v[0] = NULL; 


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
 
      /* printf("\n. n_cal : %d",n_cal); fflush(NULL); */
      /* Exit("\n"); */
     
      cal_ordering = (int *)mCalloc(n_cal,sizeof(int));
      for(i=0;i<n_cal;i++) cal_ordering[i] = i;
      
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
      
      for(i=0;i<n_cal-1;i++) assert(cal_times[cal_ordering[i]] > cal_times[cal_ordering[i+1]]);
              
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
                  assert(list_size < mixt_tree->n_otu);
                }
            }
          while(j++ < MIN(n_no_cal_tip_num-1,(int)(2*n_no_cal_tip_num/(n_cal))));
          Free(permut);

          

          // Add all the taxa that are  in the calibration set of cal
          // This should be done here so that the last node in nd_list
          // belongs to the taxa in the calibration
          clade = cal->clade_list[cal->current_clade_idx];
          for(j=0;j<clade->n_tax;j++)
            {
              nd_list[list_size] = tips[clade->tip_list[j]->num];
              list_size++;
              assert(list_size < mixt_tree->n_otu);
            }
          

          
          /* for(j=0;j<n_cal;j++) */
          /*   { */
          /*     printf("\n. %d -> %f", */
          /*            cal_ordering[j], */
          /*            cal_times[cal_ordering[j]]); */
          /*     fflush(NULL); */
          /*   } */

          /* for(j=0;j<list_size;j++) PhyML_Printf("\n. %s [%d]",nd_list[j]->name,nd_list[j]->num); */
          /* PhyML_Printf("\n. Time: %f [%f %f]\n", */
          /*              cal_times[cal_ordering[i]], */
          /*              cal->lower, */
          /*              cal->upper); */
          
          /* for(k=0;k<list_size;k++) printf("\n@ %s",nd_list[k]->name); */
          /* printf("\n"); */

          TIMES_Connect_List_Of_Taxa(nd_list, 
                                     list_size,
                                     cal_times[cal_ordering[i]], 
                                     times, 
                                     &nd_num,
                                     mixt_tree);

        }
      
      Free(cal_times);
      Free(cal_ordering);
      
      
      // Connect all remaining taxa
      for(i=0;i<mixt_tree->n_otu;i++) nd_list[i] = NULL;
      list_size = 0;
      for(i=0;i<mixt_tree->n_otu;++i)
        {
          if(tips[i]->v[0] == NULL) // Tip is not connected yet
            {
              nd_list[list_size] = tips[i];
              list_size++;
            }
        }


      
      cal = cal_list;
      do
        {
          clade = cal->clade_list[cal->current_clade_idx];
          nd_list[list_size] = clade->tip_list[0];
          list_size++;
          cal = cal->next;
        }
      while(cal);

      /* for(i=0;i<list_size;i++) printf("\n# To connect: %d lower: %f",nd_list[i]->num,time_oldest_cal); */
      
      TIMES_Connect_List_Of_Taxa(nd_list,
                                 list_size,
                                 1.5*time_oldest_cal,
                                 times,
                                 &nd_num,
                                 mixt_tree);
      
            
      // Adding root node 
      mixt_tree->n_root = mixt_tree->a_nodes[2*mixt_tree->n_otu-2];
      mixt_tree->n_root->v[1]->v[0] = mixt_tree->n_root->v[2];
      mixt_tree->n_root->v[2]->v[0] = mixt_tree->n_root->v[1];
      Update_Ancestors(mixt_tree->n_root,mixt_tree->n_root->v[2],mixt_tree);
      Update_Ancestors(mixt_tree->n_root,mixt_tree->n_root->v[1],mixt_tree);
      mixt_tree->n_root->anc = NULL;

      // Adding root edge
      mixt_tree->num_curr_branch_available = 0;
      Connect_Edges_To_Nodes_Recur(mixt_tree->a_nodes[0],mixt_tree->a_nodes[0]->v[0],mixt_tree);
      
      for(i=0;i<2*mixt_tree->n_otu-3;++i)
        {
          if(((mixt_tree->a_edges[i]->left == mixt_tree->n_root->v[1]) || (mixt_tree->a_edges[i]->rght == mixt_tree->n_root->v[1])) &&
             ((mixt_tree->a_edges[i]->left == mixt_tree->n_root->v[2]) || (mixt_tree->a_edges[i]->rght == mixt_tree->n_root->v[2])))
            {
              Add_Root(mixt_tree->a_edges[i],mixt_tree);
              break;
            }
        }
      
      
      DATE_Assign_Primary_Calibration(mixt_tree);
      DATE_Update_T_Prior_MinMax(mixt_tree);

      /* { */
      /*   Print_Node(mixt_tree->n_root,mixt_tree->n_root->v[1],mixt_tree); */
      /*   Print_Node(mixt_tree->n_root,mixt_tree->n_root->v[2],mixt_tree); */
      /*   fflush(NULL); */

      /*   int i; */
      /*   for(i=0;i<mixt_tree->rates->n_cal;i++) */
      /*     { */
      /*       PhyML_Printf("\n. Node number to which calibration [%d] applies to is [%d]",i,Find_Clade(mixt_tree->rates->a_cal[i]->target_tax, */
      /*                                                                                                mixt_tree->rates->a_cal[i]->n_target_tax, */
      /*                                                                                                mixt_tree)); */
      /*       PhyML_Printf("\n. Lower bound set to: %15f time units.",mixt_tree->rates->a_cal[i]->lower); */
      /*       PhyML_Printf("\n. Upper bound set to: %15f time units.",mixt_tree->rates->a_cal[i]->upper); */
      /*     } */
      /* } */

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
      Exit("\n");
    }

  /* for(j=0;j<mixt_tree->n_otu;j++) printf("\n. %s",mixt_tree->a_nodes[j]->name); */
        
  assert(i != 2*mixt_tree->n_otu-3);

  mixt_tree->is_mixt_tree = orig_is_mixt_tree;


  
  if(mixt_tree->is_mixt_tree == YES)
    {
      t_tree *tree;
      
      // Propagate tree topology and reorganize partial lk struct along edges
      tree = mixt_tree->next;
      do
        {
          For(i,2*tree->n_otu-1)
            {
              tree->a_nodes[i]->v[0] = tree->prev->a_nodes[i]->v[0] ? tree->a_nodes[tree->prev->a_nodes[i]->v[0]->num] : NULL;
              tree->a_nodes[i]->v[1] = tree->prev->a_nodes[i]->v[1] ? tree->a_nodes[tree->prev->a_nodes[i]->v[1]->num] : NULL;
              tree->a_nodes[i]->v[2] = tree->prev->a_nodes[i]->v[2] ? tree->a_nodes[tree->prev->a_nodes[i]->v[2]->num] : NULL;
            }
          tree->num_curr_branch_available = 0;
          Connect_Edges_To_Nodes_Recur(tree->a_nodes[0],tree->a_nodes[0]->v[0],tree);
          Add_Root(tree->a_edges[tree->prev->e_root->num],tree);
          Reorganize_Edges_Given_Lk_Struct(tree);
          Init_Partial_Lk_Tips_Int(tree);
                   
          tree = tree->next;
        }
      while(tree);
    }

  
  Free(tips);
  Free(nd_list);
  Free(no_cal_tip_num);
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
  if(tree->rates->nd_t[d->num] < tree->rates->nd_t[a->num])
    {
      PhyML_Printf("\n. a->t = %f [num:%d] d->t %f [num:%d]",
                   tree->rates->nd_t[a->num],
                   a->num,
                   tree->rates->nd_t[d->num],
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
