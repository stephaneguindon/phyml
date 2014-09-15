/*
 * author: Imran Fanaswala
 */

#ifndef BEAGLE_UTILS_CPP
#define BEAGLE_UTILS_CPP

#include  <stdio.h>
#include "beagle_utils.h"

#define CLEAN_BEAGLE_API  /* Attempting to remove unnecessary communication with  BEAGLE device */

double* int_to_double(const int* src, int num_elems)
{
    double* dest = (double*)malloc(num_elems*sizeof(double)); if (NULL==dest) Warn_And_Exit("\n. Couldn't allocate memory.\n");
    int i;
    for(i=0; i<num_elems;++i)
        dest[i] = (double)src[i];
    return dest;
}

double* short_to_double(const short* src, int num_elems)
{
    double* dest = (double*)malloc(num_elems*sizeof(double)); if (NULL==dest) Warn_And_Exit("\n. Couldn't allocate memory. \n");
    int i;
    for(i=0; i<num_elems;++i)
        dest[i] = (double)src[i];
    return dest;
}

double* float_to_double(const phydbl* src, int num_elems)
{
    double* dest = (double*)malloc(num_elems*sizeof(double)); if (NULL==dest) Warn_And_Exit("\n. Couldn't allocate memory. \n");
    int i;
    for(i=0; i<num_elems;++i)
        dest[i] = (double)src[i];
    return dest;
}

void print_beagle_flags(long inFlags) {
    if (inFlags & BEAGLE_FLAG_PROCESSOR_CPU)      fprintf(stdout, " PROCESSOR_CPU");
    if (inFlags & BEAGLE_FLAG_PROCESSOR_GPU)      fprintf(stdout, " PROCESSOR_GPU");
    if (inFlags & BEAGLE_FLAG_PROCESSOR_FPGA)     fprintf(stdout, " PROCESSOR_FPGA");
    if (inFlags & BEAGLE_FLAG_PROCESSOR_CELL)     fprintf(stdout, " PROCESSOR_CELL");
    if (inFlags & BEAGLE_FLAG_PRECISION_DOUBLE)   fprintf(stdout, " PRECISION_DOUBLE");
    if (inFlags & BEAGLE_FLAG_PRECISION_SINGLE)   fprintf(stdout, " PRECISION_SINGLE");
    if (inFlags & BEAGLE_FLAG_COMPUTATION_ASYNCH) fprintf(stdout, " COMPUTATION_ASYNCH");
    if (inFlags & BEAGLE_FLAG_COMPUTATION_SYNCH)  fprintf(stdout, " COMPUTATION_SYNCH");
    if (inFlags & BEAGLE_FLAG_EIGEN_REAL)         fprintf(stdout, " EIGEN_REAL");
    if (inFlags & BEAGLE_FLAG_EIGEN_COMPLEX)      fprintf(stdout, " EIGEN_COMPLEX");
    if (inFlags & BEAGLE_FLAG_SCALING_MANUAL)     fprintf(stdout, " SCALING_MANUAL");
    if (inFlags & BEAGLE_FLAG_SCALING_AUTO)       fprintf(stdout, " SCALING_AUTO");
    if (inFlags & BEAGLE_FLAG_SCALING_ALWAYS)     fprintf(stdout, " SCALING_ALWAYS");
    if (inFlags & BEAGLE_FLAG_SCALING_DYNAMIC)    fprintf(stdout, " SCALING_DYNAMIC");
    if (inFlags & BEAGLE_FLAG_SCALERS_RAW)        fprintf(stdout, " SCALERS_RAW");
    if (inFlags & BEAGLE_FLAG_SCALERS_LOG)        fprintf(stdout, " SCALERS_LOG");
    if (inFlags & BEAGLE_FLAG_VECTOR_NONE)        fprintf(stdout, " VECTOR_NONE");
    if (inFlags & BEAGLE_FLAG_VECTOR_SSE)         fprintf(stdout, " VECTOR_SSE");
    if (inFlags & BEAGLE_FLAG_VECTOR_AVX)         fprintf(stdout, " VECTOR_AVX");
    if (inFlags & BEAGLE_FLAG_THREADING_NONE)     fprintf(stdout, " THREADING_NONE");
    if (inFlags & BEAGLE_FLAG_THREADING_OPENMP)   fprintf(stdout, " THREADING_OPENMP");
    if (inFlags & BEAGLE_FLAG_FRAMEWORK_CPU)      fprintf(stdout, " FRAMEWORK_CPU");
    if (inFlags & BEAGLE_FLAG_FRAMEWORK_CUDA)     fprintf(stdout, " FRAMEWORK_CUDA");
    if (inFlags & BEAGLE_FLAG_FRAMEWORK_OPENCL)   fprintf(stdout, " FRAMEWORK_OPENCL");
    fflush(stdout);
}

void print_beagle_resource_list()
{
  int i;
  BeagleResourceList* rList;
  rList = beagleGetResourceList();
  fprintf(stdout, "\n\tAvailable resources:\n");
  for (i = 0; i < rList->length; i++) {
    fprintf(stdout, "\t\tResource %i:\n\t\tName : %s\n", i, rList->list[i].name);
    fprintf(stdout, "\t\t\tDesc : %s\n", rList->list[i].description);
    fprintf(stdout, "\t\t\tFlags:");
    print_beagle_flags(rList->list[i].supportFlags);
    fprintf(stdout, "\n");
  }
  fflush(stdout);
}

void print_beagle_instance_details(BeagleInstanceDetails *inst)
{
    int rNumber = inst->resourceNumber;
    fprintf(stdout, "\tUsing resource %i:\n", rNumber);
    fprintf(stdout, "\t\tRsrc Name : %s\n", inst->resourceName);
    fprintf(stdout, "\t\tImpl Name : %s\n", inst->implName);
    fprintf(stdout, "\t\tImpl Desc : %s\n", inst->implDescription);
    fprintf(stdout, "\t\tFlags:");
    fflush(stdout);
    print_beagle_flags(inst->flags);
    fflush(stdout);
}

int create_beagle_instance(t_tree *tree, int quiet, option* io)
{
    if(UNINITIALIZED != tree->b_inst){
//        fprintf(stdout,"\n\tWARNING: Creating a BEAGLE instance on a tree with an existing BEAGLE instance:%d\n",tree->b_inst);
    }
    if(!quiet){
//        print_beagle_resource_list();
    }
    int i;
    BeagleInstanceDetails inst_d;
    int num_rate_catg = tree->mod->ras->n_catg;
    int num_branches  = 2*tree->n_otu-1; //rooted tree
    //Recall that in PhyML, each edge has a "left" and "right" partial vectors. Therefore,
    //in BEAGLE we have 2*num_branches number of partials.
    //BEAGLE's partials buffer = [ tax1, tax2, ..., taxN, b1Left, b2Left, b3Left,...,bMLeft, b1Rght, b2Rght, b3Rght,...,bMRght] (N taxa, M branches)
    int num_partials  = 2*(tree->n_otu + num_branches); /* TODO: This does not seem correct; suspect poor indexing elsewhere */
    													/* In update_operations, indexes range from 0 to almost 2 (n_otu + num_branches), but not all integers are used */
        
    int resourceList[1];    
    resourceList[0] = io->beagle_resource;        
        
//    DUMP_I(tree->n_otu, num_rate_catg, num_partials, num_branches, tree->mod->ns, tree->n_pattern, tree->mod->whichmodel);
    int beagle_inst = beagleCreateInstance(
                                  tree->n_otu,                /**< Number of tip data elements (input) */
                                  num_partials,               /**< Number of partial buffer (input) */
    /* PhyML uses partials */     0,              			  /**< Number of compact state representation buffers to create (input) */
                                  tree->mod->ns,              /**< Number of states in the continuous-time Markov chain (input) */
                                  tree->n_pattern,            /**< Number of site patterns to be handled by the instance (input) */
                                  1,                          /**< Number of rate matrix eigen-decomposition,state freqs, and category weight buffers*/
                                  num_branches,               /**< Number of rate matrix buffers (input) */
                                  num_rate_catg,              /**< Number of rate categories (input) */
                                  -1,                         /**< Number of scaling buffers. Unused because we use SCALING_ALWAYS */
                                  resourceList,               /**< List of potential resource on which this instance is allowed (input, NULL implies no restriction */
                                  1,			              /**< Length of resourceList list (input) */
                                  BEAGLE_FLAG_FRAMEWORK_CPU | BEAGLE_FLAG_PROCESSOR_CPU | BEAGLE_FLAG_SCALING_ALWAYS | BEAGLE_FLAG_EIGEN_REAL | ((sizeof(float)==sizeof(phydbl)) ? BEAGLE_FLAG_PRECISION_SINGLE:BEAGLE_FLAG_PRECISION_DOUBLE),
                                  0,                		  /**< Bit-flags indicating required implementation characteristics, see BeagleFlags (input) */
                                  &inst_d);
    if (beagle_inst < 0){
        fprintf(stderr, "beagleCreateInstance() failed:%i\n\n",beagle_inst);
        return beagle_inst;
    }

    if(!quiet){
        fprintf(stdout, "\n\tUnique BEAGLE instance id:%i\n", beagle_inst);
        print_beagle_instance_details(&inst_d);
    }

    //Set the tips
    for(i=0; i<2*tree->n_otu-1; ++i) //taxa+internal nodes
      {
        //        Print_Tip_Partials(tree, tree->a_nodes[i]);
        if(tree->a_nodes[i]->tax)
          {
            assert(tree->a_nodes[i]->c_seq->len == tree->n_pattern); // number of compacts sites == number of distinct site patterns
            double* tip = short_to_double(tree->a_nodes[i]->b[0]->p_lk_tip_r, tree->n_pattern*tree->mod->ns); //The tip states are stored on the branch leading to the tip
            //Recall we store tip partials on the branch leading to the tip, rather than the tip itself.
            int ret = beagleSetTipPartials(beagle_inst, tree->a_nodes[i]->b[0]->p_lk_tip_idx, tip);
            if(ret<0){
              fprintf(stderr, "beagleSetTipPartials() on instance %i failed:%i\n\n",beagle_inst,ret);
              Free(tip);
              return ret;
            }
            Free(tip);
          }
      }
    
    //Set the pattern weights
    double* pwts = int_to_double(tree->data->wght,tree->n_pattern); //BTW, These weights are absolute counts, and not freqs
    int ret = beagleSetPatternWeights(beagle_inst, pwts);
    if(ret<0){
      fprintf(stderr, "beagleSetPatternWeights() on instance %i failed:%i\n\n",beagle_inst,ret);
      Free(pwts);
      return ret;
    }
    Free(pwts);
    
    tree->mod->b_inst=beagle_inst;
    
    update_beagle_ras(tree->mod);
    update_beagle_efrqs(tree->mod);
    
    return beagle_inst;
}

/* Update partial likelihood on edge b on the side of b where
   node d lies.
*/
void update_beagle_partials(t_tree* tree, t_edge* b, t_node* d)
{
    /*
               |
               |<- b
               |
               d
              / \
          b1 /   \ b2
            /     \
        n_v1     n_v2
    */

  if(d->tax) //Partial likelihoods are only calculated on internal nodes
    {
      PhyML_Printf("\n== t_node %d is a leaf...",d->num);
      PhyML_Printf("\n== Err. in file %s at line %d (function '%s')\n",__FILE__,__LINE__,__FUNCTION__);
      Warn_And_Exit("\n");
    }

  //Determine d's "left" and "right" neighbors.
  t_node *n_v1, *n_v2;//d's "left" and "right" neighbor nodes
  phydbl *p_lk,*p_lk_v1,*p_lk_v2;
  phydbl *Pij1,*Pij2;
  int *sum_scale, *sum_scale_v1, *sum_scale_v2;
  int *p_lk_loc;
  int dest_p_idx, child1_p_idx, child2_p_idx, Pij1_idx, Pij2_idx;
  n_v1 = n_v2                 = NULL;
  p_lk = p_lk_v1 = p_lk_v2    = NULL;
  Pij1 = Pij2                 = NULL;
  sum_scale_v1 = sum_scale_v2 = NULL;
  p_lk_loc                    = NULL;
  dest_p_idx = child1_p_idx = child2_p_idx = Pij1_idx = Pij2_idx = UNINITIALIZED;
  Set_All_P_Lk(&n_v1,&n_v2,
               &p_lk,&sum_scale,&p_lk_loc,
               &Pij1,&p_lk_v1,&sum_scale_v1,
               &Pij2,&p_lk_v2,&sum_scale_v2,
               d,b,tree,
               &dest_p_idx, &child1_p_idx, &child2_p_idx, &Pij1_idx, &Pij2_idx);
  
  
  //    fprintf(stdout, "\nUpdating partials on Branch %d (on the side where Node %d lies)\n",b->num,d->num);fflush(stdout);
  //    double* p_lk_v1_b = (double*)malloc(tree->mod->ras->n_catg*tree->mod->ns*tree->n_pattern*sizeof(double));if(NULL==p_lk_v1_b) Warn_And_Exit("Couldnt allocate memory");
  //    beagleGetPartials(tree->b_inst, child1_p_idx, BEAGLE_OP_NONE, (double*)p_lk_v1_b);
  //    double* p_lk_v2_b = (double*)malloc(tree->mod->ras->n_catg*tree->mod->ns*tree->n_pattern*sizeof(double));if(NULL==p_lk_v2_b) Warn_And_Exit("Couldnt allocate memory");
  //    beagleGetPartials(tree->b_inst, child2_p_idx, BEAGLE_OP_NONE, (double*)p_lk_v2_b);
  
  //    fprintf(stdout, "Left partials :");fflush(stdout);
  //    Dump_Arr_D(p_lk_v1_b,   tree->mod->ras->n_catg*tree->mod->ns*tree->n_pattern);
  //    fprintf(stdout, "Right partials:");fflush(stdout);
  //    Dump_Arr_D(p_lk_v2_b,   tree->mod->ras->n_catg*tree->mod->ns*tree->n_pattern);
  //    Free(p_lk_v1_b);
  //    Free(p_lk_v2_b);
  
  
  //Create the corresponding BEAGLE operation
  
  // fprintf(stderr,"%d, %d, %d, ", dest_p_idx, child1_p_idx, child2_p_idx);
  
  BeagleOperation operations[1] = {{dest_p_idx, BEAGLE_OP_NONE, BEAGLE_OP_NONE, child1_p_idx, Pij1_idx, child2_p_idx, Pij2_idx}};
  //Compute the partials
  int ret = beagleUpdatePartials(tree->b_inst, operations, 1, BEAGLE_OP_NONE);
  if(ret<0){
    fprintf(stderr, "beagleUpdatePartials() on instance %i failed:%i\n\n",tree->b_inst,ret);
    Exit("");
  }
  //Load the computed/updated partial partials
#ifndef CLEAN_BEAGLE_API
   ret = beagleGetPartials(tree->b_inst, dest_p_idx, BEAGLE_OP_NONE, (double*)p_lk);
   if(ret<0){
     fprintf(stderr, "beagleGetPartials() on instance %i failed:%i\n\n",tree->b_inst,ret);
     Exit("");
   }
#endif
  
  //    fprintf(stdout, "Updated partials:");fflush(stdout);
  //    Dump_Arr_D(p_lk, tree->mod->ras->n_catg*tree->mod->ns*tree->n_pattern);
}

int finalize_beagle_instance(t_tree *tree)
{
  if(tree->b_inst >= 0) {
    int ret = beagleFinalizeInstance(tree->b_inst);
    if(ret<0) fprintf(stderr, "\nFailed to finalize BEAGLE instance %i: %i\n\n", tree->b_inst, ret);
    return ret;
  }
  return 0;
}

void update_beagle_ras(t_mod* mod)
{
  assert(UNINITIALIZED != mod->b_inst);
  
  int ret=-1;
  if((sizeof(float)==sizeof(phydbl))) //Do we need to convert?
    {
      double* catg_rates = float_to_double(mod->ras->gamma_rr->v, mod->ras->n_catg);
      ret = beagleSetCategoryRates(mod->b_inst, catg_rates);
      Free(catg_rates);
      if(ret<0){
        fprintf(stderr, "beagleSetCategoryRates() on instance %i failed:%i\n\n",mod->b_inst,ret);
        Exit("");
      }
      double* catg_wts = float_to_double(mod->ras->gamma_r_proba->v, mod->ras->n_catg);
      if(!mod->optimizing_topology) {
        ret = beagleSetCategoryWeights(mod->b_inst, 0, catg_wts);
        Free(catg_wts);
        if(ret<0){
          fprintf(stderr, "beagleSetCategoryWeights() on instance %i failed:%i\n\n",mod->b_inst,ret);
          Exit("");
        }
      }
    }
  else
    {
      ret = beagleSetCategoryRates(mod->b_inst, mod->ras->gamma_rr->v);
      if(ret<0){
        fprintf(stderr, "beagleSetCategoryRates() on instance %i failed:%i\n\n",mod->b_inst,ret);
        Exit("");
      }
      if(!mod->optimizing_topology) {
        ret = beagleSetCategoryWeights(mod->b_inst, 0, mod->ras->gamma_r_proba->v);
        if(ret<0){
          fprintf(stderr, "beagleSetCategoryWeights() on instance %i failed:%i\n\n",mod->b_inst,ret);
          Exit("");
        }
        
      }
    }
}

void update_beagle_efrqs(t_mod* mod)
{
  assert(UNINITIALIZED != mod->b_inst);
  
  int ret=-1;
  if((sizeof(float)==sizeof(phydbl))) {
    double* efrqs = float_to_double(mod->e_frq->pi->v, mod->ns);
        ret = beagleSetStateFrequencies(mod->b_inst, 0, efrqs);
        Free(efrqs);
  } else {
    ret = beagleSetStateFrequencies(mod->b_inst, 0, mod->e_frq->pi->v);
  }
  if(ret<0){
    fprintf(stderr, "beagleSetStateFrequencies() on instance %i failed:%i\n\n",mod->b_inst,ret);
    Exit("");
  }
}

void calc_edgelks_beagle(t_edge *b, t_tree *tree)
{
  assert(UNINITIALIZED != tree->b_inst);
  
  //Compute the edge likelihood
  int parents[1]  = {b->p_lk_left_idx};
  int children[1] = {b->rght->tax?b->p_lk_tip_idx:b->p_lk_rght_idx};
  int pmats[1]    = {b->Pij_rr_idx};
  int other[1]    = {0};//Category Weights and State Frequencies both have a single buffer, hence they are both indexed at 0
  double lnL[1]   = {UNINITIALIZED};
  //    DUMP_I(parents[0], children[0], pmats[0], b->num, b->rght->tax);
  int ret=beagleCalculateEdgeLogLikelihoods(tree->b_inst, parents, children, pmats, NULL, NULL, other, other, NULL, 1, lnL, NULL, NULL);
  int i;

if(ret<0){
    fprintf(stderr, "beagleCalculateEdgeLogLikelihoods() on instance %i failed:%i\n\n",tree->b_inst,ret);
    Exit("");
  }
  
  tree->c_lnL = sizeof(phydbl)==sizeof(float)?(float)lnL[0]:lnL[0];
  
  //Retrieve the site likelihoods that were computed during the previous edge likelihood computation
  ret = beagleGetSiteLogLikelihoods(tree->b_inst,tree->cur_site_lk);//TODO: Handle when cur_site_lk is float
  if(ret<0){
    fprintf(stderr, "beagleGetSiteLogLikelihoods() on instance %i failed:%i\n\n",tree->b_inst,ret);
    Exit("");
  }
  //Transform
  for(i=0;i<tree->n_pattern;++i)
    tree->cur_site_lk[i]=EXP(tree->cur_site_lk[i]);
}

void update_beagle_eigen(t_mod* mod)
{
    assert(UNINITIALIZED != mod->b_inst);

    int whichmodel = mod->whichmodel;
    //We use Eigen Decomposition only for GTR models and AA datasets
    if((mod->io->datatype == AA || whichmodel==GTR || whichmodel==CUSTOM) && mod->use_m4mod == NO)
    {
        //Beagle expects untransformed eigen-values (i.e. recall e_val is EXP() scaled, so we undo that)
        phydbl* evals = (phydbl*)malloc(mod->ns * sizeof(phydbl));
        int i;
        for(i=0;i<mod->ns;++i)  evals[i]=LOG(mod->eigen->e_val[i]);
        int ret=-1;
        if((sizeof(float)==sizeof(phydbl)))//Need to convert to doubles?
        {
            double* eigen_vects     = float_to_double(mod->eigen->r_e_vect, mod->eigen->size*mod->eigen->size);
            double* eigen_vects_inv = float_to_double(mod->eigen->l_e_vect, mod->eigen->size*mod->eigen->size);
            double* eigen_vals      = float_to_double(evals,mod->eigen->size);
            ret = beagleSetEigenDecomposition(mod->b_inst,0,eigen_vects,eigen_vects_inv,eigen_vals);
            Free(eigen_vects);Free(eigen_vects_inv);Free(eigen_vals);
        } else {
            ret = beagleSetEigenDecomposition(mod->b_inst,0,mod->eigen->r_e_vect,mod->eigen->l_e_vect,evals);
        }
        if(ret<0){
          fprintf(stderr, "beagleSetEigenDecomposition() on instance %i failed:%i\n\n",mod->b_inst,ret);
          Free(evals);
          Exit("");
        }
        Free(evals);
    }
}

#endif // BEAGLE_UTILS_CPP

