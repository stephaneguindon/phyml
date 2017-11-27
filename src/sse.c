/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include "assert.h"
#include "sse.h"

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

#if defined(__SSE3__)

void SSE_Update_Eigen_Lr(t_edge *b, t_tree *tree)
{
  unsigned int site,catg,state;
  int is_ambigu,observed_state;
  unsigned int i,j,k,l;
  
  unsigned const int npattern = tree->n_pattern;
  unsigned const int ncatg = tree->mod->ras->n_catg;
  unsigned const int ns = tree->mod->ns;
  unsigned const int sz = (int)BYTE_ALIGN / 8;  
  unsigned const int nscatg = ncatg * ns;
  unsigned const int nblocks = ns / sz;

  __m128d _fplk[nblocks],_fplkev[sz],_colsum[nblocks];
  phydbl *ev;
  
  observed_state = -1;
  
  ev = NULL;
#ifndef WIN32
  if(posix_memalign((void **)&ev,BYTE_ALIGN,(size_t)sz*sizeof(double))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
#else
  ev = _aligned_malloc(sz * sizeof(phydbl),BYTE_ALIGN);
#endif
  
  for(site=0;site<npattern;++site)
    {
      is_ambigu = YES;
      if(b->rght->tax == YES) is_ambigu = b->rght->c_seq->is_ambigu[site];
      if(is_ambigu == NO) observed_state = b->rght->c_seq->d_state[site];

      for(catg=0;catg<ncatg;++catg)
        {
          // Dot product left partial likelihood with equilibrium freqs
          for(i=0;i<nblocks;++i) _fplk[i] = _mm_mul_pd(_mm_load_pd(b->p_lk_left + site*nscatg + catg*ns + i*sz),
                                                       _mm_load_pd(tree->mod->e_frq->pi->v + i*sz));
          
          for(i=0;i<nblocks;++i) _colsum[i] = _mm_setzero_pd();
          // Multiply by matrix of right eigen vectors
          // This matrix is made of nblocks*nblocks squares,
          // each square being of 'area' sz*sz
          for(i=0;i<nblocks;++i) // row 
            {
              for(j=0;j<nblocks;++j) // column
                {
                  for(k=0;k<sz;++k) // column in sz*sz square
                    {
                      for(l=0;l<sz;++l) // row in that same square
                        {
                          // Copy right eigen vector values column-by-column in block (i,j)
                          ev[l] = tree->mod->eigen->r_e_vect[i*nblocks*sz*sz + j*sz + l*nblocks*sz + k];
                        }
                      _fplkev[k] = _mm_mul_pd(_mm_load_pd(ev),_fplk[i]);
                    }
                  _colsum[j] = _mm_add_pd(_colsum[j],_mm_hadd_pd(_fplkev[0],_fplkev[1])); // won't work if sz != 2
                }
            }          
          for(j=0;j<nblocks;++j) _mm_store_pd(tree->eigen_lr_left + catg*ns + j*sz,_colsum[j]);


          if(b->rght->tax == YES)
            for(i=0;i<nblocks;++i) _fplk[i] = _mm_load_pd(b->p_lk_tip_r + site*ns + i*sz); 
          else
            for(i=0;i<nblocks;++i) _fplk[i] = _mm_load_pd(b->p_lk_rght + site*nscatg + catg*ns + i*sz); 

          if(is_ambigu == YES)
            {
              for(i=0;i<nblocks;++i) _colsum[i] = _mm_setzero_pd();
              // Multiply by matrix of right eigen vectors
              // This matrix is made of nblocks*nblocks squares,
              // each square being of 'area' sz*sz
              for(i=0;i<nblocks;++i) // row
                {
                  for(j=0;j<nblocks;++j) // column
                    {
                      for(l=0;l<sz;++l) // row in that same square
                        {
                          for(k=0;k<sz;++k) // column in sz*sz square
                            {
                              // Copy left eigen vector values column-by-column in block (i,j)
                              ev[k] = tree->mod->eigen->l_e_vect[i*nblocks*sz*sz + j*sz + l*nblocks*sz + k];
                            }                          
                          _fplkev[l] = _mm_mul_pd(_mm_load_pd(ev),_fplk[j]);
                        }
                      _colsum[i] = _mm_add_pd(_colsum[i],_mm_hadd_pd(_fplkev[0],_fplkev[1]));
                    }
                }
              for(j=0;j<nblocks;++j) _mm_store_pd(tree->eigen_lr_rght + catg*ns + j*sz,_colsum[j]);
            }
          else
            {
              for(state=0;state<ns;++state)
                tree->eigen_lr_rght[catg*ns + state] =
                  tree->mod->eigen->l_e_vect[state*ns + observed_state];
            }
          for(j=0;j<nblocks;++j)
            _mm_store_pd(&(tree->dot_prod[site*nscatg + catg*ns + j*sz]),
                         _mm_mul_pd(_mm_load_pd(tree->eigen_lr_rght + catg*ns + j*sz),
                                    _mm_load_pd(tree->eigen_lr_left + catg*ns + j*sz)));
        }
      
    }
  
  if(ev)  Free(ev);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl SSE_Lk_Core_One_Class_Eigen_Lr(phydbl *dot_prod, phydbl *expl, int ns)
{
  phydbl lk;
  unsigned int l;
  const unsigned sz = (int)BYTE_ALIGN / 8;
  const unsigned int nblocks = ns/sz;
  __m128d _prod[nblocks],_x;
    
  for(l=0;l<nblocks;++l) _prod[l] = _mm_load_pd(dot_prod + l*sz);
  if(expl != NULL) for(l=0;l<nblocks;++l) _prod[l] = _mm_mul_pd(_prod[l],_mm_load_pd(expl + l*sz));
  _x = _mm_setzero_pd();
  for(l=0;l<nblocks;++l) _x = _mm_add_pd(_x,_prod[l]);
  _x = _mm_hadd_pd(_x,_x);
  _mm_store_sd(&lk,_x);
  
  return lk;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl SSE_Lk_Core_One_Class_No_Eigen_Lr(phydbl *p_lk_left, phydbl *p_lk_rght, phydbl *Pij, phydbl *pi, int ns, int ambiguity_check, int state)
{
  phydbl lk = 0.0;
  phydbl dum = 0.0;
  const unsigned int sz = (int)BYTE_ALIGN / 8;
  const unsigned nblocks = ns/sz;
  unsigned int i,j,k;
  __m128d _plk;
  __m128d _plk_l[nblocks],_plk_r[nblocks];
  
  /* [ Pi . Lkr ]' x Pij x Lkl */
  
  if(ambiguity_check == NO) // tip case.
    {
      for(i=0;i<nblocks;++i) _plk_l[i] = _mm_load_pd(p_lk_left + i*sz);      
      for(i=0;i<nblocks;++i) _plk_r[i] = _mm_load_pd(Pij + state*ns + i*sz);
      for(i=0;i<nblocks;++i) _plk_r[i] = _mm_mul_pd(_plk_r[i],_mm_set1_pd(pi[state]));
      for(i=0;i<nblocks;++i) _plk_r[i] = _mm_mul_pd(_plk_r[i],_plk_l[i]);
      _plk = _mm_setzero_pd();
      for(i=0;i<nblocks;++i) _plk = _mm_add_pd(_plk,_plk_r[i]);
      _plk = _mm_hadd_pd(_plk,_plk);
      _mm_store_sd(&lk,_plk);
    }
  else
    {
      __m128d _pij[sz]; 
      __m128d _pijplk[sz];

      for(i=0;i<nblocks;++i) _plk_r[i] = _mm_load_pd(p_lk_rght + i*sz);
      for(i=0;i<nblocks;++i) _plk_r[i] = _mm_mul_pd(_plk_r[i],_mm_load_pd(pi + i*sz));
      for(i=0;i<nblocks;++i) _plk_l[i] = _mm_load_pd(p_lk_left + i*sz);

      for(j=0;j<nblocks;++j)
        {
          for(i=0;i<sz;++i) _pijplk[i] = _mm_setzero_pd();

          for(i=0;i<nblocks;++i)
            {
              for(k=0;k<sz;++k)
                {
                  _pij[k] = _mm_load_pd(Pij + j*nblocks*sz*sz + i*sz + k*ns);
                  _pij[k] = _mm_mul_pd(_pij[k],_plk_l[i]);
                  _pijplk[k] = _mm_add_pd(_pijplk[k],_pij[k]);
                }
            }
          
          _plk = _mm_setzero_pd();
          for(k=0;k<sz;++k) _plk = _mm_hadd_pd(_plk,_pijplk[k]);
          
          _plk = _mm_mul_pd(_plk,_plk_r[j]);
          _mm_store_sd(&dum,_mm_hadd_pd(_plk,_plk));
          lk += dum;
        }
    }
  return lk;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void SSE_Update_Partial_Lk(t_tree *tree, t_edge *b, t_node *d)
{
/*
           |
           |<- b
           |
           d
          / \
         /   \
        /     \
    n_v1   n_v2
*/
  t_node *n_v1, *n_v2;//d's "left" and "right" neighbor nodes
  phydbl *p_lk,*p_lk_v1,*p_lk_v2;//Partial likelihood vector of node d, d's "left" neighbor, d's "right" neighbor. We fill *p_lk, and assume *p_lk_v1 and *p_lk_v2 are already filled.
  phydbl *Pij1,*Pij2;
  phydbl *tPij1,*tPij2;
  int *sum_scale, *sum_scale_v1, *sum_scale_v2;
  int sum_scale_v1_val, sum_scale_v2_val;
  unsigned int  i,k;
  unsigned int catg/*index over the number of rate categories*/;
  unsigned int site;
  phydbl smallest_p_lk,largest_p_lk;
  phydbl p_lk_lim_inf;
  int *p_lk_loc; //Suppose site j, of a certain subtree, has "A" on one tip, and "C" on the other. If you come across this pattern again at site i<j, then you can simply copy the partial likelihoods

  const unsigned int ns = tree->mod->ns;
  const unsigned int ncatg = tree->mod->ras->n_catg;
  const unsigned int npattern = tree->n_pattern;
  
  const unsigned int sz = (int)BYTE_ALIGN / 8;
  const unsigned nblocks = ns/sz;
  
  const unsigned int ncatgns =  ncatg * ns;
  const unsigned int nsns =  ns * ns;
  
  __m128d _x1[nblocks], _x2[nblocks];
  __m128d _plk[nblocks],_plk1[nblocks],_plk2[nblocks];

    
  if(tree->n_root && tree->ignore_root == YES &&
     (d == tree->n_root->v[1] || d == tree->n_root->v[2]) &&
     (b == tree->n_root->b[1] || b == tree->n_root->b[2]))
    {
      assert(FALSE);
    }

  sum_scale_v1_val = sum_scale_v2_val = 0;

  if(d->tax)
    {
      PhyML_Printf("\n== t_node %d is a leaf...",d->num);
      PhyML_Printf("\n== Err. in file %s at line %d (function '%s')\n",__FILE__,__LINE__,__FUNCTION__);
      Exit("\n");
    }

  sum_scale_v1_val            = 0;
  sum_scale_v2_val            = 0;  
  p_lk_lim_inf                = (phydbl)P_LK_LIM_INF;
  n_v1 = n_v2                 = NULL;
  p_lk = p_lk_v1 = p_lk_v2    = NULL;
  Pij1 = Pij2                 = NULL;
  tPij1 = tPij2               = NULL;
  sum_scale_v1 = sum_scale_v2 = NULL;
  p_lk_loc                    = NULL;

  Set_All_Partial_Lk(&n_v1,&n_v2,
                     &p_lk,&sum_scale,&p_lk_loc,
                     &Pij1,&tPij1,&p_lk_v1,&sum_scale_v1,
                     &Pij2,&tPij2,&p_lk_v2,&sum_scale_v2,
                     d,b,tree);

  if(tree->mod->augmented == YES && n_v1 && n_v1->tax == NO)
    {
      PhyML_Printf("\n== SSE version of the Update_Partial_Lk function does not");
      PhyML_Printf("\n== allow augmented data.");
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
    }

  /* For every site in the alignment */
  for(site=0;site<npattern;++site)
    {
      if(n_v1->tax == YES &&
         n_v2->tax == YES &&
         n_v1->c_seq->state[site] == 'X' &&
         n_v2->c_seq->state[site] == 'X')
        {
          for(i=0;i<ncatgns;++i) p_lk[site*ncatgns + i] = 1.0;
        }                                                  
      else
        {
          /* For all rate classes */
          for(catg=0;catg<ncatg;++catg)
            {                            
              for(k=0;k<nblocks;++k) _plk1[k] =  _mm_setzero_pd();
              for(k=0;k<nblocks;++k) _plk2[k] =  _mm_setzero_pd();
              
              for(i=0;i<ns;++i) // for each *column*
                {
                  
                  for(k=0;k<nblocks;++k) _x1[k] =  _mm_setzero_pd();
                  for(k=0;k<nblocks;++k) _x2[k] =  _mm_setzero_pd();
                  
                  if(n_v1->tax == NO)
                    {
                      if(p_lk_v1[site*ncatgns + catg*ns + i] > 0.0)
                        for(k=0;k<nblocks;++k) _x1[k] = _mm_mul_pd(_mm_load_pd(tPij1 + catg*nsns + i*ns + sz*k),
                                                                   _mm_set1_pd(p_lk_v1[site*ncatgns + catg*ns + i]));
                      
                      for(k=0;k<nblocks;++k) _plk1[k] = _mm_add_pd(_plk1[k],_x1[k]);
                    }
                  else
                    {
                      if(p_lk_v1[site*ns + i] > 0.0)
                        for(k=0;k<nblocks;++k) _x1[k] = _mm_mul_pd(_mm_load_pd(tPij1 + catg*nsns + i*ns + sz*k),
                                                                   _mm_set1_pd(p_lk_v1[site*ns + i]));
                      
                      for(k=0;k<nblocks;++k) _plk1[k] = _mm_add_pd(_plk1[k],_x1[k]);
                    }
                  
                  if(n_v2->tax == NO)
                    {
                      if(p_lk_v2[site*ncatgns + catg*ns + i] > 0.0)
                        for(k=0;k<nblocks;++k) _x2[k] = _mm_mul_pd(_mm_load_pd(tPij2 + catg*nsns + i*ns + sz*k),
                                                                   _mm_set1_pd(p_lk_v2[site*ncatgns + catg*ns + i]));
                      
                      for(k=0;k<nblocks;++k) _plk2[k] = _mm_add_pd(_plk2[k],_x2[k]);
                    }
                  else
                    {
                      if(p_lk_v2[site*ns + i] > 0.0)
                        for(k=0;k<nblocks;++k) _x2[k] = _mm_mul_pd(_mm_load_pd(tPij2 + catg*nsns + i*ns + sz*k),
                                                                   _mm_set1_pd(p_lk_v2[site*ns + i]));
                      
                      for(k=0;k<nblocks;++k) _plk2[k] = _mm_add_pd(_plk2[k],_x2[k]);
                    }                  
                }
              
              for(k=0;k<nblocks;++k) _plk[k] = _mm_mul_pd(_plk1[k],_plk2[k]);
              
              for(k=0;k<nblocks;++k) _mm_store_pd(p_lk + site*ncatgns + catg*ns + sz*k,_plk[k]);
              
              
              if(tree->scaling_method == SCALE_RATE_SPECIFIC)
                {                  
                  smallest_p_lk = BIG;
                  for(i=0;i<ns;++i) 
                    if(p_lk[site*ncatgns+catg*ns+i] < smallest_p_lk) 
                      smallest_p_lk = p_lk[site*ncatgns+catg*ns+i] ;
                  
                  /* Current scaling values at that site */
                  sum_scale_v1_val = (sum_scale_v1)?(sum_scale_v1[site*ncatg+catg]):(0);
                  sum_scale_v2_val = (sum_scale_v2)?(sum_scale_v2[site*ncatg+catg]):(0);
                  
                  sum_scale[site*ncatg+catg] = sum_scale_v1_val + sum_scale_v2_val;
                  
                  /* Scaling. We have p_lk_lim_inf = 2^-500. Consider for instance that
                     smallest_p_lk = 2^-600, then curr_scaler_pow will be equal to 100, and
                     each element in the partial likelihood vector will be multiplied by
                     2^100. */
                  if(smallest_p_lk < p_lk_lim_inf &&
                     tree->mod->augmented == NO &&
                     tree->apply_lk_scaling == YES)
                    {
                      int curr_scaler_pow;
                      curr_scaler_pow = (int)(-500.*LOG2-log(smallest_p_lk))/LOG2;
                      sum_scale[site*ncatg+catg] += curr_scaler_pow;
                      for(i=0;i<ns;++i) Rate_Correction(curr_scaler_pow, p_lk + site*ncatgns + catg*ns + i);
                    }
                }
            }              
        }
            
      if(tree->scaling_method == SCALE_FAST)
        {
          sum_scale_v1_val = (sum_scale_v1)?(sum_scale_v1[site]):(0);
          sum_scale_v2_val = (sum_scale_v2)?(sum_scale_v2[site]):(0);              
          sum_scale[site] = sum_scale_v1_val + sum_scale_v2_val;
          
          largest_p_lk = -BIG; 
          for(i=0;i<ns*ncatg;++i)
            if(p_lk[site*ncatgns+i] > largest_p_lk)
              largest_p_lk = p_lk[site*ncatgns+i] ;
          
          if(largest_p_lk < INV_TWO_TO_THE_LARGE &&
             tree->mod->augmented == NO &&
             tree->apply_lk_scaling == YES)
            {
              for(i=0;i<ns*ncatg;++i) p_lk[site*ncatgns + i] *= TWO_TO_THE_LARGE;
              sum_scale[site] += LARGE;
            }
        }
    }
}

#endif
