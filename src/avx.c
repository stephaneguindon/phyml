/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#if defined(__AVX__)

#include "assert.h"
#include "avx.h"

static void AVX_Partial_Lk_Exex(__m256d *_tPij1, int state1, __m256d *_tPij2, int state2, int ns, __m256d *plk0);
static void AVX_Partial_Lk_Exin(__m256d *_tPij1, int state1, __m256d *_tPij2, __m256d *_plk2, __m256d *_pmat2plk2, int ns, __m256d *_plk0);
static void AVX_Partial_Lk_Inin(__m256d *_tPij1, __m256d *_plk1, __m256d *_pmat1plk1, __m256d *_tPij2, __m256d *_plk2, __m256d *_pmat2plk2, int ns, __m256d *_plk0);
static void AVX_Matrix_Vect_Prod(__m256d *_m_transpose,  __m256d *_v, int ns, __m256d *res);

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void AVX_Update_Eigen_Lr(t_edge *b, t_tree *tree)
{
  unsigned int site,catg,state;
  unsigned int i,j,k,l;
  
  unsigned const int npattern = tree->n_pattern;
  unsigned const int ncatg = tree->mod->ras->n_catg;
  unsigned const int ns = tree->mod->ns;
  unsigned const int sz = (int)BYTE_ALIGN / 8;
  unsigned const int ncatgns = ncatg * ns;
  unsigned const int nblocks = ns / sz;

  phydbl *p_lk_left,*p_lk_rght,*pi,*dot_prod;
  phydbl *l_ev,*r_ev;

  __m256d *_l_ev,*_r_ev,*_p_lk_left,*_p_lk_rght,*_prod_left,*_prod_rght;
  
#ifndef WIN32
  if(posix_memalign((void **)&l_ev,BYTE_ALIGN,(size_t) ns * ns * sizeof(phydbl))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  if(posix_memalign((void **)&_l_ev,BYTE_ALIGN,(size_t) ns * ns / sz * sizeof(__m256d))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  if(posix_memalign((void **)&_r_ev,BYTE_ALIGN,(size_t) ns * ns / sz * sizeof(__m256d))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  if(posix_memalign((void **)&_p_lk_left,BYTE_ALIGN,(size_t) ns / sz * sizeof(__m256d))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  if(posix_memalign((void **)&_p_lk_rght,BYTE_ALIGN,(size_t) ns / sz * sizeof(__m256d))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  if(posix_memalign((void **)&_prod_left,BYTE_ALIGN,(size_t) ns / sz * sizeof(__m256d))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  if(posix_memalign((void **)&_prod_rght,BYTE_ALIGN,(size_t) ns / sz * sizeof(__m256d))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
#else
  l_ev       = _aligned_malloc(ns * ns * sizeof(phydbl),BYTE_ALIGN);
  _l_ev      = _aligned_malloc(ns * ns / sz * sizeof(__m256d),BYTE_ALIGN);
  _r_ev      = _aligned_malloc(ns * ns / sz * sizeof(__m256d),BYTE_ALIGN);
  _p_lk_left = _aligned_malloc(ns / sz * sizeof(__m256d),BYTE_ALIGN);
  _p_lk_rght = _aligned_malloc(ns / sz * sizeof(__m256d),BYTE_ALIGN);
  _prod_left = _aligned_malloc(ns / sz * sizeof(__m256d),BYTE_ALIGN);
  _prod_rght = _aligned_malloc(ns / sz * sizeof(__m256d),BYTE_ALIGN);
#endif
  

  assert(sz == 4);
  assert(tree->update_eigen_lr == YES);
  
  r_ev = tree->mod->eigen->r_e_vect;
  
  /* Copy transpose of matrix of left eigen vectors */
  for(i=0;i<ns;++i)
    for(j=0;j<ns;++j)
      l_ev[i*ns+j] = tree->mod->eigen->l_e_vect[j*ns+i];
  
  /* Load into AVX registers */
  for(i=0;i<ns;++i)
    {
      for(j=0;j<nblocks;++j)
        {
          _r_ev[i*nblocks+j] = _mm256_load_pd(r_ev + j*sz);
          _l_ev[i*nblocks+j] = _mm256_load_pd(l_ev + j*sz);
        }
      r_ev += ns;
      l_ev += ns;
    }
          
  l_ev -= ns*ns;
  r_ev -= ns*ns;
  
  p_lk_left = b->left->tax ? b->p_lk_tip_l : b->p_lk_left;
  p_lk_rght = b->rght->tax ? b->p_lk_tip_r : b->p_lk_rght;
  pi = tree->mod->e_frq->pi->v;
  dot_prod = tree->dot_prod;
  
  for(site=0;site<npattern;++site)
    {
      for(catg=0;catg<ncatg;++catg)
        {
          for(i=0;i<nblocks;++i)
            {
              // Dot product left partial likelihood with equilibrium freqs and load...
              _p_lk_left[i] = _mm256_mul_pd(_mm256_load_pd(p_lk_left + i*sz),
                                            _mm256_load_pd(pi + i*sz));
              // ...just load.
              _p_lk_rght[i] = _mm256_load_pd(p_lk_rght + i*sz);              
            }

          AVX_Matrix_Vect_Prod(_r_ev,_p_lk_left,ns,_prod_left);
          AVX_Matrix_Vect_Prod(_l_ev,_p_lk_rght,ns,_prod_rght);
          
          for(i=0;i<nblocks;++i) _mm256_store_pd(dot_prod + i*sz,_mm256_mul_pd(_prod_left[i],_prod_rght[i]));

          dot_prod += ns;
          if(b->left->tax == NO) p_lk_left += ns;
          if(b->rght->tax == NO) p_lk_rght += ns;
        }

      if(b->left->tax == YES) p_lk_left += ns;
      if(b->rght->tax == YES) p_lk_rght += ns;
    }

  Free(l_ev);
  Free(_l_ev);
  Free(_r_ev);
  Free(_p_lk_left);
  Free(_p_lk_rght);
  Free(_prod_left);
  Free(_prod_rght);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl AVX_Lk_Core_One_Class_No_Eigen_Lr(phydbl *p_lk_left, phydbl *p_lk_rght, phydbl *Pij, phydbl *pi, int ns, int ambiguity_check, int observed_state)
{
  const unsigned int sz = (int)BYTE_ALIGN / 8;
  const unsigned nblocks = ns/sz;

  __m256d _plk_l[nblocks],_plk_r[nblocks];
  __m256d _plk; // parent partial likelihood

  /* [ Pi . Lkr ]' x Pij x Lkl */

  if(ambiguity_check == NO) // tip case.
    {
      unsigned int i;
      for(i=0;i<nblocks;++i) _plk_l[i] = _mm256_load_pd(p_lk_left + i*sz);
      for(i=0;i<nblocks;++i) _plk_r[i] = _mm256_load_pd(Pij + observed_state*ns + i*sz);
      for(i=0;i<nblocks;++i) _plk_r[i] = _mm256_mul_pd(_plk_r[i],_mm256_set1_pd(pi[observed_state]));
      for(i=0;i<nblocks;++i) _plk_r[i] = _mm256_mul_pd(_plk_r[i],_plk_l[i]);
      _plk = _mm256_setzero_pd();
      for(i=0;i<nblocks;++i) _plk = _mm256_add_pd(_plk,_plk_r[i]);
      return AVX_Vect_Norm(_plk);
    }
  else
    {
      unsigned int i,j,k;
      __m256d _pij[sz]; 
      __m256d _pijplk[sz];
      phydbl lk=.0;

      for(i=0;i<nblocks;++i) _plk_r[i] = _mm256_load_pd(p_lk_rght + i*sz);
      for(i=0;i<nblocks;++i) _plk_r[i] = _mm256_mul_pd(_plk_r[i],_mm256_load_pd(pi + i*sz));
      for(i=0;i<nblocks;++i) _plk_l[i] = _mm256_load_pd(p_lk_left + i*sz);

      for(j=0;j<nblocks;++j)
        {
          for(i=0;i<sz;++i) _pijplk[i] = _mm256_setzero_pd();

          for(i=0;i<nblocks;++i)
            {
              for(k=0;k<sz;++k)
                {
                  _pij[k] = _mm256_load_pd(Pij + j*nblocks*sz*sz + i*sz + k*ns);
                  _pij[k] = _mm256_mul_pd(_pij[k],_plk_l[i]);
                  _pijplk[k] = _mm256_add_pd(_pijplk[k],_pij[k]);
                }
            }
          
          _plk = AVX_Horizontal_Add(_pijplk);
          _plk = _mm256_mul_pd(_plk,_plk_r[j]);

          lk += AVX_Vect_Norm(_plk);
        }
      return lk;
    }
  return UNLIKELY;
}
 
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl AVX_Lk_Core_One_Class_Eigen_Lr(phydbl *dot_prod, phydbl *expl, unsigned int ns)
{
  unsigned int l;
  const unsigned int sz = (int)BYTE_ALIGN / 8;
  const unsigned nblocks = ns/sz;
  __m256d _prod[nblocks],_x;
  
  for(l=0;l<nblocks;++l) _prod[l] = _mm256_load_pd(dot_prod + l*sz);
  if(expl != NULL) for(l=0;l<nblocks;++l) _prod[l] = _mm256_mul_pd(_prod[l],_mm256_load_pd(expl + l*sz));
  _x = _mm256_setzero_pd();
  for(l=0;l<nblocks;++l) _x = _mm256_add_pd(_x,_prod[l]);
  
  return AVX_Vect_Norm(_x);
}
 
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void AVX_Lk_dLk_Core_One_Class_Eigen_Lr(phydbl *dot_prod, phydbl *expl, unsigned int ns, phydbl *lk, phydbl *dlk)
{
  unsigned int i;
  const unsigned int sz = (int)BYTE_ALIGN / 8;
  const unsigned nblocks = ns/sz*2;
  __m256d _x,_y,_z;

  _z = _mm256_setzero_pd();

  for(i=0;i<nblocks;++i)
    {
      _x = _mm256_blend_pd(_mm256_set1_pd(dot_prod[2*i]),
                           _mm256_set1_pd(dot_prod[2*i + 1]),12);
      
      _y = _mm256_load_pd(expl + 4*i);

      _z = _mm256_add_pd(_z,_mm256_mul_pd(_x,_y));
    }
  
  *lk = ((double *)&_z)[0] + ((double *)&_z)[2];
  *dlk = ((double *)&_z)[1] + ((double *)&_z)[3];

}
 
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl AVX_Vect_Norm(__m256d _z)
{
  phydbl r;
  __m256d _x = _mm256_hadd_pd(_z,_z);
  __m256d _y = _mm256_permute2f128_pd(_x,_x,0x21);
  _mm_store_sd(&r,_mm256_castpd256_pd128(_mm256_add_pd(_x,_y)));
  return r;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

__m256d AVX_Horizontal_Add(__m256d x[4])
{
  __m256d y[2],z[2];

  // y[0] = [x00+x01;x10+x11;x02+x03;x12+x13]
  y[0] = _mm256_hadd_pd(x[0], x[1]);
  // y[1] = [x20+x21;x30+x31;x22+x23;x32+x33]
  y[1] = _mm256_hadd_pd(x[2], x[3]);

  // z[0] = [x00+x01;x10+x11;x22+x23;x32+x33]
  /* z[0] = _mm256_blend_pd(y[0],y[1],0b1100); */
  z[0] = _mm256_blend_pd(y[0],y[1],12);
  // z[1] = [x02+x03;x12+x13;x20+x21;x30+x31]
  z[1] = _mm256_permute2f128_pd(y[0],y[1],0x21);

  return(_mm256_add_pd(z[0],z[1]));
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void AVX_Update_Partial_Lk(t_tree *tree, t_edge *b, t_node *d)
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
  t_node *n_v1, *n_v2;
  phydbl *plk0,*plk1,*plk2;
  phydbl *Pij1,*Pij2;
  phydbl *tPij1,*tPij2;
  int *sum_scale, *sum_scale_v1, *sum_scale_v2;
  int sum_scale_v1_val, sum_scale_v2_val;
  unsigned int i,j,k;
  unsigned int catg,site;
  unsigned int v1_idx,v2_idx,p_idx,v0_idx;
  short int state_v1,state_v2;
  short int ambiguity_check_v1,ambiguity_check_v2;
  phydbl smallest_p_lk,largest_p_lk;
  phydbl p_lk_lim_inf;
  int *p_lk_loc;
  short int p1_is_fully_ambiguous, p2_is_fully_ambiguous;
  
  const unsigned int npattern = tree->n_pattern;
  const unsigned int ns = tree->mod->ns;
  const unsigned int ncatg = tree->mod->ras->n_catg;

  const unsigned int ncatgns =  ncatg * ns;
  const unsigned int nsns =  ns * ns;
  
  const unsigned int sz = (int)BYTE_ALIGN / 8;
  const unsigned nblocks = ns/sz;

  __m256d *_tPij1,*_tPij2,*_pmat1plk1,*_pmat2plk2,*_plk1,*_plk2,*_plk0;
  
#ifndef WIN32
  if(posix_memalign((void **)&_tPij1,BYTE_ALIGN,(size_t)(ncatg * nsns / sz) * sizeof(__m256d))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  if(posix_memalign((void **)&_tPij2,BYTE_ALIGN,(size_t)(ncatg * nsns / sz) * sizeof(__m256d))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  if(posix_memalign((void **)&_pmat1plk1,BYTE_ALIGN,(size_t)ns / sz * sizeof(__m256d))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  if(posix_memalign((void **)&_pmat2plk2,BYTE_ALIGN,(size_t)ns / sz * sizeof(__m256d))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  if(posix_memalign((void **)&_plk0,BYTE_ALIGN,(size_t)(ns / sz) * sizeof(__m256d))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  if(posix_memalign((void **)&_plk1,BYTE_ALIGN,(size_t)(ns / sz) * sizeof(__m256d))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  if(posix_memalign((void **)&_plk2,BYTE_ALIGN,(size_t)(ns / sz) * sizeof(__m256d))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
#else
  _tPij1     = _aligned_malloc(ncatg * nsns / sz * sizeof(__m256d),BYTE_ALIGN);
  _tPij2     = _aligned_malloc(ncatg * nsns / sz * sizeof(__m256d),BYTE_ALIGN);
  _pmat1plk1 = _aligned_malloc(ns / sz * sizeof(__m256d),BYTE_ALIGN);
  _pmat2plk2 = _aligned_malloc(ns / sz * sizeof(__m256d),BYTE_ALIGN);
  _plk0      = _aligned_malloc(ns / sz * sizeof(__m256d),BYTE_ALIGN);
  _plk1      = _aligned_malloc(ns / sz * sizeof(__m256d),BYTE_ALIGN);
  _plk2      = _aligned_malloc(ns / sz * sizeof(__m256d),BYTE_ALIGN);
#endif

  
  
  sum_scale_v1_val            = 0;
  sum_scale_v2_val            = 0;
  p_lk_lim_inf                = (phydbl)P_LK_LIM_INF;
  n_v1 = n_v2                 = NULL;
  plk0 = plk1 = plk2          = NULL;
  Pij1 = Pij2                 = NULL;
  tPij1 = tPij2               = NULL;
  sum_scale_v1 = sum_scale_v2 = NULL;
  p_lk_loc                    = NULL;
  p1_is_fully_ambiguous       = NO;
  p2_is_fully_ambiguous       = NO;
  state_v1 = state_v2         = -1;

  
  Set_All_Partial_Lk(&n_v1,&n_v2,
                     &plk0,&sum_scale,&p_lk_loc,
                     &Pij1,&tPij1,&plk1,&sum_scale_v1,
                     &Pij2,&tPij2,&plk2,&sum_scale_v2,
                     d,b,tree);
 
  // Copy transpose of transition matrices into AVX registers
  for(i=0;i<ncatg;++i)
    {
      for(j=0;j<ns;++j)
        {
          for(k=0;k<nblocks;++k)
            {
              _tPij1[k] = _mm256_load_pd(tPij1);
              _tPij2[k] = _mm256_load_pd(tPij2);
              tPij1 += sz;
              tPij2 += sz;
            }
          _tPij1 += nblocks;
          _tPij2 += nblocks;
        }
    }
  _tPij1 = _tPij1 - (int)(ncatg * nsns / sz);
  _tPij2 = _tPij2 - (int)(ncatg * nsns / sz);
  tPij1 -= ncatg * nsns;
  tPij2 -= ncatg * nsns;

  if(tree->mod->augmented == YES)
    {
      PhyML_Printf("\n== AVX version of the Update_Partial_Lk function does not");
      PhyML_Printf("\n== allow augmented data.");
      assert(FALSE);
    }
    
  /* For every site in the alignment */
  for(site=0;site<npattern;++site)
    {
      v0_idx = site * ncatgns;

      if(n_v1->tax == YES &&
         n_v2->tax == YES &&
         n_v1->c_seq->state[site] == 'X' &&
         n_v2->c_seq->state[site] == 'X')
        {
          for(i=0;i<ncatgns;++i) plk0[v0_idx + i] = 1.0;
        }
      else
        {
          v1_idx = n_v1->tax ? site*ns : site*ncatgns;
          v2_idx = n_v2->tax ? site*ns : site*ncatgns;

          p1_is_fully_ambiguous = NO;
          for(i=0;i<ns;++i) if(Are_Equal(plk1[v1_idx+i],1.0,1.E-5) == NO) break;
          if(i==ns) p1_is_fully_ambiguous = YES;

          p2_is_fully_ambiguous = NO;
          for(i=0;i<ns;++i) if(Are_Equal(plk2[v2_idx+i],1.0,1.E-5) == NO) break;
          if(i==ns) p2_is_fully_ambiguous = YES;
          
          if(p1_is_fully_ambiguous == YES && p2_is_fully_ambiguous == YES)
            {
              for(i=0;i<ncatgns;++i) plk0[v0_idx + i] = 1.0;
            }
          else
            {
              state_v1 = state_v2 = -1;
              ambiguity_check_v1 = ambiguity_check_v2 = YES;
              
              if(n_v1 && n_v1->tax)
                {
                  ambiguity_check_v1 = n_v1->c_seq->is_ambigu[site];
                  if(ambiguity_check_v1 == NO) state_v1 = n_v1->c_seq->d_state[site];
                }
              
              if(n_v2 && n_v2->tax)
                {
                  ambiguity_check_v2 = n_v2->c_seq->is_ambigu[site];
                  if(ambiguity_check_v2 == NO) state_v2 = n_v2->c_seq->d_state[site];
                }
              
              for(catg=0;catg<ncatg;++catg)
                {
                  v1_idx = n_v1->tax ? site*ns : site*ncatgns+catg*ns;
                  v2_idx = n_v2->tax ? site*ns : site*ncatgns+catg*ns;
                  p_idx  = (int)(nsns*catg / sz);
                                        
                  for(k=0;k<nblocks;++k)
                    {
                      _plk1[k] = _mm256_load_pd(plk1 + v1_idx + k*sz);
                      _plk2[k] = _mm256_load_pd(plk2 + v2_idx + k*sz);
                    }
                  
                  if(ambiguity_check_v1 == NO && ambiguity_check_v2 == NO)
                    {
                      AVX_Partial_Lk_Exex(_tPij1 + p_idx,state_v1,
                                          _tPij2 + p_idx,state_v2,
                                          ns,_plk0);
                    }
                  else if(ambiguity_check_v1 == YES && ambiguity_check_v2 == NO)
                    {
                      AVX_Partial_Lk_Exin(_tPij2 + p_idx,state_v2,
                                          _tPij1 + p_idx,_plk1,_pmat1plk1,
                                          ns,_plk0);
                    }
                  else if(ambiguity_check_v1 == NO && ambiguity_check_v2 == YES)
                    {
                      AVX_Partial_Lk_Exin(_tPij1 + p_idx,state_v1,
                                          _tPij2 + p_idx,_plk2,_pmat2plk2,
                                          ns,_plk0);
                    }
                  else
                    {
                      AVX_Partial_Lk_Inin(_tPij1 + p_idx,_plk1,_pmat1plk1,
                                          _tPij2 + p_idx,_plk2,_pmat2plk2,
                                          ns,_plk0);
                    }
                                    
                  for(k=0;k<nblocks;++k) _mm256_store_pd(plk0+v0_idx+catg*ns+sz*k,_plk0[k]);
              

                  if(tree->scaling_method == SCALE_RATE_SPECIFIC)
                    {
                      smallest_p_lk = BIG;
                      for(i=0;i<ns;++i)
                        if(plk0[v0_idx+catg*ns+i] < smallest_p_lk)
                          smallest_p_lk = plk0[v0_idx+catg*ns+i];
                      
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
                          for(i=0;i<ns;++i) Rate_Correction(curr_scaler_pow,plk0+v0_idx+catg*ns+i);
                        }
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
          for(i=0;i<ncatgns;++i)
            if(plk0[v0_idx+i] > largest_p_lk)
              largest_p_lk = plk0[v0_idx+i];
          
          if(largest_p_lk < INV_TWO_TO_THE_LARGE &&
             tree->mod->augmented == NO &&
             tree->apply_lk_scaling == YES)
            {
              for(i=0;i<ncatgns;++i) plk0[v0_idx+i] *= TWO_TO_THE_LARGE;
              sum_scale[site] += LARGE;
            }
        }
    }

  Free(_tPij1);
  Free(_tPij2);
  Free(_pmat1plk1);
  Free(_pmat2plk2);
  Free(_plk0);
  Free(_plk1);
  Free(_plk2);
  
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

static void AVX_Partial_Lk_Exex(__m256d *_tPij1, int state1, __m256d *_tPij2, int state2, int ns, __m256d *plk0)
{
  unsigned const int sz = (int)BYTE_ALIGN / 8;
  unsigned const int nblocks = ns / sz;
  unsigned int i;

  _tPij1 = _tPij1 + state1 * nblocks;
  _tPij2 = _tPij2 + state2 * nblocks;
  for(i=0;i<nblocks;++i) plk0[i] = _mm256_mul_pd(_tPij1[i],_tPij2[i]);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

static void AVX_Partial_Lk_Exin(__m256d *_tPij1, int state1, __m256d *_tPij2, __m256d *_plk2, __m256d *_pmat2plk2, int ns, __m256d *_plk0)
{
  unsigned const int sz = (int)BYTE_ALIGN / 8;
  unsigned const int nblocks = ns / sz;
  unsigned int i;
  
  _tPij1 = _tPij1 + state1 * nblocks;
  AVX_Matrix_Vect_Prod(_tPij2,_plk2,ns,_pmat2plk2);
  
  for(i=0;i<nblocks;++i) _plk0[i] = _mm256_mul_pd(_tPij1[i],_pmat2plk2[i]);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

static void AVX_Partial_Lk_Inin(__m256d *_tPij1, __m256d *_plk1, __m256d *_pmat1plk1, __m256d *_tPij2, __m256d *_plk2, __m256d *_pmat2plk2, int ns, __m256d *_plk0)
{
  unsigned const int sz = (int)BYTE_ALIGN / 8;
  unsigned const int nblocks = ns / sz;
  unsigned int i;

  AVX_Matrix_Vect_Prod(_tPij1,_plk1,ns,_pmat1plk1);
  AVX_Matrix_Vect_Prod(_tPij2,_plk2,ns,_pmat2plk2);
  
  for(i=0;i<nblocks;++i) _plk0[i] = _mm256_mul_pd(_pmat1plk1[i],_pmat2plk2[i]);
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

static void AVX_Matrix_Vect_Prod(__m256d *_m_transpose,  __m256d *_v, int ns, __m256d *_u)
{
  unsigned const int sz = (int)BYTE_ALIGN / 8;
  unsigned const int nblocks = ns / sz;
  unsigned int i,j;  
  __m256d _x,_y;
    
  for(j=0;j<nblocks;++j) _u[j] = _mm256_setzero_pd();

  for(i=0;i<ns;++i) // for each column
    {
      for(j=0;j<nblocks;++j) // for each block of sz rows
        {
          _x = _mm256_set1_pd(((double *)&_v[(int)(i/sz)])[i]);
          _y = _mm256_mul_pd(*_m_transpose,_x);
          _u[j] = _mm256_add_pd(_u[j],_y);
          _m_transpose++;
        }
    }
}

#endif

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
