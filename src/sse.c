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

#if ((defined(__SSE__) || defined(__SSE2__) || defined(__SSE3__)) && !((defined __AVX__ || defined __AVX2__)))

void SSE_Update_Eigen_Lr(t_edge *b, t_tree *tree)
{
  unsigned int site,catg;
  unsigned int i,j;
  
  unsigned const int npattern = tree->n_pattern;
  unsigned const int ncatg = tree->mod->ras->n_catg;
  unsigned const int ns = tree->mod->ns;
  unsigned const int sz = (int)BYTE_ALIGN / 8;
  unsigned const int nblocks = ns / sz;
  unsigned const int ncatgns = ncatg*ns;

  const phydbl *p_lk_left,*p_lk_rght,*pi;
  phydbl *dot_prod,*p_lk_left_pi;
  phydbl *l_ev,*r_ev;

  __m128d *_l_ev,*_r_ev,*_prod_left,*_prod_rght;
  
  p_lk_left_pi = tree->p_lk_left_pi;
  l_ev         = tree->l_ev;
  _l_ev        = tree->_l_ev;
  _r_ev        = tree->_r_ev;
  _prod_left   = tree->_prod_left;
  _prod_rght   = tree->_prod_rght;
  
  assert(sz == 2);
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
          _r_ev[i*nblocks+j] = _mm_load_pd(r_ev + j*sz);
          _l_ev[i*nblocks+j] = _mm_load_pd(l_ev + j*sz);
        }
      r_ev += ns;
      l_ev += ns;
    }

  p_lk_left = b->left->tax ? b->p_lk_tip_l : b->p_lk_left;
  p_lk_rght = b->rght->tax ? b->p_lk_tip_r : b->p_lk_rght;
  pi = tree->mod->e_frq->pi->v;
  dot_prod = tree->dot_prod;
  
  for(site=0;site<npattern;++site)
    {
      if(tree->data->wght[site] > SMALL)
        {
          for(catg=0;catg<ncatg;++catg)
            {
              for(i=0;i<ns;++i) p_lk_left_pi[i] = p_lk_left[i] * pi[i];
              
              SSE_Matrix_Vect_Prod(_r_ev,p_lk_left_pi,ns,_prod_left);
              SSE_Matrix_Vect_Prod(_l_ev,p_lk_rght,ns,_prod_rght);
              
              for(i=0;i<nblocks;++i) _mm_store_pd(dot_prod + i*sz,_mm_mul_pd(_prod_left[i],_prod_rght[i]));
              
              dot_prod += ns;
              if(b->left->tax == NO) p_lk_left += ns;
              if(b->rght->tax == NO) p_lk_rght += ns;
            }
          
          if(b->left->tax == YES) p_lk_left += ns;
          if(b->rght->tax == YES) p_lk_rght += ns;
        }
      else
        {
          if(b->left->tax == YES) p_lk_left += ns;
          else p_lk_left += ncatgns;

          if(b->rght->tax == YES) p_lk_rght += ns;          
          else p_lk_rght += ncatgns;          

          dot_prod += ncatgns;          
        }
    }
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

void SSE_Lk_dLk_Core_One_Class_Eigen_Lr(phydbl *dot_prod, phydbl *expl, unsigned int ns, phydbl *lk, phydbl *dlk)
{
  unsigned int i;
  __m128d _x,_y,_z;

  _z = _mm_setzero_pd();

  for(i=0;i<ns;++i)
    {
      _x = _mm_set1_pd(dot_prod[i]);      
      _y = _mm_load_pd(expl + 2*i);

      _z = _mm_add_pd(_z,_mm_mul_pd(_x,_y));
    }
  
  *lk = ((double *)&_z)[0];
  *dlk = ((double *)&_z)[1];
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl SSE_Lk_Core_One_Class_No_Eigen_Lr(phydbl *p_lk_left, phydbl *p_lk_rght, phydbl *Pij, phydbl *tPij, phydbl *pi, int ns, int ambiguity_check, int state)
{
  phydbl lk,dum;
  const unsigned int sz = (int)BYTE_ALIGN / 8;
  const unsigned nblocks = ns/sz;
  unsigned int i,j;
  __m128d _plk;
  __m128d _plk_l[nblocks],_plk_r[nblocks];
  
  /* [ Pi . Lkr ]' x Pij x Lkl */
  
  if(ambiguity_check == NO) // tip case.
    {
      for(i=0;i<nblocks;++i)
        {
          _plk_l[i] = _mm_load_pd(p_lk_left);      
          _plk_r[i] = _mm_load_pd(Pij + state*ns);
          p_lk_left += sz;
          Pij += sz;
        }
      
      for(i=0;i<nblocks;++i)
        {
          _plk_r[i] = _mm_mul_pd(_plk_r[i],_mm_set1_pd(pi[state]));
          _plk_r[i] = _mm_mul_pd(_plk_r[i],_plk_l[i]);
        }
      
      _plk = _mm_setzero_pd();
      lk = 0.0;
      for(i=0;i<nblocks;++i)
        {
          _plk = _mm_hadd_pd(_plk_r[i],_plk_r[i]);
          _mm_store_sd(&dum,_plk);
          lk += dum;
        }
      return lk;
    }
  else
    {
      __m128d _pij[nblocks],_pijplk[nblocks]; 

      for(i=0;i<nblocks;++i)
        {
          _plk_r[i] = _mm_mul_pd(_mm_load_pd(p_lk_rght),_mm_load_pd(pi));
          p_lk_rght += sz;
          pi += sz;
        }

      for(i=0;i<nblocks;++i) _pijplk[i] = _mm_setzero_pd();

      for(i=0;i<ns;++i)
        {
          for(j=0;j<nblocks;++j)
            {
              _pij[j] = _mm_load_pd(tPij);
              tPij += sz;
              
              _pijplk[j] = _mm_add_pd(_pijplk[j],_mm_mul_pd(_pij[j],_mm_set1_pd(p_lk_left[i])));
            }
        }

      lk = 0.0;
      for(i=0;i<nblocks;++i)
        {
          _plk = _mm_mul_pd(_pijplk[i],_plk_r[i]);
          _plk = _mm_hadd_pd(_plk,_plk);
          _mm_store_sd(&dum,_plk);
          lk += dum;
        }
      return lk;
    }
  return UNLIKELY;
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
  t_node *n_v1, *n_v2;
  phydbl *plk0,*plk1,*plk2;
  phydbl *Pij1,*Pij2;
  phydbl *tPij1,*tPij2;
  int *sum_scale, *sum_scale_v1, *sum_scale_v2;
  int sum_scale_v1_val, sum_scale_v2_val;
  unsigned int i,j,k;
  unsigned int catg,site;
  short int state_v1,state_v2;
  short int ambiguity_check_v1,ambiguity_check_v2;
  phydbl largest_p_lk;
  int *p_lk_loc;
  
  const unsigned int npattern = tree->n_pattern;
  const unsigned int ns = tree->mod->ns;
  const unsigned int ncatg = tree->mod->ras->n_catg;

  const unsigned int ncatgns =  ncatg * ns;
  const unsigned int nsns =  ns * ns;
  
  const unsigned int sz = (int)BYTE_ALIGN / 8;
  const unsigned nblocks = ns/sz;

  __m128d *_tPij1,*_tPij2,*_pmat1plk1,*_pmat2plk2,*_plk0;
  
  _tPij1     = tree->_tPij1;
  _tPij2     = tree->_tPij2;
  _pmat1plk1 = tree->_pmat1plk1;
  _pmat2plk2 = tree->_pmat2plk2;
  _plk0      = tree->_plk0;

  sum_scale_v1_val            = 0;
  sum_scale_v2_val            = 0;
  n_v1 = n_v2                 = NULL;
  plk0 = plk1 = plk2          = NULL;
  Pij1 = Pij2                 = NULL;
  tPij1 = tPij2               = NULL;
  sum_scale_v1 = sum_scale_v2 = NULL;
  p_lk_loc                    = NULL;
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
              _tPij1[k] = _mm_load_pd(tPij1);
              _tPij2[k] = _mm_load_pd(tPij2);
              tPij1 += sz;
              tPij2 += sz;
            }
          _tPij1 += nblocks;
          _tPij2 += nblocks;
        }
    }
  _tPij1 -= ncatg*ns*nblocks;
  _tPij2 -= ncatg*ns*nblocks;

  if(tree->mod->augmented == YES)
    {
      PhyML_Printf("\n== AVX version of the Update_Partial_Lk function does not");
      PhyML_Printf("\n== allow augmented data.");
      assert(FALSE);
    }
    
  /* For every site in the alignment */
  for(site=0;site<npattern;++site)
    {
      if(tree->data->wght[site] > SMALL)
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
              if(ambiguity_check_v1 == NO && ambiguity_check_v2 == NO)
                {
                  SSE_Partial_Lk_Exex(_tPij1,state_v1,
                                      _tPij2,state_v2,
                                      ns,_plk0);
                }
              else if(ambiguity_check_v1 == YES && ambiguity_check_v2 == NO)
                {
                  SSE_Partial_Lk_Exin(_tPij2,state_v2,
                                      _tPij1,plk1,_pmat1plk1,
                                      ns,_plk0);
                }
              else if(ambiguity_check_v1 == NO && ambiguity_check_v2 == YES)
                {
                  SSE_Partial_Lk_Exin(_tPij1,state_v1,
                                      _tPij2,plk2,_pmat2plk2,
                                      ns,_plk0);
                }
              else
                {
                  SSE_Partial_Lk_Inin(_tPij1,plk1,_pmat1plk1,
                                      _tPij2,plk2,_pmat2plk2,
                                      ns,_plk0);
                }
              
              for(k=0;k<nblocks;++k) _mm_store_pd(plk0+sz*k,_plk0[k]);
              
              _tPij1 += nsns / sz;
              _tPij2 += nsns / sz;
              plk0 += ns;
              plk1 += (n_v1->tax) ? 0 : ns;
              plk2 += (n_v2->tax) ? 0 : ns;
            }
          
          _tPij1 -= ncatg * nsns / sz;
          _tPij2 -= ncatg * nsns / sz;

          plk1 += (n_v1->tax) ? ns : 0;
          plk2 += (n_v2->tax) ? ns : 0;
          
          if(tree->scaling_method == SCALE_FAST)
            {
              sum_scale_v1_val = (sum_scale_v1)?(sum_scale_v1[site]):(0);
              sum_scale_v2_val = (sum_scale_v2)?(sum_scale_v2[site]):(0);
              sum_scale[site] = sum_scale_v1_val + sum_scale_v2_val;
              
              if(sum_scale[site] >= 1024)
                {
                  /* plk0 -= ncatgns; */
                  /* plk1 -= (n_v1->tax) ? ns : ncatgns; */
                  /* plk2 -= (n_v2->tax) ? ns : ncatgns; */
                  /* PhyML_Fprintf(stderr,"\n. PARTIAL site: %d plk0: %p [%g %g %g %g] plk1: %p [%g %g %g %g] plk2: %p [%g %g %g %g]", */
                  /*               site, */
                  /*               plk0, */
                  /*               plk0[0], */
                  /*               plk0[1], */
                  /*               plk0[2], */
                  /*               plk0[3], */
                  /*               plk1, */
                  /*               plk1[0], */
                  /*               plk1[1], */
                  /*               plk1[2], */
                  /*               plk1[3], */
                  /*               plk2, */
                  /*               plk2[0], */
                  /*               plk2[1], */
                  /*               plk2[2], */
                  /*               plk2[3] */
                  /*               ); */
                  /* PhyML_Fprintf(stderr,"\n. PARTIAL site: %d d: %d n_v1: %d n_v2: %d",site,d->num,n_v1->num,n_v2->num); */
                  /* PhyML_Fprintf(stderr,"\n. PARTIAL site: %d sum n: %d sum n_v1: %d sum n_v2: %d",site,sum_scale[site],sum_scale_v1_val,sum_scale_v2_val); */
                  
                  /* plk0 += ncatgns; */
                  /* plk1 += (n_v1->tax) ? ns : ncatgns; */
                  /* plk2 += (n_v2->tax) ? ns : ncatgns; */
                  /* Exit("\n"); */
                }
              
              plk0 -= ncatgns;
              
              largest_p_lk = -BIG;
              for(i=0;i<ncatgns;++i)
                if(plk0[i] > largest_p_lk)
                  largest_p_lk = plk0[i];
              
              if(largest_p_lk < INV_TWO_TO_THE_LARGE &&
                 tree->mod->augmented == NO &&
                 tree->apply_lk_scaling == YES)
                {
                  for(i=0;i<ncatgns;++i) plk0[i] *= TWO_TO_THE_LARGE;
                  sum_scale[site] += LARGE;
                }
              
              plk0 += ncatgns;
            }
        }
      else
        {
          plk0 += ncatgns;
          plk1 += (n_v1->tax) ? ns : ncatgns;
          plk2 += (n_v2->tax) ? ns : ncatgns;        
        }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void SSE_Partial_Lk_Exex(const __m128d *_tPij1, const int state1, const __m128d *_tPij2, const int state2, const int ns, __m128d *plk0)
{
  unsigned const int sz = (int)BYTE_ALIGN / 8;
  unsigned const int nblocks = ns / sz;
  unsigned int i;

  _tPij1 = _tPij1 + state1 * nblocks;
  _tPij2 = _tPij2 + state2 * nblocks;
  for(i=0;i<nblocks;++i) plk0[i] = _mm_mul_pd(_tPij1[i],_tPij2[i]);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void SSE_Partial_Lk_Exin(const __m128d *_tPij1, const int state1, const __m128d *_tPij2, const phydbl *_plk2, __m128d *_pmat2plk2, const int ns, __m128d *_plk0)
{
  unsigned const int sz = (int)BYTE_ALIGN / 8;
  unsigned const int nblocks = ns / sz;
  unsigned int i;
  
  _tPij1 = _tPij1 + state1 * nblocks;
  SSE_Matrix_Vect_Prod(_tPij2,_plk2,ns,_pmat2plk2);
  
  for(i=0;i<nblocks;++i) _plk0[i] = _mm_mul_pd(_tPij1[i],_pmat2plk2[i]);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void SSE_Partial_Lk_Inin(const __m128d *_tPij1, const phydbl *plk1, __m128d *_pmat1plk1, const __m128d *_tPij2, const phydbl *plk2, __m128d *_pmat2plk2, const int ns, __m128d *_plk0)
{
  unsigned int i;
  unsigned const int sz = (int)BYTE_ALIGN / 8;
  unsigned const int nblocks = ns / sz;
  
  for(i=0;i<ns;++i) if(plk1[i] > 1.0 || plk1[i] < 1.0 || plk2[i] > 1.0 || plk2[i] < 1.0) break; 

  if(i != ns)
    {     
      SSE_Matrix_Vect_Prod(_tPij1,plk1,ns,_pmat1plk1);
      SSE_Matrix_Vect_Prod(_tPij2,plk2,ns,_pmat2plk2);
      
      for(i=0;i<nblocks;++i) _plk0[i] = _mm_mul_pd(_pmat1plk1[i],_pmat2plk2[i]);
    }
  else
    {
      for(i=0;i<nblocks;++i) _plk0[i] = _mm_set1_pd(1.0);      
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void SSE_Matrix_Vect_Prod(const __m128d *_m_transpose, const phydbl *_v, const int ns, __m128d *_u)
{
  unsigned const int sz = (int)BYTE_ALIGN / 8;
  unsigned const int nblocks = ns / sz;  
  unsigned int i,j;  
  __m128d _x;


  _x = _mm_set1_pd(_v[0]);
  for(j=0;j<nblocks;++j) _u[j] = _mm_mul_pd(_m_transpose[j],_x);
  _m_transpose = _m_transpose + nblocks;
  
  for(i=1;i<ns;++i)
    {
      _x = _mm_set1_pd(_v[i]);
      for(j=0;j<nblocks;++j)
        {
          _u[j] = _mm_add_pd(_u[j],_mm_mul_pd(_m_transpose[j],_x));
        }
      _m_transpose = _m_transpose + nblocks;
    }

  /* for(i=0;i<nblocks;++i) _u[i] = _mm_setzero_pd(); */

  /* for(i=0;i<ns;++i) */
  /*   { */
  /*     _x = _mm_set1_pd(_v[i]); */
  /*     for(j=0;j<nblocks;++j) _u[j] = _mm_add_pd(_u[j],_mm_mul_pd(_m_transpose[j],_x)); */
  /*     _m_transpose = _m_transpose + nblocks; */
  /*   } */
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
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

#endif
