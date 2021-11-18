/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include "make.h"

//////////////////////////////////////////////////////////////

void Make_Tree_For_Lk(t_tree *tree)
{
  int i;
  calign *cdata;


  if(tree->is_mixt_tree == YES)
    {
      MIXT_Make_Tree_For_Lk(tree);
      return;
    }
  
  const unsigned int ns = tree->mod->ns;
  const unsigned int nsns =  ns * ns;

#if (defined(__AVX__) || defined(__AVX2__) || defined(__SSE__) || defined(__SSE2__) || defined(__SSE3__))
  const unsigned int sz = (int)BYTE_ALIGN / 8;
  const unsigned int ncatg = tree->mod->ras->n_catg;
#endif
  
  cdata = tree->data;
  assert(cdata);

  
  tree->c_lnL_sorted         = (phydbl *)mCalloc(tree->n_pattern,sizeof(phydbl));
  tree->cur_site_lk          = (phydbl *)mCalloc(tree->n_pattern,sizeof(phydbl));
  tree->old_site_lk          = (phydbl *)mCalloc(tree->n_pattern,sizeof(phydbl));
  tree->site_lk_cat          = (phydbl *)mCalloc(MAX(tree->mod->ras->n_catg,tree->mod->n_mixt_classes),sizeof(phydbl));
  tree->unscaled_site_lk_cat = (phydbl *)mCalloc(MAX(tree->mod->ras->n_catg,tree->mod->n_mixt_classes)*tree->n_pattern,sizeof(phydbl));
  tree->fact_sum_scale       = (int *)mCalloc(tree->n_pattern,sizeof(int));

  
#if (defined(__AVX__) || defined(__AVX2__))
#ifndef WIN32
  if(posix_memalign((void **)&tree->dot_prod,BYTE_ALIGN,(size_t)tree->n_pattern*tree->mod->ns*MAX(tree->mod->ras->n_catg,tree->mod->n_mixt_classes)*sizeof(phydbl))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  if(posix_memalign((void **)&tree->expl,BYTE_ALIGN,(size_t)3*MAX(tree->mod->ras->n_catg,tree->mod->n_mixt_classes)*tree->mod->ns*sizeof(phydbl))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  if(posix_memalign((void **)&tree->p_lk_left_pi,BYTE_ALIGN,(size_t) ns * sizeof(phydbl))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  if(posix_memalign((void **)&tree->l_ev,BYTE_ALIGN,(size_t) nsns * sizeof(phydbl))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  if(posix_memalign((void **)&tree->_tPij1,BYTE_ALIGN,(size_t)(ncatg * nsns / sz) * sizeof(__m256d))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  if(posix_memalign((void **)&tree->_tPij2,BYTE_ALIGN,(size_t)(ncatg * nsns / sz) * sizeof(__m256d))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  if(posix_memalign((void **)&tree->_pmat1plk1,BYTE_ALIGN,(size_t)ns / sz * sizeof(__m256d))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  if(posix_memalign((void **)&tree->_pmat2plk2,BYTE_ALIGN,(size_t)ns / sz * sizeof(__m256d))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  if(posix_memalign((void **)&tree->_plk0,BYTE_ALIGN,(size_t)(ns / sz) * sizeof(__m256d))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  if(posix_memalign((void **)&tree->_l_ev,BYTE_ALIGN,(size_t) nsns / sz * sizeof(__m256d))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  if(posix_memalign((void **)&tree->_r_ev,BYTE_ALIGN,(size_t) nsns / sz * sizeof(__m256d))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  if(posix_memalign((void **)&tree->_prod_left,BYTE_ALIGN,(size_t) ns / sz * sizeof(__m256d))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  if(posix_memalign((void **)&tree->_prod_rght,BYTE_ALIGN,(size_t) ns / sz * sizeof(__m256d))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  if(posix_memalign((void **)&tree->big_lk_array,BYTE_ALIGN,(size_t) ((3*tree->n_otu-2)*tree->data->crunch_len*MAX(tree->mod->ras->n_catg,tree->mod->n_mixt_classes)*tree->mod->ns + (2*tree->n_otu-1)*tree->mod->ras->n_catg*tree->mod->ns*tree->mod->ns) * sizeof(phydbl))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
#else
  tree->dot_prod     = _aligned_malloc(tree->n_pattern*tree->mod->ns*MAX(tree->mod->ras->n_catg,tree->mod->n_mixt_classes)*sizeof(phydbl),BYTE_ALIGN);
  tree->expl         = _aligned_malloc(3*MAX(tree->mod->ras->n_catg,tree->mod->n_mixt_classes)*tree->mod->ns*sizeof(phydbl),BYTE_ALIGN);
  tree->p_lk_left_pi = _aligned_malloc(ns * sizeof(phydbl),BYTE_ALIGN);
  tree->l_ev         = _aligned_malloc(ns * ns * sizeof(phydbl),BYTE_ALIGN);
  tree->_tPij1       = _aligned_malloc(ncatg * nsns / sz * sizeof(__m256d),BYTE_ALIGN);
  tree->_tPij2       = _aligned_malloc(ncatg * nsns / sz * sizeof(__m256d),BYTE_ALIGN);
  tree->_pmat1plk1   = _aligned_malloc(ns / sz * sizeof(__m256d),BYTE_ALIGN);
  tree->_pmat2plk2   = _aligned_malloc(ns / sz * sizeof(__m256d),BYTE_ALIGN);
  tree->_plk0        = _aligned_malloc(ns / sz * sizeof(__m256d),BYTE_ALIGN);
  tree->_l_ev        = _aligned_malloc(ns * ns / sz * sizeof(__m256d),BYTE_ALIGN);
  tree->_r_ev        = _aligned_malloc(ns * ns / sz * sizeof(__m256d),BYTE_ALIGN);
  tree->_prod_left   = _aligned_malloc(ns / sz * sizeof(__m256d),BYTE_ALIGN);
  tree->_prod_rght   = _aligned_malloc(ns / sz * sizeof(__m256d),BYTE_ALIGN);
  tree->big_lk_array     = _aligned_malloc(((3*tree->n_otu-2)*tree->data->crunch_len*MAX(tree->mod->ras->n_catg,tree->mod->n_mixt_classes)*tree->mod->ns + (2*tree->n_otu-1)*tree->mod->ras->n_catg*tree->mod->ns*tree->mod->ns) * sizeof(phydbl),BYTE_ALIGN);
#endif
#elif (defined(__SSE__) || defined(__SSE2__) || defined(__SSE3__))
#ifndef WIN32
  if(posix_memalign((void **)&tree->dot_prod,BYTE_ALIGN,(size_t)tree->n_pattern*tree->mod->ns*MAX(tree->mod->ras->n_catg,tree->mod->n_mixt_classes)*sizeof(phydbl))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  if(posix_memalign((void **)&tree->expl,BYTE_ALIGN,(size_t)3*MAX(tree->mod->ras->n_catg,tree->mod->n_mixt_classes)*tree->mod->ns*sizeof(phydbl))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  if(posix_memalign((void **)&tree->p_lk_left_pi,BYTE_ALIGN,(size_t) ns * sizeof(phydbl))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  if(posix_memalign((void **)&tree->l_ev,BYTE_ALIGN,(size_t) nsns * sizeof(phydbl))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  if(posix_memalign((void **)&tree->_tPij1,BYTE_ALIGN,(size_t)(ncatg * nsns / sz) * sizeof(__m128d))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  if(posix_memalign((void **)&tree->_tPij2,BYTE_ALIGN,(size_t)(ncatg * nsns / sz) * sizeof(__m128d))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  if(posix_memalign((void **)&tree->_pmat1plk1,BYTE_ALIGN,(size_t)ns / sz * sizeof(__m128d))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  if(posix_memalign((void **)&tree->_pmat2plk2,BYTE_ALIGN,(size_t)ns / sz * sizeof(__m128d))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  if(posix_memalign((void **)&tree->_plk0,BYTE_ALIGN,(size_t)(ns / sz) * sizeof(__m128d))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  if(posix_memalign((void **)&tree->_l_ev,BYTE_ALIGN,(size_t) nsns / sz * sizeof(__m128d))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  if(posix_memalign((void **)&tree->_r_ev,BYTE_ALIGN,(size_t) nsns / sz * sizeof(__m128d))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  if(posix_memalign((void **)&tree->_prod_left,BYTE_ALIGN,(size_t) ns / sz * sizeof(__m128d))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  if(posix_memalign((void **)&tree->_prod_rght,BYTE_ALIGN,(size_t) ns / sz * sizeof(__m128d))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  if(posix_memalign((void **)&tree->big_lk_array,BYTE_ALIGN,(size_t) ((3*tree->n_otu-2)*tree->data->crunch_len*MAX(tree->mod->ras->n_catg,tree->mod->n_mixt_classes)*tree->mod->ns + (2*tree->n_otu-1)*tree->mod->ras->n_catg*tree->mod->ns*tree->mod->ns) * sizeof(phydbl))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
#else
  tree->dot_prod     = _aligned_malloc(tree->n_pattern*tree->mod->ns*MAX(tree->mod->ras->n_catg,tree->mod->n_mixt_classes)*sizeof(phydbl),BYTE_ALIGN);
  tree->expl         = _aligned_malloc(3*MAX(tree->mod->ras->n_catg,tree->mod->n_mixt_classes)*tree->mod->ns*sizeof(phydbl),BYTE_ALIGN);
  tree->p_lk_left_pi = _aligned_malloc(ns * sizeof(phydbl),BYTE_ALIGN);
  tree->l_ev         = _aligned_malloc(ns * ns * sizeof(phydbl),BYTE_ALIGN);
  tree->_tPij1       = _aligned_malloc(ncatg * nsns / sz * sizeof(__m128d),BYTE_ALIGN);
  tree->_tPij2       = _aligned_malloc(ncatg * nsns / sz * sizeof(__m128d),BYTE_ALIGN);
  tree->_pmat1plk1   = _aligned_malloc(ns / sz * sizeof(__m128d),BYTE_ALIGN);
  tree->_pmat2plk2   = _aligned_malloc(ns / sz * sizeof(__m128d),BYTE_ALIGN);
  tree->_plk0        = _aligned_malloc(ns / sz * sizeof(__m128d),BYTE_ALIGN);
  tree->_l_ev        = _aligned_malloc(ns * ns / sz * sizeof(__m128d),BYTE_ALIGN);
  tree->_r_ev        = _aligned_malloc(ns * ns / sz * sizeof(__m128d),BYTE_ALIGN);
  tree->_prod_left   = _aligned_malloc(ns / sz * sizeof(__m128d),BYTE_ALIGN);
  tree->_prod_rght   = _aligned_malloc(ns / sz * sizeof(__m128d),BYTE_ALIGN);
  tree->big_lk_array     = _aligned_malloc(((3*tree->n_otu-2)*tree->data->crunch_len*MAX(tree->mod->ras->n_catg,tree->mod->n_mixt_classes)*tree->mod->ns + (2*tree->n_otu-1)*tree->mod->ras->n_catg*tree->mod->ns*tree->mod->ns) * sizeof(phydbl),BYTE_ALIGN);
#endif
#elif (!(defined(__AVX__) || defined(__AVX2__) || defined(__SSE__) || defined(__SSE2__) || defined(__SSE3__)))
  tree->dot_prod = (phydbl *)mCalloc(tree->n_pattern*tree->mod->ns*MAX(tree->mod->ras->n_catg,tree->mod->n_mixt_classes),sizeof(phydbl));
  tree->expl = (phydbl *)mCalloc(3*MAX(tree->mod->ras->n_catg,tree->mod->n_mixt_classes)*tree->mod->ns,sizeof(phydbl));
  tree->p_lk_left_pi = (phydbl *)mCalloc(ns,sizeof(phydbl));
  tree->l_ev = (phydbl *)mCalloc(nsns,sizeof(phydbl));
  tree->big_lk_array = (phydbl *)mCalloc(((3*tree->n_otu-2)*tree->data->crunch_len*MAX(tree->mod->ras->n_catg,tree->mod->n_mixt_classes)*tree->mod->ns + (2*tree->n_otu-1)*tree->mod->ras->n_catg*tree->mod->ns*tree->mod->ns),sizeof(phydbl));
#endif
  
  tree->big_lk_array_pos = 0;
  
  tree->log_lks_aLRT = (phydbl **)mCalloc(3,sizeof(phydbl *));
  for(i=0;i<3;i++) tree->log_lks_aLRT[i] = (phydbl *)mCalloc(tree->data->init_len,sizeof(phydbl));

  for(i=0;i<2*tree->n_otu-1;++i) Make_Edge_NNI(tree->a_edges[i]);

  Make_Extra_Edge_Lk(tree);
  
  if(tree->is_mixt_tree == NO)
    {
      for(i=0;i<2*tree->n_otu-1;++i) Make_Edge_Lk(tree->a_edges[i],tree);
      for(i=0;i<2*tree->n_otu-2;++i) Make_Node_Lk(tree->a_nodes[i]);
      for(i=0;i<2*tree->n_otu-1;++i) Make_Edge_Loc(tree->a_edges[i],tree);

      Init_Partial_Lk_Tips_Double(tree);
      Init_Partial_Lk_Loc(tree);
      
      if(tree->n_root != NULL)
        {
          /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!1 */
          /* Free_Edge_Lk_Rght(tree->n_root->b[1]); */
          /* Free_Edge_Lk_Rght(tree->n_root->b[2]); */

          Free_Edge_Loc_Rght(tree->n_root->b[1]);
          Free_Edge_Loc_Rght(tree->n_root->b[2]);

          tree->n_root->b[1]->p_lk_rght = tree->e_root->p_lk_left;
          tree->n_root->b[2]->p_lk_rght = tree->e_root->p_lk_rght;
          
          tree->n_root->b[1]->p_lk_tip_r = tree->e_root->p_lk_tip_l;      
          tree->n_root->b[2]->p_lk_tip_r = tree->e_root->p_lk_tip_r;
          
          tree->n_root->b[1]->div_post_pred_rght = tree->e_root->div_post_pred_left;
          tree->n_root->b[2]->div_post_pred_rght = tree->e_root->div_post_pred_rght;
          
          tree->n_root->b[1]->sum_scale_rght = tree->e_root->sum_scale_left;
          tree->n_root->b[2]->sum_scale_rght = tree->e_root->sum_scale_rght;
          
          tree->n_root->b[1]->sum_scale_rght_cat = tree->e_root->sum_scale_left_cat;
          tree->n_root->b[2]->sum_scale_rght_cat = tree->e_root->sum_scale_rght_cat;
          
          tree->n_root->b[1]->patt_id_rght = tree->e_root->patt_id_left;
          tree->n_root->b[2]->patt_id_rght = tree->e_root->patt_id_rght;
        }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Make_Tree_For_Pars(t_tree *tree)
{
  int i;
  calign *cdata;

  if(tree->is_mixt_tree == YES)
    {
      MIXT_Make_Tree_For_Pars(tree);
      return;
    }
  
  cdata = tree->data;
  assert(cdata);

  assert(tree->mod);
  
  tree->site_pars = (int *)mCalloc(tree->n_pattern,sizeof(int));
  tree->step_mat = (int *)mCalloc(tree->mod->ns * tree->mod->ns,sizeof(int));

  for(i=0;i<2*tree->n_otu-1;++i) Make_Edge_Pars(tree->a_edges[i],tree);

  Init_Ui_Tips(tree);
  Init_Partial_Pars_Tips(tree); /* Must be called after Init_Ui_Tips is called */

  if(tree->n_root)
    {
      Free_Edge_Pars_Rght(tree->a_edges[2*tree->n_otu-3]);
      Free_Edge_Pars_Rght(tree->a_edges[2*tree->n_otu-2]);
    }

  Get_Step_Mat(tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Make_All_Edges_Lk(t_node *a, t_node *d, t_tree *tree)
{
  int i;

  for(i=0;i<3;i++) 
    if((a->v[i]) && (a->v[i] == d)) 
      Make_Edge_Lk(a->b[i],tree);

  if(d->tax) return;
  else
    {
      for(i=0;i<3;i++)
        {
          if(d->v[i] != a && d->b[i] != tree->e_root)
            Make_All_Edges_Lk(d,d->v[i],tree);
        }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

t_edge *Make_Edge_Light(t_node *a, t_node *d, int num)
{
  t_edge *b;
  
  b = (t_edge *)mCalloc(1,sizeof(t_edge));

  b->l = (scalar_dbl *)mCalloc(1,sizeof(scalar_dbl));
  Init_Scalar_Dbl(b->l);

  b->l_old = (scalar_dbl *)mCalloc(1,sizeof(scalar_dbl));
  Init_Scalar_Dbl(b->l_old);

  b->l_var = (scalar_dbl *)mCalloc(1,sizeof(scalar_dbl));
  Init_Scalar_Dbl(b->l_var);

  b->l_var_old = (scalar_dbl *)mCalloc(1,sizeof(scalar_dbl));
  Init_Scalar_Dbl(b->l_var_old);

  Init_Edge_Light(b,num);

  if(a && b)
    {
      b->left = a;
      b->rght = d;
      if(a->tax) {b->rght = a; b->left = d;} /* root */
      /* a tip is necessary on the right side of the t_edge */

      (b->left == a)?
        (Set_Edge_Dirs(b,a,d,NULL)):
        (Set_Edge_Dirs(b,d,a,NULL));

      assert(b->l_r > -1);
      assert(b->r_l > -1);
      
      b->l_old->v = b->l->v;
    }
  else
    {
      b->left = NULL;
      b->rght = NULL;
    }

  return b;

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Make_Edge_Pars(t_edge *b, t_tree *tree)
{
  assert(b);
  Make_Edge_Pars_Left(b,tree);
  Make_Edge_Pars_Rght(b,tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Make_Edge_Pars_Left(t_edge *b, t_tree *tree)
{
  b->pars_l = (int *)mCalloc(tree->data->crunch_len,sizeof(int));
  b->ui_l = (int *)mCalloc(tree->data->crunch_len,sizeof(int));
  b->p_pars_l = (int *)mCalloc(tree->data->crunch_len*tree->mod->ns,sizeof(int ));
  b->n_diff_states_l = (int *)mCalloc(tree->mod->ns,sizeof(int ));
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Make_Edge_Pars_Rght(t_edge *b, t_tree *tree)
{
  b->pars_r = (int *)mCalloc(tree->data->crunch_len,sizeof(int));
  b->ui_r = (int *)mCalloc(tree->data->crunch_len,sizeof(int));
  b->p_pars_r = (int *)mCalloc(tree->data->crunch_len*tree->mod->ns,sizeof(int ));
  b->n_diff_states_r = (int *)mCalloc(tree->mod->ns,sizeof(int ));
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Make_Edge_Lk(t_edge *b, t_tree *tree)
{
  if(tree->is_mixt_tree)
    {
      PhyML_Printf("\n== Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
      Warn_And_Exit("");
    }

  b->l_old->v = b->l->v;

/* #if (defined(__AVX__) || defined(__AVX2__) || defined(__SSE__) || defined(__SSE2__) || defined(__SSE3__)) */
/* #ifndef WIN32 */
/*   if(posix_memalign((void *)&b->Pij_rr,BYTE_ALIGN,(size_t)tree->mod->ras->n_catg*tree->mod->ns*tree->mod->ns*sizeof(phydbl))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__); */
/*   if(posix_memalign((void *)&b->tPij_rr,BYTE_ALIGN,(size_t)tree->mod->ras->n_catg*tree->mod->ns*tree->mod->ns*sizeof(phydbl))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__); */
/* #else */
/*   b->Pij_rr = _aligned_malloc(tree->mod->ras->n_catg*tree->mod->ns*tree->mod->ns*sizeof(phydbl),BYTE_ALIGN); */
/*   b->tPij_rr = _aligned_malloc(tree->mod->ras->n_catg*tree->mod->ns*tree->mod->ns*sizeof(phydbl),BYTE_ALIGN); */
/* #endif */
/* #else */
/*   b->Pij_rr = (phydbl *)mCalloc(tree->mod->ras->n_catg*tree->mod->ns*tree->mod->ns,sizeof(phydbl)); */
/*   b->tPij_rr = (phydbl *)mCalloc(tree->mod->ras->n_catg*tree->mod->ns*tree->mod->ns,sizeof(phydbl)); */
/* #endif */


  b->Pij_rr = tree->big_lk_array + tree->big_lk_array_pos;
  tree->big_lk_array_pos += tree->mod->ras->n_catg*tree->mod->ns*tree->mod->ns;

  b->tPij_rr = tree->big_lk_array + tree->big_lk_array_pos;
  tree->big_lk_array_pos += tree->mod->ras->n_catg*tree->mod->ns*tree->mod->ns;

  
  Make_Edge_Lk_Left(b,tree);
  Make_Edge_Lk_Rght(b,tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Make_Edge_Lk_Left(t_edge *b, t_tree *tree)
{
  int ns = tree->mod->ns;

  b->div_post_pred_left = (short int *)mCalloc(ns,sizeof(short int));

  b->sum_scale_left_cat = (int *)mCalloc(MAX(tree->mod->ras->n_catg,tree->mod->n_mixt_classes),sizeof(int));

  if(b->left && !b->left->tax)
    b->sum_scale_left = (int *)mCalloc(tree->data->crunch_len*MAX(tree->mod->ras->n_catg,tree->mod->n_mixt_classes),sizeof(int));
  else
    b->sum_scale_left = NULL;

  if(b->left)
    {
      if((!b->left->tax) || (tree->mod->s_opt->greedy))
        {

/* #if (defined(__AVX__) || defined(__AVX2__) || defined(__SSE__) || defined(__SSE2__) || defined(__SSE3__)) */
/* #ifndef WIN32 */
/*           if(posix_memalign((void **)&b->p_lk_left,BYTE_ALIGN,(size_t)tree->data->crunch_len*MAX(tree->mod->ras->n_catg,tree->mod->n_mixt_classes)*tree->mod->ns*sizeof(phydbl))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__); */
/* #else */
/*           b->p_lk_left = _aligned_malloc(tree->data->crunch_len*MAX(tree->mod->ras->n_catg,tree->mod->n_mixt_classes)*tree->mod->ns*sizeof(phydbl),BYTE_ALIGN); */
/* #endif */
/* #else */
/*           b->p_lk_left = (phydbl *)mCalloc(tree->data->crunch_len*MAX(tree->mod->ras->n_catg,tree->mod->n_mixt_classes)*tree->mod->ns,sizeof(phydbl)); */
/* #endif */

          b->p_lk_tip_l = NULL;

          /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
          b->p_lk_left = tree->big_lk_array + tree->big_lk_array_pos;
          tree->big_lk_array_pos += tree->data->crunch_len*MAX(tree->mod->ras->n_catg,tree->mod->n_mixt_classes)*tree->mod->ns;
        }
      else if(b->left->tax)
        {
          b->p_lk_left   = NULL;
          
#if (defined(__AVX__) || defined(__AVX2__) || defined(__SSE__) || defined(__SSE2__) || defined(__SSE3__))
#ifndef WIN32
          if(posix_memalign((void **)&b->p_lk_tip_l,BYTE_ALIGN,(size_t)tree->data->crunch_len*tree->mod->ns*sizeof(phydbl))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
#else
          b->p_lk_tip_l = _aligned_malloc(tree->data->crunch_len*tree->mod->ns*sizeof(phydbl),BYTE_ALIGN);
#endif
#else
          b->p_lk_tip_l  = (phydbl *)mCalloc(tree->data->crunch_len*tree->mod->ns,sizeof(phydbl));
#endif
        }
    }
  else
    {
      b->p_lk_left  = NULL;
      b->p_lk_tip_l = NULL;
    }

  if(b->num >= 2*tree->n_otu-3)
    {
      b->sum_scale_left = (int *)mCalloc(tree->data->crunch_len*MAX(tree->mod->ras->n_catg,tree->mod->n_mixt_classes),sizeof(int));

/* #if (defined(__AVX__) || defined(__AVX2__) || defined(__SSE__) || defined(__SSE2__) || defined(__SSE3__)) */
/* #ifndef WIN32 */
/*       if(posix_memalign((void **)&b->p_lk_left,BYTE_ALIGN,(size_t)tree->data->crunch_len*MAX(tree->mod->ras->n_catg,tree->mod->n_mixt_classes)*tree->mod->ns*sizeof(phydbl))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__); */
/* #else */
/*       b->p_lk_left = _aligned_malloc(tree->data->crunch_len*MAX(tree->mod->ras->n_catg,tree->mod->n_mixt_classes)*tree->mod->ns*sizeof(phydbl),BYTE_ALIGN); */
/* #endif */
/* #else */
/*       b->p_lk_left      = (phydbl *)mCalloc(tree->data->crunch_len*MAX(tree->mod->ras->n_catg,tree->mod->n_mixt_classes)*tree->mod->ns,sizeof(phydbl)); */
/* #endif */

      /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
      b->p_lk_left = tree->big_lk_array + tree->big_lk_array_pos;
      tree->big_lk_array_pos += tree->data->crunch_len*MAX(tree->mod->ras->n_catg,tree->mod->n_mixt_classes)*tree->mod->ns;
    }

  b->patt_id_left  = (int *)mCalloc(tree->data->crunch_len,sizeof(int));
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Make_Edge_Lk_Rght(t_edge *b, t_tree *tree)
{
  int ns = tree->mod->ns;

  b->div_post_pred_rght = (short int *)mCalloc(ns,sizeof(short int));

  b->sum_scale_rght_cat = (int *)mCalloc(MAX(tree->mod->ras->n_catg,tree->mod->n_mixt_classes),sizeof(int));

  if(b->rght && b->rght->tax == NO)
    b->sum_scale_rght = (int *)mCalloc(tree->data->crunch_len*MAX(tree->mod->ras->n_catg,tree->mod->n_mixt_classes),sizeof(int));
  else
    b->sum_scale_rght = NULL;


  if(b->rght)
    {
      if((!b->rght->tax) || (tree->mod->s_opt->greedy))
        {
/* #if (defined(__AVX__) || defined(__AVX2__) || defined(__SSE__) || defined(__SSE2__) || defined(__SSE3__)) */
/* #ifndef WIN32 */
/* 	  if(posix_memalign((void **)&b->p_lk_rght,BYTE_ALIGN,(size_t)tree->data->crunch_len*MAX(tree->mod->ras->n_catg,tree->mod->n_mixt_classes)*tree->mod->ns*sizeof(phydbl))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__); */
/* #else */
/*           b->p_lk_rght = _aligned_malloc(tree->data->crunch_len*MAX(tree->mod->ras->n_catg,tree->mod->n_mixt_classes)*tree->mod->ns*sizeof(phydbl),BYTE_ALIGN); */
/* #endif */
/* #else */
/*           b->p_lk_rght = (phydbl *)mCalloc(tree->data->crunch_len*MAX(tree->mod->ras->n_catg,tree->mod->n_mixt_classes)*tree->mod->ns,sizeof(phydbl)); */
/* #endif */

          b->p_lk_tip_r = NULL;


          /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
          b->p_lk_rght = tree->big_lk_array + tree->big_lk_array_pos;
          tree->big_lk_array_pos += tree->data->crunch_len*MAX(tree->mod->ras->n_catg,tree->mod->n_mixt_classes)*tree->mod->ns;
        }
      else if(b->rght->tax)
        {
#if (defined(__AVX__) || defined(__AVX2__) || defined(__SSE__) || defined(__SSE2__) || defined(__SSE3__))
#ifndef WIN32
          if(posix_memalign((void **)&b->p_lk_tip_r,BYTE_ALIGN,(size_t)tree->data->crunch_len*tree->mod->ns*sizeof(phydbl))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
#else
          b->p_lk_tip_r = _aligned_malloc(tree->data->crunch_len*tree->mod->ns*sizeof(phydbl),BYTE_ALIGN);
#endif
#else
          b->p_lk_tip_r  = (phydbl *)mCalloc(tree->data->crunch_len*tree->mod->ns,sizeof(phydbl));
#endif
          b->p_lk_rght = NULL;
        }
    }
  else
    {
      b->p_lk_rght  = NULL;
      b->p_lk_tip_r = NULL;
    }

  if(b->num >= 2*tree->n_otu-3)
    {
      b->sum_scale_rght = (int *)mCalloc(tree->data->crunch_len*MAX(tree->mod->ras->n_catg,tree->mod->n_mixt_classes),sizeof(int));

/* #if (defined(__AVX__) || defined(__AVX2__) || defined(__SSE__) || defined(__SSE2__) || defined(__SSE3__)) */
/* #ifndef WIN32 */
/*       if(posix_memalign((void **)&b->p_lk_rght,BYTE_ALIGN,(size_t)tree->data->crunch_len*MAX(tree->mod->ras->n_catg,tree->mod->n_mixt_classes)*tree->mod->ns*sizeof(phydbl))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__); */
/* #else */
/*       b->p_lk_rght      = _aligned_malloc(tree->data->crunch_len*MAX(tree->mod->ras->n_catg,tree->mod->n_mixt_classes)*tree->mod->ns*sizeof(phydbl),BYTE_ALIGN); */
/* #endif */

/* #else */
/*       b->p_lk_rght      = (phydbl *)mCalloc(tree->data->crunch_len*MAX(tree->mod->ras->n_catg,tree->mod->n_mixt_classes)*tree->mod->ns,sizeof(phydbl)); */
/* #endif */

      /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
      b->p_lk_rght = tree->big_lk_array + tree->big_lk_array_pos;
      tree->big_lk_array_pos += tree->data->crunch_len*MAX(tree->mod->ras->n_catg,tree->mod->n_mixt_classes)*tree->mod->ns;
    }

  b->patt_id_rght  = (int *)mCalloc(tree->data->crunch_len,sizeof(int));
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Make_Extra_Edge_Lk(t_tree *tree)
{
  int ns = tree->mod->ns;

  tree->div_post_pred_extra_0 = (short int *)mCalloc(ns,sizeof(short int));

  tree->sum_scale_cat_extra_0 = (int *)mCalloc(MAX(tree->mod->ras->n_catg,tree->mod->n_mixt_classes),sizeof(int));

  tree->sum_scale_extra_0 = (int *)mCalloc(tree->data->crunch_len*MAX(tree->mod->ras->n_catg,tree->mod->n_mixt_classes),sizeof(int));

#if (defined(__AVX__) || defined(__AVX2__) || defined(__SSE__) || defined(__SSE2__) || defined(__SSE3__))
#ifndef WIN32
  if(posix_memalign((void **)&tree->p_lk_extra_0,BYTE_ALIGN,(size_t)tree->data->crunch_len*MAX(tree->mod->ras->n_catg,tree->mod->n_mixt_classes)*tree->mod->ns*sizeof(phydbl))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
#else
  tree->p_lk_extra_0 = _aligned_malloc(tree->data->crunch_len*MAX(tree->mod->ras->n_catg,tree->mod->n_mixt_classes)*tree->mod->ns*sizeof(phydbl),BYTE_ALIGN);
#endif
#else
  tree->p_lk_extra_0 = (phydbl *)mCalloc(tree->data->crunch_len*MAX(tree->mod->ras->n_catg,tree->mod->n_mixt_classes)*tree->mod->ns,sizeof(phydbl));
#endif

  
#if (defined(__AVX__) || defined(__AVX2__) || defined(__SSE__) || defined(__SSE2__) || defined(__SSE3__))
#ifndef WIN32
  if(posix_memalign((void **)&tree->p_lk_tip_extra_0,BYTE_ALIGN,(size_t)tree->data->crunch_len*tree->mod->ns*sizeof(phydbl))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
#else
  tree->p_lk_tip_extra_0 = _aligned_malloc(tree->data->crunch_len*tree->mod->ns*sizeof(phydbl),BYTE_ALIGN);
#endif
#else
  tree->p_lk_tip_extra_0  = (phydbl *)mCalloc(tree->data->crunch_len*tree->mod->ns,sizeof(phydbl));
#endif
        
  tree->patt_id_extra_0  = (int *)mCalloc(tree->data->crunch_len,sizeof(int));

  tree->div_post_pred_extra_1 = (short int *)mCalloc(ns,sizeof(short int));

  tree->sum_scale_cat_extra_1 = (int *)mCalloc(MAX(tree->mod->ras->n_catg,tree->mod->n_mixt_classes),sizeof(int));

  tree->sum_scale_extra_1 = (int *)mCalloc(tree->data->crunch_len*MAX(tree->mod->ras->n_catg,tree->mod->n_mixt_classes),sizeof(int));

#if (defined(__AVX__) || defined(__AVX2__) || defined(__SSE__) || defined(__SSE2__) || defined(__SSE3__))
#ifndef WIN32
  if(posix_memalign((void **)&tree->p_lk_extra_1,BYTE_ALIGN,(size_t)tree->data->crunch_len*MAX(tree->mod->ras->n_catg,tree->mod->n_mixt_classes)*tree->mod->ns*sizeof(phydbl))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
#else
  tree->p_lk_extra_1 = _aligned_malloc(tree->data->crunch_len*MAX(tree->mod->ras->n_catg,tree->mod->n_mixt_classes)*tree->mod->ns*sizeof(phydbl),BYTE_ALIGN);
#endif
#else
  tree->p_lk_extra_1 = (phydbl *)mCalloc(tree->data->crunch_len*MAX(tree->mod->ras->n_catg,tree->mod->n_mixt_classes)*tree->mod->ns,sizeof(phydbl));
#endif

  
#if (defined(__AVX__) || defined(__AVX2__) || defined(__SSE__) || defined(__SSE2__) || defined(__SSE3__))
#ifndef WIN32
  if(posix_memalign((void **)&tree->p_lk_tip_extra_1,BYTE_ALIGN,(size_t)tree->data->crunch_len*tree->mod->ns*sizeof(phydbl))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
#else
  tree->p_lk_tip_extra_1 = _aligned_malloc(tree->data->crunch_len*tree->mod->ns*sizeof(phydbl),BYTE_ALIGN);
#endif
#else
  tree->p_lk_tip_extra_1  = (phydbl *)mCalloc(tree->data->crunch_len*tree->mod->ns,sizeof(phydbl));
#endif
        
  tree->patt_id_extra_1  = (int *)mCalloc(tree->data->crunch_len,sizeof(int));

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Make_Edge_Loc(t_edge *b, t_tree *tree)
{
  Make_Edge_Loc_Left(b,tree);
  Make_Edge_Loc_Rght(b,tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Make_Edge_Loc_Rght(t_edge *b, t_tree *tree)
{
  b->p_lk_loc_rght = (int *)mCalloc(tree->data->crunch_len,sizeof(int));
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Make_Edge_Loc_Left(t_edge *b, t_tree *tree)
{
  b->p_lk_loc_left = (int *)mCalloc(tree->data->crunch_len,sizeof(int));
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Make_Edge_NNI(t_edge *b)
{
  b->nni    = Make_NNI();
  b->nni->b = b;
  b->nni->left = b->left;
  b->nni->rght = b->rght;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

t_nni *Make_NNI()
{
  t_nni *a_nni;

  a_nni = (t_nni *)mCalloc(1,sizeof(t_nni));
  
  Init_NNI(a_nni);
  return a_nni;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

t_node *Make_Node_Light(int num)
{
  t_node *n;
  
  n           = (t_node *)mCalloc(1,sizeof(t_node));
  n->v        = (t_node **)mCalloc(3,sizeof(t_node *));
  n->b        = (t_edge **)mCalloc(3,sizeof(t_edge *));
  n->score    = (phydbl *)mCalloc(3,sizeof(phydbl));
  n->s_ingrp  = (int *)mCalloc(3,sizeof(int));
  n->s_outgrp = (int *)mCalloc(3,sizeof(int));
  n->cal      = (t_cal **)mCalloc(MAX_N_CAL,sizeof(t_cal *));

  Init_Node_Light(n,num);

  return n;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Make_Node_Lk(t_node *n)
{
/*   n->n_ex_nodes = (int *)mCalloc(2,sizeof(int)); */
  return;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

nexcom **Make_Nexus_Com()
{
  nexcom **com;
  int i;

  com = (nexcom **)mCalloc(N_MAX_NEX_COM,sizeof(nexcom *));

  for(i=0;i<N_MAX_NEX_COM;i++)
    {
      com[i]       = (nexcom *)mCalloc(1,sizeof(nexcom));
      com[i]->name = (char *)mCalloc(T_MAX_NEX_COM,sizeof(char));
      com[i]->parm = (nexparm **)mCalloc(N_MAX_NEX_PARM,sizeof(nexparm *));
    }

  return com;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

nexparm *Make_Nexus_Parm()
{
  nexparm *parm;

  parm        = (nexparm *)mCalloc(1,sizeof(nexparm));
  parm->name  = (char *)mCalloc(T_MAX_TOKEN,sizeof(char ));
  parm->value = (char *)mCalloc(T_MAX_TOKEN,sizeof(char ));

  return parm;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

matrix *Make_Mat(int n_otu)
{
  matrix *mat;
  int i;

  mat = (matrix *)mCalloc(1,sizeof(matrix));

  mat->n_otu = n_otu;

  mat->P        = (phydbl **)mCalloc(n_otu,sizeof(phydbl *));
  mat->Q        = (phydbl **)mCalloc(n_otu,sizeof(phydbl *));
  mat->dist     = (phydbl **)mCalloc(n_otu,sizeof(phydbl *));
  mat->on_off   = (int *)mCalloc(n_otu,sizeof(int));
  mat->name     = (char **)mCalloc(n_otu,sizeof(char *));
  mat->tip_node = (t_node **)mCalloc(n_otu,sizeof(t_node *));


  for(i=0;i<n_otu;i++)
    {
      mat->P[i]    = (phydbl *)mCalloc(n_otu,sizeof(phydbl));
      mat->Q[i]    = (phydbl *)mCalloc(n_otu,sizeof(phydbl));
      mat->dist[i] = (phydbl *)mCalloc(n_otu,sizeof(phydbl));
      mat->name[i] = (char *)mCalloc(T_MAX_NAME,sizeof(char));
    }

  return mat;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

t_tree *Make_Tree_From_Scratch(int n_otu, calign *data)
{  
  t_tree *tree;

  tree = Make_Tree(n_otu);
  Init_Tree(tree,n_otu);
  Make_All_Tree_Nodes(tree);
  Make_All_Tree_Edges(tree);
  Make_Tree_Path(tree);
  if(data)
    {
      Copy_Tax_Names_To_Tip_Labels(tree,data);
      tree->data = data;
    }
  
  
#ifdef BEAGLE
  //offset the branch's partial indices because BEAGLE insists on first storing the tips/taxa
  int num_branches = 2*tree->n_otu-1;
  int i;
  for(i=0;i<2*tree->n_otu-1;++i)
    {
      //For edgeX, its "left" partial lies at index `num_tax + edgeX->num"
      tree->a_edges[i]->p_lk_left_idx = tree->n_otu + tree->a_edges[i]->p_lk_left_idx;
      //For edgeX, its "right" partial lies at index `num_tax + edgeX->num + num_branches"
      tree->a_edges[i]->p_lk_rght_idx = tree->n_otu + tree->a_edges[i]->p_lk_left_idx + num_branches;
    }
#endif
  
  return tree;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

t_tree *Make_Tree(int n_otu)
{
  t_tree *tree;
  tree = (t_tree *)mCalloc(1,sizeof(t_tree ));
  tree->t_dir = (short int *)mCalloc((2*n_otu-2)*(2*n_otu-2),sizeof(int));
  return tree;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Make_Tree_Path(t_tree *tree)
{
  tree->curr_path = (t_node **)mCalloc(tree->n_otu,sizeof(t_node *));
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Make_All_Tree_Nodes(t_tree *tree)
{
  int i;

  tree->a_nodes = (t_node **)mCalloc(2*tree->n_otu-1,sizeof(t_node *));

  for(i=0;i<2*tree->n_otu-1;++i)
    {
      tree->a_nodes[i] = (t_node *)Make_Node_Light(i);
      if(i < tree->n_otu) tree->a_nodes[i]->tax = YES;
      else                tree->a_nodes[i]->tax = NO;
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Make_All_Tree_Edges(t_tree *tree)
{
  int i;
  tree->a_edges = (t_edge **)mCalloc(2*tree->n_otu-1,sizeof(t_edge *));
  for(i=0;i<2*tree->n_otu-1;++i) tree->a_edges[i] = (t_edge *)Make_Edge_Light(NULL,NULL,i);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

calign *Make_Calign(int n_otu, int crunch_len, int state_len, int init_len, char **sp_names_in, int n_rm, char **sp_names_out)
{
  calign *cdata;
  int j;

  cdata                = (calign *)mCalloc(1,sizeof(calign));
  cdata->c_seq         = (align **)mCalloc(n_otu,sizeof(align *));
  cdata->obs_state_frq = (phydbl *)mCalloc(T_MAX_ALPHABET,sizeof(phydbl));
  cdata->wght          = (phydbl *)mCalloc(crunch_len,sizeof(phydbl));
  cdata->ambigu        = (short int *)mCalloc(crunch_len,sizeof(short int));
  cdata->invar         = (short int *)mCalloc(crunch_len,sizeof(short int));
  cdata->sitepatt      = (int *)mCalloc(init_len,sizeof(int ));

  if(n_rm > 0) cdata->c_seq_rm = (align **)mCalloc(n_rm,sizeof(align *));

  for(j=0;j<n_otu;j++)
    {
      cdata->c_seq[j]            = (align *)mCalloc(1,sizeof(align));
      cdata->c_seq[j]->name      = (char *)mCalloc((int)(strlen(sp_names_in[j])+1),sizeof(char));
      strcpy(cdata->c_seq[j]->name,sp_names_in[j]);
      cdata->c_seq[j]->state     = (char *)mCalloc(crunch_len*state_len+1,sizeof(char));
      cdata->c_seq[j]->d_state   = (short int *)mCalloc(crunch_len*state_len,sizeof(short int));
      cdata->c_seq[j]->is_ambigu = (short int *)mCalloc(crunch_len,sizeof(short int));
    }

  for(j=0;j<n_rm;j++)
    {
      cdata->c_seq_rm[j]            = (align *)mCalloc(1,sizeof(align));
      cdata->c_seq_rm[j]->name      = (char *)mCalloc((int)(strlen(sp_names_out[j])+1),sizeof(char));
      strcpy(cdata->c_seq_rm[j]->name,sp_names_out[j]);
      cdata->c_seq_rm[j]->state     = (char *)mCalloc(crunch_len*state_len+1,sizeof(char));
      cdata->c_seq_rm[j]->d_state   = (short int *)mCalloc(crunch_len*state_len,sizeof(short int));
      cdata->c_seq_rm[j]->is_ambigu = (short int *)mCalloc(crunch_len,sizeof(short int));
    }
  
  return cdata;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


t_treelist *Make_Treelist(int list_size)
{
  t_treelist *tlist;

  tlist = (t_treelist *)mCalloc(1,sizeof(t_treelist));
  tlist->list_size = list_size;
  tlist->tree = (t_tree **)mCalloc(list_size,sizeof(t_tree *));

  return tlist;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

t_opt *Make_Optimiz()
{
  t_opt *s_opt;
  s_opt = (t_opt *)mCalloc(1,sizeof(t_opt));
  return s_opt;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Make_Custom_Model(t_mod *mod)
{
  if(!mod->r_mat)
    {
      PhyML_Printf("\n== Err. in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
  
  if(!mod->r_mat->rr->v)
    mod->r_mat->rr->v = (phydbl *)mCalloc(mod->ns*(mod->ns-1)/2,sizeof(phydbl));

  if(!mod->r_mat->rr_val->v)
    mod->r_mat->rr_val->v = (phydbl *)mCalloc(mod->ns*(mod->ns-1)/2,sizeof(phydbl));

  if(!mod->r_mat->rr_num->v)
    mod->r_mat->rr_num->v = (int *)mCalloc(mod->ns*(mod->ns-1)/2,sizeof(int *));

  if(!mod->r_mat->n_rr_per_cat->v)
    mod->r_mat->n_rr_per_cat->v = (int *)mCalloc(mod->ns*(mod->ns-1)/2,sizeof(int));
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

t_string *Make_String(int len)
{
  t_string *ts;

  ts = (t_string *)mCalloc(1,sizeof(t_string));
  ts->s = (char *)mCalloc(len,sizeof(char));

  return(ts);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

t_mod *Make_Model_Basic()
{
  t_mod *mod;

  mod = (t_mod *)mCalloc(1,sizeof(t_mod));

  mod->modelname = Make_String(T_MAX_NAME);
  Init_String(mod->modelname);

  mod->custom_mod_string = Make_String(T_MAX_NAME);
  Init_String(mod->custom_mod_string);

  mod->ras = Make_RAS_Basic();

  mod->kappa                  = (scalar_dbl *)mCalloc(1,sizeof(scalar_dbl));
  Init_Scalar_Dbl(mod->kappa);

  mod->lambda                 = (scalar_dbl *)mCalloc(1,sizeof(scalar_dbl));
  Init_Scalar_Dbl(mod->lambda);

  mod->br_len_mult      = (scalar_dbl *)mCalloc(1,sizeof(scalar_dbl));
  Init_Scalar_Dbl(mod->br_len_mult);

  mod->br_len_mult_unscaled      = (scalar_dbl *)mCalloc(1,sizeof(scalar_dbl));
  Init_Scalar_Dbl(mod->br_len_mult_unscaled);

  mod->mr                     = (scalar_dbl *)mCalloc(1,sizeof(scalar_dbl));
  Init_Scalar_Dbl(mod->mr);

  mod->e_frq_weight           = (scalar_dbl *)mCalloc(1,sizeof(scalar_dbl));
  Init_Scalar_Dbl(mod->e_frq_weight);

  mod->r_mat_weight           = (scalar_dbl *)mCalloc(1,sizeof(scalar_dbl));
  Init_Scalar_Dbl(mod->r_mat_weight);

  mod->aa_rate_mat_file       = Make_String(T_MAX_FILE);
  Init_String(mod->aa_rate_mat_file);

  return mod;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/*! Call only when the values of mod->ns & ras->n_catg is set to its final value */

void Make_Model_Complete(t_mod *mod)
{
  mod->Pij_rr = (vect_dbl *)mCalloc(1,sizeof(vect_dbl));
  Init_Vect_Dbl(0,mod->Pij_rr);

#if (defined(__AVX__) || defined(__AVX2__) || defined(__SSE__) || defined(__SSE2__) || defined(__SSE3__))
#ifndef WIN32
  if(posix_memalign((void **)&mod->Pij_rr->v,BYTE_ALIGN,(size_t)mod->ras->n_catg*mod->ns*mod->ns*sizeof(phydbl))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
#else
  mod->Pij_rr->v = _aligned_malloc(mod->ras->n_catg*mod->ns*mod->ns*sizeof(phydbl),BYTE_ALIGN);
#endif
#else
  mod->Pij_rr->v = (phydbl *)mCalloc(mod->ras->n_catg*mod->ns*mod->ns,sizeof(phydbl));
#endif

  mod->eigen     = (eigen *)Make_Eigen_Struct(mod->ns);
  
  // If r_mat (e_frq) are not NULL, then they have been created elsewhere and affected.
  if(!mod->r_mat)
    {
      mod->r_mat = (t_rmat *)Make_Rmat(mod->ns);
      Init_Rmat(mod->r_mat);
    }

  if(!mod->e_frq)
    {
      mod->e_frq = (t_efrq *)Make_Efrq(mod->ns);
      Init_Efrq(NULL,mod->e_frq);
    }

  Make_RAS_Complete(mod->ras);

  mod->e_frq->user_b_freq->len = mod->ns;

  if(mod->whichmodel < 0)
    {
      PhyML_Printf("\n== Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");
    }

  if(mod->whichmodel == CUSTOM)
    {
      Make_Custom_Model(mod);
      Translate_Custom_Mod_String(mod);
    }
  
  if(mod->io->datatype == NT && mod->whichmodel == GTR)
    {
      Make_Custom_Model(mod);
    }
  else if(mod->io->datatype == GENERIC)
    {
      Make_Custom_Model(mod);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

t_ras *Make_RAS_Basic()
{
  t_ras *ras;

  ras = (t_ras *)mCalloc(1,sizeof(t_ras));

  ras->gamma_r_proba          = (vect_dbl *)mCalloc(1,sizeof(vect_dbl));
  Init_Vect_Dbl(0,ras->gamma_r_proba);
  ras->gamma_r_proba->v = NULL;

  ras->gamma_r_proba_unscaled = (vect_dbl *)mCalloc(1,sizeof(vect_dbl));
  Init_Vect_Dbl(0,ras->gamma_r_proba_unscaled);
  ras->gamma_r_proba_unscaled->v = NULL;

  ras->gamma_rr               = (vect_dbl *)mCalloc(1,sizeof(vect_dbl));
  Init_Vect_Dbl(0,ras->gamma_rr);
  ras->gamma_rr->v = NULL;

  ras->gamma_rr_unscaled      = (vect_dbl *)mCalloc(1,sizeof(vect_dbl));
  Init_Vect_Dbl(0,ras->gamma_rr_unscaled);
  ras->gamma_rr_unscaled->v = NULL;

  ras->alpha             = (scalar_dbl *)mCalloc(1,sizeof(scalar_dbl));
  Init_Scalar_Dbl(ras->alpha);

  ras->pinvar            = (scalar_dbl *)mCalloc(1,sizeof(scalar_dbl));
  Init_Scalar_Dbl(ras->pinvar);

  ras->free_rate_mr           = (scalar_dbl *)mCalloc(1,sizeof(scalar_dbl));
  Init_Scalar_Dbl(ras->free_rate_mr);

  return(ras);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/*! Call only when the value of ras->n_catg is set to its final value */
void Make_RAS_Complete(t_ras *ras)
{
  if(!ras->gamma_r_proba->v)
    {
      ras->gamma_r_proba->v          = (phydbl *)mCalloc(ras->n_catg,sizeof(phydbl));
      ras->gamma_r_proba_unscaled->v = (phydbl *)mCalloc(ras->n_catg,sizeof(phydbl));
      ras->gamma_rr->v               = (phydbl *)mCalloc(ras->n_catg,sizeof(phydbl));
      ras->gamma_rr_unscaled->v      = (phydbl *)mCalloc(ras->n_catg,sizeof(phydbl));
      ras->skip_rate_cat             = (short int *)mCalloc(ras->n_catg,sizeof(short int));
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

t_efrq *Make_Efrq(int ns)
{
  t_efrq *e_frq;

  e_frq = (t_efrq *)mCalloc(1,sizeof(t_efrq));

  e_frq->pi               = (vect_dbl *)mCalloc(1,sizeof(vect_dbl));
  e_frq->pi->len          = ns;

#if (defined(__AVX__) || defined(__AVX2__) || defined(__SSE__) || defined(__SSE2__) || defined(__SSE3__))
#ifndef WIN32
  if(posix_memalign((void **)&e_frq->pi->v,BYTE_ALIGN,(size_t)ns*sizeof(phydbl))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
#else
  e_frq->pi->v = _aligned_malloc(ns * sizeof(phydbl),BYTE_ALIGN);
#endif
#else
  e_frq->pi->v = (phydbl *)mCalloc(ns,sizeof(phydbl));
#endif


  e_frq->pi_unscaled      = (vect_dbl *)mCalloc(1,sizeof(vect_dbl));
  e_frq->pi_unscaled->v   = (phydbl *)mCalloc(ns,sizeof(phydbl));
  e_frq->pi_unscaled->len = ns;

  e_frq->user_b_freq      = (vect_dbl *)mCalloc(1,sizeof(vect_dbl));
  Init_Vect_Dbl(0,e_frq->user_b_freq);
  e_frq->user_b_freq->v   = (phydbl *)mCalloc(T_MAX_OPTION,sizeof(phydbl));

  return e_frq;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

t_rmat *Make_Rmat(int ns)
{
  t_rmat *r_mat;

  r_mat                         = (t_rmat *)mCalloc(1,sizeof(t_rmat));

  r_mat->qmat                   = (vect_dbl *)mCalloc(1,sizeof(vect_dbl));
  Init_Vect_Dbl(0,r_mat->qmat);

  r_mat->qmat_buff              = (vect_dbl *)mCalloc(1,sizeof(vect_dbl));
  Init_Vect_Dbl(0,r_mat->qmat_buff);

  r_mat->rr                     = (vect_dbl *)mCalloc(1,sizeof(vect_dbl));
  Init_Vect_Dbl(0,r_mat->rr);

  r_mat->rr_val                 = (vect_dbl *)mCalloc(1,sizeof(vect_dbl));
  Init_Vect_Dbl(0,r_mat->rr_val);

  r_mat->rr_num                 = (vect_int *)mCalloc(1,sizeof(vect_int));
  Init_Vect_Int(0,r_mat->rr_num);

  r_mat->n_rr_per_cat           = (vect_int *)mCalloc(1,sizeof(vect_int));
  Init_Vect_Int(0,r_mat->n_rr_per_cat);

  r_mat->qmat->v                = (phydbl *)mCalloc(ns*ns,sizeof(phydbl));
  r_mat->qmat_buff->v           = (phydbl *)mCalloc(ns*ns,sizeof(phydbl));

  return(r_mat);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

option *Make_Input()
{
  int i;
  option* io                            = (option *)mCalloc(1,sizeof(option));

  io->in_align_file                     = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  io->in_tree_file                      = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  io->in_constraint_tree_file           = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  io->in_coord_file                     = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  io->out_file                          = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  io->out_tree_file                     = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  io->out_trees_file                    = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  io->out_ancestral_seq_file            = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  io->out_ancestral_tree_file           = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  io->out_boot_tree_file                = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  io->out_boot_stats_file               = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  io->out_stats_file                    = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  io->out_lk_file                       = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  io->weight_file                       = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  io->out_ps_file                       = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  io->out_summary_file                  = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  io->out_trace_file                    = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  io->out_json_trace_file               = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  io->nt_or_cd                          = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  io->run_id_string                     = (char *)mCalloc(T_MAX_OPTION,sizeof(char));
  io->clade_list_file                   = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  io->alphabet                          = (char **)mCalloc(T_MAX_ALPHABET,sizeof(char *));
  for(i=0;i<T_MAX_ALPHABET;i++) io->alphabet[i] = (char *)mCalloc(T_MAX_STATE,sizeof(char ));
  io->treelist                          = (t_treelist *)mCalloc(1,sizeof(t_treelist));
  io->mcmc                              = (t_mcmc *)MCMC_Make_MCMC_Struct();
  io->rates                             = (t_rate *)RATES_Make_Rate_Struct(-1);
  io->times                             = (t_time *)TIMES_Make_Time_Struct(-1);

  return io;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

t_mcmc *MCMC_Make_MCMC_Struct()
{
  t_mcmc *mcmc;

  mcmc               = (t_mcmc *)mCalloc(1,sizeof(t_mcmc));
  mcmc->out_filename = (char *)mCalloc(T_MAX_FILE,sizeof(char));

  return(mcmc);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

eigen *Make_Eigen_Struct(int ns)
{
  eigen *eig;

  eig              = (eigen *)mCalloc(1,sizeof(eigen));
  eig->size        = ns;
  eig->space       = (phydbl *)mCalloc(2*ns,sizeof(phydbl));
  eig->space_int   = (int *)mCalloc(2*ns,sizeof(int));
  eig->e_val       = (phydbl *)mCalloc(ns,sizeof(phydbl));
  eig->e_val_im    = (phydbl *)mCalloc(ns,sizeof(phydbl));
  eig->r_e_vect_im = (phydbl *)mCalloc(ns*ns,sizeof(phydbl));
  eig->q           = (phydbl *)mCalloc(ns*ns,sizeof(phydbl));

#if (defined(__AVX__) || defined(__AVX2__) || defined(__SSE__) || defined(__SSE2__) || defined(__SSE3__))
#ifndef WIN32
  if(posix_memalign((void **)&eig->r_e_vect,BYTE_ALIGN,ns*ns*sizeof(phydbl))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  if(posix_memalign((void **)&eig->l_e_vect,BYTE_ALIGN,ns*ns*sizeof(phydbl))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  if(posix_memalign((void **)&eig->dum,BYTE_ALIGN,ns*ns*sizeof(phydbl))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
#else
  eig->r_e_vect = _aligned_malloc(ns*ns*sizeof(phydbl),BYTE_ALIGN);
  eig->l_e_vect = _aligned_malloc(ns*ns*sizeof(phydbl),BYTE_ALIGN);
  eig->dum      = _aligned_malloc(ns*ns*sizeof(phydbl),BYTE_ALIGN);
#endif
#else
  eig->r_e_vect = (phydbl *)mCalloc(ns*ns,sizeof(phydbl));
  eig->l_e_vect = (phydbl *)mCalloc(ns*ns,sizeof(phydbl));
  eig->dum      = (phydbl *)mCalloc(ns*ns,sizeof(phydbl));
#endif
  
  Init_Eigen_Struct(eig);

  return eig;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Make_Short_L(t_tree *tree)
{
  if(!tree->short_l)
    tree->short_l = (phydbl *)mCalloc(tree->n_short_l,sizeof(phydbl));
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

xml_attr *XML_Make_Attribute(xml_attr *prev, char *attr_name, char *attr_value)
{
  xml_attr *new_attr;

  new_attr = (xml_attr *)mCalloc(1,sizeof(xml_attr));

  new_attr->prev = prev;
  new_attr->next = NULL;
  if(prev != NULL) prev->next = new_attr;

  new_attr->name = (char *)mCalloc(strlen(attr_name)+1,sizeof(char));
  strcpy(new_attr->name,attr_name);

  new_attr->value = (char *)mCalloc(strlen(attr_value)+1,sizeof(char));
  strcpy(new_attr->value,attr_value);

  return new_attr;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void XML_Make_Node_Id(xml_node *n, char *id)
{
  if(id) n->id = (char *)mCalloc(strlen(id)+1,sizeof(char));
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void XML_Make_Node_Value(xml_node *n, char *val)
{
  if(val) n->value = (char *)mCalloc(strlen(val)+1,sizeof(char));
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


xml_node *XML_Make_Node(char *name)
{
  xml_node *new_node = (xml_node *)mCalloc(1,sizeof(xml_node));
  if(name) new_node->name = (char *)mCalloc(strlen(name)+1,sizeof(char));
  new_node->ds = (t_ds *)mCalloc(1,sizeof(t_ds));
  return new_node;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Make_Best_Spr(t_tree *tree)
{
  tree->best_spr = Make_One_Spr(tree);
  Init_One_Spr(tree->best_spr);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Make_Spr(t_tree *tree)
{
  Make_Spr_List_One_Edge(tree);
  Make_Spr_List_All_Edge(tree);
  Make_Best_Spr(tree);
  if(tree->is_mixt_tree == YES) MIXT_Repeat_Task(Make_Spr,tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Make_Spr_List_One_Edge(t_tree *tree)
{
  int i;

  tree->size_spr_list_one_edge = 2*tree->n_otu-3;

  tree->spr_list_one_edge = (t_spr **)mCalloc(2*tree->n_otu-2,sizeof(t_spr *));

  for(i=0;i<2*tree->n_otu-2;++i)
    {
      tree->spr_list_one_edge[i] = Make_One_Spr(tree);
      Init_One_Spr(tree->spr_list_one_edge[i]);
    }

  tree->perform_spr_right_away = NO;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Make_Spr_List_All_Edge(t_tree *tree)
{
  int i;

  tree->size_spr_list_all_edge = 2*tree->n_otu-3;

  tree->spr_list_all_edge = (t_spr **)mCalloc(2*tree->n_otu-2,sizeof(t_spr *));

  for(i=0;i<2*tree->n_otu-2;++i)
    {
      tree->spr_list_all_edge[i] = Make_One_Spr(tree);
      Init_One_Spr(tree->spr_list_all_edge[i]);
    }

  tree->perform_spr_right_away = NO;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

t_spr *Make_One_Spr(t_tree *tree)
{
  t_spr *a_spr;
  a_spr         = (t_spr *)mCalloc(1,sizeof(t_spr));
  a_spr->path   = (t_node **)mCalloc(tree->n_otu,sizeof(t_node *));
  return a_spr;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

t_time *TIMES_Make_Time_Struct(int n_otu)
{
  t_time *times;

  times = (t_time *)mCalloc(1,sizeof(t_time));

  if(n_otu > 0)
    {
      times->calib_prob           = (phydbl *)mCalloc(2*n_otu-1,sizeof(phydbl));
      times->t_prior_min_ori      = (phydbl *)mCalloc(2*n_otu-1,sizeof(phydbl));
      times->t_prior_max_ori      = (phydbl *)mCalloc(2*n_otu-1,sizeof(phydbl));
      times->times_partial_proba  = (phydbl *)mCalloc(n_otu*n_otu,sizeof(phydbl));
      times->numb_calib_chosen    = (int *)mCalloc(n_otu*n_otu,sizeof(phydbl));
      times->a_cal                = (t_cal **)mCalloc(10*n_otu,sizeof(t_cal *));
      
      times->nd_t                 = (phydbl *)mCalloc(2*n_otu-1,sizeof(phydbl));
      times->buff_t               = (phydbl *)mCalloc(2*n_otu-1,sizeof(phydbl));
      times->true_t               = (phydbl *)mCalloc(2*n_otu-1,sizeof(phydbl));
      times->t_mean               = (phydbl *)mCalloc(2*n_otu-1,sizeof(phydbl));
      times->t_prior              = (phydbl *)mCalloc(2*n_otu-1,sizeof(phydbl));
      times->t_prior_min          = (phydbl *)mCalloc(2*n_otu-1,sizeof(phydbl));
      times->t_prior_max          = (phydbl *)mCalloc(2*n_otu-1,sizeof(phydbl));
      times->t_floor              = (phydbl *)mCalloc(2*n_otu-1,sizeof(phydbl));
      times->t_rank               = (int *)mCalloc(2*n_otu-1,sizeof(int));
      times->t_has_prior          = (short int *)mCalloc(2*n_otu-1,sizeof(short int));
      
      times->mean_t               = (phydbl *)mCalloc(2*n_otu-1,sizeof(phydbl));
      
      times->n_jps                = (int    *)mCalloc(2*n_otu-1,sizeof(int));
      times->t_jps                = (int    *)mCalloc(2*n_otu-2,sizeof(int));
      
      times->time_slice_lims      = (phydbl *)mCalloc(2*n_otu-1,sizeof(phydbl));
      times->n_time_slice_spans   = (int *)mCalloc(2*n_otu-1,sizeof(int));
      times->curr_slice           = (int *)mCalloc(2*n_otu-1,sizeof(int));
      times->has_survived         = (int    *)mCalloc(2*n_otu-1,sizeof(int));
    }
  return(times);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

t_rate *RATES_Make_Rate_Struct(int n_otu)
{
  t_rate *rates;

  rates = (t_rate  *)mCalloc(1,sizeof(t_rate));
  rates->is_allocated = NO;  

  if(n_otu > 0)
    {
      rates->is_allocated         = YES;
      rates->nd_r                 = (phydbl *)mCalloc(2*n_otu-1,sizeof(phydbl));
      rates->br_r                 = (phydbl *)mCalloc(2*n_otu-1,sizeof(phydbl));
      rates->buff_br_r            = (phydbl *)mCalloc(2*n_otu-1,sizeof(phydbl));
      rates->buff_nd_r            = (phydbl *)mCalloc(2*n_otu-1,sizeof(phydbl));
      rates->true_r               = (phydbl *)mCalloc(2*n_otu-1,sizeof(phydbl));
      rates->dens                 = (phydbl *)mCalloc(2*n_otu-2,sizeof(phydbl));
      rates->triplet              = (phydbl *)mCalloc(2*n_otu-1,sizeof(phydbl));
      rates->cov_l                = (phydbl *)mCalloc((2*n_otu-2)*(2*n_otu-2),sizeof(phydbl));
      rates->invcov               = (phydbl *)mCalloc((2*n_otu-2)*(2*n_otu-2),sizeof(phydbl));
      rates->mean_l               = (phydbl *)mCalloc(2*n_otu-2,sizeof(phydbl));
      rates->ml_l                 = (phydbl *)mCalloc(2*n_otu-2,sizeof(phydbl));
      rates->cur_l                = (phydbl *)mCalloc(2*n_otu-2,sizeof(phydbl));
      rates->u_ml_l               = (phydbl *)mCalloc(2*n_otu-1,sizeof(phydbl));
      rates->u_cur_l              = (phydbl *)mCalloc(2*n_otu-1,sizeof(phydbl));
      rates->cov_r                = (phydbl *)mCalloc((2*n_otu-2)*(2*n_otu-2),sizeof(phydbl));
      rates->cond_var             = (phydbl *)mCalloc(2*n_otu-2,sizeof(phydbl));
      rates->mean_r               = (phydbl *)mCalloc(2*n_otu-2,sizeof(phydbl));
      rates->lca                  = (t_node **)mCalloc((2*n_otu-1)*(2*n_otu-1),sizeof(t_node *));
      rates->reg_coeff            = (phydbl *)mCalloc((2*n_otu-2)*(2*n_otu-2),sizeof(phydbl));
      rates->trip_reg_coeff       = (phydbl *)mCalloc((2*n_otu-2)*(6*n_otu-9),sizeof(phydbl));
      rates->trip_cond_cov        = (phydbl *)mCalloc((2*n_otu-2)*9,sizeof(phydbl));
      rates->_2n_vect1            = (phydbl *)mCalloc(2*n_otu,sizeof(phydbl));
      rates->_2n_vect2            = (phydbl *)mCalloc(2*n_otu,sizeof(phydbl));
      rates->_2n_vect3            = (phydbl *)mCalloc(2*n_otu,sizeof(phydbl));
      rates->_2n_vect4            = (phydbl *)mCalloc(2*n_otu,sizeof(phydbl));
      rates->_2n_vect5            = (short int *)mCalloc(2*n_otu,sizeof(short int));
      rates->_2n2n_vect1          = (phydbl *)mCalloc(4*n_otu*n_otu,sizeof(phydbl));
      rates->_2n2n_vect2          = (phydbl *)mCalloc(4*n_otu*n_otu,sizeof(phydbl));
      rates->br_do_updt           = (short int *)mCalloc(2*n_otu-1,sizeof(short int));
      rates->cur_gamma_prior_mean = (phydbl *)mCalloc(2*n_otu-1,sizeof(phydbl));
      rates->cur_gamma_prior_var  = (phydbl *)mCalloc(2*n_otu-1,sizeof(phydbl));
      rates->n_tips_below         = (int *)mCalloc(2*n_otu-1,sizeof(int));
      rates->model_name           = (char *)mCalloc(T_MAX_NAME,sizeof(char));
    }

  return rates;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

t_cal *Make_Calibration()
{
  t_cal *calib;
  calib = (t_cal *)mCalloc(1, sizeof(t_cal));
  return(calib);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

t_clad *Make_Clade()
{
  t_clad *clade;
  clade = (t_clad *)mCalloc(1, sizeof(t_clad));
  return(clade);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


t_node **Make_Target_Tip(int n)
{
  t_node **this;
  this = (t_node **)mCalloc(n,sizeof(t_node *));
  return(this);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Make_All_Calibration(t_tree *tree)
{
  int i;
  t_cal **all_cal;

  assert(tree->rates);

  all_cal = (t_cal **)mCalloc(2*tree->n_otu-1,sizeof(t_cal *));

  For(i,2*tree->n_otu-1) all_cal[i] = Make_Calibration();

  tree->times->a_cal = all_cal;

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Make_Rmat_Weight(t_tree *mixt_tree)
{
  t_tree *tree, *buff_tree;
  scalar_dbl *curr_weight;

  tree = mixt_tree;
  do
    {
      if(tree->is_mixt_tree == YES) tree = tree->next;

      buff_tree = mixt_tree->next;
      do
        {
          if(buff_tree->mod->r_mat_weight == tree->mod->r_mat_weight) break;
          buff_tree = buff_tree->next;
        }
      while(buff_tree != tree);

      if(buff_tree == tree) Free(tree->mod->r_mat_weight);

      tree = tree->next;
    }
  while(tree);


  tree = mixt_tree;
  do
    {
      if(tree->is_mixt_tree == YES) tree = tree->next;
      tree->mod->r_mat_weight = NULL;
      tree = tree->next;
    }
  while(tree);

  
  tree = mixt_tree->next;
  tree->mod->r_mat_weight = (scalar_dbl *)mCalloc(1,sizeof(scalar_dbl));
  Init_Scalar_Dbl(tree->mod->r_mat_weight);
  tree->mod->r_mat_weight->v = 1.0;
  curr_weight = tree->mod->r_mat_weight;


  buff_tree = tree = mixt_tree;
  do // For each mixt_tree
    {
      if(tree->is_mixt_tree == YES) 
        {
          tree = tree->next;
        }

      buff_tree = mixt_tree->next;
      do
        {
          if(buff_tree->mod->r_mat == tree->mod->r_mat)
            {
              tree->mod->r_mat_weight = buff_tree->mod->r_mat_weight;
              break;
            }
          buff_tree = buff_tree->next;
        }
      while(buff_tree != tree);
      
      if(!tree->mod->r_mat_weight)
        {
          tree->mod->r_mat_weight = (scalar_dbl *)mCalloc(1,sizeof(scalar_dbl));
          Init_Scalar_Dbl(tree->mod->r_mat_weight);
          tree->mod->r_mat_weight->v = 1.0;
          curr_weight->next = tree->mod->r_mat_weight;
          tree->mod->r_mat_weight->prev = curr_weight;
          curr_weight = tree->mod->r_mat_weight;
        }
      
      tree = tree->next;
      
    }
  while(tree);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Make_Efrq_Weight(t_tree *mixt_tree)
{
  t_tree *tree, *buff_tree;
  scalar_dbl *curr_weight;


  tree = mixt_tree;
  do
    {
      if(tree->is_mixt_tree == YES) tree = tree->next;

      buff_tree = mixt_tree->next;
      do
        {
          if(buff_tree->mod->e_frq_weight == tree->mod->e_frq_weight) break;
          buff_tree = buff_tree->next;
        }
      while(buff_tree != tree);

      if(buff_tree == tree)
            {
          Free(tree->mod->e_frq_weight);
        }

      tree = tree->next;
    }
  while(tree);


  tree = mixt_tree;
  do
    {
      if(tree->is_mixt_tree == YES) tree = tree->next;
      tree->mod->e_frq_weight = NULL;
      tree = tree->next;
    }
  while(tree);


  tree = mixt_tree->next;
  tree->mod->e_frq_weight = (scalar_dbl *)mCalloc(1,sizeof(scalar_dbl));
  Init_Scalar_Dbl(tree->mod->e_frq_weight);
  tree->mod->e_frq_weight->v = 1.0;
  curr_weight = tree->mod->e_frq_weight;


  buff_tree = tree = mixt_tree;
  do // For each mixt_tree
    {
      if(tree->is_mixt_tree == YES) tree = tree->next;
      
      buff_tree = mixt_tree->next;
      do
        {
          if(buff_tree->mod->e_frq == tree->mod->e_frq)
            {
              tree->mod->e_frq_weight = buff_tree->mod->e_frq_weight;
              break;
            }
          buff_tree = buff_tree->next;
        }
      while(buff_tree != tree);
      
      if(!tree->mod->e_frq_weight)
        {
          tree->mod->e_frq_weight = (scalar_dbl *)mCalloc(1,sizeof(scalar_dbl));
          Init_Scalar_Dbl(tree->mod->e_frq_weight);
          tree->mod->e_frq_weight->v = 1.0;
          curr_weight->next = tree->mod->e_frq_weight;
          tree->mod->e_frq_weight->prev = curr_weight;
          curr_weight = tree->mod->e_frq_weight;
        }
      
      tree = tree->next;
      
    }
  while(tree);
  }

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

t_geo *GEO_Make_Geo_Basic()
{
  t_geo *t;
  t = (t_geo *)mCalloc(1,sizeof(t_geo));
  return(t);  
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void GEO_Make_Geo_Complete(int ldscape_sz, int n_dim, int n_tax, t_geo *t)
{
  int i;

  // F matrix
  t->f_mat = (phydbl *)mCalloc(ldscape_sz*ldscape_sz,sizeof(phydbl));

  // R matrix
  t->r_mat = (phydbl *)mCalloc(ldscape_sz*ldscape_sz,sizeof(phydbl));

  // Occupation vectors: one vector for each node
  t->occup = (int *)mCalloc((2*n_tax-1)*ldscape_sz,sizeof(int));

  // Lineage locations
  t->idx_loc = (int *)mCalloc((int)(2*n_tax-1),sizeof(int));

  // Sorted node heights
  t->sorted_nd = (t_node **)mCalloc((int)(2*n_tax-1),sizeof(t_node *));

  // Covariance matrix
  t->cov = (phydbl *)mCalloc((int)(n_dim*n_dim),sizeof(phydbl));

  // gives the location occupied beneath each node in the tree
  t->idx_loc_beneath = (int *)mCalloc((int)(2*n_tax-1)*ldscape_sz,sizeof(int));

  // Locations
  t->coord_loc = (t_geo_coord **)mCalloc(ldscape_sz,sizeof(t_geo_coord *));
  for(i=0;i<ldscape_sz;i++) t->coord_loc[i] = GEO_Make_Geo_Coord(n_dim);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

t_geo_coord *GEO_Make_Geo_Coord(int dim)
{
  t_geo_coord *t;
  t = (t_geo_coord *)mCalloc(1,sizeof(t_geo_coord));
  t->lonlat = (phydbl *)mCalloc(dim,sizeof(phydbl));
  t->id = (char *)mCalloc(T_MAX_ID_COORD,sizeof(char));

  t->cpy = (t_geo_coord *)mCalloc(1,sizeof(t_geo_coord));
  t->cpy->lonlat = (phydbl *)mCalloc(dim,sizeof(phydbl));
  t->cpy->id = (char *)mCalloc(T_MAX_ID_COORD,sizeof(char));

  return(t);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

t_phyrex_mod *PHYREX_Make_Migrep_Model(int n_otu,int dim)
{
  t_phyrex_mod *t;
  t = (t_phyrex_mod *)mCalloc(1,sizeof(t_phyrex_mod));
  t->lim_up = GEO_Make_Geo_Coord(dim);
  t->lim_do = GEO_Make_Geo_Coord(dim);
  t->sigsq = (phydbl *)mCalloc(dim,sizeof(phydbl));
  t->sigsq_scale = (phydbl *)mCalloc(2*n_otu-1,sizeof(phydbl));
  return(t);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

t_dsk *PHYREX_Make_Disk_Event(int n_dim, int n_otu)
{
  t_dsk *t;

  t         = (t_dsk *)mCalloc(1,sizeof(t_dsk));  
  t->centr  = GEO_Make_Geo_Coord(n_dim);
  t->id     = (char *)mCalloc(T_MAX_ID_DISK,sizeof(char));
  t->ldsk_a = (t_ldsk **)mCalloc(n_otu,sizeof(t_ldsk *));

  return(t);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

t_ldsk *PHYREX_Make_Lindisk_Node(int n_dim)
{
  t_ldsk *t;
  t = (t_ldsk *)mCalloc(1,sizeof(t_ldsk));
  t->coord     = GEO_Make_Geo_Coord(n_dim);
  t->cpy_coord = GEO_Make_Geo_Coord(n_dim);
  return(t);
}
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void PHYREX_Make_Lindisk_Next(t_ldsk *t)
{
  
  if(t->n_next == 0)
    t->next = (t_ldsk **)mCalloc(NEXT_BLOCK_SIZE,sizeof(t_ldsk *));
  else if(!(t->n_next%NEXT_BLOCK_SIZE))
    t->next = (t_ldsk **)mRealloc(t->next,t->n_next+NEXT_BLOCK_SIZE,sizeof(t_ldsk *));

  t->n_next++;
  /* PhyML_Printf("\n. make next for ldsk %s n_next set to %d [%p] %d",t->coord->id,t->n_next,t->next,NEXT_BLOCK_SIZE); */
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

t_poly *Make_Poly(int n)
{
  t_poly *p;
  int i;
  p = (t_poly *)mCalloc(1,sizeof(t_poly));
  p->poly_vert = (t_geo_coord **)mCalloc(n,sizeof(t_geo_coord *));
  for(i=0;i<n;i++) p->poly_vert[i] = GEO_Make_Geo_Coord(2);
  return(p);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

t_sarea *Make_Sarea(int n_poly)
{
  t_sarea *s;
  s = (t_sarea *)mCalloc(1,sizeof(t_sarea ));
  s->a_poly = (t_poly **)mCalloc(n_poly,sizeof(t_poly *));
  return(s);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

t_ll *Make_Linked_List()
{
  t_ll *list;

  // Create list if non-existing and add
  list = (t_ll *)mCalloc(1,sizeof(t_ll));

  return list;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* Matrices used in transfer bootstrap computation (tbe.c) */
void Alloc_TBE_Matrices(int n_otu, short unsigned*** i_matrix, short unsigned*** c_matrix,short unsigned*** hamming, short unsigned** min_dist, short unsigned**  min_dist_edge, int** cluster_sizes){
  int i;
  int nb_edges = 2*n_otu-3;
  (*min_dist) = (short unsigned*) malloc(nb_edges*sizeof(short unsigned)); /* array of min Hamming distances */
  (*min_dist_edge) = (short unsigned*) malloc(nb_edges*sizeof(short unsigned)); /* array of edge ids corresponding to min Hamming distances */
  (*cluster_sizes) = (int*) malloc(nb_edges*sizeof(int)); /* array of sizes of clusters associated to each branch (in the post order traversal) */
  (*c_matrix) = (short unsigned**) malloc(nb_edges*sizeof(short unsigned*)); /* matrix of cardinals of complements */
  (*i_matrix) = (short unsigned**) malloc(nb_edges*sizeof(short unsigned*)); /* matrix of cardinals of intersections */
  (*hamming) = (short unsigned**) malloc(nb_edges*sizeof(short unsigned*)); /* matrix of Hamming distances */
  for (i=0; i<nb_edges; i++){
    (*c_matrix)[i] = (short unsigned*) malloc(nb_edges*sizeof(short unsigned));
    (*i_matrix)[i] = (short unsigned*) malloc(nb_edges*sizeof(short unsigned));
    (*hamming)[i] = (short unsigned*) malloc(nb_edges*sizeof(short unsigned));
    (*min_dist)[i] = n_otu; /* initialization to the nb of taxa */
    (*cluster_sizes)[i] = 0;
  }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Make_MutMap(t_tree *tree)
{
  // (# of edges) X (# of sites) X (# mutation types: A<->C A<->G A<->T C<->G C<->T G<->T)
  do
    {
      tree->mutmap = (int *)mCalloc((2*tree->n_otu-3)*(tree->n_pattern)*6,sizeof(int));
      tree = tree->next_mixt;
    }
  while(tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

t_label *Make_Label()
{
  t_label *lab;

  lab = (t_label *)mCalloc(1,sizeof(t_label));
  lab->key = (char *)mCalloc(T_MAX_KEY,sizeof(char));
  lab->val = (char *)mCalloc(T_MAX_VAL,sizeof(char));
  lab->next = NULL;

  return(lab);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Make_Contrasts(t_tree *tree)
{
  tree->ctrst = (t_ctrst *)mCalloc(1,sizeof(t_ctrst));
  tree->ctrst->x = (phydbl *)mCalloc(2*tree->n_otu,sizeof(phydbl));
  tree->ctrst->tprime = (phydbl *)mCalloc(2*tree->n_otu,sizeof(phydbl));
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



