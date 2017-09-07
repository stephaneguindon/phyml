/*

PHYML :  a program that  computes maximum likelihood  phylogenies from
DNA or AA homologous sequences

Copyright (C) Stephane Guindon. Oct 2003 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/

#include "free.h"


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Free_All_Nodes_Light(t_tree *mixt_tree)
{
  int i;
  t_tree *tree;

  tree = mixt_tree;
  do
    {
      For(i,2*tree->n_otu-1) Free_Node(tree->a_nodes[i]);
      Free(tree->a_nodes);
      tree = tree->next;
    }
  while(tree);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Free_All_Edges_Light(t_tree *mixt_tree)
{
  int i;
  t_tree *tree;

  tree = mixt_tree;

  For(i,2*tree->n_otu-1)
    {
      Free_Scalar_Dbl(tree->a_edges[i]->l);
      Free_Scalar_Dbl(tree->a_edges[i]->l_old);
      Free_Scalar_Dbl(tree->a_edges[i]->l_var);
      Free_Scalar_Dbl(tree->a_edges[i]->l_var_old);
    }

  do
    {
      For(i,2*tree->n_otu-1) Free_Edge(tree->a_edges[i]);
      Free(tree->a_edges);
      tree = tree->next;
    }
  while(tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Free_Edge_Labels(t_edge *b)
{
  int i;

  if(b->labels)
    {
      For(i,b->n_labels-(b->n_labels%BLOCK_LABELS)+BLOCK_LABELS) Free(b->labels[i]);
      Free(b->labels);
      b->labels = NULL;
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Free_Edge(t_edge *b)
{
  Free_Edge_Labels(b);
  Free_Edge_Core(b);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Free_Edge_Core(t_edge *b)
{
  Free(b);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Free_Node(t_node *n)
{
  Free(n->b);
  Free(n->v);
  Free(n->l);
  Free(n->score);
  Free(n->s_ingrp);
  Free(n->s_outgrp);
  Free(n->cal);

  if(n->c_seq_anc != NULL) 
    {
      Free(n->c_seq_anc->state);
      Free(n->c_seq_anc);
    }
  if(n->ori_name) { Free(n->ori_name); n->ori_name = NULL; }

  /* if(n->name)     { Free(n->name);     n->name     = NULL; }  */
  /* Don't do that: see Copy_Tax_Names_To_Tip_Labels
     tree->a_nodes[i]->ori_name = tree->a_nodes[i]->name; */

  Free(n);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Free_Mat(matrix *mat)
{
  int i;

  for(i=0;i<mat->n_otu;i++)
    {
      Free(mat->P[i]);
      Free(mat->Q[i]);
      Free(mat->dist[i]);
      Free(mat->name[i]);
    }

  Free(mat->P);
  Free(mat->Q);
  Free(mat->dist);
  Free(mat->name);
  Free(mat->tip_node);

  Free(mat->on_off);
  Free(mat);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Free_Partial_Lk(phydbl *p_lk, int len, int n_catg)
{
  Free(p_lk);

/*   int i,j; */
/*   for(i=0;i<len;i++) */
/*     { */
/*       for(j=0;j<n_catg;j++) Free((*p_lk)[i][j]); */
/*       Free((*p_lk)[i]); */
/*     } */
/*   Free((*p_lk)); */
/*   (*p_lk) = NULL; */
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Free_Tree(t_tree *mixt_tree)
{
  t_tree *tree;
  t_tree *next;

  tree = mixt_tree;
  do
    {
      if(tree->mat) Free_Mat(tree->mat);
      Free(tree->t_dir);
      if(tree->short_l) Free(tree->short_l);
      if(tree->mutmap)  Free(tree->mutmap);
      Free_Bip(tree);
      Free(tree->curr_path);
      tree = tree->next;
    }
  while(tree);

  Free_All_Edges_Light(mixt_tree);
  Free_All_Nodes_Light(mixt_tree);

  tree = mixt_tree;
  next = mixt_tree->next;
  do
    {
      Free(tree);
      tree = next;
      if(!tree) break;
      next = next->next;
    }
  while(tree);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Free_Bip(t_tree *tree)
{
  int i,j;

  if(tree->has_bip)
    {
      For(i,2*tree->n_otu-2)
    {
      Free(tree->a_nodes[i]->bip_size);
      for(j=0;j<3;j++) Free(tree->a_nodes[i]->bip_node[j]);
      Free(tree->a_nodes[i]->bip_node);
    }
    }
  tree->has_bip = NO;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Free_Calign(calign *data)
{
  int i;

  if(data->io_wght) Free_Scalar_Dbl(data->io_wght);
  Free(data->invar);
  Free(data->wght);
  Free(data->ambigu);
  Free(data->b_frq);
  Free(data->sitepatt);
  for(i=0;i<data->n_otu;i++)
    {
      Free(data->c_seq[i]->name);
      if(data->c_seq[i]->state)
        {
          Free(data->c_seq[i]->state);
          Free(data->c_seq[i]->d_state);
          if(data->c_seq[i]->is_ambigu) Free(data->c_seq[i]->is_ambigu);
        }
      Free(data->c_seq[i]);
    }
  Free(data->c_seq);
  Free(data);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Free_Seq(align **d, int n_otu)
{
  int i;
  for(i=0;i<n_otu;i++)
    {
      Free(d[i]->name);
      Free(d[i]->state);
      Free(d[i]->d_state);
      if(d[i]->is_ambigu) Free(d[i]->is_ambigu);
      Free(d[i]);
    }
  Free(d);
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Free_All(align **d, calign *cdata, t_tree *tree)
{
  Free_Calign(cdata);
  Free_Seq(d,tree->n_otu);
  Free_Tree(tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Free_SubTree(t_edge *b_fcus, t_node *a, t_node *d, t_tree *tree)
{
  int i;

  if(d->tax) return;
  else
    {
      for(i=0;i<3;i++)
    {
      if(d->v[i] != a)
        {
          Free_SubTree(d->b[i],d,d->v[i],tree);
          Free_Edge(d->b[i]);
          Free_Node(d->v[i]);
        }
    }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Free_Tree_Ins_Tar(t_tree *tree)
{
  return;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Free_Tree_Pars(t_tree *mixt_tree)
{
  int i;
  t_tree *tree;

  tree = mixt_tree;
  do
    {
      Free(tree->step_mat);
      Free(tree->site_pars);
      
      For(i,2*tree->n_otu-3) 
        {
          Free_Edge_Pars(tree->a_edges[i]);
        }

      if(tree->n_root)
        {
          Free_Edge_Pars_Left(tree->n_root->b[1]);
          Free_Edge_Pars_Left(tree->n_root->b[2]);
        }
      else
        {
          Free_Edge_Pars(tree->a_edges[2*tree->n_otu-3]);
          Free_Edge_Pars(tree->a_edges[2*tree->n_otu-2]);
        }

      tree = tree->next;
    }
  while(tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Free_Edge_Pars_Left(t_edge *b)
{
  if(b->pars_l)          Free(b->pars_l);
  if(b->ui_l)            Free(b->ui_l);
  if(b->p_pars_l)        Free(b->p_pars_l);
  if(b->n_diff_states_l) Free(b->n_diff_states_l);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Free_Edge_Pars_Rght(t_edge *b)
{
  if(b->pars_r)   Free(b->pars_r);
  if(b->ui_r)     Free(b->ui_r);
  if(b->p_pars_r) Free(b->p_pars_r);
  if(b->n_diff_states_r) Free(b->n_diff_states_r);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Free_Edge_Pars(t_edge *b)
{
  Free_Edge_Pars_Left(b);
  Free_Edge_Pars_Rght(b);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Free_Tree_Lk(t_tree *mixt_tree)
{
  int i;
  t_tree *tree;

  tree = mixt_tree;
  do
    {
      Free(tree->c_lnL_sorted);
      Free(tree->cur_site_lk);
      Free(tree->old_site_lk);
      Free(tree->site_lk_cat);
      Free(tree->fact_sum_scale);
      Free(tree->eigen_lr_left);
      Free(tree->eigen_lr_rght);
      Free(tree->dot_prod);
      Free(tree->expl);
      
      for(i=0;i<3;i++) Free(tree->log_lks_aLRT[i]);
      Free(tree->log_lks_aLRT);

      Free(tree->unscaled_site_lk_cat);
      
      For(i,2*tree->n_otu-1) Free_NNI(tree->a_edges[i]->nni);
  
      if(tree->is_mixt_tree == NO)
        {
          For(i,2*tree->n_otu-3) Free_Edge_Lk(tree->a_edges[i]);
          For(i,2*tree->n_otu-3) Free_Edge_Loc(tree->a_edges[i]);
          
          if(tree->n_root != NULL)
            {
              Free(tree->n_root->b[1]->Pij_rr);
              Free(tree->n_root->b[2]->Pij_rr);
              Free(tree->n_root->b[1]->tPij_rr);
              Free(tree->n_root->b[2]->tPij_rr);
              Free_Edge_Lk_Left(tree->n_root->b[1]);
              Free_Edge_Lk_Left(tree->n_root->b[2]);
              Free_Edge_Loc_Left(tree->n_root->b[1]);
              Free_Edge_Loc_Left(tree->n_root->b[2]);
            }
          else
            {
              Free_Edge_Lk(tree->a_edges[2*tree->n_otu-3]);
              Free_Edge_Lk(tree->a_edges[2*tree->n_otu-2]);
              Free_Edge_Loc(tree->a_edges[2*tree->n_otu-3]);
              Free_Edge_Loc(tree->a_edges[2*tree->n_otu-2]);
            }
        }
      
      tree = tree->next;

    }
  while(tree);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Free_Node_Lk(t_node *n)
{
/*   Free(n->n_ex_nodes); */
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Free_Edge_Lk_Rght(t_edge *b)
{
  Free(b->div_post_pred_rght);

  if(b->p_lk_rght)
    {
      Free(b->p_lk_rght);
      if(b->sum_scale_rght) Free(b->sum_scale_rght);
    }

  if(b->p_lk_tip_r)         Free(b->p_lk_tip_r);
  if(b->sum_scale_rght_cat) Free(b->sum_scale_rght_cat);
  if(b->patt_id_rght)       Free(b->patt_id_rght);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Free_Edge_Lk_Left(t_edge *b)
{

  Free(b->div_post_pred_left);

  if(b->p_lk_left)
    {
      Free(b->p_lk_left);
      if(b->sum_scale_left) Free(b->sum_scale_left);
    }

  if(b->p_lk_tip_l)         Free(b->p_lk_tip_l);
  if(b->sum_scale_left_cat) Free(b->sum_scale_left_cat);
  if(b->patt_id_left)       Free(b->patt_id_left);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Free_Edge_Lk(t_edge *b)
{
  Free(b->tPij_rr);
  Free(b->Pij_rr);
  Free_Edge_Lk_Left(b);
  Free_Edge_Lk_Rght(b);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Free_Edge_Loc_Rght(t_edge *b)
{
  if(b->p_lk_loc_rght) Free(b->p_lk_loc_rght);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Free_Edge_Loc_Left(t_edge *b)
{
  if(b->p_lk_loc_left) Free(b->p_lk_loc_left);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Free_Edge_Loc(t_edge *b)
{
  Free_Edge_Loc_Left(b);
  Free_Edge_Loc_Rght(b);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Free_Model_Complete(t_mod *mixt_mod)
{
  Free_Eigen(mixt_mod->eigen);
  Free_Rmat(mixt_mod->r_mat);
  Free_Efrq(mixt_mod->e_frq);
  Free_Vect_Dbl(mixt_mod->Pij_rr);
  mixt_mod->r_mat = NULL;
  mixt_mod->e_frq = NULL;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Free_Rmat_Weights(t_mod *mixt_mod)
{
  t_mod *mod;

  mod = mixt_mod;

  do
    {
      Free(mod->r_mat_weight);
      mod = mod->next_mixt;
    }
  while(mod);

  if(mixt_mod->next) Free_Scalar_Dbl(mixt_mod->next->r_mat_weight);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Free_Efrq_Weights(t_mod *mixt_mod)
{
  t_mod *mod;

  mod = mixt_mod;

  do
    {
      Free(mod->e_frq_weight);
      mod = mod->next_mixt;
    }
  while(mod);

  if(mixt_mod->next) Free_Scalar_Dbl(mixt_mod->next->e_frq_weight);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Free_Model_Basic(t_mod *mixt_mod)
{
  t_mod *mod;

  Free_RAS(mixt_mod->ras);
  Free_Scalar_Dbl(mixt_mod->mr);
  Free_Scalar_Dbl(mixt_mod->kappa);
  Free_Scalar_Dbl(mixt_mod->lambda);
  Free_Scalar_Dbl(mixt_mod->br_len_mult);
  Free_Scalar_Dbl(mixt_mod->br_len_mult_unscaled);

  Free_Rmat_Weights(mixt_mod);
  Free_Efrq_Weights(mixt_mod);

  Free_String(mixt_mod->modelname);
  Free_String(mixt_mod->custom_mod_string);
  Free_String(mixt_mod->aa_rate_mat_file);

  mod = mixt_mod;
  do
    {
      if(mod->next)
        {
          mod = mod->next;
          Free(mod->prev);
        }
      else
        {
          Free(mod);
          break;
        }
    }
  while(mod);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Free_Vect_Dbl(vect_dbl *v)
{
  vect_dbl *next;

  assert(v);

  next = v->next;
  do
    {
      Free(v->v);
      Free(v);

      v = next;
      if(v) next = v->next;
    }
  while(v);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Free_Vect_Int(vect_int *v)
{
  vect_int *next;

  assert(v);

  next = v->next;
  do
    {
      Free(v->v);
      Free(v);

      v = next;
      if(v) next = v->next;
    }
  while(v);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Free_Scalar_Dbl(scalar_dbl *v)
{
  scalar_dbl *next;
  
  /* assert(v); */

  next = v->next;
  do
    {
      Free(v);
      v = next;
      if(v) next = v->next;
    }
  while(v);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Free_Scalar_Int(scalar_int *v)
{
  scalar_int *next;
  
  assert(v);

  next = v->next;
  do
    {
      Free(v);
      v = next;
      if(v) next = v->next;
    }
  while(v);
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Free_Linked_List(t_ll *t)
{
  t_ll *next,*ll;

  if(t == NULL) return;

  ll = t->head;
  next = ll->next;
  do
    {
      /* t_node *n = ll->v; */
      /* printf("\n. free node %d",n?n->num:-1); */
      /* printf(" ll: %p",ll); */
      /* printf(" hd: %p",ll?ll->head:NULL); fflush(NULL); */
      Free(ll);
      ll = next;
      if(ll) next = ll->next;
    }
  while(ll);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Free_String(t_string *ts)
{
  t_string *next;

  next = ts->next;
  do
    {
      Free(ts->s);
      Free(ts);

      ts = next;
      if(ts) next = ts->next;
    }
  while(ts);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Free_Custom_Model(t_mod *mod)
{
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Free_Efrq(t_efrq *e_frq)
{
  Free(e_frq->pi->v);
  Free(e_frq->pi);

  Free(e_frq->pi_unscaled->v);
  Free(e_frq->pi_unscaled);

  Free(e_frq->user_b_freq->v);
  Free(e_frq->user_b_freq);

  if(e_frq->next) Free_Efrq(e_frq->next);

  Free(e_frq);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Free_Rmat(t_rmat *r_mat)
{
  
  Free(r_mat->rr->v);
  Free(r_mat->rr);

  Free(r_mat->rr_num->v);

  Free(r_mat->rr_val->v);
  Free(r_mat->rr_val);

  Free(r_mat->n_rr_per_cat->v);
  Free(r_mat->n_rr_per_cat);

  Free(r_mat->rr_num);

  Free(r_mat->qmat->v);
  Free(r_mat->qmat);

  Free(r_mat->qmat_buff->v);
  Free(r_mat->qmat_buff);

  if(r_mat->next) Free_Rmat(r_mat->next);

  Free(r_mat);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Free_RAS(t_ras *ras)
{
  if(ras->gamma_r_proba->v)
    {
      Free(ras->gamma_r_proba->v);
      Free(ras->gamma_r_proba_unscaled->v);
      Free(ras->gamma_rr->v);
      Free(ras->gamma_rr_unscaled->v);
    }

  Free(ras->gamma_r_proba);
  Free(ras->skip_rate_cat);

  Free(ras->gamma_r_proba_unscaled);
  Free(ras->gamma_rr);
  Free(ras->gamma_rr_unscaled);
  Free_Scalar_Dbl(ras->pinvar);
  Free_Scalar_Dbl(ras->alpha);
  Free_Scalar_Dbl(ras->free_rate_mr);

  if(ras->next) Free_RAS(ras->next);

  Free(ras);
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Free_Model(t_mod *mod)
{
  Free_Custom_Model(mod);
  Free_Model_Complete(mod);
  if(mod->m4mod) M4_Free_M4_Model(mod->m4mod);
  Free_Model_Basic(mod);
  /* Free(mod); */
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Free(void *p)
{
#if (defined(__AVX__) || defined(__SSE3__))
#ifndef WIN32
  free(p);
#else
  _aligned_free(p);
#endif
#else  
  free(p);
#endif

  p = NULL;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Free_Input(option *io)
{
  int i;

  do
    {
      RATES_Free_Rates(io->rates);
      MCMC_Free_MCMC(io->mcmc);
      Free(io->in_align_file);
      Free(io->in_tree_file);
      Free(io->in_constraint_tree_file);
      Free(io->in_coord_file);
      Free(io->out_file);
      Free(io->out_tree_file);
      Free(io->out_trees_file);
      Free(io->out_boot_tree_file);
      Free(io->out_boot_stats_file);
      Free(io->out_stats_file);
      Free(io->weight_file);
      Free(io->out_lk_file);
      Free(io->out_summary_file);
      Free(io->out_ps_file);
      Free(io->out_trace_file);
      Free(io->out_json_trace_file);
      Free(io->out_ancestral_file);
      Free(io->nt_or_cd);
      Free(io->run_id_string);
      Free(io->clade_list_file);
      for(i=0;i<T_MAX_ALPHABET;i++) Free(io->alphabet[i]);
      Free(io->alphabet);
      if(io->short_tax_names)
        {
          for(i=0;i<io->size_tax_names;i++)
            {
              Free(io->short_tax_names[i]);
              Free(io->long_tax_names[i]);
            }
          Free(io->long_tax_names);
          Free(io->short_tax_names);
        }
      Free_Tree_List(io->treelist);

      if(io->lon) Free(io->lon);
      if(io->lat) Free(io->lat);

      if(io->next)
        {
          io = io->next;
          Free(io->prev);
        }
      else
        {
          Free(io);
          break;
        }

    }while(1);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Free_Tree_List(t_treelist *list)
{
  Free(list->tree);
  Free(list);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Free_St(supert_tree *st)
{
  int i;

  For(i,2*st->tree->n_otu-1)
    Free_NNI(st->tree->a_edges[i]->nni);

  for(i=0;i<st->n_part;i++) Free(st->match_st_node_in_gt[i]);

  Free(st->match_st_node_in_gt);

  Free_Tree(st->tree);

  Free(st);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Free_Eigen(eigen *eigen_struct)
{
  Free(eigen_struct->space_int);
  Free(eigen_struct->space);
  Free(eigen_struct->e_val);
  Free(eigen_struct->e_val_im);
  Free(eigen_struct->r_e_vect);
  Free(eigen_struct->r_e_vect_im);
  Free(eigen_struct->l_e_vect);
  Free(eigen_struct->q);

  if(eigen_struct->next) Free_Eigen(eigen_struct->next);

  Free(eigen_struct);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Free_One_Spr(t_spr *this_spr)
{
  Free(this_spr->path);
  if(this_spr->l0) Free_Scalar_Dbl(this_spr->l0);
  if(this_spr->l1) Free_Scalar_Dbl(this_spr->l1);
  if(this_spr->l2) Free_Scalar_Dbl(this_spr->l2);
  if(this_spr->v0) Free_Scalar_Dbl(this_spr->v0);
  if(this_spr->v1) Free_Scalar_Dbl(this_spr->v1);
  if(this_spr->v2) Free_Scalar_Dbl(this_spr->v2);
  
  if(this_spr->init_target_l) Free_Scalar_Dbl(this_spr->init_target_l);
  if(this_spr->init_target_v) Free_Scalar_Dbl(this_spr->init_target_v);
  
  Free(this_spr);  
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Free_Spr_List(t_tree *mixt_tree)
{
  int i;
  t_tree *tree;

  tree = mixt_tree;
  do
    {
      For(i,tree->size_spr_list+1) Free_One_Spr(tree->spr_list[i]);
      Free(tree->spr_list);
      Free_One_Spr(tree->best_spr);
      tree = tree->next;
    }
  while(tree);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Free_Triplet(triplet *t)
{
  int i,j,k;

  Free(t->F_bc);
  Free(t->F_cd);
  Free(t->F_bd);
  Free(t->pi_bc);
  Free(t->pi_cd);
  Free(t->pi_bd);

  for(k=0;k<t->mod->ras->n_catg;k++)
    {
      for(i=0;i<t->size;i++)
    {
      for(j=0;j<t->size;j++) Free(t->core[k][i][j]);
      Free(t->core[k][i]);
    }
      Free(t->core[k]);
    }
  Free(t->core);

  for(i=0;i<t->size;i++)
    {
      for(j=0;j<t->size;j++) Free(t->p_one_site[i][j]);
      Free(t->p_one_site[i]);
    }
  Free(t->p_one_site);

  for(i=0;i<t->size;i++)
    {
      for(j=0;j<t->size;j++) Free(t->sum_p_one_site[i][j]);
      Free(t->sum_p_one_site[i]);
    }
  Free(t->sum_p_one_site);

  Free_Eigen(t->eigen_struct);

  if(t->next) Free_Triplet(t->next);

  Free(t);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Free_Actual_CSeq(calign *data)
{
  int i;
  for(i=0;i<data->n_otu;i++)
    {
      Free(data->c_seq[i]->state);
      Free(data->c_seq[i]->d_state);
      data->c_seq[i]->state = NULL;
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Free_Prefix_Tree(pnode *n, int size)
{
  int i;

  for(i=0;i<size;i++)
    {
      if(n->next[i])
    {
      Free_Prefix_Tree(n->next[i],size);
    }
    }
  Free_Pnode(n);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Free_Pnode(pnode *n)
{
  Free(n->next);
  Free(n);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Free_Optimiz(t_opt *s_opt)
{
  Free(s_opt);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Free_Nexus(option *io)
{
  int i,j;

  for(i=0;i<N_MAX_NEX_COM;i++)
    {
      For(j,io->nex_com_list[i]->nparm) Free_Nexus_Parm(io->nex_com_list[i]->parm[j]);
      Free(io->nex_com_list[i]->parm);
      Free(io->nex_com_list[i]->name);
      Free(io->nex_com_list[i]);
    }
  Free(io->nex_com_list);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Free_Nexus_Com(nexcom **com)
{
  int i;

  for(i=0;i<N_MAX_NEX_COM;i++)
    {
      Free(com[i]->parm);
      Free(com[i]->name);
      Free(com[i]);
    }
  Free(com);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Free_Nexus_Parm(nexparm *parm)
{
  Free(parm->value);
  Free(parm->name);
  Free(parm);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void XML_Free_XML_Tree(xml_node *node)
{
  if(!node) return;
  if(node->child) XML_Free_XML_Tree(node->child);
  if(node->next)  XML_Free_XML_Tree(node->next);
  XML_Free_XML_Node(node);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void XML_Free_XML_Node(xml_node *node)
{
  Free(node->id);
  Free(node->name);
  Free(node->value);
  XML_Free_XML_Ds(node->ds);
  XML_Free_XML_Attr(node->attr);
  Free(node);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void XML_Free_XML_Attr(xml_attr *attr)
{
  if(attr)
    {
      Free(attr->name);
      Free(attr->value);
      if(attr->next) XML_Free_XML_Attr(attr->next);
      Free(attr);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void XML_Free_XML_Ds(t_ds *ds)
{
  if(ds->next) XML_Free_XML_Ds(ds->next);
  Free(ds);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MCMC_Free_MCMC(t_mcmc *mcmc)
{
  int i;

  Free(mcmc->move_type);
  Free(mcmc->adjust_tuning);
  Free(mcmc->out_filename);
  Free(mcmc->move_weight);
  Free(mcmc->acc_move);
  Free(mcmc->run_move);
  Free(mcmc->prev_acc_move);
  Free(mcmc->prev_run_move);
  Free(mcmc->acc_rate);
  Free(mcmc->tune_move);
  for(i=0;i<mcmc->n_moves;i++) Free(mcmc->move_name[i]);
  Free(mcmc->move_name);
  Free(mcmc->ess_run);
  Free(mcmc->start_ess);
  Free(mcmc->ess);
  Free(mcmc->sampled_val);
  Free(mcmc->mode);
  Free(mcmc);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void M4_Free_M4_Model(m4 *m4mod)
{
  int i;

  if(m4mod->o_mats)
    {
      for(i=0;i<m4mod->n_h;i++) Free(m4mod->o_mats[i]);
      Free(m4mod->o_mats);
      Free(m4mod->h_mat);
      Free(m4mod->o_rr);
      Free(m4mod->h_rr);
      Free(m4mod->o_fq);
      Free(m4mod->h_fq);
      Free(m4mod->multipl);
      Free(m4mod->multipl_unscaled);
      Free(m4mod->h_fq_unscaled);
    }

  Free(m4mod);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Free_Rates(t_rate *rates)
{
  if(rates->is_allocated == YES)
    {
      int i;
      Free(rates->nd_r);
      Free(rates->br_r);
      Free(rates->buff_r);
      Free(rates->true_r);
      Free(rates->buff_t);
      Free(rates->nd_t);
      Free(rates->true_t);
      Free(rates->t_prior);
      Free(rates->t_mean);
      Free(rates->t_prior_min);
      Free(rates->t_prior_max);
      Free(rates->t_floor);
      Free(rates->t_has_prior);
      Free(rates->t_rank);
      Free(rates->dens);
      Free(rates->triplet);
      Free(rates->n_jps);
      Free(rates->t_jps);
      Free(rates->cond_var);
      Free(rates->invcov);
      Free(rates->ml_l);
      Free(rates->cur_l);
      Free(rates->u_ml_l);
      Free(rates->u_cur_l);
      Free(rates->cov_r);
      Free(rates->lca);
      Free(rates->trip_cond_cov);
      Free(rates->trip_reg_coeff);
      Free(rates->_2n_vect1);
      Free(rates->_2n_vect2);
      Free(rates->_2n_vect3);
      Free(rates->_2n_vect4);
      Free(rates->_2n_vect5);
      Free(rates->_2n2n_vect1);
      Free(rates->_2n2n_vect2);
      Free(rates->cov_l);
      Free(rates->mean_l);
      Free(rates->mean_r);
      Free(rates->mean_t);
      Free(rates->grad_l);
      Free(rates->reg_coeff);
      Free(rates->br_do_updt);
      Free(rates->cur_gamma_prior_mean);
      Free(rates->cur_gamma_prior_var);
      Free(rates->n_tips_below);
      Free(rates->time_slice_lims);
      Free(rates->n_time_slice_spans);
      Free(rates->curr_slice);
      Free(rates->has_survived);
      Free(rates->survival_rank);
      Free(rates->survival_dur);
      Free(rates->calib_prob);
      Free(rates->t_prior_min_ori);
      Free(rates->t_prior_max_ori);
      Free(rates->times_partial_proba);
      Free(rates->numb_calib_chosen);
      for(i=0;i<rates->n_cal;i++) Free_Calib(rates->a_cal[i]);
      Free(rates->a_cal);
    }
  Free(rates);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Free_Calib(t_cal *cal)
{
  if(!cal) return;
  else 
    {    
      int i;
      for(i=0;i<cal->clade_list_size;i++) Free(cal->clade_list[i]);
      Free(cal->clade_list);      
      Free(cal->id);
      Free(cal->alpha_proba_list);
      Free(cal);
    }
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void Free_Geo(t_geo *t)
{
  int i;
  Free(t->f_mat);
  Free(t->r_mat);
  Free(t->occup);
  Free(t->idx_loc);
  Free(t->sorted_nd);
  Free(t->cov);
  Free(t->idx_loc_beneath);
  for(i=0;i<t->ldscape_sz;i++) Free_Geo_Coord(t->coord_loc[i]);
  Free(t->coord_loc);
  Free(t);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void Free_Geo_Coord(t_geo_coord *t)
{
  Free(t->cpy->lonlat); 
  Free(t->cpy->id);
  Free(t->cpy);

  Free(t->lonlat); 
  Free(t->id);
  Free(t);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void Free_Disk(t_dsk *t)
{
  Free_Geo_Coord(t->centr);
  Free(t->ldsk_a);
  Free(t->id);
  Free(t);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void Free_Ldisk(t_ldsk *t)
{
  if(t == NULL) return;
  else
    {
      Free(t->next);
      Free_Geo_Coord(t->coord);
      if(t->cpy_coord) Free_Geo_Coord(t->cpy_coord);
      Free(t);
    }

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Free_Poly(t_poly *p)
{
  int i;
  for(i=0;i<p->n_poly_vert;i++) Free_Geo_Coord(p->poly_vert[i]);
  Free(p->poly_vert);
  Free(p);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Free_Mmod(t_phyrex_mod *mmod)
{
  if(mmod == NULL) return;
  Free_Geo_Coord(mmod->lim);
  Free(mmod);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void JSON_Free_Array(json_a *a)
{
  if(a->object) JSON_Free_Object(a->object);
  Free(a);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void JSON_Free_Object(json_o *o)
{
  if(o->sv) JSON_Free_StringVal(o->sv);
  if(o->next) JSON_Free_Object(o->next);
  Free(o);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void JSON_Free_StringVal(json_sv *sv)
{
  if(sv->string) Free(sv->string);
  if(sv->value) Free(sv->value);
  if(sv->object) JSON_Free_Object(sv->object);
  if(sv->array) JSON_Free_Array(sv->array);
  if(sv->next) JSON_Free_StringVal(sv->next);
  Free(sv);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Free_NNI(t_nni *t)
{
  if(t->init_l) Free_Scalar_Dbl(t->init_l);
  if(t->init_v) Free_Scalar_Dbl(t->init_v);

  if(t->best_l) Free_Scalar_Dbl(t->best_l);
  if(t->best_v) Free_Scalar_Dbl(t->best_v);
  
  if(t->l0) Free_Scalar_Dbl(t->l0);
  if(t->v0) Free_Scalar_Dbl(t->v0);

  if(t->l1) Free_Scalar_Dbl(t->l1);
  if(t->v1) Free_Scalar_Dbl(t->v1);

  if(t->l2) Free_Scalar_Dbl(t->l2);
  if(t->v2) Free_Scalar_Dbl(t->v2);
  
  Free(t);
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
