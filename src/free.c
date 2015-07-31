/*

PHYML :  a program that  computes maximum likelihood  phyLOGenies from
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

  For(i,mat->n_otu)
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
/*   For(i,len) */
/*     { */
/*       For(j,n_catg) Free((*p_lk)[i][j]); */
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
      For(j,3) Free(tree->a_nodes[i]->bip_node[j]);
      Free(tree->a_nodes[i]->bip_node);
    }
    }
  tree->has_bip = NO;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Free_Cseq(calign *data)
{
  int i;

  Free(data->invar);
  Free(data->wght);
  Free(data->ambigu);
  Free(data->b_frq);
  Free(data->sitepatt);
  For(i,data->n_otu)
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
  For(i,n_otu)
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
  Free_Cseq(cdata);
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
      For(i,3)
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
          /* printf("\n. FREE b->num: %d %p %p",tree->a_edges[i]->num,tree->a_edges[i]->pars_l,tree->a_edges[i]->pars_r); */
          Free_Edge_Pars(tree->a_edges[i]);
        }

      if(tree->n_root)
        {
          /* printf("\n. FREE b1->num: %d %p %p",tree->n_root->b[1]->num,tree->n_root->b[1]->pars_l,tree->n_root->b[1]->pars_r); */
          Free_Edge_Pars_Left(tree->n_root->b[1]);
          Free_Edge_Pars_Left(tree->n_root->b[2]);
          /* if(tree->n_root->b[1]->pars_r) Free_Edge_Pars_Rght(tree->n_root->b[1]); */
          /* if(tree->n_root->b[2]->pars_r) Free_Edge_Pars_Rght(tree->n_root->b[2]); */
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

      For(i,3) Free(tree->log_lks_aLRT[i]);
      Free(tree->log_lks_aLRT);

      Free(tree->unscaled_site_lk_cat);

      For(i,2*tree->n_otu-3) Free_Edge_Lk(tree->a_edges[i]);

      if(tree->n_root)
        {
          Free(tree->n_root->b[1]->nni);
          Free(tree->n_root->b[2]->nni);
          Free(tree->n_root->b[1]->Pij_rr);
          Free(tree->n_root->b[2]->Pij_rr);
          Free_Edge_Lk_Left(tree->n_root->b[1]);
          Free_Edge_Lk_Left(tree->n_root->b[2]);
          /* if(tree->n_root->b[1]->p_lk_rght) Free_Edge_Lk_Rght(tree->n_root->b[1]); */
          /* if(tree->n_root->b[2]->p_lk_rght) Free_Edge_Lk_Rght(tree->n_root->b[2]); */
        }
      else
        {
          Free_Edge_Lk(tree->a_edges[2*tree->n_otu-3]);
          Free_Edge_Lk(tree->a_edges[2*tree->n_otu-2]);
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

  if(b->p_lk_tip_r) Free(b->p_lk_tip_r);

  Free(b->sum_scale_rght_cat);
  Free(b->patt_id_rght);
  Free(b->p_lk_loc_rght);

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

  if(b->p_lk_tip_l) Free(b->p_lk_tip_l);

  Free(b->sum_scale_left_cat);
  Free(b->patt_id_left);
  Free(b->p_lk_loc_left);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Free_Edge_Lk(t_edge *b)
{
  Free(b->nni);
  Free(b->Pij_rr);
  Free_Edge_Lk_Left(b);
  Free_Edge_Lk_Rght(b);
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

void Free_Model_Basic(t_mod *mixt_mod)
{
  t_mod *mod;

  Free_RAS(mixt_mod->ras);
  Free_Vect_Dbl(mixt_mod->user_b_freq);
  Free_Scalar_Dbl(mixt_mod->mr);
  Free_Scalar_Dbl(mixt_mod->kappa);
  Free_Scalar_Dbl(mixt_mod->lambda);
  Free_Scalar_Dbl(mixt_mod->br_len_mult);
  Free_Scalar_Dbl(mixt_mod->br_len_mult_unscaled);
  Free_Scalar_Dbl(mixt_mod->e_frq_weight);
  Free_Scalar_Dbl(mixt_mod->r_mat_weight);
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
  free(p);
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
      Free(io->out_lk_file);
      Free(io->out_summary_file);
      Free(io->out_ps_file);
      Free(io->out_trace_file);
      Free(io->out_ancestral_file);
      Free(io->nt_or_cd);
      Free(io->run_id_string);
      Free(io->clade_list_file);
      For(i,T_MAX_ALPHABET) Free(io->alphabet[i]);
      Free(io->alphabet);
      if(io->short_tax_names)
        {
          For(i,io->size_tax_names)
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
    Free(st->tree->a_edges[i]->nni);

  For(i,st->n_part) Free(st->match_st_node_in_gt[i]);

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

  For(k,t->mod->ras->n_catg)
    {
      For(i,t->size)
    {
      For(j,t->size) Free(t->core[k][i][j]);
      Free(t->core[k][i]);
    }
      Free(t->core[k]);
    }
  Free(t->core);

  For(i,t->size)
    {
      For(j,t->size) Free(t->p_one_site[i][j]);
      Free(t->p_one_site[i]);
    }
  Free(t->p_one_site);

  For(i,t->size)
    {
      For(j,t->size) Free(t->sum_p_one_site[i][j]);
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
  For(i,data->n_otu)
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

  For(i,size)
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

  For(i,N_MAX_NEX_COM)
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

  For(i,N_MAX_NEX_COM)
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
  For(i,mcmc->n_moves) Free(mcmc->move_name[i]);
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
      For(i,m4mod->n_h) Free(m4mod->o_mats[i]);
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
      Free(rates->t_ranked);
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
      Free(rates->node_height_dens_log_norm_const_update);
      Free(rates->curr_nd_for_cal);
      Free(rates->t_prior_min_ori);
      Free(rates->t_prior_max_ori);
      Free(rates->times_partial_proba);
      Free(rates->numb_calib_chosen);
    }
  Free_Calib(rates->calib);
  Free(rates);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Free_Calib(t_cal *cal)
{
  if(!cal) return;
  else Free_Calib(cal->next);
  Free(cal);
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
  For(i,t->ldscape_sz) Free_Geo_Coord(t->coord_loc[i]);
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
      /* if(t->min_coord) Free_Geo_Coord(t->min_coord); */
      /* if(t->max_coord) Free_Geo_Coord(t->max_coord); */
      if(t->cpy_coord) Free_Geo_Coord(t->cpy_coord);
      Free(t);
    }

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Free_Poly(t_poly *p)
{
  int i;
  For(i,p->n_poly_vert) Free_Geo_Coord(p->poly_vert[i]);
  Free(p->poly_vert);
  Free(p);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Free_Mmod(t_phyrex_mod *mmod)
{
  Free_Geo_Coord(mmod->lim);
  Free(mmod);
}
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
