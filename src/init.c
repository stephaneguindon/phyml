/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include "init.h"
#ifdef BEAGLE
#include "beagle_utils.h"
#endif

void Init_Eigen_Struct(eigen *this)
{
  this->next = NULL;
  this->prev = NULL;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Init_Scalar_Dbl(scalar_dbl *p)
{
  p->v     = -1.;
  p->onoff = ON;
  p->next  = NULL;
  p->prev  = NULL;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Init_Scalar_Int(scalar_int *p)
{
  p->v    = -1.;
  p->next = NULL;
  p->prev = NULL;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Init_Vect_Dbl(int len, vect_dbl *p)
{
  p->len  = len;
  p->next = NULL;
  p->prev = NULL;
  p->v    = NULL;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Init_Vect_Int(int len, vect_int *p)
{
  p->len  = len;
  p->next = NULL;
  p->prev = NULL;
  p->v    = NULL;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Init_String(t_string *ts)
{
  ts->len  = -1.;
  ts->next = NULL;
  ts->prev = NULL;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Init_Triplet_Struct(triplet *t)
{
  t->next = NULL;
  t->prev = NULL;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Init_Efrq(t_efrq *f)
{
  f->next = NULL;
  f->prev = NULL;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Init_Tree(t_tree *tree, int n_otu)
{
  tree->n_otu                     = n_otu;
  tree->mat                       = NULL;
  tree->n_root                    = NULL;
  tree->e_root                    = NULL;
  tree->ps_tree                   = NULL;
  tree->short_l                   = NULL;
  tree->mutmap                    = NULL;
  tree->next                      = NULL;
  tree->prev                      = NULL;
  tree->next                      = NULL;
  tree->prev                      = NULL;
  tree->mixt_tree                 = NULL;
  tree->geo                       = NULL;

  tree->is_mixt_tree              = NO;
  tree->tree_num                  = 0;
  tree->depth_curr_path           = 0;
  tree->has_bip                   = NO;
  tree->n_moves                   = 0;
  tree->n_improvements            = 0;
  tree->bl_from_node_stamps       = 0;
  tree->lock_topo                 = NO;
  tree->ps_page_number            = 0;
  tree->init_lnL                  = UNLIKELY;
  tree->best_lnL                  = UNLIKELY;
  tree->old_lnL                   = UNLIKELY;
  tree->c_lnL                     = UNLIKELY;
  tree->sum_min_sum_scale         = .0;
  tree->n_swap                    = 0;
  tree->best_pars                 = 1E+5;
  tree->n_pattern                 = -1;
  tree->n_root_pos                = -1.;
  tree->write_labels              = YES;
  tree->write_br_lens             = YES;
  tree->print_boot_val            = 0;
  tree->print_alrt_val            = 0;
  tree->num_curr_branch_available = 0;
  tree->tip_order_score           = .0;
  tree->write_tax_names           = YES;
  tree->update_alias_subpatt      = NO;
  tree->bl_ndigits                = 8;
  tree->n_short_l                 = 100;
  tree->norm_scale                = 0.0;
  tree->br_len_recorded           = NO;
  tree->max_spr_depth             = 0;
  tree->apply_lk_scaling          = YES;
  tree->dp                        = 0;
  tree->ignore_root               = YES;
  tree->annealing_temp            = 0.;
#ifdef BEAGLE
  tree->b_inst                    = UNINITIALIZED;
#endif
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Init_Edge_Light(t_edge *b, int num)
{
  b->num                  = num;
  b->bip_score            = 0;
  b->dist_btw_edges       = .0;
  b->topo_dist_btw_edges  = 0;
  b->has_zero_br_len      = NO;
  b->n_jumps              = 0;
  b->l_var->v             = -1.;
  b->does_exist           = YES;
  b->l->v                 = -1.;
  b->bin_cod_num          = -1.;
  b->l->onoff             = ON;

  b->next                 = NULL;
  b->prev                 = NULL;
  b->next_mixt            = NULL;
  b->prev_mixt            = NULL;

  b->p_lk_left            = NULL;
  b->p_lk_rght            = NULL;
  b->p_lk_loc_left        = NULL;
  b->p_lk_loc_rght        = NULL;
  b->Pij_rr               = NULL;
  b->labels               = NULL;

  b->pars_l               = NULL;
  b->pars_r               = NULL;
  b->ui_l                 = NULL;
  b->ui_r                 = NULL;
  b->p_pars_l             = NULL;
  b->p_pars_r             = NULL;
  b->n_diff_states_l      = NULL;
  b->n_diff_states_r      = NULL;

#ifdef BEAGLE
  b->p_lk_left_idx         = num;
  b->p_lk_rght_idx         = UNINITIALIZED; //Will be initialized later when the total number of branches is known (i.e. in Make_Tree_From_Scratch())
  b->Pij_rr_idx            = num;
  b->p_lk_tip_idx          = UNINITIALIZED; //Will be initialized later only if this branch is connected to a tip
#endif

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Init_Node_Light(t_node *n, int num)
{
  n->num                    = num;
  n->tax                    = -1;
  n->dist_to_root           = .0;
  n->common                 = 1;
  n->ext_node               = NULL;
  n->name                   = NULL;
  n->ori_name               = NULL;
  n->c_seq                  = NULL;
  n->c_seq_anc              = NULL;
  n->y_rank                 = 0.;
  n->y_rank_ori             = 0.;
  n->y_rank_max             = 0.;
  n->y_rank_min             = 0.;
  n->anc                    = NULL;
  n->rank                   = 0;
  n->match_node             = NULL;
  n->id_rank                = 0;
  n->next                   = NULL;
  n->prev                   = NULL;
  /* n->next                  = NULL; */
  /* n->prev                 = NULL; */
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Init_NNI(nni *a_nni)
{
  a_nni->left         = NULL;
  a_nni->rght         = NULL;
  a_nni->b            = NULL;
  a_nni->init_l       = -1.;
  a_nni->init_lk      = .0;
  a_nni->score        = +1.0;
  a_nni->best_l       = -1.;
  a_nni->swap_node_v1 = NULL;
  a_nni->swap_node_v2 = NULL;
  a_nni->swap_node_v3 = NULL;
  a_nni->swap_node_v4 = NULL;
  a_nni->lk0          = UNLIKELY;
  a_nni->lk1          = UNLIKELY;
  a_nni->lk2          = UNLIKELY;
  a_nni->l0           = -1.0;
  a_nni->l1           = -1.0;
  a_nni->l2           = -1.0;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Init_Nexus_Format(nexcom **com)
{

  /*****************************/

  strcpy(com[0]->name,"dimensions");
  com[0]->nparm = 2;
  com[0]->nxt_token_t = NEXUS_PARM;
  com[0]->cur_token_t = NEXUS_COM;

  com[0]->parm[0] = Make_Nexus_Parm();
  strcpy(com[0]->parm[0]->name,"ntax");
  com[0]->parm[0]->fp = Read_Nexus_Dimensions;
  com[0]->parm[0]->com = com[0];
  com[0]->parm[0]->nxt_token_t = NEXUS_EQUAL;
  com[0]->parm[0]->cur_token_t = NEXUS_PARM;

  com[0]->parm[1] = Make_Nexus_Parm();
  strcpy(com[0]->parm[1]->name,"nchar");
  com[0]->parm[1]->fp = Read_Nexus_Dimensions;
  com[0]->parm[1]->com = com[0];
  com[0]->parm[1]->nxt_token_t = NEXUS_EQUAL;
  com[0]->parm[1]->cur_token_t = NEXUS_PARM;

  /*****************************/

  strcpy(com[1]->name,"format");
  com[1]->nparm = 11;
  com[1]->nxt_token_t = NEXUS_PARM;
  com[1]->cur_token_t = NEXUS_COM;

  com[1]->parm[0] = Make_Nexus_Parm();
  strcpy(com[1]->parm[0]->name,"datatype");
  com[1]->parm[0]->fp = Read_Nexus_Format;
  com[1]->parm[0]->com = com[1];
  com[1]->parm[0]->nxt_token_t = NEXUS_EQUAL;
  com[1]->parm[0]->cur_token_t = NEXUS_PARM;

  com[1]->parm[1] = Make_Nexus_Parm();
  strcpy(com[1]->parm[1]->name,"respectcase");
  com[1]->parm[1]->fp = Read_Nexus_Format;
  com[1]->parm[1]->com = com[1];
  com[1]->parm[1]->nxt_token_t = NEXUS_PARM;
  com[1]->parm[1]->cur_token_t = NEXUS_VALUE;

  com[1]->parm[2] = Make_Nexus_Parm();
  strcpy(com[1]->parm[2]->name,"missing");
  com[1]->parm[2]->fp = Read_Nexus_Format;
  com[1]->parm[2]->com = com[1];
  com[1]->parm[2]->nxt_token_t = NEXUS_EQUAL;
  com[1]->parm[2]->cur_token_t = NEXUS_PARM;

  com[1]->parm[3] = Make_Nexus_Parm();
  strcpy(com[1]->parm[3]->name,"gap");
  com[1]->parm[3]->fp = Read_Nexus_Format;
  com[1]->parm[3]->com = com[1];
  com[1]->parm[3]->nxt_token_t = NEXUS_EQUAL;
  com[1]->parm[3]->cur_token_t = NEXUS_PARM;

  com[1]->parm[4] = Make_Nexus_Parm();
  strcpy(com[1]->parm[4]->name,"symbols");
  com[1]->parm[4]->fp = Read_Nexus_Format;
  com[1]->parm[4]->com = com[1];
  com[1]->parm[4]->nxt_token_t = NEXUS_EQUAL;
  com[1]->parm[4]->cur_token_t = NEXUS_PARM;

  com[1]->parm[5] = Make_Nexus_Parm();
  strcpy(com[1]->parm[5]->name,"equate");
  com[1]->parm[5]->fp = Read_Nexus_Format;
  com[1]->parm[5]->com = com[1];
  com[1]->parm[5]->nxt_token_t = NEXUS_EQUAL;
  com[1]->parm[5]->cur_token_t = NEXUS_PARM;

  com[1]->parm[6] = Make_Nexus_Parm();
  strcpy(com[1]->parm[6]->name,"matchchar");
  com[1]->parm[6]->fp = Read_Nexus_Format;
  com[1]->parm[6]->com = com[1];
  com[1]->parm[6]->nxt_token_t = NEXUS_EQUAL;
  com[1]->parm[6]->cur_token_t = NEXUS_PARM;

  com[1]->parm[7] = Make_Nexus_Parm();
  strcpy(com[1]->parm[7]->name,"transpose");
  com[1]->parm[7]->fp = Read_Nexus_Format;
  com[1]->parm[7]->com = com[1];
  com[1]->parm[7]->nxt_token_t = NEXUS_PARM;
  com[1]->parm[7]->cur_token_t = NEXUS_VALUE;

  com[1]->parm[8] = Make_Nexus_Parm();
  strcpy(com[1]->parm[8]->name,"interleave");
  com[1]->parm[8]->fp = Read_Nexus_Format;
  com[1]->parm[8]->com = com[1];
  com[1]->parm[8]->nxt_token_t = NEXUS_PARM;
  com[1]->parm[8]->cur_token_t = NEXUS_VALUE;

  com[1]->parm[9] = Make_Nexus_Parm();
  strcpy(com[1]->parm[9]->name,"items");
  com[1]->parm[9]->fp = Read_Nexus_Format;
  com[1]->parm[9]->com = com[1];
  com[1]->parm[9]->nxt_token_t = NEXUS_EQUAL;
  com[1]->parm[9]->cur_token_t = NEXUS_PARM;

  com[1]->parm[10] = Make_Nexus_Parm();
  strcpy(com[1]->parm[10]->name,"statesformat");
  com[1]->parm[10]->fp = Read_Nexus_Format;
  com[1]->parm[10]->com = com[1];
  com[1]->parm[10]->nxt_token_t = NEXUS_EQUAL;
  com[1]->parm[10]->cur_token_t = NEXUS_PARM;

  /*****************************/

  strcpy(com[2]->name,"eliminate");
  com[2]->nparm = 0;
  com[2]->nxt_token_t = NEXUS_VALUE;
  com[2]->cur_token_t = NEXUS_COM;

  /*****************************/

  strcpy(com[3]->name,"taxlabels");
  com[3]->nparm = 0;
  com[3]->nxt_token_t = -1;
  com[3]->cur_token_t = -1;

 /*****************************/

  strcpy(com[4]->name,"charstatelabels");
  com[4]->nparm = 0;
  com[4]->nxt_token_t = -1;
  com[4]->cur_token_t = -1;

  /*****************************/

  strcpy(com[5]->name,"charlabels");
  com[5]->nparm = 0;
  com[5]->nxt_token_t = -1;
  com[5]->cur_token_t = -1;

  /*****************************/

  strcpy(com[6]->name,"statelabels");
  com[6]->nparm = 0;
  com[6]->nxt_token_t = -1;
  com[6]->cur_token_t = -1;

  /*****************************/

  strcpy(com[7]->name,"matrix");
  com[7]->nparm = 1;
  com[7]->nxt_token_t = NEXUS_COM;
  com[7]->cur_token_t = NEXUS_VALUE; /* This will allow us to skip directly
                    to the matrix reading function */

  com[7]->parm[0] = Make_Nexus_Parm();
  strcpy(com[7]->parm[0]->name,"matrix");
  com[7]->parm[0]->fp = Read_Nexus_Matrix;
  com[7]->parm[0]->com = com[7];
  com[7]->parm[0]->nxt_token_t = NEXUS_COM;
  com[7]->parm[0]->cur_token_t = -1;

  /*****************************/

  strcpy(com[8]->name,"begin");
  com[8]->nparm = 3;

  com[8]->nxt_token_t = NEXUS_PARM;
  com[8]->cur_token_t = NEXUS_COM;

  com[8]->parm[0] = Make_Nexus_Parm();
  strcpy(com[8]->parm[0]->name,"data");
  com[8]->parm[0]->fp = Read_Nexus_Begin;
  com[8]->parm[0]->com = com[8];
  com[8]->parm[0]->nxt_token_t = NEXUS_COM;
  com[8]->parm[0]->cur_token_t = NEXUS_PARM;


  com[8]->parm[1] = Make_Nexus_Parm();
  strcpy(com[8]->parm[1]->name,"trees");
  com[8]->parm[1]->fp = Read_Nexus_Begin;
  com[8]->parm[1]->com = com[8];
  com[8]->parm[1]->nxt_token_t = NEXUS_COM;
  com[8]->parm[1]->cur_token_t = NEXUS_PARM;


  com[8]->parm[2] = Make_Nexus_Parm();
  strcpy(com[8]->parm[2]->name,"taxa");
  com[8]->parm[2]->fp = Read_Nexus_Taxa;
  com[8]->parm[2]->com = com[8];
  com[8]->parm[2]->nxt_token_t = NEXUS_COM;
  com[8]->parm[2]->cur_token_t = NEXUS_VALUE;


  /*****************************/

  strcpy(com[9]->name,"end");
  com[9]->nparm = 0;
  com[9]->nxt_token_t = -1;
  com[9]->cur_token_t = -1;

  /*****************************/

  strcpy(com[10]->name,"translate");
  com[10]->nparm = 1;
  com[10]->nxt_token_t = NEXUS_COM;
  com[10]->cur_token_t = NEXUS_VALUE;

  com[10]->parm[0] = Make_Nexus_Parm();
  strcpy(com[10]->parm[0]->name,"translate");
  com[10]->parm[0]->fp = Read_Nexus_Translate;
  com[10]->parm[0]->com = com[10];
  com[10]->parm[0]->nxt_token_t = NEXUS_COM;
  com[10]->parm[0]->cur_token_t = -1;

  /*****************************/

  strcpy(com[11]->name,"tree");
  com[11]->nparm = 1;
  com[11]->nxt_token_t = NEXUS_COM;
  com[11]->cur_token_t = NEXUS_VALUE;

  com[11]->parm[0] = Make_Nexus_Parm();
  strcpy(com[11]->parm[0]->name,"tree");
  com[11]->parm[0]->fp = Read_Nexus_Tree;
  com[11]->parm[0]->com = com[11];
  com[11]->parm[0]->nxt_token_t = -1;
  com[11]->parm[0]->cur_token_t = -1;


  /*****************************/

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Init_Mat(matrix *mat, calign *data)
{
  int i;

  mat->n_otu = data->n_otu;
  mat->r = mat->n_otu;
  mat->curr_int = mat->n_otu;
  mat->method = 1;

  For(i,data->n_otu)
    {
      strcpy(mat->name[i],data->c_seq[i]->name);
      mat->on_off[i] = 1;
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Set_Defaults_Input(option* io)
{
  io->fp_in_align                = NULL;
  io->fp_in_tree                 = NULL;
  io->fp_in_constraint_tree      = NULL;
  io->fp_out_tree                = NULL;
  io->fp_out_trees               = NULL;
  io->fp_out_boot_tree           = NULL;
  io->fp_out_boot_stats          = NULL;
  io->fp_out_stats               = NULL;
  io->fp_out_ancestral           = NULL;
  io->fp_in_coord                = NULL;
  io->long_tax_names             = NULL;
  io->short_tax_names            = NULL;
  io->lon                        = NULL;
  io->lat                        = NULL;
  io->z_scores                   = NULL;
  io->cstr_tree                  = NULL;
  io->next                       = NULL;
  io->prev                       = NULL;
  io->tree                       = NULL;
  io->mod                        = NULL;
  strcpy(io->nt_or_cd,"nucleotides");
  io->n_data_sets                = 1;
  io->interleaved                = 1;
  io->in_tree                    = 0;
  io->out_tree_file_open_mode    = 1;
  io->out_stats_file_open_mode   = 1;
  io->init_len                   = -1;
  io->n_otu                      = -1;
  io->n_data_set_asked           = -1;
  io->print_boot_trees           = 1;
  io->n_part                     = 1;
  io->ratio_test		 = ABAYES;
  io->multigene                  = 0;
  io->config_multigene           = 0;
  io->curr_interface             = 0;
  io->r_seed                     = -1;
  io->collapse_boot              = 0;
  io->random_boot_seq_order      = YES;
  io->print_trace                = 0;
  io->print_site_lnl             = 0;
  io->m4_model                   = NO;
  io->rm_ambigu                  = 0;
  io->append_run_ID              = 0;
  io->quiet                      = 0;
  io->datatype                   = NT;
  io->colalias                   = YES;
  io->data_file_format           = PHYLIP;
  io->tree_file_format           = PHYLIP;
  io->boot_prog_every            = 20;
  io->mem_question               = YES;
  io->do_alias_subpatt           = NO;
  io->lk_approx                  = EXACT;
  io->codpos                     = -1;
  io->mutmap                     = NO;
  io->state_len                  = 1;
  io->ancestral                  = NO;
#ifdef BEAGLE
  io->beagle_resource            = 0;
#endif

  MCMC_Init_MCMC_Struct(NULL,io,io->mcmc);
  RATES_Init_Rate_Struct(io->rates,NULL,-1);
  io->rates->model               = GUINDON;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Init_Rmat(t_rmat *rmat)
{
  rmat->n_diff_rr = 1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Set_Defaults_Model(t_mod *mod)
{
  Set_Defaults_Ras(mod->ras);

  strcpy(mod->modelname->s,"HKY85");
  strcpy(mod->custom_mod_string->s,"000000");
  mod->next                    = NULL;
  mod->prev                    = NULL;
  mod->next_mixt               = NULL;
  mod->prev_mixt               = NULL;
  mod->r_mat                   = NULL;
  mod->e_frq                   = NULL;
  mod->whichmodel              = HKY85;
  mod->n_mixt_classes          = 0;
  mod->mod_num                 = 0;
  mod->update_eigen            = NO;
  mod->is_mixt_mod             = NO;
  
  mod->kappa->v                = 4.0;
  mod->lambda->v               = 1.0;
  mod->l_var_sigma             = 1.E-2;
  mod->e_frq_weight->v         = 1.0;
  mod->r_mat_weight->v         = 1.0;

  mod->bootstrap               = 0;
  mod->ns                      = 4;
  mod->use_m4mod               = NO;
  mod->ras->gamma_median       = NO;
  mod->m4mod                   = NULL;

  /* mod->r_mat->n_diff_rr        = 0; */
  /* mod->r_mat->rr               = NULL; */
  /* mod->r_mat->rr_val           = NULL; */
  /* mod->r_mat->n_rr_per_cat     = NULL; */

  mod->io                       = NULL;
  mod->log_l                    = NO;
  mod->gamma_mgf_bl             = NO;
  mod->br_len_mult->v     = 1.0;


#if !(defined PHYTIME || defined INVITEE)
  mod->l_min = 1.E-8;
  mod->l_max = 10.0;
  /* mod->l_min = 1.E-4; */
  /* mod->l_max = 2.0; */
#else
  mod->l_min = 1.E-8;
  mod->l_max = 2.0;
#endif

  mod->l_var_min  = mod->l_min;
  mod->l_var_max  = mod->l_max;

  mod->br_len_mult->v          = 1.0;
  mod->br_len_mult_unscaled->v = 1.0;
  mod->augmented               = NO;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Set_Defaults_Ras(t_ras *ras)
{
  ras->n_catg              = 4;
  ras->normalise_rr        = YES;
  ras->pinvar->v           = 0.0;
  ras->alpha->v            = 1.0;
  ras->invar               = NO;
  ras->free_mixt_rates     = NO;
  ras->parent_class_number = 0;
  ras->init_r_proba        = YES;
  ras->init_rr             = YES;
  ras->sort_rate_classes   = NO;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Set_Defaults_Optimiz(t_opt *s_opt)
{
  s_opt->print                = YES;
  s_opt->last_opt             = YES;
  s_opt->opt_subst_param      = YES;
  s_opt->opt_alpha            = YES;
  s_opt->opt_kappa            = YES;
  s_opt->opt_bl               = YES;
  s_opt->opt_lambda           = NO;
  s_opt->opt_pinvar           = NO;
  s_opt->opt_cov_delta        = NO;
  s_opt->opt_cov_alpha        = NO;
  s_opt->opt_cov_free_rates   = NO;
  s_opt->opt_rr               = NO;
  s_opt->init_lk              = UNLIKELY;
  s_opt->n_it_max             = 1000;
  s_opt->opt_topo             = YES;
  s_opt->topo_search          = NNI_MOVE;
  s_opt->random_input_tree    = 0;
  s_opt->n_rand_starts        = 5;
  s_opt->brent_it_max         = BRENT_IT_MAX;
  s_opt->steph_spr            = YES;
  s_opt->user_state_freq      = NO;
  s_opt->opt_br_len_mult      = NO;

  /* s_opt->min_diff_lk_local    = 1.E-04; */
  /* s_opt->min_diff_lk_global   = 1.E-03; */
  /* s_opt->min_diff_lk_move     = 1.E-02; */
  s_opt->min_diff_lk_local    = 1.E-02;
  s_opt->min_diff_lk_global   = 1.E-02;
  s_opt->min_diff_lk_move     = 1.E-02;

  s_opt->p_moves_to_examine   = 0.15;
  s_opt->fast_nni             = NO;
  s_opt->greedy               = NO;
  s_opt->general_pars         = NO;
  s_opt->tree_size_mult       = 1;
  s_opt->opt_five_branch      = YES;

  s_opt->pars_thresh          = 5;

  s_opt->hybrid_thresh        = NO;
  s_opt->quickdirty           = NO;
  s_opt->spr_pars             = YES;
  s_opt->spr_lnL              = NO;
  s_opt->min_depth_path       = 0;

  s_opt->max_depth_path       = 20;
  s_opt->deepest_path         = 20;
  s_opt->max_delta_lnL_spr    = 50.;

  s_opt->br_len_in_spr        = 10;
  s_opt->opt_free_mixt_rates  = YES;
  s_opt->constrained_br_len   = NO;
  s_opt->opt_gamma_br_len     = NO;
  s_opt->first_opt_free_mixt_rates = YES;

  s_opt->wim_n_rgrft          = -1;
  s_opt->wim_n_globl          = -1;
  s_opt->wim_max_dist         = -1;
  s_opt->wim_n_optim          = -1;
  s_opt->wim_n_best           = -1;
  s_opt->wim_inside_opt       =  0;

  s_opt->opt_rmat_weight      = NO;
  s_opt->opt_efrq_weight      = NO;

  s_opt->skip_tree_traversal     = NO;
  s_opt->serial_free_rates       = YES;

  s_opt->curr_opt_free_rates     = NO;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void XML_Init_Attribute(xml_attr *attr)
{
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void XML_Init_Node(xml_node *parent, xml_node *new_node, char *name)
{
  if(name) strcpy(new_node->name,name);

  new_node->parent   = parent ? parent : NULL;
  new_node->next     = NULL;
  new_node->prev     = NULL;
  new_node->child    = NULL;
  new_node->ds->obj  = NULL;
  new_node->ds->next = NULL;

  if(parent)
    {
      if(!parent->child)
        {
          parent->child = new_node;
        }
      else
        {
          xml_node *node = parent->child;
          while(node->next) node = node->next;
          node->next = new_node;
          new_node->prev = node;
        }
    }

  new_node->attr = NULL;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Init_Rate_Struct(t_rate *rates, t_rate *existing_rates, int n_otu)
{
  int i;

  if(existing_rates && (existing_rates->model != -1))
    {
      rates->model = existing_rates->model;
    }
  else
    {
      rates->model = NONE;
    }


  if(rates->model == NONE)
    rates->model_log_rates = NO;
  else if(rates->model == THORNE)
    rates->model_log_rates = YES;
  else if(rates->model == GUINDON)
    rates->model_log_rates = YES;
  else if(rates->model == GAMMA)
    rates->model_log_rates = NO;
  else if(rates->model == STRICTCLOCK)
    rates->model_log_rates = NO;
  else
    {
      PhyML_Printf("\n== Please initialize model properly.");
      PhyML_Printf("\n== Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }


  rates->met_within_gibbs = NO;
  rates->c_lnL_rates      = UNLIKELY;
  rates->c_lnL_jps        = UNLIKELY;
  rates->adjust_rates     = 0;
  rates->use_rates        = 1;
  rates->lexp             = 1.E-3;
  rates->norm_fact        = 1.0;
  rates->inflate_var      = 1.0;
  rates->nd_t_recorded    = NO;
  rates->br_r_recorded    = NO;

  rates->birth_rate       = 1.E-2;
  rates->birth_rate_min   = 1.E-6;
  rates->birth_rate_max   = 1.E+1;

  if(rates->model_log_rates == YES)
    {
      rates->max_rate  =  LOG(10.);
      rates->min_rate  = -LOG(10.);
      /* rates->max_rate  =  MDBL_MAX; */
      /* rates->min_rate  = -MDBL_MAX; */
    }
  else
    {
      rates->max_rate  = 10.0;
      rates->min_rate  = 0.0;
    }
  /* rates->max_rate         = 6.0; */
  /* rates->min_rate         = 0.0; */


  rates->clock_r       = 1.E-4;
  rates->min_clock     = 1.E-8;
  rates->max_clock     = 1.E-3;

  /* rates->clock_r       = 3.E-4; */
  /* rates->max_clock     = 1.E-3; */
  /* rates->min_clock     = 1.E-5; */

  if(rates->model != GAMMA)
    {
      rates->nu            = 1.E-3;
      rates->min_nu        = 0.0;
      rates->max_nu        = 1.0;
    }
  else
    {
      rates->nu            = 1.E-1;
      rates->min_nu        = 1.E-2;
      rates->max_nu        = 100.0;
    }

  /* rates->max_nu        = 2.0; */

  /* rates->nu            = 1.E-4; */
  /* rates->max_nu        = 1.E-1; */
  /* rates->min_nu        = 1.E-5; */

  rates->min_dt        = 0.0;

  rates->step_rate     = 1.E-4;
  rates->approx        = 1;
  rates->bl_from_rt    = NO;

  rates->update_mean_l = NO;
  rates->update_cov_l  = NO;

  rates->p_max         = 0.01;

  rates->true_tree_size = 0.0;

  if(n_otu > 0)
    {
      For(i,(2*n_otu-2)*(2*n_otu-2)) rates->cov_l[i] = 0.0;
      
      For(i,2*n_otu-2)
        {
          rates->n_jps[i]  =  -1;
          rates->t_jps[i]  =  -1;
          rates->mean_r[i] = 1.0;
          rates->mean_l[i] = 0.0;
          rates->cur_l[i]  = 0.01;
        }
      
      For(i,2*n_otu-1)
        {
          if(rates->model_log_rates == YES)
            {
              rates->nd_r[i]   = 0.0;
              rates->br_r[i]   = 0.0;
            }
          else
            {
              rates->nd_r[i]   = 1.0;
              rates->br_r[i]   = 1.0;
            }
          
          rates->mean_t[i] = 0.0;
          rates->nd_t[i]   = 0.0;
          rates->true_t[i] = 0.0;
          if(i < n_otu)
            {
              rates->t_has_prior[i] = YES;
              rates->t_prior_max[i] = 0.0;
              rates->t_prior_min[i] = 0.0;
            }
          else
            {
              rates->t_has_prior[i] = NO;
              rates->t_prior_max[i] =  BIG;
              rates->t_prior_min[i] = -BIG;
            }
          
          rates->br_do_updt[i] = YES;
          rates->has_survived[i] = NO;
          
          rates->t_ranked[i] = i;
        }
    }
  
  rates->calib = NULL;
  rates->update_time_norm_const = NO;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Init_One_Spr(t_spr *a_spr)
{
  a_spr->lnL             = UNLIKELY;
  a_spr->pars            = 1E+5;
  a_spr->depth_path      = 0;
  a_spr->dist            = 0;
  a_spr->init_target_l   = -1.;
  a_spr->init_target_v   = -1.;
  a_spr->l0              = -1.;
  a_spr->l1              = -1.;
  a_spr->l2              = -1.;
  a_spr->v0              = -1.;
  a_spr->v1              = -1.;
  a_spr->v2              = -1.;
  a_spr->n_link          = NULL;
  a_spr->n_opp_to_link   = NULL;
  a_spr->b_opp_to_link   = NULL;
  a_spr->b_target        = NULL;
  a_spr->b_init_target   = NULL;
  a_spr->next            = NULL;
  a_spr->prev            = NULL;
  a_spr->next           = NULL;
  a_spr->prev          = NULL;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Init_Model(calign *data, t_mod *mod, option *io)
{
  int i,j;
  phydbl sum,aux;
  int result;
  phydbl *dr, *di, *space;

#ifdef BEAGLE
  mod->b_inst = UNINITIALIZED; //prevents calling an uninitialized BEAGLE instance (for ex: prevents Update_Eigen(), Update_RAS(), from calling BEAGLE)
  mod->optimizing_topology = false;
#endif

  mod->ns = io->mod->ns;

  if(io->datatype == GENERIC) mod->whichmodel = JC69;

  /* if(!mod->ras->invar) For(i,data->crunch_len) data->invar[i] = 0; */

  dr      = (phydbl *)mCalloc(  mod->ns,sizeof(phydbl));
  di      = (phydbl *)mCalloc(  mod->ns,sizeof(phydbl));
  space   = (phydbl *)mCalloc(2*mod->ns,sizeof(phydbl));

  if(mod->log_l == YES)
    {
      mod->l_min = LOG(mod->l_min);
      mod->l_max = LOG(mod->l_max);
    }

  if(mod->ras->init_r_proba == YES)
    {
      For(i,mod->ras->n_catg) mod->ras->gamma_r_proba->v[i]          = (phydbl)1./mod->ras->n_catg;
      For(i,mod->ras->n_catg) mod->ras->gamma_r_proba_unscaled->v[i] = (phydbl)(i+1);
    }
  else
    {
      mod->ras->gamma_r_proba_unscaled->v[0] = mod->ras->gamma_r_proba->v[0];
      for(i=1;i<mod->ras->n_catg;i++)
        mod->ras->gamma_r_proba_unscaled->v[i] =
          mod->ras->gamma_r_proba_unscaled->v[i-1] +
          mod->ras->gamma_r_proba->v[i];
    }

  if(mod->ras->init_rr == YES)
    {
      if(mod->ras->n_catg > 1)
        {
          For(i,mod->ras->n_catg) mod->ras->gamma_rr->v[i]          = (phydbl)i;
          For(i,mod->ras->n_catg) mod->ras->gamma_rr_unscaled->v[i] = (phydbl)i;
        }
      else
        {
          mod->ras->gamma_rr->v[0]          = 1.0;
          mod->ras->gamma_rr_unscaled->v[0] = 1.0;
        }
    }
  
  /* mod->br_len_mult->v          = 1.0; */
  /* mod->br_len_mult_unscaled->v = 1.0; */
  
  For(i,mod->ns)
    {
      mod->e_frq->pi->v[i] = data->b_frq[i];
      mod->e_frq->pi_unscaled->v[i] = mod->e_frq->pi->v[i] * 100.;
    }


  if(io->datatype == NT)
    {
      /* Set the substitution parameters to their default values
         if they are not fixed by the user */
      if(mod->s_opt->opt_kappa == YES)
        {
          mod->kappa->v  = 4.0;
          mod->lambda->v = 1.0;
        }
      
      if(mod->whichmodel == CUSTOM)
        {
          For(i,6) mod->r_mat->rr_val->v[i] = 1.0;

          /* Condition below is true for if custom model corresponds to TN93 or K80 */
          if(mod->r_mat->rr_num->v[AC] == mod->r_mat->rr_num->v[AT] &&
             mod->r_mat->rr_num->v[AT] == mod->r_mat->rr_num->v[CG] &&
             mod->r_mat->rr_num->v[CG] == mod->r_mat->rr_num->v[GT] &&
             mod->r_mat->rr_num->v[AG] != mod->r_mat->rr_num->v[AC] &&
             mod->r_mat->rr_num->v[CT] != mod->r_mat->rr_num->v[AC]) 
            {
              for(i=1;i<mod->r_mat->n_diff_rr;i++) mod->r_mat->rr_val->v[i] = 4.0;          
            }
          else if(mod->r_mat->n_diff_rr == 6) /* Custom <-> GTR model */
            {
              mod->r_mat->rr_val->v[AG] = 4.0;
              mod->r_mat->rr_val->v[CT] = 4.0;
            }
        }
      else if(mod->whichmodel == GTR)
        {
          For(i,6) mod->r_mat->rr_val->v[i] = 1.0;
          mod->r_mat->rr_val->v[AG] = 4.0;
          mod->r_mat->rr_val->v[CT] = 4.0;
        }
    }
  
  if(mod->s_opt->opt_alpha)  mod->ras->alpha->v  = 1.0;
  if(mod->s_opt->opt_pinvar) mod->ras->pinvar->v = 0.2;
  
  if(io->datatype == NT) /* Nucleotides */
    {
      /* init for nucleotides */
      mod->lambda->v    = 1.;
      /* mod->update_eigen = YES; */


      if(mod->whichmodel == JC69)
        {
          mod->e_frq->pi->v[0] = mod->e_frq->pi->v[1] = mod->e_frq->pi->v[2] = mod->e_frq->pi->v[3] = .25;
          mod->kappa->v = 1.;
          mod->s_opt->opt_state_freq = NO;
          mod->s_opt->opt_kappa      = NO;
          mod->s_opt->opt_lambda     = NO;
          mod->update_eigen          = NO;
        }
      
      if(mod->whichmodel == K80)
        {
          mod->e_frq->pi->v[0] = mod->e_frq->pi->v[1] = mod->e_frq->pi->v[2] = mod->e_frq->pi->v[3] = .25;
          mod->s_opt->opt_state_freq = NO;
          mod->s_opt->opt_lambda     = NO;
          mod->update_eigen          = NO;
        }
      
      if(mod->whichmodel == F81)
        {
          mod->kappa->v = 1.;
          mod->update_eigen = NO;
        }
      
      if(mod->whichmodel == F84)
        {
          aux = ((mod->e_frq->pi->v[0]+mod->e_frq->pi->v[2])-(mod->e_frq->pi->v[1]+mod->e_frq->pi->v[3]))/(2.*mod->kappa->v);
          mod->lambda->v = ((mod->e_frq->pi->v[1]+mod->e_frq->pi->v[3]) + aux)/((mod->e_frq->pi->v[0]+mod->e_frq->pi->v[2]) - aux);
          mod->update_eigen = NO;
        }
      
      if(mod->whichmodel == TN93)
        {
          mod->update_eigen = NO;
          if(io->mod->s_opt->opt_kappa) io->mod->s_opt->opt_lambda = YES;
        }
      
      if(mod->whichmodel == GTR)
        {
          mod->kappa->v          = 1.;
          mod->update_eigen      = YES;
          io->mod->s_opt->opt_rr = YES;
        }
      
      if(mod->whichmodel == CUSTOM)
        {
          mod->kappa->v = 1.;
          mod->update_eigen = YES;
          /* 	  io->mod->s_opt->opt_rr     = YES; */ /* What if the user decided not to optimise the rates? */
        }
      
      if(mod->whichmodel == GTR)
        {
          mod->custom_mod_string->s[0] = '0';
          mod->custom_mod_string->s[1] = '1';
          mod->custom_mod_string->s[2] = '2';
          mod->custom_mod_string->s[3] = '3';
          mod->custom_mod_string->s[4] = '4';
          mod->custom_mod_string->s[5] = '5';
          Translate_Custom_Mod_String(mod);
        }
      
      if(mod->s_opt->user_state_freq == YES && mod->whichmodel != JC69 && mod->whichmodel != K80)
        {
          For(i,4)
            {
              mod->e_frq->pi->v[i] = mod->user_b_freq->v[i];
            }
        }
      
      if(((mod->whichmodel == GTR)    ||
          (mod->whichmodel == CUSTOM) ||
          (mod->whichmodel == HKY85)) &&
         (mod->use_m4mod == NO))
        {
          mod->update_eigen = YES;
          Update_Eigen(mod);
          mod->update_eigen = NO;
        }
      /* if((mod->whichmodel != GTR)    &&  */
      /* 	 (mod->whichmodel != CUSTOM) &&  */
      /* 	 (mod->whichmodel != HKY85)) mod->update_eigen = NO; */
    }
  else if(mod->io->datatype == AA)
    {
      
      /* init for amino-acids */
      /* see comments of PMat_Empirical for details */
      /* read pi and Q from file */
      
      /* These initialisations are needed when analysing multiple
       * data sets
       */
      For(i,mod->ns*mod->ns) mod->r_mat->qmat->v[i] = .0;
      For(i,mod->ns        ) mod->e_frq->pi->v[i]   = .0;
      
      switch(mod->whichmodel)
        {
        case DAYHOFF :
          {
            Init_Qmat_Dayhoff(mod->r_mat->qmat->v,mod->e_frq->pi->v);
            if(mod->s_opt->opt_state_freq)
              For(i,mod->ns) mod->e_frq->pi->v[i] = data->b_frq[i];
            break;
          }
        case JTT :
          {
            Init_Qmat_JTT(mod->r_mat->qmat->v,mod->e_frq->pi->v);
            if(mod->s_opt->opt_state_freq)
              For(i,mod->ns) mod->e_frq->pi->v[i] = data->b_frq[i];
            break;
          }
        case MTREV :
          {
            Init_Qmat_MtREV(mod->r_mat->qmat->v,mod->e_frq->pi->v);
            if(mod->s_opt->opt_state_freq)
              For(i,mod->ns) mod->e_frq->pi->v[i] = data->b_frq[i];
            break;
          }
        case LG :
          {
            Init_Qmat_LG(mod->r_mat->qmat->v,mod->e_frq->pi->v);
            if(mod->s_opt->opt_state_freq)
              For(i,mod->ns) mod->e_frq->pi->v[i] = data->b_frq[i];
            break;
          }
        case WAG :
          {
            Init_Qmat_WAG(mod->r_mat->qmat->v,mod->e_frq->pi->v);
            if(mod->s_opt->opt_state_freq)
              For(i,mod->ns) mod->e_frq->pi->v[i] = data->b_frq[i];
            break;
          }
        case DCMUT :
          {
            Init_Qmat_DCMut(mod->r_mat->qmat->v,mod->e_frq->pi->v);
            if(mod->s_opt->opt_state_freq)
              For(i,mod->ns) mod->e_frq->pi->v[i] = data->b_frq[i];
            break;
          }
        case RTREV :
          {
            Init_Qmat_RtREV(mod->r_mat->qmat->v,mod->e_frq->pi->v);
            if(mod->s_opt->opt_state_freq)
              For(i,mod->ns) mod->e_frq->pi->v[i] = data->b_frq[i];
            break;
          }
        case CPREV :
          {
            Init_Qmat_CpREV(mod->r_mat->qmat->v,mod->e_frq->pi->v);
            if(mod->s_opt->opt_state_freq)
              For(i,mod->ns) mod->e_frq->pi->v[i] = data->b_frq[i];
            break;
          }
        case VT :
          {
            Init_Qmat_VT(mod->r_mat->qmat->v,mod->e_frq->pi->v);
            if(mod->s_opt->opt_state_freq)
              For(i,mod->ns) mod->e_frq->pi->v[i] = data->b_frq[i];
            break;
          }
        case BLOSUM62 :
          {
            Init_Qmat_Blosum62(mod->r_mat->qmat->v,mod->e_frq->pi->v);
            if(mod->s_opt->opt_state_freq)
              For(i,mod->ns) mod->e_frq->pi->v[i] = data->b_frq[i];
            break;
          }
        case MTMAM :
          {
            Init_Qmat_MtMam(mod->r_mat->qmat->v,mod->e_frq->pi->v);
            if(mod->s_opt->opt_state_freq)
              For(i,mod->ns) mod->e_frq->pi->v[i] = data->b_frq[i];
            break;
          }
        case MTART :
          {
            Init_Qmat_MtArt(mod->r_mat->qmat->v,mod->e_frq->pi->v);
            if(mod->s_opt->opt_state_freq)
              For(i,mod->ns) mod->e_frq->pi->v[i] = data->b_frq[i];
            break;
          }
        case HIVW :
          {
            Init_Qmat_HIVw(mod->r_mat->qmat->v,mod->e_frq->pi->v);
            if(mod->s_opt->opt_state_freq)
              For(i,mod->ns) mod->e_frq->pi->v[i] = data->b_frq[i];
            break;
          }
        case HIVB :
          {
            Init_Qmat_HIVb(mod->r_mat->qmat->v,mod->e_frq->pi->v);
            if(mod->s_opt->opt_state_freq)
              For(i,mod->ns) mod->e_frq->pi->v[i] = data->b_frq[i];
            break;
          }
         case AB :
          {
            Init_Qmat_AB(mod->r_mat->qmat->v,mod->e_frq->pi->v);
            if(mod->s_opt->opt_state_freq)
              For(i,mod->ns) mod->e_frq->pi->v[i] = data->b_frq[i];
            break;
          } 
        case CUSTOMAA :
          {
            Read_Qmat(mod->r_mat->qmat->v,mod->e_frq->pi->v,mod->fp_aa_rate_mat);
	    Print_Qmat_AA(mod->r_mat->qmat->v,mod->e_frq->pi->v);
            if(mod->s_opt->opt_state_freq) For(i,mod->ns) mod->e_frq->pi->v[i] = data->b_frq[i];
	    break;
	  }
	case FLU : 
	  {
	    Init_Qmat_FLU(mod->r_mat->qmat->v,mod->e_frq->pi->v);
	    if(mod->s_opt->opt_state_freq)
	      For(i,mod->ns) mod->e_frq->pi->v[i] = data->b_frq[i];
            break;
          }
        default :
          {
            Init_Qmat_LG(mod->r_mat->qmat->v,mod->e_frq->pi->v);
            if(mod->s_opt->opt_state_freq)
              For(i,mod->ns) mod->e_frq->pi->v[i] = data->b_frq[i];
            break;
          }
        }
      

      For(i,mod->ns) if(mod->e_frq->pi->v[i] < 1.E-10)
        {
          PhyML_Printf("\n. WARNING: at least one amino-acid frequency is equal to 0.0!");
          PhyML_Printf("\n. Numerical precision issues are likely to arise...");
          break;
        }
                       

      /*       /\* multiply the nth col of Q by the nth term of pi/100 just as in PAML *\/ */
      For(i,mod->ns) For(j,mod->ns) mod->r_mat->qmat->v[i*mod->ns+j] *= mod->e_frq->pi->v[j] / 100.0;
      
      /* compute diagonal terms of Q and mean rate mr = l/t */
      mod->mr->v= .0;
      For (i,mod->ns)
        {
          sum=.0;
          For(j, mod->ns) sum += mod->r_mat->qmat->v[i*mod->ns+j];
          mod->r_mat->qmat->v[i*mod->ns+i] = -sum;
          mod->mr->v += mod->e_frq->pi->v[i] * sum;
        }
      
      /* scale imod->nstantaneous rate matrix so that mu=1 */
      For (i,mod->ns*mod->ns) mod->r_mat->qmat->v[i] /= mod->mr->v;
      
      /* compute eigenvectors/values */
      result = 0;
      
      For(i,mod->ns*mod->ns) mod->r_mat->qmat_buff->v[i] = mod->r_mat->qmat->v[i];
      
      if(!Eigen(1,mod->r_mat->qmat_buff->v,mod->eigen->size,mod->eigen->e_val,
                mod->eigen->e_val_im,mod->eigen->r_e_vect,
                mod->eigen->r_e_vect_im,mod->eigen->space))
        {
          /* compute inverse(Vr) into Vi */
          For (i,mod->ns*mod->ns) mod->eigen->l_e_vect[i] = mod->eigen->r_e_vect[i];
          if(!Matinv(mod->eigen->l_e_vect,mod->eigen->size,mod->eigen->size,YES))
            {
              PhyML_Printf("\n== Err in file %s at line %d.",__FILE__,__LINE__);
              Exit("\n");
            }
          
          /* compute the diagonal terms of EXP(D) */
          For(i,mod->ns) mod->eigen->e_val[i] = (phydbl)EXP(mod->eigen->e_val[i]);
        }
      else
        {
          if (result==-1) PhyML_Printf("\n== Eigenvalues/vectors computation does not converge : computation cancelled");
          else if (result==1) PhyML_Printf("\n== Complex eigenvalues/vectors : computation cancelled");
        }
    }
  else if(mod->io->datatype == GENERIC)
    {
      /* Uniform state frequencies */
      For(i,mod->ns)  mod->e_frq->pi->v[i] = 1./(phydbl)mod->ns;
      mod->kappa->v = 1;
      mod->s_opt->opt_state_freq = NO;
      mod->s_opt->opt_kappa      = NO;
      mod->s_opt->opt_lambda     = NO;
      mod->update_eigen          = NO;
    }
  else
    {
      PhyML_Printf("\n== Datatype not recognized.\n");
      PhyML_Printf("\n== Err. in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  if(!mod->use_m4mod) Set_Model_Parameters(mod);

  Init_Eigen_Struct(mod->eigen);

  free(dr);free(di);free(space);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int Init_Qmat_Dayhoff(phydbl *daa, phydbl *pi)
{
  /* Dayhoff's model data
   * Dayhoff, M.O., Schwartz, R.M., Orcutt, B.C. (1978)
   * "A model of evolutionary change in proteins."
   * Dayhoff, M.O.(ed.) Atlas of Protein Sequence Structur., Vol5, Suppl3.
   * National Biomedical Research Foundation, Washington DC, pp.345-352.
   */
  int i,j,naa;

  naa = 20;

/*   PhyML_Printf("\n\n. REMINDER : THIS IS NOT DAYHOFF !!!\n\n"); */

/*   daa[1*20 + 0] = 0.538903;  */
/*   daa[2*20 + 0] = 0.412504;  */
/*   daa[2*20 + 1] = 0.736081;  */
/*   daa[3*20 + 0] = 0.586915;  */
/*   daa[3*20 + 1] = 0.108051;  */
/*   daa[3*20 + 2] = 5.446642;  */
/*   daa[4*20 + 0] = 2.189718;  */
/*   daa[4*20 + 1] = 0.830604;  */
/*   daa[4*20 + 2] = 0.573426;  */
/*   daa[4*20 + 3] = 0.077565;  */
/*   daa[5*20 + 0] = 1.08213;  */
/*   daa[5*20 + 1] = 2.950693;  */
/*   daa[5*20 + 2] = 1.739514;  */
/*   daa[5*20 + 3] = 0.559035;  */
/*   daa[5*20 + 4] = 0.111314;  */
/*   daa[6*20 + 0] = 1.386865;  */
/*   daa[6*20 + 1] = 0.358434;  */
/*   daa[6*20 + 2] = 0.541447;  */
/*   daa[6*20 + 3] = 5.406871;  */
/*   daa[6*20 + 4] = 0.003738;  */
/*   daa[6*20 + 5] = 4.124043;  */
/*   daa[7*20 + 0] = 2.085747;  */
/*   daa[7*20 + 1] = 0.453132;  */
/*   daa[7*20 + 2] = 1.815844;  */
/*   daa[7*20 + 3] = 1.08647;  */
/*   daa[7*20 + 4] = 0.526418;  */
/*   daa[7*20 + 5] = 0.347693;  */
/*   daa[7*20 + 6] = 0.438476;  */
/*   daa[8*20 + 0] = 0.407898;  */
/*   daa[8*20 + 1] = 2.689322;  */
/*   daa[8*20 + 2] = 5.386808;  */
/*   daa[8*20 + 3] = 0.884563;  */
/*   daa[8*20 + 4] = 0.658583;  */
/*   daa[8*20 + 5] = 4.588358;  */
/*   daa[8*20 + 6] = 0.386218;  */
/*   daa[8*20 + 7] = 0.368668;  */
/*   daa[9*20 + 0] = 0.10177;  */
/*   daa[9*20 + 1] = 0.104875;  */
/*   daa[9*20 + 2] = 0.16239;  */
/*   daa[9*20 + 3] = 0.011698;  */
/*   daa[9*20 + 4] = 0.253282;  */
/*   daa[9*20 + 5] = 0.083872;  */
/*   daa[9*20 + 6] = 0.041767;  */
/*   daa[9*20 + 7] = 0.009778;  */
/*   daa[9*20 + 8] = 0.1042;  */
/*   daa[10*20 + 0] = 0.248202;  */
/*   daa[10*20 + 1] = 0.375742;  */
/*   daa[10*20 + 2] = 0.093863;  */
/*   daa[10*20 + 3] = 0.019241;  */
/*   daa[10*20 + 4] = 0.572688;  */
/*   daa[10*20 + 5] = 0.703538;  */
/*   daa[10*20 + 6] = 0.071961;  */
/*   daa[10*20 + 7] = 0.03006;  */
/*   daa[10*20 + 8] = 0.418222;  */
/*   daa[10*20 + 9] = 3.702051;  */
/*   daa[11*20 + 0] = 0.652019;  */
/*   daa[11*20 + 1] = 5.940478;  */
/*   daa[11*20 + 2] = 2.352253;  */
/*   daa[11*20 + 3] = 0.231001;  */
/*   daa[11*20 + 4] = 0.027995;  */
/*   daa[11*20 + 5] = 3.646743;  */
/*   daa[11*20 + 6] = 1.507981;  */
/*   daa[11*20 + 7] = 0.331175;  */
/*   daa[11*20 + 8] = 0.698362;  */
/*   daa[11*20 + 9] = 0.140326;  */
/*   daa[11*20 + 10] = 0.171396;  */
/*   daa[12*20 + 0] = 0.694226;  */
/*   daa[12*20 + 1] = 0.419899;  */
/*   daa[12*20 + 2] = 0.326927;  */
/*   daa[12*20 + 3] = 0.039488;  */
/*   daa[12*20 + 4] = 0.844827;  */
/*   daa[12*20 + 5] = 1.394214;  */
/*   daa[12*20 + 6] = 0.133235;  */
/*   daa[12*20 + 7] = 0.085075;  */
/*   daa[12*20 + 8] = 0.347092;  */
/*   daa[12*20 + 9] = 4.051255;  */
/*   daa[12*20 + 10] = 6.650794;  */
/*   daa[12*20 + 11] = 0.617549;  */
/*   daa[13*20 + 0] = 0.155206;  */
/*   daa[13*20 + 1] = 0.057971;  */
/*   daa[13*20 + 2] = 0.09816;  */
/*   daa[13*20 + 3] = 0.020441;  */
/*   daa[13*20 + 4] = 0.904305;  */
/*   daa[13*20 + 5] = 0.052719;  */
/*   daa[13*20 + 6] = 0.0219;  */
/*   daa[13*20 + 7] = 0.046668;  */
/*   daa[13*20 + 8] = 0.890005;  */
/*   daa[13*20 + 9] = 0.844963;  */
/*   daa[13*20 + 10] = 2.348881;  */
/*   daa[13*20 + 11] = 0.028372;  */
/*   daa[13*20 + 12] = 1.671635;  */
/*   daa[14*20 + 0] = 1.433475;  */
/*   daa[14*20 + 1] = 0.328393;  */
/*   daa[14*20 + 2] = 0.173181;  */
/*   daa[14*20 + 3] = 0.431874;  */
/*   daa[14*20 + 4] = 0.09902;  */
/*   daa[14*20 + 5] = 0.592324;  */
/*   daa[14*20 + 6] = 0.488352;  */
/*   daa[14*20 + 7] = 0.23865;  */
/*   daa[14*20 + 8] = 0.462856;  */
/*   daa[14*20 + 9] = 0.057048;  */
/*   daa[14*20 + 10] = 0.233532;  */
/*   daa[14*20 + 11] = 0.387808;  */
/*   daa[14*20 + 12] = 0.096377;  */
/*   daa[14*20 + 13] = 0.079912;  */
/*   daa[15*20 + 0] = 4.887126;  */
/*   daa[15*20 + 1] = 0.883923;  */
/*   daa[15*20 + 2] = 4.627163;  */
/*   daa[15*20 + 3] = 1.122164;  */
/*   daa[15*20 + 4] = 3.186667;  */
/*   daa[15*20 + 5] = 1.085947;  */
/*   daa[15*20 + 6] = 0.569339;  */
/*   daa[15*20 + 7] = 1.993432;  */
/*   daa[15*20 + 8] = 0.867972;  */
/*   daa[15*20 + 9] = 0.070512;  */
/*   daa[15*20 + 10] = 0.163009;  */
/*   daa[15*20 + 11] = 0.718913;  */
/*   daa[15*20 + 12] = 0.301103;  */
/*   daa[15*20 + 13] = 0.32579;  */
/*   daa[15*20 + 14] = 1.449582;  */
/*   daa[16*20 + 0] = 2.030538;  */
/*   daa[16*20 + 1] = 0.639463;  */
/*   daa[16*20 + 2] = 2.076294;  */
/*   daa[16*20 + 3] = 0.377239;  */
/*   daa[16*20 + 4] = 1.42848;  */
/*   daa[16*20 + 5] = 0.979403;  */
/*   daa[16*20 + 6] = 0.647562;  */
/*   daa[16*20 + 7] = 0.145556;  */
/*   daa[16*20 + 8] = 0.493329;  */
/*   daa[16*20 + 9] = 0.973405;  */
/*   daa[16*20 + 10] = 0.271824;  */
/*   daa[16*20 + 11] = 1.20033;  */
/*   daa[16*20 + 12] = 1.659187;  */
/*   daa[16*20 + 13] = 0.1217;  */
/*   daa[16*20 + 14] = 0.571399;  */
/*   daa[16*20 + 15] = 6.641034;  */
/*   daa[17*20 + 0] = 0.131405;  */
/*   daa[17*20 + 1] = 0.552911;  */
/*   daa[17*20 + 2] = 0.079985;  */
/*   daa[17*20 + 3] = 0.060514;  */
/*   daa[17*20 + 4] = 0.633662;  */
/*   daa[17*20 + 5] = 0.21823;  */
/*   daa[17*20 + 6] = 0.074988;  */
/*   daa[17*20 + 7] = 0.169114;  */
/*   daa[17*20 + 8] = 0.847725;  */
/*   daa[17*20 + 9] = 0.10627;  */
/*   daa[17*20 + 10] = 0.622044;  */
/*   daa[17*20 + 11] = 0.060755;  */
/*   daa[17*20 + 12] = 0.719575;  */
/*   daa[17*20 + 13] = 3.14824;  */
/*   daa[17*20 + 14] = 0.077123;  */
/*   daa[17*20 + 15] = 0.276716;  */
/*   daa[17*20 + 16] = 0.148883;  */
/*   daa[18*20 + 0] = 0.165179;  */
/*   daa[18*20 + 1] = 0.224883;  */
/*   daa[18*20 + 2] = 0.528334;  */
/*   daa[18*20 + 3] = 0.121252;  */
/*   daa[18*20 + 4] = 1.174118;  */
/*   daa[18*20 + 5] = 0.177062;  */
/*   daa[18*20 + 6] = 0.074715;  */
/*   daa[18*20 + 7] = 0.042356;  */
/*   daa[18*20 + 8] = 5.911187;  */
/*   daa[18*20 + 9] = 0.192481;  */
/*   daa[18*20 + 10] = 0.321454;  */
/*   daa[18*20 + 11] = 0.090556;  */
/*   daa[18*20 + 12] = 0.406415;  */
/*   daa[18*20 + 13] = 10.908861;  */
/*   daa[18*20 + 14] = 0.070752;  */
/*   daa[18*20 + 15] = 0.328483;  */
/*   daa[18*20 + 16] = 0.181539;  */
/*   daa[18*20 + 17] = 3.823886;  */
/*   daa[19*20 + 0] = 1.78517;  */
/*   daa[19*20 + 1] = 0.166975;  */
/*   daa[19*20 + 2] = 0.106482;  */
/*   daa[19*20 + 3] = 0.041707;  */
/*   daa[19*20 + 4] = 1.876812;  */
/*   daa[19*20 + 5] = 0.22421;  */
/*   daa[19*20 + 6] = 0.247356;  */
/*   daa[19*20 + 7] = 0.06688;  */
/*   daa[19*20 + 8] = 0.1436;  */
/*   daa[19*20 + 9] = 9.60184;  */
/*   daa[19*20 + 10] = 1.599119;  */
/*   daa[19*20 + 11] = 0.17319;  */
/*   daa[19*20 + 12] = 1.645134;  */
/*   daa[19*20 + 13] = 0.438571;  */
/*   daa[19*20 + 14] = 0.252486;  */
/*   daa[19*20 + 15] = 0.105536;  */
/*   daa[19*20 + 16] = 1.789097;  */
/*   daa[19*20 + 17] = 0.147951;  */
/*   daa[19*20 + 18] = 0.200571; */



/*   for (i=0; i<naa; i++)  for (j=0; j<i; j++)  daa[j*naa+i] = daa[i*naa+j]; */


/*   pi[0] = 0.093862;  */
/*   pi[1] = 0.05757;  */
/*   pi[2] = 0.037017;  */
/*   pi[3] = 0.060952;  */
/*   pi[4] = 0.011004;  */
/*   pi[5] = 0.044631;  */
/*   pi[6] = 0.083648;  */
/*   pi[7] = 0.053004;  */
/*   pi[8] = 0.022528;  */
/*   pi[9] = 0.060152;  */
/*   pi[10] = 0.092314;  */
/*   pi[11] = 0.065886;  */
/*   pi[12] = 0.021075;  */
/*   pi[13] = 0.033928;  */
/*   pi[14] = 0.043333;  */
/*   pi[15] = 0.052957;  */
/*   pi[16] = 0.053465;  */
/*   pi[17] = 0.00928;  */
/*   pi[18] = 0.030025;  */
/*   pi[19] = 0.073369; */


  daa[ 1*20+ 0] =   27.00; daa[ 2*20+ 0] =   98.00; daa[ 2*20+ 1] =   32.00; daa[ 3*20+ 0] =  120.00;
  daa[ 3*20+ 1] =    0.00; daa[ 3*20+ 2] =  905.00; daa[ 4*20+ 0] =   36.00; daa[ 4*20+ 1] =   23.00;
  daa[ 4*20+ 2] =    0.00; daa[ 4*20+ 3] =    0.00; daa[ 5*20+ 0] =   89.00; daa[ 5*20+ 1] =  246.00;
  daa[ 5*20+ 2] =  103.00; daa[ 5*20+ 3] =  134.00; daa[ 5*20+ 4] =    0.00; daa[ 6*20+ 0] =  198.00;
  daa[ 6*20+ 1] =    1.00; daa[ 6*20+ 2] =  148.00; daa[ 6*20+ 3] = 1153.00; daa[ 6*20+ 4] =    0.00;
  daa[ 6*20+ 5] =  716.00; daa[ 7*20+ 0] =  240.00; daa[ 7*20+ 1] =    9.00; daa[ 7*20+ 2] =  139.00;
  daa[ 7*20+ 3] =  125.00; daa[ 7*20+ 4] =   11.00; daa[ 7*20+ 5] =   28.00; daa[ 7*20+ 6] =   81.00;
  daa[ 8*20+ 0] =   23.00; daa[ 8*20+ 1] =  240.00; daa[ 8*20+ 2] =  535.00; daa[ 8*20+ 3] =   86.00;
  daa[ 8*20+ 4] =   28.00; daa[ 8*20+ 5] =  606.00; daa[ 8*20+ 6] =   43.00; daa[ 8*20+ 7] =   10.00;
  daa[ 9*20+ 0] =   65.00; daa[ 9*20+ 1] =   64.00; daa[ 9*20+ 2] =   77.00; daa[ 9*20+ 3] =   24.00;
  daa[ 9*20+ 4] =   44.00; daa[ 9*20+ 5] =   18.00; daa[ 9*20+ 6] =   61.00; daa[ 9*20+ 7] =    0.00;
  daa[ 9*20+ 8] =    7.00; daa[10*20+ 0] =   41.00; daa[10*20+ 1] =   15.00; daa[10*20+ 2] =   34.00;
  daa[10*20+ 3] =    0.00; daa[10*20+ 4] =    0.00; daa[10*20+ 5] =   73.00; daa[10*20+ 6] =   11.00;
  daa[10*20+ 7] =    7.00; daa[10*20+ 8] =   44.00; daa[10*20+ 9] =  257.00; daa[11*20+ 0] =   26.00;
  daa[11*20+ 1] =  464.00; daa[11*20+ 2] =  318.00; daa[11*20+ 3] =   71.00; daa[11*20+ 4] =    0.00;
  daa[11*20+ 5] =  153.00; daa[11*20+ 6] =   83.00; daa[11*20+ 7] =   27.00; daa[11*20+ 8] =   26.00;
  daa[11*20+ 9] =   46.00; daa[11*20+10] =   18.00; daa[12*20+ 0] =   72.00; daa[12*20+ 1] =   90.00;
  daa[12*20+ 2] =    1.00; daa[12*20+ 3] =    0.00; daa[12*20+ 4] =    0.00; daa[12*20+ 5] =  114.00;
  daa[12*20+ 6] =   30.00; daa[12*20+ 7] =   17.00; daa[12*20+ 8] =    0.00; daa[12*20+ 9] =  336.00;
  daa[12*20+10] =  527.00; daa[12*20+11] =  243.00; daa[13*20+ 0] =   18.00; daa[13*20+ 1] =   14.00;
  daa[13*20+ 2] =   14.00; daa[13*20+ 3] =    0.00; daa[13*20+ 4] =    0.00; daa[13*20+ 5] =    0.00;
  daa[13*20+ 6] =    0.00; daa[13*20+ 7] =   15.00; daa[13*20+ 8] =   48.00; daa[13*20+ 9] =  196.00;
  daa[13*20+10] =  157.00; daa[13*20+11] =    0.00; daa[13*20+12] =   92.00; daa[14*20+ 0] =  250.00;
  daa[14*20+ 1] =  103.00; daa[14*20+ 2] =   42.00; daa[14*20+ 3] =   13.00; daa[14*20+ 4] =   19.00;
  daa[14*20+ 5] =  153.00; daa[14*20+ 6] =   51.00; daa[14*20+ 7] =   34.00; daa[14*20+ 8] =   94.00;
  daa[14*20+ 9] =   12.00; daa[14*20+10] =   32.00; daa[14*20+11] =   33.00; daa[14*20+12] =   17.00;
  daa[14*20+13] =   11.00; daa[15*20+ 0] =  409.00; daa[15*20+ 1] =  154.00; daa[15*20+ 2] =  495.00;
  daa[15*20+ 3] =   95.00; daa[15*20+ 4] =  161.00; daa[15*20+ 5] =   56.00; daa[15*20+ 6] =   79.00;
  daa[15*20+ 7] =  234.00; daa[15*20+ 8] =   35.00; daa[15*20+ 9] =   24.00; daa[15*20+10] =   17.00;
  daa[15*20+11] =   96.00; daa[15*20+12] =   62.00; daa[15*20+13] =   46.00; daa[15*20+14] =  245.00;
  daa[16*20+ 0] =  371.00; daa[16*20+ 1] =   26.00; daa[16*20+ 2] =  229.00; daa[16*20+ 3] =   66.00;
  daa[16*20+ 4] =   16.00; daa[16*20+ 5] =   53.00; daa[16*20+ 6] =   34.00; daa[16*20+ 7] =   30.00;
  daa[16*20+ 8] =   22.00; daa[16*20+ 9] =  192.00; daa[16*20+10] =   33.00; daa[16*20+11] =  136.00;
  daa[16*20+12] =  104.00; daa[16*20+13] =   13.00; daa[16*20+14] =   78.00; daa[16*20+15] =  550.00;
  daa[17*20+ 0] =    0.00; daa[17*20+ 1] =  201.00; daa[17*20+ 2] =   23.00; daa[17*20+ 3] =    0.00;
  daa[17*20+ 4] =    0.00; daa[17*20+ 5] =    0.00; daa[17*20+ 6] =    0.00; daa[17*20+ 7] =    0.00;
  daa[17*20+ 8] =   27.00; daa[17*20+ 9] =    0.00; daa[17*20+10] =   46.00; daa[17*20+11] =    0.00;
  daa[17*20+12] =    0.00; daa[17*20+13] =   76.00; daa[17*20+14] =    0.00; daa[17*20+15] =   75.00;
  daa[17*20+16] =    0.00; daa[18*20+ 0] =   24.00; daa[18*20+ 1] =    8.00; daa[18*20+ 2] =   95.00;
  daa[18*20+ 3] =    0.00; daa[18*20+ 4] =   96.00; daa[18*20+ 5] =    0.00; daa[18*20+ 6] =   22.00;
  daa[18*20+ 7] =    0.00; daa[18*20+ 8] =  127.00; daa[18*20+ 9] =   37.00; daa[18*20+10] =   28.00;
  daa[18*20+11] =   13.00; daa[18*20+12] =    0.00; daa[18*20+13] =  698.00; daa[18*20+14] =    0.00;
  daa[18*20+15] =   34.00; daa[18*20+16] =   42.00; daa[18*20+17] =   61.00; daa[19*20+ 0] =  208.00;
  daa[19*20+ 1] =   24.00; daa[19*20+ 2] =   15.00; daa[19*20+ 3] =   18.00; daa[19*20+ 4] =   49.00;
  daa[19*20+ 5] =   35.00; daa[19*20+ 6] =   37.00; daa[19*20+ 7] =   54.00; daa[19*20+ 8] =   44.00;
  daa[19*20+ 9] =  889.00; daa[19*20+10] =  175.00; daa[19*20+11] =   10.00; daa[19*20+12] =  258.00;
  daa[19*20+13] =   12.00; daa[19*20+14] =   48.00; daa[19*20+15] =   30.00; daa[19*20+16] =  157.00;
  daa[19*20+17] =    0.00; daa[19*20+18] =   28.00;

  for (i=0; i<naa; i++)  for (j=0; j<i; j++)  daa[j*naa+i] = daa[i*naa+j];

  pi[ 0] = 0.087127; pi[ 1] = 0.040904; pi[ 2] = 0.040432; pi[ 3] = 0.046872;
  pi[ 4] = 0.033474; pi[ 5] = 0.038255; pi[ 6] = 0.049530; pi[ 7] = 0.088612;
  pi[ 8] = 0.033618; pi[ 9] = 0.036886; pi[10] = 0.085357; pi[11] = 0.080482;
  pi[12] = 0.014753; pi[13] = 0.039772; pi[14] = 0.050680; pi[15] = 0.069577;
  pi[16] = 0.058542; pi[17] = 0.010494; pi[18] = 0.029916; pi[19] = 0.064718;



  return 1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


int Init_Qmat_DCMut(phydbl *daa, phydbl *pi)
{
  /*
     DCMut : new implementation based on Dayhoff et al.'s raw data and amino acid mutabilities
     C. Kosiol and N. Goldman. (2005),
     ``Different versions of the Dayhoff rate matrix'',
     Mol. Biol. Evol., 22. 193-199.
  */

  int i,j,naa;

  naa = 20;

  daa[ 1*20+ 0] =   26.78280; daa[ 2*20+ 0] =   98.44740; daa[ 2*20+ 1] =   32.70590; daa[ 3*20+ 0] =  119.98050;
  daa[ 3*20+ 1] =    0.00000; daa[ 3*20+ 2] =  893.15150; daa[ 4*20+ 0] =   36.00160; daa[ 4*20+ 1] =   23.23740;
  daa[ 4*20+ 2] =    0.00000; daa[ 4*20+ 3] =    0.00000; daa[ 5*20+ 0] =   88.77530; daa[ 5*20+ 1] =  243.99390;
  daa[ 5*20+ 2] =  102.85090; daa[ 5*20+ 3] =  134.85510; daa[ 5*20+ 4] =    0.00000; daa[ 6*20+ 0] =  196.11670;
  daa[ 6*20+ 1] =    0.00000; daa[ 6*20+ 2] =  149.34090; daa[ 6*20+ 3] = 1138.86590; daa[ 6*20+ 4] =    0.00000;
  daa[ 6*20+ 5] =  708.60220; daa[ 7*20+ 0] =  238.61110; daa[ 7*20+ 1] =    8.77910; daa[ 7*20+ 2] =  138.53520;
  daa[ 7*20+ 3] =  124.09810; daa[ 7*20+ 4] =   10.72780; daa[ 7*20+ 5] =   28.15810; daa[ 7*20+ 6] =   81.19070;
  daa[ 8*20+ 0] =   22.81160; daa[ 8*20+ 1] =  238.31480; daa[ 8*20+ 2] =  529.00240; daa[ 8*20+ 3] =   86.82410;
  daa[ 8*20+ 4] =   28.27290; daa[ 8*20+ 5] =  601.16130; daa[ 8*20+ 6] =   43.94690; daa[ 8*20+ 7] =   10.68020;
  daa[ 9*20+ 0] =   65.34160; daa[ 9*20+ 1] =   63.26290; daa[ 9*20+ 2] =   76.80240; daa[ 9*20+ 3] =   23.92480;
  daa[ 9*20+ 4] =   43.80740; daa[ 9*20+ 5] =   18.03930; daa[ 9*20+ 6] =   60.95260; daa[ 9*20+ 7] =    0.00000;
  daa[ 9*20+ 8] =    7.69810; daa[10*20+ 0] =   40.64310; daa[10*20+ 1] =   15.49240; daa[10*20+ 2] =   34.11130;
  daa[10*20+ 3] =    0.00000; daa[10*20+ 4] =    0.00000; daa[10*20+ 5] =   73.07720; daa[10*20+ 6] =   11.28800;
  daa[10*20+ 7] =    7.15140; daa[10*20+ 8] =   44.35040; daa[10*20+ 9] =  255.66850; daa[11*20+ 0] =   25.86350;
  daa[11*20+ 1] =  461.01240; daa[11*20+ 2] =  314.83710; daa[11*20+ 3] =   71.69130; daa[11*20+ 4] =    0.00000;
  daa[11*20+ 5] =  151.90780; daa[11*20+ 6] =   83.00780; daa[11*20+ 7] =   26.76830; daa[11*20+ 8] =   27.04750;
  daa[11*20+ 9] =   46.08570; daa[11*20+10] =   18.06290; daa[12*20+ 0] =   71.78400; daa[12*20+ 1] =   89.63210;
  daa[12*20+ 2] =    0.00000; daa[12*20+ 3] =    0.00000; daa[12*20+ 4] =    0.00000; daa[12*20+ 5] =  112.74990;
  daa[12*20+ 6] =   30.48030; daa[12*20+ 7] =   17.03720; daa[12*20+ 8] =    0.00000; daa[12*20+ 9] =  333.27320;
  daa[12*20+10] =  523.01150; daa[12*20+11] =  241.17390; daa[13*20+ 0] =   18.36410; daa[13*20+ 1] =   13.69060;
  daa[13*20+ 2] =   13.85030; daa[13*20+ 3] =    0.00000; daa[13*20+ 4] =    0.00000; daa[13*20+ 5] =    0.00000;
  daa[13*20+ 6] =    0.00000; daa[13*20+ 7] =   15.34780; daa[13*20+ 8] =   47.59270; daa[13*20+ 9] =  195.19510;
  daa[13*20+10] =  156.51600; daa[13*20+11] =    0.00000; daa[13*20+12] =   92.18600; daa[14*20+ 0] =  248.59200;
  daa[14*20+ 1] =  102.83130; daa[14*20+ 2] =   41.92440; daa[14*20+ 3] =   13.39400; daa[14*20+ 4] =   18.75500;
  daa[14*20+ 5] =  152.61880; daa[14*20+ 6] =   50.70030; daa[14*20+ 7] =   34.71530; daa[14*20+ 8] =   93.37090;
  daa[14*20+ 9] =   11.91520; daa[14*20+10] =   31.62580; daa[14*20+11] =   33.54190; daa[14*20+12] =   17.02050;
  daa[14*20+13] =   11.05060; daa[15*20+ 0] =  405.18700; daa[15*20+ 1] =  153.15900; daa[15*20+ 2] =  488.58920;
  daa[15*20+ 3] =   95.60970; daa[15*20+ 4] =  159.83560; daa[15*20+ 5] =   56.18280; daa[15*20+ 6] =   79.39990;
  daa[15*20+ 7] =  232.22430; daa[15*20+ 8] =   35.36430; daa[15*20+ 9] =   24.79550; daa[15*20+10] =   17.14320;
  daa[15*20+11] =   95.45570; daa[15*20+12] =   61.99510; daa[15*20+13] =   45.99010; daa[15*20+14] =  242.72020;
  daa[16*20+ 0] =  368.03650; daa[16*20+ 1] =   26.57450; daa[16*20+ 2] =  227.16970; daa[16*20+ 3] =   66.09300;
  daa[16*20+ 4] =   16.23660; daa[16*20+ 5] =   52.56510; daa[16*20+ 6] =   34.01560; daa[16*20+ 7] =   30.66620;
  daa[16*20+ 8] =   22.63330; daa[16*20+ 9] =  190.07390; daa[16*20+10] =   33.10900; daa[16*20+11] =  135.05990;
  daa[16*20+12] =  103.15340; daa[16*20+13] =   13.66550; daa[16*20+14] =   78.28570; daa[16*20+15] =  543.66740;
  daa[17*20+ 0] =    0.00000; daa[17*20+ 1] =  200.13750; daa[17*20+ 2] =   22.49680; daa[17*20+ 3] =    0.00000;
  daa[17*20+ 4] =    0.00000; daa[17*20+ 5] =    0.00000; daa[17*20+ 6] =    0.00000; daa[17*20+ 7] =    0.00000;
  daa[17*20+ 8] =   27.05640; daa[17*20+ 9] =    0.00000; daa[17*20+10] =   46.17760; daa[17*20+11] =    0.00000;
  daa[17*20+12] =    0.00000; daa[17*20+13] =   76.23540; daa[17*20+14] =    0.00000; daa[17*20+15] =   74.08190;
  daa[17*20+16] =    0.00000; daa[18*20+ 0] =   24.41390; daa[18*20+ 1] =    7.80120; daa[18*20+ 2] =   94.69400;
  daa[18*20+ 3] =    0.00000; daa[18*20+ 4] =   95.31640; daa[18*20+ 5] =    0.00000; daa[18*20+ 6] =   21.47170;
  daa[18*20+ 7] =    0.00000; daa[18*20+ 8] =  126.54000; daa[18*20+ 9] =   37.48340; daa[18*20+10] =   28.65720;
  daa[18*20+11] =   13.21420; daa[18*20+12] =    0.00000; daa[18*20+13] =  695.26290; daa[18*20+14] =    0.00000;
  daa[18*20+15] =   33.62890; daa[18*20+16] =   41.78390; daa[18*20+17] =   60.80700; daa[19*20+ 0] =  205.95640;
  daa[19*20+ 1] =   24.03680; daa[19*20+ 2] =   15.80670; daa[19*20+ 3] =   17.83160; daa[19*20+ 4] =   48.46780;
  daa[19*20+ 5] =   34.69830; daa[19*20+ 6] =   36.72500; daa[19*20+ 7] =   53.81650; daa[19*20+ 8] =   43.87150;
  daa[19*20+ 9] =  881.00380; daa[19*20+10] =  174.51560; daa[19*20+11] =   10.38500; daa[19*20+12] =  256.59550;
  daa[19*20+13] =   12.36060; daa[19*20+14] =   48.50260; daa[19*20+15] =   30.38360; daa[19*20+16] =  156.19970;
  daa[19*20+17] =    0.00000; daa[19*20+18] =   27.93790;

  for (i=0; i<naa; i++)  for (j=0; j<i; j++)  daa[j*naa+i] = daa[i*naa+j];

  pi[ 0] = 0.087127; pi[ 1] = 0.040904; pi[ 2] = 0.040432; pi[ 3] = 0.046872;
  pi[ 4] = 0.033474; pi[ 5] = 0.038255; pi[ 6] = 0.049530; pi[ 7] = 0.088612;
  pi[ 8] = 0.033619; pi[ 9] = 0.036886; pi[10] = 0.085357; pi[11] = 0.080481;
  pi[12] = 0.014753; pi[13] = 0.039772; pi[14] = 0.050680; pi[15] = 0.069577;
  pi[16] = 0.058542; pi[17] = 0.010494; pi[18] = 0.029916; pi[19] = 0.064718;

  return 1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


int Init_Qmat_MtArt(phydbl *daa, phydbl *pi) // Added by Federico Abascal
{
    /*
       Federico Abascal, April 2005 (c).

       This model has been derived from 36 artropoda mitochondrial genomes.

       Each gene of the given species was aligned individually. Then, alignments of the whole set
       of 13 genes where concatenated and passed through GBlocks (Castresana, 2000, in JME) with
       parameters and output:

            Minimum Number Of Sequences For A Conserved Position: 20
            Minimum Number Of Sequences For A Flanking Position: 32
            Maximum Number Of Contiguous Nonconserved Positions: 8
            Minimum Length Of A Block: 10
            Allowed Gap Positions: With Half
            Use Similarity Matrices: Yes

            Flank positions of the 40 selected block(s)
            Flanks: [6  22]  [26  44]  [61  70]  [77  143]  [145  185]  [208  236]  [309  640]
            [644  802]  [831  941]  [956  966]  [973  1062]  [1085  1339]  [1343  1702]
            [1754  1831]  [1840  1911]  [1916  1987]  [2011  2038]  [2097  2118]  [2125  2143]
            [2179  2215]  [2243  2268]  [2277  2288]  [2333  2347]  [2476  2518]  [2539  2558]
            [2600  2613]  [2637  2672]  [2738  2759]  [2784  2839]  [2882  2924]  [2948  3097]
            [3113  3123]  [3210  3235]  [3239  3322]  [3348  3392]  [3406  3526]  [3588  3617]
            [3660  3692]  [3803  3830]  [3909  3928]

            New number of positions in MtArt-strict.phy.fasta-gb: <b> 2664 </b> (67% of the original 3933 positions)

       The species included in the analysis were:
            Harpiosquilla harpax          [NCBI_TaxID 287944]
            Ixodes uriae                  [NCBI_TaxID 59655]
            Heptathela hangzhouensis      [NCBI_TaxID 216259]
            Triops longicaudatus          [NCBI_TaxID 58777]
            Gryllotalpa orientalis        [NCBI_TaxID 213494]
            lepidopsocid RS-2001          [NCBI_TaxID 159971]
            Locusta migratoria            [NCBI_TaxID 7004]
            Drosophila yakuba             [NCBI_TaxID 7245]
            Ostrinia furnacalis           [NCBI_TaxID 93504]
            Megabalanus volcano           [NCBI_TaxID 266495]
            Periplaneta fuliginosa        [NCBI_TaxID 36977]
            Thermobia domestica           [NCBI_TaxID 89055]
            Aleurochiton aceris           [NCBI_TaxID 266942]
            Schizaphis graminum           [NCBI_TaxID 13262]
            Pteronarcys princeps          [NCBI_TaxID 285953]
            Aleurodicus dugesii           [NCBI_TaxID 30099]
            Pollicipes polymerus          [NCBI_TaxID 36137]
            Gomphiocephalus hodgsoni      [NCBI_TaxID 221270]
            Habronattus oregonensis       [NCBI_TaxID 130930]
            Speleonectes tulumensis       [NCBI_TaxID 84346]
            Hutchinsoniella macracantha   [NCBI_TaxID 84335]
            Haemaphysalis flava           [NCBI_TaxID 181088]
            Scutigera coleoptrata         [NCBI_TaxID 29022]
            Vargula hilgendorfii          [NCBI_TaxID 6674]
            Tricholepidion gertschi       [NCBI_TaxID 89825]
            Varroa destructor             [NCBI_TaxID 109461]
            Bombyx mandarina              [NCBI_TaxID 7092]
            Thyropygus sp.                [NCBI_TaxID 174155]
            Tribolium castaneum           [NCBI_TaxID 7070]
            Pagurus longicarpus           [NCBI_TaxID 111067]
            Limulus polyphemus            [NCBI_TaxID 6850]
            Tetrodontophora bielanensis   [NCBI_TaxID 48717]
            Penaeus monodon               [NCBI_TaxID 6687]
            Daphnia pulex                 [NCBI_TaxID 6669]
            Apis mellifera                [NCBI_TaxID 7469]
            Anopheles gambiae             [NCBI_TaxID 7165]

       The topology used for inferring the model was:
            (((Daph_pulex,Trio_longi),((((((Aleu_aceri,Aleu_duges),Schi_grami),lepi_RS_20),
            ((((Ostr_furna,Bomb_manda),(Dros_yakub,Anop_gambi)),Apis_melli),Trib_casta)),
            ((Gryl_orien,Locu_migra),(Pter_princ,Peri_fulig))),(Tric_gerts,Ther_domes)),
            (Scut_coleo,Thyr_sp),Varg_hilge,Hutc_macra,((((Ixod_uriae,Haem_flava),Varr_destr),
            (Habr_orego,Hept_hangz)),Limu_polyp),(Poll_polym,Mega_volca),(Gomp_hodgs,Tetr_biela),
            ((Pagu_longi,Pena_monod),Harp_harpa),Spel_tulum));

            Note this is not the ML topoLOGy but the consensus one (based on morphoLOGical data,
            phyLOGenetic reconstruction using nuclear genes, etc). Where relationships are
            not clear, a polytomy was introduced (it contains quite a lot of polytomies!).

       The model was estimated using (the great and helpful) Ziheng Yang's Paml software package.
       A four-categorized gamma distribution was used to account for heterogeneity (alpha
       was estimated to be 0.47821). Sites with ambiguity data were taken into account.

       If you would like the data related to this matrix, please contact fabascal@uvigo.es.
       Federico Abascal (c)2005.
    */

    int i,j,naa;
    naa = 20;
    daa[1*20+ 0] = 0.2;     daa[2*20+ 0] = 0.2;     daa[2*20+ 1] = 0.2;     daa[3*20+ 0] = 0.6;
    daa[3*20+ 1] = 4.3;     daa[3*20+ 2] = 500.2;   daa[4*20+ 0] = 253.5;   daa[4*20+ 1] = 35.5;
    daa[4*20+ 2] = 98.2;    daa[4*20+ 3] = 10.6;    daa[5*20+ 0] = 0.2;     daa[5*20+ 1] = 154.0;
    daa[5*20+ 2] = 261.8;   daa[5*20+ 3] = 0.2;     daa[5*20+ 4] = 0.2;     daa[6*20+ 0] = 0.2;
    daa[6*20+ 1] = 0.2;     daa[6*20+ 2] = 183.0;   daa[6*20+ 3] = 861.8;   daa[6*20+ 4] = 0.2;
    daa[6*20+ 5] = 261.6;   daa[7*20+ 0] = 199.8;   daa[7*20+ 1] = 0.2;     daa[7*20+ 2] = 120.5;
    daa[7*20+ 3] = 12.5;    daa[7*20+ 4] = 80.5;    daa[7*20+ 5] = 2.6;     daa[7*20+ 6] = 43.9;
    daa[8*20+ 0] = 0.2;     daa[8*20+ 1] = 41.3;    daa[8*20+ 2] = 179.5;   daa[8*20+ 3] = 0.2;
    daa[8*20+ 4] = 12.4;    daa[8*20+ 5] = 313.5;   daa[8*20+ 6] = 15.2;    daa[8*20+ 7] = 0.2;
    daa[9*20+ 0] = 25.7;    daa[9*20+ 1] = 1.8;     daa[9*20+ 2] = 21.3;    daa[9*20+ 3] = 6.6;
    daa[9*20+ 4] = 63.0;    daa[9*20+ 5] = 10.5;    daa[9*20+ 6] = 6.8;     daa[9*20+ 7] = 2.7;
    daa[9*20+ 8] = 0.2;     daa[10*20+ 0] = 3.7;    daa[10*20+ 1] = 1.8;    daa[10*20+ 2] = 12.6;
    daa[10*20+ 3] = 1.2;    daa[10*20+ 4] = 78.7;   daa[10*20+ 5] = 16.3;   daa[10*20+ 6] = 1.7;
    daa[10*20+ 7] = 1.4;    daa[10*20+ 8] = 5.5;    daa[10*20+ 9] = 514.5;  daa[11*20+ 0] = 0.2;
    daa[11*20+ 1] = 208.6;  daa[11*20+ 2] = 467.3;  daa[11*20+ 3] = 1.7;    daa[11*20+ 4] = 0.2;
    daa[11*20+ 5] = 349.3;  daa[11*20+ 6] = 106.3;  daa[11*20+ 7] = 0.2;    daa[11*20+ 8] = 0.2;
    daa[11*20+ 9] = 3.5;    daa[11*20+ 10] = 3.8;   daa[12*20+ 0] = 120.6;  daa[12*20+ 1] = 5.2;
    daa[12*20+ 2] = 78.8;   daa[12*20+ 3] = 0.2;    daa[12*20+ 4] = 312.3;  daa[12*20+ 5] = 67.3;
    daa[12*20+ 6] = 0.2;    daa[12*20+ 7] = 55.7;   daa[12*20+ 8] = 0.2;    daa[12*20+ 9] = 514.8;
    daa[12*20+ 10] = 885.5; daa[12*20+ 11] = 105.6; daa[13*20+ 0] = 13.1;   daa[13*20+ 1] = 4.7;
    daa[13*20+ 2] = 19.7;   daa[13*20+ 3] = 0.2;    daa[13*20+ 4] = 184.1;  daa[13*20+ 5] = 0.2;
    daa[13*20+ 6] = 0.2;    daa[13*20+ 7] = 0.8;    daa[13*20+ 8] = 13.8;   daa[13*20+ 9] = 117.9;
    daa[13*20+ 10] = 262.6; daa[13*20+ 11] = 10.7;  daa[13*20+ 12] = 321.6; daa[14*20+ 0] = 49.3;
    daa[14*20+ 1] = 0.2;    daa[14*20+ 2] = 16.5;   daa[14*20+ 3] = 0.2;    daa[14*20+ 4] = 0.2;
    daa[14*20+ 5] = 39.3;   daa[14*20+ 6] = 7.9;    daa[14*20+ 7] = 0.2;    daa[14*20+ 8] = 0.8;
    daa[14*20+ 9] = 0.2;    daa[14*20+ 10] = 12.2;  daa[14*20+ 11] = 16.8;  daa[14*20+ 12] = 5.3;
    daa[14*20+ 13] = 14.6;  daa[15*20+ 0] = 673.0;  daa[15*20+ 1] = 2.7;    daa[15*20+ 2] = 398.4;
    daa[15*20+ 3] = 44.4;   daa[15*20+ 4] = 664.2;  daa[15*20+ 5] = 52.4;   daa[15*20+ 6] = 31.5;
    daa[15*20+ 7] = 226.0;  daa[15*20+ 8] = 10.6;   daa[15*20+ 9] = 7.2;    daa[15*20+ 10] = 8.2;
    daa[15*20+ 11] = 144.2; daa[15*20+ 12] = 111.7; daa[15*20+ 13] = 36.1;  daa[15*20+ 14] = 86.5;
    daa[16*20+ 0] = 243.9;  daa[16*20+ 1] = 0.2;    daa[16*20+ 2] = 165.9;  daa[16*20+ 3] = 0.2;
    daa[16*20+ 4] = 182.8;  daa[16*20+ 5] = 43.7;   daa[16*20+ 6] = 43.4;   daa[16*20+ 7] = 0.2;
    daa[16*20+ 8] = 18.6;   daa[16*20+ 9] = 203.7;  daa[16*20+ 10] = 47.8;  daa[16*20+ 11] = 69.5;
    daa[16*20+ 12] = 288.6; daa[16*20+ 13] = 13.5;  daa[16*20+ 14] = 46.8;  daa[16*20+ 15] = 660.4;
    daa[17*20+ 0] = 0.2;    daa[17*20+ 1] = 0.2;    daa[17*20+ 2] = 7.7;    daa[17*20+ 3] = 0.2;
    daa[17*20+ 4] = 21.6;   daa[17*20+ 5] = 6.7;    daa[17*20+ 6] = 11.0;   daa[17*20+ 7] = 1.9;
    daa[17*20+ 8] = 0.2;    daa[17*20+ 9] = 0.2;    daa[17*20+ 10] = 21.1;  daa[17*20+ 11] = 16.0;
    daa[17*20+ 12] = 70.7;  daa[17*20+ 13] = 53.7;  daa[17*20+ 14] = 0.2;   daa[17*20+ 15] = 2.4;
    daa[17*20+ 16] = 0.2;   daa[18*20+ 0] = 1.2;    daa[18*20+ 1] = 3.9;    daa[18*20+ 2] = 251.2;
    daa[18*20+ 3] = 0.2;    daa[18*20+ 4] = 72.0;   daa[18*20+ 5] = 86.7;   daa[18*20+ 6] = 7.7;
    daa[18*20+ 7] = 8.6;    daa[18*20+ 8] = 191.4;  daa[18*20+ 9] = 12.3;   daa[18*20+ 10] = 19.8;
    daa[18*20+ 11] = 117.1; daa[18*20+ 12] = 70.9;  daa[18*20+ 13] = 791.6; daa[18*20+ 14] = 18.4;
    daa[18*20+ 15] = 30.5;  daa[18*20+ 16] = 46.0;  daa[18*20+ 17] = 37.7;  daa[19*20+ 0] = 339.9;
    daa[19*20+ 1] = 0.2;    daa[19*20+ 2] = 22.6;   daa[19*20+ 3] = 0.2;    daa[19*20+ 4] = 350.4;
    daa[19*20+ 5] = 0.2;    daa[19*20+ 6] = 13.6;   daa[19*20+ 7] = 2.6;    daa[19*20+ 8] = 0.2;
    daa[19*20+ 9] = 1854.5; daa[19*20+ 10] = 84.7;  daa[19*20+ 11] = 26.1;  daa[19*20+ 12] = 281.3;
    daa[19*20+ 13] = 51.9;  daa[19*20+ 14] = 31.7;  daa[19*20+ 15] = 60.6;  daa[19*20+ 16] = 544.1;
    daa[19*20+ 17] = 0.2;   daa[19*20+ 18] = 1.6;

/*  MtArt.old: esta es la MtArt que hice con 26 secuencias (2 outgroups) con una topoLOGia incorrecta
    daa[1*20+ 0] = 0.2;     daa[2*20+ 0] = 0.2;     daa[2*20+ 1] = 8.0;     daa[3*20+ 0] = 0.2;
    daa[3*20+ 1] = 0.2;     daa[3*20+ 2] = 441.7;   daa[4*20+ 0] = 287.9;   daa[4*20+ 1] = 48.4;
    daa[4*20+ 2] = 82.4;    daa[4*20+ 3] = 0.2;     daa[5*20+ 0] = 0.2;     daa[5*20+ 1] = 149.9;
    daa[5*20+ 2] = 278.6;   daa[5*20+ 3] = 0.2;     daa[5*20+ 4] = 21.7;    daa[6*20+ 0] = 6.6;
    daa[6*20+ 1] = 0.2;     daa[6*20+ 2] = 213.9;   daa[6*20+ 3] = 760.8;   daa[6*20+ 4] = 0.2;
    daa[6*20+ 5] = 292.9;   daa[7*20+ 0] = 228.2;   daa[7*20+ 1] = 0.2;     daa[7*20+ 2] = 97.1;
    daa[7*20+ 3] = 10.4;    daa[7*20+ 4] = 98.4;    daa[7*20+ 5] = 4.0;     daa[7*20+ 6] = 48.7;
    daa[8*20+ 0] = 0.2;     daa[8*20+ 1] = 56.7;    daa[8*20+ 2] = 156.4;   daa[8*20+ 3] = 24.5;
    daa[8*20+ 4] = 15.5;    daa[8*20+ 5] = 328.6;   daa[8*20+ 6] = 7.0;     daa[8*20+ 7] = 8.4;
    daa[9*20+ 0] = 26.4;    daa[9*20+ 1] = 1.6;     daa[9*20+ 2] = 40.1;    daa[9*20+ 3] = 0.2;
    daa[9*20+ 4] = 22.1;    daa[9*20+ 5] = 13.8;    daa[9*20+ 6] = 0.2;     daa[9*20+ 7] = 3.6;
    daa[9*20+ 8] = 0.2;     daa[10*20+ 0] = 3.4;    daa[10*20+ 1] = 0.6;    daa[10*20+ 2] = 13.8;
    daa[10*20+ 3] = 0.7;    daa[10*20+ 4] = 76.9;   daa[10*20+ 5] = 12.1;   daa[10*20+ 6] = 5.4;
    daa[10*20+ 7] = 2.5;    daa[10*20+ 8] = 2.9;    daa[10*20+ 9] = 542.6;  daa[11*20+ 0] = 0.2;
    daa[11*20+ 1] = 240.2;  daa[11*20+ 2] = 602.8;  daa[11*20+ 3] = 35.5;   daa[11*20+ 4] = 0.2;
    daa[11*20+ 5] = 357.6;  daa[11*20+ 6] = 62.6;   daa[11*20+ 7] = 0.2;    daa[11*20+ 8] = 3.3;
    daa[11*20+ 9] = 0.2;    daa[11*20+ 10] = 17.5;  daa[12*20+ 0] = 119.0;  daa[12*20+ 1] = 0.2;
    daa[12*20+ 2] = 91.4;   daa[12*20+ 3] = 6.4;    daa[12*20+ 4] = 332.3;  daa[12*20+ 5] = 65.4;
    daa[12*20+ 6] = 0.2;    daa[12*20+ 7] = 60.4;   daa[12*20+ 8] = 2.4;    daa[12*20+ 9] = 492.5;
    daa[12*20+ 10] = 815.8; daa[12*20+ 11] = 67.3;  daa[13*20+ 0] = 8.2;    daa[13*20+ 1] = 6.4;
    daa[13*20+ 2] = 31.5;   daa[13*20+ 3] = 3.4;    daa[13*20+ 4] = 174.4;  daa[13*20+ 5] = 5.7;
    daa[13*20+ 6] = 5.7;    daa[13*20+ 7] = 2.1;    daa[13*20+ 8] = 11.0;   daa[13*20+ 9] = 94.4;
    daa[13*20+ 10] = 243.3; daa[13*20+ 11] = 12.3;  daa[13*20+ 12] = 357.8; daa[14*20+ 0] = 62.5;
    daa[14*20+ 1] = 0.4;    daa[14*20+ 2] = 17.5;   daa[14*20+ 3] = 0.2;    daa[14*20+ 4] = 0.2;
    daa[14*20+ 5] = 48.6;   daa[14*20+ 6] = 17.7;   daa[14*20+ 7] = 2.7;    daa[14*20+ 8] = 0.2;
    daa[14*20+ 9] = 0.2;    daa[14*20+ 10] = 11.2;  daa[14*20+ 11] = 21.7;  daa[14*20+ 12] = 5.2;
    daa[14*20+ 13] = 12.6;  daa[15*20+ 0] = 659.0;  daa[15*20+ 1] = 5.2;    daa[15*20+ 2] = 469.8;
    daa[15*20+ 3] = 52.3;   daa[15*20+ 4] = 570.7;  daa[15*20+ 5] = 47.8;   daa[15*20+ 6] = 37.3;
    daa[15*20+ 7] = 227.8;  daa[15*20+ 8] = 12.7;   daa[15*20+ 9] = 12.3;   daa[15*20+ 10] = 7.4;
    daa[15*20+ 11] = 189.0; daa[15*20+ 12] = 155.3; daa[15*20+ 13] = 43.8;  daa[15*20+ 14] = 103.4;
    daa[16*20+ 0] = 276.4;  daa[16*20+ 1] = 1.6;    daa[16*20+ 2] = 175.6;  daa[16*20+ 3] = 0.2;
    daa[16*20+ 4] = 96.2;   daa[16*20+ 5] = 71.4;   daa[16*20+ 6] = 37.4;   daa[16*20+ 7] = 0.2;
    daa[16*20+ 8] = 14.2;   daa[16*20+ 9] = 212.5;  daa[16*20+ 10] = 38.5;  daa[16*20+ 11] = 97.4;
    daa[16*20+ 12] = 254.7; daa[16*20+ 13] = 2.1;   daa[16*20+ 14] = 41.6;  daa[16*20+ 15] = 670.6;
    daa[17*20+ 0] = 6.2;    daa[17*20+ 1] = 0.2;    daa[17*20+ 2] = 0.2;    daa[17*20+ 3] = 5.6;
    daa[17*20+ 4] = 0.2;    daa[17*20+ 5] = 0.2;    daa[17*20+ 6] = 3.1;    daa[17*20+ 7] = 0.4;
    daa[17*20+ 8] = 0.2;    daa[17*20+ 9] = 15.2;   daa[17*20+ 10] = 11.5;  daa[17*20+ 11] = 32.6;
    daa[17*20+ 12] = 82.4;  daa[17*20+ 13] = 81.9;  daa[17*20+ 14] = 0.2;   daa[17*20+ 15] = 9.7;
    daa[17*20+ 16] = 0.2;   daa[18*20+ 0] = 1.6;    daa[18*20+ 1] = 7.7;    daa[18*20+ 2] = 242.5;
    daa[18*20+ 3] = 0.2;    daa[18*20+ 4] = 88.0;   daa[18*20+ 5] = 93.1;   daa[18*20+ 6] = 0.2;
    daa[18*20+ 7] = 6.0;    daa[18*20+ 8] = 113.7;  daa[18*20+ 9] = 22.1;   daa[18*20+ 10] = 17.2;
    daa[18*20+ 11] = 138.5; daa[18*20+ 12] = 37.6;  daa[18*20+ 13] = 770.2; daa[18*20+ 14] = 5.3;
    daa[18*20+ 15] = 25.0;  daa[18*20+ 16] = 55.5;  daa[18*20+ 17] = 69.3;  daa[19*20+ 0] = 307.8;
    daa[19*20+ 1] = 2.2;    daa[19*20+ 2] = 6.9;    daa[19*20+ 3] = 0.2;    daa[19*20+ 4] = 405.7;
    daa[19*20+ 5] = 0.8;    daa[19*20+ 6] = 20.2;   daa[19*20+ 7] = 5.7;    daa[19*20+ 8] = 0.2;
    daa[19*20+ 9] = 1687.9; daa[19*20+ 10] = 49.4;  daa[19*20+ 11] = 23.4;  daa[19*20+ 12] = 329.9;
    daa[19*20+ 13] = 86.3;  daa[19*20+ 14] = 27.3;  daa[19*20+ 15] = 95.0;  daa[19*20+ 16] = 443.0;
    daa[19*20+ 17] = 2.4;   daa[19*20+ 18] = 0.2;
3*/
    for (i=0; i<naa; i++) for (j=0; j<i; j++) daa[j*naa+i] = daa[i*naa+j];

    pi[0] = 0.054116;       pi[1] = 0.018227;       pi[2] = 0.039903;       pi[3] = 0.020160;       pi[4] = 0.009709;
    pi[5] = 0.018781;       pi[6] = 0.024289;       pi[7] = 0.068183;       pi[8] = 0.024518;       pi[9] = 0.092639;
    pi[10] = 0.148658;      pi[11] = 0.021718;      pi[12] = 0.061453;      pi[13] = 0.088668;      pi[14] = 0.041826;
    pi[15] = 0.091030;      pi[16] = 0.049194;      pi[17] = 0.029786;      pi[18] = 0.039443;      pi[19] = 0.057701;
/*  pi[0] = 0.055505;       pi[1] = 0.018320;       pi[2] = 0.036614;       pi[3] = 0.019517;    pi[4] = 0.009878;
    pi[5] = 0.018732;       pi[6] = 0.021591;       pi[7] = 0.068414;       pi[8] = 0.023638;    pi[9] = 0.092558;
    pi[10] = 0.154697;      pi[11] = 0.020168;      pi[12] = 0.062033;      pi[13] = 0.089433;   pi[14] = 0.040642;
    pi[15] = 0.090298;      pi[16] = 0.050134;      pi[17] = 0.027693;      pi[18] = 0.038448;   pi[19] = 0.061687;
*/
    return 1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


int Init_Qmat_HIVb(phydbl *daa, phydbl *pi) //added by FEDE
{
    /*
      Nickle DC, Heath L, Jensen MA, Gilbert PB, Mullins JI, Kosakovsky Pond SL.
      HIV-Specific Probabilistic Models of Protein Evolution.
      PLoS ONE. 2007 Jun 6;2:e503.

      [thanks to Sergei L. Kosakovsky]

      Translated from HYPHY to Phyml format by Federico Abascal.
    */

    int i,j,naa;
    naa = 20;

    daa[1*20+0]= 0.307507;        daa[2*20+0]= 0.005;           daa[2*20+1]= 0.295543;        daa[3*20+0]= 1.45504;
    daa[3*20+1]= 0.005;           daa[3*20+2]= 17.6612;         daa[4*20+0]= 0.123758;        daa[4*20+1]= 0.351721;
    daa[4*20+2]= 0.0860642;       daa[4*20+3]= 0.005;           daa[5*20+0]= 0.0551128;       daa[5*20+1]= 3.4215;
    daa[5*20+2]= 0.672052;        daa[5*20+3]= 0.005;           daa[5*20+4]= 0.005;           daa[6*20+0]= 1.48135;
    daa[6*20+1]= 0.0749218;       daa[6*20+2]= 0.0792633;       daa[6*20+3]= 10.5872;         daa[6*20+4]= 0.005;
    daa[6*20+5]= 2.5602;          daa[7*20+0]= 2.13536;         daa[7*20+1]= 3.65345;         daa[7*20+2]= 0.323401;
    daa[7*20+3]= 2.83806;         daa[7*20+4]= 0.897871;        daa[7*20+5]= 0.0619137;       daa[7*20+6]= 3.92775;
    daa[8*20+0]= 0.0847613;       daa[8*20+1]= 9.04044;         daa[8*20+2]= 7.64585;         daa[8*20+3]= 1.9169;
    daa[8*20+4]= 0.240073;        daa[8*20+5]= 7.05545;         daa[8*20+6]= 0.11974;         daa[8*20+7]= 0.005;
    daa[9*20+0]= 0.005;           daa[9*20+1]= 0.677289;        daa[9*20+2]= 0.680565;        daa[9*20+3]= 0.0176792;
    daa[9*20+4]= 0.005;           daa[9*20+5]= 0.005;           daa[9*20+6]= 0.00609079;      daa[9*20+7]= 0.005;
    daa[9*20+8]= 0.103111;        daa[10*20+0]= 0.215256;       daa[10*20+1]= 0.701427;       daa[10*20+2]= 0.005;
    daa[10*20+3]= 0.00876048;     daa[10*20+4]= 0.129777;       daa[10*20+5]= 1.49456;        daa[10*20+6]= 0.005;
    daa[10*20+7]= 0.005;          daa[10*20+8]= 1.74171;        daa[10*20+9]= 5.95879;        daa[11*20+0]= 0.005;
    daa[11*20+1]= 20.45;          daa[11*20+2]= 7.90443;        daa[11*20+3]= 0.005;          daa[11*20+4]= 0.005;
    daa[11*20+5]= 6.54737;        daa[11*20+6]= 4.61482;        daa[11*20+7]= 0.521705;       daa[11*20+8]= 0.005;
    daa[11*20+9]= 0.322319;       daa[11*20+10]= 0.0814995;     daa[12*20+0]= 0.0186643;      daa[12*20+1]= 2.51394;
    daa[12*20+2]= 0.005;          daa[12*20+3]= 0.005;          daa[12*20+4]= 0.005;          daa[12*20+5]= 0.303676;
    daa[12*20+6]= 0.175789;       daa[12*20+7]= 0.005;          daa[12*20+8]= 0.005;          daa[12*20+9]= 11.2065;
    daa[12*20+10]= 5.31961;       daa[12*20+11]= 1.28246;       daa[13*20+0]= 0.0141269;      daa[13*20+1]= 0.005;
    daa[13*20+2]= 0.005;          daa[13*20+3]= 0.005;          daa[13*20+4]= 9.29815;        daa[13*20+5]= 0.005;
    daa[13*20+6]= 0.005;          daa[13*20+7]= 0.291561;       daa[13*20+8]= 0.145558;       daa[13*20+9]= 3.39836;
    daa[13*20+10]= 8.52484;       daa[13*20+11]= 0.0342658;     daa[13*20+12]= 0.188025;      daa[14*20+0]= 2.12217;
    daa[14*20+1]= 1.28355;        daa[14*20+2]= 0.00739578;     daa[14*20+3]= 0.0342658;      daa[14*20+4]= 0.005;
    daa[14*20+5]= 4.47211;        daa[14*20+6]= 0.0120226;      daa[14*20+7]= 0.005;          daa[14*20+8]= 2.45318;
    daa[14*20+9]= 0.0410593;      daa[14*20+10]= 2.07757;       daa[14*20+11]= 0.0313862;     daa[14*20+12]= 0.005;
    daa[14*20+13]= 0.005;         daa[15*20+0]= 2.46633;        daa[15*20+1]= 3.4791;         daa[15*20+2]= 13.1447;
    daa[15*20+3]= 0.52823;        daa[15*20+4]= 4.69314;        daa[15*20+5]= 0.116311;       daa[15*20+6]= 0.005;
    daa[15*20+7]= 4.38041;        daa[15*20+8]= 0.382747;       daa[15*20+9]= 1.21803;        daa[15*20+10]= 0.927656;
    daa[15*20+11]= 0.504111;      daa[15*20+12]= 0.005;         daa[15*20+13]= 0.956472;      daa[15*20+14]= 5.37762;
    daa[16*20+0]= 15.9183;        daa[16*20+1]= 2.86868;        daa[16*20+2]= 6.88667;        daa[16*20+3]= 0.274724;
    daa[16*20+4]= 0.739969;       daa[16*20+5]= 0.243589;       daa[16*20+6]= 0.289774;       daa[16*20+7]= 0.369615;
    daa[16*20+8]= 0.711594;       daa[16*20+9]= 8.61217;        daa[16*20+10]= 0.0437673;     daa[16*20+11]= 4.67142;
    daa[16*20+12]= 4.94026;       daa[16*20+13]= 0.0141269;     daa[16*20+14]= 2.01417;       daa[16*20+15]= 8.93107;
    daa[17*20+0]= 0.005;          daa[17*20+1]= 0.991338;       daa[17*20+2]= 0.005;          daa[17*20+3]= 0.005;
    daa[17*20+4]= 2.63277;        daa[17*20+5]= 0.026656;       daa[17*20+6]= 0.005;          daa[17*20+7]= 1.21674;
    daa[17*20+8]= 0.0695179;      daa[17*20+9]= 0.005;          daa[17*20+10]= 0.748843;      daa[17*20+11]= 0.005;
    daa[17*20+12]= 0.089078;      daa[17*20+13]= 0.829343;      daa[17*20+14]= 0.0444506;     daa[17*20+15]= 0.0248728;
    daa[17*20+16]= 0.005;         daa[18*20+0]= 0.005;          daa[18*20+1]= 0.00991826;     daa[18*20+2]= 1.76417;
    daa[18*20+3]= 0.674653;       daa[18*20+4]= 7.57932;        daa[18*20+5]= 0.113033;       daa[18*20+6]= 0.0792633;
    daa[18*20+7]= 0.005;          daa[18*20+8]= 18.6943;        daa[18*20+9]= 0.148168;       daa[18*20+10]= 0.111986;
    daa[18*20+11]= 0.005;         daa[18*20+12]= 0.005;         daa[18*20+13]= 15.34;         daa[18*20+14]= 0.0304381;
    daa[18*20+15]= 0.648024;      daa[18*20+16]= 0.105652;      daa[18*20+17]= 1.28022;       daa[19*20+0]= 7.61428;
    daa[19*20+1]= 0.0812454;      daa[19*20+2]= 0.026656;       daa[19*20+3]= 1.04793;        daa[19*20+4]= 0.420027;
    daa[19*20+5]= 0.0209153;      daa[19*20+6]= 1.02847;        daa[19*20+7]= 0.953155;       daa[19*20+8]= 0.005;
    daa[19*20+9]= 17.7389;        daa[19*20+10]= 1.41036;       daa[19*20+11]= 0.265829;      daa[19*20+12]= 6.8532;
    daa[19*20+13]= 0.723274;      daa[19*20+14]= 0.005;         daa[19*20+15]= 0.0749218;     daa[19*20+16]= 0.709226;
    daa[19*20+17]= 0.005;         daa[19*20+18]= 0.0410593;
    for (i=0; i<naa; i++)  for (j=0; j<i; j++)  daa[j*naa+i] = daa[i*naa+j];

    pi[0]= 0.060490222;           pi[1]= 0.066039665;           pi[2]= 0.044127815;           pi[3]= 0.042109048;
    pi[4]= 0.020075899;           pi[5]= 0.053606488;           pi[6]= 0.071567447;           pi[7]= 0.072308239;
    pi[8]= 0.022293943;           pi[9]= 0.069730629;           pi[10]= 0.098851122;          pi[11]= 0.056968211;
    pi[12]= 0.019768318;          pi[13]= 0.028809447;          pi[14]= 0.046025282;          pi[15]= 0.05060433;
    pi[16]= 0.053636813;          pi[17]= 0.033011601;          pi[18]= 0.028350243;          pi[19]= 0.061625237;
    return 1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int Init_Qmat_HIVw(phydbl *daa, phydbl *pi)
{
    /*
      Nickle DC, Heath L, Jensen MA, Gilbert PB, Mullins JI, Kosakovsky Pond SL.
      HIV-Specific Probabilistic Models of Protein Evolution.
      PLoS ONE. 2007 Jun 6;2:e503.

      [thanks to Sergei L. Kosakovsky]

      Translated from HYPHY to Phyml format by Federico Abascal.
    */

    int i,j,naa;
    naa = 20;

    daa[1*20+0]= 0.0744808;       daa[2*20+0]= 0.617509;        daa[2*20+1]= 0.16024;         daa[3*20+0]= 4.43521;
    daa[3*20+1]= 0.0674539;       daa[3*20+2]= 29.4087;         daa[4*20+0]= 0.167653;        daa[4*20+1]= 2.86364;
    daa[4*20+2]= 0.0604932;       daa[4*20+3]= 0.005;           daa[5*20+0]= 0.005;           daa[5*20+1]= 10.6746;
    daa[5*20+2]= 0.342068;        daa[5*20+3]= 0.005;           daa[5*20+4]= 0.005;           daa[6*20+0]= 5.56325;
    daa[6*20+1]= 0.0251632;       daa[6*20+2]= 0.201526;        daa[6*20+3]= 12.1233;         daa[6*20+4]= 0.005;
    daa[6*20+5]= 3.20656;         daa[7*20+0]= 1.8685;          daa[7*20+1]= 13.4379;         daa[7*20+2]= 0.0604932;
    daa[7*20+3]= 10.3969;         daa[7*20+4]= 0.0489798;       daa[7*20+5]= 0.0604932;       daa[7*20+6]= 14.7801;
    daa[8*20+0]= 0.005;           daa[8*20+1]= 6.84405;         daa[8*20+2]= 8.59876;         daa[8*20+3]= 2.31779;
    daa[8*20+4]= 0.005;           daa[8*20+5]= 18.5465;         daa[8*20+6]= 0.005;           daa[8*20+7]= 0.005;
    daa[9*20+0]= 0.005;           daa[9*20+1]= 1.34069;         daa[9*20+2]= 0.987028;        daa[9*20+3]= 0.145124;
    daa[9*20+4]= 0.005;           daa[9*20+5]= 0.0342252;       daa[9*20+6]= 0.0390512;       daa[9*20+7]= 0.005;
    daa[9*20+8]= 0.005;           daa[10*20+0]= 0.16024;        daa[10*20+1]= 0.586757;       daa[10*20+2]= 0.005;
    daa[10*20+3]= 0.005;          daa[10*20+4]= 0.005;          daa[10*20+5]= 2.89048;        daa[10*20+6]= 0.129839;
    daa[10*20+7]= 0.0489798;      daa[10*20+8]= 1.76382;        daa[10*20+9]= 9.10246;        daa[11*20+0]= 0.592784;
    daa[11*20+1]= 39.8897;        daa[11*20+2]= 10.6655;        daa[11*20+3]= 0.894313;       daa[11*20+4]= 0.005;
    daa[11*20+5]= 13.0705;        daa[11*20+6]= 23.9626;        daa[11*20+7]= 0.279425;       daa[11*20+8]= 0.22406;
    daa[11*20+9]= 0.817481;       daa[11*20+10]= 0.005;         daa[12*20+0]= 0.005;          daa[12*20+1]= 3.28652;
    daa[12*20+2]= 0.201526;       daa[12*20+3]= 0.005;          daa[12*20+4]= 0.005;          daa[12*20+5]= 0.005;
    daa[12*20+6]= 0.005;          daa[12*20+7]= 0.0489798;      daa[12*20+8]= 0.005;          daa[12*20+9]= 17.3064;
    daa[12*20+10]= 11.3839;       daa[12*20+11]= 4.09564;       daa[13*20+0]= 0.597923;       daa[13*20+1]= 0.005;
    daa[13*20+2]= 0.005;          daa[13*20+3]= 0.005;          daa[13*20+4]= 0.362959;       daa[13*20+5]= 0.005;
    daa[13*20+6]= 0.005;          daa[13*20+7]= 0.005;          daa[13*20+8]= 0.005;          daa[13*20+9]= 1.48288;
    daa[13*20+10]= 7.48781;       daa[13*20+11]= 0.005;         daa[13*20+12]= 0.005;         daa[14*20+0]= 1.00981;
    daa[14*20+1]= 0.404723;       daa[14*20+2]= 0.344848;       daa[14*20+3]= 0.005;          daa[14*20+4]= 0.005;
    daa[14*20+5]= 3.04502;        daa[14*20+6]= 0.005;          daa[14*20+7]= 0.005;          daa[14*20+8]= 13.9444;
    daa[14*20+9]= 0.005;          daa[14*20+10]= 9.83095;       daa[14*20+11]= 0.111928;      daa[14*20+12]= 0.005;
    daa[14*20+13]= 0.0342252;     daa[15*20+0]= 8.5942;         daa[15*20+1]= 8.35024;        daa[15*20+2]= 14.5699;
    daa[15*20+3]= 0.427881;       daa[15*20+4]= 1.12195;        daa[15*20+5]= 0.16024;        daa[15*20+6]= 0.005;
    daa[15*20+7]= 6.27966;        daa[15*20+8]= 0.725157;       daa[15*20+9]= 0.740091;       daa[15*20+10]= 6.14396;
    daa[15*20+11]= 0.005;         daa[15*20+12]= 0.392575;      daa[15*20+13]= 4.27939;       daa[15*20+14]= 14.249;
    daa[16*20+0]= 24.1422;        daa[16*20+1]= 0.928203;       daa[16*20+2]= 4.54206;        daa[16*20+3]= 0.630395;
    daa[16*20+4]= 0.005;          daa[16*20+5]= 0.203091;       daa[16*20+6]= 0.458743;       daa[16*20+7]= 0.0489798;
    daa[16*20+8]= 0.95956;        daa[16*20+9]= 9.36345;        daa[16*20+10]= 0.005;         daa[16*20+11]= 4.04802;
    daa[16*20+12]= 7.41313;       daa[16*20+13]= 0.114512;      daa[16*20+14]= 4.33701;       daa[16*20+15]= 6.34079;
    daa[17*20+0]= 0.005;          daa[17*20+1]= 5.96564;        daa[17*20+2]= 0.005;          daa[17*20+3]= 0.005;
    daa[17*20+4]= 5.49894;        daa[17*20+5]= 0.0443298;      daa[17*20+6]= 0.005;          daa[17*20+7]= 2.8258;
    daa[17*20+8]= 0.005;          daa[17*20+9]= 0.005;          daa[17*20+10]= 1.37031;       daa[17*20+11]= 0.005;
    daa[17*20+12]= 0.005;         daa[17*20+13]= 0.005;         daa[17*20+14]= 0.005;         daa[17*20+15]= 1.10156;
    daa[17*20+16]= 0.005;         daa[18*20+0]= 0.005;          daa[18*20+1]= 0.005;          daa[18*20+2]= 5.06475;
    daa[18*20+3]= 2.28154;        daa[18*20+4]= 8.34835;        daa[18*20+5]= 0.005;          daa[18*20+6]= 0.005;
    daa[18*20+7]= 0.005;          daa[18*20+8]= 47.4889;        daa[18*20+9]= 0.114512;       daa[18*20+10]= 0.005;
    daa[18*20+11]= 0.005;         daa[18*20+12]= 0.579198;      daa[18*20+13]= 4.12728;       daa[18*20+14]= 0.005;
    daa[18*20+15]= 0.933142;      daa[18*20+16]= 0.490608;      daa[18*20+17]= 0.005;         daa[19*20+0]= 24.8094;
    daa[19*20+1]= 0.279425;       daa[19*20+2]= 0.0744808;      daa[19*20+3]= 2.91786;        daa[19*20+4]= 0.005;
    daa[19*20+5]= 0.005;          daa[19*20+6]= 2.19952;        daa[19*20+7]= 2.79622;        daa[19*20+8]= 0.827479;
    daa[19*20+9]= 24.8231;        daa[19*20+10]= 2.95344;       daa[19*20+11]= 0.128065;      daa[19*20+12]= 14.7683;
    daa[19*20+13]= 2.28;          daa[19*20+14]= 0.005;         daa[19*20+15]= 0.862637;      daa[19*20+16]= 0.005;
    daa[19*20+17]= 0.005;         daa[19*20+18]= 1.35482;
    for (i=0; i<naa; i++)  for (j=0; j<i; j++)  daa[j*naa+i] = daa[i*naa+j];

    pi[0]= 0.0377494;             pi[1]= 0.057321;              pi[2]= 0.0891129;             pi[3]= 0.0342034;
    pi[4]= 0.0240105;             pi[5]= 0.0437824;             pi[6]= 0.0618606;             pi[7]= 0.0838496;
    pi[8]= 0.0156076;             pi[9]= 0.0983641;             pi[10]= 0.0577867;            pi[11]= 0.0641682;
    pi[12]= 0.0158419;            pi[13]= 0.0422741;            pi[14]= 0.0458601;            pi[15]= 0.0550846;
    pi[16]= 0.0813774;            pi[17]= 0.019597;             pi[18]= 0.0205847;            pi[19]= 0.0515639;
    return 1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int Init_Qmat_JTT(phydbl *daa, phydbl *pi)
{
  int i,j,naa;

  /* JTT's model data
   * D.T.Jones, W.R.Taylor and J.M.Thornton
   * "The rapid generation of mutation data matrices from protein sequences"
   * CABIOS  vol.8 no.3 1992 pp275-282
   */

  naa = 20;

/*   PhyML_Printf("\n\n. REMINDER : THIS IS NOT JTT !!!\n\n"); */


/*   daa[1*20 + 0] = 0.592439084;  */
/*   daa[2*20 + 0] = 0.427686326;  */
/*   daa[2*20 + 1] = 0.795121783;  */
/*   daa[3*20 + 0] = 0.649119246;  */
/*   daa[3*20 + 1] = 0.116285502;  */
/*   daa[3*20 + 2] = 6.004659556;  */
/*   daa[4*20 + 0] = 2.384113318;  */
/*   daa[4*20 + 1] = 0.925080814;  */
/*   daa[4*20 + 2] = 0.648459384;  */
/*   daa[4*20 + 3] = 0.088632079;  */
/*   daa[5*20 + 0] = 1.188499461;  */
/*   daa[5*20 + 1] = 3.230029905;  */
/*   daa[5*20 + 2] = 1.910660418;  */
/*   daa[5*20 + 3] = 0.599639374;  */
/*   daa[5*20 + 4] = 0.123291416;  */
/*   daa[6*20 + 0] = 1.534092487;  */
/*   daa[6*20 + 1] = 0.390654981;  */
/*   daa[6*20 + 2] = 0.566472743;  */
/*   daa[6*20 + 3] = 5.916482256;  */
/*   daa[6*20 + 4] = 0.001457992;  */
/*   daa[6*20 + 5] = 4.564894849;  */
/*   daa[7*20 + 0] = 2.282012114;  */
/*   daa[7*20 + 1] = 0.502549742;  */
/*   daa[7*20 + 2] = 1.97013036;  */
/*   daa[7*20 + 3] = 1.198016556;  */
/*   daa[7*20 + 4] = 0.565245125;  */
/*   daa[7*20 + 5] = 0.382529852;  */
/*   daa[7*20 + 6] = 0.477515495;  */
/*   daa[8*20 + 0] = 0.447587408;  */
/*   daa[8*20 + 1] = 2.948365572;  */
/*   daa[8*20 + 2] = 5.900694893;  */
/*   daa[8*20 + 3] = 0.96157692;  */
/*   daa[8*20 + 4] = 0.757153271;  */
/*   daa[8*20 + 5] = 5.048587596;  */
/*   daa[8*20 + 6] = 0.410583482;  */
/*   daa[8*20 + 7] = 0.406115605;  */
/*   daa[9*20 + 0] = 0.105337934;  */
/*   daa[9*20 + 1] = 0.115421183;  */
/*   daa[9*20 + 2] = 0.178948502;  */
/*   daa[9*20 + 3] = 0.01319564;  */
/*   daa[9*20 + 4] = 0.269637325;  */
/*   daa[9*20 + 5] = 0.087618072;  */
/*   daa[9*20 + 6] = 0.045232757;  */
/*   daa[9*20 + 7] = 0.00958642;  */
/*   daa[9*20 + 8] = 0.116145048;  */
/*   daa[10*20 + 0] = 0.269460438;  */
/*   daa[10*20 + 1] = 0.416549984;  */
/*   daa[10*20 + 2] = 0.100490196;  */
/*   daa[10*20 + 3] = 0.021637571;  */
/*   daa[10*20 + 4] = 0.610442865;  */
/*   daa[10*20 + 5] = 0.782347175;  */
/*   daa[10*20 + 6] = 0.077565847;  */
/*   daa[10*20 + 7] = 0.034916545;  */
/*   daa[10*20 + 8] = 0.456334781;  */
/*   daa[10*20 + 9] = 4.057302736;  */
/*   daa[11*20 + 0] = 0.71706921;  */
/*   daa[11*20 + 1] = 6.469532749;  */
/*   daa[11*20 + 2] = 2.606469511;  */
/*   daa[11*20 + 3] = 0.239488567;  */
/*   daa[11*20 + 4] = 0.020594132;  */
/*   daa[11*20 + 5] = 3.978479772;  */
/*   daa[11*20 + 6] = 1.669092348;  */
/*   daa[11*20 + 7] = 0.368262509;  */
/*   daa[11*20 + 8] = 0.744835535;  */
/*   daa[11*20 + 9] = 0.153626461;  */
/*   daa[11*20 + 10] = 0.186106976;  */
/*   daa[12*20 + 0] = 0.766271778;  */
/*   daa[12*20 + 1] = 0.465902609;  */
/*   daa[12*20 + 2] = 0.352429433;  */
/*   daa[12*20 + 3] = 0.04670109;  */
/*   daa[12*20 + 4] = 0.930016402;  */
/*   daa[12*20 + 5] = 1.529072358;  */
/*   daa[12*20 + 6] = 0.149497081;  */
/*   daa[12*20 + 7] = 0.092795029;  */
/*   daa[12*20 + 8] = 0.367452788;  */
/*   daa[12*20 + 9] = 4.443975824;  */
/*   daa[12*20 + 10] = 7.305963709;  */
/*   daa[12*20 + 11] = 0.691067831;  */
/*   daa[13*20 + 0] = 0.164314785;  */
/*   daa[13*20 + 1] = 0.064046051;  */
/*   daa[13*20 + 2] = 0.108166343;  */
/*   daa[13*20 + 3] = 0.021096631;  */
/*   daa[13*20 + 4] = 0.9869692;  */
/*   daa[13*20 + 5] = 0.054959009;  */
/*   daa[13*20 + 6] = 0.025919221;  */
/*   daa[13*20 + 7] = 0.053577857;  */
/*   daa[13*20 + 8] = 0.933637734;  */
/*   daa[13*20 + 9] = 0.928915383;  */
/*   daa[13*20 + 10] = 2.549658015;  */
/*   daa[13*20 + 11] = 0.032602723;  */
/*   daa[13*20 + 12] = 1.837397986;  */
/*   daa[14*20 + 0] = 1.57937466;  */
/*   daa[14*20 + 1] = 0.360749538;  */
/*   daa[14*20 + 2] = 0.189183366;  */
/*   daa[14*20 + 3] = 0.475387263;  */
/*   daa[14*20 + 4] = 0.108255878;  */
/*   daa[14*20 + 5] = 0.645294956;  */
/*   daa[14*20 + 6] = 0.537278362;  */
/*   daa[14*20 + 7] = 0.264516985;  */
/*   daa[14*20 + 8] = 0.510062949;  */
/*   daa[14*20 + 9] = 0.063380048;  */
/*   daa[14*20 + 10] = 0.25549658;  */
/*   daa[14*20 + 11] = 0.42259906;  */
/*   daa[14*20 + 12] = 0.105330979;  */
/*   daa[14*20 + 13] = 0.092650543;  */
/*   daa[15*20 + 0] = 5.358849732;  */
/*   daa[15*20 + 1] = 0.980362668;  */
/*   daa[15*20 + 2] = 5.114953674;  */
/*   daa[15*20 + 3] = 1.242898123;  */
/*   daa[15*20 + 4] = 3.417258976;  */
/*   daa[15*20 + 5] = 1.199117804;  */
/*   daa[15*20 + 6] = 0.61024061;  */
/*   daa[15*20 + 7] = 2.200746452;  */
/*   daa[15*20 + 8] = 0.928717085;  */
/*   daa[15*20 + 9] = 0.07737779;  */
/*   daa[15*20 + 10] = 0.17657022;  */
/*   daa[15*20 + 11] = 0.768621405;  */
/*   daa[15*20 + 12] = 0.330215244;  */
/*   daa[15*20 + 13] = 0.359985962;  */
/*   daa[15*20 + 14] = 1.559718375;  */
/*   daa[16*20 + 0] = 2.21115657;  */
/*   daa[16*20 + 1] = 0.711040813;  */
/*   daa[16*20 + 2] = 2.273163111;  */
/*   daa[16*20 + 3] = 0.408851733;  */
/*   daa[16*20 + 4] = 1.57572444;  */
/*   daa[16*20 + 5] = 1.074674746;  */
/*   daa[16*20 + 6] = 0.707096016;  */
/*   daa[16*20 + 7] = 0.157418108;  */
/*   daa[16*20 + 8] = 0.543088926;  */
/*   daa[16*20 + 9] = 1.056024176;  */
/*   daa[16*20 + 10] = 0.294118282;  */
/*   daa[16*20 + 11] = 1.325676171;  */
/*   daa[16*20 + 12] = 1.834248079;  */
/*   daa[16*20 + 13] = 0.136987103;  */
/*   daa[16*20 + 14] = 0.63778855;  */
/*   daa[16*20 + 15] = 7.294506488;  */
/*   daa[17*20 + 0] = 0.142106331;  */
/*   daa[17*20 + 1] = 0.600995696;  */
/*   daa[17*20 + 2] = 0.079481793;  */
/*   daa[17*20 + 3] = 0.069125482;  */
/*   daa[17*20 + 4] = 0.701931839;  */
/*   daa[17*20 + 5] = 0.23881661;  */
/*   daa[17*20 + 6] = 0.082469159;  */
/*   daa[17*20 + 7] = 0.185722676;  */
/*   daa[17*20 + 8] = 0.911694299;  */
/*   daa[17*20 + 9] = 0.11679729;  */
/*   daa[17*20 + 10] = 0.688555894;  */
/*   daa[17*20 + 11] = 0.063328042;  */
/*   daa[17*20 + 12] = 0.780079225;  */
/*   daa[17*20 + 13] = 3.432193724;  */
/*   daa[17*20 + 14] = 0.088827331;  */
/*   daa[17*20 + 15] = 0.293602516;  */
/*   daa[17*20 + 16] = 0.158083832;  */
/*   daa[18*20 + 0] = 0.178618354;  */
/*   daa[18*20 + 1] = 0.243634289;  */
/*   daa[18*20 + 2] = 0.567657832;  */
/*   daa[18*20 + 3] = 0.134890583;  */
/*   daa[18*20 + 4] = 1.28093676;  */
/*   daa[18*20 + 5] = 0.180560506;  */
/*   daa[18*20 + 6] = 0.081426157;  */
/*   daa[18*20 + 7] = 0.044067219;  */
/*   daa[18*20 + 8] = 6.464181222;  */
/*   daa[18*20 + 9] = 0.206728215;  */
/*   daa[18*20 + 10] = 0.334379787;  */
/*   daa[18*20 + 11] = 0.102051407;  */
/*   daa[18*20 + 12] = 0.431775879;  */
/*   daa[18*20 + 13] = 11.87781125;  */
/*   daa[18*20 + 14] = 0.078070175;  */
/*   daa[18*20 + 15] = 0.359969114;  */
/*   daa[18*20 + 16] = 0.199064371;  */
/*   daa[18*20 + 17] = 4.185797039;  */
/*   daa[19*20 + 0] = 1.964727121;  */
/*   daa[19*20 + 1] = 0.18579405;  */
/*   daa[19*20 + 2] = 0.113068304;  */
/*   daa[19*20 + 3] = 0.046258503;  */
/*   daa[19*20 + 4] = 2.069983598;  */
/*   daa[19*20 + 5] = 0.247259847;  */
/*   daa[19*20 + 6] = 0.272991021;  */
/*   daa[19*20 + 7] = 0.073432734;  */
/*   daa[19*20 + 8] = 0.158879333;  */
/*   daa[19*20 + 9] = 10.52590329;  */
/*   daa[19*20 + 10] = 1.74812825;  */
/*   daa[19*20 + 11] = 0.188152512;  */
/*   daa[19*20 + 12] = 1.800124087;  */
/*   daa[19*20 + 13] = 0.478606732;  */
/*   daa[19*20 + 14] = 0.277216066;  */
/*   daa[19*20 + 15] = 0.104728903;  */
/*   daa[19*20 + 16] = 1.963379491;  */
/*   daa[19*20 + 17] = 0.16552242;  */
/*   daa[19*20 + 18] = 0.222040517; */

/*   for (i=0; i<naa; i++)  for (j=0; j<i; j++)  daa[j*naa+i] = daa[i*naa+j]; */

/*   pi[0] = 0.093862;  */
/*   pi[1] = 0.05757;  */
/*   pi[2] = 0.037017;  */
/*   pi[3] = 0.060952;  */
/*   pi[4] = 0.011004;  */
/*   pi[5] = 0.044631;  */
/*   pi[6] = 0.083648;  */
/*   pi[7] = 0.053004;  */
/*   pi[8] = 0.022528;  */
/*   pi[9] = 0.060152;  */
/*   pi[10] = 0.092314;  */
/*   pi[11] = 0.065886;  */
/*   pi[12] = 0.021075;  */
/*   pi[13] = 0.033928;  */
/*   pi[14] = 0.043333;  */
/*   pi[15] = 0.052957;  */
/*   pi[16] = 0.053465;  */
/*   pi[17] = 0.00928;  */
/*   pi[18] = 0.030025;  */
/*   pi[19] = 0.073369; */

  daa[ 1*20+ 0] =   58.00; daa[ 2*20+ 0] =   54.00; daa[ 2*20+ 1] =   45.00; daa[ 3*20+ 0] =   81.00;
  daa[ 3*20+ 1] =   16.00; daa[ 3*20+ 2] =  528.00; daa[ 4*20+ 0] =   56.00; daa[ 4*20+ 1] =  113.00;
  daa[ 4*20+ 2] =   34.00; daa[ 4*20+ 3] =   10.00; daa[ 5*20+ 0] =   57.00; daa[ 5*20+ 1] =  310.00;
  daa[ 5*20+ 2] =   86.00; daa[ 5*20+ 3] =   49.00; daa[ 5*20+ 4] =    9.00; daa[ 6*20+ 0] =  105.00;
  daa[ 6*20+ 1] =   29.00; daa[ 6*20+ 2] =   58.00; daa[ 6*20+ 3] =  767.00; daa[ 6*20+ 4] =    5.00;
  daa[ 6*20+ 5] =  323.00; daa[ 7*20+ 0] =  179.00; daa[ 7*20+ 1] =  137.00; daa[ 7*20+ 2] =   81.00;
  daa[ 7*20+ 3] =  130.00; daa[ 7*20+ 4] =   59.00; daa[ 7*20+ 5] =   26.00; daa[ 7*20+ 6] =  119.00;
  daa[ 8*20+ 0] =   27.00; daa[ 8*20+ 1] =  328.00; daa[ 8*20+ 2] =  391.00; daa[ 8*20+ 3] =  112.00;
  daa[ 8*20+ 4] =   69.00; daa[ 8*20+ 5] =  597.00; daa[ 8*20+ 6] =   26.00; daa[ 8*20+ 7] =   23.00;
  daa[ 9*20+ 0] =   36.00; daa[ 9*20+ 1] =   22.00; daa[ 9*20+ 2] =   47.00; daa[ 9*20+ 3] =   11.00;
  daa[ 9*20+ 4] =   17.00; daa[ 9*20+ 5] =    9.00; daa[ 9*20+ 6] =   12.00; daa[ 9*20+ 7] =    6.00;
  daa[ 9*20+ 8] =   16.00; daa[10*20+ 0] =   30.00; daa[10*20+ 1] =   38.00; daa[10*20+ 2] =   12.00;
  daa[10*20+ 3] =    7.00; daa[10*20+ 4] =   23.00; daa[10*20+ 5] =   72.00; daa[10*20+ 6] =    9.00;
  daa[10*20+ 7] =    6.00; daa[10*20+ 8] =   56.00; daa[10*20+ 9] =  229.00; daa[11*20+ 0] =   35.00;
  daa[11*20+ 1] =  646.00; daa[11*20+ 2] =  263.00; daa[11*20+ 3] =   26.00; daa[11*20+ 4] =    7.00;
  daa[11*20+ 5] =  292.00; daa[11*20+ 6] =  181.00; daa[11*20+ 7] =   27.00; daa[11*20+ 8] =   45.00;
  daa[11*20+ 9] =   21.00; daa[11*20+10] =   14.00; daa[12*20+ 0] =   54.00; daa[12*20+ 1] =   44.00;
  daa[12*20+ 2] =   30.00; daa[12*20+ 3] =   15.00; daa[12*20+ 4] =   31.00; daa[12*20+ 5] =   43.00;
  daa[12*20+ 6] =   18.00; daa[12*20+ 7] =   14.00; daa[12*20+ 8] =   33.00; daa[12*20+ 9] =  479.00;
  daa[12*20+10] =  388.00; daa[12*20+11] =   65.00; daa[13*20+ 0] =   15.00; daa[13*20+ 1] =    5.00;
  daa[13*20+ 2] =   10.00; daa[13*20+ 3] =    4.00; daa[13*20+ 4] =   78.00; daa[13*20+ 5] =    4.00;
  daa[13*20+ 6] =    5.00; daa[13*20+ 7] =    5.00; daa[13*20+ 8] =   40.00; daa[13*20+ 9] =   89.00;
  daa[13*20+10] =  248.00; daa[13*20+11] =    4.00; daa[13*20+12] =   43.00; daa[14*20+ 0] =  194.00;
  daa[14*20+ 1] =   74.00; daa[14*20+ 2] =   15.00; daa[14*20+ 3] =   15.00; daa[14*20+ 4] =   14.00;
  daa[14*20+ 5] =  164.00; daa[14*20+ 6] =   18.00; daa[14*20+ 7] =   24.00; daa[14*20+ 8] =  115.00;
  daa[14*20+ 9] =   10.00; daa[14*20+10] =  102.00; daa[14*20+11] =   21.00; daa[14*20+12] =   16.00;
  daa[14*20+13] =   17.00; daa[15*20+ 0] =  378.00; daa[15*20+ 1] =  101.00; daa[15*20+ 2] =  503.00;
  daa[15*20+ 3] =   59.00; daa[15*20+ 4] =  223.00; daa[15*20+ 5] =   53.00; daa[15*20+ 6] =   30.00;
  daa[15*20+ 7] =  201.00; daa[15*20+ 8] =   73.00; daa[15*20+ 9] =   40.00; daa[15*20+10] =   59.00;
  daa[15*20+11] =   47.00; daa[15*20+12] =   29.00; daa[15*20+13] =   92.00; daa[15*20+14] =  285.00;
  daa[16*20+ 0] =  475.00; daa[16*20+ 1] =   64.00; daa[16*20+ 2] =  232.00; daa[16*20+ 3] =   38.00;
  daa[16*20+ 4] =   42.00; daa[16*20+ 5] =   51.00; daa[16*20+ 6] =   32.00; daa[16*20+ 7] =   33.00;
  daa[16*20+ 8] =   46.00; daa[16*20+ 9] =  245.00; daa[16*20+10] =   25.00; daa[16*20+11] =  103.00;
  daa[16*20+12] =  226.00; daa[16*20+13] =   12.00; daa[16*20+14] =  118.00; daa[16*20+15] =  477.00;
  daa[17*20+ 0] =    9.00; daa[17*20+ 1] =  126.00; daa[17*20+ 2] =    8.00; daa[17*20+ 3] =    4.00;
  daa[17*20+ 4] =  115.00; daa[17*20+ 5] =   18.00; daa[17*20+ 6] =   10.00; daa[17*20+ 7] =   55.00;
  daa[17*20+ 8] =    8.00; daa[17*20+ 9] =    9.00; daa[17*20+10] =   52.00; daa[17*20+11] =   10.00;
  daa[17*20+12] =   24.00; daa[17*20+13] =   53.00; daa[17*20+14] =    6.00; daa[17*20+15] =   35.00;
  daa[17*20+16] =   12.00; daa[18*20+ 0] =   11.00; daa[18*20+ 1] =   20.00; daa[18*20+ 2] =   70.00;
  daa[18*20+ 3] =   46.00; daa[18*20+ 4] =  209.00; daa[18*20+ 5] =   24.00; daa[18*20+ 6] =    7.00;
  daa[18*20+ 7] =    8.00; daa[18*20+ 8] =  573.00; daa[18*20+ 9] =   32.00; daa[18*20+10] =   24.00;
  daa[18*20+11] =    8.00; daa[18*20+12] =   18.00; daa[18*20+13] =  536.00; daa[18*20+14] =   10.00;
  daa[18*20+15] =   63.00; daa[18*20+16] =   21.00; daa[18*20+17] =   71.00; daa[19*20+ 0] =  298.00;
  daa[19*20+ 1] =   17.00; daa[19*20+ 2] =   16.00; daa[19*20+ 3] =   31.00; daa[19*20+ 4] =   62.00;
  daa[19*20+ 5] =   20.00; daa[19*20+ 6] =   45.00; daa[19*20+ 7] =   47.00; daa[19*20+ 8] =   11.00;
  daa[19*20+ 9] =  961.00; daa[19*20+10] =  180.00; daa[19*20+11] =   14.00; daa[19*20+12] =  323.00;
  daa[19*20+13] =   62.00; daa[19*20+14] =   23.00; daa[19*20+15] =   38.00; daa[19*20+16] =  112.00;
  daa[19*20+17] =   25.00; daa[19*20+18] =   16.00;

  for (i=0; i<naa; i++)  for (j=0; j<i; j++)  daa[j*naa+i] = daa[i*naa+j];

  pi[ 0] = 0.076748; pi[ 1] = 0.051691; pi[ 2] = 0.042645; pi[ 3] = 0.051544;
  pi[ 4] = 0.019803; pi[ 5] = 0.040752; pi[ 6] = 0.061830; pi[ 7] = 0.073152;
  pi[ 8] = 0.022944; pi[ 9] = 0.053761; pi[10] = 0.091904; pi[11] = 0.058676;
  pi[12] = 0.023826; pi[13] = 0.040126; pi[14] = 0.050901; pi[15] = 0.068765;
  pi[16] = 0.058565; pi[17] = 0.014261; pi[18] = 0.032102; pi[19] = 0.066005;

 return 1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


int Init_Qmat_MtREV(phydbl *daa, phydbl *pi)
{
  /* J. Adachi and M. Hasegawa, ``Model of amino acid substitution in proteins
     encoded by mitochondrial DNA'' J. Mol. Evol. 42, 459 (1996) */

  int i,j,naa;

  naa = 20;


  daa[ 1*20+ 0] =   23.18; daa[ 2*20+ 0] =   26.95; daa[ 2*20+ 1] =   13.24; daa[ 3*20+ 0] =   17.67;
  daa[ 3*20+ 1] =    1.90; daa[ 3*20+ 2] =  794.38; daa[ 4*20+ 0] =   59.93; daa[ 4*20+ 1] =  103.33;
  daa[ 4*20+ 2] =   58.94; daa[ 4*20+ 3] =    1.90; daa[ 5*20+ 0] =    1.90; daa[ 5*20+ 1] =  220.99;
  daa[ 5*20+ 2] =  173.56; daa[ 5*20+ 3] =   55.28; daa[ 5*20+ 4] =   75.24; daa[ 6*20+ 0] =    9.77;
  daa[ 6*20+ 1] =    1.90; daa[ 6*20+ 2] =   63.05; daa[ 6*20+ 3] =  583.55; daa[ 6*20+ 4] =    1.90;
  daa[ 6*20+ 5] =  313.56; daa[ 7*20+ 0] =  120.71; daa[ 7*20+ 1] =   23.03; daa[ 7*20+ 2] =   53.30;
  daa[ 7*20+ 3] =   56.77; daa[ 7*20+ 4] =   30.71; daa[ 7*20+ 5] =    6.75; daa[ 7*20+ 6] =   28.28;
  daa[ 8*20+ 0] =   13.90; daa[ 8*20+ 1] =  165.23; daa[ 8*20+ 2] =  496.13; daa[ 8*20+ 3] =  113.99;
  daa[ 8*20+ 4] =  141.49; daa[ 8*20+ 5] =  582.40; daa[ 8*20+ 6] =   49.12; daa[ 8*20+ 7] =    1.90;
  daa[ 9*20+ 0] =   96.49; daa[ 9*20+ 1] =    1.90; daa[ 9*20+ 2] =   27.10; daa[ 9*20+ 3] =    4.34;
  daa[ 9*20+ 4] =   62.73; daa[ 9*20+ 5] =    8.34; daa[ 9*20+ 6] =    3.31; daa[ 9*20+ 7] =    5.98;
  daa[ 9*20+ 8] =   12.26; daa[10*20+ 0] =   25.46; daa[10*20+ 1] =   15.58; daa[10*20+ 2] =   15.16;
  daa[10*20+ 3] =    1.90; daa[10*20+ 4] =   25.65; daa[10*20+ 5] =   39.70; daa[10*20+ 6] =    1.90;
  daa[10*20+ 7] =    2.41; daa[10*20+ 8] =   11.49; daa[10*20+ 9] =  329.09; daa[11*20+ 0] =    8.36;
  daa[11*20+ 1] =  141.40; daa[11*20+ 2] =  608.70; daa[11*20+ 3] =    2.31; daa[11*20+ 4] =    1.90;
  daa[11*20+ 5] =  465.58; daa[11*20+ 6] =  313.86; daa[11*20+ 7] =   22.73; daa[11*20+ 8] =  127.67;
  daa[11*20+ 9] =   19.57; daa[11*20+10] =   14.88; daa[12*20+ 0] =  141.88; daa[12*20+ 1] =    1.90;
  daa[12*20+ 2] =   65.41; daa[12*20+ 3] =    1.90; daa[12*20+ 4] =    6.18; daa[12*20+ 5] =   47.37;
  daa[12*20+ 6] =    1.90; daa[12*20+ 7] =    1.90; daa[12*20+ 8] =   11.97; daa[12*20+ 9] =  517.98;
  daa[12*20+10] =  537.53; daa[12*20+11] =   91.37; daa[13*20+ 0] =    6.37; daa[13*20+ 1] =    4.69;
  daa[13*20+ 2] =   15.20; daa[13*20+ 3] =    4.98; daa[13*20+ 4] =   70.80; daa[13*20+ 5] =   19.11;
  daa[13*20+ 6] =    2.67; daa[13*20+ 7] =    1.90; daa[13*20+ 8] =   48.16; daa[13*20+ 9] =   84.67;
  daa[13*20+10] =  216.06; daa[13*20+11] =    6.44; daa[13*20+12] =   90.82; daa[14*20+ 0] =   54.31;
  daa[14*20+ 1] =   23.64; daa[14*20+ 2] =   73.31; daa[14*20+ 3] =   13.43; daa[14*20+ 4] =   31.26;
  daa[14*20+ 5] =  137.29; daa[14*20+ 6] =   12.83; daa[14*20+ 7] =    1.90; daa[14*20+ 8] =   60.97;
  daa[14*20+ 9] =   20.63; daa[14*20+10] =   40.10; daa[14*20+11] =   50.10; daa[14*20+12] =   18.84;
  daa[14*20+13] =   17.31; daa[15*20+ 0] =  387.86; daa[15*20+ 1] =    6.04; daa[15*20+ 2] =  494.39;
  daa[15*20+ 3] =   69.02; daa[15*20+ 4] =  277.05; daa[15*20+ 5] =   54.11; daa[15*20+ 6] =   54.71;
  daa[15*20+ 7] =  125.93; daa[15*20+ 8] =   77.46; daa[15*20+ 9] =   47.70; daa[15*20+10] =   73.61;
  daa[15*20+11] =  105.79; daa[15*20+12] =  111.16; daa[15*20+13] =   64.29; daa[15*20+14] =  169.90;
  daa[16*20+ 0] =  480.72; daa[16*20+ 1] =    2.08; daa[16*20+ 2] =  238.46; daa[16*20+ 3] =   28.01;
  daa[16*20+ 4] =  179.97; daa[16*20+ 5] =   94.93; daa[16*20+ 6] =   14.82; daa[16*20+ 7] =   11.17;
  daa[16*20+ 8] =   44.78; daa[16*20+ 9] =  368.43; daa[16*20+10] =  126.40; daa[16*20+11] =  136.33;
  daa[16*20+12] =  528.17; daa[16*20+13] =   33.85; daa[16*20+14] =  128.22; daa[16*20+15] =  597.21;
  daa[17*20+ 0] =    1.90; daa[17*20+ 1] =   21.95; daa[17*20+ 2] =   10.68; daa[17*20+ 3] =   19.86;
  daa[17*20+ 4] =   33.60; daa[17*20+ 5] =    1.90; daa[17*20+ 6] =    1.90; daa[17*20+ 7] =   10.92;
  daa[17*20+ 8] =    7.08; daa[17*20+ 9] =    1.90; daa[17*20+10] =   32.44; daa[17*20+11] =   24.00;
  daa[17*20+12] =   21.71; daa[17*20+13] =    7.84; daa[17*20+14] =    4.21; daa[17*20+15] =   38.58;
  daa[17*20+16] =    9.99; daa[18*20+ 0] =    6.48; daa[18*20+ 1] =    1.90; daa[18*20+ 2] =  191.36;
  daa[18*20+ 3] =   21.21; daa[18*20+ 4] =  254.77; daa[18*20+ 5] =   38.82; daa[18*20+ 6] =   13.12;
  daa[18*20+ 7] =    3.21; daa[18*20+ 8] =  670.14; daa[18*20+ 9] =   25.01; daa[18*20+10] =   44.15;
  daa[18*20+11] =   51.17; daa[18*20+12] =   39.96; daa[18*20+13] =  465.58; daa[18*20+14] =   16.21;
  daa[18*20+15] =   64.92; daa[18*20+16] =   38.73; daa[18*20+17] =   26.25; daa[19*20+ 0] =  195.06;
  daa[19*20+ 1] =    7.64; daa[19*20+ 2] =    1.90; daa[19*20+ 3] =    1.90; daa[19*20+ 4] =    1.90;
  daa[19*20+ 5] =   19.00; daa[19*20+ 6] =   21.14; daa[19*20+ 7] =    2.53; daa[19*20+ 8] =    1.90;
  daa[19*20+ 9] = 1222.94; daa[19*20+10] =   91.67; daa[19*20+11] =    1.90; daa[19*20+12] =  387.54;
  daa[19*20+13] =    6.35; daa[19*20+14] =    8.23; daa[19*20+15] =    1.90; daa[19*20+16] =  204.54;
  daa[19*20+17] =    5.37; daa[19*20+18] =    1.90;

  for (i=0; i<naa; i++)  for (j=0; j<i; j++)  daa[j*naa+i] = daa[i*naa+j];

  pi[ 0] = 0.072000; pi[ 1] = 0.019000; pi[ 2] = 0.039000; pi[ 3] = 0.019000;
  pi[ 4] = 0.006000; pi[ 5] = 0.025000; pi[ 6] = 0.024000; pi[ 7] = 0.056000;
  pi[ 8] = 0.028000; pi[ 9] = 0.088000; pi[10] = 0.169000; pi[11] = 0.023000;
  pi[12] = 0.054000; pi[13] = 0.061000; pi[14] = 0.054000; pi[15] = 0.072000;
  pi[16] = 0.086000; pi[17] = 0.029000; pi[18] = 0.033000; pi[19] = 0.043000;

  return 1;
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int Init_Qmat_FLU(phydbl *daa, phydbl *pi)
{
  int i,j,naa;

  /* FLU model 
   *  Cuong Cao Dang, Quang Si Le2, Olivier Gascuel and Vinh Sy Le
   * 'FLU, an amino acid substitution model for influenza proteins' 
   */

  naa = 20;
  
  daa[ 1*20+ 0] =   0.138659;  daa[ 2*20+ 0] =   0.053367;  daa[ 2*20+ 1] =   0.161001;  daa[ 3*20+ 0] =   0.584852;  
  daa[ 3*20+ 1] =   0.006772;  daa[ 3*20+ 2] =   7.737393;  daa[ 4*20+ 0] =   0.026447;  daa[ 4*20+ 1] =   0.167207;  
  daa[ 4*20+ 2] =   0.000013;  daa[ 4*20+ 3] =   0.014132;  daa[ 5*20+ 0] =   0.353754;  daa[ 5*20+ 1] =   3.292717;  
  daa[ 5*20+ 2] =   0.530643;  daa[ 5*20+ 3] =   0.145469;  daa[ 5*20+ 4] =   0.002547;  daa[ 6*20+ 0] =   1.484235;  
  daa[ 6*20+ 1] =   0.124898;  daa[ 6*20+ 2] =   0.061652;  daa[ 6*20+ 3] =   5.370511;  daa[ 6*20+ 4] =   0.000000;  
  daa[ 6*20+ 5] =   1.195629;  daa[ 7*20+ 0] =   1.132313;  daa[ 7*20+ 1] =   1.190624;  daa[ 7*20+ 2] =   0.322525;  
  daa[ 7*20+ 3] =   1.934833;  daa[ 7*20+ 4] =   0.116941;  daa[ 7*20+ 5] =   0.108051;  daa[ 7*20+ 6] =   1.593099;  
  daa[ 8*20+ 0] =   0.214758;  daa[ 8*20+ 1] =   1.879570;  daa[ 8*20+ 2] =   1.387096;  daa[ 8*20+ 3] =   0.887571;  
  daa[ 8*20+ 4] =   0.021845;  daa[ 8*20+ 5] =   5.330313;  daa[ 8*20+ 6] =   0.256492;  daa[ 8*20+ 7] =   0.058775;  
  daa[ 9*20+ 0] =   0.149927;  daa[ 9*20+ 1] =   0.246117;  daa[ 9*20+ 2] =   0.218572;  daa[ 9*20+ 3] =   0.014086;  
  daa[ 9*20+ 4] =   0.001112;  daa[ 9*20+ 5] =   0.028840;  daa[ 9*20+ 6] =   0.014211;  daa[ 9*20+ 7] =   0.000016;  
  daa[ 9*20+ 8] =   0.243190;  daa[10*20+ 0] =   0.023117;  daa[10*20+ 1] =   0.296046;  daa[10*20+ 2] =   0.000836;  
  daa[10*20+ 3] =   0.005731;  daa[10*20+ 4] =   0.005614;  daa[10*20+ 5] =   1.020367;  daa[10*20+ 6] =   0.016500;  
  daa[10*20+ 7] =   0.006516;  daa[10*20+ 8] =   0.321612;  daa[10*20+ 9] =   3.512072;  daa[11*20+ 0] =   0.474334;  
  daa[11*20+ 1] =  15.300097;  daa[11*20+ 2] =   2.646848;  daa[11*20+ 3] =   0.290043;  daa[11*20+ 4] =   0.000004;  
  daa[11*20+ 5] =   2.559587;  daa[11*20+ 6] =   3.881489;  daa[11*20+ 7] =   0.264149;  daa[11*20+ 8] =   0.347303;  
  daa[11*20+ 9] =   0.227708;  daa[11*20+10] =   0.129224;  daa[12*20+ 0] =   0.058745;  daa[12*20+ 1] =   0.890162;  
  daa[12*20+ 2] =   0.005252;  daa[12*20+ 3] =   0.041763;  daa[12*20+ 4] =   0.111457;  daa[12*20+ 5] =   0.190259;  
  daa[12*20+ 6] =   0.313974;  daa[12*20+ 7] =   0.001500;  daa[12*20+ 8] =   0.001274;  daa[12*20+ 9] =   9.017954;  
  daa[12*20+10] =   6.746936;  daa[12*20+11] =   1.331292;  daa[13*20+ 0] =   0.080491;  daa[13*20+ 1] =   0.016055;  
  daa[13*20+ 2] =   0.000836;  daa[13*20+ 3] =   0.000001;  daa[13*20+ 4] =   0.104054;  daa[13*20+ 5] =   0.032681;  
  daa[13*20+ 6] =   0.001004;  daa[13*20+ 7] =   0.001237;  daa[13*20+ 8] =   0.119029;  daa[13*20+ 9] =   1.463357;  
  daa[13*20+10] =   2.986800;  daa[13*20+11] =   0.319896;  daa[13*20+12] =   0.279911;  daa[14*20+ 0] =   0.659311;  
  daa[14*20+ 1] =   0.154027;  daa[14*20+ 2] =   0.036442;  daa[14*20+ 3] =   0.188539;  daa[14*20+ 4] =   0.000000;  
  daa[14*20+ 5] =   0.712770;  daa[14*20+ 6] =   0.319559;  daa[14*20+ 7] =   0.038632;  daa[14*20+ 8] =   0.924467;  
  daa[14*20+ 9] =   0.080543;  daa[14*20+10] =   0.634309;  daa[14*20+11] =   0.195751;  daa[14*20+12] =   0.056869;  
  daa[14*20+13] =   0.007132;  daa[15*20+ 0] =   3.011345;  daa[15*20+ 1] =   0.950138;  daa[15*20+ 2] =   3.881311;  
  daa[15*20+ 3] =   0.338372;  daa[15*20+ 4] =   0.336263;  daa[15*20+ 5] =   0.487822;  daa[15*20+ 6] =   0.307140;  
  daa[15*20+ 7] =   1.585647;  daa[15*20+ 8] =   0.580704;  daa[15*20+ 9] =   0.290381;  daa[15*20+10] =   0.570767;  
  daa[15*20+11] =   0.283808;  daa[15*20+12] =   0.007027;  daa[15*20+13] =   0.996686;  daa[15*20+14] =   2.087385;  
  daa[16*20+ 0] =   5.418298;  daa[16*20+ 1] =   0.183077;  daa[16*20+ 2] =   2.140332;  daa[16*20+ 3] =   0.135481;  
  daa[16*20+ 4] =   0.011975;  daa[16*20+ 5] =   0.602341;  daa[16*20+ 6] =   0.280125;  daa[16*20+ 7] =   0.018808;  
  daa[16*20+ 8] =   0.368714;  daa[16*20+ 9] =   2.904052;  daa[16*20+10] =   0.044926;  daa[16*20+11] =   1.526964;  
  daa[16*20+12] =   2.031511;  daa[16*20+13] =   0.000135;  daa[16*20+14] =   0.542251;  daa[16*20+15] =   2.206860;  
  daa[17*20+ 0] =   0.195966;  daa[17*20+ 1] =   1.369429;  daa[17*20+ 2] =   0.000536;  daa[17*20+ 3] =   0.000015;  
  daa[17*20+ 4] =   0.094107;  daa[17*20+ 5] =   0.044021;  daa[17*20+ 6] =   0.155245;  daa[17*20+ 7] =   0.196486;  
  daa[17*20+ 8] =   0.022373;  daa[17*20+ 9] =   0.032132;  daa[17*20+10] =   0.431278;  daa[17*20+11] =   0.000050;  
  daa[17*20+12] =   0.070460;  daa[17*20+13] =   0.814753;  daa[17*20+14] =   0.000431;  daa[17*20+15] =   0.099836;  
  daa[17*20+16] =   0.207066;  daa[18*20+ 0] =   0.018289;  daa[18*20+ 1] =   0.099855;  daa[18*20+ 2] =   0.373102;  
  daa[18*20+ 3] =   0.525399;  daa[18*20+ 4] =   0.601692;  daa[18*20+ 5] =   0.072206;  daa[18*20+ 6] =   0.104093;  
  daa[18*20+ 7] =   0.074815;  daa[18*20+ 8] =   6.448954;  daa[18*20+ 9] =   0.273934;  daa[18*20+10] =   0.340058;  
  daa[18*20+11] =   0.012416;  daa[18*20+12] =   0.874272;  daa[18*20+13] =   5.393924;  daa[18*20+14] =   0.000182;  
  daa[18*20+15] =   0.392552;  daa[18*20+16] =   0.124898;  daa[18*20+17] =   0.427755;  daa[19*20+ 0] =   3.532005;  
  daa[19*20+ 1] =   0.103964;  daa[19*20+ 2] =   0.010258;  daa[19*20+ 3] =   0.297124;  daa[19*20+ 4] =   0.054905;  
  daa[19*20+ 5] =   0.406698;  daa[19*20+ 6] =   0.285048;  daa[19*20+ 7] =   0.337230;  daa[19*20+ 8] =   0.098631;  
  daa[19*20+ 9] =  14.394052;  daa[19*20+10] =   0.890599;  daa[19*20+11] =   0.073128;  daa[19*20+12] =   4.904842;  
  daa[19*20+13] =   0.592588;  daa[19*20+14] =   0.058972;  daa[19*20+15] =   0.088256;  daa[19*20+16] =   0.654109;  
  daa[19*20+17] =   0.256900;  daa[19*20+18] =   0.167582;  
  
  for (i=0; i<naa; i++)  for (j=0; j<i; j++)  daa[j*naa+i] = daa[i*naa+j];
  
  pi[0] = 0.047072; pi[1] = 0.050910; pi[2] = 0.074214; pi[3] = 0.047860; pi[4] = 0.025022; pi[5] = 0.033304; pi[6] = 0.054587; pi[7] = 0.076373; pi[8] = 0.019964; pi[9] = 0.067134; pi[10] = 0.071498; pi[11] = 0.056785; pi[12] = 0.018151; pi[13] = 0.030496; pi[14] = 0.050656; pi[15] = 0.088409; pi[16] = 0.074339; pi[17] = 0.018524; pi[18] = 0.031474; pi[19] = 0.063229; 

  return 1;
}
  
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int Init_Qmat_LG(phydbl *daa, phydbl *pi)
{
  int i,j,naa;

  /* LG model
   * Si Quang Le & Olivier Gascuel
   * 'An improved general amino-acid replacement matrix'
   */

  naa = 20;

  daa[ 1*20+ 0] =   0.425093;  daa[ 2*20+ 0] =   0.276818;  daa[ 2*20+ 1] =   0.751878;  daa[ 3*20+ 0] =   0.395144;
  daa[ 3*20+ 1] =   0.123954;  daa[ 3*20+ 2] =   5.076149;  daa[ 4*20+ 0] =   2.489084;  daa[ 4*20+ 1] =   0.534551;
  daa[ 4*20+ 2] =   0.528768;  daa[ 4*20+ 3] =   0.062556;  daa[ 5*20+ 0] =   0.969894;  daa[ 5*20+ 1] =   2.807908;
  daa[ 5*20+ 2] =   1.695752;  daa[ 5*20+ 3] =   0.523386;  daa[ 5*20+ 4] =   0.084808;  daa[ 6*20+ 0] =   1.038545;
  daa[ 6*20+ 1] =   0.363970;  daa[ 6*20+ 2] =   0.541712;  daa[ 6*20+ 3] =   5.243870;  daa[ 6*20+ 4] =   0.003499;
  daa[ 6*20+ 5] =   4.128591;  daa[ 7*20+ 0] =   2.066040;  daa[ 7*20+ 1] =   0.390192;  daa[ 7*20+ 2] =   1.437645;
  daa[ 7*20+ 3] =   0.844926;  daa[ 7*20+ 4] =   0.569265;  daa[ 7*20+ 5] =   0.267959;  daa[ 7*20+ 6] =   0.348847;
  daa[ 8*20+ 0] =   0.358858;  daa[ 8*20+ 1] =   2.426601;  daa[ 8*20+ 2] =   4.509238;  daa[ 8*20+ 3] =   0.927114;
  daa[ 8*20+ 4] =   0.640543;  daa[ 8*20+ 5] =   4.813505;  daa[ 8*20+ 6] =   0.423881;  daa[ 8*20+ 7] =   0.311484;
  daa[ 9*20+ 0] =   0.149830;  daa[ 9*20+ 1] =   0.126991;  daa[ 9*20+ 2] =   0.191503;  daa[ 9*20+ 3] =   0.010690;
  daa[ 9*20+ 4] =   0.320627;  daa[ 9*20+ 5] =   0.072854;  daa[ 9*20+ 6] =   0.044265;  daa[ 9*20+ 7] =   0.008705;
  daa[ 9*20+ 8] =   0.108882;  daa[10*20+ 0] =   0.395337;  daa[10*20+ 1] =   0.301848;  daa[10*20+ 2] =   0.068427;
  daa[10*20+ 3] =   0.015076;  daa[10*20+ 4] =   0.594007;  daa[10*20+ 5] =   0.582457;  daa[10*20+ 6] =   0.069673;
  daa[10*20+ 7] =   0.044261;  daa[10*20+ 8] =   0.366317;  daa[10*20+ 9] =   4.145067;  daa[11*20+ 0] =   0.536518;
  daa[11*20+ 1] =   6.326067;  daa[11*20+ 2] =   2.145078;  daa[11*20+ 3] =   0.282959;  daa[11*20+ 4] =   0.013266;
  daa[11*20+ 5] =   3.234294;  daa[11*20+ 6] =   1.807177;  daa[11*20+ 7] =   0.296636;  daa[11*20+ 8] =   0.697264;
  daa[11*20+ 9] =   0.159069;  daa[11*20+10] =   0.137500;  daa[12*20+ 0] =   1.124035;  daa[12*20+ 1] =   0.484133;
  daa[12*20+ 2] =   0.371004;  daa[12*20+ 3] =   0.025548;  daa[12*20+ 4] =   0.893680;  daa[12*20+ 5] =   1.672569;
  daa[12*20+ 6] =   0.173735;  daa[12*20+ 7] =   0.139538;  daa[12*20+ 8] =   0.442472;  daa[12*20+ 9] =   4.273607;
  daa[12*20+10] =   6.312358;  daa[12*20+11] =   0.656604;  daa[13*20+ 0] =   0.253701;  daa[13*20+ 1] =   0.052722;
  daa[13*20+ 2] =   0.089525;  daa[13*20+ 3] =   0.017416;  daa[13*20+ 4] =   1.105251;  daa[13*20+ 5] =   0.035855;
  daa[13*20+ 6] =   0.018811;  daa[13*20+ 7] =   0.089586;  daa[13*20+ 8] =   0.682139;  daa[13*20+ 9] =   1.112727;
  daa[13*20+10] =   2.592692;  daa[13*20+11] =   0.023918;  daa[13*20+12] =   1.798853;  daa[14*20+ 0] =   1.177651;
  daa[14*20+ 1] =   0.332533;  daa[14*20+ 2] =   0.161787;  daa[14*20+ 3] =   0.394456;  daa[14*20+ 4] =   0.075382;
  daa[14*20+ 5] =   0.624294;  daa[14*20+ 6] =   0.419409;  daa[14*20+ 7] =   0.196961;  daa[14*20+ 8] =   0.508851;
  daa[14*20+ 9] =   0.078281;  daa[14*20+10] =   0.249060;  daa[14*20+11] =   0.390322;  daa[14*20+12] =   0.099849;
  daa[14*20+13] =   0.094464;  daa[15*20+ 0] =   4.727182;  daa[15*20+ 1] =   0.858151;  daa[15*20+ 2] =   4.008358;
  daa[15*20+ 3] =   1.240275;  daa[15*20+ 4] =   2.784478;  daa[15*20+ 5] =   1.223828;  daa[15*20+ 6] =   0.611973;
  daa[15*20+ 7] =   1.739990;  daa[15*20+ 8] =   0.990012;  daa[15*20+ 9] =   0.064105;  daa[15*20+10] =   0.182287;
  daa[15*20+11] =   0.748683;  daa[15*20+12] =   0.346960;  daa[15*20+13] =   0.361819;  daa[15*20+14] =   1.338132;
  daa[16*20+ 0] =   2.139501;  daa[16*20+ 1] =   0.578987;  daa[16*20+ 2] =   2.000679;  daa[16*20+ 3] =   0.425860;
  daa[16*20+ 4] =   1.143480;  daa[16*20+ 5] =   1.080136;  daa[16*20+ 6] =   0.604545;  daa[16*20+ 7] =   0.129836;
  daa[16*20+ 8] =   0.584262;  daa[16*20+ 9] =   1.033739;  daa[16*20+10] =   0.302936;  daa[16*20+11] =   1.136863;
  daa[16*20+12] =   2.020366;  daa[16*20+13] =   0.165001;  daa[16*20+14] =   0.571468;  daa[16*20+15] =   6.472279;
  daa[17*20+ 0] =   0.180717;  daa[17*20+ 1] =   0.593607;  daa[17*20+ 2] =   0.045376;  daa[17*20+ 3] =   0.029890;
  daa[17*20+ 4] =   0.670128;  daa[17*20+ 5] =   0.236199;  daa[17*20+ 6] =   0.077852;  daa[17*20+ 7] =   0.268491;
  daa[17*20+ 8] =   0.597054;  daa[17*20+ 9] =   0.111660;  daa[17*20+10] =   0.619632;  daa[17*20+11] =   0.049906;
  daa[17*20+12] =   0.696175;  daa[17*20+13] =   2.457121;  daa[17*20+14] =   0.095131;  daa[17*20+15] =   0.248862;
  daa[17*20+16] =   0.140825;  daa[18*20+ 0] =   0.218959;  daa[18*20+ 1] =   0.314440;  daa[18*20+ 2] =   0.612025;
  daa[18*20+ 3] =   0.135107;  daa[18*20+ 4] =   1.165532;  daa[18*20+ 5] =   0.257336;  daa[18*20+ 6] =   0.120037;
  daa[18*20+ 7] =   0.054679;  daa[18*20+ 8] =   5.306834;  daa[18*20+ 9] =   0.232523;  daa[18*20+10] =   0.299648;
  daa[18*20+11] =   0.131932;  daa[18*20+12] =   0.481306;  daa[18*20+13] =   7.803902;  daa[18*20+14] =   0.089613;
  daa[18*20+15] =   0.400547;  daa[18*20+16] =   0.245841;  daa[18*20+17] =   3.151815;  daa[19*20+ 0] =   2.547870;
  daa[19*20+ 1] =   0.170887;  daa[19*20+ 2] =   0.083688;  daa[19*20+ 3] =   0.037967;  daa[19*20+ 4] =   1.959291;
  daa[19*20+ 5] =   0.210332;  daa[19*20+ 6] =   0.245034;  daa[19*20+ 7] =   0.076701;  daa[19*20+ 8] =   0.119013;
  daa[19*20+ 9] =  10.649107;  daa[19*20+10] =   1.702745;  daa[19*20+11] =   0.185202;  daa[19*20+12] =   1.898718;
  daa[19*20+13] =   0.654683;  daa[19*20+14] =   0.296501;  daa[19*20+15] =   0.098369;  daa[19*20+16] =   2.188158;
  daa[19*20+17] =   0.189510;  daa[19*20+18] =   0.249313;


  for (i=0; i<naa; i++)  for (j=0; j<i; j++)  daa[j*naa+i] = daa[i*naa+j];

  pi[0] = 0.079066; pi[1] = 0.055941; pi[2] = 0.041977; pi[3] = 0.053052;
  pi[4] = 0.012937; pi[5] = 0.040767; pi[6] = 0.071586; pi[7] = 0.057337;
  pi[8] = 0.022355; pi[9] = 0.062157; pi[10] = 0.099081; pi[11] = 0.064600;
  pi[12] = 0.022951; pi[13] = 0.042302; pi[14] = 0.044040; pi[15] = 0.061197;
  pi[16] = 0.053287; pi[17] = 0.012066; pi[18] = 0.034155; pi[19] = 0.069147;

  return 1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int Init_Qmat_WAG(phydbl *daa, phydbl *pi)
{
  int i,j,naa;

  /* WAG's model data
   * Simon Whelan and Nick Goldman
   * 'A general empirical model of protein evolution derived from multiple
   *  protein families using a maximum-likelihood approach'
   * MBE (2001) 18:691-699
   */


  naa = 20;

  daa[ 1*20+ 0] =  55.15710; daa[ 2*20+ 0] =  50.98480; daa[ 2*20+ 1] =  63.53460;
  daa[ 3*20+ 0] =  73.89980; daa[ 3*20+ 1] =  14.73040; daa[ 3*20+ 2] = 542.94200;
  daa[ 4*20+ 0] = 102.70400; daa[ 4*20+ 1] =  52.81910; daa[ 4*20+ 2] =  26.52560;
  daa[ 4*20+ 3] =   3.02949; daa[ 5*20+ 0] =  90.85980; daa[ 5*20+ 1] = 303.55000;
  daa[ 5*20+ 2] = 154.36400; daa[ 5*20+ 3] =  61.67830; daa[ 5*20+ 4] =   9.88179;
  daa[ 6*20+ 0] = 158.28500; daa[ 6*20+ 1] =  43.91570; daa[ 6*20+ 2] =  94.71980;
  daa[ 6*20+ 3] = 617.41600; daa[ 6*20+ 4] =   2.13520; daa[ 6*20+ 5] = 546.94700;
  daa[ 7*20+ 0] = 141.67200; daa[ 7*20+ 1] =  58.46650; daa[ 7*20+ 2] = 112.55600;
  daa[ 7*20+ 3] =  86.55840; daa[ 7*20+ 4] =  30.66740; daa[ 7*20+ 5] =  33.00520;
  daa[ 7*20+ 6] =  56.77170; daa[ 8*20+ 0] =  31.69540; daa[ 8*20+ 1] = 213.71500;
  daa[ 8*20+ 2] = 395.62900; daa[ 8*20+ 3] =  93.06760; daa[ 8*20+ 4] =  24.89720;
  daa[ 8*20+ 5] = 429.41100; daa[ 8*20+ 6] =  57.00250; daa[ 8*20+ 7] =  24.94100;
  daa[ 9*20+ 0] =  19.33350; daa[ 9*20+ 1] =  18.69790; daa[ 9*20+ 2] =  55.42360;
  daa[ 9*20+ 3] =   3.94370; daa[ 9*20+ 4] =  17.01350; daa[ 9*20+ 5] =  11.39170;
  daa[ 9*20+ 6] =  12.73950; daa[ 9*20+ 7] =   3.04501; daa[ 9*20+ 8] =  13.81900;
  daa[10*20+ 0] =  39.79150; daa[10*20+ 1] =  49.76710; daa[10*20+ 2] =  13.15280;
  daa[10*20+ 3] =   8.48047; daa[10*20+ 4] =  38.42870; daa[10*20+ 5] =  86.94890;
  daa[10*20+ 6] =  15.42630; daa[10*20+ 7] =   6.13037; daa[10*20+ 8] =  49.94620;
  daa[10*20+ 9] = 317.09700; daa[11*20+ 0] =  90.62650; daa[11*20+ 1] = 535.14200;
  daa[11*20+ 2] = 301.20100; daa[11*20+ 3] =  47.98550; daa[11*20+ 4] =   7.40339;
  daa[11*20+ 5] = 389.49000; daa[11*20+ 6] = 258.44300; daa[11*20+ 7] =  37.35580;
  daa[11*20+ 8] =  89.04320; daa[11*20+ 9] =  32.38320; daa[11*20+10] =  25.75550;
  daa[12*20+ 0] =  89.34960; daa[12*20+ 1] =  68.31620; daa[12*20+ 2] =  19.82210;
  daa[12*20+ 3] =  10.37540; daa[12*20+ 4] =  39.04820; daa[12*20+ 5] = 154.52600;
  daa[12*20+ 6] =  31.51240; daa[12*20+ 7] =  17.41000; daa[12*20+ 8] =  40.41410;
  daa[12*20+ 9] = 425.74600; daa[12*20+10] = 485.40200; daa[12*20+11] =  93.42760;
  daa[13*20+ 0] =  21.04940; daa[13*20+ 1] =  10.27110; daa[13*20+ 2] =   9.61621;
  daa[13*20+ 3] =   4.67304; daa[13*20+ 4] =  39.80200; daa[13*20+ 5] =   9.99208;
  daa[13*20+ 6] =   8.11339; daa[13*20+ 7] =   4.99310; daa[13*20+ 8] =  67.93710;
  daa[13*20+ 9] = 105.94700; daa[13*20+10] = 211.51700; daa[13*20+11] =   8.88360;
  daa[13*20+12] = 119.06300; daa[14*20+ 0] = 143.85500; daa[14*20+ 1] =  67.94890;
  daa[14*20+ 2] =  19.50810; daa[14*20+ 3] =  42.39840; daa[14*20+ 4] =  10.94040;
  daa[14*20+ 5] =  93.33720; daa[14*20+ 6] =  68.23550; daa[14*20+ 7] =  24.35700;
  daa[14*20+ 8] =  69.61980; daa[14*20+ 9] =   9.99288; daa[14*20+10] =  41.58440;
  daa[14*20+11] =  55.68960; daa[14*20+12] =  17.13290; daa[14*20+13] =  16.14440;
  daa[15*20+ 0] = 337.07900; daa[15*20+ 1] = 122.41900; daa[15*20+ 2] = 397.42300;
  daa[15*20+ 3] = 107.17600; daa[15*20+ 4] = 140.76600; daa[15*20+ 5] = 102.88700;
  daa[15*20+ 6] =  70.49390; daa[15*20+ 7] = 134.18200; daa[15*20+ 8] =  74.01690;
  daa[15*20+ 9] =  31.94400; daa[15*20+10] =  34.47390; daa[15*20+11] =  96.71300;
  daa[15*20+12] =  49.39050; daa[15*20+13] =  54.59310; daa[15*20+14] = 161.32800;
  daa[16*20+ 0] = 212.11100; daa[16*20+ 1] =  55.44130; daa[16*20+ 2] = 203.00600;
  daa[16*20+ 3] =  37.48660; daa[16*20+ 4] =  51.29840; daa[16*20+ 5] =  85.79280;
  daa[16*20+ 6] =  82.27650; daa[16*20+ 7] =  22.58330; daa[16*20+ 8] =  47.33070;
  daa[16*20+ 9] = 145.81600; daa[16*20+10] =  32.66220; daa[16*20+11] = 138.69800;
  daa[16*20+12] = 151.61200; daa[16*20+13] =  17.19030; daa[16*20+14] =  79.53840;
  daa[16*20+15] = 437.80200; daa[17*20+ 0] =  11.31330; daa[17*20+ 1] = 116.39200;
  daa[17*20+ 2] =   7.19167; daa[17*20+ 3] =  12.97670; daa[17*20+ 4] =  71.70700;
  daa[17*20+ 5] =  21.57370; daa[17*20+ 6] =  15.65570; daa[17*20+ 7] =  33.69830;
  daa[17*20+ 8] =  26.25690; daa[17*20+ 9] =  21.24830; daa[17*20+10] =  66.53090;
  daa[17*20+11] =  13.75050; daa[17*20+12] =  51.57060; daa[17*20+13] = 152.96400;
  daa[17*20+14] =  13.94050; daa[17*20+15] =  52.37420; daa[17*20+16] =  11.08640;
  daa[18*20+ 0] =  24.07350; daa[18*20+ 1] =  38.15330; daa[18*20+ 2] = 108.60000;
  daa[18*20+ 3] =  32.57110; daa[18*20+ 4] =  54.38330; daa[18*20+ 5] =  22.77100;
  daa[18*20+ 6] =  19.63030; daa[18*20+ 7] =  10.36040; daa[18*20+ 8] = 387.34400;
  daa[18*20+ 9] =  42.01700; daa[18*20+10] =  39.86180; daa[18*20+11] =  13.32640;
  daa[18*20+12] =  42.84370; daa[18*20+13] = 645.42800; daa[18*20+14] =  21.60460;
  daa[18*20+15] =  78.69930; daa[18*20+16] =  29.11480; daa[18*20+17] = 248.53900;
  daa[19*20+ 0] = 200.60100; daa[19*20+ 1] =  25.18490; daa[19*20+ 2] =  19.62460;
  daa[19*20+ 3] =  15.23350; daa[19*20+ 4] = 100.21400; daa[19*20+ 5] =  30.12810;
  daa[19*20+ 6] =  58.87310; daa[19*20+ 7] =  18.72470; daa[19*20+ 8] =  11.83580;
  daa[19*20+ 9] = 782.13000; daa[19*20+10] = 180.03400; daa[19*20+11] =  30.54340;
  daa[19*20+12] = 205.84500; daa[19*20+13] =  64.98920; daa[19*20+14] =  31.48870;
  daa[19*20+15] =  23.27390; daa[19*20+16] = 138.82300; daa[19*20+17] =  36.53690;
  daa[19*20+18] =  31.47300;

  for (i=0; i<naa; i++)  for (j=0; j<i; j++)  daa[j*naa+i] = daa[i*naa+j];

  pi[0] = 0.0866279; pi[1] =  0.043972; pi[2] =  0.0390894; pi[3] =  0.0570451;
  pi[4] =  0.0193078; pi[5] =  0.0367281; pi[6] =  0.0580589; pi[7] =  0.0832518;
  pi[8] =  0.0244313; pi[9] =  0.048466; pi[10] =  0.086209; pi[11] = 0.0620286;
  pi[12] = 0.0195027; pi[13] =  0.0384319; pi[14] =  0.0457631; pi[15] = 0.0695179;
  pi[16] =  0.0610127; pi[17] =  0.0143859; pi[18] =  0.0352742; pi[19] =  0.0708956;


  return 1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


int Init_Qmat_RtREV(phydbl *daa, phydbl *pi)
{
    /*
       This model has been 'translated' from John Huelsenbeck and Fredrik Ronquist
       MrBayes program into PHYML format by Federico Abascal. Many thanks to them.
    */

    /*
    Dimmic M.W., J.S. Rest, D.P. Mindell, and D. Goldstein. 2002. RArtREV:
    An amino acid substitution matrix for inference of retrovirus and
    reverse transcriptase phyLOGeny. Journal of Molecular Evolution
    55: 65-73.
    */

    int i,j,naa;
    naa = 20;

    daa[1*20+0]= 34;         daa[2*20+0]= 51;         daa[2*20+1]= 35;         daa[3*20+0]= 10;
    daa[3*20+1]= 30;         daa[3*20+2]= 384;        daa[4*20+0]= 439;        daa[4*20+1]= 92;
    daa[4*20+2]= 128;        daa[4*20+3]= 1;          daa[5*20+0]= 32;         daa[5*20+1]= 221;
    daa[5*20+2]= 236;        daa[5*20+3]= 78;         daa[5*20+4]= 70;         daa[6*20+0]= 81;
    daa[6*20+1]= 10;         daa[6*20+2]= 79;         daa[6*20+3]= 542;        daa[6*20+4]= 1;
    daa[6*20+5]= 372;        daa[7*20+0]= 135;        daa[7*20+1]= 41;         daa[7*20+2]= 94;
    daa[7*20+3]= 61;         daa[7*20+4]= 48;         daa[7*20+5]= 18;         daa[7*20+6]= 70;
    daa[8*20+0]= 30;         daa[8*20+1]= 90;         daa[8*20+2]= 320;        daa[8*20+3]= 91;
    daa[8*20+4]= 124;        daa[8*20+5]= 387;        daa[8*20+6]= 34;         daa[8*20+7]= 68;
    daa[9*20+0]= 1;          daa[9*20+1]= 24;         daa[9*20+2]= 35;         daa[9*20+3]= 1;
    daa[9*20+4]= 104;        daa[9*20+5]= 33;         daa[9*20+6]= 1;          daa[9*20+7]= 1;
    daa[9*20+8]= 34;         daa[10*20+0]= 45;        daa[10*20+1]= 18;        daa[10*20+2]= 15;
    daa[10*20+3]= 5;         daa[10*20+4]= 110;       daa[10*20+5]= 54;        daa[10*20+6]= 21;
    daa[10*20+7]= 3;         daa[10*20+8]= 51;        daa[10*20+9]= 385;       daa[11*20+0]= 38;
    daa[11*20+1]= 593;       daa[11*20+2]= 123;       daa[11*20+3]= 20;        daa[11*20+4]= 16;
    daa[11*20+5]= 309;       daa[11*20+6]= 141;       daa[11*20+7]= 30;        daa[11*20+8]= 76;
    daa[11*20+9]= 34;        daa[11*20+10]= 23;       daa[12*20+0]= 235;       daa[12*20+1]= 57;
    daa[12*20+2]= 1;         daa[12*20+3]= 1;         daa[12*20+4]= 156;       daa[12*20+5]= 158;
    daa[12*20+6]= 1;         daa[12*20+7]= 37;        daa[12*20+8]= 116;       daa[12*20+9]= 375;
    daa[12*20+10]= 581;      daa[12*20+11]= 134;      daa[13*20+0]= 1;         daa[13*20+1]= 7;
    daa[13*20+2]= 49;        daa[13*20+3]= 1;         daa[13*20+4]= 70;        daa[13*20+5]= 1;
    daa[13*20+6]= 1;         daa[13*20+7]= 7;         daa[13*20+8]= 141;       daa[13*20+9]= 64;
    daa[13*20+10]= 179;      daa[13*20+11]= 14;       daa[13*20+12]= 247;      daa[14*20+0]= 97;
    daa[14*20+1]= 24;        daa[14*20+2]= 33;        daa[14*20+3]= 55;        daa[14*20+4]= 1;
    daa[14*20+5]= 68;        daa[14*20+6]= 52;        daa[14*20+7]= 17;        daa[14*20+8]= 44;
    daa[14*20+9]= 10;        daa[14*20+10]= 22;       daa[14*20+11]= 43;       daa[14*20+12]= 1;
    daa[14*20+13]= 11;       daa[15*20+0]= 460;       daa[15*20+1]= 102;       daa[15*20+2]= 294;
    daa[15*20+3]= 136;       daa[15*20+4]= 75;        daa[15*20+5]= 225;       daa[15*20+6]= 95;
    daa[15*20+7]= 152;       daa[15*20+8]= 183;       daa[15*20+9]= 4;         daa[15*20+10]= 24;
    daa[15*20+11]= 77;       daa[15*20+12]= 1;        daa[15*20+13]= 20;       daa[15*20+14]= 134;
    daa[16*20+0]= 258;       daa[16*20+1]= 64;        daa[16*20+2]= 148;       daa[16*20+3]= 55;
    daa[16*20+4]= 117;       daa[16*20+5]= 146;       daa[16*20+6]= 82;        daa[16*20+7]= 7;
    daa[16*20+8]= 49;        daa[16*20+9]= 72;        daa[16*20+10]= 25;       daa[16*20+11]= 110;
    daa[16*20+12]= 131;      daa[16*20+13]= 69;       daa[16*20+14]= 62;       daa[16*20+15]= 671;
    daa[17*20+0]= 5;         daa[17*20+1]= 13;        daa[17*20+2]= 16;        daa[17*20+3]= 1;
    daa[17*20+4]= 55;        daa[17*20+5]= 10;        daa[17*20+6]= 17;        daa[17*20+7]= 23;
    daa[17*20+8]= 48;        daa[17*20+9]= 39;        daa[17*20+10]= 47;       daa[17*20+11]= 6;
    daa[17*20+12]= 111;      daa[17*20+13]= 182;      daa[17*20+14]= 9;        daa[17*20+15]= 14;
    daa[17*20+16]= 1;        daa[18*20+0]= 55;        daa[18*20+1]= 47;        daa[18*20+2]= 28;
    daa[18*20+3]= 1;         daa[18*20+4]= 131;       daa[18*20+5]= 45;        daa[18*20+6]= 1;
    daa[18*20+7]= 21;        daa[18*20+8]= 307;       daa[18*20+9]= 26;        daa[18*20+10]= 64;
    daa[18*20+11]= 1;        daa[18*20+12]= 74;       daa[18*20+13]= 1017;     daa[18*20+14]= 14;
    daa[18*20+15]= 31;       daa[18*20+16]= 34;       daa[18*20+17]= 176;      daa[19*20+0]= 197;
    daa[19*20+1]= 29;        daa[19*20+2]= 21;        daa[19*20+3]= 6;         daa[19*20+4]= 295;
    daa[19*20+5]= 36;        daa[19*20+6]= 35;        daa[19*20+7]= 3;         daa[19*20+8]= 1;
    daa[19*20+9]= 1048;      daa[19*20+10]= 112;      daa[19*20+11]= 19;       daa[19*20+12]= 236;
    daa[19*20+13]= 92;       daa[19*20+14]= 25;       daa[19*20+15]= 39;       daa[19*20+16]= 196;
    daa[19*20+17]= 26;       daa[19*20+18]= 59;

    for (i=0; i<naa; i++)  for (j=0; j<i; j++)  daa[j*naa+i] = daa[i*naa+j];

    pi[0]= 0.0646;           pi[1]= 0.0453;           pi[2]= 0.0376;           pi[3]= 0.0422;
    pi[4]= 0.0114;           pi[5]= 0.0606;           pi[6]= 0.0607;           pi[7]= 0.0639;
    pi[8]= 0.0273;           pi[9]= 0.0679;           pi[10]= 0.1018;          pi[11]= 0.0751;
    pi[12]= 0.015;           pi[13]= 0.0287;          pi[14]= 0.0681;          pi[15]= 0.0488;
    pi[16]= 0.0622;          pi[17]= 0.0251;          pi[18]= 0.0318;          pi[19]= 0.0619;
    return 1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


int Init_Qmat_CpREV(phydbl *daa, phydbl *pi)
{
    /*
       This model has been 'translated' from John Huelsenbeck and Fredrik Ronquist
       MrBayes program into PHYML format by Federico Abascal. Many thanks to them.
    */

    /*
      Adachi, J., P. Waddell, W. Martin, and M. Hasegawa. 2000. Plastid
      genome phyLOGeny and a model of amino acid substitution for proteins
      encoded by chloroplast DNA. Journal of Molecular Evolution
      50:348-358.
    */

    int i,j,naa;
    naa = 20;

    daa[1*20+0]= 105;        daa[2*20+0]= 227;        daa[2*20+1]= 357;        daa[3*20+0]= 175;
    daa[3*20+1]= 43;         daa[3*20+2]= 4435;       daa[4*20+0]= 669;        daa[4*20+1]= 823;
    daa[4*20+2]= 538;        daa[4*20+3]= 10;         daa[5*20+0]= 157;        daa[5*20+1]= 1745;
    daa[5*20+2]= 768;        daa[5*20+3]= 400;        daa[5*20+4]= 10;         daa[6*20+0]= 499;
    daa[6*20+1]= 152;        daa[6*20+2]= 1055;       daa[6*20+3]= 3691;       daa[6*20+4]= 10;
    daa[6*20+5]= 3122;       daa[7*20+0]= 665;        daa[7*20+1]= 243;        daa[7*20+2]= 653;
    daa[7*20+3]= 431;        daa[7*20+4]= 303;        daa[7*20+5]= 133;        daa[7*20+6]= 379;
    daa[8*20+0]= 66;         daa[8*20+1]= 715;        daa[8*20+2]= 1405;       daa[8*20+3]= 331;
    daa[8*20+4]= 441;        daa[8*20+5]= 1269;       daa[8*20+6]= 162;        daa[8*20+7]= 19;
    daa[9*20+0]= 145;        daa[9*20+1]= 136;        daa[9*20+2]= 168;        daa[9*20+3]= 10;
    daa[9*20+4]= 280;        daa[9*20+5]= 92;         daa[9*20+6]= 148;        daa[9*20+7]= 40;
    daa[9*20+8]= 29;         daa[10*20+0]= 197;       daa[10*20+1]= 203;       daa[10*20+2]= 113;
    daa[10*20+3]= 10;        daa[10*20+4]= 396;       daa[10*20+5]= 286;       daa[10*20+6]= 82;
    daa[10*20+7]= 20;        daa[10*20+8]= 66;        daa[10*20+9]= 1745;      daa[11*20+0]= 236;
    daa[11*20+1]= 4482;      daa[11*20+2]= 2430;      daa[11*20+3]= 412;       daa[11*20+4]= 48;
    daa[11*20+5]= 3313;      daa[11*20+6]= 2629;      daa[11*20+7]= 263;       daa[11*20+8]= 305;
    daa[11*20+9]= 345;       daa[11*20+10]= 218;      daa[12*20+0]= 185;       daa[12*20+1]= 125;
    daa[12*20+2]= 61;        daa[12*20+3]= 47;        daa[12*20+4]= 159;       daa[12*20+5]= 202;
    daa[12*20+6]= 113;       daa[12*20+7]= 21;        daa[12*20+8]= 10;        daa[12*20+9]= 1772;
    daa[12*20+10]= 1351;     daa[12*20+11]= 193;      daa[13*20+0]= 68;        daa[13*20+1]= 53;
    daa[13*20+2]= 97;        daa[13*20+3]= 22;        daa[13*20+4]= 726;       daa[13*20+5]= 10;
    daa[13*20+6]= 145;       daa[13*20+7]= 25;        daa[13*20+8]= 127;       daa[13*20+9]= 454;
    daa[13*20+10]= 1268;     daa[13*20+11]= 72;       daa[13*20+12]= 327;      daa[14*20+0]= 490;
    daa[14*20+1]= 87;        daa[14*20+2]= 173;       daa[14*20+3]= 170;       daa[14*20+4]= 285;
    daa[14*20+5]= 323;       daa[14*20+6]= 185;       daa[14*20+7]= 28;        daa[14*20+8]= 152;
    daa[14*20+9]= 117;       daa[14*20+10]= 219;      daa[14*20+11]= 302;      daa[14*20+12]= 100;
    daa[14*20+13]= 43;       daa[15*20+0]= 2440;      daa[15*20+1]= 385;       daa[15*20+2]= 2085;
    daa[15*20+3]= 590;       daa[15*20+4]= 2331;      daa[15*20+5]= 396;       daa[15*20+6]= 568;
    daa[15*20+7]= 691;       daa[15*20+8]= 303;       daa[15*20+9]= 216;       daa[15*20+10]= 516;
    daa[15*20+11]= 868;      daa[15*20+12]= 93;       daa[15*20+13]= 487;      daa[15*20+14]= 1202;
    daa[16*20+0]= 1340;      daa[16*20+1]= 314;       daa[16*20+2]= 1393;      daa[16*20+3]= 266;
    daa[16*20+4]= 576;       daa[16*20+5]= 241;       daa[16*20+6]= 369;       daa[16*20+7]= 92;
    daa[16*20+8]= 32;        daa[16*20+9]= 1040;      daa[16*20+10]= 156;      daa[16*20+11]= 918;
    daa[16*20+12]= 645;      daa[16*20+13]= 148;      daa[16*20+14]= 260;      daa[16*20+15]= 2151;
    daa[17*20+0]= 14;        daa[17*20+1]= 230;       daa[17*20+2]= 40;        daa[17*20+3]= 18;
    daa[17*20+4]= 435;       daa[17*20+5]= 53;        daa[17*20+6]= 63;        daa[17*20+7]= 82;
    daa[17*20+8]= 69;        daa[17*20+9]= 42;        daa[17*20+10]= 159;      daa[17*20+11]= 10;
    daa[17*20+12]= 86;       daa[17*20+13]= 468;      daa[17*20+14]= 49;       daa[17*20+15]= 73;
    daa[17*20+16]= 29;       daa[18*20+0]= 56;        daa[18*20+1]= 323;       daa[18*20+2]= 754;
    daa[18*20+3]= 281;       daa[18*20+4]= 1466;      daa[18*20+5]= 391;       daa[18*20+6]= 142;
    daa[18*20+7]= 10;        daa[18*20+8]= 1971;      daa[18*20+9]= 89;        daa[18*20+10]= 189;
    daa[18*20+11]= 247;      daa[18*20+12]= 215;      daa[18*20+13]= 2370;     daa[18*20+14]= 97;
    daa[18*20+15]= 522;      daa[18*20+16]= 71;       daa[18*20+17]= 346;      daa[19*20+0]= 968;
    daa[19*20+1]= 92;        daa[19*20+2]= 83;        daa[19*20+3]= 75;        daa[19*20+4]= 592;
    daa[19*20+5]= 54;        daa[19*20+6]= 200;       daa[19*20+7]= 91;        daa[19*20+8]= 25;
    daa[19*20+9]= 4797;      daa[19*20+10]= 865;      daa[19*20+11]= 249;      daa[19*20+12]= 475;
    daa[19*20+13]= 317;      daa[19*20+14]= 122;      daa[19*20+15]= 167;      daa[19*20+16]= 760;
    daa[19*20+17]= 10;       daa[19*20+18]= 119;

    for (i=0; i<naa; i++)  for (j=0; j<i; j++)  daa[j*naa+i] = daa[i*naa+j];

    pi[0]= 0.076;            pi[1]= 0.062;            pi[2]= 0.041;            pi[3]= 0.037;
    pi[4]= 0.009;            pi[5]= 0.038;            pi[6]= 0.049;            pi[7]= 0.084;
    pi[8]= 0.025;            pi[9]= 0.081;            pi[10]= 0.101;           pi[11]= 0.05;
    pi[12]= 0.022;           pi[13]= 0.051;           pi[14]= 0.043;           pi[15]= 0.062;
    pi[16]= 0.054;           pi[17]= 0.018;           pi[18]= 0.031;           pi[19]= 0.066;
    return 1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


int Init_Qmat_VT(phydbl *daa, phydbl *pi)
{
    /*
       This model has been 'translated' from John Huelsenbeck and Fredrik Ronquist
       MrBayes program into PHYML format by Federico Abascal. Many thanks to them.
    */

    /*
      Muller, T., and M. Vingron. 2000. Modeling amino acid replacement.
      Journal of Computational Biology 7:761-776.
    */

    int i,j,naa;
    naa = 20;

    /* daa[1*20+0]= 0.233108;   daa[2*20+0]= 0.199097;   daa[2*20+1]= 0.210797;   daa[3*20+0]= 0.265145;    */
    /* daa[3*20+1]= 0.105191;   daa[3*20+2]= 0.883422;   daa[4*20+0]= 0.227333;   daa[4*20+1]= 0.031726;    */
    /* daa[4*20+2]= 0.027495;   daa[4*20+3]= 0.010313;   daa[5*20+0]= 0.310084;   daa[5*20+1]= 0.493763;    */
    /* daa[5*20+2]= 0.2757;     daa[5*20+3]= 0.205842;   daa[5*20+4]= 0.004315;   daa[6*20+0]= 0.567957;    */
    /* daa[6*20+1]= 0.25524;    daa[6*20+2]= 0.270417;   daa[6*20+3]= 1.599461;   daa[6*20+4]= 0.005321;    */
    /* daa[6*20+5]= 0.960976;   daa[7*20+0]= 0.876213;   daa[7*20+1]= 0.156945;   daa[7*20+2]= 0.362028;    */
    /* daa[7*20+3]= 0.311718;   daa[7*20+4]= 0.050876;   daa[7*20+5]= 0.12866;    daa[7*20+6]= 0.250447;    */
    /* daa[8*20+0]= 0.078692;   daa[8*20+1]= 0.213164;   daa[8*20+2]= 0.290006;   daa[8*20+3]= 0.134252;    */
    /* daa[8*20+4]= 0.016695;   daa[8*20+5]= 0.315521;   daa[8*20+6]= 0.104458;   daa[8*20+7]= 0.058131;    */
    /* daa[9*20+0]= 0.222972;   daa[9*20+1]= 0.08151;    daa[9*20+2]= 0.087225;   daa[9*20+3]= 0.01172;     */
    /* daa[9*20+4]= 0.046398;   daa[9*20+5]= 0.054602;   daa[9*20+6]= 0.046589;   daa[9*20+7]= 0.051089;    */
    /* daa[9*20+8]= 0.020039;   daa[10*20+0]= 0.42463;   daa[10*20+1]= 0.192364;  daa[10*20+2]= 0.069245;   */
    /* daa[10*20+3]= 0.060863;  daa[10*20+4]= 0.091709;  daa[10*20+5]= 0.24353;   daa[10*20+6]= 0.151924;   */
    /* daa[10*20+7]= 0.087056;  daa[10*20+8]= 0.103552;  daa[10*20+9]= 2.08989;   daa[11*20+0]= 0.393245;   */
    /* daa[11*20+1]= 1.755838;  daa[11*20+2]= 0.50306;   daa[11*20+3]= 0.261101;  daa[11*20+4]= 0.004067;   */
    /* daa[11*20+5]= 0.738208;  daa[11*20+6]= 0.88863;   daa[11*20+7]= 0.193243;  daa[11*20+8]= 0.153323;   */
    /* daa[11*20+9]= 0.093181;  daa[11*20+10]= 0.201204; daa[12*20+0]= 0.21155;   daa[12*20+1]= 0.08793;    */
    /* daa[12*20+2]= 0.05742;   daa[12*20+3]= 0.012182;  daa[12*20+4]= 0.02369;   daa[12*20+5]= 0.120801;   */
    /* daa[12*20+6]= 0.058643;  daa[12*20+7]= 0.04656;   daa[12*20+8]= 0.021157;  daa[12*20+9]= 0.493845;   */
    /* daa[12*20+10]= 1.105667; daa[12*20+11]= 0.096474; daa[13*20+0]= 0.116646;  daa[13*20+1]= 0.042569;   */
    /* daa[13*20+2]= 0.039769;  daa[13*20+3]= 0.016577;  daa[13*20+4]= 0.051127;  daa[13*20+5]= 0.026235;   */
    /* daa[13*20+6]= 0.028168;  daa[13*20+7]= 0.050143;  daa[13*20+8]= 0.079807;  daa[13*20+9]= 0.32102;    */
    /* daa[13*20+10]= 0.946499; daa[13*20+11]= 0.038261; daa[13*20+12]= 0.173052; daa[14*20+0]= 0.399143;   */
    /* daa[14*20+1]= 0.12848;   daa[14*20+2]= 0.083956;  daa[14*20+3]= 0.160063;  daa[14*20+4]= 0.011137;   */
    /* daa[14*20+5]= 0.15657;   daa[14*20+6]= 0.205134;  daa[14*20+7]= 0.124492;  daa[14*20+8]= 0.078892;   */
    /* daa[14*20+9]= 0.054797;  daa[14*20+10]= 0.169784; daa[14*20+11]= 0.212302; daa[14*20+12]= 0.010363;  */
    /* daa[14*20+13]= 0.042564; daa[15*20+0]= 1.817198;  daa[15*20+1]= 0.292327;  daa[15*20+2]= 0.847049;   */
    /* daa[15*20+3]= 0.461519;  daa[15*20+4]= 0.17527;   daa[15*20+5]= 0.358017;  daa[15*20+6]= 0.406035;   */
    /* daa[15*20+7]= 0.612843;  daa[15*20+8]= 0.167406;  daa[15*20+9]= 0.081567;  daa[15*20+10]= 0.214977;  */
    /* daa[15*20+11]= 0.400072; daa[15*20+12]= 0.090515; daa[15*20+13]= 0.138119; daa[15*20+14]= 0.430431;  */
    /* daa[16*20+0]= 0.877877;  daa[16*20+1]= 0.204109;  daa[16*20+2]= 0.471268;  daa[16*20+3]= 0.178197;   */
    /* daa[16*20+4]= 0.079511;  daa[16*20+5]= 0.248992;  daa[16*20+6]= 0.321028;  daa[16*20+7]= 0.136266;   */
    /* daa[16*20+8]= 0.101117;  daa[16*20+9]= 0.376588;  daa[16*20+10]= 0.243227; daa[16*20+11]= 0.446646;  */
    /* daa[16*20+12]= 0.184609; daa[16*20+13]= 0.08587;  daa[16*20+14]= 0.207143; daa[16*20+15]= 1.767766;  */
    /* daa[17*20+0]= 0.030309;  daa[17*20+1]= 0.046417;  daa[17*20+2]= 0.010459;  daa[17*20+3]= 0.011393;   */
    /* daa[17*20+4]= 0.007732;  daa[17*20+5]= 0.021248;  daa[17*20+6]= 0.018844;  daa[17*20+7]= 0.02399;    */
    /* daa[17*20+8]= 0.020009;  daa[17*20+9]= 0.034954;  daa[17*20+10]= 0.083439; daa[17*20+11]= 0.023321;  */
    /* daa[17*20+12]= 0.022019; daa[17*20+13]= 0.12805;  daa[17*20+14]= 0.014584; daa[17*20+15]= 0.035933;  */
    /* daa[17*20+16]= 0.020437; daa[18*20+0]= 0.087061;  daa[18*20+1]= 0.09701;   daa[18*20+2]= 0.093268;   */
    /* daa[18*20+3]= 0.051664;  daa[18*20+4]= 0.042823;  daa[18*20+5]= 0.062544;  daa[18*20+6]= 0.0552;     */
    /* daa[18*20+7]= 0.037568;  daa[18*20+8]= 0.286027;  daa[18*20+9]= 0.086237;  daa[18*20+10]= 0.189842;  */
    /* daa[18*20+11]= 0.068689; daa[18*20+12]= 0.073223; daa[18*20+13]= 0.898663; daa[18*20+14]= 0.032043;  */
    /* daa[18*20+15]= 0.121979; daa[18*20+16]= 0.094617; daa[18*20+17]= 0.124746; daa[19*20+0]= 1.230985;   */
    /* daa[19*20+1]= 0.113146;  daa[19*20+2]= 0.049824;  daa[19*20+3]= 0.048769;  daa[19*20+4]= 0.163831;   */
    /* daa[19*20+5]= 0.112027;  daa[19*20+6]= 0.205868;  daa[19*20+7]= 0.082579;  daa[19*20+8]= 0.068575;   */
    /* daa[19*20+9]= 3.65443;   daa[19*20+10]= 1.337571; daa[19*20+11]= 0.144587; daa[19*20+12]= 0.307309;  */
    /* daa[19*20+13]= 0.247329; daa[19*20+14]= 0.129315; daa[19*20+15]= 0.1277;   daa[19*20+16]= 0.740372;  */
    /* daa[19*20+17]= 0.022134; daa[19*20+18]= 0.125733;  */

    /* for (i=0; i<naa; i++)  for (j=0; j<i; j++)  daa[j*naa+i] = daa[i*naa+j]; */

    /* pi[0]= 0.078837;         pi[1]= 0.051238;         pi[2]= 0.042313;         pi[3]= 0.053066;          */
    /* pi[4]= 0.015175;         pi[5]= 0.036713;         pi[6]= 0.061924;         pi[7]= 0.070852;          */
    /* pi[8]= 0.023082;         pi[9]= 0.062056;         pi[10]= 0.096371;        pi[11]= 0.057324;         */
    /* pi[12]= 0.023771;        pi[13]= 0.043296;        pi[14]= 0.043911;        pi[15]= 0.063403;         */
    /* pi[16]= 0.055897;        pi[17]= 0.013272;        pi[18]= 0.034399;        pi[19]= 0.073101;         */
    /* return 1; */

    /*
       This model has been 'translated' from TREE-PUZZLE. Many thanks to  Heiko A. Schmidt,
       Korbinian Strimmer, and Arndt von Haeseler.
    */

    /*
      Muller, T., and M. Vingron. 2000. Modeling amino acid replacement.
      Journal of Computational Biology 7:761-776.
    */

    daa[ 0*20+ 1] = 1.2412691067876198;  daa[ 0*20+ 2] = 1.2184237953498958;
    daa[ 0*20+ 3] = 1.3759368509441177;  daa[ 0*20+ 4] = 2.4731223087544874;
    daa[ 0*20+ 5] = 2.2155167805137470;  daa[ 0*20+ 6] = 2.3379911207495061;
    daa[ 0*20+ 7] = 3.3386555146457697;  daa[ 0*20+ 8] = 0.9615841926910841;
    daa[ 0*20+ 9] = 0.8908203061925510;  daa[ 0*20+10] = 1.0778497408764076;
    daa[ 0*20+11] = 1.4932055816372476;  daa[ 0*20+12] = 1.9006455961717605;
    daa[ 0*20+13] = 0.6883439026872615;  daa[ 0*20+14] = 2.7355620089953550;
    daa[ 0*20+15] = 6.4208961859142883;  daa[ 0*20+16] = 5.2892514169776437;
    daa[ 0*20+17] = 0.5488578478106930;  daa[ 0*20+18] = 0.5411769916657778;
    daa[ 0*20+19] = 4.6501894691803214;

    daa[ 1*20+ 2] = 1.5720770753326880;  daa[ 1*20+ 3] = 0.7550654439001206;
    daa[ 1*20+ 4] = 1.4414262567428417;  daa[ 1*20+ 5] = 5.5120819705248678;
    daa[ 1*20+ 6] = 1.3542404860613146;  daa[ 1*20+ 7] = 1.3121700301622004;
    daa[ 1*20+ 8] = 4.9238668283945266;  daa[ 1*20+ 9] = 0.4323005487925516;
    daa[ 1*20+10] = 0.8386701149158265;  daa[ 1*20+11] = 10.0173308173660018;
    daa[ 1*20+12] = 1.2488638689609959;  daa[ 1*20+13] = 0.4224945197276290;
    daa[ 1*20+14] = 1.3091837782420783;  daa[ 1*20+15] = 1.9202994262316166;
    daa[ 1*20+16] = 1.3363401740560601;  daa[ 1*20+17] = 1.5170142153962840;
    daa[ 1*20+18] = 0.8912614404565405;  daa[ 1*20+19] = 0.7807017855806767;

    daa[ 2*20+ 3] = 7.8584219153689405;  daa[ 2*20+ 4] = 0.9784679122774127;
    daa[ 2*20+ 5] = 3.0143201670924822;  daa[ 2*20+ 6] = 2.0093434778398112;
    daa[ 2*20+ 7] = 2.4117632898861809;  daa[ 2*20+ 8] = 6.1974384977884114;
    daa[ 2*20+ 9] = 0.9179291175331520;  daa[ 2*20+10] = 0.4098311270816011;
    daa[ 2*20+11] = 4.4034547578962568;  daa[ 2*20+12] = 0.9378803706165143;
    daa[ 2*20+13] = 0.5044944273324311;  daa[ 2*20+14] = 0.7103720531974738;
    daa[ 2*20+15] = 6.1234512396801764;  daa[ 2*20+16] = 3.8852506105922231;
    daa[ 2*20+17] = 0.1808525752605976;  daa[ 2*20+18] = 1.0894926581511342;
    daa[ 2*20+19] = 0.4586061981719967;

    daa[ 3*20+ 4] = 0.2272488448121475;  daa[ 3*20+ 5] = 1.6562495638176040;
    daa[ 3*20+ 6] = 9.6883451875685065;  daa[ 3*20+ 7] = 1.9142079025990228;
    daa[ 3*20+ 8] = 2.1459640610133781;  daa[ 3*20+ 9] = 0.2161660372725585;
    daa[ 3*20+10] = 0.3574207468998517;  daa[ 3*20+11] = 1.4521790561663968;
    daa[ 3*20+12] = 0.4075239926000898;  daa[ 3*20+13] = 0.1675129724559251;
    daa[ 3*20+14] = 1.0714605979577547;  daa[ 3*20+15] = 2.2161944596741829;
    daa[ 3*20+16] = 1.5066839872944762;  daa[ 3*20+17] = 0.2496584188151770;
    daa[ 3*20+18] = 0.7447620891784513;  daa[ 3*20+19] = 0.4594535241660911;

    daa[ 4*20+ 5] = 0.4587469126746136;  daa[ 4*20+ 6] = 0.4519167943192672;
    daa[ 4*20+ 7] = 1.1034605684472507;  daa[ 4*20+ 8] = 1.5196756759380692;
    daa[ 4*20+ 9] = 0.9126668032539315;  daa[ 4*20+10] = 1.4081315998413697;
    daa[ 4*20+11] = 0.3371091785647479;  daa[ 4*20+12] = 1.2213054800811556;
    daa[ 4*20+13] = 1.6953951980808002;  daa[ 4*20+14] = 0.4326227078645523;
    daa[ 4*20+15] = 3.6366815408744255;  daa[ 4*20+16] = 1.7557065205837685;
    daa[ 4*20+17] = 1.6275179891253113;  daa[ 4*20+18] = 2.1579775140421025;
    daa[ 4*20+19] = 2.2627456996290891;

    daa[ 5*20+ 6] = 6.8124601839937675;  daa[ 5*20+ 7] = 0.8776110594765502;
    daa[ 5*20+ 8] = 7.9943228564946525;  daa[ 5*20+ 9] = 0.4882733432879921;
    daa[ 5*20+10] = 1.3318097154194044;  daa[ 5*20+11] = 6.0519085243118811;
    daa[ 5*20+12] = 1.9106190827629084;  daa[ 5*20+13] = 0.3573432522499545;
    daa[ 5*20+14] = 2.3019177728300728;  daa[ 5*20+15] = 2.3193703643237220;
    daa[ 5*20+16] = 2.1576510103471440;  daa[ 5*20+17] = 0.8959082681546182;
    daa[ 5*20+18] = 0.9183596801412757;  daa[ 5*20+19] = 0.6366932501396869;

    daa[ 6*20+ 7] = 1.3860121390169038;  daa[ 6*20+ 8] = 1.6360079688522375;
    daa[ 6*20+ 9] = 0.4035497929633328;  daa[ 6*20+10] = 0.5610717242294755;
    daa[ 6*20+11] = 4.3290086529582830;  daa[ 6*20+12] = 0.7471936218068498;
    daa[ 6*20+13] = 0.2317194387691585;  daa[ 6*20+14] = 1.5132807416252063;
    daa[ 6*20+15] = 1.8273535587773553;  daa[ 6*20+16] = 1.5839981708584689;
    daa[ 6*20+17] = 0.4198391148111098;  daa[ 6*20+18] = 0.5818111331782764;
    daa[ 6*20+19] = 0.8940572875547330;

    daa[ 7*20+ 8] = 0.8561248973045037;  daa[ 7*20+ 9] = 0.2888075033037488;
    daa[ 7*20+10] = 0.3578662395745526;  daa[ 7*20+11] = 0.8945563662345198;
    daa[ 7*20+12] = 0.5954812791740037;  daa[ 7*20+13] = 0.3693722640980460;
    daa[ 7*20+14] = 0.7744933618134962;  daa[ 7*20+15] = 3.0637776193717610;
    daa[ 7*20+16] = 0.7147489676267383;  daa[ 7*20+17] = 0.9349753595598769;
    daa[ 7*20+18] = 0.3374467649724478;  daa[ 7*20+19] = 0.6193321034173915;

    daa[ 8*20+ 9] = 0.5787937115407940;  daa[ 8*20+10] = 1.0765007949562073;
    daa[ 8*20+11] = 1.8085136096039203;  daa[ 8*20+12] = 1.3808291710019667;
    daa[ 8*20+13] = 1.3629765501081097;  daa[ 8*20+14] = 1.8370555852070649;
    daa[ 8*20+15] = 1.9699895187387506;  daa[ 8*20+16] = 1.6136654573285647;
    daa[ 8*20+17] = 0.6301954684360302;  daa[ 8*20+18] = 7.7587442309146040;
    daa[ 8*20+19] = 0.5333220944030346;

    daa[ 9*20+10] = 6.0019110258426362;  daa[ 9*20+11] = 0.6244297525127139;
    daa[ 9*20+12] = 6.7597899772045418;  daa[ 9*20+13] = 2.2864286949316077;
    daa[ 9*20+14] = 0.4811402387911145;  daa[ 9*20+15] = 0.6047491507504744;
    daa[ 9*20+16] = 2.6344778384442731;  daa[ 9*20+17] = 0.5604648274060783;
    daa[ 9*20+18] = 0.8626796044156272;  daa[ 9*20+19] = 14.8729334615190609;

    daa[10*20+11] = 0.5642322882556321;  daa[10*20+12] = 8.0327792947421148;
    daa[10*20+13] = 4.3611548063555778;  daa[10*20+14] = 1.0084320519837335;
    daa[10*20+15] = 0.8953754669269811;  daa[10*20+16] = 1.0192004372506540;
    daa[10*20+17] = 1.5183114434679339;  daa[10*20+18] = 1.2452243224541324;
    daa[10*20+19] = 3.5458093276667237;

    daa[11*20+12] = 1.7129670976916258;  daa[11*20+13] = 0.3910559903834828;
    daa[11*20+14] = 1.3918935593582853;  daa[11*20+15] = 1.9776630140912268;
    daa[11*20+16] = 2.5513781312660280;  daa[11*20+17] = 0.5851920879490173;
    daa[11*20+18] = 0.7835447533710449;  daa[11*20+19] = 0.7801080335991272;

    daa[12*20+13] = 2.3201373546296349;  daa[12*20+14] = 0.4953193808676289;
    daa[12*20+15] = 1.0657482318076852;  daa[12*20+16] = 3.3628488360462363;
    daa[12*20+17] = 1.4680478689711018;  daa[12*20+18] = 1.0899165770956820;
    daa[12*20+19] = 4.0584577156753401;

    daa[13*20+14] = 0.3746821107962129;  daa[13*20+15] = 1.1079144700606407;
    daa[13*20+16] = 0.6882725908872254;  daa[13*20+17] = 3.3448437239772266;
    daa[13*20+18] = 10.3848523331334590;  daa[13*20+19] = 1.7039730522675411;

    daa[14*20+15] = 3.5465914843628927;  daa[14*20+16] = 1.9485376673137556;
    daa[14*20+17] = 0.4326058001438786;  daa[14*20+18] = 0.4819109019647465;
    daa[14*20+19] = 0.5985498912985666;

    daa[15*20+16] = 8.8479984061248178;  daa[15*20+17] = 0.6791126595939816;
    daa[15*20+18] = 0.9547229305958682;  daa[15*20+19] = 0.9305232113028208;

    daa[16*20+17] = 0.4514203099376473;  daa[16*20+18] = 0.8564314184691215;
    daa[16*20+19] = 3.4242218450865543;

    daa[17*20+18] = 4.5377235790405388;  daa[17*20+19] = 0.5658969249032649;

    daa[18*20+19] = 1.0000000000000000;

    for (i=0; i<naa; i++)  for (j=0; j<i; j++)  daa[i*naa+j] = daa[j*naa+i];

    pi[ 0]=0.0770764620135024 ; pi[ 1]=0.0500819370772208 ;
    pi[ 2]=0.0462377395993731 ; pi[ 3]=0.0537929860758246 ;
    pi[ 4]=0.0144533387583345 ; pi[ 5]=0.0408923608974345 ;
    pi[ 6]=0.0633579339160905 ; pi[ 7]=0.0655672355884439 ;
    pi[ 8]=0.0218802687005936 ; pi[ 9]=0.0591969699027449 ;
    pi[10]=0.0976461276528445 ; pi[11]=0.0592079410822730 ;
    pi[12]=0.0220695876653368 ; pi[13]=0.0413508521834260 ;
    pi[14]=0.0476871596856874 ; pi[15]=0.0707295165111524 ;
    pi[16]=0.0567759161524817 ; pi[17]=0.0127019797647213 ;
    pi[18]=0.0323746050281867 ; pi[19]=0.0669190817443274 ;

    return 1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


int Init_Qmat_Blosum62(phydbl *daa, phydbl *pi)
{

    /*
       This model has been 'translated' from John Huelsenbeck and Fredrik Ronquist
       MrBayes program into PHYML format by Federico Abascal. Many thanks to them.
    */

    /*
      Henikoff, S., and J. G. Henikoff. 1992. Amino acid substitution
      matrices from protein blocks. Proc. Natl. Acad. Sci., U.S.A.
      89:10915-10919.
    */

    int i,j,naa;
    naa = 20;

    daa[1*20+0]= 0.735790389698;  daa[2*20+0]= 0.485391055466;  daa[2*20+1]= 1.297446705134;  daa[3*20+0]= 0.543161820899;
    daa[3*20+1]= 0.500964408555;  daa[3*20+2]= 3.180100048216;  daa[4*20+0]= 1.45999531047;   daa[4*20+1]= 0.227826574209;
    daa[4*20+2]= 0.397358949897;  daa[4*20+3]= 0.240836614802;  daa[5*20+0]= 1.199705704602;  daa[5*20+1]= 3.020833610064;
    daa[5*20+2]= 1.839216146992;  daa[5*20+3]= 1.190945703396;  daa[5*20+4]= 0.32980150463;   daa[6*20+0]= 1.1709490428;
    daa[6*20+1]= 1.36057419042;   daa[6*20+2]= 1.24048850864;   daa[6*20+3]= 3.761625208368;  daa[6*20+4]= 0.140748891814;
    daa[6*20+5]= 5.528919177928;  daa[7*20+0]= 1.95588357496;   daa[7*20+1]= 0.418763308518;  daa[7*20+2]= 1.355872344485;
    daa[7*20+3]= 0.798473248968;  daa[7*20+4]= 0.418203192284;  daa[7*20+5]= 0.609846305383;  daa[7*20+6]= 0.423579992176;
    daa[8*20+0]= 0.716241444998;  daa[8*20+1]= 1.456141166336;  daa[8*20+2]= 2.414501434208;  daa[8*20+3]= 0.778142664022;
    daa[8*20+4]= 0.354058109831;  daa[8*20+5]= 2.43534113114;   daa[8*20+6]= 1.626891056982;  daa[8*20+7]= 0.539859124954;
    daa[9*20+0]= 0.605899003687;  daa[9*20+1]= 0.232036445142;  daa[9*20+2]= 0.283017326278;  daa[9*20+3]= 0.418555732462;
    daa[9*20+4]= 0.774894022794;  daa[9*20+5]= 0.236202451204;  daa[9*20+6]= 0.186848046932;  daa[9*20+7]= 0.189296292376;
    daa[9*20+8]= 0.252718447885;  daa[10*20+0]= 0.800016530518; daa[10*20+1]= 0.622711669692; daa[10*20+2]= 0.211888159615;
    daa[10*20+3]= 0.218131577594; daa[10*20+4]= 0.831842640142; daa[10*20+5]= 0.580737093181; daa[10*20+6]= 0.372625175087;
    daa[10*20+7]= 0.217721159236; daa[10*20+8]= 0.348072209797; daa[10*20+9]= 3.890963773304; daa[11*20+0]= 1.295201266783;
    daa[11*20+1]= 5.411115141489; daa[11*20+2]= 1.593137043457; daa[11*20+3]= 1.032447924952; daa[11*20+4]= 0.285078800906;
    daa[11*20+5]= 3.945277674515; daa[11*20+6]= 2.802427151679; daa[11*20+7]= 0.752042440303; daa[11*20+8]= 1.022507035889;
    daa[11*20+9]= 0.406193586642; daa[11*20+10]= 0.445570274261;daa[12*20+0]= 1.253758266664; daa[12*20+1]= 0.983692987457;
    daa[12*20+2]= 0.648441278787; daa[12*20+3]= 0.222621897958; daa[12*20+4]= 0.76768882348;  daa[12*20+5]= 2.494896077113;
    daa[12*20+6]= 0.55541539747;  daa[12*20+7]= 0.459436173579; daa[12*20+8]= 0.984311525359; daa[12*20+9]= 3.364797763104;
    daa[12*20+10]= 6.030559379572;daa[12*20+11]= 1.073061184332;daa[13*20+0]= 0.492964679748; daa[13*20+1]= 0.371644693209;
    daa[13*20+2]= 0.354861249223; daa[13*20+3]= 0.281730694207; daa[13*20+4]= 0.441337471187; daa[13*20+5]= 0.14435695975;
    daa[13*20+6]= 0.291409084165; daa[13*20+7]= 0.368166464453; daa[13*20+8]= 0.714533703928; daa[13*20+9]= 1.517359325954;
    daa[13*20+10]= 2.064839703237;daa[13*20+11]= 0.266924750511;daa[13*20+12]= 1.77385516883; daa[14*20+0]= 1.173275900924;
    daa[14*20+1]= 0.448133661718; daa[14*20+2]= 0.494887043702; daa[14*20+3]= 0.730628272998; daa[14*20+4]= 0.356008498769;
    daa[14*20+5]= 0.858570575674; daa[14*20+6]= 0.926563934846; daa[14*20+7]= 0.504086599527; daa[14*20+8]= 0.527007339151;
    daa[14*20+9]= 0.388355409206; daa[14*20+10]= 0.374555687471;daa[14*20+11]= 1.047383450722;daa[14*20+12]= 0.454123625103;
    daa[14*20+13]= 0.233597909629;daa[15*20+0]= 4.325092687057; daa[15*20+1]= 1.12278310421;  daa[15*20+2]= 2.904101656456;
    daa[15*20+3]= 1.582754142065; daa[15*20+4]= 1.197188415094; daa[15*20+5]= 1.934870924596; daa[15*20+6]= 1.769893238937;
    daa[15*20+7]= 1.509326253224; daa[15*20+8]= 1.11702976291;  daa[15*20+9]= 0.35754441246;  daa[15*20+10]= 0.352969184527;
    daa[15*20+11]= 1.752165917819;daa[15*20+12]= 0.918723415746;daa[15*20+13]= 0.540027644824;daa[15*20+14]= 1.169129577716;
    daa[16*20+0]= 1.729178019485; daa[16*20+1]= 0.914665954563; daa[16*20+2]= 1.898173634533; daa[16*20+3]= 0.934187509431;
    daa[16*20+4]= 1.119831358516; daa[16*20+5]= 1.277480294596; daa[16*20+6]= 1.071097236007; daa[16*20+7]= 0.641436011405;
    daa[16*20+8]= 0.585407090225; daa[16*20+9]= 1.17909119726;  daa[16*20+10]= 0.915259857694;daa[16*20+11]= 1.303875200799;
    daa[16*20+12]= 1.488548053722;daa[16*20+13]= 0.488206118793;daa[16*20+14]= 1.005451683149;daa[16*20+15]= 5.15155629227;
    daa[17*20+0]= 0.465839367725; daa[17*20+1]= 0.426382310122; daa[17*20+2]= 0.191482046247; daa[17*20+3]= 0.145345046279;
    daa[17*20+4]= 0.527664418872; daa[17*20+5]= 0.758653808642; daa[17*20+6]= 0.407635648938; daa[17*20+7]= 0.508358924638;
    daa[17*20+8]= 0.30124860078;  daa[17*20+9]= 0.34198578754;  daa[17*20+10]= 0.6914746346;  daa[17*20+11]= 0.332243040634;
    daa[17*20+12]= 0.888101098152;daa[17*20+13]= 2.074324893497;daa[17*20+14]= 0.252214830027;daa[17*20+15]= 0.387925622098;
    daa[17*20+16]= 0.513128126891;daa[18*20+0]= 0.718206697586; daa[18*20+1]= 0.720517441216; daa[18*20+2]= 0.538222519037;
    daa[18*20+3]= 0.261422208965; daa[18*20+4]= 0.470237733696; daa[18*20+5]= 0.95898974285;  daa[18*20+6]= 0.596719300346;
    daa[18*20+7]= 0.308055737035; daa[18*20+8]= 4.218953969389; daa[18*20+9]= 0.674617093228; daa[18*20+10]= 0.811245856323;
    daa[18*20+11]= 0.7179934869;  daa[18*20+12]= 0.951682162246;daa[18*20+13]= 6.747260430801;daa[18*20+14]= 0.369405319355;
    daa[18*20+15]= 0.796751520761;daa[18*20+16]= 0.801010243199;daa[18*20+17]= 4.054419006558;daa[19*20+0]= 2.187774522005;
    daa[19*20+1]= 0.438388343772; daa[19*20+2]= 0.312858797993; daa[19*20+3]= 0.258129289418; daa[19*20+4]= 1.116352478606;
    daa[19*20+5]= 0.530785790125; daa[19*20+6]= 0.524253846338; daa[19*20+7]= 0.25334079019;  daa[19*20+8]= 0.20155597175;
    daa[19*20+9]= 8.311839405458; daa[19*20+10]= 2.231405688913;daa[19*20+11]= 0.498138475304;daa[19*20+12]= 2.575850755315;
    daa[19*20+13]= 0.838119610178;daa[19*20+14]= 0.496908410676;daa[19*20+15]= 0.561925457442;daa[19*20+16]= 2.253074051176;
    daa[19*20+17]= 0.266508731426;daa[19*20+18]= 1;

    for (i=0; i<naa; i++)  for (j=0; j<i; j++)  daa[j*naa+i] = daa[i*naa+j];

    pi[0]= 0.074;                 pi[1]= 0.052;                 pi[2]= 0.045;                 pi[3]= 0.054;
    pi[4]= 0.025;                 pi[5]= 0.034;                 pi[6]= 0.054;                 pi[7]= 0.074;
    pi[8]= 0.026;                 pi[9]= 0.068;                 pi[10]= 0.099;                pi[11]= 0.058;
    pi[12]= 0.025;                pi[13]= 0.047;                pi[14]= 0.039;                pi[15]= 0.057;
    pi[16]= 0.051;                pi[17]= 0.013;                pi[18]= 0.032;                pi[19]= 0.073;

    return 1;
  }

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


int Init_Qmat_MtMam(phydbl *daa, phydbl *pi)
{
    /*
       This model has been 'translated' from Ziheng Yang's PAML program
       into PHYML format by Federico Abascal. Many thanks to them.
    */

    /*
      Cao, Y. et al. 1998 Conflict amongst individual mitochondrial
      proteins in resolving the phyLOGeny of eutherian orders. Journal
      of Molecular Evolution 15:1600-1611.
    */

    int i,j,naa;
    naa = 20;

    daa[1*20+0]= 32;              daa[2*20+0]= 2;    daa[2*20+1]= 4;               daa[3*20+0]= 11;
    daa[3*20+1]= 0;               daa[3*20+2]= 864;  daa[4*20+0]= 0;               daa[4*20+1]= 186;
    daa[4*20+2]= 0;               daa[4*20+3]= 0;    daa[5*20+0]= 0;               daa[5*20+1]= 246;
    daa[5*20+2]= 8;               daa[5*20+3]= 49;   daa[5*20+4]= 0;               daa[6*20+0]= 0;
    daa[6*20+1]= 0;               daa[6*20+2]= 0;    daa[6*20+3]= 569;             daa[6*20+4]= 0;
    daa[6*20+5]= 274;             daa[7*20+0]= 78;   daa[7*20+1]= 18;              daa[7*20+2]= 47;
    daa[7*20+3]= 79;              daa[7*20+4]= 0;    daa[7*20+5]= 0;               daa[7*20+6]= 22;
    daa[8*20+0]= 8;               daa[8*20+1]= 232;  daa[8*20+2]= 458;             daa[8*20+3]= 11;
    daa[8*20+4]= 305;             daa[8*20+5]= 550;  daa[8*20+6]= 22;              daa[8*20+7]= 0;
    daa[9*20+0]= 75;              daa[9*20+1]= 0;    daa[9*20+2]= 19;              daa[9*20+3]= 0;
    daa[9*20+4]= 41;              daa[9*20+5]= 0;    daa[9*20+6]= 0;               daa[9*20+7]= 0;
    daa[9*20+8]= 0;               daa[10*20+0]= 21;  daa[10*20+1]= 6;              daa[10*20+2]= 0;
    daa[10*20+3]= 0;              daa[10*20+4]= 27;  daa[10*20+5]= 20;             daa[10*20+6]= 0;
    daa[10*20+7]= 0;              daa[10*20+8]= 26;  daa[10*20+9]= 232;            daa[11*20+0]= 0;
    daa[11*20+1]= 50;             daa[11*20+2]= 408; daa[11*20+3]= 0;              daa[11*20+4]= 0;
    daa[11*20+5]= 242;            daa[11*20+6]= 215; daa[11*20+7]= 0;              daa[11*20+8]= 0;
    daa[11*20+9]= 6;              daa[11*20+10]= 4;  daa[12*20+0]= 76;             daa[12*20+1]= 0;
    daa[12*20+2]= 21;             daa[12*20+3]= 0;   daa[12*20+4]= 0;              daa[12*20+5]= 22;
    daa[12*20+6]= 0;              daa[12*20+7]= 0;   daa[12*20+8]= 0;              daa[12*20+9]= 378;
    daa[12*20+10]= 609;           daa[12*20+11]= 59; daa[13*20+0]= 0;              daa[13*20+1]= 0;
    daa[13*20+2]= 6;              daa[13*20+3]= 5;   daa[13*20+4]= 7;              daa[13*20+5]= 0;
    daa[13*20+6]= 0;              daa[13*20+7]= 0;   daa[13*20+8]= 0;              daa[13*20+9]= 57;
    daa[13*20+10]= 246;           daa[13*20+11]= 0;  daa[13*20+12]= 11;            daa[14*20+0]= 53;
    daa[14*20+1]= 9;              daa[14*20+2]= 33;  daa[14*20+3]= 2;              daa[14*20+4]= 0;
    daa[14*20+5]= 51;             daa[14*20+6]= 0;   daa[14*20+7]= 0;              daa[14*20+8]= 53;
    daa[14*20+9]= 5;              daa[14*20+10]= 43; daa[14*20+11]= 18;            daa[14*20+12]= 0;
    daa[14*20+13]= 17;            daa[15*20+0]= 342; daa[15*20+1]= 3;              daa[15*20+2]= 446;
    daa[15*20+3]= 16;             daa[15*20+4]= 347; daa[15*20+5]= 30;             daa[15*20+6]= 21;
    daa[15*20+7]= 112;            daa[15*20+8]= 20;  daa[15*20+9]= 0;              daa[15*20+10]= 74;
    daa[15*20+11]= 65;            daa[15*20+12]= 47; daa[15*20+13]= 90;            daa[15*20+14]= 202;
    daa[16*20+0]= 681;            daa[16*20+1]= 0;   daa[16*20+2]= 110;            daa[16*20+3]= 0;
    daa[16*20+4]= 114;            daa[16*20+5]= 0;   daa[16*20+6]= 4;              daa[16*20+7]= 0;
    daa[16*20+8]= 1;              daa[16*20+9]= 360; daa[16*20+10]= 34;            daa[16*20+11]= 50;
    daa[16*20+12]= 691;           daa[16*20+13]= 8;  daa[16*20+14]= 78;            daa[16*20+15]= 614;
    daa[17*20+0]= 5;              daa[17*20+1]= 16;  daa[17*20+2]= 6;              daa[17*20+3]= 0;
    daa[17*20+4]= 65;             daa[17*20+5]= 0;   daa[17*20+6]= 0;              daa[17*20+7]= 0;
    daa[17*20+8]= 0;              daa[17*20+9]= 0;   daa[17*20+10]= 12;            daa[17*20+11]= 0;
    daa[17*20+12]= 13;            daa[17*20+13]= 0;  daa[17*20+14]= 7;             daa[17*20+15]= 17;
    daa[17*20+16]= 0;             daa[18*20+0]= 0;   daa[18*20+1]= 0;              daa[18*20+2]= 156;
    daa[18*20+3]= 0;              daa[18*20+4]= 530; daa[18*20+5]= 54;             daa[18*20+6]= 0;
    daa[18*20+7]= 1;              daa[18*20+8]= 1525;daa[18*20+9]= 16;             daa[18*20+10]= 25;
    daa[18*20+11]= 67;            daa[18*20+12]= 0;  daa[18*20+13]= 682;           daa[18*20+14]= 8;
    daa[18*20+15]= 107;           daa[18*20+16]= 0;  daa[18*20+17]= 14;            daa[19*20+0]= 398;
    daa[19*20+1]= 0;              daa[19*20+2]= 0;   daa[19*20+3]= 10;             daa[19*20+4]= 0;
    daa[19*20+5]= 33;             daa[19*20+6]= 20;  daa[19*20+7]= 5;              daa[19*20+8]= 0;
    daa[19*20+9]= 2220;           daa[19*20+10]= 100;daa[19*20+11]= 0;             daa[19*20+12]= 832;
    daa[19*20+13]= 6;             daa[19*20+14]= 0;  daa[19*20+15]= 0;             daa[19*20+16]= 237;
    daa[19*20+17]= 0;             daa[19*20+18]= 0;

    for (i=0; i<naa; i++) for (j=0; j<i; j++) daa[j*naa+i] = daa[i*naa+j];

    pi[0]= 0.0692;  pi[1]=  0.0184;  pi[2]= 0.04;    pi[3]= 0.0186;
    pi[4]= 0.0065;  pi[5]=  0.0238;  pi[6]= 0.0236;  pi[7]= 0.0557;
    pi[8]= 0.0277;  pi[9]=  0.0905;  pi[10]=0.1675;  pi[11]= 0.0221;
    pi[12]=0.0561;  pi[13]= 0.0611;  pi[14]=0.0536;  pi[15]= 0.0725;
    pi[16]=0.087;   pi[17]= 0.0293;  pi[18]=0.034;   pi[19]= 0.0428;

    return 1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


int Init_Qmat_AB(phydbl *daa, phydbl *pi)
{
    /*
       This model has been 'translated' from Alexander Mirsky,Linda Kazandjian and Maria Anisimova into PHYML format by J-E Longueville.
       Many thanks to them.
    */

    /*
    Antibody-Specific Model of Amino Acid Substitution for Immunological Inferences from Alignments of Antibody Sequences
    Alexander Mirsky, Linda Kazandjian and Maria Anisimova
    Mol Biol Evol (2015) 32 (3): 806-819. doi: 10.1093/molbev/msu340 First published online: December 21, 2014
    */

    int i,j,naa;
    naa = 20;

    daa[1*20+0]= 1.784266E-01;         daa[2*20+0]= 9.291290E-02;       daa[2*20+1]= 7.829130E-01;       daa[3*20+0]= 1.241095E+00;
    daa[3*20+1]= 5.795374E-02;         daa[3*20+2]= 7.185182E+00;       daa[4*20+0]= 8.929181E-03;       daa[4*20+1]= 1.821885E-01;
    daa[4*20+2]= 1.374268E-06;         daa[4*20+3]= 2.340019E-02;       daa[5*20+0]= 1.992269E-01;       daa[5*20+1]= 1.923901E+00;
    daa[5*20+2]= 8.705989E-02;         daa[5*20+3]= 1.843856E-01;       daa[5*20+4]= 1.046446E-08;       daa[6*20+0]= 9.521821E-01;
    daa[6*20+1]= 6.273863E-02;         daa[6*20+2]= 5.038373E-01;       daa[6*20+3]= 7.426619E+00;       daa[6*20+4]= 7.519215E-11;
    daa[6*20+5]= 3.691671E+00;         daa[7*20+0]= 1.851951E+00;       daa[7*20+1]= 1.089400E+00;       daa[7*20+2]= 4.868901E-01;
    daa[7*20+3]= 2.112400E+00;         daa[7*20+4]= 5.891123E-02;       daa[7*20+5]= 5.516340E-02;       daa[7*20+6]= 1.389370E+00;
    daa[8*20+0]= 5.241316E+00;         daa[8*20+1]= 1.049550E+01;       daa[8*20+2]= 1.405444E+01;       daa[8*20+3]= 1.126995E+01;
    daa[8*20+4]= 3.963388E+00;         daa[8*20+5]= 8.908434E+00;       daa[8*20+6]= 7.298080E+00;       daa[8*20+7]= 9.139518E+00;
    daa[9*20+0]= 1.140412E-01;          daa[9*20+1]= 3.245175E-01;         daa[9*20+2]= 1.762721E+00;         daa[9*20+3]= 3.916999E-02;
    daa[9*20+4]= 6.594967E-04;        daa[9*20+5]= 6.712736E-06;         daa[9*20+6]= 1.029959E-04;          daa[9*20+7]= 3.560482E-02;
    daa[9*20+8]= 4.706586E+00;         daa[10*20+0]= 6.969101E-02;        daa[10*20+1]= 3.932002E-01;        daa[10*20+2]= 2.769442E-02;
    daa[10*20+3]= 3.020502E-02;         daa[10*20+4]= 6.079219E-03;       daa[10*20+5]= 6.802781E-01;        daa[10*20+6]= 1.283121E-03;
    daa[10*20+7]= 2.157936E-02;         daa[10*20+8]= 5.879103E+00;        daa[10*20+9]= 1.601123E+00;       daa[11*20+0]= 7.388355E-02;
    daa[11*20+1]= 7.549240E+00;       daa[11*20+2]= 6.190318E+00;       daa[11*20+3]= 6.622772E-02;        daa[11*20+4]= 3.722878E-16;
    daa[11*20+5]= 3.030805E+00;       daa[11*20+6]= 3.608816E+00;       daa[11*20+7]= 5.504400E-02;        daa[11*20+8]= 1.455741E+00;
    daa[11*20+9]= 5.059793E-01;        daa[11*20+10]= 2.158451E-02;       daa[12*20+0]= 6.299271E-02;       daa[12*20+1]= 3.362326E-01;
    daa[12*20+2]= 3.972173E-02;         daa[12*20+3]= 3.357577E-02;         daa[12*20+4]= 7.213178E-03;       daa[12*20+5]= 1.233336E-03;
    daa[12*20+6]= 7.659566E-02;         daa[12*20+7]= 2.187264E-02;        daa[12*20+8]= 2.298295E+00;       daa[12*20+9]= 1.096748E+01;
    daa[12*20+10]= 5.647985E+00;      daa[12*20+11]= 1.238634E+00;      daa[13*20+0]= 1.130146E-01;         daa[13*20+1]= 8.208677E-02;
    daa[13*20+2]= 1.955446E-01;        daa[13*20+3]= 1.031734E-01;         daa[13*20+4]= 1.993818E-01;        daa[13*20+5]= 1.496610E-03;
    daa[13*20+6]= 5.288625E-02;         daa[13*20+7]= 1.984772E-01;         daa[13*20+8]= 5.642309E+00;       daa[13*20+9]= 2.714705E+00;
    daa[13*20+10]= 3.390618E+00;      daa[13*20+11]= 4.649035E-03;       daa[13*20+12]= 3.947940E+00;      daa[14*20+0]= 1.800713E+00;
    daa[14*20+1]= 3.498713E-01;        daa[14*20+2]= 7.342554E-03;        daa[14*20+3]= 1.509482E-01;        daa[14*20+4]= 4.878395E-03;
    daa[14*20+5]= 7.426909E-01;        daa[14*20+6]= 2.889815E-02;        daa[14*20+7]= 7.915056E-02;        daa[14*20+8]= 1.049496E+01;
    daa[14*20+9]= 5.016568E-02;        daa[14*20+10]= 1.149931E+00;       daa[14*20+11]= 9.948994E-03;       daa[14*20+12]= 7.417279E-02;
    daa[14*20+13]= 3.556198E-01;       daa[15*20+0]= 9.988358E-01;       daa[15*20+1]= 1.926435E+00;       daa[15*20+2]= 7.348346E+00;
    daa[15*20+3]= 5.822988E-01;       daa[15*20+4]= 2.639482E-01;        daa[15*20+5]= 5.906405E-04;       daa[15*20+6]= 6.776709E-02;
    daa[15*20+7]= 9.984215E-01;       daa[15*20+8]= 5.439116E+00;       daa[15*20+9]= 6.007607E-01;         daa[15*20+10]= 1.580539E-01;
    daa[15*20+11]= 8.688405E-02;       daa[15*20+12]= 1.861354E-02;        daa[15*20+13]= 9.813064E-01;       daa[15*20+14]= 1.284651E+00;
    daa[16*20+0]= 2.912317E+00;       daa[16*20+1]= 2.147175E+00;        daa[16*20+2]= 1.135258E+00;       daa[16*20+3]= 1.516881E-01;
    daa[16*20+4]= 3.225214E-06;       daa[16*20+5]= 1.202094E-01;       daa[16*20+6]= 6.016624E-02;        daa[16*20+7]= 7.862767E-02;
    daa[16*20+8]= 3.443285E+00;        daa[16*20+9]= 3.087152E+00;        daa[16*20+10]= 5.702792E-01;       daa[16*20+11]= 1.039298E+00;
    daa[16*20+12]= 1.415612E+00;      daa[16*20+13]= 3.674486E-02;       daa[16*20+14]= 9.057112E-01;       daa[16*20+15]= 3.058575E+00;
    daa[17*20+0]= 7.939549E-02;         daa[17*20+1]= 5.724286E-01;        daa[17*20+2]= 7.310937E-04;        daa[17*20+3]= 1.423897E-02;
    daa[17*20+4]= 4.440833E-01;        daa[17*20+5]= 4.332983E-05;        daa[17*20+6]= 2.252612E-02;        daa[17*20+7]= 1.386853E-01;
    daa[17*20+8]= 7.013890E+00;        daa[17*20+9]= 6.318748E-02;        daa[17*20+10]= 3.378544E-01;       daa[17*20+11]= 8.024263E-03;
    daa[17*20+12]= 1.011149E-01;      daa[17*20+13]= 2.199856E-01;      daa[17*20+14]= 5.516074E-03;        daa[17*20+15]= 1.385142E-01;
    daa[17*20+16]= 1.412361E-02;        daa[18*20+0]= 1.433528E-01;        daa[18*20+1]= 1.711315E-01;        daa[18*20+2]= 2.622763E+00;
    daa[18*20+3]= 9.078338E-01;         daa[18*20+4]= 7.741612E-01;       daa[18*20+5]= 2.737091E-02;        daa[18*20+6]= 1.240642E-01;
    daa[18*20+7]= 2.295842E-01;        daa[18*20+8]= 2.055414E+01;       daa[18*20+9]= 2.903165E-01;        daa[18*20+10]= 1.521320E-01;
    daa[18*20+11]= 7.109973E-02;        daa[18*20+12]= 2.246759E-03;       daa[18*20+13]= 7.074464E+00;     daa[18*20+14]= 1.992133E-01;
    daa[18*20+15]= 8.104751E-01;       daa[18*20+16]= 9.984255E-02;       daa[18*20+17]= 6.121284E-01;      daa[19*20+0]= 3.774477E+00;
    daa[19*20+1]= 1.366145E-01;        daa[19*20+2]= 4.931206E-02;        daa[19*20+3]= 4.076074E-01;         daa[19*20+4]= 2.243512E-02;
    daa[19*20+5]= 9.047737E-03;        daa[19*20+6]= 5.795409E-01;        daa[19*20+7]= 4.228200E-01;         daa[19*20+8]= 6.890244E+00;
    daa[19*20+9]= 7.926675E+00;      daa[19*20+10]= 3.595310E+00;      daa[19*20+11]= 3.493440E-02;       daa[19*20+12]= 4.396720E+00;
    daa[19*20+13]= 1.643946E+00;       daa[19*20+14]= 2.217442E-01;       daa[19*20+15]= 7.477041E-02;       daa[19*20+16]= 2.166054E-01;
    daa[19*20+17]= 9.663569E-02;       daa[19*20+18]= 5.010635E-01;

    for (i=0; i<naa; i++)  for (j=0; j<i; j++)  daa[j*naa+i] = daa[i*naa+j];

    pi[0]= 6.541704E-02;     pi[1]= 4.708366E-02;     pi[2]= 3.168984E-02;     pi[3]= 4.688141E-02 ;
    pi[4]= 2.150693E-02;     pi[5]= 4.240711E-02;     pi[6]= 2.842211E-02;     pi[7]= 1.005278E-01;
    pi[8]= 9.812606E-03;     pi[9]= 3.424424E-02;     pi[10]= 6.222565E-02;    pi[11]= 4.844488E-02;
    pi[12]= 1.760370E-02;    pi[13]= 3.478555E-02;    pi[14]= 3.962469E-02;    pi[15]= 1.280566E-01;
    pi[16]= 8.199314E-02;    pi[17]= 3.393045E-02;    pi[18]= 7.586119E-02;    pi[19]= 4.948141E-02;
    
    return 1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void M4_Init_Model(m4 *m4mod, calign *data, t_mod *mod)
{
  int i,j,ct;
  phydbl fq;

  if(mod->io->datatype == NT)      m4mod->n_o = 4;
  else if(mod->io->datatype == AA) m4mod->n_o = 20;
  else
    {
      PhyML_Printf("\n== Not implemented yet.");
      PhyML_Printf("\n== Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  mod->ns = m4mod->n_o * m4mod->n_h;

  For(i,m4mod->n_o) m4mod->o_fq[i] = mod->e_frq->pi->v[i]; /*! At that stage, the mod->pi vector as been initialized
                                                             under a standard non covarion type of model. Use these
                                                             frequencies as they have been set according to the
                                                             nucleotide substitution model chosen (e.g., 1/4 for JC69). !*/
  For(i,(int)(m4mod->n_h)) m4mod->multipl[i] = 1.;
  
  ct = 0;
  For(i,m4mod->n_o-1)
    {
      for(j=i+1;j<m4mod->n_o;j++)
        {
          m4mod->o_rr[ct] = MAX(mod->r_mat->qmat->v[i*m4mod->n_o+j],1.E-5);
          ct++;
        }
    }

  For(i,(int)(m4mod->n_h*(m4mod->n_h-1)/2)) m4mod->h_rr[i] = 1.;
  fq = (phydbl)(1./m4mod->n_h);

  if(mod->s_opt->opt_cov_delta) m4mod->delta = 1.0;
  if(mod->s_opt->opt_cov_alpha) m4mod->alpha = 1.0;
  For(i,m4mod->n_h) m4mod->h_fq[i] = fq;
  For(i,m4mod->n_h) m4mod->h_fq_unscaled[i] = 1.0;
  For(i,m4mod->n_h) m4mod->multipl[i] = (phydbl)i;
  For(i,m4mod->n_h) m4mod->multipl_unscaled[i] = (phydbl)i;

  Switch_Eigen(YES,mod);
  M4_Update_Qmat(m4mod,mod);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void GEO_Init_Coord(t_geo_coord *t, int n_dim)
{
  t->dim = n_dim;
  Random_String(t->id,3);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void PHYREX_Init_Disk_Event(t_dsk *t, int n_dim, t_phyrex_mod *mmod)
{
  t->prev         = NULL;
  t->next         = NULL;
  t->mmod         = NULL;
  Random_String(t->id,3);

  GEO_Init_Coord(t->centr,n_dim);

  if(mmod != NULL) t->mmod = mmod;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void PHYREX_Init_Migrep_Mod(t_phyrex_mod *t, int n_dim, phydbl max_lat, phydbl max_lon)
{
  t->name             = PHYREX_NORMAL;
  t->n_dim            = n_dim;
  
  if(n_dim != 2) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  
  t->lim->lonlat[0]   = max_lat;
  t->lim->lonlat[1]   = max_lon;

  t->lbda             = 1.E-0;
  t->min_lbda         = 1.E-3;
  t->max_lbda         = 1.E+2;
  t->prior_param_lbda = 1.0;

  t->mu               = 0.300;
  t->min_mu           = 0.001;
  t->max_mu           = 1.000;
  t->prior_param_mu   = 0.500;

  t->min_rad           = 0.0;
  t->max_rad           = 1.0*(max_lat+max_lon);
  t->rad               = .10*(max_lat+max_lon);
  t->prior_param_rad   = 0.5;
  t->update_rad        = NO;

  t->min_sigsq         = 0.0;
  t->max_sigsq         = 1.E+2;
  t->sigsq             = 1.0;
  t->prior_param_sigsq = 10.0;
  

  t->c_lnL           = UNLIKELY;
  t->c_ln_prior_rad  = UNLIKELY;
  t->c_ln_prior_lbda = UNLIKELY;
  t->c_ln_prior_mu   = UNLIKELY;

  t->soft_bound_area = 0.1;
  
  t->sampl_area      = 0.0;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void PHYREX_Init_Lindisk_Node(t_ldsk *t, t_dsk *disk, int n_dim)
{
  t->disk    = disk;
  /* disk->ldsk = t; */
  t->prev    = NULL;
  t->next    = NULL;  
  t->nd      = NULL;
  t->is_hit  = NO;
  t->n_next  = 0;
  GEO_Init_Coord(t->coord,    n_dim);
  GEO_Init_Coord(t->cpy_coord,n_dim);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MCMC_Init_MCMC_Struct(char *filename, option *io, t_mcmc *mcmc)
{
  int pid;

  mcmc->io               = io;
  mcmc->is               = NO;
  mcmc->use_data         = YES;
  mcmc->run              = 0;
  mcmc->sample_interval  = 1E+3;
  mcmc->chain_len        = 1E+7;
  mcmc->chain_len_burnin = 1E+4;
  mcmc->randomize        = YES;
  mcmc->norm_freq        = 1E+3;
  mcmc->max_tune         = 1.E+20;
  mcmc->min_tune         = 1.E-10;
  mcmc->print_every      = 2;
  mcmc->is_burnin        = NO;
  mcmc->nd_t_digits      = 4;
  mcmc->always_yes       = NO;
  mcmc->max_lag          = 1000;
  mcmc->sample_num       = 0;

  if(filename)
    {
      char *s;

      s = (char *)mCalloc(T_MAX_NAME,sizeof(char));

      strcpy(mcmc->out_filename,filename);
      pid = getpid();
      sprintf(mcmc->out_filename+strlen(mcmc->out_filename),".%d",pid);

      strcpy(s,mcmc->io->in_align_file);
      strcat(s,"_");
      strcat(s,mcmc->out_filename);
      strcat(s,".stats");
      mcmc->out_fp_stats = fopen(s,"w");

      strcpy(s,mcmc->io->in_align_file);
      strcat(s,"_");
      strcat(s,mcmc->out_filename);
      strcat(s,".trees");
      mcmc->out_fp_trees = fopen(s,"w");

      strcpy(s,mcmc->io->in_align_file);
      strcat(s,"_");
      strcat(s,mcmc->out_filename);
      strcat(s,".constree");
      mcmc->out_fp_constree = fopen(s,"w");

/*       strcpy(s,tree->mcmc->out_filename); */
/*       strcat(s,".means"); */
/*       tree->mcmc->out_fp_means = fopen(s,"w"); */

/*       strcpy(s,tree->mcmc->out_filename); */
/*       strcat(s,".lasts"); */
/*       tree->mcmc->out_fp_last  = fopen(s,"w"); */

      Free(s);
    }
  else 
    {
      mcmc->out_fp_stats = stderr;
      mcmc->out_fp_trees = stderr;
      /* tree->mcmc->out_fp_means = stderr; */
      /* tree->mcmc->out_fp_last  = stderr; */
    }
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

