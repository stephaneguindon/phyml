/*

PHYML :  a program that  computes maximum likelihood  phylogenies from
DNA or AA homologous sequences

Copyright (C) Stephane Guindon. Oct 2003 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/


#include "pars.h"


/*********************************************************/


int Pars(t_edge *b, t_tree *tree)
{
  int site,n_patterns;

  if(tree->is_mixt_tree == YES) 
    {
      MIXT_Pars(b,tree);
      return tree->c_pars;
    }

  n_patterns = tree->n_pattern;

  
  if(b == NULL)
    {
      Post_Order_Pars(tree->a_nodes[0],tree->a_nodes[0]->v[0],tree);
      if(tree->both_sides == YES) Pre_Order_Pars(tree->a_nodes[0],tree->a_nodes[0]->v[0],tree);
    }

  if(b == NULL) b = tree->a_nodes[0]->b[0];

  tree->c_pars = 0;
  for(site=0;site<n_patterns;site++)
    {
      tree->site_pars[site] = 0;
      tree->curr_site       = site;
      tree->site_pars[site] = Pars_Core(b,tree);
      tree->c_pars += tree->site_pars[site] * tree->data->wght[site];
      /* printf("\n. site %d pars: %d",site,tree->c_pars); */
    }

  return tree->c_pars;
}

/*********************************************************/

void Post_Order_Pars(t_node *a, t_node *d, t_tree *tree)
{
  int i,dir;

  dir = -1;

  if(d->tax) return;
  else
    {
      for(i=0;i<3;i++)
    {
      if(d->v[i] != a)
        Post_Order_Pars(d,d->v[i],tree);
      else dir = i;
    }
      Get_All_Partial_Pars(tree,d->b[dir],a,d);
    }
}

/*********************************************************/

void Pre_Order_Pars(t_node *a, t_node *d, t_tree *tree)
{
  int i;

  if(d->tax) return;
  else
    {
      for(i=0;i<3;i++)
        {
          if(d->v[i] != a)
            {
              Get_All_Partial_Pars(tree,d->b[i],d->v[i],d);
              Pre_Order_Pars(d,d->v[i],tree);
            }
        }
    }
}

/*********************************************************/

void Get_All_Partial_Pars(t_tree *tree, t_edge *b_fcus, t_node *a, t_node *d)
{
  Update_Partial_Pars(tree,b_fcus,d);
}

/*********************************************************/

void Site_Pars(t_tree *tree)
{
  tree->site_pars[tree->curr_site] = Pars_Core(tree->a_nodes[0]->b[0],tree);
}

/*********************************************************/

void Init_Partial_Pars_Tips(t_tree *tree)
{
  int curr_site,i,j;
  short int *state_v;
  int dim1;

  dim1 = tree->mod->ns;

  state_v = (short int *)mCalloc(tree->mod->ns,sizeof(short int));

  for(curr_site=0;curr_site<tree->data->crunch_len;curr_site++)
    {
      for(i=0;i<tree->n_otu;i++)
        {
          if(tree->a_nodes[i]->b[0]->rght->tax != 1)
            {
	      PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d\n",__FILE__,__LINE__);
              Exit("\n");
            }
          
          if(tree->io->datatype == NT)
            {
              Init_Tips_At_One_Site_Nucleotides_Int(tree->a_nodes[i]->c_seq->state[curr_site],
                                                    0,
                                                    state_v);
              for(j=0;j<tree->mod->ns;j++) tree->a_nodes[i]->b[0]->p_pars_r[curr_site*dim1+j] = MAX_PARS;
              for(j=0;j<tree->mod->ns;j++) if(state_v[j] > 0.5) tree->a_nodes[i]->b[0]->p_pars_r[curr_site*dim1+j] =  0;
            }
          else if(tree->io->datatype == AA)
            {
              Init_Tips_At_One_Site_AA_Int(tree->a_nodes[i]->c_seq->state[curr_site],
                                           0,
                                           state_v);
              for(j=0;j<tree->mod->ns;j++) tree->a_nodes[i]->b[0]->p_pars_r[curr_site*dim1+j] = MAX_PARS;
              for(j=0;j<tree->mod->ns;j++) if(state_v[j] > 0.5) tree->a_nodes[i]->b[0]->p_pars_r[curr_site*dim1+j] =  0;
            }
          else if(tree->io->datatype == GENERIC)
            {
              Init_Tips_At_One_Site_Generic_Int(tree->a_nodes[i]->c_seq->state+curr_site*tree->mod->io->state_len,
                                                tree->mod->ns,
                                                tree->mod->io->state_len,
                                                0,
                                                state_v);
              for(j=0;j<tree->mod->ns;j++) tree->a_nodes[i]->b[0]->p_pars_r[curr_site*dim1+j] = MAX_PARS;
              for(j=0;j<tree->mod->ns;j++) if(state_v[j] > 0.5) tree->a_nodes[i]->b[0]->p_pars_r[curr_site*dim1+j] =  0;
            }
        }
    }
  Free(state_v);
}

/*********************************************************/

void Init_Ui_Tips(t_tree *tree)
{
  int curr_site,i,j,br;
  short int *state_v;

  state_v = (short int *)mCalloc(tree->mod->ns,sizeof(short int));

  for(curr_site=0;curr_site<tree->data->crunch_len;curr_site++)
    {
      for(i=0;i<tree->n_otu;++i)
        {
          if(tree->io->datatype == NT)
            {
              if(tree->a_nodes[i]->b[0]->rght->tax != 1)
                {
                  PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d\n",__FILE__,__LINE__);
                  Exit("\n");
                }
              
              Init_Tips_At_One_Site_Nucleotides_Int(tree->a_nodes[i]->c_seq->state[curr_site],
                                                    0,
                                                    state_v);
              /* Init_Tips_At_One_Site_Nucleotides_Int(tree->data->c_seq[i]->state[curr_site], */
              /* 					    0, */
              /* 					    state_v);	       */
              tree->a_nodes[i]->b[0]->ui_r[curr_site] = 0;
              for(j=0;j<tree->mod->ns;j++) tree->a_nodes[i]->b[0]->ui_r[curr_site] += (int)(state_v[j] * POW(2,j));
            }
          else if(tree->io->datatype == AA)
            {
              Init_Tips_At_One_Site_AA_Int(tree->a_nodes[i]->c_seq->state[curr_site],
                                           0,
                                           state_v);
              /* Init_Tips_At_One_Site_AA_Int(tree->data->c_seq[i]->state[curr_site], */
              /* 				   0, */
              /* 				   state_v); */
              tree->a_nodes[i]->b[0]->ui_r[curr_site] = 0;
              for(j=0;j<tree->mod->ns;j++) tree->a_nodes[i]->b[0]->ui_r[curr_site] += (int)(state_v[j] * POW(2,j));
            }
          else if(tree->io->datatype == GENERIC)
            {
              Init_Tips_At_One_Site_Generic_Int(tree->a_nodes[i]->c_seq->state+curr_site*tree->mod->io->state_len,
                                                tree->mod->ns,
                                                tree->mod->io->state_len,
                                                0,
                                                state_v);
              /* Init_Tips_At_One_Site_Generic_Int(tree->data->c_seq[i]->state+curr_site*tree->mod->io->state_len, */
              /* 					tree->mod->ns, */
              /* 					tree->mod->io->state_len, */
              /* 					0, */
              /* 					state_v); */
              tree->a_nodes[i]->b[0]->ui_r[curr_site] = 0;
              for(j=0;j<tree->mod->ns;j++) tree->a_nodes[i]->b[0]->ui_r[curr_site] += (int)(state_v[j] * POW(2,j));
            }
        }
    }
  
  
  For(br,2*tree->n_otu-3)
    {
      for(curr_site=0;curr_site<tree->data->crunch_len;curr_site++)
        {
          tree->a_edges[br]->pars_r[curr_site] = 0;
          tree->a_edges[br]->pars_l[curr_site] = 0;
        }
    }
  
  
  Free(state_v);
}

/*********************************************************/

void Update_Partial_Pars(t_tree *tree, t_edge *b_fcus, t_node *n)
{
/*
           |
           |<- b_fcus
           |
           n
          / \
         /   \
        /     \
*/

  unsigned int i,j;
  unsigned int site;
  int *ui, *ui_v1, *ui_v2;
  int *p_pars_v1, *p_pars_v2, *p_pars;
  int *pars, *pars_v1, *pars_v2;
  int min_v1,min_v2;
  int v;
  
  const unsigned int ns = tree->mod->ns;
  const unsigned int n_patterns = tree->n_pattern;

  if(tree->is_mixt_tree)
    {
      MIXT_Update_Partial_Pars(tree,b_fcus,n);
      return;
    }

  if((tree->io->do_alias_subpatt == YES) &&
     (tree->update_alias_subpatt == YES))
    Alias_One_Subpatt((n==b_fcus->left)?(b_fcus->rght):(b_fcus->left),n,tree);

  if(n->tax) return;

  ui = ui_v1 = ui_v2 = NULL;
  p_pars = p_pars_v1 = p_pars_v2 = NULL;
  pars = pars_v1 = pars_v2 = NULL;

    
  if(n == b_fcus->left)
    {
      ui = b_fcus->ui_l;
      
      pars = b_fcus->pars_l;
      p_pars = b_fcus->p_pars_l;
      
      ui_v1 =
        (n == n->b[b_fcus->l_v1]->left)?
        (n->b[b_fcus->l_v1]->ui_r):
        (n->b[b_fcus->l_v1]->ui_l);
      
      ui_v2 =
        (n == n->b[b_fcus->l_v2]->left)?
        (n->b[b_fcus->l_v2]->ui_r):
        (n->b[b_fcus->l_v2]->ui_l);
      
      p_pars_v1 =
        (n == n->b[b_fcus->l_v1]->left)?
        (n->b[b_fcus->l_v1]->p_pars_r):
        (n->b[b_fcus->l_v1]->p_pars_l);
      
      p_pars_v2 =
        (n == n->b[b_fcus->l_v2]->left)?
        (n->b[b_fcus->l_v2]->p_pars_r):
        (n->b[b_fcus->l_v2]->p_pars_l);
      
      pars_v1 =
        (n == n->b[b_fcus->l_v1]->left)?
        (n->b[b_fcus->l_v1]->pars_r):
        (n->b[b_fcus->l_v1]->pars_l);
      
      pars_v2 =
        (n == n->b[b_fcus->l_v2]->left)?
        (n->b[b_fcus->l_v2]->pars_r):
        (n->b[b_fcus->l_v2]->pars_l);
    }
  else
    {
      ui = b_fcus->ui_r;
      
      pars = b_fcus->pars_r;
      p_pars = b_fcus->p_pars_r;
      
      ui_v1 =
        (n == n->b[b_fcus->r_v1]->left)?
        (n->b[b_fcus->r_v1]->ui_r):
        (n->b[b_fcus->r_v1]->ui_l);
      
      ui_v2 =
        (n == n->b[b_fcus->r_v2]->left)?
        (n->b[b_fcus->r_v2]->ui_r):
        (n->b[b_fcus->r_v2]->ui_l);
      
      p_pars_v1 =
        (n == n->b[b_fcus->r_v1]->left)?
        (n->b[b_fcus->r_v1]->p_pars_r):
        (n->b[b_fcus->r_v1]->p_pars_l);
      
      p_pars_v2 =
        (n == n->b[b_fcus->r_v2]->left)?
        (n->b[b_fcus->r_v2]->p_pars_r):
        (n->b[b_fcus->r_v2]->p_pars_l);
      
      pars_v1 =
        (n == n->b[b_fcus->r_v1]->left)?
        (n->b[b_fcus->r_v1]->pars_r):
        (n->b[b_fcus->r_v1]->pars_l);
      
      pars_v2 =
        (n == n->b[b_fcus->r_v2]->left)?
        (n->b[b_fcus->r_v2]->pars_r):
        (n->b[b_fcus->r_v2]->pars_l);
    }
  
  
  if(tree->mod->s_opt->general_pars)
    {
      for(site=0;site<n_patterns;++site)
        {
          for(i=0;i<ns;++i)
            {
              min_v1 = MAX_PARS;
              for(j=0;j<ns;++j)
                {
                  v = p_pars_v1[site*ns+j] + tree->step_mat[i*ns+j];
                  if(v < min_v1) min_v1 = v;
                }
              
              min_v2 = MAX_PARS;
              for(j=0;j<ns;++j)
                {
                  v = p_pars_v2[site*ns+j] + tree->step_mat[i*ns+j];
                  if(v < min_v2) min_v2 = v;
                }
              p_pars[site*ns+i] = min_v1 + min_v2;
            }
        }
    }
  else
    {
      for(site=0;site<n_patterns;++site)
        {
          pars[site] = pars_v1[site] + pars_v2[site];
          
          ui[site] = ui_v1[site] & ui_v2[site];
          
          if(!ui[site])
            {
              pars[site]++;
              ui[site] = ui_v1[site] | ui_v2[site];
            }          
        }
    }
}

/*********************************************************/

int Pars_Core(t_edge *b, t_tree *tree)
{
  int site;
  unsigned int i,j;
  int site_pars;
  int min_l,min_r;
  int v;
  
  const unsigned int ns = tree->mod->ns;
  site = tree->curr_site;
  site_pars = MAX_PARS;

  if(tree->mod->s_opt->general_pars)
    {
      for(i=0;i<ns;++i)
        {
          min_l = MAX_PARS;
          for(j=0;j<ns;++j)
            {
              v = b->p_pars_l[site*ns+j] + tree->step_mat[i*ns+j];
              if(v < min_l) min_l = v;
            }
          
          min_r = MAX_PARS;
          for(j=0;j<ns;++j)
            {
              v = b->p_pars_r[site*ns+j] + tree->step_mat[i*ns+j];
              if(v < min_r) min_r = v;
            }
          
          if((min_l + min_r) < site_pars) site_pars = min_l + min_r;
        }
    }
  else
    {
      site_pars = b->pars_l[site] + b->pars_r[site];
      if(!(b->ui_l[site] & b->ui_r[site])) site_pars++;
    }
  
  return site_pars;
}

/*********************************************************/
/* Is there one or more parsimoniy step(s) along this t_edge ?
   0 -> NO; 1 -> YES
*/
int One_Pars_Step(t_edge *b,t_tree *tree)
{
  int site;
  int init_general_pars;

  init_general_pars = tree->mod->s_opt->general_pars;

  tree->mod->s_opt->general_pars = 0;
  Set_Both_Sides(YES,tree);
  Pars(NULL,tree);

  for(site=0;site<tree->n_pattern;site++)
    {
      if(!(b->ui_l[site] & b->ui_r[site])) break;
    }
  tree->mod->s_opt->general_pars = init_general_pars;
  if(site == tree->n_pattern) return 0;
  else
    {
      PhyML_Printf("\n. One parsimony step ocurred at site %4d",site);
      return 1;
    }
}

/*********************************************************/
int Pars_At_Given_Edge(t_edge *b, t_tree *tree)
{
  int site,n_patterns;

/*   n_patterns = (int)FLOOR(tree->n_pattern*tree->prop_of_sites_to_consider); */
  n_patterns = tree->n_pattern;

  tree->c_pars = .0;
  for(site=0;site<n_patterns;site++)
    {
      tree->site_pars[site] = 0;
      tree->curr_site = site;
      tree->site_pars[site] = Pars_Core(b,tree);
      tree->c_pars += tree->site_pars[site] * tree->data->wght[site];
    }
  return tree->c_pars;
}

/*********************************************************/

int Update_Pars_At_Given_Edge(t_edge *b_fcus, t_tree *tree)
{
  Update_Partial_Pars(tree,b_fcus,b_fcus->left);
  Update_Partial_Pars(tree,b_fcus,b_fcus->rght);
  tree->c_pars = Pars(b_fcus,tree);
  return tree->c_pars;
}

/*********************************************************/

void Get_Step_Mat(t_tree *tree)
{
  int i;

  if(tree->io->datatype == AA)
    {
      tree->step_mat[ 0*tree->mod->ns+ 0] =    0 ;
      tree->step_mat[ 0*tree->mod->ns+ 1] =    3 ;
      tree->step_mat[ 0*tree->mod->ns+ 2] =    3 ;
      tree->step_mat[ 0*tree->mod->ns+ 3] =    2 ;
      tree->step_mat[ 0*tree->mod->ns+ 4] =    3 ;
      tree->step_mat[ 0*tree->mod->ns+ 5] =    3 ;
      tree->step_mat[ 0*tree->mod->ns+ 6] =    2 ;
      tree->step_mat[ 0*tree->mod->ns+ 7] =    2 ;
      tree->step_mat[ 0*tree->mod->ns+ 8] =    3 ;
      tree->step_mat[ 0*tree->mod->ns+ 9] =    3 ;
      tree->step_mat[ 0*tree->mod->ns+10] =    3 ;
      tree->step_mat[ 0*tree->mod->ns+11] =    3 ;
      tree->step_mat[ 0*tree->mod->ns+12] =    3 ;
      tree->step_mat[ 0*tree->mod->ns+13] =    3 ;
      tree->step_mat[ 0*tree->mod->ns+14] =    2 ;
      tree->step_mat[ 0*tree->mod->ns+15] =    2 ;
      tree->step_mat[ 0*tree->mod->ns+16] =    2 ;
      tree->step_mat[ 0*tree->mod->ns+17] =    3 ;
      tree->step_mat[ 0*tree->mod->ns+18] =    3 ;
      tree->step_mat[ 0*tree->mod->ns+19] =    2 ;
      tree->step_mat[ 1*tree->mod->ns+ 0] =    3 ;
      tree->step_mat[ 1*tree->mod->ns+ 1] =    0 ;
      tree->step_mat[ 1*tree->mod->ns+ 2] =    2 ;
      tree->step_mat[ 1*tree->mod->ns+ 3] =    3 ;
      tree->step_mat[ 1*tree->mod->ns+ 4] =    2 ;
      tree->step_mat[ 1*tree->mod->ns+ 5] =    2 ;
      tree->step_mat[ 1*tree->mod->ns+ 6] =    3 ;
      tree->step_mat[ 1*tree->mod->ns+ 7] =    2 ;
      tree->step_mat[ 1*tree->mod->ns+ 8] =    2 ;
      tree->step_mat[ 1*tree->mod->ns+ 9] =    2 ;
      tree->step_mat[ 1*tree->mod->ns+10] =    2 ;
      tree->step_mat[ 1*tree->mod->ns+11] =    2 ;
      tree->step_mat[ 1*tree->mod->ns+12] =    2 ;
      tree->step_mat[ 1*tree->mod->ns+13] =    3 ;
      tree->step_mat[ 1*tree->mod->ns+14] =    2 ;
      tree->step_mat[ 1*tree->mod->ns+15] =    2 ;
      tree->step_mat[ 1*tree->mod->ns+16] =    2 ;
      tree->step_mat[ 1*tree->mod->ns+17] =    2 ;
      tree->step_mat[ 1*tree->mod->ns+18] =    3 ;
      tree->step_mat[ 1*tree->mod->ns+19] =    3 ;
      tree->step_mat[ 2*tree->mod->ns+ 0] =    3 ;
      tree->step_mat[ 2*tree->mod->ns+ 1] =    2 ;
      tree->step_mat[ 2*tree->mod->ns+ 2] =    0 ;
      tree->step_mat[ 2*tree->mod->ns+ 3] =    2 ;
      tree->step_mat[ 2*tree->mod->ns+ 4] =    2 ;
      tree->step_mat[ 2*tree->mod->ns+ 5] =    2 ;
      tree->step_mat[ 2*tree->mod->ns+ 6] =    2 ;
      tree->step_mat[ 2*tree->mod->ns+ 7] =    3 ;
      tree->step_mat[ 2*tree->mod->ns+ 8] =    2 ;
      tree->step_mat[ 2*tree->mod->ns+ 9] =    2 ;
      tree->step_mat[ 2*tree->mod->ns+10] =    3 ;
      tree->step_mat[ 2*tree->mod->ns+11] =    1 ;
      tree->step_mat[ 2*tree->mod->ns+12] =    2 ;
      tree->step_mat[ 2*tree->mod->ns+13] =    2 ;
      tree->step_mat[ 2*tree->mod->ns+14] =    3 ;
      tree->step_mat[ 2*tree->mod->ns+15] =    2 ;
      tree->step_mat[ 2*tree->mod->ns+16] =    2 ;
      tree->step_mat[ 2*tree->mod->ns+17] =    3 ;
      tree->step_mat[ 2*tree->mod->ns+18] =    2 ;
      tree->step_mat[ 2*tree->mod->ns+19] =    3 ;
      tree->step_mat[ 3*tree->mod->ns+ 0] =    2 ;
      tree->step_mat[ 3*tree->mod->ns+ 1] =    3 ;
      tree->step_mat[ 3*tree->mod->ns+ 2] =    2 ;
      tree->step_mat[ 3*tree->mod->ns+ 3] =    0 ;
      tree->step_mat[ 3*tree->mod->ns+ 4] =    2 ;
      tree->step_mat[ 3*tree->mod->ns+ 5] =    2 ;
      tree->step_mat[ 3*tree->mod->ns+ 6] =    1 ;
      tree->step_mat[ 3*tree->mod->ns+ 7] =    2 ;
      tree->step_mat[ 3*tree->mod->ns+ 8] =    2 ;
      tree->step_mat[ 3*tree->mod->ns+ 9] =    3 ;
      tree->step_mat[ 3*tree->mod->ns+10] =    3 ;
      tree->step_mat[ 3*tree->mod->ns+11] =    2 ;
      tree->step_mat[ 3*tree->mod->ns+12] =    3 ;
      tree->step_mat[ 3*tree->mod->ns+13] =    2 ;
      tree->step_mat[ 3*tree->mod->ns+14] =    3 ;
      tree->step_mat[ 3*tree->mod->ns+15] =    3 ;
      tree->step_mat[ 3*tree->mod->ns+16] =    3 ;
      tree->step_mat[ 3*tree->mod->ns+17] =    3 ;
      tree->step_mat[ 3*tree->mod->ns+18] =    2 ;
      tree->step_mat[ 3*tree->mod->ns+19] =    2 ;
      tree->step_mat[ 4*tree->mod->ns+ 0] =    3 ;
      tree->step_mat[ 4*tree->mod->ns+ 1] =    2 ;
      tree->step_mat[ 4*tree->mod->ns+ 2] =    2 ;
      tree->step_mat[ 4*tree->mod->ns+ 3] =    2 ;
      tree->step_mat[ 4*tree->mod->ns+ 4] =    0 ;
      tree->step_mat[ 4*tree->mod->ns+ 5] =    3 ;
      tree->step_mat[ 4*tree->mod->ns+ 6] =    3 ;
      tree->step_mat[ 4*tree->mod->ns+ 7] =    2 ;
      tree->step_mat[ 4*tree->mod->ns+ 8] =    2 ;
      tree->step_mat[ 4*tree->mod->ns+ 9] =    3 ;
      tree->step_mat[ 4*tree->mod->ns+10] =    2 ;
      tree->step_mat[ 4*tree->mod->ns+11] =    3 ;
      tree->step_mat[ 4*tree->mod->ns+12] =    3 ;
      tree->step_mat[ 4*tree->mod->ns+13] =    2 ;
      tree->step_mat[ 4*tree->mod->ns+14] =    3 ;
      tree->step_mat[ 4*tree->mod->ns+15] =    2 ;
      tree->step_mat[ 4*tree->mod->ns+16] =    3 ;
      tree->step_mat[ 4*tree->mod->ns+17] =    1 ;
      tree->step_mat[ 4*tree->mod->ns+18] =    2 ;
      tree->step_mat[ 4*tree->mod->ns+19] =    3 ;
      tree->step_mat[ 5*tree->mod->ns+ 0] =    3 ;
      tree->step_mat[ 5*tree->mod->ns+ 1] =    2 ;
      tree->step_mat[ 5*tree->mod->ns+ 2] =    2 ;
      tree->step_mat[ 5*tree->mod->ns+ 3] =    2 ;
      tree->step_mat[ 5*tree->mod->ns+ 4] =    3 ;
      tree->step_mat[ 5*tree->mod->ns+ 5] =    0 ;
      tree->step_mat[ 5*tree->mod->ns+ 6] =    2 ;
      tree->step_mat[ 5*tree->mod->ns+ 7] =    3 ;
      tree->step_mat[ 5*tree->mod->ns+ 8] =    1 ;
      tree->step_mat[ 5*tree->mod->ns+ 9] =    3 ;
      tree->step_mat[ 5*tree->mod->ns+10] =    2 ;
      tree->step_mat[ 5*tree->mod->ns+11] =    2 ;
      tree->step_mat[ 5*tree->mod->ns+12] =    2 ;
      tree->step_mat[ 5*tree->mod->ns+13] =    3 ;
      tree->step_mat[ 5*tree->mod->ns+14] =    2 ;
      tree->step_mat[ 5*tree->mod->ns+15] =    3 ;
      tree->step_mat[ 5*tree->mod->ns+16] =    3 ;
      tree->step_mat[ 5*tree->mod->ns+17] =    2 ;
      tree->step_mat[ 5*tree->mod->ns+18] =    2 ;
      tree->step_mat[ 5*tree->mod->ns+19] =    3 ;
      tree->step_mat[ 6*tree->mod->ns+ 0] =    2 ;
      tree->step_mat[ 6*tree->mod->ns+ 1] =    3 ;
      tree->step_mat[ 6*tree->mod->ns+ 2] =    2 ;
      tree->step_mat[ 6*tree->mod->ns+ 3] =    1 ;
      tree->step_mat[ 6*tree->mod->ns+ 4] =    3 ;
      tree->step_mat[ 6*tree->mod->ns+ 5] =    2 ;
      tree->step_mat[ 6*tree->mod->ns+ 6] =    0 ;
      tree->step_mat[ 6*tree->mod->ns+ 7] =    2 ;
      tree->step_mat[ 6*tree->mod->ns+ 8] =    2 ;
      tree->step_mat[ 6*tree->mod->ns+ 9] =    3 ;
      tree->step_mat[ 6*tree->mod->ns+10] =    3 ;
      tree->step_mat[ 6*tree->mod->ns+11] =    2 ;
      tree->step_mat[ 6*tree->mod->ns+12] =    2 ;
      tree->step_mat[ 6*tree->mod->ns+13] =    3 ;
      tree->step_mat[ 6*tree->mod->ns+14] =    3 ;
      tree->step_mat[ 6*tree->mod->ns+15] =    3 ;
      tree->step_mat[ 6*tree->mod->ns+16] =    3 ;
      tree->step_mat[ 6*tree->mod->ns+17] =    2 ;
      tree->step_mat[ 6*tree->mod->ns+18] =    2 ;
      tree->step_mat[ 6*tree->mod->ns+19] =    2 ;
      tree->step_mat[ 7*tree->mod->ns+ 0] =    2 ;
      tree->step_mat[ 7*tree->mod->ns+ 1] =    2 ;
      tree->step_mat[ 7*tree->mod->ns+ 2] =    3 ;
      tree->step_mat[ 7*tree->mod->ns+ 3] =    2 ;
      tree->step_mat[ 7*tree->mod->ns+ 4] =    2 ;
      tree->step_mat[ 7*tree->mod->ns+ 5] =    3 ;
      tree->step_mat[ 7*tree->mod->ns+ 6] =    2 ;
      tree->step_mat[ 7*tree->mod->ns+ 7] =    0 ;
      tree->step_mat[ 7*tree->mod->ns+ 8] =    3 ;
      tree->step_mat[ 7*tree->mod->ns+ 9] =    3 ;
      tree->step_mat[ 7*tree->mod->ns+10] =    3 ;
      tree->step_mat[ 7*tree->mod->ns+11] =    3 ;
      tree->step_mat[ 7*tree->mod->ns+12] =    3 ;
      tree->step_mat[ 7*tree->mod->ns+13] =    3 ;
      tree->step_mat[ 7*tree->mod->ns+14] =    3 ;
      tree->step_mat[ 7*tree->mod->ns+15] =    2 ;
      tree->step_mat[ 7*tree->mod->ns+16] =    3 ;
      tree->step_mat[ 7*tree->mod->ns+17] =    2 ;
      tree->step_mat[ 7*tree->mod->ns+18] =    3 ;
      tree->step_mat[ 7*tree->mod->ns+19] =    2 ;
      tree->step_mat[ 8*tree->mod->ns+ 0] =    3 ;
      tree->step_mat[ 8*tree->mod->ns+ 1] =    2 ;
      tree->step_mat[ 8*tree->mod->ns+ 2] =    2 ;
      tree->step_mat[ 8*tree->mod->ns+ 3] =    2 ;
      tree->step_mat[ 8*tree->mod->ns+ 4] =    2 ;
      tree->step_mat[ 8*tree->mod->ns+ 5] =    1 ;
      tree->step_mat[ 8*tree->mod->ns+ 6] =    2 ;
      tree->step_mat[ 8*tree->mod->ns+ 7] =    3 ;
      tree->step_mat[ 8*tree->mod->ns+ 8] =    0 ;
      tree->step_mat[ 8*tree->mod->ns+ 9] =    3 ;
      tree->step_mat[ 8*tree->mod->ns+10] =    2 ;
      tree->step_mat[ 8*tree->mod->ns+11] =    2 ;
      tree->step_mat[ 8*tree->mod->ns+12] =    3 ;
      tree->step_mat[ 8*tree->mod->ns+13] =    2 ;
      tree->step_mat[ 8*tree->mod->ns+14] =    2 ;
      tree->step_mat[ 8*tree->mod->ns+15] =    3 ;
      tree->step_mat[ 8*tree->mod->ns+16] =    3 ;
      tree->step_mat[ 8*tree->mod->ns+17] =    3 ;
      tree->step_mat[ 8*tree->mod->ns+18] =    2 ;
      tree->step_mat[ 8*tree->mod->ns+19] =    3 ;
      tree->step_mat[ 9*tree->mod->ns+ 0] =    3 ;
      tree->step_mat[ 9*tree->mod->ns+ 1] =    2 ;
      tree->step_mat[ 9*tree->mod->ns+ 2] =    2 ;
      tree->step_mat[ 9*tree->mod->ns+ 3] =    3 ;
      tree->step_mat[ 9*tree->mod->ns+ 4] =    3 ;
      tree->step_mat[ 9*tree->mod->ns+ 5] =    3 ;
      tree->step_mat[ 9*tree->mod->ns+ 6] =    3 ;
      tree->step_mat[ 9*tree->mod->ns+ 7] =    3 ;
      tree->step_mat[ 9*tree->mod->ns+ 8] =    3 ;
      tree->step_mat[ 9*tree->mod->ns+ 9] =    0 ;
      tree->step_mat[ 9*tree->mod->ns+10] =    2 ;
      tree->step_mat[ 9*tree->mod->ns+11] =    2 ;
      tree->step_mat[ 9*tree->mod->ns+12] =    1 ;
      tree->step_mat[ 9*tree->mod->ns+13] =    2 ;
      tree->step_mat[ 9*tree->mod->ns+14] =    3 ;
      tree->step_mat[ 9*tree->mod->ns+15] =    2 ;
      tree->step_mat[ 9*tree->mod->ns+16] =    2 ;
      tree->step_mat[ 9*tree->mod->ns+17] =    3 ;
      tree->step_mat[ 9*tree->mod->ns+18] =    3 ;
      tree->step_mat[ 9*tree->mod->ns+19] =    2 ;
      tree->step_mat[10*tree->mod->ns+ 0] =    3 ;
      tree->step_mat[10*tree->mod->ns+ 1] =    2 ;
      tree->step_mat[10*tree->mod->ns+ 2] =    3 ;
      tree->step_mat[10*tree->mod->ns+ 3] =    3 ;
      tree->step_mat[10*tree->mod->ns+ 4] =    2 ;
      tree->step_mat[10*tree->mod->ns+ 5] =    2 ;
      tree->step_mat[10*tree->mod->ns+ 6] =    3 ;
      tree->step_mat[10*tree->mod->ns+ 7] =    3 ;
      tree->step_mat[10*tree->mod->ns+ 8] =    2 ;
      tree->step_mat[10*tree->mod->ns+ 9] =    2 ;
      tree->step_mat[10*tree->mod->ns+10] =    0 ;
      tree->step_mat[10*tree->mod->ns+11] =    3 ;
      tree->step_mat[10*tree->mod->ns+12] =    2 ;
      tree->step_mat[10*tree->mod->ns+13] =    2 ;
      tree->step_mat[10*tree->mod->ns+14] =    2 ;
      tree->step_mat[10*tree->mod->ns+15] =    3 ;
      tree->step_mat[10*tree->mod->ns+16] =    3 ;
      tree->step_mat[10*tree->mod->ns+17] =    2 ;
      tree->step_mat[10*tree->mod->ns+18] =    2 ;
      tree->step_mat[10*tree->mod->ns+19] =    2 ;
      tree->step_mat[11*tree->mod->ns+ 0] =    3 ;
      tree->step_mat[11*tree->mod->ns+ 1] =    2 ;
      tree->step_mat[11*tree->mod->ns+ 2] =    1 ;
      tree->step_mat[11*tree->mod->ns+ 3] =    2 ;
      tree->step_mat[11*tree->mod->ns+ 4] =    3 ;
      tree->step_mat[11*tree->mod->ns+ 5] =    2 ;
      tree->step_mat[11*tree->mod->ns+ 6] =    2 ;
      tree->step_mat[11*tree->mod->ns+ 7] =    3 ;
      tree->step_mat[11*tree->mod->ns+ 8] =    2 ;
      tree->step_mat[11*tree->mod->ns+ 9] =    2 ;
      tree->step_mat[11*tree->mod->ns+10] =    3 ;
      tree->step_mat[11*tree->mod->ns+11] =    0 ;
      tree->step_mat[11*tree->mod->ns+12] =    2 ;
      tree->step_mat[11*tree->mod->ns+13] =    3 ;
      tree->step_mat[11*tree->mod->ns+14] =    3 ;
      tree->step_mat[11*tree->mod->ns+15] =    2 ;
      tree->step_mat[11*tree->mod->ns+16] =    2 ;
      tree->step_mat[11*tree->mod->ns+17] =    2 ;
      tree->step_mat[11*tree->mod->ns+18] =    2 ;
      tree->step_mat[11*tree->mod->ns+19] =    3 ;
      tree->step_mat[12*tree->mod->ns+ 0] =    3 ;
      tree->step_mat[12*tree->mod->ns+ 1] =    2 ;
      tree->step_mat[12*tree->mod->ns+ 2] =    2 ;
      tree->step_mat[12*tree->mod->ns+ 3] =    3 ;
      tree->step_mat[12*tree->mod->ns+ 4] =    3 ;
      tree->step_mat[12*tree->mod->ns+ 5] =    2 ;
      tree->step_mat[12*tree->mod->ns+ 6] =    2 ;
      tree->step_mat[12*tree->mod->ns+ 7] =    3 ;
      tree->step_mat[12*tree->mod->ns+ 8] =    3 ;
      tree->step_mat[12*tree->mod->ns+ 9] =    1 ;
      tree->step_mat[12*tree->mod->ns+10] =    2 ;
      tree->step_mat[12*tree->mod->ns+11] =    2 ;
      tree->step_mat[12*tree->mod->ns+12] =    0 ;
      tree->step_mat[12*tree->mod->ns+13] =    2 ;
      tree->step_mat[12*tree->mod->ns+14] =    3 ;
      tree->step_mat[12*tree->mod->ns+15] =    2 ;
      tree->step_mat[12*tree->mod->ns+16] =    2 ;
      tree->step_mat[12*tree->mod->ns+17] =    2 ;
      tree->step_mat[12*tree->mod->ns+18] =    3 ;
      tree->step_mat[12*tree->mod->ns+19] =    2 ;
      tree->step_mat[13*tree->mod->ns+ 0] =    3 ;
      tree->step_mat[13*tree->mod->ns+ 1] =    3 ;
      tree->step_mat[13*tree->mod->ns+ 2] =    2 ;
      tree->step_mat[13*tree->mod->ns+ 3] =    2 ;
      tree->step_mat[13*tree->mod->ns+ 4] =    2 ;
      tree->step_mat[13*tree->mod->ns+ 5] =    3 ;
      tree->step_mat[13*tree->mod->ns+ 6] =    3 ;
      tree->step_mat[13*tree->mod->ns+ 7] =    3 ;
      tree->step_mat[13*tree->mod->ns+ 8] =    2 ;
      tree->step_mat[13*tree->mod->ns+ 9] =    2 ;
      tree->step_mat[13*tree->mod->ns+10] =    2 ;
      tree->step_mat[13*tree->mod->ns+11] =    3 ;
      tree->step_mat[13*tree->mod->ns+12] =    2 ;
      tree->step_mat[13*tree->mod->ns+13] =    0 ;
      tree->step_mat[13*tree->mod->ns+14] =    3 ;
      tree->step_mat[13*tree->mod->ns+15] =    2 ;
      tree->step_mat[13*tree->mod->ns+16] =    3 ;
      tree->step_mat[13*tree->mod->ns+17] =    2 ;
      tree->step_mat[13*tree->mod->ns+18] =    2 ;
      tree->step_mat[13*tree->mod->ns+19] =    2 ;
      tree->step_mat[14*tree->mod->ns+ 0] =    2 ;
      tree->step_mat[14*tree->mod->ns+ 1] =    2 ;
      tree->step_mat[14*tree->mod->ns+ 2] =    3 ;
      tree->step_mat[14*tree->mod->ns+ 3] =    3 ;
      tree->step_mat[14*tree->mod->ns+ 4] =    3 ;
      tree->step_mat[14*tree->mod->ns+ 5] =    2 ;
      tree->step_mat[14*tree->mod->ns+ 6] =    3 ;
      tree->step_mat[14*tree->mod->ns+ 7] =    3 ;
      tree->step_mat[14*tree->mod->ns+ 8] =    2 ;
      tree->step_mat[14*tree->mod->ns+ 9] =    3 ;
      tree->step_mat[14*tree->mod->ns+10] =    2 ;
      tree->step_mat[14*tree->mod->ns+11] =    3 ;
      tree->step_mat[14*tree->mod->ns+12] =    3 ;
      tree->step_mat[14*tree->mod->ns+13] =    3 ;
      tree->step_mat[14*tree->mod->ns+14] =    0 ;
      tree->step_mat[14*tree->mod->ns+15] =    2 ;
      tree->step_mat[14*tree->mod->ns+16] =    2 ;
      tree->step_mat[14*tree->mod->ns+17] =    3 ;
      tree->step_mat[14*tree->mod->ns+18] =    3 ;
      tree->step_mat[14*tree->mod->ns+19] =    3 ;
      tree->step_mat[15*tree->mod->ns+ 0] =    2 ;
      tree->step_mat[15*tree->mod->ns+ 1] =    2 ;
      tree->step_mat[15*tree->mod->ns+ 2] =    2 ;
      tree->step_mat[15*tree->mod->ns+ 3] =    3 ;
      tree->step_mat[15*tree->mod->ns+ 4] =    2 ;
      tree->step_mat[15*tree->mod->ns+ 5] =    3 ;
      tree->step_mat[15*tree->mod->ns+ 6] =    3 ;
      tree->step_mat[15*tree->mod->ns+ 7] =    2 ;
      tree->step_mat[15*tree->mod->ns+ 8] =    3 ;
      tree->step_mat[15*tree->mod->ns+ 9] =    2 ;
      tree->step_mat[15*tree->mod->ns+10] =    3 ;
      tree->step_mat[15*tree->mod->ns+11] =    2 ;
      tree->step_mat[15*tree->mod->ns+12] =    2 ;
      tree->step_mat[15*tree->mod->ns+13] =    2 ;
      tree->step_mat[15*tree->mod->ns+14] =    2 ;
      tree->step_mat[15*tree->mod->ns+15] =    0 ;
      tree->step_mat[15*tree->mod->ns+16] =    2 ;
      tree->step_mat[15*tree->mod->ns+17] =    2 ;
      tree->step_mat[15*tree->mod->ns+18] =    2 ;
      tree->step_mat[15*tree->mod->ns+19] =    3 ;
      tree->step_mat[16*tree->mod->ns+ 0] =    2 ;
      tree->step_mat[16*tree->mod->ns+ 1] =    2 ;
      tree->step_mat[16*tree->mod->ns+ 2] =    2 ;
      tree->step_mat[16*tree->mod->ns+ 3] =    3 ;
      tree->step_mat[16*tree->mod->ns+ 4] =    3 ;
      tree->step_mat[16*tree->mod->ns+ 5] =    3 ;
      tree->step_mat[16*tree->mod->ns+ 6] =    3 ;
      tree->step_mat[16*tree->mod->ns+ 7] =    3 ;
      tree->step_mat[16*tree->mod->ns+ 8] =    3 ;
      tree->step_mat[16*tree->mod->ns+ 9] =    2 ;
      tree->step_mat[16*tree->mod->ns+10] =    3 ;
      tree->step_mat[16*tree->mod->ns+11] =    2 ;
      tree->step_mat[16*tree->mod->ns+12] =    2 ;
      tree->step_mat[16*tree->mod->ns+13] =    3 ;
      tree->step_mat[16*tree->mod->ns+14] =    2 ;
      tree->step_mat[16*tree->mod->ns+15] =    2 ;
      tree->step_mat[16*tree->mod->ns+16] =    0 ;
      tree->step_mat[16*tree->mod->ns+17] =    3 ;
      tree->step_mat[16*tree->mod->ns+18] =    3 ;
      tree->step_mat[16*tree->mod->ns+19] =    3 ;
      tree->step_mat[17*tree->mod->ns+ 0] =    3 ;
      tree->step_mat[17*tree->mod->ns+ 1] =    2 ;
      tree->step_mat[17*tree->mod->ns+ 2] =    3 ;
      tree->step_mat[17*tree->mod->ns+ 3] =    3 ;
      tree->step_mat[17*tree->mod->ns+ 4] =    1 ;
      tree->step_mat[17*tree->mod->ns+ 5] =    2 ;
      tree->step_mat[17*tree->mod->ns+ 6] =    2 ;
      tree->step_mat[17*tree->mod->ns+ 7] =    2 ;
      tree->step_mat[17*tree->mod->ns+ 8] =    3 ;
      tree->step_mat[17*tree->mod->ns+ 9] =    3 ;
      tree->step_mat[17*tree->mod->ns+10] =    2 ;
      tree->step_mat[17*tree->mod->ns+11] =    2 ;
      tree->step_mat[17*tree->mod->ns+12] =    2 ;
      tree->step_mat[17*tree->mod->ns+13] =    2 ;
      tree->step_mat[17*tree->mod->ns+14] =    3 ;
      tree->step_mat[17*tree->mod->ns+15] =    2 ;
      tree->step_mat[17*tree->mod->ns+16] =    3 ;
      tree->step_mat[17*tree->mod->ns+17] =    0 ;
      tree->step_mat[17*tree->mod->ns+18] =    2 ;
      tree->step_mat[17*tree->mod->ns+19] =    3 ;
      tree->step_mat[18*tree->mod->ns+ 0] =    3 ;
      tree->step_mat[18*tree->mod->ns+ 1] =    3 ;
      tree->step_mat[18*tree->mod->ns+ 2] =    2 ;
      tree->step_mat[18*tree->mod->ns+ 3] =    2 ;
      tree->step_mat[18*tree->mod->ns+ 4] =    2 ;
      tree->step_mat[18*tree->mod->ns+ 5] =    2 ;
      tree->step_mat[18*tree->mod->ns+ 6] =    2 ;
      tree->step_mat[18*tree->mod->ns+ 7] =    3 ;
      tree->step_mat[18*tree->mod->ns+ 8] =    2 ;
      tree->step_mat[18*tree->mod->ns+ 9] =    3 ;
      tree->step_mat[18*tree->mod->ns+10] =    2 ;
      tree->step_mat[18*tree->mod->ns+11] =    2 ;
      tree->step_mat[18*tree->mod->ns+12] =    3 ;
      tree->step_mat[18*tree->mod->ns+13] =    2 ;
      tree->step_mat[18*tree->mod->ns+14] =    3 ;
      tree->step_mat[18*tree->mod->ns+15] =    2 ;
      tree->step_mat[18*tree->mod->ns+16] =    3 ;
      tree->step_mat[18*tree->mod->ns+17] =    2 ;
      tree->step_mat[18*tree->mod->ns+18] =    0 ;
      tree->step_mat[18*tree->mod->ns+19] =    3 ;
      tree->step_mat[19*tree->mod->ns+ 0] =    2 ;
      tree->step_mat[19*tree->mod->ns+ 1] =    3 ;
      tree->step_mat[19*tree->mod->ns+ 2] =    3 ;
      tree->step_mat[19*tree->mod->ns+ 3] =    2 ;
      tree->step_mat[19*tree->mod->ns+ 4] =    3 ;
      tree->step_mat[19*tree->mod->ns+ 5] =    3 ;
      tree->step_mat[19*tree->mod->ns+ 6] =    2 ;
      tree->step_mat[19*tree->mod->ns+ 7] =    2 ;
      tree->step_mat[19*tree->mod->ns+ 8] =    3 ;
      tree->step_mat[19*tree->mod->ns+ 9] =    2 ;
      tree->step_mat[19*tree->mod->ns+10] =    2 ;
      tree->step_mat[19*tree->mod->ns+11] =    3 ;
      tree->step_mat[19*tree->mod->ns+12] =    2 ;
      tree->step_mat[19*tree->mod->ns+13] =    2 ;
      tree->step_mat[19*tree->mod->ns+14] =    3 ;
      tree->step_mat[19*tree->mod->ns+15] =    3 ;
      tree->step_mat[19*tree->mod->ns+16] =    3 ;
      tree->step_mat[19*tree->mod->ns+17] =    3 ;
      tree->step_mat[19*tree->mod->ns+18] =    3 ;
      tree->step_mat[19*tree->mod->ns+19] =    0 ;
    }
  else
    {
      tree->step_mat[0*tree->mod->ns+0] = 0;
      tree->step_mat[0*tree->mod->ns+1] = 2;
      tree->step_mat[0*tree->mod->ns+2] = 1;
      tree->step_mat[0*tree->mod->ns+3] = 2;

      tree->step_mat[1*tree->mod->ns+0] = 2;
      tree->step_mat[1*tree->mod->ns+1] = 0;
      tree->step_mat[1*tree->mod->ns+2] = 2;
      tree->step_mat[1*tree->mod->ns+3] = 1;

      tree->step_mat[2*tree->mod->ns+0] = 1;
      tree->step_mat[2*tree->mod->ns+1] = 2;
      tree->step_mat[2*tree->mod->ns+2] = 0;
      tree->step_mat[2*tree->mod->ns+3] = 2;

      tree->step_mat[3*tree->mod->ns+0] = 2;
      tree->step_mat[3*tree->mod->ns+1] = 1;
      tree->step_mat[3*tree->mod->ns+2] = 2;
      tree->step_mat[3*tree->mod->ns+3] = 0;

/*       tree->step_mat[0*tree->mod->ns+0] = 0; */
/*       tree->step_mat[0*tree->mod->ns+1] = 1; */
/*       tree->step_mat[0*tree->mod->ns+2] = 1; */
/*       tree->step_mat[0*tree->mod->ns+3] = 1; */

/*       tree->step_mat[1*tree->mod->ns+0] = 1; */
/*       tree->step_mat[1*tree->mod->ns+1] = 0; */
/*       tree->step_mat[1*tree->mod->ns+2] = 1; */
/*       tree->step_mat[1*tree->mod->ns+3] = 1; */

/*       tree->step_mat[2*tree->mod->ns+0] = 1; */
/*       tree->step_mat[2*tree->mod->ns+1] = 1; */
/*       tree->step_mat[2*tree->mod->ns+2] = 0; */
/*       tree->step_mat[2*tree->mod->ns+3] = 1; */

/*       tree->step_mat[3*tree->mod->ns+0] = 1; */
/*       tree->step_mat[3*tree->mod->ns+1] = 1; */
/*       tree->step_mat[3*tree->mod->ns+2] = 1; */
/*       tree->step_mat[3*tree->mod->ns+3] = 0; */

    }

  for(i=0;i<tree->mod->ns;i++) tree->step_mat[i*tree->mod->ns+i] = 0;
}

/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/

// Tree should be ready for likelihood analysis when calling
// this function.
void Stepwise_Add_Pars(t_tree *tree)
{
  t_edge **residuals,**targets,*best_target;
  int *nd_idx,i,j,n_targets,*tg_idx;

  residuals   = (t_edge **)mCalloc(tree->n_otu-3,sizeof(t_edge *));
  targets     = (t_edge **)mCalloc(2*tree->n_otu-3,sizeof(t_edge *));
  best_target = NULL;
  nd_idx      = Permutate(tree->n_otu-3);

  // Remove all tips except that corresponding to a_nodes[0], 
  // a_nodes[1] and a_nodes[2].  
  for(i=0;i<tree->n_otu-3;++i)
    {
      Prune_Subtree(tree->a_nodes[i+3]->v[0],                   
                    tree->a_nodes[i+3],
                    NULL,
                    residuals+i,
                    tree);
    }

  // Initial targets
  n_targets = 3;
  for(i=0;i<n_targets;i++) targets[i] = tree->a_nodes[i]->b[0];

  // Regraft each tip on the tree at most parsimonious position
  for(i=0;i<tree->n_otu-3;++i)
    {
      Set_Both_Sides(YES,tree);
      Pars(NULL,tree);

      /* printf("\n. [%d/%d]",i,tree->n_otu-3); */

      tree->best_pars = 1E+8;
      best_target = NULL;
      tg_idx      = Permutate(n_targets);

      for(j=0;j<n_targets;j++)
        {
          Graft_Subtree(targets[tg_idx[j]],
                        tree->a_nodes[nd_idx[i]+3]->v[0],
                        NULL,
                        residuals[i],
                        NULL,
                        tree);
          
          Update_Partial_Pars(tree,
                        tree->a_nodes[nd_idx[i]+3]->b[0],
                        tree->a_nodes[nd_idx[i]+3]->v[0]);
          Pars(tree->a_nodes[nd_idx[i]+3]->b[0],tree);
          

          if(tree->c_pars < tree->best_pars)
            {
              tree->best_pars = tree->c_pars;
              best_target = targets[tg_idx[j]];
            }
          
          Prune_Subtree(tree->a_nodes[nd_idx[i]+3]->v[0],                        
                        tree->a_nodes[nd_idx[i]+3],
                        NULL,
                        residuals+i,
                        tree);
        }

      assert(best_target);
            
      Graft_Subtree(best_target,
                    tree->a_nodes[nd_idx[i]+3]->v[0],
                    NULL,
                    residuals[i],
                    NULL,
                    tree);
      
      targets[n_targets]   = residuals[i]; 
      targets[n_targets+1] = tree->a_nodes[nd_idx[i]+3]->b[0];
      
      Free(tg_idx);
      n_targets+=2;
    }

  /* printf("\n. pars: %d",tree->c_pars); */
  /* Exit("\n"); */

  Free(nd_idx);
  Free(residuals);
  Free(targets);
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

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
