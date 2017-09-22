/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include "assert.h"
#include "lk.h"
#ifdef BEAGLE
#include "beagle_utils.h"
#endif
#include <stdint.h>

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Init_Tips_At_One_Site_Nucleotides_Float(char state, int pos, phydbl *p_lk)
{
  switch(state)
    {
    case 'A' : p_lk[pos+0]=1.; p_lk[pos+1]=p_lk[pos+2]=p_lk[pos+3]=.0;
      break;
    case 'C' : p_lk[pos+1]=1.; p_lk[pos+0]=p_lk[pos+2]=p_lk[pos+3]=.0;
      break;
    case 'G' : p_lk[pos+2]=1.; p_lk[pos+1]=p_lk[pos+0]=p_lk[pos+3]=.0;
      break;
    case 'T' : p_lk[pos+3]=1.; p_lk[pos+1]=p_lk[pos+2]=p_lk[pos+0]=.0;
      break;
    case 'U' : p_lk[pos+3]=1.; p_lk[pos+1]=p_lk[pos+2]=p_lk[pos+0]=.0;
      break;
    case 'M' : p_lk[pos+0]=p_lk[pos+1]=1.; p_lk[pos+2]=p_lk[pos+3]=.0;
      break;
    case 'R' : p_lk[pos+0]=p_lk[pos+2]=1.; p_lk[pos+1]=p_lk[pos+3]=.0;
      break;
    case 'W' : p_lk[pos+0]=p_lk[pos+3]=1.; p_lk[pos+1]=p_lk[pos+2]=.0;
      break;
    case 'S' : p_lk[pos+1]=p_lk[pos+2]=1.; p_lk[pos+0]=p_lk[pos+3]=.0;
      break;
    case 'Y' : p_lk[pos+1]=p_lk[pos+3]=1.; p_lk[pos+0]=p_lk[pos+2]=.0;
      break;
    case 'K' : p_lk[pos+2]=p_lk[pos+3]=1.; p_lk[pos+0]=p_lk[pos+1]=.0;
      break;
    case 'B' : p_lk[pos+1]=p_lk[pos+2]=p_lk[pos+3]=1.; p_lk[pos+0]=.0;
      break;
    case 'D' : p_lk[pos+0]=p_lk[pos+2]=p_lk[pos+3]=1.; p_lk[pos+1]=.0;
      break;
    case 'H' : p_lk[pos+0]=p_lk[pos+1]=p_lk[pos+3]=1.; p_lk[pos+2]=.0;
      break;
    case 'V' : p_lk[pos+0]=p_lk[pos+1]=p_lk[pos+2]=1.; p_lk[pos+3]=.0;
      break;
    case 'N' : case 'X' : case '?' : case 'O' : case '-' :
      p_lk[pos+0]=p_lk[pos+1]=p_lk[pos+2]=p_lk[pos+3]=1.;break;
    default :
      {
        PhyML_Fprintf(stderr,"\n. Unknown character state : %c\n",state);
        Exit("\n. Init failed (data type supposed to be DNA)\n");
        break;
      }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Init_Tips_At_One_Site_Nucleotides_Int(char state, int pos, short int *p_pars)
{
  switch(state)
    {
    case 'A' : p_pars[pos+0]=1; p_pars[pos+1]=p_pars[pos+2]=p_pars[pos+3]=0;
      break;
    case 'C' : p_pars[pos+1]=1; p_pars[pos+0]=p_pars[pos+2]=p_pars[pos+3]=0;
      break;
    case 'G' : p_pars[pos+2]=1; p_pars[pos+1]=p_pars[pos+0]=p_pars[pos+3]=0;
      break;
    case 'T' : p_pars[pos+3]=1; p_pars[pos+1]=p_pars[pos+2]=p_pars[pos+0]=0;
      break;
    case 'U' : p_pars[pos+3]=1; p_pars[pos+1]=p_pars[pos+2]=p_pars[pos+0]=0;
      break;
    case 'M' : p_pars[pos+0]=p_pars[pos+1]=1; p_pars[pos+2]=p_pars[pos+3]=0;
      break;
    case 'R' : p_pars[pos+0]=p_pars[pos+2]=1; p_pars[pos+1]=p_pars[pos+3]=0;
      break;
    case 'W' : p_pars[pos+0]=p_pars[pos+3]=1; p_pars[pos+1]=p_pars[pos+2]=0;
      break;
    case 'S' : p_pars[pos+1]=p_pars[pos+2]=1; p_pars[pos+0]=p_pars[pos+3]=0;
      break;
    case 'Y' : p_pars[pos+1]=p_pars[pos+3]=1; p_pars[pos+0]=p_pars[pos+2]=0;
      break;
    case 'K' : p_pars[pos+2]=p_pars[pos+3]=1; p_pars[pos+0]=p_pars[pos+1]=0;
      break;
    case 'B' : p_pars[pos+1]=p_pars[pos+2]=p_pars[pos+3]=1; p_pars[pos+0]=0;
      break;
    case 'D' : p_pars[pos+0]=p_pars[pos+2]=p_pars[pos+3]=1; p_pars[pos+1]=0;
      break;
    case 'H' : p_pars[pos+0]=p_pars[pos+1]=p_pars[pos+3]=1; p_pars[pos+2]=0;
      break;
    case 'V' : p_pars[pos+0]=p_pars[pos+1]=p_pars[pos+2]=1; p_pars[pos+3]=0;
      break;
    case 'N' : case 'X' : case '?' : case 'O' : case '-' :
      p_pars[pos+0]=p_pars[pos+1]=p_pars[pos+2]=p_pars[pos+3]=1;break;
    default :
      {
        PhyML_Fprintf(stderr,"\n. Unknown character state : %c\n",state);
        Exit("\n. Init failed (data type supposed to be DNA)\n");
        break;
      }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Init_Tips_At_One_Site_AA_Float(char aa, int pos, phydbl *p_lk)
{
  int i;

  for(i=0;i<20;i++) p_lk[pos+i] = .0;

  switch(aa){
  case 'A' : p_lk[pos+0]= 1.; break;/* Alanine */
  case 'R' : p_lk[pos+1]= 1.; break;/* Arginine */
  case 'N' : p_lk[pos+2]= 1.; break;/* Asparagine */
  case 'D' : p_lk[pos+3]= 1.; break;/* Aspartic acid */
  case 'C' : p_lk[pos+4]= 1.; break;/* Cysteine */
  case 'Q' : p_lk[pos+5]= 1.; break;/* Glutamine */
  case 'E' : p_lk[pos+6]= 1.; break;/* Glutamic acid */
  case 'G' : p_lk[pos+7]= 1.; break;/* Glycine */
  case 'H' : p_lk[pos+8]= 1.; break;/* Histidine */
  case 'I' : p_lk[pos+9]= 1.; break;/* Isoleucine */
  case 'L' : p_lk[pos+10]=1.; break;/* Leucine */
  case 'K' : p_lk[pos+11]=1.; break;/* Lysine */
  case 'M' : p_lk[pos+12]=1.; break;/* Methionine */
  case 'F' : p_lk[pos+13]=1.; break;/* Phenylalanin */
  case 'P' : p_lk[pos+14]=1.; break;/* Proline */
  case 'S' : p_lk[pos+15]=1.; break;/* Serine */
  case 'T' : p_lk[pos+16]=1.; break;/* Threonine */
  case 'W' : p_lk[pos+17]=1.; break;/* Tryptophan */
  case 'Y' : p_lk[pos+18]=1.; break;/* Tyrosine */
  case 'V' : p_lk[pos+19]=1.; break;/* Valine */

  case 'B' : p_lk[pos+2]= 1.; break;/* Asparagine */
  case 'Z' : p_lk[pos+5]= 1.; break;/* Glutamine */

  case 'X' : case '?' : case '-' : for(i=0;i<20;i++) p_lk[pos+i] = 1.; break;
  default :
    {
      PhyML_Fprintf(stderr,"\n. Unknown character state : %c\n",aa);
      Exit("\n. Init failed (data type supposed to be amino-acids)\n");
      break;
    }
  }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Init_Tips_At_One_Site_AA_Int(char aa, int pos, short int *p_pars)
{
  int i;

  for(i=0;i<20;i++) p_pars[pos+i] = .0;

  switch(aa){
  case 'A' : p_pars[pos+0]  = 1; break;/* Alanine */
  case 'R' : p_pars[pos+1]  = 1; break;/* Arginine */
  case 'N' : p_pars[pos+2]  = 1; break;/* Asparagine */
  case 'D' : p_pars[pos+3]  = 1; break;/* Aspartic acid */
  case 'C' : p_pars[pos+4]  = 1; break;/* Cysteine */
  case 'Q' : p_pars[pos+5]  = 1; break;/* Glutamine */
  case 'E' : p_pars[pos+6]  = 1; break;/* Glutamic acid */
  case 'G' : p_pars[pos+7]  = 1; break;/* Glycine */
  case 'H' : p_pars[pos+8]  = 1; break;/* Histidine */
  case 'I' : p_pars[pos+9]  = 1; break;/* Isoleucine */
  case 'L' : p_pars[pos+10] = 1; break;/* Leucine */
  case 'K' : p_pars[pos+11] = 1; break;/* Lysine */
  case 'M' : p_pars[pos+12] = 1; break;/* Methionine */
  case 'F' : p_pars[pos+13] = 1; break;/* Phenylalanin */
  case 'P' : p_pars[pos+14] = 1; break;/* Proline */
  case 'S' : p_pars[pos+15] = 1; break;/* Serine */
  case 'T' : p_pars[pos+16] = 1; break;/* Threonine */
  case 'W' : p_pars[pos+17] = 1; break;/* Tryptophan */
  case 'Y' : p_pars[pos+18] = 1; break;/* Tyrosine */
  case 'V' : p_pars[pos+19] = 1; break;/* Valine */

  case 'B' : p_pars[pos+2]  = 1; break;/* Asparagine */
  case 'Z' : p_pars[pos+5]  = 1; break;/* Glutamine */

  case 'X' : case '?' : case '-' : for(i=0;i<20;i++) p_pars[pos+i] = 1; break;
  default :
    {
      PhyML_Fprintf(stderr,"\n. Unknown character state : %c\n",aa);
      Exit("\n. Init failed (data type supposed to be amino-acids)\n");
      break;
    }
  }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Init_Tips_At_One_Site_Generic_Float(char *state, int ns, int state_len, int pos, phydbl *p_lk)
{
  int i;
  int state_int;

  for(i=0;i<ns;i++) p_lk[pos+i] = 0.;

  if(Is_Ambigu(state,GENERIC,state_len)) for(i=0;i<ns;i++) p_lk[pos+i] = 1.;
  else
    {
      char format[6];
      sprintf(format,"%%%dd",state_len);
      if(!sscanf(state,format,&state_int))
    {
      PhyML_Fprintf(stderr,"\n. state='%c'",state);
      PhyML_Fprintf(stderr,"\n. Err in file %s at line %d (function '%s')\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
      if(state_int > ns)
    {
      PhyML_Fprintf(stderr,"\n. %s %d cstate: %.2s istate: %d state_len: %d.\n",__FILE__,__LINE__,state,state_int,state_len);
      PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
      Warn_And_Exit("");
    }
      p_lk[pos+state_int] = 1.;
      /*       PhyML_Printf("\n. %s %d cstate: %.2s istate: %d state_len: %d ns: %d pos: %d",__FILE__,__LINE__,state,state_int,state_len,ns,pos); */
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Init_Tips_At_One_Site_Generic_Int(char *state, int ns, int state_len, int pos, short int *p_pars)
{
  int i;
  int state_int;

  for(i=0;i<ns;i++) p_pars[pos+i] = 0;

  if(Is_Ambigu(state,GENERIC,state_len)) for(i=0;i<ns;i++) p_pars[pos+i] = 1;
  else
    {
      char format[6];
      sprintf(format,"%%%dd",state_len);
      if(!sscanf(state,format,&state_int))
        {
          PhyML_Fprintf(stderr,"\n. state='%c'",state);
          PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
          Warn_And_Exit("");
        }
      if(state_int > ns)
        {
          PhyML_Fprintf(stderr,"\n. %s %d cstate: %.2s istate: %d state_len: %d.\n",__FILE__,__LINE__,state,state_int,state_len);
          PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
          Warn_And_Exit("");
        }
      p_pars[pos+state_int] = 1;
/*       PhyML_Printf("\n* %s %d cstate: %.2s istate: %d state_len: %d ns: %d pos: %d",__FILE__,__LINE__,state,state_int,state_len,ns,pos); */
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Get_All_Partial_Lk_Scale(t_tree *tree, t_edge *b_fcus, t_node *a, t_node *d)
{
  Update_Partial_Lk(tree,b_fcus,d);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Post_Order_Lk(t_node *a, t_node *d, t_tree *tree)
{
  int i,dir;

  dir = -1;
  
  if(d->tax) return;
  else
    {
      if(tree->is_mixt_tree)
        {
          MIXT_Post_Order_Lk(a,d,tree);
          return;
        }

      if(tree->n_root)
        {
          for(i=0;i<3;i++)
            {
              if(d->v[i] != a && d->b[i] != tree->e_root)
                Post_Order_Lk(d,d->v[i],tree);
              else dir = i;
            }
        }
      else
        {
          for(i=0;i<3;i++)
            {
              if(d->v[i] != a)
                Post_Order_Lk(d,d->v[i],tree);
              else dir = i;
            }
        }

      if(dir < 0)
        {
          PhyML_Printf("\n. a: %d d: %d [%d %d %d] [%d %d %d]",
                       a->num,
                       d->num,
                       d->v[0] ? d->v[0]->num : -1,
                       d->v[1] ? d->v[1]->num : -1,
                       d->v[2] ? d->v[2]->num : -1,
                       d->b[0] ? d->b[0]->num : -1,
                       d->b[1] ? d->b[1]->num : -1,
                       d->b[2] ? d->b[2]->num : -1);
          fflush(NULL);
        }
      assert(dir > -1);

      if(tree->ignore_root == NO && d->b[dir] == tree->e_root)
        {
          if(d == tree->n_root->v[1]) Update_Partial_Lk(tree,tree->n_root->b[1],d);
          else                        Update_Partial_Lk(tree,tree->n_root->b[2],d);
        }
      else
        {
          Update_Partial_Lk(tree,d->b[dir],d);
        }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Pre_Order_Lk(t_node *a, t_node *d, t_tree *tree)
{
  int i;

  if(d->tax) return;
  else
    {
      if(tree->is_mixt_tree)
        {
          MIXT_Pre_Order_Lk(a,d,tree);
          return;
        }

      if(tree->n_root)
        {
          for(i=0;i<3;++i)
            {
              if(d->v[i] != a && d->b[i] != tree->e_root)
                {
                  Update_Partial_Lk(tree,d->b[i],d);
                  Pre_Order_Lk(d,d->v[i],tree);
                }
            }
        }
      else
        {
          for(i=0;i<3;++i)
            {
              if(d->v[i] != a)
                {
                  Update_Partial_Lk(tree,d->b[i],d);
                  Pre_Order_Lk(d,d->v[i],tree);
                }
            }
        }
    }      
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

// Updates all partial likelihood vectors. Depending on whether 
// both_sides = YES or NO, only 'up' or 'up'&'down' partials will
// be updates
void Update_All_Partial_Lk(t_tree *tree)
{
  if(tree->n_root)
    {
      if(tree->ignore_root == NO)
        {          
          Post_Order_Lk(tree->n_root,tree->n_root->v[1],tree);
          Post_Order_Lk(tree->n_root,tree->n_root->v[2],tree);
          
          Update_Partial_Lk(tree,tree->n_root->b[1],tree->n_root);
          Update_Partial_Lk(tree,tree->n_root->b[2],tree->n_root);
          
          if(tree->both_sides == YES)
            {
              Pre_Order_Lk(tree->n_root,tree->n_root->v[2],tree);
              Pre_Order_Lk(tree->n_root,tree->n_root->v[1],tree);
            }
        }
      else
        {

          Post_Order_Lk(tree->e_root->rght,tree->e_root->left,tree);
          Post_Order_Lk(tree->e_root->left,tree->e_root->rght,tree);
          
          if(tree->both_sides == YES)
            {
              Pre_Order_Lk(tree->e_root->rght,tree->e_root->left,tree);
              Pre_Order_Lk(tree->e_root->left,tree->e_root->rght,tree);
            }
        }
    }
  else
    {
      Post_Order_Lk(tree->a_nodes[0],tree->a_nodes[0]->v[0],tree);
      if(tree->both_sides == YES)
        Pre_Order_Lk(tree->a_nodes[0],tree->a_nodes[0]->v[0],tree);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl Lk(t_edge *b, t_tree *tree)
{
  unsigned int br,catg,state,ambiguity_check,site;
  phydbl len,*expl,*dot_prod;  
  const unsigned int ns = tree->mod->ns;
  const unsigned int ncatg = tree->mod->ras->n_catg;
  const unsigned int npatterns = tree->n_pattern;

  
  if(tree->eval_alnL == NO) return UNLIKELY;

  if(b == NULL && tree->mod->s_opt->curr_opt_free_rates == YES)
    {
      tree->mod->s_opt->curr_opt_free_rates = NO;
      Optimize_Free_Rate_Weights(tree,YES,YES);
      tree->mod->s_opt->curr_opt_free_rates = YES;
    }
  
  
  if(tree->is_mixt_tree == YES) 
    {
#ifdef BEAGLE
      Warn_And_Exit(TODO_BEAGLE);
#endif
      return MIXT_Lk(b,tree);
    }
  
  tree->old_lnL = tree->c_lnL;
  
#ifdef PHYREX
  PHYREX_Ldsk_To_Tree(tree);
#endif


#if (defined PHYTIME || defined INVITEE || defined PHYREX)
  if((tree->rates) && (tree->rates->bl_from_rt)) RATES_Update_Cur_Bl(tree);
#endif

if(tree->rates && tree->io->lk_approx == NORMAL)
    {
#ifdef BEAGLE
      Warn_And_Exit(TODO_BEAGLE);
#endif
      tree->c_lnL = Lk_Normal_Approx(tree);
      return tree->c_lnL;
    }
  
 expl     = tree->expl;
 dot_prod = tree->dot_prod;
 
 
 if(b == NULL)
   {
     Update_Boundaries(tree->mod);
     Update_RAS(tree->mod);
     Update_Efrq(tree->mod);
     Update_Eigen(tree->mod);
   }
  
 if(tree->mod->s_opt->skip_tree_traversal == NO)
    {
      if(!b) //Update PMat for all edges
        {
          for(br=0;br<2*tree->n_otu-3;++br) Update_PMat_At_Given_Edge(tree->a_edges[br],tree);
          if(tree->n_root && tree->ignore_root == NO)
            {
              Update_PMat_At_Given_Edge(tree->n_root->b[1],tree);
              Update_PMat_At_Given_Edge(tree->n_root->b[2],tree);
            }
        }
      else//Update PMat for a specific edge
        {
          if(tree->use_eigen_lr == NO) Update_PMat_At_Given_Edge(b,tree);
        }
      
      if(!b)
        {
          if(tree->n_root != NULL)
            {
              if(tree->ignore_root == NO)
                {
                  Post_Order_Lk(tree->n_root,tree->n_root->v[1],tree);
                  Post_Order_Lk(tree->n_root,tree->n_root->v[2],tree);
                  
                  Update_Partial_Lk(tree,tree->n_root->b[1],tree->n_root);
                  Update_Partial_Lk(tree,tree->n_root->b[2],tree->n_root);
                  
                  if(tree->both_sides == YES)
                    {
                      Pre_Order_Lk(tree->n_root,tree->n_root->v[2],tree);
                      Pre_Order_Lk(tree->n_root,tree->n_root->v[1],tree);
                    }
                }
              else
                {
                  Post_Order_Lk(tree->e_root->rght,tree->e_root->left,tree);
                  Post_Order_Lk(tree->e_root->left,tree->e_root->rght,tree);
                  
                  if(tree->both_sides == YES)
                    {
                      Pre_Order_Lk(tree->e_root->rght,tree->e_root->left,tree);
                      Pre_Order_Lk(tree->e_root->left,tree->e_root->rght,tree);
                    }
                }
            }
          else
            {
              Post_Order_Lk(tree->a_nodes[0],tree->a_nodes[0]->v[0],tree);
              if(tree->both_sides == YES)
                Pre_Order_Lk(tree->a_nodes[0],tree->a_nodes[0]->v[0],tree);
            }
        }
    }

  if(!b)
    {
      if(tree->n_root) 
        {
          if(tree->ignore_root == NO)
            b = (tree->n_root->v[1]->tax == NO)?(tree->n_root->b[2]):(tree->n_root->b[1]);
          else
            b = tree->e_root;
        }
      else                                        
        b = tree->a_nodes[0]->b[0];
    }

  tree->c_lnL             = .0;
  tree->sum_min_sum_scale = .0;

#ifdef BEAGLE
  calc_edgelks_beagle(b, tree);
#else

  if(tree->update_eigen_lr == YES) Update_Eigen_Lr(b,tree);
  
  if(tree->use_eigen_lr == YES)
    {  
      for(catg=0;catg<ncatg;++catg)
        {
          len = MAX(0.0,b->l->v)*tree->mod->ras->gamma_rr->v[catg];
          len *= tree->mod->br_len_mult->v;
          if(tree->mixt_tree != NULL)     len *= tree->mixt_tree->mod->ras->gamma_rr->v[tree->mod->ras->parent_class_number];
          if(len < tree->mod->l_min)      len = tree->mod->l_min;
          else if(len > tree->mod->l_max) len = tree->mod->l_max;
          for(state=0;state<ns;++state) expl[catg*ns+state] = (phydbl)POW(tree->mod->eigen->e_val[state],len);
        }
    }


  for(site=0;site<npatterns;++site)
    {
      ambiguity_check = -1;
      state           = -1;
      tree->curr_site = site;
      
      if(tree->data->wght[tree->curr_site] > SMALL)
        {
          if((b->rght->tax) && (tree->mod->s_opt->greedy == NO))
            {
              ambiguity_check = b->rght->c_seq->is_ambigu[tree->curr_site];
              if(ambiguity_check == NO)
                {
                  state = b->rght->c_seq->d_state[tree->curr_site];
                }
            }

          if(tree->mod->use_m4mod) ambiguity_check = YES;
          
          if(tree->use_eigen_lr == YES)
            {
              Lk_Core_Eigen_Lr(expl,dot_prod + site*ns*ncatg,NO,b,tree);
            }
          else
            {
              if(b->rght->tax == YES)
                Lk_Core(state,ambiguity_check,NO,
                        b->p_lk_left + site*ns*ncatg,
                        b->p_lk_tip_r + site*ns,
                        b->Pij_rr,b,tree);
              else
                Lk_Core(state,ambiguity_check,NO,
                        b->p_lk_left + site*ns*ncatg,
                        b->p_lk_rght  + site*ns*ncatg,
                        b->Pij_rr,b,tree);                
            }
        }
    }
#endif

  return tree->c_lnL;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// First and second derivative of the log-likelihood with respect
// to the length of edge b
phydbl dLk(phydbl *l, t_edge *b, t_tree *tree)
{
  unsigned int catg,state,site;
  phydbl len,rr;
  phydbl lk,dlk,d2lk,dlnlk,d2lnlk;
  phydbl ev,expevlen;
   
  const unsigned int ns = tree->mod->ns;
  const unsigned int ncatg = tree->mod->ras->n_catg;
  const unsigned int npattern = tree->n_pattern;

  phydbl *dot_prod = tree->dot_prod;
  phydbl *expl = tree->expl;
  phydbl *expld = expl + 1*ncatg*ns;
  phydbl *expld2 = expl + 2*ncatg*ns;
  
  assert(isnan(*l) == FALSE);

  if(*l < tree->mod->l_min)      *l = tree->mod->l_min;
  else if(*l > tree->mod->l_max) *l = tree->mod->l_max;      


  assert(b != NULL);
  
  if(tree->is_mixt_tree == YES)
    {
#ifdef BEAGLE
      Warn_And_Exit(TODO_BEAGLE);
#endif
      return MIXT_dLk(l,b,tree);
    }
    
  if(tree->update_eigen_lr == YES) Update_Eigen_Lr(b,tree);
  
  for(catg=0;catg<ncatg;catg++)
    {
      rr = tree->mod->ras->gamma_rr->v[catg];
      rr *=  tree->mod->br_len_mult->v;      

      len = (*l) * rr;
            
      if(isinf(len) || isnan(len)) 
        {
          PhyML_Fprintf(stderr,"\n. len=%f rr=%f l=%f",len,rr,*l);
          Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
        }

      if(len < tree->mod->l_min)      len = tree->mod->l_min;
      else if(len > tree->mod->l_max) len = tree->mod->l_max;      
      // value of rr should be corrected too if any of these two conditions
      // is true. Leads to numerical precision issues though...
      
      
      for(state=0;state<ns;++state) 
        {
          ev = log(tree->mod->eigen->e_val[state]);
          expevlen = exp(ev*len);

          expl[0*ncatg*ns + catg*ns + state] = expevlen;
          expl[1*ncatg*ns + catg*ns + state] = expevlen*ev*rr;
          expl[2*ncatg*ns + catg*ns + state] = expevlen*ev*ev*rr*rr;
        }
    }
  

  dlnlk       = 0.0;
  d2lnlk      = 0.0;
  tree->c_lnL = .0;
  
  for(site=0;site<npattern;++site)
    {
      tree->curr_site = site;
      
      dlk   = Lk_Core_Eigen_Lr(expld ,dot_prod + site*ns*ncatg,YES,b,tree);
      d2lk  = Lk_Core_Eigen_Lr(expld2,dot_prod + site*ns*ncatg,YES,b,tree);
      lk    = Lk_Core_Eigen_Lr(expl  ,dot_prod + site*ns*ncatg,NO ,b,tree);

      /* printf("\n. site: %4d l:%15G lk:%15G dlk:%15G d2lk:%15G",site,b->l->v,lk,dlk,d2lk); */
      
      if(!(lk > 0.0))
        {
          /* PhyML_Printf("\n. dlk: %G d2lk: %G lk: %G",dlk,d2lk,lk); */
          int l;
          int i,j;
          phydbl *U,*V,*Q;
          printf("\n. l: %G dlk: %G d2lk: %G lk: %G",b->l->v,dlk,d2lk,lk);
          lk = 0.0;
          for(l=0;l<ns;++l)
            {
              lk += dot_prod[site*ns*ncatg+l] * expl[l];
              printf("\n. lk: %12G dot_prod: %12G expl: %12G left: %12G rght: %12G",
                     lk,dot_prod[site*ns*ncatg+l],expl[l],
                     b->p_lk_left[site*ns*ncatg + catg*ns + l],
                     b->rght->tax ? b->p_lk_tip_r[site*ns + l] : b->p_lk_rght[site*ns*ncatg + catg*ns + l]);
            }
          printf("\n. U\n");
          U = tree->mod->eigen->r_e_vect;
          for(i=0;i<ns;i++) { for(j=0;j<ns;j++) printf("%12G ",U[i*ns+j]); printf("\n"); }
          printf("\n");
          printf("\n. V\n");
          V = tree->mod->eigen->l_e_vect;
          for(i=0;i<ns;i++) { for(j=0;j<ns;j++) printf("%12G ",V[i*ns+j]); printf("\n"); }
          printf("\n");
          printf("\n. exp(D)\n");
          for(i=0;i<ns;i++) printf("%12G\n",tree->mod->eigen->e_val[i]);
          printf("\n");
          printf("\n. Q\n");
          Q = tree->mod->r_mat->qmat->v;
          for(i=0;i<ns;i++) { for(j=0;j<ns;j++) printf("%12G ",Q[i*ns+j]); printf("\n"); }
          printf("\n");
          printf("\n. F\n");
          for(i=0;i<ns;i++) printf("f(%3d)=%12G\n",i,tree->mod->e_frq->pi->v[i]); 
          printf("\n\n. tstv : %G",tree->mod->kappa->v);
          
          for(catg=0;catg<ncatg;++catg)
            {
              for(state=0;state<ns;++state) 
                {
                  printf("\n. catg: %3d state: %3d -- expl: %12G expld: %12G expld2: %12G dot_prod: %12G",
                         catg,
                         state,
                         expl[0*ncatg*ns + catg*ns + state], 
                         expl[1*ncatg*ns + catg*ns + state], 
                         expl[2*ncatg*ns + catg*ns + state],
                         dot_prod[site*ns*ncatg]);
                }
            }
          
          fflush(NULL);
          assert(FALSE);
        }

      dlk /= lk;
      d2lk /= lk;

      dlnlk  += tree->data->wght[site] * dlk;
      d2lnlk += tree->data->wght[site] * d2lk - dlk*dlk;     
    }

  tree->c_dlnL  = dlnlk;
  tree->c_d2lnL = d2lnlk;

  return tree->c_lnL;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* Core of the likelihood calculation. Assume that the partial likelihoods on both
   sides of t_edge *b are up-to-date. Calculate the log-likelihood at one site.
   Note: this function can be used to evaluate first or second derivative of the 
   likelihood function with respect to the length of b, at a given site. Hence, be
   careful with the meaning of 'site_lk'. If 'derivative=TRUE', then site_lk is 
   either the first or second derivative of the likelihood at that site, given the
   length of edge 'b'.
*/

phydbl Lk_Core(int state, int ambiguity_check, short int derivative,
               phydbl *p_lk_left, phydbl *p_lk_rght,
               phydbl *Pij_rr,
               t_edge *b,
               t_tree *tree)
{
  phydbl site_lk,res,*pi;
  unsigned int catg;
  int num_prec_issue;
  
  const unsigned int ns    = tree->mod->ns;
  const unsigned int ncatg = tree->mod->ras->n_catg;
  const unsigned int site  = tree->curr_site;  
  const unsigned nsns      = ns*ns;

  pi = tree->mod->e_frq->pi->v;
  
  if(tree->mod->s_opt->skip_tree_traversal == NO)
    {
      for(catg=0;catg<ncatg;++catg)
        {
#if (defined(__AVX__))
          tree->site_lk_cat[catg] = AVX_Lk_Core_One_Class_No_Eigen_Lr(p_lk_left + catg*ns,
                                                                      (b->rght->tax) ? (p_lk_rght) : (p_lk_rght + catg*ns),
                                                                      Pij_rr + catg*nsns,
                                                                      pi,ns, ambiguity_check, state);
          
#elif (defined(__SSE3__))
          tree->site_lk_cat[catg] = SSE_Lk_Core_One_Class_No_Eigen_Lr(p_lk_left + catg*ns,
                                                                      (b->rght->tax) ? (p_lk_rght) : (p_lk_rght + catg*ns),
                                                                      Pij_rr + catg*nsns,
                                                                      pi,ns, ambiguity_check, state);
#else
          tree->site_lk_cat[catg] = Lk_Core_One_Class_No_Eigen_Lr(p_lk_left + catg*ns,
                                                                  (b->rght->tax) ? (p_lk_rght) : (p_lk_rght + catg*ns),
                                                                  Pij_rr + catg*nsns,
                                                                  pi,ns, ambiguity_check, state);
#endif


          /* printf("\n. site: %4d lk: %15G class: %3d", */
          /*        site, */
          /*        tree->site_lk_cat[catg], */
          /*        catg); */

        }

      Pull_Scaling_Factors(site,b,tree);

    }

  site_lk = .0;
  for(catg=0;catg<ncatg;++catg) site_lk += tree->unscaled_site_lk_cat[site*ncatg + catg] * tree->mod->ras->gamma_r_proba->v[catg];

  if(tree->mod->ras->invar == YES && derivative == NO)
    {
      num_prec_issue = NO;
      phydbl inv_site_lk = Invariant_Lk(tree->fact_sum_scale[site],site,&num_prec_issue,tree);

      switch(num_prec_issue)
        {
        case YES :         
          {
            site_lk = inv_site_lk * tree->mod->ras->pinvar->v;
            break;
          }
        case NO : 
          {
            site_lk = site_lk * (1. - tree->mod->ras->pinvar->v) + inv_site_lk * tree->mod->ras->pinvar->v;
            break;
          }
        }
    }
  
  // likelihood (or 1st, 2nd derivative) not rescaled here. Valid only if all partial likelihoods
  // were scaled using the same factor, i.e., when scaling_method == SCALE_FAST. In this case, the
  // scaling factors will cancel out in dlk/lk and d2lk/lk
  res = site_lk;
  
  if(derivative == NO)
    {
      phydbl log_site_lk;
      log_site_lk = log(site_lk) - (phydbl)LOG2 * tree->fact_sum_scale[site]; // log_site_lk =  log(site_lk_scaled / 2^(left_subtree+right_subtree))
      tree->c_lnL_sorted[site] = log_site_lk;
      tree->c_lnL += tree->data->wght[site] * log_site_lk;
      tree->cur_site_lk[site] = exp(log_site_lk); // note to self : add opt out option to avoid calculating this if not necessary
    }

  /* printf("\n. clnL: %f",tree->c_lnL); */
  return res;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl Lk_Core_Eigen_Lr(phydbl *expl, phydbl *dot_prod, short int derivative, t_edge *b, t_tree *tree)
{
  phydbl site_lk,res;
  unsigned int catg;
  int num_prec_issue;
  
  const unsigned int ns    = tree->mod->ns;
  const unsigned int ncatg = tree->mod->ras->n_catg;
  const unsigned int site  = tree->curr_site;
  
  if(tree->mod->s_opt->skip_tree_traversal == NO)
    {
      for(catg=0;catg<ncatg;++catg)
        {
#if (defined(__AVX__))
          tree->site_lk_cat[catg] = AVX_Lk_Core_One_Class_Eigen_Lr(dot_prod + catg*ns,
                                                                   expl ? (expl + catg*ns) : NULL,
                                                                   ns);
#elif (defined(__SSE3__))
          tree->site_lk_cat[catg] = SSE_Lk_Core_One_Class_Eigen_Lr(dot_prod + catg*ns,
                                                       expl ? (expl + catg*ns) : NULL,
                                                       ns);
#else
          tree->site_lk_cat[catg] = Lk_Core_One_Class_Eigen_Lr(dot_prod + catg*ns,
                                                               expl ? (expl + catg*ns) : NULL,
                                                               ns);
#endif
        }

      Pull_Scaling_Factors(site,b,tree);

    }
  
  
  site_lk = .0;
  for(catg=0;catg<ncatg;++catg) site_lk += tree->unscaled_site_lk_cat[site*ncatg+catg] * tree->mod->ras->gamma_r_proba->v[catg];
  
  if(tree->mod->ras->invar == YES && derivative == NO)
    {
      num_prec_issue = NO;
      phydbl inv_site_lk = Invariant_Lk(tree->fact_sum_scale[site],site,&num_prec_issue,tree);

      switch(num_prec_issue)
        {
        case YES :
          {
            site_lk = inv_site_lk * tree->mod->ras->pinvar->v;
            break;
          }
        case NO :
          {
            site_lk = site_lk * (1. - tree->mod->ras->pinvar->v) + inv_site_lk * tree->mod->ras->pinvar->v;
            break;
          }
        }
    }

  // likelihood (or 1st, 2nd derivative) not rescaled here. Valid only if all partial likelihoods
  // were scaled using the same factor, i.e., when scaling_method == SCALE_FAST. In this case, the
  // scaling factors will cancel out in dlk/lk and d2lk/lk
  res = site_lk;
  
  if(derivative == NO)
    {
      phydbl log_site_lk;
      log_site_lk = log(site_lk) - (phydbl)LOG2 * tree->fact_sum_scale[site]; // log_site_lk =  log(site_lk_scaled / 2^(left_subtree+right_subtree))
      tree->c_lnL_sorted[site] = log_site_lk;
      tree->c_lnL += tree->data->wght[site] * log_site_lk;
      tree->cur_site_lk[site] = exp(log_site_lk); // note to self : add opt out option to avoid calculating this if not necessary
    }
  
  return res;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

// Multiply matrices of eigenvectors with partial likelihoods
// to the left and right of edge b and take the dot-product of
// the two resulting vectors

void Update_Eigen_Lr(t_edge *b, t_tree *tree)
{
  unsigned int site,catg,state,i;
  int is_ambigu,observed_state;
  phydbl *dum,*dumdum;

  const unsigned int npattern = tree->n_pattern;
  const unsigned int ns = tree->mod->ns;
  const unsigned int ncatg = tree->mod->ras->n_catg;

  const unsigned int dim1 = ncatg * ns;
  const unsigned int dim2 = ns;

  if(tree->is_mixt_tree == YES)
    {      
      MIXT_Update_Eigen_Lr(b,tree);
      return;
    }

#if (defined(__AVX__))
  AVX_Update_Eigen_Lr(b,tree);
  return;
#elif (defined(__SSE3__))
  SSE_Update_Eigen_Lr(b,tree);
  return;
#endif

  dum    = (phydbl *)mCalloc(tree->mod->ns,sizeof(phydbl));
  dumdum = (phydbl *)mCalloc(tree->mod->ns,sizeof(phydbl));
  
  assert(tree->update_eigen_lr == YES);
  
  observed_state = -1;

  for(site=0;site<npattern;++site)
    {
      is_ambigu = YES;
      if(b->rght->tax == YES) is_ambigu = b->rght->c_seq->is_ambigu[site];
      if(is_ambigu == NO) observed_state = b->rght->c_seq->d_state[site];

      for(catg=0;catg<ncatg;++catg)
        {
          // Dot product left partial likelihood with equilibrium freqs
          for(state=0;state<ns;++state) dum[state] = b->p_lk_left[site*dim1 + catg*dim2 + state] * tree->mod->e_frq->pi->v[state];

          // Multiply by matrix of right eigen vectors
          for(state=0;state<ns;++state)
            {
              dumdum[state] = 0.0;
              for(i=0;i<ns;++i)
                {
                  dumdum[state] += 
                    dum[i] * 
                    tree->mod->eigen->r_e_vect[i*dim2 + state];
                }
            }

          for(state=0;state<ns;++state) tree->eigen_lr_left[catg*dim2 + state] = dumdum[state];

          if(b->rght->tax == YES)
            for(state=0;state<ns;++state) dum[state] = (phydbl)b->p_lk_tip_r[site*dim2 + state]; 
          else
            for(state=0;state<ns;++state) dum[state] = b->p_lk_rght[site*dim1 + catg*dim2 + state]; 

          
          if(is_ambigu == YES)
            {
              /* Multiply  matrix of left eigen vectors by vector of partial likelihoods on  */
              /* the righthand side of b */
              for(state=0;state<ns;++state)
                {
                  dumdum[state] = 0.0;
                  for(i=0;i<ns;i++)
                    {
                      dumdum[state] += 
                        tree->mod->eigen->l_e_vect[state*dim2 + i] *
                        dum[i];
                    }
                }
            }
          else
            {
              for(state=0;state<ns;++state)
                {
                  dumdum[state] = tree->mod->eigen->l_e_vect[state*dim2 + observed_state] * dum[observed_state];
                }
            }

          for(state=0;state<ns;++state) tree->eigen_lr_rght[catg*dim2 + state] = dumdum[state];
          
          for(state=0;state<ns;++state)
            {
              tree->dot_prod[site*dim1 + catg*dim2 + state] =
                tree->eigen_lr_rght[catg*dim2 + state] *
                tree->eigen_lr_left[catg*dim2 + state];
            }
        }
    }
  Free(dum);
  Free(dumdum);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Rate_Correction(int exponent, phydbl *site_lk_cat)
{
  int piecewise_exponent;
  phydbl multiplier,dum;
  unsigned long long int one = 1;
  
  dum = *site_lk_cat;
  if(exponent >= 0)
    {
      /* Multiply by 2^exponent */
      do
        {
          piecewise_exponent = MIN(exponent,63);
          multiplier = (phydbl)(one << piecewise_exponent);
          dum = dum * multiplier;
          exponent = exponent - piecewise_exponent;
        }
      while(exponent != 0);
    }
  else
    {
      /* Divide by 2^exponent */
      do
        {
          piecewise_exponent = MAX(exponent,-63);
          multiplier = 1. / (phydbl)(one << -piecewise_exponent);
          dum = dum * multiplier;
          exponent = exponent - piecewise_exponent;
        }
      while(exponent != 0);
    }

  *site_lk_cat = dum;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl Lk_Core_One_Class_Eigen_Lr(phydbl *dot_prod, phydbl *expl, int ns)
{
  unsigned int l;
  phydbl lk = 0.0;
  if(expl != NULL) for(l=0;l<ns;++l) lk += dot_prod[l] * expl[l]; 
  else for(l=0;l<ns;++l) lk += dot_prod[l];
  return lk;
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl Lk_Core_One_Class_No_Eigen_Lr(phydbl *p_lk_left, phydbl *p_lk_rght, phydbl *Pij,phydbl *pi, int ns, int ambiguity_check, int state)
{
  unsigned int l,k;
  phydbl lk = 0.0;
  phydbl sum;

    
  if(ambiguity_check == NO)/* tip case */
    {      
      sum = .0;
      for(l=0;l<ns;++l) sum += Pij[state*ns+l] * p_lk_left[l];               
      lk += sum * pi[state];
    }
  else /* If the character observed at the tip is ambiguous: ns x ns terms to consider */
    {
      for(k=0;k<ns;++k)
        {
          if(p_lk_rght[k] > .0) /* Only bother ascending into the subtrees if the likelihood of state k, at site "site*dim2" is > 0 */
            {
              sum = .0;
              for(l=0;l<ns;l++) sum += Pij[k*ns+l] * p_lk_left[l];                       
              lk += sum * pi[k] * p_lk_rght[k];
            }
        }
    } 

  return lk;
}

/* #endif */

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

// Returns the scaled likelihood for invariable sites
phydbl Invariant_Lk(int fact_sum_scale, int site, int *num_prec_issue, t_tree *tree)
{
  int exponent,piecewise_exponent;
  phydbl multiplier;
  phydbl inv_site_lk = 0.;
  
  (*num_prec_issue) = NO;
  
  /* The substitution model does include invariable sites */
  if(tree->mod->ras->invar == YES)
    {
      /* The site is invariant */
      if(tree->data->invar[site] > -0.5)
        {
          inv_site_lk = tree->mod->e_frq->pi->v[tree->data->invar[site]];
          
          /* printf("\n. inv_site_lk = %f [%c] [%d]",inv_site_lk,tree->data->c_seq[0]->state[site],tree->data->invar[site]); */

          if(tree->apply_lk_scaling == YES)
            {
              exponent = fact_sum_scale;              
              do
                {
                  piecewise_exponent = MIN(exponent,63);
                  multiplier = (phydbl)((unsigned long long)(1) << piecewise_exponent);
                  inv_site_lk *= multiplier;
                  exponent -= piecewise_exponent;
                }
              while(exponent != 0);
            }
          
          /* Update the value of site_lk */
          if(isinf(inv_site_lk)) // P(D|r=0) >> P(D|r>0) => assume P(D) = P(D|r=0)P(r=0)
            {
              int i;
              PhyML_Printf("\n. fact_sum_scale: %d",fact_sum_scale);              
              PhyML_Printf("\n. pi: %f",tree->mod->e_frq->pi->v[tree->data->invar[site]]);              
              for(i=0;i<tree->mod->ns;i++) PhyML_Printf("\n. pi %d: %f",i,tree->mod->e_frq->pi->v[i]);
              PhyML_Printf("\n. Numerical precision issue alert.");
              PhyML_Printf("\n. File %s at line %d (function '%s')\n",__FILE__,__LINE__,__FUNCTION__);
              (*num_prec_issue) = YES;
            }
        }
    }

  return inv_site_lk;

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* Update partial likelihood on edge b on the side of b where
   node d lies.
*/

void Update_Partial_Lk(t_tree *tree, t_edge *b, t_node *d)
{
  if(tree->eval_alnL == NO) return;
  
  if(tree->is_mixt_tree)
    {
      MIXT_Update_Partial_Lk(tree,b,d);
      return;
    }
  
  if((tree->io->do_alias_subpatt == YES) &&
     (tree->update_alias_subpatt == YES))
    Alias_One_Subpatt((d==b->left)?(b->rght):(b->left),d,tree);

  if(d->tax) return;

#ifdef BEAGLE
  update_beagle_partials(tree, b, d);
#else
  if(tree->mod->use_m4mod == NO)
    {
      if(tree->io->datatype == NT ||
         tree->io->datatype == AA)
        {
#if (defined(__AVX__))
          AVX_Update_Partial_Lk(tree,b,d);
#elif (defined(__SSE3__))
          SSE_Update_Partial_Lk(tree,b,d);
#else
          Default_Update_Partial_Lk(tree,b,d);
#endif
        }
      else
        {
          Update_Partial_Lk_Generic(tree,b,d);
        }
    }
  else
    {
      Update_Partial_Lk_Generic(tree,b,d);
    }
#endif


  //  Print_Edge_Likelihoods(tree, b, false);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

#ifndef BEAGLE

void Update_Partial_Lk_Generic(t_tree *tree, t_edge *b, t_node *d)
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
  phydbl p1_lk1,p2_lk2;
  phydbl *p_lk,*p_lk_v1,*p_lk_v2;
  phydbl *Pij1,*Pij2;
  phydbl *tPij1,*tPij2;
  int *sum_scale, *sum_scale_v1, *sum_scale_v2;
  int sum_scale_v1_val, sum_scale_v2_val;
  int i,j;
  int catg,site;
  int n_patterns;
  short int ambiguity_check_v1,ambiguity_check_v2;
  int state_v1,state_v2;
  phydbl smallest_p_lk,largest_p_lk;
  int *p_lk_loc;

  unsigned const int ncatg = tree->mod->ras->n_catg;
  unsigned const int ns = tree->mod->ns;
  unsigned const int ncatgns = ncatg * ns;
  unsigned const int nsns = ns * ns;
  

  if(tree->n_root && tree->ignore_root == YES &&
     (d == tree->n_root->v[1] || d == tree->n_root->v[2]) &&
     (b == tree->n_root->b[1] || b == tree->n_root->b[2]))
    {
      assert(FALSE);
    }


  state_v1 = state_v2 = -1;
  ambiguity_check_v1 = ambiguity_check_v2 = NO;
  sum_scale_v1_val = sum_scale_v2_val = 0;
  p1_lk1 = p2_lk2 = .0;

  if(d->tax)
    {
      PhyML_Fprintf(stderr,"\n. t_node %d is a leaf...",d->num);
      PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d (function '%s').\n",__FILE__,__LINE__,__FUNCTION__);
      Exit("\n");
    }

  n_patterns = tree->n_pattern;

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

  /* For every site in the alignment */
  for(site=0;site<n_patterns;site++)
    {
      state_v1 = state_v2 = -1;
      ambiguity_check_v1 = ambiguity_check_v2 = NO;

      if(tree->mod->s_opt->greedy == NO)
        {
          /* n_v1 and n_v2 are tip nodes */
          if(n_v1 && n_v1->tax)
            {
              /* Is the state at this tip ambiguous? */
              ambiguity_check_v1 = n_v1->c_seq->is_ambigu[site];
              /* if(ambiguity_check_v1 == NO) state_v1 = Get_State_From_Partial_Pars(n_v1->b[0]->p_lk_tip_r,site*ns,tree); */
              if(ambiguity_check_v1 == NO) state_v1 = n_v1->c_seq->d_state[site];
            }

          if(n_v2 && n_v2->tax)
            {
              /* Is the state at this tip ambiguous? */
              ambiguity_check_v2 = n_v2->c_seq->is_ambigu[site];
              /* ambiguity_check_v2 = tree->data->c_seq[n_v2->num]->is_ambigu[site]; */
              /* if(ambiguity_check_v2 == NO) state_v2 = Get_State_From_Partial_Pars(n_v2->b[0]->p_lk_tip_r,site*ns,tree); */
              if(ambiguity_check_v2 == NO) state_v2 = n_v2->c_seq->d_state[site];
            }
        }

      if(tree->mod->use_m4mod)
        {
          ambiguity_check_v1 = YES;
          ambiguity_check_v2 = YES;
        }

      /* For all the rate classes */
      for(catg=0;catg<ncatg;catg++)
        {
          if(tree->mod->ras->skip_rate_cat[catg] == YES) continue;
          
          smallest_p_lk  =  BIG;
          
          /* For all the states at node d */
          for(i=0;i<tree->mod->ns;i++)
            {
              p1_lk1 = .0;
              
              if(n_v1)
                {
                  /* n_v1 is a tip */
                  if((n_v1->tax) && (!tree->mod->s_opt->greedy))
                    {
                      if(ambiguity_check_v1 == NO)
                        {
                          /* For the (non-ambiguous) state at node n_v1 */
                          p1_lk1 = Pij1[catg*nsns+i*ns+state_v1];
                        }
                      else
                        {
                          /* For all the states at node n_v1 */
                          for(j=0;j<tree->mod->ns;j++)
                            {
                              p1_lk1 += Pij1[catg*nsns+i*ns+j] * (phydbl)n_v1->b[0]->p_lk_tip_r[site*ns+j];
                            }
                        }
                    }
                  /* n_v1 is an internal node */
                  else
                    {
                      /* For the states at node n_v1 */
                      for(j=0;j<tree->mod->ns;j++)
                        {
                          p1_lk1 += Pij1[catg*nsns+i*ns+j] * p_lk_v1[site*ncatgns+catg*ns+j];
                        }
                    }
                }
              else
                {
                  p1_lk1 = 1.0;
                }
              
              p2_lk2 = .0;
              
              /* We do exactly the same as for node n_v1 but for node n_v2 this time.*/
              if(n_v2)
                {
                  /* n_v2 is a tip */
                  if((n_v2->tax) && (!tree->mod->s_opt->greedy))
                    {
                      if(ambiguity_check_v2 == NO)
                        {
                          /* For the (non-ambiguous) state at node n_v2 */
                          p2_lk2 = Pij2[catg*nsns+i*ns+state_v2];
                        }
                      else
                        {
                          /* For all the states at node n_v2 */
                          for(j=0;j<tree->mod->ns;j++)
                            {
                              p2_lk2 += Pij2[catg*nsns+i*ns+j] * (phydbl)n_v2->b[0]->p_lk_tip_r[site*ns+j];
                            }
                        }
                    }
                  /* n_v2 is an internal node */
                  else
                    {
                      /* For all the states at node n_v2 */
                      for(j=0;j<tree->mod->ns;j++)
                        {
                          p2_lk2 += Pij2[catg*nsns+i*ns+j] * p_lk_v2[site*ncatgns+catg*ns+j];
                        }
                    }
                }
              else
                {
                  p2_lk2 = 1.0;
                }
              
              p_lk[site*ncatgns+catg*ns+i] = p1_lk1 * p2_lk2;
              
              /* 	      PhyML_Printf("\n+ %G",p_lk[site*ncatgns+catg*ns+i]); */
              
            }
          
          if(tree->scaling_method == SCALE_RATE_SPECIFIC)
            {
              smallest_p_lk = BIG;
              for(i=0;i<ns;++i)
                if(p_lk[site*ncatgns+catg*ns+i] < smallest_p_lk)
                  smallest_p_lk = p_lk[site*ncatgns+catg*ns+i];
              
              /* Current scaling values at that site */
              sum_scale_v1_val = (sum_scale_v1)?(sum_scale_v1[site*ncatg+catg]):(0);
              sum_scale_v2_val = (sum_scale_v2)?(sum_scale_v2[site*ncatg+catg]):(0);
              
              sum_scale[site*ncatg+catg] = sum_scale_v1_val + sum_scale_v2_val;
              
              /* Scaling. We have p_lk_lim_inf = 2^-500. Consider for instance that 
                 smallest_p_lk = 2^-600, then curr_scaler_pow will be equal to 100, and
                 each element in the partial likelihood vector will be multiplied by
                 2^100. */
              if(smallest_p_lk < (phydbl)P_LK_LIM_INF &&
                 tree->mod->augmented == NO &&
                 tree->apply_lk_scaling == YES &&
                 (n_v1->tax == NO || n_v2->tax == NO))
                
                {
                  int curr_scaler_pow;
                  curr_scaler_pow = (int)(-500.*LOG2-log(smallest_p_lk))/LOG2;
                  sum_scale[site*ncatg+catg] += curr_scaler_pow;
                  for(i=0;i<ns;++i) Rate_Correction(curr_scaler_pow, p_lk + site*nsns + catg*ns + i);
                  
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

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Default_Update_Partial_Lk(t_tree *tree, t_edge *b, t_node *d)
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
  phydbl p1_lk1,p2_lk2;//Partial likelihood at d's "left" neighbor, d's "right" neighbor
  phydbl *p_lk;
  phydbl *p_lk_v1,*p_lk_v2;//Partial likelihood vector of node d, d's "left" neighbor, d's "right" neighbor. We fill *p_lk, and assume *p_lk_v1 and *p_lk_v2 are already filled.
  phydbl *Pij1,*Pij2;
  phydbl *tPij1,*tPij2;
  int *sum_scale, *sum_scale_v1, *sum_scale_v2;
  int sum_scale_v1_val, sum_scale_v2_val;
  unsigned int i,j;
  unsigned int catg;
  unsigned int site;
  int ambiguity_check_v1,ambiguity_check_v2;
  int state_v1,state_v2;
  phydbl smallest_p_lk,largest_p_lk;
  phydbl p_lk_lim_inf;
  int *p_lk_loc;//Suppose site j, of a certain subtree, has "A" on one tip, and "C" on the other. If you come across this pattern again at site i<j, then you can simply copy the partial likelihoods
  
  const unsigned int ncatg = tree->mod->ras->n_catg;
  const unsigned int ns = tree->mod->ns;
  const unsigned int ncatgns = ncatg * ns;//Dimension of a matrix L that holds rate-specific character likelihoods. IOW, L[ij] is the likelihood of character j in rate class i
  const unsigned int nsns = ns * ns;//Dimensions of the transition prob. matrix
  const unsigned int n_patterns = tree->n_pattern;


  if(tree->n_root && tree->ignore_root == YES &&
     (d == tree->n_root->v[1] || d == tree->n_root->v[2]) &&
     (b == tree->n_root->b[1] || b == tree->n_root->b[2]))
    {
      assert(FALSE);
    }

  state_v1 = state_v2 = -1;
  ambiguity_check_v1 = ambiguity_check_v2 = NO;
  sum_scale_v1_val = sum_scale_v2_val = 0;
  p1_lk1 = p2_lk2 = .0;
  
  if(d->tax)
    {
      PhyML_Fprintf(stderr,"\n. t_node %d is a leaf...",d->num);
      PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d (function '%s')\n",__FILE__,__LINE__,__FUNCTION__);
      Exit("\n");
    }


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


  /* For every site in the alignment */
  for(site=0;site<n_patterns;site++)
    {
      state_v1 = state_v2 = -1;
      ambiguity_check_v1 = ambiguity_check_v2 = NO;

      if(!tree->mod->s_opt->greedy)
        {
          /* n_v1 and n_v2 are tip nodes */
          if(n_v1 && n_v1->tax)
            {
              /* Is the state at this tip ambiguous? */
              ambiguity_check_v1 = n_v1->c_seq->is_ambigu[site];
              if(ambiguity_check_v1 == NO) state_v1 = n_v1->c_seq->d_state[site];
            }
          
          if(n_v2 && n_v2->tax)
            {
              /* Is the state at this tip ambiguous? */
              ambiguity_check_v2 = n_v2->c_seq->is_ambigu[site];
              if(ambiguity_check_v2 == NO) state_v2 = n_v2->c_seq->d_state[site];
            }
          
          if(tree->mod->augmented == YES && n_v1 && n_v1->tax == NO)
            {
              state_v1 = Assign_State(n_v1->c_seq_anc->state+site*tree->mod->io->state_len,
                                      tree->io->datatype,
                                      tree->mod->io->state_len);
            }
          
          if(tree->mod->augmented == YES && n_v2 && n_v2->tax == NO)
            {
              state_v2 = Assign_State(n_v2->c_seq_anc->state+site*tree->mod->io->state_len,
                                      tree->io->datatype,
                                      tree->mod->io->state_len);
            }
        }
      
      if(tree->mod->use_m4mod)
        {
          ambiguity_check_v1 = YES;
          ambiguity_check_v2 = YES;
        }


      /* For all the rate classes */
      for(catg=0;catg<ncatg;++catg)
        {
          smallest_p_lk  =  BIG;
          
          /* For all states at node d */
          for(i=0;i<ns;++i)
            {
              if(tree->mod->augmented == YES)
                {
                  for(i=0;i<tree->mod->ns;i++) p_lk[site*ncatgns+catg*ns+i] = 0.0;
                  i = Assign_State(d->c_seq_anc->state+site*tree->mod->io->state_len,
                                   tree->io->datatype,
                                   tree->mod->io->state_len);
                }
              
              p1_lk1 = .0;
              
              if(n_v1)
                {
                  /* n_v1 is a tip */
                  if((n_v1->tax) && (!tree->mod->s_opt->greedy))
                    {
                      if(ambiguity_check_v1 == NO)
                        {
                          /* For the (non-ambiguous) state at node n_v1 */
                          p1_lk1 = Pij1[catg*nsns+i*ns+state_v1];
                        }
                      else
                        {
                          /* For all the states at node n_v1 */
                          p1_lk1 = .0;
                          for(j=0;j<ns;++j) p1_lk1 += Pij1[catg*nsns+i*ns+j] * p_lk_v1[site*ns+j];
                        }
                    }
                  /* n_v1 is an internal node */
                  else
                    {
                      //"catg*nsns" offsets into a rate-specific flat 4*4 DNA matrix P. Then, P[i,j] is given by "i*4+l"
                      //"site*ncatgns" offsets into a rate-specific flat num_rates*4 matrix L. Then, L[i,j] is given by "catg*ns+l"
                      if(tree->mod->augmented == YES)
                        {
                          p1_lk1 = Pij1[catg*nsns+i*ns+state_v1] * p_lk_v1[site*ncatgns+catg*ns+state_v1];
                        }
                      else
                        {
                          p1_lk1 = 0.0;
                          for(j=0;j<ns;++j) p1_lk1 += Pij1[catg*nsns+i*ns+j] * p_lk_v1[site*ncatgns+catg*ns+j];
                        }
                    }
                }
              else
                {
                  p1_lk1 = 1.0;
                }
              
              p2_lk2 = .0;
              
              /* We do exactly the same as for node n_v1 but for node n_v2 this time.*/
              if(n_v2)
                {
                  /* n_v2 is a tip */
                  if((n_v2->tax) && (!tree->mod->s_opt->greedy))
                    {
                      if(ambiguity_check_v2 == NO)
                        {
                          /* For the (non-ambiguous) state at node n_v2 */
                          p2_lk2 = Pij2[catg*nsns+i*ns+state_v2];
                        }
                      else
                        {
                          /* For all the states at node n_v2 */
                          p2_lk2 = 0.0;
                          for(j=0;j<ns;++j) p2_lk2 += Pij2[catg*nsns+i*ns+j] * p_lk_v2[site*ns+j];
                        }
                    }
                  /* n_v2 is an internal node */
                  else
                    {
                      if(tree->mod->augmented == YES)
                        {
                          p2_lk2 = Pij2[catg*nsns+i*ns+state_v2] * p_lk_v2[site*ncatgns+catg*ns+state_v2];
                        }
                      else
                        {
                          /* For all the states at node n_v2 */
                          p2_lk2 = 0.0;
                          for(j=0;j<ns;j++) p2_lk2 += Pij2[catg*nsns+i*ns+j] * p_lk_v2[site*ncatgns+catg*ns+j];
                        }
                    }
                }
              else
                {
                  p2_lk2 = 1.0;
                }
              
              /* Partial likelihood of character "i" at site "site" under rate "catg" */
              p_lk[site*ncatgns+catg*ns+i] = p1_lk1 * p2_lk2;
              
              if(tree->mod->augmented == YES) break;
            }
          
          if(tree->scaling_method == SCALE_RATE_SPECIFIC)
            {
              smallest_p_lk = BIG;
              for(i=0;i<ns;++i)
                if(p_lk[site*ncatgns+catg*ns+i] < smallest_p_lk)
                  smallest_p_lk = p_lk[site*ncatgns+catg*ns+i];
              
              /* Current scaling values at that site */
              sum_scale_v1_val = (sum_scale_v1)?(sum_scale_v1[site*ncatg+catg]):(0);
              sum_scale_v2_val = (sum_scale_v2)?(sum_scale_v2[site*ncatg+catg]):(0);
              
              sum_scale[catg*n_patterns+site] = sum_scale_v1_val + sum_scale_v2_val;
              
              
              /* Scaling. We have p_lk_lim_inf = 2^-500. Consider for instance that
                 smallest_p_lk = 2^-600, then curr_scaler_pow will be equal to 100, and
                 each element in the partial likelihood vector will be multiplied by
                 2^100. */
              if(smallest_p_lk < p_lk_lim_inf &&
                 tree->mod->augmented == NO &&
                 tree->apply_lk_scaling == YES &&
                 (n_v1->tax == NO || n_v2->tax == NO))
                {
                  int curr_scaler_pow;
                  curr_scaler_pow = (int)(-500.*LOG2-log(smallest_p_lk))/LOG2;
                  sum_scale[site*ncatg+catg] += curr_scaler_pow;
                  for(i=0;i<ns;++i) Rate_Correction(curr_scaler_pow, p_lk + site*nsns + catg*ns + i);
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

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl Return_Abs_Lk(t_tree *tree)
{
  Lk(NULL,tree);
  return FABS(tree->c_lnL);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

matrix *ML_Dist(calign *data, t_mod *mod)
{
  int i,j,k,l;
  phydbl init;
  int n_catg;
  phydbl d_max,sum;
  matrix *mat;
  calign *twodata,*tmpdata;
  int state0, state1;
  phydbl *F,len;
  eigen *eigen_struct;

  tmpdata         = (calign *)mCalloc(1,sizeof(calign));
  tmpdata->c_seq  = (align **)mCalloc(2,sizeof(align *));
  tmpdata->b_frq  = (phydbl *)mCalloc(mod->ns,sizeof(phydbl));
  tmpdata->ambigu = (short int *)mCalloc(data->crunch_len,sizeof(short int));
  F               = (phydbl *)mCalloc(mod->ns*mod->ns,sizeof(phydbl ));
  eigen_struct    = (eigen *)Make_Eigen_Struct(mod->ns);

  tmpdata->n_otu  = 2;

  tmpdata->crunch_len = data->crunch_len;
  tmpdata->init_len   = data->init_len;

  mat = NULL;
  if(mod->io->datatype == NT)           mat = (mod->whichmodel < 10)?(K80_dist(data,1E+6)):(JC69_Dist(data,mod));
  else if(mod->io->datatype == AA)      mat = JC69_Dist(data,mod);
  else if(mod->io->datatype == GENERIC) mat = JC69_Dist(data,mod);

  for(i=0;i<mod->ras->n_catg;i++) /* Don't use the discrete gamma distribution */
    {
      mod->ras->gamma_rr->v[i]      = 1.0;
      mod->ras->gamma_r_proba->v[i] = 1.0;
    }

  n_catg = mod->ras->n_catg;
  mod->ras->n_catg = 1;

  for(j=0;j<data->n_otu-1;j++)
    {
      tmpdata->c_seq[0]       = data->c_seq[j];
      tmpdata->c_seq[0]->name = data->c_seq[j]->name;
      tmpdata->wght           = data->wght;

      for(k=j+1;k<data->n_otu;k++)
        {
          tmpdata->c_seq[1]       = data->c_seq[k];
          tmpdata->c_seq[1]->name = data->c_seq[k]->name;
          
          twodata = Compact_Cdata(tmpdata,mod->io);
          
          Check_Ambiguities(twodata,mod->io->datatype,mod->io->state_len);
          
          Hide_Ambiguities(twodata);
          
          init = mat->dist[j][k];
          
          if((init > DIST_MAX-SMALL) || (init < .0)) init = 0.1;
          
          d_max = init;
          
          For(i,mod->ns*mod->ns) F[i]=.0;
          len = 0.0;
          For(l,twodata->c_seq[0]->len)
            {
              state0 = Assign_State(twodata->c_seq[0]->state+l*mod->io->state_len,mod->io->datatype,mod->io->state_len);
              state1 = Assign_State(twodata->c_seq[1]->state+l*mod->io->state_len,mod->io->datatype,mod->io->state_len);
              
              if((state0 > -1) && (state1 > -1))
                {
                  F[mod->ns*state0+state1] += twodata->wght[l];
                  len += twodata->wght[l];
                }
            }
          
          if(len > .0) 
            {
              For(i,mod->ns*mod->ns) F[i] /= len;
            }
          
          sum = 0.;
          For(i,mod->ns*mod->ns) sum += F[i];
          
          /* if(sum < .001) d_max = -1.; */
          if(sum < .001) d_max = init;
          else if((sum > 1. - .001) && (sum < 1. + .001)) Opt_Dist_F(&(d_max),F,mod);
          else
            {
              PhyML_Fprintf(stderr,"\n\n. sum = %f",sum);
              PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
              Exit("");
            }
          
          if(d_max >= DIST_MAX) d_max = DIST_MAX;
          
          
          /* Do not correct for dist < BL_MIN, otherwise Fill_Missing_Dist
           *  will not be called
           */
          mat->dist[j][k] = d_max;
          mat->dist[k][j] = mat->dist[j][k];
          Free_Calign(twodata);
        }
    }
  
  mod->ras->n_catg = n_catg;
  
  
  Free(tmpdata->ambigu);
  Free(tmpdata->b_frq);
  Free(tmpdata->c_seq);
  free(tmpdata);
  Free_Eigen(eigen_struct);
  Free(F);

  return mat;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Lk_Given_Two_Seq(calign *data, int numseq1, int numseq2, phydbl dist, t_mod *mod, phydbl *loglk)
{
  align *seq1,*seq2;
  phydbl site_lk,log_site_lk;
  int i,j,k,l;
/*   phydbl **p_lk_l,**p_lk_r; */
  phydbl *p_lk_l,*p_lk_r;
  phydbl len;
  int dim1,dim2;

  dim1 = mod->ns;
  dim2 = mod->ns * mod->ns;

  DiscreteGamma(mod->ras->gamma_r_proba->v, mod->ras->gamma_rr->v, mod->ras->alpha->v,
        mod->ras->alpha->v,mod->ras->n_catg,mod->ras->gamma_median);

  seq1 = data->c_seq[numseq1];
  seq2 = data->c_seq[numseq2];


  p_lk_l = (phydbl *)mCalloc(data->c_seq[0]->len * mod->ns,sizeof(phydbl));
  p_lk_r = (phydbl *)mCalloc(data->c_seq[0]->len * mod->ns,sizeof(phydbl));


  for(i=0;i<mod->ras->n_catg;i++)
    {
      len = dist*mod->ras->gamma_rr->v[i];
      if(len < mod->l_min) len = mod->l_min;
      else if(len > mod->l_max) len = mod->l_max;
      PMat(len,mod,dim2*i,mod->Pij_rr->v,NULL);
    }
  
  if(mod->io->datatype == NT)
    {
      For(i,data->c_seq[0]->len)
        {
          Init_Tips_At_One_Site_Nucleotides_Float(seq1->state[i],i*mod->ns,p_lk_l);
          Init_Tips_At_One_Site_Nucleotides_Float(seq2->state[i],i*mod->ns,p_lk_r);
        }
    }
  else if(mod->io->datatype == AA)
    {
      For(i,data->c_seq[0]->len)
        {
          Init_Tips_At_One_Site_AA_Float(seq1->state[i],i*mod->ns,p_lk_l);
          Init_Tips_At_One_Site_AA_Float(seq2->state[i],i*mod->ns,p_lk_r);
        }
    }
  else
    {
      PhyML_Fprintf(stderr,"\n\n. Not implemented yet...");
      PhyML_Fprintf(stderr,"\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
      Warn_And_Exit("\n");
    }
  
  
  site_lk = .0;
  *loglk = 0;
  
  For(i,data->c_seq[0]->len)
    {
      if(data->wght[i] > 0.0)
        {
          site_lk = log_site_lk = .0;
          if(!data->ambigu[i])
            {
              for(k=0;k<mod->ns;k++) {if(p_lk_l[i*mod->ns+k] > .0001) break;}
              for(l=0;l<mod->ns;l++) {if(p_lk_r[i*mod->ns+l] > .0001) break;}
              for(j=0;j<mod->ras->n_catg;j++)
                {
                  site_lk +=
                    mod->ras->gamma_r_proba->v[j] *
                    mod->e_frq->pi->v[k] *
                    p_lk_l[i*dim1+k] *
                    mod->Pij_rr->v[j*dim2+k*dim1+l] *
                    p_lk_r[i*dim1+l];
                }
            }
          else
            {
              for(j=0;j<mod->ras->n_catg;j++)
                {
                  for(k=0;k<mod->ns;k++) /*sort sum terms ? No global effect*/
                    {
                      for(l=0;l<mod->ns;l++)
                        {
                          site_lk +=
                            mod->ras->gamma_r_proba->v[j] *
                            mod->e_frq->pi->v[k] *
                            p_lk_l[i*dim1+k] *
                            mod->Pij_rr->v[j*dim2+k*dim1+l] *
                            p_lk_r[i*dim1+l];
                        }
                    }
                }
            }
          
          if(site_lk <= .0)
            {
              PhyML_Fprintf(stderr,"\n\n. '%c' '%c'\n",seq1->state[i],seq2->state[i]);
              Exit("\n. Err: site lk <= 0\n");
            }
          
          log_site_lk += (phydbl)log(site_lk);
          
          *loglk += data->wght[i] * log_site_lk;/* sort sum terms ? No global effect*/
    }
    }

/*   For(i,data->c_seq[0]->len) */
/*     { */
/*       Free(p_lk_l[i]); */
/*       Free(p_lk_r[i]); */
/*     } */

  Free(p_lk_l); Free(p_lk_r);
  return *loglk;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Unconstraint_Lk(t_tree *tree)
{
  int i;

  tree->unconstraint_lk = .0;
  for(i=0;i<tree->data->crunch_len;i++) tree->unconstraint_lk += tree->data->wght[i]*(phydbl)log(tree->data->wght[i]);
  tree->unconstraint_lk -= tree->data->init_len*(phydbl)log(tree->data->init_len);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Init_Partial_Lk_Tips_Double(t_tree *tree)
{
  unsigned int curr_site,i,dim1;

  dim1 = tree->mod->ns;


  for(i=0;i<tree->n_otu;i++)
    {
      if(!tree->a_nodes[i]->c_seq || 
	 strcmp(tree->a_nodes[i]->c_seq->name,tree->a_nodes[i]->name))
        {
          PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
          Exit("");
        }
    }

  for(curr_site=0;curr_site<tree->data->crunch_len;curr_site++)
    {
      for(i=0;i<tree->n_otu;i++)
        {
          if (tree->io->datatype == NT)
            Init_Tips_At_One_Site_Nucleotides_Float(tree->a_nodes[i]->c_seq->state[curr_site],
                                                    curr_site*dim1,
                                                    tree->a_nodes[i]->b[0]->p_lk_tip_r);

          else if(tree->io->datatype == AA)
            Init_Tips_At_One_Site_AA_Float(tree->a_nodes[i]->c_seq->state[curr_site],
                                           curr_site*dim1,
                                           tree->a_nodes[i]->b[0]->p_lk_tip_r);

          else if(tree->io->datatype == GENERIC)
            Init_Tips_At_One_Site_Generic_Float(tree->a_nodes[i]->c_seq->state+curr_site*tree->mod->io->state_len,
                                                tree->mod->ns,
                                                tree->mod->io->state_len,
                                                curr_site*dim1,
                                                tree->a_nodes[i]->b[0]->p_lk_tip_r);
        }
    }

  if(tree->mod->use_m4mod)
    {
      M4_Init_Partial_Lk_Tips_Double(tree);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Init_Partial_Lk_Tips_Int(t_tree *tree)
{
  int curr_site,i,dim1;

  dim1 = tree->mod->ns;

  for(i=0;i<tree->n_otu;i++)
    {
      if(!tree->a_nodes[i]->c_seq || 
	 strcmp(tree->a_nodes[i]->c_seq->name,tree->a_nodes[i]->name))
        {
          PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
          Exit("");
        }
    }

  for(curr_site=0;curr_site<tree->data->crunch_len;curr_site++)
    {
      for(i=0;i<tree->n_otu;i++)
        {
          /* printf("\n. site: %3d %c",curr_site,tree->a_nodes[i]->c_seq->state[curr_site]); */
          /* printf("\n. init at %s %p",tree->a_nodes[i]->name,tree->a_nodes[i]->b[0]->p_lk_tip_r); fflush(NULL); */
          if(tree->io->datatype == NT)
            {
              Init_Tips_At_One_Site_Nucleotides_Float(tree->a_nodes[i]->c_seq->state[curr_site],
                                                    curr_site*dim1,
                                                    tree->a_nodes[i]->b[0]->p_lk_tip_r);
              /* Init_Tips_At_One_Site_Nucleotides_Int(tree->data->c_seq[i]->state[curr_site], */
              /* 					    curr_site*dim1, */
              /* 					    tree->a_nodes[i]->b[0]->p_lk_tip_r); */
            }
          else if(tree->io->datatype == AA)
            Init_Tips_At_One_Site_AA_Float(tree->a_nodes[i]->c_seq->state[curr_site],
                                           curr_site*dim1,
                                           tree->a_nodes[i]->b[0]->p_lk_tip_r);
          /* Init_Tips_At_One_Site_AA_Int(tree->data->c_seq[i]->state[curr_site], */
          /* 				 curr_site*dim1,					    */
          /* 				 tree->a_nodes[i]->b[0]->p_lk_tip_r); */
          
          else if(tree->io->datatype == GENERIC)
            {
              Init_Tips_At_One_Site_Generic_Float(tree->a_nodes[i]->c_seq->state+curr_site*tree->mod->io->state_len,
                                                  tree->mod->ns,
                                                  tree->mod->io->state_len,
                                                  curr_site*dim1,
                                                  tree->a_nodes[i]->b[0]->p_lk_tip_r);
              
              /* Init_Tips_At_One_Site_Generic_Int(tree->data->c_seq[i]->state+curr_site*tree->mod->io->state_len, */
              /* 					tree->mod->ns, */
              /* 					tree->mod->io->state_len, */
              /* 					curr_site*dim1, */
              /* 					tree->a_nodes[i]->b[0]->p_lk_tip_r); */
            }
#ifdef BEAGLE
          //Recall that tip partials are stored on the branch leading
          //to the tip, rather than on the tip itself (hence `p_lk_tip_idx`
          //is a field of the branch (i.e. b[0]) rather than the node.
          //Secondly, the BEAGLE's partial buffers are laid out as
          //BEAGLE's partials buffer = [ tax1, tax2, ..., taxN, b1Left, b2Left, b3Left,...,bMLeft, b1Rght, b2Rght, b3Rght,...,bMRght] (N taxa, M branches)
        tree->a_nodes[i]->b[0]->p_lk_tip_idx = i;
#endif
      }
    }
   if(tree->mod->use_m4mod)
     {
       M4_Init_Partial_Lk_Tips_Int(tree);
     }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Update_PMat_At_Given_Edge(t_edge *b_fcus, t_tree *tree)
{
  int i;
  phydbl len;
  phydbl l_min, l_max;
  phydbl shape, scale, mean, var;

  assert(b_fcus);
  assert(tree);

  
  if(tree->is_mixt_tree == YES)
    {
      MIXT_Update_PMat_At_Given_Edge(b_fcus,tree);
      return;
    }

  if(tree->io->mod->gamma_mgf_bl == YES) Set_Br_Len_Var(b_fcus,tree);

  l_min = tree->mod->l_min;
  l_max = tree->mod->l_max;

  len = -1.0;

  if(tree->mod->log_l == YES) b_fcus->l->v = exp(b_fcus->l->v);

  if(b_fcus->l->v < l_min) b_fcus->l->v = l_min;
  if(b_fcus->l->v > l_max) b_fcus->l->v = l_max;

  for(i=0;i<tree->mod->ras->n_catg;i++)
    {
      if(tree->mod->ras->skip_rate_cat[i] == YES) continue;
      
      //Update the branch length
      if(b_fcus->has_zero_br_len == YES)
        {
#ifdef BEAGLE
          Warn_And_Exit(TODO_BEAGLE);
#endif
          len = -1.0;
          mean = -1.0;
          var  = -1.0;
        }
      else
        {
          len = MAX(0.0,b_fcus->l->v)*tree->mod->ras->gamma_rr->v[i];//branch_len * rate
          len *= tree->mod->br_len_mult->v;
          if(tree->mixt_tree)  len *= tree->mixt_tree->mod->ras->gamma_rr->v[tree->mod->ras->parent_class_number];
          if(len < l_min)      len = l_min;
          else if(len > l_max) len = l_max;
          
          mean = len;
          var  = MAX(0.0,b_fcus->l_var->v) * POW(tree->mod->ras->gamma_rr->v[i]*tree->mod->br_len_mult->v,2);
          if(tree->mixt_tree) var *= POW(tree->mixt_tree->mod->ras->gamma_rr->v[tree->mod->ras->parent_class_number],2);

          /* var = 1.E-10; */

          if(var > tree->mod->l_var_max) var = tree->mod->l_var_max;
          if(var < tree->mod->l_var_min) var = tree->mod->l_var_min;
        }

      //Update the transition prob. matrix
      if(tree->mod->gamma_mgf_bl == NO)
          {
#ifdef BEAGLE
            assert(UNINITIALIZED != tree->mod->b_inst);
#endif

            PMat(len,tree->mod,i*tree->mod->ns*tree->mod->ns,b_fcus->Pij_rr,b_fcus->tPij_rr);
          }
      else
          {
#ifdef BEAGLE
            Warn_And_Exit(TODO_BEAGLE);
#endif
            
            shape = mean*mean/var;
            scale = var/mean;
            PMat_MGF_Gamma(b_fcus->Pij_rr+tree->mod->ns*tree->mod->ns*i,shape,scale,1.0,tree->mod);
          }
    }

#ifdef BEAGLE
  int whichmodel = tree->mod->whichmodel;
  //Only for some models we use Beagle to compute/update the P-matrices, for other models
  //we compute them in PhyML and explicitly set the P-matrices in BEAGLE
  if((tree->mod->io->datatype == AA || whichmodel==GTR || whichmodel==CUSTOM) && tree->mod->use_m4mod == NO)
    {
      if(b_fcus->has_zero_br_len == YES)
        Warn_And_Exit(TODO_BEAGLE);
      
      //
      update_beagle_eigen(tree->mod);
      update_beagle_ras(tree->mod);
      
      //
      len = MAX(0.0, b_fcus->l->v) * tree->mod->br_len_mult->v;
      int p_matrices[1]     = b_fcus->Pij_rr_idx;
      double branch_lens[1] = len;
      int ret = beagleUpdateTransitionMatrices(tree->b_inst,0,p_matrices,NULL,NULL,branch_lens,1);
      if(ret<0)
        {
          PhyML_Fprintf(stderr, "beagleUpdateTransitionMatrices() on instance %i failed:%i\n\n",tree->b_inst,ret);
          Exit("");
        }
      //Retrieve a "local" copy of the P-matrix
      ret = beagleGetTransitionMatrix(tree->b_inst, b_fcus->Pij_rr_idx, b_fcus->Pij_rr);
      if(ret<0)
        {
          PhyML_Fprintf(stderr, "beagleGetTransitionMatrix() on instance %i failed:%i\n\n",tree->b_inst,ret);
          Exit("");
        }
    }
  else
    {
      int ret = beagleSetTransitionMatrix(tree->b_inst, b_fcus->Pij_rr_idx, b_fcus->Pij_rr, -1);
      if(ret<0)
        {
          PhyML_Fprintf(stderr, "beagleSetTransitionMatrix() on instance %i failed:%i\n\n",tree->b_inst,ret);
          Exit("");
        }
  }
#endif
    if(tree->mod->log_l == YES) b_fcus->l->v = log(b_fcus->l->v);

//      Print_Model(tree->mod);
//      Dump_Arr_D(tree->cur_site_lk, tree->n_pattern);
//      Print_Edge_PMats(tree, b_fcus);
//      Print_Edge_Likelihoods(tree,b_fcus,true);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


/* void Update_Partial_Lk_On_A_Path(t_node *a, t_node *d, t_edge *b, t_node *target_one_side, t_node *target_other_side, t_tree *tree) */
/* { */


/*   /\* */
/*                 \               / */
/* 	   target\___________b_/ */
/* 		 /  \	\  \   \ */
/* 		/    \	 \  \	\ */

/*     target is the root of the subtree at which we want */
/*     the likelihood to be updated */
/*   *\/ */



/* /\*   PhyML_Printf("Update_p_lk on (%d %d) at %d (target=%d %d)\n", *\/ */
/* /\* 	 b->left->num, *\/ */
/* /\* 	 b->rght->num, *\/ */
/* /\* 	 a->num, *\/ */
/* /\* 	 target_one_side->num, *\/ */
/* /\* 	 target_other_side->num); *\/ */

/*   Update_Partial_Lk(tree,b,a); */
/*   if((a == target_one_side) && (d == target_other_side))  */
/*     return; */
/*   else */
/*     { */
/*       Update_Partial_Lk_On_A_Path(d, */
/* 			    d->v[tree->t_dir[d->num][target_other_side->num]], */
/* 			    d->b[tree->t_dir[d->num][target_other_side->num]], */
/* 			    target_one_side, */
/* 			    target_other_side, */
/* 			    tree);  */
/*     } */
/* } */

void Update_Partial_Lk_Along_A_Path(t_node **path, int path_length, t_tree *tree)
{
  int i,j;
  
  for(i=0;i<path_length-1;++i)
    {
      for(j=0;j<3;++j)
        if(path[i]->v[j] == path[i+1])
          {
            if(path[i] == path[i]->b[j]->left)
              {
                Update_Partial_Lk(tree,path[i]->b[j],path[i]->b[j]->left);
              }
            else if(path[i] == path[i]->b[j]->rght)
              {
                Update_Partial_Lk(tree,path[i]->b[j],path[i]->b[j]->rght);
              }
            else
              {
                PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d. \n",__FILE__,__LINE__);
                Exit("");
              }
            break;
          }
#ifdef DEBUG
      if(j == 3)
        {
          PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d.\n",__FILE__,__LINE__);
          Exit("");
        }
#endif
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl Lk_Dist(phydbl *F, phydbl dist, t_mod *mod)
{
  int i,j,k;
  phydbl lnL,len;
  int dim1,dim2;
  phydbl pi, pijk;

  // Compute likelihood of the model made of the
  // first class of the mixture.
  if(mod->is_mixt_mod == YES) mod = mod->next;
  assert(mod);

  if(mod->log_l == YES) dist = exp(dist);

  for(k=0;k<mod->ras->n_catg;k++)
    {
      len = dist*mod->ras->gamma_rr->v[k];
      if(len < mod->l_min)      len = mod->l_min;
      else if(len > mod->l_max) len = mod->l_max;
      PMat(len,mod,mod->ns*mod->ns*k,mod->Pij_rr->v,NULL);
    }

  dim1 = mod->ns*mod->ns;
  dim2 = mod->ns;
  lnL  = .0;
  pi   = -1.;
  pijk = -1.;

  for(i=0;i<mod->ns-1;i++)
    {
      pi = mod->e_frq->pi->v[i];

      for(j=i+1;j<mod->ns;j++)
        {
          for(k=0;k<mod->ras->n_catg;k++)
            {
              pijk = mod->Pij_rr->v[dim1*k+dim2*i+j];              
              lnL += (F[dim1*k+dim2*i+j] + F[dim1*k+dim2*j+i]) * log(pi * pijk);

              /* printf("\nXXXXXX i: %d j: %d f: %G Pij:%G F:%G %G %G", */
              /*        i,j,mod->e_frq->pi->v[i], */
              /*        mod->Pij_rr->v[dim1*k+dim2*i+j], */
              /*        F[dim1*k+dim2*j+i],lnL, */
              /*        (F[dim1*k+dim2*i+j] + F[dim1*k+dim2*j+i])*log(pi * pijk)); */
            }
        }
    }
  
  for(i=0;i<mod->ns;i++) 
    {
      pi = mod->e_frq->pi->v[i];
      
      for(k=0;k<mod->ras->n_catg;k++) 
        {
          pijk = mod->Pij_rr->v[dim1*k+dim2*i+i];
          lnL += F[dim1*k+dim2*i+i]* log(pi * pijk);

          /* printf("\nYYYYYY i: %d j: %d f: %G Pij:%G F:%G lnL:%G pi:%G pijk:%G", */
          /*        i,j,mod->e_frq->pi->v[i], */
          /*        mod->Pij_rr->v[dim1*k+dim2*i+j], */
          /*        F[dim1*k+dim2*j+i],lnL, */
          /*        pi,pijk); */
        }
    }

  return lnL;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl Update_Lk_At_Given_Edge(t_edge *b_fcus, t_tree *tree)
{
  Update_Partial_Lk(tree,b_fcus,b_fcus->left);
  Update_Partial_Lk(tree,b_fcus,b_fcus->rght);
  tree->c_lnL = Lk(b_fcus,tree);
  return tree->c_lnL;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Print_Lk_Given_Edge_Recurr(t_node *a, t_node *d, t_edge *b, t_tree *tree)
{
  PhyML_Printf("\n___ Edge %3d (left=%3d rght=%3d) lnL=%f",
     b->num,
     b->left->num,
     b->rght->num,
     Lk(b,tree));

  if(d->tax) return;
  else
    {
      int i;
      for(i=0;i<3;i++)
    if(d->v[i] != a)
      Print_Lk_Given_Edge_Recurr(d,d->v[i],d->b[i],tree);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Alias_Subpatt(t_tree *tree)
{

  if(tree->n_root && tree->ignore_root == NO)
    {
      Alias_Subpatt_Post(tree->n_root,tree->n_root->v[2],tree);
      Alias_Subpatt_Post(tree->n_root,tree->n_root->v[1],tree);
    }
  else
    {
      Alias_Subpatt_Post(tree->a_nodes[0],tree->a_nodes[0]->v[0],tree);
      /* if(tree->both_sides)  */
      Alias_Subpatt_Pre(tree->a_nodes[0],tree->a_nodes[0]->v[0],tree);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Alias_One_Subpatt(t_node *a, t_node *d, t_tree *tree)
{
  int i,j;
  int *patt_id_v1, *patt_id_v2, *patt_id_d;
  int *p_lk_loc_d, *p_lk_loc_v1, *p_lk_loc_v2;
  t_node *v1, *v2;
  t_edge *b0, *b1, *b2;
  int curr_patt_id_v1, curr_patt_id_v2;
  int curr_p_lk_loc_v1, curr_p_lk_loc_v2;
  int num_subpatt;

  b0 = b1 = b2 = NULL;

  if(d->tax)
    {
      patt_id_d  = (d == d->b[0]->left)?(d->b[0]->patt_id_left):(d->b[0]->patt_id_rght);
      p_lk_loc_d = (d == d->b[0]->left)?(d->b[0]->p_lk_loc_left):(d->b[0]->p_lk_loc_rght);

      for(i=0;i<tree->n_pattern;i++)
    {
      for(j=0;j<tree->n_pattern;j++)
        {
          if(patt_id_d[i] == patt_id_d[j])
        {
          p_lk_loc_d[i] = j;
          break;
        }
          if(j > i)
        {
          PhyML_Fprintf(stderr,"\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
          Warn_And_Exit("");
        }
        }
    }
      return;
    }
  else
    {
      v1 = v2 = NULL;
      for(i=0;i<3;i++)
    {
      if(d->v[i] != a && d->b[i] != tree->e_root)
        {
          if(!v1) { v1=d->v[i]; b1=d->b[i];}
          else    { v2=d->v[i]; b2=d->b[i];}
        }
      else
        {
          b0 = d->b[i];
        }
    }


      patt_id_v1  = (v1 == b1->left)?(b1->patt_id_left):(b1->patt_id_rght);
      patt_id_v2  = (v2 == b2->left)?(b2->patt_id_left):(b2->patt_id_rght);
      patt_id_d   = (d  == b0->left)?(b0->patt_id_left):(b0->patt_id_rght);
      p_lk_loc_d  = (d  == b0->left)?(b0->p_lk_loc_left):(b0->p_lk_loc_rght);
      p_lk_loc_v1 = (v1 == b1->left)?(b1->p_lk_loc_left):(b1->p_lk_loc_rght);
      p_lk_loc_v2 = (v2 == b2->left)?(b2->p_lk_loc_left):(b2->p_lk_loc_rght);

      num_subpatt = 0;
      for(i=0;i<tree->n_pattern;i++)
    {
      curr_patt_id_v1  = patt_id_v1[i];
      curr_patt_id_v2  = patt_id_v2[i];
      curr_p_lk_loc_v1 = p_lk_loc_v1[i];
      curr_p_lk_loc_v2 = p_lk_loc_v2[i];

      p_lk_loc_d[i] = i;

      if((curr_p_lk_loc_v1 == i) || (curr_p_lk_loc_v2 == i))
        {
          p_lk_loc_d[i] = i;
          patt_id_d[i] = num_subpatt;
          num_subpatt++;
        }
      else
        if(curr_p_lk_loc_v1 == curr_p_lk_loc_v2)
          {
        p_lk_loc_d[i] = curr_p_lk_loc_v1;
        patt_id_d[i] = patt_id_d[curr_p_lk_loc_v1];
          }
        else
          {
        for(j=MAX(curr_p_lk_loc_v1,curr_p_lk_loc_v2);j<tree->n_pattern;j++)
          {
            if((patt_id_v1[j] == curr_patt_id_v1) &&
               (patt_id_v2[j] == curr_patt_id_v2))
              {
            p_lk_loc_d[i] = j;

            if(j == i)
              {
                patt_id_d[i] = num_subpatt;
                num_subpatt++;
              }
            else patt_id_d[i] = patt_id_d[j];
            break;
              }
            if(j > i)
              {
                PhyML_Fprintf(stderr,"\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
                Warn_And_Exit("");
              }
          }
          }
    }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Alias_Subpatt_Post(t_node *a, t_node *d, t_tree *tree)
{

  if(d->tax) return;
  else
    {
      int i;

      for(i=0;i<3;i++)
    {
      if(d->v[i] != a && d->b[i] != tree->e_root)
        {
          Alias_Subpatt_Post(d,d->v[i],tree);
        }
    }
      Alias_One_Subpatt(a, d, tree);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Alias_Subpatt_Pre(t_node *a, t_node *d, t_tree *tree)
{
  if(d->tax) return;
  else
    {
      int i;
      
      for(i=0;i<3;++i)
        {
          if(d->v[i] != a && d->b[i] != tree->e_root)
            {
              Alias_One_Subpatt(d->v[i],d,tree);
              Alias_Subpatt_Pre(d,d->v[i],tree);
            }
        }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Init_Partial_Lk_Loc(t_tree *tree)
{
  int i,j;
  t_node *d;
  int *patt_id_d;

  For(i,2*tree->n_otu-1)
    {
      for(j=0;j<tree->n_pattern;j++)
        {
          tree->a_edges[i]->p_lk_loc_left[j] = j;
          tree->a_edges[i]->p_lk_loc_rght[j] = j;
        }
    }
  
  for(i=0;i<tree->n_otu;i++)
    {
      d = tree->a_nodes[i];
      patt_id_d = (d == d->b[0]->left)?(d->b[0]->patt_id_left):(d->b[0]->patt_id_rght);
      for(j=0;j<tree->n_pattern;j++)
        {
          patt_id_d[j] = (int)tree->a_nodes[d->num]->c_seq->state[j];
        }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Lk_Normal_Approx(t_tree *tree)
{
  phydbl lnL;
  int i;
  int dim;
  phydbl first_order;

  dim = 2*tree->n_otu-3;

  lnL = Dnorm_Multi_Given_InvCov_Det(tree->rates->u_cur_l,
                     tree->rates->mean_l,
                     tree->rates->invcov,
                     tree->rates->covdet,
                     2*tree->n_otu-3,YES);

  first_order = 0.0;
  for(i=0;i<dim;i++) first_order += (tree->rates->u_cur_l[i] - tree->rates->mean_l[i]) * tree->rates->grad_l[i];

  lnL += first_order;

  /* printf("\n"); */
  /* for(i=0;i<dim;i++) printf("%f\t",tree->rates->u_cur_l[i]); */
  /* printf("\n. Lk=%f %f",lnL,tree->mod->l_min); */


/*   err = NO; */
/*   dim = 2*tree->n_otu-3; */
/*   for(i=0;i<dim;i++) */
/*     { */
/*       if((tree->rates->mean_l[i] / tree->mod->l_min < 1.1) &&  */
/* 	 (tree->rates->mean_l[i] / tree->mod->l_min > 0.9)) */
/* 	{	   */
/* 	  lnL -= Log_Dnorm(tree->a_edges[i]->l->v,tree->rates->mean_l[i],SQRT(tree->rates->cov_l[i*dim+i]),&err); */
/* 	  if(err) */
/* 	    { */
/* 	      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__); */
/* 	      Warn_And_Exit(""); */
/* 	    } */
/* 	  lambda = 1./SQRT(tree->rates->cov_l[i*dim+i]); */
/* 	  l = (tree->mod->log_l == YES)?(exp(tree->a_edges[i]->l->v)):(tree->a_edges[i]->l->v); */
/* 	  lnL += log(Dexp(l,lambda)); */
/* /\* 	  printf("\n. lambda = %f",lambda); *\/ */
/* 	} */
/*     } */

  return(lnL);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Wrap_Part_Lk_At_Given_Edge(t_edge *b, t_tree *tree, supert_tree *stree)
{
  return -1.0;
  /* return PART_Lk_At_Given_Edge(b,stree);; */
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl Wrap_Part_Lk(t_edge *b, t_tree *tree, supert_tree *stree)
{
  return -1.0;
  /* return PART_Lk(stree); */
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Wrap_Lk(t_edge *b, t_tree *tree, supert_tree *stree)
{
  Lk(NULL,tree);
  return tree->c_lnL;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Wrap_Geo_Lk(t_edge *b, t_tree *tree, supert_tree *stree)
{
  TIPO_Lk(tree);
  return tree->geo_lnL;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Wrap_Lk_At_Given_Edge(t_edge *b, t_tree *tree, supert_tree *stree)
{
  Lk(b,tree);
  return tree->c_lnL;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Wrap_Lk_Rates(t_edge *b, t_tree *tree, supert_tree *stree)
{
  RATES_Lk_Rates(tree);
  return tree->rates->c_lnL_rates;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl Wrap_Lk_Times(t_edge *b, t_tree *tree, supert_tree *stree)
{
  TIMES_Lk_Times(NO,tree);
  return tree->rates->c_lnL_times;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl Wrap_Lk_Linreg(t_edge *b, t_tree *tree, supert_tree *stree)
{
  /* RATES_Lk_Linreg(tree); */
  return -1.;
  /* return tree->rates->c_lnL_linreg; */
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Wrap_Diff_Lk_Norm_At_Given_Edge(t_edge *b, t_tree *tree, supert_tree *stree)
{
  phydbl diff;
  diff = Diff_Lk_Norm_At_Given_Edge(b,tree);
  return(-diff);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Sample_Ancestral_Seq(int mutmap, int fromprior, t_tree *tree)
{
  int rate_cat;
  int i,j,k,l;
  phydbl *probs;
  phydbl sum;
  phydbl u;
  int n_mut;
  FILE *fp;
  phydbl *muttime;
  int *muttype;
  char *s;
  int *ordering;

  probs = (phydbl *)mCalloc(tree->mod->ras->n_catg,sizeof(phydbl));

  muttime = (phydbl *)mCalloc((2*tree->n_otu-3)*5, // At most 5 mutations per branch on average
               sizeof(phydbl));

  muttype = (int *)mCalloc((2*tree->n_otu-3)*5, // At most 5 mutations per branch on average
               sizeof(int));

  ordering = (int *)mCalloc((2*tree->n_otu-3)*5, // At most 5 mutations per branch on average
                sizeof(int));

  s = (char *)mCalloc(T_MAX_NAME,sizeof(char));

  For(i,2*tree->n_otu-1) 
    if(tree->a_nodes[i]->tax == NO)
      {
        tree->a_nodes[i]->c_seq_anc = (align *)mCalloc(1,sizeof(align));;
        tree->a_nodes[i]->c_seq_anc->state = (char *)mCalloc(tree->n_pattern,sizeof(char));
      }

  if(fromprior == YES)
    {
      /* Update P(D_x|X=i) for each state i and node X */
      Set_Both_Sides(YES,tree);
      Lk(NULL,tree);
    }

  for(i=0;i<tree->n_pattern;i++)
    {
      /* Sample the rate class from its posterior density */
      for(j=0;j<tree->mod->ras->n_catg;j++)
        {
          if(fromprior == NO)
            probs[j] =
              tree->unscaled_site_lk_cat[i*tree->mod->ras->n_catg+j]*
              tree->mod->ras->gamma_r_proba->v[j];
          else
            probs[j] = tree->mod->ras->gamma_r_proba->v[j];
        }


      /* Scale probas. */
      sum = .0;
      for(j=0;j<tree->mod->ras->n_catg;j++) sum += probs[j];
      for(j=0;j<tree->mod->ras->n_catg;j++) probs[j]/=sum;

      /* CDF */
      for(j=1;j<tree->mod->ras->n_catg;j++) probs[j] += probs[j-1];

      /* Sample rate */
      u = Uni();
      rate_cat = -1;
      for(j=0;j<tree->mod->ras->n_catg;j++)
        if(probs[j] > u)
          {
            rate_cat = j;
            break;
          }
      
      n_mut = 0;
      Sample_Ancestral_Seq_Pre(tree->a_nodes[0],tree->a_nodes[0]->v[0],tree->a_nodes[0]->b[0],
                               i,rate_cat,
                               muttype,muttime,&n_mut,
                               mutmap,fromprior,tree);
      

      for(j=0;j<n_mut;j++) ordering[j] = 0;
      
      for(j=0;j<n_mut-1;j++)
        {
          for(k=j+1;k<n_mut;k++)
            {
              if(muttime[k] > muttime[j]) ordering[k]++;
              else ordering[j]++;
            }
        }
      
      strcpy(s,"rosetta.");
      sprintf(s+strlen(s),"%d",i);
      fp = fopen(s,"a");
      PhyML_Fprintf(fp,"\n-1 -1 -1.0 -1");

      for(j=0;j<n_mut;j++)
    {
      for(k=0;k<n_mut;k++)
        {
          if(ordering[k] == j)
        {
          for(l=0;l<tree->data->init_len;l++) if(tree->data->sitepatt[l] == i) break;
          PhyML_Fprintf(fp,"\n%4d %4d %12f %4d",j,muttype[k],muttime[k],l);
          /* PhyML_Fprintf(fp,"\n%d",muttype[ordering[j]]); */
          break;
        }
        }
    }


      for(j=0;j<n_mut;j++)
        {
          muttype[j] = -2;
          muttime[j] = +1.;
        }
      
      fclose(fp);
    }

  Free(s);
  Free(muttype);
  Free(muttime);
  Free(ordering);
  Free(probs);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Sample_Ancestral_Seq_Pre(t_node *a, t_node *d, t_edge *b,
                              int site, int rate_cat,
                              int *muttype, phydbl *muttime, int *n_mut,
                              int mutmap, int fromprior, t_tree *tree)
{

  int i,j;
  int sa,sd;
  phydbl *Pij;
  phydbl *p_lk;
  int dim1, dim2, dim3;
  phydbl sum;
  phydbl u;
  char *c;
  phydbl *probs;

  probs = (phydbl *)mCalloc(tree->mod->ns,sizeof(phydbl));

  if(a->tax)
    c = a->c_seq->state+site*tree->mod->io->state_len;
  /* c = tree->data->c_seq[a->num]->state+site*tree->mod->io->state_len; */
  else
    c = a->c_seq_anc->state+site*tree->mod->io->state_len;
  /* c = tree->anc_data->c_seq[a->num-tree->n_otu]->state+site*tree->mod->io->state_len; */
  
  sa = Assign_State(c,
                    tree->mod->io->datatype,
                    tree->mod->io->state_len);
  
  if(sa == -1) /* c is an indel */
    {
      for(j=0;j<tree->mod->ns;j++) probs[j] = tree->mod->e_frq->pi->v[j];
      
      for(j=1;j<tree->mod->ns;j++) probs[j] += probs[j-1];
      
      u = Uni();
      for(j=0;j<tree->mod->ns;j++)
        if(probs[j] > u)
          {
            sa = j;
            break;
          }
    }
  
  if(d->tax == NO) // Need to sample state at node d
    {
      
      dim1 = tree->mod->ras->n_catg * tree->mod->ns;
      dim2 = tree->mod->ns;
      dim3 = tree->mod->ns * tree->mod->ns;
      sum  = 0.0;
            
      Pij  = b->Pij_rr;
      
      if(d == b->left)
        p_lk = b->p_lk_left;
      else
        p_lk = b->p_lk_rght;
      
      for(j=0;j<tree->mod->ns;j++) probs[j] = 0.0;
      
      /* Formula (10) in Nielsen's Mutation Maping paper, e.g. */
      for(j=0;j<tree->mod->ns;j++)
        {
          if(fromprior == NO)
            probs[j] =
              p_lk[site*dim1+rate_cat*dim2+j] *
              Pij[rate_cat*dim3+sa*dim2+j];
          else
            probs[j] = Pij[rate_cat*dim3+sa*dim2+j];
        }
      
      /* Scale the probabilities */
      sum = 0.0;
      for(j=0;j<tree->mod->ns;j++) sum += probs[j];
      for(j=0;j<tree->mod->ns;j++) probs[j] /= sum;
      
      /* CDF */
      for(j=1;j<tree->mod->ns;j++) probs[j] += probs[j-1];
      
      /* Sample state according to their posterior probas. */
      sd = -1;
      u = Uni();
      for(j=0;j<tree->mod->ns;j++)
        if(probs[j] > u)
          {
            sd = j;
            break;
          }
      
      /* Assign state */
      /* tree->anc_data->c_seq[d->num-tree->n_otu]->state[site] = Reciproc_Assign_State(sd,tree->io->datatype); */
      /* printf("\n<> %p",d->c_seq_anc); fflush(NULL); */
      d->c_seq_anc->state[site] = Reciproc_Assign_State(sd,tree->io->datatype);
    }
  else
    {
      c = d->c_seq->state+site*tree->mod->io->state_len;
      /* c = tree->data->c_seq[d->num]->state+site*tree->mod->io->state_len; */
      
      sd = Assign_State(c,
                        tree->mod->io->datatype,
                        tree->mod->io->state_len);
      
      if(sd == -1) // c is an indel
        {
          for(j=0;j<tree->mod->ns;j++) probs[j] = tree->mod->e_frq->pi->v[j];
          
          for(j=1;j<tree->mod->ns;j++) probs[j] += probs[j-1];
          
          u = Uni();
          for(j=0;j<tree->mod->ns;j++)
            if(probs[j] > u)
              {
                sd = j;
                break;
              }
        }
    }
  
  /* if(site == 92) */
  /*   { */
  /*     printf("\n. sa=%d (%s,%d) sd=%d (%s,%d) b->l->v=%f", */
  /* 	     sa,a->tax?a->name:"",a->num,sd,d->tax?d->name:"",d->num,b->l->v); */
  /*   } */
  
  if(mutmap == YES) Map_Mutations(a,d,sa,sd,b,site,rate_cat,muttype,muttime,n_mut,tree);
  
  Free(probs);
  
  if(d->tax) return;
  else
    {
      for(i=0;i<3;i++)
        {
          if(d->v[i] != a)
            {
              Sample_Ancestral_Seq_Pre(d,d->v[i],d->b[i],site,rate_cat,muttype,muttime,n_mut,mutmap,fromprior,tree);
            }
        }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Map_Mutations(t_node *a, t_node *d, int sa, int sd, t_edge *b, int site, int rate_cat, int *muttype, phydbl *muttime, int *n_mut, t_tree *tree)
{
  int i,j;
  phydbl *probs,*all_probs;
  int slast; // Last state visited
  phydbl tlast, td;
  phydbl *Q;
  phydbl u,sum;
  int *mut; // Array of mutations
  int thismut;
  int n_mut_branch;
  phydbl br,cr,ta,gr;
  int n_iter;
  int first_mut;

  all_probs = (phydbl *)mCalloc(tree->mod->ns*tree->mod->ns,sizeof(phydbl));
  mut = (int *)mCalloc(tree->mod->ns*tree->mod->ns,sizeof(int));

  // Edge rate
  br = tree->rates->br_r[d->num];

  // Clock (i.e., overall) rate
  cr = tree->rates->clock_r;

  // Age of node a
  ta = tree->rates->nd_t[a->num];

  // Site relative rate
  gr = tree->mod->ras->gamma_rr->v[rate_cat];


  // Rate matrix
  Q = tree->mod->r_mat->qmat->v;

  // Length of the 'time' interval considered: product of the branch length by the
  // current relative rate at that site (set when sampling ancestral sequences)
  td = b->l->v * gr;

  // Matrix of change probabilities
  for(i=0;i<tree->mod->ns;i++)
    {
      // We only care about the non-diagonal elements here
      for(j=0;j<tree->mod->ns;j++) all_probs[i*tree->mod->ns+j] = -Q[i*tree->mod->ns+j] / Q[i*tree->mod->ns+i];

      // Set the diagonal to 0
      all_probs[i*tree->mod->ns+i] = 0.0;

      // \sum_{j != i} -q_{ij}/q_{ii}
      sum = 0;
      for(j=0;j<tree->mod->ns;j++) sum += all_probs[i*tree->mod->ns+j];

      // Diagonal: 1 - \sum_{j != i} -q_{ij}/q_{ii}
      all_probs[i*tree->mod->ns+i] = 1.-sum;

      // Get the cumulative probas
      for(j=1;j<tree->mod->ns;j++) all_probs[i*tree->mod->ns+j] += all_probs[i*tree->mod->ns+j-1];
    }

  For(i,tree->mod->ns*tree->mod->ns) mut[i] = 0;
  tlast = .0;
  slast = sa;
  probs = NULL;
  n_mut_branch = 0;
  n_iter = 0;
  first_mut = YES;

  do
    {
      if((sa != sd) && (first_mut == YES)) // ancestral and descendant states are distinct
    {
      // Sample a time for the first mutation conditional on at least one mutation
      // occurring (see formula A2 in Nielsen's Genetics paper (2001))
      u = Uni();
      tlast = -log(1. - u*(1.-exp(Q[sa*tree->mod->ns+sa]*td)))/-Q[sa*tree->mod->ns+sa];
    }
      else
    {
      // Sample a time for the next mutation
      tlast = tlast + Rexp(-Q[slast*tree->mod->ns+slast]);
    }

      // Select the appropriate vector of change probabilities
      probs = all_probs+slast*tree->mod->ns;

      /* printf("\n. slast=%2d sd=%2d tlast=%12G td=%12G p=%12f rcat=%12f site=%4d", */
      /* 	 slast,sd,tlast,td,-Q[slast*tree->mod->ns+slast], */
      /* 	 tree->mod->ras->gamma_rr->v[rate_cat],site); */

      // The time for the next mutation does not exceed the length
      // of the time interval -> sample a new mutation event
      if(tlast < td)
    {
      first_mut = NO;

      n_mut_branch++;

      u = Uni();
      for(i=0;i<tree->mod->ns;i++)
        if(probs[i] > u)
          {
        // Record mutation type
        mut[slast*tree->mod->ns+i]++;

        // Record mutation type in the site mutation array
        thismut = MIN(i,slast) * tree->mod->ns + MAX(i,slast) - (MIN(i,slast)+1+(int)POW(MIN(i,slast)+1,2))/2;
        muttype[(*n_mut)+n_mut_branch-1] = thismut;

        if((thismut > 5) || (thismut < 0))
          {
            PhyML_Fprintf(stderr,"\n. thismut = %d",thismut);
            PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d\n\n",__FILE__,__LINE__);
            Warn_And_Exit("");
          }

        // Record time of mutation
        muttime[(*n_mut)+n_mut_branch-1] = ta + br*cr*gr;

        // Update the last state
        slast = i;
        break;
          }
    }
      else
    {
      if(slast == sd) break;
      else
        {
          // Restart from the beginning
          For(i,tree->mod->ns*tree->mod->ns) mut[i] = 0;
          for(i=0;i<n_mut_branch;i++) muttype[(*n_mut)+n_mut_branch-1] = -2;
          for(i=0;i<n_mut_branch;i++) muttime[(*n_mut)+n_mut_branch-1] = +1.;
          tlast = 0.0;
          slast = sa;
          n_mut_branch = 0;
          first_mut = YES;
          n_iter++;
        }
    }
    }
  while(1);

  (*n_mut) += n_mut_branch;


  for(i=0;i<tree->mod->ns;i++)
    {
      for(j=i+1;j<tree->mod->ns;j++)
    {
      if(mut[i*tree->mod->ns+j] + mut[j*tree->mod->ns+i] > 0)
        {
          thismut = MIN(i,j) * tree->mod->ns + MAX(i,j) - (MIN(i,j)+1+(int)POW(MIN(i,j)+1,2))/2;
          tree->mutmap[thismut*(tree->n_pattern)*(2*tree->n_otu-3) + b->num*(tree->n_pattern) + site]++;
          /* if(site == 92) */
          /* 	{ */
          /* 	  printf("\nx sa=%d sd=%d mut=%d",sa,sd,thismut); */
          /* 	} */
        }
    }
    }

  Free(all_probs);
  Free(mut);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int Check_Lk_At_Given_Edge(int verbose, t_tree *tree)
{
  int res;
  int i;
  phydbl *lk;

  lk = (phydbl *)mCalloc(2*tree->n_otu-3,sizeof(phydbl));

  res = 0;
  For(i,2*tree->n_otu-3)
    {
      lk[i] = Lk(tree->a_edges[i],tree);
      if(verbose == YES) PhyML_Printf("\n. Edge %3d %13G %f %13G",
                                      tree->a_edges[i]->num,tree->a_edges[i]->l->v,lk[i],
                                      tree->a_edges[i]->l_var->v);
            
    }

  if(tree->n_root && tree->ignore_root == NO)
    {
      Lk(tree->n_root->b[1],tree);
      if(verbose == YES) PhyML_Printf("\nx Edge %3d %13G %f %13G",
                                      tree->n_root->b[1]->num,tree->n_root->b[1]->l->v,tree->c_lnL,
                                      tree->n_root->b[1]->l_var->v
                                      );

      Lk(tree->n_root->b[2],tree);
      if(verbose == YES) PhyML_Printf("\nx Edge %3d %13G %f %13G",
                                      tree->n_root->b[2]->num,tree->n_root->b[2]->l->v,tree->c_lnL,
                                      tree->n_root->b[2]->l_var->v
                                      );

    }

  res=1;
  for(i=1;i<2*tree->n_otu-3;i++)
    {
      if(FABS(lk[i]-lk[i-1]) > 1.E-2) res=0;
    }
  Free(lk);

  return res;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Ancestral_Sequences(t_tree *tree, int print)
{
  int i;

  if(print == YES)
    {
      PhyML_Printf("\n\n. Estimating ancestral sequences...");

      strcpy(tree->io->out_ancestral_file,tree->io->out_file);
      if(tree->io->append_run_ID) { strcat(tree->io->out_ancestral_file,"_"); strcat(tree->io->out_ancestral_file,tree->io->run_id_string); }
      strcat(tree->io->out_ancestral_file,"_phyml_ancestral_seq.txt");
      tree->io->fp_out_ancestral = Openfile(tree->io->out_ancestral_file,1);
      

      char *s = (char *)mCalloc((int)strlen(tree->io->in_align_file)+50,sizeof(char));
      strcpy(s,tree->io->in_align_file);
      strcat(s,"_phyml_ancestral_tree.txt");
      FILE *fp = Openfile(s,WRITE);
      
      PhyML_Fprintf(tree->io->fp_out_ancestral,"\n\n\n");
      PhyML_Fprintf(tree->io->fp_out_ancestral,"\n. Printing marginal probabilities of ancestral sequences at each site");
      PhyML_Fprintf(tree->io->fp_out_ancestral,"\n. of the alignment and each node of the tree. The tree in Newick format");
      PhyML_Fprintf(tree->io->fp_out_ancestral,"\n. with internal nodes labels corresponding to those given below can be");
      PhyML_Fprintf(tree->io->fp_out_ancestral,"\n. found in the file '%s'.",s);
      Free(s);
      PhyML_Fprintf(tree->io->fp_out_ancestral,"\n\n");
      PhyML_Fprintf(tree->io->fp_out_ancestral,"Site\tNodeLabel\t");
      for(i=0;i<tree->mod->ns;i++) PhyML_Fprintf(tree->io->fp_out_ancestral,"%c\t",Reciproc_Assign_State(i,tree->io->datatype));
      PhyML_Fprintf(tree->io->fp_out_ancestral,"\n");

      
      short int bck_boot_val = tree->print_boot_val;
      short int bck_alrt_val = tree->print_alrt_val;
      
      tree->print_node_num = YES;
      tree->print_boot_val = NO;
      tree->print_alrt_val = NO;
      s = Write_Tree(tree,NO);
      PhyML_Fprintf(fp,"%s",s);
      tree->print_node_num = NO;
      tree->print_boot_val = bck_boot_val;
      tree->print_alrt_val = bck_alrt_val;

      Free(s);
      fclose(fp);
    }

  for(i=0;i<2*tree->n_otu-2;i++)
    if(tree->a_nodes[i]->tax == NO)
      Ancestral_Sequences_One_Node(tree->a_nodes[i],tree,print);

  if(tree->n_root) Ancestral_Sequences_One_Node(tree->n_root,tree,print);


  fclose(tree->io->fp_out_ancestral);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Ancestral_Sequences_One_Node(t_node *d, t_tree *tree, int print)
{
  if(d->tax) return;
  else
    {
      if(tree->is_mixt_tree) 
        {
          MIXT_Ancestral_Sequences_One_Node(d,tree,print);
        }
      else
        {
          t_node *v0,*v1,*v2; // three neighbours of d
          t_edge *b0,*b1,*b2;
          int i,j;
          int catg;
          phydbl p0, p1, p2;
          phydbl *p;
          int site,csite;
          phydbl *p_lk0, *p_lk1, *p_lk2;
          int *sum_scale0, *sum_scale1, *sum_scale2;
          phydbl sum_probas;
          phydbl *Pij0, *Pij1, *Pij2;          
          phydbl inc,sum_scale;
          FILE *fp;

          unsigned const int ncatg = tree->mod->ras->n_catg;
          unsigned const int ns = tree->mod->ns;
          unsigned const int nsns = ns*ns;
          unsigned const int ncatgns = ns*ncatg;


          if(tree->scaling_method == SCALE_RATE_SPECIFIC)
            {
              PhyML_Fprintf(stderr,"\n. Likelihood rescaling method not compatible with the calculation of ancestral state probabilities.");
              Exit("\n");
            }

          
          if(!d) return;
          
          fp = tree->io->fp_out_ancestral;
          assert(fp != NULL);

          
          p = (phydbl *)mCalloc(ns,sizeof(phydbl));
              
          for(site=0;site<tree->data->init_len;site++) // For each site in the current partition element
            {
              csite = tree->data->sitepatt[site];
                                    
              for(i=0;i<ns;i++) p[i] = .0;
                  
              v0 = d->v[0];
              v1 = d->v[1];
              v2 = d->v[2];
              
              b0 = d->b[0];
              b1 = d->b[1];
              b2 = d->b[2];
              
              Pij0 = b0->Pij_rr;
              Pij1 = b1->Pij_rr;
              Pij2 = b2->Pij_rr;

              sum_scale = 0.0;
              
              if(v0 == b0->left)
                {
                  p_lk0 = b0->p_lk_left;
                  sum_scale0 = b0->sum_scale_left;
                }
              else
                {
                  p_lk0 = b0->p_lk_rght;
                  sum_scale0 = b0->sum_scale_rght;
                }
              
              if(v1 == b1->left)
                {
                  p_lk1 = b1->p_lk_left;
                  sum_scale1 = b1->sum_scale_left;
                }
              else
                {
                  p_lk1 = b1->p_lk_rght;
                  sum_scale1 = b1->sum_scale_rght;
                }
              
              if(v2 == b2->left)
                {
                  p_lk2 = b2->p_lk_left;
                  sum_scale2 = b2->sum_scale_left;
                }
              else
                {
                  p_lk2 = b2->p_lk_rght;
                  sum_scale2 = b2->sum_scale_rght;
                }
              

              for(catg=0;catg<ncatg;++catg)
                {
                  for(i=0;i<ns;++i)
                    {
                      p0 = .0;
                      if(v0->tax)
                        {
                          for(j=0;j<ns;++j)
                            {
                              p0 += v0->b[0]->p_lk_tip_r[csite*ns+j] * Pij0[catg*nsns+i*ns+j];
                              /* if(isinf(p0) || isnan(p0) || p0 < SMALL) */
                              if(isinf(p0) || isnan(p0)) 
                               {
                                  PhyML_Fprintf(stderr,"\n. p0: %G v0->b[0]->p_lk_tip_r[csite*ns+j]: %G Pij0[catg*nsns+i*ns+j]: %G\n",
                                                p0,
                                                v0->b[0]->p_lk_tip_r[csite*ns+j],
                                                Pij0[catg*nsns+i*ns+j]);
                                  Exit("\n");
                                }
                            }
                        }
                      else
                        {
                          for(j=0;j<ns;j++)
                            {
                              /* p0 += p_lk0[csite*ncatgns+catg*ns+j] * Pij0[catg*nsns+i*ns+j] / (phydbl)POW(2,sum_scale0[csite]); */
                              p0 += p_lk0[csite*ncatgns+catg*ns+j] * Pij0[catg*nsns+i*ns+j];
                              if(isinf(p0) || isnan(p0))
                                {
                                  PhyML_Fprintf(stderr,"\n. p0: %G p_lk0[csite*ncatgns+catg*ns+j]: %G Pij0[catg*nsns+i*ns+j]: %G (phydbl)POW(2,sum_scale0[csite*ncatg+catg]): %G sum_scale0: %d\n",
                                                p0,
                                                p_lk0[csite*ncatgns+catg*ns+j],
                                                Pij0[catg*nsns+i*ns+j],
                                                (phydbl)POW(2,sum_scale0[csite]),
                                                sum_scale0[csite]);
                                  Exit("\n");
                                }
                            }
                          if(catg == 0 && i == 0) sum_scale += sum_scale0[csite];
                        }
                      
                      p1 = .0;
                      if(v1->tax)
                        {
                          for(j=0;j<ns;j++)
                            {
                              p1 += v1->b[0]->p_lk_tip_r[csite*ns+j] * Pij1[catg*nsns+i*ns+j];
                              if(isinf(p1) || isnan(p1))
                                {
                                  PhyML_Fprintf(stderr,"\n. p1: %G v1->b[0]->p_lk_tip_r[csite*ns+j]: %G Pij1[catg*nsns+i*ns+j]: %G\n",
                                                p1,
                                                v1->b[0]->p_lk_tip_r[csite*ns+j],
                                                Pij1[catg*nsns+i*ns+j]);
                                  Exit("\n");
                                }
                            }
                        }
                      else
                        {
                          for(j=0;j<ns;j++)
                            {
                              p1 += p_lk1[csite*ncatgns+catg*ns+j] * Pij1[catg*nsns+i*ns+j];
                              if(isinf(p1) || isnan(p1))
                                {
                                  PhyML_Fprintf(stderr,"\n. p1: %G p_lk1[csite*ncatgns+catg*ns+j]: %G Pij1[catg*nsns+i*ns+j]: %G (phydbl)POW(2,sum_scale1[csite*ncatg+catg]): %G\n",
                                                p1,
                                                p_lk1[csite*ncatgns+catg*ns+j],
                                                Pij1[catg*nsns+i*ns+j],
                                                (phydbl)POW(2,sum_scale1[csite]));
                                  Exit("\n");
                                }
                            }
                          if(catg == 0 && i == 0) sum_scale += sum_scale1[csite];
                        }
                      
                      p2 = .0;
                      if(v2->tax)
                        {
                          for(j=0;j<ns;j++)
                            {
                              p2 += v2->b[0]->p_lk_tip_r[csite*ns+j] * Pij2[catg*nsns+i*ns+j];
                            }
                        }
                      else
                        {
                          for(j=0;j<ns;j++)
                            {
                              p2 += p_lk2[csite*ncatgns+catg*ns+j] * Pij2[catg*nsns+i*ns+j];
                              if(isinf(p2) || isnan(p2))
                                {
                                  PhyML_Fprintf(stderr,"\n. p2: %G p_lk2[csite*ncatgns+catg*ns+j]: %G Pij2[catg*nsns+i*ns+j]: %G (phydbl)POW(2,sum_scale2[csite]): %G\n",
                                                p2,
                                                p_lk2[csite*ncatgns+catg*ns+j],
                                                Pij2[catg*nsns+i*ns+j],
                                                (phydbl)POW(2,sum_scale2[csite]));
                                  Exit("\n");
                                }
                            }
                          if(catg == 0 && i == 0) sum_scale += sum_scale2[csite];
                        }
                          
                      inc =
                        p0*p1*p2*
                        tree->mod->e_frq->pi->v[i] *
                        tree->mod->ras->gamma_r_proba->v[catg];

                      p[i] += inc;                      

                      
                      if(isinf(p[i]) || isnan(p[i]))
                        {
                          PhyML_Fprintf(stderr,"\n. site: %4d p0: %G p1: %G p2: %G tree->mod->e_frq->pi->v[i]: %G tree->cur_site_lk[csite]: %G  tree->mod->ras->gamma_r_proba->v[catg]: %G tree->c_lnL_sorted[csite]: %G",
                                       csite,
                                       p0,p1,p2,
                                       tree->mod->e_frq->pi->v[i] ,
                                       tree->cur_site_lk[csite] ,
                                       tree->mod->ras->gamma_r_proba->v[catg],
                                       tree->c_lnL_sorted[csite]);
                          Exit("\n");
                        }
                    }

                  /* printf("\n. site: %d || %d %d %d", */
                  /*        csite, */
                  /*        v0->tax ? -1 : sum_scale0[csite*ncatg+catg], */
                  /*        v1->tax ? -1 : sum_scale1[csite*ncatg+catg], */
                  /*        v2->tax ? -1 : sum_scale2[csite*ncatg+catg]); */

                }

              for(i=0;i<ns;i++) p[i] = log(p[i]) - (phydbl)LOG2 * sum_scale;
              for(i=0;i<ns;i++) p[i] -= tree->c_lnL_sorted[csite];
              for(i=0;i<ns;i++) p[i] = exp(p[i]);
              
              /* sum_probas = 0.0; */
              /* for(i=0;i<ns;i++) sum_probas += p[i]; */
              /* for(i=0;i<ns;i++) p[i]/=sum_probas; */
              
              
              if(print == YES)
                {
                  PhyML_Fprintf(fp,"%4d\t%9d\t",site+1,d->num);
                  sum_probas = .0;
                  for(i=0;i<ns;i++)
                    {
                      PhyML_Fprintf(fp,"%.4f\t",p[i]);
                      sum_probas += p[i];
                    }
                  PhyML_Fprintf(fp,"\n");
                  fflush(NULL);
                  if(Are_Equal(sum_probas,1.0,0.01) == NO)
                    {
                      PhyML_Fprintf(stderr,"\n. Probabilities do not sum to 1.0! Aborting.");
                      for(i=0;i<ns;++i) PhyML_Fprintf(stderr,"\n. p[%2d]=%G",i,p[i]);
                      Exit("\n");
                    }
                }
            }
          Free(p);
        }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// Computes the value of fact_sum_scale, the part of the scaling factors
// that is common to all classes of the mixture and scale the 
// likelihood for each mixture using the part of the scaling
// factors that is class-specific
 
void Pull_Scaling_Factors(int site, t_edge *b, t_tree *tree)
{
  unsigned int catg;
  const unsigned int ncatg = tree->mod->ras->n_catg;

  if(tree->apply_lk_scaling == NO)
    {
      tree->fact_sum_scale[site] = 0;
      for(catg=0;catg<ncatg;++catg) tree->unscaled_site_lk_cat[site*ncatg+catg] = tree->site_lk_cat[catg];
      return;
    }
  else
    {
      switch(tree->scaling_method)
        {
        case SCALE_RATE_SPECIFIC : 
          {
            int *sum_scale_left_cat,*sum_scale_rght_cat;
            int exponent;
            phydbl max_sum_scale,min_sum_scale;
            phydbl sum,tmp,dum;
            
            sum_scale_left_cat = b->sum_scale_left_cat;
            sum_scale_rght_cat = b->sum_scale_rght_cat;
            
            max_sum_scale =   (phydbl)BIG;
            min_sum_scale =  -(phydbl)BIG;
            
            for(catg=0;catg<ncatg;++catg)
              {
                sum_scale_left_cat[catg] =
                  (b->sum_scale_left)?
                  (b->sum_scale_left[site*ncatg+catg]):
                  (0.0);
                
                sum_scale_rght_cat[catg] =
                  (b->sum_scale_rght)?
                  (b->sum_scale_rght[site*ncatg+catg]):
                  (0.0);
                
                sum = sum_scale_left_cat[catg] + sum_scale_rght_cat[catg];
                
                if(sum < .0)
                  {
                    PhyML_Fprintf(stderr,"\n. tree: %s\n",Write_Tree(tree,NO));
                    PhyML_Fprintf(stderr,"\n. b->num = %d  sum = %G root ? %d",sum,b->num,b == tree->e_root);
                    PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d.\n",__FILE__,__LINE__);
                    Exit("\n");
                  }
                
                dum = log(FABS(tree->site_lk_cat[catg]));
                
                tmp = sum + ((phydbl)LOGBIG - dum)/(phydbl)LOG2;
                if(tmp < max_sum_scale) max_sum_scale = tmp; /* min of the maxs */
                
                tmp = sum + ((phydbl)LOGSMALL - dum)/(phydbl)LOG2;
                if(tmp > min_sum_scale) min_sum_scale = tmp; /* max of the mins */
                
                assert(isnan(tmp) == NO);
              }
            
            if(min_sum_scale > max_sum_scale)
              {
#ifdef SAFEMODE
                PhyML_Printf("\n. Numerical precision issue alert.");
                PhyML_Printf("\n. min_sum_scale = %G max_sum_scale = %G",min_sum_scale,max_sum_scale);
#endif
                min_sum_scale = max_sum_scale;
              }
            
            tree->fact_sum_scale[site] = (int)((max_sum_scale + min_sum_scale) / 2);
                        
            /* Apply scaling factors */
            for(catg=0;catg<ncatg;++catg)
              {
                exponent = -(sum_scale_left_cat[catg]+sum_scale_rght_cat[catg])+tree->fact_sum_scale[site];
                Rate_Correction(exponent,tree->site_lk_cat + catg);
              }
            
            break;
          }
        case SCALE_FAST :
          {
            int sum_scale_left,sum_scale_rght;
            
            sum_scale_left =
              (b->sum_scale_left)?
              (b->sum_scale_left[site]):
              (0.0);
            
            sum_scale_rght =
              (b->sum_scale_rght)?
              (b->sum_scale_rght[site]):
              (0.0);
            
            tree->fact_sum_scale[site] = sum_scale_left + sum_scale_rght;

            break;
          }
        }
      for(catg=0;catg<ncatg;++catg) tree->unscaled_site_lk_cat[site*ncatg+catg] = tree->site_lk_cat[catg]; 
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

// Tree should be ready for likelihood analysis when calling
// this function.
 void Stepwise_Add_Lk(t_tree *tree)
 {
   t_edge **residuals,**targets,*best_target;
  int *nd_idx,i,j,n_targets,*tg_idx,n_opt;

  residuals   = (t_edge **)mCalloc(tree->n_otu-3,sizeof(t_edge *));
  targets     = (t_edge **)mCalloc(2*tree->n_otu-3,sizeof(t_edge *));
  best_target = NULL;
  nd_idx      = Permutate(tree->n_otu-3);

  // Remove all tips except that corresponding to a_nodes[0], 
  // a_nodes[1] and a_nodes[2].  
  for(i=0;i<tree->n_otu-3;i++)
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
  for(i=0;i<tree->n_otu-3;i++)
    {
      Set_Both_Sides(YES,tree);
      Lk(NULL,tree);

      printf("\n. [%d/%d]",i,tree->n_otu-3);

      tree->best_lnL = UNLIKELY;
      best_target    = NULL;
      tg_idx         = Permutate(n_targets);

      for(j=0;j<n_targets;j++)
        {
          Graft_Subtree(targets[tg_idx[j]],
                        tree->a_nodes[nd_idx[i]+3]->v[0],
                        NULL,
                        residuals[i],
                        NULL,
                        tree);
          
          Update_PMat_At_Given_Edge(targets[tg_idx[j]],tree);
          Update_PMat_At_Given_Edge(tree->a_nodes[nd_idx[i]+3]->b[0],tree);
          Update_Partial_Lk(tree,residuals[i],tree->a_nodes[nd_idx[i]+3]->v[0]);
          Lk(residuals[i],tree);

          if(tree->c_lnL > tree->best_lnL)
            {
              tree->best_lnL = tree->c_lnL;
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
      
      n_opt = 0;
      do Optimize_Br_Len_Serie (tree); while(n_opt++ < 3);

      targets[n_targets]   = residuals[i]; 
      targets[n_targets+1] = tree->a_nodes[nd_idx[i]+3]->b[0];
      
      Free(tg_idx);
      n_targets+=2;
    }

  Round_Optimize(tree,5);
  PhyML_Fprintf(stderr,"\n. lk: %f",tree->c_lnL);
  Exit("\n");

  Free(nd_idx);
  Free(residuals);
  Free(targets);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


/*     |
       |
       |b
       |
       |d
      /  \
     /    \
    /      \
   /        \
  /v1        \v2

  Set p_lk and sum_scale for subtrees with d, v1 and v2 as root,
  Pij for edges b, and the two edges connecting d to v1 and d to
  v2;
  Account for rooted trees.
*/

void Set_All_Partial_Lk(t_node **n_v1, t_node **n_v2,
                        phydbl **p_lk, int **sum_scale, int **p_lk_loc,
                        phydbl **Pij1, phydbl **tPij1, phydbl **p_lk_v1, int **sum_scale_v1,
                        phydbl **Pij2, phydbl **tPij2, phydbl **p_lk_v2, int **sum_scale_v2,
                        t_node *d, t_edge *b, t_tree *tree
#ifdef BEAGLE
                        , int *dest_p_idx, int *child1_p_idx, int* child2_p_idx, int* Pij1_idx, int* Pij2_idx
#endif
                        )
{
  unsigned int i;
  
  assert(tree->is_mixt_tree == NO);
  assert(d->tax == NO);
  
  if(tree->n_root == NULL || tree->ignore_root == YES)
    {
      /* Does d lie on the "left" or "right" of the branch? */
      if(d == b->left)
        {
          *p_lk      = b->p_lk_left;
          *sum_scale = b->sum_scale_left;
          *p_lk_loc  = b->p_lk_loc_left;
#ifdef BEAGLE
          *dest_p_idx = b->p_lk_left_idx;
#endif
        }
      else
        {
          *p_lk      = b->p_lk_rght;
          *sum_scale = b->sum_scale_rght;
          *p_lk_loc  = b->p_lk_loc_rght;
#ifdef BEAGLE
          *dest_p_idx = b->p_lk_rght_idx;
#endif
        }

      *n_v1 = *n_v2 = NULL;
      for(i=0;i<3;++i)
        {
          if(d->b[i] != b)
            {
              if(!(*n_v1))
                {
                  *n_v1 = d->v[i];
#ifdef BEAGLE
                  Set_Partial_Lk_One_Side(Pij1,tPij1,p_lk_v1,sum_scale_v1,d,d->b[i],tree,child1_p_idx,Pij1_idx);
#else
                  Set_Partial_Lk_One_Side(Pij1,tPij1,p_lk_v1,sum_scale_v1,d,d->b[i],tree);
#endif
                }
              else if(!(*n_v2))
                {
                  *n_v2 = d->v[i];
#ifdef BEAGLE
                  Set_Partial_Lk_One_Side(Pij2,tPij2,p_lk_v2,sum_scale_v2,d,d->b[i],tree,child2_p_idx,Pij2_idx);
#else
                  Set_Partial_Lk_One_Side(Pij2,tPij2,p_lk_v2,sum_scale_v2,d,d->b[i],tree);
#endif
                }
              else assert(FALSE);
            }
        }
    }
  else
    {
      if(b == tree->e_root)
        {
          if(d == tree->n_root->v[1])      b = tree->n_root->b[1];
          else if(d == tree->n_root->v[2]) b = tree->n_root->b[2];
          else assert(FALSE);
        }

      if(d == tree->n_root)
        {
          if(b == tree->n_root->b[1])
            {
              *p_lk      = tree->n_root->b[1]->p_lk_left;
              *sum_scale = tree->n_root->b[1]->sum_scale_left;
              *p_lk_loc  = tree->n_root->b[1]->p_lk_loc_left;
#ifdef BEAGLE
              *dest_p_idx = tree->n_root->b[1]->p_lk_left_idx;
#endif
            }
          else
            {
              *p_lk      = tree->n_root->b[2]->p_lk_left;
              *sum_scale = tree->n_root->b[2]->sum_scale_left;
              *p_lk_loc  = tree->n_root->b[2]->p_lk_loc_left;
#ifdef BEAGLE
              *dest_p_idx = tree->n_root->b[2]->p_lk_left_idx;
#endif
            }

          *n_v1         = NULL;
          *Pij1         = NULL;
          *tPij1        = NULL;
          *p_lk_v1      = NULL;
          *sum_scale_v1 = NULL;

          if(b == tree->n_root->b[1])
            {
              *n_v2         = tree->n_root->v[2];
              *Pij2         = tree->n_root->b[2]->Pij_rr;
              *tPij2        = tree->n_root->b[2]->tPij_rr;
              *p_lk_v2      = tree->n_root->b[2]->p_lk_rght;
              *sum_scale_v2 = tree->n_root->b[2]->sum_scale_rght;
#ifdef BEAGLE
              *child2_p_idx = tree->n_root->b[2]->p_lk_rght_idx;
              *Pij2_idx     = tree->n_root->b[2]->Pij_rr_idx;
#endif
            }
          else if(b == tree->n_root->b[2])
            {
              *n_v2         = tree->n_root->v[1];
              *Pij2         = tree->n_root->b[1]->Pij_rr;
              *tPij2        = tree->n_root->b[1]->tPij_rr;
              *p_lk_v2      = tree->n_root->b[1]->p_lk_rght;
              *sum_scale_v2 = tree->n_root->b[1]->sum_scale_rght;
#ifdef BEAGLE
              *child2_p_idx = tree->n_root->b[1]->p_lk_rght_idx;
              *Pij2_idx     = tree->n_root->b[1]->Pij_rr_idx;
#endif
            }
          else assert(FALSE);
        }
      else if(d == tree->n_root->v[1] || d == tree->n_root->v[2])
        {
          if(b == tree->n_root->b[1] || b == tree->n_root->b[2])
            {
              if(b == tree->n_root->b[1])
                {
                  *p_lk      = tree->n_root->b[1]->p_lk_rght;
                  *sum_scale = tree->n_root->b[1]->sum_scale_rght;
                  *p_lk_loc  = tree->n_root->b[1]->p_lk_loc_left;
#ifdef BEAGLE
                  *dest_p_idx = tree->n_root->b[1]->p_lk_rght_idx;
#endif
                }
              else
                {
                  *p_lk      = tree->n_root->b[2]->p_lk_rght;
                  *sum_scale = tree->n_root->b[2]->sum_scale_rght;
                  *p_lk_loc  = tree->n_root->b[2]->p_lk_loc_rght;
#ifdef BEAGLE
                  *dest_p_idx = tree->n_root->b[2]->p_lk_rght_idx;
#endif
                }

              *n_v1 = *n_v2 = NULL;
              for(i=0;i<3;++i)
                {
                  if(d->b[i] != tree->e_root)
                    {
                      if(!(*n_v1))
                        {
                          *n_v1 = d->v[i];
#ifdef BEAGLE
                          Set_Partial_Lk_One_Side(Pij1,tPij1,p_lk_v1,sum_scale_v1,d,d->b[i],tree,child1_p_idx,Pij1_idx);
#else
                          Set_Partial_Lk_One_Side(Pij1,tPij1,p_lk_v1,sum_scale_v1,d,d->b[i],tree);
#endif
                        }
                      else
                        {
                          *n_v2 = d->v[i];
#ifdef BEAGLE
                          Set_Partial_Lk_One_Side(Pij2,tPij2,p_lk_v2,sum_scale_v2,d,d->b[i],tree,child2_p_idx,Pij2_idx);
#else
                          Set_Partial_Lk_One_Side(Pij2,tPij2,p_lk_v2,sum_scale_v2,d,d->b[i],tree);
#endif
                        }
                    }
                }
            }
          else
            {
              if(d == b->left)
                {
                  *p_lk      = b->p_lk_left;
                  *sum_scale = b->sum_scale_left;
                  *p_lk_loc  = b->p_lk_loc_left;
#ifdef BEAGLE
                  *dest_p_idx = b->p_lk_left_idx;
#endif
                }
              else
                {
                  *p_lk      = b->p_lk_rght;
                  *sum_scale = b->sum_scale_rght;
                  *p_lk_loc  = b->p_lk_loc_rght;
#ifdef BEAGLE
                  *dest_p_idx = b->p_lk_rght_idx;
#endif
                }


              *n_v1 = tree->n_root;
#ifdef BEAGLE
              Set_Partial_Lk_One_Side(Pij1,tPij1,p_lk_v1,sum_scale_v1,d,
                                (d == tree->n_root->v[1])?
                                (tree->n_root->b[1]):
                                (tree->n_root->b[2]),
                                tree,child1_p_idx,Pij1_idx);
#else
              Set_Partial_Lk_One_Side(Pij1,tPij1,p_lk_v1,sum_scale_v1,d,
                                (d == tree->n_root->v[1])?
                                (tree->n_root->b[1]):
                                (tree->n_root->b[2]),
                                tree);
#endif
              for(i=0;i<3;i++)
                {
                  if(d->b[i] != tree->e_root && d->b[i] != b)
                    {
                      *n_v2 = d->v[i];
#ifdef BEAGLE
                      Set_Partial_Lk_One_Side(Pij2,tPij2,p_lk_v2,sum_scale_v2,d,d->b[i],tree,child2_p_idx,Pij2_idx);
#else
                      Set_Partial_Lk_One_Side(Pij2,tPij2,p_lk_v2,sum_scale_v2,d,d->b[i],tree);
#endif
                      break;
                    }
                }

            }
        }
      else
        {
          if(d == b->left)
            {
              *p_lk      = b->p_lk_left;
              *sum_scale = b->sum_scale_left;
              *p_lk_loc  = b->p_lk_loc_left;
#ifdef BEAGLE
              *dest_p_idx = b->p_lk_left_idx;
#endif
            }
          else
            {
              *p_lk      = b->p_lk_rght;
              *sum_scale = b->sum_scale_rght;
              *p_lk_loc  = b->p_lk_loc_rght;
#ifdef BEAGLE
              *dest_p_idx = b->p_lk_rght_idx;
#endif
            }

          *n_v1 = *n_v2 = NULL;
          for(i=0;i<3;i++)
            {
              if(d->b[i] != b)
                {
                  if(!(*n_v1))
                    {
                      *n_v1 = d->v[i];
#ifdef BEAGLE
                      Set_Partial_Lk_One_Side(Pij1,tPij1,p_lk_v1,sum_scale_v1,d,d->b[i],tree,child1_p_idx,Pij1_idx);
#else
                      Set_Partial_Lk_One_Side(Pij1,tPij1,p_lk_v1,sum_scale_v1,d,d->b[i],tree);
#endif
                    }
                  else
                    {
                      *n_v2 = d->v[i];
#ifdef BEAGLE
                      Set_Partial_Lk_One_Side(Pij2,tPij2,p_lk_v2,sum_scale_v2,d,d->b[i],tree,child2_p_idx,Pij2_idx);
#else
                      Set_Partial_Lk_One_Side(Pij2,tPij2,p_lk_v2,sum_scale_v2,d,d->b[i],tree);
#endif
                    }
                }
            }
        }
    }
}

/*     |
       |
       |d
      /  \
     /    \
    /      \b
   /        \
  /          \x (either n_v1 or n_v2)

  Returns p_lk and sum_scale for subtree with x as root, Pij for edge b
*/

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Set_Partial_Lk_One_Side(phydbl **Pij, phydbl **tPij, phydbl **p_lk,  int **sum_scale, t_node *d, t_edge *b, t_tree *tree
#ifdef BEAGLE
                                     , int* child_p_idx, int* Pij_idx
#endif
                                     )
{

  if(Pij != NULL)
    {
      *Pij  = b->Pij_rr;
      *tPij = b->tPij_rr;
#ifdef BEAGLE
      *Pij_idx = b->Pij_rr_idx;
#endif
    }

  if(d->tax == NO)
    {
      if(d == b->left) // if d is on the left of b, then d's neighbor is on the right
        {
          *p_lk      = (b->rght->tax == YES) ? b->p_lk_tip_r : b->p_lk_rght;
          *sum_scale = b->sum_scale_rght;
#ifdef BEAGLE
          *child_p_idx = b->rght->tax? b->p_lk_tip_idx: b->p_lk_rght_idx;
#endif
        }
      else
        {
          *p_lk      = b->p_lk_left;
          *sum_scale = b->sum_scale_left;
#ifdef BEAGLE
          *child_p_idx   = b->rght->tax? b->p_lk_tip_idx: b->p_lk_left_idx;
#endif
        }
    }
  else
    {
#ifdef BEAGLE
      Warn_And_Exit(TODO_BEAGLE);
#endif
      *p_lk        = NULL;
      *sum_scale   = NULL;
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
