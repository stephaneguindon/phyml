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
        PhyML_Fprintf(stderr,"\n. Unknown character state : '%c'.\n",state);
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
        PhyML_Fprintf(stderr,"\n. Unknown character state : '%c'.\n",state);
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
      PhyML_Fprintf(stderr,"\n. Unknown character state : '%c'.\n",aa);
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
      PhyML_Fprintf(stderr,"\n. Unknown character state : '%c'.\n",aa);
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

      if(tree->n_root != NULL)
        {
          for(i=0;i<3;++i)
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
          PhyML_Printf("\n. a->num: %d d->num: %d d->v[0]->num: %d d->v[1]->num: %d d->v[2]->num: %d d->b[0]->num: %d d->b[1]->num: %d d->b[2]->num: %d root ? %d e_root ? %d\n",
                       a?a->num:-1,
                       d?d->num:1,
                       d->v[0]?d->v[0]->num:-1,
                       d->v[1]?d->v[1]->num:-1,
                       d->v[2]?d->v[2]->num:-1,
                       d->b[0]?d->b[0]->num:-1,
                       d->b[1]?d->b[1]->num:-1,
                       d->b[2]?d->b[2]->num:-1,
                       tree->n_root?tree->n_root->num:-1,
                       tree->e_root?tree->e_root->num:-1);
          assert(FALSE);
        }

      /* PhyML_Printf("\n. a:%d [%d] d:%d dir:%d [%p %p] [%p %p %p] [%p %p]", */
      /*              a->num, */
      /*              a == tree->n_root, */
      /*              d->num, */
      /*              dir, */
      /*              d->b[dir], */
      /*              tree->e_root, */
      /*              d->b[0],d->b[1],d->b[2], */
      /*              tree->n_root->b[1],tree->n_root->b[2]); */
      
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
      Post_Order_Lk(tree->a_nodes[tree->tip_root],tree->a_nodes[tree->tip_root]->v[0],tree);
      if(tree->both_sides == YES)
        Pre_Order_Lk(tree->a_nodes[tree->tip_root],tree->a_nodes[tree->tip_root]->v[0],tree);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl Lk(t_edge *b, t_tree *tree)
{
  unsigned int br,catg,state,ambiguity_check,site;
  phydbl len,*expl,*dot_prod,*p_lk_left,*p_lk_rght;  

  const unsigned int ns = tree->mod->ns;
  const unsigned int ncatg = tree->mod->ras->n_catg;
  const unsigned int npatterns = tree->n_pattern;
  const unsigned int nsncatg = ns * ncatg;

  
  tree->numerical_warning = NO;
  
  /* if(tree->eval_alnL == NO) return UNLIKELY; */
  
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
      MIXT_Lk(b,tree);
      return tree->c_lnL;
    }
  
  tree->old_lnL = tree->c_lnL;
  

  if(tree->rates && tree->io && tree->io->lk_approx == NORMAL)
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
          for(br=0;br<2*tree->n_otu-3;++br)
            {
              Update_PMat_At_Given_Edge(tree->a_edges[br],tree);
            }
          
          if(tree->n_root && tree->ignore_root == NO)
            {
              Update_PMat_At_Given_Edge(tree->n_root->b[1],tree);
              Update_PMat_At_Given_Edge(tree->n_root->b[2],tree);
            }
        }
      else //Update PMat for a specific edge
        {
          if(tree->use_eigen_lr == NO)
            {
              if(tree->n_root &&
                 (b == tree->n_root->b[1] || b == tree->n_root->b[2]) &&
                 tree->ignore_root == YES)
                {
                  Update_PMat_At_Given_Edge(tree->e_root,tree);
                }
              else
                {
                  Update_PMat_At_Given_Edge(b,tree);
                }
            }
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
              Post_Order_Lk(tree->a_nodes[tree->tip_root],tree->a_nodes[tree->tip_root]->v[0],tree);
              if(tree->both_sides == YES)
                Pre_Order_Lk(tree->a_nodes[tree->tip_root],tree->a_nodes[tree->tip_root]->v[0],tree);
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
        b = tree->a_nodes[tree->tip_root]->b[0];
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
          for(state=0;state<ns;++state) expl[catg*ns+state] = (phydbl)pow(tree->mod->eigen->e_val[state],len);
        }
    }

  p_lk_left = b->p_lk_left;
  p_lk_rght = b->rght->tax ? b->p_lk_tip_r : b->p_lk_rght;

  for(site=0;site<npatterns;++site)
    {
      ambiguity_check = -1;
      state           = -1;
      tree->curr_site = site;
      
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
          if(tree->data->wght[site] > SMALL) Lk_Core_Eigen_Lr(expl,dot_prod,b,tree);
          dot_prod += nsncatg;
        }
      else
        {
          if(tree->data->wght[site] > SMALL) Lk_Core(state,ambiguity_check,p_lk_left,p_lk_rght,b->Pij_rr,b->tPij_rr,b,tree);

          if(b->rght->tax == YES)
            {
              p_lk_left += nsncatg;
              p_lk_rght += ns;
            }
          else
            {
              p_lk_left += nsncatg;
              p_lk_rght += nsncatg;
            }
        }
    }
#endif

  return tree->c_lnL;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// First derivative of the log-likelihood with respect
// to the length of edge b
phydbl dLk(phydbl *l, t_edge *b, t_tree *tree)
{
  unsigned int catg,state,site;
  phydbl len,rr;
  phydbl lk,dlk,dlnlk,lnlk;
  phydbl ev,expevlen;
   
  const unsigned int ns = tree->mod->ns;
  const unsigned int ncatg = tree->mod->ras->n_catg;
  const unsigned int npattern = tree->n_pattern;

  phydbl *dot_prod = tree->dot_prod;
  phydbl *expl = tree->expl;

  tree->numerical_warning = NO;
  
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
          assert(FALSE);
        }

      if(len < tree->mod->l_min)      len = tree->mod->l_min;
      else if(len > tree->mod->l_max) len = tree->mod->l_max;      
      // value of rr should be corrected too if any of these two conditions
      // is true. Leads to numerical precision issues though...
      
      
      for(state=0;state<ns;++state) 
        {
          ev = tree->mod->eigen->e_val[state];
          expevlen = exp(ev*len);

          expl[catg*2*ns + 2*state] = expevlen;
          expl[catg*2*ns + 2*state + 1] = expevlen*ev*rr;
        }
    }
    
  dlnlk  = 0.0;
  lnlk   = 0.0;
  
  for(site=0;site<npattern;++site)
    {
      if(tree->data->wght[site] > SMALL) 
        {
          tree->curr_site = site;
          
          Lk_dLk_Core_Eigen_Lr(expl,dot_prod+site*ns*ncatg,b,&lk,&dlk,tree);

          assert(lk > .0);
          
          dlk /= lk;
          dlnlk += tree->data->wght[site] * dlk;
          lnlk += tree->data->wght[site] * (log(lk) - (phydbl)LOG2 * tree->fact_sum_scale[site]);
        }
    }

  tree->c_dlnL = dlnlk;
  tree->c_lnL  = lnlk;

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

phydbl Lk_Core(int state, int ambiguity_check,
               phydbl *p_lk_left, phydbl *p_lk_rght,
               phydbl *Pij_rr,
               phydbl *tPij_rr,
               t_edge *b,
               t_tree *tree)
{
  phydbl site_lk,res,*pi,*site_lk_cat,log_site_lk;
  unsigned int catg;
  
  const unsigned int ns    = tree->mod->ns;
  const unsigned int ncatg = tree->mod->ras->n_catg;
  const unsigned int site  = tree->curr_site;  
  const unsigned nsns      = ns*ns;

  
  assert(tree->data->wght[site] > SMALL);
  
  pi = tree->mod->e_frq->pi->v;
  
  if(tree->mod->s_opt->skip_tree_traversal == NO)
    {
      for(catg=0;catg<ncatg;++catg)
        {
          if(ns == 4 || ns == 20)
            {
#if (defined(__AVX__) || defined(__AVX2__))
              tree->site_lk_cat[catg] = AVX_Lk_Core_One_Class_No_Eigen_Lr(p_lk_left,p_lk_rght,Pij_rr,tPij_rr,pi,ns,ambiguity_check,state);
#elif (defined(__SSE__) || defined(__SSE2__) || defined(__SSE3__))
              tree->site_lk_cat[catg] = SSE_Lk_Core_One_Class_No_Eigen_Lr(p_lk_left,p_lk_rght,Pij_rr,tPij_rr,pi,ns,ambiguity_check,state);
#else
              tree->site_lk_cat[catg] = Lk_Core_One_Class_No_Eigen_Lr(p_lk_left,p_lk_rght,Pij_rr,pi,ns,ambiguity_check,state);
#endif
            }
          else
            {
              tree->site_lk_cat[catg] = Lk_Core_One_Class_No_Eigen_Lr(p_lk_left,p_lk_rght,Pij_rr,pi,ns, ambiguity_check, state);
            }
          
          Pij_rr += nsns;
          tPij_rr += nsns;
          if(b->left->tax == NO) p_lk_left += ns;
          if(b->rght->tax == NO) p_lk_rght += ns;
        }

      Pull_Scaling_Factors(site,b,tree);

    }

  site_lk = .0;
  site_lk_cat = tree->unscaled_site_lk_cat + site*ncatg;
  for(catg=0;catg<ncatg;++catg) site_lk += site_lk_cat[catg] * tree->mod->ras->gamma_r_proba->v[catg];

  if(tree->mod->ras->invar == YES)
    {
      int num_prec_issue = NO;
      phydbl inv_site_lk = Invariant_Lk(tree->fact_sum_scale[site],site,&num_prec_issue,tree);

      switch(num_prec_issue)
        {
        case YES :         
          {
            assert(isinf(inv_site_lk));
            tree->fact_sum_scale[site] = 0;
            inv_site_lk = Invariant_Lk(0,site,&num_prec_issue,tree);
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

  if(tree->apply_lk_scaling == YES) res = site_lk / pow(2,tree->fact_sum_scale[site]);
  else                              res = site_lk;

  if(site_lk < SMALL)
    {
      site_lk = SMALL;
      tree->numerical_warning = YES;
    }

  
  log_site_lk = log(site_lk) - (phydbl)LOG2 * tree->fact_sum_scale[site]; // log_site_lk =  log(site_lk_scaled / 2^(left_subtree+right_subtree))
  tree->c_lnL_sorted[site] = log_site_lk;
  tree->c_lnL += tree->data->wght[site] * log_site_lk;
  tree->cur_site_lk[site] = exp(log_site_lk); // note to self : add opt out option to avoid calculating this if not necessary
    
  /* printf("\n. clnL: %f",tree->c_lnL); */
  return res;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl Lk_Core_Eigen_Lr(phydbl *expl, phydbl *dot_prod, t_edge *b, t_tree *tree)
{
  phydbl site_lk,res,*site_lk_cat,log_site_lk;
  unsigned int catg;
  int num_prec_issue;
  
  const unsigned int ns    = tree->mod->ns;
  const unsigned int ncatg = tree->mod->ras->n_catg;
  const unsigned int site  = tree->curr_site;

  assert(tree->data->wght[site] > SMALL);
  
  if(tree->mod->s_opt->skip_tree_traversal == NO)
    {
      for(catg=0;catg<ncatg;++catg)
        {
          if(tree->mod->io->datatype == NT || tree->mod->io->datatype == AA)
            {
#if (defined(__AVX__) || defined(__AVX2__))
              tree->site_lk_cat[catg] = AVX_Lk_Core_One_Class_Eigen_Lr(dot_prod,expl ? expl : NULL,ns);
#elif (defined(__SSE__) || defined(__SSE2__) || defined(__SSE3__))
              tree->site_lk_cat[catg] = SSE_Lk_Core_One_Class_Eigen_Lr(dot_prod,expl ? expl : NULL,ns);
#else
              tree->site_lk_cat[catg] = Lk_Core_One_Class_Eigen_Lr(dot_prod,expl ? expl : NULL,ns);
#endif
            }
          else
            {
              tree->site_lk_cat[catg] = Lk_Core_One_Class_Eigen_Lr(dot_prod,expl ? expl : NULL,ns);
            }

          dot_prod += ns;
          if(expl) expl += ns;
        }

      Pull_Scaling_Factors(site,b,tree);

    }
  
  
  site_lk_cat = tree->unscaled_site_lk_cat + site*ncatg;
  site_lk = .0;
  for(catg=0;catg<ncatg;++catg) site_lk += site_lk_cat[catg] * tree->mod->ras->gamma_r_proba->v[catg];
  
  if(tree->mod->ras->invar == YES)
    {
      num_prec_issue = NO;
      phydbl inv_site_lk = Invariant_Lk(tree->fact_sum_scale[site],site,&num_prec_issue,tree);

      switch(num_prec_issue)
        {
        case YES :
          {
            assert(isinf(inv_site_lk));
            tree->fact_sum_scale[site] = 0;
            inv_site_lk = Invariant_Lk(0,site,&num_prec_issue,tree);
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

  if(site_lk < SMALL)
    {
      site_lk = SMALL;
      tree->numerical_warning = YES;
    }

  // likelihood (or 1st, 2nd derivative) not rescaled here. Valid only if all partial likelihoods
  // were scaled using the same factor, i.e., when scaling_method == SCALE_FAST. In this case, the
  // scaling factors will cancel out in dlk/lk and d2lk/lk
  res = site_lk;
  
  log_site_lk = log(site_lk) - (phydbl)LOG2 * tree->fact_sum_scale[site]; // log_site_lk =  log(site_lk_scaled / 2^(left_subtree+right_subtree))
  tree->c_lnL_sorted[site] = log_site_lk;
  tree->c_lnL += tree->data->wght[site] * log_site_lk;
  tree->cur_site_lk[site] = exp(log_site_lk); // note to self : add opt out option to avoid calculating this if not necessary
      
  return res;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// Compute likelihood and first derivative of likelihood with respect to the length of edge b *unscaled* 
void Lk_dLk_Core_Eigen_Lr(phydbl *expl, phydbl *dot_prod, t_edge *b, phydbl *lk, phydbl *dlk, t_tree *tree)
{
  phydbl core_lk,core_dlk;
  unsigned int catg;
  int num_prec_issue;
  
  const unsigned int ns    = tree->mod->ns;
  const unsigned int ncatg = tree->mod->ras->n_catg;
  const unsigned int site  = tree->curr_site;
  
  *lk = *dlk = 0.0;

  assert(tree->data->wght[site] > SMALL);
  
  if(tree->mod->s_opt->skip_tree_traversal == NO)
    {
      for(catg=0;catg<ncatg;++catg)
        {
          if(tree->mod->io->datatype == NT || tree->mod->io->datatype == AA)
            {
#if (defined(__AVX__) || defined(__AVX2__))
              AVX_Lk_dLk_Core_One_Class_Eigen_Lr(dot_prod,
                                                 expl ? expl : NULL,
                                                 ns,&core_lk,&core_dlk);
#elif (defined(__SSE__) || defined(__SSE2__) ||  defined(__SSE3__))
              SSE_Lk_dLk_Core_One_Class_Eigen_Lr(dot_prod,
                                                 expl ? expl : NULL,
                                                 ns,&core_lk,&core_dlk);
#else
              Lk_dLk_Core_One_Class_Eigen_Lr(dot_prod,
                                             expl ? expl : NULL,
                                             ns,&core_lk,&core_dlk);
#endif
            }
          else
            {
              Lk_dLk_Core_One_Class_Eigen_Lr(dot_prod,
                                             expl ? expl : NULL,
                                             ns,&core_lk,&core_dlk);
            }

          *lk  += core_lk * tree->mod->ras->gamma_r_proba->v[catg];
          *dlk += core_dlk * tree->mod->ras->gamma_r_proba->v[catg];
          
          dot_prod += ns;
          if(expl) expl += 2*ns;
        }
      Pull_Scaling_Factors(site,b,tree);
    }
    
  if(tree->mod->ras->invar == YES)
    {
      num_prec_issue = NO;
      phydbl inv_site_lk = Invariant_Lk(tree->fact_sum_scale[site],site,&num_prec_issue,tree);

      switch(num_prec_issue)
        {
        case YES :
          {
            *lk = inv_site_lk * tree->mod->ras->pinvar->v;
            *dlk = 0.0;
            break;
          }
        case NO :
          {
            *lk = *lk * (1. - tree->mod->ras->pinvar->v) + inv_site_lk * tree->mod->ras->pinvar->v;
            *dlk = *dlk * (1. - tree->mod->ras->pinvar->v);
            break;
          }
        }
    }

  if(*lk < SMALL)
    {
      *lk = SMALL;
      tree->numerical_warning = YES;
    }
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Update_Eigen_Lr(t_edge *b, t_tree *tree)
{
  unsigned int site,catg,i,j;
  phydbl *dot_prod,*r_e_vect,*l_e_vect,*p_lk_left,*p_lk_rght,*pi;
  phydbl left,rght;
  
  const unsigned int npattern = tree->n_pattern;
  const unsigned int ns = tree->mod->ns;
  const unsigned int ncatg = tree->mod->ras->n_catg;
  const unsigned int nsncatg = ns*ncatg;

  
  if(tree->is_mixt_tree == YES)
    {
      MIXT_Update_Eigen_Lr(b,tree);
      return;
    }

  if(tree->mod->ns == 4 || tree->mod->ns == 20)
    {
#if (defined(__AVX__) || defined(__AVX2__))
      AVX_Update_Eigen_Lr(b,tree);
      return;
#elif (defined(__SSE__) || defined(__SSE2__) ||  defined(__SSE3__))
      SSE_Update_Eigen_Lr(b,tree);
      return;
#endif
    }
    
  assert(tree->update_eigen_lr == YES);
  
  dot_prod = tree->dot_prod;
  r_e_vect = tree->mod->eigen->r_e_vect;
  l_e_vect = tree->mod->eigen->l_e_vect;
  pi       = tree->mod->e_frq->pi->v;
    
  if(b->left->tax == YES) p_lk_left = b->p_lk_tip_l;
  else                    p_lk_left = b->p_lk_left;

  if(b->rght->tax == YES) p_lk_rght = b->p_lk_tip_r;
  else                    p_lk_rght = b->p_lk_rght;

  for(site=0;site<npattern;++site)
    {
      if(tree->data->wght[site] > SMALL)
        {
          for(catg=0;catg<ncatg;++catg)
            {
              for(i=0;i<ns;++i)
                {
                  left = rght = 0.0;
                  for(j=0;j<ns;++j)
                    {
                      left += r_e_vect[j*ns + i] * p_lk_left[j] * pi[j];
                      rght += l_e_vect[i*ns + j] * p_lk_rght[j];
                    }
                  dot_prod[i] = left*rght;
                }
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
          else p_lk_left += nsncatg;
          
          if(b->rght->tax == YES) p_lk_rght += ns;
          else p_lk_rght += nsncatg;

          dot_prod += nsncatg;
        }
    }
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

void Lk_dLk_Core_One_Class_Eigen_Lr(phydbl *dot_prod, phydbl *expl, unsigned int ns, phydbl *lk, phydbl *dlk)
{
  unsigned int i;

  *lk = *dlk = 0.0;
  for(i=0;i<ns;++i)
    {
      *lk  += dot_prod[i] * expl[2*i];
      *dlk += dot_prod[i] * expl[2*i+1];      
    }
}
 
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl Lk_Core_One_Class_No_Eigen_Lr(phydbl *p_lk_left, phydbl *p_lk_rght, phydbl *Pij, phydbl *pi, int ns, int ambiguity_check, int state)
{
  unsigned int l,k;
  phydbl lk = 0.0;
  phydbl sum;

    
  if(ambiguity_check == NO)/* tip case */
    {      
      sum = .0;
      Pij += state*ns;
      for(l=0;l<ns;++l) sum += Pij[l] * p_lk_left[l];               
      lk += sum * pi[state];
    }
  else /* If the character observed at the tip is ambiguous: ns x ns terms to consider */
    {
      for(k=0;k<ns;++k)
        {
          if(p_lk_rght[k] > .0) /* Only bother ascending into the subtrees if the likelihood of state k, at site "site*dim2" is > 0 */
            {
              sum = .0;
              for(l=0;l<ns;l++)
                {
                  sum += Pij[l] * p_lk_left[l];
                }
              lk += sum * pi[k] * p_lk_rght[k];
            }
          Pij += ns;
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
          
          /* printf("\n. inv_site_lk = %f [%c] [%d] invar: %d",inv_site_lk,tree->data->c_seq[0]->state[site],tree->data->invar[site],tree->data->invar[site]); */

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
              PhyML_Fprintf(stderr,"\n. fact_sum_scale: %d",fact_sum_scale);              
              PhyML_Fprintf(stderr,"\n. pi: %f",tree->mod->e_frq->pi->v[tree->data->invar[site]]);              
              for(i=0;i<tree->mod->ns;i++) PhyML_Fprintf(stderr,"\n. pi %d: %f",i,tree->mod->e_frq->pi->v[i]);
              PhyML_Fprintf(stderr,"\n. Numerical precision issue alert.");
              PhyML_Fprintf(stderr,"\n. File %s at line %d (function '%s')\n",__FILE__,__LINE__,__FUNCTION__);
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
  /* if(tree->eval_alnL == NO) return; */
  if(b->left == d && b->update_partial_lk_left == NO) return;
  if(b->rght == d && b->update_partial_lk_rght == NO) return;
  
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
      if(tree->mod->ns == 4 || tree->mod->ns == 20)
        {
#if (defined(__AVX__) || defined(__AVX2__))
          AVX_Update_Partial_Lk(tree,b,d);
#elif (defined(__SSE__) || defined(__SSE2__) || defined(__SSE3__))
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
  smallest_p_lk               = BIG;

  Set_All_Partial_Lk(&n_v1,&n_v2,
                     &p_lk,&sum_scale,&p_lk_loc,
                     &Pij1,&tPij1,&p_lk_v1,&sum_scale_v1,
                     &Pij2,&tPij2,&p_lk_v2,&sum_scale_v2,
                     d,b,tree);

  /* For every site in the alignment */
  for(site=0;site<n_patterns;site++)
    {
      if(tree->data->wght[site] > SMALL)
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
              
              smallest_p_lk = BIG;
              
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
                  
                  /* if(site == 0) PhyML_Printf("\n+ site: %d %G",site,p_lk[site*ncatgns+catg*ns+i]); */
                  
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
              
              sum_scale[site*ncatg+catg] = sum_scale_v1_val + sum_scale_v2_val;
              
              assert(sum_scale[site] < 1024);

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
      else
        {
          for(i=0;i<ns*ncatg;++i) p_lk[site*ncatgns + i] = 0.0;
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
  phydbl *p_lk;
  phydbl *p_lk_v1,*p_lk_v2;//Partial likelihood vector of node d, d's "left" neighbor, d's "right" neighbor. We fill *p_lk, and assume *p_lk_v1 and *p_lk_v2 are already filled.
  phydbl *Pij1,*Pij2;
  phydbl *tPij1,*tPij2;
  int *sum_scale, *sum_scale_v1, *sum_scale_v2;
  int *p_lk_loc;//Suppose site j, of a certain subtree, has "A" on one tip, and "C" on the other. If you come across this pattern again at site i<j, then you can simply copy the partial likelihoods
  
  const unsigned int ncatg = tree->mod->ras->n_catg;
  const unsigned int ns = tree->mod->ns;
  const unsigned int n_patterns = tree->n_pattern;


  if(tree->n_root && tree->ignore_root == YES &&
     (d == tree->n_root->v[1] || d == tree->n_root->v[2]) &&
     (b == tree->n_root->b[1] || b == tree->n_root->b[2]))
    {
      assert(FALSE);
    }
  
  if(d->tax)
    {
      PhyML_Fprintf(stderr,"\n. t_node %d is a leaf...",d->num);
      PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d (function '%s')\n",__FILE__,__LINE__,__FUNCTION__);
      Exit("\n");
    }


  n_v1 = n_v2                 = NULL;
  p_lk = p_lk_v1 = p_lk_v2    = NULL;
  Pij1 = Pij2                 = NULL;
  tPij1 = tPij2               = NULL;
  p_lk_loc                    = NULL;
  sum_scale_v1                = NULL;
  sum_scale_v2                = NULL;
  
  Set_All_Partial_Lk(&n_v1,&n_v2,
                     &p_lk,&sum_scale,&p_lk_loc,
                     &Pij1,&tPij1,&p_lk_v1,&sum_scale_v1,
                     &Pij2,&tPij2,&p_lk_v2,&sum_scale_v2,
                     d,b,tree);
  
  Core_Default_Update_Partial_Lk(n_v1,n_v2,
                                 p_lk,p_lk_v1,p_lk_v2,
                                 Pij1,Pij2,
                                 sum_scale,sum_scale_v1,sum_scale_v2,
                                 ns,ncatg,n_patterns,
                                 tree->apply_lk_scaling,
                                 tree->data->wght);
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Core_Default_Update_Partial_Lk(const t_node *n_v1, const t_node *n_v2,
                                    phydbl *plk0, const phydbl *plk1, const phydbl *plk2,
                                    const phydbl *Pij1, const phydbl *Pij2,
                                    int *sum_scale0, const int *sum_scale1, const int *sum_scale2,
                                    const int ns, const int ncatg, const int npatterns, const int apply_scaling,
                                    const phydbl *wght)
{
  
  unsigned int i,site,ncatgns,catg,nsns;
  int state_v1,state_v2;
  int ambiguity_check_v1,ambiguity_check_v2;
  int sum_scale_v1_val, sum_scale_v2_val;
  phydbl largest_p_lk;
  const phydbl *init_Pij1, *init_Pij2;
  
  ncatgns = ncatg*ns;
  nsns = ns*ns;
  init_Pij1 = Pij1;
  init_Pij2 = Pij2;
  
  /* For every site in the alignment */
  for(site=0;site<npatterns;site++)
    {
      if(wght[site] > SMALL)
        {
          state_v1 = state_v2 = -1;
          ambiguity_check_v1 = ambiguity_check_v2 = YES;
          
          /* n_v1 and n_v2 are tip nodes */
          if(n_v1->tax)
            {
              /* Is the state at this tip ambiguous? */
              ambiguity_check_v1 = n_v1->c_seq->is_ambigu[site];
              if(ambiguity_check_v1 == NO) state_v1 = n_v1->c_seq->d_state[site];
            }
          
          if(n_v2->tax)
            {
              /* Is the state at this tip ambiguous? */
              ambiguity_check_v2 = n_v2->c_seq->is_ambigu[site];
              if(ambiguity_check_v2 == NO) state_v2 = n_v2->c_seq->d_state[site];
            }
          
          Pij1 = init_Pij1;
          Pij2 = init_Pij2;
          
          /* For all the rate classes */
          for(catg=0;catg<ncatg;++catg)
            {
              if(ambiguity_check_v1 == NO && ambiguity_check_v2 == NO)
                {
                  Partial_Lk_Exex(Pij1,state_v1,
                                  Pij2,state_v2,
                                  ns,plk0);
                }
              else if(ambiguity_check_v1 == YES && ambiguity_check_v2 == NO)
                {
                  Partial_Lk_Exin(Pij2,state_v2,
                                  Pij1,plk1,
                                  ns,plk0);
                }
              else if(ambiguity_check_v1 == NO && ambiguity_check_v2 == YES)
                {
                  Partial_Lk_Exin(Pij1,state_v1,
                                  Pij2,plk2,
                                  ns,plk0);
                }
              else
                {
                  Partial_Lk_Inin(Pij1,plk1,
                                  Pij2,plk2,
                                  ns,plk0);
                }
              
              Pij1 += nsns;
              Pij2 += nsns;
              
              plk1 += (n_v1->tax) ? 0 : ns;
              plk2 += (n_v2->tax) ? 0 : ns;
              plk0 += ns;
            }
          
          plk1 += (n_v1->tax) ? ns : 0;
          plk2 += (n_v2->tax) ? ns : 0;
          
          sum_scale_v1_val = (sum_scale1)?(sum_scale1[site]):(0);
          sum_scale_v2_val = (sum_scale2)?(sum_scale2[site]):(0);
          sum_scale0[site] = sum_scale_v1_val + sum_scale_v2_val;
          
          plk0 -= ncatgns;
          largest_p_lk = -BIG;
          for(i=0;i<ncatgns;++i)
            if(plk0[i] > largest_p_lk)
              largest_p_lk = plk0[i];
          
          if(largest_p_lk < INV_TWO_TO_THE_LARGE && apply_scaling == YES)
            {
              for(i=0;i<ncatgns;++i) plk0[i] *= TWO_TO_THE_LARGE;
              sum_scale0[site] += LARGE;
            }
          plk0 += ncatgns;
        }
      else
        {
          plk0 += ncatgns;
          plk1 += (n_v1->tax) ? ns : ncatgns;
          plk2 += (n_v2->tax) ? ns : ncatgns;          
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
  tmpdata->obs_state_frq  = (phydbl *)mCalloc(mod->ns,sizeof(phydbl));
  tmpdata->ambigu = (short int *)mCalloc(data->crunch_len,sizeof(short int));
  F               = (phydbl *)mCalloc(mod->ns*mod->ns,sizeof(phydbl ));
  eigen_struct    = (eigen *)Make_Eigen_Struct(mod->ns);

  Set_Update_Eigen(YES,mod);
  Update_Boundaries(mod);
  Update_Eigen(mod);

  tmpdata->n_otu = 2;

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
          
          for(i=0;i<mod->ns*mod->ns;++i) F[i]=.0;
          len = 0.0;
          for(l=0;l<twodata->c_seq[0]->len;++l)
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
              for(i=0;i<mod->ns*mod->ns;++i) F[i] /= len;
            }
          
          sum = 0.;
          for(i=0;i<mod->ns*mod->ns;++i) sum += F[i];
          
          /* for(i=0;i<mod->ns*mod->ns;++i) PhyML_Printf("\n. %g",F[i]); */

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
  Free(tmpdata->obs_state_frq);
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
// Multinomial log likelihood
void Unconstraint_Lk(t_tree *tree)
{
  int i;

  tree->unconstraint_lk = .0;
  for(i=0;i<tree->data->crunch_len;i++) tree->unconstraint_lk += tree->data->wght[i]*(phydbl)log(tree->data->wght[i]);
  tree->unconstraint_lk -= tree->data->init_len*(phydbl)log(tree->data->init_len);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// Log-likelihood assuming a tree with edge of infinite lengths
void Composite_Lk(t_tree *tree)
{
  int i;
  tree->composite_lk = 0.0;
  for(i=0;i<tree->mod->ns;++i)
    tree->composite_lk +=
      tree->data->obs_state_frq[i]*
      tree->data->init_len*
      tree->n_otu*
      log(tree->mod->e_frq->pi->v[i]);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Init_Partial_Lk_Tips_Double(t_tree *tree)
{
  unsigned int curr_site,i,dim1;

  if(tree->is_mixt_tree == YES) return;
  
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
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Init_Partial_Lk_Tips_Int(t_tree *tree)
{
  int curr_site,i,dim1;

  if(tree->is_mixt_tree == YES) return;

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
  assert(tree->eval_alnL == YES);
  
  if(tree->is_mixt_tree == YES)
    {
      MIXT_Update_PMat_At_Given_Edge(b_fcus,tree);
      return;
    }
  
  if(b_fcus->Pij_rr == NULL)
    {
      PhyML_Printf("\n. b_fcus is e_root ? %d node left: %d node rght: %d left is root ? %d right is root ? %d [%p] [%d] [%d]",
                   (b_fcus == tree->e_root) ? 1 : 0,
                   b_fcus->left->num,
                   b_fcus->rght->num,
                   (b_fcus->left == tree->n_root) ? 1 : 0,
                   (b_fcus->rght == tree->n_root) ? 1 : 0,
                   tree->aux_tree,
                   tree->eval_alnL,
                   tree->is_mixt_tree);
      assert(false);
    }

  if(tree->mixt_tree != NULL) assert(tree->mod->ras->n_catg == 1);
  
  if(tree->io->mod->gamma_mgf_bl == YES) Set_Br_Len_Var(b_fcus,tree);

  l_min = tree->mod->l_min;
  l_max = tree->mod->l_max;

  len = -1.0;

  if(tree->mod->log_l == YES) b_fcus->l->v = exp(b_fcus->l->v);

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
                assert(FALSE);
              }
            break;
          }
#ifdef DEBUG
      if(j == 3)
        {
          PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d.\n",__FILE__,__LINE__);
          assert(FALSE);
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
  /* if(mod->is_mixt_mod == YES) mod = mod->next; */
  /* assert(mod); */

  if(mod->log_l == YES) dist = exp(dist);

  for(k=0;k<mod->ras->n_catg;k++)
    {
      len = dist*mod->ras->gamma_rr->v[k];
      if(len < mod->l_min)      len = mod->l_min;
      else if(len > mod->l_max) len = mod->l_max;
      PMat(len,mod,mod->ns*mod->ns*k,mod->Pij_rr->v,NULL);
      /* PhyML_Printf("\n. p: %g len: %g",mod->Pij_rr->v[0],len); */
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
              /* PhyML_Printf("\n. pijk: %g lnL: %g",pijk,lnL); */
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
          /* PhyML_Printf("\n. pijk: %g lnL: %g",pijk,lnL); */
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

void Init_Partial_Lk_Loc(t_tree *tree)
{
  int i,j;
  t_node *d;
  int *patt_id_d;

  if(tree->is_mixt_tree == YES) return;
  
  for(i=0;i<2*tree->n_otu-1;++i)
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
          assert(tree->a_nodes[d->num]->c_seq);
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

#if (defined GEO)
phydbl Wrap_Geo_Lk(t_edge *b, t_tree *tree, supert_tree *stree)
{
  TIPO_Lk(tree);
  return tree->geo_lnL;
}
#endif

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Wrap_Lk_At_Given_Edge(t_edge *b, t_tree *tree, supert_tree *stree)
{
  Lk(b,tree);
  return tree->c_lnL;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

#if(defined PHYREX || PHYTIME)
phydbl Wrap_Lk_Rates(t_edge *b, t_tree *tree, supert_tree *stree)
{
  RATES_Lk(tree);
  return tree->rates->c_lnL;
}
#endif

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

#if(defined PHYREX || PHYTIME)
phydbl Wrap_Lk_Times(t_edge *b, t_tree *tree, supert_tree *stree)
{
  TIMES_Lk(tree);
  return tree->times->c_lnL;
}
#endif

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

int Check_Lk_At_Given_Edge(int verbose, t_tree *tree)
{
  int res;
  int i;
  phydbl *lk;

  lk = (phydbl *)mCalloc(2*tree->n_otu-3,sizeof(phydbl));

  res = 0;
  for(i=0;i<2*tree->n_otu-3;++i)
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
                    PhyML_Fprintf(stderr,"\n. tree: %s\n",Write_Tree(tree));
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
      do Optimize_Br_Len_Serie (2,tree); while(n_opt++ < 3);

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
#ifdef BEAGLE
          *dest_p_idx = b->p_lk_left_idx;
#endif
        }
      else
        {
          *p_lk      = b->p_lk_rght;
          *sum_scale = b->sum_scale_rght;
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
              else
                {
                  PhyML_Printf("\n. Issue detected with node %d.\n",d->num);
                  assert(FALSE);
                }
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
#ifdef BEAGLE
              *dest_p_idx = tree->n_root->b[1]->p_lk_left_idx;
#endif
            }
          else
            {
              *p_lk      = tree->n_root->b[2]->p_lk_left;
              *sum_scale = tree->n_root->b[2]->sum_scale_left;
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
#ifdef BEAGLE
                  *dest_p_idx = tree->n_root->b[1]->p_lk_rght_idx;
#endif
                }
              else
                {
                  *p_lk      = tree->n_root->b[2]->p_lk_rght;
                  *sum_scale = tree->n_root->b[2]->sum_scale_rght;
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

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

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

          if(*p_lk == NULL) PhyML_Printf("\n. b:%d b->left:%d b->rght:%d d:%d",
                                         b->num,
                                         b->left->num,
                                         b->rght->num,
                                         d->num);
          assert(*p_lk);
        }
      else
        {
          *p_lk      = b->p_lk_left;
          *sum_scale = b->sum_scale_left;
#ifdef BEAGLE
          *child_p_idx   = b->rght->tax? b->p_lk_tip_idx: b->p_lk_left_idx;
#endif

          if(*p_lk == NULL) PhyML_Printf("\n. b:%d b->left:%d b->rght:%d d:%d",
                                         b->num,
                                         b->left->num,
                                         b->rght->num,
                                         d->num);
          assert(*p_lk);
        }
    }
  else
    {
#ifdef BEAGLE
      Warn_And_Exit(TODO_BEAGLE);
#endif
      *p_lk        = NULL;
      *sum_scale   = NULL;
      PhyML_Printf("\n. WARNING. p_lk set to NULL. d->num: %d b->num: %d",d->num,b->num);
        
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Switch_Partial_Lk_Post(t_node *a, t_node *d, t_edge *b, short int yesno, t_tree *tree)
{
  if(a == tree->n_root) assert(FALSE);
  if(d->tax == NO)
    {
      int i;

      for(i=0;i<3;++i)
        if(d->v[i] != a)
          Switch_Partial_Lk_Post(d,d->v[i],d->b[i],yesno,tree);
    } 

  if(b->left == d) b->update_partial_lk_left = yesno;
  else if(b->rght == d) b->update_partial_lk_rght = yesno;
  else assert(FALSE);

  return;
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Switch_Partial_Lk_Pre(t_node *a, t_node *d, t_edge *b, short int yesno, t_tree *tree)
{  
  if(a == tree->n_root) assert(FALSE);
  if(d->tax == YES) return;
  else
    {
      int i;
      
      for(i=0;i<3;++i)
        {
          if(d->v[i] != a)
            {
              if(d->b[i]->left == d) d->b[i]->update_partial_lk_left = yesno;
              else if(d->b[i]->rght == d) d->b[i]->update_partial_lk_rght = yesno;
              else assert(FALSE);
            }
        }
      
      for(i=0;i<3;++i)
        if(d->v[i] != a)
          Switch_Partial_Lk_Pre(d,d->v[i],d->b[i],yesno,tree);
    }
  
  return;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Partial_Lk_Inin(const phydbl *Pij1, const phydbl *plk1, const phydbl *Pij2, const phydbl *plk2, const int ns, phydbl *plk0)
{
  unsigned int i,j;

  
  for(i=0;i<ns;++i) if(plk1[i] > 1.0 || plk1[i] < 1.0 || plk2[i] > 1.0 || plk2[i] < 1.0) break; 

  if(i != ns)
    {
      for(i=0;i<ns;++i)
        {
          phydbl u1 = 0.0;
          phydbl u2 = 0.0;

          for(j=0;j<ns;++j)
            {
              u1 += Pij1[j] * plk1[j];
              u2 += Pij2[j] * plk2[j];
            }
          
          Pij1 += ns;
          Pij2 += ns;
          plk0[i] = u1*u2;
        }
    }
  else
    {
      for(i=0;i<ns;++i) plk0[i] = 1.0;
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Partial_Lk_Exex(const phydbl *Pij1, const int state1, const phydbl *Pij2, const int state2, const int ns, phydbl *plk0)
{
  unsigned int i;
  for(i=0;i<ns;++i)
    {
      plk0[i] = Pij1[state1]*Pij2[state2];
      Pij1 += ns;
      Pij2 += ns;      
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Partial_Lk_Exin(const phydbl *Pij1, const int state1, const phydbl *Pij2, const phydbl *plk2, const int ns, phydbl *plk0)
{
  unsigned int i,j;
  
  for(i=0;i<ns;++i)
    {
      phydbl u2 = 0.0;
      for(j=0;j<ns;++j) u2 += Pij2[j] * plk2[j];
      plk0[i] = Pij1[state1]*u2;
      Pij1 += ns;
      Pij2 += ns;
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
