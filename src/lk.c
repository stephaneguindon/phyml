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
    PhyML_Printf("\n== Unknown character state : %c\n",state);
    Exit("\n== Init failed (data type supposed to be DNA)\n");
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
        PhyML_Printf("\n. Unknown character state : %c\n",state);
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

  For(i,20) p_lk[pos+i] = .0;

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

  case 'X' : case '?' : case '-' : For(i,20) p_lk[pos+i] = 1.; break;
  default :
    {
      PhyML_Printf("\n. Unknown character state : %c\n",aa);
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

  For(i,20) p_pars[pos+i] = .0;

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

  case 'X' : case '?' : case '-' : For(i,20) p_pars[pos+i] = 1; break;
  default :
    {
      PhyML_Printf("\n. Unknown character state : %c\n",aa);
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

  For(i,ns) p_lk[pos+i] = 0.;

  if(Is_Ambigu(state,GENERIC,state_len)) For(i,ns) p_lk[pos+i] = 1.;
  else
    {
      char format[6];
      sprintf(format,"%%%dd",state_len);
      if(!sscanf(state,format,&state_int))
    {
      PhyML_Printf("\n== state='%c'",state);
      PhyML_Printf("\n== Err in file %s at line %d (function '%s')\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
      if(state_int > ns)
    {
      PhyML_Printf("\n. %s %d cstate: %.2s istate: %d state_len: %d.\n",__FILE__,__LINE__,state,state_int,state_len);
      PhyML_Printf("\n== Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
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

  For(i,ns) p_pars[pos+i] = 0;

  if(Is_Ambigu(state,GENERIC,state_len)) For(i,ns) p_pars[pos+i] = 1;
  else
    {
      char format[6];
      sprintf(format,"%%%dd",state_len);
      if(!sscanf(state,format,&state_int))
    {
      PhyML_Printf("\n. state='%c'",state);
      PhyML_Printf("\n== Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
      Warn_And_Exit("");
    }
      if(state_int > ns)
    {
      PhyML_Printf("\n. %s %d cstate: %.2s istate: %d state_len: %d.\n",__FILE__,__LINE__,state,state_int,state_len);
      PhyML_Printf("\n== Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
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
  Update_P_Lk(tree,b_fcus,d);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* void Post_Order_Lk(t_node *a, t_node *d, t_tree *tree) */
/* { */
/*   int i,dir; */

/*   if(tree->is_mixt_tree == YES) */
/*     { */
/*       MIXT_Post_Order_Lk(a,d,tree); */
/*       return; */
/*     } */

/*   dir = -1; */

/*   if(d->tax)  */
/*     { */
/*       Get_All_Partial_Lk_Scale(tree,d->b[0],a,d); */
/*       return; */
/*     } */
/*   else */
/*     { */
/*       For(i,3) */
/* 	{ */
/* 	  if(d->v[i] != a) */
/* 	    Post_Order_Lk(d,d->v[i],tree); */
/* 	  else dir = i; */
/* 	}       */
/*       Get_All_Partial_Lk_Scale(tree,d->b[dir],a,d); */
/*     } */
/* } */

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
          For(i,3)
            {
              if(d->v[i] != a && d->b[i] != tree->e_root)
                Post_Order_Lk(d,d->v[i],tree);
              else dir = i;
            }
        }
      else
        {
          For(i,3)
            {
              if(d->v[i] != a)
                Post_Order_Lk(d,d->v[i],tree);
              else dir = i;
            }
        }
      
      if(tree->ignore_root == NO && d->b[dir] == tree->e_root)
        {
          if(d == tree->n_root->v[1]) Get_All_Partial_Lk_Scale(tree,tree->n_root->b[1],tree->n_root,d);
          else                        Get_All_Partial_Lk_Scale(tree,tree->n_root->b[2],tree->n_root,d);
        }
      else
        {
          Get_All_Partial_Lk_Scale(tree,d->b[dir],a,d);
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
          For(i,3)
            {
              if(d->v[i] != a && d->b[i] != tree->e_root)
                {
                  Get_All_Partial_Lk_Scale(tree,d->b[i],d->v[i],d);
                  Pre_Order_Lk(d,d->v[i],tree);
                }
            }
        }
      else
        {
          For(i,3)
            {
              if(d->v[i] != a)
                {
                  Get_All_Partial_Lk_Scale(tree,d->b[i],d->v[i],d);
                  Pre_Order_Lk(d,d->v[i],tree);
                }
            }
        }
    }      
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl Lk(t_edge *b, t_tree *tree)
{
  int br;
  
  if(b == NULL && tree->mod->s_opt->curr_opt_free_rates == YES)
    {
      tree->mod->s_opt->curr_opt_free_rates = NO;
      Optimize_Free_Rate_Weights(tree,YES,YES);
      tree->mod->s_opt->curr_opt_free_rates = YES;
    }
  
  
  if(tree->is_mixt_tree) {
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
  
  if(!b) Set_Model_Parameters(tree->mod);
  
  Set_Br_Len_Var(tree);
  
  if(tree->mod->s_opt->skip_tree_traversal == NO)
    {
      if(!b)//Update the PMat for all edges
        {
          For(br,2*tree->n_otu-3)
            {
              Update_PMat_At_Given_Edge(tree->a_edges[br],tree);
            }
          if(tree->n_root && tree->ignore_root == NO)
            {
              Update_PMat_At_Given_Edge(tree->n_root->b[1],tree);
              Update_PMat_At_Given_Edge(tree->n_root->b[2],tree);
            }
        }
      else//Update the PMat for a specific edge
        {
          Update_PMat_At_Given_Edge(b,tree);
        }
      
      if(!b)
        {
          if(tree->n_root)
            {
              if(tree->ignore_root == NO)
                {
                  Post_Order_Lk(tree->n_root,tree->n_root->v[1],tree);
                  Post_Order_Lk(tree->n_root,tree->n_root->v[2],tree);
                  
                  Update_P_Lk(tree,tree->n_root->b[1],tree->n_root);
                  Update_P_Lk(tree,tree->n_root->b[2],tree->n_root);
                  
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

  int n_patterns,ambiguity_check,state;
  n_patterns = tree->n_pattern;
  For(tree->curr_site,n_patterns)
    {
      ambiguity_check = -1;
      state           = -1;

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

          Lk_Core(state,ambiguity_check,b,tree);

          /* fprintf(stdout,"site %d lnL: %g:%d %g\n",tree->curr_site,tree->c_lnL_sorted[tree->curr_site],tree->data->wght[tree->curr_site],tree->c_lnL); */
          /* fflush(stdout); */
        }
    }
#endif

/*   Qksort(tree->c_lnL_sorted,NULL,0,n_patterns-1); */

/*   tree->c_lnL = .0; */
/*   For(tree->curr_site,n_patterns) */
/*     { */
/*       tree->c_lnL += tree->c_lnL_sorted[tree->curr_site]; */
/*     } */


  Adjust_Min_Diff_Lk(tree);

//  Print_All_Edge_PMats(tree);
//  Print_All_Edge_Likelihoods(tree);
//  DUMP_D(tree->c_lnL);

  return tree->c_lnL;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* Core of the likelihood calculcation. Assume that the partial likelihoods on both
   sides of t_edge *b are up-to-date. Calculate the log-likelihood at one site.
*/

phydbl Lk_Core(int state, int ambiguity_check, t_edge *b, t_tree *tree)
{
  phydbl log_site_lk;
  phydbl site_lk_cat/*rate specific site lk*/, site_lk, inv_site_lk;
  int fact_sum_scale;
  phydbl sum;
  int catg/*Index of the current rate classes*/,ns/*size of the state-space*/,k,l,site/*index of the current site in the MSA*/;
  int dim1,dim2,dim3;
  int exponent;
  int num_prec_issue;

#ifdef BEAGLE
  dim1 = tree->n_pattern * tree->mod->ns;
#else
  dim1 = tree->mod->ras->n_catg * tree->mod->ns;//Dimension of a matrix L that holds rate-specific character likelihoods. IOW, L[ij] is the likelihood of character j in rate class i
#endif
  dim2 = tree->mod->ns;
  dim3 = tree->mod->ns * tree->mod->ns;//Dimensions of the transition prob. matrix
  
  log_site_lk     = .0;
  site_lk         = .0;
  site_lk_cat     = .0;
  site            = tree->curr_site;
  ns              = tree->mod->ns;
  
  
  /* Skip this if no tree traveral was required, i.e. likelihood in each class of the mixture is already up to date */
  if(tree->mod->s_opt->skip_tree_traversal == NO)
    {
      /* Actual likelihood calculation */
      /* For all classes of rates */
      For(catg,tree->mod->ras->n_catg)
        {
          site_lk_cat = .0;

          /* b is an external edge */
          if((b->rght->tax) && (!tree->mod->s_opt->greedy)) /* By convention, tips are always on the right of an external edge */
            {
              if(!ambiguity_check)/* If the character observed at the tip is NOT ambiguous: ns x 1 terms to consider */
                {
                  sum = .0;
                  For(l,ns)
                    {
                      sum +=
                        b->Pij_rr[(catg*dim3)+(state*dim2)+l] *
                        b->p_lk_left[(site*dim1)+(catg*dim2)+l];
                    }
                  
                  site_lk_cat += sum * tree->mod->e_frq->pi->v[state];
                }
              else /* If the character observed at the tip is ambiguous: ns x ns terms to consider */
                {
                  For(k,ns)
                    {
                      sum = .0;
                      if(b->p_lk_tip_r[site*dim2+k] > .0) /* Only bother ascending into the subtrees if the likelihood of state k, at site "site*dim2" is > 0 */
                        {
                          For(l,ns)
                            {
                              sum +=
                                b->Pij_rr[(catg*dim3)+(k*dim2)+l] *
                                b->p_lk_left[(site*dim1)+(catg*dim2)+l];
                            }
                          
                          site_lk_cat +=
                            sum *
                            tree->mod->e_frq->pi->v[k] *
                            b->p_lk_tip_r[site*dim2+k];
                        }
                    }
                }
            }
          else /* b is an internal edge: ns x ns terms to consider */
            {
              For(k,ns)
                {
                  sum = .0;
                  if(b->p_lk_rght[(site*dim1)+(catg*dim2)+k] > .0) //Only bother descending into the subtrees if the likelihood of state k, at site "site*dim2" is > 0
                    {
                      For(l,ns)
                        {
                          sum +=
                            b->Pij_rr[(catg*dim3)+(k*dim2)+l] *     //"catg*dim3" offsets into a rate-specific flat 4*4 DNA matrix P. Then, P[i,j] is given by "k*4+l"
                            b->p_lk_left[(site*dim1)+(catg*dim2)+l];//"site*dim1" offsets into a rate-specific flat num_rates*4 matrix L. Then, L[i,j] is given by "catg*dim2+l"
                        }
                      
                      site_lk_cat +=
                        sum *
                        tree->mod->e_frq->pi->v[k] *
                        b->p_lk_rght[(site*dim1)+(catg*dim2)+k];
                    }
                }
            }
          
          tree->site_lk_cat[catg] = site_lk_cat;
          
        } /* site likelihood for all rate classes */
      Pull_Scaling_Factors(site,b,tree);
    }

  fact_sum_scale = tree->fact_sum_scale[site];


  //Likelihood of the site; is the sum of the individual rate specific likelihoods
  site_lk = .0;
  For(catg,tree->mod->ras->n_catg)
    {
      site_lk +=
        tree->unscaled_site_lk_cat[catg*tree->n_pattern + site]* 
        tree->mod->ras->gamma_r_proba->v[catg]; //density
    }
  
  if(tree->mod->ras->invar == YES)
    {
      num_prec_issue = NO;
      inv_site_lk = Invariant_Lk(fact_sum_scale,site,&num_prec_issue,tree);
      
      if(num_prec_issue == YES) // inv_site_lk >> site_lk
        {
          site_lk = inv_site_lk * tree->mod->ras->pinvar->v;
        }
      else
        {
          site_lk = site_lk * (1. - tree->mod->ras->pinvar->v) + inv_site_lk * tree->mod->ras->pinvar->v;
        }          
    }
  
  log_site_lk = LOG(site_lk) - (phydbl)LOG2 * fact_sum_scale; // log_site_lk =  log(site_lk_scaled / 2^(left_subtree+right_subtree))      
  
  // Calculation of the site likelihood (using scaling factors)...
  int piecewise_exponent;
  phydbl multiplier;
  if(fact_sum_scale >= 0)
    {
      tree->cur_site_lk[site] = site_lk;
      exponent = -fact_sum_scale;
      do
        {
          piecewise_exponent = MAX(exponent,-63);
          multiplier = 1. / (phydbl)((unsigned long long)(1) << -piecewise_exponent);
          tree->cur_site_lk[site] *= multiplier;
          exponent -= piecewise_exponent;
        }
      while(exponent != 0);
    }
  else
    {
      //In some cases fact_sum_scale can be negative. If you rescale the partials of two independent subtrees and make some of
      //these numbers large in order to avoid underflow, then there is a chance that when multiplied them together you will
      //get an overflow, in which case fact_sum_scale can become negative.
      
      tree->cur_site_lk[site] = site_lk;
      exponent = fact_sum_scale;
      do
        {
          piecewise_exponent = MIN(exponent,63);
          multiplier = (phydbl)((unsigned long long)(1) << piecewise_exponent);
          tree->cur_site_lk[site] *= multiplier;
          exponent -= piecewise_exponent;
        }
      while(exponent != 0);
    }
  
  // ... or using the log-likelihood
  if(isinf(site_lk) || isnan(site_lk))
    {
      tree->cur_site_lk[site] = EXP(log_site_lk);
    }
  
  if(isinf(log_site_lk) || isnan(log_site_lk))
    {
      PhyML_Printf("\n== Site = %d",site);
      PhyML_Printf("\n== Invar = %d",tree->data->invar[site]);
      PhyML_Printf("\n== Mixt = %d",tree->is_mixt_tree);
      PhyML_Printf("\n== Lk = %G log(Lk) = %f < %G",site_lk,log_site_lk,-BIG);
      For(catg,tree->mod->ras->n_catg) PhyML_Printf("\n== rr=%f p=%f",tree->mod->ras->gamma_rr->v[catg],tree->mod->ras->gamma_r_proba->v[catg]);
      PhyML_Printf("\n== Pinv = %G",tree->mod->ras->pinvar->v);
      PhyML_Printf("\n== Bl mult = %G",tree->mod->br_len_mult->v);
      PhyML_Printf("\n== fact_sum_scale = %d",fact_sum_scale);
      PhyML_Printf("\n== n_catg: %d",tree->mod->ras->n_catg);
      
      if(tree->mod->whichmodel == GTR || tree->mod->whichmodel == CUSTOM)
        {
          int i,j;
          PhyML_Printf("\n== Rate matrix\n");
          For(i,tree->mod->ns)
            {
              For(j,tree->mod->ns)
                {
                  PhyML_Printf("%12G ",i,tree->mod->r_mat->qmat->v[i*4+j]); 
                }
              PhyML_Printf("\n");
            }
          fflush(NULL);
          
          PhyML_Printf("\n== Relative rates\n");
          For(i,tree->mod->ns*(tree->mod->ns-1)/2)
            {
              PhyML_Printf("\n== rr[%3d]: %12G %12G",i,tree->mod->r_mat->rr->v[i],tree->mod->r_mat->rr_val->v[i]); 
            }
          fflush(NULL);
          
        }
      
      /* int i; */
      /* For(i,2*tree->n_otu-3) */
      /* 	{ */
      /* 	  PhyML_Printf("\n. b%3d->l->v = %f %f [%G] %f %f %f %f [%s]",i, */
      /* 		       tree->a_edges[i]->l->v, */
      /* 		       tree->a_edges[i]->gamma_prior_mean, */
      /* 		       tree->a_edges[i]->gamma_prior_var, */
      /* 		       tree->rates->nd_t[tree->a_edges[i]->left->num], */
      /* 		       tree->rates->nd_t[tree->a_edges[i]->rght->num], */
      /* 		       tree->rates->br_r[tree->a_edges[i]->left->num], */
      /* 		       tree->rates->br_r[tree->a_edges[i]->rght->num], */
      /* 		       tree->a_edges[i] == tree->e_root ? "YES" : "NO"); */
      /* 	  fflush(NULL); */
      
      /* 	} */
      
      Print_Site(tree->data,site,tree->n_otu,"\n",tree->mod->io->state_len,stdout);
      PhyML_Printf("\n== Err. in file %s at line %d (function '%s')",__FILE__,__LINE__,__FUNCTION__);
      Exit("\n");
    }

/* Multiply log likelihood by the number of times this site pattern is found in the data */
  tree->c_lnL_sorted[site] = tree->data->wght[site]*log_site_lk;
  
  tree->c_lnL += tree->data->wght[site]*log_site_lk;


  return log_site_lk;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* phydbl Lk_At_Given_Edge(t_edge *b_fcus, t_tree *tree) */
/* { */
/*   int n_patterns; */

/*   if(tree->is_mixt_tree) return MIXT_Lk_At_Given_Edge(tree); */

/*   n_patterns = tree->n_pattern; */

/* #ifdef PHYTIME */
/*   if((tree->rates) && (tree->rates->bl_from_rt)) RATES_Update_Cur_Bl(tree); */
/* #endif */

/*   Check_Br_Len_Bounds(tree); */

/*   if(tree->rates && tree->io->lk_approx == NORMAL) */
/*     { */
/*       tree->c_lnL = Lk_Normal_Approx(tree); */
/*       return tree->c_lnL; */
/*     } */

/*   Update_PMat_At_Given_Edge(b_fcus,tree); */

/*   if(b_fcus->left->tax) */
/*     { */
/*       PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__); */
/*       Warn_And_Exit(""); */
/*     } */

/*   tree->c_lnL             = .0; */
/*   tree->sum_min_sum_scale = .0; */
/*   For(tree->curr_site,n_patterns) */
/*     if(tree->data->wght[tree->curr_site] > SMALL)  */
/*       Lk_Core(b_fcus,tree); */

/* /\*   Qksort(tree->c_lnL_sorted,NULL,0,n_patterns-1); *\/ */

/* /\*   tree->c_lnL = .0; *\/ */
/* /\*   For(tree->curr_site,n_patterns) *\/ */
/* /\*     { *\/ */
/* /\*       tree->c_lnL += tree->c_lnL_sorted[tree->curr_site]; *\/ */
/* /\*     } *\/ */


/*   Adjust_Min_Diff_Lk(tree); */

/*   return tree->c_lnL; */
/* } */

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Rate_Correction(int exponent, phydbl *site_lk_cat, t_tree *tree)
{
  int piecewise_exponent;
  phydbl multiplier;

  if(exponent >= 0)
    {
      /*! Multiply by 2^exponent */
      do
        {
          piecewise_exponent = MIN(exponent,63);
          multiplier = (phydbl)((unsigned long long)(1) << piecewise_exponent);
          (*site_lk_cat) *= multiplier;
          exponent -= piecewise_exponent;
        }
      while(exponent != 0);
    }
  else
    {
      do
        {
          piecewise_exponent = MAX(exponent,-63);
          multiplier = 1. / (phydbl)((unsigned long long)(1) << -piecewise_exponent);
          (*site_lk_cat) *= multiplier;
          exponent -= piecewise_exponent;
        }
      while(exponent != 0);
    }

  if(isinf(*site_lk_cat))
    {
      PhyML_Printf("\n== Numerical precision issue alert.");
      PhyML_Printf("\n== exponent: %d site_lk_cat: %f", exponent,*site_lk_cat);
      PhyML_Printf("\n== File %s at line %d\n\n",__FILE__,__LINE__);
      (*site_lk_cat) = BIG / 10;
    }
  
  if(*site_lk_cat < SMALL)
    {
      if(tree->mod->ras->n_catg == 1)
        {
          PhyML_Printf("\n== site_lk_cat = %G",*site_lk_cat);
          PhyML_Printf("\n== exponent: %d", exponent);
          PhyML_Printf("\n== Numerical precision issue alert.");
          PhyML_Printf("\n== File %s at line %d (funtction '%s')\n",__FILE__,__LINE__,__FUNCTION__);
	  PhyML_Printf("\n== %s",Write_Tree(tree,NO));
          Exit("\n");
        }
      (*site_lk_cat) = .0;
    }
}

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
              PhyML_Printf("\n== fact_sum_scale: %d",fact_sum_scale);              
              PhyML_Printf("\n== pi: %f",tree->mod->e_frq->pi->v[tree->data->invar[site]]);              
              For(i,tree->mod->ns) PhyML_Printf("\n== pi %d: %f",i,tree->mod->e_frq->pi->v[i]);
              PhyML_Printf("\n== Numerical precision issue alert.");
              PhyML_Printf("\n== File %s at line %d (function '%s')\n",__FILE__,__LINE__,__FUNCTION__);
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

void Update_P_Lk(t_tree *tree, t_edge *b, t_node *d)
{

  if(tree->is_mixt_tree)
    {
      MIXT_Update_P_Lk(tree,b,d);
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
      if(tree->io->datatype == NT)
        {
          Update_P_Lk_Nucl(tree,b,d);
        }
      else if(tree->io->datatype == AA)
        {
          Update_P_Lk_AA(tree,b,d);
        }
      else
        {
          Update_P_Lk_Generic(tree,b,d);
        }
    }
  else
    {
      Update_P_Lk_Generic(tree,b,d);
    }
#endif
//  Print_Edge_Likelihoods(tree, b, false);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

#ifndef BEAGLE
void Update_P_Lk_Generic(t_tree *tree, t_edge *b, t_node *d)
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
  int *sum_scale, *sum_scale_v1, *sum_scale_v2;
  int sum_scale_v1_val, sum_scale_v2_val;
  int i,j;
  int catg,site;
  int n_patterns;
  short int ambiguity_check_v1,ambiguity_check_v2;
  int state_v1,state_v2;
  int NsNg, Ns, NsNs;
  phydbl curr_scaler;
  int curr_scaler_pow, piecewise_scaler_pow;
  phydbl p_lk_lim_inf;
  phydbl smallest_p_lk;
  int *p_lk_loc;

  p_lk_lim_inf = (phydbl)P_LK_LIM_INF;

  NsNg = tree->mod->ras->n_catg * tree->mod->ns;
  Ns = tree->mod->ns;
  NsNs = tree->mod->ns * tree->mod->ns;

  state_v1 = state_v2 = -1;
  ambiguity_check_v1 = ambiguity_check_v2 = NO;
  sum_scale_v1_val = sum_scale_v2_val = 0;
  p1_lk1 = p2_lk2 = .0;

  if(d->tax)
    {
      PhyML_Printf("\n== t_node %d is a leaf...",d->num);
      PhyML_Printf("\n== Err. in file %s at line %d (function '%s').\n",__FILE__,__LINE__,__FUNCTION__);
      Exit("\n");
    }

  n_patterns = tree->n_pattern;

  n_v1 = n_v2                 = NULL;
  p_lk = p_lk_v1 = p_lk_v2    = NULL;
  Pij1 = Pij2                 = NULL;
  sum_scale_v1 = sum_scale_v2 = NULL;
  p_lk_loc                    = NULL;

  Set_All_P_Lk(&n_v1,&n_v2,
               &p_lk,&sum_scale,&p_lk_loc,
               &Pij1,&p_lk_v1,&sum_scale_v1,
               &Pij2,&p_lk_v2,&sum_scale_v2,
               d,b,tree);

  /* For every site in the alignment */
  For(site,n_patterns)
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
              /* if(ambiguity_check_v1 == NO) state_v1 = Get_State_From_P_Pars(n_v1->b[0]->p_lk_tip_r,site*Ns,tree); */
              if(ambiguity_check_v1 == NO) state_v1 = n_v1->c_seq->d_state[site];
            }

          if(n_v2 && n_v2->tax)
            {
              /* Is the state at this tip ambiguous? */
              ambiguity_check_v2 = n_v2->c_seq->is_ambigu[site];
              /* ambiguity_check_v2 = tree->data->c_seq[n_v2->num]->is_ambigu[site]; */
              /* if(ambiguity_check_v2 == NO) state_v2 = Get_State_From_P_Pars(n_v2->b[0]->p_lk_tip_r,site*Ns,tree); */
              if(ambiguity_check_v2 == NO) state_v2 = n_v2->c_seq->d_state[site];
            }
        }

      if(tree->mod->use_m4mod)
        {
          ambiguity_check_v1 = YES;
          ambiguity_check_v2 = YES;
        }

      if(p_lk_loc[site] < site) /* Have we seen this pattern before? */
        {
          Copy_P_Lk(p_lk,p_lk_loc[site],site,tree);
          Copy_Scale(sum_scale,p_lk_loc[site],site,tree);
        }
      else
        {
          /* For all the rate classes */
          For(catg,tree->mod->ras->n_catg)
            {
              if(tree->mod->ras->skip_rate_cat[catg] == YES) continue;

              smallest_p_lk  =  BIG;

              /* For all the states at node d */
              For(i,tree->mod->ns)
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
                              p1_lk1 = Pij1[catg*NsNs+i*Ns+state_v1];
                            }
                          else
                            {
                              /* For all the states at node n_v1 */
                              For(j,tree->mod->ns)
                                {
                                  p1_lk1 += Pij1[catg*NsNs+i*Ns+j] * (phydbl)n_v1->b[0]->p_lk_tip_r[site*Ns+j];
                                }
                            }
                        }
                      /* n_v1 is an internal node */
                      else
                        {
                          /* For the states at node n_v1 */
                          For(j,tree->mod->ns)
                            {
                              p1_lk1 += Pij1[catg*NsNs+i*Ns+j] * p_lk_v1[site*NsNg+catg*Ns+j];
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
                              p2_lk2 = Pij2[catg*NsNs+i*Ns+state_v2];
                            }
                          else
                            {
                              /* For all the states at node n_v2 */
                              For(j,tree->mod->ns)
                                {
                                  p2_lk2 += Pij2[catg*NsNs+i*Ns+j] * (phydbl)n_v2->b[0]->p_lk_tip_r[site*Ns+j];
                                }
                            }
                        }
                      /* n_v2 is an internal node */
                      else
                        {
                          /* For all the states at node n_v2 */
                          For(j,tree->mod->ns)
                            {
                              p2_lk2 += Pij2[catg*NsNs+i*Ns+j] * p_lk_v2[site*NsNg+catg*Ns+j];
                            }
                        }
                    }
                  else
                    {
                      p2_lk2 = 1.0;
                    }

                  p_lk[site*NsNg+catg*Ns+i] = p1_lk1 * p2_lk2;

                  /* 	      PhyML_Printf("\n+ %G",p_lk[site*NsNg+catg*Ns+i]); */

                  if(p_lk[site*NsNg+catg*Ns+i] < smallest_p_lk) smallest_p_lk = p_lk[site*NsNg+catg*Ns+i] ;
                }

              /* Current scaling values at that site */
              sum_scale_v1_val = (sum_scale_v1)?(sum_scale_v1[catg*n_patterns+site]):(0);
              sum_scale_v2_val = (sum_scale_v2)?(sum_scale_v2[catg*n_patterns+site]):(0);

              sum_scale[catg*n_patterns+site] = sum_scale_v1_val + sum_scale_v2_val;

              /* 	  PhyML_Printf("\n@ %d",sum_scale[catg*n_patterns+site]); */

              /* Scaling. For details read utilities.h */
              if(smallest_p_lk < p_lk_lim_inf)
                {
                  curr_scaler_pow = (int)(LOG(p_lk_lim_inf)-LOG(smallest_p_lk))/LOG2;
                  curr_scaler     = (phydbl)((unsigned long long)(1) << curr_scaler_pow);

                  /* 	      if(fabs(curr_scaler_pow) > 63 || fabs(curr_scaler_pow) > 63) */
                  /* 		{ */
                  /* 		  PhyML_Printf("\n. p_lk_lim_inf = %G smallest_p_lk = %G",p_lk_lim_inf,smallest_p_lk); */
                  /* 		  PhyML_Printf("\n. curr_scaler_pow = %d",curr_scaler_pow); */
                  /* 		  PhyML_Printf("\n. Err in file %s at line %d.",__FILE__,__LINE__); */
                  /* 		  Warn_And_Exit("\n"); */
                  /* 		} */

                  sum_scale[catg*n_patterns+site] += curr_scaler_pow;

                  do
                    {
                      piecewise_scaler_pow = MIN(curr_scaler_pow,63);
                      curr_scaler = (phydbl)((unsigned long long)(1) << piecewise_scaler_pow);
                      For(i,tree->mod->ns)
                        {
                          p_lk[site*NsNg+catg*Ns+i] *= curr_scaler;

                          if(p_lk[site*NsNg+catg*Ns+i] > BIG)
                            {
                              PhyML_Printf("\n. p_lk_lim_inf = %G smallest_p_lk = %G",p_lk_lim_inf,smallest_p_lk);
                              PhyML_Printf("\n. curr_scaler_pow = %d",curr_scaler_pow);
                              PhyML_Printf("\n. Err. in file %s at line %d (function '%s').",__FILE__,__LINE__,__FUNCTION__);
                              Warn_And_Exit("\n");
                            }
                        }
                      curr_scaler_pow -= piecewise_scaler_pow;
                    }
                  while(curr_scaler_pow != 0);
                }
            }
        }
    }
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Update_P_Lk_Nucl(t_tree *tree, t_edge *b, t_node *d)
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
  phydbl *p_lk,*p_lk_v1,*p_lk_v2;//Partial likelihood vector of node d, d's "left" neighbor, d's "right" neighbor. We fill *p_lk, and assume *p_lk_v1 and *p_lk_v2 are already filled.
  phydbl *Pij1,*Pij2;
  int *sum_scale, *sum_scale_v1, *sum_scale_v2;
  int sum_scale_v1_val, sum_scale_v2_val;
  int i;//index over the number of states
  int catg/*index over the number of rate categories*/,site;
  int n_patterns;//Number of distinct site patterns
  short int ambiguity_check_v1,ambiguity_check_v2;
  int state_v1,state_v2;
  int dim1, dim2, dim3;
  phydbl curr_scaler;
  int curr_scaler_pow, piecewise_scaler_pow;
  phydbl p_lk_lim_inf;
  phydbl smallest_p_lk;
  phydbl p0,p1,p2,p3;
  int *p_lk_loc;//Suppose site j, of a certain subtree, has "A" on one tip, and "C" on the other. If you come across this pattern again at site i<j, then you can simply copy the partial likelihoods


  p_lk_lim_inf = (phydbl)P_LK_LIM_INF;

  dim1 = tree->mod->ras->n_catg * tree->mod->ns;//Dimension of a matrix L that holds rate-specific character likelihoods. IOW, L[ij] is the likelihood of character j in rate class i
  dim2 = tree->mod->ns;
  dim3 = tree->mod->ns * tree->mod->ns;//Dimensions of the transition prob. matrix

  state_v1 = state_v2 = -1;
  ambiguity_check_v1 = ambiguity_check_v2 = NO;
  sum_scale_v1_val = sum_scale_v2_val = 0;
  p1_lk1 = p2_lk2 = .0;
  curr_scaler = .0;
  curr_scaler_pow = piecewise_scaler_pow = 0;

  if(d->tax)
    {
      PhyML_Printf("\n== t_node %d is a leaf...",d->num);
      PhyML_Printf("\n== Err. in file %s at line %d (function '%s')\n",__FILE__,__LINE__,__FUNCTION__);
      Exit("\n");
    }

  n_patterns = tree->n_pattern;

  n_v1 = n_v2                 = NULL;
  p_lk = p_lk_v1 = p_lk_v2    = NULL;
  Pij1 = Pij2                 = NULL;
  sum_scale_v1 = sum_scale_v2 = NULL;
  p_lk_loc                    = NULL;
  Set_All_P_Lk(&n_v1,&n_v2,
               &p_lk,&sum_scale,&p_lk_loc,
               &Pij1,&p_lk_v1,&sum_scale_v1,
               &Pij2,&p_lk_v2,&sum_scale_v2,
               d,b,tree);


  /* For every site in the alignment */
  For(site,n_patterns)
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


      if(p_lk_loc[site] < site)
        {
          Copy_P_Lk(p_lk,p_lk_loc[site],site,tree);
          Copy_Scale(sum_scale,p_lk_loc[site],site,tree);
        }
      else
        {
          /* For all the rate classes */
          For(catg,tree->mod->ras->n_catg)
            {
              smallest_p_lk  =  BIG;

              /* For all states at node d */
              For(i,tree->mod->ns)
                {
                  if(tree->mod->augmented == YES)
                    {
                      For(i,tree->mod->ns) p_lk[site*dim1+catg*dim2+i] = 0.0;
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
                              p1_lk1 = Pij1[catg*dim3+i*dim2+state_v1];

                              assert(!isnan(p1_lk1));

                              /* if(isnan(p1_lk1)) */
                              /*   { */
                              /*     PhyML_Printf("\n== Tree %d",tree->tree_num); */
                              /*     PhyML_Printf("\n== catg=%d dim3=%d dim2=%d i=%d state_v1=%d",catg,dim3,dim2,i,state_v1); */
                              /*     PhyML_Printf("\n== Pij1[0] = %G l = %G",Pij1[0],b->l->v); */
                              /*     Print_Model(tree->mod); */
                              /*     Generic_Exit(__FILE__,__LINE__,__FUNCTION__); */
                              /*   } */
                            }
                          else
                            {
                              /* For all the states at node n_v1 */
                              p0=Pij1[catg*dim3+i*dim2+0] * (phydbl)n_v1->b[0]->p_lk_tip_r[site*dim2+0];
                              p1=Pij1[catg*dim3+i*dim2+1] * (phydbl)n_v1->b[0]->p_lk_tip_r[site*dim2+1];
                              p2=Pij1[catg*dim3+i*dim2+2] * (phydbl)n_v1->b[0]->p_lk_tip_r[site*dim2+2];
                              p3=Pij1[catg*dim3+i*dim2+3] * (phydbl)n_v1->b[0]->p_lk_tip_r[site*dim2+3];
                              p1_lk1 = p0+p1+p2+p3;

                              assert(!isnan(p1_lk1));

                              /* if(isnan(p1_lk1)) */
                              /*   { */
                              /*     PhyML_Printf("\n== p0=%f p1=%f p2=%f p3=%f",p0,p1,p2,p3); */
                              /*     PhyML_Printf("\n== Err. in file %s at line %d (function '%s').",__FILE__,__LINE__,__FUNCTION__); */
                              /*     Exit("\n"); */
                              /*   } */
                            }
                        }
                      /* n_v1 is an internal node */
                      else
                        {
                          //"catg*dim3" offsets into a rate-specific flat 4*4 DNA matrix P. Then, P[i,j] is given by "i*4+l"
                          //"site*dim1" offsets into a rate-specific flat num_rates*4 matrix L. Then, L[i,j] is given by "catg*dim2+l"
                          if(tree->mod->augmented == YES)
                            {
                              p1_lk1 = Pij1[catg*dim3+i*dim2+state_v1] * p_lk_v1[site*dim1+catg*dim2+state_v1];
                            }
                          else
                            {
                              p0=Pij1[catg*dim3+i*dim2+0] * p_lk_v1[site*dim1+catg*dim2+0];
                              p1=Pij1[catg*dim3+i*dim2+1] * p_lk_v1[site*dim1+catg*dim2+1];
                              p2=Pij1[catg*dim3+i*dim2+2] * p_lk_v1[site*dim1+catg*dim2+2];
                              p3=Pij1[catg*dim3+i*dim2+3] * p_lk_v1[site*dim1+catg*dim2+3];
                              p1_lk1 = p0+p1+p2+p3;
                            }
                          if(isnan(p1_lk1))
                            {
                              /* PhyML_Printf("\n== p0=%f p1=%f p2=%f p3=%f",p0,p1,p2,p3); */
                              PhyML_Printf("\n== Err. in file %s at line %d.",__FILE__,__LINE__);
                              Warn_And_Exit("\n");
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
                              p2_lk2 = Pij2[catg*dim3+i*dim2+state_v2];

                              assert(!isnan(p2_lk2));

                              /* if(isnan(p2_lk2)) */
                              /*   { */
                              /*     PhyML_Printf("\n== Tree %d",tree->tree_num); */
                              /*     PhyML_Printf("\n== catg=%d dim3=%d dim2=%d i=%d state_v2=%d",catg,dim3,dim2,i,state_v2); */
                              /*     PhyML_Printf("\n== Pij2[0] = %G l = %G",Pij2[0],b->l->v); */
                              /*     Print_Model(tree->mod); */
                              /*     Generic_Exit(__FILE__,__LINE__,__FUNCTION__); */
                              /*   } */
                            }
                          else
                            {
                              /* For all the states at node n_v2 */
                              p0=Pij2[catg*dim3+i*dim2+0] * (phydbl)n_v2->b[0]->p_lk_tip_r[site*dim2+0];
                              p1=Pij2[catg*dim3+i*dim2+1] * (phydbl)n_v2->b[0]->p_lk_tip_r[site*dim2+1];
                              p2=Pij2[catg*dim3+i*dim2+2] * (phydbl)n_v2->b[0]->p_lk_tip_r[site*dim2+2];
                              p3=Pij2[catg*dim3+i*dim2+3] * (phydbl)n_v2->b[0]->p_lk_tip_r[site*dim2+3];
                              p2_lk2 = p0+p1+p2+p3;

                              assert(!isnan(p2_lk2));

                              /* if(isnan(p2_lk2)) */
                              /*   { */
                              /*     PhyML_Printf("\n== p0=%f p1=%f p2=%f p3=%f",p0,p1,p2,p3); */
                              /*     PhyML_Printf("\n== Err. in file %s at line %d.",__FILE__,__LINE__); */
                              /*     Warn_And_Exit("\n"); */
                              /*   }                              */
                            }
                        }
                      /* n_v2 is an internal node */
                      else
                        {
                          if(tree->mod->augmented == YES)
                            {
                              p2_lk2 = Pij2[catg*dim3+i*dim2+state_v2] * p_lk_v2[site*dim1+catg*dim2+state_v2];
                            }
                          else
                            {
                              /* For all the states at node n_v2 */
                              p0=Pij2[catg*dim3+i*dim2+0] * p_lk_v2[site*dim1+catg*dim2+0];
                              p1=Pij2[catg*dim3+i*dim2+1] * p_lk_v2[site*dim1+catg*dim2+1];
                              p2=Pij2[catg*dim3+i*dim2+2] * p_lk_v2[site*dim1+catg*dim2+2];
                              p3=Pij2[catg*dim3+i*dim2+3] * p_lk_v2[site*dim1+catg*dim2+3];
                              p2_lk2 = p0+p1+p2+p3;
                            }

                          if(isnan(p2_lk2))
                            {
                              PhyML_Printf("\n== %f %f",b->l->v,b->l->v*b->l->v*tree->mod->l_var_sigma);
                              /* PhyML_Printf("\n== p0=%f p1=%f p2=%f p3=%f",p0,p1,p2,p3); */
                              PhyML_Printf("\n== Pij2[0]=%f Pij2[1]=%f Pij2[2]=%f Pij2[3]=%f",Pij2[catg*dim3+i*dim2+0],Pij2[catg*dim3+i*dim2+1],Pij2[catg*dim3+i*dim2+2],Pij2[catg*dim3+i*dim2+3]);
                              PhyML_Printf("\n== Err. in file %s at line %d (function '%s').",__FILE__,__LINE__,__FUNCTION__);
                              Warn_And_Exit("\n");
                            }
                        }
                    }
                  else
                    {
                      p2_lk2 = 1.0;
                    }
                  
                  /* Partial likelihood of character "i" at site "site" under rate "catg" */
                  p_lk[site*dim1+catg*dim2+i] = p1_lk1 * p2_lk2;
                  
                  if(p_lk[site*dim1+catg*dim2+i] < smallest_p_lk) smallest_p_lk = p_lk[site*dim1+catg*dim2+i] ;

                  if(tree->mod->augmented == YES) break;

                }
              
              /* Current scaling values at that site */
              sum_scale_v1_val = (sum_scale_v1)?(sum_scale_v1[catg*n_patterns+site]):(0);
              sum_scale_v2_val = (sum_scale_v2)?(sum_scale_v2[catg*n_patterns+site]):(0);
              
              sum_scale[catg*n_patterns+site] = sum_scale_v1_val + sum_scale_v2_val;
              
              /* Scaling */
              if(smallest_p_lk < p_lk_lim_inf && tree->mod->augmented == NO)
                {
                  curr_scaler_pow = (int)(LOG(p_lk_lim_inf)-LOG(smallest_p_lk))/LOG2;
                  curr_scaler     = (phydbl)((unsigned long long)(1) << curr_scaler_pow);
                  
                  sum_scale[catg*n_patterns+site] += curr_scaler_pow;
                  //                  fprintf(stderr,"\n%d",sum_scale[catg*n_patterns+site]);
                                        
                  do
                    {
                      piecewise_scaler_pow = MIN(curr_scaler_pow,63);
                      curr_scaler = (phydbl)((unsigned long long)(1) << piecewise_scaler_pow);
                      For(i,tree->mod->ns)
                        {
                          p_lk[site*dim1+catg*dim2+i] *= curr_scaler;

                          if(p_lk[site*dim1+catg*dim2+i] > BIG)
                            {
                              PhyML_Printf("\n== p_lk_lim_inf = %G smallest_p_lk = %G",p_lk_lim_inf,smallest_p_lk);
                              PhyML_Printf("\n== curr_scaler_pow = %d",curr_scaler_pow);
                              PhyML_Printf("\n== Err. in file %s at line %d (function '%s').",__FILE__,__LINE__,__FUNCTION__);
                              Exit("\n");
                            }
//                              fprintf(stderr,"\n%e",p_lk[site*dim1+catg*dim2+i]);
                        }
                      curr_scaler_pow -= piecewise_scaler_pow;
                    }
                  while(curr_scaler_pow != 0);
                }
            }
        }
    }

//  fprintf(stdout, "Updated partials:");fflush(stdout);
//  Dump_Arr_D(p_lk, tree->mod->ras->n_catg*tree->mod->ns*tree->n_pattern);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Update_P_Lk_AA(t_tree *tree, t_edge *b, t_node *d)
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
  int *sum_scale, *sum_scale_v1, *sum_scale_v2;
  int sum_scale_v1_val, sum_scale_v2_val;
  int i;
  int catg,site;
  int n_patterns;
  short int ambiguity_check_v1,ambiguity_check_v2;
  int state_v1,state_v2;
  int dim1, dim2, dim3;
  phydbl curr_scaler;
  int curr_scaler_pow, piecewise_scaler_pow;
  phydbl p_lk_lim_inf;
  phydbl smallest_p_lk;
  phydbl p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19;
  int *p_lk_loc;

  p_lk_lim_inf = (phydbl)P_LK_LIM_INF;

  dim1 = tree->mod->ras->n_catg * tree->mod->ns;
  dim2 = tree->mod->ns;
  dim3 = tree->mod->ns * tree->mod->ns;

  state_v1 = state_v2 = -1;
  ambiguity_check_v1 = ambiguity_check_v2 = NO;
  sum_scale_v1_val = sum_scale_v2_val = 0;
  p1_lk1 = p2_lk2 = .0;

  if(d->tax)
    {
      PhyML_Printf("\n== t_node %d is a leaf...",d->num);
      PhyML_Printf("\n== Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
      Warn_And_Exit("\n");
    }

  n_patterns = tree->n_pattern;

  n_v1 = n_v2                 = NULL;
  p_lk = p_lk_v1 = p_lk_v2    = NULL;
  Pij1 = Pij2                 = NULL;
  sum_scale_v1 = sum_scale_v2 = NULL;
  p_lk_loc                    = NULL;

  Set_All_P_Lk(&n_v1,&n_v2,
               &p_lk,&sum_scale,&p_lk_loc,
               &Pij1,&p_lk_v1,&sum_scale_v1,
               &Pij2,&p_lk_v2,&sum_scale_v2,
               d,b,tree);


  /* For every site in the alignment */
  For(site,n_patterns)
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
          /* if(ambiguity_check_v1 == NO) state_v1 = Get_State_From_P_Pars(n_v1->b[0]->p_lk_tip_r,site*dim2,tree); */
              if(ambiguity_check_v1 == NO) state_v1 = n_v1->c_seq->d_state[site];
        }

      if(n_v2 && n_v2->tax)
        {
          /* Is the state at this tip ambiguous? */
          ambiguity_check_v2 = n_v2->c_seq->is_ambigu[site];
          /* if(ambiguity_check_v2 == NO) state_v2 = Get_State_From_P_Pars(n_v2->b[0]->p_lk_tip_r,site*dim2,tree); */
              if(ambiguity_check_v2 == NO) state_v2 = n_v2->c_seq->d_state[site];
        }
    }

      if(tree->mod->use_m4mod)
        {
          ambiguity_check_v1 = YES;
          ambiguity_check_v2 = YES;
        }

      if(p_lk_loc[site] < site)
        {
          Copy_P_Lk(p_lk,p_lk_loc[site],site,tree);
          Copy_Scale(sum_scale,p_lk_loc[site],site,tree);
        }
      else
        {
          /* For all the rate classes */
          For(catg,tree->mod->ras->n_catg)
            {
              smallest_p_lk  =  BIG;
              
              /* For all the state at node d */
              For(i,tree->mod->ns)
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
                              p1_lk1 = Pij1[catg*dim3+i*dim2+state_v1];
                            }
                          else
                            {
                              /* For all the states at node n_v1 */
                              p0  = Pij1[catg*dim3+i*dim2+ 0] * (phydbl)n_v1->b[0]->p_lk_tip_r[site*dim2+ 0];
                              p1  = Pij1[catg*dim3+i*dim2+ 1] * (phydbl)n_v1->b[0]->p_lk_tip_r[site*dim2+ 1];
                              p2  = Pij1[catg*dim3+i*dim2+ 2] * (phydbl)n_v1->b[0]->p_lk_tip_r[site*dim2+ 2];
                              p3  = Pij1[catg*dim3+i*dim2+ 3] * (phydbl)n_v1->b[0]->p_lk_tip_r[site*dim2+ 3];
                              p4  = Pij1[catg*dim3+i*dim2+ 4] * (phydbl)n_v1->b[0]->p_lk_tip_r[site*dim2+ 4];
                              p5  = Pij1[catg*dim3+i*dim2+ 5] * (phydbl)n_v1->b[0]->p_lk_tip_r[site*dim2+ 5];
                              p6  = Pij1[catg*dim3+i*dim2+ 6] * (phydbl)n_v1->b[0]->p_lk_tip_r[site*dim2+ 6];
                              p7  = Pij1[catg*dim3+i*dim2+ 7] * (phydbl)n_v1->b[0]->p_lk_tip_r[site*dim2+ 7];
                              p8  = Pij1[catg*dim3+i*dim2+ 8] * (phydbl)n_v1->b[0]->p_lk_tip_r[site*dim2+ 8];
                              p9  = Pij1[catg*dim3+i*dim2+ 9] * (phydbl)n_v1->b[0]->p_lk_tip_r[site*dim2+ 9];
                              p10 = Pij1[catg*dim3+i*dim2+10] * (phydbl)n_v1->b[0]->p_lk_tip_r[site*dim2+10];
                              p11 = Pij1[catg*dim3+i*dim2+11] * (phydbl)n_v1->b[0]->p_lk_tip_r[site*dim2+11];
                              p12 = Pij1[catg*dim3+i*dim2+12] * (phydbl)n_v1->b[0]->p_lk_tip_r[site*dim2+12];
                              p13 = Pij1[catg*dim3+i*dim2+13] * (phydbl)n_v1->b[0]->p_lk_tip_r[site*dim2+13];
                              p14 = Pij1[catg*dim3+i*dim2+14] * (phydbl)n_v1->b[0]->p_lk_tip_r[site*dim2+14];
                              p15 = Pij1[catg*dim3+i*dim2+15] * (phydbl)n_v1->b[0]->p_lk_tip_r[site*dim2+15];
                              p16 = Pij1[catg*dim3+i*dim2+16] * (phydbl)n_v1->b[0]->p_lk_tip_r[site*dim2+16];
                              p17 = Pij1[catg*dim3+i*dim2+17] * (phydbl)n_v1->b[0]->p_lk_tip_r[site*dim2+17];
                              p18 = Pij1[catg*dim3+i*dim2+18] * (phydbl)n_v1->b[0]->p_lk_tip_r[site*dim2+18];
                              p19 = Pij1[catg*dim3+i*dim2+19] * (phydbl)n_v1->b[0]->p_lk_tip_r[site*dim2+19];
                              p1_lk1 = p0+p1+p2+p3+p4+p5+p6+p7+p8+p9+p10+p11+p12+p13+p14+p15+p16+p17+p18+p19;
                            }
                        }
                      /* n_v1 is an internal node */
                      else
                        {
                          /* For the states at node n_v1 */
                          p0  = Pij1[catg*dim3+i*dim2+ 0] * (phydbl)p_lk_v1[site*dim1+catg*dim2+ 0];
                          p1  = Pij1[catg*dim3+i*dim2+ 1] * (phydbl)p_lk_v1[site*dim1+catg*dim2+ 1];
                          p2  = Pij1[catg*dim3+i*dim2+ 2] * (phydbl)p_lk_v1[site*dim1+catg*dim2+ 2];
                          p3  = Pij1[catg*dim3+i*dim2+ 3] * (phydbl)p_lk_v1[site*dim1+catg*dim2+ 3];
                          p4  = Pij1[catg*dim3+i*dim2+ 4] * (phydbl)p_lk_v1[site*dim1+catg*dim2+ 4];
                          p5  = Pij1[catg*dim3+i*dim2+ 5] * (phydbl)p_lk_v1[site*dim1+catg*dim2+ 5];
                          p6  = Pij1[catg*dim3+i*dim2+ 6] * (phydbl)p_lk_v1[site*dim1+catg*dim2+ 6];
                          p7  = Pij1[catg*dim3+i*dim2+ 7] * (phydbl)p_lk_v1[site*dim1+catg*dim2+ 7];
                          p8  = Pij1[catg*dim3+i*dim2+ 8] * (phydbl)p_lk_v1[site*dim1+catg*dim2+ 8];
                          p9  = Pij1[catg*dim3+i*dim2+ 9] * (phydbl)p_lk_v1[site*dim1+catg*dim2+ 9];
                          p10 = Pij1[catg*dim3+i*dim2+10] * (phydbl)p_lk_v1[site*dim1+catg*dim2+10];
                          p11 = Pij1[catg*dim3+i*dim2+11] * (phydbl)p_lk_v1[site*dim1+catg*dim2+11];
                          p12 = Pij1[catg*dim3+i*dim2+12] * (phydbl)p_lk_v1[site*dim1+catg*dim2+12];
                          p13 = Pij1[catg*dim3+i*dim2+13] * (phydbl)p_lk_v1[site*dim1+catg*dim2+13];
                          p14 = Pij1[catg*dim3+i*dim2+14] * (phydbl)p_lk_v1[site*dim1+catg*dim2+14];
                          p15 = Pij1[catg*dim3+i*dim2+15] * (phydbl)p_lk_v1[site*dim1+catg*dim2+15];
                          p16 = Pij1[catg*dim3+i*dim2+16] * (phydbl)p_lk_v1[site*dim1+catg*dim2+16];
                          p17 = Pij1[catg*dim3+i*dim2+17] * (phydbl)p_lk_v1[site*dim1+catg*dim2+17];
                          p18 = Pij1[catg*dim3+i*dim2+18] * (phydbl)p_lk_v1[site*dim1+catg*dim2+18];
                          p19 = Pij1[catg*dim3+i*dim2+19] * (phydbl)p_lk_v1[site*dim1+catg*dim2+19];
                          p1_lk1 = p0+p1+p2+p3+p4+p5+p6+p7+p8+p9+p10+p11+p12+p13+p14+p15+p16+p17+p18+p19;
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
              p2_lk2 = Pij2[catg*dim3+i*dim2+state_v2];
            }
              else
            {
              /* For all the states at node n_v2 */
              p0  = Pij2[catg*dim3+i*dim2+ 0] * (phydbl)n_v2->b[0]->p_lk_tip_r[site*dim2+ 0];
              p1  = Pij2[catg*dim3+i*dim2+ 1] * (phydbl)n_v2->b[0]->p_lk_tip_r[site*dim2+ 1];
              p2  = Pij2[catg*dim3+i*dim2+ 2] * (phydbl)n_v2->b[0]->p_lk_tip_r[site*dim2+ 2];
              p3  = Pij2[catg*dim3+i*dim2+ 3] * (phydbl)n_v2->b[0]->p_lk_tip_r[site*dim2+ 3];
              p4  = Pij2[catg*dim3+i*dim2+ 4] * (phydbl)n_v2->b[0]->p_lk_tip_r[site*dim2+ 4];
              p5  = Pij2[catg*dim3+i*dim2+ 5] * (phydbl)n_v2->b[0]->p_lk_tip_r[site*dim2+ 5];
              p6  = Pij2[catg*dim3+i*dim2+ 6] * (phydbl)n_v2->b[0]->p_lk_tip_r[site*dim2+ 6];
              p7  = Pij2[catg*dim3+i*dim2+ 7] * (phydbl)n_v2->b[0]->p_lk_tip_r[site*dim2+ 7];
              p8  = Pij2[catg*dim3+i*dim2+ 8] * (phydbl)n_v2->b[0]->p_lk_tip_r[site*dim2+ 8];
              p9  = Pij2[catg*dim3+i*dim2+ 9] * (phydbl)n_v2->b[0]->p_lk_tip_r[site*dim2+ 9];
              p10 = Pij2[catg*dim3+i*dim2+10] * (phydbl)n_v2->b[0]->p_lk_tip_r[site*dim2+10];
              p11 = Pij2[catg*dim3+i*dim2+11] * (phydbl)n_v2->b[0]->p_lk_tip_r[site*dim2+11];
              p12 = Pij2[catg*dim3+i*dim2+12] * (phydbl)n_v2->b[0]->p_lk_tip_r[site*dim2+12];
              p13 = Pij2[catg*dim3+i*dim2+13] * (phydbl)n_v2->b[0]->p_lk_tip_r[site*dim2+13];
              p14 = Pij2[catg*dim3+i*dim2+14] * (phydbl)n_v2->b[0]->p_lk_tip_r[site*dim2+14];
              p15 = Pij2[catg*dim3+i*dim2+15] * (phydbl)n_v2->b[0]->p_lk_tip_r[site*dim2+15];
              p16 = Pij2[catg*dim3+i*dim2+16] * (phydbl)n_v2->b[0]->p_lk_tip_r[site*dim2+16];
              p17 = Pij2[catg*dim3+i*dim2+17] * (phydbl)n_v2->b[0]->p_lk_tip_r[site*dim2+17];
              p18 = Pij2[catg*dim3+i*dim2+18] * (phydbl)n_v2->b[0]->p_lk_tip_r[site*dim2+18];
              p19 = Pij2[catg*dim3+i*dim2+19] * (phydbl)n_v2->b[0]->p_lk_tip_r[site*dim2+19];
              p2_lk2 = p0+p1+p2+p3+p4+p5+p6+p7+p8+p9+p10+p11+p12+p13+p14+p15+p16+p17+p18+p19;

            }
            }
          /* n_v2 is an internal node */
          else
            {
              /* For all the states at node n_v2 */
              p0  = Pij2[catg*dim3+i*dim2+ 0] * (phydbl)p_lk_v2[site*dim1+catg*dim2+ 0];
              p1  = Pij2[catg*dim3+i*dim2+ 1] * (phydbl)p_lk_v2[site*dim1+catg*dim2+ 1];
              p2  = Pij2[catg*dim3+i*dim2+ 2] * (phydbl)p_lk_v2[site*dim1+catg*dim2+ 2];
              p3  = Pij2[catg*dim3+i*dim2+ 3] * (phydbl)p_lk_v2[site*dim1+catg*dim2+ 3];
              p4  = Pij2[catg*dim3+i*dim2+ 4] * (phydbl)p_lk_v2[site*dim1+catg*dim2+ 4];
              p5  = Pij2[catg*dim3+i*dim2+ 5] * (phydbl)p_lk_v2[site*dim1+catg*dim2+ 5];
              p6  = Pij2[catg*dim3+i*dim2+ 6] * (phydbl)p_lk_v2[site*dim1+catg*dim2+ 6];
              p7  = Pij2[catg*dim3+i*dim2+ 7] * (phydbl)p_lk_v2[site*dim1+catg*dim2+ 7];
              p8  = Pij2[catg*dim3+i*dim2+ 8] * (phydbl)p_lk_v2[site*dim1+catg*dim2+ 8];
              p9  = Pij2[catg*dim3+i*dim2+ 9] * (phydbl)p_lk_v2[site*dim1+catg*dim2+ 9];
              p10 = Pij2[catg*dim3+i*dim2+10] * (phydbl)p_lk_v2[site*dim1+catg*dim2+10];
              p11 = Pij2[catg*dim3+i*dim2+11] * (phydbl)p_lk_v2[site*dim1+catg*dim2+11];
              p12 = Pij2[catg*dim3+i*dim2+12] * (phydbl)p_lk_v2[site*dim1+catg*dim2+12];
              p13 = Pij2[catg*dim3+i*dim2+13] * (phydbl)p_lk_v2[site*dim1+catg*dim2+13];
              p14 = Pij2[catg*dim3+i*dim2+14] * (phydbl)p_lk_v2[site*dim1+catg*dim2+14];
              p15 = Pij2[catg*dim3+i*dim2+15] * (phydbl)p_lk_v2[site*dim1+catg*dim2+15];
              p16 = Pij2[catg*dim3+i*dim2+16] * (phydbl)p_lk_v2[site*dim1+catg*dim2+16];
              p17 = Pij2[catg*dim3+i*dim2+17] * (phydbl)p_lk_v2[site*dim1+catg*dim2+17];
              p18 = Pij2[catg*dim3+i*dim2+18] * (phydbl)p_lk_v2[site*dim1+catg*dim2+18];
              p19 = Pij2[catg*dim3+i*dim2+19] * (phydbl)p_lk_v2[site*dim1+catg*dim2+19];
              p2_lk2 = p0+p1+p2+p3+p4+p5+p6+p7+p8+p9+p10+p11+p12+p13+p14+p15+p16+p17+p18+p19;
                        }
                    }
                  else
                    {
                      p2_lk2 = 1.0;
            }

          p_lk[site*dim1+catg*dim2+i] = p1_lk1 * p2_lk2;

          if(p_lk[site*dim1+catg*dim2+i] < smallest_p_lk) smallest_p_lk = p_lk[site*dim1+catg*dim2+i] ;
        }

          /* Current scaling values at that site */
          sum_scale_v1_val = (sum_scale_v1)?(sum_scale_v1[catg*n_patterns+site]):(0);
          sum_scale_v2_val = (sum_scale_v2)?(sum_scale_v2[catg*n_patterns+site]):(0);

          sum_scale[catg*n_patterns+site] = sum_scale_v1_val + sum_scale_v2_val;

          /* Scaling */
          if(smallest_p_lk < p_lk_lim_inf)
            {
              curr_scaler_pow = (int)(LOG(p_lk_lim_inf)-LOG(smallest_p_lk))/LOG2;
              curr_scaler     = (phydbl)((unsigned long long)(1) << curr_scaler_pow);

              /* 	      if(fabs(curr_scaler_pow) > 63 || fabs(curr_scaler_pow) > 63) */
              /* 		{ */
              /* 		  PhyML_Printf("\n. p_lk_lim_inf = %G smallest_p_lk = %G",p_lk_lim_inf,smallest_p_lk); */
              /* 		  PhyML_Printf("\n. curr_scaler_pow = %d",curr_scaler_pow); */
              /* 		  PhyML_Printf("\n. Err in file %s at line %d.",__FILE__,__LINE__); */
              /* 		  Warn_And_Exit("\n"); */
              /* 		} */

              sum_scale[catg*n_patterns+site] += curr_scaler_pow;

              do
                {
                  piecewise_scaler_pow = MIN(curr_scaler_pow,63);
                  curr_scaler = (phydbl)((unsigned long long)(1) << piecewise_scaler_pow);
                  For(i,tree->mod->ns)
                {
                  p_lk[site*dim1+catg*dim2+i] *= curr_scaler;

                  if(p_lk[site*dim1+catg*dim2+i] > BIG)
                    {
                      PhyML_Printf("\n. p_lk_lim_inf = %G smallest_p_lk = %G",p_lk_lim_inf,smallest_p_lk);
                      PhyML_Printf("\n. curr_scaler_pow = %d",curr_scaler_pow);
                      PhyML_Printf("\n. Err in file %s at line %d.",__FILE__,__LINE__);
                      Warn_And_Exit("\n");
                    }
                }
                  curr_scaler_pow -= piecewise_scaler_pow;
                }
              while(curr_scaler_pow != 0);
            }
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
  int state0, state1,len;
  phydbl *F;
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

  For(i,mod->ras->n_catg) /* Don't use the discrete gamma distribution */
    {
      mod->ras->gamma_rr->v[i]      = 1.0;
      mod->ras->gamma_r_proba->v[i] = 1.0;
    }

  n_catg = mod->ras->n_catg;
  mod->ras->n_catg = 1;

  For(j,data->n_otu-1)
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
          len = 0;
          For(l,twodata->c_seq[0]->len)
            {
              state0 = Assign_State(twodata->c_seq[0]->state+l*mod->io->state_len,mod->io->datatype,mod->io->state_len);
              state1 = Assign_State(twodata->c_seq[1]->state+l*mod->io->state_len,mod->io->datatype,mod->io->state_len);
              
              if((state0 > -1) && (state1 > -1))
                {
                  F[mod->ns*state0+state1] += twodata->wght[l];
                  len += (int)twodata->wght[l];
                }
            }
          if(len > .0) {For(i,mod->ns*mod->ns) F[i] /= (phydbl)len;}
          
          sum = 0.;
          For(i,mod->ns*mod->ns) sum += F[i];
          
          
          /* if(sum < .001) d_max = -1.; */
          if(sum < .001) d_max = init;
          else if((sum > 1. - .001) && (sum < 1. + .001)) Opt_Dist_F(&(d_max),F,mod);
          else
            {
              PhyML_Printf("\n== sum = %f\n",sum);
              PhyML_Printf("\n== Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
              Exit("");
            }
          
          if(d_max >= DIST_MAX) d_max = DIST_MAX;
          
          
          /* Do not correct for dist < BL_MIN, otherwise Fill_Missing_Dist
           *  will not be called
           */
          mat->dist[j][k] = d_max;
          mat->dist[k][j] = mat->dist[j][k];
          Free_Cseq(twodata);
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


  For(i,mod->ras->n_catg)
    {
      len = dist*mod->ras->gamma_rr->v[i];
      if(len < mod->l_min) len = mod->l_min;
      else if(len > mod->l_max) len = mod->l_max;
      PMat(len,mod,dim2*i,mod->Pij_rr->v);
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
      PhyML_Printf("\n. Not implemented yet...");
      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
      Warn_And_Exit("\n");
    }


  site_lk = .0;
  *loglk = 0;

  For(i,data->c_seq[0]->len)
    {
      if(data->wght[i])
    {
      site_lk = log_site_lk = .0;
      if(!data->ambigu[i])
        {
          For(k,mod->ns) {if(p_lk_l[i*mod->ns+k] > .0001) break;}
          For(l,mod->ns) {if(p_lk_r[i*mod->ns+l] > .0001) break;}
          For(j,mod->ras->n_catg)
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
          For(j,mod->ras->n_catg)
        {
          For(k,mod->ns) /*sort sum terms ? No global effect*/
            {
              For(l,mod->ns)
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
          PhyML_Printf("'%c' '%c'\n",seq1->state[i],seq2->state[i]);
          Exit("\n. Err: site lk <= 0\n");
        }

      log_site_lk += (phydbl)LOG(site_lk);

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

  For(i,tree->data->crunch_len)
    {
      tree->unconstraint_lk +=
    tree->data->wght[i]*(phydbl)LOG(tree->data->wght[i]);
    }
  tree->unconstraint_lk -=
    tree->data->init_len*(phydbl)LOG(tree->data->init_len);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Init_P_Lk_Tips_Double(t_tree *tree)
{
  int curr_site,i,j,k,dim1,dim2;

  dim1 = tree->mod->ras->n_catg * tree->mod->ns;
  dim2 = tree->mod->ns;


  For(i,tree->n_otu)
    {
      if(!tree->a_nodes[i]->c_seq || 
	 strcmp(tree->a_nodes[i]->c_seq->name,tree->a_nodes[i]->name))
        {
          PhyML_Printf("\n== Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
          Exit("");
        }
    }

  For(curr_site,tree->data->crunch_len)
    {
      For(i,tree->n_otu)
        {
          if (tree->io->datatype == NT)
        /* Init_Tips_At_One_Site_Nucleotides_Float(tree->data->c_seq[i]->state[curr_site], */
        /* 					    curr_site*dim1+0*dim2, */
        /* 					    tree->a_nodes[i]->b[0]->p_lk_rght); */
            Init_Tips_At_One_Site_Nucleotides_Float(tree->a_nodes[i]->c_seq->state[curr_site],
                                                    curr_site*dim1+0*dim2,
                                                    tree->a_nodes[i]->b[0]->p_lk_rght);

          else if(tree->io->datatype == AA)
            Init_Tips_At_One_Site_AA_Float(tree->a_nodes[i]->c_seq->state[curr_site],
                                           curr_site*dim1+0*dim2,
                                           tree->a_nodes[i]->b[0]->p_lk_rght);
        /* Init_Tips_At_One_Site_AA_Float(tree->data->c_seq[i]->state[curr_site], */
        /* 				   curr_site*dim1+0*dim2, */
        /* 				   tree->a_nodes[i]->b[0]->p_lk_rght); */

          else if(tree->io->datatype == GENERIC)
            Init_Tips_At_One_Site_Generic_Float(tree->a_nodes[i]->c_seq->state+curr_site*tree->mod->io->state_len,
                                                tree->mod->ns,
                                                tree->mod->io->state_len,
                                                curr_site*dim1+0*dim2,
                                                tree->a_nodes[i]->b[0]->p_lk_rght);
        /* Init_Tips_At_One_Site_Generic_Float(tree->data->c_seq[i]->state+curr_site*tree->mod->io->state_len, */
        /* 					tree->mod->ns, */
        /* 					tree->mod->io->state_len, */
        /* 					curr_site*dim1+0*dim2, */
        /* 					tree->a_nodes[i]->b[0]->p_lk_rght); */

          for(j=1;j<tree->mod->ras->n_catg;j++)
            {
              For(k,tree->mod->ns)
                {
                  tree->a_nodes[i]->b[0]->p_lk_rght[curr_site*dim1+j*dim2+k] =
                    tree->a_nodes[i]->b[0]->p_lk_rght[curr_site*dim1+0*dim2+k];
                }
            }
        }
    }

  if(tree->mod->use_m4mod)
    {
      M4_Init_P_Lk_Tips_Double(tree);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Init_P_Lk_Tips_Int(t_tree *tree)
{
  int curr_site,i,dim1;

  dim1 = tree->mod->ns;

  For(i,tree->n_otu)
    {
      if(!tree->a_nodes[i]->c_seq || 
	 strcmp(tree->a_nodes[i]->c_seq->name,tree->a_nodes[i]->name))
        {
          PhyML_Printf("\n== Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
          Exit("");
        }
    }

  For(curr_site,tree->data->crunch_len)
    {
      For(i,tree->n_otu)
        {
          /* printf("\n. site: %3d %c",curr_site,tree->a_nodes[i]->c_seq->state[curr_site]); */
          /* printf("\n. init at %s %p",tree->a_nodes[i]->name,tree->a_nodes[i]->b[0]->p_lk_tip_r); fflush(NULL); */
          if(tree->io->datatype == NT)
            {
              Init_Tips_At_One_Site_Nucleotides_Int(tree->a_nodes[i]->c_seq->state[curr_site],
                                                    curr_site*dim1,
                                                    tree->a_nodes[i]->b[0]->p_lk_tip_r);
          /* Init_Tips_At_One_Site_Nucleotides_Int(tree->data->c_seq[i]->state[curr_site], */
          /* 					    curr_site*dim1, */
          /* 					    tree->a_nodes[i]->b[0]->p_lk_tip_r); */
            }
          else if(tree->io->datatype == AA)
            Init_Tips_At_One_Site_AA_Int(tree->a_nodes[i]->c_seq->state[curr_site],
                         curr_site*dim1,
                         tree->a_nodes[i]->b[0]->p_lk_tip_r);
        /* Init_Tips_At_One_Site_AA_Int(tree->data->c_seq[i]->state[curr_site], */
        /* 				 curr_site*dim1,					    */
        /* 				 tree->a_nodes[i]->b[0]->p_lk_tip_r); */

          else if(tree->io->datatype == GENERIC)
            {
              Init_Tips_At_One_Site_Generic_Int(tree->a_nodes[i]->c_seq->state+curr_site*tree->mod->io->state_len,
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
       M4_Init_P_Lk_Tips_Int(tree);
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

  if(tree->is_mixt_tree)
    {
      MIXT_Update_PMat_At_Given_Edge(b_fcus,tree);
      return;
    }

  l_min = tree->mod->l_min;
  l_max = tree->mod->l_max;

  len = -1.0;

  if(tree->mod->log_l == YES) {b_fcus->l->v = EXP(b_fcus->l->v);}

  For(i,tree->mod->ras->n_catg)
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
          if(tree->mixt_tree)  var *= POW(tree->mixt_tree->mod->ras->gamma_rr->v[tree->mod->ras->parent_class_number],2);

          if(var > tree->mod->l_var_max) var = tree->mod->l_var_max;
          if(var < tree->mod->l_var_min) var = tree->mod->l_var_min;
        }

      //Update the transition prob. matrix
      if(tree->mod->gamma_mgf_bl == NO)
          {
#ifdef BEAGLE
            assert(UNINITIALIZED != tree->mod->b_inst);
#endif
            PMat(len,tree->mod,i*tree->mod->ns*tree->mod->ns,b_fcus->Pij_rr);
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
      int p_matrices[1]     = {b_fcus->Pij_rr_idx};
      double branch_lens[1] = {len};
      int ret = beagleUpdateTransitionMatrices(tree->b_inst,0,p_matrices,NULL,NULL,branch_lens,1);
      if(ret<0){
          fprintf(stderr, "beagleUpdateTransitionMatrices() on instance %i failed:%i\n\n",tree->b_inst,ret);
          Exit("");
      }
      //Retrieve a "local" copy of the P-matrix
      ret = beagleGetTransitionMatrix(tree->b_inst, b_fcus->Pij_rr_idx, b_fcus->Pij_rr);
      if(ret<0){
          fprintf(stderr, "beagleGetTransitionMatrix() on instance %i failed:%i\n\n",tree->b_inst,ret);
          Exit("");
      }
  }
  else
  {
      int ret = beagleSetTransitionMatrix(tree->b_inst, b_fcus->Pij_rr_idx, b_fcus->Pij_rr, -1);
      if(ret<0){
          fprintf(stderr, "beagleSetTransitionMatrix() on instance %i failed:%i\n\n",tree->b_inst,ret);
          Exit("");
      }
  }
#endif
    if(tree->mod->log_l == YES) {b_fcus->l->v = LOG(b_fcus->l->v);}

//      Print_Model(tree->mod);
//      Dump_Arr_D(tree->cur_site_lk, tree->n_pattern);
//      Print_Edge_PMats(tree, b_fcus);
//      Print_Edge_Likelihoods(tree,b_fcus,true);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


/* void Update_P_Lk_On_A_Path(t_node *a, t_node *d, t_edge *b, t_node *target_one_side, t_node *target_other_side, t_tree *tree) */
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

/*   Update_P_Lk(tree,b,a); */
/*   if((a == target_one_side) && (d == target_other_side))  */
/*     return; */
/*   else */
/*     { */
/*       Update_P_Lk_On_A_Path(d, */
/* 			    d->v[tree->t_dir[d->num][target_other_side->num]], */
/* 			    d->b[tree->t_dir[d->num][target_other_side->num]], */
/* 			    target_one_side, */
/* 			    target_other_side, */
/* 			    tree);  */
/*     } */
/* } */

void Update_P_Lk_Along_A_Path(t_node **path, int path_length, t_tree *tree)
{
  int i,j;

  For(i,path_length-1)
    {
      For(j,3)
    if(path[i]->v[j] == path[i+1])
      {
        if(path[i] == path[i]->b[j]->left)
          {
        Update_P_Lk(tree,path[i]->b[j],path[i]->b[j]->left);
          }
        else if(path[i] == path[i]->b[j]->rght)
          {
        Update_P_Lk(tree,path[i]->b[j],path[i]->b[j]->rght);
          }
        else
          {
        PhyML_Printf("\n== Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
        Exit("");
          }
        break;
      }
#ifdef DEBUG
      if(j == 3)
    {
      PhyML_Printf("\n== Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
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

  if(mod->log_l == YES) dist = EXP(dist);

  For(k,mod->ras->n_catg)
    {
      len = dist*mod->ras->gamma_rr->v[k];
      if(len < mod->l_min)      len = mod->l_min;
      else if(len > mod->l_max) len = mod->l_max;
      PMat(len,mod,mod->ns*mod->ns*k,mod->Pij_rr->v);
    }

  dim1 = mod->ns*mod->ns;
  dim2 = mod->ns;
  lnL = .0;

/*   For(i,mod->ns) */
/*     { */
/*       For(j,mod->ns) */
/* 	{ */
/* 	  For(k,mod->ras->n_catg) */
/* 	    { */
/*  	      lnL += */
/* 		F[dim1*k+dim2*i+j] * */
/* 		LOG(mod->pi[i] * mod->Pij_rr[dim1*k+dim2*i+j]); */
/* 	    } */
/* 	} */
/*     } */

  For(i,mod->ns-1)
    {
      for(j=i+1;j<mod->ns;j++)
        {
          For(k,mod->ras->n_catg)
            {
              lnL +=
                (F[dim1*k+dim2*i+j] + F[dim1*k+dim2*j+i])*
                LOG(mod->e_frq->pi->v[i] * mod->Pij_rr->v[dim1*k+dim2*i+j]);
              /* printf("\n. f: %f Pij:%f F:%f", */
              /*        mod->e_frq->pi->v[i],  */
              /*        mod->Pij_rr->v[dim1*k+dim2*i+j], */
              /*        F[dim1*k+dim2*j+i]); */
            }
        }
    }
  
  For(i,mod->ns) For(k,mod->ras->n_catg) lnL += F[dim1*k+dim2*i+i]* LOG(mod->e_frq->pi->v[i] * mod->Pij_rr->v[dim1*k+dim2*i+i]);
  
  return lnL;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Update_Lk_At_Given_Edge(t_edge *b_fcus, t_tree *tree)
{
  Update_P_Lk(tree,b_fcus,b_fcus->left);
  Update_P_Lk(tree,b_fcus,b_fcus->rght);
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
      For(i,3)
    if(d->v[i] != a)
      Print_Lk_Given_Edge_Recurr(d,d->v[i],d->b[i],tree);
    }
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

      For(i,tree->n_pattern)
    {
      For(j,tree->n_pattern)
        {
          if(patt_id_d[i] == patt_id_d[j])
        {
          p_lk_loc_d[i] = j;
          break;
        }
          if(j > i)
        {
          PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
          Warn_And_Exit("");
        }
        }
    }
      return;
    }
  else
    {
      v1 = v2 = NULL;
      For(i,3)
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
      For(i,tree->n_pattern)
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
            PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
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

      For(i,3)
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

      For(i,3)
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

void Copy_P_Lk(phydbl *p_lk, int site_from, int site_to, t_tree *tree)
{
  int i,j;
  int dim1,dim2;


  dim1 = tree->mod->ras->n_catg * tree->mod->ns;
  dim2 = tree->mod->ns;

/*   PhyML_Printf("\n# %d %d",site_to,site_from); */

  For(i,tree->mod->ns) For(j,tree->mod->ras->n_catg)
    {
      p_lk[site_to*dim1+j*dim2+i] = p_lk[site_from*dim1+j*dim2+i];
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Copy_Scale(int *scale, int site_from, int site_to, t_tree *tree)
{
  int i;

  For(i,tree->mod->ras->n_catg)
    {
      scale[i*tree->n_pattern+site_to] = scale[i*tree->n_pattern+site_from];
/*       PhyML_Printf("\n. %d",scale[i*tree->n_pattern+site_to]); */
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Init_P_Lk_Loc(t_tree *tree)
{
  int i,j;
  t_node *d;
  int *patt_id_d;

  For(i,2*tree->n_otu-1)
    {
      For(j,tree->n_pattern)
        {
          tree->a_edges[i]->p_lk_loc_left[j] = j;
          tree->a_edges[i]->p_lk_loc_rght[j] = j;
        }
    }
  
  For(i,tree->n_otu)
    {
      d = tree->a_nodes[i];
      patt_id_d = (d == d->b[0]->left)?(d->b[0]->patt_id_left):(d->b[0]->patt_id_rght);
      For(j,tree->n_pattern)
    {
      patt_id_d[j] = (int)tree->a_nodes[d->num]->c_seq->state[j];
      /* patt_id_d[j] = (int)tree->data->c_seq[d->num]->state[j]; */
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
  For(i,dim) first_order += (tree->rates->u_cur_l[i] - tree->rates->mean_l[i]) * tree->rates->grad_l[i];

  lnL += first_order;

  /* printf("\n"); */
  /* For(i,dim) printf("%f\t",tree->rates->u_cur_l[i]); */
  /* printf("\n. Lk=%f %f",lnL,tree->mod->l_min); */


/*   err = NO; */
/*   dim = 2*tree->n_otu-3; */
/*   For(i,dim) */
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
/* 	  l = (tree->mod->log_l == YES)?(EXP(tree->a_edges[i]->l->v)):(tree->a_edges[i]->l->v); */
/* 	  lnL += LOG(Dexp(l,lambda)); */
/* /\* 	  printf("\n. lambda = %f",lambda); *\/ */
/* 	} */
/*     } */

  return(lnL);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Wrap_Part_Lk_At_Given_Edge(t_edge *b, t_tree *tree, supert_tree *stree)
{
  return PART_Lk_At_Given_Edge(b,stree);;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Wrap_Part_Lk(t_edge *b, t_tree *tree, supert_tree *stree)
{
  return PART_Lk(stree);
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
  TIMES_Lk_Times(tree);
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

  For(i,tree->n_pattern)
    {
      /* Sample the rate class from its posterior density */
      For(j,tree->mod->ras->n_catg)
        {
          if(fromprior == NO)
            probs[j] =
              tree->unscaled_site_lk_cat[j*tree->n_pattern+i]*
              tree->mod->ras->gamma_r_proba->v[j];
          else
            probs[j] = tree->mod->ras->gamma_r_proba->v[j];
        }


      /* Scale probas. */
      sum = .0;
      For(j,tree->mod->ras->n_catg) sum += probs[j];
      For(j,tree->mod->ras->n_catg) probs[j]/=sum;

      /* CDF */
      for(j=1;j<tree->mod->ras->n_catg;j++) probs[j] += probs[j-1];

      /* Sample rate */
      u = Uni();
      rate_cat = -1;
      For(j,tree->mod->ras->n_catg)
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
      

      For(j,n_mut) ordering[j] = 0;
      
      For(j,n_mut-1)
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

      For(j,n_mut)
    {
      For(k,n_mut)
        {
          if(ordering[k] == j)
        {
          For(l,tree->data->init_len) if(tree->data->sitepatt[l] == i) break;
          PhyML_Fprintf(fp,"\n%4d %4d %12f %4d",j,muttype[k],muttime[k],l);
          /* PhyML_Fprintf(fp,"\n%d",muttype[ordering[j]]); */
          break;
        }
        }
    }


      For(j,n_mut)
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
      For(j,tree->mod->ns) probs[j] = tree->mod->e_frq->pi->v[j];
      
      for(j=1;j<tree->mod->ns;j++) probs[j] += probs[j-1];
      
      u = Uni();
      For(j,tree->mod->ns)
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
      
      For(j,tree->mod->ns) probs[j] = 0.0;
      
      /* Formula (10) in Nielsen's Mutation Maping paper, e.g. */
      For(j,tree->mod->ns)
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
      For(j,tree->mod->ns) sum += probs[j];
      For(j,tree->mod->ns) probs[j] /= sum;
      
      /* CDF */
      for(j=1;j<tree->mod->ns;j++) probs[j] += probs[j-1];
      
      /* Sample state according to their posterior probas. */
      sd = -1;
      u = Uni();
      For(j,tree->mod->ns)
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
          For(j,tree->mod->ns) probs[j] = tree->mod->e_frq->pi->v[j];
          
          for(j=1;j<tree->mod->ns;j++) probs[j] += probs[j-1];
          
          u = Uni();
          For(j,tree->mod->ns)
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
      For(i,3)
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
  br =
    (tree->rates->model_log_rates == YES)?
    EXP(tree->rates->br_r[d->num]):
    tree->rates->br_r[d->num];

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
  For(i,tree->mod->ns)
    {
      // We only care about the non-diagonal elements here
      For(j,tree->mod->ns) all_probs[i*tree->mod->ns+j] = -Q[i*tree->mod->ns+j] / Q[i*tree->mod->ns+i];

      // Set the diagonal to 0
      all_probs[i*tree->mod->ns+i] = 0.0;

      // \sum_{j != i} -q_{ij}/q_{ii}
      sum = 0;
      For(j,tree->mod->ns) sum += all_probs[i*tree->mod->ns+j];

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
      tlast = -LOG(1. - u*(1.-EXP(Q[sa*tree->mod->ns+sa]*td)))/-Q[sa*tree->mod->ns+sa];
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
      For(i,tree->mod->ns)
        if(probs[i] > u)
          {
        // Record mutation type
        mut[slast*tree->mod->ns+i]++;

        // Record mutation type in the site mutation array
        thismut = MIN(i,slast) * tree->mod->ns + MAX(i,slast) - (MIN(i,slast)+1+(int)POW(MIN(i,slast)+1,2))/2;
        muttype[(*n_mut)+n_mut_branch-1] = thismut;

        if((thismut > 5) || (thismut < 0))
          {
            PhyML_Printf("\n. thismut = %d",thismut);
            PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
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
          For(i,n_mut_branch) muttype[(*n_mut)+n_mut_branch-1] = -2;
          For(i,n_mut_branch) muttime[(*n_mut)+n_mut_branch-1] = +1.;
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


  For(i,tree->mod->ns)
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

void ML_Ancestral_Sequences(t_tree *tree)
{
  int i;

  /* PhyML_Printf("\n. %s \n",Write_Tree(tree,NO)); */

  PhyML_Printf("\n\n. Estimating ancestral sequences...");

  strcpy(tree->io->out_ancestral_file,tree->io->out_file);
  if(tree->io->append_run_ID) { strcat(tree->io->out_ancestral_file,"_"); strcat(tree->io->out_ancestral_file,tree->io->run_id_string); }
  strcat(tree->io->out_ancestral_file,"_phyml_ancestral_seq");
  tree->io->fp_out_ancestral = Openfile(tree->io->out_ancestral_file,1);

  if(tree->n_root)
    {
      PhyML_Fprintf(tree->io->fp_out_ancestral,"\n== Printing the tree structure. Starting from the root node");
      PhyML_Fprintf(tree->io->fp_out_ancestral,"\n== and displaying the nodes underneath recursively. Edge numbers");
      PhyML_Fprintf(tree->io->fp_out_ancestral,"\n== and their lengths are also provided.\n\n");
      Print_Node_Brief(tree->n_root,tree->n_root->v[0],tree,tree->io->fp_out_ancestral);
      Print_Node_Brief(tree->n_root,tree->n_root->v[1],tree,tree->io->fp_out_ancestral);
    }
  else
    {
      PhyML_Fprintf(tree->io->fp_out_ancestral,"\n== Printing the tree structure. Starting from node 0 (taxon %s)",tree->a_nodes[0]->name);
      PhyML_Fprintf(tree->io->fp_out_ancestral,"\n== and displaying the nodes underneath recursively. Edge numbers");
      PhyML_Fprintf(tree->io->fp_out_ancestral,"\n== and their lengths are also provided.\n\n");
      Print_Node_Brief(tree->a_nodes[0],tree->a_nodes[0]->v[0],tree,tree->io->fp_out_ancestral);
    }
  
  PhyML_Fprintf(tree->io->fp_out_ancestral,"\n\n\n");
  PhyML_Fprintf(tree->io->fp_out_ancestral,"\n== Printing marginal probabilities of ancestral sequences at each site");
  PhyML_Fprintf(tree->io->fp_out_ancestral,"\n== of the alignment and each node of the tree.");
  PhyML_Fprintf(tree->io->fp_out_ancestral,"\n\n");
  PhyML_Fprintf(tree->io->fp_out_ancestral,"Site\tNode\t");
  For(i,tree->mod->ns) PhyML_Fprintf(tree->io->fp_out_ancestral,"%c\t",Reciproc_Assign_State(i,tree->io->datatype));
  PhyML_Fprintf(tree->io->fp_out_ancestral,"\n");

  For(i,2*tree->n_otu-2)
    if(tree->a_nodes[i]->tax == NO)
      ML_Ancestral_Sequences_One_Node(tree->a_nodes[i],tree);

  if(tree->n_root) ML_Ancestral_Sequences_One_Node(tree->n_root,tree);


  fclose(tree->io->fp_out_ancestral);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void ML_Ancestral_Sequences_One_Node(t_node *mixt_d, t_tree *mixt_tree)
{
  if(mixt_d->tax) return;
  else
    {
      t_node *v0,*v1,*v2; // three neighbours of d
      t_edge *b0,*b1,*b2;
      int i,j;
      int catg;
      phydbl p0, p1, p2;
      phydbl *p;
      t_node *d,*curr_mixt_d;
      t_tree *tree, *curr_mixt_tree;
      int site,csite;
      phydbl *p_lk0, *p_lk1, *p_lk2;
      int *sum_scale0, *sum_scale1, *sum_scale2;
      phydbl r_mat_weight_sum, e_frq_weight_sum, sum_probas;
      phydbl *Pij0, *Pij1, *Pij2;
      int NsNs, Ns, NsNg;
      FILE *fp;

      if(!mixt_d) return;

      curr_mixt_tree = mixt_tree;
      curr_mixt_d    = mixt_d;
      fp             = mixt_tree->io->fp_out_ancestral;


      do /* For each partition element */
        {
          if(curr_mixt_tree->next)
            {
              r_mat_weight_sum = MIXT_Get_Sum_Chained_Scalar_Dbl(curr_mixt_tree->next->mod->r_mat_weight);
              e_frq_weight_sum = MIXT_Get_Sum_Chained_Scalar_Dbl(curr_mixt_tree->next->mod->e_frq_weight);
              sum_probas       = MIXT_Get_Sum_Of_Probas_Across_Mixtures(r_mat_weight_sum, e_frq_weight_sum, curr_mixt_tree);
            }
          else
            {
              r_mat_weight_sum = 1.;
              e_frq_weight_sum = 1.;
              sum_probas       = 1.;
            }

          Ns   = curr_mixt_tree->next ? curr_mixt_tree->next->mod->ns : curr_mixt_tree->mod->ns;
          NsNs = Ns*Ns;
          NsNg = Ns*curr_mixt_tree->mod->ras->n_catg;

          p = (phydbl *)mCalloc(Ns,sizeof(phydbl));

          /* For(site,curr_mixt_tree->n_pattern) // For each site in the current partition element */
          For(site,curr_mixt_tree->data->init_len) // For each site in the current partition element
            {
              csite = curr_mixt_tree->data->sitepatt[site];

              d    = curr_mixt_d->next ? curr_mixt_d->next : curr_mixt_d;
              tree = curr_mixt_tree->next ? curr_mixt_tree->next : curr_mixt_tree;

              For(i,tree->mod->ns) p[i] = .0;

              do // For each class of the mixture model that applies to the current partition element
                {
                  if(tree->is_mixt_tree == YES)
                    {
                      tree = tree->next;
                      d    = d->next;
                    }

                  v0 = d->v[0];
                  v1 = d->v[1];
                  v2 = d->v[2];

                  b0 = d->b[0];
                  b1 = d->b[1];
                  b2 = d->b[2];

                  Pij0 = b0->Pij_rr;
                  Pij1 = b1->Pij_rr;
                  Pij2 = b2->Pij_rr;

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


                  For(catg,tree->mod->ras->n_catg)
                    {
                      For(i,Ns)
                        {
                          p0 = .0;
                          if(v0->tax)
                            For(j,tree->mod->ns)
                              {
                                p0 += v0->b[0]->p_lk_tip_r[csite*Ns+j] * Pij0[catg*NsNs+i*Ns+j];

                                /* printf("\n. p0 %d %f", */
                                /*        v0->b[0]->p_lk_tip_r[site*Ns+j], */
                                /*        Pij0[catg*NsNs+i*Ns+j]); */
                              }
                          else
                            For(j,tree->mod->ns)
                              {
                                p0 += p_lk0[csite*NsNg+catg*Ns+j] * Pij0[catg*NsNs+i*Ns+j] / (phydbl)POW(2,sum_scale0[catg*curr_mixt_tree->n_pattern+csite]);

                                /* p0 += p_lk0[site*NsNg+catg*Ns+j] * Pij0[catg*NsNs+i*Ns+j]; */

                                /* printf("\n. p0 %f %f", */
                                /*        p_lk0[site*NsNg+catg*Ns+j], */
                                /*        Pij0[catg*NsNs+i*Ns+j]); */
                              }
                          p1 = .0;
                          if(v1->tax)
                            For(j,tree->mod->ns)
                              {
                                p1 += v1->b[0]->p_lk_tip_r[csite*Ns+j] * Pij1[catg*NsNs+i*Ns+j];

                                /* printf("\n. p1 %d %f", */
                                /*        v1->b[0]->p_lk_tip_r[site*Ns+j], */
                                /*        Pij1[catg*NsNs+i*Ns+j]); */
                              }

                          else
                            For(j,tree->mod->ns)
                              {
                                p1 += p_lk1[csite*NsNg+catg*Ns+j] * Pij1[catg*NsNs+i*Ns+j] / (phydbl)POW(2,sum_scale1[catg*curr_mixt_tree->n_pattern+csite]);

                                /* p1 += p_lk1[site*NsNg+catg*Ns+j] * Pij1[catg*NsNs+i*Ns+j];  */

                                /* printf("\n. p1 %f %f", */
                                /*        p_lk1[site*NsNg+catg*Ns+j], */
                                /*        Pij1[catg*NsNs+i*Ns+j]); */
                             }


                          p2 = .0;
                          if(v2->tax)
                            For(j,tree->mod->ns)
                              {
                                p2 += v2->b[0]->p_lk_tip_r[csite*Ns+j] * Pij2[catg*NsNs+i*Ns+j];
                                /* printf("\n. p2 %d %f", */
                                /*        v2->b[0]->p_lk_tip_r[site*Ns+j], */
                                /*        Pij2[catg*NsNs+i*Ns+j]); */
                              }
                          else
                            For(j,tree->mod->ns)
                              {
                                p2 += p_lk2[csite*NsNg+catg*Ns+j] * Pij2[catg*NsNs+i*Ns+j] / (phydbl)POW(2,sum_scale2[catg*curr_mixt_tree->n_pattern+csite]);

                                /* p2 += p_lk2[site*NsNg+catg*Ns+j] * Pij2[catg*NsNs+i*Ns+j]; */

                                /* printf("\n. p2 %f %f", */
                                /*        p_lk2[site*NsNg+catg*Ns+j], */
                                /*        Pij2[catg*NsNs+i*Ns+j]);  */
                             }

                          p[i] +=
                            p0*p1*p2*
                            tree->mod->e_frq->pi->v[i] /
                            tree->cur_site_lk[csite] *
                            curr_mixt_tree->mod->ras->gamma_r_proba->v[tree->mod->ras->parent_class_number] *
                            tree->mod->r_mat_weight->v / r_mat_weight_sum *
                            tree->mod->e_frq_weight->v / e_frq_weight_sum /
                            sum_probas;

                        }
                    }

                  PhyML_Fprintf(fp,"%4d\t%4d\t",site+1,d->num);
                  For(i,Ns)
                    {
                      PhyML_Fprintf(fp,"%.4f\t",p[i]);
                    }
                  PhyML_Fprintf(fp,"\n");
                  fflush(NULL);
                  /* Exit("\n"); */
                  

                  tree = tree->next;
                  d    = d->next;


                }
              while(tree && d && tree->is_mixt_tree == NO);
            }

          Free(p);
          curr_mixt_tree = curr_mixt_tree->next_mixt;
          curr_mixt_d    = curr_mixt_d->next_mixt;
        }
      while(curr_mixt_tree != NULL);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// Returns the value of fact_sum_scale, the part of the scaling factors
// that is common to all classes of the mixture and scale the 
// likelihood for each mixture (using the part of the scaling
// factors that is class-specific)

void Pull_Scaling_Factors(int site,
                          t_edge *b,
                          t_tree *tree)
{
  int catg;
  int *sum_scale_left_cat,*sum_scale_rght_cat;
  int exponent;
  phydbl max_sum_scale,min_sum_scale;
  phydbl sum,tmp;
  phydbl site_lk_cat;

  sum_scale_left_cat = b->sum_scale_left_cat;
  sum_scale_rght_cat = b->sum_scale_rght_cat;

  if(tree->apply_lk_scaling == YES)
    {
      max_sum_scale =   (phydbl)BIG;
      min_sum_scale =  -(phydbl)BIG;
      
      For(catg,tree->mod->ras->n_catg)
        {
          sum_scale_left_cat[catg] =
            (b->sum_scale_left)?
            (b->sum_scale_left[catg*tree->n_pattern+site]):
            (0.0);
          
          sum_scale_rght_cat[catg] =
            (b->sum_scale_rght)?
            (b->sum_scale_rght[catg*tree->n_pattern+site]):
            (0.0);
          
          sum = sum_scale_left_cat[catg] + sum_scale_rght_cat[catg];
          
          if(sum < .0)
            {
              printf("\n== tree: %s\n",Write_Tree(tree,NO));
              PhyML_Printf("\n== b->num = %d  sum = %G root ? %d",sum,b->num,b == tree->e_root);
              PhyML_Printf("\n== Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
              Exit("\n");
            }
          
          tmp = sum + ((phydbl)LOGBIG - LOG(tree->site_lk_cat[catg]))/(phydbl)LOG2;
          if(tmp < max_sum_scale) max_sum_scale = tmp; /* min of the maxs */
          
          tmp = sum + ((phydbl)LOGSMALL - LOG(tree->site_lk_cat[catg]))/(phydbl)LOG2;
          if(tmp > min_sum_scale) min_sum_scale = tmp; /* max of the mins */
          
        }
      
      if(min_sum_scale > max_sum_scale)
        {
          /* PhyML_Printf("\n== Numerical precision issue alert."); */
          /* PhyML_Printf("\n== min_sum_scale = %G max_sum_scale = %G",min_sum_scale,max_sum_scale); */
          /* PhyML_Printf("\n== Err in file %s at line %d\n\n",__FILE__,__LINE__); */
          /* Warn_And_Exit("\n"); */
          min_sum_scale = max_sum_scale;
        }
      
      tree->fact_sum_scale[site] = (int)((max_sum_scale + min_sum_scale) / 2);
            
      /* fact_sum_scale = (int)(max_sum_scale / 2); */
      
      /* Apply scaling factors */
      For(catg,tree->mod->ras->n_catg)
        {
          exponent = -(sum_scale_left_cat[catg]+sum_scale_rght_cat[catg])+tree->fact_sum_scale[site];
          site_lk_cat = tree->site_lk_cat[catg];
          Rate_Correction(exponent,&site_lk_cat,tree);
          tree->site_lk_cat[catg] = site_lk_cat;
        }
    }
  else // No scaling of lk
    {
      tree->fact_sum_scale[site] = 0;
    }
  
  For(catg,tree->mod->ras->n_catg) 
    {
      tree->unscaled_site_lk_cat[catg*tree->n_pattern + site] = tree->site_lk_cat[catg];
      
      if(isinf(tree->unscaled_site_lk_cat[catg*tree->n_pattern + site]) || 
         isnan(tree->unscaled_site_lk_cat[catg*tree->n_pattern + site]))
        {
          PhyML_Printf("\n== Err. in file %s at line %d (function '%s')",__FILE__,__LINE__,__FUNCTION__);
          Exit("\n");
        }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

