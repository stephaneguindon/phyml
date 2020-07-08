/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include "ancestral.h"

void Sample_Ancestral_Seq(int fullmutmap, int fromprior, t_tree *tree)
{
  int rate_cat;
  int i,j,k,l;
  phydbl *probs;
  phydbl sum;
  int n_mut;
  FILE *fp;
  phydbl *muttime;
  int *muttype,*muttax;
  char *s,*mut_string;
  int *ordering;

  
  if(tree->is_mixt_tree == YES)
    {
      MIXT_Sample_Ancestral_Seq(fullmutmap,fromprior,tree);
      return;
    }

  probs = (phydbl *)mCalloc(tree->mod->ras->n_catg,sizeof(phydbl));

  muttime = (phydbl *)mCalloc((2*tree->n_otu-3)*5, // At most 5 mutations per branch on average
               sizeof(phydbl));

  muttype = (int *)mCalloc((2*tree->n_otu-3)*5, // At most 5 mutations per branch on average
               sizeof(int));

  muttax = (int *)mCalloc((2*tree->n_otu-3)*5, // At most 5 mutations per branch on average
                          sizeof(int));

  ordering = (int *)mCalloc((2*tree->n_otu-3)*5, // At most 5 mutations per branch on average
                sizeof(int));

  s = (char *)mCalloc(T_MAX_NAME,sizeof(char));

  for(i=0;i<2*tree->n_otu-1;++i)
    if(tree->a_nodes[i]->tax == NO)
      {
        tree->a_nodes[i]->c_seq_anc = (align *)mCalloc(1,sizeof(align));;
        tree->a_nodes[i]->c_seq_anc->state = (char *)mCalloc(tree->n_pattern,sizeof(char));
      }

   
  /* Update P(D_x|X=i) for each state i and node X */
  Set_Both_Sides(YES,tree);
  Lk(NULL,tree);


  for(i=0;i<tree->n_pattern;++i)
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

      rate_cat = Sample_i_With_Proba_pi(probs,tree->mod->ras->n_catg);
           
      n_mut = 0;
      if(tree->n_root != NULL)
        {          
          Sample_Ancestral_Seq_Pre(tree->n_root,tree->n_root->v[1],tree->n_root->b[1],
                                   i,rate_cat,
                                   muttype,muttime,muttax,&n_mut,
                                   fullmutmap,fromprior,tree);

          
          Sample_Ancestral_Seq_Pre(tree->n_root,tree->n_root->v[2],tree->n_root->b[2],
                                   i,rate_cat,
                                   muttype,muttime,muttax,&n_mut,
                                   fullmutmap,fromprior,tree);


        }
      else
        {
          Sample_Ancestral_Seq_Pre(tree->a_nodes[0],tree->a_nodes[0]->v[0],tree->a_nodes[0]->b[0],
                                   i,rate_cat,
                                   muttype,muttime,muttax,&n_mut,
                                   fullmutmap,fromprior,tree);
        }
      
      for(j=0;j<n_mut;j++) ordering[j] = 0;
      
      for(j=0;j<n_mut-1;j++)
        {
          for(k=j+1;k<n_mut;k++)
            {
              if(muttime[k] > muttime[j]) ordering[k]++;
              else ordering[j]++;
            }
        }
      
      strcpy(s,"mutmap.");
      sprintf(s+strlen(s),"%d.",tree->io->r_seed);
      sprintf(s+strlen(s),"%d",i);
      fp = fopen(s,"a");
      PhyML_Fprintf(fp,"\n-1 -1 -1.0 -1 %d",tree->mixt_tree ? tree->mixt_tree->mcmc->run : tree->mcmc->run);
      
      Reverse_Muttime(muttime,n_mut,tree);

      for(j=0;j<n_mut;j++)
        {
          for(k=0;k<n_mut;k++)
            {
              if(ordering[k] == j)
                {
                  for(l=0;l<tree->data->init_len;l++) if(tree->data->sitepatt[l] == i) break;
                  mut_string = Mutation_Id(muttype[k],tree);
                  PhyML_Fprintf(fp,"\n%4d %s %g %4d %s",j,mut_string,muttime[k],l,tree->a_nodes[muttax[k]]->name);
                  Free(mut_string);
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
  Free(muttax);
  Free(ordering);
  Free(probs);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Sample_Ancestral_Seq_Pre(t_node *a, t_node *d, t_edge *b,
                              int site, int r_cat,
                              int *muttype, phydbl *muttime, int *muttax, int *n_mut,
                              int fullmutmap, int fromprior, t_tree *tree)
{

  int i;
  int sa, sd;
  int ns;
  phydbl *probs,sum;

  /* PhyML_Printf("\n>> a: %d d: %d b->left: %d b->rght: %d",a?a->num:-1,d?d->num:-1,b?b->left->num:-1,b?b->rght->num:-1); */

  ns = tree->mod->ns;
  
  probs = (phydbl *)mCalloc(ns,sizeof(phydbl));

  if(a->tax == TRUE) // Sample state at tip if observed state is ambiguous
    {
      assert(b);
      if(a == b->left) for(i=0;i<ns;++i) probs[i] = b->p_lk_tip_l[site*ns+i];
      else             for(i=0;i<ns;++i) probs[i] = b->p_lk_tip_r[site*ns+i];
      for(i=0;i<ns;++i) probs[i] *= tree->mod->e_frq->pi->v[i];
      sum = 0.0;
      for(i=0;i<ns;++i) sum += probs[i];
      for(i=0;i<ns;++i) probs[i] /= sum;
      sa = Sample_i_With_Proba_pi(probs,ns);
    }
  else if(a == tree->n_root)
    {
      sa = Sample_Ancestral_Seq_Core(NULL,tree->n_root,NULL,r_cat,site,tree);
    }
  else
    {
      sa = Assign_State(a->c_seq_anc->state + site,
                        tree->mod->io->datatype,
                        tree->mod->io->state_len);
    }

  sd = Sample_Ancestral_Seq_Core(a,d,b,r_cat,site,tree);
  
  if(fullmutmap == YES) Map_Mutations(a,d,sa,sd,b,site,r_cat,muttype,muttime,muttax,n_mut,tree);
    
  if(d->tax) return;
  else
    {
      for(i=0;i<3;++i)
        {
          if(d->v[i] != a && d->b[i] != tree->e_root)
            {
              Sample_Ancestral_Seq_Pre(d,d->v[i],d->b[i],site,r_cat,muttype,muttime,muttax,n_mut,fullmutmap,fromprior,tree);
            }
        }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int Sample_Ancestral_Seq_Core(t_node *a, t_node *d, t_edge *b, int r_cat, int site, t_tree *tree)
{
  int i,j;
  int ns;

  ns = tree->mod->ns;
  
  int dim1 = tree->mod->ras->n_catg * ns;
  int dim2 = ns;
  int dim3 = ns * ns;
  
  phydbl *probs,*Pij,sum;
  
  int state;

  
  probs = (phydbl *)mCalloc(ns,sizeof(phydbl));
  state = -1;

  if(d->tax == YES)
    {
      assert(b);
      if(d == b->left)      for(i=0;i<ns;++i) probs[i] = b->p_lk_tip_l[site*ns+i];
      else if(d == b->rght) for(i=0;i<ns;++i) probs[i] = b->p_lk_tip_r[site*ns+i];
      else assert(FALSE);
      for(i=0;i<ns;++i) probs[i] *= tree->mod->e_frq->pi->v[i];
      sum = 0.0;
      for(i=0;i<ns;++i) sum += probs[i];
      for(i=0;i<ns;++i) probs[i] /= sum;
      state = Sample_i_With_Proba_pi(probs,ns);
    }
  else
    {
      if(a == NULL) // d is root node
        {
          assert(d == tree->n_root);
          
          phydbl r,l;
          
          Update_PMat_At_Given_Edge(tree->n_root->b[1],tree);
          Update_PMat_At_Given_Edge(tree->n_root->b[2],tree);
          
          for(i=0;i<ns;++i)
            {
              if(tree->e_root->left == tree->n_root->v[1])
                Pij = tree->n_root->b[1]->Pij_rr;
              else
                Pij = tree->n_root->b[2]->Pij_rr;
              
              l = 0.0;
              for(j=0;j<ns;++j)
                {
                  if(tree->e_root->left->tax == NO)
                    l += tree->e_root->p_lk_left[site*dim1+r_cat*dim2+j] * Pij[r_cat*dim3+i*dim2+j];
                  else
                    l += tree->e_root->p_lk_tip_l[site*dim2+j] * Pij[r_cat*dim3+i*dim2+j];
                }
              
              
              if(tree->e_root->rght == tree->n_root->v[1])
                Pij = tree->n_root->b[1]->Pij_rr;
              else
                Pij = tree->n_root->b[2]->Pij_rr;
              
              r = 0.0;
              for(j=0;j<ns;++j)
                {
                  if(tree->e_root->rght->tax == NO)
                    r += tree->e_root->p_lk_rght[site*dim1+r_cat*dim2+j] * Pij[r_cat*dim3+i*dim2+j];
                  else
                    r += tree->e_root->p_lk_tip_r[site*dim2+j] * Pij[r_cat*dim3+i*dim2+j];
                }
              
              probs[i] = r*l*tree->mod->e_frq->pi->v[i];
            }
          
          sum = 0.0;
          for(i=0;i<ns;++i) sum += probs[i];
          for(i=0;i<ns;++i) probs[i] /= sum;
          
          state = Sample_i_With_Proba_pi(probs,ns);
          
          d->c_seq_anc->state[site] = Reciproc_Assign_State(state,tree->io->datatype);
        }
      else
        {
          assert(b);
          
          // State (already) sampled at node a
          state = Assign_State(a->c_seq_anc->state + site,
                               tree->mod->io->datatype,
                               tree->mod->io->state_len);
          
          Pij = b->Pij_rr;
                              
          for(i=0;i<ns;++i)
            {
              if(d == b->left)
                probs[i] = b->p_lk_left[site*dim1+r_cat*dim2+i] * Pij[r_cat*dim3+state*dim2+i];
              else if(d == b->rght)
                probs[i] = b->p_lk_rght[site*dim1+r_cat*dim2+i] * Pij[r_cat*dim3+state*dim2+i];
              else assert(FALSE);
            }
          
          sum = 0.0;
          for(i=0;i<ns;++i) sum += probs[i];
          for(i=0;i<ns;++i) probs[i] /= sum;
          
          state = Sample_i_With_Proba_pi(probs,ns);
          
          d->c_seq_anc->state[site] = Reciproc_Assign_State(state,tree->io->datatype);
        }      
    }

  Free(probs);
  return state;
 
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Map_Mutations(t_node *a, t_node *d, int sa, int sd, t_edge *b, int site, int rcat, int *muttype, phydbl *muttime, int *muttax, int *n_mut, t_tree *tree)
{
  int i,j;
  phydbl *probs,*all_probs;
  int slast,snew; // Last state visited
  phydbl tlast;
  phydbl *Q;
  phydbl u,sum;
  int *mut; // Array of mutations
  int thismut;
  int n_mut_branch;
  int n_iter;
  int first_mut;
  phydbl T;
  int ns,tax_idx;
  phydbl rr;
  
  ns = tree->mod->ns;
  
  all_probs = (phydbl *)mCalloc(ns*ns,sizeof(phydbl));
  mut = (int *)mCalloc(ns*ns,sizeof(int));
  
  // Site relative rate
  rr = tree->mod->ras->gamma_rr->v[rcat];
      
  // Rate matrix
  Q = tree->mod->r_mat->qmat->v;
  
  // Length of the 'time' interval considered.
  T = 0.0;
#ifdef PHYTIME
  T = tree->rates->cur_l[d->num]*rr;
#else
  T = b->l->v*rr;
#endif
  
  
  /* PhyML_Printf("\n. Mutmap: a:%d d:%d ta:%G td:%G cr:%G rr:%G l:%G", */
  /*              a?a->num:-1,d?d->num:-1, */
  /*              tree->times->nd_t[a->num], */
  /*              tree->times->nd_t[d->num], */
  /*              cr,rr, */
  /*              fabs(tree->times->nd_t[a->num]-tree->times->nd_t[d->num])*cr*rr); */

  // Matrix of change probabilities
  for(i=0;i<ns;i++)
    {
      // We only care about the non-diagonal elements here
      for(j=0;j<ns;j++) all_probs[i*ns+j] = Q[i*ns+j];
      
      // Set the diagonal to 0 so that p(i->i)=0.0;
      all_probs[i*ns+i] = 0.0;

      // Normalise so that \sum_j p(i->j) = 1.0;
      sum = 0;
      for(j=0;j<ns;j++) sum += all_probs[i*ns+j];
      for(j=0;j<ns;j++) all_probs[i*ns+j] /= sum;
    }
  
  for(i=0;i<ns*ns;++i) mut[i] = 0;
  tlast = .0;
  slast = sa;
  snew  = sa;
  probs = NULL;
  n_mut_branch = 0;
  n_iter = 0;
  first_mut = YES;
  
  do
    {
      if((sa != sd) && (first_mut == YES)) // ancestral and descendant states are distinct
        {
          // Sample a time for the first mutation conditional on at least one mutation
          // occurring (see formula 2.1 in Hobolth and Stone, 2009).
          u = Uni();
          tlast = -log(1. - u*(1.-exp(Q[sa*ns+sa]*T)))/-Q[sa*ns+sa];
        }
      else
        {
          // Sample a time for the next mutation
          tlast = tlast + Rexp(-Q[slast*ns+slast]);
        }
      
      // Select the appropriate vector of change probabilities
      probs = all_probs+slast*ns;
            

      /* printf("\n. sa=%2d sd=%2d slast=%2d tlast=%12G T=%12G  rcat=%2d site=%4d",sa,sd,slast,tlast,T,rcat,site); */
      
      // The time for the next mutation does not exceed the length
      // of the time interval -> sample a new mutation event
      if(tlast < T)
        {
          first_mut = NO;
          
          n_mut_branch++;

          if(n_mut_branch > 5 && n_mut_branch%10 == 0)
            PhyML_Printf("\n. # of mutations on edge %d (length: %g) exceeds 5 (%d) ! The program will probably crash soon...:(",
                         b->num,
                         tree->rates ? tree->rates->cur_l[d->num] : b->l->v,
                         n_mut_branch);
          
          snew = Sample_i_With_Proba_pi(probs,ns);
          
          // Record mutation type
          mut[slast*ns+snew]++;
          
          // Record mutation type in the site mutation array
          thismut = slast*ns+snew;
          
          muttype[(*n_mut)+n_mut_branch-1] = thismut;

          // Record time of mutation
          muttime[(*n_mut)+n_mut_branch-1] = tlast;

#ifdef PHYTIME
          // Transform into time in calendar units
          muttime[(*n_mut)+n_mut_branch-1] /= tree->rates->cur_l[d->num];
          muttime[(*n_mut)+n_mut_branch-1] *= fabs(tree->times->nd_t[a->num]-tree->times->nd_t[d->num]);
#endif

          tax_idx = -1;
          Random_Tax_Idx(a,d,&tax_idx,tree);
          muttax[(*n_mut)+n_mut_branch-1] = tax_idx;

          // Update the last state
          slast = snew;
        }
      else
        {
          if(slast == sd) break;
          else
            {
              // Restart from the beginning
              for(i=0;i<ns*ns;++i) mut[i] = 0;
              for(i=0;i<n_mut_branch;i++) muttype[(*n_mut)+n_mut_branch-1] = -2;
              for(i=0;i<n_mut_branch;i++) muttime[(*n_mut)+n_mut_branch-1] = +1.;
              tlast = 0.0;
              slast = sa;
              snew = sa;
              n_mut_branch = 0;
              first_mut = YES;
              n_iter++;
            }
        }
    }
  while(++n_iter < 10000);


  if(n_iter == 10000)
    {
      PhyML_Printf("\n. sa=%2d sd=%2d slast=%2d tlast=%12G T=%12G  rcat=%12f site=%4d",sa,sd,slast,tlast,T,rcat,site);
      assert(FALSE);
    }
  
  (*n_mut) += n_mut_branch;
  
  
  /* for(i=0;i<ns;i++) */
  /*   { */
  /*     for(j=i+1;j<ns;j++) */
  /*       { */
  /*         if(mut[i*ns+j] + mut[j*ns+i] > 0) */
  /*           { */
  /*             thismut = MIN(i,j) * ns + MAX(i,j) - (MIN(i,j)+1+(int)POW(MIN(i,j)+1,2))/2; */

  /*             if(tree->mixt_tree != NULL) */
  /*               tree->mixt_tree->mutmap[thismut*(tree->n_pattern)*(2*tree->n_otu-3) + b->num*(tree->n_pattern) + site]++; */
  /*             else */
  /*               tree->mutmap[thismut*(tree->n_pattern)*(2*tree->n_otu-3) + b->num*(tree->n_pattern) + site]++; */
  /*           } */
  /*       } */
  /*   } */
  
  Free(all_probs);
  Free(mut);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void Ancestral_Sequences(t_tree *tree, int print)
{
  int i;

  if(print == YES)
    {

      PhyML_Printf("\n\n. Estimating ancestral sequences...");

      strcpy(tree->io->out_ancestral_seq_file,tree->io->out_file);
      strcat(tree->io->out_ancestral_seq_file,"_phyml_ancestral_");
      if(tree->io->append_run_ID)
        {
          strcat(tree->io->out_ancestral_seq_file,tree->io->run_id_string);
          strcat(tree->io->out_ancestral_seq_file,"_");
        }
      strcat(tree->io->out_ancestral_seq_file,"seq.txt");
      tree->io->fp_out_ancestral_seq = Openfile(tree->io->out_ancestral_seq_file,1);
      

      strcpy(tree->io->out_ancestral_tree_file,tree->io->out_file);
      strcat(tree->io->out_ancestral_tree_file,"_phyml_ancestral_");
      if(tree->io->append_run_ID)
        {
          strcat(tree->io->out_ancestral_tree_file,tree->io->run_id_string);
          strcat(tree->io->out_ancestral_tree_file,"_");
        }
      strcat(tree->io->out_ancestral_tree_file,"tree.txt");
      tree->io->fp_out_ancestral_tree = Openfile(tree->io->out_ancestral_tree_file,1);

      
      PhyML_Fprintf(tree->io->fp_out_ancestral_seq,"\n\n\n");
      PhyML_Fprintf(tree->io->fp_out_ancestral_seq,"\n. Printing marginal probabilities of ancestral sequences at each site");
      PhyML_Fprintf(tree->io->fp_out_ancestral_seq,"\n. of the alignment and each node of the tree. The tree in Newick format");
      PhyML_Fprintf(tree->io->fp_out_ancestral_seq,"\n. with internal nodes labels corresponding to those given below can be");
      PhyML_Fprintf(tree->io->fp_out_ancestral_seq,"\n. found in the file '%s'.",tree->io->out_ancestral_tree_file);
      PhyML_Fprintf(tree->io->fp_out_ancestral_seq,"\n");
      PhyML_Fprintf(tree->io->fp_out_ancestral_seq,"\n");
      PhyML_Fprintf(tree->io->fp_out_ancestral_seq,"\n. Ancestral reconstruction is conducted based on the \"Minimum Posterior Expected");
      PhyML_Fprintf(tree->io->fp_out_ancestral_seq,"\n. Error\" (MPEE) criterion, which accomodates for uncertainty in the selection of ");
      PhyML_Fprintf(tree->io->fp_out_ancestral_seq,"\n. ancestral character states.");
      PhyML_Fprintf(tree->io->fp_out_ancestral_seq,"\n");
      PhyML_Fprintf(tree->io->fp_out_ancestral_seq,"\n");
      PhyML_Fprintf(tree->io->fp_out_ancestral_seq,"\n. Recommended citation:");
      PhyML_Fprintf(tree->io->fp_out_ancestral_seq,"\n. Oliva A., Pulicani S., Lefort V., Brehelin L., Gascuel O. and S. Guindon,");
      PhyML_Fprintf(tree->io->fp_out_ancestral_seq,"\n. \"Accounting for ambiguity in ancestral sequence reconstruction\",");
      PhyML_Fprintf(tree->io->fp_out_ancestral_seq,"\n. Bioinformatics, Volume 35, Issue 21, 1 November 2019.");
      PhyML_Fprintf(tree->io->fp_out_ancestral_seq,"\n");
      
      PhyML_Fprintf(tree->io->fp_out_ancestral_seq,"\n\n");
      PhyML_Fprintf(tree->io->fp_out_ancestral_seq,"Site\tNodeLabel\t");
      for(i=0;i<tree->mod->ns;i++) PhyML_Fprintf(tree->io->fp_out_ancestral_seq,"%10c\t",Reciproc_Assign_State(i,tree->io->datatype));
      PhyML_Fprintf(tree->io->fp_out_ancestral_seq,"MPEE\t");
      PhyML_Fprintf(tree->io->fp_out_ancestral_seq,"\n");

      
      short int bck_support = tree->io->print_support_val;
      
      tree->io->print_node_num = YES;
      tree->io->print_support_val = NO;
      char *s = Write_Tree(tree);
      PhyML_Fprintf(tree->io->fp_out_ancestral_tree,"%s",s);
      tree->io->print_node_num = NO;
      tree->io->print_support_val = bck_support;

      Free(s);
      fclose(tree->io->fp_out_ancestral_tree);
    }

  for(i=0;i<2*tree->n_otu-2;i++)
    if(tree->a_nodes[i]->tax == NO)
      Ancestral_Sequences_One_Node(tree->a_nodes[i],tree,print);

  if(tree->n_root) Ancestral_Sequences_One_Node(tree->n_root,tree,print);


  fclose(tree->io->fp_out_ancestral_seq);
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
          
          fp = tree->io->fp_out_ancestral_seq;
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

              if(tree->mod->ras->invar == YES)
                {
                  int num_prec_issue = NO;
                  phydbl inv_site_lk = Invariant_Lk(sum_scale,csite,&num_prec_issue,tree);


                  switch(num_prec_issue)
                    {
                    case YES :         
                      {
                        assert(isinf(inv_site_lk));
                        inv_site_lk = Invariant_Lk(0,csite,&num_prec_issue,tree);
                        for(i=0;i<ns;i++) p[i] = inv_site_lk * tree->mod->ras->pinvar->v * tree->mod->e_frq->pi->v[i];
                        break;
                      }
                    case NO : 
                      {
                        for(i=0;i<ns;i++) p[i] = p[i] * (1.-tree->mod->ras->pinvar->v) + inv_site_lk * tree->mod->ras->pinvar->v * tree->mod->e_frq->pi->v[i];
                        break;
                      }
                    }

                }

              
              for(i=0;i<ns;i++) p[i] = log(p[i]) - (phydbl)LOG2 * sum_scale;
              for(i=0;i<ns;i++) p[i] -= tree->c_lnL_sorted[csite];
              for(i=0;i<ns;i++) p[i] = exp(p[i]);


              
              /* sum_probas = 0.0; */
              /* for(i=0;i<ns;i++) sum_probas += p[i]; */
              /* for(i=0;i<ns;i++) p[i]/=sum_probas; */
              
              sum_probas = 0.0;
              for(i=0;i<ns;i++) sum_probas += p[i];
              if(Are_Equal(sum_probas,1.0,0.01) == NO)
                {
                  PhyML_Fprintf(stderr,"\n. Probabilities do not sum to 1.0! Aborting.");
                  for(i=0;i<ns;++i) PhyML_Fprintf(stderr,"\n. p[%2d]=%G",i,p[i]);
                  Exit("\n");
                }
              
              if(print == YES)
                {
                  PhyML_Fprintf(fp,"%4d\t%9d\t",site+1,d->num);                  
                  for(i=0;i<ns;i++) PhyML_Fprintf(fp,"%10g\t",p[i]);
                  PhyML_Fprintf(fp,"%s",Bit_To_Character_String(Integer_To_Bit(MPEE_Infer(p,ns),ns),ns));
                  PhyML_Fprintf(fp,"\n");
                }

              
              
            }
          Free(p);
        }
    }
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

int MPEE_Infer(const phydbl *p, const int ns)
{
  unsigned int i,j,sorted;
  phydbl *alpha;
  const unsigned mesh = 50;
  int *best_state,*idx,*count,highest_count,highest_count_idx,most_frequent_state;
  
  idx = (int *)mCalloc(ns,sizeof(int));
  
  alpha = (phydbl *)mCalloc(mesh,sizeof(phydbl));
  best_state = (int *)mCalloc(mesh+1,sizeof(int));
  count = (int *)mCalloc(mesh+1,sizeof(int));

  
  // idx[0] gives the index of the largest value in p.
  // idx[1] gives the index of the second largest value in p
  // etc.
  for(i=0;i<ns;++i) idx[i] = i;
  do
    {
      sorted = YES;
      for(i=0;i<ns-1;++i)
        {
          if(p[idx[i]] < p[idx[i+1]])
            {
              sorted = NO;

              j        = idx[i];
              idx[i]   = idx[i+1];
              idx[i+1] = j;
            }
        }
    }
  while(sorted == NO);

  /* for(i=0;i<ns;++i) */
  /*   { */
  /*     PhyML_Printf("\n. p[%d]=%g",i,p[i]); */
  /*   } */

  /* for(i=0;i<ns;++i) */
  /*   { */
  /*     PhyML_Printf("\n. idx[%d]=%d",i,idx[i]); */
  /*   } */

  for(i=0;i<mesh+1;++i)
    {
      for(j=0;j<ns;++j) alpha[j] = i*((phydbl)j/(j+1))/mesh;
      best_state[i] = MPEE_Score(alpha,idx,p,ns);
    }

  
  /* for(i=0;i<mesh;++i) */
  /*   { */
  /*     PhyML_Printf("\n. best_state[%d]=%d",i,best_state[i]); */
  /*   } */

  for(i=0;i<mesh+1;++i)
    {
      count[i] = 0;
      for(j=0;j<mesh+1;++j)
        {
          if(best_state[i] == best_state[j]) count[i]++;
        }
    }

  highest_count = 0;
  highest_count_idx = -1;
  for(i=0;i<mesh+1;++i)
    {
      if(count[i] > highest_count)
        {
          highest_count = count[i];
          highest_count_idx = i;
        }
    }
  
  most_frequent_state = best_state[highest_count_idx];
  
  Free(alpha);
  Free(best_state);
  Free(idx);
  Free(count);
  
  return(most_frequent_state);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

int MPEE_Score(const phydbl *alpha, int *idx, const phydbl *p, const int ns)
{
  unsigned int i;
  phydbl *mpee_score,*cdf_sorted;
  int min_idx,res;
  phydbl a,b;
  phydbl min;
    
  mpee_score = (phydbl *)mCalloc(ns,sizeof(phydbl));
  cdf_sorted = (phydbl *)mCalloc(ns,sizeof(phydbl));

  
  cdf_sorted[0] = p[idx[0]];
  for(i=1;i<ns;++i) cdf_sorted[i] = cdf_sorted[i-1] + p[idx[i]];

  for(i=0;i<ns-1;++i) // for the different levels of ambiguity
    {
      a = alpha[i];
      b = (phydbl)(ns-1.-a*(i+1))/(ns-i-1);

      mpee_score[i] = a + (b-a)*(1-cdf_sorted[i]);
    }
  mpee_score[ns-1] = (phydbl)(ns-1.)/ns;
  
  min = mpee_score[0];
  min_idx = 0;
  for(i=1;i<ns;++i)
    if(mpee_score[i] < min)
      {
        min = mpee_score[i];
        min_idx = i;
      }
  
  res = 0;
  for(i=0;i<min_idx+1;++i)
    {
      res += (int)pow(2,ns-1-idx[i]);
    }
  /* If best score ambiguity level is say 2 and A and G have highest
     posterior probs, then res = 2^(4-1-0) + 2^(4-1-2) = 10, which is
     bit representation is [1,0,1,0]
  */
  
  Free(mpee_score);
  Free(cdf_sorted);
  
  return(res);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void Reverse_Muttime(phydbl *muttime, int n_mut, t_tree *tree)
{
  phydbl h;
  int i;  
  h = Tree_Height(tree);
  for(i=0;i<n_mut;++i) muttime[i] = h - muttime[i];
  
}


/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
