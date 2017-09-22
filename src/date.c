/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/


/* Routines for molecular dating */


#include "date.h"

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int DATE_Main(int argc, char **argv)
{
  option *io;

  io = Get_Input(argc,argv);
  Free(io);
  return(0);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// Update t_prior_min and t_prior_max on a given ranked tree
// given (primary and secondary) calibration information.
// Make sure secondary and primary calibration are up-to-date
void DATE_Update_T_Prior_MinMax(t_tree *tree)
{
  int i,j;

  /* printf("\n"); */
  
  for(i=tree->n_otu;i<2*tree->n_otu-1;++i) // All internal nodes 
    {
      tree->rates->t_prior_max[i] = 0.0;
      tree->rates->t_prior_min[i] = -INFINITY;

      if(tree->a_nodes[i]->n_cal > 0) // Primary calibration found on that node
        {
          for(j=0;j<tree->a_nodes[i]->n_cal;++j)
            {
              tree->rates->t_prior_max[i] = MIN(tree->rates->t_prior_max[i],tree->a_nodes[i]->cal[j]->upper);
              tree->rates->t_prior_min[i] = MAX(tree->rates->t_prior_min[i],tree->a_nodes[i]->cal[j]->lower);
            }         
        }
    }

  /* // TO DO: chain t_rank */
  /* TIMES_Update_Node_Ordering(tree); */
  /* rk = tree->rates->t_rank; */

  /* for(i=0;i<tree->n_otu-2;i++) */
  /*   { */
  /*     /\* printf("\n> [%3d] t:%f|%f -- min:%f|%f", *\/ */
  /*     /\*        tree->a_nodes[rk[i+1]]->num, *\/ */
  /*     /\*        tree->rates->nd_t[rk[i]], *\/ */
  /*     /\*        tree->rates->nd_t[rk[i+1]], *\/ */
  /*     /\*        tree->rates->t_prior_min[rk[i]], *\/ */
  /*     /\*        tree->rates->t_prior_min[rk[i+1]]); *\/ */

  /*     if(tree->rates->t_prior_min[rk[i+1]] < tree->rates->t_prior_min[rk[i]]) */
  /*       { */
  /*         /\* printf("  --> %f",tree->rates->t_prior_min[rk[i]]); *\/ */
  /*         tree->rates->t_prior_min[rk[i+1]] = tree->rates->t_prior_min[rk[i]]; */
  /*       } */
  /*   } */
  

  /* for(i=tree->n_otu-2;i>0;i--) */
  /*   { */
  /*     /\* printf("\n< [%3d] t:%f|%f -- max:%f|%f", *\/ */
  /*     /\*        tree->a_nodes[rk[i-1]]->num, *\/ */
  /*     /\*        tree->rates->nd_t[rk[i]], *\/ */
  /*     /\*        tree->rates->nd_t[rk[i-1]], *\/ */
  /*     /\*        tree->rates->t_prior_max[rk[i]], *\/ */
  /*     /\*        tree->rates->t_prior_max[rk[i-1]]); *\/ */
  /*     if(tree->rates->t_prior_max[rk[i-1]] > tree->rates->t_prior_max[rk[i]]) */
  /*       { */
  /*         /\* printf("  --> %f ",tree->rates->t_prior_max[rk[i]]); *\/ */
  /*         tree->rates->t_prior_max[rk[i-1]] = tree->rates->t_prior_max[rk[i]]; */
  /*       } */
  /*   } */
    
  /* /\* for(i=tree->n_otu;i<2*tree->n_otu-1;i++) *\/ */
  /* /\*   { *\/ */
  /* /\*     if(tree->rates->t_prior_min[i] > tree->rates->t_prior_max[i]) *\/ */
  /* /\*       { *\/ */
  /* /\*         PhyML_Printf("\n. i: %d t_prior_min: %f t_prior_max: %f", *\/ */
  /* /\*                      i, *\/ */
  /* /\*                      tree->rates->t_prior_min[i], *\/ */
  /* /\*                      tree->rates->t_prior_max[i]); *\/ */
  /* /\*         Generic_Exit(__FILE__,__LINE__,__FUNCTION__); *\/ */
  /* /\*       } *\/ */
  /* /\*   } *\/ */
  
  /* /\* printf("\n. min:%f max:%f rk: %d n_cal: %d low: %f up: %f", *\/ */
  /* /\*        tree->rates->t_prior_min[tree->n_root->num], *\/ */
  /* /\*        tree->rates->t_prior_max[tree->n_root->num], *\/ */
  /* /\*        tree->rates->t_rank[tree->n_root->num], *\/ */
  /* /\*        tree->a_nodes[tree->n_root->num]->n_cal, *\/ */
  /* /\*        tree->a_nodes[tree->n_root->num]->n_cal > 0 ? tree->a_nodes[tree->n_root->num]->cal[0]->lower : -1, *\/ */
  /* /\*        tree->a_nodes[tree->n_root->num]->n_cal > 0 ? tree->a_nodes[tree->n_root->num]->cal[0]->upper : -1); *\/ */
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void DATE_Assign_Primary_Calibration(t_tree *tree)
{
  int i,j,idx,node_num;
  t_clad *clade;
  t_cal *cal;
  
  clade = NULL;
  cal  = NULL;
  
  for(i=0;i<tree->rates->n_cal;++i)
    {
      cal = tree->rates->a_cal[i];
      clade = cal->clade_list[cal->current_clade_idx];
      clade->target_nd = NULL;
    }
  
  for(i=0;i<2*tree->n_otu-1;++i) 
    for(j=0;j<MAX_N_CAL;++j) 
    {
      tree->a_nodes[i]->cal[j] = NULL;
      tree->a_nodes[i]->n_cal  = 0;
    }
  
  for(i=0;i<tree->rates->n_cal;++i)
    {
      cal = tree->rates->a_cal[i];
      clade = cal->clade_list[cal->current_clade_idx];

      node_num = Find_Clade(clade->tax_list,
                            clade->n_tax,
                            tree);

      clade->target_nd = tree->a_nodes[node_num];
      
      idx = tree->a_nodes[node_num]->n_cal;
      tree->a_nodes[node_num]->cal[idx] = tree->rates->a_cal[i];
      tree->a_nodes[node_num]->n_cal++;
      

      if(tree->a_nodes[node_num]->n_cal == MAX_N_CAL)
        {
          PhyML_Fprintf(stderr,"\n. A node cannot have more than %d calibration",MAX_N_CAL); 
          PhyML_Fprintf(stderr,"\n. constraints attached to it. Feel free to increase the"); 
          PhyML_Fprintf(stderr,"\n. value of the variable MAX_N_CAL in utilities.h if");
          PhyML_Fprintf(stderr,"\n. necessary.");
          Exit("\n");
        }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// Return splitted calibration intervals. Make sure all primary
// and secondary calibration intervals are up-to-date.
phydbl *DATE_Splitted_Calibration(t_tree *tree)
{
  phydbl *minmax,*splitted_cal,buff;
  int i,len,done;

  // One t_prior_min and one t_prior_max per internal nodes except root, so 
  // 2 x # of internal nodes boundaries in total at most.
  minmax = (phydbl *)mCalloc(2*(tree->n_otu-2),sizeof(phydbl));
  For(i,2*(tree->n_otu-2)) minmax[i] = +INFINITY;
  splitted_cal = (phydbl *)mCalloc((int)(4*tree->n_otu-10),sizeof(phydbl));


  len = 0;
  for(i = tree->n_otu; i < 2*tree->n_otu-1; i++)
    {
      if(tree->a_nodes[i] != tree->n_root)
        {
          minmax[len]   = MAX(tree->rates->t_prior_min[i],tree->rates->nd_t[tree->n_root->num]);
          minmax[len+1] = tree->rates->t_prior_max[i];        
          len+=2;
        }
    }

  
  // Bubble sort of all these times in increasing order
  do
    {
      done = YES;
      for(i=0;i<len-1;i++) 
        {
          if(minmax[i] > minmax[i+1])
            {
              buff        = minmax[i];
              minmax[i]   = minmax[i+1];
              minmax[i+1] = buff;
              done = NO;
            }
        }
    }
  while(done == NO);

  for(i=0;i<len-1;i++) assert(!(minmax[i] > minmax[i+1]));
  
  // Remove ties
  for(i=0;i<len-1;i++) 
    if(Are_Equal(minmax[i],minmax[i+1],1.E-6) == YES) 
      minmax[i] = 0.0;

  // Sort again to effectively remove ties
  do
    {
      done = YES;
      for(i=0;i<len-1;i++) 
        {
          if(minmax[i] > minmax[i+1])
            {
              buff        = minmax[i];
              minmax[i]   = minmax[i+1];
              minmax[i+1] = buff;
              done = NO;
            }
        }
    }
  while(done == NO);

  splitted_cal[0] = minmax[0];
  len = 1;
  for(i = 1; i < 2*(tree->n_otu-2); i++)                                        
    {
      splitted_cal[len]   = minmax[i];
      if(len+1 < 4*tree->n_otu-10) splitted_cal[len+1] = minmax[i];
      len+=2;
    }

  /* For(i,4*tree->n_otu-10) PhyML_Printf("\n. split -- %3d %12f",i,splitted_cal[i]); */

  Free(minmax);
   
  return splitted_cal;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl DATE_J_Sum_Product(t_tree *tree)
{
  phydbl prod, total,*splitted_cal;
  int fact,idx,ans,rk;
  
  DATE_Assign_Primary_Calibration(tree);
  DATE_Update_T_Prior_MinMax(tree);

  splitted_cal = DATE_Splitted_Calibration(tree);
  
  ans   = 0;
  total = 0.0;
  idx   = 0;
  rk    = 1;
  do
    {
      prod = 1.0;
      fact = 1;
      ans = DATE_J_Sum_Product_Pre(tree->a_nodes[tree->rates->t_rank[rk]], // Oldest node after root (as rk=1)
                                   idx,
                                   -1,
                                   prod,fact,&total,splitted_cal,rk,tree);
      idx+=2;
    }
  while(ans != 1);

  Free(splitted_cal);

  return(total);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int DATE_J_Sum_Product_Pre(t_node *d, int split_idx_d, int split_idx_a, phydbl prod, int fact, phydbl *total, phydbl *splitted_cal, int rk, t_tree *tree)
{
  int ans,idx;

  ans = DATE_Is_Split_Accessible(d,split_idx_d,splitted_cal,tree);

  switch(ans)
    {
    case 1 : // split interval is younger than t_prior_max. No need to go further.
      {
        return ans;
        break;
      }
    case 0 : // split interval is within [t_prior_min,t_prior_max]
      {
        int local_ans;

        // Calculate J for this time interval
        prod *= DATE_J(tree->rates->birth_rate, 
                       tree->rates->death_rate,
                       FABS(splitted_cal[split_idx_d+1]),
                       FABS(splitted_cal[split_idx_d]));
        
        // Remove factorial term from current product
        prod *= fact;
        
        if(split_idx_d == split_idx_a) fact++;
        else fact = 1;
        
        prod /= fact;
                
        if(tree->rates->t_rank[tree->n_otu-2] == d->num) // Youngest internal node 
          {
            (*total) += prod;
            return 0;
          }

        idx = split_idx_d;
        do
          {
            local_ans = DATE_J_Sum_Product_Pre(tree->a_nodes[tree->rates->t_rank[rk+1]],
                                               idx,
                                               split_idx_d,
                                               prod,fact,total,splitted_cal,rk+1,tree);
            idx+=2;
          }
        while(local_ans == 0);

        break;
      }
    case -1 : // split interval is older than t_prior_min. Move forward.
      {
        int local_ans;
        // Advance to younger split intervals and stop once you're in
        idx = split_idx_d+2;
        do
          {
            local_ans = DATE_J_Sum_Product_Pre(d,
                                               idx,
                                               idx+1,
                                               prod,fact,total,splitted_cal,rk,tree);
            idx+=2;
          }
        while(local_ans == -1);
        break;
      }
    }

  return ans;
}
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int DATE_Is_Split_Accessible(t_node *d, int which, phydbl *splitted_cal, t_tree *tree)
{
  phydbl eps;

  assert(d->tax == NO);

  eps = FABS(tree->rates->t_prior_min[d->num]) / 1.E+6;

  assert(eps > MDBL_MIN);

  // Upper and lower bound of splitted calibration interval are equal to zero
  if(Are_Equal(splitted_cal[which],0.0,eps) && Are_Equal(splitted_cal[which+1],0.0,eps)) return +1;


  if(Are_Equal(tree->rates->t_prior_min[d->num],splitted_cal[which],eps)   ||
     Are_Equal(tree->rates->t_prior_max[d->num],splitted_cal[which+1],eps) ||
     (tree->rates->t_prior_min[d->num] < splitted_cal[which] && 
      tree->rates->t_prior_max[d->num] > splitted_cal[which+1]))                   return  0; // splitted interval is within [t_prior_min,t_prior_max]
  else if(Are_Equal(tree->rates->t_prior_max[d->num],splitted_cal[which],eps) ||
          splitted_cal[which] > tree->rates->t_prior_max[d->num])                  return +1; // splitted interval is younger than [t_prior_min,t_prior_max]
  else if(Are_Equal(tree->rates->t_prior_min[d->num],splitted_cal[which+1],eps) ||
          splitted_cal[which+1] < tree->rates->t_prior_min[d->num])                return -1; // splitted interval is older than [t_prior_min,t_prior_max]
  else
    {
      PhyML_Printf("\n. d->num: %d d->tax: %d",d->num,d->tax);
      PhyML_Printf("\n. t_prior_min: %f t_prior_max: %f",
                   tree->rates->t_prior_min[d->num],
                   tree->rates->t_prior_max[d->num]);
      PhyML_Printf("\n. splitted_cal_min: %f splitted_cal_max: %f",
                   splitted_cal[which],
                   splitted_cal[which+1]);      
      PhyML_Printf("\n");
      assert(FALSE); // splitted interval cannot be partially overlapping [t_prior_min,t_prior_max]
    }
  return(0);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl DATE_J(phydbl birth_r, phydbl death_r, phydbl t_min, phydbl t_pls)
{
  phydbl d,b,J;
  assert(t_pls > t_min);
  d = death_r;
  b = birth_r;
  J = (b-d)*(exp(t_min*d+t_pls*b) - exp(t_min*b+t_pls*d));
  J /= ((b*exp(t_min*b)-d*exp(t_min*d)) * (b*exp(t_pls*b)-d*exp(t_pls*d)));
  /* printf("  J : %f",J); */
  return(J);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int DATE_Check_Calibration_Constraints(t_tree *tree)
{
  int i,j;
  phydbl lower,upper;

  lower = upper = -1.;

  For(i,2*tree->n_otu-1)
    {
      if(tree->a_nodes[i]->n_cal > 1)
        {
          lower = tree->a_nodes[i]->cal[0]->lower;
          upper = tree->a_nodes[i]->cal[0]->upper;
          for(j=1; j < tree->a_nodes[i]->n_cal; j++)
            {
              lower = MAX(lower,tree->a_nodes[i]->cal[j]->lower);
              upper = MIN(upper,tree->a_nodes[i]->cal[j]->upper);
              if(upper < lower) 
                {
                  /* PhyML_Printf("\n. Inconsistency detected on node %d",i); */
                  return 0; 
                }
            }
        }
    }
  return 1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int DATE_Check_Time_Constraints(t_tree *tree)
{
  int i;

  For(i,2*tree->n_otu-1)
    {
      if(tree->a_nodes[i] != tree->n_root && tree->a_nodes[i]->tax == NO)
        {
          if(tree->rates->nd_t[i] > tree->rates->t_prior_max[i] ||
             tree->rates->nd_t[i] < tree->rates->t_prior_min[i])
            {
              /* PhyML_Printf("\n. Node %d t: %f min:%f max:%f",i,tree->rates->nd_t[i],tree->rates->t_prior_min[i],tree->rates->t_prior_max[i]); */
              return 0;
            }
        }
    }
  return 1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl *DATE_MCMC(t_tree *tree)
{
  t_mcmc *mcmc;
  char *s_tree;
  int move, n_vars, i, j, adjust_len;
  phydbl u;
  phydbl *res;
  FILE *fp_stats,*fp_tree;
  int t_beg;

  t_beg = (int)time(NULL);

  fp_stats = tree->io->fp_out_stats;
  fp_tree = tree->io->fp_out_tree;

  TIMES_Randomize_Tree_With_Time_Constraints(tree->rates->a_cal[0],tree);
  
  mcmc = MCMC_Make_MCMC_Struct();
  tree->mcmc = mcmc;
  MCMC_Init_MCMC_Struct(NULL,NULL,mcmc);
  MCMC_Complete_MCMC(mcmc,tree);
  
  MCMC_Randomize_Birth(tree);
  MCMC_Randomize_Death(tree);
  MCMC_Randomize_Clock_Rate(tree);
  MCMC_Randomize_Rate_Across_Sites(tree);
  MCMC_Randomize_Rates(tree);
  
  n_vars                  = 10;
  adjust_len              = tree->io->mcmc->chain_len_burnin;
  mcmc->sample_interval   = tree->io->mcmc->sample_interval;
  mcmc->print_every       = tree->io->mcmc->print_every;
  mcmc->chain_len         = tree->io->mcmc->chain_len;
  mcmc->chain_len_burnin  = tree->io->mcmc->chain_len_burnin;
  tree->rates->bl_from_rt = YES;

  res = (phydbl *)mCalloc(tree->mcmc->chain_len / tree->mcmc->sample_interval * n_vars,sizeof(phydbl));
  
  Set_Both_Sides(YES,tree);
  Switch_Eigen(YES,tree->mod);
  Lk(NULL,tree);
  Switch_Eigen(NO,tree->mod);
  RATES_Lk_Rates(tree);
  DATE_Assign_Primary_Calibration(tree);
  TIMES_Lk_Times(NO,tree);

  /* Time_To_Branch(tree); */
  /* tree->bl_ndigits = 1; */
  /* printf("\n. Random init tree: %s",Write_Tree(tree,NO)); */
  /* tree->bl_ndigits = 7; */
  RATES_Update_Cur_Bl(tree);

  
  PhyML_Printf("\n. log(Pr(Seq|Tree)) = %f",tree->c_lnL);
  PhyML_Printf("\n. log(Pr(Tree)) = %f",tree->rates->c_lnL_times);
    
  tree->extra_tree = Make_Tree_From_Scratch(tree->n_otu,tree->data);
  tree->extra_tree->mod = tree->mod;
  Copy_Tree(tree,tree->extra_tree);
  tree->extra_tree->rates = RATES_Make_Rate_Struct(tree->n_otu);
  RATES_Init_Rate_Struct(tree->extra_tree->rates,NULL,tree->n_otu);
  RATES_Copy_Rate_Struct(tree->rates,tree->extra_tree->rates,tree->n_otu);
  RATES_Duplicate_Calib_Struct(tree,tree->extra_tree);
  MIXT_Chain_Cal(tree->extra_tree);
  DATE_Assign_Primary_Calibration(tree->extra_tree);
  TIMES_Randomize_Tree_With_Time_Constraints(tree->extra_tree->rates->a_cal[0],tree->extra_tree);
  TIMES_Lk_Times(NO,tree->extra_tree);
  PhyML_Printf("\n. log(Pr(extra tree)) = %f",tree->extra_tree->rates->c_lnL_times);
  mcmc = MCMC_Make_MCMC_Struct();
  tree->extra_tree->mcmc = mcmc;
  MCMC_Init_MCMC_Struct(NULL,NULL,mcmc);
  MCMC_Complete_MCMC(mcmc,tree->extra_tree);

 
  PhyML_Fprintf(fp_stats,"\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t",
                "sample",
                "lnL(posterior)",
                "lnL(seq)",
                "lnL(times)",
                "lnL(rates)",
                "birth",
                "death",
                "clock",
                "root",
                "tstv",
                "nu");
  for(i=0;i<tree->mod->ras->n_catg;i++) PhyML_Fprintf(fp_stats,"rr%d\t",i);
  for(i=0;i<tree->mod->ras->n_catg;i++) PhyML_Fprintf(fp_stats,"pr%d\t",i);


  for(i=0;i<tree->rates->n_cal;++i)
    {
      t_cal *cal = tree->rates->a_cal[i];
      for(j=0;j<cal->clade_list_size;++j)
        {
          t_clad *clade = cal->clade_list[j];
          PhyML_Fprintf(fp_stats,"t(calib:%s_clade:%s)\t",cal->id,clade->id);
        }
    }

  for(i=0;i<tree->rates->n_cal;++i)
    {
      t_cal *cal = tree->rates->a_cal[i];
      PhyML_Fprintf(fp_stats,"clade(calib:%s)\t",cal->id);
    }


  
  for(i=0;i<2*tree->n_otu-2;++i) PhyML_Fprintf(fp_stats,"br%d\t",i);
  for(i=0;i<2*tree->n_otu-1;++i) if(tree->a_nodes[i]->tax == NO) PhyML_Fprintf(fp_stats,"t%d\t",i);
  
  PhyML_Fprintf(fp_stats,"accRT\t");
  PhyML_Fprintf(fp_stats,"tuneRT\t");

  PhyML_Fprintf(fp_stats,"accT\t");
  PhyML_Fprintf(fp_stats,"tuneT\t");

  PhyML_Fprintf(fp_stats,"accClock\t");
  PhyML_Fprintf(fp_stats,"tuneClock\t");

  PhyML_Fprintf(fp_stats,"accTCr\t");
  PhyML_Fprintf(fp_stats,"tuneTCr\t");

  PhyML_Fprintf(fp_stats,"accTreeRates\t");
  PhyML_Fprintf(fp_stats,"tuneTreeRates\t");

  PhyML_Fprintf(fp_stats,"accSprW\t");
  PhyML_Fprintf(fp_stats,"tuneSprW\t");

  PhyML_Fprintf(fp_stats,"accSpr\t");
  PhyML_Fprintf(fp_stats,"tuneSpr\t");

  PhyML_Fprintf(fp_stats,"accSprLoc\t");
  PhyML_Fprintf(fp_stats,"tuneSprLoc\t");

  fflush(NULL);

  PhyML_Printf("\n\n");
  PhyML_Printf("\n. MCMC settings. Chain length:  %g steps.",(phydbl)tree->mcmc->chain_len);
  PhyML_Printf("\n. MCMC settings. Burnin:  %g steps.",(phydbl)tree->mcmc->chain_len_burnin);
  PhyML_Printf("\n. MCMC settings. Sample every %g steps.",(phydbl)tree->mcmc->sample_interval);
  PhyML_Printf("\n. MCMC settings. Print out every %g steps.",(phydbl)tree->mcmc->print_every);
  PhyML_Printf("\n\n");


  for(i=0;i<tree->mcmc->n_moves;i++) tree->mcmc->start_ess[i] = YES;
  Set_Both_Sides(NO,tree);
  tree->mcmc->always_yes = NO;
  move                   = -1;
  
  // Get in the range of sensible values for clock_r
  i = 0;
  do { MCMC_Clock_R(tree); } while(i++ < 100);


  
  
  do
    {
      if(tree->mcmc->run > adjust_len) for(i=0;i<tree->mcmc->n_moves;i++) tree->mcmc->adjust_tuning[i] = NO;

      if(tree->c_lnL < UNLIKELY + 0.1)
        {
          PhyML_Printf("\n. Move '%s' failed\n",tree->mcmc->move_name[move]);
          assert(FALSE);
        }

      u = Uni();
      for(move=0;move<tree->mcmc->n_moves;move++) if(tree->mcmc->move_weight[move] > u-1.E-10) break;
      
      assert(!(move == tree->mcmc->n_moves));      
      
      if(!strcmp(tree->mcmc->move_name[move],"clock"))                   MCMC_Clock_R(tree);
      else if(!strcmp(tree->mcmc->move_name[move],"birth_rate"))         MCMC_Birth_Rate(tree);
      else if(!strcmp(tree->mcmc->move_name[move],"death_rate"))         MCMC_Death_Rate(tree);
      else if(!strcmp(tree->mcmc->move_name[move],"birth_death_updown")) MCMC_Birth_Death_Updown(tree);
      else if(!strcmp(tree->mcmc->move_name[move],"tree_height"))        MCMC_Tree_Height(tree);
      else if(!strcmp(tree->mcmc->move_name[move],"times"))              MCMC_Times_All(tree);
      else if(!strcmp(tree->mcmc->move_name[move],"times_and_rates"))    MCMC_Times_And_Rates_All(tree);
      else if(!strcmp(tree->mcmc->move_name[move],"spr"))                MCMC_Prune_Regraft(tree);
      else if(!strcmp(tree->mcmc->move_name[move],"spr_local"))          MCMC_Prune_Regraft_Local(tree);
      else if(!strcmp(tree->mcmc->move_name[move],"spr_weighted"))       MCMC_Prune_Regraft_Weighted(tree);
      else if(!strcmp(tree->mcmc->move_name[move],"updown_t_cr"))        MCMC_Updown_T_Cr(tree);
      else if(!strcmp(tree->mcmc->move_name[move],"kappa"))              MCMC_Kappa(tree);
      else if(!strcmp(tree->mcmc->move_name[move],"ras"))                MCMC_Rate_Across_Sites(tree);
      else if(!strcmp(tree->mcmc->move_name[move],"nu"))                 MCMC_Nu(tree);
      else if(!strcmp(tree->mcmc->move_name[move],"subtree_height"))     MCMC_Subtree_Height(tree);
      else if(!strcmp(tree->mcmc->move_name[move],"time_slice"))         MCMC_Time_Slice(tree);
      else if(!strcmp(tree->mcmc->move_name[move],"br_rate"))            MCMC_Rates_All(tree);
      else if(!strcmp(tree->mcmc->move_name[move],"tree_rates"))         MCMC_Tree_Rates(tree);
      else if(!strcmp(tree->mcmc->move_name[move],"clade_change"))       MCMC_Clade_Change(tree);
      else continue;
      
      if(!TIMES_Check_Node_Height_Ordering(tree))
        {
          PhyML_Fprintf(stderr,"\n. move: %s",tree->mcmc->move_name[move]);
          Exit("\n");
        }
      
      
      if(!(tree->rates->c_lnL_times > UNLIKELY))
        {
          PhyML_Fprintf(stderr,"\n. move: %s",tree->mcmc->move_name[move]);
          PhyML_Fprintf(stderr,"\n. glnL=%f",tree->rates->c_lnL_times);
          assert(FALSE);
        }

      (void)signal(SIGINT,MCMC_Terminate);
      
      if(!(tree->mcmc->run%tree->mcmc->print_every))
        {
          phydbl mean_r,post;

          mean_r = RATES_Average_Substitution_Rate(tree);
          post = Get_Lk(tree) + tree->rates->c_lnL_times + tree->rates->c_lnL_rates;

          if(tree->mcmc->run < adjust_len) PhyML_Printf("\n.");
          else PhyML_Printf("\n.");
          PhyML_Printf(" %10d lnL: [%12.2f -- %12.2f -- %12.2f -- %12.2f] clock: %12f root age: %12f [time: %7d sec run: %8d] mean:%15f %20s",
                       tree->mcmc->run,
                       post,
                       Get_Lk(tree),
                       tree->rates->c_lnL_times,
                       tree->rates->c_lnL_rates,
                       tree->rates->clock_r,
                       tree->rates->nd_t[tree->n_root->num],                       
                       (int)time(NULL) - t_beg,
                       tree->mcmc->run,
                       mean_r,
                       tree->mcmc->move_name[move]);
        }
      
      if(!(tree->mcmc->run%tree->mcmc->sample_interval))
        {
          phydbl mean_r,post;
          
          mean_r = RATES_Average_Substitution_Rate(tree);
          post = Get_Lk(tree) + tree->rates->c_lnL_times + tree->rates->c_lnL_rates;
          
          PhyML_Fprintf(fp_stats,"\n%6d\t%9.1f\t%9.1f\t%9.1f\t%9.1f\t%12G\t%12G\t%12G\t%12G\t%12G\t%12G\t",
                        tree->mcmc->run,
                        post,
                        Get_Lk(tree),
                        tree->rates->c_lnL_times,
                        tree->rates->c_lnL_rates,
                        tree->rates->birth_rate,
                        tree->rates->death_rate,
                        mean_r,
                        tree->rates->nd_t[tree->n_root->num],
                        tree->next->mod->kappa->v,
                        tree->rates->nu);
          
          for(i=0;i<tree->mod->ras->n_catg;i++) PhyML_Fprintf(fp_stats,"%G\t",tree->mod->ras->gamma_rr->v[i]);
          for(i=0;i<tree->mod->ras->n_catg;i++) PhyML_Fprintf(fp_stats,"%G\t",tree->mod->ras->gamma_r_proba->v[i]);


          for(i=0;i<tree->rates->n_cal;i++)
            {
              t_cal *cal = tree->rates->a_cal[i];
              for(j=0;j<cal->clade_list_size;++j)
                {
                  t_clad *clade = cal->clade_list[j];
                  PhyML_Fprintf(fp_stats,"%G\t",tree->rates->nd_t[clade->target_nd->num]);
                }
            }
          
          for(i=0;i<tree->rates->n_cal;++i)
            {
              t_cal *cal = tree->rates->a_cal[i];
              /* PhyML_Fprintf(fp_stats,"%d\t",cal->current_clade_idx); */
              PhyML_Fprintf(fp_stats,"%s\t",cal->clade_list[cal->current_clade_idx]->id);
            }


          for(i=0;i<2*tree->n_otu-2;++i) PhyML_Fprintf(fp_stats,"%G\t",tree->rates->br_r[i]);

          for(i=0;i<2*tree->n_otu-1;++i) if(tree->a_nodes[i]->tax == NO) PhyML_Fprintf(fp_stats,"%G\t",tree->rates->nd_t[i]);

          PhyML_Fprintf(fp_stats,"%G\t",tree->mcmc->acc_rate[tree->mcmc->num_move_times_and_rates_root]);
          PhyML_Fprintf(fp_stats,"%G\t",tree->mcmc->tune_move[tree->mcmc->num_move_times_and_rates_root]);

          PhyML_Fprintf(fp_stats,"%G\t",tree->mcmc->acc_rate[tree->mcmc->num_move_root_time]);
          PhyML_Fprintf(fp_stats,"%G\t",tree->mcmc->tune_move[tree->mcmc->num_move_root_time]);

          PhyML_Fprintf(fp_stats,"%G\t",tree->mcmc->acc_rate[tree->mcmc->num_move_clock_r]);
          PhyML_Fprintf(fp_stats,"%G\t",tree->mcmc->tune_move[tree->mcmc->num_move_clock_r]);

          PhyML_Fprintf(fp_stats,"%G\t",tree->mcmc->acc_rate[tree->mcmc->num_move_updown_t_cr]);
          PhyML_Fprintf(fp_stats,"%G\t",tree->mcmc->tune_move[tree->mcmc->num_move_updown_t_cr]);

          PhyML_Fprintf(fp_stats,"%G\t",tree->mcmc->acc_rate[tree->mcmc->num_move_tree_rates]);
          PhyML_Fprintf(fp_stats,"%G\t",tree->mcmc->tune_move[tree->mcmc->num_move_tree_rates]);

          PhyML_Fprintf(fp_stats,"%G\t",tree->mcmc->acc_rate[tree->mcmc->num_move_spr_weighted]);
          PhyML_Fprintf(fp_stats,"%G\t",tree->mcmc->tune_move[tree->mcmc->num_move_spr_weighted]);

          PhyML_Fprintf(fp_stats,"%G\t",tree->mcmc->acc_rate[tree->mcmc->num_move_spr]);
          PhyML_Fprintf(fp_stats,"%G\t",tree->mcmc->tune_move[tree->mcmc->num_move_spr]);

          PhyML_Fprintf(fp_stats,"%G\t",tree->mcmc->acc_rate[tree->mcmc->num_move_spr_local]);
          PhyML_Fprintf(fp_stats,"%G\t",tree->mcmc->tune_move[tree->mcmc->num_move_spr_local]);
          
          if(tree->mcmc->sample_num == 0)
            {
              PhyML_Fprintf(fp_tree,"\n#NEXUS");
              PhyML_Fprintf(fp_tree,"\nBEGIN TREES;");
            }
          else
            {
              fseek(fp_tree,-11,SEEK_CUR);
            }
          
          
          Time_To_Branch(tree);
          tree->bl_ndigits = 3;
          /* RATES_Update_Cur_Bl(tree); */
          s_tree = Write_Tree(tree,NO);
          tree->bl_ndigits = 7;
          PhyML_Fprintf(fp_tree,"\ntree %f [&lnP=-1] = [&R] %s",s_tree,tree->c_lnL);          
          Free(s_tree);
          PhyML_Fprintf(fp_tree,"\nEND;");          
          fflush(NULL);
          RATES_Update_Cur_Bl(tree);

          tree->mcmc->sample_num++;

        }

      tree->mcmc->run++;
      MCMC_Get_Acc_Rates(tree->mcmc);
    }
  while(tree->mcmc->run < tree->mcmc->chain_len);


  return(res);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// Update the list of nodes that are younger than lim
void DATE_List_Of_Nodes_Younger_Than(t_node *a, t_node *d, phydbl lim, t_ll **list, t_tree *tree)
{
  if(tree->rates->nd_t[d->num] > lim) Push_Bottom_Linked_List(d,list,YES);
  
  if(d->tax == YES) return;
  else
    {
      int i;

      if(d == tree->n_root)
        {
          DATE_List_Of_Nodes_Younger_Than(d,d->v[1],lim,list,tree);
          DATE_List_Of_Nodes_Younger_Than(d,d->v[2],lim,list,tree);
        }
      else
        {
          for(i=0;i<3;i++)
            if(d->v[i] != a && d->b[i] != tree->e_root)
              DATE_List_Of_Nodes_Younger_Than(d,d->v[i],lim,list,tree);            
        }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// Update the list of nodes that are younger than lim with direct ancestors 
// not younger than lim
void DATE_List_Of_Nodes_And_Ancestors_Younger_Than(t_node *a, t_node *d, phydbl lim, t_ll **list, t_tree *tree)
{
  if(tree->rates->nd_t[d->num] > lim && a != NULL && tree->rates->nd_t[a->num] > lim) Push_Bottom_Linked_List(d,list,YES);
  
  if(d->tax == YES) return;
  else
    {
      int i;
      
      if(d == tree->n_root)
        {
          DATE_List_Of_Nodes_And_Ancestors_Younger_Than(d,d->v[1],lim,list,tree);
          DATE_List_Of_Nodes_And_Ancestors_Younger_Than(d,d->v[2],lim,list,tree);
        }
      else
        {      
          for(i=0;i<3;i++)
            if(d->v[i] != a && d->b[i] != tree->e_root)
              DATE_List_Of_Nodes_And_Ancestors_Younger_Than(d,d->v[i],lim,list,tree);            
        }
    }
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// List of valid regraft nodes, taking into account calibration
// constraints. The subtree (defined by prune and prune_daughter
// will be re-attached *on top of* one of the nodes in this list
// (as opposed to *on one of the sister edges below*).
t_ll *DATE_List_Of_Regraft_Nodes(t_node *prune, t_node *prune_daughter, phydbl *t_min, phydbl *t_max, int verbose, t_tree *tree)
{
  t_node *n,*m;
  int i,j;
  t_ll *out,*in,*ll;
  int is_clade_affected;
  t_clad *clade;
  t_cal *cal;

  cal = NULL;
  clade = NULL;
  n = NULL;
  m = NULL;
  *t_min = -INFINITY;
  in = NULL;
  out = NULL;
  is_clade_affected = NO;


  // Find the oldest LCA of calibrated sets among the nodes between
  // prune and root. These clades might see the position of their
  // LCA change.

  if(prune != tree->n_root)
    {
      n = prune;
      while(n)
        {
          for(i=0;i<tree->rates->n_cal;++i)
            {
              // That node is the LCA of calibration a_cal[i]
              cal = tree->rates->a_cal[i];
              clade = cal->clade_list[cal->current_clade_idx];
              
              if(n == clade->target_nd)
                {
                  is_clade_affected = NO;

                  for(j=0;j<clade->n_tax;++j)
                    {
                      m = clade->tip_list[j];
                      do
                        {
                          if(m == prune_daughter)
                            {
                              is_clade_affected = YES;
                              break;
                            }
                          m = m->anc;
                        }
                      while(m);
                      
                      if(is_clade_affected == YES) break;
                    }

                  
                  // Maximum of the lower bounds for calibration intervals
                  /* if(is_clade_affected == YES) *t_min = MAX(*t_min,tree->rates->a_cal[i]->lower); */
                  if(is_clade_affected == YES) *t_min = MAX(*t_min,tree->rates->t_prior_min[n->num]);
                }
            }
          n = n->anc;
        }
    }

  // Find the oldest internal node within intervals defined by
  // calibrations affected by the pruning.
  n = prune_daughter;
  while(n->anc && !(tree->rates->nd_t[n->anc->num] < *t_min))
    {
      n = n->anc;
      assert(n);
    }
  

  if(verbose) printf("\n. Apical: %d @ time %f min: %f",n->num,tree->rates->nd_t[n->num],*t_min); fflush(NULL);
     
  // List all nodes younger than this apical node
  DATE_List_Of_Nodes_Younger_Than(n->anc,n,-INFINITY,&in,tree);

  assert(in != NULL);
  
  if(verbose)
    {
      ll = in->head;
      t_node *x;
      do
        {
          x = (t_node *)ll->v;
          PhyML_Printf("\nx Inlist %d @ %f",x->num,tree->rates->nd_t[x->num]);
          ll = ll->next;
        }
      while(ll != NULL);
    }
  
  // Remove from that list the nodes that are too young to be suitable regraft points.
  
  
  n = prune_daughter;
  out = NULL;
  while(n)
    {
      for(i=0;i<tree->rates->n_cal;i++)
        {
          cal = tree->rates->a_cal[i];
          clade = cal->clade_list[cal->current_clade_idx];

          if(n->anc && n->anc == clade->target_nd)
            {
              for(j=0;j<clade->n_tax;++j)
                {
                  m = clade->tip_list[j];
                  do
                    {
                      if(m == prune_daughter) break;
                      if(m == prune) break;
                      m = m->anc;
                    }
                  while(m);
                  if(m == prune) break; // Prune-regraft anywhere below calibrated node will not change that node.
                }
              
              if(m != prune)
                {
                  for(j=0;j<3;++j)
                    {
                      if(n->anc->v[j] != n->anc->anc && n->anc->b[j] != tree->e_root && n->anc->v[j] != n)
                        {
                          DATE_List_Of_Nodes_And_Ancestors_Younger_Than(n->anc,
                                                                        n->anc->v[j],
                                                                        tree->rates->a_cal[i]->upper,
                                                                        &out,
                                                                        tree);
                          break;
                        }
                    }
                }
            }
        }
      n = n->anc;
    }

  // Remove nodes that are `strictly' younger than prune_daughter
  DATE_List_Of_Nodes_And_Ancestors_Younger_Than(tree->n_root,tree->n_root->v[1],tree->rates->nd_t[prune_daughter->num],&out,tree);
  DATE_List_Of_Nodes_And_Ancestors_Younger_Than(tree->n_root,tree->n_root->v[2],tree->rates->nd_t[prune_daughter->num],&out,tree);

  // Remove nodes that are below prune_daughter (prune_daughter included)
  DATE_List_Of_Nodes_Younger_Than(prune,prune_daughter,-INFINITY,&out,tree);

  // Add prune node to the list of node that can't be targeted for regraft
  Push_Bottom_Linked_List(prune,&out,YES);

  // Add root node as one cannot regraft above it
  /* Push_Bottom_Linked_List(tree->n_root,&out); */

  if(verbose)
    {
      printf("\nx outlist: %p",(void *)out); fflush(NULL);
      ll = out->head;
      do
        {
          t_node *x = (t_node *)ll->v;
          PhyML_Printf("\nx Outlist %d @ %f",x->num,tree->rates->nd_t[x->num]);
          ll = ll->next;
        }
      while(ll != NULL);
    }
  

  /* Print_List(in); */
  ll = out->head;
  do
    {
      if(verbose)
        {
          t_node *x = (t_node *)ll->v;
          printf("\nx Remove %d",x->num);
        }

      Remove_From_Linked_List(NULL,ll->v,&in);

      if(verbose) PhyML_Printf("\n. List in (in->head:%p in->tail:%p):",in?in->head:NULL,in?in->tail:NULL);

      /* Print_List(in); */
      ll = ll->next;
    }
  while(ll != NULL);
  
  Free_Linked_List(out);

  if(verbose)
    {
      ll = in->head;
      do
        {
          t_node *x;
          x = (t_node *)ll->v;
          printf("\n. In1: %d",x->num); fflush(NULL);
          ll = ll->next;
        }
      while(ll != NULL);
    }
  
  return(in);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl DATE_Lk_Calib(t_tree *tree)
{
  phydbl lnL;
  int i;
  t_cal *cal;

  lnL = 0.0;
  for(i=0;i<tree->rates->n_cal;++i)
    {
      cal = tree->rates->a_cal[i];
      lnL += LOG(cal->alpha_proba_list[cal->current_clade_idx]);
    }
  
  return lnL;
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

