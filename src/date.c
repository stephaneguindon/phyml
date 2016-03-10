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
  int seed;

  seed = getpid();

  /* seed = 1; */
  /* seed = 8596; */
  /* seed = 23595; */
  /* seed = 10868; */

  printf("\n. seed: %d",seed);
  srand(seed);

  io = Get_Input(argc,argv);
  Free(io);
  return(0);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void DATE_XML(char *xml_filename)
{
  FILE *fp_xml_in;
  xml_node *xnd,*xnd_dum,*xnd_cal,*xroot;
  t_tree *mixt_tree;
  phydbl low,up,*res;
  char *clade_name;

  mixt_tree = XML_Process_Base(xml_filename);
  assert(mixt_tree);

  mixt_tree->rates = RATES_Make_Rate_Struct(mixt_tree->n_otu);
  RATES_Init_Rate_Struct(mixt_tree->rates,NULL,mixt_tree->n_otu);

  fp_xml_in = fopen(xml_filename,"r");
  if(!fp_xml_in)
    {
      PhyML_Printf("\n== Could not find the XML file '%s'.\n",xml_filename);
      Exit("\n");
    }

  xroot = XML_Load_File(fp_xml_in);

  if(xroot == NULL)
    {
      PhyML_Printf("\n== Encountered an issue while loading the XML file.\n");
      Exit("\n");
    }
  
  // Looking for calibration node(s)
  xnd = XML_Search_Node_Name("calibration",YES,xroot);

  if(xnd == NULL)
    {
      PhyML_Printf("\n== No calibration information seems to be provided.");
      PhyML_Printf("\n== Please amend your XML file. \n");
      Exit("\n");
    }
  else
    {
      if(XML_Search_Node_Name("upper",NO,xnd->child) == NULL && XML_Search_Node_Name("lower",NO,xnd->child) == NULL)
	{
	  PhyML_Printf("\n== There is no calibration information provided. \n");
	  PhyML_Printf("\n== Please check your data. \n");
	  Exit("\n");
	}
    }


  MIXT_Prepare_All(-1,mixt_tree);
  if(mixt_tree->n_root == NULL) Add_Root(mixt_tree->a_edges[0],mixt_tree);

  clade_name = (char *)mCalloc(T_MAX_NAME,sizeof(char));

  xnd = xroot->child;
  assert(xnd);
  do
    {
      if(!strcmp(xnd->name,"calibration")) // Found a XML node <calibration>.
	{
          xnd_cal = xnd;

          // TO DO: make sure calibs are shared across partition elements -> need to write chain function to
          // call once the calib struct on the first mixt_tree is initialized.
          /* mixt_tree->rates->tot_num_cal++; */
	  /* if (mixt_tree->rates->calib == NULL) mixt_tree->rates->calib = Make_Calib(mixt_tree->n_otu); */

	  low = -BIG;
	  up  = BIG;

	  xnd_dum = XML_Search_Node_Name("lower",YES,xnd_cal);
	  if(xnd_dum != NULL) low = String_To_Dbl(xnd_dum->value); 

	  xnd_dum = XML_Search_Node_Name("upper",YES,xnd_cal);
	  if(xnd_dum != NULL) up = String_To_Dbl(xnd_dum->value);
          
          do
            {
              if(!strcmp("appliesto",xnd_cal->child->name)) 
                {
                  clade_name = XML_Get_Attribute_Value(xnd_cal->child,"clade.id");
                  
                  if(!clade_name)
                    {
                      PhyML_Printf("\n== Attribute 'value=CLADE_NAME' is mandatory");
                      PhyML_Printf("\n== Please amend your XML file accordingly.");
                      Exit("\n");
                    }
                  
                  if(strcmp("root",clade_name))
                    {
                      xml_node *xnd_clade;
                        
                      xnd_clade = XML_Search_Node_Generic("clade","id",clade_name,YES,xroot);
                      
                      if(xnd_clade != NULL) // found clade with a given name
                        {
                          char **clade;
                          int clade_size,nd_num;
                          t_cal *cal;

                          clade      = XML_Read_Clade(xnd_clade->child,mixt_tree);
                          clade_size = XML_Number_Of_Taxa_In_Clade(xnd_clade->child);
                          // TO DO: chain all calibrations
                          cal        = Make_Calibration();
                          
                          Init_Calibration(cal);

                          mixt_tree->rates->a_cal[mixt_tree->rates->n_cal] = cal;
                          mixt_tree->rates->n_cal++;

                          cal->is_primary   = YES;
                          cal->target_tax   = clade;
                          cal->n_target_tax = clade_size;
                          cal->lower        = low;
                          cal->upper        = up;

                          nd_num = Find_Clade(clade,clade_size,mixt_tree);

                          PhyML_Printf("\n. Node number to which calibration [%s] applies to is [%d]",clade_name,nd_num);                          
                          PhyML_Printf("\n. Lower bound set to: %15f time units.",low);
                          PhyML_Printf("\n. Upper bound set to: %15f time units.",up);
                          PhyML_Printf("\n. .......................................................................");
                        }
                      else
                        {
                          PhyML_Printf("\n== Calibration information for clade [%s] was not found.", clade_name);
                          PhyML_Printf("\n== Err. in file %s at line %d\n",__FILE__,__LINE__);
                          Exit("\n");
                        }                      
                    }
                }
              xnd_cal->child = xnd_cal->child->next;
            }
          while(xnd_cal->child != NULL);
        }
          /* mixt_tree->rates->calib = mixt_tree->rates->calib->next;	    */
      xnd = xnd->next;
    }
  while(xnd != NULL);

  DATE_Assign_Primary_Calibration(mixt_tree);
  DATE_Update_T_Prior_MinMax(mixt_tree);
  DATE_Chain_Cal(mixt_tree);

  Update_Ancestors(mixt_tree->n_root,mixt_tree->n_root->v[2],mixt_tree);
  Update_Ancestors(mixt_tree->n_root,mixt_tree->n_root->v[1],mixt_tree);		  
 
  TIMES_Randomize_Tree_With_Time_Constraints(mixt_tree->rates->a_cal[0], mixt_tree);
  
  /* { */
  /*   Print_Node(mixt_tree->n_root,mixt_tree->n_root->v[1],mixt_tree); */
  /*   Print_Node(mixt_tree->n_root,mixt_tree->n_root->v[2],mixt_tree); */
  /* } */

  mixt_tree->rates->birth_rate = 0.10;
  mixt_tree->rates->death_rate = 0.05;
  mixt_tree->rates->bl_from_rt = YES;
  mixt_tree->rates->clock_r    = 0.01 / FABS(mixt_tree->rates->nd_t[mixt_tree->n_root->num]);
  mixt_tree->rates->model      = BIRTHDEATH;

  RATES_Update_Cur_Bl(mixt_tree);

  res = DATE_MCMC(mixt_tree);

  Free(res);
  Free(clade_name);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// Update t_prior_min and t_prior_max on a given ranked tree
// given (primary and secondary) calibration information.
// Make sure secondary and primary calibration are up-to-date
void DATE_Update_T_Prior_MinMax(t_tree *tree)
{
  int i,j,*rk;

  rk = tree->rates->t_rank;

  for(i=tree->n_otu;i<2*tree->n_otu-1;i++) // All internal nodes except the root 
    {
      if(tree->a_nodes[i] != tree->n_root)
        {          
          if(tree->a_nodes[i]->n_cal > 0) // Primary calibration found on that node
            {
              tree->rates->t_prior_max[i] = 0.0;
              tree->rates->t_prior_min[i] = -INFINITY;
              For(j,tree->a_nodes[i]->n_cal)
                {
                  tree->rates->t_prior_max[i] = MIN(tree->rates->t_prior_max[i],tree->a_nodes[i]->cal[j]->upper);
                  tree->rates->t_prior_min[i] = MAX(tree->rates->t_prior_min[i],MAX(tree->a_nodes[i]->cal[j]->lower,tree->rates->nd_t[tree->n_root->num]));
                }
            }
          else
            {
              tree->rates->t_prior_max[i] = 0.0;
              tree->rates->t_prior_min[i] = tree->rates->nd_t[tree->n_root->num];
            }
        }
    }

  // TO DO: chain t_rank
  TIMES_Update_Node_Ordering(tree);

  For(i,tree->n_otu-2)
    {
      for(j=i+1;j<tree->n_otu-1;j++)
        {
          if(tree->rates->t_prior_min[rk[j]] < tree->rates->t_prior_min[rk[i]])
            tree->rates->t_prior_min[rk[j]] = tree->rates->t_prior_min[rk[i]];

          if(tree->rates->t_prior_max[rk[i]] > tree->rates->t_prior_max[rk[j]])
            tree->rates->t_prior_max[rk[i]] = tree->rates->t_prior_max[rk[j]];
        }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void DATE_Assign_Primary_Calibration(t_tree *tree)
{
  int i,j,idx,node_num;
  
  For(i,2*tree->n_otu-1) 
    For(j,MAX_N_CAL) 
    {
      tree->a_nodes[i]->cal[j] = NULL;
      tree->a_nodes[i]->n_cal  = 0;
    }

  For(i,tree->rates->n_cal)
    {
      node_num = Find_Clade(tree->rates->a_cal[i]->target_tax,
                            tree->rates->a_cal[i]->n_target_tax,
                            tree);
      
      idx = tree->a_nodes[node_num]->n_cal;
      tree->a_nodes[node_num]->cal[idx] = tree->rates->a_cal[i];
      tree->a_nodes[node_num]->n_cal++;

      if(tree->a_nodes[node_num]->n_cal == MAX_N_CAL)
        {
          PhyML_Printf("\n== A node cannot have more than %d calibration",MAX_N_CAL); 
          PhyML_Printf("\n== constraints attached to it. Feel free to increase the"); 
          PhyML_Printf("\n== value of the variable MAX_N_CAL in utilities.h if");
          PhyML_Printf("\n== necessary.");
          Exit("\n");
        }

      /* printf("\n. Assign cal [%f %f] to %d (%d)", */
      /*        tree->rates->a_cal[i]->lower, */
      /*        tree->rates->a_cal[i]->upper, */
      /*        node_num, */
      /*        tree->a_nodes[node_num]->n_cal); */
      /* int j; */
      /* For(j,tree->rates->a_cal[i]->n_target_tax) */
      /*   printf("\n> %s",tree->rates->a_cal[i]->target_tax[j]); */
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
      For(i,len-1) 
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

  For(i,len-1) assert(!(minmax[i] > minmax[i+1]));
  
  // Remove ties
  For(i,len-1) 
    if(Are_Equal(minmax[i],minmax[i+1],1.E-6) == YES) 
      minmax[i] = 0.0;

  // Sort again to effectively remove ties
  do
    {
      done = YES;
      For(i,len-1) 
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

  return(total);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int DATE_J_Sum_Product_Pre(t_node *d, int split_idx_d, int split_idx_a, phydbl prod, int fact, phydbl *total, phydbl *splitted_cal, int rk, t_tree *tree)
{
  int ans,idx;

  ans = DATE_Is_Split_Accessible(d,split_idx_d,splitted_cal,tree);

  /* printf("\n. IN d: %d [%12f %12f] %d ans: %d prod: %f fact: %d", */
  /*        d->num, */
  /*        splitted_cal[split_idx_d], */
  /*        splitted_cal[split_idx_d+1], */
  /*        tree->rates->t_rank[d->num],ans,prod,fact); */

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
        
        /* PhyML_Printf("\n. Node: %d [%12f %12f] %4d %4d [%12G %12G] %4d", */
        /*              d->num, */
        /*              splitted_cal[split_idx_d], */
        /*              splitted_cal[split_idx_d+1], */
        /*              split_idx_a, */
        /*              split_idx_d, */
        /*              prod,*total,fact); */
        
        /* fflush(NULL); */
        
        if(tree->rates->t_rank[tree->n_otu-2] == d->num) // Youngest internal node 
          {
            (*total) += prod;
            /* printf(" == total: %f",*total); */
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

  /* printf("\n. OUT d: %d [%12f %12f] %d ans: %d", */
  /*        d->num, */
  /*        splitted_cal[split_idx_d], */
  /*        splitted_cal[split_idx_d+1], */
  /*        tree->rates->t_rank[d->num],ans); */

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
      PhyML_Printf("\n== d->num: %d d->tax: %d",d->num,d->tax);
      PhyML_Printf("\n== t_prior_min: %f t_prior_max: %f",
                   tree->rates->t_prior_min[d->num],
                   tree->rates->t_prior_max[d->num]);
      PhyML_Printf("\n== splitted_cal_min: %f splitted_cal_max: %f",
                   splitted_cal[which],
                   splitted_cal[which+1]);      
      PhyML_Printf("\n");
      assert(FALSE); // splitted interval cannot be partially overlapping [t_prior_min,t_prior_max]
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl DATE_J(phydbl birth_r, phydbl death_r, phydbl t_min, phydbl t_pls)
{
  phydbl d,b,J;
  assert(t_pls > t_min);
  d = death_r;
  b = birth_r;
  J = (b-d)*(EXP(t_min*d+t_pls*b) - EXP(t_min*b+t_pls*d));
  J /= ((b*EXP(t_min*b)-d*EXP(t_min*d)) * (b*EXP(t_pls*b)-d*EXP(t_pls*d)));
  /* printf("  J : %f",J); */
  return(J);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void DATE_Chain_Cal(t_tree *mixt_tree)
{
  int i;
  For(i,mixt_tree->rates->n_cal-1) 
    {
      mixt_tree->rates->a_cal[i]->next   = mixt_tree->rates->a_cal[i+1];
      mixt_tree->rates->a_cal[i+1]->prev = mixt_tree->rates->a_cal[i];
    }
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
  int move, n_vars, i, adjust_len;
  phydbl u;
  phydbl *res;
  FILE *fp_stats;


  fp_stats = tree->io->fp_out_stats;

  mcmc = MCMC_Make_MCMC_Struct();

  tree->mcmc = mcmc;

  MCMC_Init_MCMC_Struct(NULL,NULL,mcmc);

  MCMC_Complete_MCMC(mcmc,tree);

  n_vars      = 10;
  adjust_len  = 1E+6;

  res = (phydbl *)mCalloc(tree->mcmc->chain_len / tree->mcmc->sample_interval * n_vars,sizeof(phydbl));
  
  Lk(NULL,tree);
  TIMES_Lk_Birth_Death(tree);
  
  PhyML_Printf("\n. log(Pr(Seq|Tree)) = %f",tree->c_lnL);
  PhyML_Printf("\n. log(Pr(Tree)) = %f",tree->rates->c_lnL_times);
  

  For(i,mcmc->n_moves) tree->mcmc->start_ess[i] = YES;

  Set_Both_Sides(NO,tree);
  mcmc->use_data   = YES; 
  mcmc->always_yes = NO;
  move             = -1;
  do
    {

      /* tree->mcmc->adjust_tuning[i] = NO; */
      if(mcmc->run > adjust_len) For(i,mcmc->n_moves) tree->mcmc->adjust_tuning[i] = NO;

      if(tree->c_lnL < UNLIKELY + 0.1)
        {
          PhyML_Printf("\n== Move '%s' failed\n",tree->mcmc->move_name[move]);
          assert(FALSE);
        }

      u = Uni();

      For(move,tree->mcmc->n_moves) if(tree->mcmc->move_weight[move] > u-1.E-10) break;

      assert(!(move == tree->mcmc->n_moves));
      
      if(!strcmp(tree->mcmc->move_name[move],"birth_rate"))
        MCMC_Birth_Rate(tree);

      if(!strcmp(tree->mcmc->move_name[move],"death_rate"))
        MCMC_Death_Rate(tree);
      
      tree->mcmc->run++;
      MCMC_Get_Acc_Rates(tree->mcmc);

      if(!(tree->mcmc->run%tree->mcmc->sample_interval))
        {
          PhyML_Fprintf(fp_stats,"\n%6d\t%9.1f\t%9.1f\t%10G\t%10G",
                        tree->mcmc->run,
                        tree->c_lnL,
                        tree->rates->c_lnL_times,
                        tree->rates->birth_rate,
                        tree->rates->death_rate);
          tree->mcmc->sample_num++;
        }
    }
  while(tree->mcmc->run < tree->mcmc->chain_len);


  return(res);
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

