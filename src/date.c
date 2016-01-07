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

void DATE_XML(char *xml_filename)
{
  FILE *fp_xml_in;
  xml_node *xnd,*xnd_dum,*xnd_cal,*xroot;
  t_tree *mixt_tree;
  phydbl low,up;
  int node_num;
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
                  
                  if(!strcmp("root",clade_name))
                    {
                      node_num = mixt_tree->n_root->num;
                    }
                  else
                    {
                      xml_node *xnd_clade;
                        
                      node_num = -1;
                      xnd_clade = XML_Search_Node_Generic("clade","id",clade_name,YES,xroot);
                      
                      if(xnd_clade != NULL) // found clade with a given name
                        {
                          char **clade;
                          int clade_size;
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


                          PhyML_Printf("\n. Node number to which calibration [%s] applies to is [%d]",
                                       clade_name,
                                       Find_Clade(clade,
                                                  clade_size,
                                                  mixt_tree));
                          
                          PhyML_Printf("\n. Lower bound set to: %15f time units.", low);
                          PhyML_Printf("\n. Upper bound set to: %15f time units.", up);
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

  fflush(NULL);

  DATE_Assign_Primary_Calibration(mixt_tree);
  DATE_Update_T_Prior_MinMax(mixt_tree);

  Update_Ancestors(mixt_tree->n_root,mixt_tree->n_root->v[2],mixt_tree);
  Update_Ancestors(mixt_tree->n_root,mixt_tree->n_root->v[1],mixt_tree);		  
  MCMC_Randomize_Node_Times(mixt_tree);
  

  PhyML_Printf("\n. K=%f",DATE_J_Sum_Product(mixt_tree));

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
          if(tree->a_nodes[i]->cal != NULL) // Primary calibration found on that node
            {
              tree->rates->t_prior_max[i] =     tree->a_nodes[i]->cal->upper;
              tree->rates->t_prior_min[i] = MAX(tree->a_nodes[i]->cal->lower,tree->rates->nd_t[tree->n_root->num]);
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


  For(i,tree->n_otu-1)
    {
      PhyML_Printf("\n. Node %d age: %f t_prior_min: %f t_prior_max: %f",
                   rk[i],
                   tree->rates->nd_t[rk[i]],
                   tree->rates->t_prior_min[rk[i]],
                   tree->rates->t_prior_max[rk[i]]);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void DATE_Assign_Primary_Calibration(t_tree *tree)
{
  int i,node_num;
  
  For(i,2*tree->n_otu-1) tree->a_nodes[i]->cal = NULL;

  For(i,tree->rates->n_cal)
    {
      node_num = Find_Clade(tree->rates->a_cal[i]->target_tax,
                            tree->rates->a_cal[i]->n_target_tax,
                            tree);

      tree->a_nodes[node_num]->cal = tree->rates->a_cal[i];
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// Return splitted calibration intervals. Make sure all primary
// and secondary calibration intervals are up-to-date.
phydbl *DATE_Splitted_Calibration(t_tree *tree)
{
  phydbl *splitted_cal,buff;
  int i,len,done;

  // One t_prior_min and one t_prior_max per internal nodes, so 
  // 2 x # of internal nodes boundaries in total at most.
  splitted_cal = (phydbl *)mCalloc(2*(tree->n_otu-1),sizeof(phydbl));
  For(i,2*(tree->n_otu-1)) splitted_cal[i] = -INFINITY;

  len = 0;
  for(i = tree->n_otu; i < 2*tree->n_otu-1; i++)
    {
      splitted_cal[len]   = tree->rates->t_prior_min[i];
      splitted_cal[len+1] = tree->rates->t_prior_max[i];        
      len++;
    }

  
  // Bubble sort of all these times in increasing order
  done = YES;
  do
    {
      done = YES;
      For(i,len-1) 
        {
          if(splitted_cal[i] > splitted_cal[i+1])
            {
              buff              = splitted_cal[i];
              splitted_cal[i]   = splitted_cal[i+1];
              splitted_cal[i+1] = buff;
              done = NO;
            }
        }
    }
  while(done == NO);

  For(i,len-1) assert(!(splitted_cal[len] > splitted_cal[len+1]));

  return splitted_cal;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl DATE_J_Sum_Product(t_tree *tree)
{
  phydbl prod, total,*splitted_cal;
  int fact,idx,ans;
  
  DATE_Assign_Primary_Calibration(tree);
  DATE_Update_T_Prior_MinMax(tree);

  splitted_cal = DATE_Splitted_Calibration(tree);
  
  ans   = 0;
  total = 0.0;
  idx   = 0;
  do
    {
      prod = 0.0;
      fact = 0;
      ans = DATE_J_Sum_Product_Pre(tree->n_root,
                                   idx,
                                   -1,
                                   &prod,&fact,&total,splitted_cal,tree);
      idx++;
      if(ans == 1) break;
    }
  while(ans <= 0);

  return(total);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int DATE_J_Sum_Product_Pre(t_node *d, int split_idx_d, int split_idx_a, phydbl *prod, int *fact, phydbl *total, phydbl *splitted_cal, t_tree *tree)
{
  int ans,idx;

  ans = DATE_Is_Split_Accessible(d,split_idx_d,splitted_cal,tree);
  if(ans != 0) return ans;

  // Calulate J for this time interval
  (*prod) *= DATE_J(tree->rates->birth_rate, 
                    tree->rates->death_rate,
                    splitted_cal[split_idx_d],
                    splitted_cal[split_idx_d+1]);
  
  // Remove factorial term from current product
  (*prod) *= (*fact);

  if(split_idx_d == split_idx_a) (*fact)++;
  else *fact = 1;

  (*prod) /= (*fact);

  if(d->tax) 
    {
      (*total) += (*prod);
      return 0;
    }
  else
    {
      idx = split_idx_d;
      do
        {
          ans = DATE_J_Sum_Product_Pre(tree->a_nodes[tree->rates->t_rank[d->num+1]],
                                       idx,
                                       split_idx_d,
                                       prod,fact,total,splitted_cal,tree);
          idx++;
          if(ans == 1) return ans;
        }
      while(ans <= 0);
    }
  assert(FALSE);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int DATE_Is_Split_Accessible(t_node *d, int which, phydbl *splitted_cal, t_tree *tree)
{
  if(tree->rates->t_prior_min[d->num] < splitted_cal[which] &&
     tree->rates->t_prior_max[d->num] > splitted_cal[which+1]) return 0; // splitted interval is within [t_prior_min,t_prior_max]
  else if(splitted_cal[which] > tree->rates->t_prior_max[d->num]) return 1; // splitted interval is younger than [t_prior_min,t_prior_max]
  else if(splitted_cal[which+1] < tree->rates->t_prior_min[d->num]) return -1; // splitted interval is older than [t_prior_min,t_prior_max]
  else assert(FALSE); // splitted interval cannot be partially overlapping [t_prior_min,t_prior_max]
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl DATE_J(phydbl birth_r, phydbl death_r, phydbl t_min, phydbl t_pls)
{
  phydbl d,b,J;
  assert(t_pls > t_min);
  d = death_r;
  b = birth_r;
  J = (d-b)*(EXP(t_pls*d + t_min*b) - EXP(t_pls*b + t_min*d));
  J /= (b*EXP(t_min*b) - d*EXP(t_min*d)) * (b*EXP(t_pls*b) - d*EXP(t_pls*d));
  return(J);
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

