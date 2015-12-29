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

                          clade             = XML_Read_Clade(xnd_clade->child,mixt_tree);
                          clade_size        = XML_Number_Of_Taxa_In_Clade(xnd_clade->child);
                          // TO DO: chain all calibrations
                          cal               = Make_Calibration();

                          mixt_tree->rates->a_cal[mixt_tree->rates->n_cal] = cal;
                          mixt_tree->rates->n_cal++;

                          cal->is_primary                   = YES;
                          cal->target_tax                   = clade;
                          cal->n_target_tax                 = clade_size;
                        }
                      else
                        {
                          PhyML_Printf("\n== Calibration information for clade [%s] was not found.", clade_name);
                          PhyML_Printf("\n== Err. in file %s at line %d\n",__FILE__,__LINE__);
                          Exit("\n");
                        }
                      
                      PhyML_Printf("\n. Node number to which calibration [%s] applies to is [%d]",clade_name,node_num);
                      PhyML_Printf("\n. Lower bound set to: %15f time units.", low);
                      PhyML_Printf("\n. Upper bound set to: %15f time units.", up);
                      PhyML_Printf("\n. .......................................................................");
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
  DATE_Update_Secondary_Cal(mixt_tree);
  DATE_Update_T_Prior_MinMax(mixt_tree);

  Free(clade_name);
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// Update secondary calibration on all internal nodes. Requires two tree traversals.
void DATE_Update_Secondary_Cal(t_tree *tree)
{
  DATE_Update_Secondary_Cal_Post(tree->n_root,tree->n_root->v[1],tree);
  DATE_Update_Secondary_Cal_Post(tree->n_root,tree->n_root->v[2],tree);
  DATE_Update_Secondary_Cal_Pre(tree->n_root,tree->n_root->v[1],tree);
  DATE_Update_Secondary_Cal_Pre(tree->n_root,tree->n_root->v[2],tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// Updates upper bounds
void DATE_Update_Secondary_Cal_Post(t_node *a, t_node *d, t_tree *tree)
{
  if(d->tax == YES) return;
  else
    {
      int i;
      phydbl max;

      For(i,3)
        {
          if(d->v[i] != a && d->b[i] != tree->e_root)
            {
              DATE_Update_Secondary_Cal_Post(d,d->v[i],tree);
            }
        }

      if(d->cal->is_primary == NO)
        {
          max = -DBL_MIN;
          For(i,3) 
            if(d->v[i] != a && d->v[i]->cal->upper > max) 
              max = d->v[i]->cal->upper;
          d->cal->upper = max;
        }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// Updates lower bounds
void DATE_Update_Secondary_Cal_Pre(t_node *a, t_node *d, t_tree *tree)
{
  if(d->tax == YES) return;
  else
    {
      int i;
      phydbl t_root;

      t_root = tree->rates->nd_t[tree->n_root->num];
      
      if(d->cal->is_primary == NO)
        For(i,3) 
          if(d->v[i] == a) 
            d->cal->lower = MAX(t_root,a->cal->lower);

      For(i,3)
        {
          if(d->v[i] != a && d->b[i] != tree->e_root)
            {
              DATE_Update_Secondary_Cal_Pre(d,d->v[i],tree);
            }
        }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// Update t_prior_min and t_prior_max on a given ranked tree
// given (primary and secondary) calibration information.
// Make sure secondary and primary calibration is up-to-date
void DATE_Update_T_Prior_MinMax(t_tree *tree)
{
  int i,j;

  for(i=tree->n_otu;i<2*tree->n_otu-1;i++)
    {
      tree->rates->t_prior_max[i] = tree->a_nodes[i]->cal->upper;;
      tree->rates->t_prior_min[i] = tree->a_nodes[i]->cal->lower;
    }

  TIMES_Update_Node_Ordering(tree);

  for(i=tree->n_otu;i<2*tree->n_otu-2;i++)
    {
      for(j=i+1;j<2*tree->n_otu-1;j++)
        {
          if(tree->rates->t_prior_min[j] < tree->rates->t_prior_min[i])
            tree->rates->t_prior_min[j] = tree->rates->t_prior_min[i];

          if(tree->rates->t_prior_max[i] > tree->rates->t_prior_max[j])
            tree->rates->t_prior_max[i] = tree->rates->t_prior_max[j];
        }
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

