/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

/* Routines that implement Etheridge and Barton's model of continuous-space
   coalescent.
*/

#include "phyrex.h"

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int PHYREX_Main(int argc, char *argv[])
{
#if (defined PHYREXSIM)
  PHYREX_Main_Simulate(argc,argv);
#elif (defined PHYREX)
  option *io;
  io = Get_Input(argc,argv);
  Free(io);
#endif
  return(0);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

#if (defined PHYREX)
void PHYREX_XML(char *xml_filename)
{
  FILE *fp_xml_in;
  xml_node *xnd,*xroot;
  t_tree *mixt_tree,*tree;
  phydbl *res;
  int seed;
  char *dum_string;
 
  mixt_tree = XML_Process_Base(xml_filename);
  assert(mixt_tree);

  mixt_tree->rates = RATES_Make_Rate_Struct(mixt_tree->n_otu);
  RATES_Init_Rate_Struct(mixt_tree->rates,NULL,mixt_tree->n_otu);
  
  mixt_tree->times = TIMES_Make_Time_Struct(mixt_tree->n_otu);
  TIMES_Init_Time_Struct(mixt_tree->times,NULL,mixt_tree->n_otu);

  mixt_tree->mmod = PHYREX_Make_Migrep_Model(mixt_tree->n_otu,2);
  mixt_tree->mmod->n_dim = 2;
  PHYREX_Set_Default_Migrep_Mod(mixt_tree->n_otu,mixt_tree->mmod);
  
  tree = mixt_tree;
  do
    {
      // All rate stuctures point to the same object
      tree->rates = mixt_tree->rates;
      tree->times = mixt_tree->times;
      tree = tree->next;
    }
  while(tree);

  
  fp_xml_in = fopen(xml_filename,"r");
  if(!fp_xml_in)
    {
      PhyML_Fprintf(stderr,"\n. Could not find the XML file '%s'.\n",xml_filename);
      Exit("\n");
    }

  /* xroot = XML_Load_File(fp_xml_in); */
  xroot = mixt_tree->xml_root;

  if(xroot == NULL)
    {
      PhyML_Fprintf(stderr,"\n. Encountered an issue while loading the XML file.\n");
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
    }

  xnd = XML_Search_Node_Name("phyrex",NO,xroot);

  if(xnd == NULL)
    {
      PhyML_Fprintf(stderr,"\n. Cound not find the \"root\" of the XML file (it should have \'phyrex\' as tag name).\n");
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
    }

  dum_string = XML_Get_Attribute_Value(xnd,"mcmc.chain.len");
  if(dum_string != NULL) mixt_tree->io->mcmc->chain_len = (int)String_To_Dbl(dum_string);
  
  dum_string = XML_Get_Attribute_Value(xnd,"mcmc.sample.every");
  if(dum_string != NULL) mixt_tree->io->mcmc->sample_interval = (int)String_To_Dbl(dum_string);

  dum_string = XML_Get_Attribute_Value(xnd,"mcmc.print.every");
  if(dum_string != NULL) mixt_tree->io->mcmc->print_every = (int)String_To_Dbl(dum_string);

  dum_string = XML_Get_Attribute_Value(xnd,"mcmc.burnin");
  if(dum_string != NULL) mixt_tree->io->mcmc->chain_len_burnin = (int)String_To_Dbl(dum_string);


  dum_string = XML_Get_Attribute_Value(xnd,"ignore.sequences");
  if(!dum_string) dum_string = XML_Get_Attribute_Value(xnd,"ignore.seq");
  if(!dum_string) dum_string = XML_Get_Attribute_Value(xnd,"ignore.data");

  if(dum_string != NULL)
    {
      int select = XML_Validate_Attr_Int(dum_string,6,
                                         "true","yes","y",
                                         "false","no","n");
      if(select < 3) mixt_tree->eval_alnL = NO;
      else mixt_tree->eval_alnL = YES;
    }
  
  dum_string = XML_Get_Attribute_Value(xnd,"mutmap");
  if(dum_string != NULL)
    {
      int select = XML_Validate_Attr_Int(dum_string,6,
                                         "true","yes","y",
                                         "false","no","n");
      if(select < 3) mixt_tree->io->mutmap = YES;
      else mixt_tree->io->mutmap = NO;
    }

    
  // Looking for XML node with rate-across-lineage info
  xnd = XML_Search_Node_Name("lineagerates",YES,xroot);
  
  if(xnd == NULL)
    {
      PhyML_Fprintf(stdout,"\n. The model of rate variation across lineages is not specified.");
      PhyML_Fprintf(stdout,"\n. Using the geometric Brownian model (see Guindon, 2012, Syst. Biol.).\n");
      mixt_tree->rates->model      = GUINDON;
      mixt_tree->mod->gamma_mgf_bl = YES;
      strcpy(mixt_tree->rates->model_name,"geometric Brownian"); 
    }
  else
    {
      char *model_name;
      model_name = XML_Get_Attribute_Value(xnd,"model");

      if(model_name == NULL)
        {
          PhyML_Fprintf(stderr,"\n. Please specify a model of rate variation across lineages,");
          PhyML_Fprintf(stderr,"\n. e.g., <lineagerates model=\"geometricbrownian\"/>.");
          PhyML_Fprintf(stderr,"\n. See the manual for more options.");
          assert(FALSE);
        }
      else
        {
          if(!strcmp(model_name,"geometricbrownian"))
            {
              mixt_tree->rates->model      = GUINDON;
              mixt_tree->mod->gamma_mgf_bl = YES;
              strcpy(mixt_tree->rates->model_name,"geometric Brownian"); 
            }
          else if(!strcmp(model_name,"geometric"))
            {
              mixt_tree->rates->model      = GUINDON;
              mixt_tree->mod->gamma_mgf_bl = YES;
              strcpy(mixt_tree->rates->model_name,"geometric Brownian"); 
            }
          else if(!strcmp(model_name,"brownian"))
            {
              mixt_tree->rates->model      = GUINDON;
              mixt_tree->mod->gamma_mgf_bl = YES;
              strcpy(mixt_tree->rates->model_name,"geometric Brownian"); 
            }
          else if(!strcmp(model_name,"geo"))
            {
              mixt_tree->rates->model      = GUINDON;
              mixt_tree->mod->gamma_mgf_bl = YES;
              strcpy(mixt_tree->rates->model_name,"geometric Brownian"); 
            }
          else if(!strcmp(model_name,"lognormal"))
            {
              mixt_tree->rates->model      = LOGNORMAL;
              mixt_tree->mod->gamma_mgf_bl = NO;
              strcpy(mixt_tree->rates->model_name,"lognormal (uncorrelated)"); 
            }
          else if(!strcmp(model_name,"normal"))
            {
              mixt_tree->rates->model      = LOGNORMAL;
              mixt_tree->mod->gamma_mgf_bl = NO;
              strcpy(mixt_tree->rates->model_name,"lognormal (uncorrelated)"); 
            }
          else if(!strcmp(model_name,"strictclock"))
            {
              mixt_tree->rates->model      = STRICTCLOCK;
              mixt_tree->mod->gamma_mgf_bl = NO;
              strcpy(mixt_tree->rates->model_name,"strict clock"); 
            }
          else if(!strcmp(model_name,"clock"))
            {
              mixt_tree->rates->model      = STRICTCLOCK;
              mixt_tree->mod->gamma_mgf_bl = NO;
              strcpy(mixt_tree->rates->model_name,"strict clock"); 
            }
          else
            {
              assert(FALSE);
            }
        }      
    }

  // Spatial model
  xnd = XML_Search_Node_Name("spatialmodel",YES,xroot);

  if(xnd != NULL)
    {
      char *modname;
      modname = XML_Get_Attribute_Value(xnd,"name");
      XML_Set_Attribute_Value(xnd,"name",modname);

      if(modname != NULL)
        {
          if(!strcmp(xnd->attr->value,"slfv")) mixt_tree->mmod->id = SLFV_GAUSSIAN;
          else if(!strcmp(xnd->attr->value,"rw")) mixt_tree->mmod->id = RW;
          else if(!strcmp(xnd->attr->value,"rrw")) mixt_tree->mmod->id = RRW;
          else assert(FALSE);
        }
    }
  

  // Looking for XML node with rate-across-lineage info
  xnd = XML_Search_Node_Name("clockrate",YES,xroot);
  
  if(xnd != NULL)
    {
      char *clock_r;
      clock_r = XML_Get_Attribute_Value(xnd,"value");
      if(clock_r == NULL) clock_r = XML_Get_Attribute_Value(xnd,"clock.val");
      if(clock_r == NULL) clock_r = XML_Get_Attribute_Value(xnd,"val");
      if(clock_r != NULL)
        {
          mixt_tree->rates->clock_r = String_To_Dbl(clock_r);
        }
      
      char *opt_clock;
      opt_clock = XML_Get_Attribute_Value(xnd,"optimise.clock");
      if(opt_clock == NULL) opt_clock = XML_Get_Attribute_Value(xnd,"optimize.clock");
      if(opt_clock == NULL) opt_clock = XML_Get_Attribute_Value(xnd,"optimize.rate");
      if(opt_clock == NULL) opt_clock = XML_Get_Attribute_Value(xnd,"opt.clock");
      
      if(opt_clock != NULL)
        {
          int select = XML_Validate_Attr_Int(opt_clock,6,
                                             "true","yes","y",
                                             "false","no","n");
          if(select < 3)  mixt_tree->mod->s_opt->opt_clock_r = YES;
          else mixt_tree->mod->s_opt->opt_clock_r = NO;
        }
    }

  
  // Looking for calibration info
  xnd = XML_Search_Node_Name("calibration",YES,xroot);

  if(xnd == NULL)
    {
      PhyML_Fprintf(stderr,"\n. No calibration information seems to be provided.");
      PhyML_Fprintf(stderr,"\n. Please amend your XML file. \n");
      assert(FALSE);
    }
  else
    {
      assert(xnd->child);
      if(XML_Search_Node_Name("upper",NO,xnd->child) == NULL && XML_Search_Node_Name("lower",NO,xnd->child) == NULL)
        {
          PhyML_Fprintf(stderr,"\n. There is no calibration information provided. \n");
          PhyML_Fprintf(stderr,"\n. Please check your data. \n");
          assert(FALSE);
        }
    }

  

  // Looking for coordinate file
  xnd = XML_Search_Node_Name("coordinates",YES,xroot);

  if(xnd == NULL)
    {
      PhyML_Fprintf(stderr,"\n. No spatial information (i.e., coordinates) seems to be provided.");
      PhyML_Fprintf(stderr,"\n. Please amend your XML file. \n");
      Exit("\n");
    }
  else
    {
      char *coord_file;
      coord_file = XML_Get_Attribute_Value(xnd,"file.name");
      
      strcpy(mixt_tree->io->in_coord_file,coord_file);
      mixt_tree->io->fp_in_coord = Openfile(mixt_tree->io->in_coord_file,READ);
    }

  seed = (mixt_tree->io->r_seed < 0)?(time(NULL)):(mixt_tree->io->r_seed);
  srand(seed);
  mixt_tree->io->r_seed = seed;

  
  MIXT_Check_Model_Validity(mixt_tree);
  MIXT_Init_Model(mixt_tree);
  Print_Data_Structure(NO,stdout,mixt_tree);
  tree = MIXT_Starting_Tree(mixt_tree);
  if(mixt_tree->io->in_tree < 2) Add_Root(tree->a_edges[0],tree);
  Copy_Tree(tree,mixt_tree);
  Free_Tree(tree);
  MIXT_Connect_Cseqs_To_Nodes(mixt_tree);
  MIXT_Init_T_Beg(mixt_tree);  
  MIXT_Make_Tree_For_Lk(mixt_tree);
  MIXT_Make_Tree_For_Pars(mixt_tree);
  MIXT_Make_Spr(mixt_tree);  
  MIXT_Chain_All(mixt_tree);
  MIXT_Check_Edge_Lens_In_All_Elem(mixt_tree);
  MIXT_Turn_Branches_OnOff_In_All_Elem(ON,mixt_tree);
  MIXT_Check_Invar_Struct_In_Each_Partition_Elem(mixt_tree);
  MIXT_Check_RAS_Struct_In_Each_Partition_Elem(mixt_tree);

  
  XML_Read_Calibration(xroot,mixt_tree);
  MIXT_Chain_Cal(mixt_tree);


  if(TIMES_Calibrations_Apply_To_Tips_Only(mixt_tree) == YES &&
     mixt_tree->mod->s_opt->opt_topo == NO)
    {
      TIMES_Randomize_Tip_Times_Given_Calibrations(mixt_tree); // Topology is unchanged
      TIMES_Bl_To_Times(mixt_tree);
      Update_Ancestors(mixt_tree->n_root,mixt_tree->n_root->v[2],mixt_tree);
      Update_Ancestors(mixt_tree->n_root,mixt_tree->n_root->v[1],mixt_tree);
    }
  else
    {
      TIMES_Randomize_Tree_With_Time_Constraints(mixt_tree->times->a_cal[0],mixt_tree);
    }
  
  MIXT_Propagate_Tree_Update(mixt_tree);

  /* Create ldsks and connect tree tips to them */
  /* once tip dates have been set properly (in */
  /* TIMES_Randomize_Tree_With_Time_Constraints) */
  PHYREX_Make_And_Connect_Tip_Disks(mixt_tree);

  
  /* Read spatial coordinates */
  PHYREX_Read_Tip_Coordinates(mixt_tree);
  PHYREX_Init_Migrep_Mod(mixt_tree->mmod,2,
                         mixt_tree->mmod->lim_do->lonlat[0],
                         mixt_tree->mmod->lim_do->lonlat[1],
                         mixt_tree->mmod->lim_up->lonlat[0],
                         mixt_tree->mmod->lim_up->lonlat[1]);

  

  /* Initialize parameters of migrep model */
  mixt_tree->mmod->lbda = 1.;
  mixt_tree->mmod->mu   = 0.8;

  if(mixt_tree->mmod->id == SLFV_GAUSSIAN || mixt_tree->mmod->id == SLFV_UNIFORM)
    mixt_tree->mmod->rad  = 0.1*((mixt_tree->mmod->lim_up->lonlat[0]-mixt_tree->mmod->lim_do->lonlat[0])+
                                 (mixt_tree->mmod->lim_up->lonlat[1]-mixt_tree->mmod->lim_do->lonlat[1]));
  else if(mixt_tree->mmod->id == RW || mixt_tree->mmod->id == RRW)
    mixt_tree->mmod->rad  = 2.0*((mixt_tree->mmod->lim_up->lonlat[0]-mixt_tree->mmod->lim_do->lonlat[0])+
                                 (mixt_tree->mmod->lim_up->lonlat[1]-mixt_tree->mmod->lim_do->lonlat[1]));
  else assert(FALSE);

  mixt_tree->mmod->sigsq = PHYREX_Update_Sigsq(mixt_tree);

  
  /* Random genealogy or user-defined tree */
  switch(mixt_tree->io->in_tree)
    {
    case 0 : case 1 :
      {
        PHYREX_Simulate_Backward_Core(mixt_tree->young_disk,YES,mixt_tree);
        PHYREX_Ldsk_To_Tree(mixt_tree);
        break;
      }
    case 2:
      {
        PHYREX_Tree_To_Ldsk(mixt_tree);
        break;
      }
    }

  Update_Ancestors(mixt_tree->n_root,mixt_tree->n_root->v[2],mixt_tree);
  Update_Ancestors(mixt_tree->n_root,mixt_tree->n_root->v[1],mixt_tree);  
  
  MIXT_Set_Ignore_Root(YES,mixt_tree);
  MIXT_Set_Bl_From_Rt(YES,mixt_tree);

  PHYREX_Oldest_Sampled_Disk(mixt_tree);

  if(mixt_tree->mmod->id == RW || mixt_tree->mmod->id == RRW)
    {
      Make_Contrasts(mixt_tree);
      mixt_tree->mmod->max_num_of_intervals = 3*tree->n_otu-2;
  }

  assert(PHYREX_Check_Struct(mixt_tree));
  PHYREX_Lk(mixt_tree);
  Set_Update_Eigen(YES,mixt_tree->mod);
  Lk(NULL,mixt_tree);
  Set_Update_Eigen(NO,mixt_tree->mod);
  PhyML_Printf("\n. Init lnPr(seq|phylo): %f lnPr(coor|phylo): %f",mixt_tree->c_lnL,mixt_tree->mmod->c_lnL);
  PhyML_Printf("\n. Random seed: %d",mixt_tree->io->r_seed);
    
  res = PHYREX_MCMC(mixt_tree);
  Free(res);  
  
  // Cleaning up...
  RATES_Free_Rates(mixt_tree->rates);
  RATES_Free_Rates(mixt_tree->extra_tree->rates);
  TIMES_Free_Times(mixt_tree->times);
  TIMES_Free_Times(mixt_tree->extra_tree->times);
  MCMC_Free_MCMC(mixt_tree->mcmc);
  MCMC_Free_MCMC(mixt_tree->extra_tree->mcmc);
  Free_Mmod(mixt_tree->mmod);
  Free_Spr_List_One_Edge(mixt_tree);
  Free_Tree_Pars(mixt_tree);
  Free_Tree_Lk(mixt_tree);

  if(mixt_tree->io->fp_out_trees)      fclose(mixt_tree->io->fp_out_trees);
  if(mixt_tree->io->fp_out_tree)       fclose(mixt_tree->io->fp_out_tree);
  if(mixt_tree->io->fp_out_stats)      fclose(mixt_tree->io->fp_out_stats);
  if(mixt_tree->io->fp_out_json_trace) fclose(mixt_tree->io->fp_out_json_trace);
  Free_Input(mixt_tree->io);


  tree = mixt_tree;
  do
    {
      Free_Calign(tree->data);
      tree = tree->next_mixt;
    }
  while(tree);

  tree = mixt_tree;
  do
    {
      Free_Optimiz(tree->mod->s_opt);
      tree = tree->next;
    }
  while(tree);

  
  Free_Model_Complete(mixt_tree->mod);
  Free_Model_Basic(mixt_tree->mod);
  Free_Tree(mixt_tree->extra_tree);  
  Free_Tree(mixt_tree);  
  Free(res);
  XML_Free_XML_Tree(xroot);
  fclose(fp_xml_in);
}
#endif

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int PHYREX_Main_Simulate(int argc, char *argv[])
{
  t_tree *tree;
  int seed,i,modid;
  char *s;
  t_dsk *disk;
  int n_sites,n_otus;
  phydbl width,height;
  phydbl lbda, rad,mu;
  
  s = (char *)mCalloc(T_MAX_NAME,sizeof(char));

  n_otus  = (int)atoi(argv[1]);
  n_sites = 1;
  width   = (phydbl)atof(argv[2]); 
  height  = (phydbl)atof(argv[3]); 
  lbda    = (phydbl)atof(argv[4]);
  rad     = (phydbl)atof(argv[5]);
  mu      = (phydbl)atof(argv[6]);
  seed    = (int)atoi(argv[7]);
  modid   = (int)atoi(argv[8]);
  
  printf("\n. seed: %d",seed);
  srand(seed);
  
  tree = PHYREX_Simulate(n_otus,n_sites,width,height,lbda,rad,mu,seed,modid);
  /* tree = PHYREX_Simulate_Independent_Loci(n_otus,500,20.,20.,seed); */

  disk = tree->young_disk;
  for(i=0;i<disk->n_ldsk_a;i++) Free_Ldisk(disk->ldsk_a[i]);
  while(disk->prev)
    {
      disk = disk->prev;
      if(disk->next->ldsk != NULL) Free_Ldisk(disk->next->ldsk);
      Free_Disk(disk->next);
    }
  
  /* Root */
  Free_Ldisk(disk->ldsk);
  Free_Disk(disk);

  RATES_Free_Rates(tree->rates);
  /* MCMC_Free_MCMC(tree->mcmc); */
  Free_Mmod(tree->mmod);
  Free_Spr_List_One_Edge(tree);
  Free_Spr_List_All_Edge(tree);
  Free_Tree_Pars(tree);
  Free_Tree_Lk(tree);
  Free_Input(tree->io);
  Free_Optimiz(tree->mod->s_opt);
  Free_Model_Complete(tree->mod);
  Free_Model_Basic(tree->mod);
  Free_Calign(tree->data);
  Free_Tree(tree);
  Free(s);

  return 0;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
/* Test whether coord is in disk. Will actually only works in disk */
/* is a rectangle... */

int PHYREX_Is_In_Disk(t_geo_coord *coord, t_dsk *disk, t_phyrex_mod *mmod)
{
  int i;

  assert(disk->centr->dim);

  if(mmod->id == SLFV_UNIFORM)
    {
      for(i=0;i<disk->centr->dim;i++)
        {
          if(fabs(coord->lonlat[i] - disk->centr->lonlat[i]) > disk->mmod->rad + 1.E-20)
            {
              return(NO);
            }
        }
      return(YES);
    }
  else if(mmod->id == SLFV_GAUSSIAN)
    {
      return(YES);
    }

  return(-1);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl PHYREX_Lk(t_tree *tree)
{
  switch(tree->mmod->id)
    {
    case SLFV_GAUSSIAN : case SLFV_UNIFORM :
      {
        return(SLFV_Lk_Gaussian(tree));
        break;
      }
    case RW :
      {
        return(RW_Lk(tree));
        break;
      }
    case RRW :
      {
        return(RRW_Lk(tree));
        break;
      }
    default : assert(FALSE);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl PHYREX_Lk_Core(t_dsk *disk, t_tree *tree)
{
  switch(tree->mmod->id)
    {
    case SLFV_GAUSSIAN : case SLFV_UNIFORM :
      {
        return(SLFV_Lk_Gaussian_Core(disk,tree));
        break;
      }
    case RW :
        {
          return(RW_Lk(tree));
          break;
        }
    case RRW :
        {
          return(RRW_Lk(tree));
          break;
        }
    default : assert(FALSE);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// Warning: the calculation below does not incoporate time information
// since things get messy when considering serial samples
phydbl PHYREX_Lk_Range(t_dsk *young, t_dsk *old, t_tree *tree)
{
  switch(tree->mmod->id)
    {
    case SLFV_GAUSSIAN : case SLFV_UNIFORM :
      {
        return(SLFV_Lk_Gaussian_Range(young,old,tree));
        break;
      }
    case RW :
        {
          return(RW_Lk(tree));
          break;
        }
    case RRW :
        {
          return(RRW_Lk(tree));
          break;
        }
    default : assert(FALSE);
    }
  
  return(-1.);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

#if (defined PHYREX)
phydbl *PHYREX_MCMC(t_tree *tree)
{
  t_mcmc *mcmc;
  int move,i,n_vars;
  phydbl u;
  FILE *fp_tree,*fp_stats;
  phydbl *res;
  t_dsk *disk;

  fp_tree    = tree->io->fp_out_tree;
  fp_stats   = tree->io->fp_out_stats;

  if(tree->io->mcmc == NULL)
    {
      mcmc = MCMC_Make_MCMC_Struct();
      tree->mcmc = mcmc;
    }
  else
    {
      tree->mcmc = tree->io->mcmc;
      mcmc = tree->mcmc;
    }

  
  mcmc->io               = NULL;
  mcmc->is               = NO;
  mcmc->run              = 0;
  mcmc->randomize        = YES;
  mcmc->norm_freq        = 1E+3;
  mcmc->max_tune         = 1.E+20;
  mcmc->is_burnin        = YES;
  mcmc->nd_t_digits      = 1;
  mcmc->max_lag          = 1000;
  mcmc->sample_size      = mcmc->chain_len/mcmc->sample_interval;
  mcmc->sample_num       = 0;
  disk                   = NULL;
    
  PhyML_Fprintf(fp_tree,"#NEXUS");
  PhyML_Fprintf(fp_tree,"\n\nBegin taxa;");
  PhyML_Fprintf(fp_tree,"\n\tDimensions ntax=%d;",tree->n_otu);
  PhyML_Fprintf(fp_tree,"\n\tTaxlabels");
  for(int i=0;i<tree->n_otu;++i)
    {
      PhyML_Fprintf(fp_tree,"\n\t\t'%s'",tree->a_nodes[i]->name);
    }
  PhyML_Fprintf(fp_tree,"\n\t;");
  PhyML_Fprintf(fp_tree,"\n\tend;");
  PhyML_Fprintf(fp_tree,"\n\n");
  PhyML_Fprintf(fp_tree,"\nBegin trees;");
  PhyML_Fprintf(fp_tree,"\n\tTranslate");
  for(int i=0;i<tree->n_otu;++i)
    {
      PhyML_Fprintf(fp_tree,"\n\t%d '%s'",i+1,tree->a_nodes[i]->name);
      if(i<tree->n_otu-1) PhyML_Fprintf(fp_tree,",");
    }
  PhyML_Fprintf(fp_tree,"\n;");
  
  
  MCMC_Complete_MCMC(mcmc,tree);

  n_vars = 12;

  res = (phydbl *)mCalloc(tree->mcmc->chain_len / tree->mcmc->sample_interval * n_vars,sizeof(phydbl));

  
  PHYREX_Update_Sigsq(tree);
  
  MIXT_Set_Bl_From_Rt(YES,tree);

  if(XML_Search_Node_Attribute_Value("add","true",NO,tree->xml_root) == NULL)
    {
      PhyML_Fprintf(fp_stats,"\n# before rand glnL: %f alnL: %f",tree->mmod->c_lnL,tree->c_lnL);
      PhyML_Fprintf(fp_stats,"\n# ninter: %d",PHYREX_Total_Number_Of_Intervals(tree));
      PhyML_Fprintf(fp_stats,"\n# ncoal: %d",PHYREX_Total_Number_Of_Coal_Disks(tree));
      PhyML_Fprintf(fp_stats,"\n# nhits: %d",PHYREX_Total_Number_Of_Hit_Disks(tree));
      PhyML_Fprintf(fp_stats,"\n# true lbda: %f",tree->mmod->lbda);
      PhyML_Fprintf(fp_stats,"\n# true mu: %f",tree->mmod->mu);
      PhyML_Fprintf(fp_stats,"\n# true rad: %f",SLFV_Update_Radius(tree));
      PhyML_Fprintf(fp_stats,"\n# true sigsq: %f",tree->mmod->sigsq);
      PhyML_Fprintf(fp_stats,"\n# true neigh. size: %f",SLFV_Neighborhood_Size(tree));
      PhyML_Fprintf(fp_stats,"\n# fst-based estimate of neighborhood size: %f",SLFV_Neighborhood_Size_Regression(tree));
      PhyML_Fprintf(fp_stats,"\n# nucleotide diversity: %f",Nucleotide_Diversity(tree->data));
      PhyML_Fprintf(fp_stats,"\n# length of a generation: %G time units",SLFV_Generation_Length(tree));
      PhyML_Fprintf(fp_stats,"\n# clock rate: %G subst. per time unit",tree->rates->clock_r);      
      PhyML_Fprintf(fp_stats,"\n# after rand glnL: %f alnL: %f",tree->mmod->c_lnL,tree->c_lnL);
      PhyML_Fprintf(fp_stats,"\n# ninter: %d",PHYREX_Total_Number_Of_Intervals(tree));
      PhyML_Fprintf(fp_stats,"\n# ncoal: %d",PHYREX_Total_Number_Of_Coal_Disks(tree));
      PhyML_Fprintf(fp_stats,"\n# nhits: %d",PHYREX_Total_Number_Of_Hit_Disks(tree));
      PhyML_Fprintf(fp_stats,"\n# start lbda: %f",tree->mmod->lbda);
      PhyML_Fprintf(fp_stats,"\n# start mu: %f",tree->mmod->mu);
      PhyML_Fprintf(fp_stats,"\n# start rad: %f",tree->mmod->rad);
      PhyML_Fprintf(fp_stats,"\n# dist. in tree: ");
      for(int i=0;i<tree->n_otu-1;++i) for(int j=i+1;j<tree->n_otu;++j) PhyML_Fprintf(fp_stats,"%G ",PHYREX_Dist_Between_Two_Ldsk(tree->a_nodes[i]->ldsk,tree->a_nodes[j]->ldsk,tree));
      PhyML_Fprintf(fp_stats,"\n# dist. in space: ");
      for(int i=0;i<tree->n_otu-1;++i) for(int j=i+1;j<tree->n_otu;++j) PhyML_Fprintf(fp_stats,"%G ",Euclidean_Dist(tree->a_nodes[i]->ldsk->coord,tree->a_nodes[j]->ldsk->coord));
      PhyML_Printf("\n");
      
      fflush(NULL);      
      
      PhyML_Fprintf(fp_stats,"\n");
      PhyML_Fprintf(fp_stats,"%s\t","sample");
      PhyML_Fprintf(fp_stats,"%s\t","lnP");
      PhyML_Fprintf(fp_stats,"%s\t","alnL");
      PhyML_Fprintf(fp_stats,"%s\t","glnL");
      PhyML_Fprintf(fp_stats,"%s\t","clock");
      PhyML_Fprintf(fp_stats,"%s\t","evolrate");
      PhyML_Fprintf(fp_stats,"%s\t","lbda");
      PhyML_Fprintf(fp_stats,"%s\t","mu");
      PhyML_Fprintf(fp_stats,"%s\t","rad");
      PhyML_Fprintf(fp_stats,"%s\t","sigsq");
      PhyML_Fprintf(fp_stats,"%s\t","neff");
      PhyML_Fprintf(fp_stats,"%s\t","neigh");
      PhyML_Fprintf(fp_stats,"%s\t","rhoe");
      PhyML_Fprintf(fp_stats,"%s\t","realsigsqroot");
      PhyML_Fprintf(fp_stats,"%s\t","realsigsqtips");
      PhyML_Fprintf(fp_stats,"%s\t","realsigsqtipsbis");
      PhyML_Fprintf(fp_stats,"%s\t","realsigsqtipster");
      PhyML_Fprintf(fp_stats,"%s\t","dispdist");
      PhyML_Fprintf(fp_stats,"%s\t","nInt");
      PhyML_Fprintf(fp_stats,"%s\t","nCoal");
      PhyML_Fprintf(fp_stats,"%s\t","nHit");
      PhyML_Fprintf(fp_stats,"%s\t","rootTime");
      PhyML_Fprintf(fp_stats,"%s\t","rootLon");
      PhyML_Fprintf(fp_stats,"%s\t","rootLat");
      for(int i=0;i<Scalar_Len(tree->mod->kappa);++i) PhyML_Fprintf(fp_stats,"tstv%d\t",i);
      if(tree->mod->ras->free_mixt_rates == NO) PhyML_Fprintf(fp_stats,"alpha\t");
      else
        {
          for(int i=0;i<tree->mod->ras->n_catg;++i) PhyML_Fprintf(fp_stats,"p(%d)\t",i+1);
          for(int i=0;i<tree->mod->ras->n_catg;++i) PhyML_Fprintf(fp_stats,"rr(%d)\t",i+1);
        }

      PhyML_Fprintf(fp_stats,"%s\t","MeanBr");
      PhyML_Fprintf(fp_stats,"%s\t","TreeLen");
      PhyML_Fprintf(fp_stats,"%s\t","accLbda");
      PhyML_Fprintf(fp_stats,"%s\t","accMu");
      PhyML_Fprintf(fp_stats,"%s\t","accRad");
      PhyML_Fprintf(fp_stats,"%s\t","accInDelDisk");
      PhyML_Fprintf(fp_stats,"%s\t","accInDelHit");
      PhyML_Fprintf(fp_stats,"%s\t","accScaleTime");
      PhyML_Fprintf(fp_stats,"%s\t","accSPR");
      PhyML_Fprintf(fp_stats,"%s\t","accSPRlocal");
      PhyML_Fprintf(fp_stats,"%s\t","accPath");
      PhyML_Fprintf(fp_stats,"%s\t","accIndelDiskSerial");
      PhyML_Fprintf(fp_stats,"%s\t","accIndelHitSerial");
      PhyML_Fprintf(fp_stats,"%s\t","accLdskGivenDisk");
      PhyML_Fprintf(fp_stats,"%s\t","accDiskGivenLdsk");
      PhyML_Fprintf(fp_stats,"%s\t","accDiskAndLdsk");
      PhyML_Fprintf(fp_stats,"%s\t","accLdskMulti");
      PhyML_Fprintf(fp_stats,"%s\t","accDiskMulti");
      PhyML_Fprintf(fp_stats,"%s\t","accMoveDiskUD");
      PhyML_Fprintf(fp_stats,"%s\t","accAddRemoveJump");
      PhyML_Fprintf(fp_stats,"%s\t","accLdskTipToRoot");
      PhyML_Fprintf(fp_stats,"%s\t","accSigsqScale");
      PhyML_Fprintf(fp_stats,"%s\t","tuneLbda");
      PhyML_Fprintf(fp_stats,"%s\t","tuneRad");
      PhyML_Fprintf(fp_stats,"%s\t","tuneMu");
      PhyML_Fprintf(fp_stats,"%s\t","tuneIndelDisk");
      PhyML_Fprintf(fp_stats,"%s\t","tuneIndelHit");
      PhyML_Fprintf(fp_stats,"%s\t","tuneLdskGivenDisk");
      PhyML_Fprintf(fp_stats,"%s\t","tuneIndelDiskSerial");
      PhyML_Fprintf(fp_stats,"%s\t","tuneIndelHitSerial");
      for(int i=0;i<2*tree->n_otu-2;++i) PhyML_Fprintf(fp_stats,"s%d\t",i);
    }
  

  for(i=0;i<mcmc->n_moves;i++) tree->mcmc->start_ess[i] = YES;

  PHYREX_Lk(tree);        
  Set_Update_Eigen(YES,tree->mod);
  Lk(NULL,tree);
  Set_Update_Eigen(NO,tree->mod);
  RATES_Lk_Rates(tree);
  
  if(isnan(tree->c_lnL) || isinf(tree->c_lnL))
    {
      PhyML_Fprintf(stderr,"\n. Cannot compute sequence log-likelihood. Aborting.");
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);            
    }

  if(isnan(tree->mmod->c_lnL) || isinf(tree->mmod->c_lnL))
    {
      PhyML_Fprintf(stderr,"\n. Cannot compute location log-likelihood. Aborting.");
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);            
    }
  

  Set_Both_Sides(NO,tree);
  mcmc->always_yes = NO;
  move             = -1;
  do
    {
      MIXT_Propagate_Tree_Update(tree);
      assert(PHYREX_Check_Struct(tree));

      if(mcmc->run > mcmc->chain_len_burnin)
        for(i=0;i<mcmc->n_moves;i++) tree->mcmc->adjust_tuning[i] = NO;
      else
        {
          for(i=0;i<mcmc->n_moves;i++) tree->mcmc->adjust_tuning[i] = YES;
          tree->mcmc->is_burnin = NO;
        }
            
      u = Uni();

      for(move=0;move<tree->mcmc->n_moves;move++) if(tree->mcmc->move_weight[move] > u-1.E-10) break;

      tree->mcmc->move_idx = move;
      
      assert(!(move == tree->mcmc->n_moves));
      
      if(!(tree->mcmc->run%tree->mcmc->print_every))
        {
          PhyML_Fprintf(stdout,"\n. %10d %30s %20f %20f %20f",
                        tree->mcmc->run,
                        tree->mcmc->move_name[move],
                        tree->mmod->c_lnL,
                        tree->c_lnL,
                        tree->mmod->c_lnL+tree->c_lnL+tree->rates->c_lnL_rates);
          if(tree->numerical_warning == YES) PhyML_Fprintf(stdout," -- WARNING: numerical precision issue detected...");
        }

            
      if(!(tree->mcmc->run%tree->mcmc->sample_interval))
        {
          disk = tree->young_disk;
          while(disk->prev) disk = disk->prev;

          PhyML_Fprintf(fp_stats,"\n");
          PhyML_Fprintf(fp_stats,"%6d\t",tree->mcmc->run);
          PhyML_Fprintf(fp_stats,"%g\t",tree->c_lnL+tree->mmod->c_lnL+tree->rates->c_lnL_rates);
          PhyML_Fprintf(fp_stats,"%g\t",tree->c_lnL);
          PhyML_Fprintf(fp_stats,"%g\t",tree->mmod->c_lnL);
          PhyML_Fprintf(fp_stats,"%g\t",tree->rates->clock_r);
          PhyML_Fprintf(fp_stats,"%g\t",RATES_Realized_Substitution_Rate(tree));
          PhyML_Fprintf(fp_stats,"%g\t",tree->mmod->lbda);
          PhyML_Fprintf(fp_stats,"%g\t",tree->mmod->mu);
          PhyML_Fprintf(fp_stats,"%g\t",tree->mmod->rad);
          PhyML_Fprintf(fp_stats,"%g\t",tree->mmod->sigsq);
          PhyML_Fprintf(fp_stats,"%g\t",tree->times->scaled_pop_size);
          PhyML_Fprintf(fp_stats,"%g\t",SLFV_Neighborhood_Size(tree));
          PhyML_Fprintf(fp_stats,"%g\t",SLFV_Effective_Density(tree));
          PhyML_Fprintf(fp_stats,"%g\t",PHYREX_Root_To_Tip_Realized_Sigsq(tree));
          PhyML_Fprintf(fp_stats,"%g\t",PHYREX_Tip_To_Root_Realized_Sigsq(tree));
          PhyML_Fprintf(fp_stats,"%g\t",PHYREX_Tip_To_Root_Realized_Bis_Sigsq(tree));
          PhyML_Fprintf(fp_stats,"%g\t",PHYREX_Tip_To_Root_Realized_Ter_Sigsq(tree));
          PhyML_Fprintf(fp_stats,"%g\t",PHYREX_Realized_Dispersal_Dist(tree));
          PhyML_Fprintf(fp_stats,"%d\t",PHYREX_Total_Number_Of_Intervals(tree));
          PhyML_Fprintf(fp_stats,"%d\t",PHYREX_Total_Number_Of_Coal_Disks(tree));
          PhyML_Fprintf(fp_stats,"%d\t",PHYREX_Total_Number_Of_Hit_Disks(tree));
          PhyML_Fprintf(fp_stats,"%g\t",disk->time);
          PhyML_Fprintf(fp_stats,"%g\t",disk->ldsk->coord->lonlat[0]);
          PhyML_Fprintf(fp_stats,"%g\t",disk->ldsk->coord->lonlat[1]);
          Output_Scalar_Dbl(tree->mod->kappa,"\t",fp_stats);
          if(tree->mod->ras->free_mixt_rates == NO) PhyML_Fprintf(fp_stats,"%g\t",tree->mod->ras->alpha->v);
          else
            {
              for(int i=0;i<tree->mod->ras->n_catg;++i) PhyML_Fprintf(fp_stats,"%g\t",tree->mod->ras->gamma_r_proba->v[i]);
              for(int i=0;i<tree->mod->ras->n_catg;++i) PhyML_Fprintf(fp_stats,"%g\t",tree->mod->ras->gamma_rr->v[i]);
            }
          PhyML_Fprintf(fp_stats,"%g\t",RATES_Get_Mean_Rate_In_Subtree(tree->n_root,tree));
          PhyML_Fprintf(fp_stats,"%g\t",Tree_Length(tree));
          PhyML_Fprintf(fp_stats,"%g\t",tree->mcmc->acc_rate[tree->mcmc->num_move_phyrex_lbda]);
          PhyML_Fprintf(fp_stats,"%g\t",tree->mcmc->acc_rate[tree->mcmc->num_move_phyrex_mu]);
          PhyML_Fprintf(fp_stats,"%g\t",tree->mcmc->acc_rate[tree->mcmc->num_move_phyrex_rad]);
          PhyML_Fprintf(fp_stats,"%g\t",tree->mcmc->acc_rate[tree->mcmc->num_move_phyrex_indel_disk]);
          PhyML_Fprintf(fp_stats,"%g\t",tree->mcmc->acc_rate[tree->mcmc->num_move_phyrex_indel_hit]);
          PhyML_Fprintf(fp_stats,"%g\t",tree->mcmc->acc_rate[tree->mcmc->num_move_phyrex_scale_times]);
          PhyML_Fprintf(fp_stats,"%g\t",tree->mcmc->acc_rate[tree->mcmc->num_move_phyrex_spr]);
          PhyML_Fprintf(fp_stats,"%g\t",tree->mcmc->acc_rate[tree->mcmc->num_move_phyrex_spr_local]);
          PhyML_Fprintf(fp_stats,"%g\t",tree->mcmc->acc_rate[tree->mcmc->num_move_phyrex_traj]);
          PhyML_Fprintf(fp_stats,"%g\t",tree->mcmc->acc_rate[tree->mcmc->num_move_phyrex_indel_disk_serial]);
          PhyML_Fprintf(fp_stats,"%g\t",tree->mcmc->acc_rate[tree->mcmc->num_move_phyrex_indel_hit_serial]);
          PhyML_Fprintf(fp_stats,"%g\t",tree->mcmc->acc_rate[tree->mcmc->num_move_phyrex_ldsk_given_disk]);
          PhyML_Fprintf(fp_stats,"%g\t",tree->mcmc->acc_rate[tree->mcmc->num_move_phyrex_disk_given_ldsk]);
          PhyML_Fprintf(fp_stats,"%g\t",tree->mcmc->acc_rate[tree->mcmc->num_move_phyrex_ldsk_and_disk]);
          PhyML_Fprintf(fp_stats,"%g\t",tree->mcmc->acc_rate[tree->mcmc->num_move_phyrex_ldsk_multi]);
          PhyML_Fprintf(fp_stats,"%g\t",tree->mcmc->acc_rate[tree->mcmc->num_move_phyrex_disk_multi]);
          PhyML_Fprintf(fp_stats,"%g\t",tree->mcmc->acc_rate[tree->mcmc->num_move_phyrex_move_disk_ud]);
          PhyML_Fprintf(fp_stats,"%g\t",tree->mcmc->acc_rate[tree->mcmc->num_move_phyrex_add_remove_jump]);
          PhyML_Fprintf(fp_stats,"%g\t",tree->mcmc->acc_rate[tree->mcmc->num_move_phyrex_ldsk_tip_to_root]);
          PhyML_Fprintf(fp_stats,"%g\t",tree->mcmc->acc_rate[tree->mcmc->num_move_phyrex_sigsq_scale]);
          PhyML_Fprintf(fp_stats,"%g\t",tree->mcmc->tune_move[tree->mcmc->num_move_phyrex_lbda]);
          PhyML_Fprintf(fp_stats,"%g\t",tree->mcmc->tune_move[tree->mcmc->num_move_phyrex_rad]);
          PhyML_Fprintf(fp_stats,"%g\t",tree->mcmc->tune_move[tree->mcmc->num_move_phyrex_mu]);
          PhyML_Fprintf(fp_stats,"%g\t",tree->mcmc->tune_move[tree->mcmc->num_move_phyrex_indel_disk]);
          PhyML_Fprintf(fp_stats,"%g\t",tree->mcmc->tune_move[tree->mcmc->num_move_phyrex_indel_hit]);
          PhyML_Fprintf(fp_stats,"%g\t",tree->mcmc->tune_move[tree->mcmc->num_move_phyrex_ldsk_given_disk]);
          PhyML_Fprintf(fp_stats,"%g\t",tree->mcmc->tune_move[tree->mcmc->num_move_phyrex_indel_disk_serial]);
          PhyML_Fprintf(fp_stats,"%g\t",tree->mcmc->tune_move[tree->mcmc->num_move_phyrex_indel_hit_serial]);
          for(int i=0;i<2*tree->n_otu-2;++i) PhyML_Fprintf(fp_stats,"%g\t",tree->mmod->sigsq_scale[i]);
          
          /* res[0 * tree->mcmc->chain_len / tree->mcmc->sample_interval +  tree->mcmc->run / tree->mcmc->sample_interval] = tree->mmod->lbda;  */
          /* res[1 * tree->mcmc->chain_len / tree->mcmc->sample_interval +  tree->mcmc->run / tree->mcmc->sample_interval] = tree->mmod->mu;  */
          /* res[2 * tree->mcmc->chain_len / tree->mcmc->sample_interval +  tree->mcmc->run / tree->mcmc->sample_interval] = PHYREX_Update_Sigsq(tree);  */
          /* res[3 * tree->mcmc->chain_len / tree->mcmc->sample_interval +  tree->mcmc->run / tree->mcmc->sample_interval] = PHYREX_Neighborhood_Size(tree); */
          /* res[4 * tree->mcmc->chain_len / tree->mcmc->sample_interval +  tree->mcmc->run / tree->mcmc->sample_interval] = tree->mmod->rad; */
          /* res[5 * tree->mcmc->chain_len / tree->mcmc->sample_interval +  tree->mcmc->run / tree->mcmc->sample_interval] = PHYREX_Total_Number_Of_Intervals(tree); */
          /* res[6 * tree->mcmc->chain_len / tree->mcmc->sample_interval +  tree->mcmc->run / tree->mcmc->sample_interval] = PHYREX_Total_Number_Of_Coal_Disks(tree); */
          /* res[7 * tree->mcmc->chain_len / tree->mcmc->sample_interval +  tree->mcmc->run / tree->mcmc->sample_interval] = PHYREX_Total_Number_Of_Hit_Disks(tree); */
          /* res[8 * tree->mcmc->chain_len / tree->mcmc->sample_interval +  tree->mcmc->run / tree->mcmc->sample_interval] = PHYREX_Effective_Density(tree); */
          /* res[9 * tree->mcmc->chain_len / tree->mcmc->sample_interval +  tree->mcmc->run / tree->mcmc->sample_interval] = PHYREX_Coalescence_Rate(tree); */

          /* MCMC_Copy_To_New_Param_Val(tree->mcmc,tree); */
          
          /* for(i=0;i<tree->mcmc->n_moves;i++) if(tree->mcmc->start_ess[i] == YES) MCMC_Update_Effective_Sample_Size(i,tree->mcmc,tree); */
          /* for(i=0;i<tree->mcmc->n_moves;i++) MCMC_Update_Mode(i,tree->mcmc,tree); */
          
          
          PHYREX_Ldsk_To_Tree(tree);
          TIMES_Time_To_Bl(tree);
          tree->bl_ndigits = 3;
          tree->write_tax_names = NO;
          PHYREX_Label_Nodes_With_Locations(tree);
          PHYREX_Label_Edges(tree);
          char *s = Write_Tree(tree);
          PhyML_Fprintf(fp_tree,"\ntree %d [&lnP=%f,precision={1.e-1,1e-02,1e-01}] = [&R] %s",tree->mcmc->sample_num,tree->c_lnL,s);
          PhyML_Fprintf(fp_tree,"\nend;");
          fseek(fp_tree,-5,SEEK_END);
          tree->write_tax_names = YES;
          Free(s);

          fflush(NULL);
                
          tree->mcmc->sample_num++;
        }

      if(!strcmp(tree->mcmc->move_name[move],"phyrex_lbda"))
        MCMC_PHYREX_Lbda(tree);
      
      if(!strcmp(tree->mcmc->move_name[move],"phyrex_mu"))
        MCMC_PHYREX_Mu(tree);
      
      if(!strcmp(tree->mcmc->move_name[move],"phyrex_rad"))
        MCMC_PHYREX_Radius(tree);

      if(!strcmp(tree->mcmc->move_name[move],"phyrex_sigsq"))
        MCMC_PHYREX_Sigsq(tree);

      if(!strcmp(tree->mcmc->move_name[move],"phyrex_neff"))
        MCMC_PHYREX_Neff(tree);

      if(!strcmp(tree->mcmc->move_name[move],"phyrex_indel_disk"))
        MCMC_PHYREX_Indel_Disk(tree);

      if(!strcmp(tree->mcmc->move_name[move],"phyrex_indel_hit"))
        MCMC_PHYREX_Indel_Hit(tree);

      if(!strcmp(tree->mcmc->move_name[move],"phyrex_move_disk_ud"))
        MCMC_PHYREX_Move_Disk_Updown(tree);

      if(!strcmp(tree->mcmc->move_name[move],"phyrex_swap_disk"))
        MCMC_PHYREX_Swap_Disk(tree);

      if(!strcmp(tree->mcmc->move_name[move],"phyrex_scale_times"))
        MCMC_PHYREX_Scale_Times(tree);

      if(!strcmp(tree->mcmc->move_name[move],"phyrex_spr"))
        MCMC_PHYREX_Prune_Regraft(tree);

      if(!strcmp(tree->mcmc->move_name[move],"phyrex_spr_local"))
        MCMC_PHYREX_Prune_Regraft_Local(tree);

      if(!strcmp(tree->mcmc->move_name[move],"phyrex_traj"))
        MCMC_PHYREX_Lineage_Traj(tree);
      
      if(!strcmp(tree->mcmc->move_name[move],"phyrex_disk_multi"))
        MCMC_PHYREX_Disk_Multi(tree);

      if(!strcmp(tree->mcmc->move_name[move],"phyrex_ldsk_multi"))
        MCMC_PHYREX_Ldsk_Multi(tree);

      if(!strcmp(tree->mcmc->move_name[move],"phyrex_ldsk_and_disk"))
        MCMC_PHYREX_Ldsk_And_Disk(tree);

      if(!strcmp(tree->mcmc->move_name[move],"phyrex_ldsk_tip_to_root"))
        MCMC_PHYREX_Ldsk_Tip_To_Root(tree);

      if(!strcmp(tree->mcmc->move_name[move],"phyrex_ldsk_given_disk"))
        MCMC_PHYREX_Ldsk_Given_Disk(tree);

      if(!strcmp(tree->mcmc->move_name[move],"phyrex_disk_given_ldsk"))
        MCMC_PHYREX_Disk_Given_Ldsk(tree);

      if(!strcmp(tree->mcmc->move_name[move],"phyrex_indel_disk_serial"))
        MCMC_PHYREX_Indel_Disk_Serial(tree);
      
      if(!strcmp(tree->mcmc->move_name[move],"phyrex_indel_hit_serial"))
        MCMC_PHYREX_Indel_Hit_Serial(tree);
      
      if(!strcmp(tree->mcmc->move_name[move],"phyrex_add_remove_jump"))
        MCMC_PHYREX_Add_Remove_Jump(tree);

      if(!strcmp(tree->mcmc->move_name[move],"phyrex_sigsq_scale"))
        MCMC_PHYREX_Sigsq_Scale(tree);

      if(!strcmp(tree->mcmc->move_name[move],"kappa"))
        MCMC_Kappa(tree);
      
      if(!strcmp(tree->mcmc->move_name[move],"rr"))
        MCMC_RR(tree);
            
      if(!strcmp(tree->mcmc->move_name[move],"ras"))
        MCMC_Rate_Across_Sites(tree);

      if(!strcmp(tree->mcmc->move_name[move],"br_rate"))
        MCMC_Rates_All(tree);
      
      if(!strcmp(tree->mcmc->move_name[move],"tree_rates"))
        MCMC_Tree_Rates(tree);

      if(!strcmp(tree->mcmc->move_name[move],"clock"))
        MCMC_Clock_R(tree);


      /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
      /* PHYREX_Lk(tree); */
      /* Lk(NULL,tree); */
      
      /* PhyML_Printf("\n. Move: %s tree->mmod->c_lnL: %f kappa: %f", */
      /*              tree->mcmc->move_name[move], */
      /*              tree->mmod->c_lnL, */
      /*              tree->mod->kappa->v); */

      if(tree->mmod->c_lnL < UNLIKELY || tree->c_lnL < UNLIKELY)
        {
          PhyML_Printf("\n. Move: %s tree->mmod->c_lnL: %f tree->c_lnL: %f",
                       tree->mcmc->move_name[move],
                       tree->mmod->c_lnL,tree->c_lnL);
          assert(FALSE);
        }

      if(tree->mmod->safe_phyrex == YES)
        {
          /* phydbl c_lnL = tree->c_lnL; */
          /* Lk(NULL,tree); */
          /* if(Are_Equal(c_lnL,tree->c_lnL,1.E-5) == NO) */
          /*   { */
          /*     PhyML_Fprintf(stderr,"\n. Problem detected with move %s",tree->mcmc->move_name[move]); */
          /*     PhyML_Fprintf(stderr,"\n. c_lnL: %f -> %f",c_lnL,tree->c_lnL); */
          /*     Generic_Exit(__FILE__,__LINE__,__FUNCTION__); */
          /*   } */

          phydbl g_lnL = tree->mmod->c_lnL;
          PHYREX_Lk(tree);
          if(Are_Equal(g_lnL,tree->mmod->c_lnL,1.E-5) == NO)
            {
              PhyML_Fprintf(stderr,"\n. Problem detected with move %s. Iteration %d",tree->mcmc->move_name[move],tree->mcmc->run);
              PhyML_Fprintf(stderr,"\n. g_lnL: %f -> %f [%g]",g_lnL,tree->mmod->c_lnL,g_lnL-tree->mmod->c_lnL);
              Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
            }

          /* phydbl r_lnL = tree->rates->c_lnL_rates; */
          /* RATES_Lk_Rates(tree); */
          /* if(Are_Equal(r_lnL,tree->rates->c_lnL_rates,1.E-5) == NO) */
          /*   { */
          /*     PhyML_Fprintf(stderr,"\n. Problem detected with move %s",tree->mcmc->move_name[move]); */
          /*     PhyML_Fprintf(stderr,"\n. r_lnL: %f -> %f [%g]",r_lnL,tree->rates->c_lnL_rates,r_lnL-tree->rates->c_lnL_rates); */
          /*     Generic_Exit(__FILE__,__LINE__,__FUNCTION__); */
          /*   } */
        }
      
      tree->mcmc->run++;
      MCMC_Get_Acc_Rates(tree->mcmc);
      
      
      
      (void)signal(SIGINT,MCMC_Terminate);
    }
  while(tree->mcmc->run < tree->mcmc->chain_len);

  PhyML_Fprintf(stdout,"\n. The analysis completed !");
  
  fclose(fp_tree);
  fclose(fp_stats);
  /* fclose(fp_summary); */

  return(res);
}
#endif

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl PHYREX_Wrap_Lk(t_edge *b, t_tree *tree, supert_tree *stree)
{
  return PHYREX_Lk(tree);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void PHYREX_Remove_Disk(t_dsk *disk)
{
  t_dsk *prev;
  t_dsk *next;

  prev = disk->prev;
  next = disk->next;

  assert(next != NULL);
  
  if(prev != NULL) prev->next = next;
  next->prev = prev;
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

/* Insert disk event based on its time. Insertion above the root is permitted  
*/
void PHYREX_Insert_Disk(t_dsk *ins, t_tree *tree)
{
  t_dsk *disk;

  assert(ins != NULL);
  
  disk = tree->young_disk;
  while(disk->prev != NULL && disk->prev->time > ins->time) disk = disk->prev;
  
  ins->prev       = disk->prev;
  ins->next       = disk;
  disk->prev      = ins;
  if(ins->prev != NULL) ins->prev->next = ins;
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

t_ldsk *PHYREX_Prev_Coal_Lindisk(t_ldsk *t)
{
  if(t == NULL) return NULL;

  if(t->n_next > 1) 
    {
      return t;
    }
  else
    {
      return PHYREX_Prev_Coal_Lindisk(t->prev);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

t_ldsk *PHYREX_Next_Coal_Lindisk(t_ldsk *t)
{
  assert(!(t == NULL)); 

  if(t->n_next > 1 || t->next == NULL) return t;
  else
    {
      if(t->n_next > 1) // Should have t->is_coal = YES
        {
          PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
          Warn_And_Exit("");
        }
      return PHYREX_Next_Coal_Lindisk(t->next[0]);
    }
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

/*  Generate a new trajectory, including disk event centers, between  
    'y_ldsk' a ``young'' lindisk event and 'o_ldsk' an old one. 'y_ldsk 
    and 'o_ldsk' remain unaffected. No disk events should be present    
    between y_ldsk and o_ldsk: we need to generate some first. n_cur_disk
    is the current number of disks between y_ldsk and o_ldsk.
*/
int PHYREX_One_New_Traj(t_ldsk *y_ldsk, t_ldsk *o_ldsk, int dir_o_y, t_dsk *xtra_dsk, int n_cur_disk, t_tree *tree)
{
  t_phyrex_mod *mmod;
  t_dsk *disk,**disk_new;
  int i,n,K;
  int min_n_disk,n_new_disk,n_disk;

  mmod     = tree->mmod;
  disk     = NULL;
  disk_new = NULL;
  K        = 2;

  /* printf("\n# New traj from %s to %s",y_ldsk->coord->id,o_ldsk->coord->id); */
  /* fflush(NULL); */

  /* Minimum number of disks between y_ldsk and o_ldsk */
  min_n_disk = 0;
  for(i=0;i<mmod->n_dim;i++)
    {
      /* PhyML_Printf("\n# y_ldsk %s : %f o_ldsk->disk->centr: %f rad: %f", */
      /*              y_ldsk->coord->id, */
      /*              y_ldsk->coord->lonlat[i], */
      /*              o_ldsk->disk->centr->lonlat[i], */
      /*              mmod->rad); */
      if(y_ldsk->coord->lonlat[i] < o_ldsk->disk->centr->lonlat[i])
        {
          n_disk = 0;
          while(y_ldsk->coord->lonlat[i] + (2*n_disk-1)*mmod->rad < o_ldsk->disk->centr->lonlat[i] - mmod->rad) n_disk++;
          if(n_disk > min_n_disk) min_n_disk = n_disk;
        }
      else
        {
          n_disk = 0;
          while(y_ldsk->coord->lonlat[i] - (2*n_disk-1)*mmod->rad > o_ldsk->disk->centr->lonlat[i] + mmod->rad) n_disk++;
          if(n_disk > min_n_disk) min_n_disk = n_disk;
        }
      /* printf("  -- min_n_disk: %d",min_n_disk); */
    }
  
  /* printf("\n# min_n_disk: %d cur_n_disk: %d",min_n_disk,n_cur_disk); */
  /* fflush(NULL); */
  
  /* How many disks along the new path between y_ldsk and o_ldsk */
  n_new_disk = Rand_Int(n_cur_disk-K,n_cur_disk+K);
  if(n_new_disk < min_n_disk) n_new_disk = min_n_disk;

  if(xtra_dsk != NULL) n_new_disk++;

  /* printf("\n# Add n_new_disk: %d",n_new_disk); fflush(NULL); */
  
  if(n_new_disk > 0)
    {
      /* Make new disks to create a new path between ldsk_left and ldsk_up */
      disk_new = (t_dsk **)mCalloc(n_new_disk,sizeof(t_dsk *));
      for(i=0;i<n_new_disk-1;i++)  disk_new[i] = PHYREX_Make_Disk_Event(mmod->n_dim,tree->n_otu);
      if(xtra_dsk != NULL) disk_new[n_new_disk-1] = xtra_dsk;
      else                 disk_new[n_new_disk-1] = PHYREX_Make_Disk_Event(mmod->n_dim,tree->n_otu);

      for(i=0;i<n_new_disk;i++)  PHYREX_Init_Disk_Event(disk_new[i],mmod->n_dim,mmod);

      /* Times of these new disks. If xtra_dsk != NULL, then make sure you do not */
      /* reset the time of that disk  */
      n = (xtra_dsk != NULL) ? (n_new_disk-1) : (n_new_disk);
      for(i=0;i<n;i++)
        disk_new[i]->time =
        Uni()*(y_ldsk->disk->time - o_ldsk->disk->time) + o_ldsk->disk->time;       
     
      /* Insert these events */
      for(i=0;i<n_new_disk;i++)
        {
          assert(!tree->young_disk->next);
          disk = tree->young_disk;
          while(disk->time > disk_new[i]->time) disk = disk->prev;
          PHYREX_Insert_Disk(disk_new[i],tree);
        }
            
      /* for(i=0;i<n_new_disk;i++) */
      /*   { */
      /*     printf("\n> disk_new: %f [%s]",disk_new[i]->time,disk_new[i]->id); fflush(NULL); */
      /*   } */
      
      /* Add new lindisks to the new disk events */
      for(i=0;i<n_new_disk;i++)
        {
          disk_new[i]->ldsk = PHYREX_Make_Lindisk_Node(tree->mmod->n_dim);
          PHYREX_Init_Lindisk_Node(disk_new[i]->ldsk,disk_new[i],tree->mmod->n_dim);
          PHYREX_Make_Lindisk_Next(disk_new[i]->ldsk);
          /* printf("\n# Add ldsk %s to %s",disk_new[i]->ldsk->coord->id,disk_new[i]->id); fflush(NULL); */
        }
      
      /* Connect them */
      PHYREX_Connect_Ldsk_Given_Disk(disk_new,n_new_disk,y_ldsk,o_ldsk,dir_o_y);

      Free(disk_new);
    }
  else
    {
      o_ldsk->next[dir_o_y] = y_ldsk;
      y_ldsk->prev          = o_ldsk;
    }
 
  /* Generate a trajectory */
  PHYREX_One_New_Traj_Given_Disk(y_ldsk,o_ldsk,tree);  

  return(n_new_disk);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

/* Generate a new trajectory, including disk event centers, between 
  'y_ldsk' a ``young'' lindisk event and 'o_ldsk' an old one. 'y_ldsk 
  and 'o_ldsk' remain unaffected. Disk events between these two ldsk 
  should already be set. 
*/
void PHYREX_One_New_Traj_Given_Disk(t_ldsk *y_ldsk, t_ldsk *o_ldsk, t_tree *tree)
{
  int n_disk_btw;
  t_ldsk *ldsk;
  phydbl min, max;
  phydbl *min_disk_coord, *max_disk_coord;
  int i,k;
  phydbl rad;


  /* Number of disks between y_ldsk and o_ldsk */
  ldsk = y_ldsk;
  n_disk_btw = 0;
  while(ldsk->prev != o_ldsk)
    {
      n_disk_btw++;
      ldsk = ldsk->prev;
    }

  if(n_disk_btw == 0) return;

  ldsk = y_ldsk;
  rad  = ldsk->disk->mmod->rad;
  k    = 0;
  while(ldsk != o_ldsk)
    {
      if(!ldsk->disk->next) /* Don't change location at tip node */
        {
          ldsk = ldsk->prev;
          continue;
        }
      
      PHYREX_Store_Geo_Coord(ldsk->coord);
      PHYREX_Store_Geo_Coord(ldsk->disk->centr);

      for(i=0;i<tree->mmod->n_dim;i++)
        {
          min = 
            MAX(ldsk->disk->mmod->lim_do->lonlat[i],
                MAX(ldsk->coord->lonlat[i] - 2.*rad,
                    o_ldsk->disk->centr->lonlat[i] - rad*(2.*(n_disk_btw-k)-1.)));

          max = 
            MIN(ldsk->disk->mmod->lim_up->lonlat[i],
                MIN(ldsk->coord->lonlat[i] + 2.*rad,
                    o_ldsk->disk->centr->lonlat[i] + rad*(2.*(n_disk_btw-k)-1.)));
          
          assert(!(max < min));
                
          /* New coordinates for the lindisk */
          ldsk->coord->lonlat[i] = Uni()*(max - min) + min;
        }

      /* New coordinate for the centre of the corresponding disk event */
      PHYREX_Get_Min_Max_Disk_Given_Ldsk(ldsk->disk,&min_disk_coord,&max_disk_coord,tree);
      for(i=0;i<tree->mmod->n_dim;i++) ldsk->disk->centr->lonlat[i] = Uni()*(max_disk_coord[i] - min_disk_coord[i]) + min_disk_coord[i];
      Free(min_disk_coord);
      Free(max_disk_coord);

      ldsk = ldsk->prev;
      k++;
    }
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl PHYREX_Uniform_Path_Density(t_ldsk *y_ldsk, t_ldsk *o_ldsk, t_tree *tree)
{
  int n_disk_btw;
  t_ldsk *ldsk;
  phydbl min, max;
  int i,k;
  phydbl rad;
  phydbl log_dens;
  phydbl *min_disk_coord, *max_disk_coord;

  if(y_ldsk == o_ldsk) return .0;

  /* Number of disks between y_ldsk and o_ldsk */
  ldsk = y_ldsk;
  n_disk_btw = 0;
  while(ldsk->prev != o_ldsk)
    {
      n_disk_btw++;
      ldsk = ldsk->prev;
    }

  if(n_disk_btw == 0) return .0;

  log_dens = 0.0;
  ldsk     = y_ldsk;
  rad      = ldsk->disk->mmod->rad;
  k        = 0;
  while(ldsk != o_ldsk)
    {
      if(!ldsk->disk->next)
        {
          ldsk = ldsk->prev;
          continue;
        }
      
      for(i=0;i<ldsk->disk->mmod->n_dim;i++)
        {
          min = 
            MAX(ldsk->disk->mmod->lim_do->lonlat[i],
                MAX(ldsk->coord->lonlat[i] - 2.*rad,
                    o_ldsk->disk->centr->lonlat[i] - rad*(2.*(n_disk_btw-k)-1.)));

          max = 
            MIN(ldsk->disk->mmod->lim_up->lonlat[i],
                MIN(ldsk->coord->lonlat[i] + 2.*rad,
                    o_ldsk->disk->centr->lonlat[i] + rad*(2.*(n_disk_btw-k)-1.)));

          assert(!(max < min));

          log_dens -= log(max - min);          
        }

      PHYREX_Get_Min_Max_Disk_Given_Ldsk(ldsk->disk,&min_disk_coord,&max_disk_coord,tree);
      for(i=0;i<tree->mmod->n_dim;i++) log_dens -= log(max_disk_coord[i] - min_disk_coord[i]);
      Free(min_disk_coord);
      Free(max_disk_coord);

      ldsk = ldsk->prev;
      k++;
    }

  return(log_dens);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

/* Return the index of the 'next' element of 'old' that should be
   used in order to reach 'young'. 
*/
int PHYREX_Get_Next_Direction(t_ldsk *young, t_ldsk *old)
{
  if(young->disk->time < old->disk->time)
    {
      PhyML_Printf("\n. young (%s) @ time %f; old (%s) @ time %f",
                   young->coord->id,young->disk->time,
                   old->coord->id,old->disk->time);
      assert(FALSE);
    }
  
  assert(!(young == NULL));
 
  if(young->prev == old)
    {
      int i;
      for(i=0;i<old->n_next;i++) if(old->next[i] == young) return i;
    }
  else
    {
      return PHYREX_Get_Next_Direction(young->prev,old);
    }
  return(-1);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void  PHYREX_Update_Lindisk_List_Range(t_dsk *young, t_dsk *old, t_tree *tree)
{
  t_dsk *disk;

  disk = young;
  do
    {
      assert(disk);
      PHYREX_Update_Lindisk_List_Core(disk,tree);
      if(disk == old) break;      
      disk = disk->prev;
    }
  while(1);  
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void PHYREX_Update_Lindisk_List(t_tree *tree)
{
  PHYREX_Update_Lindisk_List_Pre(tree->young_disk->prev,tree);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void PHYREX_Update_Lindisk_List_Pre(t_dsk *disk, t_tree *tree)
{
  if(!disk) return;
  else
    {
      PHYREX_Update_Lindisk_List_Core(disk,tree);
      PHYREX_Update_Lindisk_List_Pre(disk->prev,tree);    
    }
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void PHYREX_Update_Lindisk_List_Core(t_dsk *disk, t_tree *tree)
{
  int i;

  if(!disk) return;
  if(!disk->next) return;
  
  assert(disk->ldsk_a);
  
  // Set ldsk_a[i] to NULL if it does not point to a tip node
  for(i=0;i<tree->n_otu;++i)
    if(disk->ldsk_a[i] &&
       !(disk->ldsk_a[i]->nd != NULL &&
         disk->ldsk_a[i]->nd->tax == YES &&
         disk->age_fixed == YES &&
         disk->ldsk_a[i]->disk == disk))      
      disk->ldsk_a[i] = NULL;
    
  disk->n_ldsk_a = 0;
  for(i=0;i<tree->n_otu;++i)
    if(disk->ldsk_a[i] != NULL)
      disk->n_ldsk_a++;
  
  // Make sure the tip nodes are all at the top of ldsk_a
  for(i=0;i<disk->n_ldsk_a;++i) assert(disk->ldsk_a[i] != NULL);
  for(i=disk->n_ldsk_a;i<tree->n_otu;++i) assert(disk->ldsk_a[i] == NULL);

  for(i=0;i<disk->next->n_ldsk_a;++i)
    {      
      // disk->next->ldsk_a[i] does not coalesce or jump on disk->next
      // --> add it disk->ldsk_a array
      if((disk->next->ldsk_a[i]->prev != NULL &&
          disk->next->ldsk_a[i]->prev != disk->next->ldsk) ||
         (disk->next->ldsk_a[i]->prev == NULL))
        {
          disk->ldsk_a[disk->n_ldsk_a] = disk->next->ldsk_a[i];
          disk->n_ldsk_a++;
        }
    }

  // A jump or coalescence has occurred on disk->next
  // --> add the lineage to disk->ldsk_a array
  if(disk->next->ldsk != NULL)
    {
      disk->ldsk_a[disk->n_ldsk_a] = disk->next->ldsk;
      disk->n_ldsk_a++;
    }

  /* PhyML_Printf("\n"); */
  /* for(i=0;i<disk->n_ldsk_a;++i) PhyML_Printf("\n. Disk %s [%12f] ldsk_a: %s (%3d)",disk->id,disk->time,disk->ldsk_a[i]->coord->id,disk->ldsk_a[i]->nd ? disk->ldsk_a[i]->nd->tax : -1); */
  /* if(disk->ldsk != NULL) PhyML_Printf("\n* Has %s on it (next: %s @ %12f %s)", */
  /*                                     disk->ldsk->coord->id, */
  /*                                     disk->ldsk->next[0]->coord->id, */
  /*                                     disk->ldsk->next[0]->disk->time, */
  /*                                     disk->ldsk->next[0]->prev->coord->id); */
  /* PhyML_Printf("\n. n_ldsk_a: %d",disk->n_ldsk_a); */
  
  if(disk->n_ldsk_a == 0 || disk->n_ldsk_a > tree->n_otu) 
    {
      PhyML_Fprintf(stderr,"\n. disk: %s (%p,%d) time: %f next: %s (%f,%d,%c,%d) prev: %s (%f,%d) disk->n_ldsk_a: %d coord: %s n_otu: %d",
                    disk->id,
                    disk,
                    disk->age_fixed,
                    disk->time,
                    disk->next ? disk->next->id : "??",
                    disk->next ? disk->next->time : 0.0,
                    disk->next ? disk->next->n_ldsk_a : -1,
                    disk->next->ldsk ? 'y' : 'n',
                    disk->next->ldsk ?  disk->next->ldsk->n_next : -1,
                    disk->prev ? disk->prev->id : "??",
                    disk->prev ? disk->prev->time : 0.0,
                    disk->prev ? disk->prev->n_ldsk_a : -1,
                    disk->n_ldsk_a,
                    disk->ldsk?disk->ldsk->coord->id:"??",
                    tree->n_otu);
      assert(FALSE);
    }
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

/* Connect all the ldsk between y_ldsk (young ldsk) and o_ldsk (old ldsk).
   The disk between y_ldsk and o_ldsk should have all been set already
   Note: the disks in **disk are sorted in ascending order of their 
   times
*/

void PHYREX_Connect_Ldsk_Given_Disk(t_dsk **disk, int n_disk, t_ldsk *y_ldsk, t_ldsk *o_ldsk, int dir_o_y)
{
  int i,j;
  t_dsk *disk_tmp;

  /* Sort these events by ascending order of their times */
  for(i=0;i<n_disk-1;i++)
    {
      for(j=i+1;j<n_disk;j++)
        {
          if(disk[j]->time > disk[i]->time)
            {
              disk_tmp = disk[i];
              disk[i]  = disk[j];
              disk[j]  = disk_tmp;
            }
        }
    }


  for(i=0;i<n_disk;i++)
    {
      if(!i)
        {
          disk[i]->ldsk->next[0] = y_ldsk;
          y_ldsk->prev           = disk[i]->ldsk;
          /* printf("\n. connect %s to %s",disk[i]->ldsk->coord->id,y_ldsk->coord->id); fflush(NULL); */
        }
      else
        {
          disk[i]->ldsk->next[0] = disk[i-1]->ldsk;
          disk[i-1]->ldsk->prev  = disk[i]->ldsk;
          /* printf("\n. connect %s to %s",disk[i]->ldsk->coord->id,disk[i-1]->ldsk->coord->id); fflush(NULL); */
        }
    }
  
  /* printf("\n. connect %s next dir: %d [%d] to %s",o_ldsk->coord->id,dir_o_y,o_ldsk->n_next,disk[n_disk-1]->ldsk->coord->id); fflush(NULL); */
  o_ldsk->next[dir_o_y]      = disk[n_disk-1]->ldsk;
  disk[n_disk-1]->ldsk->prev = o_ldsk;

}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void PHYREX_Print_Struct(char sign, t_tree *tree)
{
  t_dsk *disk;
  int i;
  t_ldsk *ldisk;

  PHYREX_Update_Lindisk_List(tree);

  disk = tree->young_disk;
  /* while(disk->prev) disk = disk->prev; */
  do
    {
      PhyML_Printf("\n%c Disk: %s @ %7.3f has %3s on it is_coal? %2d rad: %f age fixed? %d #out: %d",
                   sign,
                   disk->id,
                   disk->time,disk->ldsk?disk->ldsk->coord->id:NULL,
                   disk->ldsk?disk->ldsk->n_next:-1,
                   tree->mmod->rad,
                   disk->age_fixed,
                   PHYREX_Number_Of_Outgoing_Ldsks(disk)); 
      /* for(j=0;j<tree->mmod->n_dim;j++) PhyML_Printf(" %f",disk->centr->lonlat[j]); */
      /* fflush(NULL); */


      for(i=0;i<disk->n_ldsk_a;i++)
        {
          ldisk = disk->ldsk_a[i];

          PhyML_Printf("\n%c ldisk: %s%c prev: %s",
                       sign,
                       ldisk->coord->id,
                       (ldisk->nd && ldisk->nd->tax == YES && ldisk->disk == disk) ? '*' : ' ',
                       ldisk->prev ? ldisk->prev->coord->id : NULL);
          fflush(NULL);
          
          /* for(j=0;j<tree->mmod->n_dim;j++) */
          /*   { */
          /*     PhyML_Printf(" %f",ldisk->coord->lonlat[j]); */
          /*     fflush(NULL); */

          /*     if(FABS(ldisk->coord->lonlat[j] - ldisk->disk->centr->lonlat[j]) > 2.*tree->mmod->rad && */
          /*        ldisk->disk->ldsk == ldisk) PhyML_Printf(" ! "); */

          /*     if(ldisk->prev) */
          /*       { */
          /*         if(ldisk->coord->lonlat[j] < ldisk->prev->disk->centr->lonlat[j] - tree->mmod->rad) PhyML_Printf(" #a "); */
          /*         if(ldisk->coord->lonlat[j] > ldisk->prev->disk->centr->lonlat[j] + tree->mmod->rad) PhyML_Printf(" #b "); */
          /*       } */

          /*     if(ldisk->next) */
          /*       { */
          /*         if(ldisk->coord->lonlat[j] < ldisk->disk->centr->lonlat[j] - tree->mmod->rad) PhyML_Printf(" $a "); */
          /*         if(ldisk->coord->lonlat[j] > ldisk->disk->centr->lonlat[j] + tree->mmod->rad) PhyML_Printf(" $b "); */
          /*       } */
    }


      disk = disk->prev;      
    }
  while(disk);
  PhyML_Printf("\n. End of struct");
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

int PHYREX_Check_Struct(t_tree *tree)
{
  int i;
  t_ldsk *ldisk;
  t_dsk *disk;

  disk = tree->young_disk;

  do
    {
      if(disk->n_ldsk_a == 0)
        {
          PhyML_Printf("\n. disk %s time %f has n_ldsk_a=0 (next has %d prev has %d)",
                       disk->id,
                       disk->time,
                       disk->next ? disk->next->n_ldsk_a : -1,
                       disk->prev ? disk->prev->n_ldsk_a : -1);
          assert(FALSE);
          return 0;
        }

      disk = disk->prev;
    }
  while(disk);
  
  // Check times
  for(i=0;i<tree->n_otu;++i)
    {
      ldisk = tree->a_nodes[i]->ldsk;
      assert(ldisk);
      do
        {          
          if(ldisk->prev->disk->time > ldisk->disk->time)
            {
              PhyML_Printf("\n. ldisk->id: %s ldisk->prev->id: %s ldsk->disk->time: %f  ldsk->prev->disk->time: %f diff: %g ldisk->prev->disk: %s ldisk->disk: %s",
                           ldisk->coord->id,
                           ldisk->prev->coord->id,
                           ldisk->disk->time,
                           ldisk->prev->disk->time,
                           ldisk->prev->disk->time-ldisk->disk->time,
                           ldisk->prev->disk->id,
                           ldisk->disk->id);
              assert(FALSE);
              return 0;
            }
          ldisk = ldisk->prev;
        }
      while(ldisk->prev);
    }
  return 1;
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void PHYREX_Store_Geo_Coord(t_geo_coord *t)
{
  int i;

  for(i=0;i<t->dim;i++) t->cpy->lonlat[i] = t->lonlat[i];  
  t->cpy->dim = t->dim;
  strcpy(t->cpy->id,t->id);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/


void PHYREX_Restore_Geo_Coord(t_geo_coord *t)
{
  int i;

  for(i=0;i<t->dim;i++) t->lonlat[i] = t->cpy->lonlat[i];  
  t->dim = t->cpy->dim;
  strcpy(t->id,t->cpy->id);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

int PHYREX_Total_Number_Of_Intervals(t_tree *tree)
{
  t_dsk *disk;
  int n_intervals;

  assert(!(tree->young_disk->next));

  disk = tree->young_disk;
  n_intervals = 0;
  while(disk->prev)
    {
      n_intervals++;
      disk = disk->prev;
    }
  return(n_intervals);

}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

int PHYREX_Number_Of_Intervals_Range(t_dsk *young, t_dsk *old, t_tree *tree)
{
  t_dsk *disk;
  int n_intervals;

  disk = young;
  n_intervals = 0;
  do
    {
      n_intervals++;
      disk = disk->prev;
      assert(disk);
    }
  while(disk != old);
  
  return(n_intervals);
}


/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

int PHYREX_Total_Number_Of_Hit_Disks(t_tree *tree)
{
  t_dsk *disk;
  int n_hit_disks;

  assert(!(tree->young_disk->next));

  disk = tree->young_disk;
  n_hit_disks = 0;
  while(disk)
    {
      if(disk->ldsk) n_hit_disks++;
      disk = disk->prev;
    }
  return(n_hit_disks);

}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

int PHYREX_Total_Number_Of_Coal_Disks(t_tree *tree)
{
  t_dsk *disk;
  int n_coal_disks;

  assert(!(tree->young_disk->next));

  disk = tree->young_disk;
  n_coal_disks = 0;
  while(disk)
    {
      if(disk->ldsk && disk->ldsk->n_next > 1) n_coal_disks++;
      disk = disk->prev;
    }

  return(n_coal_disks);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl PHYREX_Log_Dunif_Rectangle_Overlap(t_ldsk *ldsk, t_dsk *disk, t_phyrex_mod *mmod)
{
  phydbl up, down, left, rght;
  phydbl log_dens;

  log_dens  = 0.0;

  up   = MIN(disk->centr->lonlat[0] + mmod->rad, mmod->lim_up->lonlat[0]);
  down = MAX(disk->centr->lonlat[0] - mmod->rad, mmod->lim_do->lonlat[0]);
  rght = MIN(disk->centr->lonlat[1] + mmod->rad, mmod->lim_up->lonlat[1]);
  left = MAX(disk->centr->lonlat[1] - mmod->rad, mmod->lim_do->lonlat[1]);


  if(ldsk->coord->lonlat[0] < down)  return UNLIKELY;
  if(ldsk->coord->lonlat[0] > up)    return UNLIKELY;
  if(ldsk->coord->lonlat[1] < left)  return UNLIKELY;
  if(ldsk->coord->lonlat[1] > rght)  return UNLIKELY;

  log_dens = -log(up-down)-log(rght-left);

  return(log_dens);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
/* Samples uniformly within a rectangle (with truncation for border) */
/* and returns the corresponding log density */
phydbl PHYREX_Runif_Rectangle_Overlap(t_ldsk *ldsk, t_dsk *disk, t_phyrex_mod *mmod)
{
  phydbl up, down, left, rght;

  up   = MIN(disk->centr->lonlat[0] + mmod->rad, mmod->lim_up->lonlat[0]);
  down = MAX(disk->centr->lonlat[0] - mmod->rad, mmod->lim_do->lonlat[0]);
  rght = MIN(disk->centr->lonlat[1] + mmod->rad, mmod->lim_up->lonlat[1]);
  left = MAX(disk->centr->lonlat[1] - mmod->rad, mmod->lim_do->lonlat[1]);

  ldsk->coord->lonlat[0] = Uni()*(up - down) + down;
  ldsk->coord->lonlat[1] = Uni()*(rght - left) + left;

  /* printf("\n. disk %s (%f %f) rad: %f up: %f down: %f rght: %f left: %f", */
  /*        disk->id, */
  /*        disk->centr->lonlat[0], */
  /*        disk->centr->lonlat[1], */
  /*        mmod->rad, */
  /*        up,down,rght,left); */
  return(log(up-down)+log(rght-left));
}
 
/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
 
phydbl PHYREX_Rnorm_Trunc(t_ldsk *ldsk, t_dsk *disk, t_phyrex_mod *mmod)
{
  int i,err;

  err = NO;
  for(i=0;i<mmod->n_dim;++i)
    {
      ldsk->coord->lonlat[i] = Rnorm_Trunc(disk->centr->lonlat[i],
                                           mmod->rad,
                                           mmod->lim_do->lonlat[i],
                                           mmod->lim_up->lonlat[i],&err);
      assert(err != YES);
    }
  

  return(0.0);
}
 
/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl PHYREX_Wrap_Prior_Radius(t_edge *e, t_tree *tree, supert_tree *st)
{
  return PHYREX_LnPrior_Radius(tree);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl PHYREX_LnPrior_Lbda(t_tree *tree)
{
  if(tree->mmod->id == RW)  return(0.0);
  if(tree->mmod->id == RRW) return(0.0);
  
  if(tree->mmod->lbda < tree->mmod->min_lbda) return UNLIKELY;
  if(tree->mmod->lbda > tree->mmod->max_lbda) return UNLIKELY;

  /* tree->mmod->c_ln_prior_lbda = */
  /*   log(tree->mmod->prior_param_lbda) - */
  /*   tree->mmod->prior_param_lbda*tree->mmod->lbda; */

  /* tree->mmod->c_ln_prior_lbda -= log(exp(-tree->mmod->prior_param_lbda*tree->mmod->min_lbda)- */
  /*                                    exp(-tree->mmod->prior_param_lbda*tree->mmod->max_lbda)); */

  tree->mmod->c_ln_prior_lbda = -log(tree->mmod->max_lbda - tree->mmod->min_lbda);;

  return(tree->mmod->c_ln_prior_lbda);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl PHYREX_LnPrior_Mu(t_tree *tree)
{
  if(tree->mmod->id == RW)  return(0.0);
  if(tree->mmod->id == RRW) return(0.0);

  if(tree->mmod->mu < tree->mmod->min_mu) return UNLIKELY;
  if(tree->mmod->mu > tree->mmod->max_mu) return UNLIKELY;

  tree->mmod->c_ln_prior_mu = -log(tree->mmod->max_mu - tree->mmod->min_mu);

  /* tree->mmod->c_ln_prior_mu = -2.*log(tree->mmod->mu); */

  return(tree->mmod->c_ln_prior_mu);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl PHYREX_LnPrior_Radius(t_tree *tree)
{
  if(tree->mmod->id == RW)  return(0.0);
  if(tree->mmod->id == RRW) return(0.0);

  if(tree->mmod->rad < tree->mmod->min_rad) return UNLIKELY;
  if(tree->mmod->rad > tree->mmod->max_rad) return UNLIKELY;

  tree->mmod->c_ln_prior_rad =
    log(tree->mmod->prior_param_rad) -
    tree->mmod->prior_param_rad*tree->mmod->rad;

  tree->mmod->c_ln_prior_rad -= log(exp(-tree->mmod->prior_param_rad*tree->mmod->min_rad)-
                                    exp(-tree->mmod->prior_param_rad*tree->mmod->max_rad));

  return(tree->mmod->c_ln_prior_rad);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl PHYREX_LnPrior_Sigsq(t_tree *tree)
{
  return(0.0);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void PHYREX_Initial_Ldsk_Pos(t_tree *tree)
{
  t_dsk *disk;
  int i,j;
  phydbl mean;

  disk = tree->young_disk->prev;

  do
    {
      if(disk->ldsk)
        {
          for(i=0;i<tree->mmod->n_dim;i++)
            {
              mean = 0.0;
              for(j=0;j<disk->ldsk->n_next;j++) mean += disk->ldsk->next[j]->coord->lonlat[i];
              disk->ldsk->coord->lonlat[i] = mean / (phydbl)disk->ldsk->n_next;
            }
        }
      disk = disk->prev;
    }
  while(disk);


}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl PHYREX_Min_Radius(t_tree *tree)
{
  phydbl ori_rad, min_rad;

  ori_rad = tree->mmod->rad;
  tree->mmod->rad = tree->mmod->max_rad;
  do
    {
      PHYREX_Lk(tree);
      tree->mmod->rad -= 1.0;
    }
  while(tree->mmod->c_lnL > UNLIKELY + 0.1);

  min_rad = tree->mmod->rad + 2.0;
  tree->mmod->rad = ori_rad;
  PHYREX_Lk(tree);
  return(min_rad);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
/* Get the minimum and maximum values a ldsk can take, given the position of the disk centre */
void PHYREX_Get_Min_Max_Ldsk_Given_Disk(t_ldsk *ldsk, phydbl **min, phydbl **max, t_tree *tree)
{
  phydbl *loc_min,*loc_max;
  int i;

  if(!ldsk->disk->next) return;

  loc_min = (phydbl *)mCalloc(tree->mmod->n_dim, sizeof(phydbl));
  loc_max = (phydbl *)mCalloc(tree->mmod->n_dim, sizeof(phydbl));

  for(i=0;i<tree->mmod->n_dim;i++)
    {
      loc_min[i] = ldsk->disk->centr->lonlat[i] - tree->mmod->rad;
      loc_max[i] = ldsk->disk->centr->lonlat[i] + tree->mmod->rad;     
      
      if(ldsk->prev)
        {
          loc_min[i] = MAX(loc_min[i],ldsk->prev->disk->centr->lonlat[i] - tree->mmod->rad);
          loc_max[i] = MIN(loc_max[i],ldsk->prev->disk->centr->lonlat[i] + tree->mmod->rad);
        }

      loc_min[i] = MAX(tree->mmod->lim_do->lonlat[i],loc_min[i]);
      loc_max[i] = MIN(tree->mmod->lim_up->lonlat[i],loc_max[i]);
    }

  (*min) = loc_min;
  (*max) = loc_max;
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
/* Get the minimum and maximum values a disk centre can take, given the position of the ldsk */
void PHYREX_Get_Min_Max_Disk_Given_Ldsk(t_dsk *disk, phydbl **min, phydbl **max, t_tree *tree)
{
  phydbl *loc_min,*loc_max;
  int i,j;
  phydbl tmp_min, tmp_max;

  if(!disk->next) return;

  loc_min = (phydbl *)mCalloc(tree->mmod->n_dim,sizeof(phydbl));
  loc_max = (phydbl *)mCalloc(tree->mmod->n_dim,sizeof(phydbl));

  if(!disk->ldsk || tree->mmod->id == SLFV_GAUSSIAN)
    {
      for(i=0;i<tree->mmod->n_dim;i++)
        {
          loc_min[i] = tree->mmod->lim_do->lonlat[i];
          loc_max[i] = tree->mmod->lim_up->lonlat[i];
        }
    }
  else
    {
      for(i=0;i<tree->mmod->n_dim;i++)
        {
          tmp_min = +INFINITY;
          tmp_max = -INFINITY;
          for(j=0;j<disk->ldsk->n_next;j++)
            {
              if(disk->ldsk->next[j]->coord->lonlat[i] < tmp_min) tmp_min = disk->ldsk->next[j]->coord->lonlat[i];
              if(disk->ldsk->next[j]->coord->lonlat[i] > tmp_max) tmp_max = disk->ldsk->next[j]->coord->lonlat[i];
            }

          if(disk->ldsk->coord->lonlat[i] < tmp_min) tmp_min = disk->ldsk->coord->lonlat[i];
          if(disk->ldsk->coord->lonlat[i] > tmp_max) tmp_max = disk->ldsk->coord->lonlat[i];

          loc_min[i] = MAX(tree->mmod->lim_do->lonlat[i],
                           tmp_max - tree->mmod->rad);
          
          loc_max[i] = MIN(tree->mmod->lim_up->lonlat[i],
                           tmp_min + tree->mmod->rad);
          
          assert(!(loc_max[i] < loc_min[i]));
        }
    }

  (*min) = loc_min;
  (*max) = loc_max;
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void PHYREX_Update_Disk_Ldsk_Subtree(t_ldsk *root_ldsk, t_tree *tree)
{
  int i;  
  for(i=0;i<root_ldsk->n_next;i++) PHYREX_Update_Disk_Ldsk_Subtree_Pre(root_ldsk,root_ldsk->next[i],root_ldsk,tree);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void PHYREX_Update_Disk_Ldsk_Subtree_Pre(t_ldsk *old_ldsk, t_ldsk *young_ldsk, t_ldsk *root_ldsk, t_tree *tree)
{  
  if(!young_ldsk->disk->next) 
    {
      PHYREX_One_New_Traj_Given_Disk(young_ldsk,root_ldsk,tree);
      return;
    }
  else
    {
      int i;
      PHYREX_Update_Disk_Ldsk_Subtree_Pre(young_ldsk,young_ldsk->next[0],root_ldsk,tree);
      if(young_ldsk->n_next > 1) 
        {
          for(i=1;i<young_ldsk->n_next;i++) PHYREX_Update_Disk_Ldsk_Subtree_Pre(young_ldsk,young_ldsk->next[i],young_ldsk,tree);
        }
    }
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void PHYREX_Restore_Disk_Ldsk_Subtree(t_ldsk *root_ldsk, t_tree *tree)
{
  int i;  
  for(i=0;i<root_ldsk->n_next;i++) PHYREX_Restore_Disk_Ldsk_Subtree_Pre(root_ldsk,root_ldsk->next[i],tree);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void PHYREX_Restore_Disk_Ldsk_Subtree_Pre(t_ldsk *old_ldsk, t_ldsk *young_ldsk, t_tree *tree)
{  

  if(!young_ldsk->disk->next) return;
  else
    {
      int i;

      PHYREX_Restore_Geo_Coord(young_ldsk->coord);
      PHYREX_Restore_Geo_Coord(young_ldsk->disk->centr);

      for(i=0;i<young_ldsk->n_next;i++)
        {
          PHYREX_Restore_Disk_Ldsk_Subtree_Pre(young_ldsk,young_ldsk->next[i],tree);
        }
    }
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void PHYREX_Proposal_Disk_Ldsk_Subtree(t_ldsk *root_ldsk, phydbl *logdens, t_tree *tree)
{
  int i;  
  (*logdens) = 0.0;
  for(i=0;i<root_ldsk->n_next;i++) PHYREX_Proposal_Disk_Ldsk_Subtree_Pre(root_ldsk,root_ldsk->next[i],root_ldsk,logdens,tree);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void PHYREX_Proposal_Disk_Ldsk_Subtree_Pre(t_ldsk *old_ldsk, t_ldsk *young_ldsk, t_ldsk *root_ldsk, phydbl *logdens, t_tree *tree)
{  

  if(!young_ldsk->disk->next) 
    {
      (*logdens) += PHYREX_Uniform_Path_Density(young_ldsk,root_ldsk,tree);  
      return;
    }
  else
    {
      int i;
      for(i=0;i<young_ldsk->n_next;i++)
        {
          PHYREX_Proposal_Disk_Ldsk_Subtree_Pre(young_ldsk,young_ldsk->next[i],root_ldsk,logdens,tree);
        }
    }
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
/* Update the tree structure given the whole set of ldsk events */
/* Coalescent events involving multiple lineages are resolved using */
/* very short internal edges. Tip nodes in the tree are always connected */
/* to the corresponding ldsks. */
void PHYREX_Ldsk_To_Tree(t_tree *tree)
{
  int i,j;
  t_dsk *disk;
  t_ldsk *root_ldsk;
  
  /* Reset */
  for(i=0;i<2*tree->n_otu-1;++i) 
    {
      for(j=0;j<3;++j) 
        {
          tree->a_nodes[i]->v[j] = NULL;
          tree->a_nodes[i]->b[j] = NULL;
        }
    }

  // Erase all connections to internal nodes
  disk = tree->young_disk;
  do
    {
      if(disk->ldsk)
        {
          disk->ldsk->nd = NULL;
          assert(disk->age_fixed == NO);
        }
      disk = disk->prev;
    }
  while(disk->prev);
  
  
  // Make sure oldest disk has a ldsk on it
  assert(disk->ldsk);
  root_ldsk = disk->ldsk;
  
  if(tree->n_root == NULL) tree->n_root = tree->a_nodes[2*tree->n_otu-2];
  assert(tree->n_root);
  
  i = 2*tree->n_otu-3;
  tree->num_curr_branch_available = 0;
  PHYREX_Ldsk_To_Tree_Post(tree->n_root,root_ldsk,&i,tree);


  for(i=0;i<tree->n_otu;++i) assert(tree->a_nodes[i]->v[0]);

  for(i=0;i<3;++i) 
    if(tree->n_root->v[2]->v[i] == tree->n_root) 
      { 
        tree->n_root->v[2]->v[i] = tree->n_root->v[1]; 
        break; 
      }

  for(i=0;i<3;++i) 
    if(tree->n_root->v[1]->v[i] == tree->n_root) 
      { 
        tree->n_root->v[1]->v[i] = tree->n_root->v[2]; 
        break; 
      }

  Connect_Edges_To_Nodes_Serial(tree);


  tree->e_root = NULL;
  for(i=0;i<2*tree->n_otu-3;++i)
    {
      if((tree->a_edges[i]->left == tree->n_root->v[1] && tree->a_edges[i]->rght == tree->n_root->v[2]) ||
         (tree->a_edges[i]->left == tree->n_root->v[2] && tree->a_edges[i]->rght == tree->n_root->v[1]))
        {        
          tree->e_root = tree->a_edges[i];
          break;
        }
    }
  assert(!(tree->e_root == NULL));
  
  tree->n_root->b[1] = tree->a_edges[2*tree->n_otu-3];
  tree->n_root->b[2] = tree->a_edges[2*tree->n_otu-2];

  tree->n_root->b[1]->left = tree->n_root;
  tree->n_root->b[1]->rght = tree->n_root->v[1];

  tree->n_root->b[2]->left = tree->n_root;
  tree->n_root->b[2]->rght = tree->n_root->v[2];
  
  Update_Ancestors(tree->n_root,tree->n_root->v[2],tree);
  Update_Ancestors(tree->n_root,tree->n_root->v[1],tree);
  
  MIXT_Propagate_Tree_Update(tree);  
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void PHYREX_Ldsk_To_Tree_Post(t_node *a, t_ldsk *ldsk, int *available, t_tree *tree)
{
  assert(ldsk);
  assert(a);

  ldsk->nd = a;  
  a->ldsk = ldsk;
  tree->times->nd_t[a->num] = ldsk->disk->time;

  if(!ldsk->next) return; // Tip node
  else
    {
      t_node *parent,*son;
      int idx_next;
      t_ldsk *t;

      parent       = a;
      parent->v[1] = NULL;
      parent->v[2] = NULL;
      idx_next     = 0;
      do
        {
          t = ldsk->next[idx_next]; 

          // Descend along that lineage as long as
          // one has not reached a tip (t->next == NULL)
          // or a coalescent event (t->n_next > 1)
          while(t->next && t->n_next <= 1) t = t->next[0];
          
          
          if(t->nd == NULL) // t->disk is not a sampled disk
            {
              assert(t->disk->age_fixed == NO);
              son = tree->a_nodes[*available];
              (*available) = (*available)-1;
            }
          else
            {
              assert(t->disk->age_fixed == YES);
              son = t->nd;
            }
                    
          PHYREX_Ldsk_To_Tree_Post(son,t,available,tree);          

          // Resolve multifurcation
          if(parent->v[2] != NULL && idx_next >= 2) 
            {
              t_node *new_parent;

              new_parent = tree->a_nodes[*available];
              new_parent->ldsk = ldsk;
              
              (*available) = (*available)-1;              

              new_parent->v[0]       = parent;
              new_parent->v[1]       = parent->v[2];
              new_parent->v[2]       = son;

              parent->v[2]           = new_parent;
              son->v[0]              = new_parent;
              new_parent->v[1]->v[0] = new_parent;
              
              /* PhyML_Printf("\n[] connect %d to %d",parent->num,new_parent->num); */
              /* PhyML_Printf("\n[] connect %d to %d",new_parent->num,new_parent->v[1]->num); */
              /* PhyML_Printf("\n[] connect %d to %d",new_parent->num,new_parent->v[2]->num); */
              
              tree->times->nd_t[new_parent->num] = ldsk->disk->time;

              parent = new_parent;
            }
          else
            {
              son->v[0] = parent;          
              if(!parent->v[1]) parent->v[1] = son;
              else              parent->v[2] = son;
              
              /* printf("\n[] connect %d to %d",parent->num,son->num); */
            }

          idx_next++;
        }
      while(idx_next != ldsk->n_next);
    }
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
// Make sure PHYREX_Make_And_Connect_Tip_Disks was called beforehand
void PHYREX_Tree_To_Ldsk(t_tree *tree)
{
  t_dsk *a_disk,*disk;
  t_node *n;
  
  assert(tree->n_root);
  assert(tree->young_disk);
  
  // Initialise root disk
  a_disk = PHYREX_Make_Disk_Event(tree->mmod->n_dim,tree->n_otu);
  PHYREX_Init_Disk_Event(a_disk,tree->mmod->n_dim,NULL);

  a_disk->prev = NULL; // last (i.e., oldest) disk
  
  a_disk->ldsk = PHYREX_Make_Lindisk_Node(tree->mmod->n_dim);
  PHYREX_Init_Lindisk_Node(a_disk->ldsk,a_disk,tree->mmod->n_dim);

  tree->n_root->ldsk = a_disk->ldsk;
  a_disk->ldsk->nd = tree->n_root;
  
  // Initialize centre of event on the root disk
  a_disk->centr->lonlat[0] = Uni()*(tree->mmod->lim_up->lonlat[0]-tree->mmod->lim_do->lonlat[0])+tree->mmod->lim_do->lonlat[0];
  a_disk->centr->lonlat[1] = Uni()*(tree->mmod->lim_up->lonlat[1]-tree->mmod->lim_do->lonlat[1])+tree->mmod->lim_do->lonlat[1];      
  
  /* Its location */
  switch(tree->mmod->id)
    {
    case SLFV_UNIFORM:
      {
        PHYREX_Runif_Rectangle_Overlap(a_disk->ldsk,a_disk,tree->mmod);
        break;
      }
    case SLFV_GAUSSIAN:
      {
        PHYREX_Rnorm_Trunc(a_disk->ldsk,a_disk,tree->mmod);
        break;
      }
    }
  
  a_disk->ldsk->nd = tree->n_root;
  a_disk->time = tree->times->nd_t[tree->n_root->num];

  
  /* Inflate_Times_To_Get_Reasonnable_Edge_Lengths(1.E-3,tree); */
  Get_Node_Ranks_From_Times(tree);

  PHYREX_Tree_To_Ldsk_Post(tree->n_root,tree->n_root->v[1],a_disk,tree);
  PHYREX_Tree_To_Ldsk_Post(tree->n_root,tree->n_root->v[2],a_disk,tree);

  // Create a doubly-chained list of disks
  disk = a_disk;
  disk->time = tree->times->nd_t[tree->n_root->num];
  n = tree->n_root;
  while(n->rk_next)
    {
      // Only jump to next disk if n and n->rk_next are on distinct disks
      // If n->ldsk == n->rk_next->ldsk then n->ldsk->disk and n->rk_next->ldsk->disk
      // are the same disk (multiple merger)
      if((n->ldsk->disk != n->rk_next->ldsk->disk) && (n->ldsk != n->rk_next->ldsk))
        {
          disk->next = n->rk_next->ldsk->disk;          
          disk->next->prev = disk;
          disk->next->time = tree->times->nd_t[n->rk_next->num];
          disk = disk->next;
        }
      
      n = n->rk_next;
    }

  
  // Fill in ldsk_a arrays throughout the tree
  PHYREX_Update_Lindisk_List(tree);

  assert(tree->young_disk->prev);

}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void PHYREX_Tree_To_Ldsk_Post(t_node *a, t_node *d, t_dsk *a_disk, t_tree *tree)
{
  int i;
  
  assert(a);
  assert(d);
  assert(a_disk);
  
  
  if(d->tax)
    {
      assert(d->ldsk);
      PHYREX_Make_Lindisk_Next(a_disk->ldsk);
      d->ldsk->prev = a_disk->ldsk;
      a_disk->ldsk->next[a_disk->ldsk->n_next-1] = d->ldsk;
      a_disk->ldsk->next[a_disk->ldsk->n_next-1]->nd = d;      
      return;
    }
  else
    {
      if(tree->times->nd_t[d->num] > tree->times->nd_t[a->num])
        {
          t_dsk *d_disk;
          
          PHYREX_Make_Lindisk_Next(a_disk->ldsk);

          // Make and initialize descendent disk
          d_disk = PHYREX_Make_Disk_Event(tree->mmod->n_dim,tree->n_otu);
          assert(d_disk);
          PHYREX_Init_Disk_Event(d_disk,tree->mmod->n_dim,NULL);
          
          d_disk->time = tree->times->nd_t[d->num];
          
          d_disk->ldsk = PHYREX_Make_Lindisk_Node(tree->mmod->n_dim);
          PHYREX_Init_Lindisk_Node(d_disk->ldsk,d_disk,tree->mmod->n_dim);
          
          // Initialize centre of event on the root disk
          d_disk->centr->lonlat[0] = Uni()*(tree->mmod->lim_up->lonlat[0]-tree->mmod->lim_do->lonlat[0])+tree->mmod->lim_do->lonlat[0];
          d_disk->centr->lonlat[1] = Uni()*(tree->mmod->lim_up->lonlat[1]-tree->mmod->lim_do->lonlat[1])+tree->mmod->lim_do->lonlat[1];      
          
          /* Its location */
          switch(tree->mmod->id)
            {
            case SLFV_UNIFORM:
              {
                PHYREX_Runif_Rectangle_Overlap(d_disk->ldsk,d_disk,tree->mmod);
                break;
              }
            case SLFV_GAUSSIAN:
              {
                PHYREX_Rnorm_Trunc(d_disk->ldsk,d_disk,tree->mmod);
                break;
              }
            }
          
          d_disk->ldsk->nd = d;
          d->ldsk = d_disk->ldsk;
          
          a_disk->ldsk->next[a_disk->ldsk->n_next-1] = d_disk->ldsk;
          d_disk->ldsk->prev = a_disk->ldsk;
          
          for(i=0;i<3;++i)
            {
              if(d->v[i] != a && d->b[i] != tree->e_root)
                {
                  PHYREX_Tree_To_Ldsk_Post(d,d->v[i],d_disk,tree);
                }
            }

        }
      else // time[d] == time[a] -> add lineage to a_disk instead of creating d_disk
        {
          d->ldsk = a_disk->ldsk;

          for(i=0;i<3;++i)
            {
              if(d->v[i] != a && d->b[i] != tree->e_root)
                {
                  PHYREX_Tree_To_Ldsk_Post(d,d->v[i],a_disk,tree);
                }
            }
        } 
    }
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void PHYREX_Simulate_Disk_And_Node_Times(t_tree *tree)
{
  t_dsk *disk;

  disk = tree->young_disk;
  assert(disk->age_fixed == YES); // age of youngest disk should be fixed
  
  do
    {
      // disk->prev is not a sample --> simulate its age by sampling in exp distribution
      if(disk->prev->age_fixed == NO)
        {
          disk->prev->time = disk->time - Rexp(tree->mmod->lbda);

          // set time of internal node sitting on disk->prev
          if(disk->prev->ldsk != NULL)
            tree->times->nd_t[disk->prev->ldsk->nd->num] =
              disk->prev->time;
        }
      PhyML_Printf("\n. Simulate times disk %s time: %f",disk->id,disk->time);
      disk = disk->prev;
    }
  while(disk->prev);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
// Connections between tip nodes and corresponding ldsks (plus the
// associated disks) should be made early on and never modified
// after that.
void PHYREX_Make_And_Connect_Tip_Disks(t_tree *tree)
{
  t_dsk *disk,*new_disk;
  t_node *n;
  
  Get_Node_Ranks_From_Tip_Times(tree);

  // Find most recent tip
  n = tree->a_nodes[0];
  while(n->rk_next) n = n->rk_next;
  
  // Create most recent disk
  disk = PHYREX_Make_Disk_Event(tree->mmod->n_dim,tree->n_otu);
  PHYREX_Init_Disk_Event(disk,tree->mmod->n_dim,NULL);
  disk->age_fixed = YES; // Sample disks have their ages fixed
  disk->time = tree->times->nd_t[n->num]; // Set time of youngest disk
  PhyML_Fprintf(stdout,"\n. Youngest sampled disk time set to %f (disk id: %s)",disk->time,disk->id);
  
  // ldsk_a[0] is connected to youngest tip
  disk->ldsk_a[0] = PHYREX_Make_Lindisk_Node(tree->mmod->n_dim);
  PHYREX_Init_Lindisk_Node(disk->ldsk_a[0],disk,tree->mmod->n_dim);
  disk->n_ldsk_a = 1;
  disk->ldsk_a[0]->nd = n;
  n->ldsk = disk->ldsk_a[0];
  n->ldsk->disk = disk;

  
  // Set pointer to young_disk here and not elsewhere!
  tree->young_disk = disk;
  new_disk = NULL;
  do
    {
      assert(n->rk_prev);
      
      // n and n->rk_prev have distinct time stamps -> they should be on two distinct disks
      if(Are_Equal(tree->times->nd_t[n->num],tree->times->nd_t[n->rk_prev->num],SMALL) == NO)
        {
          new_disk = PHYREX_Make_Disk_Event(tree->mmod->n_dim,tree->n_otu);
          PHYREX_Init_Disk_Event(new_disk,tree->mmod->n_dim,NULL);
          new_disk->age_fixed = YES;
          
          new_disk->ldsk_a[0] = PHYREX_Make_Lindisk_Node(tree->mmod->n_dim);
          PHYREX_Init_Lindisk_Node(new_disk->ldsk_a[0],new_disk,tree->mmod->n_dim);
          new_disk->ldsk_a[0]->nd = n->rk_prev;
          new_disk->ldsk_a[0]->disk = new_disk;
          new_disk->n_ldsk_a = 1;

          new_disk->next = disk;
          disk->prev = new_disk;

          disk = new_disk;
        }
      // n and n->rk_prev have the same time stamp -> they sit on the same disk
      else
        {
          disk->ldsk_a[disk->n_ldsk_a] = PHYREX_Make_Lindisk_Node(tree->mmod->n_dim);
          PHYREX_Init_Lindisk_Node(disk->ldsk_a[disk->n_ldsk_a],disk,tree->mmod->n_dim);
          disk->ldsk_a[disk->n_ldsk_a]->nd = n->rk_prev;
          disk->ldsk_a[disk->n_ldsk_a]->disk = disk;
          disk->n_ldsk_a++;
        }

      // Set sampled disk time
      disk->time = tree->times->nd_t[n->rk_prev->num];
      PhyML_Fprintf(stdout,"\n. Set sampled disk (id: %s) time to %15f",disk->id,disk->time);
      
      n->rk_prev->ldsk = disk->ldsk_a[disk->n_ldsk_a-1];
      n = n->rk_prev;
    }
  while(n->rk_prev);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void PHYREX_Remove_Lindisk_Next(t_ldsk *ldsk, t_ldsk *rm)
{
  t_ldsk **new_next;
  int i,pos;

  new_next = (t_ldsk **)mCalloc(ldsk->n_next-1+NEXT_BLOCK_SIZE,sizeof(t_ldsk *));

  pos = 0;
  for(i=0;i<ldsk->n_next;i++)
    {
      if(ldsk->next[i] != rm)
        {
          new_next[pos] = ldsk->next[i];
          pos++;
        }
    }

  ldsk->n_next--;
  Free(ldsk->next);
  ldsk->next = new_next;
  /* printf("\n. remove next for ldsk %s n_next set to %d",ldsk->coord->id,ldsk->n_next); */
  /* fflush(NULL); */
}  

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
/* Returns the vector of average pairwise distances between ldsk on each disk */
phydbl *PHYREX_Mean_Pairwise_Distance_Between_Lineage_Locations(t_tree *tree)
{
  phydbl *dist;
  int block,n_disks,i,j, k;
  t_dsk *disk;
  
  PHYREX_Update_Lindisk_List(tree);

  dist = NULL;
  block   = 100;
  disk    = tree->young_disk;
  n_disks = 0;
  do
    {
      if(!n_disks) dist = (phydbl *)mCalloc(block,sizeof(phydbl));
      else if(!(n_disks%block)) dist = (phydbl *)mRealloc(dist,n_disks+block,sizeof(phydbl));
      
      dist[n_disks] = 0.0;
      for(i=0;i<disk->n_ldsk_a-1;i++)
        {
          for(j=i+1;j<disk->n_ldsk_a;j++)
            {
              for(k=0;k<tree->mmod->n_dim;k++)
                {
                  dist[n_disks] += FABS(disk->ldsk_a[i]->coord->lonlat[k] - disk->ldsk_a[j]->coord->lonlat[k]);
                  printf("\n * %d %f %f %f",
                         k,
                         disk->ldsk_a[i]->coord->lonlat[k],
                         disk->ldsk_a[j]->coord->lonlat[k],
                         tree->mmod->lim_up->lonlat[k]);
                }
            }
        }
      dist[n_disks] /= (phydbl)(disk->n_ldsk_a * (disk->n_ldsk_a-1) / 2.);

      printf("\n %d %f %f",disk->n_ldsk_a,disk->time,dist[n_disks]);

      n_disks++;
      disk = disk->prev;
    }
  while(disk->prev);

  return(dist);

}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl PHYREX_Random_Select_Time_Between_Jumps(t_tree *tree)
{
  t_dsk  *disk,**valid_disks;
  int n_valid_disks,block,select;
  phydbl time;

  valid_disks   = NULL;
  disk          = NULL;
  block         = 100;

  assert(!(tree->young_disk->next));

  disk = tree->young_disk->prev;
  n_valid_disks = 0;
  do
    {
      if(disk->ldsk != NULL && disk->prev != NULL)
        {
          if(!n_valid_disks) valid_disks = (t_dsk **)mCalloc(block,sizeof(t_dsk *));
          else if(!(n_valid_disks%block)) valid_disks = (t_dsk **)mRealloc(valid_disks,n_valid_disks+block,sizeof(t_dsk *));
          valid_disks[n_valid_disks] = disk;
          n_valid_disks++;
        }
      disk = disk->prev;
    }
  while(disk->prev);

  if(!n_valid_disks) return -1.0;

  select = Rand_Int(0,n_valid_disks-1);

  time = valid_disks[select]->time - valid_disks[select]->ldsk->prev->disk->time;
  
  Free(valid_disks);

  return(time);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

int PHYREX_Is_In_Ldscape(t_ldsk *ldsk, t_phyrex_mod *mmod)
{
  int j;
  for(j=0;j<mmod->n_dim;j++)
    if(ldsk->coord->lonlat[j] > mmod->lim_up->lonlat[j] ||
       ldsk->coord->lonlat[j] < mmod->lim_do->lonlat[j]) return NO;
  return YES;
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl PHYREX_Mean_Time_Between_Events(t_tree *tree)
{
  int n_inter;
  phydbl T;

  n_inter = PHYREX_Total_Number_Of_Intervals(tree);  
  T = -tree->times->nd_t[tree->n_root->num];
  return((phydbl)(T/n_inter));
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void PHYREX_Rand_Pairs_Coal_Times_Dist(t_tree *tree)
{
  t_node *anc;
  phydbl dist;
  int i, j;

  i = Rand_Int(0,tree->n_otu-1);
  j = Rand_Int(0,tree->n_otu-1);

  if(i == j) PhyML_Printf("\nxxWxx 0.0 0.0");
  else
    {

      anc = Find_Lca_Pair_Of_Nodes(tree->a_nodes[i],tree->a_nodes[j],tree);
      if(anc == NULL) 
        {
          PhyML_Fprintf(stderr,"\n. %s",Write_Tree(tree));
          PhyML_Fprintf(stderr,"\n. %s %s",tree->a_nodes[i]->name,tree->a_nodes[j]->name);
          Generic_Exit(__FILE__,__LINE__,__FUNCTION__);            
        }
      
      PhyML_Printf("\nxxWxx %12f",tree->times->nd_t[anc->num]);
      dist = Euclidean_Dist(tree->a_nodes[i]->ldsk->coord,tree->a_nodes[j]->ldsk->coord);
      PhyML_Printf(" %f",dist);
    }
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void PHYREX_Read_Tip_Coordinates(t_tree *tree)
{
  char *s;
  FILE *fp;
  int i,*done,found_sw,found_ne;
  phydbl sw_lon, sw_lat,  ne_lon, ne_lat;
  t_node *n;
  
  assert(tree->young_disk);
  
  s    = (char *)mCalloc(T_MAX_LINE,sizeof(char));
  fp   = tree->io->fp_in_coord;
  done = (int *)mCalloc(tree->n_otu,sizeof(int));

  found_sw = NO;
  found_ne = NO;

  for(i=0;i<tree->n_otu;i++) done[i] = NO;

  do
    {
      if(fscanf(fp,"%s",s) == EOF) break;
      for(i=0;i<strlen(s);++i) if(s[i] == '#') break; /* skip comment */
      if(i != strlen(s)) continue;
      
      for(i=0;i<tree->n_otu;i++) if(strstr(tree->a_nodes[i]->name,s)) break;

      if(i != tree->n_otu) /* Found a match */
        {
          assert(tree->a_nodes[i]->ldsk);
          // First column: latitude, second one: longitude
          if(fscanf(fp,"%lf",&(tree->a_nodes[i]->ldsk->coord->lonlat[1])) == EOF) break;
          if(fscanf(fp,"%lf",&(tree->a_nodes[i]->ldsk->coord->lonlat[0])) == EOF) break;          
          done[i] = YES;
        }
      else
        {
          if(!strcmp(s,"|SouthWest|") || !strcmp(s,"|southwest|") || !strcmp(s,"|Southwest|"))
            {
              found_sw = YES;
              if(fscanf(fp,"%lf",&(sw_lat)) == EOF) break;              
              if(fscanf(fp,"%lf",&(sw_lon)) == EOF) break;
            }
          else if(!strcmp(s,"|NorthEast|") || !strcmp(s,"|northeast|") || !strcmp(s,"|Northeast|"))
            {
              found_ne = YES;
              if(fscanf(fp,"%lf",&(ne_lat)) == EOF) break;
              if(fscanf(fp,"%lf",&(ne_lon)) == EOF) break;
            }
          else /* Haven't found any match but still need to skip long and lat for unsampled location */
            {
              phydbl dum;
              if(fscanf(fp,"%lf",&dum) == EOF) break;
              if(fscanf(fp,"%lf",&dum) == EOF) break;
            }
        }
    }
  while(1);
  
  if(found_ne == NO)
    {
      PhyML_Fprintf(stderr,"\n. Could not find coordinates for northernmost  point.");
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
    }

  if(found_sw == NO)
    {
      PhyML_Fprintf(stderr,"\n. Could not find coordinates for southernmost point.");
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
    }

  for(i=0;i<tree->n_otu;i++) 
    if(done[i] == NO) 
      {
        PhyML_Fprintf(stderr,"\n. Could not find coordinates for '%s'.",tree->a_nodes[i]->name);
        Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
      }

  n = NULL;
  for(i=0;i<tree->n_otu;i++) 
    {
      n = tree->a_nodes[i];
      
      PhyML_Printf("\n. Coordinates of '%-50s': %12f\t %12f",
                   tree->a_nodes[i]->name,
                   n->ldsk->coord->lonlat[0],
                   n->ldsk->coord->lonlat[1]);
    }

  tree->mmod->lim_up->lonlat[0] = ne_lon;
  tree->mmod->lim_up->lonlat[1] = ne_lat;

  tree->mmod->lim_do->lonlat[0] = sw_lon;
  tree->mmod->lim_do->lonlat[1] = sw_lat;
  

  Free(s);
  Free(done);

}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl PHYREX_Tree_Height(t_tree *tree)
{
  t_dsk *disk;
  
  disk = tree->young_disk;  
  while(disk && disk->prev) disk = disk->prev;
  
  return(disk->time);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

int PHYREX_Random_Insert_Ldsk_In_Next_List(t_ldsk *ins, t_ldsk *where)
{
  int size, pos, i, *rk;
  t_ldsk **next_cpy;

  size = where->n_next;

  rk = (int *)mCalloc(size,sizeof(int));
  next_cpy = (t_ldsk **)mCalloc(size,sizeof(t_ldsk *));

  for(i=0;i<size;i++) next_cpy[i] = where->next[i];

  pos  = Rand_Int(0,size);

  for(i=0;i<size;i++)
    {
      if(i < pos) rk[i] = i;
      else        rk[i] = i+1;
    }

  PHYREX_Make_Lindisk_Next(where);
  
  for(i=0;i<size;i++) where->next[rk[i]] = next_cpy[i];
  
  where->next[pos]= ins;


  Free(rk);
  Free(next_cpy);

  return pos;
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void PHYREX_Insert_Ldsk_In_Next_List(t_ldsk *ins, int pos, t_ldsk *where)
{
  int size, i, *rk;
  t_ldsk **next_cpy;

  size = where->n_next;

  rk = (int *)mCalloc(size,sizeof(int));
  next_cpy = (t_ldsk **)mCalloc(size,sizeof(t_ldsk *));

  for(i=0;i<size;i++) next_cpy[i] = where->next[i];

  for(i=0;i<size;i++)
    {
      if(i < pos) rk[i] = i;
      else        rk[i] = i+1;
    }

  PHYREX_Make_Lindisk_Next(where);
  
  for(i=0;i<size;i++) where->next[rk[i]] = next_cpy[i];
  
  where->next[pos]= ins;

  Free(rk);
  Free(next_cpy);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
/* Remove path between young ldsk and old ldsk. If young->prev == old, */
/* then young->prev set to NULL and old->next[dir_to_young] set */
/* to NULL as well. */

t_ldsk *PHYREX_Remove_Path(t_ldsk *young, t_ldsk *old, int *pos_old, t_tree *tree)
{
  t_ldsk *path,*ldsk;
  int dir_old_young,jumps;

  path = NULL;
    
  dir_old_young = PHYREX_Get_Next_Direction(young,old);
  assert(dir_old_young >= 0);
  *pos_old = dir_old_young;
  PHYREX_Remove_Lindisk_Next(old,old->next[dir_old_young]);


  jumps = 1;
  ldsk = young->prev;
  while(ldsk != old)
    {
      if(jumps == 1) path = ldsk;
      PHYREX_Remove_Disk(ldsk->disk);
      ldsk = ldsk->prev;
      jumps++;
    }

  if(jumps == 1)
    path = NULL;
  else
    {
      /* Set end of path to NULL */
      ldsk = path;
      while(ldsk->prev != old) { ldsk = ldsk->prev; assert(ldsk); }
      ldsk->prev = NULL;
    }

  
  /* path == NULL if young->prev = old, otherwise path points to the first ldsk */
  /* after jump away from young towards the past (i.e., towards old) */
  return(path);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void PHYREX_Insert_Path(t_ldsk *young, t_ldsk *old, t_ldsk *path, int pos, t_tree *tree)
{
  t_ldsk *ldsk;

  assert(path != young);

  if(path == NULL)
    {
      young->prev = old;
      PHYREX_Insert_Ldsk_In_Next_List(young,pos,old);      
    }
  else
    {
      /* Attach path to the young ldsk */
      path->next[0] = young;
      young->prev = path;
      
      ldsk = path;
      do
        {
          PHYREX_Insert_Disk(ldsk->disk,tree);
          ldsk = ldsk->prev;
        }
      while(ldsk);

      /* Get to the end of path */
      ldsk = path;
      while(ldsk->prev != NULL) { ldsk = ldsk->prev; assert(ldsk); }
      
      /* Attach it to old ldsk (both ways)*/
      PHYREX_Insert_Ldsk_In_Next_List(ldsk,pos,old);            
      ldsk->prev = old;
    }
      
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl PHYREX_Time_Tree_Length(t_tree *tree)
{
  phydbl len;
  int i;
  t_dsk *disk;

  disk = tree->young_disk;
  while(disk->prev) disk = disk->prev;

  len = 0.0;
  for(i=0;i<disk->ldsk->n_next;i++) PHYREX_Time_Tree_Length_Pre(disk->ldsk,disk->ldsk->next[i],&len,tree);
  
  assert(!(isnan(len) || isinf(len)));

  return(len);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void PHYREX_Time_Tree_Length_Pre(t_ldsk *a, t_ldsk *d, phydbl *len, t_tree *tree)
{
  int i;

  (*len) += FABS(a->disk->time - d->disk->time);
  
  if(d->disk->next == NULL) return;
  else
    for(i=0;i<d->n_next;i++) PHYREX_Time_Tree_Length_Pre(d,d->next[i],len,tree);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

int PHYREX_Is_On_Path(t_ldsk *target, t_ldsk *young, t_ldsk *old)
{
  t_ldsk *ldsk;

  if(target == young || target == old) return NO;

  assert(!(young->disk->time < old->disk->time));

  ldsk = young->prev; 
  while(ldsk != old)
    {
      if(ldsk == target) return YES;
      ldsk = ldsk->prev;
      assert(!(ldsk == NULL));
    }
  return NO;
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

int PHYREX_Path_Len(t_ldsk *young, t_ldsk *old)
{
  t_ldsk *ldsk;
  int len;

  len = 0;
  ldsk = young;
  while(ldsk != old)
    {
      len++;
      ldsk = ldsk->prev;
    }
  len++;
  return len;
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void PHYREX_Print_Disk_Lk(t_tree *tree)
{
  t_dsk *disk;

  PHYREX_Update_Lindisk_List(tree);

  disk = tree->young_disk->prev;

  do
    {
      PhyML_Printf("\n. Disk: %s time: %12f lk: %12f cumlk: %12f",
                   disk->id,
                   disk->time,
                   PHYREX_Lk_Core(disk,tree),
                   -1.);
      
      disk = disk->prev;
    }
  while(disk->prev);
  
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

t_ldsk *PHYREX_Find_Lca_Pair_Of_Ldsk(t_ldsk *n1, t_ldsk *n2, t_tree *tree)
{
  t_ldsk **list1, **list2, *lca;
  int len1, len2;

  assert(n1);
  assert(n2);

  if(n1 == n2) return(n1);

  PHYREX_Get_List_Of_Ancestors(n1,&list1,&len1,tree);
  PHYREX_Get_List_Of_Ancestors(n2,&list2,&len2,tree);
  
  len1 = len1-1;
  len2 = len2-1;
         
  assert(list1[len1] == list2[len2]);

  do
    {
      if((len1 < 0 || len2 < 0) || (list1[len1] != list2[len2])) break;
      len1--;
      len2--;
    }
  while(len1 && len2);
  

  lca = list1[len1+1];

  Free(list1);
  Free(list2);
  
  assert(lca);

  return(lca);
  
} 

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void PHYREX_Get_List_Of_Ancestors(t_ldsk *start, t_ldsk ***list, int *len, t_tree *tree)
{
  int block,i;
  t_ldsk *ldsk;

  assert(start);
  block = 100;

  *list = (t_ldsk **)mCalloc(block,sizeof(t_ldsk *));
  
  ldsk = start;
  i = 0;
  do
    {
      (*list)[i] = ldsk;
      if(!(i%block)) *list = (t_ldsk **)mRealloc(*list,i+block,sizeof(t_ldsk *));
      i++;
      ldsk = ldsk->prev;
    }
  while(ldsk);

  (*len) = i;
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl PHYREX_Dist_To_Lca(t_ldsk *d, t_ldsk *lca)
{
  return(fabs(d->disk->time - lca->disk->time));      
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl PHYREX_Dist_Between_Two_Ldsk(t_ldsk *n1, t_ldsk *n2, t_tree *tree)
{
  t_ldsk *lca;

  lca = PHYREX_Find_Lca_Pair_Of_Ldsk(n1,n2,tree);
  
  return(PHYREX_Dist_To_Lca(n1,lca)+PHYREX_Dist_To_Lca(n2,lca));
}


/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void PHYREX_Print_MultiTypeTree_Config_File(int n_sites, char *filename, t_tree *tree)
{
  int i, j, n_demes;
  char *s,**deme_names;
  FILE *fp;


  fp = Openfile(filename,WRITE);
  assert(fp);

  deme_names = (char **)mCalloc(n_sites,sizeof(char *));

  n_demes = 0;
  for(i=0;i<tree->n_otu;i++)
    {
      s = strrchr(tree->a_nodes[i]->ldsk->coord->id,'_');
      for(j=0;j<n_demes;j++) if(!strcmp(s+1,deme_names[j])) break;
      if(j == n_demes)
        {
          deme_names[n_demes] = (char *)mCalloc(strlen(s+1)+1,sizeof(char));
          strcpy(deme_names[n_demes],s+1);
          n_demes++;
        }
    }

  // n_demes is the number of non-empty sampling sites 


  PhyML_Fprintf(fp,"<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>");
  PhyML_Fprintf(fp,"\n<beast beautitemplate='MultiTypeTree' beautistatus='' namespace=\"beast.core:beast.evolution.alignment:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.nuc:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood\" version=\"2.0\">");

  PhyML_Fprintf(fp,"\n<data id=\"data\" name=\"alignment\">");


  for(i=0;i<tree->n_otu;i++)
    {
      PhyML_Fprintf(fp,"\n<sequence id=\"%s\" taxon=\"%s\" totalcount=\"4\" value=\"%s\"/>",
                   tree->a_nodes[i]->ldsk->coord->id,
                   /* tree->a_nodes[i]->coord->id, */
                   tree->a_nodes[i]->name,
                   tree->a_nodes[i]->c_seq->state);
    }

  PhyML_Fprintf(fp,"\n</data>");

  PhyML_Fprintf(fp,"\n<map name=\"Uniform\" >beast.math.distributions.Uniform</map>");
  PhyML_Fprintf(fp,"\n<map name=\"Exponential\" >beast.math.distributions.Exponential</map>");
  PhyML_Fprintf(fp,"\n<map name=\"LogNormal\" >beast.math.distributions.LogNormalDistributionModel</map>");
  PhyML_Fprintf(fp,"\n<map name=\"Normal\" >beast.math.distributions.Normal</map>");
  PhyML_Fprintf(fp,"\n<map name=\"Beta\" >beast.math.distributions.Beta</map>");
  PhyML_Fprintf(fp,"\n<map name=\"Gamma\" >beast.math.distributions.Gamma</map>");
  PhyML_Fprintf(fp,"\n<map name=\"LaplaceDistribution\" >beast.math.distributions.LaplaceDistribution</map>");
  PhyML_Fprintf(fp,"\n<map name=\"prior\" >beast.math.distributions.Prior</map>");
  PhyML_Fprintf(fp,"\n<map name=\"InverseGamma\" >beast.math.distributions.InverseGamma</map>");
  PhyML_Fprintf(fp,"\n<map name=\"OneOnX\" >beast.math.distributions.OneOnX</map>");


  PhyML_Fprintf(fp,"\n<run id=\"mcmc\" spec=\"MCMC\" chainLength=\"1000000000\">");
  PhyML_Fprintf(fp,"\n<state id=\"state\" storeEvery=\"10000\">");
  PhyML_Fprintf(fp,"\n<stateNode id=\"Tree.t:data\" spec=\"beast.evolution.tree.StructuredCoalescentMultiTypeTree\">");
  PhyML_Fprintf(fp,"\n<migrationModel id=\"migModelInit.t:data\" spec=\"beast.evolution.tree.MigrationModel\">");


  s = (char *)mCalloc(T_MAX_LINE,sizeof(char));
  For(i,n_demes*(n_demes-1)) strcat(s,"1.0 ");
  PhyML_Fprintf(fp,"\n<parameter id=\"RealParameter.0\" dimension=\"%d\" estimate=\"false\" name=\"rateMatrix\">%s</parameter>",n_demes*(n_demes-1),s);
  Free(s);

  s = (char *)mCalloc(T_MAX_LINE,sizeof(char));
  for(i=0;i<n_demes;i++) strcat(s,"1.0 ");
  PhyML_Fprintf(fp,"\n<parameter id=\"RealParameter.01\" dimension=\"%d\" estimate=\"false\" name=\"popSizes\">%s</parameter>",n_demes,s);
  Free(s);

  PhyML_Fprintf(fp,"\n</migrationModel>");

  PhyML_Fprintf(fp,"\n<typeTrait id=\"typeTraitSet.t:data\" spec=\"beast.evolution.tree.TraitSet\" traitname=\"type\" value=\"");

  for(i=0;i<tree->n_otu;i++)
    {
      s = strrchr(tree->a_nodes[i]->ldsk->coord->id,'_');
      PhyML_Fprintf(fp,"%s=%s",
                   tree->a_nodes[i]->ldsk->coord->id,
                   s+1);

      if(i < tree->n_otu-1) PhyML_Fprintf(fp,",");
      else PhyML_Fprintf(fp,"\">");
    }

  PhyML_Fprintf(fp,"\n<taxa id=\"TaxonSet.0\" spec=\"TaxonSet\">");
  PhyML_Fprintf(fp,"\n<alignment idref=\"data\"/>");
  PhyML_Fprintf(fp,"\n</taxa>");
  PhyML_Fprintf(fp,"\n</typeTrait>");
  PhyML_Fprintf(fp,"\n<taxonset idref=\"TaxonSet.0\"/>");
  PhyML_Fprintf(fp,"\n</stateNode>");
  PhyML_Fprintf(fp,"\n<parameter id=\"kappa.s:data\" lower=\"0.0\" name=\"stateNode\">2.0</parameter>");

  s = (char *)mCalloc(T_MAX_LINE,sizeof(char));
  for(i=0;i<n_demes;i++) strcat(s,"1.0 ");
  PhyML_Fprintf(fp,"\n<parameter id=\"popSizes.t:data\" dimension=\"%d\" name=\"stateNode\">%s</parameter>",n_demes,s);
  Free(s);

  s = (char *)mCalloc(T_MAX_LINE,sizeof(char));
  For(i,n_demes*(n_demes-1)) strcat(s,"1.0 ");
  PhyML_Fprintf(fp,"\n<parameter id=\"rateMatrix.t:data\" dimension=\"%d\" name=\"stateNode\">%s</parameter>",n_demes*(n_demes-1),s);
  Free(s);


  PhyML_Fprintf(fp,"\n<parameter id=\"freqParameter.s:data\" dimension=\"4\" lower=\"0.0\" name=\"stateNode\" upper=\"1.0\">0.25</parameter>");
  PhyML_Fprintf(fp,"\n</state>");


  PhyML_Fprintf(fp,"\n<distribution id=\"posterior\" spec=\"util.CompoundDistribution\">");
  PhyML_Fprintf(fp,"\n<distribution id=\"prior\" spec=\"util.CompoundDistribution\">");
  PhyML_Fprintf(fp,"\n<prior id=\"KappaPrior.s:data\" name=\"distribution\" x=\"@kappa.s:data\">");
  PhyML_Fprintf(fp,"\n<LogNormal id=\"LogNormalDistributionModel.0\" name=\"distr\">");
  PhyML_Fprintf(fp,"\n<parameter id=\"RealParameter.02\" estimate=\"false\" name=\"M\">1.0</parameter>");
  PhyML_Fprintf(fp,"\n<parameter id=\"RealParameter.03\" estimate=\"false\" name=\"S\">1.25</parameter>");
  PhyML_Fprintf(fp,"\n</LogNormal>");
  PhyML_Fprintf(fp,"\n</prior>");

  PhyML_Fprintf(fp,"\n<prior id=\"popSizesPrior.t:data\" name=\"distribution\" x=\"@popSizes.t:data\">");
  PhyML_Fprintf(fp,"\n<LogNormal id=\"LogNormalDistributionModel.01\" name=\"distr\">");
  PhyML_Fprintf(fp,"\n<parameter id=\"RealParameter.04\" estimate=\"false\" name=\"M\">1.0</parameter>");
  PhyML_Fprintf(fp,"\n<parameter id=\"RealParameter.05\" estimate=\"false\" lower=\"0.0\" name=\"S\" upper=\"5.0\">1.25</parameter>");
  PhyML_Fprintf(fp,"\n</LogNormal>");
  PhyML_Fprintf(fp,"\n</prior>");

  PhyML_Fprintf(fp,"\n<prior id=\"rateMatrixPrior.t:data\" name=\"distribution\" x=\"@rateMatrix.t:data\">");
  PhyML_Fprintf(fp,"\n<LogNormal id=\"LogNormalDistributionModel.02\" name=\"distr\">");
  PhyML_Fprintf(fp,"\n<parameter id=\"RealParameter.06\" estimate=\"false\" name=\"M\">1.0</parameter>");
  PhyML_Fprintf(fp,"\n<parameter id=\"RealParameter.07\" estimate=\"false\" lower=\"0.0\" name=\"S\" upper=\"5.0\">1.25</parameter>");
  PhyML_Fprintf(fp,"\n</LogNormal>");
  PhyML_Fprintf(fp,"\n</prior>");

  PhyML_Fprintf(fp,"\n<distribution id=\"structuredCoalescent.t:data\" spec=\"multitypetree.distributions.StructuredCoalescentTreeDensity\" multiTypeTree=\"@Tree.t:data\">");
  PhyML_Fprintf(fp,"\n<migrationModel id=\"migModel.t:data\" spec=\"beast.evolution.tree.MigrationModel\" popSizes=\"@popSizes.t:data\" rateMatrix=\"@rateMatrix.t:data\">");
  PhyML_Fprintf(fp,"\n</migrationModel>");
  PhyML_Fprintf(fp,"\n</distribution>");

  PhyML_Fprintf(fp,"\n<distribution id=\"likelihood\" spec=\"util.CompoundDistribution\">");
  PhyML_Fprintf(fp,"\n<distribution id=\"treeLikelihood.data\" spec=\"TreeLikelihood\" data=\"@data\" tree=\"@Tree.t:data\">");
  PhyML_Fprintf(fp,"\n<siteModel id=\"SiteModel.s:data\" spec=\"SiteModel\">");
  PhyML_Fprintf(fp,"\n<parameter id=\"mutationRate.s:data\" estimate=\"false\" name=\"mutationRate\">1.0</parameter>");
  PhyML_Fprintf(fp,"\n<parameter id=\"gammaShape.s:data\" estimate=\"false\" name=\"shape\">1.0</parameter>");
  PhyML_Fprintf(fp,"\n<parameter id=\"proportionInvariant.s:data\" estimate=\"false\" lower=\"0.0\" name=\"proportionInvariant\" upper=\"1.0\">0.0</parameter>");
  PhyML_Fprintf(fp,"\n<substModel id=\"hky.s:data\" spec=\"HKY\" kappa=\"@kappa.s:data\">");
  PhyML_Fprintf(fp,"\n<frequencies id=\"estimatedFreqs.s:data\" spec=\"Frequencies\" frequencies=\"@freqParameter.s:data\"/>");
  PhyML_Fprintf(fp,"\n</substModel>");
  PhyML_Fprintf(fp,"\n</siteModel>");
  PhyML_Fprintf(fp,"\n<branchRateModel id=\"StrictClock.c:data\" spec=\"beast.evolution.branchratemodel.StrictClockModel\">");
  PhyML_Fprintf(fp,"\n<parameter id=\"clockRate.c:data\" estimate=\"false\" name=\"clock.rate\">%G</parameter>",1.0);
  PhyML_Fprintf(fp,"\n</branchRateModel>");
  PhyML_Fprintf(fp,"\n</distribution>");
  PhyML_Fprintf(fp,"\n</distribution>");
  PhyML_Fprintf(fp,"\n</distribution>");
  PhyML_Fprintf(fp,"\n</distribution>");
  PhyML_Fprintf(fp,"\n");
  PhyML_Fprintf(fp,"\n<operator id=\"STX.t:data\" spec=\"multitypetree.operators.TypedSubtreeExchange\" migrationModel=\"@migModel.t:data\" multiTypeTree=\"@Tree.t:data\" weight=\"10.0\"/>");
  PhyML_Fprintf(fp,"\n<operator id=\"TWB.t:data\" spec=\"multitypetree.operators.TypedWilsonBalding\" alpha=\"0.2\" migrationModel=\"@migModel.t:data\" multiTypeTree=\"@Tree.t:data\" weight=\"10.0\"/>");
  PhyML_Fprintf(fp,"\n<operator id=\"NR.t:data\" spec=\"multitypetree.operators.NodeRetype\" migrationModel=\"@migModel.t:data\" multiTypeTree=\"@Tree.t:data\" weight=\"10.0\"/>");
  PhyML_Fprintf(fp,"\n<operator id=\"NSR1.t:data\" spec=\"multitypetree.operators.NodeShiftRetype\" migrationModel=\"@migModel.t:data\" multiTypeTree=\"@Tree.t:data\" rootOnly=\"true\" weight=\"10.0\"/>");
  PhyML_Fprintf(fp,"\n<operator id=\"NSR2.t:data\" spec=\"multitypetree.operators.NodeShiftRetype\" migrationModel=\"@migModel.t:data\" multiTypeTree=\"@Tree.t:data\" noRoot=\"true\" weight=\"10.0\"/>");
  PhyML_Fprintf(fp,"<operator id=\"MTU.t:data\" spec=\"multitypetree.operators.MultiTypeUniform\" includeRoot=\"true\" migrationModel=\"@migModel.t:data\" multiTypeTree=\"@Tree.t:data\" weight=\"10.0\"/>\n");
  PhyML_Fprintf(fp,"<operator id=\"MTTS.t:data\" spec=\"multitypetree.operators.MultiTypeTreeScale\" migrationModel=\"@migModel.t:data\" multiTypeTree=\"@Tree.t:data\" scaleFactor=\"0.98\" useOldTreeScaler=\"true\" weight=\"10.0\"/>\n");
  PhyML_Fprintf(fp,"\n<operator id=\"MTTUpDown.t:data\" spec=\"multitypetree.operators.MultiTypeTreeScale\" migrationModel=\"@migModel.t:data\" multiTypeTree=\"@Tree.t:data\" scaleFactor=\"0.98\" useOldTreeScaler=\"true\" weight=\"10.0\">");
  PhyML_Fprintf(fp,"\n<parameter idref=\"popSizes.t:data\"/>");
  PhyML_Fprintf(fp,"\n</operator>");
  PhyML_Fprintf(fp,"\n<operator id=\"KappaScaler.s:data\" spec=\"ScaleOperator\" parameter=\"@kappa.s:data\" scaleFactor=\"0.5\" weight=\"0.1\"/>");
  PhyML_Fprintf(fp,"\n<operator id=\"popSizesScaler.t:data\" spec=\"ScaleOperator\" parameter=\"@popSizes.t:data\" scaleFactor=\"0.8\" weight=\"1.0\"/>");
  PhyML_Fprintf(fp,"\n<operator id=\"rateMatrixScaler.t:data\" spec=\"ScaleOperator\" parameter=\"@rateMatrix.t:data\" scaleFactor=\"0.8\" weight=\"1.0\"/>");
  PhyML_Fprintf(fp,"\n<operator id=\"FrequenciesExchanger.s:data\" spec=\"DeltaExchangeOperator\" delta=\"0.01\" weight=\"0.1\">");
  PhyML_Fprintf(fp,"\n<parameter idref=\"freqParameter.s:data\"/>");
  PhyML_Fprintf(fp,"\n</operator>");
  PhyML_Fprintf(fp,"\n");
  PhyML_Fprintf(fp,"\n<logger id=\"tracelog\" fileName=\"$(filebase).log\" logEvery=\"10000\">");
  PhyML_Fprintf(fp,"\n<log idref=\"likelihood\"/>");
  PhyML_Fprintf(fp,"\n<log idref=\"prior\"/>");
  PhyML_Fprintf(fp,"\n<log idref=\"treeLikelihood.data\"/>");
  PhyML_Fprintf(fp,"\n<log id=\"treeHeight.t:data\" spec=\"beast.evolution.tree.TreeHeightLogger\" tree=\"@Tree.t:data\"/>");
  /* PhyML_Fprintf(fp,"\n<log id=\"treeLength.t:data\" spec=\"multitypetree.util.TreeLengthLogger\" tree=\"@Tree.t:data\"/>"); */
  /* PhyML_Fprintf(fp,"\n<log id=\"changeCounts.t:data\" spec=\"multitypetree.util.TypeChangeCounts\" migrationModel=\"@migModel.t:data\" multiTypeTree=\"@Tree.t:data\"/>"); */
  /* PhyML_Fprintf(fp,"\n<log id=\"rootTypeLogger.t:data\" spec=\"multitypetree.util.TreeRootTypeLogger\" multiTypeTree=\"@Tree.t:data\"/>"); */
  PhyML_Fprintf(fp,"\n<log id=\"migModelLogger.t:data\" spec=\"multitypetree.util.MigrationModelLogger\" migrationModel=\"@migModel.t:data\" multiTypeTree=\"@Tree.t:data\"/>");
  PhyML_Fprintf(fp,"\n</logger>");
  PhyML_Fprintf(fp,"\n");
  PhyML_Fprintf(fp,"\n<logger id=\"screenlog\" logEvery=\"50000\">");
  PhyML_Fprintf(fp,"\n<log id=\"ESS.0\" spec=\"util.ESS\" arg=\"@posterior\"/>");
  PhyML_Fprintf(fp,"\n<log idref=\"likelihood\"/>");
  PhyML_Fprintf(fp,"\n</logger>");
  PhyML_Fprintf(fp,"\n");
  PhyML_Fprintf(fp,"\n</run>");
  PhyML_Fprintf(fp,"\n</beast>");
  
  for(i=0;i<n_demes;i++) Free(deme_names[i]);
  Free(deme_names);

  fclose(fp);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

int PHYREX_Number_Of_Sampled_Demes(t_tree *tree)
{
  int i,j,n_demes;
  t_dsk *disk;
  char **deme_list;

  deme_list = (char **)mCalloc(tree->n_otu,sizeof(char *));

  disk = tree->young_disk;

  n_demes = 0;
  for(i=0;i<tree->n_otu;i++)
    {
      for(j=0;j<n_demes;j++)
        {
          if(deme_list[j] != NULL && !strcmp(strstr(disk->ldsk_a[i]->coord->id,"_deme"),deme_list[j]))
            {
              break;
            }
        }

      if(j == n_demes)
        {
          deme_list[j] = strstr(disk->ldsk_a[i]->coord->id,"_deme");
          /* printf("\n. deme_list[%d]: %s",j,deme_list[j]); */
          n_demes++;
        }
    }

  Free(deme_list);


  return(n_demes);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
/* Returns the number of ldsks going from disk towards the past */
int PHYREX_Number_Of_Outgoing_Ldsks(t_dsk *disk)
{
  int i;

  if(disk->ldsk == NULL) return disk->n_ldsk_a;
  else
    {
      int n_out;

      n_out = 0;
      for(i=0;i<disk->n_ldsk_a;++i)
        {
          if(disk->ldsk_a[i] && disk->ldsk_a[i]->prev != disk->ldsk)
            {
              n_out++;
            }
        }
      return(n_out+1); /* +1 so as to count disk->ldsk in */
    }
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
/* Select uniformly at random a lineage going "out of" disk (towards the past) */
t_ldsk *PHYREX_Random_Select_Outgoing_Ldsk(t_dsk *disk)
{
  int i,*permut,n_ldsk_a_out;
  t_ldsk **ldsk_a_out,*target_ldsk;

  ldsk_a_out = (t_ldsk **)mCalloc(disk->n_ldsk_a,sizeof(t_ldsk *));

  n_ldsk_a_out = 0;
  if(disk->ldsk != NULL)
    {
      ldsk_a_out[0] = disk->ldsk;
      n_ldsk_a_out = 1;
    }
  
  for(i=0;i<disk->n_ldsk_a;++i)
    {
      if(disk->ldsk_a[i] && disk->ldsk_a[i]->prev != disk->ldsk)
        {
          ldsk_a_out[n_ldsk_a_out] = disk->ldsk_a[i];
          n_ldsk_a_out++;
        }
    }
  
  permut = Permutate(n_ldsk_a_out);
  target_ldsk = ldsk_a_out[permut[0]];

  Free(permut);
  Free(ldsk_a_out);
  
  return(target_ldsk);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
/* Remove and free all disks and ldsks except for sampled disks.
   Connect the sampled disks together */
 
void PHYREX_Strip_And_Reconnect_Tree(t_tree *tree)
{
  t_dsk *disk,*prev;
  int i,orig_size;
  
  disk = tree->young_disk;
  do
    {
      if(disk->age_fixed == YES)
        {
          orig_size = disk->n_ldsk_a;
          for(i=0;i<orig_size;++i)
            {
              if(disk->ldsk_a[i] &&
                 disk->ldsk_a[i]->disk != disk)
                {
                  disk->ldsk_a[i] = NULL;
                  disk->n_ldsk_a--;
                }
            }
          for(i=0;i<disk->n_ldsk_a;++i) disk->ldsk_a[i]->prev = NULL;
        }
      disk = disk->prev;
    }
  while(disk);

  disk = tree->young_disk;
  do
    {
      prev = disk->prev;
      if(disk->age_fixed == NO)
        {
          PHYREX_Remove_Disk(disk);
          Free_Ldisk(disk->ldsk);
          Free_Disk(disk);
        }
      disk = prev;
    }
  while(disk);
  
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

int PHYREX_Scale_All(phydbl scale, t_dsk *start_disk, t_tree *tree)
{
  t_dsk *disk;
  t_dsk **sorted_disk;
  int n_disk,n_disks_scaled,n_nodes_scaled,sorted,i;

  n_disk         = 0;
  n_disks_scaled = 0;
  n_nodes_scaled = 0;

  disk = start_disk->prev;
  assert(disk);
  
  do
    {      
      if(disk->age_fixed == NO)
        {
          disk->time = disk->time * scale + start_disk->time * (1.-scale);
          n_disks_scaled++;
          if(disk->ldsk && disk->ldsk->n_next > 1) n_nodes_scaled++;
        }
      
      n_disk++;
      disk = disk->prev;
    }
  while(disk);

  sorted_disk = (t_dsk **)mCalloc(n_disk,sizeof(t_dsk *));
  disk = start_disk->prev;
  n_disk = 0;
  do
    {
      sorted_disk[n_disk] = disk;
      n_disk++;
      disk = disk->prev;      
    }
  while(disk);
  

  do
    {
      sorted = YES;
      for(i=0;i<n_disk-1;++i)
        {
          if(sorted_disk[i]->time < sorted_disk[i+1]->time)
            {
              sorted = NO;
              disk             = sorted_disk[i];
              sorted_disk[i]   = sorted_disk[i+1];
              sorted_disk[i+1] = disk;
            }
        }
    }
  while(sorted == NO);

  
  for(i=0;i<n_disk;++i)
    {
      PHYREX_Remove_Disk(sorted_disk[i]);
      PHYREX_Insert_Disk(sorted_disk[i],tree);
    }

  Free(sorted_disk);

  if(tree->mmod->id == SLFV_GAUSSIAN || tree->mmod->id == SLFV_UNIFORM) return(n_disks_scaled);
  else if(tree->mmod->id == RW || tree->mmod->id == RRW) return(n_nodes_scaled);
  else return(-1);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void PHYREX_Oldest_Sampled_Disk(t_tree *tree)
{
  t_dsk *disk;

  disk = tree->young_disk;
  do
    {
      if(disk->age_fixed == YES) tree->old_samp_disk = disk;
      disk = disk->prev;      
    }
  while(disk);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl PHYREX_Time_Of_Descendants(t_ldsk *ldsk, t_tree *tree)
{
  if(ldsk->disk->age_fixed == YES || ldsk->disk == tree->young_disk) return ldsk->disk->time;
  else
    {
      int i;
      phydbl t,min;
      
      min = tree->young_disk->time;
      
      for(i=0;i<ldsk->n_next;++i)
        {
          t = PHYREX_Time_Of_Descendants(ldsk->next[i],tree);
          if(t < min) min = t;
        }
      return(min);
    }
  return(1.);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl PHYREX_Time_Of_Prev_Sampled_Disk(t_dsk *disk, t_tree *tree)
{

  if(disk->prev == NULL) return -INFINITY;
  else
    {      
      disk = disk->prev;
      while(disk && disk->age_fixed == NO) disk = disk->prev;
      if(!disk) return -INFINITY;
      else return(disk->time);
    }
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl PHYREX_Time_Of_Next_Sampled_Disk(t_dsk *disk, t_tree *tree)
{

  if(disk->next == NULL) return tree->young_disk->time;
  else
    {      
      disk = disk->next;
      while(disk && disk->age_fixed == NO) disk = disk->next;
      if(!disk) return tree->young_disk->time;
      else return(disk->time);
    }
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl PHYREX_Lk_Time_Component(t_tree *tree)
{
  phydbl lnL;
  t_dsk *disk;
  phydbl dt;
  int n_evt;

  n_evt = 0;
  dt    = 0.0;
  lnL   = 0.0;
  disk  = tree->young_disk->prev;
  do
    {

      if(disk->age_fixed == NO)
        {
          dt += fabs(disk->next->time - disk->time);
          n_evt++;
        }
      
      disk = disk->prev;

    }
  while(disk);

  lnL += n_evt * log(tree->mmod->lbda) - tree->mmod->lbda * dt;

  return(lnL);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

t_ldsk *PHYREX_Find_Ldsk_From_Id(char *id, t_ldsk *root)
{
  t_ldsk *ldsk;

  ldsk = NULL;
  
  if(!strcmp(root->coord->id,id)) return(root);

  if(root->n_next == 0) return(NULL);
  
  for(int i=0;i<root->n_next;++i)
    {
      ldsk = PHYREX_Find_Ldsk_From_Id(id,root->next[i]);
      if(ldsk != NULL) break;
    }
  
  return(ldsk);
}


/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

t_geo_coord *PHYREX_Mean_Next_Loc(t_ldsk *ldsk, t_tree *tree)
{
  t_geo_coord *mean;
  int i,j;
  
  mean = GEO_Make_Geo_Coord(tree->mmod->n_dim);

  assert(ldsk->n_next > 0);
  
  for(i=0;i<ldsk->n_next;++i)
    {
      for(j=0;j<tree->mmod->n_dim;++j)
        {
          mean->lonlat[j] += ldsk->next[i]->coord->lonlat[j];
        }            
    }
  
  for(j=0;j<tree->mmod->n_dim;++j)
    {
      mean->lonlat[j] /= (phydbl)ldsk->n_next;
    }            
  
  return(mean);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
// realized sigsq = (1/d) E(\sum_x D_x^2), where $D_x$ is the
// square euclidean distance along the x-axis between the
// location of a lineage at time $t$ and time $t+1$.
phydbl PHYREX_Root_To_Tip_Realized_Sigsq(t_tree *tree)
{
  t_dsk *disk,*root_disk;
  int i;
  phydbl res,*sigsq;
  phydbl t_root, t_tip;

  sigsq = (phydbl *)mCalloc(tree->young_disk->n_ldsk_a,sizeof(phydbl));
  
  disk = tree->young_disk;
  while(disk->prev) disk = disk->prev;
  root_disk = disk;
  
  t_root = root_disk->time;
  t_tip = tree->young_disk->time;
  assert(t_tip - t_root > 0.0);
  
  disk = tree->young_disk;
  for(i=0;i<disk->n_ldsk_a;++i)
    sigsq[i] = pow(Euclidean_Dist(root_disk->ldsk->coord,disk->ldsk_a[i]->coord) / (t_tip - t_root),2);

  res = Quantile(sigsq,disk->n_ldsk_a,0.5);
  
  Free(sigsq);
  
  return(res/tree->mmod->n_dim);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
// Mean Haversine distance (in km) per year (measured on path between root and tips on the most recent disk) 
phydbl PHYREX_Realized_Dispersal_Dist(t_tree *tree)
{
  t_dsk *disk,*root_disk;
  phydbl dist,t_tip,t_root;
  int i;
  
  disk = tree->young_disk;
  while(disk->prev) disk = disk->prev;
  root_disk = disk;

  t_root = root_disk->time;
  t_tip = tree->young_disk->time;
  assert(t_tip - t_root > 0.0);
  \
  disk = tree->young_disk;
  dist = 0.0;
  for(i=0;i<disk->n_ldsk_a;++i)
    {
      dist += Haversine_Distance(root_disk->ldsk->coord->lonlat[0],
                                 root_disk->ldsk->coord->lonlat[1],
                                 disk->ldsk_a[i]->coord->lonlat[0],
                                 disk->ldsk_a[i]->coord->lonlat[1])/(t_tip - t_root);
        
    }
  
  return(dist/disk->n_ldsk_a);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl PHYREX_Tip_To_Root_Realized_Sigsq(t_tree *tree)
{
  t_dsk *disk;
  phydbl *sigsq,res;
  int i;


  sigsq = (phydbl *)mCalloc(tree->young_disk->n_ldsk_a,sizeof(phydbl));
  
  disk = tree->young_disk;
  for(i=0;i<disk->n_ldsk_a;++i)
    {
      sigsq[i] =
        pow(Euclidean_Dist(disk->ldsk_a[i]->coord,
                           disk->ldsk_a[i]->prev->coord)/fabs(disk->time - disk->ldsk_a[i]->prev->disk->time),2);
                   
    }

  res = Quantile(sigsq,disk->n_ldsk_a,0.5);

  Free(sigsq);
  
  return(res/tree->mmod->n_dim);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl PHYREX_Tip_To_Root_Realized_Bis_Sigsq(t_tree *tree)
{
  t_dsk *disk;
  phydbl sumdist,sumt;
  int i;

  sumdist = 0.0;
  sumt = 0.0;
  disk = tree->young_disk;
  for(i=0;i<disk->n_ldsk_a;++i)
    {
      sumdist += pow(Euclidean_Dist(disk->ldsk_a[i]->coord,
                                    disk->ldsk_a[i]->prev->coord),2);
      sumt += tree->mmod->n_dim * fabs(disk->time - disk->ldsk_a[i]->prev->disk->time);                   
    }
  
  return(sumdist/sumt);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl PHYREX_Tip_To_Root_Realized_Ter_Sigsq(t_tree *tree)
{
  t_ldsk *ldsk;
  phydbl mean,dist;
  int i;

  mean = 0.0;
  for(i=0;i<tree->young_disk->n_ldsk_a;++i)
    {
      ldsk = tree->young_disk->ldsk_a[i];
      while(ldsk && ldsk->prev && fabs(ldsk->prev->disk->time-tree->young_disk->ldsk_a[i]->disk->time) < 1.0)
        {
          ldsk = ldsk->prev;
          if(ldsk == NULL) return(-1.);
        }
      dist = Euclidean_Dist(ldsk->coord,tree->young_disk->ldsk_a[i]->coord);
      mean += pow(dist,2);
      
      /* PhyML_Printf("\n. %3d (%12f %12f) (%12f %12f) dist: %12f mean: %12f", */
      /*              i, */
      /*              tree->young_disk->ldsk_a[i]->coord->lonlat[0], */
      /*              tree->young_disk->ldsk_a[i]->coord->lonlat[1], */
      /*              ldsk->coord->lonlat[0],                   */
      /*              ldsk->coord->lonlat[1], */
      /*              Euclidean_Dist(ldsk->coord,tree->young_disk->ldsk_a[i]->coord), */
      /*              mean); */
                   
    }
  
  return(mean/(tree->young_disk->n_ldsk_a*tree->mmod->n_dim));
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void PHYREX_Label_Nodes_With_Locations(t_tree *tree)
{
  t_dsk *disk;
  t_node *n;
  phydbl lon,lat;


  lon = lat = -1.;
  
  for(int i=0;i<tree->n_otu;++i)
    {
      if(tree->a_nodes[i]->label == NULL)
        {
          tree->a_nodes[i]->label = Make_Label();
          tree->a_nodes[i]->label->next = Make_Label();
        }

      lat = tree->a_nodes[i]->ldsk->coord->lonlat[1];
      lon = tree->a_nodes[i]->ldsk->coord->lonlat[0];

      sprintf(tree->a_nodes[i]->label->key,"&location");
      sprintf(tree->a_nodes[i]->label->val,"{%f,%f}",lat,lon);
      sprintf(tree->a_nodes[i]->label->next->key,"location");
      sprintf(tree->a_nodes[i]->label->next->val,"{%f,%f}",lat,lon);
    }

  disk = tree->young_disk->prev;

  do
    {
      if(disk->ldsk && disk->ldsk->nd != NULL)
        {
          n = disk->ldsk->nd;
          if(n->label == NULL)
            {
              n->label = Make_Label();
              n->label->next = Make_Label();
            }

          lat = disk->ldsk->coord->lonlat[1];
          lon = disk->ldsk->coord->lonlat[0];

          sprintf(n->label->key,"&location");
          sprintf(n->label->val,"{%f,%f}",lat,lon);
          sprintf(n->label->next->key,"location");
          sprintf(n->label->next->val,"{%f,%f}",lat,lon);

          /* Print same label on all internal nodes with exactly the */
          /* same coalescence time. */
          for(int i=tree->n_otu;i<2*tree->n_otu-1;++i)
            {
              if(tree->a_nodes[i] != n &&
                 Are_Equal(tree->times->nd_t[i],
                           tree->times->nd_t[n->num],
                           1.E-10) == YES)
                {
                  n = tree->a_nodes[i];
                  if(n->label == NULL)
                    {
                      n->label = Make_Label();
                      n->label->next = Make_Label();
                    }

                  lat = disk->ldsk->coord->lonlat[1];
                  lon = disk->ldsk->coord->lonlat[0];
                  
                  sprintf(n->label->key,"&location");
                  sprintf(n->label->val,"{%f,%f}",lat,lon);
                  sprintf(n->label->next->key,"location");
                  sprintf(n->label->next->val,"{%f,%f}",lat,lon);
                }
            }
        }
      
      disk = disk->prev;
    }
  while(disk);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void PHYREX_Label_Edges(t_tree *tree)
{
  for(int i=0;i<2*tree->n_otu-1;++i)
    {
      if(tree->a_edges[i]->label == NULL)
        {
          tree->a_edges[i]->label = Make_Label();
          tree->a_edges[i]->label->next = Make_Label();
        }
      
      sprintf(tree->a_edges[i]->label->key,"&rate");
      sprintf(tree->a_edges[i]->label->val,"0.0");
      
      sprintf(tree->a_edges[i]->label->next->key,"location.rate");
      sprintf(tree->a_edges[i]->label->next->val,"0.0");
    }


}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void PHYREX_Remove_All_Disks_Except_Coal_And_Tips(t_tree *tree)
{
  t_dsk *disk,*dum;

  disk = tree->young_disk;
  dum = NULL;
  
  do
    {
      dum = disk->prev;

      if(disk->age_fixed == NO)
        {
          if(disk->ldsk == NULL) 
            {
              PHYREX_Remove_Disk(disk);
              Free_Disk(disk);
            }
          else if(disk->ldsk != NULL && disk->ldsk->nd == NULL)
            {
              t_ldsk *old_ldsk,*young_ldsk;
              int dir_old_young;
              
              assert(disk->ldsk->n_next == 1);

              old_ldsk   = disk->ldsk->prev;
              young_ldsk = disk->ldsk->next[0];

              dir_old_young = PHYREX_Get_Next_Direction(disk->ldsk,old_ldsk);
              assert(dir_old_young != -1);

              /* New connections between old_ldsk and young_ldsk */
              old_ldsk->next[dir_old_young] = young_ldsk;
              young_ldsk->prev              = old_ldsk;

              PHYREX_Remove_Disk(disk);

              Free_Ldisk(disk->ldsk);
              Free_Disk(disk);

            }
        }

      disk = dum;
    }
  while(disk);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
  
phydbl PHYREX_Path_Logdensity(t_ldsk *young, t_ldsk *old, phydbl sd, t_tree *tree)
{
  switch(tree->mmod->id)
    {
    case SLFV_UNIFORM : case SLFV_GAUSSIAN : 
      {
        return(SLFV_Path_Logdensity(young,old,sd,tree));
        break;
      }      
    case RW : case RRW :
      {
        return(0.0);
        break;
      }
    default : assert(FALSE);
    }
  return(-1.);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void PHYREX_Sample_Path(t_ldsk *young, t_ldsk *old, phydbl sd, phydbl *global_hr, t_tree *tree)
{
  switch(tree->mmod->id)
    {
    case SLFV_UNIFORM : case SLFV_GAUSSIAN : 
      {
        SLFV_Sample_Path(young,old,sd,global_hr,tree);
        break;
      }
      
    case RW : case RRW :
      {
        assert(FALSE);
        break;
      }
    default : assert(FALSE);
    }
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

t_ldsk *PHYREX_Generate_Path(t_ldsk *young, t_ldsk *old, phydbl n_evt, phydbl sd, t_tree *tree)
{
  switch(tree->mmod->id)
    {
    case SLFV_UNIFORM : case SLFV_GAUSSIAN : 
      {
        return(SLFV_Generate_Path(young,old,n_evt,sd,tree));
        break;
      }      
    case RW : case RRW : 
      {
        return(SLFV_Generate_Path(young,old,0,sd,tree));
        break;
      }
    default : assert(FALSE);
    }
  return(NULL);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

t_tree *PHYREX_Simulate(int n_otu, int n_sites, phydbl w, phydbl h, phydbl  lbda, phydbl rad, phydbl mu, int r_seed, int modid)
{
  switch(modid)
    {
    case SLFV_UNIFORM : case SLFV_GAUSSIAN : 
      {
        return(SLFV_Simulate(n_otu,n_sites,w,h,lbda,rad,mu,r_seed));
        break;
      }
      
    case RW : case RRW : 
      {
        assert(FALSE);
        break;
      }
    default : assert(FALSE);
    }
  return(NULL);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void PHYREX_Simulate_Backward_Core(t_dsk *disk,int avoid_multiple_mergers, t_tree *tree)
{
  switch(tree->mmod->id)
    {
    case SLFV_UNIFORM : case SLFV_GAUSSIAN : 
      {
        SLFV_Simulate_Backward_Core(disk,avoid_multiple_mergers,tree);
        break;
      }
      
    case RW : case RRW : 
      {
        SLFV_Simulate_Backward_Core(disk,avoid_multiple_mergers,tree);
        break;
      }

    default : assert(FALSE);

    }
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl PHYREX_Update_Sigsq(t_tree *tree)
{  
  switch(tree->mmod->id)
    {
    case SLFV_UNIFORM : case SLFV_GAUSSIAN : 
      {
        return(SLFV_Update_Sigsq(tree));
        break; 
      }

    case RW : case RRW :
      {
        return(tree->mmod->sigsq);
        break;
      }

    default : assert(FALSE);

    }

  return(-1.);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
  
