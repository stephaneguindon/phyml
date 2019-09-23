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
  /* PHYREX_Main_Simulate(argc,argv); */
  option *io;
  io = Get_Input(argc,argv);
  Free(io);
  return(0);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

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


  tree = mixt_tree;
  do
    {
      // All rate stuctures point to the same object
      tree->rates = mixt_tree->rates;
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



    // Looking for XML node with rate-across-lineage info
  xnd = XML_Search_Node_Name("clockrate",YES,xroot);
  
  if(xnd != NULL)
    {
      char *value = XML_Get_Attribute_Value(xnd,"value");
      
      if(value == NULL)
        {
          PhyML_Fprintf(stderr,"\n. Please specify a value for the average rate of substitution (clockrate),");
          PhyML_Fprintf(stderr,"\n. e.g., <clockrate value=\"1E-5\"/>.");
          PhyML_Fprintf(stderr,"\n. See the manual for more options.");
          assert(FALSE);
        }
      else
        {
          mixt_tree->rates->clock_r       = atof(value);
          mixt_tree->rates->clock_r_fixed = YES;
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

  mixt_tree->mmod = PHYREX_Make_Migrep_Model(2);
  mixt_tree->mmod->n_dim = 2;

  
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
      TIMES_Randomize_Tree_With_Time_Constraints(mixt_tree->rates->a_cal[0],mixt_tree);
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
  /* mixt_tree->mmod->lbda = Uni()*(mixt_tree->mmod->max_lbda - mixt_tree->mmod->min_lbda) + mixt_tree->mmod->min_lbda; */
  mixt_tree->mmod->lbda = 1.0;
  mixt_tree->mmod->mu   = Uni()*(mixt_tree->mmod->max_mu - mixt_tree->mmod->min_mu) + mixt_tree->mmod->min_mu;
  mixt_tree->mmod->rad  = Uni()*(mixt_tree->mmod->max_rad - mixt_tree->mmod->min_rad) + mixt_tree->mmod->min_rad;
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

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int PHYREX_Main_Simulate(int argc, char *argv[])
{
  t_tree *tree;
  int seed,pid,i;
  char *s;
  t_dsk *disk;
  int n_sites,n_otus;
  phydbl width,height;
  phydbl lbda, rad,mu;
  
  s = (char *)mCalloc(T_MAX_NAME,sizeof(char));

  pid     = getpid();
  seed    = pid;
  n_otus  = (int)atoi(argv[1]);
  /* n_sites = (int)atoi(argv[2]); */
  n_sites = 1;
  width   = (phydbl)atof(argv[2]); 
  height  = (phydbl)atof(argv[3]); 
  lbda    = (phydbl)atof(argv[4]);
  rad     = (phydbl)atof(argv[5]);
  mu      = (phydbl)atof(argv[6]);

  
  printf("\n. seed: %d",seed);
  srand(seed);
  
  tree = PHYREX_Simulate(n_otus,n_sites,width,height,lbda,rad,mu,seed);
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

t_tree *PHYREX_Simulate_Independent_Loci(int n_otu, int n_loci, phydbl w, phydbl h, int r_seed)
{
  t_tree *tree;
  int n_dim,i;
  t_phyrex_mod *mmod;
  option *io;
  t_mod *mod;
  t_opt *s_opt;
  calign *cdata;
  phydbl area;
  phydbl subst_rate;
  int locus_idx;

  n_dim = 2; // 2-dimensional landscape
  area  = w * h;

  io    = (option *)Make_Input();
  mod   = (t_mod *)Make_Model_Basic();
  s_opt = (t_opt *)Make_Optimiz();

  Set_Defaults_Input(io);
  Set_Defaults_Model(mod);
  Set_Defaults_Optimiz(s_opt);


  // io stuff
  io->mod      = mod;
  mod->io      = io;
  mod->s_opt   = s_opt;
  io->r_seed   = r_seed;
  io->n_otu    = n_otu;
  io->init_len = n_loci; /* sequence length */
  io->data     = Make_Empty_Alignment(io);
  io->colalias = NO;

  Make_Model_Complete(io->mod);
  Set_Model_Name(io->mod);

  cdata = Compact_Data(io->data,io);
  Free_Seq(io->data,io->n_otu);
  
  Print_Settings(io);


  // tree stuff
  tree = Make_Tree_From_Scratch(n_otu,cdata);
  Random_Tree(tree);
  Connect_CSeqs_To_Nodes(cdata,io,tree);
  tree->rates = RATES_Make_Rate_Struct(tree->n_otu);
  io->rates->model = STRICTCLOCK;
  RATES_Init_Rate_Struct(tree->rates,io->rates,tree->n_otu);
    
  tree->data      = cdata;
  tree->mod       = mod;
  tree->io        = io;
  tree->n_pattern = tree->data->crunch_len;

  Init_Model(cdata,mod,io);
      
  tree->mod->ras->n_catg = 1;
  tree->mod->whichmodel  = HKY85;
  tree->mod->kappa->v    = 4.0;
    
  Make_Tree_For_Pars(tree);
  Make_Tree_For_Lk(tree);
  Make_Spr(tree);

  // migrep model stuff */
  mmod = PHYREX_Make_Migrep_Model(n_dim);
  tree->mmod = mmod;
  PHYREX_Init_Migrep_Mod(mmod,n_dim,0.0,0.0,w,h);


  /* /\* Death proba param *\/ */
  /* mmod->mu = Uni(); */
  /* /\* Theta (radius) *\/ */
  /* mmod->rad = Uni()*(1.0 - 0.2) + 0.2; */
  /* /\* Rate of events *\/ */
  /* mmod->lbda = 1.0; */
  /* /\* Population density *\/ */
  /* mmod->rho = Uni()*(3.0 - 0.1) + 0.1; */


  /* Death proba param */
  mmod->mu = 0.8;
  /* Theta (radius) */
  mmod->rad = 0.3;
  /* Rate of events */
  mmod->lbda = 10.*w*h;
  /* (Actual, not effective) population density */
  mmod->rho = 1.;

  
  
  // Duration of a generation in number of events
  mmod->gen_cal_time = 1./(2.*mmod->mu*pow(mmod->rad,2)/pow(w*h,2)*
                           (sqrt(2.)*mmod->rad*exp(-.5*pow(h/mmod->rad,2)) + h*sqrt(PI)*erf(sqrt(2.)*h/(2.*mmod->rad)) - sqrt(2.)*mmod->rad)*
                           (sqrt(2.)*mmod->rad*exp(-.5*pow(w/mmod->rad,2)) + w*sqrt(PI)*erf(sqrt(2.)*w/(2.*mmod->rad)) - sqrt(2.)*mmod->rad));
  // Divide by rate of events per calendar time unit to get duration of gen. in calendar time unit
  mmod->gen_cal_time *= 1./mmod->lbda;

  /* Dispersal parameter (per generation) */
  mmod->sigsq = 4.*pow(mmod->rad,4)*mmod->lbda*PI*mmod->mu/area * mmod->gen_cal_time;

 
  phydbl p = Prob_Two_Lineages_Coal_One_Event(w,h,mmod->mu,mmod->rad);
  printf("\n. p.sim: %G p.appx: %G",p,4.*pow(mmod->mu,2)*pow(PI,2)*pow(mmod->rad,4)/(pow(w*h,2)));
  phydbl rhoe = 1./(1.-exp(-mmod->lbda*p*mmod->gen_cal_time));
  rhoe /= area;
  
  PhyML_Printf("\n. lambda=%G; mu=%G; rad=%G; clockr=%G; sigsq=%G; neigh=%G; rhoe(small rad)=%G rhoe(numerical)=%G rhoe(large rad)=%G one_gen=%G",
               mmod->lbda,
               mmod->mu,
               mmod->rad,
               tree->rates->clock_r,
               mmod->sigsq,
               2./mmod->mu,
               (1./(1.-exp(-4.*pow(mmod->mu,2)*pow(PI,2)*pow(mmod->rad,4)/(pow(w*h,2))*mmod->lbda*mmod->gen_cal_time)))/area,
               rhoe,
               /* (1./(1.-exp(-mmod->lbda*mmod->mu/(w*h))))/area, */
               (1./(1.-exp(-mmod->mu)))/area,
               mmod->gen_cal_time);
  fflush(NULL);
  /* Exit("\n"); */

  
  // Initialize position of sampled individuals
  tree->young_disk = PHYREX_Make_Disk_Event(mmod->n_dim,tree->n_otu);
  PHYREX_Init_Disk_Event(tree->young_disk,mmod->n_dim,NULL);
  tree->young_disk->time             = 0.0;
  tree->young_disk->mmod             = mmod;
  tree->young_disk->centr->lonlat[0] = .5*(mmod->lim_up->lonlat[0]-mmod->lim_do->lonlat[0]);
  tree->young_disk->centr->lonlat[1] = .5*(mmod->lim_up->lonlat[1]-mmod->lim_do->lonlat[1]);      
  tree->young_disk->n_ldsk_a         = tree->n_otu;  
  
  for(i=0;i<tree->n_otu;++i) 
    {
      char *s;
      tree->young_disk->ldsk_a[i] = PHYREX_Make_Lindisk_Node(mmod->n_dim);
      PHYREX_Init_Lindisk_Node(tree->young_disk->ldsk_a[i],tree->young_disk,mmod->n_dim);
      s = (char *)mCalloc(strlen(tree->young_disk->ldsk_a[i]->coord->id)+1+20,sizeof(char));
      strcpy(s,tree->young_disk->ldsk_a[i]->coord->id);
      strcat(s,"_deme0\0");
      Free(tree->young_disk->ldsk_a[i]->coord->id);
      tree->young_disk->ldsk_a[i]->coord->id = s;
    }
      
  for(i=0;i<tree->n_otu;i++)
    {
      tree->young_disk->ldsk_a[i]->coord->lonlat[0] = Uni()*(tree->mmod->lim_up->lonlat[0]-tree->mmod->lim_do->lonlat[0]) + tree->mmod->lim_do->lonlat[0]; // longitude
      tree->young_disk->ldsk_a[i]->coord->lonlat[1] = Uni()*(tree->mmod->lim_up->lonlat[1]-tree->mmod->lim_do->lonlat[1]) + tree->mmod->lim_do->lonlat[1]; // latitude
    }



  // Generate sequences
  subst_rate = .0;
  locus_idx = 0;
  do
    {
      PHYREX_Simulate_Backward_Core(tree->young_disk,NO,tree);
  
      // Random selection of a pair. Print out physical distances
      // between tips and time of coalescence for this pair, at this
      // particular locus
      {
        int *permut = Permutate(tree->n_otu);
        t_ldsk *lin1 = tree->young_disk->ldsk_a[permut[0]];
        t_ldsk *lin2 = tree->young_disk->ldsk_a[permut[1]];
        phydbl dist =
          sqrt(pow(lin1->coord->lonlat[0]-lin2->coord->lonlat[0],2) +
               pow(lin1->coord->lonlat[1]-lin2->coord->lonlat[1],2));

        t_dsk *disk = tree->young_disk;
        do
          {
            if(disk->ldsk == lin1->prev) lin1 = lin1->prev;
            if(disk->ldsk == lin2->prev) lin2 = lin2->prev;
            if(lin1 == lin2) break; // found MRCA
            disk = disk->prev;
          }
        while(disk);

        assert(lin1 && lin2);
        
        PhyML_Printf("\n. #$# locus %4d ; physical distance: %G time to coalescence: %G",locus_idx,dist,lin1->disk->time);
      }
      
      PHYREX_Ldsk_To_Tree(tree);  

      
      Update_Ancestors(tree->n_root,tree->n_root->v[2],tree);
      Update_Ancestors(tree->n_root,tree->n_root->v[1],tree);
      RATES_Fill_Lca_Table(tree);


      tree->rates->bl_from_rt = YES;
      /* tree->rates->clock_r    = 1.E-3/FABS(T); // slow */
      /* tree->rates->clock_r    = 1.E-0/FABS(T); // fast */
      /* tree->rates->clock_r    = 1.E-1/FABS(T); // medium */
      /* tree->rates->clock_r = 1.E-7; // slow; */
      tree->rates->clock_r = 1.E-4; // medium;
      /* tree->rates->clock_r = 1.E-2; // fast; */

      
      PhyML_Printf("\n. #!# mutation rate at that locus: %G subst. per time unit",tree->rates->clock_r);
      subst_rate += tree->rates->clock_r;
      
      RATES_Update_Cur_Bl(tree);

      char *s = Write_Tree(tree,NO);
      PhyML_Printf("\n. #@# tree: %s",s);
      Free(s);

      Evolve(tree->data,tree->mod,locus_idx,tree);
      
      t_dsk *disk = tree->young_disk->prev;
      while(disk->prev)
        {
          disk = disk->prev;
          if(disk->next->ldsk != NULL) Free_Ldisk(disk->next->ldsk);
          Free_Disk(disk->next);
        }

      locus_idx++;
    }
  while(locus_idx < n_loci);

  PhyML_Printf("\n. #s# average rate of mutation : %G",subst_rate/(phydbl)n_loci); 
  
  tree->data->format = IBDSIM;
  PhyML_Printf("\n\n");
  Print_CSeq(stdout,NO,tree->data,tree);


  Exit("\n");

  return tree;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// Simulate Etheridge-Barton model backwards in time, following n_otu lineages
// on a rectangle of dimension width w x h
// See Kelleher, Barton & Etheridge, Bioinformatics, 2013.
t_tree *PHYREX_Simulate(int n_otu, int n_sites, phydbl w, phydbl h, phydbl  lbda, phydbl rad, phydbl mu, int r_seed)
{  
  t_tree *tree;
  int n_dim,i;
  t_phyrex_mod *mmod;
  t_dsk *disk;
  option *io;
  t_mod *mod;
  t_opt *s_opt;
  calign *cdata;
  /* phydbl T; */

  n_dim = 2; // 2-dimensional landscape

  io    = (option *)Make_Input();
  mod   = (t_mod *)Make_Model_Basic();
  s_opt = (t_opt *)Make_Optimiz();

  Set_Defaults_Input(io);
  Set_Defaults_Model(mod);
  Set_Defaults_Optimiz(s_opt);

  io->mod      = mod;
  mod->io      = io;
  mod->s_opt   = s_opt;
  io->r_seed   = r_seed;

  io->n_otu    = n_otu;
  io->init_len = 1000; /* sequence length */

  io->data = Make_Empty_Alignment(io);

  Make_Model_Complete(io->mod);
  Set_Model_Name(io->mod);

  io->colalias = NO;
  cdata = Compact_Data(io->data,io);
  Free_Seq(io->data,io->n_otu);

  tree = Make_Tree_From_Scratch(n_otu,cdata);

  Connect_CSeqs_To_Nodes(cdata,io,tree);
  
  tree->rates = RATES_Make_Rate_Struct(tree->n_otu);
  RATES_Init_Rate_Struct(tree->rates,io->rates,tree->n_otu);
  
  tree->data      = cdata;
  tree->mod       = mod;
  tree->io        = io;
  tree->n_pattern = tree->data->crunch_len;


  /* Allocate migrep model */
  mmod = PHYREX_Make_Migrep_Model(n_dim);
  tree->mmod = mmod;
  PHYREX_Init_Migrep_Mod(mmod,n_dim,0.0,0.0,w,h);

  mmod->lbda = lbda;
  mmod->rad = rad;
  mmod->mu = mu;
  mmod->sigsq = PHYREX_Update_Sigsq(tree);


  
  // Duration of a generation in number of events
  mmod->gen_cal_time = 1./(2.*mmod->mu*pow(mmod->rad,2)/pow(w*h,2)*
                           (sqrt(2.)*mmod->rad*exp(-.5*pow(h/mmod->rad,2)) + h*sqrt(PI)*erf(sqrt(2.)*h/(2.*mmod->rad)) - sqrt(2.)*mmod->rad)*
                           (sqrt(2.)*mmod->rad*exp(-.5*pow(w/mmod->rad,2)) + w*sqrt(PI)*erf(sqrt(2.)*w/(2.*mmod->rad)) - sqrt(2.)*mmod->rad));
  // Divide by rate of events per calendar time unit to get duration of gen. in calendar time unit
  mmod->gen_cal_time *= 1./mmod->lbda;


  /* Prepare start disk */
  tree->young_disk = PHYREX_Make_Disk_Event(tree->mmod->n_dim,tree->n_otu);
  tree->young_disk->time = 0.0;
  tree->young_disk->age_fixed = YES;
  tree->young_disk->centr->lonlat[0] = .5*(tree->mmod->lim_up->lonlat[0]-tree->mmod->lim_do->lonlat[0]);
  tree->young_disk->centr->lonlat[1] = .5*(tree->mmod->lim_up->lonlat[1]-tree->mmod->lim_do->lonlat[1]);   
  tree->young_disk->n_ldsk_a         = tree->n_otu;

  for(i=0;i<tree->n_otu;++i) 
    {
      char *s;
      tree->young_disk->ldsk_a[i] = PHYREX_Make_Lindisk_Node(tree->mmod->n_dim);
      PHYREX_Init_Lindisk_Node(tree->young_disk->ldsk_a[i],tree->young_disk,tree->mmod->n_dim);
      s = (char *)mCalloc(strlen(tree->young_disk->ldsk_a[i]->coord->id)+1+20,sizeof(char));
      strcpy(s,tree->young_disk->ldsk_a[i]->coord->id);
      strcat(s,"_deme0\0");
      Free(tree->young_disk->ldsk_a[i]->coord->id);
      tree->young_disk->ldsk_a[i]->coord->id = s;
    }
  
  /* Generate coordinates for the tip nodes (uniform distribution on the rectangle) */
  for(i=0;i<tree->n_otu;i++)
    {
      tree->young_disk->ldsk_a[i]->coord->lonlat[0] = Uni()*(tree->mmod->lim_up->lonlat[0]-tree->mmod->lim_do->lonlat[0])+tree->mmod->lim_do->lonlat[0]; // longitude
      tree->young_disk->ldsk_a[i]->coord->lonlat[1] = Uni()*(tree->mmod->lim_up->lonlat[1]-tree->mmod->lim_do->lonlat[1])+tree->mmod->lim_do->lonlat[1]; // latitude
    }

  for(i=0;i<tree->n_otu;i++)
    {
      tree->young_disk->ldsk_a[i]->nd = tree->a_nodes[i];
      tree->a_nodes[i]->ldsk = tree->young_disk->ldsk_a[i];
    }
  
  PHYREX_Simulate_Backward_Core(tree->young_disk,NO,tree);

  PHYREX_Ldsk_To_Tree(tree);
  Update_Ancestors(tree->n_root,tree->n_root->v[2],tree);
  Update_Ancestors(tree->n_root,tree->n_root->v[1],tree);
  RATES_Fill_Lca_Table(tree);
  
  
  /* min_rate = 1.E-5; */
  /* max_rate = 1.E-4; */

  /* T = PHYREX_Tree_Height(tree); */

  tree->rates->bl_from_rt = YES;
  tree->rates->clock_r    = 1.0;
  tree->rates->model      = STRICTCLOCK;

  RATES_Update_Cur_Bl(tree);

  Init_Model(cdata,mod,io);

  tree->mod->ras->n_catg   = 1;
  tree->mod->whichmodel    = HKY85;
  tree->mod->kappa->v      = 4.0;
    
  Make_Tree_For_Pars(tree);
  Make_Tree_For_Lk(tree);
  Make_Spr(tree);
  Evolve(tree->data,tree->mod,0,tree);

  /* PhyML_Printf("@@@ %G %G %G : ",mmod->lbda,mmod->mu,mmod->rad); */
  /* for(int i=0;i<tree->n_otu-1;++i) for(int j=i+1;j<tree->n_otu;++j) PhyML_Printf("%G ",PHYREX_Dist_Between_Two_Ldsk(tree->a_nodes[i]->ldsk,tree->a_nodes[j]->ldsk,tree)); */
  /* PhyML_Printf(" : "); */
  /* for(int i=0;i<tree->n_otu-1;++i) for(int j=i+1;j<tree->n_otu;++j) PhyML_Printf("%G ",Euclidean_Dist(tree->a_nodes[i]->ldsk->coord,tree->a_nodes[j]->ldsk->coord)); */
  /* PhyML_Printf("\n"); */
  /* for(int i=0;i<tree->n_otu-1;++i) PhyML_Printf("\n%s %G %G", */
  /*                                               tree->a_nodes[i]->name, */
  /*                                               tree->a_nodes[i]->ldsk->coord->lonlat[0], */
  /*                                               tree->a_nodes[i]->ldsk->coord->lonlat[1]); */
  /* PhyML_Printf("\n\n"); */
  /* Print_CSeq(stdout,NO,tree->data,tree); */
  /* Exit("\n"); */
  
  /* Init_Partial_Lk_Tips_Double(tree); */
  /* Init_Partial_Lk_Loc(tree); */

  disk = tree->young_disk->prev;
  while(disk->prev) disk = disk->prev;
  
  PhyML_Printf("\n. Useful parameters: lambda=%G; mu=%G; rad=%G; clockr=%G; sigsq=%G",
               mmod->lbda,
               mmod->mu,
               mmod->rad,
               tree->rates->clock_r,
               mmod->sigsq);

  PhyML_Printf("\n. Useful statistics: t.root=%f n.int=%d n.coal=%d n.hit=%d root.x=%f root.y=%f nt.div=%f\n",
               disk->time,
               PHYREX_Total_Number_Of_Intervals(tree),
               PHYREX_Total_Number_Of_Coal_Disks(tree),
               PHYREX_Total_Number_Of_Hit_Disks(tree),
               disk->ldsk->coord->lonlat[0],
               disk->ldsk->coord->lonlat[1],
               Nucleotide_Diversity(tree->data));
  
  
  PhyML_Printf("\n. Tree: ");
  PhyML_Printf("\n. %s \n",Write_Tree(tree,NO));

  PhyML_Printf("\n. Spatial coordinates: ");
  for(i=0;i<tree->n_otu;i++)
    {
      PhyML_Printf("\n. %15s: [%12f ; %12f]",
                   tree->young_disk->ldsk_a[i]->nd->name,
                   tree->young_disk->ldsk_a[i]->coord->lonlat[0],                   
                   tree->young_disk->ldsk_a[i]->coord->lonlat[1]);
    }
  PhyML_Printf("\n");

  return(tree);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
// Simulate Etheridge-Barton model backwards in time, following n_otu lineages
// on a rectangle.
phydbl PHYREX_Simulate_Backward_Core(t_dsk *init_disk, int avoid_multiple_mergers, t_tree *tree)
{
  t_dsk *disk,*new_disk,*oldest_disk;
  int i,j,reached_oldest_disk,err;
  phydbl prob_hit,u,new_time,lnL;
  t_phyrex_mod *mmod;
  t_node *n;
  
  mmod = tree->mmod;
  
  Get_Node_Ranks_From_Tip_Times(tree);
  
  // Get to the youngest node
  n = tree->a_nodes[0];
  while(n->rk_next) n = n->rk_next;

  disk = init_disk;  
  // Make sure at least one of the youngest nodes is on init_disk
  for(i=0;i<init_disk->n_ldsk_a;++i) if(init_disk->ldsk_a[i]->nd == n) break;
  assert(i != init_disk->n_ldsk_a);
  
  // Get to the oldest sampled disk
  disk = init_disk;
  while(disk->prev) disk = disk->prev;
  oldest_disk = disk;
  reached_oldest_disk = NO;
  if(init_disk->prev == NULL) reached_oldest_disk = YES; // Only one sampled disk
  

  lnL = 0.0;
  disk = init_disk;
  do
    {      
      /* MRCA reached */
      if(PHYREX_Number_Of_Outgoing_Ldsks(disk) == 1 && reached_oldest_disk == YES) break;

      /* Proposed new time */
      new_time = disk->time - Rexp(mmod->lbda);

      lnL += log(mmod->lbda) - mmod->lbda * fabs(disk->time - new_time);
      
      /* New time is older than previous sampled disk (disk->prev) */
      if(disk->prev && new_time < disk->prev->time)
        {
          new_disk = disk->prev;
          new_time = disk->prev->time;

          // Reached oldest disk -> set it to NULL
          // so as to indicate that it is ok to stop
          // the simulation, continue otherwise
          if(new_disk == oldest_disk) reached_oldest_disk = YES;
          
          /* Connect disk and new_disk */
          disk->prev = new_disk;      
          new_disk->next = disk;
          
          /* Time of next event */
          new_disk->time = new_time;
          
          /* Populate new_disk->ldsk_a array */
          PHYREX_Update_Lindisk_List_Core(new_disk,tree);
        }
      else
        {
          /* Create new disk */
          new_disk = PHYREX_Make_Disk_Event(mmod->n_dim,tree->n_otu);
          PHYREX_Init_Disk_Event(new_disk,mmod->n_dim,NULL);

          /* new_disk->prev now points to previous sampled disk */
          new_disk->prev = disk->prev;
          if(disk->prev) disk->prev->next = new_disk;
          
          /* Connect disk and new_disk */
          disk->prev = new_disk;      
          new_disk->next = disk;
          
          /* Time of next event */
          new_disk->time = new_time;
          
          /* Coordinates of new event */
          for(j=0;j<mmod->n_dim;++j) new_disk->centr->lonlat[j] = Uni()*(mmod->lim_up->lonlat[j]-mmod->lim_do->lonlat[j])+mmod->lim_do->lonlat[j];

          for(j=0;j<mmod->n_dim;++j) lnL += log(1./(mmod->lim_up->lonlat[j]-mmod->lim_do->lonlat[j]));

          
          /* Populate new_disk->ldsk_a array */
          PHYREX_Update_Lindisk_List_Core(new_disk,tree);
          
          new_disk->ldsk = NULL;
          /* Which lineages in new_disk->ldsk_a are hit? */
          for(i=0;i<new_disk->n_ldsk_a;++i)
            {
              prob_hit = -1.;
              switch(mmod->name)
                {
                case PHYREX_UNIFORM: 
                  { 
                    prob_hit = mmod->mu; 
                    break; 
                  }
                case PHYREX_NORMAL:  
                  { 
                    prob_hit = log(mmod->mu);
                    for(j=0;j<mmod->n_dim;++j)
                      {                    
                        prob_hit +=
                          -pow(new_disk->ldsk_a[i]->coord->lonlat[j] -
                               new_disk->centr->lonlat[j],2)/(2.*pow(mmod->rad,2));
                      }
                    prob_hit = exp(prob_hit);                    
                    break; 
                  }
                }
              
              u = Uni();
              if(!(u > prob_hit)) // disk->ldsk_a[i] is  hit
                {
                  lnL += log(prob_hit); 
                  // new_disk->ldsk_a[i] is hit -> coalesce (or just jump) to parent (i.e., new_disk->ldsk)
                  if(new_disk->ldsk == NULL)
                    {
                      new_disk->ldsk = PHYREX_Make_Lindisk_Node(mmod->n_dim);
                      PHYREX_Init_Lindisk_Node(new_disk->ldsk,new_disk,mmod->n_dim);
                      
                      // Generate coordinates for newly created ldsk
                      switch(tree->mmod->name)
                        {
                        case PHYREX_UNIFORM:
                          {
                            PHYREX_Runif_Rectangle_Overlap(new_disk->ldsk,new_disk,tree->mmod);
                            break;
                          }
                        case PHYREX_NORMAL:
                          {
                            PHYREX_Rnorm_Trunc(new_disk->ldsk,new_disk,tree->mmod);
                            for(j=0;j<tree->mmod->n_dim;++j)
                              lnL += Log_Dnorm_Trunc(new_disk->ldsk->coord->lonlat[j],
                                                     new_disk->centr->lonlat[j],
                                                     tree->mmod->rad,
                                                     tree->mmod->lim_do->lonlat[j],
                                                     tree->mmod->lim_up->lonlat[j],&err);
                            break;
                          }
                        }
                    }
                  new_disk->ldsk_a[i]->prev = new_disk->ldsk;
                  PHYREX_Make_Lindisk_Next(new_disk->ldsk);
                  new_disk->ldsk->next[new_disk->ldsk->n_next-1] = new_disk->ldsk_a[i];
                }
              else
                {
                  lnL += log(1. - prob_hit);
                }
            }
        }

      disk = new_disk;
    }
  while(1);  
  disk->prev = NULL;

  return(lnL);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
// Simulate Etheridge-Barton model forwards in time, following n_otu lineages
// on a rectangle of dimension width x h
t_sarea *PHYREX_Simulate_Forward_Core(int n_sites, t_tree *tree)
{
  t_dsk *disk,*new_disk;
  t_ldsk *ldsk,**ldsk_a_pop,**ldsk_a_samp,**ldsk_a_tmp,**ldsk_a_tips,*new_ldsk;
  t_ll *ldsk_list,*dum_ll;
  int i,j,n_disk,n_dim,n_otu,init_pop_size,curr_pop_size,parent_id,n_lineages,sample_size,n_poly,*permut,n_sampled_demes;
  phydbl sum,*parent_prob,prob_death,tree_height,max_x,max_y,trans_x,trans_y,t_sim, one_gen;
  short int dies,n_remain;
  t_phyrex_mod *mmod;
  t_poly **poly;
  t_sarea *area;
  short int *is_sampled;
  phydbl w,h;
  phydbl cx,cy;
  phydbl lx,ly;
  int n_survive,n_offspring, n_gen;
  
  mmod          = tree->mmod;
  n_dim         = 2;
  n_otu         = tree->n_otu;
  w             = (mmod->lim_up->lonlat[0]-mmod->lim_do->lonlat[0]);
  h             = (mmod->lim_up->lonlat[1]-mmod->lim_do->lonlat[1]);
  init_pop_size = mmod->rho*w*h;
  n_gen         = 10; // number of generations to simulate
  one_gen       = mmod->gen_cal_time; // time duration of one generation (calendar unit)
  t_sim         = n_gen * one_gen; // total simulation duration (calendar unit)

 
  /* Allocate and initialise first disk event */
  disk = PHYREX_Make_Disk_Event(n_dim,init_pop_size);
  PHYREX_Init_Disk_Event(disk,n_dim,NULL);
  disk->time = 0.0;
  disk->prev = NULL;
  disk->n_ldsk_a = init_pop_size;
  
  /* Allocate coordinates for all individuals in the starting population */
  ldsk_list = NULL;
  i = 0;
  do
    {
      ldsk = PHYREX_Make_Lindisk_Node(n_dim);
      PHYREX_Init_Lindisk_Node(ldsk,disk,n_dim);
      Push_Bottom_Linked_List(ldsk,&ldsk_list,NO);
      i++;
    }
  while(i != init_pop_size);

  
  /* Fill in array of lineages on the first disk */  
  dum_ll = ldsk_list->head;
  i = 0;
  do
    {
      disk->ldsk_a[i] = (t_ldsk *)dum_ll->v;
      dum_ll = dum_ll->next;
      i++;
    }
  while(dum_ll != NULL);

  
  
  /* Generate coordinates for all individuals in starting population */
  for(i=0;i<init_pop_size;i++)
    {
      disk->ldsk_a[i]->coord->lonlat[0] = Uni()*w; // longitude
      disk->ldsk_a[i]->coord->lonlat[1] = Uni()*h; // latitude
    }

  
  n_disk = 0;
  do
    {
      /* Create new disk */
      new_disk = PHYREX_Make_Disk_Event(n_dim,1);
      Free(new_disk->ldsk_a);
      new_disk->ldsk_a = NULL;
      PHYREX_Init_Disk_Event(new_disk,n_dim,NULL);
      new_disk->prev = disk;      
      disk->next = new_disk;
      n_disk++;

            
      /* Coordinates of event */
      new_disk->centr->lonlat[0] = Uni()*(mmod->lim_up->lonlat[0]-mmod->lim_do->lonlat[0])+mmod->lim_do->lonlat[0];
      new_disk->centr->lonlat[1] = Uni()*(mmod->lim_up->lonlat[1]-mmod->lim_do->lonlat[1])+mmod->lim_do->lonlat[1];

      cx = new_disk->centr->lonlat[0];
      cy = new_disk->centr->lonlat[1];
      

      /* Time of new disk */
      new_disk->time = disk->time + Rexp(mmod->lbda);

      /* PhyML_Printf("\n. new_disk->time: %f t_sim: %f",new_disk->time,t_sim); */
      
      /* Size of current population */
      curr_pop_size = disk->n_ldsk_a;

      if(curr_pop_size == 0)
        {
          PhyML_Fprintf(stderr,"\n== Population went extinct after %d events",n_disk);
          Exit("\n");
        }
      
      
      if(new_disk->time < t_sim)
        {
          /* Select one parent */
          parent_prob = (phydbl *)mCalloc(curr_pop_size,sizeof(phydbl));
          for(i=0;i<curr_pop_size;++i)
            {
              lx = disk->ldsk_a[i]->coord->lonlat[0];
              ly = disk->ldsk_a[i]->coord->lonlat[1];
              
              switch(mmod->name)
                {
                case PHYREX_UNIFORM:
                  {
                    if(PHYREX_Is_In_Disk(disk->ldsk_a[i]->coord,new_disk,mmod) == YES)
                      parent_prob[i] = 1.0;
                    else
                      parent_prob[i] = 0.0;
                    break;
                  }
                case PHYREX_NORMAL:
                  {
                    parent_prob[i] = 0.0;
                    parent_prob[i] += -pow(lx - cx,2)/(2.*pow(mmod->rad,2));
                    parent_prob[i] += -pow(ly - cy,2)/(2.*pow(mmod->rad,2));
                    parent_prob[i] = exp(parent_prob[i]);
                    break;
                  }
                }
            }
          
          sum = 0.0;
          for(i=0;i<curr_pop_size;++i) sum += parent_prob[i];
          
          if(sum < 1.E-100)
            {
              sum = curr_pop_size;
              for(i=0;i<curr_pop_size;++i) parent_prob[i] = 1.;
            }
          
          for(i=0;i<curr_pop_size;++i) parent_prob[i] /= sum;
          
          
          
          parent_id = Sample_i_With_Proba_pi(parent_prob,curr_pop_size);
          new_disk->ldsk = disk->ldsk_a[parent_id];
                    
          Free(parent_prob);
          
          
          /* Which lineages survive that event? */
          n_survive = 0;
          for(i=0;i<curr_pop_size;++i)
            {
              lx = disk->ldsk_a[i]->coord->lonlat[0];
              ly = disk->ldsk_a[i]->coord->lonlat[1];
              
              prob_death = 0.0;
              switch(mmod->name)
                {
                case PHYREX_UNIFORM:
                  {
                    if(PHYREX_Is_In_Disk(disk->ldsk_a[i]->coord,new_disk,mmod) == YES)
                      prob_death = mmod->mu;
                    break;
                  }
                case PHYREX_NORMAL:
                  {
                    prob_death = log(mmod->mu);
                    prob_death += -pow(lx - cx,2)/(2.*pow(mmod->rad,2));
                    prob_death += -pow(ly - cy,2)/(2.*pow(mmod->rad,2));
                    prob_death = exp(prob_death);
                    break;
                  }
                }
              
              dies = NO;
              if(Uni() < prob_death) dies = YES;
              
              if(dies == NO)
                {
                  if(n_survive == 0) new_disk->ldsk_a = (t_ldsk **)mCalloc(1,sizeof(t_ldsk *));
                  else new_disk->ldsk_a = (t_ldsk **)mRealloc(new_disk->ldsk_a,n_survive+1,sizeof(t_ldsk *));              
                  new_disk->ldsk_a[n_survive] = disk->ldsk_a[i];
                  n_survive++;
                }
            }
          
          /* PhyML_Printf("\n. n_survive: %d",n_survive); */
          
          /* Offspring */
          /* How many? */
          phydbl r = mmod->rad;
          phydbl u = mmod->mu;
          phydbl rho = mmod->rho;
          phydbl mean =
            0.5*PI*rho*u*r*r*
            (erf(0.5*sqrt(2.)/r*(cy-0.0))*erf(0.5*sqrt(2.)/r*(cx-0.0)) +
             erf(0.5*sqrt(2.)/r*(cy-h))*erf(0.5*sqrt(2.)/r*(cx-w)) -
             erf(0.5*sqrt(2.)/r*(cy-0.0))*erf(0.5*sqrt(2.)/r*(cx-w)) -
             erf(0.5*sqrt(2.)/r*(cy-h))*erf(0.5*sqrt(2.)/r*(cx-0.0)));
          
          n_offspring = Rpois(mean);
          
          new_disk->n_ldsk_a = n_survive + n_offspring;
          
          /* PhyML_Printf("\n. n_offspring: %d",n_offspring); */
          
          /* Where */
          for(i=0;i<n_offspring;++i)
            {
              /* New lindisk */
              new_ldsk = PHYREX_Make_Lindisk_Node(n_dim);
              PHYREX_Init_Lindisk_Node(new_ldsk,new_disk,n_dim);
              Push_Bottom_Linked_List(new_ldsk,&ldsk_list,NO);
              
              if(new_disk->ldsk_a == NULL) new_disk->ldsk_a = (t_ldsk **)mCalloc(n_survive+1,sizeof(t_ldsk *));
              else new_disk->ldsk_a = (t_ldsk **)mRealloc(new_disk->ldsk_a,n_survive+i+1,sizeof(t_ldsk *));              
              
              new_disk->ldsk_a[n_survive+i] = new_ldsk;
              
              
              /* Generate new location */
              switch(mmod->name)
                {
                case PHYREX_UNIFORM: { PHYREX_Runif_Rectangle_Overlap(new_ldsk,new_disk,mmod); break; }
                case PHYREX_NORMAL:  { PHYREX_Rnorm_Trunc(new_ldsk,new_disk,mmod); break; }
                }
              
              /* Connect to parent */
              new_ldsk->prev = disk->ldsk_a[parent_id];          
            }
        }
      else
        {
          new_disk->time = t_sim;
          new_disk->n_ldsk_a = disk->n_ldsk_a;
          new_disk->ldsk_a = (t_ldsk **)mCalloc(curr_pop_size,sizeof(t_ldsk *));
          for(i=0;i<curr_pop_size;++i) new_disk->ldsk_a[i] = disk->ldsk_a[i];
        }
      

      disk = new_disk;
          
      /* printf("\n. pop size: %6d # of events: %6d",disk->n_ldsk_a,n_disk); */
    }
  while(disk->time < t_sim);

  
  /* Dispersal stuff */
  /* { */
  /*   phydbl T = disk->time; // total simulation time (in calendar unit) */
  /*   t_ldsk *dum_ldsk = disk->ldsk_a[Rand_Int(0,disk->n_ldsk_a-1)]; // random selection of a lineage among all the lineages available at present time */
  /*   t_dsk *dum_dsk = disk; */
  /*   phydbl s = w*h; // area */
  /*   printf("\n. s: %f",s); fflush(NULL); */
  /*   phydbl gentime = 1./(2.*PI*mmod->rad*mmod->rad*mmod->mu*mmod->lbda/s); */
  /*   phydbl ssq = 0.0; // sum of squared difference between current and previous position of lineage, when previous pos != current pos. */
  /*   phydbl curr_pos,prev_pos; */
  /*   int curr_gen,prev_gen,nhits=0; */
  /*   prev_pos = dum_ldsk->coord->lonlat[0]; */
  /*   curr_pos = prev_pos; */
  /*   curr_gen = prev_gen = 1; */
  /*   // Compute sum of squared difference of positions */
  /*   do */
  /*     { */
  /*       curr_gen  = 1 + (int)(dum_dsk->time-T)/gentime; // increment the number of generations (generation time measured in calendar time units) */
  /*       curr_pos  = dum_ldsk->coord->lonlat[0]; // current position of lineage */
        
  /*       if(dum_ldsk->disk == dum_dsk) // lineage was born at that time */
  /*         { */
  /*           dum_ldsk = dum_ldsk->prev; // jump to parent */
  /*           nhits++; */
  /*         } */
        
  /*       if(curr_gen != prev_gen) */
  /*         { */
  /*           ssq += pow(curr_pos-prev_pos,2); */
  /*           prev_pos = curr_pos; */
  /*           prev_gen = curr_gen; */
  /*         } */
        
  /*       dum_dsk = dum_dsk->prev; */
  /*     } */
  /*   while(dum_dsk); */
    
  /*   PhyML_Printf("\n # var T nhits T/nhits gentime sigsq theta u lambda"); */
  /*   PhyML_Printf("\n a@z %G %f %d %f %f %f %f %f %f", */
  /*                (1./(T/gentime))*ssq, */
  /*                T, */
  /*                nhits, */
  /*                T/nhits, */
  /*                gentime, */
  /*                4*pow(mmod->rad,4)*mmod->lbda/s*PI*mmod->mu*gentime, */
  /*                mmod->rad, */
  /*                mmod->mu, */
  /*                mmod->lbda); */
    
  /*   Exit("\n"); */
  /* } */
  

  /* /\* Coalescence stuff *\/ */
  /* { */
  /*   t_ldsk *lin1, *lin2; */
  /*   int n_evts,have_coal,coal_evt; */
  /*   phydbl T,t; */

  /*   /\* Selection of two lineages at random *\/ */
  /*   while(disk->next) disk = disk->next; */
  /*   permut = Permutate(disk->n_ldsk_a); */
  /*   lin1 = disk->ldsk_a[permut[0]]; */
  /*   lin2 = disk->ldsk_a[permut[1]]; */
  /*   Free(permut); */
    
  /*   printf("\n. disk->time-lin1->prev->time: %f",disk->time-lin1->prev->disk->time); */
  /*   printf("\n. disk->time-lin2->prev->time: %f",disk->time-lin2->prev->disk->time); */
    
  /*   /\* Go back in time *\/ */
  /*   curr_t = 0.0; */
  /*   n_evts = 0; */
  /*   coal_evt = 0; */
  /*   have_coal = 0; */
  /*   T = disk->time; */
  /*   do */
  /*     { */
  /*       if(disk->ldsk) n_evts++; */

  /*       if(disk->ldsk && disk->ldsk == lin1->prev) lin1 = lin1->prev; */
  /*       if(disk->ldsk && disk->ldsk == lin2->prev) lin2 = lin2->prev; */
        
  /*       if(lin1 == lin2) */
  /*         { */
  /*           have_coal = 1; */
  /*           coal_evt  = n_evts; */
  /*           break; */
  /*         } */
        
  /*       disk = disk->prev; */
  /*       t=disk?disk->time:0.0; */
  /*       curr_t = (T-t); */
               
  /*       /\* printf("\n>>> coal %d n_evts: %d curr_t: %f disk->time: %f T: %f one_gen: %f", *\/ */
  /*       /\*        have_coal, *\/ */
  /*       /\*        n_evts, *\/ */
  /*       /\*        curr_t, *\/ */
  /*       /\*        t, *\/ */
  /*       /\*        T, *\/ */
  /*       /\*        one_gen); *\/ */
  /*     } */
  /*   /\* while(disk && n_evts < 2); *\/ */
  /*   while(disk && curr_t < one_gen); */

  /*   PhyML_Printf("\n. @ coal : %d %d %f",have_coal,coal_evt,curr_t); */
  /*   /\* Exit("\n"); *\/ */
  /* } */


  
  
  /* Allocate coordinates for all the tips first (will grow afterwards) */
  ldsk_a_samp = (t_ldsk **)mCalloc(n_otu,sizeof(t_ldsk *));
  ldsk_a_tips = (t_ldsk **)mCalloc(n_otu,sizeof(t_ldsk *));
  ldsk_a_tmp  = (t_ldsk **)mCalloc(n_otu,sizeof(t_ldsk *));

  while(disk->next) disk = disk->next;
  ldsk_a_pop = disk->ldsk_a;

  
  /* /\* Sample n_otu individuals uniformly at random... *\/   */
  /* while(disk->next) disk = disk->next;  */
  /* permut = Permutate(disk->n_ldsk_a); */
  /* for(i=0;i<n_otu;i++) ldsk_a_samp[i] = ldsk_a_pop[permut[i]]; */
  /* Free(permut); */
  /* n_sampled_demes = 1; */
  
  
  /* ... or take samples in random polygons */
  n_poly = n_sites;
  is_sampled = (short int *)mCalloc(n_poly,sizeof(short int));

  do
    {
      poly = (t_poly **)mCalloc(n_poly,sizeof(t_poly *));
      for(i=0;i<n_poly;++i) poly[i] = Rpoly(3); /* triangles */
      for(i=0;i<n_poly;++i)
        {
          for(j=0;j<poly[i]->n_poly_vert;++j)
            {
              poly[i]->poly_vert[j]->lonlat[0] *= (mmod->lim_up->lonlat[0]-mmod->lim_do->lonlat[0])*0.5;
              poly[i]->poly_vert[j]->lonlat[1] *= (mmod->lim_up->lonlat[1]-mmod->lim_do->lonlat[1])*0.5;
            }
          
          max_x = 0.0;
          max_y = 0.0;
          for(j=0;j<poly[i]->n_poly_vert;++j)
            {
              if(poly[i]->poly_vert[j]->lonlat[0] > max_x) max_x = poly[i]->poly_vert[j]->lonlat[0];
              if(poly[i]->poly_vert[j]->lonlat[1] > max_y) max_y = poly[i]->poly_vert[j]->lonlat[1];
            }
          
          trans_x = Uni()*(mmod->lim_up->lonlat[0] - max_x);
          trans_y = Uni()*(mmod->lim_up->lonlat[1] - max_y);
          
          for(j=0;j<poly[i]->n_poly_vert;++j)
            {
              poly[i]->poly_vert[j]->lonlat[0] += trans_x;
              poly[i]->poly_vert[j]->lonlat[1] += trans_y;
            }
        }
      
      for(i=0;i<n_otu;++i) ldsk_a_samp[i] = NULL;

      for(i=0;i<n_poly;++i) is_sampled[i] = NO;

      permut = Permutate(n_poly);
  
      sample_size = 0;
      for(i=0;i<curr_pop_size;++i)
        {
          for(j=0;j<n_poly;++j)
            {
              if(Is_In_Polygon(ldsk_a_pop[i]->coord,poly[permut[j]]) == YES)
                {
                  char *s;
                  int k;
                                   
                  s = (char *)mCalloc((int)strlen(ldsk_a_pop[i]->coord->id)+1+50,sizeof(char));
                  For(k,(int)strlen(ldsk_a_pop[i]->coord->id)+1+20) s[k]='\0';
                  sprintf(s,"%d_",i);
                  strcat(s,ldsk_a_pop[i]->coord->id);
                  Free(ldsk_a_pop[i]->coord->id);
                  strcat(s,"_deme");
                  sprintf(s+strlen(s),"%d",permut[j]);
                  ldsk_a_pop[i]->coord->id = s;


                  ldsk_a_samp[sample_size] = ldsk_a_pop[i];
                  sample_size++;
                  PhyML_Printf("\n@ Coord: %f %f %s %p",
                               ldsk_a_samp[sample_size-1]->coord->lonlat[0],
                               ldsk_a_samp[sample_size-1]->coord->lonlat[1],
                               ldsk_a_pop[i]->coord->id,ldsk_a_pop[i]);

                  is_sampled[permut[j]] = YES;
                  break;
                }
            }
          if(sample_size == n_otu) break;
        }

      Free(permut);
      
      if(i == curr_pop_size)
        {
          for(j=0;j<n_poly;j++) Free_Poly(poly[j]);
          Free(poly);
          /* PhyML_Printf("\n. Not enough individuals in polygon(s) (only %d found).",sample_size); */
          /* Generic_Exit(__FILE__,__LINE__,__FUNCTION__);       */
        }
      else break;
    }
  while(1);

  n_sampled_demes = 0;
  for(i=0;i<n_poly;++i) if(is_sampled[i] == YES) n_sampled_demes++;

  for(i=0;i<n_otu;++i) ldsk_a_tips[i] = ldsk_a_samp[i];

  area = Make_Sarea(n_sampled_demes);
  area->n_poly = n_sampled_demes;
  n_sampled_demes = 0;
  for(i=0;i<n_poly;++i) if(is_sampled[i] == YES) area->a_poly[n_sampled_demes++] = poly[i];


  for(i=0;i<area->n_poly;++i) 
    {
      /* PhyML_Printf("\n@ Poly %3d area = %f",i,Area_Of_Poly_Monte_Carlo(area->a_poly[i],mmod->lim)); */

      for(j=0;j<area->a_poly[i]->n_poly_vert;++j)
        {
          PhyML_Printf("\n@ Poly %3d point %d (x,y) = (%f,%f)",
                       i,
                       j,
                       area->a_poly[i]->poly_vert[j]->lonlat[0],
                       area->a_poly[i]->poly_vert[j]->lonlat[1]);
        }
    }

  tree->young_disk       = disk;
  disk->ldsk_a           = ldsk_a_tips;
  disk->mmod             = tree->mmod;
  disk->centr->lonlat[0] = .5*(tree->mmod->lim_up->lonlat[0]-tree->mmod->lim_do->lonlat[0]);
  disk->centr->lonlat[1] = .5*(tree->mmod->lim_up->lonlat[1]-tree->mmod->lim_do->lonlat[1]);      
  disk->n_ldsk_a         = n_otu;  

  int n_discs = 0;
  
  n_lineages = n_otu;
  do
    {
      n_remain = 0;
      for(i=0;i<n_lineages;++i) 
        {

          if((disk->prev->ldsk != NULL) && (disk->prev->ldsk == ldsk_a_samp[i]->prev)) /* Coalescent event is sampled */
            {
              PHYREX_Make_Lindisk_Next(disk->prev->ldsk);
              disk->prev->ldsk->next[disk->prev->ldsk->n_next-1] = ldsk_a_samp[i];
            }
          else
            {
              ldsk_a_tmp[n_remain] = ldsk_a_samp[i];
              n_remain++;
            }
        }

      for(i=0;i<n_remain;i++) ldsk_a_samp[i] = ldsk_a_tmp[i];
      if((disk->prev->ldsk != NULL) && (disk->prev->ldsk->n_next > 0)) ldsk_a_samp[i] = disk->prev->ldsk;

      if((disk->prev->ldsk != NULL) && (disk->prev->ldsk->n_next > 0))
        {
          n_lineages -= (disk->prev->ldsk->n_next);
          n_lineages += 1;
        }

      if(n_lineages != n_remain+(disk->prev->ldsk && disk->prev->ldsk->n_next>0)?1:0) 
        {
          PhyML_Fprintf(stderr,"\n. n_lineages: %d n_remain: %d n_next: %d",
                        n_lineages,
                        n_remain,
                        disk->prev->ldsk->n_next);
          Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
        }

      /* None of the sampled lineages was hit */
      if((disk->prev->ldsk != NULL) && (disk->prev->ldsk->n_next == 0))
        {
          /* printf("\n. free %s @ %f",disk->ldsk->coord->id,disk->time); */
          Free_Ldisk(disk->prev->ldsk);
          disk->prev->ldsk = NULL;
        }

      disk = disk->prev;
      
      if(disk->prev == NULL)
        {
          PhyML_Fprintf(stderr,"\n. # lineages left: %d",n_remain);
          PhyML_Fprintf(stderr,"\n. Sample has not coalesced completely.");
          fflush(NULL);
          Exit("\n");
        }
    }
  while(n_lineages > 1);

  
  /* for(i=0;i<n_otu;i++) printf("\n> %s",tree->young_disk->ldsk_a[i]->coord->id); */

  disk->prev = NULL;

  disk = tree->young_disk;
  tree_height = disk->time;
  n_discs = 0;
  while(disk)
    {
      disk->time -= tree_height;
      disk = disk->prev;
      n_discs++;
    }

  Free(ldsk_a_tmp);
  Free(ldsk_a_samp);
  Free(ldsk_a_pop);
  Free(is_sampled);
  
  /* PHYREX_Print_Struct('#',tree); */
  /* Exit("\n"); */

  return(area);
 
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
/* Test whether coord is in disk. Will actually only works in disk */
/* is a rectangle... */

int PHYREX_Is_In_Disk(t_geo_coord *coord, t_dsk *disk, t_phyrex_mod *mmod)
{
  int i;

  assert(disk->centr->dim);

  if(mmod->name == PHYREX_UNIFORM)
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
  else if(mmod->name == PHYREX_NORMAL)
    {
      return(YES);
    }

  return(-1);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl PHYREX_Lk(t_tree *tree)
{
  t_dsk *disk;
  int n_evt;
  
  assert(!tree->young_disk->next);
  assert(tree->young_disk->prev);
  
  tree->mmod->c_lnL = 0.0;

  tree->mmod->c_lnL += PHYREX_LnPrior_Radius(tree);
  
  if(isinf(tree->mmod->c_lnL) || isnan(tree->mmod->c_lnL)) 
    {
      tree->mmod->c_lnL = UNLIKELY;
      return tree->mmod->c_lnL;
    }

  PHYREX_Update_Lindisk_List(tree);
  tree->mmod->c_lnL += PHYREX_Lk_Core(tree->young_disk,tree);
  
  n_evt = 0;
  disk = tree->young_disk->prev;
  do
    {
      if(disk->age_fixed == NO) n_evt++;

      assert(disk);
      if(disk->time > disk->next->time)
        {
          PhyML_Printf("\n. Anomaly in ordering of disk times (disk->time: %f disk->next->time: %f)",disk->time,disk->next->time);
          PhyML_Printf("\n. Currently doing following move: %s",tree->mcmc->move_name[tree->mcmc->move_idx]);
          tree->mmod->c_lnL = UNLIKELY;
          return tree->mmod->c_lnL;
        }

      tree->mmod->c_lnL += PHYREX_Lk_Core(disk,tree);
      if(disk->prev == NULL) break;
      disk = disk->prev;      
    }
  while(1);

  /* tree->mmod->c_lnL += PHYREX_Lk_Time_Component(tree); */

  tree->mmod->c_lnL += n_evt * log(tree->mmod->lbda) - tree->mmod->lbda * fabs(tree->young_disk->time - disk->time);
    
  if(isinf(tree->mmod->c_lnL) || isnan(tree->mmod->c_lnL)) tree->mmod->c_lnL = UNLIKELY;

  return(tree->mmod->c_lnL);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl PHYREX_Lk_Core(t_dsk *disk, t_tree *tree)
{
  phydbl lnL,log_prob_hit,log_mu;
  int was_hit,i,j,err;
  phydbl two_theta_two;

  two_theta_two = 2.*pow(tree->mmod->rad,2);
  lnL           = 0.0;
  log_mu        = log(tree->mmod->mu);
  was_hit       = (disk->ldsk != NULL);

  if(disk == tree->young_disk) return 0.0;
  if(disk->age_fixed == YES) return 0.0;
  
  for(i=0;i<disk->n_ldsk_a;++i)
    if(disk->ldsk_a[i]->prev->disk->time > disk->time)
      return UNLIKELY;
  
  if(disk->ldsk != NULL)
    {
      for(i=0;i<disk->n_ldsk_a;++i) if(disk->ldsk_a[i]->prev == disk->ldsk) break;
      if(i == disk->n_ldsk_a) return UNLIKELY;
      /* assert(i != disk->n_ldsk_a); */
    }

  for(i=0;i<disk->n_ldsk_a;++i)
    {
      if(PHYREX_Is_In_Ldscape(disk->ldsk_a[i],tree->mmod) == NO)
        {
          for(j=0;j<tree->mmod->n_dim;j++) PhyML_Fprintf(stdout,"\n. Found a lineage outside habitat (%20f)",disk->ldsk_a[i]->coord->lonlat[j]);
          return(UNLIKELY);     
        }
            
      log_prob_hit = log_mu;
      for(j=0;j<tree->mmod->n_dim;j++)
        log_prob_hit += -pow(disk->ldsk_a[i]->coord->lonlat[j] - disk->centr->lonlat[j],2)/two_theta_two;

      
      if(disk->ldsk_a[i]->prev == disk->ldsk) // disk->ldsk_a[i] was hit
        {
          lnL += log_prob_hit;
        }
      else // disk->ldsk_a[i] was not hit
        {
          lnL += log(1. - exp(log_prob_hit));
        }
    }

  /* a hit occurred */
  if(was_hit == TRUE)
    {
      err = NO;
      for(j=0;j<tree->mmod->n_dim;j++) lnL += Log_Dnorm_Trunc(disk->ldsk->coord->lonlat[j],
                                                              disk->centr->lonlat[j],
                                                              tree->mmod->rad,
                                                              tree->mmod->lim_do->lonlat[j],
                                                              tree->mmod->lim_up->lonlat[j],&err);
    }

  /* Likelihood for the disk center */
  for(j=0;j<tree->mmod->n_dim;j++) lnL += log(1./(tree->mmod->lim_up->lonlat[j]-tree->mmod->lim_do->lonlat[j]));
  
  return(lnL);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// Warning: the calculation below does not incoporate time information
// since things get messy when considering serial samples
phydbl PHYREX_Lk_Range(t_dsk *young, t_dsk *old, t_tree *tree)
{
  t_dsk *disk;
  phydbl lnL;
  int n_evt;
  
  return(PHYREX_Lk(tree));
  
  assert(young);

  lnL  = 0.0;

  PHYREX_Update_Lindisk_List_Core(young,tree);
  lnL += PHYREX_Lk_Core(young,tree);

  n_evt = 0;
  disk = young->prev;
  do
    {
      if(disk->age_fixed == NO) n_evt++;      
      assert(disk);
      if(disk->time > disk->next->time) return UNLIKELY;
      PHYREX_Update_Lindisk_List_Core(disk,tree);
      lnL += PHYREX_Lk_Core(disk,tree);
      if(disk == old) break;
      disk = disk->prev;
    }
  while(disk);

  /* lnL += PHYREX_Lk_Time_Component(tree); */
  lnL += n_evt * log(tree->mmod->lbda) - tree->mmod->lbda * fabs(young->time - old->time);
 
  return(lnL);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl PHYREX_Lk_Core_Range(t_dsk *young, t_dsk *old, t_tree *tree)
{
  t_dsk *disk;
  phydbl lnL;
  
  assert(young);
  
  lnL  = 0.0;

  PHYREX_Update_Lindisk_List_Core(young,tree);
  lnL += PHYREX_Lk_Core(young,tree);

  disk = young->prev;
  do
    {
      assert(disk);
      PHYREX_Update_Lindisk_List_Core(disk,tree);
      lnL += PHYREX_Lk_Core(disk,tree);
      if(disk == old) break;
      disk = disk->prev;
    }
  while(1);
      

  return(lnL);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

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
  
  MCMC_Complete_MCMC(mcmc,tree);

  n_vars = 12;

  res = (phydbl *)mCalloc(tree->mcmc->chain_len / tree->mcmc->sample_interval * n_vars,sizeof(phydbl));


  /* Starting parameter values */
  
  /* tree->mmod->lbda = Uni()*(tree->mmod->max_lbda - tree->mmod->min_lbda) + tree->mmod->min_lbda; */
  tree->mmod->lbda = (phydbl)(tree->n_otu) / fabs(PHYREX_Tree_Height(tree));
  tree->mmod->mu   = Uni()*(tree->mmod->max_mu - tree->mmod->min_mu) + tree->mmod->min_mu;
  tree->mmod->rad  = Uni()*(tree->mmod->max_rad - tree->mmod->min_rad) + tree->mmod->min_rad;

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
      PhyML_Fprintf(fp_stats,"\n# true rad: %f",PHYREX_Update_Radius(tree));
      PhyML_Fprintf(fp_stats,"\n# true sigsq: %f",tree->mmod->sigsq);
      PhyML_Fprintf(fp_stats,"\n# true neigh. size: %f",PHYREX_Neighborhood_Size(tree));
      PhyML_Fprintf(fp_stats,"\n# fst-based estimate of neighborhood size: %f",PHYREX_Neighborhood_Size_Regression(tree));
      PhyML_Fprintf(fp_stats,"\n# nucleotide diversity: %f",Nucleotide_Diversity(tree->data));
      PhyML_Fprintf(fp_stats,"\n# length of a generation: %G time units",PHYREX_Generation_Length(tree));
      PhyML_Fprintf(fp_stats,"\n# clock rate: %G subst. per time unit",tree->rates->clock_r);
      
      
      PhyML_Fprintf(fp_stats,"\n# after rand glnL: %f alnL: %f",tree->mmod->c_lnL,tree->c_lnL);
      PhyML_Fprintf(fp_stats,"\n# ninter: %d",PHYREX_Total_Number_Of_Intervals(tree));
      PhyML_Fprintf(fp_stats,"\n# ncoal: %d",PHYREX_Total_Number_Of_Coal_Disks(tree));
      PhyML_Fprintf(fp_stats,"\n# nhits: %d",PHYREX_Total_Number_Of_Hit_Disks(tree));
      PhyML_Fprintf(fp_stats,"\n# start lbda: %f",tree->mmod->lbda);
      PhyML_Fprintf(fp_stats,"\n# start mu: %f",tree->mmod->mu);
      PhyML_Fprintf(fp_stats,"\n# start rad: %f",tree->mmod->rad);
      fflush(NULL);
      
      
      PhyML_Fprintf(fp_stats,"\n");
      PhyML_Fprintf(fp_stats,"%s\t","sample");
      PhyML_Fprintf(fp_stats,"%s\t","lnP");
      PhyML_Fprintf(fp_stats,"%s\t","alnL");
      PhyML_Fprintf(fp_stats,"%s\t","glnL");
      PhyML_Fprintf(fp_stats,"%s\t","clock");
      PhyML_Fprintf(fp_stats,"%s\t","lbda");
      PhyML_Fprintf(fp_stats,"%s\t","mu");
      PhyML_Fprintf(fp_stats,"%s\t","rad");

      PhyML_Fprintf(fp_stats,"%s\t","neigh");
      PhyML_Fprintf(fp_stats,"%s\t","rhoe");

      PhyML_Fprintf(fp_stats,"%s\t","sigsq");
      PhyML_Fprintf(fp_stats,"%s\t","sigsqobs");
      PhyML_Fprintf(fp_stats,"%s\t","dispdist");
      PhyML_Fprintf(fp_stats,"%s\t","nInt");
      PhyML_Fprintf(fp_stats,"%s\t","nCoal");
      PhyML_Fprintf(fp_stats,"%s\t","nHit");
      PhyML_Fprintf(fp_stats,"%s\t","rootTime");
      PhyML_Fprintf(fp_stats,"%s\t","rootLon");
      PhyML_Fprintf(fp_stats,"%s\t","rootLat");
      for(int i=0;i<Scalar_Len(tree->mod->kappa);++i) PhyML_Fprintf(fp_stats,"tstv%d\t",i);
      PhyML_Fprintf(fp_stats,"%s\t","alpha");
      PhyML_Fprintf(fp_stats,"%s\t","MeanBr");
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
      PhyML_Fprintf(fp_stats,"%s\t","tuneLbda");
      PhyML_Fprintf(fp_stats,"%s\t","tuneRad");
      PhyML_Fprintf(fp_stats,"%s\t","tuneMu");
      PhyML_Fprintf(fp_stats,"%s\t","tuneIndelDisk");
      PhyML_Fprintf(fp_stats,"%s\t","tuneIndelHit");
      PhyML_Fprintf(fp_stats,"%s\t","tuneLdskGivenDisk");
      PhyML_Fprintf(fp_stats,"%s\t","tuneIndelDiskSerial");
      PhyML_Fprintf(fp_stats,"%s\t","tuneIndelHitSerial");
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


      
      /* tree->mmod->lbda = 1.0; */
      /* tree->mmod->mu   = 0.5; */
      /* tree->mmod->rad  = 1.5; */

      /* if(tree->mcmc->run == 0) */
        /* { */
          /* PHYREX_Strip_And_Reconnect_Tree(tree); */
          /* PHYREX_Simulate_Backward_Core(tree->young_disk,NO,tree); */
          /* PHYREX_Ldsk_To_Tree(tree); */
        /* } */
            
      if(!strcmp(tree->mcmc->move_name[move],"phyrex_lbda"))
        MCMC_PHYREX_Lbda(tree);
      
      if(!strcmp(tree->mcmc->move_name[move],"phyrex_mu"))
        MCMC_PHYREX_Mu(tree);
      
      if(!strcmp(tree->mcmc->move_name[move],"phyrex_rad"))
        MCMC_PHYREX_Radius(tree);
    
      
      if(!strcmp(tree->mcmc->move_name[move],"phyrex_indel_disk"))
        MCMC_PHYREX_Indel_Disk(tree);

      /* if(!strcmp(tree->mcmc->move_name[move],"phyrex_indel_hit")) */
      /*   MCMC_PHYREX_Indel_Hit(tree); */

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
      
      
      if(!(tree->mcmc->run%tree->mcmc->sample_interval))
        {

          /* FILE *fp; */
          /* fp = fopen("file.txt","w"); */
          /* PHYREX_Check_Point(fp,tree); */
          /* fclose(fp); */
          /* Exit("\n"); */

          
          disk = tree->young_disk;
          while(disk->prev) disk = disk->prev;

          PhyML_Fprintf(fp_stats,"\n");
          PhyML_Fprintf(fp_stats,"%6d\t",tree->mcmc->run);
          PhyML_Fprintf(fp_stats,"%g\t",tree->c_lnL+tree->mmod->c_lnL+tree->rates->c_lnL_rates);
          PhyML_Fprintf(fp_stats,"%g\t",tree->c_lnL);
          PhyML_Fprintf(fp_stats,"%g\t",tree->mmod->c_lnL);
          PhyML_Fprintf(fp_stats,"%g\t",tree->rates->clock_r);
          PhyML_Fprintf(fp_stats,"%g\t",tree->mmod->lbda);
          PhyML_Fprintf(fp_stats,"%g\t",tree->mmod->mu);
          PhyML_Fprintf(fp_stats,"%g\t",tree->mmod->rad);
          PhyML_Fprintf(fp_stats,"%g\t",PHYREX_Neighborhood_Size(tree));
          PhyML_Fprintf(fp_stats,"%g\t",PHYREX_Effective_Density(tree));
          PhyML_Fprintf(fp_stats,"%g\t",PHYREX_Update_Sigsq(tree));
          PhyML_Fprintf(fp_stats,"%g\t",PHYREX_Realized_Sigsq(tree));
          PhyML_Fprintf(fp_stats,"%g\t",PHYREX_Realized_Dispersal_Dist(tree));
          PhyML_Fprintf(fp_stats,"%d\t",PHYREX_Total_Number_Of_Intervals(tree));
          PhyML_Fprintf(fp_stats,"%d\t",PHYREX_Total_Number_Of_Coal_Disks(tree));
          PhyML_Fprintf(fp_stats,"%d\t",PHYREX_Total_Number_Of_Hit_Disks(tree));
          PhyML_Fprintf(fp_stats,"%g\t",disk->time);
          PhyML_Fprintf(fp_stats,"%g\t",disk->ldsk->coord->lonlat[0]);
          PhyML_Fprintf(fp_stats,"%g\t",disk->ldsk->coord->lonlat[1]);
          Output_Scalar_Dbl(tree->mod->kappa,"\t",fp_stats);
          PhyML_Fprintf(fp_stats,"%g\t",tree->mod->ras->alpha->v);
          /* for(int i=0;i<2*tree->n_otu-1;++i) PhyML_Fprintf(fp_stats,"%g\t",tree->rates->br_r[i]); */
          PhyML_Fprintf(fp_stats,"%g\t",RATES_Get_Mean_Rate_In_Subtree(tree->n_root,tree));
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
          PhyML_Fprintf(fp_stats,"%g\t",tree->mcmc->tune_move[tree->mcmc->num_move_phyrex_lbda]);
          PhyML_Fprintf(fp_stats,"%g\t",tree->mcmc->tune_move[tree->mcmc->num_move_phyrex_rad]);
          PhyML_Fprintf(fp_stats,"%g\t",tree->mcmc->tune_move[tree->mcmc->num_move_phyrex_mu]);
          PhyML_Fprintf(fp_stats,"%g\t",tree->mcmc->tune_move[tree->mcmc->num_move_phyrex_indel_disk]);
          PhyML_Fprintf(fp_stats,"%g\t",tree->mcmc->tune_move[tree->mcmc->num_move_phyrex_indel_hit]);
          PhyML_Fprintf(fp_stats,"%g\t",tree->mcmc->tune_move[tree->mcmc->num_move_phyrex_ldsk_given_disk]);
          PhyML_Fprintf(fp_stats,"%g\t",tree->mcmc->tune_move[tree->mcmc->num_move_phyrex_indel_disk_serial]);
          PhyML_Fprintf(fp_stats,"%g\t",tree->mcmc->tune_move[tree->mcmc->num_move_phyrex_indel_hit_serial]);

          
          
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
          
          /* if(tree->mcmc->sample_num == 0) */
          /*   { */
          /*     PhyML_Fprintf(fp_tree,"\n#NEXUS"); */
          /*     PhyML_Fprintf(fp_tree,"\nBEGIN TREES;"); */
          /*   } */
          /* else */
          /*   { */
          /*     fseek(fp_tree,-5,SEEK_CUR); */
          /*   } */

          
          PHYREX_Ldsk_To_Tree(tree);
          TIMES_Time_To_Bl(tree);
          tree->bl_ndigits = 3;
          /* RATES_Update_Cur_Bl(tree); */
          char *s = Write_Tree(tree,NO);
          PhyML_Fprintf(fp_tree,"\ntree %d [&lnP=%f] = [&R]  %s",tree->mcmc->sample_num,tree->c_lnL,s);
          Free(s);
          PhyML_Fprintf(fp_tree,"\nEND;");

          fflush(NULL);
                
          tree->mcmc->sample_num++;
        }
      
      (void)signal(SIGINT,MCMC_Terminate);
    }
  while(tree->mcmc->run < tree->mcmc->chain_len);

  fclose(fp_tree);
  fclose(fp_stats);
  /* fclose(fp_summary); */

  return(res);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl PHYREX_Wrap_Lk(t_edge *b, t_tree *tree, supert_tree *stree)
{
  return PHYREX_Lk(tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void PHYREX_New_Traj(t_dsk *start, t_dsk *end, t_tree *tree)
{
  t_dsk *disk;
  int i,j;
  int n_hit_up,n_hit_tot;
  /* phydbl min_up, min_do; */
  /* phydbl max_up, max_do; */



  for(i=0;i<start->n_ldsk_a;i++)
    {      
      /* printf("\n<><><>"); */

      disk = end;
      while(disk != start)
        {
          /* disk->ldsk_a[i]->min_coord = GEO_Make_Geo_Coord(disk->mmod->n_dim); */
          /* disk->ldsk_a[i]->max_coord = GEO_Make_Geo_Coord(disk->mmod->n_dim);                     */
          disk = disk->prev;
        }

      disk = end;
      n_hit_tot = 0;
      while(disk != start)
        {
          if(disk->ldsk_a[i]->is_hit == YES) n_hit_tot++;
          disk = disk->prev;
        }


      /* printf("\n. n_hit_tot: %d",n_hit_tot); */
      /* fflush(NULL); */

      /* printf("\n. Start disk name: %s lindisk %s at %.2f [%.2f %.2f] end disk: %s lindisk %s at %.2f [%.2f %.2f]", */
      /*        start->id, */
      /*        start->ldsk_a[i]->coord->id, */
      /*        start->time, */
      /*        start->ldsk_a[i]->coord->lonlat[0], */
      /*        start->ldsk_a[i]->coord->lonlat[1], */
      /*        end->id, */
      /*        end->ldsk_a[i]->coord->id, */
      /*        end->time, */
      /*        end->ldsk_a[i]->coord->lonlat[0], */
      /*        end->ldsk_a[i]->coord->lonlat[1] */
      /*        ); */
      /* disk = end; */
      /* while(disk != start) */
      /*   { */
      /*     printf("\n. %d %s %.2f %.2f", */
      /*            i, */
      /*            disk->ldsk_a[i]->coord->id, */
      /*            disk->ldsk_a[i]->coord->lonlat[0], */
      /*            disk->ldsk_a[i]->coord->lonlat[1]); */
      /*     disk = disk->prev; */
      /*   } */
      
      n_hit_up = 0;
      disk = end;
      while(disk->prev != start)
        {
          if(disk->ldsk_a[i]->is_hit == YES)
            {              
              n_hit_up++;

              for(j=0;j<disk->mmod->n_dim;j++)
                {
                  /* min_up = end->ldsk_a[i]->coord->lonlat[j] - n_hit_up * 2. * disk->mmod->rad; */
                  /* min_do = start->ldsk_a[i]->coord->lonlat[j] - (n_hit_tot - n_hit_up) * 2. * disk->mmod->rad; */
                  /* disk->ldsk_a[i]->prev->min_coord->lonlat[j] =  */
                  /*   MAX(MAX(MAX(min_up,min_do),0.0),disk->prev->centr->lonlat[j] - disk->mmod->rad); */

                  /* max_up = end->ldsk_a[i]->coord->lonlat[j] + n_hit_up * 2. * disk->mmod->rad; */
                  /* max_do = start->ldsk_a[i]->coord->lonlat[j] + (n_hit_tot - n_hit_up) * 2. * disk->mmod->rad; */
                  /* disk->ldsk_a[i]->prev->max_coord->lonlat[j] =  */
                  /*   MIN(MIN(MIN(max_up,max_do),disk->mmod->lim->lonlat[j]),disk->prev->centr->lonlat[j] + disk->mmod->rad); */

                  /* printf("\n. curr: %s min_up: %.2f min_do: %.2f max_up: %.2f max_do: %.2f %d %d start: %.2f [%.2f %.2f]", */
                  /*        disk->ldsk_a[i]->prev->coord->id, */
                  /*        min_up,min_do,max_up,max_do, */
                  /*        n_hit_tot,n_hit_up, */
                  /*        start->ldsk_a[i]->coord->lonlat[j], */
                  /*        disk->ldsk_a[i]->prev->min_coord->lonlat[j], */
                  /*        disk->ldsk_a[i]->prev->max_coord->lonlat[j]); */
                }
            }
          disk = disk->prev;
        }
      
      disk = end;
      while(disk->prev != start)
        {
          if(disk->ldsk_a[i]->is_hit == YES)
            {
              /* for(j=0;j<disk->mmod->n_dim;j++) */
                /* disk->ldsk_a[i]->prev->coord->lonlat[j] =  */
                /* Uni()* */
                /* (disk->ldsk_a[i]->prev->max_coord->lonlat[j]  - */
                /*  disk->ldsk_a[i]->prev->min_coord->lonlat[j]) + */
                /* disk->ldsk_a[i]->prev->min_coord->lonlat[j]; */
            }
          else
            {
              for(j=0;j<disk->mmod->n_dim;j++)
                disk->ldsk_a[i]->prev->coord->lonlat[j] = 
                disk->ldsk_a[i]->coord->lonlat[j];
            }

          disk = disk->prev;
        }
      
      disk = end;
      while(disk != start)
        {
          /* Free_Geo_Coord(disk->ldsk_a[i]->min_coord); */
          /* Free_Geo_Coord(disk->ldsk_a[i]->max_coord); */
          disk = disk->prev;
        }
    }

  /* gtk_widget_queue_draw(tree->draw_area); */
  
  /* for(i=0;i<start->n_ldsk_a;i++) */
  /*   { */
  /*     printf("\n<><><>"); */
  /*     disk = end; */
  /*     while(disk != start) */
  /*       { */
  /*         printf("\nx %s %.2f %.2f centr: %.2f %.2f", */
  /*                disk->ldsk_a[i]->coord->id, */
  /*                disk->ldsk_a[i]->coord->lonlat[0], */
  /*                disk->ldsk_a[i]->coord->lonlat[1], */
  /*                disk->centr->lonlat[0], */
  /*                disk->centr->lonlat[1]); */
  /*         disk = disk->prev; */
  /*       } */
  /*     fflush(NULL); */
  /*   } */
  
  /* PHYREX_Lk(tree); */
  /* printf("\n. >> Lk: %f rad: %f",tree->mmod->c_lnL,tree->mmod->rad); */
  /* sleep(5); */
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

/* Insert disk event based on its time. Insertion above the root is not permitted  
*/
void PHYREX_Insert_Disk(t_dsk *ins, t_tree *tree)
{
  t_dsk *disk;

  assert(!(ins == NULL));
  
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
  /* phydbl main_mass,l; */

  log_dens  = 0.0;
  /* main_mass = 1. - mmod->soft_bound_area; */

  up   = MIN(disk->centr->lonlat[0] + mmod->rad, mmod->lim_up->lonlat[0]);
  down = MAX(disk->centr->lonlat[0] - mmod->rad, mmod->lim_do->lonlat[0]);
  rght = MIN(disk->centr->lonlat[1] + mmod->rad, mmod->lim_up->lonlat[1]);
  left = MAX(disk->centr->lonlat[1] - mmod->rad, mmod->lim_do->lonlat[1]);


  /* l = main_mass / (.5*(1.-main_mass)*(up - down)); */
  /* if(ldsk->coord->lonlat[0] < down)    log_dens += log(.5*(1.-main_mass)) + log(l) - l*(down - ldsk->coord->lonlat[0]); */
  /* else if(ldsk->coord->lonlat[0] > up) log_dens += log(.5*(1.-main_mass)) + log(l) - l*(ldsk->coord->lonlat[0] - up); */
  /* else log_dens += log(main_mass) - log(up - down); */

  /* l = main_mass / (.5*(1.-main_mass)*(rght - left)); */
  /* if(ldsk->coord->lonlat[1] < left)      log_dens += log(.5*(1.-main_mass)) + log(l) - l*(left - ldsk->coord->lonlat[1]); */
  /* else if(ldsk->coord->lonlat[1] > rght) log_dens += log(.5*(1.-main_mass)) + log(l) - l*(ldsk->coord->lonlat[1] - rght); */
  /* else log_dens += log(main_mass) - log(rght - left); */


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
      ldsk->coord->lonlat[i] = Rnorm_Trunc(disk->centr->lonlat[i],mmod->rad,mmod->lim_do->lonlat[i],mmod->lim_up->lonlat[i],&err);
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
  tree->mmod->c_ln_prior_sigsq = 
    log(tree->mmod->prior_param_sigsq) - 
    tree->mmod->prior_param_sigsq*tree->mmod->sigsq; 
  return(tree->mmod->c_ln_prior_sigsq);
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

  if(!disk->ldsk || tree->mmod->name == PHYREX_NORMAL)
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
  tree->rates->nd_t[a->num] = ldsk->disk->time;

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
              
              tree->rates->nd_t[new_parent->num] = ldsk->disk->time;

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
  switch(tree->mmod->name)
    {
    case PHYREX_UNIFORM:
      {
        PHYREX_Runif_Rectangle_Overlap(a_disk->ldsk,a_disk,tree->mmod);
        break;
      }
    case PHYREX_NORMAL:
      {
        PHYREX_Rnorm_Trunc(a_disk->ldsk,a_disk,tree->mmod);
        break;
      }
    }
  
  a_disk->ldsk->nd = tree->n_root;
  a_disk->time = tree->rates->nd_t[tree->n_root->num];

  /* Inflate_Times_To_Get_Reasonnable_Edge_Lengths(1.E-3,tree); */
  Get_Node_Ranks_From_Times(tree);

  PHYREX_Tree_To_Ldsk_Post(tree->n_root,tree->n_root->v[1],a_disk,tree);
  PHYREX_Tree_To_Ldsk_Post(tree->n_root,tree->n_root->v[2],a_disk,tree);

  // Create a doubly-chained list of disks
  disk = a_disk;
  disk->time = tree->rates->nd_t[tree->n_root->num];
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
          disk->next->time = tree->rates->nd_t[n->rk_next->num];
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
      if(tree->rates->nd_t[d->num] > tree->rates->nd_t[a->num])
        {
          t_dsk *d_disk;
          
          PHYREX_Make_Lindisk_Next(a_disk->ldsk);

          // Make and initialize descendent disk
          d_disk = PHYREX_Make_Disk_Event(tree->mmod->n_dim,tree->n_otu);
          assert(d_disk);
          PHYREX_Init_Disk_Event(d_disk,tree->mmod->n_dim,NULL);
          
          d_disk->time = tree->rates->nd_t[d->num];
          
          d_disk->ldsk = PHYREX_Make_Lindisk_Node(tree->mmod->n_dim);
          PHYREX_Init_Lindisk_Node(d_disk->ldsk,d_disk,tree->mmod->n_dim);
          
          // Initialize centre of event on the root disk
          d_disk->centr->lonlat[0] = Uni()*(tree->mmod->lim_up->lonlat[0]-tree->mmod->lim_do->lonlat[0])+tree->mmod->lim_do->lonlat[0];
          d_disk->centr->lonlat[1] = Uni()*(tree->mmod->lim_up->lonlat[1]-tree->mmod->lim_do->lonlat[1])+tree->mmod->lim_do->lonlat[1];      
          
          /* Its location */
          switch(tree->mmod->name)
            {
            case PHYREX_UNIFORM:
              {
                PHYREX_Runif_Rectangle_Overlap(d_disk->ldsk,d_disk,tree->mmod);
                break;
              }
            case PHYREX_NORMAL:
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
            tree->rates->nd_t[disk->prev->ldsk->nd->num] =
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
  disk->time = tree->rates->nd_t[n->num]; // Set time of youngest disk
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
      if(Are_Equal(tree->rates->nd_t[n->num],tree->rates->nd_t[n->rk_prev->num],SMALL) == NO)
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
      disk->time = tree->rates->nd_t[n->rk_prev->num];
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
  T = -tree->rates->nd_t[tree->n_root->num];
  return((phydbl)(T/n_inter));
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

/* Uses the regression technique described in Barton et al. TPB (2013)
   to estimate the size of the neighborhood when the population evolve
   according to a spatial Lambda-Fleming-Viot process.
*/
phydbl PHYREX_Neighborhood_Size_Regression(t_tree *tree)
{
  int i,j,pair;
  t_node *anc;
  phydbl *dist,min_dist;
  phydbl QA,Qr,*fst,fst0,fst_min_dist;
  phydbl cov_fst_dist,var_dist,slope;
  phydbl eps;

  eps = 1.E-10;

  QA = Mean_Identity(tree->data);

  fst  = (phydbl *)mCalloc(tree->n_otu*(tree->n_otu-1)/2,sizeof(phydbl));
  dist = (phydbl *)mCalloc(tree->n_otu*(tree->n_otu-1)/2,sizeof(phydbl));

  pair = 0;
  fst0 = 0.0;
  for(i=0;i<tree->n_otu-1;i++)
    {
      fst_min_dist = 0.0;
      min_dist = MDBL_MAX;
      for(j=i+1;j<tree->n_otu;j++)
        {
          anc = Find_Lca_Pair_Of_Nodes(tree->a_nodes[i],tree->a_nodes[j],tree);
          if(anc == NULL) 
            {
              PhyML_Fprintf(stderr,"\n. %s",Write_Tree(tree,NO));
              PhyML_Fprintf(stderr,"\n. %s %s",tree->a_nodes[i]->name,tree->a_nodes[j]->name);
              Generic_Exit(__FILE__,__LINE__,__FUNCTION__);            
            }
          
          dist[pair] = Euclidean_Dist(tree->a_nodes[i]->ldsk->coord,tree->a_nodes[j]->ldsk->coord);
          dist[pair] = log(dist[pair]);

          Qr = Pairwise_Identity(i,j,tree->data);
          fst[pair] = (Qr-QA)/(1.-QA);

          if(dist[pair] < min_dist)
            {
              min_dist     = dist[pair];
              fst_min_dist = fst[pair];
            }

          pair++;
        }
      fst0 += fst_min_dist;
    }

  fst0 /= (phydbl)(tree->n_otu);

  cov_fst_dist = Covariance(dist,fst,pair);
  var_dist = Variance(dist,pair);

  slope = cov_fst_dist / var_dist;

  Free(dist);
  Free(fst);

  return((fst0-1.+eps)/(slope+eps));
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
          PhyML_Fprintf(stderr,"\n. %s",Write_Tree(tree,NO));
          PhyML_Fprintf(stderr,"\n. %s %s",tree->a_nodes[i]->name,tree->a_nodes[j]->name);
          Generic_Exit(__FILE__,__LINE__,__FUNCTION__);            
        }
      
      PhyML_Printf("\nxxWxx %12f",tree->rates->nd_t[anc->num]);
      dist = Euclidean_Dist(tree->a_nodes[i]->ldsk->coord,tree->a_nodes[j]->ldsk->coord);
      PhyML_Printf(" %f",dist);
    }
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl PHYREX_Generation_Length(t_tree *tree)
{
  return(1./(2.*tree->mmod->mu*POW(tree->mmod->rad,2)*PI*PHYREX_Rate_Per_Unit_Area(tree)));
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl PHYREX_Neighborhood_Size(t_tree *tree)
{
  switch(tree->mmod->name)
    {
    case PHYREX_UNIFORM: { return(1./tree->mmod->mu); break; }
    case PHYREX_NORMAL:  { return(2./tree->mmod->mu); break; }
    }
  return(-1.);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl PHYREX_Update_Sigsq(t_tree *tree)
{  
  switch(tree->mmod->name)
    {
    case PHYREX_UNIFORM: { return(-1.0); break;}
    case PHYREX_NORMAL:  
      {
        return(4.*PI*
               PHYREX_Rate_Per_Unit_Area(tree) *
               pow(tree->mmod->rad,4)*
               tree->mmod->mu); 
        break; 
      }
    }
  return(-1.);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl PHYREX_Update_Radius(t_tree *tree)
{
  switch(tree->mmod->name)
    {
    case PHYREX_UNIFORM: { return(-1.0); break;}
    case PHYREX_NORMAL:  
      { 
        return(POW(tree->mmod->sigsq/(PHYREX_Rate_Per_Unit_Area(tree)*4.*PI*tree->mmod->mu),0.25));
        break; 
      }
    }
  return(-1.);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl PHYREX_Sample_Rad_From_Prior(t_tree *tree)
{
  phydbl u,h,w,lbda,mu,b,a;

  h    = (tree->mmod->lim_up->lonlat[0]-tree->mmod->lim_do->lonlat[0]);
  w    = (tree->mmod->lim_up->lonlat[1]-tree->mmod->lim_do->lonlat[1]);
  lbda = tree->mmod->lbda;
  mu   = tree->mmod->mu;
  a    = tree->mmod->min_sigsq;
  b    = tree->mmod->max_sigsq;

  u = Uni();

  return(POW(((b-a)*u+a)*h*w/(4.*PI*lbda*mu),0.25));

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
          if(fscanf(fp,"%lf",&(tree->a_nodes[i]->ldsk->coord->lonlat[0])) == EOF) break;
          if(fscanf(fp,"%lf",&(tree->a_nodes[i]->ldsk->coord->lonlat[1])) == EOF) break;          
          done[i] = YES;
        }
      else
        {
          if(!strcmp(s,"|SouthWest|") || !strcmp(s,"|southwest|") || !strcmp(s,"|Southwest|"))
            {
              found_sw = YES;
              if(fscanf(fp,"%lf",&(sw_lon)) == EOF) break;
              if(fscanf(fp,"%lf",&(sw_lat)) == EOF) break;              
            }
          else if(!strcmp(s,"|NorthEast|") || !strcmp(s,"|northeast|") || !strcmp(s,"|Northeast|"))
            {
              found_ne = YES;
              if(fscanf(fp,"%lf",&(ne_lon)) == EOF) break;
              if(fscanf(fp,"%lf",&(ne_lat)) == EOF) break;
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

phydbl PHYREX_Rate_Per_Unit_Area(t_tree *tree)
{
  int i;
  phydbl denom;

  denom = (tree->mmod->lim_up->lonlat[0]-tree->mmod->lim_do->lonlat[0]);  
  for(i=1;i<tree->mmod->n_dim;i++) denom *= (tree->mmod->lim_up->lonlat[i]-tree->mmod->lim_do->lonlat[i]);
  
  return(tree->mmod->lbda / denom);

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
/* end is older ldsk. beg is younger ldsk. Remove path between */
/* beg ldsk and end ldsk. If beg->prev == end, then beg->prev set */
/* to NULL and old->next[dir_to_beg] set to NULL as well. */

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

phydbl PHYREX_Effective_Density(t_tree *tree)
{
  return(PHYREX_Neighborhood_Size(tree)/(4.*PI*PHYREX_Update_Sigsq(tree)));
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

t_ldsk *PHYREX_Generate_Path(t_ldsk *young, t_ldsk *old, phydbl n_evt, phydbl sd, t_tree *tree)
{
  int i,j,swap,err;
  phydbl *time,dum,mode;
  t_ldsk *path,**ldsk_a;
  t_dsk *disk;
  phydbl X,Y,Xp,Yp;
  phydbl slope,inter;
  
  path        = NULL;
  
  if(n_evt == 0) return(NULL); // path is set to NULL

  time   = (phydbl *)mCalloc(n_evt,sizeof(phydbl));
  ldsk_a = (t_ldsk **)mCalloc(n_evt,sizeof(t_ldsk *));

  
  for(i=0;i<n_evt;i++) time[i] =  Uni()*(young->disk->time - old->disk->time) + old->disk->time;
  
  /* Bubble sort time in decreasing order */
  do
    {
      swap = NO;
      for(i=0;i<n_evt-1;i++) 
        {
          if(time[i+1] > time[i])
            {
              swap = YES;
              dum       = time[i+1];
              time[i+1] = time[i];
              time[i]   = dum;
            }
        }
    }
  while(swap == YES);
  
  for(i=0;i<n_evt;i++)
    {
      ldsk_a[i] = PHYREX_Make_Lindisk_Node(tree->mmod->n_dim);
      disk = PHYREX_Make_Disk_Event(tree->mmod->n_dim,tree->n_otu);
      PHYREX_Init_Lindisk_Node(ldsk_a[i],disk,tree->mmod->n_dim);
      PHYREX_Make_Lindisk_Next(ldsk_a[i]);
      PHYREX_Init_Disk_Event(disk,tree->mmod->n_dim,tree->mmod);
      disk->ldsk = ldsk_a[i];
      disk->time = time[i];      
    }

  for(i=0;i<n_evt-1;i++) ldsk_a[i]->prev = ldsk_a[i+1];
  ldsk_a[i]->prev = NULL;

  for(i=1;i<n_evt;i++) ldsk_a[i]->next[0] = ldsk_a[i-1];
  ldsk_a[0]->next[0] = NULL;

  path = ldsk_a[0];

  /* Instantiate path */
  for(i=0;i<tree->mmod->n_dim;i++)
    {
      X = old->coord->lonlat[i];
      Y = old->disk->time;

      for(j=0;j<n_evt;j++)
        {

          if(j == 0)
            {
              Xp = young->coord->lonlat[i];
              Yp = young->disk->time;
            }
          else
            {
              assert(ldsk_a[j]->n_next == 1);
              Xp = ldsk_a[j]->next[0]->coord->lonlat[i];
              Yp = ldsk_a[j]->next[0]->disk->time;
            }
          
          slope = (Yp-Y)/(Xp-X);
          inter = Y - slope*X;

          if(j == 0)
            {
              mode = (young->disk->time - inter)/slope;
            }
          else
            {
              mode = (ldsk_a[j]->next[0]->disk->time - inter)/slope;
            }

          ldsk_a[j]->coord->lonlat[i] = Rnorm_Trunc(mode,
                                                    sd,
                                                    tree->mmod->lim_do->lonlat[i],
                                                    tree->mmod->lim_up->lonlat[i],&err);
          
          
          /* ldsk_a[j]->disk->centr->lonlat[i] = Rnorm_Trunc(mode, */
          /*                                                 sd, */
          /*                                                 0.0, */
          /*                                                 tree->mmod->lim->lonlat[i],&err); */

          
          ldsk_a[j]->disk->centr->lonlat[i] = Rnorm_Trunc(ldsk_a[j]->coord->lonlat[i],
                                                          sd,
                                                          tree->mmod->lim_do->lonlat[i],
                                                          tree->mmod->lim_up->lonlat[i],&err);
        }
    }
  
  Free(ldsk_a);
  Free(time);

  return(path);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl PHYREX_Path_Logdensity(t_ldsk *young, t_ldsk *old, phydbl sd, t_tree *tree)
{
  int i,j,err;
  t_ldsk *ldsk;
  phydbl lnDens,mode;
  phydbl X,Y,Xp,Yp;
  phydbl slope,inter;
  int dir_to_young;
  
  lnDens = 0.0;
  mode   = 0.0;
    
  for(i=0;i<tree->mmod->n_dim;i++)
    {     

      j    = 0;
      ldsk = young->prev;

      X = old->coord->lonlat[i];
      Y = old->disk->time;

      // Density up to the last jump (to old ldsk) which should
      // not be accounted for (it is not part of the simulation
      // when randomly generating a path using PHYREX_Generate_Path
      while(ldsk != old)
        {

          if(ldsk->n_next > 1)
            dir_to_young = PHYREX_Get_Next_Direction(young,ldsk);
          else
            dir_to_young = 0;
          
          Xp = ldsk->next[dir_to_young]->coord->lonlat[i];
          Yp = ldsk->next[dir_to_young]->disk->time;
          
          slope = (Yp-Y)/(Xp-X);
          inter = Y - slope*X;

          assert(ldsk != NULL);
          
          /* mode = (old->coord->lonlat[i] - young->coord->lonlat[i])/(n_evt+1.-j) + young->coord->lonlat[i]; */

          mode = (ldsk->disk->time - inter)/slope;
          
          lnDens += Log_Dnorm_Trunc(ldsk->coord->lonlat[i],
                                    mode,
                                    sd,
                                    tree->mmod->lim_do->lonlat[i],
                                    tree->mmod->lim_up->lonlat[i],&err);

          /* lnDens += Log_Dnorm_Trunc(ldsk->disk->centr->lonlat[i], */
          /*                           mode, */
          /*                           sd, */
          /*                           0.0, */
          /*                           tree->mmod->lim->lonlat[i],&err); */

          lnDens += Log_Dnorm_Trunc(ldsk->disk->centr->lonlat[i],
                                    ldsk->coord->lonlat[i],
                                    sd,
                                    tree->mmod->lim_do->lonlat[i],
                                    tree->mmod->lim_up->lonlat[i],&err);
          /* lnDens += log(1./tree->mmod->lim->lonlat[i]); */

          ldsk = ldsk->prev;
          j++;
        }
    }
  
  return(lnDens);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void PHYREX_Sample_Path(t_ldsk *young, t_ldsk *old, phydbl sd, phydbl *global_hr, t_tree *tree)
{
  phydbl new_ldsk_pos,cur_ldsk_pos;
  phydbl new_centr_pos,cur_centr_pos;
  phydbl new_glnL,cur_glnL;
  phydbl u,alpha,ratio,hr;
  int n_mcmc_iter;
  int i,err;
  t_ldsk *ldsk;
  int accept,reject;
  int path_len;
  

  path_len = PHYREX_Path_Len(young,old)-2;
  if(path_len == 0) return;

  assert(young != old);
  new_glnL = cur_glnL = tree->mmod->c_lnL;
  
  n_mcmc_iter = 0;
  accept = reject = 0;
  do
    {
          
      ldsk = young->prev;
      do
        {
          assert(ldsk->prev);
          
          PHYREX_Store_Geo_Coord(ldsk->coord);
          PHYREX_Store_Geo_Coord(ldsk->disk->centr);
          
          hr = 0.0;
          ratio = 0.0;
          
          /* cur_glnL = PHYREX_Lk_Range(young->disk,old->disk,tree); */
          cur_glnL = PHYREX_Lk_Range(ldsk->disk,ldsk->prev->disk,tree);
          new_glnL = cur_glnL;
          
          for(i=0;i<tree->mmod->n_dim;i++)    
            {
              cur_ldsk_pos  = ldsk->coord->lonlat[i];
              cur_centr_pos = ldsk->disk->centr->lonlat[i];
              
              /* new_centr_pos = Rnorm_Trunc(cur_centr_pos, */
              /*                             0.5*sd, */
              /*                             0.0, */
              /*                             tree->mmod->lim->lonlat[i],&err); */
              
              new_centr_pos = Uni() * (tree->mmod->lim_up->lonlat[i]-tree->mmod->lim_do->lonlat[i])+tree->mmod->lim_do->lonlat[i];
              
              new_ldsk_pos = Rnorm_Trunc(new_centr_pos,
                                         sd,
                                         tree->mmod->lim_do->lonlat[i],
                                         tree->mmod->lim_up->lonlat[i],&err);
              
              
              /* PhyML_Printf("\n. cur_centr: %12f new_centr: %12f",cur_centr_pos,new_centr_pos); */
              /* PhyML_Printf(" d(new|cur)=%12f d(cur|new): %12f", */
              /*              Log_Dnorm_Trunc(new_centr_pos, */
              /*                              cur_centr_pos, */
              /*                              0.5*sd, */
              /*                              0.0, */
              /*                              tree->mmod->lim->lonlat[i],&err), */
              /*              Log_Dnorm_Trunc(cur_centr_pos, */
              /*                              new_centr_pos, */
              /*                              0.5*sd, */
              /*                              0.0, */
              /*                              tree->mmod->lim->lonlat[i],&err)); */
              
              
              /* PhyML_Printf(" cur_ldsk: %12f new_ldsk: %12f",cur_ldsk_pos,new_ldsk_pos); */
              /* PhyML_Printf(" d(new|cur)=%12f d(cur|new): %12f", */
              /*              Log_Dnorm_Trunc(new_ldsk_pos, */
              /*                              new_centr_pos, */
              /*                              0.5*sd, */
              /*                              0.0, */
              /*                              tree->mmod->lim->lonlat[i],&err), */
              /*              Log_Dnorm_Trunc(cur_ldsk_pos, */
              /*                              cur_centr_pos, */
              /*                              0.5*sd, */
              /*                              0.0, */
              /*                              tree->mmod->lim->lonlat[i],&err)); */
              
              
              hr -= Log_Dnorm_Trunc(new_ldsk_pos,
                                    new_centr_pos,
                                    sd,
                                    tree->mmod->lim_do->lonlat[i],
                                    tree->mmod->lim_up->lonlat[i],&err);
              
              /* hr -= Log_Dnorm_Trunc(new_centr_pos, */
              /*                       cur_centr_pos, */
              /*                       0.5*sd, */
              /*                       0.0, */
              /*                       tree->mmod->lim->lonlat[i],&err); */
              
              
              hr += Log_Dnorm_Trunc(cur_ldsk_pos,
                                    cur_centr_pos,
                                    sd,
                                    tree->mmod->lim_do->lonlat[i],
                                    tree->mmod->lim_up->lonlat[i],&err);
              
              /* hr += Log_Dnorm_Trunc(cur_centr_pos, */
              /*                       new_centr_pos, */
              /*                       0.5*sd, */
              /*                       0.0, */
              /*                       tree->mmod->lim->lonlat[i],&err); */
              
              
              ldsk->coord->lonlat[i] = new_ldsk_pos;
              ldsk->disk->centr->lonlat[i]  = new_centr_pos;                            
            }
          
          /* new_glnL = PHYREX_Lk_Range(young->disk,old->disk,tree); */
          /* new_glnL = PHYREX_Lk_Core_Range(young->disk,old->disk,tree); */
          new_glnL = PHYREX_Lk_Range(ldsk->disk,ldsk->prev->disk,tree);
          /* new_glnL = PHYREX_Lk(tree); */
          
          
          /* PhyML_Printf("\n. k: %4d lnL: %f",k,cur_glnL); */
          
          ratio += (new_glnL - cur_glnL);
          ratio += hr;
          
          /* (*global_hr) -= ratio; */
                      
          ratio = exp(ratio);
          alpha = MIN(1.,ratio);
          
          u = Uni();
          
          if(u > alpha) /* Reject */
            {
              PHYREX_Restore_Geo_Coord(ldsk->coord);
              PHYREX_Restore_Geo_Coord(ldsk->disk->centr);
              reject++;
            }
          else
            {
              cur_glnL = new_glnL;
              accept++;
            }
          
          ldsk = ldsk->prev;
        }
      while(ldsk != old);
      
      n_mcmc_iter++;
    }
  while(n_mcmc_iter < 100 * path_len);
  /* while(n_mcmc_iter < 1); */
  /* PhyML_Printf("\n. ratio: %f %d %d",(phydbl)(accept)/(phydbl)(accept+reject),accept,reject); */
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

  /* printf("\n. len1: %d len2: %d",len1,len2); fflush(NULL); */
  /* printf("\n. %f %f %d %d [%f %f]", */
  /*        list1[1]->disk->time, */
  /*        list2[1]->disk->time, */
  /*        len1,len2, */
  /*        n1->disk->time, */
  /*        n2->disk->time); */
         
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
// Coalescence rate with time expressed in calendar unit
phydbl PHYREX_Coalescence_Rate(t_tree *tree)
{
  phydbl mu,theta,lbda,w,h;

  mu    = tree->mmod->mu;
  theta = tree->mmod->rad;
  lbda  = tree->mmod->lbda;
  w     = (tree->mmod->lim_up->lonlat[0]-tree->mmod->lim_do->lonlat[0]);
  h     = (tree->mmod->lim_up->lonlat[1]-tree->mmod->lim_do->lonlat[1]);

  return(4.*POW(PI,2)*POW(theta,4)*POW(mu,2)*lbda / (POW(w,2)*POW(h,2)));
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/


phydbl Prob_Two_Lineages_Coal_One_Event(phydbl w, phydbl h, phydbl mu, phydbl rad)
{
  phydbl cx,cy;
  phydbl l1x,l1y;
  phydbl l2x,l2y;
  phydbl prob_hit;
  int n_hit,n_trials,trial;

  n_trials = 10000000;
  trial = 0;
  n_hit = 0;
  do
    {
      cx = Uni()*w;
      cy = Uni()*h;
      
      l1x = Uni()*w;
      l1y = Uni()*h;

      l2x = Uni()*w;
      l2y = Uni()*h;

      prob_hit = log(mu);
      prob_hit += -POW(l1x - cx,2)/(2.*pow(rad,2));
      prob_hit += -POW(l1y - cy,2)/(2.*pow(rad,2));

      prob_hit += log(mu);
      prob_hit += -POW(l2x - cx,2)/(2.*pow(rad,2));
      prob_hit += -POW(l2y - cy,2)/(2.*pow(rad,2));

      prob_hit = exp(prob_hit);
      
      if(!(Uni() > prob_hit)) n_hit++;

      trial++;
    }
  while(trial != n_trials);
  
  return((phydbl)n_hit/n_trials);
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


  /* disk = tree->young_disk; */
  /* do */
  /*   { */
  /*     PhyML_Printf("\n. %s %d",disk->id,disk->age_fixed); */
  /*     if(disk->age_fixed == YES) */
  /*       { */
  /*         for(i=0;i<disk->n_ldsk_a;++i) */
  /*           { */
  /*             PhyML_Printf("\n. %s %f %s",disk->id,disk->time,disk->ldsk_a[i] ? disk->ldsk_a[i]->coord->id : "?"); */
  /*           } */
  /*       } */
  /*     disk = disk->prev; */
  /*   } */
  /* while(disk); */
  
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

int PHYREX_Scale_All(phydbl scale, t_dsk *start_disk, t_tree *tree)
{
  t_dsk *disk;
  t_dsk **sorted_disk;
  int n_disk,n_disk_scaled,sorted,i;

  n_disk        = 0;
  n_disk_scaled = 0;

  disk = start_disk->prev;
  assert(disk);
  
  do
    {      
      if(disk->age_fixed == NO)
        {
          disk->time = disk->time * scale + start_disk->time * (1.-scale);
          n_disk_scaled++;
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

  return(n_disk_scaled);
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
      dt += fabs(disk->next->time - disk->time);

      if(disk->age_fixed == NO)
        {
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
phydbl PHYREX_Realized_Sigsq(t_tree *tree)
{
  t_dsk *disk,*root_disk;
  int i;
  phydbl sum;
  phydbl t_root, t_tip;
  
  disk = tree->young_disk;
  while(disk->prev) disk = disk->prev;
  root_disk = disk;
  
  t_root = root_disk->time;
  t_tip = tree->young_disk->time;
  assert(t_tip - t_root > 0.0);
  
  disk = tree->young_disk;
  sum = 0.0;
  for(i=0;i<disk->n_ldsk_a;++i)
    sum += pow(Euclidean_Dist(root_disk->ldsk->coord,disk->ldsk_a[i]->coord) / (t_tip - t_root),2);
  
  return(sum/(tree->mmod->n_dim * disk->n_ldsk_a));
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
// Mean Haversine distance (in km) per year between root and tips on the most recent disk 
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
/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
