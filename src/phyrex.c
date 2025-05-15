/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

/* Routines that implement routines that deal with PhyREX data structure
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
 
 res = NULL;
 
  mixt_tree = XML_Process_Base(xml_filename);
  assert(mixt_tree);

  mixt_tree->rates = RATES_Make_Rate_Struct(mixt_tree->n_otu);
  RATES_Init_Rate_Struct(mixt_tree->rates,NULL,mixt_tree->n_otu);
  
  mixt_tree->times = TIMES_Make_Time_Struct(mixt_tree->n_otu);
  TIMES_Init_Time_Struct(mixt_tree->times,NULL,mixt_tree->n_otu);

  mixt_tree->mmod = PHYREX_Make_Migrep_Model(mixt_tree->n_otu,2);
  mixt_tree->mmod->n_dim = 2;
  PHYREX_Set_Default_Migrep_Mod(mixt_tree->n_otu,mixt_tree->mmod);
  
  Make_Contrasts(mixt_tree);
  Make_Contmod(mixt_tree);
  Set_Defaults_Contmod(mixt_tree);

  
  tree = mixt_tree;
  do
    {
      // All rate stuctures point to the same object
      tree->rates   = mixt_tree->rates;
      tree->times   = mixt_tree->times;
      tree->contmod = mixt_tree->contmod;
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


  dum_string = XML_Get_Attribute_Value(xnd,"mcmc.output.times");
  if(dum_string != NULL)
    {
      int select = XML_Validate_Attr_Int(dum_string,6,
                                         "true","yes","y",
                                         "false","no","n");
      if(select < 3) mixt_tree->io->mcmc_output_times = YES;
      else mixt_tree->io->mcmc_output_times = NO;
    }

  dum_string = XML_Get_Attribute_Value(xnd,"mcmc.output.trees");
  if(dum_string != NULL)
    {
      int select = XML_Validate_Attr_Int(dum_string,6,
                                         "true","yes","y",
                                         "false","no","n");
      if(select < 3) mixt_tree->io->mcmc_output_trees = YES;
      else mixt_tree->io->mcmc_output_trees = NO;
    }

  dum_string = XML_Get_Attribute_Value(xnd,"mcmc.verbose");
  if(dum_string != NULL)
    {
      int select = XML_Validate_Attr_Int(dum_string,6,
                                         "true","yes","y",
                                         "false","no","n");
      if(select < 3) mixt_tree->io->mcmc->out_verbose = 1;
      else mixt_tree->io->mcmc->out_verbose = 0;
    }




  
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
      mixt_tree->rates->model_id = GUINDON;
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
              mixt_tree->rates->model_id = GUINDON;
              mixt_tree->mod->gamma_mgf_bl = YES;
              strcpy(mixt_tree->rates->model_name,"geometric Brownian"); 
            }
          else if(!strcmp(model_name,"geometric"))
            {
              mixt_tree->rates->model_id = GUINDON;
              mixt_tree->mod->gamma_mgf_bl = YES;
              strcpy(mixt_tree->rates->model_name,"geometric Brownian"); 
            }
          else if(!strcmp(model_name,"brownian"))
            {
              mixt_tree->rates->model_id = GUINDON;
              mixt_tree->mod->gamma_mgf_bl = YES;
              strcpy(mixt_tree->rates->model_name,"geometric Brownian"); 
            }
          else if(!strcmp(model_name,"integrated"))
            {
              mixt_tree->rates->model_id = GUINDON;
              mixt_tree->mod->gamma_mgf_bl = YES;
              strcpy(mixt_tree->rates->model_name,"geometric Brownian"); 
            }
          else if(!strcmp(model_name,"geo"))
            {
              mixt_tree->rates->model_id = GUINDON;
              mixt_tree->mod->gamma_mgf_bl = YES;
              strcpy(mixt_tree->rates->model_name,"geometric Brownian"); 
            }
          else if(!strcmp(model_name,"lognormal"))
            {
              mixt_tree->rates->model_id = LOGNORMAL;
              mixt_tree->mod->gamma_mgf_bl = NO;
              strcpy(mixt_tree->rates->model_name,"lognormal (uncorrelated)"); 
            }
          else if(!strcmp(model_name,"uncorrelated"))
            {
              mixt_tree->rates->model_id = LOGNORMAL;
              mixt_tree->mod->gamma_mgf_bl = NO;
              strcpy(mixt_tree->rates->model_name,"lognormal (uncorrelated)"); 
            }
          else if(!strcmp(model_name,"uncorr"))
            {
              mixt_tree->rates->model_id = LOGNORMAL;
              mixt_tree->mod->gamma_mgf_bl = NO;
              strcpy(mixt_tree->rates->model_name,"lognormal (uncorrelated)"); 
            }
          else if(!strcmp(model_name,"normal"))
            {
              mixt_tree->rates->model_id = LOGNORMAL;
              mixt_tree->mod->gamma_mgf_bl = NO;
              strcpy(mixt_tree->rates->model_name,"lognormal (uncorrelated)"); 
            }
          else if(!strcmp(model_name,"gamma"))
            {
              mixt_tree->rates->model_id = GAMMA;
              mixt_tree->mod->gamma_mgf_bl = NO;
              strcpy(mixt_tree->rates->model_name,"gamma (uncorrelated)"); 
            }
          else if(!strcmp(model_name,"strictclock"))
            {
              mixt_tree->rates->model_id = STRICTCLOCK;
              mixt_tree->mod->gamma_mgf_bl = NO;
              strcpy(mixt_tree->rates->model_name,"strict clock"); 
            }
          else if(!strcmp(model_name,"clock"))
            {
              mixt_tree->rates->model_id = STRICTCLOCK;
              mixt_tree->mod->gamma_mgf_bl = NO;
              strcpy(mixt_tree->rates->model_name,"strict clock"); 
            }
          else if(!strcmp(model_name,"thorne"))
            {
              mixt_tree->rates->model_id = THORNE;
              mixt_tree->mod->gamma_mgf_bl = NO;
              strcpy(mixt_tree->rates->model_name,"autocorrelated"); 
            }
          else if(!strcmp(model_name,"autocorrelated"))
            {
              mixt_tree->rates->model_id = THORNE;
              mixt_tree->mod->gamma_mgf_bl = NO;
              strcpy(mixt_tree->rates->model_name,"autocorrelated"); 
            }
          else if(!strcmp(model_name,"autocorr"))
            {
              mixt_tree->rates->model_id = THORNE;
              mixt_tree->mod->gamma_mgf_bl = NO;
              strcpy(mixt_tree->rates->model_name,"autocorrelated"); 
            }
          else
            {
              assert(FALSE);
            }
        }

      char *autocor_prior;
      autocor_prior = XML_Get_Attribute_Value(xnd,"autocor.prior.rate");
      if(autocor_prior != NULL)
        {
          mixt_tree->rates->autocor_rate_prior = String_To_Dbl(autocor_prior);
          if(mixt_tree->rates->autocor_rate_prior < 0.0)
            {
              PhyML_Printf("\n. Autocorrelation prior needs to be stricly positive.");
              assert(false);
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

      mixt_tree->mmod->use_locations = YES;

      if(modname != NULL)
        {
          if(!strcmp(modname,"slfv")) mixt_tree->mmod->model_id = SLFV_GAUSSIAN;
          else if(!strcmp(modname,"rw")) mixt_tree->mmod->model_id = RW;
          else if(!strcmp(modname,"rrw+gamma")) mixt_tree->mmod->model_id = RRW_GAMMA;
          else if(!strcmp(modname,"rrw+lognormal")) mixt_tree->mmod->model_id = RRW_LOGNORMAL;
          else if(!strcmp(modname,"ibm")) mixt_tree->mmod->model_id = IBM;
          else if(!strcmp(modname,"ribm")) mixt_tree->mmod->model_id = RIBM;
          else if(!strcmp(modname,"iwn")) mixt_tree->mmod->model_id = IWNc;
          else if(!strcmp(modname,"riwn")) mixt_tree->mmod->model_id = RIWNc;
          else if(!strcmp(modname,"iwnu")) mixt_tree->mmod->model_id = IWNu;
          else if(!strcmp(modname,"riwnu")) mixt_tree->mmod->model_id = RIWNu;
          else if(!strcmp(modname,"iou")) mixt_tree->mmod->model_id = IOU;
          else
            {
              PhyML_Printf("\n. Unknown spatial model name '%s'. Aborting. ",modname);
              assert(FALSE);
            }
        }

      char *dist_type;
      dist_type = XML_Get_Attribute_Value(xnd,"distance.type");
      XML_Set_Attribute_Value(xnd,"distance.type",dist_type);
      
      if(dist_type != NULL)
        {
          if(!strcmp(dist_type,"great circle"))   mixt_tree->mmod->dist_type = HAVERSINE;
          if(!strcmp(dist_type,"greatcircle"))    mixt_tree->mmod->dist_type = HAVERSINE;
          else if(!strcmp(dist_type,"haversine")) mixt_tree->mmod->dist_type = HAVERSINE;
          else if(!strcmp(dist_type,"manhattan")) mixt_tree->mmod->dist_type = MANHATTAN;
          else if(!strcmp(dist_type,"euclidean")) mixt_tree->mmod->dist_type = EUCLIDEAN;
          else
            {
              PhyML_Printf("\n. Unknown distance type '%s'. Aborting. ",dist_type);
              assert(FALSE);
            }
        }

      char *sampling;
      sampling = XML_Get_Attribute_Value(xnd,"sampling");
      XML_Set_Attribute_Value(xnd,"sampling",sampling);

      if(sampling != NULL)
        {
          if(!strcmp(sampling,"detection")) mixt_tree->mmod->sampling_scheme = SPATIAL_SAMPLING_DETECTION;
          else if(!strcmp(sampling,"survey")) mixt_tree->mmod->sampling_scheme = SPATIAL_SAMPLING_SURVEY;
          else
            {
              PhyML_Printf("\n. Unknown sampling scheme '%s'. Aborting. ",sampling);
              assert(FALSE);
            }
        }
      
      char *integrateAncestralLocations;
      integrateAncestralLocations = XML_Get_Attribute_Value(xnd,"integrateAncestralLocations");

      if(integrateAncestralLocations != NULL)
        {
          int select = XML_Validate_Attr_Int(integrateAncestralLocations,6,
                                             "true","yes","y",
                                             "false","no","n");
          if(select < 3) mixt_tree->mmod->integrateAncestralLocations = YES;
          else mixt_tree->mmod->integrateAncestralLocations = NO;
        }
      
      char *rrw_prior_sd;
      rrw_prior_sd = XML_Get_Attribute_Value(xnd,"rrw.prior.sd");

      if(rrw_prior_sd != NULL)
        {
          char *val;
          val = XML_Get_Attribute_Value(xnd,"rrw.prior.sd");
          if(val != NULL) mixt_tree->mmod->rrw_prior_sd = String_To_Dbl(val);
        }

      char *rw_prior_sd;
      rw_prior_sd = XML_Get_Attribute_Value(xnd,"rw.prior.sd");

      if(rw_prior_sd != NULL)
        {
          mixt_tree->mmod->rw_prior_sd = String_To_Dbl(rw_prior_sd);
        }

      char *rw_prior_mean;
      rw_prior_mean = XML_Get_Attribute_Value(xnd,"rw.prior.mean");

      if(rw_prior_mean != NULL)
        {
          mixt_tree->mmod->rw_prior_mean = String_To_Dbl(rw_prior_mean);
        }

      char *prior_distrib;
      prior_distrib = XML_Get_Attribute_Value(xnd,"rw.prior.distrib");
      if(prior_distrib != NULL)
        {
          if(!strcmp(prior_distrib,"exponential")) mixt_tree->mmod->rw_prior_distrib = EXPONENTIAL_PRIOR;
          else if(!strcmp(prior_distrib,"normal")) mixt_tree->mmod->rw_prior_distrib = NORMAL_PRIOR;
          else if(!strcmp(prior_distrib,"flat")) mixt_tree->mmod->rw_prior_distrib = FLAT_PRIOR;
          else
            {
              PhyML_Printf("\n. Unknown prior distribution '%s'. Aborting. ",prior_distrib);
              assert(FALSE);
            }
        }

      char *rw_root_mean;
      char *rw_root_var;

      rw_root_mean = XML_Get_Attribute_Value(xnd,"rw.root.mean.location");

      if(rw_root_mean != NULL)
        {
          mixt_tree->mmod->rw_root_mean[LOCATION] = String_To_Dbl(rw_root_mean);
        }

      rw_root_var = XML_Get_Attribute_Value(xnd,"rw.root.var.location");

      if(rw_root_var != NULL)
        {
          mixt_tree->mmod->rw_root_var[LOCATION] = String_To_Dbl(rw_root_var);
        }

      rw_root_mean = XML_Get_Attribute_Value(xnd,"rw.root.mean.velocity");

      if(rw_root_mean != NULL)
        {
          mixt_tree->mmod->rw_root_mean[VELOCITY] = String_To_Dbl(rw_root_mean);
        }

      rw_root_var = XML_Get_Attribute_Value(xnd,"rw.root.var.velocity");

      if(rw_root_var != NULL)
        {
          mixt_tree->mmod->rw_root_var[VELOCITY] = String_To_Dbl(rw_root_var);
        }



      char *observational_model;
      observational_model = XML_Get_Attribute_Value(xnd,"observational.model");

      if(observational_model != NULL)
        {
          int select = XML_Validate_Attr_Int(observational_model,6,
                                             "true","yes","y",
                                             "false","no","n");
          if(select < 3) mixt_tree->contmod->obs_model = YES;
          else mixt_tree->contmod->obs_model = NO;
        }

      char *observational_var;
      observational_var = XML_Get_Attribute_Value(xnd,"observational.var");

      if(observational_var != NULL)
        {
          phydbl v = String_To_Dbl(observational_var);
          mixt_tree->contmod->obs_var[0] = v;
          mixt_tree->contmod->obs_var[1] = v;
        }

      char *observational_est;
      observational_est = XML_Get_Attribute_Value(xnd,"observational.var.est");

      if(observational_est != NULL)
        {
          int select = XML_Validate_Attr_Int(observational_est,6,
                                             "true","yes","y",
                                             "false","no","n");
          if(select < 3) mixt_tree->contmod->obs_model_est = YES;
          else mixt_tree->contmod->obs_model_est = NO;
        }

    }
  else
    {
      mixt_tree->mmod->use_locations = NO;
      mixt_tree->mmod->sampling_scheme = SPATIAL_SAMPLING_DETECTION;
      mixt_tree->mmod->model_id = RRW_GAMMA;
    }


  // Looking for XML node with tree generating model info
  xnd = XML_Search_Node_Name("treegenerating",YES,xroot);

  if(xnd != NULL)
    {
      char *treemodel;
      treemodel = XML_Get_Attribute_Value(xnd,"model");

      if(!strcmp(treemodel,"none")) mixt_tree->eval_tlnL = NO;
      else mixt_tree->eval_tlnL = YES;
      
      char *prior_mean;
      prior_mean = XML_Get_Attribute_Value(xnd,"neff.prior.mean");
      if(prior_mean != NULL) mixt_tree->times->neff_prior_mean = String_To_Dbl(prior_mean);
      
      char *prior_var;
      prior_var = XML_Get_Attribute_Value(xnd,"neff.prior.var");
      if(prior_var != NULL) mixt_tree->times->neff_prior_var = String_To_Dbl(prior_var);

      char *prior_distrib;
      prior_distrib = XML_Get_Attribute_Value(xnd,"neff.prior.distrib");
      if(prior_distrib != NULL)
        {
          if(!strcmp(prior_distrib,"exponential")) mixt_tree->times->neff_prior_distrib = EXPONENTIAL_PRIOR;
          else if(!strcmp(prior_distrib,"normal")) mixt_tree->times->neff_prior_distrib = NORMAL_PRIOR;
          else if(!strcmp(prior_distrib,"flat")) mixt_tree->times->neff_prior_distrib = FLAT_PRIOR;
          else
            {
              PhyML_Printf("\n. Unknown prior distribution '%s'. Aborting. ",prior_distrib);
              assert(FALSE);
            }
        }

      
      char *expgrowth;
      expgrowth = XML_Get_Attribute_Value(xnd,"expgrowth");
      if(expgrowth != NULL)
        {
          int select = XML_Validate_Attr_Int(expgrowth,6,
                                             "true","yes","y",
                                             "false","no","n");
          if(select < 3) mixt_tree->times->coalescent_model_id = EXPCOALESCENT;
          else mixt_tree->times->coalescent_model_id = STRICTCOALESCENT;
        }

      char *powgrowth;
      powgrowth = XML_Get_Attribute_Value(xnd,"powgrowth");
      if(powgrowth != NULL)
        {
          int select = XML_Validate_Attr_Int(powgrowth,6,
                                             "true","yes","y",
                                             "false","no","n");
          if(select < 3) mixt_tree->times->coalescent_model_id = POWLAW;
          else mixt_tree->times->coalescent_model_id = STRICTCOALESCENT;
        }

      char *fix_node_ages;
      fix_node_ages = XML_Get_Attribute_Value(xnd,"fix.node.ages");
      if(fix_node_ages != NULL)
        {
          int select = XML_Validate_Attr_Int(fix_node_ages,6,
                                             "true","yes","y",
                                             "false","no","n");
          if(select < 3) mixt_tree->mod->s_opt->opt_node_ages = NO;
          else mixt_tree->mod->s_opt->opt_node_ages = YES;
        }

      char *fix_neff;
      fix_neff = XML_Get_Attribute_Value(xnd,"fix.neff");
      if(fix_neff != NULL)
        {
          int select = XML_Validate_Attr_Int(fix_neff,6,
                                             "true","yes","y",
                                             "false","no","n");
          if(select < 3) mixt_tree->mod->s_opt->opt_neff = NO;
          else mixt_tree->mod->s_opt->opt_neff = YES;
        }
      
    }

  // Cross validation
  xnd = XML_Search_Node_Name("phyrex",YES,xroot);
  short int do_cv = YES;
  if(xnd != NULL)
    {
      char *crossvald;
      crossvald = XML_Get_Attribute_Value(xnd,"cross.validation");
      if(crossvald != NULL)
        {
          int select = XML_Validate_Attr_Int(crossvald,6,
                                             "true","yes","y",
                                             "false","no","n");
          if(select < 3)  do_cv = YES;
          else do_cv = NO;
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

      char *prior_mean;
      prior_mean = XML_Get_Attribute_Value(xnd,"prior.mean");
      if(prior_mean != NULL)
        {
          mixt_tree->rates->clock_r_has_prior = YES;
          mixt_tree->rates->clock_r_prior_mean = String_To_Dbl(prior_mean);
          if(mixt_tree->rates->clock_r_prior_mean < 0.0 || mixt_tree->rates->clock_r_prior_mean > 1000.)
            {
              PhyML_Printf("\n. Prior mean for clock rate should be set to a value in [0,1000.]");
              assert(false);
            }
        }

      char *prior_var;
      prior_var = XML_Get_Attribute_Value(xnd,"prior.var");
      if(prior_var != NULL)
        {
          mixt_tree->rates->clock_r_has_prior = YES;
          mixt_tree->rates->clock_r_prior_var = String_To_Dbl(prior_var);
          if(mixt_tree->rates->clock_r_prior_var < 0.0 || mixt_tree->rates->clock_r_prior_var > 1000.)
            {
              PhyML_Printf("\n. Prior variance for clock rate should be set to a value in [0,1000.]");
              assert(false);
            }
        }
      
      if((prior_var && !prior_mean) ||
         (!prior_var && prior_mean))
        {
          PhyML_Printf("\n. Could not read the prior mean or variance for the clock rate in your XML file.");
          PhyML_Printf("\n. Please use ``prior.mean='value', prior.var='value'.");
          assert(false);
        }

      char *init_value;
      init_value = XML_Get_Attribute_Value(xnd,"init.value");
      if(init_value != NULL)
        {
          mixt_tree->rates->init_clock_r = YES;
          mixt_tree->rates->clock_r = String_To_Dbl(init_value);
          if(mixt_tree->rates->clock_r < 0.0)
            {
              PhyML_Printf("\n. Initial value for clock rate should be positive.");
              assert(false);
            }
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
  MIXT_Chain_Models(mixt_tree);
  Init_Model(mixt_tree->mod->io->cdata,mixt_tree->mod,mixt_tree->mod->io);
  Set_Model_Parameters(mixt_tree->mod);
  Print_Data_Structure(NO,stdout,mixt_tree);
  tree = MIXT_Starting_Tree(mixt_tree);
  if(mixt_tree->io->in_tree < 2) Add_Root(tree->a_edges[0],tree);
  Copy_Tree(tree,mixt_tree);
  Free_Tree(tree);
  if(mixt_tree->mod->io->edge_len_unit == CALENDAR) Convert_Lengths_From_Calendar_To_Substitutions(mixt_tree);
  Copy_Tree(mixt_tree,mixt_tree->next);
  Connect_CSeqs_To_Nodes(mixt_tree->mod->io->cdata,mixt_tree->mod->io,mixt_tree);
  Init_T_Beg(mixt_tree);  
  Make_Tree_For_Lk(mixt_tree);
  Make_Tree_For_Pars(mixt_tree);
  Make_Spr(mixt_tree);  

  
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
        TIMES_Bl_To_Times(1, mixt_tree);

        // TIMES_Bl_To_Times(0, mixt_tree);
        // for (int i = 0; i < mixt_tree->n_otu; ++i)
        // {

        //     PhyML_Printf("\n<clade id=\"clad%d\">", i + 1);
        //     PhyML_Printf("\n\t<taxon value=\"%s\"/>", mixt_tree->a_nodes[i]->name);
        //     PhyML_Printf("\n</clade>");
        //     PhyML_Printf("\n<calibration id=\"cal%d\">", i + 1);
        //     PhyML_Printf("\n\t<lower>%f</lower>", mixt_tree->times->nd_t[i]);
        //     PhyML_Printf("\n\t<upper>%f</upper>", mixt_tree->times->nd_t[i]);
        //     PhyML_Printf("\n\t<appliesto clade.id=\"clad%d\"/>", i + 1);
        //     PhyML_Printf("\n</calibration>");
        // }

        // Exit("\n");

        Update_Ancestors(mixt_tree->n_root, mixt_tree->n_root->v[2], mixt_tree->n_root->b[2], mixt_tree);
        Update_Ancestors(mixt_tree->n_root, mixt_tree->n_root->v[1], mixt_tree->n_root->b[1], mixt_tree);
    }
  else
    {
      TIMES_Randomize_Tree_With_Time_Constraints(mixt_tree->times->a_cal[0],mixt_tree);
    }
  
  MIXT_Propagate_Tree_Update(mixt_tree);

  if(RRW_Is_Rw(mixt_tree->mmod) == YES || VELOC_Is_Integrated_Velocity(mixt_tree->mmod) == YES)
    {
      mixt_tree->aux_tree = (t_tree **)mCalloc(3,sizeof(t_tree *));


      /* Auxilliary tree for updating dispersal parameter in RRW and RW models first (i=0) and pop. size after (i=1) */
      for(int i=0;i<3;++i)
        {
          t_tree *aux_tree;
          
          mixt_tree->aux_tree[i] = Make_Tree_From_Scratch(mixt_tree->n_otu,mixt_tree->data);
          aux_tree = mixt_tree->aux_tree[i];
          
          aux_tree->mod = mixt_tree->mod;
          aux_tree->io = mixt_tree->io;          
          
          aux_tree->mmod = PHYREX_Make_Migrep_Model(mixt_tree->n_otu,2);
          aux_tree->mmod->n_dim = 2;
          PHYREX_Set_Default_Migrep_Mod(mixt_tree->n_otu,aux_tree->mmod);
          aux_tree->mmod->model_id = mixt_tree->mmod->model_id; 
          aux_tree->mmod->use_locations = mixt_tree->mmod->use_locations;
          aux_tree->mmod->integrateAncestralLocations = mixt_tree->mmod->integrateAncestralLocations;
          
          Copy_Tree(mixt_tree,aux_tree);

          aux_tree->rates = RATES_Make_Rate_Struct(mixt_tree->n_otu);
          RATES_Init_Rate_Struct(aux_tree->rates,NULL,mixt_tree->n_otu);
          RATES_Copy_Rate_Struct(mixt_tree->rates,aux_tree->rates,mixt_tree->n_otu);
          aux_tree->rates->model_id = LOGNORMAL;
      
          aux_tree->times = TIMES_Make_Time_Struct(mixt_tree->n_otu);
          TIMES_Init_Time_Struct(aux_tree->times,NULL,mixt_tree->n_otu);
          TIMES_Copy_Time_Struct(mixt_tree->times,aux_tree->times,mixt_tree->n_otu);
          
          RATES_Duplicate_Calib_Struct(mixt_tree,aux_tree);
          MIXT_Chain_Cal(aux_tree);  
          DATE_Assign_Primary_Calibration(aux_tree);
          
          aux_tree->data = mixt_tree->data;
          Connect_CSeqs_To_Nodes(mixt_tree->data,mixt_tree->io,aux_tree);
          Share_Lk_Struct(mixt_tree->next,aux_tree);
          
          aux_tree->mcmc = MCMC_Make_MCMC_Struct();
          MCMC_Init_MCMC_Struct(NULL,NULL,aux_tree->mcmc);
          MCMC_Complete_MCMC(aux_tree->mcmc,aux_tree);
          aux_tree->mcmc->chain_len_burnin = mixt_tree->io->mcmc->chain_len_burnin;
          aux_tree->mcmc->run = 1;


          Make_Contrasts(aux_tree);
          Make_Contmod(aux_tree);

          aux_tree->contmod->obs_var = mixt_tree->contmod->obs_var;
          aux_tree->contmod->obs_model = mixt_tree->contmod->obs_model;
          aux_tree->contmod->obs_model_est = mixt_tree->contmod->obs_model_est;
        }
    }

  
  
  /* Create ldsks and connect tree tips to them */
  /* once tip dates have been set properly (in */
  /* TIMES_Randomize_Tree_With_Time_Constraints) */
  PHYREX_Make_And_Connect_Tip_Disks(mixt_tree);
  for(int i=0;i<3;++i) if(mixt_tree->aux_tree && mixt_tree->aux_tree[i]) PHYREX_Make_And_Connect_Tip_Disks(mixt_tree->aux_tree[i]);
  
  /* Read spatial coordinates */
  PHYREX_Read_Tip_Coordinates(mixt_tree);
  for(int i=0;i<3;++i) if(mixt_tree->aux_tree && mixt_tree->aux_tree[i]) PHYREX_Read_Tip_Coordinates(mixt_tree->aux_tree[i]);

  PHYREX_Init_Migrep_Mod(mixt_tree->mmod,2,
                         mixt_tree->mmod->lim_do->lonlat[0],
                         mixt_tree->mmod->lim_do->lonlat[1],
                         mixt_tree->mmod->lim_up->lonlat[0],
                         mixt_tree->mmod->lim_up->lonlat[1]);


  /* Initialize parameters of migrep model */
  mixt_tree->mmod->lbda = 0.1;
  mixt_tree->mmod->mu   = 1.0;
  
  if(mixt_tree->mmod->model_id == SLFV_GAUSSIAN || mixt_tree->mmod->model_id == SLFV_UNIFORM)
    mixt_tree->mmod->rad  = 0.05*((mixt_tree->mmod->lim_up->lonlat[0]-mixt_tree->mmod->lim_do->lonlat[0])+
                                  (mixt_tree->mmod->lim_up->lonlat[1]-mixt_tree->mmod->lim_do->lonlat[1]));
  else if(RRW_Is_Rw(mixt_tree->mmod) == YES || VELOC_Is_Integrated_Velocity(mixt_tree->mmod) == YES)
    mixt_tree->mmod->rad = mixt_tree->mmod->sigsq[0];
  
  mixt_tree->mmod->sigsq[0] = PHYREX_Update_Sigsq(mixt_tree);
  for(int i=1;i<mixt_tree->mmod->n_dim;++i) mixt_tree->mmod->sigsq[i] = mixt_tree->mmod->sigsq[0];

  for(int i=0;i<3;++i)
    {
      if(mixt_tree->aux_tree && mixt_tree->aux_tree[i])
        {
          mixt_tree->aux_tree[i]->mmod->lbda     = mixt_tree->mmod->lbda;
          mixt_tree->aux_tree[i]->mmod->mu       = mixt_tree->mmod->mu;
          mixt_tree->aux_tree[i]->mmod->rad      = mixt_tree->mmod->rad;
          mixt_tree->aux_tree[i]->mmod->sigsq[0] = mixt_tree->mmod->sigsq[0];
          for(int j=1;j<mixt_tree->mmod->n_dim;++j) mixt_tree->aux_tree[i]->mmod->sigsq[j] = mixt_tree->mmod->sigsq[0];
        }
    }
  
  
  /* Random genealogy or user-defined tree */
  switch(mixt_tree->io->in_tree)
    {
    case 0 : case 1 :
      {
        if(RRW_Is_Rw(mixt_tree->mmod) == YES || VELOC_Is_Integrated_Velocity(mixt_tree->mmod) == YES)
          PHYREX_Simulate_Backward_Core(mixt_tree->young_disk,YES,mixt_tree);
        else
          PHYREX_Simulate_Backward_Core(mixt_tree->young_disk,NO,mixt_tree);

        PHYREX_Ldsk_To_Tree(mixt_tree);
        break;
      }
    case 2:
      {
        PHYREX_Tree_To_Ldsk(mixt_tree);
        break;
      }
    }



  Update_Ancestors(mixt_tree->n_root,mixt_tree->n_root->v[2],mixt_tree->n_root->b[2],mixt_tree);
  Update_Ancestors(mixt_tree->n_root,mixt_tree->n_root->v[1],mixt_tree->n_root->b[1],mixt_tree);  

  MCMC_Randomize_Rate_Across_Sites(mixt_tree);
  MCMC_Randomize_Rates(mixt_tree);
  MCMC_Randomize_Clock_Rate(mixt_tree);
  MCMC_Randomize_Sigsq_Scale(mixt_tree);
  MCMC_Randomize_Sigsq(mixt_tree);
  MCMC_Randomize_Locations(mixt_tree);
  MCMC_Randomize_Veloc(mixt_tree);
  MCMC_Randomize_Contmod(mixt_tree);

  Set_Ignore_Root(YES,mixt_tree);
  Set_Bl_From_Rt(YES,mixt_tree);


  PHYREX_Oldest_Sampled_Disk(mixt_tree);
  for(int i=0;i<3;++i) if(mixt_tree->aux_tree && mixt_tree->aux_tree[i]) PHYREX_Oldest_Sampled_Disk(mixt_tree->aux_tree[i]);

  if(RRW_Is_Rw(mixt_tree->mmod) == YES || VELOC_Is_Integrated_Velocity(mixt_tree->mmod) == YES)
    {
      
      mixt_tree->times->augmented_coalescent = NO;
      for(int i=0;i<3;++i) mixt_tree->aux_tree[i]->times->augmented_coalescent = NO;

      PHYREX_Remove_All_Disks_Except_Coal_And_Tips(mixt_tree);
      
      for(int i=0;i<3;++i)
        {
          PHYREX_Simulate_Backward_Core(mixt_tree->aux_tree[i]->young_disk,YES,mixt_tree->aux_tree[i]);
          PHYREX_Ldsk_To_Tree(mixt_tree->aux_tree[i]);
          PHYREX_Remove_All_Disks_Except_Coal_And_Tips(mixt_tree->aux_tree[i]);
          Update_Ancestors(mixt_tree->aux_tree[i]->n_root,mixt_tree->aux_tree[i]->n_root->v[2],mixt_tree->aux_tree[i]->n_root->b[2],mixt_tree->aux_tree[i]);
          Update_Ancestors(mixt_tree->aux_tree[i]->n_root,mixt_tree->aux_tree[i]->n_root->v[1],mixt_tree->aux_tree[i]->n_root->b[1],mixt_tree->aux_tree[i]);  
        }
    }
          
  assert(PHYREX_Check_Struct(mixt_tree,YES));
  TIMES_Lk(mixt_tree);
  RATES_Lk(mixt_tree);
  LOCATION_Prior(mixt_tree);
  TIMES_Prior(mixt_tree);
  RATES_Prior(mixt_tree);
  Set_Update_Eigen(YES,mixt_tree->mod);
  RATES_Update_Edge_Lengths(mixt_tree);

  mixt_tree->contmod->both_sides[LOCATION] = NO;
  mixt_tree->contmod->both_sides[VELOCITY] = NO;
  LOCATION_Lk(NULL,mixt_tree);

  Lk(NULL,mixt_tree);
  Set_Update_Eigen(NO,mixt_tree->mod);
  char *s = LOCATION_Model_Id(mixt_tree->mmod);
  PhyML_Printf("\n. Spatial diffusion model: %s",s);
  Free(s);
  PhyML_Printf("\n. Ne prior mean: %f",mixt_tree->times->neff_prior_mean);
  PhyML_Printf("\n. Exponential growth of Ne: %s",(mixt_tree->times->coalescent_model_id == EXPCOALESCENT)?"yes":"no");
  PhyML_Printf("\n. Init lnL(seq|phylo): %f",mixt_tree->c_lnL);
  PhyML_Printf("\n. Init lnL(rates|phylo): %f",mixt_tree->rates->c_lnL);
  PhyML_Printf("\n. Init lnL(coord|phylo): %f",mixt_tree->mmod->c_lnL);
  PhyML_Printf("\n. Init lnL(phylo): %f",mixt_tree->times->c_lnL);
  PhyML_Printf("\n. Init lnP(rates|phylo): %f",mixt_tree->rates->c_lnP);
  PhyML_Printf("\n. Init lnP(coord|phylo): %f",mixt_tree->mmod->c_lnP);
  PhyML_Printf("\n. Init lnP(phylo): %f",mixt_tree->times->c_lnP);
  PhyML_Printf("\n. Init clock rate: %f",mixt_tree->rates->clock_r);
  PhyML_Printf("\n. Random seed: %d",mixt_tree->io->r_seed);

  if (do_cv == YES)
      MCMC_Crossvalidate_Locations(mixt_tree);
  else
  {
      res = PHYREX_MCMC(mixt_tree);
      Free(res);
  }

  // Cleaning up...
  PHYREX_Free_Ldsk_Struct(mixt_tree);
  for(int i=0;i<3;++i) if(mixt_tree->aux_tree && mixt_tree->aux_tree[i]) PHYREX_Free_Ldsk_Struct(mixt_tree->aux_tree[i]);
  RATES_Free_Rates(mixt_tree->rates);
  for(int i=0;i<3;++i) if(mixt_tree->aux_tree && mixt_tree->aux_tree[i]) RATES_Free_Rates(mixt_tree->aux_tree[i]->rates);


  TIMES_Free_Times(mixt_tree->times);
  for(int i=0;i<3;++i) if(mixt_tree->aux_tree && mixt_tree->aux_tree[i]) TIMES_Free_Times(mixt_tree->aux_tree[i]->times);
  /* MCMC_Free_MCMC(mixt_tree->mcmc); */
  for(int i=0;i<3;++i) if(mixt_tree->aux_tree && mixt_tree->aux_tree[i]) MCMC_Free_MCMC(mixt_tree->aux_tree[i]->mcmc);
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
  for(int i=0;i<3;++i) if(mixt_tree->aux_tree && mixt_tree->aux_tree[i]) Free_Tree(mixt_tree->aux_tree[i]);  
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
  char *datadir;
  
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
  datadir = argv[9];
  
  printf("\n. seed: %d",seed);
  srand(seed);
  
  tree = PHYREX_Simulate(n_otus,n_sites,width,height,lbda,rad,mu,seed,modid,datadir);
  /* tree = PHYREX_Simulate_Independent_Loci(n_otus,500,20.,20.,seed); */
  
  disk = tree->young_disk;
  for(i=0;i<disk->n_ldsk_a;i++) PHYREX_Free_Ldisk(disk->ldsk_a[i]);
  while(disk->prev)
    {
      disk = disk->prev;
      if(disk->next->ldsk != NULL) PHYREX_Free_Ldisk(disk->next->ldsk);
      PHYREX_Free_Disk(disk->next);
    }
  
  /* Root */
  PHYREX_Free_Ldisk(disk->ldsk);
  PHYREX_Free_Disk(disk);

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

  if(mmod->model_id == SLFV_UNIFORM)
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
  else if(mmod->model_id == SLFV_GAUSSIAN)
    {
      return(YES);
    }

  return(-1);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl PHYREX_Lk(t_tree *tree)
{
  phydbl lnP;
  
  switch(tree->mmod->model_id)
    {
    case SLFV_GAUSSIAN : case SLFV_UNIFORM :
      {
        lnP = SLFV_Lk_Gaussian(tree) + TIMES_Lk_SLFV(tree);
        break;
      }
    case RW :
      {
        lnP = RW_Lk(tree) + TIMES_Lk_Coalescent(tree);
        break;
      }
    case RRW_GAMMA : case RRW_LOGNORMAL : 
      {
        lnP = RRW_Lk(tree) + TIMES_Lk_Coalescent(tree);
        break;
      }
    case IBM : case RIBM : case IWNc : case IWNu : case RIWNc : case RIWNu : case IOU :  
      {
        lnP = VELOC_Lk(NULL,tree) + TIMES_Lk_Coalescent(tree);
        break;
      }
    default : assert(FALSE);
    }

  return(lnP);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl PHYREX_Lk_Core(t_dsk *disk, t_tree *tree)
{
  phydbl lnL;
  
  lnL = UNLIKELY;
  
  switch(tree->mmod->model_id)
    {
    case SLFV_GAUSSIAN : case SLFV_UNIFORM :
      {
        lnL = SLFV_Lk_Gaussian_Core(disk,tree);
        break;
      }
    case RW : case RRW_GAMMA : case RRW_LOGNORMAL :  
        {
          lnL = RRW_Lk_Core(disk,tree);
          break;
        }
    case IBM : case RIBM : case IWNc : case IWNu : case RIWNc : case RIWNu : case IOU : 
      {
        PhyML_Printf("\n. PHYREX_Lk_Core function not implemented for IBM not IWN model");
        Exit("\n");
        break;
      }
    default : assert(FALSE);
    }
  
  return lnL;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// Warning: the calculation below does not incoporate time information
// since things get messy when considering serial samples
phydbl PHYREX_Lk_Range(t_dsk *young, t_dsk *old, t_tree *tree)
{    
  switch(tree->mmod->model_id)
    {
    case SLFV_GAUSSIAN : case SLFV_UNIFORM :
      {
        return(SLFV_Lk_Gaussian_Range(young,old,tree) + TIMES_Lk_SLFV_Range(young,old,tree));
        break;
      }
    case RW :
      {
        return(RW_Lk_Range(young,old,tree) + TIMES_Lk_Coalescent_Range(young,old,tree));
        break;
      }
    case RRW_GAMMA : case RRW_LOGNORMAL :
      {
        return(RRW_Lk_Range(young,old,tree) + TIMES_Lk_Coalescent_Range(young,old,tree));
        break;
      }
    case IBM : case RIBM : case IWNc : case IWNu : case RIWNc : case RIWNu : case IOU : 
      {
        PhyML_Printf("\n. PHYREX_Lk_Core function not implemented for IBM nor IWN models");
        Exit("\n");
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
  MCMC_Run(tree);
  return(NULL);
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

  if(prev != NULL) prev->next = next;
  if(next != NULL) next->prev = prev;

  disk->next = NULL;
  disk->prev = NULL;
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
  while(disk->prev != NULL && disk->prev->time > ins->time)
    {
      /* PhyML_Printf("\n disk->prev->time: %f ins->time: %f",disk->prev->time,ins->time); */
      disk = disk->prev;
    }
  PHYREX_Insert_Disk_At(ins,disk);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void PHYREX_Move_Disk_Updown(t_dsk *this, phydbl target_time, t_tree *tree)
{
  t_dsk *disk;

  if(target_time < this->time)
    {
      disk = this->prev;
      while(disk && target_time < disk->time)
        {
          if((disk->prev != NULL && disk->prev->time < target_time) || (disk->prev == NULL))
            {
              PHYREX_Remove_Disk(this);
              PHYREX_Insert_Disk_At(this,disk);
              break;
            }
          disk = disk->prev;
        }
    }
  else
    {
      disk = this->next;
      while(disk && target_time > disk->time)
        {
          if(disk->next != NULL && disk->next->time > target_time)
            {
              PHYREX_Remove_Disk(this);
              PHYREX_Insert_Disk_At(this,disk->next);
              break;
            }
          disk = disk->next;
          assert(disk);
        }
    }
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void PHYREX_Insert_Disk_At(t_dsk *ins, t_dsk *disk)
{
  ins->prev = disk->prev;
  ins->next = disk;
  disk->prev = ins;
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
      PhyML_Printf("\n. young (%s) @ time: %f; old (%s) @ time: %f delta: %G",
                   young->coord->id,young->disk->time,
                   old->coord->id,old->disk->time,
                   old->disk->time-young->disk->time);
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
  
  assert(young->time > old->time);
  
  disk = young;
  do
    {
      assert(disk);
      PHYREX_Update_Lindisk_List_Core(disk,tree);
      if(disk == old) break;      
      disk = disk->prev;
    }
  while(disk);  
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

  /* PhyML_Printf("\n. Time: %f",disk->time); */
  
  assert(disk->ldsk_a);
  
  // Set ldsk_a[i] to NULL if it does not point to a tip node
  for(i=0;i<tree->n_otu;++i)
    if(disk->ldsk_a[i] &&
       !(disk->ldsk_a[i]->nd != NULL &&
         disk->ldsk_a[i]->nd->tax == YES &&
         disk->age_fixed == YES &&
         disk->ldsk_a[i]->disk == disk))      
      {
        disk->ldsk_a[i] = NULL;
      }
  
  disk->n_ldsk_a = 0;
  for(i=0;i<tree->n_otu;++i)
    if(disk->ldsk_a[i] != NULL)
      {
        disk->n_ldsk_a++;
      }

  /* PhyML_Printf(" [%d]",disk->n_ldsk_a); */
  
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

  /* PhyML_Printf(" <%d>",disk->n_ldsk_a); */

  
  // A jump or coalescence has occurred on disk->next
  // --> add the lineage to disk->ldsk_a array
  if(disk->next->ldsk != NULL)
    {
      disk->ldsk_a[disk->n_ldsk_a] = disk->next->ldsk;
      disk->n_ldsk_a++;
    }

  /* PhyML_Printf(" !%d!",disk->n_ldsk_a); */

  /* PhyML_Printf("\n"); */
  /* for(i=0;i<disk->n_ldsk_a;++i) PhyML_Printf("\n. Disk %s [%12f] ldsk_a: %s (%3d)",disk->id,disk->time,disk->ldsk_a[i]->coord->id,disk->ldsk_a[i]->nd ? disk->ldsk_a[i]->nd->tax : -1); */
  /* if(disk->ldsk != NULL) PhyML_Printf("\n* Has %s on it (next: %s @ %12f %s)", */
  /*                                     disk->ldsk->coord->id, */
  /*                                     disk->ldsk->next[0]->coord->id, */
  /*                                     disk->ldsk->next[0]->disk->time, */
  /*                                     disk->ldsk->next[0]->prev->coord->id); */

  
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

  /* PHYREX_Update_Lindisk_List(tree); */

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

int PHYREX_Check_Struct(t_tree *tree, int exit)
{
  int i;
  t_dsk *disk;
  
  disk = tree->young_disk;

  do
    {
      if(disk->n_ldsk_a == 0)
        {
          if(exit == YES)
            {
              PhyML_Printf("\n. disk %s time %f has n_ldsk_a=0 (next has %d prev has %d)",
                       disk->id,
                       disk->time,
                       disk->next ? disk->next->n_ldsk_a : -1,
                       disk->prev ? disk->prev->n_ldsk_a : -1);
              assert(FALSE);
            }
          return 0;
        }

      disk = disk->prev;
    }
  while(disk);
  
  // Check times
  /* for(i=0;i<tree->n_otu;++i) */
  /*   { */
  /*     ldisk = tree->a_nodes[i]->ldsk; */
  /*     assert(ldisk); */
  /*     do */
  /*       { */
  /*         if(ldisk->prev->disk->time > ldisk->disk->time) */
  /*           { */
  /*             if(exit == YES) */
  /*               { */
  /*                 PhyML_Printf("\n. ldisk->id: %s ldisk->prev->id: %s ldsk->disk->time: %f  ldsk->prev->disk->time: %f diff: %g ldisk->prev->disk: %s ldisk->disk: %s", */
  /*                              ldisk->coord->id, */
  /*                              ldisk->prev->coord->id, */
  /*                              ldisk->disk->time, */
  /*                              ldisk->prev->disk->time, */
  /*                              ldisk->prev->disk->time-ldisk->disk->time, */
  /*                              ldisk->prev->disk->id, */
  /*                              ldisk->disk->id); */
  /*                 assert(FALSE); */
  /*               } */
  /*             return 0; */
  /*           } */
  /*         ldisk = ldisk->prev; */
  /*       } */
  /*     while(ldisk->prev); */
  /*   } */

  disk = tree->young_disk;
  do
    {
      if(disk->ldsk != NULL)
        {
          for(i=0;i<disk->ldsk->n_next;++i)
            {
              if(disk->ldsk->next[i]->disk->time < disk->time)
                {
                  if(exit == YES)
                    {
                      PhyML_Printf("\n. disk->time: %f [%f] disk->ldsk->next[i]->disk->time: %f [%f]",
                                   disk->time,
                                   disk->ldsk->base_line,
                                   disk->ldsk->next[i]->disk->time,
                                   disk->ldsk->next[i]->disk->ldsk->base_line);
                      assert(FALSE);
                    }
                  return(0);
                }
            }
        }
      disk = disk->prev;
    }
  while(disk);

  return(1);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void PHYREX_Store_All_Geo_Coord(t_tree *tree)
{
  int i;

  for(i=0;i<2*tree->n_otu-1;++i)
    {
      PHYREX_Store_Geo_Coord(tree->a_nodes[i]->ldsk->coord);
    }
}


/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void PHYREX_Restore_All_Geo_Coord(t_tree *tree)
{
  int i;

  for(i=0;i<2*tree->n_otu-1;++i)
    {
      PHYREX_Restore_Geo_Coord(tree->a_nodes[i]->ldsk->coord);
    }
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

void PHYREX_Store_All_Veloc(t_tree *tree)
{
  int i;

  for(i=0;i<2*tree->n_otu-1;++i)
    {
      PHYREX_Store_Veloc(tree->a_nodes[i]->ldsk->veloc);
    }
}


/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void PHYREX_Restore_All_Veloc(t_tree *tree)
{
  int i;

  for(i=0;i<2*tree->n_otu-1;++i)
    {
      PHYREX_Restore_Veloc(tree->a_nodes[i]->ldsk->veloc);
    }
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void PHYREX_Store_Veloc(t_geo_veloc *v)
{
  int i;

  for(i=0;i<v->dim;i++) v->cpy->deriv[i] = v->deriv[i];  
  v->cpy->dim = v->dim;
  strcpy(v->cpy->id,v->id);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void PHYREX_Restore_Veloc(t_geo_veloc *v)
{
  int i;
  for(i=0;i<v->dim;i++) v->deriv[i] = v->cpy->deriv[i];  
  v->dim = v->cpy->dim;
  strcpy(v->id,v->cpy->id);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

int PHYREX_Total_Number_Of_Floating_Disks(t_tree *tree)
{
  t_dsk *disk;
  int n_disks;

  disk = tree->young_disk->prev;
  n_disks = 0;
  do
    {
      if(disk->age_fixed == NO) n_disks++;
      disk = disk->prev;
    }
  while(disk);
  
  return(n_disks);
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

int PHYREX_Total_Number_Of_Single_Hit_Disks(t_tree *tree)
{
  t_dsk *disk;
  int n_hit_disks;

  assert(!(tree->young_disk->next));

  disk = tree->young_disk;
  n_hit_disks = 0;
  while(disk)
    {
      if(disk->ldsk && disk->ldsk->n_next == 1) n_hit_disks++;
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

  if(!disk->ldsk || tree->mmod->model_id == SLFV_GAUSSIAN)
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
    if(tree->n_root->v[1]->v[i] == tree->n_root)
      {
        tree->n_root->v[1]->v[i] = tree->n_root->v[2];
        break;
      }

  for(i=0;i<3;++i)
    if(tree->n_root->v[2]->v[i] == tree->n_root)
      {
        tree->n_root->v[2]->v[i] = tree->n_root->v[1];
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
  
  Update_Ancestors(tree->n_root,tree->n_root->v[2],tree->n_root->b[2],tree);
  Update_Ancestors(tree->n_root,tree->n_root->v[1],tree->n_root->b[1],tree);
  
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
  switch(tree->mmod->model_id)
    {
    case SLFV_UNIFORM:
      {
        PHYREX_Runif_Rectangle_Overlap(a_disk->ldsk,a_disk,tree->mmod);
        break;
      }
    case SLFV_GAUSSIAN :
      {
        PHYREX_Rnorm_Trunc(a_disk->ldsk,a_disk,tree->mmod);
        break;
      }
    default :
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
      /* if(tree->times->nd_t[d->num] > tree->times->nd_t[a->num]) */
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
          switch(tree->mmod->model_id)
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
            default :
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
              if(d->v[i] != a && !(a == tree->n_root && d->b[i] == tree->e_root))
                {
                  PHYREX_Tree_To_Ldsk_Post(d,d->v[i],d_disk,tree);
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
  assert(tree->mmod);
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
      else if(!((n_disks+1)%block)) dist = (phydbl *)mRealloc(dist,n_disks+block,sizeof(phydbl));
      
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
          else if(!((n_valid_disks+1)%block)) valid_disks = (t_dsk **)mRealloc(valid_disks,n_valid_disks+block,sizeof(t_dsk *));
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
  int i,j,ij_dist;

  i = Rand_Int(0,tree->n_otu-1);
  j = Rand_Int(0,tree->n_otu-1);
  ij_dist = 0;

  if(i == j) PhyML_Printf("\nxxWxx 0.0 0.0");
  else
    {

      anc = Find_Lca_Pair_Of_Nodes(tree->a_nodes[i],tree->a_nodes[j],&ij_dist,tree);
      if(anc == NULL) 
        {
          PhyML_Fprintf(stderr,"\n. %s",Write_Tree(tree));
          PhyML_Fprintf(stderr,"\n. %s %s",tree->a_nodes[i]->name,tree->a_nodes[j]->name);
          Generic_Exit(__FILE__,__LINE__,__FUNCTION__);            
        }
      
      PhyML_Printf("\nxxWxx %12f",tree->times->nd_t[anc->num]);
      dist = Euclidean_Distance(tree->a_nodes[i]->ldsk->coord,tree->a_nodes[j]->ldsk->coord);
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

  sw_lon = sw_lat = 0.0;
  ne_lon = ne_lat = 1.0;
    
  rewind(fp);
  
  found_sw = NO;
  found_ne = NO;

  for(i=0;i<tree->n_otu;i++) done[i] = NO;

  do
    {
      if(fscanf(fp,"%s",s) == EOF) break;
      
      for(i=0;i<strlen(s);++i) if(s[i] == '#') break; /* skip comment */
      if(i != strlen(s)) continue;
      
      for(i=0;i<tree->n_otu;i++) if(!strcmp(tree->a_nodes[i]->name,s)) break;

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
        }
    }
  while(1);
  
  if(found_ne == NO && (tree->mmod->model_id == SLFV_UNIFORM || tree->mmod->model_id == SLFV_GAUSSIAN))
    {
      PhyML_Fprintf(stderr,"\n. Could not find coordinates for northernmost  point.");
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
    }

  if(found_sw == NO && (tree->mmod->model_id == SLFV_UNIFORM || tree->mmod->model_id == SLFV_GAUSSIAN))
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

  if(target == young || target == old) return YES;

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

void PHYREX_Print_Disk(t_tree *tree)
{
  t_dsk *disk;


  disk = tree->n_root->ldsk->disk;

  do
    {
      PhyML_Printf("\n. Disk: %s time: %12f %2d",
                   disk->id,
                   disk->time,
                   disk->age_fixed ? 1 : 0);
      
      disk = disk->next;
    }
  while(disk);
  
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
      if(!((i+1)%block)) *list = (t_ldsk **)mRealloc(*list,i+block,sizeof(t_ldsk *));
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
          PHYREX_Free_Ldisk(disk->ldsk);
          PHYREX_Free_Disk(disk);
        }
      disk = prev;
    }
  while(disk);
  
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl PHYREX_Get_Baseline_Times(t_ldsk *ldsk, t_tree *tree)
{
  if(ldsk->next == NULL)
    {
      ldsk->base_line = ldsk->disk->time;
      return(ldsk->disk->time);
    }
  else
    {
      phydbl base_line, min_base_line;

      base_line = +INFINITY;
      min_base_line = +INFINITY;

      for(int i=0;i<ldsk->n_next;++i)
        {
          base_line = PHYREX_Get_Baseline_Times(ldsk->next[i],tree);
          if(base_line < min_base_line) min_base_line = base_line;
        }
      
      ldsk->base_line = min_base_line;
      return(min_base_line);
    }
  return(+INFINITY);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

int PHYREX_Scale_All(phydbl scale, t_dsk *start_disk, t_tree *tree)
{
  t_dsk *disk;
  int n_disks_scaled,n_nodes_scaled,n_hits_scaled,i,n_samp_disks;
  t_dsk **samp_disks_a;
  
  n_disks_scaled    = 0;
  n_nodes_scaled    = 0;
  n_hits_scaled     = 0;
  n_samp_disks      = 0;
  samp_disks_a      = NULL;
  
  disk = start_disk->prev;
  assert(disk);

  disk = start_disk->prev;
  do
    {
      if(disk->age_fixed == NO)
        {
          n_disks_scaled++;
          if(disk->ldsk && disk->ldsk->n_next > 1) n_nodes_scaled++;
          if(disk->ldsk && disk->ldsk->n_next == 1) n_hits_scaled++;
        }
      else
        {
          n_samp_disks++;
        }
      disk = disk->prev;
    }
  while(disk);


  samp_disks_a = (t_dsk **)mCalloc(n_samp_disks,sizeof(t_dsk *));

  i = 0;
  disk = start_disk->prev;
  do
    {
      if(disk->age_fixed == YES)
        {
          samp_disks_a[i] = disk;
          i++;
        }
      disk = disk->prev;
    }
  while(disk);

  disk = start_disk->prev;
  do
    {
      if(disk->age_fixed == NO)
        {
          disk->time = disk->time * scale + start_disk->time * (1.-scale);
        }
      disk = disk->prev;
    }
  while(disk);

  for(i=0;i<n_samp_disks;++i) PHYREX_Remove_Disk(samp_disks_a[i]);

  for(i=0;i<n_samp_disks;++i) PHYREX_Insert_Disk(samp_disks_a[i],tree);
  
  Free(samp_disks_a);
  
  if(tree->mmod->model_id == SLFV_GAUSSIAN || tree->mmod->model_id == SLFV_UNIFORM) return(n_disks_scaled);
  else
    {
      if(tree->times->model_id == COALESCENT)
        {
          if(tree->times->augmented_coalescent == NO)
            return(n_nodes_scaled);
          else
            return(n_hits_scaled + n_nodes_scaled);
        }
    }
  return(-1);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

t_dsk *PHYREX_Next_Floating_Disk(t_dsk *disk)
{
  t_dsk *next;
  next = disk->next;
  while(next)
    {
      if(next->age_fixed == NO) return(next);
      next = next->next;
    }
  return(NULL);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

// For a given ldsk, its age cannot be younger than the age of its baseline disk.
// The time of the baseline is the age of the oldest tip that descends from the
// node under scrutiny.
void PHYREX_Get_Baselines(t_tree *tree)
{
  PHYREX_Get_Baselines_Post(tree->n_root->ldsk,tree);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void PHYREX_Get_Baselines_Post(t_ldsk *ldsk, t_tree *tree)
{
  if(ldsk->next == NULL) ldsk->baseline = ldsk->disk;
  else
    {
      phydbl min_t;
      int i,min_idx;

      min_t = +INFINITY;
      min_idx = -1;
      
      for(i=0;i<ldsk->n_next;++i) PHYREX_Get_Baselines_Post(ldsk->next[i],tree);

      for(i=0;i<ldsk->n_next;++i)
        {
          if(ldsk->next[i]->baseline->time < min_t)
            {
              min_t = ldsk->next[i]->baseline->time;
              min_idx = i;
            }
          ldsk->baseline = ldsk->next[min_idx]->baseline;
          assert(ldsk->baseline->time > ldsk->disk->time);
        }
    }
  return;
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
  GEO_Init_Coord(mean,tree->mmod->n_dim);
  
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
    sigsq[i] = pow(Euclidean_Distance(root_disk->ldsk->coord,disk->ldsk_a[i]->coord) / (t_tip - t_root),2);

  res = Quantile(sigsq,disk->n_ldsk_a,0.5);
  
  Free(sigsq);
  
  return(res/tree->mmod->n_dim);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

/* Sum of displacements as derived from the comparison between coordinates at the two ends of each branch */
phydbl PHYREX_Realized_Displacement_Dist(short int dist_type, t_tree *tree)
{
  int i;
  phydbl disp;

  if (RRW_Is_Rw(tree->mmod) == YES && tree->mmod->integrateAncestralLocations == YES)
    RRW_Sample_Node_Locations_Joint(tree);
  if (VELOC_Is_Integrated_Velocity(tree->mmod) == YES && tree->mmod->integrateAncestralLocations == YES)
    VELOC_Sample_Node_Locations_Joint(tree);

  disp = 0.0;
  for(i=0;i<2*tree->n_otu-1;++i)
    {
      if(tree->a_nodes[i] != tree->n_root)
        {
          switch(dist_type)
            {
            case HAVERSINE :
              {
                disp += Haversine_Distance(tree->a_nodes[i]->anc->ldsk->coord,tree->a_nodes[i]->ldsk->coord);
                break;
              }
            case MANHATTAN :
              {
                disp += Manhattan_Distance(tree->a_nodes[i]->anc->ldsk->coord,tree->a_nodes[i]->ldsk->coord);
                break;
              }
            case EUCLIDEAN :
              {
                disp += Euclidean_Distance(tree->a_nodes[i]->anc->ldsk->coord,tree->a_nodes[i]->ldsk->coord);
                break;
              }
            default : assert(false);
            }
        }
    }
  
  return(disp);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

/* Sum of actual distances travelled between the two ends of each branch. May differ from dispersal distance when "jumps" take place along edges, such as with the SLFV model */
phydbl PHYREX_Realized_Dispersal_Dist(short int dist_type, t_tree *tree)
{
  t_dsk *disk,*root_disk;
  phydbl dist;
  int i;

  if (RRW_Is_Rw(tree->mmod) == YES && tree->mmod->integrateAncestralLocations == YES)
    RRW_Sample_Node_Locations_Joint(tree);
  if (VELOC_Is_Integrated_Velocity(tree->mmod) == YES && tree->mmod->integrateAncestralLocations == YES)
    VELOC_Sample_Node_Locations_Joint(tree);

  disk = tree->young_disk;
  while(disk->prev) disk = disk->prev;
  root_disk = disk;

  dist = 0.0;
  disk = root_disk;
  do
    {
      if(disk->ldsk != NULL)
        {
          for(i=0;i<disk->ldsk->n_next;++i)
            {
              switch(dist_type)
                {
                case HAVERSINE :
                  {
                    dist += Haversine_Distance(disk->ldsk->coord,disk->ldsk->next[i]->coord);
                    break;
                  }
                case MANHATTAN :
                  {
                    dist += Manhattan_Distance(disk->ldsk->coord,disk->ldsk->next[i]->coord);
                    break;
                  }
                case EUCLIDEAN :
                  {
                    dist += Euclidean_Distance(disk->ldsk->coord,disk->ldsk->next[i]->coord);
                    break;
                  }
                default : assert(false);
                }
            }
        }
      disk = disk->next;
    }
  while(disk);
  
  return(dist);
}

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

// Mean Haversine distance (in km) per year (measured on all paths between each internal node and its children) 
phydbl PHYREX_Realized_Dispersal_Dist_Alt(t_tree *tree)
{
  t_dsk *disk,*root_disk;
  phydbl dt,sum_dist;
  phydbl mean,var;
  int n,i,j;
  t_geo_coord *new_coord;

  
  new_coord = GEO_Make_Geo_Coord(tree->mmod->n_dim);
  GEO_Init_Coord(new_coord,tree->mmod->n_dim);
  
  disk = tree->young_disk;
  while(disk->prev) disk = disk->prev;
  root_disk = disk;


  n = 0;
  sum_dist = 0.0;
  disk = root_disk;
  do
    {
      if(disk->ldsk != NULL)
        {
          for(i=0;i<disk->ldsk->n_next;++i)
            {
              dt = disk->ldsk->next[i]->disk->time - disk->time;

              if(dt > 1.0)
                {
                  for(j=0;j<tree->mmod->n_dim;++j)
                    {
                      mean = disk->ldsk->coord->lonlat[j] + (disk->ldsk->next[i]->coord->lonlat[j] - disk->ldsk->coord->lonlat[j]) * 1.0 / dt;
                      var  = tree->mmod->sigsq[j] * (dt - 1.)/dt;
                      new_coord->lonlat[j] = Rnorm(mean,sqrt(var));
                    }
                  sum_dist += Haversine_Distance(disk->ldsk->coord,new_coord);
                  n++;
                }
            }
        }
      disk = disk->next;
    }
  while(disk);

  Free_Geo_Coord(new_coord);

  if(n == 0) return(0.0);
  else return(sum_dist / n);
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
        pow(Euclidean_Distance(disk->ldsk_a[i]->coord,
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
      sumdist += pow(Euclidean_Distance(disk->ldsk_a[i]->coord,
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
      dist = Euclidean_Distance(ldsk->coord,tree->young_disk->ldsk_a[i]->coord);
      mean += pow(dist,2);
      
      /* PhyML_Printf("\n. %3d (%12f %12f) (%12f %12f) dist: %12f mean: %12f", */
      /*              i, */
      /*              tree->young_disk->ldsk_a[i]->coord->lonlat[0], */
      /*              tree->young_disk->ldsk_a[i]->coord->lonlat[1], */
      /*              ldsk->coord->lonlat[0],                   */
      /*              ldsk->coord->lonlat[1], */
      /*              Euclidean_Distance(ldsk->coord,tree->young_disk->ldsk_a[i]->coord), */
      /*              mean); */
                   
    }
  
  return(mean/(tree->young_disk->n_ldsk_a*tree->mmod->n_dim));
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
              t_dsk *prev,*next;
              prev = disk->prev;
              next = disk->next;
              assert(next && prev);
              PHYREX_Remove_Disk(disk);
              PHYREX_Update_Lindisk_List_Range(next,prev,tree);
              PHYREX_Free_Disk(disk);
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

              PHYREX_Update_Lindisk_List_Range(young_ldsk->disk,old_ldsk->disk,tree);

              PHYREX_Free_Ldisk(disk->ldsk);
              PHYREX_Free_Disk(disk);
            }
        }

      disk = dum;
    }
  while(disk);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
  
phydbl PHYREX_Path_Logdensity(t_ldsk *young, t_ldsk *old, phydbl *sd, t_tree *tree)
{
  switch(tree->mmod->model_id)
    {
    case SLFV_UNIFORM : case SLFV_GAUSSIAN : 
      {
        return(SLFV_Path_Logdensity(young,old,sd,tree));
        break;
      }      
    case RW : case RRW_GAMMA : case RRW_LOGNORMAL : 
    case IBM : case RIBM : case IWNc : case IWNu : case RIWNc : case RIWNu : case IOU : 
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

void PHYREX_Sample_Path(t_ldsk *young, t_ldsk *old, phydbl *sd, phydbl *global_hr, t_tree *tree)
{
  switch(tree->mmod->model_id)
    {
    case SLFV_UNIFORM : case SLFV_GAUSSIAN : 
      {
        SLFV_Sample_Path(young,old,sd,global_hr,tree);
        break;
      }
      
    case RW : case RRW_GAMMA : case RRW_LOGNORMAL : 
    case IBM : case RIBM : case IWNc : case IWNu : case RIWNc : case RIWNu : case IOU :
      {
        assert(FALSE);
        break;
      }
    default : assert(FALSE);
    }
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

t_ldsk *PHYREX_Generate_Path(t_ldsk *young, t_ldsk *old, phydbl n_evt, phydbl *sd, t_tree *tree)
{
  switch(tree->mmod->model_id)
    {
    case SLFV_UNIFORM : case SLFV_GAUSSIAN : 
      {
        return(SLFV_Generate_Path(young,old,n_evt,sd,tree));
        break;
      }      
    case RW : case RRW_GAMMA : case RRW_LOGNORMAL :  
    case IBM : case RIBM : case IWNc : case IWNu : case RIWNc : case RIWNu : case IOU : 
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

t_tree *PHYREX_Simulate(int n_otu, int n_sites, phydbl w, phydbl h, phydbl  lbda, phydbl rad, phydbl mu, int r_seed, int modid, char *datadir)
{
  switch(modid)
    {
    case SLFV_UNIFORM : case SLFV_GAUSSIAN : 
      {
        return(SLFV_Simulate(n_otu,n_sites,w,h,lbda,rad,mu,r_seed,datadir));
        break;
      }
      
    case RW : case RRW_GAMMA : case RRW_LOGNORMAL :
    case IBM : case RIBM : case IWNc : case IWNu : case RIWNc : case RIWNu : case IOU : 
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
  switch(tree->mmod->model_id)
    {
    case SLFV_UNIFORM : case SLFV_GAUSSIAN : 
      {
        SLFV_Simulate_Backward_Core(disk,avoid_multiple_mergers,tree);
        break;
      }
      
    case RW : case RRW_GAMMA : case RRW_LOGNORMAL :
    case IBM : case RIBM : case IWNc : case IWNu : case RIWNc : case RIWNu : case IOU : 
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
  switch(tree->mmod->model_id)
    {
    case SLFV_UNIFORM : case SLFV_GAUSSIAN : 
      {
        return(SLFV_Update_Sigsq(tree));
        break;
      }

    case RW : case RRW_GAMMA : case RRW_LOGNORMAL :
    case IBM : case RIBM : case IWNc : case IWNu : case RIWNc : case RIWNu : case IOU : 
      {
        return(tree->mmod->sigsq[0]);
        break;
      }
      
    default : assert(FALSE);

    }

  return(-1.);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

/* Return the edge sitting between two ldsks. Should work for multifurcating trees */
/* Make sure node ranks are up-to-date */
t_edge *PHYREX_Edge_Between_Two_Ldsks(t_ldsk *a, t_ldsk *d, t_tree *tree)
{
  if(d->n_next == 2) // d is has out-degree 2
    {
      assert(d->nd->anc->ldsk == a);
      return(d->nd->b_anc);
    }
  else
    {
      t_node *n;
      n = d->nd;
      while(n->ldsk == d)
        {
          if(n->prev != d->nd) break;
          n = n->prev;
        }
      assert(n->anc->ldsk == a);
      return(n->b_anc);
    }
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void PHYREX_Update_Node_Times_Given_Disks(t_tree *tree)
{
  for(int i=0;i<2*tree->n_otu-1;++i) tree->times->nd_t[i] = tree->a_nodes[i]->ldsk->disk->time;
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void PHYREX_Update_Disk_Times_Given_Nodes(t_tree *tree)
{
  for(int i=0;i<2*tree->n_otu-1;++i) tree->a_nodes[i]->ldsk->disk->time = tree->times->nd_t[i];
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
/* Make sure node times and tree topology are in sync with ldsk structure */
void PHYREX_Update_Edge_Rates_Given_Ldsks(t_tree *tree)
{
  t_dsk *disk;
  int i;
  t_node *n;
  
  disk = tree->young_disk;
  if(disk == NULL) return;

  Get_Node_Ranks_From_Times(tree);

  do
    {
      if(disk->age_fixed == YES)
        {
          for(i=0;i<disk->n_ldsk_a;++i)
            {
              if(disk->ldsk_a[i]->nd != NULL && disk->ldsk_a[i]->nd->tax == YES)
                {
                  tree->rates->br_r[disk->ldsk_a[i]->nd->num] = disk->ldsk_a[i]->rr; 
                }
            }
        }
      else if(disk->ldsk->n_next > 1)
        {
          n = disk->ldsk->nd;
          while(n && n->ldsk == disk->ldsk) { tree->rates->br_r[n->num] = disk->ldsk->rr; n = n->rk_next; }
          n = disk->ldsk->nd;
          while(n && n->ldsk == disk->ldsk) { tree->rates->br_r[n->num] = disk->ldsk->rr; n = n->rk_prev; }          
        }    
      disk = disk->prev;
    }
  while(disk);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
/* Make sure node times and tree topology are in sync with ldsk structure */
void PHYREX_Update_Edge_Sigsq_Given_Ldsks(t_tree *tree)
{
  t_dsk *disk;
  int i;
  t_node *n;
  
  disk = tree->young_disk;
  if(disk == NULL) return;

  Get_Node_Ranks_From_Times(tree);

  do
    {
      if(disk->age_fixed == YES)
        {
          for(i=0;i<disk->n_ldsk_a;++i)
            {
              if(disk->ldsk_a[i]->nd != NULL && disk->ldsk_a[i]->nd->tax == YES)
                {
                  tree->mmod->sigsq_scale[disk->ldsk_a[i]->nd->num] = disk->ldsk_a[i]->sigsq; 
                }
            }
        }
      else if(disk->ldsk->n_next > 1)
        {
          n = disk->ldsk->nd;
          while(n && n->ldsk == disk->ldsk) { tree->mmod->sigsq_scale[n->num] = disk->ldsk->sigsq; n = n->rk_next; }
          n = disk->ldsk->nd;
          while(n && n->ldsk == disk->ldsk) { tree->mmod->sigsq_scale[n->num] = disk->ldsk->sigsq; n = n->rk_prev; }          
        }    
      disk = disk->prev;
    }
  while(disk);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
/* Make sure node times and tree topology are in sync with ldsk structure */
void PHYREX_Update_Ldsk_Rates_Given_Edges(t_tree *tree)
{
  int i;
  
  if(tree->young_disk == NULL) return;

  for(i=0;i<2*tree->n_otu-1;++i)
    {
      PHYREX_Update_Ldsk_Rates_Given_One_Edge(tree->a_nodes[i],tree);
    }
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
/* Make sure node times and tree topology are in sync with ldsk structure */
void PHYREX_Update_Ldsk_Sigsq_Given_Edges(t_tree *tree)
{
  int i;
  
  if(tree->young_disk == NULL) return;

  for(i=0;i<2*tree->n_otu-1;++i)
    {
      PHYREX_Update_Ldsk_Sigsq_Given_One_Edge(tree->a_nodes[i],tree);
    }
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void PHYREX_Update_Ldsk_Rates_Given_One_Edge(t_node *d, t_tree *tree)
{
  t_ldsk *ldsk;
  t_node *a;

  if(tree->young_disk == NULL) return;
  
  a = d->anc;
  
  ldsk = d->ldsk;
  assert(ldsk);
  
  while(ldsk && ldsk->nd != a)
    {
      ldsk->rr = tree->rates->br_r[d->num];
      ldsk = ldsk->prev;
    }
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void PHYREX_Update_Ldsk_Sigsq_Given_One_Edge(t_node *d, t_tree *tree)
{
  t_ldsk *ldsk;
  t_node *a;

  if(tree->young_disk == NULL) return;
  
  a = d->anc;
  
  ldsk = d->ldsk;
  assert(ldsk);

  while(ldsk && ldsk->nd != a)
    {
      ldsk->sigsq = tree->mmod->sigsq_scale[d->num];
      ldsk = ldsk->prev;
    }
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void PHYREX_Duplicate_Ldsk_Struct(t_tree *from, t_tree *where)  
{
  t_dsk *disk;
  int i,j;
  
  disk = from->n_root->ldsk->disk;

  do
    {
      disk->img = PHYREX_Make_Disk_Event(from->mmod->n_dim,from->n_otu);
      PHYREX_Init_Disk_Event(disk->img,from->mmod->n_dim,from->mmod);

      disk->img->n_ldsk_a = disk->n_ldsk_a;
      disk->img->age_fixed = disk->age_fixed;
      disk->img->time = disk->time;

      for(i=0;i<from->mmod->n_dim;++i) disk->img->centr->lonlat[i] = disk->centr->lonlat[i]; 

      if(disk->ldsk != NULL)
        {
          disk->img->ldsk = PHYREX_Make_Lindisk_Node(from->mmod->n_dim);
          PHYREX_Init_Lindisk_Node(disk->img->ldsk,disk->img,from->mmod->n_dim);
          
          disk->ldsk->img = disk->img->ldsk;
          disk->ldsk->img->disk = disk->img;

          for(i=0;i<from->mmod->n_dim;++i) disk->ldsk->img->coord->lonlat[i] = disk->ldsk->coord->lonlat[i];
          for(i=0;i<disk->ldsk->n_next;++i) PHYREX_Make_Lindisk_Next(disk->ldsk->img);
          
          if(disk->ldsk->prev != NULL) disk->ldsk->img->prev = disk->ldsk->prev->img;
          else disk->ldsk->img->prev = NULL;
          
          if(disk->ldsk->nd != NULL)
            {
              disk->ldsk->img->nd = where->a_nodes[disk->ldsk->nd->num];
              where->a_nodes[disk->ldsk->nd->num]->ldsk = disk->ldsk->img;
            }
          
          disk->ldsk->img->rr = disk->ldsk->rr;
          disk->ldsk->img->sigsq = disk->ldsk->sigsq;
        }
      
      if(disk->age_fixed == YES)
        {
          for(i=0;i<from->n_otu;++i)
            {
              if(disk->ldsk_a[i] && disk->ldsk_a[i]->n_next == 0 && disk->ldsk_a[i]->disk == disk)
                {
                  assert(disk->ldsk_a[i]->nd != NULL);

                  disk->img->ldsk_a[i] = PHYREX_Make_Lindisk_Node(from->mmod->n_dim);
                  PHYREX_Init_Lindisk_Node(disk->img->ldsk_a[i],disk->img,from->mmod->n_dim);
                  
                  disk->ldsk_a[i]->img = disk->img->ldsk_a[i];
                  
                  disk->img->ldsk_a[i]->prev = disk->ldsk_a[i]->prev->img;
                  disk->img->ldsk_a[i]->next = NULL;
                  disk->img->ldsk_a[i]->n_next = 0;
                  disk->img->ldsk_a[i]->disk = disk->img;
                  
                  for(j=0;j<from->mmod->n_dim;++j) disk->img->ldsk_a[i]->coord->lonlat[j] = disk->ldsk_a[i]->coord->lonlat[j];

                  disk->img->ldsk_a[i]->nd = where->a_nodes[disk->ldsk_a[i]->nd->num];
                  where->a_nodes[disk->ldsk_a[i]->nd->num]->ldsk = disk->img->ldsk_a[i];

                  disk->ldsk_a[i]->img->rr = disk->ldsk_a[i]->rr;
                  disk->ldsk_a[i]->img->sigsq = disk->ldsk_a[i]->sigsq;
                }
            }
        }
        
      if(disk->prev != NULL)
        {
          disk->img->prev = disk->prev->img;
          disk->prev->img->next = disk->img;
        }
      else disk->img->prev = NULL;

      disk = disk->next;
    }
  while(disk != NULL);
  
  
  disk = from->n_root->ldsk->disk;  
  do
    {
      if(disk->ldsk != NULL)
        {
          assert(disk->ldsk->img->n_next == disk->ldsk->n_next); 
          for(i=0;i<disk->ldsk->n_next;++i)
            {
              assert(disk->ldsk->next[i]->img != NULL);
              disk->ldsk->img->next[i] = disk->ldsk->next[i]->img;
            }
        }
      disk = disk->next;
    }
  while(disk);
  
  where->old_samp_disk = from->old_samp_disk->img;
  where->young_disk = from->young_disk->img;
  PHYREX_Update_Lindisk_List(where);
  PHYREX_Ldsk_To_Tree(where);  
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl PHYREX_Get_Posterior(t_tree *tree)
{
  phydbl lnP;

  lnP = 0.0;

  // Likelihoods
  lnP += tree->c_lnL;
  lnP += tree->rates->c_lnL;
  lnP += tree->times->c_lnL;
  lnP += tree->mmod->c_lnL;

  // Priors
  lnP += tree->rates->c_lnP;
  lnP += tree->times->c_lnP;
  lnP += tree->mmod->c_lnP;

  
  return(lnP);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void PHYREX_Evolve_All(t_tree *tree)
{
  int i;
  RRW_Generate(tree);
  PhyML_Printf("\n. Done with RRW_Generate");

  for(i=0;i<tree->n_otu;++i)
    {
      PhyML_Printf("\n<clade id=\"clad%d\">",i+1);
      PhyML_Printf("\n\t<taxon value=\"%s\"/>",tree->a_nodes[i]->name);
      PhyML_Printf("\n</clade>");

      PhyML_Printf("\n<calibration id=\"cal%d\">",i+1);
      PhyML_Printf("\n\t<lower>%f</lower>",tree->times->nd_t[tree->a_nodes[i]->num]);
      PhyML_Printf("\n\t<upper>%f</upper>",tree->times->nd_t[tree->a_nodes[i]->num]);
      PhyML_Printf("\n\t<appliesto clade.id=\"clad%d\"/>",i+1);
      PhyML_Printf("\n</calibration>");
    }

  PhyML_Printf("#traits\t latMg\t longMg");
  PhyML_Printf("\n|SouthWest| -26 +40");
  PhyML_Printf("\n|NorthEast| -11 +53");
  
  for(i=0;i<tree->n_otu;++i)
    {
      PhyML_Printf("\n%s\t%f\t%f",
                   tree->a_nodes[i]->name,
                   tree->a_nodes[i]->ldsk->coord->lonlat[1],
                   tree->a_nodes[i]->ldsk->coord->lonlat[0]);
    }

  PhyML_Printf("\n. dispDist: %f",PHYREX_Realized_Dispersal_Dist(HAVERSINE,tree));
  PhyML_Printf("\n. sigSqLon: %f",tree->mmod->sigsq[0]);
  PhyML_Printf("\n. sigSqLat: %f",tree->mmod->sigsq[1]);
  PhyML_Printf("\n. neff: %f",tree->times->scaled_pop_size);
  PhyML_Printf("\n. growth: %f",tree->times->neff_growth);
  PhyML_Printf("\n. root time: %f",tree->times->nd_t[tree->n_root->num]);
  PhyML_Printf("\n. root lon: %f",tree->n_root->ldsk->coord->lonlat[0]);
  PhyML_Printf("\n. root lat: %f",tree->n_root->ldsk->coord->lonlat[1]);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void PHYREX_Record_Disk_Times(t_tree *tree)
{
  t_dsk *disk;

  disk = tree->n_root->ldsk->disk;

  do
    {
      disk->time_bkp = disk->time;
      disk = disk->next;
    }
  while(disk);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void PHYREX_Restore_Disk_Times(t_tree *tree)
{
  t_dsk *disk;
  int i,n_samp_disks;
  t_dsk **samp_disks_a;
    
  n_samp_disks      = 0;
  samp_disks_a      = NULL;
  
  disk = tree->young_disk->prev;
  assert(disk);

  disk = tree->young_disk->prev;
  do
    {
      if(disk->age_fixed == YES)
        {
          n_samp_disks++;
        }
      disk = disk->prev;
    }
  while(disk);


  samp_disks_a = (t_dsk **)mCalloc(n_samp_disks,sizeof(t_dsk *));

  i = 0;
  disk = tree->young_disk->prev;
  do
    {
      if(disk->age_fixed == YES)
        {
          samp_disks_a[i] = disk;
          i++;
        }
      disk = disk->prev;
    }
  while(disk);

  disk = tree->young_disk->prev;
  do
    {
      if(disk->age_fixed == NO)
        {
          disk->time = disk->time_bkp;
        }
      disk = disk->prev;
    }
  while(disk);

  for(i=0;i<n_samp_disks;++i) PHYREX_Remove_Disk(samp_disks_a[i]);

  for(i=0;i<n_samp_disks;++i) PHYREX_Insert_Disk(samp_disks_a[i],tree);
  
  Free(samp_disks_a);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void PHYREX_Swap_Coords(t_ldsk *a, t_ldsk *b, t_tree *tree)
{
  phydbl buff;
  
  for(int i=0;i<tree->mmod->n_dim;++i)
    {
      buff = a->coord->lonlat[i];
      a->coord->lonlat[i] = b->coord->lonlat[i];
      b->coord->lonlat[i] = buff;

      buff = a->veloc->deriv[i];
      a->veloc->deriv[i] = b->veloc->deriv[i];
      b->veloc->deriv[i] = buff;
    }
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void PHYREX_Check_Disk_Times(t_tree *tree)
{
  t_dsk *disk;

  disk = tree->young_disk;

  do
    {
      if(Are_Equal(disk->time,disk->time_bkp,1.E-10) == NO)
        {
          PhyML_Printf("\n. Disk %s time inconsistency detected (t: %f t.bkup: %f -- diff: %g).",disk->id,disk->time,disk->time_bkp,disk->id,disk->time-disk->time_bkp);
        }
      disk = disk->prev;
    }
  while(disk);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void PHYREX_Exchange_Ldsk(t_ldsk *a, t_ldsk *d, t_ldsk *w, t_ldsk *v, int aw, int dv, t_tree *tree)
{
  t_ldsk *zv,*zw;
  
  assert(!(a->disk->time > d->disk->time));

  zv = v;
  while(zv->n_next == 1) zv = zv->next[0]; 
  assert(zv->nd);
  
  zw = w;
  while(zw->n_next == 1) zw = zw->next[0]; 
  assert(zw->nd);
    
  /* Connect ldsks */
  v->prev     = a;
  a->next[aw] = v;
  w->prev     = d;
  d->next[dv] = w;

  Exchange_Nodes(a->nd,d->nd,zw->nd,zv->nd,tree);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

short int PHYREX_Check_Lk(t_tree *tree)
{
  if(tree->mmod->safe_phyrex == YES)
    {
      short int failed = NO;
      
      /* if(Are_Equal(RATES_Realized_Substitution_Rate(tree),tree->rates->clock_r,1.E-1)== NO) */
      /*   { */
      /*     PhyML_Fprintf(stderr,"\n. Problem detected with move %s",tree->mcmc->move_name[move]); */
      /*     assert(false); */
      /*   } */
      
      phydbl c_lnL = tree->c_lnL;
      RATES_Update_Edge_Lengths(tree);
      if(tree->eval_alnL == YES) Lk(NULL,tree);

      if(Are_Equal(c_lnL,tree->c_lnL,1.E-3) == NO)
        {
          PhyML_Fprintf(stdout,"\n. a_lnL: %f -> %f [%g]",c_lnL,tree->c_lnL,c_lnL-tree->c_lnL);
          failed = YES;
        }
      
      if(tree->mmod->use_locations == YES)
        {
          phydbl g_lnL = tree->mmod->c_lnL;
          LOCATION_Lk(NULL,tree);
          
          if(Are_Equal(g_lnL,tree->mmod->c_lnL,1.E-3) == NO)
            {
              PhyML_Fprintf(stdout,"\n. g_lnL: %f -> %f [%g]",g_lnL,tree->mmod->c_lnL,g_lnL-tree->mmod->c_lnL);
              failed = YES;
            }
        }
      
      phydbl r_lnL = tree->rates->c_lnL;
      RATES_Lk(tree);
      if(Are_Equal(r_lnL,tree->rates->c_lnL,1.E-3) == NO)
        {
          PhyML_Fprintf(stderr,"\n. r_lnL: %f -> %f [%g]",r_lnL,tree->rates->c_lnL,r_lnL-tree->rates->c_lnL);
          failed = YES;
        }
      
      phydbl t_lnL = tree->times->c_lnL;
      TIMES_Lk(tree);
      if(Are_Equal(t_lnL,tree->times->c_lnL,1.E-3) == NO)
        {
          PhyML_Fprintf(stdout,"\n. t_lnL: %f -> %f [%g]",t_lnL,tree->times->c_lnL,t_lnL-tree->times->c_lnL);
          failed = YES;
        }
      
      
      /* RATES_Check_Rates(tree); */
      
      /* phydbl subsrate = RATES_Realized_Substitution_Rate(tree); */
      /* if(Are_Equal(subsrate,tree->rates->clock_r,1.E-8) == NO) */
      /*   { */
      /*     PhyML_Fprintf(stderr,"\n. Problem detected with move %s",tree->mcmc->move_name[move]); */
      /*     PhyML_Fprintf(stderr,"\n. rate : %f -> %f",subsrate,tree->rates->clock_r); */
      /*     Generic_Exit(__FILE__,__LINE__,__FUNCTION__); */
      /*   } */
      
      
      for(int i=0;i<2*tree->n_otu-1;++i)
        {
          if(Are_Equal(tree->a_nodes[i]->ldsk->rr,tree->rates->br_r[tree->a_nodes[i]->num],1.E-6) == NO)
            {
              PhyML_Fprintf(stderr,"\n. Problem detected with move %s",tree->mcmc->move_name[tree->mcmc->move_idx]);
              PhyML_Fprintf(stderr,"\n. rate : %f -> %f",tree->a_nodes[i]->ldsk->rr,tree->rates->br_r[tree->a_nodes[i]->num]);
              failed = YES;
            }
        }
      
      for(int i=0;i<2*tree->n_otu-1;++i)
        {
          if(Are_Equal(tree->a_nodes[i]->ldsk->disk->time,
                       tree->times->nd_t[tree->a_nodes[i]->num],1.E-6) == NO)
            {
              PhyML_Fprintf(stderr,"\n. Problem detected with move %s",tree->mcmc->move_name[tree->mcmc->move_idx]);
              PhyML_Fprintf(stderr,"\n. time : %f -> %f",
                            tree->a_nodes[i]->ldsk->disk->time,
                            tree->times->nd_t[tree->a_nodes[i]->num]);
              failed = YES;
            }
        }

      assert(PHYREX_Check_Struct(tree,YES));

      if(failed == YES)
        {
          PhyML_Fprintf(stderr,"\n. Problem detected with move %s; mcmc step # %d",tree->mcmc->move_name[tree->mcmc->move_idx],tree->mcmc->run);          
          Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
          return(0);
        }
    }

  return(1);
}
  
/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void PHYREX_Check_Lk_Nodes(t_tree *tree)
{
  tree->contmod->both_sides[LOCATION] = YES;
  tree->contmod->both_sides[VELOCITY] = YES;
  LOCATION_Lk(NULL,tree);

  PhyML_Printf("\n. Reference lk: %f",tree->mmod->c_lnL);
  
  for(int i=0; i < 2*tree->n_otu-1; ++i)
  {
    LOCATION_Lk(tree->a_nodes[i],tree);
    PhyML_Printf("\n. Node %d lk: %f",i,LOCATION_Lk(tree->a_nodes[i],tree));
  }
  
  tree->contmod->both_sides[LOCATION] = NO;
  tree->contmod->both_sides[VELOCITY] = NO;
}
