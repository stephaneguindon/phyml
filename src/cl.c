/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include "cl.h"


/*********************************************************/
/**
* Fill the Option fields, with the argc array
*/
int Read_Command_Line(option *io, int argc, char **argv)
{
  int c;
  int idx;
  int i;
  int writemode;
  short int opt_m;
  
  writemode = WRITE;
  opt_m = 0;
  
  if(argc == 1) Exit("\n. No argument was passed to the program. Please check the documentation. \n");
  PhyML_Printf("",writemode);
  
  struct option longopts[] =
    {
      {"n_rgrft",           required_argument,NULL,0},
      {"n_globl",           required_argument,NULL,1},
      {"max_dist",          required_argument,NULL,2},
      {"n_optim",           required_argument,NULL,3},
      {"n_best",            required_argument,NULL,4},
      {"model",             required_argument,NULL,5},
      {"search",            required_argument,NULL,6},
      {"datatype",          required_argument,NULL,7},
      {"multiple",          required_argument,NULL,8},
      {"input",             required_argument,NULL,9},
      {"bootstrap",         required_argument,NULL,10},
      {"ts/tv",             required_argument,NULL,11},
      {"nclasses",          required_argument,NULL,12},
      {"pinv",              required_argument,NULL,13},
      {"alpha",             required_argument,NULL,14},
      {"inputtree",         required_argument,NULL,15},
      {"min_diff_lk_local", required_argument,NULL,16},
      {"min_diff_lk_global",required_argument,NULL,17},
      {"steph_spr",         no_argument,NULL,18},
      {"brent_it_max",      required_argument,NULL,19},
      {"rand_start",        no_argument,NULL,20},
      {"n_rand_starts",     required_argument,NULL,21},
      {"sequential",        no_argument,NULL,22},
      {"inside_opt",        no_argument,NULL,23},
      {"p_moves",           required_argument,NULL,24},
      {"fast_nni",          no_argument,NULL,25},
      {"g_pars",            no_argument,NULL,26},
      {"r_seed",            required_argument,NULL,27},
      {"collapse_boot",     required_argument,NULL,28},
      {"random_boot",       required_argument,NULL,29},
      {"print_trace",       no_argument,NULL,30},
      {"print_site_lnl",    no_argument,NULL,31},
      {"print_site_lk",    no_argument,NULL,31},
      {"cov",               no_argument,NULL,32},
      {"cov_delta",         required_argument,NULL,33},
      {"cov_alpha",         required_argument,NULL,34},
      {"cov_ncats",         required_argument,NULL,35},
      {"ps",                no_argument,NULL,36},
      {"cov_free",          no_argument,NULL,37},
      {"no_gap",            no_argument,NULL,38},
      {"n_rr_branch",       required_argument,NULL,39},
      {"append",            no_argument,NULL,40},
      {"no_five_branch",    no_argument,NULL,41},
      {"pars_thresh",       required_argument,NULL,42},
      {"min_diff_lk_move",  required_argument,NULL,43},
      {"hybrid",            no_argument,NULL,44},
      {"use_median",        no_argument,NULL,45},
      {"run_id",            required_argument,NULL,46},
      {"pars",              no_argument,NULL,47},
      {"quiet",             no_argument,NULL,48},
      {"version",           no_argument,NULL,49},
      {"calibration_file",    required_argument,NULL,50},
      {"calibration",         required_argument,NULL,50},
      {"clade_file",          required_argument,NULL,50},
      {"boot_progress_every", required_argument,NULL,51},
      {"aa_rate_file",        required_argument,NULL,52},
      {"chain_len",           required_argument,NULL,53},
      {"sample_freq",         required_argument,NULL,54},
      {"burnin",              required_argument,NULL,55},
      {"no_memory_check",     no_argument,NULL,56},
      {"no_colalias",         no_argument,NULL,57},
      {"alias_subpatt",       no_argument,NULL,58},      
      {"no_sequences",        no_argument,NULL,59},      
      {"prior",               no_argument,NULL,59},      
      {"fastlk",              no_argument,NULL,60},      
      {"free_rates",          no_argument,NULL,61},
      {"freerates",           no_argument,NULL,61},
      {"freerate",            no_argument,NULL,61},
      {"free_rate",            no_argument,NULL,61},
      {"is",                  no_argument,NULL,62},
      // no 63 since it corresponds to character '?' 
      {"rate_model",          required_argument,NULL,64},
      {"ratemodel",           required_argument,NULL,64},
      {"log_l",               no_argument,NULL,65},
      {"gamma_lens",          no_argument,NULL,66},
      {"il",                  no_argument,NULL,66},
      {"codpos",              required_argument,NULL,67},
      {"constraint_file",     required_argument,NULL,68},
      {"constraint_tree",     required_argument,NULL,68},
      {"help",                no_argument,NULL,69},
      {"mutmap",              no_argument,NULL,70},
      {"parvals",             required_argument,NULL,71},
      {"constrained_lens",    no_argument,NULL,72},
      {"xml",                 required_argument,NULL,73},
      {"l_var",               required_argument,NULL,74},
#ifdef BEAGLE
      {"beagle_resource",     required_argument,NULL,75},
#endif
      {"ancestral",           no_argument,NULL,76},
      {"anc",                 no_argument,NULL,76},
      {"coord_file",          required_argument,NULL,77},
      {"json_trace",          no_argument,NULL,78},
      {"weights",             required_argument,NULL,79},
      {"tbe",                 no_argument,NULL,80},
      {"leave_duplicates",    no_argument,NULL,81},
      {"precision",           required_argument,NULL,82},
      {"l_min",               required_argument,NULL,83},
      {"print_mat_and_exit",  no_argument,NULL,84},
      {0,0,0,0}
    };

  io->datatype = UNDEFINED;

  #ifndef PHYML
  int open_ps_file = 0;
  #endif

  idx=-1;

    do
    {     
      c = getopt_long(argc,argv,"qi:d:t:m:b:n:f:zk:v:c:a:u:ho:s:x:g:l:ep",longopts,&idx);

      switch(c)
	{
        case 84 :
          {
            io->print_mat_and_exit = YES;
            break;
          }
        case 83 :
          {
            io->mod->l_min = String_To_Dbl(optarg);
            break;
          }
	case 82 :
	  {
	    if ((!atoi(optarg)) || (atoi(optarg) < 1) || (atoi(optarg) > DECIMAL_DIG -3))
	      {
		    PhyML_Printf("\n. The number of digits must be [1 - %d]\n", DECIMAL_DIG -3);
		    Exit("\n");
	      }
	    else 
	      {
		io->precision = atoi(optarg);
	      }
	    break;
	  }
	case 81 :
	  {
	    io->leave_duplicates = YES;
	    break;
	  }
	case 80 :
	  {
	    io->do_tbe = YES;
	    io->do_boot = NO;
            io->do_alrt = NO;
	    break;
	  }
        case 79:
          {
            io->has_io_weights = YES;
            strcpy(io->weight_file, optarg);
            break;
          }
        case 78:
          {
	    io->print_json_trace = YES;
	    break;            
          }

        case 77:
          {
	    char *tmp;
	    tmp = (char *)mCalloc(T_MAX_FILE, sizeof(char));

	    if(strlen(optarg) > T_MAX_FILE -11)
	      {
		char choix;
		strcpy (tmp, "\n. The file name'");
		strcat (tmp, optarg);
		strcat (tmp, "' is too long.\n");
		PhyML_Printf("%s",tmp);
		PhyML_Printf("\n. Type any key to exit.\n");
		if(!scanf("%c",&choix)) Exit("\n");
		Exit("\n");
	      }
	    else if (!Filexists (optarg))
	      {
		char choix;
		strcpy (tmp, "\n. The file '");
		strcat (tmp, optarg);
		strcat (tmp, "' doesn't exist.\n");
		PhyML_Printf("%s",tmp);
		PhyML_Printf("\n. Type any key to exit.\n");
		if(!scanf("%c",&choix)) Exit("\n");
		Exit("\n");
	      }
	    else
	      {
                strcpy(io->in_coord_file, optarg);
                io->fp_in_coord = Openfile(io->in_coord_file,READ);
	      }
	    Free(tmp);
            break;
          }

	case 76:
          {
            io->ancestral = YES;
            break;
          }
#ifdef BEAGLE
	case 75:
          {
            io->beagle_resource = (int)atoi(optarg);
            break;
          }
#endif
	case 74:
          {
            io->mod->l_var_sigma = String_To_Dbl(optarg);
            break;
          }
	case 73:
	  {
#ifdef INVITEE
            
            Free_Optimiz(io->mod->s_opt);
            Free_Model_Basic(io->mod);
            Free_Input(io);
            PhyTime_XML(optarg);
            return 0;

#elif defined(PHYML)
           
            Free_Optimiz(io->mod->s_opt);
            Free_Model_Basic(io->mod);
            Free_Input(io);
            io = PhyML_XML(optarg);
            Free(io);
            return 0;

#elif defined(PHYTIME)

            Free_Optimiz(io->mod->s_opt);
            Free_Model_Basic(io->mod);
            Free_Input(io);
            DATE_XML(optarg);
            return 0;

#elif defined(PHYREX)

            Free_Optimiz(io->mod->s_opt);
            Free_Model_Basic(io->mod);
            Free_Input(io);
            PHYREX_XML(optarg);
            return 0;
            
#endif
            break;
	  }
	case 72:
	  {
	    io->mod->s_opt->constrained_br_len = YES;
	    break;
	  }
	case 71:
	  {
	    io->mcmc->in_fp_par = fopen(optarg,"r");
	    io->mcmc->randomize = NO;
	    break;
	  }
	case 70:
	  {
	    io->mutmap = YES;
	    break;
	  }
	case 68:
	  {
	    char *tmp;
	    tmp = (char *)mCalloc(T_MAX_FILE, sizeof(char));

	    if(strlen(optarg) > T_MAX_FILE -11)
	      {
		char choix;
		strcpy (tmp, "\n. The file name'");
		strcat (tmp, optarg);
		strcat (tmp, "' is too long.\n");
		PhyML_Printf("%s",tmp);
		PhyML_Printf("\n. Type any key to exit.\n");
		if(!scanf("%c",&choix)) Exit("\n");
		Exit("\n");
	      }
	    else if (!Filexists (optarg))
	      {
		char choix;
		strcpy (tmp, "\n. The file '");
		strcat (tmp, optarg);
		strcat (tmp, "' doesn't exist.\n");
		PhyML_Printf("%s",tmp);
		PhyML_Printf("\n. Type any key to exit.\n");
		if(!scanf("%c",&choix)) Exit("\n");
		Exit("\n");
	      }
	    else
	      {
		strcpy(io->in_constraint_tree_file, optarg);
		io->fp_in_constraint_tree = Openfile(io->in_constraint_tree_file,0);
	      }
	    Free(tmp);
	    break;
	  }
	case 67:
	  {
	    phydbl pos;
	    pos = atof(optarg);
	    io->codpos = (int)pos;
	    if(io->codpos < 1 || io->codpos > 3)
	      {
		char choix;
		PhyML_Printf("\n. Coding position must be set to 1, 2 or 3.\n");
		PhyML_Printf("\n. Type any key to exit.\n");
		if(!scanf("%c",&choix)) Exit("\n");
		Exit("\n");
	      }
	    break;
	    
	    break;
	  }
	case 66:
	  {
	    io->mod->gamma_mgf_bl = YES;
	    io->mod->s_opt->opt_gamma_br_len = YES;
	    break;
	  }
	case 65:
	  {
	    io->mod->log_l = YES;
	    break;
	  }
	case 64:
	  {
	    char *s;
	    int i;
	    s = (char *)mCalloc(T_MAX_NAME,sizeof(char));
	    i = 0;
	    while(optarg[i++]) s[i]=tolower(optarg[i]);
	    if(!strcmp(optarg,"gbd")) io->rates->model_id = THORNE;
	    else if(!strcmp(optarg,"autocorrelated")) io->rates->model_id = THORNE;
	    else if(!strcmp(optarg,"autocorr")) io->rates->model_id = THORNE;
	    else if(!strcmp(optarg,"geom")) io->rates->model_id = GUINDON;
	    else if(!strcmp(optarg,"integrated")) io->rates->model_id = GUINDON;
	    else if(!strcmp(optarg,"uncorrelated")) io->rates->model_id = LOGNORMAL;
	    else if(!strcmp(optarg,"uncorr")) io->rates->model_id = LOGNORMAL;
	    else if(!strcmp(optarg,"lognorm")) io->rates->model_id = LOGNORMAL;
	    else if(!strcmp(optarg,"clock")) io->rates->model_id = STRICTCLOCK;
	    else if(!strcmp(optarg,"strictclock")) io->rates->model_id = STRICTCLOCK;
	    else if(!strcmp(optarg,"strict_clock")) io->rates->model_id = STRICTCLOCK;
	    else 
	      {
		PhyML_Printf("\n. rate_model should be 'autocorrelated', 'uncorrelated', 'integrated' or 'clock'.");
		Exit("\n");
	      }
	    Free(s);
	    break;
	  }
	
	
	case 62:
	  {
	    io->mcmc->is = YES;
	    break;
	  }
	case 61:
	  {
	    io->mod->ras->free_mixt_rates            = YES;
	    io->mod->s_opt->opt_free_mixt_rates = YES;
	    break;
	  }
	case 60:
	  {
	    io->lk_approx = NORMAL;		
	    break;
	  }
	case 59:
	  {
	    break;
	  }
	case 58:
	  {
	    io->do_alias_subpatt = YES;
	    break;
	  }
	case 57:
	  {	    
	    io->colalias = NO;
	    break;
	  }
	case 56:
	  {
	    io->mem_question = NO;
	    break;
	  }
	case 55:
	  {
	    phydbl len;
	    len = atof(optarg);
	    io->mcmc->chain_len_burnin = (int)len;
	    if(io->mcmc->chain_len_burnin < 1)
	      {
		char choix;
		PhyML_Printf("\n. chain_len_burnin must be an integer greater than 0.\n");
		PhyML_Printf("\n. Type any key to exit.\n");
		if(!scanf("%c",&choix)) Exit("\n");
		Exit("\n");
	      }
	    break;
	  }	  
	case 54:
	  {
	    phydbl len;
	    len = atof(optarg);
	    io->mcmc->sample_interval = (int)len;
	    if(io->mcmc->sample_interval < 1)
	      {
		char choix;
		PhyML_Printf("\n. sample_interval must be an integer greater than 0.\n");
		PhyML_Printf("\n. Type any key to exit.\n");
		if(!scanf("%c",&choix)) Exit("\n");
		Exit("\n");
	      }
	    break;
	  }	  
	case 53:
	  {
	    phydbl len;
	    len = atof(optarg);
	    io->mcmc->chain_len = (int)len;
	    if(io->mcmc->chain_len < 1)
	      {
		char choix;
		PhyML_Printf("\n. chain_len must be an integer greater than 0.\n");
		PhyML_Printf("\n. Type any key to exit.\n");
		if(!scanf("%c",&choix)) Exit("\n");
		Exit("\n");
	      }
	    break;
	  }	  
	case 52:
	  {
	    char *s;
	    s = (char *)mCalloc(T_MAX_FILE, sizeof(char));
	    strcpy(s,optarg);
	    io->mod->fp_aa_rate_mat = Openfile(s,0);
	    strcpy(io->mod->aa_rate_mat_file->s,s);
	    Free(s);
	    break;
	  }
	case 51:
	  {
	    io->boot_prog_every = atoi(optarg);
	    if(io->boot_prog_every < 1)
	      {
		char choix;
		PhyML_Printf("\n. boot_progress_every must be an integer greater than 0.\n");
		PhyML_Printf("\n. Type any key to exit.\n");
		if(!scanf("%c",&choix)) Exit("\n");
		Exit("\n");
	      }
	    break;
	  }
	case 50:
	  {
	    strcpy(io->clade_list_file,optarg);
	    break;
	  }
	case 49:
	  {
	    PhyML_Printf("\n. This is PhyML version %s.\n\n",VERSION);
	    Exit("");
	    break;
	  }	  
	case 48 : 
	  {
	    io->quiet = 1;
	    break;
	  }
	case 'p' : case 47 : 
	  {
	    io->in_tree = 1;
	    break;
	  }
	case 46 : 
	  {
            if(strlen(optarg) > 0)
              {
                io->append_run_ID = YES;
                strcpy(io->run_id_string,optarg);
              }
            break;
	  }
	case 45 : 
	  {
	    io->mod->ras->gamma_median = 1;
	    break;
	  }
	case 44 :
	  {
	    io->mod->s_opt->hybrid_thresh = 0;
	    break;
	  }
	case 43 :
	  {
	    io->mod->s_opt->min_diff_lk_move = atof(optarg);
	    if(io->mod->s_opt->min_diff_lk_move < 0)
	      {
		char choix;
		PhyML_Printf("\n. Min_diff_lk_move must be a double greater than 0.\n");
		PhyML_Printf("\n. Type any key to exit.\n");
		if(!scanf("%c",&choix)) Exit("\n");
		Exit("\n");
	      }
	    break;
	  }
	case 42 :
	  {
	    io->mod->s_opt->pars_thresh = (int)atoi(optarg);
	    if(io->mod->s_opt->pars_thresh < 0)
	      {
		PhyML_Printf("\n. The parsimony threshold must be an integer greater than 0.\n");
		PhyML_Printf("\n. Type any key to exit.\n");
		Exit("\n");
	      }
	    break;
	  }
	case 41 :
	  {
	    io->mod->s_opt->opt_five_branch = 0;
	    break;
	  }
	case 40 :
	  {
	    writemode = APPEND;
	    break;
	  }
	case 39 :
	  {
	    break;
	  }
	case 38 :
	  {
	    io->rm_ambigu = 1;
	    break;
	  }
	case 37 :
	  {
	    io->mod->s_opt->opt_cov_free_rates = YES;
#ifdef M4
	    io->mod->m4mod->use_cov_alpha      = NO;
	    io->mod->m4mod->use_cov_free       = YES;
#endif
	    break;
	  }
	case 36 :
	  {
#ifndef PHYML
            open_ps_file = 1;
#endif
	    break;
	  }
	case 35 :
	  {
#ifdef M4
	    io->mod->m4mod->n_h = (int)atoi(optarg);

	    if(io->mod->m4mod->n_h < 1)
	      {
		char choix;
		PhyML_Printf("\n. The number of classes must be greater than 0.\n");
		PhyML_Printf("\n. Type any key to exit.\n");
		if(!scanf("%c",&choix)) Exit("\n");
		Exit("\n");
	      }
#endif
	    break;
	  }
	case 34 :
	  {
#ifdef M4
	    io->mod->m4mod->use_cov_alpha = YES;
	    io->mod->m4mod->use_cov_free  = NO;
	    
	    if(!strcmp(optarg,"e") || !strcmp(optarg,"E") ||
	       !strcmp(optarg,"estimated") || !strcmp(optarg,"ESTIMATED"))
	      {
		io->mod->s_opt->opt_cov_alpha = YES;
		io->mod->m4mod->alpha         = 1.0;
	      }
	    else
	      {
		io->mod->m4mod->alpha = (phydbl)atof(optarg);
		
		if(io->mod->m4mod->alpha < 1.E-5)
		  {
		    char choix;
		    PhyML_Printf("\n. The value of alpha must be greater than 1.E-5.\n");
		    PhyML_Printf("\n. Type any key to exit.\n");
		    if(!scanf("%c",&choix)) Exit("\n");
		    Exit("\n");
		  }
	      }
#endif
	    break;
	  }
	case 33 :
	  {
#ifdef M4
	    if(!strcmp(optarg,"e") || !strcmp(optarg,"E") ||
	       !strcmp(optarg,"estimated") || !strcmp(optarg,"ESTIMATED"))
	      {
		io->mod->s_opt->opt_cov_delta = YES;
		io->mod->m4mod->delta         = 1.0;
	      }
	    else
	      {
		io->mod->m4mod->delta = (phydbl)atof(optarg);
		
		if(atof(optarg) < 1.E-10)
		  {
		    char choix;
		    PhyML_Printf("\n. The value of delta must be larger than 1.E-10.\n");
		    PhyML_Printf("\n. Type any key to exit.\n");
		    if(!scanf("%c",&choix)) Exit("\n");
		    Exit("\n");
		  }
	      }
#endif
	    break;
	  }
	case 32 :
	  {
	    io->mod->use_m4mod = YES;
	    break;
	  }
	case 31 :
	  {
	    io->print_site_lnl = YES;
	    break;
	  }
	case 30 :
	  {
	    io->print_trace = YES;
	    break;
	  }
	case 29 :
	  {
	    io->random_boot_seq_order = (int)atoi(optarg);
	    break;
	  }
	case 28 :
	  {
	    io->collapse_boot = (int)atoi(optarg);
	    break;
	  }
	case 27 :
	  {
	    io->r_seed = (int)atoi(optarg);
	    break;
	  }
	case 26 :
	  {
	    io->mod->s_opt->general_pars = YES;
	    break;
	  }
	case 25 :
	  {
	    io->mod->s_opt->fast_nni = YES;
	    break;
	  }
	case 24 :
	  {
	    io->mod->s_opt->p_moves_to_examine = (phydbl)atof(optarg);
	    break;
	  }
	case 23 :
	  {
	    io->mod->s_opt->wim_inside_opt = 1;
	    break;
	  }
	case 0 :
	  {
	    io->mod->s_opt->wim_n_rgrft = atoi(optarg);
	    break;
	  }
	case 1 :
	  {
	    io->mod->s_opt->wim_n_globl = atoi(optarg);
	    break;
	  }
	case 2 :
	  {
	    io->mod->s_opt->wim_max_dist = atoi(optarg);
	    break;
	  }
	case 3 :
	  {
	    io->mod->s_opt->wim_n_optim = atoi(optarg);
	    break;
	  }
	case 4 :
	  {
	    io->mod->s_opt->wim_n_best = atoi(optarg);
	    break;
	  }
	case 16 :
	  {
	    io->mod->s_opt->min_diff_lk_local = atof(optarg);
	    break;
	  }
	case 17 :
	  {
	    io->mod->s_opt->min_diff_lk_global = atof(optarg);
	    break;
	  }
	case 18 :
	  {
	    io->mod->s_opt->steph_spr = NO;
	    io->mod->s_opt->greedy    = YES;
	    break;
	  }
	case 19 :
	  {
	    io->mod->s_opt->brent_it_max = atoi(optarg);
	    break;
	  }
	case 20 :
	  {
	    io->mod->s_opt->random_input_tree = YES;
	    break;
	  }
	case 21 :
	  {
	    io->mod->s_opt->random_input_tree = YES;
	    io->mod->s_opt->n_rand_starts = atoi(optarg);
	    if(io->mod->s_opt->n_rand_starts < 1) Exit("\n. Number of random starting trees must be > 0.\n\n");
	  }
	case 's':case 6:
	  {
	    if((!strcmp(optarg,"spr")) || (!strcmp(optarg,"SPR")))
	      {
		io->mod->s_opt->topo_search = SPR_MOVE;
		io->mod->s_opt->greedy      = (io->mod->s_opt->steph_spr)?(0):(1);
	      }
	    else if((!strcmp(optarg,"nni")) || (!strcmp(optarg,"NNI")))
	      {
                PhyML_Printf("\n. The NNI option is deprecated. PhyML now uses a mix of SPRs and NNIs.");
		io->mod->s_opt->topo_search         = NNI_MOVE;
		io->mod->s_opt->random_input_tree   = 0;
	      }
	    else if((!strcmp(optarg,"best")) || (!strcmp(optarg,"BEST")))
	      {
                PhyML_Printf("\n. The BEST option is deprecated. PhyML now uses a mix of SPRs and NNIs.");
		io->mod->s_opt->topo_search = BEST_OF_NNI_AND_SPR;
		io->mod->s_opt->greedy      = (io->mod->s_opt->steph_spr)?(0):(1);
	      }
	    break;
	  }
	  
	case 'd':case 7:
	  {
	    if(!strcmp(optarg,"nt"))
	      {
		io->datatype        = NT;
		io->mod->ns         = 4;
#ifdef M4
                io->mod->m4mod->n_o = 4;
#endif
                
		if((io->mod->whichmodel == LG)        ||
		   (io->mod->whichmodel == WAG)       ||
		   (io->mod->whichmodel == DAYHOFF)   ||
		   (io->mod->whichmodel == JTT)       ||
		   (io->mod->whichmodel == BLOSUM62)  ||
		   (io->mod->whichmodel == MTREV)     ||
		   (io->mod->whichmodel == RTREV)     ||
		   (io->mod->whichmodel == CPREV)     ||
		   (io->mod->whichmodel == DCMUT)     ||
		   (io->mod->whichmodel == VT)        ||
		   (io->mod->whichmodel == MTMAM)     ||
		   (io->mod->whichmodel == MTART)     ||
		   (io->mod->whichmodel == HIVW)      ||
		   (io->mod->whichmodel == HIVB)      ||
		   (io->mod->whichmodel == AB)        ||
		   (io->mod->whichmodel == CUSTOMAA)
		   )
		  {
		    io->mod->whichmodel = HKY85;
		    strcpy(io->mod->modelname->s, "HKY85\0");
		  }
	      }
	    else if (!strcmp(optarg,"aa"))
	      {
		io->datatype              = AA;
		io->mod->s_opt->opt_kappa = NO;
		io->mod->ns               = 20;
#ifdef M4
		io->mod->m4mod->n_o       = 20;
#endif
                
		if(
		   (io->mod->whichmodel == JC69)   ||
		   (io->mod->whichmodel == K80)    ||
		   (io->mod->whichmodel == F81)    ||
		   (io->mod->whichmodel == HKY85)  ||
		   (io->mod->whichmodel == F84)    ||
		   (io->mod->whichmodel == TN93)   ||
		   (io->mod->whichmodel == GTR)    ||
		   (io->mod->whichmodel == CUSTOM)
		   )
		  {
		    io->mod->whichmodel = LG;
		    strcpy(io->mod->modelname->s, "LG\0");
		  }
	      }
	    else if ((!strcmp(optarg,"generic")) || (!strcmp(optarg,"gen")))
	      {
		io->datatype = GENERIC;
	      }
	    else
	      {
		char choix;
		PhyML_Printf("\n. Unknown argument to -d option: please use `nt' for DNA or `aa' for Amino-Acids\n");
		PhyML_Printf("\n. Type any key to exit.\n");
		if(!scanf("%c",&choix)) Exit("\n");
		Exit("\n");
	      }
	    
	    break;
	  }
	case 'm': case 5 :
	  {
            opt_m = 1;
	    if (!isalpha(optarg[0]))
              {
                if(strchr(optarg,',') == NULL)
                  {
                    strcpy(io->mod->custom_mod_string->s,optarg);
                    if (strlen(io->mod->custom_mod_string->s) != 6)
                      {
                        Warn_And_Exit("\n. The custom model string should be of length 6.\n");
                      }
                    
                    io->datatype              = NT;
                    io->mod->whichmodel       = CUSTOM;
                    strcpy(io->mod->modelname->s, "custom");
                    io->mod->s_opt->opt_kappa = NO;
                    io->mod->s_opt->opt_rr    = YES;
                  }
                else
                  {
                    phydbl v;
                    int n_rr;
                    const char *d = ",";
                    char *tok = strtok(optarg,d);
                    
                    io->datatype           = NT;
                    io->mod->ns            = 4;
                    io->mod->whichmodel    = GTR;
                    io->mod->s_opt->opt_rr = NO;                    

                    io->mod->r_mat = (t_rmat *)Make_Rmat(io->mod->ns);
                    Init_Rmat(io->mod->r_mat);
                    Make_Custom_Model(io->mod);

                    n_rr = 0;
                    while(tok && n_rr < 6)
                      {
                        v = strtod(tok,NULL);
                        if (v != 0)
                          {
                            io->mod->r_mat->rr->v[n_rr] = v;
                            io->mod->r_mat->rr_val->v[n_rr] = log(v);
                          }
                        else
                          {
                            PhyML_Printf("\n. Invalid relative rate parameter value: '%s'.\n", tok);
                            Exit("\n");
                          }
                        tok = strtok (NULL,d);
                        n_rr++;
                      }
                    assert(n_rr <= 6);
                  }                
              }         
            else
              {
                char *s = To_Upper_String(optarg);

                if (strcmp(s, "JC69") == 0)
                  {
                    io->datatype              = NT;
                    io->mod->ns               = 4;
                    io->mod->whichmodel       = JC69;
                  }
                else if(strcmp(s, "K80") == 0)
                  {
                    io->datatype              = NT;
                    io->mod->ns               = 4;
                    io->mod->whichmodel       = K80;
                  }
                else if(strcmp(s, "F81") == 0)
                  {
                    io->datatype              = NT;
                    io->mod->ns               = 4;
                    io->mod->whichmodel       = F81;
                  }
                else if (strcmp(s, "HKY85") == 0)
                  {
                    io->datatype              = NT;
                    io->mod->ns               = 4;
                    io->mod->whichmodel       = HKY85;
                  }
                else if(strcmp(s, "F84") == 0)
                  {
                    io->datatype              = NT;
                    io->mod->ns               = 4;
                    io->mod->whichmodel       = F84;
                  }
                else if (strcmp(s,"TN93") == 0)
                  {
                    io->datatype              = NT;
                    io->mod->ns               = 4;
                    io->mod->whichmodel       = TN93;
                  }
                else if(strncmp(s, "GTR", 3) == 0)
                  {
                    io->datatype              = NT;
                    io->mod->ns               = 4;
                    io->mod->whichmodel       = GTR;
                    io->mod->s_opt->opt_rr    = YES;                
                  }
                else if(strcmp(s, "DAYHOFF") == 0)
                  {
                    io->datatype              = AA;
                    io->mod->ns               = 20;
                    io->mod->whichmodel       = DAYHOFF;
                  }
                else if(strcmp(s, "JTT") == 0)
                  {
                    io->datatype              = AA;
                    io->mod->ns               = 20;
                    io->mod->whichmodel       = JTT;
                  }
                else if(strcmp(s, "MTREV") == 0)
                  {
                    io->datatype             = AA;
                    io->mod->ns              = 20;
                    io->mod->whichmodel      = MTREV;
                  }
                else if(strcmp(s, "LG") == 0)
                  {
                    io->datatype              = AA;
                    io->mod->ns               = 20;
                    io->mod->whichmodel       = LG;
                  }
                else if(strcmp(s, "WAG") == 0)
                  {
                    io->datatype              = AA;
                    io->mod->ns               = 20;
                    io->mod->whichmodel       = WAG;
                  }
                else if(strcmp(s, "DCMUT") == 0)
                  {
                    io->datatype              = AA;
                    io->mod->ns               = 20;
                    io->mod->whichmodel       = DCMUT;
                  }
                else if(strcmp(s, "RTREV") == 0)
                  {
                    io->datatype              = AA;
                    io->mod->ns               = 20;
                    io->mod->whichmodel       = RTREV;
                  }
                else if(strcmp(s, "CPREV") == 0)
                  {
                    io->datatype              = AA;
                    io->mod->ns               = 20;
                    io->mod->whichmodel       = CPREV;
                  }
                else if(strcmp(s, "VT") == 0)
                  {
                    io->datatype              = AA;
                    io->mod->ns               = 20;
                    io->mod->whichmodel       = VT;
                  }
                else if(strcmp(s, "BLOSUM62") == 0)
                  {
                    io->datatype              = AA;
                    io->mod->ns               = 20;
                    io->mod->whichmodel       = BLOSUM62;
                  }
                else if(strcmp(s, "MTMAM") == 0)
                  {
                    io->datatype              = AA;
                    io->mod->ns               = 20;
                    io->mod->whichmodel       = MTMAM;
                  }
                else if (strcmp(s,"MTART") == 0)
                  {
                    io->datatype              = AA;
                    io->mod->ns               = 20;
                    io->mod->whichmodel       = MTART;
                  }
                else if (strcmp(s,"HIVW") == 0)
                  {
                    io->datatype              = AA;
                    io->mod->ns               = 20;
                    io->mod->whichmodel       = HIVW;
                  }
                else if(strcmp(s, "HIVB") == 0)
                  {
                    io->datatype              = AA;
                    io->mod->ns               = 20;
                    io->mod->whichmodel       = HIVB;
                  }
                else if(strcmp(s, "AB") == 0)
                  {
                    io->datatype              = AA;
                    io->mod->ns               = 20;
                    io->mod->whichmodel       = AB;
                  }
                else if (strcmp(s, "CUSTOM") == 0)
                  {
                    io->datatype              = AA;
                    io->mod->ns               = 20;
                    io->mod->whichmodel       = CUSTOMAA;
                  }
                else if(strcmp(s, "FLU") == 0)
                  {
                    io->datatype              = AA;
                    io->mod->ns               = 20;
                    io->mod->whichmodel       = FLU;
                  }
                else
                  {
                    PhyML_Printf("\n. The model name is incorrect. Please see the documentation.\n");
                    Exit("\n");
                  }
                Free(s);
              }
          
            
            Set_Model_Name(io->mod);
	    
	    break;
	  }
	  
	case 'a':case 14 :
	  {
	    if ((strcmp (optarg, "e") == 0) ||
		(strcmp (optarg, "E") == 0) ||
		(strcmp (optarg, "estimated") == 0) ||
		(strcmp (optarg, "ESTIMATED") == 0))
	      {
		io->mod->s_opt->opt_alpha = YES;
	      }
	    else if (atof(optarg) < 1.E-10)
	      {
		char choix;
		PhyML_Printf("\n. Alpha must be > 1.E-10.\n");
		PhyML_Printf("\n. Type any key to exit.\n");
		if(!scanf("%c",&choix)) Exit("\n");
		Exit("\n");
	      }
	    else
	      {
		io->mod->ras->alpha->v = (phydbl)atof(optarg);
		io->mod->s_opt->opt_alpha  = NO;
	      }
	    break;
	  }
	case 'b':case 10:
	  {
	    if ((int)String_To_Dbl(optarg) < -5)
	      {
		char choix;
		PhyML_Printf("\n. Branch test value must be a positive integer for bootstrap, or between -1 and -4 for aLRT branch test\n");
		PhyML_Printf("\n. Type any key to exit.\n");
		if(!scanf("%c",&choix)) Exit("\n");
		Exit("\n");
	      }
	    else
	      {
		if((int)String_To_Dbl(optarg) > 0)
		  {
                    io->do_alrt           = NO;
		    io->ratio_test        = 0;
		    io->n_boot_replicates = (int)atoi(optarg);
		    io->print_boot_trees  = YES;

		    if(io->n_data_sets > 1)
		      {
			char choix;
			PhyML_Printf("\n. Bootstrap option is not allowed with multiple data sets\n");
			PhyML_Printf("\n. Type any key to exit.\n");
			if(!scanf("%c",&choix)) Exit("\n");
			Exit("\n");
		      }
		  }
		else if (atoi(optarg)==0)
		  {
		    io->do_alrt    = NO;
                    io->do_tbe     = NO;
                    io->do_boot    = NO;
		    io->ratio_test = 0;
		  }
		else
		  {
                    io->do_alrt = YES;
                    io->do_tbe  = NO;
                    io->do_boot = NO;
		    io->ratio_test = -(int)atoi(optarg);
		  }
	      }
	    break;
	  }
	case 'c':case 12:
	  {
	    if ((!atoi(optarg)) || (atoi(optarg) < 0))
	      {
		char choix;
		PhyML_Printf("\n. Unknown argument to -c option: the number of rate categories must be a positive integer\n");
		PhyML_Printf("\n. Type any key to exit.\n");
		if(!scanf("%c",&choix)) Exit("\n");		  
		Exit("\n");
	      }
	    else 
	      {
		io->mod->ras->n_catg = atoi(optarg);
		if(io->mod->ras->n_catg < 1) 
		  {
		    PhyML_Printf("\n. The number of rate categories must be a positive integer\n");
		    Exit("\n");
		  }
	      }
	    break;
	  }
	case 'f':
	  {
	    if(!strcmp(optarg,"e"))
	      {
	        if(io->datatype == NT)
		  io->mod->s_opt->opt_state_freq = NO;
		else if (io->datatype == AA)
		  io->mod->s_opt->opt_state_freq = YES;
		else
		  {
		    PhyML_Printf("\n. Please define the data type (nt or aa) before setting the -f option\n");
		    Exit("\n");
		  }
	      }
	    else if(!strcmp(optarg,"m"))
	      {
	        if (io->datatype == NT)
		  io->mod->s_opt->opt_state_freq = YES;
		else if (io->datatype == AA)
		  io->mod->s_opt->opt_state_freq = NO;
		else
		  {
		    PhyML_Printf("\n. Please define the data type (nt or aa) before setting the -f option\n");
		    Exit("\n");
		  }
	      }
	    else if(!isalpha(optarg[0]))
	      {
		phydbl sum;
		double val1,val2,val3,val4;
		
                io->mod->e_frq = (t_efrq *)Make_Efrq(4);
                Init_Efrq(NULL,io->mod->e_frq);

		io->mod->s_opt->opt_state_freq  = NO;
                io->mod->e_frq->user_state_freq = YES;

		sscanf(optarg,"%lf,%lf,%lf,%lf",&val1,&val2,&val3,&val4);
		io->mod->e_frq->user_b_freq->v[0] = (phydbl)val1;
		io->mod->e_frq->user_b_freq->v[1] = (phydbl)val2;
		io->mod->e_frq->user_b_freq->v[2] = (phydbl)val3;
		io->mod->e_frq->user_b_freq->v[3] = (phydbl)val4;
		
		sum =
		  (io->mod->e_frq->user_b_freq->v[0] +
		   io->mod->e_frq->user_b_freq->v[1] +
		   io->mod->e_frq->user_b_freq->v[2] +
		   io->mod->e_frq->user_b_freq->v[3]);
		
		io->mod->e_frq->user_b_freq->v[0] /= sum;
		io->mod->e_frq->user_b_freq->v[1] /= sum;
		io->mod->e_frq->user_b_freq->v[2] /= sum;
		io->mod->e_frq->user_b_freq->v[3] /= sum;
		
                
		if(io->mod->e_frq->user_b_freq->v[0] < .0 ||
		   io->mod->e_frq->user_b_freq->v[1] < .0 ||
		   io->mod->e_frq->user_b_freq->v[2] < .0 ||
		   io->mod->e_frq->user_b_freq->v[3] < .0 ||
		   io->mod->e_frq->user_b_freq->v[0] > 1. ||
		   io->mod->e_frq->user_b_freq->v[1] > 1. ||
		   io->mod->e_frq->user_b_freq->v[2] > 1. ||
		   io->mod->e_frq->user_b_freq->v[3] > 1.)
		  {
		    Warn_And_Exit("\n. Invalid base frequencies.\n");
		  }
	      }
	    break;
	  }
	  
	case 'h':case 69:
	  {
	    Usage();
	    break;
	  }
	  
	case 'i':case 9:
	  {
	    char *tmp;
	    tmp = (char *) mCalloc (T_MAX_FILE, sizeof(char));
	    if (strlen (optarg) > T_MAX_FILE -16)
	      {
		char choix;
		strcpy (tmp, "\n. The file name'");
		strcat (tmp, optarg);
		strcat (tmp, "' is too long.\n");
		PhyML_Printf("%s",tmp);
		PhyML_Printf("\n. Type any key to exit.\n");
		if(!scanf("%c",&choix)) Exit("\n");
		Exit("\n");
	      }
	    
	    else if (!Filexists (optarg))
	      {
		char choix;
		strcpy (tmp, "\n. The file '");
		strcat (tmp, optarg);
		strcat (tmp, "' does not exist.\n");
		PhyML_Printf("%s",tmp);
		PhyML_Printf("\n. Type any key to exit.\n");
		if(!scanf("%c",&choix)) Exit("\n");
		Exit("\n");
	      }
	    else
	      {
		strcpy(io->in_align_file, optarg);
		io->fp_in_align = Openfile(io->in_align_file,0);
                
		strcpy(io->out_file, optarg);
		strcpy(io->out_tree_file,optarg);
#ifdef PHYML
		strcat(io->out_tree_file,"_phyml_tree");
#elif M4
		strcat(io->out_tree_file,"_m4_tree");
#elif PHYREX
		strcat(io->out_tree_file,"_phyrex_tree");
#endif
                
		strcpy(io->out_stats_file,optarg);
#ifdef PHYML
		strcat(io->out_stats_file,"_phyml_stats");
#elif M4
		strcat(io->out_stats_file,"_m4_stats");
#elif PHYREX
		strcat(io->out_stats_file,"_phyrex_stats");
#endif


#ifdef PHYREX
		strcpy(io->out_summary_file,optarg);
		strcat(io->out_summary_file,"_phyrex_summary");
#endif


	      }
	    Free (tmp);
	    break;
	  }
	  
	case 't':case 11:
	  {
            if(opt_m == 0)
              {
                PhyML_Fprintf(stderr,"\n. Please use the -m option before -t in the command line.");
                Exit("\n");
              }
            
            if((io->mod->whichmodel != JC69) && (io->mod->whichmodel != F81) && (io->mod->whichmodel != GTR))
              {
                if ((strcmp(optarg, "e") == 0) ||
                    (strcmp(optarg, "E") == 0) ||
                    (strcmp(optarg, "estimated") == 0) ||
                    (strcmp(optarg, "ESTIMATED") == 0))
                  {
                    io->mod->kappa->v = 4.0;
                    io->mod->s_opt->opt_kappa = YES;
                    if(io->mod->whichmodel == TN93) io->mod->s_opt->opt_lambda = YES;
                  }
                else
                  {
                    io->mod->s_opt->opt_kappa  = NO;
                    io->mod->s_opt->opt_lambda = NO;

                    // Added the 2 TsTv ratios for TN93
                    // lambda is the ratio of both TsTv ratios
                    // kappa is the mean of both TsTv ratios
                    if (io->mod->whichmodel == TN93)
                      {
                        double TsTvPur, TsTvPyr;
                        TsTvPur = TsTvPyr = -1.;
                        if(!isalpha(optarg[0]))
                          {
                            sscanf(optarg,"%lf,%lf",&TsTvPur,&TsTvPyr);
                          }
                        else
                          {
                            PhyML_Fprintf(stderr,"\n. The TN93 model requires two ts/tv ratios.\n");
                            Exit("\n");
                          }
                        if ( (TsTvPur < .0) || (TsTvPyr < .0) )
                          {
                            PhyML_Fprintf(stderr,"\n. ts/tv for purines: %f",TsTvPur);
                            PhyML_Fprintf(stderr,"\n. ts/tv for pyrimidines: %f",TsTvPyr);
                            PhyML_Fprintf(stderr,"\n. The TN93 model requires two ts/tv ratios");
                            PhyML_Fprintf(stderr,"\n. The command-line option should look as follows: ... -t 3.5,4.3 ... ");
                            PhyML_Fprintf(stderr,"\n. The ts/tv ratio must be a positive number.\n");
                            Exit("\n");
                          }
                        io->mod->lambda->v = (phydbl)(TsTvPur / TsTvPyr);
                        io->mod->kappa->v = (phydbl)((TsTvPur + TsTvPyr)/2.);
                      }
                    else
                      {
                        // -- End TN93 rates for purines & pyrimidines
                        if (atof(optarg) < .0)
                          {
                            char choix;
                            PhyML_Printf("\n. The ts/tv ratio must be a positive number\n");
                            PhyML_Printf("\n. Type any key to exit.\n");
                            if(!scanf("%c",&choix)) Exit("\n");
                            Exit("\n");
                          }
                        else
                          {
                            io->mod->kappa->v = (phydbl)atof(optarg);
                          }
                      }
                  }
	      }
	    break;
	  }
	case 'n':case 8:
	  {
	    if ((!atoi(optarg)) || (atoi(optarg) < 0))
	      {
		char choix;
		PhyML_Printf("\n. The number of alignments must be a positive integer\n");
		PhyML_Printf("\n. Type any key to exit.\n");
		if(!scanf("%c",&choix)) Exit("\n");
		Exit("\n");
	      }
	    else io->n_data_sets = atoi (optarg);
	    break;
	  }
	case 'q':case 22:
	  {
	    io->interleaved = NO;
	    break;
	  }
	case 'u':case 15:
	  {
	    char *tmp;
	    tmp = (char *)mCalloc(T_MAX_FILE, sizeof(char));
	    if(strlen(optarg) > T_MAX_FILE -11)
	      {
		char choix;
		strcpy (tmp, "\n. The file name'");
		strcat (tmp, optarg);
		strcat (tmp, "' is too long.\n");
		PhyML_Printf("%s",tmp);
		PhyML_Printf("\n. Type any key to exit.\n");
		if(!scanf("%c",&choix)) Exit("\n");
		Exit("\n");
	      }
	    else if (! Filexists (optarg))
	      {
		char choix;
		strcpy (tmp, "\n. The file '");
		strcat (tmp, optarg);
		strcat (tmp, "' doesn't exist.\n");
		PhyML_Printf("%s",tmp);
		PhyML_Printf("\n. Type any key to exit.\n");
		if(!scanf("%c",&choix)) Exit("\n");
		Exit("\n");
	      }
	    else
	      {
		strcpy(io->in_tree_file, optarg);
		io->in_tree = 2;
		io->fp_in_tree = Openfile(io->in_tree_file,READ);
	      }
	    Free(tmp);
	    break;
	  }
	  
	case 'v':case 13:
	  {
	    if ((strcmp (optarg, "e") == 0) ||
		(strcmp (optarg, "E") == 0) ||
		(strcmp (optarg, "estimated") == 0) ||
		(strcmp (optarg, "ESTIMATED") == 0)) 
	      {
		io->mod->s_opt->opt_pinvar = YES;
		io->mod->ras->invar        = YES;
	      }
	    
	    else if ((atof(optarg) < 0.0) || (atof(optarg) > 1.0))
	      {
		char choix;
		PhyML_Printf("\n. The proportion of invariable site must be a number between 0.0 and 1.0\n");
		PhyML_Printf("\n. Type any key to exit.");
		if(!scanf("%c",&choix)) Exit("\n");
		Exit("\n");
	      }
	    else
	      {
		io->mod->ras->pinvar->v = (phydbl)atof(optarg);
		if (io->mod->ras->pinvar->v > 0.0+SMALL)
		  io->mod->ras->invar = 1;
		else
		  io->mod->ras->invar = 0;
		io->mod->s_opt->opt_pinvar = 0;
	      }
	    break;
	  }
	case 'o':
	  {
	    if(!strcmp(optarg,"tlr"))
	      {
		io->mod->s_opt->opt_topo        = YES;
		io->mod->s_opt->opt_bl          = YES;
		io->mod->s_opt->opt_subst_param = YES;
	      }
	    else if(!strcmp(optarg,"tl"))
	      {
		io->mod->s_opt->opt_topo        = YES;
		io->mod->s_opt->opt_bl          = YES;
		io->mod->s_opt->opt_subst_param = NO;
	      }
	    else if(!strcmp(optarg,"t"))
	      {
		Warn_And_Exit("\n. You can't optimize the topology without adjusting branch length too...\n");
	      }
	    else if(!strcmp(optarg,"lr"))
	      {
		io->mod->s_opt->opt_topo        = NO;
		io->mod->s_opt->opt_bl          = YES;
		io->mod->s_opt->opt_subst_param = YES;
	      }
	    else if(!strcmp(optarg,"l"))
	      {
		io->mod->s_opt->opt_topo        = NO;
		io->mod->s_opt->opt_bl          = YES;
		io->mod->s_opt->opt_subst_param = NO;
	      }
	    else if(!strcmp(optarg,"r"))
	      {
		io->mod->s_opt->opt_topo        = NO;
		io->mod->s_opt->opt_bl          = NO;
		io->mod->s_opt->opt_subst_param = YES;
	      }
	    else if(!strcmp(optarg,"none") || !strcmp(optarg,"n"))
	      {
		io->mod->s_opt->opt_topo        = NO;
		io->mod->s_opt->opt_bl          = NO;
		io->mod->s_opt->opt_subst_param = NO;
	      }
	    else
	      {
		char choix;
		PhyML_Printf ("\n. The optimization parameter must be 'tlr' or 'tl' or 'lr' or 'l' or 'r' or 'n'.");
		PhyML_Printf("\n. Type any key to exit.\n");
		if(!scanf("%c",&choix)) Exit("\n");
		Exit("\n");
	      }
	    break;
	  }
	
	case '?':
	  {      
	    Exit("\n");
	    break;
	  }
	  
	case -1:
	  {      
	    break;
	  }

	default:
	  {
	    Usage();
	    break;
	  }
	}
    }while(c != -1);

  
  /*   if((io->mod->whichmodel == K80) || (io->mod->whichmodel == JC69)) */
  /*     { */
  /*       if(io->mod->s_opt->opt_state_freq) */
  /* 	{ */
  /* 	  char c; */
  /* 	  PhyML_Printf("\n. WARNING: nucleotide frequencies must be set to 1/4 with this model.\n"); */
  /* 	  PhyML_Printf("\n. Type the enter key to resume the analysis.\n"); */
  /* 	  scanf("%c",&c); */
  /* 	} */
  /*       io->mod->s_opt->opt_state_freq = 0; */
  /*     } */
  
  
  if(io->mod->s_opt->constrained_br_len == YES)
    {
      io->mod->s_opt->opt_topo = NO;
      /* io->mod->s_opt->opt_bl   = NO; */
    }

  
  if(io->do_tbe == YES)
    {
      io->do_alrt = NO;
      io->do_boot = NO;
    }
  else
    {
      if(io->do_alrt == NO && io->n_boot_replicates > 0) io->do_boot = YES;        
    }
  
#ifndef PHYML
  if((open_ps_file) || (io->m4_model == YES))
    {
      strcpy(io->out_ps_file,io->in_align_file);
      strcat(io->out_ps_file, "_mc_tree.ps");
      io->fp_out_ps = Openfile(io->out_ps_file,WRITE);
    }
#endif 
  
  
  if(io->datatype == UNDEFINED) io->datatype = NT;
  
  if ((io->datatype == NT) && (io->mod->whichmodel > 10))
    {
      char choix;
      PhyML_Printf("\n. Err.: model incompatible with the data type. Please use JC69, K80, F81, HKY, F84, TN93 or GTR\n");
      PhyML_Printf("\n. Type any key to exit.\n");
      if(!scanf("%c",&choix)) Exit("\n");
      Warn_And_Exit("\n");
    }
  else if ((io->datatype == AA) && (io->mod->whichmodel < 11))
    {
      char choix;
      PhyML_Printf("\n. Err.: model incompatible with the data type. Please use LG, Dayhoff, JTT, MtREV, WAG, DCMut, RtREV, CpREV, VT, Blosum62, MtMam, MtArt, HIVw, HIVb or AB.\n");
      PhyML_Printf("\n. Type any key to exit.\n");
      if(!scanf("%c",&choix)) Exit("\n");
      Exit("\n");
    }
  
  if(io->mod->use_m4mod == NO)
    {
      io->mod->s_opt->opt_cov_delta      = 0;
      io->mod->s_opt->opt_cov_alpha      = 0;
      io->mod->s_opt->opt_cov_free_rates = 0;
    }
  
  if((io->mod->s_opt->opt_cov_free_rates) && (io->mod->s_opt->opt_cov_alpha))
    {
      io->mod->s_opt->opt_cov_free_rates = 1;
#ifdef M4
      io->mod->m4mod->use_cov_alpha      = 0;
      io->mod->m4mod->use_cov_free       = 1;
#endif
    }
  
  if(io->print_site_lnl && io->fp_in_align != NULL)
    {
      strcpy(io->out_lk_file,io->in_align_file);
      strcat(io->out_lk_file, "_phyml_lk");
      if(io->append_run_ID) { strcat(io->out_lk_file,"_"); strcat(io->out_lk_file,io->run_id_string); }
      strcat(io->out_lk_file, ".txt");
      io->fp_out_lk = Openfile(io->out_lk_file,1);
    }
  
  if(io->print_trace && io->fp_in_align != NULL)
    {
      strcpy(io->out_trace_file,io->in_align_file);
      strcat(io->out_trace_file,"_phyml_trace");
      if(io->append_run_ID) { strcat(io->out_trace_file,"_"); strcat(io->out_trace_file,io->run_id_string); }
      strcat(io->out_trace_file,".txt");
      io->fp_out_trace = Openfile(io->out_trace_file,WRITE);
    }

  if(io->print_json_trace && io->fp_in_align != NULL)
    {
      strcpy(io->out_json_trace_file,io->in_align_file);
      strcat(io->out_json_trace_file,"_phyml_trace");
      if(io->append_run_ID) { strcat(io->out_json_trace_file,"_"); strcat(io->out_json_trace_file,io->run_id_string); }
      strcat(io->out_json_trace_file,".json");
      io->fp_out_json_trace = Openfile(io->out_json_trace_file,READWRITE);
    }

  
  if(io->mod->s_opt->random_input_tree && io->fp_in_align != NULL)
    {
      strcpy(io->out_trees_file,io->in_align_file);
      strcat(io->out_trees_file,"_phyml_rand_trees");
      if(io->append_run_ID) { strcat(io->out_trees_file,"_"); strcat(io->out_trees_file,io->run_id_string); }
      strcat(io->out_trees_file,".txt");
      io->fp_out_trees = Openfile(io->out_trees_file,1);
    }
  
  if((io->print_boot_trees) && (io->do_boot == YES || io->do_tbe == YES) && (io->fp_in_align != NULL))
    {
      strcpy(io->out_boot_tree_file,io->in_align_file);
      strcat(io->out_boot_tree_file,"_phyml_boot_trees");
      if(io->append_run_ID) { strcat(io->out_boot_tree_file,"_"); strcat(io->out_boot_tree_file,io->run_id_string); }
      strcat(io->out_boot_tree_file,".txt");
      io->fp_out_boot_tree = Openfile(io->out_boot_tree_file,1);
      
      
      strcpy(io->out_boot_stats_file,io->in_align_file);
      strcat(io->out_boot_stats_file,"_phyml_boot_stats");
      if(io->append_run_ID) { strcat(io->out_boot_stats_file,"_"); strcat(io->out_boot_stats_file,io->run_id_string); }
      strcat(io->out_boot_stats_file,".txt");
      io->fp_out_boot_stats = Openfile(io->out_boot_stats_file,1);
    }
  
  if(io->append_run_ID && io->fp_in_align != NULL)
    {
      strcat(io->out_tree_file,"_");
      strcat(io->out_stats_file,"_");
      strcat(io->out_tree_file,io->run_id_string);
      strcat(io->out_stats_file,io->run_id_string);
    }

  if(io->fp_in_align != NULL)
    {
      strcat(io->out_tree_file,".txt");
      strcat(io->out_stats_file,".txt");
    }


#ifdef PHYREX
  strcat(io->out_summary_file,".txt");
#endif

  
  if(io->mod->ras->n_catg == 1) io->mod->s_opt->opt_alpha = 0;
  
  if(io->mod->s_opt->opt_subst_param == NO)
    {
      io->mod->s_opt->opt_alpha  = NO;
      io->mod->s_opt->opt_kappa  = NO;
      io->mod->s_opt->opt_lambda = NO;
      io->mod->s_opt->opt_pinvar = NO;
      io->mod->s_opt->opt_rr     = NO;	
    }
  
  if(io->mod->whichmodel != K80 && 
     io->mod->whichmodel != HKY85 && 
     io->mod->whichmodel != F84 &&
     io->mod->whichmodel != TN93)
    {
      io->mod->s_opt->opt_kappa = NO;
    }
  
  if(io->datatype == AA && io->mod->whichmodel == CUSTOMAA && !io->mod->fp_aa_rate_mat)
    {
      PhyML_Printf("\n. Custom model option with amino-acid requires you to specify a rate matrix file through the '--aa_rate_file' option.\n");
      Exit("\n");
    }
  
#if !defined(PHYTIME) 
  // Make sure you don't erase the input file...
  if(!strcmp(io->out_tree_file,io->in_align_file) ||
     !strcmp(io->out_stats_file,io->in_align_file)) 
    {
      PhyML_Fprintf(stderr,"\n. The alignment file '%s' does not seem to exist...",io->in_align_file);
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);    
    }
  
  if(io->fp_in_align != NULL)
    {
      io->fp_out_tree  = Openfile(io->out_tree_file,writemode);
      io->fp_out_stats = Openfile(io->out_stats_file,writemode);
    }
#endif

#if defined(PHYREX)
  if(io->fp_in_align != NULL) io->fp_out_summary = Openfile(io->out_summary_file,writemode);
#endif


    if(io->quiet == NO)
    {
      PhyML_Printf("\n\n. Command line: ");
      for(i=0;i<argc;i++) PhyML_Printf("%s ",argv[i]);
    }

  return 1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

