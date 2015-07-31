/*

PhyML:  a program that  computes maximum likelihood phyLOGenies from
DNA or AA homoLOGous sequences.

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

  PhyML_Printf("\n. command-line: ");
  For(i,argc) PhyML_Printf("%s ",argv[i]);


  if(argc == 1) Exit("\n. No argument was passed to the program. Please check the documentation. \n");

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
      {0,0,0,0}
    };

  io->datatype = UNDEFINED;

  #ifndef PHYML
  int open_ps_file = 0;
  #endif

  idx=-1;

    do
    {     
      c = getopt_long(argc,argv,"qi:d:m:b:n:t:f:zk:v:c:a:u:ho:s:x:g:l:ep",longopts,&idx);

      switch(c)
	{

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
#ifndef INVITEE
            Free_Optimiz(io->mod->s_opt);
            M4_Free_M4_Model(io->mod->m4mod);
            Free_Model_Basic(io->mod);
            Free_Input(io);
            PhyML_XML(optarg);
            return 0;
#endif

#ifdef INVITEE           
            Free_Optimiz(io->mod->s_opt);
            M4_Free_M4_Model(io->mod->m4mod);
            Free_Model_Basic(io->mod);
            Free_Input(io);
            PhyTime_XML(optarg);
            return 0;
#endif
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
	    if(!strcmp(optarg,"gbd")) io->rates->model               = THORNE;
	    else if(!strcmp(optarg,"gbs")) io->rates->model          = GUINDON;
	    else if(!strcmp(optarg,"gamma")) io->rates->model        = GAMMA;
	    else if(!strcmp(optarg,"clock")) io->rates->model        = STRICTCLOCK;
	    else if(!strcmp(optarg,"strictclock")) io->rates->model  = STRICTCLOCK;
	    else if(!strcmp(optarg,"strict_clock")) io->rates->model = STRICTCLOCK;
	    else 
	      {
		PhyML_Printf("\n. rate_model should be 'gbs', 'gbd', 'gamma' or 'clock'.");
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
	    io->mcmc->use_data = NO;
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
	    io->append_run_ID = YES;
	    strcpy(io->run_id_string,optarg);
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
	    writemode = 2;
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
	    io->mod->m4mod->use_cov_alpha      = NO;
	    io->mod->m4mod->use_cov_free       = YES;
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
	    io->mod->m4mod->n_h = (int)atoi(optarg);
	    
	    if(io->mod->m4mod->n_h < 1)
	      {
		char choix;
		PhyML_Printf("\n. The number of classes must be greater than 0.\n");
		PhyML_Printf("\n. Type any key to exit.\n");
		if(!scanf("%c",&choix)) Exit("\n");
		Exit("\n");
	      }
	    break;
	  }
	case 34 :
	  {
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
	    break;
	  }
	case 33 :
	  {
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
	    io->mod->s_opt->steph_spr = 0;
	    io->mod->s_opt->greedy    = 1;
	    break;
	  }
	case 19 :
	  {
	    io->mod->s_opt->brent_it_max = atoi(optarg);
	    break;
	  }
	case 20 :
	  {
	    io->mod->s_opt->random_input_tree = 1;
	    break;
	  }
	case 21 :
	  {
	    io->mod->s_opt->random_input_tree = 1;
	    io->mod->s_opt->n_rand_starts = atoi(optarg);
	    if(io->mod->s_opt->n_rand_starts < 1) Exit("\n== Number of random starting trees must be > 0.\n\n");
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
		io->mod->s_opt->topo_search         = NNI_MOVE;
		io->mod->s_opt->random_input_tree   = 0;
	      }
	    else if((!strcmp(optarg,"best")) || (!strcmp(optarg,"BEST")))
	      {
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
		io->mod->m4mod->n_o = 4;
		
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
		io->mod->m4mod->n_o       = 20;

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
	    int i;
	    
	    For(i,strlen(optarg)) Uppercase(optarg+i);
	    
	    if(!isalpha(optarg[0]))
	      {
		strcpy(io->mod->custom_mod_string->s,optarg);
		
		if(strlen(io->mod->custom_mod_string->s) != 6)
		  {
		    Warn_And_Exit("\n. The string should be of length 6.\n");
		  }
		else
		  {
		    /* Make_Custom_Model(io->mod); */
		    /* Translate_Custom_Mod_String(io->mod); */
		  }
		
		io->datatype              = NT;
		io->mod->whichmodel       = CUSTOM;
		strcpy(io->mod->modelname->s, "custom");
		io->mod->s_opt->opt_kappa = NO;
		io->mod->s_opt->opt_rr    = YES;
	      }
	    
	    else if (strcmp(optarg, "JC69") == 0)
	      {
		io->datatype              = NT;
		io->mod->whichmodel       = JC69;
	      }
	    else if(strcmp(optarg, "K80") == 0)
	      {
		io->datatype              = NT;
		io->mod->whichmodel       = K80;
	      }
	    else if(strcmp(optarg, "F81") == 0)
	      {
		io->datatype              = NT;
		io->mod->whichmodel       = F81;
	      }
	    else if (strcmp(optarg, "HKY85") == 0)
	      {
		io->datatype              = NT;
		io->mod->whichmodel       = HKY85;
	      }
	    else if(strcmp(optarg, "F84") == 0)
	      {
		io->datatype              = NT;
		io->mod->whichmodel       = F84;
	      }
	    else if (strcmp (optarg,"TN93") == 0)
	      {
		io->datatype              = NT;
		io->mod->whichmodel       = TN93;
	      }
	    else if(strcmp (optarg, "GTR") == 0)
	      {
		io->datatype              = NT;
		io->mod->whichmodel       = GTR;
	      }
	    else if(strcmp(optarg, "DAYHOFF") == 0)
	      {
		io->datatype              = AA;
		io->mod->whichmodel       = DAYHOFF;
	      }
	    else if(strcmp (optarg, "JTT") == 0)
	      {
		io->datatype              = AA;
		io->mod->whichmodel       = JTT;
	      }
	    else if(strcmp(optarg, "MTREV") == 0)
	      {
		io->datatype             = AA;
		io->mod->whichmodel      = MTREV;
	      }
	    else if(strcmp (optarg, "LG") == 0)
	      {
		io->datatype              = AA;
		io->mod->whichmodel       = LG;
	      }
	    else if(strcmp (optarg, "WAG") == 0)
	      {
		io->datatype              = AA;
		io->mod->whichmodel       = WAG;
	      }
	    else if(strcmp(optarg, "DCMUT") == 0)
	      {
		io->datatype              = AA;
		io->mod->whichmodel       = DCMUT;
	      }
	    else if(strcmp (optarg, "RTREV") == 0)
	      {
		io->datatype              = AA;
		io->mod->whichmodel       = RTREV;
	      }
	    else if(strcmp(optarg, "CPREV") == 0)
	      {
		io->datatype              = AA;
		io->mod->whichmodel       = CPREV;
	      }
	    else if(strcmp(optarg, "VT") == 0)
	      {
		io->datatype              = AA;
		io->mod->whichmodel       = VT;
	      }
	    else if(strcmp(optarg, "BLOSUM62") == 0)
	      {
		io->datatype              = AA;
		io->mod->whichmodel       = BLOSUM62;
	      }
	    else if(strcmp(optarg, "MTMAM") == 0)
	      {
		io->datatype              = AA;
		io->mod->whichmodel       = MTMAM;
	      }
	    else if (strcmp(optarg,"MTART") == 0)
	      {
		io->datatype              = AA;
		io->mod->whichmodel       = MTART;
	      }
	    else if (strcmp(optarg,"HIVW") == 0)
	      {
		io->datatype              = AA;
		io->mod->whichmodel       = HIVW;
	      }
	    else if(strcmp(optarg, "HIVB") == 0)
	      {
		io->datatype              = AA;
		io->mod->whichmodel       = HIVB;
	      }
	    else if(strcmp(optarg, "AB") == 0)
	      {
		io->datatype              = AA;
		io->mod->whichmodel       = AB;
	      }
	    else if (strcmp(optarg, "CUSTOM") == 0)
	      {
		io->datatype              = AA;
		io->mod->whichmodel       = CUSTOMAA;
	      }
	    else
	      {
		PhyML_Printf("\n. The model name is incorrect. Please see the documentation.\n");
		Exit("\n");
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
		io->mod->s_opt->opt_alpha  = 0;
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
		    io->ratio_test       = 0;
		    io->mod->bootstrap   = (int)atoi(optarg);
		    io->print_boot_trees = 1;
		    
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
		    io->mod->bootstrap = 0;
		    io->ratio_test     = 0;
		  }
		else
		  {
		    io->mod->bootstrap = 0;
		    io->ratio_test     = -(int)atoi(optarg);
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
		
		io->mod->s_opt->opt_state_freq  = 0;
		io->mod->s_opt->user_state_freq = 1;
		
		/* 		sscanf(optarg,"%lf,%lf,%lf,%lf", */
		/* 		       io->mod->user_b_freq, */
		/* 		       io->mod->user_b_freq+1, */
		/* 		       io->mod->user_b_freq+2, */
		/* 		       io->mod->user_b_freq+3); */
		sscanf(optarg,"%lf,%lf,%lf,%lf",&val1,&val2,&val3,&val4);
		io->mod->user_b_freq->v[0] = (phydbl)val1;
		io->mod->user_b_freq->v[1] = (phydbl)val2;
		io->mod->user_b_freq->v[2] = (phydbl)val3;
		io->mod->user_b_freq->v[3] = (phydbl)val4;
		
		sum =
		  (io->mod->user_b_freq->v[0] +
		   io->mod->user_b_freq->v[1] +
		   io->mod->user_b_freq->v[2] +
		   io->mod->user_b_freq->v[3]);
		
		io->mod->user_b_freq->v[0] /= sum;
		io->mod->user_b_freq->v[1] /= sum;
		io->mod->user_b_freq->v[2] /= sum;
		io->mod->user_b_freq->v[3] /= sum;
		
		
		if(io->mod->user_b_freq->v[0] < .0 ||
		   io->mod->user_b_freq->v[1] < .0 ||
		   io->mod->user_b_freq->v[2] < .0 ||
		   io->mod->user_b_freq->v[3] < .0 ||
		   io->mod->user_b_freq->v[0] > 1. ||
		   io->mod->user_b_freq->v[1] > 1. ||
		   io->mod->user_b_freq->v[2] > 1. ||
		   io->mod->user_b_freq->v[3] > 1.)
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
	    if ((io->mod->whichmodel != JC69) &&
		(io->mod->whichmodel != F81)  &&
		(io->mod->whichmodel != GTR))
	      {
		if ((strcmp(optarg, "e") == 0) ||
		    (strcmp(optarg, "E") == 0) ||
		    (strcmp(optarg, "estimated") == 0) ||
		    (strcmp(optarg, "ESTIMATED") == 0))
		  {
		    io->mod->kappa->v              = 4.0;
		    io->mod->s_opt->opt_kappa      = YES;
		    if (io->mod->whichmodel == TN93)
		      io->mod->s_opt->opt_lambda   = YES;
		  }
		else
		  {
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
			io->mod->s_opt->opt_kappa  = 0;
			io->mod->s_opt->opt_lambda = 0;
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
		PhyML_Printf ("\n. The optimization parameter must be 'tlr' or 'tl' or 'lr' or 'l' or 'r' or ''.");
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

#ifndef PHYML
  if((open_ps_file) || (io->m4_model == YES))
    {
      strcpy(io->out_ps_file,io->in_align_file);
      strcat(io->out_ps_file, "_mc_tree.ps");
      io->fp_out_ps = Openfile(io->out_ps_file,WRITE);
    }
#endif 
  
  
  if(io->datatype == UNDEFINED) io->datatype = NT;


  if((io->mod->s_opt->n_rand_starts)           && 
     (io->mod->s_opt->topo_search == NNI_MOVE) && 
     (io->mod->s_opt->random_input_tree))
    {
      Warn_And_Exit("\n== The random starting tree option is only compatible with SPR based search options.\n"); 
    }
  
  if ((io->datatype == NT) && (io->mod->whichmodel > 10))
    {
      char choix;
      PhyML_Printf("\n== Err.: model incompatible with the data type. Please use JC69, K80, F81, HKY, F84, TN93 or GTR\n");
      PhyML_Printf("\n== Type any key to exit.\n");
      if(!scanf("%c",&choix)) Exit("\n");
      Warn_And_Exit("\n");
    }
  else if ((io->datatype == AA) && (io->mod->whichmodel < 11))
    {
      char choix;
      PhyML_Printf("\n== Err.: model incompatible with the data type. Please use LG, Dayhoff, JTT, MtREV, WAG, DCMut, RtREV, CpREV, VT, Blosum62, MtMam, MtArt, HIVw, HIVb or AB.\n");
      PhyML_Printf("\n== Type any key to exit.\n");
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
      io->mod->m4mod->use_cov_alpha      = 0;
      io->mod->m4mod->use_cov_free       = 1;
    }
  
  if(io->print_site_lnl)
    {
      strcpy(io->out_lk_file,io->in_align_file);
      strcat(io->out_lk_file, "_phyml_lk");
      if(io->append_run_ID) { strcat(io->out_lk_file,"_"); strcat(io->out_lk_file,io->run_id_string); }
      io->fp_out_lk = Openfile(io->out_lk_file,1);
    }
  
  if(io->print_trace)
    {
      strcpy(io->out_trace_file,io->in_align_file);
      strcat(io->out_trace_file,"_phyml_trace");
      if(io->append_run_ID) { strcat(io->out_trace_file,"_"); strcat(io->out_trace_file,io->run_id_string); }
      io->fp_out_trace = Openfile(io->out_trace_file,1);
    }
  
  if(io->mod->s_opt->random_input_tree)
    {
      strcpy(io->out_trees_file,io->in_align_file);
      strcat(io->out_trees_file,"_phyml_rand_trees");
      if(io->append_run_ID) { strcat(io->out_trees_file,"_"); strcat(io->out_trees_file,io->run_id_string); }
      io->fp_out_trees = Openfile(io->out_trees_file,1);
    }
  
  if((io->print_boot_trees) && (io->mod->bootstrap > 0))
    {
      strcpy(io->out_boot_tree_file,io->in_align_file);
      strcat(io->out_boot_tree_file,"_phyml_boot_trees");
      if(io->append_run_ID) { strcat(io->out_boot_tree_file,"_"); strcat(io->out_boot_tree_file,io->run_id_string); }
      io->fp_out_boot_tree = Openfile(io->out_boot_tree_file,1);
      
      strcpy(io->out_boot_stats_file,io->in_align_file);
      strcat(io->out_boot_stats_file,"_phyml_boot_stats");
      if(io->append_run_ID) { strcat(io->out_boot_stats_file,"_"); strcat(io->out_boot_stats_file,io->run_id_string); }
      io->fp_out_boot_stats = Openfile(io->out_boot_stats_file,1);
    }
  
  if(io->append_run_ID)
    {
      strcat(io->out_tree_file,"_");
      strcat(io->out_stats_file,"_");
      strcat(io->out_tree_file,io->run_id_string);
      strcat(io->out_stats_file,io->run_id_string);
    }
  
  if(io->mod->ras->n_catg == 1) io->mod->s_opt->opt_alpha = 0;
  
  if(!io->mod->s_opt->opt_subst_param)
    {
      io->mod->s_opt->opt_alpha  = 0;
      io->mod->s_opt->opt_kappa  = 0;
      io->mod->s_opt->opt_lambda = 0;
      io->mod->s_opt->opt_pinvar = 0;
      io->mod->s_opt->opt_rr     = 0;	
    }
  
  if(io->mod->whichmodel != K80 && 
     io->mod->whichmodel != HKY85 && 
     io->mod->whichmodel != F84 &&
     io->mod->whichmodel != TN93)
    {
      io->mod->s_opt->opt_kappa = 0;
    }
  
  if(io->datatype == AA && io->mod->whichmodel == CUSTOMAA && !io->mod->fp_aa_rate_mat)
    {
      PhyML_Printf("\n== Custom model option with amino-acid requires you to specify a rate matrix file through the '--aa_rate_file' option.\n");
      Exit("\n");
    }
  
#if !defined(PHYTIME) 
  // Make sure you don't erase the input file...
  if(!strcmp(io->out_tree_file,io->in_align_file) ||
     !strcmp(io->out_stats_file,io->in_align_file)) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);    

  writemode = WRITE;

  io->fp_out_tree  = Openfile(io->out_tree_file,writemode);
  io->fp_out_stats = Openfile(io->out_stats_file,writemode);
#endif
  

#if defined(PHYREX)
  io->fp_out_summary = Openfile(io->out_summary_file,writemode);
#endif

  writemode++; // just to silence a warning message at compilation
  
  if(io->mod->whichmodel == GTR) 
    {
      /* Make_Custom_Model(io->mod); */
      io->mod->s_opt->opt_rr = 1;
    }
  

  return 1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

