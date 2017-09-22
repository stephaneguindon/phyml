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
  /* return(PHYREX_Main_Estimate(argc,argv)); */
  return(PHYREX_Main_Simulate(argc,argv));
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int PHYREX_Main_Estimate(int argc, char *argv[])
{
  t_tree *tree;
  phydbl *res;
  int n_dim,i;
  t_dsk *disk;
  option *io;
  calign *cdata;
  t_ldsk **ldsk_a;

  n_dim = 2;

  io = (option *)Get_Input(argc,argv);
  assert(io);

  Get_Seq(io);
  assert(io->data);

  Make_Model_Complete(io->mod);
  Set_Model_Name(io->mod);
  Print_Settings(io);

  cdata = Compact_Data(io->data,io);
  Free_Seq(io->data,cdata->n_otu);

  tree = Make_Tree_From_Scratch(cdata->n_otu,cdata);
  Connect_CSeqs_To_Nodes(cdata,io,tree);

  tree->rates = RATES_Make_Rate_Struct(tree->n_otu);
  RATES_Init_Rate_Struct(tree->rates,io->rates,tree->n_otu);
  
  /* Allocate migrep model */
  tree->mmod = PHYREX_Make_Migrep_Model(n_dim);
  PHYREX_Init_Migrep_Mod(tree->mmod,n_dim,10.0,10.0);

  tree->data      = cdata;
  tree->mod       = io->mod;
  tree->io        = io;
  tree->n_pattern = tree->data->crunch_len;

  /* Allocate and initialise first disk event */
  disk = PHYREX_Make_Disk_Event(n_dim,tree->n_otu);
  PHYREX_Init_Disk_Event(disk,n_dim,NULL);
  disk->time             = 0.0;
  disk->mmod             = tree->mmod;
  disk->n_ldsk_a         = tree->n_otu;  
  tree->disk             = disk;

  /* Allocate coordinates for all the tips first (will grow afterwards) */
  ldsk_a = (t_ldsk **)mCalloc(tree->n_otu,sizeof(t_ldsk *));
  for(i=0;i<tree->n_otu;i++) 
    {
      ldsk_a[i] = PHYREX_Make_Lindisk_Node(n_dim);
      PHYREX_Init_Lindisk_Node(ldsk_a[i],disk,n_dim);
    }
  
  PHYREX_Read_Tip_Coordinates(ldsk_a,tree);

  tree->disk->ldsk_a = ldsk_a;

  /* Initialize parameters of migrep model */
  tree->mmod->lbda  = Uni()*(0.3 - 0.05) + 0.05;
  tree->mmod->mu    = Uni()*(1.0 - 0.3)  + 0.3;
  tree->mmod->rad   = Uni()*(5.0 - 1.5)  + 1.5;
  tree->mmod->sigsq = PHYREX_Update_Sigsq(tree);

  /* Random genealogy */
  PHYREX_Simulate_Backward_Core(NO,tree->disk,tree);

  PHYREX_Ldsk_To_Tree(tree);  

  Update_Ancestors(tree->n_root,tree->n_root->v[2],tree);
  Update_Ancestors(tree->n_root,tree->n_root->v[1],tree);
  RATES_Fill_Lca_Table(tree);

  disk = tree->disk;
  while(disk->prev) disk = disk->prev;

  tree->rates->bl_from_rt = YES;
  tree->rates->clock_r    = 0.01 / FABS(disk->time);
  tree->rates->model      = STRICTCLOCK;
  RATES_Update_Cur_Bl(tree);

  Init_Model(tree->data,io->mod,io);
  Prepare_Tree_For_Lk(tree);

  res = PHYREX_MCMC(tree);

  Free(res);  

  return 0;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int PHYREX_Main_Simulate(int argc, char *argv[])
{
  t_tree *tree;
  phydbl *res;
  int seed,pid,i;
  char *s;
  t_dsk *disk;
  int n_sites,n_otus;

  s = (char *)mCalloc(T_MAX_NAME,sizeof(char));

  pid     = getpid();
  seed    = pid;
  n_otus  = (int)atoi(argv[1]);
  n_sites = (int)atoi(argv[2]);

  printf("\n. seed: %d",seed);
  /* seed = 32076; /\* !!!!!!!!!! *\/ */
  srand(seed);
  
  /* tree = PHYREX_Simulate(n_otus,n_sites,10.,10.,seed); */
  tree = PHYREX_Simulate_Independent_Loci(n_otus,500,20.,20.,seed);

  disk = tree->disk;
  while(disk->prev) disk = disk->prev;

  strcpy(s,"phyrex_trees");
  sprintf(s+strlen(s),".%d",tree->mod->io->r_seed);
  tree->io->fp_out_tree = Openfile(s,WRITE);
  strcpy(s,"phyrex_stats");
  sprintf(s+strlen(s),".%d",tree->mod->io->r_seed);
  tree->io->fp_out_stats = Openfile(s,WRITE);
  strcpy(s,"phyrex_summary");
  sprintf(s+strlen(s),".%d",tree->mod->io->r_seed);
  tree->io->fp_out_summary = Openfile(s,WRITE);
  strcpy(s,"phyrex_mtt");
  sprintf(s+strlen(s),".%d.xml",tree->mod->io->r_seed);
  PHYREX_Print_MultiTypeTree_Config_File(n_sites,s,tree);

  res = PHYREX_MCMC(tree);

  disk = tree->disk;
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
  MCMC_Free_MCMC(tree->mcmc);
  Free_Mmod(tree->mmod);
  Free_Spr_List(tree);
  Free_Triplet(tree->triplet_struct);
  Free_Tree_Pars(tree);
  Free_Tree_Lk(tree);
  Free_Input(tree->io);
  Free_Optimiz(tree->mod->s_opt);
  Free_Model_Complete(tree->mod);
  Free_Model_Basic(tree->mod);
  Free_Calign(tree->data);
  Free_Tree(tree);
  Free(res);  
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
  phydbl min_neigh, max_neigh;
  phydbl area, neigh;
  phydbl T;
  phydbl Ne,maxNe,minNe;
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
    
  Prepare_Tree_For_Lk(tree);



  // migrep model stuff */
  mmod = PHYREX_Make_Migrep_Model(n_dim);
  tree->mmod = mmod;
  PHYREX_Init_Migrep_Mod(mmod,n_dim,w,h);


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
  tree->disk = PHYREX_Make_Disk_Event(mmod->n_dim,tree->n_otu);
  PHYREX_Init_Disk_Event(tree->disk,mmod->n_dim,NULL);
  tree->disk->time             = 0.0;
  tree->disk->mmod             = mmod;
  tree->disk->centr->lonlat[0] = .5*mmod->lim->lonlat[0];
  tree->disk->centr->lonlat[1] = .5*mmod->lim->lonlat[1];      
  tree->disk->n_ldsk_a         = tree->n_otu;  
  
  for(i=0;i<tree->n_otu;++i) 
    {
      char *s;
      tree->disk->ldsk_a[i] = PHYREX_Make_Lindisk_Node(mmod->n_dim);
      PHYREX_Init_Lindisk_Node(tree->disk->ldsk_a[i],tree->disk,mmod->n_dim);
      s = (char *)mCalloc(strlen(tree->disk->ldsk_a[i]->coord->id)+1+20,sizeof(char));
      strcpy(s,tree->disk->ldsk_a[i]->coord->id);
      strcat(s,"_deme0\0");
      Free(tree->disk->ldsk_a[i]->coord->id);
      tree->disk->ldsk_a[i]->coord->id = s;
    }
      
  for(i=0;i<tree->n_otu;i++)
    {
      tree->disk->ldsk_a[i]->coord->lonlat[0] = Uni()*tree->mmod->lim->lonlat[0]; // longitude
      tree->disk->ldsk_a[i]->coord->lonlat[1] = Uni()*tree->mmod->lim->lonlat[1]; // latitude
    }



  // Generate sequences
  subst_rate = .0;
  locus_idx = 0;
  do
    {
      PHYREX_Simulate_Backward_Core(NO,tree->disk,tree);
  
      // Random selection of a pair. Print out physical distances
      // between tips and time of coalescence for this pair, at this
      // particular locus
      {
        int *permut = Permutate(tree->n_otu);
        t_ldsk *lin1 = tree->disk->ldsk_a[permut[0]];
        t_ldsk *lin2 = tree->disk->ldsk_a[permut[1]];
        phydbl dist =
          sqrt(pow(lin1->coord->lonlat[0]-lin2->coord->lonlat[0],2) +
               pow(lin1->coord->lonlat[1]-lin2->coord->lonlat[1],2));

        t_dsk *disk = tree->disk;
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

      T = PHYREX_Tree_Height(tree);

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
      
      t_dsk *disk = tree->disk->prev;
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
t_tree *PHYREX_Simulate(int n_otu, int n_sites, phydbl w, phydbl h, int r_seed)
{  
  t_tree *tree;
  int n_dim,i;
  t_phyrex_mod *mmod;
  t_dsk *disk;
  option *io;
  t_mod *mod;
  t_opt *s_opt;
  calign *cdata;
  /* phydbl max_mu, min_mu; */
  /* phydbl max_rad, min_rad; */
  /* phydbl min_rate, max_rate; */
  /* phydbl min_lbda, max_lbda; */
  phydbl min_neigh, max_neigh;
  /* phydbl min_sigsq, max_sigsq; */
  phydbl area, neigh;
  phydbl T;
  phydbl Ne,maxNe,minNe;

  n_dim = 2; // 2-dimensional landscape
  area  = w * h;

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
  io->init_len = 10; /* sequence length */

  io->data = Make_Empty_Alignment(io);

  Make_Model_Complete(io->mod);
  Set_Model_Name(io->mod);

  /* Print_Settings(io); */

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
  PHYREX_Init_Migrep_Mod(mmod,n_dim,w,h);


  do
    {
      /* Effective population size */
      minNe = 100.; maxNe = 5000.;
      Ne = Uni() * (maxNe - minNe) + minNe;
      
      /* Neighborhood size */
      max_neigh = 0.01*Ne; min_neigh = 0.001*Ne;
      neigh = Uni()*(max_neigh - min_neigh)  + min_neigh;
    }
  while(neigh < 2.0);
  
  /* Death parameter */
  mmod->mu = 2./neigh;

  /* Theta (radius) */
  mmod->rad = Uni()*(3.0 - 0.5) + 0.5;
  mmod->sigsq = neigh / (4.*PI*Ne/area);
  mmod->lbda = area * mmod->sigsq / (4.*PI*mmod->mu*pow(mmod->rad,4));
  mmod->rho = Uni()*(10.0 - 0.1) + 0.1;
  // Duration of a generation in number of events
  mmod->gen_cal_time = 1./(2.*mmod->mu*pow(mmod->rad,2)/pow(w*h,2)*
                           (sqrt(2.)*mmod->rad*exp(-.5*pow(h/mmod->rad,2)) + h*sqrt(PI)*erf(sqrt(2.)*h/(2.*mmod->rad)) - sqrt(2.)*mmod->rad)*
                           (sqrt(2.)*mmod->rad*exp(-.5*pow(w/mmod->rad,2)) + w*sqrt(PI)*erf(sqrt(2.)*w/(2.*mmod->rad)) - sqrt(2.)*mmod->rad));
  // Divide by rate of events per calendar time unit to get duration of gen. in calendar time unit
  mmod->gen_cal_time *= 1./mmod->lbda;
  
  /* /\* !!!!!!!!!!!!! *\/ */
  /* mmod->lbda  = 1.0; */
  /* mmod->mu    = 0.5; */
  /* mmod->rad   = 1.E+10; */
  /* neigh       = 2./mmod->mu; */
  /* mmod->sigsq = PHYREX_Update_Sigsq(tree); */

  PHYREX_Simulate_Backward_Core(YES,tree->disk,tree);
  /* mmod->samp_area = PHYREX_Simulate_Forward_Core(n_sites,tree); */
  
  PHYREX_Ldsk_To_Tree(tree);  

  Update_Ancestors(tree->n_root,tree->n_root->v[2],tree);
  Update_Ancestors(tree->n_root,tree->n_root->v[1],tree);
  RATES_Fill_Lca_Table(tree);

  /* min_rate = 1.E-5; */
  /* max_rate = 1.E-4; */

  T = PHYREX_Tree_Height(tree);

  tree->rates->bl_from_rt = YES;
  /* tree->rates->clock_r    = Uni()*(max_rate - min_rate) + min_rate; */
  tree->rates->clock_r    = 0.01/FABS(T);
  tree->rates->model      = STRICTCLOCK;

  RATES_Update_Cur_Bl(tree);

  Init_Model(cdata,mod,io);

  tree->mod->ras->n_catg   = 1;
  tree->mod->whichmodel    = HKY85;
  tree->mod->kappa->v      = 4.0;
    
  Prepare_Tree_For_Lk(tree);
  Evolve(tree->data,tree->mod,0,tree);


  if(tree->mod->s_opt->greedy) Init_Partial_Lk_Tips_Double(tree);
  else                         Init_Partial_Lk_Tips_Int(tree);
  Init_Partial_Lk_Loc(tree);

  disk = tree->disk->prev;
  while(disk->prev) disk = disk->prev;
  
  PhyML_Printf("\n. Useful parameters: lambda=%G; mu=%G; rad=%G; clockr=%G; sigsq=%G; neigh=%G; N=%G; rhoe=%G",
               mmod->lbda,
               mmod->mu,
               mmod->rad,
               tree->rates->clock_r,
               mmod->sigsq,
               neigh,
               area*neigh/(4*PI*mmod->sigsq),
               neigh/(4.*PI*mmod->sigsq));
  fflush(NULL);

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
                   tree->disk->ldsk_a[i]->nd->name,
                   tree->disk->ldsk_a[i]->coord->lonlat[0],                   
                   tree->disk->ldsk_a[i]->coord->lonlat[1]);
    }
  PhyML_Printf("\n");

  return(tree);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
// Simulate Etheridge-Barton model backwards in time, following n_otu lineages
// on a rectangle of dimension width x h
phydbl PHYREX_Simulate_Backward_Core(int new_loc, t_dsk *init_disk, t_tree *tree)
{
  t_dsk *disk,*new_disk;
  int i,j;
  phydbl prob_hit,u;
  t_phyrex_mod *mmod;

  mmod  = tree->mmod;

  if(new_loc == YES)
    {
      init_disk = PHYREX_Make_Disk_Event(mmod->n_dim,tree->n_otu);
      PHYREX_Init_Disk_Event(init_disk,mmod->n_dim,NULL);
      init_disk->time             = 0.0;
      init_disk->mmod             = mmod;
      init_disk->centr->lonlat[0] = .5*mmod->lim->lonlat[0];
      init_disk->centr->lonlat[1] = .5*mmod->lim->lonlat[1];      
      init_disk->n_ldsk_a         = tree->n_otu;  
      tree->disk                  = init_disk;

      for(i=0;i<tree->n_otu;++i) 
        {
          char *s;
          init_disk->ldsk_a[i] = PHYREX_Make_Lindisk_Node(mmod->n_dim);
          PHYREX_Init_Lindisk_Node(init_disk->ldsk_a[i],init_disk,mmod->n_dim);
          s = (char *)mCalloc(strlen(init_disk->ldsk_a[i]->coord->id)+1+20,sizeof(char));
          strcpy(s,init_disk->ldsk_a[i]->coord->id);
          strcat(s,"_deme0\0");
          Free(init_disk->ldsk_a[i]->coord->id);
          init_disk->ldsk_a[i]->coord->id = s;
        }
      
      /* PhyML_Printf("\n. WARNING: position of samples are not random."); */
      /* Generate coordinates for the tip nodes (uniform distribution on the rectangle) */
      for(i=0;i<tree->n_otu;i++)
        {
          init_disk->ldsk_a[i]->coord->lonlat[0] = Uni()*tree->mmod->lim->lonlat[0]; // longitude
          init_disk->ldsk_a[i]->coord->lonlat[1] = Uni()*tree->mmod->lim->lonlat[1]; // latitude
          /* init_disk->ldsk_a[i]->coord->lonlat[0] = (i/(int)SQRT(tree->n_otu)+1)*tree->mmod->lim->lonlat[0]/(SQRT(tree->n_otu)+1); // longitude */
          /* init_disk->ldsk_a[i]->coord->lonlat[1] = (i%(int)SQRT(tree->n_otu)+1)*tree->mmod->lim->lonlat[1]/(SQRT(tree->n_otu)+1); // latitude */
        }
    }
 
  disk = init_disk;
  tree->disk = init_disk;

  /* Create new disk */
  new_disk = PHYREX_Make_Disk_Event(mmod->n_dim,tree->n_otu);
  PHYREX_Init_Disk_Event(new_disk,mmod->n_dim,NULL);
  new_disk->time = disk->time + Rexp(mmod->lbda);
  new_disk->n_ldsk_a = tree->n_otu;
  for(i=0;i<new_disk->n_ldsk_a;++i) 
    {
      new_disk->ldsk_a[i] = PHYREX_Make_Lindisk_Node(mmod->n_dim);
      PHYREX_Init_Lindisk_Node(new_disk->ldsk_a[i],new_disk,mmod->n_dim);
      new_disk->ldsk_a[i] = disk->ldsk_a[i];
    }
  disk->prev = new_disk;
  new_disk->next = disk;

  /* Move to it */
  disk = disk->prev;


  do
    {
      /* Create new disk */
      new_disk = PHYREX_Make_Disk_Event(mmod->n_dim,1);
      Free(new_disk->ldsk_a);
      new_disk->ldsk_a = NULL;
      new_disk->n_ldsk_a = 0;
      PHYREX_Init_Disk_Event(new_disk,mmod->n_dim,NULL);
      disk->prev = new_disk;      
      new_disk->next = disk;


      /* Time of next event */
      new_disk->time = disk->time + Rexp(mmod->lbda);

      /* Coordinates of current event */
      disk->centr->lonlat[0] = Uni()*mmod->lim->lonlat[0];
      disk->centr->lonlat[1] = Uni()*mmod->lim->lonlat[1];      
      
      /* Parent ?*/
      disk->ldsk = PHYREX_Make_Lindisk_Node(mmod->n_dim);
      PHYREX_Init_Lindisk_Node(disk->ldsk,disk,mmod->n_dim);
            
      /* Its location */
      switch(mmod->name)
        {
        case PHYREX_UNIFORM: { PHYREX_Runif_Rectangle_Overlap(disk->ldsk,disk,mmod); break; }
        case PHYREX_NORMAL:  { PHYREX_Rnorm_Trunc(disk->ldsk,disk,mmod); break; }
        }

      /* Which lineages in disk->ldsk_a are hit? -> populate new_disk->ldsk_a */
      for(i=0;i<disk->n_ldsk_a;++i)
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
                for(j=0;j<mmod->n_dim;++j) prob_hit += -POW(disk->ldsk_a[i]->coord->lonlat[j] - disk->centr->lonlat[j],2)/(2.*POW(mmod->rad,2));
                prob_hit = exp(prob_hit);
                break; 
              }
            }
          
          if(PHYREX_Is_In_Disk(disk->ldsk_a[i]->coord,disk,mmod) == YES)
            {
              u = Uni();
              if(u > prob_hit) // disk->ldsk_a[i] is not hit -> becomes a member of new_disk->ldsk_a
                {
                  if(new_disk->n_ldsk_a == 0) new_disk->ldsk_a = (t_ldsk **)mCalloc(1,sizeof(t_ldsk *));
                  else new_disk->ldsk_a = (t_ldsk **)mRealloc(new_disk->ldsk_a,new_disk->n_ldsk_a+1,sizeof(t_ldsk *));
                  new_disk->ldsk_a[new_disk->n_ldsk_a] = disk->ldsk_a[i];
                  new_disk->n_ldsk_a++;
                }
              else
                {
                  // disk->ldsk_a[i] is hit -> coalesce (or just jump) to parent (i.e., disk->ldsk)
                  disk->ldsk_a[i]->prev = disk->ldsk;

                  PHYREX_Make_Lindisk_Next(disk->ldsk);
                  disk->ldsk->next[disk->ldsk->n_next-1] = disk->ldsk_a[i];
                }
            }
                 
        }
          
      if(disk->n_ldsk_a == new_disk->n_ldsk_a) // No hit
        {          
          Free_Ldisk(disk->ldsk); // Free parent 
          disk->ldsk = NULL;
        }
      else
        {
          new_disk->ldsk_a = (t_ldsk **)mRealloc(new_disk->ldsk_a,new_disk->n_ldsk_a+1,sizeof(t_ldsk *));
          new_disk->ldsk_a[new_disk->n_ldsk_a] = disk->ldsk;
          new_disk->n_ldsk_a++;
        }
      
      if(new_disk->n_ldsk_a == 1) break;
      
      disk = new_disk;
    }
  while(1);

  disk->prev = NULL;


  /* do */
  /*   { */
  /*     PhyML_Printf("\n. disk %s n_next: %d",disk->id,disk->ldsk?disk->ldsk->n_next:-1); */
  /*     disk = disk->next; */
  /*   } */
  /* while(disk); */


  /* PHYREX_Print_Struct('#',tree); */
  /* Exit("\n"); */

  return(-1.);
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
  w             = mmod->lim->lonlat[0];
  h             = mmod->lim->lonlat[1];
  init_pop_size = mmod->rho*w*h;
  n_gen         = 10; // number of generations to simulate
  one_gen       = mmod->gen_cal_time; // time duration of one generation (calendar unit)
  t_sim         = n_gen * one_gen; // total simulation duration (calendar unit)



  printf("\n. one_gen: %f",one_gen);
  
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
      new_disk->centr->lonlat[0] = Uni()*mmod->lim->lonlat[0];
      new_disk->centr->lonlat[1] = Uni()*mmod->lim->lonlat[1];

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
  /* {     */
  /*   phydbl T = disk->time; // total simulation time (in calendar unit) */
  /*   t_ldsk *dum_ldsk = disk->ldsk_a[Rand_Int(0,disk->n_ldsk_a-1)]; */
  /*   t_dsk *dum_dsk = disk; */
  /*   phydbl s = w*h; // area */
  /*   printf("\n. s: %f",s); fflush(NULL); */
  /*   phydbl gentime = 1./(2.*PI*mmod->rad*mmod->rad*mmod->mu*mmod->lbda/s); */
  /*   phydbl ssq = 0.0; */
  /*   phydbl curr_pos,prev_pos; */
  /*   int curr_gen,prev_gen,nhits=0; */
  /*   prev_pos = dum_ldsk->coord->lonlat[0]; */
  /*   curr_pos = prev_pos; */
  /*   curr_gen = prev_gen = 1; */
  /*   do */
  /*     {       */
  /*       curr_gen  = 1 + (int)(dum_dsk->time-T)/gentime; */
  /*       curr_pos  = dum_ldsk->coord->lonlat[0]; */
        
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
              poly[i]->poly_vert[j]->lonlat[0] *= mmod->lim->lonlat[0]*0.5;
              poly[i]->poly_vert[j]->lonlat[1] *= mmod->lim->lonlat[1]*0.5;
            }
          
          max_x = 0.0;
          max_y = 0.0;
          for(j=0;j<poly[i]->n_poly_vert;++j)
            {
              if(poly[i]->poly_vert[j]->lonlat[0] > max_x) max_x = poly[i]->poly_vert[j]->lonlat[0];
              if(poly[i]->poly_vert[j]->lonlat[1] > max_y) max_y = poly[i]->poly_vert[j]->lonlat[1];
            }
          
          trans_x = Uni()*(mmod->lim->lonlat[0] - max_x);
          trans_y = Uni()*(mmod->lim->lonlat[1] - max_y);
          
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

  tree->disk             = disk;
  disk->ldsk_a           = ldsk_a_tips;
  disk->mmod             = tree->mmod;
  disk->centr->lonlat[0] = .5*tree->mmod->lim->lonlat[0];
  disk->centr->lonlat[1] = .5*tree->mmod->lim->lonlat[1];      
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

  
  /* for(i=0;i<n_otu;i++) printf("\n> %s",tree->disk->ldsk_a[i]->coord->id); */

  disk->prev = NULL;

  disk = tree->disk;
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
          if(FABS(coord->lonlat[i] - disk->centr->lonlat[i]) > disk->mmod->rad + 1.E-20)
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
  phydbl lnL;
  phydbl log_lbda;
  t_dsk *disk;
  
  assert(!tree->disk->next);
  assert(tree->disk->prev);
  
  tree->mmod->c_lnL = 0.0;
  log_lbda          = log(tree->mmod->lbda);

  /* TO DO: create a proper PHYREX_LogPost() function */
  tree->mmod->c_lnL += PHYREX_LnPrior_Radius(tree);
  tree->mmod->c_lnL += PHYREX_LnPrior_Mu(tree);
  tree->mmod->c_lnL += PHYREX_LnPrior_Lbda(tree);
 
  if(isinf(tree->mmod->c_lnL) || isnan(tree->mmod->c_lnL)) 
    {
      tree->mmod->c_lnL = UNLIKELY;
      return tree->mmod->c_lnL;
    }

  PHYREX_Update_Lindisk_List(tree);

  disk = tree->disk->prev;
  assert(disk);
  do
    {
      lnL = PHYREX_Lk_Core(disk,tree);
      lnL += log_lbda - tree->mmod->lbda * FABS(disk->time - disk->next->time);
      tree->mmod->c_lnL += lnL;
      disk->c_lnL = tree->mmod->c_lnL;
      disk = disk->prev;
    }
  while(disk);

  if(isinf(tree->mmod->c_lnL) || isnan(tree->mmod->c_lnL)) tree->mmod->c_lnL = UNLIKELY;

  return(tree->mmod->c_lnL);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl PHYREX_Lk_Core(t_dsk *disk, t_tree *tree)
{
  phydbl lnL,log_prob_hit,log_mu,log_dens_coal;
  int was_hit,i,j,k,err;
  phydbl two_theta_two;

  two_theta_two = 2.*POW(tree->mmod->rad,2);
  lnL           = 0.0;
  log_mu        = log(tree->mmod->mu);
  was_hit       = (disk->ldsk != NULL);

  for(i=0;i<disk->n_ldsk_a;i++)
    {
      if(PHYREX_Is_In_Ldscape(disk->ldsk_a[i],tree->mmod) == NO) return(UNLIKELY);     

      if(was_hit && disk->ldsk_a[i] == disk->ldsk)
        {
          for(k=0;k<disk->ldsk->n_next;k++)
            {              
              log_prob_hit = log_mu;
              for(j=0;j<tree->mmod->n_dim;j++)
                log_prob_hit += -POW(disk->ldsk->next[k]->coord->lonlat[j] - disk->centr->lonlat[j],2)/two_theta_two;

              lnL += log_prob_hit;
            }
        }
      else
        {
          log_prob_hit = log_mu;
          for(j=0;j<tree->mmod->n_dim;j++)
            log_prob_hit += -POW(disk->ldsk_a[i]->coord->lonlat[j] - disk->centr->lonlat[j],2)/two_theta_two;
          
          lnL += log(1. - exp(log_prob_hit));
        }
    }

  /* a hit occurred */
  if(was_hit == TRUE)
    {
      err = NO;
      log_dens_coal = 0.0;
      for(j=0;j<tree->mmod->n_dim;j++) log_dens_coal += Log_Dnorm_Trunc(disk->ldsk->coord->lonlat[j],
                                                                disk->centr->lonlat[j],
                                                                tree->mmod->rad,
                                                                0.0,
                                                                tree->mmod->lim->lonlat[j],&err);
      lnL += log_dens_coal;
    }

  /* Likelihood for the disk center */
  for(i=0;i<disk->mmod->n_dim;i++) lnL -= log(tree->mmod->lim->lonlat[i]);
  
  return(lnL);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl PHYREX_Lk_Range(t_dsk *young, t_dsk *old, t_tree *tree)
{
  t_dsk *disk;
  phydbl lnL,log_lbda;
  
  assert(young);
  assert(young->next);
  
  log_lbda = log(tree->mmod->lbda);

  lnL  = 0.0;
  disk = young;
  do
    {
      PHYREX_Update_Lindisk_List_Core(disk,tree);
      lnL += PHYREX_Lk_Core(disk,tree);
      lnL += log_lbda - tree->mmod->lbda * FABS(disk->time - disk->next->time);
      if(disk == old) break;
      disk = disk->prev;
    }
  while(disk);

  return(lnL);
}
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl *PHYREX_MCMC(t_tree *tree)
{
  t_mcmc *mcmc;
  int move,i,n_vars,burnin,true_ncoal,true_nint,true_nhits,n_demes;
  phydbl u;
  t_dsk *disk;
  FILE *fp_tree,*fp_stats,*fp_summary;
  phydbl *res;
  phydbl true_root_x, true_root_y,true_lbda,true_mu,true_sigsq,true_neigh,fst_neigh,diversity,true_rad,true_height,true_rhoe,tot_samp_area;
  int adjust_len;

  fp_tree    = tree->io->fp_out_tree;
  fp_stats   = tree->io->fp_out_stats;
  fp_summary = tree->io->fp_out_summary;

  mcmc = MCMC_Make_MCMC_Struct();

  disk = tree->disk;
  while(disk->prev) disk = disk->prev;

  tree->mcmc = mcmc;

  mcmc->io               = NULL;
  mcmc->is               = NO;
  mcmc->run              = 0;
  mcmc->chain_len_burnin = 1E+5;
  mcmc->randomize        = YES;
  mcmc->norm_freq        = 1E+3;
  mcmc->max_tune         = 1.E+20;
  mcmc->min_tune         = 1.E-10;
  mcmc->print_every      = 2;
  mcmc->is_burnin        = NO;
  mcmc->nd_t_digits      = 1;
  mcmc->chain_len        = 1E+8;
  mcmc->sample_interval  = 1E+3;  
  mcmc->max_lag          = 1000;
  mcmc->sample_size      = mcmc->chain_len/mcmc->sample_interval;
  mcmc->sample_num       = 0;
  adjust_len             = 1E+6;

  
  tot_samp_area = 0.0;
  /* if(tree->mmod->samp_area != NULL) */
  /*   for(i=0;i<tree->mmod->samp_area->n_poly;i++) tot_samp_area += Area_Of_Poly_Monte_Carlo(tree->mmod->samp_area->a_poly[i],tree->mmod->lim); */


  MCMC_Complete_MCMC(mcmc,tree);

  n_vars                 = 12;
  true_root_x            = disk->ldsk->coord->lonlat[0];
  true_root_y            = disk->ldsk->coord->lonlat[1];

  res = (phydbl *)mCalloc(tree->mcmc->chain_len / tree->mcmc->sample_interval * n_vars,sizeof(phydbl));

  PHYREX_Lk(tree);
  Lk(NULL,tree);

  true_lbda   = tree->mmod->lbda;
  true_mu     = tree->mmod->mu;
  true_sigsq  = tree->mmod->sigsq;
  true_rad    = tree->mmod->rad;
  true_neigh  = PHYREX_Neighborhood_Size(tree);
  fst_neigh   = PHYREX_Neighborhood_Size_Regression(tree);
  diversity   = Nucleotide_Diversity(tree->data);
  true_ncoal  = PHYREX_Total_Number_Of_Coal_Disks(tree);
  true_nint   = PHYREX_Total_Number_Of_Intervals(tree);
  true_nhits  = PHYREX_Total_Number_Of_Hit_Disks(tree);
  true_height = PHYREX_Tree_Height(tree);
  true_rhoe   = PHYREX_Effective_Density(tree);
  n_demes     = 0;
  /* n_demes     = tree->mmod->samp_area->n_poly; */
  
  PhyML_Fprintf(fp_stats,"\n# before rand glnL: %f alnL: %f",tree->mmod->c_lnL,tree->c_lnL);
  PhyML_Fprintf(fp_stats,"\n# ninter: %d",PHYREX_Total_Number_Of_Intervals(tree));
  PhyML_Fprintf(fp_stats,"\n# ncoal: %d",PHYREX_Total_Number_Of_Coal_Disks(tree));
  PhyML_Fprintf(fp_stats,"\n# nhits: %d",PHYREX_Total_Number_Of_Hit_Disks(tree));
  PhyML_Fprintf(fp_stats,"\n# root pos: %f %f",true_root_x,true_root_y);
  PhyML_Fprintf(fp_stats,"\n# root time: %f",disk->time);
  PhyML_Fprintf(fp_stats,"\n# true lbda: %f",tree->mmod->lbda);
  PhyML_Fprintf(fp_stats,"\n# true mu: %f",tree->mmod->mu);
  PhyML_Fprintf(fp_stats,"\n# true rad: %f",PHYREX_Update_Radius(tree));
  PhyML_Fprintf(fp_stats,"\n# true sigsq: %f",tree->mmod->sigsq);
  PhyML_Fprintf(fp_stats,"\n# true neigh. size: %f",PHYREX_Neighborhood_Size(tree));
  PhyML_Fprintf(fp_stats,"\n# fst-based estimate of neighborhood size: %f",PHYREX_Neighborhood_Size_Regression(tree));
  PhyML_Fprintf(fp_stats,"\n# true pop. density: %f",true_rhoe);
  PhyML_Fprintf(fp_stats,"\n# nucleotide diversity: %f",Nucleotide_Diversity(tree->data));
  PhyML_Fprintf(fp_stats,"\n# length of a generation: %G time units",PHYREX_Generation_Length(tree));
  PhyML_Fprintf(fp_stats,"\n# clock rate: %G subst. per time unit",tree->rates->clock_r);
  /* PhyML_Fprintf(fp_stats,"\n# of sampled demes: %d",n_demes); */
  /* if(tree->mmod->samp_area != NULL) */
  /*   for(i=0;i<tree->mmod->samp_area->n_poly;i++) PhyML_Fprintf(fp_stats,"\n# area of deme%d: %f", */
  /*                                                      i, */
  /*                                                      Area_Of_Poly_Monte_Carlo(tree->mmod->samp_area->a_poly[i],tree->mmod->lim)); */
 
  /* Starting parameter values */
  tree->mmod->lbda = Uni()*(0.5 - 0.2) + 0.2;
  tree->mmod->mu   = Uni()*(0.6 - 0.3) + 0.3;
  tree->mmod->rad  = Uni()*(4.0 - 2.0) + 2.0;
  PHYREX_Update_Sigsq(tree);

  /* tree->mmod->lbda = Uni()*(0.50 - 0.20) + 0.20; */
  /* tree->mmod->mu   = Uni()*(0.30 - 0.05) + 0.05; */
  /* tree->mmod->rad  = Uni()*(3.00 - 1.00) + 1.00; */
  /* PHYREX_Update_Sigsq(tree); */

  /* MCMC_Randomize_Rate_Across_Sites(tree); */
  MCMC_Randomize_Kappa(tree);

  /* Random genealogy */
  PHYREX_Simulate_Backward_Core(NO,tree->disk,tree);

  PHYREX_Lk(tree);

  Switch_Eigen(YES,tree->mod);
  Lk(NULL,tree);
  Switch_Eigen(NO,tree->mod);


  disk = tree->disk;
  while(disk->prev) disk = disk->prev;

  PhyML_Fprintf(fp_stats,"\n# after rand glnL: %f alnL: %f",tree->mmod->c_lnL,tree->c_lnL);
  PhyML_Fprintf(fp_stats,"\n# ninter: %d",PHYREX_Total_Number_Of_Intervals(tree));
  PhyML_Fprintf(fp_stats,"\n# ncoal: %d",PHYREX_Total_Number_Of_Coal_Disks(tree));
  PhyML_Fprintf(fp_stats,"\n# nhits: %d",PHYREX_Total_Number_Of_Hit_Disks(tree));
  PhyML_Fprintf(fp_stats,"\n# root pos: %f %f",true_root_x,true_root_y);
  PhyML_Fprintf(fp_stats,"\n# root time: %f",disk->time);
  PhyML_Fprintf(fp_stats,"\n# start lbda: %f",tree->mmod->lbda);
  PhyML_Fprintf(fp_stats,"\n# start mu: %f",tree->mmod->mu);
  PhyML_Fprintf(fp_stats,"\n# start rad: %f",tree->mmod->rad);
  fflush(NULL);


  PhyML_Fprintf(fp_stats,"\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s",
                "sample",
                "lnP",
                "alnL",
                "glnL",
                "lbda",
                "modLbda",
                "mu",
                "modNeigh",
                "sigsq",
                "modSigsq",
                "rad",
                "neigh",
                "rhoe",
                "nInt",
                "nCoal",
                "nHit",
                "rootTime",
                "tstv",
                "alpha",
                "accLbda",
                "accMu",
                "accRad",
                "accInDelDisk",
                "accInDelHit",
                "accScaleTime",
                "accSPR",
                "accPath",
                "accSim",
                "accLbdaTimes",
                "accSimPlus",
                "accIndelSerial",

                "accLdskGivenDisk",
                "accDiskGivenLdsk",
                "accDiskAndLdsk",
                "accLdskMulti",
                "accDiskMulti",

                "tuneLbda",
                "tuneRad",
                "tuneMu");

  for(i=0;i<mcmc->n_moves;i++) tree->mcmc->start_ess[i] = YES;

  Set_Both_Sides(NO,tree);
  mcmc->always_yes = NO;
  move             = -1;
  do
    {

      /* tree->mcmc->adjust_tuning[i] = NO; */
      if(mcmc->run > adjust_len) for(i=0;i<mcmc->n_moves;i++) tree->mcmc->adjust_tuning[i] = NO;

      if(tree->mmod->c_lnL < UNLIKELY + 0.1)
        {
          PhyML_Printf("\n. Move '%s' failed\n",tree->mcmc->move_name[move]);
          assert(FALSE);
        }

      u = Uni();

      for(move=0;move<tree->mcmc->n_moves;move++) if(tree->mcmc->move_weight[move] > u-1.E-10) break;

      assert(!(move == tree->mcmc->n_moves));

      /* printf("\n. %10d %30s %f",tree->mcmc->run,tree->mcmc->move_name[move],tree->mmod->c_lnL); fflush(NULL); */
      /* printf("\n. %10d %30s %f",tree->mcmc->run,tree->mcmc->move_name[move],PHYREX_Lk(tree)); */
      
      if(!strcmp(tree->mcmc->move_name[move],"phyrex_lbda"))
        MCMC_PHYREX_Lbda(tree);

      if(!strcmp(tree->mcmc->move_name[move],"phyrex_mu"))
        MCMC_PHYREX_Mu(tree);

      if(!strcmp(tree->mcmc->move_name[move],"phyrex_rad"))
        MCMC_PHYREX_Radius(tree);

      /* /\* if(!strcmp(tree->mcmc->move_name[move],"phyrex_sigsq")) *\/ */
      /* /\*   MCMC_PHYREX_Sigsq(tree); *\/ */

      if(!strcmp(tree->mcmc->move_name[move],"phyrex_indel_disk"))
        MCMC_PHYREX_Indel_Disk(tree);

      if(!strcmp(tree->mcmc->move_name[move],"phyrex_indel_hit"))
        MCMC_PHYREX_Indel_Hit(tree);

      if(!strcmp(tree->mcmc->move_name[move],"phyrex_move_disk_ud"))
        MCMC_PHYREX_Move_Disk_Updown(tree);

      if(!strcmp(tree->mcmc->move_name[move],"phyrex_swap_disk"))
        MCMC_PHYREX_Swap_Disk(tree);

      if(!strcmp(tree->mcmc->move_name[move],"phyrex_spr"))
        MCMC_PHYREX_Prune_Regraft(tree);

      if(!strcmp(tree->mcmc->move_name[move],"phyrex_scale_times"))
        MCMC_PHYREX_Scale_Times(tree);

      if(!strcmp(tree->mcmc->move_name[move],"phyrex_sim"))
        MCMC_PHYREX_Simulate_Backward(tree);

      if(!strcmp(tree->mcmc->move_name[move],"phyrex_sim_plus"))
        MCMC_PHYREX_Simulate_Backward_Plus(tree);

      if(!strcmp(tree->mcmc->move_name[move],"phyrex_traj"))
        MCMC_PHYREX_Lineage_Traj(tree);

      if(!strcmp(tree->mcmc->move_name[move],"phyrex_lbda_times"))
        MCMC_PHYREX_Lbda_Times(tree);

      if(!strcmp(tree->mcmc->move_name[move],"phyrex_disk_multi"))
        MCMC_PHYREX_Disk_Multi(tree);

      if(!strcmp(tree->mcmc->move_name[move],"phyrex_ldsk_multi"))
        MCMC_PHYREX_Ldsk_Multi(tree);

      if(!strcmp(tree->mcmc->move_name[move],"phyrex_ldsk_and_disk"))
        MCMC_PHYREX_Ldsk_And_Disk(tree);

      if(!strcmp(tree->mcmc->move_name[move],"phyrex_ldsk_given_disk"))
        MCMC_PHYREX_Ldsk_Given_Disk(tree);

      if(!strcmp(tree->mcmc->move_name[move],"phyrex_disk_given_ldsk"))
        MCMC_PHYREX_Disk_Given_Ldsk(tree);

      if(!strcmp(tree->mcmc->move_name[move],"phyrex_indel_disk_serial"))
        MCMC_PHYREX_Indel_Disk_Serial(tree);

      if(!strcmp(tree->mcmc->move_name[move],"phyrex_indel_hit_serial"))
        MCMC_PHYREX_Indel_Hit_Serial(tree);

      if(!strcmp(tree->mcmc->move_name[move],"kappa"))
        MCMC_Kappa(tree);

      if(!strcmp(tree->mcmc->move_name[move],"ras"))
        MCMC_Rate_Across_Sites(tree);

      /* /\* if(!strcmp(tree->mcmc->move_name[move],"phyrex_ldscape_lim")) *\/ */
      /* /\*   MCMC_PHYREX_Ldscape_Limits(tree); *\/ */

      tree->mcmc->run++;
      MCMC_Get_Acc_Rates(tree->mcmc);
      
      if(!(tree->mcmc->run%tree->mcmc->sample_interval))
        {
          /* Lk(NULL,tree); */

          /* RATES_Update_Cur_Bl(tree); */
          /* char *s = Write_Tree(tree,NO); */
          /* PhyML_Fprintf(fp_tree,"\n[%f %f] %s",s,tree->rates->nd_t[tree->n_root->num],tree->c_lnL); */
          /* Free(s); */
          /* fflush(NULL); */

          disk = tree->disk;
          while(disk->prev) disk = disk->prev;

          PhyML_Fprintf(fp_stats,"\n%6d\t%9.1f\t%9.1f\t%9.1f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%G\t%6d\t%6d\t%6d\t%8.1f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\t%G\t%G\t%G",
                        tree->mcmc->run,
                        tree->c_lnL+tree->mmod->c_lnL,
                        tree->c_lnL,
                        tree->mmod->c_lnL,
                        tree->mmod->lbda,
                        tree->mcmc->mode[tree->mcmc->num_move_phyrex_lbda],
                        tree->mmod->mu,
                        tree->mcmc->mode[tree->mcmc->num_move_phyrex_mu],
                        PHYREX_Update_Sigsq(tree),
                        tree->mcmc->mode[tree->mcmc->num_move_phyrex_sigsq],
                        tree->mmod->rad,
                        PHYREX_Neighborhood_Size(tree),
                        PHYREX_Effective_Density(tree),
                        PHYREX_Total_Number_Of_Intervals(tree),
                        PHYREX_Total_Number_Of_Coal_Disks(tree),
                        PHYREX_Total_Number_Of_Hit_Disks(tree),
                        disk->time,
                        tree->mod->kappa->v,
                        tree->mod->ras->alpha->v,
                        tree->mcmc->acc_rate[tree->mcmc->num_move_phyrex_lbda],
                        tree->mcmc->acc_rate[tree->mcmc->num_move_phyrex_mu],
                        tree->mcmc->acc_rate[tree->mcmc->num_move_phyrex_rad],
                        tree->mcmc->acc_rate[tree->mcmc->num_move_phyrex_indel_disk],
                        tree->mcmc->acc_rate[tree->mcmc->num_move_phyrex_indel_hit],
                        tree->mcmc->acc_rate[tree->mcmc->num_move_phyrex_scale_times],
                        tree->mcmc->acc_rate[tree->mcmc->num_move_phyrex_spr],
                        tree->mcmc->acc_rate[tree->mcmc->num_move_phyrex_traj],
                        tree->mcmc->acc_rate[tree->mcmc->num_move_phyrex_sim],
                        tree->mcmc->acc_rate[tree->mcmc->num_move_phyrex_lbda_times],
                        tree->mcmc->acc_rate[tree->mcmc->num_move_phyrex_sim_plus],
                        tree->mcmc->acc_rate[tree->mcmc->num_move_phyrex_indel_hit_serial],

                        tree->mcmc->acc_rate[tree->mcmc->num_move_phyrex_ldsk_given_disk],
                        tree->mcmc->acc_rate[tree->mcmc->num_move_phyrex_disk_given_ldsk],
                        tree->mcmc->acc_rate[tree->mcmc->num_move_phyrex_ldsk_and_disk],
                        tree->mcmc->acc_rate[tree->mcmc->num_move_phyrex_ldsk_multi],
                        tree->mcmc->acc_rate[tree->mcmc->num_move_phyrex_disk_multi],

                        tree->mcmc->tune_move[tree->mcmc->num_move_phyrex_lbda],
                        tree->mcmc->tune_move[tree->mcmc->num_move_phyrex_rad],
                        tree->mcmc->tune_move[tree->mcmc->num_move_phyrex_mu]);

          fflush(fp_stats);
          
          res[0 * tree->mcmc->chain_len / tree->mcmc->sample_interval +  tree->mcmc->run / tree->mcmc->sample_interval] = tree->mmod->lbda; 
          res[1 * tree->mcmc->chain_len / tree->mcmc->sample_interval +  tree->mcmc->run / tree->mcmc->sample_interval] = tree->mmod->mu; 
          res[2 * tree->mcmc->chain_len / tree->mcmc->sample_interval +  tree->mcmc->run / tree->mcmc->sample_interval] = PHYREX_Update_Sigsq(tree); 
          res[3 * tree->mcmc->chain_len / tree->mcmc->sample_interval +  tree->mcmc->run / tree->mcmc->sample_interval] = PHYREX_Neighborhood_Size(tree);
          res[4 * tree->mcmc->chain_len / tree->mcmc->sample_interval +  tree->mcmc->run / tree->mcmc->sample_interval] = tree->mmod->rad;
          res[5 * tree->mcmc->chain_len / tree->mcmc->sample_interval +  tree->mcmc->run / tree->mcmc->sample_interval] = PHYREX_Total_Number_Of_Intervals(tree);
          res[6 * tree->mcmc->chain_len / tree->mcmc->sample_interval +  tree->mcmc->run / tree->mcmc->sample_interval] = PHYREX_Total_Number_Of_Coal_Disks(tree);
          res[7 * tree->mcmc->chain_len / tree->mcmc->sample_interval +  tree->mcmc->run / tree->mcmc->sample_interval] = PHYREX_Total_Number_Of_Hit_Disks(tree);
          res[8 * tree->mcmc->chain_len / tree->mcmc->sample_interval +  tree->mcmc->run / tree->mcmc->sample_interval] = PHYREX_Effective_Density(tree);
          res[9 * tree->mcmc->chain_len / tree->mcmc->sample_interval +  tree->mcmc->run / tree->mcmc->sample_interval] = PHYREX_Coalescence_Rate(tree);

          MCMC_Copy_To_New_Param_Val(tree->mcmc,tree);
          
          for(i=0;i<tree->mcmc->n_moves;i++) if(tree->mcmc->start_ess[i] == YES) MCMC_Update_Effective_Sample_Size(i,tree->mcmc,tree);
          for(i=0;i<tree->mcmc->n_moves;i++) MCMC_Update_Mode(i,tree->mcmc,tree);


          burnin = (int)(0.5*(tree->mcmc->run / tree->mcmc->sample_interval));
          
          rewind(fp_summary);

          PhyML_Fprintf(fp_summary,"\n# SampArea\t NDemes\t TrueLbda\t TrueMu\t TrueSig\t TrueRad\t TrueNeigh\t TrueRhoe\t \t ClockRate\t Diversity\t TrueInt\t TrueCoal\t TrueHits\t RegNeigh\t TrueXroot\t TrueYroot\t TrueHeight\t Lbda5\t Lbda50\t Lbda95\t LbdaMod \t Mu5\t Mu50\t Mu95\t  MuMod \t Sig5\t Sig50\t Sig95\t SigMod \t Neigh5\t Neigh50\t Neigh95\t NeighMod \t Rad5\t Rad50\t Rad95\t Int5\t Int50\t Int95\t Coal5\t Coal50\t Coal95\t Hit5\t Hit50\t Hit95\t Rhoe5\t Rhoe50\t Rhoe95\t CoalRate5\t CoalRate50\t CoalRate95\t ESSLbda \t ESSMu \t ESSSig \t Run");
          
          PhyML_Fprintf(fp_summary,"\n %G\t %d\t %G\t %G\t %G\t %G\t %G\t %G\t %G\t %G\t %d\t %d\t %d\t %G\t %G\t %G\t %G\t ",
                        tot_samp_area,
                        n_demes,
                        true_lbda,
                        true_mu,
                        true_sigsq,
                        true_rad,
                        true_neigh,
                        true_rhoe,
                        tree->rates->clock_r,
                        diversity,
                        true_nint,
                        true_ncoal,
                        true_nhits,
                        fst_neigh,
                        true_root_x,
                        true_root_y,
                        true_height);
          
          PhyML_Fprintf(fp_summary,"%G\t %G\t %G\t %G\t",
                        /* Lbda5 */  Quantile(res+0*tree->mcmc->chain_len / tree->mcmc->sample_interval+burnin,tree->mcmc->run / tree->mcmc->sample_interval+1-burnin,0.025),
                        /* Lbda50 */ Quantile(res+0*tree->mcmc->chain_len / tree->mcmc->sample_interval+burnin,tree->mcmc->run / tree->mcmc->sample_interval+1-burnin,0.50),
                        /* Lbda95 */ Quantile(res+0*tree->mcmc->chain_len / tree->mcmc->sample_interval+burnin,tree->mcmc->run / tree->mcmc->sample_interval+1-burnin,0.975),
                        /* LbdaMod*/ tree->mcmc->mode[tree->mcmc->num_move_phyrex_lbda]);
          
          PhyML_Fprintf(fp_summary,"%G\t %G\t %G\t %G\t",
                        /* mu5 */   Quantile(res+1*tree->mcmc->chain_len / tree->mcmc->sample_interval+burnin,tree->mcmc->run / tree->mcmc->sample_interval+1-burnin,0.025),
                        /* mu50 */  Quantile(res+1*tree->mcmc->chain_len / tree->mcmc->sample_interval+burnin,tree->mcmc->run / tree->mcmc->sample_interval+1-burnin,0.50),
                        /* mu95 */  Quantile(res+1*tree->mcmc->chain_len / tree->mcmc->sample_interval+burnin,tree->mcmc->run / tree->mcmc->sample_interval+1-burnin,0.975),
                        /* muMod */ 2./tree->mcmc->mode[tree->mcmc->num_move_phyrex_mu]);
                    
          PhyML_Fprintf(fp_summary,"%G\t %G\t %G\t %G\t",
                        /* sig5 */   Quantile(res+2*tree->mcmc->chain_len / tree->mcmc->sample_interval+burnin,tree->mcmc->run / tree->mcmc->sample_interval+1-burnin,0.025),
                        /* sig50*/   Quantile(res+2*tree->mcmc->chain_len / tree->mcmc->sample_interval+burnin,tree->mcmc->run / tree->mcmc->sample_interval+1-burnin,0.50),
                        /* sig95*/   Quantile(res+2*tree->mcmc->chain_len / tree->mcmc->sample_interval+burnin,tree->mcmc->run / tree->mcmc->sample_interval+1-burnin,0.975),
                        /* sigMod */ tree->mcmc->mode[tree->mcmc->num_move_phyrex_sigsq]);
          
          PhyML_Fprintf(fp_summary,"%G\t %G\t %G\t %G\t",
                        /* Neigh5 */   Quantile(res+3*tree->mcmc->chain_len / tree->mcmc->sample_interval+burnin,tree->mcmc->run / tree->mcmc->sample_interval+1-burnin,0.025),
                        /* Neigh50*/   Quantile(res+3*tree->mcmc->chain_len / tree->mcmc->sample_interval+burnin,tree->mcmc->run / tree->mcmc->sample_interval+1-burnin,0.50),
                        /* Neigh95*/   Quantile(res+3*tree->mcmc->chain_len / tree->mcmc->sample_interval+burnin,tree->mcmc->run / tree->mcmc->sample_interval+1-burnin,0.975),
                        /* NeighMod */ tree->mcmc->mode[tree->mcmc->num_move_phyrex_mu]);

          
          PhyML_Fprintf(fp_summary,"%G\t %G\t %G\t",
                        /* Rad5 */  Quantile(res+4*tree->mcmc->chain_len / tree->mcmc->sample_interval+burnin,tree->mcmc->run / tree->mcmc->sample_interval+1-burnin,0.025),
                        /* Rad50 */ Quantile(res+4*tree->mcmc->chain_len / tree->mcmc->sample_interval+burnin,tree->mcmc->run / tree->mcmc->sample_interval+1-burnin,0.50),
                        /* Rad95 */ Quantile(res+4*tree->mcmc->chain_len / tree->mcmc->sample_interval+burnin,tree->mcmc->run / tree->mcmc->sample_interval+1-burnin,0.975));


          PhyML_Fprintf(fp_summary,"%G\t %G\t %G\t",
                        /* Int5 */  Quantile(res+5*tree->mcmc->chain_len / tree->mcmc->sample_interval+burnin,tree->mcmc->run / tree->mcmc->sample_interval+1-burnin,0.025),
                        /* Int50 */ Quantile(res+5*tree->mcmc->chain_len / tree->mcmc->sample_interval+burnin,tree->mcmc->run / tree->mcmc->sample_interval+1-burnin,0.50),
                        /* Int95 */ Quantile(res+5*tree->mcmc->chain_len / tree->mcmc->sample_interval+burnin,tree->mcmc->run / tree->mcmc->sample_interval+1-burnin,0.975));

          PhyML_Fprintf(fp_summary,"%G\t %G\t %G\t",
                        /* Coal5 */  Quantile(res+6*tree->mcmc->chain_len / tree->mcmc->sample_interval+burnin,tree->mcmc->run / tree->mcmc->sample_interval+1-burnin,0.025),
                        /* Coal50 */ Quantile(res+6*tree->mcmc->chain_len / tree->mcmc->sample_interval+burnin,tree->mcmc->run / tree->mcmc->sample_interval+1-burnin,0.50),
                        /* Coal95 */ Quantile(res+6*tree->mcmc->chain_len / tree->mcmc->sample_interval+burnin,tree->mcmc->run / tree->mcmc->sample_interval+1-burnin,0.975));

          PhyML_Fprintf(fp_summary,"%G\t %G\t %G\t",
                        /* Hit5 */  Quantile(res+7*tree->mcmc->chain_len / tree->mcmc->sample_interval+burnin,tree->mcmc->run / tree->mcmc->sample_interval+1-burnin,0.025),
                        /* Hit50 */ Quantile(res+7*tree->mcmc->chain_len / tree->mcmc->sample_interval+burnin,tree->mcmc->run / tree->mcmc->sample_interval+1-burnin,0.50),
                        /* Hit95 */ Quantile(res+7*tree->mcmc->chain_len / tree->mcmc->sample_interval+burnin,tree->mcmc->run / tree->mcmc->sample_interval+1-burnin,0.975));
                    
          PhyML_Fprintf(fp_summary,"%G\t %G\t %G\t",
                        /* rhoe5 */  Quantile(res+8*tree->mcmc->chain_len / tree->mcmc->sample_interval+burnin,tree->mcmc->run / tree->mcmc->sample_interval+1-burnin,0.025),
                        /* rhoe50 */ Quantile(res+8*tree->mcmc->chain_len / tree->mcmc->sample_interval+burnin,tree->mcmc->run / tree->mcmc->sample_interval+1-burnin,0.50),
                        /* rhoe95 */ Quantile(res+8*tree->mcmc->chain_len / tree->mcmc->sample_interval+burnin,tree->mcmc->run / tree->mcmc->sample_interval+1-burnin,0.975));

          PhyML_Fprintf(fp_summary,"%G\t %G\t %G\t",
                        /* CoalRate5 */  Quantile(res+9*tree->mcmc->chain_len / tree->mcmc->sample_interval+burnin,tree->mcmc->run / tree->mcmc->sample_interval+1-burnin,0.025),
                        /* CoalRate50 */ Quantile(res+9*tree->mcmc->chain_len / tree->mcmc->sample_interval+burnin,tree->mcmc->run / tree->mcmc->sample_interval+1-burnin,0.50),
                        /* CoalRate95 */ Quantile(res+9*tree->mcmc->chain_len / tree->mcmc->sample_interval+burnin,tree->mcmc->run / tree->mcmc->sample_interval+1-burnin,0.975));

          PhyML_Fprintf(fp_summary,"%G\t %G\t %G\t",
                        tree->mcmc->ess[tree->mcmc->num_move_phyrex_lbda],
                        tree->mcmc->ess[tree->mcmc->num_move_phyrex_mu],  
                        tree->mcmc->ess[tree->mcmc->num_move_phyrex_sigsq]);

          PhyML_Fprintf(fp_summary,"%d\t",tree->mcmc->run);

          PhyML_Fprintf(fp_summary,"\n\n");

          fflush(NULL);
                
          tree->mcmc->sample_num++;
        }



      if(tree->mcmc->run > 2*adjust_len                            &&
         tree->mcmc->sample_num > 1E+2                             &&
         tree->mcmc->ess[tree->mcmc->num_move_phyrex_lbda]  > 10000. &&
         tree->mcmc->ess[tree->mcmc->num_move_phyrex_mu]    > 10000. &&
         tree->mcmc->ess[tree->mcmc->num_move_phyrex_sigsq] > 10000.) break;

      /* if(tree->mcmc->run > tree->mcmc->sample_interval           &&  */
      /*    tree->mcmc->ess[tree->mcmc->num_move_phyrex_lbda]  > 1. && */
      /*    tree->mcmc->ess[tree->mcmc->num_move_phyrex_mu]    > 1. && */
      /*    tree->mcmc->ess[tree->mcmc->num_move_phyrex_sigsq] > 1.) break; */

      (void)signal(SIGINT,MCMC_Terminate);
    }
  while(tree->mcmc->run < tree->mcmc->chain_len);

  fclose(fp_tree);
  fclose(fp_stats);
  fclose(fp_summary);

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

  assert(!(prev == NULL));
  assert(!(next == NULL));    
  
  prev->next = next;
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
  
  disk = tree->disk->prev;  
  while(disk != NULL && disk->time > ins->time) disk = disk->prev;

  assert(!(disk == NULL));

  ins->prev       = disk;
  ins->next       = disk->next;
  disk->next      = ins;
  ins->next->prev = ins;
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
          assert(!tree->disk->next);
          disk = tree->disk;
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
            MAX(0.0,
                MAX(ldsk->coord->lonlat[i] - 2.*rad,
                    o_ldsk->disk->centr->lonlat[i] - rad*(2.*(n_disk_btw-k)-1.)));

          max = 
            MIN(ldsk->disk->mmod->lim->lonlat[i],
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
            MAX(0.0,
                MAX(ldsk->coord->lonlat[i] - 2.*rad,
                    o_ldsk->disk->centr->lonlat[i] - rad*(2.*(n_disk_btw-k)-1.)));

          max = 
            MIN(ldsk->disk->mmod->lim->lonlat[i],
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
      fflush(NULL);
      return(-1);   
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

/* void PHYREX_Update_Lindisk_List(phydbl time, t_ldsk **list, int *pos, t_dsk *disk) */
/* { */
/*   t_dsk *root_dsk; */

/*   *pos = 0; */
/*   root_dsk = disk; */
/*   while(root_dsk->prev) root_dsk = root_dsk->prev; */
/*   /\* printf("\n. root_dsk: %s",root_dsk?root_dsk->id:"xx"); *\/ */
/*   PHYREX_Update_Lindisk_List_Pre(root_dsk->ldsk,time,list,pos); */
/* } */

/* /\*\//////////////////////////////////////////////////////////// */
/* ////////////////////////////////////////////////////////////\*\/ */

/* void PHYREX_Update_Lindisk_List_Pre(t_ldsk *ldsk, phydbl time, t_ldsk **list, int *pos) */
/* { */
/*   /\* printf("\n. time: %f pos: %d ldsk: %s disk: %s n_next: %d ldsk->disk->time: %f disk: %s", *\/ */
/*   /\*        time, *\/ */
/*   /\*        *pos, *\/ */
/*   /\*        ldsk ? ldsk->coord->id : "xx", *\/ */
/*   /\*        ldsk ? ldsk->disk->id : "zz", *\/ */
/*   /\*        ldsk ? ldsk->n_next : -1, *\/ */
/*   /\*        ldsk ? ldsk->disk->time : -1., *\/ */
/*   /\*        ldsk ? ldsk->disk->id : "yy"); fflush(NULL); *\/ */

/*   if(ldsk == NULL) Generic_Exit(__FILE__,__LINE__,__FUNCTION__); */

/*   if((ldsk->prev != NULL) && (ldsk->disk->time > time) && (ldsk->prev->disk->time < time)) */
/*     { */
/*       list[*pos] = ldsk; */
/*       *pos = *pos + 1; */
/*     } */
/*   else if(Are_Equal(ldsk->disk->time,time,1.E-10)) */
/*     { */
/*       list[*pos] = ldsk; */
/*       *pos = *pos + 1;       */
/*     } */
/*   else if(ldsk->disk->time < time) */
/*     { */
/*       int i; */
/*       for(i=0;i<ldsk->n_next;i++) */
/*         PHYREX_Update_Lindisk_List_Pre(ldsk->next[i],time,list,pos);     */
/*     } */
/*   else Generic_Exit(__FILE__,__LINE__,__FUNCTION__); */
    
/* } */

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void PHYREX_Update_Lindisk_List(t_tree *tree)
{
  PHYREX_Update_Lindisk_List_Pre(tree->disk->prev,tree);
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

  if(!disk->next) return;

  for(i=0;i<tree->n_otu;++i) disk->ldsk_a[i] = NULL;

  disk->n_ldsk_a = 0;
  for(i=0;i<disk->next->n_ldsk_a;++i)
    {
      if(disk->next->ldsk_a[i]->prev != disk->ldsk)
        {
          disk->ldsk_a[disk->n_ldsk_a] = disk->next->ldsk_a[i];
          disk->n_ldsk_a++;
        }
    }
  
  /* if(disk->ldsk)  */
  /*   { */
  /*     disk->ldsk_a[disk->n_ldsk_a] = disk->ldsk; */
  /*     disk->n_ldsk_a++;           */
  /*   } */
  
  if(disk->n_ldsk_a > tree->n_otu) 
    {
      PhyML_Fprintf(stderr,"\n. disk: %s next: %s disk->n_ldsk_a: %d coord: %s",disk->id, disk->next->id, disk->n_ldsk_a, disk->ldsk?disk->ldsk->coord->id:"xx");
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
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
  int i,j;
  t_ldsk *ldisk;

  PHYREX_Update_Lindisk_List(tree);

  disk = tree->disk;
  while(disk->prev) disk = disk->prev;
  do
    {
      PhyML_Printf("\n%c Disk: %s @ %7.3f has %3s on it is_coal? %2d rad: %f coord ",
             sign,
             disk->id,
             disk->time,disk->ldsk?disk->ldsk->coord->id:NULL,
             disk->ldsk?disk->ldsk->n_next:-1,
             tree->mmod->rad); fflush(NULL);
      for(j=0;j<tree->mmod->n_dim;j++) PhyML_Printf(" %f",disk->centr->lonlat[j]);
      fflush(NULL);


      for(i=0;i<disk->n_ldsk_a;i++)
        {
          ldisk = disk->ldsk_a[i];

          PhyML_Printf("\n%c ldisk: %s prev: %s",
                       sign,
                       ldisk->coord->id,
                       ldisk->prev ? ldisk->prev->coord->id : NULL);
          fflush(NULL);
          
          for(j=0;j<tree->mmod->n_dim;j++) 
            {
              PhyML_Printf(" %f",ldisk->coord->lonlat[j]);
              fflush(NULL);

              if(FABS(ldisk->coord->lonlat[j] - ldisk->disk->centr->lonlat[j]) > 2.*tree->mmod->rad &&
                 ldisk->disk->ldsk == ldisk) PhyML_Printf(" ! ");

              if(ldisk->prev)
                {
                  if(ldisk->coord->lonlat[j] < ldisk->prev->disk->centr->lonlat[j] - tree->mmod->rad) PhyML_Printf(" #a ");
                  if(ldisk->coord->lonlat[j] > ldisk->prev->disk->centr->lonlat[j] + tree->mmod->rad) PhyML_Printf(" #b ");
                }

              if(ldisk->next)
                {
                  if(ldisk->coord->lonlat[j] < ldisk->disk->centr->lonlat[j] - tree->mmod->rad) PhyML_Printf(" $a ");
                  if(ldisk->coord->lonlat[j] > ldisk->disk->centr->lonlat[j] + tree->mmod->rad) PhyML_Printf(" $b ");
                }
            }
        }

      fflush(NULL);
      disk = disk->next;      
    }
  while(disk);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void PHYREX_Check_Struct(t_tree *tree)
{
  t_dsk *disk;
  int i,j;
  t_ldsk *ldisk;

  disk = tree->disk;
  while(disk->prev) disk = disk->prev;
  do
    {
      /* PHYREX_Update_Lindisk_List(disk->time,disk->ldsk_a,&(disk->n_ldsk_a),disk); */
      PHYREX_Update_Lindisk_List(tree);

      for(i=0;i<disk->n_ldsk_a;i++)
        {
          ldisk = disk->ldsk_a[i];
          if(ldisk->prev != NULL)
            {
              for(j=0;j<tree->mmod->n_dim;j++)
                {
                  if(FABS(ldisk->coord->lonlat[j] - 
                          ldisk->prev->coord->lonlat[j]) > 2.*tree->mmod->rad)
                    {
                      PHYREX_Print_Struct('=',tree);
                      PhyML_Fprintf(stderr,"\n. %f %f %f",
                                    ldisk->coord->lonlat[j], 
                                    ldisk->prev->coord->lonlat[j],
                                    2.*tree->mmod->rad);
                      PhyML_Fprintf(stderr,"\n. Radius: %f",tree->mmod->rad);
                      PhyML_Fprintf(stderr,"\n. Check ldsk %s",ldisk->coord->id);
                      PhyML_Fprintf(stderr,"\n. Centr: %f",ldisk->prev->disk->centr->lonlat[j]);
                      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
                    }
                }
            }
        }

      disk = disk->next;      
    }
  while(disk);
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

  assert(!(tree->disk->next));

  disk = tree->disk;
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

int PHYREX_Total_Number_Of_Hit_Disks(t_tree *tree)
{
  t_dsk *disk;
  int n_hit_disks;

  assert(!(tree->disk->next));

  disk = tree->disk;
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

  assert(!(tree->disk->next));

  disk = tree->disk;
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

  up   = MIN(disk->centr->lonlat[0] + mmod->rad, mmod->lim->lonlat[0]);
  down = MAX(disk->centr->lonlat[0] - mmod->rad, 0.0);
  rght = MIN(disk->centr->lonlat[1] + mmod->rad, mmod->lim->lonlat[1]);
  left = MAX(disk->centr->lonlat[1] - mmod->rad, 0.0);


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

  up   = MIN(disk->centr->lonlat[0] + mmod->rad, mmod->lim->lonlat[0]);
  down = MAX(disk->centr->lonlat[0] - mmod->rad, 0.0);
  rght = MIN(disk->centr->lonlat[1] + mmod->rad, mmod->lim->lonlat[1]);
  left = MAX(disk->centr->lonlat[1] - mmod->rad, 0.0);

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
  phydbl up, down, left, rght;
  int err;

  up   = mmod->lim->lonlat[0];
  down = 0.0;
  rght = mmod->lim->lonlat[1];
  left = 0.0;

  err = NO;
  ldsk->coord->lonlat[0] = Rnorm_Trunc(disk->centr->lonlat[0],mmod->rad,down,up,&err);
  ldsk->coord->lonlat[1] = Rnorm_Trunc(disk->centr->lonlat[1],mmod->rad,left,rght,&err);

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

  /* tree->mmod->c_ln_prior_rad = */
  /*   log(tree->mmod->prior_param_rad) - */
  /*   tree->mmod->prior_param_rad*tree->mmod->rad; */

  /* tree->mmod->c_ln_prior_rad -= log(exp(-tree->mmod->prior_param_lbda*tree->mmod->min_rad)- */
  /*                                   exp(-tree->mmod->prior_param_lbda*tree->mmod->max_rad)); */

  tree->mmod->c_ln_prior_rad = -log(tree->mmod->max_rad - tree->mmod->min_rad);

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

  disk = tree->disk->prev;

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

      loc_min[i] = MAX(0.0,loc_min[i]);
      loc_max[i] = MIN(tree->mmod->lim->lonlat[i],loc_max[i]);
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
          loc_min[i] = 0.0;
          loc_max[i] = tree->mmod->lim->lonlat[i];
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

          loc_min[i] = MAX(0.0,
                           tmp_max - tree->mmod->rad);

          loc_max[i] = MIN(tree->mmod->lim->lonlat[i],
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
/* very short internal edges */
void PHYREX_Ldsk_To_Tree(t_tree *tree)
{
  int i,j;
  t_dsk *disk;

  /* Reset */
  for(i=0;i<2*tree->n_otu-1;++i) 
    {
      for(j=0;j<3;++j) 
        {
          tree->a_nodes[i]->v[j] = NULL;
          tree->a_nodes[i]->b[j] = NULL;
        }
    }

  disk = tree->disk->prev;
  while(disk) 
    {
      if(disk->ldsk) disk->ldsk->nd = NULL;
      disk = disk->prev;
    }

  /* Connect tips */
  for(i=0;i<tree->n_otu;++i) 
    {
      tree->disk->ldsk_a[i]->nd = tree->a_nodes[i];
      tree->a_nodes[i]->coord = tree->disk->ldsk_a[i]->coord;
    }

  disk = tree->disk;
  while(disk->prev) disk = disk->prev;
  
  tree->n_root = tree->a_nodes[2*tree->n_otu-2];
  i = 2*tree->n_otu-3;
  tree->num_curr_branch_available = 0;
  PHYREX_Ldsk_To_Tree_Post(tree->n_root,disk->ldsk,&i,tree);

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
  /* tree->num_curr_branch_available = 0; */
  /* Connect_Edges_To_Nodes_Recur(tree->a_nodes[0],tree->a_nodes[0]->v[0],tree); */

  tree->e_root = NULL;
  for(i=0;i<2*tree->n_otu-3;++i)
    {
      /* printf("\n %d %d",tree->a_edges[i]->left->num,tree->a_edges[i]->rght->num); */
      if((tree->a_edges[i]->left == tree->n_root->v[1] && tree->a_edges[i]->rght == tree->n_root->v[2]) ||
         (tree->a_edges[i]->left == tree->n_root->v[2] && tree->a_edges[i]->rght == tree->n_root->v[1]))
        {        
          tree->e_root = tree->a_edges[i];
          break;
        }
    }
  assert(!(tree->e_root == NULL));
  
  tree->n_root->b[1]  = tree->a_edges[2*tree->n_otu-3];
  tree->n_root->b[2]  = tree->a_edges[2*tree->n_otu-2];

  /* For(i,2*tree->n_otu-1) */
  /*   { */
  /*     printf("\n. * Edge %d %p", */
  /*            tree->a_edges[i]->num, */
  /*            tree->a_edges[i]); */
  /*   } */
  /* PhyML_Printf("\n. tree->n_root->b[1]: %p",tree->n_root->b[1]);  */
  /* PhyML_Printf("\n. tree->n_root->b[2]: %p",tree->n_root->b[2]);  */
  /* fflush(NULL); */
  /* Exit("\n"); */
}
/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void PHYREX_Ldsk_To_Tree_Post(t_node *a, t_ldsk *ldsk, int *available, t_tree *tree)
{
  assert(ldsk);
  assert(a);

  ldsk->nd = a;  
  tree->rates->nd_t[a->num] = ldsk->disk->time;
  a->coord = ldsk->coord;

  if(!ldsk->next) return;
  else
    {
      t_node *parent,*son;
      int n_next;
      t_ldsk *t;

      parent       = a;
      parent->v[1] = NULL;
      parent->v[2] = NULL;
      n_next       = 0;
      do
        {
          t = ldsk->next[n_next]; 

          /* if(t == NULL) */
          /*   { */
          /*     PhyML_Printf("\n. ldsk:%p ldsk->next:%p n_next:%d", */
          /*                  ldsk,ldsk?ldsk->next:NULL,n_next); */
          /*     Generic_Exit(__FILE__,__LINE__,__FUNCTION__); */
          /*   } */


          while(t->next && t->n_next == 1) t = t->next[0];
         
          if(t->nd == NULL) 
            {
              son = tree->a_nodes[*available];
              (*available) = (*available)-1;
            }
          else
            {
              son = t->nd;
            }
                    
          /* printf("\n. a: %d son: %d n_next: %d",a->num,son->num,ldsk->n_next); */

          PHYREX_Ldsk_To_Tree_Post(son,t,available,tree);          

          if(parent->v[2] != NULL && n_next >= 2) 
            {
              t_node *new_parent;
              /* phydbl orig_l2; */

              new_parent = tree->a_nodes[*available];
              (*available) = (*available)-1;              

              new_parent->v[0]       = parent;
              new_parent->v[1]       = parent->v[2];
              new_parent->v[2]       = son;

              parent->v[2]           = new_parent;
              son->v[0]              = new_parent;
              new_parent->v[1]->v[0] = new_parent;
              
              /* printf("\n# connect %d to %d",parent->num,new_parent->num); */
              /* printf("\n# connect %d to %d",new_parent->num,new_parent->v[1]->num); */
              /* printf("\n# connect %d to %d",new_parent->num,new_parent->v[2]->num);               */
              /* fflush(NULL); */
              
              tree->rates->nd_t[new_parent->num] = ldsk->disk->time;

              parent = new_parent;
            }
          else
            {
              son->v[0] = parent;          
              if(!parent->v[1]) parent->v[1] = son;
              else              parent->v[2] = son;
              
              /* printf("\n. connect %d to %d",parent->num,son->num); */
            }

          n_next++;
        }
      while(n_next != ldsk->n_next);
    }
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
  disk    = tree->disk;
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
                         tree->mmod->lim->lonlat[k]);
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

  assert(!(tree->disk->next));

  disk = tree->disk->prev;
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
  for(j=0;j<mmod->n_dim;j++) if(ldsk->coord->lonlat[j] > mmod->lim->lonlat[j] ||
                        ldsk->coord->lonlat[j] < 0.0) return NO;
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
          
          dist[pair] = Euclidean_Dist(tree->a_nodes[i]->coord,tree->a_nodes[j]->coord);
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
      dist = Euclidean_Dist(tree->a_nodes[i]->coord,tree->a_nodes[j]->coord);
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
               POW(tree->mmod->rad,4)*
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

  h    = tree->mmod->lim->lonlat[0];
  w    = tree->mmod->lim->lonlat[1];
  lbda = tree->mmod->lbda;
  mu   = tree->mmod->mu;
  a    = tree->mmod->min_sigsq;
  b    = tree->mmod->max_sigsq;

  u = Uni();

  return(POW(((b-a)*u+a)*h*w/(4.*PI*lbda*mu),0.25));

}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void PHYREX_Read_Tip_Coordinates(t_ldsk **ldsk_a, t_tree *tree)
{
  char *s;
  FILE *fp;
  int i,*done,found_sw,found_ne;
  phydbl sw_lon, sw_lat,  ne_lon, ne_lat;

  s    = (char *)mCalloc(T_MAX_LINE,sizeof(char));
  fp   = tree->io->fp_in_coord;
  done = (int *)mCalloc(tree->n_otu,sizeof(int));

  found_sw = NO;
  found_ne = NO;

  for(i=0;i<tree->n_otu;i++) done[i] = NO;

  do
    {
      if(fscanf(fp,"%s",s) == EOF) break;
      For(i,strlen(s)) if(s[i] == '#') break; /* skip comment */
      if(i != strlen(s)) continue;
      
      for(i=0;i<tree->n_otu;i++) if(strstr(tree->a_nodes[i]->name,s)) break;

      if(i != tree->n_otu) /* Found a match */
        {
          if(fscanf(fp,"%lf",&(ldsk_a[i]->coord->lonlat[0])) == EOF) break;
          if(fscanf(fp,"%lf",&(ldsk_a[i]->coord->lonlat[1])) == EOF) break;          
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

  for(i=0;i<tree->n_otu;i++) 
    {
      ldsk_a[i]->coord->lonlat[0] -= sw_lon;
      ldsk_a[i]->coord->lonlat[1] -= sw_lat;

      ldsk_a[i]->coord->lonlat[0] /= (ne_lon - sw_lon);
      ldsk_a[i]->coord->lonlat[1] /= (ne_lat - sw_lat);

      ldsk_a[i]->coord->lonlat[0] *= tree->mmod->lim->lonlat[0];
      ldsk_a[i]->coord->lonlat[1] *= tree->mmod->lim->lonlat[1];

      PhyML_Printf("\n. Scaled coordinates of '%-50s': %12f\t %12f",
                   tree->a_nodes[i]->name,
                   ldsk_a[i]->coord->lonlat[0],
                   ldsk_a[i]->coord->lonlat[1]);
    }

  Free(s);
  Free(done);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl PHYREX_Rate_Per_Unit_Area(t_tree *tree)
{
  int i;
  phydbl denom;

  denom = tree->mmod->lim->lonlat[0];
  
  for(i=1;i<tree->mmod->n_dim;i++) denom *= tree->mmod->lim->lonlat[i];
  
  return(tree->mmod->lbda / denom);

}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl PHYREX_Tree_Height(t_tree *tree)
{
  t_dsk *disk;

  disk = tree->disk;
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

t_ldsk *PHYREX_Remove_Path(t_ldsk *beg, t_ldsk *end, int *pos_end, t_tree *tree)
{
  t_ldsk *ldsk;
  int dir_end_beg;

  dir_end_beg = PHYREX_Get_Next_Direction(beg,end);
  *pos_end = dir_end_beg;
  PHYREX_Remove_Lindisk_Next(end,end->next[dir_end_beg]);

  if(beg->prev == end) return NULL;
  
  /* PhyML_Printf("\n- rm beg: %12f %12f %12f %s", */
  /*              beg->coord->lonlat[0],beg->coord->lonlat[1],beg->disk->time, */
  /*              beg->coord->id); */
  /* PhyML_Printf("\n- rm end: %12f %12f %12f %s", */
  /*              end->coord->lonlat[0],end->coord->lonlat[1],end->disk->time, */
  /*              end->coord->id); */

  ldsk = beg->prev;
  while(1)
    {
      PHYREX_Remove_Disk(ldsk->disk);
      /* PhyML_Printf("\n- rm %12f %12f %s",ldsk->coord->lonlat[0],ldsk->coord->lonlat[1],ldsk->coord->id); */
      if(ldsk->prev == end) break;
      ldsk = ldsk->prev;
    }
  
  ldsk->prev = NULL;
  ldsk = beg->prev;
  beg->prev = NULL;

  return(ldsk);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void PHYREX_Insert_Path(t_ldsk *beg, t_ldsk *end, t_ldsk *path, int pos, t_tree *tree)
{
  t_ldsk *ldsk;

  assert(!(path == beg || path == end));

  if(path == NULL) 
    {
      beg->prev = end;      

      /* printf("\n. in (beg-end) %12f %12f %12f %12f %s %s", */
      /*        beg->coord->lonlat[0], */
      /*        beg->coord->lonlat[1], */
      /*        end->coord->lonlat[0], */
      /*        end->coord->lonlat[1], */
      /*        beg->coord->id,end->coord->id); */

      /* PHYREX_Insert_Ldsk_In_Next_List(beg,end->n_next,end); */
      PHYREX_Insert_Ldsk_In_Next_List(beg,pos,end);
    }
  else
    {

      /* printf("\n. in beg %12f %12f %12f %s", */
      /*        beg->coord->lonlat[0], */
      /*        beg->coord->lonlat[1], */
      /*        beg->disk->time, */
      /*        beg->coord->id); */

      /* printf("\n. in end %12f %12f %12f %s", */
      /*        end->coord->lonlat[0], */
      /*        end->coord->lonlat[1], */
      /*        end->disk->time, */
      /*        end->coord->id); */
      
      ldsk = path;
      while(1)
        {
          PHYREX_Insert_Disk(ldsk->disk,tree);
          /* printf("\n+ in %12f %12f %12f %s", */
          /*        ldsk->coord->lonlat[0], */
          /*        ldsk->coord->lonlat[1], */
          /*        ldsk->disk->time, */
          /*        ldsk->coord->id); */
          if(ldsk->prev == NULL) break;
          ldsk = ldsk->prev;
        }
      
      beg->prev = path;
      ldsk->prev = end;
      path->next[0] = beg;

      /* PHYREX_Insert_Ldsk_In_Next_List(ldsk,end->n_next,end); */
      PHYREX_Insert_Ldsk_In_Next_List(ldsk,pos,end);
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

t_ldsk *PHYREX_Generate_Path(t_ldsk *beg, t_ldsk *end, phydbl cur_n_evt, phydbl sd, t_tree *tree)
{
  int n_evt,i,j,swap,err;
  phydbl dt,*time,dum,mode;
  t_ldsk *path,**ldsk_a;
  t_dsk *disk;

  dt = FABS(beg->disk->time - end->disk->time);

  /* How many hit events ? */
  if(cur_n_evt < .0)
    /* Not sure that rate is ok when considering landscape with boundaries... */
    n_evt = Rpois(dt*2.*PHYREX_Rate_Per_Unit_Area(tree)*PI*POW(sd,2)*tree->mmod->mu); 
  else
    n_evt = Rpois(cur_n_evt);


  if(n_evt <= 0) return(NULL);

  time   = (phydbl *)mCalloc(n_evt,sizeof(phydbl));
  ldsk_a = (t_ldsk **)mCalloc(n_evt,sizeof(t_ldsk *));

  for(i=0;i<n_evt;i++) time[i] =  beg->disk->time - FABS(Uni()*(end->disk->time - beg->disk->time));

  /* Invert time direction */
  for(i=0;i<n_evt;i++) time[i] = -time[i];
  
  /* Bubble sort time in ascending order */
  do
    {
      swap = NO;
      for(i=0;i<n_evt-1;i++) 
        {
          if(time[i+1] < time[i])
            {
              swap = YES;
              dum       = time[i+1];
              time[i+1] = time[i];
              time[i]   = dum;
            }
        }
    }while(swap == YES);

  
  for(i=0;i<n_evt;i++)
    {
      ldsk_a[i] = PHYREX_Make_Lindisk_Node(tree->mmod->n_dim);
      disk = PHYREX_Make_Disk_Event(tree->mmod->n_dim,tree->n_otu);
      PHYREX_Init_Lindisk_Node(ldsk_a[i],disk,tree->mmod->n_dim);
      PHYREX_Make_Lindisk_Next(ldsk_a[i]);
      PHYREX_Init_Disk_Event(disk,tree->mmod->n_dim,tree->mmod);
      disk->ldsk = ldsk_a[i];
      disk->time = -time[i];      
      /* printf("\n. Generate ldsk %s",ldsk_a[i]->coord->id); */
    }

  for(i=0;i<n_evt-1;i++) ldsk_a[i]->prev = ldsk_a[i+1];
  ldsk_a[i]->prev = NULL;

  for(i=1;i<n_evt;i++) ldsk_a[i]->next[0] = ldsk_a[i-1];
  ldsk_a[0]->next[0] = NULL;

  path = ldsk_a[0];

  /* Generate path */
  for(i=0;i<tree->mmod->n_dim;i++)
    {
      for(j=0;j<n_evt;j++)
        {
          if(j == 0)
            mode = (end->coord->lonlat[i] - beg->coord->lonlat[i])/(n_evt+1.) + beg->coord->lonlat[i];          
          else
            mode = (end->coord->lonlat[i] - ldsk_a[j-1]->coord->lonlat[i])/(n_evt+1.-j) + ldsk_a[j-1]->coord->lonlat[i];          

          ldsk_a[j]->coord->lonlat[i] = Rnorm_Trunc(mode,
                                                    sd,
                                                    0.0,
                                                    tree->mmod->lim->lonlat[i],&err);
          

          /* ldsk_a[j]->coord->lonlat[i] = Uni()*tree->mmod->lim->lonlat[i]; */
          /* ldsk_a[j]->coord->lonlat[i] = mode; */

          ldsk_a[j]->disk->centr->lonlat[i] = Rnorm_Trunc(ldsk_a[j]->coord->lonlat[i],
                                                          sd,
                                                          0.0,
                                                          tree->mmod->lim->lonlat[i],&err);

          /* ldsk_a[j]->disk->centr->lonlat[i] = Uni()*tree->mmod->lim->lonlat[i]; */
          /* ldsk_a[j]->disk->centr->lonlat[i] = ldsk_a[j]->coord->lonlat[i]; */
        }
    }


  /* for(j=0;j<n_evt;j++) PhyML_Printf("\n. in %12f %12f",ldsk_a[j]->coord->lonlat[0],ldsk_a[j]->coord->lonlat[1]); */


  Free(ldsk_a);
  Free(time);

  return(path);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl PHYREX_Path_Logdensity(t_ldsk *beg, t_ldsk *end, phydbl cur_n_evt, phydbl sd, t_tree *tree)
{
  int i,j,err,n_evt;
  t_ldsk *ldsk;
  phydbl lnDens,mode,rate;

  lnDens = 0.0;
  mode   = 0.0;

  n_evt = 0;
  ldsk = beg->prev;
  while(ldsk != end)
    {
      n_evt++;
      ldsk = ldsk->prev;
      assert(ldsk != NULL);
    }

  for(i=0;i<tree->mmod->n_dim;i++)
    {     
      j    = 0;
      ldsk = beg;
      while(ldsk->prev != end)
        {
          assert(!(ldsk == NULL));
          
          mode = (end->coord->lonlat[i] - ldsk->coord->lonlat[i])/(n_evt+1.-j) + ldsk->coord->lonlat[i];          

          lnDens += Log_Dnorm_Trunc(ldsk->prev->coord->lonlat[i],
                                    mode,
                                    sd,
                                    0.0,
                                    tree->mmod->lim->lonlat[i],&err);

          /* lnDens += log(1./tree->mmod->lim->lonlat[i]); */

          lnDens += Log_Dnorm_Trunc(ldsk->prev->disk->centr->lonlat[i],
                                    ldsk->prev->coord->lonlat[i],
                                    sd,
                                    0.0,
                                    tree->mmod->lim->lonlat[i],&err);

          /* lnDens += log(1./tree->mmod->lim->lonlat[i]); */

          ldsk = ldsk->prev;
          j++;
        }
    }

  if(cur_n_evt < 0)
    rate = 2.*PHYREX_Rate_Per_Unit_Area(tree)*PI*POW(sd,2)*tree->mmod->mu*FABS(end->disk->time - beg->disk->time);
  else
    rate = cur_n_evt;

  lnDens += Dpois(n_evt,rate,YES);
  lnDens += (n_evt) * log(1./FABS(end->disk->time - beg->disk->time));
  lnDens += LnFact(n_evt);

  return(lnDens);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl PHYREX_Time_Tree_Length(t_tree *tree)
{
  phydbl len;
  int i;
  t_dsk *disk;

  disk = tree->disk;
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

int PHYREX_Is_On_Path(t_ldsk *target, t_ldsk *beg, t_ldsk *end)
{
  t_ldsk *ldsk;

  if(target == beg || target == end) return NO;

  assert(!(beg->disk->time < end->disk->time));

  ldsk = beg->prev; 
  while(ldsk != end)
    {
      if(ldsk == target) return YES;
      ldsk = ldsk->prev;
      assert(!(ldsk == NULL));
    }
  return NO;
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

int PHYREX_Path_Len(t_ldsk *beg, t_ldsk *end)
{
  t_ldsk *ldsk;
  int len;

  len = 0;
  ldsk = beg;
  while(ldsk != end)
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

  disk = tree->disk->prev;

  do
    {
      PhyML_Printf("\n. Disk: %p time: %12f lk: %12f cumlk: %12f",
                   disk,
                   disk->time,
                   PHYREX_Lk_Core(disk,tree),
                   disk->c_lnL);
      
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
  return(FABS(d->disk->time - lca->disk->time));      
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl PHYREX_Dist_Between_Two_Ldsk(t_ldsk *n1,  t_ldsk *n2, t_tree *tree)
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
      s = strrchr(tree->a_nodes[i]->coord->id,'_');
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
                   tree->a_nodes[i]->coord->id,
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
      s = strrchr(tree->a_nodes[i]->coord->id,'_');
      PhyML_Fprintf(fp,"%s=%s",
                   tree->a_nodes[i]->coord->id,
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

  disk = tree->disk;

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
  w     = tree->mmod->lim->lonlat[0];
  h     = tree->mmod->lim->lonlat[1];

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
/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
