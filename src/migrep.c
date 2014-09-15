/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homoLOGous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

/* Routines that implement Etheridge and Barton's model of continuous-space
   coalescent.
*/

#include "migrep.h"

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


int MIGREP_Draw(GtkWidget *widget, cairo_t *cr, gpointer *data)
{   
  t_tree *tree;
  phydbl h, w;
  t_disk_evt *devt;
  t_migrep_mod *mmod;
  int i;

  tree = (t_tree *)(data);
  
  mmod = tree->mmod;

  devt = tree->devt;
  while(devt->prev) devt = devt->prev;

  h = FABS(devt->time);
  w = tree->mmod->lim->lonlat[1];

  cairo_select_font_face(cr, "Courier",
                         CAIRO_FONT_SLANT_NORMAL,
                         CAIRO_FONT_WEIGHT_NORMAL);
  
  cairo_set_font_size(cr,0.03);
  cairo_scale (cr,WINDOW_WIDTH/1.2,WINDOW_HEIGHT/1.2);
  cairo_translate(cr,0.05,0.05);
  cairo_set_source_rgb(cr,0,0,0);
  cairo_set_line_width(cr,0.002);

  while(devt->next) devt = devt->next;
  do
    {      
      cairo_move_to(cr,0,FABS(devt->time)/h);
      cairo_show_text(cr,devt->id); 
      cairo_stroke(cr);

      cairo_move_to(cr,devt->centr->lonlat[1]/w - mmod->rad/w,FABS(devt->time)/h);
      cairo_set_line_width(cr,0.004);
      cairo_set_source_rgba(cr, 1, 0, 0, 0.3);
      cairo_set_line_width(cr,0.002);
      cairo_line_to(cr,devt->centr->lonlat[1]/w + mmod->rad/w,FABS(devt->time)/h);
      cairo_stroke(cr);
      cairo_arc(cr,devt->centr->lonlat[1]/w,FABS(devt->time)/h,0.005,0.,2.*PI);
      cairo_stroke(cr);
      cairo_set_source_rgb(cr, 0, 0, 0);
          
      For(i,devt->n_ldsk_a)
        {
          cairo_move_to(cr,devt->ldsk_a[i]->coord->lonlat[1]/w,FABS(devt->time)/h);
          cairo_show_text(cr,devt->ldsk_a[i]->coord->id); 
          cairo_stroke(cr);
          cairo_move_to(cr,devt->ldsk_a[i]->coord->lonlat[1]/w,FABS(devt->time)/h);
          cairo_line_to(cr,devt->ldsk_a[i]->coord->lonlat[1]/w,FABS(devt->prev->time)/h);
          cairo_stroke(cr);

          if(devt->ldsk_a[i]->prev->devt == devt->prev)
            {
              cairo_line_to(cr,devt->ldsk_a[i]->coord->lonlat[1]/w,FABS(devt->prev->time)/h);
              cairo_line_to(cr,devt->ldsk_a[i]->prev->coord->lonlat[1]/w,FABS(devt->prev->time)/h);
              cairo_stroke(cr);
            }

          if(devt->ldsk_a[i]->is_coal && devt->ldsk_a[i]->devt == devt) 
            {
              cairo_set_source_rgba(cr, 0, 0.2, 0.8, 0.5);
              cairo_arc(cr,devt->ldsk_a[i]->coord->lonlat[1]/w,FABS(devt->time)/h,0.01,0.,2.*PI);
              cairo_fill (cr);
              cairo_stroke (cr);
              cairo_set_source_rgb (cr, 0, 0, 0);
            }
          else
            {
              cairo_set_source_rgba(cr, 0.2, 0.2, 0.2, 0.8);
              cairo_arc(cr,devt->ldsk_a[i]->coord->lonlat[1]/w,FABS(devt->time)/h,0.005,0.,2.*PI);
              cairo_fill (cr);
              cairo_stroke (cr);
              cairo_set_source_rgb (cr, 0, 0, 0);
            }
        }
      devt = devt->prev;
      if(devt->prev == NULL) break;
    }while(1);

  cairo_move_to(cr,devt->ldsk_a[0]->coord->lonlat[1]/w,FABS(devt->time)/h);
  cairo_show_text(cr,devt->ldsk_a[0]->coord->id); 
  cairo_arc(cr,devt->ldsk_a[0]->coord->lonlat[1]/w,FABS(devt->time)/h,0.005,0.,2.*PI);
  cairo_stroke (cr);

  /* cairo_scale (cr,WINDOW_WIDTH,WINDOW_HEIGHT); */
  /* cairo_set_source_rgb (cr, 0, 0, 0); */
  /* cairo_move_to (cr, 0, 0); */
  /* cairo_line_to (cr, 1, 1); */
  /* cairo_move_to (cr, 1, 0); */
  /* cairo_line_to (cr, 0, 1); */
  /* cairo_set_line_width (cr, 0.2); */
  /* cairo_stroke (cr); */
  
  /* cairo_rectangle (cr, 0, 0, 0.5, 0.5); */
  /* cairo_set_source_rgba (cr, 1, 0, 0, 0.80); */
  /* cairo_fill (cr); */
  
  /* cairo_rectangle (cr, 0, 0.5, 0.5, 0.5); */
  /* cairo_set_source_rgba (cr, 0, 1, 0, 0.60); */
  /* cairo_fill (cr); */
  
  /* cairo_rectangle (cr, 0.5, 0, 0.5, 0.5); */
  /* cairo_set_source_rgba (cr, 0, 0, 1, 0.40); */
  /* cairo_fill (cr); */
  
  /* /\* Set color for background *\/ */
  /* cairo_set_source_rgb(cr, 1, 1, 1); */
  /* /\* fill in the background color*\/ */
  /* cairo_paint(cr); */
  
  /* /\* set color for rectangle *\/ */
  /* cairo_set_source_rgb(cr, 0.42, 0.65, 0.80); */
  /* /\* set the line width *\/ */
  /* cairo_set_line_width(cr,6); */
  /* /\* draw the rectangle's path beginning at 3,3 *\/ */
  /* cairo_rectangle (cr, 3, 3, 100, 100); */
  /* /\* stroke the rectangle's path with the chosen color so it's actually visible *\/ */
  /* cairo_stroke(cr); */
  
  /* /\* draw circle *\/ */
  /* cairo_set_source_rgb(cr, 0.17, 0.63, 0.12); */
  /* cairo_set_line_width(cr,2); */
  /* cairo_arc(cr, 150, 210, 20, 0, 2*G_PI); */
  /* cairo_stroke(cr); */
  
  /* /\* draw horizontal line *\/ */
  /* cairo_set_source_rgb(cr, 0.77, 0.16, 0.13); */
  /* cairo_set_line_width(cr, 6); */
  /* cairo_move_to(cr, 80,160); */
  /* cairo_line_to(cr, 200, 160); */
  /* cairo_stroke(cr); */
  
  /* /\* cairo_set_source_rgb(cr,1,1,0); *\/ */
  /* /\* cairo_paint(cr); *\/ */
  
  return FALSE;
}

int MIGREP_Main(int argc, char *argv[])
{
  GtkWidget *window;
  GtkWidget *draw_area;
  pthread_t calc_tid;
  t_tree *tree;
  int seed;

  seed = time(NULL);
  /* seed = 1409782620; */
  /* seed = 32680; */
  /* seed = 32544; */
  /* seed = 32615; */
  /* seed = 3; */
  /* seed = 1410475473; */
  seed = 1410476172;
  printf("\n. Seed: %d",seed);
  srand(seed);
  tree = MIGREP_Simulate_Backward((int)atoi(argv[1]),10.,10.);

  printf("\n. seed = %d",seed);

  /* gdk_threads_init(); */
  /* gdk_threads_enter(); */
  
  gtk_init(&argc,&argv);
  
  window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
  g_signal_connect(window,"destroy",G_CALLBACK(gtk_main_quit),NULL);
  gtk_container_set_border_width (GTK_CONTAINER (window),0);
  
  draw_area= gtk_drawing_area_new();
  gtk_widget_set_size_request(draw_area,WINDOW_WIDTH,WINDOW_HEIGHT);
  g_signal_connect(draw_area,"draw",G_CALLBACK(MIGREP_Draw),(gpointer)tree);
  
  gtk_container_add(GTK_CONTAINER(window),draw_area);

  gtk_widget_show_all (window);

  tree->draw_area = draw_area;

  pthread_create(&calc_tid,NULL,(void *(*)(void *))MIGREP_MCMC,tree);
    
  gtk_main();
  /* gdk_threads_leave(); */
  
  return 0;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// Simulate Etheridge-Barton model backward in time, following n_otu lineages
// on a rectangle of dimension width x height
// See Kelleher, Barton & Etheridge, Bioinformatics, 2013.
t_tree *MIGREP_Simulate_Backward(int n_otu, phydbl width, phydbl height)
{  
  int i;
  t_tree *tree;
  int n_dim,n_lineages,n_lineages_new,n_devt,n_hit;
  phydbl curr_t,dt_disk;
  t_migrep_mod *mmod;
  t_disk_evt *devt;
  t_lindisk_nd **ldsk_a,**ldsk_a_tmp,*new_ldsk;
  t_node *nd;

  n_dim = 2; // 2-dimensional landscape

  tree = Make_Tree_From_Scratch(n_otu,NULL);
  nd   = Make_Node_Light(0);
  devt = MIGREP_Make_Disk_Event(n_dim);
  MIGREP_Init_Disk_Event(devt);
  
  // Allocate coordinates for all the tips first (will grow afterwards)
  ldsk_a = (t_lindisk_nd **)mCalloc(n_otu,sizeof(t_lindisk_nd *));
  For(i,n_otu) 
    {
      ldsk_a[i] = MIGREP_Make_Lindisk_Node(n_dim);
      MIGREP_Init_Lindisk_Node(ldsk_a[i],devt,n_dim);
    }

  ldsk_a_tmp = (t_lindisk_nd **)mCalloc(n_otu,sizeof(t_lindisk_nd *));

  // Generate coordinates for the tip nodes (uniform distribution on the rectangle)
  For(i,n_otu)
    {
      /* ldsk_a[i]->coord->lonlat[0] = Uni()*width;  // longitude */
      /* ldsk_a[i]->coord->lonlat[1] = Uni()*height; // latitude */
      ldsk_a[i]->coord->lonlat[0] = (i+1.)/(n_otu+1.)*width;  // longitude
      ldsk_a[i]->coord->lonlat[1] = (i+1.)/(n_otu+1.)*height; // latitude
    }

  // Allocate migrep model
  mmod = MIGREP_Make_Migrep_Model(n_dim);
  MIGREP_Init_Migrep_Mod(mmod,n_dim);
  mmod->lim->lonlat[0] = width;
  mmod->lim->lonlat[1] = height;
  
  // First disk event (at time 0)
  devt->time             = 0.0;
  devt->ldsk_a           = ldsk_a;
  devt->mmod             = mmod;
  devt->n_ldsk_a         = n_otu;
  devt->centr->lonlat[0] = .5*width;
  devt->centr->lonlat[1] = .5*height;      
  
  // Allocate and initialise for next event
  devt->prev = MIGREP_Make_Disk_Event(n_dim);
  MIGREP_Init_Disk_Event(devt->prev);
  devt->prev->next = devt;
  
  For(i,n_otu)
    {
      printf("\nx disk %s [%15f] %3d %s %15f %15f",
             devt?devt->id:"xxx",
             devt?devt->time:0.0,
             i,
             ldsk_a[i]->coord->id,
             ldsk_a[i]->coord->lonlat[0],
             ldsk_a[i]->coord->lonlat[1]);
    }

  // Move to it
  devt = devt->prev;

  // Initialize parameters of migrep model
  mmod->lbda = 0.2;
  mmod->mu   = 1.0;
  mmod->rad  = 3.0;
  
  curr_t      = 0.0;
  dt_disk     = 0.0;
  n_lineages  = n_otu;
  n_devt      = 0;
  do
    {
      // Time of next event
      dt_disk = Rexp(mmod->lbda);
      curr_t -= dt_disk;
      
      // Coordinates of next event
      GEO_Init_Coord(devt->centr,n_dim);
      devt->centr->lonlat[0] = Uni()*width;
      devt->centr->lonlat[1] = Uni()*height;      

      devt->time = curr_t;
      devt->mmod = mmod;

      /* printf("\n. Disk %s has %d lindisk nodes and %d disks under",devt->id,devt->n_ldsk_a,devt->n_disk_under); */

      // New lindisk (will not be used is no lineage is hit) 
      new_ldsk = MIGREP_Make_Lindisk_Node(n_dim);
      MIGREP_Init_Lindisk_Node(new_ldsk,devt,n_dim);
      do
        {
          Runif_Disk(new_ldsk->coord->lonlat,
                     new_ldsk->coord->lonlat+1,
                     devt->centr->lonlat[0],
                     devt->centr->lonlat[1],
                     devt->mmod->rad);
          // Check that the next lindisk is within boundaries
          if(new_ldsk->coord->lonlat[0] > 0.0 &&
             new_ldsk->coord->lonlat[0] < devt->mmod->lim->lonlat[0] &&
             new_ldsk->coord->lonlat[1] > 0.0 &&
             new_ldsk->coord->lonlat[1] < devt->mmod->lim->lonlat[0]) break;
        }
      while(1);

      n_hit          = 0;
      n_lineages_new = 0;
      For(i,n_lineages)
        {
          ldsk_a[i]->prev = NULL;

          if(MIGREP_Is_In_Disk(ldsk_a[i]->coord,devt) == YES) 
            {
              if(Uni() < mmod->mu)
                {
                  printf("\n. Hit and die %s @ time %f. Center: %f %f Go to %f %f",
                         ldsk_a[i]->coord->id,
                         devt->time,
                         devt->centr->lonlat[0],
                         devt->centr->lonlat[1],
                         new_ldsk->coord->lonlat[0],
                         new_ldsk->coord->lonlat[1]);
                  
                  MIGREP_Make_Lindisk_Next(new_ldsk);

                  ldsk_a[i]->prev                    = new_ldsk;
                  ldsk_a[i]->prev->is_hit            = YES;
                  new_ldsk->next[new_ldsk->n_next-1] = ldsk_a[i]; 
                  devt->ldsk                         = new_ldsk;

                  n_hit++;

                  if(n_hit == 1)
                    {
                      ldsk_a_tmp[n_lineages_new] = new_ldsk;
                      n_lineages_new++;
                    }
                }
            }
          
          if(ldsk_a[i]->prev == NULL)
            {
              ldsk_a[i]->prev = ldsk_a[i];
              /* ldsk_a[i]->prev = MIGREP_Make_Lindisk_Node(n_dim); */
              /* MIGREP_Init_Lindisk_Node(ldsk_a[i]->prev,n_dim); */
              /* For(j,mmod->n_dim) ldsk_a[i]->prev->coord->lonlat[j] = ldsk_a[i]->coord->lonlat[j]; */
              /* ldsk_a[i]->prev->next = ldsk_a[i]; */
              ldsk_a_tmp[n_lineages_new] = ldsk_a[i]->prev;
              n_lineages_new++;
            }          
        }
      
      if((n_hit > 0) && (n_lineages_new != n_lineages - n_hit + 1))
        {
          PhyML_Printf("\n== Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
          Warn_And_Exit("");
        }

      if(n_hit > 0) n_lineages -= (n_hit-1);

      if(n_hit > 1) // a coalescent event occurred
        {
          devt->nd          = nd;
          new_ldsk->is_coal = YES;
          PhyML_Printf("\n. Coalescent @ disk %s on ldsk %s",devt->id,devt->ldsk->coord->id);
        }

 
      ldsk_a = (t_lindisk_nd **)mCalloc(n_lineages,sizeof(t_lindisk_nd *));
      For(i,n_lineages) ldsk_a[i] = ldsk_a_tmp[i];
      
      devt->ldsk_a = ldsk_a;
      devt->n_ldsk_a = n_lineages;

      For(i,n_lineages)
        {
          printf("\n. disk %s [%15f] %3d %s %15f %15f",
                 devt?devt->id:"xxx",
                 devt?devt->time:0.0,
                 i,
                 ldsk_a[i]->coord->id,
                 ldsk_a[i]->coord->lonlat[0],
                 ldsk_a[i]->coord->lonlat[1]);
        }

      n_devt++;

      if(n_lineages == 1) break;

      devt->prev = MIGREP_Make_Disk_Event(n_dim);
      MIGREP_Init_Disk_Event(devt->prev);
      devt->prev->next = devt;
      
      devt = devt->prev;          
    }
  while(1);

  while(devt->next) devt = devt->next;

  tree->devt = devt;
  tree->mmod = mmod;


  /* while(devt->prev) */
  /*   { */
  /*     printf("\n<><><>"); */
  /*     printf("\n devt %f %s",devt->time,devt->id); */
  /*     For(i,devt->n_ldsk_a) */
  /*       { */
  /*         printf("\n. %s %f %f", */
  /*            devt->ldsk_a[i]->coord->id, */
  /*            devt->ldsk_a[i]->coord->lonlat[0], */
  /*            devt->ldsk_a[i]->coord->lonlat[1]); */
  /*       } */
  /*     devt = devt->prev; */
  /*   } */


  MIGREP_Lk(devt,mmod);
  printf("\n. LK: %f",mmod->c_lnL);
  /* Exit("\n"); */
  /* MIGREP_MCMC(tree); */

  return(tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int MIGREP_Is_In_Disk(t_geo_coord *coord, t_disk_evt *devt)
{
  int i;

  For(i,devt->centr->dim)
    {
      if(FABS(coord->lonlat[i] - devt->centr->lonlat[i]) > devt->mmod->rad)
        {
          return(NO);
        }
    }
  return(YES);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl MIGREP_Lk(t_disk_evt *devt, t_migrep_mod *mmod)
{
  phydbl lnL;
  phydbl log_lbda,log_pi_rad,log_mu,log_one_mu;
  t_lindisk_nd *lindisk_nd;
  int i;
  short int was_hit, n_hit;


  // Rewind if necessary
  while(devt->next) { devt = devt->next; }  

  if(!devt->prev)
    {
      PhyML_Printf("\n== Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
      Warn_And_Exit("");
    }

  log_lbda   = LOG(mmod->lbda);
  log_mu     = LOG(mmod->mu);
  log_pi_rad = LOG(PI*mmod->rad*mmod->rad);
  log_one_mu = LOG(1. - mmod->mu);  
  lnL        = 0.0;

  do
    {
      /* PhyML_Printf("\n. Likelihood - disk %s has %d lindisk nodes [%f]",devt->id,devt->n_ldsk_a,lnL); */

      lnL += log_lbda - mmod->lbda * (devt->time - devt->prev->time);

      was_hit = NO;
      n_hit   = 0;

      For(i,devt->n_ldsk_a)
        {
          lindisk_nd = devt->ldsk_a[i];
          was_hit    = lindisk_nd->prev->devt == devt->prev;

          if(was_hit == YES) n_hit++;

          if(MIGREP_Is_In_Disk(lindisk_nd->coord,devt->prev) == YES)
            {
              if(was_hit == YES) // was hit
                {
                  /* PhyML_Printf("\n. lindisk %s was hit and gave %s",lindisk_nd->coord->id,lindisk_nd->prev->coord->id); */
                  lnL += log_mu;
                }
              else // was not hit
                {
                  /* PhyML_Printf("\n. lindisk %s was not hit",lindisk_nd->coord->id); */
                  lnL += log_one_mu;
                }
            }
          else 
            { 
              // has changed spatial coordinates but was outside disk 
              if(was_hit == YES) 
                {
                  /* printf("\n. FAIL lindisk: %s [%.3f %.3f] prev: %s [%.3f %.3f] centr: [%.3f %.3f] rad: %.3f", */
                  /*        lindisk_nd->coord->id, */
                  /*        lindisk_nd->coord->lonlat[0], */
                  /*        lindisk_nd->coord->lonlat[1], */
                  /*        lindisk_nd->prev->coord->id, */
                  /*        lindisk_nd->prev->coord->lonlat[0], */
                  /*        lindisk_nd->prev->coord->lonlat[1], */
                  /*        devt->centr->lonlat[0], */
                  /*        devt->centr->lonlat[1], */
                  /*        devt->mmod->rad); */
                  /* fflush(NULL); */
                  mmod->c_lnL = UNLIKELY;
                  return UNLIKELY;
                }
            }
        }
      
      if(n_hit >= 2) // a coalescent event occurred
        {
          lnL -= log_pi_rad;
        }
    
      devt = devt->prev;

      if(!devt->prev) break;

    }
  while(1);
  
  mmod->c_lnL = lnL;
  return(lnL);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIGREP_MCMC(t_tree *tree)
{
  t_mcmc *mcmc;

  mcmc = MCMC_Make_MCMC_Struct();
  MCMC_Complete_MCMC(mcmc,tree);

  tree->mcmc = mcmc;

  mcmc->io               = NULL;
  mcmc->is               = NO;
  mcmc->use_data         = YES;
  mcmc->run              = 0;
  mcmc->sample_interval  = 1E+3;
  mcmc->chain_len        = 1E+6;
  mcmc->chain_len_burnin = 1E+5;
  mcmc->randomize        = YES;
  mcmc->norm_freq        = 1E+3;
  mcmc->max_tune         = 1.E+20;
  mcmc->min_tune         = 1.E-10;
  mcmc->print_every      = 2;
  mcmc->is_burnin        = NO;
  mcmc->nd_t_digits      = 1;
  mcmc->chain_len        = 1.E+8;
  mcmc->sample_interval  = 50;
  
  tree->mmod->lbda = Uni()*(tree->mmod->max_lbda - tree->mmod->min_lbda) + tree->mmod->min_lbda;
  tree->mmod->mu   = Uni()*(tree->mmod->max_mu - tree->mmod->min_mu) + tree->mmod->min_mu;
  tree->mmod->rad  = tree->mmod->max_rad;

  MIGREP_Lk(tree->devt,tree->mmod);
  printf("\n. LK: %f",tree->mmod->c_lnL);

  gtk_widget_queue_draw(tree->draw_area);
  sleep(2);

  PhyML_Printf("\n %10s %10s %10s %10s",
               "lnL",
               "lbda",
               "mu",
               "rad");
  mcmc->sample_interval = 100000;
  mcmc->run = 0;
  do
    {
      /* MCMC_Migrep_Lbda(tree); */
      /* MCMC_Migrep_Mu(tree); */
      /* MCMC_Migrep_Radius(tree); */
      MCMC_Migrep_Triplet(tree);
      /* MCMC_Migrep_Slice(tree); */
      Exit("\n");

      if(mcmc->run%mcmc->sample_interval == 0)
      /* if(mcmc->run == 0) */
        {
          PhyML_Printf("\n %10f %10f %10f %10f",
                       tree->mmod->c_lnL,
                       tree->mmod->lbda,
                       tree->mmod->mu,
                       tree->mmod->rad);

          /* gdk_threads_enter(); */
          gtk_widget_queue_draw(tree->draw_area);
          sleep(1);
          /* gdk_threads_leave(); */
        }

      mcmc->run++;
    }
  while(tree->mcmc->run < tree->mcmc->chain_len);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl MIGREP_Wrap_Lk(t_edge *b, t_tree *tree, supert_tree *stree)
{
  return MIGREP_Lk(tree->devt,tree->mmod);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIGREP_New_Traj(t_disk_evt *start, t_disk_evt *end, t_tree *tree)
{
  t_disk_evt *devt;
  int i,j;
  int n_hit_up,n_hit_tot;
  phydbl min_up, min_do;
  phydbl max_up, max_do;



  For(i,start->n_ldsk_a)
    {      
      /* printf("\n<><><>"); */

      devt = end;
      while(devt != start)
        {
          devt->ldsk_a[i]->min_coord = GEO_Make_Geo_Coord(devt->mmod->n_dim);
          devt->ldsk_a[i]->max_coord = GEO_Make_Geo_Coord(devt->mmod->n_dim);                    
          devt = devt->prev;
        }

      devt = end;
      n_hit_tot = 0;
      while(devt != start)
        {
          if(devt->ldsk_a[i]->is_hit == YES) n_hit_tot++;
          devt = devt->prev;
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
      /* devt = end; */
      /* while(devt != start) */
      /*   { */
      /*     printf("\n. %d %s %.2f %.2f", */
      /*            i, */
      /*            devt->ldsk_a[i]->coord->id, */
      /*            devt->ldsk_a[i]->coord->lonlat[0], */
      /*            devt->ldsk_a[i]->coord->lonlat[1]); */
      /*     devt = devt->prev; */
      /*   } */
      
      n_hit_up = 0;
      devt = end;
      while(devt->prev != start)
        {
          if(devt->ldsk_a[i]->is_hit == YES)
            {              
              n_hit_up++;

              For(j,devt->mmod->n_dim)
                {
                  min_up = end->ldsk_a[i]->coord->lonlat[j] - n_hit_up * 2. * devt->mmod->rad;
                  min_do = start->ldsk_a[i]->coord->lonlat[j] - (n_hit_tot - n_hit_up) * 2. * devt->mmod->rad;
                  devt->ldsk_a[i]->prev->min_coord->lonlat[j] = 
                    MAX(MAX(MAX(min_up,min_do),0.0),devt->prev->centr->lonlat[j] - devt->mmod->rad);

                  max_up = end->ldsk_a[i]->coord->lonlat[j] + n_hit_up * 2. * devt->mmod->rad;
                  max_do = start->ldsk_a[i]->coord->lonlat[j] + (n_hit_tot - n_hit_up) * 2. * devt->mmod->rad;
                  devt->ldsk_a[i]->prev->max_coord->lonlat[j] = 
                    MIN(MIN(MIN(max_up,max_do),devt->mmod->lim->lonlat[j]),devt->prev->centr->lonlat[j] + devt->mmod->rad);

                  /* printf("\n. curr: %s min_up: %.2f min_do: %.2f max_up: %.2f max_do: %.2f %d %d start: %.2f [%.2f %.2f]", */
                  /*        devt->ldsk_a[i]->prev->coord->id, */
                  /*        min_up,min_do,max_up,max_do, */
                  /*        n_hit_tot,n_hit_up, */
                  /*        start->ldsk_a[i]->coord->lonlat[j], */
                  /*        devt->ldsk_a[i]->prev->min_coord->lonlat[j], */
                  /*        devt->ldsk_a[i]->prev->max_coord->lonlat[j]); */
                }
            }
          devt = devt->prev;
        }
      
      devt = end;
      while(devt->prev != start)
        {
          if(devt->ldsk_a[i]->is_hit == YES)
            {
              For(j,devt->mmod->n_dim)
                devt->ldsk_a[i]->prev->coord->lonlat[j] = 
                Uni()*
                (devt->ldsk_a[i]->prev->max_coord->lonlat[j]  -
                 devt->ldsk_a[i]->prev->min_coord->lonlat[j]) +
                devt->ldsk_a[i]->prev->min_coord->lonlat[j];
            }
          else
            {
              For(j,devt->mmod->n_dim)
                devt->ldsk_a[i]->prev->coord->lonlat[j] = 
                devt->ldsk_a[i]->coord->lonlat[j];
            }

          devt = devt->prev;
        }
      
      devt = end;
      while(devt != start)
        {
          Free_Geo_Coord(devt->ldsk_a[i]->min_coord);
          Free_Geo_Coord(devt->ldsk_a[i]->max_coord);
          devt = devt->prev;
        }
    }

  /* gtk_widget_queue_draw(tree->draw_area); */
  
  /* For(i,start->n_ldsk_a) */
  /*   { */
  /*     printf("\n<><><>"); */
  /*     devt = end; */
  /*     while(devt != start) */
  /*       { */
  /*         printf("\nx %s %.2f %.2f centr: %.2f %.2f", */
  /*                devt->ldsk_a[i]->coord->id, */
  /*                devt->ldsk_a[i]->coord->lonlat[0], */
  /*                devt->ldsk_a[i]->coord->lonlat[1], */
  /*                devt->centr->lonlat[0], */
  /*                devt->centr->lonlat[1]); */
  /*         devt = devt->prev; */
  /*       } */
  /*     fflush(NULL); */
  /*   } */
  
  /* MIGREP_Lk(tree->devt,tree->mmod); */
  /* printf("\n. >> Lk: %f rad: %f",tree->mmod->c_lnL,tree->mmod->rad); */
  /* sleep(5); */
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIGREP_Copy_Coord(t_geo_coord *ori, t_geo_coord *cpy)
{
  int i;
  For(i,ori->dim) cpy->lonlat[i] = ori->lonlat[i];
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIGREP_Remove_Devt(t_disk_evt *devt)
{
  t_disk_evt *prev;
  t_disk_evt *next;

  prev = devt->prev;
  next = devt->next;

  if(prev == NULL)
    {
      PhyML_Printf("\n== Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
      Warn_And_Exit("");
    }

  if(next == NULL)
    {
      PhyML_Printf("\n== Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
      Warn_And_Exit("");
    }
  
  printf("\n. Remove disk %s",devt->id); fflush(NULL);
  
  prev->next = next;
  next->prev = prev;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// Insert disk event. devt->prev and devt->next need to be set
// accordingly.
void MIGREP_Insert_Devt(t_disk_evt *devt)
{
  if(devt->prev == NULL)
    {
      PhyML_Printf("\n== Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
      Warn_And_Exit("");
    }

  if(devt->next == NULL)
    {
      PhyML_Printf("\n== Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
      Warn_And_Exit("");
    }

  printf("\n. Insert disk %s @ %f",devt->id,devt->time); fflush(NULL);

  devt->prev->next = devt;
  devt->next->prev = devt;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

t_lindisk_nd *MIGREP_Prev_Coal_Lindisk(t_lindisk_nd *t)
{
  if(t == NULL)
    {
      PhyML_Printf("\n== Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
      Warn_And_Exit("");
    }

  if(t->is_coal == YES) 
    {
      return t;
    }
  else
    {
      return MIGREP_Prev_Coal_Lindisk(t->prev);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

t_lindisk_nd *MIGREP_Next_Coal_Lindisk(t_lindisk_nd *t)
{
  if(t == NULL)
    {
      PhyML_Printf("\n== Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
      Warn_And_Exit("");
    }

  if(t->is_coal == YES || t->next == NULL) return t;
  else
    {
      if(t->n_next > 1) // Should have t->is_coal = YES
        {
          PhyML_Printf("\n== Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
          Warn_And_Exit("");
        }

      return MIGREP_Next_Coal_Lindisk(t->next[0]);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

// Generate a new trajectory, including disk event centers, between
// 'y_ldsk' a ``young'' lindisk event and 'o_ldsk' an old one. 'y_ldsk
// and 'o_ldsk' remain unaffected.
void MIGREP_One_New_Traj(t_lindisk_nd *y_ldsk, t_lindisk_nd *o_ldsk)
{
  int n_disk_btw;
  t_lindisk_nd *ldsk;
  phydbl min, max;
  int i;
  phydbl rad;

  // Number of disks between y_ldsk and o_ldsk
  ldsk = y_ldsk;
  n_disk_btw = 0;
  while(ldsk->prev != o_ldsk)
    {
      n_disk_btw++;
      ldsk = ldsk->prev;
    }

  printf("\n. Number of disks between %s and %s: %d",y_ldsk->coord->id,o_ldsk->coord->id,n_disk_btw);
  fflush(NULL);

  ldsk = y_ldsk;
  rad  = ldsk->devt->mmod->rad;
  while(ldsk->prev != o_ldsk)
    {
      For(i,ldsk->devt->mmod->n_dim)
        {
          min = 
            MIN(0,
                MIN(ldsk->coord->lonlat[i] - 2.*rad,
                    o_ldsk->coord->lonlat[i] - 2.*rad*n_disk_btw));

          max = 
            MAX(ldsk->devt->mmod->lim->lonlat[i],
                MAX(ldsk->coord->lonlat[i] + 2.*rad,
                    o_ldsk->coord->lonlat[i] + 2.*rad*n_disk_btw));
          
          printf("\n. ldsk %s at %f min: %f max: %f",ldsk->coord->id,ldsk->coord->lonlat[i],min,max);

          if(max < min)
            {
              PhyML_Printf("\n== Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
              Warn_And_Exit("");              
            }
                
          // New coordinate for the lindisk
          ldsk->prev->coord->lonlat[i] = Uni()*(max - min) + min;

          printf("\n. NEW TRAJ ldsk %s at %f [devt: %s] located at %f",
                 ldsk->prev->coord->id,
                 ldsk->prev->devt->time,
                 ldsk->prev->devt->id,
                 ldsk->prev->coord->lonlat[0]);
          fflush(NULL);

          // New coordinate for the centre of the corresponding disk event
          max = ldsk->prev->coord->lonlat[i] + rad;
          min = ldsk->prev->coord->lonlat[i] - rad;
          ldsk->prev->devt->centr->lonlat[i] = Uni()*(max - min) + min;
        }
      ldsk = ldsk->prev;
      n_disk_btw--;
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

// Return the index of the 'next' element of 'old' that should be
// used in order to reach 'young'.
int MIGREP_Get_Next_Direction(t_lindisk_nd *young, t_lindisk_nd *old)
{
  if(young->devt->time < old->devt->time)
    {
      PhyML_Printf("\n== Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
      Warn_And_Exit("");
    }

  if(young == NULL)
    {
      PhyML_Printf("\n== Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
      Warn_And_Exit("");
    }

  if(young->prev == old)
    {
      int i;
      For(i,old->n_next) if(old->next[i] == young) return i;
    }
  else
    {
      return MIGREP_Get_Next_Direction(young->prev,old);
    }
  return(-1);
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
