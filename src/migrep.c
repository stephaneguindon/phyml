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
  t_dsk *disk;
  t_migrep_mod *mmod;
  int i;

  tree = (t_tree *)(data);
  
  mmod = tree->mmod;

  disk = tree->disk;
  while(disk->prev) disk = disk->prev;

  h = FABS(disk->time);
  w = tree->mmod->lim->lonlat[1];

  cairo_select_font_face(cr, "Courier",
                         CAIRO_FONT_SLANT_NORMAL,
                         CAIRO_FONT_WEIGHT_NORMAL);
  
  cairo_set_font_size(cr,0.03);
  cairo_scale (cr,WINDOW_WIDTH/1.2,WINDOW_HEIGHT/1.2);
  cairo_translate(cr,0.05,0.05);
  cairo_set_source_rgb(cr,0,0,0);
  cairo_set_line_width(cr,0.002);

  while(disk->next) disk = disk->next;
  do
    {      
      MIGREP_Update_Lindisk_List(disk->time,disk->ldsk_a,&(disk->n_ldsk_a),disk);

      cairo_move_to(cr,0,FABS(disk->time)/h);
      cairo_show_text(cr,disk->id); 
      cairo_stroke(cr);

      cairo_move_to(cr,disk->centr->lonlat[1]/w - mmod->rad/w,FABS(disk->time)/h);
      cairo_set_line_width(cr,0.004);
      cairo_set_source_rgba(cr, 1, 0, 0, 0.3);
      cairo_set_line_width(cr,0.002);
      cairo_line_to(cr,disk->centr->lonlat[1]/w + mmod->rad/w,FABS(disk->time)/h);
      cairo_stroke(cr);
      cairo_arc(cr,disk->centr->lonlat[1]/w,FABS(disk->time)/h,0.005,0.,2.*PI);
      cairo_stroke(cr);
      cairo_set_source_rgb(cr, 0, 0, 0);
          
      For(i,disk->n_ldsk_a)
        {
          cairo_move_to(cr,disk->ldsk_a[i]->coord->lonlat[1]/w,FABS(disk->time)/h);
          /* cairo_show_text(cr,disk->ldsk_a[i]->coord->id);  */
          cairo_stroke(cr);
          cairo_move_to(cr,disk->ldsk_a[i]->coord->lonlat[1]/w,FABS(disk->time)/h);
          cairo_line_to(cr,disk->ldsk_a[i]->coord->lonlat[1]/w,FABS(disk->prev->time)/h);
          cairo_stroke(cr);

          if(disk->ldsk_a[i]->prev->disk == disk->prev)
            {
              cairo_line_to(cr,disk->ldsk_a[i]->coord->lonlat[1]/w,FABS(disk->prev->time)/h);
              cairo_line_to(cr,disk->ldsk_a[i]->prev->coord->lonlat[1]/w,FABS(disk->prev->time)/h);
              cairo_stroke(cr);
            }

          if(disk->ldsk_a[i]->is_coal && disk->ldsk_a[i]->disk == disk) 
            {
              cairo_set_source_rgba(cr, 0, 0.2, 0.8, 0.5);
              cairo_arc(cr,disk->ldsk_a[i]->coord->lonlat[1]/w,FABS(disk->time)/h,0.01,0.,2.*PI);
              cairo_fill(cr);
              cairo_stroke(cr);
              cairo_set_source_rgb(cr, 0, 0, 0);
            }
          else
            {
              cairo_set_source_rgba(cr, 0.2, 0.2, 0.2, 0.8);
              cairo_arc(cr,disk->ldsk_a[i]->coord->lonlat[1]/w,FABS(disk->time)/h,0.005,0.,2.*PI);
              cairo_fill(cr);
              cairo_stroke(cr);
              cairo_set_source_rgb(cr, 0, 0, 0);
            }
        }
      disk = disk->prev;
      if(disk->prev == NULL) break;
    }while(1);

  // Root disk
  cairo_move_to(cr,disk->ldsk->coord->lonlat[1]/w,FABS(disk->time)/h);
  cairo_show_text(cr,disk->ldsk->coord->id); 
  cairo_arc(cr,disk->ldsk->coord->lonlat[1]/w,FABS(disk->time)/h,0.005,0.,2.*PI);
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
  /* seed = 1411708990; */
  /* seed = 1411709351; */
  /* seed = 1412025046; */
  /* seed = 1412197537; */
  seed = 1412214370;
  printf("\n. Seed: %d",seed);
  srand(seed);
  tree = MIGREP_Simulate_Backward((int)atoi(argv[1]),10.,10.);

  MIGREP_MCMC(tree);
  Exit("\n");

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
  int n_dim,n_lineages,n_lineages_new,n_disk,n_hit;
  phydbl curr_t,dt_dsk;
  t_migrep_mod *mmod;
  t_dsk *disk;
  t_ldsk **ldsk_a,**ldsk_a_tmp,*new_ldsk;
  t_node *nd;

  n_dim = 2; // 2-dimensional landscape

  tree = Make_Tree_From_Scratch(n_otu,NULL);
  nd   = Make_Node_Light(0);
  disk = MIGREP_Make_Disk_Event(n_dim,n_otu);
  MIGREP_Init_Disk_Event(disk,NULL);
  
  // Allocate coordinates for all the tips first (will grow afterwards)
  ldsk_a = (t_ldsk **)mCalloc(n_otu,sizeof(t_ldsk *));
  For(i,n_otu) 
    {
      ldsk_a[i] = MIGREP_Make_Lindisk_Node(n_dim);
      MIGREP_Init_Lindisk_Node(ldsk_a[i],disk,n_dim);
    }

  ldsk_a_tmp = (t_ldsk **)mCalloc(n_otu,sizeof(t_ldsk *));

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
  disk->time             = 0.0;
  disk->mmod             = mmod;
  disk->centr->lonlat[0] = .5*width;
  disk->centr->lonlat[1] = .5*height;      
  
  // Allocate and initialise for next event
  disk->prev = MIGREP_Make_Disk_Event(n_dim,n_otu);
  MIGREP_Init_Disk_Event(disk->prev,NULL);
  disk->prev->next = disk;
  
  For(i,n_otu)
    {
      printf("\nx disk %s [%15f] %3d %s %15f %15f",
             disk?disk->id:"xxx",
             disk?disk->time:0.0,
             i,
             ldsk_a[i]->coord->id,
             ldsk_a[i]->coord->lonlat[0],
             ldsk_a[i]->coord->lonlat[1]);
    }

  // Move to it
  disk = disk->prev;

  // Initialize parameters of migrep model
  mmod->lbda = 0.2;
  mmod->mu   = 1.0;
  mmod->rad  = 3.0;
  
  curr_t      = 0.0;
  dt_dsk     = 0.0;
  n_lineages  = n_otu;
  n_disk      = 0;
  do
    {
      // Time of next event
      dt_dsk = Rexp(mmod->lbda);
      curr_t -= dt_dsk;
      
      // Coordinates of next event
      GEO_Init_Coord(disk->centr,n_dim);
      disk->centr->lonlat[0] = Uni()*width;
      disk->centr->lonlat[1] = Uni()*height;      

      disk->time = curr_t;
      disk->mmod = mmod;

      /* printf("\n. Disk %s has %d lindisk nodes and %d disks under",disk->id,disk->n_ldsk_a,disk->n_disk_under); */

      // New lindisk (will not be used is no lineage is hit) 
      new_ldsk = MIGREP_Make_Lindisk_Node(n_dim);
      MIGREP_Init_Lindisk_Node(new_ldsk,disk,n_dim);
      do
        {
          Runif_Disk(new_ldsk->coord->lonlat,
                     new_ldsk->coord->lonlat+1,
                     disk->centr->lonlat[0],
                     disk->centr->lonlat[1],
                     disk->mmod->rad);
          // Check that the next lindisk is within boundaries
          if(new_ldsk->coord->lonlat[0] > 0.0 &&
             new_ldsk->coord->lonlat[0] < disk->mmod->lim->lonlat[0] &&
             new_ldsk->coord->lonlat[1] > 0.0 &&
             new_ldsk->coord->lonlat[1] < disk->mmod->lim->lonlat[0]) break;
        }
      while(1);

      n_hit          = 0;
      n_lineages_new = 0;
      For(i,n_lineages)
        {
          ldsk_a[i]->prev = NULL;

          if(MIGREP_Is_In_Disk(ldsk_a[i]->coord,disk) == YES) 
            {
              if(Uni() < mmod->mu)
                {
                  printf("\n. Hit and die %s @ time %f. Center: %f %f Go to %f %f",
                         ldsk_a[i]->coord->id,
                         disk->time,
                         disk->centr->lonlat[0],
                         disk->centr->lonlat[1],
                         new_ldsk->coord->lonlat[0],
                         new_ldsk->coord->lonlat[1]);
                  
                  MIGREP_Make_Lindisk_Next(new_ldsk);

                  ldsk_a[i]->prev                    = new_ldsk;
                  ldsk_a[i]->prev->is_hit            = YES;
                  new_ldsk->next[new_ldsk->n_next-1] = ldsk_a[i]; 
                  disk->ldsk                         = new_ldsk;

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
          disk->nd          = nd;
          new_ldsk->is_coal = YES;
          PhyML_Printf("\n. Coalescent @ disk %s on ldsk %s",disk->id,disk->ldsk->coord->id);
        }
 
      ldsk_a = (t_ldsk **)mCalloc(n_otu,sizeof(t_ldsk *));
      For(i,n_lineages) ldsk_a[i] = ldsk_a_tmp[i];
      
      For(i,n_lineages)
        {
          printf("\n. disk %s [%15f] %3d %s %15f %15f",
                 disk?disk->id:"xxx",
                 disk?disk->time:0.0,
                 i,
                 ldsk_a[i]->coord->id,
                 ldsk_a[i]->coord->lonlat[0],
                 ldsk_a[i]->coord->lonlat[1]);
        }

      n_disk++;


      if(n_lineages == 1) break;

      disk->prev = MIGREP_Make_Disk_Event(n_dim,n_otu);
      MIGREP_Init_Disk_Event(disk->prev,NULL);
      disk->prev->next = disk;
      
      disk = disk->prev;          
    }
  while(1);

  while(disk->next) disk = disk->next;

  tree->disk = disk;
  tree->mmod = mmod;


  /* while(disk->prev) */
  /*   { */
  /*     printf("\n<><><>"); */
  /*     printf("\n disk %f %s",disk->time,disk->id); */
  /*     For(i,disk->n_ldsk_a) */
  /*       { */
  /*         printf("\n. %s %f %f", */
  /*            disk->ldsk_a[i]->coord->id, */
  /*            disk->ldsk_a[i]->coord->lonlat[0], */
  /*            disk->ldsk_a[i]->coord->lonlat[1]); */
  /*       } */
  /*     disk = disk->prev; */
  /*   } */


  MIGREP_Lk(disk,mmod);
  printf("\n. LK: %f",mmod->c_lnL);
  /* Exit("\n"); */
  /* MIGREP_MCMC(tree); */

  return(tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int MIGREP_Is_In_Disk(t_geo_coord *coord, t_dsk *disk)
{
  int i;

  For(i,disk->centr->dim)
    {
      if(FABS(coord->lonlat[i] - disk->centr->lonlat[i]) > disk->mmod->rad)
        {
          return(NO);
        }
    }
  return(YES);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl MIGREP_Lk(t_dsk *disk, t_migrep_mod *mmod)
{
  phydbl lnL;
  phydbl log_lbda,log_pi_rad,log_mu,log_one_mu;
  t_ldsk *lindisk_nd;
  int i;
  short int was_hit, n_hit;

  // Rewind if necessary
  while(disk->next) { disk = disk->next; }  

  if(!disk->prev) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);    

  log_lbda   = LOG(mmod->lbda);
  log_mu     = LOG(mmod->mu);
  log_pi_rad = LOG(PI*mmod->rad*mmod->rad);
  log_one_mu = LOG(1. - mmod->mu);  
  lnL        = 0.0;

  do
    {
      /* PhyML_Printf("\n. Likelihood - disk %s has %d lindisk nodes [%f]",disk->id,disk->n_ldsk_a,lnL); */

      lnL += log_lbda - mmod->lbda * (disk->time - disk->prev->time);

      MIGREP_Update_Lindisk_List(disk->time,disk->ldsk_a,&(disk->n_ldsk_a),disk);

      was_hit = NO;
      n_hit   = 0;

      For(i,disk->n_ldsk_a)
        {
          lindisk_nd = disk->ldsk_a[i];
          was_hit    = lindisk_nd->prev->disk == disk->prev;

          if(was_hit == YES) n_hit++;

          if(MIGREP_Is_In_Disk(lindisk_nd->coord,disk->prev) == YES)
            {
              if(was_hit == YES) // was hit
                {
                  /* PhyML_Printf("\n. %d/%d lindisk %s was hit and gave %s",i,disk->n_ldsk_a,lindisk_nd->coord->id,lindisk_nd->prev->coord->id); */
                  lnL += log_mu;
                }
              else // was not hit
                {
                  /* PhyML_Printf("\n. %d/%d lindisk %s was not hit",i,disk->n_ldsk_a,lindisk_nd->coord->id); */
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
                  /*        disk->centr->lonlat[0], */
                  /*        disk->centr->lonlat[1], */
                  /*        disk->mmod->rad); */
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
    
      disk = disk->prev;

      if(!disk->prev) break;

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

  MIGREP_Lk(tree->disk,tree->mmod);
  printf("\n. LK: %f",tree->mmod->c_lnL);

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

      /* gtk_widget_queue_draw(tree->draw_area); */
      /* sleep(5); */

      MCMC_Migrep_Triplet(tree);

      /* MCMC_Migrep_Slice(tree); */
      /* Exit("\n"); */

      if(mcmc->run%mcmc->sample_interval == 0)
      /* if(mcmc->run == 0) */
        {
          PhyML_Printf("\n %10f %10f %10f %10f",
                       tree->mmod->c_lnL,
                       tree->mmod->lbda,
                       tree->mmod->mu,
                       tree->mmod->rad);

          /* gdk_threads_enter(); */

          /* gtk_widget_queue_draw(tree->draw_area); */
          /* sleep(1); */

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
  return MIGREP_Lk(tree->disk,tree->mmod);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIGREP_New_Traj(t_dsk *start, t_dsk *end, t_tree *tree)
{
  t_dsk *disk;
  int i,j;
  int n_hit_up,n_hit_tot;
  phydbl min_up, min_do;
  phydbl max_up, max_do;



  For(i,start->n_ldsk_a)
    {      
      /* printf("\n<><><>"); */

      disk = end;
      while(disk != start)
        {
          disk->ldsk_a[i]->min_coord = GEO_Make_Geo_Coord(disk->mmod->n_dim);
          disk->ldsk_a[i]->max_coord = GEO_Make_Geo_Coord(disk->mmod->n_dim);                    
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

              For(j,disk->mmod->n_dim)
                {
                  min_up = end->ldsk_a[i]->coord->lonlat[j] - n_hit_up * 2. * disk->mmod->rad;
                  min_do = start->ldsk_a[i]->coord->lonlat[j] - (n_hit_tot - n_hit_up) * 2. * disk->mmod->rad;
                  disk->ldsk_a[i]->prev->min_coord->lonlat[j] = 
                    MAX(MAX(MAX(min_up,min_do),0.0),disk->prev->centr->lonlat[j] - disk->mmod->rad);

                  max_up = end->ldsk_a[i]->coord->lonlat[j] + n_hit_up * 2. * disk->mmod->rad;
                  max_do = start->ldsk_a[i]->coord->lonlat[j] + (n_hit_tot - n_hit_up) * 2. * disk->mmod->rad;
                  disk->ldsk_a[i]->prev->max_coord->lonlat[j] = 
                    MIN(MIN(MIN(max_up,max_do),disk->mmod->lim->lonlat[j]),disk->prev->centr->lonlat[j] + disk->mmod->rad);

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
              For(j,disk->mmod->n_dim)
                disk->ldsk_a[i]->prev->coord->lonlat[j] = 
                Uni()*
                (disk->ldsk_a[i]->prev->max_coord->lonlat[j]  -
                 disk->ldsk_a[i]->prev->min_coord->lonlat[j]) +
                disk->ldsk_a[i]->prev->min_coord->lonlat[j];
            }
          else
            {
              For(j,disk->mmod->n_dim)
                disk->ldsk_a[i]->prev->coord->lonlat[j] = 
                disk->ldsk_a[i]->coord->lonlat[j];
            }

          disk = disk->prev;
        }
      
      disk = end;
      while(disk != start)
        {
          Free_Geo_Coord(disk->ldsk_a[i]->min_coord);
          Free_Geo_Coord(disk->ldsk_a[i]->max_coord);
          disk = disk->prev;
        }
    }

  /* gtk_widget_queue_draw(tree->draw_area); */
  
  /* For(i,start->n_ldsk_a) */
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
  
  /* MIGREP_Lk(tree->disk,tree->mmod); */
  /* printf("\n. >> Lk: %f rad: %f",tree->mmod->c_lnL,tree->mmod->rad); */
  /* sleep(5); */
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void MIGREP_Copy_Coord(t_geo_coord *ori, t_geo_coord *cpy)
{
  int i;
  For(i,ori->dim) cpy->lonlat[i] = ori->lonlat[i];
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void MIGREP_Remove_Disk(t_dsk *disk)
{
  t_dsk *prev;
  t_dsk *next;

  prev = disk->prev;
  next = disk->next;

  /* printf("\n. Remove disk %s (prev: %s next: %s)", */
  /*        disk->id, */
  /*        prev?prev->id:NULL, */
  /*        next?next->id:NULL); fflush(NULL); */

  if(prev == NULL) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  if(next == NULL) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);    
  
  
  prev->next = next;
  next->prev = prev;

  /* printf(" set %s->next to %s and %s->prev to %s", */
  /*        prev->id,prev->next->id, */
  /*        next->id,next->prev->id); fflush(NULL); */
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

/* Insert disk event. disk->prev and disk->next need to be set 
   accordingly. 
*/
void MIGREP_Insert_Disk(t_dsk *disk)
{

  if(disk->prev == NULL) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);    
  if(disk->next == NULL) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);    
  if(disk->prev->time > disk->time) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  if(disk->next->time < disk->time) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);

  /* printf("\n. Insert disk %s @ %f between %s and %s", */
  /*        disk->id,disk->time, */
  /*        disk->prev->id, */
  /*        disk->next->id); */
  /* fflush(NULL); */

  disk->prev->next = disk;
  disk->next->prev = disk;
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

t_ldsk *MIGREP_Prev_Coal_Lindisk(t_ldsk *t)
{
  if(t == NULL) return NULL;

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

t_ldsk *MIGREP_Next_Coal_Lindisk(t_ldsk *t)
{
  if(t == NULL) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);    

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

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

/*  Generate a new trajectory, including disk event centers, between  
    'y_ldsk' a ``young'' lindisk event and 'o_ldsk' an old one. 'y_ldsk 
    and 'o_ldsk' remain unaffected. No disk events should be present    
    between y_ldsk and o_ldsk: we need to generate some first.          
*/
void MIGREP_One_New_Traj(t_ldsk *y_ldsk, t_ldsk *o_ldsk, int dir_o_y, t_tree *tree)
{
  t_migrep_mod *mmod;
  t_dsk *disk,**disk_new;
  phydbl max_dist;
  int i;
  int min_n_disk,n_new_disk;

  mmod     = tree->mmod;
  disk     = NULL;
  disk_new = NULL;

  /* printf("\n. New traj from %s to %s",y_ldsk->coord->id,o_ldsk->coord->id); */
  /* fflush(NULL); */

  /* Minimum number of disks between y_ldsk and o_ldsk */
  max_dist = -1.;
  For(i,mmod->n_dim)
    {
      if(FABS(y_ldsk->coord->lonlat[i] - o_ldsk->coord->lonlat[i]) > max_dist)
        max_dist = FABS(y_ldsk->coord->lonlat[i] - o_ldsk->coord->lonlat[i]);
    }
  min_n_disk = (int)(max_dist / (2. * mmod->rad));
  
  /* printf("\n. min_n_disk: %d [%f]",min_n_disk,max_dist); */
  /* fflush(NULL); */
  
  /* How many disks along the new path between y_ldsk and o_ldsk */
  n_new_disk = Rand_Int(min_n_disk,min_n_disk+5);
  /* printf("\n. n_new_disk: %d",n_new_disk); fflush(NULL); */
  
  if(n_new_disk > 0)
    {
      /* Make new disks to create a new path between ldsk_left and ldsk_up */
      disk_new = (t_dsk **)mCalloc(n_new_disk,sizeof(t_dsk *));
      For(i,n_new_disk) disk_new[i] = MIGREP_Make_Disk_Event(mmod->n_dim,tree->n_otu);
      For(i,n_new_disk) MIGREP_Init_Disk_Event(disk_new[i],mmod);
      
      /* Times of these new disks */
      For(i,n_new_disk)
        disk_new[i]->time =
        Uni()*(y_ldsk->disk->time - o_ldsk->disk->time) + o_ldsk->disk->time;
      
      disk = tree->disk;
      
      /* Insert these events */
      For(i,n_new_disk)
        {
          if(tree->disk->next) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
          disk = tree->disk;
          while(disk->time > disk_new[i]->time) disk = disk->prev;
          disk_new[i]->prev = disk;
          disk_new[i]->next = disk->next;
          MIGREP_Insert_Disk(disk_new[i]);
        }
            
      /* For(i,n_new_disk) */
      /*   { */
      /*     printf("\n> disk_new: %f [%s]",disk_new[i]->time,disk_new[i]->id); fflush(NULL); */
      /*   } */
      
      /* Add new lindisks to the new disk events */
      For(i,n_new_disk)
        {
          disk_new[i]->ldsk = MIGREP_Make_Lindisk_Node(tree->mmod->n_dim);
          MIGREP_Init_Lindisk_Node(disk_new[i]->ldsk,disk_new[i],tree->mmod->n_dim);
          MIGREP_Make_Lindisk_Next(disk_new[i]->ldsk);
          /* printf("\n. Add ldsk %s to %s",disk_new[i]->ldsk->coord->id,disk_new[i]->id); fflush(NULL); */
        }
      
      /* Connect them */
      MIGREP_Connect_Ldsk_Given_Disk(disk_new,n_new_disk,y_ldsk,o_ldsk,dir_o_y);

      Free(disk_new);
    }
  else
    {
      o_ldsk->next[dir_o_y] = y_ldsk;
      y_ldsk->prev          = o_ldsk;
    }
 
  /* Generate a trajectory */
  MIGREP_One_New_Traj_Given_Disk(y_ldsk,o_ldsk);  
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

/* Generate a new trajectory, including disk event centers, between 
  'y_ldsk' a ``young'' lindisk event and 'o_ldsk' an old one. 'y_ldsk 
  and 'o_ldsk' remain unaffected. Disk events between these two ldsk 
  should already be set. 
*/
void MIGREP_One_New_Traj_Given_Disk(t_ldsk *y_ldsk, t_ldsk *o_ldsk)
{
  int n_disk_btw;
  t_ldsk *ldsk;
  phydbl min, max;
  int i;
  phydbl rad;

  /* Number of disks between y_ldsk and o_ldsk */
  ldsk = y_ldsk;
  n_disk_btw = 0;
  while(ldsk->prev != o_ldsk)
    {
      n_disk_btw++;
      ldsk = ldsk->prev;
    }

  /* printf("\n. Number of disks between %s and %s: %d",y_ldsk->coord->id,o_ldsk->coord->id,n_disk_btw); */
  /* fflush(NULL); */

  ldsk = y_ldsk;
  rad  = ldsk->disk->mmod->rad;
  while(ldsk->prev != o_ldsk)
    {
      /* printf("\n. ldsk %s at %f disk: %s ",ldsk->coord->id,ldsk->disk->time,ldsk->disk->id); */
      /* fflush(NULL); */

      For(i,ldsk->disk->mmod->n_dim)
        {
          min = 
            MAX(0,
                MAX(ldsk->coord->lonlat[i] - 2.*rad,
                    o_ldsk->coord->lonlat[i] - 2.*rad*n_disk_btw));

          max = 
            MIN(ldsk->disk->mmod->lim->lonlat[i],
                MIN(ldsk->coord->lonlat[i] + 2.*rad,
                    o_ldsk->coord->lonlat[i] + 2.*rad*n_disk_btw));
          

          if(max < min) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
                
          // New coordinate for the lindisk
          ldsk->prev->coord->lonlat[i] = Uni()*(max - min) + min;

          /* printf("\n. NEW TRAJ ldsk %s at %f [disk: %s] located at %f", */
          /*        ldsk->prev->coord->id, */
          /*        ldsk->prev->disk->time, */
          /*        ldsk->prev->disk->id, */
          /*        ldsk->prev->coord->lonlat[i]); */
          /* fflush(NULL); */

          // New coordinate for the centre of the corresponding disk event
          max = ldsk->prev->coord->lonlat[i] + rad;
          min = ldsk->prev->coord->lonlat[i] - rad;
          ldsk->prev->disk->centr->lonlat[i] = Uni()*(max - min) + min;
        }
      ldsk = ldsk->prev;
      n_disk_btw--;
    }
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

/* Return the index of the 'next' element of 'old' that should be
   used in order to reach 'young'. 
*/
int MIGREP_Get_Next_Direction(t_ldsk *young, t_ldsk *old)
{
  if(young->disk->time < old->disk->time)
    {
      PhyML_Printf("\n== young (%s) @ time %f; old (%s) @ time %f",
                   young->coord->id,young->disk->time,
                   old->coord->id,old->disk->time);
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);    
    }

  if(young == NULL) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);    

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

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void MIGREP_Update_Lindisk_List(phydbl time, t_ldsk **list, int *pos, t_dsk *disk)
{
  t_dsk *root_dsk;

  *pos = 0;
  root_dsk = disk;
  while(root_dsk->prev) root_dsk = root_dsk->prev;
  /* printf("\n. root_dsk: %s",root_dsk?root_dsk->id:"xx"); */
  MIGREP_Update_Lindisk_List_Pre(root_dsk->ldsk,time,list,pos);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void MIGREP_Update_Lindisk_List_Pre(t_ldsk *ldsk, phydbl time, t_ldsk **list, int *pos)
{
  /* printf("\n. time: %f pos: %d ldsk: %s disk: %s n_next: %d ldsk->disk->time: %f disk: %s", */
  /*        time, */
  /*        *pos, */
  /*        ldsk ? ldsk->coord->id : "xx", */
  /*        ldsk ? ldsk->disk->id : "zz", */
  /*        ldsk ? ldsk->n_next : -1, */
  /*        ldsk ? ldsk->disk->time : -1., */
  /*        ldsk ? ldsk->disk->id : "yy"); fflush(NULL); */

  if(ldsk == NULL) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);

  if((ldsk->prev != NULL) && (ldsk->disk->time > time) && (ldsk->prev->disk->time < time))
    {
      list[*pos] = ldsk;
      *pos = *pos + 1;
    }
  else if(Are_Equal(ldsk->disk->time,time,SMALL_DBL))
    {
      list[*pos] = ldsk;
      *pos = *pos + 1;      
    }
  else if(ldsk->disk->time < time)
    {
      int i;
      For(i,ldsk->n_next)
        MIGREP_Update_Lindisk_List_Pre(ldsk->next[i],time,list,pos);    
    }
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

/* Connect all the ldsk between y_ldsk (young ldsk) and o_ldsk (old ldsk).
   The disk between y_ldsk and o_ldsk should have all been set already
   Note: the disks in **disk are sorted in ascending order of their 
   times
*/

void MIGREP_Connect_Ldsk_Given_Disk(t_dsk **disk, int n_disk, t_ldsk *y_ldsk, t_ldsk *o_ldsk, int dir_o_y)
{
  int i,j;
  t_dsk *disk_tmp;

  /* Sort these events by ascending order of their times */
  For(i,n_disk-1)
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


  For(i,n_disk)
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

void MIGREP_Print_Struct(char sign, t_tree *tree)
{
  t_dsk *disk;
  int i;
  t_ldsk *ldisk;

  disk = tree->disk;
  while(disk->prev) disk = disk->prev;
  do
    {
      printf("\n%c Disk: %s @ %7.3f has %3s on it is_coal? %2d",
             sign,
             disk->id,
             disk->time,disk->ldsk?disk->ldsk->coord->id:NULL,
             disk->ldsk?disk->ldsk->is_coal:-1); fflush(NULL);

      MIGREP_Update_Lindisk_List(disk->time,disk->ldsk_a,&(disk->n_ldsk_a),disk);

      For(i,disk->n_ldsk_a)
        {
          ldisk = disk->ldsk_a[i];
          printf("\n%c ldisk: %s prev: %s",
                 sign,
                 ldisk->coord->id,
                 ldisk->prev ? ldisk->prev->coord->id : NULL);
        }

      disk = disk->next;      
    }
  while(disk);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
