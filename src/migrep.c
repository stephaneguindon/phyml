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
  int i;

  tree = (t_tree *)(data);
  
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
      /* MIGREP_Update_Lindisk_List(disk->time,disk->ldsk_a,&(disk->n_ldsk_a),disk); */
      MIGREP_Update_Lindisk_List(tree);

      /* cairo_move_to(cr,0,FABS(disk->time)/h); */
      /* cairo_show_text(cr,disk->id); */
      /* cairo_stroke(cr); */

      /* cairo_move_to(cr,disk->centr->lonlat[1]/w - mmod->rad/w,FABS(disk->time)/h); */
      /* cairo_set_line_width(cr,0.004); */
      /* cairo_set_source_rgba(cr, 1, 0, 0, 0.3); */
      /* cairo_set_line_width(cr,0.002); */
      /* cairo_line_to(cr,disk->centr->lonlat[1]/w + mmod->rad/w,FABS(disk->time)/h); */
      /* cairo_stroke(cr); */
      /* cairo_arc(cr,disk->centr->lonlat[1]/w,FABS(disk->time)/h,0.005,0.,2.*PI); */
      /* cairo_stroke(cr); */
      /* cairo_set_source_rgb(cr, 0, 0, 0); */
          
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
  /* seed = 1413426621; */
  /* seed = 1413492392; */
  /* seed = 1413497813; */
  /* seed = 1414890813; */
  printf("\n# Seed: %d",seed); fflush(NULL);
  /* seed = 1412394873; */
  srand(seed);
  tree = MIGREP_Simulate_Backward((int)atoi(argv[1]),10.,10.);

  MIGREP_MCMC(tree);
  Exit("\n");

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
  MIGREP_Init_Disk_Event(disk,n_dim,NULL);
  
  // Allocate coordinates for all the tips first (will grow afterwards)
  ldsk_a = (t_ldsk **)mCalloc(n_otu,sizeof(t_ldsk *));
  For(i,n_otu) 
    {
      ldsk_a[i] = MIGREP_Make_Lindisk_Node(n_dim);
      MIGREP_Init_Lindisk_Node(ldsk_a[i],disk,n_dim);
    }

  ldsk_a_tmp = (t_ldsk **)mCalloc(n_otu,sizeof(t_ldsk *));

  /* Generate coordinates for the tip nodes (uniform distribution on the rectangle) */
  For(i,n_otu)
    {
      ldsk_a[i]->coord->lonlat[0] = Uni()*width;  // longitude
      ldsk_a[i]->coord->lonlat[1] = Uni()*height; // latitude
      /* ldsk_a[i]->coord->lonlat[0] = (i+1.)/(n_otu+1.)*width;  // longitude */
      /* ldsk_a[i]->coord->lonlat[1] = (i+1.)/(n_otu+1.)*height; // latitude */
    }

  /* Allocate migrep model */
  mmod = MIGREP_Make_Migrep_Model(n_dim);
  MIGREP_Init_Migrep_Mod(mmod,n_dim);
  mmod->lim->lonlat[0] = width;
  mmod->lim->lonlat[1] = height;
  
  /* First disk event (at time 0) */
  disk->time             = 0.0;
  disk->mmod             = mmod;
  disk->centr->lonlat[0] = .5*width;
  disk->centr->lonlat[1] = .5*height;      
  disk->ldsk_a           = ldsk_a;
  disk->n_ldsk_a         = n_otu;

  /* Allocate and initialise for next event */
  disk->prev = MIGREP_Make_Disk_Event(n_dim,n_otu);
  MIGREP_Init_Disk_Event(disk->prev,n_dim,NULL);
  disk->prev->next = disk;
  
  /* For(i,n_otu) */
  /*   { */
  /*     printf("\nx disk %s [%15f] %3d %s %15f %15f", */
  /*            disk?disk->id:"xxx", */
  /*            disk?disk->time:0.0, */
  /*            i, */
  /*            ldsk_a[i]->coord->id, */
  /*            ldsk_a[i]->coord->lonlat[0], */
  /*            ldsk_a[i]->coord->lonlat[1]); */
  /*   } */

  /* Move to it */
  disk = disk->prev;

  /* Initialize parameters of migrep model */
  mmod->lbda = 0.3;
  mmod->mu   = 0.8;
  mmod->rad  = 3.0;
  
  curr_t      = 0.0;
  dt_dsk     = 0.0;
  n_lineages  = n_otu;
  n_disk      = 0;
  do
    {
      /* Time of next event */
      dt_dsk = Rexp(mmod->lbda);
      curr_t -= dt_dsk;
      
      /* Coordinates of next event */
      disk->centr->lonlat[0] = Uni()*width;
      disk->centr->lonlat[1] = Uni()*height;      

      disk->time = curr_t;
      disk->mmod = mmod;

      /* printf("\n. Disk %s has %d lindisk nodes and %d disks under",disk->id,disk->n_ldsk_a,disk->n_disk_under); */

      /* New lindisk (will not be used if no lineage is hit)  */
      new_ldsk = MIGREP_Make_Lindisk_Node(n_dim);
      MIGREP_Init_Lindisk_Node(new_ldsk,disk,n_dim);

      /* Sample the location of new_ldsk uniformly in the disk. */
      /* Takes into account the limits of the landscape */
      MIGREP_Runif_Rectangle_Overlap(new_ldsk,disk,mmod);

      n_hit          = 0;
      n_lineages_new = 0;
      For(i,n_lineages)
        {
          ldsk_a[i]->prev = NULL;

          if(MIGREP_Is_In_Disk(ldsk_a[i]->coord,disk) == YES) 
            {
              if(Uni() < mmod->mu)
                {
                  /* printf("\n. Hit and die %s @ time %f. Center: %f %f Go to %f %f", */
                  /*        ldsk_a[i]->coord->id, */
                  /*        disk->time, */
                  /*        disk->centr->lonlat[0], */
                  /*        disk->centr->lonlat[1], */
                  /*        new_ldsk->coord->lonlat[0], */
                  /*        new_ldsk->coord->lonlat[1]); */
                  
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
          /* PhyML_Printf("\n. Coalescent @ disk %s on ldsk %s",disk->id,disk->ldsk->coord->id); */
        }
 
      ldsk_a = (t_ldsk **)mCalloc(n_otu,sizeof(t_ldsk *));
      For(i,n_lineages) ldsk_a[i] = ldsk_a_tmp[i];
      
      /* For(i,n_lineages) */
      /*   { */
      /*     printf("\n. disk %s [%15f] %3d %s %15f %15f", */
      /*            disk?disk->id:"xxx", */
      /*            disk?disk->time:0.0, */
      /*            i, */
      /*            ldsk_a[i]->coord->id, */
      /*            ldsk_a[i]->coord->lonlat[0], */
      /*            ldsk_a[i]->coord->lonlat[1]); */
      /*   } */

      n_disk++;


      if(n_lineages == 1) break;

      disk->prev = MIGREP_Make_Disk_Event(n_dim,n_otu);
      MIGREP_Init_Disk_Event(disk->prev,n_dim,NULL);
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


  MIGREP_Lk(tree);
  /* Exit("\n"); */
  /* MIGREP_MCMC(tree); */

  return(tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
/* Test whether coord is in disk. Will actually only works in disk */
/* is a rectangle... */

int MIGREP_Is_In_Disk(t_geo_coord *coord, t_dsk *disk)
{
  int i;
  /* PhyML_Printf("\n<> disk %s %d",disk->id,disk->centr->dim); */

  if(!disk->centr->dim) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);

  For(i,disk->centr->dim)
    {
      /* PhyML_Printf("\n<> disk %s ldsk %s %f ctr: %f",disk->id,coord->id,coord->lonlat[i],disk->centr->lonlat[i]); */
      if(FABS(coord->lonlat[i] - disk->centr->lonlat[i]) > disk->mmod->rad + 1.E-20)
        {
          /* PhyML_Printf("  NO"); */
          return(NO);
        }
    }
  /* PhyML_Printf("  YES"); */
  return(YES);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl MIGREP_Lk(t_tree *tree)
{
  phydbl lnL;
  phydbl log_mu,log_one_mu;
  int i,n_inter;
  short int was_hit, n_hit;
  t_migrep_mod *mmod;
  t_dsk *disk;

  mmod = tree->mmod;
  disk = tree->disk;

  if(disk->next)  Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  if(!disk->prev) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);    

  log_mu     = LOG(mmod->mu);
  log_one_mu = LOG(1. - mmod->mu);  
  lnL        = 0.0;

  MIGREP_Update_Lindisk_List(tree);

  /* PhyML_Printf("\n\n. New likelihood call"); */
  do
    {
      /* Likelihood for the disk center */
      For(i,disk->mmod->n_dim) lnL -= LOG(disk->mmod->lim->lonlat[i]);
      
      /* PhyML_Printf("\n. Likelihood - disk %s has %d lindisk nodes [%f] rad: %f",disk->id,disk->n_ldsk_a,lnL,mmod->rad); */
      /* fflush(NULL); */

      was_hit = NO;
      n_hit   = 0;
      
      For(i,disk->n_ldsk_a)
        {
          was_hit = disk->ldsk_a[i]->prev->disk == disk->prev;
          
          if(was_hit == YES) n_hit++;
          
          if(was_hit == YES && !disk->prev->ldsk)
            {
              PhyML_Printf("\n== disk: %s disk->prev: %s",disk->id,disk->prev->id);
              Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
            }

          if(MIGREP_Is_In_Disk(disk->ldsk_a[i]->coord,disk->prev) == YES) /* 'Departure' point is in disk */
            {
              if(was_hit == YES) /* was hit */
                {
                  /* PhyML_Printf("\n. %d/%d lindisk %s (@ %f) was hit and gave %s (@ %f)", */
                  /*              i,disk->n_ldsk_a, */
                  /*              disk->ldsk_a[i]->coord->id, */
                  /*              disk->ldsk_a[i]->prev->coord->id, */
                  /*              disk->ldsk_a[i]->coord->lonlat[0], */
                  /*              disk->ldsk_a[i]->prev->coord->lonlat[0]); */
                  if(MIGREP_Is_In_Disk(disk->prev->ldsk->coord,disk->prev) == YES) /* 'Arrival' point is in disk */
                    {
                      lnL += log_mu;
                    }
                  else /* Landed outside the disk */
                    {
                      /* PhyML_Printf("\n. Landed outside"); */
                      /* PhyML_Printf("\n. ldsk: %s %f %f centr: %f %f rad: %f", */
                      /*              disk->ldsk_a[i]->prev->coord->id, */
                      /*              disk->ldsk_a[i]->prev->coord->lonlat[0], */
                      /*              disk->ldsk_a[i]->prev->coord->lonlat[1], */
                      /*              disk->prev->centr->lonlat[0], */
                      /*              disk->prev->centr->lonlat[1], */
                      /*              tree->mmod->rad); */
                      /* fflush(NULL); */

                      /* mmod->c_lnL = UNLIKELY; */
                      /* return UNLIKELY; */
                    }
                }
              else /* was not hit */
                {
                  lnL += log_one_mu;
                  /* PhyML_Printf("\n. %d/%d lindisk %s was not hit, lnL: %f [%f %f]",i,disk->n_ldsk_a,disk->ldsk_a[i]->coord->id,lnL,log_one_mu,mmod->mu); */
                }
            }
          else
            {
              /* PhyML_Printf("\n. %s was hit: %d %s %s %s lnL:%f",disk->ldsk_a[i]->coord->id,was_hit,disk->ldsk_a[i]->prev->disk->id,disk->prev->id,disk->id,lnL); */
              /* ldsk has changed spatial coordinates but was outside disk, so it could not have been hit  */
              if(was_hit == YES)
                {
                  /* printf("\n. FAIL lindisk: %s [%.3f %.3f] prev: %s [%.3f %.3f] centr: [%.3f %.3f] rad: %.3f", */
                  /*        disk->ldsk_a[i]->coord->id, */
                  /*        disk->ldsk_a[i]->coord->lonlat[0], */
                  /*        disk->ldsk_a[i]->coord->lonlat[1], */
                  /*        disk->ldsk_a[i]->prev->coord->id, */
                  /*        disk->ldsk_a[i]->prev->coord->lonlat[0], */
                  /*        disk->ldsk_a[i]->prev->coord->lonlat[1], */
                  /*        disk->prev->centr->lonlat[0], */
                  /*        disk->prev->centr->lonlat[1], */
                  /*        disk->mmod->rad); */
                  /* fflush(NULL); */

                  mmod->c_lnL = UNLIKELY;
                  return UNLIKELY;
                }
            }
        }
      
      /* a hit occurred */
      if(n_hit >= 1) 
        {
          lnL += MIGREP_Log_Dunif_Rectangle_Overlap(disk->prev->ldsk,disk->prev,mmod);

          /* For(i,disk->prev->ldsk->n_next)  */
          /*   lnL += MIGREP_Log_Dunif_Rectangle_Overlap(disk->prev->ldsk->next[i],disk->prev,mmod); */
        }

      disk = disk->prev;

      if(!disk->prev) break;
    }
  while(1);
  
  n_inter = MIGREP_Total_Number_Of_Intervals(tree);

  lnL += (n_inter)*LOG(mmod->lbda) + mmod->lbda*disk->time;
  /* lnL -= LnGamma((phydbl)(n_inter)); */
  
  /* lnL += Dpois((phydbl)n_inter,-mmod->lbda*disk->time,YES); */

  if(isinf(lnL) || isnan(lnL)) lnL = UNLIKELY;

  mmod->c_lnL = lnL;

  return(lnL);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIGREP_MCMC(t_tree *tree)
{
  t_mcmc *mcmc;
  int move,i;
  phydbl u;
  t_dsk *disk;

  mcmc = MCMC_Make_MCMC_Struct();
  MCMC_Complete_MCMC(mcmc,tree);

  tree->mcmc = mcmc;

  mcmc->io               = NULL;
  mcmc->is               = NO;
  mcmc->use_data         = YES;
  mcmc->run              = 0;
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
  mcmc->sample_interval  = 1E+3;
  
  MIGREP_Lk(tree);
  MIGREP_LnPrior_Radius(tree);
  MIGREP_LnPrior_Mu(tree);
  MIGREP_LnPrior_Lbda(tree);

  printf("\n# before rand lnL: %f",tree->mmod->c_lnL);
  printf("\n# ninter: %d",MIGREP_Total_Number_Of_Intervals(tree));
  printf("\n# ncoal: %d",MIGREP_Total_Number_Of_Coal_Disks(tree));
  printf("\n# nhits: %d",MIGREP_Total_Number_Of_Hit_Disks(tree));

  /* gtk_widget_queue_draw(tree->draw_area); */
  /* sleep(3); */

  disk = tree->disk;
  while(disk->prev) disk = disk->prev;

  tree->mmod->lbda = Uni()*(tree->mmod->max_lbda - tree->mmod->min_lbda) + tree->mmod->min_lbda;
  tree->mmod->mu   = Uni()*(tree->mmod->max_mu - tree->mmod->min_mu) + tree->mmod->min_mu;
  /* tree->mmod->rad  = Uni()*(tree->mmod->max_rad - tree->mmod->min_rad) + tree->mmod->min_rad; */
  tree->mmod->rad  = tree->mmod->max_rad;

  mcmc->always_yes = YES;
  mcmc->run        = 0;
  do
    {
        MCMC_MIGREP_Delete_Disk(tree);
        MCMC_MIGREP_Insert_Disk(tree);
        MCMC_MIGREP_Move_Disk_Centre(tree);
        MCMC_MIGREP_Move_Disk_Updown(tree);
        MCMC_MIGREP_Swap_Disk(tree);
        MCMC_MIGREP_Delete_Hit(tree);
        MCMC_MIGREP_Insert_Hit(tree);
        MCMC_MIGREP_Move_Ldsk(tree);
        MCMC_MIGREP_Shift_Ldsk_To_Centre(tree);
        MCMC_MIGREP_Shift_Centre_To_Median(tree);
        mcmc->run++;
        
        if(!(mcmc->run%mcmc->sample_interval))
          {
            PhyML_Printf("\n %13d %13f %13f %13f %13f %13d %13d %13d %13f %13f",
                         mcmc->run,
                         tree->mmod->c_lnL,
                         tree->mmod->lbda,
                         tree->mmod->mu,
                         tree->mmod->rad,
                         MIGREP_Total_Number_Of_Intervals(tree),
                         MIGREP_Total_Number_Of_Coal_Disks(tree),
                         MIGREP_Total_Number_Of_Hit_Disks(tree),
                         disk->ldsk->coord->lonlat[0],
                         disk->ldsk->coord->lonlat[1]);
          }
    }
  while(mcmc->run < 10.*mcmc->sample_interval);

  MIGREP_Lk(tree);
  printf("\n# after rand lnL: %f",tree->mmod->c_lnL);
  printf("\n# ninter: %d",MIGREP_Total_Number_Of_Intervals(tree));
  printf("\n# ncoal: %d",MIGREP_Total_Number_Of_Coal_Disks(tree));
  printf("\n# nhits: %d",MIGREP_Total_Number_Of_Hit_Disks(tree));
  printf("\n# root pos: %f %f",disk->ldsk->coord->lonlat[0],disk->ldsk->coord->lonlat[1]);

  PhyML_Printf("\n %13s %13s %13s %13s %13s %13s %13s %13s %13s %13s",
               "run",
               "lnL",
               "lbda",
               "mu",
               "rad",
               "nInt",
               "nCoal",
               "nHit",
               "xRootLon",
               "xRootLat");

  /* gtk_widget_queue_draw(tree->draw_area); */
  /* sleep(2); */

  PhyML_Printf("\n *** \n");

  mcmc->always_yes = NO;
  do
    {
      if(mcmc->run == (int)mcmc->chain_len * 0.01) 
        {
          For(i,mcmc->n_moves) tree->mcmc->adjust_tuning[i] = NO;
        }

      u = Uni();

      For(move,tree->mcmc->n_moves) if(tree->mcmc->move_weight[move] > u) break;
      
      if(!strcmp(tree->mcmc->move_name[move],"migrep_lbda"))
        MCMC_MIGREP_Lbda(tree);

      if(!strcmp(tree->mcmc->move_name[move],"migrep_mu"))
        MCMC_MIGREP_Mu(tree);

      if(!strcmp(tree->mcmc->move_name[move],"migrep_rad"))
        MCMC_MIGREP_Radius(tree);

      if(!strcmp(tree->mcmc->move_name[move],"migrep_delete_disk"))
        MCMC_MIGREP_Delete_Disk(tree);

      if(!strcmp(tree->mcmc->move_name[move],"migrep_insert_disk"))
        MCMC_MIGREP_Insert_Disk(tree);

      if(!strcmp(tree->mcmc->move_name[move],"migrep_move_disk_ct"))
        MCMC_MIGREP_Move_Disk_Centre(tree);

      if(!strcmp(tree->mcmc->move_name[move],"migrep_move_disk_ud"))
        MCMC_MIGREP_Move_Disk_Updown(tree);

      if(!strcmp(tree->mcmc->move_name[move],"migrep_swap_disk"))
        MCMC_MIGREP_Swap_Disk(tree);

      if(!strcmp(tree->mcmc->move_name[move],"migrep_delete_hit"))
        MCMC_MIGREP_Delete_Hit(tree);

      if(!strcmp(tree->mcmc->move_name[move],"migrep_insert_hit"))
        MCMC_MIGREP_Insert_Hit(tree);

      if(!strcmp(tree->mcmc->move_name[move],"migrep_move_ldsk"))
        MCMC_MIGREP_Move_Ldsk(tree);

      if(!strcmp(tree->mcmc->move_name[move],"migrep_shift_ldsk_to_centre"))
        MCMC_MIGREP_Shift_Ldsk_To_Centre(tree);

      if(!strcmp(tree->mcmc->move_name[move],"migrep_shift_centr_to_median"))
        MCMC_MIGREP_Shift_Centre_To_Median(tree);

      /* if(!strcmp(tree->mcmc->move_name[move],"migrep_shift_ldsk_path")) */
      /*   MCMC_MIGREP_Shift_Ldsk_Path(tree); */

      /* if(!strcmp(tree->mcmc->move_name[move],"migrep_shift_disk_path")) */
      /*   MCMC_MIGREP_Shift_Disk_Path(tree); */

      /* if(!strcmp(tree->mcmc->move_name[move],"migrep_subtree")) */
      /*   MCMC_MIGREP_New_Subtree(tree); */
      
      /* if(!strcmp(tree->mcmc->move_name[move],"migrep_rad_mu")) */
      /*   MCMC_MIGREP_Rad_Mu(tree); */

      tree->mcmc->run++;
      MCMC_Get_Acc_Rates(tree->mcmc);
      
      if(!(tree->mcmc->run%tree->mcmc->sample_interval))
        {
          PhyML_Printf("\n %13d %13f %13f %13f %13f %13d %13d %13d %13f %13f",
                       tree->mcmc->run,
                       tree->mmod->c_lnL,
                       tree->mmod->lbda,
                       tree->mmod->mu,
                       tree->mmod->rad,
                       MIGREP_Total_Number_Of_Intervals(tree),
                       MIGREP_Total_Number_Of_Coal_Disks(tree),
                       MIGREP_Total_Number_Of_Hit_Disks(tree),
                       disk->ldsk->coord->lonlat[0],
                       disk->ldsk->coord->lonlat[1]);
          /* gtk_widget_queue_draw(tree->draw_area); */
          /* sleep(1); */
        }
    }
  while(tree->mcmc->run < tree->mcmc->chain_len);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl MIGREP_Wrap_Lk(t_edge *b, t_tree *tree, supert_tree *stree)
{
  return MIGREP_Lk(tree);
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
  
  /* MIGREP_Lk(tree); */
  /* printf("\n. >> Lk: %f rad: %f",tree->mmod->c_lnL,tree->mmod->rad); */
  /* sleep(5); */
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
    between y_ldsk and o_ldsk: we need to generate some first. n_cur_disk
    is the current number of disks between y_ldsk and o_ldsk.
*/
int MIGREP_One_New_Traj(t_ldsk *y_ldsk, t_ldsk *o_ldsk, int dir_o_y, t_dsk *xtra_dsk, int n_cur_disk, t_tree *tree)
{
  t_migrep_mod *mmod;
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
  For(i,mmod->n_dim)
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
      For(i,n_new_disk-1)  disk_new[i] = MIGREP_Make_Disk_Event(mmod->n_dim,tree->n_otu);
      if(xtra_dsk != NULL) disk_new[n_new_disk-1] = xtra_dsk;
      else                 disk_new[n_new_disk-1] = MIGREP_Make_Disk_Event(mmod->n_dim,tree->n_otu);

      For(i,n_new_disk)  MIGREP_Init_Disk_Event(disk_new[i],mmod->n_dim,mmod);

      /* Times of these new disks. If xtra_dsk != NULL, then make sure you do not */
      /* reset the time of that disk  */
      n = (xtra_dsk != NULL) ? (n_new_disk-1) : (n_new_disk);
      For(i,n)
        disk_new[i]->time =
        Uni()*(y_ldsk->disk->time - o_ldsk->disk->time) + o_ldsk->disk->time;       
     
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
          /* printf("\n# Add ldsk %s to %s",disk_new[i]->ldsk->coord->id,disk_new[i]->id); fflush(NULL); */
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
  MIGREP_One_New_Traj_Given_Disk(y_ldsk,o_ldsk,tree);  

  return(n_new_disk);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

/* Generate a new trajectory, including disk event centers, between 
  'y_ldsk' a ``young'' lindisk event and 'o_ldsk' an old one. 'y_ldsk 
  and 'o_ldsk' remain unaffected. Disk events between these two ldsk 
  should already be set. 
*/
void MIGREP_One_New_Traj_Given_Disk(t_ldsk *y_ldsk, t_ldsk *o_ldsk, t_tree *tree)
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
      
      MIGREP_Copy_Geo_Coord(ldsk->coord);
      MIGREP_Copy_Geo_Coord(ldsk->disk->centr);

      For(i,tree->mmod->n_dim)
        {
          min = 
            MAX(0.0,
                MAX(ldsk->coord->lonlat[i] - 2.*rad,
                    o_ldsk->disk->centr->lonlat[i] - rad*(2.*(n_disk_btw-k)-1.)));

          max = 
            MIN(ldsk->disk->mmod->lim->lonlat[i],
                MIN(ldsk->coord->lonlat[i] + 2.*rad,
                    o_ldsk->disk->centr->lonlat[i] + rad*(2.*(n_disk_btw-k)-1.)));
          
          if(max < min) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
                
          /* New coordinates for the lindisk */
          ldsk->coord->lonlat[i] = Uni()*(max - min) + min;
        }

      /* New coordinate for the centre of the corresponding disk event */
      MIGREP_Get_Min_Max_Disk_Given_Ldsk(ldsk->disk,&min_disk_coord,&max_disk_coord,tree);
      For(i,tree->mmod->n_dim) ldsk->disk->centr->lonlat[i] = Uni()*(max_disk_coord[i] - min_disk_coord[i]) + min_disk_coord[i];
      Free(min_disk_coord);
      Free(max_disk_coord);

      ldsk = ldsk->prev;
      k++;
    }
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl MIGREP_Uniform_Path_Density(t_ldsk *y_ldsk, t_ldsk *o_ldsk, t_tree *tree)
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
      
      For(i,ldsk->disk->mmod->n_dim)
        {
          min = 
            MAX(0.0,
                MAX(ldsk->coord->lonlat[i] - 2.*rad,
                    o_ldsk->disk->centr->lonlat[i] - rad*(2.*(n_disk_btw-k)-1.)));

          max = 
            MIN(ldsk->disk->mmod->lim->lonlat[i],
                MIN(ldsk->coord->lonlat[i] + 2.*rad,
                    o_ldsk->disk->centr->lonlat[i] + rad*(2.*(n_disk_btw-k)-1.)));

          if(max < min) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);

          log_dens -= LOG(max - min);          
        }

      MIGREP_Get_Min_Max_Disk_Given_Ldsk(ldsk->disk,&min_disk_coord,&max_disk_coord,tree);
      For(i,tree->mmod->n_dim) log_dens -= LOG(max_disk_coord[i] - min_disk_coord[i]);
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

/* void MIGREP_Update_Lindisk_List(phydbl time, t_ldsk **list, int *pos, t_dsk *disk) */
/* { */
/*   t_dsk *root_dsk; */

/*   *pos = 0; */
/*   root_dsk = disk; */
/*   while(root_dsk->prev) root_dsk = root_dsk->prev; */
/*   /\* printf("\n. root_dsk: %s",root_dsk?root_dsk->id:"xx"); *\/ */
/*   MIGREP_Update_Lindisk_List_Pre(root_dsk->ldsk,time,list,pos); */
/* } */

/* /\*\//////////////////////////////////////////////////////////// */
/* ////////////////////////////////////////////////////////////\*\/ */

/* void MIGREP_Update_Lindisk_List_Pre(t_ldsk *ldsk, phydbl time, t_ldsk **list, int *pos) */
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
/*       For(i,ldsk->n_next) */
/*         MIGREP_Update_Lindisk_List_Pre(ldsk->next[i],time,list,pos);     */
/*     } */
/*   else Generic_Exit(__FILE__,__LINE__,__FUNCTION__); */
    
/* } */

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void MIGREP_Update_Lindisk_List(t_tree *tree)
{
  MIGREP_Update_Lindisk_List_Pre(tree->disk->prev);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void MIGREP_Update_Lindisk_List_Pre(t_dsk *disk)
{
  if(!disk) return;
  else
    {
      int i;
      disk->n_ldsk_a = 0;
      For(i,disk->next->n_ldsk_a)
        {
          if(disk->next->ldsk_a[i]->prev != disk->ldsk)
            {
              disk->ldsk_a[disk->n_ldsk_a] = disk->next->ldsk_a[i];
              disk->n_ldsk_a++;
            }
        }
      
      if(disk->ldsk) 
        {
          disk->ldsk_a[disk->n_ldsk_a] = disk->ldsk;
          disk->n_ldsk_a++;          
        }

      if(disk->n_ldsk_a == 0) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
      MIGREP_Update_Lindisk_List_Pre(disk->prev);    
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
  int i,j;
  t_ldsk *ldisk;

  disk = tree->disk;
  while(disk->prev) disk = disk->prev;
  do
    {
      printf("\n%c Disk: %s @ %7.3f has %3s on it is_coal? %2d rad: %f coord ",
             sign,
             disk->id,
             disk->time,disk->ldsk?disk->ldsk->coord->id:NULL,
             disk->ldsk?disk->ldsk->is_coal:-1,
             tree->mmod->rad); fflush(NULL);
      For(j,tree->mmod->n_dim) PhyML_Printf(" %f",disk->centr->lonlat[j]);

      MIGREP_Update_Lindisk_List(tree);

      For(i,disk->n_ldsk_a)
        {
          ldisk = disk->ldsk_a[i];

          printf("\n%c ldisk: %s prev: %s",
                 sign,
                 ldisk->coord->id,
                 ldisk->prev ? ldisk->prev->coord->id : NULL);

          For(j,tree->mmod->n_dim) 
            {
              PhyML_Printf(" %f",ldisk->coord->lonlat[j]);

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

      disk = disk->next;      
    }
  while(disk);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void MIGREP_Check_Struct(t_tree *tree)
{
  t_dsk *disk;
  int i,j;
  t_ldsk *ldisk;

  disk = tree->disk;
  while(disk->prev) disk = disk->prev;
  do
    {
      /* MIGREP_Update_Lindisk_List(disk->time,disk->ldsk_a,&(disk->n_ldsk_a),disk); */
      MIGREP_Update_Lindisk_List(tree);

      For(i,disk->n_ldsk_a)
        {
          ldisk = disk->ldsk_a[i];
          if(ldisk->prev != NULL)
            {
              For(j,tree->mmod->n_dim)
                {
                  if(FABS(ldisk->coord->lonlat[j] - 
                          ldisk->prev->coord->lonlat[j]) > 2.*tree->mmod->rad)
                    {
                      MIGREP_Print_Struct('=',tree);
                      PhyML_Printf("\n== %f %f %f",
                                   ldisk->coord->lonlat[j], 
                                   ldisk->prev->coord->lonlat[j],
                                   2.*tree->mmod->rad);
                      PhyML_Printf("\n== Radius: %f",tree->mmod->rad);
                      PhyML_Printf("\n== Check ldsk %s",ldisk->coord->id);
                      PhyML_Printf("\n== Centr: %f",ldisk->prev->disk->centr->lonlat[j]);
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

void MIGREP_Copy_Geo_Coord(t_geo_coord *t)
{
  int i;

  For(i,t->dim) t->cpy->lonlat[i] = t->lonlat[i];  
  t->cpy->dim = t->dim;
  strcpy(t->cpy->id,t->id);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void MIGREP_Restore_Geo_Coord(t_geo_coord *t)
{
  int i;

  For(i,t->dim) t->lonlat[i] = t->cpy->lonlat[i];  
  t->dim = t->cpy->dim;
  strcpy(t->id,t->cpy->id);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

int MIGREP_Total_Number_Of_Intervals(t_tree *tree)
{
  t_dsk *disk;
  int n_intervals;

  if(tree->disk->next) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);

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

int MIGREP_Total_Number_Of_Hit_Disks(t_tree *tree)
{
  t_dsk *disk;
  int n_hit_disks;

  if(tree->disk->next) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);

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

int MIGREP_Total_Number_Of_Coal_Disks(t_tree *tree)
{
  t_dsk *disk;
  int n_coal_disks;

  if(tree->disk->next) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  disk = tree->disk;
  n_coal_disks = 0;
  while(disk)
    {
      if(disk->ldsk && disk->ldsk->is_coal == YES) n_coal_disks++;
      disk = disk->prev;
    }

  return(n_coal_disks);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl MIGREP_Log_Dunif_Rectangle_Overlap(t_ldsk *ldsk, t_dsk *disk, t_migrep_mod *mmod)
{
  phydbl up, down, left, rght;
  phydbl log_dens,l;
  phydbl main_mass;

  log_dens  = 0.0;
  main_mass = 1. - mmod->soft_bound_area;

  up   = MIN(disk->centr->lonlat[0] + mmod->rad, mmod->lim->lonlat[0]);
  down = MAX(disk->centr->lonlat[0] - mmod->rad, 0.0);
  rght = MIN(disk->centr->lonlat[1] + mmod->rad, mmod->lim->lonlat[1]);
  left = MAX(disk->centr->lonlat[1] - mmod->rad, 0.0);


  /* l = main_mass / (.5*(1.-main_mass)*(up - down)); */
  /* if(ldsk->coord->lonlat[0] < down)    log_dens += LOG(.5*(1.-main_mass)) + LOG(l) - l*(down - ldsk->coord->lonlat[0]); */
  /* else if(ldsk->coord->lonlat[0] > up) log_dens += LOG(.5*(1.-main_mass)) + LOG(l) - l*(ldsk->coord->lonlat[0] - up); */
  /* else log_dens += LOG(main_mass) - LOG(up - down); */

  /* l = main_mass / (.5*(1.-main_mass)*(rght - left)); */
  /* if(ldsk->coord->lonlat[1] < left)      log_dens += LOG(.5*(1.-main_mass)) + LOG(l) - l*(left - ldsk->coord->lonlat[1]); */
  /* else if(ldsk->coord->lonlat[1] > rght) log_dens += LOG(.5*(1.-main_mass)) + LOG(l) - l*(ldsk->coord->lonlat[1] - rght); */
  /* else log_dens += LOG(main_mass) - LOG(rght - left); */


  if(ldsk->coord->lonlat[0] < down)  return UNLIKELY;
  if(ldsk->coord->lonlat[0] > up)    return UNLIKELY;
  if(ldsk->coord->lonlat[1] < left)  return UNLIKELY;
  if(ldsk->coord->lonlat[1] > rght)  return UNLIKELY;

  log_dens = -LOG(up-down)-LOG(rght-left);

  return(log_dens);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
/* Samples uniformly within a rectangle (with truncation for border) */
/* and returns the corresponding log density */
phydbl MIGREP_Runif_Rectangle_Overlap(t_ldsk *ldsk, t_dsk *disk, t_migrep_mod *mmod)
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
  return(LOG(up-down)+LOG(rght-left));
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl MIGREP_Wrap_Prior_Radius(t_edge *e, t_tree *tree, supert_tree *st)
{
  return MIGREP_LnPrior_Radius(tree);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl MIGREP_LnPrior_Lbda(t_tree *tree)
{
  tree->mmod->c_ln_prior_lbda = LOG(tree->mmod->prior_param_lbda) - tree->mmod->prior_param_lbda * tree->mmod->lbda;
  return(tree->mmod->c_ln_prior_lbda);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl MIGREP_LnPrior_Mu(t_tree *tree)
{
  tree->mmod->c_ln_prior_mu = LOG(tree->mmod->prior_param_mu) - tree->mmod->prior_param_mu * (1. - tree->mmod->mu) - LOG(1.-EXP(-tree->mmod->prior_param_mu));
  return(tree->mmod->c_ln_prior_mu);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl MIGREP_LnPrior_Radius(t_tree *tree)
{
  tree->mmod->c_ln_prior_rad = LOG(tree->mmod->prior_param_rad) - tree->mmod->prior_param_rad * tree->mmod->rad;
  return(tree->mmod->c_ln_prior_rad);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void MIGREP_Initial_Ldsk_Pos(t_tree *tree)
{
  t_dsk *disk;
  int i,j;
  phydbl mean;

  disk = tree->disk->prev;

  do
    {
      if(disk->ldsk)
        {
          For(i,tree->mmod->n_dim)
            {
              mean = 0.0;
              For(j,disk->ldsk->n_next) mean += disk->ldsk->next[j]->coord->lonlat[i];
              disk->ldsk->coord->lonlat[i] = mean / (phydbl)disk->ldsk->n_next;
            }
        }
      disk = disk->prev;
    }
  while(disk);


}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl MIGREP_Min_Radius(t_tree *tree)
{
  phydbl ori_rad, min_rad;

  ori_rad = tree->mmod->rad;
  tree->mmod->rad = tree->mmod->max_rad;
  do
    {
      MIGREP_Lk(tree);
      tree->mmod->rad -= 1.0;
    }
  while(tree->mmod->c_lnL > UNLIKELY + 0.1);

  min_rad = tree->mmod->rad + 2.0;
  tree->mmod->rad = ori_rad;
  MIGREP_Lk(tree);
  return(min_rad);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
/* Get the minimum and maximum values a ldsk can take, given the position of the disk centre */
void MIGREP_Get_Min_Max_Ldsk_Given_Disk(t_ldsk *ldsk, phydbl **min, phydbl **max, t_tree *tree)
{
  phydbl *loc_min,*loc_max;
  int i;

  if(!ldsk->disk->next) return;

  loc_min = (phydbl *)mCalloc(tree->mmod->n_dim, sizeof(phydbl));
  loc_max = (phydbl *)mCalloc(tree->mmod->n_dim, sizeof(phydbl));

  For(i,tree->mmod->n_dim)
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
void MIGREP_Get_Min_Max_Disk_Given_Ldsk(t_dsk *disk, phydbl **min, phydbl **max, t_tree *tree)
{
  phydbl *loc_min,*loc_max;
  int i,j;
  phydbl tmp_min, tmp_max;

  if(!disk->next) return;

  loc_min = (phydbl *)mCalloc(tree->mmod->n_dim,sizeof(phydbl));
  loc_max = (phydbl *)mCalloc(tree->mmod->n_dim,sizeof(phydbl));

  if(!disk->ldsk)
    {
      For(i,tree->mmod->n_dim)
        {
          loc_min[i] = 0.0;
          loc_max[i] = tree->mmod->lim->lonlat[i];
        }
    }
  else
    {
      For(i,tree->mmod->n_dim)
        {
          tmp_min = +INFINITY;
          tmp_max = -INFINITY;
          For(j,disk->ldsk->n_next)
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
                    
          if(loc_max[i] < loc_min[i]) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
        }
    }

  (*min) = loc_min;
  (*max) = loc_max;
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void MIGREP_Update_Disk_Ldsk_Subtree(t_ldsk *root_ldsk, t_tree *tree)
{
  int i;  
  For(i,root_ldsk->n_next) MIGREP_Update_Disk_Ldsk_Subtree_Pre(root_ldsk,root_ldsk->next[i],root_ldsk,tree);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void MIGREP_Update_Disk_Ldsk_Subtree_Pre(t_ldsk *old_ldsk, t_ldsk *young_ldsk, t_ldsk *root_ldsk, t_tree *tree)
{  
  if(!young_ldsk->disk->next) 
    {
      MIGREP_One_New_Traj_Given_Disk(young_ldsk,root_ldsk,tree);
      return;
    }
  else
    {
      int i;
      MIGREP_Update_Disk_Ldsk_Subtree_Pre(young_ldsk,young_ldsk->next[0],root_ldsk,tree);
      if(young_ldsk->n_next > 1) 
        {
          for(i=1;i<young_ldsk->n_next;i++) MIGREP_Update_Disk_Ldsk_Subtree_Pre(young_ldsk,young_ldsk->next[i],young_ldsk,tree);
        }
    }
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void MIGREP_Restore_Disk_Ldsk_Subtree(t_ldsk *root_ldsk, t_tree *tree)
{
  int i;  
  For(i,root_ldsk->n_next) MIGREP_Restore_Disk_Ldsk_Subtree_Pre(root_ldsk,root_ldsk->next[i],tree);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void MIGREP_Restore_Disk_Ldsk_Subtree_Pre(t_ldsk *old_ldsk, t_ldsk *young_ldsk, t_tree *tree)
{  

  if(!young_ldsk->disk->next) return;
  else
    {
      int i;

      MIGREP_Restore_Geo_Coord(young_ldsk->coord);
      MIGREP_Restore_Geo_Coord(young_ldsk->disk->centr);

      For(i,young_ldsk->n_next)
        {
          MIGREP_Restore_Disk_Ldsk_Subtree_Pre(young_ldsk,young_ldsk->next[i],tree);
        }
    }
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void MIGREP_Proposal_Disk_Ldsk_Subtree(t_ldsk *root_ldsk, phydbl *logdens, t_tree *tree)
{
  int i;  
  (*logdens) = 0.0;
  For(i,root_ldsk->n_next) MIGREP_Proposal_Disk_Ldsk_Subtree_Pre(root_ldsk,root_ldsk->next[i],root_ldsk,logdens,tree);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void MIGREP_Proposal_Disk_Ldsk_Subtree_Pre(t_ldsk *old_ldsk, t_ldsk *young_ldsk, t_ldsk *root_ldsk, phydbl *logdens, t_tree *tree)
{  

  if(!young_ldsk->disk->next) 
    {
      (*logdens) += MIGREP_Uniform_Path_Density(young_ldsk,root_ldsk,tree);  
      return;
    }
  else
    {
      int i;
      For(i,young_ldsk->n_next)
        {
          MIGREP_Proposal_Disk_Ldsk_Subtree_Pre(young_ldsk,young_ldsk->next[i],root_ldsk,logdens,tree);
        }
    }
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
/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
