/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include "slfv.h"


/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl SLFV_Prob_Two_Lineages_Coal(t_ldsk *l0, t_ldsk *l1, t_tree *tree)
{
  phydbl eucl,prob,area,W,H;

  assert(tree->mmod->n_dim == 2);

  W = (tree->mmod->lim_up->lonlat[0]-tree->mmod->lim_do->lonlat[0]);
  H = (tree->mmod->lim_up->lonlat[1]-tree->mmod->lim_do->lonlat[1]);
  
  area = W*H;  

  eucl = Euclidean_Dist(l0->coord,l1->coord);
  
  prob = 0.0;
  prob += 2.*log(tree->mmod->mu);
  prob += 2.*log(tree->mmod->rad);
  prob += log(PI);
  prob -= log(4.*area);
  prob -= 0.25*eucl*eucl/(tree->mmod->rad*tree->mmod->rad);


  prob +=
    log(erf((l0->coord->lonlat[0] + l1->coord->lonlat[0])/(2.*tree->mmod->rad)) +
        erf((2.*W - l0->coord->lonlat[0] - l1->coord->lonlat[0])/(2.*tree->mmod->rad)));

  prob +=
    log(erf((l0->coord->lonlat[1] + l1->coord->lonlat[1])/(2.*tree->mmod->rad)) +
        erf((2.*H - l0->coord->lonlat[1] - l1->coord->lonlat[1])/(2.*tree->mmod->rad)));

  prob = exp(prob);
  
  return(prob);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl SLFV_Prob_Two_Random_Lineages_Coal_One_Event(phydbl w, phydbl h, phydbl mu, phydbl rad)
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
// Coalescence rate with time expressed in calendar unit
phydbl SLFV_Coalescence_Rate(t_tree *tree)
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
// Log-density of path, conditionned on its length
phydbl SLFV_Path_Logdensity(t_ldsk *young, t_ldsk *old, phydbl *sd, t_tree *tree)
{
  int i,j,err;
  t_ldsk *ldsk;
  phydbl lnP,mode,var;
  phydbl X,Y,Xp,Yp;
  /* phydbl slope,inter,mode; */
  int dir_to_young;
  phydbl dt_young,dt_old,sum;
  phydbl pos;
  
  lnP = 0.0;
  mode   = 0.0;
    
  for(i=0;i<tree->mmod->n_dim;i++)
    {
      j    = 0;
      ldsk = young->prev;

      X = old->coord->lonlat[i];
      Y = old->disk->time;

      // Density up to the last jump (to old ldsk) which should
      // not be accounted for (it is not part of the simulation
      // when randomly generating a path using SLFV_Generate_Path
      while(ldsk != old)
        {

          if(ldsk->n_next > 1)
            dir_to_young = PHYREX_Get_Next_Direction(young,ldsk);
          else
            dir_to_young = 0;
          
          Xp = ldsk->next[dir_to_young]->coord->lonlat[i];
          Yp = ldsk->next[dir_to_young]->disk->time;
          
          /* slope = (Yp-Y)/(Xp-X); */
          /* inter = Y - slope*X; */
          
          /* mode = (ldsk->disk->time - inter)/slope; */

          dt_young = fabs(Yp - ldsk->disk->time);
          dt_old = fabs(Y - ldsk->disk->time);
          sum = dt_young + dt_old;
          
          mode = (dt_young/sum)*X + (dt_old/sum)*Xp;
          var = dt_young*dt_old/sum * sd[i] * sd[i]; 
          
          assert(ldsk != NULL);
                    
          lnP += Log_Dnorm_Trunc(ldsk->coord->lonlat[i],
                                    mode,
                                    sqrt(var),
                                    tree->mmod->lim_do->lonlat[i],
                                    tree->mmod->lim_up->lonlat[i],&err);

          /* lnP += Log_Dnorm(ldsk->coord->lonlat[i], */
          /*                     mode, */
          /*                     sqrt(var), */
          /*                     &err); */

          pos = 0.0;
          if(ldsk->coord->lonlat[i] > tree->mmod->lim_up->lonlat[i]) pos = tree->mmod->lim_up->lonlat[i];
          else if(ldsk->coord->lonlat[i] < tree->mmod->lim_do->lonlat[i]) pos = tree->mmod->lim_do->lonlat[i];
          else pos = ldsk->coord->lonlat[i];
          
          lnP += Log_Dnorm_Trunc(ldsk->disk->centr->lonlat[i],
                                 pos,
                                 sd[i],
                                 tree->mmod->lim_do->lonlat[i],
                                 tree->mmod->lim_up->lonlat[i],&err);

          ldsk = ldsk->prev;
          j++;
        }
    }
  
  return(lnP);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void SLFV_Sample_Path(t_ldsk *young, t_ldsk *old, phydbl *sd, phydbl *global_hr, t_tree *tree)
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

          cur_glnL = SLFV_Lk_Range(ldsk->disk,ldsk->prev->disk,tree);
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
                                         sd[i],
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
                                    sd[i],
                                    tree->mmod->lim_do->lonlat[i],
                                    tree->mmod->lim_up->lonlat[i],&err);
              
              /* hr -= Log_Dnorm_Trunc(new_centr_pos, */
              /*                       cur_centr_pos, */
              /*                       0.5*sd, */
              /*                       0.0, */
              /*                       tree->mmod->lim->lonlat[i],&err); */
              
              
              hr += Log_Dnorm_Trunc(cur_ldsk_pos,
                                    cur_centr_pos,
                                    sd[i],
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
          
          new_glnL = SLFV_Lk_Range(ldsk->disk,ldsk->prev->disk,tree);
          if(PHYREX_Total_Number_Of_Intervals(tree) > tree->mmod->max_num_of_intervals) new_glnL = UNLIKELY;
          
          
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

t_ldsk *SLFV_Generate_Path(t_ldsk *young, t_ldsk *old, int n_evt, phydbl *sd, t_tree *tree)
{
  int i,j,swap,err;
  phydbl *time,dum,mode,var;
  t_ldsk *path,**ldsk_a;
  t_dsk *disk;
  phydbl X,Y,Xp,Yp;
  /* phydbl slope,inter; */
  phydbl dt_young,dt_old,sum;
  phydbl pos;
  
  path = NULL;
  
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
          
          /* slope = (Yp-Y)/(Xp-X); */
          /* inter = Y - slope*X; */

          /* if(j == 0) */
          /*   { */
          /*     mode = (young->disk->time - inter)/slope; */
          /*   } */
          /* else */
          /*   { */
          /*     mode = (ldsk_a[j]->next[0]->disk->time - inter)/slope; */
          /*   } */

          dt_young = fabs(Yp - ldsk_a[j]->disk->time);
          dt_old = fabs(Y - ldsk_a[j]->disk->time);
          sum = dt_young + dt_old;
          
          mode = (dt_young/sum)*X + (dt_old/sum)*Xp;
          var = dt_young*dt_old/sum * sd[i] * sd[i]; 
          
          
          /* ldsk_a[j]->coord->lonlat[i] = Rnorm(mode,sqrt(var)); */
          ldsk_a[j]->coord->lonlat[i]  = Rnorm_Trunc(mode,
                                                     sqrt(var),
                                                     tree->mmod->lim_do->lonlat[i],
                                                     tree->mmod->lim_up->lonlat[i],&err);
          
          pos = 0.0;
          if(ldsk_a[j]->coord->lonlat[i] > tree->mmod->lim_up->lonlat[i]) pos = tree->mmod->lim_up->lonlat[i];
          else if(ldsk_a[j]->coord->lonlat[i] < tree->mmod->lim_do->lonlat[i]) pos = tree->mmod->lim_do->lonlat[i];
          else pos = ldsk_a[j]->coord->lonlat[i];

          ldsk_a[j]->disk->centr->lonlat[i] = Rnorm_Trunc(pos,
                                                          sd[i],
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

phydbl SLFV_Effective_Density(t_tree *tree)
{
  return(SLFV_Neighborhood_Size(tree)/(4.*PI*SLFV_Update_Sigsq(tree)));
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl SLFV_Rate_Per_Unit_Area(t_tree *tree)
{
  int i;
  phydbl denom;

  denom = (tree->mmod->lim_up->lonlat[0]-tree->mmod->lim_do->lonlat[0]);  
  for(i=1;i<tree->mmod->n_dim;i++) denom *= (tree->mmod->lim_up->lonlat[i]-tree->mmod->lim_do->lonlat[i]);
  
  return(tree->mmod->lbda / denom);

}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl SLFV_Sample_Rad_From_Prior(t_tree *tree)
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

phydbl SLFV_Update_Radius(t_tree *tree)
{
  switch(tree->mmod->model_id)
    {
    case SLFV_UNIFORM: { return(-1.0); break;}
    case SLFV_GAUSSIAN:  
      { 
        return(POW(tree->mmod->sigsq[0]/(SLFV_Rate_Per_Unit_Area(tree)*4.*PI*tree->mmod->mu),0.25));
        break; 
      }
    }
  return(-1.);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl SLFV_Generation_Length(t_tree *tree)
{
  return(1./(2.*tree->mmod->mu*POW(tree->mmod->rad,2)*PI*SLFV_Rate_Per_Unit_Area(tree)));
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl SLFV_Neighborhood_Size(t_tree *tree)
{
  switch(tree->mmod->model_id)
    {
    case SLFV_UNIFORM: { return(1./tree->mmod->mu); break; }
    case SLFV_GAUSSIAN:  { return(2./tree->mmod->mu); break; }
    }
  return(-1.);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl SLFV_Update_Sigsq(t_tree *tree)
{  
  switch(tree->mmod->model_id)
    {
    case SLFV_UNIFORM : { return(-1.0); break;}
    case SLFV_GAUSSIAN :  
      {
        return(4.*PI*
               SLFV_Rate_Per_Unit_Area(tree) *
               pow(tree->mmod->rad,4)*
               tree->mmod->mu); 
        break; 
      }
    }
  return(-1.);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

/* Uses the regression technique described in Barton et al. TPB (2013)
   to estimate the size of the neighborhood when the population evolve
   according to a spatial Lambda-Fleming-Viot process.
*/
phydbl SLFV_Neighborhood_Size_Regression(t_tree *tree)
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
              PhyML_Fprintf(stderr,"\n. %s",Write_Tree(tree));
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

phydbl SLFV_Lk_Range(t_dsk *young, t_dsk *old, t_tree *tree)
{
  phydbl lnP;

  lnP = UNLIKELY;
  
  switch(tree->mmod->model_id)
    {
    case SLFV_GAUSSIAN :
      {
        lnP = SLFV_Lk_Gaussian_Range(young,old,tree);
        break;
      }
    case SLFV_UNIFORM :
      {
        PhyML_Fprintf(stderr,"\n. SLFV model with rectangle is not implemented. Sorry...");
        assert(false);
        break;
      }      
    default : break;
    }

  return(lnP);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl SLFV_Lk_Gaussian_Range(t_dsk *young, t_dsk *old, t_tree *tree)
{
  t_dsk *disk;
  phydbl lnL,dt,disk_lnL;
  int n_evt;
    
  assert(young);

  lnL  = 0.0;
  lnL += SLFV_Lk_Gaussian_Core(young,tree);

  disk_lnL = 0.0;
  dt = 0.0;
  n_evt = 0;
  disk = young->prev;
  do
    {
      if(disk->age_fixed == NO)
        {
          dt += fabs(disk->next->time - disk->time);
          n_evt++;
        }

      assert(disk);
      if(disk->time > disk->next->time) return UNLIKELY;
      disk_lnL = SLFV_Lk_Gaussian_Core(disk,tree);
      lnL += disk_lnL;

      if(disk->next != NULL) disk->cum_glnL = disk->next->cum_glnL + disk_lnL;      

      if(disk == old) break;
      disk = disk->prev;
    }
  while(disk);

  lnL += n_evt * log(tree->mmod->lbda) - tree->mmod->lbda * dt;
 
  return(lnL);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl SLFV_Lk_Gaussian_Core(t_dsk *disk, t_tree *tree)
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

phydbl SLFV_Lk_Gaussian(t_tree *tree)
{  
  t_dsk *disk;
  int n_evt;
    
  assert(!tree->young_disk->next);
  assert(tree->young_disk->prev);
  
  tree->mmod->c_lnL = 0.0;

  if(PHYREX_Total_Number_Of_Intervals(tree) > tree->mmod->max_num_of_intervals)
   {
     tree->mmod->c_lnL = UNLIKELY;
     return UNLIKELY;
   }
  
  if(isinf(tree->mmod->c_lnL) || isnan(tree->mmod->c_lnL)) 
    {
      tree->mmod->c_lnL = UNLIKELY;
      return tree->mmod->c_lnL;
    }

  /* tree->mmod->c_lnL += PHYREX_LnPrior_Radius(tree); */
  
  // !!!!!!! Really necessary?
  PHYREX_Update_Lindisk_List(tree);

  tree->mmod->c_lnL += SLFV_Lk_Gaussian_Core(tree->young_disk,tree);
  
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

      tree->mmod->c_lnL += SLFV_Lk_Gaussian_Core(disk,tree);
      if(disk->prev == NULL) break;
      disk = disk->prev;      
    }
  while(1);

  /* tree->mmod->c_lnL += TIMES_Lk_SLFV(tree); */
  /* tree->mmod->c_lnL += n_evt * log(tree->mmod->lbda) - tree->mmod->lbda * fabs(tree->young_disk->time - disk->time); */
    
  if(isinf(tree->mmod->c_lnL) || isnan(tree->mmod->c_lnL)) tree->mmod->c_lnL = UNLIKELY;

  return(tree->mmod->c_lnL);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
// Simulate Etheridge-Barton model forwards in time, following n_otu lineages
// on a rectangle of dimension width x h
t_sarea *SLFV_Simulate_Forward_Core(int n_sites, t_tree *tree)
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
              
              switch(mmod->model_id)
                {
                case SLFV_UNIFORM:
                  {
                    if(PHYREX_Is_In_Disk(disk->ldsk_a[i]->coord,new_disk,mmod) == YES)
                      parent_prob[i] = 1.0;
                    else
                      parent_prob[i] = 0.0;
                    break;
                  }
                case SLFV_GAUSSIAN:
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
              switch(mmod->model_id)
                {
                case SLFV_UNIFORM:
                  {
                    if(PHYREX_Is_In_Disk(disk->ldsk_a[i]->coord,new_disk,mmod) == YES)
                      prob_death = mmod->mu;
                    break;
                  }
                case SLFV_GAUSSIAN:
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
              switch(mmod->model_id)
                {
                case SLFV_UNIFORM: { PHYREX_Runif_Rectangle_Overlap(new_ldsk,new_disk,mmod); break; }
                case SLFV_GAUSSIAN:  { PHYREX_Rnorm_Trunc(new_ldsk,new_disk,mmod); break; }
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
          PHYREX_Free_Ldisk(disk->prev->ldsk);
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


/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
// Simulate Etheridge-Barton model backwards in time, following n_otu lineages
// on a rectangle.
phydbl SLFV_Simulate_Backward_Core(t_dsk *init_disk, int avoid_multiple_mergers, t_tree *tree)
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
          for(j=0;j<mmod->n_dim;++j)
            {
              switch(mmod->model_id)
                {
                case SLFV_UNIFORM : case SLFV_GAUSSIAN : 
                  { 
                    new_disk->centr->lonlat[j] = Uni()*(mmod->lim_up->lonlat[j]-mmod->lim_do->lonlat[j])+mmod->lim_do->lonlat[j];
                    break; 
                  }
                case RW : case RRW_GAMMA : case RRW_LOGNORMAL :  
                  { 
                    new_disk->centr->lonlat[j] = disk->ldsk_a[Rand_Int(0,disk->n_ldsk_a-1)]->coord->lonlat[j];
                    break; 
                  }
                default : break;
                }
            }
          
          for(j=0;j<mmod->n_dim;++j) lnL += log(1./(mmod->lim_up->lonlat[j]-mmod->lim_do->lonlat[j]));
          
          /* Populate new_disk->ldsk_a array */
          PHYREX_Update_Lindisk_List_Core(new_disk,tree);          

          new_disk->ldsk = NULL;
          /* Which lineages in new_disk->ldsk_a are hit? */
          for(i=0;i<new_disk->n_ldsk_a;++i)
            {
              prob_hit = -1.;
              switch(mmod->model_id)
                {
                case SLFV_UNIFORM : 
                  { 
                    prob_hit = mmod->mu; 
                    break; 
                  }
                case SLFV_GAUSSIAN :   
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
                  case RW : case RRW_GAMMA : case RRW_LOGNORMAL :
                    {
                      prob_hit = 1./(new_disk->n_ldsk_a+1.);
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
                      switch(tree->mmod->model_id)
                        {
                        case SLFV_UNIFORM :
                          {
                            PHYREX_Runif_Rectangle_Overlap(new_disk->ldsk,new_disk,tree->mmod);
                            break;
                          }
                        case SLFV_GAUSSIAN : case RW : case RRW_GAMMA : case RRW_LOGNORMAL :
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
              
              if(new_disk->ldsk && new_disk->ldsk->n_next == 2 && avoid_multiple_mergers == YES) break;
            }
        }

      disk = new_disk;
    }
  while(1);  
  disk->prev = NULL;

  return(lnL);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// Simulate Etheridge-Barton model backwards in time, following n_otu lineages
// on a rectangle of dimension width w x h
// See Kelleher, Barton & Etheridge, Bioinformatics, 2013.
t_tree *SLFV_Simulate(int n_otu, int n_sites, phydbl w, phydbl h, phydbl  lbda, phydbl rad, phydbl mu, int r_seed)
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
  io->init_len = 500; /* sequence length */

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

  tree->times = TIMES_Make_Time_Struct(tree->n_otu);
  TIMES_Init_Time_Struct(tree->times,io->times,tree->n_otu);

  tree->data      = cdata;
  tree->mod       = mod;
  tree->io        = io;
  tree->n_pattern = tree->data->crunch_len;


  /* Allocate migrep model */
  mmod = PHYREX_Make_Migrep_Model(tree->n_otu,n_dim);
  tree->mmod = mmod;
  PHYREX_Init_Migrep_Mod(mmod,n_dim,0.0,0.0,w,h);

  mmod->lbda = lbda;
  mmod->rad = rad;
  mmod->mu = mu;

  if(mmod->lbda > mmod->max_lbda) mmod->lbda = mmod->max_lbda;
  if(mmod->rad > mmod->max_rad) mmod->rad = mmod->max_rad;
  if(mmod->mu > mmod->max_mu) mmod->mu = mmod->max_mu;

  if(mmod->lbda < mmod->min_lbda) mmod->lbda = mmod->min_lbda;
  if(mmod->rad < mmod->min_rad) mmod->rad = mmod->min_rad;
  if(mmod->mu < mmod->min_mu) mmod->mu = mmod->min_mu;


  mmod->sigsq[0] = SLFV_Update_Sigsq(tree);
  
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

  /* /\* !!!!!!!!!!!!!!!!!!!!!!!!!! *\/ */
  /* /\* Fix coordinates of first two tips *\/ */
  /* tree->young_disk->ldsk_a[0]->coord->lonlat[0] = (tree->mmod->lim_up->lonlat[0]+tree->mmod->lim_do->lonlat[0])/2; // longitude */
  /* tree->young_disk->ldsk_a[0]->coord->lonlat[1] = (tree->mmod->lim_up->lonlat[1]+tree->mmod->lim_do->lonlat[1])/2; // latitude */
  /* tree->young_disk->ldsk_a[1]->coord->lonlat[0] = (tree->mmod->lim_up->lonlat[1]+tree->mmod->lim_do->lonlat[1])/2 + 1; // longitude */
  /* tree->young_disk->ldsk_a[1]->coord->lonlat[1] = (tree->mmod->lim_up->lonlat[1]+tree->mmod->lim_do->lonlat[1])/2 + 1; // latitude */

  
  for(i=0;i<tree->n_otu;i++)
    {
      tree->young_disk->ldsk_a[i]->nd = tree->a_nodes[i];
      tree->a_nodes[i]->ldsk = tree->young_disk->ldsk_a[i];
    }
  
  SLFV_Simulate_Backward_Core(tree->young_disk,NO,tree);

  PHYREX_Ldsk_To_Tree(tree);
  Update_Ancestors(tree->n_root,tree->n_root->v[2],tree->n_root->b[2],tree);
  Update_Ancestors(tree->n_root,tree->n_root->v[1],tree->n_root->b[1],tree);
  RATES_Fill_Lca_Table(tree);
  
  
  /* min_rate = 1.E-5; */
  /* max_rate = 1.E-4; */

  phydbl T = PHYREX_Tree_Height(tree);

  tree->rates->bl_from_rt = YES;
  tree->rates->clock_r    = 0.1/fabs(2.*T);
  tree->rates->model_id      = STRICTCLOCK;

  RATES_Update_Edge_Lengths(tree);

  Init_Model(cdata,mod,io);
  Set_Model_Parameters(mod);
  
  tree->mod->ras->n_catg   = 1;
  tree->mod->whichmodel    = HKY85;
  tree->mod->kappa->v      = 4.0;
    
  Make_Tree_For_Pars(tree);
  Make_Tree_For_Lk(tree);
  Make_Spr(tree);
  Evolve(tree->data,tree->mod,0,tree);

  PhyML_Printf("@@@ %G %G %G : ",mmod->lbda,mmod->mu,mmod->rad);
  for(int i=0;i<tree->n_otu-1;++i) for(int j=i+1;j<tree->n_otu;++j) PhyML_Printf("%G ",PHYREX_Dist_Between_Two_Ldsk(tree->a_nodes[i]->ldsk,tree->a_nodes[j]->ldsk,tree));
  PhyML_Printf(" : ");
  for(int i=0;i<tree->n_otu-1;++i) for(int j=i+1;j<tree->n_otu;++j) PhyML_Printf("%G ",Euclidean_Dist(tree->a_nodes[i]->ldsk->coord,tree->a_nodes[j]->ldsk->coord));
  PhyML_Printf("\n");
  for(int i=0;i<tree->n_otu;++i) PhyML_Printf("\n%20s %12G %12G",
                                              tree->a_nodes[i]->name,
                                              tree->a_nodes[i]->ldsk->coord->lonlat[0],
                                              tree->a_nodes[i]->ldsk->coord->lonlat[1]);
  PhyML_Printf("\n\n");
  PhyML_Printf(">> SEQUENCES\n");
  Print_CSeq(stdout,NO,tree->data,tree);
  PhyML_Printf("<< SEQUENCES");
  
  /* Init_Partial_Lk_Tips_Double(tree); */
  /* Init_Partial_Lk_Loc(tree); */

  disk = tree->young_disk->prev;
  while(disk->prev) disk = disk->prev;
  
  PhyML_Printf("\n. Parameter boundaries: lambda:[%G,%G]; mu=[%G,%G]; rad=[%G,%G]",
               mmod->min_lbda,mmod->max_lbda,
               mmod->min_mu,mmod->max_mu,
               mmod->min_rad,mmod->max_rad);
  PhyML_Printf("\n. Useful parameters: lambda=%G; mu=%G; rad=%G; clockr=%G; sigsq=%G",
               mmod->lbda,
               mmod->mu,
               mmod->rad,
               tree->rates->clock_r,
               mmod->sigsq[0]);

  PhyML_Printf("\n. Useful statistics: t.root=%f n.int=%d n.coal=%d n.hit=%d root.x=%f root.y=%f nt.div=%f\n",
               disk->time,
               PHYREX_Total_Number_Of_Intervals(tree),
               PHYREX_Total_Number_Of_Coal_Disks(tree),
               PHYREX_Total_Number_Of_Hit_Disks(tree),
               disk->ldsk->coord->lonlat[0],
               disk->ldsk->coord->lonlat[1],
               Nucleotide_Diversity(tree->data));
  
  Exit("\n");
  
  /* PhyML_Printf("\n. Tree: "); */
  /* PhyML_Printf("\n. %s \n",Write_Tree(tree)); */

  /* PhyML_Printf("\n. Spatial coordinates: "); */
  /* for(i=0;i<tree->n_otu;i++) */
  /*   { */
  /*     PhyML_Printf("\n%15s %12f %12f", */
  /*                  tree->young_disk->ldsk_a[i]->nd->name, */
  /*                  tree->young_disk->ldsk_a[i]->coord->lonlat[0],                    */
  /*                  tree->young_disk->ldsk_a[i]->coord->lonlat[1]); */
  /*   } */
  /* PhyML_Printf("\n"); */
  
  /* for(int i=0;i<io->n_otu;i++) */
  /*   { */
  /*     PhyML_Printf("\n<clade id=\"clad%d\">",i+1); */
  /*     PhyML_Printf("\n\t<taxon value=\"%s\"/>",tree->a_nodes[i]->name); */
  /*     PhyML_Printf("\n</clade>"); */
  /*     PhyML_Printf("\n<calibration id=\"cal%d\">",i+1); */
  /*     PhyML_Printf("\n\t<lower>0.0</lower>"); */
  /*     PhyML_Printf("\n\t<upper>0.0</upper>"); */
  /*     PhyML_Printf("\n\t<appliesto clade.id=\"clad%d\"/>",i+1); */
  /*     PhyML_Printf("\n</calibration>"); */
  /*   } */

  return(tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

t_tree *SLFV_Simulate_Independent_Loci(int n_otu, int n_loci, phydbl w, phydbl h, int r_seed)
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
  io->rates->model_id = STRICTCLOCK;
  RATES_Init_Rate_Struct(tree->rates,io->rates,tree->n_otu);

  tree->times = TIMES_Make_Time_Struct(tree->n_otu);
  TIMES_Init_Time_Struct(tree->times,io->times,tree->n_otu);
  
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
  mmod = PHYREX_Make_Migrep_Model(tree->n_otu,n_dim);
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
  mmod->sigsq[0] = 4.*pow(mmod->rad,4)*mmod->lbda*PI*mmod->mu/area * mmod->gen_cal_time;

 
  phydbl p = SLFV_Prob_Two_Random_Lineages_Coal_One_Event(w,h,mmod->mu,mmod->rad);
  printf("\n. p.sim: %G p.appx: %G",p,4.*pow(mmod->mu,2)*pow(PI,2)*pow(mmod->rad,4)/(pow(w*h,2)));
  phydbl rhoe = 1./(1.-exp(-mmod->lbda*p*mmod->gen_cal_time));
  rhoe /= area;
  
  PhyML_Printf("\n. lambda=%G; mu=%G; rad=%G; clockr=%G; sigsq=%G; neigh=%G; rhoe(small rad)=%G rhoe(numerical)=%G rhoe(large rad)=%G one_gen=%G",
               mmod->lbda,
               mmod->mu,
               mmod->rad,
               tree->rates->clock_r,
               mmod->sigsq[0],
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
      SLFV_Simulate_Backward_Core(tree->young_disk,NO,tree);
  
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

      
      Update_Ancestors(tree->n_root,tree->n_root->v[2],tree->n_root->b[2],tree);
      Update_Ancestors(tree->n_root,tree->n_root->v[1],tree->n_root->b[1],tree);
      RATES_Fill_Lca_Table(tree);


      tree->rates->bl_from_rt = YES;
      /* tree->rates->clock_r    = 1.E-3/FABS(T); // slow */
      /* tree->rates->clock_r    = 1.E-0/FABS(T); // fast */
      /* tree->rates->clock_r    = 1.E-1/FABS(T); // medium */
      /* tree->rates->clock_r = 1.E-7; // slow; */
      tree->rates->clock_r = 1.E-3; // medium;
      /* tree->rates->clock_r = 1.E-2; // fast; */

      
      PhyML_Printf("\n. #!# mutation rate at that locus: %G subst. per time unit",tree->rates->clock_r);
      subst_rate += tree->rates->clock_r;
      
      RATES_Update_Edge_Lengths(tree);

      char *s = Write_Tree(tree);
      PhyML_Printf("\n. #@# tree: %s",s);
      Free(s);

      Evolve(tree->data,tree->mod,locus_idx,tree);
      
      t_dsk *disk = tree->young_disk->prev;
      while(disk->prev)
        {
          disk = disk->prev;
          if(disk->next->ldsk != NULL) PHYREX_Free_Ldisk(disk->next->ldsk);
          PHYREX_Free_Disk(disk->next);
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


/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void SLFV_Integrated_Coal_Rate(t_ldsk *l0, t_ldsk *l1, phydbl T, t_tree *tree)
{
  t_ldsk *baseline;
  phydbl integrated_coal_rate,njumps,coaltime;
  int coal;
  
  PhyML_Printf("\n# ");


  coal = -1;
  baseline = l0;
  integrated_coal_rate = 0.0;
  njumps = 0.0;
  coaltime = 0.0;
  do
    {
      if(l0->prev->disk->time > l1->prev->disk->time)
        {
          assert(coal == -1);
          integrated_coal_rate +=
            tree->mmod->lbda * SLFV_Prob_Two_Lineages_Coal(l0,l1,tree) *
            fabs(baseline->disk->time - MAX(l0->prev->disk->time,l1->prev->disk->time));
          
          njumps+=1.0;
          baseline = l0->prev;

          l0 = l0->prev;
        }
      else if(l0->prev->disk->time < l1->prev->disk->time)
        {
          assert(coal == -1);
          integrated_coal_rate += 
            tree->mmod->lbda * SLFV_Prob_Two_Lineages_Coal(l0,l1,tree) *
            fabs(baseline->disk->time - MAX(l0->prev->disk->time,l1->prev->disk->time));
          
          njumps+=1.0;
          baseline = l1->prev;
          
          l1 = l1->prev;
        }
      else
        {
          assert(l0->prev);
          assert(l1->prev);
          assert(l0->prev == l1->prev);
          coaltime = l0->prev->disk->time;
          if(coaltime > T)
            {
              coal = 1;
            }
          else
            {
              coal = 0;
              integrated_coal_rate += 
                tree->mmod->lbda * SLFV_Prob_Two_Lineages_Coal(l0,l1,tree) *
                fabs(baseline->disk->time - T);              
            }
                        
          l0->prev = NULL;
          l1->prev = NULL;
        }    
    }
  while(l0->prev && l1->prev);
  
  PhyML_Printf("%d %f %g %f",
               coal,
               integrated_coal_rate,
               njumps,
               coaltime);
  PhyML_Printf("\n");
  Exit("\n");
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
// Calculation of Pr(T>n) where T is the proba of coalescence
void SLFV_Sum_Coal_Rate(t_ldsk *l0, t_ldsk *l1, int n, t_tree *tree)
{
  t_dsk *disk;
  phydbl prod_one_min_coal_prob;
  int coal,n_evts;

  
  PhyML_Printf("\n# ");

  prod_one_min_coal_prob = 0.0;
  n_evts = 1;
  coal = 0;
  disk = tree->young_disk->prev;
  do
    {

      if(n_evts <= n) prod_one_min_coal_prob += log(1. - SLFV_Prob_Two_Lineages_Coal(l0,l1,tree));
      
      if(l0->prev == disk->ldsk) l0 = disk->ldsk;
      if(l1->prev == disk->ldsk) l1 = disk->ldsk;
      if(l0 == l1 && n_evts <= n) coal = 1;
      disk = disk->prev;
      n_evts++;
    }
  while(disk);

  prod_one_min_coal_prob = exp(prod_one_min_coal_prob);
  
  PhyML_Printf("%d %d %f",coal,n_evts,prod_one_min_coal_prob);
  PhyML_Printf("\n");
  Exit("\n");
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void SLFV_Generate_Ldsk_New_Location(t_ldsk *l, t_ldsk *prev_l, phydbl rad, phydbl *hr, int dim_idx, t_tree *tree)
{
  int err;
  phydbl new_c,new_l,c;

  err = NO;

  /* Forward move */
  new_c = Uni()*(tree->mmod->lim_up->lonlat[dim_idx]-tree->mmod->lim_do->lonlat[dim_idx])+tree->mmod->lim_do->lonlat[dim_idx];
  (*hr) -= log(1./(tree->mmod->lim_up->lonlat[dim_idx]-tree->mmod->lim_do->lonlat[dim_idx]));
  c = new_c;
  new_l = Rnorm_Trunc(c,rad,tree->mmod->lim_do->lonlat[dim_idx],tree->mmod->lim_up->lonlat[dim_idx],&err);
  (*hr) -= Log_Dnorm_Trunc(new_l,c,rad,tree->mmod->lim_do->lonlat[dim_idx],tree->mmod->lim_up->lonlat[dim_idx],&err);
  
  /* Reverse move */
  c = prev_l->disk->centr->lonlat[dim_idx];
  (*hr) += Log_Dnorm_Trunc(prev_l->coord->lonlat[dim_idx],c,rad,tree->mmod->lim_do->lonlat[dim_idx],tree->mmod->lim_up->lonlat[dim_idx],&err);
  
  l->disk->centr->lonlat[dim_idx] = new_c;
  l->coord->lonlat[dim_idx] = new_l;
  
  assert(err == NO);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
