/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

/* Routines that implement the spatial part of Etheridge and Barton's model of continuous-space
   coalescent and other models of migration/spatial dispersion.
*/

#include "location.h"

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

#ifdef PHYREX

phydbl LOCATION_Lk(t_tree *tree)
{
  phydbl lnL;

  lnL = UNLIKELY;
  
  if(tree->mmod->use_locations == NO)
    {
      tree->mmod->c_lnL = 0.0;
      return(tree->mmod->c_lnL);
    }

  switch(tree->mmod->model_id)
    {
    case SLFV_GAUSSIAN :
      {
        lnL = SLFV_Lk_Gaussian(tree);
        break;
      }
    case SLFV_UNIFORM :
      {
        PhyML_Fprintf(stderr,"\n. SLFV model with rectangle is not implemented. Sorry...");
        assert(false);
        break;
      }
    case RRW_GAMMA : case RRW_LOGNORMAL :
      {
        lnL = RRW_Lk(tree);
        break;
      }
    case RW :
      {
        lnL = RW_Lk(tree);
        break;
      }      
    default : assert(false);
    }

  if(isinf(lnL) || isnan(lnL)) lnL = UNLIKELY;
  
  tree->mmod->c_lnL = lnL;
  
  return(tree->mmod->c_lnL);
      
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl LOCATION_Prior(t_tree *tree)
{
  phydbl lnP;

  lnP = 0.0;
  
  if(tree->mmod->use_locations == NO)
    {
      tree->mmod->c_lnL = 0.0;
      return(tree->mmod->c_lnL);
    }

  switch(tree->mmod->model_id)
    {
    case SLFV_GAUSSIAN :
      {
        lnP = UNLIKELY;
        break;
      }
    case SLFV_UNIFORM :
      {
        lnP = UNLIKELY;
        break;
      }
    case RRW_GAMMA : case RRW_LOGNORMAL :
      {
        lnP = RRW_Prior(tree);
        break;
      }
    case RW :
      {
        lnP = RW_Prior(tree);
        break;
      }
    default : assert(false);
    }

  if(isinf(lnP) || isnan(lnP)) lnP = UNLIKELY;
  
  tree->mmod->c_lnP = lnP;
  
  return(tree->mmod->c_lnP);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl LOCATION_Lk_Range(t_dsk *young, t_dsk *old, t_tree *tree)
{    
  if(tree->mmod->use_locations == NO) return(0.0);

  switch(tree->mmod->model_id)
    {
    case SLFV_GAUSSIAN : 
      {
        return(SLFV_Lk_Gaussian_Range(young,old,tree));
        break;
      }
    case SLFV_UNIFORM :
      {
        PhyML_Fprintf(stderr,"\n. SLFV model with rectangle is not implemented. Sorry...");
        assert(false);
        break;
      }
    case RW : case RRW_GAMMA : case RRW_LOGNORMAL : 
      {
        return(RRW_Lk_Range(young,old,tree));
        break;
      }
    default : assert(FALSE);
    }
  
  return(-1.);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl LOCATION_Lk_Path(t_dsk *young, t_dsk *old, t_tree *tree)
{    
  if(tree->mmod->use_locations == NO) return(0.0);

  switch(tree->mmod->model_id)
    {
    case SLFV_GAUSSIAN : 
      {
        return(SLFV_Lk_Gaussian_Range(young,old,tree));
        break;
      }
    case SLFV_UNIFORM :
      {
        PhyML_Fprintf(stderr,"\n. SLFV model with rectangle is not implemented. Sorry...");
        assert(false);
        break;
      }
    case RW : case RRW_GAMMA : case RRW_LOGNORMAL : 
      {
        return(RRW_Forward_Lk_Path(old->ldsk,young->ldsk,tree));
        break;
      }
    default : assert(FALSE);
    }
  
  return(-1.);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl LOCATION_Forward_Lk_Path(t_ldsk *a, t_ldsk *d, t_tree *tree)
{
  phydbl lnL;

  lnL = UNLIKELY;
  switch(tree->mmod->model_id)
    {
    case SLFV_GAUSSIAN : 
      {
        assert(false);
        break;
      }
    case SLFV_UNIFORM :
      {
        assert(false);
        break;
      }
    case RW : case RRW_GAMMA : case RRW_LOGNORMAL : 
        {
          lnL = RRW_Forward_Lk_Path(a,d,tree);
          break;
        }
    default : assert(FALSE);
    }
  
  return lnL;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl LOCATION_Lk_Core(t_dsk *disk, t_tree *tree)
{
  phydbl lnL;
  
  switch(tree->mmod->model_id)
    {
    case SLFV_GAUSSIAN : 
      {
        lnL = SLFV_Lk_Gaussian_Core(disk,tree);
        break;
      }
    case SLFV_UNIFORM :
      {
        PhyML_Fprintf(stderr,"\n. SLFV model with rectangle is not implemented. Sorry...");
        assert(false);
        break;
      }
    case RW : case RRW_GAMMA : case RRW_LOGNORMAL : 
        {
          lnL = RRW_Lk_Core(disk,tree);
          break;
        }
    default : assert(FALSE);
    }
  
  return lnL;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void LOCATION_Sample_Path(t_ldsk *young, t_ldsk *old, phydbl *sd, phydbl *global_hr, t_tree *tree)
{
  switch(tree->mmod->model_id)
    {
    case SLFV_GAUSSIAN : 
      {
        SLFV_Sample_Path(young,old,sd,global_hr,tree);
        break;
      }
    case SLFV_UNIFORM :
      {
        PhyML_Fprintf(stderr,"\n. SLFV model with rectangle is not implemented. Sorry...");
        assert(false);
        break;
      }
    case RW : case RRW_GAMMA : case RRW_LOGNORMAL : 
        {
          assert(false);
          break;
        }
    default : assert(FALSE);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

char *LOCATION_Model_Id(t_phyrex_mod *mmod)
{
  char *s;

  s = (char *)mCalloc(T_MAX_LINE,sizeof(char));
  strcpy(s,"none");
  
  switch(mmod->model_id)
    {
    case SLFV_GAUSSIAN :
      {
        strcpy(s,"spatial coalescent (Gaussian)");
        break;
      }
    case SLFV_UNIFORM :
      {
        strcpy(s,"spatial coalescent (Rectangle)");
        break;
      }
    case RW :
      {
        strcpy(s,"strict random walk");
        break;
      }
    case RRW_GAMMA :
      {
        strcpy(s,"relaxed random walk (Gamma)");
        break;
      }
    case RRW_LOGNORMAL :
      {
        strcpy(s,"relaxed random walk (Lognormal)");
        break;
      }
    default : assert(FALSE);
    }
  
  return(s);
 
}

#endif
