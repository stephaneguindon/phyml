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

phydbl LOCATION_Lk(t_tree *tree)
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
        lnP = SLFV_Lk_Gaussian(tree);
        break;
      }
    case SLFV_UNIFORM :
      {
        PhyML_Fprintf(stderr,"\n. SLFV model with rectangle is not implemented. Sorry...");
        assert(false);
        break;
      }
    case RRW :
      {
        lnP = RRW_Lk(tree);
        break;
      }
    case RW :
      {
        lnP = RW_Lk(tree);
        break;
      }      
    default : assert(false);
    }

  if(isinf(lnP) || isnan(lnP)) lnP = UNLIKELY;
  
  tree->mmod->c_lnL = lnP;
  
  return(tree->mmod->c_lnL);
      
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
    case RW : case RRW : 
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
    case RW : case RRW : 
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
    case RW : case RRW : 
        {
          assert(false);
          break;
        }
    default : assert(FALSE);
    }
}
