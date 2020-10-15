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




phydbl LOCATION_Lk(t_tree *tree)
{

  switch(tree->mmod->model_id)
    {
    case SLFV_GAUSSIAN :
      {
        tree->mmod->c_lnL = SLFV_Lk_Gaussian(tree);
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
        tree->mmod->c_lnL = RRW_Lk(tree);
        break;
      }
    case RW :
      {
        tree->mmod->c_lnL = RW_Lk(tree);
        break;
      }      
    default : assert(false);
    }

  return(tree->mmod->c_lnL);
      
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl LOCATION_Lk_Range(t_dsk *young, t_dsk *old, t_tree *tree)
{    
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
        return(RW_Lk_Range(young,old,tree));
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
