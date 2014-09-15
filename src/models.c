/*

PhyML :  a program that  computes maximum likelihood  phyLOGenies from
DNA or AA homoLOGous sequences

Copyright (C) Stephane Guindon. Oct 2003 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/

#include "models.h"

#ifdef BEAGLE
#include "beagle_utils.h"
#endif

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* Handle any number of states (>1) */
void PMat_JC69(phydbl l, int pos, phydbl *Pij, t_mod *mod)
{
  int ns;
  int i,j;

  ns = mod->ns;


  For(i,ns) Pij[pos+ ns*i+i] = 1. - ((ns - 1.)/ns)*(1. - EXP(-ns*l/(ns - 1.)));
  For(i,ns-1)
    for(j=i+1;j<ns;j++)
      {
        Pij[pos+ ns*i+j] = (1./ns)*(1. - EXP(-ns*l/(ns - 1.)));
        if(Pij[pos+ns*i+j] < SMALL_PIJ) Pij[pos+ns*i+j] = SMALL_PIJ;
        Pij[pos+ ns*j+i] = Pij[pos+ ns*i+j];
      }
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void PMat_K80(phydbl l, phydbl kappa, int pos, phydbl *Pij)
{
  phydbl Ts,Tv,e1,e2,aux;
  int i,j;
  /*0 => A*/
  /*1 => C*/
  /*2 => G*/
  /*3 => T*/

  /* Ts -> transition*/
  /* Tv -> transversion*/


  aux = -2*l/(kappa+2);
  e1 = (phydbl)EXP(aux *2);

  e2 = (phydbl)EXP(aux *(kappa+1));
  Tv = .25*(1-e1);
  Ts = .25*(1+e1-2*e2);

  Pij[pos+ 4*0+0] = Pij[pos+ 4*1+1] =
  Pij[pos+ 4*2+2] = Pij[pos+ 4*3+3] = 1.-Ts-2.*Tv;

  Pij[pos+ 4*0+1] = Pij[pos+ 4*1+0] = Tv;
  Pij[pos+ 4*0+2] = Pij[pos+ 4*2+0] = Ts;
  Pij[pos+ 4*0+3] = Pij[pos+ 4*3+0] = Tv;

  Pij[pos+ 4*1+2] = Pij[pos+ 4*2+1] = Tv;
  Pij[pos+ 4*1+3] = Pij[pos+ 4*3+1] = Ts;

  Pij[pos+ 4*2+3] = Pij[pos+ 4*3+2] = Tv;

  For(i,4) For(j,4)
    if(Pij[pos + 4*i+j] < SMALL_PIJ) Pij[pos + 4*i+j] = SMALL_PIJ;

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////



void PMat_TN93(phydbl l, t_mod *mod, int pos, phydbl *Pij)
{
  int i,j;
  phydbl e1,e2,e3;
  phydbl a1t,a2t,bt;
  phydbl A,C,G,T,R,Y;
  phydbl kappa1,kappa2;

  A = mod->e_frq->pi->v[0]; C = mod->e_frq->pi->v[1]; G = mod->e_frq->pi->v[2]; T = mod->e_frq->pi->v[3];
  R = A+G;  Y = T+C;

  if(mod->kappa->v < .0) mod->kappa->v = 1.0e-5;

  if((mod->whichmodel != F84) && (mod->whichmodel != TN93)) mod->lambda->v = 1.;
  else if(mod->whichmodel == F84)
    {
      mod->lambda->v = Get_Lambda_F84(mod->e_frq->pi->v,&mod->kappa->v);
    }

  kappa2 = mod->kappa->v*2./(1.+mod->lambda->v);
  kappa1 = kappa2 * mod->lambda->v;


  bt = l/(2.*(A*G*kappa1+C*T*kappa2+R*Y));

  a1t = kappa1;
  a2t = kappa2;
  a1t*=bt; a2t*=bt;

  e1 = (phydbl)EXP(-a1t*R-bt*Y);
  e2 = (phydbl)EXP(-a2t*Y-bt*R);
  e3 = (phydbl)EXP(-bt);


  /*A->A*/Pij[pos + 4*0+0] = A+Y*A/R*e3+G/R*e1;
  /*A->C*/Pij[pos + 4*0+1] = C*(1-e3);
  /*A->G*/Pij[pos + 4*0+2] = G+Y*G/R*e3-G/R*e1;
  /*A->T*/Pij[pos + 4*0+3] = T*(1-e3);

  /*C->A*/Pij[pos + 4*1+0] = A*(1-e3);
  /*C->C*/Pij[pos + 4*1+1] = C+R*C/Y*e3+T/Y*e2;
  /*C->G*/Pij[pos + 4*1+2] = G*(1-e3);
  /*C->T*/Pij[pos + 4*1+3] = T+R*T/Y*e3-T/Y*e2;

  /*G->A*/Pij[pos + 4*2+0] = A+Y*A/R*e3-A/R*e1;
  /*G->C*/Pij[pos + 4*2+1] = C*(1-e3);
  /*G->G*/Pij[pos + 4*2+2] = G+Y*G/R*e3+A/R*e1;
  /*G->T*/Pij[pos + 4*2+3] = T*(1-e3);

  /*T->A*/Pij[pos + 4*3+0] = A*(1-e3);
  /*T->C*/Pij[pos + 4*3+1] = C+R*C/Y*e3-C/Y*e2;
  /*T->G*/Pij[pos + 4*3+2] = G*(1-e3);
  /*T->T*/Pij[pos + 4*3+3] = T+R*T/Y*e3+C/Y*e2;

  For(i,4) For(j,4)
    if(Pij[pos + 4*i+j] < SMALL_PIJ) Pij[pos + 4*i+j] = SMALL_PIJ;

/*   /\*A->A*\/(*Pij)[0][0] = A+Y*A/R*e3+G/R*e1;  */
/*   /\*A->C*\/(*Pij)[0][1] = C*(1-e3); */
/*   /\*A->G*\/(*Pij)[0][2] = G+Y*G/R*e3-G/R*e1; */
/*   /\*A->T*\/(*Pij)[0][3] = T*(1-e3); */

/*   /\*C->A*\/(*Pij)[1][0] = A*(1-e3); */
/*   /\*C->C*\/(*Pij)[1][1] = C+R*C/Y*e3+T/Y*e2; */
/*   /\*C->G*\/(*Pij)[1][2] = G*(1-e3); */
/*   /\*C->T*\/(*Pij)[1][3] = T+R*T/Y*e3-T/Y*e2; */

/*   /\*G->A*\/(*Pij)[2][0] = A+Y*A/R*e3-A/R*e1; */
/*   /\*G->C*\/(*Pij)[2][1] = C*(1-e3); */
/*   /\*G->G*\/(*Pij)[2][2] = G+Y*G/R*e3+A/R*e1; */
/*   /\*G->T*\/(*Pij)[2][3] = T*(1-e3); */

/*   /\*T->A*\/(*Pij)[3][0] = A*(1-e3); */
/*   /\*T->C*\/(*Pij)[3][1] = C+R*C/Y*e3-C/Y*e2; */
/*   /\*T->G*\/(*Pij)[3][2] = G*(1-e3); */
/*   /\*T->T*\/(*Pij)[3][3] = T+R*T/Y*e3+C/Y*e2; */

/*   For(i,4) For(j,4) */
/*     if((*Pij)[i][j] < SMALL) (*Pij)[i][j] = SMALL; */

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Get_Lambda_F84(phydbl *pi, phydbl *kappa)
{
  phydbl A,C,G,T,R,Y,lambda;
  int kappa_has_changed;

  A = pi[0]; C = pi[1]; G = pi[2]; T = pi[3];
  R = A+G;  Y = T+C;

  if(*kappa < .0) *kappa = 1.0e-5;

  kappa_has_changed = NO;

  do
    {
      lambda = (Y+(R-Y)/(2.*(*kappa)))/(R-(R-Y)/(2.*(*kappa)));

      if(lambda < .0)
    {
      *kappa += *kappa/10.;
      kappa_has_changed = YES;
    }
    }while(lambda < .0);

  if(kappa_has_changed)
    {
      PhyML_Printf("\n. WARNING: This transition/transversion ratio\n");
      PhyML_Printf("  is impossible with these base frequencies!\n");
      PhyML_Printf("  The ratio is now set to %.3f\n",*kappa);
    }

  return lambda;
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////



/********************************************************************/
/* void PMat_Empirical(phydbl l, t_mod *mod, phydbl ***Pij)         */
/*                                                                  */
/* Computes the substitution probability matrix                     */
/* from the initial substitution rate matrix and frequency vector   */
/* and one specific branch length                                   */
/*                                                                  */
/* input : l , branch length                                        */
/* input : mod , choosen model parameters, qmat and pi              */
/* ouput : Pij , substitution probability matrix                    */
/*                                                                  */
/* matrix P(l) is computed as follows :                             */
/* P(l) = EXP(Q*t) , where :                                        */
/*                                                                  */
/*   Q = substitution rate matrix = Vr*D*inverse(Vr) , where :      */
/*                                                                  */
/*     Vr = right eigenvector matrix for Q                          */
/*     D  = diagonal matrix of eigenvalues for Q                    */
/*                                                                  */
/*   t = time interval = l / mr , where :                           */
/*                                                                  */
/*     mr = mean rate = branch length/time interval                 */
/*        = sum(i)(pi[i]*p(i->j)) , where :                         */
/*                                                                  */
/*       pi = state frequency vector                                */
/*       p(i->j) = subst. probability from i to a different state   */
/*               = -Q[ii] , as sum(j)(Q[ij]) +Q[ii] =0              */
/*                                                                  */
/* the Taylor development of EXP(Q*t) gives :                       */
/* P(l) = Vr*EXP(D*t)        *inverse(Vr)                           */
/*      = Vr*POW(EXP(D/mr),l)*inverse(Vr)                           */
/*                                                                  */
/* for performance we compute only once the following matrixes :    */
/* Vr, inverse(Vr), EXP(D/mr)                                       */
/* thus each time we compute P(l) we only have to :                 */
/* make 20 times the operation POW()                                */
/* make 2 20x20 matrix multiplications , that is :                  */
/*   16000 = 2x20x20x20 times the operation *                       */
/*   16000 = 2x20x20x20 times the operation +                       */
/*   which can be reduced to (the central matrix being diagonal) :  */
/*   8400 = 20x20 + 20x20x20 times the operation *                  */
/*   8000 = 20x20x20 times the operation +                          */
/********************************************************************/

void PMat_Empirical(phydbl l, t_mod *mod, int pos, phydbl *Pij)
{
  int n = mod->ns;
  int i, j, k;
  phydbl *U,*V,*R;
  phydbl *expt;
  phydbl *uexpt;

  expt  = mod->eigen->e_val_im;
  uexpt = mod->eigen->r_e_vect_im;
  U     = mod->eigen->r_e_vect;
  V     = mod->eigen->l_e_vect;
  R     = mod->eigen->e_val; /* exponential of the eigen value matrix */

  //Initialize a rate-specific N*N matrix
  For(i,n) For(k,n) Pij[pos+mod->ns*i+k] = .0;

  /* compute POW(EXP(D/mr),l) into mat_eDmrl */
  For(k,n) expt[k] = (phydbl)POW(R[k],l);

  /* multiply Vr*POW(EXP(D/mr),l)*Vi into Pij */
  For (i,n) For (k,n) uexpt[i*n+k] = U[i*n+k] * expt[k];

  For (i,n)
    {
      For (j,n)
    {
      For(k,n)
        {
          Pij[pos+mod->ns*i+j] += (uexpt[i*n+k] * V[k*n+j]);
        }
/* 	  if(Pij[pos+mod->ns*i+j] < SMALL) Pij[pos+mod->ns*i+j] = SMALL; */
      if(Pij[pos+mod->ns*i+j] < SMALL_PIJ) Pij[pos+mod->ns*i+j] = SMALL_PIJ;
    }

#ifndef PHYML
      phydbl sum;
      sum = .0;
      For (j,n) sum += Pij[pos+mod->ns*i+j];
      if((sum > 1.+.0001) || (sum < 1.-.0001))
    {
      PhyML_Printf("\n");
      PhyML_Printf("\n. Q\n");
      For(i,n) { For(j,n) PhyML_Printf("%7.3f ",mod->eigen->q[i*n+j]); PhyML_Printf("\n"); }
      PhyML_Printf("\n. U\n");
      For(i,n) { For(j,n) PhyML_Printf("%7.3f ",U[i*n+j]); PhyML_Printf("\n"); }
      PhyML_Printf("\n");
      PhyML_Printf("\n. V\n");
      For(i,n) { For(j,n) PhyML_Printf("%7.3f ",V[i*n+j]); PhyML_Printf("\n"); }
      PhyML_Printf("\n");
      PhyML_Printf("\n. Eigen\n");
      For(i,n)  PhyML_Printf("%E ",expt[i]);
      PhyML_Printf("\n");
      PhyML_Printf("\n. Pij\n");
      For(i,n) { For (j,n) PhyML_Printf("%f ",Pij[pos+mod->ns*i+j]); PhyML_Printf("\n"); }
      PhyML_Printf("\n. sum = %f",sum);
      if(mod->m4mod)
        {
          int i;
          PhyML_Printf("\n. mod->m4mod->alpha = %f",mod->m4mod->alpha);
          PhyML_Printf("\n. mod->m4mod->delta = %f",mod->m4mod->delta);
          For(i,mod->m4mod->n_h)
        {
          PhyML_Printf("\n. mod->m4mod->multipl[%d] = %f",i,mod->m4mod->multipl[i]);
        }
        }
      PhyML_Printf("\n. l=%f",l);
      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
#endif
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void PMat_Gamma(phydbl l, t_mod *mod, int pos, phydbl *Pij)
{
  int n;
  int i, j, k;
  phydbl *U,*V,*R;
  phydbl *expt;
  phydbl *uexpt;
  phydbl shape;


  n     = mod->ns;
  expt  = mod->eigen->e_val_im;
  uexpt = mod->eigen->r_e_vect_im;
  U     = mod->eigen->r_e_vect;
  V     = mod->eigen->l_e_vect;
  R     = mod->eigen->e_val; /* exponential of the eigen value matrix */

  if(mod->ras->n_catg == 1) shape = 1.E+4;
  else                 shape = mod->ras->alpha->v;


  For(i,n) For(k,n) Pij[pos+mod->ns*i+k] = .0;

  if(shape < 1.E-10)
    {
      PhyML_Printf("\n== Err in file %s at line %d\n\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  /* Formula 13.42, page 220 of Felsenstein's book ``Inferring Phylogenies'' */
  For(k,n) expt[k] = POW(shape/(shape-LOG(R[k])*l),shape);

  /* multiply Vr*expt*Vi into Pij */
  For(i,n) For(k,n) uexpt[i*n+k] = U[i*n+k] * expt[k];

  For (i,n)
    {
      For (j,n)
    {
      For(k,n)
        {
          Pij[pos+mod->ns*i+j] += (uexpt[i*n+k] * V[k*n+j]);
        }
      if(Pij[pos+mod->ns*i+j] < SMALL_PIJ) Pij[pos+mod->ns*i+j] = SMALL_PIJ;
    }

#ifdef DEBUG
      phydbl sum;
      sum = .0;
      For (j,n) sum += Pij[pos+mod->ns*i+j];
      if((sum > 1.+.0001) || (sum < 1.-.0001))
    {
      PhyML_Printf("\n");
      PhyML_Printf("\n. Q\n");
      For(i,n) { For(j,n) PhyML_Printf("%7.3f ",mod->eigen->q[i*n+j]); PhyML_Printf("\n"); }
      PhyML_Printf("\n. U\n");
      For(i,n) { For(j,n) PhyML_Printf("%7.3f ",U[i*n+j]); PhyML_Printf("\n"); }
      PhyML_Printf("\n");
      PhyML_Printf("\n. V\n");
      For(i,n) { For(j,n) PhyML_Printf("%7.3f ",V[i*n+j]); PhyML_Printf("\n"); }
      PhyML_Printf("\n");
      PhyML_Printf("\n. Eigen\n");
      For(i,n)  PhyML_Printf("%E ",expt[i]);
      PhyML_Printf("\n");
      PhyML_Printf("\n. Pij\n");
      For(i,n) { For (j,n) PhyML_Printf("%f ",Pij[pos+mod->ns*i+j]); PhyML_Printf("\n"); }
      PhyML_Printf("\n. sum = %f",sum);
      if(mod->m4mod)
        {
          int i;
          PhyML_Printf("\n. mod->m4mod->ras->alpha = %f",mod->m4mod->alpha);
          PhyML_Printf("\n. mod->m4mod->delta = %f",mod->m4mod->delta);
          For(i,mod->m4mod->n_h)
        {
          PhyML_Printf("\n. mod->m4mod->multipl[%d] = %f",i,mod->m4mod->multipl[i]);
        }
        }
      PhyML_Printf("\n. l=%f",l);
      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
#endif
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void PMat_Zero_Br_Len(t_mod *mod, int pos, phydbl *Pij)
{
  int n = mod->ns;
  int i, j;

  For (i,n) For (j,n) Pij[pos+mod->ns*i+j] = .0;
  For (i,n) Pij[pos+mod->ns*i+i] = 1.0;

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/*
 *Update a rate specific Transition Prob matrix for a given branch-length(already adjusted
 *with the rate prior to this function being called)
 *
 *  Pij: the P-matrix that will be adjusted
 *  l  : branch length * rate
 *  pos: offset into a specific rate-category
 *
 */
void PMat(phydbl l, t_mod *mod, int pos, phydbl *Pij)
{
  /* Warning: l is never the log of branch length here */
  if(l < 0.0)
    {
#ifdef BEAGLE
      Warn_And_Exit(TODO_BEAGLE);
#endif
      PMat_Zero_Br_Len(mod,pos,Pij);
    }
  else
    {
      switch(mod->io->datatype)
        {
        case NT :
          {
            if(mod->use_m4mod)
              {
                PMat_Empirical(l,mod,pos,Pij);
              }
            else
              {
                if((mod->whichmodel == JC69) ||
                   (mod->whichmodel == K80))
                  {
                    /* 		    PMat_JC69(l,pos,Pij,mod); */
                    PMat_K80(l,mod->kappa->v,pos,Pij);
                  }
                else
                  {
                    if(
                       (mod->whichmodel == F81)   ||
                       (mod->whichmodel == HKY85) ||
                       (mod->whichmodel == F84)   ||
                       (mod->whichmodel == TN93))
                      {
                        PMat_TN93(l,mod,pos,Pij);
                      }
                    else
                      {
#ifdef BEAGLE
                        //when there is no active instance (i.e. when we are building the initial tree)
                        if(UNINITIALIZED == mod->b_inst) {
                            PMat_Empirical(l,mod,pos,Pij);
                        }
#else
                        PMat_Empirical(l,mod,pos,Pij);
#endif
                      }
                  }
                break;
              }
          case AA :
            {
#ifdef BEAGLE
                //when there is no active instance (i.e. when we are building the initial tree)
                if(UNINITIALIZED == mod->b_inst) {
                    PMat_Empirical(l,mod,pos,Pij);
                }
#else
                PMat_Empirical(l,mod,pos,Pij);
#endif
              break;
            }
          default:
            {
              PMat_JC69(l,pos,Pij,mod);
              break;
    /* 	      PhyML_Printf("\n. Not implemented yet.\n"); */
    /* 	      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__); */
    /* 	      Warn_And_Exit(""); */
    /* 	      break; */
            }
          }
        }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


int GetDaa (phydbl *daa, phydbl *pi, char *file_name)
{
/* Get the amino acid distance (or substitution rate) matrix
   (grantham, dayhoff, jones, etc).
*/
   FILE * fdaa;
   int i,j, naa;
   phydbl dmax,dmin;
   phydbl sum;
   double val;

   naa = 20;
   dmax = .0;
   dmin = 1.E+40;

   fdaa = (FILE *)Openfile(file_name,0);

   for(i=0; i<naa; i++)
     for(j=0; j<i; j++)
       {
/* 	 if(fscanf(fdaa, "%lf", &daa[i*naa+j])) Exit("\n. err aaRatefile"); */
     if(fscanf(fdaa, "%lf", &val)) Exit("\n. err aaRatefile");
     daa[i*naa+j] = (phydbl)val;
     daa[j*naa+i]=daa[i*naa+j];
     if (dmax<daa[i*naa+j]) dmax=daa[i*naa+j];
     if (dmin>daa[i*naa+j]) dmin=daa[i*naa+j];
       }

   For(i,naa)
     {
/*        if(fscanf(fdaa,"%lf",&pi[i])!=1) Exit("\n. err aaRatefile"); */
       if(fscanf(fdaa,"%lf",&val)!=1) Exit("\n. err aaRatefile");
       pi[i] = (phydbl)val;
     }
   sum = 0.0;
   For(i, naa) sum += pi[i];
   if (FABS(1-sum)>1e-4) {
     PhyML_Printf("\nSum of freq. = %.6f != 1 in aaRateFile\n",sum);
     exit(-1);
   }

   fclose (fdaa);

   return (0);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Update_Qmat_Generic(phydbl *rr, phydbl *pi, int ns, phydbl *qmat)
{
  int i,j;
  phydbl sum,mr;

  For(i,ns*ns) qmat[i] = .0;

  if(rr[(int)(ns*(ns-1)/2)-1] < 0.00001)
    {
      PhyML_Printf("\n== rr[%d]=%f",(int)(ns*(ns-1)/2)-1,rr[(int)(ns*(ns-1)/2)-1]);
      PhyML_Printf("\n== Err in file %s at line %d\n\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  /* PhyML_Printf("\n"); */
  /* For(i,(int)(ns*(ns-1)/2)) */
  /*   { */
  /*     PhyML_Printf("\n> rr %d = %f",i,rr[i]); */
  /*   } */

  For(i,(int)(ns*(ns-1)/2))
    {
      rr[i] /= rr[(int)(ns*(ns-1)/2)-1];
    }

  /* Fill the non-diagonal parts */
  For(i,ns)
    {
      for(j=i+1;j<ns;j++)
    {
      qmat[i*ns+j] = rr[MIN(i,j) * ns + MAX(i,j) -
                (MIN(i,j)+1+(int)POW(MIN(i,j)+1,2))/2];
      qmat[j*ns+i] = qmat[i*ns+j];
    }
    }


  /* Multiply by pi */
  For(i,ns)
    {
      For(j,ns)
    {
      qmat[i*ns+j] *= pi[j];
    }
    }

  /* Compute diagonal elements */
  mr = .0;
  For(i,ns)
    {
      sum = .0;
      For(j,ns) {sum += qmat[i*ns+j];}
      qmat[i*ns+i] = -sum;
      mr += sum * pi[i];
    }

  /* For(i,ns) For(j,ns) qmat[i*ns+j] /= mr; */
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Update_Qmat_GTR(phydbl *rr, phydbl *rr_val, int *rr_num, phydbl *pi, phydbl *qmat)
{
  int i;
  phydbl mr;


  For(i,6) rr[i] = rr_val[rr_num[i]];
  For(i,6)
    if(rr[i] < 0.0)
      {
        PhyML_Printf("\n== rr%d: %f",i,rr[i]);
        PhyML_Printf("\n== Err. in file %s at line %d (function '%s').\n",__FILE__,__LINE__,__FUNCTION__);
        Exit("");
      }

  For(i,6) rr[i] /= rr[5];
  For(i,6) if(rr[i] < RR_MIN) rr[i] = RR_MIN;
  For(i,6) if(rr[i] > RR_MAX) rr[i] = RR_MAX;

  qmat[0*4+1] = (rr[0]*pi[1]);
  qmat[0*4+2] = (rr[1]*pi[2]);
  qmat[0*4+3] = (rr[2]*pi[3]);

  qmat[1*4+0] = (rr[0]*pi[0]);
  qmat[1*4+2] = (rr[3]*pi[2]);
  qmat[1*4+3] = (rr[4]*pi[3]);

  qmat[2*4+0] = (rr[1]*pi[0]);
  qmat[2*4+1] = (rr[3]*pi[1]);
  qmat[2*4+3] = (rr[5]*pi[3]);

  qmat[3*4+0] = (rr[2]*pi[0]);
  qmat[3*4+1] = (rr[4]*pi[1]);
  qmat[3*4+2] = (rr[5]*pi[2]);

  qmat[0*4+0] = -(rr[0]*pi[1]+rr[1]*pi[2]+rr[2]*pi[3]);
  qmat[1*4+1] = -(rr[0]*pi[0]+rr[3]*pi[2]+rr[4]*pi[3]);
  qmat[2*4+2] = -(rr[1]*pi[0]+rr[3]*pi[1]+rr[5]*pi[3]);
  qmat[3*4+3] = -(rr[2]*pi[0]+rr[4]*pi[1]+rr[5]*pi[2]);

  /* compute diagonal terms of Q and mean rate mr = l/t */
  mr = .0;
  For (i,4) mr += pi[i] * (-qmat[i*4+i]);
  For(i,16) qmat[i] /= mr;

  /* int j; */
  /* printf("\n"); */
  /* printf("\n. rr -- "); */
  /* For(i,5) printf(" %15f ",rr[i]); */
  /* printf("\n"); */
  /* For(i,4) */
  /*   { */
  /*     printf("\n. [%15f] \t ",pi[i]); */
  /*     For(j,4) */
  /* 	{ */
  /* 	  printf(" %15f ",qmat[i*4+j]); */
  /* 	} */
  /*   } */
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Update_Qmat_HKY(phydbl kappa, phydbl *pi, phydbl *qmat)
{
  int i;
  phydbl mr;

  /* A -> C */ qmat[0*4+1] = (phydbl)(pi[1]);
  /* A -> G */ qmat[0*4+2] = (phydbl)(kappa*pi[2]);
  /* A -> T */ qmat[0*4+3] = (phydbl)(pi[3]);

  /* C -> A */ qmat[1*4+0] = (phydbl)(pi[0]);
  /* C -> G */ qmat[1*4+2] = (phydbl)(pi[2]);
  /* C -> T */ qmat[1*4+3] = (phydbl)(kappa*pi[3]);

  /* G -> A */ qmat[2*4+0] = (phydbl)(kappa*pi[0]);
  /* G -> C */ qmat[2*4+1] = (phydbl)(pi[1]);
  /* G -> T */ qmat[2*4+3] = (phydbl)(pi[3]);

  /* T -> A */ qmat[3*4+0] = (phydbl)(pi[0]);
  /* T -> C */ qmat[3*4+1] = (phydbl)(kappa*pi[1]);
  /* T -> G */ qmat[3*4+2] = (phydbl)(pi[2]);

  qmat[0*4+0] = (phydbl)(-(qmat[0*4+1]+qmat[0*4+2]+qmat[0*4+3]));
  qmat[1*4+1] = (phydbl)(-(qmat[1*4+0]+qmat[1*4+2]+qmat[1*4+3]));
  qmat[2*4+2] = (phydbl)(-(qmat[2*4+0]+qmat[2*4+1]+qmat[2*4+3]));
  qmat[3*4+3] = (phydbl)(-(qmat[3*4+0]+qmat[3*4+1]+qmat[3*4+2]));

  /* compute diagonal terms of Q and mean rate mr = l/t */
  mr = .0;
  For (i,4) mr += pi[i] * (-qmat[i*4+i]);
  For(i,16) qmat[i] /= mr;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Translate_Custom_Mod_String(t_mod *mod)
{
  int i,j;

  For(i,6) mod->r_mat->n_rr_per_cat->v[i] = 0;

  mod->r_mat->n_diff_rr = 0;

  For(i,6)
    {
      For(j,i)
    {
      if((mod->custom_mod_string->s[i] == mod->custom_mod_string->s[j]))
        {
          break;
        }
    }

      if(i == j)
    {
      mod->r_mat->rr_num->v[i] = mod->r_mat->n_diff_rr;
      mod->r_mat->n_diff_rr++;
    }
      else
    {
      mod->r_mat->rr_num->v[i] = mod->r_mat->rr_num->v[j];
    }

      mod->r_mat->n_rr_per_cat->v[mod->r_mat->rr_num->v[j]]++;
    }

/*   PhyML_Printf("\n"); */
/*   For(i,6) PhyML_Printf("%d ",mod->rr_param_num[i]); */
/*   For(i,mod->n_diff_rr_param) PhyML_Printf("\n. Class %d size %d",i+1,mod->r_mat->n_rr_param_per_cat[i]); */
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

// Update rate across sites distribution settings.

void Update_RAS(t_mod *mod)
{
  phydbl sum;
  int i;



  if(mod->ras->free_mixt_rates == NO) DiscreteGamma(mod->ras->gamma_r_proba->v,
                            mod->ras->gamma_rr->v,
                            mod->ras->alpha->v,
                            mod->ras->alpha->v,
                            mod->ras->n_catg,
                            mod->ras->gamma_median);
  else
    {


      if(mod->ras->sort_rate_classes == YES)
        {
          Qksort(mod->ras->gamma_r_proba_unscaled->v,NULL,0,mod->ras->n_catg-1); // Unscaled class frequencies sorted in increasing order

          // Update class frequencies
          For(i,mod->ras->n_catg)
            {
              if(!i)
                mod->ras->gamma_r_proba->v[i] =
                  mod->ras->gamma_r_proba_unscaled->v[i] /  (mod->ras->gamma_r_proba_unscaled->v[mod->ras->n_catg-1]) ;
              else
                mod->ras->gamma_r_proba->v[i] =
                  (mod->ras->gamma_r_proba_unscaled->v[i] - mod->ras->gamma_r_proba_unscaled->v[i-1]) /
                  (mod->ras->gamma_r_proba_unscaled->v[mod->ras->n_catg-1]) ;
            }
        }
      else
        {
          sum = 0.0;
          For(i,mod->ras->n_catg) sum += mod->ras->gamma_r_proba_unscaled->v[i];
          For(i,mod->ras->n_catg) mod->ras->gamma_r_proba->v[i] = mod->ras->gamma_r_proba_unscaled->v[i] / sum;
        }

      do
        {
          sum = .0;
          For(i,mod->ras->n_catg)
            {
              if(mod->ras->gamma_r_proba->v[i] < 0.01) mod->ras->gamma_r_proba->v[i]=0.01;
              if(mod->ras->gamma_r_proba->v[i] > 0.99) mod->ras->gamma_r_proba->v[i]=0.99;
              sum += mod->ras->gamma_r_proba->v[i];
            }
          For(i,mod->ras->n_catg) mod->ras->gamma_r_proba->v[i]/=sum;
        }
      while((sum > 1.01) || (sum < 0.99));


      // Update class rates
      sum = .0;
      For(i,mod->ras->n_catg) sum += mod->ras->gamma_r_proba->v[i] * mod->ras->gamma_rr_unscaled->v[i];

      if(mod->ras->normalise_rr == YES)
        For(i,mod->ras->n_catg)
          mod->ras->gamma_rr->v[i] = mod->ras->gamma_rr_unscaled->v[i]/sum;
      else
        For(i,mod->ras->n_catg)
          mod->ras->gamma_rr->v[i] = mod->ras->gamma_rr_unscaled->v[i] * mod->ras->free_rate_mr->v;

      /* printf("\n"); */
      /* For(i,mod->ras->n_catg)  */
      /*   printf("\nx %3d %12f %12f xx %12f %12f", */
      /*          mod->ras->normalise_rr, */
      /*          mod->ras->gamma_r_proba->v[i], */
      /*          mod->ras->gamma_rr->v[i], */
      /*          mod->ras->gamma_r_proba_unscaled->v[i], */
      /*          mod->ras->gamma_rr_unscaled->v[i]); */
    }
#ifdef BEAGLE
  if(UNINITIALIZED != mod->b_inst)
      update_beagle_ras(mod);
#endif

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Update_Efrq(t_mod *mod)
{
  phydbl sum;
  int i;

  if((mod->io->datatype == NT) && (mod->s_opt->opt_state_freq))
    {
      sum = .0;
      For(i,mod->ns) sum += FABS(mod->e_frq->pi_unscaled->v[i]);
      For(i,mod->ns) mod->e_frq->pi->v[i] = FABS(mod->e_frq->pi_unscaled->v[i])/sum;

      do
        {
          sum = .0;
          For(i,mod->ns)
            {
              if(mod->e_frq->pi->v[i] < 0.01) mod->e_frq->pi->v[i]=0.01;
              if(mod->e_frq->pi->v[i] > 0.99) mod->e_frq->pi->v[i]=0.99;
              sum += mod->e_frq->pi->v[i];
            }
          For(i,mod->ns) mod->e_frq->pi->v[i]/=sum;
        }
      while((sum > 1.01) || (sum < 0.99));

#ifdef BEAGLE
      if(UNINITIALIZED != mod->b_inst)
        update_beagle_efrqs(mod);
#endif
    }



}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Set_Model_Parameters(t_mod *mod)
{
  Update_RAS(mod);
  Update_Efrq(mod);
  Update_Eigen(mod);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Update_Eigen(t_mod *mod)
{
  int result, n_iter;
  phydbl scalar;
  int i;


  if(mod->is_mixt_mod)
    {
      MIXT_Update_Eigen(mod);
      return;
    }

  if(mod->update_eigen)
    {
      //Update the Q-matrix first before computing the Eigen(because the Eigen is computed based on the Q-matrix)
      if(mod->use_m4mod == NO)
        {
          if(mod->io->datatype == NT)
            {
              if(mod->whichmodel == GTR)
                 Update_Qmat_GTR(mod->r_mat->rr->v, mod->r_mat->rr_val->v, mod->r_mat->rr_num->v, mod->e_frq->pi->v, mod->r_mat->qmat->v);
              else if(mod->whichmodel == CUSTOM)
                 Update_Qmat_GTR(mod->r_mat->rr->v, mod->r_mat->rr_val->v, mod->r_mat->rr_num->v, mod->e_frq->pi->v, mod->r_mat->qmat->v);
              else if(mod->whichmodel == HKY85)
                 Update_Qmat_HKY(mod->kappa->v, mod->e_frq->pi->v, mod->r_mat->qmat->v);
              else /* Any other nucleotide-based t_mod */
                 Update_Qmat_HKY(mod->kappa->v, mod->e_frq->pi->v, mod->r_mat->qmat->v);
            }
        }
      else
        {
          M4_Update_Qmat(mod->m4mod,mod);
        }

      //Now compute the Eigen
      scalar   = 1.0;
      n_iter   = 0;
      result   = 0;

      For(i,mod->ns*mod->ns) mod->r_mat->qmat_buff->v[i] = mod->r_mat->qmat->v[i];

      /* compute eigenvectors/values */
      /*       if(!EigenRealGeneral(mod->eigen->size,mod->r_mat->qmat,mod->eigen->e_val, */
      /* 			  mod->eigen->e_val_im, mod->eigen->r_e_vect, */
      /* 			  mod->eigen->space_int,mod->eigen->space)) */

      if(!Eigen(1,mod->r_mat->qmat_buff->v,mod->eigen->size,mod->eigen->e_val,
                mod->eigen->e_val_im,mod->eigen->r_e_vect,
                mod->eigen->r_e_vect_im,mod->eigen->space))
        {
          /* compute inverse(Vr) into Vi */
          For (i,mod->ns*mod->ns) mod->eigen->l_e_vect[i] = mod->eigen->r_e_vect[i];
          while(!Matinv(mod->eigen->l_e_vect, mod->eigen->size, mod->eigen->size,YES))
            {
	      PhyML_Printf("\n== Trying Q<-Q*scalar and then Root<-Root/scalar to fix this...\n");
              scalar += scalar / 3.;
              For(i,mod->eigen->size*mod->eigen->size) mod->r_mat->qmat_buff->v[i]  = mod->r_mat->qmat->v[i];
              For(i,mod->eigen->size*mod->eigen->size) mod->r_mat->qmat_buff->v[i] *= scalar;
              result = Eigen(1,mod->r_mat->qmat_buff->v,mod->eigen->size,mod->eigen->e_val,
                     mod->eigen->e_val_im,mod->eigen->r_e_vect,
                     mod->eigen->r_e_vect_im,mod->eigen->space);
              if (result == -1)
                {
                  PhyML_Printf("\n== Eigenvalues/vectors computation did not converge: computation cancelled."); 
                  Exit("\n");
                }
              else if (result == 1)
		{
                  PhyML_Printf("\n== Complex eigenvalues/vectors: computation cancelled.");
                  Exit("\n");
                }

              For (i,mod->eigen->size*mod->eigen->size) mod->eigen->l_e_vect[i] = mod->eigen->r_e_vect[i];
              n_iter++;
	      if(n_iter > 100) 
                {                  
                  PhyML_Printf("\n== Cannot work out eigen vectors.");
                  Exit("\n");
                }
            };
          For(i,mod->eigen->size) mod->eigen->e_val[i] /= scalar;

          /* compute the diagonal terms of EXP(D) */
          For(i,mod->ns) mod->eigen->e_val[i] = (phydbl)EXP(mod->eigen->e_val[i]);


      /* int j; */
      /* double *U,*V,*R; */
      /* double *expt; */
      /* double *uexpt; */
      /* int n; */

      /* expt  = mod->eigen->e_val_im; */
      /* uexpt = mod->eigen->r_e_vect_im; */
      /* U     = mod->eigen->r_e_vect; */
      /* V     = mod->eigen->l_e_vect; */
      /* R     = mod->eigen->e_val; /\* exponential of the eigen value matrix *\/ */
      /* n     = mod->ns; */

      /* PhyML_Printf("\n"); */
      /* PhyML_Printf("\n. Q\n"); */
      /* For(i,n) { For(j,n) PhyML_Printf("%7.3f ",mod->eigen->q[i*n+j]); PhyML_Printf("\n"); } */
      /* PhyML_Printf("\n. U\n"); */
      /* For(i,n) { For(j,n) PhyML_Printf("%7.3f ",U[i*n+j]); PhyML_Printf("\n"); } */
      /* PhyML_Printf("\n"); */
      /* PhyML_Printf("\n. V\n"); */
      /* For(i,n) { For(j,n) PhyML_Printf("%7.3f ",V[i*n+j]); PhyML_Printf("\n"); } */
      /* PhyML_Printf("\n"); */
      /* PhyML_Printf("\n. Eigen\n"); */
      /* For(i,n)  PhyML_Printf("%E ",mod->eigen->e_val[i]); */
      /* PhyML_Printf("\n"); */

/* 	  Exit("\n"); */
#ifdef BEAGLE
          //Recall that BEAGLE is initialized *after* all the model parameters are set
          //IOW, this function may be called before BEAGLE is initialized ("chicken-egg")
          if(UNINITIALIZED != mod->b_inst)
              update_beagle_eigen(mod);
#endif
        }
      else
        {
          PhyML_Printf("\n. Eigenvalues/vectors computation does not converge : computation cancelled");
          Warn_And_Exit("\n");
        }
    }

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Switch_From_M4mod_To_Mod(t_mod *mod)
{
  int i;

  mod->use_m4mod = 0;
  mod->ns = mod->m4mod->n_o;
  For(i,mod->ns) mod->e_frq->pi->v[i] = mod->m4mod->o_fq[i];
  mod->eigen->size = mod->ns;
  Switch_Eigen(YES,mod);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void PMat_MGF_Gamma(phydbl *Pij, phydbl shape, phydbl scale, phydbl scaling_fact, t_mod *mod)
{
  int dim;
  int i,j,k;
  phydbl *uexpt,*imbd;

  dim = mod->eigen->size;
  uexpt = mod->eigen->r_e_vect_im;
  imbd  = mod->eigen->e_val_im;

  /* Get the eigenvalues of Q (not the exponentials) */
  For(i,dim) imbd[i]  = LOG(mod->eigen->e_val[i]);

  /* Multiply them by the scaling factor */
  For(i,dim) imbd[i]  *= scaling_fact;

  For(i,dim) imbd[i] *= -scale;
  For(i,dim) imbd[i] += 1.0;
  For(i,dim) imbd[i]  = POW(imbd[i],-shape);

  For(i,dim) For(k,dim) uexpt[i*dim+k] = mod->eigen->r_e_vect[i*dim+k] * imbd[k];

  For(i,dim) For(k,dim) Pij[dim*i+k] = .0;

  For(i,dim)
    {
      For(j,dim)
	{
	  For(k,dim)
	    {
	      Pij[dim*i+j] += (uexpt[i*dim+k] * mod->eigen->l_e_vect[k*dim+j]);
	    }
	  if(Pij[dim*i+j] < SMALL_PIJ) Pij[dim*i+j] = SMALL_PIJ;
	}
    }

  /* printf("\n. shape = %G scale = %G %f",shape,scale,Pij[1]); */
  /* printf("\n. Pij: %f",Pij[1]); */

  /* printf("\n. Pmat"); */
  /* For(i,dim) */
  /*   { */
  /*     printf("\n"); */
  /*     For(j,dim) */
  /* 	{ */
  /* 	  printf("%12f ",Pij[i*dim+j]); */
  /* 	} */
  /*   } */
  /* Exit("\n"); */
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Switch_From_Mod_To_M4mod(t_mod *mod)
{
  int i;
  mod->use_m4mod = 1;
  mod->ns = mod->m4mod->n_o * mod->m4mod->n_h;
  For(i,mod->ns) mod->e_frq->pi->v[i] = mod->m4mod->o_fq[i%mod->m4mod->n_o] * mod->m4mod->h_fq[i/mod->m4mod->n_o];
  mod->eigen->size = mod->ns;
  Switch_Eigen(YES,mod);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl General_Dist(phydbl *F, t_mod *mod, eigen *eigen_struct)
{
  phydbl *pi,*mod_pi;
  int i,j,k;
  phydbl dist;
  phydbl sum;
  phydbl sum_ev;
  phydbl *F_phydbl;


  /* TO DO : call eigen decomposition function for all nt models */

  F_phydbl = (phydbl *)mCalloc(eigen_struct->size*eigen_struct->size,sizeof(phydbl));
  pi       = (phydbl *)mCalloc(eigen_struct->size,sizeof(phydbl));
  mod_pi   = (phydbl *)mCalloc(eigen_struct->size,sizeof(phydbl));

  For(i,mod->ns) mod_pi[i] = mod->e_frq->pi->v[i];

  sum = .0;
  For(i,eigen_struct->size)
    {
      For(j,eigen_struct->size)
        {
          pi[i] += (F[eigen_struct->size*i+j] + F[eigen_struct->size*j+i])/2.;
          sum += F[eigen_struct->size*i+j];
        }
    }

  Make_Symmetric(&F,eigen_struct->size);
  Divide_Mat_By_Vect(&F,mod->e_frq->pi->v,eigen_struct->size);

  /* Eigen decomposition of pi^{-1} x F */
  For(i,eigen_struct->size) For(j,eigen_struct->size) F_phydbl[eigen_struct->size*i+j] = F[eigen_struct->size*i+j];

  if(Eigen(1,F_phydbl,mod->eigen->size,mod->eigen->e_val,
       mod->eigen->e_val_im,mod->eigen->r_e_vect,
       mod->eigen->r_e_vect_im,mod->eigen->space))
    {
      For(i,mod->ns) mod->e_frq->pi->v[i] = mod_pi[i];
      Update_Qmat_GTR(mod->r_mat->rr->v, mod->r_mat->rr_val->v, mod->r_mat->rr_num->v, mod->e_frq->pi->v, mod->r_mat->qmat->v);
      Free(pi);
      Free(mod_pi);
      return -1.;
    }

  /* Get the left eigen vector of pi^{-1} x F */
  For(i,eigen_struct->size*eigen_struct->size) eigen_struct->l_e_vect[i] = eigen_struct->r_e_vect[i];
  if(!Matinv(eigen_struct->l_e_vect,eigen_struct->size,eigen_struct->size,YES)<0)
    {
      For(i,mod->ns) mod->e_frq->pi->v[i] = mod_pi[i];
      Update_Qmat_GTR(mod->r_mat->rr->v, mod->r_mat->rr_val->v, mod->r_mat->rr_num->v, mod->e_frq->pi->v, mod->r_mat->qmat->v);
      Free(pi);
      Free(mod_pi);
      return -1.;
    }

  /* LOG of eigen values */
  For(i,eigen_struct->size)
    {
/*       if(eigen_struct->e_val[i] < 0.0) eigen_struct->e_val[i] = 0.0001; */
      eigen_struct->e_val[i] = (phydbl)LOG(eigen_struct->e_val[i]);
     }

  /* Matrix multiplications LOG(pi^{-1} x F) */
  For(i,eigen_struct->size) For(j,eigen_struct->size)
    eigen_struct->r_e_vect[eigen_struct->size*i+j] =
    eigen_struct->r_e_vect[eigen_struct->size*i+j] *
    eigen_struct->e_val[j];
  For(i,eigen_struct->size) For(j,eigen_struct->size) F[eigen_struct->size*i+j] = .0;
  For(i,eigen_struct->size) For(j,eigen_struct->size) For(k,eigen_struct->size)
    F[eigen_struct->size*i+j] += eigen_struct->r_e_vect[eigen_struct->size*i+k] * eigen_struct->l_e_vect[eigen_struct->size*k+j];


  /* Trace */
  dist = .0;
  For(i,eigen_struct->size) dist+=F[eigen_struct->size*i+i];

  sum_ev = .0;
  For(i,mod->ns) sum_ev += mod->eigen->e_val[i];

/*   dist /= sum_ev; */
  dist /= -4.;


/*   For(i,mod->ns) mod->e_frq->pi->v[i] = mod_pi[i]; */
/*   Update_Qmat_GTR(mod); */
  Free(pi);
  Free(mod_pi);
  Free(F_phydbl);

  if(isnan(dist)) return -1.;
  return dist;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl GTR_Dist(phydbl *F, phydbl alpha, eigen *eigen_struct)
{
  phydbl *pi;
  int i,j,k;
  phydbl dist;
  phydbl sum;
  phydbl *F_phydbl;

  pi       = (phydbl *)mCalloc(eigen_struct->size,sizeof(phydbl));
  F_phydbl = (phydbl *)mCalloc(eigen_struct->size*eigen_struct->size,sizeof(phydbl));

/*   /\* Waddell and Steel's example *\/ */
/*   F[4*0+0] = 1415./4898.; F[4*0+1] = 8./4898.;    F[4*0+2] = 55./4898.;  F[4*0+3] = 2./4898.; */
/*   F[4*1+0] = 4./4898.;    F[4*1+1] = 1371./4898.; F[4*1+2] = 1./4898.;   F[4*1+3] = 144./4898.; */
/*   F[4*2+0] = 73./4898.;   F[4*2+1] = 0./4898.;    F[4*2+2] = 578./4898.; F[4*2+3] = 0./4898.; */
/*   F[4*3+0] = 3./4898.;    F[4*3+1] = 117./4898.;  F[4*3+2] = 1./4898.;   F[4*3+3] = 1126./4898.; */


  For(i,eigen_struct->size)
    {
      For(j,eigen_struct->size)
    {
      pi[i] += (F[eigen_struct->size*i+j] + F[eigen_struct->size*j+i])/2.;
      sum += F[eigen_struct->size*i+j];
    }
    }

/* /\*   Jukes and Cantor correction *\/ */
/*   sum = .0; */
/*   For(i,eigen_struct->size) sum += F[eigen_struct->size*i+i]; */
/*   sum = 1.-sum; */
/*   For(i,eigen_struct->size*eigen_struct->size) F[i] = sum/12.; */
/*   For(i,eigen_struct->size) F[eigen_struct->size*i+i] = (1.-sum)/4.; */
/*   For(i,eigen_struct->size) pi[i] = 1./(phydbl)eigen_struct->size; */


  Make_Symmetric(&F,eigen_struct->size);
  Divide_Mat_By_Vect(&F,pi,eigen_struct->size);


  /* Eigen decomposition of pi^{-1} x F */
  For(i,eigen_struct->size) For(j,eigen_struct->size) F_phydbl[eigen_struct->size*i+j] = F[eigen_struct->size*i+j];
  if(Eigen(1,F_phydbl,eigen_struct->size,eigen_struct->e_val,
       eigen_struct->e_val_im,eigen_struct->r_e_vect,
       eigen_struct->r_e_vect_im,eigen_struct->space))
    {
      Free(pi);
      return -1.;
    }

  /* Get the left eigen vector of pi^{-1} x F */
  For(i,eigen_struct->size*eigen_struct->size) eigen_struct->l_e_vect[i] = eigen_struct->r_e_vect[i];
  if(!Matinv(eigen_struct->l_e_vect,eigen_struct->size,eigen_struct->size,YES)<0) {Free(pi); return -1.;}

  /* Equation (3) + inverse of the moment generating function for the gamma distribution (see Waddell & Steel, 1997) */
  For(i,eigen_struct->size)
    {
      if(eigen_struct->e_val[i] < 0.0)
    {
      eigen_struct->e_val[i] = 0.0001;
    }
      if(alpha < .0)
    eigen_struct->e_val[i] = (phydbl)LOG(eigen_struct->e_val[i]);
      else
    eigen_struct->e_val[i] = alpha * (1. - (phydbl)POW(eigen_struct->e_val[i],-1./alpha));
     }

  /* Matrix multiplications pi x LOG(pi^{-1} x F) */
  For(i,eigen_struct->size) For(j,eigen_struct->size)
    eigen_struct->r_e_vect[eigen_struct->size*i+j] =
    eigen_struct->r_e_vect[eigen_struct->size*i+j] * eigen_struct->e_val[j];
  For(i,eigen_struct->size) For(j,eigen_struct->size) F[eigen_struct->size*i+j] = .0;
  For(i,eigen_struct->size) For(j,eigen_struct->size) For(k,eigen_struct->size)
    F[eigen_struct->size*i+j] += eigen_struct->r_e_vect[eigen_struct->size*i+k] * eigen_struct->l_e_vect[eigen_struct->size*k+j];
  For(i,eigen_struct->size) For(j,eigen_struct->size) F[eigen_struct->size*i+j] *= pi[i];

  /* Trace */
  dist = .0;
  For(i,eigen_struct->size) dist-=F[eigen_struct->size*i+i];

/*   PhyML_Printf("\nDIST = %f\n",dist); Exit("\n"); */

  Free(pi);
  Free(F_phydbl);

  if(isnan(dist)) return -1.;
  return dist;
}


