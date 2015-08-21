/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include "stats.h"


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* RANDOM VARIATES GENERATORS */
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


/*********************************************************************/
/* A C-function for TT800 : July 8th 1996 Version */
/* by M. Matsumoto, email: matumoto@math.keio.ac.jp */
/* tt800() generate one pseudorandom number with double precision */
/* which is uniformly distributed on [0,1]-interval */
/* for each call.  One may choose any initial 25 seeds */
/* except all zeros. */

/* See: ACM Transactions on Modelling and Computer Simulation, */
/* Vol. 4, No. 3, 1994, pages 254-266. */

phydbl tt800()
{
  int M=7;
  unsigned long y;
  static int k = 0;
  static unsigned long x[25]={ /* initial 25 seeds, change as you wish */
    0x95f24dab, 0x0b685215, 0xe76ccae7, 0xaf3ec239, 0x715fad23,
    0x24a590ad, 0x69e4b5ef, 0xbf456141, 0x96bc1b7b, 0xa7bdf825,
    0xc1de75b7, 0x8858a9c9, 0x2da87693, 0xb657f9dd, 0xffdc8a9f,
    0x8121da71, 0x8b823ecb, 0x885d05f5, 0x4e20cd47, 0x5a9ad5d9,
    0x512c0c03, 0xea857ccd, 0x4cc1d30f, 0x8891a8a1, 0xa6b7aadb
  };
  static unsigned long mag01[2]={ 
    0x0, 0x8ebfd028 /* this is magic vector `a', don't change */
  };
  if (k==25) { /* generate 25 words at one time */
    int kk;
    for (kk=0;kk<25-M;kk++) {
      x[kk] = x[kk+M] ^ (x[kk] >> 1) ^ mag01[x[kk] % 2];
    }
    for (; kk<25;kk++) {
      x[kk] = x[kk+(M-25)] ^ (x[kk] >> 1) ^ mag01[x[kk] % 2];
    }
    k=0;
  }
  y = x[k];
  y ^= (y << 7) & 0x2b5b2500; /* s and b, magic vectors */
  y ^= (y << 15) & 0xdb8b0000; /* t and c, magic vectors */
  y &= 0xffffffff; /* you may delete this line if word size = 32 */
  /* 
     the following line was added by Makoto Matsumoto in the 1996 version
     to improve lower bit's corellation.
     Delete this line to o use the code published in 1994.
  */
  y ^= (y >> 16); /* added to the 1994 version */
  k++;
  return((phydbl)y / (unsigned long) 0xffffffff);
}

/*********************************************************************/

phydbl Uni()
{
  phydbl r,mx;
  mx = (phydbl)RAND_MAX;
  r  = (phydbl)rand();
  r /= mx;
  /* r = tt800(); */
  return r;
}

/*********************************************************************/
// Return a uniform draw u s.t., min <= u <= max
int Rand_Int(int min, int max)
{
/*   phydbl u;   */
/*   u = (phydbl)rand(); */
/*   u /=  (RAND_MAX); */
/*   u *= (max - min + 1); */
/*   u += min; */
/*   return (int)FLOOR(u); */

  int u;
  /* if(max < min) Generic_Exit(__FILE__,__LINE__,__FUNCTION__); */
  u = rand();
  return (u%(max+1-min)+min);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////



/********************* random Gamma generator ************************
* Properties:
* (1) X = Gamma(alpha,lambda) = Gamma(alpha,1)/lambda
* (2) X1 = Gamma(alpha1,1), X2 = Gamma(alpha2,1) independent
*     then X = X1+X2 = Gamma(alpha1+alpha2,1)
* (3) alpha = k = integer then
*     X = Gamma(k,1) = Erlang(k,1) = -sum(LOG(Ui)) = -LOG(prod(Ui))
*     where U1,...Uk iid uniform(0,1)
*
* Decompose alpha = k+delta with k = [alpha], and 0<delta<1
* Apply (3) for Gamma(k,1)
* Apply Ahrens-Dieter algorithm for Gamma(delta,1)
*/
 
phydbl Ahrensdietergamma(phydbl alpha)
{
  phydbl x = 0.;

  if (alpha>0.) 
    {
      phydbl y = 0.;
      phydbl b = (alpha+EXP(1.))/EXP(1.);
      phydbl p = 1./alpha;
      int go = 0;
      while (go==0) 
	{
	  phydbl u = Uni();
	  phydbl w = Uni();
	  phydbl v = b*u;
	  if (v<=1.) 
	    {
	      x = POW(v,p);
	      y = EXP(-x);
	    }
	  else 
	    {
	      x = -LOG(p*(b-v));
	      y = POW(x,alpha-1.);
	    }
	  go = (w<y); // x is accepted when go=1
	}
    }
  return x;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl Rgamma(phydbl shape, phydbl scale)
{
  /* Code below is stolen from R sources. Thanks to the R team! */
  /* References:
     [1] Shape parameter a >= 1.  Algorithm GD in:

	  Ahrens, J.H. and Dieter, U. (1982).
	  Generating gamma variates by a modified
	  rejection technique.
	  Comm. ACM, 25, 47-54.


    [2] Shape parameter 0 < a < 1. Algorithm GS in:

	  Ahrens, J.H. and Dieter, U. (1974).
	  Computer methods for sampling from gamma, beta,
	  poisson and binomial distributions.
	  Computing, 12, 223-246.
  */


  double a = (double)shape;
  /* Constants : */
  const static double sqrt32 = 5.656854;
  const static double exp_m1 = 0.36787944117144232159;/* exp(-1) = 1/e */

  /* Coefficients q[k] - for q0 = sum(q[k]*a^(-k))
   * Coefficients a[k] - for q = q0+(t*t/2)*sum(a[k]*v^k)
   * Coefficients e[k] - for exp(q)-1 = sum(e[k]*q^k)
   */
  const static double q1 = 0.04166669;
  const static double q2 = 0.02083148;
  const static double q3 = 0.00801191;
  const static double q4 = 0.00144121;
  const static double q5 = -7.388e-5;
  const static double q6 = 2.4511e-4;
  const static double q7 = 2.424e-4;
  
  const static double a1 = 0.3333333;
  const static double a2 = -0.250003;
  const static double a3 = 0.2000062;
  const static double a4 = -0.1662921;
  const static double a5 = 0.1423657;
  const static double a6 = -0.1367177;
  const static double a7 = 0.1233795;
  
  /* State variables [FIXME for threading!] :*/
  static double aa = 0.;
  static double aaa = 0.;
  static double s, s2, d;    /* no. 1 (step 1) */
  static double q0, b, si, c;/* no. 2 (step 4) */
  
  double e, p, q, r, t, u, v, w, x, ret_val;
  
  if(a < 0.0 || scale <= 0.0) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  
  if (a < 1.) 
    { /* GS algorithm for parameters a < 1 */
      if(a == 0) return 0.;
      e = 1.0 + exp_m1 * a;
      for(;;) 
        {
          p = e * Uni();
          if (p >= 1.0) 
            {
              x = -LOG((e - p) / a);
              if (Rexp(1.) >= (1.0 - a) * LOG(x))
                break;
	    } 
          else 
            {
              x = EXP(LOG(p) / a);
              if (Rexp(1.) >= x)
                break;
	    }
	}
      return scale * x;
    }

    /* --- a >= 1 : GD algorithm --- */

    /* Step 1: Recalculations of s2, s, d if a has changed */
    if (a != aa) 
      {
	aa = a;
	s2 = a - 0.5;
	s = SQRT(s2);
	d = sqrt32 - s * 12.0;
      }
    /* Step 2: t = standard normal deviate,
               x = (s,1/2) -normal deviate. */

    /* immediate acceptance (i) */
    t = Rnorm(0.0,1.0);
    x = s + 0.5 * t;
    ret_val = x * x;
    if (t >= 0.0) return scale * ret_val;

    /* Step 3: u = 0,1 - uniform sample. squeeze acceptance (s) */
    u = Uni();
    if (d * u <= t * t * t) return scale * ret_val;

    /* Step 4: recalculations of q0, b, si, c if necessary */

    if (a != aaa) 
      {
	aaa = a;
	r = 1.0 / a;
	q0 = ((((((q7 * r + q6) * r + q5) * r + q4) * r + q3) * r
	       + q2) * r + q1) * r;
        
	/* Approximation depending on size of parameter a */
	/* The constants in the expressions for b, si and c */
	/* were established by numerical experiments */
        
	if (a <= 3.686) 
          {
            b = 0.463 + s + 0.178 * s2;
            si = 1.235;
            c = 0.195 / s - 0.079 + 0.16 * s;
          } 
        else if (a <= 13.022) 
          {
            b = 1.654 + 0.0076 * s2;
            si = 1.68 / s + 0.275;
	    c = 0.062 / s + 0.024;
          } 
        else 
          {
            b = 1.77;
            si = 0.75;
	    c = 0.1515 / s;
          }
      }

    /* Step 5: no quotient test if x not positive */
    if (x > 0.0) 
      {
	/* Step 6: calculation of v and quotient q */
	v = t / (s + s);
	if (FABS(v) <= 0.25)
          q = q0 + 0.5 * t * t * ((((((a7 * v + a6) * v + a5) * v + a4) * v
                                    + a3) * v + a2) * v + a1) * v;
	else
          q = q0 - s * t + 0.25 * t * t + (s2 + s2) * log(1.0 + v);
        
        
	/* Step 7: quotient acceptance (q) */
	if (LOG(1.0 - u) <= q)
          return scale * ret_val;
      }
    
    for(;;) 
      {
	/* Step 8: e = standard exponential deviate
	 *	u =  0,1 -uniform deviate
	 *	t = (b,si)-double exponential (laplace) sample */
	e = Rexp(1.0);
	u = Uni();
	u = u + u - 1.0;
	if (u < 0.0)
	    t = b - si * e;
	else
	    t = b + si * e;
	/* Step	 9:  rejection if t < tau(1) = -0.71874483771719 */
	if (t >= -0.71874483771719) {
	    /* Step 10:	 calculation of v and quotient q */
	    v = t / (s + s);
	    if (FABS(v) <= 0.25)
		q = q0 + 0.5 * t * t *
		    ((((((a7 * v + a6) * v + a5) * v + a4) * v + a3) * v
		      + a2) * v + a1) * v;
	    else
		q = q0 - s * t + 0.25 * t * t + (s2 + s2) * LOG(1.0 + v);
	    /* Step 11:	 hat acceptance (h) */
	    /* (if q not positive go to step 8) */
	    if (q > 0.0) {
		w = EXP(q)-1.0;
		/* if t is rejected sample again at step 8 */
		if (c * FABS(u) <= w * EXP(e - 0.5 * t * t))
		    break;
	    }
	}
    } /* repeat .. until  `t' is accepted */
    x = s + 0.5 * t;
    return scale * x * x;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl Rexp(phydbl lambda)
{
  return -LOG(Uni()+1.E-30)/lambda;
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

/* Returns a random deviates from an exponential distribution */
/* left-truncated at 'left' and right-truncated at 'rght' */
phydbl Rexp_Trunc(phydbl lambda, phydbl left, phydbl rght)
{
  phydbl u;
  u = Uni();
  return (left-LOG(1. - u*(1.-EXP(-lambda*(rght-left))))/lambda);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl Rnorm(phydbl mean, phydbl sd)
{
  /* Box-Muller transformation */
  phydbl u1, u2, res;
  
  /* u1=Uni(); */
  /* u2=Uni(); */
  /* u1 = SQRT(-2.*LOG(u1))*COS(6.28318530717959f*u2); */

  /* Polar */
  phydbl d,x,y;

  do
    {
      u1=Uni();
      u2=Uni();
      x = 2.*u1-1.;
      y = 2.*u2-1.;
      d = x*x + y*y;
      if(d>.0 && d<1.) break;
    }
  while(1);
  u1 = x*SQRT((-2.*LOG(d))/d);

  res = u1*sd+mean;

  if(isnan(res) || isinf(res))
    {
      printf("\n. res=%f sd=%f mean=%f u1=%f u2=%f",res,sd,mean,u1,u2);
    }
  return res;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl *Rnorm_Multid(phydbl *mu, phydbl *cov, int dim)
{
  phydbl *L,*x,*y;
  int i,j;
  
  x = (phydbl *)mCalloc(dim,sizeof(phydbl));
  y = (phydbl *)mCalloc(dim,sizeof(phydbl));

  L = (phydbl *)Cholesky_Decomp(cov,dim);

  For(i,dim) x[i]=Rnorm(0.0,1.0);
  For(i,dim) For(j,dim) y[i] += L[i*dim+j]*x[j];
  For(i,dim) y[i] += mu[i];

  Free(L);
  Free(x);

  return(y);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Rnorm_Trunc_Inverse(phydbl mean, phydbl sd, phydbl min, phydbl max, int *error)
{

  phydbl u, ret_val,eps;
  phydbl z;
  phydbl z_min,z_max;
  phydbl cdf_min, cdf_max;

  z      = 0.0;
  u      = -1.0;
  *error = 0;
  
  if(sd < 1.E-100)
    {
      PhyML_Printf("\n. Small variance detected in Rnorm_Trunc.");
      PhyML_Printf("\n. mean=%f sd=%f min=%f max=%f",mean,sd,min,max);
      *error = 1;
      return -1.0;
    }

  z_min = (min - mean)/sd;
  z_max = (max - mean)/sd;

  eps = (z_max-z_min)/1E+6;


  /*       Simple inversion method. Seems to work well. Needs more thorough testing though... */
  cdf_min = Pnorm(z_min,0.0,1.0);
  cdf_max = Pnorm(z_max,0.0,1.0);
  u = cdf_min + (cdf_max-cdf_min) * Uni();
  z = PointNormal(u);
	
  if((z < z_min-eps) || (z > z_max+eps))
    {
      *error = 1;
      PhyML_Printf("\n. Numerical precision issue detected in Rnorm_Trunc.");
      PhyML_Printf("\n. z = %f",z);
      PhyML_Printf("\n. mean=%f sd=%f z_min=%f z_max=%f min=%f max=%f",mean,sd,z_min,z_max,min,max);
      ret_val = (max - min)/2.;
      Exit("\n");
    }
  
  ret_val = z*sd+mean;
  
  return ret_val;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl Rnorm_Trunc(phydbl mean, phydbl sd, phydbl min, phydbl max, int *error)
{

  phydbl ret_val;
  phydbl z, q, u, a;
  phydbl z_min,z_max;
  int n_max_iter,n_iter;
  int algo;

  z          = 0.0;
  *error     = NO;
  ret_val    = INFINITY;
  n_max_iter = 100000;

  if(sd < 1.E-100)
    {
      PhyML_Printf("\n. Small variance detected in Rnorm_Trunc.");
      PhyML_Printf("\n. mean=%f sd=%f min=%f max=%f",mean,sd,min,max);
      *error = YES;
      return -1.0;
    }

  if(max < min)
    {
      PhyML_Printf("\n. Max < Min");
      PhyML_Printf("\n. mean=%f sd=%f min=%f max=%f",mean,sd,min,max);
      *error = YES;
      return -1.0;
    }

  z_min = (min - mean)/sd;
  z_max = (max - mean)/sd;

  if(z_min < 0.0 && z_max > 0 && (z_max - z_min > SQRT(2*PI)))
    {
      algo = 0;
    }
  if((z_min > 0.0 || Are_Equal(z_min,0.0,1.E-10)) && z_max > z_min + 2.*SQRT(EXP(1.0))/(z_min + SQRT(z_min*z_min+4.)) * EXP((z_min*z_min - z_min*SQRT(z_min*z_min+4.))/4.))
    {
      algo = 1;
    }
  else if((z_max < 0.0 || Are_Equal(z_max,0.0,1.E-10)) && -z_min > -z_max + 2.*SQRT(EXP(1.0))/(-z_max + SQRT(z_max*z_max+4.)) * EXP((z_max*z_max - (-z_max)*SQRT(z_max*z_max+4.))/4.))
    {
      algo = 2;
    }
  else
    {
      algo = 3;
    }


  switch(algo)
    {
    case 0:
      {
        n_iter = 0;
        do 
          { 
            z = Rnorm(0.0,1.0); 
            n_iter++;
            if(n_iter > n_max_iter)
              {
                PhyML_Printf("\n== Too many iterations in Rnorm_Trunc()");
                *error = YES; 
              }
          }
        while(z < z_min || z > z_max);
        break;
      }
    case 1:
      {
        n_iter = 0;
        do
          {
            a = (z_min + SQRT(z_min*z_min+4.))/2.;
            q = Rexp(a) + z_min;
            u = Uni();
            n_iter++;
            if(n_iter > n_max_iter)
              {
                PhyML_Printf("\n== Too many iterations in Rnorm_Trunc()");
                *error = YES; 
              }
          }while(u > EXP(-POW(q-a,2)/2.));
        z = q;
        break;
      }
    case 2:
      {
        n_iter = 0;
        do
          {
            a = (-z_max + SQRT(z_max*z_max+4.))/2.;
            q = Rexp(a) - z_max;
            u = Uni();
            n_iter++;
            if(n_iter > n_max_iter)
              {
                PhyML_Printf("\n== Too many iterations in Rnorm_Trunc()");
                *error = YES; 
              }
          }while(u > EXP(-POW(q-a,2)/2.));
        z = -q;
        break;
      }
    case 3:
      {
        n_iter = 0;
        do
          {
            z = Uni()*(z_max - z_min) + z_min;
            if(z_min < 0.0 && z_max > 0.0) q = EXP(-z*z/2.);
            else if (z_max < 0.0) q = EXP((z_max*z_max-z*z)/2.);
            else q = EXP((z_min*z_min-z*z)/2.);
            u = Uni();
            n_iter++;
            if(n_iter > n_max_iter)
              {
                PhyML_Printf("\n== Too many iterations in Rnorm_Trunc()");
                *error = YES; 
              }
          }while(u > q); 
        break;
      }

    default: Generic_Exit(__FILE__,__LINE__,__FUNCTION__);

    }

  ret_val = z*sd+mean;

  return ret_val;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl *Rnorm_Multid_Trunc(phydbl *mean, phydbl *cov, phydbl *min, phydbl *max, int dim)
{
  int i,j;
  phydbl *L,*x, *u;
  phydbl up, low, rec;
  int err;
  
  u = (phydbl *)mCalloc(dim,sizeof(phydbl)); 
  x = (phydbl *)mCalloc(dim,sizeof(phydbl));
 
  L = Cholesky_Decomp(cov,dim);
  
  low = (min[0]-mean[0])/L[0*dim+0];
  up  = (max[0]-mean[0])/L[0*dim+0];
  u[0] = Rnorm_Trunc(0.0,1.0,low,up,&err);

  for(i=1;i<dim;i++)
    {
      rec = .0;
      For(j,i) rec += L[i*dim+j] * u[j];
      low  = (min[i]-mean[i]-rec)/L[i*dim+i];
      up   = (max[i]-mean[i]-rec)/L[i*dim+i];
      u[i] = Rnorm_Trunc(0.0,1.0,low,up,&err);
    }

  x = Matrix_Mult(L,u,dim,dim,dim,1);

/*   PhyML_Printf("\n>>>\n"); */
/*   For(i,dim) */
/*     { */
/*       For(j,dim) */
/* 	{ */
/* 	  PhyML_Printf("%10lf ",L[i*dim+j]); */
/* 	} */
/*       PhyML_Printf("\n"); */
/*     } */
/*   PhyML_Printf("\n"); */

/*   For(i,dim) PhyML_Printf("%f ",u[i]); */
/*   PhyML_Printf("\n"); */

  
/*   PhyML_Printf("\n"); */
/*   For(i,dim) PhyML_Printf("%10lf ",x[i]); */
/*   PhyML_Printf("\n<<<\n"); */

  For(i,dim) x[i] += mean[i];

  Free(L);
  Free(u);
  
  return x;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
/* Inversion method for sampling from Geometric distibution */
phydbl Rgeom(phydbl p)
{
  phydbl x,u;

  u = Uni();
  if(u < SMALL) return(0.0);
  x = LOG(u) / LOG(1. - p);
  
  return(CEIL(x));  
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl Dgeom(phydbl k, phydbl p, int logit)
{
  phydbl prob;

  if(k < 1.) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  if(p > 1.-SMALL)
    {
      if(logit == YES) return(-INFINITY);
      else return(0.0);
    }

  if(logit == YES)
    prob = (k - 1.)*LOG(1. - p) + LOG(p);
  else
    prob = POW(1.-p,k-1.)*p;

  return(prob);
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl Pgeom(phydbl k, phydbl p)
{

  if(k < 1.) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  if(p > 1.-SMALL) return(0.0);

  return(1. - POW((1. - p),k));

}


/*
 * Random variates from the Poisson distribution. Completely stolen from R code.
 *
 * REFERENCE
 *
 * Ahrens, J.H. and Dieter, U. (1982).
 * Computer generation of Poisson deviates
 * from modified normal distributions.
 * ACM Trans. Math. Software 8, 163-179.
 */
phydbl Rpois(phydbl mmu)
{
  double mu = (double)mmu;
  
  double a0	= -0.5;
  double a1	= 0.3333333;
  double a2	= -0.2500068;
  double a3	= 0.2000118;
  double a4	= -0.1661269;
  double a5	= 0.1421878;
  double a6	= -0.1384794;
  double a7	= 0.1250060;
  double one_7	= 0.1428571428571428571;
  double one_12	= 0.0833333333333333333;
  double one_24	= 0.0416666666666666667;

  /* Factorial Table (0:9)! */
  const static double fact[10] =
    {
      1., 1., 2., 6., 24., 120., 720., 5040., 40320., 362880.
    };
  /* These are static --- persistent between calls for same mu : */
  static int l, m;
  static double b1, b2, c, c0, c1, c2, c3;
  static double pp[36], p0, p, q, s, d, omega;
  static double big_l;/* integer "w/o overflow" */
  static double muprev = 0., muprev2 = 0.;/*, muold = 0.*/
  /* Local Vars [initialize some for -Wall]: */
  double del, difmuk= 0., E= 0., fk= 0., fx, fy, g, px, py, t, u= 0., v, x;
  double pois = -1.;
  int k, kflag, big_mu, new_big_mu = FALSE;
  
  if (isnan(mu) || mu < 0.0) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  if (mu <= 0.) return 0.;

  big_mu = mu >= 10.;
  if(big_mu) new_big_mu = FALSE;
  if (!(big_mu && mu == muprev)) {/* maybe compute new persistent par.s */
    if (big_mu) 
      {
        new_big_mu = TRUE;
        /* Case A. (recalculation of s,d,l because mu has changed):
         * The poisson probabilities pk exceed the discrete normal
         * probabilities fk whenever k >= m(mu).
         */
        muprev = mu;
        s = sqrt(mu);
        d = 6. * mu * mu;
        big_l = floor(mu - 1.1484);
        /* = an upper bound to m(mu) for all mu >= 10.*/
      }
    else 
      { 
        /* Small mu ( < 10) -- not using normal approx. */
        /* Case B. (start new table and calculate p0 if necessary) */
        /*muprev = 0.;-* such that next time, mu != muprev ..*/
        if (mu != muprev) 
          {
            muprev = mu;
            m = MAX(1, (int) mu);
            l = 0; /* pp[] is already ok up to pp[l] */
            q = p0 = p = exp(-mu);
          }
        for(;;) 
          {
            /* Step U. uniform sample for inversion method */
            u = Uni();
            if (u <= p0) return 0.;
            /* Step T. table comparison until the end pp[l] of the
               pp-table of cumulative poisson probabilities
               (0.458 > ~= pp[9](= 0.45792971447) for mu=10 ) */
            if (l != 0) 
              {
                for (k = (u <= 0.458) ? 1 : MIN(l, m); k <= l; k++) if (u <= pp[k]) return((phydbl)k);
                if (l == 35) /* u > pp[35] */
                  continue;
              }
            /* Step C. creation of new poisson
               probabilities p[l..] and their cumulatives q =: pp[k] */
            l++;
            for (k = l; k <= 35; k++) 
              {
                p *= mu / k;
                q += p;
                pp[k] = q;
                if (u <= q) 
                  {
                    l = k;
                    return((phydbl)k);
                  }
              }
            l = 35;
          } /* end(repeat) */
      }/* mu < 10 */
  } /* end {initialize persistent vars} */
  /* Only if mu >= 10 : ----------------------- */
  /* Step N. normal sample */
  g = mu + s * Rnorm(0.0,1.0);/* norm_rand() ~ N(0,1), standard normal */
  if (g >= 0.) 
    {
      pois = floor(g);
      /* Step I. immediate acceptance if pois is large enough */
      if (pois >= big_l)
        return((phydbl)pois);
      /* Step S. squeeze acceptance */
      fk = pois;
      difmuk = mu - fk;
      u = Uni(); /* ~ U(0,1) - sample */
      if (d * u >= difmuk * difmuk * difmuk)
        return((phydbl)pois);
    }
  /* Step P. preparations for steps Q and H.
     (recalculations of parameters if necessary) */
  if (new_big_mu || mu != muprev2) {
    /* Careful! muprev2 is not always == muprev
       because one might have exited in step I or S
    */
    muprev2 = mu;
    omega = M_1_SQRT_2PI / s;
    /* The quantities b1, b2, c3, c2, c1, c0 are for the Hermite
     * approximations to the discrete normal probabilities fk. */
    b1 = one_24 / mu;
    b2 = 0.3 * b1 * b1;
    c3 = one_7 * b1 * b2;
    c2 = b2 - 15. * c3;
    c1 = b1 - 6. * b2 + 45. * c3;
    c0 = 1. - b1 + 3. * b2 - 15. * c3;
    c = 0.1069 / mu; /* guarantees majorization by the 'hat'-function. */
  }
  if (g >= 0.) {
    /* 'Subroutine' F is called (kflag=0 for correct return) */
    kflag = 0;
    goto Step_F;
  }
  for(;;) 
    {
      /* Step E. Exponential Sample */
      E = Rexp(1.0); /* ~ Exp(1) (standard exponential) */
      /* sample t from the laplace 'hat'
         (if t <= -0.6744 then pk < fk for all mu >= 10.) */
      u = 2. * Uni() - 1.;
      /* t = 1.8 + fsign(E, u); */
      t = 1.8 + ((u >= 0.0) ? fabs(E) : -fabs(E));
      if (t > -0.6744)
        {
          pois = floor(mu + s * t);
          fk = pois;
          difmuk = mu - fk;
          /* 'subroutine' F is called (kflag=1 for correct return) */
          kflag = 1;
        Step_F: /* 'subroutine' F : calculation of px,py,fx,fy. */
          if (pois < 10) { /* use factorials from table fact[] */
            px = -mu;
            py = pow(mu, pois) / fact[(int)pois];
          }
          else {
            /* Case pois >= 10 uses polynomial approximation
               a0-a7 for accuracy when advisable */
            del = one_12 / fk;
            del = del * (1. - 4.8 * del * del);
            v = difmuk / fk;
            if (fabs(v) <= 0.25)
              px = fk * v * v * (((((((a7 * v + a6) * v + a5) * v + a4) *
                                    v + a3) * v + a2) * v + a1) * v + a0)
                - del;
            else /* |v| > 1/4 */
              px = fk * log(1. + v) - difmuk - del;
            py = M_1_SQRT_2PI / sqrt(fk);
          }
          x = (0.5 - difmuk) / s;
          x *= x;/* x^2 */
          fx = -0.5 * x;
          fy = omega * (((c3 * x + c2) * x + c1) * x + c0);
          if (kflag > 0) {
            /* Step H. Hat acceptance (E is repeated on rejection) */
            if (c * fabs(u) <= py * exp(px + E) - fy * exp(fx + E))
              break;
          } else
            /* Step Q. Quotient acceptance (rare case) */
            if (fy - u * fy <= py * exp(px - fx))
              break;
        }/* t > -.67.. */
    }
  return((phydbl)pois);
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* DENSITIES / PROBA */
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Dnorm_Moments(phydbl x, phydbl mean, phydbl var)
{
  phydbl dens,sd,pi;

  pi = 3.141593;
  sd = SQRT(var);

  dens = 1./(SQRT(2*pi)*sd)*EXP(-((x-mean)*(x-mean)/(2.*sd*sd)));

  return dens;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Dnorm(phydbl x, phydbl mean, phydbl sd)
{
  phydbl dens;

  /* dens = -(.5*LOG2PI+LOG(sd))  - .5*POW(x-mean,2)/POW(sd,2); */
  /* return EXP(dens); */
  
  x = (x-mean)/sd;

  dens = M_1_SQRT_2_PI * EXP(-0.5*x*x);
  
  return dens / sd;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl Log_Dnorm(phydbl x, phydbl mean, phydbl sd, int *err)
{
  phydbl dens;

  *err = NO;

  x = (x-mean)/sd;
  
  dens = -(phydbl)LOG(SQRT(2.*PI)) - x*x*0.5 - LOG(sd);

  if(dens < -BIG)
    {
      PhyML_Printf("\n. dens=%f -- x=%f mean=%f sd=%f\n",dens,x,mean,sd);
      *err = 1;
    }

  return dens;
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl Log_Dnorm_Trunc(phydbl x, phydbl mean, phydbl sd, phydbl lo, phydbl up, int *err)
{
  phydbl log_dens;
  phydbl cdf_up, cdf_lo;

  if(x < lo || x > up) return -230.;

  *err = NO;
  cdf_lo = cdf_up = 0.0;

  log_dens = Log_Dnorm(x,mean,sd,err);

  if(*err == YES)
    {
      PhyML_Printf("\n== mean=%f sd=%f lo=%f up=%f cdf_lo=%G CDF_up=%G log_dens=%G",mean,sd,lo,up,cdf_lo,cdf_up,log_dens);
      PhyML_Printf("\n== Warning in file %s at line %d\n",__FILE__,__LINE__);
      *err = YES;
    }

  cdf_up = Pnorm(up,mean,sd);
  cdf_lo = Pnorm(lo,mean,sd);

  if(cdf_up - cdf_lo < 1.E-20)
    {
      log_dens = -230.; /* ~LOG(1.E-100) */
    }
  else
    {
      log_dens -= LOG(cdf_up - cdf_lo);
    }

  if(isnan(log_dens) || isinf(FABS(log_dens)))
    {
      PhyML_Printf("\n. x=%f mean=%f sd=%f lo=%f up=%f cdf_lo=%G CDF_up=%G log_dens=%G",x,mean,sd,lo,up,cdf_lo,cdf_up,log_dens);
      PhyML_Printf("\n. Warning in file %s at line %d\n",__FILE__,__LINE__);
      *err = YES;
    }

  return log_dens;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Dnorm_Trunc(phydbl x, phydbl mean, phydbl sd, phydbl lo, phydbl up)
{
  phydbl dens;
  phydbl cdf_up, cdf_lo;

  dens   = Dnorm(x,mean,sd);
  cdf_up = Pnorm(up,mean,sd);
  cdf_lo = Pnorm(lo,mean,sd);

  dens /= (cdf_up - cdf_lo);

  if(isnan(dens) || isinf(FABS(dens)))
    {
      PhyML_Printf("\n== mean=%f sd=%f lo=%f up=%f cdf_lo=%G CDF_up=%G",mean,sd,lo,up,cdf_lo,cdf_up);
      PhyML_Printf("\n== Err. in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");
    }

  return dens;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Dnorm_Multi(phydbl *x, phydbl *mu, phydbl *cov, int size, int _log)
{
  phydbl *xmmu,*invcov;
  phydbl *buff1,*buff2;
  int i;
  phydbl det,density;

  xmmu   = (phydbl *)mCalloc(size,sizeof(phydbl));
  invcov = (phydbl *)mCalloc(size*size,sizeof(phydbl));

  For(i,size) xmmu[i] = x[i] - mu[i];
  For(i,size*size) invcov[i] = cov[i];
  
  if(!Matinv(invcov,size,size,NO))
    {
      PhyML_Printf("\n== Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");      
    }

  buff1 = Matrix_Mult(xmmu,invcov,1,size,size,size);
  buff2 = Matrix_Mult(buff1,xmmu,1,size,size,1);
  
  det = Matrix_Det(cov,size,NO);
  /* det_1D(cov,size,&det); */

  density = size * LOG2PI + LOG(det) + buff2[0];
  density /= -2.;

/*   density = (1./(POW(2.*PI,size/2.)*SQRT(FABS(det)))) * EXP(-0.5*buff2[0]); */

  Free(xmmu);
  Free(invcov);
  Free(buff1);
  Free(buff2);

  return (_log)?(density):(EXP(density));
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Dnorm_Multi_Given_InvCov_Det(phydbl *x, phydbl *mu, phydbl *invcov, phydbl log_det, int size, int _log)
{
  phydbl *xmmu;
  phydbl *buff1,*buff2;
  int i;
  phydbl density;

  xmmu = (phydbl *)mCalloc(size,sizeof(phydbl));

  For(i,size) xmmu[i] = x[i] - mu[i];
  
  buff1 = Matrix_Mult(xmmu,invcov,1,size,size,size);
  buff2 = Matrix_Mult(buff1,xmmu,1,size,size,1);

  density = size * LOG2PI + log_det + buff2[0];
  density /= -2.;

  Free(xmmu);
  Free(buff1);
  Free(buff2);
  
  return (_log)?(density):(EXP(density));
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Pbinom(int N, int ni, phydbl p)
{
  return Bico(N,ni)*POW(p,ni)*POW(1-p,N-ni);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Bivariate_Normal_Density(phydbl x, phydbl y, phydbl mux, phydbl muy, phydbl sdx, phydbl sdy, phydbl rho)
{
  phydbl cx, cy;
  phydbl pi;
  phydbl dens;
  phydbl rho2;

  pi = 3.141593;

  cx = x - mux;
  cy = y - muy;

  rho2 = rho*rho;

  dens = 1./(2*pi*sdx*sdy*SQRT(1.-rho2));
  dens *= EXP((-1./(2.*(1.-rho2)))*(cx*cx/(sdx*sdx)+cy*cy/(sdy*sdy)+2*rho*cx*cy/(sdx*sdy)));
	      
  return dens;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// E(X) = n(1-p)/p ; V(X) = n(1-p)/p^2
phydbl Dnbinom(phydbl x, phydbl n, phydbl p, int logit)
{

  phydbl lnDens = LnGamma(x+n) - LnGamma(n) - LnFact(x) + n*LOG(p) + x*LOG(1.-p);

  if(logit == TRUE)
    {
      return(lnDens);
    }
  else
    {
      return(EXP(lnDens));
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl Rnbinom(phydbl n, phydbl p)
{
  phydbl y,x;
  y = Rgamma(n,(1.-p)/p);
  x = Rpois(y);
  return(x);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl Dgamma_Moments(phydbl x, phydbl mean, phydbl var)
{
  phydbl shape, scale;

  if(var < 1.E-20) 
    {
/*       var  = 1.E-20;  */
      PhyML_Printf("\n. var=%f mean=%f",var,mean);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");      
    }

  if(mean < 1.E-20) 
    { 
/*       mean = 1.E-20;  */
      PhyML_Printf("\n. var=%f mean=%f",var,mean);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");
    }


  shape = mean * mean / var;
  scale = var / mean;
  
  return(Dgamma(x,shape,scale));
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Dgamma(phydbl x, phydbl shape, phydbl scale)
{
  phydbl v;

  if(x > INFINITY) 
    {
      PhyML_Printf("\n. WARNING: huge value of x -> x = %G",x);
      x = 1.E+10;
    }

  if(x < 1.E-20)
    {
      if(x < 0.0) return 0.0;
      else
	{
	  PhyML_Printf("\n. WARNING: small value of x -> x = %G",x);
	  x = 1.E-20;
	}
    }


  if(scale < 0.0 || shape < 0.0)
    {
      PhyML_Printf("\n. scale=%f shape=%f",scale,shape);
      Exit("\n");
    }


  v = (shape-1.) * LOG(x) - shape * LOG(scale) - x / scale - LnGamma(shape);


  if(v < 500.)
    {
      v = EXP(v);
    }
  else
    {
      PhyML_Printf("\n. WARNING v=%f x=%f shape=%f scale=%f",v,x,shape,scale);
      PhyML_Printf("\n. LOG(x) = %G LnGamma(shape)=%G",LOG(x),LnGamma(shape));
      v = EXP(v);
      /* PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__); */
      /* Exit("\n"); */
    }

	 
  return v;
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl Dexp(phydbl x, phydbl param)
{
  return param * EXP(-param * x);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

/* Returns the density of an exponential distribution left-truncated */
/* at 'left' and right-truncated at 'rght' */

phydbl Dexp_Trunc(phydbl x, phydbl lambda, phydbl left, phydbl rght)
{
  return (lambda * EXP(-lambda * x))/(EXP(-lambda * left) - EXP(-lambda * rght));
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
/* Poisson probability */
phydbl Dpois(phydbl x, phydbl param, int logit)
{
  phydbl v;

  if(x < .0) 
    {
      if(logit == YES) return(-INFINITY);
      else return 0.0;
   }

  v = x * LOG(param) - param - Factln(x);
  if(logit == YES) return v;
  else
    {
      if(v < 500.)
        {
          v = EXP(v);
        }
      else
        {
          PhyML_Printf("\n. WARNING v=%f x=%f param=%f",v,x,param);
          v = EXP(500);
        }
      
      /*   PhyML_Printf("\n. Poi %f %f (x=%f param=%f)", */
      /* 	 v, */
      /* 	 POW(param,x) * EXP(-param) / EXP(LnGamma(x+1)), */
      /* 	 x,param); */
      /*   return POW(param,x) * EXP(-param) / EXP(LnGamma(x+1)); */
  
      return v;
    }
  return(-1.0);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* CDFs */
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Pnorm(phydbl x, phydbl mean, phydbl sd)
{
/*   const phydbl b1 =  0.319381530; */
/*   const phydbl b2 = -0.356563782; */
/*   const phydbl b3 =  1.781477937; */
/*   const phydbl b4 = -1.821255978; */
/*   const phydbl b5 =  1.330274429; */
/*   const phydbl p  =  0.2316419; */
/*   const phydbl c  =  0.39894228; */
  
  x = (x-mean)/sd;
  
/*   if(x >= 0.0) */
/*     { */
/*       phydbl t = 1.0 / ( 1.0 + p * x ); */
/*       return (1.0 - c * EXP( -x * x / 2.0 ) * t * */
/* 	      ( t *( t * ( t * ( t * b5 + b4 ) + b3 ) + b2 ) + b1 )); */
/*     } */
/*   else */
/*     { */
/*       phydbl t = 1.0 / ( 1.0 - p * x ); */
/*       return ( c * EXP( -x * x / 2.0 ) * t * */
/* 	       ( t *( t * ( t * ( t * b5 + b4 ) + b3 ) + b2 ) + b1 )); */
/*     } */

/* i_tail in {0,1,2} means: "lower", "upper", or "both" :
   if(lower) return  *cum := P[X <= x]
   if(upper) return *ccum := P[X >  x] = 1 - P[X <= x]
*/

/*   return Pnorm_Marsaglia(x); */
  return Pnorm_Ihaka_Derived_From_Cody(x);
}


/* G. Marsaglia. "Evaluating the Normal distribution". Journal of Statistical Software. 2004. Vol. 11. Issue 4. */
phydbl  Pnorm_Marsaglia(phydbl x)
{
  long double s=x,t=0,b=x,q=x*x,i=1;
  while(s!=t) s=(t=s)+(b*=q/(i+=2)); 
  return .5+s*exp(-.5*q-.91893853320467274178L);

}



/* Stolen from R source code */
#define SIXTEN 16

phydbl Pnorm_Ihaka_Derived_From_Cody(phydbl x)
{

    const static double a[5] = {
	2.2352520354606839287,
	161.02823106855587881,
	1067.6894854603709582,
	18154.981253343561249,
	0.065682337918207449113
    };
    const static double b[4] = {
	47.20258190468824187,
	976.09855173777669322,
	10260.932208618978205,
	45507.789335026729956
    };
    const static double c[9] = {
	0.39894151208813466764,
	8.8831497943883759412,
	93.506656132177855979,
	597.27027639480026226,
	2494.5375852903726711,
	6848.1904505362823326,
	11602.651437647350124,
	9842.7148383839780218,
	1.0765576773720192317e-8
    };
    const static double d[8] = {
	22.266688044328115691,
	235.38790178262499861,
	1519.377599407554805,
	6485.558298266760755,
	18615.571640885098091,
	34900.952721145977266,
	38912.003286093271411,
	19685.429676859990727
    };
    const static double p[6] = {
	0.21589853405795699,
	0.1274011611602473639,
	0.022235277870649807,
	0.001421619193227893466,
	2.9112874951168792e-5,
	0.02307344176494017303
    };
    const static double q[5] = {
	1.28426009614491121,
	0.468238212480865118,
	0.0659881378689285515,
	0.00378239633202758244,
	7.29751555083966205e-5
    };

    double xden, xnum, temp, del, eps, xsq, y;
    int i, lower, upper;
    double cum,ccum;
    int i_tail;
    
    i_tail = 0;
    cum = ccum = 0.0;

    if(isnan(x)) { cum = ccum = x; return (phydbl)cum; }

    /* Consider changing these : */
    eps = DBL_EPSILON * 0.5;

    /* i_tail in {0,1,2} =^= {lower, upper, both} */
    lower = i_tail != 1;
    upper = i_tail != 0;

    y = fabs(x);
    if (y <= 0.67448975) { /* qnorm(3/4) = .6744.... -- earlier had 0.66291 */
	if (y > eps) {
	    xsq = x * x;
	    xnum = a[4] * xsq;
	    xden = xsq;
	    for (i = 0; i < 3; ++i) {
		xnum = (xnum + a[i]) * xsq;
		xden = (xden + b[i]) * xsq;
	    }
	} else xnum = xden = 0.0;

	temp = x * (xnum + a[3]) / (xden + b[3]);
	if(lower)  cum = 0.5 + temp;
	if(upper) ccum = 0.5 - temp;
	}    
    else if (y <= M_SQRT_32) {

	/* Evaluate pnorm for 0.674.. = qnorm(3/4) < |x| <= SQRT(32) ~= 5.657 */

	xnum = c[8] * y;
	xden = y;
	for (i = 0; i < 7; ++i) {
	    xnum = (xnum + c[i]) * y;
	    xden = (xden + d[i]) * y;
	}
	temp = (xnum + c[7]) / (xden + d[7]);

#define do_del(X)							\
	xsq = floor(X * SIXTEN) / SIXTEN;			\
	del = (X - xsq) * (X + xsq);					\
	cum = exp(-xsq * xsq * 0.5) * exp(-del * 0.5) * temp;		\
	ccum = 1.0 - cum;						\
	
#define swap_tail						\
	if (x > 0.) {/* swap  ccum <--> cum */			\
	    temp = cum; if(lower) cum = ccum; ccum = temp;	\
	}

	do_del(y);
	swap_tail;
    }

/* else	  |x| > SQRT(32) = 5.657 :
 * the next two case differentiations were really for lower=T, log=F
 * Particularly	 *not*	for  log_p !

 * Cody had (-37.5193 < x  &&  x < 8.2924) ; R originally had y < 50
 *
 * Note that we do want symmetry(0), lower/upper -> hence use y
 */
    else if((lower && -37.5193 < x  &&  x < 8.2924) || (upper && -8.2924  < x  &&  x < 37.5193)) 
      {
	/* Evaluate pnorm for x in (-37.5, -5.657) union (5.657, 37.5) */
	xsq = 1.0 / (x * x);
	xnum = p[5] * xsq;
	xden = xsq;
	for (i = 0; i < 4; ++i) {
	    xnum = (xnum + p[i]) * xsq;
	    xden = (xden + q[i]) * xsq;
	}
	temp = xsq * (xnum + p[4]) / (xden + q[4]);
	temp = (M_1_SQRT_2_PI - temp) / y;

	do_del(x);
	swap_tail;
      }
    else 
      { /* no log_p , large x such that probs are 0 or 1 */
	if(x > 0) {	cum = 1.; ccum = 0.;	}
	else {	        cum = 0.; ccum = 1.;	}
      }

    return (phydbl)cum;


}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl Pgamma(phydbl x, phydbl shape, phydbl scale)
{
  return IncompleteGamma(x/scale,shape,LnGamma(shape));
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl Ppois(phydbl x, phydbl param)
{
  /* Press et al. (1990) approximation of the CDF for the Poisson distribution */
  if(param < SMALL || x < 0.0) 
    {
      PhyML_Printf("\n== param = %G x=%G",param,x);
      PhyML_Printf("\n== Err. in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
  return IncompleteGamma(x,param,LnGamma(param));
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* Inverse CDFs */
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl PointChi2 (phydbl prob, phydbl v)
{
/* returns z so that Prob{x<z}=prob where x is Chi2 distributed with df=v
   returns -1 if in error.   0.000002<prob<0.999998
   RATNEST FORTRAN by
       Best DJ & Roberts DE (1975) The percentage points of the
       Chi2 distribution.  Applied Statistics 24: 385-388.  (AS91)
   Converted into C by Ziheng Yang, Oct. 1993.
*/
   double aa=.6931471805, p=prob, g;
   double xx, c, ch, a=0,q=0,p1=0,p2=0,t=0,x=0,b=0,s1,s2,s3,s4,s5,s6;
   double e=.5e-6;
   
   if (p<.000002 || p>.999998 || v<=0) return ((phydbl)-1);

   g = (double)LnGamma(v/2);
   xx=v/2;   c=xx-1;
   if (v >= -1.24*log(p)) goto l1;

   ch=pow((p*xx*exp(g+xx*aa)), 1/xx);
   if (ch-e<0) return (ch);
   goto l4;
l1:
   if (v>.32) goto l3;
   ch=0.4;   a=log(1-p);
l2:
   q=ch;  p1=1+ch*(4.67+ch);  p2=ch*(6.73+ch*(6.66+ch));
   t=-0.5+(4.67+2*ch)/p1 - (6.73+ch*(13.32+3*ch))/p2;
   ch-=(1-exp(a+g+.5*ch+c*aa)*p2/p1)/t;
   if (fabs(q/ch-1)-.01 <= 0) goto l4;
   else                       goto l2;

l3:
   x=(double)PointNormal (p);
   p1=0.222222/v;   ch=v*pow((x*sqrt(p1)+1-p1), 3.0);
   if (ch>2.2*v+6)  ch=-2*(log(1-p)-c*log(.5*ch)+g);
l4:
   q=ch;   p1=.5*ch;
   if ((t=(double)IncompleteGamma (p1, xx, g))<0) {
      PhyML_Printf ("\nerr IncompleteGamma");
      return ((phydbl)-1.);
   }
   p2=p-t;
   t=p2*exp(xx*aa+g+p1-c*log(ch));
   b=t/ch;  a=0.5*t-b*c;

   s1=(210+a*(140+a*(105+a*(84+a*(70+60*a))))) / 420;
   s2=(420+a*(735+a*(966+a*(1141+1278*a))))/2520;
   s3=(210+a*(462+a*(707+932*a)))/2520;
   s4=(252+a*(672+1182*a)+c*(294+a*(889+1740*a)))/5040;
   s5=(84+264*a+c*(175+606*a))/2520;
   s6=(120+c*(346+127*c))/5040;
   ch+=t*(1+0.5*t*s1-b*c*(s1-b*(s2-b*(s3-b*(s4-b*(s5-b*s6))))));
   if (FABS(q/ch-1) > e) goto l4;

   return (phydbl)(ch);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////



/*
  The following function was extracted from the source code of R.
  It implements the methods referenced below.
 *  REFERENCE
 *
 *	Beasley, J. D. and S. G. Springer (1977).
 *	Algorithm AS 111: The percentage points of the normal distribution,
 *	Applied Statistics, 26, 118-121.
 *
 *      Wichura, M.J. (1988).
 *      Algorithm AS 241: The Percentage Points of the Normal Distribution.
 *      Applied Statistics, 37, 477-484.
 */


phydbl PointNormal (phydbl prob)
{
/* returns z so that Prob{x<z}=prob where x ~ N(0,1) and (1e-12)<prob<1-(1e-12)
   returns (-9999) if in error
   Odeh RE & Evans JO (1974) The percentage points of the normal distribution.
   Applied Statistics 22: 96-97 (AS70)

   Newer methods:
     Wichura MJ (1988) Algorithm AS 241: the percentage points of the
       normal distribution.  37: 477-484.
     Beasley JD & Springer SG  (1977).  Algorithm AS 111: the percentage 
       points of the normal distribution.  26: 118-121.
*/
   phydbl a0=-.322232431088, a1=-1, a2=-.342242088547, a3=-.0204231210245;
   phydbl a4=-.453642210148e-4, b0=.0993484626060, b1=.588581570495;
   phydbl b2=.531103462366, b3=.103537752850, b4=.0038560700634;
   phydbl y, z=0, p=prob, p1;

   p1 = (p<0.5 ? p : 1-p);
   if (p1<1e-20) z=999;
   else {
      y = SQRT (LOG(1/(p1*p1)));   
      z = y + ((((y*a4+a3)*y+a2)*y+a1)*y+a0) / ((((y*b4+b3)*y+b2)*y+b1)*y+b0);
   }
   return (p<0.5 ? -z : z);
}


/* phydbl PointNormal(phydbl p) */
/* { */
/*     double p_, q, r, val; */

/*     p_ = p; */
/*     q = p_ - 0.5; */

/*     /\*-- use AS 241 --- *\/ */
/*     /\* double ppnd16_(double *p, long *ifault)*\/ */
/*     /\*      ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3 */

/* 	    Produces the normal deviate Z corresponding to a given lower */
/* 	    tail area of P; Z is accurate to about 1 part in 10**16. */

/* 	    (original fortran code used PARAMETER(..) for the coefficients */
/* 	    and provided hash codes for checking them...) */
/* *\/ */
/*     if (fabs(q) <= .425)  */
/*       {/\* 0.075 <= p <= 0.925 *\/ */
/* 	r = .180625 - q * q; */
/* 	val = */
/* 	  q * (((((((r * 2509.0809287301226727 + */
/* 		     33430.575583588128105) * r + 67265.770927008700853) * r + */
/* 		   45921.953931549871457) * r + 13731.693765509461125) * r + */
/* 		 1971.5909503065514427) * r + 133.14166789178437745) * r + */
/* 	       3.387132872796366608) */
/* 	  / (((((((r * 5226.495278852854561 + */
/* 		   28729.085735721942674) * r + 39307.89580009271061) * r + */
/* 		 21213.794301586595867) * r + 5394.1960214247511077) * r + */
/* 	       687.1870074920579083) * r + 42.313330701600911252) * r + 1.); */
/*       } */
/*     else  */
/*       { /\* closer than 0.075 from {0,1} boundary *\/ */

/* 	/\* r = min(p, 1-p) < 0.075 *\/ */
/* 	if (q > 0) */
/* 	  r = 1-p;/\* 1-p *\/ */
/* 	else */
/* 	  r = p_;/\* = R_DT_Iv(p) ^=  p *\/ */
	
/* 	r = sqrt(-log(r)); */
/*         /\* r = sqrt(-log(r))  <==>  min(p, 1-p) = exp( - r^2 ) *\/ */
	
/*         if (r <= 5.) { /\* <==> min(p,1-p) >= exp(-25) ~= 1.3888e-11 *\/ */
/* 	  r += -1.6; */
/* 	  val = (((((((r * 7.7454501427834140764e-4 + */
/*                        .0227238449892691845833) * r + .24178072517745061177) * */
/*                      r + 1.27045825245236838258) * r + */
/*                     3.64784832476320460504) * r + 5.7694972214606914055) * */
/*                   r + 4.6303378461565452959) * r + */
/*                  1.42343711074968357734) */
/* 	    / (((((((r * */
/* 		     1.05075007164441684324e-9 + 5.475938084995344946e-4) * */
/* 		    r + .0151986665636164571966) * r + */
/* 		   .14810397642748007459) * r + .68976733498510000455) * */
/* 		 r + 1.6763848301838038494) * r + */
/* 		2.05319162663775882187) * r + 1.); */
/*         } */
/*         else  */
/* 	  { /\* very close to  0 or 1 *\/ */
/* 	    r += -5.; */
/* 	    val = (((((((r * 2.01033439929228813265e-7 + */
/* 			 2.71155556874348757815e-5) * r + */
/* 			.0012426609473880784386) * r + .026532189526576123093) * */
/* 		      r + .29656057182850489123) * r + */
/* 		     1.7848265399172913358) * r + 5.4637849111641143699) * */
/* 		   r + 6.6579046435011037772) */
/* 	      / (((((((r * */
/* 		       2.04426310338993978564e-15 + 1.4215117583164458887e-7)* */
/* 		      r + 1.8463183175100546818e-5) * r + */
/* 		     7.868691311456132591e-4) * r + .0148753612908506148525) */
/* 		   * r + .13692988092273580531) * r + */
/* 		  .59983220655588793769) * r + 1.); */
/* 	  } */
	
/* 	if(q < 0.0) */
/* 	  val = -val; */
/*         /\* return (q >= 0.)? r : -r ;*\/ */
/*       } */
/*     return (phydbl)val; */
/* } */

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* MISCs */
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Bico(int n, int k)
{
  return FLOOR(0.5+EXP(Factln(n)-Factln(k)-Factln(n-k)));
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl Factln(int n)
{
  static phydbl a[101];
  
  if (n < 0)    { Warn_And_Exit("\n== Err: negative factorial in routine FACTLN"); }
  if (n <= 1)     return 0.0;
  if (n <= 100)   return (a[n]>SMALL) ? a[n] : (a[n]=Gammln(n+1.0));
  else return     Gammln(n+1.0);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Gammln(phydbl xx)
{
  double x,tmp,ser;
  static double cof[6]={76.18009173,-86.50532033,24.01409822,
			-1.231739516,0.120858003e-2,-0.536382e-5};
  int j;
  
  x=xx-1.0;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.0;
  for (j=0;j<=5;j++) 
    {
      x += 1.0;
      ser += cof[j]/x;
    }
  return (phydbl)(-tmp+log(2.50662827465*ser));
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


/* void Plim_Binom(phydbl pH0, int N, phydbl *pinf, phydbl *psup) */
/* { */
/*   *pinf = pH0 - 1.64*SQRT(pH0*(1-pH0)/(phydbl)N); */
/*   if(*pinf < 0) *pinf = .0; */
/*   *psup = pH0 + 1.64*SQRT(pH0*(1-pH0)/(phydbl)N); */
/* } */

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl LnGamma (phydbl alpha)
{
/* returns ln(gamma(alpha)) for alpha>0, accurate to 10 decimal places.
   Stirling's formula is used for the central polynomial part of the procedure.
   Pike MC & Hill ID (1966) Algorithm 291: Logarithm of the gamma function.
   Communications of the Association for Computing Machinery, 9:684
*/
   double x=alpha, f=0, z;
   if (x<7) {
      f=1;  z=x-1;
      while (++z<7)  f*=z;
      x=z;   f=-log(f);
   }
   z = 1/(x*x);
   return (phydbl)(f + (x-0.5)*log(x) - x + .918938533204673
		   + (((-.000595238095238*z+.000793650793651)*z-.002777777777778)*z
		      +.083333333333333)/x);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl IncompleteGamma(phydbl x, phydbl alpha, phydbl ln_gamma_alpha)
{
/* returns the incomplete gamma ratio I(x,alpha) where x is the upper
	   limit of the integration and alpha is the shape parameter.
   returns (-1) if in error
   ln_gamma_alpha = ln(Gamma(alpha)), is almost redundant.
   (1) series expansion     if (alpha>x || x<=1)
   (2) continued fraction   otherwise
   RATNEST FORTRAN by
   Bhattacharjee GP (1970) The incomplete gamma integral.  Applied Statistics,
   19: 285-287 (AS32)
*/
   int i;
   double p=alpha, g=ln_gamma_alpha;
   double accurate=1e-8, overflow=1e30;
   double factor, gin=0, rn=0, a=0,b=0,an=0,dif=0, term=0, pn[6];

   if (fabs(x) < SMALL) return ((phydbl).0);
   if (x<0 || p<=0)        return ((phydbl)-1);

   factor=exp(p*log(x)-x-g);
   if (x>1 && x>=p) goto l30;
   /* (1) series expansion */
   gin=1;  term=1;  rn=p;
 l20:
   rn++;
   term*=x/rn;   gin+=term;

   if (term > accurate) goto l20;
   gin*=factor/p;
   goto l50;
 l30:
   /* (2) continued fraction */
   a=1-p;   b=a+x+1;  term=0;
   pn[0]=1;  pn[1]=x;  pn[2]=x+1;  pn[3]=x*b;
   gin=pn[2]/pn[3];
 l32:
   a++;  b+=2;  term++;   an=a*term;
   for (i=0; i<2; i++) pn[i+4]=b*pn[i+2]-an*pn[i];
   if (fabs(pn[5]) < .0) goto l35;
   rn=pn[4]/pn[5];   dif=fabs(gin-rn);
   if (dif>accurate) goto l34;
   if (dif<=accurate*rn) goto l42;
 l34:
   gin=rn;
 l35:
   for (i=0; i<4; i++) pn[i]=pn[i+2];
   if (fabs(pn[4]) < overflow) goto l32;
   for (i=0; i<4; i++) pn[i]/=overflow;
   goto l32;
 l42:
   gin=1-factor*gin;

 l50:
   return (phydbl)(gin);
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


int DiscreteGamma (phydbl freqK[], phydbl rK[],
		   phydbl alfa, phydbl beta, int K, int median)
{
  /* discretization of gamma distribution with equal proportions in each
     category
  */
   
  int i;
  phydbl gap05=1.0/(2.0*K), t, factor=alfa/beta*K, lnga1;

  if(K==1)
    {
      freqK[0] = 1.0;
      rK[0] = 1.0;
      return 0;
    }

   if (median) 
     {
       for (i=0; i<K; i++)     rK[i]=PointGamma((i*2.0+1)*gap05, alfa, beta);
       for (i=0,t=0; i<K; i++) t+=rK[i];
       for (i=0; i<K; i++)     rK[i]*=factor/t;
     }
   else {

      lnga1=LnGamma(alfa+1);
      for (i=0; i<K-1; i++)
        {
          freqK[i]=PointGamma((i+1.0)/K, alfa, beta);
        }

      for (i=0; i<K-1; i++)
        {
          freqK[i]=IncompleteGamma(freqK[i]*beta, alfa+1, lnga1);
        }

      rK[0] = freqK[0]*factor;
      rK[K-1] = (1-freqK[K-2])*factor;
      for (i=1; i<K-1; i++)  rK[i] = (freqK[i]-freqK[i-1])*factor;
   }
   for (i=0; i<K; i++) freqK[i]=1.0/K;
   return (0);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


/* Return LOG(n!) */

phydbl LnFact(int n)
{
  int i;
  phydbl res;

  res = .0;
  for(i=2;i<n+1;i++) res += LOG((phydbl)i);
  
  return(res);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int Choose(int n, int k)
{
  phydbl accum;
  int i;

  if (k > n) return(0);
  if (k > n/2) k = n-k;
  if(!k) return(1);

  accum = 1.;
  for(i=1;i<k+1;i++) accum = accum * (n-k+i) / i;

  return((int)accum);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl LnChoose(int n, int k)
{
  return(LnFact(n) - LnFact(k) - LnFact(n-k));
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl *Covariance_Matrix(t_tree *tree)
{
  phydbl *cov, *mean;
  int *ori_wght,*site_num;
  int dim,i,j,replicate,n_site,position,sample_size;

  sample_size = 100;
  dim = 2*tree->n_otu-3;

  cov      = (phydbl *)mCalloc(dim*dim,sizeof(phydbl));
  mean     = (phydbl *)mCalloc(    dim,sizeof(phydbl));
  ori_wght = (int *)mCalloc(tree->data->crunch_len,sizeof(int));
  site_num = (int *)mCalloc(tree->data->init_len,sizeof(int));
  
  For(i,tree->data->crunch_len) ori_wght[i] = tree->data->wght[i];

  n_site = 0;
  For(i,tree->data->crunch_len) For(j,tree->data->wght[i])
    {
      site_num[n_site] = i;
      n_site++;
    }

  
  tree->mod->s_opt->print = 0;
  For(replicate,sample_size)
    {
      For(i,2*tree->n_otu-3) tree->a_edges[i]->l->v = .1;

      For(i,tree->data->crunch_len) tree->data->wght[i] = 0;

      For(i,tree->data->init_len)
	{
	  position = Rand_Int(0,(int)(tree->data->init_len-1.0));
	  tree->data->wght[site_num[position]] += 1;
	}

      Round_Optimize(tree,tree->data,ROUND_MAX);
      
      For(i,2*tree->n_otu-3) For(j,2*tree->n_otu-3) cov[i*dim+j] += LOG(tree->a_edges[i]->l->v) * LOG(tree->a_edges[j]->l->v);  
      For(i,2*tree->n_otu-3) mean[i] += LOG(tree->a_edges[i]->l->v);

      PhyML_Printf("[%3d/%3d]",replicate,sample_size); fflush(NULL);
/*       PhyML_Printf("\n. %3d %12f %12f %12f ", */
/* 	     replicate, */
/*  	     cov[1*dim+1]/(replicate+1)-mean[1]*mean[1]/POW(replicate+1,2), */
/* 	     tree->a_edges[1]->l->v, */
/* 	     mean[1]/(replicate+1)); */
   }

  For(i,2*tree->n_otu-3) mean[i] /= (phydbl)sample_size;
  For(i,2*tree->n_otu-3) For(j,2*tree->n_otu-3) cov[i*dim+j] /= (phydbl)sample_size;
  For(i,2*tree->n_otu-3) For(j,2*tree->n_otu-3) cov[i*dim+j] -= mean[i]*mean[j];
/*   For(i,2*tree->n_otu-3) if(cov[i*dim+i] < var_min) cov[i*dim+i] = var_min; */
  

/*   PhyML_Printf("\n"); */
/*   For(i,2*tree->n_otu-3) PhyML_Printf("%f %f\n",mean[i],tree->a_edges[i]->l->v); */
/*   PhyML_Printf("\n"); */
/*   PhyML_Printf("\n"); */
/*   For(i,2*tree->n_otu-3) */
/*     { */
/*       For(j,2*tree->n_otu-3) */
/* 	{ */
/* 	  PhyML_Printf("%G\n",cov[i*dim+j]); */
/* 	} */
/*       PhyML_Printf("\n"); */
/*     } */

  For(i,tree->data->crunch_len) tree->data->wght[i] = ori_wght[i];

  Free(mean);
  Free(ori_wght);
  Free(site_num);

  return cov;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* Work out the Hessian for the likelihood function. Only branch lengths are considered as variable.
   This function is very much inspired from Jeff Thorne's 'hessian' function in his program 'estbranches'. */
phydbl *Hessian(t_tree *tree)
{
  phydbl *hessian;
  phydbl *plus_plus, *minus_minus, *plus_zero, *minus_zero, *plus_minus, zero_zero;
  phydbl *ori_bl,*inc,*buff;
  int *ok_edges,*is_ok;
  int dim;
  int n_ok_edges;
  int i,j;
  phydbl eps;
  phydbl lk;
  phydbl lnL,lnL1,lnL2,ori_lnL;
  phydbl l_inf;

  dim = 2*tree->n_otu-3;
  eps = (tree->mod->log_l == YES)?(0.2):(1E-4);

  hessian     = (phydbl *)mCalloc((int)dim*dim,sizeof(phydbl));
  ori_bl      = (phydbl *)mCalloc((int)dim,sizeof(phydbl));
  plus_plus   = (phydbl *)mCalloc((int)dim*dim,sizeof(phydbl));
  minus_minus = (phydbl *)mCalloc((int)dim*dim,sizeof(phydbl));
  plus_minus  = (phydbl *)mCalloc((int)dim*dim,sizeof(phydbl));
  plus_zero   = (phydbl *)mCalloc((int)dim    ,sizeof(phydbl));
  minus_zero  = (phydbl *)mCalloc((int)dim    ,sizeof(phydbl));
  inc         = (phydbl *)mCalloc((int)dim    ,sizeof(phydbl));
  buff        = (phydbl *)mCalloc((int)dim*dim,sizeof(phydbl));
  ok_edges    = (int *)mCalloc((int)dim,sizeof(int));
  is_ok       = (int *)mCalloc((int)dim,sizeof(int));

  lnL = lnL1 = lnL2 = UNLIKELY;

  Set_Both_Sides(YES,tree);
  Lk(NULL,tree);
  ori_lnL = tree->c_lnL;


  For(i,dim) ori_bl[i] = tree->a_edges[i]->l->v;

  if(tree->mod->log_l == NO)
    l_inf = MAX(tree->mod->l_min,1./(phydbl)tree->data->init_len);
  else
    l_inf = MAX(tree->mod->l_min,-LOG((phydbl)tree->data->init_len));


  n_ok_edges = 0;
  For(i,dim) 
    {
      if(tree->a_edges[i]->l->v*(1.-eps) > l_inf)
	{
	  inc[i] = eps * tree->a_edges[i]->l->v;
	  ok_edges[n_ok_edges] = i;
	  n_ok_edges++;
	  is_ok[i] = 1;
	}
      else
	{
	  inc[i] = -1.0;
	  is_ok[i] = 0;
	}
    }


  /* Fine tune the increments */
  For(i,dim)
    {
      do
	{
	  tree->a_edges[i]->l->v += inc[i];
	  lnL1 = Lk(tree->a_edges[i],tree);
	  tree->a_edges[i]->l->v = ori_bl[i];
	  inc[i] *= 1.1;
	}while((FABS(lnL1 - ori_lnL) < 1.E-1) && 
	       (tree->a_edges[i]->l->v+inc[i] < tree->mod->l_max));
      inc[i] /= 1.1;
    }



  /* zero zero */  
  zero_zero = tree->c_lnL;

  /* plus zero */  
  For(i,dim) 
    {
      if(is_ok[i])
	{
	  tree->a_edges[i]->l->v += inc[i];
	  lk = Lk(tree->a_edges[i],tree);
	  plus_zero[i] = lk;
	  tree->a_edges[i]->l->v = ori_bl[i];
	}
    }


  /* minus zero */  
  For(i,dim) 
    {
      if(is_ok[i])
	{
	  tree->a_edges[i]->l->v -= inc[i];
	  lk = Lk(tree->a_edges[i],tree);
	  minus_zero[i] = lk;
	  tree->a_edges[i]->l->v = ori_bl[i];
	}
    }


  For(i,dim) Update_PMat_At_Given_Edge(tree->a_edges[i],tree);

  /* plus plus  */  
  For(i,dim)
    {
      if(is_ok[i])
	{
	  tree->a_edges[i]->l->v += inc[i];
	  Update_PMat_At_Given_Edge(tree->a_edges[i],tree);

	  For(j,3)
	    if((!tree->a_edges[i]->left->tax) && (tree->a_edges[i]->left->v[j] != tree->a_edges[i]->rght))
	      Recurr_Hessian(tree->a_edges[i]->left,tree->a_edges[i]->left->v[j],1,inc,plus_plus+i*dim,is_ok,tree);
	  
	  For(j,3)
	    if((!tree->a_edges[i]->rght->tax) && (tree->a_edges[i]->rght->v[j] != tree->a_edges[i]->left))
	      Recurr_Hessian(tree->a_edges[i]->rght,tree->a_edges[i]->rght->v[j],1,inc,plus_plus+i*dim,is_ok,tree);
		      
	  tree->a_edges[i]->l->v = ori_bl[i];
	  Lk(NULL,tree);
	}
    }


  /* plus minus */  
  For(i,dim)
    {
      if(is_ok[i])
	{
	  tree->a_edges[i]->l->v += inc[i];
	  Update_PMat_At_Given_Edge(tree->a_edges[i],tree);
	  
	  For(j,3)
	    if((!tree->a_edges[i]->left->tax) && (tree->a_edges[i]->left->v[j] != tree->a_edges[i]->rght))
	      Recurr_Hessian(tree->a_edges[i]->left,tree->a_edges[i]->left->v[j],-1,inc,plus_minus+i*dim,is_ok,tree);
	  
	  For(j,3)
	    if((!tree->a_edges[i]->rght->tax) && (tree->a_edges[i]->rght->v[j] != tree->a_edges[i]->left))
	      Recurr_Hessian(tree->a_edges[i]->rght,tree->a_edges[i]->rght->v[j],-1,inc,plus_minus+i*dim,is_ok,tree);
	  
	  tree->a_edges[i]->l->v = ori_bl[i];
	  Lk(NULL,tree);
	}
    }


  /* minus minus */  
  For(i,dim)
    {
      if(is_ok[i])
	{
	  tree->a_edges[i]->l->v -= inc[i];
	  
	  Update_PMat_At_Given_Edge(tree->a_edges[i],tree);
	  
	  For(j,3)
	    if((!tree->a_edges[i]->left->tax) && (tree->a_edges[i]->left->v[j] != tree->a_edges[i]->rght))
	      Recurr_Hessian(tree->a_edges[i]->left,tree->a_edges[i]->left->v[j],-1,inc,minus_minus+i*dim,is_ok,tree);
	  
	  For(j,3)
	    if((!tree->a_edges[i]->rght->tax) && (tree->a_edges[i]->rght->v[j] != tree->a_edges[i]->left))
	      Recurr_Hessian(tree->a_edges[i]->rght,tree->a_edges[i]->rght->v[j],-1,inc,minus_minus+i*dim,is_ok,tree);
	  
	  tree->a_edges[i]->l->v = ori_bl[i];
	  Lk(NULL,tree);
	}
    }


  
  For(i,dim)
    {
      if(is_ok[i])
	{
	  hessian[i*dim+i] = (plus_zero[i]-2*zero_zero+minus_zero[i])/(POW(inc[i],2));

	  for(j=i+1;j<dim;j++)
	    {
	      if(is_ok[j])
		{
		  hessian[i*dim+j] = 
		    (plus_plus[i*dim+j]-plus_minus[i*dim+j]-plus_minus[j*dim+i]+minus_minus[i*dim+j])/
		    (4*inc[i]*inc[j]);
		  hessian[j*dim+i] = hessian[i*dim+j];
		}
	    }
	}
    }
        
  For(i,n_ok_edges)
    {
      For(j,n_ok_edges)
	{
	  buff[i*n_ok_edges+j] = -1.0*hessian[ok_edges[i]*dim+ok_edges[j]];
	}
    }


  if(!Matinv(buff,n_ok_edges,n_ok_edges,NO))
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");      
    }

  For(i,n_ok_edges)
    {
      For(j,n_ok_edges)
	{
	  hessian[ok_edges[i]*dim+ok_edges[j]] = buff[i*n_ok_edges+j];
	}
    }

/*   eps = 1./(phydbl)tree->data->init_len; */
  /* Approximate variance for very short branches */
  For(i,dim)
    if(inc[i] < 0.0 || hessian[i*dim+i] < MIN_VAR_BL)
      {	
	eps = 0.2 * tree->a_edges[i]->l->v;
	do
	  {
	    lnL  = Lk(tree->a_edges[i],tree);	
	    tree->a_edges[i]->l->v += eps;
	    lnL1 = Lk(tree->a_edges[i],tree);
	    tree->a_edges[i]->l->v += eps;
	    lnL2 = Lk(tree->a_edges[i],tree);
	    tree->a_edges[i]->l->v -= 2.*eps;
	    
	    hessian[i*dim+i] = (lnL2 - 2.*lnL1 + lnL) / POW(eps,2);
	
/* 	    printf("\n* l=%G eps=%f lnL=%f lnL1=%f lnL2=%f var=%f",tree->a_edges[i]->l->v,eps,lnL,lnL1,lnL2,hessian[i*dim+i]); */
	    eps *= 5.;
	  }while(FABS(lnL2 - lnL) < 1.E-3);

	hessian[i*dim+i] = -1.0 / hessian[i*dim+i];

      }
  

  /* Fit a straight line to the log-likelihood (i.e., an exponential to the likelihood) */
  /* It is only a straight line when considering branch length (rather than log(branch lengths)) */
  For(i,dim)
    if((tree->a_edges[i]->l->v / tree->mod->l_min < 1.1) &&
       (tree->a_edges[i]->l->v / tree->mod->l_min > 0.9))
      {
	phydbl *x,*y,l;
	phydbl cov,var;
	
	x=plus_plus;
	y=minus_minus;
	l=(tree->mod->log_l == YES)?(EXP(tree->a_edges[i]->l->v)):(tree->a_edges[i]->l->v); /* Get actual branch length */
	
	For(j,dim)
	  {
	    x[j] = l + (100.*l-l)*((phydbl)j/dim);
	    tree->a_edges[i]->l->v = (tree->mod->log_l)?(LOG(x[j])):(x[j]); /* Transform to log if necessary */
	    y[j] = Lk(tree->a_edges[i],tree);
	    tree->a_edges[i]->l->v = (tree->mod->log_l)?(LOG(l)):(l); /* Go back to initial edge length */
	  }
	
	cov = Covariance(x,y,dim);
	var = Covariance(x,x,dim);
	
	/* cov/var is minus the parameter of the exponential distribution.
	   The variance is therefore : */
	hessian[i*dim+i] = 1.0 / pow(cov/var,2);
	
	/* 	    printf("\n. Hessian = %G cov=%G var=%G",hessian[i*dim+i],cov,var); */
      }
  /*     } */


  For(i,dim)
    if(hessian[i*dim+i] < 0.0)
      {
	PhyML_Printf("\n. l=%G var=%G",tree->a_edges[i]->l->v,hessian[i*dim+i]);
/* 	PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__); */
/* 	Exit("\n"); */
	hessian[i*dim+i] = MIN_VAR_BL;
      }

  For(i,dim)
    {
      if(hessian[i*dim+i] < MIN_VAR_BL)
	{
	  PhyML_Printf("\n. l=%10G var(l)=%12G. WARNING: numerical precision issues may affect this analysis.",
		       tree->a_edges[i]->l->v,hessian[i*dim+i]);
	  hessian[i*dim+i] = MIN_VAR_BL;
	}
      if(hessian[i*dim+i] > MAX_VAR_BL)
	{
	  PhyML_Printf("\n. l=%10G var(l)=%12G. WARNING: numerical precision issues may affect this analysis.",
		       tree->a_edges[i]->l->v,hessian[i*dim+i]);
	  hessian[i*dim+i] = MAX_VAR_BL;
	}
    }
  
  Iter_Matinv(hessian,dim,dim,NO);

  For(i,dim*dim) hessian[i] = -1.0*hessian[i];

  For(i,dim)
    {
      For(j,dim)
	{
	  if(FABS(hessian[i*dim+j]-hessian[j*dim+i]) > 1.E-3)
	    {
	      PhyML_Printf("\n. Hessian not symmetrical.");
	      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	      Exit("\n");
	    }
	  hessian[i*dim+j] = (hessian[i*dim+j] + hessian[j*dim+i]) / 2.; 
	  hessian[j*dim+i] = hessian[i*dim+j];  
	}
    }
  
/*   printf("\n"); */
/*   printf("HESSIAN\n"); */
/*   For(i,dim) */
/*     { */
/*       PhyML_Printf("[%f] ",tree->a_edges[i]->l->v); */
/*       For(j,dim) */
/* 	{ */
/* 	  PhyML_Printf("%12lf ",hessian[i*dim+j]); */
/* 	} */
/*       PhyML_Printf("\n"); */
/*     } */

  /* Matinv(hessian,dim,dim,NO); */

  /* PhyML_Printf("\n"); */

  /* For(i,dim) */
  /*   { */
  /*     PhyML_Printf("[%f] ",tree->a_edges[i]->l->v); */
  /*     For(j,dim) */
  /* 	{ */
  /* 	  PhyML_Printf("%12G ",-hessian[i*dim+j]); */
  /* 	} */
  /*     PhyML_Printf("\n"); */
  /*   } */
  /* Exit("\n"); */


  /* Make sure to update likelihood before bailing out */
  Set_Both_Sides(YES,tree);
  Lk(NULL,tree);

  Free(ori_bl);
  Free(plus_plus);
  Free(minus_minus);
  Free(plus_zero);
  Free(minus_zero);
  Free(plus_minus);
  Free(inc);
  Free(buff);
  Free(ok_edges);
  Free(is_ok);

  return hessian;

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


/* Work out the gradient for the likelihood function. Only branch lengths are considered as variable.
 */
phydbl *Gradient(t_tree *tree)
{
  phydbl *gradient;
  phydbl *plus, *minus;
  phydbl *ori_bl,*inc;
  int *is_ok;
  int dim;
  int i;
  phydbl eps;
  phydbl lk;
  phydbl lnL,lnL1,lnL2;
  phydbl l_inf;

  dim = 2*tree->n_otu-3;
  eps = (tree->mod->log_l == YES)?(0.2):(1.E-6);

  gradient    = (phydbl *)mCalloc((int)dim,sizeof(phydbl));
  ori_bl      = (phydbl *)mCalloc((int)dim,sizeof(phydbl));
  plus        = (phydbl *)mCalloc((int)dim,sizeof(phydbl));
  minus       = (phydbl *)mCalloc((int)dim,sizeof(phydbl));
  inc         = (phydbl *)mCalloc((int)dim    ,sizeof(phydbl));
  is_ok       = (int *)mCalloc((int)dim,sizeof(int));

  lnL = lnL1 = lnL2 = UNLIKELY;

  Set_Both_Sides(YES,tree);
  Lk(NULL,tree);

  For(i,dim) ori_bl[i] = tree->a_edges[i]->l->v;

  if(tree->mod->log_l == NO)
    l_inf = MAX(tree->mod->l_min,1./(phydbl)tree->data->init_len);
  else
    l_inf = MAX(tree->mod->l_min,-LOG((phydbl)tree->data->init_len));

  For(i,dim) 
    {
      if(tree->a_edges[i]->l->v*(1.-eps) > l_inf)
	{
	  inc[i] = eps * tree->a_edges[i]->l->v;
	  is_ok[i] = YES;
	}
      else
	{
	  inc[i] = -1.0;
	  is_ok[i] = NO;
	}
    }

  /* plus */  
  For(i,dim) 
    {
      if(is_ok[i] == YES)
	{
	  tree->a_edges[i]->l->v += inc[i];
	  lk = Lk(tree->a_edges[i],tree);
	  plus[i] = lk;
	  tree->a_edges[i]->l->v = ori_bl[i];
	}
    }


  /* minus */  
  For(i,dim)
    {
      if(is_ok[i] == YES)
	{
	  tree->a_edges[i]->l->v -= inc[i];
	  lk = Lk(tree->a_edges[i],tree);
	  minus[i] = lk;
	  tree->a_edges[i]->l->v = ori_bl[i];
	}
    }


  For(i,dim)
    {
      if(is_ok[i] == YES)
	{
	  gradient[i] = (plus[i] - minus[i])/(2.*inc[i]);
	}
    }


  For(i,dim)
    {
      if(is_ok[i] == NO)
	{
	  eps = FABS(0.2 * tree->a_edges[i]->l->v);
	  lnL  = Lk(tree->a_edges[i],tree);	
	  tree->a_edges[i]->l->v += eps;
	  lnL1 = Lk(tree->a_edges[i],tree);
	  tree->a_edges[i]->l->v += eps;
	  lnL2 = Lk(tree->a_edges[i],tree);
	  tree->a_edges[i]->l->v -= eps;
	  tree->a_edges[i]->l->v -= eps;
	  gradient[i] = (4.*lnL1 - lnL2 - 3.*lnL) / (2.*eps);
	}
    }

  /* Make sure to update likelihood before bailing out */
  Set_Both_Sides(YES,tree);
  Lk(NULL,tree);

  Free(ori_bl);
  Free(plus);
  Free(minus);
  Free(inc);
  Free(is_ok);
  

/*   printf("\n"); */
/*   printf("GRADIENT\n"); */
/*   For(i,dim) */
/*     { */
/*       PhyML_Printf("[%f] ",tree->a_edges[i]->l->v); */
/*       For(j,dim) */
/* 	{ */
/* 	  printf("%12lf ",gradient[i]*gradient[j]); */
/* 	} */
/*       printf("\n"); */
/*     } */
/*   printf("\n"); */
/*   For(i,dim) */
/*     { */
/*       PhyML_Printf("[%f] [%f]\n",tree->a_edges[i]->l->v,gradient[i]); */
/*     } */

/*   Exit("\n"); */

  return gradient;

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


/* Work out the Hessian for the likelihood function using the method described by Seo et al., 2004, MBE.
   Corresponds to the outer product of the scores approach described in Porter, 2002. (matrix J1)
*/
phydbl *Hessian_Seo(t_tree *tree)
{
  phydbl *hessian,*site_hessian;
  phydbl *gradient;
  phydbl *plus, *minus, *plusplus, *zero;
  phydbl *ori_bl,*inc_plus,*inc_minus,*inc;
  int *is_ok;
  int dim;
  int i,j,k;
  phydbl eps;
  phydbl ori_lnL,lnL1,lnL2;
  phydbl l_inf;
  int l,n;
  phydbl small_var;

  dim = 2*tree->n_otu-3;
  eps = (tree->mod->log_l == YES)?(0.2):(1.E-4);

  hessian      = (phydbl *)mCalloc((int)dim*dim,sizeof(phydbl));
  site_hessian = (phydbl *)mCalloc((int)dim*dim,sizeof(phydbl));
  gradient     = (phydbl *)mCalloc((int)dim,sizeof(phydbl));
  ori_bl       = (phydbl *)mCalloc((int)dim,sizeof(phydbl));
  plus         = (phydbl *)mCalloc((int)dim*tree->n_pattern,sizeof(phydbl));
  plusplus     = (phydbl *)mCalloc((int)dim*tree->n_pattern,sizeof(phydbl));
  minus        = (phydbl *)mCalloc((int)dim*tree->n_pattern,sizeof(phydbl));
  zero         = (phydbl *)mCalloc((int)dim*tree->n_pattern,sizeof(phydbl));
  inc_plus     = (phydbl *)mCalloc((int)dim,sizeof(phydbl));
  inc_minus    = (phydbl *)mCalloc((int)dim,sizeof(phydbl));
  inc          = (phydbl *)mCalloc((int)dim,sizeof(phydbl));
  is_ok        = (int *)mCalloc((int)dim,sizeof(int));

  lnL1 = lnL2 = UNLIKELY;
  
  For(i,dim) ori_bl[i] = tree->a_edges[i]->l->v;
  
  Set_Both_Sides(YES,tree);
  Lk(NULL,tree);
  ori_lnL = tree->c_lnL;


  if(tree->mod->log_l == NO)
    l_inf = MAX(tree->mod->l_min,1./(phydbl)tree->data->init_len);
  else
    l_inf = MAX(tree->mod->l_min,-LOG((phydbl)tree->data->init_len));

  For(i,dim) 
    {
      if(tree->a_edges[i]->l->v*(1.-eps) > l_inf)
	{
	  inc_plus[i]  = FABS(eps * MAX(tree->mod->l_min,tree->a_edges[i]->l->v));
	  inc_minus[i] = FABS(eps * MAX(tree->mod->l_min,tree->a_edges[i]->l->v));
	  is_ok[i]     = YES;
	}
      else
	{
	  inc_plus[i]  = FABS(0.2 * MAX(tree->mod->l_min,tree->a_edges[i]->l->v));
	  inc_minus[i] = FABS(0.2 * MAX(tree->mod->l_min,tree->a_edges[i]->l->v));
	  is_ok[i]     = NO;
    }


    }

  /* Fine tune the increments */
  For(i,dim)
    {
      do
	{
	  tree->a_edges[i]->l->v += inc_plus[i];	  
	  lnL1 = Lk(tree->a_edges[i],tree);
	  tree->a_edges[i]->l->v = ori_bl[i];
	  inc_plus[i] *= 1.1;
	}while((FABS(lnL1 - ori_lnL) < 1.E-1) && (tree->a_edges[i]->l->v+inc_plus[i] < tree->mod->l_max));
      inc_plus[i] /= 1.1;
    }

  For(i,dim)
    {
      do
	{
	  tree->a_edges[i]->l->v -= inc_minus[i];
	  lnL1 = Lk(tree->a_edges[i],tree);
	  tree->a_edges[i]->l->v = ori_bl[i];
	  inc_minus[i] *= 1.1;
	}while((FABS(lnL1 - ori_lnL) < 1.E-1) && 
	       (tree->a_edges[i]->l->v -inc_minus[i] > tree->mod->l_min));
      inc_minus[i] /= 1.1;
    }

  For(i,dim) 
    {
      inc[i] = MIN(inc_plus[i],inc_minus[i]);
    }

  /* plus */  
  For(i,dim) 
    {
      if(is_ok[i] == YES)
	{
	  tree->a_edges[i]->l->v += inc[i];
	  Lk(tree->a_edges[i],tree);
	  For(j,tree->n_pattern) plus[i*tree->n_pattern+j] = LOG(tree->cur_site_lk[j]);
	  tree->a_edges[i]->l->v = ori_bl[i];
	}
    }


  /* minus */
  For(i,dim)
    {
      if(is_ok[i] == YES)
	{
	  tree->a_edges[i]->l->v -= inc[i];
	  Lk(tree->a_edges[i],tree);
	  For(j,tree->n_pattern) minus[i*tree->n_pattern+j] = LOG(tree->cur_site_lk[j]);
	  tree->a_edges[i]->l->v = ori_bl[i];
	}
    }


  For(i,dim)
    {
      if(is_ok[i] == NO)
	{
	  Lk(tree->a_edges[i],tree);	
	  For(j,tree->n_pattern) zero[i*tree->n_pattern+j] = LOG(tree->cur_site_lk[j]);
	  
	  tree->a_edges[i]->l->v += inc[i];
	  lnL1 = Lk(tree->a_edges[i],tree);
	  For(j,tree->n_pattern) plus[i*tree->n_pattern+j] = LOG(tree->cur_site_lk[j]);
	  
	  tree->a_edges[i]->l->v += inc[i];
	  lnL2 = Lk(tree->a_edges[i],tree);
	  For(j,tree->n_pattern) plusplus[i*tree->n_pattern+j] = LOG(tree->cur_site_lk[j]);
	  
	  tree->a_edges[i]->l->v = ori_bl[i];	

	}
    }

  For(i,dim*dim) hessian[i] = 0.0;

  For(k,tree->n_pattern)
    {
      For(i,dim) 
	{
	  if(is_ok[i] == YES)
	    gradient[i] = (plus[i*tree->n_pattern+k] - minus[i*tree->n_pattern+k])/(inc[i] + inc[i]); 
	  else
	    gradient[i] = (4.*plus[i*tree->n_pattern+k] - plusplus[i*tree->n_pattern+k] - 3.*zero[i*tree->n_pattern+k])/(inc[i] + inc[i]);
	  
	  /* if(is_ok[i] == NO) */
	  /*   printf("\n. i=%d site=%d l=%G plus=%G plusplus=%G zero=%G num=%f grad=%G", */
	  /* 	   i,k,tree->a_edges[i]->l->v, */
	  /* 	   plus[i*tree->n_pattern+k],plusplus[i*tree->n_pattern+k],zero[i*tree->n_pattern+k], */
	  /* 	   (4.*plus[i*tree->n_pattern+k] - plusplus[i*tree->n_pattern+k] - 3.*zero[i*tree->n_pattern+k]), */
	  /* 	   gradient[i]); */

	}
      For(i,dim) For(j,dim) site_hessian[i*dim+j] = gradient[i] * gradient[j];
      For(i,dim*dim) hessian[i] -= site_hessian[i] * tree->data->wght[k]; 
    }


  /* Make sure to update likelihood before bailing out */
  Set_Both_Sides(YES,tree);
  Lk(NULL,tree);

  l = tree->data->init_len;
  n = tree->mod->ns;
  /* Delta method for variance. Assumes Jukes and Cantor with p=1/n */
  small_var = (1./(l*l))*(1.-1./l)*(n-1.)*(n-1.)/(n-1.-n/l);
  For(i,dim)
    if(is_ok[i] == NO)
      {
	For(j,dim)
	  {
	    hessian[i*dim+j] = 0.;
	    hessian[j*dim+i] = 0.;
	  }
	hessian[i*dim+i] = -1./small_var;

	if(tree->mod->log_l == YES) 
	  {
	    hessian[i*dim+i] = small_var * POW(EXP(tree->a_edges[i]->l->v),-2); 
	    hessian[i*dim+i] = -1./hessian[i*dim+i];
	  }
      }

  For(i,dim)
    if(is_ok[i] == YES && hessian[i*dim+i] < -1./small_var)
      {
  	For(j,dim)
  	  {
  	    hessian[i*dim+j] = 0.;
  	    hessian[j*dim+i] = 0.;
  	  }
  	hessian[i*dim+i] = -1./small_var;

  	if(tree->mod->log_l == YES)
  	  {
  	    hessian[i*dim+i] = small_var * POW(EXP(tree->a_edges[i]->l->v),-2);
  	    hessian[i*dim+i] = -1./hessian[i*dim+i];
  	  }
      }


  For(i,dim)
    {
      For(j,dim)
	{
	  if(FABS(hessian[i*dim+j]-hessian[j*dim+i]) > 1.E-3)
	    {
	      PhyML_Printf("\n== Hessian not symmetrical.");
	      PhyML_Printf("\n== Err in file %s at line %d\n",__FILE__,__LINE__);
	      Exit("\n");
	    }
	  hessian[i*dim+j] = (hessian[i*dim+j] + hessian[j*dim+i]) / 2.; 
	  hessian[j*dim+i] = hessian[i*dim+j];  
	}
    }

  /* printf("\n"); */
  /* printf("HESSIAN SEO\n"); */
  /* For(i,dim) */
  /*   { */
  /*     PhyML_Printf("[%f] ",tree->a_edges[i]->l->v); */
  /*     For(j,dim) */
  /* 	{ */
  /* 	  PhyML_Printf("%12lf ",hessian[i*dim+j]); */
  /* 	} */
  /*     PhyML_Printf("\n"); */
  /*   } */

  Free(site_hessian);
  Free(ori_bl);
  Free(plus);
  Free(minus);
  Free(plusplus);
  Free(zero);
  Free(inc);
  Free(inc_plus);
  Free(inc_minus);
  Free(is_ok);
  Free(gradient);
  
  return hessian;

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Recurr_Hessian(t_node *a, t_node *d, int plus_minus, phydbl *inc, phydbl *res, int *is_ok, t_tree *tree)
{
  int i;
  phydbl ori_l;

  For(i,3)
    if(a->v[i] == d)
      {
	Update_P_Lk(tree,a->b[i],a);

	ori_l = a->b[i]->l->v;
	if(is_ok[a->b[i]->num])
	  {
	    if(plus_minus > 0) a->b[i]->l->v += inc[a->b[i]->num];
	    else               a->b[i]->l->v -= inc[a->b[i]->num];
	    res[a->b[i]->num] = Lk(a->b[i],tree);
	    a->b[i]->l->v = ori_l;
	    Update_PMat_At_Given_Edge(a->b[i],tree);
	  }
	break;
      }

  if(d->tax) return;
  else 
    For(i,3) 
      if(d->v[i] != a) 
	Recurr_Hessian(d,d->v[i],plus_minus,inc,res,is_ok,tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


/* Work out the Hessian for the likelihood function. Only LOGARITHM of branch lengths are considered as variable.
   This function is very much inspired from Jeff Thorne's 'hessian' function in his program 'estbranches'. */
phydbl *Hessian_Log(t_tree *tree)
{
  phydbl *hessian;
  phydbl *plus_plus, *minus_minus, *plus_zero, *minus_zero, *plus_minus, *zero_zero;
  phydbl *ori_bl,*inc,*buff;
  int *ok_edges,*is_ok;
  int dim;
  int n_ok_edges;
  int i,j;
  phydbl eps;
  phydbl lk;

  dim = 2*tree->n_otu-3;
  eps = 1.E-4;

  hessian     = (phydbl *)mCalloc((int)dim*dim,sizeof(phydbl));
  ori_bl      = (phydbl *)mCalloc((int)dim,    sizeof(phydbl));
  plus_plus   = (phydbl *)mCalloc((int)dim*dim,sizeof(phydbl));
  minus_minus = (phydbl *)mCalloc((int)dim*dim,sizeof(phydbl));
  plus_minus  = (phydbl *)mCalloc((int)dim*dim,sizeof(phydbl));
  plus_zero   = (phydbl *)mCalloc((int)dim    ,sizeof(phydbl));
  minus_zero  = (phydbl *)mCalloc((int)dim    ,sizeof(phydbl));
  zero_zero   = (phydbl *)mCalloc((int)dim    ,sizeof(phydbl));
  inc         = (phydbl *)mCalloc((int)dim    ,sizeof(phydbl));
  buff        = (phydbl *)mCalloc((int)dim*dim,sizeof(phydbl));
  ok_edges    = (int *)mCalloc((int)dim,       sizeof(int));
  is_ok       = (int *)mCalloc((int)dim,       sizeof(int));
  
  Set_Both_Sides(YES,tree);
  Lk(NULL,tree);

  For(i,dim) ori_bl[i] = tree->a_edges[i]->l->v;

  n_ok_edges = 0;
  For(i,dim) 
    {
      if(tree->a_edges[i]->l->v > 3.0/(phydbl)tree->data->init_len)
	{
	  inc[i] = FABS(eps * tree->a_edges[i]->l->v);
	  ok_edges[n_ok_edges] = i;
	  n_ok_edges++;
	  is_ok[i] = 1;
	}
      else is_ok[i] = 0;
    }

  /* zero zero */  
  lk = Log_Det(is_ok,tree);
  For(i,dim) if(is_ok[i]) zero_zero[i] = tree->c_lnL+lk;

  /* plus zero */  
  For(i,dim) 
    {
      if(is_ok[i])
	{
	  tree->a_edges[i]->l->v += inc[i];
	  lk = Lk(tree->a_edges[i],tree);
	  plus_zero[i] = lk+Log_Det(is_ok,tree);
	  tree->a_edges[i]->l->v = ori_bl[i];
	}
    }


  /* minus zero */
  For(i,dim) 
    {
      if(is_ok[i])
	{
	  tree->a_edges[i]->l->v -= inc[i];
	  lk = Lk(tree->a_edges[i],tree);
	  minus_zero[i] = lk+Log_Det(is_ok,tree);
	  tree->a_edges[i]->l->v = ori_bl[i];
	}
    }

  For(i,dim) Update_PMat_At_Given_Edge(tree->a_edges[i],tree);

  /* plus plus  */  
  For(i,dim)
    {
      if(is_ok[i])
	{
	  tree->a_edges[i]->l->v += inc[i];
	  Update_PMat_At_Given_Edge(tree->a_edges[i],tree);

	  For(j,3)
	    if((!tree->a_edges[i]->left->tax) && (tree->a_edges[i]->left->v[j] != tree->a_edges[i]->rght))
	      Recurr_Hessian_Log(tree->a_edges[i]->left,tree->a_edges[i]->left->v[j],1,inc,plus_plus+i*dim,is_ok,tree);
	  
	  For(j,3)
	    if((!tree->a_edges[i]->rght->tax) && (tree->a_edges[i]->rght->v[j] != tree->a_edges[i]->left))
	      Recurr_Hessian_Log(tree->a_edges[i]->rght,tree->a_edges[i]->rght->v[j],1,inc,plus_plus+i*dim,is_ok,tree);

/* 	  For(j,dim)  */
/* 	    if(j != i) */
/* 	      { */
/* 		if(inc[j] > 0.0) */
/* 		  { */
/* 		    tree->a_edges[j]->l->v += inc[j]; */
/* 		    Lk(tree); */
/* 		    plus_plus[i*dim+j]=tree->c_lnL; */
/* 		    tree->a_edges[j]->l->v = ori_bl[j]; */
/* 		  } */
/* 	      } */
		      
	  tree->a_edges[i]->l->v = ori_bl[i];
	  Lk(NULL,tree);
	}
    }

  /* plus minus */  
  For(i,dim)
    {
      if(is_ok[i])
	{
	  tree->a_edges[i]->l->v += inc[i];
	  Update_PMat_At_Given_Edge(tree->a_edges[i],tree);
	  
	  For(j,3)
	    if((!tree->a_edges[i]->left->tax) && (tree->a_edges[i]->left->v[j] != tree->a_edges[i]->rght))
	      Recurr_Hessian_Log(tree->a_edges[i]->left,tree->a_edges[i]->left->v[j],-1,inc,plus_minus+i*dim,is_ok,tree);
	  
	  For(j,3)
	    if((!tree->a_edges[i]->rght->tax) && (tree->a_edges[i]->rght->v[j] != tree->a_edges[i]->left))
	      Recurr_Hessian_Log(tree->a_edges[i]->rght,tree->a_edges[i]->rght->v[j],-1,inc,plus_minus+i*dim,is_ok,tree);
	  
/* 	  For(j,dim)  */
/* 	    if(j != i) */
/* 	      { */
/* 		if(inc[j] > 0.0) */
/* 		  { */
/* 		    tree->a_edges[j]->l->v -= inc[j]; */
/* 		    Lk(tree); */
/* 		    plus_minus[i*dim+j] = tree->c_lnL; */
/* 		    tree->a_edges[j]->l->v = ori_bl[j]; */
/* 		  } */
/* 	      } */

	  tree->a_edges[i]->l->v = ori_bl[i];
	  Lk(NULL,tree);
	}
    }


  /* minus minus */  
  For(i,dim)
    {
      if(is_ok[i])
	{
	  tree->a_edges[i]->l->v -= inc[i];
	  
	  Update_PMat_At_Given_Edge(tree->a_edges[i],tree);
	  
	  For(j,3)
	    if((!tree->a_edges[i]->left->tax) && (tree->a_edges[i]->left->v[j] != tree->a_edges[i]->rght))
	      Recurr_Hessian_Log(tree->a_edges[i]->left,tree->a_edges[i]->left->v[j],-1,inc,minus_minus+i*dim,is_ok,tree);
	  
	  For(j,3)
	    if((!tree->a_edges[i]->rght->tax) && (tree->a_edges[i]->rght->v[j] != tree->a_edges[i]->left))
	      Recurr_Hessian_Log(tree->a_edges[i]->rght,tree->a_edges[i]->rght->v[j],-1,inc,minus_minus+i*dim,is_ok,tree);
	  
/* 	  For(j,dim)  */
/* 	    if(j != i) */
/* 	      { */
/* 		if(inc[j] > 0.0) */
/* 		  { */
/* 		    tree->a_edges[j]->l->v -= inc[j]; */
/* 		    Lk(tree); */
/* 		    minus_minus[i*dim+j] = tree->c_lnL; */
/* 		    tree->a_edges[j]->l->v = ori_bl[j]; */
/* 		  } */
/* 	      } */

	  tree->a_edges[i]->l->v = ori_bl[i];
	  Lk(NULL,tree);
	}
    }

/*   For(i,dim) if(is_ok[i]) inc[i] = POW(tree->a_edges[i]->l->v+inc[i],2)-POW(tree->a_edges[i]->l->v,2); */
  For(i,dim) if(is_ok[i]) inc[i] = LOG(tree->a_edges[i]->l->v+inc[i])-LOG(tree->a_edges[i]->l->v);
/*   For(i,dim) inc[i] = 2.*inc[i]; */
/*   For(i,dim) if(is_ok[i]) inc[i] = SQRT(tree->a_edges[i]->l->v+inc[i])-SQRT(tree->a_edges[i]->l->v); */
  
  For(i,dim)
    {
      if(is_ok[i])
	{
	  hessian[i*dim+i] = (plus_zero[i]-2*zero_zero[i]+minus_zero[i])/(POW(inc[i],2));

	  for(j=i+1;j<dim;j++)
	    {
	      if(is_ok[j])
		{
		  hessian[i*dim+j] = 
		    (plus_plus[i*dim+j]-plus_minus[i*dim+j]-plus_minus[j*dim+i]+minus_minus[i*dim+j])/
		    (4*inc[i]*inc[i]);
		  hessian[j*dim+i] = hessian[i*dim+j];
		}
	    }
	}
    }
        

  For(i,n_ok_edges)
    {
      For(j,n_ok_edges)
	{
	  buff[i*n_ok_edges+j] = -hessian[ok_edges[i]*dim+ok_edges[j]];
	}
    }

  if(!Matinv(buff,n_ok_edges,n_ok_edges,NO))
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");      
    }

  For(i,n_ok_edges)
    {
      For(j,n_ok_edges)
	{
	  hessian[ok_edges[i]*dim+ok_edges[j]] = buff[i*n_ok_edges+j];
	}
    }

  /* Approximate variance for very short branches */
  For(i,dim)
    if(!is_ok[i])
      {
	hessian[i*dim+i] = 1./POW(tree->data->init_len,2);
      }

  if(!Matinv(hessian,dim,dim,NO))
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");      
    }

  For(i,dim*dim) hessian[i] = -1.0*hessian[i];

/*   For(i,dim) */
/*     { */
/*       PhyML_Printf("[%f] ",tree->a_edges[i]->l->v); */
/*       For(j,i+1) */
/* 	{ */
/* 	  PhyML_Printf("%12lf ",hessian[i*dim+j]); */
/* 	} */
/*       PhyML_Printf("\n"); */
/*     } */

/*   Matinv(hessian,dim,dim); */

/*   PhyML_Printf("\n"); */

  For(i,dim)
    {
      PhyML_Printf("[%f] ",tree->a_edges[i]->l->v);
      For(j,i+1)
	{
	  PhyML_Printf("%12lf ",hessian[i*dim+j]);
	}
      PhyML_Printf("\n");
    }
/*   Exit("\n"); */


  Free(ori_bl);
  Free(plus_plus);
  Free(minus_minus);
  Free(plus_zero);
  Free(minus_zero);
  Free(plus_minus);
  Free(zero_zero);
  Free(inc);
  Free(buff);
  Free(ok_edges);
  Free(is_ok);

  return hessian;

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Recurr_Hessian_Log(t_node *a, t_node *d, int plus_minus, phydbl *inc, phydbl *res, int *is_ok, t_tree *tree)
{
  int i;
  phydbl ori_l;

  For(i,3)
    if(a->v[i] == d)
      {
	Update_P_Lk(tree,a->b[i],a);

	ori_l = a->b[i]->l->v;
	if(is_ok[a->b[i]->num])
	  {
	    if(plus_minus > 0) a->b[i]->l->v += inc[a->b[i]->num];
	    else               a->b[i]->l->v -= inc[a->b[i]->num];
	    res[a->b[i]->num]  = Lk(a->b[i],tree);
	    res[a->b[i]->num] += Log_Det(is_ok,tree);
	    a->b[i]->l->v = ori_l;
	    Update_PMat_At_Given_Edge(a->b[i],tree);
	  }
	break;
      }

  if(d->tax) return;
  else 
    For(i,3) 
      if(d->v[i] != a) 
	Recurr_Hessian_Log(d,d->v[i],plus_minus,inc,res,is_ok,tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Log_Det(int *is_ok, t_tree *tree)
{
  int i;
  phydbl ldet;

  ldet = 0.0;
/*   For(i,2*tree->n_otu-3) if(is_ok[i]) ldet += LOG(2.*SQRT(tree->a_edges[i]->l->v)); */
  For(i,2*tree->n_otu-3) if(is_ok[i]) ldet += LOG(tree->a_edges[i]->l->v);
/*   For(i,2*tree->n_otu-3) if(is_ok[i]) ldet -= LOG(2*tree->a_edges[i]->l->v); */
  
  return ldet;

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Normal_Trunc_Mean(phydbl mu, phydbl sd, phydbl min, phydbl max)
{
  phydbl mean;

  mean = mu + sd * 
    (Dnorm((min-mu)/sd,0.,1.)-Dnorm((max-mu)/sd,0.,1.))/
    (Pnorm((max-mu)/sd,0.,1.)-Pnorm((min-mu)/sd,0.,1.));
  return mean;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Constraint_Normal_Trunc_Mean(phydbl wanted_mu, phydbl sd, phydbl min, phydbl max)
{
  int j;
  phydbl dx,f,fmid,xmid,rtb;
  phydbl x1, x2;

  x1 = min;
  x2 = max;

  f    = Normal_Trunc_Mean(x1,sd,min,max) - wanted_mu;
  fmid = Normal_Trunc_Mean(x2,sd,min,max) - wanted_mu;
  
  if(f*fmid >= 0.0)
    {
      PhyML_Printf("\n. Root must be bracketed for bisection!");
      PhyML_Printf("\n. f=%f fmid=%f",f,fmid);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");
    }

  rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);

  For(j,100) 
    {
      xmid=rtb+(dx *= 0.5);
      fmid=Normal_Trunc_Mean(xmid,sd,min,max)-wanted_mu;
      if(fmid <= 0.0) rtb=xmid;
      if(fmid > -1.E-10 && fmid < 1.E-10) return rtb;
    }

  Exit("Too many bisections in RTBIS");
  return(-1.);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


int Matinv(phydbl *x, int n, int m, int verbose)
{

/* x[n*m]  ... m>=n
*/

   int i,j,k;
   int *irow;
   phydbl ee, t,t1,xmax;
   phydbl det;

   ee = 1.0E-10;
   det = 1.0;
   
   irow = (int *)mCalloc(n,sizeof(int));

   For (i,n)
     {
       xmax = 0.;
       for (j=i; j<n; j++)
         if (xmax < FABS(x[j*m+i]))
	   {
	     xmax = FABS(x[j*m+i]);
	     irow[i]=j;
	   }

      det *= xmax;
      if (xmax < ee)
	{
	  Free(irow);
	  if(verbose)
	    {
	      PhyML_Printf("\n== Determinant becomes zero at %3d!\t", i+1);
	      PhyML_Printf("\n== Failed to invert the matrix.");
	    }
	  return(0);
	}
      if (irow[i] != i)
	{
	  For (j,m)
	    {
	      t = x[i*m+j];
	      x[i*m+j] = x[irow[i]*m+j];
	      x[irow[i]*m+j] = t;
	    }
	}
      t = 1./x[i*m+i];
      For (j,n)
	{
	  if (j == i) continue;
	  t1 = t*x[j*m+i];
	  For(k,m)  x[j*m+k] -= t1*x[i*m+k];
	  x[j*m+i] = -t1;
	}
      For(j,m)   x[i*m+j] *= t;
      x[i*m+i] = t;
   }                            /* i  */
   for (i=n-1; i>=0; i--)
     {
       if (irow[i] == i) continue;
       For(j,n)
	 {
	   t = x[j*m+i];
	   x[j*m+i] = x[j*m + irow[i]];
	   x[j*m + irow[i]] = t;
	 }
     }

   Free(irow);
   return (1);

/*   int i, j, k, lower, upper; */
/*   phydbl temp; */
/*   phydbl *a; */
/*   int nsize; */

/*   nsize = n; */
/*   a = x; */
  
/*   /\*Gauss-Jordan reduction -- invert matrix a in place, */
/*          overwriting previous contents of a.  On exit, matrix a */
/*          contains the inverse.*\/ */
/*   lower = 0; */
/*   upper = nsize-1; */
/*   for(i = lower; i <= upper; i++)  */
/*     { */
/*       temp = 1.0 / a[i*n+i]; */
/*       a[i*n+i] = 1.0; */
/*       for (j = lower; j <= upper; j++)  */
/* 	{ */
/* 	  a[i*n+j] *= temp; */
/* 	} */
/*       for (j = lower; j <= upper; j++)  */
/* 	{ */
/* 	  if (j != i)  */
/* 	    { */
/* 	      temp = a[j*n+i]; */
/* 	      a[j*n+i] = 0.0; */
/* 	      for (k = lower; k <= upper; k++)  */
/* 		{ */
/* 		  a[j*n+k] -= temp * a[i*n+k]; */
/* 		}	       */
/* 	    } */
/* 	} */
/*     } */

  return(1);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


int Iter_Matinv(phydbl *x, int n, int m, int verbose)
{
  phydbl *buff;
  int i,iter;
  phydbl scaler;
  int pb;

  buff = (phydbl *)mCalloc(n*m,sizeof(phydbl));

  pb = NO;
  iter   = 0;
  scaler = 1.;
  For(i,n*m) buff[i] = x[i];
  while(!Matinv(buff,n,m,verbose))
    {
      pb = YES;
      For(i,n*m) buff[i] = x[i];
      scaler *= 10.;
      For(i,n*m) buff[i] *= scaler;
      iter++;

      if(iter > 100)
	{
	  PhyML_Printf("\n== Err in file %s at line %d.",__FILE__,__LINE__);
	  return 0;
	}      
    }
  if(pb)  PhyML_Printf("\n== Managed to fix the problem by rescaling the matrix.");
  For(i,n*m) x[i] = buff[i]*scaler;
  Free(buff);
  return 1;
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl *Matrix_Mult(phydbl *A, phydbl *B, int nra, int nca, int nrb, int ncb)
{
  int i,j,k;
  phydbl *C;

  C = (phydbl *)mCalloc(nra*ncb,sizeof(phydbl));

  if(nca != nrb)
    {
      PhyML_Printf("\n. Matrices dimensions don't match.");
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");      
    }
  
  For(i,nra)
    For(j,ncb)
       For(k,nca)
         C[i*ncb+j] += A[i*nca+k] * B[k*ncb+j];
  
  return C;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl *Matrix_Transpose(phydbl *A, int dim)
{
  phydbl *tA,buff;
  int i,j;

  tA = (phydbl *)mCalloc(dim*dim,sizeof(phydbl));

  For(i,dim*dim) tA[i]=A[i];

  For(i,dim) for(j=i+1;j<dim;j++) 
    {
      buff        = tA[i*dim+j];
      tA[i*dim+j] = tA[j*dim+i];
      tA[j*dim+i]  = buff;
    }

  return tA;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Matrix_Det(phydbl *A, int size, int _log)
{
  phydbl *triA;
  int i;
  phydbl det;

  triA = Cholesky_Decomp(A,size);
  det = 0.0;
  For(i,size) det += LOG(triA[i*size+i]);
  Free(triA);
 
  if(_log == NO)
    {
      det = EXP(det);
      return det*det;
    }
  else
    {
      return 2.*det;
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


/* http://en.wikipedia.org/wiki/Multivariate_normal_distribution (Conditional distributions) */
void Normal_Conditional(phydbl *mu, phydbl *cov, phydbl *a, int n, short int *is_1, int n1, phydbl *cond_mu, phydbl *cond_cov)
{
  phydbl *mu1,*mu2;
  phydbl *sig11,*sig12,*sig21,*sig22,*sig12_invsig22,*buff;
  phydbl *ctrd_a;
  phydbl *cond_cov_norder,*cond_mu_norder;
  int    n2;
  int i,j,nr,nc;
  phydbl *buff_mat;

  n2 = n-n1;

  mu1             = (phydbl *)mCalloc(n1,   sizeof(phydbl));
  mu2             = (phydbl *)mCalloc(n2,   sizeof(phydbl));
  sig11           = (phydbl *)mCalloc(n1*n1,sizeof(phydbl));
  sig12           = (phydbl *)mCalloc(n1*n2,sizeof(phydbl));
  sig21           = (phydbl *)mCalloc(n2*n1,sizeof(phydbl));
  sig22           = (phydbl *)mCalloc(n2*n2,sizeof(phydbl));
  ctrd_a          = (phydbl *)mCalloc(n2,   sizeof(phydbl)); 
  cond_cov_norder = (phydbl *)mCalloc(n1*n1,sizeof(phydbl));
  cond_mu_norder  = (phydbl *)mCalloc(n1*n1,sizeof(phydbl));
  buff_mat        = (phydbl *)mCalloc(n2*n2,sizeof(phydbl));

  nr=0;
  For(i,n) { if(!is_1[i]) { ctrd_a[nr] = a[i]-mu[i]; nr++; } }

  nr=0;
  For(i,n) { if( is_1[i]) { mu1[nr] = mu[i]; nr++; } }

  nr=0;
  For(i,n) { if(!is_1[i]) { mu2[nr] = mu[i]; nr++; } }

  nr=0; nc=0;
  For(i,n)
    {
      if(is_1[i])
	{
	  nc = nr;
 	  for(j=i;j<n;j++)
/* 	  nc = 0; */
/* 	  For(j,n) */
	    {
	      if(is_1[j])
		{
		  sig11[nr*n1+nc] = cov[i*n+j];
		  sig11[nc*n1+nr] = cov[i*n+j];
		  nc++;
		}
	    }
	  nr++;
	}
    }


  nr=0; nc=0;
  For(i,n)
    {
      if(is_1[i])
	{
/* 	  nc = nr; */
/*  	  for(j=i;j<n;j++) */
	  nc = 0;
	  For(j,n)
	    {
	      if(!is_1[j])
		{
		  sig12[nr*n2+nc] = cov[i*n+j];
/* 		  sig12[nc*n2+nr] = cov[i*n+j]; */
		  nc++;
		}
	    }
	  nr++;
	}
    }

  nr=0; nc=0;
  For(i,n)
    {
      if(!is_1[i])
	{
/* 	  nc = nr; */
/* 	  for(j=i;j<n;j++) */
	  nc = 0;
	  For(j,n)
	    {
	      if(is_1[j])
		{
		  sig21[nr*n1+nc] = cov[i*n+j];
/* 		  sig21[nc*n1+nr] = cov[i*n+j]; */
		  nc++;
		}
	    }
	  nr++;
	}
    }


  nr=0; nc=0;
  For(i,n)
    {
      if(!is_1[i])
	{
	  nc = nr;
	  for(j=i;j<n;j++)
/* 	  nc = 0; */
/* 	  For(j,n) */
	    {
	      if(!is_1[j])
		{
		  sig22[nr*n2+nc] = cov[i*n+j];
 		  sig22[nc*n2+nr] = cov[i*n+j];
		  nc++;
		}
	    }
	  nr++;
	}
    }

  Iter_Matinv(sig22,n2,n2,NO);

  sig12_invsig22 = Matrix_Mult(sig12,sig22,n1,n2,n2,n2);

  buff = Matrix_Mult(sig12_invsig22,ctrd_a,n1,n2,n2,1);
  For(i,n1) cond_mu_norder[i] = mu1[i]+buff[i];
  Free(buff);

  buff = Matrix_Mult(sig12_invsig22,sig21,n1,n2,n2,n1);
  For(i,n1) For(j,n1) cond_cov_norder[i*n1+j] = sig11[i*n1+j] - buff[i*n1+j];
  Free(buff);

  nr = 0;
  For(i,n) if(is_1[i]) { cond_mu[i] = cond_mu_norder[nr]; nr++; }

  nr = nc = 0;
  For(i,n) 
    {
      if(is_1[i]) 
	{ 
	  nc = 0;
	  For(j,n)
	    {
	      if(is_1[j]) 
		{		  
		  cond_cov[i*n+j] = cond_cov_norder[nr*n1+nc]; 
		  nc++;
		}
	    }
	  nr++;
	}
    }

/*   For(i,n1) */
/*     { */
/*       for(j=i;j<n1;j++) */
/* 	if(FABS(cond_cov_norder[i*n1+j] - cond_cov_norder[j*n1+i]) > 1.E-3) */
/* 	  { */
/* 	    PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__); */
/* 	    Warn_And_Exit(""); */
/* 	  } */
/*     } */


  For(i,n)
    {
      for(j=i+1;j<n;j++)
	if(FABS(cond_cov[i*n+j] - cond_cov[j*n+i]) > 1.E-3)
	  {
	    PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	    Warn_And_Exit("");
	  }
    }

  Free(mu1);
  Free(mu2);
  Free(sig11);
  Free(sig12);
  Free(sig21);
  Free(sig22);
  Free(ctrd_a);
  Free(cond_cov_norder);
  Free(cond_mu_norder);
  Free(sig12_invsig22);
  Free(buff_mat);
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


/* http://en.wikipedia.org/wiki/Multivariate_normal_distribution (Conditional distributions) */
void Normal_Conditional_Unsorted(phydbl *mu, phydbl *cov, phydbl *a, int n, short int *is_1, int n1, phydbl *cond_mu, phydbl *cond_cov)
{
  phydbl *mu1,*mu2;
  phydbl *sig11,*sig12,*sig21,*sig22,*sig12_invsig22,*buff;
  phydbl *ctrd_a;
  int    n2;
  int i,j,nr,nc;

  n2 = n-n1;

  mu1             = (phydbl *)mCalloc(n1,   sizeof(phydbl));
  mu2             = (phydbl *)mCalloc(n2,   sizeof(phydbl));
  sig11           = (phydbl *)mCalloc(n1*n1,sizeof(phydbl));
  sig12           = (phydbl *)mCalloc(n1*n2,sizeof(phydbl));
  sig21           = (phydbl *)mCalloc(n2*n1,sizeof(phydbl));
  sig22           = (phydbl *)mCalloc(n2*n2,sizeof(phydbl));
  ctrd_a          = (phydbl *)mCalloc(n2,   sizeof(phydbl)); 

  nr=0;
  For(i,n) { if(!is_1[i]) { ctrd_a[nr] = a[i]-mu[i]; nr++; } }

  nr=0;
  For(i,n) { if( is_1[i]) { mu1[nr] = mu[i]; nr++; } }

  nr=0;
  For(i,n) { if(!is_1[i]) { mu2[nr] = mu[i]; nr++; } }

  nr=0; nc=0;
  For(i,n)
    {
      if(is_1[i])
	{
	  nc = nr;
 	  for(j=i;j<n;j++)
/* 	  nc = 0; */
/* 	  For(j,n) */
	    {
	      if(is_1[j])
		{
		  sig11[nr*n1+nc] = cov[i*n+j];
		  sig11[nc*n1+nr] = cov[i*n+j];
		  nc++;
		}
	    }
	  nr++;
	}
    }


  nr=0; nc=0;
  For(i,n)
    {
      if(is_1[i])
	{
/* 	  nc = nr; */
/*  	  for(j=i;j<n;j++) */
	  nc = 0;
	  For(j,n)
	    {
	      if(!is_1[j])
		{
		  sig12[nr*n2+nc] = cov[i*n+j];
/* 		  sig12[nc*n2+nr] = cov[i*n+j]; */
		  nc++;
		}
	    }
	  nr++;
	}
    }

  nr=0; nc=0;
  For(i,n)
    {
      if(!is_1[i])
	{
/* 	  nc = nr; */
/* 	  for(j=i;j<n;j++) */
	  nc = 0;
	  For(j,n)
	    {
	      if(is_1[j])
		{
		  sig21[nr*n1+nc] = cov[i*n+j];
/* 		  sig21[nc*n1+nr] = cov[i*n+j]; */
		  nc++;
		}
	    }
	  nr++;
	}
    }


  nr=0; nc=0;
  For(i,n)
    {
      if(!is_1[i])
	{
	  nc = nr;
	  for(j=i;j<n;j++)
/* 	  nc = 0; */
/* 	  For(j,n) */
	    {
	      if(!is_1[j])
		{
		  sig22[nr*n2+nc] = cov[i*n+j];
 		  sig22[nc*n2+nr] = cov[i*n+j];
		  nc++;
		}
	    }
	  nr++;
	}
    }

  if(!Matinv(sig22,n2,n2,NO))
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");      
    }
  sig12_invsig22 = Matrix_Mult(sig12,sig22,n1,n2,n2,n2);

  buff = Matrix_Mult(sig12_invsig22,ctrd_a,n1,n2,n2,1);
  For(i,n1) cond_mu[i] = mu1[i]+buff[i];
  Free(buff);

  buff = Matrix_Mult(sig12_invsig22,sig21,n1,n2,n2,n1);
  For(i,n1) For(j,n1) cond_cov[i*n1+j] = sig11[i*n1+j] - buff[i*n1+j];


  Free(mu1);
  Free(mu2);
  Free(sig11);
  Free(sig12);
  Free(sig21);
  Free(sig22);
  Free(ctrd_a);
  Free(sig12_invsig22);
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


/* http://en.wikipedia.org/wiki/Multivariate_normal_distribution (Conditional distributions) */
void Get_Reg_Coeff(phydbl *mu, phydbl *cov, phydbl *a, int n, short int *is_1, int n1, phydbl *reg_coeff)
{
  phydbl *sig12,*sig22,*sig12_invsig22;
  int    n2;
  int    i,j,nr,nc;

  n2 = n-n1;

  sig12 = (phydbl *)mCalloc(n1*n2,sizeof(phydbl));
  sig22 = (phydbl *)mCalloc(n2*n2,sizeof(phydbl));

  nr=0; nc=0;
  For(i,n)
    {
      if(is_1[i])
	{
	  nc = 0;
	  For(j,n)
	    {
	      if(!is_1[j])
		{
		  sig12[nr*n2+nc] = cov[i*n+j];
		  nc++;
		}
	    }
	  nr++;
	}
    }


  nr=0; nc=0;
  For(i,n)
    {
      if(!is_1[i])
	{
	  nc = nr;
	  for(j=i;j<n;j++)
	    {
	      if(!is_1[j])
		{
		  sig22[nr*n2+nc] = cov[i*n+j];
 		  sig22[nc*n2+nr] = cov[i*n+j];
		  nc++;
		}
	    }
	  nr++;
	}
    }


  if(!Matinv(sig22,n2,n2,NO))
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");      
    }
  sig12_invsig22 = Matrix_Mult(sig12,sig22,n1,n2,n2,n2);


  For(i,n) reg_coeff[i] = 0.0;

/*   nr = 0; */
/*   For(i,n) if(!is_1[i]) { reg_coeff[i] = sig12_invsig22[nr]; nr++; } */

  nc = 0;
  nr = 0;
  For(i,n1) 
    {
      nc = 0;
      For(j,n)
	if(!is_1[j]) 
	  { 
	    reg_coeff[i*n+j] = sig12_invsig22[nr*n2+nc]; 
	    nc++; 
	  }
      nr++;
    }


  if(nc != n2 || nr != n1)
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");
    }


  Free(sig12);
  Free(sig22);
  Free(sig12_invsig22);
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Norm_Trunc_Sd(phydbl mu, phydbl sd, phydbl a, phydbl b)
{
  phydbl pdfa, pdfb;
  phydbl cdfa, cdfb;
  phydbl ctra, ctrb;
  phydbl cond_var;
  phydbl cdfbmcdfa;

  ctra = (a - mu)/sd;
  ctrb = (b - mu)/sd;

  pdfa = Dnorm(ctra,0.0,1.0);
  pdfb = Dnorm(ctrb,0.0,1.0);

  cdfa = Pnorm(ctra,0.0,1.0);
  cdfb = Pnorm(ctrb,0.0,1.0);

  cdfbmcdfa = cdfb - cdfa;

  if(cdfbmcdfa < SMALL) 
    {
      cdfbmcdfa = SMALL;
      PhyML_Printf("\n. mu=%G sd=%G a=%G b=%G",mu,sd,a,b);
      PhyML_Printf("\n. Numerical precision issue detected.");
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
    }
	    
  cond_var = sd*sd*(1. + (ctra*pdfa - ctrb*pdfb)/cdfbmcdfa - POW((pdfa - pdfb)/cdfbmcdfa,2));

  return SQRT(cond_var);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Norm_Trunc_Mean(phydbl mu, phydbl sd, phydbl a, phydbl b)
{
  phydbl pdfa, pdfb;
  phydbl cdfa, cdfb;
  phydbl ctra, ctrb;
  phydbl cond_mu;
  phydbl cdfbmcdfa;

  ctra = (a - mu)/sd;
  ctrb = (b - mu)/sd;

  pdfa = Dnorm(ctra,0.0,1.0);
  pdfb = Dnorm(ctrb,0.0,1.0);

  cdfa = Pnorm(ctra,0.0,1.0);
  cdfb = Pnorm(ctrb,0.0,1.0);
  
  cdfbmcdfa = cdfb - cdfa;

  if(cdfbmcdfa < SMALL)
    {
      cdfbmcdfa = SMALL;
      PhyML_Printf("\n. mu=%G sd=%G a=%G b=%G",mu,sd,a,b);
      PhyML_Printf("\n. Numerical precision issue detected.");
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
    }
  
  cond_mu = mu + sd*(pdfa - pdfb)/cdfbmcdfa;

  return cond_mu;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


int Norm_Trunc_Mean_Sd(phydbl mu, phydbl sd, phydbl a, phydbl b, phydbl *trunc_mu, phydbl *trunc_sd)
{

  phydbl pdfa, pdfb;
  phydbl cdfa, cdfb;
  phydbl ctra, ctrb;
  phydbl cdfbmcdfa;

  ctra = (a - mu)/sd;
  ctrb = (b - mu)/sd;

  pdfa = Dnorm(ctra,0.0,1.0);
  pdfb = Dnorm(ctrb,0.0,1.0);

  cdfa = Pnorm(ctra,0.0,1.0);
  cdfb = Pnorm(ctrb,0.0,1.0);
  
  cdfbmcdfa = cdfb - cdfa;

  if(cdfbmcdfa < SMALL)
    {
      cdfbmcdfa = SMALL;
      PhyML_Printf("\n. mu=%G sd=%G a=%G b=%G",mu,sd,a,b);
      PhyML_Printf("\n. cdfa=%G cdfb=%G",cdfa,cdfb);
      PhyML_Printf("\n. Numerical precision issue detected.");
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      return 0;
    }
  
  *trunc_mu = mu + sd*(pdfa - pdfb)/cdfbmcdfa;
  *trunc_sd = sd*sd*(1. + (ctra*pdfa - ctrb*pdfb)/cdfbmcdfa - POW((pdfa - pdfb)/cdfbmcdfa,2));
  *trunc_sd = SQRT(*trunc_sd);
  return 1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void VarCov_Approx_Likelihood(t_tree *tree)
{
  int i,j;
  phydbl *cov;
  phydbl *mean;
  int dim;
  int iter;
  phydbl cur_mean,new_mean,diff_mean,max_diff_mean;
  phydbl cur_cov,new_cov,diff_cov,max_diff_cov;
  FILE *fp;

  
  cov = tree->rates->cov_l;
  mean = tree->rates->mean_l;
  dim = 2*tree->n_otu-3;

  fp = fopen("covariance","w");
  fprintf(fp,"\n");
  fprintf(fp,"Run\t");
  fprintf(fp,"lnL\t");
  For(i,dim) fprintf(fp,"Edge%d[%f]\t",i,tree->rates->u_ml_l[i]);

  
  For(i,dim)     mean[i] = .0;
  For(i,dim*dim) cov[i]  = .0;

  MCMC_Randomize_Branch_Lengths(tree);
  
  /* For(i,2*tree->n_otu-3) tree->a_edges[i]->l->v *= Rgamma(5.,1./5.); */
  
  Set_Both_Sides(YES,tree);
  Lk(NULL,tree);

  iter = 0;
  do
    {
      /* tree->both_sides = YES; */
      /* Lk(tree); */
      MCMC_Br_Lens(tree);
      /* MCMC_Scale_Br_Lens(tree); */


      max_diff_mean = 0.0;
      For(i,dim)
	{
	  cur_mean = mean[i];

	  mean[i] *= (phydbl)iter;
	  mean[i] += tree->a_edges[i]->l->v;
	  mean[i] /= (phydbl)(iter+1);

	  new_mean = mean[i];	  
	  diff_mean = MAX(cur_mean,new_mean)/MIN(cur_mean,new_mean);
	  if(diff_mean > max_diff_mean) max_diff_mean = diff_mean;
	  /* printf("\n. %d diff_mean = %f %f %f %f",i,diff_mean,cur_mean,new_mean,tree->a_edges[i]->l->v); */
	}

      max_diff_cov = 0.0;
      For(i,dim)
	{
	  For(j,dim)
	    {
	      cur_cov = cov[i*dim+j];

	      cov[i*dim+j] *= (phydbl)iter;
	      cov[i*dim+j] += tree->a_edges[i]->l->v * tree->a_edges[j]->l->v;
	      cov[i*dim+j] /= (phydbl)(iter+1);

	      new_cov = cov[i*dim+j];
	      diff_cov = MAX(cur_cov,new_cov)/MIN(cur_cov,new_cov);
	      if(diff_cov > max_diff_cov) max_diff_cov = diff_cov;
	    }
	}
      iter++;
      
      /* if(!(iter%10)) */
      /* printf("\n. iter=%d max_diff_mean=%f max_diff_cov=%f",iter,max_diff_mean,max_diff_cov); */

      /* if(iter && max_diff_mean < 1.01 && max_diff_cov < 1.01) break; */
      
      if(!(iter%20))
	{
	  fprintf(fp,"\n");
	  fprintf(fp,"%d\t",iter);
	  fprintf(fp,"%f\t",tree->c_lnL);
 	  For(i,dim) fprintf(fp,"%f\t",tree->a_edges[i]->l->v);
	  fflush(NULL);
	}

    }while(iter < 5000);


  For(i,dim)
    {
      For(j,dim)
	{
	  cov[i*dim+j] = cov[i*dim+j] - mean[i]*mean[j];
	  if(i == j && cov[i*dim+j] < MIN_VAR_BL) cov[i*dim+j] = MIN_VAR_BL;
	}
    }

  fclose(fp);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


/* Order statistic. x_is are uniformily distributed in [min,max] */
phydbl Dorder_Unif(phydbl x, int r, int n, phydbl min, phydbl max)
{
  phydbl cons;
  phydbl Fx;
  phydbl dens;

  if(x < min || x > max || min > max)
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");
    }

  cons = LnGamma(n+1) - LnGamma(r) - LnGamma(n-r+1);
  cons = EXP(cons);
  cons = ROUND(cons);

  Fx = (x-min)/(max-min);
  
  dens = cons * pow(Fx,r-1) * pow(1.-Fx,n-r) * (1./(max-min));

  /* printf("\n. x=%f r=%d n=%d min=%f max=%f dens=%f",x,r,n,min,max,dens); */
  /* Exit("\n"); */

  return(dens);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Covariance(phydbl *x, phydbl *y, int n)
{
  int i;
  phydbl mean_x,mean_y,mean_xy;

  mean_x = .0;
  For(i,n) mean_x += x[i];
  mean_x /= (phydbl)n;

  mean_y = .0;
  For(i,n) mean_y += y[i];
  mean_y /= (phydbl)n;

  mean_xy = .0;
  For(i,n) mean_xy += x[i]*y[i];
  mean_xy /= (phydbl)n;
  
  return (mean_xy - mean_x*mean_y);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* Sample X from a multivariate normal with mean mu and covariance cov, within
   the interval [min,max], under the linear constraint X.lambda=k 
*/
   
phydbl *Rnorm_Multid_Trunc_Constraint(phydbl *mu, phydbl *cov, phydbl *min, phydbl *max, phydbl *lambda, phydbl k, phydbl *res, int len)
{

  phydbl *loc_res;
  int i,j,cond,iter;
  phydbl *x;
  phydbl cond_mean,cond_var;
  phydbl cov_zic,cov_zii,cov_zcc;
  phydbl mean_zi, mean_zc;
  phydbl alpha;
  phydbl sum;
  int err;
  phydbl zi;


  cond    = 0;

  loc_res = NULL;
  if(!res) 
    {
      loc_res = (phydbl *)mCalloc(len,sizeof(phydbl));
      x = loc_res;
    }
  else x = res;
  


  /* zi = x[i] . lambda[i] */

  iter = 0;
  do
    {
      sum = 0.0;
      For(i,len)
	{      
	  if(i != cond)
	    {
	      cov_zic = lambda[i]    * lambda[cond] * cov[i*len+cond];
	      cov_zii = lambda[i]    * lambda[i]    * cov[i*len+i];
	      cov_zcc = lambda[cond] * lambda[cond] * cov[cond*len+cond];
	      
	      mean_zi = lambda[i];
	      mean_zc = lambda[cond];
	      
	      /* alpha = k - \sum_{j != cond, j !=i} z_j */
	      alpha = k;
	      For(j,len) if(j != cond && j != i) alpha -= lambda[j] * x[j];
	      
	      cond_mean = mean_zi + (cov_zii + cov_zic) / (cov_zii + 2.*cov_zic + cov_zcc) * (alpha - mean_zi - mean_zc);
	      cond_var  = cov_zii - POW(cov_zii + cov_zic,2)/(cov_zii + 2.*cov_zic + cov_zcc);
	      
	      if(lambda[i]*min[i] > alpha - lambda[cond]*min[i])
		{
		  PhyML_Printf("\n. Cannot satisfy the constraint.\n");
		  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		  Exit("\n");
		}

	      err = NO;
	      zi = Rnorm_Trunc(cond_mean,SQRT(cond_var),
			       MAX(lambda[i]*min[i],alpha-lambda[cond]*max[cond]),
			       MIN(lambda[i]*max[i],alpha-lambda[cond]*min[cond]),&err);
	      if(err == YES)
		{
		  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		  Exit("\n");
		}
	      sum += zi;
	      x[i] = zi / lambda[i];
	    }
	}
      
      x[cond] = (k - sum)/lambda[cond];
    
    }while(iter++ < 10);

  return(loc_res);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////




void Integrated_Brownian_Bridge_Moments(phydbl x_beg, phydbl x_end, 
					phydbl y_beg, phydbl y_end, 
					phydbl brownian_var, phydbl *mean, phydbl *var)
{
  /* phydbl *y; */
  /* phydbl *y_mean; */
  /* int n_rep; */
  int n_breaks;
  int i;
  /* int j; */
  /* phydbl traj_mean, traj_sd; */
  /* phydbl x_prev, x_curr; */
  phydbl x;
  phydbl x_step;
  phydbl sum;
  /* phydbl sumsum; */
  phydbl scaled_var;

  scaled_var = brownian_var/FABS(x_end - x_beg);

  n_breaks = 100;


  /* n_rep    = 500;   */

  /* x_step   = (x_end - x_beg)/(n_breaks+1); */

  /* y      = (phydbl *)mCalloc(n_breaks+2,sizeof(phydbl)); */
  /* y_mean = (phydbl *)mCalloc(n_rep,sizeof(phydbl)); */

  /* y[0] = y_beg; */
  /* y[n_breaks+1] = y_end; */

  /* For(i,n_rep) */
  /*   { */
  /*     for(j=1;j<n_breaks+1;j++) */
  /* 	{ */
  /* 	  x_prev = x_beg + (j-1)*x_step; */
  /* 	  x_curr = x_prev + x_step; */

  /* 	  traj_mean = y[j-1] + (y_end - y[j-1])*(x_curr - x_prev)/(x_end - x_prev); */
  /* 	  traj_sd   = SQRT(scaled_var*(x_curr - x_prev)*(x_end - x_curr)/(x_end - x_prev)); */

  /* 	  if(isnan(traj_mean) || isnan(traj_sd)) */
  /* 	    { */
  /* 	      PhyML_Printf("\n. traj_mean=%f traj_sd=%f x_end=%f x_prev=%f x_step=%f [%f %f %f %f %f %f %f] j=%d n_breaks=%d", */
  /* 			   traj_mean,traj_sd,x_end,x_prev,x_step, */
  /* 			   y[j-1],y_end,y[j-1],x_curr,x_prev,x_end,x_prev,j,n_breaks); */
  /* 	      Exit("\n"); */
  /* 	    } */

  /* 	  y[j] = Rnorm(traj_mean,traj_sd); */

  /* 	  if(isnan(y[j]) || isinf(y[j])) */
  /* 	    { */
  /* 	      printf("\n. mean=%f sd=%f %f j=%d y[j]=%f",traj_sd,traj_mean,Rnorm(traj_mean,traj_sd),j,y[j]); */
  /* 	      Exit("\n"); */
  /* 	    } */

  /* 	} */
      
  /*     sum = 0.0; */
  /*     For(j,n_breaks+2) sum += FABS(y[j]); */
  /*     y_mean[i] = sum/(n_breaks+2); */
  /*   } */

  /* sum = sumsum = 0.0; */
  /* For(i,n_rep) */
  /*   { */
  /*     sum += y_mean[i]; */
  /*     sumsum += y_mean[i] * y_mean[i]; */
  /*   } */

  /* *mean = sum/n_rep; */
  /* *var = sumsum/n_rep - (*mean) * (*mean); */

  /* if(isnan(*mean) || isnan(*var)) */
  /*   { */
  /*     PhyML_Printf("\n. sum=%f sumsum=%f n_rep=%d",sum,sumsum,n_rep); */
  /*     Exit("\n"); */
  /*   } */

  /* Free(y); */
  /* Free(y_mean); */

  /* /\* printf("\n. [%f %f]",*mean,*var); *\/ */

  phydbl mux,six;

  x_step = (x_end - x_beg)/(n_breaks+1);
  sum = y_beg;
  for(i=1;i<n_breaks+1;i++)
    {
      x = x_beg + i*x_step;

      mux = y_beg + (y_end - y_beg)*(x - x_beg)/(x_end - x_beg);
      six = SQRT(scaled_var*(x - x_beg)*(x_end - x)/(x_end - x_beg));

      sum += 
	(2.*six)/SQRT(2.*PI)*EXP(-POW(mux,2)/(2.*POW(six,2))) + 
	2.*mux*Pnorm(mux/six,.0,1.) - mux;
    }
  sum += y_end;

  (*mean) = sum / (n_breaks+2.);
  (*var)  = (1./12.)*scaled_var*(x_end - x_beg);

  /* printf(" [%f %f] -- x_beg=%f x_end=%f y_beg=%f y_end=%f sd=%f", */
  /* 	 (*mean),(*var),x_beg,x_end,y_beg,y_end,brownian_var); */
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


/* Let X'(t) = A + (B-A)t/T + X(t) and X(t) = W(t) + t/T * W(T),
i.e., X(t) is a Brownian bridge starting at 0 at t=0 and stopping
at 0 at t=T. X'(t) starts at X'(t)=A at t=0 and stops at X'(t)=B at
t=T. This function calculates the mean and variance of 
Z(T) = 1/T \int_0^T exp(X'(t)) dt. It uses a 13th order approximation
to exp(X) = 1 + X + (1/2!)X^2 + ... (1/10!)X^10
*/

void Integrated_Geometric_Brownian_Bridge_Moments(phydbl T, phydbl A, phydbl B, phydbl u, phydbl *mean, phydbl *var)
{
  Integrated_Geometric_Brownian_Bridge_Mean(T,A,B,u,mean);
  Integrated_Geometric_Brownian_Bridge_Var(T,A,B,u,var);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/*
with(CodeGeneration); C(exp(A)*(int(sum(((B-A)*s/t+(1/2)*u*(t-s)*s/t)^k/factorial(k), k = 0 .. 13), s = 0 .. t))/t, optimize)
*/

void Integrated_Geometric_Brownian_Bridge_Mean(phydbl T, phydbl A, phydbl B, phydbl u, phydbl *mean)
{
phydbl t1 = exp(A);
phydbl t2 = A * A;
phydbl t3 = t2 * t2;
phydbl t4 = t3 * t3;
phydbl t5 = t4 * t2;
phydbl t6 = B * B;
phydbl t7 = t6 * B;
phydbl t10 = t4 * A;
phydbl t11 = t6 * t6;
phydbl t14 = t2 * A;
phydbl t15 = t4 * t14;
phydbl t18 = t11 * B;
phydbl t21 = t4 * t3;
phydbl t24 = t11 * t11;
phydbl t25 = t24 * t7;
phydbl t26 = t25 * A;
phydbl t28 = t14 * t7;
phydbl t30 = t6 * A;
phydbl t32 = B * t2;
phydbl t34 = t7 * A;
phydbl t36 = B * t14;
phydbl t38 = t7 * t2;
phydbl t40 = t11 * t4;
phydbl t42 = B * t10;
phydbl t44 = t7 * t10;
phydbl t46 = t6 * t5;
phydbl t48 = B * t15;
phydbl t50 = -0.46994831769600e14 * t5 * t7 + 0.117487079424000e15 * t10 * t11 + 0.12816772300800e14 * t15 * t6 - 0.211476742963200e15 * t4 * t18 - 0.2136128716800e13 * t21 * B + 0.27605355724800e14 * t26 + 0.56844948508508160000e20 * t28 + 0.1790615878018007040000e22 * t30 - 0.1790615878018007040000e22 * t32 + 0.477497567471468544000e21 * t34 + 0.477497567471468544000e21 * t36 - 0.198957319779778560000e21 * t38 - 0.1138720923648000e16 * t40 + 0.3588696244224000e16 * t42 + 0.506098188288000e15 * t44 - 0.151829456486400e15 * t46 + 0.27605355724800e14 * t48;
phydbl t51 = u * u;
phydbl t52 = t51 * t51;
phydbl t53 = t52 * t52;
phydbl t54 = t53 * t52;
phydbl t55 = T * T;
phydbl t56 = t55 * t55;
phydbl t57 = t56 * t56;
phydbl t58 = t57 * t56;
phydbl t61 = t3 * t14;
phydbl t62 = t11 * t61;
phydbl t64 = t53 * t51;
phydbl t65 = t57 * t55;
phydbl t66 = t64 * t65;
phydbl t68 = t24 * B;
phydbl t69 = t68 * A;
phydbl t70 = t51 * u;
phydbl t71 = t55 * T;
phydbl t72 = t70 * t71;
phydbl t75 = t24 * t6;
phydbl t76 = t75 * A;
phydbl t77 = t51 * t55;
phydbl t80 = u * T;
phydbl t83 = t24 * t3;
phydbl t86 = t3 * t2;
phydbl t87 = B * t86;
phydbl t88 = t52 * t51;
phydbl t89 = t56 * t55;
phydbl t90 = t88 * t89;
phydbl t93 = t11 * t86;
phydbl t96 = t7 * t86;
phydbl t97 = t52 * t56;
phydbl t100 = t18 * t86;
phydbl t103 = t6 * t86;
phydbl t104 = t52 * u;
phydbl t105 = t56 * T;
phydbl t106 = t104 * t105;
phydbl t109 = t4 * t6;
phydbl t112 = t5 * B;
phydbl t117 = B * t61;
phydbl t122 = t11 * t6;
phydbl t123 = t122 * A;
phydbl t126 = -0.108e3 * t54 * t58 + 0.9868914671616000e16 * t62 - 0.993600e6 * t66 + 0.86387558400e11 * t69 * t72 + 0.293717698560e12 * t76 * t77 + 0.854451486720e12 * t26 * t80 - 0.35246123827200e14 * t83 * t80 - 0.795674880e9 * t87 * t90 - 0.1814138726400e13 * t93 * t72 - 0.201570969600e12 * t96 * t97 - 0.12336143339520e14 * t100 * t77 - 0.15913497600e11 * t103 * t106 - 0.388744012800e12 * t109 * t72 - 0.293717698560e12 * t112 * t77 + 0.86387558400e11 * t42 * t72 + 0.4546713600e10 * t117 * t106 + 0.854451486720e12 * t48 * t80 + 0.103520083968000e15 * t72 * t123;
phydbl t128 = t6 * t3;
phydbl t131 = t18 * t14;
phydbl t134 = t11 * t3;
phydbl t137 = t18 * A;
phydbl t140 = t11 * t7;
phydbl t141 = t140 * A;
phydbl t146 = B * A;
phydbl t151 = t24 * t14;
phydbl t154 = t24 * t2;
phydbl t159 = t24 * A;
phydbl t162 = t140 * t3;
phydbl t165 = t6 * t2;
phydbl t170 = t6 * t14;
phydbl t173 = B * t3;
phydbl t176 = -0.51760041984000e14 * t97 * t128 + 0.2898562351104000e16 * t77 * t131 - 0.3623202938880000e16 * t77 * t134 + 0.20704016793600e14 * t97 * t137 + 0.414080335872000e15 * t77 * t141 - 0.177640464089088000e18 * t32 * t72 + 0.16149133099008000e17 * t146 * t97 + 0.1184269760593920000e19 * t36 * t77 + 0.162674417664000e15 * t151 * t80 - 0.16267441766400e14 * t154 * t77 - 0.10658427845345280000e20 * t38 * t80 + 0.1016715110400e13 * t159 * t72 - 0.325348835328000e15 * t162 * t80 - 0.1776404640890880000e19 * t165 * t77 + 0.177640464089088000e18 * t30 * t72 + 0.10658427845345280000e20 * t170 * t80 - 0.5329213922672640000e19 * t173 * t80;
phydbl t177 = t11 * A;
phydbl t184 = t3 * A;
phydbl t185 = t122 * t184;
phydbl t188 = t140 * t14;
phydbl t191 = t140 * t2;
phydbl t196 = t6 * t184;
phydbl t199 = t52 * t70;
phydbl t200 = t56 * t71;
phydbl t201 = t199 * t200;
phydbl t208 = B * t184;
phydbl t213 = t7 * t184;
phydbl t224 = 0.5329213922672640000e19 * t177 * t80 + 0.1184269760593920000e19 * t34 * t77 + 0.239227084800e12 * t141 * t97 + 0.455488369459200e15 * t185 * t80 + 0.43379844710400e14 * t188 * t77 - 0.4066860441600e13 * t191 * t72 + 0.116290944000e12 * t90 * t170 + 0.7117005772800e13 * t97 * t196 - 0.9180864000e10 * t201 * t165 + 0.6120576000e10 * t201 * t34 + 0.40668604416000e14 * t77 * t159 + 0.828988832415744000e18 * t208 * t80 - 0.4710163820544000e16 * t32 * t97 + 0.75914728243200e14 * t72 * t213 + 0.4710163820544000e16 * t30 * t97 + 0.75914728243200e14 * t72 * t131 - 0.1046618496000e13 * t106 * t128 + 0.2372335257600e13 * t97 * t123;
phydbl t227 = t7 * t3;
phydbl t234 = t122 * t3;
phydbl t237 = t122 * t14;
phydbl t242 = t11 * t2;
phydbl t245 = t7 * t61;
phydbl t260 = t18 * t184;
phydbl t267 = -0.11861676288000e14 * t97 * t227 - 0.116290944000e12 * t90 * t38 + 0.1395491328000e13 * t106 * t28 - 0.75914728243200e14 * t234 * t77 + 0.9489341030400e13 * t237 * t72 - 0.414494416207872000e18 * t38 * t77 - 0.2072472081039360000e19 * t242 * t80 + 0.1518294564864000e16 * t80 * t245 + 0.10844961177600e14 * t72 * t141 + 0.418647398400e12 * t106 * t137 + 0.414494416207872000e18 * t170 * t77 - 0.94893410304000e14 * t72 * t134 - 0.2657015488512000e16 * t80 * t93 + 0.2763296108052480000e19 * t28 * t80 + 0.3188418586214400e16 * t80 * t260 + 0.37681310564352000e17 * t34 * t72 + 0.126524547072000e15 * t80 * t69;
phydbl t268 = t122 * t2;
phydbl t285 = t24 * t11;
phydbl t297 = B * t52;
phydbl t300 = t7 * t51;
phydbl t303 = t6 * u;
phydbl t309 = -0.837294796800e12 * t268 * t97 + 0.46516377600e11 * t123 * t106 + 0.1518294564864000e16 * t80 * t188 - 0.569360461824000e15 * t80 * t154 - 0.2657015488512000e16 * t80 * t234 + 0.207247208103936000e18 * t177 * t77 + 0.828988832415744000e18 * t137 * t80 - 0.2072472081039360000e19 * t128 * t80 + 0.2136128716800e13 * t285 * A - 0.119374391867867136000e21 * t11 - 0.2387487837357342720000e22 * t6 + 0.596871959339335680000e21 * t14 - 0.296067440148480000e18 * t2 * t70 * t71 + 0.1776404640890880000e19 * t14 * t51 * t55 - 0.29606744014848000e17 * t297 * t56 - 0.1776404640890880000e19 * t300 * t55 - 0.179061587801800704000e21 * t303 * T + 0.596871959339335680000e21 * A * u * T;
phydbl t311 = B * u;
phydbl t332 = B * t88;
phydbl t341 = t7 * t52;
phydbl t347 = t6 * t104;
phydbl t362 = -0.596871959339335680000e21 * t311 * T + 0.29606744014848000e17 * A * t52 * t56 - 0.7105618563563520000e19 * t3 * u * T + 0.41449441620787200e17 * t184 * t51 * t55 - 0.181160146944000e15 * t2 * t104 * t105 + 0.1570054606848000e16 * t14 * t52 * t56 - 0.138164805402624000e18 * t86 * u * T - 0.12940010496000e14 * t332 * t89 - 0.138164805402624000e18 * t122 * u * T - 0.9420327641088000e16 * t11 * t70 * t71 - 0.1570054606848000e16 * t341 * t56 - 0.41449441620787200e17 * t18 * t51 * t55 - 0.181160146944000e15 * t347 * t105 + 0.1065842784534528000e19 * t184 * u * T - 0.8074566549504000e16 * t2 * t52 * t56 + 0.672880545792000e15 * A * t104 * t105 - 0.296067440148480000e18 * t3 * t51 * t55;
phydbl t366 = B * t104;
phydbl t372 = t6 * t52;
phydbl t375 = t7 * t70;
phydbl t381 = t6 * t51;
phydbl t388 = t4 * B;
phydbl t391 = t61 * t6;
phydbl t394 = t4 * t7;
phydbl t399 = t10 * t6;
phydbl t406 = 0.59213488029696000e17 * t14 * t70 * t71 - 0.672880545792000e15 * t366 * t105 - 0.296067440148480000e18 * t11 * t51 * t55 - 0.8074566549504000e16 * t372 * t56 - 0.59213488029696000e17 * t375 * t71 - 0.1065842784534528000e19 * t18 * u * T - 0.8526742276276224000e19 * t381 * t55 - 0.7162463512072028160000e22 * B + 0.7162463512072028160000e22 * A - 0.596871959339335680000e21 * t7 - 0.2387487837357342720000e22 * t2 - 0.21596889600e11 * t388 * t97 + 0.86387558400e11 * t391 * t97 - 0.4405765478400e13 * t394 * t77 - 0.35246123827200e14 * t40 * t80 + 0.1468588492800e13 * t399 * t77 + 0.15664943923200e14 * t44 * t80 + 0.1036650700800e13 * t245 * t72;
phydbl t412 = t18 * t61;
phydbl t429 = t184 * t70;
phydbl t449 = t2 * t88;
phydbl t453 = t53 * u;
phydbl t455 = t57 * T;
phydbl t467 = t3 * t52;
phydbl t475 = 0.8811530956800e13 * t62 * t77 + 0.56393798123520e14 * t412 * t80 - 0.4699483176960e13 * t46 * t80 + 0.113667840e9 * t184 * t199 * t200 * B + 0.302356454400e12 * t184 * t52 * t56 * t11 + 0.31826995200e11 * t184 * t104 * t105 * t7 + 0.2176966471680e13 * t429 * t71 * t18 + 0.2387024640e10 * t184 * t88 * t89 * t6 - 0.70200e5 * t2 * t64 * t65 * B - 0.284169600e9 * t2 * t199 * t200 * t11 - 0.25833600e8 * t2 * t53 * t57 * t7 - 0.2387024640e10 * t449 * t89 * t18 - 0.1684800e7 * t2 * t453 * t455 * t6 - 0.284169600e9 * t3 * t199 * t200 * t6 + 0.1123200e7 * t14 * t453 * t455 * B - 0.302356454400e12 * t467 * t56 * t18 - 0.12916800e8 * t3 * t53 * t57 * B;
phydbl t496 = t14 * t104;
phydbl t500 = t53 * t70;
phydbl t502 = t57 * t71;
phydbl t510 = A * t453;
phydbl t514 = A * t199;
phydbl t522 = t105 * t3;
phydbl t527 = t184 * t11;
phydbl t536 = t11 * t14;
phydbl t539 = 0.25833600e8 * t14 * t53 * t57 * t6 - 0.3978374400e10 * t3 * t88 * t89 * t7 - 0.39783744000e11 * t3 * t104 * t105 * t11 + 0.378892800e9 * t14 * t199 * t200 * t7 + 0.3978374400e10 * t14 * t88 * t89 * t11 + 0.31826995200e11 * t496 * t105 * t18 + 0.2808e4 * A * t500 * t502 * B + 0.12916800e8 * A * t53 * t57 * t11 + 0.1123200e7 * t510 * t455 * t7 + 0.113667840e9 * t514 * t200 * t18 + 0.70200e5 * A * t64 * t65 * t6 - 0.3235002624000e13 * t366 * t522 - 0.162674417664000e15 * t77 * t191 + 0.14234011545600e14 * t527 * t72 - 0.2093236992000e13 * t134 * t97 + 0.325348835328000e15 * t62 * t80 - 0.75914728243200e14 * t93 * t77 + 0.232581888000e12 * t536 * t106;
phydbl t547 = t7 * u;
phydbl t548 = T * t4;
phydbl t561 = B * t199;
phydbl t562 = t200 * t2;
phydbl t565 = B * t51;
phydbl t566 = t55 * t61;
phydbl t571 = t18 * t2;
phydbl t582 = -0.18361728000e11 * t242 * t90 + 0.30145048451481600e17 * t208 * t77 + 0.918086400e9 * t177 * t201 - 0.162674417664000e15 * t547 * t548 + 0.966187450368000e15 * t36 * t97 - 0.103520083968000e15 * t32 * t106 - 0.110531844322099200e18 * t87 * t80 + 0.6901338931200e13 * t146 * t90 - 0.12560436854784000e17 * t38 * t72 - 0.31715712000e11 * t561 * t562 + 0.414080335872000e15 * t565 * t566 + 0.100483494838272000e18 * t28 * t77 - 0.331595532966297600e18 * t571 * t80 - 0.75362621128704000e17 * t242 * t77 + 0.966187450368000e15 * t34 * t97 - 0.1449281175552000e16 * t165 * t97 + 0.331595532966297600e18 * t196 * t80;
phydbl t601 = t56 * t184;
phydbl t604 = B * t53;
phydbl t605 = t57 * A;
phydbl t610 = t7 * t104;
phydbl t615 = t71 * t86;
phydbl t621 = 0.103520083968000e15 * t30 * t106 + 0.12560436854784000e17 * t170 * t72 - 0.75362621128704000e17 * t128 * t77 - 0.6280218427392000e16 * t173 * t72 + 0.30145048451481600e17 * t137 * t77 + 0.6280218427392000e16 * t177 * t72 + 0.552659221610496000e18 * t536 * t80 - 0.552659221610496000e18 * t227 * t80 + 0.110531844322099200e18 * t123 * t80 + 0.20704016793600e14 * t297 * t601 + 0.1669248000e10 * t604 * t605 + 0.1674589593600e13 * t341 * t601 - 0.232581888000e12 * t610 * t522 + 0.43379844710400e14 * t300 * t566 - 0.9489341030400e13 * t375 * t615 - 0.3947565868646400e16 * t68 + 0.3947565868646400e16 * t10 - 0.1345761091584000e16 * t106;
phydbl t643 = B * t70;
phydbl t649 = t6 * t70;
phydbl t660 = t53 * t57;
phydbl t665 = 0.39791463955955712000e20 * t14 * u * T + 0.1065842784534528000e19 * A * t70 * t71 - 0.8526742276276224000e19 * t2 * t51 * t55 - 0.179061587801800704000e21 * t2 * u * T - 0.29843597966966784000e20 * t565 * t55 + 0.29843597966966784000e20 * A * t51 * t55 - 0.39791463955955712000e20 * t547 * T - 0.1065842784534528000e19 * t643 * t71 - 0.7105618563563520000e19 * t11 * u * T - 0.296067440148480000e18 * t649 * t71 - 0.2842247425425408000e19 * t86 - 0.2842247425425408000e19 * t122 + 0.11861676288000e14 * t97 * t536 - 0.1046618496000e13 * t106 * t242 - 0.56521965846528000e17 * t165 * t72 + 0.459043200e9 * t660 * t30 - 0.37957364121600e14 * t72 * t103;
phydbl t696 = t68 * t2;
phydbl t703 = 0.6470005248000e13 * t106 * t170 - 0.570882816000e12 * t90 * t165 - 0.207247208103936000e18 * t173 * t77 + 0.37681310564352000e17 * t36 * t72 + 0.362320293888000e15 * t146 * t106 + 0.310560251904000e15 * t72 * t196 + 0.1345761091584000e16 * t80 * t159 + 0.358123175603601408000e21 * t146 * t80 + 0.10844961177600e14 * t76 * t80 + 0.119374391867867136000e21 * t30 * t80 + 0.17053484552552448000e20 * t146 * t77 - 0.119374391867867136000e21 * t32 * t80 + 0.380588544000e12 * t90 * t34 + 0.2898562351104000e16 * t77 * t213 + 0.28422474254254080000e20 * t36 * t80 - 0.54224805888000e14 * t696 * t80 + 0.592134880296960000e18 * t146 * t72 + 0.69013389312000e14 * t97 * t28;
phydbl t740 = -0.14324927024144056320000e23 - 0.5329213922672640000e19 * t32 * t77 - 0.6470005248000e13 * t106 * t38 - 0.517600419840000e15 * t72 * t227 + 0.3614987059200e13 * t69 * t77 - 0.42633711381381120000e20 * t165 * t80 + 0.5329213922672640000e19 * t30 * t77 + 0.28422474254254080000e20 * t34 * t80 + 0.355280928178176000e18 * t61 - 0.39475658686464000e17 * t24 - 0.39475658686464000e17 * t4 - 0.1404e4 * t6 * t500 * t502 - 0.2583360e7 * t18 * t53 * t57 - 0.280800e6 * t11 * t453 * t455 - 0.23400e5 * t7 * t64 * t65 - 0.54e2 * B * t54 * t58 - 0.21859200e8 * t11 * t53 * t57;
phydbl t765 = t453 * t455;
phydbl t772 = t500 * t502;
phydbl t783 = -0.1987200e7 * t7 * t453 * t455 - 0.129600e6 * t6 * t64 * t65 - 0.5400e4 * B * t500 * t502 + 0.985905561600e12 * t15 * u * T - 0.361498705920e12 * t77 * t5 + 0.112968345600e12 * t72 * t10 - 0.29903385600e11 * t97 * t4 - 0.1224115200e10 * t90 * t86 + 0.183617280e9 * t201 * t184 + 0.6645196800e10 * t106 * t61 + 0.1987200e7 * t765 * t14 - 0.129600e6 * t66 * t2 - 0.21859200e8 * t660 * t3 + 0.5400e4 * t772 * A + 0.6120576000e10 * t201 * t36 - 0.2372335257600e13 * t97 * t87 + 0.418647398400e12 * t106 * t208 + 0.10844961177600e14 * t72 * t117;
phydbl t794 = t18 * t3;
phydbl t823 = -0.459043200e9 * t660 * t32 - 0.58145472000e11 * t90 * t173 + 0.162674417664000e15 * t77 * t391 - 0.569360461824000e15 * t77 * t794 - 0.40668604416000e14 * t77 * t388 - 0.569360461824000e15 * t80 * t109 + 0.126524547072000e15 * t80 * t42 + 0.21859200e8 * t765 * t146 - 0.379573641216000e15 * t77 * t96 - 0.455488369459200e15 * t100 * t80 + 0.91097673891840e14 * t260 * t77 - 0.14234011545600e14 * t794 * t72 + 0.1674589593600e13 * t131 * t97 - 0.139549132800e12 * t571 * t106 + 0.7344691200e10 * t137 * t90 + 0.379573641216000e15 * t77 * t237 + 0.569360461824000e15 * t77 * t527;
phydbl t826 = t89 * t14;
phydbl t835 = t122 * t86;
phydbl t850 = t140 * t184;
phydbl t859 = t75 * t2;
phydbl t864 = -0.103520083968000e15 * t643 * t615 + 0.380588544000e12 * t332 * t826 + 0.12336143339520e14 * t185 * t77 - 0.15913497600e11 * t268 * t106 + 0.201570969600e12 * t237 * t97 - 0.65792764477440e14 * t835 * t80 - 0.388744012800e12 * t154 * t72 + 0.795674880e9 * t123 * t90 - 0.1814138726400e13 * t234 * t72 - 0.8811530956800e13 * t162 * t77 + 0.1036650700800e13 * t188 * t72 + 0.4546713600e10 * t141 * t106 + 0.56393798123520e14 * t850 * t80 - 0.86387558400e11 * t191 * t97 + 0.21596889600e11 * t159 * t97 - 0.1468588492800e13 * t696 * t77 - 0.4699483176960e13 * t859 * t80 + 0.4405765478400e13 * t151 * t77;
phydbl t866 = t68 * t14;
phydbl t905 = 0.15664943923200e14 * t866 * t80 - 0.1345761091584000e16 * t311 * t548 + 0.24482304000e11 * t7 * t88 * t826 - 0.1836172800e10 * t7 * t199 * t562 + 0.87436800e8 * t7 * t53 * t605 + 0.54224805888000e14 * t303 * T * t10 - 0.18840655282176000e17 * t80 * t794 + 0.5383044366336000e16 * t80 * t391 - 0.113043931693056000e18 * t80 * t134 + 0.12919306479206400e17 * t80 * t117 - 0.45217572677222400e17 * t80 * t103 + 0.90435145354444800e17 * t80 * t131 + 0.90435145354444800e17 * t80 * t213 - 0.45217572677222400e17 * t80 * t268 + 0.20704016793600e14 * t106 * t34 - 0.2173921763328000e16 * t72 * t242 - 0.11304393169305600e17 * t77 * t571;
phydbl t938 = t6 * t88;
phydbl t949 = 0.2898562351104000e16 * t72 * t28 - 0.310560251904000e15 * t97 * t38 + 0.114176563200e12 * t201 * t146 - 0.3768131056435200e16 * t77 * t87 - 0.1941001574400e13 * t90 * t32 + 0.20704016793600e14 * t106 * t36 + 0.869568705331200e15 * t72 * t208 + 0.12919306479206400e17 * t80 * t141 + 0.58145472000e11 * t90 * t177 + 0.1941001574400e13 * t90 * t30 + 0.11304393169305600e17 * t77 * t196 - 0.31056025190400e14 * t106 * t165 - 0.16267441766400e14 * t381 * t55 * t4 + 0.4066860441600e13 * t649 * t71 * t61 + 0.310560251904000e15 * t97 * t170 - 0.18361728000e11 * t938 * t89 * t3 + 0.1836172800e10 * t6 * t199 * t200 * t14 - 0.837294796800e12 * t372 * t56 * t86;
phydbl t979 = 0.139549132800e12 * t347 * t105 * t184 - 0.2173921763328000e16 * t72 * t128 - 0.131155200e9 * t6 * t53 * t57 * t2 + 0.5961600e7 * t30 * t765 - 0.155280125952000e15 * t97 * t173 + 0.869568705331200e15 * t72 * t137 - 0.29905802035200e14 * t25 - 0.2300446310400e13 * t285 - 0.12816772300800e14 * t25 * t2 - 0.355280928178176000e18 * t140 - 0.2131685569069056000e19 * t72 - 0.716246351207202816000e21 * t165 - 0.10800e5 * t772 - 0.59687195933933568000e20 * t77 + 0.4774975674714685440000e22 * t146 + 0.99478659889889280000e20 * t177 - 0.99478659889889280000e20 * t173;
phydbl t998 = -0.59213488029696000e17 * t97 - 0.151829456486400e15 * t859 + 0.198957319779778560000e21 * t170 + 0.17053484552552448000e20 * t137 - 0.42633711381381120000e20 * t242 + 0.17053484552552448000e20 * t208 - 0.42633711381381120000e20 * t128 - 0.12434832486236160000e20 * t227 + 0.12434832486236160000e20 * t536 + 0.7460899491741696000e19 * t196 - 0.7460899491741696000e19 * t571 + 0.2486966497247232000e19 * t123 - 0.2486966497247232000e19 * t87 - 0.25880020992000e14 * t90 - 0.431333683200e12 * t201 + 0.315805269491712000e18 * t117 + 0.43064354930688000e17 * t245 - 0.75362621128704000e17 * t234;
phydbl t1017 = -0.1138720923648000e16 * t83 - 0.2763296108052480000e19 * t134 - 0.16149133099008000e17 * t154 + 0.43064354930688000e17 * t188 - 0.1105318443220992000e19 * t268 + 0.2210636886441984000e19 * t131 + 0.3588696244224000e16 * t69 + 0.90435145354444800e17 * t260 - 0.75362621128704000e17 * t93 - 0.1105318443220992000e19 * t103 - 0.6343142400e10 * t660 + 0.1821953477836800e16 * t850 - 0.83462400e8 * t765 - 0.2125612390809600e16 * t835 - 0.1193743918678671360000e22 * t80 + 0.1821953477836800e16 * t412 + 0.315805269491712000e18 * t141;
phydbl t1038 = 0.2210636886441984000e19 * t213 - 0.331595532966297600e18 * t96 - t53 * t104 * t57 * t105 + 0.506098188288000e15 * t866 - 0.16149133099008000e17 * t109 - 0.497393299449446400e18 * t794 + 0.35528092817817600e17 * t159 - 0.142112371271270400e18 * t191 + 0.331595532966297600e18 * t237 + 0.497393299449446400e18 * t527 + 0.142112371271270400e18 * t391 - 0.35528092817817600e17 * t388 - 0.4934457335808000e16 * t394 + 0.1644819111936000e16 * t399 - 0.328963822387200e15 * t112 - 0.9868914671616000e16 * t162 + 0.328963822387200e15 * t76 - 0.1644819111936000e16 * t696;
phydbl t1070 = 0.4934457335808000e16 * t151 + 0.13816480540262400e17 * t185 - 0.13816480540262400e17 * t100 + 0.281968990617600e15 * t122 * t61 - 0.117487079424000e15 * t68 * t3 + 0.46994831769600e14 * t75 * t14 + 0.211476742963200e15 * t24 * t184 - 0.281968990617600e15 * t140 * t86 - 0.358869624422400e15 * t75 - 0.358869624422400e15 * t5 + 0.29905802035200e14 * t15 + 0.155280125952000e15 * t97 * t177 - 0.18840655282176000e17 * t77 * t227 + 0.18840655282176000e17 * t77 * t536 + 0.3768131056435200e16 * t77 * t123 - 0.12560436854784000e17 * t80 * t96 + 0.12560436854784000e17 * t80 * t237;
phydbl t1093 = B * t453;
phydbl t1119 = -0.10844961177600e14 * t311 * T * t5 + 0.3614987059200e13 * t565 * t55 * t10 - 0.1016715110400e13 * t643 * t71 * t4 + 0.239227084800e12 * t297 * t56 * t61 - 0.46516377600e11 * t366 * t105 * t86 + 0.7344691200e10 * t332 * t89 * t184 + 0.18840655282176000e17 * t80 * t527 - 0.5383044366336000e16 * t80 * t191 - 0.5961600e7 * t1093 * t455 * t2 + 0.259200e6 * B * t64 * t65 * A - 0.918086400e9 * t561 * t200 * t3 + 0.87436800e8 * t604 * t57 * t14 + 0.3235002624000e13 * t106 * t177 - 0.310560251904000e15 * t72 * t571 - 0.1449281175552000e16 * t77 * t268 + 0.517600419840000e15 * t72 * t536 - 0.51760041984000e14 * t97 * t242 + 0.31715712000e11 * t201 * t30;
phydbl t1160 = -0.1449281175552000e16 * t77 * t103 - 0.7117005772800e13 * t97 * t571 - 0.37957364121600e14 * t72 * t268 - 0.164317593600e12 * t24 * t18 + 0.164317593600e12 * t4 * t184 + 0.215666841600e12 * t514 * t200 + 0.1256043685478400e16 * t429 * t71 - 0.3450669465600e13 * t449 * t89 + 0.15790263474585600e17 * t61 * u * T - 0.241546862592000e15 * t467 * t56 + 0.34506694656000e14 * t496 * t105 - 0.5024174741913600e16 * t86 * t51 * t55 - 0.215666841600e12 * t561 * t200 - 0.34506694656000e14 * t610 * t105 - 0.5024174741913600e16 * t122 * t51 * t55 - 0.241546862592000e15 * t11 * t52 * t56 - 0.1256043685478400e16 * t18 * t70 * t71;
phydbl t1200 = -0.15790263474585600e17 * t140 * u * T - 0.3450669465600e13 * t938 * t89 + 0.12940010496000e14 * A * t88 * t89 - 0.9420327641088000e16 * t3 * t70 * t71 - 0.57088281600e11 * t201 * t6 - 0.538304436633600e15 * t77 * t140 - 0.31056025190400e14 * t97 * t18 - 0.5176004198400e13 * t106 * t11 - 0.144928117555200e15 * t72 * t122 - 0.647000524800e12 * t90 * t7 - 0.3171571200e10 * t660 * B - 0.144928117555200e15 * t72 * t86 + 0.647000524800e12 * t90 * t14 - 0.5176004198400e13 * t106 * t3 + 0.538304436633600e15 * t77 * t61 - 0.57088281600e11 * t201 * t2 + 0.31056025190400e14 * t97 * t184 + 0.3171571200e10 * t660 * A;
phydbl t1239 = -0.1614913309900800e16 * t80 * t4 - 0.1614913309900800e16 * t80 * t24 - 0.985905561600e12 * t25 * u * T - 0.361498705920e12 * t75 * t51 * t55 - 0.112968345600e12 * t68 * t70 * t71 - 0.29903385600e11 * t24 * t52 * t56 - 0.6645196800e10 * t140 * t104 * t105 - 0.12652454707200e14 * t80 * t5 - 0.1355620147200e13 * t72 * t24 - 0.12652454707200e14 * t80 * t75 - 0.1224115200e10 * t122 * t88 * t89 - 0.1530144000e10 * t201 * t11 - 0.119374391867867136000e21 * t3 - 0.19895731977977856000e20 * t18 + 0.19895731977977856000e20 * t184 - 0.2300446310400e13 * t21 - 0.338905036800e12 * t97 * t140;
phydbl t1277 = -0.153014400e9 * t660 * t7 - 0.69774566400e11 * t106 * t122 + 0.11629094400e11 * t90 * t184 - 0.1530144000e10 * t201 * t3 + 0.153014400e9 * t660 * t14 - 0.1355620147200e13 * t72 * t4 - 0.4518733824000e13 * t77 * t68 - 0.69774566400e11 * t106 * t86 - 0.11629094400e11 * t90 * t18 - 0.10929600e8 * t765 * t6 - 0.10929600e8 * t765 * t2 + 0.496800e6 * t66 * A + 0.338905036800e12 * t97 * t61 - 0.496800e6 * t66 * B + 0.4518733824000e13 * t77 * t10 - 0.51760041984000e14 * t77 * t24 - 0.10571904000e11 * t201 * t7 - 0.183617280e9 * t18 * t199 * t200;
phydbl t1314 = -0.3450669465600e13 * t97 * t122 - 0.95147136000e11 * t90 * t11 - 0.14788583424000e14 * t72 * t140 - 0.647000524800e12 * t106 * t18 - 0.834624000e9 * t660 * t6 - 0.149529010176000e15 * t80 * t68 - 0.3450669465600e13 * t97 * t86 + 0.10571904000e11 * t201 * t14 - 0.51760041984000e14 * t77 * t4 + 0.647000524800e12 * t106 * t184 - 0.95147136000e11 * t90 * t3 + 0.14788583424000e14 * t72 * t61 - 0.41731200e8 * t1093 * t455 - 0.834624000e9 * t660 * t2 + 0.41731200e8 * t510 * t455 + 0.149529010176000e15 * t80 * t10 - 0.113667840e9 * t140 * t88 * t89;
phydbl t1369 = -0.71204290560e11 * t285 * u * T - 0.8638755840e10 * t75 * t70 * t71 - 0.2399654400e10 * t68 * t52 * t56 - 0.26701608960e11 * t25 * t51 * t55 - 0.568339200e9 * t24 * t104 * t105 - 0.18944640e8 * t122 * t199 * t200 + 0.26701608960e11 * t15 * t51 * t55 - 0.568339200e9 * t4 * t104 * t105 + 0.2399654400e10 * t10 * t52 * t56 - 0.71204290560e11 * t21 * u * T + 0.113667840e9 * t61 * t88 * t89 - 0.8638755840e10 * t5 * t70 * t71 - 0.18944640e8 * t86 * t199 * t200 + 0.2583360e7 * t184 * t53 * t57 - 0.280800e6 * t3 * t453 * t455 + 0.23400e5 * t14 * t64 * t65 - 0.1404e4 * t2 * t500 * t502 + 0.54e2 * A * t54 * t58;
*mean = -t1 * (t126 + t1160 + t823 + t665 + t621 + t1369 + t267 + t1070 + t1314 + t475 + t783 + t1200 + t703 + t1119 + t740 + t582 + t406 + t1239 + t1017 + t979 + t539 + t309 + t224 + t50 + t1277 + t905 + t864 + t998 + t949 + t1038 + t362 + t176) / 0.14324927024144056320000e23;
 

 /* printf("\n. Taylor: %f",*mean); */
 
 /*  /\* C(int(exp((B-A)*s/T+(1/2)*u*(T-s)*s/T), s = 0 .. T)) *\/ */
 /*  /\* Correct but numerically unstable due to EXP(TOO BIG) *\/ */

 /* *mean = */
 /*   SQRT(0.2e1) * */
 /*   SQRT(0.3141592654e1) * */
 /*   EXP((double) ((4 * B * B - 8 * B * A + 4 * B * u * T + 4 * A * A - 4 * A * u * T + u * u * T * T) / u / T) / 0.8e1) * */
 /*   (-erf(SQRT(0.2e1) * (double) (-2 * B + 2 * A - u * T) / (double) T * pow((double) (u / T), -0.1e1 / 0.2e1) / 0.4e1) */
 /*    +erf(SQRT(0.2e1) * (double) (u *  T - 2 * B + 2 * A) / (double) T * pow((double) (u / T), -0.1e1 / 0.2e1) / 0.4e1))* */
 /*   pow((double) (u / T), -0.1e1 / 0.2e1) / 0.2e1; */
 
 /* *mean /= T; */
 /* *mean *= EXP(A); */

 /* printf("\nErf: %f [%f %f %f]", */
 /*        *mean, */
 /*        EXP((double) ((4 * B * B - 8 * B * A + 4 * B * u * T + 4 * A * A - 4 * A * u * T + u * u * T * T) / u / T) / 0.8e1), */
 /*        (double) ((4 * B * B - 8 * B * A + 4 * B * u * T + 4 * A * A - 4 * A * u * T + u * u * T * T) / u / T) / 0.8e1, */
 /*        (-erf(SQRT(0.2e1) * (double) (-2 * B + 2 * A - u * T) / (double) T * pow((double) (u / T), -0.1e1 / 0.2e1) / 0.4e1) + erf(SQRT(0.2e1) * (double) (u * T - 2 * B + 2 * A) / (double) T * pow((double) (u / T), -0.1e1 / 0.2e1) / 0.4e1)) * pow((double) (u / T), -0.1e1 / 0.2e1) / 0.2e1); */

 }
 
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Integrated_Geometric_Brownian_Bridge_Var(phydbl T, phydbl A, phydbl B, phydbl u, phydbl *var)
 { 
phydbl t2 = exp(0.2e1 * A);
phydbl t3 = A * A;
phydbl t7 = u * u;
phydbl t8 = B * t7;
phydbl t9 = T * T;
phydbl t12 = B * B;
phydbl t13 = t7 * u;
phydbl t14 = t7 * t7;
phydbl t15 = t14 * t14;
phydbl t16 = t15 * t13;
phydbl t18 = t9 * T;
phydbl t19 = t9 * t9;
phydbl t20 = t19 * t19;
phydbl t21 = t20 * t18;
phydbl t24 = t12 * t12;
phydbl t25 = t24 * B;
phydbl t29 = t15 * u;
phydbl t31 = t20 * T;
phydbl t34 = t12 * B;
phydbl t35 = t15 * t7;
phydbl t37 = t20 * t9;
phydbl t40 = t15 * t14;
phydbl t42 = t20 * t19;
phydbl t57 = t3 * A;
phydbl t58 = t3 * t3;
phydbl t59 = t58 * t58;
phydbl t60 = t59 * t57;
phydbl t64 = t7 * t9;
phydbl t65 = t59 * t3;
phydbl t68 = t13 * t18;
phydbl t69 = t59 * A;
phydbl t72 = t14 * t19;
phydbl t75 = t14 * t7;
phydbl t76 = t19 * t9;
phydbl t77 = t75 * t76;
phydbl t78 = t58 * t3;
phydbl t81 = -0.4722222240e18 * t3 * u * T - 0.2222222160e18 * t8 * t9 - 0.2370000e7 * t12 * t16 * t21 - 0.745270000e9 * t25 * t15 * t20 - 0.143603000e9 * t24 * t29 * t31 - 0.13704400e8 * t34 * t35 * t37 + 0.220000e6 * B * t40 * t42 - 0.3230000000e10 * t24 * t15 * t20 - 0.565000000e9 * t34 * t29 * t31 - 0.100000000e9 * t12 * t35 * t37 - 0.85000e5 * B * t16 * t21 + 0.6891165130e12 * t60 * u * T - 0.8770304000e12 * t64 * t65 + 0.6541049000e12 * t68 * t69 - 0.3577430000e12 * t72 * t59 - 0.5328000000e11 * t77 * t78;
phydbl t82 = t14 * t13;
phydbl t83 = t19 * t18;
phydbl t84 = t82 * t83;
phydbl t85 = t58 * A;
phydbl t88 = t14 * u;
phydbl t89 = t19 * T;
phydbl t90 = t88 * t89;
phydbl t91 = t58 * t57;
phydbl t94 = t29 * t31;
phydbl t97 = t35 * t37;
phydbl t100 = t15 * t20;
phydbl t103 = t16 * t21;
phydbl t106 = t34 * t88;
phydbl t107 = t89 * t58;
phydbl t110 = t34 * t7;
phydbl t111 = t9 * t91;
phydbl t114 = t34 * t13;
phydbl t115 = t18 * t78;
phydbl t118 = B * u;
phydbl t119 = T * t59;
phydbl t123 = t76 * t57;
phydbl t127 = t83 * t3;
phydbl t131 = t20 * A;
phydbl t134 = t12 * u;
phydbl t138 = u * T;
phydbl t139 = t25 * t58;
phydbl t142 = t91 * t12;
phydbl t145 = t24 * t58;
phydbl t148 = 0.1479000000e11 * t84 * t85 + 0.1539200000e12 * t90 * t91 + 0.587000000e9 * t94 * t57 - 0.92460000e8 * t97 * t3 - 0.3259400000e10 * t100 * t58 + 0.10000000e8 * t103 * A - 0.5386900000e13 * t106 * t107 + 0.1052436300e15 * t110 * t111 - 0.5494474000e14 * t114 * t115 - 0.2664291708e15 * t118 * t119 + 0.1065273000e13 * t34 * t75 * t123 - 0.1479180000e12 * t34 * t82 * t127 + 0.1338774000e11 * t34 * t15 * t131 + 0.3790140862e14 * t134 * T * t69 - 0.3730008441e16 * t138 * t139 + 0.1065716712e16 * t138 * t142 - 0.1197449022e17 * t138 * t145;
phydbl t150 = B * t91;
phydbl t153 = t12 * t78;
phydbl t156 = t25 * t57;
phydbl t159 = t34 * t85;
phydbl t162 = t24 * t12;
phydbl t163 = t162 * t3;
phydbl t166 = t34 * A;
phydbl t169 = t24 * t3;
phydbl t172 = t25 * t3;
phydbl t175 = t34 * t57;
phydbl t178 = t34 * t3;
phydbl t181 = B * A;
phydbl t184 = B * t78;
phydbl t187 = B * t3;
phydbl t190 = B * t57;
phydbl t193 = B * t85;
phydbl t196 = t24 * t34;
phydbl t197 = t196 * A;
phydbl t200 = t24 * A;
phydbl t203 = 0.1368513200e16 * t138 * t150 - 0.4789796130e16 * t138 * t153 + 0.9579592368e16 * t138 * t156 + 0.9579592019e16 * t138 * t159 - 0.4789796148e16 * t138 * t163 + 0.6456120000e14 * t90 * t166 - 0.1780077270e16 * t68 * t169 - 0.4005331107e16 * t64 * t172 + 0.2373436550e16 * t68 * t175 - 0.5114963006e15 * t72 * t178 + 0.1197000000e13 * t84 * t181 - 0.1335110261e16 * t64 * t184 - 0.1115603000e14 * t77 * t187 + 0.6455912700e14 * t90 * t190 + 0.7120309141e15 * t68 * t193 + 0.1368513201e16 * t138 * t197 + 0.1291400000e13 * t77 * t200;
phydbl t204 = t12 * A;
phydbl t207 = t12 * t85;
phydbl t210 = t12 * t3;
phydbl t213 = t12 * t7;
phydbl t217 = t12 * t13;
phydbl t221 = t12 * t57;
phydbl t224 = t12 * t75;
phydbl t232 = t12 * t14;
phydbl t236 = t12 * t88;
phydbl t240 = t12 * t58;
phydbl t249 = B * t58;
phydbl t252 = t25 * A;
phydbl t257 = t34 * t58;
phydbl t260 = 0.1115831100e14 * t77 * t204 + 0.4005330927e16 * t64 * t207 - 0.9683950000e14 * t90 * t210 - 0.3946635700e14 * t213 * t9 * t59 + 0.2354776600e14 * t217 * t18 * t91 + 0.5114971000e15 * t72 * t221 - 0.7991390000e12 * t224 * t76 * t58 + 0.1481000000e12 * t12 * t82 * t83 * t57 - 0.1001680000e14 * t232 * t19 * t78 + 0.3232259000e13 * t236 * t89 * t85 - 0.1780077100e16 * t68 * t240 - 0.1933540000e11 * t12 * t15 * t20 * t3 + 0.1300000000e10 * t204 * t94 - 0.2557472700e15 * t72 * t249 + 0.7120305100e15 * t68 * t252 + 0.2557476700e15 * t72 * t200 - 0.6675552600e16 * t64 * t257;
phydbl t263 = t24 * t57;
phydbl t266 = t162 * A;
phydbl t269 = t34 * t78;
phydbl t272 = t162 * t57;
phydbl t281 = B * t13;
phydbl t285 = B * t14;
phydbl t289 = B * t88;
phydbl t293 = B * t75;
phydbl t297 = t85 * t24;
phydbl t300 = t196 * t3;
phydbl t303 = B * t29;
phydbl t311 = B * t82;
phydbl t315 = B * t15;
phydbl t321 = 0.6675550630e16 * t64 * t263 + 0.1335110310e16 * t64 * t266 - 0.2486672290e16 * t138 * t269 + 0.2486672264e16 * t138 * t272 - 0.7580281740e13 * t118 * T * t65 + 0.8770300400e13 * t8 * t9 * t69 - 0.5886921000e13 * t281 * t18 * t59 + 0.2861957000e13 * t285 * t19 * t91 - 0.1077430000e13 * t289 * t89 * t78 + 0.3195000000e12 * t293 * t76 * t85 + 0.3730008390e16 * t138 * t297 - 0.1065716699e16 * t138 * t300 - 0.1395750000e10 * t303 * t31 * t3 + 0.152478000e9 * B * t35 * t37 * A - 0.7417900000e11 * t311 * t83 * t58 + 0.1310000000e11 * t315 * t20 * t57 + 0.1971235000e14 * t90 * t200;
phydbl t354 = -0.4879631000e15 * t68 * t172 - 0.9734606000e15 * t64 * t163 + 0.8132711424e15 * t68 * t263 - 0.1652248000e15 * t72 * t169 + 0.6575580000e12 * t84 * t204 - 0.9734607070e15 * t64 * t153 - 0.3000000000e10 * t12 + 0.1000000000e10 * t57 + 0.7331618000e14 * t72 * t263 - 0.1244100000e14 * t90 * t169 - 0.1253983358e17 * t210 * t68 + 0.3439800000e11 * t100 * t204 - 0.1144586120e15 * t68 * t153 + 0.3942690000e14 * t90 * t221 - 0.6457200000e13 * t77 * t210 - 0.2050264784e17 * t249 * t64 + 0.8359887000e16 * t190 * t68;
phydbl t360 = t24 * t24;
phydbl t361 = t360 * A;
phydbl t366 = t360 * t12;
phydbl t367 = t366 * A;
phydbl t374 = t360 * B;
phydbl t381 = t34 * u;
phydbl t397 = 0.2929910910e15 * t181 * t90 + 0.4879629976e15 * t68 * t207 + 0.2664291701e15 * t138 * t361 + 0.9444444500e18 * t181 * t138 + 0.7580281461e13 * t367 * t138 + 0.5833333338e18 * t204 * t138 + 0.2460317400e18 * t181 * t64 + 0.200000e6 * t374 - 0.200000e6 * t69 + 0.100000000e9 * t24 + 0.2222222150e18 * A * t7 * t9 - 0.1944444438e18 * t381 * T - 0.3224206410e17 * t281 * t18 - 0.6398809530e17 * t24 * u * T - 0.1755952430e17 * t217 * t18 - 0.1230158726e18 * t3 * t7 * t9 + 0.3224206250e17 * A * t13 * t18;
phydbl t444 = 0.1944444429e18 * t57 * u * T - 0.1755952280e17 * t3 * t13 * t18 + 0.4894179970e17 * t57 * t7 * t9 - 0.3339945030e16 * t285 * t19 - 0.4894180056e17 * t110 * t9 - 0.4722222240e18 * t134 * T + 0.8333333350e18 * A * u * T - 0.8333333400e18 * t118 * T + 0.3339946425e16 * A * t14 * t19 - 0.6398809540e17 * t58 * u * T + 0.4100529100e16 * t85 * t7 * t9 - 0.1464959090e15 * t3 * t88 * t89 + 0.6839256369e15 * t57 * t14 * t19 - 0.4238315719e16 * t78 * u * T - 0.1897556340e14 * t293 * t76 - 0.4238315750e16 * t162 * u * T - 0.2089971880e16 * t24 * t13 * t18;
phydbl t448 = t34 * t14;
phydbl t478 = t360 * t3;
phydbl t481 = t162 * t58;
phydbl t490 = -0.6839241185e15 * t448 * t19 - 0.4100529178e16 * t25 * t7 * t9 - 0.1464966800e15 * t236 * t89 + 0.1769179920e17 * t85 * u * T - 0.1797238280e16 * t3 * t14 * t19 + 0.2747601900e15 * A * t88 * t89 - 0.1547619080e17 * t58 * t7 * t9 + 0.6812170220e16 * t57 * t13 * t18 - 0.2747571750e15 * t289 * t89 - 0.1547619200e17 * t24 * t7 * t9 - 0.20000000e8 * t25 - 0.30000000e8 * t85 - 0.2114391104e15 * t138 * t478 - 0.9867158400e15 * t138 * t481 + 0.2050264696e17 * t200 * t64 + 0.2542989412e17 * t252 * t138 - 0.6357473500e17 * t240 * t138;
phydbl t507 = t59 * B;
phydbl t510 = t12 * t59;
phydbl t513 = B * t69;
phydbl t520 = t25 * t78;
phydbl t523 = t25 * t85;
phydbl t530 = 0.2508160000e12 * t84 * t190 - 0.1466339000e14 * t72 * t184 + 0.4976610000e13 * t90 * t193 + 0.3270245000e14 * t68 * t150 - 0.3450000000e11 * t100 * t187 - 0.1291300000e13 * t77 * t249 + 0.2074777207e15 * t64 * t142 - 0.7261721108e15 * t64 * t139 - 0.5186943795e14 * t64 * t507 - 0.2114391078e15 * t138 * t510 + 0.4698646924e14 * t138 * t513 + 0.3070300000e10 * t94 * t181 - 0.4841145120e15 * t64 * t269 - 0.3183718246e15 * t520 * t138 + 0.2210116600e15 * t523 * t64 - 0.8241728800e14 * t139 * t68 + 0.2003340000e14 * t156 * t72;
phydbl t553 = t19 * t85;
phydbl t559 = -0.3232400000e13 * t172 * t90 + 0.3196000000e12 * t252 * t77 + 0.4841147005e15 * t64 * t272 + 0.7261720740e15 * t64 * t297 - 0.1626544080e15 * t281 * t115 + 0.4303517000e13 * t293 * t123 - 0.3000000e7 * t360 - 0.3000000e7 * t59 - 0.17000000e8 * t196 - 0.4399012000e14 * t72 * t172 - 0.1144586200e15 * t68 * t163 - 0.48000e5 * t366 - 0.20000e5 * t65 + 0.2003375000e14 * t448 * t553 - 0.10000000e8 * t78 - 0.50000000e8 * t162 + 0.8000000e7 * t91;
phydbl t561 = t59 * t58;
phydbl t574 = A * t82;
phydbl t577 = t85 * t13;
phydbl t580 = t3 * t75;
phydbl t586 = t58 * t14;
phydbl t589 = t57 * t88;
phydbl t602 = -0.100000e6 * t367 - 0.1420e4 * t561 - 0.4325330800e13 * t510 * t68 - 0.1797236000e16 * t232 * t19 - 0.6812172380e16 * t114 * t18 - 0.1769179896e17 * t25 * u * T - 0.1230158720e18 * t213 * t9 + 0.1136992000e13 * t574 * t83 + 0.5357388956e15 * t577 * t18 - 0.1004312000e14 * t580 * t76 + 0.8983686054e15 * t91 * u * T - 0.2049036000e15 * t586 * t19 + 0.5490979267e14 * t589 * t89 - 0.9402891770e15 * t78 * t7 * t9 - 0.1136811000e13 * t311 * t83 - 0.5490960000e14 * t106 * t89 - 0.9402891530e15 * t162 * t7 * t9;
phydbl t644 = -0.2049024100e15 * t24 * t14 * t19 - 0.5357393770e15 * t25 * t13 * t18 - 0.8983686060e15 * t196 * u * T - 0.1004260240e14 * t224 * t76 + 0.1897544000e14 * A * t75 * t76 - 0.2089971420e16 * t58 * t13 * t18 - 0.5984000000e12 * t84 * t12 - 0.1907300810e15 * t64 * t196 - 0.5114985270e14 * t72 * t25 - 0.1614025000e14 * t90 * t24 - 0.1186717058e15 * t68 * t162 - 0.3719300000e13 * t77 * t34 - 0.6049800000e11 * t100 * B - 0.1186716997e15 * t68 * t78 + 0.3719150000e13 * t77 * t57 - 0.1614009620e14 * t90 * t58 + 0.1907300700e15 * t64 * t91;
phydbl t655 = t360 * t34;
phydbl t686 = -0.5983807300e12 * t84 * t3 + 0.5114970240e14 * t72 * t85 + 0.6059000000e11 * t100 * A - 0.1710641462e15 * t138 * t59 - 0.1710641438e15 * t138 * t360 - 0.6891165120e12 * t655 * u * T - 0.8770302200e12 * t366 * t7 * t9 - 0.6541040000e12 * t374 * t13 * t18 - 0.3577404000e12 * t360 * t14 * t19 - 0.1539200000e12 * t196 * t88 * t89 - 0.4698646953e13 * t138 * t65 - 0.4087803000e13 * t68 * t360 - 0.4698647004e13 * t138 * t366 - 0.5328500000e11 * t162 * t75 * t76 - 0.6280000000e11 * t84 * t24 - 0.2094778000e13 * t72 * t196 - 0.1159280000e11 * t100 * t34;
phydbl t723 = -0.8294300000e12 * t90 * t162 + 0.2583700000e12 * t77 * t85 - 0.6265560000e11 * t84 * t58 + 0.1149330000e11 * t100 * t57 - 0.4087804000e13 * t68 * t59 - 0.5763269898e13 * t64 * t374 - 0.8293904000e12 * t90 * t78 - 0.2582655100e12 * t77 * t25 - 0.1521060000e10 * t94 * t12 - 0.1510800000e10 * t94 * t3 + 0.85420000e8 * t97 * A + 0.2094757500e13 * t72 * t91 - 0.78510000e8 * t97 * B + 0.5763268400e13 * t64 * t69 - 0.3476645600e14 * t64 * t360 - 0.2192000000e12 * t84 * t34 - 0.1481620000e11 * t25 * t82 * t83;
phydbl t752 = A * t29;
phydbl t760 = -0.1101512460e14 * t72 * t162 - 0.1076151000e13 * t77 * t24 - 0.2323633800e14 * t68 * t196 - 0.3942769100e13 * t90 * t25 - 0.3179650000e11 * t100 * t12 - 0.2960324185e14 * t138 * t374 - 0.1101515600e14 * t72 * t78 + 0.2192784000e12 * t84 * t57 - 0.3476645322e14 * t64 * t59 + 0.3942680000e13 * t90 * t85 - 0.1076127000e13 * t77 * t58 + 0.2323633627e14 * t68 * t91 - 0.2880000000e10 * t303 * t31 - 0.3137962000e11 * t100 * t3 + 0.2840510000e10 * t752 * t31 + 0.2960324165e14 * t138 * t69 - 0.9682000000e10 * t196 * t75 * t76;
phydbl t765 = t360 * t24;
phydbl t814 = -0.9396930200e11 * t765 * u * T - 0.9611875000e11 * t366 * t13 * t18 - 0.5556330000e11 * t374 * t14 * t19 - 0.1234481500e12 * t655 * t7 * t9 - 0.2563440000e11 * t360 * t88 * t89 - 0.2994200000e10 * t162 * t82 * t83 + 0.1234481600e12 * t60 * t7 * t9 - 0.2563300000e11 * t59 * t88 * t89 + 0.5556248000e11 * t69 * t14 * t19 - 0.9396930200e11 * t561 * u * T + 0.9689500000e10 * t91 * t75 * t76 - 0.9611880000e11 * t65 * t13 * t18 - 0.3021580000e10 * t78 * t82 * t83 + 0.763720000e9 * t85 * t15 * t20 - 0.129241000e9 * t58 * t29 * t31 + 0.7383500e7 * t57 * t35 * t37;
phydbl t829 = t374 * t3;
phydbl t842 = t374 * A;
phydbl t853 = -0.5600000e7 * t3 * t16 * t21 - 0.3213600e7 * A * t40 * t42 - 0.5833333350e18 * t187 * t138 + 0.4304200000e13 * t77 * t166 + 0.1946921629e16 * t64 * t159 + 0.2559523820e18 * t190 * t138 - 0.3790140882e14 * t829 * t138 + 0.3511906000e17 * t181 * t68 + 0.2202988200e15 * t72 * t175 - 0.1468254070e18 * t187 * t64 - 0.3942732000e14 * t90 * t178 - 0.8132711000e15 * t68 * t257 + 0.8770302700e13 * t842 * t64 - 0.3839285591e18 * t210 * t138 + 0.1468253985e18 * t204 * t64 + 0.2559523810e18 * t166 * t138 + 0.1626542320e15 * t68 * t266;
phydbl t871 = t360 * t57;
phydbl t880 = t196 * t58;
phydbl t891 = -0.1652248000e15 * t72 * t240 + 0.1946921488e16 * t64 * t156 - 0.2433651484e16 * t64 * t145 + 0.6609041000e14 * t72 * t252 + 0.2781316600e15 * t64 * t197 - 0.2043651325e17 * t187 * t68 + 0.3594474817e16 * t181 * t72 + 0.6190476200e17 * t190 * t64 + 0.1137042247e15 * t871 * t138 - 0.3946636500e14 * t478 * t64 - 0.1769179894e18 * t178 * t138 + 0.5886927000e13 * t361 * t68 - 0.2274084525e15 * t880 * t138 - 0.9285714360e17 * t210 * t64 + 0.2043651805e17 * t204 * t68 + 0.1769179900e18 * t221 * t138 - 0.8845899000e17 * t249 * t138;
phydbl t898 = t162 * t85;
phydbl t901 = t196 * t57;
phydbl t928 = 0.8845899300e17 * t200 * t138 + 0.6190475974e17 * t166 * t64 + 0.2861967000e13 * t197 * t72 + 0.3183718369e15 * t898 * t138 + 0.1052436420e15 * t901 * t64 - 0.2354776200e14 * t300 * t68 + 0.2583300000e13 * t77 * t221 + 0.4399010000e14 * t72 * t207 - 0.3755200000e12 * t84 * t210 + 0.2509100000e12 * t84 * t166 + 0.5186945270e14 * t64 * t361 + 0.2542989379e17 * t193 * t138 - 0.2051769370e16 * t187 * t72 + 0.2289169200e15 * t68 * t159 + 0.2051768000e16 * t204 * t72 + 0.2289167900e15 * t68 * t156 - 0.1244097000e14 * t90 * t240;
phydbl t947 = t34 * t91;
phydbl t958 = t24 * t78;
phydbl t967 = 0.1466349000e14 * t72 * t266 - 0.7331650000e14 * t72 * t257 - 0.2583029000e13 * t77 * t178 + 0.1658907000e14 * t90 * t175 - 0.1841763540e15 * t481 * t64 + 0.5494474200e14 * t272 * t68 - 0.4100529010e17 * t178 * t64 - 0.6357473383e17 * t169 * t138 + 0.5638376416e15 * t138 * t947 + 0.3270242900e14 * t68 * t197 + 0.4976387700e13 * t90 * t252 + 0.4100529000e17 * t221 * t64 - 0.2861462400e15 * t68 * t145 - 0.9867158433e15 * t138 * t958 + 0.8476631585e17 * t175 * t138 + 0.1184059034e16 * t138 * t523 + 0.8359888510e16 * t166 * t68;
phydbl t978 = t65 * B;
phydbl t985 = B * t60;
phydbl t992 = t59 * t34;
phydbl t995 = t24 * t59;
phydbl t998 = t69 * t12;
phydbl t1001 = t34 * t69;
phydbl t1006 = 0.4698646863e14 * t138 * t842 - 0.1001659000e14 * t163 * t72 + 0.1077350000e13 * t266 * t90 + 0.5638376400e15 * t138 * t901 - 0.360000e6 * t880 + 0.10000e5 * t60 - 0.1357930100e13 * t978 * t64 + 0.9611860000e12 * t513 * t68 + 0.2050665000e12 * t150 * t90 + 0.1127631630e13 * t985 * t138 - 0.5000541000e12 * t507 * t72 + 0.2000230000e13 * t142 * t72 - 0.2036894700e14 * t992 * t64 - 0.4651480740e14 * t995 * t138 + 0.6789651000e13 * t998 * t64 + 0.2067324610e14 * t1001 * t138 + 0.1153425500e14 * t947 * t68;
phydbl t1008 = t91 * t24;
phydbl t1011 = t25 * t91;
phydbl t1014 = t12 * t65;
phydbl t1070 = 0.4073788600e14 * t1008 * t64 + 0.7442369200e14 * t1011 * t138 - 0.6201974310e13 * t1014 * t138 + 0.1809000000e11 * t85 * t82 * t83 * B + 0.7000817000e13 * t85 * t14 * t19 * t24 + 0.1435633000e13 * t85 * t88 * t89 * t34 + 0.2422191000e14 * t577 * t18 * t25 + 0.2034000000e12 * t85 * t75 * t76 * t12 - 0.54366000e8 * t3 * t35 * t37 * B - 0.4510900000e11 * t3 * t82 * t83 * t24 - 0.7466500000e10 * t3 * t15 * t20 * t34 - 0.2035185700e12 * t580 * t76 * t25 - 0.1040000000e10 * t3 * t29 * t31 * t12 - 0.4514750000e11 * t58 * t82 * t83 * t12 + 0.612330000e9 * t57 * t29 * t31 * B - 0.7000850000e13 * t586 * t19 * t25 - 0.3850000000e10 * t58 * t15 * t20 * B;
phydbl t1124 = 0.7537400000e10 * t57 * t15 * t20 * t12 - 0.3386100000e12 * t58 * t75 * t76 * t34 - 0.1794360000e13 * t58 * t88 * t89 * t24 + 0.5995720000e11 * t57 * t82 * t83 * t34 + 0.3388054000e12 * t57 * t75 * t76 * t24 + 0.1435490000e13 * t589 * t89 * t25 + 0.16000000e8 * A * t16 * t21 * B + 0.3773240000e10 * A * t15 * t20 * t24 + 0.692700000e9 * t752 * t31 * t34 + 0.1795620000e11 * t574 * t83 * t25 + 0.79693000e8 * A * t35 * t37 * t12 - 0.1971335000e14 * t289 * t107 - 0.2074777361e15 * t64 * t300 + 0.8241726800e14 * t297 * t68 - 0.2504208000e14 * t145 * t72 + 0.2274084510e15 * t1008 * t138 - 0.1841763400e15 * t958 * t64;
phydbl t1162 = 0.5386800000e13 * t263 * t90 - 0.7989000000e12 * t169 * t77 + 0.5641734828e16 * t193 * t64 + 0.7417430000e11 * t200 * t84 - 0.1137042238e15 * t381 * t119 + 0.8196118110e15 * t190 * t72 - 0.1647262473e15 * t187 * t90 - 0.6288580100e16 * t184 * t138 + 0.2008619000e14 * t181 * t77 - 0.5357392640e16 * t178 * t68 - 0.6579000000e12 * t311 * t127 + 0.2781316526e15 * t8 * t111 + 0.1880578242e17 * t175 * t64 - 0.1886574102e17 * t172 * t138 - 0.1410433708e17 * t169 * t64 + 0.8196108200e15 * t166 * t72 - 0.1229418200e16 * t210 * t72;
phydbl t1194 = 0.1886574129e17 * t207 * t138 + 0.1647291050e15 * t204 * t90 + 0.5357390500e16 * t221 * t68 - 0.1410433708e17 * t240 * t64 - 0.2678697000e16 * t249 * t68 + 0.5641735040e16 * t252 * t64 + 0.2678700000e16 * t200 * t68 + 0.3144290150e17 * t263 * t138 - 0.3144290124e17 * t257 * t138 + 0.6288580291e16 * t266 * t138 + 0.6609036000e14 * t285 * t553 + 0.6340000000e11 * t315 * t131 + 0.1000000e7 * t871 + 0.2300000e7 * t898 - 0.4000000e7 * t520 + 0.220000e6 * t162 * t91 - 0.140000e6 * t374 * t58;
phydbl t1224 = -0.40000e5 * t366 * t57 - 0.180000e6 * t360 * t85 + 0.1000e4 * t655 * t3 + 0.990e3 * t765 * A - 0.211100e6 * t196 * t78 - 0.30000e5 * t65 * t34 + 0.133900e6 * t69 * t24 + 0.20860e5 * t60 * t12 + 0.254000e6 * t59 * t25 - 0.1020e4 * t561 * B + 0.600000e6 * t995 - 0.20000e5 * t513 - 0.200000e6 * t1001 + 0.117000e6 * t1014 - 0.22790e5 * t985 - 0.28500e5 * t40 * t42 - 0.114000000e9 * t97;
phydbl t1235 = t162 * t78;
phydbl t1250 = t196 * t85;
phydbl t1259 = t366 * t3;
phydbl t1262 = -0.101e3 * t360 * t25 + 0.76e2 * t59 * t85 + 0.5703305000e14 * t898 * t64 - 0.7175350000e12 * t163 * t90 + 0.4667200000e13 * t272 * t72 - 0.8682763800e14 * t1235 * t138 - 0.4325349000e13 * t478 * t68 + 0.6776140000e11 * t266 * t77 - 0.2018495000e14 * t481 * t68 - 0.4073788100e14 * t880 * t64 + 0.1153425000e14 * t901 * t68 + 0.2050823000e12 * t197 * t90 + 0.7442369200e14 * t1250 * t138 - 0.2000230000e13 * t300 * t72 + 0.5000544000e12 * t361 * t72 - 0.6789650400e13 * t829 * t64 - 0.6201974000e13 * t1259 * t138;
phydbl t1277 = t360 * t58;
phydbl t1296 = 0.2036894600e14 * t871 * t64 + 0.2067324700e14 * t374 * t57 * t138 + 0.9611900000e12 * t842 * t68 + 0.1357930000e13 * t367 * t64 + 0.1127631640e13 * t655 * A * t138 - 0.4651480610e14 * t1277 * t138 - 0.6775000000e11 * t184 * t77 - 0.2018493300e14 * t958 * t68 - 0.4667170000e13 * t269 * t72 - 0.5703305700e14 * t520 * t64 - 0.7176804000e12 * t153 * t90 + 0.300000e6 * t978 - 0.10000e5 * t655 - 0.1900e4 * t765 - 0.170000e6 * t829 - 0.2747570000e15 * t90 - 0.3224206200e17 * t68;
phydbl t1314 = 0.1000000000e10 * t210 - 0.5246000e7 * t103 - 0.2222222269e18 * t64 - 0.100000000e9 * t200 + 0.400000000e9 * t249 - 0.3339948000e16 * t72 - 0.30000e5 * t1259 + 0.1000000000e10 * t221 - 0.200000000e9 * t252 - 0.200000000e9 * t169 - 0.100000000e9 * t193 + 0.270000000e9 * t240 + 0.100000000e9 * t257 - 0.100000000e9 * t263 - 0.210000000e9 * t207 + 0.40000000e8 * t172 + 0.20000000e8 * t266;
phydbl t1333 = 0.30000000e8 * t184 - 0.1897466160e14 * t77 - 0.1137580000e13 * t84 + 0.5000000e7 * t150 - 0.1000000e7 * t947 - 0.25000000e8 * t481 + 0.510000e6 * t1277 - 0.20000000e8 * t145 - 0.200000e6 * t478 + 0.1000000e7 * t901 + 0.20000000e8 * t163 - 0.130000000e9 * t156 - 0.700000e6 * t842 - 0.25000000e8 * t523 + 0.31000000e8 * t958 + 0.20000000e8 * t153 - 0.6034000000e11 * t100;
phydbl t1354 = 0.1000000e7 * t1250 - 0.2889000000e10 * t94 - 0.600000e6 * t1235 - 0.8333333300e18 * t138 + 0.700000e6 * t1011 - 0.351766e6 * t15 * t88 * t20 * t89 + 0.29000000e8 * t139 - 0.1000000e7 * t361 - 0.2000000e7 * t300 + 0.20000000e8 * t272 - 0.70000000e8 * t297 + 0.21000000e8 * t142 - 0.1000000e7 * t507 - 0.994000e6 * t992 - 0.130000e6 * t998 - 0.1000000000e10 * t34 - 0.3000000000e10 * t3;
*var = -0.1000000000e-18 * t2 * (t891 + t853 + t1296 + t81 + t1070 + t928 + t1314 + t148 + t760 + t203 + t814 + t1224 + t686 + t1354 + t1162 + t321 + t1006 + t559 + t530 + t967 + t260 + t397 + t1262 + t644 + t1333 + t1194 + t354 + t602 + t1124 + t723 + t490 + t444);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// Inverse method to sample from X where P(X=xi)=pi[i]

int Sample_i_With_Proba_pi(phydbl *pi, int len)
{
  phydbl *cum_pi;
  int i;
  phydbl u;
  
  cum_pi = (phydbl *)mCalloc(len,sizeof(phydbl));

  u = .0;
  For(i,len) u += pi[i];  
  For(i,len) cum_pi[i] = pi[i] / u;
  for(i=1;i<len;i++) cum_pi[i] += cum_pi[i-1];

  if((cum_pi[i-1] > 1. + 1.E-10) || (cum_pi[i-1] < 1. - 1.E-10))
    {
      PhyML_Printf("\n== Sum of probabilities is different from 1.0 (%f).",cum_pi[i-1]);
      PhyML_Printf("\n== Err. in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");
    }

  i = 0;
  u = Uni();
  For(i,len) if(cum_pi[i] > u) break;

  if(i == len)
    {
      PhyML_Printf("\n== Len = %d",len);
      PhyML_Printf("\n== Err. in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");
    }
  Free(cum_pi);

  return(i);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

// Return the value y such that Prob(x<y) = p
phydbl Quantile(phydbl *x, int len, phydbl p)
{
  phydbl *y,q;
  int i;
  int swap;
  phydbl buff;

  y = (phydbl *)mCalloc(len,sizeof(phydbl));
  For(i,len) y[i] = x[i];

  do
    {
      swap = NO;
      For(i,len-1) 
        {
          if(y[i+1] < y[i])
            {
              swap = YES;

              buff = y[i+1];
              y[i+1] = y[i];
              y[i] = buff;
            }
        }
    }
  while(swap == YES);
  
  q = y[(int)((len-1)*p)];

  Free(y);

  return(q);

}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

// Return p such that Prob(x<z) = p
phydbl Prob(phydbl *x, int len, phydbl z)
{
  int i;
  phydbl hit;

  hit = 0.;
  For(i,len) if(x[i] < z) hit+=1.;

  return(hit/(phydbl)len);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

// Return x  where mu is the first moment of the normal density
// and x is the value such that f(x;mu,sigma)=y

phydbl Inverse_Truncated_Normal(phydbl y, phydbl mu, phydbl sigma, phydbl lim_inf, phydbl lim_sup)
{
  phydbl p_inf, p_sup;
  
  p_inf = Pnorm(lim_inf,mu,sigma);
  p_sup = Pnorm(lim_sup,mu,sigma);

  /* return(mu + sigma * SQRT(-LOG( y * y * (p_sup - p_inf) * (p_sup - p_inf) * 2 * PI * sigma * sigma))); */
  return(mu + sigma * SQRT(-LOG( y * y * (p_sup - p_inf) * (p_sup - p_inf) * 2. * PI * sigma * sigma)));
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// Returns a vector with a permutation of all the integer from 0
// to len-1. Fisher-Yates algorithm.

int *Permutate(int len)
{
  int i,pos,tmp;
  int *x;

  x = (int *)mCalloc(len,sizeof(int));

  For(i,len) x[i] = i;

  For(i,len)
    {
      pos = Rand_Int(i,len-1);
      
      tmp    = x[i];
      x[i]   = x[pos];
      x[pos] = tmp;
    }

  return(x);
} 
 
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// Returns the p-value for the Mantel test of correlation between
// matrices x and y.

phydbl Mantel(phydbl *x, phydbl *y, int nrow, int ncol)
{ 
 
  phydbl obs_stat;  // Value of the statistic on the observed data
  phydbl mc_stat;   // Value of the statistics on the Monte Carlo generated data
  int N;
  phydbl sumx, sumy, sumxy, sumxx, sumyy;
  int i,j,k;
  int npermut;
  int *permut;
  phydbl p_val;

  N = nrow*ncol;

  sumx = .0;
  For(i,N) sumx += x[i];

  sumy = .0;
  For(i,N) sumy += y[i];

  sumxx = .0;
  For(i,N) sumxx += x[i]*x[i];

  sumyy = .0;
  For(i,N) sumyy += y[i]*y[i];

  sumxy = .0;
  For(i,N) sumxy += x[i]*y[i];

  obs_stat = (N * sumxy - sumx * sumy) / (SQRT((N-1)*sumxx - (sumx/N)*(sumx/N)) * SQRT((N-1)*sumyy - (sumy/N)*(sumy/N)));
  
  npermut = 1000;
  p_val   = 0.0;
  For(k,npermut)
    {
      permut = Permutate(nrow);

      sumxy = .0;
      For(i,nrow)
        {
          For(j,ncol)
            {
              sumxy += x[i*ncol+j] * y[permut[i]*ncol+permut[j]];
            }
        }
           
      mc_stat = (N * sumxy - sumx * sumy) / (SQRT((N-1)*sumxx - (sumx/N)*(sumx/N)) * SQRT((N-1)*sumyy - (sumy/N)*(sumy/N)));

      Free(permut);

      if(mc_stat > obs_stat) p_val += 1.;
    }
    
  return(p_val / (phydbl)npermut);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl Weighted_Mean(phydbl *x, phydbl *w, int l)
{
  int i;
  phydbl wm;
  wm = .0;
  if(w) For(i,l) wm += x[i]*w[i];
  else  
    {
      For(i,l) wm += x[i];
      wm /= (phydbl)l;
    }
  return(wm);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl Variance(phydbl *x, int l)
{
  phydbl mean,sum;
  int i;

  mean = Weighted_Mean(x,NULL,l);
  sum = 0.0;
  For(i,l) sum += x[i]*x[i];
  
  return(sum/l - mean*mean);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int Sum_Bits(int value, int range)
{
  int i;
  int sum;

  if(range > 8*(int)sizeof(int))
    {
      PhyML_Printf("\n== Err. in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");
    }

  sum = 0;
  For(i,range)
    {
      sum += (value >> i) & 1;
    }

  return(sum);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int Modulo (int a, int b)
{
   if(b < 0)
     return Modulo(-a, -b);   
   int ret = a % b;
   if(ret < 0)
     ret+=b;
   return ret;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Runif_Disk(phydbl *sampled_x, phydbl *sampled_y, phydbl centrx, phydbl centry, phydbl radius)
{
  phydbl r,theta;
  
  r     = Uni();
  theta = Uni()*2.*PI;

  (*sampled_x) = SQRT(r)*COS(theta);
  (*sampled_y) = SQRT(r)*SIN(theta);
  
  (*sampled_x) *= radius;
  (*sampled_y) *= radius;

  (*sampled_x) += centrx;
  (*sampled_y) += centry;

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Random_String(char *s, int len)
{
  int i;
  For(i,len) s[i] = Rand_Int(97,121);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int *Random_Permut(int n)
{
  int *permut;
  int i,j;
  int tmp;

  if(n < 3) 
    {
      PhyML_Printf("\n== Number of vertices in a polygon has to be at least 3.");
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
    }
    
  permut = (int *)mCalloc(n,sizeof(int));

  For(i,n) permut[i] = i;
  
  For(i,n-1)
    {
      j = Rand_Int(i,n-1);
      tmp = permut[i];
      permut[i] = permut[j];
      permut[j] = tmp;      
    }

  return(permut);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* Generate a random polygon with n vertices. Each point lies in [0,1] */
t_poly *Rpoly(int n)
{
  t_poly *p;
  int i;

  p = (t_poly *)Make_Poly(n);  
  p->n_poly_vert = n;
  
  For(i,n)
    {
      p->poly_vert[i]->lonlat[0] = Uni();
      p->poly_vert[i]->lonlat[1] = Uni();
    }
  
  return(p);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl Area_Of_Poly_Monte_Carlo(t_poly **poly, int n_poly, t_geo_coord *lim)
{
  int n_hit,n_trials,trial,i;
  t_geo_coord *point;

  point = (t_geo_coord *)GEO_Make_Geo_Coord(2);

  n_trials = 1E+7;
  trial = 0;
  n_hit = 0;
  do
    {
      point->lonlat[0] = Uni()*lim->lonlat[0];
      point->lonlat[1] = Uni()*lim->lonlat[1];
      For(i,n_poly) if(Is_In_Polygon(point,poly[i]) == YES) { n_hit++; break; }
      trial++;
    }
  while(trial < n_trials);

  Free_Geo_Coord(point);

  return((phydbl)(n_hit)/n_trials);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

int Is_In_Polygon(t_geo_coord *point, t_poly *poly)
{
  int i,j;
  phydbl x,y,x1,y1,x2,y2;
  phydbl x_intersect;
  short int is_in;

  /* Coordinates of the point to test */
  x = point->lonlat[0];
  y = point->lonlat[1];

  j = poly->n_poly_vert-1;
  is_in = NO;
  For(i,poly->n_poly_vert)
    {
      /* Edge of polygon goes from (x1,y1) to (x2,y2) */
      x1 = poly->poly_vert[i]->lonlat[0];
      y1 = poly->poly_vert[i]->lonlat[1];
      x2 = poly->poly_vert[j]->lonlat[0];
      y2 = poly->poly_vert[j]->lonlat[1];

      j = i;

      /* Shoot an horizontal ray to the right. Find out if
         this ray hits the polygon edge */
      if((y1 < y && y2 > y) || (y2 < y && y1 > y))
        {
          /* Coordinates along X-axis of the intersection between ray and edge */
          x_intersect = (y-y1)/(y1-y2)*(x1-x2)+x1;  
          if(x_intersect > x) /* Intersection is on the righthand side */
            is_in = (is_in == YES)?NO:YES;
        }
    }

  return is_in;
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

/* Modified Bessel function of the first kind. Stolen from Numerical Recipes in C. */
phydbl Bessi0(phydbl x)
{
  phydbl ax,ans;
  phydbl y;

  if ((ax=fabs(x)) < 3.75) 
    {
      y=x/3.75;
      y*=y;
      ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492+y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
    } 
  else 
    {
      y=3.75/ax;
      ans=(exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1+y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2+y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1+y*0.392377e-2))))))));
    }

  return ans;
}


/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

/* Modified Bessel function of the second kind (degree 0). Stolen from Numerical Recipes in C. */
phydbl Bessk0(phydbl x)
{
  phydbl y,ans;

  if (x <= 2.0) 
    {
      y=x*x/4.0;
      ans=(-log(x/2.0)*Bessi0(x))+(-0.57721566+y*(0.42278420+y*(0.23069756+y*(0.3488590e-1+y*(0.262698e-2+y*(0.10750e-3+y*0.74e-5))))));
    } 
  else 
    {
      y=2.0/x;
      ans=(exp(-x)/sqrt(x))*(1.25331414+y*(-0.7832358e-1+y*(0.2189568e-1+y*(-0.1062446e-1+y*(0.587872e-2+y*(-0.251540e-2+y*0.53208e-3))))));
    }

  return ans;
}


/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl Euclidean_Dist(t_geo_coord *x, t_geo_coord *y)
{
  int i;
  phydbl dist;
  
  if(x->dim != y->dim) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);    
  
  dist = 0.0;
  For(i,x->dim) dist += POW(x->lonlat[i]-y->lonlat[i],2);

  return(SQRT(dist));
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
/* Return the ranks of elements in x */
int *Ranks(phydbl *x, int len)
{
  int *rk,tmp;
  int i,swap;

  rk = (int *)mCalloc(len,sizeof(int));

  For(i,len) rk[i] = i;

  do
    {
      swap = NO;
      For(i,len-1) 
        {
          if(x[rk[i]] > x[rk[i+1]])
            {
              swap  = YES; 
              tmp   = rk[i];
              rk[i] = rk[i+1];
              rk[i+1] = tmp;
            }
        }
    }
  while(swap == YES);

  return(rk);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl *Brownian_Bridge_Generate(phydbl start, phydbl end, phydbl var, phydbl beg_time, phydbl end_time, int n_steps, phydbl *time)
{
  phydbl *state,end_brown;
  int i;

  if(n_steps == 0) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  if(beg_time > end_time) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  For(i,n_steps-1) if(!(time[i+1] > time[i])) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);

  state = Brownian_Generate(var,n_steps,beg_time,time);
  end_brown = Rnorm(state[n_steps-1],SQRT((time[n_steps-1] - end_time)*var));

  For(i,n_steps)
    {
      state[i] = state[i] - (time[i]/end_time) * end_brown;
      state[i] = start + (end - start)/end_time * time[i] + state[i];
    }

  
  return(state);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl *Brownian_Generate(phydbl var, int n_steps, phydbl beg_time, phydbl *time)
{
  phydbl *state;
  int i;

  if(n_steps == 0) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);

  state = (phydbl *)mCalloc(n_steps,sizeof(phydbl));

  state[0] = Rnorm(0.0,SQRT((time[0] - beg_time)*var));

  for(i=1;i<n_steps;i++)
    {
      state[i] = Rnorm(state[i-1],SQRT((time[i]-time[i-1])*var));
      if(time[i] < time[i-1]) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
    }

  return(state);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl *Random_Walk_Bridged_Generate(phydbl start, phydbl end, phydbl var, int n_steps)
{
  phydbl *state,end_walk;
  int i;

  if(n_steps == 0) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);

  state = Random_Walk_Generate(var,n_steps);
  end_walk = Rnorm(state[n_steps-1],SQRT(var));

  For(i,n_steps)
    {
      state[i] = state[i] - ((phydbl)(i+1.)/n_steps) * end_walk;
      state[i] = start + (end - start)/n_steps * (i+1.) + state[i];
    }
  
  return(state);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl *Random_Walk_Generate(phydbl var, int n_steps)
{
  phydbl *state;
  int i;

  if(n_steps == 0) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);

  state = (phydbl *)mCalloc(n_steps,sizeof(phydbl));

  state[0] = Rnorm(0.0,SQRT(var));

  for(i=1;i<n_steps;i++) state[i] = Rnorm(state[i-1],SQRT(var));

  return(state);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl Reflected(phydbl x, phydbl down, phydbl up)
{
  phydbl ref;
  
  ref = x;
  do
    {
      if(ref > up)   ref = up   - (ref  -  up);
      if(ref < down) ref = down + (down - ref);
    }while(!(ref < up && ref > down));

  return(ref);
}
