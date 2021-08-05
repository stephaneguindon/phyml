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
*     X = Gamma(k,1) = Erlang(k,1) = -sum(log(Ui)) = -log(prod(Ui))
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
      phydbl b = (alpha+exp(1.))/exp(1.);
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
	      y = exp(-x);
	    }
	  else 
	    {
	      x = -log(p*(b-v));
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
              x = -log((e - p) / a);
              if (Rexp(1.) >= (1.0 - a) * log(x))
                break;
	    } 
          else 
            {
              x = exp(log(p) / a);
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
	if (fabs(v) <= 0.25)
          q = q0 + 0.5 * t * t * ((((((a7 * v + a6) * v + a5) * v + a4) * v
                                    + a3) * v + a2) * v + a1) * v;
	else
          q = q0 - s * t + 0.25 * t * t + (s2 + s2) * log(1.0 + v);
        
        
	/* Step 7: quotient acceptance (q) */
	if (log(1.0 - u) <= q)
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
	    if (fabs(v) <= 0.25)
		q = q0 + 0.5 * t * t *
		    ((((((a7 * v + a6) * v + a5) * v + a4) * v + a3) * v
		      + a2) * v + a1) * v;
	    else
		q = q0 - s * t + 0.25 * t * t + (s2 + s2) * log(1.0 + v);
	    /* Step 11:	 hat acceptance (h) */
	    /* (if q not positive go to step 8) */
	    if (q > 0.0) {
		w = exp(q)-1.0;
		/* if t is rejected sample again at step 8 */
		if (c * fabs(u) <= w * exp(e - 0.5 * t * t))
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
  return -log(Uni()+SMALL)/lambda;
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

/* Returns a random deviates from an exponential distribution */
/* left-truncated at 'left' and right-truncated at 'rght' */
phydbl Rexp_Trunc(phydbl lambda, phydbl left, phydbl rght)
{
  phydbl u;
  u = Uni();
  return (left-log(1. - u*(1.-exp(-lambda*(rght-left))))/lambda);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl Rnorm(phydbl mean, phydbl sd)
{
  /* Box-Muller transformation */
  phydbl u1, u2, res;
  
  /* u1=Uni(); */
  /* u2=Uni(); */
  /* u1 = SQRT(-2.*log(u1))*COS(6.28318530717959f*u2); */

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
  u1 = x*SQRT((-2.*log(d))/d);

  res = u1*sd+mean;

  if(isnan(res) || isinf(res)) printf("\n. res=%f sd=%f mean=%f u1=%f u2=%f",res,sd,mean,u1,u2);

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

  for(i=0;i<dim;i++) x[i]=Rnorm(0.0,1.0);
  for(i=0;i<dim;i++) for(j=0;j<dim;j++) y[i] += L[i*dim+j]*x[j];
  for(i=0;i<dim;i++) y[i] += mu[i];

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
/* Borrowed from https://github.com/olafmersmann/truncnorm/blob/afc91b696db8a3feda25d39435fd979bacd962c6/src/rtruncnorm.c */

phydbl Rnorm_Trunc_Algo1(phydbl alpha, phydbl beta)
{
  phydbl z = -DBL_MAX;
  while(z < alpha || z > beta)
    {
      z = Rnorm(0.0,1.0);
    }
  return(z);
}

phydbl Rnorm_Trunc_Algo2(phydbl alpha, phydbl beta)
{
  phydbl z = 0.0;
  phydbl d_alpha = Dnorm(alpha,0.0,1.0);
  const double ub = alpha < 0.0 && beta > 0.0 ? M_1_SQRT_2PI : d_alpha;
  do
    {
      z = Uni()*(beta-alpha) + alpha;
    }
  while(Uni() * ub > Dnorm(z,0.0,1.0));
  return(z);
}

phydbl Rnorm_Trunc_Algo3(phydbl alpha, phydbl beta)
{

  phydbl z = alpha - 1.0;
  while(z < alpha || z > beta)
    {
      z = Rnorm(0,1);
      z = fabs(z);
    }
  return(z);
}

phydbl Rnorm_Trunc_Algo4(phydbl alpha, phydbl beta)
{
  phydbl z = 0.0;
  const phydbl ainv = 1.0/alpha;
  phydbl rho;
  do
    {
      z = Rexp(ainv) + alpha;
      rho = exp(-0.5 * pow((z-alpha),2));
    }
  while(Uni() > rho || z > beta);
  return(z);
}


phydbl Rnorm_Trunc(phydbl mean, phydbl sd, phydbl min, phydbl max, int *error)
{
  phydbl alpha,beta;
  phydbl d_alpha,d_beta;
  phydbl z,ret_val;
  
  alpha = (min - mean)/sd;
  beta  = (max - mean)/sd;

  d_alpha = Dnorm(alpha,0.0,1.0);
  d_beta  = Dnorm(beta,0.0,1.0);

  if(beta < alpha) return NAN;
  else
    {
      if(!(alpha > 0.0) && !(beta < 0.0))
        {
          if(!(d_alpha > 0.15) || !(d_beta > 0.15))
            {
              z = Rnorm_Trunc_Algo1(alpha,beta);
              return(mean + sd * z);
            }
          else
            {
              z = Rnorm_Trunc_Algo2(alpha,beta);
              return(mean + sd * z);
            }
        }
      else if(alpha > 0.0)
        {
          if(!(d_alpha / d_beta > 2.18))
            {
              z = Rnorm_Trunc_Algo2(alpha,beta);
              return(mean + sd * z);
            }
          else
            {
              if(!(alpha > 0.725))
                {
                  z = Rnorm_Trunc_Algo3(alpha,beta);
                  return(mean + sd * z);
                }
              else
                {
                  z = Rnorm_Trunc_Algo4(alpha,beta);
                  return(mean + sd * z);
                }
            }
        }
      else
        {
          if(!(d_beta / d_alpha > 2.18))
            {
              z = Rnorm_Trunc_Algo2(alpha,beta);
              return(mean - sd * z);
            }
          else if(beta > -0.725)
            {
              z = Rnorm_Trunc_Algo3(alpha,beta);
              return(mean - sd * z);
            }
          else
            {              
              z = Rnorm_Trunc_Algo4(alpha,beta);
              return(mean - sd * z);
            }
        }
    }
  
  ret_val = mean + sd*z;
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
      for(j=0;j<i;j++) rec += L[i*dim+j] * u[j];
      low  = (min[i]-mean[i]-rec)/L[i*dim+i];
      up   = (max[i]-mean[i]-rec)/L[i*dim+i];
      u[i] = Rnorm_Trunc(0.0,1.0,low,up,&err);
    }

  x = Matrix_Mult(L,u,dim,dim,dim,1);

/*   PhyML_Printf("\n>>>\n"); */
/*   for(i=0;i<dim;i++) */
/*     { */
/*       for(j=0;j<dim;j++) */
/* 	{ */
/* 	  PhyML_Printf("%10lf ",L[i*dim+j]); */
/* 	} */
/*       PhyML_Printf("\n"); */
/*     } */
/*   PhyML_Printf("\n"); */

/*   for(i=0;i<dim;i++) PhyML_Printf("%f ",u[i]); */
/*   PhyML_Printf("\n"); */

  
/*   PhyML_Printf("\n"); */
/*   for(i=0;i<dim;i++) PhyML_Printf("%10lf ",x[i]); */
/*   PhyML_Printf("\n<<<\n"); */

  for(i=0;i<dim;i++) x[i] += mean[i];

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
  x = log(u) / log(1. - p);
  
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
    prob = (k - 1.)*log(1. - p) + log(p);
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

  dens = 1./(SQRT(2*pi)*sd)*exp(-((x-mean)*(x-mean)/(2.*sd*sd)));

  return dens;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Dnorm(phydbl x, phydbl mean, phydbl sd)
{
  phydbl dens;

  /* dens = -(.5*LOG2PI+log(sd))  - .5*POW(x-mean,2)/POW(sd,2); */
  /* return exp(dens); */

  if(sd < SMALL && fabs(x-mean) < SMALL) return(1.0);
  
  x = (x-mean)/sd;

  dens = M_1_SQRT_2_PI * exp(-0.5*x*x);
  
  return dens / sd;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl Log_Dnorm(phydbl x, phydbl mean, phydbl sd, int *err)
{
  phydbl dens;

  *err = NO;

  if(sd < SMALL)
    {
      if(fabs(x-mean) < SMALL) return(0.0);
      else return(-INFINITY);
    }
  
  x = (x-mean)/sd;
  
  dens = -LOG_SQRT_2PI;
  dens -= pow(x,2)*.5;
  dens -= log(sd);
  
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

  if(sd < SMALL && fabs(x-mean) < SMALL) return(0.0);
  
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
      log_dens = -230.; /* ~log(1.E-100) */
    }
  else
    {
      log_dens -= log(cdf_up - cdf_lo);
    }

  if(isnan(log_dens) || isinf(fabs(log_dens)))
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

  if(sd < SMALL && fabs(x-mean) < SMALL) return(1.0);

  dens   = Dnorm(x,mean,sd);
  cdf_up = Pnorm(up,mean,sd);
  cdf_lo = Pnorm(lo,mean,sd);

  dens /= (cdf_up - cdf_lo);

  if(isnan(dens) || isinf(fabs(dens)))
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

  for(i=0;i<size;i++) xmmu[i] = x[i] - mu[i];
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

  density = size * LOG2PI + log(det) + buff2[0];
  density /= -2.;

/*   density = (1./(POW(2.*PI,size/2.)*SQRT(fabs(det)))) * exp(-0.5*buff2[0]); */

  Free(xmmu);
  Free(invcov);
  Free(buff1);
  Free(buff2);

  return (_log)?(density):(exp(density));
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

  for(i=0;i<size;i++) xmmu[i] = x[i] - mu[i];
  
  buff1 = Matrix_Mult(xmmu,invcov,1,size,size,size);
  buff2 = Matrix_Mult(buff1,xmmu,1,size,size,1);

  density = size * LOG2PI + log_det + buff2[0];
  density /= -2.;

  Free(xmmu);
  Free(buff1);
  Free(buff2);
  
  return (_log)?(density):(exp(density));
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
  dens *= exp((-1./(2.*(1.-rho2)))*(cx*cx/(sdx*sdx)+cy*cy/(sdy*sdy)+2*rho*cx*cy/(sdx*sdy)));
	      
  return dens;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// E(X) = n(1-p)/p ; V(X) = n(1-p)/p^2
phydbl Dnbinom(phydbl x, phydbl n, phydbl p, int logit)
{

  phydbl lnDens = LnGamma(x+n) - LnGamma(n) - LnFact(x) + n*log(p) + x*log(1.-p);

  if(logit == TRUE)
    {
      return(lnDens);
    }
  else
    {
      return(exp(lnDens));
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

  if(x < 1.E-200)
    {
      if(x < 0.0) return 0.0;
      else
	{
	  PhyML_Printf("\n. WARNING: Dgamma -> small value of x = %G (shape: %G scale: %G <X>: %G)",x,shape,scale,shape*scale);
	  x = 1.E-200;
	}
    }


  if(scale < 0.0 || shape < 0.0)
    {
      PhyML_Printf("\n. scale=%f shape=%f",scale,shape);
      Exit("\n");
    }


  v = (shape-1.) * log(x) - shape * log(scale) - x / scale - LnGamma(shape);


  if(v < 500.)
    {
      v = exp(v);
    }
  else
    {
      PhyML_Printf("\n. WARNING v=%f x=%f shape=%f scale=%f",v,x,shape,scale);
      PhyML_Printf("\n. log(x) = %G LnGamma(shape)=%G",log(x),LnGamma(shape));
      v = exp(v);
      /* PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__); */
      /* Exit("\n"); */
    }

	 
  return v;
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl Dexp(phydbl x, phydbl param)
{
  return param * exp(-param * x);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

/* Returns the density of an exponential distribution left-truncated */
/* at 'left' and right-truncated at 'rght' */

phydbl Dexp_Trunc(phydbl x, phydbl lambda, phydbl left, phydbl rght)
{
  return (lambda * exp(-lambda * x))/(exp(-lambda * left) - exp(-lambda * rght));
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
/* Poisson probability */
phydbl Dpois(phydbl x, phydbl param, int logit)
{
  phydbl v;

  if(param < SMALL)
    {
      if(x < SMALL)
        {
          if(logit) return 0.0;
          else return 1.0;
        }
      else
        {
          if(logit) return -INFINITY;
          else return 0.0;
        }
    }

  if(x < .0) 
    {
      if(logit == YES) return(-INFINITY);
      else return 0.0;
   }

  v = x * log(param) - param - Factln(x);
  if(logit == YES) return v;
  else
    {
      if(v < 500.)
        {
          v = exp(v);
        }
      else
        {
          PhyML_Printf("\n. WARNING v=%f x=%f param=%f",v,x,param);
          v = exp(500);
        }  
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
// Error function
phydbl Erf(phydbl x)
{
  return(2.*Pnorm(x*sqrt(2.),0.0,1.0)-1.);
}

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
/*       return (1.0 - c * exp( -x * x / 2.0 ) * t * */
/* 	      ( t *( t * ( t * ( t * b5 + b4 ) + b3 ) + b2 ) + b1 )); */
/*     } */
/*   else */
/*     { */
/*       phydbl t = 1.0 / ( 1.0 - p * x ); */
/*       return ( c * exp( -x * x / 2.0 ) * t * */
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
   if (fabs(q/ch-1) > e) goto l4;

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
      y = SQRT (log(1/(p1*p1)));   
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
  return FLOOR(0.5+exp(Factln(n)-Factln(k)-Factln(n-k)));
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl Factln(int n)
{
  /* static phydbl a[101]; */
  
  /* if (n < 0)    { Warn_And_Exit("\n== Err: negative factorial in routine FACTLN"); } */
  /* if (n <= 1)     return 0.0; */
  /* if (n <= 100)   return (a[n]>SMALL) ? a[n] : (a[n]=Gammln(n+1.0)); */
  /* else return     Gammln(n+1.0); */
  return Gammln(n+1.0);
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
   else
     {
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


/* Return log(n!) */

phydbl LnFact(int n)
{
  int i;
  phydbl res;

  res = .0;
  for(i=2;i<n+1;i++) res += log((phydbl)i);
  
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
  
  for(i=0;i<tree->data->crunch_len;i++) ori_wght[i] = tree->data->wght[i];

  n_site = 0;
  for(i=0;i<tree->data->crunch_len;i++) For(j,tree->data->wght[i])
    {
      site_num[n_site] = i;
      n_site++;
    }

  
  tree->verbose = VL0;
  for(replicate=0;replicate<sample_size;replicate++)
    {
      For(i,2*tree->n_otu-3) tree->a_edges[i]->l->v = .1;

      for(i=0;i<tree->data->crunch_len;i++) tree->data->wght[i] = 0;

      for(i=0;i<tree->data->init_len;i++)
	{
	  position = Rand_Int(0,(int)(tree->data->init_len-1.0));
	  tree->data->wght[site_num[position]] += 1;
	}

      Round_Optimize(tree,ROUND_MAX);
      
      For(i,2*tree->n_otu-3) For(j,2*tree->n_otu-3) cov[i*dim+j] += log(tree->a_edges[i]->l->v) * log(tree->a_edges[j]->l->v);  
      For(i,2*tree->n_otu-3) mean[i] += log(tree->a_edges[i]->l->v);

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

  for(i=0;i<tree->data->crunch_len;i++) tree->data->wght[i] = ori_wght[i];

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


  for(i=0;i<dim;i++) ori_bl[i] = tree->a_edges[i]->l->v;

  if(tree->mod->log_l == NO)
    l_inf = MAX(tree->mod->l_min,1./(phydbl)tree->data->init_len);
  else
    l_inf = MAX(tree->mod->l_min,-log((phydbl)tree->data->init_len));


  n_ok_edges = 0;
  for(i=0;i<dim;i++) 
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
  for(i=0;i<dim;i++)
    {
      do
	{
	  tree->a_edges[i]->l->v += inc[i];
	  lnL1 = Lk(tree->a_edges[i],tree);
	  tree->a_edges[i]->l->v = ori_bl[i];
	  inc[i] *= 1.1;
	}while((fabs(lnL1 - ori_lnL) < 1.E-1) && 
	       (tree->a_edges[i]->l->v+inc[i] < tree->mod->l_max));
      inc[i] /= 1.1;
    }



  /* zero zero */  
  zero_zero = tree->c_lnL;

  /* plus zero */  
  for(i=0;i<dim;i++) 
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
  for(i=0;i<dim;i++) 
    {
      if(is_ok[i])
	{
	  tree->a_edges[i]->l->v -= inc[i];
	  lk = Lk(tree->a_edges[i],tree);
	  minus_zero[i] = lk;
	  tree->a_edges[i]->l->v = ori_bl[i];
	}
    }


  for(i=0;i<dim;i++) Update_PMat_At_Given_Edge(tree->a_edges[i],tree);

  /* plus plus  */  
  for(i=0;i<dim;i++)
    {
      if(is_ok[i])
	{
	  tree->a_edges[i]->l->v += inc[i];
	  Update_PMat_At_Given_Edge(tree->a_edges[i],tree);

	  for(j=0;j<3;j++)
	    if((!tree->a_edges[i]->left->tax) && (tree->a_edges[i]->left->v[j] != tree->a_edges[i]->rght))
	      Recurr_Hessian(tree->a_edges[i]->left,tree->a_edges[i]->left->v[j],1,inc,plus_plus+i*dim,is_ok,tree);
	  
	  for(j=0;j<3;j++)
	    if((!tree->a_edges[i]->rght->tax) && (tree->a_edges[i]->rght->v[j] != tree->a_edges[i]->left))
	      Recurr_Hessian(tree->a_edges[i]->rght,tree->a_edges[i]->rght->v[j],1,inc,plus_plus+i*dim,is_ok,tree);
		      
	  tree->a_edges[i]->l->v = ori_bl[i];
	  Lk(NULL,tree);
	}
    }


  /* plus minus */  
  for(i=0;i<dim;i++)
    {
      if(is_ok[i])
	{
	  tree->a_edges[i]->l->v += inc[i];
	  Update_PMat_At_Given_Edge(tree->a_edges[i],tree);
	  
	  for(j=0;j<3;j++)
	    if((!tree->a_edges[i]->left->tax) && (tree->a_edges[i]->left->v[j] != tree->a_edges[i]->rght))
	      Recurr_Hessian(tree->a_edges[i]->left,tree->a_edges[i]->left->v[j],-1,inc,plus_minus+i*dim,is_ok,tree);
	  
	  for(j=0;j<3;j++)
	    if((!tree->a_edges[i]->rght->tax) && (tree->a_edges[i]->rght->v[j] != tree->a_edges[i]->left))
	      Recurr_Hessian(tree->a_edges[i]->rght,tree->a_edges[i]->rght->v[j],-1,inc,plus_minus+i*dim,is_ok,tree);
	  
	  tree->a_edges[i]->l->v = ori_bl[i];
	  Lk(NULL,tree);
	}
    }


  /* minus minus */  
  for(i=0;i<dim;i++)
    {
      if(is_ok[i])
	{
	  tree->a_edges[i]->l->v -= inc[i];
	  
	  Update_PMat_At_Given_Edge(tree->a_edges[i],tree);
	  
	  for(j=0;j<3;j++)
	    if((!tree->a_edges[i]->left->tax) && (tree->a_edges[i]->left->v[j] != tree->a_edges[i]->rght))
	      Recurr_Hessian(tree->a_edges[i]->left,tree->a_edges[i]->left->v[j],-1,inc,minus_minus+i*dim,is_ok,tree);
	  
	  for(j=0;j<3;j++)
	    if((!tree->a_edges[i]->rght->tax) && (tree->a_edges[i]->rght->v[j] != tree->a_edges[i]->left))
	      Recurr_Hessian(tree->a_edges[i]->rght,tree->a_edges[i]->rght->v[j],-1,inc,minus_minus+i*dim,is_ok,tree);
	  
	  tree->a_edges[i]->l->v = ori_bl[i];
	  Lk(NULL,tree);
	}
    }


  
  for(i=0;i<dim;i++)
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
        
  for(i=0;i<n_ok_edges;i++)
    {
      for(j=0;j<n_ok_edges;j++)
	{
	  buff[i*n_ok_edges+j] = -1.0*hessian[ok_edges[i]*dim+ok_edges[j]];
	}
    }


  if(!Matinv(buff,n_ok_edges,n_ok_edges,NO))
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");      
    }

  for(i=0;i<n_ok_edges;i++)
    {
      for(j=0;j<n_ok_edges;j++)
	{
	  hessian[ok_edges[i]*dim+ok_edges[j]] = buff[i*n_ok_edges+j];
	}
    }

/*   eps = 1./(phydbl)tree->data->init_len; */
  /* Approximate variance for very short branches */
  for(i=0;i<dim;i++)
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
	  }while(fabs(lnL2 - lnL) < 1.E-3);

	hessian[i*dim+i] = -1.0 / hessian[i*dim+i];

      }
  

  /* Fit a straight line to the log-likelihood (i.e., an exponential to the likelihood) */
  /* It is only a straight line when considering branch length (rather than log(branch lengths)) */
  for(i=0;i<dim;i++)
    if((tree->a_edges[i]->l->v / tree->mod->l_min < 1.1) &&
       (tree->a_edges[i]->l->v / tree->mod->l_min > 0.9))
      {
	phydbl *x,*y,l;
	phydbl cov,var;
	
	x=plus_plus;
	y=minus_minus;
	l=(tree->mod->log_l == YES)?(exp(tree->a_edges[i]->l->v)):(tree->a_edges[i]->l->v); /* Get actual branch length */
	
	for(j=0;j<dim;j++)
	  {
	    x[j] = l + (100.*l-l)*((phydbl)j/dim);
	    tree->a_edges[i]->l->v = (tree->mod->log_l)?(log(x[j])):(x[j]); /* Transform to log if necessary */
	    y[j] = Lk(tree->a_edges[i],tree);
	    tree->a_edges[i]->l->v = (tree->mod->log_l)?(log(l)):(l); /* Go back to initial edge length */
	  }
	
	cov = Covariance(x,y,dim);
	var = Covariance(x,x,dim);
	
	/* cov/var is minus the parameter of the exponential distribution.
	   The variance is therefore : */
	hessian[i*dim+i] = 1.0 / pow(cov/var,2);
	
	/* 	    printf("\n. Hessian = %G cov=%G var=%G",hessian[i*dim+i],cov,var); */
      }
  /*     } */


  for(i=0;i<dim;i++)
    if(hessian[i*dim+i] < 0.0)
      {
	PhyML_Printf("\n. l=%G var=%G",tree->a_edges[i]->l->v,hessian[i*dim+i]);
/* 	PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__); */
/* 	Exit("\n"); */
	hessian[i*dim+i] = MIN_VAR_BL;
      }

  for(i=0;i<dim;i++)
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

  for(i=0;i<dim;i++)
    {
      for(j=0;j<dim;j++)
	{
	  if(fabs(hessian[i*dim+j]-hessian[j*dim+i]) > 1.E-3)
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
/*   for(i=0;i<dim;i++) */
/*     { */
/*       PhyML_Printf("[%f] ",tree->a_edges[i]->l->v); */
/*       for(j=0;j<dim;j++) */
/* 	{ */
/* 	  PhyML_Printf("%12lf ",hessian[i*dim+j]); */
/* 	} */
/*       PhyML_Printf("\n"); */
/*     } */

  /* Matinv(hessian,dim,dim,NO); */

  /* PhyML_Printf("\n"); */

  /* for(i=0;i<dim;i++) */
  /*   { */
  /*     PhyML_Printf("[%f] ",tree->a_edges[i]->l->v); */
  /*     for(j=0;j<dim;j++) */
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

  for(i=0;i<dim;i++) ori_bl[i] = tree->a_edges[i]->l->v;

  if(tree->mod->log_l == NO)
    l_inf = MAX(tree->mod->l_min,1./(phydbl)tree->data->init_len);
  else
    l_inf = MAX(tree->mod->l_min,-log((phydbl)tree->data->init_len));

  for(i=0;i<dim;i++) 
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
  for(i=0;i<dim;i++) 
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
  for(i=0;i<dim;i++)
    {
      if(is_ok[i] == YES)
	{
	  tree->a_edges[i]->l->v -= inc[i];
	  lk = Lk(tree->a_edges[i],tree);
	  minus[i] = lk;
	  tree->a_edges[i]->l->v = ori_bl[i];
	}
    }


  for(i=0;i<dim;i++)
    {
      if(is_ok[i] == YES)
	{
	  gradient[i] = (plus[i] - minus[i])/(2.*inc[i]);
	}
    }


  for(i=0;i<dim;i++)
    {
      if(is_ok[i] == NO)
	{
	  eps = fabs(0.2 * tree->a_edges[i]->l->v);
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
/*   for(i=0;i<dim;i++) */
/*     { */
/*       PhyML_Printf("[%f] ",tree->a_edges[i]->l->v); */
/*       for(j=0;j<dim;j++) */
/* 	{ */
/* 	  printf("%12lf ",gradient[i]*gradient[j]); */
/* 	} */
/*       printf("\n"); */
/*     } */
/*   printf("\n"); */
/*   for(i=0;i<dim;i++) */
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
  
  for(i=0;i<dim;i++) ori_bl[i] = tree->a_edges[i]->l->v;
  
  Set_Both_Sides(YES,tree);
  Lk(NULL,tree);
  ori_lnL = tree->c_lnL;


  if(tree->mod->log_l == NO)
    l_inf = MAX(tree->mod->l_min,1./(phydbl)tree->data->init_len);
  else
    l_inf = MAX(tree->mod->l_min,-log((phydbl)tree->data->init_len));

  for(i=0;i<dim;i++) 
    {
      if(tree->a_edges[i]->l->v*(1.-eps) > l_inf)
	{
	  inc_plus[i]  = fabs(eps * MAX(tree->mod->l_min,tree->a_edges[i]->l->v));
	  inc_minus[i] = fabs(eps * MAX(tree->mod->l_min,tree->a_edges[i]->l->v));
	  is_ok[i]     = YES;
	}
      else
	{
	  inc_plus[i]  = fabs(0.2 * MAX(tree->mod->l_min,tree->a_edges[i]->l->v));
	  inc_minus[i] = fabs(0.2 * MAX(tree->mod->l_min,tree->a_edges[i]->l->v));
	  is_ok[i]     = NO;
    }


    }

  /* Fine tune the increments */
  for(i=0;i<dim;i++)
    {
      do
	{
	  tree->a_edges[i]->l->v += inc_plus[i];	  
	  lnL1 = Lk(tree->a_edges[i],tree);
	  tree->a_edges[i]->l->v = ori_bl[i];
	  inc_plus[i] *= 1.1;
	}while((fabs(lnL1 - ori_lnL) < 1.E-1) && (tree->a_edges[i]->l->v+inc_plus[i] < tree->mod->l_max));
      inc_plus[i] /= 1.1;
    }

  for(i=0;i<dim;i++)
    {
      do
	{
	  tree->a_edges[i]->l->v -= inc_minus[i];
	  lnL1 = Lk(tree->a_edges[i],tree);
	  tree->a_edges[i]->l->v = ori_bl[i];
	  inc_minus[i] *= 1.1;
	}while((fabs(lnL1 - ori_lnL) < 1.E-1) && 
	       (tree->a_edges[i]->l->v -inc_minus[i] > tree->mod->l_min));
      inc_minus[i] /= 1.1;
    }

  for(i=0;i<dim;i++) 
    {
      inc[i] = MIN(inc_plus[i],inc_minus[i]);
    }

  /* plus */  
  for(i=0;i<dim;i++) 
    {
      if(is_ok[i] == YES)
	{
	  tree->a_edges[i]->l->v += inc[i];
	  Lk(tree->a_edges[i],tree);
	  for(j=0;j<tree->n_pattern;j++) plus[i*tree->n_pattern+j] = log(tree->cur_site_lk[j]);
	  tree->a_edges[i]->l->v = ori_bl[i];
	}
    }


  /* minus */
  for(i=0;i<dim;i++)
    {
      if(is_ok[i] == YES)
	{
	  tree->a_edges[i]->l->v -= inc[i];
	  Lk(tree->a_edges[i],tree);
	  for(j=0;j<tree->n_pattern;j++) minus[i*tree->n_pattern+j] = log(tree->cur_site_lk[j]);
	  tree->a_edges[i]->l->v = ori_bl[i];
	}
    }


  for(i=0;i<dim;i++)
    {
      if(is_ok[i] == NO)
	{
	  Lk(tree->a_edges[i],tree);	
	  for(j=0;j<tree->n_pattern;j++) zero[i*tree->n_pattern+j] = log(tree->cur_site_lk[j]);
	  
	  tree->a_edges[i]->l->v += inc[i];
	  lnL1 = Lk(tree->a_edges[i],tree);
	  for(j=0;j<tree->n_pattern;j++) plus[i*tree->n_pattern+j] = log(tree->cur_site_lk[j]);
	  
	  tree->a_edges[i]->l->v += inc[i];
	  lnL2 = Lk(tree->a_edges[i],tree);
	  for(j=0;j<tree->n_pattern;j++) plusplus[i*tree->n_pattern+j] = log(tree->cur_site_lk[j]);
	  
	  tree->a_edges[i]->l->v = ori_bl[i];	

	}
    }

  For(i,dim*dim) hessian[i] = 0.0;

  for(k=0;k<tree->n_pattern;k++)
    {
      for(i=0;i<dim;i++) 
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
      for(i=0;i<dim;i++) for(j=0;j<dim;j++) site_hessian[i*dim+j] = gradient[i] * gradient[j];
      For(i,dim*dim) hessian[i] -= site_hessian[i] * tree->data->wght[k]; 
    }


  /* Make sure to update likelihood before bailing out */
  Set_Both_Sides(YES,tree);
  Lk(NULL,tree);

  l = tree->data->init_len;
  n = tree->mod->ns;
  /* Delta method for variance. Assumes Jukes and Cantor with p=1/n */
  small_var = (1./(l*l))*(1.-1./l)*(n-1.)*(n-1.)/(n-1.-n/l);
  for(i=0;i<dim;i++)
    if(is_ok[i] == NO)
      {
	for(j=0;j<dim;j++)
	  {
	    hessian[i*dim+j] = 0.;
	    hessian[j*dim+i] = 0.;
	  }
	hessian[i*dim+i] = -1./small_var;

	if(tree->mod->log_l == YES) 
	  {
	    hessian[i*dim+i] = small_var * POW(exp(tree->a_edges[i]->l->v),-2); 
	    hessian[i*dim+i] = -1./hessian[i*dim+i];
	  }
      }

  for(i=0;i<dim;i++)
    if(is_ok[i] == YES && hessian[i*dim+i] < -1./small_var)
      {
  	for(j=0;j<dim;j++)
  	  {
  	    hessian[i*dim+j] = 0.;
  	    hessian[j*dim+i] = 0.;
  	  }
  	hessian[i*dim+i] = -1./small_var;

  	if(tree->mod->log_l == YES)
  	  {
  	    hessian[i*dim+i] = small_var * POW(exp(tree->a_edges[i]->l->v),-2);
  	    hessian[i*dim+i] = -1./hessian[i*dim+i];
  	  }
      }


  for(i=0;i<dim;i++)
    {
      for(j=0;j<dim;j++)
	{
	  if(fabs(hessian[i*dim+j]-hessian[j*dim+i]) > 1.E-3)
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
  /* for(i=0;i<dim;i++) */
  /*   { */
  /*     PhyML_Printf("[%f] ",tree->a_edges[i]->l->v); */
  /*     for(j=0;j<dim;j++) */
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

  for(i=0;i<3;i++)
    if(a->v[i] == d)
      {
	Update_Partial_Lk(tree,a->b[i],a);

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
    for(i=0;i<3;i++) 
      if(d->v[i] != a) 
	Recurr_Hessian(d,d->v[i],plus_minus,inc,res,is_ok,tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


/* Work out the Hessian for the likelihood function. Only logARITHM of branch lengths are considered as variable.
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

  for(i=0;i<dim;i++) ori_bl[i] = tree->a_edges[i]->l->v;

  n_ok_edges = 0;
  for(i=0;i<dim;i++) 
    {
      if(tree->a_edges[i]->l->v > 3.0/(phydbl)tree->data->init_len)
	{
	  inc[i] = fabs(eps * tree->a_edges[i]->l->v);
	  ok_edges[n_ok_edges] = i;
	  n_ok_edges++;
	  is_ok[i] = 1;
	}
      else is_ok[i] = 0;
    }

  /* zero zero */  
  lk = Log_Det(is_ok,tree);
  for(i=0;i<dim;i++) if(is_ok[i]) zero_zero[i] = tree->c_lnL+lk;

  /* plus zero */  
  for(i=0;i<dim;i++) 
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
  for(i=0;i<dim;i++) 
    {
      if(is_ok[i])
	{
	  tree->a_edges[i]->l->v -= inc[i];
	  lk = Lk(tree->a_edges[i],tree);
	  minus_zero[i] = lk+Log_Det(is_ok,tree);
	  tree->a_edges[i]->l->v = ori_bl[i];
	}
    }

  for(i=0;i<dim;i++) Update_PMat_At_Given_Edge(tree->a_edges[i],tree);

  /* plus plus  */  
  for(i=0;i<dim;i++)
    {
      if(is_ok[i])
	{
	  tree->a_edges[i]->l->v += inc[i];
	  Update_PMat_At_Given_Edge(tree->a_edges[i],tree);

	  for(j=0;j<3;j++)
	    if((!tree->a_edges[i]->left->tax) && (tree->a_edges[i]->left->v[j] != tree->a_edges[i]->rght))
	      Recurr_Hessian_Log(tree->a_edges[i]->left,tree->a_edges[i]->left->v[j],1,inc,plus_plus+i*dim,is_ok,tree);
	  
	  for(j=0;j<3;j++)
	    if((!tree->a_edges[i]->rght->tax) && (tree->a_edges[i]->rght->v[j] != tree->a_edges[i]->left))
	      Recurr_Hessian_Log(tree->a_edges[i]->rght,tree->a_edges[i]->rght->v[j],1,inc,plus_plus+i*dim,is_ok,tree);

/* 	  for(j=0;j<dim;j++)  */
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
  for(i=0;i<dim;i++)
    {
      if(is_ok[i])
	{
	  tree->a_edges[i]->l->v += inc[i];
	  Update_PMat_At_Given_Edge(tree->a_edges[i],tree);
	  
	  for(j=0;j<3;j++)
	    if((!tree->a_edges[i]->left->tax) && (tree->a_edges[i]->left->v[j] != tree->a_edges[i]->rght))
	      Recurr_Hessian_Log(tree->a_edges[i]->left,tree->a_edges[i]->left->v[j],-1,inc,plus_minus+i*dim,is_ok,tree);
	  
	  for(j=0;j<3;j++)
	    if((!tree->a_edges[i]->rght->tax) && (tree->a_edges[i]->rght->v[j] != tree->a_edges[i]->left))
	      Recurr_Hessian_Log(tree->a_edges[i]->rght,tree->a_edges[i]->rght->v[j],-1,inc,plus_minus+i*dim,is_ok,tree);
	  
/* 	  for(j=0;j<dim;j++)  */
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
  for(i=0;i<dim;i++)
    {
      if(is_ok[i])
	{
	  tree->a_edges[i]->l->v -= inc[i];
	  
	  Update_PMat_At_Given_Edge(tree->a_edges[i],tree);
	  
	  for(j=0;j<3;j++)
	    if((!tree->a_edges[i]->left->tax) && (tree->a_edges[i]->left->v[j] != tree->a_edges[i]->rght))
	      Recurr_Hessian_Log(tree->a_edges[i]->left,tree->a_edges[i]->left->v[j],-1,inc,minus_minus+i*dim,is_ok,tree);
	  
	  for(j=0;j<3;j++)
	    if((!tree->a_edges[i]->rght->tax) && (tree->a_edges[i]->rght->v[j] != tree->a_edges[i]->left))
	      Recurr_Hessian_Log(tree->a_edges[i]->rght,tree->a_edges[i]->rght->v[j],-1,inc,minus_minus+i*dim,is_ok,tree);
	  
/* 	  for(j=0;j<dim;j++)  */
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

/*   for(i=0;i<dim;i++) if(is_ok[i]) inc[i] = POW(tree->a_edges[i]->l->v+inc[i],2)-POW(tree->a_edges[i]->l->v,2); */
  for(i=0;i<dim;i++) if(is_ok[i]) inc[i] = log(tree->a_edges[i]->l->v+inc[i])-log(tree->a_edges[i]->l->v);
/*   for(i=0;i<dim;i++) inc[i] = 2.*inc[i]; */
/*   for(i=0;i<dim;i++) if(is_ok[i]) inc[i] = SQRT(tree->a_edges[i]->l->v+inc[i])-SQRT(tree->a_edges[i]->l->v); */
  
  for(i=0;i<dim;i++)
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
        

  for(i=0;i<n_ok_edges;i++)
    {
      for(j=0;j<n_ok_edges;j++)
	{
	  buff[i*n_ok_edges+j] = -hessian[ok_edges[i]*dim+ok_edges[j]];
	}
    }

  if(!Matinv(buff,n_ok_edges,n_ok_edges,NO))
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");      
    }

  for(i=0;i<n_ok_edges;i++)
    {
      for(j=0;j<n_ok_edges;j++)
	{
	  hessian[ok_edges[i]*dim+ok_edges[j]] = buff[i*n_ok_edges+j];
	}
    }

  /* Approximate variance for very short branches */
  for(i=0;i<dim;i++)
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

/*   for(i=0;i<dim;i++) */
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

  for(i=0;i<dim;i++)
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

  for(i=0;i<3;i++)
    if(a->v[i] == d)
      {
	Update_Partial_Lk(tree,a->b[i],a);

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
    for(i=0;i<3;i++) 
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
/*   For(i,2*tree->n_otu-3) if(is_ok[i]) ldet += log(2.*SQRT(tree->a_edges[i]->l->v)); */
  For(i,2*tree->n_otu-3) if(is_ok[i]) ldet += log(tree->a_edges[i]->l->v);
/*   For(i,2*tree->n_otu-3) if(is_ok[i]) ldet -= log(2*tree->a_edges[i]->l->v); */
  
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

  for(j=0;j<100;j++) 
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


int Matinv(phydbl *x, const int n, const int m, const int verbose)
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

   for(i=0;i<n;++i)
     {
       xmax = 0.;
       for (j=i; j<n; j++)
         if (xmax < fabs(x[j*m+i]))
	   {
	     xmax = fabs(x[j*m+i]);
	     irow[i]=j;
	   }

      det *= xmax;
      if (xmax < ee)
	{
	  Free(irow);
	  if(verbose)
	    {
	      PhyML_Printf("\n. xmax=%g",xmax);
	      PhyML_Printf("\n. Determinant becomes zero at %3d!\t", i+1);
	      PhyML_Printf("\n. Failed to invert the matrix.");
	    }
	  return(0);
	}
      if (irow[i] != i)
	{
	  for (j=0;j<m;++j)
	    {
	      t = x[i*m+j];
	      x[i*m+j] = x[irow[i]*m+j];
	      x[irow[i]*m+j] = t;
	    }
	}
      t = 1./x[i*m+i];
      for (j=0;j<n;++j)
	{
	  if (j == i) continue;
	  t1 = t*x[j*m+i];
	  for(k=0;k<m;k++)  x[j*m+k] -= t1*x[i*m+k];
	  x[j*m+i] = -t1;
	}
      for(j=0;j<m;j++)   x[i*m+j] *= t;
      x[i*m+i] = t;
   }                            /* i  */
   for (i=n-1; i>=0; i--)
     {
       if (irow[i] == i) continue;
       for(j=0;j<n;j++)
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
  
  for(i=0;i<nra;i++)
    for(j=0;j<ncb;j++)
       for(k=0;k<nca;k++)
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

  for(i=0;i<dim;i++) for(j=i+1;j<dim;j++) 
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
  for(i=0;i<size;i++) det += log(triA[i*size+i]);
  Free(triA);
 
  if(_log == NO)
    {
      det = exp(det);
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
  for(i=0;i<n;i++) { if(!is_1[i]) { ctrd_a[nr] = a[i]-mu[i]; nr++; } }

  nr=0;
  for(i=0;i<n;i++) { if( is_1[i]) { mu1[nr] = mu[i]; nr++; } }

  nr=0;
  for(i=0;i<n;i++) { if(!is_1[i]) { mu2[nr] = mu[i]; nr++; } }

  nr=0; nc=0;
  for(i=0;i<n;i++)
    {
      if(is_1[i])
	{
	  nc = nr;
 	  for(j=i;j<n;j++)
/* 	  nc = 0; */
/* 	  for(j=0;j<n;j++) */
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
  for(i=0;i<n;i++)
    {
      if(is_1[i])
	{
/* 	  nc = nr; */
/*  	  for(j=i;j<n;j++) */
	  nc = 0;
	  for(j=0;j<n;j++)
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
  for(i=0;i<n;i++)
    {
      if(!is_1[i])
	{
/* 	  nc = nr; */
/* 	  for(j=i;j<n;j++) */
	  nc = 0;
	  for(j=0;j<n;j++)
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
  for(i=0;i<n;i++)
    {
      if(!is_1[i])
	{
	  nc = nr;
	  for(j=i;j<n;j++)
/* 	  nc = 0; */
/* 	  for(j=0;j<n;j++) */
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
  for(i=0;i<n1;i++) cond_mu_norder[i] = mu1[i]+buff[i];
  Free(buff);

  buff = Matrix_Mult(sig12_invsig22,sig21,n1,n2,n2,n1);
  for(i=0;i<n1;i++) for(j=0;j<n1;j++) cond_cov_norder[i*n1+j] = sig11[i*n1+j] - buff[i*n1+j];
  Free(buff);

  nr = 0;
  for(i=0;i<n;i++) if(is_1[i]) { cond_mu[i] = cond_mu_norder[nr]; nr++; }

  nr = nc = 0;
  for(i=0;i<n;i++) 
    {
      if(is_1[i]) 
	{ 
	  nc = 0;
	  for(j=0;j<n;j++)
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

/*   for(i=0;i<n1;i++) */
/*     { */
/*       for(j=i;j<n1;j++) */
/* 	if(fabs(cond_cov_norder[i*n1+j] - cond_cov_norder[j*n1+i]) > 1.E-3) */
/* 	  { */
/* 	    PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__); */
/* 	    Warn_And_Exit(""); */
/* 	  } */
/*     } */


  for(i=0;i<n;i++)
    {
      for(j=i+1;j<n;j++)
	if(fabs(cond_cov[i*n+j] - cond_cov[j*n+i]) > 1.E-3)
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
  for(i=0;i<n;i++) { if(!is_1[i]) { ctrd_a[nr] = a[i]-mu[i]; nr++; } }

  nr=0;
  for(i=0;i<n;i++) { if( is_1[i]) { mu1[nr] = mu[i]; nr++; } }

  nr=0;
  for(i=0;i<n;i++) { if(!is_1[i]) { mu2[nr] = mu[i]; nr++; } }

  nr=0; nc=0;
  for(i=0;i<n;i++)
    {
      if(is_1[i])
	{
	  nc = nr;
 	  for(j=i;j<n;j++)
/* 	  nc = 0; */
/* 	  for(j=0;j<n;j++) */
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
  for(i=0;i<n;i++)
    {
      if(is_1[i])
	{
/* 	  nc = nr; */
/*  	  for(j=i;j<n;j++) */
	  nc = 0;
	  for(j=0;j<n;j++)
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
  for(i=0;i<n;i++)
    {
      if(!is_1[i])
	{
/* 	  nc = nr; */
/* 	  for(j=i;j<n;j++) */
	  nc = 0;
	  for(j=0;j<n;j++)
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
  for(i=0;i<n;i++)
    {
      if(!is_1[i])
	{
	  nc = nr;
	  for(j=i;j<n;j++)
/* 	  nc = 0; */
/* 	  for(j=0;j<n;j++) */
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
  for(i=0;i<n1;i++) cond_mu[i] = mu1[i]+buff[i];
  Free(buff);

  buff = Matrix_Mult(sig12_invsig22,sig21,n1,n2,n2,n1);
  for(i=0;i<n1;i++) for(j=0;j<n1;j++) cond_cov[i*n1+j] = sig11[i*n1+j] - buff[i*n1+j];


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
  for(i=0;i<n;i++)
    {
      if(is_1[i])
	{
	  nc = 0;
	  for(j=0;j<n;j++)
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
  for(i=0;i<n;i++)
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


  for(i=0;i<n;i++) reg_coeff[i] = 0.0;

/*   nr = 0; */
/*   for(i=0;i<n;i++) if(!is_1[i]) { reg_coeff[i] = sig12_invsig22[nr]; nr++; } */

  nc = 0;
  nr = 0;
  for(i=0;i<n1;i++) 
    {
      nc = 0;
      for(j=0;j<n;j++)
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
  cons = exp(cons);
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
  for(i=0;i<n;i++) mean_x += x[i];
  mean_x /= (phydbl)n;

  mean_y = .0;
  for(i=0;i<n;i++) mean_y += y[i];
  mean_y /= (phydbl)n;

  mean_xy = .0;
  for(i=0;i<n;i++) mean_xy += x[i]*y[i];
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
      for(i=0;i<len;i++)
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
	      for(j=0;j<len;j++) if(j != cond && j != i) alpha -= lambda[j] * x[j];
	      
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

  scaled_var = brownian_var/fabs(x_end - x_beg);

  n_breaks = 100;


  /* n_rep    = 500;   */

  /* x_step   = (x_end - x_beg)/(n_breaks+1); */

  /* y      = (phydbl *)mCalloc(n_breaks+2,sizeof(phydbl)); */
  /* y_mean = (phydbl *)mCalloc(n_rep,sizeof(phydbl)); */

  /* y[0] = y_beg; */
  /* y[n_breaks+1] = y_end; */

  /* for(i=0;i<n_rep;i++) */
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
  /*     For(j,n_breaks+2) sum += fabs(y[j]); */
  /*     y_mean[i] = sum/(n_breaks+2); */
  /*   } */

  /* sum = sumsum = 0.0; */
  /* for(i=0;i<n_rep;i++) */
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
	(2.*six)/SQRT(2.*PI)*exp(-POW(mux,2)/(2.*POW(six,2))) + 
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
at value 1 at time t=T. X'(t) starts at X'(t)=A at t=0 and stops at X'(t)=B at
t=T. This function calculates the mean and variance of 
Z(T) = 1/T \int_0^T exp(X'(t)) dt. 
*/

void Integrated_Geometric_Brownian_Bridge_Moments(phydbl T, phydbl A, phydbl B, phydbl u, phydbl *mean, phydbl *var)
{
  Integrated_Geometric_Brownian_Bridge_Mean(T,A,B,u,mean);
  Integrated_Geometric_Brownian_Bridge_Var(T,A,B,u,var);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Integrated_Geometric_Brownian_Bridge_Mean(phydbl T, phydbl A, phydbl B, phydbl u, phydbl *mean)
{
  phydbl pnormarg1 = log(B/A)/sqrt(u*u*T) + .5*sqrt(u*u*T);
  phydbl pnormarg2 = pnormarg1 - sqrt(u*u*T);

  /* printf("\n. %G %G",pnormarg1,pnormarg2); */
  
  if((pnormarg1 >  2.0 && pnormarg2 >  2.0) ||
     (pnormarg1 < -2.0 && pnormarg2 < -2.0))
    {
      // Proposition 5.1 in Privaut & Guindon, 2015, Journal of mathematical biology 71 (6-7), 1387-1409 
      *mean = T*(B-A)/log(B/A)+u*u*T*T*((A+B)/(2.*pow(log(B/A),2)) - (B-A)/pow(log(B/A),3));
    }
  else
    {
      // Proposition 3.3 in Privaut & Guindon, 2015, Journal of mathematical biology 71 (6-7), 1387-1409
      *mean =
        (A/(u*u))*sqrt(2.*PI*u*u*T)*exp(pow(u*u*T/2.+log(B/A),2)/(2.*u*u*T))*
        (Pnorm(pnormarg1,0.,1.) -
         Pnorm(pnormarg2,0.,1.));
    }
    
  *mean = *mean/T;
}
 
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Integrated_Geometric_Brownian_Bridge_Var(phydbl T, phydbl A, phydbl B, phydbl u, phydbl *var)
{

  double U = T*u*u;
  double logz = log(B/A);
  double z = B/A;
  double m;
  
  phydbl pnormarg1 = log(B/A)/sqrt(u*u*T) + .5*sqrt(u*u*T);
  phydbl pnormarg2 = pnormarg1 - sqrt(u*u*T);

  if((pnormarg1 >  2.0 && pnormarg2 >  2.0) ||
     (pnormarg1 < -2.0 && pnormarg2 < -2.0))
    {
      // Proposition 5.1 in Privaut & Guindon, 2015, Journal of mathematical biology 71 (6-7), 1387-1409 
      *var = 0.0;
    }
  else
    {
      // Proposition 3.3 in Privaut & Guindon, 2015, Journal of mathematical biology 71 (6-7), 1387-1409
      *var =
        2.*sqrt(2.*PI*U)*exp(pow(U   +logz,2)/(2.*U))*(Pnorm(logz/sqrt(U) +    sqrt(U),0.0,1.0)-Pnorm(logz/sqrt(U) -    sqrt(U),0.0,1.0)) -
 (1.+z)*2.*sqrt(2.*PI*U)*exp(pow(U/2.+logz,2)/(2.*U))*(Pnorm(logz/sqrt(U) + .5*sqrt(U),0.0,1.0)-Pnorm(logz/sqrt(U) - .5*sqrt(U),0.0,1.0));
      *var = *var * pow(A,2) / pow(u,4);
      Integrated_Geometric_Brownian_Bridge_Mean(T,A,B,u,&m);
      *var = *var / pow(T,2);
      *var = *var - m*m;
    }
  if(*var < 0.0) *var = 0.0;
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
  for(i=0;i<len;i++) u += pi[i];
  for(i=0;i<len;i++) cum_pi[i] = pi[i] / u;
  for(i=1;i<len;i++) cum_pi[i] += cum_pi[i-1];

  if((cum_pi[i-1] > 1. + 1.E-10) || (cum_pi[i-1] < 1. - 1.E-10))
    {
      PhyML_Printf("\n== Sum of probabilities is different from 1.0 (%f).",cum_pi[i-1]);
      PhyML_Printf("\n== Err. in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");
    }

  i = 0;
  u = Uni();
  for(i=0;i<len;i++) if(cum_pi[i] > u) break;

  if(i == len)
    {
      for(i=0;i<len;i++) printf("\n== idx:%d prob:%g",i,pi[i]);
      PhyML_Printf("\n== Len = %d",len);
      PhyML_Printf("\n== Err. in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");
    }
  Free(cum_pi);

  return(i);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// Same method as previous one, but returning n_elts sampled
// elements at once, avoiding to recreate cum_pi for each sample.
// Alias Method (https://en.wikipedia.org/wiki/Alias_method)
// Inspired from:
// https://jugit.fz-juelich.de/mlz/ransampl/blob/master/lib/ransampl.c
// and
// https://possiblywrong.wordpress.com/2012/02/05/the-alias-method-and-double-precision/
int* Sample_n_i_With_Proba_pi(phydbl *pi, int len, int n_elts)
{
  int i,n;
  phydbl sum;
  int *sampled;
  int *alias;
  phydbl *prob, *p;
  int *small, *large;
  int num_small = 0, num_large = 0;
  int a, g;
  
  alias = mCalloc(len, sizeof(int));
  prob = mCalloc(len, sizeof(phydbl));
  p = mCalloc(len, sizeof(phydbl));
  small = mCalloc(len, sizeof(int));
  large = mCalloc(len, sizeof(int));
  
  sampled = (int *)mCalloc(n_elts,sizeof(int));
  
  sum = .0;
  for(i=0;i<len;i++)
    {
      if( pi[i]<0 )
	{
	  PhyML_Printf("\n== Probability < 0");
	  PhyML_Printf("\n== Err. in file %s at line %d\n",__FILE__,__LINE__);
	  Exit("\n");
	}
      sum += pi[i];
    }

  if(sum == 0. )
    {
      PhyML_Printf("\n== Sum of probabilities is 0");
      PhyML_Printf("\n== Err. in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");
    }
  
  for(i=0;i<len;i++) p[i] = pi[i] * len / sum;

  num_small = 0;
  num_large = 0;
  for ( i=len-1; i>=0; --i )
    {
      if ( p[i]<1 )
	small[num_small++] = i;
      else
	large[num_large++] = i;
    }

  while ( num_small && num_large )
    {
      a = small[--num_small];
      g = large[--num_large];
      prob[a] = p[a];
      alias[a] = g;
      p[g] = p[g] + p[a] - 1;
      if ( p[g] < 1 )
	small[num_small++] = g;
      else
	large[num_large++] = g;
    }

  while ( num_large )
    prob[ large[--num_large] ] = 1;

  while ( num_small )
      prob[ small[--num_small] ] = 1;

  Free( p );
  Free( small );
  Free( large );

  for(n=0; n<n_elts;n++)
    {
      phydbl r1 = Uni();
      phydbl r2 = Uni();
      i = (int) (len * r1);
      sampled[n] = (r2 < prob[i] ? i : alias[i]);
    }

  Free(alias );
  Free(prob );
  
  return(sampled);
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
  for(i=0;i<len;i++) y[i] = x[i];

  do
    {
      swap = NO;
      for(i=0;i<len-1;i++) 
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
  for(i=0;i<len;i++) if(x[i] < z) hit+=1.;

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

  /* return(mu + sigma * SQRT(-log( y * y * (p_sup - p_inf) * (p_sup - p_inf) * 2 * PI * sigma * sigma))); */
  return(mu + sigma * SQRT(-log( y * y * (p_sup - p_inf) * (p_sup - p_inf) * 2. * PI * sigma * sigma)));
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// Returns a vector with a permutation of all the integers from 0
// to len-1. Fisher-Yates algorithm.

int *Permutate(int len)
{
  int i,pos,tmp;
  int *x;

  x = (int *)mCalloc(len,sizeof(int));

  for(i=0;i<len;i++) x[i] = i;

  for(i=0;i<len;i++)
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
  for(i=0;i<N;i++) sumx += x[i];

  sumy = .0;
  for(i=0;i<N;i++) sumy += y[i];

  sumxx = .0;
  for(i=0;i<N;i++) sumxx += x[i]*x[i];

  sumyy = .0;
  for(i=0;i<N;i++) sumyy += y[i]*y[i];

  sumxy = .0;
  for(i=0;i<N;i++) sumxy += x[i]*y[i];

  obs_stat = (N * sumxy - sumx * sumy) / (SQRT((N-1)*sumxx - (sumx/N)*(sumx/N)) * SQRT((N-1)*sumyy - (sumy/N)*(sumy/N)));
  
  npermut = 1000;
  p_val   = 0.0;
  for(k=0;k<npermut;k++)
    {
      permut = Permutate(nrow);

      sumxy = .0;
      for(i=0;i<nrow;i++)
        {
          for(j=0;j<ncol;j++)
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
  if(w) for(i=0;i<l;i++) wm += x[i]*w[i];
  else  
    {
      for(i=0;i<l;i++) wm += x[i];
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
  for(i=0;i<l;i++) sum += x[i]*x[i];
  
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
  for(i=0;i<range;i++)
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
  for(i=0;i<len;i++) s[i] = Rand_Int(97,121);
  s[i] = '\0';
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

  for(i=0;i<n;i++) permut[i] = i;
  
  for(i=0;i<n-1;i++)
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
  
  for(i=0;i<n;i++)
    {
      p->poly_vert[i]->lonlat[0] = Uni();
      p->poly_vert[i]->lonlat[1] = Uni();
    }
  
  return(p);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl Area_Of_Poly_Monte_Carlo(t_poly *poly, t_geo_coord *lim)
{
  int n_hit,n_trials,trial;
  t_geo_coord *point;

  point = (t_geo_coord *)GEO_Make_Geo_Coord(2);

  n_trials = 1E+7;
  trial = 0;
  n_hit = 0;
  do
    {
      point->lonlat[0] = Uni()*lim->lonlat[0];
      point->lonlat[1] = Uni()*lim->lonlat[1];
      if(Is_In_Polygon(point,poly) == YES) n_hit++;
      trial++;
    }
  while(trial < n_trials);


  Free_Geo_Coord(point);

  return((phydbl)(n_hit)/n_trials*lim->lonlat[0]*lim->lonlat[1]);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

int Is_In_Polygon(t_geo_coord *point, t_poly *poly)
{
  int i,j;
  phydbl x,y,x1,y1,x2,y2;
  phydbl x_intersect;
  short int is_in;

  assert(point);
  assert(poly);

  /* Coordinates of the point to test */
  x = point->lonlat[0];
  y = point->lonlat[1];

  j = poly->n_poly_vert-1;
  is_in = NO;
  for(i=0;i<poly->n_poly_vert;i++)
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
  
  if(x->dim != y->dim)
    {
      PhyML_Printf("\n. x->dim: %d y->dim: %d",x->dim,y->dim);
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);    
    }
  
  dist = 0.0;
  for(i=0;i<x->dim;i++) dist += pow(x->lonlat[i]-y->lonlat[i],2);

  return(sqrt(dist));
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl Manhattan_Dist(t_geo_coord *x, t_geo_coord *y)
{
  int i;
  phydbl dist;
  
  if(x->dim != y->dim)
    {
      PhyML_Printf("\n. x->dim: %d y->dim: %d",x->dim,y->dim);
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);    
    }
  
  dist = 0.0;
  for(i=0;i<x->dim;i++) dist += fabs(x->lonlat[i]-y->lonlat[i]);

  return(dist);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

// Haversine distance between (lon1,lat1) and (lon2,lat2) where
// coordinates are expressed using decimals
phydbl Haversine_Distance(t_geo_coord *x, t_geo_coord *y)
{

  phydbl R = 6371.; // Earth radius, in km
  phydbl a;
  phydbl lon1,lat1,lon2,lat2;
  
  lon1 = x->lonlat[0] * PI / 180.;
  lat1 = x->lonlat[1] * PI / 180.;

  lon2 = y->lonlat[0] * PI / 180.;
  lat2 = y->lonlat[1] * PI / 180.;

  a = pow(sin(.5*(lat2-lat1)),2);
  a += cos(lat1)*cos(lat2)*pow(sin(.5*(lon2-lon1)),2);

  return(2.* R * asin(sqrt(a)));

  /* return(sqrt(pow(lon1-lon2,2)+pow(lat1-lat2,2))); */
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

/* Return a vector rk such that rk[i] gives the index of the i-th largest element in x */
int *Ranks(phydbl *x, int len)
{
  int *rk,tmp;
  int i,swap;

  rk = (int *)mCalloc(len,sizeof(int));

  for(i=0;i<len;i++) rk[i] = i;

  do
    {
      swap = NO;
      for(i=0;i<len-1;i++) 
        {
          if(x[rk[i]] < x[rk[i+1]])
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
  for(i=0;i<n_steps-1;i++) if(!(time[i+1] > time[i])) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);

  state = Brownian_Generate(var,n_steps,beg_time,time);
  end_brown = Rnorm(state[n_steps-1],SQRT((time[n_steps-1] - end_time)*var));

  for(i=0;i<n_steps;i++)
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

  for(i=0;i<n_steps;i++)
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
