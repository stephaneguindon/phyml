/*
PhyML :  a program that  computes maximum likelihood  phylogenies from
DNA or AA homologous sequences
Copyright (C) Stephane Guindon. Oct 2003 onward
All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.
*/

#include "optimiz.h"

static phydbl Br_Len_Spline(phydbl *l, t_edge *b, int n_iter_max, phydbl tol, t_tree *tree);

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Optimize_Single_Param_Generic(t_tree *tree, phydbl *param, phydbl lim_inf, phydbl lim_sup, phydbl tol, int n_max_iter, int quickdirty)
{
  phydbl lk_init;
  
  lk_init = tree->c_lnL;

  Generic_Brent_Lk(param,
                   lim_inf,
                   lim_sup,
                   tol,
                   n_max_iter,
                   quickdirty,
                   Wrap_Lk,
                   NULL,
                   tree,
                   NULL,
                   NO,NO);


  if(tree->c_lnL < lk_init - tree->mod->s_opt->min_diff_lk_global)
    {
      PhyML_Fprintf(stderr,"\n. %.10f < %.10f --> diff=%.10f param value = %f\n",tree->c_lnL,lk_init,tree->c_lnL-lk_init,*param);
      assert(FALSE);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int Generic_Brak(phydbl *param,
                 phydbl *ax, phydbl *bx, phydbl *cx,
                 phydbl *fa, phydbl *fb, phydbl *fc,
                 phydbl lim_inf, phydbl lim_sup,
                 t_tree *tree)
{
   phydbl ulim,u,r,q,fu,dum;

   u = 0.0;
   *param = *ax;

   if(*param > lim_sup) *param = lim_sup;
   if(*param < lim_inf) *param = lim_inf;
   *fa=-Lk(NULL,tree);
   *param = *bx;
   if(*param > lim_sup) *param = lim_sup;
   if(*param < lim_inf) *param = lim_inf;
   *fb=-Lk(NULL,tree);
   if (*fb > *fa) {
     SHFT(dum,*ax,*bx,dum)
       SHFT(dum,*fb,*fa,dum)
       }
   *cx=(*bx)+MNBRAK_GOLD*(*bx-*ax);
   *param = FABS(*cx);
   if(*param > lim_sup) *param = lim_sup;
   if(*param < lim_inf) *param = lim_inf;
   *fc=-Lk(NULL,tree);
   while (*fb > *fc)
     {
       
       if(*ax > lim_sup) *ax = lim_sup;
       if(*ax < lim_inf) *ax = lim_inf;
       if(*bx > lim_sup) *bx = lim_sup;
       if(*bx < lim_inf) *bx = lim_inf;
       if(*cx > lim_sup) *cx = lim_sup;
       if(*cx < lim_inf) *cx = lim_inf;
       if(u   > lim_sup) u   = lim_sup;
       if(u   < lim_inf) u   = lim_inf;
       
       r=(*bx-*ax)*(*fb-*fc);
       q=(*bx-*cx)*(*fb-*fa);
       u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
         (2.0*SIGN(MAX(FABS(q-r),MNBRAK_TINY),q-r));
       ulim=(*bx)+MNBRAK_GLIMIT*(*cx-*bx);
       
       if ((*bx-u)*(u-*cx) > lim_inf)
         {
           *param = FABS(u);
           if(*param > lim_sup) {*param = u = lim_sup;}
           if(*param < lim_inf) {*param = u = lim_inf;}
           fu=-Lk(NULL,tree);
           if (fu < *fc)
             {
               *ax=(*bx);
               *bx=u;
               *fa=(*fb);
               *fb=fu;
               (*ax)=FABS(*ax);
               (*bx)=FABS(*bx);
               (*cx)=FABS(*cx);
               return(0);
             }
           else if (fu > *fb)
             {
               *cx=u;
               *fc=fu;
               (*ax)=FABS(*ax);
               (*bx)=FABS(*bx);
               (*cx)=FABS(*cx);
               return(0);
             }
           u=(*cx)+MNBRAK_GOLD*(*cx-*bx);
           *param = FABS(u);
           if(*param > lim_sup) {*param = u = lim_sup;}
           if(*param < lim_inf) {*param = u = lim_inf;}
           fu=-Lk(NULL,tree);
         }
       else if ((*cx-u)*(u-ulim) > lim_inf)
         {
           *param = FABS(u);
           if(*param > lim_sup) {*param = u = lim_sup;}
           if(*param < lim_inf) {*param = u = lim_inf;}
           fu=-Lk(NULL,tree);
           if (fu < *fc)
             {
               SHFT(*bx,*cx,u,*cx+MNBRAK_GOLD*(*cx-*bx))
                 *param = FABS(u);
               SHFT(*fb,*fc,fu,-Lk(NULL,tree))
                 }
         }
       else if ((u-ulim)*(ulim-*cx) >= lim_inf)
         {
           u=ulim;
           *param = FABS(u);
           if(*param > lim_sup) {*param = u = lim_sup;}
           if(*param < lim_inf) {*param = u = lim_inf;}
           fu=-Lk(NULL,tree);
         }
       else
         {
           u=(*cx)+MNBRAK_GOLD*(*cx-*bx);
           *param = FABS(u);
           if(*param > lim_sup) {*param = u = lim_sup;}
           if(*param < lim_inf) {*param = u = lim_inf;}
           fu=-Lk(NULL,tree);
         }
       SHFT(*ax,*bx,*cx,u)
         SHFT(*fa,*fb,*fc,fu)
         
         
         }
   (*ax)=FABS(*ax);
   (*bx)=FABS(*bx);
   (*cx)=FABS(*cx);
   return(0);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int Generic_Brak_Lk(phydbl *param,
                    phydbl *ax, phydbl *bx, phydbl *cx,
                    phydbl *fa, phydbl *fb, phydbl *fc,
                    phydbl min, phydbl max,
                    phydbl (*obj_func)(t_edge *,t_tree *,supert_tree *),
                    t_edge *branch, t_tree *tree, supert_tree *stree)
{
   phydbl ulim,u,r,q,fu,dum;

   *param = *ax;
   if(*param < min) *param = min;
   if(*param > max) *param = max;
   *fa = -(*obj_func)(branch,tree,stree);

   *param = *bx;
   if(*param < min) *param = min;
   if(*param > max) *param = max;
   *fb = -(*obj_func)(branch,tree,stree);

   if (*fb > *fa) 
     {
       SHFT(dum,*ax,*bx,dum)
       SHFT(dum,*fb,*fa,dum)
     }

   *cx=(*bx)+MNBRAK_GOLD*(*bx-*ax);
   *param = *cx;
   if(*param < min) *param = min;
   if(*param > max) *param = max;
   *fc = -(*obj_func)(branch,tree,stree);

   /* printf("\nx la: %f lb: %f lc: %f fa: %f fb: %f fc: %f",*ax,*bx,*cx,*fa,*fb,*fc); */
   while (*fb > *fc)
     {       
       /* printf("\nx la: %f lb: %f lc: %f fa: %f fb: %f fc: %f",*ax,*bx,*cx,*fa,*fb,*fc); */
       r=(*bx-*ax)*(*fb-*fc);
       q=(*bx-*cx)*(*fb-*fa);
       u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
         (2.0*SIGN(MAX(FABS(q-r),MNBRAK_TINY),q-r));
       ulim=(*bx)+MNBRAK_GLIMIT*(*cx-*bx);
       
       if((*bx-u)*(u-*cx) > 0.0)
         {
           *param = u;
           if(*param < min) *param = min;
           if(*param > max) *param = max;
           fu = -(*obj_func)(branch,tree,stree);

           if (fu < *fc)
             {
               *ax=(*bx);
               *bx=u;
               *fa=(*fb);
               *fb=fu;
               return(0);
             }
           else if (fu > *fb)
             {
               *cx=u;
               *fc=fu;
               return(0);
             }
           u=(*cx)+MNBRAK_GOLD*(*cx-*bx);
           *param = u;
           if(*param < min) *param = min;
           if(*param > max) *param = max;
           fu = -(*obj_func)(branch,tree,stree);
         }
       else if ((*cx-u)*(u-ulim) > 0.0)
         {
           *param = u;
           if(*param < min) *param = min;
           if(*param > max) *param = max;
           fu = -(*obj_func)(branch,tree,stree);

           if (fu < *fc)
             {
               SHFT(*bx,*cx,u,*cx+MNBRAK_GOLD*(*cx-*bx))
               *param = u;
               if(*param < min) *param = min;
               if(*param > max) *param = max;
               SHFT(*fb,*fc,fu,-(*obj_func)(branch,tree,stree))
             }
         }
       else if ((u-ulim)*(ulim-*cx) >= 0.0)
         {
           u=ulim;
           *param = u;
           if(*param < min) *param = min;
           if(*param > max) *param = max;
           fu = -(*obj_func)(branch,tree,stree);
         }
       else
         {
           u=(*cx)+MNBRAK_GOLD*(*cx-*bx);
           *param = u;
           if(*param < min) *param = min;
           if(*param > max) *param = max;
           fu = -(*obj_func)(branch,tree,stree);
         }
       SHFT(*ax,*bx,*cx,u)
       SHFT(*fa,*fb,*fc,fu)
         }
   
   return(0);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl Generic_Brent(phydbl *param, phydbl ax, phydbl cx, phydbl tol,
                     int n_iter_max,
                     phydbl (*obj_func)(t_tree *),
                     t_tree *tree)
{
  int iter;
  phydbl a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  phydbl e=0.0;
  phydbl old_score;
  phydbl bx = *param;
  
  d=0.0;
  a=((ax < cx) ? ax : cx);
  b=((ax > cx) ? ax : cx);
  x=w=v=bx;
  (*param) = bx;
  fw=fv=fx=fu=(*obj_func)(tree);
  old_score = fw;
  
  /* PhyML_Printf("\n. %p init_score=%f fu=%f ax=%f cx=%f param=%f",tree,init_score,fu,ax,cx,*param); */
  
  for(iter=1;iter<=BRENT_IT_MAX;iter++)
    {
      xm=0.5*(a+b);
      tol2=2.0*(tol1=tol*x+BRENT_ZEPS);
      
      if((fabs(fu-old_score) < tol && iter > 1) || (iter > n_iter_max - 1))
        {
          (*param) = x;
          fu = (*obj_func)(tree);
          /* printf("\n. return %f [%f] %d",*param,fu,iter); */
          return fu;
        }

      if(fabs(e) > tol1)
        {
          r=(x-w)*(fx-fv);
          q=(x-v)*(fx-fw);
          p=(x-v)*q-(x-w)*r;
          q=2.0*(q-r);
          if(q > 0.0) p = -p;
          q=fabs(q);
          etemp=e;
          e=d;
          if(fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
            {
              d=BRENT_CGOLD*(e=(x >= xm ? a-x : b-x));
              /* PhyML_Printf("\n. Golden section step"); */
            }
          else
            {
              d=p/q;
              u=x+d;
              if (u-a < tol2 || b-u < tol2) d=SIGN(tol1,xm-x);
              /* PhyML_Printf("\n. Parabolic step [e=%f]",e); */
            }
        }
      else
        {
          d=BRENT_CGOLD*(e=(x >= xm ? a-x : b-x));
          /* PhyML_Printf("\n. Golden section step (default) [e=%f tol1=%f a=%f b=%f d=%f x=%f]",e,tol1,a,b,d,x); */
        }

      u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
      (*param) = u;
      old_score = fu;
      fu = (*obj_func)(tree);
      /* PhyML_Printf("\n. iter=%d/%d param=%f lnL=%f u: %f x: %f d: %f logt: %d",iter,BRENT_IT_MAX,*param,fu,u,x,d,logt); */

      if(fu <= fx)
        {
          if(u >= x) a=x; else b=x;
          SHFT(v,w,x,u)
          SHFT(fv,fw,fx,fu)
        }
      else
        {
          if (u < x) a=u; else b=u;
          if (fu < fw || fabs(w-x) < SMALL)
            {
              v=w;
              w=u;
              fv=fw;
              fw=fu;
            }
          /* 	  else if (fu <= fv || v == x || v == w) */
          else if (fu < fv || fabs(v-x) < SMALL || fabs(v-w) < SMALL)
            {
              v=u;
              fv=fu;
            }
        }
    }

  PhyML_Printf("\n. Too many iterations in Generic_Brent !");
  assert(FALSE);
  return((*obj_func)(tree));
  /* Not Reached ??  *param=x;   */
  /* Not Reached ??  return fx; */
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl RRparam_GTR_Golden(phydbl ax, phydbl bx, phydbl cx, phydbl tol,
              phydbl *xmin, t_tree *tree, calign *cdata, phydbl *param, int n_iter_max)
{
   phydbl f1,f2,x0,x1,x2,x3;
   int n_iter;


   x0=ax;
   x3=cx;
   if (FABS(cx-bx) > FABS(bx-ax))
     {
       x1=bx;
       x2=bx+GOLDEN_C*(cx-bx);
     }
   else
     {
       x2=bx;
       x1=bx-GOLDEN_C*(bx-ax);
     }
   (*param)=x1;

   Lk(NULL,tree);
   f1=-tree->c_lnL;
   (*param)=x2;

   Lk(NULL,tree);
   f2=-tree->c_lnL;

   n_iter = 0;
   while (FABS(x3-x0) > tol*(FABS(x1)+FABS(x2)))
     {

       if (f2 < f1)
     {
       SHFT3(x0,x1,x2,GOLDEN_R*x1+GOLDEN_C*x3)
       (*param)=x2;
       Lk(NULL,tree);
       SHFT2(f1,f2,-tree->c_lnL)
     }
       else
     {
       SHFT3(x3,x2,x1,GOLDEN_R*x2+GOLDEN_C*x0)
       (*param)=x1;
       Lk(NULL,tree);
       SHFT2(f2,f1,-tree->c_lnL)
     }

       if(n_iter++ > n_iter_max) break;

/*        PhyML_Printf("p=%E %f\n",(*param),tree->c_lnL); */
     }
   if (f1 < f2)
    {
       *xmin=x1;
       return f1;
     }
   else
     {
       *xmin=x2;
       return f2;
     }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Br_Len_Golden(phydbl ax, phydbl bx, phydbl cx, phydbl tol,
                     phydbl *xmin, t_edge *b_fcus, t_tree *tree)
{
   phydbl f1,f2,x0,x1,x2,x3;

   x0=ax;
   x3=cx;
   if(FABS(cx-bx) > FABS(bx-ax))
     {
       x1=bx;
       x2=bx+GOLDEN_C*(cx-bx);
     }
   else
     {
       x2=bx;
       x1=bx-GOLDEN_C*(bx-ax);
     }
   
   b_fcus->l->v = x1;
   f1 = -Lk(b_fcus,tree);
   b_fcus->l->v = x2;
   f2 = -Lk(b_fcus,tree);

   while (FABS(x3-x0) > tol*(FABS(x1)+FABS(x2)))
     {
       if (f2 < f1)
         {
           SHFT3(x0,x1,x2,GOLDEN_R*x1+GOLDEN_C*x3)
           b_fcus->l->v = x2;
           SHFT2(f1,f2,-Lk(b_fcus,tree))
         }
       else
         {
           SHFT3(x3,x2,x1,GOLDEN_R*x2+GOLDEN_C*x0)
           b_fcus->l->v = x1;
           SHFT2(f2,f1,-Lk(b_fcus,tree))
         }
     }

   if (f1 < f2)
     {
       *xmin = x1;
       return -f1;
     }
   else
     {
       *xmin = x2;
       return -f2;
     }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


int Br_Len_Brak(phydbl *ax, phydbl *bx, phydbl *cx,
        phydbl *fa, phydbl *fb, phydbl *fc,
        t_edge *b_fcus, t_tree *tree)
{
   phydbl ulim,u,r,q,fu,dum;

   b_fcus->l->v = *ax;
   *fa=-Lk(b_fcus,tree);
   b_fcus->l->v = *bx;
   *fb=-Lk(b_fcus,tree);
   if (*fb > *fa) {
      SHFT(dum,*ax,*bx,dum)
      SHFT(dum,*fb,*fa,dum)
   }
   *cx=(*bx)+MNBRAK_GOLD*(*bx-*ax);
   b_fcus->l->v = *cx;
   *fc=-Lk(b_fcus,tree);
   while (*fb > *fc + tree->mod->s_opt->min_diff_lk_local)
     {
       r=(*bx-*ax)*(*fb-*fc);
       q=(*bx-*cx)*(*fb-*fa);
       u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
               (2.0*SIGN(MAX(FABS(q-r),MNBRAK_TINY),q-r));
       ulim=(*bx)+MNBRAK_GLIMIT*(*cx-*bx);

       if ((*bx-u)*(u-*cx) > 0.0)
     {
       b_fcus->l->v = u;
       fu=-Lk(b_fcus,tree);
       if (fu < *fc)
         {
           *ax=(*bx);
           *bx=u;
           *fa=(*fb);
           *fb=fu;
/* 	       (*ax)=FABS(*ax); */
/* 	       (*bx)=FABS(*bx); */
/* 	       (*cx)=FABS(*cx); */
           return(0);
         }
       else if (fu > *fb)
         {
           *cx=u;
           *fc=fu;
/* 	       (*ax)=FABS(*ax); */
/* 	       (*bx)=FABS(*bx); */
/* 	       (*cx)=FABS(*cx); */
           return(0);
         }
       u=(*cx)+MNBRAK_GOLD*(*cx-*bx);
       b_fcus->l->v = u;
       fu=-Lk(b_fcus,tree);
     }
       else if ((*cx-u)*(u-ulim) > 0.0)
     {
       b_fcus->l->v = FABS(u);
       fu=-Lk(b_fcus,tree);
       if (fu < *fc)
         {
           SHFT(*bx,*cx,u,*cx+MNBRAK_GOLD*(*cx-*bx))
           b_fcus->l->v = u;
           SHFT(*fb,*fc,fu,-Lk(b_fcus,tree))
         }
     }
       else if ((u-ulim)*(ulim-*cx) >= 0.0)
     {
       u=ulim;
       b_fcus->l->v = u;
       fu=-Lk(b_fcus,tree);
     }
       else
     {
       u=(*cx)+MNBRAK_GOLD*(*cx-*bx);
       b_fcus->l->v = u;
       fu=-Lk(b_fcus,tree);
     }
       SHFT(*ax,*bx,*cx,u)
       SHFT(*fa,*fb,*fc,fu)
      }
   (*ax)=FABS(*ax);
   (*bx)=FABS(*bx);
   (*cx)=FABS(*cx);
   return(0);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl Fast_Br_Len(t_edge *b, t_tree *tree, int approx)
{
  phydbl init_min_diff_lk_local = tree->mod->s_opt->min_diff_lk_local;

  tree->mod->s_opt->min_diff_lk_local = 1.E-1;      

  if(tree->is_mixt_tree) MIXT_Br_Len_Opt(b,tree);   
  else Br_Len_Opt(&(b->l->v),b,tree);

  tree->mod->s_opt->min_diff_lk_local = init_min_diff_lk_local;

  return tree->c_lnL;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl Br_Len_Opt(phydbl *l, t_edge *b, t_tree *tree)
{
  phydbl lk_begin, lk_end;

  if(tree->is_mixt_tree == YES && tree->ignore_mixt_info == NO)
    {
      MIXT_Br_Len_Opt(b,tree);
      return tree->c_lnL;
    }

  if(b->l->onoff == OFF || b->l->optimize == NO) return tree->c_lnL;

  lk_begin = UNLIKELY;
  lk_end   = UNLIKELY;
  
  Set_Update_Eigen_Lr(YES,tree);
  Set_Use_Eigen_Lr(NO,tree);

  lk_begin = Lk(b,tree); /* We can't assume that the log-lk value is up-to-date */

  Set_Update_Eigen_Lr(NO,tree);
  Set_Use_Eigen_Lr(YES,tree);

  Br_Len_Spline(l,b,tree->mod->s_opt->brent_it_max,tree->mod->s_opt->min_diff_lk_local,tree);
  
  Update_PMat_At_Given_Edge(b,tree);
  
  Set_Update_Eigen_Lr(NO,tree);
  Set_Use_Eigen_Lr(NO,tree);

  /* Set_Update_Eigen_Lr(NO,tree); */
  /* Set_Use_Eigen_Lr(NO,tree); */
  /* lk_begin = Lk(b,tree); */
  /* tree->n_tot_bl_opt += Generic_Brent_Lk(l, */
  /*                                        tree->mod->l_min, */
  /*                                        tree->mod->l_max, */
  /*                                        tree->mod->s_opt->min_diff_lk_local, */
  /*                                        tree->mod->s_opt->brent_it_max, */
  /*                                        tree->mod->s_opt->quickdirty, */
  /*                                        Wrap_Lk_At_Given_Edge, */
  /*                                        b,tree,NULL,NO,NO); */

  /* lk_end = Lk(b,tree); /\* We can't assume that the log-lk value is up-to-date *\/ */

  /* PhyML_Printf("\n. b->num: %4d l=%12G lnL: %12G",b->num,b->l->v,tree->c_lnL); */
  

  lk_end = tree->c_lnL;

  if(lk_end < lk_begin - tree->mod->s_opt->min_diff_lk_local)
    {
      PhyML_Fprintf(stderr,"\n. lk_beg = %f lk_end = %f",lk_begin, lk_end);
      PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d",__FILE__,__LINE__);
      Exit("\n");
    }

  return tree->c_lnL;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Round_Optimize(t_tree *tree, int n_round_max)
{
  int n_round,each,freq;
  phydbl lk_old, lk_new;

  lk_new  = tree->c_lnL;
  lk_old  = UNLIKELY;
  n_round = 0;
  each    = 0;
  freq    = 1;
  
  while(n_round < n_round_max)
    {
      if(tree->mod->s_opt->opt_bl_one_by_one || tree->mod->s_opt->constrained_br_len) Optimize_Br_Len_Serie(n_round_max,tree);
      
      if((tree->mod->s_opt->opt_bl_one_by_one || tree->mod->s_opt->constrained_br_len) &&
         (tree->verbose > VL2) &&
         (tree->io->quiet == NO)) Print_Lk(tree,"[Branch lengths     ]");     

      if(!each)
        {
          each = freq;
          Optimiz_All_Free_Param(tree,(tree->io->quiet)?(0):(tree->verbose > VL2));
        }

      lk_new = tree->c_lnL;
      
      if(lk_new < lk_old - tree->mod->s_opt->min_diff_lk_local)
        {
          PhyML_Fprintf(stderr,"\n. lk_new = %f lk_old = %f diff = %f",lk_new,lk_old,lk_new-lk_old);
          assert(FALSE);
        }
      
      if((fabs(lk_new - lk_old) < tree->mod->s_opt->min_diff_lk_local) && (each == freq)) break;
      /* if(fabs(lk_new - lk_old) < tree->mod->s_opt->min_diff_lk_local) break; */
      lk_old  = lk_new;

      n_round++;
      each--;
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Optimize_Br_Len_Serie(int n_max_iter, t_tree *tree)
{
  phydbl lk_init,lk_end;
  int iter;
  
  Set_Both_Sides(NO,tree);
  Lk(NULL,tree);
  
  lk_init = tree->c_lnL;

  Optimize_Lvar(tree,(tree->io->quiet)?(0):(tree->verbose > VL2));
  
  iter = 0;
  do
    {      
      lk_init = tree->c_lnL;
      if(tree->n_root && tree->ignore_root == NO)
        {
          Update_Partial_Lk(tree,tree->n_root->b[1],tree->n_root);
          Optimize_Br_Len_Serie_Post(tree->n_root,tree->n_root->v[1],tree->n_root->b[1],tree);
          Update_Partial_Lk(tree,tree->n_root->b[2],tree->n_root);
          Optimize_Br_Len_Serie_Post(tree->n_root,tree->n_root->v[2],tree->n_root->b[2],tree);
        }
      else if(tree->n_root && tree->ignore_root == YES)
        {
          Optimize_Br_Len_Serie_Post(tree->e_root->rght,
                                     tree->e_root->left,
                                     tree->e_root,tree);
          
          Optimize_Br_Len_Serie_Post(tree->e_root->left,
                                     tree->e_root->rght,
                                     tree->e_root,tree);
        }
      else
        {
          Optimize_Br_Len_Serie_Post(tree->a_nodes[tree->tip_root],
                                     tree->a_nodes[tree->tip_root]->v[0],
                                     tree->a_nodes[tree->tip_root]->b[0],
                                     tree);
        }
      
      lk_end = tree->c_lnL;
      
      if(lk_end < lk_init - tree->mod->s_opt->min_diff_lk_local)
        {
          PhyML_Fprintf(stderr,"\n. lk_init: %f lk_end: %f",lk_init,lk_end);
          Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
        }

      iter++;
    }
  while(iter < 1);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void Optimize_Br_Len_Multiplier(t_tree *mixt_tree, int verbose)
{
  phydbl lk_init;
  t_tree *tree;

  /* PhyML_Printf("\n. Optimizing Br_Len_Multiplier"); */
  
  tree = mixt_tree;
  
  do
    {      
      if(tree->mod->s_opt->opt_br_len_mult == YES)
        {
          lk_init = Get_Lk(tree);
          /* Generic_Brent_Lk(&(tree->mod->br_len_mult_unscaled->v), */
          Generic_Brent_Lk(&(tree->mod->br_len_mult->v),
                           1.E-2,1.E+1,
                           tree->mod->s_opt->min_diff_lk_local,
                           tree->mod->s_opt->brent_it_max,
                           tree->mod->s_opt->quickdirty,
                           Wrap_Lk,NULL,mixt_tree,NULL,NO,NO);

          if(Get_Lk(tree) < lk_init - tree->mod->s_opt->min_diff_lk_local)
            {
              PhyML_Fprintf(stderr,"\n. %f -- %f",lk_init,tree->c_lnL);
              Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
            }
        }
      tree = tree->next_mixt;
    }
  while(tree);

  tree = mixt_tree;  
  do
    {      
      if(verbose && tree->mod->s_opt->opt_br_len_mult == YES)
        {
          Print_Lk(tree,"[Tree scale         ]");
          PhyML_Printf("[%10f]",tree->mod->br_len_mult->v);
        }

      tree = tree->next_mixt;
    }
  while(tree);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void Optimize_Br_Len_Serie_Post(t_node *a, t_node *d, t_edge *b_fcus, t_tree *tree)
{
  int i;
  phydbl lk_init;
                  
  lk_init = tree->c_lnL;
  
  if(tree->mod->s_opt->constrained_br_len == YES)
    {
      Generic_Brent_Lk(&(tree->mod->br_len_mult->v),
                       1.E-2,1.E+1,
                       tree->mod->s_opt->min_diff_lk_local,
                       tree->mod->s_opt->brent_it_max,
                       tree->mod->s_opt->quickdirty,
                       Wrap_Lk,NULL,tree,NULL,NO,NO);

      if(tree->c_lnL < lk_init - tree->mod->s_opt->min_diff_lk_local)
        {
          PhyML_Fprintf(stderr,"\n. %f -- %f",lk_init,tree->c_lnL);
          Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
        }

      return;
    }

  if(tree->io->mod->s_opt->opt_bl_one_by_one == YES) Br_Len_Opt(&(b_fcus->l->v),b_fcus,tree);

  if(tree->c_lnL < lk_init - tree->mod->s_opt->min_diff_lk_local)
    {
      PhyML_Fprintf(stderr,"\n. %f -- %f",lk_init,tree->c_lnL);
      PhyML_Fprintf(stderr,"\n. Edge: %d",b_fcus->num);
      PhyML_Fprintf(stderr,"\n. is_mixt_tree: %d",tree->is_mixt_tree);
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
    }

  if(d->tax) return;

  if(tree->n_root && tree->ignore_root == NO)
    {
      for(i=0;i<3;++i)
        {
          if(d->v[i] != a && d->b[i] != tree->e_root)
            {
              Update_Partial_Lk(tree,d->b[i],d);
              Optimize_Br_Len_Serie_Post(d,d->v[i],d->b[i],tree);
            }
        }
      
      for(i=0;i<3;++i)
        if(d->v[i] == a || d->b[i] == tree->e_root) 
          {
            Update_Partial_Lk(tree,d->b[i],d);
            if(tree->io->mod->s_opt->opt_bl_one_by_one == YES) Br_Len_Opt(&(d->b[i]->l->v),d->b[i],tree);
          }
    }
  else
    {
      // Ok if root exists but requires traversal to be initiated from a node != root
      for(i=0;i<3;++i)
        {
          if(d->v[i] != a)
            {
              Update_Partial_Lk(tree,d->b[i],d);
              Optimize_Br_Len_Serie_Post(d,d->v[i],d->b[i],tree);
            }
        }
      Update_Partial_Lk(tree,b_fcus,d);
      if(tree->io->mod->s_opt->opt_bl_one_by_one == YES) Br_Len_Opt(&(b_fcus->l->v),b_fcus,tree);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Optimiz_Ext_Br(t_tree *tree)
{
  int i;
  t_edge *b;
  phydbl lk_init;
  scalar_dbl *l_init,*v_init;

  lk_init = tree->c_lnL;

  For(i,2*tree->n_otu-3)
    {
      b = tree->a_edges[i];
      if((b->left->tax) || (b->rght->tax))
        {
          
          l_init = Duplicate_Scalar_Dbl(b->l);          
          v_init = Duplicate_Scalar_Dbl(b->l_var);          
          
          Br_Len_Opt(&(b->l->v),b,tree);
          
          if(b->nni->best_l == NULL)
            {
              b->nni->best_l = Duplicate_Scalar_Dbl(b->l);
              b->nni->best_v = Duplicate_Scalar_Dbl(b->l_var);
            }
          else
            {
              Copy_Scalar_Dbl(b->l,b->nni->best_l);
              Copy_Scalar_Dbl(b->l_var,b->nni->best_v);
            }

          if(b->nni->l0 == NULL)
            {
              b->nni->l0 = Duplicate_Scalar_Dbl(b->l);
              b->nni->v0 = Duplicate_Scalar_Dbl(b->l_var);
            }
          else
            {
              Copy_Scalar_Dbl(b->l,b->nni->l0);
              Copy_Scalar_Dbl(b->l_var,b->nni->v0);
            }

          // Revert of original edge lengths
          Copy_Scalar_Dbl(l_init,b->l);
          Copy_Scalar_Dbl(v_init,b->l_var);

          Free_Scalar_Dbl(l_init);
          Free_Scalar_Dbl(v_init);

          /* ori = b; */
          /* do */
          /*   { */
          /*     b->nni->best_l->v = b->l->v; */
          /*     b->nni->l0->v     = b->l->v; */
          /*     b->nni->best_conf = 0; */
          /*     b->l->v           = l_init; */
          /*     b = b->next; */
          /*   } */
          /* while(b); */
          /* b = ori; */
        }
    }
  tree->c_lnL = lk_init;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Optimiz_All_Free_Param(t_tree *tree, int verbose)
{
  int  init_both_sides;
  
  if(!tree) return;

  if(tree->mixt_tree && tree->mod->ras->invar == YES) return;

  init_both_sides  = tree->both_sides;


  Set_Both_Sides(NO,tree);
  Lk(NULL,tree);

  Optimize_RR_Params(tree,verbose);
  Optimize_TsTv(tree,verbose);
  Optimize_Lambda(tree,verbose);
  Optimiz_Alpha_And_Pinv(tree,verbose);
  Optimize_Pinv(tree,verbose);
  Optimize_Alpha(tree,verbose);
  Optimize_State_Freqs(tree,verbose);
  Optimize_Rmat_Weights(tree,verbose);
  Optimize_Efrq_Weights(tree,verbose);
  Optimize_Free_Rate(tree,verbose);
  Optimize_Br_Len_Multiplier(tree,verbose);
  Optimize_M4mod(tree, verbose);

  if(tree->io->print_json_trace == YES) JSON_Tree_Io(tree,tree->io->fp_out_json_trace); 
  
  Set_Both_Sides(init_both_sides,tree);

  if(tree->both_sides == YES) Lk(NULL,tree); /* Needed to update all partial likelihoods */


  /* if(tree->next) Optimiz_All_Free_Param(tree->next,verbose); */
  /* else            Optimiz_All_Free_Param(tree->next,verbose); */

  /* if(tree->nextree) Optimiz_All_Free_Param(tree->nextree,verbose); */
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Optimize_M4mod(t_tree *tree, int verbose)
{
  if (tree->mod->use_m4mod == YES)
  {
    int i;

    if (tree->mod->m4mod->delta->optimize == YES)
    {
      Set_Update_Eigen(YES,tree->mod);

      Generic_Brent_Lk(
          &(tree->mod->m4mod->delta->v), 0.01, 10.,
          tree->mod->s_opt->min_diff_lk_local, tree->mod->s_opt->brent_it_max,
          tree->mod->s_opt->quickdirty, Wrap_Lk, NULL, tree, NULL, NO, NO);

      if (verbose)
      {
        Print_Lk(tree, "[Switching param.   ]");
        PhyML_Printf("[%10f]", tree->mod->m4mod->delta->v);
      }

      Set_Update_Eigen(NO, tree->mod);
    }

    if (tree->mod->m4mod->multipl_unscaled->optimize == YES)
    {
    
      Set_Update_Eigen(YES,tree->mod);

      for (int rcat = 0; rcat < tree->mod->m4mod->n_h; rcat++)
      {
        /* 	      Optimize_Single_Param_Generic(tree,&(tree->mod->m4mod->multipl_unscaled[rcat]),
         */
        /* 					    .01,10., */
        /* 					    tree->mod->s_opt->min_diff_lk_local,
         */
        /* 					    tree->mod->s_opt->brent_it_max,
         */
        /* 					    tree->mod->s_opt->quickdirty);
         */

        Generic_Brent_Lk(
            &(tree->mod->m4mod->multipl_unscaled->v[rcat]), 0.1, 100.,
            tree->mod->s_opt->min_diff_lk_local, tree->mod->s_opt->brent_it_max,
            tree->mod->s_opt->quickdirty, Wrap_Lk, NULL, tree, NULL, NO, NO);

        if (verbose)
        {
          Print_Lk(tree, "[Rel. subst. rate   ]");
          PhyML_Printf("[%10f]", tree->mod->m4mod->multipl[rcat]);
        }
      }

      Set_Update_Eigen(NO, tree->mod);
    }

    if (tree->mod->m4mod->multipl_unscaled->optimize == YES)
    {
      Set_Update_Eigen(YES, tree->mod);

      for (int rcat = 0; rcat < tree->mod->m4mod->n_h; rcat++)
      {
        Generic_Brent_Lk(
            &(tree->mod->m4mod->h_fq_unscaled->v[rcat]), 0.1, 100.,
            tree->mod->s_opt->min_diff_lk_local, tree->mod->s_opt->brent_it_max,
            tree->mod->s_opt->quickdirty, Wrap_Lk, NULL, tree, NULL, NO, NO);

        if (verbose)
        {
          Print_Lk(tree, "[Subst. class freq  ]");
          PhyML_Printf("[%10f]", tree->mod->m4mod->h_fq[rcat]);
        }
      }

      Set_Update_Eigen(NO, tree->mod);
    }

    if (tree->mod->m4mod->alpha->optimize == YES)
    {
      Set_Update_Eigen(YES, tree->mod);

      Generic_Brent_Lk(
          &(tree->mod->m4mod->alpha->v), 0.01, 10.,
          tree->mod->s_opt->min_diff_lk_local, tree->mod->s_opt->brent_it_max,
          tree->mod->s_opt->quickdirty, Wrap_Lk, NULL, tree, NULL, NO, NO);

      if (verbose)
      {
        Print_Lk(tree, "[Alpha (covarion)   ]");
        PhyML_Printf("[%10f]", tree->mod->m4mod->alpha->v);
      }

      Set_Update_Eigen(NO, tree->mod);
    }

    /* Substitutions between nucleotides are considered to follow a
       GTR model */
    if (tree->mod->io->datatype == NT && (tree->mod->whichmodel == GTR || tree->mod->whichmodel == CUSTOM) && tree->mod->m4mod->o_rr->optimize == YES)
      {
        Set_Update_Eigen(YES, tree->mod);

        int *permut = Permutate(tree->mod->r_mat->n_diff_rr);

        for (i = 0; i < 5; i++)
          if (permut[i] != 5)
          {
            Generic_Brent_Lk(&(tree->mod->m4mod->o_rr->v[permut[i]]), 1.E-4, 1.E+4,
                             tree->mod->s_opt->min_diff_lk_local,
                             tree->mod->s_opt->brent_it_max,
                             tree->mod->s_opt->quickdirty, Wrap_Lk, NULL, tree,
                             NULL, NO, NO);
          }

        Free(permut);

        if (verbose) Print_Lk(tree, "[GTR parameters     ]");

        Set_Update_Eigen(NO, tree->mod);
      }    
  }
}

#define ITMAX 200
#define EPS   3.0e-8
#define TOLX (4*EPS)
#define STPMX 100.0
static phydbl sqrarg;
#define SQR(a) ((sqrarg=(a)) < SMALL ? 0.0 : sqrarg*sqrarg)

void BFGS(t_tree *tree,
          phydbl *p,
          int n,
          phydbl gtol,
          phydbl difff,
          phydbl step_size,
          short int logt, short int expt, short int is_positive,
          phydbl(*func)(t_tree *tree),
          int(*dfunc)(t_tree *tree,phydbl *param,int n_param,phydbl stepsize, short int logt, short int expt, phydbl(*func)(t_tree *tree),phydbl *derivatives, int is_positive),
          int(*lnsrch)(t_tree *tree, int n, phydbl *xold, phydbl fold, phydbl *g, phydbl *p, phydbl *x,phydbl *f, phydbl stpmax, int *check, short int logt, short int expt, int is_positive),
          int *failed)
{
  
  int check,i,its,j;
  phydbl den,fac,fad,fae,fp,stpmax,sum=0.0,sumdg,sumxi,temp,test,fret;
  phydbl *dg,*g,*hdg,**hessin,*pnew,*xi;
  phydbl fp_old;
  phydbl *init;

  hessin = (phydbl **)mCalloc(n,sizeof(phydbl *));
  for(i=0;i<n;i++) hessin[i] = (phydbl *)mCalloc(n,sizeof(phydbl));
  dg   = (phydbl *)mCalloc(n,sizeof(phydbl ));
  g    = (phydbl *)mCalloc(n,sizeof(phydbl ));
  pnew = (phydbl *)mCalloc(n,sizeof(phydbl ));
  hdg  = (phydbl *)mCalloc(n,sizeof(phydbl ));
  xi   = (phydbl *)mCalloc(n,sizeof(phydbl ));
  init = (phydbl *)mCalloc(n,sizeof(phydbl ));


  for(i=0;i<n;i++) init[i] = p[i];

  if(expt == YES) for(i=0;i<n;i++) p[i] = log(p[i]);
  if(logt == YES) for(i=0;i<n;i++) p[i] = exp(MIN(1.E+2,p[i]));
  fp=(*func)(tree);
  if(logt == YES) for(i=0;i<n;i++) p[i] = log(p[i]);
  if(expt == YES) for(i=0;i<n;i++) p[i] = exp(p[i]);
  
  /* PhyML_Printf("\n. ENTER BFGS WITH: %f\n",fp); */
  
  fp_old = fp;
  
  (*dfunc)(tree,p,n,step_size,logt,expt,func,g,is_positive);
  
  /* PhyML_Printf("\n. BFGS step_size: %f",step_size); */

  for (i=0;i<n;i++)
    {
      for (j=0;j<n;j++) hessin[i][j]=0.0;
      hessin[i][i]=1.0;
      xi[i] = -g[i];
      sum += p[i]*p[i];
      /* PhyML_Printf("\n. BFGS x[%d]: %f p[i]: %f",i,xi[i],p[i]); */
    }
  
  stpmax=STPMX*MAX(SQRT(sum),(phydbl)n);
  for(its=1;its<=ITMAX;its++)
    {
      
      lnsrch(tree,n,p,fp,g,xi,pnew,&fret,stpmax,&check,logt,expt,is_positive);
      
      /* PhyML_Printf("\n. BFGS: %f stpmax: %f fret: %g\n",tree->c_lnL,stpmax,fret); */

      fp_old = fp;
      fp = fret;
      
      for (i=0;i<n;i++)
        {
          xi[i]=pnew[i]-p[i];
          p[i]=pnew[i];
        }
      
      test=0.0;
      for (i=0;i<n;i++)
        {
          temp=xi[i]/MAX(p[i],1.0);
          /* printf("\n. x[i]=%f p[i]=%f",xi[i],p[i]); */
          if (temp > test) test=temp;
        }
      
      if (test < TOLX || (fabs(fp-fp_old) < difff && its > 1))
        {
          if(fp > fp_old)
            {
              for(i=0;i<n;i++) p[i] = init[i];
              *failed = YES;
            }
          
          if(expt == YES) for(i=0;i<n;i++) p[i] = log(p[i]);
          if(logt == YES) for(i=0;i<n;i++) p[i] = exp(MIN(1.E+2,p[i]));
          (*func)(tree);
          if(logt == YES) for(i=0;i<n;i++) p[i] = log(p[i]);
          if(expt == YES) for(i=0;i<n;i++) p[i] = exp(p[i]);
                    
          for(i=0;i<n;i++) Free(hessin[i]);
          free(hessin);
          free(xi);
          free(pnew);
          free(hdg);
          free(g);
          free(dg);
          free(init);
          return;
        }
      
      for (i=0;i<n;i++) dg[i]=g[i];
      
      (*dfunc)(tree,p,n,step_size,logt,expt,func,g,is_positive);
      
      test=0.0;
      den=MAX(fret,1.0);
      for (i=0;i<n;i++)
        {
          temp=g[i]*MAX(p[i],1.0)/den;
          if (temp > test) test=temp;
        }
      if (test < gtol)
        {
          *failed = NO;
          if(expt == YES) for(i=0;i<n;i++) p[i] = log(p[i]);
          if(logt == YES) for(i=0;i<n;i++) p[i] = exp(MIN(1.E+2,p[i]));
          (*func)(tree);
          if(logt == YES) for(i=0;i<n;i++) p[i] = log(p[i]);
          if(expt == YES) for(i=0;i<n;i++) p[i] = exp(p[i]);
                    
          for(i=0;i<n;i++) Free(hessin[i]);
          free(hessin);
          free(xi);
          free(pnew);
          free(hdg);
          free(g);
          free(dg);
          free(init);
          return;
        }

    for (i=0;i<n;i++) dg[i]=g[i]-dg[i];

    for (i=0;i<n;i++)
      {
        hdg[i]=0.0;
        for (j=0;j<n;j++) hdg[i] += hessin[i][j]*dg[j];
      }
    
    fac=fae=sumdg=sumxi=0.0;
    for (i=0;i<n;i++)
      {
        fac += dg[i]*xi[i];
        fae += dg[i]*hdg[i];
        sumdg += SQR(dg[i]);
        sumxi += SQR(xi[i]);
      }
    
    if(fac*fac > EPS*sumdg*sumxi)
      {
        fac=1.0/fac;
        fad=1.0/fae;
        for (i=0;i<n;i++) dg[i]=fac*xi[i]-fad*hdg[i];
        for (i=0;i<n;i++)
          {
            for (j=0;j<n;j++)
              {
                hessin[i][j] += fac*xi[i]*xi[j]
                  -fad*hdg[i]*hdg[j]+fae*dg[i]*dg[j];
              }
          }
      }
    for (i=0;i<n;i++)
      {
        xi[i]=0.0;
        for (j=0;j<n;j++) xi[i] -= hessin[i][j]*g[j];
      }
    }
  /*   PhyML_Printf("\n. Too many iterations in BFGS...\n"); */
  *failed = YES;
  for(i=0;i<n;i++) Free(hessin[i]);
  free(hessin);
  free(xi);
  free(pnew);
  free(hdg);
  free(g);
  free(dg);
}

#undef ITMAX
#undef EPS
#undef TOLX
#undef STPMX

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


#define ITMAX 2000
#define EPS   3.0e-8
#define TOLX (4*EPS)
#define STPMX 100.0
static phydbl sqrarg;
#define SQR(a) ((sqrarg=(a)) < SMALL ? 0.0 : sqrarg*sqrarg)

void BFGS_Nonaligned(t_tree *tree,
                     phydbl **p,
                     int n,
                     phydbl gtol,
                     phydbl difff,
                     phydbl step_size,
                     short int logt, short int expt, short int is_positive,
                     phydbl(*func)(t_tree *tree),
                     int(*dfunc_nonaligned)(t_tree *tree,phydbl **param,int n_param,phydbl stepsize,short int logt,short int expt,phydbl(*func)(t_tree *tree),phydbl *derivatives, int is_positive),
                     int(*lnsrch_nonaligned)(t_tree *tree, int n, phydbl **xold, phydbl fold,phydbl *g, phydbl *p, phydbl *x,phydbl *f, phydbl stpmax, int *check, short int logt, short int expt, int is_positive),
                     int *failed)
{

  int check,i,its,j;
  phydbl den,fac,fad,fae,fp,stpmax,sum=0.0,sumdg,sumxi,temp,test,fret;
  phydbl *dg,*g,*hdg,**hessin,*pnew,*xi;
  phydbl fp_old,fp_init;
  phydbl *init,*sign;

  hessin = (phydbl **)mCalloc(n,sizeof(phydbl *));
  for(i=0;i<n;i++) hessin[i] = (phydbl *)mCalloc(n,sizeof(phydbl));
  dg   = (phydbl *)mCalloc(n,sizeof(phydbl ));
  g    = (phydbl *)mCalloc(n,sizeof(phydbl ));
  pnew = (phydbl *)mCalloc(n,sizeof(phydbl ));
  hdg  = (phydbl *)mCalloc(n,sizeof(phydbl ));
  xi   = (phydbl *)mCalloc(n,sizeof(phydbl ));
  init = (phydbl *)mCalloc(n,sizeof(phydbl ));
  sign = (phydbl *)mCalloc(n,sizeof(phydbl ));

  for(i=0;i<n;i++) init[i] = (*(p[i]));

  if(expt == YES) for(i=0;i<n;i++) (*(p[i])) = log(*(p[i]));
  if(logt == YES) for(i=0;i<n;i++) (*(p[i])) = exp(MIN(1.E+2,*(p[i])));
  fp=(*func)(tree);
  fp_init = fp;
  if(logt == YES) for(i=0;i<n;i++) (*(p[i])) = log(*(p[i]));
  if(expt == YES) for(i=0;i<n;i++) (*(p[i])) = exp(*(p[i]));

  fp_old = fp;

  /* PhyML_Printf("\n- ENTER BFGS WITH: %f\n",fp); */

  (*dfunc_nonaligned)(tree,p,n,step_size,logt,expt,func,g,is_positive);

  for (i=0;i<n;i++)
    {
      for (j=0;j<n;j++) hessin[i][j]=0.0;
      hessin[i][i]=1.0;
      xi[i] = -g[i];
      sum += (*(p[i]))*(*(p[i]));
    }


  stpmax=STPMX*MAX(SQRT(sum),(phydbl)n);
  /* stpmax = 0.01*MAX(SQRT(sum),(phydbl)n); */
  for(its=1;its<=ITMAX;its++)
    {
      lnsrch_nonaligned(tree,n,p,fp,g,xi,pnew,&fret,stpmax,&check,logt,expt,is_positive);

      /* PhyML_Printf("\n. BFGS -> %f\n",tree->c_lnL); */

      fp_old = fp;
      fp = fret;

      for (i=0;i<n;i++)
        {
          xi[i]=pnew[i]-(*(p[i]));
          (*(p[i]))=pnew[i];
        }

      test=0.0;
      for (i=0;i<n;i++)
        {
          temp=xi[i]/MAX(*(p[i]),1.0);
          /* printf("\n. x[i]=%G p[i]=%f",xi[i],*(p[i])); */
          if (temp > test) test=temp;
        }

      if (test < TOLX || (FABS(fp_old-fp) < difff && its > 1))
        {          
          if(fp > fp_init)
            {
              for(i=0;i<n;i++) (*(p[i])) = init[i];
              *failed = YES;
            }
          
          if(expt == YES) for(i=0;i<n;i++) (*(p[i])) = log(*(p[i]));
          if(logt == YES) for(i=0;i<n;i++) (*(p[i])) = exp(MIN(1.E+2,*(p[i])));
          for(i=0;i<n;i++) sign[i] = *(p[i]) > .0 ? 1. : -1.;
          if(is_positive == YES) for(i=0;i<n;i++) *(p[i]) = FABS(*(p[i]));
          (*func)(tree);
          if(is_positive == YES) for(i=0;i<n;i++) *(p[i]) *= sign[i];
          if(logt == YES) for(i=0;i<n;i++) (*(p[i])) = log(*(p[i]));
          if(expt == YES) for(i=0;i<n;i++) (*(p[i])) = exp(*(p[i]));
          
          if(is_positive == YES) for(i=0;i<n;i++) *(p[i]) = FABS(*(p[i]));
          
          if(fp < fp_init) (*failed) = NO;

          /* PhyML_Printf("\n0 EXIT BFGS WITH: %f fp_old: %f\n",(*func)(tree),fp_old); */

          for(i=0;i<n;i++) Free(hessin[i]);
          free(hessin);
          free(xi);
          free(pnew);
          free(hdg);
          free(g);
          free(dg);
          free(init);
          free(sign);
          return;
        }
      
      for (i=0;i<n;i++) dg[i]=g[i];
      
      (*dfunc_nonaligned)(tree,p,n,step_size,logt,expt,func,g,is_positive);
      
      test=0.0;
      den=MAX(fret,1.0);
      for (i=0;i<n;i++)
        {
          temp=g[i]*MAX(*(p[i]),1.0)/den;
          if (temp > test) test=temp;
        }
      
      if (test < gtol)
        {
          if(expt == YES) for(i=0;i<n;i++) (*(p[i])) = log(*(p[i]));
          if(logt == YES) for(i=0;i<n;i++) (*(p[i])) = exp(MIN(1.E+2,*(p[i])));
          for(i=0;i<n;i++) sign[i] = *(p[i]) > .0 ? 1. : -1.;
          if(is_positive == YES) for(i=0;i<n;i++) *(p[i]) = FABS(*(p[i]));
          (*func)(tree);
          if(is_positive == YES) for(i=0;i<n;i++) *(p[i]) *= sign[i];
          if(logt == YES) for(i=0;i<n;i++) (*(p[i])) = log(*(p[i]));
          if(expt == YES) for(i=0;i<n;i++) (*(p[i])) = exp(*(p[i]));
          
          if(is_positive == YES) for(i=0;i<n;i++) *(p[i]) = FABS(*(p[i]));
          
          for(i=0;i<n;i++) Free(hessin[i]);
          free(hessin);
          free(xi);
          free(pnew);
          free(hdg);
          free(g);
          free(dg);
          free(init);
          free(sign);

          if(fp < fp_init) (*failed) = NO;

          /* PhyML_Printf("\n1 EXIT BFGS WITH: %f\n",(*func)(tree)); */

          return;
        }
      
      for (i=0;i<n;i++) dg[i]=g[i]-dg[i];
      
      for (i=0;i<n;i++)
        {
          hdg[i]=0.0;
          for (j=0;j<n;j++) hdg[i] += hessin[i][j]*dg[j];
        }
      
      fac=fae=sumdg=sumxi=0.0;
      for (i=0;i<n;i++)
        {
          fac += dg[i]*xi[i];
          fae += dg[i]*hdg[i];
          sumdg += SQR(dg[i]);
          sumxi += SQR(xi[i]);
        }
      
      if(fac*fac > EPS*sumdg*sumxi)
        {
          fac=1.0/fac;
          fad=1.0/fae;
          for (i=0;i<n;i++) dg[i]=fac*xi[i]-fad*hdg[i];
          for (i=0;i<n;i++)
            {
              for (j=0;j<n;j++)
                {
                  hessin[i][j] += fac*xi[i]*xi[j]
                    -fad*hdg[i]*hdg[j]+fae*dg[i]*dg[j];
                }
            }
        }
      for (i=0;i<n;i++)
        {
          xi[i]=0.0;
          for (j=0;j<n;j++) xi[i] -= hessin[i][j]*g[j];
        }
    }

  if(fp < fp_init) (*failed) = NO;

  /* PhyML_Printf("\n2 EXIT BFGS WITH: %f\n",(*func)(tree)); */

  /* PhyML_Printf("\n. Too many iterations in BFGS...\n"); */
  *failed = YES;
  for(i=0;i<n;i++) Free(hessin[i]);
  free(hessin);
  free(xi);
  free(pnew);
  free(hdg);
  free(g);
  free(dg);
  free(sign);
}

#undef ITMAX
#undef EPS
#undef TOLX
#undef STPMX

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

#define ALF 1.0e-4
#define TOLX 1.0e-7

int Lnsrch(t_tree *tree, int n, phydbl *xold, phydbl fold,
           phydbl *g, phydbl *p, phydbl *x,
           phydbl *f, phydbl stpmax, int *check, short int logt, short int expt, int is_positive)
{
  int i;
  phydbl a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp,test,tmplam;
  phydbl *local_xold,*sign;
  
  alam = alam2 = f2 = fold2 = tmplam = .0;

  local_xold = (phydbl *)mCalloc(n,sizeof(phydbl));
  sign       = (phydbl *)mCalloc(n,sizeof(phydbl));

  for(i=0;i<n;i++) local_xold[i] = xold[i];
  
  *check=0;
  for(sum=0.0,i=0;i<n;i++) sum += p[i]*p[i];
  sum=SQRT(sum);
  /* PhyML_Printf("\n. lnsrch sum: %f",sum); */
  if(sum > stpmax) for(i=0;i<n;i++) p[i] *= stpmax/sum;
  /* for(i=0;i<n;i++) PhyML_Printf("\n. lnsrch p[i]: %f",p[i]); */
  slope=0.0;
  for(i=0;i<n;i++) slope += g[i]*p[i];
  /* PhyML_Printf("\n. lnsrch slope: %f",slope); */
  test=0.0;
  for(i=0;i<n;i++)
    {
      temp=p[i]/MAX(local_xold[i],1.0);
      if (temp > test) test=temp;
    }
  alamin=TOLX/test;
  alam=1.0;
  for (;;)
    {
      for(i=0;i<n;i++)
        {
          x[i]=local_xold[i]+alam*p[i];
          xold[i] = x[i];
          /* PhyML_Printf("\n. lnsrch x[i]: %f",x[i]); */
        }
      
      /* PhyML_Printf("\n. lnsrch loop slope: %f alam: %f alam2: %f",slope,alam,alam2); */
      
      if(i==n)
        {
          if(expt == YES) for(i=0;i<n;i++) xold[i] = log(xold[i]);
          if(logt == YES) for(i=0;i<n;i++) xold[i] = exp(MIN(1.E+2,xold[i]));
          for(i=0;i<n;i++) sign[i] = xold[i] < .0 ? -1. : 1.;
          if(is_positive == YES) for(i=0;i<n;i++) xold[i] = FABS(xold[i]);
          /* for(i=0;i<n;i++) PhyML_Printf("\n. <<>> %f",xold[i]); */
          *f=Return_Abs_Lk(tree);
          if(is_positive == YES) for(i=0;i<n;i++) xold[i] *= sign[i];
          if(logt == YES) for(i=0;i<n;i++) xold[i] = log(xold[i]);
          if(expt == YES) for(i=0;i<n;i++) xold[i] = exp(xold[i]);
        }
      else *f=1.+fold+ALF*alam*slope;
      if (alam < alamin)
        {
          *check=1;
          for(i=0;i<n;i++) xold[i] = local_xold[i];
          if(is_positive == YES) for(i=0;i<n;i++) xold[i] = FABS(xold[i]);
          Free(local_xold);
          Free(sign);
          return 0;
        }
      else if (*f <= fold+ALF*alam*slope)
        {
          for(i=0;i<n;i++) xold[i] = local_xold[i];
          if(is_positive == YES) for(i=0;i<n;i++) xold[i] = FABS(xold[i]);
          Free(local_xold);
          Free(sign);
          return 0;
        }
      else
        {
          /* 	  if (alam == 1.0) */
          if ((alam < 1.0+SMALL) && (alam > 1.0-SMALL))
            tmplam = -slope/(2.0*(*f-fold-slope));
          else
            {
              rhs1 = *f-fold-alam*slope;
              rhs2=f2-fold2-alam2*slope;
              a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
              b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
              if (a < SMALL && a > -SMALL) tmplam = -slope/(2.0*b);
              else
                {
                  disc=b*b-3.0*a*slope;
                  if (disc<0.0) tmplam = 0.5*alam;
                  else if(b <= 0.0) tmplam=(-b+SQRT(disc))/(3.0*a);
                  else tmplam = -slope/(b+SQRT(disc));
                }
              if (tmplam>0.5*alam) tmplam=0.5*alam;
            }
        }
      alam2=alam;
      f2 = *f;
      fold2=fold;
      alam=MAX(tmplam,0.1*alam);
    }
  return 1;
}

#undef ALF
#undef TOLX
#undef NRANSI

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

#define ALF 1.0e-4
#define TOLX 1.0e-7

int Lnsrch_Nonaligned(t_tree *tree, int n, phydbl **xold, phydbl fold,
                      phydbl *g, phydbl *p, phydbl *x,
                      phydbl *f, phydbl stpmax, int *check, short int logt, short int expt, int is_positive)
{
  int i;
  phydbl a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp,test,tmplam;
  phydbl *local_xold,*sign;


  alam = alam2 = f2 = fold2 = tmplam = .0;

  local_xold = (phydbl *)mCalloc(n,sizeof(phydbl));
  sign       = (phydbl *)mCalloc(n,sizeof(phydbl));

  for(i=0;i<n;i++) local_xold[i] = *(xold[i]);

  *check=0;
  for(sum=0.0,i=0;i<n;i++) sum += p[i]*p[i];
  sum=SQRT(sum);
  if(sum > stpmax)
    for(i=0;i<n;i++) p[i] *= stpmax/sum;
  for(slope=0.0,i=0;i<n;i++)
    slope += g[i]*p[i];
  test=0.0;
  for(i=0;i<n;i++)
    {
      temp=p[i]/MAX(local_xold[i],1.0);
      if (temp > test) test=temp;
    }
  alamin=TOLX/test;
  alam=1.0;
  for (;;)
    {
      for(i=0;i<n;i++)
        {
          x[i]=local_xold[i]+alam*p[i];
          *(xold[i]) = x[i];
        }
      
      if(i==n)
        {
          if(expt == YES) for(i=0;i<n;i++) *(xold[i]) = log(*(xold[i]));
          if(logt == YES) for(i=0;i<n;i++) *(xold[i]) = exp(MIN(1.E+2,*(xold[i])));
          for(i=0;i<n;i++) sign[i]    = *(xold[i]) < .0 ? -1. : 1.;
          if(is_positive == YES) for(i=0;i<n;i++) *(xold[i]) = FABS(*(xold[i]));
          *f=Return_Abs_Lk(tree);
          if(is_positive == YES) for(i=0;i<n;i++) *(xold[i]) *= sign[i];
          if(logt == YES) for(i=0;i<n;i++) *(xold[i]) = log(*(xold[i]));
          if(expt == YES) for(i=0;i<n;i++) *(xold[i]) = exp(*(xold[i]));
        }
      else *f=1.+fold+ALF*alam*slope;
      
      if (alam < alamin)
        {
          *check=1;
          for(i=0;i<n;i++) *(xold[i]) = local_xold[i];
          if(is_positive == YES) for(i=0;i<n;i++) *(xold[i]) = FABS(*(xold[i]));
          Free(local_xold);
          Free(sign);
          return 0;
        }
      else if (*f <= fold+ALF*alam*slope)
        {
          for(i=0;i<n;i++) *(xold[i]) = local_xold[i];
          if(is_positive) for(i=0;i<n;i++) *(xold[i]) = FABS(*(xold[i]));
          Free(local_xold);
          Free(sign);
          return 0;
        }
      else
        {
          /* 	  if (alam == 1.0) */
          if ((alam < 1.0+SMALL) && (alam > 1.0-SMALL))
            tmplam = -slope/(2.0*(*f-fold-slope));
          else
            {
              rhs1 = *f-fold-alam*slope;
              rhs2=f2-fold2-alam2*slope;
              a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
              b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
              if (a < SMALL && a > -SMALL) tmplam = -slope/(2.0*b);
              else
                {
                  disc=b*b-3.0*a*slope;
                  if (disc<0.0) tmplam = 0.5*alam;
                  else if(b <= 0.0) tmplam=(-b+SQRT(disc))/(3.0*a);
                  else tmplam = -slope/(b+SQRT(disc));
                }
              if (tmplam>0.5*alam) tmplam=0.5*alam;
            }
        }
      alam2=alam;
      f2 = *f;
      fold2=fold;
      alam=MAX(tmplam,0.1*alam);
    }
  Free(local_xold);
  Free(sign);
  return 1;
}

#undef ALF
#undef TOLX
#undef NRANSI

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int Dist_F_Brak(phydbl *ax, phydbl *bx, phydbl *cx, phydbl *F, phydbl *param, t_mod *mod)
{
   phydbl ulim,u,r,q,dum;
   phydbl fa, fb, fc, fu;

   fa = -Lk_Dist(F,FABS(*ax),mod);
   fb = -Lk_Dist(F,FABS(*bx),mod);

   if(fb > fa)
     {
       SHFT(dum,*ax,*bx,dum)
       SHFT(dum,fb,fa,dum)
     }

   *cx=(*bx)+MNBRAK_GOLD*(*bx-*ax);
   fc = -Lk_Dist(F,FABS(*cx),mod);

   while (fb > fc)
     {
       r=(*bx-*ax)*(fb-fc);
       q=(*bx-*cx)*(fb-fa);
       u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
               (2.0*SIGN(MAX(FABS(q-r),MNBRAK_TINY),q-r));
       ulim=(*bx)+MNBRAK_GLIMIT*(*cx-*bx);

       if ((*bx-u)*(u-*cx) > 0.0)
     {
       fu = -Lk_Dist(F,FABS(u),mod);
       if (fu < fc)
         {
           *ax=(*bx);
           *bx=u;
           fa=fb;
           fb=fu;
           return(0);
         }
       else if (fu > fb)
         {
           *cx=u;
           fc=fu;
           return(0);
         }
       u=(*cx)+MNBRAK_GOLD*(*cx-*bx);
       fu = -Lk_Dist(F,FABS(u),mod);
     }
       else if ((*cx-u)*(u-ulim) > 0.0)
     {
       fu = -Lk_Dist(F,FABS(u),mod);
       if (fu < fc)
         {
           SHFT(*bx,*cx,u,*cx+MNBRAK_GOLD*(*cx-*bx))
           SHFT(fb,fc,fu,-Lk_Dist(F,FABS(u),mod))
         }
     }
       else if ((u-ulim)*(ulim-*cx) >= 0.0)
     {
       u  = ulim;
       fu = -Lk_Dist(F,FABS(u),mod);
     }
       else
     {
       u  =(*cx)+MNBRAK_GOLD*(*cx-*bx);
       fu = -Lk_Dist(F,FABS(u),mod);
     }

       SHFT(*ax,*bx,*cx,u)
       SHFT(fa,fb,fc,fu)
      }
   return(0);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl Dist_F_Brent(phydbl ax, phydbl bx, phydbl cx, phydbl tol, int n_iter_max,
                    phydbl *param, phydbl *F, t_mod *mod)
{
  int iter;
  phydbl a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  phydbl e=0.0;
  phydbl old_lnL,init_lnL, curr_lnL;
  
  d=0.0;
  a=((ax < cx) ? ax : cx);
  b=((ax > cx) ? ax : cx);
  x = w = v = bx;
  old_lnL = UNLIKELY;
  fw = fv = fx = -Lk_Dist(F,fabs(bx),mod);
  curr_lnL = init_lnL = -fw;
  
  /* PhyML_Printf("\n. bx=%f f: %f %f %f %f fx: %f",bx,mod->e_frq->pi->v[0],mod->e_frq->pi->v[1],mod->e_frq->pi->v[2],mod->e_frq->pi->v[3],fx); */
  assert(isnan(fx) == FALSE);
  assert(isinf(fx) == FALSE);

  for(iter=1;iter<=BRENT_IT_MAX;iter++)
    {
      xm=0.5*(a+b);
      
      tol2=2.0*(tol1=tol*FABS(x)+BRENT_ZEPS);
      
      if(
         ((FABS(curr_lnL-old_lnL) < mod->s_opt->min_diff_lk_local) &&
          (curr_lnL > init_lnL - mod->s_opt->min_diff_lk_local)) ||
         (iter > n_iter_max - 1)
         )
        {
          *param = x;
          curr_lnL = Lk_Dist(F,*param,mod);
          return -curr_lnL;
        }
      
      if(FABS(e) > tol1)
        {
          r=(x-w)*(fx-fv);
          q=(x-v)*(fx-fw);
          p=(x-v)*q-(x-w)*r;
          q=2.0*(q-r);
          if(q > 0.0) p = -p;
          q=FABS(q);
          etemp=e;
          e=d;
          if(FABS(p) >= FABS(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
            {
              d=BRENT_CGOLD*(e=(x >= xm ? a-x : b-x));
              /* PhyML_Printf("Golden section step\n"); */
            }
          else
            {
              d=p/q;
              u=x+d;
              if (u-a < tol2 || b-u < tol2)
                d=SIGN(tol1,xm-x);
              /* PhyML_Printf("Parabolic step\n"); */
            }
        }
      else
        {
          d=BRENT_CGOLD*(e=(x >= xm ? a-x : b-x));
          /* PhyML_Printf("Golden section step (default)\n"); */
        }
      
      u=(FABS(d) >= tol1 ? x+d : x+SIGN(tol1,d));
      (*param) = FABS(u);
      old_lnL = curr_lnL;
      fu = -Lk_Dist(F,FABS(u),mod);
      curr_lnL = -fu;
      /* PhyML_Printf("param=%f loglk=%f\n",*param,fu); */
      
      /*       if(fu <= fx)  */
      if(fu < fx)
        {
          if(iter > n_iter_max) return -fu;
          
          if(u >= x) a=x; else b=x;
          SHFT(v,w,x,u)
            SHFT(fv,fw,fx,fu)
            }
      else
        {
          if (u < x) a=u; else b=u;
          /* 	  if (fu <= fw || w == x)  */
          if (fu < fw || FABS(w-x) < SMALL)
            {
              v=w;
              w=u;
              fv=fw;
              fw=fu;
            }
          /* 	  else if (fu <= fv || v == x || v == w)  */
          else if (fu < fv || FABS(v-x) < SMALL || FABS(v-w) < SMALL)
            {
              v=u;
              fv=fu;
            }
        }
    }
  
  Exit("\n. Too many iterations in Dist_F_Brent !\n");
  return(-1);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Opt_Dist_F(phydbl *dist, phydbl *F, t_mod *mod)
{
  phydbl ax,bx,cx;

  if(*dist < mod->l_min) *dist = mod->l_min;

  ax = mod->l_min;
  bx =  (*dist);
  cx = mod->l_max;

  /* PhyML_Printf("\n. bx: %g",bx); */
  
/*   Dist_F_Brak(&ax,&bx,&cx,F,dist,mod); */
  Dist_F_Brent(ax,bx,cx,1.E-10,1000,dist,F,mod);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


int Missing_Dist_Brak(phydbl *ax, phydbl *bx, phydbl *cx, int x, int y, matrix *mat)
{
   phydbl ulim,u,r,q,dum;
   phydbl fa, fb, fc, fu;

   fa = Least_Square_Missing_Dist_XY(x,y,FABS(*ax),mat);
   fb = Least_Square_Missing_Dist_XY(x,y,FABS(*bx),mat);

   if(fb > fa)
     {
       SHFT(dum,*ax,*bx,dum)
       SHFT(dum,fb,fa,dum)
     }

   *cx=(*bx)+MNBRAK_GOLD*((*bx)-(*ax));
   fc = Least_Square_Missing_Dist_XY(x,y,FABS(*cx),mat);

   while (fb > fc)
     {
       r=((*bx)-(*ax))*(fb-fc);
       q=((*bx)-(*cx))*(fb-fa);
       u=(*bx)-(((*bx)-(*cx))*q-((*bx)-(*ax))*r)/
               (2.0*SIGN(MAX(FABS(q-r),MNBRAK_TINY),q-r));
       ulim=(*bx)+MNBRAK_GLIMIT*(*cx-*bx);

       if ((*bx-u)*(u-*cx) > 0.0)
     {
       fu = Least_Square_Missing_Dist_XY(x,y,FABS(u),mat);
       if (fu < fc)
         {
           *ax=(*bx);
           *bx=u;
           fa=fb;
           fb=fu;
           return(0);
         }
       else if (fu > fb)
         {
           *cx=u;
           fc=fu;
           return(0);
         }
       u=(*cx)+MNBRAK_GOLD*(*cx-*bx);
       fu = Least_Square_Missing_Dist_XY(x,y,FABS(u),mat);
     }
       else if ((*cx-u)*(u-ulim) > 0.0)
     {
       fu = Least_Square_Missing_Dist_XY(x,y,FABS(u),mat);
       if (fu < fc)
         {
           SHFT(*bx,*cx,u,*cx+MNBRAK_GOLD*(*cx-*bx))
         SHFT(fb,fc,fu,Least_Square_Missing_Dist_XY(x,y,FABS(u),mat))
         }
     }
       else if ((u-ulim)*(ulim-*cx) >= 0.0)
     {
       u  = ulim;
       fu = Least_Square_Missing_Dist_XY(x,y,FABS(u),mat);
     }
       else
     {
       u  =(*cx)+MNBRAK_GOLD*(*cx-*bx);
       fu = Least_Square_Missing_Dist_XY(x,y,FABS(u),mat);
     }

       SHFT(*ax,*bx,*cx,u)
       SHFT(fa,fb,fc,fu)
      }
   return(0);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Missing_Dist_Brent(phydbl ax, phydbl bx, phydbl cx, phydbl tol, int n_iter_max,
              int x, int y, matrix *mat)
{
  int iter;
  phydbl a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,xx,xm;
  phydbl e=0.0;

  d=0.0;
  a=((ax < cx) ? ax : cx);
  b=((ax > cx) ? ax : cx);
  xx=w=v=bx;
  fx=Least_Square_Missing_Dist_XY(x,y,FABS(bx),mat);
  fw=fv=-fx;

  for(iter=1;iter<=BRENT_IT_MAX;iter++)
    {
      xm=0.5*(a+b);
      tol2=2.0*(tol1=tol*FABS(xx)+BRENT_ZEPS);

      if(FABS(xx-xm) <= (tol2-0.5*(b-a)))
    {
      mat->dist[x][y] = xx;
      Least_Square_Missing_Dist_XY(x,y,mat->dist[x][y],mat);
      return -fx;
    }

      if(FABS(e) > tol1)
    {
      r=(xx-w)*(fx-fv);
      q=(xx-v)*(fx-fw);
      p=(xx-v)*q-(xx-w)*r;
      q=2.0*(q-r);
      if(q > 0.0) p = -p;
      q=FABS(q);
      etemp=e;
      e=d;
      if(FABS(p) >= FABS(0.5*q*etemp) || p <= q*(a-xx) || p >= q*(b-xx))
        {
          d=BRENT_CGOLD*(e=(xx >= xm ? a-xx : b-xx));
          /*                   PhyML_Printf("Golden section step\n"); */
        }
      else
        {
          d=p/q;
          u=xx+d;
          if (u-a < tol2 || b-u < tol2)
        d=SIGN(tol1,xm-xx);
          /*                   PhyML_Printf("Parabolic step\n"); */
        }
        }
      else
    {
      d=BRENT_CGOLD*(e=(xx >= xm ? a-xx : b-xx));
      /*               PhyML_Printf("Golden section step (default)\n"); */
    }

      u=(FABS(d) >= tol1 ? xx+d : xx+SIGN(tol1,d));
      fu = Least_Square_Missing_Dist_XY(x,y,FABS(u),mat);

/*       PhyML_Printf("param=%f loglk=%f\n",u,fu); */

/*       if(fu <= fx)  */
      if(fu < fx)
    {
      if(iter > n_iter_max) return -fu;

      if(u >= xx) a=xx; else b=xx;
      SHFT(v,w,xx,u)
      SHFT(fv,fw,fx,fu)
    }
      else
    {
      if (u < xx) a=u; else b=u;
/* 	  if (fu <= fw || w == xx)  */
      if (fu < fw || FABS(w-xx) < SMALL)
        {
          v=w;
          w=u;
          fv=fw;
          fw=fu;
        }
/* 	  else if (fu <= fv || v == xx || v == w)  */
      else if (fu < fv || FABS(v-xx) < SMALL || FABS(v-w) < SMALL)
        {
          v=u;
          fv=fu;
        }
    }
    }
  Exit("\n. Too many iterations in Missing_Dist_Brent !");
  return(-1);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Opt_Missing_Dist(int x, int y, matrix *mat)
{
  phydbl ax,bx,cx;

  ax = DIST_MAX;
  bx = DIST_MAX/4.;

  Missing_Dist_Brak(&ax,&bx,&cx,x,y,mat);
  PhyML_Printf("ax=%f bx=%f cx=%f\n",FABS(ax),FABS(bx),FABS(cx));
  Missing_Dist_Brent(FABS(ax),FABS(bx),FABS(cx),1.E-5,100,x,y,mat);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int Optimiz_Alpha_And_Pinv(t_tree *mixt_tree, int verbose)
{
  scalar_dbl **alpha;
  int n_alpha;
  t_tree *tree;
  int i;

  /* PhyML_Printf("\n. Optimizing Alpha and Inv"); */
  
  Set_Update_Eigen(NO,mixt_tree->mod);

  alpha   = NULL;
  n_alpha = 0;
  tree    = mixt_tree;

  do
    {
      if(tree->mod->ras->alpha->optimize == YES && tree->mod->ras->n_catg > 1 && tree->mod->ras->pinvar->optimize == YES)
        {
          for(i=0;i<n_alpha;i++) if(tree->mod->ras->alpha == alpha[i]) break;

          if(i == n_alpha)
            {
              if(!alpha) alpha = (scalar_dbl **)mCalloc(1,sizeof(scalar_dbl *));
              else       alpha = (scalar_dbl **)mRealloc(alpha,n_alpha+1,sizeof(scalar_dbl *));
              alpha[n_alpha] = tree->mod->ras->alpha;
              n_alpha++;

              if(tree->mod->ras->alpha->optimize == YES &&
                 tree->mod->ras->free_mixt_rates == NO)
                {
                  if(tree->mod->ras->n_catg > 1)
                    {
                      Generic_Brent_Lk(&(tree->mod->ras->alpha->v),
                                       ALPHA_MIN,ALPHA_MAX,
                                       tree->mod->s_opt->min_diff_lk_local,
                                       tree->mod->s_opt->brent_it_max,
                                       tree->mod->s_opt->quickdirty,
                                       Wrap_Lk,NULL,mixt_tree,NULL,NO,NO);
                    }
                  if(verbose == YES)
                    {
                      Print_Lk(mixt_tree,"[Alpha              ]");
                      PhyML_Printf("[%10f]",tree->mod->ras->alpha->v);
                    }
                }

              if(tree->mod->ras->pinvar->optimize == YES &&
                 tree->mod->ras->free_mixt_rates == NO)
                {
                  tree->mod->s_opt->skip_tree_traversal = YES;

                  Optimize_Single_Param_Generic(mixt_tree,&(tree->mod->ras->pinvar->v),.0001,0.9999,
                                                tree->mod->s_opt->min_diff_lk_local,
                                                tree->mod->s_opt->brent_it_max,
                                                tree->mod->s_opt->quickdirty);

                  tree->mod->s_opt->skip_tree_traversal = NO;

                  if(verbose == YES)
                    {
                      Print_Lk(mixt_tree,"[P-inv              ]");
                      PhyML_Printf("[%10f]",tree->mod->ras->pinvar->v);
                    }
                }
            }
        }
      tree = tree->next_mixt;
    }
  while(tree);

  if(alpha) Free(alpha);

  return 1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

static phydbl Br_Len_Spline(phydbl *l, t_edge *b, int n_iter_max, phydbl tol, t_tree *tree)
{
  short int converged;
  phydbl init_lnL,old_lnL;
  int iter;
  phydbl best_l, new_l, init_l, init_dl, best_lnL;
  phydbl u, v;
  phydbl fu, fv;
  phydbl dfu, dfv;
  phydbl mult;
  phydbl a_,b_,A_,B_,C_,D_,root1,root2;
  short int ok1, ok2;
  // Warning: make sure eigen_lr vectors are already up-to-date 

  Set_Use_Eigen_Lr(YES,tree);
  
  best_l = init_l = *l;
  best_lnL = old_lnL = init_lnL = tree->c_lnL;  
  mult = 1.2;
  ok1 = ok2 = NO;
  a_ = b_ = A_ = B_ = D_ = root1 = root2 = -1.;
  u = v = fu = fv = dfu = dfv = -1.;
  new_l = -1.;
  
  dLk(l,b,tree);
  init_dl = tree->c_dlnL;
  
  if(*l > tree->mod->l_max) *l = 0.5;
  if(*l < tree->mod->l_min) *l = 0.001;
  
  // Find value of l where first derivative is < 0;
  tree->c_dlnL = init_dl;
  while(tree->c_dlnL < 0.0)
    {
      *l /= mult;
      tree->n_tot_bl_opt++;
      if(*l < tree->mod->l_min)
        {
          *l = best_l;
          tree->c_lnL = best_lnL;
          return best_lnL;
        }
      dLk(l,b,tree);
      if(tree->c_lnL > best_lnL)
        {
          best_lnL = tree->c_lnL;
          best_l   = *l;
        }
    }
  u = *l;
  fu = tree->c_lnL;
  dfu = tree->c_dlnL;
  
  
  
  // Find good upper bound
  *l = init_l;
  tree->c_dlnL = init_dl;
  tree->c_lnL = init_lnL;

  while(tree->c_dlnL > 0.0)
    {
      *l *= mult;
      tree->n_tot_bl_opt++;
      if(*l > tree->mod->l_max)
        {
          *l = best_l;
          tree->c_lnL = best_lnL;
          return best_lnL;
        }
      dLk(l,b,tree);
      if(tree->c_lnL > best_lnL)
        {
          best_lnL = tree->c_lnL;
          best_l   = *l;
        }
    }
  v = *l;
  fv = tree->c_lnL;
  dfv = tree->c_dlnL;


  
  /* PhyML_Printf("\n Begin NR loop (lnL: %12G dlnL: %12G) l: %12G num: %d",tree->c_lnL,tree->c_dlnL,*l,b->num); */

  converged = NO;
  iter = 0;
  do
    {
      /* PhyML_Printf("\n. l=%12f lnL=%12f iter:%d u=%12f v=%12f root1=%12f root2=%12f dfu=%12f dfv=%12f fu=%12f fv=%12f diff=%12f tol=%12f init=%12f", */
      /*              *l,tree->c_lnL,iter, */
      /*              u,v, */
      /*              root1,root2, */
      /*              dfu,dfv, */
      /*              fu,fv, */
      /*              tree->c_lnL-old_lnL, */
      /*              tol,init_l); */

      /* PhyML_Printf("\n. l: %f dlnL: %f lnL: %f",*l,tree->c_dlnL,tree->c_lnL); */
      // Spline interpolation (https://en.wikipedia.org/wiki/Spline_interpolation)
      a_ = dfu*(v-u) - (fv-fu);
      b_ = -dfv*(v-u) + (fv-fu);
      
      A_ = 3.*a_ - 3.*b_;
      B_ = -4.*a_ + 2.*b_;
      C_ = fv-fu+a_;
      D_ = sqrt(B_*B_-4.*A_*C_);

      root1 = (-B_-D_)/(2.*A_);
      root2 = (-B_+D_)/(2.*A_);

      root1 = root1*(v-u) + u;
      root2 = root2*(v-u) + u;
      
      ok1 = NO;
      ok2 = NO;
      if(root1 > u && root1 < v) ok1 = YES;
      if(root2 > u && root2 < v) ok2 = YES;

      if(Are_Equal(root1,u,1.E-5) == YES) ok1 = YES;
      if(Are_Equal(root2,u,1.E-5) == YES) ok2 = YES;
      if(Are_Equal(root1,v,1.E-5) == YES) ok1 = YES;
      if(Are_Equal(root2,v,1.E-5) == YES) ok2 = YES;
      

      if(ok1 == YES && ok2 == YES) new_l = root1 < root2 ? root1 : root2;
      else if(ok1 == YES) new_l = root1;
      else if(ok2 == YES) new_l = root2;
      else if(u/v > 1.1 || u/v < 0.9)
        {
          PhyML_Printf("\n. iter=%4d u=%12G fu=%12G dfu=%12G v=%12G fv=%12G dfv=%12G root1=%12G root2=%12G\n",iter,u,fu,dfu,v,fv,dfv,root1,root2);
          assert(FALSE);
        }
      

      *l = new_l;
      tree->n_tot_bl_opt++;
          
      old_lnL = tree->c_lnL;
      dLk(l,b,tree);
      if(tree->c_lnL > best_lnL)
        {
          best_lnL = tree->c_lnL;
          best_l   = *l;
        }
            
      if(tree->c_dlnL > 0.0)
        {
          u = new_l;
          fu = tree->c_lnL;
          dfu = tree->c_dlnL;
        }
      else
        {
          v = new_l;
          fv = tree->c_lnL;
          dfv = tree->c_dlnL;
        }



      if(u - v < DBL_MIN) converged = YES;
      if(fabs(tree->c_lnL-old_lnL) < tol) converged = YES;
      if(++iter == n_iter_max+20) converged = YES;
      if(iter >= n_iter_max) PhyML_Fprintf(stderr,"\n. Edge length optimization took too long... l=%G lnL=%G iter:%d u=%G v=%G root1=%G root2=%G dfu=%G dfv=%G fu=%G fv=%G diff=%G tol=%G",
                                           *l,tree->c_lnL,iter,
                                           u,v,
                                           root1,root2,
                                           dfu,dfv,
                                           fu,fv,
                                           tree->c_lnL-old_lnL,
                                           tol);
      
      if(converged == NO)
        {
          if(!(u < v)) PhyML_Printf("\n. u=%g v=%g.\n",u,v);
          if(!(dfu > 0.0)) PhyML_Printf("\n. dfu=%g l=%g u=%g v=%g\n",dfu,*l,u,v);
          if(!(dfv < 0.0)) PhyML_Printf("\n. dfv=%g l=%g u=%g v=%g\n",dfv,*l);
          
          assert(u < v);
          assert(dfu > 0.0);
          assert(dfv < 0.0);
        }
    }
  while(converged == NO);

  /* PhyML_Printf("\n. l = %g l_max = %g diff=%g",*l,tree->mod->l_max,*l-tree->mod->l_max); */
  
  assert(!(*l > tree->mod->l_max)); 
  assert(!(*l < tree->mod->l_min)); 

  /* if(dfu > 1.) */
  /*   { */
  /*     PhyML_Printf("\n> [%4d] l=%f lnL=%f iter:%d u=%f v=%f root1=%f root2=%f dfu=%f dfv=%f fu=%f fv=%f diff=%f tol=%f init=%f", */
  /*                  tree->n_tot_bl_opt, */
  /*                  *l,tree->c_lnL,iter, */
  /*                  u,v, */
  /*                  root1,root2, */
  /*                  dfu,dfv, */
  /*                  fu,fv, */
  /*                  tree->c_lnL-old_lnL, */
  /*                  tol,init_l); */
  /*   } */
  /* if(dfv < -1) */
  /*   { */
  /*     PhyML_Printf("\n< [%4d] l=%f lnL=%f iter:%d u=%f v=%f root1=%f root2=%f dfu=%f dfv=%f fu=%f fv=%f diff=%f tol=%f init=%f", */
  /*                  tree->n_tot_bl_opt, */
  /*                  *l,tree->c_lnL,iter, */
  /*                  u,v, */
  /*                  root1,root2, */
  /*                  dfu,dfv, */
  /*                  fu,fv, */
  /*                  tree->c_lnL-old_lnL, */
  /*                  tol,init_l); */
  /*   } */

  *l = best_l;
  tree->c_lnL = best_lnL;
  
  if(iter == n_iter_max)
    {
      PhyML_Printf("\n. Too many iterations in edge length optimization routine (l=%G init=%G).\n",best_l,init_l);
      assert(FALSE);
    }
  
  return tree->c_lnL;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int Generic_Brent_Lk(phydbl *param, phydbl ax, phydbl cx, phydbl tol,
                     int n_iter_max, int quickdirty,
                     phydbl (*obj_func)(t_edge *,t_tree *,supert_tree *),
                     t_edge *branch, t_tree *tree, supert_tree *stree, short int logt, short int expt)
{
  int iter;
  phydbl a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  phydbl e=0.0;
  phydbl old_lnL,init_lnL;
  phydbl bx = *param;
  int n_opt_step;
  
  n_opt_step = 0;
  d=0.0;
  a=((ax < cx) ? ax : cx);
  b=((ax > cx) ? ax : cx);
  x=w=v=bx;
  (*param) = bx;
  if(expt == YES) (*param) = log(*param);
  if(logt == YES) (*param) = exp(MIN(1.E+2,*param));
  fw=fv=fx=fu=-(*obj_func)(branch,tree,stree);
  if(logt == YES) (*param) = log(*param);
  if(expt == YES) (*param) = exp(*param);
  init_lnL = old_lnL = fw;
  
  /* PhyML_Printf("\n. %p %p %p init_lnL=%f fu=%f ax=%f cx=%f param=%f",branch,tree,stree,init_lnL,fu,ax,cx,*param); */
  
  for(iter=1;iter<=BRENT_IT_MAX;iter++)
    {
      xm=0.5*(a+b);
      tol2=2.0*(tol1=tol*x+BRENT_ZEPS);
      
      if((fu < init_lnL + tol) && (quickdirty == YES) && (iter > 1))
        {
          (*param) = x;
          if(expt == YES) (*param) = log(*param);
          if(logt == YES) (*param) = exp(MIN(1.E+2,*param));
          fu = (*obj_func)(branch,tree,stree);
          if(logt == YES) (*param) = log(*param);
          if(expt == YES) (*param) = exp(*param);
          /* printf("\n. return %f [%f] %d",fu,*param,iter); */
          return n_opt_step;
        }

      if((FABS(fu-old_lnL) < tol && iter > 1) || (iter > n_iter_max - 1))
        {
          (*param) = x;
          if(expt == YES) (*param) = log(*param);
          if(logt == YES) (*param) = exp(MIN(1.E+2,*param));
          fu = (*obj_func)(branch,tree,stree);
          if(logt == YES) (*param) = log(*param);
          if(expt == YES) (*param) = exp(*param);
          /* printf("\n. return %f [%f] %d",*param,fu,iter); */
          return n_opt_step;
        }

      if(FABS(e) > tol1)
        {
          r=(x-w)*(fx-fv);
          q=(x-v)*(fx-fw);
          p=(x-v)*q-(x-w)*r;
          q=2.0*(q-r);
          if(q > 0.0) p = -p;
          q=FABS(q);
          etemp=e;
          e=d;
          if(FABS(p) >= FABS(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
            {
              d=BRENT_CGOLD*(e=(x >= xm ? a-x : b-x));
              /* PhyML_Printf("\n. Golden section step"); */
            }
          else
            {
              d=p/q;
              u=x+d;
              if (u-a < tol2 || b-u < tol2) d=SIGN(tol1,xm-x);
              /* PhyML_Printf("\n. Parabolic step [e=%f]",e); */
            }
        }
      else
        {
          d=BRENT_CGOLD*(e=(x >= xm ? a-x : b-x));
          /* PhyML_Printf("\n. Golden section step (default) [e=%f tol1=%f a=%f b=%f d=%f x=%f]",e,tol1,a,b,d,x); */
        }

      u=(FABS(d) >= tol1 ? x+d : x+SIGN(tol1,d));
      (*param) = u;
      n_opt_step++;
      old_lnL = fu;
      if(expt == YES) (*param) = log(*param);
      if(logt == YES) (*param) = exp(MIN(1.E+2,*param));
      fu = -(*obj_func)(branch,tree,stree);
      if(logt == YES) (*param) = log(*param);
      if(expt == YES) (*param) = exp(*param);
      /* PhyML_Printf("\n. iter=%d/%d param=%f lnL=%f u: %f x: %f d: %f logt: %d",iter,BRENT_IT_MAX,*param,fu,u,x,d,logt); */

      if(fu <= fx)
        {
          if(u >= x) a=x; else b=x;
          SHFT(v,w,x,u)
          SHFT(fv,fw,fx,fu)
        }
      else
        {
          if (u < x) a=u; else b=u;
          if (fu < fw || FABS(w-x) < SMALL)
            {
              v=w;
              w=u;
              fv=fw;
              fw=fu;
            }
          /* 	  else if (fu <= fv || v == x || v == w) */
          else if (fu < fv || FABS(v-x) < SMALL || FABS(v-w) < SMALL)
            {
              v=u;
              fv=fu;
            }
        }
    }

  /* PhyML_Printf("\n. Too many iterations in Generic_Brent_Lk !"); */
  assert(FALSE);
  return(n_opt_step);
  /* Not Reached ??  *param=x;   */
  /* Not Reached ??  return fx; */
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


/* find ML erstimates of node heights given fixed substitution
   rates on branches. Also optimizes the overall substitution
   rate */
void Round_Optimize_Node_Heights(t_tree *tree)
{
  phydbl cur_lnL, new_lnL;
  int n_iter;


  cur_lnL = UNLIKELY;
  new_lnL = Lk(NULL,tree);


  n_iter = 0;
  while(fabs(new_lnL - cur_lnL) > tree->mod->s_opt->min_diff_lk_local)
    {
      cur_lnL = tree->c_lnL;

      Opt_Node_Heights_Recurr(tree);

      Generic_Brent_Lk(&(tree->rates->clock_r),
                       tree->rates->min_clock,
                       tree->rates->max_clock,
                       tree->mod->s_opt->min_diff_lk_local,
                       tree->mod->s_opt->brent_it_max,
                       tree->mod->s_opt->quickdirty,
                       Wrap_Lk,NULL,tree,NULL,NO,NO);

      printf("\n. cur_lnL=%f new_lnL=%f clock_r=%G root height=%f",
         cur_lnL,new_lnL,tree->rates->clock_r,tree->times->nd_t[tree->n_root->num]);
      new_lnL = tree->c_lnL;
      n_iter++;
      if(n_iter > 100) break;
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Opt_Node_Heights_Recurr(t_tree *tree)
{
  Opt_Node_Heights_Recurr_Pre(tree->n_root,tree->n_root->v[2],tree);
  Opt_Node_Heights_Recurr_Pre(tree->n_root,tree->n_root->v[1],tree);

  Generic_Brent_Lk(&(tree->times->nd_t[tree->n_root->num]),
                   MIN(tree->times->t_prior_max[tree->n_root->num],
                       MIN(tree->times->nd_t[tree->n_root->v[2]->num],
                           tree->times->nd_t[tree->n_root->v[1]->num])),
                   tree->times->t_prior_min[tree->n_root->num],
                   tree->mod->s_opt->min_diff_lk_local,
                   tree->mod->s_opt->brent_it_max,
                   tree->mod->s_opt->quickdirty,
                   Wrap_Lk,NULL,tree,NULL,NO,NO);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Opt_Node_Heights_Recurr_Pre(t_node *a, t_node *d, t_tree *tree)
{
  if(d->tax) return;
  else
    {
      int i;
      phydbl t0,t2,t3;
      phydbl t_min,t_max;
      t_node *v2,*v3;

      v2 = v3 = NULL;
      for(i=0;i<3;i++)
    if((d->v[i] != a) && (d->b[i] != tree->e_root))
      {
        if(!v2) { v2 = d->v[i]; }
        else    { v3 = d->v[i]; }
      }

      Opt_Node_Heights_Recurr_Pre(d,v2,tree);
      Opt_Node_Heights_Recurr_Pre(d,v3,tree);

      t0 = tree->times->nd_t[a->num];
      t2 = tree->times->nd_t[v2->num];
      t3 = tree->times->nd_t[v3->num];

      t_min = t0;
      t_max = MIN(t2,t3);

      t_min = MAX(t_min,tree->times->t_prior_min[d->num]);
      t_max = MIN(t_max,tree->times->t_prior_max[d->num]);

      t_min += tree->rates->min_dt;
      t_max -= tree->rates->min_dt;

      if(t_min > t_max)
        {
          PhyML_Fprintf(stderr,"\n. Err in file %s at line %d\n",__FILE__,__LINE__);
          Exit("\n");
        }

      Generic_Brent_Lk(&(tree->times->nd_t[d->num]),
                       t_min,t_max,
                       tree->mod->s_opt->min_diff_lk_local,
                       tree->mod->s_opt->brent_it_max,
                       tree->mod->s_opt->quickdirty,
                       Wrap_Lk,NULL,tree,NULL,NO,NO);

      /* printf("\n. t%d = %f [%f;%f] lnL = %f",d->num,tree->times->nd_t[d->num],t_min,t_max,tree->c_lnL); */

    }
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Optimize_RR_Params(t_tree *mixt_tree, int verbose)
{
  t_tree *tree;
  t_rmat **r_mat;
  int *permut;
  int n_r_mat;
  int i;
  phydbl lk_new,lk_old;
  phydbl a,b;
  phydbl *opt_val;

  /* PhyML_Printf("\n. Optimizing RR params"); */

  Set_Update_Eigen(YES,mixt_tree->mod);

  n_r_mat = 0;
  tree    = mixt_tree;
  r_mat   = NULL;
  permut  = NULL;
  lk_old  = UNLIKELY;
  lk_new  = UNLIKELY;
  
  do
    {      
      if(tree->is_mixt_tree == YES) tree = tree->next;
      
      for(i=0;i<n_r_mat;i++) if(tree->mod->r_mat == r_mat[i]) break;
      
      if(i == n_r_mat) // tree->mod->r_mat was not found before
        {
          if(!r_mat) r_mat = (t_rmat **)mCalloc(1,sizeof(t_rmat *));
          else       r_mat = (t_rmat **)mRealloc(r_mat,n_r_mat+1,sizeof(t_rmat *));
          r_mat[n_r_mat] = tree->mod->r_mat;
          n_r_mat++;
          
          if((tree->mod->r_mat->optimize == YES) &&
             ((tree->mod->whichmodel == GTR) ||
             ((tree->mod->whichmodel == CUSTOM) &&
              (tree->mod->r_mat->n_diff_rr > 1))))
            {
              int i,iter;
              
              opt_val = (phydbl *)mCalloc(tree->mod->r_mat->n_diff_rr,sizeof(phydbl));
                
              iter = 0;
              do
                {
                  lk_old = mixt_tree->c_lnL;
                  int failed = NO;
                  
                  if(tree->mod->r_mat->n_diff_rr > 2)
                    {
                      for(i=0;i<tree->mod->r_mat->n_diff_rr;i++) opt_val[i] = tree->mod->r_mat->rr_val->v[i];

                      BFGS(mixt_tree,tree->mod->r_mat->rr_val->v,tree->mod->r_mat->n_diff_rr,1.e-5,tree->mod->s_opt->min_diff_lk_local,1.e-3,NO,NO,NO,
                           &Return_Abs_Lk,
                           &Num_Derivative_Several_Param,
                           &Lnsrch,&failed);
                      
                      if(failed == YES)  for(i=0;i<tree->mod->r_mat->n_diff_rr;i++) tree->mod->r_mat->rr_val->v[i] = opt_val[i];
                    }

                  permut = Permutate(tree->mod->r_mat->n_diff_rr);
                  
                  for(i=0;i<tree->mod->r_mat->n_diff_rr;i++) opt_val[i] = tree->mod->r_mat->rr_val->v[i];

                  for(i=0;i<tree->mod->r_mat->n_diff_rr;i++)
                    {
                      // Remember rr_val = log(rr) hence the upper and lower bounds below
                      a = tree->mod->r_mat->rr_val->v[permut[i]]-10.;
                      b = tree->mod->r_mat->rr_val->v[permut[i]]+10.;

                      Generic_Brent_Lk(&(tree->mod->r_mat->rr_val->v[permut[i]]),
                                       a,b,
                                       tree->mod->s_opt->min_diff_lk_local,
                                       tree->mod->s_opt->brent_it_max,
                                       tree->mod->s_opt->quickdirty,
                                       Wrap_Lk,NULL,mixt_tree,NULL,NO,NO);                        
                    }

                  
                  if(mixt_tree->c_lnL < lk_old)
                    {
                      for(i=0;i<tree->mod->r_mat->n_diff_rr;i++) tree->mod->r_mat->rr_val->v[i] = opt_val[i];
                      Lk(NULL,mixt_tree);
                    }
                  
                  
                  if(verbose) Print_Lk(tree->mixt_tree?
                                       tree->mixt_tree:
                                       tree,"[GTR parameters     ]");
                  
                  lk_new = mixt_tree->c_lnL;
                  
                  Free(permut);

                  
                  if(lk_new < lk_old - tree->mod->s_opt->min_diff_lk_global)
                    {
                      PhyML_Printf("\n. lk_new: %f lk_old: %f",lk_new,lk_old);
                      assert(FALSE);
                    }
                  
                  if(fabs(lk_new-lk_old) < tree->mod->s_opt->min_diff_lk_local) break;
                }
              /* while(++iter < tree->mod->s_opt->brent_it_max); */
              while(++iter < 1);

              Free(opt_val);
              
              if(iter == tree->mod->s_opt->brent_it_max)
                {
                  if(tree->verbose > VL0)
                    {
                      PhyML_Printf("\n. Failed to optimize GTR parameters this round...");
                    }
                }              
            }
        }
      
      tree = tree->next;
      if(!tree) break;
      
  }
  while(1);

  if(r_mat) Free(r_mat);
  
  
  Set_Update_Eigen(NO,mixt_tree->mod);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Optimize_TsTv(t_tree *mixt_tree, int verbose)
{
  scalar_dbl **tstv;
  int n_tstv;
  t_tree *tree;
  int i;

  if (mixt_tree->mod->io->datatype != NT) return;

  /* PhyML_Printf("\n. Optimizing Tstv"); */

  Set_Update_Eigen(YES,mixt_tree->mod);

  tstv   = NULL;
  n_tstv = 0;
  tree   = mixt_tree;

  do
    {
      if(tree->is_mixt_tree == YES) tree = tree->next;
      
      for(i=0;i<n_tstv;i++) if(tree->mod->kappa == tstv[i]) break;
      
      if(i == n_tstv)
        {
          if(!tstv) tstv = (scalar_dbl **)mCalloc(1,sizeof(scalar_dbl *));
          else      tstv = (scalar_dbl **)mRealloc(tstv,n_tstv+1,sizeof(scalar_dbl *));
          tstv[n_tstv] = tree->mod->kappa;
          n_tstv++;
          
          if(tree->mod->kappa->optimize == YES)
            {
              phydbl a,c;
              
              /* a = tree->mod->kappa->v * .1; */
              /* c = tree->mod->kappa->v * 10.; */
              a = TSTV_MIN;
              c = TSTV_MAX;
                            
              Generic_Brent_Lk(&(tree->mod->kappa->v),
                               a,c,
                               tree->mod->s_opt->min_diff_lk_local,
                               tree->mod->s_opt->brent_it_max,
                               tree->mod->s_opt->quickdirty,
                               Wrap_Lk,NULL,mixt_tree,NULL,NO,NO);
              
              if(verbose)
                {
                  Print_Lk(mixt_tree,"[Ts/ts ratio        ]");
                  PhyML_Printf("[%10f]",tree->mod->kappa->v);
                }              
            }
        }
      
      tree = tree->next;
      
    }
  while(tree);
  
  if(tstv) Free(tstv);
  
  Set_Update_Eigen(NO,mixt_tree->mod);
}



//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Optimize_Lvar(t_tree *mixt_tree, int verbose)
{
  t_tree *tree;
  scalar_dbl **l_var;
  int n_l_var;
  int i;
  
  tree = mixt_tree;
  n_l_var = 0;
  l_var = NULL;
  do
    {
      if(tree->is_mixt_tree == YES) tree = tree->next;
      
      for(i=0;i<n_l_var;i++) if(tree->mod->l_var_sigma == l_var[i]) break;

      if(i == n_l_var)
        {
          if(!l_var) l_var = (scalar_dbl **)mCalloc(1,sizeof(scalar_dbl *));
          else       l_var = (scalar_dbl **)mRealloc(l_var,n_l_var+1,sizeof(scalar_dbl *));
          l_var[n_l_var] = tree->mod->l_var_sigma;
          n_l_var++;

          if(tree->mod->gamma_mgf_bl == YES && tree->mod->s_opt->opt_gamma_br_len == YES)
            {
              Generic_Brent_Lk(&(tree->mod->l_var_sigma->v),
                               tree->mod->l_var_min,
                               tree->mod->l_var_max,
                               tree->mod->s_opt->min_diff_lk_local,
                               tree->mod->s_opt->brent_it_max,
                               tree->mod->s_opt->quickdirty,
                               Wrap_Lk,NULL,mixt_tree,NULL,NO,NO);
              

              if(verbose)
                {
                  Print_Lk(tree,"[Branch len. var.   ]");
                  PhyML_Printf("[%10f]",tree->mod->l_var_sigma->v);
                }
            }
        }
      tree = tree->next;
    }
  while(tree);
  
  if(l_var) Free(l_var);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Optimize_Pinv(t_tree *mixt_tree, int verbose)
{
  scalar_dbl **pinv;
  int n_pinv;
  t_tree *tree;
  int i;

  /* PhyML_Printf("\n. Optimizing Pinv"); */

  Set_Update_Eigen(NO,mixt_tree->mod);

  pinv   = NULL;
  n_pinv = 0;
  tree   = mixt_tree;

  do
    {
      for(i=0;i<n_pinv;i++) if(tree->mod->ras->pinvar == pinv[i]) break;

      if(i == n_pinv)
        {
          if(!pinv) pinv = (scalar_dbl **)mCalloc(1,sizeof(scalar_dbl *));
          else      pinv = (scalar_dbl **)mRealloc(pinv,n_pinv+1,sizeof(scalar_dbl *));
          pinv[n_pinv] = tree->mod->ras->pinvar;
          n_pinv++;
           
          if(tree->mod->ras->pinvar->optimize == YES && (tree->mod->ras->alpha->optimize == NO || tree->mod->ras->n_catg == 1))
            {
              Generic_Brent_Lk(&(tree->mod->ras->pinvar->v),
                               PINV_MIN,PINV_MAX,
                               tree->mod->s_opt->min_diff_lk_local,
                               tree->mod->s_opt->brent_it_max,
                               tree->mod->s_opt->quickdirty,
                               Wrap_Lk,NULL,mixt_tree,NULL,NO,NO);
              
              if(verbose)
                {
                  Print_Lk(tree,"[P-inv              ]");
                  PhyML_Printf("[%10f]",tree->mod->ras->pinvar->v);
                }
            }
        }
      
      tree = tree->next_mixt;

    }
  while(tree);

  if(pinv) Free(pinv);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Optimize_Alpha(t_tree *mixt_tree, int verbose)
{
  scalar_dbl **alpha;
  int n_alpha;
  t_tree *tree;
  int i;

  /* PhyML_Printf("\n. Optimizing Alpha"); */

  Set_Update_Eigen(NO,mixt_tree->mod);

  alpha   = NULL;
  n_alpha = 0;
  tree   = mixt_tree;

  do
    {

      if(tree->mod->ras->alpha->optimize == YES && tree->mod->ras->n_catg > 1)
        {
          for(i=0;i<n_alpha;i++) if(tree->mod->ras->alpha == alpha[i]) break;

          if(i == n_alpha)
            {
              if(!alpha) alpha = (scalar_dbl **)mCalloc(1,sizeof(scalar_dbl *));
              else       alpha = (scalar_dbl **)mRealloc(alpha,n_alpha+1,sizeof(scalar_dbl *));
              alpha[n_alpha] = tree->mod->ras->alpha;
              n_alpha++;

              if(tree->mod->ras->alpha->optimize == YES &&
                 tree->mod->ras->free_mixt_rates == NO &&
                 tree->mod->ras->pinvar->optimize == NO)
                {
                  
                  if(tree->mod->ras->n_catg > 1)
                    {
                      Generic_Brent_Lk(&(tree->mod->ras->alpha->v),
                                       ALPHA_MIN,ALPHA_MAX,
                                       tree->mod->s_opt->min_diff_lk_local,
                                       tree->mod->s_opt->brent_it_max,
                                       tree->mod->s_opt->quickdirty,
                                       Wrap_Lk,NULL,mixt_tree,NULL,NO,NO);
                    }
                  if(verbose)
                    {
                      Print_Lk(mixt_tree,"[Alpha              ]");
                      PhyML_Printf("[%10f]",tree->mod->ras->alpha->v);
                    }
                }
            }
        }
      tree = tree->next_mixt;
    }
  while(tree);

  if(alpha) Free(alpha);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Optimize_Free_Rate(t_tree *mixt_tree, int verbose)
{
  t_tree *tree;
  int fast;
  int i,pos,failed;
  int *permut;
  phydbl lk_before, lk_after;
  tree = mixt_tree;

  lk_before = lk_after = UNLIKELY;

  do
    {
      if((tree->mod->s_opt->opt_free_mixt_rates == YES) && (tree->mod->ras->free_mixt_rates == YES) && (tree->mod->ras->n_catg > 1))
        {
          if(tree->mod->s_opt->serial_free_rates == YES)
            {
              fast = NO;
              lk_before = mixt_tree->c_lnL;
              Optimize_Free_Rate_Weights(tree,fast,verbose);
              lk_after = mixt_tree->c_lnL;
              Optimize_Free_Rate_Rr(tree,fast,verbose);
              lk_after = mixt_tree->c_lnL;


              if (lk_after < lk_before - tree->mod->s_opt->min_diff_lk_global)
              {
                PhyML_Fprintf(stderr, "\n. lk_before: %f lk_after: %f diff: %G",
                              lk_before, lk_after, lk_before - lk_after);
                PhyML_Fprintf(stderr, "\n. Err. in file %s at line %d\n",
                              __FILE__, __LINE__);
                Exit("");
              }
            }        
          else
            {
                fast = YES;
                lk_before = mixt_tree->c_lnL;
                Optimize_Free_Rate_Weights(tree,fast,verbose);
                lk_after = mixt_tree->c_lnL;
                Optimize_Free_Rate_Rr(tree,fast,verbose);
                lk_after = mixt_tree->c_lnL;


                if (lk_after < lk_before - tree->mod->s_opt->min_diff_lk_global)
                {
                  PhyML_Fprintf(stderr,"\n. lk_before: %f lk_after: %f diff: %G",lk_before,lk_after,lk_before-lk_after);
                  PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d\n",__FILE__,__LINE__);
                  Exit("");
                }
            }


          if(FABS(lk_before - lk_after) > 0.001)
            {
              phydbl **x,*cpy;
              
              x = (phydbl **)mCalloc(2*tree->mod->ras->n_catg,sizeof(phydbl *));
              cpy = (phydbl *)mCalloc(2*tree->mod->ras->n_catg,sizeof(phydbl));

              
              lk_before = mixt_tree->c_lnL;

              pos = 0;
              for(i=0;i<tree->mod->ras->n_catg;i++) x[pos++] = tree->mod->ras->gamma_rr_unscaled->v+i;
              for(i=0;i<tree->mod->ras->n_catg;i++) x[pos++] = tree->mod->ras->gamma_r_proba_unscaled->v+i;

              pos = 0;
              for(i=0;i<tree->mod->ras->n_catg;i++) cpy[pos++] = tree->mod->ras->gamma_rr_unscaled->v[i];
              for(i=0;i<tree->mod->ras->n_catg;i++) cpy[pos++] = tree->mod->ras->gamma_r_proba_unscaled->v[i];

              
              failed = YES;
              BFGS_Nonaligned(tree,x,2*tree->mod->ras->n_catg,1.e-5,tree->mod->s_opt->min_diff_lk_global,1.e-5,NO,NO,NO,
                              &Return_Abs_Lk,
                              &Num_Derivative_Several_Param_Nonaligned,
                              &Lnsrch_Nonaligned,&failed);

              if(failed == YES)
                {
                  pos = 0;
                  for(i=0;i<tree->mod->ras->n_catg;i++) tree->mod->ras->gamma_rr_unscaled->v[i] = cpy[pos++];
                  for(i=0;i<tree->mod->ras->n_catg;i++) tree->mod->ras->gamma_r_proba_unscaled->v[i] = cpy[pos++];

                  permut = Permutate(tree->mod->ras->n_catg);
                  
                  for(i=0;i<tree->mod->ras->n_catg;++i)
                    {
                      phydbl a,c;
                      
                      a = tree->mod->ras->gamma_rr_unscaled->v[permut[i]] - 2.0;
                      c = tree->mod->ras->gamma_rr_unscaled->v[permut[i]] + 2.0;
                      
                      Generic_Brent_Lk(&(tree->mod->ras->gamma_rr_unscaled->v[permut[i]]),
                                       a,c,
                                       tree->mod->s_opt->min_diff_lk_local,
                                       tree->mod->s_opt->brent_it_max,
                                       tree->mod->s_opt->quickdirty,
                                       Wrap_Lk,NULL,mixt_tree,NULL,NO,NO);      
                    }

                  for(i=0;i<tree->mod->ras->n_catg;++i)
                    {
                      phydbl a,c;
                      
                      a = tree->mod->ras->gamma_r_proba_unscaled->v[permut[i]] - 2.0;
                      c = tree->mod->ras->gamma_r_proba_unscaled->v[permut[i]] + 2.0;
                      
                      Generic_Brent_Lk(&(tree->mod->ras->gamma_r_proba_unscaled->v[permut[i]]),
                                       a,c,
                                       tree->mod->s_opt->min_diff_lk_local,
                                       tree->mod->s_opt->brent_it_max,
                                       tree->mod->s_opt->quickdirty,
                                       Wrap_Lk,NULL,mixt_tree,NULL,NO,NO);      
                    }                  

                  Free(permut);
                                
                }
 
              lk_after = mixt_tree->c_lnL;

              if(lk_after < lk_before - tree->mod->s_opt->min_diff_lk_global)
                {
                  PhyML_Fprintf(stderr,"\n. lk_before: %f lk_after: %f diff: %G",lk_before,lk_after,lk_before-lk_after);
                  PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d\n",__FILE__,__LINE__);
                  Exit("");
                }

              Free(x);
              Free(cpy);
            }
        }
      tree = tree->next_mixt;
    }
  while(tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Optimize_Free_Rate_Rr(t_tree *tree, int fast, int verbose)
{
  phydbl lk_before, lk_after;

  lk_before = tree->c_lnL;

  if(tree->prev == NULL && tree->next == NULL)
    {
      int i;

      for(i=0;i<tree->mod->ras->n_catg;i++)
        {
          phydbl a,c;
          
          // a = tree->mod->ras->gamma_rr_unscaled->v[i] - 2.0;
          // c = tree->mod->ras->gamma_rr_unscaled->v[i] + 2.0;
          a = - 2.0;
          c = + 2.0;
                    
          Generic_Brent_Lk(&(tree->mod->ras->gamma_rr_unscaled->v[i]),
                           a,c,
                           tree->mod->s_opt->min_diff_lk_local,
                           tree->mod->s_opt->brent_it_max,
                           tree->mod->s_opt->quickdirty,
                           Wrap_Lk,NULL,tree,NULL,NO,NO);
          

          lk_after = tree->c_lnL;
          
          if(lk_after < lk_before - tree->mod->s_opt->min_diff_lk_global)
            {
              PhyML_Fprintf(stderr,"\n. lk_before: %f lk_after: %f diff: %G",lk_before,lk_after,lk_before-lk_after);
              PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d\n",__FILE__,__LINE__);
              Exit("");
            }
        }
    }
  else
    {
      int i;
      for(i=0;i<tree->mod->ras->n_catg;i++)
        {
          phydbl a,c;
          
          // a = tree->mod->ras->gamma_rr_unscaled->v[i] - 2.0;
          // c = tree->mod->ras->gamma_rr_unscaled->v[i] + 2.0;
          a = - 2.0;
          c = + 2.0;
                    
          Generic_Brent_Lk(&(tree->mod->ras->gamma_rr_unscaled->v[i]),
                           a,c,
                           tree->mod->s_opt->min_diff_lk_local,
                           tree->mod->s_opt->brent_it_max,
                           tree->mod->s_opt->quickdirty,
                           Wrap_Lk,NULL,tree,NULL,NO,NO);
        }
    }

  lk_after = tree->c_lnL;

  if(lk_after < lk_before - tree->mod->s_opt->min_diff_lk_global)
    {
      PhyML_Fprintf(stderr,"\n. lk_before: %f lk_after: %f diff: %G",lk_before,lk_after,lk_before-lk_after);
      PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d\n",__FILE__,__LINE__);
      Exit("");
    }

  if(verbose) Print_Lk(tree,"[Rate class values  ]");
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Optimize_Free_Rate_Weights(t_tree *tree, int fast, int verbose)
{
  int i,j;
  phydbl wm, lk_before, lk_after;

  lk_before = tree->c_lnL;

  if(fast == YES)
    {
      tree->mod->ras->normalise_rr = NO;
      for(i=0;i<tree->mod->ras->n_catg;i++) tree->mod->ras->gamma_rr_unscaled->v[i] = log(tree->mod->ras->gamma_rr->v[i]);
    }
  
  for(i=0;i<tree->mod->ras->n_catg;i++)
    {
      phydbl a,c;
      
      // a = tree->mod->ras->gamma_r_proba_unscaled->v[i] - 2.0;
      // c = tree->mod->ras->gamma_r_proba_unscaled->v[i] + 2.0;
      a = - 2.0;
      c = + 2.0;

      Generic_Brent_Lk(&(tree->mod->ras->gamma_r_proba_unscaled->v[i]),
                       a,c,
                       tree->mod->s_opt->min_diff_lk_local,
                       tree->mod->s_opt->brent_it_max,
                       tree->mod->s_opt->quickdirty,
                       Wrap_Lk,NULL,tree,NULL,NO,NO);      


      if(fast == YES)
        {
          wm = Weighted_Mean(tree->mod->ras->gamma_rr->v,
                             tree->mod->ras->gamma_r_proba->v,
                             tree->mod->ras->n_catg);
          
          for(j=0;j<2*tree->n_otu-1;++j) tree->a_edges[j]->l->v *= wm;
          for(j=0;j<tree->mod->ras->n_catg;j++) tree->mod->ras->gamma_rr->v[j] /= wm; 
          for(j=0;j<tree->mod->ras->n_catg;j++) tree->mod->ras->gamma_rr_unscaled->v[j] = log(tree->mod->ras->gamma_rr->v[j]);  
        }
    }

  tree->mod->ras->normalise_rr = YES;

  tree->mod->ras->gamma_rr_unscaled->v[0] = 0.0;
  for(j=1;j<tree->mod->ras->n_catg;j++) tree->mod->ras->gamma_rr_unscaled->v[j] = log(tree->mod->ras->gamma_rr->v[j]/tree->mod->ras->gamma_rr->v[0]);
          
  lk_after = tree->c_lnL;

  if(lk_after < lk_before - tree->mod->s_opt->min_diff_lk_global)
    {
      PhyML_Fprintf(stderr,"\n. lk_before: %f lk_after: %f diff: %G",lk_before,lk_after,lk_before-lk_after);
      PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d\n",__FILE__,__LINE__);
      Exit("");
    }

  if(verbose) Print_Lk(tree,"[Rate class freqs.  ]");
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Optimize_State_Freqs(t_tree *mixt_tree, int verbose)
{
  /* vect_dbl **freqs; */
  /* int n_freqs; */
  t_tree *tree;
  int i;
  phydbl lk_new,lk_old;
  int *permut;
  phydbl *opt_val;

  /* PhyML_Printf("\n. Optimizing State Freqs"); */
  
  Set_Update_Eigen(YES,mixt_tree->mod);

  /* freqs   = NULL; */
  /* n_freqs = 0; */
  tree    = mixt_tree;
  lk_old  = UNLIKELY;
  lk_new  = UNLIKELY;
  opt_val = NULL;

  do
    {
      if(tree == mixt_tree && tree->next != NULL) tree = tree->next;

      if(tree->mod->e_frq->type == ML)
        {
          int iter;
          
          opt_val = (phydbl *)mCalloc(tree->mod->ns,sizeof(phydbl));
          
          iter = 0;
          do
            {
              lk_old = mixt_tree->c_lnL;
              int failed = NO;
                  
              for(i=0;i<tree->mod->ns;i++) opt_val[i] = tree->mod->e_frq->pi_unscaled->v[i];
              
              /* PhyML_Printf("\n. -- BFGS lk=%f",mixt_tree->c_lnL); */
              BFGS(mixt_tree,tree->mod->e_frq->pi_unscaled->v,tree->mod->ns,1.e-5,tree->mod->s_opt->min_diff_lk_local,1.e-2,NO,NO,NO,
                   &Return_Abs_Lk,
                   &Num_Derivative_Several_Param,
                   &Lnsrch,&failed);
              /* PhyML_Printf("\n. ++ BFGS lk=%f",mixt_tree->c_lnL); */
              /* for(i=0;i<tree->mod->ns;i++) PhyML_Printf("\n. %p %f %f",tree,tree->mod->e_frq->pi_unscaled->v[i],tree->mod->e_frq->pi->v[i]); */
              
              if(failed == YES) for(i=0;i<tree->mod->ns;i++) tree->mod->e_frq->pi_unscaled->v[i] = opt_val[i];

                           
              permut = Permutate(tree->mod->ns);
              
              for(i=0;i<tree->mod->ns;++i)
                {
                  phydbl a,c;
                  
                  a = tree->mod->e_frq->pi_unscaled->v[permut[i]] / 2.;
                  c = tree->mod->e_frq->pi_unscaled->v[permut[i]] * 2.+1;
                  
                  /* PhyML_Printf("\n. -- lk=%f",mixt_tree->c_lnL); */
                  Generic_Brent_Lk(&(tree->mod->e_frq->pi_unscaled->v[permut[i]]),
                                   /* UNSCALED_E_FRQ_MIN,UNSCALED_E_FRQ_MAX, */
                                   a,c,
                                   tree->mod->s_opt->min_diff_lk_local,
                                   tree->mod->s_opt->brent_it_max,
                                   tree->mod->s_opt->quickdirty,
                                   Wrap_Lk,NULL,mixt_tree,NULL,NO,NO);      
                  /* PhyML_Printf("\n. ++ lk=%f",mixt_tree->c_lnL); */
                }

              Free(permut);
              
              if(verbose)
                {
                  Print_Lk(mixt_tree,"[Character freqs.   ]");
                }
              
              lk_new = mixt_tree->c_lnL;
              
              /* PhyML_Printf("\n. lk_new=%f lk_old=%f",lk_new,lk_old); */
              assert(lk_new > lk_old - tree->mod->s_opt->min_diff_lk_local);
              if(fabs(lk_new-lk_old) < tree->mod->s_opt->min_diff_lk_local) break;
            }
          /* while(++iter < tree->mod->s_opt->brent_it_max); */
          while(++iter < 1);
          
          Free(opt_val);
          
          if(iter == tree->mod->s_opt->brent_it_max)
            {
              if(tree->verbose > VL0)
                {
                  PhyML_Printf("\n. Failed to optimize frequency parameters this round...");
                }
            }
        }
    
      
      tree = tree->next;
      
    }
  while(tree);
  
  
  Set_Update_Eigen(NO,mixt_tree->mod);
  
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Optimize_Rmat_Weights(t_tree *mixt_tree, int verbose)
{
  scalar_dbl *r_mat_weight;

  /* PhyML_Printf("\n. Optimizing Rmat Weights"); */

  Set_Update_Eigen(NO,mixt_tree->mod);

  if(mixt_tree->is_mixt_tree == NO) return;

  r_mat_weight = mixt_tree->next->mod->r_mat_weight;

  if (mixt_tree->next->mod->s_opt->opt_rmat_weight == YES)
  {
    do
    {
      phydbl a, c;

      a = r_mat_weight->v * .1;
      c = r_mat_weight->v * 10.;

      Generic_Brent_Lk(&(r_mat_weight->v), a, c,
                       mixt_tree->mod->s_opt->min_diff_lk_local,
                       mixt_tree->mod->s_opt->brent_it_max,
                       mixt_tree->mod->s_opt->quickdirty, Wrap_Lk, NULL,
                       mixt_tree, NULL, NO, NO);

      if (verbose) Print_Lk(mixt_tree, "[Rate mat. weights  ]");

      r_mat_weight = r_mat_weight->next;
    } while (r_mat_weight);
  }

  Set_Update_Eigen(NO,mixt_tree->mod);
  }

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Optimize_Efrq_Weights(t_tree *mixt_tree, int verbose)
{
  scalar_dbl *e_frq_weight;

  /* PhyML_Printf("\n. Optimizing Efrq Weights"); */

  Set_Update_Eigen(NO,mixt_tree->mod);

  if(mixt_tree->is_mixt_tree == NO) return;

  e_frq_weight = mixt_tree->next->mod->e_frq_weight;


  if(mixt_tree->next->mod->s_opt->opt_efrq_weight == YES)
    {
      do
        {
          phydbl a,c;
          
          a = e_frq_weight->v * .1;
          c = e_frq_weight->v * 10.+1;
          
          Generic_Brent_Lk(&(e_frq_weight->v),
                           a,c,
                           mixt_tree->mod->s_opt->min_diff_lk_local,
                           mixt_tree->mod->s_opt->brent_it_max,
                           mixt_tree->mod->s_opt->quickdirty,
                           Wrap_Lk,NULL,mixt_tree,NULL,NO,NO);      
          
          if(verbose)
            {
              Print_Lk(mixt_tree,"[Equ. frq. weights  ]");
            }
          
          e_frq_weight = e_frq_weight->next;
        }
      while(e_frq_weight);
    }

  Set_Update_Eigen(NO,mixt_tree->mod);

}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Optimize_Lambda(t_tree *mixt_tree, int verbose)
{
  scalar_dbl **lambda;
  int n_lambda;
  t_tree *tree;
  int i;

  if (mixt_tree->mod->io->datatype != NT) return;

  /* PhyML_Printf("\n. Optimizing Lambda"); */

  Set_Update_Eigen(YES,mixt_tree->mod);

  lambda   = NULL;
  n_lambda = 0;
  tree     = mixt_tree;
  
  do
    {
      if(tree->next) tree = tree->next;
      
      for(i=0;i<n_lambda;i++) if(tree->mod->lambda == lambda[i]) break;
      
      if(i == n_lambda)
        {
          if(!lambda) lambda = (scalar_dbl **)mCalloc(1,sizeof(scalar_dbl *));
          else        lambda = (scalar_dbl **)mRealloc(lambda,n_lambda+1,sizeof(scalar_dbl *));
          lambda[n_lambda] = tree->mod->lambda;
          n_lambda++;
          
          if(tree->mod->lambda->optimize == YES)
            {
              /* PhyML_Printf("\n-- Lambda: %f",tree->c_lnL); */
              Generic_Brent_Lk(&(tree->mod->lambda->v),
                               0.001,100.,
                               tree->mod->s_opt->min_diff_lk_local,
                               tree->mod->s_opt->brent_it_max,
                               tree->mod->s_opt->quickdirty,
                               Wrap_Lk,NULL,mixt_tree,NULL,NO,NO);
              /* PhyML_Printf("\n++ Lambda: %f",tree->c_lnL); */
              
              if(verbose)
                {
                  Print_Lk(mixt_tree,"[Lambda             ]");
                  PhyML_Printf("[%10f]",tree->mod->lambda->v);
                }
            }
        }
      
      tree = tree->next;
      
    }
  while(tree);
  
  if(lambda) Free(lambda);
  
  Set_Update_Eigen(NO,mixt_tree->mod);
  
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

#if (defined PHYTIME)

void Least_Square_Node_Ages(t_tree *tree)
{

  phydbl new_error,cur_error,sum_error;
  int i,j;
  phydbl young,old;
  int dir1,dir2;
  t_node *n;
  
  Update_Ancestors(tree->n_root,tree->n_root->v[2],tree->n_root->b[2],tree);
  Update_Ancestors(tree->n_root,tree->n_root->v[1],tree->n_root->b[1],tree);  

  TIMES_Randomize_Node_Ages(tree);
  
  assert(fabs(tree->times->nd_t[tree->n_root->num]) > SMALL);

  cur_error = BIG;
  new_error = BIG;
  do
    {
      cur_error = new_error;
      
      sum_error = 0.0;
      for(i=0;i<2*tree->n_otu-1;++i)
        {
          if(tree->a_nodes[i]->tax == NO)
            {
              n = tree->a_nodes[i];
              old = -1.;
              
              dir1 = dir2 = -1;
              for(j=0;j<3;++j)
                {
                  if(n->v[j] != n->anc && n->b[j] != tree->e_root)
                    {
                      if(dir1 < 0) dir1 = j;
                      else         dir2 = j;
                    }
                  else
                    {
                      if(n != tree->n_root) old = tree->times->nd_t[n->anc->num];
                      else old = 2.*tree->times->nd_t[tree->n_root->num];
                    }
                }
              
              young = MIN(tree->times->nd_t[n->v[dir1]->num],
                          tree->times->nd_t[n->v[dir2]->num]);
              
              sum_error += Generic_Brent(tree->times->nd_t + i,
                                         young,old,1.E-10,10000,
                                         TIMES_Least_Square_Criterion,
                                         tree);
              
              /* PhyML_Printf("\n. Node %3d%c time: %15f err: %15f young: %15f old: %15f %15f %15f", */
              /*              i, */
              /*              (n == tree->n_root)?'*':' ', */
              /*              tree->times->nd_t[i], */
              /*              sum_error, */
              /*              young,old, */
              /*              tree->times->nd_t[n->v[dir1]->num], */
              /*              tree->times->nd_t[n->v[dir2]->num]); */
            }

          
#if (defined PHYREX || PHYRIME)
          if(RATES_Check_Node_Times(tree)) Exit("\n");
#endif

        }

      sum_error += Generic_Brent(&(tree->rates->clock_r),
                                 tree->rates->min_clock,
                                 tree->rates->max_clock,
                                 1.E-10,10000,
                                 TIMES_Least_Square_Criterion,
                                 tree);

      /* PhyML_Printf("\n clock_r: %g",tree->rates->clock_r); */
      
      new_error = sum_error;
      /* assert(new_error < cur_error+SMALL); */
    }
  while(fabs(new_error-cur_error) > 1.E-10);

}
#endif

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

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
