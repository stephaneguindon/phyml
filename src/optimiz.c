/*

PhyML :  a program that  computes maximum likelihood  phyLOGenies from
DNA or AA homoLOGous sequences

Copyright (C) Stephane Guindon. Oct 2003 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/

#include "optimiz.h"


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
                   NO);


  if(tree->c_lnL < lk_init - tree->mod->s_opt->min_diff_lk_global)
    {
      PhyML_Printf("\n== %.10f < %.10f --> diff=%.10f param value = %f\n",
         tree->c_lnL,lk_init,
         tree->c_lnL-lk_init,
         *param);
      Exit("\n== Optimisation failed !\n");
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


phydbl Generic_Brent(phydbl ax, phydbl bx, phydbl cx, phydbl tol,
             phydbl *xmin, t_tree *tree, int n_iter_max,
             int quickdirty)
{
  int iter;
  phydbl a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  phydbl e=0.0;
  phydbl init_lnL;


  d=0.0;
  a=((ax < cx) ? ax : cx);
  b=((ax > cx) ? ax : cx);
  x=w=v=bx;
  (*xmin) = bx;
  fw=fv=fx=-Lk(NULL,tree);
  init_lnL = -fw;

  PhyML_Printf("\n. init_lnL = %f a=%f b=%f c=%f\n",init_lnL,ax,bx,cx);

  for(iter=1;iter<=n_iter_max;iter++)
    {
      xm=0.5*(a+b);
      tol2=2.0*(tol1=tol*FABS(x)+BRENT_ZEPS);


      if(FABS(x - xm) <= (tol2 - 0.5 * (b - a)))
    {
      *xmin = x;
      Lk(NULL,tree);
      if(tree->c_lnL < init_lnL - tree->mod->s_opt->min_diff_lk_local)
        {
          PhyML_Printf("\n== Err. in file %s at line %d\n",__FILE__,__LINE__);
          Warn_And_Exit("");
        }
      return tree->c_lnL;
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
/* 	      PhyML_Printf("Golden section step\n"); */
        }
      else
        {
          d=p/q;
          u=x+d;
          if (u-a < tol2 || b-u < tol2) d=SIGN(tol1,xm-x);
/* 	      PhyML_Printf("Parabolic step [e=%f]\n",e); */
        }
        }
      else
    {
      d=BRENT_CGOLD*(e=(x >= xm ? a-x : b-x));
/* 	  PhyML_Printf("Golden section step (default) [e=%f tol1=%f a=%f b=%f d=%f]\n",e,tol1,a,b,d); */
    }

      u=(FABS(d) >= tol1 ? x+d : x+SIGN(tol1,d));
      (*xmin) = FABS(u);
      fu = -Lk(NULL,tree);

      PhyML_Printf("\n. iter=%d/%d param=%f loglk=%f",iter,BRENT_IT_MAX,*xmin,tree->c_lnL);

/*       if(fu <= fx) */
      if(fu < fx)
    {
/* 	  if(u >= x) a=x; else b=x; */
      if(u > x) a=x; else b=x;
      SHFT(v,w,x,u)
      SHFT(fv,fw,fx,fu)
    }
      else
    {
      if (u < x) a=u; else b=u;
/* 	  if (fu <= fw || w == x) */
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

  Exit("\n. Too many iterations in Generic_Brent !");
  return(-1);
  /* Not Reached ??  *xmin=x;   */
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

   b_fcus->l->v=x1;
   f1 = -Lk(b_fcus,tree);
   b_fcus->l->v=x2;
   f2 = -Lk(b_fcus,tree);
   while (FABS(x3-x0) > tol*(FABS(x1)+FABS(x2)))
     {
       if (f2 < f1)
     {
       SHFT3(x0,x1,x2,GOLDEN_R*x1+GOLDEN_C*x3)
       b_fcus->l->v=x2;
       SHFT2(f1,f2,-Lk(b_fcus,tree))
     }
       else
     {
       SHFT3(x3,x2,x1,GOLDEN_R*x2+GOLDEN_C*x0)
       b_fcus->l->v=x1;
       SHFT2(f2,f1,-Lk(b_fcus,tree))
     }
     }
   if (f1 < f2)
     {
       *xmin=FABS(x1);
       return -f1;
     }
   else
     {
       *xmin=FABS(x2);
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
       PhyML_Printf("fb=%f fc=%f\n",*fb,*fc);
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

phydbl Br_Len_Brent(phydbl prop_min, phydbl prop_max, t_edge *b_fcus, t_tree *tree)
{
  t_tree *loc_tree;
  t_edge *loc_b;
  phydbl lk_begin, lk_end;

  lk_begin = UNLIKELY;
  lk_end   = UNLIKELY;
  loc_tree = tree;
  loc_b    = b_fcus;


  /*! Rewind back to the first mixt_tree */
  while(loc_tree->prev)
    {
      loc_tree = loc_tree->prev;
      loc_b    = loc_b->prev;
    }

  if(tree->is_mixt_tree)
    {
      MIXT_Br_Len_Brent(prop_min,prop_max,b_fcus,tree);
      return loc_tree->c_lnL;
    }

  if(b_fcus->l->onoff == OFF) return loc_tree->c_lnL;

  if(isinf(prop_min) || isnan(prop_min)) prop_min = 1.E-04;
  if(isinf(prop_max) || isnan(prop_max)) prop_max = 1.E+04;

  lk_begin = Lk(loc_b,loc_tree); /*! We can't assume that the log-lk value is up-to-date */

  Generic_Brent_Lk(&(b_fcus->l->v),
                   MAX(tree->mod->l_min,MIN(tree->mod->l_max,b_fcus->l->v))*MAX(1.E-04,prop_min),
                   MAX(tree->mod->l_min,MIN(tree->mod->l_max,b_fcus->l->v))*MIN(1.E+04,prop_max),
                   tree->mod->s_opt->min_diff_lk_local,
                   tree->mod->s_opt->brent_it_max,
                   tree->mod->s_opt->quickdirty,
                   Wrap_Lk_At_Given_Edge,
                   loc_b,loc_tree,NULL,NO);

  /* if(tree->mod->gamma_mgf_bl == YES) */
  /*   { */
  /*     if(b_fcus->num == 0) */
  /*       { */
  /*         Generic_Brent_Lk(&(b_fcus->l_var), */
  /*                          1.E-4,100., */
  /*                          tree->mod->s_opt->min_diff_lk_local, */
  /*                          tree->mod->s_opt->brent_it_max, */
  /*                          tree->mod->s_opt->quickdirty, */
  /*                          Wrap_Lk_At_Given_Edge,loc_b,loc_tree,NULL); */
  /*                          /\* Wrap_Lk,NULL,loc_tree,NULL); *\/ */
  /*       } */
  /*   } */

  lk_end = loc_tree->c_lnL;


  if(lk_end < lk_begin - tree->mod->s_opt->min_diff_lk_local)
    {
      PhyML_Printf("\n== prop_min: %f prop_max: %f l: %f var:%f",prop_min,prop_max,b_fcus->l->v,b_fcus->l_var->v);
      PhyML_Printf("\n== lk_beg = %f lk_end = %f",lk_begin, lk_end);
      PhyML_Printf("\n== Err. in file %s at line %d",__FILE__,__LINE__);
      Exit("\n");
    }


  return loc_tree->c_lnL;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Round_Optimize(t_tree *tree, calign *data, int n_round_max)
{
  int n_round,each;
  phydbl lk_old, lk_new;

  lk_new = tree->c_lnL;
  lk_old = UNLIKELY;
  n_round = 0;
  each = 0;

  Set_Both_Sides(NO,tree); /* Only the down partial likelihoods need to be up-to-date here */
  Lk(NULL,tree);

  while(n_round < n_round_max)
    {

      if(tree->mod->s_opt->opt_bl || tree->mod->s_opt->constrained_br_len)
        Optimize_Br_Len_Serie(tree);

      if((tree->mod->s_opt->opt_bl || tree->mod->s_opt->constrained_br_len) &&
         (tree->mod->s_opt->print) &&
         (!tree->io->quiet)) Print_Lk(tree,"[Branch lengths     ]");

      Set_Both_Sides(YES,tree);
      Lk(NULL,tree);

      if(!each)
        {
          each = 1;
          Optimiz_All_Free_Param(tree,(tree->io->quiet)?(0):(tree->mod->s_opt->print));
          Set_Both_Sides(YES,tree);
          Lk(NULL,tree);
        }

      lk_new = tree->c_lnL;
      if(lk_new < lk_old - tree->mod->s_opt->min_diff_lk_local)
        {
          PhyML_Printf("\n== lk_new = %f lk_old = %f diff = %f",lk_new,lk_old,lk_new-lk_old);
          Exit("\n== Optimisation failed ! (Round_Optimize)\n");
        }
      if(FABS(lk_new - lk_old) < tree->mod->s_opt->min_diff_lk_local) {
//          DUMP_D(FABS(lk_new - lk_old));
//          DUMP_I(n_round);
          break;
      }
      else {
          lk_old  = lk_new;
      }
      n_round++;
      each--;
    }

  Optimiz_All_Free_Param(tree,(tree->io->quiet)?(0):(tree->mod->s_opt->print));

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Optimize_Br_Len_Serie(t_tree *tree)
{
  if(tree->mod->gamma_mgf_bl == YES)
    {
      phydbl lk_init = tree->c_lnL;

      Generic_Brent_Lk(&(tree->mod->l_var_sigma),
                       tree->mod->l_var_min,
                       tree->mod->l_var_max,
                       tree->mod->s_opt->min_diff_lk_local,
                       tree->mod->s_opt->brent_it_max,
                       tree->mod->s_opt->quickdirty,
                       Wrap_Lk,NULL,tree,NULL,NO);

      if(tree->c_lnL < lk_init - tree->mod->s_opt->min_diff_lk_local)
        {
          PhyML_Printf("\n== %f -- %f",lk_init,tree->c_lnL);
          PhyML_Printf("\n== Err. in file %s at line %d\n",__FILE__,__LINE__);
        }

      if((tree->io->quiet)?(0):(tree->mod->s_opt->print))
        {
          Print_Lk(tree,"[Branch len. var.   ]");
          PhyML_Printf("[%10f]",tree->mod->l_var_sigma);
        }
    }

  if(tree->n_root && tree->ignore_root == NO)
    {
      Update_P_Lk(tree,tree->n_root->b[1],tree->n_root);
      Optimize_Br_Len_Serie_Post(tree->n_root,tree->n_root->v[1],tree->n_root->b[1],tree);
      Update_P_Lk(tree,tree->n_root->b[2],tree->n_root);
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
      Optimize_Br_Len_Serie_Post(tree->a_nodes[0],tree->a_nodes[0]->v[0],tree->a_nodes[0]->b[0],tree);
    }
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void Optimize_Br_Len_Multiplier(t_tree *mixt_tree, int verbose)
{
  phydbl lk_init;
  t_tree *tree;



  tree = mixt_tree;
  
  do
    {      
      if(tree->mod->s_opt->opt_br_len_mult == YES)
        {
          lk_init = Get_Lk(tree);
          Generic_Brent_Lk(&(tree->mod->br_len_mult_unscaled->v),
                           1.E-2,1.E+1,
                           tree->mod->s_opt->min_diff_lk_local,
                           tree->mod->s_opt->brent_it_max,
                           tree->mod->s_opt->quickdirty,
                           Wrap_Lk,NULL,mixt_tree,NULL,NO);

          if(Get_Lk(tree) < lk_init - tree->mod->s_opt->min_diff_lk_local)
            {
              PhyML_Printf("\n== %f -- %f",lk_init,tree->c_lnL);
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
  phydbl l_infa,l_infb;
  phydbl lk_init;

  lk_init = tree->c_lnL;

  if(tree->mod->s_opt->constrained_br_len == YES)
    {
      Generic_Brent_Lk(&(tree->mod->br_len_mult->v),
                       1.E-2,1.E+1,
                       tree->mod->s_opt->min_diff_lk_local,
                       tree->mod->s_opt->brent_it_max,
                       tree->mod->s_opt->quickdirty,
                       Wrap_Lk,NULL,tree,NULL,NO);

      if(tree->c_lnL < lk_init - tree->mod->s_opt->min_diff_lk_local)
        {
          PhyML_Printf("\n== %f -- %f",lk_init,tree->c_lnL);
          Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
        }

      return;
    }

  l_infa = tree->mod->l_max/b_fcus->l->v;
  l_infb = tree->mod->l_min/b_fcus->l->v;

  Set_Both_Sides(YES,tree);

  if(tree->io->mod->s_opt->opt_bl == YES) Br_Len_Brent(l_infb,l_infa,b_fcus,tree);

  if(tree->c_lnL < lk_init - tree->mod->s_opt->min_diff_lk_local)
    {
      PhyML_Printf("\n== %f %f %G",l_infa,l_infb,b_fcus->l->v);
      PhyML_Printf("\n== %f -- %f",lk_init,tree->c_lnL);
      PhyML_Printf("\n== Edge: %d",b_fcus->num);
      PhyML_Printf("\n== is_mixt_tree: %d",tree->is_mixt_tree);
      PhyML_Printf("\n== Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
      Warn_And_Exit("");
    }

  if(d->tax) return;

  if(tree->n_root)
    {
      For(i,3)
        {
          if(d->v[i] != a && d->b[i] != tree->e_root)
            {
              Update_P_Lk(tree,d->b[i],d);
              Optimize_Br_Len_Serie_Post(d,d->v[i],d->b[i],tree);
            }
        }
      
      For(i,3)
        if(d->v[i] == a || d->b[i] == tree->e_root) 
          Update_P_Lk(tree,d->b[i],d);
    }
  else
    {
      For(i,3)
        {
          if(d->v[i] != a)
            {
              Update_P_Lk(tree,d->b[i],d);
              Optimize_Br_Len_Serie_Post(d,d->v[i],d->b[i],tree);
            }
        }
      
      For(i,3) 
        if(d->v[i] == a) 
          Update_P_Lk(tree,d->b[i],d);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Optimiz_Ext_Br(t_tree *tree)
{
  int i;
  t_edge *b,*ori;
  phydbl l_infa,l_infb;
  phydbl lk_init,l_init;

  lk_init = tree->c_lnL;

  For(i,2*tree->n_otu-3)
    {
      b = tree->a_edges[i];
      if((b->left->tax) || (b->rght->tax))
    {

      l_init = b->l->v;

/* 	  Fast_Br_Len(b,tree); */
/* 	  lk = Lk(NULL,tree,b); */

      l_infb = tree->mod->l_min/b->l->v;
      l_infa = 10.;

      Br_Len_Brent(l_infb,l_infa,b,tree);

      ori = b;
      do
        {
          b->nni->best_l    = b->l->v;
          b->nni->l0        = b->l->v;
          b->nni->best_conf = 0;
          b->l->v              = l_init;
              b = b->next;
        }
      while(b);
      b = ori;
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

  if(tree->mod->use_m4mod)
    {
      int failed,i;

      if(tree->mod->s_opt->opt_cov_delta)
        {
          Switch_Eigen(YES,tree->mod);

    /* 	  Optimize_Single_Param_Generic(tree,&(tree->mod->m4mod->delta), */
    /* 					0.01,10., */
    /* 					tree->mod->s_opt->min_diff_lk_local, */
    /* 					tree->mod->s_opt->brent_it_max, */
    /* 					tree->mod->s_opt->quickdirty); */

          Generic_Brent_Lk(&(tree->mod->m4mod->delta),
                   0.01,10.,
                   tree->mod->s_opt->min_diff_lk_local,
                   tree->mod->s_opt->brent_it_max,
                   tree->mod->s_opt->quickdirty,
                   Wrap_Lk,NULL,tree,NULL,NO);

          if(verbose)
            {
              Print_Lk(tree,"[Switching param.   ]");
              PhyML_Printf("[%10f]",tree->mod->m4mod->delta);
            }

          Switch_Eigen(NO,tree->mod);

        }

      
      if(tree->mod->s_opt->opt_cov_free_rates)
        {
          int rcat;
          
          Switch_Eigen(YES,tree->mod);
          
          For(rcat,tree->mod->m4mod->n_h)
            {
              /* 	      Optimize_Single_Param_Generic(tree,&(tree->mod->m4mod->multipl_unscaled[rcat]), */
              /* 					    .01,10., */
              /* 					    tree->mod->s_opt->min_diff_lk_local, */
              /* 					    tree->mod->s_opt->brent_it_max, */
              /* 					    tree->mod->s_opt->quickdirty); */
              
              Generic_Brent_Lk(&(tree->mod->m4mod->multipl_unscaled[rcat]),
                               0.1,100.,
                               tree->mod->s_opt->min_diff_lk_local,
                               tree->mod->s_opt->brent_it_max,
                               tree->mod->s_opt->quickdirty,
                               Wrap_Lk,NULL,tree,NULL,NO);
              
              if(verbose)
                {
                  Print_Lk(tree,"[Rel. subst. rate   ]");
                  PhyML_Printf("[%10f]",tree->mod->m4mod->multipl[rcat]);
                }
            }
          
          For(rcat,tree->mod->m4mod->n_h)
            {
              
              /*  	      Optimize_Single_Param_Generic(tree,&(tree->mod->m4mod->h_fq_unscaled[rcat]), */
              /* 					    .01,100., */
              /* 					    tree->mod->s_opt->min_diff_lk_local, */
              /* 					    tree->mod->s_opt->brent_it_max, */
              /* 					    tree->mod->s_opt->quickdirty); */
              
              Generic_Brent_Lk(&(tree->mod->m4mod->h_fq_unscaled[rcat]),
                               0.1,100.,
                               tree->mod->s_opt->min_diff_lk_local,
                               tree->mod->s_opt->brent_it_max,
                               tree->mod->s_opt->quickdirty,
                               Wrap_Lk,NULL,tree,NULL,NO);
              
              
              if(verbose)
                {
                  Print_Lk(tree,"[Subst. class freq  ]");
                  PhyML_Printf("[%10f]",tree->mod->m4mod->h_fq[rcat]);
                }
            }
          
          Switch_Eigen(NO,tree->mod);
          
        }
      
      if(tree->mod->s_opt->opt_cov_alpha)
        {
          
          Switch_Eigen(YES,tree->mod);
          
          /* 	  Optimize_Single_Param_Generic(tree,&(tree->mod->m4mod->ras->alpha), */
          /* 					.01,10., */
          /* 					tree->mod->s_opt->min_diff_lk_local, */
          /* 					tree->mod->s_opt->brent_it_max, */
          /* 					tree->mod->s_opt->quickdirty); */
          
          Generic_Brent_Lk(&(tree->mod->m4mod->alpha),
                           0.01,10.,
                           tree->mod->s_opt->min_diff_lk_local,
                           tree->mod->s_opt->brent_it_max,
                           tree->mod->s_opt->quickdirty,
                           Wrap_Lk,NULL,tree,NULL,NO);
          
          
          if(verbose)
            {
              Print_Lk(tree,"[Alpha (covarion)   ]");
              PhyML_Printf("[%10f]",tree->mod->m4mod->alpha);
            }
          
          Switch_Eigen(NO,tree->mod);
          
        }
      
      
      /* Substitutions between nucleotides are considered to follow a
         GTR model */
      
      if(tree->mod->io->datatype == NT)
        {
          if(tree->mod->whichmodel == GTR || tree->mod->whichmodel == CUSTOM)
            {
              
              Switch_Eigen(YES,tree->mod);
              
              For(i,5) tree->mod->m4mod->o_rr[i] = LOG(tree->mod->m4mod->o_rr[i]);
              
              failed = YES;
              
              /* BFGS(tree,tree->mod->m4mod->o_rr,5,1.e-5,tree->mod->s_opt->min_diff_lk_local,1.e-5,NO,YES, */
              BFGS(tree,tree->mod->m4mod->o_rr,5,1.e-5,tree->mod->s_opt->min_diff_lk_local,1.e-5,YES,NO,
                   &Return_Abs_Lk,
                   &Num_Derivative_Several_Param,
                   &Lnsrch,&failed);
              
              For(i,5) tree->mod->m4mod->o_rr[i] = EXP(tree->mod->m4mod->o_rr[i]);
              
              For(i,5)
                {
                  /* 	      Optimize_Single_Param_Generic(tree,&(tree->mod->m4mod->o_rr[i]), */
                  /* 					    1.E-20,1.E+10, */
                  /* 					    tree->mod->s_opt->min_diff_lk_local, */
                  /* 					    tree->mod->s_opt->brent_it_max, */
                  /* 					    tree->mod->s_opt->quickdirty); */
                  
                  Generic_Brent_Lk(&(tree->mod->m4mod->o_rr[i]),
                                   1.E-4,1.E+4,
                                   tree->mod->s_opt->min_diff_lk_local,
                                   tree->mod->s_opt->brent_it_max,
                                   tree->mod->s_opt->quickdirty,
                                   Wrap_Lk,NULL,tree,NULL,NO);
                  
                }
              
              if(verbose) Print_Lk(tree,"[GTR parameters     ]");
              
              Switch_Eigen(NO,tree->mod);
            }
        }
    }
  
  Set_Both_Sides(init_both_sides,tree);

  if(tree->both_sides == YES) Lk(NULL,tree); /* Needed to update all partial likelihoods */


  /* if(tree->next) Optimiz_All_Free_Param(tree->next,verbose); */
  /* else            Optimiz_All_Free_Param(tree->next,verbose); */

  /* if(tree->nextree) Optimiz_All_Free_Param(tree->nextree,verbose); */
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
          int logt,
          int is_positive,
          phydbl(*func)(t_tree *tree),
          int(*dfunc)(t_tree *tree,phydbl *param,int n_param,phydbl stepsize,int logt, phydbl(*func)(t_tree *tree),phydbl *derivatives, int is_positive),
          int(*lnsrch)(t_tree *tree, int n, phydbl *xold, phydbl fold, phydbl *g, phydbl *p, phydbl *x,phydbl *f, phydbl stpmax, int *check, int logt, int is_positive),
          int *failed)
{
  
  int check,i,its,j;
  phydbl den,fac,fad,fae,fp,stpmax,sum=0.0,sumdg,sumxi,temp,test,fret;
  phydbl *dg,*g,*hdg,**hessin,*pnew,*xi;
  phydbl fp_old;
  phydbl *init,*sign;

  hessin = (phydbl **)mCalloc(n,sizeof(phydbl *));
  For(i,n) hessin[i] = (phydbl *)mCalloc(n,sizeof(phydbl));
  dg   = (phydbl *)mCalloc(n,sizeof(phydbl ));
  g    = (phydbl *)mCalloc(n,sizeof(phydbl ));
  pnew = (phydbl *)mCalloc(n,sizeof(phydbl ));
  hdg  = (phydbl *)mCalloc(n,sizeof(phydbl ));
  xi   = (phydbl *)mCalloc(n,sizeof(phydbl ));
  init = (phydbl *)mCalloc(n,sizeof(phydbl ));
  sign = (phydbl *)mCalloc(n,sizeof(phydbl ));


  For(i,n) init[i] = p[i];


  /*! p is log transformed */
  if(logt == YES) For(i,n) p[i] = EXP(MIN(1.E+2,p[i]));
  fp=(*func)(tree);
  if(logt == YES) For(i,n) p[i] = LOG(p[i]);
  
  /* PhyML_Printf("\n. ENTER BFGS WITH: %f\n",fp); */
  
  fp_old = fp;
  
  (*dfunc)(tree,p,n,step_size,logt,func,g,is_positive);
  
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
      /* PhyML_Printf("\n. BFGS -> %f stpmax: %f\n",tree->c_lnL,stpmax); */
      
      lnsrch(tree,n,p,fp,g,xi,pnew,&fret,stpmax,&check,logt,is_positive);
      
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
      if (test < TOLX || (FABS(fp-fp_old) < difff && its > 1))
        {
          if(fp > fp_old)
            {
              For(i,n) p[i] = init[i];
              *failed = YES;
            }
          
          if(logt == YES) For(i,n) p[i] = EXP(MIN(1.E+2,p[i]));
          For(i,n) sign[i] = p[i] > .0 ? 1. : -1.;
          if(is_positive == YES) For(i,n) p[i] = FABS(p[i]);
          (*func)(tree);
          if(is_positive == YES) For(i,n) p[i] *= sign[i];
          if(logt == YES) For(i,n) p[i] = LOG(p[i]);
          
          if(is_positive == YES) For(i,n) p[i] = FABS(p[i]);
          
          For(i,n) Free(hessin[i]);
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
      
      (*dfunc)(tree,p,n,step_size,logt,func,g, is_positive);
      
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
          if(logt == YES) For(i,n) p[i] = EXP(MIN(1.E+2,p[i]));
          For(i,n) sign[i] = p[i] > .0 ? 1. : -1.;
          if(is_positive == YES) For(i,n) p[i] = FABS(p[i]);
          (*func)(tree);
          if(is_positive == YES) For(i,n) p[i] *= sign[i];
          if(logt == YES) For(i,n) p[i] = LOG(p[i]);
          
          if(is_positive == YES) For(i,n) p[i] = FABS(p[i]);
          
          For(i,n) Free(hessin[i]);
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
  For(i,n) Free(hessin[i]);
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
                     int logt,
                     int is_positive,
                     phydbl(*func)(t_tree *tree),
                     int(*dfunc_nonaligned)(t_tree *tree,phydbl **param,int n_param,phydbl stepsize,int logt,phydbl(*func)(t_tree *tree),phydbl *derivatives, int is_positive),
                     int(*lnsrch_nonaligned)(t_tree *tree, int n, phydbl **xold, phydbl fold,phydbl *g, phydbl *p, phydbl *x,phydbl *f, phydbl stpmax, int *check, int logt, int is_positive),
                     int *failed)
{

  int check,i,its,j;
  phydbl den,fac,fad,fae,fp,stpmax,sum=0.0,sumdg,sumxi,temp,test,fret;
  phydbl *dg,*g,*hdg,**hessin,*pnew,*xi;
  phydbl fp_old;
  phydbl *init,*sign;

  hessin = (phydbl **)mCalloc(n,sizeof(phydbl *));
  For(i,n) hessin[i] = (phydbl *)mCalloc(n,sizeof(phydbl));
  dg   = (phydbl *)mCalloc(n,sizeof(phydbl ));
  g    = (phydbl *)mCalloc(n,sizeof(phydbl ));
  pnew = (phydbl *)mCalloc(n,sizeof(phydbl ));
  hdg  = (phydbl *)mCalloc(n,sizeof(phydbl ));
  xi   = (phydbl *)mCalloc(n,sizeof(phydbl ));
  init = (phydbl *)mCalloc(n,sizeof(phydbl ));
  sign = (phydbl *)mCalloc(n,sizeof(phydbl ));

  For(i,n) init[i] = (*(p[i]));

  if(logt == YES) For(i,n) (*(p[i])) = EXP(MIN(1.E+2,*(p[i])));
  fp=(*func)(tree);
  if(logt == YES) For(i,n) (*(p[i])) = LOG(*(p[i]));

  fp_old = fp;

  /* PhyML_Printf("\n- ENTER BFGS WITH: %f\n",fp); */

  (*dfunc_nonaligned)(tree,p,n,step_size,logt,func,g,is_positive);

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
      lnsrch_nonaligned(tree,n,p,fp,g,xi,pnew,&fret,stpmax,&check,logt,is_positive);

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

      if(fp > fp_old)
        {
              For(i,n) (*(p[i])) = init[i];
          *failed = 1;
        }

          if(logt == YES) For(i,n) (*(p[i])) = EXP(MIN(1.E+2,*(p[i])));
          For(i,n) sign[i] = *(p[i]) > .0 ? 1. : -1.;
          if(is_positive == YES) For(i,n) *(p[i]) = FABS(*(p[i]));
      (*func)(tree);
          if(is_positive == YES) For(i,n) *(p[i]) *= sign[i];
          if(logt == YES) For(i,n) (*(p[i])) = LOG(*(p[i]));

          if(is_positive == YES) For(i,n) *(p[i]) = FABS(*(p[i]));

      For(i,n) Free(hessin[i]);
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

      (*dfunc_nonaligned)(tree,p,n,step_size,logt,func,g,is_positive);

      test=0.0;
      den=MAX(fret,1.0);
      for (i=0;i<n;i++)
    {
      temp=g[i]*MAX(*(p[i]),1.0)/den;
      if (temp > test) test=temp;
    }

      if (test < gtol)
    {
          if(logt == YES) For(i,n) (*(p[i])) = EXP(MIN(1.E+2,*(p[i])));
          For(i,n) sign[i] = *(p[i]) > .0 ? 1. : -1.;
          if(is_positive == YES) For(i,n) *(p[i]) = FABS(*(p[i]));
      (*func)(tree);
          if(is_positive == YES) For(i,n) *(p[i]) *= sign[i];
          if(logt == YES) For(i,n) (*(p[i])) = LOG(*(p[i]));

          if(is_positive == YES) For(i,n) *(p[i]) = FABS(*(p[i]));

      For(i,n) Free(hessin[i]);
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
  PhyML_Printf("\n. Too many iterations in BFGS...\n");
  *failed = YES;
  For(i,n) Free(hessin[i]);
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
           phydbl *f, phydbl stpmax, int *check, int logt, int is_positive)
{
  int i;
  phydbl a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp,test,tmplam;
  phydbl *local_xold,*sign;
  
  alam = alam2 = f2 = fold2 = tmplam = .0;

  local_xold = (phydbl *)mCalloc(n,sizeof(phydbl));
  sign       = (phydbl *)mCalloc(n,sizeof(phydbl));

  For(i,n) local_xold[i] = xold[i];
  
  *check=0;
  for(sum=0.0,i=0;i<n;i++) sum += p[i]*p[i];
  sum=SQRT(sum);
  /* PhyML_Printf("\n. lnsrch sum: %f",sum); */
  if(sum > stpmax) For(i,n) p[i] *= stpmax/sum;
  /* For(i,n) PhyML_Printf("\n. lnsrch p[i]: %f",p[i]); */
  slope=0.0;
  For(i,n) slope += g[i]*p[i];
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
      For(i,n)
        {
          x[i]=local_xold[i]+alam*p[i];
          xold[i] = x[i];
          /* PhyML_Printf("\n. lnsrch x[i]: %f",x[i]); */
        }
      
      /* PhyML_Printf("\n. lnsrch loop slope: %f alam: %f alam2: %f",slope,alam,alam2); */
      
      if(i==n)
        {
          if(logt == YES) For(i,n) xold[i] = EXP(MIN(1.E+2,xold[i]));
          For(i,n) sign[i] = xold[i] < .0 ? -1. : 1.;
          if(is_positive == YES) For(i,n) xold[i] = FABS(xold[i]);
          /* For(i,n) PhyML_Printf("\n. <<>> %f",xold[i]); */
          *f=Return_Abs_Lk(tree);
          if(is_positive == YES) For(i,n) xold[i] *= sign[i];
          if(logt == YES) For(i,n) xold[i] = LOG(xold[i]);
        }
      else *f=1.+fold+ALF*alam*slope;
      if (alam < alamin)
        {
          *check=1;
          For(i,n) xold[i] = local_xold[i];
          if(is_positive == YES) For(i,n) xold[i] = FABS(xold[i]);
          Free(local_xold);
          Free(sign);
          return 0;
        }
      else if (*f <= fold+ALF*alam*slope)
        {
          For(i,n) xold[i] = local_xold[i];
          if(is_positive == YES) For(i,n) xold[i] = FABS(xold[i]);
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
  Free(sign);
  Free(local_xold);
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
                      phydbl *f, phydbl stpmax, int *check, int logt, int is_positive)
{
  int i;
  phydbl a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp,test,tmplam;
  phydbl *local_xold,*sign;


  alam = alam2 = f2 = fold2 = tmplam = .0;

  local_xold = (phydbl *)mCalloc(n,sizeof(phydbl));
  sign       = (phydbl *)mCalloc(n,sizeof(phydbl));

  For(i,n) local_xold[i] = *(xold[i]);

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
          if(logt == YES) For(i,n) *(xold[i]) = EXP(MIN(1.E+2,*(xold[i])));
          For(i,n) sign[i]    = *(xold[i]) < .0 ? -1. : 1.;
          if(is_positive == YES) For(i,n) *(xold[i]) = FABS(*(xold[i]));
      *f=Return_Abs_Lk(tree);
          if(is_positive == YES) For(i,n) *(xold[i]) *= sign[i];
          if(logt == YES) For(i,n) *(xold[i]) = LOG(*(xold[i]));
    }
      else *f=1.+fold+ALF*alam*slope;

      if (alam < alamin)
    {
      *check=1;
      For(i,n) *(xold[i]) = local_xold[i];
          if(is_positive == YES) For(i,n) *(xold[i]) = FABS(*(xold[i]));
      Free(local_xold);
          Free(sign);
      return 0;
    }
      else if (*f <= fold+ALF*alam*slope)
    {
      For(i,n) *(xold[i]) = local_xold[i];
          if(is_positive) For(i,n) *(xold[i]) = FABS(*(xold[i]));
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
  fw = fv = fx = -Lk_Dist(F,FABS(bx),mod);
  curr_lnL = init_lnL = -fw;

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
          /*                   PhyML_Printf("Golden section step\n"); */
        }
      else
        {
          d=p/q;
          u=x+d;
          if (u-a < tol2 || b-u < tol2)
        d=SIGN(tol1,xm-x);
          /*                   PhyML_Printf("Parabolic step\n"); */
        }
        }
      else
    {
      d=BRENT_CGOLD*(e=(x >= xm ? a-x : b-x));
      /*               PhyML_Printf("Golden section step (default)\n"); */
    }

      u=(FABS(d) >= tol1 ? x+d : x+SIGN(tol1,d));
      (*param) = FABS(u);
      old_lnL = curr_lnL;
      fu = -Lk_Dist(F,FABS(u),mod);
      curr_lnL = -fu;
/*       PhyML_Printf("param=%f loglk=%f\n",*param,fu); */

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
  Exit("\n. Too many iterations in Dist_F_Brent !");
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
  t_tree *tree;
  int    iter;
  phydbl best_alpha, best_pinv, best_mult;
  phydbl slope, intercept;
  phydbl lk_b, lk_a;
  phydbl f0,f1,f2,x0,x1,x2,x3;
  phydbl pinv0, pinv1;
  phydbl a, b, c;
  phydbl fa, fb, fc;
  phydbl K;
  phydbl alpha0, alpha1;
  phydbl best_lnL;
  scalar_dbl **alpha;
  int n_alpha;
  int i;

  Switch_Eigen(NO,mixt_tree->mod);

  alpha    = NULL;
  n_alpha  = 0;
  tree     = mixt_tree;

  do
    {
      For(i,n_alpha) if(tree->mod->ras->alpha == alpha[i]) break;

      if(i == n_alpha)
    {
      if(!alpha) alpha = (scalar_dbl **)mCalloc(1,sizeof(scalar_dbl *));
      else       alpha = (scalar_dbl **)mRealloc(alpha,n_alpha+1,sizeof(scalar_dbl *));

          alpha[n_alpha] = tree->mod->ras->alpha;
      n_alpha++;

      if((tree->mod->s_opt->opt_pinvar) && (tree->mod->s_opt->opt_alpha) && (tree->mod->ras->n_catg > 1))
        {

          lk_b     = UNLIKELY;
          lk_a     = UNLIKELY;

          /* PhyML_Printf("\n\n. %p Init lnL = %f alpha=%f pinv=%f", */
              /*              tree, */
              /*              mixt_tree->c_lnL, */
              /*              tree->mod->ras->alpha, */
              /*              tree->mod->ras->pinvar->v); */

          /* Two (full) steps to compute  pinv_alpha_slope & pinv_alpha_intercept */

          Set_Both_Sides(YES,mixt_tree);
          Lk(NULL,mixt_tree);
          lk_b = mixt_tree->c_lnL;

          Optimize_Br_Len_Serie(mixt_tree);

          Set_Both_Sides(NO,mixt_tree);

          Optimize_Single_Param_Generic(mixt_tree,&(tree->mod->ras->alpha->v),0.01,100.,
                        mixt_tree->mod->s_opt->min_diff_lk_local,
                        mixt_tree->mod->s_opt->brent_it_max,
                        mixt_tree->mod->s_opt->quickdirty);

          Optimize_Single_Param_Generic(mixt_tree,&(tree->mod->ras->pinvar->v),.0001,0.9999,
                        tree->mod->s_opt->min_diff_lk_local,
                        tree->mod->s_opt->brent_it_max,
                        tree->mod->s_opt->quickdirty);

          pinv0  = tree->mod->ras->pinvar->v;
          alpha0 = tree->mod->ras->alpha->v;
          f0 = mixt_tree->c_lnL;

          Set_Both_Sides(YES,mixt_tree);
          Lk(NULL,mixt_tree);

          Optimize_Br_Len_Serie(mixt_tree);

          Set_Both_Sides(NO,mixt_tree);
          Optimize_Single_Param_Generic(mixt_tree,&(tree->mod->ras->alpha->v),0.01,100.,
                        tree->mod->s_opt->min_diff_lk_local,
                        tree->mod->s_opt->brent_it_max,
                        tree->mod->s_opt->quickdirty);

          Optimize_Single_Param_Generic(mixt_tree,&(tree->mod->ras->pinvar->v),.0001,0.9999,
                        tree->mod->s_opt->min_diff_lk_local,
                        tree->mod->s_opt->brent_it_max,
                        tree->mod->s_opt->quickdirty);

          lk_a = mixt_tree->c_lnL;

          pinv1  = tree->mod->ras->pinvar->v;
          alpha1 = tree->mod->ras->alpha->v;
          f1 = mixt_tree->c_lnL;
          best_lnL = f1;

          if(lk_a < lk_b - mixt_tree->mod->s_opt->min_diff_lk_local)
            {
              PhyML_Printf("\n== Err. in file %s at line %d\n",__FILE__,__LINE__);
              Exit("\n");
            }
          else if(FABS(lk_a - lk_b) < mixt_tree->mod->s_opt->min_diff_lk_local)
            {
              if(alpha) Free(alpha);
              return 1;
            }

          Record_Br_Len(mixt_tree);
          best_alpha = tree->mod->ras->alpha->v;
          best_pinv  = tree->mod->ras->pinvar->v;
          best_mult  = tree->mod->br_len_mult->v;

          /* PhyML_Printf("\n\n. Init lnL after std opt = %f [%f] best_alpha=%f best_pinv=%f",mixt_tree->c_lnL,Lk(NULL,mixt_tree),best_alpha,best_pinv); */
          /* PhyML_Printf("\n. Best_lnL = %f %d",best_lnL,tree->mod->ras->invar); */

          slope     = (pinv1 - pinv0)/(alpha1 - alpha0);
          intercept = pinv1 - slope * alpha1;


              /* printf("\n. slope = %f pinv1=%f pinv0=%f alpha1=%f alpha0=%f", */
              /*        slope,pinv1,pinv0,alpha1,alpha0); */


          if((slope > 0.001) && (slope < 1./0.001))
        {
          /* PhyML_Printf("\n. pinv0 = %f, pinv1 = %f, alpha0 = %f, alpha1 = %f",pinv0,pinv1,alpha0,alpha1); */
          /* PhyML_Printf("\n. slope = %f intercept = %f",slope,intercept); */

          K = 0.381966;

          if(alpha1 < alpha0)
            {
              c  = alpha0;
              b  = alpha1;
              fc = f0;
              fb = f1;

              a = (0.1 < alpha1)?(0.1):(0.5*alpha1);
              tree->mod->ras->alpha->v = a;
              tree->mod->ras->pinvar->v = slope * tree->mod->ras->alpha->v + intercept;
              if(tree->mod->ras->pinvar->v > 1.0) tree->mod->ras->pinvar->v = 0.9;
              if(tree->mod->ras->pinvar->v < 0.0) tree->mod->ras->pinvar->v = 0.001;
              Set_Both_Sides(YES,mixt_tree);
              Lk(NULL,mixt_tree);

              Optimize_Br_Len_Serie(mixt_tree);

              fa = mixt_tree->c_lnL;

              iter = 0;

              /* PhyML_Printf("\n. a=%f, b=%f, c=%f, fa=%f, fb=%f, fc=%f (alpha=%f pinv=%f)",a,b,c,fa,fb,fc,tree->mod->ras->alpha->v,tree->mod->ras->pinvar->v); */

              while(fa > fb)
            {
              a = a/5.;
              tree->mod->ras->alpha->v = a;
              tree->mod->ras->pinvar->v = slope * tree->mod->ras->alpha->v + intercept;
              if(tree->mod->ras->pinvar->v > 1.0) tree->mod->ras->pinvar->v = 0.9;
              if(tree->mod->ras->pinvar->v < 0.0) tree->mod->ras->pinvar->v = 0.001;
              Set_Both_Sides(YES,mixt_tree);
              Lk(NULL,mixt_tree);
              Optimize_Br_Len_Serie(mixt_tree);
              fa = mixt_tree->c_lnL;
              /* PhyML_Printf("\n1 a=%f, b=%f, c=%f, fa=%f, fb=%f, fc=%f",a,b,c,fa,fb,fc); */
              if(iter++ > 10)
                            {
                              if(alpha) Free(alpha);
                              return 0;
                            }
            }
            }
          else
            {
              a  = alpha0;
              b  = alpha1;
              fa = f0;
              fb = f1;

              c = (alpha1 < 2.)?(2.0):(2.*alpha1);
              tree->mod->ras->alpha->v = c;
              tree->mod->ras->pinvar->v = slope * tree->mod->ras->alpha->v + intercept;
              if(tree->mod->ras->pinvar->v > 1.0) tree->mod->ras->pinvar->v = 0.9;
              if(tree->mod->ras->pinvar->v < 0.0) tree->mod->ras->pinvar->v = 0.001;
              Set_Both_Sides(YES,mixt_tree);
              Lk(NULL,mixt_tree);
              Optimize_Br_Len_Serie(mixt_tree);
              fc = mixt_tree->c_lnL;

              /* PhyML_Printf("\n. a=%f, b=%f, c=%f, fa=%f, fb=%f, fc=%f (alpha=%f pinv=%f)",a,b,c,fa,fb,fc,tree->mod->ras->alpha->v,tree->mod->ras->pinvar->v); */

              iter = 0;
              while(fc > fb)
            {
              c = c*2.;
              tree->mod->ras->alpha->v = c;
              tree->mod->ras->pinvar->v = slope * tree->mod->ras->alpha->v + intercept;
              if(tree->mod->ras->pinvar->v > 1.0) tree->mod->ras->pinvar->v = 0.9;
              if(tree->mod->ras->pinvar->v < 0.0) tree->mod->ras->pinvar->v = 0.001;
              Set_Both_Sides(YES,mixt_tree);
              Lk(NULL,mixt_tree);
              Optimize_Br_Len_Serie(mixt_tree);
              fc = mixt_tree->c_lnL;
              /* PhyML_Printf("\n2 a=%f, b=%f, c=%f, fa=%f, fb=%f, fc=%f",a,b,c,fa,fb,fc); */
              if(iter++ > 10)
                            {
                              if(alpha) Free(alpha);
                              return 0;
                            }
            }
            }


          if(FABS(b - c) > FABS(a - b))
            {
              x0 = a; x1 = b; x3 = c;
              x2 = b + K * FABS(b - c);

              f0 = fa;
              f1 = fb;
              tree->mod->ras->alpha->v = x2;
              tree->mod->ras->pinvar->v = slope * tree->mod->ras->alpha->v + intercept;
              if(tree->mod->ras->pinvar->v > 1.0) tree->mod->ras->pinvar->v = 0.9;
              if(tree->mod->ras->pinvar->v < 0.0) tree->mod->ras->pinvar->v = 0.001;
              Set_Both_Sides(YES,mixt_tree);
              Lk(NULL,mixt_tree);
              Optimize_Br_Len_Serie(mixt_tree);
              f2 = mixt_tree->c_lnL;
            }
          else /* |b -c| < |a - b| */
            {
              x0 = a; x2 = b; x3 = c;
              x1 = b - K * FABS(b - a);

              f0 = fa;
              f2 = fb;
              tree->mod->ras->alpha->v = x1;
              tree->mod->ras->pinvar->v = slope * tree->mod->ras->alpha->v + intercept;
              if(tree->mod->ras->pinvar->v > 1.0) tree->mod->ras->pinvar->v = 0.9;
              if(tree->mod->ras->pinvar->v < 0.0) tree->mod->ras->pinvar->v = 0.001;
              Set_Both_Sides(YES,mixt_tree);
              Lk(NULL,mixt_tree);
              Optimize_Br_Len_Serie(mixt_tree);
              f1 = mixt_tree->c_lnL;
            }

          iter = 0;
          do
            {
              /* PhyML_Printf("\n. x0=%f, x1=%f, x2=%f, x3=%f, f0=%f, f1=%f, f2=%f, f3=%f", */
              /* 	 x0,x1,x2,x3,f0,f1,f2,f3); */

              if(f1 > f2)
            {
              x3 = x2;
              x2 = x1;
              x1 = x2 - K * FABS(x2 - x0);

              f2 = f1;

              tree->mod->ras->alpha->v = x1;
              tree->mod->ras->pinvar->v = slope * tree->mod->ras->alpha->v + intercept;
              if(tree->mod->ras->pinvar->v > 1.0) tree->mod->ras->pinvar->v = 0.9;
              if(tree->mod->ras->pinvar->v < 0.0) tree->mod->ras->pinvar->v = 0.001;
              Set_Both_Sides(YES,mixt_tree);
              Lk(NULL,mixt_tree);
              Optimize_Br_Len_Serie(mixt_tree);
              f1 = mixt_tree->c_lnL;
              if(f1 > best_lnL)
                {
                  Record_Br_Len(mixt_tree);
                  best_alpha = tree->mod->ras->alpha->v;
                  best_pinv  = tree->mod->ras->pinvar->v;
                  best_mult  = tree->mod->br_len_mult->v;
                  /* PhyML_Printf("\n>.< New alpha=%f pinv=%f",best_alpha,best_pinv); */
                }
              /* PhyML_Printf("\n> f1=%f",f1); */
            }
              else /* f1 < f2 */
            {
              x0 = x1;
              x1 = x2;
              x2 = x2 + K * FABS(x3 - x2);

              f0 = f1;
              f1 = f2;

              tree->mod->ras->alpha->v = x2;
              tree->mod->ras->pinvar->v = slope * tree->mod->ras->alpha->v + intercept;
              if(tree->mod->ras->pinvar->v > 1.0) tree->mod->ras->pinvar->v = 0.9;
              if(tree->mod->ras->pinvar->v < 0.0) tree->mod->ras->pinvar->v = 0.001;
              Set_Both_Sides(YES,mixt_tree);
              Lk(NULL,mixt_tree);
              Optimize_Br_Len_Serie(mixt_tree);
              f2 = mixt_tree->c_lnL;
              if(f2 > best_lnL)
                {
                  Record_Br_Len(mixt_tree);
                  best_alpha = tree->mod->ras->alpha->v;
                  best_pinv  = tree->mod->ras->pinvar->v;
                  best_mult  = tree->mod->br_len_mult->v;
                  /* PhyML_Printf("\n>o< New alpha=%f pinv=%f",best_alpha,best_pinv); */
                }
              /* PhyML_Printf("\n> f2=%f",f2); */
            }

              if(FABS(f1 - f2) < 0.01) break;

              iter++;

            }while(iter < 100);
        }
          
          tree->mod->ras->alpha->v  = best_alpha;
          tree->mod->ras->pinvar->v = best_pinv;
          tree->mod->br_len_mult->v = best_mult;
          Restore_Br_Len(mixt_tree);
          Set_Both_Sides(YES,mixt_tree);
          Lk(NULL,mixt_tree);

          if(verbose)
            {
              Print_Lk(mixt_tree,"[Alpha              ]");
              PhyML_Printf("[%10f]",tree->mod->ras->alpha->v);
              Print_Lk(mixt_tree,"[P-inv              ]");
              PhyML_Printf("[%10f]",tree->mod->ras->pinvar->v);
            }
        }
    }

      tree = tree->next_mixt;

    }
  while(tree);

  /* PhyML_Printf("\n\n. Init lnL after golden opt = %f [%f] best_alpha=%f best_pinv=%f",tree->c_lnL,Lk(tree),best_alpha,best_pinv); */

  if(alpha) Free(alpha);

  return 1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Generic_Brent_Lk(phydbl *param, phydbl ax, phydbl cx, phydbl tol,
                        int n_iter_max, int quickdirty,
                        phydbl (*obj_func)(t_edge *,t_tree *,supert_tree *),
                        t_edge *branch, t_tree *tree, supert_tree *stree, int logt)
{
  int iter;
  phydbl a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  phydbl e=0.0;
  phydbl old_lnL,init_lnL;
  phydbl bx = *param;

  d=0.0;
  a=((ax < cx) ? ax : cx);
  b=((ax > cx) ? ax : cx);
  x=w=v=bx;
  old_lnL = UNLIKELY;
  (*param) = bx;
  if(logt == YES) (*param) = EXP(MIN(1.E+2,*param));
  fw=fv=fx=fu=-(*obj_func)(branch,tree,stree);
  if(logt == YES) (*param) = LOG(*param);
  init_lnL = -fw;

  /* PhyML_Printf("\n. %p %p %p init_lnL = %f a=%f b=%f c=%f",branch,tree,stree,init_lnL,ax,bx,cx); */
  
  for(iter=1;iter<=BRENT_IT_MAX;iter++)
    {
      xm=0.5*(a+b);
      tol2=2.0*(tol1=tol*x+BRENT_ZEPS);
      
      if((fu > init_lnL + tol) && (quickdirty))
        {
          (*param) = x;
          if(logt == YES) (*param) = EXP(MIN(1.E+2,*param));
          fu = (*obj_func)(branch,tree,stree);
          if(logt == YES) (*param) = LOG(*param);
          return fu;
        }

/*       if(((FABS(cur_lnL-old_lnL) < tol) && (cur_lnL > init_lnL - tol)) || (iter > n_iter_max - 1)) */
      if((FABS(fu-old_lnL) < tol) || (iter > n_iter_max - 1))
        {
          (*param) = x;
          if(logt == YES) (*param) = EXP(MIN(1.E+2,*param));
          fu = (*obj_func)(branch,tree,stree);
          if(logt == YES) (*param) = LOG(*param);
    /* 	  Exit("\n"); */
          return fu;
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
              /* PhyML_Printf(" Golden section step\n"); */
            }
          else
            {
              d=p/q;
              u=x+d;
              if (u-a < tol2 || b-u < tol2) d=SIGN(tol1,xm-x);
              /* PhyML_Printf(" Parabolic step [e=%f]\n",e); */
            }
            }
      else
        {
          d=BRENT_CGOLD*(e=(x >= xm ? a-x : b-x));
          /* PhyML_Printf(" Golden section step (default) [e=%f tol1=%f a=%f b=%f d=%f]\n",e,tol1,a,b,d); */
        }

      u=(FABS(d) >= tol1 ? x+d : x+SIGN(tol1,d));
      (*param) = u;
      old_lnL = fu;
      if(logt == YES) (*param) = EXP(MIN(1.E+2,*param));
      fu = -(*obj_func)(branch,tree,stree);
      if(logt == YES) (*param) = LOG(*param);

      /* PhyML_Printf("\n. iter=%d/%d param=%f lnL=%f",iter,BRENT_IT_MAX,*param,fu); */

      if(fu <= fx)
        {
          if(u >= x) a=x; else b=x;
          SHFT(v,w,x,u)
          SHFT(fv,fw,fx,fu)
        }
      else
        {
          if (u < x) a=u; else b=u;
    /* 	  if (fu <= fw || w == x) */
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

  Exit("\n. Too many iterations in Generic_Brent_Lk !");
  return(-1);
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


  printf("\n. cur_lnL = %f new_lnL=%f",cur_lnL,new_lnL);


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
                   Wrap_Lk,NULL,tree,NULL,NO);

      printf("\n. cur_lnL=%f new_lnL=%f clock_r=%G root height=%f",
         cur_lnL,new_lnL,tree->rates->clock_r,tree->rates->nd_t[tree->n_root->num]);
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

  Generic_Brent_Lk(&(tree->rates->nd_t[tree->n_root->num]),
           MIN(tree->rates->t_prior_max[tree->n_root->num],
               MIN(tree->rates->nd_t[tree->n_root->v[2]->num],
               tree->rates->nd_t[tree->n_root->v[1]->num])),
           tree->rates->t_prior_min[tree->n_root->num],
           tree->mod->s_opt->min_diff_lk_local,
           tree->mod->s_opt->brent_it_max,
           tree->mod->s_opt->quickdirty,
           Wrap_Lk,NULL,tree,NULL,NO);
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
      For(i,3)
    if((d->v[i] != a) && (d->b[i] != tree->e_root))
      {
        if(!v2) { v2 = d->v[i]; }
        else    { v3 = d->v[i]; }
      }

      Opt_Node_Heights_Recurr_Pre(d,v2,tree);
      Opt_Node_Heights_Recurr_Pre(d,v3,tree);

      t0 = tree->rates->nd_t[a->num];
      t2 = tree->rates->nd_t[v2->num];
      t3 = tree->rates->nd_t[v3->num];

      t_min = t0;
      t_max = MIN(t2,t3);

      t_min = MAX(t_min,tree->rates->t_prior_min[d->num]);
      t_max = MIN(t_max,tree->rates->t_prior_max[d->num]);

      t_min += tree->rates->min_dt;
      t_max -= tree->rates->min_dt;

      if(t_min > t_max)
        {
          PhyML_Printf("\n== Err in file %s at line %d\n",__FILE__,__LINE__);
          Exit("\n");
        }

      Generic_Brent_Lk(&(tree->rates->nd_t[d->num]),
               t_min,t_max,
               tree->mod->s_opt->min_diff_lk_local,
               tree->mod->s_opt->brent_it_max,
               tree->mod->s_opt->quickdirty,
               Wrap_Lk,NULL,tree,NULL,NO);

      /* printf("\n. t%d = %f [%f;%f] lnL = %f",d->num,tree->rates->nd_t[d->num],t_min,t_max,tree->c_lnL); */

    }
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Optimize_RR_Params(t_tree *mixt_tree, int verbose)
{
  t_tree *tree;
  t_rmat **r_mat;
  int n_r_mat;
  int i;

  Switch_Eigen(YES,mixt_tree->mod);

  n_r_mat = 0;
  tree    = mixt_tree;
  r_mat   = NULL;

  do
    {
      if(tree->next) tree = tree->next;
      
      For(i,n_r_mat) if(tree->mod->r_mat == r_mat[i]) break;
      
      if(i == n_r_mat) // tree->mod->r_mat was not found before
        {
          if(!r_mat) r_mat = (t_rmat **)mCalloc(1,sizeof(t_rmat *));
          else       r_mat = (t_rmat **)mRealloc(r_mat,n_r_mat+1,sizeof(t_rmat *));
          r_mat[n_r_mat] = tree->mod->r_mat;
          n_r_mat++;
          
          if((tree->mod->whichmodel == GTR) ||
             ((tree->mod->whichmodel == CUSTOM) &&
              (tree->mod->s_opt->opt_rr) &&
              (tree->mod->r_mat->n_diff_rr > 1)))
            {
              int failed,i;
              
              For(i,tree->mod->r_mat->n_diff_rr) tree->mod->r_mat->rr_val->v[i] = LOG(tree->mod->r_mat->rr_val->v[i]);
              
              failed = YES;
              
              /* BFGS(mixt_tree,tree->mod->r_mat->rr_val->v,tree->mod->r_mat->n_diff_rr,1.e-5,tree->mod->s_opt->min_diff_lk_local,1.e-5,NO,YES, */
              if(tree->mod->r_mat->n_diff_rr > 2)
                {
                  BFGS(mixt_tree,tree->mod->r_mat->rr_val->v,tree->mod->r_mat->n_diff_rr,1.e-5,tree->mod->s_opt->min_diff_lk_local,1.e-5,YES,NO,
                       &Return_Abs_Lk,
                       &Num_Derivative_Several_Param,
                       &Lnsrch,&failed);
                }

              For(i,tree->mod->r_mat->n_diff_rr) tree->mod->r_mat->rr_val->v[i] = EXP(tree->mod->r_mat->rr_val->v[i]);
              

              if(failed == YES)
                {
                  For(i,tree->mod->r_mat->n_diff_rr)
                    if(i != 5)
                      {
                        Generic_Brent_Lk(&(tree->mod->r_mat->rr_val->v[i]),
                                         1.E-2,1.E+2,
                                         tree->mod->s_opt->min_diff_lk_local,
                                         tree->mod->s_opt->brent_it_max,
                                         tree->mod->s_opt->quickdirty,
                                         Wrap_Lk,NULL,mixt_tree,NULL,NO);
                      }
                }
              
              if(verbose) Print_Lk(tree->mixt_tree?
                                   tree->mixt_tree:
                                   tree,"[GTR parameters     ]");
              
            }
    }

      tree = tree->next;
      if(!tree) break;

  }
  while(1);

  if(r_mat) Free(r_mat);

  Switch_Eigen(NO,mixt_tree->mod);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Optimize_TsTv(t_tree *mixt_tree, int verbose)
{
  scalar_dbl **tstv;
  int n_tstv;
  t_tree *tree;
  int i;

  Switch_Eigen(YES,mixt_tree->mod);

  tstv   = NULL;
  n_tstv = 0;
  tree   = mixt_tree;

  do
    {
      if(tree->next) tree = tree->next;
      
      For(i,n_tstv) if(tree->mod->kappa == tstv[i]) break;
      
      if(i == n_tstv)
        {
          if(!tstv) tstv = (scalar_dbl **)mCalloc(1,sizeof(scalar_dbl *));
          else      tstv = (scalar_dbl **)mRealloc(tstv,n_tstv+1,sizeof(scalar_dbl *));
          tstv[n_tstv] = tree->mod->kappa;
          n_tstv++;
          
          if(tree->mod->s_opt->opt_kappa)
            {
              Generic_Brent_Lk(&(tree->mod->kappa->v),
                               0.1,100.,
                               tree->mod->s_opt->min_diff_lk_local,
                               tree->mod->s_opt->brent_it_max,
                               tree->mod->s_opt->quickdirty,
                               Wrap_Lk,NULL,mixt_tree,NULL,NO);
              
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
  
  Switch_Eigen(NO,mixt_tree->mod);
  
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Optimize_Pinv(t_tree *mixt_tree, int verbose)
{
  scalar_dbl **pinv;
  int n_pinv;
  t_tree *tree;
  int i;

  Switch_Eigen(NO,mixt_tree->mod);

  pinv   = NULL;
  n_pinv = 0;
  tree   = mixt_tree;

  do
    {

      For(i,n_pinv) if(tree->mod->ras->pinvar == pinv[i]) break;

      if(i == n_pinv)
    {
      if(!pinv) pinv = (scalar_dbl **)mCalloc(1,sizeof(scalar_dbl *));
      else      pinv = (scalar_dbl **)mRealloc(pinv,n_pinv+1,sizeof(scalar_dbl *));
      pinv[n_pinv] = tree->mod->ras->pinvar;
      n_pinv++;

      if(tree->mod->s_opt->opt_pinvar == YES && tree->mod->s_opt->opt_alpha == NO)
        {
          Generic_Brent_Lk(&(tree->mod->ras->pinvar->v),
                   0.0001,0.9999,
                   tree->mod->s_opt->min_diff_lk_local,
                   tree->mod->s_opt->brent_it_max,
                   tree->mod->s_opt->quickdirty,
                   Wrap_Lk,NULL,mixt_tree,NULL,NO);

          if(verbose)
        {
          Print_Lk(mixt_tree,"[P-inv              ]");
          PhyML_Printf("[%10f]",tree->mod->ras->pinvar->v);
        }
        }
    }

      tree = tree->next;

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

  Switch_Eigen(NO,mixt_tree->mod);

  alpha   = NULL;
  n_alpha = 0;
  tree   = mixt_tree;

  do
    {
      if(tree->mod->s_opt->opt_alpha == YES && tree->mod->ras->n_catg > 1)
        {
          For(i,n_alpha) if(tree->mod->ras->alpha == alpha[i]) break;

          if(i == n_alpha)
            {
              if(!alpha) alpha = (scalar_dbl **)mCalloc(1,sizeof(scalar_dbl *));
              else       alpha = (scalar_dbl **)mRealloc(alpha,n_alpha+1,sizeof(scalar_dbl *));
              alpha[n_alpha] = tree->mod->ras->alpha;
              n_alpha++;

              if(tree->mod->s_opt->opt_alpha == YES &&
                 tree->mod->ras->free_mixt_rates == NO &&
                 tree->mod->s_opt->opt_pinvar == NO)
                {
                  if(tree->mod->ras->n_catg > 1)
                    {
                      Generic_Brent_Lk(&(tree->mod->ras->alpha->v),
                                       tree->mod->ras->alpha->v/2.,100.,
                                       tree->mod->s_opt->min_diff_lk_local,
                                       tree->mod->s_opt->brent_it_max,
                                       tree->mod->s_opt->quickdirty,
                                       Wrap_Lk,NULL,mixt_tree,NULL,NO);
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
  phydbl lk_before, lk_after;
  tree = mixt_tree;

  lk_before = lk_after = UNLIKELY;

  do
    {
      if((tree->mod->s_opt->opt_free_mixt_rates) && (tree->mod->ras->free_mixt_rates == YES) && (tree->mod->ras->n_catg > 1))
        {
          if(tree->mod->s_opt->serial_free_rates == YES)
            {
              fast = YES;
              lk_before = tree->c_lnL;
              Optimize_Free_Rate_Weights(tree,fast,verbose);
              Optimize_Free_Rate_Rr(tree,fast,verbose);
              lk_after = tree->c_lnL;

              if(lk_after < lk_before - tree->mod->s_opt->min_diff_lk_global)
                {
                  PhyML_Printf("\n== lk_before: %f lk_after: %f diff: %G",lk_before,lk_after,lk_before-lk_after);
                  PhyML_Printf("\n== Err. in file %s at line %d\n",__FILE__,__LINE__);
                  Exit("");
                }
            }

          if(FABS(lk_before - lk_after) > 0.001)
            {
              phydbl **x;
              x = (phydbl **)mCalloc(2*tree->n_otu-3 + 2*tree->mod->ras->n_catg,sizeof(phydbl *));
              pos = 0;
              
              lk_before = tree->c_lnL;

              /* For(i,2*tree->n_otu-3) x[pos++] = &(tree->a_edges[i]->l->v); */
              For(i,tree->mod->ras->n_catg) x[pos++] = tree->mod->ras->gamma_rr_unscaled->v+i;
              For(i,tree->mod->ras->n_catg) x[pos++] = tree->mod->ras->gamma_r_proba_unscaled->v+i;

              /* For(i,2*tree->n_otu-3 + 2*tree->mod->ras->n_catg) *(x[i]) = LOG(MAX(1.E-10,*(x[i]))); */
              /* For(i,2*tree->mod->ras->n_catg) printf("\n:: %12f",*(x[i])); fflush(NULL); */
              For(i,2*tree->mod->ras->n_catg) *(x[i]) = LOG(MAX(1.E-10,*(x[i])));
              /* For(i,2*tree->n_otu-3 + 2*tree->mod->ras->n_catg) printf("\n<> %12f",*(x[i])); */
              /* For(i,2*tree->mod->ras->n_catg) printf("\n<> %12f",*(x[i])); fflush(NULL); */

              failed = YES;
              /* BFGS_Nonaligned(tree,x,2*tree->n_otu-3 + 2*tree->mod->ras->n_catg,1.e-5,tree->mod->s_opt->min_diff_lk_global,1.e-5,YES, */
              BFGS_Nonaligned(tree,x,2*tree->mod->ras->n_catg,1.e-5,tree->mod->s_opt->min_diff_lk_global,1.e-5,YES,NO,
                              &Return_Abs_Lk,
                              &Num_Derivative_Several_Param_Nonaligned,
                              &Lnsrch_Nonaligned,&failed);


              /* For(i,2*tree->n_otu-3 + 2*tree->mod->ras->n_catg) *(x[i]) = EXP(*(x[i])); */
              For(i,2*tree->mod->ras->n_catg) *(x[i]) = EXP(MIN(1.E+2,*(x[i])));

              lk_after = tree->c_lnL;

              /* For(i,2*tree->mod->ras->n_catg) printf("\n>< %12f",*(x[i])); fflush(NULL); */

              if(lk_after < lk_before - tree->mod->s_opt->min_diff_lk_global)
                {
                  PhyML_Printf("\n== lk_before: %f lk_after: %f diff: %G",lk_before,lk_after,lk_before-lk_after);
                  PhyML_Printf("\n== Err. in file %s at line %d\n",__FILE__,__LINE__);
                  Exit("");
                }

              Free(x);
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
      phydbl wm;

      /* tree->mod->s_opt->curr_opt_free_rates = YES; */

      if(fast == YES)
        {
          For(i,tree->mod->ras->n_catg) tree->mod->ras->skip_rate_cat[i] = YES;
          tree->mod->ras->normalise_rr                                   = NO;

          wm = Weighted_Mean(tree->mod->ras->gamma_rr_unscaled->v,
                             tree->mod->ras->gamma_r_proba->v,
                             tree->mod->ras->n_catg);

          tree->mod->ras->free_rate_mr->v = 100.;
          For(i,2*tree->n_otu-1) tree->a_edges[i]->l->v /= (wm * tree->mod->ras->free_rate_mr->v);
        }


      For(i,tree->mod->ras->n_catg-1)
        {
          if(fast == YES) tree->mod->ras->skip_rate_cat[i] = NO;

          Generic_Brent_Lk(&(tree->mod->ras->gamma_rr_unscaled->v[i]),
                           0.0,
                           100.,
                           tree->mod->s_opt->min_diff_lk_local,
                           tree->mod->s_opt->brent_it_max,
                           tree->mod->s_opt->quickdirty,
                           Wrap_Lk,NULL,tree,NULL,NO);

          if(fast == YES) tree->mod->ras->skip_rate_cat[i] = YES;

          lk_after = tree->c_lnL;

          if(lk_after < lk_before - tree->mod->s_opt->min_diff_lk_global)
            {
              PhyML_Printf("\n== lk_before: %f lk_after: %f diff: %G",lk_before,lk_after,lk_before-lk_after);
              PhyML_Printf("\n== Err. in file %s at line %d\n",__FILE__,__LINE__);
              Exit("");
            }
        }

      if(fast == YES)
        {
          For(i,tree->mod->ras->n_catg) tree->mod->ras->skip_rate_cat[i] = NO;
          tree->mod->ras->normalise_rr                                   = YES;

          wm = Weighted_Mean(tree->mod->ras->gamma_rr_unscaled->v,
                             tree->mod->ras->gamma_r_proba->v,
                             tree->mod->ras->n_catg);

          For(i,2*tree->n_otu-1) tree->a_edges[i]->l->v *= (wm * tree->mod->ras->free_rate_mr->v);
        }

      /* tree->mod->s_opt->curr_opt_free_rates = NO; */

    }
  else
    {
      int i;
      For(i,tree->mod->ras->n_catg-1)
        {
          Generic_Brent_Lk(&(tree->mod->ras->gamma_rr_unscaled->v[i]),
                           1.E-2,100.,
                           tree->mod->s_opt->min_diff_lk_local,
                           tree->mod->s_opt->brent_it_max,
                           tree->mod->s_opt->quickdirty,
                           Wrap_Lk,NULL,tree,NULL,NO);
        }
    }

  lk_after = tree->c_lnL;

  if(lk_after < lk_before - tree->mod->s_opt->min_diff_lk_global)
    {
      PhyML_Printf("\n== lk_before: %f lk_after: %f diff: %G",lk_before,lk_after,lk_before-lk_after);
      PhyML_Printf("\n== Err. in file %s at line %d\n",__FILE__,__LINE__);
      Exit("");
    }

  if(verbose) Print_Lk(tree,"[Rate class values  ]");

  /* int i; */
  /* For(i,tree->mod->ras->n_catg) */
  /*   { */
  /*     printf("\n+ c %2d p: %15f r: %15f up: %15f ur: %5f", */
  /*            i+1, */
  /*            tree->mod->ras->gamma_r_proba->v[i], */
  /*            tree->mod->ras->gamma_rr->v[i], */
  /*            tree->mod->ras->gamma_r_proba_unscaled->v[i], */
  /*            tree->mod->ras->gamma_rr_unscaled->v[i]); */
  /*   } */
  /* fflush(NULL); */

  /* printf("\n. LK: %f",Lk(NULL,tree)); */

  /* int i; */
  /* printf("\n"); */
  /* printf("X*X %f ",tree->c_lnL); */
  /* /\* For(i,tree->mod->ras->n_catg) printf("%f ",tree->mod->ras->gamma_rr_unscaled->v[i]); *\/ */
  /* For(i,tree->mod->ras->n_catg) printf("%f ",tree->mod->ras->gamma_rr->v[i]); */
  /* For(i,tree->mod->ras->n_catg) printf("%f ",tree->mod->ras->gamma_r_proba->v[i]); */
  /* For(i,2*tree->n_otu-3) printf("%f ",tree->a_edges[i]->l->v); */
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Optimize_Free_Rate_Weights(t_tree *tree, int fast, int verbose)
{
  int i;
  phydbl wm;
  phydbl lk_before, lk_after;


  if(tree->mod->s_opt->first_opt_free_mixt_rates == YES)
    {
      /* tree->mod->s_opt->opt_alpha     = YES; */
      /* tree->mod->s_opt->opt_pinvar    = NO; */
      /* tree->mod->ras->free_mixt_rates = NO; */

      /* Optimize_Alpha(tree,YES); */

      /* For(i,tree->mod->ras->n_catg) */
      /*   { */
      /*     tree->mod->ras->gamma_r_proba_unscaled->v[i] = 1./(phydbl)tree->mod->ras->n_catg; */
      /*     tree->mod->ras->gamma_rr_unscaled->v[i]      = tree->mod->ras->gamma_rr->v[i]; */
      /*   } */

      /* tree->mod->s_opt->opt_alpha                 = NO; */
      /* tree->mod->ras->free_mixt_rates             = YES; */
      /* tree->mod->s_opt->first_opt_free_mixt_rates = NO; */

      /* Lk(NULL,tree); */
    }

  lk_before = tree->c_lnL;

  /*! Only skip tree traversal when data is not partitionned */
  if(tree->prev == NULL && tree->next == NULL && fast == YES)
    {
      tree->mod->s_opt->skip_tree_traversal = YES;
      tree->mod->ras->normalise_rr          = NO;

      wm = Weighted_Mean(tree->mod->ras->gamma_rr_unscaled->v,
                         tree->mod->ras->gamma_r_proba->v,
                         tree->mod->ras->n_catg);

      tree->mod->ras->free_rate_mr->v = 100.;
      For(i,2*tree->n_otu-1) tree->a_edges[i]->l->v /= (wm * tree->mod->ras->free_rate_mr->v);
    }

  /* BFGS returns negative values sometimes: need to log-transform?... to do... */
  /* int failed = NO; */
  /* BFGS(tree,tree->mod->ras->gamma_r_proba_unscaled->v,tree->mod->ras->n_catg,1.e-5,tree->mod->s_opt->min_diff_lk_local,1.e-5,NO, */
  /*      &Return_Abs_Lk, */
  /*      &Num_Derivative_Several_Param, */
  /*      &Lnsrch,&failed); */


  For(i,tree->mod->ras->n_catg-1)
    {
      Generic_Brent_Lk(&(tree->mod->ras->gamma_r_proba_unscaled->v[i]),
                       0.0,
                       100.,
                       tree->mod->s_opt->min_diff_lk_local,
                       tree->mod->s_opt->brent_it_max,
                       tree->mod->s_opt->quickdirty,
                       Wrap_Lk,NULL,tree,NULL,NO);
    }

  if(tree->mod->s_opt->skip_tree_traversal == YES && fast == YES)
    {
      tree->mod->s_opt->skip_tree_traversal = NO;
      tree->mod->ras->normalise_rr          = YES;

      wm = Weighted_Mean(tree->mod->ras->gamma_rr_unscaled->v,
                         tree->mod->ras->gamma_r_proba->v,
                         tree->mod->ras->n_catg);

      For(i,2*tree->n_otu-1) tree->a_edges[i]->l->v *= (wm * tree->mod->ras->free_rate_mr->v);
    }

  lk_after = tree->c_lnL;

  if(lk_after < lk_before - tree->mod->s_opt->min_diff_lk_global)
    {
      PhyML_Printf("\n== lk_before: %f lk_after: %f diff: %G",lk_before,lk_after,lk_before-lk_after);
      PhyML_Printf("\n== Err. in file %s at line %d\n",__FILE__,__LINE__);
      Exit("");
    }

  if(verbose) Print_Lk(tree,"[Rate class freqs.  ]");

  /* printf("\n. Before: %f after: %f",lk_before,lk_after); */
  /* For(i,tree->mod->ras->n_catg) */
  /*   { */
  /*     printf("\n+ c %2d p: %15f r: %15f up: %15f ur: %5f", */
  /*            i+1, */
  /*            tree->mod->ras->gamma_r_proba->v[i], */
  /*            tree->mod->ras->gamma_rr->v[i], */
  /*            tree->mod->ras->gamma_r_proba_unscaled->v[i], */
  /*            tree->mod->ras->gamma_rr_unscaled->v[i]); */
  /*   } */
  /* fflush(NULL); */

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Optimize_State_Freqs(t_tree *mixt_tree, int verbose)
{
  vect_dbl **freqs;
  int n_freqs;
  t_tree *tree;
  int i;
  int failed;

  Switch_Eigen(YES,mixt_tree->mod);

  freqs   = NULL;
  n_freqs = 0;
  tree    = mixt_tree;

  do
    {
      if(tree->next) tree = tree->next;

      For(i,n_freqs) if(tree->mod->e_frq->pi_unscaled == freqs[i]) break;

      if(i == n_freqs)
    {
      if(!freqs) freqs = (vect_dbl **)mCalloc(1,sizeof(vect_dbl *));
      else       freqs = (vect_dbl **)mRealloc(freqs,n_freqs+1,sizeof(vect_dbl *));
      freqs[n_freqs] = tree->mod->e_frq->pi_unscaled;
      n_freqs++;

      if((tree->mod->s_opt->opt_state_freq) && (tree->io->datatype == NT))
        {
          failed = YES;

          BFGS(mixt_tree,tree->mod->e_frq->pi_unscaled->v,tree->mod->ns,1.e-5,tree->mod->s_opt->min_diff_lk_local,1.e-5,NO,YES,
               &Return_Abs_Lk,
               &Num_Derivative_Several_Param,
               &Lnsrch,&failed);
          
          if(failed == YES)
        {
          For(i,tree->mod->ns)
            {
              Generic_Brent_Lk(&(tree->mod->e_frq->pi_unscaled->v[i]),
                       0.,100.,
                       tree->mod->s_opt->min_diff_lk_local,
                       tree->mod->s_opt->brent_it_max,
                       tree->mod->s_opt->quickdirty,
                       Wrap_Lk,NULL,mixt_tree,NULL,NO);
            }
        }

          if(verbose)
        {
          Print_Lk(mixt_tree,"[Nucleotide freqs.  ]");
        }
        }
    }

      tree = tree->next;

    }
  while(tree);

  if(freqs) Free(freqs);

  Switch_Eigen(NO,mixt_tree->mod);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Optimize_Rmat_Weights(t_tree *mixt_tree, int verbose)
{
  scalar_dbl *r_mat_weight;

  Switch_Eigen(NO,mixt_tree->mod);

  if(mixt_tree->is_mixt_tree == NO) return;

  r_mat_weight = mixt_tree->next->mod->r_mat_weight;
  
  if(mixt_tree->next->mod->s_opt->opt_rmat_weight == YES)
    {
      do
        {
          Generic_Brent_Lk(&(r_mat_weight->v),
                           0.,100.,
                           mixt_tree->mod->s_opt->min_diff_lk_local,
                           mixt_tree->mod->s_opt->brent_it_max,
                           mixt_tree->mod->s_opt->quickdirty,
                           Wrap_Lk,NULL,mixt_tree,NULL,NO);
          
          if(verbose)
            {
              Print_Lk(mixt_tree,"[Rate mat. weights  ]");
            }
          
          r_mat_weight = r_mat_weight->next;
        }
      while(r_mat_weight);
    }

  Switch_Eigen(NO,mixt_tree->mod);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Optimize_Efrq_Weights(t_tree *mixt_tree, int verbose)
{
  scalar_dbl *e_frq_weight;

  Switch_Eigen(NO,mixt_tree->mod);

  if(mixt_tree->is_mixt_tree == NO) return;

  e_frq_weight = mixt_tree->next->mod->e_frq_weight;


  if(mixt_tree->next->mod->s_opt->opt_efrq_weight == YES)
    {
      do
        {
          Generic_Brent_Lk(&(e_frq_weight->v),
                           0.,100.,
                           mixt_tree->mod->s_opt->min_diff_lk_local,
                           mixt_tree->mod->s_opt->brent_it_max,
                           mixt_tree->mod->s_opt->quickdirty,
                           Wrap_Lk,NULL,mixt_tree,NULL,NO);
          
          if(verbose)
            {
              Print_Lk(mixt_tree,"[Equ. frq. weights  ]");
            }
          
          e_frq_weight = e_frq_weight->next;
        }
      while(e_frq_weight);
    }

  Switch_Eigen(NO,mixt_tree->mod);

}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Optimize_Lambda(t_tree *mixt_tree, int verbose)
{
  scalar_dbl **lambda;
  int n_lambda;
  t_tree *tree;
  int i;

  Switch_Eigen(YES,mixt_tree->mod);

  lambda   = NULL;
  n_lambda = 0;
  tree     = mixt_tree;

  do
    {
      if(tree->next) tree = tree->next;

      For(i,n_lambda) if(tree->mod->lambda == lambda[i]) break;

      if(i == n_lambda)
    {
      if(!lambda) lambda = (scalar_dbl **)mCalloc(1,sizeof(scalar_dbl *));
      else        lambda = (scalar_dbl **)mRealloc(lambda,n_lambda+1,sizeof(scalar_dbl *));
      lambda[n_lambda] = tree->mod->lambda;
      n_lambda++;

      if(tree->mod->s_opt->opt_lambda)
        {
          Generic_Brent_Lk(&(tree->mod->lambda->v),
                   0.001,100.,
                   tree->mod->s_opt->min_diff_lk_local,
                   tree->mod->s_opt->brent_it_max,
                   tree->mod->s_opt->quickdirty,
                   Wrap_Lk,NULL,mixt_tree,NULL,NO);

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

  Switch_Eigen(NO,mixt_tree->mod);

}

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

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

