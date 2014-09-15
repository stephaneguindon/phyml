/*

PHYML :  a program that  computes maximum likelihood  phyLOGenies from
DNA or AA homoLOGous sequences 

Copyright (C) Stephane Guindon. Oct 2003 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/

/*

The code below is an implementation of the building tree algorithm
described in "BIONJ: an improved version of the NJ algorithm based 
on a simple model of sequence data." (1997) O. Gascuel. Mol Biol Evol. 
14:685-95.  

*/

#include "bionj.h"


void Bionj(matrix *mat)
{
  int x,y,i;
  phydbl vxy,lx,ly,lamda,score;

  Clean_Tree_Connections(mat->tree);  
  For(i,mat->tree->n_otu) mat->tip_node[i] = mat->tree->a_nodes[i];
  mat->tree->num_curr_branch_available = 0;

  while(mat->r > 3)
    {
      x = y =  0;
      vxy   = .0;
      score = .0;
      Compute_Sx(mat);
      Best_Pair(mat,&x,&y,&score);
      vxy=Variance(mat,x,y);
      lx=Br_Length(mat,x,y);    
      ly=Br_Length(mat,y,x);
      lamda=Lamda(mat,x,y,vxy); 
      Update_Mat(mat,x,y,lx,ly,vxy,lamda);
      Update_Tree(mat,x,y,lx,ly,score);      
    }
  Finish(mat);
}
  
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Finish(matrix *mat)
{
  phydbl dxy,dxz,dyz;
  int x,y,z;
  t_node *nx,*ny,*nz,*new;
  int i;

  dxy = dxz = dyz = -1.;
  x = y = z = -1;

  For(i,mat->n_otu)
    {
      if(mat->on_off[i])
	{
	  if(x < 0) x=i;
	  else if(y < 0) y = i;
	  else if(z < 0) z = i;
	}
    }

  dxy = Dist(mat,x,y);
  dxz = Dist(mat,x,z);
  dyz = Dist(mat,y,z);

  nx = mat->tip_node[x];
  ny = mat->tip_node[y];
  nz = mat->tip_node[z];

  new = mat->tree->a_nodes[mat->curr_int];
  new->num = mat->curr_int;
  new->v[0] = nx;
  new->v[1] = ny;
  new->v[2] = nz;

  nx->v[0] = new;
  ny->v[0] = new;
  nz->v[0] = new;

  
  Connect_One_Edge_To_Two_Nodes(new,nx,mat->tree->a_edges[mat->tree->num_curr_branch_available],mat->tree);
  Connect_One_Edge_To_Two_Nodes(new,ny,mat->tree->a_edges[mat->tree->num_curr_branch_available],mat->tree);
  Connect_One_Edge_To_Two_Nodes(new,nz,mat->tree->a_edges[mat->tree->num_curr_branch_available],mat->tree);
  
 
  nx->b[0]->l->v = .5*(dxy-dyz+dxz);
  ny->b[0]->l->v = .5*(dyz-dxz+dxy);
  nz->b[0]->l->v = .5*(dxz-dxy+dyz);
   
  new->b[0]->l->v = nx->b[0]->l->v;
  new->b[1]->l->v = ny->b[0]->l->v;
  new->b[2]->l->v = nz->b[0]->l->v;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Update_Mat(matrix *mat, int x, int y, phydbl lx, phydbl ly, phydbl vxy, phydbl lamda)
{
  int i;
  int a,b;
  
  a = b = -1;
  For(i,mat->n_otu)
    {
      if((mat->on_off[i]) && (i != x) && (i != y))
	{
	  if(x > i)
	    {
	      a=x;
	      b=i;
	    }
	  else
	    {
	      a=i;
	      b=x;
	    }
	  mat->dist[a][b]=Dist_Red(mat,x,lx,y,ly,i,lamda);
	  mat->dist[b][a]=Var_Red(mat,x,y,i,lamda,vxy);
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Update_Tree(matrix *mat, int x, int y, phydbl lx, phydbl ly, phydbl score)
{
  t_node *new, *nx, *ny;

  nx            = mat->tip_node[x];
  ny            = mat->tip_node[y];
  new           = mat->tree->a_nodes[mat->curr_int];
  nx->v[0]      = new;
  ny->v[0]      = new;
  new->v[1]     = nx;
  new->v[2]     = ny;
  new->num      = mat->curr_int;

  Connect_One_Edge_To_Two_Nodes(new,nx,mat->tree->a_edges[mat->tree->num_curr_branch_available],mat->tree);
  Connect_One_Edge_To_Two_Nodes(new,ny,mat->tree->a_edges[mat->tree->num_curr_branch_available],mat->tree);

  nx->b[0]->l->v   = lx;
  ny->b[0]->l->v   = ly;
  
  new->b[1]->l->v  = lx;
  new->b[2]->l->v  = ly;
  new->score[0] = score;

  nx->l[0]      = lx;
  ny->l[0]      = ly;
  
  new->l[1]     = lx;
  new->l[2]     = ly;
 
  mat->tip_node[x] = new;
  mat->on_off[y] = 0;
  mat->curr_int++;
  mat->r--;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Best_Pair(matrix *mat, int *x, int *y,phydbl *score)
{
  int i,j/* ,n_ties */;
  phydbl Qmin,Qmin2;
  phydbl *t_Qij;
/*   int *ties; */

  t_Qij = (phydbl *)mCalloc(mat->n_otu * mat->n_otu,sizeof(phydbl ));
/*   ties  = (int *)mCalloc(mat->n_otu * mat->n_otu,sizeof(int )); */

  Qmin = 1.e+10;
  
  For(i,mat->n_otu)
    {
      if(mat->on_off[i])
	{
	  for(j=0;j<i;j++)
	    {
	      if(mat->on_off[j])
		{
		  t_Qij[mat->n_otu*i+j] = Q_Agglo(mat,i,j);

		  if(t_Qij[mat->n_otu*i+j] < Qmin - 1.E-05)
		    {
		      *x = i;
		      *y = j;
		      Qmin = t_Qij[mat->n_otu*i+j];
		    }
		}
	    }
	}
    }
  
/*   n_ties = 0; */
/*   For(i,mat->n_otu) */
/*     { */
/*       if(mat->on_off[i]) */
/* 	{ */
/* 	  for(j=0;j<i;j++) */
/* 	    { */
/* 	      if(mat->on_off[j]) */
/* 		{ */
/* 		  if((t_Qij[mat->n_otu*i+j] < Qmin + 1.E-05) && (t_Qij[mat->n_otu*i+j] > Qmin - 1.E-05)) */
/* 		    { */
/* 		      ties[n_ties] = mat->n_otu*i+j; */
/* 		      n_ties++; */
/* 		    } */
/* 		} */
/* 	    } */
/* 	} */
/*     } */
   
/*   if(!n_ties) */
/*     { */
/*       PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__); */
/*       Warn_And_Exit(""); */
/*     } */


/*   /\* Useful if some pairwise distances are null *\/ */
/*   if(n_ties > 1) */
/*     { */
/*       int cand; */
/*       *x = *y = -1; */
/*       cand = (int)RINT(rand()/(phydbl)(RAND_MAX) * (n_ties-1)); */
/*       *x = (int)(ties[cand] / mat->n_otu); */
/*       *y = (int)(ties[cand] % mat->n_otu); */
/*     } */



  Qmin2 = 1e+10;

  For(i,mat->n_otu)
    {
      if((i != *y) && (i != *x) && (t_Qij[mat->n_otu*(*x)+i] < Qmin2)) Qmin2 = t_Qij[mat->n_otu*(*x)+i];
    }

  For(i,mat->n_otu)
    {
      if((i != *y) && (i != *x) && (t_Qij[mat->n_otu*i+(*y)] < Qmin2)) Qmin2 = t_Qij[mat->n_otu*i+(*y)];
    }

  *score = FABS(Qmin2 - Qmin)/FABS(Qmin);

  Free(t_Qij);
/*   Free(ties); */
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Compute_Sx(matrix *mat)
{
  int i,j;
  
  For(i,mat->n_otu)
    {
      mat->dist[i][i] = .0;
      if(mat->on_off[i])
	{
	  For(j,mat->n_otu)
	    {
	      if((i != j) && (mat->on_off[j]))
		{
		  mat->dist[i][i] += Dist(mat,i,j);
		}
	    }
	}
    }
}
	      
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Sum_S(matrix *mat, int i)
{
  return mat->dist[i][i];
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Dist(matrix *mat, int x, int y)
{
    if(x > y)
      return(mat->dist[x][y]);
    else
      return(mat->dist[y][x]);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Variance(matrix *mat, int x, int y)
{
    if(x > y)
      {
	return(mat->dist[y][x]);
      }
    else
      {
	return(mat->dist[x][y]);
      }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Br_Length(matrix *mat, int x, int y)
{
    return .5*(Dist(mat,x,y)+
	      (Sum_S(mat,x)-Sum_S(mat,y))/(phydbl)(mat->r-2.)); 
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Dist_Red(matrix *mat, int x, phydbl lx, int y, phydbl ly, int i, phydbl lamda)
{
  phydbl Dui;
  Dui=lamda*(Dist(mat,x,i)-lx)
     +(1.-lamda)*(Dist(mat,y,i)-ly);
  return(Dui);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Var_Red(matrix *mat, int x, int y, int i, phydbl lamda, phydbl vxy)
{
  phydbl Vui;
  Vui=lamda*(Variance(mat,x,i))
     +(1.-lamda)*(Variance(mat,y,i))
    -lamda*(1.-lamda)*vxy;
  return(Vui);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Lamda(matrix *mat, int x, int y, phydbl vxy)
{
    phydbl lamda=0.0;
    int i;
    
    if(mat->method == 0) /* NJ (Saitou & Nei, 1987) */
      lamda = 0.5;
    else /* BioNJ (Gascuel, 1997) */
      {
	if(vxy < SMALL && vxy > -SMALL)
	  lamda=0.5;
	else
	  {
	    For(i,mat->n_otu)
	      {
		if((x != i) && (y != i) && (mat->on_off[i]))
		  lamda = lamda + Variance(mat,y,i) - Variance(mat,x,i);
	      }
	    lamda = 0.5 + lamda/(2.*(mat->r-2)*vxy);
	  }
	
	if(lamda > 1.0)
	  lamda = 0.5;/*1.0;*/
	else if(lamda < 0.0)
	  lamda = 0.5;/*0.0;*/
      }

    return(lamda);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Q_Agglo(matrix *mat, int x, int y)
{
  phydbl Qxy;

  Qxy = .0;
  Qxy=(mat->r-2.)*Dist(mat,x,y)-Sum_S(mat,x)-Sum_S(mat,y); 
  return(Qxy);                       
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Bionj_Br_Length(matrix *mat)
{
  int x;

  x = Bionj_Br_Length_Post(mat->tree->a_nodes[0],
			   mat->tree->a_nodes[0]->v[0],
			   mat);
  mat->tree->a_nodes[0]->b[0]->l->v = Dist(mat,0,x);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


int Bionj_Br_Length_Post(t_node *a, t_node *d, matrix *mat)
{
  int i;

  if(d->tax)
    {
      return d->num;
    }
  else
    {
      int d_v1, d_v2;
      phydbl lx, ly, vxy,lamda;
      int x,y;

      d_v1 = d_v2 = -1;
      For(i,3)
	if(d->v[i] != a) {(d_v1 < 0)?(d_v1 = i):(d_v2 = i);}
      

      x = Bionj_Br_Length_Post(d,d->v[d_v1],mat);
      y = Bionj_Br_Length_Post(d,d->v[d_v2],mat);

      vxy = .0;
      Compute_Sx(mat);
      vxy=Variance(mat,(x),(y));
      lx=Br_Length(mat,(x),(y));    
      ly=Br_Length(mat,(y),(x));
      lamda=Lamda(mat,(x),(y),vxy); 
      Update_Mat(mat,(x),(y),lx,ly,vxy,lamda);

      d->b[d_v1]->l->v = lx;
      d->b[d_v2]->l->v = ly;
      
      mat->on_off[y] = 0;
      mat->r--;

      return x;
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////







