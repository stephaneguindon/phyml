/*

PHYML :  a program that  computes maximum likelihood  phyLOGenies from
DNA or AA homoLOGous sequences

Copyright (C) Stephane Guindon. Oct 2003 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/

#include "draw.h"


void DR_Draw_Tree(char *file_name, t_tree *tree)
{
  FILE *ps_tree;
  
  ps_tree = (FILE *)fopen(file_name,"w");  
  DR_Print_Postscript_Header(1,ps_tree);
  tree->ps_tree = DR_Make_Tdraw_Struct(tree);
  DR_Init_Tdraw_Struct(tree->ps_tree);
  DR_Get_Tree_Box_Width(tree->ps_tree,tree);
  Dist_To_Root(tree->n_root,tree);
  tree->ps_tree->max_dist_to_root = DR_Get_Max_Dist_To_Root(tree);
  DR_Get_X_Coord(NO,tree->ps_tree,tree);
  DR_Get_Y_Coord(NO,tree->ps_tree,tree);
  DR_Print_Tree_Postscript(1,NO,ps_tree,tree);
  DR_Print_Postscript_EOF(ps_tree);
  fclose(ps_tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void DR_Get_Tree_Coord(t_tree *tree)
{
  DR_Init_Tdraw_Struct(tree->ps_tree);
  DR_Get_Tree_Box_Width(tree->ps_tree,tree);
  if(!tree->n_root) 
    {
      PhyML_Printf("\n. Adding root before rendering the tree.");
      Add_Root(tree->a_edges[0],tree);
    }
  Dist_To_Root(tree->n_root,tree);
  tree->ps_tree->max_dist_to_root = DR_Get_Max_Dist_To_Root(tree);
  DR_Get_X_Coord(NO,tree->ps_tree,tree);
  DR_Get_Y_Coord(NO,tree->ps_tree,tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void DR_Print_Postscript_Header(int n_pages, FILE *fp)
{
  if(!fp)
    {
      PhyML_Printf("\n== Failed to open the postscript file.");
      PhyML_Printf("\n== Did you forget the '--ps' option ?.");
      Exit("\n");
    }

  PhyML_Fprintf(fp,"%%!PS-Adobe-3.0\n");
  PhyML_Fprintf(fp,"%%%%BoundingBox: 0 0 595.28 841.89\n");
  PhyML_Fprintf(fp,"%%%%DocumentFonts: Times-Roman Times-Roman\n");
  PhyML_Fprintf(fp,"%%%%Creator: Stephane Guindon\n");
  PhyML_Fprintf(fp,"%%%%Title: tree\n");
  PhyML_Fprintf(fp,"%%%%EndComments\n");
  PhyML_Fprintf(fp,"%%%%Pages: %d\n",n_pages);

  PhyML_Fprintf(fp,"/lt {lineto} bind def\n");
  PhyML_Fprintf(fp,"/mt {moveto} bind def\n");
  PhyML_Fprintf(fp,"/sc {setrgbcolor} bind def\n");
  PhyML_Fprintf(fp,"/ct {curveto} bind def\n");
  PhyML_Fprintf(fp,"/np {newpath} bind def\n");
  PhyML_Fprintf(fp,"/cp {closepath} bind def\n");
  PhyML_Fprintf(fp,"/gs {gsave} bind def\n");
  PhyML_Fprintf(fp,"/gr {grestore} bind def\n");
  
  PhyML_Fprintf(fp,"/Times-Roman findfont\n");
  PhyML_Fprintf(fp,"12 scalefont\n");
  PhyML_Fprintf(fp,"setfont\n");

  PhyML_Fprintf(fp,"/clipbox\n");
  PhyML_Fprintf(fp,"{\n");
  PhyML_Fprintf(fp,"newpath\n");
  PhyML_Fprintf(fp,"20 20 mt\n");
  PhyML_Fprintf(fp,"580 20 lt\n");
  PhyML_Fprintf(fp,"580 820 lt\n");
  PhyML_Fprintf(fp,"20 820 lt\n");
  PhyML_Fprintf(fp,"20 20 lt\n");
  PhyML_Fprintf(fp,"closepath\n");
  PhyML_Fprintf(fp,"clip\n");
  PhyML_Fprintf(fp,"} bind def\n");


  /* PhyML_Fprintf(fp,"gs\n"); */
  /* PhyML_Fprintf(fp,"newpath\n"); */
  /* PhyML_Fprintf(fp,"20 20 mt\n"); */
  /* PhyML_Fprintf(fp,"580 20 lt\n"); */
  /* PhyML_Fprintf(fp,"580 820 lt\n"); */
  /* PhyML_Fprintf(fp,"20 820 lt\n"); */
  /* PhyML_Fprintf(fp,"20 20 lt\n"); */
  /* PhyML_Fprintf(fp,"closepath\n"); */
  /* PhyML_Fprintf(fp,"stroke\n"); */
  /* PhyML_Fprintf(fp,"gr\n"); */
  /* PhyML_Fprintf(fp,"0 0 0 sc\n"); */


}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void DR_Print_Postscript_EOF(FILE *fp)
{
  PhyML_Fprintf(fp,"%%%%Trailer\n");
  PhyML_Fprintf(fp,"%%%%EOF\n");
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void DR_Print_Tree_Postscript(int page_num, int render_name, FILE *fp, t_tree *tree)
{
  tdraw *draw;
  t_node *n_root;
  
  draw = tree->ps_tree;
/*   DR_Get_Tree_Coord(tree); */
  n_root = tree->n_root;

/*   PhyML_Fprintf(fp,"%%%%Page: %d %d\n",page_num,page_num); */
  /* PhyML_Fprintf(fp,"0.001 setlinewidth\n"); */
/*   PhyML_Fprintf(fp,"0.5 0.5 0.4 sc\n"); */
  /* PhyML_Fprintf(fp,"0 0 0 sc\n"); */
/*   PhyML_Fprintf(fp,"clipbox\n"); */
/*   PhyML_Fprintf(fp,"stroke\n"); */
  PhyML_Fprintf(fp,"20 20 translate\n");
  PhyML_Fprintf(fp,"newpath\n");


  draw->ycoord[n_root->num] = (draw->ycoord[n_root->v[2]->num] + draw->ycoord[n_root->v[1]->num])/2. + 20; 
  draw->xcoord[n_root->num] = 0.0; 
  DR_Print_Tree_Postscript_Pre(n_root,n_root->v[2],n_root->b[2],render_name,fp,draw,tree);
  DR_Print_Tree_Postscript_Pre(n_root,n_root->v[1],n_root->b[1],render_name,fp,draw,tree);


  PhyML_Fprintf(fp,"closepath\n");
  PhyML_Fprintf(fp,"0 0 translate\n");
  PhyML_Fprintf(fp,"stroke\n");
  PhyML_Fprintf(fp,"showpage\n");
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void DR_Print_Tree_Postscript_Pre(t_node *a, t_node *d, t_edge *b, int render_name, FILE *fp, tdraw *w, t_tree *tree)
{
  int i;
  phydbl R, G, B;
  
  R = G = B = 0.0;

  PhyML_Fprintf(fp,"gs\n");
  
  PhyML_Fprintf(fp,"%.1f %.1f mt\n",w->xcoord[a->num],w->ycoord[a->num]);

/*   PhyML_Fprintf(fp,"%.1f %.1f lt\n",w->xcoord[a->num],w->ycoord[d->num]); */
/*   PhyML_Fprintf(fp,"%.1f %.1f lt\n",w->xcoord[d->num],w->ycoord[d->num]); */

  phydbl min,max,step,val;
  
  min = 0.0;
  max = 5.;
  
  step = (max-min)/13.;
  
  /* val = tree->rates->mean_r[d->num] / (phydbl)(tree->mcmc->run/tree->mcmc->sample_interval+1.); */
  /* val = tree->rates->mean_r[d->num]; */
  /* val = tree->rates->has_survived[d->num]; */
  /* if(val > 0.5) {R=1.; G=.0; B=0.;} */
  /* val = tree->geo->ldscape[tree->geo->idx_loc[d->num]*tree->geo->n_dim+0] + 2.5; */
  val = tree->geo->coord_loc[tree->geo->idx_loc[d->num]]->lonlat[0] + 2.5;
  /* val = 0.; */
  
  if(val <= min+1.*step)
    {R=.0; G=1.; B=1.;}
  else if(val > min+1.*step && val <= min+2.*step)
    {R=.0; G=1.; B=.8;}
  else if(val > min+2.*step && val <= min+3.*step)
    {R=.0; G=1.; B=.5;}
  else if(val > min+3.*step && val <= min+4.*step)
    {R=.0; G=1.; B=.3;}
  else if(val > min+4.*step && val <= min+5.*step)
    {R=.0; G=1.; B=.0;}
  else if(val > min+5.*step && val <= min+6.*step)
    {R=.25; G=1.; B=0.;}
  else if(val > min+6.*step && val <= min+7.*step)
    {R=.5; G=1.; B=0.;}
  else if(val > min+7.*step && val <= min+8.*step)
    {R=.75; G=1.; B=.0;}
  else if(val > min+8.*step && val <= min+9.*step)
    {R=1.; G=1.; B=.0;}
  else if(val > min+9.*step && val <= min+10.*step)
    {R=1.; G=.75; B=.0;}
  else if(val > min+10.*step && val <= min+11.*step)
    {R=1.; G=.5; B=.0;}
  else if(val > min+11.*step && val <= min+12.*step)
    {R=1.; G=.25; B=.0;}
  else if(val > min+12.*step)
    {R=1.; G=.0; B=0.;}


  /* R = 0.; G = 0.; B = 0.; */

  PhyML_Fprintf(fp,"2 setlinewidth\n");
  /* PhyML_Fprintf(fp,"%.1f %.1f lt\n",w->xcoord[a->num],w->ycoord[d->num]); */
  /* PhyML_Fprintf(fp,"%.1f %.1f lt\n",w->xcoord[d->num],w->ycoord[d->num]); */

  phydbl xa = w->xcoord[a->num];
  phydbl xd = MIN(w->xcoord[a->num] + 5,w->xcoord[d->num]);
  phydbl ya = w->ycoord[a->num];
  phydbl yd = w->ycoord[d->num];
  
  PhyML_Fprintf(fp,"%.1f %.1f %.1f %.1f %.1f %.1f ct\n",
                xa + (xd-xa)/2.,(ya+yd)/2.,
                xd - (xd-xa)/5.,yd,
                xd,yd);
  
  PhyML_Fprintf(fp,"%.1f %.1f lt\n",w->xcoord[d->num],w->ycoord[d->num]);

  /* PhyML_Fprintf(fp,"%.1f %.1f %.1f %.1f %.1f %.1f ct\n", */
  /* 		w->xcoord[a->num], */
  /* 		w->ycoord[d->num], */
  /* 		w->xcoord[d->num], */
  /* 		w->ycoord[d->num], */
  /* 		w->xcoord[d->num], */
  /* 		w->ycoord[d->num]); */


  if(tree->rates && tree->rates->has_survived[d->num] == YES)
    {
      PhyML_Fprintf(fp," /Helvetica findfont 16 scalefont\n");
      PhyML_Fprintf(fp,"setfont\n");
      PhyML_Fprintf(fp,"%.1f %.1f mt\n",w->xcoord[d->num]-5,w->ycoord[d->num]);
      PhyML_Fprintf(fp,"0 0 0 sc\n");
      PhyML_Fprintf(fp,"(*) show \n");      
    }

  PhyML_Fprintf(fp,"%f %f %f sc\n",R,G,B);

  if(d->tax) 
    {
      PhyML_Fprintf(fp,"stroke\n");
      PhyML_Fprintf(fp,"0 setgray\n");
      PhyML_Fprintf(fp,"2 setlinewidth\n");
      PhyML_Fprintf(fp,"np %.1f %.1f 1 0 360 arc cp\n",w->xcoord[d->num],w->ycoord[d->num]);
      PhyML_Fprintf(fp,"%.1f %.1f %.1f sc fill\n",R,G,B);
/*       PhyML_Fprintf(fp,"%f setgray fill\n",greylevel); */
      PhyML_Fprintf(fp,"0 0 0 sc\n");


      PhyML_Fprintf(fp," /Helvetica findfont 10 scalefont\n");
      PhyML_Fprintf(fp,"setfont\n");
      PhyML_Fprintf(fp,"%.1f %.1f mt\n",w->xcoord[d->num]+2,w->ycoord[d->num]-6);
      /* PhyML_Fprintf(fp,"(%d) show \n",d->num); */
      /* PhyML_Fprintf(fp,"(%s) show \n",d->name); */
      PhyML_Fprintf(fp," /Helvetica findfont 14 scalefont\n");
      PhyML_Fprintf(fp,"setfont\n");


      PhyML_Fprintf(fp,"%.1f %.1f mt\n",w->xcoord[d->num] - (w->xcoord[d->num] - w->xcoord[a->num])/2.,w->ycoord[d->num]);
      PhyML_Fprintf(fp," /Helvetica findfont 10 scalefont\n");
      PhyML_Fprintf(fp,"setfont\n");

#if (defined GEO)
      /* PhyML_Fprintf(fp,"([%4.4f,%4.4f]) show \n", */
      /*               tree->geo->ldscape[tree->geo->loc[d->num]*tree->geo->n_dim+0], */
      /*               tree->geo->ldscape[tree->geo->loc[d->num]*tree->geo->n_dim+1]); */
#endif

     
      /* PhyML_Fprintf(fp,"%.1f %.1f mt\n",w->xcoord[d->num]+5,w->ycoord[d->num]); */
      /* PhyML_Fprintf(fp,"(%.10s) show \n",d->name); */

      
/*       if(render_name) */
/* 	{ */
/* 	  if(tree->io->long_tax_names) */
/* 	    PhyML_Fprintf(fp,"(%s) show \n",tree->io->long_tax_names[d->num]); */
/* 	  else */
/* 	    PhyML_Fprintf(fp,"(%s) show \n",d->name); */
/* 	} */

      PhyML_Fprintf(fp,"stroke\n");
      PhyML_Fprintf(fp,"gr\n");
      PhyML_Fprintf(fp,"0 0 0 sc\n");
      return;
    }
  else
    {
      PhyML_Fprintf(fp,"stroke\n");
      PhyML_Fprintf(fp,"0 setgray\n");
      PhyML_Fprintf(fp,"2 setlinewidth\n");
      PhyML_Fprintf(fp,"np %.1f %.1f 1 0 360 arc cp\n",w->xcoord[d->num],w->ycoord[d->num]);
      PhyML_Fprintf(fp,"%.1f %.1f %.1f sc fill\n",R,G,B);
/*       PhyML_Fprintf(fp,"%f setgray fill\n",greylevel); */
      PhyML_Fprintf(fp,"0 0 0 sc\n");

      PhyML_Fprintf(fp," /Helvetica findfont 10 scalefont\n");
      PhyML_Fprintf(fp,"setfont\n");
      PhyML_Fprintf(fp,"%.1f %.1f mt\n",w->xcoord[d->num]+2,w->ycoord[d->num]);

      /* PhyML_Fprintf(fp,"(%d) show \n",b->num); */

      PhyML_Fprintf(fp,"%.1f %.1f mt\n",w->xcoord[d->num],w->ycoord[d->num]);
      PhyML_Fprintf(fp," /Helvetica findfont 14 scalefont\n");
      PhyML_Fprintf(fp,"setfont\n");


      PhyML_Fprintf(fp,"%.1f %.1f mt\n",w->xcoord[d->num] - (w->xcoord[d->num] - w->xcoord[a->num])/2.,w->ycoord[d->num]);
      PhyML_Fprintf(fp," /Helvetica findfont 10 scalefont\n");
      PhyML_Fprintf(fp,"setfont\n");

#if (defined GEO)
      /* PhyML_Fprintf(fp,"([%4.4f,%4.4f]) show \n", */
      /*               tree->geo->ldscape[tree->geo->loc[d->num]*tree->geo->n_dim+0], */
      /*               tree->geo->ldscape[tree->geo->loc[d->num]*tree->geo->n_dim+1]); */
#endif      


      PhyML_Fprintf(fp,"stroke\n");
      PhyML_Fprintf(fp,"gr\n");
      PhyML_Fprintf(fp,"0 0 0 sc\n");
      For(i,3)
	if(d->v[i] != a && d->b[i] != tree->e_root) DR_Print_Tree_Postscript_Pre(d,d->v[i],d->b[i],render_name,fp,w,tree);
    }




  return;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void DR_Get_X_Coord_Pre(t_node *a, t_node *d, t_edge *b, tdraw *w, int fixed_tips, t_tree *tree)
{
  int i;

  if(!(d->tax && fixed_tips == YES)) w->xcoord[d->num] =  d->dist_to_root * (phydbl)w->tree_box_width/w->max_dist_to_root;

  if(d->tax) return;
  else
    {
      For(i,3)
	if((d->v[i] != a) && (d->b[i] != tree->e_root)) 
	  DR_Get_X_Coord_Pre(d,d->v[i],d->b[i],w,fixed_tips,tree);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void DR_Get_X_Coord(int fixed_tips, tdraw *w, t_tree *tree)
{
  if(!(tree->n_root->v[2]->tax && fixed_tips == YES)) w->xcoord[tree->n_root->v[2]->num] = tree->n_root->v[2]->dist_to_root * (phydbl)w->tree_box_width/w->max_dist_to_root;
  if(!(tree->n_root->v[1]->tax && fixed_tips == YES)) w->xcoord[tree->n_root->v[1]->num] = tree->n_root->v[1]->dist_to_root * (phydbl)w->tree_box_width/w->max_dist_to_root;
  DR_Get_X_Coord_Pre(tree->n_root,tree->n_root->v[2],NULL,w,fixed_tips,tree);
  DR_Get_X_Coord_Pre(tree->n_root,tree->n_root->v[1],NULL,w,fixed_tips,tree);
  w->xcoord[tree->n_root->num] = 0;
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void DR_Get_Y_Coord(int fixed_tips, tdraw *w, t_tree *tree)
{
  int next_y_slot;
  next_y_slot = 0;
  DR_Get_Y_Coord_Post(tree->n_root,tree->n_root->v[2],NULL,&next_y_slot,fixed_tips,w,tree);
  DR_Get_Y_Coord_Post(tree->n_root,tree->n_root->v[1],NULL,&next_y_slot,fixed_tips,w,tree);
  w->ycoord[tree->n_root->num] = (int)((w->ycoord[tree->n_root->v[2]->num] + w->ycoord[tree->n_root->v[2]->num]) / 2.) + 20;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void DR_Get_Y_Coord_Post(t_node *a, t_node *d, t_edge *b, int *next_y_slot, int fixed_tips, tdraw *w, t_tree *tree)
{
  int i;

  if(d->tax) 
    {
      if(!fixed_tips)
	{
/* 	  w->ycoord[d->num] = *next_y_slot + (int)(w->page_height / (2.*tree->n_otu)); */
	  w->ycoord[d->num] = *next_y_slot + 20;
	  (*next_y_slot) += (int)(w->page_height / (tree->n_otu-1));
          printf("\n. %s %f",d->name,w->ycoord[d->num]);
	}
    }
  else
    {
      int d1, d2;

      d1 = d2 = -1;
      For(i,3)
	{
	  if(d->v[i] != a && d->b[i] != tree->e_root)
	    {
	      DR_Get_Y_Coord_Post(d,d->v[i],d->b[i],next_y_slot,fixed_tips,w,tree);
	      if(d1<0) d1 = i;
	      else     d2 = i;
	    }
	}
      w->ycoord[d->num] = (w->ycoord[d->v[d1]->num] + w->ycoord[d->v[d2]->num])/2.; 
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


tdraw *DR_Make_Tdraw_Struct(t_tree *tree)
{
  tdraw *w;

  w = (tdraw *)mCalloc(1,sizeof(tdraw));
  w->xcoord = (phydbl *)mCalloc(2*tree->n_otu-1,sizeof(phydbl));
  w->ycoord = (phydbl *)mCalloc(2*tree->n_otu-1,sizeof(phydbl));
  w->xcoord_s = (phydbl *)mCalloc(2*tree->n_otu-1,sizeof(phydbl));
  w->ycoord_s = (phydbl *)mCalloc(2*tree->n_otu-1,sizeof(phydbl));
  w->cdf_mat  = (int *)mCalloc((2*tree->n_otu-2)*(2*tree->n_otu-2),sizeof(int));
  w->cdf_mat_x  = (phydbl *)mCalloc(2*tree->n_otu-1,sizeof(phydbl));
  w->cdf_mat_y  = (phydbl *)mCalloc(2*tree->n_otu-1,sizeof(phydbl));

  return w;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void DR_Init_Tdraw_Struct(tdraw *w)
{
  w->page_width  = 580-20;
  w->page_height = 820-20;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void DR_Get_Tree_Box_Width(tdraw *w, t_tree *tree)
{
  int i;
  int max_name_len, curr_len;

  max_name_len = curr_len = 0;
  For(i,tree->n_otu)
    {
      curr_len = (int)strlen(tree->a_nodes[i]->name);
      if(curr_len > max_name_len) max_name_len = curr_len;
    }

  w->tree_box_width = w->page_width - max_name_len * 8.66667;
  /* w->tree_box_width = w->page_width - max_name_len * 10.; */
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl DR_Get_Max_Dist_To_Root(t_tree *tree)
{
  phydbl mx;
  int i;

  mx = .0;
  For(i,tree->n_otu)
    {
      if(tree->a_nodes[i]->dist_to_root > mx)
	{
	  mx = tree->a_nodes[i]->dist_to_root;
	}
    }

  return mx;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void DR_Get_Tree_Coord_Scaled(tdraw *w, t_tree *tree)
{
  int i;
  int max_x,min_x;
  int max_y,min_y;
 
  max_x = -INT_MAX;
  min_x =  INT_MAX;

  For(i,2*tree->n_otu-1)
    {
      if(w->xcoord[i] > max_x) max_x = w->xcoord[i];
      if(w->xcoord[i] < min_x) min_x = w->xcoord[i];
    }

  max_y = -INT_MAX;
  min_y =  INT_MAX;

  For(i,2*tree->n_otu-1)
    {
      if(w->ycoord[i] > max_y) max_y = w->ycoord[i];
      if(w->ycoord[i] < min_y) min_y = w->ycoord[i];
    }


  For(i,2*tree->n_otu-1)
    {
      w->xcoord_s[i] = (phydbl)(w->xcoord[i] - min_x) / (max_x - min_x);
      w->ycoord_s[i] = (phydbl)(w->ycoord[i] - min_y) / (max_y - min_y);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void DR_Get_Cdf_Mat(t_tree *tree)
{
  int i,j,k;
  phydbl min_x,max_x,y;
  phydbl x_mat,y_mat;
  t_node *d, *a;
  phydbl eps;

  eps = 1.E-6;

  For(i,2*tree->n_otu-1)  tree->ps_tree->cdf_mat_x[i] = tree->ps_tree->xcoord_s[i];
  For(i,2*tree->n_otu-1)  tree->ps_tree->cdf_mat_y[i] = tree->ps_tree->ycoord_s[i];

  Qksort(tree->ps_tree->cdf_mat_x,NULL,0,2*tree->n_otu-2);
  Qksort(tree->ps_tree->cdf_mat_y,NULL,0,2*tree->n_otu-2);

  For(i,2*tree->n_otu-2) /* x coordinates */
    {
      For(j,2*tree->n_otu-2) /* y coordinates */
	{
	  For(k,2*tree->n_otu-2) /* all nodes in the tree */
	    {
	      d = tree->a_nodes[k];
	      a = tree->a_nodes[k]->anc;

	      min_x = tree->ps_tree->xcoord_s[a->num];
	      max_x = tree->ps_tree->xcoord_s[d->num];

	      y = tree->ps_tree->ycoord_s[d->num];

	      x_mat = tree->ps_tree->cdf_mat_x[i];
	      y_mat = tree->ps_tree->cdf_mat_y[j];
	
/* 	      printf("\n. x_mat=%.1f ymat=%.1f min=%.1f max=%.1f y=%.1f", */
/* 		     x_mat,y_mat,min_x,max_x,y); */
	      
	      if((min_x < x_mat + eps) && (max_x > x_mat) && (y > y_mat))
		{
		  tree->ps_tree->cdf_mat[j*(2*tree->n_otu-2)+i] += 1;
/* 		  PhyML_Printf("\n. Add 1 to [%.1f,%.1f]", */
/* 			       tree->ps_tree->cdf_mat_x[i], */
/* 			       tree->ps_tree->cdf_mat_y[j]); */
		}
	    }
	}
    }
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

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

