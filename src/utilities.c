/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include "utilities.h"
#include "assert.h"

#ifdef BEAGLE
#include "beagle_utils.h"
#endif

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl String_To_Dbl(char *string)
{
  phydbl buff;
  char *endptr;

  if(!string)
    {
      PhyML_Fprintf(stderr,"\n. String object empty.");
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);    
    }

  buff = strtod(string,&endptr);

  if(string == endptr || errno == ERANGE)
    {
      PhyML_Printf("\n. Error in translating string '%s' to double.",string);
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);    
    }
  return buff;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int String_To_Int(char *string)
{
  int buff;
  char *endptr;

  if(!string)
    {
      PhyML_Fprintf(stderr,"\n. String object empty.");
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);    
    }

  buff = (int)strtol(string,&endptr, 10);

  if(string == endptr || errno == ERANGE)
    {
      PhyML_Fprintf(stderr,"\n. Error in translating string '%s' to integer.",string);
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);    
    }

  return buff;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Unroot_Tree(char **subtrees)
{
  char **tmp_sub;
  int degree,i,j;

  PhyML_Printf("\n. Removing the root...\n");

  tmp_sub = Sub_Trees(subtrees[0],&degree);
  if(degree >= 2)
    {
      strcpy(subtrees[2],subtrees[1]);
      Clean_Multifurcation(tmp_sub,degree,2);
      for(j=0;j<2;j++) strcpy(subtrees[j],tmp_sub[j]);
    }
  else
    {
      tmp_sub = Sub_Trees(subtrees[1],&degree);
      strcpy(subtrees[2],subtrees[0]);
      Clean_Multifurcation(tmp_sub,degree,2);
      for(j=0;j<2;j++) strcpy(subtrees[j],tmp_sub[j]);
    }

  for(i=0;i<degree;i++) Free(tmp_sub[i]);
  Free(tmp_sub);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Set_Edge_Dirs(t_edge *b, t_node *a, t_node *d, t_tree *tree)
{
  int i;

  if(a == b->rght)
    {
      PhyML_Fprintf(stderr,"\n. a->num = %d ; d->num = %d",a->num,d->num);
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);    
    }

  if(d == b->left)
    {
      PhyML_Fprintf(stderr,"\n. a->num = %d ; d->num = %d",a->num,d->num);
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);    
    }

  b->l_r = b->r_l = -1;
  for(i=0;i<3;i++)
    {
      /* if((a->v[i]) && ((a->v[i] == d) || (e_root && a->b[i] == e_root))) */
      if((a->v[i]) && ((a->v[i] == d)))
        {
          b->l_r  = i; /* we consider here that 'a' is on the left handside of 'b'*/
          a->b[i] = b;
        }
      /* if((d->v[i]) && ((d->v[i] == a) || (e_root && d->b[i] == e_root))) */
      if((d->v[i]) && ((d->v[i] == a)))
        {
          b->r_l = i; /* we consider here that 'd' is on the right handside of 'b'*/
          d->b[i] = b;
        }
    }

  if(a->tax) {b->r_l = 0; for(i=0;i<3;i++) if(d->v[i]==a) {b->l_r = i; break;}}

  b->l_v1 = b->l_v2 = b->r_v1 = b->r_v2 = -1;
  for(i=0;i<3;i++)
    {
      if(b->left->v[i] != b->rght)
        {
          if(b->l_v1 < 0) b->l_v1 = i;
          else            b->l_v2 = i;
        }

      if(b->rght->v[i] != b->left)
        {
          if(b->r_v1 < 0) b->r_v1 = i;
          else            b->r_v2 = i;
        }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Restrict_To_Coding_Position(align **data, option *io)
{
  int i,j,curr_pos;

  if(io->codpos != -1)
    {
      for(i=0;i<io->n_otu;i++)
    {
      curr_pos = 0;
      for(j=io->codpos-1;j<data[i]->len;j+=3)
        {
          data[i]->state[curr_pos] = data[i]->state[j];
          curr_pos++;
        }
      data[i]->len /= 3;
    }
    }
}



void Uppercase(char *ch)
{
  /* convert ch to upper case -- either ASCII or EBCDIC */
   *ch = isupper((int)*ch) ? *ch : toupper((int)*ch);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Lowercase(char *ch)
{
  /* convert ch to upper case -- either ASCII or EBCDIC */
  *ch = isupper((int)*ch) ? tolower((int)*ch) : *ch;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

calign *Compact_Data(align **data, option *io)
{
  calign *cdata_tmp,*cdata;
  int i,j,k,site;
  int n_patt,which_patt;
  char **sp_names;
  int n_otu;
  pnode *proot;
  int compress;
  int n_ambigu,is_ambigu;
  scalar_dbl *io_wght;
  phydbl len,inc,n_invar;


  n_otu      = io->n_otu;
  n_patt     = 0;
  which_patt = 0;

  sp_names = (char **)mCalloc(n_otu,sizeof(char *));
  for(i=0;i<n_otu;i++)
    {
      sp_names[i] = (char *)mCalloc(T_MAX_NAME,sizeof(char));
      strcpy(sp_names[i],data[i]->name);
    }

  cdata_tmp = Make_Calign(n_otu,data[0]->len,io->state_len,data[0]->len,sp_names);
  Init_Calign(n_otu,data[0]->len,data[0]->len,cdata_tmp);

  proot = (pnode *)Create_Pnode(T_MAX_ALPHABET);

  for(i=0;i<n_otu;i++) Free(sp_names[i]);
  Free(sp_names);

  if(data[0]->len%io->state_len)
    {
      PhyML_Fprintf(stderr,"\n. Sequence length is not a multiple of %d\n",io->state_len);
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
    }

  // Read in weights given in input file
  io_wght = NULL;
  if(io->has_io_weights == YES) 
    {
      io_wght = Read_Io_Weights(io);
      if(Scalar_Len(io_wght) - data[0]->len > 0)
        {
          PhyML_Fprintf(stderr,"\n. Sequence length (%d) differs from number of weights (%d).\n",
                        data[0]->len,
                        Scalar_Len(io_wght));
          Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
        }
    }

  compress = io->colalias;
  n_ambigu = 0;
  is_ambigu = NO;

  if(!io->quiet && !compress) PhyML_Printf("\n. WARNING: sequences are not compressed !\n");

  inc = -1.0;
  len = 0.0;
  Fors(site,data[0]->len,io->state_len)
    {      
      if(io->has_io_weights == YES)
        inc = Scalar_Elem(site,io_wght);
      else                
        inc = 1.;
      
      // Sequence length taking into account input weights, if any
      len += inc;

      if(io->rm_ambigu == YES)
        {
          is_ambigu = NO;
          for(j=0;j<n_otu;j++) if(Is_Ambigu(data[j]->state+site,io->datatype,io->state_len)) break;
          if(j != n_otu)
            {
              is_ambigu = YES;
              n_ambigu++;
            }
        }
      
      if(!is_ambigu)
        {
          if(compress)
            {
              which_patt = -1;
              
              Traverse_Prefix_Tree(site,-1,&which_patt,&n_patt,data,io,proot);
              if(which_patt == n_patt-1) /* New pattern found */
                {
                  n_patt--;
                  k=n_patt;
                }
              else
                {
                  k = n_patt-10;
                }
            }
          else
            {
              k = n_patt;
            }
          
          if(k == n_patt) /* add a new site pattern */
            {
              for(j=0;j<n_otu;j++)
                Copy_One_State(data[j]->state+site,
                               cdata_tmp->c_seq[j]->state+n_patt*io->state_len,
                               io->state_len);
              
              for(i=0;i<n_otu;i++)
                {
                  for(j=0;j<n_otu;j++)
                    {
                      if(!(Are_Compatible(cdata_tmp->c_seq[i]->state+n_patt*io->state_len,
                                          cdata_tmp->c_seq[j]->state+n_patt*io->state_len,
                                          io->state_len,
                                          io->datatype))) break;
                    }
                  if(j != n_otu) break;
                }
              
              if((j == n_otu) && (i == n_otu)) /* all characters at that site are compatible with one another:
                                                  the site may be invariant */
                {
                  for(j=0;j<n_otu;j++)
                    {
                      cdata_tmp->invar[n_patt] = Assign_State(cdata_tmp->c_seq[j]->state+n_patt*io->state_len,
                                                              io->datatype,
                                                              io->state_len);
                      
                      if(cdata_tmp->invar[n_patt] > -1.) break;
                    }
                }
              else cdata_tmp->invar[n_patt] = -1;
              
              cdata_tmp->sitepatt[site] = n_patt;
              cdata_tmp->wght[n_patt] += inc;
              n_patt += 1;
            }
          else
            {
              cdata_tmp->sitepatt[site] = which_patt;
              cdata_tmp->wght[which_patt] += inc;
            }
        }
    }
  
  data[0]->len -= n_ambigu;
  
  cdata_tmp->init_len                   = data[0]->len;
  cdata_tmp->crunch_len                 = n_patt;
  for(i=0;i<n_otu;i++) cdata_tmp->c_seq[i]->len = n_patt;
  for(i=0;i<n_otu;i++) cdata_tmp->c_seq[i]->num = i;
  
  if(!io->quiet) PhyML_Printf("\n. %d patterns found (out of a total of %d sites). \n",n_patt,data[0]->len);

  if((io->rm_ambigu == YES) && (n_ambigu > 0)) PhyML_Printf("\n. Removed %d columns of the alignment as they contain ambiguous characters (e.g., gaps) \n",n_ambigu);

  n_invar=0.0;
  for(i=0;i<cdata_tmp->crunch_len;i++) if(cdata_tmp->invar[i] > -1.) n_invar+=cdata_tmp->wght[i];

  if(io->quiet == NO) 
    {
      if((n_invar - ceil(n_invar)) < 1.E-10)     
        PhyML_Printf("\n. %d sites without polymorphism (%.2f%c).\n",(int)n_invar,100.*(phydbl)n_invar/len,'%');
      else
        PhyML_Printf("\n. %f sites without polymorphism (%.2f%c).\n",n_invar,100.*(phydbl)n_invar/len,'%');
    }

  cdata_tmp->obs_pinvar = (phydbl)n_invar/len;

  cdata_tmp->io = io;

  
  if(io->datatype == NT)      Get_Base_Freqs(cdata_tmp);
  else if(io->datatype == AA) Get_AA_Freqs(cdata_tmp);
  else {/* Uniform state frequency distribution.*/}

  cdata = Copy_Cseq(cdata_tmp,io);

  Free_Calign(cdata_tmp);
  Free_Prefix_Tree(proot,T_MAX_ALPHABET);

  Check_Ambiguities(cdata,io->datatype,io->state_len);
  Set_D_States(cdata,io->datatype,io->state_len);

  if(io_wght != NULL) Free_Scalar_Dbl(io_wght);

  return cdata;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


calign *Compact_Cdata(calign *data, option *io)
{
  calign *cdata;
  int i,j,k,site;
  int n_patt,which_patt;
  int n_otu;

  n_otu = data->n_otu;

  cdata         = (calign *)mCalloc(1,sizeof(calign));
  cdata->n_otu  = n_otu;
  cdata->c_seq  = (align **)mCalloc(n_otu,sizeof(align *));
  cdata->wght   = (phydbl *)mCalloc(data->crunch_len,sizeof(phydbl));
  cdata->b_frq  = (phydbl *)mCalloc(io->mod->ns,sizeof(phydbl));
  cdata->ambigu = (short int *)mCalloc(data->crunch_len,sizeof(short int));
  cdata->invar  = (short int *)mCalloc(data->crunch_len,sizeof(short int));

  cdata->crunch_len = cdata->init_len = -1;
  for(j=0;j<n_otu;j++)
    {
      cdata->c_seq[j]            = (align *)mCalloc(1,sizeof(align));
      cdata->c_seq[j]->name      = (char *)mCalloc(T_MAX_NAME,sizeof(char));
      strcpy(cdata->c_seq[j]->name,data->c_seq[j]->name);
      cdata->c_seq[j]->state     = (char *)mCalloc(data->crunch_len,sizeof(char));
      cdata->c_seq[j]->is_ambigu = (short int *)mCalloc(data->crunch_len,sizeof(short int));
      cdata->c_seq[j]->state[0]  = data->c_seq[j]->state[0];
    }


  n_patt = which_patt =  0;

  for(site=0;site<data->crunch_len;site++)
    {
      if(data->wght[site] > 0.0)
        {
          for(k=0;k<n_patt;k++)
            {
              for(j=0;j<n_otu;j++)
                {
                  if(strncmp(cdata->c_seq[j]->state+k*io->state_len,
                             data->c_seq[j]->state+site*io->state_len,
                             io->state_len))
                    break;
                }
              
              if(j == n_otu)
                {
                  which_patt = k;
                  break;
                }
            }
                    
          if(k == n_patt)
            {
              for(j=0;j<n_otu;j++) Copy_One_State(data->c_seq[j]->state+site*io->state_len,
                                          cdata->c_seq[j]->state+n_patt*io->state_len,
                                          io->state_len);
              
              for(i=0;i<n_otu;i++)
                {
                  for(j=0;j<n_otu;j++)
                    {
                      if(!(Are_Compatible(cdata->c_seq[i]->state+n_patt*io->state_len,
                                          cdata->c_seq[j]->state+n_patt*io->state_len,
                                          io->state_len,
                                          io->datatype))) break;
                    }
                  if(j != n_otu) break;
                }
              
              if((j == n_otu) && (i == n_otu))
                {
                  for(j=0;j<n_otu;j++)
                    {
                      cdata->invar[n_patt] = Assign_State(cdata->c_seq[j]->state+n_patt*io->state_len,
                                                          io->datatype,
                                                          io->state_len);
                      
                      if(cdata->invar[n_patt] > -1.) break;
                    }
                }
              else cdata->invar[n_patt] = -1;
              
              cdata->wght[n_patt] += data->wght[site];
              n_patt+=1;
            }
          else cdata->wght[which_patt] += data->wght[site];
          
          /*       Print_Site(cdata,k,n_otu,"\n",io->stepsize); */
        }
    }
  
  cdata->init_len   = data->crunch_len;
  cdata->crunch_len = n_patt;
  for(i=0;i<n_otu;i++) cdata->c_seq[i]->len = n_patt;

  if(io->datatype == NT)      Get_Base_Freqs(cdata);
  else if(io->datatype == AA) Get_AA_Freqs(cdata);
  else {/* Not implemented yet */}

  return cdata;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Traverse_Prefix_Tree(int site, int seqnum, int *patt_num, int *n_patt, align **data, option *io, pnode *n)
{
  if(seqnum == io->n_otu-1)
    {
      n->weight++;
      if(n->weight == 1)
        {
          n->num = *n_patt;
          (*n_patt) += 1;
        }
      (*patt_num) = n->num;
      return;
    }
  else
    {
      int next_state;
      
      next_state = -1;
      next_state = Assign_State_With_Ambiguity(data[seqnum+1]->state+site,
                                               io->datatype,
                                               io->state_len);
      
      if(!n->next[next_state]) n->next[next_state] = Create_Pnode(T_MAX_ALPHABET);
      Traverse_Prefix_Tree(site,seqnum+1,patt_num,n_patt,data,io,n->next[next_state]);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


pnode *Create_Pnode(int size)
{
  pnode *n;
  int i;

  n = (pnode *)mCalloc(1,sizeof(pnode ));
  n->next = (pnode **)mCalloc(size,sizeof(pnode *));
  for(i=0;i<size;i++) n->next[i] = NULL;
  n->weight = 0;
  n->num = -1;
  return n;
}
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Get_Base_Freqs(calign *data)
{
  int i,j,k;
  phydbl A,C,G,T;
  phydbl fA,fC,fG,fT;
  phydbl w;

  fA = fC = fG = fT = .25;

  for(k=0;k<8;k++)
    {
      A = C = G = T = .0;
      for(i=0;i<data->n_otu;i++)
        {
          for(j=0;j<data->crunch_len;j++)
            {
              w = data->wght[j];
              if(w)
                {
                  switch(data->c_seq[i]->state[j])
                    {
                    case 'A' : A+=w;
                      break;
                    case 'C' : C+=w;
                      break;
                    case 'G' : G+=w;
                      break;
                    case 'T' : T+=w;
                      break;
                    case 'U' : T+=w;
                      break;
                    case 'M' : C+=w*fC/(fC+fA); A+=w*fA/(fA+fC);
                      break;
                    case 'R' : G+=w*fG/(fA+fG); A+=w*fA/(fA+fG);
                      break;
                    case 'W' : T+=w*fT/(fA+fT); A+=w*fA/(fA+fT);
                      break;
                    case 'S' : C+=w*fC/(fC+fG); G+=w*fG/(fC+fG);
                      break;
                    case 'Y' : C+=w*fC/(fC+fT); T+=w*fT/(fT+fC);
                      break;
                    case 'K' : G+=w*fG/(fG+fT); T+=w*fT/(fT+fG);
                      break;
                    case 'B' : C+=w*fC/(fC+fG+fT); G+=w*fG/(fC+fG+fT); T+=w*fT/(fC+fG+fT);
                      break;
                    case 'D' : A+=w*fA/(fA+fG+fT); G+=w*fG/(fA+fG+fT); T+=w*fT/(fA+fG+fT);
                      break;
                    case 'H' : A+=w*fA/(fA+fC+fT); C+=w*fC/(fA+fC+fT); T+=w*fT/(fA+fC+fT);
                      break;
                    case 'V' : A+=w*fA/(fA+fC+fG); C+=w*fC/(fA+fC+fG); G+=w*fG/(fA+fC+fG);
                      break;
                    case 'N' : case 'X' : case '?' : case 'O' : case '-' :
                      A+=w*fA; C+=w*fC; G+=w*fG; T+=w*fT; break;
                    default : break;
                    }
                }
            }
        }
      fA = A/(A+C+G+T);
      fC = C/(A+C+G+T);
      fG = G/(A+C+G+T);
      fT = T/(A+C+G+T);
    }
  
  data->b_frq[0] = fA;
  data->b_frq[1] = fC;
  data->b_frq[2] = fG;
  data->b_frq[3] = fT;  
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Get_AA_Freqs(calign *data)
{
  int i,j,k;
  phydbl A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y;
  phydbl fA,fC,fD,fE,fF,fG,fH,fI,fK,fL,fM,fN,fP,fQ,fR,fS,fT,fV,fW,fY;
  int w;
  phydbl sum;

  fA = fC = fD = fE = fF = fG = fH = fI = fK = fL =
  fM = fN = fP = fQ = fR = fS = fT = fV = fW = fY = 1./20.;

  for(k=0;k<8;k++)
    {
      A = C = D = E = F = G = H = I = K = L =
        M = N = P = Q = R = S = T = V = W = Y = .0;
      
      for(i=0;i<data->n_otu;i++)
        {
          for(j=0;j<data->crunch_len;j++)
            {
              w = data->wght[j];
              if(w)
                {
                  switch(data->c_seq[i]->state[j])
                    {
                    case 'A' : A+=w;		break;
                    case 'C' : C+=w;		break;
                    case 'D' : D+=w;		break;
                    case 'E' : E+=w;		break;
                    case 'F' : F+=w;		break;
                    case 'G' : G+=w;		break;
                    case 'H' : H+=w;		break;
                    case 'I' : I+=w;		break;
                    case 'K' : K+=w;		break;
                    case 'L' : L+=w;		break;
                    case 'M' : M+=w;		break;
                    case 'N' : N+=w;		break;
                    case 'P' : P+=w;		break;
                    case 'Q' : Q+=w;		break;
                    case 'R' : R+=w;		break;
                    case 'S' : S+=w;		break;
                    case 'T' : T+=w;		break;
                    case 'V' : V+=w;		break;
                    case 'W' : W+=w;		break;
                    case 'Y' : Y+=w;		break;
                    case 'Z' : Q+=w;		break;
                    case 'X' : case '?' : case 'O' : case '-' :
                      A+=w*fA;
                      C+=w*fC;
                      D+=w*fD;
                      E+=w*fE;
                      F+=w*fF;
                      G+=w*fG;
                      H+=w*fH;
                      I+=w*fI;
                      K+=w*fK;
                      L+=w*fL;
                      M+=w*fM;
                      N+=w*fN;
                      P+=w*fP;
                      Q+=w*fQ;
                      R+=w*fR;
                      S+=w*fS;
                      T+=w*fT;
                      V+=w*fV;
                      W+=w*fW;
                      Y+=w*fY;
                      break;
                    default : break;
                    }
                }
            }
        }
      sum = (A+C+D+E+F+G+H+I+K+L+M+N+P+Q+R+S+T+V+W+Y);
      fA = A/sum;      fC = C/sum;      fD = D/sum;      fE = E/sum;
      fF = F/sum;      fG = G/sum;      fH = H/sum;      fI = I/sum;
      fK = K/sum;      fL = L/sum;      fM = M/sum;      fN = N/sum;
      fP = P/sum;      fQ = Q/sum;      fR = R/sum;      fS = S/sum;
      fT = T/sum;      fV = V/sum;      fW = W/sum;      fY = Y/sum;
    }
  
  data->b_frq[0]  = fA;  data->b_frq[1]  = fR;  data->b_frq[2]  = fN;  data->b_frq[3]  = fD;
  data->b_frq[4]  = fC;  data->b_frq[5]  = fQ;  data->b_frq[6]  = fE;  data->b_frq[7]  = fG;
  data->b_frq[8]  = fH;  data->b_frq[9]  = fI;  data->b_frq[10] = fL;  data->b_frq[11] = fK;
  data->b_frq[12] = fM;  data->b_frq[13] = fF;  data->b_frq[14] = fP;  data->b_frq[15] = fS;
  data->b_frq[16] = fT;  data->b_frq[17] = fW;  data->b_frq[18] = fY;  data->b_frq[19] = fV; 
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

// Swap the nodes on the left and right of e1 with the nodes
// on the left and right of e2 respectively, or on the
// right and left of e2 if swap == YES

void Swap_Nodes_On_Edges(t_edge *e1, t_edge *e2, int swap, t_tree *tree)
{
  t_node *buff;

  e1->left->l[e1->l_r] = e1->l->v;
  e1->rght->l[e1->r_l] = e1->l->v;

  e2->left->l[e2->l_r] = e2->l->v;
  e2->rght->l[e2->r_l] = e2->l->v;

  printf("\n. Swap edge %d (%d %d) with %d (%d %d)",e1->num,e1->left->num,e1->rght->num,e2->num,e2->left->num,e2->rght->num);

  if(swap == NO)
    {
      buff = e1->left;
      e1->left = e2->left;
      e2->left = buff;

      buff = e1->rght;
      e1->rght = e2->rght;
      e2->rght = buff;

    }
  else
    {
      buff = e1->left;
      e1->left = e2->rght;
      e2->rght = buff;

      buff = e1->rght;
      e1->rght = e2->left;
      e2->left = buff;
    }

  Connect_One_Edge_To_Two_Nodes(e1->left,e1->rght,e1,tree);
  Connect_One_Edge_To_Two_Nodes(e2->left,e2->rght,e2,tree);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
/* As opposed to Connect_Edges_To_Nodes_Recur, the ordering of 
   edges connected to tips does not depend on the topology
*/
void Connect_Edges_To_Nodes_Serial(t_tree *tree)
{
  int i,j;

  /* Reset */
  For(i,2*tree->n_otu-1) for(j=0;j<3;j++) tree->a_nodes[i]->b[j] = NULL;

  tree->num_curr_branch_available = 0;
  
  for(i=0;i<tree->n_otu;i++)
    {
      if(!tree->a_nodes[i]->tax) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);

      Connect_One_Edge_To_Two_Nodes(tree->a_nodes[i],
                                    tree->a_nodes[i]->v[0],
                                    tree->a_edges[tree->num_curr_branch_available],
                                    tree);
      if(tree->n_otu == 2) break;
    }

  for(i=tree->n_otu;i<2*tree->n_otu-2;i++)
    {
      if(tree->a_nodes[i] == tree->n_root) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);

      for(j=0;j<3;j++) 
        if(!tree->a_nodes[i]->b[j])
          Connect_One_Edge_To_Two_Nodes(tree->a_nodes[i],
                                        tree->a_nodes[i]->v[j],
                                        tree->a_edges[tree->num_curr_branch_available],
                                        tree);
    }
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void Connect_Edges_To_Nodes_Recur(t_node *a, t_node *d, t_tree *tree)
{
  int i;

  Connect_One_Edge_To_Two_Nodes(a,d,tree->a_edges[tree->num_curr_branch_available],tree);

  if(d->tax) return;
  else 
    for(i=0;i<3;i++) 
      if(d->v[i] != a) /* Don't add d->b[i] != tree->e_root condition here since tree is not wired yet... */ 
        Connect_Edges_To_Nodes_Recur(d,d->v[i],tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Connect_One_Edge_To_Two_Nodes(t_node *a, t_node *d, t_edge *b, t_tree *tree)
{
  int i,dir_a_d,dir_d_a;

  if(a == NULL || d == NULL) 
    {
      PhyML_Fprintf(stderr,"\n. a: %d d: %d",a?a->num:-1,d?d->num:-1);
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
    }

  dir_a_d = -1;
  for(i=0;i<3;i++) if(a->v[i] == d) { dir_a_d = i; break; }

  dir_d_a = -1;
  for(i=0;i<3;i++) if(d->v[i] == a) { dir_d_a = i; break; }

  if(dir_a_d == -1) Generic_Exit(__FILE__,__LINE__,__FUNCTION__); 
  if(dir_d_a == -1) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);    

  a->b[dir_a_d] = b;
  d->b[dir_d_a] = b;
  b->num        = tree->num_curr_branch_available;
  b->left       = a;
  b->rght       = d;
  if(a->tax) {b->rght = a; b->left = d;} /* root */
  /* a tip is necessarily on the righthand side of the t_edge */

  tree->num_curr_branch_available += 1;

  (b->left == a)?
    (Set_Edge_Dirs(b,a,d,tree)):
    (Set_Edge_Dirs(b,d,a,tree));

  b->l->v            = a->l[b->l_r];
  if(a->tax) b->l->v = a->l[b->r_l];
  b->l_old->v        = b->l->v;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Update_Dirs(t_tree *tree)
{
  int i;
  int buff;
  t_edge *b;

  b = NULL;
  buff = -1;
  For(i,2*tree->n_otu-3)
    {
      b = tree->a_edges[i];

      if((!b->left->tax) && (b->left->v[b->l_v1]->num < b->left->v[b->l_v2]->num))
        {
          buff    = b->l_v1;
          b->l_v1 = b->l_v2;
          b->l_v2 = buff;
        }
      if((!b->rght->tax) && (b->rght->v[b->r_v1]->num < b->rght->v[b->r_v2]->num))
        {
          buff    = b->r_v1;
          b->r_v1 = b->r_v2;
          b->r_v2 = buff;
        }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Exit(char *message)
{
  fflush(NULL);
  PhyML_Fprintf(stderr,"%s",message);
  exit(1);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void *mCalloc(int nb, size_t size)
{
  void *allocated;

  if((allocated = calloc((size_t)nb,size)) != NULL)
/*   if((allocated = malloc((size_t)nb*(size_t)size)) != NULL) */
    return allocated;
  else
    Generic_Exit(__FILE__,__LINE__,__FUNCTION__);    

  return NULL;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void *mRealloc(void *p,int nb, size_t size)
{
  if((p = realloc(p,(size_t)nb*size)) != NULL)
    return p;
  else
    Exit("\n. Err.: low memory\n");

  return NULL;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


/* t_tree *Make_Light_Tree_Struct(int n_otu) */
/* { */
/*   t_tree *tree; */
/*   int i; */

/*   tree          = (t_tree *)mCalloc(1,sizeof(t_tree )); */
/*   tree->a_edges = (t_edge **)mCalloc(2*n_otu-3,sizeof(t_edge *)); */
/*   tree->a_nodes   = (t_node **)mCalloc(2*n_otu-2,sizeof(t_node *)); */
/*   tree->n_otu   = n_otu; */

/*   For(i,2*n_otu-3) */
/*     tree->a_edges[i] = Make_Edge_Light(NULL,NULL,i); */

/*   For(i,2*n_otu-2) */
/*     tree->a_nodes[i] = Make_Node_Light(i); */

/*   return tree; */
/* } */

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


int Sort_Phydbl_Decrease(const void *a, const void *b)
{
    if((*(phydbl *)(a)) >= (*(phydbl *)(b))) return -1;
    else return 1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Qksort_Int(int *A, int *B, int ilo, int ihi)
{
    phydbl pivot;	// pivot value for partitioning array
    int ulo, uhi;	// indices at ends of unpartitioned region
    int ieq;		// least index of array entry with value equal to pivot
    int tempEntry;	// temporary entry used for swapping

    if (ilo >= ihi) {
    return;
    }
    // Select a pivot value.
    pivot = A[(ilo + ihi)/2];
    // Initialize ends of unpartitioned region and least index of entry
    // with value equal to pivot.
    ieq = ulo = ilo;
    uhi = ihi;
    // While the unpartitioned region is not empty, try to reduce its size.
    while (ulo <= uhi) {
      if (A[uhi] > pivot) {
        // Here, we can reduce the size of the unpartitioned region and
        // try again.
        uhi--;
    } else {
        // Here, A[uhi] <= pivot, so swap entries at indices ulo and
        // uhi.
        tempEntry = A[ulo];
        A[ulo]    = A[uhi];
        A[uhi]    = tempEntry;

        if(B)
          {
            tempEntry = B[ulo];
            B[ulo]    = B[uhi];
            B[uhi]    = tempEntry;
          }
        
        // After the swap, A[ulo] <= pivot.
        if (A[ulo] < pivot) 
          {
            // Swap entries at indices ieq and ulo.
            tempEntry = A[ieq];
            A[ieq] = A[ulo];
            A[ulo] = tempEntry;
            

            if(B)
              {
                tempEntry = B[ieq];
                B[ieq] = B[ulo];
                B[ulo] = tempEntry;
              }
            
            
            // After the swap, A[ieq] < pivot, so we need to change
            // ieq.
            ieq++;
            // We also need to change ulo, but we also need to do
            // that when A[ulo] = pivot, so we do it after this if
            // statement.
          }
        // Once again, we can reduce the size of the unpartitioned
        // region and try again.
        ulo++;
      }
    }
    // Now, all entries from index ilo to ieq - 1 are less than the pivot
    // and all entries from index uhi to ihi + 1 are greater than the
    // pivot.  So we have two regions of the array that can be sorted
    // recursively to put all of the entries in order.
    Qksort_Int(A, B, ilo, ieq - 1);
    Qksort_Int(A, B, uhi + 1, ihi);    
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


/* Sort in ascending order. Elements in B (if provided) are also re-ordered according to the ordering of A  */
void Qksort(phydbl *A, phydbl *B, int ilo, int ihi)
{
    phydbl pivot;	// pivot value for partitioning array
    int ulo, uhi;	// indices at ends of unpartitioned region
    int ieq;		// least index of array entry with value equal to pivot
    phydbl tempEntry;	// temporary entry used for swapping

    if (ilo >= ihi) {
    return;
    }
    // Select a pivot value.
    pivot = A[(ilo + ihi)/2];
    // Initialize ends of unpartitioned region and least index of entry
    // with value equal to pivot.
    ieq = ulo = ilo;
    uhi = ihi;
    // While the unpartitioned region is not empty, try to reduce its size.
    while (ulo <= uhi) {
      if (A[uhi] > pivot) {
        // Here, we can reduce the size of the unpartitioned region and
        // try again.
        uhi--;
    } else {
        // Here, A[uhi] <= pivot, so swap entries at indices ulo and
        // uhi.
        tempEntry = A[ulo];
        A[ulo]    = A[uhi];
        A[uhi]    = tempEntry;

        if(B)
          {
        tempEntry = B[ulo];
        B[ulo]    = B[uhi];
        B[uhi]    = tempEntry;
          }



        // After the swap, A[ulo] <= pivot.
        if (A[ulo] < pivot) {
        // Swap entries at indices ieq and ulo.
        tempEntry = A[ieq];
        A[ieq] = A[ulo];
        A[ulo] = tempEntry;


        if(B)
          {
            tempEntry = B[ieq];
            B[ieq] = B[ulo];
            B[ulo] = tempEntry;
          }


        // After the swap, A[ieq] < pivot, so we need to change
        // ieq.
        ieq++;
        // We also need to change ulo, but we also need to do
        // that when A[ulo] = pivot, so we do it after this if
        // statement.
        }
        // Once again, we can reduce the size of the unpartitioned
        // region and try again.
        ulo++;
    }
    }
    // Now, all entries from index ilo to ieq - 1 are less than the pivot
    // and all entries from index uhi to ihi + 1 are greater than the
    // pivot.  So we have two regions of the array that can be sorted
    // recursively to put all of the entries in order.
    Qksort(A, B, ilo, ieq - 1);
    Qksort(A, B, uhi + 1, ihi);
}

/********************************************************/

void Qksort_Matrix(phydbl **A, int col, int ilo, int ihi)
{
    phydbl pivot;	// pivot value for partitioning array
    int ulo, uhi;	// indices at ends of unpartitioned region
    int ieq;		// least index of array entry with value equal to pivot
    phydbl *tempEntry;	// temporary entry used for swapping

    tempEntry = NULL;

    if (ilo >= ihi) {
    return;
    }
    // Select a pivot value.
    pivot = A[(ilo + ihi)/2][col];
    // Initialize ends of unpartitioned region and least index of entry
    // with value equal to pivot.
    ieq = ulo = ilo;
    uhi = ihi;
    // While the unpartitioned region is not empty, try to reduce its size.
    while (ulo <= uhi) {
    if (A[uhi][col] > pivot) {
        // Here, we can reduce the size of the unpartitioned region and
        // try again.
        uhi--;
    } else {
        // Here, A[uhi] <= pivot, so swap entries at indices ulo and
        // uhi.
        tempEntry = A[ulo];
        A[ulo] = A[uhi];
        A[uhi] = tempEntry;
        // After the swap, A[ulo] <= pivot.
        if (A[ulo][col] < pivot) {
        // Swap entries at indices ieq and ulo.
        tempEntry = A[ieq];
        A[ieq] = A[ulo];
        A[ulo] = tempEntry;
        // After the swap, A[ieq] < pivot, so we need to change
        // ieq.
        ieq++;
        // We also need to change ulo, but we also need to do
        // that when A[ulo] = pivot, so we do it after this if
        // statement.
        }
        // Once again, we can reduce the size of the unpartitioned
        // region and try again.
        ulo++;
    }
    }
    // Now, all entries from index ilo to ieq - 1 are less than the pivot
    // and all entries from index uhi to ihi + 1 are greater than the
    // pivot.  So we have two regions of the array that can be sorted
    // recursively to put all of the entries in order.
    Qksort_Matrix(A, col, ilo, ieq - 1);
    Qksort_Matrix(A, col, uhi + 1, ihi);
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

char *Add_Taxa_To_Constraint_Tree(FILE *fp, calign *cdata)
{
  char *line,*long_line;
  t_tree *tree;
  int i,j,open;

  rewind(fp);

  line = Return_Tree_String_Phylip(fp);
  tree = Read_Tree(&line);

  long_line = (char *)mCalloc(T_MAX_LINE,sizeof(char));
  strcpy(long_line,line);
  i = 1;
  open = 1;
  while(open)
    {
      if(line[i]=='(') open++;
      if(line[i]==')') open--;      
      if(i > T_MAX_LINE) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
      i++;
    }
  long_line[i-1] = '\0';
  
  for(i=0;i<cdata->n_otu;i++)
    {
      for(j=0;j<tree->n_otu;j++)
        {
          if(!strcmp(tree->a_nodes[j]->name,cdata->c_seq[i]->name))
            break;
        }
      
      if(j == tree->n_otu)
        {
          strcat(long_line,",");
          strcat(long_line,cdata->c_seq[i]->name);
        }
      
    }
  
  strcat(long_line,");");

  Free_Tree(tree);
  Free(line);


  return long_line;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Check_Constraint_Tree_Taxa_Names(t_tree *tree, calign *cdata)
{
  int i,j,n_otu_tree,n_otu_cdata;

  n_otu_tree  = tree->n_otu;
  n_otu_cdata = cdata->n_otu;

  for(i=0;i<n_otu_tree;i++)
    {
      for(j=0;j<n_otu_cdata;j++)
        {
          if(!strcmp(tree->a_nodes[i]->name,cdata->c_seq[j]->name))
            break;
        }
      
      if(j==n_otu_cdata)
        {
          PhyML_Fprintf(stderr,"\n. '%s' was not found in sequence data set\n",tree->a_nodes[i]->name);
          Warn_And_Exit("\n. PhyML finished prematurely.");
        }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Copy_Tax_Names_To_Tip_Labels(t_tree *tree, calign *data)
{
  int i;

  for(i=0;i<tree->n_otu;i++)
    {
      tree->a_nodes[i]->name = (char *)mCalloc((int)strlen(data->c_seq[i]->name)+1,sizeof(char));
      tree->a_nodes[i]->ori_name = tree->a_nodes[i]->name;
      strcpy(tree->a_nodes[i]->name,data->c_seq[i]->name);
      tree->a_nodes[i]->tax = 1;
      tree->a_nodes[i]->num = i;
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Share_Lk_Struct(t_tree *t_full, t_tree *t_empt)
{
  int i,j,n_otu;
  t_edge *b_e,*b_f;
  t_node *n_e, *n_f;

  n_otu                        = t_full->n_otu;
  t_empt->c_lnL_sorted         = t_full->c_lnL_sorted;
  t_empt->unscaled_site_lk_cat = t_full->unscaled_site_lk_cat;
  t_empt->cur_site_lk          = t_full->cur_site_lk;
  t_empt->old_site_lk          = t_full->old_site_lk;
  t_empt->triplet_struct       = t_full->triplet_struct;
  t_empt->log_lks_aLRT         = t_full->log_lks_aLRT;
  t_empt->site_lk_cat          = t_full->site_lk_cat;
  t_empt->fact_sum_scale       = t_full->fact_sum_scale;
  t_empt->eigen_lr_left        = t_full->eigen_lr_left;
  t_empt->eigen_lr_rght        = t_full->eigen_lr_rght;
  t_empt->dot_prod             = t_full->dot_prod;
  t_empt->expl                 = t_full->expl;

  For(i,2*n_otu-3)
    {
      b_f = t_full->a_edges[i];
      b_e = t_empt->a_edges[i];

      b_e->Pij_rr  = b_f->Pij_rr;
      b_e->tPij_rr = b_f->tPij_rr;

      b_e->nni = b_f->nni;
    }


  for(i=n_otu;i<2*n_otu-2;i++)
    {
      n_f = t_full->a_nodes[i];
      n_e = t_empt->a_nodes[i];

      for(j=0;j<3;j++)
    {
      if(n_f->b[j]->left == n_f)
        {
          if(n_e->b[j]->left == n_e)
        {
          n_e->b[j]->p_lk_left          = n_f->b[j]->p_lk_left;
          n_e->b[j]->p_lk_loc_left      = n_f->b[j]->p_lk_loc_left;
          n_e->b[j]->patt_id_left       = n_f->b[j]->patt_id_left;
          n_e->b[j]->sum_scale_left     = n_f->b[j]->sum_scale_left;
          n_e->b[j]->sum_scale_left_cat = n_f->b[j]->sum_scale_left_cat;
          n_e->b[j]->p_lk_tip_l         = n_f->b[j]->p_lk_tip_l;
        }
          else
        {
          n_e->b[j]->p_lk_rght          = n_f->b[j]->p_lk_left;
          n_e->b[j]->p_lk_loc_rght      = n_f->b[j]->p_lk_loc_left;
          n_e->b[j]->patt_id_rght       = n_f->b[j]->patt_id_left;
          n_e->b[j]->sum_scale_rght     = n_f->b[j]->sum_scale_left;
          n_e->b[j]->sum_scale_rght_cat = n_f->b[j]->sum_scale_left_cat;
          n_e->b[j]->p_lk_tip_r         = n_f->b[j]->p_lk_tip_l;
        }
        }
      else
        {
          if(n_e->b[j]->rght == n_e)
        {
          n_e->b[j]->p_lk_rght          = n_f->b[j]->p_lk_rght;
          n_e->b[j]->p_lk_loc_rght      = n_f->b[j]->p_lk_loc_rght;
          n_e->b[j]->patt_id_rght       = n_f->b[j]->patt_id_rght;
          n_e->b[j]->sum_scale_rght     = n_f->b[j]->sum_scale_rght;
          n_e->b[j]->sum_scale_rght_cat = n_f->b[j]->sum_scale_rght_cat;
          n_e->b[j]->p_lk_tip_r         = n_f->b[j]->p_lk_tip_r;
        }
          else
        {
          n_e->b[j]->p_lk_left          = n_f->b[j]->p_lk_rght;
          n_e->b[j]->p_lk_loc_left      = n_f->b[j]->p_lk_loc_rght;
          n_e->b[j]->patt_id_left       = n_f->b[j]->patt_id_rght;
          n_e->b[j]->sum_scale_left     = n_f->b[j]->sum_scale_rght;
          n_e->b[j]->sum_scale_left_cat = n_f->b[j]->sum_scale_rght_cat;
          n_e->b[j]->p_lk_tip_l         = n_f->b[j]->p_lk_tip_r;
        }
        }
    }
    }

  for(i=0;i<n_otu;i++)
    {
      n_f = t_full->a_nodes[i];
      n_e = t_empt->a_nodes[i];

      if(n_f->b[0]->rght == n_f)
        {
          n_e->b[0]->p_lk_rght          = n_f->b[0]->p_lk_rght;
          n_e->b[0]->p_lk_loc_rght      = n_f->b[0]->p_lk_loc_rght;
          n_e->b[0]->patt_id_rght       = n_f->b[0]->patt_id_rght;
          n_e->b[0]->sum_scale_rght     = n_f->b[0]->sum_scale_rght;
          n_e->b[0]->sum_scale_rght_cat = n_f->b[0]->sum_scale_rght_cat;
          n_e->b[0]->p_lk_tip_r         = n_f->b[0]->p_lk_tip_r;
        }
      else
        {
          Generic_Exit(__FILE__,__LINE__,__FUNCTION__);    
        }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Share_Spr_Struct(t_tree *t_full, t_tree *t_empt)
{
  t_empt->size_spr_list = t_full->size_spr_list;
  t_empt->spr_list      = t_full->spr_list;
  t_empt->best_spr      = t_full->best_spr;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Share_Pars_Struct(t_tree *t_full, t_tree *t_empt)
{
  int i;

  t_empt->site_pars = t_full->site_pars;
  t_empt->step_mat  = t_full->step_mat;

  For(i,2*t_full->n_otu-3)
    {
      t_empt->a_edges[i]->ui_l     = t_full->a_edges[i]->ui_l;
      t_empt->a_edges[i]->ui_r     = t_full->a_edges[i]->ui_r;

      t_empt->a_edges[i]->pars_l   = t_full->a_edges[i]->pars_l;
      t_empt->a_edges[i]->pars_r   = t_full->a_edges[i]->pars_r;

      t_empt->a_edges[i]->p_pars_l = t_full->a_edges[i]->p_pars_l;
      t_empt->a_edges[i]->p_pars_r = t_full->a_edges[i]->p_pars_r;
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int Sort_Edges_NNI_Score(t_tree *tree, t_edge **sorted_edges, int n_elem)
{
  int i,j,done;
  t_edge *buff;

  do
    {
      done = YES;
      for(i=0;i<n_elem-1;i++)
        {
          for(j=i+1;j<n_elem;j++)
            {
              if(sorted_edges[j]->nni->score  < sorted_edges[i]->nni->score)
                {
                  done = NO;
                  buff = sorted_edges[j];
                  sorted_edges[j] = sorted_edges[i];
                  sorted_edges[i] = buff;
                }
            }
        }
    }
  while(done == NO);

  return 1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int Sort_Edges_Depth(t_tree *tree, t_edge **sorted_edges, int n_elem)
{
  int i,j;
  t_edge *buff;
  phydbl *depth,buff_depth;

  depth = (phydbl *)mCalloc(n_elem,sizeof(phydbl));

  for(i=0;i<n_elem;i++)
    depth[i] =
    sorted_edges[i]->left->bip_size[sorted_edges[i]->l_r] *
    sorted_edges[i]->rght->bip_size[sorted_edges[i]->r_l] ;


  for(i=0;i<n_elem-1;i++)
    {
      for(j=i+1;j<n_elem;j++)
    {
      if(depth[i] > depth[j])
        {
          buff = sorted_edges[i];
          sorted_edges[i] = sorted_edges[j];
          sorted_edges[j] = buff;

          buff_depth = depth[i];
          depth[i] = depth[j];
          depth[j] = buff_depth;
        }
    }
    }

  Free(depth);

  return 1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void NNI(t_tree *tree, t_edge *b_fcus, int do_swap)
{
  t_node *v1,*v2,*v3,*v4;
  phydbl lk0, lk1, lk2;
  phydbl lk0_init, lk1_init, lk2_init;
  scalar_dbl *len0,*len1,*len2;
  scalar_dbl *var0,*var1,*var2;
  phydbl l_infa, l_infb;
  phydbl lk_init;
  
  if(tree->prev) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);    

  lk_init                = tree->c_lnL;
  b_fcus->nni->init_lk   = tree->c_lnL;;
  b_fcus->nni->best_conf = 0;
  b_fcus->nni->score     = +1.0;
  lk0 = lk1 = lk2        = UNLIKELY;
  v1 = v2 = v3 = v4      = NULL;

  if(b_fcus->nni->init_l != NULL)  Copy_Scalar_Dbl(b_fcus->l,b_fcus->nni->init_l);
  else                             b_fcus->nni->init_l = Duplicate_Scalar_Dbl(b_fcus->l);

  v1 = b_fcus->left->v[b_fcus->l_v1];
  v2 = b_fcus->left->v[b_fcus->l_v2];
  v3 = b_fcus->rght->v[b_fcus->r_v1];
  v4 = b_fcus->rght->v[b_fcus->r_v2];

  Record_Br_Len(tree);

  if(v1->num < v2->num) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
  if(v3->num < v4->num) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);

  /************************************************************/
  Swap(v2,b_fcus->left,b_fcus->rght,v3,tree);
  Set_Both_Sides(YES,tree);

  MIXT_Set_Alias_Subpatt(YES,tree);
  lk1_init = Update_Lk_At_Given_Edge(b_fcus,tree);
  MIXT_Set_Alias_Subpatt(NO,tree);

  l_infa = 1.;
  l_infb = 1.E-4;
  lk1    = lk1_init;

  if(tree->mod->s_opt->nni_br_len_opt == YES)
    {
      if(tree->mod->s_opt->fast_nni)
        {
          lk1 = Fast_Br_Len(b_fcus,tree,YES);
        }
      else 
        {
          lk1 = Br_Len_Brent(b_fcus,tree);
        }
    }

  if(lk1 < lk1_init - tree->mod->s_opt->min_diff_lk_local)
    {
      PhyML_Printf("\n. %f %f %G",l_infa,l_infb,b_fcus->l->v);
      PhyML_Printf("\n. %f -- %f",lk1_init,lk1);
      PhyML_Printf("\n. Err. in NNI (1)");
    }

  len1 = Duplicate_Scalar_Dbl(b_fcus->l);
  var1 = Duplicate_Scalar_Dbl(b_fcus->l_var);
  Swap(v3,b_fcus->left,b_fcus->rght,v2,tree);
  /************************************************************/


  /************************************************************/
  Swap(v2,b_fcus->left,b_fcus->rght,v4,tree);
  Restore_Br_Len(tree);
  Set_Both_Sides(YES,tree);

  MIXT_Set_Alias_Subpatt(YES,tree);
  lk2_init = Update_Lk_At_Given_Edge(b_fcus,tree);
  MIXT_Set_Alias_Subpatt(NO,tree);

  l_infa = 1.;
  l_infb = 1.E-4;
 
  lk2    = lk2_init;

  if(tree->mod->s_opt->nni_br_len_opt == YES)
    {
      if(tree->mod->s_opt->fast_nni)
        {
          lk2 = Fast_Br_Len(b_fcus,tree,YES);
        }
      else 
        {
          lk2 = Br_Len_Brent(b_fcus,tree);
        }
    }

  if(lk2 < lk2_init - tree->mod->s_opt->min_diff_lk_local)
    {
      PhyML_Printf("\n. %f %f %G",l_infa,l_infb,b_fcus->l->v);
      PhyML_Printf("\n. %f -- %f",lk2_init,lk2);
      PhyML_Printf("\n. Err. in NNI (2)");
   }

  len2 = Duplicate_Scalar_Dbl(b_fcus->l);
  var2 = Duplicate_Scalar_Dbl(b_fcus->l_var);
  Swap(v4,b_fcus->left,b_fcus->rght,v2,tree);
  /************************************************************/



  /************************************************************/
  Restore_Br_Len(tree);
  Set_Both_Sides(YES,tree);

  MIXT_Set_Alias_Subpatt(YES,tree);
  lk0_init = Update_Lk_At_Given_Edge(b_fcus,tree);
  MIXT_Set_Alias_Subpatt(NO,tree);


  if(FABS(lk0_init - lk_init) > tree->mod->s_opt->min_diff_lk_local)
    {
      PhyML_Fprintf(stderr,"\n. lk_init = %f; lk = %f diff = %f l = %G",
                    lk_init,
                    lk0_init,
                    lk_init-lk0_init,
                    b_fcus->l->v);
      PhyML_Fprintf(stderr,"\n. Curr_lnL = %f",Lk(NULL,tree));
      Exit("\n. Err. in NNI (3)");
    }

  l_infa = 1.;
  l_infb = 1.E-4;
  lk0    = lk0_init;

  if(tree->mod->s_opt->nni_br_len_opt == YES)
    {
      if(tree->mod->s_opt->fast_nni)
        {
          lk0 = Fast_Br_Len(b_fcus,tree,YES);
        }
      else 
        {
          lk0 = Br_Len_Brent(b_fcus,tree);
        }
    }

  if(lk0 < lk_init - tree->mod->s_opt->min_diff_lk_local)
    {
      PhyML_Printf("\n. %f %f %f",l_infa,l_infb,b_fcus->l->v);
      PhyML_Printf("\n. %f -- %f",lk0_init,lk0);
      PhyML_Printf("\n. Err. in NNI (3)\n");
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
    }

  len0 = Duplicate_Scalar_Dbl(b_fcus->l);
  var0 = Duplicate_Scalar_Dbl(b_fcus->l_var);
  /************************************************************/

  b_fcus->nni->lk0 = lk0;
  b_fcus->nni->lk1 = lk1;
  b_fcus->nni->lk2 = lk2;
  
  if(b_fcus->nni->l0 == NULL) b_fcus->nni->l0 = Duplicate_Scalar_Dbl(len0);
  else                        Copy_Scalar_Dbl(len0,b_fcus->nni->l0);

  if(b_fcus->nni->l1 == NULL) b_fcus->nni->l1 = Duplicate_Scalar_Dbl(len1);
  else                        Copy_Scalar_Dbl(len1,b_fcus->nni->l1);

  if(b_fcus->nni->l2 == NULL) b_fcus->nni->l2 = Duplicate_Scalar_Dbl(len2);
  else                        Copy_Scalar_Dbl(len2,b_fcus->nni->l2);

  if(b_fcus->nni->v0 == NULL) b_fcus->nni->v0 = Duplicate_Scalar_Dbl(var0);
  else                        Copy_Scalar_Dbl(var0,b_fcus->nni->v0);

  if(b_fcus->nni->v1 == NULL) b_fcus->nni->v1 = Duplicate_Scalar_Dbl(var1);
  else                        Copy_Scalar_Dbl(var1,b_fcus->nni->v1);

  if(b_fcus->nni->v2 == NULL) b_fcus->nni->v2 = Duplicate_Scalar_Dbl(var2);
  else                        Copy_Scalar_Dbl(var2,b_fcus->nni->v2);

  b_fcus->nni->score = lk0 - MAX(lk1,lk2);

  if((b_fcus->nni->score <  tree->mod->s_opt->min_diff_lk_local) &&
     (b_fcus->nni->score > -tree->mod->s_opt->min_diff_lk_local))
    {
      b_fcus->nni->score = .0;
      b_fcus->nni->lk1 = b_fcus->nni->lk0;
      b_fcus->nni->lk2 = b_fcus->nni->lk0;
    }

  if(lk0 > MAX(lk1,lk2))
    {
      b_fcus->nni->best_conf    = 0;
      b_fcus->nni->swap_node_v1 = NULL;
      b_fcus->nni->swap_node_v2 = NULL;
      b_fcus->nni->swap_node_v3 = NULL;
      b_fcus->nni->swap_node_v4 = NULL;
      
      if(b_fcus->nni->best_l == NULL) b_fcus->nni->best_l = Duplicate_Scalar_Dbl(len0);
      else                            Copy_Scalar_Dbl(len0,b_fcus->nni->best_l);

      if(b_fcus->nni->best_v == NULL) b_fcus->nni->best_v = Duplicate_Scalar_Dbl(var0);
      else                            Copy_Scalar_Dbl(var0,b_fcus->nni->best_v);
    }
  else if(lk1 > MAX(lk0,lk2))
    {
      b_fcus->nni->best_conf    = 1;
      b_fcus->nni->swap_node_v1 = v2;
      b_fcus->nni->swap_node_v2 = b_fcus->left;
      b_fcus->nni->swap_node_v3 = b_fcus->rght;
      b_fcus->nni->swap_node_v4 = v3;
      
      if(b_fcus->nni->best_l == NULL) b_fcus->nni->best_l = Duplicate_Scalar_Dbl(len1);
      else                            Copy_Scalar_Dbl(len1,b_fcus->nni->best_l);

      if(b_fcus->nni->best_v == NULL) b_fcus->nni->best_v = Duplicate_Scalar_Dbl(var1);
      else                            Copy_Scalar_Dbl(var1,b_fcus->nni->best_v);
    }
  else if(lk2 > MAX(lk0,lk1))
    {
      b_fcus->nni->best_conf    = 2;
      b_fcus->nni->swap_node_v1 = v2;
      b_fcus->nni->swap_node_v2 = b_fcus->left;
      b_fcus->nni->swap_node_v3 = b_fcus->rght;
      b_fcus->nni->swap_node_v4 = v4;

      if(b_fcus->nni->best_l == NULL) b_fcus->nni->best_l = Duplicate_Scalar_Dbl(len2);
      else                            Copy_Scalar_Dbl(len2,b_fcus->nni->best_l);

      if(b_fcus->nni->best_v == NULL) b_fcus->nni->best_v = Duplicate_Scalar_Dbl(var2);
      else                            Copy_Scalar_Dbl(var2,b_fcus->nni->best_v);
    }
  else
    {
      b_fcus->nni->score        = +1.0;
      b_fcus->nni->best_conf    = 0;
      b_fcus->nni->swap_node_v1 = NULL;
      b_fcus->nni->swap_node_v2 = NULL;
      b_fcus->nni->swap_node_v3 = NULL;
      b_fcus->nni->swap_node_v4 = NULL;

      if(b_fcus->nni->best_l == NULL) b_fcus->nni->best_l = Duplicate_Scalar_Dbl(len0);
      else                            Copy_Scalar_Dbl(len0,b_fcus->nni->best_l);

      if(b_fcus->nni->best_v == NULL) b_fcus->nni->best_v = Duplicate_Scalar_Dbl(var0);
      else                            Copy_Scalar_Dbl(var0,b_fcus->nni->best_v);
    }


  if(do_swap == YES)
    {
      if((lk1 > lk0) || (lk2 > lk0))
        {
          tree->n_swap++;

          if(lk1 > lk2)
            {
              Swap(v2,b_fcus->left,b_fcus->rght,v3,tree);
              
              if(b_fcus->nni->best_l == NULL) b_fcus->nni->best_l = Duplicate_Scalar_Dbl(len1);
              else                            Copy_Scalar_Dbl(len1,b_fcus->nni->best_l);
              
              if(b_fcus->nni->best_v == NULL) b_fcus->nni->best_v = Duplicate_Scalar_Dbl(var1);
              else                            Copy_Scalar_Dbl(var1,b_fcus->nni->best_v);
            }
          else
            {
              Swap(v2,b_fcus->left,b_fcus->rght,v4,tree);
              
              if(b_fcus->nni->best_l == NULL) b_fcus->nni->best_l = Duplicate_Scalar_Dbl(len2);
              else                            Copy_Scalar_Dbl(len2,b_fcus->nni->best_l);
              
              if(b_fcus->nni->best_v == NULL) b_fcus->nni->best_v = Duplicate_Scalar_Dbl(var2);
              else                            Copy_Scalar_Dbl(var2,b_fcus->nni->best_v);
            }
        }
    }
  else
    {
      Restore_Br_Len(tree);
      Update_PMat_At_Given_Edge(b_fcus,tree);
      tree->c_lnL = lk_init;
    }

  Free_Scalar_Dbl(len0);
  Free_Scalar_Dbl(len1);
  Free_Scalar_Dbl(len2);
  Free_Scalar_Dbl(var0);
  Free_Scalar_Dbl(var1);
  Free_Scalar_Dbl(var2);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void NNI_Pars(t_tree *tree, t_edge *b_fcus, int do_swap)
{
  t_node *v1,*v2,*v3,*v4;
  int pars0, pars1, pars2;
  int pars_init;

  pars_init              = tree->c_pars;
  b_fcus->nni->best_conf = 0;
  b_fcus->nni->score     = +1.0;

  pars0 = pars1 = pars2  = 0;
  v1 = v2 = v3 = v4      = NULL;


  v1                     = b_fcus->left->v[b_fcus->l_v1];
  v2                     = b_fcus->left->v[b_fcus->l_v2];
  v3                     = b_fcus->rght->v[b_fcus->r_v1];
  v4                     = b_fcus->rght->v[b_fcus->r_v2];

  if(v1->num < v2->num) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);    
  if(v3->num < v4->num) Generic_Exit(__FILE__,__LINE__,__FUNCTION__); 


  /***********/
  Swap(v2,b_fcus->left,b_fcus->rght,v3,tree);
  Set_Both_Sides(YES,tree);
  pars1 = Update_Pars_At_Given_Edge(b_fcus,tree);
  Swap(v3,b_fcus->left,b_fcus->rght,v2,tree);
  /***********/

  /***********/
  Swap(v2,b_fcus->left,b_fcus->rght,v4,tree);
  Set_Both_Sides(YES,tree);
  pars2 = Update_Pars_At_Given_Edge(b_fcus,tree);
  Swap(v4,b_fcus->left,b_fcus->rght,v2,tree);
  /***********/


  /***********/
  Set_Both_Sides(YES,tree);
  pars0 = Update_Pars_At_Given_Edge(b_fcus,tree);

  if(pars0 != pars_init)
     {
       PhyML_Fprintf(stderr,"\n. pars_init = %d; pars0 = %d\n",
          pars_init,
          pars0);
       Warn_And_Exit("\n. Err. in NNI (3)\n");
     }
   /***********/

   tree->c_pars = pars0;

   b_fcus->nni->score = MIN(pars1,pars2) - pars0;

   if(pars0 < MIN(pars1,pars2))
     {
       b_fcus->nni->best_conf    = 0;
       b_fcus->nni->swap_node_v1 = NULL;
       b_fcus->nni->swap_node_v2 = NULL;
       b_fcus->nni->swap_node_v3 = NULL;
       b_fcus->nni->swap_node_v4 = NULL;
      }
   else if(pars1 < MIN(pars0,pars2))
     {
       b_fcus->nni->best_conf    = 1;
       b_fcus->nni->swap_node_v1 = v2;
       b_fcus->nni->swap_node_v2 = b_fcus->left;
       b_fcus->nni->swap_node_v3 = b_fcus->rght;
       b_fcus->nni->swap_node_v4 = v3;
     }
   else if(pars2 > MIN(pars0,pars1))
     {
       b_fcus->nni->best_conf    = 2;
       b_fcus->nni->swap_node_v1 = v2;
       b_fcus->nni->swap_node_v2 = b_fcus->left;
       b_fcus->nni->swap_node_v3 = b_fcus->rght;
       b_fcus->nni->swap_node_v4 = v4;
     }
   else
     {
       b_fcus->nni->score        = +1.0;
       b_fcus->nni->swap_node_v1 = NULL;
       b_fcus->nni->swap_node_v2 = NULL;
       b_fcus->nni->swap_node_v3 = NULL;
       b_fcus->nni->swap_node_v4 = NULL;
     }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Swap(t_node *a, t_node *b, t_node *c, t_node *d, t_tree *tree)
{
  int ab, ba, cd, dc, bc;
  int i;

  /* \             /d      \             /a
   *  \           /         \           /
   *   \b__...__c/    ->     \b__...__c/
   *   /         \	     /         \
   *  /           \	    /   	\
   * /a            \  	   /d            \
   *
   * nodes b and c are not necessarily on the same branch
   */

  if(!tree) return;

  if(!a || !b || !c || !d) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);    

  ab = ba = cd = dc = bc = -1;
  
  for(i=0;i<3;i++) if(a->v[i] == b) { ab = i; break; }
  for(i=0;i<3;i++) if(b->v[i] == a) { ba = i; break; }
  for(i=0;i<3;i++) if(c->v[i] == d) { cd = i; break; }
  for(i=0;i<3;i++) if(d->v[i] == c) { dc = i; break; }
  for(i=0;i<3;i++) if(b->v[i] == c) { bc = i; break; }
  
  if(ab < 0 || ba < 0 || cd < 0 || dc < 0)
    {
      PhyML_Fprintf(stderr,"\n. ab=%d ba=%d cd=%d dc=%d bc=%d",ab,ba,cd,dc,bc);
      PhyML_Fprintf(stderr,"\n. Nodes %d %d %d %d.",a->num,b->num,c->num,d->num);
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);    
    }
  
  a->v[ab] = c;
  d->v[dc] = b;
  b->v[ba] = d;
  c->v[cd] = a;
  b->b[ba] = d->b[dc];
  c->b[cd] = a->b[ab];
  
  (a->b[ab]->left == b)?
    (a->b[ab]->left = c):
    (a->b[ab]->rght = c);
  
  (d->b[dc]->left == c)?
    (d->b[dc]->left = b):
    (d->b[dc]->rght = b);
  
  for(i=0;i<3;i++)
    {
      if(a->b[ab]->left->v[i] == a->b[ab]->rght) a->b[ab]->l_r = i;
      if(a->b[ab]->rght->v[i] == a->b[ab]->left) a->b[ab]->r_l = i;
      if(d->b[dc]->left->v[i] == d->b[dc]->rght) d->b[dc]->l_r = i;
      if(d->b[dc]->rght->v[i] == d->b[dc]->left) d->b[dc]->r_l = i;
    }
  
  a->b[ab]->l_v1 = a->b[ab]->l_v2 =
  a->b[ab]->r_v1 = a->b[ab]->r_v2 =
  d->b[dc]->l_v1 = d->b[dc]->l_v2 =
  d->b[dc]->r_v1 = d->b[dc]->r_v2 = -1;
  
  for(i=0;i<3;i++)
    {
      if(i != a->b[ab]->l_r)
        {
          if(a->b[ab]->l_v1 < 0) a->b[ab]->l_v1 = i;
          else a->b[ab]->l_v2 = i;
        }
      if(i != a->b[ab]->r_l)
        {
          if(a->b[ab]->r_v1 < 0) a->b[ab]->r_v1 = i;
          else a->b[ab]->r_v2 = i;
        }
      if(i != d->b[dc]->l_r)
        {
          if(d->b[dc]->l_v1 < 0) d->b[dc]->l_v1 = i;
          else d->b[dc]->l_v2 = i;
        }
      if(i != d->b[dc]->r_l)
        {
          if(d->b[dc]->r_v1 < 0) d->b[dc]->r_v1 = i;
          else d->b[dc]->r_v2 = i;
        }
    }

  Update_Dirs(tree);
  
  if(tree->n_root != NULL)
    {
      tree->n_root->v[1] = tree->e_root->left;
      tree->n_root->v[2] = tree->e_root->rght;
      tree->n_root->b[1]->rght = tree->e_root->left;
      tree->n_root->b[2]->rght = tree->e_root->rght;
    }

  if(tree->next)
    Swap(a->next,b->next,c->next,d->next,tree->next);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Update_SubTree_Partial_Lk(t_edge *b_fcus, t_node *a, t_node *d, t_tree *tree)
{
  int i;
  Update_Partial_Lk(tree,b_fcus,a);
  if(d->tax) return;
  else for(i=0;i<3;++i) if(d->v[i] != a) Update_SubTree_Partial_Lk(d->b[i],d,d->v[i],tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Copy_Seq_Names_To_Tip_Labels(t_tree *tree, calign *data)
{
  int i;
  for(i=0;i<tree->n_otu;++i) strcpy(tree->a_nodes[i]->name,data->c_seq[i]->name);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

calign *Copy_Cseq(calign *ori, option *io)
{
  calign *new;
  int i,j,k,n_otu,c_len;
  char **sp_names;

  n_otu = ori->n_otu;
  c_len = ori->crunch_len;

  sp_names = (char **)mCalloc(n_otu,sizeof(char *));
  for(i=0;i<n_otu;i++)
    {
      sp_names[i] = (char *)mCalloc(strlen(ori->c_seq[i]->name)+1,sizeof(char));
      strcpy(sp_names[i],ori->c_seq[i]->name);
    }

  new = Make_Calign(n_otu,c_len+1,io->state_len,ori->init_len,sp_names);
  Init_Calign(n_otu,c_len+1,ori->init_len,new);

  new->obs_pinvar = ori->obs_pinvar;

  for(i=0;i<n_otu;i++) new->c_seq[i]->num = ori->c_seq[i]->num;

  for(i=0;i<ori->init_len;i++) new->sitepatt[i] = ori->sitepatt[i];

  for(j=0;j<ori->crunch_len;j++)
    {
      for(i=0;i<ori->n_otu;i++)
        {
          for(k=0;k<io->state_len;k++)
            {
              new->c_seq[i]->state[j*io->state_len+k] =
                ori->c_seq[i]->state[j*io->state_len+k];
            }
          new->c_seq[i]->is_ambigu[j] = ori->c_seq[i]->is_ambigu[j];
        }
      
      new->wght[j]   = ori->wght[j];
      new->ambigu[j] = ori->ambigu[j];
      new->invar[j]  = ori->invar[j];
    }
  
  for(i=0;i<ori->n_otu;i++)
    {
      new->c_seq[i]->len = ori->c_seq[i]->len;
      strcpy(new->c_seq[i]->name,ori->c_seq[i]->name);
    }
  
  for(i=0;i<ori->n_otu;i++) new->c_seq[i]->state[c_len*io->state_len] = '\0';

  for(i=0;i<T_MAX_ALPHABET;i++) new->b_frq[i] = ori->b_frq[i];

  new->init_len           = ori->init_len;
  new->clean_len          = ori->clean_len;
  new->crunch_len         = ori->crunch_len;
  new->n_otu              = ori->n_otu;
  new->io                 = ori->io;

  for(i=0;i<n_otu;i++) Free(sp_names[i]);
  Free(sp_names);

  return new;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


int Filexists(char *filename)
{
  FILE *fp;
  fp =fopen(filename,"r");
  if (fp) {
    fclose(fp);
    return 1;
  } else
    return 0;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

matrix *K80_dist(calign *data, phydbl g_shape)
{
  int i,j,k;
  int diff;
  matrix *mat;
  phydbl **len;

  len = (phydbl **)mCalloc(data->n_otu,sizeof(phydbl *));
  for(i=0;i<data->n_otu;i++)
    len[i] = (phydbl *)mCalloc(data->n_otu,sizeof(phydbl));

  mat = Make_Mat(data->n_otu);

  Init_Mat(mat,data);

  diff = 0;
  
  For(i,data->c_seq[0]->len)
    {
      for(j=0;j<data->n_otu-1;j++)
        {
          for(k=j+1;k<data->n_otu;k++)
            {
              if(((data->c_seq[j]->state[i] == 'A' || data->c_seq[j]->state[i] == 'G') &&
                  (data->c_seq[k]->state[i] == 'C' || data->c_seq[k]->state[i] == 'T'))||
                 ((data->c_seq[j]->state[i] == 'C' || data->c_seq[j]->state[i] == 'T') &&
                  (data->c_seq[k]->state[i] == 'A' || data->c_seq[k]->state[i] == 'G')))
                {
                  diff++;
                  mat->Q[j][k]+=data->wght[i];
                  len[j][k]+=data->wght[i];
                  len[k][j]=len[j][k];
                }
              
              else
                if(((data->c_seq[j]->state[i] == 'A' && data->c_seq[k]->state[i] == 'G') ||
                    (data->c_seq[j]->state[i] == 'G' && data->c_seq[k]->state[i] == 'A'))||
                   ((data->c_seq[j]->state[i] == 'C' && data->c_seq[k]->state[i] == 'T') ||
                    (data->c_seq[j]->state[i] == 'T' && data->c_seq[k]->state[i] == 'C')))
                  {
                    diff++;
                    mat->P[j][k]+=data->wght[i];
                    len[j][k]+=data->wght[i];
                    len[k][j]=len[j][k];
                  }
                else
                  if((data->c_seq[j]->state[i] == 'A' ||
                      data->c_seq[j]->state[i] == 'C' ||
                      data->c_seq[j]->state[i] == 'G' ||
                      data->c_seq[j]->state[i] == 'T')&&
                     (data->c_seq[k]->state[i] == 'A' ||
                      data->c_seq[k]->state[i] == 'C' ||
                      data->c_seq[k]->state[i] == 'G' ||
                      data->c_seq[k]->state[i] == 'T'))
                    {
                      len[j][k]+=data->wght[i];
                      len[k][j]=len[j][k];
                    }
            }
        }
    }
  
  
  for(i=0;i<data->n_otu-1;i++)
    for(j=i+1;j<data->n_otu;j++)
      {
    if(len[i][j] > .0)
      {
        mat->P[i][j] /= len[i][j];
        mat->Q[i][j] /= len[i][j];
      }
    else
      {
        mat->P[i][j] = .5;
        mat->Q[i][j] = .5;
      }

    mat->P[j][i] = mat->P[i][j];
    mat->Q[j][i] = mat->Q[i][j];


    if((1-2*mat->P[i][j]-mat->Q[i][j] <= .0) || (1-2*mat->Q[i][j] <= .0))
      {
        mat->dist[i][j] = -1.;
        mat->dist[j][i] = -1.;
        continue;
      }

    mat->dist[i][j] = (g_shape/2)*
      (POW(1-2*mat->P[i][j]-mat->Q[i][j],-1./g_shape) +
       0.5*POW(1-2*mat->Q[i][j],-1./g_shape) - 1.5);

    if(mat->dist[i][j] > DIST_MAX) mat->dist[i][j] = DIST_MAX;

    mat->dist[j][i] = mat->dist[i][j];
      }

  for(i=0;i<data->n_otu;i++) free(len[i]);
  free(len);
  return mat;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


matrix *JC69_Dist(calign *data, t_mod *mod)
{
  int site,i,j,k;
  matrix *mat;
  phydbl **len;
  int datatype;


  len = (phydbl **)mCalloc(data->n_otu,sizeof(phydbl *));
  for(i=0;i<data->n_otu;i++)
    len[i] = (phydbl *)mCalloc(data->n_otu,sizeof(phydbl));

  mat = Make_Mat(data->n_otu);
  Init_Mat(mat,data);

  datatype = mod->io->datatype;

  For(site,data->c_seq[0]->len)
    {
      for(j=0;j<data->n_otu-1;j++)
    {
      for(k=j+1;k<data->n_otu;k++)
        {
          if((!Is_Ambigu(data->c_seq[j]->state+site*mod->io->state_len,datatype,mod->io->state_len)) &&
         (!Is_Ambigu(data->c_seq[k]->state+site*mod->io->state_len,datatype,mod->io->state_len)))
        {
          len[j][k]+=data->wght[site];
          len[k][j]=len[j][k];


          if(strncmp(data->c_seq[j]->state+site*mod->io->state_len,
                 data->c_seq[k]->state+site*mod->io->state_len,mod->io->state_len))
/* 		  if(!Are_Compatible(data->c_seq[j]->state+site*mod->io->state_len, */
/* 				     data->c_seq[k]->state+site*mod->io->state_len, */
/* 				     mod->io->state_len, */
/* 				     mod->io->datatype)) */
            mat->P[j][k]+=data->wght[site];
        }
        }
    }
    }


  for(i=0;i<data->n_otu-1;i++)
    for(j=i+1;j<data->n_otu;j++)
      {
    if(len[i][j] > .0) mat->P[i][j] /= len[i][j];
    else               mat->P[i][j] = 1.;

    mat->P[j][i] = mat->P[i][j];

    if((1.-(mod->ns)/(mod->ns-1.)*mat->P[i][j]) < .0) mat->dist[i][j] = -1.;
    else
      mat->dist[i][j] = -(mod->ns-1.)/(mod->ns)*(phydbl)log(1.-(mod->ns)/(mod->ns-1.)*mat->P[i][j]);

/* 	PhyML_Printf("\n. Incorrect JC distances"); */
/* 	mat->dist[i][j] = len[i][j]; */

    if(mat->dist[i][j] > DIST_MAX) mat->dist[i][j] = DIST_MAX;

    mat->dist[j][i] = mat->dist[i][j];
      }

  for(i=0;i<data->n_otu;i++) free(len[i]);
  free(len);

  return mat;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


matrix *Hamming_Dist(calign *data, t_mod *mod)
{
  int i,j,k;
  matrix *mat;
  phydbl **len;
  int datatype;

  len = (phydbl **)mCalloc(data->n_otu,sizeof(phydbl *));
  for(i=0;i<data->n_otu;i++)
    len[i] = (phydbl *)mCalloc(data->n_otu,sizeof(phydbl));

  mat = Make_Mat(data->n_otu);
  Init_Mat(mat,data);

  datatype = mod->io->datatype;

  for(i=0;i<data->crunch_len;i++)
    {
      for(j=0;j<data->n_otu-1;j++)
    {
      for(k=j+1;k<data->n_otu;k++)
        {
          if((!Is_Ambigu(data->c_seq[j]->state+i*mod->io->state_len,datatype,mod->io->state_len)) &&
         (!Is_Ambigu(data->c_seq[k]->state+i*mod->io->state_len,datatype,mod->io->state_len)))
        {
          len[j][k]+=data->wght[i];
          len[k][j]=len[j][k];
/* 		  if(data->c_seq[j]->state[i] != data->c_seq[k]->state[i]) */
          if(!Are_Compatible(data->c_seq[j]->state+i*mod->io->state_len,
                     data->c_seq[k]->state+i*mod->io->state_len,
                     mod->io->state_len,
                     mod->io->datatype))
            {
              mat->P[j][k]+=data->wght[i];
            }
        }
        }
    }
    }

  for(i=0;i<data->n_otu-1;i++)
    for(j=i+1;j<data->n_otu;j++)
      {
    if(len[i][j] > .0)
      {
        mat->P[i][j] /= len[i][j];
      }
    else
      {
        mat->P[i][j] = 1.;
      }

    mat->P[j][i] = mat->P[i][j];

    mat->dist[i][j] = mat->P[i][j];


    if(mat->dist[i][j] > DIST_MAX)
      {
        mat->dist[i][j] = DIST_MAX;
      }
    mat->dist[j][i] = mat->dist[i][j];
      }

  for(i=0;i<data->n_otu;i++) free(len[i]);
  free(len);

  return mat;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* Test if the given site pattern is invariant. Does not handle ambiguities */

int Is_Invar(int patt_num, int stepsize, int datatype, calign *data)
{
  int i, j;

  for(i=0;i<data->n_otu;i++)
    {
      for(j=0;j<data->n_otu;j++)
    {
      if(!(Are_Compatible(data->c_seq[i]->state+patt_num,
                  data->c_seq[j]->state+patt_num,
                  stepsize,
                  datatype)))
        {
          break;
        }
    }
      if(j != data->n_otu) break;
    }

  if(i == data->n_otu) return 1;
  else                 return 0;
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int Is_Ambigu(char *state, int datatype, int stepsize)
{
  int val,i;

  val = -1;
  if(datatype == NT)
    {
      for(i=0;i<stepsize;i++)
        {
          switch(state[i])
            {
            case 'A' : case 'C' : case 'G' : case 'T' : case 'U' : { val=NO; break; }
            default : { val=YES; break; }
            }
          if(val == YES) break;
        }
    }
  else if(datatype == AA)
    {
      switch(state[0])
        {
        case 'X' : case '?' : case '-' : case '.' : {val=YES; break; }
        default : { val=NO; break; }
        }
    }
  else if(datatype == GENERIC)
    {
      int i;
      for(i=0;i<stepsize;i++) if(!isdigit(state[i])) break;
      if(i == stepsize) val = NO;
      else              val = YES;
    }
  
  return val;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Check_Ambiguities(calign *data, int datatype, int stepsize)
{
  int i,j;

  for(j=0;j<data->crunch_len;j++)
    {
      data->ambigu[j] = NO;
      for(i=0;i<data->n_otu;i++)
        {
          data->c_seq[i]->is_ambigu[j] = NO;
        }
      
      for(i=0;i<data->n_otu;i++)
        {
          if(Is_Ambigu(data->c_seq[i]->state+j*stepsize,
                       datatype,
                       stepsize))
            {
              data->ambigu[j]              = YES;
              data->c_seq[i]->is_ambigu[j] = YES;
            }
        }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Set_D_States(calign *data, int datatype, int stepsize)
{
  int i,j;

  for(j=0;j<data->crunch_len;j++)
    {
      for(i=0;i<data->n_otu;i++)
        {
          if(data->c_seq[i]->is_ambigu[j] == NO)
            {
              data->c_seq[i]->d_state[j] = Assign_State(data->c_seq[i]->state+j,
                                                        datatype,
                                                        stepsize);
            }
        }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int Get_State_From_Ui(int ui, int datatype)
{
  if(datatype == NT)
    {
      switch(ui)
    {
    case 1 : {return 0; break;}
    case 2 : {return 1; break;}
    case 4 : {return 2; break;}
    case 8 : {return 3; break;}
    default :
      {
        PhyML_Fprintf(stderr,"\n. ui=%d",ui);
        PhyML_Fprintf(stderr,"\n. Err in file %s at line %d\n",__FILE__,__LINE__);
        Warn_And_Exit("\n. PhyML finished prematurely.");
        break;
      }
    }
    }
  else if(datatype == AA)
    {
      switch(ui)
    {
    case 1 :      {return 0;  break;}
    case 2 :      {return 1;  break;}
    case 4 :      {return 2;  break;}
    case 8 :      {return 3;  break;}
    case 16 :     {return 4;  break;}
    case 32 :     {return 5;  break;}
    case 64 :     {return 6;  break;}
    case 128 :    {return 7;  break;}
    case 256 :    {return 8;  break;}
    case 512 :    {return 9;  break;}
    case 1024 :   {return 10; break;}
    case 2048 :   {return 11; break;}
    case 4096 :   {return 12; break;}
    case 8192 :   {return 13; break;}
    case 16384 :  {return 14; break;}
    case 32768 :  {return 15; break;}
    case 65536 :  {return 16; break;}
    case 131072 : {return 17; break;}
    case 262144 : {return 18; break;}
    case 524288 : {return 19; break;}
    default :
      {
        PhyML_Fprintf(stderr,"\n. ui=%d",ui);
        PhyML_Fprintf(stderr,"\n. Err in file %s at line %d\n",__FILE__,__LINE__);
        Warn_And_Exit("\n. PhyML finished prematurely.");
      }
    }
    }
  else Generic_Exit(__FILE__,__LINE__,__FUNCTION__);    
  return -1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


int Assign_State(char *c, int datatype, int stepsize)
{
  int state[3];
  int i;

  state[0] = state[1] = state[2] = -1;
  if(datatype == NT)
    {
      for(i=0;i<stepsize;i++)
        {
          switch(c[i])
            {
            case 'A' : {state[i]=0;  break;}
            case 'C' : {state[i]=1;  break;}
            case 'G' : {state[i]=2;  break;}
            case 'T' : {state[i]=3;  break;}
            case 'U' : {state[i]=3;  break;}
            default  : {state[i]=-1; break;}
            }
        }
      return (stepsize>1)?(state[0]*16+state[1]*4+state[2]):(state[0]);
    }
  else if(datatype == AA)
    {
      switch(c[0])
    {
    case 'A' : {state[0]=0 ; break;}
    case 'R' : {state[0]=1 ; break;}
    case 'N' : {state[0]=2 ; break;}
    case 'D' : {state[0]=3 ; break;}
    case 'C' : {state[0]=4 ; break;}
    case 'Q' : {state[0]=5 ; break;}
    case 'E' : {state[0]=6 ; break;}
    case 'G' : {state[0]=7 ; break;}
    case 'H' : {state[0]=8 ; break;}
    case 'I' : {state[0]=9 ; break;}
    case 'L' : {state[0]=10; break;}
    case 'K' : {state[0]=11; break;}
    case 'M' : {state[0]=12; break;}
    case 'F' : {state[0]=13; break;}
    case 'P' : {state[0]=14; break;}
    case 'S' : {state[0]=15; break;}
    case 'T' : {state[0]=16; break;}
    case 'W' : {state[0]=17; break;}
    case 'Y' : {state[0]=18; break;}
    case 'V' : {state[0]=19; break;}

    case 'B' : {state[0] = 2; break;}
    case 'Z' : {state[0] = 5; break;}
    default  : {state[0]=-1;  break;}
    }
      return state[0];
    }
  else if(datatype == GENERIC)
    {
      char format[6];
      int ret;

      sprintf(format,"%%%dd",stepsize);
      ret = sscanf(c,format,state);
      if(!ret) state[0] = -1;
      return state[0];
    }
  else
    {
      PhyML_Printf("\n. Not implemented yet.\n");
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);    
    }

  return -1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

char Reciproc_Assign_State(int i_state, int datatype)
{
  if(datatype == NT)
    {
      i_state = i_state%4;
      switch(i_state)
        {
        case 0 :   {return 'A';  break;}
        case 1 :   {return 'C';  break;}
        case 2 :   {return 'G';  break;}
        case 3 :   {return 'T';  break;}
        default  :
          {
            PhyML_Printf("\n. i_state = %d",i_state);
            PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
            Warn_And_Exit("\n. PhyML finished prematurely.");
            break;
          }
        }
    }
  else if(datatype == AA)
    {
      i_state = i_state%20;
      switch(i_state)
        {
        case 0  : {return 'A' ; break;}
        case 1  : {return 'R' ; break;}
        case 2  : {return 'N' ; break;}
        case 3  : {return 'D' ; break;}
        case 4  : {return 'C' ; break;}
        case 5  : {return 'Q' ; break;}
        case 6  : {return 'E' ; break;}
        case 7  : {return 'G' ; break;}
        case 8  : {return 'H' ; break;}
        case 9  : {return 'I' ; break;}
        case 10 : {return 'L';  break;}
        case 11 : {return 'K';  break;}
        case 12 : {return 'M';  break;}
        case 13 : {return 'F';  break;}
        case 14 : {return 'P';  break;}
        case 15 : {return 'S';  break;}
        case 16 : {return 'T';  break;}
        case 17 : {return 'W';  break;}
        case 18 : {return 'Y';  break;}
        case 19 : {return 'V';  break;}
        default  :
          {
            PhyML_Printf("\n. Err. in file %s at line %d\n",__FILE__,__LINE__);
            Warn_And_Exit("\n. PhyML finished prematurely.");
            break;
          }
        }
    }
  else if(datatype == GENERIC)
    {
      return (char)i_state;
    }
  return -1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


int Assign_State_With_Ambiguity(char *c, int datatype, int stepsize)
{
  int state[3];
  int i;

  state[0] = state[1] = state[2] = -1;
  if(datatype == NT)
    {
      for(i=0;i<stepsize;i++)
        {
          switch(c[i])
            {
            case 'A' : {state[i]= 0;  break;}
            case 'C' : {state[i]= 1;  break;}
            case 'G' : {state[i]= 2;  break;}
            case 'T' : {state[i]= 3;  break;}
            case 'U' : {state[i]= 3;  break;}
            case 'M' : {state[i]= 4;  break;}
            case 'R' : {state[i]= 5;  break;}
            case 'W' : {state[i]= 6;  break;}
            case 'S' : {state[i]= 7;  break;}
            case 'Y' : {state[i]= 8;  break;}
            case 'K' : {state[i]= 9;  break;}
            case 'B' : {state[i]=10;  break;}
            case 'D' : {state[i]=11;  break;}
            case 'H' : {state[i]=12;  break;}
            case 'V' : {state[i]=13;  break;}
            case 'N' : case 'X' : case '?' : case 'O' : case '-' : {state[i]=T_MAX_ALPHABET-1;  break;}
            default :
              {
                PhyML_Printf("\n. Unknown character state : '%c'\n",c[i]);
                Warn_And_Exit("\n. Init failed (data type supposed to be DNA)\n");
                break;
              }
            }
          return (stepsize>1)?(state[0]*16+state[1]*4+state[2]):(state[0]);
        }
    }
  else if(datatype == AA)
    {
      switch(c[0])
        {
        case 'A' : {state[0]= 0; break;}
        case 'R' : {state[0]= 1; break;}
        case 'N' : {state[0]= 2; break;}
        case 'D' : {state[0]= 3; break;}
        case 'C' : {state[0]= 4; break;}
        case 'Q' : {state[0]= 5; break;}
        case 'E' : {state[0]= 6; break;}
        case 'G' : {state[0]= 7; break;}
        case 'H' : {state[0]= 8; break;}
        case 'I' : {state[0]= 9; break;}
        case 'L' : {state[0]=10; break;}
        case 'K' : {state[0]=11; break;}
        case 'M' : {state[0]=12; break;}
        case 'F' : {state[0]=13; break;}
        case 'P' : {state[0]=14; break;}
        case 'S' : {state[0]=15; break;}
        case 'T' : {state[0]=16; break;}
        case 'W' : {state[0]=17; break;}
        case 'Y' : {state[0]=18; break;}
        case 'V' : {state[0]=19; break;}
        case 'B' : {state[0]= 2; break;}
        case 'Z' : {state[0]= 5; break;}
        case 'X' : case '?' : case '-' : {state[0]=T_MAX_ALPHABET-1; break;}
        default  :
          {
            PhyML_Printf("\n. Unknown character state : %c\n",state[0]);
            Warn_And_Exit("\n. Init failed (data type supposed to be amino-acids)\n");
            break;
          }
        }
      return state[0];
    }
  else if(datatype == GENERIC)
    {
      if(Is_Ambigu(c,GENERIC,stepsize)) state[0] = T_MAX_ALPHABET-1;
      else
        {
          char format[6];
          sprintf(format,"%%%dd",stepsize);
          if(!sscanf(c,format,state))
            {
              PhyML_Printf("\n. Error reading character. Was expecting an integer, got '%c' instead.\n",c[0]);
              PhyML_Printf("\n. Err. in file %s at line %d (function '%s')\n",__FILE__,__LINE__,__FUNCTION__);
              Warn_And_Exit("\n. PhyML finished prematurely.");
            }
        }
      return state[0];
    }
  
  return -1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Clean_Tree_Connections(t_tree *tree)
{

  int i;
  For(i,2*tree->n_otu-2)
    {
      tree->a_nodes[i]->v[0] = NULL;
      tree->a_nodes[i]->v[1] = NULL;
      tree->a_nodes[i]->v[2] = NULL;
      tree->a_nodes[i]->b[0] = NULL;
      tree->a_nodes[i]->b[1] = NULL;
      tree->a_nodes[i]->b[2] = NULL;
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Bootstrap(t_tree *tree)
{
  int *site_num, n_site;
  int replicate,j,k;
  int position,init_len;
  calign *boot_data;
  t_tree *boot_tree;
  t_mod *boot_mod;
  matrix *boot_mat;
  char *s;
/*   phydbl rf; */


  if(tree->is_mixt_tree == YES)
    {
      PhyML_Printf("\n. Bootstrap option not yet available for partition/mixture analysis.");
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);    
    }

  tree->print_boot_val = 1;
  tree->print_alrt_val = 0;
  boot_tree            = NULL;

  site_num = (int *)mCalloc(tree->data->init_len,sizeof(int));

  Free_Bip(tree);
  Alloc_Bip(tree);
  Get_Bip(tree->a_nodes[0],tree->a_nodes[0]->v[0],tree);

  n_site = 0;
  for(j=0;j<tree->data->crunch_len;j++) For(k,tree->data->wght[j])
    {
      site_num[n_site] = j;
      n_site++;
    }

  boot_data = Copy_Cseq(tree->data,tree->io);

  PhyML_Printf("\n\n. Non parametric bootstrap analysis \n\n");
  PhyML_Printf("  [");

  for(replicate=0;replicate<tree->mod->bootstrap;replicate++)
    {
      for(j=0;j<boot_data->crunch_len;j++) boot_data->wght[j] = 0;

      init_len = 0;
      for(j=0;j<boot_data->init_len;j++)
        {
          position = Rand_Int(0,(int)(tree->data->init_len-1.0));
          boot_data->wght[site_num[position]] += 1;
          init_len++;
        }

      if(init_len != tree->data->init_len) Exit("\n. Pb. when copying sequences\n");

      init_len = 0;
      for(j=0;j<boot_data->crunch_len;j++) init_len += boot_data->wght[j];

      if(init_len != tree->data->init_len) Exit("\n. Pb. when copying sequences\n");

      if(tree->io->datatype == NT)      Get_Base_Freqs(boot_data);
      else if(tree->io->datatype == AA) Get_AA_Freqs(boot_data);

      if(tree->io->random_boot_seq_order) Randomize_Sequence_Order(boot_data);

      Set_D_States(boot_data,tree->io->datatype,tree->io->state_len);

      boot_mod        = Copy_Model(tree->mod);

      boot_mod->s_opt = tree->mod->s_opt; /* WARNING: re-using the same address here instead of creating a copying
                                             requires to leave the value of s_opt unchanged during the boostrap. */
      boot_mod->io    = tree->io; /* WARNING: re-using the same address here instead of creating a copying
                                     requires to leave the value of io unchanged during the boostrap. */

      Init_Model(boot_data,boot_mod,tree->io);

      if(tree->io->mod->use_m4mod) M4_Init_Model(boot_mod->m4mod,boot_data,boot_mod);

      if(tree->io->in_tree == 2)
        {
          rewind(tree->io->fp_in_tree);
          boot_tree = Read_Tree_File_Phylip(tree->io->fp_in_tree);
        }
      else
        {
          boot_mat = ML_Dist(boot_data,boot_mod);
          boot_mat->tree = Make_Tree_From_Scratch(boot_data->n_otu,boot_data);
          Fill_Missing_Dist(boot_mat);
          Bionj(boot_mat);
          boot_tree = boot_mat->tree;
          boot_tree->mat = boot_mat;
        }

      boot_tree->mod                  = boot_mod;
      boot_tree->io                   = tree->io;
      boot_tree->data                 = boot_data;
      boot_tree->verbose              = VL0;
      /* boot_tree->verbose              = VL3; */
      boot_tree->n_pattern            = boot_tree->data->crunch_len;
      boot_tree->io->print_site_lnl   = NO;
      boot_tree->io->print_trace      = NO;
      boot_tree->io->print_json_trace = NO;
      boot_tree->n_root               = NULL;
      boot_tree->e_root               = NULL;

      Set_Both_Sides(YES,boot_tree);

      if((boot_tree->mod->s_opt->random_input_tree) && (boot_tree->mod->s_opt->topo_search == SPR_MOVE)) Random_Tree(boot_tree);

      Connect_CSeqs_To_Nodes(boot_data,tree->io,boot_tree);
      Check_Br_Lens(boot_tree);
      Share_Lk_Struct(tree,boot_tree);
      Share_Spr_Struct(tree,boot_tree);
      Share_Pars_Struct(tree,boot_tree);
      Update_Dirs(boot_tree);

      Init_Partial_Lk_Tips_Double(boot_tree);
      Init_Ui_Tips(boot_tree);
      Init_Partial_Pars_Tips(boot_tree);
      Br_Len_Not_Involving_Invar(boot_tree);

      if(boot_tree->io->do_alias_subpatt)
        {
          MIXT_Set_Alias_Subpatt(YES,boot_tree);
          Lk(NULL,boot_tree);
          MIXT_Set_Alias_Subpatt(NO,boot_tree);
        }

      Switch_Eigen(YES,boot_tree->mod);                 
      Lk(NULL,boot_tree);
      Switch_Eigen(NO,boot_tree->mod);

      if(boot_tree->mod->s_opt->opt_topo)
        {
          if(boot_tree->mod->s_opt->topo_search == NNI_MOVE)
            {
              Simu_Loop(boot_tree);
            }
          else if((boot_tree->mod->s_opt->topo_search == SPR_MOVE) ||
                  (boot_tree->mod->s_opt->topo_search == BEST_OF_NNI_AND_SPR))
            {
              Speed_Spr_Loop(boot_tree);
            }
        }
      else
        {
          if(boot_tree->mod->s_opt->opt_subst_param || boot_tree->mod->s_opt->opt_bl)
            Round_Optimize(boot_tree,ROUND_MAX);
          else
            Lk(NULL,boot_tree);
        }

      Free_Bip(boot_tree);

      Alloc_Bip(boot_tree);

      Match_Tip_Numbers(tree,boot_tree);

      Get_Bip(boot_tree->a_nodes[0],
              boot_tree->a_nodes[0]->v[0],
              boot_tree);

      Compare_Bip(tree,boot_tree,NO);

      Check_Br_Lens(boot_tree);
      Br_Len_Involving_Invar(boot_tree);

      if(tree->io->print_boot_trees)
        {
          s = Write_Tree(boot_tree,NO);
          PhyML_Fprintf(tree->io->fp_out_boot_tree,"%s\n",s);
          Free(s);
          Print_Fp_Out_Lines(tree->io->fp_out_boot_stats,0,0,boot_tree,tree->io,replicate+1);
        }

      /*       rf = .0; */
      /*       For(j,2*tree->n_otu-3)  */
      /* 	rf += tree->a_edges[j]->bip_score; */


      PhyML_Printf(".");
#ifndef QUIET
      fflush(stdout);
#endif
      if(!((replicate+1)%tree->io->boot_prog_every))
        {
          PhyML_Printf("] %4d/%4d\n  ",replicate+1,tree->mod->bootstrap);
          if(replicate != tree->mod->bootstrap-1) PhyML_Printf("[");
        }

      Free_Tree(boot_tree);
      Free_Model(boot_mod);
    }

  if(((replicate)%tree->io->boot_prog_every)) PhyML_Printf("] %4d/%4d\n ",replicate,tree->mod->bootstrap);

  tree->lock_topo = YES; /* Topology should not be modified afterwards */

  if(tree->io->print_boot_trees)
    {
      fclose(tree->io->fp_out_boot_tree);
      fclose(tree->io->fp_out_boot_stats);
    }

  Free_Calign(boot_data);
  Free(site_num);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Br_Len_Involving_Invar(t_tree *tree)
{
  int i;

  if(tree->is_mixt_tree)
    {
      MIXT_Br_Len_Involving_Invar(tree);
      return;
    }

  For(i,2*tree->n_otu-1) tree->a_edges[i]->l->v *= (1.0-tree->mod->ras->pinvar->v);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Br_Len_Not_Involving_Invar(t_tree *tree)
{
  int i;

  if(tree->is_mixt_tree)
    {
      MIXT_Br_Len_Not_Involving_Invar(tree);
      return;
    }

  For(i,2*tree->n_otu-1) tree->a_edges[i]->l->v /= (1.0-tree->mod->ras->pinvar->v);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Getstring_Stdin(char *s)
{
  if(!fgets(s,T_MAX_LINE,stdin)) Exit("");
  if (strchr(s, '\n') != NULL)
    *strchr(s, '\n') = '\0';
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl Num_Derivatives_One_Param(phydbl (*func)(t_tree *tree), t_tree *tree,
                                 phydbl f0, phydbl *param, int which, int n_param, phydbl stepsize, int logt,
                                 phydbl *err, int precise, int is_positive)
{
  int i,j;
  phydbl errt,fac,hh,**a,ans,*sign;
  int n_iter;

  sign = (phydbl *)mCalloc(n_param,sizeof(phydbl));

  a = (phydbl **)mCalloc(11,sizeof(phydbl *));
  for(i=0;i<11;i++) a[i] = (phydbl *)mCalloc(11,sizeof(phydbl));

  n_iter = 10; /* */

  ans  = .0;

  if(stepsize < SMALL) Warn_And_Exit("\n. h must be nonzero in Dfridr.");

  hh=stepsize;

  if(!precise)
    {
      param[which]   = param[which]+hh;

      if(logt == YES) for(i=0;i<n_param;i++) param[i] = exp(MIN(1.E+2,param[i]));
      for(i=0;i<n_param;i++) sign[i] = param[i] > .0 ? 1. : -1.;
      if(is_positive == YES) for(i=0;i<n_param;i++) param[i] = FABS(param[i]);
      a[0][0]  = (*func)(tree);
      if(is_positive == YES) for(i=0;i<n_param;i++) param[i] *= sign[i];
      if(logt == YES) for(i=0;i<n_param;i++) param[i] = log(param[i]);

      /* printf("\n. f0=%f f1=%f hh=%G %f",f0,a[0][0],hh,param[which]); */

      a[0][0]  -= f0;
      a[0][0]  /= hh;
      param[which]   = param[which]-hh;

      ans =  a[0][0];
    }
  else
    {
      param[which]   = param[which]+hh;

      if(logt == YES) for(i=0;i<n_param;i++) param[i] = exp(MIN(1.E+2,param[i]));
      for(i=0;i<n_param;i++) sign[i] = param[i] > .0 ? 1. : -1.;
      if(is_positive == YES) for(i=0;i<n_param;i++) param[i] = FABS(param[i]);
      a[0][0]  = (*func)(tree);
      if(is_positive == YES) for(i=0;i<n_param;i++) param[i] *= sign[i];
      if(logt == YES) for(i=0;i<n_param;i++) param[i] = log(param[i]);

      /*   param[which]   = param[which]-2*hh; */
      /*   a[0][0] -= (*func)(tree); */
      /*   a[0][0] /= (2.0*hh); */
      /*   param[which]   = param[which]+hh; */
      a[0][0]  -= f0;
      a[0][0]  /= hh;
      param[which]   = param[which]-hh;

      *err=1e30;
      for(i=1;i<n_iter;i++)
        {
          hh /= 1.4;
          
          /*       param[which]   = param[which]+hh; */
          /*       a[0][i]  = (*func)(tree); */
          /*       param[which]   = param[which]-2*hh; */
          /*       a[0][i] -= (*func)(tree); */
          /*       a[0][i] /= (2.0*hh); */
          /*       param[which]   = param[which]+hh; */
          
          param[which]   = param[which]+hh;
          
          if(logt == YES) for(j=0;j<n_param;j++) param[j] = exp(MIN(1.E+2,param[j]));
          for(i=0;i<n_param;i++) sign[i] = param[i] > .0 ? 1. : -1.;
          if(is_positive == YES) for(i=0;i<n_param;i++) param[i] = FABS(param[i]);
          a[0][i]  = (*func)(tree);
          if(is_positive == YES) for(i=0;i<n_param;i++) param[i] *= sign[i];
          if(logt == YES) for(j=0;j<n_param;j++) param[j] = log(param[j]);
          
          /*   param[which]   = param[which]-2*hh; */
          /*   a[0][i] -= (*func)(tree); */
          /*   a[0][i] /= (2.0*hh); */
          /*   param[which]   = param[which]+hh; */
          a[0][i]  -= f0;
          a[0][i]  /= hh;
          param[which]   = param[which]-hh;
          
          
          fac=1.4*1.4;
          for (j=1;j<=i;j++)
            {
              a[j][i]=(a[j-1][i]*fac-a[j-1][i-1])/(fac-1.0);
              fac=1.4*1.4*fac;
              
              errt=MAX(FABS(a[j][i]-a[j-1][i]),FABS(a[j][i]-a[j-1][i-1]));
              
              if (errt <= *err)
                {
                  *err=errt;
                  ans=a[j][i];
                }
            }
          
          if(FABS(a[i][i]-a[i-1][i-1]) >= 2.0*(*err)) break;
        }
    }
  for(i=0;i<11;i++) Free(a[i]);
  Free(a);
  Free(sign);
  
  return ans;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl Num_Derivatives_One_Param_Nonaligned(phydbl (*func)(t_tree *tree), t_tree *tree,
                                            phydbl f0, phydbl **param, int which, int n_param, phydbl stepsize, int logt,
                                            phydbl *err, int precise, int is_positive)
{
  int i,j;
  phydbl errt,fac,hh,**a,ans,*sign;
  int n_iter;

  sign = (phydbl *)mCalloc(n_param,sizeof(phydbl));

  a = (phydbl **)mCalloc(11,sizeof(phydbl *));
  for(i=0;i<11;i++) a[i] = (phydbl *)mCalloc(11,sizeof(phydbl));

  n_iter = 10; /* */

  ans  = .0;

  if(stepsize < SMALL) Warn_And_Exit("\n. h must be nonzero in Dfridr.");

  hh=stepsize;

  if(!precise)
    {
      *(param[which])   = *(param[which])+hh;

      if(logt == YES) for(i=0;i<n_param;i++) *(param[i]) = exp(MIN(1.E+2,*(param[i])));
      for(i=0;i<n_param;i++) sign[i] = (*(param[i])) > .0 ? 1. : -1.;
      if(is_positive == YES) for(i=0;i<n_param;i++) (*(param[i])) = FABS(*(param[i]));
      a[0][0]  = (*func)(tree);
      if(is_positive == YES) for(i=0;i<n_param;i++) (*(param[i])) *= sign[i];
      if(logt == YES) for(i=0;i<n_param;i++) *(param[i]) = log(*(param[i]));

      /* printf("\n. f0=%f f1=%f hh=%G %f",f0,a[0][0],hh,*(param[which])); */

      a[0][0]  -= f0;
      a[0][0]  /= hh;
      *(param[which]) = *(param[which])-hh;

      ans =  a[0][0];
    }
  else
    {
      *(param[which])   = *(param[which])+hh;

      if(logt == YES) for(i=0;i<n_param;i++) *(param[i]) = exp(MIN(1.E+2,*(param[i])));
      for(i=0;i<n_param;i++) sign[i] = (*(param[i])) > .0 ? 1. : -1.;
      if(is_positive == YES) for(i=0;i<n_param;i++) (*(param[i])) = FABS(*(param[i]));
      a[0][0]  = (*func)(tree);
      if(is_positive == YES) for(i=0;i<n_param;i++) (*(param[i])) *= sign[i];
      if(logt == YES) for(i=0;i<n_param;i++) *(param[i]) = log(*(param[i]));

      /*   *(param[which]   = *(param[which]-2*hh; */
      /*   a[0][0] -= (*func)(tree); */
      /*   a[0][0] /= (2.0*hh); */
      /*   *(param[which]   = *(param[which]+hh; */
      a[0][0]  -= f0;
      a[0][0]  /= hh;
      *(param[which])   = *(param[which])-hh;

      *err=1e30;
      for(i=1;i<n_iter;i++)
    {
      hh /= 1.4;

      /*       *(param[which]   = *(param[which]+hh; */
      /*       a[0][i]  = (*func)(tree); */
      /*       *(param[which]   = *(param[which]-2*hh; */
      /*       a[0][i] -= (*func)(tree); */
      /*       a[0][i] /= (2.0*hh); */
      /*       *(param[which]   = *(param[which]+hh; */

      *(param[which])   = *(param[which])+hh;

          if(logt == YES) for(j=0;j<n_param;j++) *(param[j]) = exp(MIN(1.E+2,*(param[j])));
          for(i=0;i<n_param;i++) sign[i] = (*(param[i])) > .0 ? 1. : -1.;
          if(is_positive == YES) for(i=0;i<n_param;i++) (*(param[i])) = FABS(*(param[i]));
      a[0][i]  = (*func)(tree);
          if(is_positive == YES) for(i=0;i<n_param;i++) (*(param[i])) *= sign[i];
          if(logt == YES) for(j=0;j<n_param;j++) *(param[j]) = log(*(param[j]));

      /*   *(param[which]   = *(param[which]-2*hh; */
      /*   a[0][i] -= (*func)(tree); */
      /*   a[0][i] /= (2.0*hh); */
      /*   *(param[which]   = *(param[which]+hh; */
      a[0][i]  -= f0;
      a[0][i]  /= hh;
      *(param[which])   = *(param[which])-hh;


      fac=1.4*1.4;
      for (j=1;j<=i;j++)
        {
          a[j][i]=(a[j-1][i]*fac-a[j-1][i-1])/(fac-1.0);
          fac=1.4*1.4*fac;

          errt=MAX(FABS(a[j][i]-a[j-1][i]),FABS(a[j][i]-a[j-1][i-1]));

          if (errt <= *err)
        {
          *err=errt;
          ans=a[j][i];
        }
        }

      if(FABS(a[i][i]-a[i-1][i-1]) >= 2.0*(*err)) break;
    }
    }
  for(i=0;i<11;i++) Free(a[i]);
  Free(a);
  Free(sign);
  return ans;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int Num_Derivative_Several_Param(t_tree *tree, phydbl *param, int n_param, phydbl stepsize, int logt,
                                 phydbl (*func)(t_tree *tree), phydbl *derivatives, int is_positive)
{
  int i;
  phydbl err,f0,*sign;

  sign = (phydbl *)mCalloc(n_param,sizeof(phydbl));

  if(logt == YES)   for(i=0;i<n_param;i++) param[i] = exp(MIN(1.E+2,param[i]));
  for(i=0;i<n_param;i++) sign[i] = (param[i]) > .0 ? 1. : -1.;
  if(is_positive == YES) for(i=0;i<n_param;i++) param[i] = FABS(param[i]);
  f0 = (*func)(tree);
  if(is_positive == YES) for(i=0;i<n_param;i++) param[i] *= sign[i];
  if(logt == YES)   for(i=0;i<n_param;i++) param[i] = log(param[i]);

  for(i=0;i<n_param;i++)
    {
      /* for(j=0;j<tree->mod->r_mat->n_diff_rr;j++) PhyML_Printf("\n. 00%d %f",i,tree->mod->r_mat->rr_val->v[j]); */
      derivatives[i] = Num_Derivatives_One_Param(func,
                                                 tree,
                                                 f0,
                                                 param,
                                                 i,
                                                 n_param,
                                                 stepsize,
                                                 logt,
                                                 &err,
                                                 0,
                                                 is_positive
                                                 );
    }
  
  Free(sign);

  return 1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int Num_Derivative_Several_Param_Nonaligned(t_tree *tree, phydbl **param, int n_param, phydbl stepsize, int logt,
                                            phydbl (*func)(t_tree *tree), phydbl *derivatives, int is_positive)
{
  int i;
  phydbl err,f0,*sign;

  sign = (phydbl *)mCalloc(n_param,sizeof(phydbl));

  if(logt == YES)   for(i=0;i<n_param;i++) (*(param[i])) = exp(MIN(1.E+2,*(param[i])));
  for(i=0;i<n_param;i++) sign[i] = (*(param[i])) > .0 ? 1. : -1.;
  if(is_positive == YES) for(i=0;i<n_param;i++) *(param[i]) = FABS(*(param[i]));
  f0 = (*func)(tree);
  if(is_positive == YES) for(i=0;i<n_param;i++) *(param[i]) *= sign[i];
  if(logt == YES)   for(i=0;i<n_param;i++) (*(param[i])) = log(*(param[i]));


  for(i=0;i<n_param;i++)
    {

      derivatives[i] = Num_Derivatives_One_Param_Nonaligned(func,
                         tree,
                         f0,
                                                            param,
                                                            i,
                                                            n_param,
                         stepsize,
                                                            logt,
                         &err,
                                                            0,
                                                            is_positive
                         );
    }

  Free(sign);

  return 1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


int Compare_Two_States(char *state1, char *state2, int state_size)
{

  /* 1 the two states are identical */
  /* 0 the two states are different */
  int i;

  for(i=0;i<state_size;i++) if(state1[i] != state2[i]) break;

  return (i==state_size)?(1):(0);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Copy_One_State(char *from, char *to, int state_size)
{
  int i;
  for(i=0;i<state_size;i++) to[i] = from[i];
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Copy_Dist(phydbl **cpy, phydbl **orig, int n)
{
  int i,j;
  for(i=0;i<n;i++) for(j=0;j<n;j++) cpy[i][j] = orig[i][j];
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

t_mod *Copy_Model(t_mod *ori)
{
  t_mod *cpy;

  cpy             = Make_Model_Basic();

  cpy->ns          = ori->ns;
  cpy->ras->n_catg = ori->ras->n_catg;
  cpy->whichmodel  = ori->whichmodel;

  Make_Model_Complete(cpy);
  Record_Model(ori,cpy);
  /* cpy->m4mod = M4_Copy_M4_Model(ori, ori->m4mod); */

#ifdef BEAGLE
  cpy->b_inst              = ori->b_inst;
  cpy->optimizing_topology = ori->optimizing_topology;
#endif

  return cpy;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Record_Model(t_mod *ori, t_mod *cpy)
{
  int i;

  cpy->ns                   = ori->ns;
  cpy->ras->n_catg          = ori->ras->n_catg;
  cpy->ras->normalise_rr    = ori->ras->normalise_rr;
  cpy->l_var_sigma          = ori->l_var_sigma;

  cpy->kappa->v             = ori->kappa->v;
  cpy->ras->alpha->v        = ori->ras->alpha->v;
  cpy->lambda->v            = ori->lambda->v;
  cpy->ras->pinvar->v       = ori->ras->pinvar->v;
  cpy->br_len_mult->v       = ori->br_len_mult->v;

  strcpy(cpy->modelname->s,ori->modelname->s);
  strcpy(cpy->custom_mod_string->s,ori->custom_mod_string->s);

  cpy->mod_num              = ori->mod_num;
  cpy->whichmodel           = ori->whichmodel;
  cpy->update_eigen         = ori->update_eigen;
  cpy->bootstrap            = ori->bootstrap;
  cpy->ras->invar           = ori->ras->invar;
  cpy->r_mat->n_diff_rr     = ori->r_mat->n_diff_rr;
  cpy->l_min                = ori->l_min;
  cpy->l_max                = ori->l_max;
  cpy->log_l                = ori->log_l;
  cpy->ras->free_mixt_rates = ori->ras->free_mixt_rates;
  cpy->ras->gamma_median    = ori->ras->gamma_median;


  if((ori->whichmodel == CUSTOM) || (ori->whichmodel == GTR))
    {
      For(i,ori->ns*(ori->ns-1)/2)
        {
          cpy->r_mat->rr_num->v[i]       = ori->r_mat->rr_num->v[i];
          cpy->r_mat->rr_val->v[i]       = ori->r_mat->rr_val->v[i];
          cpy->r_mat->rr->v[i]           = ori->r_mat->rr->v[i];
          cpy->r_mat->n_rr_per_cat->v[i] = ori->r_mat->n_rr_per_cat->v[i];
        }
    }
  
  for(i=0;i<cpy->ns;i++)
    {
      cpy->e_frq->pi->v[i]          = ori->e_frq->pi->v[i];
      cpy->e_frq->pi_unscaled->v[i] = ori->e_frq->pi_unscaled->v[i];
      cpy->e_frq->user_b_freq->v[i] = ori->e_frq->user_b_freq->v[i];
    }
  
  For(i,cpy->ns*cpy->ns) cpy->r_mat->qmat->v[i] = ori->r_mat->qmat->v[i];
  
  for(i=0;i<cpy->ras->n_catg;i++)
    {
      cpy->ras->gamma_r_proba->v[i]          = ori->ras->gamma_r_proba->v[i];
      cpy->ras->gamma_rr->v[i]               = ori->ras->gamma_rr->v[i];
      cpy->ras->gamma_r_proba_unscaled->v[i] = ori->ras->gamma_r_proba_unscaled->v[i];
      cpy->ras->gamma_rr_unscaled->v[i]      = ori->ras->gamma_rr_unscaled->v[i];
    }
  
  cpy->use_m4mod = ori->use_m4mod;
  
  cpy->eigen->size = ori->eigen->size;
  For(i,2*ori->ns)       cpy->eigen->space[i]       = ori->eigen->space[i];
  For(i,2*ori->ns)       cpy->eigen->space_int[i]   = ori->eigen->space_int[i];
  for(i=0;i<ori->ns;i++)         cpy->eigen->e_val[i]       = ori->eigen->e_val[i];
  for(i=0;i<ori->ns;i++)         cpy->eigen->e_val_im[i]    = ori->eigen->e_val_im[i];
  For(i,ori->ns*ori->ns) cpy->eigen->r_e_vect[i]    = ori->eigen->r_e_vect[i];
  For(i,ori->ns*ori->ns) cpy->eigen->r_e_vect[i]    = ori->eigen->r_e_vect[i];
  For(i,ori->ns*ori->ns) cpy->eigen->r_e_vect_im[i] = ori->eigen->r_e_vect_im[i];
  For(i,ori->ns*ori->ns) cpy->eigen->l_e_vect[i]    = ori->eigen->l_e_vect[i];
  For(i,ori->ns*ori->ns) cpy->eigen->q[i]           = ori->eigen->q[i];
  
#ifdef BEAGLE
  cpy->b_inst              = ori->b_inst;
  cpy->optimizing_topology = ori->optimizing_topology;
#endif
  
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Test_Node_Table_Consistency(t_tree *tree)
{
  int i;

  For(i,2*tree->n_otu-2)
    {
      if(tree->a_nodes[i]->num != i)
        {
          PhyML_Printf("\n. Node table is not consistent with node numbers.");
          Generic_Exit(__FILE__,__LINE__,__FUNCTION__);    
        }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Get_Bip(t_node *a, t_node *d, t_tree *tree)
{
  int i,j;
  t_node *tmp;
  int swapped;

  if(!d || !a || !tree)
    {
      PhyML_Printf("\n. d: %p a: %p tree: %p",d,a,tree);
      PhyML_Printf("\n. Err. in file %s at line %d (function '%s').\n",__FILE__,__LINE__,__FUNCTION__);
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
    }

  if(d->tax)
    {
      if(d->common)
        {
          d->bip_node[0] = (t_node **)mCalloc(1,sizeof(t_node *));
          d->bip_node[0][0] = d;
          d->bip_size[0]    = 1;
          d->bip_size[1]    = -1;
          d->bip_size[2]    = -1;
          
          for(i=0;i<3;i++)
            {
              if(a->v[i] == d)
                {
                  a->bip_size[i] = 0;
                  for(j=0;j<tree->n_otu;j++)
                    {
                      if(strcmp(tree->a_nodes[j]->name,d->name))
                        {
                          a->bip_node[i] = (t_node **)realloc(a->bip_node[i],(a->bip_size[i]+1)*sizeof(t_node *));
                          a->bip_node[i][a->bip_size[i]] = tree->a_nodes[j];
                          a->bip_size[i]++;
                        }
                    }
                  
                  /* Sort bipartition */
                  do
                    {
                      swapped = NO;
                      For(j,a->bip_size[i]-1)
                        {
                          if(a->bip_node[i][j]->num > a->bip_node[i][j+1]->num)
                            {
                              swapped = YES;
                              tmp                 = a->bip_node[i][j];
                              a->bip_node[i][j]   = a->bip_node[i][j+1];
                              a->bip_node[i][j+1] = tmp;
                            }
                        }
                    }while(swapped == YES);
                  

                  break;
                  
                }
            }
        }
      return;
    }
  else
    {
      int k;
      int d_a;
      
      d_a = -1;
      
      for(i=0;i<3;i++)
        {
          if(d->v[i] != a) Get_Bip(d,d->v[i],tree);
          else if(d->v[i] == a) d_a = i;
        }
      
      d->bip_size[d_a] = 0;
      for(i=0;i<3;i++)
        if(d->v[i] != a)
          {
            for(j=0;j<3;j++)
              {
                if(d->v[i]->v[j] == d)
                  {
                    For(k,d->v[i]->bip_size[j])
                      {
                        d->bip_node[d_a] = (t_node **)realloc(d->bip_node[d_a],(d->bip_size[d_a]+1)*sizeof(t_node *));
                        d->bip_node[d_a][d->bip_size[d_a]] = d->v[i]->bip_node[j][k];
                        d->bip_size[d_a]++;
                      }
                    break;
                  }
              }
          }
      
      do
        {
          swapped = NO;
          For(j,d->bip_size[d_a]-1)
            {
              if(d->bip_node[d_a][j]->num > d->bip_node[d_a][j+1]->num)
                {
                  swapped = YES;
                  tmp                   = d->bip_node[d_a][j];
                  d->bip_node[d_a][j]   = d->bip_node[d_a][j+1];
                  d->bip_node[d_a][j+1] = tmp;
                }
            }
        }while(swapped == YES);
      
      
      for(i=0;i<3;i++)
        if(a->v[i] == d)
          {
            a->bip_size[i] = 0;
            for(j=0;j<tree->n_otu;j++)
              {
                For(k,d->bip_size[d_a])
                  {
                    if(d->bip_node[d_a][k] == tree->a_nodes[j])
                      break;
                  }
                
                if((k == d->bip_size[d_a]) && (tree->a_nodes[j]->common))
                  {
                    a->bip_node[i] = (t_node **)realloc(a->bip_node[i],(a->bip_size[i]+1)*sizeof(t_node *));
                    a->bip_node[i][a->bip_size[i]] = tree->a_nodes[j];
                    a->bip_size[i]++;
                  }
              }
            
            do
              {
                swapped = NO;
                For(j,a->bip_size[i]-1)
                  {
                    if(a->bip_node[i][j]->num > a->bip_node[i][j+1]->num)
                      {
                        swapped = YES;
                        tmp                 = a->bip_node[i][j];
                        a->bip_node[i][j]   = a->bip_node[i][j+1];
                        a->bip_node[i][j+1] = tmp;
                      }
                  }
              }while(swapped == YES);
            
            if(a->bip_size[i] != tree->n_otu - d->bip_size[d_a])
              {
                PhyML_Printf("%d %d \n",a->bip_size[i],tree->n_otu - d->bip_size[d_a]);
                Warn_And_Exit("\n. Problem in counting bipartitions \n");
              }
            break;
          }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Alloc_Bip(t_tree *tree)
{
  int i;

  if(tree->has_bip) return;

  tree->has_bip = YES;

  For(i,2*tree->n_otu-2)
    {
      tree->a_nodes[i]->bip_size = (int *)mCalloc(3,sizeof(int));
      tree->a_nodes[i]->bip_node = (t_node ***)mCalloc(3,sizeof(t_node **));
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int Sort_Phydbl_Increase(const void *a, const void *b)
{
  if((*(phydbl *)(a)) <= (*(phydbl *)(b))) return -1;
  else return 1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int Sort_String(const void *a, const void *b)
{
  return(strcmp((*(const char **)(a)), (*(const char **)(b))));
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int Compare_Bip(t_tree *tree1, t_tree *tree2, int on_existing_edges_only)
{
  int i,j,k;
  t_edge *b1,*b2;
/*   char **bip1,**bip2; */
/*   int *bip1,*bip2; */
  t_node **bip1, **bip2;
  int bip_size1, bip_size2, bip_size;
  int different,identical;
  int n_edges;

  /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
  /* WARNING: call Match_Tip_Numbers and Get_Bip before using this function. */
  /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

  if(on_existing_edges_only == YES)
    {
      n_edges = 0;
      For(i,2*tree1->n_otu-3)
        if(tree1->a_edges[i]->does_exist && tree2->a_edges[i]->does_exist) n_edges++;
      n_edges -= tree1->n_otu;
    }
  else
    {
      n_edges = tree1->n_otu-3;
    }


  identical = 0;
  different = 0;
  For(i,2*tree1->n_otu-3)
    {
      b1 = tree1->a_edges[i];
      bip_size1 = MIN(b1->left->bip_size[b1->l_r],b1->rght->bip_size[b1->r_l]);
      
      if(bip_size1 > 1 && ((on_existing_edges_only == YES && b1->does_exist) || (on_existing_edges_only == NO)))
        {
          For(j,2*tree2->n_otu-3)
            {
              b2 = tree2->a_edges[j];
              bip_size2 = MIN(b2->left->bip_size[b2->l_r],b2->rght->bip_size[b2->r_l]);
              
              if(bip_size2 > 1 && ((on_existing_edges_only == YES && b2->does_exist) || (on_existing_edges_only == NO)))
                {
                  if(bip_size1 == bip_size2)
                    {
                      bip_size = bip_size1;
                      
                      if(b1->left->bip_size[b1->l_r] == b1->rght->bip_size[b1->r_l])
                        {
                          /* 			  if(b1->left->bip_name[b1->l_r][0][0] < b1->rght->bip_name[b1->r_l][0][0]) */
                          if(b1->left->bip_node[b1->l_r][0]->num < b1->rght->bip_node[b1->r_l][0]->num)
                            {
                              /* 			      bip1 = b1->left->bip_name[b1->l_r]; */
                              bip1 = b1->left->bip_node[b1->l_r];
                            }
                          else
                            {
                              /* 			      bip1 = b1->rght->bip_name[b1->r_l]; */
                              bip1 = b1->rght->bip_node[b1->r_l];
                            }
                        }
                      else if(b1->left->bip_size[b1->l_r] < b1->rght->bip_size[b1->r_l])
                        {
                          /* 			  bip1 = b1->left->bip_name[b1->l_r]; */
                          bip1 = b1->left->bip_node[b1->l_r];
                        }
                      else
                        {
                          /* 			  bip1 = b1->rght->bip_name[b1->r_l]; */
                          bip1 = b1->rght->bip_node[b1->r_l];
                        }
                      
                      
                      if(b2->left->bip_size[b2->l_r] == b2->rght->bip_size[b2->r_l])
                        {
                          /* 			  if(b2->left->bip_name[b2->l_r][0][0] < b2->rght->bip_name[b2->r_l][0][0]) */
                          if(b2->left->bip_node[b2->l_r][0]->num < b2->rght->bip_node[b2->r_l][0]->num)
                            {
                              /* 			      bip2 = b2->left->bip_name[b2->l_r]; */
                              bip2 = b2->left->bip_node[b2->l_r];
                            }
                          else
                            {
                              /* 			      bip2 = b2->rght->bip_name[b2->r_l]; */
                              bip2 = b2->rght->bip_node[b2->r_l];
                            }
                        }
                      else if(b2->left->bip_size[b2->l_r] < b2->rght->bip_size[b2->r_l])
                        {
                          /* 			  bip2 = b2->left->bip_name[b2->l_r]; */
                          bip2 = b2->left->bip_node[b2->l_r];
                        }
                      else
                        {
                          /* 			  bip2 = b2->rght->bip_name[b2->r_l]; */
                          bip2 = b2->rght->bip_node[b2->r_l];
                        }
                      
                      if(bip_size == 1) Warn_And_Exit("\n. Problem in Compare_Bip\n");
                      
                      for(k=0;k<bip_size;k++)
                        {
                          /* 			  if(strcmp(bip1[k],bip2[k])) break; */
                          if(bip1[k]->num != bip2[k]->num) break;
                        }
                      
                      if(k == bip_size) /* Branches b1 and b2 define the same bipartition */
                        {
                          b1->bip_score++;
                          b2->bip_score++;
                          identical++;
                          goto out;
                        }
                      else
                        {
                          different++; // Bipartitions have identical sizes but distinct elements
                        }
                    }
                  else different++; // Biparition have different sizes
                }
            }
        }
    out: ;
    }

  return n_edges - identical;
  /* return different; */
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* Modifiy the tip numbering in tree2 so that tips in
   tree1 and tree2 corresponding to the same taxon name
   also have the same tip numbering */
void Match_Tip_Numbers(t_tree *tree1, t_tree *tree2)
{
  int i,j;

  if(tree1->n_otu != tree2->n_otu)
    {
      PhyML_Printf("\n. tree1 and tree2 must have the same number of tips.");
      /* Otherwise, if tree2->n_otu < tree->n_otu, then some tips in tree2
         will have a number (->num) that is the same as the number of an
         internal node in this tree */
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);    
    }

  for(i=0;i<tree1->n_otu;i++)
    {
      for(j=0;j<tree2->n_otu;j++)
    {
      if(!strcmp(tree1->a_nodes[i]->name,tree2->a_nodes[j]->name))
        {
          tree2->a_nodes[j]->num = tree1->a_nodes[i]->num;
          break;
        }
    }
    }

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Test_Multiple_Data_Set_Format(option *io)
{
  char *line;

  line = (char *)mCalloc(T_MAX_LINE,sizeof(char));

  io->n_trees = 0;

  while(fgets(line,T_MAX_LINE,io->fp_in_tree)) if(strstr(line,";")) io->n_trees++;

  Free(line);

  if((io->mod->bootstrap > 1) && (io->n_trees > 1))
    Warn_And_Exit("\n. Bootstrap option is not allowed with multiple input trees !\n");

  rewind(io->fp_in_tree);

  return;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


int Are_Compatible(char *statea, char *stateb, int stepsize, int datatype)
{
  int i,j;
  char a,b;

  if(datatype == NT)
    {
      for(i=0;i<stepsize;i++)
        {
          a = statea[i];
          for(j=0;j<stepsize;j++)
            {
              b = stateb[j];
              
              switch(a)
                {
                case 'A':
                  {
                    switch(b)
                      {
                      case 'A' :
                      case 'M' :
                      case 'R' :
                      case 'W' :
                      case 'D' :
                      case 'H' :
                      case 'V' :
                      case 'X' : {break;}
                      default : return 0;
                      }
                    break;
                  }
                case 'G':
                  {
                    switch(b)
                      {
                      case 'G' :
                      case 'R' :
                      case 'S' :
                      case 'K' :
                      case 'B' :
                      case 'D' :
                      case 'V' :
                      case 'X' : {break;}
                      default : return 0;
                      }
                    break;
                  }
                case 'C':
                  {
                    switch(b)
                      {
                      case 'C' :
                      case 'M' :
                      case 'S' :
                      case 'Y' :
                      case 'B' :
                      case 'H' :
                      case 'V' :
                      case 'X' : { break;}
                      default : return 0;
                      }
                    break;
                  }
                case 'T':
                  {
                    switch(b)
                      {
                      case 'T' :
                      case 'W' :
                      case 'Y' :
                      case 'K' :
                      case 'B' :
                      case 'D' :
                      case 'H' :
                      case 'X' : { break;}
                      default : return 0;
                      }
                    break;
                  }
                case 'M' :
                  {
                    switch(b)
                      {
                      case 'M' :
                      case 'A' :
                      case 'C' :
                      case 'R' :
                      case 'W' :
                      case 'S' :
                      case 'Y' :
                      case 'B' :
                      case 'D' :
                      case 'H' :
                      case 'V' :
                      case 'X' : { break;}
                      default : return 0;
                      }
                    break;
                  }
                case 'R' :
                  {
                    switch(b)
                      {
                      case 'R' :
                      case 'A' :
                      case 'G' :
                      case 'M' :
                      case 'W' :
                      case 'S' :
                      case 'K' :
                      case 'B' :
                      case 'D' :
                      case 'H' :
                      case 'V' :
                      case 'X' : { break;}
                      default : return 0;
                      }
                    break;
                  }
                  
                case 'W' :
                  {
                    switch(b)
                      {
                      case 'W' :
                      case 'A' :
                      case 'T' :
                      case 'M' :
                      case 'R' :
                      case 'Y' :
                      case 'K' :
                      case 'B' :
                      case 'D' :
                      case 'H' :
                      case 'V' :
                      case 'X' : { break;}
                      default : return 0;
                      }
                    break;
                  }
                  
                case 'S' :
                  {
                    switch(b)
                      {
                      case 'S' :
                      case 'C' :
                      case 'G' :
                      case 'M' :
                      case 'R' :
                      case 'Y' :
                      case 'K' :
                      case 'B' :
                      case 'D' :
                      case 'H' :
                      case 'V' :
                      case 'X' : { break;}
                      default : return 0;
                      }
                    break;
                  }
                  
                case 'Y' :
                  {
                    switch(b)
                      {
                      case 'Y' :
                      case 'C' :
                      case 'T' :
                      case 'M' :
                      case 'W' :
                      case 'S' :
                      case 'K' :
                      case 'B' :
                      case 'D' :
                      case 'H' :
                      case 'V' :
                      case 'X' : { break;}
                      default : return 0;
                      }
                    break;
                  }
                  
                case 'K' :
                  {
                    switch(b)
                      {
                      case 'K' :
                      case 'G' :
                      case 'T' :
                      case 'R' :
                      case 'W' :
                      case 'S' :
                      case 'Y' :
                      case 'B' :
                      case 'D' :
                      case 'H' :
                      case 'V' :
                      case 'X' : { break;}
                      default : return 0;
                      }
                    break;
                  }
                case 'B' :
                  {
                    switch(b)
                      {
                      case 'B' :
                      case 'C' :
                      case 'G' :
                      case 'T' :
                      case 'M' :
                      case 'R' :
                      case 'W' :
                      case 'S' :
                      case 'Y' :
                      case 'K' :
                      case 'D' :
                      case 'H' :
                      case 'V' :
                      case 'X' : { break;}
                      default : return 0;
                      }
                    break;
                  }
                case 'D' :
                  {
                    switch(b)
                      {
                      case 'D' :
                      case 'A' :
                      case 'G' :
                      case 'T' :
                      case 'M' :
                      case 'R' :
                      case 'W' :
                      case 'S' :
                      case 'Y' :
                      case 'K' :
                      case 'B' :
                      case 'H' :
                      case 'V' :
                      case 'X' : { break;}
                      default : return 0;
                      }
                    break;
                  }
                case 'H' :
                  {
                    switch(b)
                      {
                      case 'H' :
                      case 'A' :
                      case 'C' :
                      case 'T' :
                      case 'M' :
                      case 'R' :
                      case 'W' :
                      case 'S' :
                      case 'Y' :
                      case 'K' :
                      case 'B' :
                      case 'D' :
                      case 'V' :
                      case 'X' : { break;}
                      default : return 0;
                      }
                    break;
                  }
                case 'V' :
                  {
                    switch(b)
                      {
                      case 'V' :
                      case 'A' :
                      case 'C' :
                      case 'G' :
                      case 'M' :
                      case 'R' :
                      case 'W' :
                      case 'S' :
                      case 'Y' :
                      case 'K' :
                      case 'B' :
                      case 'D' :
                      case 'H' :
                      case 'X' : { break;}
                      default : return 0;
                      }
                    break;
                  }
                case 'X' :
                  {
                    switch(b)
                      {
                      case 'X' :
                      case 'A' :
                      case 'C' :
                      case 'G' :
                      case 'T' :
                      case 'M' :
                      case 'R' :
                      case 'W' :
                      case 'S' :
                      case 'Y' :
                      case 'K' :
                      case 'B' :
                      case 'D' :
                      case 'H' :
                      case 'V' : { break;}
                      default : return 0;
                      }
                    break;
                  }
                default :
                  {
                    PhyML_Printf("\n. Err. in Are_Compatible.");
                    PhyML_Printf("\n. Please check that characters `%c` and `%c`",a,b);
                    PhyML_Printf("\n. correspond to existing nucleotides.\n");
                    Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
                    return 0;
                  }
                }
            }
        }
    }
  else if(datatype == AA)
    {
      a = statea[0]; b = stateb[0];
      switch(a)
        {
        case 'A' :
          {
            switch(b)
              {
              case 'A' :
              case 'X' : { break;}
              default : return 0;
              }
            break;
          }
        case 'R' :
          {
            switch(b)
              {
              case 'R' :
              case 'X' : { break;}
              default : return 0;
              }
            break;
          }
        case 'N' :
          {
            switch(b)
              {
              case 'N' :
              case 'B' :
              case 'X' : { break;}
              default : return 0;
              }
            break;
          }
        case 'B' :
          {
            switch(b)
              {
              case 'N' :
              case 'B' :
              case 'X' : { break;}
              default : return 0;
              }
            break;
          }
        case 'D' :
          {
            switch(b)
              {
              case 'D' :
              case 'X' : { break;}
              default : return 0;
              }
            break;
          }
        case 'C' :
          {
            switch(b)
              {
              case 'C' :
              case 'X' : { break;}
              default : return 0;
              }
            break;
          }
        case 'Q' :
          {
            switch(b)
              {
              case 'Q' :
              case 'Z' :
              case 'X' : { break;}
              default : return 0;
              }
            break;
          }
        case 'Z' :
          {
            switch(b)
              {
              case 'Q' :
              case 'Z' :
              case 'X' : { break;}
              default : return 0;
              }
            break;
          }
        case 'E' :
          {
            switch(b)
              {
              case 'E' :
              case 'X' : { break;}
              default : return 0;
              }
            break;
          }
        case 'G' :
          {
            switch(b)
              {
              case 'G' :
              case 'X' : { break;}
              default : return 0;
              }
            break;
          }
        case 'H' :
          {
            switch(b)
              {
              case 'H' :
              case 'X' : { break;}
              default : return 0;
              }
            break;
          }
        case 'I' :
          {
            switch(b)
              {
              case 'I' :
              case 'X' : { break;}
              default : return 0;
              }
            break;
          }
        case 'L' :
          {
            switch(b)
              {
              case 'L' :
              case 'X' : { break;}
              default : return 0;
              }
            break;
          }
        case 'K' :
          {
            switch(b)
              {
              case 'K' :
              case 'X' : { break;}
              default : return 0;
              }
            break;
          }
        case 'M' :
          {
            switch(b)
              {
              case 'M' :
              case 'X' : { break;}
              default : return 0;
              }
            break;
          }
        case 'F' :
          {
            switch(b)
              {
              case 'F' :
              case 'X' : { break;}
              default : return 0;
              }
            break;
          }
        case 'P' :
          {
            switch(b)
              {
              case 'P' :
              case 'X' : { break;}
              default : return 0;
              }
            break;
          }
        case 'S' :
          {
            switch(b)
              {
              case 'S' :
              case 'X' : { break;}
              default : return 0;
              }
            break;
          }
        case 'T' :
          {
            switch(b)
              {
              case 'T' :
              case 'X' : { break;}
              default : return 0;
              }
            break;
          }
        case 'W' :
          {
            switch(b)
              {
              case 'W' :
              case 'X' : { break;}
              default : return 0;
              }
            break;
          }
        case 'Y' :
          {
            switch(b)
              {
              case 'Y' :
              case 'X' : { break;}
              default : return 0;
              }
            break;
          }
        case 'V' :
          {
            switch(b)
              {
              case 'V' :
              case 'X' : { break;}
              default : return 0;
              }
            break;
          }
        case 'X' :
          {
            switch(b)
              {
              case 'A':case 'R':case 'N' :case 'B' :case 'D' :
              case 'C':case 'Q':case 'Z' :case 'E' :case 'G' :
              case 'H':case 'I':case 'L' :case 'K' :case 'M' :
              case 'F':case 'P':case 'S' :case 'T' :case 'W' :
              case 'Y':case 'V': case 'X' : { break;}
              default : return 0;
              }
            break;
          }
        default :
          {
            PhyML_Printf("\n. Err. in Are_Compatible.");
            PhyML_Printf("\n. Please check that characters `%c` and `%c`",a,b);
            PhyML_Printf("\n. correspond to existing amino-acids.\n");
            Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
            return 0;
          }
        }
    }
  else if(datatype == GENERIC)
    {
      if(Is_Ambigu(statea,GENERIC,stepsize) || Is_Ambigu(stateb,GENERIC,stepsize)) return 1;
      else
        {
          int a,b;
          char format[6];
          
          sprintf(format,"%%%dd",stepsize);
          
          if(!sscanf(statea,format,&a))
            {
              PhyML_Printf("\n. statea = %s",statea);
              PhyML_Printf("\n. Err. in file %s at line %d\n",__FILE__,__LINE__);
              Warn_And_Exit("\n. PhyML finished prematurely.");
            }
          if(!sscanf(stateb,format,&b))
            {
              PhyML_Printf("\n. statea = %s",stateb);
              PhyML_Printf("\n. Err. in file %s at line %d\n",__FILE__,__LINE__);
              Warn_And_Exit("\n. PhyML finished prematurely.");
            }
          
          /* 	  PhyML_Printf("\n. %s %d a=%d b=%d ",__FILE__,__LINE__,a,b);  */
          
          if(a == b) return 1;
        }
      return 0;
    }
  
  return 1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Hide_Ambiguities(calign *data)
{
  int i;
  for(i=0;i<data->crunch_len;i++) if(data->ambigu[i]) data->wght[i] = 0;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Copy_Tree(t_tree *ori, t_tree *cpy)
{
  int i,j;

  for(i=0;i<2*ori->n_otu-2;++i)
    {
      for(j=0;j<3;++j)
        {
          if(ori->a_nodes[i]->v[j])
            {
              cpy->a_nodes[i]->v[j] = cpy->a_nodes[ori->a_nodes[i]->v[j]->num];
              cpy->a_nodes[i]->l[j] = ori->a_nodes[i]->l[j];
              cpy->a_nodes[i]->b[j] = cpy->a_edges[ori->a_nodes[i]->b[j]->num];
            }
          else
            {
              cpy->a_nodes[i]->v[j] = NULL;
              cpy->a_nodes[i]->b[j] = NULL;
            }
        }
    }

  for(i=0;i<2*ori->n_otu-3;++i)
    {
      cpy->a_edges[i]->l->v             = ori->a_edges[i]->l->v;
      cpy->a_edges[i]->l_old->v         = ori->a_edges[i]->l_old->v;
      cpy->a_edges[i]->l_var->v         = ori->a_edges[i]->l_var->v;
      cpy->a_edges[i]->l_var_old->v     = ori->a_edges[i]->l_var_old->v;
      cpy->a_edges[i]->left             = ori->a_edges[i]->left ? cpy->a_nodes[ori->a_edges[i]->left->num] : NULL;
      cpy->a_edges[i]->rght             = ori->a_edges[i]->rght ? cpy->a_nodes[ori->a_edges[i]->rght->num] : NULL;
      cpy->a_edges[i]->l_v1             = ori->a_edges[i]->l_v1;
      cpy->a_edges[i]->l_v2             = ori->a_edges[i]->l_v2;
      cpy->a_edges[i]->r_v1             = ori->a_edges[i]->r_v1;
      cpy->a_edges[i]->r_v2             = ori->a_edges[i]->r_v2;
      cpy->a_edges[i]->l_r              = ori->a_edges[i]->l_r;
      cpy->a_edges[i]->r_l              = ori->a_edges[i]->r_l;
      cpy->a_edges[i]->does_exist       = ori->a_edges[i]->does_exist;

#ifdef BEAGLE
      cpy->a_edges[i]->p_lk_left_idx    = ori->a_edges[i]->p_lk_left_idx;
      cpy->a_edges[i]->p_lk_rght_idx    = ori->a_edges[i]->p_lk_rght_idx;
      cpy->a_edges[i]->p_lk_tip_idx     = ori->a_edges[i]->p_lk_tip_idx;
#endif
    }


  for(i=0;i<ori->n_otu;++i)
    {
      cpy->a_nodes[i]->tax = YES;

      Free(cpy->a_nodes[i]->name);

      cpy->a_nodes[i]->name = (char *)mCalloc(strlen(ori->a_nodes[i]->name)+1,sizeof(char));
      cpy->a_nodes[i]->ori_name = cpy->a_nodes[i]->name ;

      strcpy(cpy->a_nodes[i]->name,ori->a_nodes[i]->name);
    }

  if(ori->n_root)
    {
      cpy->e_root = cpy->a_edges[ori->e_root->num];
      cpy->n_root = cpy->a_nodes[ori->n_root->num];
      Add_Root(cpy->e_root,cpy);
    }
  else
    {
      cpy->a_edges[2*cpy->n_otu-3]->l->v           = ori->a_edges[2*cpy->n_otu-3]->l->v;
      cpy->a_edges[2*cpy->n_otu-3]->l_old->v       = ori->a_edges[2*cpy->n_otu-3]->l_old->v;
      cpy->a_edges[2*cpy->n_otu-3]->l_var->v       = ori->a_edges[2*cpy->n_otu-3]->l_var->v;
      cpy->a_edges[2*cpy->n_otu-3]->l_var_old->v   = ori->a_edges[2*cpy->n_otu-3]->l_var_old->v;
      cpy->a_edges[2*cpy->n_otu-2]->l->v           = ori->a_edges[2*cpy->n_otu-2]->l->v;
      cpy->a_edges[2*cpy->n_otu-2]->l_old->v       = ori->a_edges[2*cpy->n_otu-2]->l_old->v;
      cpy->a_edges[2*cpy->n_otu-2]->l_var->v       = ori->a_edges[2*cpy->n_otu-2]->l_var->v;
      cpy->a_edges[2*cpy->n_otu-2]->l_var_old->v   = ori->a_edges[2*cpy->n_otu-2]->l_var_old->v;
    }

  cpy->num_curr_branch_available = 0;
  cpy->t_beg                     = ori->t_beg;
  cpy->verbose                   = ori->verbose;


#ifdef BEAGLE
  cpy->b_inst = ori->b_inst;
#endif

  if(cpy->next) Copy_Tree(cpy,cpy->next);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Prune_Subtree(t_node *a, t_node *d, t_edge **target, t_edge **residual, t_tree *tree)
{
  t_node *v1, *v2, *buff_nd;
  t_edge *b1, *b2;
  int dir_v1, dir_v2;
  int i;
  phydbl *buff_p_lk;
  int *buff_scale;
  int *buff_p_pars;
  int *buff_pars;
  int *buff_p_lk_loc, *buff_patt_id;
  int *buff_ui;
  phydbl *buff_p_lk_tip;

  assert(a);
  assert(d);
  assert(tree);

  
  if(tree->n_root && a == tree->n_root) 
    {
      if(d == tree->e_root->left) a = tree->e_root->rght;
      else if(d == tree->e_root->rght) a = tree->e_root->left;
      else
        {
          PhyML_Printf("\n. left: %d right: %d",tree->e_root->left->num,tree->e_root->rght->num);
          assert(false);
        }
    }

  if(a->tax) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);

  dir_v1 = dir_v2 = -1;
  for(i=0;i<3;i++)
    {
      if(a->v[i] != d)
        {
          if(dir_v1 < 0) dir_v1 = i;
          else           dir_v2 = i;
        }
    }

  assert(dir_v1 > -1);
  assert(dir_v2 > -1);

  assert(a->v[dir_v1] != NULL);
  assert(a->v[dir_v2] != NULL);
  
  if(a->v[dir_v1] == a->anc) a->v[dir_v2]->anc = a->v[dir_v1];
  else                       a->v[dir_v1]->anc = a->v[dir_v2];

  if(a->v[dir_v1]->num < a->v[dir_v2]->num)
    {
      v1 = a->v[dir_v1];
      v2 = a->v[dir_v2];
      b1 = a->b[dir_v1];
      b2 = a->b[dir_v2];
    }
  else
    {
      v1 = a->v[dir_v2];
      v2 = a->v[dir_v1];
      b1 = a->b[dir_v2];
      b2 = a->b[dir_v1];
    }


  assert(NULL != b1 && NULL != b2);

  if(target)   (*target)   = b1;
  if(residual) (*residual) = b2;


  a->v[dir_v1] = NULL;
  a->v[dir_v2] = NULL;
  a->b[dir_v1] = NULL;
  a->b[dir_v2] = NULL;

#ifdef BEAGLE
  int temp;
#endif

  if(v1 == b1->left)
    {
      b1->rght = v2;

      if(v2 == b2->left)
        {
          if(tree->is_mixt_tree == NO)
            {
              buff_p_lk            = b1->p_lk_rght;
              b1->p_lk_rght        = b2->p_lk_left;
              b2->p_lk_left        = buff_p_lk;
              
              buff_p_lk_tip        = b1->p_lk_tip_r;
              b1->p_lk_tip_r       = b2->p_lk_tip_l;
              b2->p_lk_tip_l       = buff_p_lk_tip;
              
#ifdef BEAGLE
              temp                 = b1->p_lk_rght_idx;
              b1->p_lk_rght_idx    = b2->p_lk_left_idx;
              b2->p_lk_left_idx    = temp;
#endif
              buff_scale           = b1->sum_scale_rght;
              b1->sum_scale_rght   = b2->sum_scale_left;
              b2->sum_scale_left   = buff_scale;
              
              buff_scale             = b1->sum_scale_rght_cat;
              b1->sum_scale_rght_cat = b2->sum_scale_left_cat;
              b2->sum_scale_left_cat = buff_scale;
              
              buff_pars            = b1->pars_r;
              b1->pars_r           = b2->pars_l;
              b2->pars_l           = buff_pars;
              
              buff_ui              = b1->ui_r;
              b1->ui_r             = b2->ui_l;
              b2->ui_l             = buff_ui;
              
              buff_p_pars          = b1->p_pars_r;
              b1->p_pars_r         = b2->p_pars_l;
              b2->p_pars_l         = buff_p_pars;
              
              buff_p_lk_loc        = b1->p_lk_loc_rght;
              b1->p_lk_loc_rght    = b2->p_lk_loc_left;
              b2->p_lk_loc_left    = buff_p_lk_loc;
              
              buff_patt_id         = b1->patt_id_rght;
              b1->patt_id_rght     = b2->patt_id_left;
              b2->patt_id_left     = buff_patt_id;
            }
        }
      else
        {
          if(tree->is_mixt_tree == NO)
            {
              buff_p_lk            = b1->p_lk_rght; /* b1->p_lk_rght = NULL if b1->rght->tax */
              b1->p_lk_rght        = b2->p_lk_rght; /* b2->p_lk_rght = NULL if b2->rght->tax */
              b2->p_lk_rght        = buff_p_lk;
              
              buff_p_lk_tip        = b1->p_lk_tip_r;
              b1->p_lk_tip_r       = b2->p_lk_tip_r;
              b2->p_lk_tip_r       = buff_p_lk_tip;
#ifdef BEAGLE
              temp                 = b1->p_lk_rght_idx;
              b1->p_lk_rght_idx    = b2->p_lk_rght_idx;
              b2->p_lk_rght_idx    = temp;
              
              b2->p_lk_tip_idx     = b1->p_lk_tip_idx;
#endif
              buff_scale           = b1->sum_scale_rght;
              b1->sum_scale_rght   = b2->sum_scale_rght;
              b2->sum_scale_rght   = buff_scale;
              
              buff_pars            = b1->pars_r;
              b1->pars_r           = b2->pars_r;
              b2->pars_r           = buff_pars;
              
              buff_ui              = b1->ui_r;
              b1->ui_r             = b2->ui_r;
              b2->ui_r             = buff_ui;
              
              buff_p_pars          = b1->p_pars_r;
              b1->p_pars_r         = b2->p_pars_r;
              b2->p_pars_r         = buff_p_pars;
              
              buff_p_lk_loc        = b1->p_lk_loc_rght;
              b1->p_lk_loc_rght    = b2->p_lk_loc_rght;
              b2->p_lk_loc_rght    = buff_p_lk_loc;
              
              buff_patt_id         = b1->patt_id_rght;
              b1->patt_id_rght     = b2->patt_id_rght;
              b2->patt_id_rght     = buff_patt_id;
            }
        }
    }
  else
    {
      b1->left = v2;

      if(v2 == b2->left)
        {
          if(tree->is_mixt_tree == NO)
            {
              buff_p_lk            = b1->p_lk_left;
              b1->p_lk_left        = b2->p_lk_left;
              b2->p_lk_left        = buff_p_lk;
              
              buff_p_lk_tip        = b1->p_lk_tip_l;
              b1->p_lk_tip_l       = b2->p_lk_tip_l;
              b2->p_lk_tip_l       = buff_p_lk_tip;
#ifdef BEAGLE
              temp                 = b1->p_lk_left_idx;
              b1->p_lk_left_idx    = b2->p_lk_left_idx;
              b2->p_lk_left_idx    = temp;
#endif
              buff_scale           = b1->sum_scale_left;
              b1->sum_scale_left   = b2->sum_scale_left;
              b2->sum_scale_left   = buff_scale;
              
              buff_scale             = b1->sum_scale_left_cat;
              b1->sum_scale_left_cat = b2->sum_scale_left_cat;
              b2->sum_scale_left_cat = buff_scale;
              
              buff_pars            = b1->pars_l;
              b1->pars_l           = b2->pars_l;
              b2->pars_l           = buff_pars;
              
              buff_ui              = b1->ui_l;
              b1->ui_l             = b2->ui_l;
              b2->ui_l             = buff_ui;
              
              buff_p_pars          = b1->p_pars_l;
              b1->p_pars_l         = b2->p_pars_l;
              b2->p_pars_l         = buff_p_pars;
              
              buff_p_lk_loc        = b1->p_lk_loc_left;
              b1->p_lk_loc_left    = b2->p_lk_loc_left;
              b2->p_lk_loc_left    = buff_p_lk_loc;
              
              buff_patt_id         = b1->patt_id_left;
              b1->patt_id_left     = b2->patt_id_left;
              b2->patt_id_left     = buff_patt_id;
            }
        }
      else
        {
          if(tree->is_mixt_tree == NO)
            {
              buff_p_lk            = b1->p_lk_left;
              b1->p_lk_left        = b2->p_lk_rght; /* b2->p_lk_rght = NULL if b2->rght->tax */
              b2->p_lk_rght        = buff_p_lk;
              
              buff_p_lk_tip        = b1->p_lk_tip_l;
              b1->p_lk_tip_l       = b2->p_lk_tip_r;
              b2->p_lk_tip_r       = buff_p_lk_tip;
#ifdef BEAGLE
              temp                 = b1->p_lk_left_idx;
              b1->p_lk_left_idx    = b2->p_lk_rght_idx;
              b2->p_lk_rght_idx    = temp;
              
              b2->p_lk_tip_idx     = b1->p_lk_tip_idx;
#endif
              buff_scale           = b1->sum_scale_left;
              b1->sum_scale_left   = b2->sum_scale_rght;
              b2->sum_scale_rght   = buff_scale;
              
              buff_scale             = b1->sum_scale_left_cat;
              b1->sum_scale_left_cat = b2->sum_scale_rght_cat;
              b2->sum_scale_rght_cat = buff_scale;
              
              buff_pars            = b1->pars_l;
              b1->pars_l           = b2->pars_r;
              b2->pars_r           = buff_pars;
              
              buff_ui              = b1->ui_l;
              b1->ui_l             = b2->ui_r;
              b2->ui_r             = buff_ui;
              
              buff_p_pars          = b1->p_pars_l;
              b1->p_pars_l         = b2->p_pars_r;
              b2->p_pars_r         = buff_p_pars;
              
              buff_p_lk_loc        = b1->p_lk_loc_left;
              b1->p_lk_loc_left    = b2->p_lk_loc_rght;
              b2->p_lk_loc_rght    = buff_p_lk_loc;
              
              buff_patt_id         = b1->patt_id_left;
              b1->patt_id_left     = b2->patt_id_rght;
              b2->patt_id_rght     = buff_patt_id;
            }
        }
    }

  for(i=0;i<3;i++)
    if(v2->v[i] == a)
      {
        v2->v[i] = v1;
        v2->b[i] = b1;
        break;
      }

#ifdef DEBUG
  if(i == 3)
    {
      PhyML_Printf("\n. Err. in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("\n. PhyML finished prematurely.");
    }
#endif

  for(i=0;i<3;i++)
    if(v1->v[i] == a)
      {
        v1->v[i] = v2;
        break;
      }

#ifdef DEBUG
  if(i == 3)
    {
      PhyML_Printf("\n. Err. in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("\n. PhyML finished prematurely.");
    }
#endif

  if(b1->l->onoff == ON)
    {
      b1->l->v     = (b1->l->v + b2->l->v);
      b1->l_var->v = (b1->l_var->v + b2->l_var->v);
    }


  (v1 == b1->left)?
    (Set_Edge_Dirs(b1,v1,v2,tree)):
    (Set_Edge_Dirs(b1,v2,v1,tree));


  if(tree->n_root != NULL)
    {
      // Pruning one of the subtree below n_root->v[2], v2 below a
      if(tree->n_root->v[1] == v1 && tree->n_root->v[2] == a) tree->n_root->v[2] = v2;

      // Pruning one of the subtree below n_root->v[1], v2 below a
      else if(tree->n_root->v[2] == v1 && tree->n_root->v[1] == a) tree->n_root->v[1] = v2;
  
      // Pruning one of the subtree below n_root->v[1], v1 below a
      else if((tree->n_root->v[1] == v2 && tree->n_root->v[2] == a) || 
              (tree->n_root->v[2] == v2 && tree->n_root->v[1] == a))
        {
          tree->e_root = b1;
          if(tree->n_root->v[1] == v2) tree->n_root->v[2] = v1;
          if(tree->n_root->v[2] == v2) tree->n_root->v[1] = v1;
        }
      
      // Prune subtree to the left or to the right of the root node
      else if((tree->n_root->v[1] == a && tree->n_root->v[2] == d) ||
              (tree->n_root->v[1] == d && tree->n_root->v[2] == a))
        {
          tree->e_root = b1;
          tree->n_root->v[1] = v2;
          tree->n_root->v[2] = v1;
        }


      if(tree->n_root->v[1] == tree->e_root->rght)
        {
          buff_nd = tree->n_root->v[1];
          tree->n_root->v[1] = tree->n_root->v[2];
          tree->n_root->v[2] = buff_nd;
        }

      Update_Ancestors(tree->n_root,tree->n_root->v[1],tree);
      Update_Ancestors(tree->n_root,tree->n_root->v[2],tree);
      tree->n_root->anc = NULL;
    }

#ifdef DEBUG
  if(b1->left->tax == YES && b1->rght->tax == NO)
    {
      PhyML_Printf("\n. root: %d root->v1: %d root->v2: %d eroot: %d b1: %d b2: %d v1: %d v2: %d",
                   tree->n_root->num,
                   tree->n_root->v[1]->num,
                   tree->n_root->v[2]->num,
                   tree->e_root->num,
                   b1->num,b2->num,v1->num,v2->num);
      PhyML_Printf("\n. b1->left->num = %d",b1->left->num);
      PhyML_Printf("\n. b1->rght->num = %d",b1->rght->num);
      PhyML_Printf("\n. Err. in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("\n. PhyML finished prematurely.");
    }
#endif

  if(tree->is_mixt_tree == YES) MIXT_Prune_Subtree(a,d,target,residual,tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Graft_Subtree(t_edge *target, t_node *link, t_node *link_daughter, t_edge *residual, t_node *target_nd, t_tree *tree)
{
  t_node *v1, *v2;
  int i, dir_v1, dir_v2;
  phydbl *buff_p_lk;
  int *buff_scale;
  int *buff_p_pars, *buff_pars;
  int *buff_p_lk_loc, *buff_patt_id;
  phydbl *buff_p_lk_tip;
  int *buff_ui;
  t_edge *b_up;

  assert(link);
  assert(tree);
  assert(target);
    
  if(link == tree->n_root)
    {
      assert(link_daughter);
      if(link_daughter == tree->n_root->v[1]) link = tree->n_root->v[2];
      else if(link_daughter == tree->n_root->v[2]) link = tree->n_root->v[1];
      else
        {
          PhyML_Printf("\n. link: %d link_daughter: %d",
                       link->num,
                       link_daughter ? link_daughter->num : -1);
          assert(false);
        }
    }

  dir_v1 = dir_v2 = -1;
  b_up = NULL;
  for(i=0;i<3;i++)
    {
      if(link->v[i] == NULL)
        {
          if(dir_v1 < 0) dir_v1 = i;
          else           dir_v2 = i;
        }
      else b_up = link->b[i];
    }
  
  if(dir_v1 < 0 || dir_v2 < 0)
    {
      PhyML_Printf("\n. link: %d was not pruned in a clean manner...\n",link->num);
      assert(FALSE);
    }


  if(target->left == target->rght->anc)
    {
      link->anc = target->left;
      target->rght->anc = link;
    }
  else
    {
      link->anc = target->rght;
      target->left->anc = link;
    }
    

#ifdef BEAGLE
  int temp;
#endif

  if(target->left->num < target->rght->num)
    {
      v1 = target->left;
      v2 = target->rght;

      if(tree->is_mixt_tree == NO)
        {
          buff_p_lk                    = residual->p_lk_rght;
          residual->p_lk_rght          = target->p_lk_rght;
          target->p_lk_rght            = buff_p_lk;
          
          buff_p_lk_tip                = residual->p_lk_tip_r;
          residual->p_lk_tip_r         = target->p_lk_tip_r;
          target->p_lk_tip_r           = buff_p_lk_tip;
          
#ifdef BEAGLE
          temp                         = residual->p_lk_rght_idx;
          residual->p_lk_rght_idx      = target->p_lk_rght_idx;
          target->p_lk_rght_idx        = temp;
          
          temp                         = residual->p_lk_tip_idx;
          residual->p_lk_tip_idx       = target->p_lk_tip_idx;
          target->p_lk_tip_idx         = temp;
#endif
          
          buff_scale                   = residual->sum_scale_rght;
          residual->sum_scale_rght     = target->sum_scale_rght;
          target->sum_scale_rght       = buff_scale;
          
          buff_scale                   = residual->sum_scale_rght_cat;
          residual->sum_scale_rght_cat = target->sum_scale_rght_cat;
          target->sum_scale_rght_cat   = buff_scale;

          buff_pars                    = residual->pars_r;
          residual->pars_r             = target->pars_r;
          target->pars_r               = buff_pars;
          
          buff_ui                      = residual->ui_r;
          residual->ui_r               = target->ui_r;
          target->ui_r                 = buff_ui;
          
          buff_p_pars                  = residual->p_pars_r;
          residual->p_pars_r           = target->p_pars_r;
          target->p_pars_r             = buff_p_pars;
          
          buff_p_lk_loc                = residual->p_lk_loc_rght;
          residual->p_lk_loc_rght      = target->p_lk_loc_rght;
          target->p_lk_loc_rght        = buff_p_lk_loc;
          
          buff_patt_id                 = residual->patt_id_rght;
          residual->patt_id_rght       = target->patt_id_rght;
          target->patt_id_rght         = buff_patt_id;
        }
    }
  else
    {
      v1 = target->rght;
      v2 = target->left;

      if(tree->is_mixt_tree == NO)
        {
          buff_p_lk                    = residual->p_lk_rght;
          residual->p_lk_rght          = target->p_lk_left;
          target->p_lk_left            = buff_p_lk;
          
          buff_p_lk_tip                = residual->p_lk_tip_r;
          residual->p_lk_tip_r         = target->p_lk_tip_l;
          target->p_lk_tip_l           = buff_p_lk_tip;
          
#ifdef BEAGLE
          temp                         = residual->p_lk_rght_idx;
          residual->p_lk_rght_idx      = target->p_lk_left_idx;
          target->p_lk_left_idx        = temp;
#endif
          
          buff_scale                   = residual->sum_scale_rght;
          residual->sum_scale_rght     = target->sum_scale_left;
          target->sum_scale_left       = buff_scale;
          
          buff_scale                   = residual->sum_scale_rght_cat;
          residual->sum_scale_rght_cat = target->sum_scale_left_cat;
          target->sum_scale_left_cat   = buff_scale;
          
          buff_pars                    = residual->pars_r;
          residual->pars_r             = target->pars_l;
          target->pars_l               = buff_pars;
          
          buff_ui                      = residual->ui_r;
          residual->ui_r               = target->ui_l;
          target->ui_l                 = buff_ui;
          
          buff_p_pars                  = residual->p_pars_r;
          residual->p_pars_r           = target->p_pars_l;
          target->p_pars_l             = buff_p_pars;
          
          buff_p_lk_loc                = residual->p_lk_loc_rght;
          residual->p_lk_loc_rght      = target->p_lk_loc_left;
          target->p_lk_loc_left        = buff_p_lk_loc;
          
          buff_patt_id                 = residual->patt_id_rght;
          residual->patt_id_rght       = target->patt_id_left;
          target->patt_id_left         = buff_patt_id;
        }
    }
  


  for(i=0;i<3;i++)
    if(v2->b[i] == target)
      {
        v2->v[i] = link;
        v2->b[i] = residual;
        break;
      }

  link->v[dir_v2] = v2;
  link->b[dir_v2] = residual;

  residual->left  = link;
  residual->rght  = v2;

  if(v1 == target->left) target->rght = link;
  else                   target->left = link;

  link->v[dir_v1] = v1;
  link->b[dir_v1] = target;


  for(i=0;i<3;i++)
    if(v1->v[i] == v2)
      {
        v1->v[i] = link;
        break;
      }

  if(target->l->onoff == ON)
    {
      target->l->v     /= 2.0;
      target->l_var->v /= 2.0;
    }

  if(residual->l->onoff == ON)
    {
      residual->l->v     = target->l->v;
      residual->l_var->v = target->l_var->v;
    }

  Set_Edge_Dirs(target,target->left,target->rght,tree);
  Set_Edge_Dirs(residual,residual->left,residual->rght,tree);
  Set_Edge_Dirs(b_up,b_up->left,b_up->rght,tree);

  /* printf("\n here %p target: %d [%p] eroot: %d [%p] target_nd: %d [%p] root: %d [%p]", */
  /*        tree, */
  /*        target?target->num:-1,target, */
  /*        tree->e_root?tree->e_root->num:-1,tree->e_root, */
  /*        target_nd?target_nd->num:-1,target_nd, */
  /*        tree->n_root?tree->n_root->num:-1,tree->n_root); */
  /* fflush(NULL); */

  if(tree->n_root != NULL)
    {
      if(target == tree->e_root)
        {
          if(target_nd == v1)                tree->e_root = residual;
          else if(target_nd == v2)           tree->e_root = target;
          else if(target_nd == tree->n_root) tree->e_root = b_up;

          tree->n_root->v[1] = tree->e_root->left;
          tree->n_root->v[2] = tree->e_root->rght;
        }
      
      Update_Ancestors(tree->n_root,tree->n_root->v[1],tree);
      Update_Ancestors(tree->n_root,tree->n_root->v[2],tree);
      tree->n_root->anc = NULL;
    }

  if(tree->is_mixt_tree == YES) MIXT_Graft_Subtree(target,link,link_daughter,residual,target_nd,tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Reassign_Node_Nums(t_node *a, t_node *d, int *curr_ext_node, int *curr_int_node, t_tree *tree)
{
  t_node *buff;
  int i;

  if(a->tax)
    {
      buff = tree->a_nodes[*curr_ext_node];
      tree->a_nodes[*curr_ext_node] = a;
      tree->a_nodes[a->num] = buff;
      buff->num = a->num;
      a->num = *curr_ext_node;
      (*curr_ext_node)++;
    }

  if(d->tax)
    {
      buff = tree->a_nodes[*curr_ext_node];
      tree->a_nodes[*curr_ext_node] = d;
      tree->a_nodes[d->num] = buff;
      buff->num = d->num;
      d->num = *curr_ext_node;
      (*curr_ext_node)++;
      return;
    }
  else
    {
      buff = tree->a_nodes[*curr_int_node];
      tree->a_nodes[*curr_int_node] = d;
      tree->a_nodes[d->num] = buff;
      buff->num = d->num;
      d->num = *curr_int_node;
      (*curr_int_node)++;
    }

  for(i=0;i<3;i++)
    {
      if(d->v[i] != a)
        Reassign_Node_Nums(d,d->v[i],curr_ext_node,curr_int_node,tree);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Reassign_Edge_Nums(t_node *a, t_node *d, int *curr_br, t_tree *tree)
{
  t_edge *buff;
  int i,j;

  for(i=0;i<3;i++)
    if(a->v[i] == d)
      {
        buff = tree->a_edges[*curr_br];
        For(j,2*N_MAX_OTU-3) if(tree->a_edges[j] == a->b[i]) break;
        if(j == 2*N_MAX_OTU-3)
          {
            PhyML_Printf("\n. Err. in file %s at line %d (function '%s').\n",__FILE__,__LINE__,__FUNCTION__);
            Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
          }
        tree->a_edges[*curr_br] = a->b[i];
        tree->a_edges[j] = buff;
        a->b[i]->num = *curr_br;
        (*curr_br)++;
        break;
      }
  
  if(d->tax) return;
  else
    {
      for(i=0;i<3;i++)
        if(d->v[i] != a)
          Reassign_Edge_Nums(d,d->v[i],curr_br,tree);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Find_Mutual_Direction(t_node *n1, t_node *n2, short int *dir_n1_to_n2, short int *dir_n2_to_n1)
{
  int scores[3][3];
  int i,j,k,l;
  
  if(n1 == n2) return;


  for(i=0;i<3;i++)
    {
      for(j=0;j<3;j++)
    {
      scores[i][j] = 0;

      For(k,n1->bip_size[i])
        {
          For(l,n2->bip_size[j])
        {
          if(n1->bip_node[i][k] == n2->bip_node[j][l])
            {
              scores[i][j]++;
              break;
            }
        }
        }
    }
    }

  for(i=0;i<3;i++)
    {
      for(j=0;j<3;j++)
    {
      if(!scores[i][j])
        {
          *dir_n1_to_n2 = i;
          *dir_n2_to_n1 = j;
          return;
        }
    }
    }

  PhyML_Printf("\n. n1=%d n2=%d",n1->num,n2->num);
  PhyML_Printf("\n. Err. in file %s at line %d\n",__FILE__,__LINE__);
  Warn_And_Exit("\n. PhyML finished prematurely.");



/*   for(i=0;i<3;i++) */
/*     { */
/*       n_zero_line = 0; */
/*       for(j=0;j<3;j++) */
/* 	{ */
/* 	  if(!scores[i][j]) n_zero_line++; */
/* 	} */
/*       if(n_zero_line != 2) {*dir_n1_to_n2 = i; break;} */
/*     } */


/*   for(i=0;i<3;i++) */
/*     { */
/*       n_zero_col = 0; */
/*       for(j=0;j<3;j++) */
/* 	{ */
/* 	  if(!scores[j][i]) n_zero_col++; */
/* 	} */
/*       if(n_zero_col != 2) {*dir_n2_to_n1 = i; break;} */
/*     } */

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Update_Dir_To_Tips(t_node *a, t_node *d, t_tree *tree)
{
  int i,j,k;
  short int *inout;
  int d_a;
  int dim;

  dim = 2*tree->n_otu-2;

  inout = (short int *)mCalloc(tree->n_otu,sizeof(short int));
  
  for(i=0;i<3;i++)
    {
      if(a->v[i] == d)
        {
          for(j=0;j<tree->n_otu;j++) inout[j] = 1;
          For(k,a->bip_size[i]) inout[a->bip_node[i][k]->num] = 0;
          for(j=0;j<tree->n_otu;j++) if(inout[tree->a_nodes[j]->num]) tree->t_dir[a->num*dim+tree->a_nodes[j]->num] = i;
          break;
        }
    }
  
  
  if(!d->tax)
    {      
      d_a = -1;
      
      for(i=0;i<3;i++)
        {
          if(d->v[i] != a) Update_Dir_To_Tips(d,d->v[i],tree);
          else if(d->v[i] == a) d_a = i;
        }
      
      for(j=0;j<tree->n_otu;j++) inout[j] = 1;
      For(k,d->bip_size[d_a]) inout[d->bip_node[d_a][k]->num] = 0;
      for(j=0;j<tree->n_otu;j++) if(inout[tree->a_nodes[j]->num]) tree->t_dir[d->num*dim+tree->a_nodes[j]->num] = d_a;
    }
  Free(inout);  
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Fill_Dir_Table(t_tree *tree)
{
  int i,j;
  int dim;

  dim = 2*tree->n_otu-2;
  For(i,dim*dim) tree->t_dir[i] = 0;
  Free_Bip(tree);
  Alloc_Bip(tree);
  Get_Bip(tree->a_nodes[0],tree->a_nodes[0]->v[0],tree);
  Update_Dir_To_Tips(tree->a_nodes[0],tree->a_nodes[0]->v[0],tree);

  for(i=tree->n_otu;i<2*tree->n_otu-2;i++)
    for(j=i;j<2*tree->n_otu-2;j++)
      {
        Find_Mutual_Direction(tree->a_nodes[i],tree->a_nodes[j],
                              &(tree->t_dir[i*dim+j]),
                              &(tree->t_dir[j*dim+i]));
      }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


int Get_Subtree_Size(t_node *a, t_node *d)
{
  int size,i;

  if(d->tax) return 1;
  else
    {
      size = 0;
      for(i=0;i<3;i++)
    if(d->v[i] != a)
      size += Get_Subtree_Size(d,d->v[i]);
    }
  return size;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl Fast_Br_Len(t_edge *b, t_tree *tree, int approx)
{
  phydbl init_min_diff_lk_local = tree->mod->s_opt->min_diff_lk_local;

  tree->mod->s_opt->min_diff_lk_local = 1.E-1;      

  if(tree->is_mixt_tree) MIXT_Br_Len_Brent(b,tree);   
  else Br_Len_Brent(b,tree);

  tree->mod->s_opt->min_diff_lk_local = init_min_diff_lk_local;

  return tree->c_lnL;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
/*!
  Calculate the joint probability of states (nt or aa) at the
  two extremities of a given edge given the matrix of transition
  probabilities, the vector of conditional likelihoods on each
  side of the branch and the vector of equilibrium frequencies.
*/
void Joint_Proba_States_Left_Right(phydbl *Pij, phydbl *p_lk_left, phydbl *p_lk_rght,
                   vect_dbl *pi, int scale_left, int scale_rght,
                   phydbl *F, int n, int site, t_tree *tree)
{
  int i,j;
  phydbl sum = 0.0;

  for(i=0;i<n;i++) F[i] = .0;

  for(i=0;i<n;i++)
    {
      for(j=0;j<n;j++)
    {
      F[i*n+j] =
        pi->v[i] *
        Pij[i*n+j] *
        p_lk_left[i] *
        p_lk_rght[j] *
        POW(2.,-(scale_left + scale_rght));

      sum += F[i*n+j];
    }
    }

  For(i,n*n)
    {
      F[i] /= sum;
      if(isnan(F[i]) || isinf(F[i]))
        {
          for(i=0;i<n;i++) for(j=0;j<n;j++)
            PhyML_Printf("\n. %15G %15G %15G %15G %15G",
                         pi->v[i] ,
                         Pij[i*n+j] ,
                         p_lk_left[i] ,
                         p_lk_rght[j] ,
                         POW(2.,-(scale_left + scale_rght)));

          PhyML_Printf("\n. sum = %G",sum);
          Print_Site(tree->data,site,tree->n_otu,"\n",1,stderr);
          PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
          Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
        }
    }

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl Triple_Dist(t_node *a, t_tree *tree)
{
  if(a->tax) return UNLIKELY;
  else
    {
      Update_PMat_At_Given_Edge(a->b[1],tree);
      Update_PMat_At_Given_Edge(a->b[2],tree);

      Update_Partial_Lk(tree,a->b[0],a);
      Fast_Br_Len(a->b[0],tree,YES);

      Update_Partial_Lk(tree,a->b[1],a);
      Fast_Br_Len(a->b[1],tree,YES);

      Update_Partial_Lk(tree,a->b[2],a);
      Fast_Br_Len(a->b[2],tree,YES);
        
      
      Update_Partial_Lk(tree,a->b[1],a);
      Update_Partial_Lk(tree,a->b[0],a);
    }
  
  return tree->c_lnL;

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl Triple_Dist_Approx(t_node *a, t_edge *b, t_tree *tree)
{
  // !!!!!!!! NOT MIXT PROOF 
  if(a->tax) return UNLIKELY;
  else
    {
      int i;

      for(i=0;i<3;i++)
        if(a->b[i] != b)
          Update_PMat_At_Given_Edge(a->b[i],tree);

      Update_Partial_Lk(tree,b,a);
      Fast_Br_Len(b,tree,YES);

      return tree->c_lnL;

      /* t_node *v0,*v1,*v2; */
      /* phydbl d01,d02,d12; */
      /* t_ll *tips0,*tips1,*tips2; */

      /* d01 = d02 = d12 = 0.0; */

      /* v0 = a->v[0]; */
      /* v1 = a->v[1]; */
      /* v2 = a->v[2]; */

      /* tips0 = Get_List_Of_Reachable_Tips(a,v0,tree); */
      /* tips1 = Get_List_Of_Reachable_Tips(a,v1,tree); */
      /* tips2 = Get_List_Of_Reachable_Tips(a,v2,tree); */
      
      /* d01 = Length_Of_Path_Between_List_Of_Tips(tips0,tips1,tree->mat); */
      /* d02 = Length_Of_Path_Between_List_Of_Tips(tips0,tips2,tree->mat); */
      /* d12 = Length_Of_Path_Between_List_Of_Tips(tips1,tips2,tree->mat); */

      /* a->b[0]->l->v = (d01 + d02 - d12)/2.; */
      /* a->b[1]->l->v = (d01 + d12 - d02)/2.; */
      /* a->b[2]->l->v = (d02 + d12 - d01)/2.; */
      
      /* Free_Linked_List(tips0); */
      /* Free_Linked_List(tips1); */
      /* Free_Linked_List(tips2); */

    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Make_Symmetric(phydbl **F, int size)
{
  int i,j;

  for(i=0;i<size;i++)
    {
      for(j=i+1;j<size;j++)
    {
      (*F)[size*i+j] = ((*F)[size*i+j] + (*F)[size*j+i])/2.;
      (*F)[size*j+i] = (*F)[size*i+j];
    }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Round_Down_Freq_Patt(phydbl **F, t_tree *tree)
{
  int i,j;

  for(i=0;i<tree->mod->ns;i++)
    {
      for(j=0;j<tree->mod->ns;j++)
    {
      (*F)[tree->mod->ns*i+j] = RINT((*F)[tree->mod->ns*i+j]);
    }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Get_Sum_Of_Cells(phydbl *F, t_tree *tree)
{
  int i,j;
  phydbl sum = .0;

  for(i=0;i<tree->mod->ns;i++)
    for(j=0;j<tree->mod->ns;j++)
    sum += F[tree->mod->ns*i+j];

  return sum;
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Divide_Cells(phydbl **F, phydbl div, t_tree *tree)
{
  int i,j;

  for(i=0;i<tree->mod->ns;i++)
    for(j=0;j<tree->mod->ns;j++)
    (*F)[tree->mod->ns*i+j] /= div;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Divide_Mat_By_Vect(phydbl **F, phydbl *vect, int size)
{
  int i,j;
  for(i=0;i<size;i++)
    for(j=0;j<size;j++)
    (*F)[size*i+j] = (*F)[size*i+j] / vect[j];
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Multiply_Mat_By_Vect(phydbl **F, phydbl *vect, int size)
{
  int i,j;
  for(i=0;i<size;i++)
    for(j=0;j<size;j++)
    (*F)[size*i+j] = (*F)[size*i+j] * vect[j];
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Found_In_Subtree(t_node *a, t_node *d, t_node *target, int *match, t_tree *tree)
{
  if(d->tax) return;
  else
    {
      int i;
      if(d == target) *match =  1;
      for(i=0;i<3;i++)
    {
      if(d->v[i] != a)
        Found_In_Subtree(d,d->v[i],target,match,tree);
    }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Get_List_Of_Target_Edges(t_node *a, t_node *d, t_edge **list, int *list_size, t_tree *tree)
{
  int i;

  for(i=0;i<3;i++)
    {
      if(a->v[i] && a->v[i] == d)
        {
          list[*list_size] = a->b[i];
          (*list_size)++;
        }
    }
  
  if(d->tax) return;
  else
    {
      for(i=0;i<3;i++)
        {
          if(d->v[i] != a)
            Get_List_Of_Target_Edges(d,d->v[i],list,list_size,tree);
        }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Fix_All(t_tree *tree)
{
  int i;

  for(i=tree->n_otu;i<2*tree->n_otu-2;i++)
    {
      tree->a_nodes[i]->b[0]->l_old->v = tree->a_nodes[i]->b[0]->l->v;
      tree->a_nodes[i]->b[1]->l_old->v = tree->a_nodes[i]->b[1]->l->v;
      tree->a_nodes[i]->b[2]->l_old->v = tree->a_nodes[i]->b[2]->l->v;
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Record_Br_Len(t_tree *mixt_tree)
{
  int i;
  t_tree *tree;

  if(mixt_tree->br_len_recorded == YES)
    {
      PhyML_Printf("\n. Overwriting recorded edge lengths.\n");
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);    
    }

  tree = mixt_tree;

  do
    {
      For(i,2*tree->n_otu-1) tree->a_edges[i]->l_old->v     = tree->a_edges[i]->l->v;
      For(i,2*tree->n_otu-1) tree->a_edges[i]->l_var_old->v = tree->a_edges[i]->l_var->v;
      tree = tree->next;
    }
  while(tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

scalar_dbl **Copy_Br_Len(t_tree *mixt_tree)
{
  int i;
  scalar_dbl **bl, *new_l;
  t_edge *e;

  bl = (scalar_dbl **)mCalloc(2*mixt_tree->n_otu-1,sizeof(scalar_dbl *));
  
  For(i,2*mixt_tree->n_otu-1) 
    {
      e = mixt_tree->a_edges[i];
      bl[i] = (scalar_dbl *)mCalloc(1,sizeof(scalar_dbl));
      do
        {
          bl[i]->v = e->l->v;
          e = e->next;
          if(e) 
            {
              new_l = (scalar_dbl *)mCalloc(1,sizeof(scalar_dbl));
              bl[i]->next       = new_l;
              bl[i]->next->prev = bl[i];
              bl[i]             = bl[i]->next;
            }
        }
      while(e);
    }

  return(bl);
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

scalar_dbl **Copy_Br_Len_Var(t_tree *mixt_tree)
{
  int i;
  scalar_dbl **bl_var, *new_l_var;
  t_edge *e;

  bl_var = (scalar_dbl **)mCalloc(2*mixt_tree->n_otu-1,sizeof(scalar_dbl *));
  
  For(i,2*mixt_tree->n_otu-1) 
    {
      e = mixt_tree->a_edges[i];
      bl_var[i] = (scalar_dbl *)mCalloc(1,sizeof(scalar_dbl));
      do
        {
          bl_var[i]->v = e->l_var->v;
          e = e->next;
          if(e) 
            {
              new_l_var = (scalar_dbl *)mCalloc(1,sizeof(scalar_dbl));
              bl_var[i]->next       = new_l_var;
              bl_var[i]->next->prev = bl_var[i];
              bl_var[i]             = bl_var[i]->next;
            }
        }
      while(e);
    }

  return(bl_var);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Transfer_Br_Len_To_Tree(scalar_dbl **bl, t_tree *tree)
{
  int i;
  scalar_dbl *la, *lb;

  For(i,2*tree->n_otu-1) 
    {
      if(tree->a_edges[i]->l != NULL)
        {
          la = bl[i];
          lb = tree->a_edges[i]->l;
          if(lb != NULL && la != NULL)
            {
              do 
                { 
                  lb->v = la->v;
                  if(la) la = la->next;
                  if(lb) lb = lb->next;
                }
              while(la != NULL && lb != NULL);
              assert(la == NULL && lb == NULL);
            }
        }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Restore_Br_Len(t_tree *mixt_tree)
{
  int i;
  t_tree *tree;

  mixt_tree->br_len_recorded = NO;

  tree = mixt_tree;

  do
    {
      For(i,2*tree->n_otu-1) tree->a_edges[i]->l->v     = tree->a_edges[i]->l_old->v;
      For(i,2*tree->n_otu-1) tree->a_edges[i]->l_var->v = tree->a_edges[i]->l_var_old->v;
      tree = tree->next;
    }
  while(tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Get_Dist_Btw_Edges(t_node *a, t_node *d, t_tree *tree)
{
  int i;
  t_edge *b_fcus;

  b_fcus = NULL;
  for(i=0;i<3;i++) if(a->v[i] == d) {b_fcus = a->b[i]; break;}

  if(d->tax) return;
  else
    {
      for(i=0;i<3;i++)
    if(d->v[i] != a)
      {
        d->b[i]->topo_dist_btw_edges = b_fcus->topo_dist_btw_edges + 1;
        d->b[i]->dist_btw_edges      = b_fcus->dist_btw_edges + d->b[i]->l->v / 2.;
        Get_Dist_Btw_Edges(d,d->v[i],tree);
      }
    }


}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Detect_Polytomies(t_edge *b, phydbl l_thresh, t_tree *tree)
{
  if((b->l->v < l_thresh) && (!b->left->tax) && (!b->rght->tax))
    {
      b->l->v               = 0.0;
      b->has_zero_br_len = YES;
    }
  else b->has_zero_br_len = NO;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Get_List_Of_Nodes_In_Polytomy(t_node *a, t_node *d, t_node ***list, int *size_list)
{
  if(d->tax) return;
  else
    {
      int i;

      for(i=0;i<3;i++)
    {
      if(d->v[i] != a)
        {
          if(!d->b[i]->has_zero_br_len)
        {
          (*list)[*size_list] = d->v[i];
          (*size_list)++;
        }

          if(d->b[i]->has_zero_br_len)
        Get_List_Of_Nodes_In_Polytomy(d,d->v[i],list,size_list);
        }
    }
    }

}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Path_Length(t_node *dep, t_node *arr, phydbl *len, t_tree *tree)
{
  if(dep==arr) return;
  else
    {
      t_edge *next;

      next = dep->b[tree->t_dir[dep->num*(2*tree->n_otu-2)+arr->num]];

      if(next == tree->e_root)
        {
          (*len) += (tree->n_root->l[1] + tree->n_root->l[2]);
        }
      else
        {
          (*len) += next->l->v;
        }
      Path_Length(dep->v[tree->t_dir[dep->num*(2*tree->n_otu-2)+arr->num]],arr,len,tree);
      return;
    }
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Check_Path(t_node *a, t_node *d, t_node *target, t_tree *tree)
{
  PhyML_Printf("path---------\n");
  if(d==target) return;
  else Check_Path(d,d->v[tree->t_dir[d->num*(2*tree->n_otu-2)+target->num]],target,tree);
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Connect_Two_Nodes(t_node *a, t_node *d)
{
  a->v[0] = d;
  d->v[0] = a;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Get_List_Of_Adjacent_Targets(t_node *a, t_node *d, t_node ***node_list, t_edge ***edge_list, int *list_size, int curr_depth, int max_depth)
{
  int i;
  
  if(a->tax) return;
  
  for(i=0;i<3;i++)
    if(a->v[i] == d)
      {
        if(node_list != NULL) (*node_list)[*list_size] = a;
        if(edge_list != NULL) (*edge_list)[*list_size] = a->b[i];
        (*list_size)++;
      }
  if(curr_depth == max_depth) return;
  if(d->tax) return;
  else
    for(i=0;i<3;i++)
      if(d->v[i] != a) Get_List_Of_Adjacent_Targets(d,d->v[i],node_list,edge_list,list_size,curr_depth+1,max_depth);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Sort_List_Of_Adjacent_Targets(t_edge ***list, int list_size)
{
  t_edge *buff_edge;
  int i,j;

  buff_edge = NULL;

  for(i=0;i<list_size-1;i++)
    {
      for(j=i+1;j<list_size;j++)
    if((*list)[j]->topo_dist_btw_edges < (*list)[i]->topo_dist_btw_edges)
      {
        buff_edge = (*list)[j];
        (*list)[j] = (*list)[i];
        (*list)[i] = buff_edge;
      }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

t_node *Common_Nodes_Btw_Two_Edges(t_edge *a, t_edge *b)
{
  if(a->left == b->left)      return b->left;
  else if(a->left == b->rght) return b->rght;
  else if(a->rght == b->left) return b->left;
  else if(a->rght == b->rght) return b->rght;

  PhyML_Printf("\n. First t_edge = %d (%d %d); Second t_edge = %d (%d %d)\n",
     a->num,a->left->num,a->rght->num,
     b->num,b->left->num,b->rght->num);
  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
  Warn_And_Exit("\n. PhyML finished prematurely.");

  return NULL;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int KH_Test(phydbl *site_lk_M1, phydbl *site_lk_M2, t_tree *tree)
{
  phydbl *delta,mean,sd,obs_stat,threshold;
  int i;


  delta = (phydbl *)mCalloc(tree->data->init_len,sizeof(phydbl));

  threshold = .0;
  mean = .0;
  obs_stat = .0;
  for(i=0;i<tree->n_pattern;i++)
    {
      delta[i] = site_lk_M1[i] - site_lk_M2[i];
      mean += ((int)tree->data->wght[i])*delta[i];
    }

  obs_stat = mean;

  mean /= tree->data->init_len;

  for(i=0;i<tree->data->init_len;i++) delta[i] -= mean;

  sd = .0;
  for(i=0;i<tree->data->init_len;i++) sd += POW(delta[i],2);
  sd /= (phydbl)(tree->data->init_len-1.);

/*   threshold = tree->dnorm_thresh*SQRT(sd*tree->data->init_len); */


/*   PhyML_Printf("\nObs stat = %f Threshold = %f\n",obs_stat,threshold); */
  Free(delta);

  if(obs_stat > threshold) return 1;
  else                     return 0;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Random_Tree(t_tree *tree)
{
  int *is_available,*list_of_nodes;
  int i,node_num,step,n_available;
  phydbl min_edge_len;

  assert(tree);

  min_edge_len = 1.E-3;

  is_available  = (int *)mCalloc(2*tree->n_otu-2,sizeof(int));
  list_of_nodes = (int *)mCalloc(tree->n_otu,    sizeof(int));

  for(i=0;i<tree->n_otu;i++) is_available[i]  = 1;
  for(i=0;i<tree->n_otu;i++) list_of_nodes[i] = i;


  step = 0;
  do
    {
/*       node_num = (int)RINT(rand()/(phydbl)(RAND_MAX+1.0)*(tree->n_otu-1-step)); */
      node_num = Rand_Int(0,tree->n_otu-1-step);
      node_num = list_of_nodes[node_num];
      is_available[node_num] = 0;
      for(i=0;i<tree->n_otu;i++) list_of_nodes[i] = -1;
      n_available = 0;
      For(i,2*tree->n_otu-2) if(is_available[i]) {list_of_nodes[n_available++] = i;}

      tree->a_nodes[node_num]->v[0] = tree->a_nodes[tree->n_otu+step];
      tree->a_nodes[tree->n_otu+step]->v[1] = tree->a_nodes[node_num];

/*       node_num = (int)RINT(rand()/(phydbl)(RAND_MAX+1.0)*(tree->n_otu-2-step)); */
      node_num = Rand_Int(0,tree->n_otu-2-step);
      node_num = list_of_nodes[node_num];
      is_available[node_num] = 0;
      for(i=0;i<tree->n_otu;i++) list_of_nodes[i] = -1;
      n_available = 0;
      For(i,2*tree->n_otu-2) if(is_available[i]) {list_of_nodes[n_available++] = i;}

      tree->a_nodes[node_num]->v[0] = tree->a_nodes[tree->n_otu+step];
      tree->a_nodes[tree->n_otu+step]->v[2] = tree->a_nodes[node_num];

      is_available[tree->n_otu+step] = 1;
      for(i=0;i<tree->n_otu;i++) list_of_nodes[i] = -1;
      n_available = 0;
      For(i,2*tree->n_otu-2) if(is_available[i]) list_of_nodes[n_available++] = i;

      step++;
    }while(step < tree->n_otu-2);

  tree->a_nodes[list_of_nodes[0]]->v[0] = tree->a_nodes[list_of_nodes[1]];
  tree->a_nodes[list_of_nodes[1]]->v[0] = tree->a_nodes[list_of_nodes[0]];

  tree->num_curr_branch_available = 0;
  Connect_Edges_To_Nodes_Recur(tree->a_nodes[0],tree->a_nodes[0]->v[0],tree);

  For(i,2*tree->n_otu-3) if(tree->a_edges[i]->l->v < min_edge_len) tree->a_edges[i]->l->v = min_edge_len;

  Free(is_available);
  Free(list_of_nodes);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// Make sure internal edges have likelihood vectors on both
// sides and external edges have one likelihood vector on the
// lefthand side only
// Note: make sure  p_lk_tips vector are re-initialized after
// calling this function 
void Reorganize_Edges_Given_Lk_Struct(t_tree *tree)
{
  int j,i;

  For(i,2*tree->n_otu-3)
    {
      if(tree->a_edges[i]->p_lk_left && tree->a_edges[i]->left->tax == YES)
        {
          For(j,2*tree->n_otu-3)
            {
              if(!tree->a_edges[j]->p_lk_left && tree->a_edges[j]->left->tax == NO)
                {
                  Swap_Partial_Lk(tree->a_edges[i],tree->a_edges[j],LEFT,LEFT,tree);
                  break;
                }
              if(!tree->a_edges[j]->p_lk_rght && tree->a_edges[j]->rght->tax == NO)
                {
                  Swap_Partial_Lk(tree->a_edges[i],tree->a_edges[j],LEFT,RGHT,tree);
                  break;
                }
            }
        }
      
      if(tree->a_edges[i]->p_lk_rght && tree->a_edges[i]->rght->tax == YES)
        {
          For(j,2*tree->n_otu-3)
            {
              if(!tree->a_edges[j]->p_lk_left && tree->a_edges[j]->left->tax == NO)
                {
                  Swap_Partial_Lk(tree->a_edges[i],tree->a_edges[j],RGHT,LEFT,tree);
                  break;
                }
              if(!tree->a_edges[j]->p_lk_rght && tree->a_edges[j]->rght->tax == NO)
                {
                  Swap_Partial_Lk(tree->a_edges[i],tree->a_edges[j],RGHT,RGHT,tree);
                  break;
                }
            }
        }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Swap_Partial_Lk(t_edge *a, t_edge *b, int side_a, int side_b, t_tree *tree)
{
  phydbl *buff_p_lk;
  int *buff_scale;
  int *buff_p_pars;
  int *buff_pars;
  int *buff_p_lk_loc, *buff_patt_id;
  phydbl *buff_p_lk_tip;
  int *buff_ui;
  
  
  if(side_a == LEFT && side_b == LEFT)
    {
      buff_p_lk    = b->p_lk_left;
      b->p_lk_left = a->p_lk_left;
      a->p_lk_left = buff_p_lk;
      
      buff_p_lk_tip = b->p_lk_tip_l;
      b->p_lk_tip_l = a->p_lk_tip_l;
      a->p_lk_tip_l = buff_p_lk_tip;
      
      buff_scale            = b->sum_scale_left_cat;
      b->sum_scale_left_cat = a->sum_scale_left_cat;
      a->sum_scale_left_cat = buff_scale;
      
      buff_scale        = b->sum_scale_left;
      b->sum_scale_left = a->sum_scale_left;
      a->sum_scale_left = buff_scale;
      
      buff_patt_id    = b->patt_id_left;
      b->patt_id_left = a->patt_id_left;
      a->patt_id_left = buff_patt_id;
      
      buff_p_lk_loc    = b->p_lk_loc_left;
      b->p_lk_loc_left = a->p_lk_loc_left;
      a->p_lk_loc_left = buff_p_lk_loc;
      
      buff_pars = b->pars_l;
      b->pars_l = a->pars_l;
      a->pars_l = buff_pars;
      
      buff_p_pars = b->p_pars_l;
      b->p_pars_l = a->p_pars_l;
      a->p_pars_l = buff_p_pars;
      
      buff_ui = b->ui_l;
      b->ui_l = a->ui_l;
      a->ui_l = buff_ui;    
      
#ifdef BEAGLE
      temp             = b->p_lk_left_idx;
      b->p_lk_left_idx = a->p_lk_left_idx;
      a->p_lk_left_idx = temp;
      
      temp            = b->p_lk_tip_idx;
      b->p_lk_tip_idx = a->p_lk_tip_idx;
      a->p_lk_tip_idx = temp;
#endif
    }

  if(side_a == LEFT && side_b == RGHT)
    {
      buff_p_lk    = b->p_lk_rght;
      b->p_lk_rght = a->p_lk_left;
      a->p_lk_left = buff_p_lk;
      
      buff_p_lk_tip = b->p_lk_tip_r;
      b->p_lk_tip_r = a->p_lk_tip_l;
      a->p_lk_tip_l = buff_p_lk_tip;
      
      buff_scale            = b->sum_scale_rght_cat;
      b->sum_scale_rght_cat = a->sum_scale_left_cat;
      a->sum_scale_left_cat = buff_scale;
      
      buff_scale        = b->sum_scale_rght;
      b->sum_scale_rght = a->sum_scale_left;
      a->sum_scale_left = buff_scale;
      
      buff_patt_id    = b->patt_id_rght;
      b->patt_id_rght = a->patt_id_left;
      a->patt_id_left = buff_patt_id;
      
      buff_p_lk_loc    = b->p_lk_loc_rght;
      b->p_lk_loc_rght = a->p_lk_loc_left;
      a->p_lk_loc_left = buff_p_lk_loc;
      
      buff_pars = b->pars_r;
      b->pars_r = a->pars_l;
      a->pars_l = buff_pars;
      
      buff_p_pars = b->p_pars_r;
      b->p_pars_r = a->p_pars_l;
      a->p_pars_l = buff_p_pars;
      
      buff_ui = b->ui_r;
      b->ui_r = a->ui_l;
      a->ui_l = buff_ui;    
      
#ifdef BEAGLE
      temp             = b->p_lk_rght_idx;
      b->p_lk_rght_idx = a->p_lk_left_idx;
      a->p_lk_left_idx = temp;
      
      temp            = b->p_lk_tip_idx;
      b->p_lk_tip_idx = a->p_lk_tip_idx;
      a->p_lk_tip_idx = temp;
#endif
    }

  if(side_a == RGHT && side_b == LEFT)
    {
      buff_p_lk    = b->p_lk_left;
      b->p_lk_left = a->p_lk_rght;
      a->p_lk_rght = buff_p_lk;
      
      buff_p_lk_tip = b->p_lk_tip_l;
      b->p_lk_tip_l = a->p_lk_tip_r;
      a->p_lk_tip_r = buff_p_lk_tip;
      
      buff_scale            = b->sum_scale_left_cat;
      b->sum_scale_left_cat = a->sum_scale_rght_cat;
      a->sum_scale_rght_cat = buff_scale;
      
      buff_scale        = b->sum_scale_left;
      b->sum_scale_left = a->sum_scale_rght;
      a->sum_scale_rght = buff_scale;
      
      buff_patt_id    = b->patt_id_left;
      b->patt_id_left = a->patt_id_rght;
      a->patt_id_rght = buff_patt_id;
      
      buff_p_lk_loc    = b->p_lk_loc_left;
      b->p_lk_loc_left = a->p_lk_loc_rght;
      a->p_lk_loc_rght = buff_p_lk_loc;
      
      buff_pars = b->pars_l;
      b->pars_l = a->pars_r;
      a->pars_r = buff_pars;
      
      buff_p_pars = b->p_pars_l;
      b->p_pars_l = a->p_pars_r;
      a->p_pars_r = buff_p_pars;
      
      buff_ui = b->ui_l;
      b->ui_l = a->ui_r;
      a->ui_r = buff_ui;    
      
#ifdef BEAGLE
      temp             = b->p_lk_left_idx;
      b->p_lk_left_idx = a->p_lk_rght_idx;
      a->p_lk_rght_idx = temp;
      
      temp            = b->p_lk_tip_idx;
      b->p_lk_tip_idx = a->p_lk_tip_idx;
      a->p_lk_tip_idx = temp;
#endif
    }

  if(side_a == RGHT && side_b == RGHT)
    {
      buff_p_lk    = b->p_lk_rght;
      b->p_lk_rght = a->p_lk_rght;
      a->p_lk_rght = buff_p_lk;
      
      buff_p_lk_tip = b->p_lk_tip_r;
      b->p_lk_tip_r = a->p_lk_tip_r;
      a->p_lk_tip_r = buff_p_lk_tip;
      
      buff_scale            = b->sum_scale_rght_cat;
      b->sum_scale_rght_cat = a->sum_scale_rght_cat;
      a->sum_scale_rght_cat = buff_scale;
      
      buff_scale        = b->sum_scale_rght;
      b->sum_scale_rght = a->sum_scale_rght;
      a->sum_scale_rght = buff_scale;
      
      buff_patt_id    = b->patt_id_rght;
      b->patt_id_rght = a->patt_id_rght;
      a->patt_id_rght = buff_patt_id;
      
      buff_p_lk_loc    = b->p_lk_loc_rght;
      b->p_lk_loc_rght = a->p_lk_loc_rght;
      a->p_lk_loc_rght = buff_p_lk_loc;
      
      buff_pars = b->pars_r;
      b->pars_r = a->pars_r;
      a->pars_r = buff_pars;
      
      buff_p_pars = b->p_pars_r;
      b->p_pars_r = a->p_pars_r;
      a->p_pars_r = buff_p_pars;
      
      buff_ui = b->ui_r;
      b->ui_r = a->ui_r;
      a->ui_r = buff_ui;    
      
#ifdef BEAGLE
      temp             = b->p_lk_rght_idx;
      b->p_lk_rght_idx = a->p_lk_rght_idx;
      a->p_lk_rght_idx = temp;
      
      temp            = b->p_lk_tip_idx;
      b->p_lk_tip_idx = a->p_lk_tip_idx;
      a->p_lk_tip_idx = temp;
#endif
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Random_NNI(int n_moves, t_tree *tree)
{
  int i,j;
  t_edge *b;
  t_node *n1,*n2,*n_target;

  n1 = n2 = NULL;
  b = NULL;
  for(i=0;i<n_moves;++i)
    {
      n_target  = tree->a_nodes[tree->n_otu + (int)((phydbl)rand()/RAND_MAX * (2*tree->n_otu-3-tree->n_otu))];
      for(j=0;j<3;++j) if(!n_target->v[j]->tax) {b = n_target->b[j]; break;}

      for(j=0;j<3;++j) if(b->left->v[j] != b->rght) {n1 = b->left->v[j]; break;}
      for(j=0;j<3;++j) if(b->rght->v[j] != b->left) {n2 = b->rght->v[j]; break;}

      Swap(n1,b->left,b->rght,n2,tree);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Fill_Missing_Dist(matrix *mat)
{
  int i,j;
  
  for(i=0;i<mat->n_otu;i++)
    {
      for(j=i+1;j<mat->n_otu;j++)
        {
          if(i != j)
            {
              if(mat->dist[i][j] < .0)
                {
                  Fill_Missing_Dist_XY(i,j,mat);
                  mat->dist[j][i] = mat->dist[i][j];
                }
            }
        }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Fill_Missing_Dist_XY(int x, int y, matrix *mat)
{

  int i,j;
  phydbl *local_mins,**S1S2;
  int cpt;
  int pos_best_estimate;
  phydbl min_crit, curr_crit;

  local_mins = (phydbl *)mCalloc(mat->n_otu*mat->n_otu,sizeof(phydbl ));
  S1S2       = (phydbl **)mCalloc(mat->n_otu*mat->n_otu,sizeof(phydbl *));
  For(i,mat->n_otu*mat->n_otu) S1S2[i] = (phydbl *)mCalloc(2,sizeof(phydbl));

  cpt = 0;
  for(i=0;i<mat->n_otu;i++)
    {
      if((mat->dist[i][x] > .0) && (mat->dist[i][y] > .0))
        {
          for(j=0;j<mat->n_otu;j++)
            {
              if((mat->dist[j][x] > .0) && (mat->dist[j][y] > .0))
                {
                  if((i != j) && (i != x) && (i != y) && (j != x) && (j != y))
                    {
                      S1S2[cpt][0] = MIN(mat->dist[i][x] + mat->dist[j][y] - mat->dist[i][j] , mat->dist[i][y] + mat->dist[j][x] - mat->dist[i][j]);
                      S1S2[cpt][1] = MAX(mat->dist[i][x] + mat->dist[j][y] - mat->dist[i][j] , mat->dist[i][y] + mat->dist[j][x] - mat->dist[i][j]);
                      cpt++;
                    }
                }
            }
        }
    }
  
  Qksort_Matrix(S1S2,0,0,cpt-1);
  
  local_mins[0] = S1S2[0][1];
  for(i=1;i<cpt;i++) local_mins[i] = (i*local_mins[i-1] + S1S2[i][1])/(phydbl)(i+1);
  
  pos_best_estimate = 0;
  min_crit = curr_crit = BIG;
  
  for(i=0;i<cpt-1;i++)
    {
      if((local_mins[i] < S1S2[i+1][0]) && (local_mins[i] > S1S2[i][0]))
        {
          curr_crit = Least_Square_Missing_Dist_XY(x,y,local_mins[i],mat);
          if(curr_crit < min_crit)
            {
              min_crit = curr_crit;
              pos_best_estimate = i;
            }
        }
    }
  
  mat->dist[x][y] = local_mins[pos_best_estimate];
  mat->dist[y][x] = mat->dist[x][y];
  
  For(i,mat->n_otu*mat->n_otu) Free(S1S2[i]);
  Free(S1S2);
  Free(local_mins);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Least_Square_Missing_Dist_XY(int x, int y, phydbl dxy, matrix *mat)
{
  int i,j;
  phydbl fit;
  
  fit = .0;
  for(i=0;i<mat->n_otu;i++)
    {
      if((mat->dist[i][x] > .0) && (mat->dist[i][y] > .0))
        {
          for(j=0;j<mat->n_otu;j++)
            {
              if((mat->dist[j][x] > .0) && (mat->dist[j][y] > .0))
                {
                  if((i != j) && (i != x) && (i != y) && (j != x) && (j != y))
                    {
                      if(dxy < MIN(mat->dist[i][x] + mat->dist[j][y] - mat->dist[i][j] , mat->dist[i][y] + mat->dist[j][x] - mat->dist[i][j]))
                        {
                          fit += POW((mat->dist[i][x] + mat->dist[j][y]) - (mat->dist[i][y] + mat->dist[j][x]),2);
                        }
                      else if((mat->dist[i][x] + mat->dist[j][y]) < (mat->dist[i][y] + mat->dist[j][x]))
                        {
                          fit += POW(dxy - (mat->dist[i][y] + mat->dist[j][x] - mat->dist[i][j]),2);
                        }
                      else
                        {
                          fit += POW(dxy - (mat->dist[i][x] + mat->dist[j][y] - mat->dist[i][j]),2);
                        }
                    }
                }
            }
        }
    }
  return fit;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Check_Memory_Amount(t_tree *tree)
{
  /* Rough estimate of the amount of memory that has to be used */

  long int nbytes;
  int n_otu;
  t_mod *mod;

  mod    = tree->mod;
  n_otu  = tree->io->n_otu;
  nbytes = 0;

  /* Partial Pars */
  nbytes += (2*n_otu-3) * 2 * tree->data->crunch_len * sizeof(int);
  nbytes += (2*n_otu-3) * 2 * tree->data->crunch_len * sizeof(int);
  nbytes += (2*n_otu-3) * 2 * tree->data->crunch_len * mod->ns * sizeof(int);

  /* Pmat */
  nbytes += (2*n_otu-3) * mod->ras->n_catg * mod->ns * mod->ns * sizeof(phydbl);


  /* Partial Lk */
  nbytes += ((2*n_otu-3) * 2 - tree->n_otu) * tree->data->crunch_len * mod->ras->n_catg * mod->ns * sizeof(phydbl);


  /* Scaling factors */
  nbytes += ((2*n_otu-3) * 2 - tree->n_otu) * tree->data->crunch_len * sizeof(int);



  if(((phydbl)nbytes/(1.E+06)) > 256.)
/*   if(((phydbl)nbytes/(1.E+06)) > 0.) */
    {
      PhyML_Printf("\n\n. WARNING: this analysis requires at least %.0f MB of memory space.\n",(phydbl)nbytes/(1.E+06));
#ifndef BATCH

      char answer;
      if((!tree->io->quiet) && (tree->io->mem_question == YES))
        {
          PhyML_Printf("\n. Do you really want to proceed? [Y/n] ");
          if(scanf("%c", &answer))
            {
              if(answer == '\n') answer = 'Y';
              else if(answer == 'n' || answer == 'N') Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
              else getchar();
            }
          else
            {
              Warn_And_Exit("\n\n");
            }
        }
#endif
    }
  else if(((phydbl)nbytes/(1.E+06)) > 100.)
    {
      if(!tree->io->quiet) PhyML_Printf("\n\n. WARNING: this analysis will use at least %.0f MB of memory space...\n",(phydbl)nbytes/(1.E+06));
    }
  else if(((phydbl)nbytes/(1.E+06)) > 1.)
    {
      if(!tree->io->quiet) PhyML_Printf("\n\n. This analysis requires at least %.0f MB of memory space.\n",(phydbl)nbytes/(1.E+06));
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


int Get_State_From_Partial_Lk(phydbl *p_lk, int pos, t_tree *tree)
{
  int i;
  for(i=0;i<tree->mod->ns;i++) if(p_lk[pos+i] > .0) return i;
  return -1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


int Get_State_From_Partial_Pars(short int *p_pars, int pos, t_tree *tree)
{
  int i;
  for(i=0;i<tree->mod->ns;i++) if(p_pars[pos+i] > .0) return i;
  return -1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Check_Dirs(t_tree *tree)
{
  int i;

  For(i,2*tree->n_otu-3)
    {
      if(!tree->a_edges[i]->left->tax)
    {
      if(tree->a_edges[i]->left->v[tree->a_edges[i]->l_v1]->num <
         tree->a_edges[i]->left->v[tree->a_edges[i]->l_v2]->num)
        {
          PhyML_Printf("\n. Edge %d ; v1=%d v2=%d",
             tree->a_edges[i]->num,
             tree->a_edges[i]->left->v[tree->a_edges[i]->l_v1]->num,
             tree->a_edges[i]->left->v[tree->a_edges[i]->l_v2]->num);
          PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
          Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
        }
    }

      if(!tree->a_edges[i]->rght->tax)
    {
      if(tree->a_edges[i]->rght->v[tree->a_edges[i]->r_v1]->num <
         tree->a_edges[i]->rght->v[tree->a_edges[i]->r_v2]->num)
        {
          PhyML_Printf("\n. Edge %d ; v3=%d v4=%d",
             tree->a_edges[i]->num,
             tree->a_edges[i]->rght->v[tree->a_edges[i]->r_v1]->num,
             tree->a_edges[i]->rght->v[tree->a_edges[i]->r_v2]->num);
          PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
          Warn_And_Exit("\n. PhyML finished prematurely.");
        }
    }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Warn_And_Exit(const char *s)
{
  PhyML_Fprintf(stderr,"%s",s);
  fflush(NULL);
#ifndef BATCH
  PhyML_Fprintf(stderr,"\n. Type enter to exit.\n");
  Exit("");
#endif
  Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

// Apply random prune and regraft moves to an existing tree. As opposed to Random_Tree, using this 
// function does not break the likelihood structure.
void Randomize_Tree(t_tree *tree, int n_prune_regraft)
{
  t_node *rnd_node;
  t_edge *rnd_edge,*b_target,*b_residual,**target_list;
  int n_targets,n_rand,i;

  target_list = (t_edge **)mCalloc(2*tree->n_otu-3,sizeof(t_edge *));

  n_rand = n_prune_regraft;
  do
    {
      rnd_node = tree->a_nodes[Rand_Int(tree->n_otu,2*tree->n_otu-3)];
      assert(rnd_node != tree->n_root && rnd_node->tax == NO);
      
      rnd_edge = rnd_node->b[Rand_Int(0,2)];
      
      Prune_Subtree(rnd_node,
                    rnd_node == rnd_edge->left ? rnd_edge->rght : rnd_edge->left,
                    &b_target,
                    &b_residual,
                    tree);
      
      n_targets = 0;
      for(i=0;i<3;i++)
        if(b_target->left->v[i] != b_target->rght)
          Get_List_Of_Adjacent_Targets(b_target->left,b_target->left->v[i],NULL,&target_list,&n_targets,0,tree->n_otu);
      
      for(i=0;i<3;i++)
        if(b_target->rght->v[i] != b_target->left)
          Get_List_Of_Adjacent_Targets(b_target->rght,b_target->rght->v[i],NULL,&target_list,&n_targets,0,tree->n_otu);
      
      if(n_targets > 0) b_target = target_list[Rand_Int(0,n_targets-1)];
      
      assert(b_target != NULL);
      
      Graft_Subtree(b_target,rnd_node,NULL,b_residual,NULL,tree);
      
      n_rand--;
    }
  while(n_rand > 0);

  Free(target_list);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Randomize_Sequence_Order(calign *cdata)
{
  int i,exchange_with;
  phydbl buff_dbl;
  char *buff_name,*buff_state;
  short int *buff_ambigu;

  exchange_with = -1;
  for(i=0;i<cdata->n_otu;i++)
    {
      buff_dbl  = rand();
      buff_dbl /= (RAND_MAX+1.);
      buff_dbl *= cdata->n_otu;
      exchange_with = (int)FLOOR(buff_dbl);

      buff_name                         = cdata->c_seq[i]->name;
      cdata->c_seq[i]->name             = cdata->c_seq[exchange_with]->name;
      cdata->c_seq[exchange_with]->name = buff_name;

      buff_state                         = cdata->c_seq[i]->state;
      cdata->c_seq[i]->state             = cdata->c_seq[exchange_with]->state;
      cdata->c_seq[exchange_with]->state = buff_state;

      buff_ambigu                            = cdata->c_seq[i]->is_ambigu;
      cdata->c_seq[i]->is_ambigu             = cdata->c_seq[exchange_with]->is_ambigu;
      cdata->c_seq[exchange_with]->is_ambigu = buff_ambigu;
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Update_Root_Pos(t_tree *tree)
{
  if(tree->n_root_pos > -1.0)
    {
      tree->n_root->l[2] = tree->e_root->l->v * tree->n_root_pos;
      tree->n_root->l[1] = tree->e_root->l->v * (1.-tree->n_root_pos);
    }
  else
    {
/*       tree->n_root->l[0] = tree->e_root->l->v / 2.; */
/*       tree->n_root->l[1] = tree->e_root->l->v / 2.; */
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Add_Root(t_edge *target, t_tree *tree)
{
  t_edge *b1, *b2;

  assert(target);
  assert(tree);

  #ifndef PHYML
  /* PhyML_Printf("\n. Adding root on t_edge %d left = %d right = %d.",target->num,target->left ? target->left->num : -1, target->rght ? target->rght->num : -1); fflush(NULL) */;
  #endif

  tree->e_root = target;

  /* Create the root t_node if it does not exist yet */
  if(!tree->a_nodes[2*tree->n_otu-2]) tree->n_root = (t_node *)Make_Node_Light(2*tree->n_otu-2);
  else                                tree->n_root = tree->a_nodes[2*tree->n_otu-2];

  tree->a_nodes[2*tree->n_otu-2] = tree->n_root;

  tree->n_root->tax = 0;

  /* Set the position of the root */
  tree->n_root->v[0] = NULL;
  tree->n_root->v[1] = tree->e_root->left;
  tree->n_root->v[2] = tree->e_root->rght;

  /* tree->n_root->b[2] = tree->e_root; */
  /* tree->n_root->b[1] = tree->e_root; */

  if(tree->n_root_pos > -1.0)
    {
      if(tree->n_root_pos < 1.E-6 &&  tree->n_root_pos > -1.E-6)
        {
          printf("\n. WARNING: you put the root at a weird position...");
        }
/*       tree->n_root->l[0] = tree->e_root->l->v * (tree->n_root_pos/(1.+tree->n_root_pos)); */
/*       tree->n_root->l[1] = tree->e_root->l->v - tree->n_root->l[0]; */

      tree->n_root->l[2] = tree->e_root->l->v * tree->n_root_pos;
      tree->n_root->l[1] = tree->e_root->l->v * (1. - tree->n_root_pos);
    }
  else
    {
      tree->n_root->l[2] = tree->e_root->l->v / 2.;
      tree->n_root->l[1] = tree->e_root->l->v / 2.;
      tree->n_root_pos = 0.5;
    }

  b1 = tree->a_edges[2*tree->n_otu-3];
  b2 = tree->a_edges[2*tree->n_otu-2];

  tree->n_root->b[0] = NULL;
  tree->n_root->b[1] = b1;
  tree->n_root->b[2] = b2;

  b1->num  = tree->num_curr_branch_available;
  b2->num  = tree->num_curr_branch_available+1;
  b1->left = tree->n_root;
  b1->rght = tree->n_root->v[1];
  b2->left = tree->n_root;
  b2->rght = tree->n_root->v[2];

  b1->l->v     = tree->n_root->l[1];
  b2->l->v     = tree->n_root->l[2];
  b1->l_old->v = tree->n_root->l[1];
  b2->l_old->v = tree->n_root->l[2];

  b1->l_r = 1;
  b2->l_r = 2;

  b1->r_l = 0;
  b2->r_l = 0;

  b1->l_v1 = 0;
  b1->l_v2 = 2;

  b2->l_v1 = 0;
  b2->l_v2 = 1;

  b1->r_v1 = 1;
  b1->r_v2 = 2;

  b2->r_v1 = 1;
  b2->r_v2 = 2;

  /* WARNING: make sure you have freed the memory for p_lk_rght on b1 and b2 */
  if(tree->is_mixt_tree == NO)
    {
      b1->p_lk_rght = tree->e_root->p_lk_left;
      b2->p_lk_rght = tree->e_root->p_lk_rght;
      
      b1->p_lk_tip_r = tree->e_root->p_lk_tip_l;
      b2->p_lk_tip_r = tree->e_root->p_lk_tip_r;
      
      b1->sum_scale_rght = tree->e_root->sum_scale_left;
      b2->sum_scale_rght = tree->e_root->sum_scale_rght;
      
      b1->sum_scale_rght_cat = tree->e_root->sum_scale_left_cat;
      b2->sum_scale_rght_cat = tree->e_root->sum_scale_rght_cat;
      
      b1->p_lk_loc_rght = tree->e_root->p_lk_loc_left;
      b2->p_lk_loc_rght = tree->e_root->p_lk_loc_rght;
      
      b1->pars_r = tree->e_root->pars_l;
      b2->pars_r = tree->e_root->pars_r;
      
      b1->ui_r = tree->e_root->ui_l;
      b2->ui_r = tree->e_root->ui_r;
      
      b1->p_pars_r = tree->e_root->p_pars_l;
      b2->p_pars_r = tree->e_root->p_pars_r;
      
      b1->p_lk_loc_rght = tree->e_root->p_lk_loc_left;
      b2->p_lk_loc_rght = tree->e_root->p_lk_loc_rght;
      
      b1->patt_id_rght = tree->e_root->patt_id_left;
      b2->patt_id_rght = tree->e_root->patt_id_rght;
    }

  Update_Ancestors(tree->n_root,tree->n_root->v[2],tree);
  Update_Ancestors(tree->n_root,tree->n_root->v[1],tree);
  tree->n_root->anc = NULL;

  if(tree->is_mixt_tree == YES) MIXT_Add_Root(target,tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Update_Ancestors(t_node *a, t_node *d, t_tree *tree)
{
  if(d == NULL)
    {
      PhyML_Printf("\n. d is NULL; a: %d root: %d",a->num,tree->n_root->num);
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);    
    }
  d->anc = a;

  if(a == tree->n_root) a->anc = NULL;

  if(d->tax) return;
  else
    {
      int i;
      for(i=0;i<3;i++)
        if((d->v[i] != a) && (d->b[i] != tree->e_root))
          Update_Ancestors(d,d->v[i],tree);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* Generate a random unrooted tree with 'n_otu' OTUs */
#if (defined PHYTIME || defined INVITEE)
t_tree *Generate_Random_Tree_From_Scratch(int n_otu, int rooted)
{
  t_tree *tree;
  int *connected,*nonconnected,*available_nodes;
  int i,n_connected,n_nonconnected,n1,n2,new_n,n_internal,n_external,n_available;
  t_node *root,*curr_n,**internal_nodes, **external_nodes;
  phydbl *t,*tmp;

  tree = Make_Tree_From_Scratch(n_otu,NULL);

  tree->rates = RATES_Make_Rate_Struct(tree->n_otu);
  RATES_Init_Rate_Struct(tree->rates,tree->io->rates,tree->n_otu);
  
  For(i,2*tree->n_otu-2)
    {
      tree->a_nodes[i]->v[1] = NULL;
      tree->a_nodes[i]->v[2] = NULL;
    }
  
  root = (t_node *)Make_Node_Light(2*tree->n_otu-2);
  
  connected       = (int *)mCalloc(2*tree->n_otu-2,sizeof(int));
  nonconnected    = (int *)mCalloc(2*tree->n_otu-2,sizeof(int));
  available_nodes = (int *)mCalloc(2*tree->n_otu-2,sizeof(int));
  internal_nodes  = (t_node **)mCalloc(tree->n_otu-2,sizeof(t_node *));
  external_nodes  = (t_node **)mCalloc(tree->n_otu,  sizeof(t_node *));
  t               = (phydbl *)mCalloc(tree->n_otu-1,sizeof(phydbl ));
  tmp             = (phydbl *)mCalloc(2*tree->n_otu-2,sizeof(phydbl ));
  
  n_nonconnected = 2*n_otu-2;
  
  For(i,2*tree->n_otu-2) nonconnected[i] = i;
  
  available_nodes[0] = 2*n_otu-2;
  
  /* Node times are generated according to a Birth-death process.
     Formulae are as described by Yang and Rannala (1997) */
  phydbl    phi;
  phydbl    rho; /* sampling intensity */
  phydbl     mu; /* birth rate */
  phydbl lambda; /* death rate */
  phydbl      u; /* random U[0,1] */
  phydbl expval;
  
  /* rho = 1.0 and mu = 0.0 correspond to the Yule process */
  
  lambda = 6.7;
  mu     = 2.5;
  rho    = 9./150.;
  
  expval = exp(MIN(1.E+2,mu-lambda));
  phi = (rho*lambda*(expval-1.) + (mu-lambda)*expval)/(expval-1.); /* Equation 16 */
  
  for(i=0;i<tree->n_otu-1;i++)
    {
      u = rand();
      u /= RAND_MAX;
      
      if(FABS(lambda - mu) > 1.E-4)
        t[i] = (log(phi-u*rho*lambda) - log(phi-u*rho*lambda + u*(lambda-mu)))/(mu-lambda); /* Equation 15 */
      else
        t[i] = u / (1.+lambda*rho*(1-u)); /* Equation 17 */
    }
  
  Qksort(t,NULL,0,tree->n_otu-2); /* Node times ordering in ascending order */
  
  for(i=0;i<tree->n_otu-1;i++) tmp[i] =  t[tree->n_otu-2-i];
  for(i=0;i<tree->n_otu-1;i++) t[i]   = -tmp[i];
  
  
  /* Rescale t_node times such that the time at the root t_node is -100 */
  for(i=1;i<tree->n_otu-1;i++)
    {
      t[i] /= -t[0];
      t[i] *= 1.E+02;
    }
  t[0] = -1.E+02;
  
  
  n_available = 1;
  curr_n = root;
  n_connected = 0;
  do
    {
      n1 = Rand_Int(0,n_nonconnected-1);
      n1 = nonconnected[n1];
      connected[n1] = 1;
      
      n_nonconnected = 0;
      For(i,2*tree->n_otu-2) if(!connected[i]) {nonconnected[n_nonconnected++] = i;}
      
      n2 = Rand_Int(0,n_nonconnected-1);
      n2 = nonconnected[n2];
      connected[n2] = 1;
      
      n_nonconnected = 0;
      For(i,2*tree->n_otu-2) if(!connected[i]) {nonconnected[n_nonconnected++] = i;}
      
      curr_n->v[1] = tree->a_nodes[n1];
      curr_n->v[2] = tree->a_nodes[n2];
      tree->a_nodes[n1]->v[0] = curr_n;
      tree->a_nodes[n2]->v[0] = curr_n;
      
      tree->rates->nd_t[curr_n->num] = t[n_connected/2];
      
      available_nodes[n_available] = tree->a_nodes[n1]->num;
      for(i=0;i<n_available;i++)
        if(available_nodes[i] == curr_n->num)
          {
            available_nodes[i] = tree->a_nodes[n2]->num;
            break;
          }
      n_available++;
      
      new_n = Rand_Int(0,n_available-1);
      curr_n = tree->a_nodes[available_nodes[new_n]];
      
      n_connected+=2;
      
    }while(n_connected < 2*tree->n_otu-2);
  
  For(i,2*tree->n_otu-2) tmp[i] = tree->rates->nd_t[i];
  
  /* Unroot the tree */
  root->v[2]->v[0] = root->v[2];
  root->v[1]->v[0] = root->v[1];
  
  n_internal = n_external = 0;
  For(i,2*tree->n_otu-2)
    {
      if(tree->a_nodes[i]->v[1]) internal_nodes[n_internal++] = tree->a_nodes[i];
      else                       external_nodes[n_external++] = tree->a_nodes[i];
    }
  
  
  n_internal = n_external = 0;
  For(i,2*tree->n_otu-2)
    {
      if(i < tree->n_otu)
        {
          tree->a_nodes[i]      = external_nodes[n_external++];
          tree->a_nodes[i]->tax = 1;
        }
      else
        {
          tree->rates->nd_t[i]  = tmp[internal_nodes[n_internal]->num];
          tree->a_nodes[i]      = internal_nodes[n_internal++];
          tree->a_nodes[i]->tax = 0;
        }
      
      tree->a_nodes[i]->num = i;
    }
  
  for(i=0;i<tree->n_otu;i++) tree->rates->nd_t[i] = 0.0;
  
  for(i=0;i<tree->n_otu;i++)
    {
      if(!tree->a_nodes[i]->name) tree->a_nodes[i]->name = (char *)mCalloc(T_MAX_NAME,sizeof(char));
      strcpy(tree->a_nodes[i]->name,"x");
      sprintf(tree->a_nodes[i]->name+1,"%d",i);
    }
  
  
  tree->num_curr_branch_available = 0;
  Connect_Edges_To_Nodes_Recur(tree->a_nodes[0],tree->a_nodes[0]->v[0],tree);
  
  /* Add root */
  if(rooted)
    {
      For(i,2*tree->n_otu-3)
        {
          if(((tree->a_edges[i]->left == root->v[1]) || (tree->a_edges[i]->rght == root->v[1])) &&
             ((tree->a_edges[i]->left == root->v[2]) || (tree->a_edges[i]->rght == root->v[2])))
            {
              Add_Root(tree->a_edges[i],tree);
              break;
            }
        }

    }
  /* Or not... */
  else
    {
      Free_Node(root);
    }
  
  RATES_Random_Branch_Lengths(tree);
  
  Free(available_nodes);
  Free(connected);
  Free(nonconnected);
  Free(external_nodes);
  Free(internal_nodes);
  Free(t);
  Free(tmp);
  
  return tree;
}
#endif
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Random_Lineage_Rates(t_node *a, t_node *d, t_edge *b, phydbl stick_prob, phydbl *rates, int curr_rate, int n_rates, t_tree *tree)
{
  phydbl uni;
  int new_rate;
  int i;


  if(b)
    {
      uni  = rand();
      uni /= RAND_MAX;

      if(uni > stick_prob) /* Randomly pick a new rate */
    {
      uni  = rand();
      uni /= RAND_MAX;
      uni = (phydbl)(uni * (n_rates-1));
      if(uni-(int)(uni) > 0.5-BIG) new_rate = (int)(uni)+1;
      else new_rate = (int)(uni);
    }
      else
    {
      new_rate = curr_rate;
    }

      for(i=0;i<3;i++)
    if(a->v[i] == d)
      {
        a->b[i]->l->v *= rates[new_rate];
        break;
      }

      for(i=0;i<3;i++)
    if(a->v[i] == d)
      {
        if(!(a->b[i]->n_labels%BLOCK_LABELS)) Make_New_Edge_Label(a->b[i]);
        if(rates[new_rate] > 1.0)      strcpy(a->b[i]->labels[a->b[i]->n_labels],"FAST");
        else if(rates[new_rate] < 1.0) strcpy(a->b[i]->labels[a->b[i]->n_labels],"SLOW");
        else                           strcpy(a->b[i]->labels[a->b[i]->n_labels],"MEDIUM");
        a->b[i]->n_labels++;
        break;
      }
      curr_rate = new_rate;
    }

  if(d->tax) return;
  else
    {
      for(i=0;i<3;i++)
    if(d->v[i] != a)
      Random_Lineage_Rates(d,d->v[i],d->b[i],stick_prob,rates,curr_rate,n_rates,tree);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

t_edge *Find_Edge_With_Label(char *label, t_tree *tree)
{
  int i,j;

  For(i,2*tree->n_otu-3)
    {
      For(j,tree->a_edges[i]->n_labels)
    {
      if(!strcmp(tree->a_edges[i]->labels[j],label)) return tree->a_edges[i];
    }
    }
  return NULL;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Evolve(calign *data, t_mod *mod, int first_site_pos, t_tree *tree)
{
  int root_state, root_rate_class;
  int site,i;
  phydbl *orig_l;
  /* phydbl shape,scale,var,mean; */
  int switch_to_yes;
  
  orig_l = (phydbl *)mCalloc(2*tree->n_otu-3,sizeof(phydbl));
  For(i,2*tree->n_otu-3) orig_l[i] = tree->a_edges[i]->l->v;

  data->n_otu = tree->n_otu;
  data->io    = tree->io;

  if(mod->use_m4mod) tree->write_labels = YES;

  Set_Br_Len_Var(NULL,tree);

  switch_to_yes = NO;
  if(tree->mod->gamma_mgf_bl == YES) switch_to_yes = YES;

  Switch_Eigen(YES,mod);                 

  assert(first_site_pos < data->init_len);
  
  for(site=first_site_pos;site<data->init_len;++site)
    {
      Set_Model_Parameters(mod);

     /* Pick the rate class */
      root_state = root_rate_class = -1;
      root_rate_class = Pick_State(mod->ras->n_catg,mod->ras->gamma_r_proba->v);

 
      /* /\* Get the change probability matrices *\/ */
      /* For(i,2*tree->n_otu-3) */
      /*   { */
      /*     var   = MAX(0.0,tree->a_edges[i]->l_var->v) * POW(tree->mod->ras->gamma_rr->v[root_rate_class],2); */
      /*     mean  = orig_l[i]; */

      /*     shape = mean * mean / var; */
      /*     scale = var / mean; */

      /*     tree->a_edges[i]->l->v = Rgamma(shape,scale); */

      /*   } */

      for(i=0;i<2*tree->n_otu-3;++i) Update_PMat_At_Given_Edge(tree->a_edges[i],tree);
      
      /* Pick the root nucleotide/aa */
      root_state = Pick_State(mod->ns,mod->e_frq->pi->v);
      data->c_seq[0]->state[site] = Reciproc_Assign_State(root_state,tree->io->datatype);

      /* printf("\n. root_state: %d root_rate_class: %d [%f %f %f %f]", */
      /*        root_state, */
      /*        root_rate_class, */
      /*        mod->e_frq->pi->v[0], */
      /*        mod->e_frq->pi->v[1], */
      /*        mod->e_frq->pi->v[2], */
      /*        mod->e_frq->pi->v[3]); */
      /* Generic_Exit(__FILE__,__LINE__,__FUNCTION__); */

      /* tree->a_nodes[0] is considered as the root t_node */
      Evolve_Recur(tree->a_nodes[0],
                   tree->a_nodes[0]->v[0],
                   tree->a_nodes[0]->b[0],
                   root_state,
                   root_rate_class,
                   site,
                   data,
                   mod,
                   tree);

/*       PhyML_Printf("%s\n",Write_Tree(tree,NO)); */

      data->wght[site] = 1;
    }
  data->crunch_len = data->init_len;
  /* Print_CSeq(stdout,NO,data); */
  For(i,2*tree->n_otu-3) tree->a_edges[i]->l->v = orig_l[i];
  Free(orig_l);

  if(switch_to_yes == YES) tree->mod->gamma_mgf_bl = YES;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int Pick_State(int n, phydbl *prob)
{
  int pos;
  phydbl uni;

  do
    {
      pos  = rand();
      pos  = (pos % n);
      uni  = (phydbl)rand();
      uni /= (phydbl)RAND_MAX;
      if(uni < prob[pos]) break;
    }
  while(1);

  return (int)pos;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Evolve_Recur(t_node *a, t_node *d, t_edge *b, int a_state, int r_class, int site_num, calign *gen_data, t_mod *mod, t_tree *tree)
{
  int d_state;
  int dim1,dim2;

  dim1 = tree->mod->ns * tree->mod->ns;
  dim2 = tree->mod->ns;

  d_state = Pick_State(mod->ns,b->Pij_rr+r_class*dim1+a_state*dim2);

  /* PhyML_Printf("\n>> %c (%d) L:%G %G %G %G %G",Reciproc_Assign_State(d_state,mod->io->datatype),d_state, */
  /*              b->l->v, */
  /*              b->Pij_rr[r_class*dim1+a_state*dim2+0], */
  /*              b->Pij_rr[r_class*dim1+a_state*dim2+1], */
  /*              b->Pij_rr[r_class*dim1+a_state*dim2+2], */
  /*              b->Pij_rr[r_class*dim1+a_state*dim2+3]); */

  if(mod->use_m4mod)
    {
      phydbl rrate; /* relative rate of substitutions */

      rrate = mod->m4mod->multipl[(int)d_state/mod->m4mod->n_o];
      if(!(b->n_labels%BLOCK_LABELS)) Make_New_Edge_Label(b);
      if(rrate > 1.0) strcpy(b->labels[b->n_labels],"FASTER");
      else strcpy(b->labels[b->n_labels],"SLOWER");
      b->n_labels++;
    }

  if(d->tax)
    {
      gen_data->c_seq[d->num]->state[site_num] = Reciproc_Assign_State(d_state,tree->io->datatype);
      return;
    }
  else
    {
      int i;
      for(i=0;i<3;i++)
        if(d->v[i] != a)
          Evolve_Recur(d,d->v[i],d->b[i],
                       d_state,r_class,site_num,gen_data,
                       mod,tree);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Site_Diversity(t_tree *tree)
{
  int i,j,k,ns;
  int *div,sum;
  
  ns = tree->mod->ns;
  
  div = (int *)mCalloc(ns,sizeof(int));
  
  Site_Diversity_Post(tree->a_nodes[0],tree->a_nodes[0]->v[0],tree->a_nodes[0]->b[0],tree);
  Site_Diversity_Pre (tree->a_nodes[0],tree->a_nodes[0]->v[0],tree->a_nodes[0]->b[0],tree);
  
  For(i,2*tree->n_otu-3)
    {
      for(j=0;j<ns;j++)
        {
          tree->a_edges[i]->div_post_pred_left[j] = 0;
          tree->a_edges[i]->div_post_pred_rght[j] = 0;
        }
    }
  
  for(i=0;i<tree->n_pattern;i++)
    {
      For(j,2*tree->n_otu-3)
        {
          Binary_Decomposition(tree->a_edges[j]->ui_l[i],div,ns);
          sum = 0;
          for(k=0;k<ns;k++) sum += div[k];
          tree->a_edges[j]->div_post_pred_left[sum-1] += tree->data->wght[i];
          
          Binary_Decomposition(tree->a_edges[j]->ui_r[i],div,ns);
          sum = 0;
          for(k=0;k<ns;k++) sum += div[k];
          tree->a_edges[j]->div_post_pred_rght[sum-1] += tree->data->wght[i];
        }
    }
  
  /*   For(j,2*tree->n_otu-3) */
  /*     { */
  /*       PhyML_Printf("\n. Edge %4d   div_left = %4d %4d %4d %4d -- div_rght = %4d %4d %4d %4d", */
  /* 	     j, */
  /* 	     tree->a_edges[j]->div_post_pred_left[0], */
  /* 	     tree->a_edges[j]->div_post_pred_left[1], */
  /* 	     tree->a_edges[j]->div_post_pred_left[2], */
  /* 	     tree->a_edges[j]->div_post_pred_left[3], */
  /* 	     tree->a_edges[j]->div_post_pred_rght[0], */
  /* 	     tree->a_edges[j]->div_post_pred_rght[1], */
  /* 	     tree->a_edges[j]->div_post_pred_rght[2], */
  /* 	     tree->a_edges[j]->div_post_pred_rght[3]); */
  /*     } */
  
  Free(div);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Site_Diversity_Post(t_node *a, t_node *d, t_edge *b, t_tree *tree)
{
  if(d->tax) return;
  else
    {
      int i;

      for(i=0;i<3;i++)
        if(d->v[i] != a)
          Site_Diversity_Post(d,d->v[i],d->b[i],tree);
      
      Subtree_Union(d,b,tree);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Site_Diversity_Pre(t_node *a, t_node *d, t_edge *b, t_tree *tree)
{
  if(d->tax) return;
  else
    {
      int i;

      for(i=0;i<3;i++)
    if(d->v[i] != a)
      {
        Subtree_Union(d,d->b[i],tree);
        Site_Diversity_Pre(d,d->v[i],d->b[i],tree);
      }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Subtree_Union(t_node *n, t_edge *b_fcus, t_tree *tree)
{
/*
           |
       |<- b_cus
       |
       n
          / \
         /   \
        /     \
*/

  int site;
  int *ui, *ui_v1, *ui_v2;

  ui = ui_v1 = ui_v2 = NULL;

  if(n == b_fcus->left)
    {
      ui = b_fcus->ui_l;

      ui_v1 =
      (n == n->b[b_fcus->l_v1]->left)?
      (n->b[b_fcus->l_v1]->ui_r):
      (n->b[b_fcus->l_v1]->ui_l);

      ui_v2 =
      (n == n->b[b_fcus->l_v2]->left)?
      (n->b[b_fcus->l_v2]->ui_r):
      (n->b[b_fcus->l_v2]->ui_l);
    }
  else
    {
      ui = b_fcus->ui_r;

      ui_v1 =
      (n == n->b[b_fcus->r_v1]->left)?
      (n->b[b_fcus->r_v1]->ui_r):
      (n->b[b_fcus->r_v1]->ui_l);

      ui_v2 =
      (n == n->b[b_fcus->r_v2]->left)?
      (n->b[b_fcus->r_v2]->ui_r):
      (n->b[b_fcus->r_v2]->ui_l);
    }

  for(site=0;site<tree->n_pattern;site++) ui[site] = ui_v1[site] | ui_v2[site];

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Binary_Decomposition(int value, int *bit_vect, int size)
{
  int i,cumul;

  for(i=0;i<size;i++) bit_vect[i] = 0;

  cumul = 0;
  for(i=size-1;i>=0;i--)
    {
      if(value - cumul < (int)POW(2,i))
        {
          bit_vect[i] = 0;
        }
      else
        {
          bit_vect[i] = 1;
          cumul += (int)POW(2,i);
        }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Print_Diversity_Header(FILE *fp, t_tree *tree)
{
  /*   PhyML_Fprintf(fp,"t_edge side mean\n");  */
  PhyML_Fprintf(fp,"t_edge side diversity count\n");
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* Estimation of density using kernel smoothing.
- where : point where I want to estimate the density,
- x : data vector,
- sample_size :  number of data points in x
*/
phydbl Univariate_Kernel_Density_Estimate(phydbl where, phydbl *x, int sample_size)
{
  phydbl sd,h;
  phydbl density,sqrt2pi,cons;
  int i;

  sqrt2pi = 2.506628;

  sd = SQRT(Var(x,sample_size));
  h = 1.06 * sd * POW(sample_size,-1./5.); /* Quick and dirty way to set the bandwidth */

  cons = (1./sample_size) * (1./h) * (1./sqrt2pi);

  density = .0;
  for(i=0;i<sample_size;i++) density += exp(-0.5 * POW((x[i] - where)/h,2));
  density *= cons;

  return density;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


/* Estimation of a multivariate density using kernel smoothing.

- where : vector where I want to estimate the density,
- x : data matrix, i.e., sample of vectors,
- sample_size : number of vectors,
- vect_size : vector length.

See "Multivariate Density Estimation" by David Scott. pp 150.
*/
phydbl Multivariate_Kernel_Density_Estimate(phydbl *where, phydbl **x, int sample_size, int vect_size)
{
  phydbl sd,*h,cons,density,tmp;
  phydbl _2pi;
  int i,j;

  h = (phydbl *)mCalloc(vect_size,sizeof(phydbl));

  _2pi = 6.283185;

  for(i=0;i<vect_size;i++)
    {
      sd = SQRT(Var(x[i],sample_size));
/*       h[i] = POW(4./(vect_size+2.),1./(vect_size+4)) * sd * POW(sample_size,-1./(vect_size+4)); */
      h[i] = sd * POW(sample_size,-1./(vect_size+4));
/*       PhyML_Printf("\n. sd = %f, h[i] = %f",sd,h[i]); */
    }

  cons = sample_size;
  for(i=0;i<vect_size;i++) cons *= h[i];
  cons *= POW(_2pi,vect_size/2.);
  cons = 1./cons;

  density = .0;
  for(i=0;i<sample_size;i++)
    {
      tmp = 1.0;
      for(j=0;j<vect_size;j++)
    {
      tmp *= exp(-0.5 * POW((x[j][i] - where[j])/h[j],2));
    }
      density += tmp;
    }

  density *= cons;

  Free(h);

  return density;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Var(phydbl *x, int n)
{
  phydbl mean, sum2;
  int i;

  mean = Mean(x,n);

  sum2 = .0;
  for(i=0;i<n;i++) sum2 += x[i] * x[i];

  return (1./n) * (sum2 - n * POW(mean,2));
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Mean(phydbl *x, int n)
{
  int i;
  phydbl sum;

  sum = .0;

  for(i=0;i<n;i++) sum += x[i];

  return sum / n;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Best_Of_NNI_And_SPR(t_tree *tree)
{
  if(tree->mod->s_opt->random_input_tree)
      Speed_Spr_Loop(tree); /* Don't do simultaneous NNIs if starting tree is random */
  else
    {
      t_tree *ori_tree,*best_tree;
      t_mod *ori_mod,*best_mod;
      scalar_dbl  **ori_bl,**best_bl;
      phydbl best_lnL,ori_lnL,nni_lnL,spr_lnL;
      int i;
#ifdef BEAGLE
      tree->b_inst = create_beagle_instance(tree, tree->io->quiet, tree->io);
#endif

      ori_mod   = Copy_Model(tree->mod);
      best_mod  = Copy_Model(tree->mod);

      ori_tree  = Make_Tree_From_Scratch(tree->n_otu,tree->data);
      best_tree = Make_Tree_From_Scratch(tree->n_otu,tree->data);

      Copy_Tree(tree,ori_tree);//Save a backup of the original tree in ori_tree
      Record_Br_Len(tree);
      ori_bl = Copy_Br_Len(tree);


      best_lnL = UNLIKELY;
      Lk(NULL,tree);
      ori_lnL = tree->c_lnL; /* Record likelihood of the starting tree */

      // ****** Perform NNI ******
      Simu_Loop(tree); /* Perform simultaneous NNIs */
      best_lnL = tree->c_lnL; /* Record the likelihood */
      nni_lnL = tree->c_lnL;

      //Mark the NNI tree as the "best" tree
      Copy_Tree(tree,best_tree); /* Record the tree topology and branch lengths */
      Record_Br_Len(tree);
      best_bl = Copy_Br_Len(tree);
      Transfer_Br_Len_To_Tree(best_bl,best_tree);
      Record_Model(tree->mod,best_mod);

      Copy_Tree(ori_tree,tree); /* Back to the original tree topology */
      Transfer_Br_Len_To_Tree(ori_bl,tree); /* Back to the original branch lengths */
      Record_Model(ori_mod,tree->mod); /* Back to the original model */

      /* Make sure the tree is in its original form */
      Lk(NULL,tree);
      if(FABS(tree->c_lnL - ori_lnL) > tree->mod->s_opt->min_diff_lk_local)
        {
          PhyML_Printf("\n. ori_lnL = %f, c_lnL = %f",ori_lnL,tree->c_lnL);
          PhyML_Printf("\n. Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
          Warn_And_Exit("\n. PhyML finished prematurely.");
        }


      // ****** Perform SPR ******
      Speed_Spr_Loop(tree);
      spr_lnL = tree->c_lnL;


      //Did SPR perform better than NNI?
      if(tree->c_lnL > best_lnL)
        {
#ifdef BEAGLE
          finalize_beagle_instance(best_tree);//Free the old BEAGLE instance associated with the NNI tree (since SPR is better)
#endif
          best_lnL = spr_lnL;
          Copy_Tree(tree,best_tree); /* Record tree topology, branch lengths and model parameters */
          Record_Br_Len(tree);
          For(i,2*tree->n_otu-1) Free_Scalar_Dbl(best_bl[i]);
          Free(best_bl);
          best_bl = Copy_Br_Len(tree);
          Transfer_Br_Len_To_Tree(best_bl,best_tree);
          Record_Model(tree->mod,best_mod);
        }

      Copy_Tree(best_tree,tree);
      Init_Partial_Lk_Tips_Double(tree);
      Init_Ui_Tips(tree);
      Init_Partial_Pars_Tips(tree);
      Transfer_Br_Len_To_Tree(best_bl,tree);
      Record_Model(best_mod,tree->mod);

      /* Make sure the current tree has the best topology, branch lengths and model parameters */
      Lk(NULL,tree);
      if(FABS(tree->c_lnL - best_lnL) > tree->mod->s_opt->min_diff_lk_local)
        {
          PhyML_Fprintf(stderr,"\n. best_lnL = %f, c_lnL = %f",best_lnL,tree->c_lnL);
          PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d\n",__FILE__,__LINE__);
          Warn_And_Exit("\n. PhyML finished prematurely.");
        }

      if(tree->verbose > VL0)
        {
          PhyML_Printf("\n\n. Log likelihood obtained after NNI moves : %f",nni_lnL);
          PhyML_Printf("\n. Log likelihood obtained after SPR moves : %f",spr_lnL);
        }

      For(i,2*tree->n_otu-1) Free_Scalar_Dbl(ori_bl[i]);
      Free(ori_bl);

      For(i,2*tree->n_otu-1) Free_Scalar_Dbl(best_bl[i]);
      Free(best_bl);

      Free_Tree(ori_tree);
      Free_Tree(best_tree);

      Free_Model_Complete(ori_mod);
      Free_Model_Complete(best_mod);

      M4_Free_M4_Model(ori_mod->m4mod);
      M4_Free_M4_Model(best_mod->m4mod);

      Free_Model_Basic(ori_mod);
      Free_Model_Basic(best_mod);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


/* Polynomial interpolation. Adapted from "Numerical Recipes in C".
Press, Flannery, Teukolsky, Vetterling, 1988.
*/
int Polint(phydbl *xa, phydbl *ya, int n, phydbl x, phydbl *y, phydbl *dy)
{
   int i,m,ns=1;
   phydbl den,dif,dift,ho,hp,w;
   phydbl *c,*d;

   dif=FABS(x-xa[1]);

   c = (phydbl *)mCalloc(n,sizeof(phydbl));
   d = (phydbl *)mCalloc(n,sizeof(phydbl));

   for(i=1;i<=n;i++)
     {
       if((dift=FABS(x-xa[i])) < dif)
     {
       ns=i;
       dif=dift;
     }
       c[i]=ya[i];
       d[i]=ya[i];
     }

   *y=ya[ns--];

   for (m=1;m<n;m++)
     {
       for (i=1;i<=n-m;i++)
     {
       ho=xa[i]-x;
       hp=xa[i+m]-x;
       w=c[i+1]-d[i];
       if((den=ho-hp) < SMALL && (den=ho-hp) > -SMALL )
         {
/* 	       Rprintf("\n. Error in routine POLINT.\n"); */
           Exit("\n. Error in routine POLINT.\n");
           return(-1);
         }
       den=w/den;
       d[i]=hp*den;
       c[i]=ho*den;
     }
       *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
     }

   Free(d);
   Free(c);
   return(0);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void JF(t_tree *tree)
{
  //printing loglk for each site, to compute SH-like tests */
  phydbl sum=0.0;
  PhyML_Printf("\n\nSITES LKS:\n");
  int n_patterns = (int)FLOOR(tree->n_pattern);
  int site=0;
  for(site=0;site<n_patterns;site++) {
    int wei=0;
    For(wei,tree->data->wght[site]) {
      PhyML_Printf("%f\n",tree->c_lnL_sorted[site]);
      sum+=tree->c_lnL_sorted[site];
    }
  }

  PhyML_Printf("\n\nsum=%f\n\n",sum);
  int i=0;
  For(i,2*tree->n_otu-3)
    {
      if((!tree->a_edges[i]->left->tax) && (!tree->a_edges[i]->rght->tax))
        {
          PhyML_Printf("%3d %f %f %f\n",
                       tree->a_edges[i]->bip_score,tree->a_edges[i]->alrt_statistic, tree->a_edges[i]->ratio_test,tree->a_edges[i]->l->v);
        }
    }


/*   //printing loglk for each site, to compute SH-like tests */
/*   phydbl sum=0.0; */
/*   PhyML_Printf("\n\nSITES LKS:\n"); */
/*   int n_patterns = (int)FLOOR(tree->n_pattern); */
/*   int site=0; */
/*   for(site=0;site<n_patterns;site++) { */
/*     int wei=0; */
/*     For(wei,tree->data->wght[site]) { */
/*       PhyML_Printf("%f\n",tree->c_lnL_sorted[site] / tree->data->wght[site]); */
/*       sum+=tree->c_lnL_sorted[site] / tree->data->wght[site]; */
/*     } */
/*   } */

/*   PhyML_Printf("\n\nsum=%f\n\n",sum); */

/*   int i=0; */
/*   For(i,2*tree->n_otu-3) */
/*     { */
/*       if((!tree->a_edges[i]->left->tax) && (!tree->a_edges[i]->rght->tax)) */
/* 	{ */
/* 	  PhyML_Printf("%3d %f %f %f\n", */
/* 		 tree->a_edges[i]->bip_score,tree->a_edges[i]->alrt_statistic, tree->a_edges[i]->ratio_test,tree->a_edges[i]->l->v); */
/* 	} */
/*     } */
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

t_tree *Dist_And_BioNJ(calign *cdata, t_mod *mod, option *io)
{
  t_tree *tree;
  matrix *mat;

  if(mod->s_opt->random_input_tree == NO)
    {
      if(!io->quiet) PhyML_Printf("\n. Computing pairwise distances...");

      mat = ML_Dist(cdata,mod);
      Fill_Missing_Dist(mat);

      if(!io->quiet) PhyML_Printf("\n\n. Building BioNJ tree...");
      mat->tree = Make_Tree_From_Scratch(cdata->n_otu,cdata);      
      Bionj(mat);
      
      tree      = mat->tree;
      tree->mat = mat;
    }
  else
    {
      tree = Make_Tree_From_Scratch(cdata->n_otu,cdata);
      tree->mat = NULL;
    }

  return tree;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Add_BioNJ_Branch_Lengths(t_tree *tree, calign *cdata, t_mod *mod, matrix *mat)
{
  short unsigned int freemat = NO;
  if(mat == NULL) freemat = YES;
  Connect_CSeqs_To_Nodes(cdata,mod->io,tree);
  if(mat == NULL) mat = ML_Dist(cdata,mod);
  mat->tree = tree;
  mat->method = 0;
  Bionj_Br_Length(mat);
  if(freemat == YES) Free_Mat(mat);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

char *Bootstrap_From_String(char *s_tree, calign *cdata, t_mod *mod, option *io)
{
  t_tree *tree;

  tree = Read_Tree(&s_tree);

  tree->n_root = NULL;
  tree->e_root = NULL;

  if(!tree)
    {
      PhyML_Printf("\n. Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
      Exit("");
    }

  tree->mod         = mod;
  tree->io          = io;
  tree->data        = cdata;
  tree->n_pattern   = tree->data->crunch_len;

  Connect_CSeqs_To_Nodes(cdata,io,tree);
  if(tree->mod->s_opt->random_input_tree) Random_Tree(tree);
  Make_Tree_4_Pars(tree,cdata,cdata->init_len);
  Make_Tree_4_Lk(tree,cdata,cdata->init_len);
  tree->triplet_struct = Make_Triplet_Struct(mod);
  Init_Triplet_Struct(tree->triplet_struct);
  Unscale_Br_Len_Multiplier_Tree(tree);
  Br_Len_Not_Involving_Invar(tree);
  Make_Spr_List(tree);
  Make_Best_Spr(tree);

  Set_Both_Sides(YES,tree);
  Lk(NULL,tree);

#ifdef MPI
  Bootstrap_MPI(tree);
#else
  Bootstrap(tree);
#endif

  Free(s_tree);

  Rescale_Br_Len_Multiplier_Tree(tree);
  Br_Len_Involving_Invar(tree);

  s_tree = Write_Tree(tree,NO);

  Free_Spr_List(tree);
  Free_Triplet(tree->triplet_struct);
  Free_Tree_Pars(tree);
  Free_Tree_Lk(tree);
  Free_Tree(tree);

  return s_tree;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

char *aLRT_From_String(char *s_tree, calign *cdata, t_mod *mod, option *io)
{
  t_tree *tree;

  tree = Read_Tree(&s_tree);

  tree->n_root = NULL;
  tree->e_root = NULL;

  if(!tree)
    {
      PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
      Warn_And_Exit("\n. PhyML finished prematurely.");
    }

  tree->mod         = mod;
  tree->io          = io;
  tree->data        = cdata;
  tree->n_pattern   = tree->data->crunch_len;

  Connect_CSeqs_To_Nodes(cdata,io,tree);
  if(tree->mod->s_opt->random_input_tree) Random_Tree(tree);
  Make_Tree_4_Pars(tree,cdata,cdata->init_len);
  Make_Tree_4_Lk(tree,cdata,cdata->init_len);
  tree->triplet_struct = Make_Triplet_Struct(mod);
  Init_Triplet_Struct(tree->triplet_struct);

  Unscale_Br_Len_Multiplier_Tree(tree);
  Br_Len_Not_Involving_Invar(tree);

  Make_Spr_List(tree);
  Make_Best_Spr(tree);

#ifdef BEAGLE
  tree->b_inst = create_beagle_instance(tree, io->quiet, io);
#endif

  Set_Both_Sides(YES,tree);
  Lk(NULL,tree);
  
  aLRT(tree);

  Free(s_tree);

  Rescale_Br_Len_Multiplier_Tree(tree);
  Br_Len_Involving_Invar(tree);

  s_tree = Write_Tree(tree,NO);

#ifdef BEAGLE
  finalize_beagle_instance(tree);
#endif

  Free_Spr_List(tree);
  Free_Triplet(tree->triplet_struct);
  Free_Tree_Pars(tree);
  Free_Tree_Lk(tree);
  Free_Tree(tree);

  return s_tree;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Prepare_Tree_For_Lk(t_tree *tree)
{
  Connect_CSeqs_To_Nodes(tree->data,tree->io,tree);
  Make_Tree_4_Pars(tree,tree->data,tree->data->init_len);
  Make_Tree_4_Lk(tree,tree->data,tree->data->init_len);
  tree->triplet_struct = Make_Triplet_Struct(tree->mod);
  Init_Triplet_Struct(tree->triplet_struct);
  Make_Spr_List(tree);
  Make_Best_Spr(tree);

  if(tree->is_mixt_tree == YES) MIXT_Prepare_Tree_For_Lk(tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Find_Common_Tips(t_tree *tree1, t_tree *tree2)
{
  int i,j;

  for(i=0;i<tree1->n_otu;i++) tree1->a_nodes[i]->common = 0;
  for(i=0;i<tree2->n_otu;i++) tree2->a_nodes[i]->common = 0;

  for(i=0;i<tree1->n_otu;i++)
    {
      for(j=0;j<tree2->n_otu;j++)
    {
      if(!strcmp(tree1->a_nodes[i]->name,tree2->a_nodes[j]->name))
        {
          tree1->a_nodes[i]->common = 1;
          tree2->a_nodes[j]->common = 1;
          break;
        }
    }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl Get_Tree_Size(t_tree *tree)
{
  int i;
  phydbl tree_size;

  tree_size = 0.0;
  For(i,2*tree->n_otu-3) tree_size += tree->a_edges[i]->l->v;

  if(tree->n_root != NULL)
    {
      tree_size += tree->n_root->b[1]->l->v;
      tree_size += tree->n_root->b[2]->l->v;
    }


/*   For(i,2*tree->n_otu-3)  */
/*     tree_size +=  */
/*     FABS(tree->rates->nd_t[tree->a_edges[i]->left->num] -  */
/* 	 tree->rates->nd_t[tree->a_edges[i]->rght->num]); */

  tree->size = tree_size;
  return tree_size;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Dist_To_Root_Pre(t_node *a, t_node *d, t_edge *b, t_tree *tree)
{
  int i;

  if(b) d->dist_to_root = a->dist_to_root + b->l->v;


  if(d->tax) return;
  else
    {
      for(i=0;i<3;i++)
    if((d->v[i] != a) && (d->b[i] != tree->e_root))
      Dist_To_Root_Pre(d,d->v[i],d->b[i],tree);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Dist_To_Root(t_node *n_root, t_tree *tree)
{
  /* n_root->v[2]->dist_to_root = tree->rates->cur_l[n_root->v[2]->num]; */
  /* n_root->v[1]->dist_to_root = tree->rates->cur_l[n_root->v[1]->num]; */
  n_root->v[2]->dist_to_root = tree->e_root->l->v * tree->n_root_pos;
  n_root->v[1]->dist_to_root = tree->e_root->l->v * (1. - tree->n_root_pos);
  Dist_To_Root_Pre(n_root,n_root->v[2],NULL,tree);
  Dist_To_Root_Pre(n_root,n_root->v[1],NULL,tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* 'Borrowed' fromn libgen */
char *Basename(char *path)
{
  char *p;

  if( path == NULL || *path == '\0' ) return ".";

  p = path + strlen(path) - 1;

  while( *p == '/' )
    {
      if( p == path ) return path;
      *p-- = '\0';
    }

  while( p >= path && *p != '/' ) p--;

  return p + 1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* Find the Last Common Ancestor of n1 and n2 */
t_node *Find_Lca_Pair_Of_Nodes(t_node *n1, t_node *n2, t_tree *tree)
{
  t_node **list1, **list2, *lca;
  int size1, size2;

  if(n1 == n2) return(n1);

  if(!tree->n_root)
    {
      PhyML_Printf("\n. The tree must be rooted in this function.");
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);    
    }

  list1 = (t_node **)mCalloc(2*tree->n_otu-1,sizeof(t_node *));
  list2 = (t_node **)mCalloc(2*tree->n_otu-1,sizeof(t_node *));

  Get_List_Of_Ancestors(n1,list1,&size1,tree);
  Get_List_Of_Ancestors(n2,list2,&size2,tree);

  while(list1[size1] == list2[size2])
    {
      size1--;
      size2--;

      if(size1 < 0 || size2 < 0) break;
    }

  lca = list1[size1+1];

  Free(list1);
  Free(list2);

  if(lca == NULL)
    {
      PhyML_Printf("\n. %s",Write_Tree(tree,NO));
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
    }
  return lca;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* Find the Last Common Ancestor of all the nodes in node_list */
t_node *Find_Lca_Clade(t_node **node_list, int node_list_size, t_tree *tree)
{
  t_node ***list, *lca;
  int *size;
  int i;

  if(!tree->n_root)
    {
      PhyML_Printf("\n. The tree must be rooted in this function.");
      PhyML_Printf("\n. Err. in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("\n. PhyML finished prematurely.");
    }

  list = (t_node ***)mCalloc(node_list_size,sizeof(t_node **));
  for(i=0;i<node_list_size;i++) list[i] = (t_node **)mCalloc(2*tree->n_otu-1,sizeof(t_node *));
  size = (int *)mCalloc(node_list_size,sizeof(int));
  
  for(i=0;i<node_list_size;i++)
    {
      if(!Get_List_Of_Ancestors(node_list[i],list[i],size+i,tree))
        {
          for(i=0;i<node_list_size;i++) PhyML_Printf("\n. %s",node_list[i]->name);
          Generic_Exit(__FILE__,__LINE__,__FUNCTION__);    
        }
    }

  /* for(i=0;i<node_list_size;i++) */
  /*   { */
  /*     int j; */
  /*     PhyML_Printf("\n. Listing all ancestors of node number %d [%s]", */
  /*                  node_list[i]->num, */
  /*                  node_list[i]->tax ? node_list[i]->name : NULL); */
  /*     For(j,size[i]) PhyML_Printf("\n.  > %d <",list[i][j]->num); */
  /*   } */

  if(node_list_size > 1)
    {
      do
        {
          for(i=0;i<node_list_size-1;i++)
            {
              assert(list[i][size[i]-1]);
              assert(list[i+1][size[i+1]-1]);
              /* PhyML_Printf("\n. %d %d %d %d",list[i][size[i]-1]->num,size[i],list[i+1][size[i+1]-1]->num,size[i+1]); */
              if(list[i][size[i]-1] != list[i+1][size[i+1]-1])
                {
                  /* PhyML_Printf("\n. Break at %d %d",list[i][size[i]]->num,list[i+1][size[i+1]]->num); */
                  break;
                }
            }

          if(i != node_list_size-1) break;
          
          for(i=0;i<node_list_size;i++)
            {
              size[i]--;
              assert(size[i] > 0);
            }
          
          if(node_list_size == 1) break;
          
        }while(1);
      lca = list[0][size[0]];
    }
  else
    {
      lca = node_list[0];
    }
  
  for(i=0;i<node_list_size;i++) Free(list[i]);
  Free(list);
  Free(size);

  /* PhyML_Printf("\n. LCA: %d",lca->num); */

  return lca;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* Returns the list of the ancestors of ref_t_node from ref_t_node to the root included */
int Get_List_Of_Ancestors(t_node *ref_node, t_node **list, int *size, t_tree *tree)
{
  t_node *n;
  
  n = ref_node;
  list[0] = n;  
  *size = 1;

  if(!n)
    {
      PhyML_Printf("\n. There seems to be a problem with the calibration file.\n");
      return 0;
    }
  
  while(n != tree->n_root)
    {
      n = n->anc;
      if(!n)
        {
          PhyML_Printf("\n. n->anc has not been set properly (call Update_Ancestors first...)\n");
          return 0;
        }
      list[*size] = n;
      *size = *size+1;
    }
  return 1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int Edge_Num_To_Node_Num(int edge_num, t_tree *tree)
{
  int node_num;
  t_edge *b;

  b = tree->a_edges[edge_num];

  node_num = (b->left == b->rght->anc)?(b->rght->num):(b->left->num);

  return node_num;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Time_To_Branch(t_tree *tree)
{
  Time_To_Branch_Pre(tree->n_root,tree->n_root->v[2],tree);
  Time_To_Branch_Pre(tree->n_root,tree->n_root->v[1],tree);
  tree->n_root->l[1] = tree->rates->nd_t[tree->n_root->v[1]->num] - tree->rates->nd_t[tree->n_root->num];
  tree->n_root->l[2] = tree->rates->nd_t[tree->n_root->v[1]->num] - tree->rates->nd_t[tree->n_root->num];
  tree->e_root->l->v = tree->n_root->l[1] + tree->n_root->l[2];
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Time_To_Branch_Pre(t_node *a, t_node *d, t_tree *tree)
{
  int i;
  
  /* tree->rates->cur_l[d->num] = FABS(tree->rates->nd_t[d->num] - tree->rates->nd_t[a->num]); */
  tree->rates->cur_l[d->num] = tree->rates->nd_t[d->num] - tree->rates->nd_t[a->num];
  
  if(d->tax) return;
  else
    {
      for(i=0;i<3;i++)
        if((d->v[i] != a) && (d->b[i] != tree->e_root))
          Time_To_Branch_Pre(d,d->v[i],tree);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

// Assume an ultrametric tree.
void Branch_To_Time(t_tree *tree)
{
  Branch_To_Time_Pre(tree->n_root,tree->n_root->v[2],tree);
  Branch_To_Time_Pre(tree->n_root,tree->n_root->v[1],tree);

  tree->rates->nd_t[tree->n_root->num] =
    tree->rates->nd_t[tree->n_root->v[1]->num] -
    tree->n_root->l[1];
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Branch_To_Time_Pre(t_node *a, t_node *d, t_tree *tree)
{
  int i;

  if(d->tax)
    {
      tree->rates->nd_t[d->num] = 0.0;
      return;
    }
  else
    {
      for(i=0;i<3;i++)
    if((d->v[i] != a) && (d->b[i] != tree->e_root))
      {
            Branch_To_Time_Pre(d,d->v[i],tree);
          }

      for(i=0;i<3;i++)
    if((d->v[i] != a) && (d->b[i] != tree->e_root))
      {
            tree->rates->nd_t[d->num] = tree->rates->nd_t[d->v[i]->num] - d->b[i]->l->v;
            break;
          }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Branch_Lengths_To_Rate_Lengths(t_tree *tree)
{
  Branch_Lengths_To_Rate_Lengths_Pre(tree->n_root,tree->n_root->v[2],tree);
  Branch_Lengths_To_Rate_Lengths_Pre(tree->n_root,tree->n_root->v[1],tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Branch_Lengths_To_Rate_Lengths_Pre(t_node *a, t_node *d, t_tree *tree)
{
  int i;

  tree->rates->cur_l[d->num] =
    tree->rates->br_r[d->num] *
    tree->rates->clock_r *
    tree->rates->norm_fact;

  if(d->tax) return;
  else
    {
      for(i=0;i<3;i++)
    if((d->v[i] != a) && (d->b[i] != tree->e_root))
      Branch_Lengths_To_Rate_Lengths_Pre(d,d->v[i],tree);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int Find_Clade(char **tax_name_list, int list_size, t_tree *tree)
{
  int *tax_num_list;
  t_node **tax_node_list;
  int i,j;
  int n_matches;
  t_node *lca;

  tax_num_list = (int *)mCalloc(list_size,sizeof(int));
  tax_node_list = (t_node **)mCalloc(list_size,sizeof(t_node *));

  for(i=0;i<list_size;i++) tax_num_list[i] = -1;

  n_matches = 0;

  for(i=0;i<list_size;i++)
    {
      for(j=0;j<tree->n_otu;j++)
        {
          if(!strcmp(tax_name_list[i],tree->a_nodes[j]->name))
            {
              tax_num_list[i] = tree->a_nodes[j]->num;
              tax_node_list[i] = tree->a_nodes[j];
              n_matches++;
              break;
            }
        }
      
      if(j == tree->n_otu)
        {
          PhyML_Printf("\n. Problem with the calibration file.");
          PhyML_Printf("\n. Could not find taxon with name '%s' in the sequence or tree file.",tax_name_list[i]);
          Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
        }
    }

  lca = Find_Lca_Clade(tax_node_list,n_matches,tree);

  Free(tax_num_list);
  Free(tax_node_list);

  if(lca) return lca->num;
  else    return -1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Find_Clade_Pre(t_node *a, t_node *d, int *tax_num_list, int list_size, int *num, t_tree *tree)
{
  int i,j,k;
  int score;


  for(i=0;i<3;i++)
    if((d->v[i] == a) || (d->b[i] == tree->e_root))
      {
    if(list_size == d->bip_size[i])
      {
        score = 0;
        For(j,d->bip_size[i])
          {
        for(k=0;k<list_size;k++)
          {
            if(tax_num_list[k] == d->bip_node[i][j]->num)
              {
            score++;
            break;
              }
          }
          }
        if(score == list_size) *num = d->num;
      }
    break;
      }

  if(d->tax) return;
  else
    for(i=0;i<3;i++)
      if((d->v[i] != a) && (d->b[i] != tree->e_root))
        Find_Clade_Pre(d,d->v[i],tax_num_list,list_size,num,tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

t_edge *Find_Root_Edge(FILE *fp_input_tree, t_tree *tree)
{
  char **subs;
  int degree;
  int i,j;
  t_node *left, *rght;
  int l_r, r_l;
  int score;
  char *line;
  char c;
  t_edge *root_edge;

  line = (char *)mCalloc(T_MAX_LINE,sizeof(char));

  rewind(fp_input_tree);

  do c=fgetc(fp_input_tree);
  while((c != '(') && (c != EOF));

  if(c==EOF)
    {
      Free(line);
      return NULL;
    }

  i=0;
  for(;;)
    {
      if((c == ' ') || (c == '\n'))
    {
      c=fgetc(fp_input_tree);
      if(c==EOF) break;
      else continue;
    }

      line[i]=c;
      i++;
      c=fgetc(fp_input_tree);
      if(c==EOF || c==';') break;
    }


  Free_Bip(tree);
  Alloc_Bip(tree);
  Get_Bip(tree->a_nodes[0],tree->a_nodes[0]->v[0],tree);

  subs = Sub_Trees(line,&degree);
  Clean_Multifurcation(subs,degree,3);
  if(degree != 2)
    {
      PhyML_Printf("\n. The tree does not seem to be rooted...");
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("\n. PhyML finished prematurely.");
    }

  left = rght = NULL;
  l_r = r_l = -1;

  For(i,2*tree->n_otu-3)
    {
      left = tree->a_edges[i]->left;
      rght = tree->a_edges[i]->rght;
      l_r  = tree->a_edges[i]->l_r;
      r_l  = tree->a_edges[i]->r_l;

      score = 0;
      For(j,left->bip_size[l_r]) if(strstr(subs[1],left->bip_node[l_r][j]->name)) score++;
      if(score == left->bip_size[l_r]) break;

      score = 0;
      For(j,rght->bip_size[r_l]) if(strstr(subs[1],rght->bip_node[r_l][j]->name)) score++;
      if(score == rght->bip_size[r_l]) break;
    }

  root_edge = tree->a_edges[i];

  i = 0;
  while(subs[i] != NULL) Free(subs[i++]);
  Free(subs);
  Free(line);

  if(i == 2*tree->n_otu-3)
    {
      PhyML_Printf("\n. Could not find the root edge...");
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("\n. PhyML finished prematurely.");
    }

  return root_edge;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Copy_Tree_Topology_With_Labels(t_tree *ori, t_tree *cpy)
{
  int i,j;

  For(i,2*ori->n_otu-2)
    {
      for(j=0;j<3;j++)
        {
          if(ori->a_nodes[i]->v[j])
            {
              cpy->a_nodes[i]->v[j] = cpy->a_nodes[ori->a_nodes[i]->v[j]->num];
              cpy->a_nodes[i]->l[j] = ori->a_nodes[i]->l[j];
            }
          else
            cpy->a_nodes[i]->v[j] = NULL;
        }
      cpy->a_nodes[i]->num = ori->a_nodes[i]->num;
      cpy->a_nodes[i]->tax = 0;
    }

  For(i,2*ori->n_otu-3)
    {
      cpy->a_edges[i]->l->v = ori->a_edges[i]->l->v;
    }

  for(i=0;i<ori->n_otu;i++)
    {
      cpy->a_nodes[i]->tax = 1;
      strcpy(cpy->a_nodes[i]->name,ori->a_nodes[i]->name);
    }

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Set_Model_Name(t_mod *mod)
{
  switch(mod->whichmodel)
    {
    case JC69:
      {
    strcpy(mod->modelname->s, "JC69");
    break;
      }
    case K80:
      {
    strcpy(mod->modelname->s, "K80");
    break;
      }
    case F81:
      {
    strcpy(mod->modelname->s, "F81");
    break;
      }
    case HKY85:
      {
    strcpy(mod->modelname->s, "HKY85");
    break;
      }
    case F84:
      {
    strcpy(mod->modelname->s, "F84");
    break;
      }
    case TN93:
      {
    strcpy(mod->modelname->s, "TN93");
    break;
      }
    case GTR:
      {
    strcpy(mod->modelname->s, "GTR");
    break;
      }
    case CUSTOM:
      {
    strcpy(mod->modelname->s, "Custom");
    break;
      }
    case DAYHOFF:
      {
    strcpy(mod->modelname->s, "Dayhoff");
    break;
      }
    case JTT:
      {
    strcpy(mod->modelname->s, "JTT");
    break;
      }
    case MTREV:
      {
    strcpy(mod->modelname->s, "MtREV");
    break;
      }
    case LG:
      {
    strcpy(mod->modelname->s, "LG");
    break;
      }
    case WAG:
      {
    strcpy(mod->modelname->s, "WAG");
    break;
      }
    case DCMUT:
      {
    strcpy(mod->modelname->s, "DCMut");
    break;
      }
    case RTREV:
      {
    strcpy(mod->modelname->s, "RtREV");
    break;
      }
    case CPREV:
      {
    strcpy(mod->modelname->s, "CpREV");
    break;
      }
    case VT:
      {
    strcpy(mod->modelname->s, "VT");
    break;
      }
    case BLOSUM62:
      {
    strcpy(mod->modelname->s, "Blosum62");
    break;
      }
    case MTMAM:
      {
    strcpy(mod->modelname->s, "MtMam");
    break;
      }
    case MTART:
      {
    strcpy(mod->modelname->s, "MtArt");
    break;
      }
    case HIVW:
      {
    strcpy(mod->modelname->s, "HIVw");
    break;
      }
    case HIVB:
      {
    strcpy(mod->modelname->s, "HIVb");
    break;
      }
    case AB:
      {
    strcpy(mod->modelname->s, "AB");
    break;
      }
    case CUSTOMAA:
      {
    strcpy(mod->modelname->s, "Custom");
    break;
      }
    default:
      {
    PhyML_Printf("\n. Unknown model name.\n");
    PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
    Warn_And_Exit("\n. PhyML finished prematurely.");
    break;
      }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Adjust_Min_Diff_Lk(t_tree *tree)
{
  if(sizeof(phydbl) == 4)
    {
      int exponent;
      exponent = (int)FLOOR(log10(FABS(tree->c_lnL)));
      tree->mod->s_opt->min_diff_lk_local = POW(10.,exponent - FLT_DIG + 1);
      tree->mod->s_opt->min_diff_lk_move  = tree->mod->s_opt->min_diff_lk_local;
    }
/*   PhyML_Printf("\n. Exponent = %d Precision = %E DIG = %d",exponent,tree->mod->s_opt->min_diff_lk_local,FLT_DIG); */
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


/*!
  tree->a_nodes[i]->name is initially a number. It is translated into
  a string of characters using the names provided in the tax_name
  array.
 */
void Translate_Tax_Names(char **tax_names, t_tree *tree)
{
  int i;
  int tax_num;

  for(i=0;i<tree->n_otu;i++)
    {
      tax_num = strtol(tree->a_nodes[i]->name,NULL,10);
      tree->a_nodes[i]->name = tax_names[tax_num-1];
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


/*!
  Skip coment in NEXUS file.
 */
void Skip_Comment(FILE *fp)
{
  int in_comment;
  char c;

  in_comment = 1;
  do
    {
      c = fgetc(fp);
      if(c == EOF) break;
      if(c == '[')      in_comment++;
      else if(c == ']') in_comment--;
    }
  while(in_comment);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/*!
  Determine the most appropriate position of the root if outgroup taxa are specified.
 */

void Get_Best_Root_Position(t_tree *tree)
{
  int i,j;
  phydbl eps;
  phydbl s, s_max;
  t_edge *best_edge;
  int has_outgrp;

  best_edge = NULL;

  if(tree->n_root)
    {
      PhyML_Printf("\n. Tree already has a root.");
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("\n. PhyML finished prematurely.");
    }

  has_outgrp = NO;

  if(strstr(tree->a_nodes[0]->name,"*"))
    {
      /* PhyML_Printf("\n. Found outgroup taxon: %s",tree->a_nodes[0]->name); */
      tree->a_nodes[0]->s_ingrp[0]  = 0;
      tree->a_nodes[0]->s_outgrp[0] = 1;
      has_outgrp = YES;
    }
  else
    {
      tree->a_nodes[0]->s_ingrp[0]  = 1;
      tree->a_nodes[0]->s_outgrp[0] = 0;
    }

  Get_Best_Root_Position_Post(tree->a_nodes[0],tree->a_nodes[0]->v[0],&has_outgrp,tree);
  Get_Best_Root_Position_Pre(tree->a_nodes[0],tree->a_nodes[0]->v[0],tree);

  if(has_outgrp == YES)
    {

      Free_Edge_Lk_Rght(tree->a_edges[2*tree->n_otu-3]);
      Free_Edge_Lk_Rght(tree->a_edges[2*tree->n_otu-2]);
      Free_Edge_Pars_Rght(tree->a_edges[2*tree->n_otu-3]);
      Free_Edge_Pars_Rght(tree->a_edges[2*tree->n_otu-2]);

      eps = 1.E-10;
      s = s_max = 0.0;
      For(i,2*tree->n_otu-2)
        {
          for(j=0;j<3;j++)
            {
              s = (tree->a_nodes[i]->s_outgrp[j]+eps) / (tree->a_nodes[i]->s_ingrp[j] + eps) ;
              /* printf("\n. [%d %d] %d %d",i,j,tree->a_nodes[i]->s_outgrp[j],tree->a_nodes[i]->s_ingrp[j]); */
              if(s > s_max)
                {
                  s_max = s;
                  best_edge = tree->a_nodes[i]->b[j];
                }
            }
        }      
      Add_Root(best_edge,tree);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


/*!
  Determine the most appropriate position of the root if outgroup taxa are specified.
  Post-traversal.
 */
void Get_Best_Root_Position_Post(t_node *a, t_node *d, int *has_outgrp, t_tree *tree)
{
  if(d->tax)
    {
      if(strstr(d->name,"*"))
        {
          *has_outgrp = YES;
          /* PhyML_Printf("\n. Found outgroup taxon: %s",d->name); */
          d->s_ingrp[0]  = NO;
          d->s_outgrp[0] = YES;
        }
      else
        {
          d->s_ingrp[0]  = YES;
          d->s_outgrp[0] = NO;
        }
      return;
    }
  else
    {
      int i;
      
      for(i=0;i<3;i++)
        if(d->v[i] != a && (d->b[i] != tree->e_root))
          Get_Best_Root_Position_Post(d,d->v[i],has_outgrp,tree);
      
      Get_OutIn_Scores(a,d);
      
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


/*!
  Determine the most appropriate position of the root if outgroup taxa are specified.
  Pre-traversal.
 */
void Get_Best_Root_Position_Pre(t_node *a, t_node *d, t_tree *tree)
{
  if(d->tax)
    {
      return;
    }
  else
    {
      int i;

      for(i=0;i<3;i++)
    if(d->v[i] != a && (d->b[i] != tree->e_root))
      {
        Get_OutIn_Scores(d->v[i],d);
        Get_Best_Root_Position_Pre(d,d->v[i],tree);
      }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


/*!
  Determine the most appropriate position of the root if outgroup taxa are specified.
  Core.
 */
void Get_OutIn_Scores(t_node *a, t_node *d)
{
  int i,d_v1,d_v2,v1_d,v2_d,d_a;

  d_a = v1_d = v2_d = -1;
  d_v1 = d_v2 = -1;
  for(i=0;i<3;i++)
    {
      if(d->v[i] != a)
    {
      if(d_v1 < 0) d_v1 = i;
      else         d_v2 = i;
    }
    }

  for(i=0;i<3;i++) if(d->v[i] == a) { d_a = i; break; }
  for(i=0;i<3;i++) if(d->v[d_v1]->v[i] == d) { v1_d = i; break; }
  for(i=0;i<3;i++) if(d->v[d_v2]->v[i] == d) { v2_d = i; break; }

  d->s_ingrp[d_a] =
    d->v[d_v1]->s_ingrp[v1_d] +
    d->v[d_v2]->s_ingrp[v2_d] ;

  d->s_outgrp[d_a] =
    d->v[d_v1]->s_outgrp[v1_d] +
    d->v[d_v2]->s_outgrp[v2_d] ;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int Check_Sequence_Name(char *s)
{
  int i;
  /* if(rindex(s,':')) */
  For(i,strlen(s))
    {
      if(s[i] == ':')
    {
      PhyML_Printf("\n. Character ':' is not permitted in sequence name (%s).",s);
      PhyML_Printf("\n. Err. in file %s at line %d",__FILE__,__LINE__);
      Warn_And_Exit("\n. PhyML finished prematurely.");
    }
    }
  /* if(rindex(s,',')) */
  For(i,strlen(s))
    {
      if(s[i] == ',')
    {
      PhyML_Printf("\n. Character ',' is not permitted in sequence name (%s).",s);
      PhyML_Printf("\n. Err in file %s at line %d",__FILE__,__LINE__);
      Warn_And_Exit("\n. PhyML finished prematurely.");
    }
    }
  /* if(rindex(s,' ')) */
  For(i,strlen(s))
    {
      if(s[i] == ' ')
    {
      PhyML_Printf("\n. Character ' ' is not permitted in sequence name (%s).",s);
      PhyML_Printf("\n. Err in file %s at line %d",__FILE__,__LINE__);
      Warn_And_Exit("\n. PhyML finished prematurely.");
    }
    }

  return 1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int Scale_Subtree_Height(t_node *a, phydbl K, phydbl floor, int *n_nodes, t_tree *tree)
{
  phydbl new_height;

  if(a->tax == YES) return 0;
  
  *n_nodes = 0;

  new_height = .0;

  if(!(tree->rates->nd_t[a->num] > floor)) new_height = K*(tree->rates->nd_t[a->num]-floor)+floor;

  if(a == tree->n_root)
    {
      tree->rates->nd_t[tree->n_root->num] = new_height;
      *n_nodes = 1;
      
      Scale_Node_Heights_Post(tree->n_root,tree->n_root->v[2],K,floor,n_nodes,tree);
      Scale_Node_Heights_Post(tree->n_root,tree->n_root->v[1],K,floor,n_nodes,tree);
    }
  else
    {
      int i;
      
      if(new_height < tree->rates->nd_t[a->anc->num]) return 0;
      else
        {
          tree->rates->nd_t[a->num] = new_height;
          *n_nodes = 1;
        }
      
      for(i=0;i<3;i++)
        if(a->v[i] != a->anc && a->b[i] != tree->e_root)
          {
            Scale_Node_Heights_Post(a,a->v[i],K,floor,n_nodes,tree);
          }
    }
  
  return 1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Scale_Node_Heights_Post(t_node *a, t_node *d, phydbl K, phydbl floor, int *n_nodes, t_tree *tree)
{
  if(d == tree->n_root)
    {
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);    
    }

  if(d->tax) return;
  else
    {
      int i;
      
      /* It is tempting to set floor = tree->rates->t_prior_max[d->num]; but
         it then becomes possible for nodes with different floor values
         to have their orders interverted (i.e., ancestor below descendant)
      */
      if((tree->rates->nd_t[d->num] > floor) == NO) // If node is strictly older than floor
        {
          tree->rates->nd_t[d->num] = K*(tree->rates->nd_t[d->num]-floor)+floor;
          *n_nodes = *n_nodes+1;
        }
      
      if(tree->rates->nd_t[d->num] < tree->rates->nd_t[a->num])
        {
          PhyML_Printf("\n. K = %f floor = %f t_prior_max(a) = %f t_prior_max(d) = %f a->t = %f d->t %f",
                       K,floor,tree->rates->t_prior_max[a->num],tree->rates->t_prior_max[d->num],
                       tree->rates->nd_t[a->num],tree->rates->nd_t[d->num]);
          PhyML_Printf("\n. Err. in file %s at line %d\n",__FILE__,__LINE__);
          Warn_And_Exit("\n. PhyML finished prematurely.");
        }
      
      for(i=0;i<3;i++)
        if(d->v[i] != a && d->b[i] != tree->e_root)
          Scale_Node_Heights_Post(d,d->v[i],K,floor,n_nodes,tree);
      
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int Scale_Subtree_Rates(t_node *a, phydbl mult, int *n_nodes, t_tree *tree)
{
  int res;
  int i;
    
  *n_nodes = 0;
  res      = 1;

  if(a == tree->n_root)
    {
      res = Scale_Subtree_Rates_Post(a,a->v[2],mult,n_nodes,tree);
      if(res) res = Scale_Subtree_Rates_Post(a,a->v[1],mult,n_nodes,tree);
      return res;
    }
  else
    {
      for(i=0;i<3;i++) if((a->v[i] != a->anc) &&
          (a->b[i] != tree->e_root) &&
          (res == 1)) res = Scale_Subtree_Rates_Post(a,a->v[i],mult,n_nodes,tree);
      return res;
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int Scale_Subtree_Rates_Post(t_node *a, t_node *d, phydbl mult, int *n_nodes, t_tree *tree)
{

  tree->rates->br_r[d->num] *= mult;

  *n_nodes = *n_nodes+1;
  
  if(tree->rates->br_r[d->num] < tree->rates->min_rate) return 0;
  if(tree->rates->br_r[d->num] > tree->rates->max_rate) return 0;
  
  if(d->tax) return 1;
  else
    {
      int i,res;
      
      res = 1;
      for(i=0;i<3;i++)
        {
          if((d->v[i] != a) &&
             (d->b[i] != tree->e_root) &&
             (res == 1))
            {
              res = Scale_Subtree_Rates_Post(d,d->v[i],mult,n_nodes,tree);
            }
        }
      return res;
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Get_Node_Ranks(t_tree *tree)
{
  tree->n_root->rank = 1;
  Get_Node_Ranks_Pre(tree->n_root,tree->n_root->v[2],tree);
  Get_Node_Ranks_Pre(tree->n_root,tree->n_root->v[1],tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Get_Node_Ranks_Pre(t_node *a, t_node *d,t_tree *tree)
{
  d->rank = a->rank+1;

  if(d->tax) return;
  else
    {
      int i;

      for(i=0;i<3;i++)
    {
      if(d->v[i] != a && d->b[i] != tree->e_root)
        {
          Get_Node_Ranks_Pre(d,d->v[i],tree);
        }
    }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Log_Br_Len(t_tree *tree)
{
  int i;
  For(i,2*tree->n_otu-3) tree->a_edges[i]->l->v = log(tree->a_edges[i]->l->v);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl Diff_Lk_Norm_At_Given_Edge(t_edge *b, t_tree *tree)
{
  int i,dim,err;
  phydbl lk_exact,lk_norm,sum;

  Record_Br_Len(tree);

  dim = 2*tree->n_otu-3;
  sum = 0.0;

  for(i=0;i<tree->n_short_l;i++)
    {
      b->l->v = tree->short_l[i];

      lk_exact = Lk(b,tree);
      lk_norm = tree->norm_scale + Log_Dnorm(b->l->v,tree->rates->mean_l[b->num],
                         tree->rates->cov_l[b->num*dim+b->num],&err);

      if(err)
        {
          Generic_Exit(__FILE__,__LINE__,__FUNCTION__);    
        }

      sum += pow(lk_exact - lk_norm,2);
    }

  Restore_Br_Len(tree);
  Lk(b,tree);

  return(sum);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Adjust_Variances(t_tree *tree)
{
  int i;
  phydbl new_diff,curr_diff;

  Make_Short_L(tree);
  for(i=0;i<tree->n_short_l;i++)
    {
      tree->short_l[i] = tree->mod->l_min + i*(0.1 - tree->mod->l_min)/tree->n_short_l;
    }


  For(i,2*tree->n_otu-3)
    {
      if(tree->a_edges[i]->l->v < 1.1*tree->mod->l_min)
        {
          tree->rates->mean_l[i]                     = -1.00;
          tree->rates->cov_l[i*(2*tree->n_otu-3)+i]  =  0.1;
          tree->norm_scale                           = -100;
          
          
          new_diff = curr_diff = 10.0;
          do
            {
              curr_diff = new_diff;
              
              Generic_Brent_Lk(&(tree->norm_scale),
                               -1E+6,
                               0.0,
                               1.E-10,
                               10000,
                               NO,
                               Wrap_Diff_Lk_Norm_At_Given_Edge,tree->a_edges[i],tree,NULL,NO);
              
              /* 		      Generic_Brent_Lk(&(tree->rates->mean_l[0]), */
              /* 				       -100., */
              /* 				       10*tree->mod->l_min, */
              /* 				       1.E-3, */
              /* 				       10000, */
              /* 				       NO, */
              /* 				       Wrap_Diff_Lk_Norm_At_Given_Edge,tree->a_edges[0],tree,NULL); */
              
              Generic_Brent_Lk(&(tree->rates->cov_l[i*(2*tree->n_otu-3)+i]),
                               0.0,
                               10.0,
                               1.E-10,
                               10000,
                               NO,
                               Wrap_Diff_Lk_Norm_At_Given_Edge,tree->a_edges[i],tree,NULL,NO);
              
              new_diff = Diff_Lk_Norm_At_Given_Edge(tree->a_edges[i],tree);
            }while(FABS(new_diff-curr_diff) > 1.E-3);
        }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Effective_Sample_Size(phydbl first_val, phydbl last_val, phydbl sum, phydbl sumsq, phydbl sumcurnext, int n)
{
  phydbl numerator,denom;
  phydbl mean;
  phydbl r;

  mean = sum / n;
  denom = sumsq - n * POW(mean,2);
  numerator = sumcurnext - (n+1.)*POW(mean,2) + (first_val+last_val)*mean;

  r = numerator/denom;

  return (phydbl)n * (1.-r)/(1.+r);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl Rescale_Br_Len_Multiplier_Tree(t_tree *tree)
{
  int i;

  if(tree->is_mixt_tree)
    {
      MIXT_Rescale_Br_Len_Multiplier_Tree(tree);
      return(-1.);
    }

  For(i,2*tree->n_otu-1) tree->a_edges[i]->l->v *= tree->mod->br_len_mult->v;
  return(-1.);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl Unscale_Br_Len_Multiplier_Tree(t_tree *tree)
{
  int i;

  if(tree->is_mixt_tree)
    {
      MIXT_Unscale_Br_Len_Multiplier_Tree(tree);
      return(-1.);
    }

  For(i,2*tree->n_otu-1) tree->a_edges[i]->l->v /= tree->mod->br_len_mult->v;
  return(-1.);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Reflect(phydbl x, phydbl l, phydbl u)
{
  int rounds;
  phydbl tmp;
  int k;

  if(u < l)
    {
      tmp = u;
      u   = l;
      l   = tmp;
    }

  if(x < l) x = x + 2.*(l - x);

  if(((x-u) > (u-l)) && (x > u))
    {
      k = (x - (2.*u-l))/(2.*(u-l));
      x = x - 2.*k*(u-l);
    }

  rounds = 0;
  do
    {
      rounds++;
      /* printf("\n. l=%f u=%f x=%f",l,u,x); */
      if(x > u || x < l)
    {
      if(x > u) x = x - 2.*(x - u);
      else      x = x + 2.*(l - x);
    }
      else break;
      /* printf(" x'=%f",x); */
    }
  while(rounds < 100);

  if(rounds == 100 && (x > u || x < l))
    {
      PhyML_Printf("\n. u=%f l=%f x=%f",u,l,x);
      PhyML_Printf("\n. Err. in file %s at line %d\n",__FILE__,__LINE__);
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
    }

  return x;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


int Are_Equal(phydbl a, phydbl b, phydbl eps)
{
  if(FABS(a-b) < eps) return TRUE; /* a==b */
  else return FALSE;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* Returns 1 if small_tree is displayed by big_tree, 0 otherwise
   Does not account for the root positions, if any.
*/
int Check_Topo_Constraints(t_tree *big_tree, t_tree *small_tree)
{
  if(!small_tree) return 1;

  if(small_tree->n_otu < 4) return 1;

  if(small_tree->n_otu > big_tree->n_otu)
    {
      PhyML_Printf("\n");
      PhyML_Printf("\n. The tree that defines the topological constraints can not");
      PhyML_Printf("\n. display more taxa than %d",big_tree->n_otu);
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
    }

  t_tree *big_tree_cpy;
  int diffs,i;

  big_tree_cpy = Make_Tree_From_Scratch(big_tree->n_otu,NULL);
  Copy_Tree(big_tree,big_tree_cpy);


  Prune_Tree(big_tree_cpy,small_tree);

  /* For(i,2*small_tree->n_otu-3) printf("\nz %d . %d . %d", */
  /* 				      big_tree->a_edges[i]->does_exist, */
  /* 				      big_tree_cpy->a_edges[i]->does_exist, */
  /* 				      small_tree->a_edges[i]->does_exist); */

  Free_Bip(small_tree);
  Alloc_Bip(small_tree);
  Get_Bip(small_tree->a_nodes[0],small_tree->a_nodes[0]->v[0],small_tree);

  Free_Bip(big_tree_cpy);
  Alloc_Bip(big_tree_cpy);
  Match_Tip_Numbers(small_tree,big_tree_cpy);
  Get_Bip(big_tree_cpy->a_nodes[0],big_tree_cpy->a_nodes[0]->v[0],big_tree_cpy);

  For(i,2*big_tree_cpy->n_otu-3) big_tree_cpy->a_edges[i]->bip_score = 0;
  For(i,2*small_tree->n_otu-3)   small_tree->a_edges[i]->bip_score = 0;

  diffs = Compare_Bip(small_tree,big_tree_cpy,NO);

  /* printf("\n"); */
  /* printf("\n. %s",Write_Tree(big_tree_cpy,NO)); */
  /* printf("\n. %s",Write_Tree(small_tree,NO)); */
  /* printf("\n. diffs=%d",diffs); */


  Free_Tree(big_tree_cpy);

  t_tree *big_tree_cpy_bis;
  big_tree_cpy_bis = Make_Tree_From_Scratch(big_tree->n_otu,NULL);
  Copy_Tree(big_tree,big_tree_cpy_bis);
  Free_Tree(big_tree_cpy_bis);

  if(diffs == 0) return 1; /* Constraint is satisfied */
  else           return 0;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Prune_Tree(t_tree *big_tree, t_tree *small_tree)
{
  int i,j;
  int curr_ext_node, curr_int_node, curr_br, n_pruned_nodes;;
  t_node **pruned_nodes;
  t_edge **residual_edges;


  pruned_nodes   = (t_node **)mCalloc(big_tree->n_otu,sizeof(t_node *));
  residual_edges = (t_edge **)mCalloc(big_tree->n_otu,sizeof(t_edge *));

  n_pruned_nodes = 0;
  for(i=0;i<big_tree->n_otu;i++)
    {
      for(j=0;j<small_tree->n_otu;j++)
        if(!strcmp(small_tree->a_nodes[j]->name,big_tree->a_nodes[i]->name))
          break;
      
      if(j == small_tree->n_otu)
        {
          Prune_Subtree(big_tree->a_nodes[i]->v[0],
                        big_tree->a_nodes[i],
                        NULL,&(residual_edges[n_pruned_nodes]),
                        big_tree);
          
          pruned_nodes[n_pruned_nodes] = big_tree->a_nodes[i];
          n_pruned_nodes++;
        }
    }
  
  if(!n_pruned_nodes)
    {
      Free(pruned_nodes);
      Free(residual_edges);
      return;
    }
  
  Free(big_tree->t_dir);
  
  big_tree->n_otu -= n_pruned_nodes;
  

  curr_ext_node = 0;
  curr_int_node = big_tree->n_otu;
  curr_br = 0;
  For(i,big_tree->n_otu+n_pruned_nodes)
    {
      for(j=0;j<n_pruned_nodes;j++)
        if(!strcmp(pruned_nodes[j]->name,big_tree->a_nodes[i]->name))
          break;
      
      if(j == n_pruned_nodes) /* That t_node still belongs to the tree */
        {
          Reassign_Node_Nums(big_tree->a_nodes[i],big_tree->a_nodes[i]->v[0],
                             &curr_ext_node,&curr_int_node,big_tree);
          break;
        }
    }
  

  Reassign_Edge_Nums(big_tree->a_nodes[0],big_tree->a_nodes[0]->v[0],&curr_br,big_tree);
  

  big_tree->t_dir = (short int *)mCalloc((2*big_tree->n_otu-2)*(2*big_tree->n_otu-2),sizeof(short int));
  
  for(i=0;i<n_pruned_nodes;i++)
    {
      Free_Edge(residual_edges[i]);
      Free_Edge(pruned_nodes[i]->b[0]);
      Free_Node(pruned_nodes[i]->v[0]);
      Free_Node(pruned_nodes[i]);
    }
  
  Free(pruned_nodes);
  Free(residual_edges);

  big_tree->a_edges[2*big_tree->n_otu-3] = big_tree->a_edges[2*(big_tree->n_otu+n_pruned_nodes)-3];
  big_tree->a_edges[2*big_tree->n_otu-2] = big_tree->a_edges[2*(big_tree->n_otu+n_pruned_nodes)-2];
  big_tree->a_nodes[2*big_tree->n_otu-2] = big_tree->a_nodes[2*(big_tree->n_otu+n_pruned_nodes)-2];

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* For every node in small_tree, find which node in big_tree
   it corresponds to and initialize the variable match_node
   accordingly (vice versa for big_tree)
*/
void Match_Nodes_In_Small_Tree(t_tree *small_tree, t_tree *big_tree)
{
  int i,j,k,l,m,n,identical;
  int *score;

  if(small_tree->n_otu > big_tree->n_otu)
    {
      PhyML_Printf("\n. small_tree->n_otu=%d big_tree->n_otu=%d",small_tree->n_otu,big_tree->n_otu);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
    }

  Free_Bip(big_tree);
  Alloc_Bip(big_tree);
  Get_Bip(big_tree->a_nodes[0],big_tree->a_nodes[0]->v[0],big_tree);

  Free_Bip(small_tree);
  Alloc_Bip(small_tree);
  Get_Bip(small_tree->a_nodes[0],small_tree->a_nodes[0]->v[0],small_tree);

  if(!Check_Topo_Constraints(big_tree,small_tree))
    {
      PhyML_Printf("\n. small_tree and big_tree cannot have distinct topologies.");
      PhyML_Printf("\n. Err. in file %s at line %d\n",__FILE__,__LINE__);
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
    }

  For(i,2*small_tree->n_otu-1) small_tree->a_nodes[i]->match_node = NULL;
  For(i,2*big_tree->n_otu-1)   big_tree->a_nodes[i]->match_node   = NULL;

  score = (int *)mCalloc(3,sizeof(int));

  for(i=0;i<small_tree->n_otu;i++)
    {
      for(j=0;j<big_tree->n_otu;j++)
    {
      if(!strcmp(small_tree->a_nodes[i]->name,big_tree->a_nodes[j]->name))
        {
          small_tree->a_nodes[i]->match_node = big_tree->a_nodes[j];
          big_tree->a_nodes[j]->match_node   = small_tree->a_nodes[i];
          break;
        }
    }
    }

  For(i,2*small_tree->n_otu-2)
    {
      if(small_tree->a_nodes[i]->tax == NO)
    {
      For(j,2*big_tree->n_otu-2)
        {
          if(big_tree->a_nodes[j]->tax == NO)
        {
          for(k=0;k<3;k++) score[k] = 0;

          for(k=0;k<3;k++)
            {
              for(l=0;l<3;l++)
            {
              identical = 0;
              For(m,small_tree->a_nodes[i]->bip_size[k])
                {
                  For(n,big_tree->a_nodes[j]->bip_size[l])
                {
                  if(!strcmp(small_tree->a_nodes[i]->bip_node[k][m]->name,big_tree->a_nodes[j]->bip_node[l][n]->name))
                    {
                      identical++;
                      break;
                    }
                }
                }
              if(identical == small_tree->a_nodes[i]->bip_size[k])
                {
                  score[k]++;
                }
            }
            }

          /* printf("\n. [%d] [%d] %d %d %d -- %d %d %d",i,j, */
          /* 	 score[0],score[1],score[2], */
          /* 	 small_tree->a_nodes[i]->bip_size[0], */
          /* 	 small_tree->a_nodes[i]->bip_size[1], */
          /* 	 small_tree->a_nodes[i]->bip_size[2]); */

          if(
             score[0] == 1 &&
             score[1] == 1 &&
             score[2] == 1
             )
            {
              small_tree->a_nodes[i]->match_node = big_tree->a_nodes[j];
              big_tree->a_nodes[j]->match_node   = small_tree->a_nodes[i];
              break;
            }
        }
        }
    }
    }

  Free(score);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Find_Surviving_Edges_In_Small_Tree(t_tree *small_tree, t_tree *big_tree)
{
  int i;

  Match_Nodes_In_Small_Tree(small_tree,big_tree);

  For(i,2*small_tree->n_otu-1) small_tree->rates->has_survived[i] = NO;

  Find_Surviving_Edges_In_Small_Tree_Post(big_tree->n_root,big_tree->n_root->v[2],small_tree,big_tree);
  Find_Surviving_Edges_In_Small_Tree_Post(big_tree->n_root,big_tree->n_root->v[1],small_tree,big_tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Find_Surviving_Edges_In_Small_Tree_Post(t_node *a, t_node *d, t_tree *small_tree, t_tree *big_tree)
{
  if(d->match_node && !a->match_node)
    {
      small_tree->rates->has_survived[d->match_node->num] = YES;
    }

  if(d->tax == YES) return;
  else
    {
      int i;

      for(i=0;i<3;i++)
    {
      if(d->v[i] != a && d->b[i] != big_tree->e_root)
        {
          Find_Surviving_Edges_In_Small_Tree_Post(d,d->v[i],small_tree,big_tree);
        }
    }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Set_Taxa_Id_Ranking(t_tree *tree)
{
  int i,j;

  for(i=0;i<tree->n_otu;i++) tree->a_nodes[i]->id_rank = 0;

  for(i=0;i<tree->n_otu;i++)
    {
      for(j=i+1;j<tree->n_otu;j++)
    {
      if(strcmp(tree->a_nodes[i]->name,tree->a_nodes[j]->name) > 0)
        tree->a_nodes[i]->id_rank++;
      else
        tree->a_nodes[j]->id_rank++;
    }
    }
  /* for(i=0;i<tree->n_otu;i++) PhyML_Printf("\n. %20s %4d",tree->a_nodes[i]->name,tree->a_nodes[i]->id_rank); */
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Get_Edge_Binary_Coding_Number(t_tree *tree)
{
  int i,j;
  int list_size;
  t_node **list;
  t_edge *b;
  int max_left,max_rght;

  if(tree->n_otu > 1000)
    {
      PhyML_Printf("\n. Can't work out edge binary code if the number of taxa >1000.");
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
    }

  Free_Bip(tree);
  Alloc_Bip(tree);
  Get_Bip(tree->a_nodes[0],tree->a_nodes[0]->v[0],tree);

  Set_Taxa_Id_Ranking(tree);

  b = NULL;
  For(i,2*tree->n_otu-3)
    {
      b = tree->a_edges[i];

      max_left = 0;
      For(j,b->left->bip_size[b->l_r])
    if(b->left->bip_node[b->l_r][j]->id_rank > max_left)
      max_left = b->left->bip_node[b->l_r][j]->id_rank;

      max_rght = 0;
      For(j,b->rght->bip_size[b->r_l])
    if(b->rght->bip_node[b->r_l][j]->id_rank > max_rght)
      max_rght = b->rght->bip_node[b->r_l][j]->id_rank;


      if(max_left < max_rght)
    {
      list = b->left->bip_node[b->l_r];
      list_size = b->left->bip_size[b->l_r];
    }
      else
    {
      list = b->rght->bip_node[b->r_l];
      list_size = b->rght->bip_size[b->r_l];
    }

      b->bin_cod_num = 0.;
      for(j=0;j<list_size;j++) b->bin_cod_num += POW(2,list[j]->id_rank);
      /* printf("\n. %f",b->bin_cod_num); */
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Make_MutMap(t_tree *tree)
{
  // (# of edges) X (# of sites) X (# mutation types: A<->C A<->G A<->T C<->G C<->T G<->T)
  tree->mutmap = (int *)mCalloc((2*tree->n_otu-3)*(tree->n_pattern)*6,sizeof(int));
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int Get_Mutmap_Val(int edge, int site, int mut, t_tree *tree)
{
  int dim1,dim2;

  dim1 = (tree->n_pattern)*(2*tree->n_otu-3);
  dim2 = (tree->n_pattern);

  return tree->mutmap[mut*dim1 + edge*dim2 + site];
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Get_Mutmap_Coord(int idx, int *edge, int *site, int *mut, t_tree *tree)
{
  int dim1,dim2;

  dim1 = (tree->n_pattern)*(2*tree->n_otu-3);
  dim2 = (tree->n_pattern);

  (*mut) =  (int)idx/dim1;
  (*edge) = (int)(idx - (*mut)*dim1)/dim2;
  (*site) = (int)(idx - (*mut)*dim1 - (*edge)*dim2);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Copy_Edge_Lengths(t_tree *to, t_tree *from)
{
  int i;
  For(i,2*from->n_otu-1) to->a_edges[i]->l->v = from->a_edges[i]->l->v;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

char *To_Lower_String(char *in)
{
  char *out;
  int i;
  int len;

  len = (int)strlen(in);

  out = (char *)mCalloc(len+1,sizeof(char));

  for(i=0;i<len;i++) out[i] = (char)tolower(in[i]);

  out[len] = '\0';
  return(out);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

char *To_Upper_String(char *in)
{
  char *out;
  int i;
  int len;

  len = (int)strlen(in);

  out = (char *)mCalloc(len+1,sizeof(char));

  for(i=0;i<len;i++)
    {
      out[i] = (char)toupper(in[i]);
    }

  out[len] = '\0';
  return(out);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Connect_CSeqs_To_Nodes(calign *cdata, option *io, t_tree *tree)
{
  int i,j,n_otu_tree,n_otu_cdata;
  
  n_otu_tree  = tree->n_otu;
  n_otu_cdata = cdata->n_otu;
  
  if((n_otu_tree != n_otu_cdata) && (io->fp_in_constraint_tree == NULL))
    {
      PhyML_Printf("\n. Number of taxa in the tree: %d, number of sequences: %d.",n_otu_tree,n_otu_cdata);
      Warn_And_Exit("\n. The number of tips in the tree is not the same as the number of sequences\n");
    }
  
  for(i=0;i<n_otu_tree;i++)
    {
      for(j=0;j<n_otu_cdata;j++)
        {          
          if(!strcmp(tree->a_nodes[i]->name,cdata->c_seq[j]->name)) break;
        }
      
      if(j==n_otu_cdata)
        {
          PhyML_Printf("\n. Taxon '%s' was not found in sequence file '%s'.\n",
                       tree->a_nodes[i]->name,
                       io->in_align_file);
          Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
        }
      tree->a_nodes[i]->c_seq = cdata->c_seq[j];
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Switch_Eigen(int state, t_mod *mod)
{
  t_mod *buff;

  buff = mod;
  do
    {
      buff->update_eigen = state;
      buff = buff->next;
      if(!buff) break;
    }
  while(1);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Set_Both_Sides(int yesno, t_tree *tree)
{
  tree->both_sides = yesno;
  if(tree->is_mixt_tree == YES) MIXT_Set_Both_Sides(yesno,tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Set_Use_Eigen_Lr(int yesno, t_tree *tree)
{
  tree->use_eigen_lr = yesno;
  if(tree->is_mixt_tree == YES) MIXT_Set_Use_Eigen_Lr(yesno,tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Set_Update_Eigen_Lr(int yesno, t_tree *tree)
{
  tree->update_eigen_lr = yesno;
  if(tree->is_mixt_tree == YES) MIXT_Set_Update_Eigen_Lr(yesno,tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// Returns the matrix of pairwise distances between tips
phydbl *Dist_Btw_Tips(t_tree *tree)
{
  int i,j;
  phydbl *dist;

  dist = (phydbl *)mCalloc(tree->n_otu*tree->n_otu,sizeof(phydbl));

  for(i=0;i<tree->n_otu-1;i++)
    {
      for(j=i+1;j<tree->n_otu;j++)
        {
          Path_Length(tree->a_nodes[i],tree->a_nodes[j],dist+i*tree->n_otu+j,tree);
          dist[j*tree->n_otu+i] = dist[i*tree->n_otu+j];
        }
    }


  return(dist);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Random_SPRs_On_Rooted_Tree(t_tree *tree)
{
  int i,j;
  t_node *a, *d, *v1, *v2,*target_nd;
  phydbl u;
  phydbl lim_sup,lim_inf;
  t_edge *residual,*nouse;
  phydbl new_t;
  int err;
  int ok_dir,dir12,dir21;

  /* printf("\n. >>>>>>>>>>>>>>.\n"); */
  /* Print_Node(tree->n_root,tree->n_root->v[1],tree); */
  /* Print_Node(tree->n_root,tree->n_root->v[2],tree); */
  /* printf("\n. >>>>>>>>>>>>>>."); */
  /* GEO_Update_Occup(tree->geo,tree); */
  /* GEO_Lk(tree->geo,tree); */
  /* PhyML_Printf("\n>> Init loglk: %f",tree->geo->c_lnL); */


  target_nd = NULL;
  residual = NULL;
  nouse = NULL;
  a = d = v1 = v2 = NULL;
  for(i=0;i<tree->n_otu;i++) // We will perform tree->n_otu random SPRs
    {
      a = tree->a_nodes[Rand_Int(tree->n_otu,2*tree->n_otu-3)];

      if(a == tree->n_root)
        {
          PhyML_Printf("\n. Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
          Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
        }
      /* if(a == tree->n_root->v[1] || a == tree->n_root->v[2]) continue; */

      v1 = v2 = NULL;
      for(j=0;j<3;j++)
        {
          if(a->v[j] != a->anc && a->b[j] != tree->e_root)
            {
              if(!v1) v1 = a->v[j];
              else    v2 = a->v[j];
            }
        }

      u = Uni();
      if(u < .5) d = v1;
      else       d = v2;


      for(j=0;j<3;j++) if(a->v[j] == d && a->b[j] == tree->e_root) break;
      if(j != 3) continue;

      lim_sup = tree->rates->nd_t[d->num];
      lim_inf = tree->rates->nd_t[tree->n_root->num];

      phydbl scale = Uni()*(50.-5.)+5.;
      scale = 50;
      err = NO;
      new_t = Rnorm_Trunc(tree->rates->nd_t[a->num],
                          FABS(tree->rates->nd_t[tree->n_root->num]/scale),
                          lim_inf,
                          lim_sup,
                          &err);


      if(err == YES) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);

      if(new_t < tree->rates->nd_t[a->num])
        {
          target_nd = a;
          do
            {
              if(tree->rates->nd_t[target_nd->anc->num] < new_t)
                {
                  if(target_nd == a)
                    {
                      if(v1 == d) target_nd = v2;
                      else        target_nd = v1;
                    }
                  break;
                }
              else
                {
                  target_nd = target_nd->anc;
                }
            }
          while(1);
        }
      else
        {
          if(v1 == d) target_nd = v2;
          else        target_nd = v1;

          do
            {
              if(tree->rates->nd_t[target_nd->num] > new_t)
                {
                  break;
                }
              else
                {
                  v1 = v2 = NULL;
                  for(j=0;j<3;j++)
                    {
                      if(target_nd->v[j] != target_nd->anc)
                        {
                          if(!v1) v1 = target_nd->v[j];
                          else    v2 = target_nd->v[j];
                        }
                    }
                  u = Uni();
                  if(u < 0.5) target_nd = v1;
                  else        target_nd = v2;
                }
            }while(1);
        }

      ok_dir = -1;
      for(j=0;j<3;j++) if(target_nd->v[j] == target_nd->anc || target_nd->b[j] == tree->e_root) ok_dir = j;

      dir12 = -1;
      for(j=0;j<3;j++) if(tree->n_root->v[1]->v[j] == tree->n_root->v[2]) dir12 = j;

      if(dir12 == -1) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);        

      dir21 = -1;
      for(j=0;j<3;j++) if(tree->n_root->v[2]->v[j] == tree->n_root->v[1]) dir21 = j;

      if(dir21 == -1) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);



      /* a = tree->a_nodes[101]; */
      /* d = tree->a_nodes[105]; */
      /* target_nd = tree->a_nodes[105]; */
      /* printf("\n pull %d %d target %d",a->num,d->num,target_nd->num); */
      /* printf("\n. %d %d %d",tree->n_root->num,tree->n_root->v[1]->num,tree->n_root->v[2]->num); */
      /* printf("\n. %d %d",tree->e_root->num,target_nd->b[ok_dir]->num); */

      Prune_Subtree(a,d,&nouse,&residual,tree);

      Graft_Subtree(target_nd->b[ok_dir],a,NULL,residual,NULL,tree);

      tree->rates->nd_t[a->num] = new_t;

      if(target_nd == tree->n_root->v[1])
        {
          tree->n_root->v[1] = a;
          tree->e_root = tree->n_root->v[2]->b[dir21];
        }

      if(target_nd == tree->n_root->v[2])
        {
          tree->n_root->v[2] = a;
          tree->e_root = tree->n_root->v[1]->b[dir12];
        }

      if(tree->n_root->v[1] == NULL)
        {
          tree->n_root->v[1] = tree->n_root->v[2]->v[dir21];
          tree->e_root = tree->n_root->v[2]->b[dir21];
        }

      if(tree->n_root->v[2] == NULL)
        {
          tree->n_root->v[2] = tree->n_root->v[1]->v[dir12];
          tree->e_root = tree->n_root->v[1]->b[dir12];
        }

      Update_Ancestors(tree->n_root,tree->n_root->v[2],tree);
      Update_Ancestors(tree->n_root,tree->n_root->v[1],tree);
      tree->n_root->anc = NULL;

      /* Print_Node(tree->n_root,tree->n_root->v[1],tree); */
      /* Print_Node(tree->n_root,tree->n_root->v[2],tree); */

      Time_To_Branch(tree);

      /* GEO_Update_Occup(tree->geo,tree); */
      /* GEO_Lk(tree->geo,tree); */
      /* PhyML_Printf("\nxx Init loglk: %f",tree->geo->c_lnL); */

    }


}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Best_Root_Position_IL_Model(t_tree *tree)
{

  if(tree->n_root)
    {
      PhyML_Printf("\n. The tree already has a root node");
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
    }
  else
    {
      int i;
      t_edge *best_edge;
      phydbl best_lnL;

      Free_Edge_Lk_Rght(tree->a_edges[2*tree->n_otu-3]);
      Free_Edge_Lk_Rght(tree->a_edges[2*tree->n_otu-2]);
      Free_Edge_Pars_Rght(tree->a_edges[2*tree->n_otu-3]);
      Free_Edge_Pars_Rght(tree->a_edges[2*tree->n_otu-2]);


      best_edge = NULL;
      best_lnL  = UNLIKELY;
      For(i,2*tree->n_otu-3)
        {
          PhyML_Printf("\n. Positionning root node on edge %4d",tree->a_edges[i]->num);
          Add_Root(tree->a_edges[i],tree);
          tree->ignore_root = NO;
          Set_Both_Sides(YES,tree);
          Lk(NULL,tree);

          /* Optimize_Br_Len_Serie(tree); */

          Update_Partial_Lk(tree,tree->n_root->b[1],tree->n_root);
          Br_Len_Brent(tree->n_root->b[1],tree);
          Update_Partial_Lk(tree,tree->n_root->b[2],tree->n_root);
          Br_Len_Brent(tree->n_root->b[2],tree);

          PhyML_Printf(" -- lnL: %20f",tree->c_lnL);
          if(tree->c_lnL > best_lnL)
            {
              best_lnL  = tree->c_lnL;
              best_edge = tree->a_edges[i];
            }
        }

      Add_Root(best_edge,tree);
      Set_Both_Sides(YES,tree);
      Lk(NULL,tree);
      Update_Partial_Lk(tree,tree->n_root->b[1],tree->n_root);
      Br_Len_Brent(tree->n_root->b[1],tree);
      Update_Partial_Lk(tree,tree->n_root->b[2],tree->n_root);
      Br_Len_Brent(tree->n_root->b[2],tree);
      tree->ignore_root = YES;
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Set_Br_Len_Var(t_edge *b, t_tree *tree)
{
  if(tree->is_mixt_tree)
    {
      MIXT_Set_Br_Len_Var(b,tree);
      return;
    }

  if(tree->rates == NO && tree->mod->gamma_mgf_bl == YES)
    {

      phydbl len;
      if(b == NULL)
        {
          int i;
          
          For(i,2*tree->n_otu-1)
            {
              len = MAX(0.0,tree->a_edges[i]->l->v);
              tree->a_edges[i]->l_var->v = POW(len,2)*tree->mod->l_var_sigma;
            }
        }
      else
        {
          len = MAX(0.0,b->l->v);
          b->l_var->v = POW(len,2)*tree->mod->l_var_sigma;
        }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Check_Br_Lens(t_tree *tree)
{
  int i;
  scalar_dbl *l;

  For(i,2*tree->n_otu-1)
    {
      l = tree->a_edges[i]->l;
      do
        {
          /* if(l->v < tree->mod->l_min) l->v = tree->mod->l_min; */
          /* if(l->v > tree->mod->l_max) l->v = tree->mod->l_max; */
          if(l->v < 0.0) l->v = 0.0;
          l = l->next;
        }
      while(l);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Build_Distrib_Number_Of_Diff_States_Under_Model(t_tree *tree)
{
  calign *orig_data;
  t_mod *orig_mod;
  int iter,n_iter_tot,i,j;
  phydbl *n_diff_states_all_l,*n_diff_states_all_r;

  Calculate_Number_Of_Diff_States(tree);
  
  PhyML_Printf("\n TRUE     edge    side    states val");
  For(i,2*tree->n_otu-3)
    {
      if(tree->a_edges[i]->left->tax == NO && tree->a_edges[i]->rght->tax == NO)
        {
          for(j=0;j<tree->mod->ns;j++)
            {
              PhyML_Printf("\n TRUE %3d 0 %3d %d",
                           i,
                           j+1,
                           tree->a_edges[i]->n_diff_states_l[j]);
              
              PhyML_Printf("\n TRUE %3d 1 %3d %d",
                           i,
                           j+1,
                           tree->a_edges[i]->n_diff_states_r[j]);
            }
        }
    }



  n_iter_tot = 100;

  n_diff_states_all_l = (phydbl *)mCalloc((n_iter_tot) * (tree->mod->ns) * (2*tree->n_otu-3) * 2, sizeof(phydbl));
  n_diff_states_all_r = (phydbl *)mCalloc((n_iter_tot) * (tree->mod->ns) * (2*tree->n_otu-3) * 2, sizeof(phydbl));

  orig_mod  = Copy_Model(tree->mod);
  orig_data = Copy_Cseq(tree->data,tree->io);

  orig_mod->io = tree->io;
  orig_mod->s_opt = tree->mod->s_opt;

  iter = 0;

  do
    {
      Evolve(tree->data,tree->mod,0,tree);

      Calculate_Number_Of_Diff_States(tree);

      For(i,2*tree->n_otu-3)
        {
          for(j=0;j<tree->mod->ns;j++)
            {
              n_diff_states_all_l[j*(2*tree->n_otu-3)*(n_iter_tot) + i*(n_iter_tot) + iter] = tree->a_edges[i]->n_diff_states_l[j];
              n_diff_states_all_r[j*(2*tree->n_otu-3)*(n_iter_tot) + i*(n_iter_tot) + iter] = tree->a_edges[i]->n_diff_states_r[j];
            }
        }


      Free_Calign(tree->data);
      Free_Model_Complete(tree->mod);
      Free_Model_Basic(tree->mod);
      
      tree->mod  = Copy_Model(orig_mod);
      tree->data = Copy_Cseq(orig_data,tree->io);

      tree->mod->io    = orig_mod->io;
      tree->mod->s_opt = orig_mod->s_opt;

      Connect_CSeqs_To_Nodes(tree->data,tree->io,tree);

      iter++;
    }
  while(iter < n_iter_tot);


  PhyML_Printf("\n SIM     edge    side    states low      up");
  For(i,2*tree->n_otu-3)
    {
      if(tree->a_edges[i]->left->tax == NO && tree->a_edges[i]->rght->tax == NO)
        {
          for(j=0;j<tree->mod->ns;j++)
            {              
              PhyML_Printf("\n SIM %3d 0 %3d %.0f %.0f",
                     i,
                     j+1,
                     Quantile(n_diff_states_all_l + j*(2*tree->n_otu-3)*(n_iter_tot) + i*(n_iter_tot), n_iter_tot, 0.10),
                     Quantile(n_diff_states_all_l + j*(2*tree->n_otu-3)*(n_iter_tot) + i*(n_iter_tot), n_iter_tot, 0.90));
              
              PhyML_Printf("\n SIM %3d 1 %3d %.0f %.0f",
                           i,
                           j+1,
                           Quantile(n_diff_states_all_r + j*(2*tree->n_otu-3)*(n_iter_tot) + i*(n_iter_tot), n_iter_tot, 0.10),
                           Quantile(n_diff_states_all_r + j*(2*tree->n_otu-3)*(n_iter_tot) + i*(n_iter_tot), n_iter_tot, 0.90));
            }
        }
    }


  Add_Root(tree->a_edges[0],tree);
  DR_Draw_Tree("treefile",tree);

   
  Free(n_diff_states_all_l);
  Free(n_diff_states_all_r);
  
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* Calculate the number of sites at which 1,...,n states (n: 4 or 20) */
/* are observed, for every subtree */
 
void Calculate_Number_Of_Diff_States(t_tree *tree)
{  
  Init_Ui_Tips(tree);
  Calculate_Number_Of_Diff_States_Post(tree->a_nodes[0],tree->a_nodes[0]->v[0],tree->a_nodes[0]->b[0],tree);
  Calculate_Number_Of_Diff_States_Pre(tree->a_nodes[0],tree->a_nodes[0]->v[0],tree->a_nodes[0]->b[0],tree);

  /* int i; */

  /* For(i,2*tree->n_otu-3) */
  /*   { */
  /*     if(tree->a_edges[i]->left->tax == NO && tree->a_edges[i]->rght->tax == NO) */
        /* printf("\n. Edge %d left : %d %d %d %d right: %d %d %d %d", */
        /*        i, */
        /*        tree->a_edges[i]->n_diff_states_l[0], */
        /*        tree->a_edges[i]->n_diff_states_l[1], */
        /*        tree->a_edges[i]->n_diff_states_l[2], */
        /*        tree->a_edges[i]->n_diff_states_l[3], */
        /*        tree->a_edges[i]->n_diff_states_r[0], */
        /*        tree->a_edges[i]->n_diff_states_r[1], */
        /*        tree->a_edges[i]->n_diff_states_r[2], */
        /*        tree->a_edges[i]->n_diff_states_r[3]); */
    /* } */
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Calculate_Number_Of_Diff_States_Post(t_node *a, t_node *d, t_edge *b, t_tree *tree)
{
  if(d->tax) return;
  else
    {
      int i;

      for(i=0;i<3;i++)
        if(d->v[i] != a)
          Calculate_Number_Of_Diff_States_Post(d,d->v[i],d->b[i],tree);

      Calculate_Number_Of_Diff_States_Core(a,d,b,tree);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Calculate_Number_Of_Diff_States_Pre(t_node *a, t_node *d, t_edge *b, t_tree *tree)
{

  if(d->tax) return;
  else
    {
      int i;

      for(i=0;i<3;i++)
        if(d->v[i] != a)
          {
            Calculate_Number_Of_Diff_States_Core(d->v[i],d,d->b[i],tree);
            Calculate_Number_Of_Diff_States_Pre(d,d->v[i],d->b[i],tree);
          }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Calculate_Number_Of_Diff_States_Core(t_node *a, t_node *d, t_edge *b, t_tree *tree)
{
  int *ui, *ui_v1, *ui_v2;
  int sum,site,state;
  int *diff;
  t_node *v1, *v2;

  ui = ui_v1 = ui_v2 = NULL;
  v1 = v2 = NULL;

  if(d == b->left)
    {
      v1 = (d == d->b[b->l_v1]->left)?
        (d->b[b->l_v1]->rght):
        (d->b[b->l_v1]->left);

      v2 = (d == d->b[b->l_v2]->left)?
        (d->b[b->l_v2]->rght):
        (d->b[b->l_v2]->left);
       
      ui = b->ui_l;
      diff = b->n_diff_states_l;

      ui_v1 = 
      (d == d->b[b->l_v1]->left)?
      (d->b[b->l_v1]->ui_r):
      (d->b[b->l_v1]->ui_l);

      ui_v2 = 
      (d == d->b[b->l_v2]->left)?
      (d->b[b->l_v2]->ui_r):
      (d->b[b->l_v2]->ui_l);

    }
  else
    {
      v1 = (d == d->b[b->r_v1]->left)?
        (d->b[b->r_v1]->rght):
        (d->b[b->r_v1]->left);

      v2 = (d == d->b[b->r_v2]->left)?
        (d->b[b->r_v2]->rght):
        (d->b[b->r_v2]->left);

      ui = b->ui_r;
      diff = b->n_diff_states_r;
      
      ui_v1 = 
        (d == d->b[b->r_v1]->left)?
        (d->b[b->r_v1]->ui_r):
        (d->b[b->r_v1]->ui_l);
      
      ui_v2 = 
      (d == d->b[b->r_v2]->left)?
      (d->b[b->r_v2]->ui_r):
      (d->b[b->r_v2]->ui_l);

    }

  for(state=0;state<tree->mod->ns;state++) diff[state] = 0;
  
  for(site=0;site<tree->n_pattern;site++)
    {
      if(v1->tax == YES)
        {
          int sum;
          sum = Sum_Bits(ui_v1[site],tree->mod->ns);
          if(sum > 1)
            {
              int val = ui_v1[site];
              int pos, iter;
              phydbl u = Uni();

              iter = 0;
              do
                {
                  pos = Rand_Int(0,tree->mod->ns-1);
                  if(((val >> pos) & 1) && (u > 1./sum)) break;
                }
              while(iter++ < 1000);


              if(iter == 1000) 
                {
                  Generic_Exit(__FILE__,__LINE__,__FUNCTION__);                
                }

              ui_v1[site] = POW(2,pos);
            }
        }

      if(v2->tax == YES)
        {
          int sum;
          sum = Sum_Bits(ui_v2[site],tree->mod->ns);
          if(sum > 1)
            {
              int val = ui_v2[site];
              int pos, iter;
              phydbl u = Uni();
              
              iter = 0;
              do
                {
                  pos = Rand_Int(0,tree->mod->ns-1);
                  if(((val >> pos) & 1) && (u > 1./sum)) break;
                }
              while(iter++ < 1000);

              if(iter == 1000) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);

              ui_v2[site] = POW(2,pos);
            }
        }

      ui[site] = ui_v1[site] | ui_v2[site];
      
      sum = Sum_Bits(ui[site],tree->mod->ns);

      /* printf("\n. ui_v1: %d ui_v2: %d ui: %d sum: %d",ui_v1[site],ui_v2[site],ui[site],sum); fflush(NULL); */

      if(sum-1 > tree->mod->ns-1) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);

      diff[sum-1]++;

    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
/* Returns the number of distinct states observed at a particular site */

int Number_Of_Diff_States_One_Site(int site, t_tree *tree)
{
  int n_states;


  Number_Of_Diff_States_One_Site_Post(tree->a_nodes[0],
                                      tree->a_nodes[0]->v[0],
                                      tree->a_nodes[0]->b[0],
                                      site,tree);
  
  n_states = Sum_Bits(tree->a_nodes[0]->b[0]->ui_r[site] | tree->a_nodes[0]->b[0]->ui_l[site],tree->mod->ns);


  return(n_states);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Number_Of_Diff_States_One_Site_Post(t_node *a, t_node *d, t_edge *b, int site, t_tree *tree)
{
  if(d->tax) return;
  else
    {
      int i;

      for(i=0;i<3;i++)
        if(d->v[i] != a && d->b[i] != tree->e_root)
          Number_Of_Diff_States_One_Site_Post(d,d->v[i],d->b[i],site,tree);

      Number_Of_Diff_States_One_Site_Core(a,d,b,site,tree);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int Number_Of_Diff_States_One_Site_Core(t_node *a, t_node *d, t_edge *b, int site, t_tree *tree)
{
  int *ui, *ui_v1, *ui_v2;
  int sum;
  t_node *v1, *v2;
  
  ui = ui_v1 = ui_v2 = NULL;
  v1 = v2 = NULL;
  
  if(d == b->left)
    {
      v1 = (d == d->b[b->l_v1]->left)?
        (d->b[b->l_v1]->rght):
        (d->b[b->l_v1]->left);
      
      v2 = (d == d->b[b->l_v2]->left)?
        (d->b[b->l_v2]->rght):
        (d->b[b->l_v2]->left);
      
      ui = b->ui_l;
      
      ui_v1 = 
        (d == d->b[b->l_v1]->left)?
        (d->b[b->l_v1]->ui_r):
        (d->b[b->l_v1]->ui_l);
      
      ui_v2 = 
        (d == d->b[b->l_v2]->left)?
        (d->b[b->l_v2]->ui_r):
        (d->b[b->l_v2]->ui_l);      
    }
  else
    {
      v1 = (d == d->b[b->r_v1]->left)?
        (d->b[b->r_v1]->rght):
        (d->b[b->r_v1]->left);
      
      v2 = (d == d->b[b->r_v2]->left)?
        (d->b[b->r_v2]->rght):
        (d->b[b->r_v2]->left);
      
      ui = b->ui_r;
      
      ui_v1 = 
        (d == d->b[b->r_v1]->left)?
        (d->b[b->r_v1]->ui_r):
        (d->b[b->r_v1]->ui_l);
      
      ui_v2 = 
        (d == d->b[b->r_v2]->left)?
        (d->b[b->r_v2]->ui_r):
        (d->b[b->r_v2]->ui_l);      
    }
  
  if(v1->tax == YES) // Check for ambiguous character state at this tip
    {
      sum = Sum_Bits(ui_v1[site],tree->mod->ns);
      
      if(sum > 1)
        {
          int val = ui_v1[site];
          int pos, iter;
          phydbl u = Uni();

          // Select a state uniformly at random
          iter = 0;
          do
            {
              pos = Rand_Int(0,tree->mod->ns-1);
              if(((val >> pos) & 1) && (u > 1./sum)) break;
            }
          while(iter++ < 1000);


          if(iter == 1000) 
            {
              Generic_Exit(__FILE__,__LINE__,__FUNCTION__);            
            }
          
          ui_v1[site] = POW(2,pos);
        }
    }

  if(v2->tax == YES)
    {
      sum = Sum_Bits(ui_v2[site],tree->mod->ns);
      if(sum > 1)
        {
          int val = ui_v2[site];
          int pos, iter;
          phydbl u = Uni();
          
          iter = 0;
          do
            {
              pos = Rand_Int(0,tree->mod->ns-1);
              if(((val >> pos) & 1) && (u > 1./sum)) break;
            }
          while(iter++ < 1000);
          
          if(iter == 1000) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
          
          ui_v2[site] = POW(2,pos);
        }
    }
  
  ui[site] = ui_v1[site] | ui_v2[site];
      
  sum = Sum_Bits(ui[site],tree->mod->ns);

  /* printf("\n. ui_v1: %d ui_v2: %d ui: %d sum: %d",ui_v1[site],ui_v2[site],ui[site],sum); fflush(NULL); */
  
  if(sum-1 > tree->mod->ns-1)
    {
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
    }

  return sum;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl Get_Lk(t_tree *tree)
{
  t_tree *loc_tree;

  loc_tree = tree;
  /*! Rewind back to the first mixt_tree */
  while(loc_tree->prev) loc_tree = loc_tree->prev;

  return loc_tree->c_lnL;

}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl Get_dLk(t_tree *tree)
{
  t_tree *loc_tree;

  loc_tree = tree;
  /*! Rewind back to the first mixt_tree */
  while(loc_tree->prev) loc_tree = loc_tree->prev;

  return loc_tree->c_dlnL;

}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl Get_d2Lk(t_tree *tree)
{
  t_tree *loc_tree;

  loc_tree = tree;
  /*! Rewind back to the first mixt_tree */
  while(loc_tree->prev) loc_tree = loc_tree->prev;

  return loc_tree->c_d2lnL;

}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

align **Make_Empty_Alignment(option *io)
{
  int i;
  char *line;
  align **data;

  line   = (char *)mCalloc(T_MAX_LINE,sizeof(char));
  data   = (align **)mCalloc(io->n_otu,sizeof(align *));
  
  for(i=0;i<io->n_otu;i++)
    {
      data[i]        = (align *)mCalloc(1,sizeof(align));
      data[i]->name  = (char *)mCalloc(T_MAX_NAME,sizeof(char));
      data[i]->state = (char *)mCalloc(io->init_len*io->state_len+1,sizeof(char));

      data[i]->is_ambigu = NULL;
      data[i]->len = 0;
      
      Random_String(data[i]->name,5);

      while(data[i]->len < io->init_len * io->state_len)
        {
          data[i]->state[data[i]->len] = 'X';
          data[i]->len++;
        }
    }

  for(i=0;i<io->n_otu;i++) data[i]->state[data[i]->len] = '\0';

  Free(line);

  return data;
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
/* Mean observed frequency of difference between the n(n-1)/2 pairs of sequences */
phydbl Mean_Identity(calign *data)
{
  int i,j,n;
  phydbl tot_idt;
  
  n = data->n_otu;
  
  tot_idt = 0.0;
  for(i=0;i<n-1;i++)
    {
      for(j=i+1; j<n; j++)
        {
          tot_idt += Pairwise_Identity(i,j,data);
        }
    }

  return(tot_idt / (phydbl)(n*(n-1.)/2.));
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

/* Observed frequency of difference for the (i,j)-pair of sequences */
phydbl Pairwise_Identity(int i, int j, calign *data)
{
  int k;
  phydbl div,p,d;
  
  div = 0.0;
  for(k=0;k<data->crunch_len;k++) if(data->c_seq[i]->state[k] == data->c_seq[j]->state[k]) div += (phydbl)data->wght[k];

  /* observed proportion of identity */
  p = 1. - div / (phydbl)data->init_len;

  d = 0.0;
  if(data->io->datatype == NT)
    {
      if(p > 3./4.) return 0.25;
      else
        {
          /* Jukes & Cantor distance */
          d = -(3./4.)*log(1. - 4./3.*p);
        }
    }
  else if(data->io->datatype == AA)
    {
      if(p > 19./20.) return 1./20.;
      else
        {
          /* Jukes & Cantor distance */
          d = -(19./20.)*log(1. - 20./19.*p);
        }
    }
  else Generic_Exit(__FILE__,__LINE__,__FUNCTION__);

  return(exp(-d));
}


/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl Fst(int i, int j, calign *data)
{
  phydbl FA, Fr;

  FA = Mean_Identity(data);
  Fr = Pairwise_Identity(i,j,data);

  return((Fr-FA)/(1-FA));
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl Nucleotide_Diversity(calign *data)
{

  int i,j,n;
  phydbl pair_div;
  
  n = data->n_otu;
  
  pair_div = 0.0;
  for(i=0;i<n-1;i++)
    {
      for(j=i+1; j<n; j++)
        {
          pair_div += 1. - Pairwise_Identity(i,j,data);
        }
    }

  return(pair_div / (phydbl)(n*(n-1.)/2.));
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void Copy_Scalar_Dbl(scalar_dbl *from, scalar_dbl *to)
{
  scalar_dbl *f,*t;
  f = from;
  t = to;
  do
    {
      assert(t);
      assert(f);
      t->v = f->v;
      t = t->next;
      f = f->next;
    }
  while(f);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

scalar_dbl *Duplicate_Scalar_Dbl(scalar_dbl *from)
{
  scalar_dbl *to,*t,*f;

  to = (scalar_dbl *)mCalloc(1,sizeof(scalar_dbl));

  t = to;
  f = from;
  do
    {

      t->v = f->v;
      f = f->next;
      if(f) t->next = (scalar_dbl *)mCalloc(1,sizeof(scalar_dbl));
      t = t->next;
    }
  while(f);

  return(to);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

phydbl Scalar_Elem(int pos, scalar_dbl *scl)
{
  scalar_dbl *loc;
  loc = scl;
  while(--pos >= 0) loc = loc->next;
  assert(loc);
  return(loc->v);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void *Linked_List_Elem(int pos, t_ll *ll)
{
  t_ll *loc;
  
  if(ll == NULL) return NULL;
  
  loc = ll->head;
  while(--pos >= 0) 
    {
      assert(loc);
      loc = loc->next;
    }
  assert(loc);
  return(loc->v);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

int Scalar_Len(scalar_dbl *scl)
{
  int len;
  scalar_dbl *loc;

  if(!scl) return 0;

  loc = scl;
  len = 0;
  do 
    {
      len++;
      loc = loc->next;
    }
  while(loc != NULL);

  return(len);

}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

int Linked_List_Len(t_ll *list)
{
  int len;
  t_ll *loc;

  if(list == NULL) return 0;

  loc = list->head;
  len = 0;
  do 
    {
      len++;
      loc = loc->next;
    }
  while(loc != NULL);

  return(len);

}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void Push_Bottom_Linked_List(void *what, t_ll **list, bool remove_duplicates)
{
  t_ll *new,*ll;

  /* { */
  /*   if(*list) */
  /*     { */
  /*       t_node *n = what; */
  /*       t_node *m = (*list)?(*list)->tail->v:NULL; */
  /*       t_node *o = (*list)?(*list)->head->v:NULL; */
  /*       printf("\n. before: push node %d bot %d head: %d --",n?n->num:-1,m?m->num:-1,o?o->num:-1); */
  /*     } */
  /*   else */
  /*     { */
  /*       printf("\n. "); */
  /*     } */
  /* } */

  /* if(*list == NULL) printf("\n"); */
  
  new = (t_ll *)mCalloc(1,sizeof(t_ll));  
  new->v = (void *)what;
  
  /* t_node *n = what; */
  /* printf("\n push node %p list: %p %p",new,*list,(*list)?(*list)->head:NULL); */


  // First elem of list
  if(*list == NULL) 
    {
      *list     = new;
      new->tail = new;
      new->head = new;
      new->next = NULL;
      new->prev = NULL;
    }
  else
    {
      ll = (*list)->head;

      if(remove_duplicates == YES)
        {
          do
            {
              if(ll->v == what)
                {
                  Free(new);
                  return; // 'what' already in list 
                }
              ll = ll->next;
            }
          while(ll);
        }

      new->prev = (*list)->tail;
      (*list)->tail->next = new;
      new->next = NULL;
      new->head = (*list)->head;      

      ll = (*list)->head;
      do
        {
          ll->tail = new;
          ll = ll->next;
        }
      while(ll);
    }
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void Remove_From_Linked_List(t_ll *elem, void *val, t_ll **list)
{
  t_ll *ll;
  
  if(*list == NULL) return;
  
  ll = (*list)->head;

  /* t_node *n = elem ? elem->v : val; */

  do
    {
      if((elem && ll == elem) || (val && ll->v == val))
        {
          if(ll == (*list)->head && ll != (*list)->tail)
            {
              // Re-initialise head of list
              t_ll *mm,*newhead;
              mm = (*list);
              newhead = (*list)->head->next;
              do { mm->head = newhead; mm = mm->next; } while(mm);              
              (*list) = (*list)->head;
              (*list)->head->prev = NULL;
            }
          else if(ll != (*list)->head && ll == (*list)->tail) 
            {
              // Re-initialise tail of list
              t_ll *mm,*newtail;
              mm = (*list);
              newtail = (*list)->tail->prev;
              do { mm->tail = newtail; mm = mm->next; } while(mm);              
              (*list)->tail->next = NULL;
            }
          else if(ll == (*list)->head && ll == (*list)->tail) 
            {
              (*list)->tail = NULL;
              (*list)->head = NULL;
              (*list)->next = NULL;
              (*list)->prev = NULL;
              (*list) = NULL;
              /* printf("\n. free %p",ll); */
              Free(ll);
              return;
            }
          else
            {
              ll->prev->next = ll->next;
              ll->next->prev = ll->prev;
            }

          ll->next = NULL;
          ll->prev = NULL;
          ll->head = NULL;
          ll->tail = NULL;
          /* printf("\n. free %p",ll); */
          Free(ll);
          return;
        }
      ll = ll->next;
    }
  while(ll != NULL);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

t_ll *Get_List_Of_Reachable_Tips(t_node *a, t_node *d, t_tree *tree)
{
  t_ll *list;
  Get_List_Of_Reachable_Tips_Post(a,d,&list,tree);
  return list;
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void Get_List_Of_Reachable_Tips_Post(t_node *a, t_node *d, t_ll **list, t_tree *tree)
{
  if(d->tax)
    {
      /* printf("\n. push %d list: %p",d->num,list->head); */
      Push_Bottom_Linked_List(d,list,YES);
      return;
    }
  else
    {
      int i;
      for(i=0;i<3;i++)
        {
          if(d->v[i] != a && d->b[i] != tree->e_root)
            Get_List_Of_Reachable_Tips_Post(d,d->v[i],list,tree);
        }
    }  
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
// tips0: first list of tips. tips1: second list of tips
phydbl Length_Of_Path_Between_List_Of_Tips(t_ll *tips0, t_ll *tips1, matrix *mat)
{
  phydbl d,n;
  t_ll *x, *y;
  t_node *nx, *ny;

  d = 0.0;
  n = 0;

  /* Print_Mat(mat); */

  // Add all distances between tips in distinct lists
  x = tips0->head;
  do
    {
      y = tips1->head;
      do
        {
          nx = (t_node *)x->v;
          ny = (t_node *)y->v;

          d += mat->dist[nx->c_seq->num][ny->c_seq->num];

          /* printf("\n. nx: %d ny: %d [%p %p] d: %G", */
          /*        nx->c_seq->num, */
          /*        ny->c_seq->num, */
          /*        x->head, */
          /*        y->head, */
          /*        mat->dist[nx->c_seq->num][ny->c_seq->num]); fflush(NULL); */

          n++; //  number of pairs of tips, each tip in different lists
          y = y->next;
        }
      while(y);
      x = x->next;
    }
  while(x);
  
  // Remove distances between tips in tips0
  x = tips0->head;
  do
    {
      y = tips0->head;
      do
        {
          nx = (t_node *)x->v;
          ny = (t_node *)y->v;
          d -= mat->dist[nx->c_seq->num][ny->c_seq->num];
          y = y->next;
        }
      while(y);
      x = x->next;
    }
  while(x);


  // Remove distances between tips in tips1
  x = tips1->head;
  do
    {
      y = tips1->head;
      do
        {
          nx = (t_node *)x->v;
          ny = (t_node *)y->v;
          d -= mat->dist[nx->c_seq->num][ny->c_seq->num];
          y = y->next;
        }
      while(y);
      x = x->next;
    }
  while(x);

  return d/(phydbl)n;
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void Random_Walk_Along_Tree_On_Radius(t_node *a, t_node *d, t_edge *b, phydbl *radius, t_edge **target_edge, t_node **target_nd, phydbl *target_time, t_tree *tree)
{
  assert(tree->rates);
  
  phydbl delta,ta,td,u;

  /* if(tree->rates->nd_t[a->num] < tree->rates->nd_t[d->num]) */
  /*   { */
  /*     printf("\n. a: %d d: %d radius: %f l: %f [%f] -- %f | %f [%d]", */
  /*            a->num, */
  /*            d->num, */
  /*            *radius, */
  /*            b->l->v, */
  /*            FABS(tree->rates->nd_t[a->num] - tree->rates->nd_t[d->num]) * tree->rates->clock_r * tree->rates->br_r[d->num], */
  /*            tree->rates->cur_l[d->num], */
  /*            tree->rates->cur_l[a->num], */
  /*            b == tree->e_root); fflush(NULL); */
  /*   } */
  /* else */
  /*   { */
  /*     printf("\n. a: %d d: %d radius: %f l: %f [%f] -- %f | %f [%d]", */
  /*            a->num, */
  /*            d->num, */
  /*            *radius, */
  /*            b->l->v, */
  /*            FABS(tree->rates->nd_t[a->num] - tree->rates->nd_t[d->num]) * tree->rates->clock_r * tree->rates->br_r[a->num], */
  /*            tree->rates->cur_l[d->num], */
  /*            tree->rates->cur_l[a->num], */
  /*            b == tree->e_root); fflush(NULL); */
  /*   } */
  
  delta = *radius;

  if(!(delta > 0.0))
    {
      PhyML_Fprintf(stderr,"\n. delta=%G",delta);
      Generic_Exit(__FILE__,__LINE__,__FUNCTION__);    
    }
  
  (*radius) -= b->l->v;
  if(*radius < 0.0)
    {
      *target_edge = b;
      ta = tree->rates->nd_t[a->num];
      td = tree->rates->nd_t[d->num];

      if(b != tree->e_root)
        {
          if(ta < td)
            {
              *target_time = ta + delta / (tree->rates->clock_r * tree->rates->br_r[d->num]);
              *target_nd = d;
              /* printf("\n< ta: %f td: %f new_time: %f delta: %f c: %f r: %f",ta,td,*target_time,delta,tree->rates->clock_r,tree->rates->br_r[d->num]); */
            }
          else
            {
              *target_time = ta - delta / (tree->rates->clock_r * tree->rates->br_r[a->num]);
              *target_nd = a;
              /* printf("\n> ta: %f td: %f new_time: %f delta: %f c: %f r: %f",ta,td,*target_time,delta,tree->rates->clock_r,tree->rates->br_r[a->num]); */
            }
        }
      else
        {
          // target falls on edge below root leading to node a
          if(delta < tree->rates->cur_l[a->num])
            {
              *target_time = ta - delta / (tree->rates->clock_r * tree->rates->br_r[a->num]);
              *target_nd = a;
              /* printf("\n.. delta: %f l(a): %f time: %f ta: %f",delta,tree->rates->cur_l[a->num],*target_time,ta); */
            }
          else
            {
              u = Uni();
              if(u < .5)
                {
                  // target falls on edge below root leading to node d
                  *target_time = tree->rates->nd_t[tree->n_root->num] + (delta - tree->rates->cur_l[a->num])/(tree->rates->clock_r * tree->rates->br_r[d->num]);
                  *target_nd = d;
                  /* printf("\n<< ta: %f td: %f new_time: %f delta: %f c: %f",ta,td,*target_time,delta,tree->rates->clock_r); */
                }
              else
                {
                  // target falls above root
                  *target_time = tree->rates->nd_t[tree->n_root->num] - (delta - tree->rates->cur_l[a->num])/tree->rates->clock_r;
                  *target_nd = tree->n_root;                  
                  /* printf("\n>> ta: %f td: %f new_time: %f delta: %f c: %f",ta,td,*target_time,delta,tree->rates->clock_r); */
                }
            }
        }
      
      return;
    }
  
  if(d->tax == YES) return;
  else
    {      
      int i,dir1,dir2;
      
      dir1 = dir2 = -1;
      for(i=0;i<3;i++)
        if(d->v[i] != a)
          {
            if(dir1 < 0) dir1 = i;
            else         dir2 = i;
          }
             
      u = Uni();
      if(u < .5)
        Random_Walk_Along_Tree_On_Radius(d,d->v[dir1],d->b[dir1],radius,target_edge,target_nd,target_time,tree);
      else
        Random_Walk_Along_Tree_On_Radius(d,d->v[dir2],d->b[dir2],radius,target_edge,target_nd,target_time,tree);
    }
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void Table_Top(unsigned int width)
{
  unsigned int i;
  PhyML_Printf("\n\t \u256D");
  for(i=0;i<width;++i) PhyML_Printf("\u2500");
  PhyML_Printf("\u256E ");
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void Table_Row(unsigned int width)
{
  unsigned int i;
  /* PhyML_Printf("\n\t \u2502"); */
  PhyML_Printf("\n\t \u251C");
  for(i=0;i<width;++i) PhyML_Printf("\u2500");
  /* PhyML_Printf("\u2502"); */
  PhyML_Printf("\u2524");
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void Table_Bottom(unsigned int width)
{
  unsigned int i;
  PhyML_Printf("\n\t \u2570");
  for(i=0;i<width;++i) PhyML_Printf("\u2500");
  PhyML_Printf("\u256F ");
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

t_cal *Duplicate_Calib(t_cal *from)
{
  int i;
  t_cal *to;
  
  to = Make_Calibration();

  to->current_clade_idx = from->current_clade_idx;
  to->lower = from->lower;
  to->upper = from->upper;
  to->is_primary = from->is_primary;
  to->clade_list_size = from->clade_list_size;

  to->id = (char *)mCalloc(strlen(from->id)+1,sizeof(char)); 
  strcpy(to->id,from->id);

  to->alpha_proba_list = (phydbl *)mCalloc(from->clade_list_size,sizeof(phydbl));
  to->clade_list = (t_clad **)mCalloc(from->clade_list_size,sizeof(t_clad *));

  for(i=0;i<from->clade_list_size;++i)
    {
      to->alpha_proba_list[i] = from->alpha_proba_list[i];
      to->clade_list[i] = Duplicate_Clade(from->clade_list[i]);
    }
  
  return(to);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

t_clad *Duplicate_Clade(t_clad *from)
{
  int i;
  t_clad *to;

  to = Make_Clade();

  to->id = (char *)mCalloc(strlen(from->id)+1,sizeof(char));
  strcpy(to->id,from->id);
  
  to->n_tax = from->n_tax;
  to->target_nd = from->target_nd;

  to->tax_list = (char **)mCalloc(from->n_tax,sizeof(char *));
  to->tip_list = (t_node **)mCalloc(from->n_tax,sizeof(t_node *));
  
  for(i=0;i<from->n_tax;++i)
    {
      to->tax_list[i] = (char *)mCalloc(strlen(from->tax_list[i])+1,sizeof(char));
      strcpy(to->tax_list[i],from->tax_list[i]);
      to->tip_list[i] = from->tip_list[i];
    }
  
  return to;
  
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
