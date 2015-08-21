/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include "nexus.h"

void Find_Nexus_Com(char *token, nexcom **found_com, nexparm **default_parm, nexcom **com_list)
{
  int i,j,tokenlen,ndiff;
    
  For(i,N_MAX_NEX_COM) 
    {
      tokenlen = strlen(token);
      ndiff = -1;
      if(tokenlen && (tokenlen == strlen(com_list[i]->name)))
	{
	  ndiff = 0;
	  For(j,tokenlen)
	    {
	      Lowercase(token+j);
	      Lowercase(com_list[i]->name+j);
	      if(token[j] != com_list[i]->name[j]) ndiff++;
	    }
	}
      if(!ndiff) { *found_com = com_list[i]; break; }
    }

  if(*found_com && (*found_com)->nparm) *default_parm = (*found_com)->parm[0];

  /* if(*found_com) PhyML_Printf("\n. Found command '%s'.\n",(*found_com)->name); */
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Find_Nexus_Parm(char *token, nexparm **found_parm, nexcom *curr_com)
{
  int i,j;
  int tokenlen;
  int ndiff;

  if(!curr_com)
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("");
    }

  For(i,curr_com->nparm)
    {
      tokenlen = strlen(token);
      ndiff = -1;
      if(tokenlen == strlen(curr_com->parm[i]->name))
	{
	  ndiff = 0;
	  For(j,tokenlen)
	    {
	      Lowercase(token+j);
	      Lowercase(curr_com->parm[i]->name+j);
	      if(token[j] != curr_com->parm[i]->name[j]) ndiff++;
	    }
	}
      if(!ndiff) { *found_parm = curr_com->parm[i]; break; }
    }

  /* if(*found_parm) PhyML_Printf("\n. Found parameter '%s'.\n",(*found_parm)->name); */
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


int Read_Nexus_Taxa(char *token, nexparm *curr_parm, option *io)
{

  PhyML_Printf("\n. Skipping 'taxa' block");

  do
    {
      Get_Token(io->fp_in_align,token);
      if(token[0] == ';') break;
    }while(strlen(token) > 0);
  
  fseek(io->fp_in_align,-1*sizeof(char),SEEK_CUR);

  return 1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


int Read_Nexus_Translate(char *token, nexparm *curr_parm, option *io)
{
  int tax_num;
  char *end;

  PhyML_Printf("\n. Reading 'translate' block");
  io->size_tax_names = 0;

  do
    {
      Get_Token(io->fp_in_tree,token);
      if(token[0] == ';') break;
      tax_num = (int)strtol(token,&end,10);
      if(*end =='\0' && token[0])
	{
	  io->size_tax_names++;

	  io->short_tax_names = (char **)realloc(io->short_tax_names,io->size_tax_names*sizeof(char *));
	  io->short_tax_names[io->size_tax_names-1] = (char *)mCalloc(strlen(token)+1,sizeof(char));
	  sprintf(io->short_tax_names[io->size_tax_names-1],"%d",tax_num);

	  Get_Token(io->fp_in_tree,token);

	  io->long_tax_names = (char **)realloc(io->long_tax_names,io->size_tax_names*sizeof(char *));
	  io->long_tax_names[io->size_tax_names-1] = (char *)mCalloc(strlen(token)+1,sizeof(char));
	  strcpy(io->long_tax_names[io->size_tax_names-1],token);

/* 	  printf("\n. Copying %s number %d",io->long_tax_names[io->size_long_tax_names-1],tax_num-1); */
	}
    }while(strlen(token) > 0);
  
  fseek(io->fp_in_tree,-1*sizeof(char),SEEK_CUR);

  return 1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


int Read_Nexus_Matrix(char *token, nexparm *curr_parm, option *io)
{

  if(io->interleaved) io->data = Read_Seq_Interleaved(io);
  else                io->data = Read_Seq_Sequential(io);

  fseek(io->fp_in_align,-1*sizeof(char),SEEK_CUR);

  return 1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


int Read_Nexus_Tree(char *token, nexparm *curr_parm, option *io)
{
  io->treelist->tree = (t_tree **)realloc(io->treelist->tree,(io->treelist->list_size+1)*sizeof(t_tree *));
  io->tree = Read_Tree_File_Phylip(io->fp_in_tree);
  if(!(io->treelist->list_size%10) && io->treelist->list_size > 1) 
    {
      PhyML_Printf("\n. Reading tree %d",io->treelist->list_size);
      if(io->tree->n_root) PhyML_Printf(" (that is a rooted tree)");
      else                 PhyML_Printf(" (that is an unrooted tree)");
    }
  io->treelist->tree[io->treelist->list_size] = io->tree;
  io->treelist->list_size++;
  fseek(io->fp_in_tree,-1*sizeof(char),SEEK_CUR);
  return 1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


int Read_Nexus_Begin(char *token, nexparm *curr_parm, option *io)
{
  if(token[0] == '=') return 0;

  if(!curr_parm)
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("");
    }

  if(!strcmp(curr_parm->name,"data") || !strcmp(curr_parm->name,"trees")) 
    PhyML_Printf("\n. Reading '%s' block.\n",curr_parm->value);
  else
    {
      PhyML_Printf("\n. The '%s' block type is not supported by PhyML. Sorry.\n",curr_parm->name);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("");
    }

  return 1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


int Read_Nexus_Dimensions(char *token, nexparm *curr_parm, option *io)
{
  if(token[0] == '=') return 0;

  if(!curr_parm)
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("");
    }
  
  strcpy(curr_parm->value,token);

  if(!strcmp(curr_parm->name,"ntax"))
    {
      sscanf(curr_parm->value,"%d",&(io->n_otu));
    }

  if(!strcmp(curr_parm->name,"nchar"))
    {
      sscanf(curr_parm->value,"%d",&(io->init_len));
    }
  return 1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


int Read_Nexus_Format(char *token, nexparm *curr_parm, option *io)
{
  int i;  
  
  if(token[0] == '=') return 0;
  
  if(!curr_parm)
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("");
    }

  For(i,strlen(token)) Lowercase(token+i);

  strcpy(curr_parm->value,token);


  /* printf("\n. >> %s",curr_parm->value); */
    
  if(!strcmp(curr_parm->name,"datatype"))
    {
      if(!strcmp(curr_parm->value,"standard"))
	{
	  io->datatype = GENERIC;
	  io->mod->whichmodel = JC69;
	  io->mod->s_opt->opt_kappa  = NO;
	  io->mod->s_opt->opt_lambda = NO;
	  io->mod->ns = 2;
	  io->alphabet[0][0] = '0'; io->alphabet[0][1] = '\0';
	  io->alphabet[1][0] = '1'; io->alphabet[1][1] = '\0';
	}

      else if(!strcmp(curr_parm->value,"dna"))
	{
	  io->datatype = NT;
	  io->mod->ns = 4;
	}

      else if(!strcmp(curr_parm->value,"rna"))
	{
	  io->datatype = NT;
	  io->mod->ns = 4;
	}

      else if(!strcmp(curr_parm->value,"nucleotide"))
	{
	  io->datatype = NT;
	  io->mod->ns = 4;
	}

      else if(!strcmp(curr_parm->value,"protein"))
	{
	  io->datatype = AA;
	  io->mod->ns = 20;
	}
      
      else if(!strcmp(curr_parm->value,"continuous"))
	{
	  PhyML_Printf("\n== The 'continuous' format is not supported by PhyML. Sorry.\n");
	  PhyML_Printf("\n== Err. in file %s at line %d\n",__FILE__,__LINE__);
	  Exit("");
	}
    }

  else if(!strcmp(curr_parm->name,"missing"))
    {
      PhyML_Printf("\n== The 'missing' subcommand is not supported by PhyML. Please remove it from the NEXUS file.");
      PhyML_Printf("\n== Note that the characters 'X', '?' and '-' will be considered as indels by default."); 
      PhyML_Printf("\n== Err. in file %s at line %d\n",__FILE__,__LINE__);
      Exit("");
    }

  else if(!strcmp(curr_parm->name,"gap"))
    {
      PhyML_Printf("\n== The 'gap' subcommand is not supported by PhyML. Please remove it from the NEXUS file.");
      PhyML_Printf("\n== Note that the characters 'X', '?' and '-' will be considered as indels by default."); 
      PhyML_Printf("\n== Err. in file %s at line %d\n",__FILE__,__LINE__);
      Exit("");
    }

  else if(!strcmp(curr_parm->name,"symbols"))
    {
      if(*token != '"' || *(token+strlen(token)-1) != '"')
	{
	  PhyML_Printf("\n== Symbols list is supposed to be displayed between quotation marks (e.g., \"ACTG\").\n");
	  PhyML_Printf("\n== Err. in file %s at line %d\n",__FILE__,__LINE__);
	  Exit("");
	}


      int i,has_spaces,state_len;

      i          = 0;
      has_spaces = 0;     
      token++; /* Get rid of the first '"' character */
      while(token[i] != '"')  { if(token[i] == ' ') { has_spaces = 1; break; } i++; }

      io->mod->ns = 0;
      if(!has_spaces)
	{
	  while(token[i] != '"') 
	    { 
	      io->alphabet[io->mod->ns][0] = token[i]; 
	      io->alphabet[io->mod->ns][1] = '\0'; 
	      io->mod->ns++;
	      i++;
	      if(io->mod->ns > T_MAX_ALPHABET)
		{
		  PhyML_Printf("\n== The alphabet cannot contain more than %d characters. Sorry.",T_MAX_ALPHABET);
		  PhyML_Printf("\n== Err. in file %s at line %d\n",__FILE__,__LINE__);
		  Exit("");
		}
	    }
	}
      else
	{
	  i = 0;
	  do
	    {
	      state_len = 0;
	      while(token[i] != ' ' && token[i] != '"') 
		{ 
		  io->alphabet[io->mod->ns][state_len] = token[i];
		  state_len++;
		  i++;
		  if(state_len > T_MAX_STATE)
		    {
		      PhyML_Printf("\n== A state cannot contain more than %d characters. Sorry.\n",T_MAX_STATE);
		      PhyML_Printf("\n== Err. in file %s at line %d\n",__FILE__,__LINE__);
		      Exit("");
		    }
		}
	      
	      io->alphabet[io->mod->ns][state_len] = '\0';
	      io->mod->ns++;
	      if(token[i] != '"') i++;
	    }
	  while(token[i] != '"');
	}

      int len;
      len = strlen(io->alphabet[0]);
      For(i,io->mod->ns)
	{
	  if(strlen(io->alphabet[i]) != len)
	    {
	      PhyML_Printf("\n== All character states defined in the symbol list are supposed to have the same length.\n");
	      PhyML_Printf("\n== Er.r in file %s at line %d\n",__FILE__,__LINE__);
	      Exit("");
	    }
	}
      io->state_len = len;      

/*       For(i,io->mod->ns) PhyML_Printf("\n. '%s'",io->alphabet[i]); */
    }
  
  else if(!strcmp(curr_parm->name,"equate"))
    {
      PhyML_Printf("\n== PhyML does not recognize the command '%s' yet. Sorry.",curr_parm->name);
      PhyML_Printf("\n== Err. in file %s at line %d\n",__FILE__,__LINE__);
      Exit("");
    }
  
  else if(!strcmp(curr_parm->name,"matchchar"))
    {
      PhyML_Printf("\n== PhyML does not recognize the command '%s' yet. Sorry.",curr_parm->name);
      PhyML_Printf("\n== Err. in file %s at line %d\n",__FILE__,__LINE__);
      Exit("");
    }

  else if(!strcmp(curr_parm->name,"items"))
    {
      PhyML_Printf("\n== PhyML does not recognize the command '%s' yet. Sorry.",curr_parm->name);
      PhyML_Printf("\n== Err. in file %s at line %d\n",__FILE__,__LINE__);
      Exit("");
    }

  else if(!strcmp(curr_parm->name,"interleave"))
    {
      io->interleaved = YES;
    }

  return 1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


int Read_Nexus_Eliminate(char *token, nexparm *curr_parm, option *io)
{
  if(token[0] == '=') return 0;

  PhyML_Printf("\n== 'Eliminate' command is not supported by PhyML. Sorry.");
  PhyML_Printf("\n== Err. in file %s at line %d\n",__FILE__,__LINE__);
  Exit("");

  return 1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


int Read_Nexus_Taxlabel(char *token, nexparm *curr_parm, option *io)
{
  if(token[0] == '=') return 0;

  PhyML_Printf("\n== 'Taxlabels' command is not supported by PhyML. Sorry.");
  PhyML_Printf("\n== Err. in file %s at line %d\n",__FILE__,__LINE__);
  Exit("");

  return 1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


int Read_Nexus_Charstatelabels(char *token, nexparm *curr_parm, option *io)
{

  if(token[0] == '=') return 0;

  PhyML_Printf("\n== 'CharStateLabels' command is not supported by PhyML. Sorry.");
  PhyML_Printf("\n== Err. in file %s at line %d\n",__FILE__,__LINE__);
  Exit("");

  return 1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


int Read_Nexus_Charlabels(char *token, nexparm *curr_parm, option *io)
{
  if(token[0] == '=') return 0;

  PhyML_Printf("\n== 'CharLabels' command is not supported by PhyML. Sorry.");
  PhyML_Printf("\n== Err. in file %s at line %d\n",__FILE__,__LINE__);
  Exit("");

  return 1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


int Read_Nexus_Statelabels(char *token, nexparm *curr_parm, option *io)
{
  if(token[0] == '=') return 0;

  PhyML_Printf("\n== 'StateLabels' command is not supported by PhyML. Sorry.");
  PhyML_Printf("\n== Err. in file %s at line %d\n",__FILE__,__LINE__);
  Exit("");

  return 1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

