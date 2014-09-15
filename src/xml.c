#include "xml.h"

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

xml_node *XML_Load_File(FILE *fp)
{
  int c;
  char *buffer,*bufptr;
  int bufsize;
  xml_node *parent,*node;

  buffer = (char *)mCalloc(T_MAX_XML_TAG,sizeof(char));

  bufsize = T_MAX_XML_TAG;
  bufptr  = buffer;
  parent  = NULL;
  node    = NULL;

  while((c = fgetc(fp)) != EOF)
    {
      if(c == '<' && bufptr > buffer) 
        {
          *bufptr = '\0';

          /* PhyML_Printf("\n. Read value '%s' for node '%s'",buffer,node->name); */
          /* fflush(NULL); */

          XML_Set_Node_Value(node,buffer);
          bufptr = buffer;
        }
              
      if(c == '<')
        {
          bufptr = buffer;

          while((c = fgetc(fp)) != EOF)
            {
              if(isspace(c) != NO || c == '>' || (c == '/' && bufptr > buffer)) break; // End of open or close tag
              else if(c == '<')
                {
                  Exit("\n== Bare < in element!");
                }             
              else if(XML_Add_Character(c,&bufptr,&buffer,&bufsize))
                {
                  PhyML_Printf("\n== Err in file %s at line %d\n",__FILE__,__LINE__);
                  Exit("\n");     
                }
            }

          *bufptr = '\0';
          
          if(!strcmp(buffer,"!--")) // Get the rest of the comment
            {
              while((c = fgetc(fp)) != EOF)
                {
                  
                  if(c == '>' && bufptr > (buffer + 4) && bufptr[-3] != '-' &&
                     bufptr[-2] == '-' && bufptr[-1] == '-') break;
                  else if(XML_Add_Character(c,&bufptr,&buffer,&bufsize))
                    {
                      PhyML_Printf("\n== Err in file %s at line %d\n",__FILE__,__LINE__);
                      Exit("\n");         
                    }
                }
              *bufptr = '\0';

              if(c != '>')
                {
                  PhyML_Printf("\n== Early EOF in comment node.");
                  Exit("\n");     
                }             
            }     
          else if(buffer[0] == '/') // Close tag
            {
              if(strcmp(buffer+1,parent->name))
                {
                  PhyML_Printf("\n== Opened tag with name '%s' and closed it with '%s'...",node->name,buffer+1);
                  PhyML_Printf("\n== Err in file %s at line %d\n",__FILE__,__LINE__);
                  Exit("\n");
                }

              /* printf("\n. Closing node with name '%s'",node->name); */

              if(node->parent)
                {
                  parent = parent->parent;
                  node   = parent;
                }
            }
          else if(buffer[0] == '?')
            {
              while((c = fgetc(fp)) != EOF)
                {
                  if (c == '>' && bufptr > buffer && bufptr[-1] == '?')
                    break;
                  else if (XML_Add_Character(c, &bufptr, &buffer, &bufsize))
                    {
                      PhyML_Printf("\n== Err in file %s at line %d\n",__FILE__,__LINE__);
                      Exit("\n");         
                    }
                }

              if(c != '>')
                {
                  PhyML_Printf("\n== An error occurred when reading the processing instruction.");
                  PhyML_Printf("\n== Err in file %s at line %d\n",__FILE__,__LINE__);
                  Exit("\n");     
                }

              *bufptr = '\0';

            }
          else // Open tag
            {
              node = XML_Make_Node(buffer);
              XML_Init_Node(parent,node,buffer);
              if(!parent) parent = node;

              if(isspace(c) != NO) c=XML_Parse_Element(fp,node);
              else if(c == '/')
                {
                  if((c=fgetc(fp)) != '>')
                    {
                      PhyML_Printf("\n== Expected '>' but read '%c' instead",c);
                      Exit("\n");
                    }
                  c = '/';
                }

              if(c != '/') parent = node;

              buffer[0] = '\0';
            }     
          bufptr = buffer;
        }
      else if(isspace(c) == NO)
        {
          if(XML_Add_Character(c,&bufptr,&buffer,&bufsize))
            {
              PhyML_Printf("\n== Err in file %s at line %d\n",__FILE__,__LINE__);
              Exit("\n");
            }
        }
    }
  Free(buffer);
  return node;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int XML_Add_Character(int c, char  **bufptr, char **buffer, int *bufsize)
{
  char *newbuffer;

  if(*bufptr >= (*buffer + *bufsize - 4))
    {
      // Increase the size of the buffer...
      
      if (*bufsize < 1024)
        (*bufsize) *= 2;
      else
        (*bufsize) += 1024;

    if((newbuffer = realloc(*buffer, *bufsize)) == NULL)
      {
        Free(*buffer);
        PhyML_Printf("\n== Unable to expand string buffer to %d bytes!", *bufsize);
        Exit("\n");
      }
    
    *bufptr = newbuffer + (*bufptr - *buffer);
    *buffer = newbuffer;
  }

  /* *(*bufptr)++ = tolower(c); */
  *(*bufptr)++ = c;
  return 0;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int XML_Parse_Element(FILE *fp, xml_node *n)
{
  int c;
  int quote;
  char *name, *value, *ptr;
  int namesize, valsize;

  name  = (char *)mCalloc(64,sizeof(char));
  value = (char *)mCalloc(64,sizeof(char));
  
  namesize = 64;
  valsize  = 64;
  
  while((c = fgetc(fp)) != EOF)
    {

      if(isspace(c) != NO) continue;

      if(c == '/') // End of tag
        {
          /* printf("\n. Closing node '%s'.",n->name); */

          quote = fgetc(fp);
          if(quote != '>')
            {
              PhyML_Printf("\n== Expected '>' after '%c' but read '%c' instead",c,quote);
              Exit("\n");
            }
          break;
        }
      else if(c == '<')
        {
          Exit("\n== Bare < in element!");        
        }
      else if(c == '>') // End of tag
        {
          break;
        }

      name[0] = c;
      ptr     = name + 1;

      if(c == '\"' || c == '\'') // Name is in quotes
        {
          quote = c;

          while((c = fgetc(fp)) != EOF)
            {
              if(XML_Add_Character(c,&ptr,&name,&namesize))
                {
                  PhyML_Printf("\n== Err in file %s at line %d\n",__FILE__,__LINE__);
                  Exit("\n");
                }
              if(c == quote) break;
            }     
        }
      else // Name not in quotes
        {
          while((c = fgetc(fp)) != EOF)
            {
              if(isspace(c) != NO || c == '=' || c == '/' || c == '>' || c == '?')
                break;
              else
                {
                  if(XML_Add_Character(c,&ptr,&name,&namesize))
                    {
                      PhyML_Printf("\n== Err in file %s at line %d\n",__FILE__,__LINE__);
                      Exit("\n");
                    }
                }         
            }
        }
      
      *ptr = '\0';
            
      while(c != EOF && isspace(c) != NO) c = fgetc(fp);

      if(c == '=') // Read the attribute value
        {
          while((c = fgetc(fp)) != EOF && isspace(c) != NO);

          if(c == EOF)
            {
              PhyML_Printf("\n== Missing value in attribute.");
              PhyML_Printf("\n== Err in file %s at line %d\n",__FILE__,__LINE__);
              Exit("\n");
            }

          if(c == '\'' || c == '\"')
            {
              quote = c;
              ptr   = value;

              while((c = fgetc(fp)) != EOF)
                {
                  if(c == quote) break;
                  else
                    {
                      if(XML_Add_Character(c,&ptr,&value,&valsize))
                        {
                          PhyML_Printf("\n== Err in file %s at line %d\n",__FILE__,__LINE__);
                          Exit("\n");
                        }
                    }
                }
              *ptr = '\0';
            }
          else
            {
              value[0] = c;
              ptr      = value + 1;
              
              while((c = fgetc(fp)) != EOF)
                {
                  if(isspace(c) != NO || c == '=' || c == '/' || c == '>')
                    break;
                  else
                    {
                      if(XML_Add_Character(c,&ptr,&value,&valsize))
                        {
                          PhyML_Printf("\n== Err in file %s at line %d\n",__FILE__,__LINE__);
                          Exit("\n");
                        }                     
                    }
                }             
            }
        }

      /* printf("\n. Setting attribute '%s=%s' to node '%s'",name,value,n->name); */
      XML_Set_Attribute(n,name,value);

      if(c == '>') break;


    }
  Free(name);
  Free(value);

  /* printf("\n. Return '%c'\n",c); */
  return(c);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

xml_attr *XML_Search_Attribute(xml_node *n, char *target_attr_name)
{
  xml_attr *attr;

  attr = n->attr;
  do
    {
      if(!strcmp(attr->name,target_attr_name)) return attr;
      attr = attr->next;
    }
  while(attr);

  return(NULL);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int XML_Set_Attribute(xml_node *n, char *attr_name, char *attr_value)
{
  xml_attr *prev;
  char *s;

  prev = NULL;
  while(n->attr != NULL) 
    {
      prev    = n->attr;
      n->attr = n->attr->next;
    }

  n->attr = XML_Make_Attribute(prev,attr_name,attr_value);
  XML_Init_Attribute(n->attr);
  n->n_attr++;

  // rewind
  while(n->attr->prev != NULL) n->attr = n->attr->prev; 

  s = To_Lower_String(attr_name);
  if(!strcmp(s,"id"))
    {
      XML_Set_Node_Id(n,attr_value);
      /* printf("\n. Node '%s' id is '%s'",n->name,n->id); */
    }
  Free(s);

  return(0);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int XML_Set_Node_Id(xml_node *n, char *id)
{
  XML_Make_Node_Id(n,id);
  strcpy(n->id,id);
  return(0);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


int XML_Set_Node_Value(xml_node *n, char *val)
{
  XML_Make_Node_Value(n,val);
  strcpy(n->value,val);
  return(0);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

xml_node *XML_Search_Node_Generic(char *nd_name, char *attr_name, char *attr_val, int skip, xml_node *node)
{

  xml_node *match;

  /* if(nd_name) printf("\n. [1] nd_name:%s attr_name:%s attr_val:%s \n", nd_name, attr_name, attr_val); */
  /* else  printf("\n. attr_name:%s attr_val:%s \n", attr_name, attr_val); */
  
  /* printf("\n. name:%s child:%s next:%s ", */
  /*     node?node->name:"xx", */
  /*     node->child?node->child->name:"xx", */
  /*     node->next?node->next->name:"xx"); fflush(NULL); */


  match = NULL;
  if(skip == NO && nd_name && attr_name && attr_val)
    {
      if(!strcmp(nd_name, node -> name))
        {
          xml_attr *attr = XML_Search_Attribute(node, attr_name);
          if(attr && !strcmp(attr -> value, attr_val)) match = node;
        }
    }
  else if(skip == NO && !nd_name && attr_name && attr_val)
    {
      xml_attr *attr = XML_Search_Attribute(node, attr_name);
      if(attr && !strcmp(attr -> value, attr_val)) match = node;
    }
  else if(skip == NO && nd_name && !attr_name && attr_val)
    {
      if(!strcmp(nd_name, node -> name))
        {
          do
            {
              if(!strcmp(node -> attr -> value, attr_val)) 
                {
                  match = node;
                  break;
                }
              node -> attr = node -> attr -> next;
              if(!node -> attr) break;
            }
          while(1);
        }
    }
  else if(skip == NO && nd_name && attr_name && !attr_val)
    {
      if(!strcmp(nd_name, node -> name))
        {
          do
            {
              if(!strcmp(node -> attr -> name, attr_name)) 
                {
                  match = node;
                  break;
                }
              node -> attr = node -> attr -> next;
              if(!node -> attr) break;
            }
          while(1);
        }
    }
  else if(skip == NO && nd_name && !attr_name && !attr_val)
    {
      if(!strcmp(nd_name, node -> name)) match = node;
    }
  else if(skip == NO && !nd_name && attr_name && !attr_val)
    {
      xml_attr *attr = XML_Search_Attribute(node, attr_name);
      if(attr) match = node;
    }
  else if(skip == NO && !nd_name && !attr_name && attr_val)
    {
      do
        {
          if(!strcmp(node -> attr -> value, attr_val)) 
            {
              match = node;
              break;
            }
          node -> attr = node -> attr -> next;
          if(!node -> attr) break;
        }
      while(1);
    }

  // If node has a child, node = child, else if node has next, node = next, else if node
  // has parent, node = parent->next else node = NULL

  if(match) return(match);
  if(node -> child)
    {
      match = XML_Search_Node_Generic(nd_name, attr_name, attr_val, NO, node -> child);
    }
  if(!match && node -> next)
    {
      match = XML_Search_Node_Generic(nd_name, attr_name, attr_val, NO, node -> next);
    }
  if(match == NULL && node -> parent)
    {
      if(node -> parent == NULL) // Reached the root
        {
          PhyML_Printf("\n== Could not find a node with name '%s'.", attr_name);
          Exit("\n");
        }
      return NULL;
    }
  return match;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

xml_node *XML_Search_Node_Name(char *name, int skip, xml_node *node)
{

  xml_node *match;
  
  /* printf("\n. name:%s child:%s next:%s ", */
  /*     node?node->name:"xx", */
  /*     node->child?node->child->name:"xx", */
  /*     node->next?node->next->name:"xx"); fflush(NULL); */


  match = NULL;
  if(skip == NO && !strcmp(node->name,name)) match = node;
  else
    {
      // If node has a child, node = child, else if node has next, node = next, else if node
      // has parent, node = parent->next else node = NULL
      if(node->child)
        {
          match = XML_Search_Node_Name(name,NO,node->child);
        }
      if(match == NULL && node->next)
        {
          match = XML_Search_Node_Name(name,NO,node->next);
        }
      if(match == NULL && node->parent)
        {
          if(node->parent == NULL) // Reached the root
            {
              PhyML_Printf("\n== Could not find a node with name '%s'.",name);
              Exit("\n");
            }
          return NULL;
        }
    }
  return match;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

xml_node *XML_Search_Node_ID(char *id, int skip, xml_node *node)
{
  xml_node *match;
  
  if(!node)
    {
      PhyML_Printf("\n== Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");         
    }
      

  match = NULL;
  if(skip == NO && node->id && !strcmp(node->id,id)) match = node;
  else
    {
      // If node has a child, node = child, else if node has next, node = next, else if node
      // has parent, node = parent->next else node = NULL
      if(node->child)
        {
          match = XML_Search_Node_ID(id,NO,node->child);
        }
      if(match == NULL && node->next)
        {
          match = XML_Search_Node_ID(id,NO,node->next);
        }
      if(match == NULL && node->parent)
        {
          if(node->parent == NULL) // Reached the root
            {
              PhyML_Printf("\n== Could not find a node with id '%s'.",id);
              Exit("\n");
            }
          
          return NULL;
        }
    }
  return match;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

xml_node *XML_Search_Node_Attribute_Value(char *attr_name, char *value, int skip, xml_node *node)
{
  xml_node *match;

  
  if(!node)
    {
      PhyML_Printf("\n== node: %p attr: %p",node,node?node->attr:NULL);
      PhyML_Printf("\n== Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");         
    }
  
  match = NULL;

  if(skip) 
    {
      match = XML_Search_Node_Attribute_Value(attr_name, value, NO, node->child);
      return match;
    }

  if(skip == NO && node->attr)
    {
      xml_attr *attr;

      attr = node->attr;
      do
        {
          if(!strcmp(attr->name,attr_name) && 
             !strcmp(attr->value,value)) 
            {
              match = node;
              break;
            }
          attr = attr->next;
          if(!attr) break;
        }
      while(1);
    }

  if(match) return(match);

  if(node->child)
    {
      match = XML_Search_Node_Attribute_Value(attr_name,value,NO,node->child);
      return match;
    }
  if(node->next && !match)
    {
      match = XML_Search_Node_Attribute_Value(attr_name,value,NO,node->next);
      return match;
    }
  return NULL;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

char *XML_Get_Attribute_Value(xml_node *node, char *attr_name)
{
  xml_attr *attr;
  
  attr = node->attr;

  while(attr && strcmp(attr->name,attr_name))
    {
      attr = attr->next;
      
      /* if(attr == NULL) */
      /*        { */
      /*          PhyML_Printf("\n== Could not find an attribute with name '%s' in node with name '%s'.",id,node->name); */
      /*          Exit("\n"); */
      /*        } */
    }

  return(attr?attr->value:NULL);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int XML_Validate_Attr_Int(char *target, int num, ...)
{
  va_list args;                     
  int i;
  char *s,*sc_s;
  char *sc_target;
  
  sc_target = To_Lower_String(target);

  va_start(args,num);           
  For(i,num+1)
    {
      s = va_arg(args, char *); 
      sc_s = To_Lower_String(s);
      if(!strcmp(sc_s,sc_target)) 
        {
          Free(sc_s);
          break;
        }
      Free(sc_s);
    }
  va_end(args);

  if(i == num+1) 
    {
      i = -1;
      PhyML_Printf("\n== Attribute value '%s' is not valid",target);
      Exit("\n");
    }

  Free(sc_target);

  return(i);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void XML_Check_Siterates_Node(xml_node *parent)
{
  xml_node *n;
  int n_weights_nodes;
  char *rate_value = NULL;
  phydbl buff;
  int n_zeros;
  char *endptr;

  if(!parent)
    {
      PhyML_Printf("\n== Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");         
    }

  if(strcmp(parent->name,"siterates"))
    {
      PhyML_Printf("\n== Node name '%s' should be 'siterates'",parent->name);
      Exit("\n");
    }
  
  // Check that only one 'weights' node is present
  n_weights_nodes = 0;
  n = parent->child;
  do
    {
      if(!strcmp(n->name,"weights")) n_weights_nodes++;
      if(n_weights_nodes > 1)
        {
          PhyML_Printf("\n== Only one distribution is authorized for 'siterates' nodes.");
          Exit("\n");
        }
      n = n->next;
      if(!n) break;
    }
  while(1);

  // Check that one rate value is set to zero if gamma+inv model is used
  n = XML_Search_Node_Attribute_Value("family","gamma+inv",YES,parent);
  if(!n) return;
  else
    {
      n_zeros = 0;
      n = parent->child;
      do
        {
          if(!strcmp(n->name,"instance"))
            {
              rate_value = NULL;
              rate_value = XML_Get_Attribute_Value(n,"init.value");  

              if(rate_value)
                {
                  errno = 0;
                  buff = strtod(rate_value,&endptr);
                  
                  if(rate_value == endptr || errno == ERANGE)
                    {
                      PhyML_Printf("\n== value: %s",rate_value);
                      PhyML_Printf("\n== Error in reading attribute 'init.value' in node 'instance'.");
                      Exit("\n");
                    }
                  
                  if(buff < 1.E-20) n_zeros++;           
                }
            }
          n = n->next;
          if(!n) break;
        }
      while(1);
      
      if(n_zeros != 1)
        {
          PhyML_Printf("\n== # of zero-rates: %d",n_zeros);
          PhyML_Printf("\n== Exactly one rate value has to be set to zero when using the 'gamma+inv' model.");
          PhyML_Printf("\n== Component id: %s",parent->id);
          Exit("\n");
        }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int XML_Get_Number_Of_Classes_Siterates(xml_node *parent)
{
  xml_node *n;
  int n_classes;

  if(!parent)
    {
      PhyML_Printf("\n== Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");         
    }

  n_classes = 0;
  n = parent->child;
  do
    {
      if(!strcmp(n->name,"instance")) n_classes++;        
      n = n->next;
      if(!n) break;
    }
  while(1);
  
  n = NULL;
  n = XML_Search_Node_Attribute_Value("family","gamma+inv",YES,parent);

  if(!n) return n_classes;
  else return n_classes-1;
}

//////////////////////////////////////////////////////////////

int XML_Siterates_Has_Invariants(xml_node *parent)
{
  xml_node *n;

  if(!parent)
    {
      PhyML_Printf("\n== Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");         
    }
      
  n = NULL;
  n = XML_Search_Node_Attribute_Value("family","gamma+inv",YES,parent);

  if(!n) return NO;
  else return YES;
}

//////////////////////////////////////////////////////////////

void XML_Count_Number_Of_Node_With_ID(char *id, int *count, xml_node *n)
{
  if(!id) return;
  if(n->id && !strcmp(n->id,id)) (*count)++;
  
  if(n->child) XML_Count_Number_Of_Node_With_ID(id,count,n->child);
  if(n->next)  XML_Count_Number_Of_Node_With_ID(id,count,n->next);
    
}

//////////////////////////////////////////////////////////////

void XML_Count_Number_Of_Node_With_Name(char *name, int *count, xml_node *n)
{
  if(!name) return;
  if(n->name && !strcmp(n->name,name)) (*count)++;
  
  if(n->child) XML_Count_Number_Of_Node_With_Name(name,count,n->child);
  if(n->next)  XML_Count_Number_Of_Node_With_Name(name,count,n->next);
    
}

//////////////////////////////////////////////////////////////

void XML_Check_Duplicate_ID(xml_node *n)
{
  int count;
  
  count = 0;
  XML_Count_Number_Of_Node_With_ID(n->id,&count,n);
  
  if(count > 1)
    {
      PhyML_Printf("\n== Node ID'%s' was found in more than once.",n->id);
      PhyML_Printf("\n== Each ID must be unique. Please amend your XML");
      PhyML_Printf("\n== file accordingly.");
      Exit("\n");
    }

  if(n->child) XML_Check_Duplicate_ID(n->child);
  if(n->next) XML_Check_Duplicate_ID(n->next);
}

//////////////////////////////////////////////////////////////

xml_node *XML_Copy_XML_Graph(xml_node *root)
{
  xml_node *cpy_root;

  cpy_root = XML_Make_Node(root->name);
  XML_Copy_XML_Node(cpy_root,root);

  return(cpy_root);
}

//////////////////////////////////////////////////////////////

void XML_Copy_XML_Node(xml_node *cpy_root, xml_node *root)
{
  xml_attr *attr,*cpy_attr;
  
  strcpy(cpy_root->name,root->name);

  XML_Make_Node_Id(cpy_root,root->id);  
  if(root->id) strcpy(cpy_root->id,root->id);

  XML_Make_Node_Value(cpy_root,root->value);
  if(root->value) strcpy(cpy_root->value,root->value);

  cpy_root->n_attr = root->n_attr;

  if(root->attr)
    {
      cpy_root->attr = XML_Make_Attribute(NULL,root->attr->name,root->attr->value);
      XML_Init_Attribute(cpy_root->attr);
      attr           = root->attr;
      cpy_attr       = cpy_root->attr;
      while(attr->next)
        {
          fflush(NULL);
          cpy_attr->next = XML_Make_Attribute(cpy_attr,attr->next->name,attr->next->value);
          XML_Init_Attribute(cpy_attr->next);
          attr           = attr->next;
          cpy_attr       = cpy_attr->next;
        }
    }
   
  if(root->child)
    {
      cpy_root->child = XML_Make_Node(root->child->name);
      cpy_root->child->parent = cpy_root;
      XML_Copy_XML_Node(cpy_root->child,root->child);
    }

  if(root->next)
    {
      cpy_root->next = XML_Make_Node(root->next->name);
      cpy_root->next->prev = cpy_root;
      XML_Copy_XML_Node(cpy_root->next,root->next);
    }
}

//////////////////////////////////////////////////////////////

void XML_Write_XML_Graph(FILE *fp, xml_node *root)
{
  int indent;
  indent = 0;
  XML_Write_XML_Node(fp,&indent,root);
}

//////////////////////////////////////////////////////////////

void XML_Write_XML_Node(FILE *fp, int *indent, xml_node *root)
{
  xml_node *n;
  xml_attr *attr;
  char *s;
  int i;

  s = (char *)mCalloc((*indent)+1,sizeof(char));
  For(i,(*indent)) s[i]='\t';
  s[i]='\0';


  PhyML_Fprintf(fp,"\n%s",s);
  
  n = root;
  
  PhyML_Fprintf(fp,"<%s",n->name);

  attr = n->attr;  
  while(attr)
    {
      PhyML_Fprintf(fp," %s=\"%s\"",attr->name,attr->value);
      fflush(NULL);
      attr = attr->next;
    }

  
  if(n->child)
    {
      (*indent)++;
      PhyML_Fprintf(fp,">");
      XML_Write_XML_Node(fp,indent,n->child);
      PhyML_Fprintf(fp,"\n%s</%s>\n",s,n->name);
      (*indent)--;
    }
  else
    {
      PhyML_Fprintf(fp,"/>");
    }

  if(n->next)
    XML_Write_XML_Node(fp,indent,n->next);
  

  Free(s);
}

//////////////////////////////////////////////////////////////

void Check_Mandatory_XML_Node(xml_node *root, char *name)
{
  if(!XML_Search_Node_Name(name,NO,root))
    {
      PhyML_Printf("\n== Could not find mandatory XML node with name '%s'.",name);
      PhyML_Printf("\n== Please amend your XML file.");
      Exit("\n");
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
