/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include "io.h"
#include "assert.h"

#ifdef BEAGLE
#include "libhmsbeagle/beagle.h"
#endif

/* Tree parser function. We need to pass a pointer to the string of characters
   since this string might be freed and then re-allocated by that function
   (i.e., its address in memory might change)
*/
t_tree *Read_Tree(char **s_tree)
{
  char  **subs;
  int     i, n_ext, n_int, n_otu;
  t_tree *tree;
  int     degree, len;
  t_node *root_node;

  /* PhyML_Printf("\n. Reading tree %s",(*s_tree)); */

  n_int = n_ext = 0;

  n_otu = 0;
  for (i = 0; i < (int)strlen((*s_tree)); ++i)
  {
    if ((*s_tree)[i] == '[') do
        ++i;
      while ((*s_tree)[i] != ']'); // Skip labels
    if ((*s_tree)[i] == ',') n_otu++;
  }
  n_otu += 1;

  tree = Make_Tree_From_Scratch(n_otu, NULL);
  subs = Sub_Trees((*s_tree), &degree);
  Clean_Multifurcation(subs, degree, 3);

  if (degree == 2)
  {
    root_node      = tree->a_nodes[2 * n_otu - 2];
    root_node->num = 2 * n_otu - 2;
    tree->n_root   = root_node;
    n_int -= 1;
  }
  else
  {
    root_node      = tree->a_nodes[n_otu];
    root_node->num = n_otu;
    tree->n_root   = NULL;
  }

  if (degree > 3) /* Multifurcation at the root. Need to re-assemble the
                     subtrees since Clean_Multifurcation added sets of
                     parenthesis and the corresponding NULL edges */
  {
    degree = 3;
    Free((*s_tree));
    len = 0;
    for (i = 0; i < degree; i++) len += (strlen(subs[i]) + 1);
    len += 5;

    (*s_tree) = (char *)mCalloc(len, sizeof(char));

    (*s_tree)[0] = '(';
    (*s_tree)[1] = '\0';
    for (i = 0; i < degree; i++)
    {
      strcat((*s_tree), subs[i]);
      strcat((*s_tree), ",\0");
    }

    sprintf((*s_tree) + strlen((*s_tree)) - 1, "%s", ");\0");

    i = 0;
    while (subs[i] != NULL) Free(subs[i++]);
    Free(subs);

    subs = Sub_Trees((*s_tree), &degree);
  }

  root_node->tax = 0;

  tree->has_branch_lengths        = 0;
  tree->num_curr_branch_available = tree->n_otu;
  for (i = 0; i < degree; i++)
    R_rtree((*s_tree), subs[i], root_node, tree, &n_int, &n_ext);

  i = degree;
  while (subs[i] != NULL) Free(subs[i++]);
  Free(subs);

  if (tree->n_root)
  {
    if (tree->n_root->v[1]->tax == NO && tree->n_root->v[2]->tax == NO)
    {
      tree->e_root       = tree->a_edges[tree->num_curr_branch_available];
      tree->n_root->b[1] = tree->a_edges[tree->num_curr_branch_available + 1];
      tree->n_root->b[2] = tree->a_edges[tree->num_curr_branch_available + 2];
    }
    else
    {
      for (i = 0; i < tree->n_otu; ++i)
        if (tree->a_edges[i]->left == NULL && tree->a_edges[i]->rght == NULL)
          break;
      assert(i != tree->n_otu);
      tree->e_root = tree->a_edges[i];

      tree->n_root->b[1] = tree->a_edges[tree->num_curr_branch_available];
      tree->n_root->b[2] = tree->a_edges[tree->num_curr_branch_available + 1];
    }

    tree->n_root->v[2]->v[0] = tree->n_root->v[1];
    tree->n_root->v[1]->v[0] = tree->n_root->v[2];

    tree->n_root->b[1]->left = tree->n_root;
    tree->n_root->b[2]->left = tree->n_root;
    tree->n_root->b[1]->rght = tree->n_root->v[1];
    tree->n_root->b[2]->rght = tree->n_root->v[2];

    // Reading edge lengths
    subs = Sub_Trees((*s_tree), &degree);
    Read_Branch_Length(subs[0], (*s_tree), tree->n_root->b[1], tree);
    Read_Branch_Length(subs[1], (*s_tree), tree->n_root->b[2], tree);
    // Edge labels
    Read_Edge_Label(subs[0], (*s_tree), tree->n_root->b[1]);
    Read_Edge_Label(subs[1], (*s_tree), tree->n_root->b[2]);
    // Node labels
    Read_Node_Label(subs[0], (*s_tree), tree->n_root);
    Read_Node_Label(subs[1], (*s_tree), tree->n_root);

    Free(subs);

    Connect_One_Edge_To_Two_Nodes(tree->n_root->v[2], tree->n_root->v[1],
                                  tree->e_root, tree);

    tree->e_root->l->v = tree->n_root->b[2]->l->v + tree->n_root->b[1]->l->v;
    if (tree->e_root->l->v > 0.0)
      tree->n_root_pos = tree->n_root->b[2]->l->v / tree->e_root->l->v;
    else
      tree->n_root_pos = .5;

    Update_Ancestors(tree->n_root, tree->n_root->v[2], tree->n_root->b[2],
                     tree);
    Update_Ancestors(tree->n_root, tree->n_root->v[1], tree->n_root->b[1],
                     tree);
  }

  return tree;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* 'a' in t_node a stands for ancestor. 'd' stands for descendant */
void R_rtree(char *s_tree_a, char *s_tree_d, t_node *a, t_tree *tree,
             int *n_int, int *n_ext)
{
  int     i;
  t_node *d;
  int     n_otu = tree->n_otu;

  if (strstr(s_tree_a, " "))
  {
    PhyML_Fprintf(stderr, "\n. [%s]", s_tree_a);
    Generic_Exit(__FILE__, __LINE__, __FUNCTION__);
  }

  // Internal edge
  if (s_tree_d[0] == '(')
  {
    char  **subs;
    int     degree;
    t_edge *b;

    (*n_int) += 1;

    if ((*n_int + n_otu) == (2 * n_otu - 1))
    {
      PhyML_Fprintf(stderr, "\n. The number of internal nodes in the tree "
                            "exceeds the number of taxa minus one.");
      PhyML_Fprintf(
          stderr,
          "\n. There probably is a formating problem in the input tree.");
      Generic_Exit(__FILE__, __LINE__, __FUNCTION__);
    }

    d      = tree->a_nodes[n_otu + *n_int];
    d->num = n_otu + *n_int;
    d->tax = 0;
    b      = tree->a_edges[tree->num_curr_branch_available];

    Read_Edge_Label(s_tree_d, s_tree_a, b);
    Read_Branch_Support(s_tree_d, s_tree_a, b, tree);
    Read_Branch_Length(s_tree_d, s_tree_a, b, tree);
    Read_Node_Label(s_tree_d, s_tree_a, d);

    if (tree->n_root && a == tree->n_root)
    {
      if (!a->v[1])
        a->v[1] = d;
      else
        a->v[2] = d;
    }
    else
    {
      for (i = 0; i < 3; i++)
      {
        if (!a->v[i])
        {
          a->v[i] = d;
          break;
        }
      }
    }
    d->v[0] = a;

    if (!(tree->n_root && a == tree->n_root))
      Connect_One_Edge_To_Two_Nodes(
          a, d, tree->a_edges[tree->num_curr_branch_available], tree);

    subs = Sub_Trees(s_tree_d, &degree);

    if (degree < 2)
    {
      PhyML_Fprintf(stderr,
                    "\n. A problem was detected in the following subtree:");
      PhyML_Fprintf(stderr, "\n. %s", s_tree_d);
      Generic_Exit(__FILE__, __LINE__, __FUNCTION__);
    }

    if (degree > 2)
    {
      Clean_Multifurcation(subs, degree, 2);

      Free(s_tree_d);

      s_tree_d =
          (char *)mCalloc(strlen(subs[0]) + strlen(subs[1]) + 5, sizeof(char));

      i = 0;
      while (subs[i] != NULL) Free(subs[i++]);
      Free(subs);

      subs = Sub_Trees(s_tree_d, &degree);
    }

    R_rtree(s_tree_d, subs[0], d, tree, n_int, n_ext);
    R_rtree(s_tree_d, subs[1], d, tree, n_int, n_ext);

    i = 2;
    while (subs[i] != NULL) Free(subs[i++]);
    Free(subs);
  }

  // External edge
  else
  {
    int i;

    d      = tree->a_nodes[*n_ext];
    d->tax = 1;

    Read_Edge_Label(s_tree_d, s_tree_a, tree->a_edges[*n_ext]);
    Read_Node_Name(s_tree_d, s_tree_a, d, tree);
    Read_Branch_Length(s_tree_d, s_tree_a, tree->a_edges[*n_ext], tree);
    Read_Node_Label(s_tree_d, s_tree_a, d);

    if (tree->n_root && a == tree->n_root)
    {
      if (!a->v[1])
        a->v[1] = d;
      else
        a->v[2] = d;
    }
    else
    {
      for (i = 0; i < 3; i++)
      {
        if (!a->v[i])
        {
          a->v[i] = d;
          break;
        }
      }
    }
    d->v[0] = a;

    if (!(tree->n_root && a == tree->n_root))
    {
      Connect_One_Edge_To_Two_Nodes(a, d, tree->a_edges[*n_ext], tree);
    }

    d->num = *n_ext;

    (*n_ext) += 1;
  }

  Free(s_tree_d);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Read_Node_Name(char *s_d, char *s_a, t_node *d, t_tree *tree)
{
  char *p;
  int   i;

  assert(s_d[0] != '(');

  d->name = (char *)mCalloc(strlen(s_d) + 1, sizeof(char));

  p = strstr(s_a, s_d);
  i = 0;
  while (p[i] != ':' && p[i] != '[')
  {
    d->name[i] = p[i];
    i++;
  }
  d->name[i]  = '\0';
  d->ori_name = d->name;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Read_Branch_Support(char *s_d, char *s_a, t_edge *b, t_tree *tree)
{
  char *sub_tp;
  char *p;
  int   i;

  if (s_d[0] != '(') return; // b is an external edge -> no edge support on it

  sub_tp = (char *)mCalloc(10 + strlen(s_d) + 1, sizeof(char));

  sub_tp[0] = '(';
  sub_tp[1] = '\0';
  strcat(sub_tp, s_d);
  p = strstr(s_a, sub_tp);

  if (!p)
  {
    sub_tp[0] = ',';
    sub_tp[1] = '\0';
    strcat(sub_tp, s_d);
    p = strstr(s_a, sub_tp);
  }
  assert(p);

  p++;

  i = 0;
  if (p[i] == '(')
  {
    i = Next_Matching_Char(p, '(', ')', i);
    i++;
  }
  if (p[i] == '[')
  {
    i = Next_Matching_Char(p, '[', ']', i);
    i++;
  }

  if (p[i] != ':')
  {
    b->support_val = atof((char *)(p + i));
  }

  Free(sub_tp);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Read_Branch_Length(char *s_d, char *s_a, t_edge *b, t_tree *tree)
{
  char *sub_tp;
  char *p;
  int   i;

  b->l->v = -1.;

  sub_tp = (char *)mCalloc(10 + strlen(s_d) + 1, sizeof(char));

  sub_tp[0] = '(';
  sub_tp[1] = '\0';
  strcat(sub_tp, s_d);
  p = strstr(s_a, sub_tp);

  if (!p)
  {
    sub_tp[0] = ',';
    sub_tp[1] = '\0';
    strcat(sub_tp, s_d);
    p = strstr(s_a, sub_tp);
  }
  assert(p);

  p++;

  i = 0;
  while (p[i] != ':' && p[i] != '\0')
  {
    if (p[i] == '(') i = Next_Matching_Char(p, '(', ')', i);
    if (p[i] == '[') i = Next_Matching_Char(p, '[', ']', i);
    i++;
  }
  if (p[i] == ':')
  {
    i++;
    if (p[i] == '[')
    {
      i = Next_Matching_Char(p, '[', ']', i);
      i++;
    }
    b->l->v                  = atof(p + i);
    tree->has_branch_lengths = YES;
    b->does_exist            = YES;
    /* PhyML_Printf("\n. READ LENGTH for s_d: %s b: %f",s_d,b->l->v); */
  }

  Free(sub_tp);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Read_Edge_Label(char *s_d, char *s_a, t_edge *b)
{
  char *sub_tp;
  char *p;
  int   i;

  sub_tp = (char *)mCalloc(10 + strlen(s_d) + 1, sizeof(char));

  sub_tp[0] = '(';
  sub_tp[1] = '\0';
  strcat(sub_tp, s_d);
  p = strstr(s_a, sub_tp);

  if (!p)
  {
    sub_tp[0] = ',';
    sub_tp[1] = '\0';
    strcat(sub_tp, s_d);
    p = strstr(s_a, sub_tp);
  }
  assert(p);

  p++;

  i = 0;
  while (p[i] != ':' && p[i] != '\0')
  {
    if (p[i] == '(') i = Next_Matching_Char(p, '(', ')', i);
    if (p[i] == '[') i = Next_Matching_Char(p, '[', ']', i);
    i++;
  }
  if (p[i] == ':')
  {
    i++;
    if (p[i] == '[')
    {
      char *s_lab;

      s_lab    = (char *)mCalloc(strlen(s_a) + 1, sizeof(char));
      s_lab[0] = '[';

      if (sscanf(p + i, "[%[^]]]", s_lab + 1) != 1)
      {
        PhyML_Fprintf(
            stderr, "\n. Edge label is in wrong format. A proper label should");
        PhyML_Fprintf(stderr,
                      "\n. look as follows: \"[xxx={yyy},xxxx={yy},...]\"");
        assert(FALSE);
      }

      s_lab[strlen(s_lab)] = ']';
      s_lab[strlen(s_lab)] = '\0';

      b->label = Read_Labels(s_lab);

      /* PhyML_Printf("\n. READ EDGE LABEL for s_d: %s lab.keyval:
       * %s:%s",s_d,b->label->key,b->label->val); */
      Free(s_lab);
      i++;
    }
  }

  Free(sub_tp);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Read_Node_Label(char *s_d, char *s_a, t_node *n)
{
  char *sub_tp;
  char *p;
  int   i;

  sub_tp = (char *)mCalloc(10 + strlen(s_d) + 1, sizeof(char));

  sub_tp[0] = '(';
  sub_tp[1] = '\0';
  strcat(sub_tp, s_d);
  p = strstr(s_a, sub_tp);

  if (!p)
  {
    sub_tp[0] = ',';
    sub_tp[1] = '\0';
    strcat(sub_tp, s_d);
    p = strstr(s_a, sub_tp);
  }
  assert(p);

  p++;

  i = 0;
  if (p[0] == '(') i = Next_Matching_Char(p, '(', ')', i);
  while (p[i] != '[' && p[i] != '\0') i++;

  if (p[i] == '[')
  {

    char *s_lab;

    s_lab    = (char *)mCalloc(strlen(s_a) + 1, sizeof(char));
    s_lab[0] = '[';

    if (sscanf(p + i, "[%[^]]]", s_lab + 1) != 1)
    {
      PhyML_Fprintf(stderr,
                    "\n. Node label is in wrong format. A proper label should");
      PhyML_Fprintf(stderr,
                    "\n. look as follows: \"[xxx={yyy},xxxx={yy},...]\"");
      assert(FALSE);
    }

    s_lab[strlen(s_lab)] = ']';
    s_lab[strlen(s_lab)] = '\0';

    n->label = Read_Labels(s_lab);

    /* PhyML_Printf("\n. READ NODE LABEL for s_d: %s %s lab.keyval: %s:%s
     * %s:%s", */
    /*              s_d, */
    /*              s_lab, */
    /*              n->label->key,n->label->val, */
    /*              n->label->next->key,n->label->next->val); */
    Free(s_lab);
  }
  Free(sub_tp);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Clean_Multifurcation(char **subtrees, int current_deg, int end_deg)
{

  if (current_deg <= end_deg)
    return;
  else
  {
    char *s_tmp;
    int   i;

    /* s_tmp = (char *)mCalloc(T_MAX_LINE,sizeof(char)); */
    s_tmp = (char *)mCalloc(10 + (int)strlen(subtrees[0]) + 1 +
                                (int)strlen(subtrees[1]) + 1,
                            sizeof(char));

    strcat(s_tmp, "(\0");
    strcat(s_tmp, subtrees[0]);
    strcat(s_tmp, ",\0");
    strcat(s_tmp, subtrees[1]);
    strcat(
        s_tmp,
        ")#NULL\0"); /* Add the label 'NULL' to identify a non-existing edge */
    Free(subtrees[0]);
    subtrees[0] = s_tmp;

    for (i = 1; i < current_deg - 1; i++) strcpy(subtrees[i], subtrees[i + 1]);

    Clean_Multifurcation(subtrees, current_deg - 1, end_deg);
  }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

char **Sub_Trees(char *tree, int *degree)
{
  char **subs;
  int    posbeg, posend;
  int    i;

  subs = NULL;

  if (tree[0] != '(')
  {
    *degree = 1;
    return NULL;
  }

  posbeg = posend = 1;
  (*degree)       = 0;
  do
  {
    posbeg = posend;
    while (tree[posend] != ',' && tree[posend] != ')')
    {
      if (tree[posend] == '(')
        posend = Next_Matching_Char(tree, '(', ')', posend);
      if (tree[posend] == '[')
        posend = Next_Matching_Char(tree, '[', ']', posend);
      posend++;
    }

    posend -= 1;

    if (*degree == 0)
      subs = (char **)mCalloc(1, sizeof(char *));
    else
      subs = (char **)mRealloc(subs, *degree + 1, sizeof(char *));

    subs[(*degree)] = (char *)mCalloc(strlen(tree) + 1, sizeof(char));
    strncpy(subs[(*degree)], tree + posbeg, posend - posbeg + 1);
    subs[(*degree)][posend - posbeg + 1] =
        '\0'; /* Thanks to Jean-Baka Domelevo-Entfellner */

    posend += 2;

    (*degree)++;
    if ((*degree) == NODE_DEG_MAX)
    {
      for (i = 0; i < (*degree); ++i)
        PhyML_Fprintf(stderr, "\n. Subtree %d : %s\n", i + 1, subs[i]);
      PhyML_Fprintf(stderr,
                    "\n. The degree of a t_node cannot be greater than %d\n",
                    NODE_DEG_MAX);
      Warn_And_Exit("\n");
    }
  } while (tree[posend - 1] != ')');

  subs            = (char **)mRealloc(subs, *degree + 1, sizeof(char *));
  subs[(*degree)] = NULL;

  return subs;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int Next_Matching_Char(char *s, char o, char c, int pos)
{
  int curr;

  curr = pos + 1;

  while (*(s + curr) != c)
  {
    if (*(s + curr) == o) curr = Next_Matching_Char(s, o, c, curr);
    curr++;
  }

  return curr;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Print_Tree(FILE *fp, t_tree *tree)
{
  char *s_tree;
  int   i;

  s_tree = (char *)Write_Tree(tree);

  if (OUTPUT_TREE_FORMAT == NEWICK)
    PhyML_Fprintf(fp, "%s\n", s_tree);
  else if (OUTPUT_TREE_FORMAT == NEXUS)
  {
    PhyML_Fprintf(fp, "#NEXUS\n");
    PhyML_Fprintf(fp, "BEGIN TREES;\n");
    PhyML_Fprintf(fp, "\tTRANSLATE\n");
    for (i = 0; i < tree->n_otu; ++i)
      PhyML_Fprintf(fp, "\t%3d\t%s,\n", i + 1, tree->a_nodes[i]->name);
    PhyML_Fprintf(fp, "\tUTREE PAUP_1=\n");
    PhyML_Fprintf(fp, "%s\n", s_tree);
    PhyML_Fprintf(fp, "ENDBLOCK;");
  }
  Free(s_tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

char *Write_Tree(t_tree *tree)
{
  char *s;
  int   i, available;

#ifndef MPI
  int init_len;
  init_len  = 3 * (int)T_MAX_NAME;
  s         = (char *)mCalloc(init_len, sizeof(char));
  available = init_len;
#elif defined MPI
  s = (char *)mCalloc(T_MAX_LINE, sizeof(char));
#endif

  i    = -1;
  s[0] = '(';

  if (tree->n_root == NULL)
  {
    i = 0;
    while ((!tree->a_nodes[tree->n_otu + i]->v[0]) ||
           (!tree->a_nodes[tree->n_otu + i]->v[1]) ||
           (!tree->a_nodes[tree->n_otu + i]->v[2]))
      i++;

    R_wtree(tree->a_nodes[tree->n_otu + i],
            tree->a_nodes[tree->n_otu + i]->v[0],
            tree->a_nodes[tree->n_otu + i]->b[0], &available, &s, tree);
    R_wtree(tree->a_nodes[tree->n_otu + i],
            tree->a_nodes[tree->n_otu + i]->v[1],
            tree->a_nodes[tree->n_otu + i]->b[1], &available, &s, tree);
    R_wtree(tree->a_nodes[tree->n_otu + i],
            tree->a_nodes[tree->n_otu + i]->v[2],
            tree->a_nodes[tree->n_otu + i]->b[2], &available, &s, tree);
  }
  else
  {
    R_wtree(tree->n_root, tree->n_root->v[1], tree->n_root->b[1], &available,
            &s, tree);
    R_wtree(tree->n_root, tree->n_root->v[2], tree->n_root->b[2], &available,
            &s, tree);
  }

  s[(int)strlen(s) - 1] = ')';

  if (tree->n_root != NULL)
    if (tree->print_labels == YES)
      Print_Labels(NULL, s + (int)strlen(s), tree->n_root->label);

  if (tree->io && tree->io->print_node_num == YES)
  {
    if (!tree->n_root)
      sprintf(s + (int)strlen(s), "%d", tree->a_nodes[tree->n_otu + i]->num);
    else
      sprintf(s + (int)strlen(s), "%d", tree->n_root->num);
  }
  s[(int)strlen(s)] = ';';

  return s;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void R_wtree(t_node *pere, t_node *fils, t_edge *b, int *available,
             char **s_tree, t_tree *tree)
{
  int    i, p;
  char  *format;
  int    last_len;
  phydbl mean_len;

  format = (char *)mCalloc(100, sizeof(char));

  sprintf(format, "%%.%df", tree->bl_ndigits);

  p = -1;
  if (fils->tax == YES) // Tip node
  {
    last_len = (int)strlen(*s_tree);

    if (OUTPUT_TREE_FORMAT == NEWICK)
    {
      if (tree->write_tax_names == YES)
      {
        if (tree->io && tree->io->long_tax_names)
        {
          strcat(*s_tree, tree->io->long_tax_names[fils->num]);
        }
        else
        {
          if (fils->name && strlen(fils->name) > 0)
            strcat(*s_tree, fils->name);
          else
            sprintf(*s_tree + (int)strlen(*s_tree), "%d", fils->num + 1);
        }
      }
      else if (tree->write_tax_names == NO)
      {
        sprintf(*s_tree + (int)strlen(*s_tree), "%d", fils->num + 1);
      }
    }
    else if (OUTPUT_TREE_FORMAT == NEXUS)
    {
      sprintf(*s_tree + (int)strlen(*s_tree), "%d", fils->num + 1);
    }
    else
    {
      PhyML_Printf("\n. Unknown tree format.");
      PhyML_Printf("\n. Err. in file %s at line %d\n", __FILE__, __LINE__);
      PhyML_Printf("\n. s=%s\n", *s_tree);
    }

    if ((fils->b) && (fils->b[0]) && (tree->write_br_lens == YES))
    {
      if (tree->print_labels == YES)
        Print_Labels(NULL, *s_tree + (int)strlen(*s_tree), fils->label);

      (*s_tree)[(int)strlen(*s_tree)] = ':';

      if (tree->print_labels == YES)
        Print_Labels(NULL, *s_tree + (int)strlen(*s_tree), b->label);

      if (tree->is_mixt_tree == NO)
        mean_len = b->l->v;
      else
        mean_len = MIXT_Get_Mean_Edge_Len(b, tree);
      sprintf(*s_tree + (int)strlen(*s_tree), format, MAX(0.0, mean_len));
    }

    (*s_tree)[(int)strlen(*s_tree)] = ',';

#ifndef MPI
    (*available) -= ((int)strlen(*s_tree) - last_len);

    /* printf("\n0 Available = %d [%d
     * %d]",(*available),(int)strlen(*s_tree),last_len); */
    /* printf("\n0 %s [%d,%d]",*s_tree,(int)(int)strlen(*s_tree),*available); */

    if (*available < 0)
    {
      PhyML_Fprintf(stderr, "\n. s=%s\n", *s_tree);
      PhyML_Fprintf(stderr, "\n. len=%d\n", (int)strlen(*s_tree));
      PhyML_Fprintf(
          stderr,
          "\n. The sequence names in your input file might be too long.");
      PhyML_Fprintf(stderr, "\n. Err. in file %s at line %d\n", __FILE__,
                    __LINE__);
      Warn_And_Exit("");
    }

    if (*available < (int)S_TREE_CHUNK)
    {
      (*s_tree) = (char *)mRealloc(
          *s_tree, (int)strlen(*s_tree) + 3 * (int)S_TREE_CHUNK, sizeof(char));
      for (i = 0; i < 3 * (int)S_TREE_CHUNK; ++i)
        (*s_tree)[(int)strlen(*s_tree) + i] = '\0';
      (*available) = 3 * (int)S_TREE_CHUNK;
    }
#endif
  }
  else // Internal node
  {
    (*s_tree)[(int)strlen(*s_tree)] = '(';

#ifndef MPI
    (*available) -= 1;

    if (*available < (int)S_TREE_CHUNK)
    {
      (*s_tree) = (char *)mRealloc(
          *s_tree, (int)strlen(*s_tree) + 3 * (int)S_TREE_CHUNK, sizeof(char));
      for (i = 0; i < 3 * (int)S_TREE_CHUNK; ++i)
        (*s_tree)[(int)strlen(*s_tree) + i] = '\0';
      (*available) = 3 * (int)S_TREE_CHUNK;
    }
#endif

    for (i = 0; i < 3; i++)
    {
      if ((fils->v[i] != pere) && (fils->b[i] != tree->e_root))
        R_wtree(fils, fils->v[i], fils->b[i], available, s_tree, tree);
      else
        p = i;
    }

    if (p < 0) assert(false);

    last_len = (int)strlen(*s_tree);

    (*s_tree)[last_len - 1] = ')';

    if ((fils->b) && (tree->write_br_lens == YES))
    {
      if (tree->print_labels == YES)
        Print_Labels(NULL, *s_tree + (int)strlen(*s_tree), fils->label);

      if (tree->io && tree->io->print_support_val == YES)
      {
        if (tree->io->do_boot == YES || tree->io->do_bayesboot == YES)
        {
          sprintf(*s_tree + (int)strlen(*s_tree), "%.0f",
                  fils->b[p]->support_val);
        }
        else
        {
          sprintf(*s_tree + (int)strlen(*s_tree), "%f",
                  fils->b[p]->support_val);
        }
      }
      if (tree->io && tree->io->print_node_num == YES)
      {
        sprintf(*s_tree + (int)strlen(*s_tree), "%d", fils->num);
      }

      fflush(NULL);

      (*s_tree)[(int)strlen(*s_tree)] = ':';

      if (tree->print_labels == YES)
        Print_Labels(NULL, *s_tree + (int)strlen(*s_tree), b->label);

      if (tree->is_mixt_tree == NO)
        mean_len = b->l->v;
      else
        mean_len = MIXT_Get_Mean_Edge_Len(b, tree);
      sprintf(*s_tree + (int)strlen(*s_tree), format, MAX(0.0, mean_len));
    }

    (*s_tree)[(int)strlen(*s_tree)] = ',';

#ifndef MPI
    (*available) -= ((int)strlen(*s_tree) - last_len);

    if (*available < 0)
    {
      PhyML_Fprintf(stderr, "\n. available = %d", *available);
      PhyML_Fprintf(stderr, "\n. Err. in file %s at line %d\n", __FILE__,
                    __LINE__);
      Warn_And_Exit("");
    }

    if (*available < (int)S_TREE_CHUNK)
    {
      (*s_tree) = (char *)mRealloc(
          *s_tree, (int)strlen(*s_tree) + 3 * (int)S_TREE_CHUNK, sizeof(char));
      for (i = 0; i < 3 * (int)S_TREE_CHUNK; ++i)
        (*s_tree)[(int)strlen(*s_tree) + i] = '\0';
      (*available) = 3 * (int)S_TREE_CHUNK;
    }
#endif
  }

  Free(format);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Detect_Align_File_Format(option *io)
{
  int    c;
  fpos_t curr_pos;

  fgetpos(io->fp_in_align, &curr_pos);

  errno = 0;

  while ((c = fgetc(io->fp_in_align)) != EOF)
  {
    if (errno)
      io->data_file_format = PHYLIP;
    else if (c == '#')
    {
      char s[10], t[6] = "NEXUS";
      if (!fgets(s, 6, io->fp_in_align))
      {
        PhyML_Fprintf(stderr, "\n. Err. in file %s at line %d\n", __FILE__,
                      __LINE__);
        Exit("\n");
      }
      if (!strcmp(t, s))
      {
        fsetpos(io->fp_in_align, &curr_pos);
        io->data_file_format = NEXUS;
        return;
      }
    }
  }
  fsetpos(io->fp_in_align, &curr_pos);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Detect_Tree_File_Format(option *io)
{
  int    c;
  fpos_t curr_pos;

  assert(io->fp_in_tree != NULL);

  fgetpos(io->fp_in_tree, &curr_pos);

  errno = 0;

  while ((c = fgetc(io->fp_in_tree)) != EOF)
  {
    if (errno)
    {
      io->tree_file_format = PHYLIP;
      PhyML_Printf("\n. Detected PHYLIP tree file format.");
    }
    else if (c == '#')
    {
      char s[10], t[6] = "NEXUS";
      if (!fgets(s, 6, io->fp_in_tree))
      {
        PhyML_Fprintf(stderr, "\n. Err. in file %s at line %d\n", __FILE__,
                      __LINE__);
        Warn_And_Exit("");
      }
      if (!strcmp(t, s))
      {
        fsetpos(io->fp_in_tree, &curr_pos);
        io->tree_file_format = NEXUS;
        PhyML_Printf("\n. Detected NEXUS tree file format.");
        return;
      }
    }
  }

  fsetpos(io->fp_in_tree, &curr_pos);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

align **Get_Seq(option *io)
{
  io->data = NULL;

  if (!io->fp_in_align)
  {
    PhyML_Fprintf(stderr, "\n. Filehandle to '%s' seems to be closed.",
                  io->in_align_file);
    Exit("\n");
  }

  Detect_Align_File_Format(io);

  switch (io->data_file_format)
  {
  case PHYLIP:
  {
    io->data = Get_Seq_Phylip(io);
    break;
  }
  case NEXUS:
  {
    io->nex_com_list = Make_Nexus_Com();
    Init_Nexus_Format(io->nex_com_list, io->fp_in_align);
    Get_Nexus_Data(io);
    Free_Nexus(io);
    break;
  }
  default:
  {
    PhyML_Fprintf(stderr, "\n. Err. in file %s at line %d (function '%s')\n",
                  __FILE__, __LINE__, __FUNCTION__);
    Exit("\n");
    break;
  }
  }

  if (!io->data)
  {
    PhyML_Fprintf(stderr, "\n. Err. in file %s at line %d (function '%s')\n",
                  __FILE__, __LINE__, __FUNCTION__);
    Exit("\n");
  }
  else
  {
    Post_Process_Data(io);
  }

  if (io->n_otu < 3)
  {
    PhyML_Fprintf(
        stderr,
        "\n. PhyML needs at least three sequences to perform an analysis.");
    assert(FALSE);
  }

  return io->data;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Post_Process_Data(option *io)
{
  int    i, j, swap;
  align *data_buff;

  for (i = 0; i < io->data[0]->len; ++i)
  {
    for (j = 0; j < io->n_otu; ++j)
    {
      if ((io->data[j]->state[i] == '*') || (io->data[j]->state[i] == '?') ||
          (io->data[j]->state[i] == '-'))
        io->data[j]->state[i] = 'X';
      if ((io->datatype == NT) && (io->data[j]->state[i] == 'N'))
        io->data[j]->state[i] = 'X';
      if (io->data[j]->state[i] == 'U') io->data[j]->state[i] = 'T';
    }
  }

  for (i = 0; i < io->n_otu; ++i) io->data[i]->len = io->data[0]->len;

  /* Sequences are ordered alphabetically */
  data_buff = NULL;
  swap      = TRUE;
  /* swap = FALSE; */
  while (swap == TRUE)
  {
    swap = FALSE;
    for (i = 0; i < io->n_otu - 1; i++)
    {
      for (j = i + 1; j < io->n_otu; j++)
      {
        if (strcmp(io->data[i]->name, io->data[j]->name) < 0)
        {
          swap        = TRUE;
          data_buff   = io->data[i];
          io->data[i] = io->data[j];
          io->data[j] = data_buff;
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Get_Nexus_Data(option *io)
{
  char    *token;
  nexcom  *curr_com;
  nexparm *curr_parm;
  int      nxt_token_t, cur_token_t;
  FILE    *fp;

  assert(io->nex_com_list[0] != NULL);
  fp = io->nex_com_list[0]->fp;

  token = (char *)mCalloc(T_MAX_TOKEN, sizeof(char));

  curr_com    = NULL;
  curr_parm   = NULL;
  nxt_token_t = NEXUS_COM;
  cur_token_t = -1;

  do
  {
    /* PhyML_Printf("\n+ Token: '%s' next_token=%d
     * cur_token=%d",token,nxt_token_t,cur_token_t); */
    if (!Get_Token(fp, token)) break;

    if (token[0] == ';')
    {
      curr_com    = NULL;
      curr_parm   = NULL;
      nxt_token_t = NEXUS_COM;
      cur_token_t = -1;
    }

    if (nxt_token_t == NEXUS_EQUAL)
    {
      cur_token_t = NEXUS_VALUE;
      nxt_token_t = NEXUS_PARM;
      continue;
    }

    if ((nxt_token_t == NEXUS_COM) && (cur_token_t != NEXUS_VALUE))
    {
      Find_Nexus_Com(token, &curr_com, &curr_parm, io->nex_com_list);
      if (curr_com)
      {
        nxt_token_t = curr_com->nxt_token_t;
        cur_token_t = curr_com->cur_token_t;
      }
      if (cur_token_t != NEXUS_VALUE) continue;
    }

    if ((nxt_token_t == NEXUS_PARM) && (cur_token_t != NEXUS_VALUE))
    {
      Find_Nexus_Parm(token, &curr_parm, curr_com);
      if (curr_parm)
      {
        nxt_token_t = curr_parm->nxt_token_t;
        cur_token_t = curr_parm->cur_token_t;
      }
      if (cur_token_t != NEXUS_VALUE) continue;
    }

    if (cur_token_t == NEXUS_VALUE)
    {
      if ((curr_parm->func)(token, curr_parm, io)) /* Read in parameter value */
      {
        nxt_token_t = NEXUS_PARM;
        cur_token_t = -1;
      }
    }
  } while (strlen(token) > 0);

  Free(token);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int Get_Token(FILE *fp, char *token)
{
  char c;

  assert(fp != NULL);

  c = ' ';

  while (c == ' ' || c == '\t' || c == '\n' || c == '\r')
  {
    c = fgetc(fp);
    if (c == EOF) return 0;
  }

  if (c == '"')
  {
    do
    {
      *token = c;
      token++;
      c = fgetc(fp);
      if (c == EOF) return 0;
    } while (c != '"');
    *token = c;
    /* c = fgetc(fp); */
    if (c == EOF) return 0;
    *(token + 1) = '\0';
    return 1;
  }

  if (c == '[')
  {
    Skip_Comment(fp);
    c      = fgetc(fp);
    *token = c;
    token++;
    if (c == EOF) return 0;
    return 1;
  }

  if (c == '#')
  {
    *token = c;
    token++;
  }
  else if (c == ';')
  {
    *token = c;
    token++;
  }
  else if (c == ',')
  {
    *token = c;
    token++;
  }
  else if (c == '.')
  {
    *token = c;
    token++;
  }
  else if (c == '=')
  {
    *token = c;
    token++;
  }
  else if (c == '(')
  {
    *token = c;
    token++;
  }
  else if (c == ')')
  {
    *token = c;
    token++;
  }
  else if (c == '{')
  {
    *token = c;
    token++;
  }
  else if (c == '}')
  {
    *token = c;
    token++;
  }
  else if (c == '?')
  {
    *token = c;
    token++;
  }
  else if (c == '-')
  {
    *token = c;
    token++;
  }
  else
  {
    /* while(isgraph(c) && c != ';' && c != '-' && c != ',' && c != '=') */
    while (isgraph(c) && c != ';' && c != ',' && c != '=')
    {
      *(token++) = c;
      c          = fgetc(fp);
      if (c == EOF) return 0;
    }

    fseek(fp, -1 * sizeof(char), SEEK_CUR);
  }
  *token = '\0';
  return 1;
}

/* void Get_Token(char *line, char *token) */
/* { */
/*   while(**line == ' ' || **line == '\t') (*line)++; */

/*   if(**line == '"')  */
/*     { */
/*       do { *token = **line; (*line)++; token++; } while(**line != '"'); */
/*       *token = **line; */
/*       (*line)++; */
/*       *(token+1) = '\0'; */
/*       return; */
/*     } */

/*   if(**line == '[')  */
/*     { */
/*       int in_comment; */

/*       in_comment = 1; */
/*       do  */
/* 	{  */
/* 	  (*line)++;  */
/* 	  if(**line == '[')  */
/* 	    { */
/* 	      in_comment++; */
/* 	    } */
/* 	  else if(**line == ']') in_comment--;	   */
/* 	} */
/*       while(in_comment); */
/*       (*line)++; */
/*       return; */
/*     } */

/*   if(**line == '#')      {*token = **line; (*line)++; token++; } */
/*   else if(**line == ';') {*token = **line; (*line)++; token++; } */
/*   else if(**line == ',') {*token = **line; (*line)++; token++; } */
/*   else if(**line == '.') {*token = **line; (*line)++; token++; } */
/*   else if(**line == '=') {*token = **line; (*line)++; token++; } */
/*   else if(**line == '(') {*token = **line; (*line)++; token++; } */
/*   else if(**line == ')') {*token = **line; (*line)++; token++; } */
/*   else if(**line == '{') {*token = **line; (*line)++; token++; } */
/*   else if(**line == '}') {*token = **line; (*line)++; token++; } */
/*   else if(**line == '?') {*token = **line; (*line)++; token++; } */
/*   else if(**line == '-') {*token = **line; (*line)++; token++; } */
/*   else */
/*     { */
/*       while(isgraph(**line) && **line != ';' && **line != '=' && **line !=
 * ',')  */
/* 	{ */
/* 	  *(token++) = **line; */
/* 	  (*line)++;  */
/* 	} */
/*     } */
/*   *token = '\0'; */
/* } */

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

align **Get_Seq_Phylip(option *io)
{
  Read_Ntax_Len_Phylip(io->fp_in_align, &io->n_otu, &io->init_len);

  if (io->n_otu > N_MAX_OTU)
  {
    PhyML_Fprintf(stderr, "\n. The number of taxa should not exceed %d",
                  N_MAX_OTU);
    assert(FALSE);
  }

  if (io->interleaved == YES)
    io->data = Read_Seq_Interleaved(io);
  else
    io->data = Read_Seq_Sequential(io);

  return io->data;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Read_Ntax_Len_Phylip(FILE *fp, int *n_otu, int *n_tax)
{
  char *line;
  int   readok;

  line = (char *)mCalloc(T_MAX_LINE, sizeof(char));

  readok = 0;
  do
  {
    if (fscanf(fp, "%s", line) == EOF)
    {
      Free(line);
      PhyML_Fprintf(stderr, "\n. PhyML can't read in this alignment.");
      PhyML_Fprintf(stderr, "\n. Could it be that sequence file is empty?");
      PhyML_Fprintf(stderr, "\n. Err. in file %s at line %d\n", __FILE__,
                    __LINE__);
      Exit("\n");
    }
    else
    {
      if (strcmp(line, "\n") && strcmp(line, "\r") && strcmp(line, "\t"))
      {
        sscanf(line, "%d", n_otu);
        if (*n_otu <= 0)
          Warn_And_Exit("\n. The number of taxa cannot be negative.\n");

        if (!fscanf(fp, "%s", line)) Exit("\n");
        sscanf(line, "%d", n_tax);
        if (*n_tax <= 0)
          Warn_And_Exit("\n. The sequence length cannot be negative.\n");
        else
          readok = 1;
      }
    }
  } while (!readok);

  Free(line);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

align **Read_Seq_Sequential(option *io)
{
  int     i;
  char   *line;
  align **data;
  /*   char c; */
  char *format;

  format = (char *)mCalloc(T_MAX_NAME, sizeof(char));
  line   = (char *)mCalloc(T_MAX_LINE, sizeof(char));
  data   = (align **)mCalloc(io->n_otu, sizeof(align *));

  /*   while((c=fgetc(in))!='\n'); */
  /*  while(((c=fgetc(io->fp_in_align))!='\n') && (c != ' ') && (c != '\r') &&
   * (c != '\t')); */

  sprintf(format, "%%%ds", T_MAX_NAME);

  for (i = 0; i < io->n_otu; i++)
  {
    data[i]       = (align *)mCalloc(1, sizeof(align));
    data[i]->name = (char *)mCalloc(T_MAX_NAME, sizeof(char));
    data[i]->state =
        (char *)mCalloc(io->init_len * io->state_len + 1, sizeof(char));

    data[i]->is_ambigu = NULL;
    data[i]->len       = 0;

    if (!fscanf(io->fp_in_align, format, data[i]->name)) Exit("\n");

    Check_Sequence_Name(data[i]->name);

    while (data[i]->len < io->init_len * io->state_len)
      assert(Read_One_Line_Seq(&data, i, io->fp_in_align));

    if (data[i]->len != io->init_len * io->state_len)
    {
      PhyML_Fprintf(
          stderr,
          "\n. Err. Problem with species %s's sequence (check the format).\n",
          data[i]->name);
      PhyML_Fprintf(stderr,
                    "\n. Observed sequence length: %d, expected length: %d\n",
                    data[i]->len, io->init_len * io->state_len);
      Warn_And_Exit("");
    }
  }

  for (i = 0; i < io->n_otu; i++) data[i]->state[data[i]->len] = '\0';

  Restrict_To_Coding_Position(data, io);

  Free(format);
  Free(line);

  return data;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

align **Read_Seq_Interleaved(option *io)
{
  int     i, end, num_block;
  char   *line;
  align **data;
  /*   char c; */
  char  *format;
  fpos_t curr_pos;

  line   = (char *)mCalloc(T_MAX_LINE, sizeof(char));
  format = (char *)mCalloc(T_MAX_NAME, sizeof(char));
  data   = (align **)mCalloc(io->n_otu, sizeof(align *));

  /*   while(((c=fgetc(io->fp_in_align))!='\n') && (c != ' ') && (c != '\r') &&
   * (c != '\t')); */

  sprintf(format, "%%%ds", T_MAX_NAME);

  end = 0;
  for (i = 0; i < io->n_otu; i++)
  {
    data[i]       = (align *)mCalloc(1, sizeof(align));
    data[i]->name = (char *)mCalloc(T_MAX_NAME, sizeof(char));
    data[i]->state =
        (char *)mCalloc(io->init_len * io->state_len + 1, sizeof(char));

    data[i]->len       = 0;
    data[i]->is_ambigu = NULL;

    if (!fscanf(io->fp_in_align, format, data[i]->name)) Exit("\n");

    Check_Sequence_Name(data[i]->name);

    if (!Read_One_Line_Seq(&data, i, io->fp_in_align))
    {
      end = 1;
      if ((i != io->n_otu) && (i != io->n_otu - 1))
      {
        PhyML_Fprintf(stderr, "\n. Err.: problem with species %s's sequence.\n",
                      data[i]->name);
        PhyML_Fprintf(stderr,
                      "\n. Observed sequence length: %d, expected length: %d\n",
                      data[i]->len, io->init_len * io->state_len);
        Exit("");
      }
      break;
    }
  }

  if (data[0]->len == io->init_len * io->state_len) end = 1;

  /*   if(end) printf("\n. finished yet '%c'\n",fgetc(io->fp_in_align)); */
  if (!end)
  {
    end       = 0;
    num_block = 1;
    do
    {
      num_block++;

      /* interblock */
      if (!fgets(line, T_MAX_LINE, io->fp_in_align)) break;

      if (line[0] != 13 && line[0] != 10)
      {
        PhyML_Fprintf(stderr,
                      "\n. Err.: one or more missing sequences in block %d.\n",
                      num_block - 1);
        Exit("");
      }

      for (i = 0; i < io->n_otu; i++)
        if (data[i]->len != io->init_len * io->state_len) break;

      if (i == io->n_otu) break;

      for (i = 0; i < io->n_otu; i++)
      {
        /* Skip the taxon name, if any, in this interleaved block */
        fgetpos(io->fp_in_align, &curr_pos);
        if (!fscanf(io->fp_in_align, format, line))
        {
          PhyML_Fprintf(stderr, "\n. Err. in file %s at line %d\n", __FILE__,
                        __LINE__);
          Warn_And_Exit("");
        }
        if (line && strcmp(line, data[i]->name))
          fsetpos(io->fp_in_align, &curr_pos);

        if (data[i]->len > io->init_len * io->state_len)
        {
          PhyML_Fprintf(stderr,
                        "\n. Observed sequence length=%d expected length=%d.\n",
                        data[i]->len, io->init_len * io->state_len);
          PhyML_Fprintf(stderr,
                        "\n. Err.: Problem with species %s's sequence.\n",
                        data[i]->name);
          Exit("");
        }
        else if (!Read_One_Line_Seq(&data, i, io->fp_in_align))
        {
          end = 1;
          if ((i != io->n_otu) && (i != io->n_otu - 1))
          {
            PhyML_Fprintf(stderr,
                          "\n. Err.: Problem with species %s's sequence.\n",
                          data[i]->name);
            PhyML_Fprintf(
                stderr,
                "\n. Observed sequence length: %d, expected length: %d.\n",
                data[i]->len, io->init_len * io->state_len);
            Exit("");
          }
          break;
        }
      }
    } while (!end);
  }

  for (i = 0; i < io->n_otu; i++) data[i]->state[data[i]->len] = '\0';

  for (i = 0; i < io->n_otu; i++)
  {
    if (data[i]->len != io->init_len * io->state_len)
    {
      PhyML_Fprintf(stderr,
                    "\n. Check sequence '%s' length (expected length: %d, "
                    "observed length: %d) [OTU %d].\n",
                    data[i]->name, io->init_len, data[i]->len, i + 1);
      Exit("");
    }
  }

  Restrict_To_Coding_Position(data, io);

  Free(format);
  Free(line);

  return data;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int Read_One_Line_Seq(align ***data, int num_otu, FILE *in)
{
  char c     = ' ';
  int  nchar = 0;

  while (1)
  {
    /*       if((c == EOF) || (c == '\n') || (c == '\r')) break; */

    if ((c == 13) || (c == 10))
    {
      /* PhyML_Printf("[%d %d]\n",c,nchar); fflush(NULL); */
      if (!nchar)
      {
        c = (char)fgetc(in);
        continue;
      }
      else
      {
        /* PhyML_Printf("break\n"); */
        break;
      }
    }
    else if (c == EOF)
    {
      /* PhyML_Printf("EOL\n"); */
      break;
    }
    else if ((c == ' ') || (c == '\t') || (c == 32))
    {
      /* PhyML_Printf("[%d]",c); */
      c = (char)fgetc(in);
      continue;
    }

    nchar++;
    Uppercase(&c);

    if (c == '.')
    {
      c = (*data)[0]->state[(*data)[num_otu]->len];
      if (!num_otu)
        Warn_And_Exit(
            "\n. Err: Symbol \".\" should not appear in the first sequence\n");
    }
    (*data)[num_otu]->state[(*data)[num_otu]->len] = c;
    (*data)[num_otu]->len++;
    c = (char)fgetc(in);
    /* PhyML_Printf("[%c %d]",c,c); fflush(NULL); */
    if (c == ';') break;
  }

  /* printf("\n. Exit nchar: %d [%d]\n",nchar,c==EOF); */
  if (c == EOF)
    return 0;
  else
    return 1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

scalar_dbl *Read_Io_Weights(option *io)
{
  scalar_dbl *w, *ori, *prev = NULL;
  double      val;

  assert(io->weight_file);

  io->fp_weight_file = Openfile(io->weight_file, READ);

  w   = (scalar_dbl *)mCalloc(1, sizeof(scalar_dbl));
  ori = w;

  do
  {
    if (fscanf(io->fp_weight_file, "%lf,", &val) == EOF) break;
    w->v    = (phydbl)val;
    w->next = (scalar_dbl *)mCalloc(1, sizeof(scalar_dbl));
    prev    = w;
    w       = w->next;
  } while (1);

  /* Remove the last allocated and empty element of the list */
  if (prev != NULL && prev->next != NULL)
  {
    Free(prev->next);
    prev->next = NULL;
  }

  fclose(io->fp_weight_file);

  return (ori);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

char *Return_Tree_String_Phylip(FILE *fp_input_tree)
{
  char *line;
  int   i;
  char  c;
  int   open, maxopen;

  if (fp_input_tree == NULL)
  {
    PhyML_Fprintf(stderr, "\n. Err. in file %s at line %d\n", __FILE__,
                  __LINE__);
    Warn_And_Exit("");
  }

  do
  {
    c = fgetc(fp_input_tree);
  } while ((c != '(') && (c != EOF));

  if (c == EOF) return NULL;

  line    = (char *)mCalloc(1, sizeof(char));
  open    = 1;
  maxopen = open;

  i = 0;
  for (;;)
  {
    if ((c == ' ') || (c == '\n'))
    {
      c = fgetc(fp_input_tree);
      if (c == EOF || c == ';')
        break;
      else
        continue;
    }

    /* if(c == '[') */
    /*   { */
    /*     Skip_Comment(fp_input_tree); */
    /*     c = fgetc(fp_input_tree); */
    /*     if(c == EOF || c == ';') break; */
    /*   } */

    line = (char *)mRealloc(line, i + 2, sizeof(char));

    line[i] = c;
    i++;
    c = fgetc(fp_input_tree);
    if (c == EOF || c == ';') break;
    if (c == '(') open++;
    if (c == ')') open--;
    if (open > maxopen) maxopen = open;
  }
  line[i] = '\0';

  /* if(maxopen == 1) return NULL; */
  return line;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

t_tree *Read_Tree_File_Phylip(FILE *fp_input_tree)
{
  char   *line;
  t_tree *tree;

  line = Return_Tree_String_Phylip(fp_input_tree);
  tree = Read_Tree(&line);
  Free(line);

  return tree;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Print_Site(calign *cdata, int num, int n_otu, char *sep, int stepsize,
                FILE *fp)
{
  int i, j;
  PhyML_Fprintf(fp, "\n");
  for (i = 0; i < n_otu; i++)
  {
    PhyML_Fprintf(fp, "%20s ", cdata->c_seq[i]->name);
    for (j = 0; j < stepsize; j++)
      PhyML_Fprintf(fp, "%c", cdata->c_seq[i]->state[num + j]);
    PhyML_Fprintf(fp, "%s", sep);
  }
  PhyML_Fprintf(fp, "%s", sep);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Print_Site_Lk(t_tree *tree, FILE *fp)
{
  int    site;
  int    catg;
  char  *s;
  phydbl postmean, sum;

  assert(fp);
  rewind(fp);

  if (tree->is_mixt_tree == YES)
  {
    MIXT_Print_Site_Lk(tree, fp);
    return;
  }

  assert(tree->io->print_site_lnl == YES);

  if (!tree->io->print_trace)
  {
    s = (char *)mCalloc(T_MAX_LINE, sizeof(char));

    PhyML_Fprintf(fp, "# Note : P(D|M) is the probability of site D given the "
                      "model M (i.e., the site likelihood)\n");
    if (tree->mod->ras->n_catg > 1 || tree->mod->ras->invar)
    {
      PhyML_Fprintf(fp, "# P*(D|M,rr[x]) is the scaled probability of site D "
                        "given the model M and the relative rate\n");
      PhyML_Fprintf(fp, "# of evolution rr[x], where x is the class of rate to "
                        "be considered.\n");
      PhyML_Fprintf(fp, "# The actual conditional probability is given by "
                        "P*(D|M,rr[x])/2^F, where\n");
      PhyML_Fprintf(fp, "# F is the scaling factor (see column 'Scaler').\n");
      PhyML_Fprintf(fp, "# For invariant sites, P(D|M,rr[0]=0) is the actual "
                        "conditional probability\n");
      PhyML_Fprintf(fp, "# (i.e., it is not scaled).\n");
    }

    PhyML_Fprintf(fp, "\n\n");

    sprintf(s, "Site");
    PhyML_Fprintf(fp, "%-12s", s);

    sprintf(s, "P(D|M)");
    PhyML_Fprintf(fp, "%-15s", s);

    sprintf(s, "Scaler");
    PhyML_Fprintf(fp, "%-7s", s);

    sprintf(s, "Pattern");
    PhyML_Fprintf(fp, "%-9s", s);

    if (tree->mod->ras->n_catg > 1)
    {
      for (catg = 0; catg < tree->mod->ras->n_catg; catg++)
      {
        sprintf(s, "P*(D|M,rr[%d]=%5.4f)", catg + 1,
                tree->mod->ras->gamma_rr->v[catg]);
        PhyML_Fprintf(fp, "%-23s", s);
      }

      sprintf(s, "Posterior mean");
      PhyML_Fprintf(fp, "%-22s", s);
    }

    if (tree->mod->ras->invar)
    {
      sprintf(s, "P(D|M,rr[0]=0)");
      PhyML_Fprintf(fp, "%-16s", s);
    }

    sprintf(s, "NDistinctStates");
    PhyML_Fprintf(fp, "%-16s", s);

    PhyML_Fprintf(fp, "\n");

    Init_Ui_Tips(tree);

    for (site = 0; site < tree->data->init_len; site++)
    {
      PhyML_Fprintf(fp, "%-12d", site + 1);

      PhyML_Fprintf(fp, "%-15g", tree->cur_site_lk[tree->data->sitepatt[site]]);

      PhyML_Fprintf(fp, "%-7d",
                    tree->fact_sum_scale[tree->data->sitepatt[site]]);

      PhyML_Fprintf(fp, "%-9d", tree->data->sitepatt[site]);

      if (tree->mod->ras->n_catg > 1)
      {
        for (catg = 0; catg < tree->mod->ras->n_catg; catg++)
        {
          PhyML_Fprintf(fp, "%-23g",
                        tree->unscaled_site_lk_cat[tree->data->sitepatt[site] *
                                                       tree->mod->ras->n_catg +
                                                   catg]);
        }

        postmean = .0;
        for (catg = 0; catg < tree->mod->ras->n_catg; catg++)
          postmean += tree->mod->ras->gamma_rr->v[catg] *
                      tree->unscaled_site_lk_cat[tree->data->sitepatt[site] *
                                                     tree->mod->ras->n_catg +
                                                 catg] *
                      tree->mod->ras->gamma_r_proba->v[catg];

        sum = .0;
        for (catg = 0; catg < tree->mod->ras->n_catg; catg++)
        {
          sum += tree->unscaled_site_lk_cat[tree->data->sitepatt[site] *
                                                tree->mod->ras->n_catg +
                                            catg] *
                 tree->mod->ras->gamma_r_proba->v[catg];
        }

        postmean /= sum;

        PhyML_Fprintf(fp, "%-22g", postmean);
      }

      if (tree->mod->ras->invar)
      {
        if ((phydbl)tree->data->invar[tree->data->sitepatt[site]] > -0.5)
          PhyML_Fprintf(fp, "%-16g",
                        tree->mod->e_frq->pi
                            ->v[tree->data->invar[tree->data->sitepatt[site]]]);
        else
          PhyML_Fprintf(fp, "%-16g", 0.0);
      }

      PhyML_Fprintf(
          fp, "%-16d",
          Number_Of_Diff_States_One_Site(tree->data->sitepatt[site], tree));

      PhyML_Fprintf(fp, "\n");
    }
    Free(s);
  }
  else
  {
    for (site = 0; site < tree->data->init_len; site++)
      PhyML_Fprintf(fp, "%.2f\t",
                    log(tree->cur_site_lk[tree->data->sitepatt[site]]));
    PhyML_Fprintf(fp, "\n");
  }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Print_Seq(FILE *fp, align **data, int n_otu)
{
  int i, j;

  PhyML_Fprintf(fp, "%d\t%d\n", n_otu, data[0]->len);
  for (i = 0; i < n_otu; i++)
  {
    PhyML_Fprintf(fp, "%s\t", data[i]->name);
    /* for(j=0;j<20;j++) */
    /*   { */
    /*     if(j<(int)strlen(data[i]->name)) */
    /*       fputc(data[i]->name[j],fp); */
    /*     else fputc(' ',fp); */
    /*   } */
    /*       PhyML_Printf("%10d  ",i); */
    For(j, data[i]->len) { PhyML_Fprintf(fp, "%c", data[i]->state[j]); }
    PhyML_Fprintf(fp, "\n");
  }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Print_CSeq(FILE *fp, int compressed, calign *cdata, t_tree *tree)
{
  int i, j;
  int n_otu;

  n_otu = cdata->n_otu;
  if (cdata->format == PHYLIP)
  {
    PhyML_Fprintf(fp, "%d\t%d\n", n_otu, cdata->init_len);
  }
  else if (cdata->format == NEXUS)
  {
    PhyML_Fprintf(fp, "#NEXUS\n");
    PhyML_Fprintf(fp, "begin data\n");
    PhyML_Fprintf(fp, "dimensions ntax=%d nchar=%d;\n", n_otu, cdata->init_len);
    PhyML_Fprintf(fp, "format sequential datatype=dna;\n");
    PhyML_Fprintf(fp, "matrix\n");
  }
  else
  {
    PhyML_Fprintf(fp,
                  "This sample file is to be processed by IBDSim "
                  "(http://www1.montpellier.inra.fr/CBGP/software/ibdsim/)");
  }

  if (cdata->format == PHYLIP || cdata->format == NEXUS)
  {
    for (i = 0; i < n_otu; i++)
    {
      for (j = 0; j < 50; j++)
      {
        if (j < (int)strlen(cdata->c_seq[i]->name))
          fputc(cdata->c_seq[i]->name[j], fp);
        else
          fputc(' ', fp);
      }

      if (compressed == YES) /* Print out compressed sequences */
        PhyML_Fprintf(fp, "%s", cdata->c_seq[i]->state);
      else /* Print out uncompressed sequences */
      {
        for (j = 0; j < cdata->init_len; j++)
        {
          PhyML_Fprintf(fp, "%c", cdata->c_seq[i]->state[cdata->sitepatt[j]]);
        }
      }
      PhyML_Fprintf(fp, "\n");
    }
    PhyML_Fprintf(fp, "\n");

    if (cdata->format == NEXUS)
    {
      PhyML_Fprintf(fp, ";\n");
      PhyML_Fprintf(fp, "END;\n");
    }
  }
  else if (cdata->format == IBDSIM)
  {
    for (i = 0; i < cdata->init_len; i++) PhyML_Fprintf(fp, "\nlocus %6d", i);
    for (i = 0; i < n_otu; i++)
    {
      PhyML_Fprintf(fp, "\npop");
      PhyML_Fprintf(fp, "%12f  %12f , ",
                    tree ? tree->a_nodes[i]->coord->lonlat[0] : -1.,
                    tree ? tree->a_nodes[i]->coord->lonlat[1] : -1.);

      if (tree != NULL)
      {
        for (j = 0; j < cdata->init_len; j++)
        {
          switch (tree->a_nodes[i]->c_seq->state[j])
          {
          case 'A':
          {
            PhyML_Fprintf(fp, "001 ");
            break;
          }
          case 'C':
          {
            PhyML_Fprintf(fp, "002 ");
            break;
          }
          case 'G':
          {
            PhyML_Fprintf(fp, "003 ");
            break;
          }
          case 'T':
          {
            PhyML_Fprintf(fp, "004 ");
            break;
          }
          }
        }
      }
    }
  }
  /*   PhyML_Printf("\t"); */
  /*   for(j=0;j<cdata->crunch_len;j++) */
  /*     PhyML_Printf("%.0f ",cdata->wght[j]); */
  /*   PhyML_Printf("\n"); */
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Print_CSeq_Select(FILE *fp, int compressed, calign *cdata, t_tree *tree)
{
  int    i, j;
  int    n_otu;
  phydbl eps;

  int slice = 14;

  eps   = 1.E-6;
  n_otu = 0;
  for (i = 0; i < cdata->n_otu; i++)
    if (tree->times->nd_t[i] < tree->times->time_slice_lims[slice] + eps)
      n_otu++;

  PhyML_Fprintf(fp, "%d\t%d\n", n_otu, cdata->init_len);

  n_otu = cdata->n_otu;

  for (i = 0; i < n_otu; i++)
  {
    if (tree->times->nd_t[i] < tree->times->time_slice_lims[slice] + eps)
    {
      for (j = 0; j < 50; j++)
      {
        if (j < (int)strlen(cdata->c_seq[i]->name))
          fputc(cdata->c_seq[i]->name[j], fp);
        else
          fputc(' ', fp);
      }

      if (compressed == YES) /* Print out compressed sequences */
        PhyML_Fprintf(fp, "%s", cdata->c_seq[i]->state);
      else /* Print out uncompressed sequences */
      {
        for (j = 0; j < cdata->init_len; j++)
        {
          PhyML_Fprintf(fp, "%c", cdata->c_seq[i]->state[cdata->sitepatt[j]]);
        }
      }
      PhyML_Fprintf(fp, "\n");
    }
  }

  if (cdata->format == 1)
  {
    PhyML_Fprintf(fp, ";\n");
    PhyML_Fprintf(fp, "END;\n");
  }

  /*   PhyML_Printf("\t"); */
  /*   for(j=0;j<cdata->crunch_len;j++) */
  /*     PhyML_Printf("%.0f ",cdata->wght[j]); */
  /*   PhyML_Printf("\n"); */
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Print_Dist(matrix *mat)
{
  int i, j;

  for (i = 0; i < mat->n_otu; i++)
  {
    PhyML_Printf("%s ", mat->name[i]);

    for (j = 0; j < mat->n_otu; j++) PhyML_Printf("%9.6f ", mat->dist[i][j]);
    PhyML_Printf("\n");
  }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Print_Node(t_node *a, t_node *d, t_tree *tree)
{
  int i;
  int dir;
  dir = -1;
  for (i = 0; i < 3; i++)
    if (a->v[i] == d)
    {
      dir = i;
      break;
    }
  PhyML_Printf("Node nums: %3d %3d  (dir:%3d) (a->anc:%3d) (d->anc:%3d) "
               "ta:%8.4f td:%8.4f t_min:%6.2f t_max:%6.2f",
               a->num, d->num, dir, a->anc ? a->anc->num : (-1),
               d->anc ? d->anc->num : (-1),
               tree->rates ? tree->times->nd_t[a->num] : -1.,
               tree->rates ? tree->times->nd_t[d->num] : -1.,
               tree->rates ? tree->times->t_prior_min[a->num] : -1.,
               tree->rates ? tree->times->t_prior_max[a->num] : -1.);

  PhyML_Printf(" names = '%10s' '%10s' ; ", a->name, d->name);
  for (i = 0; i < 3; i++)
    if (a->v[i] == d)
    {
      if (a->b[i])
      {
        PhyML_Printf("Branch num = %3d%c (%3d %3d) %f", a->b[i]->num,
                     a->b[i] == tree->e_root ? '*' : ' ', a->b[i]->left->num,
                     a->b[i]->rght->num, a->b[i]->l->v);
        if (a->b[i]->left->tax) PhyML_Printf(" WARNING LEFT->TAX!");
        break;
      }
    }
  PhyML_Printf("\n");

  if (d->tax)
    return;
  else
    for (i = 0; i < 3; i++)
      if (d->v[i] != a && d->b[i] != tree->e_root) Print_Node(d, d->v[i], tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Print_Node_Brief(t_node *a, t_node *d, t_tree *tree, FILE *fp)
{
  int i;
  int dir;

  dir = -1;
  for (i = 0; i < 3; i++)
    if (a->v[i] == d)
    {
      dir = i;
      break;
    }

  PhyML_Fprintf(fp, "\n");
  PhyML_Fprintf(fp, "Node nums: %3d %3d  (dir:%3d)", a->num, d->num, dir);

  PhyML_Fprintf(fp, "\tnames = '%10s' '%10s' ; ", a->name, d->name);
  for (i = 0; i < 3; i++)
    if (a->v[i] == d)
    {
      if (a->b[i])
      {
        PhyML_Fprintf(fp, "Branch num = %3d%c (%3d %3d) length:%10f",
                      a->b[i]->num, a->b[i] == tree->e_root ? '*' : ' ',
                      a->b[i]->left->num, a->b[i]->rght->num, a->b[i]->l->v);
        if (a->b[i]->left->tax) PhyML_Printf(" WARNING LEFT->TAX!");
        break;
      }
    }

  if (d->tax)
    return;
  else
    for (i = 0; i < 3; i++)
      if (d->v[i] != a && d->b[i] != tree->e_root)
        Print_Node_Brief(d, d->v[i], tree, fp);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Print_Model(t_mod *mod)
{
  int i, j, k;

  PhyML_Printf("\n. name=%s", mod->modelname->s);
  PhyML_Printf("\n. string=%s", mod->custom_mod_string);
  PhyML_Printf("\n. mod_num=%d", mod->mod_num);
  PhyML_Printf("\n. ns=%d", mod->ns);
  PhyML_Printf("\n. n_catg=%d", mod->ras->n_catg);
  PhyML_Printf("\n. kappa=%f", mod->kappa->v);
  PhyML_Printf("\n. alpha=%f", mod->ras->alpha->v);
  PhyML_Printf("\n. lambda=%f", mod->lambda->v);
  PhyML_Printf("\n. pinvar=%f", mod->ras->pinvar->v);
  PhyML_Printf("\n. br_len_mult=%f", mod->br_len_mult->v);
  PhyML_Printf("\n. whichmodel=%d", mod->whichmodel);
  PhyML_Printf("\n. update_eigen=%d", mod->update_eigen);
  PhyML_Printf("\n. bootstrap=%d", mod->io->n_boot_replicates);
  PhyML_Printf("\n. n_diff_rr=%d", mod->r_mat->n_diff_rr);
  PhyML_Printf("\n. invar=%d", mod->ras->invar);
  PhyML_Printf("\n. use_m4mod=%d", mod->use_m4mod);
  PhyML_Printf("\n. gamma_median=%d", mod->ras->gamma_median);
  PhyML_Printf("\n. state_len=%d", mod->io->state_len);
  PhyML_Printf("\n. log_l=%d", mod->log_l);
  PhyML_Printf("\n. l_min=%f", mod->l_min);
  PhyML_Printf("\n. l_max=%f", mod->l_max);
  PhyML_Printf("\n. free_mixt_rates=%d", mod->ras->free_mixt_rates);
  PhyML_Printf("\n. gamma_mgf_bl=%d", mod->gamma_mgf_bl);

  PhyML_Printf("\n. Pi\n");
  for (i = 0; i < mod->ns; i++) PhyML_Printf(" %f ", mod->e_frq->pi->v[i]);
  PhyML_Printf("\n");
  for (i = 0; i < mod->ns; i++)
    PhyML_Printf(" %f ", mod->e_frq->pi_unscaled->v[i]);

  PhyML_Printf("\n. Rates\n");
  for (i = 0; i < mod->ras->n_catg; i++)
    PhyML_Printf(" %f ", mod->ras->gamma_r_proba->v[i]);
  PhyML_Printf("\n");
  for (i = 0; i < mod->ras->n_catg; i++)
    PhyML_Printf(" %f ", mod->ras->gamma_r_proba_unscaled->v[i]);
  PhyML_Printf("\n");
  for (i = 0; i < mod->ras->n_catg; i++)
    PhyML_Printf(" %f ", mod->ras->gamma_rr->v[i]);
  PhyML_Printf("\n");
  for (i = 0; i < mod->ras->n_catg; i++)
    PhyML_Printf(" %f ", mod->ras->gamma_rr_unscaled->v[i]);

  PhyML_Printf("\n. Qmat \n");
  if (mod->whichmodel == CUSTOM)
  {
    fflush(NULL);
    for (i = 0; i < 6; i++)
    {
      PhyML_Printf(" %12f ", mod->r_mat->rr->v[i]);
      fflush(NULL);
    }
    for (i = 0; i < 6; i++)
    {
      PhyML_Printf(" %12f ", mod->r_mat->rr_val->v[i]);
      fflush(NULL);
    }
    for (i = 0; i < 6; i++)
    {
      PhyML_Printf(" %12d ", mod->r_mat->rr_num->v[i]);
      fflush(NULL);
    }
    for (i = 0; i < 6; i++)
    {
      PhyML_Printf(" %12d ", mod->r_mat->n_rr_per_cat->v[i]);
      fflush(NULL);
    }
  }
  for (i = 0; i < mod->ns; i++)
  {
    PhyML_Printf("  ");
    for (j = 0; j < 4; j++)
      PhyML_Printf("%8.5f  ", mod->r_mat->qmat->v[i * 4 + j]);
    PhyML_Printf("\n");
  }

  PhyML_Printf("\n. Freqs");
  PhyML_Printf("\n");
  for (i = 0; i < mod->ns; i++)
    PhyML_Printf(" %12f ", mod->e_frq->user_b_freq->v[i]);
  PhyML_Printf("\n");
  for (i = 0; i < mod->ns; i++) PhyML_Printf(" %12f ", mod->e_frq->pi->v[i]);
  PhyML_Printf("\n");
  for (i = 0; i < mod->ns; i++)
    PhyML_Printf(" %12f ", mod->e_frq->pi_unscaled->v[i]);

  PhyML_Printf("\n. Eigen\n");
  For(i, 2 * mod->ns) PhyML_Printf(" %f ", mod->eigen->space[i]);
  /* PhyML_Printf("\n"); */
  /* For(i,2*mod->ns)       PhyML_Printf(" %f ",mod->eigen->space_int[i]); */
  PhyML_Printf("\n");
  for (i = 0; i < mod->ns; i++) PhyML_Printf(" %f ", mod->eigen->e_val[i]);
  PhyML_Printf("\n");
  for (i = 0; i < mod->ns; i++) PhyML_Printf(" %f ", mod->eigen->e_val_im[i]);
  PhyML_Printf("\n");
  For(i, mod->ns * mod->ns) PhyML_Printf(" %f ", mod->eigen->r_e_vect[i]);
  PhyML_Printf("\n");
  For(i, mod->ns * mod->ns) PhyML_Printf(" %f ", mod->eigen->r_e_vect_im[i]);
  PhyML_Printf("\n");
  For(i, mod->ns * mod->ns) PhyML_Printf(" %f ", mod->eigen->l_e_vect[i]);
  PhyML_Printf("\n");
  For(i, mod->ns * mod->ns) PhyML_Printf(" %f ", mod->eigen->q[i]);
  PhyML_Printf("\n");

  PhyML_Printf("\n. Pij");
  for (k = 0; k < mod->ras->n_catg; k++)
  {
    PMat(0.01 * mod->ras->gamma_rr->v[k], mod, mod->ns * mod->ns * k,
         mod->Pij_rr->v, NULL);
    PhyML_Printf("\n. l=%f\n", 0.01 * mod->ras->gamma_rr->v[k]);
    for (i = 0; i < mod->ns; i++)
    {
      PhyML_Printf("  ");
      for (j = 0; j < mod->ns; j++)
        PhyML_Printf("%8.5f  ",
                     mod->Pij_rr->v[k * mod->ns * mod->ns + i * mod->ns + j]);
      PhyML_Printf("\n");
    }
  }

  PhyML_Printf("\n");

  fflush(NULL);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Print_Mat(matrix *mat)
{
  int i, j;

  PhyML_Printf("\n\n");
  PhyML_Printf("%d", mat->n_otu);
  PhyML_Printf("\n");

  for (i = 0; i < mat->n_otu; i++)
  {
    for (j = 0; j < 13; j++)
    {
      if (j >= (int)strlen(mat->name[i]))
        putchar(' ');
      else
        putchar(mat->name[i][j]);
    }

    for (j = 0; j < mat->n_otu; j++)
    {
      char s[2] = "-";
      if (mat->dist[i][j] < .0)
        PhyML_Printf("%12s", s);
      else
        PhyML_Printf("%12f", mat->dist[i][j]);
    }
    PhyML_Printf("\n");
  }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

FILE *Openfile(char *filename, int mode)
{
  FILE *fp;
  char *s;
  int   open_test = 0;

  s = filename;

  fp = NULL;

  switch (mode)
  {
  case READ:
  {
    while (!(fp = (FILE *)fopen(s, "r")) && ++open_test < 10)
    {
      PhyML_Printf("\n. Can't open file '%s', enter a new name : ", s);
      Getstring_Stdin(s);
    }
    break;
  }
  case WRITE:
  {
    fp = (FILE *)fopen(s, "w");
    break;
  }
  case APPEND:
  {
    fp = (FILE *)fopen(s, "a");
    break;
  }
  case READWRITE:
  {
    fp = (FILE *)fopen(s, "w+");
    break;
  }

  default:
    break;
  }

  /*   Free(s); */

  return fp;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Print_Fp_Out(FILE *fp_out, time_t t_beg, time_t t_end, t_tree *tree,
                  option *io, int n_data_set, int num_tree, int add_citation,
                  int precision)
{
  char *s;
  char  format[8];
  div_t hour, min;
  int   i, j;

  if (precision > 0) snprintf(format, 8, "%%.%huf", (unsigned short)precision);

  if (n_data_set == 1)
  {
    rewind(fp_out);
    Print_Banner_Small(fp_out);
  }

  PhyML_Fprintf(fp_out, "\n. Sequence filename: \t\t\t%s",
                Basename(io->in_align_file));
  PhyML_Fprintf(fp_out, "\n. Data set: \t\t\t\t#%d", n_data_set);

  if (io->mod->s_opt->random_input_tree)
    PhyML_Fprintf(fp_out, "\n. Random init tree: \t\t\t#%d", num_tree + 1);
  else if (io->n_trees > 1)
    PhyML_Fprintf(fp_out, "\n. Starting tree number: \t\t#%d", num_tree + 1);

  /* if(io->mod->s_opt->opt_topo) */
  /*   PhyML_Fprintf(fp_out,"\n. Tree topology search: \t\tSPRs"); */
  /* else */
  /*   PhyML_Fprintf(fp_out,"\n. Tree topology: \t\t\tfixed"); */

  /* was after Sequence file ; moved here FLT */
  s = (char *)mCalloc(T_MAX_LINE, sizeof(char));
  if (io->in_tree == 2)
  {
    strcat(strcat(strcat(s, "user tree ("), io->in_tree_file), ")");
  }
  else
  {
    if (!io->mod->s_opt->random_input_tree)
    {
      if (io->in_tree == 0) strcat(s, "BioNJ");
      if (io->in_tree == 1) strcat(s, "parsimony");
    }
    else
      strcat(s, "random tree");
  }

  PhyML_Fprintf(fp_out, "\n. Initial tree: \t\t\t%s", s);
  Free(s);

  if (tree->io->datatype == NT)
  {
    PhyML_Fprintf(fp_out, "\n. Model of nucleotides substitution: \t%s",
                  tree->mod->modelname->s);
    if (io->mod->whichmodel == CUSTOM)
      PhyML_Fprintf(fp_out, " (%s)", io->mod->custom_mod_string);
  }
  else if (tree->io->datatype == AA)
  {
    PhyML_Fprintf(fp_out, "\n. Model of amino acids substitution: \t%s",
                  tree->mod->modelname->s);
    if (io->mod->whichmodel == CUSTOMAA)
      PhyML_Fprintf(fp_out, " (%s)", tree->mod->aa_rate_mat_file->s);
  }
  else
  {
    PhyML_Fprintf(fp_out, "\n. Substitution model: \t\t%s",
                  tree->mod->modelname->s);
  }

  PhyML_Fprintf(fp_out, "\n. Integrated length (IL) model: \t%s",
                (tree->mod->gamma_mgf_bl == YES) ? "yes" : "no");

  if (tree->mod->gamma_mgf_bl == YES)
    PhyML_Fprintf(fp_out, "\n. Variance of IL model: \t\t%f",
                  tree->mod->l_var_sigma->v);

  PhyML_Fprintf(fp_out, "\n. Number of taxa: \t\t\t%d",
                tree->n_otu); /*added FLT*/

  PhyML_Fprintf(fp_out, "\n. Log-likelihood: \t\t\t%.5f",
                tree->c_lnL); /*was last ; moved here FLT*/

  Unconstraint_Lk(tree);
  PhyML_Fprintf(fp_out, "\n. Unconstrained log-likelihood: \t%.5f",
                tree->unconstraint_lk);

  Composite_Lk(tree);
  PhyML_Fprintf(fp_out, "\n. Composite log-likelihood: \t\t%.5f",
                tree->composite_lk);

  PhyML_Fprintf(fp_out, "\n. Parsimony: \t\t\t\t%d", tree->c_pars);

  PhyML_Fprintf(fp_out, "\n. Tree size: \t\t\t\t%.5f", Get_Tree_Size(tree));

  /* if(tree->mod->ras->n_catg > 1 && tree->mod->ras->free_mixt_rates == NO) */
  if (tree->mod->ras->free_mixt_rates == NO)
  {
    PhyML_Fprintf(fp_out, "\n. Discrete gamma model: \t\t%s", "Yes");
    PhyML_Fprintf(fp_out, "\n  - Number of classes: \t\t\t%d",
                  tree->mod->ras->n_catg);
    PhyML_Fprintf(fp_out, "\n  - Gamma shape parameter: \t\t%.3f",
                  tree->mod->ras->alpha->v);
    for (i = 0; i < tree->mod->ras->n_catg; i++)
    {
      PhyML_Fprintf(fp_out,
                    "\n  - Relative rate in class %d: \t\t%.5f [freq=%4f] \t\t",
                    i + 1, tree->mod->ras->gamma_rr->v[i],
                    tree->mod->ras->gamma_r_proba->v[i]);
    }
  }
  else if (tree->mod->ras->free_mixt_rates == YES)
  {
    int *rk;
    rk = Ranks(tree->mod->ras->gamma_rr->v, tree->mod->ras->n_catg);
    PhyML_Fprintf(fp_out, "\n. FreeRate model: \t\t\t%s", "Yes");
    PhyML_Fprintf(fp_out, "\n  - Number of classes: \t\t\t%d",
                  tree->mod->ras->n_catg);
    for (i = 0; i < tree->mod->ras->n_catg; i++)
    {
      PhyML_Fprintf(fp_out,
                    "\n  - Relative rate in class %d: \t\t%.5f [freq=%4f] \t\t",
                    i + 1, tree->mod->ras->gamma_rr->v[rk[i]],
                    tree->mod->ras->gamma_r_proba->v[rk[i]]);
    }
    Free(rk);
  }

  if (tree->mod->ras->invar)
    PhyML_Fprintf(fp_out, "\n. Proportion of invariant: \t\t%.3f",
                  tree->mod->ras->pinvar->v);

  /*was before Discrete gamma model ; moved here FLT*/
  if ((tree->mod->whichmodel == K80) || (tree->mod->whichmodel == HKY85) ||
      (tree->mod->whichmodel == F84))
  {
    PhyML_Fprintf(fp_out, "\n. Transition/transversion ratio: \t");
    if (precision > 0)
      PhyML_Fprintf(fp_out, format, tree->mod->kappa->v);
    else
      PhyML_Fprintf(fp_out, "%f", tree->mod->kappa->v);
  }
  else if (tree->mod->whichmodel == TN93)
  {
    PhyML_Fprintf(fp_out,
                  "\n. Transition/transversion ratio for purines: \t\t");
    if (precision > 0)
      PhyML_Fprintf(fp_out, format,
                    tree->mod->kappa->v * 2. * tree->mod->lambda->v /
                        (1. + tree->mod->lambda->v));
    else
      PhyML_Fprintf(fp_out, "%f",
                    tree->mod->kappa->v * 2. * tree->mod->lambda->v /
                        (1. + tree->mod->lambda->v));

    PhyML_Fprintf(fp_out,
                  "\n. Transition/transversion ratio for pyrimidines: \t");
    if (precision > 0)
      PhyML_Fprintf(fp_out, format,
                    tree->mod->kappa->v * 2. / (1. + tree->mod->lambda->v));
    else
      PhyML_Fprintf(fp_out, "%f",
                    tree->mod->kappa->v * 2. / (1. + tree->mod->lambda->v));
  }

  if (tree->io->datatype == NT)
  {
    PhyML_Fprintf(fp_out, "\n. Nucleotides frequencies:");
    if (precision > 0)
    {
      PhyML_Fprintf(fp_out, "\n  - f(A)=  ");
      PhyML_Fprintf(fp_out, format, tree->mod->e_frq->pi->v[0]);
      PhyML_Fprintf(fp_out, "\n  - f(C)=  ");
      PhyML_Fprintf(fp_out, format, tree->mod->e_frq->pi->v[1]);
      PhyML_Fprintf(fp_out, "\n  - f(G)=  ");
      PhyML_Fprintf(fp_out, format, tree->mod->e_frq->pi->v[2]);
      PhyML_Fprintf(fp_out, "\n  - f(T)=  ");
      PhyML_Fprintf(fp_out, format, tree->mod->e_frq->pi->v[3]);
    }
    else
    {
      PhyML_Fprintf(fp_out, "\n  - f(A)= %8.5f", tree->mod->e_frq->pi->v[0]);
      PhyML_Fprintf(fp_out, "\n  - f(C)= %8.5f", tree->mod->e_frq->pi->v[1]);
      PhyML_Fprintf(fp_out, "\n  - f(G)= %8.5f", tree->mod->e_frq->pi->v[2]);
      PhyML_Fprintf(fp_out, "\n  - f(T)= %8.5f", tree->mod->e_frq->pi->v[3]);
    }
  }

  /*****************************************/
  if ((tree->mod->whichmodel == GTR) || (tree->mod->whichmodel == CUSTOM))
  {
    Update_Qmat_GTR(tree->mod->r_mat->rr->v, tree->mod->r_mat->rr_val->v,
                    tree->mod->r_mat->rr_num->v, tree->mod->e_frq->pi->v,
                    tree->mod->r_mat->qmat->v);

    PhyML_Fprintf(fp_out, "\n");
    PhyML_Fprintf(fp_out, ". GTR relative rate parameters :");
    if (precision > 0)
    {
      PhyML_Fprintf(fp_out, "\n  A <-> C   ");
      if (tree->mod->r_mat->rr->v[0] < 10) PhyML_Fprintf(fp_out, " ");
      PhyML_Fprintf(fp_out, format, tree->mod->r_mat->rr->v[0]);
      PhyML_Fprintf(fp_out, "\n  A <-> G   ");
      if (tree->mod->r_mat->rr->v[1] < 10) PhyML_Fprintf(fp_out, " ");
      PhyML_Fprintf(fp_out, format, tree->mod->r_mat->rr->v[1]);
      PhyML_Fprintf(fp_out, "\n  A <-> T   ");
      if (tree->mod->r_mat->rr->v[2] < 10) PhyML_Fprintf(fp_out, " ");
      PhyML_Fprintf(fp_out, format, tree->mod->r_mat->rr->v[2]);
      PhyML_Fprintf(fp_out, "\n  C <-> G   ");
      if (tree->mod->r_mat->rr->v[3] < 10) PhyML_Fprintf(fp_out, " ");
      PhyML_Fprintf(fp_out, format, tree->mod->r_mat->rr->v[3]);
      PhyML_Fprintf(fp_out, "\n  C <-> T   ");
      if (tree->mod->r_mat->rr->v[4] < 10) PhyML_Fprintf(fp_out, " ");
      PhyML_Fprintf(fp_out, format, tree->mod->r_mat->rr->v[4]);
      PhyML_Fprintf(fp_out, "\n  G <-> T   ");
      if (tree->mod->r_mat->rr->v[5] < 10) PhyML_Fprintf(fp_out, " ");
      PhyML_Fprintf(fp_out, format, tree->mod->r_mat->rr->v[5]);
    }
    else
    {
      PhyML_Fprintf(fp_out, "\n  A <-> C   %8.5f", tree->mod->r_mat->rr->v[0]);
      PhyML_Fprintf(fp_out, "\n  A <-> G   %8.5f", tree->mod->r_mat->rr->v[1]);
      PhyML_Fprintf(fp_out, "\n  A <-> T   %8.5f", tree->mod->r_mat->rr->v[2]);
      PhyML_Fprintf(fp_out, "\n  C <-> G   %8.5f", tree->mod->r_mat->rr->v[3]);
      PhyML_Fprintf(fp_out, "\n  C <-> T   %8.5f", tree->mod->r_mat->rr->v[4]);
      PhyML_Fprintf(fp_out, "\n  G <-> T   %8.5f", tree->mod->r_mat->rr->v[5]);
    }

    PhyML_Fprintf(fp_out, "\n. Instantaneous rate matrix : ");
    if (precision > 0)
    {
      PhyML_Fprintf(fp_out, "\n  [A");
      for (i = 0; i < precision + 4; i++) PhyML_Fprintf(fp_out, "-");
      PhyML_Fprintf(fp_out, "C");
      for (i = 0; i < precision + 4; i++) PhyML_Fprintf(fp_out, "-");
      PhyML_Fprintf(fp_out, "G");
      for (i = 0; i < precision + 4; i++) PhyML_Fprintf(fp_out, "-");
      PhyML_Fprintf(fp_out, "T");
      for (i = 0; i < precision + 1; i++) PhyML_Fprintf(fp_out, "-");
      PhyML_Fprintf(fp_out, "]\n");
      for (i = 0; i < 4; i++)
      {
        for (j = 0; j < 4; j++)
        {
          if (i == j)
            PhyML_Fprintf(fp_out, "  ");
          else
            PhyML_Fprintf(fp_out, "   ");
          PhyML_Fprintf(fp_out, format, tree->mod->r_mat->qmat->v[i * 4 + j]);
        }
        PhyML_Fprintf(fp_out, "\n");
      }
    }
    else
    {
      PhyML_Fprintf(fp_out, "\n  [A---------C---------G---------T------]\n");
      for (i = 0; i < 4; i++)
      {
        PhyML_Fprintf(fp_out, "  ");
        for (j = 0; j < 4; j++)
          PhyML_Fprintf(fp_out, "%8.5f  ",
                        tree->mod->r_mat->qmat->v[i * 4 + j]);
        PhyML_Fprintf(fp_out, "\n");
      }
    }
    // PhyML_Fprintf(fp_out,"\n");
  }

  if (tree->io->datatype == AA &&
      (tree->mod->e_frq->type == ML || tree->mod->e_frq->type == EMPIRICAL))
  {
    PhyML_Fprintf(fp_out, "\n. Amino-acid frequencies");
    if (precision > 0)
    {
      PhyML_Fprintf(fp_out, "\n- f(Ala)= ");
      PhyML_Fprintf(fp_out, format, tree->mod->e_frq->pi->v[0]);
      PhyML_Fprintf(fp_out, " f(Arg)= ");
      PhyML_Fprintf(fp_out, format, tree->mod->e_frq->pi->v[1]);
      PhyML_Fprintf(fp_out, " f(Asn)= ");
      PhyML_Fprintf(fp_out, format, tree->mod->e_frq->pi->v[2]);
      PhyML_Fprintf(fp_out, "\n- f(Asp)= ");
      PhyML_Fprintf(fp_out, format, tree->mod->e_frq->pi->v[3]);
      PhyML_Fprintf(fp_out, " f(Cys)= ");
      PhyML_Fprintf(fp_out, format, tree->mod->e_frq->pi->v[4]);
      PhyML_Fprintf(fp_out, " f(Gln)= ");
      PhyML_Fprintf(fp_out, format, tree->mod->e_frq->pi->v[5]);
      PhyML_Fprintf(fp_out, "\n- f(Glu)= ");
      PhyML_Fprintf(fp_out, format, tree->mod->e_frq->pi->v[6]);
      PhyML_Fprintf(fp_out, " f(Gly)= ");
      PhyML_Fprintf(fp_out, format, tree->mod->e_frq->pi->v[7]);
      PhyML_Fprintf(fp_out, " f(His)= ");
      PhyML_Fprintf(fp_out, format, tree->mod->e_frq->pi->v[8]);
      PhyML_Fprintf(fp_out, "\n- f(Ile)= ");
      PhyML_Fprintf(fp_out, format, tree->mod->e_frq->pi->v[9]);
      PhyML_Fprintf(fp_out, " f(Leu)= ");
      PhyML_Fprintf(fp_out, format, tree->mod->e_frq->pi->v[10]);
      PhyML_Fprintf(fp_out, " f(Lys)= ");
      PhyML_Fprintf(fp_out, format, tree->mod->e_frq->pi->v[11]);
      PhyML_Fprintf(fp_out, "\n- f(Met)= ");
      PhyML_Fprintf(fp_out, format, tree->mod->e_frq->pi->v[12]);
      PhyML_Fprintf(fp_out, " f(Phe)= ");
      PhyML_Fprintf(fp_out, format, tree->mod->e_frq->pi->v[13]);
      PhyML_Fprintf(fp_out, " f(Pro)= ");
      PhyML_Fprintf(fp_out, format, tree->mod->e_frq->pi->v[14]);
      PhyML_Fprintf(fp_out, "\n- f(Ser)= ");
      PhyML_Fprintf(fp_out, format, tree->mod->e_frq->pi->v[15]);
      PhyML_Fprintf(fp_out, " f(Thr)= ");
      PhyML_Fprintf(fp_out, format, tree->mod->e_frq->pi->v[16]);
      PhyML_Fprintf(fp_out, " f(Trp)= ");
      PhyML_Fprintf(fp_out, format, tree->mod->e_frq->pi->v[17]);
      PhyML_Fprintf(fp_out, "\n- f(Tyr)= ");
      PhyML_Fprintf(fp_out, format, tree->mod->e_frq->pi->v[18]);
      PhyML_Fprintf(fp_out, " f(Val)= ");
      PhyML_Fprintf(fp_out, format, tree->mod->e_frq->pi->v[19]);
    }
    else
    {
      PhyML_Fprintf(fp_out, "\n- f(Ala)= ");
      PhyML_Fprintf(fp_out, "%f", tree->mod->e_frq->pi->v[0]);
      PhyML_Fprintf(fp_out, " f(Arg)= ");
      PhyML_Fprintf(fp_out, "%f", tree->mod->e_frq->pi->v[1]);
      PhyML_Fprintf(fp_out, " f(Asn)= ");
      PhyML_Fprintf(fp_out, "%f", tree->mod->e_frq->pi->v[2]);
      PhyML_Fprintf(fp_out, "\n- f(Asp)= ");
      PhyML_Fprintf(fp_out, "%f", tree->mod->e_frq->pi->v[3]);
      PhyML_Fprintf(fp_out, " f(Cys)= ");
      PhyML_Fprintf(fp_out, "%f", tree->mod->e_frq->pi->v[4]);
      PhyML_Fprintf(fp_out, " f(Gln)= ");
      PhyML_Fprintf(fp_out, "%f", tree->mod->e_frq->pi->v[5]);
      PhyML_Fprintf(fp_out, "\n- f(Glu)= ");
      PhyML_Fprintf(fp_out, "%f", tree->mod->e_frq->pi->v[6]);
      PhyML_Fprintf(fp_out, " f(Gly)= ");
      PhyML_Fprintf(fp_out, "%f", tree->mod->e_frq->pi->v[7]);
      PhyML_Fprintf(fp_out, " f(His)= ");
      PhyML_Fprintf(fp_out, "%f", tree->mod->e_frq->pi->v[8]);
      PhyML_Fprintf(fp_out, "\n- f(Ile)= ");
      PhyML_Fprintf(fp_out, "%f", tree->mod->e_frq->pi->v[9]);
      PhyML_Fprintf(fp_out, " f(Leu)= ");
      PhyML_Fprintf(fp_out, "%f", tree->mod->e_frq->pi->v[10]);
      PhyML_Fprintf(fp_out, " f(Lys)= ");
      PhyML_Fprintf(fp_out, "%f", tree->mod->e_frq->pi->v[11]);
      PhyML_Fprintf(fp_out, "\n- f(Met)= ");
      PhyML_Fprintf(fp_out, "%f", tree->mod->e_frq->pi->v[12]);
      PhyML_Fprintf(fp_out, " f(Phe)= ");
      PhyML_Fprintf(fp_out, "%f", tree->mod->e_frq->pi->v[13]);
      PhyML_Fprintf(fp_out, " f(Pro)= ");
      PhyML_Fprintf(fp_out, "%f", tree->mod->e_frq->pi->v[14]);
      PhyML_Fprintf(fp_out, "\n- f(Ser)= ");
      PhyML_Fprintf(fp_out, "%f", tree->mod->e_frq->pi->v[15]);
      PhyML_Fprintf(fp_out, " f(Thr)= ");
      PhyML_Fprintf(fp_out, "%f", tree->mod->e_frq->pi->v[16]);
      PhyML_Fprintf(fp_out, " f(Trp)= ");
      PhyML_Fprintf(fp_out, "%f", tree->mod->e_frq->pi->v[17]);
      PhyML_Fprintf(fp_out, "\n- f(Tyr)= ");
      PhyML_Fprintf(fp_out, "%f", tree->mod->e_frq->pi->v[18]);
      PhyML_Fprintf(fp_out, " f(Val)= ");
      PhyML_Fprintf(fp_out, "%f", tree->mod->e_frq->pi->v[19]);
    }
  }

  /*****************************************/
  if (io->ratio_test == 1)
  {
    PhyML_Fprintf(fp_out, ". aLRT statistics to test branches");
  }
  else if (io->ratio_test == 2)
  {
    PhyML_Fprintf(fp_out, ". aLRT branch supports (cubic approximation, "
                          "mixture of Chi2s distribution)");
  }

  PhyML_Fprintf(fp_out, "\n");
  PhyML_Fprintf(fp_out, "\n. Run ID:\t\t\t\t%s",
                (io->append_run_ID) ? (io->run_id_string) : ("none"));
  PhyML_Fprintf(fp_out, "\n. Random seed:\t\t\t\t%d", io->r_seed);
  PhyML_Fprintf(fp_out, "\n. Subtree patterns aliasing:\t\t%s",
                io->do_alias_subpatt ? "yes" : "no");
  PhyML_Fprintf(fp_out, "\n. Version:\t\t\t\t%s", VERSION);

  hour = div(t_end - t_beg, 3600);
  min  = div(t_end - t_beg, 60);
  min.quot -= hour.quot * 60;

  PhyML_Fprintf(fp_out, "\n. Time used:\t\t\t\t%dh%dm%ds (%d seconds)",
                hour.quot, min.quot, (int)(t_end - t_beg) % 60,
                (int)(t_end - t_beg));

  if (add_citation == YES)
  {
    PhyML_Fprintf(fp_out, "\n\n");
    PhyML_Fprintf(fp_out, " ooooooooooooooooooooooooooooooooooooooooooooooooooo"
                          "ooooooooooooooooooooooooooooooooooooooooooooo\n");
    PhyML_Fprintf(fp_out, " Suggested citations:\n");
    PhyML_Fprintf(fp_out, " S. Guindon, JF. Dufayard, V. Lefort, M. Anisimova, "
                          "W. Hordijk, O. Gascuel\n");
    PhyML_Fprintf(
        fp_out, " \"New algorithms and methods to estimate maximum-likelihood "
                "phylogenies: assessing the performance of PhyML 3.0.\"\n");
    PhyML_Fprintf(fp_out, " Systematic Biology. 2010. 59(3):307-321.\n");
    PhyML_Fprintf(fp_out, "\n");
    PhyML_Fprintf(fp_out, " S. Guindon & O. Gascuel\n");
    PhyML_Fprintf(fp_out,
                  " \"A simple, fast, and accurate algorithm to estimate large "
                  "phylogenies by maximum likelihood\"\n");
    PhyML_Fprintf(fp_out, " Systematic Biology. 2003. 52(5):696-704.\n");
    PhyML_Fprintf(fp_out, " ooooooooooooooooooooooooooooooooooooooooooooooooooo"
                          "ooooooooooooooooooooooooooooooooooooooooooooo\n");

    PhyML_Fprintf(fp_out,
                  "\n\n Please consider making a donation to support PhyML "
                  "development by clicking on the following link:\n");
    PhyML_Fprintf(fp_out,
                  "\n          https://fondation-cnrs.org/faire-un-don/");
    PhyML_Fprintf(fp_out, "\n\n              ~ Thank you ! ~\n\n");
  }
  else
  {
    PhyML_Fprintf(fp_out, "\n\n");
    PhyML_Fprintf(fp_out, " ooooooooooooooooooooooooooooooooooooooooooooooooooo"
                          "ooooooooooooooooooooooooooooooooooooooooooooo\n");
    PhyML_Fprintf(fp_out, "\n");
  }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/*FLT wrote this function*/
void Print_Fp_Out_Lines(FILE *fp_out, time_t t_beg, time_t t_end, t_tree *tree,
                        option *io, int n_data_set)
{
  char *s;
  /*div_t hour,min;*/

  if (n_data_set == 1)
  {

    PhyML_Fprintf(fp_out, ". Sequence file : [%s]\n\n",
                  Basename(io->in_align_file));

    if ((tree->io->datatype == NT) || (tree->io->datatype == AA))
    {
      (tree->io->datatype == NT)
          ? (PhyML_Fprintf(fp_out,
                           ". Model of nucleotides substitution : %s\n\n",
                           io->mod->modelname->s))
          : (PhyML_Fprintf(fp_out,
                           ". Model of amino acids substitution : %s\n\n",
                           io->mod->modelname->s));
    }

    s = (char *)mCalloc(100, sizeof(char));

    switch (io->in_tree)
    {
    case 0:
    {
      strcpy(s, "BioNJ");
      break;
    }
    case 1:
    {
      strcpy(s, "parsimony");
      break;
    }
    case 2:
    {
      strcpy(s, "user tree (");
      strcat(s, io->in_tree_file);
      strcat(s, ")");
      break;
    }
    }

    PhyML_Fprintf(fp_out, ". Initial tree : [%s]\n\n", s);

    Free(s);

    PhyML_Fprintf(fp_out, "\n");

    /*headline 1*/
    PhyML_Fprintf(fp_out, ". Data\t");

    PhyML_Fprintf(fp_out, "Nb of \t");

    PhyML_Fprintf(fp_out, "Likelihood\t");

    PhyML_Fprintf(fp_out, "Discrete   \t");

    if (tree->mod->ras->n_catg > 1)
      PhyML_Fprintf(fp_out, "Number of \tGamma shape\t");

    PhyML_Fprintf(fp_out, "Proportion of\t");

    if (tree->mod->whichmodel <= 6) PhyML_Fprintf(fp_out, "Transition/ \t");

    PhyML_Fprintf(fp_out, "Nucleotides frequencies               \t");

    if ((tree->mod->whichmodel == GTR) || (tree->mod->whichmodel == CUSTOM))
      PhyML_Fprintf(fp_out, "Instantaneous rate matrix              \t");

    /*    PhyML_Fprintf(fp_out,"Time\t");*/

    PhyML_Fprintf(fp_out, "\n");

    /*headline 2*/
    PhyML_Fprintf(fp_out, "  set\t");

    PhyML_Fprintf(fp_out, "taxa\t");

    PhyML_Fprintf(fp_out, "loglk     \t");

    PhyML_Fprintf(fp_out, "gamma model\t");

    if (tree->mod->ras->n_catg > 1)
      PhyML_Fprintf(fp_out, "categories\tparameter  \t");

    PhyML_Fprintf(fp_out, "invariant    \t");

    if (tree->mod->whichmodel <= 6) PhyML_Fprintf(fp_out, "transversion\t");

    PhyML_Fprintf(fp_out, "f(A)      f(C)      f(G)      f(T)    \t");

    if ((tree->mod->whichmodel == GTR) || (tree->mod->whichmodel == CUSTOM))
      PhyML_Fprintf(fp_out, "[A---------C---------G---------T------]\t");

    /*    PhyML_PhyML_Fprintf(fp_out,"used\t");*/

    PhyML_Fprintf(fp_out, "\n");

    /*headline 3*/
    if (tree->mod->whichmodel == TN93)
    {
      PhyML_Fprintf(fp_out, "    \t      \t          \t           \t");
      if (tree->mod->ras->n_catg > 1)
        PhyML_Fprintf(fp_out, "         \t         \t");
      PhyML_Fprintf(fp_out, "             \t");
      PhyML_Fprintf(fp_out, "purines pyrimid.\t");
      PhyML_Fprintf(fp_out, "\n");
    }

    PhyML_Fprintf(fp_out, "\n");
  }

  /*line items*/

  PhyML_Fprintf(fp_out, "  #%d\t", n_data_set);

  PhyML_Fprintf(fp_out, "%d   \t", tree->n_otu);

  PhyML_Fprintf(fp_out, "%.5f\t", tree->c_lnL);

  PhyML_Fprintf(fp_out, "%s        \t",
                (tree->mod->ras->n_catg > 1) ? ("Yes") : ("No "));
  if (tree->mod->ras->n_catg > 1)
  {
    PhyML_Fprintf(fp_out, "%d        \t", tree->mod->ras->n_catg);
    PhyML_Fprintf(fp_out, "%.3f    \t", tree->mod->ras->alpha->v);
  }

  /*if(tree->mod->ras->invar)*/
  PhyML_Fprintf(fp_out, "%.3f    \t", tree->mod->ras->pinvar->v);

  if (tree->mod->whichmodel <= 5)
  {
    PhyML_Fprintf(fp_out, "%.3f     \t", tree->mod->kappa->v);
  }
  else if (tree->mod->whichmodel == TN93)
  {
    PhyML_Fprintf(fp_out, "%.3f   ",
                  tree->mod->kappa->v * 2. * tree->mod->lambda->v /
                      (1. + tree->mod->lambda->v));
    PhyML_Fprintf(fp_out, "%.3f\t",
                  tree->mod->kappa->v * 2. / (1. + tree->mod->lambda->v));
  }

  if (tree->io->datatype == NT)
  {
    PhyML_Fprintf(fp_out, "%8.5f  ", tree->mod->e_frq->pi->v[0]);
    PhyML_Fprintf(fp_out, "%8.5f  ", tree->mod->e_frq->pi->v[1]);
    PhyML_Fprintf(fp_out, "%8.5f  ", tree->mod->e_frq->pi->v[2]);
    PhyML_Fprintf(fp_out, "%8.5f\t", tree->mod->e_frq->pi->v[3]);
  }
  /*
  hour = div(t_end-t_beg,3600);
  min  = div(t_end-t_beg,60  );

  min.quot -= hour.quot*60;

  PhyML_Fprintf(fp_out,"%dh%dm%ds\t", hour.quot,min.quot,(int)(t_end-t_beg)%60);
  if(t_end-t_beg > 60)
    PhyML_Fprintf(fp_out,". -> %d seconds\t",(int)(t_end-t_beg));
  */

  /*****************************************/
  if ((tree->mod->whichmodel == GTR) || (tree->mod->whichmodel == CUSTOM))
  {
    int i, j;

    for (i = 0; i < 4; i++)
    {
      if (i != 0)
      {
        /*format*/
        PhyML_Fprintf(fp_out, "      \t     \t          \t           \t");
        if (tree->mod->ras->n_catg > 1)
          PhyML_Fprintf(fp_out, "          \t           \t");
        PhyML_Fprintf(
            fp_out, "             \t                                      \t");
      }
      for (j = 0; j < 4; j++)
        PhyML_Fprintf(fp_out, "%8.5f  ", tree->mod->r_mat->qmat->v[i * 4 + j]);
      if (i < 3) PhyML_Fprintf(fp_out, "\n");
    }
  }
  /*****************************************/

  PhyML_Fprintf(fp_out, "\n\n");
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Print_Freq(t_tree *tree)
{

  switch (tree->io->datatype)
  {
  case NT:
  {
    PhyML_Printf("A : %f\n", tree->mod->e_frq->pi->v[0]);
    PhyML_Printf("C : %f\n", tree->mod->e_frq->pi->v[1]);
    PhyML_Printf("G : %f\n", tree->mod->e_frq->pi->v[2]);
    PhyML_Printf("T : %f\n", tree->mod->e_frq->pi->v[3]);
    break;
  }
  case AA:
  {
    PhyML_Printf("A : %f\n", tree->mod->e_frq->pi->v[0]);
    PhyML_Printf("R : %f\n", tree->mod->e_frq->pi->v[1]);
    PhyML_Printf("N : %f\n", tree->mod->e_frq->pi->v[2]);
    PhyML_Printf("D : %f\n", tree->mod->e_frq->pi->v[3]);
    PhyML_Printf("C : %f\n", tree->mod->e_frq->pi->v[4]);
    PhyML_Printf("Q : %f\n", tree->mod->e_frq->pi->v[5]);
    PhyML_Printf("E : %f\n", tree->mod->e_frq->pi->v[6]);
    PhyML_Printf("G : %f\n", tree->mod->e_frq->pi->v[7]);
    PhyML_Printf("H : %f\n", tree->mod->e_frq->pi->v[8]);
    PhyML_Printf("I : %f\n", tree->mod->e_frq->pi->v[9]);
    PhyML_Printf("L : %f\n", tree->mod->e_frq->pi->v[10]);
    PhyML_Printf("K : %f\n", tree->mod->e_frq->pi->v[11]);
    PhyML_Printf("M : %f\n", tree->mod->e_frq->pi->v[12]);
    PhyML_Printf("F : %f\n", tree->mod->e_frq->pi->v[13]);
    PhyML_Printf("P : %f\n", tree->mod->e_frq->pi->v[14]);
    PhyML_Printf("S : %f\n", tree->mod->e_frq->pi->v[15]);
    PhyML_Printf("T : %f\n", tree->mod->e_frq->pi->v[16]);
    PhyML_Printf("W : %f\n", tree->mod->e_frq->pi->v[17]);
    PhyML_Printf("Y : %f\n", tree->mod->e_frq->pi->v[18]);
    PhyML_Printf("V : %f\n", tree->mod->e_frq->pi->v[19]);
    PhyML_Printf("N : %f\n", tree->mod->e_frq->pi->v[20]);
    break;
  }
  default:
  {
    break;
  }
  }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Print_Settings(option *io)
{
  int   answer;
  char *s;

  s = (char *)mCalloc(100, sizeof(char));

  PhyML_Printf("\n\n\n");
  PhyML_Printf("\n\n");

  PhyML_Printf("\n  "
               "////////////////////////////////////"
               ".\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\"
               "\\\\\\\\\\\\\\\\\\\\\\");
  PhyML_Printf("\n  "
               "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\"
               "\\\\\\\\\\.//////////////////////////////////////////\n");

  PhyML_Printf("\n        . Sequence filename:\t\t\t\t %s",
               Basename(io->in_align_file));

  if (io->datatype == NT)
    strcpy(s, "dna");
  else if (io->datatype == AA)
    strcpy(s, "aa");
  else
    strcpy(s, "generic");

  PhyML_Printf("\n        . Data type:\t\t\t\t\t %s", s);
  PhyML_Printf("\n        . Alphabet size:\t\t\t\t %d", io->mod->ns);

  PhyML_Printf("\n        . Sequence format:\t\t\t\t %s",
               io->interleaved ? "interleaved" : "sequential");
  PhyML_Printf("\n        . Number of data sets:\t\t\t\t %d", io->n_data_sets);

  PhyML_Printf("\n        . Nb of bootstrapped data sets:\t\t\t %d",
               io->n_boot_replicates);

  if (io->do_bayesboot)
    PhyML_Printf("\n        . Compute Bayesian Bootstrap:\t\t\t yes");
  if (io->do_tbe)
    PhyML_Printf("\n        . Compute TBE:\t\t\t yes");

  if (io->n_boot_replicates > 0)
    PhyML_Printf("\n        . Compute approximate likelihood ratio test:\t no");
  else
  {
    if (io->ratio_test == 1)
      PhyML_Printf("\n        . Compute approximate likelihood ratio test:\t "
                   "yes (aLRT statistics)");
    else if (io->ratio_test == 2)
      PhyML_Printf("\n        . Compute approximate likelihood ratio test:\t "
                   "yes (Chi2-based parametric branch supports)");
    else if (io->ratio_test == 3)
      PhyML_Printf("\n        . Compute approximate likelihood ratio test:\t "
                   "yes (Minimum of SH-like and Chi2-based branch supports)");
    else if (io->ratio_test == 4)
      PhyML_Printf("\n        . Compute approximate likelihood ratio test:\t "
                   "yes (SH-like branch supports)");
    else if (io->ratio_test == 5)
      PhyML_Printf("\n        . Compute approximate likelihood ratio test:\t "
                   "yes (aBayes branch supports)");
  }

  PhyML_Printf("\n        . Model name:\t\t\t\t\t %s", io->mod->modelname->s);

  if (io->datatype == AA && io->mod->whichmodel == CUSTOMAA)
    PhyML_Printf(" (%s)", io->mod->aa_rate_mat_file->s);

  if (io->datatype == NT)
  {
    if ((io->mod->whichmodel == K80) || (io->mod->whichmodel == HKY85) ||
        (io->mod->whichmodel == F84) || (io->mod->whichmodel == TN93))
    {
      if (io->mod->s_opt && io->mod->kappa->optimize)
        PhyML_Printf("\n        . Ts/tv ratio:\t\t\t\t\t estimated");
      else
        PhyML_Printf("\n        . Ts/tv ratio:\t\t\t\t\t %f",
                     io->mod->kappa->v);
    }
  }

  if (io->mod->s_opt && io->mod->ras->pinvar->optimize)
    PhyML_Printf("\n        . Proportion of invariable sites:\t\t estimated");
  else
    PhyML_Printf("\n        . Proportion of invariable sites:\t\t %f",
                 io->mod->ras->pinvar->v);

  if (io->mod->ras->free_mixt_rates == NO)
    PhyML_Printf("\n        . RAS model:\t\t\t\t\t discrete Gamma");
  else
    PhyML_Printf("\n        . RAS model:\t\t\t\t\t FreeRate");

  PhyML_Printf("\n        . Number of subst. rate catgs:\t\t\t %d",
               io->mod->ras->n_catg);

  if (io->mod->ras->n_catg > 1)
  {

    if (io->mod->ras->free_mixt_rates == NO)
    {
      if (io->mod->s_opt && io->mod->ras->alpha->optimize)
        PhyML_Printf(
            "\n        . Gamma distribution parameter:\t\t\t estimated");
      else
        PhyML_Printf("\n        . Gamma distribution parameter:\t\t\t %f",
                     io->mod->ras->alpha->v);
      PhyML_Printf("\n        . 'Middle' of each rate class:\t\t\t %s",
                   (io->mod->ras->gamma_median) ? ("median") : ("mean"));
    }
  }

  if (io->datatype == AA)
  {
    if (io->mod->e_frq->type == ML)
      PhyML_Printf(
          "\n        . Amino-acid equilibrium frequencies:\t\t optimized");
    else if (io->mod->e_frq->type == EMPIRICAL)
      PhyML_Printf(
          "\n        . Amino-acid equilibrium frequencies:\t\t empirical");
    else if (io->mod->e_frq->type == USER)
      PhyML_Printf(
          "\n        . Amino-acid equilibrium frequencies:\t\t user-defined");
    else if (io->mod->e_frq->type == MODEL)
      PhyML_Printf(
          "\n        . Amino-acid equilibrium frequencies:\t\t model-defined");
  }
  else if (io->datatype == NT)
  {
    if ((io->mod->whichmodel != JC69) && (io->mod->whichmodel != K80) &&
        (io->mod->whichmodel != F81))
    {
      if (io->mod->e_frq->type == ML)
        PhyML_Printf(
            "\n        . Nucleotide equilibrium frequencies:\t\t optimized");
      else if (io->mod->e_frq->type == EMPIRICAL)
        PhyML_Printf(
            "\n        . Nucleotide equilibrium frequencies:\t\t empirical");
      else if (io->mod->e_frq->type == USER)
        PhyML_Printf(
            "\n        . Nucleotide equilibrium frequencies:\t\t user-defined");
      else if (io->mod->e_frq->type == MODEL)
        PhyML_Printf("\n        . Nucleotide equilibrium frequencies:\t\t "
                     "model-defined");
    }
  }

  PhyML_Printf("\n        . Optimise tree topology:\t\t\t %s",
               (io->mod->s_opt && io->mod->s_opt->opt_topo) ? "yes" : "no");

  switch (io->in_tree)
  {
  case 0:
  {
    strcpy(s, "BioNJ");
    break;
  }
  case 1:
  {
    strcpy(s, "parsimony");
    break;
  }
  case 2:
  {
    strcpy(s, "user tree (");
    strcat(s, Basename(io->in_tree_file));
    strcat(s, ")");
    break;
  }
  }

  if (io->mod->s_opt && io->mod->s_opt->opt_topo)
  {
    /* if(io->mod->s_opt->topo_search == NNI_MOVE) PhyML_Printf("\n        .
     * Tree topology search:\t\t\t\t NNIs"); */
    /* else if(io->mod->s_opt->topo_search == SPR_MOVE) PhyML_Printf("\n . Tree
     * topology search:\t\t\t\t SPRs"); */
    /* else if(io->mod->s_opt->topo_search == BEST_OF_NNI_AND_SPR)
     * PhyML_Printf("\n        . Tree topology search:\t\t\t\t Best of NNIs and
     * SPRs"); */

    PhyML_Printf("\n        . Starting tree:\t\t\t\t %s", s);

    PhyML_Printf("\n        . Add random input tree:\t\t\t %s",
                 (io->mod->s_opt->random_input_tree) ? "yes" : "no");
    if (io->mod->s_opt->random_input_tree)
      PhyML_Printf("\n        . Number of random starting trees:\t\t %d",
                   io->mod->s_opt->n_rand_starts);
  }
  else if (io->mod->s_opt && !io->mod->s_opt->random_input_tree)
    PhyML_Printf("\n        . Evaluated tree:\t\t\t\t \"%s\"", s);

  PhyML_Printf("\n        . Optimise branch lengths:\t\t\t %s",
               (io->mod->s_opt && io->mod->s_opt->opt_bl_one_by_one) ? "yes"
                                                                     : "no");

  PhyML_Printf("\n        . Minimum length of an edge:\t\t\t %g",
               io->mod->l_min);

  answer = NO;
  if (io->mod->ras->alpha->optimize || io->mod->kappa->optimize ||
      io->mod->lambda->optimize || io->mod->ras->pinvar->optimize ||
      io->mod->r_mat->optimize)
    answer = YES;

  PhyML_Printf("\n        . Optimise substitution model parameters:\t %s",
               (answer) ? "yes" : "no");

  PhyML_Printf("\n        . Run ID:\t\t\t\t\t %s",
               (io->append_run_ID) ? (io->run_id_string) : ("none"));
  PhyML_Printf("\n        . Random seed:\t\t\t\t\t %d", io->r_seed);
  PhyML_Printf("\n        . Subtree patterns aliasing:\t\t\t %s",
               io->do_alias_subpatt ? "yes" : "no");
  PhyML_Printf("\n        . Version:\t\t\t\t\t %s", VERSION);
  PhyML_Printf("\n        . Byte alignment:\t\t\t\t %d", BYTE_ALIGN);
  PhyML_Printf("\n        . AVX enabled:\t\t\t\t\t %s",
#if defined(__AVX__)
               "yes"
#else
               "no"
#endif
  );
  PhyML_Printf("\n        . SSE enabled:\t\t\t\t\t %s",
#if defined(__SSE3__)
               "yes"
#else
               "no"
#endif
  );

  /* PhyML_Printf("\n\n                         \u205C \u205C \u205C \u205C
   * \u205C \u205C \u205C \u205C \u205C \u205C \u205C \u205C  \n"); */

  PhyML_Printf("\n");
  PhyML_Printf("\n  "
               "////////////////////////////////////"
               ".\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\"
               "\\\\\\\\\\\\\\\\\\\\\\");
  PhyML_Printf("\n  "
               "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\"
               "\\\\\\\\\\.//////////////////////////////////////////\n");

  PhyML_Printf("\n\n");
  fflush(NULL);

  Free(s);
}

void Print_Banner(FILE *fp)
{
  PhyML_Fprintf(fp, "\n");
  PhyML_Fprintf(fp, " ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo"
                    "ooooooooooooooooooooooooooooooooooooooo\n");
  PhyML_Fprintf(fp, "                                                          "
                    "                                        \n");
  PhyML_Fprintf(fp,
                "                                 ---  PhyML %s  ---           "
                "                                  \n",
                VERSION);
  PhyML_Fprintf(fp, "                                                          "
                    "                                        \n");
  PhyML_Fprintf(fp, "    A simple, fast, and accurate algorithm to estimate "
                    "large phylogenies by maximum likelihood    \n");
  PhyML_Fprintf(fp, "                            Stephane Guindon & Olivier "
                    "Gascuel                                      \n");
  PhyML_Fprintf(fp, "                                                          "
                    "                                        \n");
  PhyML_Fprintf(
      fp,
      "                           https://github.com/stephaneguindon/phyml/    "
      "                      \n");
  PhyML_Fprintf(fp, "                                                          "
                    "                                        \n");
  PhyML_Fprintf(fp, "                            Copyright CNRS - Universite "
                    "Montpellier                                  \n");
  PhyML_Fprintf(fp, "                                                          "
                    "                                        \n");
  PhyML_Fprintf(fp, " ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo"
                    "ooooooooooooooooooooooooooooooooooooooo\n");
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Print_Banner_Small(FILE *fp)
{
  PhyML_Fprintf(fp, "\n");
  PhyML_Fprintf(fp, " ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo"
                    "ooooooooooooooooooooooooooooooooooooooo\n");
  PhyML_Fprintf(fp,
                "                                  ---  PhyML %s  ---          "
                "                                   \n",
                VERSION);
  PhyML_Fprintf(
      fp, "                            "
          "https://github.com/stephaneguindon/phyml/                      \n");
  PhyML_Fprintf(fp, "                             Copyright CNRS - Universite "
                    "Montpellier                                 \n");
  PhyML_Fprintf(fp, " ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo"
                    "ooooooooooooooooooooooooooooooooooooooo\n");
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Print_Data_Set_Number(option *io, FILE *fp)
{
  PhyML_Fprintf(fp, "\n");
  PhyML_Fprintf(fp, "                                                          "
                    "                                        \n");
  PhyML_Fprintf(fp,
                "                                 [ Data set number %3d ]      "
                "                                     \n",
                io->curr_gt + 1);
  PhyML_Fprintf(fp, "                                                          "
                    "                                        \n");
}
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Print_Lk(t_tree *tree, char *string)
{
  t_tree *loc_tree;

  loc_tree = tree;
  /*! Rewind back to the first mixt_tree */
  while (loc_tree->prev) loc_tree = loc_tree->prev;

  time(&(loc_tree->t_current));
  PhyML_Printf("\n. (%5d sec) [%15.4f] %s",
               (int)(loc_tree->t_current - loc_tree->t_beg), Get_Lk(tree),
               string);
#ifndef QUIET
  fflush(NULL);
#endif
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Print_List(t_ll *list)
{
  t_ll *ll;
  ll = list->head;
  do
  {
    t_node *n;
    n = (t_node *)ll->v;
    PhyML_Printf("\n. list elem: %p next: %p prev: %p [%d] %p %p", (void *)ll,
                 (void *)ll->next, (void *)ll->prev, n->num, (void *)ll->head,
                 (void *)ll->tail);
    ll = ll->next;
  } while (ll != NULL);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Print_Pars(t_tree *tree)
{
  time(&(tree->t_current));
  PhyML_Printf("\n. (%5d sec) [%5d]", (int)(tree->t_current - tree->t_beg),
               tree->c_pars);
#ifndef QUIET
  fflush(NULL);
#endif
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Print_Lk_And_Pars(t_tree *tree)
{
  time(&(tree->t_current));

  PhyML_Printf("\n. (%5d sec) [%15.4f] [%5d]",
               (int)(tree->t_current - tree->t_beg), tree->c_lnL, tree->c_pars);

#ifndef QUIET
  fflush(NULL);
#endif
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Read_Qmat(phydbl *daa, phydbl *pi, FILE *fp)
{
  int    i, j;
  phydbl sum;
  double val;

  assert(fp);
  rewind(fp);

  for (i = 1; i < 20; i++)
  {
    for (j = 0; j < 19; j++)
    {
      /* 	  if(!fscanf(fp,"%lf",&(daa[i*20+j]))) Exit("\n"); */
      if (!fscanf(fp, "%lf", &val))
      {
        PhyML_Fprintf(stderr,
                      "\n. Rate matrix file does not appear to have a proper "
                      "format. Please refer to the documentation.");
        Exit("\n");
      }
      daa[i * 20 + j] = (phydbl)val;
      daa[j * 20 + i] = daa[i * 20 + j];
      if (j == i - 1) break;
    }
  }

  for (i = 0; i < 20; i++)
  {
    if (!fscanf(fp, "%lf", &val)) Exit("\n");
    pi[i] = (phydbl)val;
  }
  sum = .0;
  for (i = 0; i < 20; i++) sum += pi[i];
  if (FABS(sum - 1.) > 1.E-06)
  {
    PhyML_Printf("\n. Sum of amino-acid frequencies: %f", sum);
    PhyML_Printf("\n. Scaling amino-acid frequencies...\n");
    for (i = 0; i < 20; i++) pi[i] /= sum;
  }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Print_Qmat_AA(phydbl *daa, phydbl *pi)
{
  int i, j, cpt;

  cpt = 0;
  for (i = 0; i < 20; i++)
  {
    for (j = 0; j < i; j++)
    {
      PhyML_Printf("daa[%2d*20+%2d] = %10f;  ", i, j, daa[i * 20 + j]);
      cpt++;
      if (!(cpt % 4)) PhyML_Printf("\n");
    }
  }

  PhyML_Printf("\n\n");
  PhyML_Printf("for (i=0; i<naa; i++)  for (j=0; j<i; j++)  daa[j*naa+i] = "
               "daa[i*naa+j];\n\n");
  for (i = 0; i < 20; i++) PhyML_Printf("pi[%d] = %f; ", i, pi[i]);
  PhyML_Printf("\n");
  PhyML_Printf("Ala\tArg\tAsn\tAsp\tCys\tGln\tGlu\tGly\tHis\tIle\tLeu\tLys\tMet"
               "\tPhe\tPro\tSer\tThr\tTrp\tTyr\tVal\n");
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Print_Square_Matrix_Generic(int n, phydbl *mat)
{
  int i, j;

  PhyML_Printf("\n");
  for (i = 0; i < n; i++)
  {
    PhyML_Printf("[%3d]", i);
    for (j = 0; j < n; j++)
    {
      PhyML_Printf("%12.5f ", mat[i * n + j]);
    }
    PhyML_Printf("\n");
  }
  PhyML_Printf("\n");
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Print_Diversity(FILE *fp, t_tree *tree)
{

  Print_Diversity_Pre(tree->a_nodes[0], tree->a_nodes[0]->v[0],
                      tree->a_nodes[0]->b[0], fp, tree);

  /*       mean_div_left = .0; */
  /*       for(k=0;k<ns;k++)  */
  /* 	{ */
  /* 	  mean_div_left += (k+1) * tree->a_edges[j]->div_post_pred_left[k]; */
  /* 	} */
  /*       mean_div_rght = .0; */
  /*       for(k=0;k<ns;k++) mean_div_rght += (k+1) *
   * tree->a_edges[j]->div_post_pred_rght[k]; */

  /*       mean_div_left /= (phydbl)tree->data->init_len; */
  /*       mean_div_rght /= (phydbl)tree->data->init_len; */

  /*       PhyML_Fprintf(fp,"%4d 0 %f\n",j,mean_div_left); */
  /*       PhyML_Fprintf(fp,"%4d 1 %f\n",j,mean_div_rght); */

  /*       mean_div_left = .0; */
  /*       for(k=0;k<ns;k++) mean_div_left +=
   * tree->a_edges[j]->div_post_pred_left[k]; */

  /*       mean_div_rght = .0; */
  /*       for(k=0;k<ns;k++)  */
  /* 	{ */
  /* 	  mean_div_rght += tree->a_edges[j]->div_post_pred_rght[k]; */
  /* 	} */

  /*       if((mean_div_left != tree->data->init_len) || (mean_div_rght !=
   * tree->data->init_len)) */
  /* 	{ */
  /* 	  PhyML_Printf("\n. mean_div_left = %f mean_div_rght = %f init_len =
   * %d", */
  /* 		 mean_div_left,mean_div_rght,tree->data->init_len); */
  /* 	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__); */
  /* 	  Warn_And_Exit(""); */
  /* 	} */
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Print_Diversity_Pre(t_node *a, t_node *d, t_edge *b, FILE *fp,
                         t_tree *tree)
{
  int k, ns;

  ns = -1;

  if (d->tax)
    return;
  else
  {

    if (tree->io->datatype == NT)
      ns = 4;
    else if (tree->io->datatype == AA)
      ns = 20;

    if (d == b->left)
      for (k = 0; k < ns; k++)
        PhyML_Fprintf(fp, "%4d 0 %2d %4d\n", b->num, k,
                      b->div_post_pred_left[k]);
    else
      for (k = 0; k < ns; k++)
        PhyML_Fprintf(fp, "%4d 1 %2d %4d\n", b->num, k,
                      b->div_post_pred_rght[k]);

    for (k = 0; k < 3; k++)
      if (d->v[k] != a) Print_Diversity_Pre(d, d->v[k], d->b[k], fp, tree);
  }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Print_Tip_Partials(t_tree *tree, t_node *d)
{
  if (!d->tax)
  {
    fprintf(stdout, "Node %d is not a Taxa\n", d->num);
    fflush(stdout);
    return;
  }

  assert(d->b[0]->rght == d);
  assert(d->b[0]->rght->tax);
  fprintf(stdout, "Taxa/Node %d\n", d->num);
  int site, i;
  for (site = 0; site < tree->data->n_pattern; ++site)
  {
    fprintf(stdout, "[%d: ", site);
    for (i = 0; i < tree->mod->ns; ++i)
    {
      int tip_partial = d->b[0]->p_lk_tip_r[site * tree->mod->ns + i];
      fprintf(stdout, "%d", tip_partial);
      fflush(stdout);
    }
    fprintf(stdout, "] ");
    fflush(stdout);
  }
  fprintf(stdout, "\n");
  fflush(stdout);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Print_Edge_PMats(t_tree *tree, t_edge *b)
{
  phydbl *Pij;
#ifdef BEAGLE
  Pij = (phydbl *)malloc(tree->mod->ns * tree->mod->ns *
                         tree->mod->ras->n_catg * sizeof(phydbl));
  if (NULL == Pij) Warn_And_Exit(__PRETTY_FUNCTION__);
  int ret = beagleGetTransitionMatrix(tree->b_inst, b->Pij_rr_idx, Pij);
  if (ret < 0)
  {
    fprintf(stderr, "beagleGetTransitionMatrix() on instance %i failed:%i\n\n",
            tree->b_inst, ret);
    Free(Pij);
    Exit("");
  }
#else
  Pij = b->Pij_rr;
#endif
  fprintf(stdout,
          "\nflattened P-Matrices (for each rate category) "
          "state*state*num_rates[%d*%d*%d] for branch num:%i\n",
          tree->mod->ns, tree->mod->ns, tree->mod->ras->n_catg, b->num);
  int i;
  for (i = 0; i < tree->mod->ns * tree->mod->ns * tree->mod->ras->n_catg; ++i)
  {
    fprintf(stdout, "%f,", Pij[i]);
    fflush(stdout);
  }
  fprintf(stdout, "\n");
  fflush(stdout);
#ifdef BEAGLE
  Free(Pij);
#endif
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Print_All_Edge_PMats(t_tree *tree)
{
  int i;
  for (i = 0; i < 2 * tree->n_otu - 3; ++i)
    Print_Edge_PMats(tree, tree->a_edges[i]);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Print_Edge_Likelihoods(t_tree *tree, t_edge *b,
                            bool scientific /*Print in scientific notation?*/)
{
  int   catg, site, j;
  char *fmt = scientific
                  ? "[%d,%d,%d]%e "
                  : "[%d,%d,%d]%f "; // rate category, site, state, likelihood

  phydbl *lk_left  = b->p_lk_left;
  phydbl *lk_right = b->p_lk_rght;

  fprintf(stdout, "\n");
  fflush(stdout);
  if (NULL != lk_left) // not a tip?
  {
    fprintf(
        stdout,
        "Partial Likelihoods on LEFT subtree of Branch %d [rate,site,state]:\n",
        b->num);
    for (catg = 0; catg < tree->mod->ras->n_catg; ++catg)
      for (site = 0; site < tree->data->n_pattern; ++site)
        for (j = 0; j < tree->mod->ns; ++j)
#ifdef BEAGLE
          fprintf(stdout, fmt, catg, site, j,
                  lk_left[catg * tree->n_pattern * tree->mod->ns +
                          site * tree->mod->ns + j]);
#else
          fprintf(stdout, fmt, catg, site, j,
                  lk_left[catg * tree->mod->ns +
                          site * tree->mod->ras->n_catg * tree->mod->ns + j]);
#endif
    fflush(stdout);
  }
  else // is a tip
  {
    fprintf(stdout, "Likelihoods on LEFT tip of Branch %d [site,state]:\n",
            b->num);
    for (site = 0; site < tree->data->n_pattern; ++site)
      for (j = 0; j < tree->mod->ns; ++j)
        fprintf(stdout, "[%d,%d]%.1f ", site, j,
                b->p_lk_tip_l[site * tree->mod->ns + j]);
    fflush(stdout);
  }
  fprintf(stdout, "\n");
  fflush(stdout);
  if (NULL != lk_right) // not a tip?
  {
    fprintf(stdout,
            "Partial Likelihoods on RIGHT subtree of Branch %d "
            "[rate,site,state]:\n",
            b->num);
    for (catg = 0; catg < tree->mod->ras->n_catg; ++catg)
      for (site = 0; site < tree->data->n_pattern; ++site)
        for (j = 0; j < tree->mod->ns; ++j)
#ifdef BEAGLE
          fprintf(stdout, fmt, catg, site, j,
                  lk_right[catg * tree->n_pattern * tree->mod->ns +
                           site * tree->mod->ns + j]);
#else
          fprintf(stdout, fmt, catg, site, j,
                  lk_right[catg * tree->mod->ns +
                           site * tree->mod->ras->n_catg * tree->mod->ns + j]);
#endif
    fflush(stdout);
  }
  else // is a tip
  {
    fprintf(stdout, "Likelihoods on RIGHT tip of Branch %d [site,state]:\n",
            b->num);
    for (site = 0; site < tree->data->n_pattern; ++site)
      for (j = 0; j < tree->mod->ns; ++j)
        fprintf(stdout, "[%d,%d]%.1f ", site, j,
                b->p_lk_tip_r[site * tree->mod->ns + j]);
    fflush(stdout);
  }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Print_All_Edge_Likelihoods(t_tree *tree)
{
  int i;
  for (i = 0; i < 2 * tree->n_otu - 3; ++i)
    Print_Edge_Likelihoods(tree, tree->a_edges[i], false);
  fflush(NULL);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Dump_Arr_S(short int *arr, int num)
{
  int i;

  if (NULL == arr)
  {
    PhyML_Fprintf(stderr, "\n. Err. in file %s at line %d (function '%s')\n",
                  __FILE__, __LINE__, __FUNCTION__);
    Exit("\n. Trying to print NULL array");
    return;
  }
  fprintf(stdout, "[");
  for (i = 0; i < num; ++i)
  {
    fprintf(stdout, "%d,", arr[i]);
    fflush(stdout);
  }
  fprintf(stdout, "]\n");
  fflush(stdout);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Dump_Arr_D(phydbl *arr, int num)
{
  int i;

  if (NULL == arr)
  {
    PhyML_Fprintf(stderr, "\n. Err. in file %s at line %d (function '%s')\n",
                  __FILE__, __LINE__, __FUNCTION__);
    Exit("\n. Trying to print NULL array");
    return;
  }
  fprintf(stdout, "[");
  for (i = 0; i < num; ++i)
  {
    fprintf(stdout, "%g,", arr[i]);
    fflush(stdout);
  }
  fprintf(stdout, "]\n");
  fflush(stdout);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Dump_Arr_I(int *arr, int num)
{
  int i;

  if (NULL == arr)
  {
    PhyML_Fprintf(stderr, "\n. Err. in file %s at line %d (function '%s')\n",
                  __FILE__, __LINE__, __FUNCTION__);
    Exit("\n. Trying to print NULL array");
    return;
  }
  fprintf(stdout, "[");
  for (i = 0; i < num; ++i)
  {
    fprintf(stdout, "%d,", arr[i]);
    fflush(stdout);
  }
  fprintf(stdout, "]\n");
  fflush(stdout);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Print_Tree_Structure(t_tree *tree)
{
  int i;
  PhyML_Fprintf(stdout, "\n. n_otu: %d", tree->n_otu);

  for (i = 0; i < 2 * tree->n_otu - 1; ++i)
  {
    if (tree->a_edges[i])
      PhyML_Fprintf(
          stdout,
          "\n. Edge %p %3d, Length: %f LeftNode %3d [%s], RightNode %3d [%s]",
          tree->a_edges[i], tree->a_edges[i]->num, tree->a_edges[i]->l->v,
          tree->a_edges[i]->left ? tree->a_edges[i]->left->num : -1,
          tree->a_edges[i]->left
              ? tree->a_edges[i]->left->tax ? tree->a_edges[i]->left->name : ""
              : "",
          tree->a_edges[i]->rght ? tree->a_edges[i]->rght->num : -1,
          tree->a_edges[i]->left
              ? tree->a_edges[i]->rght->tax ? tree->a_edges[i]->rght->name : ""
              : "");
    else
      PhyML_Fprintf(stdout, "\n. NULL");
  }

  for (i = 0; i < 2 * tree->n_otu - 1; ++i)
  {
    if (tree->a_nodes[i])
      PhyML_Fprintf(
          stdout,
          "\n. Node %p %3d v0: %3d v1: %3d v2: %3d b0: %3d b1: %3d b2: %3d",
          tree->a_nodes[i], tree->a_nodes[i]->num,
          tree->a_nodes[i]->v[0] ? tree->a_nodes[i]->v[0]->num : -1,
          tree->a_nodes[i]->v[1] ? tree->a_nodes[i]->v[1]->num : -1,
          tree->a_nodes[i]->v[2] ? tree->a_nodes[i]->v[2]->num : -1,
          tree->a_nodes[i]->v[0] ? tree->a_nodes[i]->b[0]->num : -1,
          tree->a_nodes[i]->v[1] ? tree->a_nodes[i]->b[1]->num : -1,
          tree->a_nodes[i]->v[2] ? tree->a_nodes[i]->b[2]->num : -1);
    else
      PhyML_Fprintf(stdout, "\n. NULL");
  }
}
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

t_tree *Read_User_Tree(calign *cdata, t_mod *mod, option *io)
{
  t_tree *tree;

  PhyML_Printf("\n. Reading tree...");
  fflush(NULL);
  if (io->n_trees == 1) rewind(io->fp_in_tree);

  Detect_Tree_File_Format(io);

  tree = NULL;

  switch (io->tree_file_format)
  {
  case PHYLIP:
  {
    tree = Read_Tree_File_Phylip(io->fp_in_tree);
    break;
  }
  case NEXUS:
  {
    io->nex_com_list = Make_Nexus_Com();
    Init_Nexus_Format(io->nex_com_list, io->fp_in_tree);
    Get_Nexus_Data(io);
    tree = io->tree;
    Free_Nexus(io);
    break;
  }

  default:
  {
    PhyML_Fprintf(stderr, "\n. Problem detected with tree file format");
    Exit("\n");
    break;
  }
  }

  if (tree == NULL) Exit("\n. Input tree not found...");
  /* Add branch lengths if necessary */
  if (tree->has_branch_lengths == NO)
  {
    assert(cdata != NULL);
    assert(mod != NULL);
    Add_BioNJ_Branch_Lengths(tree, cdata, mod, NULL);
  }
  return tree;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Print_Time_Info(time_t t_beg, time_t t_end)
{
  div_t hour, min;

  hour = div(t_end - t_beg, 3600);
  min  = div(t_end - t_beg, 60);
  min.quot -= hour.quot * 60;

  PhyML_Printf("\n\n. Time used %dh%dm%ds\n", hour.quot, min.quot,
               (int)(t_end - t_beg) % 60);
  PhyML_Printf("\nooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo"
               "ooooooooooooooooooooooooooooooooooo\n");
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Print_Time_Info_Brief(time_t t_beg, time_t t_end)
{
  div_t hour, min;

  hour = div(t_end - t_beg, 3600);
  min  = div(t_end - t_beg, 60);
  min.quot -= hour.quot * 60;

  PhyML_Printf("\n. Time used %dh%dm%ds\n", hour.quot, min.quot,
               (int)(t_end - t_beg) % 60);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void PhyML_Printf(char *format, ...)
{
  va_list ptr;

#ifdef MPI
  if (Global_myRank == 0)
  {
    va_start(ptr, format);
    vprintf(format, ptr);
    va_end(ptr);
  }
#else
  va_start(ptr, format);
  vprintf(format, ptr);
  va_end(ptr);
#endif

  fflush(NULL);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void PhyML_Fprintf(FILE *fp, char *format, ...)
{
  va_list ptr;

#ifdef MPI
  if (Global_myRank == 0)
  {
    va_start(ptr, format);
    vfprintf(fp, format, ptr);
    va_end(ptr);
  }
#else
  va_start(ptr, format);
  vfprintf(fp, format, ptr);
  va_end(ptr);
#endif

  fflush(NULL);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int PhyML_Fscanf(FILE *fp, char *format, ...)
{
  va_list ptr;
  int     rv;

  rv = -1;

#ifdef MPI
  if (Global_myRank == 0)
  {
    va_start(ptr, format);
    rv = vfscanf(fp, format, ptr);
    if (!rv) Generic_Exit(__FILE__, __LINE__, __FUNCTION__);
    va_end(ptr);
  }
#else
  va_start(ptr, format);
  rv = vfscanf(fp, format, ptr);
  if (!rv) Generic_Exit(__FILE__, __LINE__, __FUNCTION__);
  va_end(ptr);
#endif

  return (rv);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Read_Clade_Priors(char *file_name, t_tree *tree)
{
  FILE  *fp;
  char  *s, *line;
  int    n_clade_priors;
  int    clade_size;
  char **clade_list;
  int    i, pos;
  phydbl prior_low, prior_up;
  int    node_num;

  PhyML_Printf("\n");
  PhyML_Printf("\n. Reading prior on node ages.\n");

  line = (char *)mCalloc(T_MAX_LINE, sizeof(char));
  s    = (char *)mCalloc(T_MAX_LINE, sizeof(char));

  clade_list = (char **)mCalloc(tree->n_otu, sizeof(char *));
  for (i = 0; i < tree->n_otu; i++)
    clade_list[i] = (char *)mCalloc(T_MAX_NAME, sizeof(char));

  fp = Openfile(file_name, 0);

  n_clade_priors = 0;
  do
  {
    if (!fgets(line, T_MAX_LINE, fp)) break;

    clade_size = 0;
    pos        = 0;
    do
    {
      i = 0;

      while (line[pos] == ' ') pos++;

      while ((line[pos] != ' ') && (line[pos] != '\n') && line[pos] != '#')
      {
        s[i] = line[pos];
        i++;
        pos++;
      }
      s[i] = '\0';

      /* PhyML_Printf("\n. s = %s\n",s); */

      if (line[pos] == '\n' || line[pos] == '#') break;
      pos++;

      if (strcmp(s, "|"))
      {
        strcpy(clade_list[clade_size], s);
        clade_size++;
      }
      else
        break;
    } while (1);

    if (line[pos] != '#' && line[pos] != '\n')
    {
      phydbl val1, val2;
      /* 	  sscanf(line+pos,"%lf %lf",&prior_up,&prior_low); */
      sscanf(line + pos, "%lf %lf", &val1, &val2);
      prior_up  = (phydbl)val1;
      prior_low = (phydbl)val2;
      node_num  = -1;
      if (!strcmp("@root@", clade_list[0]))
        node_num = tree->n_root->num;
      else
        node_num = Find_Clade(clade_list, clade_size, tree);

      n_clade_priors++;

      if (node_num < 0)
      {
        PhyML_Printf("\n");
        PhyML_Printf("\n");
        PhyML_Printf("\n. "
                     "........................................................."
                     "........");
        PhyML_Printf("\n. WARNING: could not find any clade in the tree "
                     "referred to with the following taxon names:");
        for (i = 0; i < clade_size; i++)
          PhyML_Printf("\n. \"%s\"", clade_list[i]);
        PhyML_Printf("\n. "
                     "........................................................."
                     "........");
      }
      else
      {
        tree->times->t_has_prior[node_num] = YES;

        tree->times->t_prior_min[node_num] = -MAX(prior_low, prior_up);
        tree->times->t_prior_max[node_num] = -MIN(prior_low, prior_up);

        if (fabs(prior_low - prior_up) < 1.E-6 &&
            tree->a_nodes[node_num]->tax == YES)
          tree->times->nd_t[node_num] = prior_low;

        PhyML_Printf("\n");
        PhyML_Printf("\n. "
                     "[%3d]...................................................."
                     "..............",
                     n_clade_priors);
        PhyML_Printf("\n. Node %4d matches the clade referred to with the "
                     "following taxon names:",
                     node_num);
        for (i = 0; i < clade_size; i++)
          PhyML_Printf("\n. - \"%s\"", clade_list[i]);
        PhyML_Printf("\n. Lower bound set to: %15f time units.",
                     MIN(prior_low, prior_up));
        PhyML_Printf("\n. Upper bound set to: %15f time units.",
                     MAX(prior_low, prior_up));
        PhyML_Printf("\n. "
                     "........................................................."
                     "..............");
      }
    }
  } while (1);

  PhyML_Printf("\n. Read prior information on %d %s.", n_clade_priors,
               n_clade_priors > 1 ? "clades" : "clade");

  if (!n_clade_priors)
  {
    PhyML_Fprintf(stderr, "\n. PhyTime could not find any prior on node age.");
    PhyML_Fprintf(stderr,
                  "\n. This is likely due to a problem in the calibration ");
    PhyML_Fprintf(stderr,
                  "\n. file format. Make sure, for instance, that there is ");
    PhyML_Fprintf(stderr,
                  "\n. a blank character between the end of the last name");
    PhyML_Fprintf(stderr,
                  "\n. of each clade and the character `|'. Otherwise, ");
    PhyML_Fprintf(stderr, "\n. please refer to the example file.\n");
    PhyML_Fprintf(stderr, "\n. Err. in file %s at line %d (function '%s')\n",
                  __FILE__, __LINE__, __FUNCTION__);
    Warn_And_Exit("");
  }

  for (i = 0; i < tree->n_otu; i++) Free(clade_list[i]);
  Free(clade_list);
  Free(line);
  Free(s);
  fclose(fp);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

option *Get_Input(int argc, char **argv)
{
  option *io;
  t_mod  *mod;
  t_opt  *s_opt;
  int     rv;

  rv = 1;

  io    = (option *)Make_Input();
  mod   = (t_mod *)Make_Model_Basic();
  s_opt = (t_opt *)Make_Optimiz();

  Set_Defaults_Input(io);
  Set_Defaults_Model(mod);
  Set_Defaults_Optimiz(s_opt);

  io->mod    = mod;
  mod->io    = io;
  mod->s_opt = s_opt;

#ifdef MPI
  rv = Read_Command_Line(io, argc, argv);
#elif (defined PHYTIME || defined INVITEE || defined PHYREX || defined TEST)
  rv = Read_Command_Line(io, argc, argv);
#else
  switch (argc)
  {
  case 1:
  {
    Launch_Interface(io);
    break;
  }
  default:
  {
    rv = Read_Command_Line(io, argc, argv);
  }
  }
#endif

  if (rv)
    return io;
  else
    return NULL;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int Set_Whichmodel(int select)
{
  int wm;

  wm = -1;

  switch (select)
  {
  case 1:
  {
    wm = JC69;
    break;
  }
  case 2:
  {
    wm = K80;
    break;
  }
  case 3:
  {
    wm = F81;
    break;
  }
  case 4:
  {
    wm = HKY85;
    break;
  }
  case 5:
  {
    wm = F84;
    break;
  }
  case 6:
  {
    wm = TN93;
    break;
  }
  case 7:
  {
    wm = GTR;
    break;
  }
  case 8:
  {
    wm = CUSTOM;
    break;
  }
  case 11:
  {
    wm = WAG;
    break;
  }
  case 12:
  {
    wm = DAYHOFF;
    break;
  }
  case 13:
  {
    wm = JTT;
    break;
  }
  case 14:
  {
    wm = BLOSUM62;
    break;
  }
  case 15:
  {
    wm = MTREV;
    break;
  }
  case 16:
  {
    wm = RTREV;
    break;
  }
  case 17:
  {
    wm = CPREV;
    break;
  }
  case 18:
  {
    wm = DCMUT;
    break;
  }
  case 19:
  {
    wm = VT;
    break;
  }
  case 20:
  {
    wm = MTMAM;
    break;
  }
  case 21:
  {
    wm = MTART;
    break;
  }
  case 22:
  {
    wm = HIVW;
    break;
  }
  case 23:
  {
    wm = HIVB;
    break;
  }
  case 24:
  {
    wm = FLU;
    break;
  }
  case 25:
  {
    wm = CUSTOMAA;
    break;
  }
  case 26:
  {
    wm = LG;
    break;
  }
  case 27:
  {
    wm = AB;
    break;
  }
  default:
  {
    PhyML_Fprintf(
        stderr, "\n. Model number %d is unknown. Please use a valid model name",
        select);
    Exit("\n");
    break;
  }
  }

  return wm;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Print_Data_Structure(int final, FILE *fp, t_tree *mixt_tree)
{
  int     n_partition_elem;
  char   *s;
  t_tree *tree, *cpy_mixt_tree;
  int     c, cc, cc_efrq, cc_rmat, cc_lens;
  char   *param;
  int    *link_efrq, *link_rmat, *link_lens;
  phydbl  r_mat_weight_sum, e_frq_weight_sum;

  PhyML_Fprintf(fp, "\n. Random seed: %d",mixt_tree->io->r_seed);

  PhyML_Fprintf(fp, "\n\n. Starting tree: %s",
                mixt_tree->io->in_tree == 2 ? mixt_tree->io->in_tree_file
                                            : "BioNJ");

  cpy_mixt_tree = mixt_tree;

  n_partition_elem = 1;
  tree             = mixt_tree;
  do
  {
    tree = tree->next_mixt;
    if (!tree) break;
    n_partition_elem++;
  } while (1);

  s    = (char *)mCalloc(2, sizeof(char));
  s[0] = ' ';
  s[1] = '\0';
  tree = mixt_tree;
  do
  {
    s = (char *)mRealloc(
        s, (int)(strlen(s) + strlen(tree->io->in_align_file) + 2 + 2),
        sizeof(char));
    strcat(s, tree->io->in_align_file);
    strcat(s, ", ");
    tree = tree->next_mixt;
    if (!tree) break;
  } while (1);
  s[(int)strlen(s) - 2] = ' ';
  s[(int)strlen(s) - 1] = '\0';

  if (final == NO)
    PhyML_Fprintf(fp, "\n\n. Processing %d data %s (%s)", n_partition_elem,
                  n_partition_elem > 1 ? "sets" : "set", s);
  if (final == YES)
    PhyML_Fprintf(fp, "\n\n. Processed %d data %s (%s)", n_partition_elem,
                  n_partition_elem > 1 ? "sets" : "set", s);
  Free(s);

  if (final == YES)
    PhyML_Fprintf(fp, "\n\n. Final log-likelihood: %f", mixt_tree->c_lnL);

  r_mat_weight_sum =
      MIXT_Get_Sum_Chained_Scalar_Dbl(mixt_tree->next->mod->r_mat_weight);
  e_frq_weight_sum =
      MIXT_Get_Sum_Chained_Scalar_Dbl(mixt_tree->next->mod->e_frq_weight);

  do
  {
    int class = 0;

    PhyML_Fprintf(fp, "\n\n");
    PhyML_Fprintf(fp, "\n "
                      "________________________________________________________"
                      "_______________ ");
    PhyML_Fprintf(fp, "\n|                                                     "
                      "                  |");
    PhyML_Fprintf(fp, "\n| %40s      (partition element %2d)  |",
                  mixt_tree->io->in_align_file, mixt_tree->dp);
    PhyML_Fprintf(fp, "\n|_____________________________________________________"
                      "__________________|");
    PhyML_Fprintf(fp, "\n");

    PhyML_Fprintf(fp, "\n. Number of rate classes:\t\t%20d",
                  mixt_tree->mod->ras->n_catg +
                      (mixt_tree->mod->ras->invar ? 1 : 0));
    if (mixt_tree->mod->ras->n_catg > 1)
    {
      PhyML_Fprintf(fp, "\n. Model of rate variation:\t\t%20s",
                    mixt_tree->mod->ras->free_mixt_rates ? "FreeRates"
                    : mixt_tree->mod->ras->invar         ? "Gamma+Inv"
                                                         : "Gamma");
      if (mixt_tree->mod->ras->free_mixt_rates == NO)
      {
        PhyML_Fprintf(fp, "\n. Gamma shape parameter value:\t\t%20.2f",
                      mixt_tree->mod->ras->alpha->v);
        PhyML_Fprintf(fp, "\n   Optimise: \t\t\t\t%20s",
                      mixt_tree->mod->ras->alpha->optimize == YES ? "yes"
                                                                  : "no");
      }
      if (mixt_tree->mod->ras->invar == YES)
      {
        PhyML_Fprintf(fp, "\n. Proportion of invariable sites:\t%20.2f",
                      mixt_tree->mod->ras->pinvar->v);
        PhyML_Fprintf(fp, "\n   Optimise: \t\t\t\t%20s",
                      mixt_tree->mod->ras->pinvar->optimize == YES ? "yes"
                                                                   : "no");
      }
    }
    PhyML_Fprintf(fp, "\n. Relative average rate:\t\t%20f",
                  mixt_tree->mod->br_len_mult->v);

    tree = mixt_tree;
    do
    {
      if (tree->is_mixt_tree) tree = tree->next;

      PhyML_Fprintf(fp, "\n");
      PhyML_Fprintf(fp, "\n. Mixture class %d", class + 1);

      if (mixt_tree->mod->ras->n_catg > 1)
      {
        if (tree->mod->ras->invar == NO)
        {
          PhyML_Fprintf(fp, "\n   Relative substitution rate:\t%20f",
                        mixt_tree->mod->ras->gamma_rr
                            ->v[tree->mod->ras->parent_class_number]);
          PhyML_Fprintf(fp, "\n   Rel. rate freq. (> 0 rates):\t%20f",
                        mixt_tree->mod->ras->gamma_r_proba
                            ->v[tree->mod->ras->parent_class_number]);
          PhyML_Fprintf(fp, "\n   Rate class number:\t\t%20d",
                        tree->mod->ras->parent_class_number);
        }
        else
        {
          PhyML_Fprintf(fp, "\n   Relative substitution rate:\t%20f", 0.0);
          PhyML_Fprintf(fp, "\n   Relative rate freq.:\t\t%20f",
                        mixt_tree->mod->ras->pinvar->v);
        }
      }

      PhyML_Fprintf(fp, "\n   Substitution model:\t\t%20s",
                    tree->mod->modelname->s);

      if (tree->mod->whichmodel == CUSTOM)
        PhyML_Fprintf(fp, "\n   Substitution model code:\t%20s",
                      tree->mod->custom_mod_string->s);

      if (tree->mod->whichmodel == CUSTOMAA)
        PhyML_Fprintf(fp, "\n   Rate matrix file name:\t%20s",
                      tree->mod->aa_rate_mat_file->s);

      if (tree->mod->whichmodel == K80 || tree->mod->whichmodel == HKY85 ||
          tree->mod->whichmodel == TN93)
      {
        PhyML_Fprintf(fp, "\n   Value of the ts/tv ratio:\t%20f",
                      tree->mod->kappa->v);
        PhyML_Fprintf(fp, "\n   Optimise ts/tv ratio:\t%20s",
                      tree->mod->kappa->optimize ? "yes" : "no");
      }
      else if (tree->mod->whichmodel == GTR || tree->mod->whichmodel == CUSTOM)
      {
        PhyML_Fprintf(fp, "\n   Optimise subst. rates:\t%20s",
                      tree->mod->r_mat->optimize ? "yes" : "no");
        if (final == YES)
        {
          PhyML_Fprintf(fp, "\n   Subst. rate A<->C:\t\t%20.2f",
                        tree->mod->r_mat->rr->v[0]);
          PhyML_Fprintf(fp, "\n   Subst. rate A<->G:\t\t%20.2f",
                        tree->mod->r_mat->rr->v[1]);
          PhyML_Fprintf(fp, "\n   Subst. rate A<->T:\t\t%20.2f",
                        tree->mod->r_mat->rr->v[2]);
          PhyML_Fprintf(fp, "\n   Subst. rate C<->G:\t\t%20.2f",
                        tree->mod->r_mat->rr->v[3]);
          PhyML_Fprintf(fp, "\n   Subst. rate C<->T:\t\t%20.2f",
                        tree->mod->r_mat->rr->v[4]);
          PhyML_Fprintf(fp, "\n   Subst. rate G<->T:\t\t%20.2f",
                        tree->mod->r_mat->rr->v[5]);
        }
      }

      PhyML_Fprintf(fp, "\n   Rate matrix weight:\t\t%20f",
                    tree->mod->r_mat_weight->v / r_mat_weight_sum);

      if (tree->io->datatype == NT && tree->mod->whichmodel != JC69 &&
          tree->mod->whichmodel != K80)
      {
        PhyML_Fprintf(fp, "\n   Optimise nucleotide freq.:\t%20s",
                      tree->mod->e_frq->type == ML ? "yes" : "no");
        if (final == YES)
        {
          PhyML_Fprintf(fp, "\n   Freq(A):\t\t\t%20.2f",
                        tree->mod->e_frq->pi->v[0]);
          PhyML_Fprintf(fp, "\n   Freq(C):\t\t\t%20.2f",
                        tree->mod->e_frq->pi->v[1]);
          PhyML_Fprintf(fp, "\n   Freq(G):\t\t\t%20.2f",
                        tree->mod->e_frq->pi->v[2]);
          PhyML_Fprintf(fp, "\n   Freq(T):\t\t\t%20.2f",
                        tree->mod->e_frq->pi->v[3]);
        }
      }
      else if (tree->io->datatype == AA)
      {
        char *s;

        s = (char *)mCalloc(50, sizeof(char));

        if (tree->mod->e_frq->type == EMPIRICAL)
        {
          strcpy(s, "Empirical");
        }
        else if (tree->mod->e_frq->type == MODEL)
        {
          strcpy(s, "Model");
        }
        else if (tree->mod->e_frq->type == ML)
        {
          strcpy(s, "Optimized");
        }

        PhyML_Fprintf(fp, "\n   Amino-acid freq.:\t\t%20s", s);

        Free(s);
      }

      PhyML_Fprintf(fp, "\n   Equ. freq. weight:\t\t%20f",
                    tree->mod->e_frq_weight->v / e_frq_weight_sum);

      PhyML_Fprintf(fp, "\n   Integrated length (IL) model:%20s",
                    (tree->mod->gamma_mgf_bl == YES) ? "yes" : "no");

      if (tree->mod->gamma_mgf_bl == YES)
        PhyML_Fprintf(fp, "\n   Variance of IL model:\t%20.2g",
                      tree->mod->l_var_sigma->v);

      

      class ++;

      tree = tree->next;

      if (tree && tree->is_mixt_tree == YES) break;
    } while (tree);

    PhyML_Fprintf(fp, "\n");
    AIC(mixt_tree);
    BIC(mixt_tree);
    if (mixt_tree->mod->aic->print == YES)
      PhyML_Fprintf(fp, "\n   AIC: %20.2f", mixt_tree->mod->aic->v);
    if (mixt_tree->mod->bic->print == YES)
      PhyML_Fprintf(fp, "\n   BIC: %20.2f", mixt_tree->mod->bic->v);

    mixt_tree = mixt_tree->next_mixt;
    if (!mixt_tree) break;
  } while (1);

  tree = cpy_mixt_tree;
  c    = 0;
  do
  {
    if (tree->is_mixt_tree) tree = tree->next;
    c++;
    tree = tree->next;
  } while (tree);

  link_efrq = (int *)mCalloc(c, sizeof(int));
  link_lens = (int *)mCalloc(c, sizeof(int));
  link_rmat = (int *)mCalloc(c, sizeof(int));

  PhyML_Fprintf(fp, "\n");
  PhyML_Fprintf(fp, "\n");
  PhyML_Fprintf(fp, "\n "
                    "__________________________________________________________"
                    "_____________ ");
  PhyML_Fprintf(fp, "\n|                                                       "
                    "                |");
  PhyML_Fprintf(fp, "\n|                        Model summary table            "
                    "                |");
  PhyML_Fprintf(fp, "\n|_______________________________________________________"
                    "________________|");
  PhyML_Fprintf(fp, "\n");
  PhyML_Fprintf(fp, "\n");
  /* PhyML_Fprintf(fp,"                  "); */
  PhyML_Fprintf(fp, "  ------------------");
  tree = cpy_mixt_tree;
  c    = 0;
  do
  {
    if (tree->is_mixt_tree) tree = tree->next;
    PhyML_Fprintf(fp, "---");
    tree = tree->next;
  } while (tree);

  param = (char *)mCalloc(30, sizeof(char));
  PhyML_Fprintf(fp, "\n");
  strcpy(param, "Partition element ");
  PhyML_Fprintf(fp, "  %18s", param);
  tree = cpy_mixt_tree;
  c    = 0;
  do
  {
    if (tree->is_mixt_tree) tree = tree->next;
    PhyML_Fprintf(fp, "%2d ", tree->mixt_tree->dp);
    tree = tree->next;
  } while (tree);

  PhyML_Fprintf(fp, "\n");
  PhyML_Fprintf(fp, "  ------------------");
  tree = cpy_mixt_tree;
  c    = 0;
  do
  {
    if (tree->is_mixt_tree) tree = tree->next;
    PhyML_Fprintf(fp, "---");
    tree = tree->next;
  } while (tree);

  tree = cpy_mixt_tree;
  c    = 0;
  do
  {
    if (tree->is_mixt_tree) tree = tree->next;
    link_rmat[c] = -1;
    link_lens[c] = -1;
    link_efrq[c] = -1;
    tree         = tree->next;
    c++;
  } while (tree);

  mixt_tree = cpy_mixt_tree;
  cc_efrq   = 97;
  cc_rmat   = 97;
  cc_lens   = 97;
  cc        = 0;
  do
  {
    if (mixt_tree->is_mixt_tree) mixt_tree = mixt_tree->next;

    if (link_efrq[cc] < 0)
    {
      link_efrq[cc] = cc_efrq;
      tree          = mixt_tree->next;
      c             = cc + 1;
      if (tree)
      {
        do
        {
          if (tree->is_mixt_tree) tree = tree->next;

          if (mixt_tree->mod->e_frq == tree->mod->e_frq) link_efrq[c] = cc_efrq;

          tree = tree->next;
          c++;
        } while (tree);
      }
      cc_efrq++;
    }

    if (link_lens[cc] < 0)
    {
      link_lens[cc] = cc_lens;
      tree          = mixt_tree->next;
      c             = cc + 1;
      if (tree)
      {
        do
        {
          if (tree->is_mixt_tree) tree = tree->next;

          if (mixt_tree->a_edges[0]->l == tree->a_edges[0]->l)
            link_lens[c] = cc_lens;

          tree = tree->next;
          c++;
        } while (tree);
      }
      cc_lens++;
    }

    if (link_rmat[cc] < 0)
    {
      link_rmat[cc] = cc_rmat;
      tree          = mixt_tree->next;
      c             = cc + 1;
      if (tree)
      {
        do
        {
          if (tree->is_mixt_tree) tree = tree->next;

          if (mixt_tree->mod->r_mat == tree->mod->r_mat &&
              mixt_tree->mod->whichmodel == tree->mod->whichmodel &&
              !strcmp(mixt_tree->mod->custom_mod_string->s,
                      tree->mod->custom_mod_string->s) &&
              !strcmp(mixt_tree->mod->aa_rate_mat_file->s,
                      tree->mod->aa_rate_mat_file->s))
            link_rmat[c] = cc_rmat;

          tree = tree->next;
          c++;
        } while (tree);
      }
      cc_rmat++;
    }

    cc++;
    mixt_tree = mixt_tree->next;
  } while (mixt_tree);

  PhyML_Fprintf(fp, "\n");
  strcpy(param, "State frequencies ");
  PhyML_Fprintf(fp, "  %18s", param);

  tree = cpy_mixt_tree;
  c    = 0;
  do
  {
    if (tree->is_mixt_tree) tree = tree->next;
    PhyML_Fprintf(fp, "%2c ", link_efrq[c]);
    tree = tree->next;
    c++;
  } while (tree);

  PhyML_Fprintf(fp, "\n");
  strcpy(param, "Branch lengths ");
  PhyML_Fprintf(fp, "  %18s", param);

  tree = cpy_mixt_tree;
  c    = 0;
  do
  {
    if (tree->is_mixt_tree) tree = tree->next;
    PhyML_Fprintf(fp, "%2c ", link_lens[c]);
    tree = tree->next;
    c++;
  } while (tree);

  PhyML_Fprintf(fp, "\n");
  strcpy(param, "Rate matrix ");
  PhyML_Fprintf(fp, "  %18s", param);

  tree = cpy_mixt_tree;
  c    = 0;
  do
  {
    if (tree->is_mixt_tree) tree = tree->next;
    PhyML_Fprintf(fp, "%2c ", link_rmat[c]);
    tree = tree->next;
    c++;
  } while (tree);

  PhyML_Fprintf(fp, "\n");
  PhyML_Fprintf(fp, "  ------------------");
  tree = cpy_mixt_tree;
  c    = 0;
  do
  {
    if (tree->is_mixt_tree) tree = tree->next;
    PhyML_Fprintf(fp, "---");
    tree = tree->next;
  } while (tree);
  PhyML_Fprintf(fp, "\n");

  if (final == YES)
  {
    tree = cpy_mixt_tree;
    c    = 0;
    do
    {
      PhyML_Fprintf(fp, "\n");
      PhyML_Fprintf(fp, "\n. Tree estimated from data partition %d", c++);
      Br_Len_Involving_Invar(tree->next);
      Rescale_Br_Len_Multiplier_Tree(tree->next);
      s = Write_Tree(
          tree->next); /*! tree->next is not a mixt_tree so edge lengths
                            are not averaged over when writing the tree out. */
      PhyML_Fprintf(fp, "\n %s", s);
      Br_Len_Not_Involving_Invar(tree->next);
      Unscale_Br_Len_Multiplier_Tree(tree->next);
      Free(s);
      tree = tree->next_mixt;
    } while (tree);
  }

  Free(param);
  Free(link_efrq);
  Free(link_rmat);
  Free(link_lens);

  mixt_tree = cpy_mixt_tree;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

option *PhyML_XML(char *xml_filename)
{
  char     *most_likely_tree;
  phydbl    best_lnL;
  int       num_rand_tree;
  t_tree   *mixt_tree, *tree;
  option   *io, *dum;
  xml_node *xml_root;

  mixt_tree        = XML_Process_Base(xml_filename);
  io               = mixt_tree->io;
  most_likely_tree = NULL;
  xml_root         = mixt_tree->xml_root;
  best_lnL         = UNLIKELY;
  dum              = NULL;

  dum          = (option *)mCalloc(1, sizeof(option));
  dum->use_xml = YES;
  
  for (num_rand_tree = 0; num_rand_tree < io->mod->s_opt->n_rand_starts;
       num_rand_tree++)
  {
    MIXT_Check_Model_Validity(mixt_tree);
    MIXT_Chain_Models(mixt_tree);
    Init_Model(mixt_tree->mod->io->cdata, mixt_tree->mod, mixt_tree->mod->io);
    Set_Model_Parameters(mixt_tree->mod);

    Print_Data_Structure(NO, stdout, mixt_tree);
    tree = MIXT_Starting_Tree(mixt_tree);
    Copy_Tree(tree, mixt_tree);
    Free_Tree(tree);


    if (mixt_tree->io->mod->s_opt->random_input_tree)
    {
      PhyML_Printf("\n\n. [%3d/%3d]", num_rand_tree + 1,
                   mixt_tree->io->mod->s_opt->n_rand_starts);
      Random_Tree(mixt_tree);
    }

    Copy_Tree(mixt_tree, mixt_tree->next);
    Connect_CSeqs_To_Nodes(mixt_tree->mod->io->cdata, mixt_tree->mod->io,
                           mixt_tree);
    Init_T_Beg(mixt_tree);

    Make_Tree_For_Pars(mixt_tree);
    Make_Tree_For_Lk(mixt_tree);
    Make_Spr(mixt_tree);

    MIXT_Chain_All(mixt_tree);
    MIXT_Check_Edge_Lens_In_All_Elem(mixt_tree);
    MIXT_Turn_Branches_OnOff_In_All_Elem(ON, mixt_tree);
    MIXT_Check_Invar_Struct_In_Each_Partition_Elem(mixt_tree);
    MIXT_Check_RAS_Struct_In_Each_Partition_Elem(mixt_tree);
    Br_Len_Not_Involving_Invar(mixt_tree);
    Unscale_Br_Len_Multiplier_Tree(mixt_tree);
    Set_Both_Sides(YES, mixt_tree);

    Set_Update_Eigen(YES, mixt_tree->mod);
    Lk(NULL, mixt_tree);
    Set_Update_Eigen(NO, mixt_tree->mod);

    if (mixt_tree->mod->s_opt->opt_topo)
    {
      Global_Spr_Search(mixt_tree);
    }
    else
    {
      Round_Optimize(mixt_tree, ROUND_MAX);
    }

    PhyML_Printf("\n\n. Log-likelihood = %f", mixt_tree->c_lnL);

    if ((num_rand_tree == io->mod->s_opt->n_rand_starts - 1) &&
        (io->mod->s_opt->random_input_tree))
    {
      num_rand_tree--;
      io->mod->s_opt->random_input_tree = NO;
    }

    Br_Len_Involving_Invar(mixt_tree);
    Rescale_Br_Len_Multiplier_Tree(mixt_tree);

    /*! Print the tree estimated using the current random (or BioNJ) starting
     * tree */
    if (io->mod->s_opt->n_rand_starts > 1)
    {
      Print_Tree(io->fp_out_trees, mixt_tree);
      fflush(NULL);
    }

    if (mixt_tree->c_lnL > best_lnL)
    {
      best_lnL = mixt_tree->c_lnL;
      if (most_likely_tree) Free(most_likely_tree);
      if (io->ratio_test) 
      {
        aLRT(mixt_tree);
        Collect_Edge_Support_Values(mixt_tree);
      }
      
      most_likely_tree     = Write_Tree(mixt_tree);
      mixt_tree->lock_topo = NO;
    }

    tree = mixt_tree;
    do
    {
      tree->mod->aic->print = YES;
      tree->mod->bic->print = YES;
      tree                  = tree->next_mixt;
    } while (tree);

    tree = mixt_tree;
    do
    {
      if (tree->io->print_site_lnl == YES)
        Print_Site_Lk(tree, tree->io->fp_out_lk);
      tree = tree->next_mixt;
    } while (tree);

    MIXT_Init_T_End(mixt_tree);
    Print_Data_Structure(YES, mixt_tree->io->fp_out_stats, mixt_tree);

    MIXT_Cv(mixt_tree);

    Free_Best_Spr(mixt_tree);
    Free_Spr_List_One_Edge(mixt_tree);
    Free_Spr_List_All_Edge(mixt_tree);
    Free_Tree_Pars(mixt_tree);
    Free_Tree_Lk(mixt_tree);
  }

  /*! Print the most likely tree in the output file */
  if (!mixt_tree->io->quiet)
    PhyML_Printf("\n\n. Printing the most likely tree in file '%s'...\n",
                 Basename(mixt_tree->io->out_tree_file));
  PhyML_Fprintf(mixt_tree->io->fp_out_tree, "%s\n", most_likely_tree);

  /*! Bootstrap analysis */
  MIXT_Bootstrap(most_likely_tree, xml_root);

  while (io->prev != NULL) io = io->prev;

  Free(most_likely_tree);

  tree = mixt_tree;
  do
  {
    Free_Calign(tree->data);
    tree = tree->next_mixt;
  } while (tree);

  tree = mixt_tree;
  do
  {
    Free_Optimiz(tree->mod->s_opt);
    tree = tree->next;
  } while (tree);

  Free_Model_Complete(mixt_tree->mod);
  Free_Model_Basic(mixt_tree->mod);
  Free_Tree(mixt_tree);

  if (io->fp_out_trees) fclose(io->fp_out_trees);
  if (io->fp_out_tree) fclose(io->fp_out_tree);
  if (io->fp_out_stats) fclose(io->fp_out_stats);
  if (io->fp_out_json_trace) fclose(io->fp_out_json_trace);
  Free_Input(io);

  XML_Free_XML_Tree(xml_root);

  return (dum);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/*! Check that the same nodes in the different mixt_trees are
  connected to the same taxa
*/
void Check_Taxa_Sets(t_tree *mixt_tree)
{
  t_tree *tree;
  int     i;

  tree = mixt_tree;
  do
  {
    if (tree->next)
    {
      for (i = 0; i < tree->n_otu; i++)
      {
        if (strcmp(tree->a_nodes[i]->name, tree->next->a_nodes[i]->name))
        {
          PhyML_Fprintf(
              stderr,
              "\n. There seems to be a problem in one (or more) of your");
          PhyML_Fprintf(stderr,
                        "\n. sequence alignments. PhyML could not match taxon");
          PhyML_Fprintf(stderr,
                        "\n. '%s' found in file '%s' with any of the taxa",
                        tree->a_nodes[i]->name, tree->io->in_align_file);
          PhyML_Fprintf(stderr, "\n. listed in file '%s'.",
                        tree->next->io->in_align_file);
          Exit("\n");
        }
      }
    }
    tree = tree->next;
  } while (tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Make_Ratematrix_From_XML_Node(xml_node *instance, option *io, t_mod *mod)
{
  char *model = NULL;
  int   select;

  model = XML_Get_Attribute_Value(instance, "model");

  if (model == NULL)
  {
    PhyML_Fprintf(stderr, "\n. Poorly formated XML file.");
    PhyML_Fprintf(stderr,
                  "\n. Attribute 'model' is mandatory in a <ratematrix> node.");
    Exit("\n");
  }

  select = XML_Validate_Attr_Int(model, 27,
                                 "xxxxx",    // 0
                                 "JC69",     // 1
                                 "K80",      // 2
                                 "F81",      // 3
                                 "HKY85",    // 4
                                 "F84",      // 5
                                 "TN93",     // 6
                                 "GTR",      // 7
                                 "CUSTOM",   // 8
                                 "xxxxx",    // 9
                                 "xxxxx",    // 10
                                 "WAG",      // 11
                                 "DAYHOFF",  // 12
                                 "JTT",      // 13
                                 "BLOSUM62", // 14
                                 "MTREV",    // 15
                                 "RTREV",    // 16
                                 "CPREV",    // 17
                                 "DCMUT",    // 18
                                 "VT",       // 19
                                 "MTMAM",    // 20
                                 "MTART",    // 21
                                 "HIVW",     // 22
                                 "HIVB",     // 23
                                 "FLU",      // 24
                                 "CUSTOMAA", // 25
                                 "LG",       // 26
                                 "AB");      // 27

  if (select < 9)
  {
    if (io->datatype != NT)
    {
      PhyML_Fprintf(stderr,
                    "\n. Data type and selected model are incompatible");
      Exit("\n");
    }
  }
  else
  {
    if (io->datatype != AA)
    {
      PhyML_Fprintf(stderr,
                    "\n. Data type and selected model are incompatible");
      Exit("\n");
    }
  }

  mod->r_mat = (t_rmat *)Make_Rmat(mod->ns);
  Init_Rmat(mod->r_mat);
  Make_Custom_Model(mod);

  /*! Set model number & name */
  mod->whichmodel = Set_Whichmodel(select);
  Set_Model_Name(mod);

  if (mod->whichmodel == K80 || mod->whichmodel == HKY85 ||
      mod->whichmodel == TN93)
  {
    char *tstv, *opt_tstv;

    tstv = XML_Get_Attribute_Value(instance, "tstv");

    if (tstv)
    {
      mod->kappa->v        = String_To_Dbl(tstv);
      mod->kappa->optimize = NO;
    }
    else
    {
      mod->kappa->optimize = YES;
    }

    opt_tstv = XML_Get_Attribute_Value(instance, "optimise.tstv");

    if (opt_tstv)
    {
      if (!strcmp(opt_tstv, "true") || !strcmp(opt_tstv, "yes"))
      {
        mod->kappa->optimize        = YES;
        mod->s_opt->opt_subst_param = YES;
      }
      else
      {
        mod->kappa->optimize        = NO;
        mod->s_opt->opt_subst_param = NO;
      }
    }
  }
  else
  {
    mod->kappa->optimize = NO;
  }

  if (mod->whichmodel == GTR || mod->whichmodel == CUSTOM)
  {
    xml_node *rr;
    char     *val;
    phydbl    v;

    rr = XML_Search_Node_Name("rr", YES, instance);

    if (rr != NULL)
    {
      // A<->C
      val = XML_Get_Attribute_Value(rr, "AC");
      if (val == NULL)
      {
        PhyML_Printf("\n. Please specify the relative rate of substitution "
                     "between A and T");
        Generic_Exit(__FILE__, __LINE__, __FUNCTION__);
      }
      else
      {
        v = strtod(val, NULL);
        if (v != 0 && v > 0.0)
        {
          mod->r_mat->rr->v[0]     = v;
          mod->r_mat->rr_val->v[0] = log(v);
        }
        else
        {
          PhyML_Printf("\n. Invalid relative rate parameter value: '%s'.\n",
                       val);
          Generic_Exit(__FILE__, __LINE__, __FUNCTION__);
        }
      }

      // A<->G
      val = XML_Get_Attribute_Value(rr, "AG");
      if (val == NULL)
      {
        PhyML_Printf("\n. Please specify the relative rate of substitution "
                     "between A and T");
        Generic_Exit(__FILE__, __LINE__, __FUNCTION__);
      }
      else
      {
        v = strtod(val, NULL);
        if (v != 0 && v > 0.0)
        {
          mod->r_mat->rr->v[1]     = v;
          mod->r_mat->rr_val->v[1] = log(v);
        }
        else
        {
          PhyML_Printf("\n. Invalid relative rate parameter value: '%s'.\n",
                       val);
          Generic_Exit(__FILE__, __LINE__, __FUNCTION__);
        }
      }

      // A<->T
      val = XML_Get_Attribute_Value(rr, "AT");
      if (val == NULL)
      {
        PhyML_Printf("\n. Please specify the relative rate of substitution "
                     "between A and T");
        Generic_Exit(__FILE__, __LINE__, __FUNCTION__);
      }
      else
      {
        v = strtod(val, NULL);
        if (v != 0 && v > 0.0)
        {
          mod->r_mat->rr->v[2]     = v;
          mod->r_mat->rr_val->v[2] = log(v);
        }
        else
        {
          PhyML_Printf("\n. Invalid relative rate parameter value: '%s'.\n",
                       val);
          Generic_Exit(__FILE__, __LINE__, __FUNCTION__);
        }
      }

      // C<->G
      val = XML_Get_Attribute_Value(rr, "CG");
      if (val == NULL)
      {
        PhyML_Printf("\n. Please specify the relative rate of substitution "
                     "between A and T");
        Generic_Exit(__FILE__, __LINE__, __FUNCTION__);
      }
      else
      {
        v = strtod(val, NULL);
        if (v != 0 && v > 0.0)
        {
          mod->r_mat->rr->v[3]     = v;
          mod->r_mat->rr_val->v[3] = log(v);
        }
        else
        {
          PhyML_Printf("\n. Invalid relative rate parameter value: '%s'.\n",
                       val);
          Generic_Exit(__FILE__, __LINE__, __FUNCTION__);
        }
      }

      // C<->T
      val = XML_Get_Attribute_Value(rr, "CT");
      if (val == NULL)
      {
        PhyML_Printf("\n. Please specify the relative rate of substitution "
                     "between A and T");
        Generic_Exit(__FILE__, __LINE__, __FUNCTION__);
      }
      else
      {
        v = strtod(val, NULL);
        if (v != 0 && v > 0.0)
        {
          mod->r_mat->rr->v[4]     = v;
          mod->r_mat->rr_val->v[4] = log(v);
        }
        else
        {
          PhyML_Printf("\n. Invalid relative rate parameter value: '%s'.\n",
                       val);
          Generic_Exit(__FILE__, __LINE__, __FUNCTION__);
        }
      }

      // G<->T
      val = XML_Get_Attribute_Value(rr, "GT");

      if (val != NULL)
      {
        v = strtod(val, NULL);
        if (v != 0 && v > 0.0)
        {
          mod->r_mat->rr->v[5]     = v;
          mod->r_mat->rr_val->v[5] = log(v);
        }
        else
        {
          PhyML_Printf("\n. Invalid relative rate parameter value: '%s'.\n",
                       val);
          Generic_Exit(__FILE__, __LINE__, __FUNCTION__);
        }
      }
    }
  }

  char *opt_rr;

  opt_rr = XML_Get_Attribute_Value(instance, "optimise.rr");

  if (opt_rr != NULL)
  {
    assert(mod->r_mat != NULL);

    if (!strcmp(opt_rr, "yes") || !strcmp(opt_rr, "true"))
    {
      mod->r_mat->optimize        = YES;
      mod->s_opt->opt_subst_param = YES;
    }
    else
    {
      mod->r_mat->optimize        = NO;
      mod->s_opt->opt_subst_param = NO;
    }
  }

  /*! Custom model for nucleotide sequences. Read the corresponding
    code. */
  if (mod->whichmodel == CUSTOM)
  {
    char *model_code = NULL;

    model_code = XML_Get_Attribute_Value(instance, "model.code");

    if (!model_code)
    {
      PhyML_Fprintf(stderr,
                    "\n. No valid 'model.code' attribute could be found.\n");
      PhyML_Fprintf(stderr, "\n. Please fix your XML file.\n");
      Exit("\n");
    }
    else
    {
      strcpy(mod->custom_mod_string->s, model_code);
      Make_Custom_Model(mod);
      Translate_Custom_Mod_String(mod);
    }
  }

  /*! Custom model for amino-acids. Read in the rate matrix file */
  if (mod->whichmodel == CUSTOMAA)
  {
    char *r_mat_file;

    r_mat_file = XML_Get_Attribute_Value(instance, "ratematrix.file");

    if (!r_mat_file)
    {
      PhyML_Fprintf(stderr,
                    "\n. No valid 'ratematrix.file' attribute could be found.");
      PhyML_Fprintf(stderr, "\n. Please fix your XML file.\n");
      Exit("\n");
    }
    else
    {
      strcpy(mod->aa_rate_mat_file->s, r_mat_file);
    }

    /* Free(r_mat_file); */
  }

  char *buff;

  buff = XML_Get_Attribute_Value(instance->parent, "optimise.weights");
  if (buff && (!strcmp(buff, "yes") || !strcmp(buff, "true")))
  {
    mod->s_opt->opt_rmat_weight = YES;
  }
  else
  {
    mod->s_opt->opt_rmat_weight = NO;
  }

  char *il_model = NULL;
  il_model       = XML_Get_Attribute_Value(instance, "integrated.lens");

  if (il_model)
  {
    if (!strcmp(il_model, "yes") || !strcmp(il_model, "true"))
    {
      mod->gamma_mgf_bl = YES;
    }
    else
    {
      mod->gamma_mgf_bl = NO;
    }
  }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Make_Efrq_From_XML_Node(xml_node *instance, option *io, t_mod *mod)
{
  char *buff;

  mod->e_frq = (t_efrq *)Make_Efrq(mod->ns);
  Init_Efrq(NULL, mod->e_frq);

  buff = XML_Get_Attribute_Value(instance, "freqs");

  if (buff)
  {
    if (!strcmp(buff, "empirical"))
    {
      mod->e_frq->type         = EMPIRICAL;
      mod->e_frq->pi->optimize = NO;
    }
    if (!strcmp(buff, "model"))
    {
      mod->e_frq->type         = MODEL;
      mod->e_frq->pi->optimize = NO;
    }
    else if (!strcmp(buff, "optimized"))
    {
      mod->e_frq->type         = ML;
      mod->e_frq->pi->optimize = YES;
    }
    else if (!strcmp(buff, "user"))
    {
      mod->e_frq->type         = ML;
      mod->e_frq->pi->optimize = NO;
    }
  }

  buff = XML_Get_Attribute_Value(instance, "base.freqs");

  if (buff)
  {
    if (io->datatype == AA)
    {
      PhyML_Fprintf(
          stderr,
          "\n. Option 'base.freqs' is not allowed with amino-acid data.");
      Exit("\n");
    }
    else
    {
      phydbl A, C, G, T;

      if (mod->e_frq->type != USER)
      {
        PhyML_Fprintf(stderr, "\n. Please set 'freqs=user'.");
        Exit("\n");
      }

      sscanf(buff, "%lf,%lf,%lf,%lf", &A, &C, &G, &T);
      mod->e_frq->user_b_freq->v[0] = (phydbl)A;
      mod->e_frq->user_b_freq->v[1] = (phydbl)C;
      mod->e_frq->user_b_freq->v[2] = (phydbl)G;
      mod->e_frq->user_b_freq->v[3] = (phydbl)T;
      mod->e_frq->type              = USER;
    }
  }

  buff = XML_Get_Attribute_Value(instance->parent, "optimise.weights");
  if (buff && (!strcmp(buff, "yes") || !strcmp(buff, "true")))
  {
    mod->s_opt->opt_efrq_weight = YES;
  }
  else
  {
    mod->s_opt->opt_efrq_weight = NO;
  }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Make_Topology_From_XML_Node(xml_node *instance, option *io, t_mod *mod)
{
  // Starting tree
  char *init_tree;

  init_tree = XML_Get_Attribute_Value(instance, "init.tree");

  if (!init_tree)
  {
    PhyML_Fprintf(stderr,
                  "\n. Attribute 'init.tree=bionj|user|random' is mandatory");
    PhyML_Fprintf(stderr, "\n. Please amend your XML file accordingly.");
    Exit("\n");
  }

  if (!strcmp(init_tree, "user") || !strcmp(init_tree, "User"))
  {
    char *starting_tree = NULL;
    starting_tree       = XML_Get_Attribute_Value(instance, "file.name");

    if (!Filexists(starting_tree))
    {
      PhyML_Fprintf(stderr, "\n. The tree file '%s' could not be found.",
                    starting_tree);
      Exit("\n");
    }
    else
    {
      strcpy(io->in_tree_file, starting_tree);
      io->in_tree    = 2;
      io->fp_in_tree = Openfile(io->in_tree_file, 0);
    }
  }
  else if (!strcmp(init_tree, "random") || !strcmp(init_tree, "Random"))
  {
    char *n_rand_starts = NULL;

    io->mod->s_opt->random_input_tree = YES;

    n_rand_starts = XML_Get_Attribute_Value(instance, "n.rand.starts");

    if (n_rand_starts)
    {
      mod->s_opt->n_rand_starts = atoi(n_rand_starts);
      if (mod->s_opt->n_rand_starts < 1)
        Exit("\n. Number of random starting trees must be > 0.\n\n");
    }

    strcpy(io->out_trees_file, io->in_align_file);
    strcat(io->out_trees_file, "_phyml_rand_trees");
    if (io->append_run_ID)
    {
      strcat(io->out_trees_file, "_");
      strcat(io->out_trees_file, io->run_id_string);
    }
    strcat(io->out_trees_file, ".txt");
    io->fp_out_trees = Openfile(io->out_trees_file, 1);
  }
  else if (!strcmp(init_tree, "parsimony") || !strcmp(init_tree, "Parsimony"))
  {
    io->in_tree = 1;
  }

  // Estimate tree topology
  char *optimise = NULL;

  optimise = XML_Get_Attribute_Value(instance, "optimise.topology");

  if (optimise)
  {
    int select;

    select = XML_Validate_Attr_Int(optimise, 6, "true", "yes", "y", "false",
                                   "no", "n");

    if (select > 2)
      io->mod->s_opt->opt_topo = NO;
    else
    {
      char *search;
      int   select;

      search = XML_Get_Attribute_Value(instance, "search");

      if (search == NULL)
      {
        io->mod->s_opt->topo_search = SPR_MOVE;
        io->mod->s_opt->opt_topo    = YES;
      }
      else
      {
        select = XML_Validate_Attr_Int(search, 4, "spr", "nni", "best", "none");

        switch (select)
        {
        case 0:
        {
          io->mod->s_opt->topo_search = SPR_MOVE;
          io->mod->s_opt->opt_topo    = YES;
          break;
        }
        case 1:
        {
          io->mod->s_opt->topo_search = NNI_MOVE;
          io->mod->s_opt->opt_topo    = YES;
          break;
        }
        case 2:
        {
          io->mod->s_opt->topo_search = BEST_OF_NNI_AND_SPR;
          io->mod->s_opt->opt_topo    = YES;
          break;
        }
        case 3:
        {
          io->mod->s_opt->opt_topo = NO;
          break;
        }
        default:
        {
          PhyML_Fprintf(stderr, "\n. Topology search option '%s' is not valid.",
                        search);
          Exit("\n");
          break;
        }
        }
      }
    }
  }

  // Edge length unit
  char *edge_len;

  edge_len = XML_Get_Attribute_Value(instance, "edge.lengths");

  if (edge_len != NULL)
  {
    if (!strcmp(edge_len, "calendar") || !strcmp(edge_len, "cal"))
    {
      io->edge_len_unit = CALENDAR;
    }
    else if (!strcmp(edge_len, "mutations") ||
             !strcmp(edge_len, "substitutions") || !strcmp(edge_len, "mut") ||
             !strcmp(edge_len, "subst"))
    {
      io->edge_len_unit = SUBSTITUTIONS;
    }
  }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Make_RAS_From_XML_Node(xml_node *parent, t_mod *mod)
{
  xml_node *w;
  char     *family;
  int       select;

  family           = NULL;
  mod->ras->n_catg = 0;

  XML_Check_Siterates_Node(parent);

  w = XML_Search_Node_Name("weights", YES, parent);
  if (w)
  {
    family = XML_Get_Attribute_Value(w, "family");
    if (family == NULL)
      select = -1;
    else
      select =
          XML_Validate_Attr_Int(family, 3, "gamma", "gamma+inv", "freerates");
    switch (select)
    {
    case 0: // Gamma model
    {
      char *alpha, *alpha_opt;

      mod->ras->pinvar->optimize = NO;
      mod->ras->invar            = NO;

      alpha = XML_Get_Attribute_Value(w, "alpha");

      if (alpha)
      {
        if (!strcmp(alpha, "estimate") || !strcmp(alpha, "estimated") ||
            !strcmp(alpha, "optimise") || !strcmp(alpha, "optimised"))
        {
          mod->ras->alpha->optimize = YES;
        }
        else
        {
          mod->ras->alpha->optimize = NO;
          mod->ras->alpha->v        = String_To_Dbl(alpha);
        }
      }

      alpha_opt = XML_Get_Attribute_Value(w, "optimise.alpha");

      if (alpha_opt)
      {
        if (!strcmp(alpha_opt, "yes") || !strcmp(alpha_opt, "true"))
        {
          mod->ras->alpha->optimize = YES;
        }
        else
        {
          mod->ras->alpha->optimize = NO;
        }
      }

      mod->ras->n_catg = XML_Get_Number_Of_Classes_Siterates(parent);

      Make_RAS_Complete(mod->ras);

      break;
    }
    case 1: // Gamma+Inv model
    {
      char *alpha, *pinv, *alpha_opt, *pinv_opt;

      mod->ras->invar            = YES;
      mod->ras->pinvar->optimize = YES;

      alpha = XML_Get_Attribute_Value(w, "alpha");

      if (alpha)
      {
        if (!strcmp(alpha, "estimate") || !strcmp(alpha, "estimated") ||
            !strcmp(alpha, "optimise") || !strcmp(alpha, "optimised"))
        {
          mod->ras->alpha->optimize = YES;
        }
        else
        {
          mod->ras->alpha->optimize = NO;
          mod->ras->alpha->v        = String_To_Dbl(alpha);
          ;
        }
      }

      alpha_opt = XML_Get_Attribute_Value(w, "optimise.alpha");

      if (alpha_opt)
      {
        if (!strcmp(alpha_opt, "yes") || !strcmp(alpha_opt, "true"))
        {
          mod->ras->alpha->optimize = YES;
        }
        else
        {
          mod->ras->alpha->optimize = NO;
        }
      }

      pinv = XML_Get_Attribute_Value(w, "pinv");

      if (pinv)
      {
        if (!strcmp(pinv, "estimate") || !strcmp(pinv, "estimated") ||
            !strcmp(pinv, "optimise") || !strcmp(pinv, "optimised"))
        {
          mod->ras->pinvar->optimize = YES;
        }
        else
        {
          mod->ras->pinvar->optimize = NO;
          mod->ras->pinvar->v        = String_To_Dbl(pinv);
          ;
        }
      }

      pinv_opt = XML_Get_Attribute_Value(w, "optimise.pinv");

      if (pinv_opt)
      {
        if (!strcmp(pinv_opt, "yes") || !strcmp(pinv_opt, "true"))
        {
          mod->ras->pinvar->optimize = YES;
        }
        else
        {
          mod->ras->pinvar->optimize = NO;
        }
      }

      mod->ras->n_catg = XML_Get_Number_Of_Classes_Siterates(parent);
      break;
    }
    case 2: // FreeRate model
    {
      char *opt_freerates;

      mod->ras->free_mixt_rates = YES;

      mod->s_opt->opt_free_mixt_rates = YES;

      opt_freerates = XML_Get_Attribute_Value(w, "optimise.freerates");

      if (opt_freerates)
      {
        if (!strcmp(opt_freerates, "yes") || !strcmp(opt_freerates, "true"))
        {
          mod->s_opt->opt_free_mixt_rates = YES;
        }
        else
        {
          mod->s_opt->opt_free_mixt_rates = NO;
        }
      }

      mod->ras->n_catg = XML_Get_Number_Of_Classes_Siterates(parent);
      break;
    }
    default:
    {
      if (family != NULL)
        PhyML_Printf("\n. family: %s", family);
      else
      {
        PhyML_Printf("\n. Please specify a model family (\"gamma\", "
                     "\"gamma+inv\" or \"freerate\")");
        PhyML_Printf("\n. for every 'weights' element in every 'siterate' "
                     "element. Note that ");
        PhyML_Printf("\n. a \"gamma\" or a \"freerate\" model with a single "
                     "rate class amounts");
        PhyML_Printf("\n. to no variation of rates across sites, if this is "
                     "the model you'd");
        PhyML_Printf("\n. like to implement...");
      }

      PhyML_Printf("\n. Err. in file %s at line %d\n", __FILE__, __LINE__);
      Exit("\n");
    }
    }
  }

  int nc = XML_Get_Number_Of_Classes_Siterates(parent);

  if (w && nc != mod->ras->n_catg)
  {
    PhyML_Printf("\n. <siterates> component '%s'. The number of classes ",
                 parent->id);
    PhyML_Printf("\n. specified in the <weight> component should be ");
    PhyML_Printf("\n. the same as that of <instance> nodes. Please fix");
    PhyML_Printf("\n. your XML file accordingly.");
    Exit("\n");
  }

  if (!w) mod->ras->n_catg = nc;

  Make_RAS_Complete(mod->ras);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void Generic_Exit(const char *file, int line, const char *function)
{
  PhyML_Fprintf(stderr, "\n. Err. in file '%s' (line %d)", file, line);
  if (function != NULL) PhyML_Printf(", function '%s'", function);

#if defined(PHYTIME)
  PhyML_Fprintf(stderr, "\n. PhyTime finished prematurely.");
#elif defined(PHYREX)
  PhyML_Fprintf(stderr, "\n. PhyREX finished prematurely.");
#elif defined(PHYML)
  PhyML_Fprintf(stderr, "\n. PhyML finished prematurely.");
#else
  PhyML_Fprintf(stderr, "\n. The execution finished prematurely.");
#endif

  Exit("\n");
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void JSON_Write_Object(json_o *obj, FILE *where)
{
  json_kv *kv;

  assert(obj);
  assert(where);

  kv = obj->kv;
  assert(kv);

  PhyML_Fprintf(where, "{");

  do
  {
    PhyML_Fprintf(where, "\"%s\":", kv->key);
    if (kv->value != NULL)
      PhyML_Fprintf(where, "\"%s\"", kv->value);
    else if (kv->array != NULL)
      JSON_Write_Array(kv->array, where);
    else if (kv->object != NULL)
      JSON_Write_Object(kv->object, where);
    kv = kv->next;
    if (kv) PhyML_Fprintf(where, ",");
  } while (kv);

  PhyML_Fprintf(where, "}");
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void JSON_Write_Array(json_a *array, FILE *where)
{
  json_o *o;

  assert(where);
  assert(array);

  o = array->object;

  assert(o);

  PhyML_Fprintf(where, "[");

  do
  {
    JSON_Write_Object(o, where);
    o = o->next;
    if (o) PhyML_Fprintf(where, ",");
  } while (o);

  PhyML_Fprintf(where, "]\n");
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void JSON_Write_All(json_a *array, FILE *where)
{
  JSON_Write_Array(array, where);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

json_o *JSON_Tree_To_Object(t_tree *mixt_tree)
{
  t_tree  *tree;
  json_kv *state, *kv;
  json_o  *ret, *o;
  json_a  *a;
  int      string_len;
  char    *s;
  time_t   t_end;

  string_len = 20;

  ret       = (json_o *)mCalloc(1, sizeof(json_o));
  ret->next = NULL;

  state   = (json_kv *)mCalloc(1, sizeof(json_kv));
  ret->kv = state;

  state->key = (char *)mCalloc(string_len, sizeof(char));
  strcpy(state->key, "state");

  state->value  = NULL;
  state->array  = NULL;
  state->object = NULL;
  state->next   = NULL;

  state->array = (json_a *)mCalloc(1, sizeof(json_a));
  a            = state->array;

  a->object = (json_o *)mCalloc(1, sizeof(json_o));
  o         = a->object;

  // State num
  o->kv      = (json_kv *)mCalloc(1, sizeof(json_kv));
  kv         = o->kv;
  kv->value  = NULL;
  kv->array  = NULL;
  kv->object = NULL;
  kv->next   = NULL;
  kv->key    = (char *)mCalloc(string_len, sizeof(char));
  strcpy(kv->key, "type");
  kv->value = (char *)mCalloc(string_len, sizeof(char));
  strcpy(kv->value, "t_num");

  kv->next = (json_kv *)mCalloc(1, sizeof(json_kv));
  kv       = kv->next;

  kv->value  = NULL;
  kv->array  = NULL;
  kv->object = NULL;
  kv->next   = NULL;
  kv->key    = (char *)mCalloc(string_len, sizeof(char));
  strcpy(kv->key, "id");
  kv->value = (char *)mCalloc(string_len, sizeof(char));
  strcpy(kv->value, "state_num");

  kv->next = (json_kv *)mCalloc(1, sizeof(json_kv));
  kv       = kv->next;

  kv->value  = NULL;
  kv->array  = NULL;
  kv->object = NULL;
  kv->next   = NULL;
  kv->key    = (char *)mCalloc(string_len, sizeof(char));
  strcpy(kv->key, "value");
  kv->value = (char *)mCalloc(string_len, sizeof(char));
  sprintf(kv->value, "%d", mixt_tree->json_num);
  mixt_tree->json_num++;

  kv->next = NULL;

  // Time
  o->next = (json_o *)mCalloc(1, sizeof(json_o));
  o       = o->next;

  o->kv = (json_kv *)mCalloc(1, sizeof(json_kv));
  kv    = o->kv;

  kv->value  = NULL;
  kv->array  = NULL;
  kv->object = NULL;
  kv->next   = NULL;
  kv->key    = (char *)mCalloc(string_len, sizeof(char));
  strcpy(kv->key, "type");
  kv->value = (char *)mCalloc(string_len, sizeof(char));
  strcpy(kv->value, "t_time");

  kv->next   = (json_kv *)mCalloc(1, sizeof(json_kv));
  kv         = kv->next;
  kv->value  = NULL;
  kv->array  = NULL;
  kv->object = NULL;
  kv->next   = NULL;
  kv->key    = (char *)mCalloc(string_len, sizeof(char));
  strcpy(kv->key, "id");
  kv->value = (char *)mCalloc(string_len, sizeof(char));
  strcpy(kv->value, "time");

  kv->next   = (json_kv *)mCalloc(1, sizeof(json_kv));
  kv         = kv->next;
  kv->value  = NULL;
  kv->array  = NULL;
  kv->object = NULL;
  kv->next   = NULL;
  kv->key    = (char *)mCalloc(string_len, sizeof(char));
  strcpy(kv->key, "value");
  kv->value = (char *)mCalloc(string_len, sizeof(char));
  time(&t_end);
  sprintf(kv->value, "%d", (int)(t_end - mixt_tree->t_beg));

  kv->next = NULL;

  // Tree
  o->next = (json_o *)mCalloc(1, sizeof(json_o));
  o       = o->next;

  o->kv = (json_kv *)mCalloc(1, sizeof(json_kv));
  kv    = o->kv;

  kv->value  = NULL;
  kv->array  = NULL;
  kv->object = NULL;
  kv->next   = NULL;
  kv->key    = (char *)mCalloc(string_len, sizeof(char));
  strcpy(kv->key, "type");
  kv->value = (char *)mCalloc(string_len, sizeof(char));
  strcpy(kv->value, "t_tree");

  kv->next   = (json_kv *)mCalloc(1, sizeof(json_kv));
  kv         = kv->next;
  kv->value  = NULL;
  kv->array  = NULL;
  kv->object = NULL;
  kv->next   = NULL;
  kv->key    = (char *)mCalloc(string_len, sizeof(char));
  strcpy(kv->key, "id");
  kv->value = (char *)mCalloc(string_len, sizeof(char));
  strcpy(kv->value, "tree");

  kv->next   = (json_kv *)mCalloc(1, sizeof(json_kv));
  kv         = kv->next;
  kv->value  = NULL;
  kv->array  = NULL;
  kv->object = NULL;
  kv->next   = NULL;
  kv->key    = (char *)mCalloc(string_len, sizeof(char));
  strcpy(kv->key, "value");
  s         = Write_Tree(mixt_tree);
  kv->value = (char *)mCalloc((int)strlen(s) + 1, sizeof(char));
  sprintf(kv->value, "%s", s);
  Free(s);

  kv->next = NULL;

  // Likelihood
  o->next = (json_o *)mCalloc(1, sizeof(json_o));
  o       = o->next;

  o->kv = (json_kv *)mCalloc(1, sizeof(json_kv));
  kv    = o->kv;

  kv->value  = NULL;
  kv->array  = NULL;
  kv->object = NULL;
  o->next    = NULL;
  kv->key    = (char *)mCalloc(string_len, sizeof(char));
  strcpy(kv->key, "type");
  kv->value = (char *)mCalloc(string_len, sizeof(char));
  strcpy(kv->value, "t_param");

  kv->next   = (json_kv *)mCalloc(1, sizeof(json_kv));
  kv         = kv->next;
  kv->value  = NULL;
  kv->array  = NULL;
  kv->object = NULL;
  kv->next   = NULL;
  kv->key    = (char *)mCalloc(string_len, sizeof(char));
  strcpy(kv->key, "id");
  kv->value = (char *)mCalloc(string_len, sizeof(char));
  strcpy(kv->value, "likelihood");

  kv->next   = (json_kv *)mCalloc(1, sizeof(json_kv));
  kv         = kv->next;
  kv->value  = NULL;
  kv->array  = NULL;
  kv->object = NULL;
  kv->next   = NULL;
  kv->key    = (char *)mCalloc(string_len, sizeof(char));
  strcpy(kv->key, "value");
  kv->value = (char *)mCalloc(string_len, sizeof(char));
  sprintf(kv->value, "%G", mixt_tree->c_lnL);

  kv->next = NULL;

  // TsTv
  {
    int          n_tstv, i;
    scalar_dbl **tstv;

    n_tstv = 0;
    tstv   = NULL;
    tree   = mixt_tree;

    do
    {
      if (tree->is_mixt_tree == YES) tree = tree->next;

      for (i = 0; i < n_tstv; ++i)
        if (tree->mod->kappa == tstv[i]) break;

      if (i == n_tstv)
      {
        if (!tstv)
          tstv = (scalar_dbl **)mCalloc(1, sizeof(scalar_dbl *));
        else
          tstv =
              (scalar_dbl **)mRealloc(tstv, n_tstv + 1, sizeof(scalar_dbl *));
        tstv[n_tstv] = tree->mod->kappa;
        n_tstv++;

        if (tree->mod->kappa->optimize == YES)
        {
          o->next = (json_o *)mCalloc(1, sizeof(json_o));
          o       = o->next;

          o->kv = (json_kv *)mCalloc(1, sizeof(json_kv));
          kv    = o->kv;

          kv->value  = NULL;
          kv->array  = NULL;
          kv->object = NULL;
          o->next    = NULL;
          kv->key    = (char *)mCalloc(string_len, sizeof(char));
          strcpy(kv->key, "type");
          kv->value = (char *)mCalloc(string_len, sizeof(char));
          strcpy(kv->value, "t_param");

          kv->next   = (json_kv *)mCalloc(1, sizeof(json_kv));
          kv         = kv->next;
          kv->value  = NULL;
          kv->array  = NULL;
          kv->object = NULL;
          kv->next   = NULL;
          kv->key    = (char *)mCalloc(string_len, sizeof(char));
          strcpy(kv->key, "id");
          kv->value = (char *)mCalloc(string_len, sizeof(char));
          strcpy(kv->value, "tstv");
          sprintf(kv->value + strlen(kv->value), "%d", n_tstv);

          kv->next   = (json_kv *)mCalloc(1, sizeof(json_kv));
          kv         = kv->next;
          kv->value  = NULL;
          kv->array  = NULL;
          kv->object = NULL;
          kv->next   = NULL;
          kv->key    = (char *)mCalloc(string_len, sizeof(char));
          strcpy(kv->key, "value");
          kv->value = (char *)mCalloc(string_len, sizeof(char));
          sprintf(kv->value, "%G", tree->mod->kappa->v);

          kv->next = NULL;
        }
      }
      tree = tree->next;
    } while (tree);

    if (tstv) Free(tstv);
  }

  // Alpha
  {
    int          n_alpha, i;
    scalar_dbl **alpha;

    n_alpha = 0;
    alpha   = NULL;
    tree    = mixt_tree;

    do
    {

      for (i = 0; i < n_alpha; i++)
        if (tree->mod->ras->alpha == alpha[i]) break;

      if (i == n_alpha)
      {
        if (!alpha)
          alpha = (scalar_dbl **)mCalloc(1, sizeof(scalar_dbl *));
        else
          alpha =
              (scalar_dbl **)mRealloc(alpha, n_alpha + 1, sizeof(scalar_dbl *));
        alpha[n_alpha] = tree->mod->ras->alpha;
        n_alpha++;

        if (tree->mod->ras->alpha->optimize == YES)
        {
          o->next = (json_o *)mCalloc(1, sizeof(json_o));
          o       = o->next;

          o->kv = (json_kv *)mCalloc(1, sizeof(json_kv));
          kv    = o->kv;

          kv->value  = NULL;
          kv->array  = NULL;
          kv->object = NULL;
          o->next    = NULL;
          kv->key    = (char *)mCalloc(string_len, sizeof(char));
          strcpy(kv->key, "type");
          kv->value = (char *)mCalloc(string_len, sizeof(char));
          strcpy(kv->value, "t_param");

          kv->next   = (json_kv *)mCalloc(1, sizeof(json_kv));
          kv         = kv->next;
          kv->value  = NULL;
          kv->array  = NULL;
          kv->object = NULL;
          kv->next   = NULL;
          kv->key    = (char *)mCalloc(string_len, sizeof(char));
          strcpy(kv->key, "id");
          kv->value = (char *)mCalloc(string_len, sizeof(char));
          strcpy(kv->value, "alpha");
          sprintf(kv->value + strlen(kv->value), "%d", n_alpha);

          kv->next   = (json_kv *)mCalloc(1, sizeof(json_kv));
          kv         = kv->next;
          kv->value  = NULL;
          kv->array  = NULL;
          kv->object = NULL;
          kv->next   = NULL;
          kv->key    = (char *)mCalloc(string_len, sizeof(char));
          strcpy(kv->key, "value");
          kv->value = (char *)mCalloc(string_len, sizeof(char));
          sprintf(kv->value, "%G", tree->mod->ras->alpha->v);

          kv->next = NULL;
        }
      }
      tree = tree->next_mixt;
    } while (tree);

    if (alpha) Free(alpha);
  }

  // Tree size
  {
    tree         = mixt_tree;
    int tree_num = 1;
    do
    {
      o->next = (json_o *)mCalloc(1, sizeof(json_o));
      o       = o->next;

      o->kv = (json_kv *)mCalloc(1, sizeof(json_kv));
      kv    = o->kv;

      kv->value  = NULL;
      kv->array  = NULL;
      kv->object = NULL;
      o->next    = NULL;
      kv->key    = (char *)mCalloc(string_len, sizeof(char));
      strcpy(kv->key, "type");
      kv->value = (char *)mCalloc(string_len, sizeof(char));
      strcpy(kv->value, "t_param");

      kv->next   = (json_kv *)mCalloc(1, sizeof(json_kv));
      kv         = kv->next;
      kv->value  = NULL;
      kv->array  = NULL;
      kv->object = NULL;
      kv->next   = NULL;
      kv->key    = (char *)mCalloc(string_len, sizeof(char));
      strcpy(kv->key, "id");
      kv->value = (char *)mCalloc(string_len, sizeof(char));
      strcpy(kv->value, "tree_size");
      sprintf(kv->value + strlen(kv->value), "%d", tree_num);

      kv->next   = (json_kv *)mCalloc(1, sizeof(json_kv));
      kv         = kv->next;
      kv->value  = NULL;
      kv->array  = NULL;
      kv->object = NULL;
      kv->next   = NULL;
      kv->key    = (char *)mCalloc(string_len, sizeof(char));
      strcpy(kv->key, "value");
      kv->value = (char *)mCalloc(string_len, sizeof(char));
      sprintf(kv->value, "%G", Get_Tree_Size(tree));

      kv->next = NULL;

      tree = tree->next_mixt;
    } while (tree);
  }

  return (ret);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

json_o *JSON_Tree_To_Object_Light(t_tree *mixt_tree)
{
  t_tree  *tree;
  json_kv *state, *kv;
  json_o  *ret, *o;
  int      string_len;
  char    *s;
  time_t   t_end;

  string_len = 20;

  ret       = (json_o *)mCalloc(1, sizeof(json_o));
  ret->next = NULL;

  state   = (json_kv *)mCalloc(1, sizeof(json_kv));
  ret->kv = state;

  state->key = (char *)mCalloc(string_len, sizeof(char));
  strcpy(state->key, "state");

  state->value  = NULL;
  state->array  = NULL;
  state->object = NULL;
  state->next   = NULL;

  state->object = (json_o *)mCalloc(1, sizeof(json_o));
  o             = state->object;

  o->kv = (json_kv *)mCalloc(1, sizeof(json_kv));
  kv    = o->kv;

  // State num
  kv->value  = NULL;
  kv->array  = NULL;
  kv->object = NULL;
  kv->next   = NULL;
  kv->key    = (char *)mCalloc(string_len, sizeof(char));
  strcpy(kv->key, "state_num");
  kv->value = (char *)mCalloc(string_len, sizeof(char));
  sprintf(kv->value, "%d", mixt_tree->json_num);

  // Time
  kv->next   = (json_kv *)mCalloc(1, sizeof(json_kv));
  kv         = kv->next;
  kv->value  = NULL;
  kv->array  = NULL;
  kv->object = NULL;
  kv->next   = NULL;
  kv->key    = (char *)mCalloc(string_len, sizeof(char));
  strcpy(kv->key, "time");
  kv->value = (char *)mCalloc(string_len, sizeof(char));
  time(&t_end);
  sprintf(kv->value, "%d", (int)(t_end - mixt_tree->t_beg));

  // Tree
  kv->next = (json_kv *)mCalloc(1, sizeof(json_kv));
  kv       = kv->next;

  kv->value  = NULL;
  kv->array  = NULL;
  kv->object = NULL;
  kv->next   = NULL;
  kv->key    = (char *)mCalloc(string_len, sizeof(char));
  strcpy(kv->key, "tree");
  s         = Write_Tree(mixt_tree);
  kv->value = (char *)mCalloc((int)strlen(s) + 1, sizeof(char));
  sprintf(kv->value, "%s", s);
  Free(s);

  // Likelihood
  kv->next   = (json_kv *)mCalloc(1, sizeof(json_kv));
  kv         = kv->next;
  kv->value  = NULL;
  kv->array  = NULL;
  kv->object = NULL;
  o->next    = NULL;
  kv->key    = (char *)mCalloc(string_len, sizeof(char));
  strcpy(kv->key, "likelihood");
  kv->value = (char *)mCalloc(string_len, sizeof(char));
  sprintf(kv->value, "%G", mixt_tree->c_lnL);

  // TsTv
  {
    int          n_tstv, i;
    scalar_dbl **tstv;

    n_tstv = 0;
    tstv   = NULL;
    tree   = mixt_tree;

    do
    {
      if (tree->is_mixt_tree == YES) tree = tree->next;

      for (i = 0; i < n_tstv; ++i)
        if (tree->mod->kappa == tstv[i]) break;

      if (i == n_tstv)
      {
        if (!tstv)
          tstv = (scalar_dbl **)mCalloc(1, sizeof(scalar_dbl *));
        else
          tstv =
              (scalar_dbl **)mRealloc(tstv, n_tstv + 1, sizeof(scalar_dbl *));
        tstv[n_tstv] = tree->mod->kappa;
        n_tstv++;

        if (tree->mod->kappa->optimize == YES)
        {

          kv->next   = (json_kv *)mCalloc(1, sizeof(json_kv));
          kv         = kv->next;
          kv->value  = NULL;
          kv->array  = NULL;
          kv->object = NULL;
          o->next    = NULL;
          kv->key    = (char *)mCalloc(string_len, sizeof(char));
          strcpy(kv->key, "tstv");
          sprintf(kv->key + strlen(kv->key), "%d", n_tstv);
          kv->value = (char *)mCalloc(string_len, sizeof(char));
          sprintf(kv->value, "%G", tree->mod->kappa->v);
        }
      }
      tree = tree->next;
    } while (tree);

    if (tstv) Free(tstv);
  }

  // Alpha
  {
    int          n_alpha, i;
    scalar_dbl **alpha;

    n_alpha = 0;
    alpha   = NULL;
    tree    = mixt_tree;

    do
    {

      for (i = 0; i < n_alpha; i++)
        if (tree->mod->ras->alpha == alpha[i]) break;

      if (i == n_alpha)
      {
        if (!alpha)
          alpha = (scalar_dbl **)mCalloc(1, sizeof(scalar_dbl *));
        else
          alpha =
              (scalar_dbl **)mRealloc(alpha, n_alpha + 1, sizeof(scalar_dbl *));
        alpha[n_alpha] = tree->mod->ras->alpha;
        n_alpha++;

        if (tree->mod->ras->alpha->optimize == YES)
        {

          kv->next   = (json_kv *)mCalloc(1, sizeof(json_kv));
          kv         = kv->next;
          kv->value  = NULL;
          kv->array  = NULL;
          kv->object = NULL;
          o->next    = NULL;
          kv->key    = (char *)mCalloc(string_len, sizeof(char));
          strcpy(kv->key, "alpha");
          sprintf(kv->key + strlen(kv->key), "%d", n_alpha);
          kv->value = (char *)mCalloc(string_len, sizeof(char));
          sprintf(kv->value, "%G", tree->mod->ras->alpha->v);
        }
      }
      tree = tree->next_mixt;
    } while (tree);

    if (alpha) Free(alpha);
  }

  // Tree size
  {
    tree         = mixt_tree;
    int tree_num = 1;
    do
    {

      kv->next = (json_kv *)mCalloc(1, sizeof(json_kv));
      kv       = kv->next;

      kv->value  = NULL;
      kv->array  = NULL;
      kv->object = NULL;
      o->next    = NULL;
      kv->key    = (char *)mCalloc(string_len, sizeof(char));
      strcpy(kv->key, "tree_size");
      sprintf(kv->key + strlen(kv->key), "%d", tree_num);
      kv->value = (char *)mCalloc(string_len, sizeof(char));
      sprintf(kv->value, "%G", Get_Tree_Size(tree));

      tree = tree->next_mixt;
    } while (tree);
  }

  kv->next = NULL;

  return (ret);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void JSON_Tree_Io(t_tree *tree, FILE *where)
{
  // Append
  json_o *o;
  fpos_t  pos;
  char    c;

  fgetpos(where, &pos);

  rewind(where);
  c = fgetc(where);

  if (c != '[')
  {
    PhyML_Fprintf(where, "[");
  }
  else
  {
    fsetpos(where, &pos);
    fseek(where, -1, SEEK_CUR);
    PhyML_Fprintf(where, ",");
  }

  PhyML_Fprintf(where, "\n");
  /* o = JSON_Tree_To_Object(tree); */
  o = JSON_Tree_To_Object_Light(tree);
  JSON_Write_Object(o, where);
  JSON_Free_Object(o);
  PhyML_Fprintf(where, "]");
  fflush(where);

  // Print out latest tree + info
  /* json_o *o; */

  /* rewind(where); */
  /* /\* PhyML_Fprintf(where,"["); *\/ */
  /* o = JSON_Tree_To_Object(tree); */
  /* JSON_Write_Object(o,where); */
  /* JSON_Free_Object(o); */
  /* /\* PhyML_Fprintf(where,"]"); *\/ */
  /* fflush(where); */
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void Print_Lk_Given_Edge_Recurr(t_node *a, t_node *d, t_edge *b, t_tree *tree)
{
  PhyML_Printf("\n___ Edge %3d (left=%3d rght=%3d) lnL=%f", b->num,
               b->left->num, b->rght->num, Lk(b, tree));

  if (d->tax)
    return;
  else
  {
    int i;
    for (i = 0; i < 3; i++)
      if (d->v[i] != a) Print_Lk_Given_Edge_Recurr(d, d->v[i], d->b[i], tree);
  }
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void Collect_Edge_Support_Values(t_tree *tree)
{
  int i;

  for (i = 0; i < 2 * tree->n_otu - 3; ++i)
  {
    if (tree->io->do_boot == YES || tree->io->do_bayesboot == YES)
    {
      tree->a_edges[i]->support_val = (phydbl)tree->a_edges[i]->bip_score;
    }
    else if (tree->io->do_tbe == YES)
    {
      int pminus1 =
          MIN(tree->a_edges[i]->left->bip_size[tree->a_edges[i]->l_r],
              tree->a_edges[i]->rght->bip_size[tree->a_edges[i]->r_l]) -
          1;
      tree->a_edges[i]->support_val =
          1. - ((tree->a_edges[i]->tdist_score) /
                (tree->io->n_boot_replicates * 1.0)) /
                   pminus1;
    }
    else if (tree->io->do_alrt == YES)
    {
      tree->a_edges[i]->support_val = tree->a_edges[i]->ratio_test;
    }
  }
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

#if defined(PHYREX)
void PHYREX_Output_Tree_Structure(FILE *fp, t_tree *tree)
{
  char *s;
  s = PHYREX_Print_Tree_Structure(tree);
  PhyML_Fprintf(fp, "%s", s);
  Free(s);
}
#endif

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

#if defined(PHYREX)
char *PHYREX_Print_Tree_Structure(t_tree *tree)
{
  t_dsk *disk;
  char  *s, *buff;
  FILE  *fp;
  fpos_t pos;

  fp = tmpfile();
  assert(fp);

  buff = (char *)mCalloc(T_MAX_LINE, sizeof(char));

  disk = tree->young_disk;
  while (disk->prev != NULL) disk = disk->prev;

  assert(disk->ldsk);
  assert(disk->ldsk->n_next > 0);

  fgetpos(fp, &pos);
  PhyML_Fprintf(fp, "begin\n");
  fsetpos(fp, &pos);
  if (fgets(buff, T_MAX_LINE, fp) == NULL)
    Generic_Exit(__FILE__, __LINE__, __FUNCTION__);
  s = (char *)mCalloc((int)strlen(buff) + 1, sizeof(char));
  sprintf(s + strlen(s), "%s", buff);

  fgetpos(fp, &pos);
  PhyML_Fprintf(fp, "%d %d\n", tree->n_otu, tree->mmod->n_dim);
  fsetpos(fp, &pos);
  if (fgets(buff, T_MAX_LINE, fp) == NULL)
    Generic_Exit(__FILE__, __LINE__, __FUNCTION__);
  s = (char *)mRealloc(s, (int)(strlen(s) + strlen(buff) + 1), sizeof(char));
  sprintf(s + strlen(s), "%s", buff);

  do
  {
    fgetpos(fp, &pos);
    PhyML_Fprintf(fp, "%s ", disk->id);
    PhyML_Fprintf(fp, "%f ", disk->time);
    PhyML_Fprintf(fp, "%f %f ", disk->centr->lonlat[0], disk->centr->lonlat[1]);

    if (disk->ldsk != NULL)
    {
      s = (char *)mRealloc(s, (int)strlen(s) + 9, sizeof(char));
      PhyML_Fprintf(fp, "*%s ", disk->ldsk->coord->id);
      PhyML_Fprintf(fp, "%f %f ", disk->ldsk->coord->lonlat[0],
                    disk->ldsk->coord->lonlat[1]);
    }
    for (int i = 0; i < disk->n_ldsk_a; ++i)
    {
      PhyML_Fprintf(fp, "%s ", disk->ldsk_a[i]->coord->id);
    }
    PhyML_Fprintf(fp, "\n");

    fsetpos(fp, &pos);
    if (fgets(buff, T_MAX_LINE, fp) == NULL)
      Generic_Exit(__FILE__, __LINE__, __FUNCTION__);
    s = (char *)mRealloc(s, (int)(strlen(s) + strlen(buff) + 1), sizeof(char));
    sprintf(s + strlen(s), "%s", buff);

    disk = disk->next;
  } while (disk);

  fgetpos(fp, &pos);
  PhyML_Fprintf(fp, "end");
  fsetpos(fp, &pos);
  if (fgets(buff, T_MAX_LINE, fp) == NULL)
    Generic_Exit(__FILE__, __LINE__, __FUNCTION__);
  s = (char *)mRealloc(s, (int)(strlen(s) + strlen(buff) + 1), sizeof(char));
  sprintf(s + strlen(s), "%s", buff);

  fclose(fp);

  Free(buff);

  return (s);
}
#endif

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

#if defined(PHYREX)
void PHYREX_Check_Point(FILE *fp, t_tree *tree)
{
  xml_node *n;
  char     *s;

  if (XML_Search_Node_Attribute_Value("add", "true", NO, tree->xml_root) ==
      NULL)
  {
    XML_Add_Attribute(tree->xml_root, "add", "true");
  }

  n = XML_Search_Node_Name("slfv", YES, tree->xml_root);
  if (n == NULL)
  {
    n = XML_Add_Node(tree->xml_root, "slfv");
    XML_Set_Node_Id(n, "SLFV1");
  }

  s = PHYREX_Print_Tree_Structure(tree);
  XML_Set_Node_Value(n, s);
  XML_Update_XML_Struct_Given_Model_Params(tree);
  /* XML_Update_XML_Struct_Given_MCMC_Params(tree); */
  XML_Write_XML_Graph(fp, tree->xml_root);
}
#endif

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

#if defined(PHYREX)
void PHYREX_Input_Tree_Structure(FILE *fp)
{
  char   *tk, *s, delim[4] = "\n\r ";
  int     n_otu, n_dim, idx;
  t_dsk  *disk;
  t_ldsk *new_ldsk, *root_ldsk, *ldsk;

  s = (char *)mCalloc(T_MAX_LINE, sizeof(char));

  do
  {
    if (fgets(s, T_MAX_LINE, fp) == NULL)
      Generic_Exit(__FILE__, __LINE__, __FUNCTION__);
  } while (strstr(s, "begin") == NULL);

  if (fgets(s, T_MAX_LINE, fp) == NULL)
    Generic_Exit(__FILE__, __LINE__, __FUNCTION__);

  tk    = strtok(s, delim);
  n_otu = strtod(tk, NULL);

  tk    = strtok(NULL, delim);
  n_dim = strtod(tk, NULL);

  root_ldsk = NULL;
  new_ldsk  = NULL;
  while (fgets(s, T_MAX_LINE, fp) != NULL)
  {
    tk = strtok(s, delim);
    if (strcmp(tk, "end") && tk != NULL)
    {
      disk = PHYREX_Make_Disk_Event(n_dim, n_otu);
      strcpy(disk->id, tk);

      tk         = strtok(NULL, delim);
      disk->time = strtold(tk, NULL);

      tk                     = strtok(NULL, delim);
      disk->centr->lonlat[0] = strtold(tk, NULL);
      tk                     = strtok(NULL, delim);
      disk->centr->lonlat[1] = strtold(tk, NULL);

      tk = strtok(NULL, delim);
      if (tk[0] == '*') // disk has ldsk on it
      {
        if (root_ldsk == NULL)
        {
          disk->ldsk       = PHYREX_Make_Lindisk_Node(n_dim);
          disk->ldsk->disk = disk;
          root_ldsk        = disk->ldsk;
        }
        else
        {
          disk->ldsk = PHYREX_Find_Ldsk_From_Id(tk + 1, root_ldsk);
          assert(disk->ldsk != NULL);
        }

        strcpy(disk->ldsk->coord->id, tk + 1);

        tk                           = strtok(NULL, delim);
        disk->ldsk->coord->lonlat[0] = strtold(tk, NULL);
        tk                           = strtok(NULL, delim);
        disk->ldsk->coord->lonlat[1] = strtold(tk, NULL);

        do
        {
          tk = strtok(NULL, delim);
          if (tk == NULL) break;
          new_ldsk = PHYREX_Make_Lindisk_Node(n_dim);
          strcpy(new_ldsk->coord->id, tk);
          PHYREX_Make_Lindisk_Next(disk->ldsk);
          disk->ldsk->next[disk->ldsk->n_next - 1] = new_ldsk;
          disk->ldsk_a[disk->ldsk->n_next - 1]     = new_ldsk;
          new_ldsk->prev                           = disk->ldsk;
        } while (tk);
      }
      else
      {
        idx = 0;
        do
        {
          tk = strtok(NULL, delim);
          if (tk == NULL) break;
          ldsk = PHYREX_Find_Ldsk_From_Id(tk, root_ldsk);
          assert(ldsk);
          disk->ldsk_a[idx] = ldsk;
          idx++;
        } while (tk);
      }
    }
  }

  Free(s);
}
#endif

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void Output_Scalar_Dbl(scalar_dbl *t, char *sep, FILE *fp)
{
  scalar_dbl *l;
  l = t;
  do
  {
    PhyML_Fprintf(fp, "%g%s", l->v, sep);
    l = l->next;
  } while (l);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

// s should look like that: "xxx={yyy},xxxx={yy},..."
t_label *Read_Labels(char *s)
{
  t_label  *lab, *init_lab;
  char     *s_cpy, *s_small;
  char     *label;
  char     *key, *val;
  int       cur, beg, op;
  short int finished;

  if (s[0] != '[' || s[(int)strlen(s) - 1] != ']')
  {
    PhyML_Fprintf(stderr,
                  "\n. Label is in wrong format. A proper label should");
    PhyML_Fprintf(stderr, "\n. look as follows: \"[xxx={yyy},xxxx={yy},...]\"");
    assert(FALSE);
  }

  lab      = Make_Label();
  init_lab = lab;

  s_cpy = (char *)mCalloc((int)strlen(s) + 1, sizeof(char));
  strcpy(s_cpy, s);

  label = (char *)mCalloc(strlen(s), sizeof(char));

  cur = beg = 0;
  finished  = NO;
  do
  {
    cur++;
    beg = cur;
    op  = 0;
    do
    {
      if (s_cpy[cur] == '{') op++;
      if (s_cpy[cur] == '}') op--;

      cur++;
    } while (!(s_cpy[cur] == ',' && op == 0) && s_cpy[cur] != ']');

    if (s_cpy[cur] == ']')
      finished = YES;
    else
      finished = NO;
    s_cpy[cur] = '\0';

    strcpy(label, s_cpy + beg);

    key = strtok_r(label, "=", &s_small);
    val = strtok_r(NULL, "=", &s_small);

    Free(lab->key);
    lab->key = (char *)mCalloc(strlen(key) + 1, sizeof(char));
    strcpy(lab->key, key);

    Free(lab->val);
    lab->val = (char *)mCalloc(strlen(val) + 1, sizeof(char));
    strcpy(lab->val, val);

    if (finished == NO)
    {
      lab->sep  = ',';
      lab->next = Make_Label();
      lab       = lab->next;
    }
  } while (finished == NO);

  lab = init_lab;

  Free(label);
  Free(s_cpy);

  return (lab);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void Print_Labels(FILE *fp_where, char *s_where, t_label *label)
{
  if (label == NULL)
    return;
  else
  {
    t_label *lab;

    lab = label;

    if (fp_where != NULL)
    {
      PhyML_Fprintf(fp_where, "[");
    }
    else
    {
      sprintf(s_where + (int)strlen(s_where), "[");
    }

    while (lab)
    {
      if (fp_where != NULL)
      {
        if (lab->prev == NULL)
          PhyML_Fprintf(fp_where, "&%s=%s", lab->key, lab->val);
        else
          PhyML_Fprintf(fp_where, "%s=%s", lab->key, lab->val);
        if (lab->next != NULL) PhyML_Fprintf(fp_where, ",");
      }
      else
      {
        if (lab->prev == NULL)
          sprintf(s_where + (int)strlen(s_where), "&%s=%s", lab->key, lab->val);
        else
          sprintf(s_where + (int)strlen(s_where), "%s=%s", lab->key, lab->val);
        if (lab->next != NULL) sprintf(s_where + (int)strlen(s_where), ",");
      }

      lab = lab->next;
    }

    if (fp_where != NULL)
    {
      PhyML_Fprintf(fp_where, "]");
    }
    else
    {
      sprintf(s_where + (int)strlen(s_where), "]");
    }
  }
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

#if (defined PHYREX)
void PHYREX_Print_MultiTypeTree_Config_File(int n_sites, char *filename,
                                            t_tree *tree)
{
  int   i, j, n_demes;
  char *s, **deme_names;
  FILE *fp;

  fp = Openfile(filename, WRITE);
  assert(fp);

  deme_names = (char **)mCalloc(n_sites, sizeof(char *));

  n_demes = 0;
  for (i = 0; i < tree->n_otu; i++)
  {
    s = strrchr(tree->a_nodes[i]->ldsk->coord->id, '_');
    for (j = 0; j < n_demes; j++)
      if (!strcmp(s + 1, deme_names[j])) break;
    if (j == n_demes)
    {
      deme_names[n_demes] = (char *)mCalloc(strlen(s + 1) + 1, sizeof(char));
      strcpy(deme_names[n_demes], s + 1);
      n_demes++;
    }
  }

  // n_demes is the number of non-empty sampling sites

  PhyML_Fprintf(fp,
                "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>");
  PhyML_Fprintf(
      fp, "\n<beast beautitemplate='MultiTypeTree' beautistatus='' "
          "namespace=\"beast.core:beast.evolution.alignment:beast.evolution."
          "tree.coalescent:beast.core.util:beast.evolution.nuc:beast.evolution."
          "operators:beast.evolution.sitemodel:beast.evolution."
          "substitutionmodel:beast.evolution.likelihood\" version=\"2.0\">");

  PhyML_Fprintf(fp, "\n<data id=\"data\" name=\"alignment\">");

  for (i = 0; i < tree->n_otu; i++)
  {
    PhyML_Fprintf(
        fp,
        "\n<sequence id=\"%s\" taxon=\"%s\" totalcount=\"4\" value=\"%s\"/>",
        tree->a_nodes[i]->ldsk->coord->id,
        /* tree->a_nodes[i]->coord->id, */
        tree->a_nodes[i]->name, tree->a_nodes[i]->c_seq->state);
  }

  PhyML_Fprintf(fp, "\n</data>");

  PhyML_Fprintf(
      fp, "\n<map name=\"Uniform\" >beast.math.distributions.Uniform</map>");
  PhyML_Fprintf(fp, "\n<map name=\"Exponential\" "
                    ">beast.math.distributions.Exponential</map>");
  PhyML_Fprintf(fp,
                "\n<map name=\"LogNormal\" "
                ">beast.math.distributions.LogNormalDistributionModel</map>");
  PhyML_Fprintf(
      fp, "\n<map name=\"Normal\" >beast.math.distributions.Normal</map>");
  PhyML_Fprintf(fp,
                "\n<map name=\"Beta\" >beast.math.distributions.Beta</map>");
  PhyML_Fprintf(fp,
                "\n<map name=\"Gamma\" >beast.math.distributions.Gamma</map>");
  PhyML_Fprintf(fp, "\n<map name=\"LaplaceDistribution\" "
                    ">beast.math.distributions.LaplaceDistribution</map>");
  PhyML_Fprintf(fp,
                "\n<map name=\"prior\" >beast.math.distributions.Prior</map>");
  PhyML_Fprintf(fp, "\n<map name=\"InverseGamma\" "
                    ">beast.math.distributions.InverseGamma</map>");
  PhyML_Fprintf(
      fp, "\n<map name=\"OneOnX\" >beast.math.distributions.OneOnX</map>");

  PhyML_Fprintf(fp,
                "\n<run id=\"mcmc\" spec=\"MCMC\" chainLength=\"1000000000\">");
  PhyML_Fprintf(fp, "\n<state id=\"state\" storeEvery=\"10000\">");
  PhyML_Fprintf(
      fp, "\n<stateNode id=\"Tree.t:data\" "
          "spec=\"beast.evolution.tree.StructuredCoalescentMultiTypeTree\">");
  PhyML_Fprintf(fp, "\n<migrationModel id=\"migModelInit.t:data\" "
                    "spec=\"beast.evolution.tree.MigrationModel\">");

  s = (char *)mCalloc(T_MAX_LINE, sizeof(char));
  For(i, n_demes * (n_demes - 1)) strcat(s, "1.0 ");
  PhyML_Fprintf(fp,
                "\n<parameter id=\"RealParameter.0\" dimension=\"%d\" "
                "estimate=\"false\" name=\"rateMatrix\">%s</parameter>",
                n_demes * (n_demes - 1), s);
  Free(s);

  s = (char *)mCalloc(T_MAX_LINE, sizeof(char));
  for (i = 0; i < n_demes; i++) strcat(s, "1.0 ");
  PhyML_Fprintf(fp,
                "\n<parameter id=\"RealParameter.01\" dimension=\"%d\" "
                "estimate=\"false\" name=\"popSizes\">%s</parameter>",
                n_demes, s);
  Free(s);

  PhyML_Fprintf(fp, "\n</migrationModel>");

  PhyML_Fprintf(
      fp, "\n<typeTrait id=\"typeTraitSet.t:data\" "
          "spec=\"beast.evolution.tree.TraitSet\" traitname=\"type\" value=\"");

  for (i = 0; i < tree->n_otu; i++)
  {
    s = strrchr(tree->a_nodes[i]->ldsk->coord->id, '_');
    PhyML_Fprintf(fp, "%s=%s", tree->a_nodes[i]->ldsk->coord->id, s + 1);

    if (i < tree->n_otu - 1)
      PhyML_Fprintf(fp, ",");
    else
      PhyML_Fprintf(fp, "\">");
  }

  PhyML_Fprintf(fp, "\n<taxa id=\"TaxonSet.0\" spec=\"TaxonSet\">");
  PhyML_Fprintf(fp, "\n<alignment idref=\"data\"/>");
  PhyML_Fprintf(fp, "\n</taxa>");
  PhyML_Fprintf(fp, "\n</typeTrait>");
  PhyML_Fprintf(fp, "\n<taxonset idref=\"TaxonSet.0\"/>");
  PhyML_Fprintf(fp, "\n</stateNode>");
  PhyML_Fprintf(fp, "\n<parameter id=\"kappa.s:data\" lower=\"0.0\" "
                    "name=\"stateNode\">2.0</parameter>");

  s = (char *)mCalloc(T_MAX_LINE, sizeof(char));
  for (i = 0; i < n_demes; i++) strcat(s, "1.0 ");
  PhyML_Fprintf(fp,
                "\n<parameter id=\"popSizes.t:data\" dimension=\"%d\" "
                "name=\"stateNode\">%s</parameter>",
                n_demes, s);
  Free(s);

  s = (char *)mCalloc(T_MAX_LINE, sizeof(char));
  For(i, n_demes * (n_demes - 1)) strcat(s, "1.0 ");
  PhyML_Fprintf(fp,
                "\n<parameter id=\"rateMatrix.t:data\" dimension=\"%d\" "
                "name=\"stateNode\">%s</parameter>",
                n_demes * (n_demes - 1), s);
  Free(s);

  PhyML_Fprintf(
      fp, "\n<parameter id=\"freqParameter.s:data\" dimension=\"4\" "
          "lower=\"0.0\" name=\"stateNode\" upper=\"1.0\">0.25</parameter>");
  PhyML_Fprintf(fp, "\n</state>");

  PhyML_Fprintf(
      fp,
      "\n<distribution id=\"posterior\" spec=\"util.CompoundDistribution\">");
  PhyML_Fprintf(
      fp, "\n<distribution id=\"prior\" spec=\"util.CompoundDistribution\">");
  PhyML_Fprintf(fp, "\n<prior id=\"KappaPrior.s:data\" name=\"distribution\" "
                    "x=\"@kappa.s:data\">");
  PhyML_Fprintf(
      fp, "\n<LogNormal id=\"LogNormalDistributionModel.0\" name=\"distr\">");
  PhyML_Fprintf(fp, "\n<parameter id=\"RealParameter.02\" estimate=\"false\" "
                    "name=\"M\">1.0</parameter>");
  PhyML_Fprintf(fp, "\n<parameter id=\"RealParameter.03\" estimate=\"false\" "
                    "name=\"S\">1.25</parameter>");
  PhyML_Fprintf(fp, "\n</LogNormal>");
  PhyML_Fprintf(fp, "\n</prior>");

  PhyML_Fprintf(fp, "\n<prior id=\"popSizesPrior.t:data\" "
                    "name=\"distribution\" x=\"@popSizes.t:data\">");
  PhyML_Fprintf(
      fp, "\n<LogNormal id=\"LogNormalDistributionModel.01\" name=\"distr\">");
  PhyML_Fprintf(fp, "\n<parameter id=\"RealParameter.04\" estimate=\"false\" "
                    "name=\"M\">1.0</parameter>");
  PhyML_Fprintf(fp, "\n<parameter id=\"RealParameter.05\" estimate=\"false\" "
                    "lower=\"0.0\" name=\"S\" upper=\"5.0\">1.25</parameter>");
  PhyML_Fprintf(fp, "\n</LogNormal>");
  PhyML_Fprintf(fp, "\n</prior>");

  PhyML_Fprintf(fp, "\n<prior id=\"rateMatrixPrior.t:data\" "
                    "name=\"distribution\" x=\"@rateMatrix.t:data\">");
  PhyML_Fprintf(
      fp, "\n<LogNormal id=\"LogNormalDistributionModel.02\" name=\"distr\">");
  PhyML_Fprintf(fp, "\n<parameter id=\"RealParameter.06\" estimate=\"false\" "
                    "name=\"M\">1.0</parameter>");
  PhyML_Fprintf(fp, "\n<parameter id=\"RealParameter.07\" estimate=\"false\" "
                    "lower=\"0.0\" name=\"S\" upper=\"5.0\">1.25</parameter>");
  PhyML_Fprintf(fp, "\n</LogNormal>");
  PhyML_Fprintf(fp, "\n</prior>");

  PhyML_Fprintf(
      fp, "\n<distribution id=\"structuredCoalescent.t:data\" "
          "spec=\"multitypetree.distributions."
          "StructuredCoalescentTreeDensity\" multiTypeTree=\"@Tree.t:data\">");
  PhyML_Fprintf(
      fp, "\n<migrationModel id=\"migModel.t:data\" "
          "spec=\"beast.evolution.tree.MigrationModel\" "
          "popSizes=\"@popSizes.t:data\" rateMatrix=\"@rateMatrix.t:data\">");
  PhyML_Fprintf(fp, "\n</migrationModel>");
  PhyML_Fprintf(fp, "\n</distribution>");

  PhyML_Fprintf(
      fp,
      "\n<distribution id=\"likelihood\" spec=\"util.CompoundDistribution\">");
  PhyML_Fprintf(
      fp, "\n<distribution id=\"treeLikelihood.data\" spec=\"TreeLikelihood\" "
          "data=\"@data\" tree=\"@Tree.t:data\">");
  PhyML_Fprintf(fp, "\n<siteModel id=\"SiteModel.s:data\" spec=\"SiteModel\">");
  PhyML_Fprintf(fp, "\n<parameter id=\"mutationRate.s:data\" "
                    "estimate=\"false\" name=\"mutationRate\">1.0</parameter>");
  PhyML_Fprintf(fp, "\n<parameter id=\"gammaShape.s:data\" estimate=\"false\" "
                    "name=\"shape\">1.0</parameter>");
  PhyML_Fprintf(fp,
                "\n<parameter id=\"proportionInvariant.s:data\" "
                "estimate=\"false\" lower=\"0.0\" name=\"proportionInvariant\" "
                "upper=\"1.0\">0.0</parameter>");
  PhyML_Fprintf(
      fp,
      "\n<substModel id=\"hky.s:data\" spec=\"HKY\" kappa=\"@kappa.s:data\">");
  PhyML_Fprintf(fp,
                "\n<frequencies id=\"estimatedFreqs.s:data\" "
                "spec=\"Frequencies\" frequencies=\"@freqParameter.s:data\"/>");
  PhyML_Fprintf(fp, "\n</substModel>");
  PhyML_Fprintf(fp, "\n</siteModel>");
  PhyML_Fprintf(fp,
                "\n<branchRateModel id=\"StrictClock.c:data\" "
                "spec=\"beast.evolution.branchratemodel.StrictClockModel\">");
  PhyML_Fprintf(fp,
                "\n<parameter id=\"clockRate.c:data\" estimate=\"false\" "
                "name=\"clock.rate\">%G</parameter>",
                1.0);
  PhyML_Fprintf(fp, "\n</branchRateModel>");
  PhyML_Fprintf(fp, "\n</distribution>");
  PhyML_Fprintf(fp, "\n</distribution>");
  PhyML_Fprintf(fp, "\n</distribution>");
  PhyML_Fprintf(fp, "\n</distribution>");
  PhyML_Fprintf(fp, "\n");
  PhyML_Fprintf(fp, "\n<operator id=\"STX.t:data\" "
                    "spec=\"multitypetree.operators.TypedSubtreeExchange\" "
                    "migrationModel=\"@migModel.t:data\" "
                    "multiTypeTree=\"@Tree.t:data\" weight=\"10.0\"/>");
  PhyML_Fprintf(fp, "\n<operator id=\"TWB.t:data\" "
                    "spec=\"multitypetree.operators.TypedWilsonBalding\" "
                    "alpha=\"0.2\" migrationModel=\"@migModel.t:data\" "
                    "multiTypeTree=\"@Tree.t:data\" weight=\"10.0\"/>");
  PhyML_Fprintf(fp, "\n<operator id=\"NR.t:data\" "
                    "spec=\"multitypetree.operators.NodeRetype\" "
                    "migrationModel=\"@migModel.t:data\" "
                    "multiTypeTree=\"@Tree.t:data\" weight=\"10.0\"/>");
  PhyML_Fprintf(
      fp, "\n<operator id=\"NSR1.t:data\" "
          "spec=\"multitypetree.operators.NodeShiftRetype\" "
          "migrationModel=\"@migModel.t:data\" multiTypeTree=\"@Tree.t:data\" "
          "rootOnly=\"true\" weight=\"10.0\"/>");
  PhyML_Fprintf(
      fp, "\n<operator id=\"NSR2.t:data\" "
          "spec=\"multitypetree.operators.NodeShiftRetype\" "
          "migrationModel=\"@migModel.t:data\" multiTypeTree=\"@Tree.t:data\" "
          "noRoot=\"true\" weight=\"10.0\"/>");
  PhyML_Fprintf(fp, "<operator id=\"MTU.t:data\" "
                    "spec=\"multitypetree.operators.MultiTypeUniform\" "
                    "includeRoot=\"true\" migrationModel=\"@migModel.t:data\" "
                    "multiTypeTree=\"@Tree.t:data\" weight=\"10.0\"/>\n");
  PhyML_Fprintf(
      fp, "<operator id=\"MTTS.t:data\" "
          "spec=\"multitypetree.operators.MultiTypeTreeScale\" "
          "migrationModel=\"@migModel.t:data\" multiTypeTree=\"@Tree.t:data\" "
          "scaleFactor=\"0.98\" useOldTreeScaler=\"true\" weight=\"10.0\"/>\n");
  PhyML_Fprintf(
      fp, "\n<operator id=\"MTTUpDown.t:data\" "
          "spec=\"multitypetree.operators.MultiTypeTreeScale\" "
          "migrationModel=\"@migModel.t:data\" multiTypeTree=\"@Tree.t:data\" "
          "scaleFactor=\"0.98\" useOldTreeScaler=\"true\" weight=\"10.0\">");
  PhyML_Fprintf(fp, "\n<parameter idref=\"popSizes.t:data\"/>");
  PhyML_Fprintf(fp, "\n</operator>");
  PhyML_Fprintf(
      fp, "\n<operator id=\"KappaScaler.s:data\" spec=\"ScaleOperator\" "
          "parameter=\"@kappa.s:data\" scaleFactor=\"0.5\" weight=\"0.1\"/>");
  PhyML_Fprintf(
      fp,
      "\n<operator id=\"popSizesScaler.t:data\" spec=\"ScaleOperator\" "
      "parameter=\"@popSizes.t:data\" scaleFactor=\"0.8\" weight=\"1.0\"/>");
  PhyML_Fprintf(
      fp,
      "\n<operator id=\"rateMatrixScaler.t:data\" spec=\"ScaleOperator\" "
      "parameter=\"@rateMatrix.t:data\" scaleFactor=\"0.8\" weight=\"1.0\"/>");
  PhyML_Fprintf(
      fp, "\n<operator id=\"FrequenciesExchanger.s:data\" "
          "spec=\"DeltaExchangeOperator\" delta=\"0.01\" weight=\"0.1\">");
  PhyML_Fprintf(fp, "\n<parameter idref=\"freqParameter.s:data\"/>");
  PhyML_Fprintf(fp, "\n</operator>");
  PhyML_Fprintf(fp, "\n");
  PhyML_Fprintf(fp, "\n<logger id=\"tracelog\" fileName=\"$(filebase).log\" "
                    "logEvery=\"10000\">");
  PhyML_Fprintf(fp, "\n<log idref=\"likelihood\"/>");
  PhyML_Fprintf(fp, "\n<log idref=\"prior\"/>");
  PhyML_Fprintf(fp, "\n<log idref=\"treeLikelihood.data\"/>");
  PhyML_Fprintf(
      fp,
      "\n<log id=\"treeHeight.t:data\" "
      "spec=\"beast.evolution.tree.TreeHeightLogger\" tree=\"@Tree.t:data\"/>");
  /* PhyML_Fprintf(fp,"\n<log id=\"treeLength.t:data\"
   * spec=\"multitypetree.util.TreeLengthLogger\" tree=\"@Tree.t:data\"/>"); */
  /* PhyML_Fprintf(fp,"\n<log id=\"changeCounts.t:data\"
   * spec=\"multitypetree.util.TypeChangeCounts\"
   * migrationModel=\"@migModel.t:data\" multiTypeTree=\"@Tree.t:data\"/>"); */
  /* PhyML_Fprintf(fp,"\n<log id=\"rootTypeLogger.t:data\"
   * spec=\"multitypetree.util.TreeRootTypeLogger\"
   * multiTypeTree=\"@Tree.t:data\"/>"); */
  PhyML_Fprintf(
      fp,
      "\n<log id=\"migModelLogger.t:data\" "
      "spec=\"multitypetree.util.MigrationModelLogger\" "
      "migrationModel=\"@migModel.t:data\" multiTypeTree=\"@Tree.t:data\"/>");
  PhyML_Fprintf(fp, "\n</logger>");
  PhyML_Fprintf(fp, "\n");
  PhyML_Fprintf(fp, "\n<logger id=\"screenlog\" logEvery=\"50000\">");
  PhyML_Fprintf(fp,
                "\n<log id=\"ESS.0\" spec=\"util.ESS\" arg=\"@posterior\"/>");
  PhyML_Fprintf(fp, "\n<log idref=\"likelihood\"/>");
  PhyML_Fprintf(fp, "\n</logger>");
  PhyML_Fprintf(fp, "\n");
  PhyML_Fprintf(fp, "\n</run>");
  PhyML_Fprintf(fp, "\n</beast>");

  for (i = 0; i < n_demes; i++) Free(deme_names[i]);
  Free(deme_names);

  fclose(fp);
}
#endif

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

#if (defined PHYREX || PHYTIME)
void PHYREX_Print_MCMC_Stats(t_tree *tree)
{
  FILE *fp_stats;

  fp_stats = tree->io->fp_out_stats;

  if (tree->mcmc->run == 0)
  {
    if (XML_Search_Node_Attribute_Value("add", "true", NO, tree->xml_root) ==
        NULL)
    {
      PhyML_Fprintf(fp_stats, "\n# before rand glnL: %f alnL: %f",
                    tree->mmod->c_lnL, tree->c_lnL);
      PhyML_Fprintf(fp_stats, "\n# ninter: %d",
                    PHYREX_Total_Number_Of_Intervals(tree));
      PhyML_Fprintf(fp_stats, "\n# ncoal: %d",
                    PHYREX_Total_Number_Of_Coal_Disks(tree));
      PhyML_Fprintf(fp_stats, "\n# nhits: %d",
                    PHYREX_Total_Number_Of_Hit_Disks(tree));
      PhyML_Fprintf(fp_stats, "\n# true lbda: %f", tree->mmod->lbda);
      PhyML_Fprintf(fp_stats, "\n# true mu: %f", tree->mmod->mu);
      PhyML_Fprintf(fp_stats, "\n# true rad: %f", SLFV_Update_Radius(tree));
      PhyML_Fprintf(fp_stats, "\n# true sigsq: %f", tree->mmod->sigsq[0]);
      PhyML_Fprintf(fp_stats, "\n# true neigh. size: %f",
                    SLFV_Neighborhood_Size(tree));
      PhyML_Fprintf(fp_stats, "\n# fst-based estimate of neighborhood size: %f",
                    SLFV_Neighborhood_Size_Regression(tree));
      PhyML_Fprintf(fp_stats, "\n# nucleotide diversity: %f",
                    Nucleotide_Diversity(tree->data));
      PhyML_Fprintf(fp_stats, "\n# length of a generation: %G time units",
                    SLFV_Generation_Length(tree));
      PhyML_Fprintf(fp_stats, "\n# clock rate: %G subst. per time unit",
                    tree->rates->clock_r);
      PhyML_Fprintf(fp_stats, "\n# after rand glnL: %f alnL: %f",
                    tree->mmod->c_lnL, tree->c_lnL);
      PhyML_Fprintf(fp_stats, "\n# ninter: %d",
                    PHYREX_Total_Number_Of_Intervals(tree));
      PhyML_Fprintf(fp_stats, "\n# ncoal: %d",
                    PHYREX_Total_Number_Of_Coal_Disks(tree));
      PhyML_Fprintf(fp_stats, "\n# nhits: %d",
                    PHYREX_Total_Number_Of_Hit_Disks(tree));
      PhyML_Fprintf(fp_stats, "\n# start lbda: %f", tree->mmod->lbda);
      PhyML_Fprintf(fp_stats, "\n# start mu: %f", tree->mmod->mu);
      PhyML_Fprintf(fp_stats, "\n# start rad: %f", tree->mmod->rad);
      PhyML_Fprintf(fp_stats, "\n# dist. in tree: ");
      for (int i = 0; i < tree->n_otu - 1; ++i)
        for (int j = i + 1; j < tree->n_otu; ++j)
          PhyML_Fprintf(fp_stats, "%G ",
                        PHYREX_Dist_Between_Two_Ldsk(tree->a_nodes[i]->ldsk,
                                                     tree->a_nodes[j]->ldsk,
                                                     tree));
      PhyML_Fprintf(fp_stats, "\n# dist. in space: ");
      for (int i = 0; i < tree->n_otu - 1; ++i)
        for (int j = i + 1; j < tree->n_otu; ++j)
          PhyML_Fprintf(fp_stats, "%G ",
                        Euclidean_Distance(tree->a_nodes[i]->ldsk->coord,
                                           tree->a_nodes[j]->ldsk->coord));

      PhyML_Fprintf(fp_stats, "\n");
      PhyML_Fprintf(fp_stats, "%s\t", "sample");
      PhyML_Fprintf(fp_stats, "%s\t", "lnPost");
      PhyML_Fprintf(fp_stats, "%s\t", "lnLAlgn");
      if (VELOC_Is_Integrated_Velocity(tree->mmod) == YES)
      {
        PhyML_Fprintf(fp_stats, "%s\t", "lnLPIV");
        PhyML_Fprintf(fp_stats, "%s\t", "lnLVeloc");
        PhyML_Fprintf(fp_stats, "%s\t", "lnLSpac");
      }
      else
        PhyML_Fprintf(fp_stats, "%s\t", "lnLSpac");

      PhyML_Fprintf(fp_stats, "%s\t", "lnLRate");
      PhyML_Fprintf(fp_stats, "%s\t", "lnLTime");
      PhyML_Fprintf(fp_stats, "%s\t", "lnPSpac");
      PhyML_Fprintf(fp_stats, "%s\t", "lnPRate");
      PhyML_Fprintf(fp_stats, "%s\t", "lnPTime");
      PhyML_Fprintf(fp_stats, "%s\t", "substRate");
      PhyML_Fprintf(fp_stats, "%s\t", "clockRate");
      for (int i = 0; i < tree->mmod->n_dim; ++i)
        PhyML_Fprintf(fp_stats, "%s%s\t", "sigSq",
                      (i == 0) ? ("Lon") : ((i == 1) ? ("Lat") : ("xx")));
      PhyML_Fprintf(fp_stats, "%s\t", "nEff");
      PhyML_Fprintf(fp_stats, "%s\t", "growth");

      PhyML_Fprintf(fp_stats, "%s\t", "speed");
      PhyML_Fprintf(fp_stats, "%s\t", "displacement");
      PhyML_Fprintf(fp_stats, "%s\t", "speedfromveloc");

      if (IWN_Is_Iwn(tree->mmod) == YES)
        PhyML_Fprintf(fp_stats, "%s\t", "stickiness");
      if (IOU_Is_Iou(tree->mmod) == YES)
        PhyML_Fprintf(fp_stats, "%s\t", "stickiness");
      if (IOU_Is_Iou(tree->mmod) == YES)
        for (int i = 0; i < tree->mmod->n_dim; ++i)
          PhyML_Fprintf(fp_stats, "%s%s\t", "tren",
                        (i == 0) ? ("Lon") : ((i == 1) ? ("Lat") : ("xx")));

      if (tree->contmod->obs_model == YES)
        PhyML_Fprintf(fp_stats, "%s\t", "obsVarLon");
      if (tree->contmod->obs_model == YES)
        PhyML_Fprintf(fp_stats, "%s\t", "obsVarLat");

      PhyML_Fprintf(fp_stats, "%s\t", "rootTime");
      PhyML_Fprintf(fp_stats, "%s\t", "rootLon");
      PhyML_Fprintf(fp_stats, "%s\t", "rootLat");
      PhyML_Fprintf(fp_stats, "%s\t", "nu");
      PhyML_Fprintf(fp_stats, "%s\t", "rrNormFact");
      PhyML_Fprintf(fp_stats, "%s\t", "rrRateMult");
      PhyML_Fprintf(fp_stats, "%s\t", "rrwNormFact");

      PhyML_Fprintf(fp_stats, "%s\t", "accNarrowExchange");
      PhyML_Fprintf(fp_stats, "%s\t", "accWideExchange");
      PhyML_Fprintf(fp_stats, "%s\t", "accMoveRootTime");
      PhyML_Fprintf(fp_stats, "%s\t", "accMoveScaleTime");
      PhyML_Fprintf(fp_stats, "%s\t", "accNodeTimes");
      PhyML_Fprintf(fp_stats, "%s\t", "accRatesShrink");
      PhyML_Fprintf(fp_stats, "%s\t", "accSigsq");
      PhyML_Fprintf(fp_stats, "%s\t", "accSigsqScale");
      PhyML_Fprintf(fp_stats, "%s\t", "accNodeVeloc");
      PhyML_Fprintf(fp_stats, "%s\t", "accScaleVeloc");
      PhyML_Fprintf(fp_stats, "%s\t", "tuneRatesShrink");
      PhyML_Fprintf(fp_stats, "%s\t", "tuneVeloc");

      if (tree->mcmc->out_verbose == 1)
        if (tree->mmod->model_id == RRW_GAMMA ||
            tree->mmod->model_id == RRW_LOGNORMAL)
          for (int i = 0; i < 2 * tree->n_otu - 1; ++i)
            PhyML_Fprintf(fp_stats, "sigSq%d\t", i);

      PhyML_Fprintf(fp_stats, "sigSqSonLeft\t");
      PhyML_Fprintf(fp_stats, "sigSqSonRght\t");

      PhyML_Fprintf(fp_stats, "rrSonLeft\t");
      PhyML_Fprintf(fp_stats, "rrSonRght\t");

      PhyML_Fprintf(fp_stats, "root_VelocLon\t");
      PhyML_Fprintf(fp_stats, "root_VelocLat\t");

      if (tree->mcmc->out_verbose == 1)
      {
        for (int i = 0; i < tree->n_otu; ++i)
        {
          PhyML_Fprintf(fp_stats, "%s_VelocLon\t", tree->a_nodes[i]->name);
          PhyML_Fprintf(fp_stats, "%s_VelocLat\t", tree->a_nodes[i]->name);
        }

        if (VELOC_Is_Integrated_Velocity(tree->mmod) == YES)
        {
          for (int i = tree->n_otu; i < 2 * tree->n_otu - 1; ++i)
          {
            PhyML_Fprintf(fp_stats, "%d_VelocLon%s\t", i,
                          (tree->a_nodes[i] == tree->n_root->v[1] ||
                           tree->a_nodes[i] == tree->n_root->v[2])
                              ? "*"
                              : "");
            PhyML_Fprintf(fp_stats, "%d_VelocLat%s\t", i,
                          (tree->a_nodes[i] == tree->n_root->v[1] ||
                           tree->a_nodes[i] == tree->n_root->v[2])
                              ? "*"
                              : "");
          }
        }
      }

      if (tree->mcmc->out_verbose == 1)
      {
        for (int i = tree->n_otu; i < 2 * tree->n_otu - 1; ++i)
        {
          PhyML_Fprintf(fp_stats, "%d_Lon%s\t", i,
                        (tree->a_nodes[i] == tree->n_root->v[1] ||
                         tree->a_nodes[i] == tree->n_root->v[2])
                            ? "*"
                            : "");
          PhyML_Fprintf(fp_stats, "%d_Lat%s\t", i,
                        (tree->a_nodes[i] == tree->n_root->v[1] ||
                         tree->a_nodes[i] == tree->n_root->v[2])
                            ? "*"
                            : "");
        }
      }
    }
  }

  if (!(tree->mcmc->run % tree->mcmc->sample_interval) &&
      tree->mcmc->sample_interval > 0)
  {

    if (tree->eval_glnL == YES)
      LOCATION_Lk(NULL,
                  tree); // Required in order to have contmod->lnL up-to-date
    if (tree->eval_tlnL == YES) TIMES_Lk(tree);
    if (tree->eval_alnL == YES) Lk(NULL, tree);
    if (tree->eval_rlnL == YES) RATES_Lk(tree);

    time(&(tree->mcmc->time_end));

    PhyML_Fprintf(fp_stats, "\n");
    PhyML_Fprintf(fp_stats, "%6d\t", tree->mcmc->run);
    PhyML_Fprintf(fp_stats, "%.2f\t", PHYREX_Get_Posterior(tree));
    PhyML_Fprintf(fp_stats, "%.2f\t", tree->c_lnL);
    if (VELOC_Is_Integrated_Velocity(tree->mmod) == YES)
    {
      PhyML_Fprintf(fp_stats, "%.2f\t", tree->mmod->c_lnL);
      PhyML_Fprintf(fp_stats, "%.2f\t",
                    tree->contmod->lnL[VELOCITY * tree->mmod->n_dim + 0] +
                        tree->contmod->lnL[VELOCITY * tree->mmod->n_dim + 1]);
      PhyML_Fprintf(fp_stats, "%.2f\t",
                    tree->contmod->lnL[LOCATION * tree->mmod->n_dim + 0] +
                        tree->contmod->lnL[LOCATION * tree->mmod->n_dim + 1]);
    }
    else
      PhyML_Fprintf(fp_stats, "%.2f\t", tree->mmod->c_lnL);

    PhyML_Fprintf(fp_stats, "%.2f\t", tree->rates->c_lnL);
    PhyML_Fprintf(fp_stats, "%.2f\t", tree->times->c_lnL);
    PhyML_Fprintf(fp_stats, "%f\t", tree->mmod->c_lnP);
    PhyML_Fprintf(fp_stats, "%.2f\t", tree->rates->c_lnP);
    PhyML_Fprintf(fp_stats, "%.2f\t", tree->times->c_lnP);
    PhyML_Fprintf(fp_stats, "%g\t", RATES_Realized_Substitution_Rate(tree));
    PhyML_Fprintf(fp_stats, "%g\t", tree->rates->clock_r);
    for (int i = 0; i < tree->mmod->n_dim; ++i)
      PhyML_Fprintf(fp_stats, "%g\t", tree->mmod->sigsq[i]);
    PhyML_Fprintf(fp_stats, "%g\t", tree->times->scaled_pop_size);
    PhyML_Fprintf(fp_stats, "%g\t", tree->times->neff_growth);

    PhyML_Fprintf(fp_stats, "%g\t",
                  PHYREX_Realized_Dispersal_Dist(tree->mmod->dist_type, tree) /
                      PHYREX_Time_Tree_Length(tree));
    PhyML_Fprintf(
        fp_stats, "%g\t",
        PHYREX_Realized_Displacement_Dist(tree->mmod->dist_type, tree) /
            PHYREX_Time_Tree_Length(tree));
    if (VELOC_Is_Integrated_Velocity(tree->mmod) == YES)
      PhyML_Fprintf(fp_stats, "%g\t", VELOC_Mean_Speed(tree));
    else
      PhyML_Fprintf(fp_stats, "-1.\t");

    if (IWN_Is_Iwn(tree->mmod) == YES)
      PhyML_Fprintf(fp_stats, "%g\t", tree->mmod->omega);
    if (IOU_Is_Iou(tree->mmod) == YES)
      PhyML_Fprintf(fp_stats, "%g\t", tree->mmod->ou_theta);
    if (IOU_Is_Iou(tree->mmod) == YES)
      for (int i = 0; i < tree->mmod->n_dim; ++i)
        PhyML_Fprintf(fp_stats, "%g\t", tree->mmod->ou_mu[i]);

    if (tree->contmod->obs_model == YES)
      PhyML_Fprintf(fp_stats, "%g\t", tree->contmod->obs_var[0]);
    if (tree->contmod->obs_model == YES)
      PhyML_Fprintf(fp_stats, "%g\t", tree->contmod->obs_var[1]);

    PhyML_Fprintf(fp_stats, "%.2f\t", tree->n_root->ldsk->disk->time);

    if (RRW_Is_Rw(tree->mmod) == YES &&
        tree->mmod->integrateAncestralLocations == YES)
      RRW_Sample_Node_Locations_Joint(tree);
    if (VELOC_Is_Integrated_Velocity(tree->mmod) == YES &&
        tree->mmod->integrateAncestralLocations == YES)
      VELOC_Sample_Node_Locations_Joint(tree);

    PhyML_Fprintf(fp_stats, "%g\t", tree->n_root->ldsk->coord->lonlat[0]);
    PhyML_Fprintf(fp_stats, "%g\t", tree->n_root->ldsk->coord->lonlat[1]);

    PhyML_Fprintf(fp_stats, "%g\t", tree->rates->nu);
    PhyML_Fprintf(fp_stats, "%g\t", tree->rates->norm_fact);
    PhyML_Fprintf(fp_stats, "%g\t", RATES_Mean_Rate_Multiplier(tree));
    PhyML_Fprintf(fp_stats, "%g\t", tree->mmod->sigsq_scale_norm_fact);
    PhyML_Fprintf(
        fp_stats, "%g\t",
        tree->mcmc->acc_rate[tree->mcmc->num_move_phyrex_narrow_exchange]);
    PhyML_Fprintf(
        fp_stats, "%g\t",
        tree->mcmc->acc_rate[tree->mcmc->num_move_phyrex_wide_exchange]);
    PhyML_Fprintf(fp_stats, "%g\t",
                  tree->mcmc->acc_rate[tree->mcmc->num_move_root_time]);
    PhyML_Fprintf(
        fp_stats, "%g\t",
        tree->mcmc->acc_rate[tree->mcmc->num_move_phyrex_scale_times]);
    PhyML_Fprintf(fp_stats, "%g\t",
                  tree->mcmc->acc_rate[tree->mcmc->num_move_phyrex_node_times]);
    PhyML_Fprintf(fp_stats, "%g\t",
                  tree->mcmc->acc_rate[tree->mcmc->num_move_rates_shrink]);
    PhyML_Fprintf(fp_stats, "%g\t",
                  tree->mcmc->acc_rate[tree->mcmc->num_move_phyrex_sigsq]);
    PhyML_Fprintf(
        fp_stats, "%g\t",
        tree->mcmc->acc_rate[tree->mcmc->num_move_phyrex_sigsq_scale]);
    PhyML_Fprintf(fp_stats, "%g\t",
                  tree->mcmc->acc_rate[tree->mcmc->num_move_phyrex_node_veloc]);
    PhyML_Fprintf(fp_stats, "%g\t",
                  tree->mcmc->acc_rate[tree->mcmc->num_move_phyrex_all_veloc]);
    PhyML_Fprintf(fp_stats, "%g\t",
                  tree->mcmc->tune_move[tree->mcmc->num_move_rates_shrink]);
    PhyML_Fprintf(
        fp_stats, "%g\t",
        tree->mcmc->tune_move[tree->mcmc->num_move_phyrex_node_veloc]);

    if (tree->mcmc->out_verbose == 1)
      if (tree->mmod->model_id == RRW_GAMMA ||
          tree->mmod->model_id == RRW_LOGNORMAL)
        for (int i = 0; i < 2 * tree->n_otu - 1; ++i)
          PhyML_Fprintf(fp_stats, "%g\t", tree->mmod->sigsq_scale[i]);

    PhyML_Fprintf(fp_stats, "%f\t",
                  tree->mmod->sigsq_scale[tree->n_root->v[1]->num]);
    PhyML_Fprintf(fp_stats, "%f\t",
                  tree->mmod->sigsq_scale[tree->n_root->v[2]->num]);

    PhyML_Fprintf(fp_stats, "%f\t", tree->rates->br_r[tree->n_root->v[1]->num]);
    PhyML_Fprintf(fp_stats, "%f\t", tree->rates->br_r[tree->n_root->v[2]->num]);

    if (VELOC_Is_Integrated_Velocity(tree->mmod) == YES)
    {
      PhyML_Fprintf(fp_stats, "%g\t", tree->n_root->ldsk->veloc->deriv[0]);
      PhyML_Fprintf(fp_stats, "%g\t", tree->n_root->ldsk->veloc->deriv[1]);

      if (tree->mcmc->out_verbose == 1)
        for (int i = 0; i < 2 * tree->n_otu - 1; ++i)
        {
          PhyML_Fprintf(fp_stats, "%g\t",
                        tree->a_nodes[i]->ldsk->veloc->deriv[0]);
          PhyML_Fprintf(fp_stats, "%g\t",
                        tree->a_nodes[i]->ldsk->veloc->deriv[1]);
        }
    }
    else if (RRW_Is_Rw(tree->mmod) == YES)
    {
      RRW_Tip_Velocities(tree);

      PhyML_Fprintf(fp_stats, "%g\t", -1.);
      PhyML_Fprintf(fp_stats, "%g\t", -1.);

      if (tree->mcmc->out_verbose == 1)
      {
        for (int i = 0; i < tree->n_otu; ++i)
        {
          PhyML_Fprintf(fp_stats, "%g\t",
                        tree->a_nodes[i]->ldsk->veloc->deriv[0]);
          PhyML_Fprintf(fp_stats, "%g\t",
                        tree->a_nodes[i]->ldsk->veloc->deriv[1]);
        }
      }
    }

    if (tree->mcmc->out_verbose == 1)
      for (int i = tree->n_otu; i < 2 * tree->n_otu - 1; ++i)
      {
        PhyML_Fprintf(fp_stats, "%g\t",
                      tree->a_nodes[i]->ldsk->coord->lonlat[0]);
        PhyML_Fprintf(fp_stats, "%g\t",
                      tree->a_nodes[i]->ldsk->coord->lonlat[1]);
      }

    /* for(int i=0;i<tree->n_otu;++i) */
    /*   { */
    /*     PhyML_Fprintf(fp_stats,"%g\t%g\t", */
    /*                   (tree->a_nodes[i]->ldsk->coord->lonlat[0] -
     * tree->a_nodes[i]->anc->ldsk->coord->lonlat[0])/ */
    /*                   fabs(tree->times->nd_t[tree->a_nodes[i]->num] -
     * tree->times->nd_t[tree->a_nodes[i]->anc->num]), */
    /*                   (tree->a_nodes[i]->ldsk->coord->lonlat[1] -
     * tree->a_nodes[i]->anc->ldsk->coord->lonlat[1])/ */
    /*                   fabs(tree->times->nd_t[tree->a_nodes[i]->num] -
     * tree->times->nd_t[tree->a_nodes[i]->anc->num])); */
    /*   } */

    /* for(int i=0;i<tree->n_otu;++i) */
    /*   { */
    /*     PhyML_Fprintf(fp_stats,"%g\t", */
    /*                   Euclidean_Distance(tree->a_nodes[i]->ldsk->coord,tree->a_nodes[i]->anc->ldsk->coord)/
     */
    /*                   fabs(tree->times->nd_t[tree->a_nodes[i]->num] -
     * tree->times->nd_t[tree->a_nodes[i]->anc->num])); */

    /*   } */

    /* for(int i=0;i<tree->n_otu-1;++i) */
    /*   { */
    /*     PhyML_Fprintf(fp_stats,"%g\t",tree->a_nodes[tree->n_otu+i]->ldsk->coord->lonlat[0]);
     */
    /*     PhyML_Fprintf(fp_stats,"%g\t",tree->a_nodes[tree->n_otu+i]->ldsk->coord->lonlat[1]);
     */
    /*   } */

    /* if(tree->io->mcmc_output_times == YES) */
    /*   { */
    /*     for(int i=tree->n_otu;i<2*tree->n_otu-1;++i)
     * PhyML_Fprintf(fp_stats,"%f\t",tree->times->nd_t[tree->a_nodes[i]->num]);
     */
    /*   } */

    fflush(NULL);

    time(&(tree->mcmc->time_beg));

    tree->mcmc->sample_num++;
  }
}
#endif

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

#if (defined PHYREX || defined TEST)
void PHYREX_Print_MCMC_Tree(t_tree *tree)
{
  FILE *fp_tree;

  if (tree->io->mcmc_output_trees == NO) return;

  fp_tree = tree->io->fp_out_tree;

  if (tree->mcmc->run == 0)
  {
    PhyML_Fprintf(fp_tree, "#NEXUS");
    PhyML_Fprintf(fp_tree, "\n\nBegin taxa;");
    PhyML_Fprintf(fp_tree, "\n\tDimensions ntax=%d;", tree->n_otu);
    PhyML_Fprintf(fp_tree, "\n\tTaxlabels");
    for (int i = 0; i < tree->n_otu; ++i)
      PhyML_Fprintf(fp_tree, "\n\t\t'%s'", tree->a_nodes[i]->name);
    PhyML_Fprintf(fp_tree, "\n\t;");
    PhyML_Fprintf(fp_tree, "\n\tend;");
    PhyML_Fprintf(fp_tree, "\n\n");
    PhyML_Fprintf(fp_tree, "\nBegin trees;");
    PhyML_Fprintf(fp_tree, "\n\tTranslate");
    for (int i = 0; i < tree->n_otu; ++i)
    {
      PhyML_Fprintf(fp_tree, "\n\t%d '%s'", i + 1, tree->a_nodes[i]->name);
      if (i < tree->n_otu - 1) PhyML_Fprintf(fp_tree, ",");
    }
    PhyML_Fprintf(fp_tree, "\n;");
  }

  if (!(tree->mcmc->run % tree->mcmc->sample_interval) &&
      tree->mcmc->sample_interval > 0)
  {
    /* if(RRW_Is_Rw(tree->mmod) == YES &&
     * tree->mmod->integrateAncestralLocations == YES)
     * RRW_Sample_Node_Locations_Joint(tree); */
    /* if(VELOC_Is_Integrated_Velocity(tree->mmod) == YES &&
     * tree->mmod->integrateAncestralLocations == YES)
     * VELOC_Sample_Node_Locations_Joint(tree); */
    Record_Br_Len(tree);
    TIMES_Time_To_Bl(tree);
    tree->bl_ndigits = 3;
    /* tree->bl_ndigits = 7; */
    tree->write_tax_names = NO;
    Label_Nodes_With_Locations(tree);
    if (VELOC_Is_Integrated_Velocity(tree->mmod) == YES)
      Label_Nodes_With_Velocities(tree);
    else if (RRW_Is_Rw(tree->mmod) == YES)
    {
      RRW_Tip_Velocities(tree);
      Label_Nodes_With_Velocities(tree);
    }
    Label_Edges(tree);
    char *s = Write_Tree(tree);
    Free_All_Node_Labels(tree);
    PhyML_Fprintf(fp_tree,
                  "\ntree %d [&lnP=%f,precision={1.e-1,1e-02,1e-01}] = [&R] %s",
                  tree->mcmc->run, tree->c_lnL, s);
    PhyML_Fprintf(fp_tree, "\nend;");
    fseek(fp_tree, -5, SEEK_END);
    tree->write_tax_names = YES;
    Free(s);
    Restore_Br_Len(tree);
    fflush(NULL);
  }
}
#endif

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

#if (defined PHYREX || defined TEST)
void PHYREX_Print_MCMC_Summary(t_tree *tree)
{
  char *s;

  s = (char *)mCalloc(100, sizeof(char));

  strcpy(s, "");

  if (tree->mcmc->run == 0)
  {
    if (VELOC_Is_Integrated_Velocity(tree->mmod) == YES ||
        RRW_Is_Rw(tree->mmod) == YES)
      sprintf(s, "%s", "speed\0");

    PhyML_Fprintf(
        stdout, "\n. %23s\t%25s\t%15s\t%15s\t%15s\t%13s\t%10s\t%19s\t%13s",
        "run (aux.run)", "operator", "lnSpac",
        /* "lnTime", */
        "lnAlgn", "lnPost", "rootTime", "substRate", "(rootLon;rootLat)", s);
  }

  if (!(tree->mcmc->run % tree->mcmc->print_every) &&
      tree->mcmc->print_every > 0)
  {
    if (RRW_Is_Rw(tree->mmod) == YES &&
        tree->mmod->integrateAncestralLocations == YES)
      RRW_Sample_Node_Locations_Joint(tree);
    if (VELOC_Is_Integrated_Velocity(tree->mmod) == YES &&
        tree->mmod->integrateAncestralLocations == YES)
      VELOC_Sample_Node_Locations_Joint(tree);

    if (VELOC_Is_Integrated_Velocity(tree->mmod) == YES)
      sprintf(s, "%13f", VELOC_Mean_Speed(tree));
    else if (RRW_Is_Rw(tree->mmod) == YES)
      sprintf(s, "%13f",
              PHYREX_Realized_Dispersal_Dist(tree->mmod->dist_type, tree) /
                  PHYREX_Time_Tree_Length(tree));

    PhyML_Fprintf(stdout,
                  "\n. %10d "
                  "(%10d)\t%25s\t%15.2f\t%15.2f\t%15.2f\t%13.1f\t%10f\t(%8.2f;%"
                  "8.2f)\t%13s",
                  tree->mcmc->run, tree->aux_tree[0]->mcmc->run,
                  tree->mcmc->move_idx > -1
                      ? tree->mcmc->move_name[tree->mcmc->move_idx]
                      : "",
                  tree->mmod->c_lnL,
                  /* tree->times->c_lnL, */
                  tree->c_lnL, PHYREX_Get_Posterior(tree),
                  tree->times->nd_t[tree->n_root->num], tree->rates->clock_r,
                  tree->n_root->ldsk->coord->lonlat[0],
                  tree->n_root->ldsk->coord->lonlat[1], s);

    if (tree->numerical_warning == YES)
      PhyML_Fprintf(stdout,
                    " -- WARNING: numerical precision issue detected...");
  }

  Free(s);
}
#endif

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

xml_node *Generate_PhyREX_XMLObj(char *run_id, char *out_file_name,
                                 char *seq_file_name, char *coord_file_name)
{
  xml_node *root, *nd, *ndnd, *ndndnd;

  root = XML_Make_Node("phyrex");
  XML_Init_Node(NULL, root, "phyrex");
  root->attr = XML_Make_Attribute(NULL, "run.id", run_id);

  XML_Add_Attribute(root, "output.file", out_file_name);
  XML_Add_Attribute(root, "mcmc.chain.len", "2E+8");
  XML_Add_Attribute(root, "mcmc.sample.every", "1E+4");
  XML_Add_Attribute(root, "mcmc.print.every", "1E+4");
  XML_Add_Attribute(root, "mcmc.burnin", "1E+6");
  XML_Add_Attribute(root, "mutmap", "no");
  XML_Add_Attribute(root, "ignore.sequences", "no");
  XML_Add_Attribute(root, "mcmc.output.trees", "yes");

  nd       = XML_Add_Node(root, "spatialmodel");
  nd->attr = XML_Make_Attribute(NULL, "name", "ibm");
  XML_Add_Attribute(nd, "rw.prior.distrib", "flat");
  XML_Add_Attribute(nd, "rw.prior.mean", "1.");
  XML_Add_Attribute(nd, "sampling", "detection");
  XML_Add_Attribute(nd, "integrateAncestralLocations", "true");
  XML_Add_Attribute(nd, "distance.type", "greatcircle");
  XML_Add_Attribute(nd, "observational.model", "yes");

  nd       = XML_Add_Node(root, "lineagerates");
  nd->attr = XML_Make_Attribute(NULL, "model", "lognormal");
  XML_Add_Attribute(nd, "autocor.prior.rate", "100");

  nd       = XML_Add_Node(root, "treegenerating");
  nd->attr = XML_Make_Attribute(NULL, "model", "coalescent");
  XML_Add_Attribute(nd, "neff.prior.distrib", "flat");
  XML_Add_Attribute(nd, "fix.node.ages", "no");

  nd       = XML_Add_Node(root, "clockrate");
  nd->attr = XML_Make_Attribute(NULL, "prior.mean", "0.00058");
  XML_Add_Attribute(nd, "prior.var", "1.E-2");

  nd         = XML_Add_Node(root, "topology");
  ndnd       = XML_Add_Node(nd, "instance");
  ndnd->attr = XML_Make_Attribute(NULL, "id", "T1");
  XML_Add_Attribute(ndnd, "init.tree", "BioNJ");

  nd         = XML_Add_Node(root, "ratematrices");
  nd->attr   = XML_Make_Attribute(NULL, "id", "RM1");
  ndnd       = XML_Add_Node(nd, "instance");
  ndnd->attr = XML_Make_Attribute(NULL, "id", "M1");
  XML_Add_Attribute(ndnd, "model", "HKY85");
  XML_Add_Attribute(ndnd, "optimise.tstv", "yes");

  nd         = XML_Add_Node(root, "siterates");
  nd->attr   = XML_Make_Attribute(NULL, "id", "SR1");
  ndnd       = XML_Add_Node(nd, "instance");
  ndnd->attr = XML_Make_Attribute(NULL, "id", "R1");
  XML_Add_Attribute(ndnd, "init.value", "1.0");
  ndnd       = XML_Add_Node(nd, "weights");
  ndnd->attr = XML_Make_Attribute(NULL, "id", "D1");
  XML_Add_Attribute(ndnd, "family", "freerates");
  ndndnd       = XML_Add_Node(ndnd, "instance");
  ndndnd->attr = XML_Make_Attribute(NULL, "appliesto", "R1");
  XML_Add_Attribute(ndndnd, "value", "0.25");

  nd         = XML_Add_Node(root, "equfreqs");
  nd->attr   = XML_Make_Attribute(NULL, "id", "EF1");
  ndnd       = XML_Add_Node(nd, "instance");
  ndnd->attr = XML_Make_Attribute(NULL, "id", "F1");
  XML_Add_Attribute(ndnd, "optimise.freqs", "no");

  nd         = XML_Add_Node(root, "branchlengths");
  nd->attr   = XML_Make_Attribute(NULL, "id", "BL1");
  ndnd       = XML_Add_Node(nd, "instance");
  ndnd->attr = XML_Make_Attribute(NULL, "id", "L1");
  XML_Add_Attribute(ndnd, "optimise.lens", "yes");

  nd       = XML_Add_Node(root, "partitionelem");
  nd->attr = XML_Make_Attribute(NULL, "id", "partition1");
  XML_Add_Attribute(nd, "file.name", seq_file_name);
  XML_Add_Attribute(nd, "data.type", "nt");
  XML_Add_Attribute(nd, "interleaved", "no");

  ndnd       = XML_Add_Node(nd, "mixtureelem");
  ndnd->attr = XML_Make_Attribute(NULL, "list", "T1");
  ndnd       = XML_Add_Node(nd, "mixtureelem");
  ndnd->attr = XML_Make_Attribute(NULL, "list", "M1");
  ndnd       = XML_Add_Node(nd, "mixtureelem");
  ndnd->attr = XML_Make_Attribute(NULL, "list", "F1");
  ndnd       = XML_Add_Node(nd, "mixtureelem");
  ndnd->attr = XML_Make_Attribute(NULL, "list", "R1");
  ndnd       = XML_Add_Node(nd, "mixtureelem");
  ndnd->attr = XML_Make_Attribute(NULL, "list", "L1");

  nd       = XML_Add_Node(root, "coordinates");
  nd->attr = XML_Make_Attribute(NULL, "id", "coordinates");
  XML_Add_Attribute(nd, "file.name", coord_file_name);

  return (root);
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void Print_Trace(void)
{
  void  *array[10];
  char **strings;
  int    size, i;

  size    = backtrace(array, 10);
  strings = backtrace_symbols(array, size);

  if (strings != NULL)
  {
    PhyML_Printf("\n\n. Obtained %d stack frames.\n", size);
    for (i = 0; i < size; ++i) PhyML_Printf("\n. %s", strings[i]);

    Free(strings);
  }
}

/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/
