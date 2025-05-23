#include "xml.h"

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

t_tree *XML_Process_Base(char *xml_filename)
{
  FILE        *fp;
  xml_node    *root, *p_elem, *m_elem, *parent, *instance;
  option      *io;
  void        *buff;
  t_mod       *mod, *iomod;
  t_tree      *tree, *mixt_tree, *root_tree;
  int          select;
  char        *component;
  int          i, j, n_components;
  int          first_m_elem;
  int          class_number;
  scalar_dbl **lens, **lens_var, **ori_lens, **lens_old, **lens_var_old,
      **ori_lens_old, **ori_lens_var, **ori_lens_var_old;
  t_ds *ds;
  char *outputfile;
  char *alignment;
  char *s;
  int   lens_size;
  int   first;
  int  *class_num;
  int   n_partitionelem;

  fp = fopen(xml_filename, "r");
  if (!fp)
  {
    PhyML_Fprintf(stderr, "\n. Could not find the XML file '%s'.\n",
                  xml_filename);
    Exit("\n");
  }

  root = XML_Load_File(fp);

  if (!root)
  {
    PhyML_Fprintf(stderr,
                  "\n. Encountered an issue while loading the XML file.\n");
    Exit("\n");
  }

  class_num = (int *)mCalloc(N_MAX_MIXT_CLASSES, sizeof(int));
  for (i = 0; i < N_MAX_MIXT_CLASSES; i++) class_num[i] = i;

  component = (char *)mCalloc(T_MAX_NAME, sizeof(char));

  m_elem           = NULL;
  p_elem           = root;
  io               = NULL;
  mixt_tree        = NULL;
  root_tree        = NULL;
  mod              = NULL;
  tree             = NULL;
  lens             = NULL;
  lens_var         = NULL;
  lens_old         = NULL;
  lens_var_old     = NULL;
  select           = -1;
  class_number     = -1;
  ds               = NULL;
  lens_size        = 0;
  ori_lens         = NULL;
  ori_lens_old     = NULL;
  ori_lens_var     = NULL;
  ori_lens_var_old = NULL;
  first            = YES;

  // Make sure there are no duplicates in node's IDs
  XML_Check_Duplicate_ID(root);

  int count = 0;
  XML_Count_Number_Of_Node_With_Name("topology", &count, root);

  if (count > 1)
  {
    PhyML_Fprintf(stderr,
                  "\n. There should not more than one 'topology' node.");
    PhyML_Fprintf(stderr, "\n. Found %d. Please fix your XML file", count);
    Exit("\n");
  }
  else if (count < 1)
  {
    PhyML_Fprintf(stderr, "\n. There should be at least one 'topology' node.");
    PhyML_Fprintf(stderr, "\n. Found none. Please fix your XML file");
    Exit("\n");
  }

#if defined PHYML
  p_elem = XML_Search_Node_Name("phyml", NO, p_elem);
#elif defined PHYTIME
  p_elem = XML_Search_Node_Name("phytime", NO, p_elem);
#endif

  if (!p_elem)
  {
#if defined PHYML
    PhyML_Fprintf(stderr, "\n. The 'phyml' tag is mandatory.");
#elif defined PHYTIME
    PhyML_Fprintf(stderr, "\n. The 'phytime' tag is mandatory.");
#endif
    PhyML_Fprintf(stderr, "\n. Please amend your XML file accordingly.");
    assert(false);
  }

  /*! Input file
   */
  outputfile = XML_Get_Attribute_Value(p_elem, "output.file");

  if (!outputfile || !(strlen(outputfile) > 0))
  {
#if defined PHYML
    PhyML_Fprintf(
        stderr, "\n. The 'outputfile' attribute in 'phyml' tag is mandatory.");
#elif defined PHYTIME
    PhyML_Fprintf(
        stderr,
        "\n. The 'outputfile' attribute in 'phytime' tag is mandatory.");
#endif
    PhyML_Fprintf(stderr, "\n. Please amend your XML file accordingly.");
    assert(false);
  }

  io = (option *)Make_Input();
  Set_Defaults_Input(io);

  s = XML_Get_Attribute_Value(p_elem, "run.id");
  if (s && strlen(s) > 0)
  {
    io->append_run_ID = YES;
    strcpy(io->run_id_string, s);
  }

  s = XML_Get_Attribute_Value(p_elem, "r.seed");
  if (s && strlen(s) > 0) io->r_seed = (int)atoi(s);

  strcpy(io->out_file, outputfile);
  strcpy(io->out_tree_file, outputfile);
  strcpy(io->out_stats_file, outputfile);

#if defined(PHYTIME)
  strcat(io->out_tree_file, "_phytime_tree");
#elif defined(PHYREX)
  strcat(io->out_tree_file, "_phyrex_tree");
#elif defined(PHYML)
  strcat(io->out_tree_file, "_phyml_tree");
#endif

#if defined(PHYTIME)
  strcat(io->out_stats_file, "_phytime_stats");
#elif defined(PHYREX)
  strcat(io->out_stats_file, "_phyrex_stats");
#elif defined(PHYML)
  strcat(io->out_stats_file, "_phyml_stats");
#endif

  if (io->append_run_ID)
  {
    strcat(io->out_tree_file, "_");
    strcat(io->out_tree_file, io->run_id_string);
  }
  if (io->append_run_ID)
  {
    strcat(io->out_stats_file, "_");
    strcat(io->out_stats_file, io->run_id_string);
  }

  strcat(io->out_stats_file, ".txt");
  strcat(io->out_tree_file, ".txt");

  if (XML_Search_Node_Attribute_Value("add", "true", NO, root) != NULL)
  {
    io->fp_out_tree  = Openfile(io->out_tree_file, APPEND);
    io->fp_out_stats = Openfile(io->out_stats_file, APPEND);
  }
  else
  {
    io->fp_out_tree  = Openfile(io->out_tree_file, WRITE);
    io->fp_out_stats = Openfile(io->out_stats_file, WRITE);
  }

  s = XML_Get_Attribute_Value(p_elem, "print.trace");

  if (s)
  {
    select =
        XML_Validate_Attr_Int(s, 6, "true", "yes", "y", "false", "no", "n");

    if (select < 3)
    {
      io->print_trace = YES;
      strcpy(io->out_trace_file, outputfile);
      strcat(io->out_trace_file, "_phyml_trace");
      if (io->append_run_ID)
      {
        strcat(io->out_trace_file, "_");
        strcat(io->out_trace_file, io->run_id_string);
      }
      strcat(io->out_trace_file, ".txt");
      io->fp_out_trace = Openfile(io->out_trace_file, 1);
    }
  }

  s = XML_Get_Attribute_Value(p_elem, "print.json.trace");

  if (s)
  {
    select =
        XML_Validate_Attr_Int(s, 6, "true", "yes", "y", "false", "no", "n");

    if (select < 3)
    {
      io->print_json_trace = YES;
      strcpy(io->out_json_trace_file, outputfile);
      strcat(io->out_json_trace_file, "_phyml_trace");
      if (io->append_run_ID)
      {
        strcat(io->out_json_trace_file, "_");
        strcat(io->out_json_trace_file, io->run_id_string);
      }
      strcat(io->out_json_trace_file, ".json");
      io->fp_out_json_trace = Openfile(io->out_json_trace_file, READWRITE);
    }
  }

  s = XML_Get_Attribute_Value(p_elem, "print.site.lk");

  if (s)
  {
    select =
        XML_Validate_Attr_Int(s, 6, "true", "yes", "y", "false", "no", "n");
    if (select < 3) io->print_site_lnl = YES;

    strcpy(io->out_lk_file, outputfile);
    strcat(io->out_lk_file, "_phyml_lk");
    if (io->append_run_ID)
    {
      strcat(io->out_lk_file, "_");
      strcat(io->out_lk_file, io->run_id_string);
    }
    strcat(io->out_lk_file, ".txt");
    io->fp_out_lk = Openfile(io->out_lk_file, WRITE);
  }

  s = XML_Get_Attribute_Value(p_elem, "branch.test");
  if (s)
  {
    if (!strcmp(s, "aLRT"))
    {
      io->ratio_test = ALRTSTAT;
    }
    else if (!strcmp(s, "aBayes"))
    {
      io->ratio_test = ABAYES;
    }
    else if (!strcmp(s, "SH"))
    {
      io->ratio_test = SH;
    }
    else if (!strcmp(s, "no") || !strcmp(s, "none"))
    {
      io->ratio_test = 0;
    }
    else
    {
      PhyML_Fprintf(stderr, "\n. '%s' is not a valid option for 'branch.test'.",
                    s);
      PhyML_Fprintf(stderr, "\n. Please use 'aLRT' or 'aBayes' or 'SH'.");
      Exit("\n");
    }
  }

  s = XML_Get_Attribute_Value(p_elem, "quiet");
  if (s)
  {
    select =
        XML_Validate_Attr_Int(s, 6, "true", "yes", "y", "false", "no", "n");
    if (select < 3) io->quiet = YES;
  }

  s = XML_Get_Attribute_Value(p_elem, "memory.check");
  if (s)
  {
    select =
        XML_Validate_Attr_Int(s, 6, "true", "yes", "y", "false", "no", "n");
    if (select >= 3) io->mem_question = NO;
  }

  n_partitionelem = 0;
  /*! Read all partitionelem nodes and mixturelem nodes in each of them
   */
  do
  {
    p_elem = XML_Search_Node_Name("partitionelem", YES, p_elem);

    if (p_elem == NULL)
      break;
    else
      n_partitionelem++;

    buff = (option *)Make_Input();
    Set_Defaults_Input(buff);
    io->next       = buff;
    io->next->prev = io;

    io = io->next;
    if (first == YES)
    {
      io       = io->prev;
      io->next = NULL;
      Free_Input(buff);
      first = NO;
    }

    /*! Set the datatype (required when compressing data)
     */
    char *dt = NULL;
    dt       = XML_Get_Attribute_Value(p_elem, "data.type");
    if (!dt)
    {
      PhyML_Fprintf(stderr,
                    "\n. Please specify the type of data ('aa' or 'nt') for "
                    "partition element '%s'",
                    XML_Get_Attribute_Value(p_elem, "id"));
      PhyML_Fprintf(stderr,
                    "\n. Syntax: 'data.type=\"aa\"' or 'data.type=\"nt\"'");
      Exit("\n");
    }

    select = XML_Validate_Attr_Int(dt, 2, "aa", "nt");
    switch (select)
    {
    case 0:
    {
      io->datatype = AA;
      break;
    }
    case 1:
    {
      io->datatype = NT;
      break;
    }
    default:
    {
      PhyML_Fprintf(stderr,
                    "\n. Unknown data type. Must be either 'aa' or 'nt'.");
      Exit("\n");
    }
    }

    char *format = NULL;
    format       = XML_Get_Attribute_Value(p_elem, "format");
    if (format)
    {
      if (!strcmp(format, "interleave") || !strcmp(format, "interleaved"))
      {
        io->interleaved = YES;
      }
      else if (!strcmp(format, "sequential"))
      {
        io->interleaved = NO;
      }
    }

    /*! Attach a model to this io struct
     */
    io->mod = (t_mod *)Make_Model_Basic();
    Set_Defaults_Model(io->mod);
    io->mod->ras->n_catg = 1;
    io->mod->io          = io;
    iomod                = io->mod;

    if (io->datatype == AA)
      io->mod->ns = 20;
    else if (io->datatype == NT)
      io->mod->ns = 4;
    else
      assert(FALSE);

    /*! Attach an optimization structure to this model
     */
    iomod->s_opt = (t_opt *)Make_Optimiz();
    Set_Defaults_Optimiz(iomod->s_opt);

    /*! Input file
     */
    alignment = XML_Get_Attribute_Value(p_elem, "file.name");

    if (!alignment)
    {
      PhyML_Fprintf(stderr,
                    "\n. 'file.name' tag is mandatory. Please amend your");
      PhyML_Fprintf(stderr, "\n. XML file accordingly.");
      Exit("\n");
    }

    strcpy(io->in_align_file, alignment);

    /*! Open pointer to alignment
     */
    io->fp_in_align = Openfile(io->in_align_file, READ);

    /*! Load sequence file
     */
    io->data = Get_Seq(io);

    /*! Close pointer to alignment
     */
    fclose(io->fp_in_align);

    /*! Compress alignment
     */
    io->cdata = Compact_Data(io->data, io);

    /*! Free uncompressed alignment
     */
    Free_Seq(io->data, io->n_otu);

    /*! Create new mixture tree
     */
    buff = (t_tree *)Make_Tree_From_Scratch(io->cdata->n_otu, io->cdata);

    if (mixt_tree)
    {
      mixt_tree->next_mixt            = buff;
      mixt_tree->next_mixt->prev_mixt = mixt_tree;
      mixt_tree                       = mixt_tree->next_mixt;
      mixt_tree->dp                   = mixt_tree->prev_mixt->dp + 1;
    }
    else
      mixt_tree = buff;

    /*! Connect mixt_tree to xml struct
     */
    mixt_tree->xml_root = root;

    /*! Connect mixt_tree to io struct
     */
    mixt_tree->io = io;

    /*! Connect mixt_tree to model struct
     */
    mixt_tree->mod = iomod;

    /*! mixt_tree is a mixture tree
     */
    mixt_tree->is_mixt_tree = YES;
    mixt_tree->ignore_mixt_info = NO;

    /*! mixt_tree is a mixture tree
     */
    mixt_tree->mod->is_mixt_mod = YES;

    /*! Connect mixt_tree to compressed data
     */
    mixt_tree->data = io->cdata;

    /*! Remove branch lengths from mixt_tree */
    for (i = 0; i < 2 * mixt_tree->n_otu - 1; ++i)
    {
      Free_Scalar_Dbl(mixt_tree->a_edges[i]->l);
      Free_Scalar_Dbl(mixt_tree->a_edges[i]->l_old);
      Free_Scalar_Dbl(mixt_tree->a_edges[i]->l_var);
      Free_Scalar_Dbl(mixt_tree->a_edges[i]->l_var_old);
    }

    /*! Connect last tree of the mixture for the
      previous partition element to the next mixture tree
    */
    if (tree)
    {
      tree->next      = mixt_tree;
      mixt_tree->prev = tree;
    }

    if (!root_tree) root_tree = mixt_tree;

    /*! Tree size scaling factor */
    char *scale_tree = NULL;
    scale_tree       = XML_Get_Attribute_Value(p_elem, "optimise.tree.scale");

    if (scale_tree)
    {
      int select;

      select = XML_Validate_Attr_Int(scale_tree, 6, "true", "yes", "y", "false",
                                     "no", "n");

      if (select < 3) mixt_tree->mod->s_opt->opt_br_len_mult = YES;
    }

    scale_tree = NULL;
    scale_tree = XML_Get_Attribute_Value(p_elem, "tree.scale");

    if (scale_tree)
    {
      mixt_tree->mod->br_len_mult->v          = String_To_Dbl(scale_tree);
      mixt_tree->mod->br_len_mult_unscaled->v = String_To_Dbl(scale_tree);
      Free(scale_tree);
    }

    char *cv = NULL;
    cv       = XML_Get_Attribute_Value(root, "cv.type");

    if (cv)
    {
      if (!strcmp(cv, "kfold.pos"))
        mixt_tree->mod->cv_type = KFOLD_POS;
      else if (!strcmp(cv, "kfold.col"))
        mixt_tree->mod->cv_type = KFOLD_COL;
      else if (!strcmp(cv, "maxfold"))
        mixt_tree->mod->cv_type = MAXFOLD;
        else mixt_tree->mod->cv_type = -1;
    }
    else
      mixt_tree->mod->cv_type = -1;

    /*! Process all the mixtureelem tags in this partition element
     */
    n_components = 0;
    m_elem       = p_elem;
    first_m_elem = 0;
    mod          = NULL;
    tree         = NULL;
    class_number = 0;
    do
    {
      m_elem = XML_Search_Node_Name("mixtureelem", YES, m_elem);
      if (m_elem == NULL) break;

      if (!strcmp(m_elem->name, "mixtureelem"))
      {
        first_m_elem++;

        /*! Rewind tree and model when processing a new mixtureelem node
         */
        if (first_m_elem > 1)
        {
          while (tree->prev && tree->prev->is_mixt_tree == NO)
          {
            tree = tree->prev;
          }
          while (mod->prev && mod->prev->is_mixt_mod == NO)
          {
            mod = mod->prev;
          }
        }

        /*! Read and process model components
         */
        char *list;
        list = XML_Get_Attribute_Value(m_elem, "list");

        j = 0;
        for (i = 0; i < (int)strlen(list); ++i)
          if (list[i] == ',') j++;

        if (j != n_components && first_m_elem > 1)
        {
          PhyML_Fprintf(stderr,
                        "\n. Discrepancy in the number of elements in nodes "
                        "'mixtureelem' partitionelem id '%s'",
                        p_elem->id);
          PhyML_Fprintf(stderr, "\n. Check 'mixturelem' node with list '%s'",
                        list);
          Exit("\n");
        }
        n_components = j;

        i = j        = 0;
        component[0] = '\0';
        while (1)
        {
          if (list[j] == ',' || j == (int)strlen(list))
          {
            /*! Reading a new component
             */

            if (first_m_elem ==
                YES) /* Only true when processing the first mixtureelem node */
            {
              t_tree *this_tree;
              t_mod  *this_mod;

              /*! Create new tree
               */
              this_tree =
                  (t_tree *)Make_Tree_From_Scratch(io->cdata->n_otu, io->cdata);
              
              /*! Update the number of mixtures */
              iomod->n_mixt_classes++;

              if (tree)
              {
                tree->next       = this_tree;
                tree->next->prev = tree;
              }
              else
              {
                mixt_tree->next       = this_tree;
                mixt_tree->next->prev = mixt_tree;
              }

              tree            = this_tree;
              tree->mixt_tree = mixt_tree;

              /*! Create a new model
               */
              this_mod = (t_mod *)Make_Model_Basic();
              Set_Defaults_Model(this_mod);
              this_mod->ras->n_catg = 1;
              this_mod->ns          = iomod->ns;
              this_mod->r_mat       = (t_rmat *)Make_Rmat(this_mod->ns);
              Init_Rmat(this_mod->r_mat);

              /*! All br_len_multiplier point to the corresponding */
              /*! parameter in the relevant mixt_tree */
              Free_Scalar_Dbl(this_mod->br_len_mult);
              this_mod->br_len_mult = iomod->br_len_mult;

              Free_Scalar_Dbl(this_mod->br_len_mult_unscaled);
              this_mod->br_len_mult_unscaled = iomod->br_len_mult_unscaled;

              if (mod)
              {
                mod->next       = this_mod;
                mod->next->prev = mod;
              }
              else
              {
                this_mod->prev = iomod;
              }

              mod = this_mod;
              if (!iomod->next) iomod->next = mod;
              mod->io = io;

              mod->s_opt = (t_opt *)Make_Optimiz();
              Set_Defaults_Optimiz(mod->s_opt);

              mod->ras->alpha->optimize  = NO;
              mod->ras->pinvar->optimize = NO;

              tree->data = io->cdata;
              tree->io   = io;
              tree->mod  = mod;

              if (tree->data->n_pattern != tree->prev->data->n_pattern)
                Generic_Exit(__FILE__, __LINE__, __FUNCTION__);
            }

            /*! Read a component */
            component[i] = '\0';
            if (j != (int)strlen(list) - 1) i = 0;

            /*! Find which node this ID corresponds to
             */
            instance = XML_Search_Node_ID(component, YES, root);

            if (!instance)
            {
              PhyML_Fprintf(stderr, "\n. Could not find a node with id: '%s'.",
                            component);
              PhyML_Fprintf(stderr,
                            "\n. Problem with 'mixtureelem' node, list '%s'.",
                            list);
              Exit("\n");
            }

            if (!instance->parent)
            {
              PhyML_Fprintf(stderr, "\n. Node '%s' with id:'%s' has no parent.",
                            instance->name, component);
              Exit("\n");
            }

            parent = instance->parent;

            ////////////////////////////////////////
            //        SUBSTITUTION MODEL          //
            ////////////////////////////////////////

            if (!strcmp(parent->name, "ratematrices"))
            {
              /* ! First time we process this 'instance' node which has this
               * 'ratematrices' parent */
              if (instance->ds->obj == NULL)
              {
                Make_Ratematrix_From_XML_Node(instance, io, mod);

                ds = instance->ds;

                /*! Connect the data structure n->ds to mod->r_mat */
                ds->obj = (t_rmat *)(mod->r_mat);

                /*! Create and connect the data structure n->ds->next to
                 * mod->kappa */
                ds->next = (t_ds *)mCalloc(1, sizeof(t_ds));
                ds       = ds->next;
                ds->obj  = (scalar_dbl *)(mod->kappa);

                /*! Create and connect the data structure n->ds->next to
                 * mod->whichmodel */
                ds->next = (t_ds *)mCalloc(1, sizeof(t_ds));
                ds       = ds->next;
                ds->obj  = (int *)(&mod->whichmodel);

                /*! Create and connect the data structure n->ds->next to
                 * mod->modelname */
                ds->next = (t_ds *)mCalloc(1, sizeof(t_ds));
                ds       = ds->next;
                ds->obj  = (t_string *)(mod->modelname);

                /*! Create and connect the data structure n->ds->next to mod->ns
                 */
                ds->next = (t_ds *)mCalloc(1, sizeof(t_ds));
                ds       = ds->next;
                ds->obj  = (int *)(&mod->ns);

                /*! Create and connect the data structure n->ds->next to
                 * mod->modelname */
                ds->next = (t_ds *)mCalloc(1, sizeof(t_ds));
                ds       = ds->next;
                ds->obj  = (t_string *)(mod->custom_mod_string);

                /*! Create and connect the data structure n->ds->next to
                 * mod->fp_aa_rate_mat */
                ds->next = (t_ds *)mCalloc(1, sizeof(t_ds));
                ds       = ds->next;
                ds->obj  = (FILE *)(mod->fp_aa_rate_mat);

                /*! Create and connect the data structure n->ds->next to
                 * mod->aa_rate_mat_file */
                ds->next = (t_ds *)mCalloc(1, sizeof(t_ds));
                ds       = ds->next;
                ds->obj  = (t_string *)mod->aa_rate_mat_file;

                /*! Create and connect the data structure n->ds->next to
                 * mod->s_opt->opt_rmat_weight */
                ds->next = (t_ds *)mCalloc(1, sizeof(t_ds));
                ds       = ds->next;
                ds->obj  = (int *)(&mod->s_opt->opt_rmat_weight);

                /*! Create and connect the data structure n->ds->next to
                 * mod->gamma_mgf_bl */
                ds->next = (t_ds *)mCalloc(1, sizeof(t_ds));
                ds       = ds->next;
                ds->obj  = (int *)(&mod->gamma_mgf_bl);

                /*! Create and connect the data structure n->ds->next to
                 * mod->l_var_sigma */
                ds->next = (t_ds *)mCalloc(1, sizeof(t_ds));
                ds       = ds->next;
                ds->obj  = (scalar_dbl *)(mod->l_var_sigma);
              }
              else
              {
                /*! Connect to already existing r_mat & kappa
                 * structs. */
                t_ds *ds;

                Free(mod->r_mat);
                ds         = instance->ds;
                mod->r_mat = (t_rmat *)ds->obj;

                Free_Scalar_Dbl(mod->kappa);
                ds         = ds->next;
                mod->kappa = (scalar_dbl *)ds->obj;

                ds              = ds->next;
                mod->whichmodel = *((int *)ds->obj);

                Free_String(mod->modelname);
                ds             = ds->next;
                mod->modelname = (t_string *)ds->obj;

                ds      = ds->next;
                mod->ns = *((int *)ds->obj);

                Free_String(mod->custom_mod_string);
                ds                     = ds->next;
                mod->custom_mod_string = (t_string *)ds->obj;

                ds                  = ds->next;
                mod->fp_aa_rate_mat = (FILE *)ds->obj;

                Free_String(mod->aa_rate_mat_file);
                ds                    = ds->next;
                mod->aa_rate_mat_file = (t_string *)ds->obj;

                ds                          = ds->next;
                mod->s_opt->opt_rmat_weight = *((int *)ds->obj);

                ds                = ds->next;
                mod->gamma_mgf_bl = *((int *)ds->obj);

                Free_Scalar_Dbl(mod->l_var_sigma);
                ds               = ds->next;
                mod->l_var_sigma = (scalar_dbl *)ds->obj;
              }
            }

            ////////////////////////////////////////
            //           STATE FREQS              //
            ////////////////////////////////////////

            else if (!strcmp(parent->name, "equfreqs"))
            {
              /* If n->ds == NULL, the corrresponding node data structure,
               * n->ds, has not */
              /* been initialized. If not, do nothing. */
              if (instance->ds->obj == NULL)
              {
                Make_Efrq_From_XML_Node(instance, io, mod);

                t_ds *ds;

                ds      = instance->ds;
                ds->obj = (t_efrq *)mod->e_frq;

                ds->next = (t_ds *)mCalloc(1, sizeof(t_ds));
                ds       = ds->next;
                ds->obj  = (int *)(&mod->e_frq->type);

                ds->next = (t_ds *)mCalloc(1, sizeof(t_ds));
                ds       = ds->next;
                ds->obj  = (vect_dbl *)(mod->e_frq->user_b_freq);
              }
              else
              {
                /* Connect the data structure n->ds to mod->e_frq */

                ds         = instance->ds;
                mod->e_frq = (t_efrq *)ds->obj;

                ds               = ds->next;
                mod->e_frq->type = *((int *)ds->obj);

                ds                      = ds->next;
                mod->e_frq->user_b_freq = (vect_dbl *)ds->obj;
              }
            }

            //////////////////////////////////////////
            //             TOPOLOGY                 //
            //////////////////////////////////////////

            else if (!strcmp(parent->name, "topology"))
            {
              if (parent->ds->obj == NULL)
                Make_Topology_From_XML_Node(instance, io, iomod);

              ds = parent->ds;

              int buff;
              ds->obj = (int *)(&buff);
            }

            //////////////////////////////////////////
            //                RAS                   //
            //////////////////////////////////////////

            else if (!strcmp(parent->name, "siterates"))
            {
              char *rate_value = NULL;
              /* scalar_dbl *r; */
              phydbl val;

              /*! First time we process this 'siterates' node, check that its
                format is valid. and process it afterwards.
              */
              if (parent->ds->obj == NULL)
              {
                class_number = 0;

                Make_RAS_From_XML_Node(parent, iomod);

                ds = parent->ds;

                ds->obj = (t_ras *)iomod->ras;

                ds->next = (t_ds *)mCalloc(1, sizeof(t_ds));
                ds       = ds->next;
                ds->obj  = (int *)(&iomod->s_opt->opt_free_mixt_rates);
              }
              else /*! Connect ras struct to an already defined one. Same for
                      opt_alpha & opt_free_mixt_rates */
              {
                ds = parent->ds;

                if (iomod->ras != (t_ras *)ds->obj) Free_RAS(iomod->ras);
                iomod->ras = (t_ras *)ds->obj;

                ds                                = ds->next;
                iomod->s_opt->opt_free_mixt_rates = *((int *)ds->obj);
              }

              rate_value = XML_Get_Attribute_Value(instance, "init.value");

              val = 1.;
              if (rate_value) val = String_To_Dbl(rate_value);

              if (instance->ds->obj == NULL)
              {
                instance->ds->obj       = (phydbl *)(&val);
                instance->ds->next      = (t_ds *)mCalloc(1, sizeof(t_ds));
                instance->ds->next->obj = (int *)(class_num + class_number);

                iomod->ras->gamma_rr->v[class_number] = val;
                iomod->ras->init_rr                   = NO;

                if (Are_Equal(val, 0.0, 1E-20) == NO) class_number++;
              }

              /*! Note: ras is already connected to the relevant t_ds stucture.
                No need to connect ras->gamma_rr or ras->p_invar */

              /*! Invariant */
              if (Are_Equal(val, 0.0, 1E-20))
              {
                mod->ras->invar = YES;
              }
              else
              {
                mod->ras->parent_class_number =
                    *((int *)instance->ds->next->obj);
              }

              xml_node *orig_w = NULL;
              orig_w           = XML_Search_Node_Attribute_Value(
                  "appliesto", instance->id, YES, instance->parent);

              if (orig_w)
              {
                char *weight;
                weight = XML_Get_Attribute_Value(orig_w, "value");
                if (mod->ras->invar == YES)
                {
                  iomod->ras->pinvar->v = String_To_Dbl(weight);
                }
                else
                {
                  int class;
                  class = *((int *)instance->ds->next->obj);
                  iomod->ras->gamma_r_proba->v[class] = String_To_Dbl(weight);
                  iomod->ras->init_r_proba            = NO;
                }
              }
            }

            //////////////////////////////////////////////
            //           BRANCH LENGTHS                 //
            //////////////////////////////////////////////

            else if (!strcmp(parent->name, "branchlengths"))
            {
              int i;
              int n_otu;

              n_otu = tree->n_otu;

              if (instance->ds->obj == NULL)
              {
                if (!lens)
                {
                  ori_lens     = (scalar_dbl **)mCalloc(2 * tree->n_otu - 1,
                                                        sizeof(scalar_dbl *));
                  ori_lens_old = (scalar_dbl **)mCalloc(2 * tree->n_otu - 1,
                                                        sizeof(scalar_dbl *));

                  ori_lens_var     = (scalar_dbl **)mCalloc(2 * tree->n_otu - 1,
                                                            sizeof(scalar_dbl *));
                  ori_lens_var_old = (scalar_dbl **)mCalloc(
                      2 * tree->n_otu - 1, sizeof(scalar_dbl *));

                  lens     = ori_lens;
                  lens_old = ori_lens_old;

                  lens_var     = ori_lens_var;
                  lens_var_old = ori_lens_var_old;

                  lens_size = 2 * tree->n_otu - 1;
                }
                else
                {
                  ori_lens = (scalar_dbl **)mRealloc(
                      ori_lens, 2 * tree->n_otu - 1 + lens_size,
                      sizeof(scalar_dbl *));
                  ori_lens_old = (scalar_dbl **)mRealloc(
                      ori_lens_old, 2 * tree->n_otu - 1 + lens_size,
                      sizeof(scalar_dbl *));

                  ori_lens_var = (scalar_dbl **)mRealloc(
                      ori_lens_var, 2 * tree->n_otu - 1 + lens_size,
                      sizeof(scalar_dbl *));
                  ori_lens_var_old = (scalar_dbl **)mRealloc(
                      ori_lens_var_old, 2 * tree->n_otu - 1 + lens_size,
                      sizeof(scalar_dbl *));

                  lens = ori_lens + lens_size;
                  ;
                  lens_old = ori_lens_old + lens_size;
                  ;

                  lens_var = ori_lens_var + lens_size;
                  ;
                  lens_var_old = ori_lens_var_old + lens_size;
                  ;

                  lens_size += 2 * tree->n_otu - 1;
                }

                for (i = 0; i < 2 * tree->n_otu - 1; ++i)
                {
                  lens[i] = (scalar_dbl *)mCalloc(1, sizeof(scalar_dbl));
                  Init_Scalar_Dbl(lens[i]);

                  lens_old[i] = (scalar_dbl *)mCalloc(1, sizeof(scalar_dbl));
                  Init_Scalar_Dbl(lens_old[i]);

                  lens_var[i] = (scalar_dbl *)mCalloc(1, sizeof(scalar_dbl));
                  Init_Scalar_Dbl(lens_var[i]);

                  lens_var_old[i] =
                      (scalar_dbl *)mCalloc(1, sizeof(scalar_dbl));
                  Init_Scalar_Dbl(lens_var_old[i]);

                  Free_Scalar_Dbl(tree->a_edges[i]->l);
                  Free_Scalar_Dbl(tree->a_edges[i]->l_old);

                  Free_Scalar_Dbl(tree->a_edges[i]->l_var);
                  Free_Scalar_Dbl(tree->a_edges[i]->l_var_old);

                  if (tree->prev &&
                      tree->prev->a_edges[i]->l == mixt_tree->a_edges[i]->l &&
                      tree->prev->is_mixt_tree == NO)
                  {
                    PhyML_Fprintf(stderr, "\n. %p %p", tree->a_edges[i]->l,
                                  mixt_tree->a_edges[i]->l);
                    PhyML_Fprintf(
                        stderr, "\n. Only one set of edge lengths is allowed ");
                    PhyML_Fprintf(
                        stderr,
                        "\n. in a 'partitionelem'. Please fix your XML file.");
                    Exit("\n");
                  }
                }

                char *l_min = NULL;
                l_min       = XML_Get_Attribute_Value(instance, "l.min");

                if (l_min)
                {
                  iomod->l_min = String_To_Dbl(l_min);
                }

                char *opt_bl = NULL;
                opt_bl = XML_Get_Attribute_Value(instance, "optimise.lens");

                if (opt_bl)
                {
                  if (!strcmp(opt_bl, "yes") || !strcmp(opt_bl, "true"))
                  {
                    iomod->s_opt->opt_bl_one_by_one = YES;
                  }
                  else
                  {
                    iomod->s_opt->opt_bl_one_by_one = NO;
                  }
                }

                ds = instance->ds;

                ds->obj = (scalar_dbl **)lens;

                ds->next = (t_ds *)mCalloc(1, sizeof(t_ds));
                ds       = ds->next;
                ds->obj  = (scalar_dbl **)lens_old;

                ds->next = (t_ds *)mCalloc(1, sizeof(t_ds));
                ds       = ds->next;
                ds->obj  = (scalar_dbl **)lens_var;

                ds->next = (t_ds *)mCalloc(1, sizeof(t_ds));
                ds       = ds->next;
                ds->obj  = (scalar_dbl **)lens_var_old;

                ds->next = (t_ds *)mCalloc(1, sizeof(t_ds));
                ds       = ds->next;
                ds->obj  = (int *)(&iomod->s_opt->opt_bl_one_by_one);
              }
              else
              {
                /* Free existing edge lengths, which will be replaced with those
                 * read from <branchlengths> !*/
                for (i = 0; i < 2 * tree->n_otu - 1; ++i)
                {
                  Free_Scalar_Dbl(tree->a_edges[i]->l);
                  Free_Scalar_Dbl(tree->a_edges[i]->l_old);
                  Free_Scalar_Dbl(tree->a_edges[i]->l_var);
                  Free_Scalar_Dbl(tree->a_edges[i]->l_var_old);
                }

                ds = instance->ds;

                lens = (scalar_dbl **)ds->obj;

                ds       = ds->next;
                lens_old = (scalar_dbl **)ds->obj;

                ds       = ds->next;
                lens_var = (scalar_dbl **)ds->obj;

                ds           = ds->next;
                lens_var_old = (scalar_dbl **)ds->obj;

                ds                              = ds->next;
                iomod->s_opt->opt_bl_one_by_one = *((int *)ds->obj);
              }

              if (n_otu != tree->n_otu)
              {
                PhyML_Fprintf(stderr, "\n. All the data sets should display "
                                      "the same number of sequences.");
                PhyML_Fprintf(stderr,
                              "\n. Found at least one data set with %d "
                              "sequences and one with %d sequences.",
                              n_otu, tree->n_otu);
                Exit("\n");
              }

              for (i = 0; i < 2 * tree->n_otu - 1; ++i)
              {
                tree->a_edges[i]->l          = lens[i];
                mixt_tree->a_edges[i]->l     = lens[i];
                tree->a_edges[i]->l_old      = lens_old[i];
                mixt_tree->a_edges[i]->l_old = lens_old[i];

                tree->a_edges[i]->l_var          = lens_var[i];
                mixt_tree->a_edges[i]->l_var     = lens_var[i];
                tree->a_edges[i]->l_var_old      = lens_var_old[i];
                mixt_tree->a_edges[i]->l_var_old = lens_var_old[i];
              }
            }

            ///////////////////////////////////////////////
            ///////////////////////////////////////////////
            ///////////////////////////////////////////////

            if (first_m_elem >
                1) // Done with this component, move to the next tree and model
            {
              if (tree->next) tree = tree->next;
              if (mod->next) mod = mod->next;
            }
          }
          else if (list[j] != ' ')
          {
            component[i] = list[j];
            i++;
          }
          j++;
          if (j == (int)strlen(list) + 1) break;

        } // end of mixtureelem processing
      } // end of partitionelem processing
    } while (1);
  } while (1);

  if (n_partitionelem == 0)
  {
    PhyML_Fprintf(stderr, "\n. No <partitionelem> found. Please amend your XML "
                          "file accordingly.");
    Exit("\n");
  }

  if (ori_lens) Free(ori_lens);
  if (ori_lens_old) Free(ori_lens_old);
  if (ori_lens_var) Free(ori_lens_var);
  if (ori_lens_var_old) Free(ori_lens_var_old);

  while (io->prev != NULL) io = io->prev;
  while (mixt_tree->prev != NULL) mixt_tree = mixt_tree->prev;

  /*! Finish making the models */
  tree = mixt_tree;
  do
  {
    Make_Model_Complete(tree->mod);
    tree = tree->next;
  } while (tree);

  Check_Taxa_Sets(mixt_tree);

#if defined PHYML
  Check_Mandatory_XML_Node(root, "phyml");
#elif defined PHYTIME
  Check_Mandatory_XML_Node(root, "phytime");
#endif

  Check_Mandatory_XML_Node(root, "topology");
  Check_Mandatory_XML_Node(root, "branchlengths");
  Check_Mandatory_XML_Node(root, "ratematrices");
  Check_Mandatory_XML_Node(root, "equfreqs");
  Check_Mandatory_XML_Node(root, "siterates");
  Check_Mandatory_XML_Node(root, "partitionelem");
  Check_Mandatory_XML_Node(root, "mixtureelem");

  if (!io->mod->s_opt->random_input_tree) io->mod->s_opt->n_rand_starts = 1;

  int r_seed = (io->r_seed < 0) ? (time(NULL)) : (io->r_seed);
  srand(r_seed);
  io->r_seed = r_seed;

  Free(component);
  Free(class_num);

  fclose(fp);
  return mixt_tree;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

xml_node *XML_Load_File(FILE *fp)
{
  int       c;
  char     *buffer, *bufptr;
  int       bufsize;
  xml_node *parent, *node;

  buffer = (char *)mCalloc(T_MAX_XML_TAG, sizeof(char));

  bufsize = T_MAX_XML_TAG;
  bufptr  = buffer;
  parent  = NULL;
  node    = NULL;

  while ((c = fgetc(fp)) != EOF)
  {
    if (c == '<' && bufptr > buffer)
    {
      *bufptr = '\0';

      /* PhyML_Printf("\n. Read value '%s' for node '%s'",buffer,node->name); */
      /* fflush(NULL); */

      XML_Set_Node_Value(node, buffer);
      bufptr = buffer;
    }

    if (c == '<')
    {
      bufptr = buffer;

      while ((c = fgetc(fp)) != EOF)
      {
        if (isspace(c) != NO || c == '>' || (c == '/' && bufptr > buffer))
          break; // End of open or close tag
        else if (c == '<')
        {
          Exit("\n. Bare < in element!");
        }
        else if (XML_Add_Character(c, &bufptr, &buffer, &bufsize))
        {
          PhyML_Printf("\n. Err in file %s at line %d\n", __FILE__, __LINE__);
          Exit("\n");
        }
      }

      *bufptr = '\0';

      if (!strcmp(buffer, "!--")) // Get the rest of the comment
      {
        while ((c = fgetc(fp)) != EOF)
        {

          if (c == '>' && bufptr > (buffer + 4) && bufptr[-3] != '-' &&
              bufptr[-2] == '-' && bufptr[-1] == '-')
            break;
          else if (XML_Add_Character(c, &bufptr, &buffer, &bufsize))
          {
            PhyML_Fprintf(stderr, "\n. Err in file %s at line %d\n", __FILE__,
                          __LINE__);
            Exit("\n");
          }
        }
        *bufptr = '\0';

        if (c != '>')
        {
          PhyML_Fprintf(stderr, "\n. Early EOF in comment node.");
          Exit("\n");
        }
      }
      else if (buffer[0] == '/') // Close tag
      {
        if (strcmp(buffer + 1, parent->name))
        {
          PhyML_Fprintf(
              stderr,
              "\n. Opened tag with name '%s' and closed it with '%s'...",
              node->name, buffer + 1);
          PhyML_Fprintf(stderr, "\n. Err. in file %s at line %d\n", __FILE__,
                        __LINE__);
          Exit("\n");
        }

        /* printf("\n. Closing node with name '%s'",node->name); */

        if (node->parent)
        {
          parent = parent->parent;
          node   = parent;
        }
      }
      else if (buffer[0] == '?')
      {
        while ((c = fgetc(fp)) != EOF)
        {
          if (c == '>' && bufptr > buffer && bufptr[-1] == '?')
            break;
          else if (XML_Add_Character(c, &bufptr, &buffer, &bufsize))
          {
            PhyML_Fprintf(stderr, "\n. Err in file %s at line %d\n", __FILE__,
                          __LINE__);
            Exit("\n");
          }
        }

        if (c != '>')
        {
          PhyML_Fprintf(
              stderr,
              "\n. An error occurred when reading the processing instruction.");
          PhyML_Fprintf(stderr, "\n. Err in file %s at line %d\n", __FILE__,
                        __LINE__);
          Exit("\n");
        }

        *bufptr = '\0';
      }
      else // Open tag
      {
        node = XML_Make_Node(buffer);
        XML_Init_Node(parent, node, buffer);
        if (!parent) parent = node;

        if (isspace(c) != NO)
          c = XML_Parse_Element(fp, node);
        else if (c == '/')
        {
          if ((c = fgetc(fp)) != '>')
          {
            PhyML_Fprintf(stderr, "\n. Expected '>' but read '%c' instead", c);
            Exit("\n");
          }
          c = '/';
        }

        if (c != '/') parent = node;

        buffer[0] = '\0';
      }
      bufptr = buffer;
    }
    else if (isspace(c) == NO)
    {
      if (XML_Add_Character(c, &bufptr, &buffer, &bufsize))
      {
        PhyML_Fprintf(stderr, "\n. Err in file %s at line %d\n", __FILE__,
                      __LINE__);
        Exit("\n");
      }
    }
  }
  Free(buffer);
  return node;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int XML_Add_Character(int c, char **bufptr, char **buffer, int *bufsize)
{
  char *newbuffer;

  if (*bufptr >= (*buffer + *bufsize - 4))
  {
    // Increase the size of the buffer...
    int delta;

    delta = (*bufptr - *buffer);

    if (*bufsize < 1024)
      (*bufsize) *= 2;
    else
      (*bufsize) += 1024;

    if ((newbuffer = realloc(*buffer, *bufsize)) == NULL)
    {
      Free(*buffer);
      PhyML_Fprintf(stderr, "\n. Unable to expand string buffer to %d bytes!",
                    *bufsize);
      Exit("\n");
    }

    *bufptr = newbuffer + delta;
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
  int   c;
  int   quote;
  char *name, *value, *ptr;
  int   namesize, valsize;

  name  = (char *)mCalloc(64, sizeof(char));
  value = (char *)mCalloc(64, sizeof(char));

  namesize = 64;
  valsize  = 64;

  while ((c = fgetc(fp)) != EOF)
  {

    if (isspace(c) != NO) continue;

    if (c == '/') // End of tag
    {
      /* printf("\n. Closing node '%s'.",n->name); */

      quote = fgetc(fp);
      if (quote != '>')
      {
        PhyML_Fprintf(stderr,
                      "\n. Expected '>' after '%c' but read '%c' instead", c,
                      quote);
        Exit("\n");
      }
      break;
    }
    else if (c == '<')
    {
      Exit("\n. Bare < in element!");
    }
    else if (c == '>') // End of tag
    {
      break;
    }

    name[0] = c;
    ptr     = name + 1;

    if (c == '\"' || c == '\'') // Name is in quotes
    {
      quote = c;

      while ((c = fgetc(fp)) != EOF)
      {
        if (XML_Add_Character(c, &ptr, &name, &namesize))
        {
          PhyML_Fprintf(stderr, "\n. Err in file %s at line %d\n", __FILE__,
                        __LINE__);
          Exit("\n");
        }
        if (c == quote) break;
      }
    }
    else // Name not in quotes
    {
      while ((c = fgetc(fp)) != EOF)
      {
        if (isspace(c) != NO || c == '=' || c == '/' || c == '>' || c == '?')
          break;
        else
        {
          if (XML_Add_Character(c, &ptr, &name, &namesize))
          {
            PhyML_Fprintf(stderr, "\n. Err in file %s at line %d\n", __FILE__,
                          __LINE__);
            Exit("\n");
          }
        }
      }
    }

    *ptr = '\0';

    while (c != EOF && isspace(c) != NO) c = fgetc(fp);

    if (c == '=') // Read the attribute value
    {
      while ((c = fgetc(fp)) != EOF && isspace(c) != NO);

      if (c == EOF)
      {
        PhyML_Fprintf(stderr, "\n. Missing value in attribute.");
        PhyML_Fprintf(stderr, "\n. Err in file %s at line %d\n", __FILE__,
                      __LINE__);
        Exit("\n");
      }

      if (c == '\'' || c == '\"')
      {
        quote = c;
        ptr   = value;

        while ((c = fgetc(fp)) != EOF)
        {
          if (c == quote)
            break;
          else
          {
            if (XML_Add_Character(c, &ptr, &value, &valsize))
            {
              PhyML_Fprintf(stderr, "\n. Err in file %s at line %d\n", __FILE__,
                            __LINE__);
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

        while ((c = fgetc(fp)) != EOF)
        {
          if (isspace(c) != NO || c == '=' || c == '/' || c == '>')
            break;
          else
          {
            if (XML_Add_Character(c, &ptr, &value, &valsize))
            {
              PhyML_Fprintf(stderr, "\n. Err in file %s at line %d\n", __FILE__,
                            __LINE__);
              Exit("\n");
            }
          }
        }
      }
    }

    /* printf("\n. Setting attribute '%s=%s' to node '%s'",name,value,n->name);
     */
    XML_Make_And_Set_Attribute(n, name, value);

    if (c == '>') break;
  }
  Free(name);
  Free(value);

  /* printf("\n. Return '%c'\n",c); */
  return (c);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

xml_attr *XML_Search_Attribute(xml_node *n, char *target_attr_name)
{
  xml_attr *attr;

  attr = n->attr;
  do
  {
    if (!strcmp(attr->name, target_attr_name)) return attr;
    attr = attr->next;
  } while (attr);

  return (NULL);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int XML_Make_And_Set_Attribute(xml_node *n, char *attr_name, char *attr_value)
{
  xml_attr *prev;
  char     *s;

  prev = NULL;
  while (n->attr != NULL)
  {
    prev    = n->attr;
    n->attr = n->attr->next;
  }

  n->attr = XML_Make_Attribute(prev, attr_name, attr_value);
  XML_Init_Attribute(n->attr);
  n->n_attr++;

  // rewind
  while (n->attr->prev != NULL) n->attr = n->attr->prev;

  s = To_Lower_String(attr_name);
  if (!strcmp(s, "id"))
  {
    XML_Set_Node_Id(n, attr_value);
    /* printf("\n. Node '%s' id is '%s'",n->name,n->id); */
  }
  Free(s);

  return (0);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int XML_Set_Attribute_Value(xml_node *n, char *attr_name, char *attr_value)
{
  xml_attr *attr;
  char     *s;

  attr = XML_Search_Attribute(n, attr_name);
  if (attr == NULL) return (-1);

  s = To_Lower_String(attr_value);
  strcpy(attr->value, s);

  Free(s);

  return (0);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int XML_Add_Attribute(xml_node *n, char *attr_name, char *attr_value)
{
  xml_attr *attr;

  attr = n->attr;
  while (attr->next != NULL) attr = attr->next;
  assert(attr);

  attr->next = XML_Make_Attribute(attr, attr_name, attr_value);
  XML_Init_Attribute(attr);
  n->n_attr++;

  return (0);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

xml_node *XML_Add_Node(xml_node *parent, char *nd_name)
{
  xml_node *n;
  n = XML_Make_Node(nd_name);
  XML_Init_Node(parent, n, nd_name);
  return (n);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int XML_Set_Node_Id(xml_node *n, char *id)
{
  XML_Make_Node_Id(n, id);
  strcpy(n->id, id);
  return (0);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int XML_Set_Node_Value(xml_node *n, char *val)
{
  XML_Make_Node_Value(n, val);
  strcpy(n->value, val);
  return (0);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

xml_node *XML_Search_Node_Generic(char *nd_name, char *attr_name,
                                  char *attr_val, int skip, xml_node *node)
{

  xml_node *match;

  /* if(nd_name) printf("\n. [1] nd_name:%s attr_name:%s attr_val:%s \n",
   * nd_name, attr_name, attr_val); */
  /* else  printf("\n. attr_name:%s attr_val:%s \n", attr_name, attr_val); */

  /* printf("\n. name:%s child:%s next:%s ", */
  /*     node?node->name:"xx", */
  /*     node->child?node->child->name:"xx", */
  /*     node->next?node->next->name:"xx"); fflush(NULL); */

  match = NULL;
  if (skip == NO && nd_name && attr_name && attr_val)
  {
    if (!strcmp(nd_name, node->name))
    {
      xml_attr *attr = XML_Search_Attribute(node, attr_name);
      if (attr && !strcmp(attr->value, attr_val)) match = node;
    }
  }
  else if (skip == NO && !nd_name && attr_name && attr_val)
  {
    xml_attr *attr = XML_Search_Attribute(node, attr_name);
    if (attr && !strcmp(attr->value, attr_val)) match = node;
  }
  else if (skip == NO && nd_name && !attr_name && attr_val)
  {
    if (!strcmp(nd_name, node->name))
    {
      do
      {
        if (!strcmp(node->attr->value, attr_val))
        {
          match = node;
          break;
        }
        node->attr = node->attr->next;
        if (!node->attr) break;
      } while (1);
    }
  }
  else if (skip == NO && nd_name && attr_name && !attr_val)
  {
    if (!strcmp(nd_name, node->name))
    {
      do
      {
        if (!strcmp(node->attr->name, attr_name))
        {
          match = node;
          break;
        }
        node->attr = node->attr->next;
        if (!node->attr) break;
      } while (1);
    }
  }
  else if (skip == NO && nd_name && !attr_name && !attr_val)
  {
    if (!strcmp(nd_name, node->name)) match = node;
  }
  else if (skip == NO && !nd_name && attr_name && !attr_val)
  {
    xml_attr *attr = XML_Search_Attribute(node, attr_name);
    if (attr) match = node;
  }
  else if (skip == NO && !nd_name && !attr_name && attr_val)
  {
    do
    {
      if (!strcmp(node->attr->value, attr_val))
      {
        match = node;
        break;
      }
      node->attr = node->attr->next;
      if (!node->attr) break;
    } while (1);
  }

  // If node has a child, node = child, else if node has next, node = next, else
  // if node has parent, node = parent->next else node = NULL

  if (match) return (match);
  if (node->child)
  {
    match =
        XML_Search_Node_Generic(nd_name, attr_name, attr_val, NO, node->child);
  }
  if (!match && node->next)
  {
    match =
        XML_Search_Node_Generic(nd_name, attr_name, attr_val, NO, node->next);
  }
  if (match == NULL && node->parent)
  {
    if (node->parent == NULL) // Reached the root
    {
      PhyML_Fprintf(stderr, "\n. Could not find a node with name '%s'.",
                    attr_name);
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
  if (skip == NO && !strcmp(node->name, name))
    match = node;
  else
  {
    // If node has a child, node = child, else if node has next, node = next,
    // else if node has parent, node = parent->next else node = NULL
    if (node->child)
    {
      match = XML_Search_Node_Name(name, NO, node->child);
    }
    if (match == NULL && node->next)
    {
      match = XML_Search_Node_Name(name, NO, node->next);
    }
    if (match == NULL && node->parent)
    {
      if (node->parent == NULL) // Reached the root
      {
        PhyML_Fprintf(stderr, "\n. Could not find a node with name '%s'.",
                      name);
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

  if (!node)
  {
    PhyML_Fprintf(stderr, "\n. Err. in file %s at line %d\n", __FILE__,
                  __LINE__);
    Exit("\n");
  }

  match = NULL;
  if (skip == NO && node->id && !strcmp(node->id, id))
    match = node;
  else
  {
    // If node has a child, node = child, else if node has next, node = next,
    // else if node has parent, node = parent->next else node = NULL
    if (node->child)
    {
      match = XML_Search_Node_ID(id, NO, node->child);
    }
    if (match == NULL && node->next)
    {
      match = XML_Search_Node_ID(id, NO, node->next);
    }
    if (match == NULL && node->parent)
    {
      if (node->parent == NULL) // Reached the root
      {
        PhyML_Fprintf(stderr, "\n. Could not find a node with id '%s'.", id);
        Exit("\n");
      }

      return NULL;
    }
  }
  return match;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

xml_node *XML_Search_Node_Attribute_Value(char *attr_name, char *value,
                                          int skip, xml_node *node)
{
  xml_node *match;

  if (!node)
  {
    PhyML_Fprintf(stderr, "\n. node: %p attr: %p", node,
                  node ? node->attr : NULL);
    PhyML_Fprintf(stderr, "\n. Err in file %s at line %d\n", __FILE__,
                  __LINE__);
    Exit("\n");
  }

  match = NULL;

  if (skip)
  {
    match = XML_Search_Node_Attribute_Value(attr_name, value, NO, node->child);
    return match;
  }

  if (skip == NO && node->attr)
  {
    xml_attr *attr;
    char     *sname, *sval;
    char     *qname, *qval;

    qname = To_Lower_String(attr_name);
    qval  = To_Lower_String(value);

    attr = node->attr;
    do
    {
      sname = To_Lower_String(attr->name);
      sval  = To_Lower_String(attr->value);

      if (!strcmp(sname, qname) && !strcmp(sval, qval))
      {
        match = node;
        Free(sname);
        Free(sval);
        break;
      }

      Free(sname);
      Free(sval);

      attr = attr->next;

      if (!attr) break;
    } while (1);

    Free(qval);
    Free(qname);
  }

  if (match)
    return (match);
  else
  {
    if (node->child)
    {
      match =
          XML_Search_Node_Attribute_Value(attr_name, value, NO, node->child);
    }
    if (node->next && !match)
    {
      match = XML_Search_Node_Attribute_Value(attr_name, value, NO, node->next);
    }
  }
  return match;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

char *XML_Get_Attribute_Value(xml_node *node, char *attr_name)
{
  xml_attr *attr;

  attr = node->attr;

  while (attr && strcmp(attr->name, attr_name)) attr = attr->next;

  return (attr ? attr->value : NULL);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int XML_Validate_Attr_Int(char *target, int num, ...)
{
  va_list args;
  int     i;
  char   *s, *sc_s;
  char   *sc_target;

  sc_target = To_Lower_String(target);

  va_start(args, num);
  for (i = 0; i < num; i++)
  {
    s    = va_arg(args, char *);
    sc_s = To_Lower_String(s);
    if (!strcmp(sc_s, sc_target))
    {
      Free(sc_s);
      break;
    }
    Free(sc_s);
  }
  va_end(args);

  if (i == num)
  {
    i = -1;
    PhyML_Fprintf(stderr, "\n. Attribute value '%s' is not valid", target);
    Exit("\n");
  }

  Free(sc_target);

  return (i);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void XML_Check_Siterates_Node(xml_node *parent)
{
  xml_node *n;
  int       n_weights_nodes;
  char     *rate_value = NULL;
  phydbl    buff;
  int       n_zeros;
  char     *endptr;

  if (!parent)
  {
    PhyML_Fprintf(stderr, "\n. Err in file %s at line %d\n", __FILE__,
                  __LINE__);
    Exit("\n");
  }

  if (strcmp(parent->name, "siterates"))
  {
    PhyML_Fprintf(stderr, "\n. Node name '%s' should be 'siterates'",
                  parent->name);
    Exit("\n");
  }

  // Check that only one 'weights' node is present
  n_weights_nodes = 0;
  n               = parent->child;
  do
  {
    if (!strcmp(n->name, "weights")) n_weights_nodes++;
    if (n_weights_nodes > 1)
    {
      PhyML_Fprintf(
          stderr,
          "\n. Only one distribution is authorized for 'siterates' nodes.");
      Exit("\n");
    }
    n = n->next;
    if (!n) break;
  } while (1);

  // Check that one rate value is set to zero if gamma+inv model is used
  n = XML_Search_Node_Attribute_Value("family", "gamma+inv", YES, parent);
  if (!n)
    return;
  else
  {
    n_zeros = 0;
    n       = parent->child;
    do
    {
      if (!strcmp(n->name, "instance"))
      {
        rate_value = NULL;
        rate_value = XML_Get_Attribute_Value(n, "init.value");

        if (rate_value)
        {
          errno = 0;
          buff  = strtod(rate_value, &endptr);

          if (rate_value == endptr || errno == ERANGE)
          {
            PhyML_Fprintf(stderr, "\n. value: %s", rate_value);
            PhyML_Fprintf(stderr, "\n. Error in reading attribute 'init.value' "
                                  "in node 'instance'.");
            Exit("\n");
          }

          if (buff < 1.E-20) n_zeros++;
        }
      }
      n = n->next;
      if (!n) break;
    } while (1);

    if (n_zeros != 1)
    {
      PhyML_Fprintf(stderr, "\n. # of zero-rates: %d", n_zeros);
      PhyML_Fprintf(stderr, "\n. Exactly one rate value has to be set to zero "
                            "when using the 'gamma+inv' model.");
      PhyML_Fprintf(stderr, "\n. Component id: %s", parent->id);
      Exit("\n");
    }
  }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int XML_Get_Number_Of_Classes_Siterates(xml_node *parent)
{
  xml_node *n;
  int       n_classes;

  if (!parent)
  {
    PhyML_Printf("\n. Err. in file %s at line %d\n", __FILE__, __LINE__);
    Generic_Exit(__FILE__, __LINE__, __FUNCTION__);
  }

  n_classes = 0;
  n         = parent->child;
  do
  {
    if (!strcmp(n->name, "instance")) n_classes++;
    n = n->next;
    if (!n) break;
  } while (1);

  n = NULL;
  n = XML_Search_Node_Attribute_Value("family", "gamma+inv", YES, parent);

  if (!n)
    return n_classes;
  else
    return n_classes - 1;
}

//////////////////////////////////////////////////////////////

int XML_Siterates_Has_Invariants(xml_node *parent)
{
  xml_node *n;

  if (!parent)
  {
    PhyML_Fprintf(stderr, "\n. Err in file %s at line %d\n", __FILE__,
                  __LINE__);
    Exit("\n");
  }

  n = NULL;
  n = XML_Search_Node_Attribute_Value("family", "gamma+inv", YES, parent);

  if (!n)
    return NO;
  else
    return YES;
}

//////////////////////////////////////////////////////////////

void XML_Count_Number_Of_Node_With_ID(char *id, int *count, xml_node *n)
{
  if (!id) return;
  if (n->id && !strcmp(n->id, id)) (*count)++;

  if (n->child) XML_Count_Number_Of_Node_With_ID(id, count, n->child);
  if (n->next) XML_Count_Number_Of_Node_With_ID(id, count, n->next);
}

//////////////////////////////////////////////////////////////

void XML_Count_Number_Of_Node_With_Name(char *name, int *count, xml_node *n)
{
  if (!name) return;
  if (n->name && !strcmp(n->name, name)) (*count)++;

  if (n->child) XML_Count_Number_Of_Node_With_Name(name, count, n->child);
  if (n->next) XML_Count_Number_Of_Node_With_Name(name, count, n->next);
}

//////////////////////////////////////////////////////////////

void XML_Check_Duplicate_ID(xml_node *n)
{
  int count;

  count = 0;
  XML_Count_Number_Of_Node_With_ID(n->id, &count, n);

  if (count > 1)
  {
    PhyML_Fprintf(stderr, "\n. Node ID '%s' was found more than once.", n->id);
    PhyML_Fprintf(stderr, "\n. Each ID must be unique. Please amend your XML");
    PhyML_Fprintf(stderr, "\n. file accordingly.");
    Exit("\n");
  }

  if (n->child) XML_Check_Duplicate_ID(n->child);
  if (n->next) XML_Check_Duplicate_ID(n->next);
}

//////////////////////////////////////////////////////////////

xml_node *XML_Copy_XML_Graph(xml_node *root)
{
  xml_node *cpy_root;

  cpy_root = XML_Make_Node(root->name);
  XML_Copy_XML_Node(cpy_root, root);

  return (cpy_root);
}

//////////////////////////////////////////////////////////////

void XML_Copy_XML_Node(xml_node *cpy_root, xml_node *root)
{
  xml_attr *attr, *cpy_attr;

  strcpy(cpy_root->name, root->name);

  XML_Make_Node_Id(cpy_root, root->id);
  if (root->id) strcpy(cpy_root->id, root->id);

  XML_Make_Node_Value(cpy_root, root->value);
  if (root->value) strcpy(cpy_root->value, root->value);

  cpy_root->n_attr = root->n_attr;

  if (root->attr)
  {
    cpy_root->attr =
        XML_Make_Attribute(NULL, root->attr->name, root->attr->value);
    XML_Init_Attribute(cpy_root->attr);
    attr     = root->attr;
    cpy_attr = cpy_root->attr;
    while (attr->next)
    {
      fflush(NULL);
      cpy_attr->next =
          XML_Make_Attribute(cpy_attr, attr->next->name, attr->next->value);
      XML_Init_Attribute(cpy_attr->next);
      attr     = attr->next;
      cpy_attr = cpy_attr->next;
    }
  }

  if (root->child)
  {
    cpy_root->child         = XML_Make_Node(root->child->name);
    cpy_root->child->parent = cpy_root;
    XML_Copy_XML_Node(cpy_root->child, root->child);
  }

  if (root->next)
  {
    cpy_root->next       = XML_Make_Node(root->next->name);
    cpy_root->next->prev = cpy_root;
    XML_Copy_XML_Node(cpy_root->next, root->next);
  }
}

//////////////////////////////////////////////////////////////

void XML_Write_XML_Graph(FILE *fp, xml_node *root)
{
  int indent;
  indent = 0;
  XML_Write_XML_Node(fp, &indent, YES, YES, root);
}

//////////////////////////////////////////////////////////////

void XML_Write_XML_Node(FILE *fp, int *indent, short int recurr_next,
                        short int recurr_child, xml_node *root)
{
  xml_node *n;
  xml_attr *attr;
  char     *s;
  int       i;

  assert(root);
  
  s = (char *)mCalloc((*indent) + 1, sizeof(char));
  for (i = 0; i < (*indent); ++i) s[i] = '\t';
  s[i] = '\0';

  PhyML_Fprintf(fp, "\n%s", s);

  n = root;

  PhyML_Fprintf(fp, "<%s", n->name);

  attr = n->attr;
  while (attr)
  {
    PhyML_Fprintf(fp, " %s=\"%s\"", attr->name, attr->value);
    attr = attr->next;
  }
  PhyML_Fprintf(fp, ">");

  if (n->value != NULL) XML_Write_Node_Value(fp, s, n);

  if (n->child && recurr_child == YES)
  {
    (*indent)++;
    XML_Write_XML_Node(fp, indent, recurr_next, recurr_child, n->child);
    (*indent)--;
  }

  PhyML_Fprintf(fp, "\n%s</%s>", s, n->name);

  if (n->next && recurr_next == YES)
    XML_Write_XML_Node(fp, indent, recurr_next, recurr_child, n->next);

  Free(s);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void XML_Write_Node_Value(FILE *fp, char *indent, xml_node *n)
{
  PhyML_Fprintf(fp, "\n");

  char *tk = strtok(n->value, "\n");

  do
  {
    PhyML_Fprintf(fp, "%s%s", indent, tk);
    tk = strtok(NULL, "\n");
    if (tk != NULL) PhyML_Fprintf(fp, "\n");
  } while (tk != NULL);
}
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Check_Mandatory_XML_Node(xml_node *root, char *name)
{
  if (!XML_Search_Node_Name(name, NO, root))
  {
    PhyML_Fprintf(
        stderr, "\n. Could not find mandatory XML node with name '%s'.", name);
    PhyML_Fprintf(stderr, "\n. Please amend your XML file.");
    Exit("\n");
  }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int XML_Number_Of_Taxa_In_Clade(xml_node *n_clade)
{
  int clade_size = 0;
  if (n_clade)
  {
    do
    {
      clade_size++;
      if (n_clade->next)
        n_clade = n_clade->next;
      else
        break;
    } while (n_clade);
  }
  else
  {
    PhyML_Fprintf(stderr, "\n. Clade is empty.");
    PhyML_Fprintf(stderr, "\n. Err. in file %s at line %d\n", __FILE__,
                  __LINE__);
    Exit("\n");
  }
  return (clade_size);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

char **XML_Read_Clade(xml_node *xnd_clade, t_tree *tree)
{
  int    i;
  char **clade;

  clade = (char **)mCalloc(tree->n_otu, sizeof(char *));

  if (xnd_clade)
  {
    i = 0;
    do
    {
      clade[i] = xnd_clade->attr->value;
      i++;
      if (xnd_clade->next)
        xnd_clade = xnd_clade->next;
      else
        break;
    } while (xnd_clade);
  }
  else
  {
    PhyML_Fprintf(stderr, "== Clade is empty. \n");
    PhyML_Fprintf(stderr, "\n. Err. in file %s at line %d\n", __FILE__,
                  __LINE__);
    Exit("\n");
  }

  return (clade);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void XML_Read_Calibration(xml_node *xroot, t_tree *tree)
{
  xml_node *xnd_clade, *xnd_cal, *xnd, *xnd_dum;
  phydbl    low, up, alpha_proba_dbl, t, t_ref;
  t_cal    *cal;
  t_clad   *clade;
  char     *calib_id, *clade_id, *clade_name, *alpha_proba_string;
  int       i, j;

  clade = NULL;
  cal   = NULL;

  xnd = xroot->child;
  assert(xnd);

  do
  {
    if (!strcmp(xnd->name, "calibration")) // Found a XML node <calibration>.
    {
      // TO DO: make sure calibs are shared across partition elements -> need to
      // write chain function to call once the calib struct on the first
      // mixt_tree is initialized.
      /* mixt_tree->rates->tot_num_cal++; */
      /* if (mixt_tree->rates->calib == NULL) mixt_tree->rates->calib =
       * Make_Calib(mixt_tree->n_otu); */

      xnd_cal = xnd;

      low = 0.0;
      up  = BIG;

      cal = Make_Calibration();
      Init_Calibration(cal);

      xnd_dum = XML_Search_Node_Name("lower", YES, xnd_cal);
      if (xnd_dum != NULL) low = String_To_Dbl(xnd_dum->value);

      xnd_dum = XML_Search_Node_Name("upper", YES, xnd_cal);
      if (xnd_dum != NULL) up = String_To_Dbl(xnd_dum->value);

      calib_id = XML_Get_Attribute_Value(xnd_cal, "id");
      cal->id  = (char *)mCalloc(strlen(calib_id) + 1, sizeof(char));
      strcpy(cal->id, calib_id);

      cal->clade_list_size   = 0;
      cal->current_clade_idx = 0;
      cal->lower             = low;
      cal->upper             = up;
      cal->is_primary        = YES;

      tree->times->a_cal[tree->times->n_cal] = cal;
      tree->times->n_cal++;

      if (tree->times->n_cal > 10 * tree->n_otu)
      {
        PhyML_Fprintf(stderr, "\n. There are too many clades defined that are "
                              "not found in the tree...");
        PhyML_Fprintf(stderr, "\n. Please remove some of them.");
        Generic_Exit(__FILE__, __LINE__, __FUNCTION__);
      }

      xnd_dum = xnd_cal->child;
      do
      {
        clade_name = NULL;
        if (!strcmp("appliesto", xnd_dum->name))
        {
          clade_name = XML_Get_Attribute_Value(xnd_dum, "clade.id");

          if (!clade_name)
          {
            PhyML_Fprintf(stderr,
                          "\n. Attribute 'value=CLADE_NAME' is mandatory");
            PhyML_Fprintf(stderr,
                          "\n. Please amend your XML file accordingly.");
            Exit("\n");
          }

          alpha_proba_string = XML_Get_Attribute_Value(xnd_dum, "probability");

          alpha_proba_dbl = -1;
          if (!alpha_proba_string)
            alpha_proba_dbl = 1.0;
          else
            alpha_proba_dbl = String_To_Dbl(alpha_proba_string);
          assert(!(alpha_proba_dbl < 0. || alpha_proba_dbl > 1.));

          if (strcmp("root", clade_name))
          {
            xnd_clade =
                XML_Search_Node_Generic("clade", "id", clade_name, YES, xroot);

            if (xnd_clade != NULL) // found clade with a given name
            {
              char **xclade;
              int    clade_size, nd_num;
              int    i;

              xclade     = XML_Read_Clade(xnd_clade->child, tree);
              clade_size = XML_Number_Of_Taxa_In_Clade(xnd_clade->child);

              clade = Make_Clade();

              clade->n_tax    = clade_size;
              clade->tax_list = (char **)mCalloc(clade_size, sizeof(char *));
              for (i = 0; i < clade_size; ++i)
                clade->tax_list[i] =
                    (char *)mCalloc(strlen(xclade[i]) + 1, sizeof(char));
              for (i = 0; i < clade_size; ++i)
                strcpy(clade->tax_list[i], xclade[i]);

              clade_id = XML_Get_Attribute_Value(xnd_dum, "clade.id");
              if (clade_id == NULL)
              {
                PhyML_Fprintf(stderr, "\n. Attribute \"clade.id\" is missing "
                                      "in \"appliesto\" tag.");
                Generic_Exit(__FILE__, __LINE__, __FUNCTION__);
              }
              clade->id = (char *)mCalloc(strlen(clade_id) + 1, sizeof(char));
              strcpy(clade->id, clade_id);

              nd_num = Find_Clade(clade->tax_list, clade_size, tree);

              if (nd_num > -1)
              {
                clade->target_nd = tree->a_nodes[nd_num];
                clade->tip_list  = Make_Target_Tip(clade->n_tax);
                Init_Target_Tip(clade, tree);

                if (cal->clade_list_size == 0)
                  cal->clade_list = (t_clad **)mCalloc(1, sizeof(t_clad *));
                else
                  cal->clade_list = (t_clad **)mRealloc(
                      cal->clade_list, cal->clade_list_size + 1,
                      sizeof(t_clad *));
                if (cal->clade_list_size == 0)
                  cal->alpha_proba_list = (phydbl *)mCalloc(1, sizeof(phydbl));
                else
                  cal->alpha_proba_list = (phydbl *)mRealloc(
                      cal->alpha_proba_list, cal->clade_list_size + 1,
                      sizeof(phydbl));

                cal->clade_list[cal->clade_list_size]       = clade;
                cal->alpha_proba_list[cal->clade_list_size] = alpha_proba_dbl;
                cal->clade_list_size++;
              }
              else
              {
                Free_Clade(clade);
              }

              Free(xclade);
            }
            else
            {
              PhyML_Fprintf(
                  stderr,
                  "\n. Calibration information for clade [%s] was not found.",
                  clade_name);
              PhyML_Fprintf(stderr, "\n. Err. in file %s at line %d\n",
                            __FILE__, __LINE__);
              Exit("\n");
            }
          }
        }
        xnd_dum = xnd_dum->next;
      } while (xnd_dum != NULL);
      if (clade_name == NULL)
      {
        PhyML_Fprintf(
            stderr,
            "\n. Could not find calibration information for calibration `%s'",
            cal->id);
        PhyML_Fprintf(stderr, "\n. Please amend your XML file.");
        Exit("\n");
      }
    }
    xnd = xnd->next;
  } while (xnd != NULL);

  PhyML_Printf("\n\n..........................................................."
               "............");
  for (i = 0; i < tree->times->n_cal; ++i)
  {
    cal = tree->times->a_cal[i];

    phydbl sum = 0.0;
    for (j = 0; j < cal->clade_list_size; j++) sum += cal->alpha_proba_list[j];
    for (j = 0; j < cal->clade_list_size; j++) cal->alpha_proba_list[j] /= sum;

    for (j = 0; j < cal->clade_list_size; j++)
    {
      clade = cal->clade_list[j];
      PhyML_Printf("\n Calibration id: %s.", cal->id);
      PhyML_Printf("\n Lower bound set to: %15f time units.", cal->lower);
      PhyML_Printf("\n Upper bound set to: %15f time units.", cal->upper);
      PhyML_Printf(
          "\n This calibration applies to node %d with probability %G.",
          clade->target_nd->num, cal->alpha_proba_list[j]);
      PhyML_Printf("\n........................................................."
                   "..............");
    }
  }

  t_ref = t = -1.;

  for (i = 0; i < tree->times->n_cal; ++i)
  {
    cal = tree->times->a_cal[i];
    for (j = 0; j < cal->clade_list_size; j++)
    {
      clade = cal->clade_list[j];
      if (clade->target_nd->tax == YES)
      {
        if (t_ref < -.5)
        {
          t_ref = Uni() * (cal->upper - cal->lower) + cal->lower;
          t_ref = fabs(t_ref);
        }
        else
        {
          t = Uni() * (cal->upper - cal->lower) + cal->lower;
          t = fabs(t);

          if (Are_Equal(t, t_ref, 1.E-10) == NO)
          {
            tree->times->is_asynchronous = YES;
            break;
          }
        }
      }
    }
  }

  PhyML_Printf("\n\n");
  PhyML_Printf("\n. Number of calibrations: %d", tree->times->n_cal);
  PhyML_Printf("\n. Is asynchronous: %s",
               tree->times->is_asynchronous ? "yes" : "no");
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void XML_Update_XML_Struct_Given_Model_Params(t_tree *tree)
{
  xml_node *parent, *child, *x;
  char     *s, *val;
  t_ds     *ds;
  phydbl    v;

  val = (char *)mCalloc(T_MAX_LINE, sizeof(char));

  parent = XML_Search_Node_Name("ratematrices", YES, tree->xml_root);
  if (parent != NULL)
  {
    child = parent->child;

    assert(child != NULL);

    do
    {
      ds = child->ds;

      assert(ds);

      x = XML_Search_Node_Name("rr", YES, child);

      if (x != NULL)
      {
        s = XML_Get_Attribute_Value(x, "AC");
        if (s != NULL)
        {
          sprintf(val, "%f", ((t_rmat *)(ds->obj))->rr->v[0]);
          XML_Set_Attribute_Value(x, "AC", val);
        }

        s = XML_Get_Attribute_Value(x, "AG");
        if (s != NULL)
        {
          sprintf(val, "%f", ((t_rmat *)(ds->obj))->rr->v[1]);
          XML_Set_Attribute_Value(x, "AC", val);
        }

        s = XML_Get_Attribute_Value(x, "AT");
        if (s != NULL)
        {
          sprintf(val, "%f", ((t_rmat *)(ds->obj))->rr->v[2]);
          XML_Set_Attribute_Value(x, "AT", val);
        }

        s = XML_Get_Attribute_Value(x, "CG");
        if (s != NULL)
        {
          sprintf(val, "%f", ((t_rmat *)(ds->obj))->rr->v[3]);
          XML_Set_Attribute_Value(x, "CG", val);
        }

        s = XML_Get_Attribute_Value(x, "CT");
        if (s != NULL)
        {
          sprintf(val, "%f", ((t_rmat *)(ds->obj))->rr->v[4]);
          XML_Set_Attribute_Value(x, "CT", val);
        }
      }

      ds = ds->next;
      assert(ds);

      s = XML_Get_Attribute_Value(child, "tstv");
      if (s != NULL)
      {
        sprintf(val, "%f", ((scalar_dbl *)(ds->obj))->v);
        XML_Set_Attribute_Value(child, "tstv", val);
      }

      child = child->next;
    } while (child);
  }

  parent = XML_Search_Node_Name("siterates", YES, tree->xml_root);

  if (parent != NULL)
  {
    child = parent->child;
    assert(child != NULL);

    int class = 0;
    do
    {
      s = XML_Get_Attribute_Value(child, "init.value");
      if (s != NULL)
      {
        v = atof(s);
        if (Are_Equal(v, 0.0, 1E-20) == NO)
        {
          sprintf(val, "%f", ((t_ras *)(parent->ds->obj))->gamma_rr->v[class]);
          XML_Set_Attribute_Value(child, "init.value", val);
          class ++;
        }
      }

      child = child->next;
    } while (child);

    child = XML_Search_Node_Name("weights", YES, parent);
    assert(child != NULL);
    s = XML_Get_Attribute_Value(child, "family");
    assert(s != NULL);
    if (!strcmp(s, "gamma"))
    {
      s = XML_Get_Attribute_Value(child, "alpha");
      if (s != NULL)
      {
        sprintf(val, "%f", ((t_ras *)(parent->ds->obj))->alpha->v);
        XML_Set_Attribute_Value(child, "alpha", val);
      }
    }
    else if (!strcmp(s, "gamma+inv"))
    {
      s = XML_Get_Attribute_Value(child, "alpha");
      if (s != NULL)
      {
        sprintf(val, "%f", ((t_ras *)(parent->ds->obj))->alpha->v);
        XML_Set_Attribute_Value(child, "alpha", val);

        sprintf(val, "%f", ((t_ras *)(parent->ds->obj))->pinvar->v);
        XML_Set_Attribute_Value(child, "pinv", val);
      }
    }

    child = XML_Search_Node_Name("weights", YES, parent);
    child = child->child;
    assert(child);

    class = 0;
    do
    {
      // Find the rate class that this weight instance points to
      s = XML_Get_Attribute_Value(child, "appliesto");
      assert(s);
      x = XML_Search_Node_Attribute_Value("id", s, YES, parent);
      assert(x);

      s = XML_Get_Attribute_Value(x, "init.value");
      if (s != NULL)
      {
        v = atof(s);
        if (Are_Equal(v, 0.0, 1E-20) == NO)
        {
          sprintf(val, "%f",
                  ((t_ras *)(parent->ds->obj))->gamma_r_proba->v[class]);
          XML_Set_Attribute_Value(child, "value", val);
          class ++;
        }
        else
        {
          sprintf(val, "%f", 0.0);
          XML_Set_Attribute_Value(child, "value", val);
        }
      }
      child = child->next;
    } while (child);
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
