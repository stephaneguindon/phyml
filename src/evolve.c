/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include "evolve.h"


int EVOLVE_Main(int argc, char **argv)
{
  option *io;
  t_tree *tree;
  calign *cdata;
  char   *dum;
  FILE   *fp;
  int     r_seed;

  dum = (char *)mCalloc(100, sizeof(char));

  io  = NULL;

  io = (option *)Get_Input(argc, argv);
  if (!io)
    return (0);
  else if (io->use_xml == YES)
  {
    Free(io);
    return (0);
  }

  r_seed = (io->r_seed < 0) ? (time(NULL)) : (io->r_seed);
  srand(r_seed);
  io->r_seed = r_seed;

  io->data = Make_Empty_Alignment(io);

  Make_Model_Complete(io->mod);
  Set_Model_Name(io->mod);
  Print_Settings(io);

  io->colalias = NO;
  cdata        = Compact_Data(io->data, io);
  Free_Seq(io->data, io->n_otu);

  tree = Make_Tree_From_Scratch(io->n_otu, cdata);

  Connect_CSeqs_To_Nodes(cdata, io, tree);

  tree->rates = RATES_Make_Rate_Struct(tree->n_otu);
  RATES_Init_Rate_Struct(tree->rates, io->rates, tree->n_otu);

  tree->times = TIMES_Make_Time_Struct(tree->n_otu);
  TIMES_Init_Time_Struct(tree->times, io->times, tree->n_otu);

  tree->data                   = cdata;
  tree->mod                    = io->mod;
  tree->io                     = io;
  tree->times->scaled_pop_size = 1.E+3;
  tree->times->neff_growth     = 0.5;

  EVOLVE_Coalescent(tree);

  Init_Model(tree->data, tree->mod, io);
  Set_Model_Parameters(tree->mod);
  Set_Model_Name(tree->mod);

  Make_Tree_For_Pars(tree);
  Make_Tree_For_Lk(tree);
  Make_Spr(tree);

  sprintf(dum, "%s%d%s", "./", io->r_seed, "_evolve_stats.txt");
  fp = Openfile(dum, WRITE);
  EVOLVE_Seq(tree->data, tree->mod, fp, tree);
  fclose(fp);
  
  sprintf(dum, "%s%d%s", "./", io->r_seed, "_evolve_data.txt");
  fp = Openfile(dum, WRITE);
  Print_CSeq(fp, NO, tree->data, tree);
  fclose(fp);

  sprintf(dum, "%s%d%s", "./", io->r_seed, "_evolve_tree.txt");
  fp = Openfile(dum, WRITE);
  PhyML_Fprintf(fp, "%s\n", Write_Tree(tree));
  fclose(fp);

      

  // Write XML configuration file for downstream PhyML analyses

  xml_node *root, *nd, *ndnd, *ndndnd;
  char     *model, *cv_type;

  model = (char *)mCalloc(100, sizeof(char));
  cv_type = (char *)mCalloc(100, sizeof(char));

  strcpy(model, "GTR");
  strcpy(cv_type, "kfold.pos");

  root = XML_Make_Node("phyml");
  XML_Init_Node(NULL, root, "phyml");
  sprintf(dum, "%s_%s",model, cv_type);
  root->attr = XML_Make_Attribute(NULL, "run.id",dum);
  sprintf(dum, "%d_%s_%s", io->r_seed, model, cv_type);
  XML_Add_Attribute(root, "output.file",dum);
  XML_Add_Attribute(root, "cv.type",cv_type);

  nd         = XML_Add_Node(root, "topology");
  ndnd       = XML_Add_Node(nd, "instance");
  ndnd->attr = XML_Make_Attribute(NULL, "id", "T1");
  XML_Add_Attribute(ndnd, "init.tree", "BioNJ");

  nd         = XML_Add_Node(root, "ratematrices");
  nd->attr   = XML_Make_Attribute(NULL, "id", "RM1");
  ndnd       = XML_Add_Node(nd, "instance");
  ndnd->attr = XML_Make_Attribute(NULL, "id", "M1");
  XML_Add_Attribute(ndnd, "model", model);
  XML_Add_Attribute(ndnd, "optimise.rr","yes");

  nd         = XML_Add_Node(root, "siterates");
  nd->attr   = XML_Make_Attribute(NULL, "id", "SR1");
  ndnd       = XML_Add_Node(nd, "instance");
  ndnd->attr = XML_Make_Attribute(NULL, "id", "R1");
  XML_Add_Attribute(ndnd, "init.value", "1.0");
  // ndnd       = XML_Add_Node(nd, "instance");
  // ndnd->attr = XML_Make_Attribute(NULL, "id", "R2");
  // XML_Add_Attribute(ndnd, "init.value", "1.0");
  // ndnd       = XML_Add_Node(nd, "instance");
  // ndnd->attr = XML_Make_Attribute(NULL, "id", "R3");
  // XML_Add_Attribute(ndnd, "init.value", "1.0");
  ndnd       = XML_Add_Node(nd, "weights");
  ndnd->attr = XML_Make_Attribute(NULL, "id", "D1");
  XML_Add_Attribute(ndnd, "family", "freerates");
  ndndnd       = XML_Add_Node(ndnd, "instance");
  ndndnd->attr = XML_Make_Attribute(NULL, "appliesto", "R1");
  XML_Add_Attribute(ndndnd, "value", "1.0");
  // ndndnd       = XML_Add_Node(ndnd, "instance");
  // ndndnd->attr = XML_Make_Attribute(NULL, "appliesto", "R2");
  // XML_Add_Attribute(ndndnd, "value", "1.0");
  // ndndnd       = XML_Add_Node(ndnd, "instance");
  // ndndnd->attr = XML_Make_Attribute(NULL, "appliesto", "R3");
  // XML_Add_Attribute(ndndnd, "value", "1.0");

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
  dum      = (char *)mCalloc(100, sizeof(char));
  sprintf(dum, "%s%d%s", "./", io->r_seed, "_evolve_data.txt");
  XML_Add_Attribute(nd, "file.name", dum);
  XML_Add_Attribute(nd, "data.type", "nt");
  XML_Add_Attribute(nd, "interleaved", "no");

  // ndnd       = XML_Add_Node(nd, "mixtureelem");
  // ndnd->attr = XML_Make_Attribute(NULL, "list", "T1,T1,T1");
  // ndnd       = XML_Add_Node(nd, "mixtureelem");
  // ndnd->attr = XML_Make_Attribute(NULL, "list", "M1,M1,M1");
  // ndnd       = XML_Add_Node(nd, "mixtureelem");
  // ndnd->attr = XML_Make_Attribute(NULL, "list", "F1,F1,F1");
  // ndnd       = XML_Add_Node(nd, "mixtureelem");
  // ndnd->attr = XML_Make_Attribute(NULL, "list", "R1,R2,R3");
  // ndnd       = XML_Add_Node(nd, "mixtureelem");
  // ndnd->attr = XML_Make_Attribute(NULL, "list", "L1,L1,L1");

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

  sprintf(dum, "%s%d_%s_%s_%s", "./",io->r_seed, model, cv_type, "evolve_config.xml");
  fp = Openfile(dum, WRITE);
  XML_Write_XML_Graph(fp, root);
  fclose(fp);


///
///

  strcpy(model, "GTR");
  strcpy(cv_type, "kfold.col");

  sprintf(dum, "%s_%s", model, cv_type);
  XML_Set_Attribute_Value(root, "run.id", dum);
  sprintf(dum, "%d_%s_%s", io->r_seed, model, cv_type);
  XML_Set_Attribute_Value(root, "output.file", dum);
  XML_Set_Attribute_Value(root, "cv.type", cv_type);

  nd         = XML_Search_Node_Name("ratematrices", NO, root);
  ndnd       = XML_Search_Node_Name("instance" , NO, nd);
  XML_Set_Attribute_Value(ndnd, "model", model);

  sprintf(dum, "%s%d_%s_%s_%s", "./",io->r_seed, model, cv_type,
          "evolve_config.xml");
  fp = Openfile(dum, WRITE);
  XML_Write_XML_Graph(fp, root);
  fclose(fp);

  ///
  ///

  strcpy(model, "GTR");
  strcpy(cv_type, "maxfold");

  sprintf(dum, "%s_%s", model, cv_type);
  XML_Set_Attribute_Value(root, "run.id", dum);
  sprintf(dum, "%d_%s_%s", io->r_seed, model, cv_type);
  XML_Set_Attribute_Value(root, "output.file", dum);
  XML_Set_Attribute_Value(root, "cv.type", cv_type);

  nd   = XML_Search_Node_Name("ratematrices", NO, root);
  ndnd = XML_Search_Node_Name("instance", NO, nd);
  XML_Set_Attribute_Value(ndnd, "model", model);

  sprintf(dum, "%s%d_%s_%s_%s", "./", io->r_seed, model, cv_type,
          "evolve_config.xml");
  fp = Openfile(dum, WRITE);
  XML_Write_XML_Graph(fp, root);
  fclose(fp);

  ///
  ///

  strcpy(model, "HKY85");
  strcpy(cv_type, "kfold.pos");

  sprintf(dum, "%s_%s", model, cv_type);
  XML_Set_Attribute_Value(root, "run.id", dum);
  sprintf(dum, "%d_%s_%s", io->r_seed, model, cv_type);
  XML_Set_Attribute_Value(root, "output.file", dum);
  XML_Set_Attribute_Value(root, "cv.type", cv_type);

  nd   = XML_Search_Node_Name("ratematrices", NO, root);
  ndnd = XML_Search_Node_Name("instance", NO, nd);
  XML_Set_Attribute_Value(ndnd, "model", model);

  sprintf(dum, "%s%d_%s_%s_%s", "./", io->r_seed, model, cv_type,
          "evolve_config.xml");
  fp = Openfile(dum, WRITE);
  XML_Write_XML_Graph(fp, root);
  fclose(fp);

  ///
  ///

  strcpy(model, "HKY85");
  strcpy(cv_type, "kfold.col");

  sprintf(dum, "%s_%s", model, cv_type);
  XML_Set_Attribute_Value(root, "run.id", dum);
  sprintf(dum, "%d_%s_%s", io->r_seed, model, cv_type);
  XML_Set_Attribute_Value(root, "output.file", dum);
  XML_Set_Attribute_Value(root, "cv.type", cv_type);

  nd   = XML_Search_Node_Name("ratematrices", NO, root);
  ndnd = XML_Search_Node_Name("instance", NO, nd);
  XML_Set_Attribute_Value(ndnd, "model", model);

  sprintf(dum, "%s%d_%s_%s_%s", "./", io->r_seed, model, cv_type,
          "evolve_config.xml");
  fp = Openfile(dum, WRITE);
  XML_Write_XML_Graph(fp, root);
  fclose(fp);

  ///
  ///

  strcpy(model, "HKY85");
  strcpy(cv_type, "maxfold");

  sprintf(dum, "%s_%s", model, cv_type);
  XML_Set_Attribute_Value(root, "run.id", dum);
  sprintf(dum, "%d_%s_%s", io->r_seed, model, cv_type);
  XML_Set_Attribute_Value(root, "output.file", dum);
  XML_Set_Attribute_Value(root, "cv.type", cv_type);

  nd   = XML_Search_Node_Name("ratematrices", NO, root);
  ndnd = XML_Search_Node_Name("instance", NO, nd);
  XML_Set_Attribute_Value(ndnd, "model", model);

  sprintf(dum, "%s%d_%s_%s_%s", "./", io->r_seed, model, cv_type,
          "evolve_config.xml");
  fp = Openfile(dum, WRITE);
  XML_Write_XML_Graph(fp, root);
  fclose(fp);

  ///
  ///

  strcpy(model, "K80");
  strcpy(cv_type, "kfold.pos");

  sprintf(dum, "%s_%s", model, cv_type);
  XML_Set_Attribute_Value(root, "run.id", dum);
  sprintf(dum, "%d_%s_%s", io->r_seed, model, cv_type);
  XML_Set_Attribute_Value(root, "output.file", dum);
  XML_Set_Attribute_Value(root, "cv.type", cv_type);

  nd   = XML_Search_Node_Name("ratematrices", NO, root);
  ndnd = XML_Search_Node_Name("instance", NO, nd);
  XML_Set_Attribute_Value(ndnd, "model", model);

  sprintf(dum, "%s%d_%s_%s_%s", "./", io->r_seed, model, cv_type,
          "evolve_config.xml");
  fp = Openfile(dum, WRITE);
  XML_Write_XML_Graph(fp, root);
  fclose(fp);

  ///
  ///

  strcpy(model, "K80");
  strcpy(cv_type, "kfold.col");

  sprintf(dum, "%s_%s", model, cv_type);
  XML_Set_Attribute_Value(root, "run.id", dum);
  sprintf(dum, "%d_%s_%s", io->r_seed, model, cv_type);
  XML_Set_Attribute_Value(root, "output.file", dum);
  XML_Set_Attribute_Value(root, "cv.type", cv_type);

  nd   = XML_Search_Node_Name("ratematrices", NO, root);
  ndnd = XML_Search_Node_Name("instance", NO, nd);
  XML_Set_Attribute_Value(ndnd, "model", model);

  sprintf(dum, "%s%d_%s_%s_%s", "./", io->r_seed, model, cv_type,
          "evolve_config.xml");
  fp = Openfile(dum, WRITE);
  XML_Write_XML_Graph(fp, root);
  fclose(fp);

  ///
  ///

  strcpy(model, "K80");
  strcpy(cv_type, "maxfold");

  sprintf(dum, "%s_%s", model, cv_type);
  XML_Set_Attribute_Value(root, "run.id", dum);
  sprintf(dum, "%d_%s_%s", io->r_seed, model, cv_type);
  XML_Set_Attribute_Value(root, "output.file", dum);
  XML_Set_Attribute_Value(root, "cv.type", cv_type);

  nd   = XML_Search_Node_Name("ratematrices", NO, root);
  ndnd = XML_Search_Node_Name("instance", NO, nd);
  XML_Set_Attribute_Value(ndnd, "model", model);

  sprintf(dum, "%s%d_%s_%s_%s", "./", io->r_seed, model, cv_type,
          "evolve_config.xml");
  fp = Openfile(dum, WRITE);
  XML_Write_XML_Graph(fp, root);
  fclose(fp);

  ///
  ///

  strcpy(model, "JC69");
  strcpy(cv_type, "kfold.pos");

  sprintf(dum, "%s_%s", model, cv_type);
  XML_Set_Attribute_Value(root, "run.id", dum);
  sprintf(dum, "%d_%s_%s", io->r_seed, model, cv_type);
  XML_Set_Attribute_Value(root, "output.file", dum);
  XML_Set_Attribute_Value(root, "cv.type", cv_type);

  nd   = XML_Search_Node_Name("ratematrices", NO, root);
  ndnd = XML_Search_Node_Name("instance", NO, nd);
  XML_Set_Attribute_Value(ndnd, "model", model);

  sprintf(dum, "%s%d_%s_%s_%s", "./", io->r_seed, model, cv_type,
          "evolve_config.xml");
  fp = Openfile(dum, WRITE);
  XML_Write_XML_Graph(fp, root);
  fclose(fp);

  ///
  ///

  strcpy(model, "JC69");
  strcpy(cv_type, "kfold.col");

  sprintf(dum, "%s_%s", model, cv_type);
  XML_Set_Attribute_Value(root, "run.id", dum);
  sprintf(dum, "%d_%s_%s", io->r_seed, model, cv_type);
  XML_Set_Attribute_Value(root, "output.file", dum);
  XML_Set_Attribute_Value(root, "cv.type", cv_type);

  nd   = XML_Search_Node_Name("ratematrices", NO, root);
  ndnd = XML_Search_Node_Name("instance", NO, nd);
  XML_Set_Attribute_Value(ndnd, "model", model);

  sprintf(dum, "%s%d_%s_%s_%s", "./", io->r_seed, model, cv_type,
          "evolve_config.xml");
  fp = Openfile(dum, WRITE);
  XML_Write_XML_Graph(fp, root);
  fclose(fp);

  ///
  ///

  strcpy(model, "JC69");
  strcpy(cv_type, "maxfold");

  sprintf(dum, "%s_%s", model, cv_type);
  XML_Set_Attribute_Value(root, "run.id", dum);
  sprintf(dum, "%d_%s_%s", io->r_seed, model, cv_type);
  XML_Set_Attribute_Value(root, "output.file", dum);
  XML_Set_Attribute_Value(root, "cv.type", cv_type);

  nd   = XML_Search_Node_Name("ratematrices", NO, root);
  ndnd = XML_Search_Node_Name("instance", NO, nd);
  XML_Set_Attribute_Value(ndnd, "model", model);

  sprintf(dum, "%s%d_%s_%s_%s", "./", io->r_seed, model, cv_type,
          "evolve_config.xml");
  fp = Openfile(dum, WRITE);
  XML_Write_XML_Graph(fp, root);
  fclose(fp);


  Free(dum);

  return (-1);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void EVOLVE_Coalescent(t_tree *tree)
{
  t_node   *n, *new_n, **avail;
  int       n_lineages, idx_ancestor, *rand_idx, idx_edge;
  short int no_tip_left;
  phydbl    Ne, dt;
  phydbl    T;

  Ne = tree->times->scaled_pop_size;

  avail = (t_node **)mCalloc(tree->n_otu, sizeof(t_node *));

  if (Ne > tree->times->scaled_pop_size_max ||
      Ne < tree->times->scaled_pop_size_min)
    assert(FALSE);

  Get_Node_Ranks_From_Tip_Times(tree);

  n = tree->a_nodes[0];
  while (n->rk_next) n = n->rk_next;

  for (int i = 0; i < tree->n_otu; ++i) avail[i] = NULL;
  avail[0] = n;
  avail[1] = n->rk_prev;

  // PhyML_Printf("\n. n->time = %f n->rk_prev->time: %f",
  // tree->times->nd_t[n->num], tree->times->nd_t[n->rk_prev->num]);
  idx_ancestor = tree->n_otu;
  idx_edge     = 0;
  new_n        = NULL;
  n            = n->rk_prev;
  n_lineages   = 2;
  no_tip_left  = FALSE;
  T            = tree->times->nd_t[n->num];
  do
  {
    dt = Rexp(Bico(n_lineages, 2) / Ne);
    T  = T - dt;

    // PhyML_Printf("\n. T: %f", T);

    if (n->rk_prev == NULL) no_tip_left = TRUE;

    if (n->rk_prev != NULL && T < tree->times->nd_t[n->rk_prev->num])
    {
      T = tree->times->nd_t[n->rk_prev->num];
      if (n->rk_prev->tax == YES) n = n->rk_prev;
      avail[n_lineages] = n;
      n_lineages++;
      // PhyML_Printf("\n. Found tip %s at time %f", n->name,
      // tree->times->nd_t[n->num]);
    }
    else
    {
      new_n                           = tree->a_nodes[idx_ancestor];
      tree->times->nd_t[idx_ancestor] = T;

      idx_ancestor++;

      rand_idx = Permutate(n_lineages);

      avail[rand_idx[0]]->v[0] = new_n;
      avail[rand_idx[1]]->v[0] = new_n;

      new_n->v[1] = avail[rand_idx[0]];
      new_n->v[2] = avail[rand_idx[1]];

      avail[rand_idx[0]]    = new_n;
      avail[rand_idx[1]]    = avail[n_lineages - 1];
      avail[n_lineages - 1] = NULL;

      Connect_One_Edge_To_Two_Nodes(new_n, new_n->v[1], tree->a_edges[idx_edge],
                                    tree);
      Connect_One_Edge_To_Two_Nodes(new_n, new_n->v[2],
                                    tree->a_edges[idx_edge + 1], tree);

      // PhyML_Printf("\n. Added coalescent node %d at time %f descendants %d %d
      // # lineages: %d",
      //              idx_ancestor - 1, tree->times->nd_t[idx_ancestor - 1],
      //              new_n->v[1]->num,
      //              new_n->v[2]->num,
      //              n_lineages);

      n_lineages--;
      idx_edge += 2;
      Free(rand_idx);
      // for(int i=0; i < n_lineages; ++i) PhyML_Printf("\n. Avail:
      // %d",avail[i]->num);
    }
  } while (n_lineages > 1 || no_tip_left == FALSE);

  tree->n_root = new_n;

  tree->n_root->v[1]->v[0] = tree->n_root->v[2];
  tree->n_root->v[2]->v[0] = tree->n_root->v[1];
  Connect_One_Edge_To_Two_Nodes(tree->n_root->v[1], tree->n_root->v[2],
                                tree->a_edges[idx_edge], tree);
  tree->e_root = tree->a_edges[idx_edge];

  // Transform of times for exponentially growning or declining pop size
  phydbl g = tree->times->neff_growth;
  for (int i = 0; i < 2 * tree->n_otu - 1; ++i)
  {
    tree->times->nd_t[i] = -((1. / g) * log(1. + g * -tree->times->nd_t[i]));
    assert(!isnan(tree->times->nd_t[i]));
  }

  Update_Ancestors(tree->n_root, tree->n_root->v[2], tree->n_root->b[2], tree);
  Update_Ancestors(tree->n_root, tree->n_root->v[1], tree->n_root->b[1], tree);
  RATES_Fill_Lca_Table(tree);

  phydbl L                = TIMES_Tree_Length(tree);
  tree->rates->bl_from_rt = YES;
  tree->rates->clock_r    = 0.1 / L * (2 * tree->n_otu - 2);
  tree->rates->model_id   = STRICTCLOCK;

  RATES_Update_Edge_Lengths(tree);

  Free(avail);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void EVOLVE_Seq(calign *data, t_mod *mod, FILE *fp_stats, t_tree *tree)
{
  int        root_state, root_rate_class;
  int        site, n_otu, ns;
  phydbl    *orig_l,*weights;
  phydbl     shape, scale, var, mean, r_mult, sum;
  phydbl    *state_probs_one_site, *state_probs_all_sites;
  int        switch_to_yes;
  short int *truth;
  
  orig_l = (phydbl *)mCalloc(2 * tree->n_otu - 1, sizeof(phydbl));
  for (int i = 0; i < 2 * tree->n_otu - 1; ++i)
    orig_l[i] = tree->a_edges[i]->l->v;

  data->n_otu  = tree->n_otu;
  data->io     = tree->io;
  n_otu        = tree->n_otu;
  ns           = tree->mod->ns;

  state_probs_all_sites =
      (phydbl *)mCalloc(ns * n_otu * data->init_len, sizeof(phydbl));
  truth = (short int *)mCalloc(ns * n_otu * data->init_len, sizeof(short int));
  weights = (phydbl *)mCalloc(n_otu * data->init_len, sizeof(phydbl));

  for (int i = 0; i < n_otu * data->init_len; ++i) weights[i] = 1.;

  if (mod->use_m4mod) tree->print_labels = YES;

  Set_Br_Len_Var(NULL, tree);

  switch_to_yes = NO;
  if (tree->mod->gamma_mgf_bl == YES) switch_to_yes = YES;

  Set_Update_Eigen(YES, mod);

  if (!Set_Model_Parameters(mod))
    Generic_Exit(__FILE__, __LINE__, __FUNCTION__);

  if (tree->mod->whichmodel == GTR)
  {
    for (int i = 0; i < mod->ns; ++i)
      tree->mod->e_frq->pi_unscaled->v[i] = Uni();

    sum = 0.0;
    for (int i = 0; i < mod->ns; ++i)
      sum += tree->mod->e_frq->pi_unscaled->v[i];

    for (int i = 0; i < mod->ns; ++i)
      tree->mod->e_frq->pi->v[i] = tree->mod->e_frq->pi_unscaled->v[i] / sum;

    tree->mod->r_mat->rr_val->v[AC] = log(Uni() * (2. - .5) + .5);
    tree->mod->r_mat->rr_val->v[AG] = log(Uni() * (8. - 2.) + 2.);
    tree->mod->r_mat->rr_val->v[AT] = log(Uni() * (2. - .5) + .5);
    tree->mod->r_mat->rr_val->v[CG] = log(Uni() * (2. - .5) + .5);
    tree->mod->r_mat->rr_val->v[CT] = log(Uni() * (8. - 2.) + 2.);
    tree->mod->r_mat->rr_val->v[GT] = log(Uni() * (2. - .5) + .5);

    tree->mod->r_mat->n_diff_rr = 6;
    tree->mod->whichmodel       = GTR;
    tree->mod->update_eigen     = YES;
    Update_Eigen(tree->mod);
  }
  else if (tree->mod->whichmodel == HKY85)
  {
    for (int i = 0; i < mod->ns; ++i)
      tree->mod->e_frq->pi_unscaled->v[i] = Uni();

    sum = 0.0;
    for (int i = 0; i < mod->ns; ++i)
      sum += tree->mod->e_frq->pi_unscaled->v[i];

    for (int i = 0; i < mod->ns; ++i)
      tree->mod->e_frq->pi->v[i] = tree->mod->e_frq->pi_unscaled->v[i] / sum;

    tree->mod->kappa->v = Uni() * (8. - 2.) + 2.;

    tree->mod->custom_mod_string->s[0] = '0';
    tree->mod->custom_mod_string->s[1] = '1';
    tree->mod->custom_mod_string->s[2] = '2';
    tree->mod->custom_mod_string->s[3] = '3';
    tree->mod->custom_mod_string->s[4] = '4';
    tree->mod->custom_mod_string->s[5] = '5';
    Translate_Custom_Mod_String(tree->mod);

    tree->mod->r_mat->rr_val->v[AC] = log(1.0);
    tree->mod->r_mat->rr_val->v[AG] = log(tree->mod->kappa->v);
    tree->mod->r_mat->rr_val->v[AT] = log(1.0);
    tree->mod->r_mat->rr_val->v[CG] = log(1.0);
    tree->mod->r_mat->rr_val->v[CT] = log(tree->mod->kappa->v);
    tree->mod->r_mat->rr_val->v[GT] = log(1.0);

    tree->mod->r_mat->n_diff_rr = 6;
    tree->mod->whichmodel       = GTR;
    tree->mod->update_eigen     = YES;
    Update_Eigen(tree->mod);
  }
  else if (tree->mod->whichmodel == K80)
  {
    for (int i = 0; i < mod->ns; ++i)
      tree->mod->e_frq->pi_unscaled->v[i] = 1.0;

    sum = 0.0;
    for (int i = 0; i < mod->ns; ++i)
      sum += tree->mod->e_frq->pi_unscaled->v[i];

    for (int i = 0; i < mod->ns; ++i)
      tree->mod->e_frq->pi->v[i] = tree->mod->e_frq->pi_unscaled->v[i] / sum;


    tree->mod->kappa->v = Uni() * (8. - 2.) + 2.;

    tree->mod->custom_mod_string->s[0] = '0';
    tree->mod->custom_mod_string->s[1] = '1';
    tree->mod->custom_mod_string->s[2] = '2';
    tree->mod->custom_mod_string->s[3] = '3';
    tree->mod->custom_mod_string->s[4] = '4';
    tree->mod->custom_mod_string->s[5] = '5';
    Translate_Custom_Mod_String(tree->mod);

    tree->mod->r_mat->rr_val->v[AC] = log(1.0);
    tree->mod->r_mat->rr_val->v[AG] = log(tree->mod->kappa->v);
    tree->mod->r_mat->rr_val->v[AT] = log(1.0);
    tree->mod->r_mat->rr_val->v[CG] = log(1.0);
    tree->mod->r_mat->rr_val->v[CT] = log(tree->mod->kappa->v);
    tree->mod->r_mat->rr_val->v[GT] = log(1.0);

    tree->mod->r_mat->n_diff_rr = 6;
    tree->mod->whichmodel       = GTR;
    tree->mod->update_eigen     = YES;
    Update_Eigen(tree->mod);
  }
  else if (tree->mod->whichmodel == JC69)
  {
    for (int i = 0; i < mod->ns; ++i)
      tree->mod->e_frq->pi_unscaled->v[i] = 1.0;

    sum = 0.0;
    for (int i = 0; i < mod->ns; ++i)
      sum += tree->mod->e_frq->pi_unscaled->v[i];

    for (int i = 0; i < mod->ns; ++i)
      tree->mod->e_frq->pi->v[i] = tree->mod->e_frq->pi_unscaled->v[i] / sum;

    tree->mod->custom_mod_string->s[0] = '0';
    tree->mod->custom_mod_string->s[1] = '1';
    tree->mod->custom_mod_string->s[2] = '2';
    tree->mod->custom_mod_string->s[3] = '3';
    tree->mod->custom_mod_string->s[4] = '4';
    tree->mod->custom_mod_string->s[5] = '5';
    Translate_Custom_Mod_String(tree->mod);

    tree->mod->r_mat->rr_val->v[AC] = log(1.0);
    tree->mod->r_mat->rr_val->v[AG] = log(1.0);
    tree->mod->r_mat->rr_val->v[AT] = log(1.0);
    tree->mod->r_mat->rr_val->v[CG] = log(1.0);
    tree->mod->r_mat->rr_val->v[CT] = log(1.0);
    tree->mod->r_mat->rr_val->v[GT] = log(1.0);

    tree->mod->r_mat->n_diff_rr = 6;
    tree->mod->whichmodel       = GTR;
    tree->mod->update_eigen     = YES;
    Update_Eigen(tree->mod);
  }

  PhyML_Fprintf(fp_stats ? fp_stats : stdout,
                "\n. pi: %7f %7f %7f %7f rr: %7f %7f %7f %7f %7f %7f",
                mod->e_frq->pi->v[0], mod->e_frq->pi->v[1],
                mod->e_frq->pi->v[2], mod->e_frq->pi->v[3],
                mod->r_mat->rr->v[AC], mod->r_mat->rr->v[AG],
                mod->r_mat->rr->v[AT], mod->r_mat->rr->v[CG],
                mod->r_mat->rr->v[CT], mod->r_mat->rr->v[GT]);

  for (site = 0; site < data->init_len; ++site)
  {

    /* Pick the rate class */
    root_state = root_rate_class = 0;

    // root_rate_class =
    //     EVOLVE_Pick_State(mod->ras->n_catg, mod->ras->gamma_r_proba->v);

    var  = 5.0;
    mean = 1.0;

    // Set edge lengths
    shape = mean * mean / var;
    scale = var / mean;

    // r_mult = Rgamma(shape, scale);
    r_mult = 1.0; // No scaling of edge lengths -> no RAS

    // RAS : all edges multiplied by gamma distributed scaling factor
    for (int i = 0; i < 2 * tree->n_otu - 1; ++i)
      tree->a_edges[i]->l->v = orig_l[i] * r_mult;


    for (int i = 0; i < 2 * tree->n_otu - 1; ++i)
      Update_PMat_At_Given_Edge(tree->a_edges[i], tree);

    // Pick the root nucleotide/aa
    root_state = EVOLVE_Pick_State(mod->ns, mod->e_frq->pi->v);
    tree->a_nodes[0]->c_seq->state[site] =
        Reciproc_Assign_State(root_state, tree->io->datatype);
    tree->a_nodes[0]->c_seq->d_state[site] = root_state;


    // PhyML_Printf("\n. root_state: %d root_rate_class: %d [%f %f %f %f]",
    //              root_state,
    //              root_rate_class,
    //              mod->e_frq->pi->v[0],
    //              mod->e_frq->pi->v[1],
    //              mod->e_frq->pi->v[2],
    //              mod->e_frq->pi->v[3]);
    /* Generic_Exit(__FILE__,__LINE__,__FUNCTION__); */

    /* tree->a_nodes[0] is considered as the root t_node */
    EVOLVE_Seq_Recur(tree->a_nodes[0], tree->a_nodes[0]->v[0],
                     tree->a_nodes[0]->b[0], root_state, root_rate_class, site,
                     data, mod, tree);

    data->wght[site] = 1;

    state_probs_one_site = EVOLVE_Site_Lk(site, root_rate_class, data, tree);

    for (int tax_id = 0; tax_id < n_otu; ++tax_id)
    {
      for (int state = 0; state < ns; ++state)
      {
        state_probs_all_sites[site * n_otu * ns + tax_id * ns + state] =
            state_probs_one_site[tax_id * ns + state];
        truth[site * n_otu * ns + tax_id * ns + state] = 0;                  
      }
      truth[site * n_otu * ns + tax_id * ns +
            tree->a_nodes[tax_id]->c_seq->d_state[site]] = 1;
    }

    Free(state_probs_one_site);
  }

  data->n_pattern = data->init_len;
  /* Print_CSeq(stdout,NO,data); */
  for (int i = 0; i < 2 * tree->n_otu - 3; ++i)
    tree->a_edges[i]->l->v = orig_l[i];
  Free(orig_l);

  ROC(state_probs_all_sites, truth, ns, n_otu * data->init_len, weights,
      "SIM",fp_stats ? fp_stats : stdout);

  Free(state_probs_all_sites);
  Free(truth);
  Free(weights);

  if (switch_to_yes == YES) tree->mod->gamma_mgf_bl = YES;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int EVOLVE_Pick_State(int n, phydbl *prob)
{
  int    pos;
  phydbl uni;

  do
  {
    pos = rand();
    pos = (pos % n);
    uni = (phydbl)rand();
    uni /= (phydbl)RAND_MAX;
    if (uni < prob[pos]) break;
  } while (1);

  return (int)pos;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void EVOLVE_Seq_Recur(t_node *a, t_node *d, t_edge *b, int a_state, int r_class,
                      int site_num, calign *gen_data, t_mod *mod, t_tree *tree)
{
  int d_state;
  int dim1, dim2;

  dim1 = tree->mod->ns * tree->mod->ns;
  dim2 = tree->mod->ns;

  //   PhyML_Printf("\n## L:%G %G %G %G %G",
  //                b->l->v,
  //                b->Pij_rr[r_class * dim1 + a_state * dim2 + 0],
  //                b->Pij_rr[r_class * dim1 + a_state * dim2 + 1],
  //                b->Pij_rr[r_class * dim1 + a_state * dim2 + 2],
  //                b->Pij_rr[r_class * dim1 + a_state * dim2 + 3]);

  d_state =
      EVOLVE_Pick_State(mod->ns, b->Pij_rr + r_class * dim1 + a_state * dim2);

  //   PhyML_Printf("\n>> %c->%c L:%G %G %G %G %G",
  //                Reciproc_Assign_State(a_state, mod->io->datatype),
  //                Reciproc_Assign_State(d_state, mod->io->datatype),
  //                b->l->v,
  //                b->Pij_rr[r_class * dim1 + a_state * dim2 + 0],
  //                b->Pij_rr[r_class * dim1 + a_state * dim2 + 1],
  //                b->Pij_rr[r_class * dim1 + a_state * dim2 + 2],
  //                b->Pij_rr[r_class * dim1 + a_state * dim2 + 3]);

  if (d->tax)
  {
    d->c_seq->state[site_num] =
        Reciproc_Assign_State(d_state, tree->io->datatype);
    d->c_seq->d_state[site_num] = d_state;
    return;
  }
  else
  {
    int i;
    for (i = 0; i < 3; i++)
      if (d->v[i] != a && !(a == tree->n_root && d->b[i] == tree->e_root))
      {
        EVOLVE_Seq_Recur(d, d->v[i], d->b[i], d_state, r_class, site_num,
                         gen_data, mod, tree);
      }
  }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// Evaluate likelihood at site_idx of data under the very same model
// that was used to generate that site, i.e., the same set of
// transition probability matrices.
phydbl *EVOLVE_Site_Lk(int site_idx, int rate_class, calign *data, t_tree *tree)
{
  int     ori_len, ns, n_otu;
  align **dum_data;
  calign *ori_data;
  t_edge *b;
  phydbl *Pij, *p_lk_left, sum, *state_probs;

  ori_len            = tree->io->init_len;
  tree->io->init_len = 1;
  ns                 = tree->mod->ns;
  ori_data           = tree->data;
  n_otu              = tree->n_otu;

  state_probs = (phydbl *)mCalloc(ns * n_otu, sizeof(phydbl));

  dum_data = Make_Empty_Alignment(tree->io);
  for (int i = 0; i < tree->n_otu; ++i)
  {
    dum_data[i]->state[0] = data->c_seq[i]->state[site_idx];
    strcpy(dum_data[i]->name, data->c_seq[i]->name);
  }

  tree->data = Compact_Data(dum_data, tree->io);
  Free_Seq(dum_data, tree->n_otu);

  Connect_CSeqs_To_Nodes(tree->data, tree->io, tree);

  for (int i = 0; i < tree->n_otu; ++i)
    Init_Partial_Lk_Tips_Double_One_Character(i, 0, tree);

  Post_Order_Lk(tree->n_root, tree->n_root->v[1], tree);
  Post_Order_Lk(tree->n_root, tree->n_root->v[2], tree);

  Pre_Order_Lk(tree->n_root, tree->n_root->v[2], tree);
  Pre_Order_Lk(tree->n_root, tree->n_root->v[1], tree);

  for (int tax_id = 0; tax_id < tree->n_otu; ++tax_id)
  {
    b = tree->a_nodes[tax_id]->b[0];

    p_lk_left = b->p_lk_left;
    Pij       = b->Pij_rr;

    for (int tip_state = 0; tip_state < ns; ++tip_state)
    {
      state_probs[tax_id * ns + tip_state] = 0.0;
      for (int int_state = 0; int_state < ns; ++int_state)
      {
        state_probs[tax_id * ns + tip_state] +=
            tree->mod->ras->gamma_r_proba->v[rate_class] *
            tree->mod->e_frq->pi->v[tip_state] *
            Pij[rate_class * ns * ns + tip_state * ns + int_state] *
            p_lk_left[rate_class * ns + int_state];
      }
    }

    sum = 0.0;
    for (int state = 0; state < ns; ++state)
      sum += state_probs[tax_id * ns + state];
    for (int state = 0; state < ns; ++state)
      state_probs[tax_id * ns + state] /= sum;

    // for (int state = 0; state < ns; ++state)
    //     PhyML_Printf("\n. Site: %4d tax: %20s state: %3d prob: %15G
    //     obs_state: %3d",
    //                  site_idx,
    //                  data->c_seq[tax_id]->name,
    //                  state,
    //                  state_prob[tax_id * ns + state],
    //                  tree->a_nodes[tax_id]->c_seq->d_state[0]);

    // PhyML_Printf("\n###%s,%s,%d,%g,%g", tree->mod->modelname->s,
    //              tree->a_nodes[tax_id]->name, site_idx,
    //              log(state_probs[tax_id * ns +
    //                              tree->a_nodes[tax_id]->c_seq->d_state[0]]),
    //              tree->a_nodes[tax_id]->b[0]->l->v);
  }

  tree->io->init_len = ori_len;

  Free_Actual_CSeq(tree->data);

  tree->data = ori_data;

  Connect_CSeqs_To_Nodes(tree->data, tree->io, tree);

  return (state_probs);
}
