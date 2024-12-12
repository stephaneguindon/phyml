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
    char *dum;
    FILE *fp;

    io = NULL;

    io = (option *)Get_Input(argc, argv);
    if (!io)
        return (0);
    else if (io->use_xml == YES)
    {
        Free(io);
        return (0);
    }

    io->r_seed = (io->r_seed < 0) ? (time(NULL)) : (io->r_seed);
    srand(io->r_seed);

    io->data = Make_Empty_Alignment(io);

    Make_Model_Complete(io->mod);
    Set_Model_Name(io->mod);
    Print_Settings(io);

    io->colalias = NO;
    cdata = Compact_Data(io->data, io);
    Free_Seq(io->data, io->n_otu);

    tree = Make_Tree_From_Scratch(io->n_otu, cdata);

    Connect_CSeqs_To_Nodes(cdata, io, tree);

    tree = Make_Tree_From_Scratch(io->n_otu, cdata);

    Connect_CSeqs_To_Nodes(cdata, io, tree);

    tree->rates = RATES_Make_Rate_Struct(tree->n_otu);
    RATES_Init_Rate_Struct(tree->rates, io->rates, tree->n_otu);

    tree->times = TIMES_Make_Time_Struct(tree->n_otu);
    TIMES_Init_Time_Struct(tree->times, io->times, tree->n_otu);

    tree->data = cdata;
    tree->mod = io->mod;
    tree->io = io;
    tree->n_pattern = tree->data->crunch_len;
    tree->times->scaled_pop_size = 1.E-0;

    EVOLVE_Coalescent(tree);

    Make_Tree_For_Pars(tree);
    Make_Tree_For_Lk(tree);
    Make_Spr(tree);

    Evolve(tree->data, tree->mod, 0, tree);

    dum = (char *)mCalloc(100, sizeof(char));
    sprintf(dum, "%s%d%s", "./",io->r_seed, "_evolve_data.txt");
    fp = Openfile(dum, WRITE);
    Print_CSeq(fp, NO, tree->data, tree);
    fclose(fp);
    Free(dum);

    dum = (char *)mCalloc(100, sizeof(char));
    sprintf(dum, "%s%d%s", "./",io->r_seed, "_evolve_tree.txt");
    fp = Openfile(dum, WRITE);
    PhyML_Fprintf(fp,"%s\n",Write_Tree(tree));
    fclose(fp);
    Free(dum);

    return (-1);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void EVOLVE_Coalescent(t_tree *tree)
{
    t_node *n, *new_n, **avail;
    int n_lineages, idx_ancestor, *rand_idx, idx_edge;
    short int no_tip_left;
    phydbl Ne, dt, exp_g;
    phydbl T;

    Ne = tree->times->scaled_pop_size;
    exp_g = tree->times->neff_growth;

    avail = (t_node **)mCalloc(tree->n_otu, sizeof(t_node *));

    if (Ne > tree->times->scaled_pop_size_max || Ne < tree->times->scaled_pop_size_min)
        assert(FALSE);

    Get_Node_Ranks_From_Tip_Times(tree);

    n = tree->a_nodes[0];
    while (n->rk_next)
        n = n->rk_next;

    for (int i = 0; i < tree->n_otu; ++i)
        avail[i] = NULL;
    avail[0] = n;
    avail[1] = n->rk_prev;

    // PhyML_Printf("\n. n->time = %f n->rk_prev->time: %f", tree->times->nd_t[n->num], tree->times->nd_t[n->rk_prev->num]);
    idx_ancestor = tree->n_otu;
    idx_edge = 0;
    new_n = NULL;
    n = n->rk_prev;
    n_lineages = 2;
    no_tip_left = FALSE;
    T = tree->times->nd_t[n->num];
    do
    {
        dt = Rexp(Bico(n_lineages, 2)/Ne);
        T = T - dt;

        // PhyML_Printf("\n. T: %f", T);

        if (n->rk_prev == NULL)
            no_tip_left = TRUE;

        if (n->rk_prev != NULL && T < tree->times->nd_t[n->rk_prev->num])
        {
            T = tree->times->nd_t[n->rk_prev->num];
            if (n->rk_prev->tax == YES)
                n = n->rk_prev;
            avail[n_lineages] = n;
            n_lineages++;
            // PhyML_Printf("\n. Found tip %s at time %f", n->name, tree->times->nd_t[n->num]);
        }
        else
        {
            new_n = tree->a_nodes[idx_ancestor];
            tree->times->nd_t[idx_ancestor] = T;

            idx_ancestor++;

            rand_idx = Permutate(n_lineages);

            avail[rand_idx[0]]->v[0] = new_n;
            avail[rand_idx[1]]->v[0] = new_n;

            new_n->v[1] = avail[rand_idx[0]];
            new_n->v[2] = avail[rand_idx[1]];

            avail[rand_idx[0]] = new_n;
            avail[rand_idx[1]] = avail[n_lineages - 1];
            avail[n_lineages - 1] = NULL;

            Connect_One_Edge_To_Two_Nodes(new_n, new_n->v[1], tree->a_edges[idx_edge], tree);
            Connect_One_Edge_To_Two_Nodes(new_n, new_n->v[2], tree->a_edges[idx_edge + 1], tree);

            // PhyML_Printf("\n. Added coalescent node %d at time %f descendants %d %d # lineages: %d",
            //              idx_ancestor - 1, tree->times->nd_t[idx_ancestor - 1],
            //              new_n->v[1]->num,
            //              new_n->v[2]->num,
            //              n_lineages);

            n_lineages--;
            idx_edge += 2;
            Free(rand_idx);
            // for(int i=0; i < n_lineages; ++i) PhyML_Printf("\n. Avail: %d",avail[i]->num);
        }
    } while (n_lineages > 1 || no_tip_left == FALSE);

    tree->n_root = new_n;

    tree->n_root->v[1]->v[0] = tree->n_root->v[2];
    tree->n_root->v[2]->v[0] = tree->n_root->v[1];
    Connect_One_Edge_To_Two_Nodes(tree->n_root->v[1], tree->n_root->v[2], tree->a_edges[idx_edge], tree);
    tree->e_root = tree->a_edges[idx_edge];

    Update_Ancestors(tree->n_root, tree->n_root->v[2], tree->n_root->b[2], tree);
    Update_Ancestors(tree->n_root, tree->n_root->v[1], tree->n_root->b[1], tree);
    RATES_Fill_Lca_Table(tree);

    phydbl L = TIMES_Tree_Length(tree);
    tree->rates->bl_from_rt = YES;
    tree->rates->clock_r = 0.1 / L * (2 * tree->n_otu - 2);
    tree->rates->model_id = STRICTCLOCK;

    RATES_Update_Edge_Lengths(tree);

    Free(avail);
}