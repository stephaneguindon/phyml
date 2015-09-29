/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include <config.h>

#ifndef UTILITIES_H
#define UTILITIES_H


#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <time.h>
#include <limits.h>
#include <errno.h>
#include <float.h>
#include <assert.h>
#include <stdbool.h>

extern int n_sec1;
extern int n_sec2;

#define For(i,n)                     for(i=0; i<n; i++)
#define Fors(i,n,s)                  for(i=0; i<n; i+=s)
#define PointGamma(prob,alpha,beta)  PointChi2(prob,2.0*(alpha))/(2.0*(beta))
#define SHFT2(a,b,c)                 (a)=(b);(b)=(c);
#define SHFT3(a,b,c,d)               (a)=(b);(b)=(c);(c)=(d);
#define SIGN(a,b)                    ((b) > 0.0 ? fabs(a) : -fabs(a))
#define SHFT(a,b,c,d)                (a)=(b);(b)=(c);(c)=(d);

#ifndef MAX
#define MAX(a,b)                     ((a)>(b)?(a):(b))
#endif

#ifndef MIN
#define MIN(a,b)                     ((a)<(b)?(a):(b))
#endif

#define READ  0
#define WRITE 1

#ifndef isnan
# define isnan(x)						 \
  (sizeof (x) == sizeof (long double) ? isnan_ld (x)		 \
   : sizeof (x) == sizeof (double) ? isnan_d (x)		 \
   : isnan_f (x))
static inline int isnan_f  (float       x) { return x != x; }
static inline int isnan_d  (double      x) { return x != x; }
static inline int isnan_ld (long double x) { return x != x; }
#endif

#ifndef isinf
# define isinf(x)						 \
  (sizeof (x) == sizeof (long double) ? isinf_ld (x)		 \
   : sizeof (x) == sizeof (double) ? isinf_d (x)		 \
   : isinf_f (x))
static inline int isinf_f  (float       x) { return isnan (x - x); }
static inline int isinf_d  (double      x) { return isnan (x - x); }
static inline int isinf_ld (long double x) { return isnan (x - x); }
#endif

#define AC 0
#define AG 1
#define AT 2
#define CG 3
#define CT 4
#define GT 5

#ifndef M_1_SQRT_2PI
#define M_1_SQRT_2PI	0.398942280401432677939946059934	/* 1/sqrt(2pi) */
#endif



#define T_MAX_MCMC_MOVE_NAME 500

#define WINDOW_WIDTH  800
#define WINDOW_HEIGHT 800

#define RR_MIN 0.01
#define RR_MAX 200.0

#define PHYREX_UNIFORM 0
#define PHYREX_NORMAL  1

#define MCMC_MOVE_RANDWALK_UNIFORM       0
#define MCMC_MOVE_LOG_RANDWALK_UNIFORM   1
#define MCMC_MOVE_RANDWALK_NORMAL        2
#define MCMC_MOVE_LOG_RANDWALK_NORMAL    3
#define MCMC_MOVE_SCALE_THORNE           4
#define MCMC_MOVE_SCALE_GAMMA            5

#define N_MAX_MOVES     50

#define N_MAX_NEX_COM   20
#define T_MAX_NEX_COM   100
#define N_MAX_NEX_PARM  50
#define T_MAX_TOKEN     200
#define T_MAX_ID_COORD  5
#define T_MAX_ID_DISK   5

#define N_MAX_MIXT_CLASSES 1000

#define NEXUS_COM    0
#define NEXUS_PARM   1
#define NEXUS_EQUAL  2
#define NEXUS_VALUE  3
#define NEXUS_SPACE  4

#define  NNI_MOVE            0
#define  SPR_MOVE            1
#define  BEST_OF_NNI_AND_SPR 2

#define  M_1_SQRT_2_PI   0.398942280401432677939946059934
#define  M_SQRT_32       5.656854249492380195206754896838
#define  PI              3.14159265358979311600
#define  SQRT2PI         2.50662827463100024161
#define  LOG2PI          1.83787706640934533908
#define  LOG2            0.69314718055994528623
#define  LOG_SQRT_2_PI   0.918938533204672741780329736406

#define  NORMAL 1
#define  EXACT  2

#define  PHYLIP 0
#define  NEXUS  1

#ifndef YES
#define  YES 1
#endif

#ifndef NO
#define  NO  0
#endif

#ifndef TRUE
#define  TRUE  1
#endif

#ifndef FALSE
#define  FALSE 0
#endif

#define  ON  1
#define  OFF 0

#define SMALL_DBL 1.E-20

#define  NT 0 /*! nucleotides */
#define  AA 1 /*! amino acids */
#define  GENERIC 2 /*! custom alphabet */
#define  UNDEFINED -1

#define  ACGT 0 /*! A,G,G,T encoding */
#define  RY   1 /*! R,Y     encoding */

#define INTERFACE_DATA_TYPE      0
#define INTERFACE_MULTIGENE      1
#define INTERFACE_MODEL          2
#define INTERFACE_TOPO_SEARCH    3
#define INTERFACE_BRANCH_SUPPORT 4

#define LEFT 0
#define RGHT 1


#ifndef INFINITY
#define INFINITY HUGE
#endif

#define  N_MAX_OPTIONS        100

#define NEXT_BLOCK_SIZE        50

#define  T_MAX_FILE           200
#define  T_MAX_LINE       2000000
#define  T_MAX_NAME           100
#define  T_MAX_ID              20
#define  T_MAX_SEQ        2000000
#define  T_MAX_OPTION         100
#define  T_MAX_LABEL           10
#define  T_MAX_STATE            5
#define  N_MAX_LABEL           10
#define  BLOCK_LABELS         100

#define  NODE_DEG_MAX         500
#define  BRENT_IT_MAX         500
#define  BRENT_CGOLD    0.3819660
#define  BRENT_ZEPS        1.e-10
#define  MNBRAK_GOLD     1.618034
#define  MNBRAK_GLIMIT      100.0
#define  MNBRAK_TINY       1.e-20
#define  ALPHA_MIN           0.04
#define  ALPHA_MAX            100
#define  BL_START          1.e-04
#define  GOLDEN_R      0.61803399
#define  GOLDEN_C  (1.0-GOLDEN_R)
#define  N_MAX_INSERT          20
#define  N_MAX_OTU           4000
#define  UNLIKELY          -1.e10
#define  NJ_SEUIL             0.1
#define  ROUND_MAX            100
#define  DIST_MAX            2.00
#define  AROUND_LK           50.0
#define  PROP_STEP            1.0
#define  T_MAX_ALPHABET        22
#define  MDBL_MIN         FLT_MIN
#define  MDBL_MAX         FLT_MAX
#define  POWELL_ITMAX         200
#define  LINMIN_TOL       2.0E-04
#define  SCALE_POW             10    /*! Scaling factor will be 2^SCALE_POW or 2^(-SCALE_POW) [[ WARNING: SCALE_POW < 31 ]]*/
#define  DEFAULT_SIZE_SPR_LIST 20
#define  NEWICK                 0
#define  NEXUS                  1
#define  OUTPUT_TREE_FORMAT NEWICK
#define  MAX_PARS      1000000000

#define  LIM_SCALE_VAL     1.E-50 /*! Scaling limit (deprecated) */

#define  MIN_CLOCK_RATE   1.E-10

#define  MIN_VAR_BL       1.E-8
#define  MAX_VAR_BL       1.E+3

#define JC69       1
#define K80        2
#define F81        3
#define HKY85      4
#define F84        5
#define TN93       6
#define GTR        7
#define CUSTOM     8

#define WAG       11
#define DAYHOFF   12
#define JTT       13
#define BLOSUM62  14
#define MTREV     15
#define RTREV     16
#define CPREV     17
#define DCMUT     18
#define VT        19
#define MTMAM     20
#define MTART     21
#define HIVW      22
#define HIVB      23
#define FLU       24
#define CUSTOMAA  25
#define LG        26
#define AB        27
// Amino acid ordering:
// Ala Arg Asn Asp Cys Gln Glu Gly His Ile Leu Lys Met Phe Pro Ser Thr Trp Tyr Val

#define COMPOUND_COR   0
#define COMPOUND_NOCOR 1
#define EXPONENTIAL    2
#define GAMMA          3
#define THORNE         4
#define GUINDON        5
#define STRICTCLOCK    6
#define NONE          -1

#define ALRTSTAT       1
#define ALRTCHI2       2
#define MINALRTCHI2SH  3
#define SH             4
#define ABAYES         5


/*  /\* Uncomment the lines below to switch to single precision *\/ */
/*  typedef	float phydbl; */
/*  #define LOG logf */
/*  #define POW powf */
/*  #define EXP expf */
/*  #define FABS fabsf */
/*  #define SQRT sqrtf */
/*  #define CEIL ceilf */
/*  #define FLOOR floorf */
/*  #define RINT rintf */
/*  #define ROUND roundf */
/*  #define TRUNC truncf */
/*  #define COS cosf */
/*  #define SIN sinf */
/*  #define TAN tanf */
/*  #define SMALL FLT_MIN */
/*  #define BIG  FLT_MAX */
/*  #define SMALL_PIJ 1.E-10 */
/*  #define BL_MIN 1.E-5 */
/*  #define  P_LK_LIM_INF   2.168404e-19 /\* 2^-62 *\/ */
/*  #define  P_LK_LIM_SUP   4.611686e+18 /\* 2^62 *\/ */

/* Uncomment the line below to switch to double precision */
typedef	double phydbl;
#define LOG log
#define POW pow
#define EXP exp
#define FABS fabs
#define SQRT sqrt
#define CEIL ceil
#define FLOOR floor
#define RINT rint
#define ROUND round
#define TRUNC trunc
#define COS cos
#define SIN sin
#define TAN tan
#define SMALL DBL_MIN
#define BIG  DBL_MAX
#define SMALL_PIJ 1.E-20
#define LOGBIG 690.
#define LOGSMALL -690.


#if !(defined PHYTIME || defined SERGEII)
#define BL_MIN 1.E-8
#define BL_MAX 100.
#else
#define BL_MIN 1.E-6
#define BL_MAX 1.
#endif

/* #define P_LK_LIM_INF 7.888609052e-31 */
/* #define P_LK_LIM_MAX 1.267650600e+30 */
/* #define P_LK_LIM_INF 4.909093465e-91 /\* R: format(2^(-300),digits=10) *\/ */
/* #define P_LK_LIM_SUP 2.037035976e+90 /\* R: format(2^(+300),digits=10) *\/ */
#define  P_LK_LIM_INF   3.054936e-151 /* 2^-500 */
#define  P_LK_LIM_SUP   3.273391e+150 /* 2^500 */
//#define  P_LK_LIM_INF   4.656612873e-10 /*2^-31 */

#define T_MAX_XML_TAG 64

#define NARGS_SEQ(_1,_2,_3,_4,_5,_6,_7,_8,N,...) N
#define NARGS(...) NARGS_SEQ(__VA_ARGS__, 8, 7, 6, 5, 4, 3, 2, 1)
#define PRIMITIVE_CAT(x, y) x ## y
#define CAT(x, y) PRIMITIVE_CAT(x, y)
#define FOR_EACH(macro, ...) CAT(FOR_EACH_, NARGS(__VA_ARGS__))(macro, __VA_ARGS__)
#define FOR_EACH_1(m, x1) m(x1)
#define FOR_EACH_2(m, x1, x2) m(x1) m(x2)
#define FOR_EACH_3(m, x1, x2, x3) m(x1) m(x2) m(x3)
#define FOR_EACH_4(m, x1, x2, x3, x4) m(x1) m(x2) m(x3) m(x4)
#define FOR_EACH_5(m, x1, x2, x3, x4, x5) m(x1) m(x2) m(x3) m(x4) m(x5)
#define FOR_EACH_6(m, x1, x2, x3, x4, x5, x6) m(x1) m(x2) m(x3) m(x4) m(x5) m(x6)
#define FOR_EACH_7(m, x1, x2, x3, x4, x5, x6, x7) m(x1) m(x2) m(x3) m(x4) m(x5) m(x6) m(x7)
#define FOR_EACH_8(m, x1, x2, x3, x4, x5, x6, x7, x8) m(x1) m(x2) m(x3) m(x4) m(x5) m(x6) m(x7) m(x8)
#define DUMP_EACH_INT(v)     fprintf(stderr,"\n\t\tDEBUG:%s:\t\t%s--->%i",__PRETTY_FUNCTION__,#v,(v));fflush(stderr);
#define DUMP_EACH_STRING(v)  fprintf(stderr,"\n\t\tDEBUG:%s:\t\t%s--->%s",__PRETTY_FUNCTION__,#v,(v));fflush(stderr);
#define DUMP_EACH_DECIMAL(v) fprintf(stderr,"\n\t\tDEBUG:%s:\t\t%s--->%f",__PRETTY_FUNCTION__,#v,(v));fflush(stderr);
#define DUMP_EACH_SCI(v)     fprintf(stderr,"\n\t\tDEBUG:%s:\t\t%s--->%e",__PRETTY_FUNCTION__,#v,(v));fflush(stderr);
#define DUMP_I(...) FOR_EACH(DUMP_EACH_INT, __VA_ARGS__)
#define DUMP_S(...) FOR_EACH(DUMP_EACH_STRING, __VA_ARGS__)
#define DUMP_D(...) FOR_EACH(DUMP_EACH_DECIMAL, __VA_ARGS__)
#define DUMP_E(...) FOR_EACH(DUMP_EACH_SCI, __VA_ARGS__)

/*!********************************************************/

typedef struct __Scalar_Int {
  int                      v;
  struct __Scalar_Int  *next;
  struct __Scalar_Int  *prev;
}scalar_int;


/*!********************************************************/

typedef struct __Scalar_Dbl {
  phydbl                  v;
  bool                onoff;
  struct __Scalar_Dbl *next;
  struct __Scalar_Dbl *prev;
}scalar_dbl;

/*!********************************************************/

typedef struct __Vect_Int {
  int                  *v;
  int                 len;
  struct __Vect_Int *next;
  struct __Vect_Int *prev;
}vect_int;


/*!********************************************************/

typedef struct __Vect_Dbl {
  phydbl               *v;
  int                 len;
  struct __Vect_Dbl *next;
  struct __Vect_Dbl *prev;
}vect_dbl;

/*!********************************************************/

typedef struct __String {
  char                 *s;
  int                 len;
  struct __String   *next;
  struct __String   *prev;
}t_string;

/*!********************************************************/

typedef struct __Node {
  struct __Node                       **v; /*! table of pointers to neighbor nodes. Dimension = 3 */
  struct __Node               ***bip_node; /*! three lists of pointer to tip nodes. One list for each direction */
  struct __Edge                       **b; /*! table of pointers to neighbor branches */
  struct __Node                      *anc; /*! direct ancestor t_node (for rooted tree only) */
  struct __Node                 *ext_node;
  struct __Node               *match_node;
  struct __Align                   *c_seq; /*! corresponding compressed sequence */
  struct __Align               *c_seq_anc; /*! corresponding compressed ancestral sequence */


  struct __Node                     *next; /*! tree->a_nodes[i]->next <=> tree->next->a_nodes[i] */
  struct __Node                   *prev; /*! See above */
  struct __Node                *next_mixt; /*! Next mixture tree*/
  struct __Node                *prev_mixt; /*! Parent mixture tree */

  struct __Calibration            **calib;
  short int             *calib_applies_to;
  int                             n_calib;

  int                           *bip_size; /*! Size of each of the three lists from bip_node */
  int                                 num; /*! t_node number */
  int                                 tax; /*! tax = 1 -> external node, else -> internal t_node */
  int                        check_branch; /*! check_branch=1 is the corresponding branch is labelled with '*' */
  char                              *name; /*! taxon name (if exists) */
  char                          *ori_name; /*! taxon name (if exists) */

  phydbl                           *score; /*! score used in BioNJ to determine the best pair of nodes to agglomerate */
  phydbl                               *l; /*! lengths of the (three or one) branche(s) connected this t_node */
  phydbl                     dist_to_root; /*! distance to the root t_node */

  short int                        common;
  phydbl                           y_rank;
  phydbl                       y_rank_ori;
  phydbl                       y_rank_min;
  phydbl                       y_rank_max;

  int                            *s_ingrp; /*! does the subtree beneath belong to the ingroup */
  int                           *s_outgrp; /*! does the subtree beneath belong to the outgroup */

  int                             id_rank; /*! order taxa alphabetically and use id_rank to store the ranks */
  int                                rank;
  int                            rank_max;

  struct __Geo_Coord               *coord;

}t_node;


/*!********************************************************/

typedef struct __Edge {
  /*!
    syntax :  (node) [edge]
(left_1) .                   .(right_1)
          \ (left)  (right) /
           \._____________./
           /    [b_fcus]   \
          /                 \
(left_2) .                   .(right_2)

  */

  struct __Node               *left,*rght; /*! t_node on the left/right side of the t_edge */
  short int         l_r,r_l,l_v1,l_v2,r_v1,r_v2;
  /*! these are directions (i.e., 0, 1 or 2): */
  /*! l_r (left to right) -> left[b_fcus->l_r] = right */
  /*! r_l (right to left) -> right[b_fcus->r_l] = left */
  /*! l_v1 (left t_node to first t_node != from right) -> left[b_fcus->l_v1] = left_1 */
  /*! l_v2 (left t_node to secnd t_node != from right) -> left[b_fcus->l_v2] = left_2 */
  /*! r_v1 (right t_node to first t_node != from left) -> right[b_fcus->r_v1] = right_1 */
  /*! r_v2 (right t_node to secnd t_node != from left) -> right[b_fcus->r_v2] = right_2 */

  struct __NNI                       *nni;
  struct __Edge                     *next;
  struct __Edge                     *prev;
  struct __Edge                *next_mixt;
  struct __Edge                *prev_mixt;

  int                                 num; /*! branch number */
  scalar_dbl                           *l; /*! branch length */
  scalar_dbl                      *best_l; /*! best branch length found so far */
  scalar_dbl                       *l_old; /*! old branch length */
  scalar_dbl                       *l_var; /*! variance of edge length */
  scalar_dbl                   *l_var_old; /*! variance of edge length (previous value) */


  int                           bip_score; /*! score of the bipartition generated by the corresponding edge
                          bip_score = 1 iif the branch is found in both trees to be compared,
                          bip_score = 0 otherwise. */

  phydbl            *p_lk_left,*p_lk_rght; /*! likelihoods of the subtree on the left and right side (for each site and each relative rate category) */
  short int      *p_lk_tip_r, *p_lk_tip_l;
#ifdef BEAGLE
  int        p_lk_left_idx, p_lk_rght_idx;
  int                        p_lk_tip_idx;
#endif

  short int           *div_post_pred_left; /*! posterior prediction of nucleotide/aa diversity (left-hand subtree) */
  short int           *div_post_pred_rght; /*! posterior prediction of nucleotide/aa diversity (rght-hand subtree) */
  short int                    does_exist;


  int                       *patt_id_left;
  int                       *patt_id_rght;
  int                      *p_lk_loc_left;
  int                      *p_lk_loc_rght;


  phydbl                          *Pij_rr; /*! matrix of change probabilities and its first and secnd derivates (rate*state*state) */
#ifdef BEAGLE
  int                          Pij_rr_idx;
#endif
  int                     *pars_l,*pars_r; /*! parsimony of the subtree on the left and right sides (for each site) */
  unsigned int               *ui_l, *ui_r; /*! union - intersection vectors used in Fitch's parsimony algorithm */
  int                *p_pars_l, *p_pars_r; /*! conditional parsimony vectors */

  int                         num_st_left; /*! number of the subtree on the left side */
  int                         num_st_rght; /*! number of the subtree on the right side */


  /*! Below are the likelihood scaling factors (used in functions
     `Get_All_Partial_Lk_Scale' in lk.c. */
  /*
For every site, every subtree and every rate class, PhyML maintains a`sum_scale_pow' value
where sum_scale_pow = sum_scale_pow_v1 + sum_scale_pow_v2 + curr_scale_pow'
sum_scale_pow_v1 and sum_scale_pow_v2 are sum_scale_pow of the left and
right subtrees. curr_scale_pow is an integer greater than one.
The smaller the partials, the larger curr_scale_pow.

Now the partials for this subtree are scaled by *multiplying* each of
them by 2^curr_scale_pow. The reason for doing the scaling this way is
that multiplications by 2^x (x an integer) can be done in an 'exact'
manner (i.e., there is no loss of numerical precision)

At the root edge, the log-likelihood is then
    logL = logL' - (sum_scale_pow_left + sum_scale_pow_right)log(2),
where L' is the scaled likelihood.
*/
  int                 *sum_scale_left_cat;
  int                 *sum_scale_rght_cat;
  int                     *sum_scale_left;
  int                     *sum_scale_rght;

  phydbl                          bootval; /*! bootstrap value (if exists) */

  short int                      is_alive; /*! is_alive = 1 if this t_edge is used in a tree */

  phydbl                   dist_btw_edges;
  int                 topo_dist_btw_edges;

  int                     has_zero_br_len;

  phydbl                       ratio_test; /*! approximate likelihood ratio test */
  phydbl                   alrt_statistic; /*! aLRT statistic */

  char                           **labels; /*! string of characters that labels the corresponding t_edge */
  int                            n_labels; /*! number of labels */
  int                             n_jumps; /*! number of jumps of substitution rates */

  int                    *n_diff_states_l; /*! Number of different states found in the subtree on the left of this edge */
  int                    *n_diff_states_r; /*! Number of different states found in the subtree on the right of this edge */

  phydbl                      bin_cod_num;
}t_edge;

/*!********************************************************/

typedef struct __Tree{

  struct __Node                       *n_root; /*! root t_node */
  struct __Edge                       *e_root; /*! t_edge on which lies the root */
  struct __Node                     **a_nodes; /*! array of nodes that defines the tree topology */
  struct __Edge                     **a_edges; /*! array of edges */
  struct __Model                         *mod; /*! substitution model */
  struct __Calign                       *data; /*! sequences */
  struct __Tree                         *next; /*! set to NULL by default. Used for mixture models */
  struct __Tree                         *prev; /*! set to NULL by default. Used for mixture models */
  struct __Tree                    *next_mixt; /*! set to NULL by default. Used for mixture models */
  struct __Tree                    *prev_mixt; /*! set to NULL by default. Used for mixture models */
  struct __Tree                    *mixt_tree; /*! set to NULL by default. Used for mixture models */
  struct __Option                         *io; /*! input/output */
  struct __Matrix                        *mat; /*! pairwise distance matrix */
  struct __Node                   **curr_path; /*! list of nodes that form a path in the tree */
  struct __SPR                     **spr_list;
  struct __SPR                      *best_spr;
  struct __Tdraw                     *ps_tree; /*! structure for drawing trees in postscript format */
  struct __T_Rate                       *rates; /*! structure for handling rates of evolution */
  struct __Tmcmc                        *mcmc;
  struct __Triplet            *triplet_struct;
  struct __Phylogeo                      *geo;
  struct __Migrep_Model                 *mmod;
  struct __Disk_Event                   *disk;

  int                            is_mixt_tree;
  int                                tree_num; /*! tree number. Used for mixture models */
  int                          ps_page_number; /*! when multiple trees are printed, this variable give the current page number */
  int                         depth_curr_path; /*! depth of the t_node path defined by curr_path */
  int                                 has_bip; /*!if has_bip=1, then the structure to compare
                         tree topologies is allocated, has_bip=0 otherwise */
  int                                   n_otu; /*! number of taxa */
  int                               curr_site; /*! current site of the alignment to be processed */
  int                               curr_catg; /*! current class of the discrete gamma rate distribution */
  int                                  n_swap; /*! number of NNIs performed */
  int                               n_pattern; /*! number of distinct site patterns */
  int                      has_branch_lengths; /*! =1 iff input tree displays branch lengths */
  int                          print_boot_val; /*! if print_boot_val=1, the bootstrap values are printed */
  int                          print_alrt_val; /*! if print_boot_val=1, the aLRT values are printed */
  int                              both_sides; /*! both_sides=1 -> a pre-order and a post-order tree
                          traversals are required to compute the likelihood
                          of every subtree in the phylogeny*/
  int               num_curr_branch_available; /*!gives the number of the next cell in a_edges that is free to receive a pointer to a branch */
  short int                            *t_dir;
  int                          n_improvements;
  int                                 n_moves;

  int                                      dp; /*! Data partition */
  int                               s_mod_num; /*! Substitution model number */
  int                               lock_topo; /*! = 1 any subsequent topological modification will be banished */
  int                            write_labels;
  int                           write_br_lens;
  int                                 *mutmap; /*! Mutational map */

  phydbl                             init_lnL;
  phydbl                             best_lnL; /*! highest value of the loglikelihood found so far */
  int                               best_pars; /*! highest value of the parsimony found so far */
  phydbl                                c_lnL; /*! loglikelihood */
  phydbl                              old_lnL; /*! old loglikelihood */
  phydbl                    sum_min_sum_scale; /*! common factor of scaling factors */
  phydbl                        *c_lnL_sorted; /*! used to compute c_lnL by adding sorted terms to minimize CPU errors */
  phydbl                         *cur_site_lk; /*! vector of loglikelihoods at individual sites */
  phydbl                         *old_site_lk; /*! vector of likelihoods at individual sites */
  phydbl                       annealing_temp; /*! annealing temperature in simulated annealing optimization algo */

  phydbl                *unscaled_site_lk_cat; /*! partially scaled site likelihood at individual sites */

  phydbl                         *site_lk_cat; /*! loglikelihood at a single site and for each class of rate*/
  phydbl                      unconstraint_lk; /*! unconstrained (or multinomial) likelihood  */
  int                         *fact_sum_scale;
  phydbl                       **log_lks_aLRT; /*! used to compute several branch supports */
  phydbl                           n_root_pos; /*! position of the root on its t_edge */
  phydbl                                 size; /*! tree size */
  int                              *site_pars;
  int                                  c_pars;
  int                               *step_mat;

  int                           size_spr_list;
  int                  perform_spr_right_away;

  time_t                                t_beg;
  time_t                            t_current;

  int                     bl_from_node_stamps; /*! == 1 -> Branch lengths are determined by t_node times */
  phydbl                        sum_y_dist_sq;
  phydbl                           sum_y_dist;
  phydbl                      tip_order_score;
  int                         write_tax_names;
  int                    update_alias_subpatt;

  phydbl                           geo_mig_sd; /*! standard deviation of the migration step random variable */
  phydbl                              geo_lnL; /*! log likelihood of the phylo-geography model */

  int                              bl_ndigits;

  phydbl                             *short_l; /*! Vector of short branch length values */
  int                               n_short_l; /*! Length of short_l */
  phydbl                           norm_scale;

  short int                   br_len_recorded;
  int                           max_spr_depth;

  short int                  apply_lk_scaling; /*! Applying scaling of likelihoods. YES/NO */

  phydbl                                   *K; /*! a vector of the norm.constants for the node times prior. */

  short int                       ignore_root;
#ifdef BEAGLE
  int                                  b_inst; /*! The BEAGLE instance id associated with this tree. */
#endif

}t_tree;

/*!********************************************************/

typedef struct __Super_Tree {
  struct __Tree                           *tree;
  struct __List_Tree                  *treelist; /*! list of trees. One tree for each data set to be processed */
  struct __Calign                    *curr_cdata;
  struct __Option                   **optionlist; /*! list of pointers to input structures (used in supertrees) */

  struct __Node           ***match_st_node_in_gt;
  /*!  match_st_in_gt_node[subdataset number][supertree t_node number]
   *  gives the t_node in tree estimated from 'subdataset number' that corresponds
   *  to 'supertree t_node number' in the supertree
   */

  struct __Node           *****map_st_node_in_gt;
  /*!  mat_st_gt_node[gt_num][st_node_num][direction] gives the
   *  t_node in gt gt_num that maps t_node st_node_num in st.
   */

  struct __Edge             ***map_st_edge_in_gt;
  /*!  map_st_gt_br[gt_num][st_branch_num] gives the
   *  branch in gt gt_num that maps branch st_branch_num
   *  in st.
   */

  struct __Edge            ****map_gt_edge_in_st;
  /*!  mat_gt_st_br[gt_num][gt_branch_num][] is the list of
   *  branches in st that map branch gt_branch_num
   *  in gt gt_num.
   */

  int                   **size_map_gt_edge_in_st;
  /*!  size_map_gt_st_br[gt_num][gt_branch_num] gives the
   *  size of the list map_gt_st_br[gt_num][gt_branch_num][]
   */


  struct __Edge             ***match_st_edge_in_gt;
  /*! match_st_edge_in_gt[gt_num][st_branch_num] gives the
   * branch in gt gt_num that matches branch st_branch_num
   */

  struct __Edge             ***match_gt_edge_in_st;
  /*! match_gt_edge_in_st[gt_num][gt_branch_num] gives the
   * branch in st that matches branch gt_branch_num
   */

  struct __Node                  ****closest_match;
  /*! closest_match[gt_num][st_node_num][dir] gives the
   * closest t_node in st that matches a t_node in gt gt_num
   */

  int                              ***closest_dist;
  /*! closest_dist[gt_num][st_node_num][dir] gives the
   * number of edges to traverse to get to node
   * closest_match[gt_num][st_node_num][dir]
   */

  int                                         n_part;
  /*! number of trees */

  phydbl                                      **bl;
  /*! bl[gt_num][gt_branch_num] gives the length of
   * branch gt_branch_num
   */

  phydbl                                  **bl_cpy;
  /*! copy of bl */

  phydbl                                     **bl0;
  /*! bl estimated during NNI (original topo)
   * See Mg_NNI.
   */

  phydbl                                     **bl1;
  /*! bl estimated during NNI (topo conf 1)
   * See Mg_NNI.
   */

  phydbl                                     **bl2;
  /*! bl estimated during NNI (topo conf 2)
   * See Mg_NNI.
   */

  int                                *bl_partition;
  /*! partition[gt_num] gives the t_edge partition number
   * gt_num belongs to.
   */
  int                                   n_bl_part;

  struct __Model                          **s_mod; /*! substitution model */

  int                                     n_s_mod;
  int                                 lock_br_len;

}supert_tree;

/*!********************************************************/

typedef struct __List_Tree { /*! a list of trees */
  struct __Tree   **tree;
  int           list_size;                /*! number of trees in the list */
}t_treelist;

/*!********************************************************/

typedef struct __Align {
  char           *name; /*! sequence name */
  int              len; /*! sequence length */
  char          *state; /*! sequence itself */
  int         *d_state; /*! sequence itself (digits) */
  short int *is_ambigu; /*! is_ambigu[site] = 1 if state[site] is an ambiguous character. 0 otherwise */
}align;

/*!********************************************************/


typedef struct __Calign {
  struct __Align **c_seq;             /*! compressed sequences      */
  struct __Option    *io;             /*! input/output */
  phydbl          *b_frq;             /*! observed state frequencies */
  short int       *invar;             /*! < 0 -> polymorphism observed */
  int              *wght;             /*! # of each site in c_align */
  short int      *ambigu;             /*! ambigu[i]=1 is one or more of the sequences at site
                                        i display an ambiguous character */
  phydbl      obs_pinvar;
  int              n_otu;             /*! number of taxa */
  int          clean_len;             /*! uncrunched sequences lenghts without gaps */
  int         crunch_len;             /*! crunched sequences lengths */
  int           init_len;             /*! length of the uncompressed sequences */
  int          *sitepatt;             /*! this array maps the position of the patterns in the
                                        compressed alignment to the positions in the uncompressed
                                        one */
  int             format;             /*! 0 (default): PHYLIP. 1: NEXUS. */
}calign;

/*!********************************************************/

typedef struct __Matrix { /*! mostly used in BIONJ */
  phydbl    **P,**Q,**dist; /*! observed proportions of transition, transverion and  distances
                   between pairs of  sequences */

  t_tree             *tree; /*! tree... */
  int              *on_off; /*! on_off[i]=1 if column/line i corresponds to a t_node that has not
                   been agglomerated yet */
  int                n_otu; /*! number of taxa */
  char              **name; /*! sequence names */
  int                    r; /*! number of nodes that have not been agglomerated yet */
  struct __Node **tip_node; /*! array of pointer to the leaves of the tree */
  int             curr_int; /*! used in the NJ/BIONJ algorithms */
  int               method; /*! if method=1->NJ method is used, BIONJ otherwise */
}matrix;

/*!********************************************************/

typedef struct __RateMatrix {
  int                    n_diff_rr; /*! number of different relative substitution rates in the custom model */
  vect_dbl                     *rr; /*! relative rate parameters of the GTR or custom model (given by rr_val[rr_num[i]]) */
  vect_dbl                 *rr_val; /*! relative rate parameters of the GTR or custom model */
  vect_int                 *rr_num;
  vect_int           *n_rr_per_cat; /*! number of rate parameters in each category */
  vect_dbl                   *qmat;
  vect_dbl              *qmat_buff;

  struct __RateMatrix        *next;
  struct __RateMatrix        *prev;
}t_rmat;

/*!********************************************************/

typedef struct __RAS {
  /*! Rate across sites */
  int                       n_catg; /*! number of categories in the discrete gamma distribution */
  int                        invar; /*! =1 iff the substitution model takes into account invariable sites */
  int                 gamma_median; /*! 1: use the median of each bin in the discrete gamma distribution. 0: the mean is used */
  vect_dbl          *gamma_r_proba; /*! probabilities of the substitution rates defined by the discrete gamma distribution */
  vect_dbl *gamma_r_proba_unscaled;
  vect_dbl               *gamma_rr; /*! substitution rates defined by the RAS distribution */
  vect_dbl      *gamma_rr_unscaled; /*! substitution rates defined by the RAS distribution (unscaled) */
  scalar_dbl                *alpha; /*! gamma shapa parameter */
  int              free_mixt_rates;
  scalar_dbl         *free_rate_mr; /*! mean relative rate as given by the FreeRate model */


  int          parent_class_number;
  scalar_dbl               *pinvar; /*! proportion of invariable sites */

  short int                init_rr;
  short int           init_r_proba;
  short int           normalise_rr;

  short int         *skip_rate_cat; /*! indicates whether of not the the likelihood for a given rate class shall be calculated.
                                        The default model in PhyML has four rate classes and the likelihood at a given site is
                                        then \sum_{i=1}^{4} \Pr(D|R=i) \Pr(R=i). Now, when optimising
                                        the rate for say the first rate class (i=1) one does not need to
                                        re-compute \Pr(D|R=2), \Pr(D|R=3) and \Pr(D|R=4) (which is the time
                                        consuming part of the likelihood calculation). This is where
                                        'skip_rate_category' comes in handy. */

  short int      sort_rate_classes; /*! When doing MCMC moves on rate classes, one needs to have the rate classes sorted */

  struct __RAS               *next;
  struct __RAS               *prev;

}t_ras;

/*!********************************************************/

typedef struct __EquFreq {
  /*! Equilibrium frequencies */
  vect_dbl             *pi; /*! states frequencies */
  vect_dbl    *pi_unscaled; /*! states frequencies (unscaled) */

  struct __EquFreq   *next;
  struct __EquFreq   *prev;

}t_efrq;

/*!********************************************************/

typedef struct __Model {
  struct __Optimiz          *s_opt; /*! pointer to parameters to optimize */
  struct __Eigen            *eigen;
  struct __M4               *m4mod;
  struct __Option              *io;
  struct __Model             *next;
  struct __Model             *prev;
  struct __Model        *next_mixt;
  struct __Model        *prev_mixt;
  struct __RateMatrix       *r_mat;
  struct __EquFreq          *e_frq;
  struct __RAS                *ras;

  t_string       *aa_rate_mat_file;
  FILE             *fp_aa_rate_mat;

  vect_dbl            *user_b_freq; /*! user-defined nucleotide frequencies */
  t_string              *modelname;
  t_string      *custom_mod_string; /*! string of characters used to define custom models of substitution */

  int                      mod_num; /*! model number */

  int                 update_eigen; /*! update_eigen=1-> eigen values/vectors need to be updated */

  int                   whichmodel;
  int                  is_mixt_mod;
  int                    augmented;
  int                           ns; /*! number of states (4 for ADN, 20 for AA) */

  int                    bootstrap; /*! Number of bootstrap replicates (0 : no bootstrap analysis is launched) */
  int                    use_m4mod; /*! Use a Makrkov modulated Markov model ? */

  scalar_dbl                *kappa; /*! transition/transversion rate */
  scalar_dbl               *lambda; /*! parameter used to define the ts/tv ratios in the F84 and TN93 models */
  scalar_dbl          *br_len_mult; /*! when users want to fix the relative length of edges and simply estimate the total length of the tree. This multiplier does the trick */
  scalar_dbl *br_len_mult_unscaled;

  vect_dbl                 *Pij_rr; /*! matrix of change probabilities */
  scalar_dbl                   *mr; /*! mean rate = branch length/time interval  mr = -sum(i)(vct_pi[i].mat_Q[ii]) */

  short int                 log_l; /*! Edge lengths are actually log(Edge lengths) if log_l == YES !*/
  phydbl                    l_min; /*! Minimum branch length !*/
  phydbl                    l_max; /*! Maximum branch length !*/

  phydbl              l_var_sigma; /*! For any edge b we have b->l_var->v = l_var_sigma * (b->l->v)^2 */
  phydbl                l_var_min; /*! Min of variance of branch lengths (used in conjunction with gamma_mgf_bl == YES) */
  phydbl                l_var_max; /*! Max of variance of branch lengths (used in conjunction with gamma_mgf_bl == YES) */

  int                gamma_mgf_bl; /*! P = \int_0^inf exp(QL) p(L) where L=\int_0^t R(s) ds and p(L) is the gamma density. Set to NO by default !*/

  int              n_mixt_classes; /* Number of classes in the mixture model. */

  scalar_dbl        *r_mat_weight;
  scalar_dbl        *e_frq_weight;
#ifdef BEAGLE
  int                      b_inst;
  bool        optimizing_topology; /*! This is a flag that prevents the resetting of category weights. Why? Read */
                                   /*  Recall that while optimizing the toplogy, PhyML temporarily only uses 2
                                    *  rate categories. Recall also that a BEAGLE instance is created with all the
                                    *  required categories, but we temporarily assign 0 weight to the other categories
                                    *  thus effectively using only 2 categories. However, subsequent calls to
                                    *  update the rates (i.e. update_beagle_ras()) will reset the weights. This flag
                                    *  prevents this resetting from happening                                    */
#endif

}t_mod;



/*!********************************************************/

typedef struct __Eigen{
  int              size; /*! matrix is size * size */
  phydbl             *q; /*! matrix for which eigen values and vectors are computed */
  phydbl         *space;
  int        *space_int;
  phydbl         *e_val; /*! eigen values (vector), real part. */
  phydbl      *e_val_im; /*! eigen values (vector), imaginary part */
  phydbl      *r_e_vect; /*! right eigen vector (matrix), real part */
  phydbl   *r_e_vect_im; /*! right eigen vector (matrix), imaginary part */
  phydbl      *l_e_vect; /*! left eigen vector (matrix), real part */

  struct __Eigen  *prev;
  struct __Eigen  *next;
}eigen;

/*!********************************************************/

typedef struct __Option { /*! mostly used in 'help.c' */
  struct __Model                *mod; /*! pointer to a substitution model */
  struct __Tree                *tree; /*! pointer to the current tree */
  struct __Align              **data; /*! pointer to the uncompressed sequences */
  struct __Tree           *cstr_tree; /*! pointer to a constraint tree (can be a multifurcating one) */
  struct __Calign             *cdata; /*! pointer to the compressed sequences */
  struct __Super_Tree            *st; /*! pointer to supertree */
  struct __Tnexcom    **nex_com_list;
  struct __List_Tree       *treelist; /*! list of trees. */
  struct __Option              *next;
  struct __Option              *prev;

  int                    interleaved; /*! interleaved or sequential sequence file format ? */
  int                        in_tree; /*! =1 iff a user input tree is used as input */

  char                *in_align_file; /*! alignment file name */
  FILE                  *fp_in_align; /*! pointer to the alignment file */

  char                 *in_tree_file; /*! input tree file name */
  FILE                   *fp_in_tree; /*! pointer to the input tree file */

  char      *in_constraint_tree_file; /*! input constraint tree file name */
  FILE        *fp_in_constraint_tree; /*! pointer to the input constraint tree file */

  char                *out_tree_file; /*! name of the tree file */
  FILE                  *fp_out_tree;

  char               *out_trees_file; /*! name of the tree file */
  FILE                 *fp_out_trees; /*! pointer to the tree file containing all the trees estimated using random starting trees */

  char           *out_boot_tree_file; /*! name of the tree file */
  FILE             *fp_out_boot_tree; /*! pointer to the bootstrap tree file */

  char          *out_boot_stats_file; /*! name of the tree file */
  FILE            *fp_out_boot_stats; /*! pointer to the statistics file */

  char               *out_stats_file; /*! name of the statistics file */
  FILE                 *fp_out_stats;

  char               *out_trace_file; /*! name of the file in which the likelihood of the model is written */
  FILE                 *fp_out_trace;

  char                  *out_lk_file; /*! name of the file in which the likelihood of the model is written */
  FILE                    *fp_out_lk;

  char             *out_summary_file; /*! name of the file in which summary statistics are written */
  FILE               *fp_out_summary;

  char                  *out_ps_file; /*! name of the file in which tree(s) is(are) written */
  FILE                    *fp_out_ps;

  char           *out_ancestral_file; /*! name of the file containing the ancestral sequences */
  FILE             *fp_out_ancestral; /*! pointer to the file containing the ancestral sequences */

  char                *in_coord_file; /*! name of input file containing coordinates */
  FILE                  *fp_in_coord; /*! pointer to the file containing coordinates */

  char                     *out_file; /*! name of the output file */

  char              *clade_list_file;

  int                       datatype; /*! 0->DNA, 1->AA */
  int               print_boot_trees; /*! =1 if the bootstrapped trees are printed in output */
  int       out_stats_file_open_mode; /*! opening file mode for statistics file */
  int        out_tree_file_open_mode; /*! opening file mode for tree file */
  int                    n_data_sets; /*! number of data sets to be analysed */
  int                        n_trees; /*! number of trees */
  int                       init_len; /*! sequence length */
  int                          n_otu; /*! number of taxa */
  int               n_data_set_asked; /*! number of bootstrap replicates */
  char                     *nt_or_cd; /*! nucleotide or codon data ? (not used) */
  int                      multigene; /*! if=1 -> analyse several partitions. */
  int               config_multigene;
  int                         n_part; /*! number of data partitions */
  int                        curr_gt;
  int                     ratio_test; /*! from 1 to 4 for specific branch supports, 0 of not */
  int                    ready_to_go;
  int                data_file_format; /*! Data format: Phylip or Nexus */
  int                tree_file_format; /*! Tree format: Phylip or Nexus */
  int                      state_len;

  int                 curr_interface;
  int                         r_seed; /*! random seed */
  int                  collapse_boot; /*! 0 -> branch length on bootstrap trees are not collapsed if too small */
  int          random_boot_seq_order; /*! !0 -> sequence order in bootstrapped data set is random */
  int                    print_trace;
  int                 print_site_lnl;
  int                       m4_model;
  int                      rm_ambigu; /*! 0 is the default. 1: columns with ambiguous characters are discarded prior further analysis */
  int                       colalias;
  int                  append_run_ID;
  char                *run_id_string;
  int                          quiet; /*! 0 is the default. 1: no interactive question (for batch mode) */
  int                      lk_approx; /* EXACT or NORMAL */
  char                    **alphabet;
  int                         codpos;
  int                         mutmap;
  int                        use_xml;

  char              **long_tax_names;
  char             **short_tax_names;
  int                 size_tax_names;

  phydbl                    *z_scores;
  phydbl                         *lat;
  phydbl                         *lon;

  int                 boot_prog_every;

  int                    mem_question;
  int                do_alias_subpatt;
  struct __Tmcmc                *mcmc;
  struct __T_Rate              *rates;
  
#ifdef BEAGLE
  int                 beagle_resource;
#endif

  int                       ancestral;
}option;

/*!********************************************************/

typedef struct __Optimiz { /*! parameters to be optimised (mostly used in 'optimiz.c') */
  int                    print; /*! =1 -> verbose mode  */

  int                opt_alpha; /*! =1 -> the gamma shape parameter is optimised */
  int                opt_kappa; /*! =1 -> the ts/tv ratio parameter is optimised */
  int               opt_lambda; /*! =1 -> the F84|TN93 model specific parameter is optimised */
  int               opt_pinvar; /*! =1 -> the proportion of invariants is optimised */
  int           opt_state_freq; /*! =1 -> the nucleotide frequencies are optimised */
  int                   opt_rr; /*! =1 -> the relative rate parameters of the GTR or the customn model are optimised */
  int          opt_subst_param; /*! if opt_topo=0 and opt_subst_param=1 -> the numerical parameters of the
                   model are optimised. if opt_topo=0 and opt_free_param=0 -> no parameter is
                   optimised */
  int            opt_cov_delta;
  int            opt_cov_alpha;
  int       opt_cov_free_rates;


  int                   opt_bl; /*! =1 -> the branch lengths are optimised */
  int                 opt_topo; /*! =1 -> the tree topology is optimised */
  int              topo_search;
  phydbl               init_lk; /*! initial loglikelihood value */
  int                 n_it_max; /*! maximum bnumber of iteration during an optimisation step */
  int                 last_opt; /*! =1 -> the numerical parameters are optimised further while the
                                       tree topology remains fixed */
  int        random_input_tree; /*! boolean */
  int            n_rand_starts; /*! number of random starting points */
  int             brent_it_max;
  int                steph_spr;
  int          user_state_freq;
  int          opt_five_branch;
  int              pars_thresh;
  int            hybrid_thresh;
  int          opt_br_len_mult;

  phydbl        tree_size_mult; /*! tree size multiplier */
  phydbl     min_diff_lk_local;
  phydbl    min_diff_lk_global;
  phydbl      min_diff_lk_move;
  phydbl    p_moves_to_examine;
  int                 fast_nni;
  int                   greedy;
  int             general_pars;
  int               quickdirty;
  int                 spr_pars;
  int                  spr_lnL;
  int           max_depth_path;
  int           min_depth_path;
  int             deepest_path;
  phydbl     max_delta_lnL_spr;
  int            br_len_in_spr;
  int      opt_free_mixt_rates;
  int       constrained_br_len;
  int         opt_gamma_br_len;
  int first_opt_free_mixt_rates;
  int              wim_n_rgrft;
  int              wim_n_globl;
  int             wim_max_dist;
  int              wim_n_optim;
  int               wim_n_best;
  int           wim_inside_opt;

  int          opt_rmat_weight;
  int          opt_efrq_weight;

  int      skip_tree_traversal;
  int        serial_free_rates;

  int      curr_opt_free_rates;
}t_opt;

/*!********************************************************/

typedef struct __NNI{

  struct __Node         *left;
  struct __Node         *rght;
  struct __Edge            *b;

  phydbl                score;
  phydbl               init_l;
  phydbl              init_lk;
  phydbl               best_l;
  phydbl          lk0,lk1,lk2;
  phydbl             l0,l1,l2;

  struct __Node *swap_node_v1;
  struct __Node *swap_node_v2;
  struct __Node *swap_node_v3;
  struct __Node *swap_node_v4;

  int       best_conf;   /*! best topological configuration :
                ((left_1,left_2),right_1,right_2) or
                ((left_1,right_2),right_1,left_2) or
                ((left_1,right_1),right_1,left_2)  */
}nni;

/*!********************************************************/

typedef struct __SPR{
  struct __Node         *n_link;
  struct __Node  *n_opp_to_link;
  struct __Edge  *b_opp_to_link;
  struct __Edge       *b_target;
  struct __Edge  *b_init_target;
  struct __Node          **path;
  phydbl          init_target_l;
  phydbl          init_target_v;
  phydbl               l0,l1,l2;
  phydbl               v0,v1,v2;
  phydbl                    lnL;
  int                depth_path;
  int                      pars;
  int                      dist;

  struct __SPR            *next;
  struct __SPR            *prev;
  struct __SPR       *next_mixt;
  struct __SPR       *prev_mixt;

}t_spr;

/*!********************************************************/

typedef struct __Triplet{
  int    size;
  phydbl *F_bc;
  phydbl *F_cd;
  phydbl *F_bd;
  phydbl ****core;
  phydbl ***p_one_site;
  phydbl ***sum_p_one_site;
  phydbl *pi_bc;
  phydbl *pi_cd;
  phydbl *pi_bd;
  struct __Eigen *eigen_struct;
  struct __Model *mod;

  struct __Triplet *next;
  struct __Triplet *prev;
}triplet;

/*!********************************************************/

typedef struct __Pnode{
  struct __Pnode **next;
  int weight;
  int num;
}pnode;

/*!********************************************************/

typedef struct __M4 {
  int                  n_h; /*! number of hidden states */
  int                  n_o; /*! number of observable states  */
  int        use_cov_alpha;
  int         use_cov_free;

  phydbl          **o_mats; /*! set of matrices of substitution rates across observable states */
  phydbl          *multipl; /*! vector of values that multiply each o_mats matrix */
  phydbl             *o_rr; /*! relative rates (symmetric) of substitution between observable states */
  phydbl             *h_rr; /*! relative rates (symmetric) of substitution between hidden states */
  phydbl            *h_mat; /*! matrix that describes the substitutions between hidden states (aka switches) */
  phydbl             *o_fq; /*! equilibrium frequencies for the observable states */
  phydbl             *h_fq; /*! equilibrium frequencies for the hidden states */
  phydbl    *h_fq_unscaled; /*! unscaled equilibrium frequencies for the hidden states */
  phydbl *multipl_unscaled; /*! unscaled  vector of values that multiply each o_mats matrix */

  phydbl             delta; /*! switching rate */
  phydbl             alpha; /*! gamma shape parameter */
}m4;

/*!********************************************************/

typedef struct __Tdraw {
  phydbl          *xcoord; /*! t_node coordinates on the x axis */
  phydbl          *ycoord; /*! t_node coordinates on the y axis */
  phydbl        *xcoord_s; /*! t_node coordinates on the x axis (scaled) */
  phydbl        *ycoord_s; /*! t_node coordinates on the y axis (scaled) */
  int          page_width;
  int         page_height;
  int      tree_box_width;

  int         *cdf_mat;
  phydbl       *cdf_mat_x;
  phydbl       *cdf_mat_y;


  phydbl max_dist_to_root;
}tdraw;

/*!********************************************************/

typedef struct __T_Rate {
  phydbl lexp; /*! Parameter of the exponential distribution that governs the rate at which substitution between rate classes ocur */
  phydbl alpha;
  phydbl less_likely;
  phydbl birth_rate;
  phydbl birth_rate_min;
  phydbl birth_rate_max;
  phydbl min_rate;
  phydbl max_rate;
  phydbl c_lnL1;
  phydbl c_lnL2;
  phydbl c_lnL_rates; /*! Prob(Br len | time stamps, model of rate evolution) */
  phydbl c_lnL_times; /*! Prob(time stamps) */
  phydbl c_lnL_jps; /*! Prob(# Jumps | time stamps, rates, model of rate evolution) */
  phydbl clock_r; /*! Mean substitution rate, i.e., 'molecular clock' rate */
  phydbl min_clock;
  phydbl max_clock;
  phydbl lbda_nu;
  phydbl min_dt;
  phydbl step_rate;
  phydbl true_tree_size;
  phydbl p_max;
  phydbl norm_fact;
  phydbl nu; /*! Parameter of the Exponential distribution for the corresponding model */
  phydbl min_nu;
  phydbl max_nu;
  phydbl covdet;
  phydbl sum_invalid_areas;


  phydbl     *nd_r;  /*! Current rates at nodes */
  phydbl     *br_r;  /*! Current rates along edges */
  phydbl     *nd_t; /*! Current t_node times */
  phydbl     *triplet;
  phydbl     *true_t; /*! true t_node times (including root node) */
  phydbl     *true_r; /*! true t_edge rates (on rooted tree) */
  phydbl     *buff_t;
  phydbl     *buff_r;
  phydbl     *dens; /*! Probability densities of mean substitution rates at the nodes */
  phydbl     *ml_l; /*! ML t_edge lengths (rooted) */
  phydbl     *cur_l; /*! Current t_edge lengths (rooted) */
  phydbl     *u_ml_l; /*! ML t_edge lengths (unrooted) */
  phydbl     *u_cur_l; /*! Current t_edge lengths (unrooted) */
  phydbl     *invcov;
  phydbl     *cov_r;
  phydbl     *mean_r; /*! average values of br_r taken across the sampled values during the MCMC */
  phydbl     *mean_t; /*! average values of nd_t taken across the sampled values during the MCMC */
  phydbl     *_2n_vect1;
  phydbl     *_2n_vect2;
  phydbl     *_2n_vect3;
  phydbl     *_2n_vect4;
  short int  *_2n_vect5;
  phydbl     *_2n2n_vect1;
  phydbl     *_2n2n_vect2;
  phydbl     *trip_cond_cov;
  phydbl     *trip_reg_coeff;
  phydbl     *cond_var;
  phydbl     *reg_coeff;
  phydbl     *t_prior;
  phydbl     *t_prior_min;
  phydbl     *t_prior_max;
  phydbl     *t_floor;
  phydbl     *t_mean;
  int        *t_ranked;
  phydbl     *mean_l;
  phydbl     *cov_l;
  phydbl     *grad_l; /* gradient */
  phydbl     inflate_var;
  phydbl     *time_slice_lims;
  phydbl     *survival_rank;
  phydbl     *survival_dur;
  phydbl     *calib_prob;

  int adjust_rates; /*! if = 1, branch rates are adjusted such that a modification of a given t_node time
               does not modify any branch lengths */
  int use_rates; /*! if = 0, branch lengths are expressed as differences between t_node times */
  int bl_from_rt; /*! if =1, branch lengths are obtained as the product of cur_r and t */
  int approx;
  int model; /*! Model number */
  int is_allocated;
  int met_within_gibbs;

  int update_mean_l;
  int update_cov_l;

  int *n_jps;
  int *t_jps;
  int n_time_slices;
  int *n_time_slice_spans;
  int *curr_slice;
  int *n_tips_below;

  short int *t_has_prior;
  struct __Node **lca; /*! 2-way table of common ancestral nodes for each pari of nodes */

  short int *br_do_updt;
  phydbl *cur_gamma_prior_mean;
  phydbl *cur_gamma_prior_var;

  int model_log_rates;

  short int nd_t_recorded;
  short int br_r_recorded;

  int *has_survived;

  struct __Calibration *calib;
  int tot_num_cal;
  int *curr_nd_for_cal;
  phydbl c_lnL_Hastings_ratio;

  phydbl     *t_prior_min_ori;
  phydbl     *t_prior_max_ori;
  phydbl     *times_partial_proba;

  phydbl log_K_cur;
  int cur_comb_numb;
  int *numb_calib_chosen;
  phydbl *node_height_dens_log_norm_const_update;
  int update_time_norm_const;

}t_rate;

/*!********************************************************/

typedef struct __Tmcmc {
  struct __Option *io;

  phydbl *tune_move;
  phydbl *move_weight;

  phydbl *acc_rate;
  int *acc_move;
  int *run_move;
  int *prev_acc_move;
  int *prev_run_move;
  int *num_move;
  int *move_type;
  char **move_name;

  int num_move_nd_r;
  int num_move_br_r;
  int num_move_nd_t;
  int num_move_nu;
  int num_move_clock_r;
  int num_move_tree_height;
  int num_move_subtree_height;
  int num_move_kappa;
  int num_move_tree_rates;
  int num_move_subtree_rates;
  int num_move_updown_nu_cr;
  int num_move_updown_t_cr;
  int num_move_updown_t_br;
  int num_move_ras;
  int num_move_cov_rates;
  int num_move_cov_switch;
  int num_move_birth_rate;
  int num_move_jump_calibration;
  int num_move_geo_lambda;
  int num_move_geo_sigma;
  int num_move_geo_tau;
  int num_move_geo_dum;
  int num_move_geo_updown_tau_lbda;
  int num_move_geo_updown_lbda_sigma;
  int num_move_phyrex_lbda;
  int num_move_phyrex_mu;
  int num_move_phyrex_rad;
  int num_move_phyrex_indel_disk;
  int num_move_phyrex_move_disk_ct;
  int num_move_phyrex_move_disk_ud;
  int num_move_phyrex_swap_disk;
  int num_move_phyrex_indel_hit;
  int num_move_phyrex_move_ldsk;
  int num_move_phyrex_spr;
  int num_move_phyrex_scale_times;
  int num_move_phyrex_ldscape_lim;
  int num_move_phyrex_sigsq;
  int num_move_phyrex_sim;
  int num_move_phyrex_traj;
  int num_move_phyrex_lbda_times;
  int num_move_phyrex_ldsk_given_disk;
  int num_move_phyrex_indel_disk_serial;

  int nd_t_digits;
  int *monitor;

  char *out_filename;

  time_t t_beg;
  time_t t_cur;
  time_t t_last_print;

  FILE *out_fp_stats;
  FILE *out_fp_trees;
  FILE *out_fp_means;
  FILE *out_fp_last;
  FILE *out_fp_constree;
  FILE *in_fp_par;

  int *adjust_tuning;
  int n_moves;
  int use_data;
  int randomize;
  int norm_freq;
  int run;
  int chain_len;
  int sample_interval;
  int chain_len_burnin;
  int print_every;
  int is_burnin;
  int max_lag;

  phydbl max_tune;
  phydbl min_tune;
  
  phydbl *sampled_val;
  int sample_size;
  int sample_num;
  phydbl *ess;
  int    *ess_run;
  int    *start_ess;
  phydbl *mode;
  int always_yes; /* Always accept proposed move (as long as log-likelihood > UNLIKELY) */
  int is; /* Importance sampling? Yes or NO */
}t_mcmc;

/*!********************************************************/

typedef struct __Tpart {
  int *ns;         /*! number of states for each partition (e.g., 2, 4, 3) */
  int *cum_ns;     /*! cumulative number of states (e.g., 0, 2, 6) */
  int ns_max;      /*! maximum number of states */
  int ns_min;      /*! minimum number of states */
  int n_partitions; /*! number of partitions */
  struct __Eigen    *eigen;
}part;

/*!********************************************************/

typedef struct __Tnexcom {
  char *name;
  int nparm;
  int nxt_token_t;
  int cur_token_t;
  struct __Tnexparm **parm;
}nexcom;

/*!********************************************************/

typedef struct __Tnexparm {
  char *name;
  char *value;
  int nxt_token_t;
  int cur_token_t;
  int (*fp)(char *, struct __Tnexparm *, struct __Option *);
  struct __Tnexcom *com;
}nexparm;

/*!********************************************************/

typedef struct __ParamInt {
  int val;
}t_param_int;

/*!********************************************************/

typedef struct __ParamDbl {
  phydbl val;
}t_param_dbl;

/*!********************************************************/

typedef struct __XML_node {

  struct __XML_attr *attr;   // Pointer to the first element of a list of attributes
  int n_attr;                // Number of attributes
  struct __XML_node      *next;   // Next sibling
  struct __XML_node      *prev;   // Previous sibling
  struct __XML_node    *parent; // Parent of this node
  struct __XML_node     *child;  // Child of this node
  char *id;
  char *name;
  char *value;
  struct __Generic_Data_Structure *ds; // Pointer to a data strucuture. Can be a scalar, a vector, anything.
}xml_node;

/*!********************************************************/

typedef struct __Generic_Data_Structure {
  void *obj;
  struct __Generic_Data_Structure *next;
}t_ds;


/*!********************************************************/

typedef struct __XML_attr {
  char *name;
  char *value;
  struct __XML_attr *next; // Next attribute
  struct __XML_attr *prev; // Previous attribute
}xml_attr;

/*!********************************************************/

typedef struct __Calibration {
  phydbl *proba; // Probability of this calibration (set by the user and fixed throughout)
  phydbl lower; // lower bound
  phydbl upper; // upper bound
  int cur_applies_to;
  phydbl calib_proba;
  struct __Node **all_applies_to;
  int n_all_applies_to;
  struct __Calibration *next;
  struct __Calibration *prev;
}t_cal;

/*!********************************************************/

typedef struct __Phylogeo{
  phydbl                    *cov; // Covariance of migrations (n_dim x n_dim)
  phydbl                  *r_mat; // R matrix. Gives the rates of migrations between locations. See article.
  phydbl                  *f_mat; // F matrix. See article.
  int                     *occup; // Vector giving the number of lineages that occupy each location
  int                   *idx_loc; // Index of location for each lineage
  int           *idx_loc_beneath; // Gives the index of location occupied beneath each node in the tree
  int                 ldscape_sz; // Landscape size: number of locations
  int                      n_dim; // Dimension of the data (e.g., longitude + lattitude -> n_dim = 2)
  int                update_fmat;
  struct __Geo_Coord **coord_loc; // Coordinates of the observed locations

  phydbl              sigma; // Dispersal parameter
  phydbl          min_sigma;
  phydbl          max_sigma;
  phydbl       sigma_thresh; // beyond sigma_thresh, there is no dispersal bias.

  phydbl               lbda; // Competition parameter
  phydbl           min_lbda;
  phydbl           max_lbda;

  phydbl              c_lnL;

  struct __Node **sorted_nd; // Table of nodes sorted wrt their heights.

  phydbl                tau; // overall migration rate parameter
  phydbl            min_tau;
  phydbl            max_tau;

  phydbl                dum; // dummy parameter use to assess non-identifiability issues
  phydbl            min_dum;
  phydbl            max_dum;

  
}t_geo;

/*!********************************************************/
// Structure for the Etheridge-Barton migration/reproduction model
typedef struct __Migrep_Model{
  int                           name;
  int                          n_dim;

  phydbl                        lbda; // rate at which events occur
  phydbl                    min_lbda; // min of rate at which events occur
  phydbl                    max_lbda; // max of rate at which events occur
  phydbl            prior_param_lbda; // parameter of the parameter for the prior on lbda

  phydbl                          mu; // per-capita and per event death probability
  phydbl                      min_mu; // min of per-capita and per event death probability
  phydbl                      max_mu; // max of per-capita and per event death probability
  phydbl              prior_param_mu; // parameter of the parameter for the prior on mu

  phydbl                         rad; // radius of the migrep disk 
  phydbl                     min_rad; // min of radius of the migrep disk 
  phydbl                     max_rad; // max of radius of the migrep disk 
  phydbl             prior_param_rad; // parameter of the parameter for the prior on radius
  int                     update_rad;

  phydbl                       sigsq; // parent to offspring distance variance (i.e., gene flow) parameter. 
  phydbl                   min_sigsq; // min 
  phydbl                   max_sigsq; // max  
  phydbl           prior_param_sigsq; // parameter of the parameter for the prior 

  phydbl                         rho; // intensity parameter of the Poisson point processs

  phydbl                       c_lnL; // current value of log-likelihood 
  phydbl              c_ln_prior_rad; // current value of log prior for the prior on radius
  phydbl             c_ln_prior_lbda; // current value of log prior for the prior on lbda
  phydbl               c_ln_prior_mu; // current value of log prior for the prior on mu
  phydbl            c_ln_prior_sigsq; // current value of log prior for the prior on sigsq=4.pi.lbda.mu.rad^4

  int                    safe_phyrex;
  phydbl             soft_bound_area;

  struct __Geo_Coord            *lim; // max longitude and lattitude (the min are both set to zero)                       

  phydbl                  sampl_area;
}t_phyrex_mod;

/*!********************************************************/

typedef struct __Disk_Event{
  struct __Geo_Coord           *centr;
  phydbl                         time;
  struct __Disk_Event           *next;
  struct __Disk_Event           *prev;
  struct __Lindisk_Node      **ldsk_a; // array of lindisk nodes corresponding to this disk event
  int                        n_ldsk_a; // size of ldsk_a
  struct __Lindisk_Node         *ldsk;
  struct __Migrep_Model         *mmod;
  char                            *id;
  phydbl                        c_lnL;
}t_dsk;

/*!********************************************************/

typedef struct __Geo_Coord{
  phydbl             *lonlat; /* longitude-latitude vector */
  int                    dim;
  char                   *id;
  struct __Geo_Coord    *cpy; /* keep a copy of this coordinate */

}t_geo_coord;

/*!********************************************************/

typedef struct __Lindisk_Node{
  struct __Disk_Event     *disk;
  struct __Lindisk_Node  **next;
  struct __Lindisk_Node   *prev;
  struct __Geo_Coord     *coord; 
  struct __Geo_Coord *cpy_coord; 
  short int              is_hit;
  int                    n_next;
  struct __Node             *nd;
}t_ldsk;

/*!********************************************************/

typedef struct __Polygon{
  struct __Geo_Coord **poly_vert; /* array of polygon vertex coordinates */
  int n_poly_vert; /* number of vertices */
}t_poly;

/*!********************************************************/
/*!********************************************************/
/*!********************************************************/
/*!********************************************************/
/*!********************************************************/
/*!********************************************************/





void Unroot_Tree(char **subtrees);
void Init_Tree(t_tree *tree,int n_otu);
void Init_Edge_Light(t_edge *b,int num);
void Init_Node_Light(t_node *n,int num);
void Set_Edge_Dirs(t_edge *b,t_node *a,t_node *d,t_tree *tree);
void Init_NNI(nni *a_nni);
void Init_Nexus_Format(nexcom **com);
void Restrict_To_Coding_Position(align **data,option *io);
void Uppercase(char *ch);
void Lowercase(char *ch);
calign *Compact_Data(align **data,option *io);
calign *Compact_Cdata(calign *data,option *io);
void Traverse_Prefix_Tree(int site,int seqnum,int *patt_num,int *n_patt,align **data,option *io,pnode *n);
pnode *Create_Pnode(int size);
void Get_Base_Freqs(calign *data);
void Get_AA_Freqs(calign *data);
void Swap_Nodes_On_Edges(t_edge *e1,t_edge *e2,int swap,t_tree *tree);
void Connect_Edges_To_Nodes_Recur(t_node *a,t_node *d,t_tree *tree);
void Connect_One_Edge_To_Two_Nodes(t_node *a,t_node *d,t_edge *b,t_tree *tree);
void Update_Dirs(t_tree *tree);
void Exit(char *message);
void *mCalloc(int nb,size_t size);
void *mRealloc(void *p,int nb,size_t size);
int Sort_Phydbl_Decrease(const void *a,const void *b);
void Qksort_Int(int *A,int *B,int ilo,int ihi);
void Qksort(phydbl *A,phydbl *B,int ilo,int ihi);
void Qksort_Matrix(phydbl **A,int col,int ilo,int ihi);
void Order_Tree_Seq(t_tree *tree,align **data);
char *Add_Taxa_To_Constraint_Tree(FILE *fp,calign *cdata);
void Check_Constraint_Tree_Taxa_Names(t_tree *tree,calign *cdata);
void Order_Tree_CSeq(t_tree *tree,calign *cdata);
void Init_Mat(matrix *mat,calign *data);
void Copy_Tax_Names_To_Tip_Labels(t_tree *tree,calign *data);
void Share_Lk_Struct(t_tree *t_full,t_tree *t_empt);
void Share_Spr_Struct(t_tree *t_full,t_tree *t_empt);
void Share_Pars_Struct(t_tree *t_full,t_tree *t_empt);
int Sort_Edges_NNI_Score(t_tree *tree,t_edge **sorted_edges,int n_elem);
int Sort_Edges_Depth(t_tree *tree,t_edge **sorted_edges,int n_elem);
void NNI(t_tree *tree,t_edge *b_fcus,int do_swap);
void NNI_Pars(t_tree *tree,t_edge *b_fcus,int do_swap);
void Swap(t_node *a,t_node *b,t_node *c,t_node *d,t_tree *tree);
void Update_All_Partial_Lk(t_edge *b_fcus,t_tree *tree);
void Update_SubTree_Partial_Lk(t_edge *b_fcus,t_node *a,t_node *d,t_tree *tree);
void Copy_Seq_Names_To_Tip_Labels(t_tree *tree,calign *data);
calign *Copy_Cseq(calign *ori,option *io);
int Filexists(char *filename);
matrix *K80_dist(calign *data,phydbl g_shape);
matrix *JC69_Dist(calign *data,t_mod *mod);
matrix *Hamming_Dist(calign *data,t_mod *mod);
int Is_Invar(int patt_num,int stepsize,int datatype,calign *data);
int Is_Ambigu(char *state,int datatype,int stepsize);
void Check_Ambiguities(calign *data,int datatype,int stepsize);
int Get_State_From_Ui(int ui,int datatype);
int Assign_State(char *c,int datatype,int stepsize);
char Reciproc_Assign_State(int i_state,int datatype);
int Assign_State_With_Ambiguity(char *c,int datatype,int stepsize);
void Clean_Tree_Connections(t_tree *tree);
void Bootstrap(t_tree *tree);
void Br_Len_Involving_Invar(t_tree *tree);
void Br_Len_Not_Involving_Invar(t_tree *tree);
void Getstring_Stdin(char *s);
phydbl Num_Derivatives_One_Param(phydbl (*func)(t_tree *tree), t_tree *tree,
                                 phydbl f0, phydbl *param, int which, int n_param, phydbl stepsize, int logt,
                                 phydbl *err, int precise, int is_positive);
phydbl Num_Derivatives_One_Param_Nonaligned(phydbl (*func)(t_tree *tree), t_tree *tree,
                                            phydbl f0, phydbl **param, int which, int n_param, phydbl stepsize, int logt,
                                            phydbl *err, int precise, int is_positive);
int Num_Derivative_Several_Param(t_tree *tree,phydbl *param,int n_param,phydbl stepsize,int logt,phydbl(*func)(t_tree *tree),phydbl *derivatives, int is_positive);
int Num_Derivative_Several_Param_Nonaligned(t_tree *tree, phydbl **param, int n_param, phydbl stepsize, int logt,
                                            phydbl (*func)(t_tree *tree), phydbl *derivatives, int is_positive);
int Compare_Two_States(char *state1,char *state2,int state_size);
void Copy_One_State(char *from,char *to,int state_size);
void Copy_Dist(phydbl **cpy,phydbl **orig,int n);
t_mod *Copy_Model(t_mod *ori);
void Record_Model(t_mod *ori,t_mod *cpy);
void Set_Defaults_Input(option *io);
void Set_Defaults_Model(t_mod *mod);
void Set_Defaults_Optimiz(t_opt *s_opt);
void Test_Node_Table_Consistency(t_tree *tree);
void Get_Bip(t_node *a,t_node *d,t_tree *tree);
void Alloc_Bip(t_tree *tree);
int Sort_Phydbl_Increase(const void *a,const void *b);
int Sort_String(const void *a,const void *b);
int Compare_Bip(t_tree *tree1,t_tree *tree2,int on_existing_edges_only);
void Match_Tip_Numbers(t_tree *tree1,t_tree *tree2);
void Test_Multiple_Data_Set_Format(option *io);
int Are_Compatible(char *statea,char *stateb,int stepsize,int datatype);
void Hide_Ambiguities(calign *data);
void Copy_Tree(t_tree *ori,t_tree *cpy);
void Prune_Subtree(t_node *a,t_node *d,t_edge **target,t_edge **residual,t_tree *tree);
void Graft_Subtree(t_edge *target,t_node *link,t_edge *residual,t_tree *tree);
void Reassign_Node_Nums(t_node *a,t_node *d,int *curr_ext_node,int *curr_int_node,t_tree *tree);
void Reassign_Edge_Nums(t_node *a,t_node *d,int *curr_br,t_tree *tree);
void Find_Mutual_Direction(t_node *n1,t_node *n2,short int *dir_n1_to_n2,short int *dir_n2_to_n1);
void Update_Dir_To_Tips(t_node *a,t_node *d,t_tree *tree);
void Fill_Dir_Table(t_tree *tree);
int Get_Subtree_Size(t_node *a,t_node *d);
void Fast_Br_Len(t_edge *b,t_tree *tree,int approx);
void Init_Eigen_Struct(eigen *this);
phydbl Triple_Dist(t_node *a,t_tree *tree,int approx);
void Make_Symmetric(phydbl **F,int size);
void Round_Down_Freq_Patt(phydbl **F,t_tree *tree);
phydbl Get_Sum_Of_Cells(phydbl *F,t_tree *tree);
void Divide_Cells(phydbl **F,phydbl div,t_tree *tree);
void Divide_Mat_By_Vect(phydbl **F,phydbl *vect,int size);
void Multiply_Mat_By_Vect(phydbl **F,phydbl *vect,int size);
void Found_In_Subtree(t_node *a,t_node *d,t_node *target,int *match,t_tree *tree);
void Get_List_Of_Target_Edges(t_node *a,t_node *d,t_edge **list,int *list_size,t_tree *tree);
void Fix_All(t_tree *tree);
void Record_Br_Len(t_tree *tree);
void Restore_Br_Len(t_tree *tree);
void Get_Dist_Btw_Edges(t_node *a,t_node *d,t_tree *tree);
void Detect_Polytomies(t_edge *b,phydbl l_thresh,t_tree *tree);
void Get_List_Of_Nodes_In_Polytomy(t_node *a,t_node *d,t_node ***list,int *size_list);
void Check_Path(t_node *a,t_node *d,t_node *target,t_tree *tree);
void Connect_Two_Nodes(t_node *a,t_node *d);
void Get_List_Of_Adjacent_Targets(t_node *a,t_node *d,t_node ***node_list,t_edge ***edge_list,int *list_size);
void Sort_List_Of_Adjacent_Targets(t_edge ***list,int list_size);
t_node *Common_Nodes_Btw_Two_Edges(t_edge *a,t_edge *b);
int KH_Test(phydbl *site_lk_M1,phydbl *site_lk_M2,t_tree *tree);
void Random_Tree(t_tree *tree);
void Reorganize_Edges_Given_Lk_Struct(t_tree *tree);
void Random_NNI(int n_moves,t_tree *tree);
void Fill_Missing_Dist(matrix *mat);
void Fill_Missing_Dist_XY(int x,int y,matrix *mat);
phydbl Least_Square_Missing_Dist_XY(int x,int y,phydbl dxy,matrix *mat);
void Check_Memory_Amount(t_tree *tree);
int Get_State_From_P_Lk(phydbl *p_lk,int pos,t_tree *tree);
int Get_State_From_P_Pars(short int *p_pars,int pos,t_tree *tree);
void Check_Dirs(t_tree *tree);
void Warn_And_Exit(const char *s);
void Randomize_Sequence_Order(calign *cdata);
void Update_Root_Pos(t_tree *tree);
void Add_Root(t_edge *target,t_tree *tree);
void Update_Ancestors(t_node *a,t_node *d,t_tree *tree);
#if (defined PHYTIME || defined SERGEII)
t_tree *Generate_Random_Tree_From_Scratch(int n_otu,int rooted);
#endif
void Random_Lineage_Rates(t_node *a,t_node *d,t_edge *b,phydbl stick_prob,phydbl *rates,int curr_rate,int n_rates,t_tree *tree);
t_edge *Find_Edge_With_Label(char *label,t_tree *tree);
void Evolve(calign *data,t_mod *mod,t_tree *tree);
int Pick_State(int n,phydbl *prob);
void Evolve_Recur(t_node *a,t_node *d,t_edge *b,int a_state,int r_class,int site_num,calign *gen_data,t_mod *mod,t_tree *tree);
void Site_Diversity(t_tree *tree);
void Site_Diversity_Post(t_node *a,t_node *d,t_edge *b,t_tree *tree);
void Site_Diversity_Pre(t_node *a,t_node *d,t_edge *b,t_tree *tree);
void Subtree_Union(t_node *n,t_edge *b_fcus,t_tree *tree);
void Binary_Decomposition(int value,int *bit_vect,int size);
void Print_Diversity_Header(FILE *fp,t_tree *tree);
phydbl Univariate_Kernel_Density_Estimate(phydbl where,phydbl *x,int sample_size);
phydbl Multivariate_Kernel_Density_Estimate(phydbl *where,phydbl **x,int sample_size,int vect_size);
phydbl Var(phydbl *x,int n);
phydbl Mean(phydbl *x,int n);
void Best_Of_NNI_And_SPR(t_tree *tree);
int Polint(phydbl *xa,phydbl *ya,int n,phydbl x,phydbl *y,phydbl *dy);
void JF(t_tree *tree);
t_tree *Dist_And_BioNJ(calign *cdata,t_mod *mod,option *io);
void Add_BioNJ_Branch_Lengths(t_tree *tree,calign *cdata,t_mod *mod);
char *Bootstrap_From_String(char *s_tree,calign *cdata,t_mod *mod,option *io);
char *aLRT_From_String(char *s_tree,calign *cdata,t_mod *mod,option *io);
void Prepare_Tree_For_Lk(t_tree *tree);
void Find_Common_Tips(t_tree *tree1,t_tree *tree2);
phydbl Get_Tree_Size(t_tree *tree);
void Dist_To_Root_Pre(t_node *a,t_node *d,t_edge *b,t_tree *tree);
void Dist_To_Root(t_node *n_root,t_tree *tree);
char *Basename(char *path);
t_node *Find_Lca_Pair_Of_Nodes(t_node *n1,t_node *n2,t_tree *tree);
t_node *Find_Lca_Clade(t_node **node_list,int node_list_size,t_tree *tree);
int Get_List_Of_Ancestors(t_node *ref_node,t_node **list,int *size,t_tree *tree);
int Edge_Num_To_Node_Num(int edge_num,t_tree *tree);
void Time_To_Branch(t_tree *tree);
void Time_To_Branch_Pre(t_node *a,t_node *d,t_tree *tree);
void Branch_Lengths_To_Rate_Lengths(t_tree *tree);
void Branch_Lengths_To_Rate_Lengths_Pre(t_node *a,t_node *d,t_tree *tree);
int Find_Clade(char **tax_name_list,int list_size,t_tree *tree);
void Find_Clade_Pre(t_node *a,t_node *d,int *tax_num_list,int list_size,int *num,t_tree *tree);
t_edge *Find_Root_Edge(FILE *fp_input_tree,t_tree *tree);
void Copy_Tree_Topology_With_Labels(t_tree *ori,t_tree *cpy);
void Set_Model_Name(t_mod *mod);
void Adjust_Min_Diff_Lk(t_tree *tree);
void Translate_Tax_Names(char **tax_names,t_tree *tree);
void Skip_Comment(FILE *fp);
void Get_Best_Root_Position(t_tree *tree);
void Get_Best_Root_Position_Post(t_node *a,t_node *d,int *has_outgrp,t_tree *tree);
void Get_Best_Root_Position_Pre(t_node *a,t_node *d,t_tree *tree);
void Get_OutIn_Scores(t_node *a,t_node *d);
int Check_Sequence_Name(char *s);
int Scale_Subtree_Height(t_node *a,phydbl K,phydbl floor,int *n_nodes,t_tree *tree);
void Scale_Node_Heights_Post(t_node *a,t_node *d,phydbl K,phydbl floor,int *n_nodes,t_tree *tree);
int Scale_Subtree_Rates(t_node *a,phydbl mult,int *n_nodes,t_tree *tree);
void Check_Br_Len_Bounds(t_tree *tree);
int Scale_Subtree_Rates_Post(t_node *a,t_node *d,phydbl mult,int *n_nodes,t_tree *tree);
void Get_Node_Ranks(t_tree *tree);
void Get_Node_Ranks_Pre(t_node *a,t_node *d,t_tree *tree);
void Log_Br_Len(t_tree *tree);
phydbl Diff_Lk_Norm_At_Given_Edge(t_edge *b,t_tree *tree);
void Adjust_Variances(t_tree *tree);
phydbl Effective_Sample_Size(phydbl first_val,phydbl last_val,phydbl sum,phydbl sumsq,phydbl sumcurnext,int n);
void Rescale_Free_Rate_Tree(t_tree *tree);
phydbl Rescale_Br_Len_Multiplier_Tree(t_tree *tree);
phydbl Unscale_Br_Len_Multiplier_Tree(t_tree *tree);
phydbl Reflect(phydbl x,phydbl l,phydbl u);
int Are_Equal(phydbl a,phydbl b,phydbl eps);
int Check_Topo_Constraints(t_tree *big_tree,t_tree *small_tree);
void Prune_Tree(t_tree *big_tree,t_tree *small_tree);
void Match_Nodes_In_Small_Tree(t_tree *small_tree,t_tree *big_tree);
void Find_Surviving_Edges_In_Small_Tree(t_tree *small_tree,t_tree *big_tree);
void Find_Surviving_Edges_In_Small_Tree_Post(t_node *a,t_node *d,t_tree *small_tree,t_tree *big_tree);
void Set_Taxa_Id_Ranking(t_tree *tree);
void Get_Edge_Binary_Coding_Number(t_tree *tree);
void Make_Ancestral_Seq(t_tree *tree);
void Make_MutMap(t_tree *tree);
int Get_Mutmap_Val(int edge,int site,int mut,t_tree *tree);
void Get_Mutmap_Coord(int idx,int *edge,int *site,int *mut,t_tree *tree);
void Copy_Edge_Lengths(t_tree *to,t_tree *from);
void Init_Scalar_Dbl(scalar_dbl *p);
void Init_Scalar_Int(scalar_int *p);
void Init_Vect_Dbl(int len,vect_dbl *p);
void Init_Vect_Int(int len,vect_int *p);
char *To_Lower_String(char *in);
phydbl String_To_Dbl(char *string);
char *To_Upper_String(char *in);
void Connect_CSeqs_To_Nodes(calign *cdata, option *io, t_tree *tree);
void Switch_Eigen(int state, t_mod *mod);
void Joint_Proba_States_Left_Right(phydbl *Pij, phydbl *p_lk_left, phydbl *p_lk_rght,
                   vect_dbl *pi, int scale_left, int scale_rght,
                   phydbl *F, int n, int site, t_tree *tree);
void Set_Both_Sides(int yesno, t_tree *tree);
void Set_D_States(calign *data, int datatype, int stepsize);
void Branch_To_Time(t_tree *tree);
void Branch_To_Time_Pre(t_node *a, t_node *d, t_tree *tree);
void Path_Length(t_node *dep, t_node *arr, phydbl *len, t_tree *tree);
phydbl *Dist_Btw_Tips(t_tree *tree);
void Random_SPRs_On_Rooted_Tree(t_tree *tree);


void Set_P_Lk_One_Side(phydbl **Pij, phydbl **p_lk,  int **sum_scale, t_node *d, t_edge *b, t_tree *tree
#ifdef BEAGLE
                       , int* child_p_idx, int* Pij_idx
#endif
                       );


void Set_All_P_Lk(t_node **n_v1, t_node **n_v2,
                                 phydbl **p_lk , int **sum_scale , int **p_lk_loc,
                  phydbl **Pij1, phydbl **p_lk1, int **sum_scale1,
                  phydbl **Pij2, phydbl **p_lk2, int **sum_scale2,
                  t_node *d, t_edge *b, t_tree *tree
#ifdef BEAGLE
                  , int *dest_p_idx, int *child1_p_idx, int* child2_p_idx, int* Pij1_idx, int* Pij2_idx
#endif
                  );

void Optimum_Root_Position_IL_Model(t_tree *tree);
void Set_Br_Len_Var(t_tree *tree);
void Check_Br_Lens(t_tree *tree);
void Calculate_Number_Of_Diff_States_Post(t_node *a, t_node *d, t_edge *b, t_tree *tree);
void Calculate_Number_Of_Diff_States_Pre(t_node *a, t_node *d, t_edge *b, t_tree *tree);
void Calculate_Number_Of_Diff_States_Core(t_node *a, t_node *d, t_edge *b, t_tree *tree);
void Calculate_Number_Of_Diff_States(t_tree *tree);
void Build_Distrib_Number_Of_Diff_States_Under_Model(t_tree *tree);
int Number_Of_Diff_States_One_Site(int site, t_tree *tree);
void Number_Of_Diff_States_One_Site_Post(t_node *a, t_node *d, t_edge *b, int site, t_tree *tree);
int Number_Of_Diff_States_One_Site_Core(t_node *a, t_node *d, t_edge *b, int site, t_tree *tree);
phydbl Get_Lk(t_tree *tree);
align **Make_Empty_Alignment(option *io);
void Connect_Edges_To_Nodes_Serial(t_tree *tree);
phydbl Mean_Identity(calign *data);
phydbl Pairwise_Identity(int i, int j, calign *data);
phydbl Fst(int i, int j, calign *data);
phydbl Nucleotide_Diversity(calign *data);


#include "xml.h"
#include "free.h"
#include "spr.h"
#include "lk.h"
#include "optimiz.h"
#include "models.h"
#include "bionj.h"
#include "simu.h"
#include "eigen.h"
#include "pars.h"
#include "alrt.h"
#include "stats.h"
#include "help.h"
#include "io.h"
#include "make.h"
#include "nexus.h"
#include "init.h"
#include "mcmc.h"

#ifdef GEO
#include "geo.h"
#endif

#ifdef PHYREX
#include "phyrex.h"
#endif

#ifdef MPI
#include "mpi_boot.h"
#endif

#ifdef MG
#include "mg.h"
#endif

#ifdef TIME
#include "times.h"
#include "rates.h"
#endif

#ifdef _NOT_NEEDED_A_PRIORI
#include "m4.h"
#endif


#endif
