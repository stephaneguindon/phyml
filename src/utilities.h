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

// #define _POSIX_C_SOURCE 200112L

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <execinfo.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <signal.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
/* #include <malloc/malloc.h> */
/* #include <malloc.h> */

#if (defined(__AVX__) || defined(__AVX2__))
#include <immintrin.h>
#include <pmmintrin.h>
#include <xmmintrin.h>
#elif defined(__ARM_NEON)
#include "sse2neon.h"
#elif (defined(__SSE__) || defined(__SSE2__) || defined(__SSE3__))
#include <emmintrin.h>
#include <pmmintrin.h>
#endif

extern int n_sec1;
extern int n_sec2;

#define TREE_COMP_RF_PLUS_LENGTH 0
#define TREE_COMP_RF_PLAIN 1
#define TREE_COMP_LENGTH_SHARED_EDGES 2
#define TREE_COMP_LENGTH_EXTERNAL_EDGES 3

#define HAVERSINE 0
#define MANHATTAN 1
#define EUCLIDEAN 2

#define __FUNCTION__ NULL

#define EMPIRICAL 0
#define ML 1
#define USER 2
#define MODEL 3

#define For(i, n) for (i = 0; i < n; i++)
#define Fors(i, n, s) for (i = 0; i < n; i += s)
#define PointGamma(prob, alpha, beta)                                          \
  PointChi2(prob, 2.0 * (alpha)) / (2.0 * (beta))
#define SHFT2(a, b, c)                                                         \
  (a) = (b);                                                                   \
  (b) = (c);
#define SHFT3(a, b, c, d)                                                      \
  (a) = (b);                                                                   \
  (b) = (c);                                                                   \
  (c) = (d);
#define SIGN(a, b) ((b) > 0.0 ? fabs(a) : -fabs(a))
#define SHFT(a, b, c, d)                                                       \
  (a) = (b);                                                                   \
  (b) = (c);                                                                   \
  (c) = (d);

#ifndef MAX
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#endif

#ifndef MIN
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#endif

#define READ 0
#define WRITE 1
#define APPEND 2
#define READWRITE 3

#ifndef isnan
#define isnan(x)                                                               \
  (sizeof(x) == sizeof(long double) ? isnan_ld(x)                              \
   : sizeof(x) == sizeof(double)    ? isnan_d(x)                               \
                                    : isnan_f(x))
static inline int isnan_f(float x) { return x != x; }
static inline int isnan_d(double x) { return x != x; }
static inline int isnan_ld(long double x) { return x != x; }
#endif

#ifndef isinf
#define isinf(x)                                                               \
  (sizeof(x) == sizeof(long double) ? isinf_ld(x)                              \
   : sizeof(x) == sizeof(double)    ? isinf_d(x)                               \
                                    : isinf_f(x))
static inline int isinf_f(float x) { return isnan(x - x); }
static inline int isinf_d(double x) { return isnan(x - x); }
static inline int isinf_ld(long double x) { return isnan(x - x); }
#endif

#define LOG12 2.4849066497880003546
#define LOG4 1.3862943611198905725
#define LOG3 1.0986122886681097821

#define BIRTHDEATH 0
#define COALESCENT 1

#define STRICTCOALESCENT 0
#define EXPCOALESCENT 1
#define POWLAW 2

#define SLFV_GAUSSIAN 2 /* Spatial Lambda-Fleming-Viot model (Gaussian) */
#define SLFV_UNIFORM 3  /* Spatial Lambda-Fleming-Viot model (Uniform) */
#define RW 4 /* standard Brownian diffusion model in phylogeography */
#define RRW_GAMMA                                                              \
  5 /* Lemey's relaxed random walk (Gamma distributed relative diffusion       \
       rates) */
#define RRW_LOGNORMAL                                                          \
  6 /* Lemey's relaxed random walk (Lognormal distributed relative diffusion   \
       rates) */
#define IBM 7    /* Integrated Brownian Motion */
#define RIBM 8   /* Relaxed Integrated Brownian Motion */
#define IWNc 9   /* Integrated white noise */
#define RIWNc 10 /* Relaxed integrated white noise */
#define IWNu 11  /* Integrated white noise uncorrelated */
#define RIWNu 12 /* Relaxed integrated white noise uncorrelated */
#define IOU 13   /* Integrated Ornstein-Uhlenbeck */

#define AC 0
#define AG 1
#define AT 2
#define CG 3
#define CT 4
#define GT 5

#if (defined __AVX__ || defined __AVX2__)
#define BYTE_ALIGN 32
#elif (defined(__SSE__) || defined(__SSE2__) || defined(__SSE3__) ||           \
       defined(__ARM_NEON))
#define BYTE_ALIGN 16
#else
#define BYTE_ALIGN 1
#endif

#ifndef M_1_SQRT_2PI
#define M_1_SQRT_2PI 0.398942280401432677939946059934 /* 1/sqrt(2*pi) */
#endif

#ifndef LOG_SQRT_2PI
#define LOG_SQRT_2PI 0.91893853320467266954 /* log(sqrt(2*pi)) */
#endif

// verbose levels
#define VL0 0
#define VL1 1
#define VL2 2
#define VL3 3
#define VL4 4

#define T_MAX_MCMC_MOVE_NAME 500

#define WINDOW_WIDTH 800
#define WINDOW_HEIGHT 800

#define ALPHA_MIN 0.01
#define ALPHA_MAX 1000.

#define E_FRQ_WEIGHT_MIN 0.01
#define E_FRQ_WEIGHT_MAX 100.

#define R_MAT_WEIGHT_MIN 0.01
#define R_MAT_WEIGHT_MAX 100.

#define E_FRQ_MIN 0.001
#define E_FRQ_MAX 0.999

#define UNSCALED_E_FRQ_MIN -100.
#define UNSCALED_E_FRQ_MAX +100.

#define TSTV_MIN 0.05
#define TSTV_MAX 100.0

#define PINV_MIN 0.00001
#define PINV_MAX 0.99999

#define RR_MIN 0.0001
#define RR_MAX 10000.0

#define UNSCALED_RR_MIN log(RR_MIN)
#define UNSCALED_RR_MAX log(RR_MAX)

#define LOCATION 0
#define VELOCITY 1

#ifdef PHYML
#define GAMMA_RR_UNSCALED_MIN -1000.
#define GAMMA_RR_UNSCALED_MAX 1000.
/* #define GAMMA_RR_UNSCALED_MIN 0.01 */
/* #define GAMMA_RR_UNSCALED_MAX 100. */
#else
#define GAMMA_RR_UNSCALED_MIN 0.01
#define GAMMA_RR_UNSCALED_MAX 200.
#endif

#ifdef PHYML
#define GAMMA_R_PROBA_UNSCALED_MIN -5.
#define GAMMA_R_PROBA_UNSCALED_MAX 5.
#else
#define GAMMA_R_PROBA_UNSCALED_MIN 0.01
#define GAMMA_R_PROBA_UNSCALED_MAX 200.
#endif

#define PHYREX_UNIFORM 0
#define PHYREX_NORMAL 1

#define MCMC_MOVE_RANDWALK_UNIFORM 0
#define MCMC_MOVE_LOG_RANDWALK_UNIFORM 1
#define MCMC_MOVE_RANDWALK_NORMAL 2
#define MCMC_MOVE_LOG_RANDWALK_NORMAL 3
#define MCMC_MOVE_SCALE_THORNE 4
#define MCMC_MOVE_SCALE_GAMMA 5

#define N_MAX_MOVES 50

#define N_MAX_NEX_COM 20
#define T_MAX_NEX_COM 100
#define N_MAX_NEX_PARM 50
#define T_MAX_TOKEN 200
#define T_MAX_ID_COORD 10
#define T_MAX_ID_DISK 10
#define T_MAX_KEY 200
#define T_MAX_VAL 1000

#define N_MAX_MIXT_CLASSES 1000

#define NEXUS_COM 0
#define NEXUS_PARM 1
#define NEXUS_EQUAL 2
#define NEXUS_VALUE 3
#define NEXUS_SPACE 4

#define NNI_MOVE 0
#define SPR_MOVE 1
#define BEST_OF_NNI_AND_SPR 2

#define M_1_SQRT_2_PI 0.398942280401432677939946059934
#define M_SQRT_32 5.656854249492380195206754896838
#define PI 3.14159265358979311600
#define SQRT2PI 2.50662827463100024161
#define LOG2PI 1.83787706640934533908
#define LOG2 0.69314718055994528623
#define LOG_SQRT_2_PI 0.918938533204672741780329736406

#define NORMAL 1
#define EXACT 2

#define PHYLIP 0
#define NEXUS 1
#define IBDSIM 2

#ifndef YES
#define YES 1
#endif

#ifndef NO
#define NO 0
#endif

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#define SPATIAL_SAMPLING_DETECTION 0
#define SPATIAL_SAMPLING_SURVEY 1

#define ON 1
#define OFF 0

#define SMALL_DBL 1.E-20

#define NT 0      /*! nucleotides */
#define AA 1      /*! amino acids */
#define GENERIC 2 /*! custom alphabet */
#define UNDEFINED -1

#define ACGT 0 /*! A,G,G,T encoding */
#define RY 1   /*! R,Y     encoding */

#define INTERFACE_DATA_TYPE 0
#define INTERFACE_MULTIGENE 1
#define INTERFACE_MODEL 2
#define INTERFACE_TOPO_SEARCH 3
#define INTERFACE_BRANCH_SUPPORT 4

#define LEFT 0
#define RGHT 1

#ifndef INFINITY
#define INFINITY HUGE
#endif

#define MAX_N_CAL 100
#define N_MAX_OPTIONS 100

#define NEXT_BLOCK_SIZE 50

#define T_MAX_FILE 200
#define T_MAX_LINE 2000000
#define T_MAX_NAME 1000
#define T_MAX_ID 20
#define T_MAX_SEQ 2000000
#define T_MAX_OPTION 100
#define T_MAX_STATE 5
#define S_TREE_CHUNK 10000

#define NODE_DEG_MAX 2000
#define BRENT_IT_MAX 1000
#define BRENT_CGOLD 0.3819660
#define BRENT_ZEPS 1.e-10
#define MNBRAK_GOLD 1.618034
#define MNBRAK_GLIMIT 100.0
#define MNBRAK_TINY 1.e-20
#define BL_START 1.e-04
#define GOLDEN_R 0.61803399
#define GOLDEN_C (1.0 - GOLDEN_R)
#define N_MAX_INSERT 20
#define N_MAX_OTU 4000
#define UNLIKELY -1.e20
#define NJ_SEUIL 0.1
#define ROUND_MAX 100
#define DIST_MAX 2.00
#define DIST_MIN 1.e-10
#define AROUND_LK 50.0
#define PROP_STEP 1.0
#define T_MAX_ALPHABET 22
#define MDBL_MIN FLT_MIN
#define MDBL_MAX FLT_MAX
#define POWELL_ITMAX 200
#define LINMIN_TOL 2.0E-04
#define SCALE_POW                                                              \
  10 /*! Scaling factor will be 2^SCALE_POW or 2^(-SCALE_POW) [[ WARNING:      \
        SCALE_POW < 31 ]]*/
#define DEFAULT_SIZE_SPR_LIST 20
#define NEWICK 0
#define OUTPUT_TREE_FORMAT NEWICK
#define MAX_PARS 1000000000

#define LIM_SCALE_VAL 1.E-50 /*! Scaling limit (deprecated) */

#define MIN_CLOCK_RATE 1.E-10

#define MIN_VAR_BL 1.E-8
#define MAX_VAR_BL 1.E+3

#define KFOLD_POS 0
#define KFOLD_COL 1
#define MAXFOLD 2

#define MASK_TYPE_POSITION                                                     \
  0 // Individual positions (i.e., site x taxon) are masked
#define MASK_TYPE_COLUMN 1 // Columns (i.e., site x .) are masked

#define JC69 1
#define K80 2
#define F81 3
#define HKY85 4
#define F84 5
#define TN93 6
#define GTR 7
#define CUSTOM 8

#define WAG 11
#define DAYHOFF 12
#define JTT 13
#define BLOSUM62 14
#define MTREV 15
#define RTREV 16
#define CPREV 17
#define DCMUT 18
#define VT 19
#define MTMAM 20
#define MTART 21
#define HIVW 22
#define HIVB 23
#define FLU 24
#define CUSTOMAA 25
#define LG 26
#define AB 27
// Amino acid ordering:
// Ala Arg Asn Asp Cys Gln Glu Gly His Ile Leu Lys Met Phe Pro Ser Thr Trp Tyr
// Val
#define OTHER 30

#define EXPONENTIAL_PRIOR 0
#define NORMAL_PRIOR 1
#define FLAT_PRIOR 2

#define COMPOUND_COR 0
#define COMPOUND_NOCOR 1
#define EXPONENTIAL 2
#define LOGNORMAL 3
#define THORNE 4
#define GUINDON 5
#define STRICTCLOCK 6
#define GAMMA 7
#define NONE -1

#define ALRTSTAT 1
#define ALRTCHI2 2
#define MINALRTCHI2SH 3
#define SH 4
#define ABAYES 5

#define CALENDAR 0
#define SUBSTITUTIONS 1

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
typedef double phydbl;
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
#define BIG DBL_MAX
#define SMALL_PIJ 1.E-100
#define LOGBIG 690.
#define LOGSMALL -690.

#if !(defined PHYTIME)
#define BL_MIN 1.E-8
#define BL_MAX 100.
#else
#define BL_MIN 1.E-6
#define BL_MAX 1.
#endif

// Do *not* change the values below and leave the lines with
// curr_scaler_pow = (int)(-XXX.-LOG(smallest_p_lk))/LOG2;
// as XXX depends on what the value of P_LK_LIM_INF is
#define P_LK_LIM_INF 3.054936e-151 /* 2^-500 */
#define P_LK_LIM_SUP 3.273391e+150 /* 2^500 */
/* #define  P_LK_LIM_INF   1.499697e-241 /\* 2^-800 *\/ */
/* #define  P_LK_LIM_SUP   6.668014e+240 /\* 2^800 *\/ */

/* #define LARGE 100 */
/* #define TWO_TO_THE_LARGE 1267650600228229401496703205376.0 /\* 2^100 In R :
 * sprintf("%.100f", 2^100)*\/ */

/* #define LARGE 200 */
/* #define TWO_TO_THE_LARGE
 * 1606938044258990275541962092341162602522202993782792835301376.0 /\* 2^200 In
 * R : sprintf("%.100f", 2^256)*\/ */

#define LARGE 256
#define TWO_TO_THE_LARGE                                                       \
  115792089237316195423570985008687907853269984665640564039457584007913129639936.0 /* 2^256 In R : sprintf("%.100f", 2^256)*/

/* #define LARGE 128 */
/* #define TWO_TO_THE_LARGE 340282366920938463463374607431768211456.0 /\* 2^128
 * In R : sprintf("%.100f", 2^128)*\/ */

/* #define LARGE 320 */
/* #define TWO_TO_THE_LARGE
 * 2135987035920910082395021706169552114602704522356652769947041607822219725780640550022962086936576.0
 */

#define INV_TWO_TO_THE_LARGE (1. / TWO_TO_THE_LARGE)

#define SCALE_RATE_SPECIFIC 1
#define SCALE_FAST 2

#define SCALENO 0
#define SCALEYES 1

#define T_MAX_XML_TAG 64

#define NARGS_SEQ(_1, _2, _3, _4, _5, _6, _7, _8, N, ...) N
#define NARGS(...) NARGS_SEQ(__VA_ARGS__, 8, 7, 6, 5, 4, 3, 2, 1)
#define PRIMITIVE_CAT(x, y) x##y
#define CAT(x, y) PRIMITIVE_CAT(x, y)
#define FOR_EACH(macro, ...)                                                   \
  CAT(FOR_EACH_, NARGS(__VA_ARGS__))(macro, __VA_ARGS__)
#define FOR_EACH_1(m, x1) m(x1)
#define FOR_EACH_2(m, x1, x2) m(x1) m(x2)
#define FOR_EACH_3(m, x1, x2, x3) m(x1) m(x2) m(x3)
#define FOR_EACH_4(m, x1, x2, x3, x4) m(x1) m(x2) m(x3) m(x4)
#define FOR_EACH_5(m, x1, x2, x3, x4, x5) m(x1) m(x2) m(x3) m(x4) m(x5)
#define FOR_EACH_6(m, x1, x2, x3, x4, x5, x6)                                  \
  m(x1) m(x2) m(x3) m(x4) m(x5) m(x6)
#define FOR_EACH_7(m, x1, x2, x3, x4, x5, x6, x7)                              \
  m(x1) m(x2) m(x3) m(x4) m(x5) m(x6) m(x7)
#define FOR_EACH_8(m, x1, x2, x3, x4, x5, x6, x7, x8)                          \
  m(x1) m(x2) m(x3) m(x4) m(x5) m(x6) m(x7) m(x8)
#define DUMP_EACH_INT(v)                                                       \
  fprintf(stderr, "\n\t\tDEBUG:%s:\t\t%s--->%i", __PRETTY_FUNCTION__, #v,      \
          (v));                                                                \
  fflush(stderr);
#define DUMP_EACH_STRING(v)                                                    \
  fprintf(stderr, "\n\t\tDEBUG:%s:\t\t%s--->%s", __PRETTY_FUNCTION__, #v,      \
          (v));                                                                \
  fflush(stderr);
#define DUMP_EACH_DECIMAL(v)                                                   \
  fprintf(stderr, "\n\t\tDEBUG:%s:\t\t%s--->%f", __PRETTY_FUNCTION__, #v,      \
          (v));                                                                \
  fflush(stderr);
#define DUMP_EACH_SCI(v)                                                       \
  fprintf(stderr, "\n\t\tDEBUG:%s:\t\t%s--->%e", __PRETTY_FUNCTION__, #v,      \
          (v));                                                                \
  fflush(stderr);
#define DUMP_I(...) FOR_EACH(DUMP_EACH_INT, __VA_ARGS__)
#define DUMP_S(...) FOR_EACH(DUMP_EACH_STRING, __VA_ARGS__)
#define DUMP_D(...) FOR_EACH(DUMP_EACH_DECIMAL, __VA_ARGS__)
#define DUMP_E(...) FOR_EACH(DUMP_EACH_SCI, __VA_ARGS__)

/*!********************************************************/
// Generic linked list
typedef struct __Generic_LL
{
  void                *v;
  struct __Generic_LL *next;
  struct __Generic_LL *prev;
  struct __Generic_LL *tail;
  struct __Generic_LL *head;
} t_ll;

/*!********************************************************/

typedef struct __Scalar_Int
{
  int                  v;
  bool                 optimize;
  struct __Scalar_Int *next;
  struct __Scalar_Int *prev;
} scalar_int;

/*!********************************************************/

typedef struct __Scalar_Dbl
{
  phydbl               v;
  bool                 onoff;
  bool                 optimize;
  bool                 print;
  struct __Scalar_Dbl *next;
  struct __Scalar_Dbl *prev;
} scalar_dbl;

/*!********************************************************/

typedef struct __Vect_Int
{
  int               *v;
  int                len;
  bool               optimize;
  struct __Vect_Int *next;
  struct __Vect_Int *prev;
} vect_int;

/*!********************************************************/

typedef struct __Vect_Dbl
{
  phydbl            *v;
  int                len;
  bool               optimize;
  struct __Vect_Dbl *next;
  struct __Vect_Dbl *prev;
} vect_dbl;

/*!********************************************************/

typedef struct __String
{
  char            *s;
  int              len;
  struct __String *next;
  struct __String *prev;
} t_string;

/*!********************************************************/

typedef struct __Node
{
  struct __Node  **v; /*! table of pointers to neighbor nodes. Dimension = 3 */
  struct __Node ***bip_node; /*! three lists of pointer to tip nodes. One list
                                for each direction */
  struct __Edge **b;         /*! table of pointers to neighbor branches */
  struct __Node  *anc;   /*! direct ancestor t_node (for rooted tree only) */
  struct __Edge  *b_anc; /*! edge between this node and its direct ancestor (for
                            rooted tree only) */
  struct __Node  *ext_node;
  struct __Node  *match_node;
  struct __Align *c_seq;     /*! corresponding compressed sequence */
  struct __Align *c_seq_anc; /*! corresponding compressed ancestral sequence */

  struct __Node *next; /*! tree->a_nodes[i]->next <=> tree->next->a_nodes[i] */
  struct __Node *prev; /*! See above */
  struct __Node *next_mixt; /*! Next mixture tree*/
  struct __Node *prev_mixt; /*! Parent mixture tree */
  struct __Calibration *
      *cal; /*! List of calibration constraints attached to this node */
  struct __Lindisk_Node
      *ldsk; /*! Used in PhyREX. Lineage/Disk this node corresponds to */
  struct __Node *rk_next; /*! Next node in the list of ranked nodes (from oldest
                             to youngest) */
  struct __Node *rk_prev; /*! Previous node in the list of ranked nodes (from
                             oldest to youngest) */
  struct __Label *label;

  int *bip_size; /*! Size of each of the three lists from bip_node */
  int  num;      /*! t_node number */
  int  tax;      /*! tax = 1 -> external node, else -> internal t_node */

  char *name;     /*! taxon name (if exists) */
  char *ori_name; /*! taxon name (if exists) */
  int   n_cal;    /*! Number of calibration constraints */

  phydbl *score; /*! score used in BioNJ to determine the best pair of nodes to
                    agglomerate */
  phydbl dist_to_root; /*! distance to the root t_node */

  short int common;
  phydbl    y_rank;
  phydbl    y_rank_ori;
  phydbl    y_rank_min;
  phydbl    y_rank_max;

  int *s_ingrp;  /*! does the subtree beneath belong to the ingroup */
  int *s_outgrp; /*! does the subtree beneath belong to the outgroup */

  int id_rank; /*! order taxa alphabetically and use id_rank to store the ranks
                */
  int rank;
  int rank_max;

  struct __Geo_Coord *coord;
  struct __Geo_Veloc *veloc;

} t_node;

/*!********************************************************/

typedef struct __Edge
{
  /*!
    syntax :  (node) [edge]
(left_1) .                   .(right_1)
          \ (left)  (right) /
           \._____________./
           /    [b_fcus]   \
          /                 \
(left_2) .                   .(right_2)

  */

  struct __Node *left, *rght; /*! t_node on the left/right side of the t_edge */
  short int      l_r, r_l, l_v1, l_v2, r_v1, r_v2;
  /*! these are directions (i.e., 0, 1 or 2): */
  /*! l_r (left to right) -> left[b_fcus->l_r] = right */
  /*! r_l (right to left) -> right[b_fcus->r_l] = left */
  /*! l_v1 (left t_node to first t_node != from right) -> left[b_fcus->l_v1] =
   * left_1 */
  /*! l_v2 (left t_node to secnd t_node != from right) -> left[b_fcus->l_v2] =
   * left_2 */
  /*! r_v1 (right t_node to first t_node != from left) -> right[b_fcus->r_v1] =
   * right_1 */
  /*! r_v2 (right t_node to secnd t_node != from left) -> right[b_fcus->r_v2] =
   * right_2 */

  struct __NNI   *nni;
  struct __Edge  *next;
  struct __Edge  *prev;
  struct __Edge  *next_mixt;
  struct __Edge  *prev_mixt;
  struct __Label *label;

  int         num;       /*! branch number */
  scalar_dbl *l;         /*! branch length */
  scalar_dbl *best_l;    /*! best branch length found so far */
  scalar_dbl *l_old;     /*! old branch length */
  scalar_dbl *l_var;     /*! variance of edge length */
  scalar_dbl *l_var_old; /*! variance of edge length (previous value) */

  int bip_score; /*! score of the bipartition generated by the corresponding
edge bip_score = 1 iif the branch is found in both trees to be compared,
bip_score = 0 otherwise. */
  phydbl tdist_score; /*! average transfer distance over bootstrap trees */

  int num_st_left; /*! number of the subtree on the left side */
  int num_st_rght; /*! number of the subtree on the right side */

  phydbl *Pij_rr;  /*! matrix of change probabilities and its first and secnd
                      derivates (rate*state*state) */
  phydbl *tPij_rr; /*! transpose matrix of change probabilities and its first
                      and secnd derivates (rate*state*state) */
#ifdef BEAGLE
  int Pij_rr_idx;
#endif

  short int *div_post_pred_left; /*! posterior prediction of nucleotide/aa
                                    diversity (left-hand subtree) */
  short int *div_post_pred_rght; /*! posterior prediction of nucleotide/aa
                                    diversity (rght-hand subtree) */
  short int does_exist;

  phydbl *p_lk_left,
      *p_lk_rght; /*! likelihoods of the subtree on the left and right side (for
                     each site and each relative rate category) */
  phydbl *p_lk_tip_r, *p_lk_tip_l;

#ifdef BEAGLE
  int p_lk_left_idx, p_lk_rght_idx;
  int p_lk_tip_idx;
#endif

  int *patt_id_left;
  int *patt_id_rght;
  int *p_lk_loc_left;
  int *p_lk_loc_rght;

  int *pars_l, *pars_r; /*! parsimony of the subtree on the left and right sides
                           (for each site) */
  int *ui_l, *ui_r; /*! union - intersection vectors used in Fitch's parsimony
                       algorithm */
  int *p_pars_l, *p_pars_r; /*! conditional parsimony vectors */

  /*! Below are the likelihood scaling factors (used in functions
    `Get_All_Partial_Lk_Scale' in lk.c. */
  /*
    For every site, every subtree and every rate class, PhyML maintains
    a`sum_scale_pow' value where sum_scale_pow = sum_scale_pow_v1 +
    sum_scale_pow_v2 + curr_scale_pow' sum_scale_pow_v1 and sum_scale_pow_v2 are
    sum_scale_pow of the left and right subtrees. curr_scale_pow is an integer
    greater than one. The smaller the partials, the larger curr_scale_pow.

    Now the partials for this subtree are scaled by *multiplying* each of
    them by 2^curr_scale_pow. The reason for doing the scaling this way is
    that multiplications by 2^x (x an integer) can be done in an 'exact'
    manner (i.e., there is no loss of numerical precision)

    At the root edge, the log-likelihood is then
    logL = logL' - (sum_scale_pow_left + sum_scale_pow_right)log(2),
    where L' is the scaled likelihood.
  */

  int *sum_scale_left_cat;
  int *sum_scale_rght_cat;
  int *sum_scale_left;
  int *sum_scale_rght;

  phydbl bootval; /*! bootstrap value (if exists) */

  short int is_alive; /*! is_alive = 1 if this t_edge is used in a tree */

  phydbl dist_btw_edges;
  int    topo_dist_btw_edges;

  int has_zero_br_len;

  phydbl ratio_test;     /*! approximate likelihood ratio test */
  phydbl alrt_statistic; /*! aLRT statistic */
  phydbl support_val;

  int n_jumps; /*! number of jumps of substitution rates */

  int *n_diff_states_l; /*! Number of different states found in the subtree on
                           the left of this edge */
  int *n_diff_states_r; /*! Number of different states found in the subtree on
                           the right of this edge */

  phydbl bin_cod_num;

  short int update_partial_lk_left;
  short int update_partial_lk_rght;

} t_edge;

/*!********************************************************/

typedef struct __Tree
{

  struct __Node   *n_root;  /*! root t_node */
  struct __Edge   *e_root;  /*! t_edge on which lies the root */
  struct __Node  **a_nodes; /*! array of nodes that defines the tree topology */
  struct __Edge  **a_edges; /*! array of edges */
  struct __Model  *mod;     /*! substitution model */
  struct __Calign *data;    /*! sequences */
  struct __Tree   *next; /*! set to NULL by default. Used for mixture models */
  struct __Tree   *prev; /*! set to NULL by default. Used for mixture models */
  struct __Tree
      *next_mixt; /*! set to NULL by default. Used for mixture models */
  struct __Tree
      *prev_mixt; /*! set to NULL by default. Used for mixture models */
  struct __Tree
      *mixt_tree; /*! set to NULL by default. Used for mixture models */
  struct __Tree **aux_tree;   /*! set to NULL by default. Used as a latent
                                 variable in molecular dating */
  struct __Option *io;        /*! input/output */
  struct __Matrix *mat;       /*! pairwise distance matrix */
  struct __Node  **curr_path; /*! list of nodes that form a path in the tree */
  struct __SPR   **spr_list_one_edge;
  struct __SPR   **spr_list_all_edge;
  struct __SPR    *best_spr;
  struct __Tdraw
      *ps_tree; /*! structure for drawing trees in postscript format */
  struct __T_Rate       *rates; /*! structure for handling rates of evolution */
  struct __T_Time       *times; /*! structure for handling node ages */
  struct __Tmcmc        *mcmc;
  struct __Phylogeo     *geo;
  struct __Migrep_Model *mmod;
  struct __Disk_Event   *young_disk;  /*! Youngest disk (i.e., disk which age is
                                         the closest to present). Used in PhyREX */
  struct __Disk_Event *old_samp_disk; /*! Oldest sampled disk. Used in PhyREX */
  struct __XML_node   *xml_root;
  struct __Generic_LL *edge_list;
  struct __Generic_LL *node_list;
  struct __Independent_Contrasts
      *ctrst; /*! Pointer to data structure used for independent contrasts */
  struct __Continuous_Model *contmod;

#if (defined(__AVX__) || defined(__AVX2__))
  __m256d *_tPij1, *_tPij2, *_pmat1plk1, *_pmat2plk2, *_plk0, *_l_ev, *_r_ev,
      *_prod_left, *_prod_rght;
#elif (defined(__SSE__) || defined(__SSE2__) || defined(__SSE3__) ||           \
       defined(__ARM_NEON))
  __m128d *_tPij1, *_tPij2, *_pmat1plk1, *_pmat2plk2, *_plk0, *_l_ev, *_r_ev,
      *_prod_left, *_prod_rght;
#endif

  phydbl *p_lk_left_pi, *l_ev;
  phydbl *big_lk_array;
  int     big_lk_array_pos;

  short int eval_alnL; /*! Evaluate likelihood for genetic data */
  short int eval_rlnL; /*! Evaluate likelihood for rates along the tree */
  short int eval_glnL; /*! Evaluate likelihood under phylogeo model */
  short int eval_tlnL; /*! Evaluate likelihood under tree-generating process */
  short int scaling_method;
  short int fully_nni_opt;
  short int numerical_warning;

  short int use_eigen_lr;
  int       is_mixt_tree;
  int       tree_num;  /*! tree number. Used for mixture models */
  int ps_page_number;  /*! when multiple trees are printed, this variable give
                          the current page number */
  int depth_curr_path; /*! depth of the t_node path defined by curr_path */
  int has_bip;         /*!if has_bip=1, then the structure to compare
 tree topologies is allocated, has_bip=0 otherwise */
  int n_otu;           /*! number of taxa */
  int curr_site;       /*! current site of the alignment to be processed */
  int curr_catg; /*! current class of the discrete gamma rate distribution */
  int n_swap;    /*! number of NNIs performed */
  int has_branch_lengths; /*! =1 iff input tree displays branch lengths */
  short int both_sides;   /*! both_sides=1 -> a pre-order and a post-order tree
     traversals are required to compute the likelihood
     of every subtree in the phylogeny*/
  int num_curr_branch_available; /*!gives the number of the next cell in a_edges
                                    that is free to receive a pointer to a
                                    branch */
  short int *t_dir;
  int        n_moves;
  int        verbose;

  int dp;        /*! Data partition */
  int s_mod_num; /*! Substitution model number */
  int lock_topo; /*! = 1 any subsequent topological modification will be
                    banished */
  int       print_labels;
  int       write_br_lens;
  int      *mutmap; /*! Mutational map */
  int       json_num;
  short int update_eigen_lr;
  int       tip_root; /*! Index of tip node used as the root */
  phydbl   *dot_prod;

  phydbl *expl;

  phydbl  init_lnL;
  phydbl  best_lnL;  /*! highest value of the loglikelihood found so far */
  int     best_pars; /*! highest value of the parsimony found so far */
  phydbl  c_lnL;     /*! loglikelihood */
  phydbl  p_lnL;     /*! loglikelihood (previous value) */
  phydbl  old_lnL;   /*! old loglikelihood */
  phydbl  sum_min_sum_scale; /*! common factor of scaling factors */
  phydbl *c_lnL_sorted;      /*! used to compute c_lnL by adding sorted terms to
                                minimize CPU errors */
  phydbl *cur_site_lk;       /*! vector of loglikelihoods at individual sites */
  phydbl *old_site_lk;       /*! vector of likelihoods at individual sites */
  phydbl  annealing_temp;    /*! annealing temperature in simulated annealing
                                optimization algo */
  phydbl c_dlnL;  /*! First derivative of the log-likelihood with respect to the
                     length of a branch */
  phydbl c_d2lnL; /*! Second derivative of the log-likelihood with respect to
                     the length of a branch */

  phydbl *unscaled_site_lk_cat; /*! partially scaled site likelihood at
                                   individual sites */

  phydbl *site_lk_cat; /*! loglikelihood at a single site and for each class of
                          rate*/
  phydbl unconstraint_lk; /*! unconstrained (or multinomial) log-likelihood  */
  phydbl composite_lk;    /*! composite log-likelihood  */
  int   *fact_sum_scale;
  phydbl **log_lks_aLRT; /*! used to compute several branch supports */
  phydbl   n_root_pos;   /*! position of the root on its t_edge */
  phydbl   size;         /*! tree size */
  int     *site_pars;
  int      c_pars;
  int     *step_mat;

  int size_spr_list_one_edge;
  int size_spr_list_all_edge;
  int perform_spr_right_away;

  time_t t_beg;
  time_t t_current;

  int bl_from_node_stamps; /*! == 1 -> Branch lengths are determined by t_node
                              times */
  phydbl sum_y_dist_sq;
  phydbl sum_y_dist;
  phydbl tip_order_score;
  int    write_tax_names;
  int    update_alias_subpatt;

  phydbl geo_mig_sd; /*! standard deviation of the migration step random
                        variable */
  phydbl geo_lnL;    /*! log likelihood of the phylo-geography model */

  int bl_ndigits;

  phydbl *short_l;   /*! Vector of short branch length values */
  int     n_short_l; /*! Length of short_l */
  phydbl  norm_scale;

  short int br_len_recorded;

  short int apply_lk_scaling; /*! Applying scaling of likelihoods. YES/NO */

  phydbl *K; /*! a vector of the norm.constants for the node times prior. */

  short int ignore_root;
  short int ignore_mixt_info;
#ifdef BEAGLE
  int b_inst; /*! The BEAGLE instance id associated with this tree. */
#endif

  // Extra partial lk structure for bookkeeping
  short int *div_post_pred_extra_0;
  int       *sum_scale_cat_extra_0;
  int       *sum_scale_extra_0;
  phydbl    *p_lk_extra_0;
  phydbl    *p_lk_tip_extra_0;
  int       *patt_id_extra_0;

  short int *div_post_pred_extra_1;
  int       *sum_scale_cat_extra_1;
  int       *sum_scale_extra_1;
  phydbl    *p_lk_extra_1;
  phydbl    *p_lk_tip_extra_1;
  int       *patt_id_extra_1;

  int n_edges_traversed;
  int n_tot_bl_opt;

  bool opt_topo;

} t_tree;

/*!********************************************************/

typedef struct __Super_Tree
{
  struct __Tree *tree;
  struct __List_Tree *
      treelist; /*! list of trees. One tree for each data set to be processed */
  struct __Calign  *curr_cdata;
  struct __Option **optionlist; /*! list of pointers to input structures (used
                                   in supertrees) */

  struct __Node ***match_st_node_in_gt;
  /*!  match_st_in_gt_node[subdataset number][supertree t_node number]
   *  gives the t_node in tree estimated from 'subdataset number' that
   * corresponds to 'supertree t_node number' in the supertree
   */

  struct __Node *****map_st_node_in_gt;
  /*!  mat_st_gt_node[gt_num][st_node_num][direction] gives the
   *  t_node in gt gt_num that maps t_node st_node_num in st.
   */

  struct __Edge ***map_st_edge_in_gt;
  /*!  map_st_gt_br[gt_num][st_branch_num] gives the
   *  branch in gt gt_num that maps branch st_branch_num
   *  in st.
   */

  struct __Edge ****map_gt_edge_in_st;
  /*!  mat_gt_st_br[gt_num][gt_branch_num][] is the list of
   *  branches in st that map branch gt_branch_num
   *  in gt gt_num.
   */

  int **size_map_gt_edge_in_st;
  /*!  size_map_gt_st_br[gt_num][gt_branch_num] gives the
   *  size of the list map_gt_st_br[gt_num][gt_branch_num][]
   */

  struct __Edge ***match_st_edge_in_gt;
  /*! match_st_edge_in_gt[gt_num][st_branch_num] gives the
   * branch in gt gt_num that matches branch st_branch_num
   */

  struct __Edge ***match_gt_edge_in_st;
  /*! match_gt_edge_in_st[gt_num][gt_branch_num] gives the
   * branch in st that matches branch gt_branch_num
   */

  struct __Node ****closest_match;
  /*! closest_match[gt_num][st_node_num][dir] gives the
   * closest t_node in st that matches a t_node in gt gt_num
   */

  int ***closest_dist;
  /*! closest_dist[gt_num][st_node_num][dir] gives the
   * number of edges to traverse to get to node
   * closest_match[gt_num][st_node_num][dir]
   */

  int n_part;
  /*! number of trees */

  phydbl **bl;
  /*! bl[gt_num][gt_branch_num] gives the length of
   * branch gt_branch_num
   */

  phydbl **bl_cpy;
  /*! copy of bl */

  phydbl **bl0;
  /*! bl estimated during NNI (original topo)
   * See Mg_NNI.
   */

  phydbl **bl1;
  /*! bl estimated during NNI (topo conf 1)
   * See Mg_NNI.
   */

  phydbl **bl2;
  /*! bl estimated during NNI (topo conf 2)
   * See Mg_NNI.
   */

  int *bl_partition;
  /*! partition[gt_num] gives the t_edge partition number
   * gt_num belongs to.
   */
  int n_bl_part;

  struct __Model **s_mod; /*! substitution model */

  int n_s_mod;
  int lock_br_len;

} supert_tree;

/*!********************************************************/

typedef struct __List_Tree
{ /*! a list of trees */
  struct __Tree **tree;
  int             list_size; /*! number of trees in the list */
} t_treelist;

/*!********************************************************/

typedef struct __Align
{
  char      *name;      /*! sequence name */
  int        len;       /*! sequence length */
  char      *state;     /*! sequence itself */
  short int *d_state;   /*! sequence itself (digits) */
  short int *is_ambigu; /*! is_ambigu[site] = 1 if state[site] is an ambiguous
                           character. 0 otherwise */
  short int is_duplicate;
  int       num;
} align;

/*!********************************************************/

typedef struct __Calign
{
  struct __Align **c_seq;         /*! compressed sequences      */
  struct __Align **c_seq_rm;      /*! removed sequences      */
  struct __Option *io;            /*! input/output */
  phydbl          *obs_state_frq; /*! observed state frequencies */
  short int       *invar;         /*! < 0 -> polymorphism observed */
  phydbl          *wght;          /*! # of each site in c_align */
  short int *ambigu; /*! ambigu[i]=1 is one or more of the sequences at site
                    i display an ambiguous character */
  phydbl obs_pinvar;
  int    n_otu;        /*! number of taxa */
  int    n_rm;         /*! number of taxa removed */
  int    clean_len;    /*! uncrunched sequences lenghts without gaps */
  int    n_pattern;    /*! crunched sequences lengths */
  int    init_len;     /*! length of the uncompressed sequences */
  int   *sitepatt;     /*! this array maps the position of the patterns in the
                      compressed alignment to the positions in the uncompressed
                      one */
  int         format;  /*! 0 (default): PHYLIP. 1: NEXUS. */
  scalar_dbl *io_wght; /*! weight of each *site* (not pattern) given as input */

  int n_masked;         /* Number of masked positions (or columns, depending on
                           mask_type) */
  short int mask_type;  /* MASK_TYPE_POSITION or MASK_TYPE_COLUMN */
  int      *masked_pos; /* Vector of masked positions/columns */

} calign;

/*!********************************************************/

typedef struct __Matrix
{                          /*! mostly used in BIONJ */
  phydbl **P, **Q, **dist; /*! observed proportions of transition, transverion
                  and  distances between pairs of  sequences */

  t_tree *tree; /*! tree... */
  int *on_off; /*! on_off[i]=1 if column/line i corresponds to a t_node that has
      not been agglomerated yet */
  int             n_otu; /*! number of taxa */
  char          **name;  /*! sequence names */
  int             r; /*! number of nodes that have not been agglomerated yet */
  struct __Node **tip_node; /*! array of pointer to the leaves of the tree */
  int             curr_int; /*! used in the NJ/BIONJ algorithms */
  int             method; /*! if method=1->NJ method is used, BIONJ otherwise */
} matrix;

/*!********************************************************/

typedef struct __RateMatrix
{
  int n_diff_rr; /*! number of different relative substitution rates in the
                    custom model */
  vect_dbl
      *rr; /*! relative rate parameters of the GTR or custom model (rescaled) */
  vect_dbl *rr_val; /*! log of relative rate parameters of the GTR or custom
                       model (unscaled) */
  vect_int *rr_num;
  vect_int *n_rr_per_cat; /*! number of rate parameters in each category */
  vect_dbl *qmat;
  vect_dbl *qmat_buff;

  bool optimize;

  struct __RateMatrix *next;
  struct __RateMatrix *prev;
} t_rmat;

/*!********************************************************/

typedef struct __RAS
{
  /*! Rate across sites */
  int n_catg; /*! number of categories in the discrete gamma distribution */
  int invar;  /*! =1 iff the substitution model takes into account invariable
                 sites */
  int gamma_median; /*! 1: use the median of each bin in the discrete gamma
                       distribution. 0: the mean is used */
  vect_dbl *gamma_r_proba; /*! probabilities of the substitution rates defined
                              by the discrete gamma distribution */
  vect_dbl *gamma_r_proba_unscaled;
  vect_dbl *gamma_rr; /*! substitution rates defined by the RAS distribution */
  vect_dbl *gamma_rr_unscaled; /*! substitution rates defined by the RAS
                                  distribution (unscaled) */
  scalar_dbl *alpha;           /*! gamma shapa parameter */
  int         free_mixt_rates;
  scalar_dbl
      *free_rate_mr; /*! mean relative rate as given by the FreeRate model */

  int         parent_class_number;
  scalar_dbl *pinvar; /*! proportion of invariable sites */

  short int init_rr;
  short int init_r_proba;
  short int normalise_rr;

  short int
      *skip_rate_cat; /*! indicates whether of not the the likelihood for a
                         given rate class shall be calculated. The default model
                         in PhyML has four rate classes and the likelihood at a
                         given site is then \sum_{i=1}^{4} \Pr(D|R=i) \Pr(R=i).
                         Now, when optimising the rate for say the first rate
                         class (i=1) one does not need to re-compute \Pr(D|R=2),
                         \Pr(D|R=3) and \Pr(D|R=4) (which is the time consuming
                         part of the likelihood calculation). This is where
                          'skip_rate_category' comes in handy. */

  short int sort_rate_classes; /*! When doing MCMC moves on rate classes, one
                                  needs to have the rate classes sorted */

  bool optimize;

  struct __RAS *next;
  struct __RAS *prev;

} t_ras;

/*!********************************************************/

typedef struct __EquFreq
{
  /*! Equilibrium frequencies */
  vect_dbl *pi;          /*! states frequencies */
  vect_dbl *pi_unscaled; /*! states frequencies (unscaled) */

  phydbl *b_frq; /*! vector of empirical state frequencies */

  vect_dbl *user_b_freq; /*! user-defined nucleotide frequencies */

  short int type;

  struct __EquFreq *next;
  struct __EquFreq *prev;

} t_efrq;

/*!********************************************************/

typedef struct __Model
{
  struct __Optimiz    *s_opt; /*! pointer to parameters to optimize */
  struct __Eigen      *eigen;
  struct __M4         *m4mod;
  struct __Option     *io;
  struct __Model      *next;
  struct __Model      *prev;
  struct __Model      *next_mixt;
  struct __Model      *prev_mixt;
  struct __RateMatrix *r_mat;
  struct __EquFreq    *e_frq;
  struct __RAS        *ras;

  t_string *aa_rate_mat_file;
  FILE     *fp_aa_rate_mat;

  t_string *modelname;
  t_string *custom_mod_string; /*! string of characters used to define custom
                                  models of substitution */

  int mod_num; /*! model number */

  int update_eigen; /*! update_eigen=1-> eigen values/vectors need to be updated
                     */
  int cv_type;      /* Type of cross-validation method */

  int whichmodel;
  int is_mixt_mod;
  int augmented;
  int ns; /*! number of states (4 for ADN, 20 for AA) */

  int use_m4mod; /*! Use a Markov modulated Markov model ? */

  scalar_dbl *kappa;  /*! transition/transversion rate */
  scalar_dbl *lambda; /*! parameter used to define the ts/tv ratios in the F84
                         and TN93 models */
  scalar_dbl *br_len_mult; /*! when users want to fix the relative length of
                              edges and simply estimate the total length of the
                              tree. This multiplier does the trick */
  scalar_dbl *br_len_mult_unscaled;

  vect_dbl   *Pij_rr; /*! matrix of change probabilities */
  scalar_dbl *mr;     /*! mean rate = branch length/time interval  mr =
                         -sum(i)(vct_pi[i].mat_Q[ii]) */
  scalar_dbl *aic;
  scalar_dbl *bic;

  short int
      log_l; /*! Edge lengths are actually log(Edge lengths) if log_l == YES !*/
  phydbl l_min; /*! Minimum branch length !*/
  phydbl l_max; /*! Maximum branch length !*/

  scalar_dbl *l_var_sigma; /*! For any edge b we have b->l_var->v = l_var_sigma
                            * (b->l->v)^2 */
  phydbl l_var_min; /*! Min of variance of branch lengths (used in conjunction
                       with gamma_mgf_bl == YES) */
  phydbl l_var_max; /*! Max of variance of branch lengths (used in conjunction
                       with gamma_mgf_bl == YES) */

  int gamma_mgf_bl; /*! P = \int_0^inf exp(QL) p(L) where L=\int_0^t R(s) ds and
                       p(L) is the gamma density. Set to NO by default !*/

  int n_mixt_classes; /* Number of classes in the mixture model. */

  scalar_dbl *r_mat_weight;
  scalar_dbl *e_frq_weight;
#ifdef BEAGLE
  int  b_inst;
  bool optimizing_topology; /*! This is a flag that prevents the resetting of
                               category weights. Why? Read */
  /*  Recall that while optimizing the topology, PhyML temporarily only uses 2
   *  rate categories. Recall also that a BEAGLE instance is created with all
   * the required categories, but we temporarily assign 0 weight to the other
   * categories thus effectively using only 2 categories. However, subsequent
   * calls to update the rates (i.e. update_beagle_ras()) will reset the
   * weights. This flag prevents this resetting from happening */
#endif

} t_mod;

/*!********************************************************/

typedef struct __Eigen
{
  int     size; /*! matrix is size * size */
  phydbl *q;    /*! matrix for which eigen values and vectors are computed */
  phydbl *space;
  int    *space_int;
  phydbl *e_val;       /*! eigen values (vector), real part. */
  phydbl *e_val_im;    /*! eigen values (vector), imaginary part */
  phydbl *r_e_vect_im; /*! right eigen vector (matrix), imaginary part */
  phydbl *l_e_vect;    /*! left eigen vector (matrix), real part */
  phydbl *r_e_vect;    /*! right eigen vector (matrix), real part */
  phydbl *dum;         /*! dummy ns*ns space */

  struct __Eigen *prev;
  struct __Eigen *next;
} eigen;

/*!********************************************************/

typedef struct __Option
{                             /*! mostly used in 'help.c' */
  struct __Model  *mod;       /*! pointer to a substitution model */
  struct __Tree   *tree;      /*! pointer to the current tree */
  struct __Align **data;      /*! pointer to the uncompressed sequences */
  struct __Tree   *cstr_tree; /*! pointer to a constraint tree (can be a
                                 multifurcating one) */
  struct __Calign     *cdata; /*! pointer to the compressed sequences */
  struct __Super_Tree *st;    /*! pointer to supertree */
  struct __Tnexcom   **nex_com_list;
  struct __List_Tree  *treelist; /*! list of trees. */
  struct __Option     *next;
  struct __Option     *prev;
  struct __Tmcmc      *mcmc;
  struct __T_Rate     *rates;
  struct __T_Time     *times;

  int interleaved; /*! interleaved or sequential sequence file format ? */
  int in_tree;     /*! =1 iff a user input tree is used as input */

  char *in_align_file; /*! alignment file name */
  FILE *fp_in_align;   /*! pointer to the alignment file */

  char *in_tree_file; /*! input tree file name */
  FILE *fp_in_tree;   /*! pointer to the input tree file */

  char *in_constraint_tree_file; /*! input constraint tree file name */
  FILE *fp_in_constraint_tree; /*! pointer to the input constraint tree file */

  char *out_tree_file; /*! name of the tree file */
  FILE *fp_out_tree;

  char *weight_file; /*! name of the file containing site weights */
  FILE *fp_weight_file;

  char *out_trees_file; /*! name of the tree file */
  FILE *fp_out_trees;   /*! pointer to the tree file containing all the trees
                           estimated using random starting trees */

  char *out_boot_tree_file; /*! name of the tree file */
  FILE *fp_out_boot_tree;   /*! pointer to the bootstrap tree file */

  char *out_boot_stats_file; /*! name of the tree file */
  FILE *fp_out_boot_stats;   /*! pointer to the statistics file */

  char *out_stats_file; /*! name of the statistics file */
  FILE *fp_out_stats;

  char *out_trace_file; /*! name of the file in which the trace is written */
  FILE *fp_out_trace;

  char *out_json_trace_file; /*! name of the file in which json trace is written
                              */
  FILE *fp_out_json_trace;

  char *out_lk_file; /*! name of the file in which the likelihood of the model
                        is written */
  FILE *fp_out_lk;

  char *out_summary_file; /*! name of the file in which summary statistics are
                             written */
  FILE *fp_out_summary;

  char *out_ps_file; /*! name of the file in which tree(s) is(are) written */
  FILE *fp_out_ps;

  char *out_ancestral_seq_file;  /*! name of the file containing the ancestral
                                    sequences */
  char *out_ancestral_tree_file; /*! name of the file containing the tree with
                                    internal node labelled according to refs in
                                    ancestral_seq_file */

  FILE *fp_out_ancestral_seq;  /*! pointer to the file containing the ancestral
                                  sequences */
  FILE *fp_out_ancestral_tree; /*! pointer to the file containing the tree with
                                  labels on internal nodes  */

  char *in_xml_file;
  FILE *fp_in_xml; /*! pointer to the file containing XML formatted data */

  char *in_coord_file; /*! name of input file containing coordinates */
  FILE *fp_in_coord;   /*! pointer to the file containing coordinates */

  char *out_file; /*! name of the output file */

  char *clade_list_file;

  int datatype;         /*! 0->DNA, 1->AA */
  int print_boot_trees; /*! =1 if the bootstrapped trees are printed in output
                         */
  int   out_stats_file_open_mode; /*! opening file mode for statistics file */
  int   out_tree_file_open_mode;  /*! opening file mode for tree file */
  int   n_data_sets;              /*! number of data sets to be analysed */
  int   n_trees;                  /*! number of trees */
  int   init_len;                 /*! sequence length */
  int   n_otu;                    /*! number of taxa */
  int   n_data_set_asked;         /*! number of bootstrap replicates */
  char *nt_or_cd;                 /*! nucleotide or codon data ? (not used) */
  int   multigene;                /*! if=1 -> analyse several partitions. */
  int   config_multigene;
  int   n_part; /*! number of data partitions */
  int   curr_gt;
  int   ratio_test; /*! from 1 to 4 for specific branch supports, 0 of not */
  int   ready_to_go;
  int   data_file_format; /*! Data format: Phylip or Nexus */
  int   tree_file_format; /*! Tree format: Phylip or Nexus */
  int   state_len;

  int curr_interface;
  int r_seed;        /*! random seed */
  int collapse_boot; /*! 0 -> branch length on bootstrap trees are not collapsed
                        if too small */
  int random_boot_seq_order; /*! !0 -> sequence order in bootstrapped data set
                                is random */
  int print_trace;
  int print_json_trace;
  int print_site_lnl;
  int m4_model;
  int rm_ambigu; /*! 0 is the default. 1: columns with ambiguous characters are
                    discarded prior further analysis */
  int   colalias;
  int   append_run_ID;
  char *run_id_string;
  int   quiet; /*! 0 is the default. 1: no interactive question (for batch mode)
                */
  int    lk_approx; /* EXACT or NORMAL */
  char **alphabet;
  int    codpos;
  int    mutmap;
  int    use_xml;

  char **long_tax_names;
  char **short_tax_names;
  int    size_tax_names;

  phydbl *z_scores;
  phydbl *lat;
  phydbl *lon;

  int boot_prog_every;

  int mem_question;
  int do_alias_subpatt;

#ifdef BEAGLE
  int beagle_resource;
#endif

  int ancestral;
  int has_io_weights;
  int tbe_bootstrap; /* Replace standard bootstrap with tbe bootstrap (only when
                        b>0) */

  int leave_duplicates; /* Leave duplicated sequences */
  int precision;        /* Decimal output precision for values in stats file */
  int n_boot_replicates;

  short int print_mat_and_exit;
  short int print_node_num; /*! print node numbers if print_node_num=1 */
  short int print_support_val;

  short int do_tbe;
  short int do_boot;
  short int do_bayesboot;
  short int do_alrt;

  short int edge_len_unit;

  short int mcmc_output_times;
  short int mcmc_output_trees;

} option;

/*!********************************************************/

typedef struct __Optimiz
{
  short int opt_subst_param; /*! if opt_topo=0 and opt_subst_param=1 -> the
          numerical parameters of the model are optimised. if opt_topo=0 and
          opt_free_param=0 -> no parameter is optimised */
  short int opt_clock_r;
  short int opt_bl_one_by_one; /*! =1 -> the branch lengths are optimised */
  short int opt_topo;          /*! =1 -> the tree topology is optimised */
  short int topo_search;
  short int opt_node_ages;
  short int opt_neff;

  phydbl init_lk; /*! initial loglikelihood value */
  int n_it_max; /*! maximum bnumber of iteration during an optimisation step */
  int last_opt; /*! =1 -> the numerical parameters are optimised further while
                   the tree topology remains fixed */
  int random_input_tree; /*! boolean */
  int n_rand_starts;     /*! number of random starting points */
  int brent_it_max;
  int steph_spr;
  int opt_five_branch;
  int pars_thresh;
  int hybrid_thresh;
  int opt_br_len_mult;
  int min_n_triple_moves;
  int max_rank_triple_move;
  int n_improvements;
  int max_spr_depth;
  int max_no_better_tree_found;

  phydbl tree_size_mult; /*! tree size multiplier */
  phydbl min_diff_lk_local;
  phydbl min_diff_lk_global;
  phydbl min_diff_lk_move;
  phydbl p_moves_to_examine;
  int    fast_nni;
  int    greedy;
  int    general_pars;
  int    quickdirty;
  int    spr_pars;
  int    spr_lnL;
  int    max_depth_path;
  int    min_depth_path;
  int    deepest_path;
  int    eval_list_regraft;
  phydbl max_delta_lnL_spr;
  phydbl max_delta_lnL_spr_current;
  phydbl worst_lnL_spr;
  int    br_len_in_spr;
  int    opt_free_mixt_rates;
  int    constrained_br_len;
  int    opt_gamma_br_len;
  int    first_opt_free_mixt_rates;
  int    wim_n_rgrft;
  int    wim_n_globl;
  int    wim_max_dist;
  int    wim_n_optim;
  int    wim_n_best;
  int    wim_inside_opt;

  int opt_rmat_weight;
  int opt_efrq_weight;

  int skip_tree_traversal;
  int serial_free_rates;

  int curr_opt_free_rates;

  int nni_br_len_opt;

  int apply_spr_right_away;
  int apply_spr;

  phydbl l_min_spr;
} t_opt;

/*!********************************************************/

typedef struct __NNI
{

  struct __Node *left;
  struct __Node *rght;
  struct __Edge *b;

  phydbl      score;
  scalar_dbl *init_l;
  scalar_dbl *init_v;
  phydbl      init_lk;
  scalar_dbl *best_l;
  scalar_dbl *best_v;
  phydbl      lk0, lk1, lk2;
  scalar_dbl *l0, *l1, *l2;
  scalar_dbl *v0, *v1, *v2;

  struct __Node *swap_node_v1;
  struct __Node *swap_node_v2;
  struct __Node *swap_node_v3;
  struct __Node *swap_node_v4;

  int best_conf; /*! best topological configuration :
        ((left_1,left_2),right_1,right_2) or
        ((left_1,right_2),right_1,left_2) or
        ((left_1,right_1),right_1,left_2)  */
} t_nni;

/*!********************************************************/

typedef struct __SPR
{
  struct __Node  *n_link;
  struct __Node  *n_opp_to_link;
  struct __Edge  *b_opp_to_link;
  struct __Edge  *b_target;
  struct __Edge  *b_init_target;
  struct __Node **path;

  scalar_dbl *init_target_l;
  scalar_dbl *init_target_v;

  scalar_dbl *l0, *l1, *l2;
  scalar_dbl *v0, *v1, *v2;

  phydbl lnL;
  int    depth_path;
  int    pars;
  int    dist;

  struct __SPR *next;
  struct __SPR *prev;
  struct __SPR *next_mixt;
  struct __SPR *prev_mixt;
  struct __SPR *path_prev;
  struct __SPR *path_next;

} t_spr;

/*!********************************************************/

typedef struct __Pnode
{
  struct __Pnode **next;
  int              weight;
  int              num;
} pnode;

/*!********************************************************/

typedef struct __M4
{
  int n_h; /*! number of hidden states */
  int n_o; /*! number of observable states  */
  int use_cov_alpha;
  int use_cov_free;

  phydbl **o_mats; /*! set of matrices of substitution rates across observable
                      states */
  vect_dbl *o_rr;  /*! relative rates (symmetric) of substitution between
                    observable states */
  phydbl *h_rr;    /*! relative rates (symmetric) of substitution between hidden
                      states */
  phydbl *h_mat;   /*! matrix that describes the substitutions between hidden
                      states (aka switches) */
  phydbl *o_fq;    /*! equilibrium frequencies for the observable states */

  vect_dbl *h_fq;          /*! equilibrium frequencies for the hidden states */
  vect_dbl *h_fq_unscaled; /*! unscaled equilibrium frequencies for the hidden
                            states */

  vect_dbl *multipl_unscaled; /*! unscaled  vector of values that multiply each
                                o_mats matrix */
  vect_dbl *multipl; /*! vector of values that multiply each o_mats matrix */

  scalar_dbl *delta; /*! switching rate */
  scalar_dbl *alpha; /*! gamma shape parameter */
} m4;

/*!********************************************************/

typedef struct __Tdraw
{
  phydbl *xcoord;   /*! t_node coordinates on the x axis */
  phydbl *ycoord;   /*! t_node coordinates on the y axis */
  phydbl *xcoord_s; /*! t_node coordinates on the x axis (scaled) */
  phydbl *ycoord_s; /*! t_node coordinates on the y axis (scaled) */
  int     page_width;
  int     page_height;
  int     tree_box_width;

  int    *cdf_mat;
  phydbl *cdf_mat_x;
  phydbl *cdf_mat_y;

  phydbl max_dist_to_root;
} tdraw;

/*!********************************************************/

typedef struct __T_Rate
{
  phydbl lexp; /*! Parameter of the exponential distribution that governs the
                  rate at which substitution between rate classes ocur */
  phydbl alpha;
  phydbl less_likely;
  phydbl c_lnL1;
  phydbl c_lnL2;

  phydbl c_lnL; /*! Current value of ProbDens(Br len | time stamps, model of
                   rate evolution) */
  phydbl p_lnL; /*! Previous value of c_lnL */

  phydbl c_lnP; /*! Current value of ProbDens(time stamps, model of rate
                   evolution) */
  phydbl p_lnP; /*! Previous value of c_lnP */

  phydbl    clock_r; /*! Mean substitution rate, i.e., 'molecular clock' rate */
  phydbl    clock_r_prior_mean;
  phydbl    clock_r_prior_var;
  short int clock_r_has_prior;
  short int
      init_clock_r; /* Set to YES if user gives initial value of clock rate */
  phydbl min_clock;
  phydbl max_clock;

  phydbl lbda_nu;
  phydbl min_dt;
  phydbl step_rate;
  phydbl true_tree_size;
  phydbl p_max;
  phydbl norm_fact;
  phydbl nu; /*! Parameter of the Exponential distribution for the corresponding
                model */
  phydbl min_nu;
  phydbl max_nu;
  phydbl autocor_rate_prior;
  phydbl covdet;
  phydbl sum_invalid_areas;

  phydbl min_rate;
  phydbl max_rate;

  phydbl *nd_r; /*! Current rates at nodes */
  phydbl *br_r; /*! Current rates along edges */
  phydbl *triplet;
  phydbl *true_r; /*! true t_edge rates (on rooted tree) */
  phydbl *buff_t;
  phydbl *buff_br_r;
  phydbl *buff_nd_r;
  phydbl *
      dens; /*! Probability densities of mean substitution rates at the nodes */
  phydbl *ml_l;    /*! ML t_edge lengths (rooted) */
  phydbl *cur_l;   /*! Current t_edge lengths (rooted) */
  phydbl *u_ml_l;  /*! ML t_edge lengths (unrooted) */
  phydbl *u_cur_l; /*! Current t_edge lengths (unrooted) */
  phydbl *invcov;
  phydbl *cov_r;
  phydbl *mean_r; /*! average values of br_r taken across the sampled values
                     during the MCMC */
  phydbl    *_2n_vect1;
  phydbl    *_2n_vect2;
  phydbl    *_2n_vect3;
  phydbl    *_2n_vect4;
  short int *_2n_vect5;
  phydbl    *_2n2n_vect1;
  phydbl    *_2n2n_vect2;
  phydbl    *trip_cond_cov;
  phydbl    *trip_reg_coeff;
  phydbl    *cond_var;
  phydbl    *reg_coeff;
  phydbl    *mean_l;
  phydbl    *cov_l;
  phydbl    *grad_l; /* gradient */
  phydbl     inflate_var;

  int adjust_rates; /*! if = 1, branch rates are adjusted such that a
               modification of a given t_node time does not modify any branch
               lengths */
  int use_rates; /*! if = 0, branch lengths are expressed as differences between
                    t_node times */
  int bl_from_rt; /*! if =1, branch lengths are obtained as the product of cur_r
                     and t */
  int   approx;
  int   model_id; /*! Model number */
  char *model_name;
  int   is_allocated;
  int   met_within_gibbs;

  int update_mean_l;
  int update_cov_l;

  int *n_tips_below;

  struct __Node *
      *lca; /*! 2-way table of common ancestral nodes for each pair of nodes */

  short int *br_do_updt;
  phydbl    *cur_gamma_prior_mean;
  phydbl    *cur_gamma_prior_var;

  short int br_r_recorded;

  phydbl log_K_cur;
  int    cur_comb_numb;

  struct __T_Rate *next;
  struct __T_Rate *prev;

} t_rate;

/*!********************************************************/

typedef struct __T_Time
{
  phydbl *nd_t;   /*! Current t_node times */
  phydbl *buff_t; /*! Current t_node times */
  phydbl *true_t; /*! Current t_node times */

  phydbl c_lnL; /*! ProbDens(time stamps|param) */
  phydbl c_lnP; /*! ProbDens(param) */
  phydbl p_lnL; /*! Previous value of Prob(time stamps|param) */
  phydbl p_lnP; /*! Previous value of Prob(param) */

  phydbl c_lnL_jps; /*! Prob(# Jumps | time stamps, rates, model of rate
                       evolution) */

  phydbl scaled_pop_size; // Product of effective population size with length of
                          // a generation in calendar time unit
  phydbl scaled_pop_size_min;
  phydbl scaled_pop_size_max;

  phydbl neff_growth;
  phydbl neff_growth_min;
  phydbl neff_growth_max;

  phydbl birth_rate;
  phydbl birth_rate_min;
  phydbl birth_rate_max;
  phydbl birth_rate_pivot;

  phydbl death_rate;
  phydbl death_rate_min;
  phydbl death_rate_max;
  phydbl death_rate_pivot;

  phydbl *t_prior;
  phydbl *t_prior_min;
  phydbl *t_prior_max;
  phydbl *t_floor;
  phydbl *t_mean;
  int *t_rank; /* rank of nodes, e.g., tree->nd_a[tree->rates->t_rank[0]] is the
                  oldest node */

  short int nd_t_recorded;

  short int is_asynchronous;

  phydbl *time_slice_lims;
  phydbl *calib_prob;

  struct __Calibration **a_cal; /* array of calibration data */
  int                    n_cal; /* number of elements in a_cal */

  phydbl *t_prior_min_ori;
  phydbl *t_prior_max_ori;
  phydbl *times_partial_proba;

  int *has_survived;

  int *n_jps;
  int *t_jps;

  int *numb_calib_chosen;

  int  n_time_slices;
  int *n_time_slice_spans;
  int *curr_slice;

  short int *t_has_prior;

  phydbl *mean_t; /*! average values of nd_t taken across the sampled values
                     during the MCMC */

  short int model_id;
  short int coalescent_model_id;

  int update_time_norm_const;

  short int augmented_coalescent;

  phydbl    neff_prior_mean;
  phydbl    neff_prior_var;
  short int neff_prior_distrib;

} t_time;

/*!********************************************************/

typedef struct __Tmcmc
{
  struct __Option *io;

  phydbl *tune_move;
  phydbl *move_weight;
  phydbl *move_prob;

  phydbl *acc_rate;
  int    *acc_move;
  int    *run_move;
  int    *prev_acc_move;
  int    *prev_run_move;
  int    *num_move;
  int    *move_type;
  char  **move_name;

  time_t time_beg;
  time_t time_end;

  int num_move_nd_r;
  int num_move_br_r;
  int num_move_times;
  int num_move_times_and_rates;
  int num_move_times_and_rates_root;
  int num_move_root_time;
  int num_move_nu;
  int num_move_clock_r;
  int num_move_tree_height;
  int num_move_time_slice;
  int num_move_subtree_height;
  int num_move_kappa;
  int num_move_rr;
  int num_move_spr;
  int num_move_spr_weighted;
  int num_move_spr_local;
  int num_move_spr_root;
  int num_move_tree_rates;
  int num_move_subtree_rates;
  int num_move_updown_nu_cr;
  int num_move_updown_t_cr;
  int num_move_updown_t_br;
  int num_move_ras;
  int num_move_rates_shrink;
  int num_move_cov_rates;
  int num_move_cov_switch;
  int num_move_birth_rate;
  int num_move_death_rate;
  int num_move_birth_death_updown;
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
  int num_move_phyrex_spr;
  int num_move_phyrex_spr_slide;
  int num_move_phyrex_narrow_exchange;
  int num_move_phyrex_wide_exchange;
  int num_move_phyrex_scale_times;
  int num_move_phyrex_ldscape_lim;
  int num_move_phyrex_sigsq;
  int num_move_time_neff;
  int num_move_time_neff_growth;
  int num_move_phyrex_sim;
  int num_move_phyrex_traj;
  int num_move_phyrex_indel_disk_serial;
  int num_move_phyrex_sim_plus;
  int num_move_phyrex_indel_hit_serial;
  int num_move_phyrex_ldsk_and_disk;
  int num_move_phyrex_ldsk_multi;
  int num_move_phyrex_disk_multi;
  int num_move_phyrex_ldsk_given_disk;
  int num_move_phyrex_disk_given_ldsk;
  int num_move_phyrex_add_remove_jump;
  int num_move_clade_change;
  int num_move_phyrex_ldsk_tip_to_root;
  int num_move_phyrex_sigsq_scale;
  int num_move_phyrex_ldsk_tips;
  int num_move_phyrex_node_times;
  int num_move_phyrex_node_veloc;
  int num_move_phyrex_correlated_node_veloc;
  int num_move_phyrex_all_veloc;
  int num_move_phyrex_shuffle_node_times;
  int num_move_phyrex_iwn_omega;
  int num_move_phyrex_iou_theta;
  int num_move_phyrex_iou_theta_sigsq;
  int num_move_phyrex_iou_mu;
  int num_move_obs_var;
  int num_move_phyrex_tip_loc;

  int  nd_t_digits;
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
  int  n_moves;
  int  move_idx;
  int  randomize;
  int  norm_freq;
  int  run;
  int  chain_len;
  int  sample_interval;
  int  chain_len_burnin;
  int  print_every;
  int  is_burnin;
  int  max_lag;

  phydbl max_tune;
  phydbl min_tune;

  phydbl *sampled_val;
  int     sample_size;
  int     sample_num;
  phydbl *ess;
  int    *ess_run;
  int    *start_ess;
  phydbl *mode;
  int always_yes; /* Always accept proposed move (as long as log-likelihood >
                     UNLIKELY) */
  int is;         /* Importance sampling? Yes or NO */

  short int out_verbose;

} t_mcmc;

/*!********************************************************/

typedef struct __Tpart
{
  int            *ns; /*! number of states for each partition (e.g., 2, 4, 3) */
  int            *cum_ns; /*! cumulative number of states (e.g., 0, 2, 6) */
  int             ns_max; /*! maximum number of states */
  int             ns_min; /*! minimum number of states */
  int             n_partitions; /*! number of partitions */
  struct __Eigen *eigen;
} part;

/*!********************************************************/

typedef struct __Tnexcom
{ /*! Nexus command */
  char               *name;
  int                 nparm;
  int                 nxt_token_t;
  int                 cur_token_t;
  FILE               *fp;
  struct __Tnexparm **parm;
} nexcom;

/*!********************************************************/

typedef struct __Tnexparm
{ /*! Nexus parameter value */
  char *name;
  char *value;
  int   nxt_token_t;
  int   cur_token_t;
  FILE *fp;
  int (*func)(char *, struct __Tnexparm *, struct __Option *);
  struct __Tnexcom *com;
} nexparm;

/*!********************************************************/

typedef struct __ParamInt
{
  int val;
} t_param_int;

/*!********************************************************/

typedef struct __ParamDbl
{
  phydbl val;
} t_param_dbl;

/*!********************************************************/

typedef struct __XML_node
{

  struct __XML_attr
      *attr;   // Pointer to the first element of a list of attributes
  int  n_attr; // Number of attributes
  struct __XML_node *next;   // Next sibling
  struct __XML_node *prev;   // Previous sibling
  struct __XML_node *parent; // Parent of this node
  struct __XML_node *child;  // Child of this node
  char              *id;
  char              *name;
  char              *value;
  struct __Generic_Data_Structure
      *ds; // Pointer to a data strucuture. Can be a scalar, a vector, anything.
} xml_node;

/*!********************************************************/

typedef struct __Generic_Data_Structure
{
  void                            *obj;
  struct __Generic_Data_Structure *next;
} t_ds;

/*!********************************************************/

typedef struct __XML_attr
{
  char              *name;
  char              *value;
  struct __XML_attr *next; // Next attribute
  struct __XML_attr *prev; // Previous attribute
} xml_attr;

/*!********************************************************/

typedef struct __Calibration
{
  struct __Calibration *next; // Next calibration
  struct __Calibration *prev; // Previous calibration
  struct __Clade      **clade_list;

  phydbl *
      alpha_proba_list; // list of alpha proba, one for each clade in clade_list

  int current_clade_idx; // index of the clade the calibration time interval
                         // currently applies to
  int clade_list_size;

  phydbl lower; // lower bound
  phydbl upper; // upper bound

  short int is_primary; // Is it a primary or secondary calibration interval?

  char *id; // calibration ID
} t_cal;

/*!********************************************************/

typedef struct __Clade
{
  char           *id;
  struct __Node **tip_list;  // list of tips defining the clade
  char          **tax_list;  // list of names of tips defining the clade
  int             n_tax;     // number of taxa in the clade
  struct __Node  *target_nd; // The node this calibration applies to
} t_clad;

/*!********************************************************/

typedef struct __Phylogeo
{
  phydbl *cov;   // Covariance of migrations (n_dim x n_dim)
  phydbl *r_mat; // R matrix. Gives the rates of migrations between locations.
                 // See article.
  phydbl *f_mat; // F matrix. See article.
  int *occup; // Vector giving the number of lineages that occupy each location
  int *idx_loc;         // Index of location for each lineage
  int *idx_loc_beneath; // Gives the index of location occupied beneath each
                        // node in the tree
  int ldscape_sz;       // Landscape size: number of locations
  int n_dim; // Dimension of the data (e.g., longitude + lattitude -> n_dim = 2)
  int update_fmat;
  struct __Geo_Coord **coord_loc; // Coordinates of the observed locations

  phydbl sigma; // Dispersal parameter
  phydbl min_sigma;
  phydbl max_sigma;
  phydbl sigma_thresh; // beyond sigma_thresh, there is no dispersal bias.

  phydbl lbda; // Competition parameter
  phydbl min_lbda;
  phydbl max_lbda;

  phydbl c_lnL; // ProbDens(location|param)
  phydbl p_lnL; // Previous value of c_lnL

  phydbl c_lnP; // ProbDens(param)
  phydbl p_lnP; // Previous value of p_lnP

  struct __Node **sorted_nd; // Table of nodes sorted wrt their heights.

  phydbl tau; // overall migration rate parameter
  phydbl min_tau;
  phydbl max_tau;

  phydbl dum; // dummy parameter use to assess non-identifiability issues
  phydbl min_dum;
  phydbl max_dum;

} t_geo;

/*!********************************************************/
// Data structure for the migration/reproduction model
typedef struct __Migrep_Model
{
  struct __Geo_Coord *lim_up; // max longitude and lattitude
  struct __Geo_Coord *lim_do; // min longitude and lattitude
  phydbl
      *sigsq_scale; // Scaling factors for the variance parameter in RRW model

  phydbl sigsq_scale_min;
  phydbl sigsq_scale_max;

  phydbl sigsq_scale_norm_fact;

  short int model_id;
  short int dist_type;
  int       n_dim;
  int       safe_phyrex;
  int       max_num_of_intervals;
  int       update_rad;

  short int integrateAncestralLocations;

  phydbl lbda;             // rate at which events occur
  phydbl min_lbda;         // min of rate at which events occur
  phydbl max_lbda;         // max of rate at which events occur
  phydbl prior_param_lbda; // parameter of the parameter for the prior on lbda

  phydbl mu;             // per-capita and per event death probability
  phydbl min_mu;         // min of per-capita and per event death probability
  phydbl max_mu;         // max of per-capita and per event death probability
  phydbl prior_param_mu; // parameter of the parameter for the prior on mu

  phydbl rad;             // radius of the migrep disk
  phydbl min_rad;         // min of radius of the migrep disk
  phydbl max_rad;         // max of radius of the migrep disk
  phydbl prior_param_rad; // parameter of the parameter for the prior on radius

  phydbl
      *sigsq; // parent to offspring distance variance (i.e., gene flow)
              // parameter. First elem is for latitude, second is for longitude
  phydbl min_sigsq;         // min
  phydbl max_sigsq;         // max
  phydbl prior_param_sigsq; // parameter of the parameter for the prior

  phydbl rho;          // intensity parameter of the Poisson point processs
  phydbl gen_cal_time; // duration of one generation in calendar time unit
  phydbl nu; // parameter of hyperprior on sigsq_scale (see Eq. (1) in Lemey et
             // al., 2010).

  phydbl omega; // step size in white noise jump model
  phydbl min_omega;
  phydbl max_omega;

  phydbl ou_theta; // "stickiness" parameter in OU model
  phydbl min_ou_theta;
  phydbl max_ou_theta;

  phydbl *ou_mu; // drift parameter in OU model
  phydbl  min_ou_mu;
  phydbl  max_ou_mu;

  phydbl min_veloc;
  phydbl max_veloc;

  phydbl c_lnL; // current value of log-likelihood
  phydbl p_lnL; // previous value of log-likelihood

  phydbl c_lnP; // current value of log-prior
  phydbl p_lnP; // previous value of log-prior

  phydbl c_ln_prior_rad;   // current value of log prior for the prior on radius
  phydbl c_ln_prior_lbda;  // current value of log prior for the prior on lbda
  phydbl c_ln_prior_mu;    // current value of log prior for the prior on mu
  phydbl c_ln_prior_sigsq; // current value of log prior for the prior on
                           // sigsq=4.pi.lbda.mu.rad^4

  phydbl soft_bound_area;

  short int sampling_scheme;
  short int use_locations;

  // Prior mean, sd and distribution for sigsq parameter, controlling either the
  // variation of coordinates or velocities (and coordinates)
  phydbl    rw_prior_mean;
  phydbl    rw_prior_sd;
  short int rw_prior_distrib;

  phydbl *rw_root_mean;
  phydbl *rw_root_var;

  phydbl rrw_prior_sd; // value of parameter governing the variance of diffusion
                       // coefficient across edges

  short int print_lk;

  short int do_location_sampling;
} t_phyrex_mod;

/*!********************************************************/

typedef struct __Disk_Event
{
  struct __Geo_Coord  *centr;
  struct __Disk_Event *next;
  struct __Disk_Event *prev;
  struct __Disk_Event *img;
  struct __Lindisk_Node *
      *ldsk_a; // array of lindisk nodes corresponding to this disk event.
  struct __Lindisk_Node *ldsk;
  struct __Migrep_Model *mmod;

  phydbl cum_glnL; // cumulative log-likelihood (locations)
  phydbl cum_tlnL; // cumulative log-likelihood (times);
  phydbl time;
  phydbl time_bkp;

  short int age_fixed; // time is fixed for disks corresponding to samples.

  int   n_ldsk_a; // size of ldsk_a.
  char *id;

} t_dsk;

/*!********************************************************/

// Data structure for implementing Felsenstein's independent contrasts
// likelihood calculation Use the same notation as that in the 1973 article (Am
// J Hum Genet 25:471-492)
typedef struct __Independent_Contrasts
{
  phydbl *x;
  phydbl *tprime;
} t_ctrst;

/*!********************************************************/

typedef struct __Continuous_Model
{
  phydbl *mu_down;
  phydbl *var_down;
  phydbl *logrem_down;

  phydbl *mu_up;
  phydbl *var_up;
  phydbl *logrem_up;

  phydbl *lnL_up;
  phydbl *lnL_down;
  phydbl *lnL;

  phydbl combined_lnL;

  phydbl *obs_var; // Variance of observational model (one param for each
                   // spatial dimension)
  phydbl    obs_var_min;
  phydbl    obs_var_max;
  short int obs_model;     // Observational model 0/1
  short int obs_model_est; // Estimate variance of observational model 0/1
  phydbl    obs_var_prior_mean;

  short int *both_sides;

} t_contmod;

/*!********************************************************/

typedef struct __Geo_Coord
{
  phydbl             *lonlat; /* longitude-latitude vector */
  int                 dim;
  char               *id;
  struct __Geo_Coord *cpy; /* keep a copy of this coordinate */

} t_geo_coord;

/*!********************************************************/

typedef struct __Geo_Velocity
{
  phydbl                *deriv; /* derivative of longitude-latitude vector */
  int                    dim;
  char                  *id;
  struct __Geo_Velocity *cpy; /* keep a copy of this velocity */

} t_geo_veloc;

/*!********************************************************/

typedef struct __Lindisk_Node
{
  struct __Disk_Event    *disk;
  struct __Disk_Event    *baseline;
  struct __Lindisk_Node **next;
  struct __Lindisk_Node  *prev;
  struct __Lindisk_Node  *img;
  struct __Geo_Coord     *coord;
  struct __Geo_Coord     *cpy_coord;
  struct __Geo_Velocity  *veloc;
  struct __Node          *nd;

  short int is_hit;

  int n_next;

  phydbl rr;
  phydbl sigsq;
  phydbl base_line;
} t_ldsk;

/*!********************************************************/

typedef struct __Polygon
{
  struct __Geo_Coord **poly_vert;   /* array of polygon vertex coordinates */
  int                  n_poly_vert; /* number of vertices */
} t_poly;

/*!********************************************************/

typedef struct __SampArea
{
  int      n_poly; /* Number of polygons making the sampling area */
  t_poly **a_poly; /* Polygons making the sampling area */
} t_sarea;

/*!********************************************************/

typedef struct __JSON_KeyVal
{
  char                 *key;
  char                 *value;
  struct __JSON_Object *object;
  struct __JSON_Array  *array;
  struct __JSON_KeyVal *next;
} json_kv;

/*!********************************************************/

typedef struct __JSON_Object
{
  struct __JSON_KeyVal *kv;
  struct __JSON_Object *next;
} json_o;

/*!********************************************************/

typedef struct __JSON_Array
{
  struct __JSON_Object *object;
} json_a;

/*!********************************************************/

typedef struct __Label
{
  char           *key;
  char           *val;
  char            sep;
  struct __Label *next;
  struct __Label *prev;
} t_label;

/*!********************************************************/
/*!********************************************************/
/*!********************************************************/

void    Unroot_Tree(char **subtrees);
void    Set_Edge_Dirs(t_edge *b, t_node *a, t_node *d, t_tree *tree);
void    Restrict_To_Coding_Position(align **data, option *io);
void    Uppercase(char *ch);
void    Lowercase(char *ch);
calign *Compact_Data(align **data, option *io);
calign *Compact_Cdata(calign *data, option *io);
void    Traverse_Prefix_Tree(int site, int seqnum, int *patt_num, int *n_patt,
                             align **data, option *io, pnode *n);
pnode  *Create_Pnode(int size);
void    Get_Base_Freqs(calign *data);
void    Get_AA_Freqs(calign *data);
void    Swap_Nodes_On_Edges(t_edge *e1, t_edge *e2, int swap, t_tree *tree);
void    Connect_Edges_To_Nodes_Recur(t_node *a, t_node *d, t_tree *tree);
void    Connect_One_Edge_To_Two_Nodes(t_node *a, t_node *d, t_edge *b,
                                      t_tree *tree);
void    Update_Dirs(t_tree *tree);
void    Exit(char *message);
void   *mCalloc(int nb, size_t size);
void   *mRealloc(void *p, int nb, size_t size);
int     Sort_Phydbl_Decrease(const void *a, const void *b);
void    Qksort_Int(int *A, int *B, int ilo, int ihi);
void    Qksort(phydbl *A, phydbl *B, int ilo, int ihi);
void    Qksort_Matrix(phydbl **A, int col, int ilo, int ihi);
void    Order_Tree_Seq(t_tree *tree, align **data);
char   *Add_Taxa_To_Constraint_Tree(FILE *fp, calign *cdata);
void    Check_Constraint_Tree_Taxa_Names(t_tree *tree, calign *cdata);
void    Order_Tree_CSeq(t_tree *tree, calign *cdata);
void    Init_Mat(matrix *mat, calign *data);
void    Copy_Tax_Names_To_Tip_Labels(t_tree *tree, calign *data);
void    Share_Lk_Struct(t_tree *t_full, t_tree *t_empt);
void    Share_Spr_Struct(t_tree *t_full, t_tree *t_empt);
void    Share_Pars_Struct(t_tree *t_full, t_tree *t_empt);
int     Sort_Edges_NNI_Score(t_tree *tree, t_edge **sorted_edges, int n_elem);
int     Sort_Edges_Depth(t_tree *tree, t_edge **sorted_edges, int n_elem);
void    NNI(t_tree *tree, t_edge *b_fcus, int do_swap);
void    NNI_Pars(t_tree *tree, t_edge *b_fcus, int do_swap);
void    Swap(t_node *a, t_node *b, t_node *c, t_node *d, t_tree *tree);
void    Update_SubTree_Partial_Lk(t_edge *b_fcus, t_node *a, t_node *d,
                                  t_tree *tree);
void    Copy_Seq_Names_To_Tip_Labels(t_tree *tree, calign *data);
calign *Copy_Cseq(calign *ori, option *io, t_tree *tree);
int     Filexists(char *filename);
int     Is_Invar(int patt_num, int stepsize, int datatype, calign *data);
int     Is_Ambigu(char *state, int datatype, int stepsize);
void    Check_Ambiguities(calign *data, int datatype, int stepsize);
int     Get_State_From_Ui(int ui, int datatype);
int     Assign_State(char *c, int datatype, int stepsize);
char    Reciproc_Assign_State(int i_state, int datatype);
int     Assign_State_With_Ambiguity(char *c, int datatype, int stepsize);
void    Clean_Tree_Connections(t_tree *tree);
void    Bootstrap(t_tree *tree);
void    Br_Len_Involving_Invar(t_tree *tree);
void    Br_Len_Not_Involving_Invar(t_tree *tree);
void    Getstring_Stdin(char *s);
phydbl  Num_Derivatives_One_Param(phydbl (*func)(t_tree *tree), t_tree *tree,
                                  phydbl f0, phydbl *param, int which,
                                  int n_param, phydbl stepsize, short int logt,
                                  short int expt, phydbl *err, int precise,
                                  int is_positive);
phydbl  Num_Derivatives_One_Param_Nonaligned(
     phydbl (*func)(t_tree *tree), t_tree *tree, phydbl f0, phydbl **param,
     int which, int n_param, phydbl stepsize, short int logt, short int expt,
     phydbl *err, int precise, int is_positive);
int     Num_Derivative_Several_Param(t_tree *tree, phydbl *param, int n_param,
                                     phydbl stepsize, short int logt,
                                     short int expt, phydbl (*func)(t_tree *tree),
                                     phydbl *derivatives, int is_positive);
int     Num_Derivative_Several_Param_Nonaligned(t_tree *tree, phydbl **param,
                                                int n_param, phydbl stepsize,
                                                short int logt, short int expt,
                                                phydbl (*func)(t_tree *tree),
                                                phydbl *derivatives,
                                                int     is_positive);
int     Compare_Two_States(char *state1, char *state2, int state_size);
void    Copy_One_State(char *from, char *to, int state_size);
void    Copy_Dist(phydbl **cpy, phydbl **orig, int n);
t_mod  *Copy_Model(t_mod *ori);
void    Record_Model(t_mod *ori, t_mod *cpy);
void    Set_Defaults_Input(option *io);
void    Set_Defaults_Model(t_mod *mod);
void    Set_Defaults_Optimiz(t_opt *s_opt);
void    Test_Node_Table_Consistency(t_tree *tree);
void    Get_Bip(t_node *a, t_node *d, t_tree *tree);
void    Alloc_Bip(t_tree *tree);
int     Sort_Phydbl_Increase(const void *a, const void *b);
int     Sort_String(const void *a, const void *b);
phydbl  Compare_Bip(t_tree *tree1, t_tree *tree2, int on_existing_edges_only, int comparison_criterion, phydbl tree2_brlen_cutoff);
void    Compare_Bip_Distance(t_tree *tree1, t_tree *tree2);
void    Match_Tip_Numbers(t_tree *tree1, t_tree *tree2);
void    Test_Multiple_Data_Set_Format(option *io);
int     Are_Compatible(char *statea, char *stateb, int stepsize, int datatype);
void    Hide_Ambiguities(calign *data);
void    Copy_Tree(t_tree *ori, t_tree *cpy);
void    Prune_Subtree(t_node *a, t_node *d, t_edge **target, t_edge **residual,
                      t_tree *tree);
void    Graft_Subtree(t_edge *target, t_node *link, t_node *link_daughter,
                      t_edge *residual, t_node *target_nd, t_tree *tree);
void    Reassign_Node_Nums(t_node *a, t_node *d, unsigned int *curr_ext_node,
                           unsigned int *curr_int_node, t_tree *tree);
void    Reassign_Edge_Nums(t_node *a, t_node *d, int *curr_br, t_tree *tree);
void    Find_Mutual_Direction(t_node *n1, t_node *n2, short int *dir_n1_to_n2,
                              short int *dir_n2_to_n1);
void    Update_Dir_To_Tips(t_node *a, t_node *d, t_tree *tree);
void    Fill_Dir_Table(t_tree *tree);
int     Get_Subtree_Size(t_node *a, t_node *d);
void    Init_Eigen_Struct(eigen *this);
phydbl  Triple_Dist(t_node *a, t_tree *tree);
phydbl  Triple_Dist_Approx(t_node *a, t_edge *b, t_tree *tree);
void    Make_Symmetric(phydbl **F, int size);
void    Divide_Mat_By_Vect(phydbl **F, phydbl *vect, int size);
void    Found_In_Subtree(t_node *a, t_node *d, t_node *target, int *match,
                         t_tree *tree);
void    Get_List_Of_Target_Edges(t_node *a, t_node *d, t_edge **list,
                                 int *list_size, t_tree *tree);
void    Fix_All(t_tree *tree);
void    Record_Br_Len(t_tree *tree);
void    Restore_Br_Len(t_tree *tree);
void    Get_Dist_Btw_Edges(t_node *a, t_node *d, t_tree *tree);
void    Detect_Polytomies(t_edge *b, phydbl l_thresh, t_tree *tree);
void    Get_List_Of_Nodes_In_Polytomy(t_node *a, t_node *d, t_node ***list,
                                      int *size_list);
void    Check_Path(t_node *a, t_node *d, t_node *target, t_tree *tree);
void    Connect_Two_Nodes(t_node *a, t_node *d);
void    Get_List_Of_Adjacent_Targets(t_node *a, t_node *d, t_node ***node_list,
                                     t_edge ***edge_list, int *list_size,
                                     int curr_depth, int max_depth);
void    Sort_List_Of_Adjacent_Targets(t_edge ***list, int list_size);
t_node *Common_Nodes_Btw_Two_Edges(t_edge *a, t_edge *b);
void    Random_Tree(t_tree *tree);
void    Reorganize_Edges_Given_Lk_Struct(t_tree *tree);
void    Random_NNI(int n_moves, t_tree *tree);
void    Fill_Missing_Dist(matrix *mat);
void    Fill_Missing_Dist_XY(int x, int y, matrix *mat);
phydbl  Least_Square_Missing_Dist_XY(int x, int y, phydbl dxy, matrix *mat);
void    Check_Memory_Amount(t_tree *tree);
int     Get_State_From_P_Lk(phydbl *p_lk, int pos, t_tree *tree);
int     Get_State_From_P_Pars(short int *p_pars, int pos, t_tree *tree);
void    Check_Dirs(t_tree *tree);
void    Warn_And_Exit(const char *s);
void    Randomize_Sequence_Order(calign *cdata);
void    Update_Root_Pos(t_tree *tree);
void    Add_Root(t_edge *target, t_tree *tree);
void    Update_Ancestors(t_node *a, t_node *d, t_edge *b, t_tree *tree);
#if (defined PHYTIME || defined SERGEII)
t_tree *Generate_Random_Tree_From_Scratch(int n_otu, int rooted);
#endif
void    Random_Lineage_Rates(t_node *a, t_node *d, t_edge *b, phydbl stick_prob,
                             phydbl *rates, int curr_rate, int n_rates,
                             t_tree *tree);
t_edge *Find_Edge_With_Label(char *label, t_tree *tree);
void    Site_Diversity(t_tree *tree);
void    Site_Diversity_Post(t_node *a, t_node *d, t_edge *b, t_tree *tree);
void    Site_Diversity_Pre(t_node *a, t_node *d, t_edge *b, t_tree *tree);
void    Subtree_Union(t_node *n, t_edge *b_fcus, t_tree *tree);
void    Binary_Decomposition(int value, int *bit_vect, int size);
void    Print_Diversity_Header(FILE *fp, t_tree *tree);
void    Best_Of_NNI_And_SPR(t_tree *tree);
int     Polint(phydbl *xa, phydbl *ya, int n, phydbl x, phydbl *y, phydbl *dy);
t_tree *Dist_And_BioNJ(calign *cdata, t_mod *mod, option *io);
void    Add_BioNJ_Branch_Lengths(t_tree *tree, calign *cdata, t_mod *mod,
                                 matrix *mat);
char   *Bootstrap_From_String(char *s_tree, calign *cdata, t_mod *mod,
                              option *io);
char   *aLRT_From_String(char *s_tree, calign *cdata, t_mod *mod, option *io);
void    Find_Common_Tips(t_tree *tree1, t_tree *tree2);
phydbl  Get_Tree_Size(t_tree *tree);
void    Dist_To_Root_Pre(t_node *a, t_node *d, t_edge *b, t_tree *tree);
void    Dist_To_Root(t_tree *tree);
char   *Basename(char *path);
t_node *Find_Lca_Pair_Of_Nodes(t_node *n1, t_node *n2, int *dist, t_tree *tree);
t_node *Find_Lca_Clade(t_node **node_list, int node_list_size, t_tree *tree);
int     Get_List_Of_Ancestors(t_node *ref_node, t_node **list, int *size,
                              t_tree *tree);
int     Edge_Num_To_Node_Num(int edge_num, t_tree *tree);
void    Branch_Lengths_To_Rate_Lengths(t_tree *tree);
void    Branch_Lengths_To_Rate_Lengths_Pre(t_node *a, t_node *d, t_tree *tree);
int     Find_Clade(char **tax_name_list, int list_size, t_tree *tree);
void    Find_Clade_Pre(t_node *a, t_node *d, int *tax_num_list, int list_size,
                       int *num, t_tree *tree);
t_edge *Find_Root_Edge(FILE *fp_input_tree, t_tree *tree);
void    Copy_Tree_Topology_With_Labels(t_tree *ori, t_tree *cpy);
void    Set_Model_Name(t_mod *mod);
void    Adjust_Min_Diff_Lk(t_tree *tree);
void    Translate_Tax_Names(char **tax_names, t_tree *tree);
void    Skip_Comment(FILE *fp);
void    Get_Best_Root_Position(t_tree *tree);
void    Get_Best_Root_Position_Post(t_node *a, t_node *d, int *has_outgrp,
                                    t_tree *tree);
void    Get_Best_Root_Position_Pre(t_node *a, t_node *d, t_tree *tree);
void    Get_OutIn_Scores(t_node *a, t_node *d);
int     Check_Sequence_Name(char *s);
int     Scale_Subtree_Height(t_node *a, phydbl K, phydbl floor, int *n_nodes,
                             t_tree *tree);
void    Scale_Node_Heights_Post(t_node *a, t_node *d, phydbl K, phydbl floor,
                                int *n_nodes, t_tree *tree);
int     Scale_Subtree_Rates(t_node *a, phydbl mult, int *n_nodes, t_tree *tree);
void    Check_Br_Len_Bounds(t_tree *tree);
int    Scale_Subtree_Rates_Post(t_node *a, t_node *d, phydbl mult, int *n_nodes,
                                t_tree *tree);
void   Get_Node_Ranks(t_tree *tree);
void   Get_Node_Ranks_Pre(t_node *a, t_node *d, t_tree *tree);
void   Log_Br_Len(t_tree *tree);
phydbl Diff_Lk_Norm_At_Given_Edge(t_edge *b, t_tree *tree);
void   Adjust_Variances(t_tree *tree);
phydbl Effective_Sample_Size(phydbl first_val, phydbl last_val, phydbl sum,
                             phydbl sumsq, phydbl sumcurnext, int n);
void   Rescale_Free_Rate_Tree(t_tree *tree);
phydbl Rescale_Br_Len_Multiplier_Tree(t_tree *tree);
phydbl Unscale_Br_Len_Multiplier_Tree(t_tree *tree);
phydbl Reflect(phydbl x, phydbl l, phydbl u);
int    Are_Equal(phydbl a, phydbl b, phydbl eps);
int    Check_Topo_Constraints(t_tree *big_tree, t_tree *small_tree);
void   Prune_Tree(t_tree *big_tree, t_tree *small_tree);
void   Match_Nodes_In_Small_Tree(t_tree *small_tree, t_tree *big_tree);
void   Find_Surviving_Edges_In_Small_Tree(t_tree *small_tree, t_tree *big_tree);
void   Find_Surviving_Edges_In_Small_Tree_Post(t_node *a, t_node *d,
                                               t_tree *small_tree,
                                               t_tree *big_tree);
void   Set_Taxa_Id_Ranking(t_tree *tree);
void   Get_Edge_Binary_Coding_Number(t_tree *tree);
void   Make_Ancestral_Seq(t_tree *tree);
void   Make_MutMap(t_tree *tree);
int    Get_Mutmap_Val(int edge, int site, int mut, t_tree *tree);
void   Get_Mutmap_Coord(int idx, int *edge, int *site, int *mut, t_tree *tree);
void   Copy_Edge_Lengths(t_tree *to, t_tree *from);
void   Init_Scalar_Dbl(scalar_dbl *p);
void   Init_Scalar_Int(scalar_int *p);
void   Init_Vect_Dbl(int len, vect_dbl *p);
void   Init_Vect_Int(int len, vect_int *p);
char  *To_Lower_String(char *in);
phydbl String_To_Dbl(char *string);
int    String_To_Int(char *string);
char  *To_Upper_String(char *in);
void   Connect_CSeqs_To_Nodes(calign *cdata, option *io, t_tree *tree);
void   Joint_Proba_States_Left_Right(phydbl *Pij, phydbl *p_lk_left,
                                     phydbl *p_lk_rght, vect_dbl *pi,
                                     int scale_left, int scale_rght, phydbl *F,
                                     int n, int site, t_tree *tree);
void   Set_Both_Sides(int yesno, t_tree *tree);
void   Set_D_States(calign *data, int datatype, int stepsize);
void   Path_Length(t_node *dep, t_node *arr, phydbl *len, t_tree *tree);
phydbl *Dist_Btw_Tips(t_tree *tree);
void    Best_Root_Position_IL_Model(t_tree *tree);
void    Set_Br_Len_Var(t_edge *b, t_tree *tree);
void    Check_Br_Lens(t_tree *tree);
void    Calculate_Number_Of_Diff_States_Post(t_node *a, t_node *d, t_edge *b,
                                             t_tree *tree);
void    Calculate_Number_Of_Diff_States_Pre(t_node *a, t_node *d, t_edge *b,
                                            t_tree *tree);
void    Calculate_Number_Of_Diff_States_Core(t_node *a, t_node *d, t_edge *b,
                                             t_tree *tree);
void    Calculate_Number_Of_Diff_States(t_tree *tree);
void    Build_Distrib_Number_Of_Diff_States_Under_Model(t_tree *tree);
int     Number_Of_Diff_States_One_Site(int site, t_tree *tree);
void    Number_Of_Diff_States_One_Site_Post(t_node *a, t_node *d, t_edge *b,
                                            int site, t_tree *tree);
int     Number_Of_Diff_States_One_Site_Core(t_node *a, t_node *d, t_edge *b,
                                            int site, t_tree *tree);
phydbl  Get_Lk(t_tree *tree);
phydbl  Get_d2Lk(t_tree *tree);
phydbl  Get_dLk(t_tree *tree);
void    Connect_Edges_To_Nodes_Serial(t_tree *tree);
phydbl  Mean_Identity(calign *data);
phydbl  Pairwise_Identity(int i, int j, calign *data);
phydbl  Fst(int i, int j, calign *data);
phydbl  Nucleotide_Diversity(calign *data);
void    Swap_Partial_Lk(t_edge *a, t_edge *b, int side_a, int side_b,
                        t_tree *tree);
scalar_dbl **Copy_Br_Len_Var(t_tree *mixt_tree);
scalar_dbl **Copy_Br_Len(t_tree *mixt_tree);
void         Transfer_Br_Len_To_Tree(scalar_dbl **bl, t_tree *tree);
void         Copy_Scalar_Dbl(scalar_dbl *from, scalar_dbl *to);
scalar_dbl  *Duplicate_Scalar_Dbl(scalar_dbl *from);
scalar_dbl  *Read_Weights(option *io);
phydbl       Scalar_Elem(int pos, scalar_dbl *scl);
int          Scalar_Len(scalar_dbl *scl);
int          Vect_Len(vect_dbl *scl);
void         Set_Scalar_Dbl(phydbl val, scalar_dbl *from);
void         Set_Scalar_Dbl_Min_Thresh(phydbl thresh, scalar_dbl *from);
void         Set_Scalar_Dbl_Max_Thresh(phydbl thresh, scalar_dbl *from);
void List_Of_Regraft_Nodes(t_node *a, t_node *d, phydbl time_thresh, t_ll *list,
                           t_tree *tree);
void Push_Bottom_Linked_List(void *what, t_ll **list, bool remove_duplicates);
void Remove_From_Linked_List(t_ll *elem, void *val, t_ll **list);
int  Linked_List_Len(t_ll *list);
void   *Linked_List_Elem(int pos, t_ll *ll);
void    Randomize_Tree(t_tree *tree, int n_prune_regraft);
t_ll   *Get_List_Of_Reachable_Tips(t_node *a, t_node *d, t_tree *tree);
void    Get_List_Of_Reachable_Tips_Post(t_node *a, t_node *d, t_ll **list,
                                        t_tree *tree);
phydbl  Length_Of_Path_Between_List_Of_Tips(t_ll *tips0, t_ll *tips1,
                                            matrix *mat);
void    Set_Update_Eigen_Lr(int yn, t_tree *tree);
void    Set_Use_Eigen_Lr(int yn, t_tree *tree);
void    Random_Walk_Along_Tree_On_Radius(t_node *a, t_node *d, t_edge *b,
                                         phydbl *radius, t_edge **target_edge,
                                         t_node **target_nd, phydbl *target_time,
                                         t_tree *tree);
void    Table_Top(unsigned int width);
void    Table_Row(unsigned int width);
void    Table_Bottom(unsigned int width);
t_cal  *Duplicate_Calib(t_cal *from);
t_clad *Duplicate_Clade(t_clad *from);
void    Swap_Partial_Lk_Extra(t_edge *b, t_node *d, int whichone, t_tree *tree);
void    Remove_Duplicates(calign *data, option *io, t_tree *tree);
short int Are_Sequences_Identical(align *seq1, align *seq2);
char     *Mutation_Id(int mut_idx, t_tree *tree);
void      Random_Tax_Idx(t_node *a, t_node *d, int *idx, t_tree *tree);
void      List_Taxa_In_Clade(t_node *a, t_node *d, t_tree *tree);
void      Alias_Subpatt_Pre(t_node *a, t_node *d, t_tree *tree);
void      Alias_Subpatt_Post(t_node *a, t_node *d, t_tree *tree);
void      Alias_One_Subpatt(t_node *a, t_node *d, t_tree *tree);
void      Alias_Subpatt(t_tree *tree);
void   Map_Mutations(t_node *a, t_node *d, int sa, int sd, t_edge *b, int site,
                     int rcat, int *muttype, phydbl *muttime, int *muttax,
                     int *n_mut, t_tree *tree);
void   Set_Update_Eigen(int yesno, t_mod *mod);
int   *Order_Int(const int *u, const int n);
int   *Order_Dbl(const phydbl *u, const int n);
char   Integer_To_IUPAC_Code(int x);
void   Shuffle_Sites(const phydbl prop, align **data, const int n_otu);
void   Multiply_Scalar_Dbl(phydbl mult, scalar_dbl *x);
void   Insert_Duplicates(t_tree *tree);
void   Get_Node_Ranks_From_Dist_To_Root(t_tree *tree);
void   Get_Node_Ranks_From_Times(t_tree *tree);
void   Get_Node_Ranks_From_Tip_Times(t_tree *tree);
phydbl Tree_Height(t_tree *tree);
void   Post_Inflate_Times_To_Get_Reasonnable_Edge_Lengths(t_node *a, t_node *d,
                                                          t_edge *b, phydbl min_l,
                                                          t_tree *tree);
void  Inflate_Times_To_Get_Reasonnable_Edge_Lengths(phydbl min_l, t_tree *tree);
void  Refactor_Tree(t_tree *tree);
void  Refactor_External(t_node *a, t_node *d, int *idx, t_tree *tree);
void  Refactor_Internal(t_node *a, t_node *d, t_edge *b, int *idx_nd,
                        int *idx_br, t_tree *tree);
int  *Integer_To_Bit(int val, const int ns);
char *Bit_To_Character_String(int *bit, int ns);
t_tree *Duplicate_Tree(t_tree *ori);
matrix *K80_dist(calign *data, phydbl g_shape);
matrix *JC69_Dist(calign *data, t_mod *mod);
matrix *Hamming_Dist(calign *data, t_mod *mod);
phydbl  Tree_Length(t_tree *tree);
void    Remove_Duplicates_From_Tree(calign *data, t_tree *tree);
void    Reset_Lk(t_tree *tree);
void    Set_Lk(t_tree *tree);
void    Reset_Prior(t_tree *tree);
void    Set_Prior(t_tree *tree);
t_edge *Get_Edge(t_node *a, t_node *d, t_tree *tree);
void Exchange_Nodes(t_node *a, t_node *d, t_node *w, t_node *v, t_tree *tree);
void Init_T_Beg(t_tree *tree);
void Set_Ignore_Root(int yesno, t_tree *tree);
void Set_Bl_From_Rt(int yesno, t_tree *tree);
void Replace_Short_With_Long_Tax_Names(t_tree *tree, option *io);
void Convert_Lengths_From_Calendar_To_Substitutions(t_tree *tree);
void Convert_Lengths_From_Calendar_To_Substitutions_Post(t_node *a, t_node *d,
                                                         t_tree *tree);
t_label *Get_Next_Label(t_label *curr_lab);
int      Scale_Subtree_Veloc(t_node *a, phydbl mult, int *n_nodes, int dim,
                             t_tree *tree);
int   Scale_Subtree_Veloc_Post(t_node *a, t_node *d, phydbl mult, int *n_nodes,
                               int dim, t_tree *tree);
void  Label_Edges(t_tree *tree);
void  Label_Nodes_With_Velocities(t_tree *tree);
void  Label_Nodes_With_Locations(t_tree *tree);
void  Edge_Labels_To_Rates(t_tree *tree);
void  Node_Labels_To_Velocities(t_tree *tree);
void  Node_Labels_To_Locations(t_tree *tree);
int   Add_Subtree_Veloc(t_node *a, phydbl add, int *n_nodes, int dim,
                        t_tree *tree);
int   Add_Subtree_Veloc_Post(t_node *a, t_node *d, phydbl add, int *n_nodes,
                             int dim, t_tree *tree);
int   Contmod_Start(short int datatype, short int dim_idx, t_tree *tree);
t_ll *Get_Velocity_Targets(t_node *a, t_node *d, t_tree *tree);
void Get_Velocity_Targets_Post(t_node *a, t_node *d, t_ll **list, t_tree *tree);
char D_State_To_Character(int d_state, t_tree *tree);
void ROC(phydbl *probs, short int *truth, int nclasses, int n_obs,
         phydbl *weights, char *tag, FILE *fp);
phydbl AIC(t_tree *tree);
phydbl BIC(t_tree *tree);
void   Set_Edge_Length_Optimizer(t_tree *tree);
int    Number_Of_Free_Params(t_tree *mixt_tree);

#include "alrt.h"
#include "ancestral.h"
#include "bionj.h"
#include "cv.h"
#include "eigen.h"
#include "evolve.h"
#include "free.h"
#include "help.h"
#include "init.h"
#include "io.h"
#include "lk.h"
#include "make.h"
#include "mcmc.h"
#include "models.h"
#include "nexus.h"
#include "optimiz.h"
#include "pars.h"
#include "simu.h"
#include "spr.h"
#include "stats.h"
#include "xml.h"
#include "m4.h"

#ifdef GEO
#include "geo.h"
#endif

#if (defined PHYREX || PHYREXSIM || TEST)
#include "ibm.h"
#include "iou.h"
#include "iwn.h"
#include "location.h"
#include "phyrex.h"
#include "rrw.h"
#include "rw.h"
#include "slfv.h"
#include "velocity.h"
#endif

#ifdef MPI
#include "mpi_boot.h"
#endif

#ifdef MG
#include "mg.h"
#endif

#ifdef TIME
#include "rates.h"
#include "times.h"
#endif

#ifdef _NOT_NEEDED_A_PRIORI
#include "m4.h"
#endif

#if (defined(__AVX__) || defined(__AVX2__))
#include "avx.h"
#elif (defined(__SSE__) || defined(__SSE2__) || defined(__SSE3__) ||           \
       defined(__ARM_NEON))
#include "sse.h"
#endif

#endif
