/*
    Copyright (C) 2015 Tomas Flouri, Diego Darriba

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: Tomas Flouri <Tomas.Flouri@h-its.org>,
    Exelixis Lab, Heidelberg Instutute for Theoretical Studies
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/

#ifndef PLL_H
#define PLL_H
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <ctype.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#if (!defined(__clang__) && defined(__GNUC__) && (__GNUC__ < 4 || \
     (__GNUC__ == 4 && __GNUC_MINOR__ < 7)))
  #if ((__GNUC__ == 4) && (__GNUC_MINOR__ == 6))
    #if (defined(HAVE_AVX2))
      #error "GCC 4.6.x. Please run ./configure --disable-avx2"
    #endif
  #else
    #if (defined(HAVE_AVX2) || defined(HAVE_AVX))
      #error "GCC < 4.6. Please run ./configure --disable-avx --disable-avx2"
    #endif
  #endif
#endif

#ifdef HAVE_X86INTRIN_H
#include <x86intrin.h>
#endif

/* platform specific */

#if (!defined(__APPLE__) && !defined(__WIN32__) && !defined(__WIN64__))
#include <sys/sysinfo.h>
#endif

#if (defined(__WIN32__) || defined(__WIN64__))
#define PLL_EXPORT __declspec(dllexport)
#else
#define PLL_EXPORT
#endif

/* macros */

#define PLL_MIN(a,b) ((a) < (b) ? (a) : (b))
#define PLL_MAX(a,b) ((a) > (b) ? (a) : (b))
#define PLL_SWAP(x,y) do { __typeof__ (x) _t = x; x = y; y = _t; } while(0)
#define PLL_STAT(x) ((pll_hardware.init || pll_hardware_probe()) \
                     && pll_hardware.x)

/* constants */

#define PLL_FAILURE  0
#define PLL_SUCCESS  1

#define PLL_FALSE  0
#define PLL_TRUE   1

#define PLL_ALIGNMENT_CPU   8
#define PLL_ALIGNMENT_SSE  16
#define PLL_ALIGNMENT_AVX  32

#define PLL_LINEALLOC 2048

#define PLL_ASCII_SIZE 256

#define PLL_SCALE_FACTOR 115792089237316195423570985008687907853269984665640564039457584007913129639936.0  /*  2**256 (exactly)  */
#define PLL_SCALE_THRESHOLD (1.0/PLL_SCALE_FACTOR)
#define PLL_SCALE_FACTOR_SQRT 340282366920938463463374607431768211456.0 /* 2**128 */
#define PLL_SCALE_THRESHOLD_SQRT (1.0/PLL_SCALE_FACTOR_SQRT)
#define PLL_SCALE_BUFFER_NONE -1

/* in per-rate scaling mode, maximum difference between scalers
 * please see https://github.com/xflouris/libpll/issues/44  */
#define PLL_SCALE_RATE_MAXDIFF 4

#define PLL_MISC_EPSILON 1e-8
#define PLL_ONE_EPSILON 1e-15
#define PLL_ONE_MIN (1-PLL_ONE_EPSILON)
#define PLL_ONE_MAX (1+PLL_ONE_EPSILON)

/* attribute flags */

#define PLL_ATTRIB_ARCH_CPU            0
#define PLL_ATTRIB_ARCH_SSE       (1 << 0)
#define PLL_ATTRIB_ARCH_AVX       (1 << 1)
#define PLL_ATTRIB_ARCH_AVX2      (1 << 2)
#define PLL_ATTRIB_ARCH_AVX512    (1 << 3)
#define PLL_ATTRIB_ARCH_MASK         0xF

#define PLL_ATTRIB_PATTERN_TIP    (1 << 4)

/* ascertainment bias correction */
#define PLL_ATTRIB_AB_LEWIS        (1 << 5)
#define PLL_ATTRIB_AB_FELSENSTEIN  (2 << 5)
#define PLL_ATTRIB_AB_STAMATAKIS   (3 << 5)
#define PLL_ATTRIB_AB_MASK         (7 << 5)
#define PLL_ATTRIB_AB_FLAG         (1 << 8)

#define PLL_ATTRIB_RATE_SCALERS    (1 << 9)

/* topological rearrangements */

#define PLL_UTREE_MOVE_SPR                  1
#define PLL_UTREE_MOVE_NNI                  2

#define PLL_UTREE_MOVE_NNI_LEFT             1
#define PLL_UTREE_MOVE_NNI_RIGHT            2

#define PLL_TREE_TRAVERSE_POSTORDER         1
#define PLL_TREE_TRAVERSE_PREORDER          2

/* error codes */

#define PLL_ERROR_FILE_OPEN                100
#define PLL_ERROR_FILE_SEEK                101
#define PLL_ERROR_FILE_EOF                 102
#define PLL_ERROR_FASTA_ILLEGALCHAR        103
#define PLL_ERROR_FASTA_UNPRINTABLECHAR    104
#define PLL_ERROR_FASTA_INVALIDHEADER      105
#define PLL_ERROR_PHYLIP_SYNTAX            106
#define PLL_ERROR_PHYLIP_LONGSEQ           107
#define PLL_ERROR_PHYLIP_NONALIGNED        108
#define PLL_ERROR_PHYLIP_ILLEGALCHAR       109
#define PLL_ERROR_PHYLIP_UNPRINTABLECHAR   110
#define PLL_ERROR_NEWICK_SYNTAX            111
#define PLL_ERROR_MEM_ALLOC                112
#define PLL_ERROR_PARAM_INVALID            113
#define PLL_ERROR_TIPDATA_ILLEGALSTATE     114
#define PLL_ERROR_TIPDATA_ILLEGALFUNCTION  115
#define PLL_ERROR_TREE_CONVERSION          116
#define PLL_ERROR_INVAR_INCOMPAT           117
#define PLL_ERROR_INVAR_PROPORTION         118
#define PLL_ERROR_INVAR_PARAMINDEX         119
#define PLL_ERROR_INVAR_NONEFOUND          120
#define PLL_ERROR_AB_INVALIDMETHOD         121
#define PLL_ERROR_AB_NOSUPPORT             122
#define PLL_ERROR_SPR_TERMINALBRANCH       123
#define PLL_ERROR_SPR_NOCHANGE             124
#define PLL_ERROR_NNI_INVALIDMOVE          125
#define PLL_ERROR_NNI_TERMINALBRANCH       126
#define PLL_ERROR_STEPWISE_STRUCT          127
#define PLL_ERROR_STEPWISE_TIPS            128
#define PLL_ERROR_STEPWISE_UNSUPPORTED     129
#define PLL_ERROR_EINVAL                   130

/* utree specific */

#define PLL_UTREE_SHOW_LABEL             (1 << 0)
#define PLL_UTREE_SHOW_BRANCH_LENGTH     (1 << 1)
#define PLL_UTREE_SHOW_CLV_INDEX         (1 << 2)
#define PLL_UTREE_SHOW_SCALER_INDEX      (1 << 3)
#define PLL_UTREE_SHOW_PMATRIX_INDEX     (1 << 4)

/* GAMMA discretization modes */
#define PLL_GAMMA_RATES_MEAN             0
#define PLL_GAMMA_RATES_MEDIAN           1

/* structures and data types */

typedef struct pll_hardware_s
{
  int init;
  /* cpu features */
  int altivec_present;
  int mmx_present;
  int sse_present;
  int sse2_present;
  int sse3_present;
  int ssse3_present;
  int sse41_present;
  int sse42_present;
  int popcnt_present;
  int avx_present;
  int avx2_present;

  /* TODO: add chip,core,mem info */
} pll_hardware_t;

typedef struct pll_partition
{
  unsigned int tips;
  unsigned int clv_buffers;
  unsigned int states;
  unsigned int sites;
  unsigned int pattern_weight_sum;
  unsigned int rate_matrices;
  unsigned int prob_matrices;
  unsigned int rate_cats;
  unsigned int scale_buffers;
  unsigned int attributes;

  const unsigned int * map;

  /* vectorization options */
  size_t alignment;
  unsigned int states_padded;

  double ** clv;
  double ** pmatrix;
  double * rates;
  double * rate_weights;
  double ** subst_params;
  unsigned int ** scale_buffer;
  double ** frequencies;
  double * prop_invar;
  int * invariant;
  unsigned int * pattern_weights;

  int * eigen_decomp_valid;
  double ** eigenvecs;
  double ** inv_eigenvecs;
  double ** eigenvals;

  /* tip-tip precomputation data */
  unsigned int maxstates;
  unsigned char ** tipchars;
  unsigned char * charmap;
  double * ttlookup;
  unsigned int * tipmap;

  /* ascertainment bias correction */
  int asc_bias_alloc;
} pll_partition_t;


/* Structure for driving likelihood operations */

typedef struct pll_operation
{
  unsigned int parent_clv_index;
  int parent_scaler_index;
  unsigned int child1_clv_index;
  unsigned int child1_matrix_index;
  int child1_scaler_index;
  unsigned int child2_clv_index;
  unsigned int child2_matrix_index;
  int child2_scaler_index;
} pll_operation_t;

/* Doubly-linked list */

typedef struct pll_dlist
{
  struct pll_dlist * next;
  struct pll_dlist * prev;
  void * data;
} pll_dlist_t;

/* multiple sequence alignment */
typedef struct pll_msa_s
{
  int count;
  int length;

  char ** sequence;
  char ** label;
} pll_msa_t;

/* Simple structure for handling FASTA parsing */

typedef struct pll_fasta
{
  FILE * fp;
  char line[PLL_LINEALLOC];
  const unsigned int * chrstatus;
  long no;
  long filesize;
  long lineno;
  long stripped_count;
  long stripped[256];
} pll_fasta_t;

/* Simple structure for handling PHYLIP parsing */
typedef struct pll_phylip_s
{
  FILE * fp;
  char * line;
  size_t line_size;
  size_t line_maxsize;
  char buffer[PLL_LINEALLOC];
  const unsigned int * chrstatus;
  long no;
  long filesize;
  long lineno;
  long stripped_count;
  long stripped[256];
} pll_phylip_t;

/* Simple unrooted and rooted tree structure for parsing newick */

typedef struct pll_unode_s
{
  char * label;
  double length;
  unsigned int node_index;
  unsigned int clv_index;
  int scaler_index;
  unsigned int pmatrix_index;
  struct pll_unode_s * next;
  struct pll_unode_s * back;

  void * data;
} pll_unode_t;

typedef struct pll_utree_s
{
  unsigned int tip_count;
  unsigned int inner_count;
  unsigned int edge_count;

  pll_unode_t ** nodes;

} pll_utree_t;

typedef struct pll_rnode_s
{
  char * label;
  double length;
  unsigned int node_index;
  unsigned int clv_index;
  int scaler_index;
  unsigned int pmatrix_index;
  struct pll_rnode_s * left;
  struct pll_rnode_s * right;
  struct pll_rnode_s * parent;

  void * data;
} pll_rnode_t;

typedef struct pll_rtree_s
{
  unsigned int tip_count;
  unsigned int inner_count;
  unsigned int edge_count;

  pll_rnode_t ** nodes;

  pll_rnode_t * root;

} pll_rtree_t;

/* structures for handling topological rearrangement move rollbacks */

typedef struct pll_utree_rb_s
{
  int move_type;
  union
  {
    struct
    {
      pll_unode_t * p;
      pll_unode_t * r;
      pll_unode_t * rb;
      pll_unode_t * pnb;
      pll_unode_t * pnnb;
      double r_len;
      double pnb_len;
      double pnnb_len;
    } spr;
    struct
    {
      pll_unode_t * p;
      int nni_type;
    } nni;
  };
} pll_utree_rb_t;

/* structures for parsimony */

typedef struct pll_parsimony_s
{
  /* common information */
  unsigned int tips;
  unsigned int inner_nodes;
  unsigned int sites;
  unsigned int states;
  unsigned int attributes;
  size_t alignment;

  /* fast unweighted parsimony */
  unsigned int ** packedvector;
  unsigned int * node_cost;
  unsigned int packedvector_count;
  unsigned int const_cost;
  int * informative;
  unsigned int informative_count;

  /* weighted parsimony */
  unsigned int score_buffers;
  unsigned int ancestral_buffers;
  double * score_matrix;
  double ** sbuffer;
  unsigned int ** anc_states;
} pll_parsimony_t;


typedef struct pll_pars_buildop_s
{
  unsigned int parent_score_index;
  unsigned int child1_score_index;
  unsigned int child2_score_index;
} pll_pars_buildop_t;

typedef struct pll_pars_recop_s
{
  unsigned int node_score_index;
  unsigned int node_ancestral_index;
  unsigned int parent_score_index;
  unsigned int parent_ancestral_index;
} pll_pars_recop_t;

/* structures for SVG visualization */

typedef struct pll_svg_attrib_s
{
  int precision;
  long width;
  long font_size;
  long tip_spacing;
  long stroke_width;
  long legend_show;
  long legend_spacing;
  long margin_left;
  long margin_right;
  long margin_bottom;
  long margin_top;
  long node_radius;
  double legend_ratio;
} pll_svg_attrib_t;

/* Reentrant versions of the `random' family of functions.
   These functions all use the following data structure to contain
   state, rather than global state variables. Taken and modified from
   glibc 2.23 */

struct pll_random_data
{
  int32_t *fptr;        /* Front pointer.  */
  int32_t *rptr;        /* Rear pointer.  */
  int32_t *state;       /* Array of state values.  */
  int rand_type;        /* Type of random number generator.  */
  int rand_deg;         /* Degree of random number generator.  */
  int rand_sep;         /* Distance between front and rear.  */
  int32_t *end_ptr;     /* Pointer behind state table.  */
};

/* common data */

PLL_EXPORT extern __thread int pll_errno;
PLL_EXPORT extern __thread char pll_errmsg[200];
PLL_EXPORT extern pll_hardware_t pll_hardware;

PLL_EXPORT extern const unsigned int pll_map_bin[256];
PLL_EXPORT extern const unsigned int pll_map_nt[256];
PLL_EXPORT extern const unsigned int pll_map_aa[256];
PLL_EXPORT extern const unsigned int pll_map_fasta[256];
PLL_EXPORT extern const unsigned int pll_map_phylip[256];

PLL_EXPORT extern const double pll_aa_rates_dayhoff[190];
PLL_EXPORT extern const double pll_aa_rates_lg[190];
PLL_EXPORT extern const double pll_aa_rates_dcmut[190];
PLL_EXPORT extern const double pll_aa_rates_jtt[190];
PLL_EXPORT extern const double pll_aa_rates_mtrev[190];
PLL_EXPORT extern const double pll_aa_rates_wag[190];
PLL_EXPORT extern const double pll_aa_rates_rtrev[190];
PLL_EXPORT extern const double pll_aa_rates_cprev[190];
PLL_EXPORT extern const double pll_aa_rates_vt[190];
PLL_EXPORT extern const double pll_aa_rates_blosum62[190];
PLL_EXPORT extern const double pll_aa_rates_mtmam[190];
PLL_EXPORT extern const double pll_aa_rates_mtart[190];
PLL_EXPORT extern const double pll_aa_rates_mtzoa[190];
PLL_EXPORT extern const double pll_aa_rates_pmb[190];
PLL_EXPORT extern const double pll_aa_rates_hivb[190];
PLL_EXPORT extern const double pll_aa_rates_hivw[190];
PLL_EXPORT extern const double pll_aa_rates_jttdcmut[190];
PLL_EXPORT extern const double pll_aa_rates_flu[190];
PLL_EXPORT extern const double pll_aa_rates_stmtrev[190];
PLL_EXPORT extern const double pll_aa_rates_lg4m[4][190];
PLL_EXPORT extern const double pll_aa_rates_lg4x[4][190];

PLL_EXPORT extern const double pll_aa_freqs_dayhoff[20];
PLL_EXPORT extern const double pll_aa_freqs_lg[20];
PLL_EXPORT extern const double pll_aa_freqs_dcmut[20];
PLL_EXPORT extern const double pll_aa_freqs_jtt[20];
PLL_EXPORT extern const double pll_aa_freqs_mtrev[20];
PLL_EXPORT extern const double pll_aa_freqs_wag[20];
PLL_EXPORT extern const double pll_aa_freqs_rtrev[20];
PLL_EXPORT extern const double pll_aa_freqs_cprev[20];
PLL_EXPORT extern const double pll_aa_freqs_vt[20];
PLL_EXPORT extern const double pll_aa_freqs_blosum62[20];
PLL_EXPORT extern const double pll_aa_freqs_mtmam[20];
PLL_EXPORT extern const double pll_aa_freqs_mtart[20];
PLL_EXPORT extern const double pll_aa_freqs_mtzoa[20];
PLL_EXPORT extern const double pll_aa_freqs_pmb[20];
PLL_EXPORT extern const double pll_aa_freqs_hivb[20];
PLL_EXPORT extern const double pll_aa_freqs_hivw[20];
PLL_EXPORT extern const double pll_aa_freqs_jttdcmut[20];
PLL_EXPORT extern const double pll_aa_freqs_flu[20];
PLL_EXPORT extern const double pll_aa_freqs_stmtrev[20];
PLL_EXPORT extern const double pll_aa_freqs_lg4m[4][20];
PLL_EXPORT extern const double pll_aa_freqs_lg4x[4][20];

#ifdef __cplusplus
extern "C" {
#endif

/* functions in pll.c */

PLL_EXPORT pll_partition_t * pll_partition_create(unsigned int tips,
                                                  unsigned int clv_buffers,
                                                  unsigned int states,
                                                  unsigned int sites,
                                                  unsigned int rate_matrices,
                                                  unsigned int prob_matrices,
                                                  unsigned int rate_cats,
                                                  unsigned int scale_buffers,
                                                  unsigned int attributes);

PLL_EXPORT void pll_partition_destroy(pll_partition_t * partition);

PLL_EXPORT int pll_set_tip_states(pll_partition_t * partition,
                                  unsigned int tip_index,
                                  const unsigned int * map,
                                  const char * sequence);

PLL_EXPORT int pll_set_tip_clv(pll_partition_t * partition,
                               unsigned int tip_index,
                               const double * clv,
                               int padding);

PLL_EXPORT void pll_set_pattern_weights(pll_partition_t * partition,
                                        const unsigned int * pattern_weights);

PLL_EXPORT int pll_set_asc_bias_type(pll_partition_t * partition,
                                     int asc_bias_type);

PLL_EXPORT void pll_set_asc_state_weights(pll_partition_t * partition,
                                          const unsigned int * state_weights);

/* functions in list.c */

PLL_EXPORT int pll_dlist_append(pll_dlist_t ** dlist, void * data);
PLL_EXPORT int pll_dlist_remove(pll_dlist_t ** dlist, void * data);
PLL_EXPORT int pll_dlist_prepend(pll_dlist_t ** dlist, void * data);

/* functions in models.c */

PLL_EXPORT void pll_set_subst_params(pll_partition_t * partition,
                                     unsigned int params_index,
                                     const double * params);

PLL_EXPORT void pll_set_frequencies(pll_partition_t * partition,
                                    unsigned int params_index,
                                    const double * frequencies);

PLL_EXPORT void pll_set_category_rates(pll_partition_t * partition,
                                       const double * rates);

PLL_EXPORT void pll_set_category_weights(pll_partition_t * partition,
                                         const double * rate_weights);

PLL_EXPORT int pll_update_eigen(pll_partition_t * partition,
                                unsigned int params_index);

PLL_EXPORT int pll_update_prob_matrices(pll_partition_t * partition,
                                        const unsigned int * params_index,
                                        const unsigned int * matrix_indices,
                                        const double * branch_lengths,
                                        unsigned int count);

PLL_EXPORT unsigned int pll_count_invariant_sites(pll_partition_t * partition,
                                                  unsigned int * state_inv_count);

PLL_EXPORT int pll_update_invariant_sites(pll_partition_t * partition);

PLL_EXPORT int pll_update_invariant_sites_proportion(pll_partition_t * partition,
                                                     unsigned int params_index,
                                                     double prop_invar);

PLL_EXPORT void * pll_aligned_alloc(size_t size, size_t alignment);

PLL_EXPORT void pll_aligned_free(void * ptr);

/* functions in likelihood.c */

PLL_EXPORT double pll_compute_root_loglikelihood(pll_partition_t * partition,
                                                 unsigned int clv_index,
                                                 int scaler_index,
                                                 const unsigned int * freqs_indices,
                                                 double * persite_lnl);

PLL_EXPORT double pll_compute_edge_loglikelihood(pll_partition_t * partition,
                                                 unsigned int parent_clv_index,
                                                 int parent_scaler_index,
                                                 unsigned int child_clv_index,
                                                 int child_scaler_index,
                                                 unsigned int matrix_index,
                                                 const unsigned int * freqs_indices,
                                                 double * persite_lnl);

/* functions in partials.c */

PLL_EXPORT void pll_update_partials(pll_partition_t * partition,
                                    const pll_operation_t * operations,
                                    unsigned int count);


/* functions in derivatives.c */

PLL_EXPORT int pll_update_sumtable(pll_partition_t * partition,
                                      unsigned int parent_clv_index,
                                      unsigned int child_clv_index,
                                      int parent_scaler_index,
                                      int child_scaler_index,
                                      const unsigned int * params_indices,
                                      double *sumtable);

PLL_EXPORT int pll_compute_likelihood_derivatives(pll_partition_t * partition,
                                                  int parent_scaler_index,
                                                  int child_scaler_index,
                                                  double branch_length,
                                                  const unsigned int * params_indices,
                                                  const double * sumtable,
                                                  double * d_f,
                                                  double * dd_f);

/* functions in gamma.c */

PLL_EXPORT int pll_compute_gamma_cats(double alpha,
                                      unsigned int categories,
                                      double * output_rates,
                                      int rates_mode);

/* functions in output.c */

PLL_EXPORT void pll_show_pmatrix(const pll_partition_t * partition,
                                 unsigned int index,
                                 unsigned int float_precision);

PLL_EXPORT void pll_show_clv(const pll_partition_t * partition,
                             unsigned int clv_index,
                             int scaler_index,
                             unsigned int float_precision);

/* functions in fasta.c */

PLL_EXPORT pll_fasta_t * pll_fasta_open(const char * filename,
                                        const unsigned int * map);

PLL_EXPORT int pll_fasta_getnext(pll_fasta_t * fd, char ** head,
                                 long * head_len,  char ** seq,
                                 long * seq_len, long * seqno);

PLL_EXPORT void pll_fasta_close(pll_fasta_t * fd);

PLL_EXPORT long pll_fasta_getfilesize(const pll_fasta_t * fd);

PLL_EXPORT long pll_fasta_getfilepos(pll_fasta_t * fd);

PLL_EXPORT int pll_fasta_rewind(pll_fasta_t * fd);

/* functions in parse_rtree.y */

PLL_EXPORT pll_rtree_t * pll_rtree_parse_newick(const char * filename);

PLL_EXPORT pll_rtree_t * pll_rtree_parse_newick_string(const char * s);

PLL_EXPORT void pll_rtree_destroy(pll_rtree_t * root,
                                  void (*cb_destroy)(void *));

PLL_EXPORT void pll_rtree_reset_template_indices(pll_rnode_t * node,
                                                 unsigned int tip_count);

PLL_EXPORT void pll_rtree_graph_destroy(pll_rnode_t * root,
                                        void (*cb_destroy)(void *));

PLL_EXPORT pll_rtree_t * pll_rtree_wraptree(pll_rnode_t * root,
                                            unsigned int tip_count);
/* functions in parse_utree.y */

PLL_EXPORT pll_utree_t * pll_utree_parse_newick(const char * filename);

PLL_EXPORT pll_utree_t * pll_utree_parse_newick_string(const char * s);

PLL_EXPORT void pll_utree_destroy(pll_utree_t * tree,
                                  void (*cb_destroy)(void *));

PLL_EXPORT void pll_utree_reset_template_indices(pll_unode_t * node,
                                                 unsigned int tip_count);

PLL_EXPORT void pll_utree_graph_destroy(pll_unode_t * root,
                                        void (*cb_destroy)(void *));

PLL_EXPORT pll_utree_t * pll_utree_wraptree(pll_unode_t * root,
                                            unsigned int tip_count);

/* functions in utree.c */

PLL_EXPORT void pll_utree_show_ascii(const pll_unode_t * tree, int options);

PLL_EXPORT char * pll_utree_export_newick(const pll_unode_t * root,
                                   char * (*cb_serialize)(const pll_unode_t *));

PLL_EXPORT int pll_utree_traverse(pll_unode_t * root,
                                  int traversal,
                                  int (*cbtrav)(pll_unode_t *),
                                  pll_unode_t ** outbuffer,
                                  unsigned int * trav_size);

#if 0
PLL_EXPORT unsigned int pll_utree_query_tipnodes(pll_utree_t * root,
                                                 pll_utree_t ** node_list);

PLL_EXPORT unsigned int pll_utree_query_innernodes(pll_utree_t * root,
                                                   pll_utree_t ** node_list);
#endif

PLL_EXPORT void pll_utree_create_operations(pll_unode_t * const* trav_buffer,
                                            unsigned int trav_buffer_size,
                                            double * branches,
                                            unsigned int * pmatrix_indices,
                                            pll_operation_t * ops,
                                            unsigned int * matrix_count,
                                            unsigned int * ops_count);

PLL_EXPORT int pll_utree_check_integrity(const pll_utree_t * root);

PLL_EXPORT pll_unode_t * pll_utree_graph_clone(const pll_unode_t * root);

PLL_EXPORT pll_utree_t * pll_utree_clone(const pll_utree_t * root);

PLL_EXPORT pll_utree_t * pll_rtree_unroot(pll_rtree_t * tree);

PLL_EXPORT int pll_utree_every(pll_utree_t * tree,
                               int (*cb)(pll_unode_t *));

PLL_EXPORT int pll_utree_every_const(const pll_utree_t * tree,
                                     int (*cb)(const pll_unode_t *));

PLL_EXPORT void pll_utree_create_pars_buildops(pll_unode_t * const* trav_buffer,
                                               unsigned int trav_buffer_size,
                                               pll_pars_buildop_t * ops,
                                               unsigned int * ops_count);

/* functions in phylip.c */

PLL_EXPORT void pll_msa_destroy(pll_msa_t * msa);

PLL_EXPORT pll_phylip_t * pll_phylip_open(const char * filename,
                                          const unsigned int * map);

PLL_EXPORT int pll_phylip_rewind(pll_phylip_t * fd);

PLL_EXPORT void pll_phylip_close(pll_phylip_t * fd);

PLL_EXPORT pll_msa_t * pll_phylip_parse_interleaved(pll_phylip_t * fd);

PLL_EXPORT pll_msa_t * pll_phylip_parse_sequential(pll_phylip_t * fd);

/* functions in rtree.c */

PLL_EXPORT void pll_rtree_show_ascii(const pll_rnode_t * root, int options);

PLL_EXPORT char * pll_rtree_export_newick(const pll_rnode_t * root,
                                   char * (*cb_serialize)(const pll_rnode_t *));

PLL_EXPORT int pll_rtree_traverse(pll_rnode_t * root,
                                  int traversal,
                                  int (*cbtrav)(pll_rnode_t *),
                                  pll_rnode_t ** outbuffer,
                                  unsigned int * trav_size);

#if 0
PLL_EXPORT unsigned int pll_rtree_query_tipnodes(pll_rtree_t * root,
                                                 pll_rtree_t ** node_list);

PLL_EXPORT unsigned int pll_rtree_query_innernodes(pll_rtree_t * root,
                                                   pll_rtree_t ** node_list);
#endif

PLL_EXPORT void pll_rtree_create_operations(pll_rnode_t * const* trav_buffer,
                                            unsigned int trav_buffer_size,
                                            double * branches,
                                            unsigned int * pmatrix_indices,
                                            pll_operation_t * ops,
                                            unsigned int * matrix_count,
                                            unsigned int * ops_count);

#if 0
PLL_EXPORT int pll_rtree_traverse_preorder(pll_rtree_t * root,
                                           int (*cbtrav)(pll_rtree_t *),
                                           pll_rtree_t ** outbuffer,
                                           unsigned int * trav_size);
#endif

PLL_EXPORT void pll_rtree_create_pars_buildops(pll_rnode_t * const* trav_buffer,
                                               unsigned int trav_buffer_size,
                                               pll_pars_buildop_t * ops,
                                               unsigned int * ops_count);

PLL_EXPORT void pll_rtree_create_pars_recops(pll_rnode_t * const* trav_buffer,
                                             unsigned int trav_buffer_size,
                                             pll_pars_recop_t * ops,
                                             unsigned int * ops_count);

/* functions in core_partials.c */

PLL_EXPORT void pll_core_create_lookup(unsigned int states,
                                       unsigned int rate_cats,
                                       double * lookup,
                                       const double * left_matrix,
                                       const double * right_matrix,
                                       const unsigned int * tipmap,
                                       unsigned int tipmap_size,
                                       unsigned int attrib);

PLL_EXPORT void pll_core_update_partial_tt(unsigned int states,
                                           unsigned int sites,
                                           unsigned int rate_cats,
                                           double * parent_clv,
                                           unsigned int * parent_scaler,
                                           const unsigned char * left_tipchars,
                                           const unsigned char * right_tipchars,
                                           const unsigned int * tipmap,
                                           unsigned int tipmap_size,
                                           const double * lookup,
                                           unsigned int attrib);

PLL_EXPORT void pll_core_update_partial_ti(unsigned int states,
                                           unsigned int sites,
                                           unsigned int rate_cats,
                                           double * parent_clv,
                                           unsigned int * parent_scaler,
                                           const unsigned char * left_tipchars,
                                           const double * right_clv,
                                           const double * left_matrix,
                                           const double * right_matrix,
                                           const unsigned int * right_scaler,
                                           const unsigned int * tipmap,
                                           unsigned int tipmap_size,
                                           unsigned int attrib);

PLL_EXPORT void pll_core_update_partial_ii(unsigned int states,
                                           unsigned int sites,
                                           unsigned int rate_cats,
                                           double * parent_clv,
                                           unsigned int * parent_scaler,
                                           const double * left_clv,
                                           const double * right_clv,
                                           const double * left_matrix,
                                           const double * right_matrix,
                                           const unsigned int * left_scaler,
                                           const unsigned int * right_scaler,
                                           unsigned int attrib);

PLL_EXPORT void pll_core_create_lookup_4x4(unsigned int rate_cats,
                                           double * lookup,
                                           const double * left_matrix,
                                           const double * right_matrix);

PLL_EXPORT void pll_core_update_partial_tt_4x4(unsigned int sites,
                                               unsigned int rate_cats,
                                               double * parent_clv,
                                               unsigned int * parent_scaler,
                                               const unsigned char * left_tipchars,
                                               const unsigned char * right_tipchars,
                                               const double * lookup,
                                               unsigned int attrib);

PLL_EXPORT void pll_core_update_partial_ti_4x4(unsigned int sites,
                                               unsigned int rate_cats,
                                               double * parent_clv,
                                               unsigned int * parent_scaler,
                                               const unsigned char * left_tipchars,
                                               const double * right_clv,
                                               const double * left_matrix,
                                               const double * right_matrix,
                                               const unsigned int * right_scaler,
                                               unsigned int attrib);

/* functions in core_derivatives.c */

PLL_EXPORT int pll_core_update_sumtable_ti_4x4(unsigned int sites,
                                               unsigned int rate_cats,
                                               const double * parent_clv,
                                               const unsigned char * left_tipchars,
                                               const unsigned int * parent_scaler,
                                               double * const * eigenvecs,
                                               double * const * inv_eigenvecs,
                                               double * const * freqs,
                                               const unsigned int * tipmap,
                                               double * sumtable,
                                               unsigned int attrib);

PLL_EXPORT int pll_core_update_sumtable_ii(unsigned int states,
                                           unsigned int sites,
                                           unsigned int rate_cats,
                                           const double * parent_clv,
                                           const double * child_clv,
                                           const unsigned int * parent_scaler,
                                           const unsigned int * child_scaler,
                                           double * const * eigenvecs,
                                           double * const * inv_eigenvecs,
                                           double * const * freqs,
                                           double * sumtable,
                                           unsigned int attrib);

PLL_EXPORT int pll_core_update_sumtable_ti(unsigned int states,
                                           unsigned int sites,
                                           unsigned int rate_cats,
                                           const double * parent_clv,
                                           const unsigned char * left_tipchars,
                                           const unsigned int * parent_scaler,
                                           double * const * eigenvecs,
                                           double * const * inv_eigenvecs,
                                           double * const * freqs,
                                           const unsigned int * tipmap,
                                           unsigned int tipmap_size,
                                           double * sumtable,
                                           unsigned int attrib);

PLL_EXPORT int pll_core_likelihood_derivatives(unsigned int states,
                                               unsigned int sites,
                                               unsigned int rate_cats,
                                               const double * rate_weights,
                                               const unsigned int * parent_scaler,
                                               const unsigned int * child_scaler,
                                               const int * invariant,
                                               const unsigned int * pattern_weights,
                                               double branch_length,
                                               const double * prop_invar,
                                               double * const * freqs,
                                               const double * rates,
                                               double * const * eigenvals,
                                               const double * sumtable,
                                               double * d_f,
                                               double * dd_f,
                                               unsigned int attrib);

/* functions in core_likelihood.c */

PLL_EXPORT double pll_core_edge_loglikelihood_ii(unsigned int states,
                                                 unsigned int sites,
                                                 unsigned int rate_cats,
                                                 const double * parent_clv,
                                                 const unsigned int * parent_scaler,
                                                 const double * child_clv,
                                                 const unsigned int * child_scaler,
                                                 const double * pmatrix,
                                                 double * const * frequencies,
                                                 const double * rate_weights,
                                                 const unsigned int * pattern_weights,
                                                 const double * invar_proportion,
                                                 const int * invar_indices,
                                                 const unsigned int * freqs_indices,
                                                 double * persite_lnl,
                                                 unsigned int attrib);

PLL_EXPORT double pll_core_edge_loglikelihood_ti(unsigned int states,
                                                 unsigned int sites,
                                                 unsigned int rate_cats,
                                                 const double * parent_clv,
                                                 const unsigned int * parent_scaler,
                                                 const unsigned char * tipchars,
                                                 const unsigned int * tipmap,
                                                 unsigned int tipmap_size,
                                                 const double * pmatrix,
                                                 double * const * frequencies,
                                                 const double * rate_weights,
                                                 const unsigned int * pattern_weights,
                                                 const double * invar_proportion,
                                                 const int * invar_indices,
                                                 const unsigned int * freqs_indices,
                                                 double * persite_lnl,
                                                 unsigned int attrib);

PLL_EXPORT double pll_core_edge_loglikelihood_ti_4x4(unsigned int sites,
                                                     unsigned int rate_cats,
                                                     const double * parent_clv,
                                                     const unsigned int * parent_scaler,
                                                     const unsigned char * tipchars,
                                                     const double * pmatrix,
                                                     double * const * frequencies,
                                                     const double * rate_weights,
                                                     const unsigned int * pattern_weights,
                                                     const double * invar_proportion,
                                                     const int * invar_indices,
                                                     const unsigned int * freqs_indices,
                                                     double * persite_lnl,
                                                     unsigned int attrib);

PLL_EXPORT double pll_core_root_loglikelihood(unsigned int states,
                                              unsigned int sites,
                                              unsigned int rate_cats,
                                              const double * clv,
                                              const unsigned int * scaler,
                                              double * const * frequencies,
                                              const double * rate_weights,
                                              const unsigned int * pattern_weights,
                                              const double * invar_proportion,
                                              const int * invar_indices,
                                              const unsigned int * freqs_indices,
                                              double * persite_lnl,
                                              unsigned int attrib);

/* functions in core_partials_sse.c */

#ifdef HAVE_SSE3
PLL_EXPORT void pll_core_create_lookup_sse(unsigned int states,
                                           unsigned int rate_cats,
                                           double * ttlookup,
                                           const double * left_matrix,
                                           const double * right_matrix,
                                           const unsigned int * tipmap,
                                           unsigned int tipmap_size);

PLL_EXPORT void pll_core_create_lookup_4x4_sse(unsigned int rate_cats,
                                               double * lookup,
                                               const double * left_matrix,
                                               const double * right_matrix);

PLL_EXPORT void pll_core_update_partial_tt_sse(unsigned int states,
                                               unsigned int sites,
                                               unsigned int rate_cats,
                                               double * parent_clv,
                                               unsigned int * parent_scaler,
                                               const unsigned char * left_tipchars,
                                               const unsigned char * right_tipchars,
                                               const double * lookup,
                                               unsigned int tipstates_count,
                                               unsigned int attrib);

PLL_EXPORT void pll_core_update_partial_tt_4x4_sse(unsigned int sites,
                                                   unsigned int rate_cats,
                                                   double * parent_clv,
                                                   unsigned int * parent_scaler,
                                                   const unsigned char * left_tipchars,
                                                   const unsigned char * right_tipchars,
                                                   const double * lookup,
                                                   unsigned int attrib);

PLL_EXPORT void pll_core_update_partial_ti_sse(unsigned int states,
                                               unsigned int sites,
                                               unsigned int rate_cats,
                                               double * parent_clv,
                                               unsigned int * parent_scaler,
                                               const unsigned char * left_tipchars,
                                               const double * right_clv,
                                               const double * left_matrix,
                                               const double * right_matrix,
                                               const unsigned int * right_scaler,
                                               const unsigned int * tipmap,
                                               unsigned int tipmap_size,
                                               unsigned int attrib);


PLL_EXPORT void pll_core_update_partial_ti_4x4_sse(unsigned int sites,
                                                   unsigned int rate_cats,
                                                   double * parent_clv,
                                                   unsigned int * parent_scaler,
                                                   const unsigned char * left_tipchar,
                                                   const double * right_clv,
                                                   const double * left_matrix,
                                                   const double * right_matrix,
                                                   const unsigned int * right_scaler,
                                                   unsigned int attrib);

PLL_EXPORT void pll_core_update_partial_ii_sse(unsigned int states,
                                               unsigned int sites,
                                               unsigned int rate_cats,
                                               double * parent_clv,
                                               unsigned int * parent_scaler,
                                               const double * left_clv,
                                               const double * right_clv,
                                               const double * left_matrix,
                                               const double * right_matrix,
                                               const unsigned int * left_scaler,
                                               const unsigned int * right_scaler,
                                               unsigned int attrib);

PLL_EXPORT void pll_core_update_partial_ii_4x4_sse(unsigned int sites,
                                                   unsigned int rate_cats,
                                                   double * parent_clv,
                                                   unsigned int * parent_scaler,
                                                   const double * left_clv,
                                                   const double * right_clv,
                                                   const double * left_matrix,
                                                   const double * right_matrix,
                                                   const unsigned int * left_scaler,
                                                   const unsigned int * right_scaler,
                                                   unsigned int attrib);
#endif

/* functions in core_partials_avx.c */

#ifdef HAVE_AVX
PLL_EXPORT void pll_core_create_lookup_avx(unsigned int states,
                                           unsigned int rate_cats,
                                           double * lookup,
                                           const double * left_matrix,
                                           const double * right_matrix,
                                           const unsigned int * tipmap,
                                           unsigned int tipmap_size);

PLL_EXPORT void pll_core_create_lookup_4x4_avx(unsigned int rate_cats,
                                               double * lookup,
                                               const double * left_matrix,
                                               const double * right_matrix);

PLL_EXPORT void pll_core_create_lookup_20x20_avx(unsigned int rate_cats,
                                               double * ttlookup,
                                               const double * left_matrix,
                                               const double * right_matrix,
                                               const unsigned int * tipmap,
                                               unsigned int tipmap_size);

PLL_EXPORT void pll_core_update_partial_tt_avx(unsigned int states,
                                               unsigned int sites,
                                               unsigned int rate_cats,
                                               double * parent_clv,
                                               unsigned int * parent_scaler,
                                               const unsigned char * left_tipchars,
                                               const unsigned char * right_tipchars,
                                               const double * lookup,
                                               unsigned int tipstates_count,
                                               unsigned int attrib);

PLL_EXPORT void pll_core_update_partial_tt_4x4_avx(unsigned int sites,
                                                   unsigned int rate_cats,
                                                   double * parent_clv,
                                                   unsigned int * parent_scaler,
                                                   const unsigned char * left_tipchars,
                                                   const unsigned char * right_tipchars,
                                                   const double * lookup,
                                                   unsigned int attrib);

PLL_EXPORT void pll_core_update_partial_ti_avx(unsigned int states,
                                               unsigned int sites,
                                               unsigned int rate_cats,
                                               double * parent_clv,
                                               unsigned int * parent_scaler,
                                               const unsigned char * left_tipchars,
                                               const double * right_clv,
                                               const double * left_matrix,
                                               const double * right_matrix,
                                               const unsigned int * right_scaler,
                                               const unsigned int * tipmap,
                                               unsigned int tipmap_size,
                                               unsigned int attrib);

PLL_EXPORT void pll_core_update_partial_ti_4x4_avx(unsigned int sites,
                                                   unsigned int rate_cats,
                                                   double * parent_clv,
                                                   unsigned int * parent_scaler,
                                                   const unsigned char * left_tipchar,
                                                   const double * right_clv,
                                                   const double * left_matrix,
                                                   const double * right_matrix,
                                                   const unsigned int * right_scaler,
                                                   unsigned int attrib);

PLL_EXPORT void pll_core_update_partial_ti_20x20_avx(unsigned int sites,
                                                     unsigned int rate_cats,
                                                     double * parent_clv,
                                                     unsigned int * parent_scaler,
                                                     const unsigned char * left_tipchar,
                                                     const double * right_clv,
                                                     const double * left_matrix,
                                                     const double * right_matrix,
                                                     const unsigned int * right_scaler,
                                                     const unsigned int * tipmap,
                                                     unsigned int tipmap_size);

PLL_EXPORT void pll_core_update_partial_ii_avx(unsigned int states,
                                               unsigned int sites,
                                               unsigned int rate_cats,
                                               double * parent_clv,
                                               unsigned int * parent_scaler,
                                               const double * left_clv,
                                               const double * right_clv,
                                               const double * left_matrix,
                                               const double * right_matrix,
                                               const unsigned int * left_scaler,
                                               const unsigned int * right_scaler,
                                               unsigned int attrib);

PLL_EXPORT void pll_core_update_partial_ii_4x4_avx(unsigned int sites,
                                                   unsigned int rate_cats,
                                                   double * parent_clv,
                                                   unsigned int * parent_scaler,
                                                   const double * left_clv,
                                                   const double * right_clv,
                                                   const double * left_matrix,
                                                   const double * right_matrix,
                                                   const unsigned int * left_scaler,
                                                   const unsigned int * right_scaler,
                                                   unsigned int attrib);
#endif

/* functions in core_partials_avx2.c */

#ifdef HAVE_AVX2
PLL_EXPORT void pll_core_update_partial_ti_avx2(unsigned int states,
                                                unsigned int sites,
                                                unsigned int rate_cats,
                                                double * parent_clv,
                                                unsigned int * parent_scaler,
                                                const unsigned char * left_tipchars,
                                                const double * right_clv,
                                                const double * left_matrix,
                                                const double * right_matrix,
                                                const unsigned int * right_scaler,
                                                const unsigned int * tipmap,
                                                unsigned int tipmap_size,
                                                unsigned int attrib);

PLL_EXPORT
void pll_core_update_partial_ti_20x20_avx2(unsigned int sites,
                                           unsigned int rate_cats,
                                           double * parent_clv,
                                           unsigned int * parent_scaler,
                                           const unsigned char * left_tipchar,
                                           const double * right_clv,
                                           const double * left_matrix,
                                           const double * right_matrix,
                                           const unsigned int * right_scaler,
                                           const unsigned int * tipmap,
                                           unsigned int tipmap_size);

PLL_EXPORT void pll_core_update_partial_ii_avx2(unsigned int states,
                                                unsigned int sites,
                                                unsigned int rate_cats,
                                                double * parent_clv,
                                                unsigned int * parent_scaler,
                                                const double * left_clv,
                                                const double * right_clv,
                                                const double * left_matrix,
                                                const double * right_matrix,
                                                const unsigned int * left_scaler,
                                                const unsigned int * right_scaler,
                                                unsigned int attrib);
#endif


/* functions in core_derivatives_sse.c */

#ifdef HAVE_SSE3
PLL_EXPORT int pll_core_update_sumtable_ii_sse(unsigned int states,
                                               unsigned int sites,
                                               unsigned int rate_cats,
                                               const double * parent_clv,
                                               const double * child_clv,
                                               const unsigned int * parent_scaler,
                                               const unsigned int * child_scaler,
                                               double * const * eigenvecs,
                                               double * const * inv_eigenvecs,
                                               double * const * freqs,
                                               double * sumtable,
                                               unsigned int attrib);

PLL_EXPORT int pll_core_update_sumtable_ti_sse(unsigned int states,
                                               unsigned int sites,
                                               unsigned int rate_cats,
                                               const double * parent_clv,
                                               const unsigned char * left_tipchars,
                                               const unsigned int * parent_scaler,
                                               double * const * eigenvecs,
                                               double * const * inv_eigenvecs,
                                               double * const * freqs,
                                               const unsigned int * tipmap,
                                               double * sumtable,
                                               unsigned int attrib);

#endif

/* functions in core_derivatives_avx.c */

#ifdef HAVE_AVX

PLL_EXPORT int pll_core_update_sumtable_ii_avx(unsigned int states,
                                               unsigned int sites,
                                               unsigned int rate_cats,
                                               const double * clvp,
                                               const double * clvc,
                                               const unsigned int * parent_scaler,
                                               const unsigned int * child_scaler,
                                               double * const * eigenvecs,
                                               double * const * inv_eigenvecs,
                                               double * const * freqs,
                                               double * sumtable,
                                               unsigned int attrib);

PLL_EXPORT int pll_core_update_sumtable_ti_avx(unsigned int states,
                                               unsigned int sites,
                                               unsigned int rate_cats,
                                               const double * parent_clv,
                                               const unsigned char * left_tipchars,
                                               const unsigned int * parent_scaler,
                                               double * const * eigenvecs,
                                               double * const * inv_eigenvecs,
                                               double * const * freqs,
                                               const unsigned int * tipmap,
                                               unsigned int tipmap_size,
                                               double * sumtable,
                                               unsigned int attrib);

PLL_EXPORT int pll_core_likelihood_derivatives_avx(unsigned int states,
                                                   unsigned int states_padded,
                                                   unsigned int rate_cats,
                                                   unsigned int ef_sites,
                                                   const unsigned int * pattern_weights,
                                                   const double * rate_weights,
                                                   const int * invariant,
                                                   const double * prop_invar,
                                                   double * const * freqs,
                                                   const double * sumtable,
                                                   const double * diagptable,
                                                   double * d_f,
                                                   double * dd_f);
#endif

/* functions in core_derivatives_avx2.c */

#ifdef HAVE_AVX2

PLL_EXPORT int pll_core_update_sumtable_ii_avx2(unsigned int states,
                                                unsigned int sites,
                                                unsigned int rate_cats,
                                                const double * clvp,
                                                const double * clvc,
                                                const unsigned int * parent_scaler,
                                                const unsigned int * child_scaler,
                                                double * const * eigenvecs,
                                                double * const * inv_eigenvecs,
                                                double * const * freqs,
                                                double * sumtable,
                                                unsigned int attrib);

PLL_EXPORT int pll_core_update_sumtable_ti_avx2(unsigned int states,
                                                unsigned int sites,
                                                unsigned int rate_cats,
                                                const double * parent_clv,
                                                const unsigned char * left_tipchars,
                                                const unsigned int * parent_scaler,
                                                double * const * eigenvecs,
                                                double * const * inv_eigenvecs,
                                                double * const * freqs,
                                                const unsigned int * tipmap,
                                                unsigned int tipmap_size,
                                                double * sumtable,
                                                unsigned int attrib);

PLL_EXPORT
int pll_core_likelihood_derivatives_avx2(unsigned int states,
                                         unsigned int states_padded,
                                         unsigned int rate_cats,
                                         unsigned int ef_sites,
                                         const unsigned int * pattern_weights,
                                         const double * rate_weights,
                                         const int * invariant,
                                         const double * prop_invar,
                                         double * const * freqs,
                                         const double * sumtable,
                                         const double * diagptable,
                                         double * d_f,
                                         double * dd_f);
#endif

/* functions in core_likelihood_sse.c */

#ifdef HAVE_SSE3
PLL_EXPORT
double pll_core_edge_loglikelihood_ii_sse(unsigned int states,
                                          unsigned int sites,
                                          unsigned int rate_cats,
                                          const double * parent_clv,
                                          const unsigned int * parent_scaler,
                                          const double * child_clv,
                                          const unsigned int * child_scaler,
                                          const double * pmatrix,
                                          double * const * frequencies,
                                          const double * rate_weights,
                                          const unsigned int * pattern_weights,
                                          const double * invar_proportion,
                                          const int * invar_indices,
                                          const unsigned int * freqs_indices,
                                          double * persite_lnl);

PLL_EXPORT
double pll_core_edge_loglikelihood_ii_4x4_sse(unsigned int sites,
                                              unsigned int rate_cats,
                                              const double * parent_clv,
                                              const unsigned int * parent_scaler,
                                              const double * child_clv,
                                              const unsigned int * child_scaler,
                                              const double * pmatrix,
                                              double * const * frequencies,
                                              const double * rate_weights,
                                              const unsigned int * pattern_weights,
                                              const double * invar_proportion,
                                              const int * invar_indices,
                                              const unsigned int * freqs_indices,
                                              double * persite_lnl,
                                              unsigned int attrib);

PLL_EXPORT
double pll_core_edge_loglikelihood_ti_sse(unsigned int states,
                                          unsigned int sites,
                                          unsigned int rate_cats,
                                          const double * parent_clv,
                                          const unsigned int * parent_scaler,
                                          const unsigned char * tipchars,
                                          const unsigned int * tipmap,
                                          const double * pmatrix,
                                          double * const * frequencies,
                                          const double * rate_weights,
                                          const unsigned int * pattern_weights,
                                          const double * invar_proportion,
                                          const int * invar_indices,
                                          const unsigned int * freqs_indices,
                                          double * persite_lnl);

PLL_EXPORT
double pll_core_edge_loglikelihood_ti_4x4_sse(unsigned int sites,
                                              unsigned int rate_cats,
                                              const double * parent_clv,
                                              const unsigned int * parent_scaler,
                                              const unsigned char * tipchars,
                                              const double * pmatrix,
                                              double * const * frequencies,
                                              const double * rate_weights,
                                              const unsigned int * pattern_weights,
                                              const double * invar_proportion,
                                              const int * invar_indices,
                                              const unsigned int * freqs_indices,
                                              double * persite_lnl,
                                              unsigned int attrib);

PLL_EXPORT double pll_core_root_loglikelihood_4x4_sse(unsigned int sites,
                                                      unsigned int rate_cats,
                                                      const double * clv,
                                                      const unsigned int * scaler,
                                                      double * const * frequencies,
                                                      const double * rate_weights,
                                                      const unsigned int * pattern_weights,
                                                      const double * invar_proportion,
                                                      const int * invar_indices,
                                                      const unsigned int * freqs_indices,
                                                      double * persite_lnl);

PLL_EXPORT double pll_core_root_loglikelihood_sse(unsigned int states,
                                                  unsigned int sites,
                                                  unsigned int rate_cats,
                                                  const double * clv,
                                                  const unsigned int * scaler,
                                                  double * const * frequencies,
                                                  const double * rate_weights,
                                                  const unsigned int * pattern_weights,
                                                  const double * invar_proportion,
                                                  const int * invar_indices,
                                                  const unsigned int * freqs_indices,
                                                  double * persite_lnl);
#endif

/* functions in core_likelihood_avx.c */

#ifdef HAVE_AVX
PLL_EXPORT double pll_core_edge_loglikelihood_ii_avx(unsigned int states,
                                                     unsigned int sites,
                                                     unsigned int rate_cats,
                                                     const double * parent_clv,
                                                     const unsigned int * parent_scaler,
                                                     const double * child_clv,
                                                     const unsigned int * child_scaler,
                                                     const double * pmatrix,
                                                     double * const * frequencies,
                                                     const double * rate_weights,
                                                     const unsigned int * pattern_weights,
                                                     const double * invar_proportion,
                                                     const int * invar_indices,
                                                     const unsigned int * freqs_indices,
                                                     double * persite_lnl);

PLL_EXPORT double pll_core_edge_loglikelihood_ii_4x4_avx(unsigned int sites,
                                                         unsigned int rate_cats,
                                                         const double * parent_clv,
                                                         const unsigned int * parent_scaler,
                                                         const double * child_clv,
                                                         const unsigned int * child_scaler,
                                                         const double * pmatrix,
                                                         double * const * frequencies,
                                                         const double * rate_weights,
                                                         const unsigned int * pattern_weights,
                                                         const double * invar_proportion,
                                                         const int * invar_indices,
                                                         const unsigned int * freqs_indices,
                                                         double * persite_lnl,
                                                         unsigned int attrib);

PLL_EXPORT double pll_core_edge_loglikelihood_ti_4x4_avx(unsigned int sites,
                                                         unsigned int rate_cats,
                                                         const double * parent_clv,
                                                         const unsigned int * parent_scaler,
                                                         const unsigned char * tipchars,
                                                         const double * pmatrix,
                                                         double * const * frequencies,
                                                         const double * rate_weights,
                                                         const unsigned int * pattern_weights,
                                                         const double * invar_proportion,
                                                         const int * invar_indices,
                                                         const unsigned int * freqs_indices,
                                                         double * persite_lnl,
                                                         unsigned int attrib);

PLL_EXPORT double pll_core_edge_loglikelihood_ti_20x20_avx(unsigned int sites,
                                                           unsigned int rate_cats,
                                                           const double * parent_clv,
                                                           const unsigned int * parent_scaler,
                                                           const unsigned char * tipchars,
                                                           const unsigned int * tipmap,
                                                           unsigned int tipmap_size,
                                                           const double * pmatrix,
                                                           double * const * frequencies,
                                                           const double * rate_weights,
                                                           const unsigned int * pattern_weights,
                                                           const double * invar_proportion,
                                                           const int * invar_indices,
                                                           const unsigned int * freqs_indices,
                                                           double * persite_lnl);

PLL_EXPORT double pll_core_edge_loglikelihood_ti_avx(unsigned int states,
                                                     unsigned int sites,
                                                     unsigned int rate_cats,
                                                     const double * parent_clv,
                                                     const unsigned int * parent_scaler,
                                                     const unsigned char * tipchars,
                                                     const unsigned int * tipmap,
                                                     const double * pmatrix,
                                                     double * const * frequencies,
                                                     const double * rate_weights,
                                                     const unsigned int * pattern_weights,
                                                     const double * invar_proportion,
                                                     const int * invar_indices,
                                                     const unsigned int * freqs_indices,
                                                     double * persite_lnl);

PLL_EXPORT double pll_core_root_loglikelihood_4x4_avx(unsigned int sites,
                                                      unsigned int rate_cats,
                                                      const double * clv,
                                                      const unsigned int * scaler,
                                                      double * const * frequencies,
                                                      const double * rate_weights,
                                                      const unsigned int * pattern_weights,
                                                      const double * invar_proportion,
                                                      const int * invar_indices,
                                                      const unsigned int * freqs_indices,
                                                      double * persite_lnl);

PLL_EXPORT double pll_core_root_loglikelihood_avx(unsigned int states,
                                                  unsigned int sites,
                                                  unsigned int rate_cats,
                                                  const double * clv,
                                                  const unsigned int * scaler,
                                                  double * const * frequencies,
                                                  const double * rate_weights,
                                                  const unsigned int * pattern_weights,
                                                  const double * invar_proportion,
                                                  const int * invar_indices,
                                                  const unsigned int * freqs_indices,
                                                  double * persite_lnl);
#endif


/* functions in core_likelihood_avx2.c */

#ifdef HAVE_AVX2
PLL_EXPORT
double pll_core_root_loglikelihood_avx2(unsigned int states,
                                        unsigned int sites,
                                        unsigned int rate_cats,
                                        const double * clv,
                                        const unsigned int * scaler,
                                        double * const * frequencies,
                                        const double * rate_weights,
                                        const unsigned int * pattern_weights,
                                        const double * invar_proportion,
                                        const int * invar_indices,
                                        const unsigned int * freqs_indices,
                                        double * persite_lnl);

PLL_EXPORT
double pll_core_edge_loglikelihood_ti_20x20_avx2(unsigned int sites,
                                                 unsigned int rate_cats,
                                                 const double * parent_clv,
                                                 const unsigned int * parent_scaler,
                                                 const unsigned char * tipchars,
                                                 const unsigned int * tipmap,
                                                 unsigned int tipmap_size,
                                                 const double * pmatrix,
                                                 double * const * frequencies,
                                                 const double * rate_weights,
                                                 const unsigned int * pattern_weights,
                                                 const double * invar_proportion,
                                                 const int * invar_indices,
                                                 const unsigned int * freqs_indices,
                                                 double * persite_lnl);


PLL_EXPORT
double pll_core_edge_loglikelihood_ii_avx2(unsigned int states,
                                           unsigned int sites,
                                           unsigned int rate_cats,
                                           const double * parent_clv,
                                           const unsigned int * parent_scaler,
                                           const double * child_clv,
                                           const unsigned int * child_scaler,
                                           const double * pmatrix,
                                           double * const * frequencies,
                                           const double * rate_weights,
                                           const unsigned int * pattern_weights,
                                           const double * invar_proportion,
                                           const int * invar_indices,
                                           const unsigned int * freqs_indices,
                                           double * persite_lnl);
#endif

/* functions in core_pmatrix.c */

PLL_EXPORT int pll_core_update_pmatrix(double ** pmatrix,
                                       unsigned int states,
                                       unsigned int rate_cats,
                                       const double * rates,
                                       const double * branch_lengths,
                                       const unsigned int * matrix_indices,
                                       const unsigned int * params_indices,
                                       const double * prop_invar,
                                       double * const * eigenvals,
                                       double * const * eigenvecs,
                                       double * const * inv_eigenvecs,
                                       unsigned int count,
                                       unsigned int attrib);

/* functions in core_pmatrix_avx2.c */

#ifdef HAVE_AVX2
PLL_EXPORT int pll_core_update_pmatrix_20x20_avx2(double ** pmatrix,
                                                  unsigned int rate_cats,
                                                  const double * rates,
                                                  const double * branch_lengths,
                                                  const unsigned int * matrix_indices,
                                                  const unsigned int * params_indices,
                                                  const double * prop_invar,
                                                  double * const * eigenvals,
                                                  double * const * eigenvecs,
                                                  double * const * inv_eigenvecs,
                                                  unsigned int count);
#endif

/* functions in core_pmatrix_avx.c */

#ifdef HAVE_AVX
PLL_EXPORT int pll_core_update_pmatrix_4x4_avx(double ** pmatrix,
                                               unsigned int rate_cats,
                                               const double * rates,
                                               const double * branch_lengths,
                                               const unsigned int * matrix_indices,
                                               const unsigned int * params_indices,
                                               const double * prop_invar,
                                               double * const * eigenvals,
                                               double * const * eigenvecs,
                                               double * const * inv_eigenvecs,
                                               unsigned int count);

PLL_EXPORT int pll_core_update_pmatrix_20x20_avx(double ** pmatrix,
                                                 unsigned int rate_cats,
                                                 const double * rates,
                                                 const double * branch_lengths,
                                                 const unsigned int * matrix_indices,
                                                 const unsigned int * params_indices,
                                                 const double * prop_invar,
                                                 double * const * eigenvals,
                                                 double * const * eigenvecs,
                                                 double * const * inv_eigenvecs,
                                                 unsigned int count);
#endif

/* functions in core_pmatrix_sse.c */

#ifdef HAVE_SSE3
PLL_EXPORT int pll_core_update_pmatrix_4x4_sse(double ** pmatrix,
                                               unsigned int rate_cats,
                                               const double * rates,
                                               const double * branch_lengths,
                                               const unsigned int * matrix_indices,
                                               const unsigned int * params_indices,
                                               const double * prop_invar,
                                               double * const * eigenvals,
                                               double * const * eigenvecs,
                                               double * const * inv_eigenvecs,
                                               unsigned int count);
#endif

/* functions in compress.c */

PLL_EXPORT unsigned int * pll_compress_site_patterns(char ** sequence,
                                                     const unsigned int * map,
                                                     int count,
                                                     int * length);

/* functions in utree_moves.c */

PLL_EXPORT int pll_utree_spr(pll_unode_t * p,
                             pll_unode_t * r,
                             pll_utree_rb_t * rb,
                             double * branch_lengths,
                             unsigned int * matrix_indices);

PLL_EXPORT int pll_utree_spr_safe(pll_unode_t * p,
                                  pll_unode_t * r,
                                  pll_utree_rb_t * rb,
                                  double * branch_lengths,
                                  unsigned int * matrix_indices);

PLL_EXPORT int pll_utree_nni(pll_unode_t * p,
                             int type,
                             pll_utree_rb_t * rb);

PLL_EXPORT int pll_utree_rollback(pll_utree_rb_t * rollback,
                                  double * branch_lengths,
                                  unsigned int * matrix_indices);

/* functions in parsimony.c */

PLL_EXPORT int pll_set_parsimony_sequence(pll_parsimony_t * pars,
                                          unsigned int tip_index,
                                          const unsigned int * map,
                                          const char * sequence);

PLL_EXPORT pll_parsimony_t * pll_parsimony_create(unsigned int tips,
                                                  unsigned int states,
                                                  unsigned int sites,
                                                  const double * score_matrix,
                                                  unsigned int score_buffers,
                                                  unsigned int ancestral_buffers);

PLL_EXPORT double pll_parsimony_build(pll_parsimony_t * pars,
                                      const pll_pars_buildop_t * operations,
                                      unsigned int count);

PLL_EXPORT void pll_parsimony_reconstruct(pll_parsimony_t * pars,
                                          const unsigned int * map,
                                          const pll_pars_recop_t * operations,
                                          unsigned int count);

PLL_EXPORT double pll_parsimony_score(pll_parsimony_t * pars,
                                      unsigned int score_buffer_index);

PLL_EXPORT void pll_parsimony_destroy(pll_parsimony_t * pars);

/* functions in utree_svg.c */

PLL_EXPORT pll_svg_attrib_t * pll_svg_attrib_create(void);

PLL_EXPORT void pll_svg_attrib_destroy(pll_svg_attrib_t * attrib);

PLL_EXPORT int pll_utree_export_svg(pll_utree_t * tree,
                                    pll_unode_t * root,
                                    const pll_svg_attrib_t * attribs,
                                    const char * filename);

/* functions in fast_parsimony.c */

PLL_EXPORT pll_parsimony_t * pll_fastparsimony_init(const pll_partition_t * partition);

PLL_EXPORT void pll_fastparsimony_update_vectors(pll_parsimony_t * parsimony,
                                                 const pll_pars_buildop_t * ops,
                                                 unsigned int count);

PLL_EXPORT unsigned int pll_fastparsimony_root_score(const pll_parsimony_t * parsimony,
                                                     unsigned int root_index);

PLL_EXPORT unsigned int pll_fastparsimony_edge_score(const pll_parsimony_t * parsimony,
                                                     unsigned int node1_score_index,
                                                     unsigned int node2_score_index);

PLL_EXPORT void pll_fastparsimony_update_vector_4x4(pll_parsimony_t * parsimony,
                                                    const pll_pars_buildop_t * op);

PLL_EXPORT unsigned int pll_fastparsimony_edge_score_4x4(const pll_parsimony_t * parsimony,
                                                         unsigned int node1_score_index,
                                                         unsigned int node2_score_index);

PLL_EXPORT void pll_fastparsimony_update_vector(pll_parsimony_t * parsimony,
                                                const pll_pars_buildop_t * op);

/* functions in fast_parsimony_sse.c */

PLL_EXPORT void pll_fastparsimony_update_vector_4x4_sse(pll_parsimony_t * parsimony,
                                                         const pll_pars_buildop_t * op);

PLL_EXPORT unsigned int pll_fastparsimony_edge_score_4x4_sse(const pll_parsimony_t * parsimony,
                                                             unsigned int node1_score_index,
                                                             unsigned int node2_score_index);

PLL_EXPORT unsigned int pll_fastparsimony_edge_score_sse(const pll_parsimony_t * parsimony,
                                                         unsigned int node1_score_index,
                                                         unsigned int node2_score_index);

PLL_EXPORT void pll_fastparsimony_update_vector_sse(pll_parsimony_t * parsimony,
                                                    const pll_pars_buildop_t * op);

/* functions in fast_parsimony_avx.c */

PLL_EXPORT void pll_fastparsimony_update_vector_4x4_avx(pll_parsimony_t * parsimony,
                                                         const pll_pars_buildop_t * op);

PLL_EXPORT unsigned int pll_fastparsimony_edge_score_4x4_avx(const pll_parsimony_t * parsimony,
                                                             unsigned int node1_score_index,
                                                             unsigned int node2_score_index);

PLL_EXPORT void pll_fastparsimony_update_vector_avx(pll_parsimony_t * parsimony,
                                                    const pll_pars_buildop_t * op);


PLL_EXPORT unsigned int pll_fastparsimony_edge_score_avx(const pll_parsimony_t * parsimony,
                                                         unsigned int node1_score_index,
                                                         unsigned int node2_score_index);

/* functions in fast_parsimony_avx2.c */

PLL_EXPORT void pll_fastparsimony_update_vector_4x4_avx2(pll_parsimony_t * parsimony,
                                                          const pll_pars_buildop_t * op);

PLL_EXPORT unsigned int pll_fastparsimony_edge_score_4x4_avx2(const pll_parsimony_t * parsimony,
                                                              unsigned int node1_score_index,
                                                              unsigned int node2_score_index);

PLL_EXPORT void pll_fastparsimony_update_vector_avx2(pll_parsimony_t * parsimony,
                                                     const pll_pars_buildop_t * op);

PLL_EXPORT unsigned int pll_fastparsimony_edge_score_avx2(const pll_parsimony_t * parsimony,
                                                          unsigned int node1_score_index,
                                                          unsigned int node2_score_index);

/* functions in stepwise.c */

PLL_EXPORT pll_utree_t * pll_fastparsimony_stepwise(pll_parsimony_t ** list,
                                                    char * const * labels,
                                                    unsigned int * score,
                                                    unsigned int count,
                                                    unsigned int seed);

/* functions in random.c */

PLL_EXPORT extern int pll_random_r(struct pll_random_data * __buf,
                                   int32_t * __result);

PLL_EXPORT extern int pll_srandom_r(unsigned int __seed,
                                    struct pll_random_data * __buf);

PLL_EXPORT extern int pll_initstate_r(unsigned int __seed,
                                      char * __statebuf,
                                      size_t __statelen,
                                      struct pll_random_data * __buf);

PLL_EXPORT extern int pll_setstate_r(char * __statebuf,
                                     struct pll_random_data * __buf);

/* functions in hardware.c */

PLL_EXPORT int pll_hardware_probe(void);

PLL_EXPORT void pll_hardware_dump();

PLL_EXPORT void pll_hardware_ignore();

#ifdef __cplusplus
} /* extern "C" */
#endif
#endif
