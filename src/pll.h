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

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <x86intrin.h>

#define LIB_NAME "libpll"
#define LIB_VERSION "v0.0.1"

/* platform specific */

#if (!defined(__APPLE__) && !defined(__WIN32__) && !defined(__WIN64__))
#include <sys/sysinfo.h>
#endif

#if (defined(__WIN32__) || defined(__WIN64__))
#define PLL_EXPORT __declspec(dllexport)
#else
#define PLL_EXPORT
#endif

/* constants */

#define PLL_FAILURE  0
#define PLL_SUCCESS  1

#define PLL_ALIGNMENT_CPU   8
#define PLL_ALIGNMENT_SSE  16
#define PLL_ALIGNMENT_AVX  32

#define PLL_LINEALLOC 2048

#define PLL_SCALE_FACTOR 115792089237316195423570985008687907853269984665640564039457584007913129639936.0  /*  2**256 (exactly)  */
#define PLL_SCALE_THRESHOLD (1.0/PLL_SCALE_FACTOR)
#define PLL_SCALE_BUFFER_NONE -1

/* attribute flags */

#define PLL_ATTRIB_ARCH_CPU            0
#define PLL_ATTRIB_ARCH_SSE       1 << 0
#define PLL_ATTRIB_ARCH_AVX       1 << 1
#define PLL_ATTRIB_ARCH_AVX2      1 << 2
#define PLL_ATTRIB_ARCH_AVX512    1 << 3
#define PLL_ATTRIB_ARCH_MASK         0xF
#define PLL_ATTRIB_MIXT_LINKED    1 << 4  /** Q matrices linked to rate categories */
#define PLL_ATTRIB_MIXT_UNLINKED  1 << 5  /** Q matrices unlinked */
#define PLL_ATTRIB_MIXT_MASK        0xF0

/* error codes */

#define PLL_ERROR_FILE_OPEN                1
#define PLL_ERROR_FILE_SEEK                2
#define PLL_ERROR_FILE_EOF                 3
#define PLL_ERROR_FASTA_ILLEGALCHAR        4
#define PLL_ERROR_FASTA_UNPRINTABLECHAR    5
#define PLL_ERROR_FASTA_INVALIDHEADER      6
#define PLL_ERROR_MEM_ALLOC                7
#define PLL_ERROR_NEWICK_SYNTAX            8
#define PLL_ERROR_TIP_DATA_ILLEGAL_STATE   9
#define PLL_ERROR_MULTIPLE_ARCH           10

#define PLL_ERROR_ALPHA                  101
#define PLL_ERROR_PINV                   102
#define PLL_ERROR_MIXTURE                103

/* utree specific */

#define PLL_UTREE_SHOW_LABEL             1 << 0
#define PLL_UTREE_SHOW_BRANCH_LENGTH     1 << 1
#define PLL_UTREE_SHOW_CLV_INDEX         1 << 2
#define PLL_UTREE_SHOW_SCALER_INDEX      1 << 3
#define PLL_UTREE_SHOW_PMATRIX_INDEX     1 << 4



/* structures and data types */

typedef struct pll_partition
{
  unsigned int tips;
  unsigned int clv_buffers;
  unsigned int states;
  unsigned int sites;
  unsigned int rate_matrices;
  unsigned int prob_matrices;
  unsigned int rate_cats;
  unsigned int scale_buffers;
  unsigned int attributes;

  /* vectorization options */
  size_t alignment;
  unsigned int states_padded;

  unsigned int mixture;
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

/* Simple unrooted tree structure for parsing newick */

typedef struct pll_utree
{
  char * label;
  double length;
  unsigned int clv_index;
  int scaler_index;
  unsigned int pmatrix_index;
  struct pll_utree * next;
  struct pll_utree * back;

  void * data;
} pll_utree_t;

typedef struct pll_rtree
{
  char * label;
  double length;
  unsigned int clv_index;
  int scaler_index;
  unsigned int pmatrix_index;
  struct pll_rtree * left;
  struct pll_rtree * right;
  struct pll_rtree * parent;

  void * data;
} pll_rtree_t;

/* common data */

PLL_EXPORT extern int pll_errno;
PLL_EXPORT extern char pll_errmsg[200];
PLL_EXPORT extern const unsigned int pll_map_bin[256];
PLL_EXPORT extern const unsigned int pll_map_nt[256];
PLL_EXPORT extern const unsigned int pll_map_aa[256];
PLL_EXPORT extern const unsigned int pll_map_fasta[256];

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
                                                  unsigned int mixture,
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

PLL_EXPORT void pll_set_tip_clv(pll_partition_t * partition,
                                unsigned int tip_index,
                                const double * clv);

PLL_EXPORT void pll_set_pattern_weights(pll_partition_t * partition,
                                        const unsigned int * pattern_weights);

/* functions in list.c */

PLL_EXPORT int pll_dlist_append(pll_dlist_t ** dlist, void * data);
PLL_EXPORT int pll_dlist_remove(pll_dlist_t ** dlist, void * data);
PLL_EXPORT int pll_dlist_prepend(pll_dlist_t ** dlist, void * data);

/* functions in models.c */

PLL_EXPORT void pll_set_subst_params(pll_partition_t * partition, 
                                     unsigned int params_index,
                                     unsigned int mixture_index,
                                     const double * params);

PLL_EXPORT void pll_set_frequencies(pll_partition_t * partition, 
                                    unsigned int params_index,
                                    unsigned int mixture_index,
                                    const double * frequencies);

PLL_EXPORT void pll_set_category_rates(pll_partition_t * partition,
                                       const double * rates);

PLL_EXPORT void pll_set_category_weights(pll_partition_t * partition,
                                         const double * rate_weights);

PLL_EXPORT void pll_update_prob_matrices(pll_partition_t * partition, 
                                         unsigned int params_index, 
                                         unsigned int * matrix_indices, 
                                         double * branch_lengths, 
                                         unsigned int count);

PLL_EXPORT int pll_update_invariant_sites(pll_partition_t * partition);

PLL_EXPORT int pll_update_invariant_sites_proportion(pll_partition_t * partition,
                                                     unsigned int params_index,
                                                     double prop_invar);

PLL_EXPORT void * pll_aligned_alloc(size_t size, size_t alignment);

PLL_EXPORT void pll_aligned_free(void * ptr);

/* functions in likelihood.c */

PLL_EXPORT void pll_update_partials(pll_partition_t * partition,
                                    const pll_operation_t * operations,
                                    unsigned int count);

PLL_EXPORT double pll_compute_root_loglikelihood(pll_partition_t * partition,
                                                 unsigned int clv_index,
                                                 int scaler_index,
                                                 unsigned int freqs_index);

PLL_EXPORT double pll_compute_edge_loglikelihood(pll_partition_t * partition,
                                                 unsigned int parent_clv_index,
                                                 int parent_scaler_index,
                                                 unsigned int child_clv_index,
                                                 int child_scaler_index,
                                                 unsigned int matrix_index,
                                                 unsigned int freqs_index);

/* functions in gamma.c */

PLL_EXPORT int pll_compute_gamma_cats(double alpha,
                                      unsigned int categories,
                                      double * output_rates);

/* functions in output.c */

PLL_EXPORT void pll_show_pmatrix(pll_partition_t * partition,
                                 unsigned int index,
                                 unsigned int float_precision);

PLL_EXPORT void pll_show_clv(pll_partition_t * partition,
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

PLL_EXPORT long pll_fasta_getfilesize(pll_fasta_t * fd);

PLL_EXPORT long pll_fasta_getfilepos(pll_fasta_t * fd);

PLL_EXPORT int pll_fasta_rewind(pll_fasta_t * fd);

/* functions in parse_rtree.y */

PLL_EXPORT pll_rtree_t * pll_rtree_parse_newick(const char * filename,
                                                unsigned int * tip_count);

PLL_EXPORT void pll_rtree_destroy(pll_rtree_t * root);

/* functions in parse_utree.y */

PLL_EXPORT pll_utree_t * pll_utree_parse_newick(const char * filename,
                                                unsigned int * tip_count);

PLL_EXPORT void pll_utree_destroy(pll_utree_t * root);

/* functions in utree.c */

PLL_EXPORT void pll_utree_show_ascii(pll_utree_t * tree, int options);

PLL_EXPORT char * pll_utree_export_newick(pll_utree_t * root);

PLL_EXPORT int pll_utree_traverse(pll_utree_t * root,
                                  int (*cbtrav)(pll_utree_t *),
                                  pll_utree_t ** outbuffer,
                                  unsigned int * trav_size);

PLL_EXPORT unsigned int pll_utree_query_tipnodes(pll_utree_t * root,
                                                 pll_utree_t ** node_list);

PLL_EXPORT unsigned int pll_utree_query_innernodes(pll_utree_t * root,
                                                   pll_utree_t ** node_list);

PLL_EXPORT void pll_utree_create_operations(pll_utree_t ** trav_buffer,
                                            unsigned int trav_buffer_size,
                                            double * branches,
                                            unsigned int * pmatrix_indices,
                                            pll_operation_t * ops,
                                            unsigned int * matrix_count,
                                            unsigned int * ops_count);

/* functions in rtree.c */

PLL_EXPORT void pll_rtree_show_ascii(pll_rtree_t * tree, int options);

PLL_EXPORT char * pll_rtree_export_newick(pll_rtree_t * root);

PLL_EXPORT int pll_rtree_traverse(pll_rtree_t * root,
                                  int (*cbtrav)(pll_rtree_t *),
                                  pll_rtree_t ** outbuffer,
                                  unsigned int * trav_size);

PLL_EXPORT unsigned int pll_rtree_query_tipnodes(pll_rtree_t * root,
                                                 pll_rtree_t ** node_list);

PLL_EXPORT unsigned int pll_rtree_query_innernodes(pll_rtree_t * root,
                                                   pll_rtree_t ** node_list);

PLL_EXPORT void pll_rtree_create_operations(pll_rtree_t ** trav_buffer,
                                            unsigned int trav_buffer_size,
                                            double * branches,
                                            unsigned int * pmatrix_indices,
                                            pll_operation_t * ops,
                                            unsigned int * matrix_count,
                                            unsigned int * ops_count);

/* functions in core_likelihood.c */

PLL_EXPORT void pll_core_update_partial(unsigned int states,
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

/* functions in likelihood_avx.c */

void pll_update_partials_avx(pll_partition_t * partition,
                             const pll_operation_t * op);
#ifdef __cplusplus
} /* extern "C" */
#endif
