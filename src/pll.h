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

#ifndef __APPLE__
#include <sys/sysinfo.h>
#endif

#ifdef _WIN32
#define PLL_EXPORT __declspec(dllexport)
#else
#define PLL_EXPORT
#endif

/* macros */

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* constants */

#define PLL_FAILURE  0
#define PLL_SUCCESS  1

#define PLL_ALIGNMENT_CPU   8
#define PLL_ALIGNMENT_SSE  16
#define PLL_ALIGNMENT_AVX  32

#define PLL_LINEALLOC 2048

/* attribute flags */

#define PLL_ATTRIB_ARCH_SSE       1 << 0
#define PLL_ATTRIB_ARCH_AVX       1 << 1
#define PLL_ATTRIB_ARCH_AVX2      1 << 2
#define PLL_ATTRIB_ARCH_AVX512    1 << 3

/* error codes */

#define PLL_ERROR_FILE_OPEN              1 
#define PLL_ERROR_FILE_SEEK              2
#define PLL_ERROR_FILE_EOF               3
#define PLL_ERROR_FASTA_ILLEGALCHAR      4
#define PLL_ERROR_FASTA_UNPRINTABLECHAR  5
#define PLL_ERROR_FASTA_INVALIDHEADER    6
#define PLL_ERROR_MEM_ALLOC              7
#define PLL_ERROR_NEWICK_SYNTAX          8


/* structures and data types */

typedef struct pll_partition
{
  int tips;
  int clv_buffers;
  int states;
  int sites;
  int rate_matrices;
  int prob_matrices;
  int rate_cats;
  int scale_buffers;
  int attributes;

  size_t alignment;

  double ** clv;
  double ** pmatrix;
  double * rates;
  double ** subst_params;
  double * scale_buffer;
  double ** frequencies;
  double * prop_invar;
  int * invariant;

  int * eigen_decomp_valid;
  double ** eigenvecs;
  double ** inv_eigenvecs;
  double ** eigenvals;
} pll_partition_t;


/* Structure for driving likelihood operations */

typedef struct pll_operation
{
  int parent_clv_index;
  int child1_clv_index;
  int child1_matrix_index;
  int child2_clv_index;
  int child2_matrix_index;
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
  unsigned int * chrstatus;
  long no;
  long filesize;
  long lineno;
  long stripped_count;
  long stripped[256];
} pll_fasta_t;

/* Simple unrooted tree structure for parsing newick */

typedef struct tree_noderec
{
  char * label;
  double length;
  struct tree_noderec * next;
  struct tree_noderec * back;

  void * data;
} pll_utree_t;

/* common data */

extern int pll_errno;
extern char pll_errmsg[200];
extern unsigned int pll_map_bin[256];
extern unsigned int pll_map_nt[256];
extern unsigned int pll_map_aa[256];
extern unsigned int pll_map_fasta[256];

extern double pll_aa_rates_dayhoff[190];
extern double pll_aa_rates_lg[190];
extern double pll_aa_rates_dcmut[190];
extern double pll_aa_rates_jtt[190];
extern double pll_aa_rates_mtrev[190];
extern double pll_aa_rates_wag[190];
extern double pll_aa_rates_rtrev[190];
extern double pll_aa_rates_cprev[190];
extern double pll_aa_rates_vt[190];
extern double pll_aa_rates_blosum62[190];
extern double pll_aa_rates_mtmam[190];
extern double pll_aa_rates_mtart[190];
extern double pll_aa_rates_mtzoa[190];
extern double pll_aa_rates_pmb[190];
extern double pll_aa_rates_hivb[190];
extern double pll_aa_rates_hivw[190];
extern double pll_aa_rates_jttdcmut[190];
extern double pll_aa_rates_flu[190];
extern double pll_aa_rates_stmtrev[190];

extern double pll_aa_freqs_dayhoff[20];
extern double pll_aa_freqs_lg[20];
extern double pll_aa_freqs_dcmut[20];
extern double pll_aa_freqs_jtt[20];
extern double pll_aa_freqs_mtrev[20];
extern double pll_aa_freqs_wag[20];
extern double pll_aa_freqs_rtrev[20];
extern double pll_aa_freqs_cprev[20];
extern double pll_aa_freqs_vt[20];
extern double pll_aa_freqs_blosum62[20];
extern double pll_aa_freqs_mtmam[20];
extern double pll_aa_freqs_mtart[20];
extern double pll_aa_freqs_mtzoa[20];
extern double pll_aa_freqs_pmb[20];
extern double pll_aa_freqs_hivb[20];
extern double pll_aa_freqs_hivw[20];
extern double pll_aa_freqs_jttdcmut[20];
extern double pll_aa_freqs_flu[20];
extern double pll_aa_freqs_stmtrev[20];

#ifdef __cplusplus
extern "C" {
#endif

/* functions in pll.c */

PLL_EXPORT pll_partition_t * pll_create_partition(int tips,
                                                  int clv_buffers,
                                                  int states,
                                                  int sites,
                                                  int rate_matrices,
                                                  int prob_matrices,
                                                  int rate_cats,
                                                  int scale_buffers,
                                                  int attributes);

PLL_EXPORT void pll_destroy_partition(pll_partition_t * partition);

PLL_EXPORT int pll_set_tip_states(pll_partition_t * partition, 
                                  int tip_index, 
                                  const unsigned int * map,
                                  const char * sequence);

PLL_EXPORT void pll_set_tip_clv(pll_partition_t * partition,
                                int tip_index,
                                const double * clv);


/* functions in dlist.c */

PLL_EXPORT int pll_dlist_append(pll_dlist_t ** dlist, void * data);
PLL_EXPORT int pll_dlist_remove(pll_dlist_t ** dlist, void * data);
PLL_EXPORT int pll_dlist_prepend(pll_dlist_t ** dlist, void * data);

/* functions in models.c */

PLL_EXPORT void pll_set_subst_params(pll_partition_t * partition, 
                                     int params_index, 
                                     const double * params);

PLL_EXPORT void pll_set_frequencies(pll_partition_t * partition, 
                                    int params_index,
                                    const double * frequencies);

PLL_EXPORT void pll_set_category_rates(pll_partition_t * partition,
                                       const double * rates);

PLL_EXPORT void pll_update_prob_matrices(pll_partition_t * partition, 
                                         int params_index, 
                                         int * matrix_indices, 
                                         double * branch_lengths, 
                                         int count);

PLL_EXPORT int pll_update_invariant_sites(pll_partition_t * partition);

PLL_EXPORT int pll_update_invariant_sites_proportion(pll_partition_t * partition, 
                                                     int params_index,
                                                     double prop_invar);

/* functions in likelihood.c */

PLL_EXPORT void pll_update_partials(pll_partition_t * partition, 
                                    pll_operation_t * operations, 
                                    int count);

PLL_EXPORT double pll_compute_root_loglikelihood(pll_partition_t * partition, 
                                                 int clv_index, 
                                                 int freqs_index);

PLL_EXPORT double pll_compute_edge_loglikelihood(pll_partition_t * partition, 
                                                 int parent_clv_index, 
                                                 int child_clv_index, 
                                                 int matrix_index,
                                                 int freqs_index);

/* functions in gamma.c */

PLL_EXPORT int pll_compute_gamma_cats(double alpha, 
                                      int categories, 
                                      double * output_rates);

/* functions in output.c */

PLL_EXPORT void pll_show_pmatrix(pll_partition_t * partition, 
                                 int index, 
                                 int float_precision);

PLL_EXPORT void pll_show_clv(pll_partition_t * partition, 
                             int index, 
                             int float_precision);

/* functions in fasta.c */

PLL_EXPORT pll_fasta_t * pll_fasta_open(const char * filename,
                                        unsigned int * map);

PLL_EXPORT int pll_fasta_getnext(pll_fasta_t * fd, char ** head,
                                 long * head_len,  char ** seq,
                                 long * seq_len, long * seqno);

PLL_EXPORT void pll_fasta_close(pll_fasta_t * fd);

PLL_EXPORT long pll_fasta_getfilesize(pll_fasta_t * fd);

PLL_EXPORT long pll_fasta_getfilepos(pll_fasta_t * fd);

/* functions in unrooted.y */

PLL_EXPORT pll_utree_t * pll_parse_newick_utree(const char * filename, 
                                                int * tip_count);

PLL_EXPORT void pll_destroy_utree(pll_utree_t * root);

/* functions in tree.c */

PLL_EXPORT void pll_show_ascii_utree(pll_utree_t * tree);

PLL_EXPORT char * pll_write_newick_utree(pll_utree_t * root);

PLL_EXPORT void pll_traverse_utree(pll_utree_t * tree, 
                                   int tips, 
                                   double ** branch_lengths, 
                                   int ** indices,
                                   pll_operation_t ** ops,
                                   int * edge_pmatrix_index,
                                   int * edge_node1_clv_index,
                                   int * edge_node2_clv_index);

PLL_EXPORT char ** pll_query_utree_tipnames(pll_utree_t * tree,
                                            int tips);

#ifdef __cplusplus
} /* extern "C" */
#endif
