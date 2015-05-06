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

/* constants */

#define PLL_FAILURE  0
#define PLL_SUCCESS  1

#define PLL_ALIGNMENT_CPU   8
#define PLL_ALIGNMENT_SSE  16
#define PLL_ALIGNMENT_AVX  32

/* attribute flags */

#define PLL_ATTRIB_ARCH_SSE       1 << 0
#define PLL_ATTRIB_ARCH_AVX       1 << 1
#define PLL_ATTRIB_ARCH_AVX2      1 << 2
#define PLL_ATTRIB_ARCH_AVX512    1 << 3

/* error codes */

#define PLL_ERROR_UNKNOWN_PARTITION   1

/* structures and data types */

typedef struct
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
  double prop_invar;

  /* multiple to which memory is aligned */
  size_t alignment;

  double ** clv;
  double ** pmatrix;
  double * rates;
  double ** subst_params;
  double * scale_buffer;
  double ** frequencies;
  int * invariant;

  /* eigen decomposition */
  int * eigen_decomp_valid;
  double ** eigenvecs;
  double ** inv_eigenvecs;
  double ** eigenvals;
} pll_partition_t;

typedef struct
{
  int parent_clv_index;
  int child1_clv_index;
  int child1_matrix_index;
  int child2_clv_index;
  int child2_matrix_index;
} pll_operation_t;

typedef struct pll_dlist
{
  struct pll_dlist * next;
  struct pll_dlist * prev;
  void * data;
} pll_dlist_t;


/* common data */

extern int pll_errno;
extern unsigned int pll_map_bin[256];
extern unsigned int pll_map_nt[256];
extern unsigned int pll_map_aa[256];

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

PLL_EXPORT int pll_destroy_partition(pll_partition_t * partition);

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
                                     double * params, 
                                     int count);

PLL_EXPORT void pll_set_frequencies(pll_partition_t * partition, 
                                    int params_index,
                                    double * frequencies);

PLL_EXPORT void pll_set_category_rates(pll_partition_t * partition, 
                                       double * rates);

PLL_EXPORT void pll_update_prob_matrices(pll_partition_t * partition, 
                                         int params_index, 
                                         int * matrix_indices, 
                                         double * branch_lengths, 
                                         int count);

PLL_EXPORT int pll_update_invariant_sites(pll_partition_t * partition);

PLL_EXPORT int pll_update_invariant_sites_proportion(pll_partition_t * partition, 
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

PLL_EXPORT void pll_show_pmatrix(pll_partition_t * partition, int index);
PLL_EXPORT void pll_show_clv(pll_partition_t * partition, int index);

#ifdef __cplusplus
} /* extern "C" */
#endif
