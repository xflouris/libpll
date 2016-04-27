/*
 * utils.h
 *
 *  Created on: Sep 28, 2015
 *      Author: diego
 */

#ifndef UTILS_H_
#define UTILS_H_

#include "pll_optimize.h"

#include <assert.h>
#include <stdarg.h>
#include <search.h>
#include <time.h>

#define STATES    4
#define RATE_CATS 4

#define SUBST_PARAMS (STATES*(STATES-1)/2)

typedef struct
{
  pll_partition_t * partition;
  unsigned int parent_clv_index;
  unsigned int parent_scaler_index;
  unsigned int child_clv_index;
  unsigned int child_scaler_index;
  unsigned int edge_pmatrix_index;
  unsigned int freqs_indices[RATE_CATS];
  unsigned int params_index;

  /* traverse */
  unsigned int * matrix_indices;
  double * branch_lengths;
  pll_operation_t * operations;

  /* subst rate specific */
  int * symmetries;
  int n_subst_params;

  /* frequency specific */
  unsigned int highest_freq_state;
} my_params_t;

void fatal (const char * format, ...) __attribute__ ((noreturn));

pll_partition_t * partition_fasta_create (const char *file,
                                          unsigned int states,
                                          unsigned int n_rate_matrices,
                                          unsigned int n_rate_cats,
                                          int attributes,
                                          int rooted,
                                          unsigned int tip_count,
                                          const char **tipnames);

void set_missing_branch_length (pll_utree_t * tree, double length);

int * build_model_symmetries (const char * modelmatrix);

/* a callback function for performing a full traversal */
int cb_full_traversal(pll_utree_t * node);

/* a callback function for rates optimization */
double target_rates_opt(void * params, double * x);
/* a callback function for frequencies optimization */
double target_freqs_opt (void * p, double * x);

#endif /* UTILS_H_ */
