/*
 * optimize.h
 *
 *  Created on: May 12, 2015
 *      Author: diego
 */
#include "pll.h"

#ifndef OPTIMIZE_H_
#define OPTIMIZE_H_

typedef enum {
  PARAM_SUBST_RATES,
  PARAM_PINV,
  PARAM_ALPHA,
  PARAM_FREQUENCIES
} pll_parameter;

typedef struct {
  pll_partition_t * partition;
  pll_operation_t * operations;
  int params_index;
  int * matrix_indices;
  double * branch_lengths;
  double alpha;
  int n_cat_gamma;
  int n_prob_matrices;
  int n_operations;
  int is_rooted;
  int freqs_index;
  int clv_index1;
  /* additional stuff for unrooted trees */
  int clv_index2;
  int branch_matrix_index;

  pll_parameter which_parameter;

  double * min_params;
  double * max_params;
} lk_stuff;

double optimize(lk_stuff * tree, double epsilon);

#endif /* OPTIMIZE_H_ */
