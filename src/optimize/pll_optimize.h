/*
    Copyright (C) 2015 Diego Darriba

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

    Contact: Diego Darriba <Diego.Darriba@h-its.org>,
    Exelixis Lab, Heidelberg Instutute for Theoretical Studies
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/
#ifndef PLL_OPTIMIZE_H_
#define PLL_OPTIMIZE_H_

#include <pll.h>
#include <search.h>

/* Parameters mask */

#define PLL_PARAMETER_SUBST_RATES      1
#define PLL_PARAMETER_ALPHA            2
#define PLL_PARAMETER_PINV             4
#define PLL_PARAMETER_FREQUENCIES      8
#define PLL_PARAMETER_BRANCH_LENGTHS  16
#define PLL_PARAMETER_TOPOLOGY        32

/* L-BFGS-B bound type */

#define PLL_LBFGSB_BOUND_NONE  0
#define PLL_LBFGSB_BOUND_LOWER 1
#define PLL_LBFGSB_BOUND_BOTH  2
#define PLL_LBFGSB_BOUND_UPPER 3

/* error codes */

#define PLL_ERROR_TAXA_MISMATCH       101
#define PLL_ERROR_SEQLEN_MISMATCH     102
#define PLL_ERROR_ALIGN_UNREADABLE    103

/* Structure with information necessary for evaluating the likelihood */

typedef struct
{
  pll_partition_t * partition;
  pll_operation_t * operations;
  double * branch_lengths;
  int * matrix_indices;
  int rooted;
  int freqs_index;
  union {
      struct {
        int root_clv_index;
      } rooted_t;
      struct {
        int parent_clv_index;
        int child_clv_index;
        int edge_pmatrix_index;
      } unrooted_t;
    } where;

  double alpha_value;
} pll_likelihood_info_t;

/* Structure with information for parameter optimization */

typedef struct
{
  pll_likelihood_info_t lk_params;
  int params_index;
  unsigned int which_parameters;
  int * subst_params_symmetries;
  double factr;
  double pgtol;
} pll_optimize_options_t;

/* functions in pll_optimize.c */

PLL_EXPORT double pll_optimize_parameters_lbfgsb(pll_optimize_options_t * p);
PLL_EXPORT pll_partition_t * pll_create_partition_fasta (char *file,
                                                         int states,
                                                         int n_rate_matrices,
                                                         int n_rate_cats,
                                                         int attributes,
                                                         int rooted,
                                                         int tip_count,
                                                         char **tipnames);
#endif /* PLL_OPTIMIZE_H_ */
