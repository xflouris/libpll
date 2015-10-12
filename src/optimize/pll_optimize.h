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
//#include <search.h>

/* Parameters mask */
#define PLL_PARAMETER_SUBST_RATES          1
#define PLL_PARAMETER_ALPHA                2
#define PLL_PARAMETER_PINV                 4
#define PLL_PARAMETER_FREQUENCIES          8
#define PLL_PARAMETER_BRANCHES_SINGLE     16
#define PLL_PARAMETER_BRANCHES_ALL        32
#define PLL_PARAMETER_BRANCHES_ITERATIVE  64
#define PLL_PARAMETER_TOPOLOGY           128

/* L-BFGS-B bound type */
#define PLL_LBFGSB_BOUND_NONE  0
#define PLL_LBFGSB_BOUND_LOWER 1
#define PLL_LBFGSB_BOUND_BOTH  2
#define PLL_LBFGSB_BOUND_UPPER 3

/* Parameter defaults and limits */
#define PLL_OPT_MIN_BRANCH_LEN       1.0e-4
#define PLL_OPT_MAX_BRANCH_LEN          100
#define PLL_OPT_TOL_BRANCH_LEN       1.0e-4
#define PLL_OPT_DEFAULT_RATE_RATIO        1
#define PLL_OPT_DEFAULT_FREQ_RATIO        1
#define PLL_OPT_DEFAULT_PINV            0.5
#define PLL_OPT_DEFAULT_ALPHA           0.5
#define PLL_OPT_DEFAULT_BRANCH_LEN      0.1
#define PLL_OPT_LNL_UNLIKELY         -1e+80

/* Branch lengths optimization algorithm */
#define PLL_BRANCH_OPT_NEWTON 1
#define PLL_BRANCH_OPT_BRENT  2
#define PLL_BRANCH_OPT_LBFGSB 3

/* error codes */
#define PLL_ERROR_PARAMETER           100
#define PLL_ERROR_TAXA_MISMATCH       101
#define PLL_ERROR_SEQLEN_MISMATCH     102
#define PLL_ERROR_ALIGN_UNREADABLE    103
#define PLL_ERROR_LBFGSB_UNKNOWN      110
#define PLL_ERROR_NEWTON_DERIV        120
#define PLL_ERROR_NEWTON_LIMIT        130
#define PLL_ERROR_NEWTON_UNKNOWN      140

/* Structure with information necessary for evaluating the likelihood */

typedef struct
{
  pll_partition_t * partition;
  pll_operation_t * operations;
  double * branch_lengths;
  unsigned int * matrix_indices;
  int rooted;
  unsigned int freqs_index;
  union {
      struct {
        unsigned int root_clv_index;
        int scaler_index;
      } rooted_t;
      struct {
        unsigned int parent_clv_index;
        int parent_scaler_index;
        unsigned int child_clv_index;
        int child_scaler_index;
        unsigned int edge_pmatrix_index;
      } unrooted_t;
    } where;

  char __padding__[4];
  double alpha_value;
} pll_likelihood_info_t;

/* Structure with information for parameter optimization */

typedef struct
{
  pll_likelihood_info_t lk_params;
  unsigned int highest_freq_state;
  unsigned int params_index;
  unsigned int mixture_index;
  unsigned int which_parameters;
  int * subst_params_symmetries;
  double factr;
  double pgtol;

  double * sumtable;
} pll_optimize_options_t;

 typedef struct
{
  pll_partition_t * partition;
  pll_utree_t * tree;
  unsigned int params_index;
  unsigned int freqs_index;
  double * sumtable;
} pll_newton_tree_params_t;

/* functions in pll_optimize.c */
PLL_EXPORT double * pll_compute_empirical_frequencies(pll_partition_t * partition);
PLL_EXPORT double * pll_compute_empirical_subst_rates(pll_partition_t * partition);

PLL_EXPORT double pll_derivative_func(void * parameters,
                                      double proposal,
                                      double *df, double *ddf);

/* core Newton-Raphson optimization function */
PLL_EXPORT double pll_minimize_newton(double x1,
                                      double xguess,
                                      double x2,
                                      unsigned int max_iters,
                                      double *score,
                                      void * params,
                                      double (deriv_func)(void *,
                                                          double,
                                                          double *, double *));

/* core Brent optimization function */
PLL_EXPORT double pll_minimize_brent(double xmin,
                                     double xguess,
                                     double xmax,
                                     double xtol,
                                     double *fx,
                                     double *f2x,
                                     void * params,
                                     double (*target_funk)(
                                         void *,
                                         double));

/* optimization functions */
PLL_EXPORT double pll_optimize_parameters_brent(pll_optimize_options_t * p);
PLL_EXPORT double pll_optimize_parameters_lbfgsb(pll_optimize_options_t * p);
PLL_EXPORT double pll_optimize_branch_lengths_iterative (
                                              pll_partition_t * partition,
                                              pll_utree_t * tree,
                                              unsigned int params_index,
                                              unsigned int freqs_index,
                                              double tolerance,
                                              int smoothings);

#endif /* PLL_OPTIMIZE_H_ */
