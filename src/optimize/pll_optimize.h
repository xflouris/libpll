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
#ifndef OPTIMIZE_H_
#define OPTIMIZE_H_

#include <pll.h>

#define PLL_PARAMETER_SUBST_RATES      1
#define PLL_PARAMETER_ALPHA            2
#define PLL_PARAMETER_PINV             4
#define PLL_PARAMETER_FREQUENCIES      8
#define PLL_PARAMETER_BRANCH_LENGTHS  16
#define PLL_PARAMETER_TOPOLOGY        32

#define PLL_LBFGSB_BOUND_NONE  0
#define PLL_LBFGSB_BOUND_LOWER 1
#define PLL_LBFGSB_BOUND_BOTH  2
#define PLL_LBFGSB_BOUND_UPPER 3

typedef struct
{
  /* pll stuff */
  pll_partition_t * partition;
  pll_operation_t * operations;
  double * branch_lengths;
  int * matrix_indices;
  int clv1;
  int clv2;
  int edge_pmatrix_index;
  int num_gamma_cats;
  int params_index;
  double alpha_value;

  /* optimization parameters */
  unsigned int which_parameters;

  /* substitution matrix symmetries and parameters */
  int * subst_params_symmetries;

  /* tolerances / stopping criteria. */
  double factr;
  double pgtol;
} pll_optimize_options;

/* functions in pll_optimize.c */

PLL_EXPORT double pll_optimize_parameters_lbfgsb(pll_optimize_options * params);

#endif /* OPTIMIZE_H_ */
