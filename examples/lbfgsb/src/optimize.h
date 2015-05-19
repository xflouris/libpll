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

#include "pll.h"

enum
{
  PARAM_SUBST_RATES = 1,
  PARAM_PINV = 2,
  PARAM_ALPHA = 4,
  PARAM_FREQUENCIES = 8,
  PARAM_BRANCH_LENGTHS = 16
};

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
  /* tolerances in the stopping criteria. */
  double factr;
  double pgtol;
} opt_params;

double optimize_parameters(opt_params * params);

#endif /* OPTIMIZE_H_ */
