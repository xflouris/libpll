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
#include "optimize.h"
#include "lbfgsb.h"

#define MAX_VARIABLES 1024

static double compute_lnl_unrooted (opt_params * params, double *x)
{
  pll_partition_t * partition = params->partition;
  pll_operation_t * operations = params->operations;
  double * branch_lengths = params->branch_lengths;
  int * matrix_indices = params->matrix_indices;
  int clv1 = params->clv1;
  int clv2 = params->clv2;
  int edge_pmatrix_index = params->edge_pmatrix_index;
  int num_gamma_cats = params->num_gamma_cats;
  int params_index = params->params_index;

  double score;
  double * xptr = x;

  if (params->which_parameters & PARAM_SUBST_RATES)
  {
    double subsparams[6];
    memcpy (subsparams, xptr, 5 * sizeof(double));
    subsparams[5] = 1.0;
    pll_set_subst_params (partition, 0, subsparams);
    xptr += 5;
  }
  if (params->which_parameters & PARAM_PINV)
  {
    pll_update_invariant_sites_proportion (partition, params_index, xptr[0]);
    xptr++;
  }
  if (params->which_parameters & PARAM_ALPHA)
  {
    double * rate_cats;
    rate_cats = alloca(params->num_gamma_cats);
    params->alpha_value = xptr[0];
    pll_compute_gamma_cats (xptr[0], num_gamma_cats, rate_cats);
    pll_set_category_rates (partition, rate_cats);
    xptr++;
  }
  if (params->which_parameters & PARAM_FREQUENCIES)
  {
    printf ("We do not know how to do this!");
    assert(0);
  }
  if (params->which_parameters & PARAM_BRANCH_LENGTHS)
  {
    int num_branch_lengths = 2 * partition->tips - 3;
    memcpy (branch_lengths, xptr, num_branch_lengths * sizeof(double));
    xptr += num_branch_lengths;
  }

  pll_update_prob_matrices (partition, 0, matrix_indices, branch_lengths,
                            2 * partition->tips - 3);

  pll_update_partials (partition, operations, partition->tips - 2);
  score = -1
      * pll_compute_edge_loglikelihood (partition, clv1, clv2,
                                        edge_pmatrix_index, 0);
  return score;
}

double optimize_parameters (opt_params * params)
{
  int i;

  /* L-BFGS-B */

  /* Local variables */
  double score, g[MAX_VARIABLES];
  double lower_bounds[MAX_VARIABLES];
  integer m, num_variables;
  double upper_bounds[MAX_VARIABLES], x[MAX_VARIABLES], wa[43251];
  integer bound_type[MAX_VARIABLES], iwa[3072];

  /*     static char task[60]; */
  integer taskValue;
  integer *task = &taskValue; /* must initialize !! */

  /*     static char csave[60]; */
  integer csaveValue;
  integer *csave = &csaveValue;
  double dsave[29];
  integer isave[44];
  logical lsave[4];

  /*     We wish to have output at every iteration. */
  integer iprint = 0;

  /*     We specify the dimension n of the sample problem and the number */
  /*        m of limited memory corrections stored.  (n and m should not */
  /*        exceed the limits nmax and mmax respectively.) */

  m = 5;
  num_variables = 0;

  /*     We now provide nbd which defines the bounds on the variables: */
  /*                    l   specifies the lower bounds, */
  /*                    u   specifies the upper bounds. */
  /*         nbd(i)=0 if x(i) is unbounded, */
  /*                1 if x(i) has only a lower bound, */
  /*                2 if x(i) has both lower and upper bounds, and */
  /*                3 if x(i) has only an upper bound. */
  integer * nbd_ptr = bound_type;
  double * l_ptr = lower_bounds, *u_ptr = upper_bounds;
  if (params->which_parameters & PARAM_SUBST_RATES)
  {
    num_variables = 5;
    for (i = 0; i < num_variables; i++)
    {
      nbd_ptr[i] = 2;
      l_ptr[i] = 0.0001;
      u_ptr[i] = 1000.;
      x[i] = params->partition->subst_params[params->params_index][i];
    }
    nbd_ptr += 5;
    l_ptr += 5;
    u_ptr += 5;
  }
  if (params->which_parameters & PARAM_PINV)
  {
    *nbd_ptr = 2;
    *l_ptr = 0.0;
    *u_ptr = 0.99;
    x[num_variables] = params->partition->prop_invar[params->params_index];
    num_variables += 1;
    nbd_ptr++;
    l_ptr++;
    u_ptr++;
  }
  if (params->which_parameters & PARAM_ALPHA)
  {
    *nbd_ptr = 2;
    *l_ptr = 0.02;
    *u_ptr = 100.0;
    x[num_variables] = params->alpha_value;
    num_variables += 1;
    nbd_ptr++;
    l_ptr++;
    u_ptr++;
  }
  if (params->which_parameters & PARAM_FREQUENCIES)
  {
    printf ("We do not know how to do this!");
    assert(0);
  }
  if (params->which_parameters & PARAM_BRANCH_LENGTHS)
  {
    int num_branch_lengths = 2 * params->partition->tips - 3;
    for (i = 0; i < num_branch_lengths; i++)
    {
      nbd_ptr[i] = 1; /* lower bound only */
      l_ptr[i] = 0.00001;
      x[num_variables + i] = params->branch_lengths[i];
    }
    num_variables += num_branch_lengths;
    nbd_ptr += num_branch_lengths;
    l_ptr += num_branch_lengths;
    u_ptr += num_branch_lengths;
  }

  /*     We start the iteration by initializing task. */

  *task = (integer) START;

  int continue_opt = 1;
  while (continue_opt)
  {
    /*     This is the call to the L-BFGS-B code. */
    setulb (&num_variables, &m, x, lower_bounds, upper_bounds, bound_type,
            &score, g, &(params->factr), &(params->pgtol), wa, iwa, task,
            &iprint, csave, lsave, isave, dsave);
    if (IS_FG(*task))
    {
      /*        the minimization routine has returned to request the */
      /*        function f and gradient g values at the current x. */
      /*        Compute function value f for the sample problem. */

      score = compute_lnl_unrooted (params, x);

      double h, temp;
      for (i = 0; i < num_variables; i++)
      {
        double ERROR_X = 1.0e-4;
        temp = x[i];
        h = ERROR_X * fabs (temp);
        if (h < 1e-12)
          h = ERROR_X;

        x[i] = temp + h;
        h = x[i] - temp;
        double lnderiv = compute_lnl_unrooted (params, x);

        g[i] = (lnderiv - score) / h;

        /* reset variable */
        x[i] = temp;
      }
    }
    else if (*task != NEW_X)
    {
      continue_opt = 0;
    }
  }

  return score;
}
