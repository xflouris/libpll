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
#include "pll_optimize.h"
#include "lbfgsb/lbfgsb.h"

#define MAX_VARIABLES 1024

static int v_int_max (int * v, int n)
{
  int i, max = v[0];
  for (i = 1; i < n; i++)
    if (v[i] > max)
      max = v[i];
  return max;
}

static double compute_lnl_unrooted (pll_optimize_options_t * params,
                                    double *x)
{
  pll_partition_t * partition = params->lk_params.partition;
  pll_operation_t * operations = params->lk_params.operations;
  double * branch_lengths = params->lk_params.branch_lengths;
  int * matrix_indices = params->lk_params.matrix_indices;
  int params_index = params->params_index;
  int num_branch_lengths;
  int num_inner_nodes;
  double score;
  double * xptr = x;

  if (params->lk_params.rooted)
  {
    num_branch_lengths = 2 * partition->tips - 2;
    num_inner_nodes = partition->tips - 1;
  }
  else
  {
    num_branch_lengths = 2 * partition->tips - 3;
    num_inner_nodes = partition->tips - 2;
  }

  /* update substitution rate parameters */
  if (params->which_parameters & PLL_PARAMETER_SUBST_RATES)
  {
    int * symm;
    int n_subst_rates;
    double * subst_rates;

    symm = params->subst_params_symmetries;
    n_subst_rates = partition->states * (partition->states - 1) / 2;
    subst_rates = (double *) malloc ((size_t) n_subst_rates * sizeof(double));

    assert(subst_rates);

    if (symm)
    {
      int i, j, k;
      int n_subst_free_params = 0;

      /* compute the number of free parameters */
      n_subst_free_params = v_int_max (symm, n_subst_rates);

      /* assign values to the substitution rates */
      k=0;
      for (i = 0; i <= n_subst_free_params; i++)
      {
        double next_value = (i == symm[n_subst_rates - 1]) ? 1.0 : xptr[k++];
        for (j = 0; j < n_subst_rates; j++)
          if (symm[j] == i)
          {
            subst_rates[j] = next_value;
          }
      }
      pll_set_subst_params (partition, 0, subst_rates);
      xptr += n_subst_free_params;
    }
    else
    {
      memcpy (subst_rates, xptr, (n_subst_rates - 1) * sizeof(double));
      subst_rates[n_subst_rates - 1] = 1.0;
    }
    free (subst_rates);
  }

  /* update proportion of invariant parameters */
  if (params->which_parameters & PLL_PARAMETER_PINV)
  {
    pll_update_invariant_sites_proportion (partition, params_index, xptr[0]);
    xptr++;
  }
  if (params->which_parameters & PLL_PARAMETER_ALPHA)
  {
    /* assign discrete rates */
    double * rate_cats;
    rate_cats = alloca(partition->rate_cats);
    params->lk_params.alpha_value = xptr[0];
    pll_compute_gamma_cats (xptr[0], partition->rate_cats, rate_cats);
    pll_set_category_rates (partition, rate_cats);
    xptr++;
  }
  if (params->which_parameters & PLL_PARAMETER_BRANCH_LENGTHS)
  {
    /* assign branch lengths */
    memcpy (branch_lengths, xptr, num_branch_lengths * sizeof(double));
    xptr += num_branch_lengths;
  }

  pll_update_prob_matrices (partition, 0,
                            matrix_indices,
                            branch_lengths,
                            num_branch_lengths);

  pll_update_partials (partition,
                       operations,
                       num_inner_nodes);

  if (params->lk_params.rooted)
  {
    score = -1
        * pll_compute_root_loglikelihood (
            partition, params->lk_params.where.rooted_t.root_clv_index,
            params->lk_params.freqs_index);
  }
  else
  {
    score = -1
        * pll_compute_edge_loglikelihood (
            partition, params->lk_params.where.unrooted_t.parent_clv_index,
            params->lk_params.where.unrooted_t.child_clv_index,
            params->lk_params.where.unrooted_t.edge_pmatrix_index,
            params->lk_params.freqs_index);
  }
  return score;
}

PLL_EXPORT double pll_optimize_parameters_lbfgsb (
    pll_optimize_options_t * params)
{
  int i;
  pll_partition_t * partition = params->lk_params.partition;
  /* L-BFGS-B */

  /* Local variables */
  double score, g[MAX_VARIABLES];
  double lower_bounds[MAX_VARIABLES];
  int m, num_variables;
  double upper_bounds[MAX_VARIABLES], x[MAX_VARIABLES], wa[43251];
  int bound_type[MAX_VARIABLES], iwa[3072];

  /*     static char task[60]; */
  int taskValue;
  int *task = &taskValue; /* must initialize !! */

  /*     static char csave[60]; */
  int csaveValue;
  int *csave = &csaveValue;
  double dsave[29];
  int isave[44];
  logical lsave[4];

  int iprint = -1;

  /*     We specify the dimension n of the sample problem and the number */
  /*        m of limited memory corrections stored.  (n and m should not */
  /*        exceed the limits nmax and mmax respectively.) */

  m = 5;
  num_variables = 0;

  /*     We now provide nbd which defines the bounds on the variables: */
  /*                    l   specifies the lower bounds, */
  /*                    u   specifies the upper bounds. */
  int * nbd_ptr = bound_type;
  double * l_ptr = lower_bounds, *u_ptr = upper_bounds;
  if (params->which_parameters & PLL_PARAMETER_SUBST_RATES)
  {
    int n_subst_rates;
    int n_subst_free_params;

    n_subst_rates = partition->states * (partition->states - 1) / 2;
    if (params->subst_params_symmetries)
      n_subst_free_params = v_int_max (params->subst_params_symmetries,
                                       n_subst_rates);
    else
      n_subst_free_params = n_subst_rates;

    for (i = 0; i < n_subst_free_params; i++)
    {
      nbd_ptr[i] = PLL_LBFGSB_BOUND_BOTH;
      l_ptr[i] = 0.001;
      u_ptr[i] = 1000.;
      x[i] = partition->subst_params[params->params_index][i];
    }
    nbd_ptr += n_subst_free_params;
    l_ptr += n_subst_free_params;
    u_ptr += n_subst_free_params;
    num_variables += n_subst_free_params;
  }
  if (params->which_parameters & PLL_PARAMETER_PINV)
  {
    *nbd_ptr = PLL_LBFGSB_BOUND_BOTH;
    *l_ptr = 0.0;
    *u_ptr = 0.99;
    x[num_variables] = partition->prop_invar[params->params_index];
    num_variables += 1;
    nbd_ptr++;
    l_ptr++;
    u_ptr++;
  }
  if (params->which_parameters & PLL_PARAMETER_ALPHA)
  {
    *nbd_ptr = PLL_LBFGSB_BOUND_BOTH;
    *l_ptr = 0.02;
    *u_ptr = 100.0;
    x[num_variables] = params->lk_params.alpha_value;
    num_variables += 1;
    nbd_ptr++;
    l_ptr++;
    u_ptr++;
  }
  if (params->which_parameters
      & (PLL_PARAMETER_FREQUENCIES | PLL_PARAMETER_TOPOLOGY))
  {
    return PLL_FAILURE;
  }
  if (params->which_parameters & PLL_PARAMETER_BRANCH_LENGTHS)
  {
    int num_branch_lengths = params->lk_params.rooted?
        (2 * partition->tips - 3) :
        (2 * partition->tips - 2);
    for (i = 0; i < num_branch_lengths; i++)
    {
      nbd_ptr[i] = PLL_LBFGSB_BOUND_LOWER;
      l_ptr[i] = 0.00001;
      x[num_variables + i] = params->lk_params.branch_lengths[i];
    }
    num_variables += num_branch_lengths;
    nbd_ptr += num_branch_lengths;
    l_ptr += num_branch_lengths;
    u_ptr += num_branch_lengths;
  }

  /*     We start the iteration by initializing task. */
  *task = (int) START;

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
