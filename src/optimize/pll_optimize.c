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

#define UPDATE_SCALERS 1
#define BL_OPT_METHOD PLL_BRANCH_OPT_NEWTON

static int is_nan(double v)
{
  return v!=v;
}

static int v_int_max (int * v, int n)
{
  int i, max = v[0];
  for (i = 1; i < n; i++)
    if (v[i] > max)
      max = v[i];
  return max;
}

static int d_equals(double a, double b)
{
  return (fabs(a-b) < 1e-10);
}

static int set_x_to_parameters(pll_optimize_options_t * params,
                               double *x)
{
  pll_partition_t * partition = params->lk_params.partition;
  pll_operation_t * operations = params->lk_params.operations;
  double * branch_lengths = params->lk_params.branch_lengths;
  unsigned int * matrix_indices = params->lk_params.matrix_indices;
  unsigned int params_index = params->params_index;
  unsigned int n_branches, n_inner_nodes;
  double * xptr = x;

  if (params->lk_params.rooted)
  {
    n_branches = 2 * partition->tips - 2;
    n_inner_nodes = partition->tips - 1;
  }
  else
  {
    n_branches = 2 * partition->tips - 3;
    n_inner_nodes = partition->tips - 2;
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
      k = 0;
      for (i = 0; i <= n_subst_free_params; i++)
      {
        double next_value = (i == symm[n_subst_rates - 1]) ? 1.0 : xptr[k++];
        for (j = 0; j < n_subst_rates; j++)
          if (symm[j] == i)
          {
            subst_rates[j] = next_value;
          }
      }
      xptr += n_subst_free_params;
    }
    else
    {
      memcpy (subst_rates, xptr, ((size_t)n_subst_rates - 1) * sizeof(double));
      subst_rates[n_subst_rates - 1] = 1.0;
      xptr += n_subst_rates-1;
    }

    pll_set_subst_params (partition, 0, subst_rates);
    free (subst_rates);
  }
  /* update stationary frequencies */
  if (params->which_parameters & PLL_PARAMETER_FREQUENCIES)
  {
    unsigned int i;
    unsigned int n_states = partition->states;
    double sum_ratios = 1.0;
    double *freqs = (double *) malloc ((size_t) n_states * sizeof(double));
    for (i = 0; i < (n_states - 1); ++i)
    {
      assert(!is_nan(xptr[i]));
      sum_ratios += xptr[i];
    }
    for (i = 0; i < (n_states - 1); ++i)
      freqs[i] = xptr[i] / sum_ratios;
    freqs[n_states - 1] = 1.0 / sum_ratios;
    pll_set_frequencies (partition, params->params_index, freqs);
    free (freqs);
    xptr += (n_states - 1);
  }
  /* update proportion of invariant sites */
  if (params->which_parameters & PLL_PARAMETER_PINV)
  {
    assert(!is_nan(xptr[0]));
    if (!pll_update_invariant_sites_proportion (partition,
                                                params_index,
                                                xptr[0]))
    {
      return PLL_FAILURE;
    }
    xptr++;
  }
  /* update gamma shape parameter */
  if (params->which_parameters & PLL_PARAMETER_ALPHA)
  {
    assert(!is_nan(xptr[0]));
    /* assign discrete rates */
    double * rate_cats;
    rate_cats = malloc((size_t)partition->rate_cats * sizeof(double));
    params->lk_params.alpha_value = xptr[0];
    if (!pll_compute_gamma_cats (xptr[0], partition->rate_cats, rate_cats))
    {
      return PLL_FAILURE;
    }
    pll_set_category_rates (partition, rate_cats);

    free(rate_cats);
    xptr++;
  }
  /* update all branch lengths */
  if (params->which_parameters & PLL_PARAMETER_BRANCHES_ALL)
  {
    /* assign branch lengths */
    memcpy (branch_lengths, xptr, (size_t)n_branches * sizeof(double));
    xptr += n_branches;
  }
  /* update single branch */
  if (params->which_parameters & PLL_PARAMETER_BRANCHES_SINGLE)
   {
    assert(!is_nan(xptr[0]));
     /* assign branch length */
     *branch_lengths = *xptr;
     pll_update_prob_matrices (partition,
                         0,
                         &params->lk_params.where.unrooted_t.edge_pmatrix_index,
                         xptr,
                         1);
     xptr += 1;
   }
  else
   {
       pll_update_prob_matrices (partition, 0, matrix_indices, branch_lengths,
                                 n_branches);

       pll_update_partials (partition, operations, n_inner_nodes);
   }
  return PLL_SUCCESS;
}

static double compute_lnl_unrooted (pll_optimize_options_t * params, double *x)
{
  pll_partition_t * partition = params->lk_params.partition;
  double score;

  if (x && !set_x_to_parameters(params, x))
    return -INFINITY;

  if (params->lk_params.rooted)
  {
    score = -1
        * pll_compute_root_loglikelihood (
            partition, params->lk_params.where.rooted_t.root_clv_index,
            params->lk_params.where.rooted_t.scaler_index,
            params->lk_params.freqs_index);
  }
  else
  {
    score = -1
        * pll_compute_edge_loglikelihood (
            partition, params->lk_params.where.unrooted_t.parent_clv_index,
            params->lk_params.where.unrooted_t.parent_scaler_index,
            params->lk_params.where.unrooted_t.child_clv_index,
            params->lk_params.where.unrooted_t.child_scaler_index,
            params->lk_params.where.unrooted_t.edge_pmatrix_index,
            params->lk_params.freqs_index);
  }

  return score;
} /* compute_lnl_unrooted */

static unsigned int count_n_free_variables (pll_optimize_options_t * params)
{
  unsigned int num_variables = 0;
  pll_partition_t * partition = params->lk_params.partition;

  /* count number of variables for dynamic allocation */
  if (params->which_parameters & PLL_PARAMETER_SUBST_RATES)
  {
    int n_subst_rates = partition->states * (partition->states - 1) / 2;
    num_variables +=
        params->subst_params_symmetries ?
            (unsigned int) v_int_max (params->subst_params_symmetries, n_subst_rates) :
            (unsigned int) n_subst_rates - 1;
  }
  if (params->which_parameters & PLL_PARAMETER_FREQUENCIES)
    num_variables += partition->states - 1;
  num_variables += (params->which_parameters & PLL_PARAMETER_PINV) != 0;
  num_variables += (params->which_parameters & PLL_PARAMETER_ALPHA) != 0;
  num_variables += (params->which_parameters & PLL_PARAMETER_BRANCHES_SINGLE)
      != 0;
  if (params->which_parameters & PLL_PARAMETER_BRANCHES_ALL)
  {
    unsigned int num_branch_lengths =
        params->lk_params.rooted ?
            (2 * partition->tips - 3) : (2 * partition->tips - 2);
    num_variables += num_branch_lengths;
  }
  return num_variables;
} /* count_n_free_variables */


/******************************************************************************/
/* BRENT'S OPTIMIZATION */
/******************************************************************************/

#define ITMAX 100
#define CGOLD 0.3819660
#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define ZEPS 1.0e-10
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

static double brent_opt (double ax, double bx, double cx, double tol,
                             double *foptx, double *f2optx, double fax,
                             double fbx, double fcx,
                             pll_optimize_options_t * params,
                             double (*target_funk)(
                                 pll_optimize_options_t *,
                                 double))
{
  int iter;
  double a, b, d = 0, etemp, fu, fv, fw, fx, p, q, r, tol1, tol2, u, v, w, x,
      xm;
  double xw, wv, vx;
  double e = 0.0;

  a = (ax < cx ? ax : cx);
  b = (ax > cx ? ax : cx);
  x = bx;
  fx = fbx;
  if (fax < fcx)
  {
    w = ax;
    fw = fax;
    v = cx;
    fv = fcx;
  }
  else
  {
    w = cx;
    fw = fcx;
    v = ax;
    fv = fax;
  }

  for (iter = 1; iter <= ITMAX; iter++)
  {
    xm = 0.5 * (a + b);
    tol2 = 2.0 * (tol1 = tol * fabs (x) + ZEPS);
    if (fabs (x - xm) <= (tol2 - 0.5 * (b - a)))
    {
      *foptx = fx;
      xw = x - w;
      wv = w - v;
      vx = v - x;
      *f2optx = 2.0 * (fv * xw + fx * wv + fw * vx)
          / (v * v * xw + x * x * wv + w * w * vx);
      return x;
    }

    if (fabs (e) > tol1)
    {
      r = (x - w) * (fx - fv);
      q = (x - v) * (fx - fw);
      p = (x - v) * q - (x - w) * r;
      q = 2.0 * (q - r);
      if (q > 0.0)
        p = -p;
      q = fabs (q);
      etemp = e;
      e = d;
      if (fabs (p) >= fabs (0.5 * q * etemp) || p <= q * (a - x)
          || p >= q * (b - x))
        d = CGOLD * (e = (x >= xm ? a - x : b - x));
      else
      {
        d = p / q;
        u = x + d;
        if (u - a < tol2 || b - u < tol2)
          d = SIGN(tol1, xm - x);
      }
    }
    else
    {
      d = CGOLD * (e = (x >= xm ? a - x : b - x));
    }

    u = (fabs (d) >= tol1 ? x + d : x + SIGN(tol1, d));
    fu = target_funk (params, u);
    if (fu <= fx)
    {
      if (u >= x)
        a = x;
      else
        b = x;

      SHFT(v, w, x, u)
      SHFT(fv, fw, fx, fu)
    }
    else
    {
      if (u < x)
        a = u;
      else
        b = u;
      if (fu <= fw || w == x)
      {
        v = w;
        w = u;
        fv = fw;
        fw = fu;
      }
      else if (fu <= fv || v == x || v == w)
      {
        v = u;
        fv = fu;
      }
    }
  }

  *foptx = fx;
  xw = x - w;
  wv = w - v;
  vx = v - x;
  *f2optx = 2.0 * (fv * xw + fx * wv + fw * vx)
      / (v * v * xw + x * x * wv + w * w * vx);

  return x;
}

static double brent_target(pll_optimize_options_t * p, double x)
{
  return compute_lnl_unrooted(p, &x);
}

/* most of the code for Brent optimization taken from IQ-Tree
 * http://www.cibiv.at/software/iqtree
 * --------------------------------------------------------------------------
 * Bui Quang Minh, Minh Anh Thi Nguyen, and Arndt von Haeseler (2013)
 * Ultrafast approximation for phylogenetic bootstrap.
 * Mol. Biol. Evol., 30:1188-1195. (free reprint, DOI: 10.1093/molbev/mst024) */
PLL_EXPORT double pll_minimize_brent(pll_optimize_options_t * params,
                                     double xmin,
                                     double xguess,
                                     double xmax,
                                     double *fx,
                                     double *f2x,
                                     double (*target_funk)(
                                         pll_optimize_options_t *,
                                         double))
{
  double eps, optx, ax, bx, cx, fa, fb, fc;
  int outbounds_ax, outbounds_cx;
  double tolerance = params->pgtol;

  /* first attempt to bracketize minimum */
  if (xguess < xmin)
    xguess = xmin;
  if (xguess > xmax)
    xguess = xmax;
  eps = xguess * tolerance * 50.0;
  ax = xguess - eps;
  outbounds_ax = ax < xmin;
  if (outbounds_ax)
    ax = xmin;
  bx = xguess;
  cx = xguess + eps;
  outbounds_cx = cx > xmax;
  if (outbounds_cx)
    cx = xmax;

  /* check if this works */
  fa = target_funk (params, ax);
  fb = target_funk (params, bx);
  fc = target_funk (params, cx);

  /* if it works use these borders else be conservative */
  if ((fa < fb) || (fc < fb))
  {
    if (!outbounds_ax)
      fa = target_funk (params, xmin);
    if (!outbounds_cx)
      fc = target_funk (params, xmax);
    ax = xmin;
    cx = xmax;
  }

  optx = brent_opt (ax, bx, cx, tolerance, fx, f2x, fa, fb, fc, params,
                    target_funk);
  if (*fx > fb) // if worse, return initial value
  {
    *fx = target_funk (params, bx);
    return bx;
  }

  return optx; /* return optimal x */
}

PLL_EXPORT double pll_optimize_parameters_brent(pll_optimize_options_t * params)
{
  double score = 0;

  /* Brent parameters */
  double xmin;
  double xguess;
  double xmax;
  double f2x;

  switch (params->which_parameters)
  {
    case PLL_PARAMETER_ALPHA:
      xmin = 0.0201 + PLL_LBFGSB_ERROR;
      xmax = 100.0;
      xguess = params->lk_params.alpha_value;
      if (xguess < xmin || xguess > xmax)
         xguess = PLL_OPT_DEFAULT_ALPHA;
      break;
    case PLL_PARAMETER_PINV:
      xmin = PLL_LBFGSB_ERROR;
      xmax = 0.99;
      xguess = params->lk_params.partition->prop_invar[params->params_index];
      if (xguess < xmin || xguess > xmax)
         xguess = PLL_OPT_DEFAULT_PINV;
      break;
    case PLL_PARAMETER_BRANCHES_SINGLE:
      xmin = PLL_OPT_MIN_BRANCH_LEN + PLL_LBFGSB_ERROR;
      xmax = PLL_OPT_MAX_BRANCH_LEN;
      xguess = params->lk_params.branch_lengths[0];
      if (xguess < xmin || xguess > xmax)
        xguess = PLL_OPT_DEFAULT_BRANCH_LEN;
//      //TODO: Temporarily unavailable
//      assert(0);
      break;
    default:
      /* unavailable or multiple parameter */
      return -INFINITY;
  }

  double xres = pll_minimize_brent(params,
                                   xmin, xguess, xmax,
                                   &score, &f2x,
                                   &brent_target);
  score = brent_target(params, xres);
  return score;
}

PLL_EXPORT double pll_optimize_parameters_newton(pll_optimize_options_t * params)
{
  double score = 0.0;

  /* Brent parameters */
  double xmin;
  double xguess;
  double xmax;

  assert(params->which_parameters == PLL_PARAMETER_BRANCHES_SINGLE);

  xmin = PLL_OPT_MIN_BRANCH_LEN + PLL_LBFGSB_ERROR;
  xmax = PLL_OPT_MAX_BRANCH_LEN;
  xguess = params->lk_params.branch_lengths[0];
  if (xguess < xmin || xguess > xmax)
    xguess = PLL_OPT_DEFAULT_BRANCH_LEN;

  double xres = pll_minimize_newton(params,
                                   xmin, xguess, xmax,
                                   200,
                                   &score);

  if (pll_errno)
  {
    printf("ERROR: %s\n", pll_errmsg);
    exit(1);
  }

  params->lk_params.branch_lengths[0] = xres;
  pll_update_prob_matrices(params->lk_params.partition,
                           params->params_index,
                           &(params->lk_params.where.unrooted_t.edge_pmatrix_index),
                           &xres,1);

  if (score < 0)
    score *= -1;
  else
  {
   score = brent_target(params, xres);
  }
  return score;
}

/******************************************************************************/
/* L-BFGS-B OPTIMIZATION */
/******************************************************************************/

PLL_EXPORT double pll_optimize_parameters_lbfgsb (
                                              pll_optimize_options_t * params)
{
  unsigned int i;
  pll_partition_t * partition = params->lk_params.partition;

  /* L-BFGS-B parameters */
  int max_corrections;
  unsigned int num_variables;
  double score = 0;
  double *x, *g, *lower_bounds, *upper_bounds, *wa;
  int *bound_type, *iwa;

  int taskValue;
  int *task = &taskValue;

  int csaveValue;
  int *csave = &csaveValue;
  double dsave[29];
  int isave[44];
  logical lsave[4];

  int iprint = -1;

  pll_errno = 0;

  /* ensure that the 2 branch optimization modes are not set together */
  assert(!((params->which_parameters & PLL_PARAMETER_BRANCHES_ALL)
      && (params->which_parameters & PLL_PARAMETER_BRANCHES_SINGLE)));

  max_corrections = 5;
  num_variables = count_n_free_variables (params);

  x = (double *) calloc ((size_t) num_variables, sizeof(double));
  g = (double *) calloc ((size_t) num_variables, sizeof(double));
  lower_bounds = (double *) calloc ((size_t) num_variables, sizeof(double));
  upper_bounds = (double *) calloc ((size_t) num_variables, sizeof(double));
  bound_type = (int *) calloc ((size_t) num_variables, sizeof(int));

  {
    int * nbd_ptr = bound_type;
    double * l_ptr = lower_bounds, *u_ptr = upper_bounds;
    unsigned int check_n = 0;

    /* substitution rate parameters */
    if (params->which_parameters & PLL_PARAMETER_SUBST_RATES)
    {
      unsigned int n_subst_rates;
      unsigned int n_subst_free_params;

      n_subst_rates = partition->states * (partition->states - 1) / 2;
      if (params->subst_params_symmetries)
      {
        n_subst_free_params =(unsigned int) v_int_max (
                                         params->subst_params_symmetries,
                                         (int) n_subst_rates);
      }
      else
      {
        n_subst_free_params = n_subst_rates;
      }

      int current_rate = 0;
      for (i = 0; i < n_subst_free_params; i++)
      {
        nbd_ptr[i] = PLL_LBFGSB_BOUND_BOTH;
        l_ptr[i] = 0.001;
        u_ptr[i] = 100.;
        unsigned int j = i;
        if (params->subst_params_symmetries)
        {
          if (params->subst_params_symmetries[n_subst_rates-1] == current_rate)
            current_rate++;
          for (j=0; j<n_subst_rates; j++)
          {
            if (params->subst_params_symmetries[j] == current_rate)
              break;
          }
          current_rate++;
        }
        x[check_n + i] = partition->subst_params[params->params_index][j];

        if (x[check_n + i] < *l_ptr || x[check_n + i] > *u_ptr)
                                   x[check_n + 1] = PLL_OPT_DEFAULT_RATE_RATIO;
      }
      nbd_ptr += n_subst_free_params;
      l_ptr += n_subst_free_params;
      u_ptr += n_subst_free_params;
      check_n += n_subst_free_params;
    }

    /* stationary frequency parameters */
    if (params->which_parameters & PLL_PARAMETER_FREQUENCIES)
    {
      unsigned int n_freqs_free_params;

      n_freqs_free_params = params->lk_params.partition->states - 1;
      for (i = 0; i < n_freqs_free_params; i++)
      {
        nbd_ptr[i] = PLL_LBFGSB_BOUND_BOTH;
        l_ptr[i] = PLL_LBFGSB_ERROR;
        u_ptr[i] = 1000;
        x[check_n + i] = params->freq_ratios[i];

        if (x[check_n + i] < *l_ptr || x[check_n + i] > *u_ptr)
                              x[check_n + 1] = PLL_OPT_DEFAULT_FREQ_RATIO;
      }
      check_n += n_freqs_free_params;
      nbd_ptr += n_freqs_free_params;
      l_ptr += n_freqs_free_params;
      u_ptr += n_freqs_free_params;
    }

    /* proportion of invariant sites */
    if (params->which_parameters & PLL_PARAMETER_PINV)
    {
      *nbd_ptr = PLL_LBFGSB_BOUND_BOTH;
      *l_ptr = PLL_LBFGSB_ERROR;
      *u_ptr = 0.99;
      x[check_n] = partition->prop_invar[params->params_index];

      if (x[check_n] < *l_ptr || x[check_n] > *u_ptr)
                      x[check_n] = PLL_OPT_DEFAULT_PINV;

      check_n += 1;
      nbd_ptr++;
      l_ptr++;
      u_ptr++;
    }

    /* gamma shape parameter */
    if (params->which_parameters & PLL_PARAMETER_ALPHA)
    {
      *nbd_ptr = PLL_LBFGSB_BOUND_BOTH;
      /* minimum alpha + error offset */
      *l_ptr = 0.0201 + PLL_LBFGSB_ERROR;
      *u_ptr = 100.0;
      x[check_n] = params->lk_params.alpha_value;

      if (x[check_n] < *l_ptr || x[check_n] > *u_ptr)
                x[check_n] = PLL_OPT_DEFAULT_ALPHA;

      check_n += 1;
      nbd_ptr++;
      l_ptr++;
      u_ptr++;
    }

    /* topology (UNIMPLEMENTED) */
    if (params->which_parameters & PLL_PARAMETER_TOPOLOGY)
    {
      return PLL_FAILURE;
    }

    /* single branch length */
    if (params->which_parameters & PLL_PARAMETER_BRANCHES_SINGLE)
    {
        *nbd_ptr = PLL_LBFGSB_BOUND_LOWER;
        *l_ptr = PLL_OPT_MIN_BRANCH_LEN + PLL_LBFGSB_ERROR;
        x[check_n] = params->lk_params.branch_lengths[0];

        if (x[check_n] < *l_ptr)
          x[check_n] = PLL_OPT_DEFAULT_BRANCH_LEN;

        check_n += 1;
        nbd_ptr++;
        l_ptr++;
        u_ptr++;
    }

    /* all branches */
    if (params->which_parameters & PLL_PARAMETER_BRANCHES_ALL)
    {
      unsigned int num_branch_lengths =
          params->lk_params.rooted ?
              (2 * partition->tips - 3) : (2 * partition->tips - 2);
      for (i = 0; i < num_branch_lengths; i++)
      {
        nbd_ptr[i] = PLL_LBFGSB_BOUND_LOWER;
        l_ptr[i] = PLL_OPT_MIN_BRANCH_LEN + PLL_LBFGSB_ERROR;
        x[check_n + i] =
            params->lk_params.branch_lengths[i] > l_ptr[i]?
                params->lk_params.branch_lengths[i] :
                PLL_OPT_DEFAULT_BRANCH_LEN;
      }
      check_n += num_branch_lengths;
      nbd_ptr += num_branch_lengths;
      l_ptr += num_branch_lengths;
      u_ptr += num_branch_lengths;
    }
    assert(check_n == num_variables);
  }

  /*     We start the iteration by initializing task. */
  *task = (int) START;

  iwa = (int *) calloc (3 * (size_t)num_variables, sizeof(int));
  wa = (double *) calloc (
      (2 * (size_t)max_corrections + 5) * (size_t)num_variables
          + 12 * (size_t)max_corrections * ((size_t)max_corrections + 1),
      sizeof(double));

  int continue_opt = 1;
  while (continue_opt)
  {
    /*     This is the call to the L-BFGS-B code. */
    setulb ((int *)&num_variables, &max_corrections, x, lower_bounds, upper_bounds,
            bound_type, &score, g, &(params->factr), &(params->pgtol), wa, iwa,
            task, &iprint, csave, lsave, isave, dsave);
    if (IS_FG(*task))
    {
      /*
       * the minimization routine has returned to request the
       * function f and gradient g values at the current x.
       * Compute function value f for the sample problem.
       */

      score = compute_lnl_unrooted (params, x);
      if (is_nan(score) || d_equals(score, -INFINITY))
        break;

      double h, temp;
      for (i = 0; i < num_variables; i++)
      {
        temp = x[i];
        h = PLL_LBFGSB_ERROR * fabs (temp);
        if (h < 1e-12)
          h = PLL_LBFGSB_ERROR;

        x[i] = temp + h;
        h = x[i] - temp;
        double lnderiv = compute_lnl_unrooted (params, x);

        g[i] = (lnderiv - score) / h;

        /* reset variable */
        x[i] = temp;
      }

      if (!set_x_to_parameters (params, x))
        return -INFINITY;
    }
    else if (*task != NEW_X)
      continue_opt = 0;
  }

  free (iwa);
  free (wa);
  free (x);
  free (g);
  free (lower_bounds);
  free (upper_bounds);
  free (bound_type);

  if (is_nan(score))
  {
    score = -INFINITY;
    if (!pll_errno)
    {
      pll_errno = PLL_ERROR_LBFGSB_UNKNOWN;
      snprintf(pll_errmsg, 200, "Unknown LBFGSB error");
    }
  }

  return score;
} /* pll_optimize_parameters_lbfgsb */

/******************************************************************************/
/* GENERIC */
/******************************************************************************/

static double recomp_iterative (pll_optimize_options_t * params,
                                pll_utree_t * tree,
                                double prev_lnl)
{
  int i;
  double lnl, new_lnl;
  pll_utree_t *tr_p, *tr_q, *tr_z;

  DBG("Optimizing branch %3d - %3d (%.6f) [%.4f]\n",
      tree->clv_index, tree->back->clv_index, tree->length, prev_lnl);

  lnl = prev_lnl;
  tr_p = tree;
  tr_q = tree->next;
  tr_z = tree->next?tree->next->next:NULL;

  /* set Branch Length */
  assert(d_equals(tree->length, tr_p->back->length));

  params->lk_params.branch_lengths[0]                    = tr_p->length;
  params->lk_params.where.unrooted_t.child_clv_index     = tr_p->back->clv_index;
  params->lk_params.where.unrooted_t.child_scaler_index  = tr_p->back->scaler_index;
  params->lk_params.where.unrooted_t.parent_clv_index    = tr_p->clv_index;
  params->lk_params.where.unrooted_t.parent_scaler_index = tr_p->scaler_index;
  params->lk_params.where.unrooted_t.edge_pmatrix_index  = tr_p->pmatrix_index;

#if(BL_OPT_METHOD == PLL_BRANCH_OPT_LBFGSB)
  new_lnl = -1 * pll_optimize_parameters_lbfgsb(params);
#elif(BL_OPT_METHOD == PLL_BRANCH_OPT_BRENT)
  new_lnl = -1 * pll_optimize_parameters_brent(params);
#else
  new_lnl = -1 * pll_optimize_parameters_newton(params);
#endif

  /* ensure that new_lnl is not NaN */
  assert (new_lnl == new_lnl);

  if (new_lnl < lnl)
  {
    /* revert */
    pll_update_prob_matrices(params->lk_params.partition,
                             params->params_index,
                             &(tr_p->pmatrix_index),
                             &(tr_p->length),
                             1);
//    new_lnl = pll_compute_edge_loglikelihood (params->lk_params.partition,
//                                          tree->back->clv_index, tree->back->scaler_index,
//                                          tree->clv_index, tree->scaler_index,
//                                          tree->pmatrix_index,
//                                          params->lk_params.freqs_index);
    DBG("Revert branch %.14f\n", tree->length);
  }
  else
  {
    /* consolidate */
    DBG("Consolidate branch %.14f  to %.14f\n", tr_p->length, params->lk_params.branch_lengths[0]);
    lnl = new_lnl;
    tree->length = params->lk_params.branch_lengths[0];
    tree->back->length = params->lk_params.branch_lengths[0];
  }

  DBG(" Optimized branch %3d - %3d (%.6f) [%.4f]\n",
      tr_p->clv_index, tr_p->back->clv_index, tr_p->length, lnl);

  /* update children */
  if (tr_p && tr_z)
  {
    /* update children 'Q'
     * CLV at P is recomputed with children P->back and Z->back
     * Scaler is updated by subtracting Q->back and adding P->back
     */
    pll_operation_t new_op;

    /* set CLV */
    new_op.parent_clv_index    = tr_p->clv_index;
    new_op.parent_scaler_index = tr_p->scaler_index;
    new_op.child1_clv_index    = tr_p->back->clv_index;
    new_op.child1_matrix_index = tr_p->back->pmatrix_index;
    new_op.child1_scaler_index = tr_p->back->scaler_index;
    new_op.child2_clv_index    = tr_z->back->clv_index;
    new_op.child2_matrix_index = tr_z->back->pmatrix_index;
    new_op.child2_scaler_index = tr_z->back->scaler_index;
#if(UPDATE_SCALERS)
    /* update scalers */
    if (tr_p->scaler_index != PLL_SCALE_BUFFER_NONE)
    {
      int n = params->lk_params.partition->sites;
      for (i=0; i<n; i++)
      {
        params->lk_params.partition->scale_buffer[tr_p->scaler_index][i] =
            params->lk_params.partition->scale_buffer[tr_p->scaler_index][i]
                + ((tr_p->back->scaler_index != PLL_SCALE_BUFFER_NONE) ?
                    params->lk_params.partition->scale_buffer[tr_p->back->scaler_index][i] :
                    0)
                - ((tr_q->back->scaler_index != PLL_SCALE_BUFFER_NONE) ?
                    params->lk_params.partition->scale_buffer[tr_q->back->scaler_index][i] :
                    0);
      }
    }
#endif
    pll_update_partials (params->lk_params.partition, &new_op, 1);

    /* eval */
    lnl = recomp_iterative (params, tr_q->back, lnl);

    /* update children 'Z'
     * CLV at P is recomputed with children P->back and Q->back
     * Scaler is updated by subtracting Z->back and adding Q->back
     */

    /* set CLV */
    new_op.parent_clv_index    = tr_p->clv_index;
    new_op.parent_scaler_index = tr_p->scaler_index;
    new_op.child1_clv_index    = tr_p->back->clv_index;
    new_op.child1_matrix_index = tr_p->back->pmatrix_index;
    new_op.child1_scaler_index = tr_p->back->scaler_index;
    new_op.child2_clv_index    = tr_q->back->clv_index;
    new_op.child2_matrix_index = tr_q->back->pmatrix_index;
    new_op.child2_scaler_index = tr_q->back->scaler_index;
#if(UPDATE_SCALERS)
    /* update scalers */
    if (tr_p->scaler_index != PLL_SCALE_BUFFER_NONE)
    {
      int n = params->lk_params.partition->sites;
      for (i = 0; i < n; i++)
        params->lk_params.partition->scale_buffer[tr_p->scaler_index][i] =
            params->lk_params.partition->scale_buffer[tr_p->scaler_index][i]
                + ((tr_q->back->scaler_index != PLL_SCALE_BUFFER_NONE) ?
                    params->lk_params.partition->scale_buffer[tr_q->back->scaler_index][i] :
                    0)
                - ((tr_z->back->scaler_index != PLL_SCALE_BUFFER_NONE) ?
                    params->lk_params.partition->scale_buffer[tr_z->back->scaler_index][i] :
                    0);
    }
#endif
    pll_update_partials (params->lk_params.partition, &new_op, 1);

   /* eval */
    lnl = recomp_iterative (params, tr_z->back, lnl);

    /* reset to initial state
     * CLV at P is recomputed with children Q->back and Z->back
     * Scaler is updated by subtracting P->back and adding Z->back
     */

    /* reset CLV */
    new_op.parent_clv_index    = tr_p->clv_index;
    new_op.parent_scaler_index = tr_p->scaler_index;
    new_op.child1_clv_index    = tr_q->back->clv_index;
    new_op.child1_matrix_index = tr_q->back->pmatrix_index;
    new_op.child1_scaler_index = tr_q->back->scaler_index;
    new_op.child2_clv_index    = tr_z->back->clv_index;
    new_op.child2_matrix_index = tr_z->back->pmatrix_index;
    new_op.child2_scaler_index = tr_z->back->scaler_index;
#if(UPDATE_SCALERS)
    /* update scalers */
    if (tr_p->scaler_index != PLL_SCALE_BUFFER_NONE)
    {
      int n = params->lk_params.partition->sites;
      for (i = 0; i < n; i++)
        params->lk_params.partition->scale_buffer[tr_p->scaler_index][i] =
            params->lk_params.partition->scale_buffer[tr_p->scaler_index][i]
                + ((tr_q->back->scaler_index != PLL_SCALE_BUFFER_NONE) ?
                    params->lk_params.partition->scale_buffer[tr_q->back->scaler_index][i] :
                    0)
                - ((tr_p->back->scaler_index != PLL_SCALE_BUFFER_NONE) ?
                    params->lk_params.partition->scale_buffer[tr_p->back->scaler_index][i] :
                    0);
    }
#endif
    pll_update_partials (params->lk_params.partition, &new_op, 1);
  }

  return lnl;
} /* recomp_iterative */

PLL_EXPORT double pll_optimize_branch_lengths_iterative (
                                              pll_optimize_options_t * params,
                                              pll_utree_t * tree,
                                              int smoothings)
{
  int i;
  double lnl = 0.0;

  params->which_parameters = PLL_PARAMETER_BRANCHES_SINGLE;

  lnl = pll_compute_edge_loglikelihood (params->lk_params.partition,
                                            tree->back->clv_index, tree->back->scaler_index,
                                              tree->clv_index, tree->scaler_index,
                                              tree->pmatrix_index,
                                              params->lk_params.freqs_index);

  for (i=0; i<smoothings; i++)
  {
      double new_lnl = recomp_iterative (params, tree, lnl);
      new_lnl = recomp_iterative (params, tree->back, new_lnl);
      if (fabs(new_lnl - lnl) < params->pgtol)
        break;
      lnl = new_lnl;
  }

  return -1*lnl;
} /* pll_optimize_branch_lengths_iterative */




/******************************************************************************/
/* NEWTON-RAPHSON OPTIMIZATION */
/******************************************************************************/

PLL_EXPORT double pll_minimize_newton(pll_optimize_options_t * params,
                                      double x1,
                                      double xguess,
                                      double x2,
                                      unsigned int max_iters,
                                      double *score)
{
  unsigned int i;
  double df, dx, dxold, f;
  double temp, xh, xl, rts, rts_old;

  double tolerance = 1e-4; //params->pgtol
  if (params->which_parameters != PLL_PARAMETER_BRANCHES_SINGLE)
  {
    snprintf(pll_errmsg, 200,
             "Newton-Raphson defined only for single branches");
    pll_errno = PLL_ERROR_PARAMETER;
    return 0.0;
  }
  pll_errno = 0;

  rts = xguess;
  if (rts < x1)
    rts = x1;
  if (rts > x2)
    rts = x2;

  pll_compute_likelihood_derivatives (
      params->lk_params.partition,
      params->lk_params.where.unrooted_t.parent_clv_index,
      params->lk_params.where.unrooted_t.parent_scaler_index,
      params->lk_params.where.unrooted_t.child_clv_index,
      params->lk_params.where.unrooted_t.child_scaler_index,
      rts,
      params->params_index,
      params->lk_params.freqs_index,
      &f, &df);
  DBG("[NR deriv] BL=%f   f=%f  df=%f  nextBL=%f\n", rts, f, df, rts-f/df);
  if (!isfinite(f) || !isfinite(df))
  {
    snprintf (pll_errmsg, 200, "wrong likelihood derivatives");
    pll_errno = PLL_ERROR_NEWTON_DERIV;
    return 0.0;
  }
  if (df >= 0.0 && fabs (f) < tolerance)
    return rts;
  if (f < 0.0)
  {
    xl = rts;
    xh = x2;
  }
  else
  {
    xh = rts;
    xl = x1;
  }

  dx = dxold = fabs (xh - xl);
  for (i = 1; i <= max_iters; i++)
  {
    rts_old = rts;
    if ((df <= 0.0) // function is concave
    || (((rts - xh) * df - f) * ((rts - xl) * df - f) >= 0.0) // out of bound
        )
    {
      dxold = dx;
      dx = 0.5 * (xh - xl);
      rts = xl + dx;
      if (xl == rts)
        return rts;
    }
    else
    {
      dxold = dx;
      dx = f / df;
      temp = rts;
      rts -= dx;
      if (temp == rts)
        return rts;
    }
    if (fabs (dx) < tolerance || (i == max_iters))
      return rts_old;

    if (rts < x1) rts = x1;
    *score = pll_compute_likelihood_derivatives (
          params->lk_params.partition,
          params->lk_params.where.unrooted_t.parent_clv_index,
          params->lk_params.where.unrooted_t.parent_scaler_index,
          params->lk_params.where.unrooted_t.child_clv_index,
          params->lk_params.where.unrooted_t.child_scaler_index,
          rts,
          params->params_index,
          params->lk_params.freqs_index,
          &f, &df);

    if (!isfinite(f) || !isfinite(df))
    {
      snprintf (pll_errmsg, 200, "wrong likelihood derivatives [it]");
      pll_errno = PLL_ERROR_NEWTON_DERIV;
      return 0.0;;
    }

    if (df > 0.0 && fabs (f) < tolerance)
      return rts;

    if (f < 0.0)
      xl = rts;
    else
      xh = rts;
  }

  snprintf(pll_errmsg, 200, "Exceeded maximum number of iterations");
  pll_errno = PLL_ERROR_NEWTON_LIMIT;
  return 0.0;
}
