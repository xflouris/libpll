/*
 * optimize.c
 *
 *  Created on: May 12, 2015
 *      Author: diego
 */
#include "optimize.h"

#include <lbfgs.h>
#include <stdio.h>

static lbfgsfloatval_t targetFunk(double *x,
                                  lk_stuff * tree) {
  /* compute positive logLK */
  switch (tree->which_parameter)
    {
    case PARAM_SUBST_RATES:
      pll_set_subst_params (tree->partition, tree->params_index, x, 5);
      break;
    case PARAM_ALPHA:
      pll_compute_gamma_cats(x[0], tree->n_cat_gamma, tree->partition->rates);
      break;
    case PARAM_PINV:
      pll_update_invariant_sites_proportion(tree->partition, x[0]);
      break;
    case PARAM_FREQUENCIES:
      pll_set_frequencies(tree->partition, tree->params_index, x);
      break;
    default:
      assert(0);
      break;
    }

  pll_update_prob_matrices (tree->partition, 0, tree->matrix_indices,
                            tree->branch_lengths, tree->n_prob_matrices);
  pll_update_partials (tree->partition, tree->operations, tree->n_operations);

  lbfgsfloatval_t logl;
  if (tree->is_rooted)
  {
    logl = pll_compute_root_loglikelihood (tree->partition, tree->clv_index1,
                                           tree->freqs_index);
  }
  else
  {
    logl = pll_compute_edge_loglikelihood (tree->partition, tree->clv_index1,
                                           tree->clv_index2,
                                           tree->branch_matrix_index,
                                           tree->freqs_index);
  }
  return -logl;
}

const double ERROR_X = 1.0e-4;

static lbfgsfloatval_t derivativeFunk (lbfgsfloatval_t * x,
                                       lbfgsfloatval_t * dfx,
                                       int ndim,
                                       lk_stuff * tree)
{
  int dim;

  lbfgsfloatval_t fx = targetFunk(x, tree);
  double h, temp;
  for (dim = 0; dim < ndim; dim++)
  {

    temp = x[dim];
    h = ERROR_X * fabs (temp);
    if (h < 1e-12)
      h = ERROR_X;

    if ((temp+h) < (tree->max_params[dim])) {
      /* search normal */
      x[dim] = temp + h;
      h = x[dim] - temp;
      dfx[dim] = (targetFunk(x, tree) - fx) / h;
    } else if ((temp-h) > (tree->min_params[dim])) {
      /* search reverse */
      x[dim] = temp - h;
      dfx[dim] = (targetFunk(x, tree) - fx) / h;
    }

    /* reset variable */
    x[dim] = temp;

  }

  return fx;
}

static lbfgsfloatval_t evaluate(void *instance,
                                const lbfgsfloatval_t *x,
                                lbfgsfloatval_t *g,
                                const int n,
                                const lbfgsfloatval_t step
                                )
{
  lbfgsfloatval_t fx = 0.0;

  lk_stuff * tree = (lk_stuff *) instance;

  double * xtest = (double *) malloc (n * sizeof(double));
  memcpy (xtest, x, n * sizeof(double));

  fx = derivativeFunk (xtest, g, n, tree);

    free(xtest);

    return fx;
}

static int progress(
    void *instance,
    const lbfgsfloatval_t *x,
    const lbfgsfloatval_t *g,
    const lbfgsfloatval_t fx,
    const lbfgsfloatval_t xnorm,
    const lbfgsfloatval_t gnorm,
    const lbfgsfloatval_t step,
    int n,
    int k,
    int ls
    )
{

#ifndef DEBUG
  return 0;
#endif

  int dim;
    printf("Iteration %d:\n", k);
    printf("  fx = %.8f, ", fx);
    for (dim=0; dim<n; ++dim)
      printf("%.4f ", x[dim]);
    printf("\n");
    for (dim=0; dim<n; ++dim)
          printf("%.4f ", g[dim]);
        printf("\n");
    printf("xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
    printf("\n");
    return 0;
}

double optimize(lk_stuff * tree, double epsilon)
{
  int i, ret;

  lbfgsfloatval_t fx;
  lbfgsfloatval_t * x;
  lbfgs_parameter_t param;

  int ndim;
  double min_val, max_val; /* assuming homogeneity */
  double * params_buf;

  switch (tree->which_parameter)
    {
    case PARAM_SUBST_RATES:
      ndim = 5;
      min_val = 0.001;
      max_val = 500.0;
      params_buf = tree->partition->subst_params[tree->params_index];
      break;
    case PARAM_ALPHA:
      ndim = 1;
      min_val = 0.002;
      max_val = 1000;
      params_buf = &(tree->alpha);
      break;
    case PARAM_PINV:
      ndim = 1;
      min_val = 0.0;
      max_val = 0.999;
      if (tree->partition->prop_invar[tree->params_index] < 1e-12)
               pll_update_invariant_sites_proportion(tree->partition, 0.5);
      params_buf = &(tree->partition->prop_invar[tree->params_index]);
      break;
    case PARAM_FREQUENCIES:
      ndim = 3;
      min_val = 0.0001;
      max_val = 0.9996;
      params_buf = tree->partition->frequencies[tree->params_index];
      break;
    default:
      assert(0);
      break;
    }
  x = lbfgs_malloc(ndim);
  if (!x)
    {
      printf ("ERROR: Failed to allocate a memory block for variables.\n");
      return 1;
    }

  memcpy(x, params_buf, ndim*sizeof(double));

  tree->min_params = (double *) malloc (ndim * sizeof(double));
  tree->max_params = (double *) malloc (ndim * sizeof(double));
  for (i=0; i<ndim; i++) {
    tree->min_params[i] = min_val;
    tree->max_params[i] = max_val;
  }

//  printf("   INI: ");
//  for (i=0;i<ndim;i++)
//    printf("%f ", x[i]);
//  printf("\n");

  /* Initialize the parameters for the L-BFGS optimization. */
  lbfgs_parameter_init(&param);
param.epsilon = epsilon;
  /*
      Start the L-BFGS optimization; this will invoke the callback functions
      evaluate() and progress() when necessary.
   */
  assert(params_buf);
  ret = lbfgs(ndim, x, tree->min_params, tree->max_params, &fx, evaluate, progress, tree, &param);

  free (tree->min_params);
  free (tree->max_params);

  switch (tree->which_parameter)
    {
    case PARAM_SUBST_RATES:
      break;
    case PARAM_ALPHA:
      tree->alpha = x[0];
      break;
    case PARAM_PINV:
      break;
    case PARAM_FREQUENCIES:
      break;
    default:
      assert(0);
      break;
    }

  /* Report the result. */
//  printf ("  L-BFGS optimization terminated with status code = %d\n", ret);
//  printf ("  final = %.8f\n", (double) fx);
//  printf ("   END: ");
//  for (i = 0; i < ndim; i++)
//    printf ("%f ", x[i]);
//  printf ("\n\n");

  lbfgs_free (x);

  return fx;
}
