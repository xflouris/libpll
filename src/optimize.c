#include "pll.h"

/*****************************************************
 Multi dimensional optimization with BFGS method
 *****************************************************/

typedef int bool;
#define true  1
#define false 0

typedef struct {
  pll_partition_t * partition;
  int params_index;
  pll_operation_t * operations;
  int * matrix_indices;
  double * branch_lengths;
  int n_prob_matrices;
  int n_operations;
  int is_rooted;
  int freqs_index;
  int clv_index1;
  /* additional stuff for unrooted trees */
  int clv_index2;
  int branch_matrix_index;
  int ndim;
} lk_stuff;

static double derivativeFunk (double x[], double dfx[], lk_stuff tree);
static void dfpmin (double p[], int n, double lower[], double upper[], double gtol, int *iter, double *fret, lk_stuff tree);
static void lnsrch (int n, double xold[], double fold, double g[], double p[], double x[],
                    double *f, double stpmax, int *check, double lower[], double upper[],
                    lk_stuff tree);
static void fixBound (double x[], double lower[], double upper[], int n);
static double targetFunk (double x[], lk_stuff tree);
static void getVariables(double *variables, pll_partition_t * partition, int params_index, int ndim);
static void setVariables(double *variables, pll_partition_t * partition, int params_index, int ndim);
double minimizeMultiDimen (double guess[], int ndim, double lower[], double upper[],
                    bool bound_check[], double gtol, lk_stuff tree);

const double ERROR_X = 1.0e-4;

#define ALF 1.0e-4
#define TOLX 1.0e-7

static double targetFunk (double x[], lk_stuff tree)
{
  /* compute positive logLK */
  pll_set_subst_params (tree.partition, 0, x + 1, tree.ndim);

  pll_update_prob_matrices (tree.partition, 0, tree.matrix_indices,
                            tree.branch_lengths, tree.n_prob_matrices);
  pll_update_partials (tree.partition, tree.operations, tree.n_operations);

  double logl;
  if (tree.is_rooted)
  {
    logl = pll_compute_root_loglikelihood (tree.partition, tree.clv_index1,
                                           tree.freqs_index);
  }
  else
  {
    logl = pll_compute_edge_loglikelihood (tree.partition,
                                           tree.clv_index1, tree.clv_index2,
                                           tree.branch_matrix_index,
                                           tree.freqs_index);
  }

  return -logl;
  //setVariables (x, partition, params_index);
}

static double random_double()
{
  return ((double) rand()) / ((double) RAND_MAX + 1);
}

static double inline
FMAX (double a, double b)
{
  if (a > b)
    return a;
  return b;
}

static void
fixBound (double x[], double lower[], double upper[], int n)
{
  int i;
  for (i = 1; i <= n; i++)
  {
    if (x[i] < lower[i])
      x[i] = lower[i];
    else if (x[i] > upper[i])
      x[i] = upper[i];
  }
}

static void
lnsrch (int n, double xold[], double fold, double g[], double p[], double x[],
        double *f, double stpmax, int *check, double lower[], double upper[],
        lk_stuff tree)
{
  int i;
  double a, alam, alam2 = 0, alamin, b, disc, f2 = 0, fold2 = 0, rhs1, rhs2,
      slope, sum, temp, test, tmplam;

  *check = 0;
  for (sum = 0.0, i = 1; i <= n; i++)
    sum += p[i] * p[i];
  sum = sqrt (sum);
  if (sum > stpmax)
    for (i = 1; i <= n; i++)
      p[i] *= stpmax / sum;
  for (slope = 0.0, i = 1; i <= n; i++)
    slope += g[i] * p[i];
  test = 0.0;
  for (i = 1; i <= n; i++)
  {
    temp = fabs (p[i]) / FMAX (fabs (xold[i]), 1.0);
    if (temp > test)
      test = temp;
  }
  alamin = TOLX / test;
  alam = 1.0;
  /*
   int rep = 0;
   do {
   for (i=1;i<=n;i++) x[i]=xold[i]+alam*p[i];
   if (!checkRange(x))
   alam *= 0.5;
   else
   break;
   rep++;
   } while (rep < 10);
   */
  bool first_time = true;
  for (;;)
  {
    for (i = 1; i <= n; i++)
      x[i] = xold[i] + alam * p[i];
    fixBound (x, lower, upper, n);
    //checkRange(x);
    *f = targetFunk (x, tree);
    if (alam < alamin)
    {
      for (i = 1; i <= n; i++)
        x[i] = xold[i];
      *check = 1;
      return;
    }
    else if (*f <= fold + ALF * alam * slope)
      return;
    else
    {
      if (first_time)
        tmplam = -slope / (2.0 * (*f - fold - slope));
      else
      {
        rhs1 = *f - fold - alam * slope;
        rhs2 = f2 - fold2 - alam2 * slope;
        a = (rhs1 / (alam * alam) - rhs2 / (alam2 * alam2)) / (alam - alam2);
        b = (-alam2 * rhs1 / (alam * alam) + alam * rhs2 / (alam2 * alam2))
            / (alam - alam2);
        if (a == 0.0)
          tmplam = -slope / (2.0 * b);
        else
        {
          disc = b * b - 3.0 * a * slope;
          if (disc < 0.0) //nrerror("Roundoff problem in lnsrch.");
            tmplam = 0.5 * alam;
          else if (b <= 0.0)
            tmplam = (-b + sqrt (disc)) / (3.0 * a);
          else
            tmplam = -slope / (b + sqrt (disc));
        }
        if (tmplam > 0.5 * alam)
          tmplam = 0.5 * alam;
      }
    }
    alam2 = alam;
    f2 = *f;
    fold2 = fold;
    alam = FMAX (tmplam, 0.1 * alam);
    first_time = false;
  }
}
#undef ALF
#undef TOLX

const int MAX_ITER = 3;

double
minimizeMultiDimen (double guess[], int ndim, double lower[], double upper[],
                    bool bound_check[], double gtol, lk_stuff tree)
{
  int i, iter;
  double fret, minf = 10000000.0;
  double *minx = (double *) malloc ((size_t) (ndim + 1) * sizeof(double));
  int count = 0;
  bool restart;
  do
  {
    dfpmin (guess, ndim, lower, upper, gtol, &iter, &fret, tree);
    if (fret < minf)
    {
      minf = fret;
      for (i = 1; i <= ndim; i++)
        minx[i] = guess[i];
    }
    count++;
    // restart the search if at the boundary
    // it's likely to end at a local optimum at the boundary
    restart = false;

    for (i = 1; i <= ndim; i++)
      if (bound_check[i])
        if (fabs (guess[i] - lower[i]) < 1e-4
            || fabs (guess[i] - upper[i]) < 1e-4)
        {
          restart = true;
          break;
        }

    if (!restart)
      break;

    if (count == MAX_ITER)
      break;

    do
    {
      for (i = 1; i <= ndim; i++)
      {
        guess[i] = random_double () * (upper[i] - lower[i]) / 3 + lower[i];
      }
    }
    while (false);
    printf ("Restart estimation at the boundary... \n");
  }
  while (count < MAX_ITER);
  if (count > 1)
  {
    for (i = 1; i <= ndim; i++)
      guess[i] = minx[i];
    fret = minf;
  }
  free (minx);

  return fret;
}

#define ITMAX 200
#define EPS 3.0e-8
#define TOLX (4*EPS)
#define STPMX 100.0

static void
free_all (double *xi, double *pnew, double **hessin, double *hdg, double *g,
          double *dg, int n)
{
  int i;

  free (xi);
  free (pnew);
  free (hdg);
  free (g);
  free (dg);

  for (i = 0; i <= n; i++)
    free (hessin[i]);

  free (hessin);
}

static void
dfpmin (double p[], int n, double lower[], double upper[], double gtol,
        int *iter, double *fret, lk_stuff tree)
{
  int check, i, its, j;
  double den, fac, fad, fae, fp, stpmax, sum = 0.0, sumdg, sumxi, temp, test;
  double *dg, *g, *hdg, **hessin, *pnew, *xi;

  dg = (double *) malloc (((size_t) n + 1) * sizeof(double));
  g = (double *) malloc (((size_t) n + 1) * sizeof(double));
  hdg = (double *) malloc (((size_t) n + 1) * sizeof(double));
  hessin = (double **) malloc (((size_t) n + 1) * sizeof(double *));
  for (i = 0; i <= n; i++)
    hessin[i] = (double*) malloc (((size_t) n + 1) * sizeof(double));
  pnew = (double *) malloc (((size_t) n + 1) * sizeof(double));
  xi = (double *) malloc (((size_t) n + 1) * sizeof(double));
  fp = derivativeFunk (p, g, tree);

  for (i = 1; i <= n; i++)
  {
    for (j = 1; j <= n; j++)
      hessin[i][j] = 0.0;
    hessin[i][i] = 1.0;
    xi[i] = -g[i];
    sum += p[i] * p[i];
  }
  //checkBound(p, xi, lower, upper, n);
  //checkDirection(p, xi);

  stpmax = STPMX * FMAX (sqrt (sum), (double) n);
  for (its = 1; its <= ITMAX; its++)
  {
    *iter = its;
    lnsrch (n, p, fp, g, xi, pnew, fret, stpmax, &check, lower, upper, tree);
    fp = *fret;
    for (i = 1; i <= n; i++)
    {
      xi[i] = pnew[i] - p[i];
      p[i] = pnew[i];
    }
    test = 0.0;
    for (i = 1; i <= n; i++)
    {
      temp = fabs (xi[i]) / FMAX (fabs (p[i]), 1.0);
      if (temp > test)
        test = temp;
    }
    if (test < TOLX)
    {
      free_all (xi, pnew, hessin, hdg, g, dg, n);
      return;
    }
    for (i = 1; i <= n; i++)
      dg[i] = g[i];
    derivativeFunk (p, g, tree);
    test = 0.0;
    den = FMAX (fabs (*fret), 1.0); // fix bug found by Tung, as also suggested by NR author
    for (i = 1; i <= n; i++)
    {
      temp = fabs (g[i]) * FMAX (fabs (p[i]), 1.0) / den;
      if (temp > test)
        test = temp;
    }
    if (test < gtol)
    {
      free_all (xi, pnew, hessin, hdg, g, dg, n);
      return;
    }
    for (i = 1; i <= n; i++)
      dg[i] = g[i] - dg[i];
    for (i = 1; i <= n; i++)
    {
      hdg[i] = 0.0;
      for (j = 1; j <= n; j++)
        hdg[i] += hessin[i][j] * dg[j];
    }
    fac = fae = sumdg = sumxi = 0.0;
    for (i = 1; i <= n; i++)
    {
      fac += dg[i] * xi[i];
      fae += dg[i] * hdg[i];
      sumdg += (dg[i]*dg[i]);
      sumxi += (xi[i]*xi[i]);
    }
    if (fac * fac > EPS * sumdg * sumxi)
    {
      fac = 1.0 / fac;
      fad = 1.0 / fae;
      for (i = 1; i <= n; i++)
        dg[i] = fac * xi[i] - fad * hdg[i];
      for (i = 1; i <= n; i++)
      {
        for (j = 1; j <= n; j++)
        {
          hessin[i][j] += fac * xi[i] * xi[j] - fad * hdg[i] * hdg[j]
              + fae * dg[i] * dg[j];
        }
      }
    }
    for (i = 1; i <= n; i++)
    {
      xi[i] = 0.0;
      for (j = 1; j <= n; j++)
        xi[i] -= hessin[i][j] * g[j];
    }
    //checkBound(p, xi, lower, upper, n);
    //checkDirection(p, xi);
    //if (*iter > 200) cout << "iteration=" << *iter << endl;
  }
  // BQM: TODO disable this message!
  //nrerror("too many iterations in dfpmin");
  free_all (xi, pnew, hessin, hdg, g, dg, n);
}
#undef ITMAX
#undef SQR
#undef EPS
#undef TOLX
#undef STPMX

/**
 the approximated derivative function
 @param x the input vector x
 @param dfx the derivative at x
 @return the function value at x
 */
static double
derivativeFunk (double x[], double dfx[], lk_stuff tree)
{
  int dim;
  /*
   if (!checkRange(x))
   return INFINITIVE;
   */
  double fx = targetFunk(x, tree);
  int ndim = tree.ndim;
  double h, temp;
  for (dim = 1; dim <= ndim; dim++)
  {
    temp = x[dim];
    h = ERROR_X * fabs (temp);
    if (h == 0.0)
      h = ERROR_X;
    x[dim] = temp + h;
    h = x[dim] - temp;
    dfx[dim] = (targetFunk(x, tree) - fx) / h;
    x[dim] = temp;
  }
  return fx;
}

static void setVariables(double *variables, pll_partition_t * partition, int params_index, int ndim) {
    memcpy(variables+1, partition->subst_params[params_index], ndim * sizeof(double));
    //memcpy(variables+ncategory+1, proportion, (ncategory-1)*sizeof(double));
}

static void getVariables(double *variables, pll_partition_t * partition, int params_index, int ndim) {
    memcpy(partition->subst_params[params_index], variables+1, ndim * sizeof(double));
//    memcpy(proportion, variables+ncategory+1, (ncategory-1)*sizeof(double));
//    double sum = 0.0;
//    for (int i = 0; i < ncategory-1; i++)
//        sum += proportion[i];
//    proportion[ncategory-1] = 1.0 - sum;
}

double
optimizeParameters (double epsilon, int num_cats, int num_subst_rates,
                    pll_partition_t * partition, int params_index,
                    pll_operation_t * operations, double * branch_lengths,
                    int * matrix_indices, int n_prob_matrices, int n_operations,
                    int is_rooted, int clv_index1, int clv_index2,
                    int branch_matrix_index, int freqs_index
                    )
{

  // return if nothing to be optimized
  if (num_subst_rates == 0)
    return 0.0;

  double *variables = (double *) malloc (
      (size_t) (num_subst_rates + 1) * sizeof(double));
  double *upper_bound = (double *) malloc (
      (size_t) (num_subst_rates + 1) * sizeof(double));
  double *lower_bound = (double *) malloc (
      (size_t) (num_subst_rates + 1) * sizeof(double));
  bool *bound_check = (bool *) malloc (
      (size_t) (num_subst_rates + 1) * sizeof(bool));
  int i;
  double score;

  // by BFGS algorithm
  setVariables (variables, partition, params_index, num_subst_rates);
  for (i = 1; i <= num_subst_rates; i++)
  {
    lower_bound[i] = 1e-4;
    upper_bound[i] = 100.0;
    bound_check[i] = false;
  }
  for (i = num_subst_rates - num_cats + 2; i <= num_subst_rates; i++)
    upper_bound[i] = 1.0;

  lk_stuff tree;
  tree.ndim = num_subst_rates;
  tree.partition = partition;
  tree.params_index = params_index;
  tree.operations = operations;
  tree.branch_lengths = branch_lengths;
  tree.matrix_indices = matrix_indices;
  tree.n_prob_matrices = n_prob_matrices;
  tree.n_operations = n_operations;

  tree.is_rooted = is_rooted;
  tree.clv_index1 = clv_index1;
  tree.clv_index2 = clv_index2;
  tree.branch_matrix_index = branch_matrix_index;
  tree.freqs_index = freqs_index;

  score = -minimizeMultiDimen (variables, num_subst_rates, lower_bound,
                               upper_bound, bound_check, epsilon,
                               tree);

  getVariables (variables, partition, params_index, num_subst_rates);

  free (bound_check);
  free (lower_bound);
  free (upper_bound);
  free (variables);

  return score;
}
