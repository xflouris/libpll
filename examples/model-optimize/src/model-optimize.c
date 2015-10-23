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
#include "utils.h"

#define RATE_CATS 4

/* parameters to optimize */
#define OPTIMIZE_BRANCHES     1
#define OPTIMIZE_SUBST_PARAMS 1
#define OPTIMIZE_ALPHA        1
#define OPTIMIZE_FREQS        1
#define OPTIMIZE_PINV         0

/* tolerances */
#define OPT_EPSILON       1e-1
#define OPT_PARAM_EPSILON 1e-1

/* use L-BFGS-B instead of Brent for single parameters (alpha / pInv) */
#define USE_LBFGSB 0

/* if CHECK_LOCAL_CONVERGENCE, the parameters are no longer optimize when
 * they do not improve the likelihood at one iteration */
#define CHECK_LOCAL_CONVERGENCE 1

/* if KEEP_UPDATE, branch lengths are iteratively updated during the
 * optimization process */
#define KEEP_UPDATE 1

static int v_int_max (int * v, int n)
{
  int i, max = v[0];
  for (i = 1; i < n; i++)
    if (v[i] > max)
      max = v[i];
  return max;
}

int main (int argc, char * argv[])
{
  /* iterators */
  unsigned int i;

  /* tree properties */
  pll_utree_t     * tree           = NULL;
  pll_partition_t * partition      = NULL;
  pll_operation_t * operations     = NULL;
  double          * branch_lengths = NULL;
  unsigned int    * matrix_indices = NULL;
  pll_utree_t    ** travbuffer     = NULL;
  unsigned int    * data           = NULL;
  unsigned int      matrix_count,
                    ops_count;
  unsigned int tip_count,
               nodes_count,
               branch_count,
               inner_nodes_count;

  /* optimization parameters */
  pll_optimize_options_t params;
  int parameters_to_optimize;
  int * subst_params_symmetries;
  int n_subst_free_params;

  time_t start_time, end_time;

  /*
   * The initialization part is similar to the newick-fasta-unrooted example.
   * It is vaguely commented here.
   * Please refer to the newick-fasta-unrooted example for details.
   */

  if (argc != 4)
    fatal (" syntax: %s [newick] [fasta] [model]", argv[0]);

  printf ("Parsing tree\n");

  /* parse the input tree */
  tree = pll_utree_parse_newick (argv[1], &tip_count);
  nodes_count = 2 * tip_count - 2;
  branch_count = 2 * tip_count - 3;
  inner_nodes_count = tip_count - 2;

  printf ("Traversing branches\n");

  /* fix missing branch lengths */
  set_missing_branch_length (tree, 0.1);

  /*  obtain an array of pointers to tip nodes */
  pll_utree_t ** tipnodes = (pll_utree_t **) calloc ((size_t) tip_count,
                                                     sizeof(pll_utree_t *));
  pll_utree_query_tipnodes (tree, tipnodes);
  pll_utree_t ** innernodes = (pll_utree_t **) calloc (
      (size_t) inner_nodes_count, sizeof(pll_utree_t *));
  pll_utree_query_innernodes (tree, innernodes);

  /* create and populate the list of tipnames */
  char **tipnames = (char **) malloc (tip_count * sizeof(char *));
  data = (unsigned int *) malloc ((size_t) tip_count * sizeof(unsigned int));
  for (i = 0; i < tip_count; ++i)
    tipnames[i] = tipnodes[i]->label;

  printf ("Reading fasta file\n");

  /* create the partition instance */
  partition = partition_fasta_create (argv[2],
                                      STATES,
                                      1,
                                      RATE_CATS,
                                      PLL_ATTRIB_ARCH_AVX, 0, tip_count,
                                      (const char **) tipnames);

  if (!partition)
    fatal ("Error %d: %s\n", pll_errno, pll_errmsg);

  printf ("  Tips:   %u\n", partition->tips);
  printf ("  Length: %u\n", partition->sites);

  /* we no longer need these arrays */
  free (data);
  free (tipnodes);
  free (tipnames);

  travbuffer = (pll_utree_t **) malloc (
      (size_t) nodes_count * sizeof(pll_utree_t *));
  branch_lengths = (double *) malloc ((size_t) branch_count * sizeof(double));
  matrix_indices = (unsigned int *) malloc (
      (size_t) branch_count * sizeof(unsigned int));
  operations = (pll_operation_t *) malloc (
      (size_t) inner_nodes_count * sizeof(pll_operation_t));

  /* perform a postorder traversal of the unrooted tree */
  unsigned int traversal_size;
  if (!pll_utree_traverse (tree, cb_full_traversal, travbuffer,
                           &traversal_size))
    fatal ("Function pll_utree_traverse() requires inner nodes as parameters");

  /* given the computed traversal descriptor, generate the operations
   structure, and the corresponding probability matrix indices that
   may need recomputing */
  pll_utree_create_operations (travbuffer, traversal_size, branch_lengths,
                               matrix_indices, operations, &matrix_count,
                               &ops_count);

  /* we'll use 4 rate categories, and currently initialize them to 0 */
  double rate_cats[RATE_CATS] =
    { 0 };

  /* compute the discretized category rates from a gamma distribution
   with alpha shape 1 and store them in rate_cats  */
  pll_compute_gamma_cats (1, RATE_CATS, rate_cats);

  /* set frequencies at model with index 0 (we currently have only one model) */
  double * empirical_freqs = pll_compute_empirical_frequencies (partition);
  printf ("Empirical frequencies: ");
  for (i = 0; i < STATES; i++)
    printf ("%.4f ", empirical_freqs[i]);
  printf ("\n");
  pll_set_frequencies (partition, 0, 0, empirical_freqs);
  free (empirical_freqs);

  /* Additionally, for these examples we define the substitution rate matrix
   * symmetries. This a six-integer code that defines how the substitution rate
   * parameters are linked to each other. For example, JC/F81 model is 000000,
   * HKY/K80 = 010010 and GTR/SYM = 012345. There is no need to start in 0.
   * For example, 174673 would work exactly as 010010.
   */
  subst_params_symmetries = build_model_symmetries (argv[3]);
  if (!subst_params_symmetries)
    fatal ("Error: Invalid matrix symmetries %s\n", argv[3]);

  printf ("Model: ");
  for (i = 0; i < SUBST_PARAMS; i++)
    printf ("%d", subst_params_symmetries[i]);
  printf ("\n");

  n_subst_free_params = v_int_max (
                          subst_params_symmetries,
                          (int) SUBST_PARAMS);

  /* set empirical substitution rates for GTR model */
  if (n_subst_free_params == 5)
  {
    double * empirical_subst_rates = pll_compute_empirical_subst_rates (
        partition);
    printf ("Empirical rates: ");
    for (i = 0; i < SUBST_PARAMS; i++)
      printf ("%.4f ", empirical_subst_rates[i]);
    printf ("\n");
    pll_set_subst_params (partition, 0, 0, empirical_subst_rates);
    free (empirical_subst_rates);
  }
  else
  {
    /* if there are symmetries in the rate matrix, we initialize it to
     * {1,1,1,1,1,1} */
    double start_subst_rates[SUBST_PARAMS] = {1,1,1,1,1,1};
    pll_set_subst_params (partition, 0, 0, start_subst_rates);
  }

  /* set rate categories */
  pll_set_category_rates (partition, rate_cats);

  printf ("Updating probability matrices\n");

  pll_update_prob_matrices (partition, 0, matrix_indices, branch_lengths,
                            branch_count);

  printf ("Updating partials\n");

  pll_update_partials (partition, operations, inner_nodes_count);

  double logl = pll_compute_edge_loglikelihood (partition, tree->clv_index,
                                                tree->scaler_index,
                                                tree->back->clv_index,
                                                tree->back->scaler_index,
                                                tree->pmatrix_index, 0);

  char * newick = pll_utree_export_newick (tree);
  printf ("Starting tree: %s\n", newick);
  free (newick);
  printf ("Log-L: %f\n", logl);

  /* pll stuff */
  params.lk_params.partition = partition;
  params.lk_params.operations = operations;
  params.lk_params.branch_lengths = branch_lengths;
  params.lk_params.matrix_indices = matrix_indices;
  params.lk_params.alpha_value = 0;
  params.lk_params.freqs_index = 0;
  params.lk_params.rooted = 0;
  params.lk_params.where.unrooted_t.parent_clv_index = tree->clv_index;
  params.lk_params.where.unrooted_t.parent_scaler_index = tree->scaler_index;
  params.lk_params.where.unrooted_t.child_clv_index = tree->back->clv_index;
  params.lk_params.where.unrooted_t.child_scaler_index =
      tree->back->scaler_index;
  params.lk_params.where.unrooted_t.edge_pmatrix_index = tree->pmatrix_index;

  /* optimization parameters */
  params.params_index = 0;
  params.mixture_index = 0;
  params.subst_params_symmetries = subst_params_symmetries;
  params.factr = 1e7;
  params.pgtol = OPT_PARAM_EPSILON;

  parameters_to_optimize =
      (OPTIMIZE_SUBST_PARAMS * PLL_PARAMETER_SUBST_RATES * n_subst_free_params > 0)
      | (OPTIMIZE_ALPHA * PLL_PARAMETER_ALPHA)
      | (OPTIMIZE_BRANCHES * PLL_PARAMETER_BRANCHES_ALL)
      | (OPTIMIZE_PINV * PLL_PARAMETER_PINV)
      | (OPTIMIZE_FREQS * PLL_PARAMETER_FREQUENCIES);

  start_time = time (NULL);

  double cur_logl = logl - 10;
  int smoothings = 1;

  printf ("\nParameter optimization start\n");
  double lnl_monitor = logl;
  while (fabs (cur_logl - logl) > OPT_EPSILON)
  {
    logl = cur_logl;

    if (parameters_to_optimize & PLL_PARAMETER_BRANCHES_ALL)
    {
      /* move to random node */
      int inner_index = rand () % inner_nodes_count;

      tree = innernodes[inner_index];
      pll_utree_traverse (tree, cb_full_traversal, travbuffer, &traversal_size);
      pll_utree_create_operations (travbuffer, traversal_size, branch_lengths,
                                   matrix_indices, operations, &matrix_count,
                                   &ops_count);
      pll_update_prob_matrices (partition, 0, matrix_indices, branch_lengths,
                                branch_count);
      pll_update_partials (partition, operations, tip_count - 2);

      params.which_parameters = PLL_PARAMETER_BRANCHES_SINGLE;

      cur_logl = -1 * pll_optimize_branch_lengths_iterative (
          partition, tree, params.params_index, params.lk_params.freqs_index,
          params.pgtol, smoothings++, KEEP_UPDATE);

      pll_utree_traverse (tree, cb_full_traversal, travbuffer, &traversal_size);
      pll_utree_create_operations (travbuffer, traversal_size, branch_lengths,
                                   matrix_indices, operations, &matrix_count,
                                   &ops_count);
      params.lk_params.where.unrooted_t.parent_clv_index = tree->clv_index;
      params.lk_params.where.unrooted_t.parent_scaler_index =
          tree->scaler_index;
      params.lk_params.where.unrooted_t.child_clv_index = tree->back->clv_index;
      params.lk_params.where.unrooted_t.child_scaler_index =
          tree->back->scaler_index;
      params.lk_params.where.unrooted_t.edge_pmatrix_index =
          tree->pmatrix_index;

      /* if the branch lengths are not updated during the optimization,
       * we need to update them together and recompute the CLVs */
      if (!KEEP_UPDATE)
      {
        pll_update_prob_matrices (partition, 0, matrix_indices, branch_lengths,
                                  2 * tip_count - 3);
        pll_update_partials (partition, operations, tip_count - 2);
        cur_logl = pll_compute_edge_loglikelihood (partition, tree->clv_index,
                                                   tree->scaler_index,
                                                   tree->back->clv_index,
                                                   tree->back->scaler_index,
                                                   tree->pmatrix_index, 0);
      }

      printf ("  %5ld s [branches]: %f\n", time (NULL) - start_time, cur_logl);

      assert(cur_logl >= lnl_monitor);

      if (CHECK_LOCAL_CONVERGENCE && cur_logl - lnl_monitor < OPT_EPSILON)
      {
        parameters_to_optimize &= ~PLL_PARAMETER_BRANCHES_ALL;
      }

      lnl_monitor = cur_logl;
    }

    if (parameters_to_optimize & PLL_PARAMETER_SUBST_RATES)
    {
      params.which_parameters = PLL_PARAMETER_SUBST_RATES;
      cur_logl = -1 * pll_optimize_parameters_lbfgsb (&params);

      printf ("  %5ld s [s_rates]: %f\n", time (NULL) - start_time, cur_logl);
      printf ("             %f %f %f %f %f %f\n", partition->subst_params[0][0],
              partition->subst_params[0][1], partition->subst_params[0][2],
              partition->subst_params[0][3], partition->subst_params[0][4],
              partition->subst_params[0][5]);

      assert(cur_logl >= lnl_monitor);

      if (CHECK_LOCAL_CONVERGENCE && cur_logl - lnl_monitor < OPT_EPSILON)
              parameters_to_optimize &= ~PLL_PARAMETER_SUBST_RATES;

      lnl_monitor = cur_logl;
    }

    if (parameters_to_optimize & PLL_PARAMETER_ALPHA)
    {
      params.which_parameters = PLL_PARAMETER_ALPHA;
#if(USE_LBFGSB)
      cur_logl = -1 * pll_optimize_parameters_lbfgsb (&params);
#else
      cur_logl = -1 * pll_optimize_parameters_brent (&params);
#endif
      printf ("  %5ld s [alpha]: %f\n", time (NULL) - start_time, cur_logl);
      printf ("             %f\n", params.lk_params.alpha_value);

      assert(cur_logl >= lnl_monitor);

      if (CHECK_LOCAL_CONVERGENCE && cur_logl - lnl_monitor < OPT_EPSILON)
        parameters_to_optimize &= ~PLL_PARAMETER_ALPHA;

      lnl_monitor = cur_logl;
    }

    if (parameters_to_optimize & PLL_PARAMETER_PINV)
    {
      params.which_parameters = PLL_PARAMETER_PINV;
#if(USE_LBFGSB)
      cur_logl = -1 * pll_optimize_parameters_lbfgsb (&params);
#else
      cur_logl = -1 * pll_optimize_parameters_brent (&params);
#endif
      printf ("  %5ld s [p-inv]: %f\n", time (NULL) - start_time, cur_logl);
      printf ("             %f\n", partition->prop_invar[0]);

      assert(cur_logl >= lnl_monitor);

      if (CHECK_LOCAL_CONVERGENCE && cur_logl - lnl_monitor < OPT_EPSILON)
        parameters_to_optimize &= ~PLL_PARAMETER_PINV;

      lnl_monitor = cur_logl;
    }

    if (parameters_to_optimize & PLL_PARAMETER_FREQUENCIES)
    {
      params.which_parameters = PLL_PARAMETER_FREQUENCIES;
      cur_logl = -1 * pll_optimize_parameters_lbfgsb (&params);
      printf ("  %5ld s [freqs]: %f\n", time (NULL) - start_time, cur_logl);
      printf ("             ");
      for (i = 0; i < partition->states; i++)
        printf ("%f ", partition->frequencies[0][i]);
      printf ("\n");

      assert(cur_logl >= lnl_monitor);

      if (CHECK_LOCAL_CONVERGENCE && cur_logl - lnl_monitor < OPT_EPSILON)
        parameters_to_optimize &= ~PLL_PARAMETER_FREQUENCIES;

      lnl_monitor = cur_logl;
    }

    printf ("Iteration: %5ld s. : %f\n", time (NULL) - start_time, cur_logl);
  }
  end_time = time (NULL);

  printf ("Final Log-L: %f\n", cur_logl);
  printf ("Time:  %ld s.\n", end_time - start_time);

  printf ("Alpha: %f\n", params.lk_params.alpha_value);
  printf ("P-inv: %f\n", partition->prop_invar[0]);
  printf ("Rates: %f %f %f %f %f %f\n", partition->subst_params[0][0],
          partition->subst_params[0][1], partition->subst_params[0][2],
          partition->subst_params[0][3], partition->subst_params[0][4],
          partition->subst_params[0][5]);

  newick = pll_utree_export_newick (tree);
  printf ("\nFinal tree: %s\n", newick);
  free (newick);
  printf ("Final Log-L: %f\n", logl);

  /* clean */

  free (innernodes);
  pll_utree_destroy (tree);
  pll_partition_destroy (partition);

  free (travbuffer);
  free (subst_params_symmetries);
  free (branch_lengths);
  free (matrix_indices);
  free (operations);

  return (0);
}
